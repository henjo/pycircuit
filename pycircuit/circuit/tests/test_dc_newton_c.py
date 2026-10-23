"""The DC operating point's base Newton solve in C (`_dc_newton_c`,
`csrc/dcnewton.c`): every base solve of a gold deck's operating point --
its solution's bytes and iteration count, or its exception's type and
text, in order -- StructuralKLU's counters, the operating point and the
warnings, the same with the C as without; its declines; a planted fault
seen.  Each deck skips where the benchmark data is not fetched."""
import hashlib
import json
import warnings

import numpy as np
import pytest

from pycircuit._testing import benchdata
from pycircuit.circuit import (
    _dc_newton_c,
    _paths,
    _tran_initial,
    circuit,
    linearsolver,
    nrsolver,
    spice_import,
)

_DECKS = {name: 'Netlists/CircuitSim90/' + path for name, path in (
    ('gm6', 'MOS2/gm6.cir'), ('mike2', 'MOS3/mike2.cir'), ('gm3', 'MOS3/gm3.cir'),
    ('gm19', 'MOS3/gm19.cir'))}


## the chain's tests: DC's operating point by the ordinary Newton chain
pytestmark = pytest.mark.usefixtures('ordinary_dc')


def _imported(name):
    text = benchdata.spice_data(_DECKS[name])
    if text is None:
        pytest.skip('benchmark data not fetched (benchmarks/fetch_spice_suite.py)')
    ## (without SPICE's junction gmin: with it mike2 and gm3 converge on the
    ## plain solve at once, and there is no ladder for the C to serve)
    return spice_import.import_netlist(text, gmin=0.0)


def _record(name, monkeypatch, on, x0=None):
    """The deck's operating point, every base solve recorded in order, and
    what the run leaves: `(record, the C's _paths counts)`."""
    circuit.default_toolkit = circuit.numeric
    monkeypatch.setattr(_dc_newton_c, 'ENABLED', on)
    log, sks, out = [], [], {}
    orig = nrsolver.StandardNewton.solve_system

    def solve_system(self, *a, **k):
        try:
            r = orig(self, *a, **k)
        except Exception as e:
            log.append(('raise', type(e).__name__, str(e)))
            raise
        log.append(('ok', hashlib.sha256(np.asarray(r[0]).tobytes()).hexdigest(), r[1]))
        return r
    monkeypatch.setattr(nrsolver.StandardNewton, 'solve_system', solve_system)
    init = linearsolver.StructuralKLU.__init__

    def klu_init(self, *a, **k):
        init(self, *a, **k)
        sks.append(self)
    monkeypatch.setattr(linearsolver.StructuralKLU, '__init__', klu_init)
    op = _tran_initial._InitialState._solve_operating_point

    def operating_point(self, *a, **k):
        x = op(self, *a, **k)
        out['x0'] = np.asarray(x).tobytes()
        return x
    monkeypatch.setattr(_tran_initial._InitialState, '_solve_operating_point', operating_point)
    imp = _imported(name)
    tr, kw = imp.transient()
    kw = dict(kw, tend=kw['timestep'])
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        if x0 is None:
            tr.solve(**kw)
        else:
            from pycircuit.circuit.dcanalysis import DC
            dc = DC(tr.cir, toolkit=circuit.numeric, epar=tr.epar)
            try:
                dc.solve(x0=x0(tr.cir.n))
            except Exception as e:                             # noqa: BLE001
                log.append(('dc raised', type(e).__name__, str(e)))
    counts = {k: v for k, v in _paths.since(before).items() if k.startswith('dcnewton_c:')}
    klu = [(s.refactors, s.factors, s.pivot_fallbacks, s.residual_fallbacks, s.foreign)
           for s in sks]
    rec = json.dumps({'log': log, 'klu': klu, 'x0': out.get('x0', b'').hex(),
                      'warnings': [(w.category.__name__, str(w.message)) for w in W]})
    return rec, counts


@pytest.mark.parametrize('name', ['mike2', 'gm3', 'gm19'])
def test_every_base_solve_is_the_python_one(name, monkeypatch):
    """mike2 and gm3 dense (every layer: junction gmin, gshunt, source
    stepping; failing solves, singular ones), gm19 on StructuralKLU (its
    refactors, zero-pivot and residual fallbacks, singular and flag
    hand-backs): every solve's outcome, the counters, x0 and the warnings
    the same -- and the C served."""
    py, _ = _record(name, monkeypatch, False)
    c, counts = _record(name, monkeypatch, True)
    assert c == py
    assert counts.get('dcnewton_c:served', 0) > 0, counts


@pytest.mark.parametrize('name', ['mike2'])
def test_drawn_starting_points(name, monkeypatch):
    """The operating point from drawn starts -- other rungs, other
    failures -- every solve the same with the C as without."""
    rng = np.random.default_rng(len(name))
    served = 0
    for _ in range(3):
        ## (rough starts: the plain Newton fails and the ladders run -- a
        ## first solve converging is declined, AutoSolver's choice not made)
        scale = 10.0 ** rng.uniform(0, 1)
        seed = int(rng.integers(1 << 30))

        def x0(n, scale=scale, seed=seed):
            return np.random.default_rng(seed).standard_normal(n) * scale
        py, _ = _record(name, monkeypatch, False, x0)
        c, counts = _record(name, monkeypatch, True, x0)
        assert c == py
        served += counts.get('dcnewton_c:served', 0)
    assert served > 0


def test_the_declines(monkeypatch):
    """A circuit below `nrsolver.DC_C_MIN_N` (gm6: 20 unknowns), numpy reporting
    underflow, the switch off: declined and counted, the answer the
    Python's; an evaluation that is not a DC composition is not one."""
    py, _ = _record('gm6', monkeypatch, False)
    got, counts = _record('gm6', monkeypatch, True)
    assert counts.get('dcnewton_c:small', 0) > 0 and 'dcnewton_c:served' not in counts
    assert got == py
    py, _ = _record('mike2', monkeypatch, False)
    with np.errstate(under='warn'):
        got, counts = _record('mike2', monkeypatch, True)
    assert counts.get('dcnewton_c:underflow', 0) > 0 and 'dcnewton_c:served' not in counts
    assert got == py
    _, counts = _record('mike2', monkeypatch, False)
    assert counts.get('dcnewton_c:off', 0) > 0
    assert not _dc_newton_c.composition(lambda x: (x, x))


def test_a_planted_fault_is_seen(monkeypatch):
    """Alive: a layer's conductance one ulp off in the C's inputs (mike2's:
    junction gmin's; gshunt's where it runs) moves the record."""
    py, _ = _record('mike2', monkeypatch, False)
    load = _dc_newton_c._load

    def planted(rec, dc, x0, lam, gs, jg, *a):
        if gs is not None:
            gs = np.nextafter(gs, np.inf)
        if jg is not None:
            jg = (np.nextafter(jg[0], np.inf), jg[1])
        return load(rec, dc, x0, lam, gs, jg, *a)
    monkeypatch.setattr(_dc_newton_c, '_load', planted)
    c, counts = _record('mike2', monkeypatch, True)
    assert counts.get('dcnewton_c:served', 0) > 0
    assert c != py


## -- the rounding floor (roadmap 197) ----------------------------------------

def _companion_dc(cells=16, g=3.5e5):
    """The transient's stalled step as an operating point: per cell a VCVS
    sets `b{k}` to 0.7 V(a) (`a` a divider's node), a behavioural g V(b{k})
    to ground stands in for a capacitor's companion C/h v and a current
    source g (0.7 V(a) + 6e-10 (k + 1)) for its history -- the VCVS's
    current their ~2e-4 A difference, its update at the rounding of its
    residual (`nrsolver._within_rounding`); 16 cells: past `DC_C_MIN_N`."""
    from pycircuit.circuit.elements import IS, VCVS, VS, R, SubCircuit, gnd
    from pycircuit.circuit.elements_hdl import BSourceHdl
    c = SubCircuit()
    c['vs'] = VS('s', gnd, v=3.3)
    c['r1'] = R('s', 'a', r=1e3)
    c['r2'] = R('a', gnd, r=3e3)
    for k in range(cells):
        c[f'e{k}'] = VCVS('a', gnd, f'b{k}', gnd, g=0.7)
        c[f'gq{k}'] = BSourceHdl(f'b{k}', gnd, f'b{k}', gnd, i_func=lambda v: g * v)
        c[f'ih{k}'] = IS(gnd, f'b{k}', i=g * (0.7 * 2.475 + 6e-10 * (k + 1)))
    return c


def _floor_dc(cells, on, within=None):
    """`_companion_dc`'s operating point twice -- the second from zeros, a
    solve the C may serve (the first chose the solver) -- its bytes or its
    exception's type, the second's counts, and the solver chosen; `within`
    stands in for `nrsolver._within_rounding`."""
    from pycircuit.circuit.dcanalysis import DC
    kept = _dc_newton_c.ENABLED, nrsolver._within_rounding
    _dc_newton_c.ENABLED = on
    if within is not None:
        nrsolver._within_rounding = within
    out = []
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dc = DC(_companion_dc(cells), toolkit=circuit.numeric, reltol=1e-9, iabstol=1e-13,
                    vabstol=1e-12)
            for k in range(2):
                before = _paths.snapshot()
                try:
                    out.append(np.asarray(dc.solve(x0=np.zeros(dc.cir.n) if k else None).x,
                                          float).tobytes())
                except nrsolver.NoConvergenceError as e:
                    out.append(type(e).__name__)
    finally:
        _dc_newton_c.ENABLED, nrsolver._within_rounding = kept
    return out, _paths.since(before), type(dc._get_linearsolver()._choice).__name__


@pytest.mark.parametrize('cells, solver', [(16, 'DenseSolver'), (130, 'StructuralKLU')])
def test_the_rounding_floor_is_the_python_one(cells, solver):
    """dcnewton.c's rounding floor (`floor_holds`: the rows of J^-1 by
    numpy's dgesv of J^T or klu_tsolve, a block at a time -- 130 failing
    unknowns, five blocks -- the flags read before the verdict): the second
    operating point's bytes the Python's, the C serving it -- and the
    Python's converged by the floor, without which neither converges at all
    (every rung of the ladder stalls at the rounding)."""
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    seen, orig = [], nrsolver._within_rounding

    def within(*a):
        seen.append(orig(*a))
        return seen[-1]
    py, _, kind = _floor_dc(cells, False, within)
    got, counts, _ = _floor_dc(cells, True)
    assert got == py and all(type(x) is bytes for x in py) and kind == solver
    assert sum(seen) >= 2 and counts.get('dcnewton_c:served', 0) == 1, (seen, counts)
    refused = _floor_dc(cells, False, lambda *a: False)[0]
    assert refused == ['NoConvergenceError'] * 2
