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
