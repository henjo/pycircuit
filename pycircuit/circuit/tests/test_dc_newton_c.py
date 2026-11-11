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
    _tran_core,
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


def _ladder(cells):
    """`cells` resistors in a row from a 1 V source, each node to ground
    through a cubic conductance (a behavioural current, no limiter): a
    nonlinear operating point of cells + 2 unknowns."""
    from pycircuit.circuit.elements import VS, R, SubCircuit, gnd
    from pycircuit.circuit.elements_hdl import BSourceHdl
    c = SubCircuit()
    c['vs'] = VS('n0', gnd, v=1.0)
    for k in range(cells):
        c[f'r{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
        c[f'g{k}'] = BSourceHdl(f'n{k + 1}', gnd, f'n{k + 1}', gnd,
                                i_func=lambda v: 1e-4 * v + 1e-3 * v ** 3)
    return c


def _first_dc(cells, on):
    """`_ladder(cells)`'s operating point, the analysis's first solve, the
    C on or off: `(x bytes, the solver's class, its counts, the run's
    counts, the analysis)`."""
    from pycircuit.circuit.dcanalysis import DC
    kept = _dc_newton_c.ENABLED
    _dc_newton_c.ENABLED = on
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dc = DC(_ladder(cells), toolkit=circuit.numeric)
            before = _paths.snapshot()
            x = np.asarray(dc.solve().x, float).tobytes()
            counts = _paths.since(before)
    finally:
        _dc_newton_c.ENABLED = kept
    ch = dc._get_linearsolver()._choice
    klu = tuple(getattr(ch, k, None) for k in ('analyses', 'factors', 'refactors',
                                                'pivot_fallbacks', 'residual_fallbacks'))
    return x, type(ch).__name__, klu, counts, dc


def test_the_first_solve_is_the_cs_where_the_structure_chooses():
    """From the AutoSolver's structural floor (150) the circuit's structure
    chooses the StructuralKLU before the first matrix does
    (`_dc_newton_c.first_choice`), so the C serves the operating point's
    first solve too -- the factor its sequence makes with no numeric yet, as
    StructuralKLU's first solve -- x, the solver and its counts the
    Python's; below the floor the first solve is still the Python's
    (`dcnewton_c:linsolver`)."""
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    py = _first_dc(160, False)
    got = _first_dc(160, True)
    assert got[:3] == py[:3] and got[1] == 'StructuralKLU', (got[1:3], py[1:3])
    assert got[3].get('dcnewton_c:served', 0) >= 1, got[3]
    assert 'dcnewton_c:linsolver' not in got[3], got[3]
    small = _first_dc(40, True)
    assert small[1] == 'DenseSolver' and small[3].get('dcnewton_c:linsolver', 0) >= 1, small[3]


def test_a_klu_record_lays_out_klus_entries_alone(monkeypatch):
    """The record bound to KLU makes its active entries KLU's -- never the
    m^2 index arrays the dense solve alone lays out -- and keeps its R,
    |J| and C on small pages (`_tran_core.pattern_zeros`: resident where
    written; under huge pages one written entry made 2 MB resident, the
    whole m^2 of an IBM power grid's 25 850 unknowns)."""
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    laid, real = [], _dc_newton_c._Rec._dense_active

    def dense_active(self):
        laid.append(self.n)
        return real(self)
    monkeypatch.setattr(_dc_newton_c._Rec, '_dense_active', dense_active)
    *_, counts, dc = _first_dc(160, True)
    assert not laid, laid
    rec = dc.__dict__['_dc_newton_c']
    sk = dc._get_linearsolver()._choice
    assert rec.klu is sk and len(rec.act[0]) == sk.nnz < (rec.n - 1) ** 2, counts
    for k in ('R', 'absJ', 'Cs'):
        assert _tran_core._mapping(rec.buf[k]) is not None, k


def test_the_convergence_and_refactor_tests_are_decided_on_the_pattern(monkeypatch):
    """On KLU's pattern the C's convergence test (`aF < reltol * (|J| |x| +
    aF) + abstol`) and the refactor's residual test are decided by the
    sparse sums' bound (pcx_verdict.h) -- no dense product over the m^2
    |J| or R, whose untouched pages cost seconds at size -- and the
    operating point is the Python's, bit for bit; `Is`, which the rounding
    floor and the iteration limit read, made only there."""
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    py = _first_dc(160, False)
    got = _first_dc(160, True)
    assert got[:3] == py[:3]
    assert got[3].get('dcnewton_c:served', 0) >= 1, got[3]
    assert 'dcnewton_c:dense_conv' not in got[3] and 'dcnewton_c:dense_resid' not in got[3], got[3]
    rec = got[4].__dict__['_dc_newton_c']
    assert int(rec.s.is_made) == 0


def test_an_undecided_refactor_residual_is_the_dense_test_on_r(monkeypatch):
    """The refactor's residual verdict where its bound cannot decide (the
    solver's tolerance at the refactor's own residual, measured on the
    Python's solve): the dense test on R itself (pcx_kluseq.h's
    PCX_A_ON_PATTERN: R zero off KLU's pattern), counted -- the operating
    point the Python's, bit for bit.  Without the sparse verdict first
    (R as A built) nothing is counted."""
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    seen = []
    real = linearsolver.StructuralKLU._solve

    def solve(self, A, b, toolkit, kept):
        reused = self._numeric is not None
        x = real(self, A, b, toolkit, kept)
        if reused:
            arr, rhs = np.asarray(A), np.asarray(b, float)
            seen.append(float(np.abs(arr @ np.asarray(x) - rhs).max())
                        / max(float(np.abs(rhs).max()), 1e-300))
        return x
    with monkeypatch.context() as mp:
        mp.setattr(linearsolver.StructuralKLU, '_solve', solve)
        _first_dc(160, False)
    assert seen, 'no refactor to validate'
    monkeypatch.setattr(linearsolver.StructuralKLU, 'REFACTOR_RESIDUAL_TOL', seen[0])
    py = _first_dc(160, False)
    got = _first_dc(160, True)
    assert got[:3] == py[:3]
    assert got[3].get('dcnewton_c:served', 0) >= 1 and got[3].get('dcnewton_c:dense_resid', 0) >= 1, got[3]


def _base_solves(build, on, **kw):
    """`build()`'s operating point (the DC's parameters `kw`), the C on or
    off: every base solve's outcome in order -- x's bytes or the exception's
    type and text -- and the C's counts."""
    from pycircuit.circuit.dcanalysis import DC
    log = []
    real = nrsolver.StandardNewton.solve_system

    def solve_system(self, *a, **k):
        try:
            r = real(self, *a, **k)
        except Exception as e:
            log.append(('raise', type(e).__name__, str(e)))
            raise
        log.append(('ok', np.asarray(r[0]).tobytes(), r[1]))
        return r
    kept = _dc_newton_c.ENABLED, nrsolver.StandardNewton.solve_system
    _dc_newton_c.ENABLED = on
    nrsolver.StandardNewton.solve_system = solve_system
    before = _paths.snapshot()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dc = DC(build(), toolkit=circuit.numeric, **kw)
            try:
                dc.solve()
            except Exception as e:                             # noqa: BLE001
                log.append(('dc raised', type(e).__name__, str(e)))
    finally:
        _dc_newton_c.ENABLED, nrsolver.StandardNewton.solve_system = kept
    return log, _paths.since(before)


@pytest.mark.parametrize('cells, klu', [(160, True), (40, False)])
def test_the_iteration_limit_is_the_pythons(cells, klu):
    """An operating point whose Newton meets its iteration limit: on KLU's
    pattern, the residual's test decided there, the C makes `Is` for the
    exception's scale; on the dense record its test made it already.  Every
    base solve -- the limited ones' exceptions, their text, the ladders'
    rungs after them -- is the Python's."""
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    py, _ = _base_solves(lambda: _ladder(cells), False, maxiter=1)
    got, counts = _base_solves(lambda: _ladder(cells), True, maxiter=1)
    assert got == py
    assert counts.get('dcnewton_c:maxiter', 0) >= 1, counts
    assert any(e[0] == 'raise' for e in py), py
    assert (_first_dc(cells, True)[1] == 'StructuralKLU') is klu


def test_an_undecided_convergence_test_is_the_dense_one(monkeypatch):
    """A KLU operating point whose convergence test the sparse verdict
    cannot decide (`test_pcnr_c._stiff_ladder`: a row's sum |J| |x| past
    1e300) makes the dense test there, the Python's -- which holds, so the
    Newton stops where the Python's does (an undecided verdict read as a
    failure ran on to the iteration limit): every base solve the same with
    the C on and off, the undecided tests counted."""
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.tests.test_pcnr_c import _stiff_ladder
    ## (StructuralKLU from 40 unknowns: the C's KLU on the pattern)
    monkeypatch.setattr(DC, '_STRUCTURAL_KLU_MIN_N', 40)
    py, _ = _base_solves(lambda: _stiff_ladder(60), False)
    got, counts = _base_solves(lambda: _stiff_ladder(60), True)
    assert got == py
    assert counts.get('dcnewton_c:served', 0) >= 1, counts
    assert counts.get('dcnewton_c:dense_conv', 0) >= 1, counts
