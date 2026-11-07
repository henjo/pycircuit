"""The chord's fallback in the C Newton (newton.h's THE CHORD'S FALLBACK,
`_tran_newton_c.FALLBACK`; roadmap 268).  `ChordNewton.solve_system` falls
back to the full Newton from the same seed where its chord breaks or meets
its limit; the C Newton handed those attempts back, and the Python repeated
the whole attempt on dense state -- every Gear-2 hand-back measured was one
(a 4 000-stage chain's 1 us run: 12, about 70 % of its time).  Now the C
runs the full Newton there: its first evaluation the chord's first C (what
`_C_lookup` finds), no early check past the chord's, its iterations added,
`chord_fallbacks` counted, a sparse run's numeric refactored as the
Python's full Newton refactors it -- every result the Python's."""
import warnings

import numpy as np
import pytest

from pycircuit._testing import benchdata
from pycircuit.circuit import (
    _paths,
    _tran_loop_c,
    _tran_lte_c,
    _tran_newton_c,
    circuit,
    spice_import,
)
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.transient import Transient

PATHS = {'python': (False, False, False), 'python loop, C Newton': (False, True, True),
         'C loop': (True, True, True)}


def _deck(monkeypatch, name, path, fallback=True, dense=False):
    """todd3 / slowlatch (CircuitSim90, the suite's way) on `path`: the run's
    bytes and statistics, its linear solver's refactors, its counts."""
    data = benchdata.spice_data(f'Netlists/CircuitSim90/{name}.cir')
    if data is None:
        pytest.skip('benchmark data not fetched (benchmarks/fetch_spice_suite.py)')
    loop, newton, lte = PATHS[path]
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', loop)
    monkeypatch.setattr(_tran_newton_c, 'ENABLED', newton)
    monkeypatch.setattr(_tran_lte_c, 'ENABLED', lte)
    monkeypatch.setattr(_tran_newton_c, 'FALLBACK', fallback)
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(data)
    tr, kw = imp.transient()
    if dense:
        ## (numpy's dense solve whatever the size: the full Newton's dgesv)
        ls = tr._get_linearsolver()
        ls.min_n = ls.structural_min_n = 10 ** 9
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(**kw)
    st = tr.statistics
    ch = tr._get_linearsolver()._choice
    run = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps, st.newton_iterations, st.chord_fallbacks,
           st.early_aborts)
    klu = (type(ch).__name__, getattr(ch, 'refactors', None), getattr(ch, 'factors', None))
    return run, klu, _paths.since(before)


@pytest.mark.parametrize('dense', [False, True], ids=['sparse', 'dense'])
@pytest.mark.parametrize('name', ['MOS3/todd3', 'MOS2/slowlatch'])
def test_the_c_falls_back_as_the_python_does(monkeypatch, name, dense):
    """todd3 (chord breaks and limits) and slowlatch (a limit), on KLU and on
    numpy's dense solve: the C loop's run with the fallback is the Python's
    -- bytes, steps, rejections, Newton iterations, chord fallbacks, early
    aborts, the KLU's refactors and factors -- and the C loop's without it
    (the attempts handed back); no chord hands an attempt back."""
    want, wklu, wcounts = _deck(monkeypatch, name, 'python', dense=dense)
    got, gklu, counts = _deck(monkeypatch, name, 'C loop', dense=dense)
    off, oklu, ocounts = _deck(monkeypatch, name, 'C loop', fallback=False, dense=dense)
    assert want[5] > 0 and wklu[0] == ('DenseSolver' if dense else 'StructuralKLU')
    assert got == want == off
    assert gklu == wklu == oklu
    assert counts['newton_c:fallbacks'] == want[5]
    assert counts.get('early:checks', 0) == wcounts.get('early:checks', 0) > 0
    assert counts.get('loop_c:bail:newton', 0) == 0
    assert not any(k.startswith('newton_c:bail:') for k in counts)
    assert ocounts['loop_c:bail:newton'] > 0 and 'newton_c:fallbacks' not in ocounts


def test_the_python_loops_c_newton_falls_back_too(monkeypatch):
    """The Python loop on the C Newton (`_tran_newton_c.solve`): the same
    run, the fallbacks served, none handed back."""
    want, wklu, wcounts = _deck(monkeypatch, 'MOS3/todd3', 'python')
    got, gklu, counts = _deck(monkeypatch, 'MOS3/todd3', 'python loop, C Newton')
    assert got == want and gklu == wklu
    assert counts.get('early:checks', 0) == wcounts.get('early:checks', 0) > 0
    assert counts['newton_c:fallbacks'] == want[5] > 0
    assert not any(k.startswith('newton_c:bail:') for k in counts)


def _active(stages=260, vto=0.5, kp=2e-4):
    """A chain whose stage gain is above one -- a pulse travels it -- with 264
    unknowns: KLU and the chord ('auto')."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=5e-9, tr=2e-9, tf=2e-9, pw=4e-8, per=1e-7)
    for k in range(stages):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd, vto=vto,
                                     kp=kp)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
    return c


def _limited(monkeypatch, path, maxiter=4, armed=False, plant=None, mos=(), **par):
    """The active chain with few iterations, on `path`: chord and full Newton
    both stop at the limit on the steps over the input's edge.  `mos` its
    devices' parameters where others, `par` the transient's, `plant` called on
    each C loop's record as it is loaded."""
    loop, newton, lte = PATHS[path]
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', loop)
    monkeypatch.setattr(_tran_newton_c, 'ENABLED', newton)
    monkeypatch.setattr(_tran_lte_c, 'ENABLED', lte)
    if plant is not None:
        load = _tran_loop_c._load

        def planted(loop_, tr_, td, run_, *a, **k):
            r = load(loop_, tr_, td, run_, *a, **k)
            plant(run_)
            return r
        monkeypatch.setattr(_tran_loop_c, '_load', planted)
    tr = Transient(_active(**dict(mos)), toolkit=circuit.numeric, **par)
    tr.par.maxiter = maxiter
    if armed:
        tr._damped_last_resort = True
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(tend=6e-8, timestep=1e-9)
    st = tr.statistics
    ch = tr._get_linearsolver()._choice
    run = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps, st.newton_iterations, st.chord_fallbacks,
           st.early_aborts)
    klu = (type(ch).__name__, ch.refactors, ch.factors, ch.pivot_fallbacks, ch.residual_fallbacks,
           ch.chord_refactors, ch.chord_factors, ch.chord_pivot_fallbacks,
           ch.chord_residual_fallbacks)
    return run, klu, _paths.since(before)


def test_the_fallback_s_limit_is_taken_without_a_repeat(monkeypatch):
    """A chord whose full Newton fallback meets the limit: the C Newton raises
    its exception and the C loop rejects the attempt, where the Python's repeat
    could only reproduce it -- nothing handed back, the run the Python's:
    bytes, statistics, the KLU's every count, the chord and early-check counts."""
    want, wklu, wc = _limited(monkeypatch, 'python')
    assert wklu[0] == 'StructuralKLU'
    for path in ('python loop, C Newton', 'C loop'):
        got, gklu, c = _limited(monkeypatch, path)
        assert got == want and gklu == wklu, path
        assert c['newton_c:maxiter'] > 0, c
        assert not any(k.startswith('newton_c:bail:') for k in c), c
        assert c.get('loop_c:bail:newton', 0) == 0, c
        assert c.get('chord:chosen') == wc.get('chord:chosen') > 0
        assert c.get('early:checks', 0) == wc.get('early:checks', 0) > 0
    assert c['loop_c:newton_rejects'] > 0, c


def test_an_armed_retry_still_hands_the_limit_back(monkeypatch):
    """With the damped retry armed the limit is not the end of the attempt:
    the C hands it back, its numeric and counts as they were, and the Python's
    retry runs -- the run the Python's."""
    want, wklu, _wc = _limited(monkeypatch, 'python', armed=True)
    got, gklu, c = _limited(monkeypatch, 'python loop, C Newton', armed=True)
    assert got == want and gklu == wklu
    assert c.get('newton_c:bail:maxiter', 0) > 0 and not c.get('newton_c:maxiter'), c


#: the first fallback meets its limit here, at 4 iterations; on the active
#: chain at 6 it converges -- neither run's numeric made before it
_GAIN_BELOW_ONE = (('vto', 0.0), ('kp', 2e-5))
_FIRST = {'converges': (6, ()), 'limit': (4, _GAIN_BELOW_ONE)}


def _seeded(monkeypatch, path, case, **k):
    maxiter, mos = _FIRST[case]
    return _limited(monkeypatch, path, maxiter, mos=mos, lte_jacobian='seed', **k)


@pytest.mark.parametrize('case', _FIRST)
def test_the_fallback_makes_the_run_s_first_numeric(monkeypatch, case):
    """The chord's first fallback with no run numeric yet (the seed's error
    test): its first solve factors one in the C, as `StructuralKLU._solve`, and
    a full error test after it reads it -- nothing handed back, the run the
    Python's: bytes, statistics, the KLU's every count (a factor, then
    refactors) -- whether that fallback converges or meets its limit."""
    want, wklu, _c = _seeded(monkeypatch, 'python', case)
    assert wklu[0] == 'StructuralKLU' and wklu[1] > 0 and wklu[2] == 1
    assert (want[5] > 0) == (case == 'converges')
    for path in ('python loop, C Newton', 'C loop'):
        run, klu, c = _seeded(monkeypatch, path, case)
        assert run == want and klu == wklu, path
        assert c['newton_c:first_numeric'] == 1, (path, c)
        assert not any(k.startswith(('newton_c:bail:', 'loop_c:bail:')) for k in c), (path, c)
    assert c['loop_c:attempts'] > 20


def test_a_hand_back_frees_the_first_numeric(monkeypatch):
    """A fallback handed back after making the numeric (the damped retry armed:
    the limit is not the attempt's end): the C frees it, and the Python's
    repeat factors its own from the same values -- a factor, not a refactor of
    the C's -- the run and its counts the Python's."""
    want, wklu, _c = _seeded(monkeypatch, 'python', 'limit', armed=True)
    run, klu, c = _seeded(monkeypatch, 'python loop, C Newton', 'limit', armed=True)
    assert run == want and klu == wklu and wklu[2] == 1
    assert c['newton_c:bail:maxiter'] > 0 and 'newton_c:first_numeric' not in c, c


def test_an_attempt_the_loop_hands_back_frees_its_first_numeric(monkeypatch):
    """A fallback that converged and made the numeric, then an error test that
    hands the attempt back (its validation failing, planted on the C loop's
    test alone): the attempt is the Python's, so the numeric goes with it and
    the repeat factors its own -- the run and its counts the Python's."""
    def failing(run_):
        if run_.sk is not None:
            run_.lte.s.resid_tol = 0.0
    want, wklu, _c = _seeded(monkeypatch, 'python', 'converges')
    run, klu, c = _seeded(monkeypatch, 'C loop', 'converges', plant=failing)
    assert run == want and klu == wklu and wklu[2] == 1
    assert c['loop_c:bail:lte'] > 0, c
