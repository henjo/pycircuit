"""THE EARLY ABORT (`Transient.early_abort`, 2026-10-14; Andreas: "Build it
now"; doc/hdl_roadmap_260824.md, 132).

A multistep attempt is judged by its error estimate at its converged
point; a rejected one has paid for every Newton iteration first.  Its
estimate at the FIRST iterate -- the same chain on x1 and its charge, J^-1
by the seed's factorisation the Newton already holds -- told the 21 863
rejections of the 27 gold decks apart from the accepted attempts (over 2:
18 063 of them, 9 accepted steps), so an attempt over the threshold there
is abandoned and rejected at that error.  Measured: 9 % fewer Newton
iterations, ~6 % less solve time, no forced step, every gold column
passing as before.  The default threshold 1.25 since roadmap 134 (was 2):
5.8 % fewer Newton iterations again, 3.3-3.7 % less solve time.

What these hold: the rule is ONE rule on the three paths that make a
multistep attempt -- the Python loop and Newton (`_newton`'s check), the
Python loop on the C Newton (newton.c's, on `_tran_lte_c.early_context`)
and the C loop (loop.c) -- bit for bit, for the full Newton and the chord,
dense and sparse; a check that never aborts changes nothing; the loop's
rejection is the controller's law at the estimate; the arm keeps a full
attempt in the point's budget."""
import warnings
from types import SimpleNamespace

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
from pycircuit.circuit.stepcontroller import MIN_SHRINK_RATIO, IntegralController
from pycircuit.circuit.transient import Transient, _LMMSteps, _SteppingLoop

PATHS = {'python': (False, False, False), 'python loop, C Newton': (False, True, True),
         'C loop': (True, True, True)}


def _gm2(monkeypatch, path='C loop', sparse=False, **kw):
    """gm2 (CircuitSim90, MOS level 3; the suite's way) on `path`: `(the
    transient, its result, the `_paths` counts of the run)`."""
    deck = benchdata.spice_data('Netlists/CircuitSim90/MOS3/gm2.cir')
    if deck is None:
        pytest.skip('benchmark data not fetched (benchmarks/fetch_spice_suite.py)')
    loop, newton, lte = PATHS[path]
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', loop)
    monkeypatch.setattr(_tran_newton_c, 'ENABLED', newton)
    monkeypatch.setattr(_tran_lte_c, 'ENABLED', lte)
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(deck)
    tr, skw = imp.transient(**kw)
    if sparse:
        ## (KLU on the circuit's pattern however small: StructuralKLU)
        ls = tr._get_linearsolver()
        ls.min_n, ls.max_fill = 2, 1.0
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(**skw)
    return tr, res, _paths.since(before)


def _record(tr, res):
    st = tr.statistics
    return (np.asarray(res.x, float).tobytes(), st.accepted_steps, st.rejected_steps,
            st.newton_iterations, st.force_accepts, st.early_aborts)


def test_an_attempt_far_over_at_its_first_iterate_is_abandoned(monkeypatch):
    """gm2, default against `early_abort=0` (measured at the default 1.25:
    834 / 257 / 2061 Newton iterations, 205 abandoned; at 2, 835 / 263 / 2165
    and 148; off, 822 / 256 / 2376): every abandoned attempt is a rejection,
    none is forced, and the Newton work falls."""
    tr, _res, counts = _gm2(monkeypatch)
    off, _r, off_counts = _gm2(monkeypatch, early_abort=0)
    st, so = tr.statistics, off.statistics
    assert st.early_aborts == counts['early:aborts'] > 100, (st.early_aborts, counts)
    assert counts['early:checks'] > st.early_aborts
    assert so.early_aborts == 0 and 'early:checks' not in off_counts
    assert st.force_accepts == so.force_accepts == 0
    assert st.rejected_steps >= st.early_aborts
    assert st.newton_iterations < 0.95 * so.newton_iterations, (st.newton_iterations,
                                                                so.newton_iterations)


@pytest.mark.parametrize('chord', [False, True])
@pytest.mark.parametrize('sparse', [False, True])
def test_the_early_abort_is_one_rule_on_every_path(monkeypatch, chord, sparse):
    """The Python Newton's check (`_newton`), newton.c's on the Python loop
    and the C loop's: the same steps, rejections, Newton iterations and
    abandoned attempts, the waveform byte for byte -- the full Newton and
    the chord, dense and on KLU.  The full Newton on a sparse matrix holds
    no factors: no check there, on any path."""
    got = {p: _gm2(monkeypatch, p, sparse, chord_jacobian=chord) for p in PATHS}
    recs = {p: _record(tr, res) for p, (tr, res, _c) in got.items()}
    ref = recs['python']
    for p, rec in recs.items():
        assert rec == ref, (p, rec[1:], ref[1:])
    for p, (_tr, _res, counts) in got.items():
        assert counts.get('early:checks', 0) == got['python'][2].get('early:checks', 0), p
    if sparse and not chord:
        assert ref[5] == 0
    else:
        assert ref[5] > 50, ref[1:]


def test_a_check_that_never_aborts_changes_nothing(monkeypatch):
    """Every armed attempt's estimate made, none over the threshold: the run
    is the run with the abort off, byte for byte -- the check reads the
    attempt and the controller and writes neither."""
    tr, res, counts = _gm2(monkeypatch, early_abort=1e300)
    off, res_off, _c = _gm2(monkeypatch, early_abort=0)
    assert counts['early:checks'] > 100 and tr.statistics.early_aborts == 0
    assert _record(tr, res) == _record(off, res_off)


def test_an_abandoned_attempt_is_rejected_by_the_controllers_law():
    """`IntegralController.reject_early`: the shrink a judged rejection at
    the same error gets (`_shrink`: the first retry at a point at 1/p, a
    second at the order the two measured, MIN_SHRINK_RATIO-bounded) and
    the same retry record."""
    ctrl = IntegralController()
    x0 = np.zeros(2)
    s1 = SimpleNamespace(x_last=x0, h_curr=1e-9)
    h1 = ctrl.reject_early(s1, 8.0, 3.0)
    assert h1 == 1e-9 * max(MIN_SHRINK_RATIO, (0.9 ** 3.0 / 8.0) ** (1.0 / 3.0))
    assert ctrl._retry == (x0, 1e-9, 8.0)
    s2 = SimpleNamespace(x_last=x0, h_curr=h1)
    h2 = ctrl.reject_early(s2, 3.0, 3.0)
    q = min(max(np.log(8.0 / 3.0) / np.log(1e-9 / h1), 1.0), 3.0)
    assert h2 == pytest.approx(h1 * max(MIN_SHRINK_RATIO, (0.9 ** 3.0 / 3.0) ** (1.0 / q)),
                               rel=1e-15)
    assert ctrl.reject_early(SimpleNamespace(x_last=x0, h_curr=1e-9), 1e9, 3.0) == \
        1e-9 * MIN_SHRINK_RATIO


def test_the_arm_keeps_a_full_attempt_in_the_points_budget():
    """`_SteppingLoop.early_arm`: armed with history, while a rejection
    leaves a full attempt in the point's budget (`max_reject`,
    `max_retries`) and the step is above `minstep`; never on a fixed grid,
    with the abort off, under another controller or a band."""
    def loop(**over):
        ctrl = IntegralController()
        tr = SimpleNamespace(par=SimpleNamespace(early_abort=2.0, pcnr=False),
                             step_controller=ctrl)
        fam = _LMMSteps.__new__(_LMMSteps)
        fam.max_reject, fam.max_retries = 6, 10
        lp = SimpleNamespace(tr=tr, family=fam, fixed=False, X=[np.zeros(2), np.ones(2)],
                             rejects=0, point_retries=0, h=1e-9, minstep=1e-18, clamped=False,
                             run=SimpleNamespace(abstol=1e-12, max_step=1e-6))
        for k, v in over.items():
            if k in ('early_abort', 'pcnr'):
                setattr(tr.par, k, v)
            else:
                setattr(lp, k, v)
        return lp

    lp = loop()
    arm = _SteppingLoop.early_arm(lp)
    assert arm is not None and arm.th == 2.0 and arm.x_last is lp.X[-1] and arm.h == 1e-9
    assert _SteppingLoop.early_arm(loop(rejects=4)) is not None
    assert _SteppingLoop.early_arm(loop(rejects=5)) is None
    assert _SteppingLoop.early_arm(loop(point_retries=9)) is None
    assert _SteppingLoop.early_arm(loop(h=1e-18)) is None
    assert _SteppingLoop.early_arm(loop(fixed=True)) is None
    assert _SteppingLoop.early_arm(loop(X=[np.zeros(2)])) is None
    assert _SteppingLoop.early_arm(loop(early_abort=0)) is None
    assert _SteppingLoop.early_arm(loop(pcnr=True)) is None
    banded = loop()
    banded.tr.step_controller.set_lte_band(0.5, 1.0, None)
    assert _SteppingLoop.early_arm(banded) is None


def test_a_threshold_below_one_is_refused():
    """Below 1 the check would abandon attempts the error test passes."""
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import VS, C, R
    c = SubCircuit()
    c['v'] = VS(1, gnd, v=1.0)
    c['r'] = R(1, 2, r=1e3)
    c['c'] = C(2, gnd, c=1e-9)
    tr = Transient(c, toolkit=circuit.numeric, early_abort=0.5)
    with pytest.raises(ValueError, match='early_abort'):
        tr.solve(tend=1e-6, timestep=1e-8)
