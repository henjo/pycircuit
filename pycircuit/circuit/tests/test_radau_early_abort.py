"""THE RADAU EARLY ABORT (`Transient.radau_early_abort`, 2026-10-21;
Andreas: "build prototype", then "accept"; doc/hdl_roadmap_260824.md, 159).

An adaptive Radau attempt is judged by its embedded estimate once its
transform Newton has converged -- and one that does not converge runs the
transform to `maxiter` and the dense coupled Newton before the judge says
no.  Sized over Radau's 20 decks (plan, 'AN EARLY ABORT SIZED'): 736
attempts ran the transform to `maxiter` and fell back, 31 % of the time,
614 of them rejected after; the estimate at the transform's first iterate
is the final one to a few % (median ratio 1.00).  So an armed attempt runs
its first iteration alone, and the estimate there, normalised as the step
is judged, over the threshold abandons it as a rejection sized by the stage
family's retry law.

What these hold: the check is one rule on the transform's C and its Python
loop, bit for bit; a check that never aborts changes nothing (the first
iteration and the rest are the iterates one call makes); the abandoned
attempt is retried by `judge`'s law at the worst unknown; the arm keeps a
full attempt in the point's budget."""
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from pycircuit._testing import benchdata
from pycircuit.circuit import _paths, _tran_radau_tc, circuit, spice_import
from pycircuit.circuit import integrator as I
from pycircuit.circuit._tran_radau import _RadauAbort, _RadauArm
from pycircuit.circuit.stepcontroller import MIN_SHRINK_RATIO
from pycircuit.circuit.transient import _StageSteps, _SteppingLoop


def _latch(monkeypatch, c=True, **kw):
    """The MCNC BJT latch under Radau's default (adaptive, the transform),
    its loop in C or in Python (`c`): `(the transient, its result, the
    `_paths` counts of the run)`."""
    deck = benchdata.spice_data('Netlists/MCNC_BJT_LATCH/latch.cir')
    if deck is None:
        pytest.skip('benchmark data not fetched (benchmarks/fetch_spice_suite.py)')
    monkeypatch.setattr(_tran_radau_tc, 'ENABLED', c)
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(deck)
    tr, skw = imp.transient()
    tr.par.integrator = I.RadauIIA3Integrator()
    for k, v in kw.items():
        setattr(tr.par, k, v)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(**skw)
    return tr, res, _paths.since(before)


def _record(tr, res):
    st = tr.statistics
    return (np.asarray(res.x, float).tobytes(), st.accepted_steps, st.rejected_steps,
            st.force_accepts, st.early_aborts)


def test_a_radau_attempt_far_over_at_its_first_iterate_is_abandoned(monkeypatch):
    """The latch, the default (2) against `radau_early_abort=0` (measured:
    56 attempts abandoned, the transform's iterations 5 941 -> 3 288, its
    dense fallbacks 45 -> 19): every abandoned attempt a rejection, none
    forced, the fallbacks and the transform's iterations fewer."""
    tr, _res, counts = _latch(monkeypatch)
    off, _r, off_counts = _latch(monkeypatch, radau_early_abort=0)
    st, so = tr.statistics, off.statistics
    assert st.early_aborts == counts['radau.early:aborts'] > 20, (st.early_aborts, counts)
    assert counts['radau.early:checks'] > st.early_aborts
    assert so.early_aborts == 0 and 'radau.early:checks' not in off_counts
    assert st.force_accepts == so.force_accepts == 0
    assert st.rejected_steps >= st.early_aborts
    assert counts['radau.transform:fallback'] < off_counts['radau.transform:fallback']
    assert counts['radau_tc:iterations'] < 0.8 * off_counts['radau_tc:iterations']


def test_a_check_that_never_aborts_changes_nothing(monkeypatch):
    """Every armed attempt's first iteration made alone and its estimate
    checked, none over the threshold: the run is the run with the abort
    off, byte for byte -- the C's first iteration and its resumed rest are
    the iterates one call makes, and the check writes nothing the run
    reads."""
    tr, res, counts = _latch(monkeypatch, radau_early_abort=1e300)
    off, res_off, _c = _latch(monkeypatch, radau_early_abort=0)
    assert counts['radau.early:checks'] > 100 and tr.statistics.early_aborts == 0
    assert counts['radau_tc:probe'] > 100
    assert _record(tr, res) == _record(off, res_off)


def test_the_check_is_one_rule_on_the_c_and_the_python_loop(monkeypatch):
    """The transform's loop in C and in Python (`PYCIRCUIT_RADAU_TC`): the
    same attempts abandoned, the same steps and rejections, the waveform
    byte for byte."""
    tr, res, counts = _latch(monkeypatch)
    py, res_py, py_counts = _latch(monkeypatch, c=False)
    assert counts.get('radau_tc:probe', 0) > 100 and 'radau_tc:probe' not in py_counts
    assert _record(tr, res) == _record(py, res_py)
    assert counts['radau.early:checks'] == py_counts['radau.early:checks']


def _family(order='measured'):
    fam = _StageSteps.__new__(_StageSteps)
    fam.retry_order = order
    fam._retry = fam._q_carried = None
    fam.tr = SimpleNamespace(base_integrator=SimpleNamespace(EMBEDDED_ORDER=3))
    return fam


def test_an_abandoned_attempt_is_retried_by_judges_law():
    """`_StageSteps.reject_early`: the retry `judge` sizes at the same worst
    unknown (`_retry_after`: a point's first retry at the order 1/p, a
    second at the order the two tries measured, MIN_SHRINK_RATIO-bounded)
    and the same retry record."""
    fam = _family()
    x0 = np.zeros(2)
    p = 4
    h1 = fam.reject_early(_RadauAbort(3.0, 8.0, x0, 1e-9))
    assert h1 == 1e-9 * max(MIN_SHRINK_RATIO, (0.9 ** p / 8.0) ** (1.0 / p))
    assert fam._retry == (x0, 1e-9, 8.0)
    h2 = fam.reject_early(_RadauAbort(2.5, 3.0, x0, h1))
    q = min(max(np.log(8.0 / 3.0) / np.log(1e-9 / h1), 1.0), p)
    assert h2 == pytest.approx(h1 * max(MIN_SHRINK_RATIO, (0.9 ** p / 3.0) ** (1.0 / q)),
                               rel=1e-15)
    assert _family().reject_early(_RadauAbort(1e9, 1e9, x0, 1e-9)) == 1e-9 * MIN_SHRINK_RATIO


def test_the_arm_keeps_a_full_attempt_in_the_points_budget():
    """`_SteppingLoop.radau_arm`: armed on an adaptive stage step whose
    estimate is made, while a rejection would leave a full attempt in the
    point's budget, above `minstep`; never with the threshold 0."""
    fam = SimpleNamespace(max_reject=12, max_retries=12)
    x = np.zeros(2)

    def loop(**kw):
        d = {'tr': SimpleNamespace(par=SimpleNamespace(radau_early_abort=2.0),
                                   _rk_want_est=True),
             'family': fam, 'fixed': False, 'rejects': 0, 'point_retries': 0, 'h': 1e-9,
             'minstep': 1e-15, 'X': [x]}
        d.update(kw)
        return SimpleNamespace(**d)

    arm = _SteppingLoop.radau_arm(loop())
    assert type(arm) is _RadauArm and arm.th == 2.0 and arm.fam is fam
    assert arm.x_last is x and arm.h == 1e-9
    assert _SteppingLoop.radau_arm(loop(rejects=10)) is not None
    assert _SteppingLoop.radau_arm(loop(rejects=11)) is None
    assert _SteppingLoop.radau_arm(loop(point_retries=11)) is None
    assert _SteppingLoop.radau_arm(loop(h=1e-15)) is None
    assert _SteppingLoop.radau_arm(loop(fixed=True)) is None
    off = loop()
    off.tr.par.radau_early_abort = 0
    assert _SteppingLoop.radau_arm(off) is None
    noest = loop()
    noest.tr._rk_want_est = False
    assert _SteppingLoop.radau_arm(noest) is None


def test_a_threshold_under_one_is_refused(monkeypatch):
    """Under 1 the check would abandon attempts the judge accepts."""
    with pytest.raises(ValueError, match='radau_early_abort'):
        _latch(monkeypatch, radau_early_abort=0.5)
