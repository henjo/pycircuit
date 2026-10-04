"""The stage predictor's multistep fast path (`_tran_predictor.PRED_FAST`,
speed round 8, stage 2b): bit for bit the general path."""
import numpy as np
from hypothesis import given
from hypothesis import strategies as st

from pycircuit.circuit import _paths, _tran_predictor, circuit
from pycircuit.circuit.elements import C, R, SubCircuit, VSin, gnd
from pycircuit.circuit.transient import Transient


def _rc():
    c = SubCircuit()
    c['vs'] = VSin('in', gnd, va=1.0, freq=1e3)
    c['R'] = R('in', 'out', r=1e3)
    c['C'] = C('out', gnd, c=1e-7)
    return c


def _both(tr, ttarget, deg=None):
    """The prediction with the fast path on and off (the same instance,
    each with its own weights: the memo emptied before each)."""
    out = []
    for on in (True, False):
        old = _tran_predictor.PRED_FAST
        _tran_predictor.PRED_FAST = on
        try:
            tr.__dict__.pop('_pred_fast', None)
            tr.__dict__.pop('_pred_wmemo', None)
            out.append(tr._predict_state(ttarget, deg=deg))
        finally:
            _tran_predictor.PRED_FAST = old
    return out


def _same(a, b):
    if a is None or b is None:
        return a is None and b is None
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


_finite = st.floats(-1e6, 1e6, allow_nan=False, allow_infinity=False)


@given(times=st.lists(st.floats(0.0, 1e-3, allow_nan=False), min_size=2, max_size=8, unique=True),
       ahead=st.floats(1e-12, 1e-3), deg=st.sampled_from([1, 2, 3, 4]),
       data=st.data())
def test_the_fast_path_is_the_general_path_bit_for_bit(times, ahead, deg, data):
    """Any history behind the target (times close enough to merge in the
    1e-13 dedupe included), any degree: the same bytes, or both None."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    n = 4
    times = sorted(times)
    tr._pred_hist = [(float(t), np.array(data.draw(st.lists(_finite, min_size=n, max_size=n))))
                     for t in times]
    a, b = _both(tr, times[-1] + ahead, deg)
    assert _same(a, b)


def test_near_duplicate_times_and_a_wrapping_row():
    """Times inside the general path's 1e-13 dedupe (but outside the
    history's own 1e-14) and a periodic row: still the same bytes."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    t0 = 1e-3
    tr._pred_hist = [(t0 - 3e-6, np.array([0.1, -0.2, 0.3])),
                     (t0 - 2e-6, np.array([0.2, -0.1, 0.25])),
                     (t0 - 1e-6, np.array([0.3, 0.0, 0.2])),
                     (t0 - 1e-6 + 5e-17, np.array([0.31, 0.01, 0.21]))]
    tr._periodic_rows = [(1, 1.0, 0.0)]
    for deg in (1, 2, 3):
        a, b = _both(tr, t0, deg)
        assert _same(a, b)


def test_the_fast_path_serves_the_multistep_step():
    """A gear transient's predictions take it: the count says so."""
    before = _paths.snapshot()
    Transient(_rc(), toolkit=circuit.numeric).solve(tend=2e-4, timestep=1e-5,
                                                    fixed_timestep=True)
    d = _paths.since(before)
    assert d.get('pred:fast', 0) > 10 and d.get('pred:general', 0) == 0


def test_an_edited_history_takes_the_general_path():
    """A history not in ascending order (a test's edit): the general path,
    counted, and the same prediction."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    tr._pred_hist = [(2e-6, np.array([0.2, 0.1, 0.0])), (1e-6, np.array([0.1, 0.0, 0.0])),
                     (3e-6, np.array([0.3, 0.2, 0.1]))]
    before = _paths.snapshot()
    a, b = _both(tr, 4e-6, 2)
    assert _same(a, b)
    assert _paths.since(before).get('pred:general', 0) == 1


def test_a_patched_fit_takes_the_general_path(monkeypatch):
    """A `_fit` patched on the class is a caller that expects its calls:
    the fast path steps aside, counted, and the patch is reached."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    tr._pred_hist = [(1e-6, np.array([0.1, 0.0, 0.0])), (2e-6, np.array([0.2, 0.1, 0.0])),
                     (3e-6, np.array([0.3, 0.2, 0.1]))]
    seen = []
    real = Transient._fit

    def counting(self, sub, tat):
        seen.append(len(sub))
        return real(self, sub, tat)
    monkeypatch.setattr(Transient, '_fit', counting)
    before = _paths.snapshot()
    tr._predict_state(4e-6, deg=2)
    assert seen == [3] and _paths.since(before).get('pred:patched', 0) == 1
