"""The stage predictor's fast path in C (`_tran_pred_c`, the C stepping loop
plan's stage 4, step a; 2026-10-11): bit for bit
`_StagePredictor._predict_fast` -- on drawn histories (times uniform, graded,
duplicated within the dedupe, tied, out of order, at or past the target;
states of every magnitude), degrees and targets, the prediction, its absence
and the general path the same; the hand-backs only where numpy would report
a flag; and alive: the clamp's factor moved is another prediction.
"""
import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _tran_pred_c
from pycircuit.circuit._tran_predictor import _GENERAL, _StagePredictor

pytestmark = pytest.mark.skipif(_tran_pred_c.driver() is None,
                                reason='the predictor unit did not build')


class _P(_StagePredictor):
    """The predictor alone: a history, the clamp, no periodic rows."""

    _periodic_rows = None

    class base_integrator:
        order = 2


def _same(hist, ttarget, deg):
    p = _P()
    p._pred_hist = [(t, x.copy()) for t, x in hist]
    py = p._predict_fast(ttarget, deg)
    c = _tran_pred_c.predict(p, ttarget, deg)
    return py, c


def _check(py, c):
    """The C's answer is the Python's, or a hand-back to the general path
    the fast path itself did not take (a flag numpy reports)."""
    if c is _GENERAL:
        return 'handed'
    if py is None or c is None:
        assert py is None and c is None, (py, c)
        return 'none'
    assert py is not _GENERAL
    assert py.tobytes() == c.tobytes(), (py, c)
    return 'served'


_mag = st.one_of(st.floats(-2.0, 2.0), st.floats(-1e-12, 1e-12),
                 st.sampled_from([0.0, -0.0, 1e-300, -1e-300, 1e150, 5e-324]))


@settings(deadline=None)
@given(data=st.data(), hk=st.integers(0, 8), n=st.integers(1, 7), deg=st.integers(1, 4),
       h=st.floats(1e-12, 1e-3), t0=st.sampled_from([0.0, 1e-9, 3.7e-6, 1.0, 1e3]))
def test_the_prediction_is_the_python_one_bit_for_bit(data, hk, n, deg, h, t0):
    kind = data.draw(st.sampled_from(['uniform', 'graded', 'dup', 'tie', 'order']))
    ts = [t0 + k * h for k in range(hk)]
    if kind == 'graded':
        ts = [t0 + h * (1.0 - 0.7 ** k) * 3.0 for k in range(hk)]
    elif kind == 'dup' and hk >= 2:
        k = data.draw(st.integers(1, hk - 1))
        ts[k] = ts[k - 1] * (1.0 + 1e-15) if ts[k - 1] else 1e-30
    elif kind == 'tie' and hk >= 2:
        k = data.draw(st.integers(1, hk - 1))
        ts[k] = ts[k - 1]
    elif kind == 'order' and hk >= 2:
        k = data.draw(st.integers(1, hk - 1))
        ts[k], ts[k - 1] = ts[k - 1], ts[k]
    xs = [np.array(data.draw(st.lists(_mag, min_size=n, max_size=n))) for _ in range(hk)]
    reach = data.draw(st.sampled_from([0.0, 0.5, 1.0, 2.0, -0.5]))
    last = ts[-1] if ts else t0
    ttarget = last + reach * h
    py, c = _same(list(zip(ts, xs)), ttarget, deg)
    _check(py, c)


def test_the_cases_the_fast_path_hands_on_are_handed_on():
    """At or past the target, out of order, tied: `_GENERAL` both ways."""
    x = np.array([1.0, 2.0])
    for hist, tt in (([(0.0, x), (1.0, x)], 1.0), ([(1.0, x), (0.5, x)], 2.0),
                     ([(0.0, x), (0.5, x), (0.5, x)], 1.0)):
        py, c = _same(hist, tt, 2)
        assert py is _GENERAL and c is _GENERAL


def test_a_uniform_ramp_is_served_and_exact():
    """A linear history: the prediction is the line (served, both ways)."""
    hist = [(k * 1e-9, np.array([k * 0.1, -k * 0.2, 1.0])) for k in range(6)]
    py, c = _same(hist, 6e-9, 2)
    assert _check(py, c) == 'served'
    np.testing.assert_allclose(c, [0.6, -1.2, 1.0], rtol=1e-12)


def test_a_moved_clamp_is_another_prediction(monkeypatch):
    """Alive: the clamp's factor read a hair off in the C -- the prediction
    of a history whose fit the clamp bites moves."""
    hist = [(k * 1e-9, np.array([math.sin(k), (-1.0) ** k])) for k in range(5)]
    py, c = _same(hist, 5e-9, 3)
    assert _check(py, c) == 'served'
    p = _P()
    p._pred_hist = list(hist)
    p.PRED_CLAMP = math.nextafter(_P.PRED_CLAMP, 0.0)
    moved = _tran_pred_c.predict(p, 5e-9, 3)
    assert moved is not _GENERAL and moved.tobytes() != c.tobytes()


def test_the_capacities_are_one_declaration():
    """pred.h's `PRED_HIST_MAX` and `PRED_FIT_MAX` are the reference's
    history capacity, the struct's arrays' lengths and the glue's limits --
    the struct's sizes are written as numbers (cffi lays it out), so this
    holds them to the names."""
    import cffi
    ffi = cffi.FFI()
    ffi.cdef(_tran_pred_c.PRED_CDEF)
    f = dict(ffi.typeof('pred_t').fields)
    F = _tran_pred_c.FIT_MAX
    assert _tran_pred_c.HIST_MAX == _StagePredictor.PRED_HIST_MAX
    assert f['tau'].type.length == f['w'].type.length == f['ipiv'].type.length == F
    assert f['V'].type.length == F * F
    assert _tran_pred_c.MAX_DEG == F - 1
