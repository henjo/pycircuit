"""The step controller's law in C (`_tran_ctrl_c`, the C stepping loop
plan's stage 5, step a; 2026-10-11): bit for bit `IntegralController.
_evaluate` past its error test -- on drawn errors, orders, steps, caps and
retry states (the last rejection at the same point or another, larger or
not, above the target or not; a carried order), under every `retry_order`:
the verdict, the next step, `_retry` and `_q_carried` the same; and alive: a
retry at another point is another step.
"""
import math

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _tran_ctrl_c
from pycircuit.circuit import stepcontroller as sc

pytestmark = pytest.mark.skipif(_tran_ctrl_c.driver() is None,
                                reason='the controller unit did not build')


class _S:
    """The attempt as the law reads it."""

    def __init__(self, h, max_step, x_last):
        self.h_curr, self.max_step, self.x_last = h, max_step, x_last
        self.no_history, self.h_clamped = False, False


def _python(err, p, h, max_step, same, mode, retry, q_carried):
    c = sc.IntegralController()
    c.set_retry_order(mode)
    here, there = object(), object()
    c._retry = None if retry is None else ((here if same else there), retry[0], retry[1])
    c._q_carried = q_carried
    c._max_error = lambda s: (err, p)
    accepted, h_next = c._evaluate(_S(h, max_step, here))
    r = c._retry
    return (accepted, h_next, None if r is None else (r[0] is here, r[1], r[2]), c._q_carried)


def _c(err, p, h, max_step, same, mode, retry, q_carried):
    accepted, h_next, r, q = _tran_ctrl_c.evaluate(
        err, p, h, max_step, 7, mode,
        None if retry is None else ((7 if same else 3), retry[0], retry[1]), q_carried)
    return (accepted, h_next, None if r is None else (r[0] == 7, r[1], r[2]), q)


def _hex(t):
    return tuple(float(v).hex() if isinstance(v, float) else (_hex(v) if isinstance(v, tuple) else v)
                 for v in t)


_err = st.one_of(st.floats(0.0, 50.0), st.sampled_from([0.0, 1e-13, 1.0, math.nextafter(1.0, 2.0),
                                                       0.729, 1e300, 5e-324]))


@settings(deadline=None)
@given(err=_err, p=st.sampled_from([2.0, 3.0]), h=st.floats(1e-15, 1e-3),
       cap=st.floats(0.1, 10.0), same=st.booleans(), mode=st.sampled_from(sc.RETRY_ORDER_MODES),
       retry=st.one_of(st.none(), st.tuples(st.floats(0.5, 20.0), st.floats(0.0, 100.0))),
       q_carried=st.one_of(st.none(), st.floats(1.0, 3.0)))
def test_the_law_is_the_python_one_bit_for_bit(err, p, h, cap, same, mode, retry, q_carried):
    """The kept rejection's step a multiple of this one (larger or not);
    the cap a multiple of the step."""
    retry = None if retry is None else (h * retry[0], retry[1])
    a = _python(err, p, h, h * cap, same, mode, retry, q_carried)
    b = _c(err, p, h, h * cap, same, mode, retry, q_carried)
    assert _hex(a) == _hex(b), (a, b)


def test_a_second_retry_at_one_point_measures_its_order():
    """Rejected at h, then again at h/2 at the same point: the order the
    two measured is the shrink's (and is carried); at another point it is
    not."""
    first = _c(4.0, 3.0, 1e-9, 1e-6, True, 'measured', None, None)
    assert first[0] is False and first[2] is not None
    second = _c(3.0, 3.0, 0.5e-9, 1e-6, True, 'measured', (1e-9, 4.0), None)
    elsewhere = _c(3.0, 3.0, 0.5e-9, 1e-6, False, 'measured', (1e-9, 4.0), None)
    assert second[3] is not None and elsewhere[3] is None
    assert second[1] != elsewhere[1]
    assert _hex(second) == _hex(_python(3.0, 3.0, 0.5e-9, 1e-6, True, 'measured', (1e-9, 4.0), None))
