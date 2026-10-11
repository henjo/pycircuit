"""THE STEP CONTROLLER'S LAW IN C (the C stepping loop plan's stage 5, step a;
2026-10-11): `IntegralController._evaluate` past its error test -- the
default band's accept or reject, the next step (`0.9 ** p`, the shrink with
`retry_order`'s exponent from the last rejection at the same point, the
grow, the cap), `_judged`'s retry state -- as `csrc/ctrl.c`, which the
adaptive loop in C calls (`csrc/loop.c`).  Bit for bit the Python's: its float
operations in their order, `pow` and `log` the same libm's, its `min` and
`max` as Python's builtins pick.  "The same point" is the object the
controller kept (`is`) in the Python, an index of the accepted points here.

`evaluate(...)` is the twin's entry for its tests.  `ENABLED` (env
`PYCIRCUIT_CTRL_C=0`) is read on every call.
"""
import os

from pycircuit.circuit import _cbuild, _paths

_PC = _paths.COUNTS

ENABLED = os.environ.get('PYCIRCUIT_CTRL_C', '1') != '0'
STATUS = 'not loaded'

CTRL_C = _cbuild.csrc('ctrl.c')
CTRL_CDEF = _cbuild.csrc('ctrl.h')
_CODES = _cbuild.constants(CTRL_CDEF)
MODES = {'assumed': _CODES['CTRL_ASSUMED'], 'measured': _CODES['CTRL_MEASURED'],
         'carried': _CODES['CTRL_CARRIED']}

_driver = None


def driver():
    """`(ffi, cfn)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(CTRL_C, CTRL_CDEF)
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn)
        STATUS = 'c'
    return _driver or None


def new_state(ffi, mode):
    """A controller struct for `retry_order` `mode`, the module's constants
    set."""
    from pycircuit.circuit import stepcontroller as sc
    s = ffi.new('ctrl_t *')
    s.mode = MODES[mode]
    s.min_shrink, s.max_growth = sc.MIN_SHRINK_RATIO, sc.MAX_GROWTH_RATIO
    s.order_min = sc.RETRY_ORDER_MIN
    return s


def evaluate(err, p, h, max_step, point, mode, retry=None, q_carried=None):
    """The C law on one judged attempt: `(accepted, h_next, retry,
    q_carried)` -- `retry` the kept rejection `(point, h, err)` or None, as
    `_judged` leaves `_retry` (its point an index)."""
    ffi, cfn = driver()
    s = new_state(ffi, mode)
    s.err, s.p, s.h, s.max_step, s.point = err, p, h, max_step, point
    s.have_retry = retry is not None
    if retry is not None:
        s.retry_point, s.retry_h, s.retry_err = retry
    s.have_q = q_carried is not None
    s.q_carried = q_carried if q_carried is not None else 0.0
    cfn(s)
    return (bool(s.accepted), float(s.h_next),
            (int(s.retry_point), float(s.retry_h), float(s.retry_err)) if s.have_retry else None,
            float(s.q_carried) if s.have_q else None)
