"""EXACT SCALAR FAST PATHS for the numpy twins of a chained model (the speed
plan round 2, stage C; 2026-10-02).

The generated chain calls a few numpy functions on SCALARS hundreds of times
per evaluation -- PSP's `G` twin: 704 `numpy.where`, 263 `_step`, 281
`maxc` / `minc`, 82 `numpy.maximum`, ~357 comparisons -- and on scalars
numpy's dispatch is the cost: `numpy.where` 1.07 us (and it returns a 0-d
ARRAY, whose next arithmetic costs 450 ns instead of 34), `1.0 *
numpy.bool_` (`_step`) 840 ns, `numpy.maximum` 610 ns, `numpy.less` ~700
ns, against ~80 ns for the same decision taken in Python.

Each helper here answers EXACTLY what its numpy counterpart answers -- the
same type, the same bits, the same raise-mode behaviour -- on the scalar
types it takes a fast branch for, and hands everything else to the
counterpart itself:

* the decision is a Python comparison of the same operands, and the result
  is one of the operands (or a constant) as the numpy SCALAR numpy would
  return -- never a Python float or bool, which would turn a later `x / 0`
  into a `ZeroDivisionError` where numpy gives inf;
* ties and NaN, where numpy's choice follows the hardware (the operand of
  a +-0 tie, which NaN), always go to numpy;
* only `float`, `numpy.float64` (and, in comparisons, `int` within 2**53,
  where Python's exact int/float comparison and numpy's conversion agree)
  are fast; 0-d arrays, arrays, int64, sympy objects, tracers go through.

`_fwhere` returns a numpy SCALAR where `numpy.where` returns a 0-d ARRAY.
They compute the same values in every operation, and the same bits in all
but two cases, both measured (2026-10-02) and both closed:

* `**`: numpy squares / roots / reciprocates a 0-d array operand (exponent
  2, 1/2, -1) and raises a scalar through glibc `pow` -- so
  `_hdl_cse._fast_rewrite` keeps `numpy.where` for every call whose value
  reaches a `**` operand (base or exponent, through plain renames), as the
  C printer prints `_sq` for exactly those (`array_syms`);
* two NaNs meeting: which one an operation returns follows the operand
  order, which numpy's 0-d loop and its scalar arithmetic order differently
  -- so `_fwhere` is fast only for a FINITE chosen value.

The shared helpers (`hdl._maxc_numpy` and the rest) are NOT replaced: the
JAX twins, the eager path's lambdify, `CY` on frequency arrays and sympy's
evalf call them too.  These are new names, used by the rewritten twin text
only.
"""
from math import isfinite as _isfinite
from math import isnan as _isnan

import numpy as np

_f64 = np.float64
_TRUE, _FALSE = np.True_, np.False_
_ONE, _ZERO = np.float64(1.0), np.float64(0.0)
_INT_EXACT = 2 ** 53


def _num(v):
    """`v` is a scalar this module decides itself: a float or float64 that
    is not NaN, or an int Python and numpy compare alike."""
    t = type(v)
    if t is float or t is _f64:
        return not _isnan(v)
    return t is int and -_INT_EXACT <= v <= _INT_EXACT


def _fwhere(c, a, b):
    """`numpy.where(c, a, b)` for a bool condition and float arms choosing a
    FINITE value, as a float64 SCALAR (the value numpy's 0-d array holds).

    Finite only: an operation meeting two NaNs returns one of them, and
    which follows the operand order -- numpy's loop on a 0-d array and its
    scalar arithmetic order them differently (`a + b`, both NaN, gave +NaN
    on one and -NaN on the other: MosLevel3 PMOS's `i` at an all-inf
    state, 2026-10-02).  With the scalar always finite, no operation it
    enters can meet two NaNs that one path ordered differently; an inf or
    NaN arm stays numpy's 0-d array."""
    tc = type(c)
    if tc is np.bool_ or tc is bool:
        ta, tb = type(a), type(b)
        if (ta is float or ta is _f64) and (tb is float or tb is _f64):
            r = a if c else b
            if _isfinite(r):
                return r if type(r) is _f64 else _f64(r)
    return np.where(c, a, b)


def _fstep(a, b):
    """`_step`: `1.0 * (a >= b)` -- the comparison is the reference's own;
    only the multiply by a `numpy.bool_` (ufunc dispatch) is replaced by
    its result."""
    r = a >= b
    if type(r) is np.bool_:
        return _ONE if r else _ZERO
    return 1.0 * r


def _pick_max(a, b, ref):
    ta, tb = type(a), type(b)
    if (ta is float or ta is _f64) and (tb is float or tb is _f64):
        if a > b:
            return a if ta is _f64 else _f64(a)
        if b > a:
            return b if tb is _f64 else _f64(b)
    return ref(a, b)


def _pick_min(a, b, ref):
    ta, tb = type(a), type(b)
    if (ta is float or ta is _f64) and (tb is float or tb is _f64):
        if a < b:
            return a if ta is _f64 else _f64(a)
        if b < a:
            return b if tb is _f64 else _f64(b)
    return ref(a, b)


def _make(hdl):
    """The helpers, bound to `hdl`'s shared reference functions (imported
    late: hdl imports this module)."""
    maxc_ref, minc_ref = hdl._maxc_numpy, hdl._minc_numpy
    npmax, npmin = np.maximum, np.minimum

    def _fmaxc(a, b):
        """`maxc`: `numpy.maximum` behind the reference's guards."""
        return _pick_max(a, b, maxc_ref)

    def _fminc(a, b):
        """`minc`: `numpy.minimum` behind the reference's guards."""
        return _pick_min(a, b, minc_ref)

    def _fmax(a, b):
        """`numpy.maximum`."""
        return _pick_max(a, b, npmax)

    def _fmin(a, b):
        """`numpy.minimum`."""
        return _pick_min(a, b, npmin)

    def compare(op, ref):
        def f(a, b):
            if _num(a) and _num(b):
                return _TRUE if op(a, b) else _FALSE
            return ref(a, b)
        f.__doc__ = f'`numpy.{ref.__name__}`, as a numpy.bool_.'
        return f

    import operator
    return {
        '_fwhere': _fwhere, '_fstep': _fstep,
        '_fmaxc': _fmaxc, '_fminc': _fminc, '_fmax': _fmax, '_fmin': _fmin,
        '_fless': compare(operator.lt, np.less),
        '_fless_equal': compare(operator.le, np.less_equal),
        '_fgreater': compare(operator.gt, np.greater),
        '_fgreater_equal': compare(operator.ge, np.greater_equal),
        '_fequal': compare(operator.eq, np.equal),
        '_fnot_equal': compare(operator.ne, np.not_equal),
    }


_HELPERS = None


def helpers():
    """The fast helpers by name: ONE set of objects per process -- the
    chain namespace carries them (`hdl._chain_namespace`), and the compile
    cache compares a namespace's objects by identity."""
    global _HELPERS
    if _HELPERS is None:
        from pycircuit.circuit import hdl
        _HELPERS = _make(hdl)
    return _HELPERS


#: the call head each fast name replaces in the generated text
RENAMES = {
    'numpy.where': '_fwhere', '_step': '_fstep',
    'maxc': '_fmaxc', 'minc': '_fminc',
    'numpy.maximum': '_fmax', 'numpy.minimum': '_fmin',
    'numpy.less': '_fless', 'numpy.less_equal': '_fless_equal',
    'numpy.greater': '_fgreater', 'numpy.greater_equal': '_fgreater_equal',
    'numpy.equal': '_fequal', 'numpy.not_equal': '_fnot_equal',
}
