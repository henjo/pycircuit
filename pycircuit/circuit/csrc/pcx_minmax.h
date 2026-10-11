/* pcx_minmax.h -- min and max as Python and numpy pick them.
 *
 * On a tie, or a NaN, which operand comes back is part of the bits, and the
 * three rules differ:
 *   pcx_py_max(a, b)       Python's `max(a, b)`: `a` unless `b` is strictly
 *                          greater (a tie, or a NaN either side, keeps `a`);
 *   pcx_py_min(a, b)       Python's `min(a, b)`: `a` unless `b` is strictly
 *                          smaller;
 *   pcx_np_maximum(a, b)   numpy's `maximum(a, b)` for a non-NaN `a`: `a`
 *                          unless `b` is greater or NaN (numpy keeps a NaN
 *                          `a` as well: a caller passes a finite one);
 *   pcx_np_clip(x, lo, hi) numpy's `clip(x, lo, hi)` for a non-NaN `x`, as
 *                          min(max(x, lo), hi) (numpy keeps a NaN `x`).
 * A call names the Python it stands for with its operands in the Python's
 * order: `max(abs(t), 1.0)` is `pcx_py_max(fabs(t), 1.0)`.
 *
 * A unit includes it as `#include "pcx_minmax.h"` (expanded in place by
 * `_cbuild.csrc`, so a unit's pinned pair sees it).
 */

static inline double pcx_py_max(double a, double b)
{
    return (b > a) ? b : a;
}

static inline double pcx_py_min(double a, double b)
{
    return (b < a) ? b : a;
}

static inline double pcx_np_maximum(double a, double b)
{
    return (a >= b) ? a : b;
}

static inline double pcx_np_clip(double x, double lo, double hi)
{
    double m = (x > lo) ? x : lo;
    return (m < hi) ? m : hi;
}
