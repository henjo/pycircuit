/* pcx_fenv.h -- the floating-point flags a unit tests, and THE FLAG RULE.
 *
 * A unit that reproduces numpy's warnings clears the flags, makes the
 * operations numpy would make, and tests PCX_FP_FLAGS (overflow, underflow,
 * invalid, division by zero) before anything its caller reads is final: any
 * flag hands the work back, and the Python makes it with whatever warnings its
 * own error state gives.  That is more than numpy reports (its default ignores
 * underflow, and some chains run under `np.errstate(... 'ignore')`): a flag
 * costs at worst a needless hand-back, never another answer or another
 * warning.
 *
 * THE FLAG RULE -- MUST KEEP in every unit that includes this header.  A test
 * must see exactly the operations made since the last clear.  GCC does not
 * honour `#pragma STDC FENV_ACCESS`: it may move a floating-point operation
 * across `feclearexcept` / `fetestexcept` (opaque library calls) unless
 * something pins it there.  So each operation a test covers stores its result
 * through a pointer the call could read (a buffer of the unit's struct), or
 * decides a branch taken before the test.  Keeping such a result only in a
 * local, fusing two loops a test separates, or moving work into or out of a
 * cleared-and-tested region can let a flag escape its test (a warning numpy
 * gives, not given) or reach another one.
 *
 * A unit includes it as `#include "pcx_fenv.h"` (expanded in place by
 * `_cbuild.csrc`, so a unit's pinned pair sees it).
 */
#include <fenv.h>

#define PCX_FP_FLAGS (FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO)
