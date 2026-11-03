#include <math.h>
#include "libm.h"
/* THE C LIBRARY'S FUNCTIONS OVER AN ARRAY, for the libm probe: `out[i] =
   f(a[i])` (`f(a[i], b[i])` for atan2, hypot and pow) -- built with the
   kernels' flags (`_cbuild.CFLAGS`: no builtin of any of them), so each is
   the C library's own entry point, the one a kernel calls.  0, or 1 for a
   function it does not know. */
long pcx_libm(long which, const double *a, const double *b, long n, double *out)
{
    long i;
    switch (which) {
    case LIBM_EXP: for (i = 0; i < n; i++) out[i] = exp(a[i]); break;
    case LIBM_LOG: for (i = 0; i < n; i++) out[i] = log(a[i]); break;
    case LIBM_LOG1P: for (i = 0; i < n; i++) out[i] = log1p(a[i]); break;
    case LIBM_EXPM1: for (i = 0; i < n; i++) out[i] = expm1(a[i]); break;
    case LIBM_SIN: for (i = 0; i < n; i++) out[i] = sin(a[i]); break;
    case LIBM_COS: for (i = 0; i < n; i++) out[i] = cos(a[i]); break;
    case LIBM_TAN: for (i = 0; i < n; i++) out[i] = tan(a[i]); break;
    case LIBM_ASIN: for (i = 0; i < n; i++) out[i] = asin(a[i]); break;
    case LIBM_ACOS: for (i = 0; i < n; i++) out[i] = acos(a[i]); break;
    case LIBM_ATAN: for (i = 0; i < n; i++) out[i] = atan(a[i]); break;
    case LIBM_SINH: for (i = 0; i < n; i++) out[i] = sinh(a[i]); break;
    case LIBM_COSH: for (i = 0; i < n; i++) out[i] = cosh(a[i]); break;
    case LIBM_TANH: for (i = 0; i < n; i++) out[i] = tanh(a[i]); break;
    case LIBM_ATAN2: for (i = 0; i < n; i++) out[i] = atan2(a[i], b[i]); break;
    case LIBM_HYPOT: for (i = 0; i < n; i++) out[i] = hypot(a[i], b[i]); break;
    case LIBM_POW: for (i = 0; i < n; i++) out[i] = pow(a[i], b[i]); break;
    default: return 1;
    }
    return 0;
}
