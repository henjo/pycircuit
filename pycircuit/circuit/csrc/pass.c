
#include "pcx_abi.h"
#include "pass.h"
/* the class's kernel `fn` for the pass `want` (one bit) over `n` elements:
   each one's outputs in its row of OUT, handed in every pass's place (only
   the wanted pass writes) */
void pcx_pass(const void *fn, long want, const double *X, double **PR, double *OUT,
            long n, long sx, long so, long ti, double T)
{
    fused_fn_t f = (fused_fn_t) fn;
    long e;
    double *o;
    if (ti >= 0) {
        for (e = 0; e < n; e++) {
            PR[e][ti] = T;
            o = OUT + e * so;
            f(X + e * sx, PR[e], want, o, o, o, o);
        }
    } else {
        for (e = 0; e < n; e++) {
            o = OUT + e * so;
            f(X + e * sx, PR[e], want, o, o, o, o);
        }
    }
}
