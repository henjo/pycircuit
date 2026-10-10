
#include <stdint.h>
#include "pcx_abi.h"
#include "walk.h"
long hdl_fn(void **F, double **PR, const int64_t *TI, double T,
            const int64_t *NM, const int64_t *OFF, const int64_t *K,
            double *xg, const double *x0g, long start, long n)
{
    double o[64], b[64];
    long e, j, k;
    for (e = start; e < n; e++) {
        const int64_t *nm;
        if (F[e] == 0) return e;
        k = K[e];
        nm = NM + OFF[e];
        for (j = 0; j < k; j++) { o[j] = xg[nm[j]]; b[j] = x0g[nm[j]]; }
        if (TI[e] >= 0) PR[e][TI[e]] = T;
        if (((lim_fn_t) F[e])(b, PR[e], o)) return e;
        for (j = 0; j < k; j++) xg[nm[j]] = o[j];
    }
    return n;
}
