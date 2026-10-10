
#include <fenv.h>
#include <stddef.h>
#include <stdint.h>
#include "pcx_abi.h"
#include "sens.h"
#define SENS_FLAGS (FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO)
/* the full row (column) of reduced index r: the reference one skipped */
#define FULL(r) ((r) < t->iref ? (r) : (r) + 1)
long hdl_fn(sens_t *t)
{
    long s = t->s, m = t->m, k = t->k, nf = t->nf, n = s * m, nn = nf * nf;
    long i, j, r, c;
    double h = t->h;
    /* every reduced entry finite: `_stage_block` declines otherwise, and the
       loop's block makes `lu_factor` raise */
    if (!isfinite(h)) return t->status = SENS_BAIL_NONFINITE;
    for (i = 0; i < s * s; i++) if (!isfinite(t->A[i])) return t->status = SENS_BAIL_NONFINITE;
    for (i = 0; i < s; i++)
        for (r = 0; r < m; r++)
            for (c = 0; c < m; c++) {
                long f = i * nn + FULL(r) * nf + FULL(c);
                if (!isfinite(t->Cf[f]) || !isfinite(t->Gf[f])) return t->status = SENS_BAIL_NONFINITE;
            }
    /* C_n P finite: `lu_solve` refuses it otherwise */
    for (i = 0; i < m * k; i++) if (!isfinite(t->base[i])) return t->status = SENS_BAIL_RHS;
    /* the block, column-major: J[(j m + c) n + i m + r] = [i == j] C_i[r, c]
       + (h A_ij) G_j[r, c] -- `_stage_block`'s operations, under its errstate */
    feclearexcept(FE_ALL_EXCEPT);
    for (i = 0; i < s; i++)
        for (j = 0; j < s; j++) {
            double hA = h * t->A[i * s + j];
            const double *G = t->Gf + j * nn, *C = t->Cf + i * nn;
            for (c = 0; c < m; c++) {
                double *col = t->J + (j * m + c) * n + i * m;
                long fc = FULL(c);
                for (r = 0; r < m; r++) {
                    long f = FULL(r) * nf + fc;
                    double v = hA * G[f];
                    if (i == j) v = C[f] + v;
                    col[r] = v;
                }
            }
        }
    if (fetestexcept(SENS_FLAGS)) return t->status = SENS_BAIL_FLAGS;
    int n32 = (int) n, k32 = (int) k, info = 0;
    ((dgetrf32_t) t->dgetrf)(&n32, &n32, t->J, &n32, t->ipiv, &info);
    t->info = info;
    if (info != 0) return t->status = SENS_BAIL_SINGULAR;
    /* the right-hand side, column-major: every stage's block C_n P */
    for (c = 0; c < k; c++)
        for (i = 0; i < s; i++)
            for (r = 0; r < m; r++)
                t->Z[c * n + i * m + r] = t->base[r * k + c];
    char trans = 'N';
    ((dgetrs32_t) t->dgetrs)(&trans, &n32, &k32, t->J, &n32, t->ipiv, t->Z, &n32, &info, 1);
    t->info = info;
    if (info != 0) return t->status = SENS_BAIL_LAPACK;
    return t->status = SENS_DONE;
}
