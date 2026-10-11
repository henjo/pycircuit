/* sens.c -- the shooting's coupled stage step map in one call: the plain
   dense walk's step (`_PeriodWalks._stage_step` with `_stage_block`,
   SciPy's `lu_factor` and `_StageStep.solve`'s `getrs`), bit for bit.  The
   glue, `shooting/_sens_c.step` and `_map`, checks the case, fills sens_t
   (sens.h) and takes the new P from Z; on a bail the Python step runs. */
#include <fenv.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include "pcx_abi.h"
#include "pcx_fenv.h"
#include "sens.h"
/* the four flags `_stage_block`'s errstate(all='raise') raises on (inexact
   is not one) */
/* the full row (column) of reduced index r: the reference one skipped */
#define FULL(r) ((r) < t->iref ? (r) : (r) + 1)
/* the step map: J Z = B, J[i][j] = delta_ij C_i + (h A_ij) G_j over the
   reduced rows and columns, B every stage's block C_n P.
   Assumes what `_sens_c.step` checks and `_Ctx` sizes: 0 < s, the buffers
   of sens.h's sizes, n = s m small enough for int (`MAXN`), Z set.
   Returns SENS_DONE (also in `status`): J holds dgetrf's LU, ipiv its
   pivots, Z the solution (its last stage block the new P), info 0.  Else a
   bail code, also in `status`: J, ipiv and Z possibly part-written, none
   read by the glue.  The floating-point flags are left as they are, on
   every return (none cleared). */
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
       + (h A_ij) G_j[r, c] -- `_stage_block`'s operations, under its errstate.
       MUST KEEP: each result stored (col[r]) before the fetestexcept -- the
       compiler does not model the flags (no FENV_ACCESS), so a result left
       in a register may be made after the test */
    feclearexcept(FE_ALL_EXCEPT);
    for (i = 0; i < s; i++)
        for (j = 0; j < s; j++) {
            /* h A_ij rounded first (numpy's `h * A`), then its product with
               G_j -- h (A_ij G) rounds apart */
            double hA = h * t->A[i * s + j];
            const double *G = t->Gf + j * nn, *C = t->Cf + i * nn;
            for (c = 0; c < m; c++) {
                double *col = t->J + (j * m + c) * n + i * m;
                long fc = FULL(c);
                for (r = 0; r < m; r++) {
                    long f = FULL(r) * nf + fc;
                    double v = hA * G[f];
                    /* the diagonal blocks: `C3 + T`, C on the left */
                    if (i == j) v = C[f] + v;
                    col[r] = v;
                }
            }
        }
    if (fetestexcept(PCX_FP_FLAGS)) return t->status = SENS_BAIL_FLAGS;
    int n32 = (int) n, k32 = (int) k, info = 0;
    /* `lu_factor`: SciPy's own dgetrf on the column-major block, in place */
    ((dgetrf32_t) t->dgetrf)(&n32, &n32, t->J, &n32, t->ipiv, &info);
    t->info = info;
    if (info != 0) return t->status = SENS_BAIL_SINGULAR;
    /* the right-hand side, column-major: every stage's block C_n P */
    for (c = 0; c < k; c++)
        for (i = 0; i < s; i++)
            for (r = 0; r < m; r++)
                t->Z[c * n + i * m + r] = t->base[r * k + c];
    /* `_StageStep.solve`: SciPy's own dgetrs on those factors, in place */
    char trans = 'N';
    ((dgetrs32_t) t->dgetrs)(&trans, &n32, &k32, t->J, &n32, t->ipiv, t->Z, &n32, &info, 1);
    t->info = info;
    if (info != 0) return t->status = SENS_BAIL_LAPACK;
    return t->status = SENS_DONE;
}
