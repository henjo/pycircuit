/* sensklu.c -- the shooting's coupled stage step map on the walk's KLU in one
   call: `_KLUWalk.factor_values` of the C stage Newton's last assembly (the
   walk's numeric refactored, its check, a fresh factor where it fails) and
   `_StageStep.solve`'s stage system on it, bit for bit.  The glue,
   `shooting/_sensklu_c`, checks the case, fills sensklu_t (sensklu.h) and
   takes the new P from X; on a bail the Python step runs. */
#include <fenv.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include "pcx_abi.h"
#include "pcx_minmax.h"
#include "sensklu.h"

/* the flags numpy's arithmetic warns of (underflow is not one) */
#define SK_FLAGS (FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

/* `_KLUWalk._holds`: b = np.bincount(Ai, weights=Ax, minlength=M) (each
   row from +0.0, the entries added in the pattern's order), x = b solved,
   r = np.bincount(Ai, weights=Ax * x[cols]) - b, then max |r| / max(max |b|,
   1e-300) <= tol -- numpy's max propagating a NaN, Python's max(a, b) `a`
   unless `b` is greater.  1 it holds, 0 not (a failed solve: not), -1 a
   flag numpy would warn of (KLU's own flags are not numpy's: cleared) */
static int holds(sensklu_t *t)
{
    long i, p, M = t->M;
    double maxb = 0.0, maxr = 0.0;
    int nanb = 0, nanr = 0, ok;
    for (i = 0; i < M; i++) t->pb[i] = 0.0;
    for (p = 0; p < t->nnz; p++) t->pb[t->kAi[p]] += t->Ax[p];
    for (i = 0; i < M; i++) t->px[i] = t->pb[i];
    if (fetestexcept(SK_FLAGS)) return -1;
    ok = ((klu_solve_t) t->klu_solve)(t->klu_sym, t->num, (int) M, 1, t->px, t->klu_common);
    feclearexcept(FE_ALL_EXCEPT);
    if (!ok) return 0;
    for (i = 0; i < M; i++) t->pr[i] = 0.0;
    for (p = 0; p < t->nnz; p++) t->pr[t->kAi[p]] += t->Ax[p] * t->px[t->cols[p]];
    for (i = 0; i < M; i++) t->pr[i] = t->pr[i] - t->pb[i];
    for (i = 0; i < M; i++) {
        double ab = fabs(t->pb[i]), ar = fabs(t->pr[i]);
        if (isnan(ab)) nanb = 1;
        else if (ab > maxb) maxb = ab;
        if (isnan(ar)) nanr = 1;
        else if (ar > maxr) maxr = ar;
    }
    if (fetestexcept(SK_FLAGS)) return -1;
    if (nanb || nanr) return 0;
    {
        double scale = pcx_py_max(maxb, 1e-300);
        double r = maxr / scale;
        if (fetestexcept(SK_FLAGS)) return -1;
        return r <= t->resid_tol;
    }
}

/* the step: the walk's numeric made for this block (`factor_values`), the
   stage system's right-hand side -- `np.vstack` of `base` once a stage,
   copied column-major as `_KLUNumeric._run` copies it -- and klu_solve of
   its k columns at once.  Assumes the buffers of sensklu.h's sizes, Ax and
   base finite or not (the Python reads them as they are), k >= 1. */
long pcx_sensklu(sensklu_t *t)
{
    long M = t->M, m = t->m, k = t->k, s = t->s, i, r, c;
    int h;
    t->refactors = t->fallbacks = t->factors = 0;
    feclearexcept(FE_ALL_EXCEPT);
    if (t->num) {
        int ok = ((klu_refactor_t) t->klu_refactor)(t->kAp, t->kAi, t->Ax,
                                                    t->klu_sym, t->num, t->klu_common);
        feclearexcept(FE_ALL_EXCEPT);
        h = ok ? holds(t) : 0;
        if (h < 0) return t->status = SENSKLU_BAIL_FLAGS;
        t->refactors = 1;
        if (!h) {
            void *old = t->num;
            t->fallbacks = 1;
            t->num = NULL;
            ((klu_free_numeric_t) t->klu_free)(&old, t->klu_common);
        }
    }
    if (!t->num) {
        t->num = ((klu_factor_t) t->klu_factor)(t->kAp, t->kAi, t->Ax, t->klu_sym,
                                                t->klu_common);
        feclearexcept(FE_ALL_EXCEPT);
        if (!t->num) return t->status = SENSKLU_SINGULAR;
        t->factors = 1;
    }
    for (c = 0; c < k; c++)
        for (i = 0; i < s; i++)
            for (r = 0; r < m; r++) t->X[c * M + i * m + r] = t->base[r * k + c];
    if (!((klu_solve_t) t->klu_solve)(t->klu_sym, t->num, (int) M, (int) k, t->X,
                                      t->klu_common)) {
        feclearexcept(FE_ALL_EXCEPT);
        return t->status = SENSKLU_SOLVE_FAILED;
    }
    feclearexcept(FE_ALL_EXCEPT);
    return t->status = SENSKLU_DONE;
}
