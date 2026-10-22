/* lmmsens.c -- the dense multistep walk's monodromy step in one call:
   `_steps._LiveMonodromy.step` on its live columns -- `_lmm_recursion` with
   the capacitances read through their nonzero columns (`_ColumnsOf`) and
   the step's solve through numpy's `linalg.solve` -- bit for bit.  The glue,
   `shooting/_lmmsens_c.step`, checks the case, fills lmmsens_t (lmmsens.h)
   and hands Pn and Pqn to the Python step's bookkeeping; on a bail the
   Python step runs on the same state. */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include "pcx_abi.h"
#include "lmmsens.h"

/* C's nonzero columns (a NaN or inf nonzero), as `_ColumnsOf` finds them
   (``np.any(C != 0, axis=0)``), into kc; their count.  C read along its
   rows: kc first each column's flag, then (in place: nk <= j) the columns
   flagged. */
static long lmmsens_columns(long m, const double *C, int64_t *kc)
{
    long nk = 0, i, j;
    for (j = 0; j < m; j++) kc[j] = 0;
    for (i = 0; i < m; i++) {
        const double *row = C + i * m;
        for (j = 0; j < m; j++) kc[j] |= row[j] != 0.0;
    }
    for (j = 0; j < m; j++)
        if (kc[j]) kc[nk++] = j;
    return nk;
}

/* `_ColumnsOf(C) @ P` into T (m x w, row-major), kc and nk C's nonzero
   columns: numpy's matmul of C[:, kc] and P[kc], both gathered C-contiguous
   as numpy gathers them -- `cblas_dgemm` (row-major, neither transposed)
   for two columns or more, numpy's own loop for one (``0 + a b`` an entry:
   a one-column inner dimension is its special case), zeros for none (a
   zero dimension). */
static void lmmsens_product(lmmsens_t *t, const double *C, const int64_t *kc, long nk,
                            const double *P, double *T)
{
    long m = t->m, w = t->w, i, j, c;
    if (nk == 0) {
        for (i = 0; i < m * w; i++) T[i] = 0.0;
        return;
    }
    for (i = 0; i < m; i++)
        for (c = 0; c < nk; c++) t->Ck[i * nk + c] = C[i * m + kc[c]];
    for (c = 0; c < nk; c++)
        for (j = 0; j < w; j++) t->Pk[c * w + j] = P[kc[c] * w + j];
    if (nk == 1) {
        for (i = 0; i < m; i++)
            for (j = 0; j < w; j++) {
                double v = 0.0;
                v += t->Ck[i] * t->Pk[j];
                T[i * w + j] = v;
            }
        return;
    }
    ((dgemm64_t) t->dgemm)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, PCX_CBLAS_NO_TRANS,
                           (int64_t) m, (int64_t) w, (int64_t) nk, 1.0, t->Ck, (int64_t) nk,
                           t->Pk, (int64_t) w, 0.0, T, (int64_t) w);
}

/* the step: S = (b Pq where b is nonzero, else zeros) + sum_k alpha_k (C_k
   @ P_k) over k = 1 .. L-1 (C_1, P_0 and C_2, P_1), each term added as
   numpy adds it (the product scaled, then summed); Pn = -solve(Jf, S) --
   numpy's own dgesv on the column-major copies `linalg.solve` makes;
   Pqn = alpha_0 (Cn @ Pn) + S.
   Assumes what `_lmmsens_c.step` checks and `_Ctx` sizes: 2 <= m, 2 <= w
   <= wcap, 2 <= L <= 3, the inputs C-contiguous m x w or m x m float64, Pn
   and Pqn set, `reuse` true of the arrays (kcs and nks are then the last
   call's).  The three column lists are made first, whatever the return:
   kcs and nks then hold this call's.  Returns LMMSENS_DONE (also in `status`), Pn and Pqn
   written and finite, info 0.  Else a bail code, also in `status`: Pn and
   Pqn possibly part-written, not read by the glue. */
long pcx_lmmsens(lmmsens_t *t)
{
    long m = t->m, w = t->w, mw = m * w, i, r, c, k;
    double *S = t->S, *T = t->T;
    int64_t *k1 = t->kcs, *k2 = t->kcs + m, *kn = t->kcs + 2 * m;
    /* the three lists: C2's from the last C1's (before C1's moves on), C1's
       from the last Cn's, Cn's read */
    if (t->reuse & 2) {
        for (i = 0; i < t->nks[0]; i++) k2[i] = k1[i];
        t->nks[1] = t->nks[0];
    } else
        t->nks[1] = lmmsens_columns(m, t->C2, k2);
    if (t->reuse & 1) {
        for (i = 0; i < t->nks[2]; i++) k1[i] = kn[i];
        t->nks[0] = t->nks[2];
    } else
        t->nks[0] = lmmsens_columns(m, t->C1, k1);
    t->nks[2] = lmmsens_columns(m, t->Cn, kn);
    /* Python's `b * Pq if b`: a NaN is true, -0.0 false */
    if (t->b != 0.0)
        for (i = 0; i < mw; i++) S[i] = t->b * t->Pq[i];
    else
        for (i = 0; i < mw; i++) S[i] = 0.0;
    for (k = 1; k < t->L; k++) {
        double a = t->alpha[k];
        if (k == 1) lmmsens_product(t, t->C1, k1, t->nks[0], t->P0, T);
        else lmmsens_product(t, t->C2, k2, t->nks[1], t->P1, T);
        for (i = 0; i < mw; i++) S[i] = S[i] + a * T[i];
    }
    /* `linalg.solve(Jf, S)`: numpy's dgesv on Jf and S, each copied
       column-major as its gufunc copies them */
    for (c = 0; c < m; c++)
        for (r = 0; r < m; r++) t->A[c * m + r] = t->Jf[r * m + c];
    for (c = 0; c < w; c++)
        for (r = 0; r < m; r++) t->X[c * m + r] = S[r * w + c];
    int64_t m64 = m, w64 = w, info = 0;
    ((dgesv64_t) t->dgesv)(&m64, &w64, t->A, &m64, t->ipiv, t->X, &m64, &info);
    t->info = info;
    if (info != 0) return t->status = LMMSENS_BAIL_SINGULAR;
    for (r = 0; r < m; r++)
        for (c = 0; c < w; c++) t->Pn[r * w + c] = -t->X[c * m + r];
    lmmsens_product(t, t->Cn, kn, t->nks[2], t->Pn, T);
    double a0 = t->alpha[0];
    for (i = 0; i < mw; i++) t->Pqn[i] = a0 * T[i] + S[i];
    /* the Python step's own test: a result not finite runs at full width */
    for (i = 0; i < mw; i++)
        if (!isfinite(t->Pn[i]) || !isfinite(t->Pqn[i])) return t->status = LMMSENS_BAIL_NONFINITE;
    return t->status = LMMSENS_DONE;
}
