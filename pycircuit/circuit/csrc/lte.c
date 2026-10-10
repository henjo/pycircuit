
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <fenv.h>
#include "pcx_abi.h"
#include "lte.h"

#define LTE_FLAGS (FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO)

static long lte_bail(lte_t *s, long why)
{
    feclearexcept(FE_ALL_EXCEPT);
    s->status = why;
    return why;
}

long hdl_fn(lte_t *s)
{
    long n = s->n, m = n - 1, r = s->iref, i, j, ri, ci, rr, cc;
    int64_t m64 = m, one64 = 1, info64 = 0;
    double err;
    /* every input finite: anything else is the Python chain's */
    if (!(isfinite(s->h1) && isfinite(s->h2) && isfinite(s->h3) && isfinite(s->s12)
          && isfinite(s->s23) && isfinite(s->s123) && isfinite(s->coef)
          && isfinite(s->reltol) && isfinite(s->trtol) && isfinite(s->abstol)))
        return lte_bail(s, LTE_BAIL_NONFINITE);
    for (i = 0; i < n; i++)
        if (!(isfinite(s->xc[i]) && isfinite(s->xl[i]) && isfinite(s->q[i])
              && isfinite(s->q1[i]) && isfinite(s->q2[i]) && isfinite(s->q3[i])))
            return lte_bail(s, LTE_BAIL_NONFINITE);
    if (s->have_run)
        for (i = 0; i < n; i++)
            if (!isfinite(s->run[i]) || signbit(s->run[i])) return lte_bail(s, LTE_BAIL_NONFINITE);
    if (s->abs_vec)
        for (i = 0; i < n; i++) if (!isfinite(s->abst[i])) return lte_bail(s, LTE_BAIL_NONFINITE);
    feclearexcept(FE_ALL_EXCEPT);
    /* compute_lte: coef * third_divided_difference, numpy's elementwise
       operations in their order (each result stored, so each operation
       precedes the flag test); the step sums and the coefficient are the
       glue's, in the integrator's own scalar arithmetic */
    for (i = 0; i < n; i++) {
        double d1 = (s->q[i] - s->q1[i]) / s->h1;
        double d2 = (s->q1[i] - s->q2[i]) / s->h2;
        double d3 = (s->q2[i] - s->q3[i]) / s->h3;
        double dda = (d1 - d2) / s->s12;
        double ddb = (d2 - d3) / s->s23;
        s->lte[i] = s->coef * ((dda - ddb) / s->s123);
    }
    if (fetestexcept(LTE_FLAGS)) return lte_bail(s, LTE_BAIL_FLAGS);
    /* remove_row_col and numpy's solve: the reduced J column-major (numpy's
       copy for LAPACK), the reduced Eg, numpy's own dgesv */
    for (ri = 0, rr = 0; ri < n; ri++) {
        if (ri == r) continue;
        for (ci = 0, cc = 0; ci < n; ci++) {
            double v;
            if (ci == r) continue;
            v = s->J[ri * n + ci];
            if (!isfinite(v)) return lte_bail(s, LTE_BAIL_NONFINITE);
            s->A[cc * m + rr] = v;
            cc++;
        }
        s->b[rr] = s->lte[ri];
        rr++;
    }
    ((dgesv64_t) s->dgesv)(&m64, &one64, s->A, &m64, s->ipiv, s->b, &m64, &info64);
    /* (numpy's solve clears what LAPACK raised; an info is its LinAlgError) */
    feclearexcept(FE_ALL_EXCEPT);
    if (info64 != 0) return lte_bail(s, LTE_BAIL_DGESV);
    for (i = 0; i < m; i++) if (!isfinite(s->b[i])) return lte_bail(s, LTE_BAIL_LTE);
    /* the full vector, the reference row restored as +0.0 */
    for (i = 0, j = 0; i < n; i++) s->lte[i] = (i == r) ? 0.0 : s->b[j++];
    /* _reference: numpy's maximum (the first where they compare equal) of
       |x_curr| and |x_last|, folded into the running maximum */
    for (i = 0; i < n; i++) {
        double a = fabs(s->xc[i]), c = fabs(s->xl[i]);
        double loc = (a >= c) ? a : c;
        if (s->relref == LTE_RELREF_POINTLOCAL) {
            s->ref[i] = loc;
        } else if (s->have_run) {
            double rv = s->run[i];
            s->out[i] = (rv >= loc) ? rv : loc;
        } else {
            s->out[i] = loc;
        }
    }
    if (s->relref == LTE_RELREF_ALLLOCAL) {
        for (i = 0; i < n; i++) s->ref[i] = s->out[i];
    } else if (s->relref == LTE_RELREF_SIGGLOBAL) {
        /* sigglobal_reference: each unit group's maximum, broadcast back */
        long nn = s->nn, lo, hi, g;
        long cut = (nn <= 0 || nn >= n) ? n : nn;
        for (g = 0; g < 2; g++) {
            double mx;
            lo = g ? cut : 0;
            hi = g ? n : cut;
            if (lo >= hi) continue;
            mx = s->out[lo];
            for (i = lo + 1; i < hi; i++) if (s->out[i] > mx) mx = s->out[i];
            for (i = lo; i < hi; i++) s->ref[i] = mx;
        }
    }
    /* the tolerance TRTOL (reltol ref + abstol) and |lte| / tol */
    for (i = 0; i < n; i++)
        s->tq[i] = s->trtol * (s->reltol * s->ref[i] + (s->abs_vec ? s->abst[i] : s->abstol));
    if (fetestexcept(LTE_FLAGS)) return lte_bail(s, LTE_BAIL_FLAGS);
    for (i = 0; i < n; i++) if (!(s->tq[i] > 0.0 && isfinite(s->tq[i]))) return lte_bail(s, LTE_BAIL_TOLERANCE);
    for (i = 0; i < n; i++) s->tq[i] = fabs(s->lte[i]) / s->tq[i];
    if (fetestexcept(LTE_FLAGS)) return lte_bail(s, LTE_BAIL_FLAGS);
    /* `_normalised`'s mask (`Transient.lte_algebraic`): an unknown not
       judged is +0.0, as `np.where(mask, e, 0.0)` makes it */
    if (s->has_mask)
        for (i = 0; i < n; i++) if (s->mask[i] == 0.0) s->tq[i] = 0.0;
    /* np.max: no NaN, no negative zero (|lte| / a positive tolerance) */
    err = s->tq[0];
    for (i = 1; i < n; i++) if (s->tq[i] > err) err = s->tq[i];
    s->err = err;
    feclearexcept(FE_ALL_EXCEPT);
    s->status = LTE_DONE;
    return LTE_DONE;
}
