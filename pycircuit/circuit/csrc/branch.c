/* branch.c -- the branch check's cheap screen in C, bit for bit the
   verdict of `_BranchCheck._branch_screen_C` up to its SVD, the
   certificate included; its glue is `_tran_branch_c`, and the C loop
   (loop.c) calls it after each solve.
   No floating-point flag is tested: past the finiteness test the cheap
   level's every operation is a fabs, a comparison or a product of finite
   values, and the Python's reductions and products warn about none of
   them; the certificate's (numpy's own dgesv, a sum of squares, a square
   root) the Python makes under `np.errstate(all='ignore')`, deciding by
   the same comparisons on the same values. */
#include <stdint.h>
#include <float.h>
#include <math.h>
#include "pcx_abi.h"
#include "pcx_minmax.h"
#include "branch.h"

/* THE CERTIFICATE (`_tran_branch._certify`, operation for operation): 1
   where every entry of B = C[ci, cj] between two of its blocks is zero and
   every block's bound 1 / ||B_k^-1||_F -- |b| for a block of one -- is
   above tol (1 + delta_k), delta_k = growth k eps ||B_k||_F ||B_k^-1||_F
   (the inverse's rounding, at the block's own conditioning); `bound` the
   smallest.  0 otherwise: an entry between blocks, a singular block
   (dgesv's info), a bound not above (NaN included).  B_k and its inverse in
   `work` column-major, as numpy's inv hands them to LAPACK; each sum of
   squares sequential in that order (numpy's cumsum). */
static int certify(branch_t *s, const double *C, double tol)
{
    long n = s->n, r0 = s->r0, b, i, j;
    double worst = INFINITY;
    for (i = 0; i < r0; i++) {
        const double *row = C + s->ci[i] * n;
        for (j = 0; j < r0; j++)
            if (s->bid[i] != s->bid[j] && row[s->cj[j]] != 0.0) return 0;
    }
    for (b = 0; b < s->nb; b++) {
        int64_t o = s->bstart[b], k = s->bstart[b + 1] - o, info = 0, p, q;
        double bound, delta = 0.0;
        if (k == 1) {
            bound = fabs(C[s->ci[o] * n + s->cj[o]]);
        } else {
            double *A = s->work, *X = s->work + k * k, a2 = 0.0, f = 0.0;
            for (q = 0; q < k; q++)
                for (p = 0; p < k; p++) {
                    A[p + q * k] = C[s->ci[o + p] * n + s->cj[o + q]];
                    X[p + q * k] = p == q ? 1.0 : 0.0;
                }
            for (p = 0; p < k * k; p++) a2 += A[p] * A[p];
            ((dgesv64_t) s->dgesv)(&k, &k, A, &k, s->ipiv, X, &k, &info);
            if (info != 0) return 0;
            for (p = 0; p < k * k; p++) f += X[p] * X[p];
            bound = f > 0.0 ? 1.0 / sqrt(f) : INFINITY;
            delta = s->growth * (double) k * DBL_EPSILON * sqrt(a2) * sqrt(f);
        }
        if (!(bound > tol * (1.0 + delta))) return 0;
        if (bound < worst) worst = bound;
    }
    s->bound = worst;
    return 1;
}

/* `_BranchCheck._branch_screen_C` up to its SVD: |C|'s largest entry, the
   running maximum (Python's `max(cmax, nrm)`: the first where they are
   equal), the diagonal's largest and its smallest positive entry, the
   screen's two tests -- the same comparisons on the same values, so the
   same verdict.
   Assumes n >= 1 and `C` n x n; with npat > 0, every entry off the pattern
   is +0.0 (the core wrote C on its pattern: the glue's to ensure) -- only
   the pattern's entries are then read for |C|'s maximum and tested
   finite, the diagonal read whole.
   Where the cheap level fires and C has not collapsed, the certificate
   (r0 > 0): QUIET where it proves the rank (`ncert` counted, `bound`
   written), else SVD.
   QUIET, COLLAPSED, SVD: `nrm` and `ref` written.  BRANCH_BAIL_NONFINITE:
   an entry read is not finite, `nrm` and `ref` not written.  `verdict` is
   always set; the running maximum is the caller's to assign. */
long pcx_branch(branch_t *s)
{
    long n = s->n, k, i;
    const double *C = s->C;
    double nrm = 0.0, ref, tol, dmax, pmin = 0.0;
    /* `Ca.max()`: no NaN (tested) and no -0.0 (a fabs), so a running `>`
       from +0.0 is numpy's maximum in any order; the entries off the
       pattern are +0.0 and change no maximum */
    if (s->npat) {
        for (k = 0; k < s->npat; k++) {
            double v = fabs(C[s->pat[k]]);
            if (!isfinite(v)) return s->verdict = BRANCH_BAIL_NONFINITE;
            if (v > nrm) nrm = v;
        }
    } else {
        for (k = 0; k < n * n; k++) {
            double v = fabs(C[k]);
            if (!isfinite(v)) return s->verdict = BRANCH_BAIL_NONFINITE;
            if (v > nrm) nrm = v;
        }
    }
    /* `max(self._branch_cmax, nrm)`: cmax kept unless nrm is strictly
       greater (fmax would not promise which of two equal zeros) */
    ref = pcx_py_max(s->cmax, nrm);
    tol = s->tol_rel * ref;
    /* the diagonal's `d.max()` and `pos.min()` of its positive entries
       (0.0 where there are none, as `... if pos.size else 0.0`): order-free
       on finite values without -0.0 */
    dmax = fabs(C[0]);
    for (i = 0; i < n; i++) {
        double d = fabs(C[i * (n + 1)]);
        if (d > dmax) dmax = d;
        if (d > 0.0 && (pmin == 0.0 || d < pmin)) pmin = d;
    }
    s->nrm = nrm;
    s->ref = ref;
    /* the Python's two tests, in its order */
    if (nrm > tol && dmax > tol && pmin > 1e-6 * ref) return s->verdict = BRANCH_QUIET;
    if (nrm <= tol) return s->verdict = BRANCH_COLLAPSED;
    if (s->r0 > 0 && certify(s, C, tol)) {
        s->ncert++;
        return s->verdict = BRANCH_QUIET;
    }
    return s->verdict = BRANCH_SVD;
}
