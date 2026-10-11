/* branch.c -- the branch check's cheap screen in C, bit for bit the
   verdict of `_BranchCheck._branch_screen_C` up to its SVD; its glue is
   `_tran_branch_c`, and the C loop (loop.c) calls it after each solve.
   No floating-point flag is tested: past the finiteness test every
   operation is a fabs, a comparison or a product of finite values, and
   the Python's reductions and products warn about none of them. */
#include <stdint.h>
#include <math.h>
#include "branch.h"

/* `_BranchCheck._branch_screen_C` up to its SVD: |C|'s largest entry, the
   running maximum (Python's `max(cmax, nrm)`: the first where they are
   equal), the diagonal's largest and its smallest positive entry, the
   screen's two tests -- the same comparisons on the same values, so the
   same verdict.
   Assumes n >= 1 and `C` n x n; with npat > 0, every entry off the pattern
   is +0.0 (the core wrote C on its pattern: the glue's to ensure) -- only
   the pattern's entries are then read for |C|'s maximum and tested
   finite, the diagonal read whole.
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
    ref = (nrm > s->cmax) ? nrm : s->cmax;
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
    return s->verdict = BRANCH_SVD;
}
