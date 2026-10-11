#include <stdint.h>
#include <math.h>
#include "branch.h"

/* `_BranchCheck._branch_screen_C` up to its SVD: |C|'s largest entry, the
   running maximum (Python's `max(cmax, nrm)`: the first where they are
   equal), the diagonal's largest and its smallest positive entry, the
   screen's two tests -- the same comparisons on the same values, so the
   same verdict */
long hdl_fn(branch_t *s)
{
    long n = s->n, k, i;
    const double *C = s->C;
    double nrm = 0.0, ref, tol, dmax, pmin = 0.0;
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
    ref = (nrm > s->cmax) ? nrm : s->cmax;
    tol = s->tol_rel * ref;
    dmax = fabs(C[0]);
    for (i = 0; i < n; i++) {
        double d = fabs(C[i * (n + 1)]);
        if (d > dmax) dmax = d;
        if (d > 0.0 && (pmin == 0.0 || d < pmin)) pmin = d;
    }
    s->nrm = nrm;
    s->ref = ref;
    if (nrm > tol && dmax > tol && pmin > 1e-6 * ref) return s->verdict = BRANCH_QUIET;
    if (nrm <= tol) return s->verdict = BRANCH_COLLAPSED;
    return s->verdict = BRANCH_SVD;
}
