/* pcx_dense.h -- the reduced system's dense plumbing, written once (the C
 * stepping loop plan's stage 2): the Python solves with the reference row
 * and column removed (`analysis.remove_row_col`) and puts the reference back
 * as +0.0 (`insert_row`); its dense solves hand LAPACK a column-major copy
 * (numpy's, before `dgesv`; SciPy's, before `dgetrf`).
 *
 * A unit includes it as `#include "pcx_dense.h"` (expanded in place by
 * `_cbuild.csrc`, so a unit's pinned pair sees it).  `r` is the reference
 * row, 0 <= r < n; the reduced size is n - 1.
 */
#include <math.h>

/* remove_row_col on a vector: `red` (n - 1) is `full` (n) without row r */
static inline void pcx_drop_ref(const double *full, double *red, long n, long r)
{
    long i;
    for (i = 0; i < r; i++) red[i] = full[i];
    for (i = r + 1; i < n; i++) red[i - 1] = full[i];
}

/* insert_row: `full` (n) is `red` (n - 1) with +0.0 at row r */
static inline void pcx_put_ref(const double *red, double *full, long n, long r)
{
    long i;
    for (i = 0; i < r; i++) full[i] = red[i];
    full[r] = 0.0;
    for (i = r; i < n - 1; i++) full[i + 1] = red[i];
}

/* the reduced matrix of the row-major n x n `J`, column-major into `A`
   ((n - 1) x (n - 1): LAPACK's layout); nonzero where every entry copied
   is finite */
static inline int pcx_reduce_cm(const double *J, double *A, long n, long r)
{
    long m = n - 1, ri, ci, rr, cc;
    int ok = 1;
    for (ri = 0, rr = 0; ri < n; ri++) {
        if (ri == r) continue;
        for (ci = 0, cc = 0; ci < n; ci++) {
            double v;
            if (ci == r) continue;
            v = J[ri * n + ci];
            ok &= isfinite(v) != 0;
            A[cc * m + rr] = v;
            cc++;
        }
        rr++;
    }
    return ok;
}

/* `pcx_reduce_cm`, and |J| reduced, row-major, into `absJ` (the Newton's
   convergence test reads it) */
static inline int pcx_reduce_cm_abs(const double *J, double *A, double *absJ, long n, long r)
{
    long m = n - 1, ri, ci, rr, cc;
    int ok = 1;
    for (ri = 0, rr = 0; ri < n; ri++) {
        if (ri == r) continue;
        for (ci = 0, cc = 0; ci < n; ci++) {
            double v;
            if (ci == r) continue;
            v = J[ri * n + ci];
            ok &= isfinite(v) != 0;
            A[cc * m + rr] = v;
            absJ[rr * m + cc] = fabs(v);
            cc++;
        }
        rr++;
    }
    return ok;
}
