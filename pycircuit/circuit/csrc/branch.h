/* branch.h -- the branch check's cheap screen: its state and entry, the
 * interface `_tran_branch_c` loads branch.c through.  cffi lays the struct
 * out from this text and branch.c includes it, so the compiler checks the
 * entry against it.
 *
 * WHO SETS WHAT.  `_tran_branch_c.screen` keeps one struct a process and
 * sets every input at each call (`pat` re-pointed only when the pattern
 * object changes).  In the C loop, `_tran_loop_c._bind` makes the run's own
 * and sets `n`, `npat`, `pat` and `tol_rel` once a run; loop.c sets `C`
 * and `cmax` at each attempt and keeps `ref` as the new `cmax`. */
typedef struct {
    long n;                 /* C is n x n, row-major */
    long npat;              /* > 0: only the entries `pat` names can be
                               nonzero (the core wrote C on its pattern) */
    const int64_t *pat;     /* npat: the pattern's entries, as indices into
                               C's n * n values; read where npat > 0 */
    const double *C;        /* the matrix screened (`C` at the converged
                               point); only read */
    double cmax, tol_rel;   /* in: the run's largest |C| so far, the screen's
                               relative tolerance (BRANCH_SCREEN_TOL) */
    double nrm, ref;        /* out: |C|'s largest entry, the new running
                               maximum (not written on a bail) */
    long verdict;           /* out: hdl_fn's result, enum branch_verdict */
} branch_t;
long hdl_fn(branch_t *s);

/* hdl_fn's result, left in `verdict` too: `_branch_screen_C`'s cheap level */
enum branch_verdict {
    BRANCH_QUIET = 0,            /* nothing collapsed: no SVD */
    BRANCH_COLLAPSED = 1,        /* C collapsed entirely: fired, no direction */
    BRANCH_SVD = 2,              /* the SVD decides (the Python's) */
    BRANCH_BAIL_NONFINITE = 3    /* an entry not finite: the Python screens */
};
