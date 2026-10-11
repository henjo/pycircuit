/* branch.h -- the branch check's cheap screen: its state and entry, the
 * interface `_tran_branch_c` loads branch.c through.  cffi lays the struct
 * out from this text and branch.c includes it, so the compiler checks the
 * entry against it. */
typedef struct {
    long n;                 /* C is n x n, row-major */
    long npat;              /* > 0: only the entries `pat` names can be
                               nonzero (the core wrote C on its pattern) */
    const int64_t *pat;
    const double *C;
    double cmax, tol_rel;   /* in: the run's largest |C| so far, the screen's
                               relative tolerance (BRANCH_SCREEN_TOL) */
    double nrm, ref;        /* out: |C|'s largest entry, the new running
                               maximum */
    long verdict;
} branch_t;
long hdl_fn(branch_t *s);

/* hdl_fn's result, left in `verdict` too: `_branch_screen_C`'s cheap level */
enum branch_verdict {
    BRANCH_QUIET = 0,            /* nothing collapsed: no SVD */
    BRANCH_COLLAPSED = 1,        /* C collapsed entirely: fired, no direction */
    BRANCH_SVD = 2,              /* the SVD decides (the Python's) */
    BRANCH_BAIL_NONFINITE = 3    /* an entry not finite: the Python screens */
};
