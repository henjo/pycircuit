/* branch.h -- the branch check's cheap screen: its state and entry, the
 * interface `_tran_branch_c` loads branch.c through.  cffi lays the struct
 * out from this text and branch.c includes it, so the compiler checks the
 * entry against it.
 *
 * WHO SETS WHAT.  `_tran_branch_c.screen` keeps one struct a process and
 * sets every input at each call (`pat` re-pointed only when the pattern
 * object changes, the certificate's block only when the run's data or
 * growth change).  In the C loop, `_tran_loop_c._bind` makes the run's own
 * and sets `n`, `npat`, `pat`, `tol_rel` and the certificate's block once a
 * run; loop.c sets `C` and `cmax` at each attempt and keeps `ref` as the new
 * `cmax`.  Both bind the certificate by `_tran_branch_c.cert_bind`. */
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
    long verdict;           /* out: pcx_branch's result, enum branch_verdict */
    /* THE CERTIFICATE (`_BranchCheck._branch_certify`, the run's data from
       `_branch_cert_data`): where the cheap level fired, rank C >= r0 proved
       without the SVD -- B = C[ci, cj] (r0 x r0), block diagonal by `bstart`
       (every entry between two blocks tested zero), each block's
       1 / ||B_k^-1||_F above tol (1 + delta_k).  r0 = 0: none. */
    long r0, nb;            /* B's size, its blocks */
    const int64_t *ci, *cj; /* r0 each: C's rows and columns, block by block */
    const int64_t *bstart;  /* nb + 1: each block's first index into ci, cj */
    const int64_t *bid;     /* r0: each index's block */
    double growth;          /* delta_k's factor (BRANCH_CERT_GROWTH) */
    double *work;           /* 2 kmax^2: a block and its inverse */
    int64_t *ipiv;          /* kmax: dgesv's pivots */
    void *dgesv;            /* numpy's own `scipy_dgesv_64_` (dgesv64_t) */
    long ncert;             /* out, running: the screens it proved (its
                               callers zero it) */
    double bound;           /* out: the smallest block's bound, where proved */
} branch_t;
long pcx_branch(branch_t *s);

/* pcx_branch's result, left in `verdict` too: `_branch_screen_C`'s cheap level */
enum branch_verdict {
    BRANCH_QUIET = 0,            /* nothing collapsed, or the certificate proved
                                    the rank: no SVD */
    BRANCH_COLLAPSED = 1,        /* C collapsed entirely: fired, no direction */
    BRANCH_SVD = 2,              /* the SVD decides (the Python's) */
    BRANCH_BAIL_NONFINITE = 3    /* an entry not finite: the Python screens */
};
