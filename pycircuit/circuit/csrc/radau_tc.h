/* radau_tc.h -- the transformed Radau stage solve's state and entry: the interface `_tran_radau_tc` loads
 * radau_tc.c through.  cffi lays the struct out from this text and
 * radau_tc.c includes it, so the compiler checks the entry against it. */
typedef struct {
    const void *core; void *core_fn; double T;
    const double *bini; double *dummy; long bits;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv; void *dgetrs; void *klu_refactor; void *klu_solve;
    void *Symbolic; void *Numeric; void *Common;
    int *Ap; int *Ai; double *Ax;
    long n, iref, maxiter, fused;
    double h, reltol, abstol, rtol;
    double A[9];
    double p[12];
    double v[9];
    const double *qn, *u, *seed;
    double *lu_a; int64_t *lu_ipiv; int64_t lu_info;
    long klu_fresh, refactored;
    double *Y, *Yp, *q, *iv, *K, *F, *R, *b0, *x0, *r1, *x1, *y1, *ab, *dY, *stp, *tolv;
    long iters, status;
    const double *binG; long end_bits, end_cg, end_done;
    double *e_q, *e_i, *e_C, *e_G, *e_C01, *e_G01, *e_q01;
} radau_tf_t;
long hdl_fn(radau_tf_t *s);

/* hdl_fn's result, left in `status` too (`_tran_radau_tc.BAIL` names the
   bails from these): */
enum radau_tc_status {
    RADAU_TC_RUNNING = 0,
    RADAU_TC_CONVERGED = 1,
    RADAU_TC_BAIL_NONFINITE = 2, /* an input or a stage not finite */
    RADAU_TC_BAIL_FLAGS = 3,     /* a floating-point flag numpy would raise on */
    RADAU_TC_BAIL_RHS = 4,       /* a right-hand side not finite */
    RADAU_TC_BAIL_LU = 5,        /* the real block's factorisation: singular */
    RADAU_TC_BAIL_REFACTOR = 6,  /* KLU's refactorisation failed */
    RADAU_TC_BAIL_KLUSOLVE = 7,  /* KLU's solve failed */
    RADAU_TC_BAIL_RESIDUAL = 8,  /* the complex block's residual check failed */
    RADAU_TC_BAIL_STEP = 9,      /* a step not finite */
    RADAU_TC_BAIL_WALKSTOP = 10, /* a limiter declined: the walk stopped */
    RADAU_TC_BAIL_MAXITER = 11   /* no convergence in maxiter */
};
