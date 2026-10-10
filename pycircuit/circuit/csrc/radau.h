/* radau.h -- the Radau IIA stage solve's state and entry: the interface `_tran_radau_c` loads
 * radau.c through.  cffi lays the struct out from this text and
 * radau.c includes it, so the compiler checks the entry against it. */
typedef struct {
    const void *core; void *core_fn; double T;
    const double *binG, *bini; double *dummy;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv;
    long n, iref, maxiter, maxa;
    double h, reltol, abstol;
    double A[9];
    const double *qn, *u, *seed;
    double *Ya, *qa, *ia, *Ca, *Ga;
    double *K, *F, *R, *Jf, *b, *blk, *stp, *tolv;
    int64_t *ipiv;
    long assemblies, iters, status;
} radau_t;
long hdl_fn(radau_t *s);

/* hdl_fn's result, left in `status` too (`_tran_radau_c.BAIL` names the
   bails from these; 7 and 10 unused): */
enum radau_status {
    RADAU_RUNNING = 0,
    RADAU_CONVERGED = 1,
    RADAU_BAIL_NONFINITE = 2,    /* an input or a stage not finite */
    RADAU_BAIL_FLAGS = 3,        /* a floating-point flag numpy would raise on */
    RADAU_BAIL_RBOUND = 4,       /* the residual or the Jacobian out of bounds */
    RADAU_BAIL_DGESV = 5,        /* the stage system's solve: singular */
    RADAU_BAIL_ASSEMBLIES = 6,   /* more assemblies than the buffers hold */
    RADAU_BAIL_WALKSTOP = 8,     /* a limiter declined: the walk stopped */
    RADAU_BAIL_STEP = 9,         /* a step not finite */
    RADAU_BAIL_MAXITER = 11      /* no convergence in maxiter */
};
