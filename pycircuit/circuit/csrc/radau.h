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
