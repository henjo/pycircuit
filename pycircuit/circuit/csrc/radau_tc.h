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
