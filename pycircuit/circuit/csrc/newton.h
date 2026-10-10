/* newton.h -- the C Newton's state and entry: the interface `_tran_newton_c` loads
 * newton.c through.  cffi lays the struct out from this text and
 * newton.c includes it, so the compiler checks the entry against it. */
typedef struct {
    const void *core; void *core_fn;
    double T; long formula; double a0, a1, a2, h, theta;
    const double *q1, *q2, *iq1, *u;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv, *dgetrf, *dgetrs, *dgemv;
    long n, iref, maxiter, chord, first_bits;
    double reltol; const double *abstol, *xtol;
    const double *Cin0; const double *xseed; const double *cx, *cC; long have_cache, hit;
    double *C0, *q0, *iq0, *Geq0, *F0, *J0;
    double *C, *q, *iq, *Geq, *F, *J;
    double *xf, *xnf, *x0f;
    double *A, *absJ;
    double *b, *Fr, *x, *xn, *xd, *ax, *Is;
    int64_t *ipiv64; int *ipiv32;
    long iters, last_eval, status;
    long fold_j; double *xc, *Cj, *qj, *iqj, *Geqj, *Jj;
} newton_t;
long hdl_fn(newton_t *s);
