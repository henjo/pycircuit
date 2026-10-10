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

/* hdl_fn's result, left in `status` too (`_tran_newton_c.BAIL` names the
   bails from these): */
enum newton_status {
    NEWTON_RUNNING = 0,
    NEWTON_CONVERGED = 1,
    NEWTON_BAIL_X = 2,           /* an iterate not finite */
    NEWTON_BAIL_JNONFINITE = 3,  /* the chord's Jacobian not finite */
    NEWTON_BAIL_GETRF = 4,       /* the chord's factorisation: singular */
    NEWTON_BAIL_BNONFINITE = 5,  /* a chord step's right-hand side not finite */
    NEWTON_BAIL_GETRS = 6,       /* the chord's solve failed */
    NEWTON_BAIL_DGESV = 7,       /* a full step's solve: singular */
    NEWTON_BAIL_WALKSTOP = 8,    /* a limiter declined: the walk stopped */
    NEWTON_BAIL_XNEXT = 9,       /* the next iterate not finite */
    NEWTON_BAIL_CHORDBREAK = 10, /* the chord's steps stopped shrinking */
    NEWTON_BAIL_MAXITER = 11     /* no convergence in maxiter */
};
