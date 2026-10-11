/* pred.h -- the stage predictor's fast path: its state and entry, the
 * interface `_tran_pred_c` loads pred.c through.  cffi lays the struct out
 * from this text and pred.c includes it, so the compiler checks the entry
 * against it. */
typedef struct {
    long n;                 /* unknowns */
    long hk;                /* the history's nodes, oldest first (`_pred_hist`) */
    const double *ht;       /* their times, hk */
    const double *hx;       /* their states, hk x n row-major */
    double ttarget;
    long deg;               /* the fit's degree (`_pred_degree`), <= 4 */
    double clamp;           /* PRED_CLAMP */
    double *out;            /* n: the prediction */
    double *xt;             /* work: the nodes taken, newest first, 5 x n */
    double *motion;         /* work: n */
    double tau[5], w[5], V[25];   /* the fit; `w` kept for the last `tau` */
    long wk;                /* the weights `w` holds (0: none yet) */
    int64_t ipiv[5];
    void *dgesv, *dgemv;    /* numpy's own (`scipy_dgesv_64_`, `scipy_cblas_dgemv64_`) */
    long status;
} pred_t;
long hdl_fn(pred_t *s);

/* hdl_fn's result, left in `status` too */
enum pred_status {
    PRED_DONE = 1,          /* `out` is the prediction */
    PRED_NONE = 2,          /* no prediction: the seed is the last point */
    PRED_GENERAL = 3,       /* the general path's (a history not strictly
                               ascending, a tie, one at or past the target) */
    PRED_BAIL_FLAGS = 4     /* a floating-point flag numpy would report */
};
