/* pred.h -- the stage predictor's fast path: its state and entry, the
 * interface `_tran_pred_c` loads pred.c through.  cffi lays the struct out
 * from this text and pred.c includes it, so the compiler checks the entry
 * against it.
 *
 * WHO SETS WHAT.  `_tran_pred_c.new_state` makes the struct: `n`, the
 * buffers `out`, `xt`, `motion` (numpy arrays it keeps with the struct),
 * `dgesv`, `dgemv`, and `wk` = 0.  Then either `_tran_pred_c.predict` (the
 * twin's tests: a fresh struct each call) sets `hk`, `ht`, `hx`, `ttarget`,
 * `deg`, `clamp`; or, in the C loop, `_tran_loop_c._bind` points `ht` and
 * `hx` at the loop's own history buffers and sets `deg` and `clamp` once a
 * run, and loop.c sets `hk` and `ttarget` at each attempt (it keeps the
 * history's contents: the accepted points, the newest eight).
 *
 * CAPACITIES, which the glue assumes too.  The fit takes at most 5 nodes:
 * `deg` <= 4 (`_tran_pred_c.MAX_DEG`) -- `tau`, `w`, `ipiv` hold 5, `V` 25
 * and `xt` 5 x n (`new_state`'s `np.empty(5 * n)`).  The history holds at
 * most 8 nodes: loop.c keeps the newest eight in buffers of 8 and 8 x n
 * (`_tran_loop_c._bind`'s `np.empty(8)`, `np.empty(8 * n)`), `predict`
 * declines a longer one, and pred.c records the first eight nodes it keeps
 * (`nodes[PRED_HIST_MAX]`). */
typedef struct {
    long n;                 /* unknowns */
    long hk;                /* the history's nodes, oldest first (`_pred_hist`) */
    const double *ht;       /* their times, hk */
    const double *hx;       /* their states, hk x n row-major */
    double ttarget;         /* the time predicted for: the attempt's end */
    long deg;               /* the fit's degree (`_pred_degree`), <= 4 */
    double clamp;           /* PRED_CLAMP */
    double *out;            /* n: the prediction (read on PRED_DONE only) */
    double *xt;             /* work: the nodes taken, newest first, 5 x n */
    double *motion;         /* work: n */
    /* the fit, and its weights carried between calls: `w` solved for the
       normalised times `tau` (the first `wk` of each valid), reused when
       a call's `tau` has the same bits -- `_fit_fast`'s weight memo, one
       entry deep (the weights are a pure function of tau's bits); `V` the
       Vandermonde, work (dgesv leaves its LU there) */
    double tau[5], w[5], V[25];   /* the fit; `w` kept for the last `tau` */
    long wk;                /* the weights `w` holds (0: none yet) */
    int64_t ipiv[5];        /* work: dgesv's pivots */
    void *dgesv, *dgemv;    /* numpy's own (`scipy_dgesv_64_`, `scipy_cblas_dgemv64_`) */
    long status;            /* OUT: pcx_pred's result, enum pred_status */
} pred_t;
long pcx_pred(pred_t *s);

/* the capacities the struct, pred.c, the C loop and both glues share: the
   nodes a history holds (`_StagePredictor.PRED_HIST_MAX`: `ht`, `hx`) and the
   fit's (`deg` + 1 at most: `tau`, `w`, `ipiv`, `V`'s rows, `xt`'s) -- the
   arrays above are sized by them (tests/test_pred_c holds them together) */
enum pred_cap {
    PRED_HIST_MAX = 8,
    PRED_FIT_MAX = 5
};

/* pcx_pred's result, left in `status` too */
enum pred_status {
    PRED_DONE = 1,          /* `out` is the prediction */
    PRED_NONE = 2,          /* no prediction: the seed is the last point */
    PRED_GENERAL = 3,       /* the general path's (a history not strictly
                               ascending, a tie, one at or past the target) */
    PRED_BAIL_FLAGS = 4     /* a floating-point flag numpy would report */
};
