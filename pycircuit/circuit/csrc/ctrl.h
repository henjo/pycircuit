/* ctrl.h -- the step controller's law: its state and entry, the interface
 * `_tran_ctrl_c` loads ctrl.c through.  cffi lays the struct out from this
 * text and ctrl.c includes it, so the compiler checks the entry against it.
 *
 * WHO SETS WHAT.  `_tran_ctrl_c.new_state` sets `mode` and the constants
 * once a struct.  In the C loop, loop.c sets the attempt's fields before
 * each call and reads `accepted` and `h_next`; `_tran_loop_c` loads the
 * retry state from the controller at each loop entry and writes it back
 * after.  `_tran_ctrl_c.evaluate` (the twin's tests) does all of it. */
typedef struct {
    /* the attempt, each call */
    double err, p;              /* the error test's normalised error and order + 1 */
    double h, max_step;         /* the attempt's step, the run's cap */
    long point;                 /* the accepted point the attempt left from: the
                                   Python compares the object (`is`), this its index */
    /* set once (`new_state`) */
    long mode;                  /* `retry_order`: CTRL_ASSUMED / MEASURED / CARRIED */
    double min_shrink, max_growth, order_min;   /* stepcontroller's constants
                                   MIN_SHRINK_RATIO, MAX_GROWTH_RATIO, RETRY_ORDER_MIN */
    /* the controller's retry state (`_retry`, `_q_carried`), in and out:
       carried from call to call */
    long have_retry, retry_point;   /* a rejection kept (`_retry` not None) and
                                   its point's index (`_tran_loop_c` loads -2
                                   for a point that is not the loop's: no
                                   `point` equals it) */
    double retry_h, retry_err;  /* the kept rejection's step and error */
    long have_q;                /* `_q_carried` not None */
    double q_carried;           /* the last order measured (`carried` reuses it) */
    /* out */
    long accepted;              /* 1: accepted, 0: rejected */
    double h_next;              /* the next step, as `_evaluate` returns it */
    long status;                /* hdl_fn's result: CTRL_DONE */
} ctrl_t;
long hdl_fn(ctrl_t *s);

/* `retry_order` (stepcontroller.RETRY_ORDER_MODES) */
enum ctrl_mode {
    CTRL_ASSUMED = 0,
    CTRL_MEASURED = 1,
    CTRL_CARRIED = 2
};
/* hdl_fn's result, left in `status` too */
enum ctrl_status {
    CTRL_DONE = 1
};
