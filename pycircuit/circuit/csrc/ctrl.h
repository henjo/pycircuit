/* ctrl.h -- the step controller's law: its state and entry, the interface
 * `_tran_ctrl_c` loads ctrl.c through.  cffi lays the struct out from this
 * text and ctrl.c includes it, so the compiler checks the entry against it. */
typedef struct {
    double err, p;              /* the error test's normalised error and order + 1 */
    double h, max_step;         /* the attempt's step, the run's cap */
    long point;                 /* the accepted point the attempt left from: the
                                   Python compares the object (`is`), this its index */
    long mode;                  /* `retry_order`: CTRL_ASSUMED / MEASURED / CARRIED */
    double min_shrink, max_growth, order_min;   /* stepcontroller's constants */
    /* the controller's retry state (`_retry`, `_q_carried`), in and out */
    long have_retry, retry_point;
    double retry_h, retry_err;
    long have_q;
    double q_carried;
    /* out */
    long accepted;
    double h_next;
    long status;
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
