/* ctrl.c -- the step controller's law in C, bit for bit
   `IntegralController._evaluate` past its error test; its glue is
   `_tran_ctrl_c`, and the C loop (loop.c) calls it after each judged
   attempt.  No floating-point flag is tested: the Python's law is Python
   floats (`**` and `math.log` libm's pow and log, which `-fno-builtin-*`
   keeps real calls), which warn about nothing -- `**` would raise on an
   overflow, which these bases (below 1 on a rejection, at most 0.9 ** p
   / 1e-12 on an accept) and exponents (in (0, 1]) cannot reach.

   PYTHON'S min AND max.  `max(a, b)` returns a unless b > a, `min(a, b)`
   a unless b < a: written here as `(b > a) ? b : a` and `(b < a) ? b : a`,
   each naming the Python form it stands for.  fmax/fmin differ on a NaN
   (they return the other argument) and do not promise which of two equal
   zeros they keep. */
#include <stdint.h>
#include <math.h>
#include "pcx_minmax.h"
#include "ctrl.h"

/* `_retry_exponent`: the exponent after a rejection when the one before it
   was at the SAME point and larger -- 1/q, q the order the two tries
   measured, clamped to [order_min, p]; under CARRIED a first rejection
   takes the last order measured.  0: none (the caller's 1/p) -- 0.0
   stands for the Python's None, which no exponent 1/q with q in
   [order_min, p] can equal.  Reads the retry state as the attempt found
   it (the caller writes it after); writes `q_carried` and `have_q` only
   where it measures, as the Python sets `_q_carried`. */
static double retry_exponent(ctrl_t *s, double err, double p, double target)
{
    double q, m;
    if (s->mode == CTRL_ASSUMED) return 0.0;
    /* (`prev[0] is not s.x_last`: the points' indices compared) */
    if (!s->have_retry || s->retry_point != s->point || !(s->retry_h > s->h * (1.0 + 1e-9))
            || !(s->retry_err > target) || !(err > 0.0)) {
        if (s->mode == CTRL_CARRIED && s->have_q)
            return 1.0 / pcx_py_min(s->q_carried, p);
        return 0.0;
    }
    /* math.log is libm's log for a positive argument: both quotients are
       > 0 past the tests above (err > 1 here, retry_err > target > 0) */
    q = log(s->retry_err / err) / log(s->retry_h / s->h);
    m = pcx_py_max(q, s->order_min);                    /* max(q, RETRY_ORDER_MIN) */
    q = pcx_py_min(m, p);                                /* min(., p) */
    s->q_carried = q;
    s->have_q = 1;
    return 1.0 / q;
}

/* `IntegralController._evaluate` past its error test, with the default band
   (`_band_decision`: rejected above 1, no lower edge, no damper) and
   `_judged`: Python's float operations in their order, its `min` and `max`
   (the first argument unless the second is strictly beyond it).
   Assumes the default band (`_tran_loop_c` checks `lte_gamma_min` 0,
   `lte_gamma_max` 1 and no `lte_eta` before the loop runs) and `err` the
   error test's (finite, not negative).  Writes `accepted`, `h_next`, the
   retry state as `_judged` leaves `_retry` (kept on a rejection, cleared
   on an accept) and, measuring, `q_carried`; always CTRL_DONE. */
long pcx_ctrl(ctrl_t *s)
{
    /* target = `_band_target(0.9, p)` with no band: `safety ** p`, Python's
       float `**` -- libm's pow */
    double err = s->err, p = s->p, h = s->h, exponent = 1.0 / p, target = pow(0.9, p);
    if (err > 1.0) {
        /* shrink: `h * max(MIN_SHRINK_RATIO, (target / e) ** (exponent if
           retry is None else retry))` */
        double r = retry_exponent(s, err, p, target), v;
        v = pow(target / err, r == 0.0 ? exponent : r);
        s->h_next = h * pcx_py_max(s->min_shrink, v);
        s->accepted = 0;
        s->have_retry = 1;
        s->retry_point = s->point;
        s->retry_h = h;
        s->retry_err = err;
    } else {
        /* grow: `h * min(MAX_GROWTH_RATIO, (target / max(e, 1e-12)) **
           exponent)`, then `min(grown, s.max_step)` */
        double e = pcx_py_max(err, 1e-12), v, g;
        v = pow(target / e, exponent);
        g = h * pcx_py_min(s->max_growth, v);
        s->h_next = pcx_py_min(g, s->max_step);
        s->accepted = 1;
        s->have_retry = 0;
    }
    return s->status = CTRL_DONE;
}
