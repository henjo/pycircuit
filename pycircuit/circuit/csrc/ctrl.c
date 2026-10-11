#include <stdint.h>
#include <math.h>
#include "ctrl.h"

/* `_retry_exponent`: the exponent after a rejection when the one before it
   was at the SAME point and larger -- 1/q, q the order the two tries
   measured, clamped to [order_min, p]; under CARRIED a first rejection
   takes the last order measured.  0: none (the caller's 1/p) */
static double retry_exponent(ctrl_t *s, double err, double p, double target)
{
    double q, m;
    if (s->mode == CTRL_ASSUMED) return 0.0;
    if (!s->have_retry || s->retry_point != s->point || !(s->retry_h > s->h * (1.0 + 1e-9))
            || !(s->retry_err > target) || !(err > 0.0)) {
        if (s->mode == CTRL_CARRIED && s->have_q)
            return 1.0 / ((p < s->q_carried) ? p : s->q_carried);
        return 0.0;
    }
    q = log(s->retry_err / err) / log(s->retry_h / s->h);
    m = (s->order_min > q) ? s->order_min : q;          /* max(q, RETRY_ORDER_MIN) */
    q = (p < m) ? p : m;                                 /* min(., p) */
    s->q_carried = q;
    s->have_q = 1;
    return 1.0 / q;
}

/* `IntegralController._evaluate` past its error test, with the default band
   (`_band_decision`: rejected above 1, no lower edge, no damper) and
   `_judged`: Python's float operations in their order, its `min` and `max`
   (the first argument unless the second is strictly beyond it) */
long hdl_fn(ctrl_t *s)
{
    double err = s->err, p = s->p, h = s->h, exponent = 1.0 / p, target = pow(0.9, p);
    if (err > 1.0) {
        double r = retry_exponent(s, err, p, target), v;
        v = pow(target / err, r == 0.0 ? exponent : r);
        s->h_next = h * ((v > s->min_shrink) ? v : s->min_shrink);
        s->accepted = 0;
        s->have_retry = 1;
        s->retry_point = s->point;
        s->retry_h = h;
        s->retry_err = err;
    } else {
        double e = (1e-12 > err) ? 1e-12 : err, v, g;
        v = pow(target / e, exponent);
        g = h * ((v < s->max_growth) ? v : s->max_growth);
        s->h_next = (s->max_step < g) ? s->max_step : g;
        s->accepted = 1;
        s->have_retry = 0;
    }
    return s->status = CTRL_DONE;
}
