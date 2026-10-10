/* sources.c -- the independent sources' waveforms in C (the C stepping loop
   plan's stage 1, 2026-10-10): `func.py`'s `f` and `next_event` of each
   waveform the stamp plan's direct source pass serves, in the Python's own
   operations and order, and that pass's sum and scatter
   (`_stamp_plan._SourceDirect.run`): per source v = a + f(t) (a the DC
   value), in element order into bins from +0.0 -- a voltage source's -v
   into its branch row, a current source's +v and -v into its nodes.
   DECLINE, NEVER DIFFER: a value the Python would make by another route --
   a non-finite one, one past the fast add's 1e300 guard, a sine past
   `_sin_scalar`'s guards, an exponential past exp's range -- bails, and the
   Python makes it.  exp and sin are the C library's, which `_tran_sources`
   has checked against numpy's on this CPU (`_hdl_cbackend.libm_check`). */
#include <stdint.h>
#include <math.h>
#include "sources.h"

#define BIG 1e300
/* past this, j * per is no longer the Python's (its int converted) */
#define PULSE_K_MAX 4503599627370496.0

/* Python's float % (CPython's float_rem) */
static double py_mod(double vx, double wx)
{
    double mod = fmod(vx, wx);
    if (mod) {
        if ((wx < 0) != (mod < 0)) mod += wx;
    } else {
        mod = copysign(0.0, wx);
    }
    return mod;
}

/* Python's float // (CPython's _float_div_mod) */
static double py_floordiv(double vx, double wx)
{
    double mod = fmod(vx, wx), div = (vx - mod) / wx, floordiv;
    if (mod) {
        if ((wx < 0) != (mod < 0)) div -= 1.0;
    }
    if (div) {
        floordiv = floor(div);
        if (div - floordiv > 0.5) floordiv += 1.0;
    } else {
        floordiv = copysign(0.0, vx / wx);
    }
    return floordiv;
}

/* `_sin_scalar`, its guards a bail */
static int sin_f(const double *p, double t, double *f)
{
    const double td = p[0], th = p[1], om = p[2], ph = p[3], of = p[4], am = p[5];
    double dt, a, b, q;
    if (!(-BIG < t && t < BIG)) return 0;
    dt = (t > td) ? t - td : 0.0;
    a = -th * dt;
    b = om * dt + ph;
    if (!(-700.0 <= a && a <= 700.0 && -BIG < b && b < BIG)) return 0;
    q = am * exp(a);
    if (!(-BIG < q && q < BIG)) return 0;
    *f = of + q * sin(b);
    return 1;
}

/* `Pulse.f`: the period folded by Python's %, then the `where` chain -- the
   last condition that holds wins, each branch's arithmetic as written */
static int pulse_f(const double *p, double t, double *f)
{
    const double v1 = p[0], v2 = p[1], td = p[2], tr = p[3], tf = p[4], pw = p[5], per = p[6];
    double v = v1;
    if (per != 0) t = py_mod(t, per);
    if (t < ((td + tr) + pw) + tf) v = v2 + ((v1 - v2) / tf) * (t - ((td + tr) + pw));
    if (t < (td + tr) + pw) v = v2;
    if (t < td + tr) v = v1 + ((v2 - v1) / tr) * (t - td);
    if (t < td) v = v1;
    *f = v;
    return 1;
}

/* `PWL.f`: the ends held, else the first segment holding t; none (a NaN
   time): the Python's None, which raises -- a bail */
static int pwl_f(const double *p, long np_, double t, double *f, long *nan)
{
    const double *T = p, *V = p + np_;
    long i;
    if (t <= T[0]) { *f = V[0]; return 1; }
    if (t >= T[np_ - 1]) { *f = V[np_ - 1]; return 1; }
    for (i = 0; i < np_ - 1; i++) {
        if (T[i] <= t && t <= T[i + 1]) {
            if (T[i + 1] == T[i]) { *f = V[i + 1]; return 1; }
            *f = V[i] + ((V[i + 1] - V[i]) * (t - T[i])) / (T[i + 1] - T[i]);
            return 1;
        }
    }
    *nan = 1;
    return 0;
}

/* `Exp.f` (the taus nonzero: checked at the table's build) */
static int exp_f(const double *p, double t, double *f)
{
    const double v1 = p[0], v2 = p[1], td1 = p[2], tau1 = p[3], td2 = p[4], tau2 = p[5];
    double e, vt;
    if (t <= td1) { *f = v1; return 1; }
    if (t <= td2) {
        e = exp((-(t - td1)) / tau1);
        *f = v1 + (v2 - v1) * (1.0 - e);
        return isfinite(e) && isfinite(*f);
    }
    e = exp((-(td2 - td1)) / tau1);
    vt = v1 + (v2 - v1) * (1.0 - e);
    if (!(isfinite(e) && isfinite(vt))) return 0;
    e = exp((-(t - td2)) / tau2);
    *f = vt + (v1 - vt) * (1.0 - e);
    return isfinite(e) && isfinite(*f);
}

/* `AM.f` */
static int am_f(const double *p, double t, double *f)
{
    const double vo = p[0], va = p[1], fc = p[2], fm = p[3], m = p[4], pi = p[5];
    const double am = ((2.0 * pi) * fm) * t, ac = ((2.0 * pi) * fc) * t;
    double mod;
    if (!(isfinite(am) && isfinite(ac))) return 0;
    mod = 1.0 + m * sin(am);
    *f = vo + (va * mod) * sin(ac);
    return isfinite(mod) && isfinite(*f);
}

/* `SFFM.f` */
static int sffm_f(const double *p, double t, double *f)
{
    const double vo = p[0], va = p[1], fc = p[2], mdi = p[3], fm = p[4], pi = p[5];
    const double am = ((2.0 * pi) * fm) * t;
    double arg;
    if (!isfinite(am)) return 0;
    arg = ((2.0 * pi) * fc) * t + mdi * sin(am);
    if (!isfinite(arg)) return 0;
    *f = vo + va * sin(arg);
    return isfinite(*f);
}

/* the earliest of `Pulse.next_event`'s candidates strictly after t: the
   period index by Python's //, each candidate j * per + edge */
static int pulse_next(const double *p, double t, double *tn)
{
    const double td = p[2], tr = p[3], tf = p[4], pw = p[5], per = p[6];
    const double edges[5] = {td, td + tr, (td + tr) + pw, ((td + tr) + pw) + tf, per};
    double best = INFINITY, k;
    int j, e, ne = (per != 0) ? 5 : 4;
    if (per != 0) {
        k = py_floordiv(t, per);
        if (!(fabs(k) < PULSE_K_MAX)) return 0;
        for (j = -1; j <= 1; j++)
            for (e = 0; e < ne; e++) {
                double c = (k + j) * per + edges[e];
                if (c > t && c < best) best = c;
            }
    } else {
        for (e = 0; e < ne; e++)
            if (edges[e] > t && edges[e] < best) best = edges[e];
    }
    *tn = best;
    return 1;
}

static int next_event(long kind, const double *p, long np_, double t, double *tn)
{
    long i;
    switch (kind) {
    case SOURCES_SIN:
        *tn = (t < p[0]) ? p[0] : INFINITY;
        return 1;
    case SOURCES_PULSE:
        return pulse_next(p, t, tn);
    case SOURCES_PWL:
        for (i = 0; i < np_; i++)
            if (p[i] > t) { *tn = p[i]; return 1; }
        *tn = INFINITY;
        return 1;
    case SOURCES_EXP:
        *tn = (t < p[2]) ? p[2] : (t < p[4]) ? p[4] : INFINITY;
        return 1;
    default:
        *tn = INFINITY;
        return 1;
    }
}

long hdl_fn(sources_t *s)
{
    const double t = s->t;
    long k;
    s->bad = -1;
    if (s->mode == SOURCES_NEXT) {
        double best = INFINITY;
        for (k = 0; k < s->ns; k++) {
            double tn;
            if (!next_event(s->kind[k], s->par + s->poff[k], s->npts[k], t, &tn)) {
                s->bad = k;
                return s->status = SOURCES_BAIL_VALUE;
            }
            if (tn < best) best = tn;
        }
        s->t_next = best;
        return s->status = SOURCES_DONE;
    }
    for (k = 0; k < s->n; k++) s->out[k] = 0.0;
    for (k = 0; k < s->ns; k++) {
        const double *p = s->par + s->poff[k], a = s->a[k];
        double b = 0.0, v;
        long nan = 0;
        int ok;
        switch (s->kind[k]) {
        case SOURCES_DC: ok = 1; break;
        case SOURCES_SIN: ok = sin_f(p, t, &b); break;
        case SOURCES_PULSE: ok = pulse_f(p, t, &b); break;
        case SOURCES_PWL: ok = pwl_f(p, s->npts[k], t, &b, &nan); break;
        case SOURCES_EXP: ok = exp_f(p, t, &b); break;
        case SOURCES_AM: ok = am_f(p, t, &b); break;
        case SOURCES_SFFM: ok = sffm_f(p, t, &b); break;
        default: ok = 0;
        }
        /* the fast add's guard, both values: the Python adds them so */
        if (!ok || !(-BIG < a && a < BIG && -BIG < b && b < BIG)) {
            s->bad = k;
            return s->status = nan ? SOURCES_BAIL_PWL : SOURCES_BAIL_VALUE;
        }
        v = a + b;
        if (s->isv[k]) {
            s->out[s->row2[k]] += -v;
        } else {
            s->out[s->row0[k]] += v;
            s->out[s->row1[k]] += -v;
        }
    }
    return s->status = SOURCES_DONE;
}
