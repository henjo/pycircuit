#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "pcx_abi.h"
/* the units' structs -- each header declares its own unit's entry as
   `hdl_fn`, renamed here: this unit calls them through the pointers it is
   handed, as the Newton calls the core */
#define hdl_fn pcx_newton_entry_
#include "newton.h"
#undef hdl_fn
#define hdl_fn pcx_pred_entry_
#include "pred.h"
#undef hdl_fn
#define hdl_fn pcx_branch_entry_
#include "branch.h"
#undef hdl_fn
#define hdl_fn pcx_sources_entry_
#include "sources.h"
#undef hdl_fn
#define hdl_fn pcx_lte_entry_
#include "lte.h"
#undef hdl_fn
#define hdl_fn pcx_ctrl_entry_
#include "ctrl.h"
#undef hdl_fn
#include "loop.h"

#define CALL(fn, T, p) (((long (*)(T *)) (fn))(p))
/* `integrator.SHRINK_GUARD_RATIO` (test_loop_c holds the two together) */
#define SHRINK_GUARD_RATIO 0.1
#define COPY(dst, src, n) memcpy((dst), (src), (size_t) (n) * sizeof(double))

/* `_next_breakpoint` on the table: `cir.next_event` is numpy's maximum of t
   and the sources' earliest event after it (every other element's answer
   is inf: checked at the loop's entry); a corner closer than `minbreak`
   (relative) is passed over, looking again `minbreak * 1e3` on, as
   `_StepEvents._next_breakpoint` does.  `m` is Python's `max(abs(t), 1.0)`
   (the first argument kept unless the second is strictly greater); `nb` is
   numpy's `maximum(t, t_next)`, `t >= t_next ? t : t_next` for a finite t
   (a NaN `t_next` propagates, as numpy's does).  0: the table handed the
   time back, nothing written to `out`. */
static int next_breakpoint(loop_t *s, double t, double *out)
{
    sources_t *src = (sources_t *) s->src;
    double a = fabs(t), m = (1.0 > a) ? 1.0 : a, nb;
    src->mode = SOURCES_NEXT;
    src->t = t;
    if (CALL(s->src_fn, sources_t, src) != SOURCES_DONE) return 0;
    nb = (t >= src->t_next) ? t : src->t_next;
    if (nb <= t + s->minbreak * m) {
        double t2 = t + (s->minbreak * 1e3) * m;
        src->mode = SOURCES_NEXT;
        src->t = t2;
        if (CALL(s->src_fn, sources_t, src) != SOURCES_DONE) return 0;
        nb = (t2 >= src->t_next) ? t2 : src->t_next;
    }
    *out = nb;
    return 1;
}

/* What every attempt the loop makes leaves -- accepted or rejected, as the
   Python's `_tran_newton_c` write-back and error test leave it: the
   converged point, its charge and current, its C and Geq (the buffer `nxt`
   the Newton wrote), the step and companion it was made with, the screen's
   maximum, the error test's running reference and error.  The Newton's
   `xc`, `qj`, `iqj` hold only until its next call: copied here. */
static void attempt_made(loop_t *s, newton_t *nt, long nxt, long formula, double h, double h_last,
                         double tn, int screened, double ref, int judged, double err)
{
    long n = s->n;
    COPY(s->xatt, nt->xc, n);
    COPY(s->qatt, nt->qj, n);
    COPY(s->iqatt, nt->iqj, n);
    s->cur = nxt;
    s->att_h = h;
    s->att_h_last = h_last;
    s->att_t = tn;
    s->att_formula = formula;
    if (screened) {
        s->cmax = ref;
        s->screens++;
    }
    if (judged) {
        lte_t *le = (lte_t *) s->lte;
        if (s->keeps_run) {
            COPY(s->run, le->out, n);
            s->have_run = 1;
        }
        s->last_err = err;
        s->ltes++;
        s->lte_euler += formula == PCX_FORMULA_EULER;
    }
    s->iters += nt->iters;
    if (s->hk >= 2) s->pred_fast++;
    s->attempts++;
}

/* `_SteppingLoop`'s attempts -- `where` / `take` / `judge` / `accept` -- for a
   multistep family whose every piece it was checked to stand in for, fixed
   or adaptive: the steps taken into X, T and the state; it stops BEFORE an
   attempt it cannot make as the Python would, that attempt's state
   untouched.  On any status the state is the last MADE attempt's: `cur`
   names the C and Geq buffer it wrote (the other may hold a discarded
   attempt's), the `att_*` fields its step, and the rings, the predictor's
   nodes, `xlast`, X and T accepted steps only.  The counts are this call's
   (zeroed here). */
long hdl_fn(loop_t *s)
{
    newton_t *nt = (newton_t *) s->newton;
    pred_t *pr = (pred_t *) s->pred;
    sources_t *src = (sources_t *) s->src;
    branch_t *br = (branch_t *) s->branch;
    long n = s->n, r;
    s->k = 0;
    s->attempts = s->rejected = s->iters = s->euler_steps = s->breaks = s->screens = 0;
    s->pred_fast = s->ltes = s->lte_euler = 0;
    for (;;) {
        double t = s->t, h = s->h, tn, h_last, t_break, a0 = 0.0, a1 = 0.0, a2 = 0.0;
        double h_next, err = 0.0;
        long formula, st, nxt;
        int is_first, landing, on_break;
        const double *seed;
        /* `while self.t < self.tend` */
        if (!(t < s->tend)) return s->status = LOOP_DONE;
        if (s->k >= s->cap) return s->status = LOOP_FULL;
        /* 1. where: the order drop a landing or a force asks for, the
           breakpoint -- a grid crosses it, an adaptive step lands on it --
           and the run's end */
        is_first = s->is_first || s->landing || s->order_drop;
        if (!next_breakpoint(s, t, &t_break)) return s->status = LOOP_BAIL_SOURCES;
        if (s->adaptive) {
            if (t + h > t_break) {
                h = t_break - t;
                landing = 1;
            } else {
                landing = 0;
            }
            if (t + h > s->tend) h = s->tend - t;
        } else {
            /* (`_SteppingLoop.where`'s tolerances, both of them) */
            landing = t_break <= t + h * (1.0 + 1e-9);
            if (t + h > s->tend) {
                h = s->tend - t;
                /* (a uniform grid's rounding residue: the run is over -- the
                   state left as the last accept left it, `where` the
                   Python's) */
                if (h <= 1e-9 * s->timestep) return s->status = LOOP_DONE;
            }
        }
        /* 2. take: the companion (`check_order_drop`, `bdf2_alphas`; its
           shrink guard `SHRINK_GUARD_RATIO`, its growth limit `zs_ratio`,
           `ZERO_STABILITY_RATIO`) */
        tn = t + h;
        h_last = s->have_dt_last ? s->dt_last : h;
        formula = s->base;
        if (s->base == PCX_FORMULA_GEAR2) {
            if (is_first || h / h_last > s->zs_ratio
                    || (h / h_last < SHRINK_GUARD_RATIO && s->shrink_guard))
                formula = PCX_FORMULA_EULER;
        } else if (s->base == PCX_FORMULA_TRAPEZOIDAL) {
            if (is_first) formula = PCX_FORMULA_EULER;
        } else if (s->base != PCX_FORMULA_EULER) {
            return s->status = LOOP_BAIL_FORMULA;
        }
        if (formula == PCX_FORMULA_GEAR2) {
            a0 = (2.0 * h + h_last) / (h * (h + h_last));
            a1 = -(h + h_last) / (h * h_last);
            a2 = h / (h_last * (h + h_last));
        }
        /* the seed: the predictor's, or the last point */
        pr->hk = s->hk;
        pr->ttarget = tn;
        st = CALL(s->pred_fn, pred_t, pr);
        if (st == PRED_DONE) seed = pr->out;
        else if (st == PRED_NONE) seed = s->xlast;
        else return s->status = LOOP_BAIL_PRED;
        for (r = 0; r < n; r++) if (!isfinite(seed[r])) return s->status = LOOP_BAIL_NEWTON;
        /* u at the attempt's end */
        src->mode = SOURCES_U;
        src->t = tn;
        if (CALL(s->src_fn, sources_t, src) != SOURCES_DONE) return s->status = LOOP_BAIL_SOURCES;
        /* the Newton, its inputs as `_tran_newton_c.solve` sets them */
        nt->formula = formula;
        nt->a0 = a0;
        nt->a1 = a1;
        nt->a2 = a2;
        nt->h = h;
        nt->theta = 0.0;
        nt->q1 = s->qring;
        nt->q2 = (formula == PCX_FORMULA_GEAR2) ? s->qring + n : s->qring;
        nt->iq1 = (formula == PCX_FORMULA_TRAPEZOIDAL || formula == PCX_FORMULA_THETA)
                  ? s->iqring : s->qring;
        nt->u = src->out;
        nt->xseed = seed;
        nt->have_cache = s->cur >= 0;
        nxt = (s->cur == 0) ? 1 : 0;
        if (s->cur >= 0) {
            nt->cx = s->xatt;
            nt->cC = s->cC[s->cur];
        }
        nt->Cj = s->cC[nxt];
        nt->Geqj = s->Geq[nxt];
        if (CALL(s->newton_fn, newton_t, nt) != NEWTON_CONVERGED) return s->status = LOOP_BAIL_NEWTON;
        /* the branch screen, after the solve */
        if (br) {
            br->C = nt->Cj;
            br->cmax = s->cmax;
            if (CALL(s->branch_fn, branch_t, br) != BRANCH_QUIET) return s->status = LOOP_BAIL_SCREEN;
        }
        /* 3. judge: the error test, the controller's law, the budget */
        if (s->adaptive) {
            lte_t *le = (lte_t *) s->lte;
            ctrl_t *ct = (ctrl_t *) s->ctrl;
            long had_retry = ct->have_retry, kept_point = ct->retry_point, had_q = ct->have_q;
            double kept_h = ct->retry_h, kept_err = ct->retry_err, kept_q = ct->q_carried;
            double coef, hl = h_last, hl2 = 0.0, s12 = 0.0, s23 = 0.0, s123 = 0.0;
            int euler = formula == PCX_FORMULA_EULER;
            /* the estimate's scalars, as `_tran_lte_c` makes them in Python
               floats: Euler `-0.5 * (2.0 * h / (h + hl) if hl else 1.0)`, Gear
               `-h * (h + hl)` (the unary minus first), the trapezoid `-(h ** 2)`
               -- Python's `**` on floats is libm `pow`, which can differ from
               `h * h` in the last bit, so `pow` it is (`-fno-builtin-pow`
               keeps the call) */
            if (euler) {
                coef = -0.5 * ((hl != 0.0) ? (2.0 * h) / (h + hl) : 1.0);
                hl = 0.0;
            } else {
                if (!s->have_dt_last2) return s->status = LOOP_BAIL_LTE;
                hl2 = s->dt_last2;
                coef = (formula == PCX_FORMULA_GEAR2) ? (-h) * (h + hl) : -pow(h, 2.0);
                s12 = h + hl;
                s23 = hl + hl2;
                s123 = s12 + hl2;
            }
            if (!(isfinite(coef) && isfinite(s12) && isfinite(s23) && isfinite(s123)))
                return s->status = LOOP_BAIL_LTE;
            le->formula = euler ? LTE_FORMULA_EULER : LTE_FORMULA_DD3;
            le->h1 = h;
            le->h2 = hl;
            le->h3 = hl2;
            le->coef = coef;
            le->s12 = s12;
            le->s23 = s23;
            le->s123 = s123;
            le->xc = nt->xc;
            le->xl = s->xlast;
            le->q = nt->qj;
            le->q1 = s->qring;
            le->q2 = s->qring + n;
            le->q3 = s->qring + 2 * n;
            le->iq1 = s->iqring;
            le->J = nt->Jj;
            le->run = s->run;
            le->have_run = s->have_run;
            if (CALL(s->lte_fn, lte_t, le) != LTE_DONE) return s->status = LOOP_BAIL_LTE;
            err = le->err;
            ct->err = err;
            /* (`_max_error`'s exponent: the order + 1; the point the attempt
               left from as an index, where the Python keeps the object) */
            ct->p = euler ? 2.0 : 3.0;
            ct->h = h;
            ct->max_step = s->max_step;
            ct->point = s->npoint - 1;
            CALL(s->ctrl_fn, ctrl_t, ct);
            if (!ct->accepted) {
                double hn = ct->h_next;
                int growth = hn > h;
                if (s->point_retries < s->max_retries
                        && (growth || (s->rejects < s->max_reject && h > s->minstep))) {
                    /* retried at the point: the attempt stands as made */
                    attempt_made(s, nt, nxt, formula, h, h_last, tn, br != NULL,
                                 br ? br->ref : 0.0, 1, err);
                    s->rejected++;
                    s->point_retries++;
                    s->rejects += !growth;
                    s->landing = landing;
                    s->is_first = is_first;
                    s->order_drop = 0;
                    s->att_rejected = 1;
                    s->t_break = t_break;
                    /* (`_LMMSteps.judge`'s `max(h, minstep)`) */
                    s->h = (s->minstep > hn) ? s->minstep : hn;
                    continue;
                }
                /* the budget spent: the force-accept (its warning's record)
                   is the Python's -- the law's record as before */
                ct->have_retry = had_retry;
                ct->retry_point = kept_point;
                ct->retry_h = kept_h;
                ct->retry_err = kept_err;
                ct->have_q = had_q;
                ct->q_carried = kept_q;
                return s->status = LOOP_BAIL_FORCE;
            }
            s->rejects = s->point_retries = 0;
            h_next = ct->h_next;
        } else {
            h_next = s->timestep;
        }
        /* 4. events: none -- a circuit with a declared state event declines
           at the loop's entry, so no step is cut here */
        /* 5. accept (`on_break`'s tolerance is `accept`'s; the statistics'
           extremes set first, then replaced only when strictly beyond) */
        attempt_made(s, nt, nxt, formula, h, h_last, tn, br != NULL, br ? br->ref : 0.0,
                     (int) s->adaptive, err);
        s->att_rejected = 0;
        on_break = t + h >= t_break * (1.0 - 1e-12);
        s->landing = landing || on_break;
        s->is_first = 0;
        s->order_drop = 0;
        s->t = tn;
        COPY(s->X + s->k * n, nt->xc, n);
        s->T[s->k] = tn;
        if (!s->have_hm || h < s->hmin) s->hmin = h;
        if (!s->have_hm || h > s->hmax) s->hmax = h;
        s->have_hm = 1;
        s->breaks += on_break;
        s->euler_steps += (formula == PCX_FORMULA_EULER && s->base != PCX_FORMULA_EULER);
        /* the rings (`_ring_push`): the newest first, the oldest dropped */
        for (r = s->hist_len - 1; r > 0; r--) {
            COPY(s->qring + r * n, s->qring + (r - 1) * n, n);
            COPY(s->iqring + r * n, s->iqring + (r - 1) * n, n);
        }
        COPY(s->qring, nt->qj, n);
        COPY(s->iqring, nt->iqj, n);
        /* the predictor's node (`_pred_note`: appended past every node, so
           its sort keeps it last and its dedupe -- `1e-14 * max(abs(t), 1.0)`,
           `_StagePredictor`'s -- compares it with the one before alone; the
           newest eight kept: the capacity pred.c's `nodes[8]` and the glue's
           `ht` / `hx` buffers share) */
        {
            double a = fabs(tn), m = (1.0 > a) ? 1.0 : a;
            if (s->hk && fabs(tn - s->ht[s->hk - 1]) <= 1e-14 * m) {
                s->hk--;
            } else if (s->hk == PRED_HIST_MAX) {
                memmove(s->ht, s->ht + 1, (PRED_HIST_MAX - 1) * sizeof(double));
                memmove(s->hx, s->hx + n, (size_t) ((PRED_HIST_MAX - 1) * n) * sizeof(double));
                s->hk = PRED_HIST_MAX - 1;
            }
            s->ht[s->hk] = tn;
            COPY(s->hx + s->hk * n, nt->xc, n);
            s->hk++;
        }
        COPY(s->xlast, nt->xc, n);
        s->have_dt_last2 = s->have_dt_last;
        s->dt_last2 = s->dt_last;
        s->dt_last = h;
        s->have_dt_last = 1;
        s->t_break = t_break;
        /* (`accept`'s `min(h_next, max_step)`: the first argument kept unless
           the second is strictly smaller) */
        s->h = s->adaptive ? ((s->max_step < h_next) ? s->max_step : h_next) : s->timestep;
        s->npoint++;
        s->k++;
    }
}
