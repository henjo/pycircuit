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
#include "fixed.h"

#define CALL(fn, T, p) (((long (*)(T *)) (fn))(p))

/* `_next_breakpoint` on the table: `cir.next_event` is numpy's maximum of t
   and the sources' earliest event after it (every other element's answer
   is inf: checked at the loop's entry); a corner closer than `minbreak`
   (relative) is passed over */
static int next_breakpoint(fixed_t *s, double t, double *out)
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

/* `_SteppingLoop`'s fixed-step attempt, `where` / `take` / `accept`, for a
   multistep family whose every piece it was checked to stand in for: the
   steps taken into X, T and the state; it stops BEFORE a step it cannot
   take as the Python would, that step's state untouched */
long hdl_fn(fixed_t *s)
{
    newton_t *nt = (newton_t *) s->newton;
    pred_t *pr = (pred_t *) s->pred;
    sources_t *src = (sources_t *) s->src;
    branch_t *br = (branch_t *) s->branch;
    long n = s->n, r;
    s->k = 0;
    s->iters = s->euler_steps = s->breaks = s->screens = s->pred_fast = 0;
    for (;;) {
        double t = s->t, h = s->h, tn, h_last, t_break, a0 = 0.0, a1 = 0.0, a2 = 0.0;
        long formula, st, nxt;
        int is_first, landing, on_break;
        const double *seed;
        /* `while self.t < self.tend` */
        if (!(t < s->tend)) return s->status = FIXED_DONE;
        if (s->k >= s->cap) return s->status = FIXED_FULL;
        /* 1. where: the order drop a landing asks for, the breakpoint, the
           grid's end */
        is_first = s->is_first || s->landing;
        if (!next_breakpoint(s, t, &t_break)) return s->status = FIXED_BAIL_SOURCES;
        landing = t_break <= t + h * (1.0 + 1e-9);
        if (t + h > s->tend) {
            h = s->tend - t;
            /* (a uniform grid's rounding residue: the run is over -- the
               state left as the last accept left it, `where` the Python's) */
            if (h <= 1e-9 * s->timestep) return s->status = FIXED_DONE;
        }
        /* 2. take: the companion (`check_order_drop`, `bdf2_alphas`) */
        tn = t + h;
        h_last = s->have_dt_last ? s->dt_last : h;
        formula = s->base;
        if (s->base == PCX_FORMULA_GEAR2) {
            if (is_first || h / h_last > s->zs_ratio || (h / h_last < 0.1 && s->shrink_guard))
                formula = PCX_FORMULA_EULER;
        } else if (s->base == PCX_FORMULA_TRAPEZOIDAL) {
            if (is_first) formula = PCX_FORMULA_EULER;
        } else if (s->base != PCX_FORMULA_EULER) {
            return s->status = FIXED_BAIL_FORMULA;
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
        else return s->status = FIXED_BAIL_PRED;
        for (r = 0; r < n; r++) if (!isfinite(seed[r])) return s->status = FIXED_BAIL_NEWTON;
        /* u at the step's end */
        src->mode = SOURCES_U;
        src->t = tn;
        if (CALL(s->src_fn, sources_t, src) != SOURCES_DONE) return s->status = FIXED_BAIL_SOURCES;
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
            nt->cx = s->xlast;
            nt->cC = s->cC[s->cur];
        }
        nt->Cj = s->cC[nxt];
        nt->Geqj = s->Geq[nxt];
        if (CALL(s->newton_fn, newton_t, nt) != NEWTON_CONVERGED) return s->status = FIXED_BAIL_NEWTON;
        /* the branch screen, after the solve */
        if (br) {
            br->C = nt->Cj;
            br->cmax = s->cmax;
            if (CALL(s->branch_fn, branch_t, br) != BRANCH_QUIET) return s->status = FIXED_BAIL_SCREEN;
        }
        /* 5. accept */
        on_break = t + h >= t_break * (1.0 - 1e-12);
        s->landing = landing || on_break;
        s->is_first = 0;
        s->t = tn;
        memcpy(s->X + s->k * n, nt->xc, (size_t) n * sizeof(double));
        s->T[s->k] = tn;
        if (!s->have_hm || h < s->hmin) s->hmin = h;
        if (!s->have_hm || h > s->hmax) s->hmax = h;
        s->have_hm = 1;
        s->breaks += on_break;
        s->euler_steps += (formula == PCX_FORMULA_EULER && s->base != PCX_FORMULA_EULER);
        s->iters += nt->iters;
        if (br) {
            s->cmax = br->ref;
            s->screens++;
        }
        if (s->hk >= 2) s->pred_fast++;
        /* the rings (`_ring_push`): the newest first, the oldest dropped */
        for (r = s->hist_len - 1; r > 0; r--) {
            memcpy(s->qring + r * n, s->qring + (r - 1) * n, (size_t) n * sizeof(double));
            memcpy(s->iqring + r * n, s->iqring + (r - 1) * n, (size_t) n * sizeof(double));
        }
        memcpy(s->qring, nt->qj, (size_t) n * sizeof(double));
        memcpy(s->iqring, nt->iqj, (size_t) n * sizeof(double));
        /* the predictor's node (`_pred_note`: appended past every node, so
           its sort keeps it last and its dedupe compares it with the one
           before alone; the newest eight kept) */
        {
            double a = fabs(tn), m = (1.0 > a) ? 1.0 : a;
            if (s->hk && fabs(tn - s->ht[s->hk - 1]) <= 1e-14 * m) {
                s->hk--;
            } else if (s->hk == 8) {
                memmove(s->ht, s->ht + 1, 7 * sizeof(double));
                memmove(s->hx, s->hx + n, (size_t) (7 * n) * sizeof(double));
                s->hk = 7;
            }
            s->ht[s->hk] = tn;
            memcpy(s->hx + s->hk * n, nt->xc, (size_t) n * sizeof(double));
            s->hk++;
        }
        memcpy(s->xlast, nt->xc, (size_t) n * sizeof(double));
        s->cur = nxt;
        s->have_dt_last2 = s->have_dt_last;
        s->dt_last2 = s->dt_last;
        s->dt_last = h;
        s->have_dt_last = 1;
        s->last_formula = formula;
        s->last_h_last = h_last;
        s->t_break = t_break;
        s->h = s->timestep;
        s->k++;
    }
}
