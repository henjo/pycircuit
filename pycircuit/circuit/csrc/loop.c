#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#include "pcx_abi.h"
#include "pcx_minmax.h"
#include "pcx_fenv.h"
#include "pcx_dense.h"
/* the units this one calls through the pointers it is handed: their
   structs and entries */
#include "newton.h"
#include "pred.h"
#include "branch.h"
#include "sources.h"
#include "lte.h"
#include "ctrl.h"
#include "lmmsens.h"
#include "pcx_klu.h"
#include "radau.h"
#include "radau_tc.h"
#include "radau_fold.h"
#include "pcx_stagestep.h"
#include "loop.h"

/* a call to a unit through its pointer, typed by the unit's own declaration:
   a struct that is not its own does not compile */
#define CALL(fn, entry, p) (((__typeof__(&entry)) (fn))(p))
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
   (a NaN `t_next` propagates, as numpy's does).  A corner within
   `minbreak` of the run's end (relative to it) is the end, `tend`, as the
   families pass it.  0: the table handed the time back, nothing written to
   `out`. */
static int next_breakpoint(loop_t *s, double t, double *out)
{
    sources_t *src = (sources_t *) s->src;
    double m = pcx_py_max(fabs(t), 1.0), nb;
    src->mode = SOURCES_NEXT;
    src->t = t;
    if (CALL(s->src_fn, pcx_sources, src) != SOURCES_DONE) return 0;
    nb = pcx_np_maximum(t, src->t_next);
    if (nb <= t + s->minbreak * m) {
        double t2 = t + (s->minbreak * 1e3) * m;
        src->mode = SOURCES_NEXT;
        src->t = t2;
        if (CALL(s->src_fn, pcx_sources, src) != SOURCES_DONE) return 0;
        nb = pcx_np_maximum(t2, src->t_next);
    }
    if (fabs(nb - s->tend) <= s->minbreak * pcx_py_max(fabs(s->tend), 1.0)) nb = s->tend;
    *out = nb;
    return 1;
}

/* What every attempt the loop makes leaves -- accepted or rejected, as the
   Python's `_tran_newton_c` write-back and error test leave it: the
   converged point, its charge and current, its C and Geq (the buffer the
   Newton wrote), the step and companion it was made with, the screen's
   maximum, and where it was judged the error test's running reference and
   error.  The Newton's `xc`, `qj`, `iqj` hold only until its next call:
   copied here. */
static void attempt_made(loop_t *s, const loop_attempt_t *a, int judged)
{
    newton_t *nt = (newton_t *) s->newton;
    branch_t *br = (branch_t *) s->branch;
    long n = s->n;
    COPY(s->xatt, nt->xc, n);
    COPY(s->qatt, nt->qj, n);
    COPY(s->iqatt, nt->iqj, n);
    s->cur = a->nxt;
    s->att_h = a->h;
    s->att_h_last = a->h_last;
    s->att_t = a->tn;
    s->att_formula = a->formula;
    if (br) {
        s->cmax = br->ref;
        s->screens++;
    }
    if (judged) {
        lte_t *le = (lte_t *) s->lte;
        if (s->keeps_run) {
            COPY(s->run, le->out, n);
            s->have_run = 1;
        }
        s->last_err = a->err;
        s->ltes++;
        s->lte_euler += a->formula == PCX_FORMULA_EULER;
    }
    s->iters += nt->iters;
    s->early_checks += a->early_checks;
    if (s->hk >= 2) s->pred_fast++;
    s->attempts++;
}

/* 1. where (`_SteppingLoop.where`): the order drop a landing or a force
   asks for, the breakpoint -- a grid crosses it, an adaptive step lands on
   it -- and the run's end.  0, or the status to stop with */
static long where(loop_t *s, loop_attempt_t *a)
{
    double t = s->t, h = s->h;
    double mb = s->minbreak * pcx_py_max(fabs(s->tend), 1.0);
    /* (an adaptive run's residue under `minbreak` at its end is rounding,
       not a step: the run is over, nothing written -- the Python's) */
    if (s->adaptive && s->tend - t <= mb) return LOOP_DONE;
    a->t = t;
    a->have_tn = 0;
    a->is_first = s->is_first || s->landing || s->order_drop;
    a->clamped = 0;
    if (!next_breakpoint(s, t, &a->t_break)) return LOOP_BAIL_SOURCES;
    if (s->adaptive) {
        if (t + h > a->t_break) {
            h = a->t_break - t;
            a->landing = 1;
        } else {
            a->landing = 0;
        }
        a->clamped = a->landing;
        if (t + h > s->tend) {
            h = s->tend - t;
            a->clamped = 1;
        /* (a step ending within `minbreak` short of the end ends at it) */
        } else if (s->tend - (t + h) <= mb) {
            h = s->tend - t;
            a->clamped = 1;
        }
    } else {
        /* (`_SteppingLoop.where`'s tolerances, both of them) */
        a->landing = a->t_break <= t + h * (1.0 + 1e-9);
        if (t + h > s->tend) {
            h = s->tend - t;
            /* (a uniform grid's rounding residue: the run is over -- the
               state left as the last accept left it, `where` the Python's) */
            if (h <= 1e-9 * s->timestep) return LOOP_DONE;
        }
    }
    a->h = h;
    return 0;
}

/* 2. take (`_SteppingLoop.take`): the companion (`check_order_drop`,
   `bdf2_alphas`; its shrink guard `SHRINK_GUARD_RATIO`, its growth limit
   `zs_ratio`, `ZERO_STABILITY_RATIO`), the seed (the predictor's, or the
   last point), u at the attempt's end, the Newton -- its inputs as
   `_tran_newton_c.solve` sets them -- and the branch screen after it.  0,
   or the status to stop with */
static long lte_setup(loop_t *s, const loop_attempt_t *a, const double *xc, const double *q,
                      const double *J);

/* THE EARLY ABORT, armed for the attempt as `_SteppingLoop.early_arm` arms
   it: history (two accepted points), a full attempt left in the point's
   budget, above minstep -- the error test set up for the attempt, its x, q
   and J the Newton's at its first iterate (a 2nd-order step with no third
   point is not armed: `lte_setup`) */
static void arm_early(loop_t *s, const loop_attempt_t *a, newton_t *nt)
{
    nt->early_lte = 0;
    if (s->early_th > 0.0 && s->adaptive && s->npoint >= 2 && s->rejects + 1 < s->max_reject
            && s->point_retries + 1 < s->max_retries && a->h > s->minstep
            && lte_setup(s, a, s->xlast, s->qring, 0) == 0) {
        nt->early_lte = s->lte;
        nt->early_lte_fn = s->lte_fn;
        nt->early_th = s->early_th;
    }
}

static long take(loop_t *s, loop_attempt_t *a)
{
    newton_t *nt = (newton_t *) s->newton;
    pred_t *pr = (pred_t *) s->pred;
    sources_t *src = (sources_t *) s->src;
    branch_t *br = (branch_t *) s->branch;
    long n = s->n, r, st;
    double h = a->h, h_last;
    const double *seed;
    a->tn = a->have_tn ? a->tn : a->t + h;
    h_last = s->have_dt_last ? s->dt_last : h;
    a->h_last = h_last;
    a->formula = s->base;
    if (s->base == PCX_FORMULA_GEAR2) {
        if (a->is_first || h / h_last > s->zs_ratio
                || (h / h_last < SHRINK_GUARD_RATIO && s->shrink_guard))
            a->formula = PCX_FORMULA_EULER;
    } else if (s->base == PCX_FORMULA_TRAPEZOIDAL) {
        if (a->is_first) a->formula = PCX_FORMULA_EULER;
    } else if (s->base != PCX_FORMULA_EULER) {
        return LOOP_BAIL_FORMULA;
    }
    a->a0 = a->a1 = a->a2 = 0.0;
    if (a->formula == PCX_FORMULA_GEAR2) {
        a->a0 = (2.0 * h + h_last) / (h * (h + h_last));
        a->a1 = -(h + h_last) / (h * h_last);
        a->a2 = h / (h_last * (h + h_last));
    }
    pr->hk = s->hk;
    pr->ttarget = a->tn;
    st = CALL(s->pred_fn, pcx_pred, pr);
    if (st == PRED_DONE) seed = pr->out;
    else if (st == PRED_NONE) seed = s->xlast;
    else return LOOP_BAIL_PRED;
    for (r = 0; r < n; r++) if (!isfinite(seed[r])) return LOOP_BAIL_NEWTON;
    if (s->walk) {
        /* (the walk's sources: the table's row, the Python's own vector) */
        nt->u = s->wu + s->wj * n;
    } else {
        src->mode = SOURCES_U;
        src->t = a->tn;
        if (CALL(s->src_fn, pcx_sources, src) != SOURCES_DONE) return LOOP_BAIL_SOURCES;
        nt->u = src->out;
    }
    nt->formula = a->formula;
    nt->a0 = a->a0;
    nt->a1 = a->a1;
    nt->a2 = a->a2;
    nt->h = h;
    nt->theta = 0.0;
    nt->q1 = s->qring;
    nt->q2 = (a->formula == PCX_FORMULA_GEAR2) ? s->qring + n : s->qring;
    nt->iq1 = (a->formula == PCX_FORMULA_TRAPEZOIDAL || a->formula == PCX_FORMULA_THETA)
              ? s->iqring : s->qring;
    nt->xseed = seed;
    nt->have_cache = s->cur >= 0;
    a->nxt = (s->cur == 0) ? 1 : 0;
    if (s->cur >= 0) {
        nt->cx = s->xatt;
        nt->cC = s->cC[s->cur];
    }
    nt->Cj = s->cC[a->nxt];
    nt->Geqj = s->Geq[a->nxt];
    arm_early(s, a, nt);
    st = CALL(s->newton_fn, pcx_newton, nt);
    nt->early_lte = 0;
    a->early_checks = (st == NEWTON_EARLY || st == NEWTON_CONVERGED) ? nt->early_checks : 0;
    if (st == NEWTON_EARLY) return LOOP_JUDGE_EARLY;
    if (st == NEWTON_BAIL_MAXITER && s->newton_reject && h * 0.25 >= s->minstep)
        return LOOP_JUDGE_NEWTON;
    if (st != NEWTON_CONVERGED) return LOOP_BAIL_NEWTON;
    if (br) {
        br->C = nt->Cj;
        br->cmax = s->cmax;
        if (CALL(s->branch_fn, pcx_branch, br) != BRANCH_QUIET) return LOOP_BAIL_SCREEN;
    }
    return 0;
}

/* the error test on the attempt (as `_tran_lte_c.max_error` calls it): its
   scalars as `_tran_lte_c` makes them in Python floats -- Euler `-0.5 *
   (2.0 * h / (h + hl) if hl else 1.0)`, Gear `-h * (h + hl)` (the unary
   minus first), the trapezoid `-(h ** 2)`; Python's `**` on floats is libm
   `pow`, which can differ from `h * h` in the last bit, so `pow` it is
   (`-fno-builtin-pow` keeps the call) -- its vectors the attempt's and the
   rings'.  0 (`err` set), or LOOP_BAIL_LTE */
static long lte_setup(loop_t *s, const loop_attempt_t *a, const double *xc, const double *q,
                      const double *J)
{
    lte_t *le = (lte_t *) s->lte;
    long n = s->n;
    double h = a->h;
    double coef, hl = a->h_last, hl2 = 0.0, s12 = 0.0, s23 = 0.0, s123 = 0.0;
    int euler = a->formula == PCX_FORMULA_EULER;
    if (euler) {
        coef = -0.5 * ((hl != 0.0) ? (2.0 * h) / (h + hl) : 1.0);
        hl = 0.0;
    } else {
        if (!s->have_dt_last2) return LOOP_BAIL_LTE;
        hl2 = s->dt_last2;
        coef = (a->formula == PCX_FORMULA_GEAR2) ? (-h) * (h + hl) : -pow(h, 2.0);
        s12 = h + hl;
        s23 = hl + hl2;
        s123 = s12 + hl2;
    }
    if (!(isfinite(coef) && isfinite(s12) && isfinite(s23) && isfinite(s123)))
        return LOOP_BAIL_LTE;
    le->formula = euler ? LTE_FORMULA_EULER : LTE_FORMULA_DD3;
    le->h1 = h;
    le->h2 = hl;
    le->h3 = hl2;
    le->coef = coef;
    le->s12 = s12;
    le->s23 = s23;
    le->s123 = s123;
    le->xc = xc;
    le->xl = s->xlast;
    le->q = q;
    le->q1 = s->qring;
    le->q2 = s->qring + n;
    le->q3 = s->qring + 2 * n;
    le->iq1 = s->iqring;
    le->J = J;
    le->run = s->run;
    le->have_run = s->have_run;
    return 0;
}

/* the error test on the attempt's converged point (`lte_setup`, then
   pcx_lte).  0 (`err` set), or LOOP_BAIL_LTE */
static long error_test(loop_t *s, loop_attempt_t *a)
{
    newton_t *nt = (newton_t *) s->newton;
    lte_t *le = (lte_t *) s->lte;
    long st;
    if ((st = lte_setup(s, a, nt->xc, nt->qj, nt->Jj)) != 0) return st;
    if (CALL(s->lte_fn, pcx_lte, le) != LTE_DONE) return LOOP_BAIL_LTE;
    a->err = le->err;
    return 0;
}

/* a rejection the budget allows -- a growth always, a shrink while the
   point's rejections and the step allow (either family's budget at one
   point): 1, or 0 where it is spent */
static int retry_fits(const loop_t *s, const loop_attempt_t *a, double hn)
{
    double h = a->h;
    int growth = hn > h;
    return s->point_retries < s->max_retries
           && (growth || (s->rejects < s->max_reject && h > s->minstep));
}

/* the rejection retried at the point: its counts, the next step `max(hn,
   minstep)` -- either family's (`_SteppingLoop.judge`) */
static void retry_note(loop_t *s, const loop_attempt_t *a, double hn)
{
    int growth = hn > a->h;
    s->rejected++;
    s->point_retries++;
    s->rejects += !growth;
    s->imposed = 0;
    s->landing = a->landing;
    s->is_first = a->is_first;
    s->order_drop = 0;
    s->att_rejected = 1;
    s->t_break = a->t_break;
    s->h = pcx_py_max(hn, s->minstep);
}

static int retry_at_point(loop_t *s, loop_attempt_t *a)
{
    ctrl_t *ct = (ctrl_t *) s->ctrl;
    if (!retry_fits(s, a, ct->h_next)) return 0;
    attempt_made(s, a, 1);
    retry_note(s, a, ct->h_next);
    return 1;
}

/* THE EARLY ABORT's rejection (`_SteppingLoop.take`'s, `IntegralController.
   reject_early`): the law at the estimate the Newton stopped at, its retry
   record; the point's counts and the next step as `retry_at_point` keeps
   them; the Newton's iteration counted.  Nothing of the attempt is recorded
   as made: the cache, the rings and the `att_*` fields stay the last made
   attempt's, as the Python's raise leaves its state */
static void early_reject(loop_t *s, const loop_attempt_t *a)
{
    ctrl_t *ct = (ctrl_t *) s->ctrl;
    newton_t *nt = (newton_t *) s->newton;
    ct->err = nt->early_err;
    ct->p = (a->formula == PCX_FORMULA_EULER) ? 2.0 : 3.0;
    ct->h = a->h;
    ct->max_step = s->max_step;
    ct->point = s->npoint - 1;
    CALL(s->ctrl_fn, pcx_ctrl, ct);
    s->iters += nt->iters;
    s->early_checks += a->early_checks;
    if (s->hk >= 2) s->pred_fast++;
    s->rejected++;
    s->early_aborts++;
    s->point_retries++;
    s->rejects++;
    s->imposed = 0;
    s->landing = a->landing;
    s->is_first = a->is_first;
    s->order_drop = 0;
    s->t_break = a->t_break;
    s->h = pcx_py_max(ct->h_next, s->minstep);
}

/* THE NEWTON'S ITERATION LIMIT, rejected as `_SteppingLoop.take` rejects
   its exception: the attempt retried at a quarter of its step (above
   minstep: below it the rescue is the Python's), still imposed where its
   step was held; a rejection counted, its seed made; nothing of the attempt
   recorded as made, its iterations not counted (a failed solve adds none),
   the point's rejection counts and the law's record untouched */
static void newton_reject(loop_t *s, const loop_attempt_t *a)
{
    if (s->hk >= 2) s->pred_fast++;
    s->rejected++;
    s->newton_rejects++;
    s->imposed = a->clamped || s->imposed;
    s->landing = a->landing;
    s->is_first = a->is_first;
    s->order_drop = 0;
    s->t_break = a->t_break;
    s->h = a->h * 0.25;
}

/* 3. judge (`_SteppingLoop.judge`, `_LMMSteps.judge`): a fixed grid's next
   step; adaptive, the error test (`error_test`), the controller's law, the
   rejection budget (`retry_at_point`).  LOOP_JUDGE_ACCEPT (`h_next` set),
   LOOP_JUDGE_RETRY, or the status to stop with (the budget spent: the
   force-accept and its warning's record are the Python's, the law's
   record put back as it was) */
static long judge(loop_t *s, loop_attempt_t *a)
{
    ctrl_t *ct = (ctrl_t *) s->ctrl;
    long had_retry, kept_point, had_q, st;
    double kept_h, kept_err, kept_q;
    a->err = 0.0;
    if (!s->adaptive) {
        a->h_next = s->timestep;
        return LOOP_JUDGE_ACCEPT;
    }
    had_retry = ct->have_retry;
    kept_point = ct->retry_point;
    had_q = ct->have_q;
    kept_h = ct->retry_h;
    kept_err = ct->retry_err;
    kept_q = ct->q_carried;
    if ((st = error_test(s, a)) != 0) return st;
    ct->err = a->err;
    /* (`_max_error`'s exponent: the order + 1; the point the attempt left
       from as an index, where the Python keeps the object) */
    ct->p = (a->formula == PCX_FORMULA_EULER) ? 2.0 : 3.0;
    ct->h = a->h;
    ct->max_step = s->max_step;
    ct->point = s->npoint - 1;
    CALL(s->ctrl_fn, pcx_ctrl, ct);
    if (!ct->accepted) {
        if (retry_at_point(s, a)) return LOOP_JUDGE_RETRY;
        ct->have_retry = had_retry;
        ct->retry_point = kept_point;
        ct->retry_h = kept_h;
        ct->retry_err = kept_err;
        ct->have_q = had_q;
        ct->q_carried = kept_q;
        return LOOP_BAIL_FORCE;
    }
    s->rejects = s->point_retries = 0;
    a->h_next = ct->h_next;
    return LOOP_JUDGE_ACCEPT;
}

/* the accepted point's bookkeeping either family makes (`_SteppingLoop.
   accept`): its step into X and T, the statistics' extremes (set first,
   then replaced only when strictly beyond), the landing (`on_break`'s
   tolerance is `accept`'s), the point, the step history, the next step */
static void accept_point(loop_t *s, const loop_attempt_t *a, const double *x)
{
    long n = s->n;
    double h = a->h, tn = a->tn;
    int on_break = a->t + h >= a->t_break * (1.0 - 1e-12);
    s->att_rejected = 0;
    s->imposed = 0;
    s->landing = a->landing || on_break;
    s->is_first = 0;
    s->order_drop = 0;
    s->t = tn;
    COPY(s->X + s->k * n, x, n);
    s->T[s->k] = tn;
    if (!s->have_hm || h < s->hmin) s->hmin = h;
    if (!s->have_hm || h > s->hmax) s->hmax = h;
    s->have_hm = 1;
    s->breaks += on_break;
    COPY(s->xlast, x, n);
    s->have_dt_last2 = s->have_dt_last;
    s->dt_last2 = s->dt_last;
    s->dt_last = h;
    s->have_dt_last = 1;
    s->t_break = a->t_break;
    /* (`accept`'s `min(h_next, max_step)`: the first argument kept unless
       the second is strictly smaller) */
    s->h = s->adaptive ? pcx_py_min(a->h_next, s->max_step) : s->timestep;
    s->npoint++;
    s->k++;
}

/* 5. accept, the multistep family's: the attempt recorded, its order drop
   counted, the rings, the predictor's node, then the point's bookkeeping */
static void accept(loop_t *s, const loop_attempt_t *a)
{
    newton_t *nt = (newton_t *) s->newton;
    long n = s->n, r;
    double tn = a->tn;
    attempt_made(s, a, (int) s->adaptive);
    s->euler_steps += (a->formula == PCX_FORMULA_EULER && s->base != PCX_FORMULA_EULER);
    /* the rings (`_ring_push`): the newest first, the oldest dropped */
    for (r = s->hist_len - 1; r > 0; r--) {
        COPY(s->qring + r * n, s->qring + (r - 1) * n, n);
        COPY(s->iqring + r * n, s->iqring + (r - 1) * n, n);
    }
    COPY(s->qring, nt->qj, n);
    COPY(s->iqring, nt->iqj, n);
    /* the predictor's node (`_pred_note`: appended past every node, so its
       sort keeps it last and its dedupe -- `1e-14 * max(abs(t), 1.0)`,
       `_StagePredictor`'s -- compares it with the one before alone; the
       newest PRED_HIST_MAX kept, the capacity pred.c and the glue's `ht` /
       `hx` buffers share) */
    {
        double m = pcx_py_max(fabs(tn), 1.0);
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
    accept_point(s, a, nt->xc);
}
/* `remove_row_col` of the full `C` (n x n, row-major) into `dst` (m x m) */
static void reduce_mat(const double *C, double *dst, long n, long r)
{
    long i, j, o = 0;
    for (i = 0; i < n; i++) {
        if (i == r) continue;
        for (j = 0; j < n; j++) if (j != r) dst[o++] = C[i * n + j];
    }
}

/* `_LiveMonodromy._prune`'s test on the state one step on: 1 where a
   column is zero in P_n, P_{n-1} and Pq_n alike (`np.any(P != 0, axis=0)`:
   a NaN live) -- the prune then the Python's */
static int walk_dies(const loop_t *s)
{
    long m = s->m, w = s->w, r, c;
    const double *P0 = s->Pb[s->p0], *P1 = s->Pb[s->p1], *Pq = s->Pqb[s->pq];
    for (c = 0; c < w; c++) {
        int live = 0;
        for (r = 0; r < m && !live; r++)
            live = P0[r * w + c] != 0.0 || P1[r * w + c] != 0.0 || Pq[r * w + c] != 0.0;
        if (!live) return 1;
    }
    return 0;
}

/* the walk's step map after step `a` (`_walk_lmm`'s, `_lmmsens_c.step`): the
   step's J and C reduced, its coefficients (`companion_coefficients`:
   Gear-2's alphas, Euler's `1.0 / h` and `-1.0 / h`, the trapezoid's
   `2.0 / h`, `-2.0 / h` and -1), lmmsens.c on the live columns, then the
   rings moved on (`_LiveMonodromy.advance`; its prune, on its first steps
   and then every 64th, tested).  0, LOOP_WALK_MAP: the map handed back,
   the rings as they were, or LOOP_WALK_PRUNE: a column dead, the rings
   moved on and the prune the Python's */
static long walk_map(loop_t *s, const loop_attempt_t *a)
{
    newton_t *nt = (newton_t *) s->newton;
    lmmsens_t *mp = (lmmsens_t *) s->map;
    long fr = 3 - s->cn - s->cp, p2 = 3 - s->p0 - s->p1, st;
    double h = a->h;
    reduce_mat(nt->Jj, s->Jf, s->n, s->iref);
    reduce_mat(nt->Cj, s->Cr[fr], s->n, s->iref);
    mp->b = 0.0;
    if (a->formula == PCX_FORMULA_GEAR2) {
        mp->L = 3;
        mp->alpha[0] = a->a0;
        mp->alpha[1] = a->a1;
        mp->alpha[2] = a->a2;
    } else {
        mp->L = 2;
        mp->alpha[0] = (a->formula == PCX_FORMULA_EULER ? 1.0 : 2.0) / h;
        mp->alpha[1] = (a->formula == PCX_FORMULA_EULER ? -1.0 : -2.0) / h;
        if (a->formula != PCX_FORMULA_EULER) mp->b = -1.0;
    }
    mp->w = s->w;
    mp->P0 = s->Pb[s->p0];
    mp->P1 = s->Pb[s->p1];
    mp->Pq = s->Pqb[s->pq];
    mp->C1 = s->Cr[s->cn];
    mp->C2 = s->Cr[s->cp];
    mp->Cn = s->Cr[fr];
    mp->Jf = s->Jf;
    mp->Pn = s->Pb[p2];
    mp->Pqn = s->Pqb[1 - s->pq];
    /* (the capacitances' column lists the last map made, moved on one
       place: the first map of a call reads them afresh) */
    mp->reuse = s->maps ? 3 : 0;
    st = CALL(s->map_fn, pcx_lmmsens, mp);
    if (st != LMMSENS_DONE) {
        s->map_status = st;
        return LOOP_WALK_MAP;
    }
    s->maps++;
    s->p1 = s->p0;
    s->p0 = p2;
    s->pq = 1 - s->pq;
    s->cp = s->cn;
    s->cn = fr;
    s->mono_n++;
    if ((s->mono_n <= 4 || s->mono_n % 64 == 0) && walk_dies(s)) return LOOP_WALK_PRUNE;
    return 0;
}

/* the walk (`_lmmwalk_c`): the table's steps from `wj`, each made by `take`
   with its end the grid's, accepted (every one: no judge), then its map.
   It stops before a step a unit hands back, nothing of it made; or after a
   step whose map handed back or whose prune is the Python's
   (LOOP_WALK_MAP, LOOP_WALK_PRUNE: `wj` past it) */
static long walk(loop_t *s)
{
    loop_attempt_t a;
    long st;
    for (; s->wj < s->wsteps; s->wj++) {
        if (s->k >= s->cap) return LOOP_FULL;
        memset(&a, 0, sizeof a);
        a.t = s->t;
        a.h = s->wh[s->wj];
        a.tn = s->wt[s->wj];
        a.have_tn = 1;
        a.is_first = (int) s->is_first;
        a.t_break = INFINITY;
        if ((st = take(s, &a)) != 0) return st;
        accept(s, &a);
        if ((st = walk_map(s, &a)) != 0) {
            s->wj++;
            return st;
        }
    }
    return LOOP_DONE;
}

/* THE STAGE FAMILY (loop.h) */

/* the start `x`'s record in output set `o` (`A` assemblies): the stage
   whose bytes are `x`'s -- `_memo_get`'s key -- the latest put where
   several are (`_memo_put`'s update): its q, C, G, i into `a`.  1, or 0 */
static int stage_find(const loop_t *s, long o, long A, const double *x, loop_attempt_t *a)
{
    long n = s->n, nn = n * n, k;
    for (k = 3 * A - 1; k >= 0; k--) {
        if (memcmp(s->oset[o][0] + k * n, x, (size_t) n * sizeof(double))) continue;
        a->rq = s->oset[o][1] + k * n;
        a->ri = s->oset[o][2] + k * n;
        a->rC = s->oset[o][3] + k * nn;
        a->rG = s->oset[o][4] + k * nn;
        return 1;
    }
    return 0;
}

/* the core's passes at the start as `_tran_core.passes` makes them
   (formula 4; its want bits -- the pattern's where the core has one, into
   a C zeroed first -- C its own input and output, q in place, G and i from
   the core's bins) into `xr_*`, the record no generation holds: one call
   for the context's 'q' and the estimate's 'CGi', each output the same
   bits; the core's flags are not numpy's: cleared */
static void stage_pass(loop_t *s, loop_attempt_t *a)
{
    radau_t *rs = (radau_t *) s->radau;
    core_fn_t core = (core_fn_t) rs->core_fn;
    long n = s->n;
    if (s->pass_zero) memset(s->xr_C, 0, (size_t) (n * n) * sizeof(double));
    core(rs->core, s->xlast, rs->T, s->pass_bits, PCX_FORMULA_PASSES, 0.0, 0.0, 0.0, 1.0, 0.0,
         rs->dummy, rs->dummy, rs->dummy, rs->dummy, s->xr_C, s->xr_C, s->xr_q, rs->dummy,
         rs->dummy, rs->dummy, rs->dummy);
    COPY(s->xr_G, rs->binG, n * n);
    COPY(s->xr_i, rs->bini, n);
    feclearexcept(FE_ALL_EXCEPT);
    a->rq = s->xr_q;
    a->rC = s->xr_C;
    a->rG = s->xr_G;
    a->ri = s->xr_i;
}

/* the start's record as `_coupled_stage_context` reads its charge: in the
   last made attempt's generation (the first attempt: the entry's memo,
   which the glue copied into `xr_*` where it holds the start), else the
   core's passes -- counted, as `_tran_core.passes` counts them */
static void stage_start(loop_t *s, loop_attempt_t *a)
{
    if (s->made == 0 && s->first_rec) {
        a->rq = s->xr_q;
        a->rC = s->xr_C;
        a->rG = s->xr_G;
        a->ri = s->xr_i;
        return;
    }
    if (s->made > 0 && stage_find(s, s->gen0, s->A0, s->xlast, a)) return;
    stage_pass(s, a);
    a->passes++;
}

/* the start's record as the estimate reads it (`_memo_get` after the
   Newton's puts): this attempt's own set first (a stage the start's very
   bytes), else the context's -- a generation's, or the passes then made,
   made again (counted: the estimate's own; the same bits).  A record
   found is the estimate's memo hit (`radau.est:memo`) */
static void stage_est_rec(loop_t *s, loop_attempt_t *a)
{
    if (stage_find(s, a->set, a->A, s->xlast, a) || a->passes == 0) {
        a->estmemo = 1;
        return;
    }
    a->passes++;
}

/* u(tt) into row `row` of `ust`: a time this attempt asked before is that
   answer (the step's source memo, `_source_at`'s, by Python's `==`), else
   the table's -- a direct pass, counted.  0, or LOOP_BAIL_SOURCES */
static long stage_u(loop_t *s, loop_attempt_t *a, double tt, long row)
{
    sources_t *src = (sources_t *) s->src;
    long n = s->n, k;
    for (k = 0; k < a->nu; k++) {
        if (a->ut[k] != tt) continue;
        if (a->us[k] != row) COPY(s->ust + row * n, s->ust + a->us[k] * n, n);
        a->uhit++;
        return 0;
    }
    src->mode = SOURCES_U;
    src->t = tt;
    if (CALL(s->src_fn, pcx_sources, src) != SOURCES_DONE) return LOOP_BAIL_SOURCES;
    COPY(s->ust + row * n, src->out, n);
    a->ut[a->nu] = tt;
    a->us[a->nu] = row;
    a->nu++;
    a->umiss++;
    return 0;
}

/* `_radau_error_estimate` with the cost transform off (a fresh real
   factor) into the error test's `lte`: F1 = (dd1 Z1 + dd2 Z2 + dd3 Z3) / h
   (Z_i = Y_i - x_n), the right-hand side C_n F1 + f0 (numpy's matmul,
   dgemv row-major; f0 = -(i_n + u(t_n))), the real factor (gamma_r / h)
   C_n + G_n over the full width then reduced, numpy's solve (dgesv on its
   column-major copy; numpy's errstate there raises only on invalid), the
   estimate with +0.0 at the reference row -- numpy's operations in their
   order.  0, or LOOP_BAIL_ESTIMATE where numpy would warn or raise */
static long stage_estimate(loop_t *s, const loop_attempt_t *a, const double *Y)
{
    lte_t *le = (lte_t *) s->lte;
    long n = s->n, r = s->iref, i, j, ri, cj;
    int64_t m = n - 1, one = 1, info = 0;
    double h = a->h, g, *F = s->eF, *rhs = le->lte;
    volatile double sink;
    const double *x = s->xlast, *u0 = s->ust + 4 * n;
    int flags = s->under_ok ? (PCX_FP_FLAGS & ~FE_UNDERFLOW) : PCX_FP_FLAGS;
    feclearexcept(FE_ALL_EXCEPT);
    for (i = 0; i < n; i++)
        F[i] = ((s->dd1 * (Y[i] - x[i]) + s->dd2 * (Y[n + i] - x[i]))
                + s->dd3 * (Y[2 * n + i] - x[i])) / h;
    ((dgemv64_t) s->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, n, n, 1.0, a->rC, n, F, 1,
                           0.0, rhs, 1);
    for (i = 0; i < n; i++) rhs[i] = rhs[i] + (-(a->ri[i] + u0[i]));
    for (i = 0, ri = 0; i < n; i++) if (i != r) s->eb[ri++] = rhs[i];
    g = s->gamma_r / h;
    for (j = 0, cj = 0; j < n; j++) {
        for (i = 0, ri = 0; i < n; i++) {
            double v = g * a->rC[i * n + j] + a->rG[i * n + j];
            /* (the reference row and column: made, as numpy makes the full
               factor before it drops them -- their flags count) */
            if (i == r || j == r) {
                sink = v;
                continue;
            }
            s->eA[cj * m + ri] = v;
            ri++;
        }
        if (j != r) cj++;
    }
    (void) sink;
    if (fetestexcept(flags)) return LOOP_BAIL_ESTIMATE;
    feclearexcept(FE_ALL_EXCEPT);
    ((dgesv64_t) s->dgesv)(&m, &one, s->eA, &m, s->eipiv, s->eb, &m, &info);
    if (info != 0 || fetestexcept(FE_INVALID)) {
        feclearexcept(FE_ALL_EXCEPT);
        return LOOP_BAIL_ESTIMATE;
    }
    feclearexcept(FE_ALL_EXCEPT);
    pcx_put_ref(s->eb, rhs, n, r);
    return 0;
}

/* 2. take, the stage family's (`_rk_step_coupled`: the dense coupled
   Newton, the cost transform off): the stage times from the step's end as
   `_coupled_stage_context` makes them (`tn = t_end - h`, `tn + c_i h`), the
   step's nodes in merge order, the start's record, the sources at the stage
   times, the seeds, the C stage Newton writing the set the last made
   attempt's is not in, each stage's screen, the sources at the end
   (`_finish_stage_step`'s `iq`) and, adaptive, the start's record again,
   its sources and the estimate.  0, or the status to stop with: before the
   Newton nothing of the attempt is made; after it, its result is held for
   the Python's repeat (the glue's `_tran_radau_c.hold`) */
static long stage_take(loop_t *s, loop_attempt_t *a)
{
    radau_t *rs = (radau_t *) s->radau;
    long n = s->n, i, st, o;
    double h = a->h;
    const double *Y;
    a->tn = a->t + h;
    a->t0 = a->tn - h;
    for (i = 0; i < 3; i++) a->ts[i] = a->t0 + s->ct[i] * h;
    a->formula = 0;
    a->nu = a->umiss = a->uhit = a->passes = a->estmemo = 0;
    s->held = 0;
    s->cur_h = h;
    if (!pcx_stage_order(a->ts, a->tn, s->hk, s->ht, a->nodes, a->ord)) return LOOP_BAIL_NODES;
    stage_start(s, a);
    a->fast = pcx_stage_seeds(s->pred_fn, (pred_t *) s->pred, s->hk, a->ts, s->xlast, s->useed,
                              n);
    if (a->fast < 0) return LOOP_BAIL_PRED;
    for (i = 0; i < 3; i++) if ((st = stage_u(s, a, a->ts[i], i)) != 0) return st;
    o = a->set = (s->gen0 == 0) ? 1 : 0;
    COPY(s->hq, a->rq, n);
    rs->qn = s->hq;
    rs->u = s->ust;
    rs->seed = s->useed;
    rs->h = h;
    rs->Ya = s->oset[o][0];
    rs->qa = s->oset[o][1];
    rs->ia = s->oset[o][2];
    rs->Ca = s->oset[o][3];
    rs->Ga = s->oset[o][4];
    rs->klu.klu_num = s->klu_cur;
    if (CALL(s->radau_fn, pcx_radau, rs) != RADAU_CONVERGED) return LOOP_BAIL_RADAU;
    s->held = 1;
    if (rs->klu.klu_sym) s->klu_cur = rs->klu.klu_out;
    a->A = rs->assemblies;
    Y = a->Ys = s->oset[o][0] + (a->A - 1) * 3 * n;
    {
        long b = ((a->A - 1) * 3 + 2) * n;
        a->q3 = s->oset[o][1] + b;
        a->i3 = s->oset[o][2] + b;
        a->C3 = s->oset[o][3] + b * n;
        a->G3 = s->oset[o][4] + b * n;
    }
    a->cmax = s->cmax;
    if (!pcx_stage_screen(s->branch_fn, (branch_t *) s->branch,
                          s->oset[o][3] + (a->A - 1) * 3 * n * n, n, &a->cmax))
        return LOOP_BAIL_SCREEN;
    if ((st = stage_u(s, a, a->tn, 3)) != 0) return st;
    if (!s->adaptive) return 0;
    stage_est_rec(s, a);
    if ((st = stage_u(s, a, a->t0, 4)) != 0) return st;
    return stage_estimate(s, a, Y);
}

/* the coupled family's part of an attempt made: its set the last made
   attempt's (`gen0`), the one before kept by its stages (`ygen1`), its
   assemblies, the stage KLU's counts */
static void co_made(loop_t *s, const loop_attempt_t *a)
{
    radau_t *rs = (radau_t *) s->radau;
    long n = s->n;
    if (s->gen0 >= 0) {
        COPY(s->ygen1, s->oset[s->gen0][0], s->A0 * 3 * n);
        s->A1 = s->A0;
    }
    s->gen0 = a->set;
    s->A0 = a->A;
    s->st_assemblies += a->A;
    if (rs->klu.klu_sym) {
        s->st_klu_refactors += rs->klu.refactors;
        s->st_klu_factors += rs->klu.factors;
        s->st_klu_pivot += rs->klu.pivot_fallbacks;
        s->st_klu_resid += rs->klu.residual_fallbacks;
        s->st_klu_on++;
    }
}

/* ---- the transform's attempts (`tmode`: `_rk_step_transformed`) ---- */

/* (the family's retry law, below: the early abort's rejection takes it) */
static double stage_retry(loop_t *s, const loop_attempt_t *a, double worst, double p);

/* the memo's generation in slot g: its record k's state (n) */
static double *tg_state(const loop_t *s, long g, long k)
{
    return s->tg + (g * LOOP_TG_MAX + k) * s->n;
}

/* the record of the state x in slot g (`_memo_get`'s key: the state's
   bytes): its index, or -1 */
static long tg_find(const loop_t *s, long g, const double *x)
{
    long k;
    for (k = 0; k < s->tg_n[g]; k++)
        if (!memcmp(tg_state(s, g, k), x, (size_t) s->n * sizeof(double))) return k;
    return -1;
}

/* `_memo_get`'s record at x: its passes (LOOP_TG_*) -- this attempt's
   generation's, else the one before's; 0 where neither holds one */
static long tg_get(const loop_t *s, const double *x)
{
    long k = tg_find(s, s->tg_gen, x);
    if (k >= 0) return s->tg_f[s->tg_gen][k];
    k = tg_find(s, s->tg_prev, x);
    return k >= 0 ? s->tg_f[s->tg_prev][k] : 0;
}

/* `_memo_put` of the passes f at x: into this attempt's record, made from
   the one before's where that holds one (an attempt puts four at most
   into its new generation: LOOP_TG_MAX holds them) */
static void tg_put(loop_t *s, const double *x, long f)
{
    long g = s->tg_gen, k = tg_find(s, g, x), j;
    if (k < 0) {
        k = s->tg_n[g]++;
        memcpy(tg_state(s, g, k), x, (size_t) s->n * sizeof(double));
        j = tg_find(s, s->tg_prev, x);
        s->tg_f[g][k] = j >= 0 ? s->tg_f[s->tg_prev][j] : 0;
    }
    s->tg_f[g][k] |= f;
}

/* `_memo_step`: the last made attempt's generation becomes the one before
   and the attempt's own opens in the slot neither holds -- the two kept
   until it is made, so a stop takes the roll back (`tg_undo`) */
static void tg_roll(loop_t *s, loop_attempt_t *a)
{
    a->tg_gen0 = s->tg_gen;
    a->tg_prev0 = s->tg_prev;
    s->tg_prev = s->tg_gen;
    s->tg_gen = 3 - a->tg_gen0 - a->tg_prev0;
    s->tg_n[s->tg_gen] = 0;
}

/* an attempt not made, taken back: the memo's roll, and the point's f0
   where it made it (the Python's repeat makes and counts both) */
static void tc_undo(loop_t *s, const loop_attempt_t *a)
{
    s->tg_gen = a->tg_gen0;
    s->tg_prev = a->tg_prev0;
    s->have_f0 = a->had_f0;
    s->f0_point = a->f0_point0;
    if (a->repat) {
        /* (the new pattern's analysis and factor freed, the old put back:
           the Python's repeat makes its own) */
        radau_tf_t *tf = (radau_tf_t *) s->tc;
        radau_fold_t *fo = (radau_fold_t *) s->fold;
        void *sym = tf->Symbolic, *num = tf->Numeric;
        ((klu_free_numeric_t) s->klu_free_znum)(&num, tf->Common);
        ((klu_free_numeric_t) s->klu_free_sym)(&sym, tf->Common);
        tf->Symbolic = a->old_sym;
        tf->Numeric = a->old_num;
        tf->Ap = (int *) a->old_Ap;
        tf->Ai = (int *) a->old_Ai;
        fo->nz = (const uint8_t *) a->old_nz;
        fo->nnz = a->old_nnz;
        fo->nz_at = a->old_at;
        s->pcur = a->old_pcur;
    }
}

/* x_n's record's pass k (0 q, 1 i, 2 C, 3 G) as the C reads it: in the
   end set it was made in (`te_x`), else `xr_*` */
static const double *tc_xv(const loop_t *s, long k)
{
    if (s->te_x >= 0) return s->te[s->te_x][k];
    return k == 0 ? s->xr_q : k == 1 ? s->xr_i : k == 2 ? s->xr_C : s->xr_G;
}

/* one core call at x as `_tran_core.passes` makes it (formula 4, the want
   bits `bits` -- the pattern's where the core has one, C zeroed first
   then -- each pass's values whatever else is asked): q and C into the
   buffers given, G and i from the core's bins where given; the core's
   flags are not numpy's: cleared */
static void tc_pass(loop_t *s, const double *x, long bits, double *q, double *C, double *G,
                    double *iv)
{
    radau_tf_t *tf = (radau_tf_t *) s->tc;
    core_fn_t core = (core_fn_t) tf->core_fn;
    long n = s->n;
    if (s->pass_zero) memset(C, 0, (size_t) (n * n) * sizeof(double));
    core(tf->core, x, tf->T, bits, PCX_FORMULA_PASSES, 0.0, 0.0, 0.0, 1.0, 0.0,
         tf->dummy, tf->dummy, tf->dummy, tf->dummy, C, C, q, tf->dummy, tf->dummy, tf->dummy,
         tf->dummy);
    if (G) COPY(G, tf->binG, n * n);
    if (iv) COPY(iv, tf->bini, n);
    feclearexcept(FE_ALL_EXCEPT);
}

/* x_n's passes made into `xr_*` where the C holds none of `need` (one
   call makes them all: each pass's values whatever else it is asked) */
static void tc_xmake(loop_t *s, long need)
{
    if (s->te_x >= 0 || (s->xfields & need) == need) return;
    tc_pass(s, s->xlast, s->pass_bits, s->xr_q, s->xr_C, s->xr_G, s->xr_i);
    s->xfields = LOOP_TG_Q | LOOP_TG_I | LOOP_TG_C | LOOP_TG_G;
}

/* numpy's division of the complex (re, im) by the real h: its loop's form
   (Smith's, the divisor's imaginary part +0.0) */
static void tc_cdiv(double re, double im, double h, double *qr, double *qi)
{
    double rat = 0.0 / h, scl = 1.0 / (h + 0.0 * rat);
    *qr = (re + im * rat) * scl;
    *qi = (im - re * rat) * scl;
}

/* `_radau_frozen`'s scalars for the step h: P's rows 0 and 1 (the
   Newton's `p`) where h is not the step they were made for (numpy's
   `(diag(lam) @ Tinv) / h` under its errstate(all='raise'), the product
   numpy's in `tcm`), then the factors' s0 = lam0.real / h and c1 = lam1 /
   h.  0, or 1: a flag that errstate raises on, or h outside the fold's
   range (the Python's) */
static int tc_scalars(loop_t *s, double h)
{
    radau_tf_t *tf = (radau_tf_t *) s->tc;
    radau_fold_t *fo = (radau_fold_t *) s->fold;
    long k;
    if (!(1e-290 < h && h < 1e290)) return 1;
    feclearexcept(FE_ALL_EXCEPT);
    if (!s->tc_hp || h != s->tc_h) {
        for (k = 0; k < 6; k++)
            tc_cdiv(s->tcm[2 * k], s->tcm[2 * k + 1], h, &tf->p[2 * k], &tf->p[2 * k + 1]);
        if (fetestexcept(PCX_FP_FLAGS)) {
            feclearexcept(FE_ALL_EXCEPT);
            s->tc_hp = 0;
            return 1;
        }
        s->tc_h = h;
        s->tc_hp = 1;
    }
    fo->s0 = s->lam0r / h;
    tc_cdiv(s->lam1r, s->lam1i, h, &fo->c1r, &fo->c1i);
    feclearexcept(FE_ALL_EXCEPT);
    return 0;
}

/* the complex factor's new pattern (the fold's nonzeros off the last
   record's: `ComplexKLUSolver.prepare`'s new record, `_csc_of_dense`'s
   gather) into the loop's pattern buffer the current one is not in -- Ap,
   Ai and the mask laid out as A.T, its values into the fold's `ax` -- then
   `solve_prepared`'s new pattern: `klu_analyze` and a fresh
   `klu_z_factor` on it.  What it replaces is kept in the attempt for
   `tc_undo`.  0, or LOOP_BAIL_FOLD (KLU refused: the Python's) */
static long tc_repattern(loop_t *s, loop_attempt_t *a)
{
    radau_tf_t *tf = (radau_tf_t *) s->tc;
    radau_fold_t *fo = (radau_fold_t *) s->fold;
    long m = s->n - 1, i, j, k = 0, b = s->pcur == 0 ? 1 : 0;
    int *Ap = s->pAp[b], *Ai = s->pAi[b];
    uint8_t *nz = s->pnz[b];
    void *sym, *num;
    Ap[0] = 0;
    for (j = 0; j < m; j++) {
        for (i = 0; i < m; i++) {
            const double re = fo->cf[2 * (i * m + j)], im = fo->cf[2 * (i * m + j) + 1];
            const int on = re != 0.0 || im != 0.0;
            nz[j * m + i] = (uint8_t) on;
            if (on) {
                Ai[k] = (int) i;
                fo->ax[2 * k] = re;
                fo->ax[2 * k + 1] = im;
                k++;
            }
        }
        Ap[j + 1] = (int) k;
    }
    sym = ((klu_analyze_t) s->klu_analyze)((int) m, Ap, Ai, tf->Common);
    if (!sym) return LOOP_BAIL_FOLD;
    num = ((klu_factor_t) s->klu_zfactor)(Ap, Ai, fo->ax, sym, tf->Common);
    if (!num) {
        ((klu_free_numeric_t) s->klu_free_sym)(&sym, tf->Common);
        return LOOP_BAIL_FOLD;
    }
    a->repat = 1;
    a->old_sym = tf->Symbolic;
    a->old_num = tf->Numeric;
    a->old_Ap = tf->Ap;
    a->old_Ai = tf->Ai;
    a->old_nz = (void *) fo->nz;
    a->old_nnz = fo->nnz;
    a->old_at = fo->nz_at;
    a->old_pcur = s->pcur;
    tf->Symbolic = sym;
    tf->Numeric = num;
    tf->Ap = Ap;
    tf->Ai = Ai;
    fo->nz = nz;
    fo->nnz = k;
    fo->nz_at = 1;
    s->pcur = b;
    return 0;
}

/* `_radau_frozen`'s factors for the step (radau_fold.c, numpy's bytes) of
   x_n's C and G reduced: the real one into the loop's LU buffer
   column-major (`NumpyLU.make`; the Newton's first solve factors it), the
   complex one's nonzeros in the last record's pattern for KLU (the fold's
   `ax`, the Newton's `Ax`).  The fold counted where the Python takes it
   (P the step's already).  0, or LOOP_BAIL_FOLD: a flag, a non-finite
   input, nonzeros off the pattern (the Python's own path) */
static long tc_frozen(loop_t *s, loop_attempt_t *a)
{
    radau_fold_t *fo = (radau_fold_t *) s->fold;
    long n = s->n, m = n - 1, i, j;
    a->fold = s->tc_hp && a->h == s->tc_h;
    if (tc_scalars(s, a->h)) return LOOP_BAIL_FOLD;
    reduce_mat(tc_xv(s, 2), s->crr, n, s->iref);
    reduce_mat(tc_xv(s, 3), s->grr, n, s->iref);
    CALL(s->fold_fn, pcx_radau_fold, fo);
    if (fo->status != FOLD_STATUS_MADE) return LOOP_BAIL_FOLD;
    if (!fo->same) {
        long st = tc_repattern(s, a);
        if (st) return st;
    }
    for (j = 0; j < m; j++)
        for (i = 0; i < m; i++) s->lua[j * m + i] = fo->rf[i * m + j];
    return 0;
}

/* x_n's charge for the context (`_coupled_stage_context`, before the
   memo's roll: the record's, else a 'q' pass), the roll, then its C and G
   frozen (the record's, else a 'CG' pass, put) */
static void tc_start(loop_t *s, loop_attempt_t *a)
{
    const long cg = LOOP_TG_C | LOOP_TG_G;
    if (!(tg_get(s, s->xlast) & LOOP_TG_Q)) a->xpass++;
    tc_xmake(s, LOOP_TG_Q);
    tg_roll(s, a);
    if ((tg_get(s, s->xlast) & cg) != cg) {
        a->xpass++;
        tg_put(s, s->xlast, cg);
    }
    tc_xmake(s, cg);
}

/* THE EARLY ABORT's check at the probe's iterate (`_radau_early_check`):
   f0 = -(i_n + u(t_n)) once a point (`early_f0`: x_n's i -- the memo's,
   else an 'i' pass -- and the start's source), F1 at the iterate, `C_n F1 +
   f0` (numpy's dgemv), reduced, dgetrs on the probe's LU, normalised by
   lte.c's STAGE judgement (`normalised(kept=True)`: the reference row out,
   the running reference folded in as `_reference` folds it).  0 to go on,
   LOOP_JUDGE_EARLY (over the threshold: abandoned, `worst` the rejection's),
   or the status to stop with */
static long tc_check(loop_t *s, loop_attempt_t *a)
{
    radau_tf_t *tf = (radau_tf_t *) s->tc;
    lte_t *le = (lte_t *) s->lte;
    long n = s->n, r = s->iref, i, ri, st;
    int64_t m = n - 1, one = 1, info = 0;
    double h = a->h, *F = s->eF, *rhs = le->lte;
    const double *x = s->xlast, *Y = tf->Y, *Cn;
    int flags = s->under_ok ? (PCX_FP_FLAGS & ~FE_UNDERFLOW) : PCX_FP_FLAGS;
    if (!s->have_f0 || s->f0_point != s->npoint - 1) {
        if (!(tg_get(s, x) & LOOP_TG_I)) a->xpass++;
        tc_xmake(s, LOOP_TG_I);
        if ((st = stage_u(s, a, a->t0, 4)) != 0) return st;
        for (i = 0; i < n; i++) s->tcf0[i] = -(tc_xv(s, 1)[i] + s->ust[4 * n + i]);
        s->have_f0 = 1;
        s->f0_point = s->npoint - 1;
    }
    Cn = tc_xv(s, 2);
    feclearexcept(FE_ALL_EXCEPT);
    for (i = 0; i < n; i++)
        F[i] = ((s->dd1 * (Y[i] - x[i]) + s->dd2 * (Y[n + i] - x[i]))
                + s->dd3 * (Y[2 * n + i] - x[i])) / h;
    ((dgemv64_t) s->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, n, n, 1.0, Cn, n, F, 1,
                           0.0, rhs, 1);
    for (i = 0; i < n; i++) rhs[i] = rhs[i] + s->tcf0[i];
    if (fetestexcept(flags)) return LOOP_BAIL_ESTIMATE;
    for (i = 0, ri = 0; i < n; i++) if (i != r) s->eb[ri++] = rhs[i];
    ((dgetrs64_t) tf->dgetrs)("N", &m, &one, s->lua, &m, s->lupiv, s->eb, &m, &info, (size_t) 1);
    feclearexcept(FE_ALL_EXCEPT);
    if (info != 0) return LOOP_BAIL_ESTIMATE;
    pcx_put_ref(s->eb, rhs, n, r);
    le->xc = Y + 2 * n;
    le->xl = x;
    le->h1 = h;
    le->run = s->run;
    le->have_run = s->have_run;
    if (CALL(s->lte_fn, pcx_lte, le) != LTE_DONE) return LOOP_BAIL_LTE;
    if (s->keeps_run) {
        COPY(s->run, le->out, n);
        s->have_run = 1;
    }
    a->checked = 1;
    a->worst = le->worst;
    return le->err > s->tc_th ? LOOP_JUDGE_EARLY : 0;
}

/* radau_tc.c's Newton on the loop's struct for `maxit` iterations from the
   seeds in `useed`, its end passes into the end set `o`: its status, its
   iterations and refactors added to the attempt's */
static long tc_call(loop_t *s, loop_attempt_t *a, long maxit, long o)
{
    radau_tf_t *tf = (radau_tf_t *) s->tc;
    long st;
    tf->maxiter = maxit;
    tf->e_q = s->te[o][0];
    tf->e_i = s->te[o][1];
    tf->e_C = s->te[o][2];
    tf->e_G = s->te[o][3];
    st = CALL(s->tc_fn, pcx_radau_tc, tf);
    a->iters += tf->iters;
    a->refact += tf->refactored;
    return st;
}

/* on a new pattern (`tc_repattern`) the C Newton declines in the Python
   (`radau_tc:klu`) and the Python loop iterates: its first complex solve
   on the fresh factor, its second refactoring the same values (both made
   here: radau_tc.c with the numeric fresh, then with a refactor).  Armed,
   its first iteration alone, the check, then the C Newton serves the rest
   from the iterate -- the pattern now analysed -- its source reads the
   memo's hits; not armed, the Python loop makes every iteration.  The
   loop's iterations counted as its (`pyit`), the C's as the C's */
static long tc_newton_repat(loop_t *s, loop_attempt_t *a, long o)
{
    radau_tf_t *tf = (radau_tf_t *) s->tc;
    long n = s->n, st, it;
    tf->klu_fresh = 1;
    st = tc_call(s, a, 1, o);
    a->pyit = 1;
    a->iters = 0;
    if (st != RADAU_TC_MAXITER) return st;
    if (a->armed && s->tc_maxit > 1) {
        if ((st = tc_check(s, a)) != 0) return st;
        COPY(s->useed, tf->Y, 3 * n);
        tf->klu_fresh = 0;
        st = tc_call(s, a, s->tc_maxit - 1, o);
        a->tuhit = 3 + 3 * (a->iters - 1);
        a->served = st == RADAU_TC_CONVERGED;
        return st;
    }
    if (s->tc_maxit <= 1) return st;
    COPY(s->useed, tf->Y, 3 * n);
    tf->klu_fresh = 0;
    it = a->iters;
    st = tc_call(s, a, s->tc_maxit - 1, o);
    a->pyit += a->iters - it;
    a->iters = 0;
    return st;
}

/* the transform's Newton from the seeds: radau_tc.c on the loop's struct
   -- x_n's charge, the stage sources, the factors as `tc_frozen` left them
   (the LU made at the first real solve, the complex values refactored at
   the first complex one) -- with its end passes into the end set x_n's
   record is not in.  Armed, its first iteration alone (the probe), the
   check, then the rest from the iterate (`resume`: the factors kept); on a
   new pattern, `tc_newton_repat`.  0, LOOP_JUDGE_EARLY (abandoned at the
   check), or LOOP_BAIL_TRANSFORM (handed back, or not converged: the
   Python's dense fallback), LOOP_BAIL_END (no end passes: a non-finite
   stage), the check's stop */
static long tc_newton(loop_t *s, loop_attempt_t *a)
{
    radau_tf_t *tf = (radau_tf_t *) s->tc;
    long n = s->n, o = s->te_x == 0 ? 1 : 0, st, it;
    tf->qn = tc_xv(s, 0);
    tf->u = s->ust;
    tf->seed = s->useed;
    tf->h = a->h;
    tf->lu_a = s->lua;
    tf->lu_ipiv = s->lupiv;
    tf->lu_info = -1;
    tf->klu_fresh = 0;
    if (a->repat) {
        st = tc_newton_repat(s, a, o);
        if (st == LOOP_JUDGE_EARLY) return st;
    } else if (a->armed && s->tc_maxit > 1) {
        st = tc_call(s, a, 1, o);
        if (st == RADAU_TC_MAXITER) {
            a->probed = 1;
            if ((st = tc_check(s, a)) != 0) return st;
            COPY(s->useed, tf->Y, 3 * n);
            it = a->iters;
            st = tc_call(s, a, s->tc_maxit - 1, o);
            a->tuhit = 3 * (a->iters - it);
        }
        a->served = st == RADAU_TC_CONVERGED;
    } else {
        st = tc_call(s, a, s->tc_maxit, o);
        a->tuhit = 3 * (a->iters - 1);
        a->served = st == RADAU_TC_CONVERGED;
    }
    if (st != RADAU_TC_CONVERGED) return LOOP_BAIL_TRANSFORM;
    if (!tf->end_done) return LOOP_BAIL_END;
    a->pyend = !a->served;
    a->set = o;
    a->Ys = tf->Y;
    a->q3 = s->te[o][0];
    a->i3 = s->te[o][1];
    a->C3 = s->te[o][2];
    a->G3 = s->te[o][3];
    return 0;
}

/* `_stage_end_passes`: the end's record where the memo lacks passes of it
   -- the Newton's end passes handed and put whole (`radau.fuse:handed`,
   the pass they stand for counted where any was missing), or, the memo
   holding some, the missing ones by a pass (`:unhanded`) */
static void tc_end(loop_t *s, loop_attempt_t *a)
{
    const long all = LOOP_TG_Q | LOOP_TG_I | LOOP_TG_C | LOOP_TG_G;
    const double *y3 = a->Ys + 2 * s->n;
    long need = all & ~tg_get(s, y3);
    if (a->pyend) {
        /* (the Python loop converged: no hand-off -- the missing passes by
           one call, put) */
        a->fuse = need ? 4 : 5;
        if (need) tg_put(s, y3, need);
        return;
    }
    if (need == all) {
        a->fuse = 2;
        tg_put(s, y3, all);
    } else if (need) {
        a->fuse = 3;
        tg_put(s, y3, need);
    } else {
        a->fuse = 1;
    }
}

/* `_branch_after_coupled` on the stages: each one's C as `_C_at_state`
   reads it -- the memo's (`memo.C:memo`), else `cir.C`, the core's pass
   bit for bit, put -- and branch.c on it, the screen's maximum carried.
   Nothing read where nothing is screened.  0, or LOOP_BAIL_SCREEN */
static long tc_screen(loop_t *s, loop_attempt_t *a)
{
    branch_t *br = (branch_t *) s->branch;
    long n = s->n, nn = n * n, j;
    double cm = s->cmax;
    a->cmax = s->cmax;
    if (!br) return 0;
    a->screened = 1;
    for (j = 0; j < 3; j++) {
        const double *y = a->Ys + j * n, *Cj = a->C3;
        if (tg_get(s, y) & LOOP_TG_C) {
            a->mcmemo++;
        } else {
            a->cpass++;
            tg_put(s, y, LOOP_TG_C);
        }
        if (j < 2) {
            /* (C alone: `cir.C`'s, the pattern's where the core has one) */
            tc_pass(s, y, PCX_PASS_C | (s->pass_bits & PCX_WANT_PATTERN), s->tcq,
                    s->tcC + j * nn, NULL, NULL);
            Cj = s->tcC + j * nn;
        }
        br->C = Cj;
        br->cmax = cm;
        if (PCX_SCALL(s->branch_fn, pcx_branch, br) != BRANCH_QUIET) return LOOP_BAIL_SCREEN;
        cm = br->ref;
    }
    a->cmax = cm;
    return 0;
}

/* the estimate on the frozen real solve (`_radau_error_estimate` given the
   transform's `real`): x_n's record -- the memo's (`radau.est:memo`), else
   the core's 'CGi' -- F1, `C_n F1 + f0` (numpy's dgemv), reduced, dgetrs
   on the Newton's LU (`NumpyLU.solve`), +0.0 at the reference row, into
   the error test's `lte`.  0, or LOOP_BAIL_ESTIMATE where numpy would warn
   or the solve fails */
static long tc_estimate(loop_t *s, loop_attempt_t *a)
{
    lte_t *le = (lte_t *) s->lte;
    radau_tf_t *tf = (radau_tf_t *) s->tc;
    const long need = LOOP_TG_C | LOOP_TG_G | LOOP_TG_I;
    long n = s->n, r = s->iref, i, ri;
    int64_t m = n - 1, one = 1, info = 0;
    double h = a->h, *F = s->eF, *rhs = le->lte;
    const double *x = s->xlast, *u0 = s->ust + 4 * n, *Y = a->Ys, *Cn, *in;
    int flags = s->under_ok ? (PCX_FP_FLAGS & ~FE_UNDERFLOW) : PCX_FP_FLAGS;
    a->estimated = 1;
    if ((tg_get(s, x) & need) == need) a->estmemo = 1;
    else a->xpass++;
    tc_xmake(s, need);
    Cn = tc_xv(s, 2);
    in = tc_xv(s, 1);
    feclearexcept(FE_ALL_EXCEPT);
    for (i = 0; i < n; i++)
        F[i] = ((s->dd1 * (Y[i] - x[i]) + s->dd2 * (Y[n + i] - x[i]))
                + s->dd3 * (Y[2 * n + i] - x[i])) / h;
    ((dgemv64_t) s->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, n, n, 1.0, Cn, n, F, 1,
                           0.0, rhs, 1);
    for (i = 0; i < n; i++) rhs[i] = rhs[i] + (-(in[i] + u0[i]));
    if (fetestexcept(flags)) return LOOP_BAIL_ESTIMATE;
    for (i = 0, ri = 0; i < n; i++) if (i != r) s->eb[ri++] = rhs[i];
    ((dgetrs64_t) tf->dgetrs)("N", &m, &one, s->lua, &m, s->lupiv, s->eb, &m, &info, (size_t) 1);
    feclearexcept(FE_ALL_EXCEPT);
    if (info != 0) return LOOP_BAIL_ESTIMATE;
    pcx_put_ref(s->eb, rhs, n, r);
    return 0;
}

/* the transform's attempt past its start (`tc_take`): the factors, the
   seeds, the sources at the stage times, the Newton and its end passes,
   the end's record, the screens, the sources at the end and, adaptive, at
   the start, and the estimate -- `_rk_step_transformed` and
   `_finish_radau`'s order */
static long tc_body(loop_t *s, loop_attempt_t *a)
{
    long i, st;
    if ((st = tc_frozen(s, a)) != 0) return st;
    a->fast = pcx_stage_seeds(s->pred_fn, (pred_t *) s->pred, s->hk, a->ts, s->xlast, s->useed,
                              s->n);
    if (a->fast < 0) return LOOP_BAIL_PRED;
    for (i = 0; i < 3; i++) if ((st = stage_u(s, a, a->ts[i], i)) != 0) return st;
    if ((st = tc_newton(s, a)) != 0) return st;
    tc_end(s, a);
    if ((st = tc_screen(s, a)) != 0) return st;
    if ((st = stage_u(s, a, a->tn, 3)) != 0) return st;
    if (!s->adaptive) return 0;
    if ((st = stage_u(s, a, a->t0, 4)) != 0) return st;
    return tc_estimate(s, a);
}

/* 2. take, the transform's (`_rk_step_transformed`): the stage times from
   the step's end, the step's nodes in merge order, x_n's charge, the
   memo's roll and x_n's C and G frozen (`tc_start`), then `tc_body`.  0,
   or the status to stop with -- the memo's roll taken back: nothing of
   the attempt made (the Python's repeat refactors the same values, and
   the LU is the loop's own) */
static long tc_take(loop_t *s, loop_attempt_t *a)
{
    long i, st;
    double h = a->h;
    a->tn = a->t + h;
    a->t0 = a->tn - h;
    for (i = 0; i < 3; i++) a->ts[i] = a->t0 + s->ct[i] * h;
    a->formula = 0;
    a->nu = a->umiss = a->uhit = a->passes = a->estmemo = 0;
    a->iters = a->refact = a->xpass = a->cpass = a->fuse = a->mcmemo = a->fold = 0;
    a->tuhit = a->served = a->probed = a->checked = a->screened = a->estimated = 0;
    a->repat = a->pyit = a->pyend = 0;
    /* (armed as the stepping loop's `radau_arm` arms it: adaptive, a
       threshold, a full attempt left in the point's budget, above minstep) */
    a->armed = s->adaptive && s->tc_th > 0.0 && s->rejects + 1 < s->max_reject
               && s->point_retries + 1 < s->max_retries && h > s->minstep;
    s->held = 0;
    s->cur_h = h;
    if (!pcx_stage_order(a->ts, a->tn, s->hk, s->ht, a->nodes, a->ord)) return LOOP_BAIL_NODES;
    a->had_f0 = s->have_f0;
    a->f0_point0 = s->f0_point;
    tc_start(s, a);
    st = tc_body(s, a);
    /* (abandoned at the early abort, its roll and puts stand: the Python's
       raise leaves them) */
    if (st != 0 && st != LOOP_JUDGE_EARLY) tc_undo(s, a);
    return st;
}

/* an attempt's counts, made or abandoned at the early abort: its passes,
   sources, seeds and stage counts and the transform's own */
static void tc_count(loop_t *s, const loop_attempt_t *a)
{
    s->st_fast += a->fast;
    s->st_umiss += a->umiss;
    s->st_uhit += a->uhit;
    s->st_estmemo += a->estmemo;
    s->tc_frozen++;
    s->tc_iters += a->iters;
    s->tc_uhit += a->tuhit;
    s->tc_refactors += a->refact;
    s->tc_served += a->served;
    s->tc_probes += a->probed;
    s->tc_checks += a->checked;
    s->tc_xpasses += a->xpass;
    s->tc_cpasses += a->cpass;
    s->tc_ends += a->fuse != 0;
    s->tc_handed += a->fuse == 1 || a->fuse == 2;
    s->tc_endpass += a->fuse == 2 || a->fuse == 4;
    s->tc_unhanded += a->fuse == 3;
    s->tc_analyses += a->repat;
    s->tc_factors += a->repat;
    s->tc_klu += a->repat;
    s->tc_pyiters += a->pyit;
    /* (the Python loop's iterations after its first read the stage sources
       from the step's memo) */
    if (a->pyit > 1) s->tc_uhit += 3 * (a->pyit - 1);
    s->tc_pyends += a->pyend;
    s->tc_screened += a->screened;
    s->tc_mcmemo += a->mcmemo;
    s->tc_ests += a->estimated;
    s->tc_folds += a->fold;
    s->tc_foldoff += a->fold && a->repat;
}

/* a made attempt's new pattern kept: the analysis and factor it replaced
   freed, as `ComplexKLUSolver.solve_prepared` frees them at its new key
   (an attempt not made frees its own instead, `tc_undo`) */
static void tc_kept(loop_t *s, const loop_attempt_t *a)
{
    radau_tf_t *tf = (radau_tf_t *) s->tc;
    void *sym = a->old_sym, *num = a->old_num;
    if (!a->repat) return;
    if (num) {
        ((klu_free_numeric_t) s->klu_free_znum)(&num, tf->Common);
        s->tc_numfreed++;
    }
    if (sym) {
        ((klu_free_numeric_t) s->klu_free_sym)(&sym, tf->Common);
        s->tc_symfreed++;
    }
}

/* THE EARLY ABORT's rejection (`_SteppingLoop.take`'s `_RadauAbort`, the
   family's `reject_early`): the retry law at the check's worst unknown
   (`_retry_after`, its record), the point's counts and the next step as
   `take` keeps them; the attempt's counts made, its memo roll and puts
   standing; nothing of its outcome recorded as made */
static void tc_abort(loop_t *s, const loop_attempt_t *a)
{
    double hn = stage_retry(s, a, a->worst, s->pst);
    tc_kept(s, a);
    tc_count(s, a);
    s->tc_aborts++;
    s->rejected++;
    s->early_aborts++;
    s->point_retries++;
    s->rejects++;
    s->imposed = 0;
    s->landing = a->landing;
    s->is_first = a->is_first;
    s->order_drop = 0;
    s->t_break = a->t_break;
    s->h = pcx_py_max(hn, s->minstep);
}

/* the transform's part of an attempt made: its stages kept (`tcY`), the
   step `_radau_Pc` then holds, its pattern, its counts */
static void tc_made(loop_t *s, const loop_attempt_t *a)
{
    COPY(s->tcY, a->Ys, 3 * s->n);
    s->ph_made = a->h;
    s->have_ph = 1;
    tc_kept(s, a);
    tc_count(s, a);
}


/* an attempt judged (accepted, or retried at its point): the family's part
   (`co_made`, `tc_made`); the screen's maximum; the running maximum the
   judgement made and the estimate; the counts; the made attempt's `iq =
   -(i + u(t_end))`; the last attempt made's record */
static void stage_made(loop_t *s, const loop_attempt_t *a)
{
    lte_t *le = (lte_t *) s->lte;
    long n = s->n, i;
    const double *ue = s->ust + 3 * n;
    if (s->tmode) tc_made(s, a);
    else co_made(s, a);
    s->made++;
    if (s->branch) {
        s->cmax = a->cmax;
        s->screens += 3;
    }
    if (s->adaptive) {
        COPY(s->lest, le->lte, n);
        if (s->keeps_run) {
            COPY(s->run, le->out, n);
            s->have_run = 1;
        }
    }
    if (!s->tmode) {
        s->st_fast += a->fast;
        s->st_umiss += a->umiss;
        s->st_uhit += a->uhit;
        s->st_passes += a->passes;
        s->st_estmemo += a->estmemo;
    }
    for (i = 0; i < n; i++) s->liq[i] = -(a->i3[i] + ue[i]);
    s->held = 0;
    s->att_set = a->set;
    for (i = 0; i < 3; i++) s->att_ts[i] = a->ts[i];
    s->att_tend = a->tn;
    s->att_h = a->h;
    s->att_t = a->tn;
    s->attempts++;
}

/* `_StageSteps._retry_after`: the retry's order (`_retry_order`: 'assumed'
   p; from a point's second rejection -- the previous try left this point
   and was larger -- the order the two measured, `ln(e1/e2) / ln(h1/h2)`
   clamped to [RETRY_ORDER_MIN, p], kept as `_q_carried`; 'carried' takes
   it at a first rejection too), the record made, the next step
   `h max(MIN_SHRINK_RATIO, (SAFETY^p / worst)^(1/q))` -- Python's float
   arithmetic, `pow` C's (CPython's `**`) */
static double stage_retry(loop_t *s, const loop_attempt_t *a, double worst, double p)
{
    double q = p, h = a->h;
    long pt = s->npoint - 1;
    if (s->retry_order != LOOP_RETRY_ASSUMED) {
        if (s->have_sretry && s->sretry_point == pt && s->sretry_h > h * (1.0 + 1e-9)
            && s->sretry_worst > worst && worst > 0.0) {
            q = pcx_py_min(pcx_py_max(log(s->sretry_worst / worst) / log(s->sretry_h / h),
                                      s->rq_min), p);
            s->have_sq = 1;
            s->sq = q;
        } else if (s->retry_order == LOOP_RETRY_CARRIED && s->have_sq) {
            q = pcx_py_min(s->sq, p);
        }
    }
    s->have_sretry = 1;
    s->sretry_point = pt;
    s->sretry_h = h;
    s->sretry_worst = worst;
    return h * pcx_py_max(s->min_shrink, pow(pow(s->safety, p) / worst, 1.0 / q));
}

/* 3. judge, the stage family's (`_SteppingLoop.judge`, `_StageSteps.
   judge`): a fixed grid's next step; adaptive, lte.c's STAGE judgement of
   the estimate against the point, accepted at `err <= 1` (the step grown
   by `SAFETY err^(-1/p)` clamped to [K, 1/K], the retry record dropped),
   else `stage_retry` and the loop's budget.  LOOP_JUDGE_ACCEPT,
   LOOP_JUDGE_RETRY, or the status to stop with (the budget spent: the
   force-accept is the Python's, the family's retry record put back) */
static long stage_judge(loop_t *s, loop_attempt_t *a)
{
    lte_t *le = (lte_t *) s->lte;
    double h = a->h, p = s->pst, hn, kh = s->sretry_h, kw = s->sretry_worst, kq = s->sq;
    long khave = s->have_sretry, kpt = s->sretry_point, khq = s->have_sq;
    a->err = 0.0;
    if (!s->adaptive) {
        a->h_next = s->timestep;
        stage_made(s, a);
        return LOOP_JUDGE_ACCEPT;
    }
    le->xc = a->Ys + 2 * s->n;
    le->xl = s->xlast;
    le->h1 = h;
    le->run = s->run;
    le->have_run = s->have_run;
    if (CALL(s->lte_fn, pcx_lte, le) != LTE_DONE) return LOOP_BAIL_LTE;
    a->err = le->err;
    if (le->err <= 1.0) {
        double e = le->err > 1e-16 ? le->err : 1e-16;
        s->have_sretry = 0;
        a->h_next = h * pcx_py_min(1.0 / s->kclamp,
                                   pcx_py_max(s->kclamp, s->safety * pow(e, -1.0 / p)));
        stage_made(s, a);
        s->rejects = s->point_retries = 0;
        return LOOP_JUDGE_ACCEPT;
    }
    hn = stage_retry(s, a, le->worst, p);
    if (retry_fits(s, a, hn)) {
        stage_made(s, a);
        retry_note(s, a, hn);
        return LOOP_JUDGE_RETRY;
    }
    s->have_sretry = khave;
    s->sretry_point = kpt;
    s->sretry_h = kh;
    s->sretry_worst = kw;
    s->have_sq = khq;
    s->sq = kq;
    return LOOP_BAIL_FORCE;
}

/* 5. accept, the stage family's (`_roll_history` of a stage step): the
   rings -- the point's charge, its `iq` (`stage_made`'s) -- the stages and
   the point noted (`_pred_note`'s merge), then the point's bookkeeping */
static void stage_accept(loop_t *s, const loop_attempt_t *a)
{
    long n = s->n, k;
    const double *Y = a->Ys, *xv[4];
    const double *q3 = a->q3;
    for (k = s->hist_len - 1; k > 0; k--) {
        COPY(s->qring + k * n, s->qring + (k - 1) * n, n);
        COPY(s->iqring + k * n, s->iqring + (k - 1) * n, n);
    }
    COPY(s->qring, q3, n);
    COPY(s->iqring, s->liq, n);
    for (k = 0; k < 3; k++) xv[k] = Y + k * n;
    xv[3] = Y + 2 * n;
    pcx_stage_note(&s->hk, s->ht, s->hx, n, a->nodes, xv, a->ord);
    accept_point(s, a, Y + 2 * n);
    if (s->tmode) s->te_x = a->set;
}

/* the stage family's loop (`_SteppingLoop.execute` over `_StageSteps`):
   where, take, judge, accept -- the multistep family's `where` and point
   bookkeeping; its early abort and Newton rejection are not this family's */
static long stage_loop(loop_t *s)
{
    loop_attempt_t a;
    long st;
    for (;;) {
        if (!(s->t < s->tend)) return LOOP_DONE;
        if (s->k >= s->cap) return LOOP_FULL;
        if ((st = where(s, &a)) != 0) return st;
        if ((st = s->tmode ? tc_take(s, &a) : stage_take(s, &a)) != 0) {
            if (st != LOOP_JUDGE_EARLY) return st;
            tc_abort(s, &a);
            continue;
        }
        st = stage_judge(s, &a);
        if (st == LOOP_JUDGE_RETRY) continue;
        if (st != LOOP_JUDGE_ACCEPT) {
            /* (the transform's attempt not made: its memo roll taken back) */
            if (s->tmode) tc_undo(s, &a);
            return st;
        }
        stage_accept(s, &a);
    }
}

/* `_SteppingLoop.execute` for a multistep family whose every piece it was
   checked to stand in for, fixed or adaptive: attempt after attempt,
   where / take / judge / accept, the steps taken into X, T and the state.
   It stops BEFORE an attempt it cannot make as the Python would, that
   attempt's state untouched.  On any status the state is the last MADE
   attempt's: `cur` names the C and Geq buffer it wrote (the other may hold
   a discarded attempt's), the `att_*` fields its step, and the rings, the
   predictor's nodes, `xlast`, X and T accepted steps only.  The counts are
   this call's (zeroed here). */
long pcx_loop(loop_t *s)
{
    loop_attempt_t a;
    long st;
    s->k = 0;
    s->attempts = s->rejected = s->iters = s->euler_steps = s->breaks = s->screens = 0;
    s->pred_fast = s->ltes = s->lte_euler = s->early_aborts = s->early_checks = 0;
    s->newton_rejects = 0;
    s->maps = 0;
    if (s->walk) return s->status = walk(s);
    s->st_fast = s->st_umiss = s->st_uhit = s->st_passes = s->st_estmemo = s->st_assemblies = 0;
    s->st_klu_refactors = s->st_klu_factors = s->st_klu_pivot = s->st_klu_resid = 0;
    s->st_klu_on = 0;
    s->tc_frozen = s->tc_iters = s->tc_uhit = s->tc_refactors = s->tc_served = 0;
    s->tc_probes = s->tc_checks = s->tc_aborts = s->tc_xpasses = s->tc_cpasses = 0;
    s->tc_ends = s->tc_handed = s->tc_endpass = s->tc_unhanded = s->tc_screened = 0;
    s->tc_mcmemo = s->tc_ests = s->tc_folds = 0;
    s->tc_analyses = s->tc_factors = s->tc_foldoff = s->tc_klu = s->tc_pyiters = 0;
    s->tc_symfreed = s->tc_numfreed = 0;
    s->tc_pyends = 0;
    if (s->family == LOOP_FAMILY_STAGE) return s->status = stage_loop(s);
    for (;;) {
        /* `while self.t < self.tend` */
        if (!(s->t < s->tend)) return s->status = LOOP_DONE;
        if (s->k >= s->cap) return s->status = LOOP_FULL;
        if ((st = where(s, &a)) != 0) return s->status = st;
        st = take(s, &a);
        if (st == LOOP_JUDGE_EARLY) {
            early_reject(s, &a);
            continue;
        }
        if (st == LOOP_JUDGE_NEWTON) {
            newton_reject(s, &a);
            continue;
        }
        if (st != 0) return s->status = st;
        st = judge(s, &a);
        if (st == LOOP_JUDGE_RETRY) continue;
        if (st != LOOP_JUDGE_ACCEPT) return s->status = st;
        /* 4. events: none -- a circuit with a declared state event declines
           at the loop's entry, so no step is cut here */
        accept(s, &a);
    }
}
