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
    Y = s->oset[o][0] + (a->A - 1) * 3 * n;
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

/* an attempt judged (accepted, or retried at its point): its set the last
   made attempt's (`gen0`), the one before kept by its stages (`ygen1`);
   the screen's maximum; the running maximum the judgement made; the stage
   KLU's counts and the family's; the last attempt made's record */
static void stage_made(loop_t *s, const loop_attempt_t *a)
{
    radau_t *rs = (radau_t *) s->radau;
    lte_t *le = (lte_t *) s->lte;
    long n = s->n, i;
    if (s->gen0 >= 0) {
        COPY(s->ygen1, s->oset[s->gen0][0], s->A0 * 3 * n);
        s->A1 = s->A0;
    }
    s->gen0 = a->set;
    s->A0 = a->A;
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
    s->st_fast += a->fast;
    s->st_umiss += a->umiss;
    s->st_uhit += a->uhit;
    s->st_passes += a->passes;
    s->st_estmemo += a->estmemo;
    s->st_assemblies += a->A;
    if (rs->klu.klu_sym) {
        s->st_klu_refactors += rs->klu.refactors;
        s->st_klu_factors += rs->klu.factors;
        s->st_klu_pivot += rs->klu.pivot_fallbacks;
        s->st_klu_resid += rs->klu.residual_fallbacks;
        s->st_klu_on++;
    }
    {
        long b = ((a->A - 1) * 3 + 2) * n;
        const double *i3 = s->oset[a->set][2] + b, *ue = s->ust + 3 * n;
        for (i = 0; i < n; i++) s->liq[i] = -(i3[i] + ue[i]);
    }
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
    le->xc = s->oset[a->set][0] + ((a->A - 1) * 3 + 2) * s->n;
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
    long n = s->n, k, b = ((a->A - 1) * 3 + 2) * n;
    const double *Y = s->oset[a->set][0] + (a->A - 1) * 3 * n, *xv[4];
    const double *q3 = s->oset[a->set][1] + b;
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
        if ((st = stage_take(s, &a)) != 0) return st;
        st = stage_judge(s, &a);
        if (st == LOOP_JUDGE_RETRY) continue;
        if (st != LOOP_JUDGE_ACCEPT) return st;
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
