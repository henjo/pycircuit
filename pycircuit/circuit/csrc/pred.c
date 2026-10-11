/* pred.c -- the stage predictor's fast path in C, bit for bit
   `_StagePredictor._predict_fast` or handed back; its glue is
   `_tran_pred_c`, and the C loop (loop.c) calls it for every attempt's
   seed.  It tests floating-point flags: THE FLAG RULE (pcx_fenv.h) holds
   (MUST KEEP: each operation a flag test covers stores its result through
   a pointer, `motion`, `V`, `out`, before the test). */
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#include "pcx_abi.h"
/* (its min / max / clip are written out, where the other units call
   pcx_minmax.h: the named helpers moved this unit's instruction
   schedule, and a unit that tests flags keeps the code it was proven
   with -- each says the rule it follows, as pcx_minmax.h names them) */
#include "pcx_fenv.h"
#include "pred.h"

/* the flags numpy checks after its operations (inexact is not one) */

/* every return: the flags cleared, `status` set and returned */
static long done(pred_t *s, long status)
{
    feclearexcept(FE_ALL_EXCEPT);
    return s->status = status;
}

/* `_pred_fast_nodes` and the tests before the fit: the history's indices
   newest first into `nodes` (the 1e-13 duplicates dropped -- the first
   eight recorded, the history's capacity, pred.h; only the first take <= 5
   are read), `*kp` counting them all.  0, or the status to return
   (PRED_NONE, PRED_GENERAL) */
static long fast_nodes(const pred_t *s, long *nodes, long *kp)
{
    long hk = s->hk, k = 0, e;
    double prev_t = 0.0, prev_kept = 0.0, gap = INFINITY;
    int have_prev = 0;
    /* (the Python's tests in its order: too short, then not wholly behind
       the target) */
    if (hk < 2) return PRED_NONE;
    if (!(s->ht[hk - 1] < s->ttarget)) return PRED_GENERAL;
    /* `_pred_fast_nodes`: newest first; Python's `min(gap, prev - t0)` and
       `max(abs(t0), 1.0)` as ternaries -- the first argument kept unless
       the second is strictly beyond it (fmin/fmax differ on a NaN) */
    for (e = hk - 1; e >= 0; e--) {
        double t0 = s->ht[e];
        if (have_prev) {
            double a = fabs(t0), m;
            if (!(t0 < prev_t)) return PRED_GENERAL;
            if (prev_t - t0 < gap) gap = prev_t - t0;
            m = (1.0 > a) ? 1.0 : a;
            if (fabs(t0 - prev_kept) <= 1e-13 * m) { prev_t = t0; continue; }
        }
        if (k < PRED_HIST_MAX) nodes[k] = e;
        k++;
        prev_t = t0;
        prev_kept = t0;
        have_prev = 1;
    }
    /* the tie test, as `_predict_fast` writes it */
    if (!(gap > 4.5e-16 * (s->ttarget - s->ht[0]))) return PRED_GENERAL;
    if (k < 2) return PRED_NONE;
    *kp = k;
    return 0;
}

/* `_fit_fast`'s weights `w` for the `take` nodes (the memo: `tau`, `wk`).
   0, or the status to return (PRED_NONE, PRED_BAIL_FLAGS) */
static long fit(pred_t *s, const long *nodes, long take)
{
    long j, c;
    double scale, tat, d[PRED_FIT_MAX];
    /* _fit_fast: the times about the target, scaled by the largest --
       Python floats there, whose flags numpy never sees (one raised here
       reaches the Vandermonde's test: a hand-back, never a difference);
       `max(abs(v) for v in d)` keeps the first of equal values, as the
       strict `>` does */
    tat = s->ttarget;
    for (j = 0; j < take; j++) {
        d[j] = s->ht[nodes[j]] - tat;
        if (isnan(d[j])) return PRED_NONE;
    }
    scale = fabs(d[0]);
    for (j = 1; j < take; j++) if (fabs(d[j]) > scale) scale = fabs(d[j]);
    if (!(scale < INFINITY) || scale <= 0.0) return PRED_NONE;
    {
        double tau[PRED_FIT_MAX];
        for (j = 0; j < take; j++) tau[j] = d[j] / scale;
        /* the memo's key is tau's bytes (`struct.pack`): memcmp, not ==,
           so -0.0 and +0.0 are two keys there as here */
        if (!(s->wk == take && memcmp(tau, s->tau, (size_t) take * sizeof(double)) == 0)) {
            /* np.vander(tau, take, increasing=True): 1, then running
               products; its transpose column-major is it row-major */
            int64_t m64 = take, one64 = 1, info64 = 0;
            for (j = 0; j < take; j++) {
                double p = tau[j];
                s->V[j * take] = 1.0;
                if (take > 1) s->V[j * take + 1] = p;
                for (c = 2; c < take; c++) { p = p * tau[j]; s->V[j * take + c] = p; }
            }
            /* (numpy's vander is a ufunc: what it raised, it reports) */
            if (fetestexcept(PCX_FP_FLAGS)) return PRED_BAIL_FLAGS;
            for (j = 0; j < take; j++) s->w[j] = (j == 0) ? 1.0 : 0.0;
            s->wk = 0;
            ((dgesv64_t) s->dgesv)(&m64, &one64, s->V, &m64, s->ipiv, s->w, &m64, &info64);
            /* (numpy's solve: an info is its LinAlgError, which `_fit_fast`
               answers with None; an invalid operation is answered None
               here too -- numpy clears what LAPACK raises (lte.c), but from
               V's finite entries, all within [-1, 1], none arises) */
            if (info64 != 0 || fetestexcept(FE_INVALID)) return PRED_NONE;
            feclearexcept(FE_ALL_EXCEPT);
            for (j = 0; j < take; j++) if (!isfinite(s->w[j])) return PRED_NONE;
            memcpy(s->tau, tau, (size_t) take * sizeof(double));
            s->wk = take;
        }
    }
    return 0;
}

/* `w @ X` into `out`, clamped to the newest node `x0` widened by the
   motion and the reach.  0, or PRED_BAIL_FLAGS */
static long clamped_product(pred_t *s, const long *nodes, long take, const double *x0)
{
    long n = s->n, j, i;
    double ratio, wmul, dt_last;
    /* w @ X, X the nodes taken, newest first: copied into `xt` as the
       Python's fresh (take, n) array lays them out (`hx` is oldest first,
       with duplicates), then numpy's own dgemv, the call its matmul makes */
    for (j = 0; j < take; j++) memcpy(s->xt + j * n, s->hx + nodes[j] * n, (size_t) n * sizeof(double));
    ((dgemv64_t) s->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_TRANS, take, n, 1.0, s->xt, n, s->w, 1,
                           0.0, s->out, 1);
    /* the clamp: the newest node, widened by the motion and the reach --
       `ratio` and `max(ratio, 1.0)` Python floats (the first argument kept
       unless the second is strictly greater), `PRED_CLAMP * motion *
       wmul` numpy's left to right, then numpy's `clip(p, lo, hi)`, which
       is min(max(p, lo), hi) as the two ternaries write it -- for a p
       that is not NaN (numpy's clip keeps a NaN p, these do not) */
    dt_last = fabs(s->ht[nodes[0]] - s->ht[nodes[1]]);
    ratio = (dt_last > 0.0) ? fabs(s->ttarget - s->ht[nodes[0]]) / dt_last : 1.0;
    wmul = (1.0 > ratio) ? 1.0 : ratio;
    for (i = 0; i < n; i++) {
        double wc = (s->clamp * s->motion[i]) * wmul;
        double lo = x0[i] - wc, hi = x0[i] + wc, p = s->out[i];
        double m = (p > lo) ? p : lo;
        s->out[i] = (m < hi) ? m : hi;
    }
    if (fetestexcept(PCX_FP_FLAGS)) return PRED_BAIL_FLAGS;
    return 0;
}

/* `_StagePredictor._predict_fast` with `_pred_fast_nodes` and `_fit_fast`:
   the nodes newest first (the 1e-13 dedupe, the smallest gap, the tie
   test), the fit in the target's coordinate (the times in Python's float
   order, `np.vander`'s running products, numpy's own `dgesv` for the
   weights and `dgemv` for `w @ X`), the clamp to the newest node widened by
   the motion between the two newest (numpy's `clip`).
   Assumes: the struct as pred.h says it -- `deg` <= 4 (the fit's
   capacity: d, tau, w, ipiv hold 5, V 25, xt 5 x n), `ht` and `hx` hold
   hk nodes, oldest first.
   PRED_DONE: `out` the clamped prediction.  PRED_NONE: `_predict_fast`'s
   None (fewer than two nodes, the degree below 1, a NaN or degenerate
   time spread, a singular or non-finite fit).  PRED_GENERAL: a history
   the fast path does not take (the newest node not strictly before the
   target, times not strictly ascending, a tie).  PRED_BAIL_FLAGS: a flag
   numpy would report (the motion, the Vandermonde, the product and the
   clamp).  On every status but PRED_DONE `out` is not the answer (it may
   hold a product or part of the clamp); `xt`, `motion`, `V`, `ipiv` are
   work; `tau`, `w`, `wk` change only by a solve that succeeded (`wk` is 0
   while it runs).  Every return leaves the flags clear (`done`). */
long pcx_pred(pred_t *s)
{
    long n = s->n, k = 0, i, take, st, nodes[PRED_HIST_MAX];
    const double *x0, *x1;
    if ((st = fast_nodes(s, nodes, &k)) != 0) return done(s, st);
    /* min(int(deg) + 1, len(nodes)) */
    take = s->deg + 1 < k ? s->deg + 1 : k;
    if (take < 2) return done(s, PRED_NONE);
    feclearexcept(FE_ALL_EXCEPT);
    x0 = s->hx + nodes[0] * n;
    x1 = s->hx + nodes[1] * n;
    for (i = 0; i < n; i++) s->motion[i] = fabs(x0[i] - x1[i]);
    if (fetestexcept(PCX_FP_FLAGS)) return done(s, PRED_BAIL_FLAGS);
    if ((st = fit(s, nodes, take)) != 0) return done(s, st);
    if ((st = clamped_product(s, nodes, take, x0)) != 0) return done(s, st);
    return done(s, PRED_DONE);
}
