/* stagewalk.c -- the shooting's plain stage walk (contract: stagewalk.h): the
   steps of `_PeriodWalks._walk_stage` for a fully implicit tableau on the
   walk's KLU, each as the inner transient's step and the walk's step map
   make it -- the seeds (`_pred_stages`), the C stage Newton, the screen of
   each stage (`_branch_after_coupled`), the step's finish and its roll
   (`_finish_stage_step`, `_roll_history`, `_pred_note`), the walk's step map
   on its numeric (`_sensklu_c.step`) and the live columns' drop
   (`_LiveColumns.advance`) -- bit for bit, or stopped before the step it
   cannot make as the Python would.  Its glue is `shooting/_stagewalk_c`. */
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "pcx_abi.h"
#include "pcx_minmax.h"
/* the units this one calls through the pointers it is handed: their structs
   and entries (and lmmsens.h's dgemm64_t, numpy's own dgemm) */
#include "pcx_klu.h"
#include "radau.h"
#include "pred.h"
#include "branch.h"
#include "sensklu.h"
#include "lmmsens.h"
#include "stagewalk.h"

/* a call to a unit through its pointer, typed by the unit's own declaration */
#define CALL(fn, entry, p) (((__typeof__(&entry)) (fn))(p))
#define COPY(dst, src, cnt) memcpy((dst), (src), (size_t) (cnt) * sizeof(double))

/* `C_n @ P` into `base` (m x k, row-major), as numpy's matmul makes it of the
   C-ordered `C_n` and the live columns as `_LiveColumns` holds them: one
   column numpy's gemv (row-major, not transposed), more numpy's dgemm -- its
   bits the same for the column-major view and a C-ordered copy (checked on
   140 shapes): the operand as it stands */
static void product(stagewalk_t *t)
{
    long m = t->m, k = t->k;
    if (k == 1) {
        ((dgemv64_t) t->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, m, m, 1.0, t->Cn, m,
                               t->P, 1, 0.0, t->base, 1);
        return;
    }
    ((dgemm64_t) t->dgemm)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS,
                           t->trans ? PCX_CBLAS_TRANS : PCX_CBLAS_NO_TRANS, m, k, m, 1.0, t->Cn,
                           m, t->P, t->ld, 0.0, t->base, k);
}

/* `np.concatenate((x[:iref], x[iref + 1:]))` of the full `src` into `dst` */
static void reduce_vec(const double *src, double *dst, long n, long r)
{
    long i, o = 0;
    for (i = 0; i < n; i++) if (i != r) dst[o++] = src[i];
}

/* the reduced capacitance (`_reduce_ndarray`: the reference row and column
   dropped) of the full `C` (n x n, row-major) into `dst` (m x m) */
static void reduce_mat(const double *C, double *dst, long n, long r)
{
    long i, j, o = 0;
    for (i = 0; i < n; i++) {
        if (i == r) continue;
        for (j = 0; j < n; j++) if (j != r) dst[o++] = C[i * n + j];
    }
}

/* `_pred_note`'s merge of the step's nodes -- the three stages at their
   times and the point at the step's end, sorted by time as Python's stable
   sort sorts them (`ord`), each after the history's newest: appended, one
   within `1e-14 * max(abs(t), 1.0)` of the newest replacing it, the newest
   PRED_HIST_MAX kept (an append past them drops the oldest) */
static void note(stagewalk_t *t, const double *tn, const double *const *xv, const long *ord)
{
    long n = t->n, e;
    for (e = 0; e < 4; e++) {
        double te = tn[ord[e]];
        double mx = pcx_py_max(fabs(te), 1.0);
        if (t->hk && fabs(te - t->ht[t->hk - 1]) <= 1e-14 * mx) {
            t->hk--;
        } else if (t->hk == PRED_HIST_MAX) {
            memmove(t->ht, t->ht + 1, (PRED_HIST_MAX - 1) * sizeof(double));
            memmove(t->hx, t->hx + n, (size_t) ((PRED_HIST_MAX - 1) * n) * sizeof(double));
            t->hk = PRED_HIST_MAX - 1;
        }
        t->ht[t->hk] = te;
        COPY(t->hx + t->hk * n, xv[ord[e]], n);
        t->hk++;
    }
}

/* the step's nodes -- its stage times and its end -- into `tn`, and the
   order Python's stable sort takes them in into `ord`: 0, or the bail where
   the first is not after the history's newest */
static long order_nodes(stagewalk_t *t, const double *ts, double tj, double *tn, long *ord)
{
    long i, e;
    for (i = 0; i < 3; i++) tn[i] = ts[i];
    tn[3] = tj;
    for (e = 0; e < 4; e++) ord[e] = e;
    for (e = 1; e < 4; e++) {
        long v = ord[e], f = e;
        while (f > 0 && tn[ord[f - 1]] > tn[v]) {
            ord[f] = ord[f - 1];
            f--;
        }
        ord[f] = v;
    }
    if (t->hk && !(tn[ord[0]] >= t->ht[t->hk - 1])) return STAGEWALK_BAIL_NODES;
    return 0;
}

/* the seeds, `_pred_or(np.array(x), t)` at each stage time: 0, or the bail
   where the predictor's general path makes one */
static long seeds(stagewalk_t *t, const double *ts)
{
    pred_t *pr = (pred_t *) t->pred;
    long i, n = t->n, st;
    for (i = 0; i < 3; i++) {
        pr->hk = t->hk;
        pr->ttarget = ts[i];
        st = CALL(t->pred_fn, pcx_pred, pr);
        if (st == PRED_DONE) COPY(t->seed + i * n, pr->out, n);
        else if (st == PRED_NONE) COPY(t->seed + i * n, t->x, n);
        else return STAGEWALK_BAIL_PRED;
    }
    return 0;
}

/* the screen of the last assembly's stages (`A` of them made), the running
   maximum carried into `*cmax`: 0, or the bail at a stage not QUIET */
static long screen(stagewalk_t *t, const radau_t *rs, long A, double *cmax)
{
    branch_t *br = (branch_t *) t->branch;
    long i, nn = t->n * t->n;
    if (!br) return 0;
    for (i = 0; i < 3; i++) {
        br->C = rs->Ca + (A - 1) * 3 * nn + i * nn;
        br->cmax = *cmax;
        if (CALL(t->branch_fn, pcx_branch, br) != BRANCH_QUIET) return STAGEWALK_BAIL_SCREEN;
        *cmax = br->ref;
    }
    return 0;
}

/* the walk's step map into `Xw`: `C_n P`, the block refactored and solved
   (sensklu.c's result) */
static long step_map(stagewalk_t *t, const radau_t *rs, double *Xw)
{
    sensklu_t *sk = (sensklu_t *) t->sens;
    product(t);
    sk->Ax = rs->klu.kAx;
    sk->num = t->walk_num;
    sk->base = t->base;
    sk->m = t->m;
    sk->k = t->k;
    sk->s = 3;
    sk->X = Xw;
    return CALL(t->sens_fn, pcx_sensklu, sk);
}

/* the step made: the step map's numeric and counts as `block_step` commits
   them, the Newton's, the finish and the roll (the rings: iq = -(i + u(t)),
   q at the point; the predictor's nodes; the point), the records of step
   `j`: its stages `Y` of the last of `A` assemblies */
static void commit(stagewalk_t *t, const radau_t *rs, long A, double cmax, const double *Uj,
                   const double *tn, const long *ord, long j)
{
    const sensklu_t *sk = (const sensklu_t *) t->sens;
    long n = t->n, nn = n * n, i, e;
    const double *Y = rs->Ya + (A - 1) * 3 * n, *q3 = rs->qa + (A - 1) * 3 * n + 2 * n;
    const double *i3 = rs->ia + (A - 1) * 3 * n + 2 * n, *xv[4];
    t->walk_num = sk->num;
    t->walk_refactors += sk->refactors;
    t->walk_fallbacks += sk->fallbacks;
    t->walk_factors += sk->factors;
    t->klu_refactors += rs->klu.refactors;
    t->klu_factors += rs->klu.factors;
    t->klu_pivot_fallbacks += rs->klu.pivot_fallbacks;
    t->klu_residual_fallbacks += rs->klu.residual_fallbacks;
    t->assemblies += A;
    t->last_assemblies = A;
    t->cmax = cmax;
    for (e = t->hist_len - 1; e > 0; e--) {
        COPY(t->qring + e * n, t->qring + (e - 1) * n, n);
        COPY(t->iqring + e * n, t->iqring + (e - 1) * n, n);
    }
    COPY(t->qring, q3, n);
    for (i = 0; i < n; i++) t->iqring[i] = -(i3[i] + Uj[3 * n + i]);
    for (i = 0; i < 3; i++) xv[i] = Y + i * n;
    xv[3] = Y + 2 * n;
    note(t, tn, xv, ord);
    COPY(t->x, Y + 2 * n, n);
    reduce_vec(t->x, t->xr, n, t->iref);
    COPY(t->qn, q3, n);
    reduce_mat(rs->Ca + (A - 1) * 3 * nn + 2 * nn, t->Cn, n, t->iref);
    COPY(t->Xout + j * t->m, t->xr, t->m);
    COPY(t->Yout + j * 3 * n, Y, 3 * n);
}

/* any entry of the carried columns of the new `P` not finite */
static int nonfinite(const stagewalk_t *t, const double *Pn)
{
    long c, i, M = 3 * t->m;
    for (c = 0; c < t->k; c++)
        for (i = 0; i < t->m; i++)
            if (!isfinite(Pn[c * M + i])) return 1;
    return 0;
}

/* the live columns (`advance`): the new P `Pn`, then every few steps the
   columns that stayed zero dropped (`np.any(P != 0, axis=0)`: a NaN live),
   compacted into `Pc` -- in place: kk <= c; `cols[keep]`, or `keep` itself
   where every column was kept so far */
static void advance(stagewalk_t *t, double *Pn)
{
    long m = t->m, M = 3 * m, c, i, kk = 0, nadv;
    t->P = Pn;
    t->ld = M;
    t->trans = 1;
    nadv = ++t->nadv;
    if (!(t->live_on && (nadv <= 4 || nadv % 64 == 0))) return;
    for (c = 0; c < t->k; c++) {
        int live = 0;
        for (i = 0; i < m && !live; i++) live = Pn[c * M + i] != 0.0;
        if (live) {
            COPY(t->Pc + kk * m, Pn + c * M, m);
            t->cols[kk] = t->have_cols ? t->cols[c] : c;
            kk++;
        }
    }
    if (kk < t->k) {
        t->have_cols = 1;
        t->k = kk;
        t->P = t->Pc;
        t->ld = m;
    }
}

/* the walk (contract: stagewalk.h): each step, its stops in their order --
   before the Newton nothing of the step made; after it converged, its
   result left in the radau record (the stage KLU's numeric moved with it);
   after the step map the step made and `j` past it */
long pcx_stagewalk(stagewalk_t *t)
{
    radau_t *rs = (radau_t *) t->radau;
    long n = t->n, m = t->m, st;
    t->klu_refactors = t->klu_factors = t->klu_pivot_fallbacks = t->klu_residual_fallbacks = 0;
    t->walk_refactors = t->walk_fallbacks = t->walk_factors = 0;
    t->assemblies = 0;
    for (; t->j < t->steps; t->j++) {
        long j = t->j, A, ord[4];
        const double *ts = t->tst + 3 * j, *Uj = t->U + 4 * n * j, *Y;
        double tn[4], cmax = t->cmax, *Pn;
        /* (the step map's buffer: the one the carried P is not in, so a stop
           after the map leaves P as the step found it) */
        double *Xw = t->P >= t->Xs && t->P < t->Xs + 3 * m * m ? t->Xs2 : t->Xs;
        /* (a walk with no live column: `_sensklu_c.step` declines it) */
        if (t->k < 1) return t->status = STAGEWALK_BAIL_LIVE;
        if ((st = order_nodes(t, ts, t->tend[j], tn, ord))) return t->status = st;
        if ((st = seeds(t, ts))) return t->status = st;
        /* the C stage Newton, as `_tran_radau_c.solve` calls it -- a
           hand-back leaves the numeric it was given (the glue's reading) */
        rs->qn = t->qn;
        rs->u = Uj;
        rs->seed = t->seed;
        rs->h = t->hs[j];
        rs->klu.klu_num = t->klu_out;
        if (CALL(t->radau_fn, pcx_radau, rs) != RADAU_CONVERGED)
            return t->status = STAGEWALK_BAIL_RADAU;
        t->klu_out = rs->klu.klu_out;
        if (!rs->fin_ok) return t->status = STAGEWALK_BAIL_FINAL;
        A = rs->assemblies;
        Y = rs->Ya + (A - 1) * 3 * n;
        /* (the next step starts from the reduced point with +0.0 inserted:
           this one's reference entry that, bit for bit) */
        if (Y[2 * n + t->iref] != 0.0 || signbit(Y[2 * n + t->iref]))
            return t->status = STAGEWALK_BAIL_REF;
        if ((st = screen(t, rs, A, &cmax))) return t->status = st;
        st = step_map(t, rs, Xw);
        /* (a hand-back commits nothing of the step map: `block_step`'s) */
        if (st == SENSKLU_BAIL_FLAGS) return t->status = STAGEWALK_BAIL_SENS;
        commit(t, rs, A, cmax, Uj, tn, ord, j);
        if (st == SENSKLU_SINGULAR || st == SENSKLU_SOLVE_FAILED) {
            /* (the Python's LinAlgError, raised by the glue) */
            t->j++;
            return t->status = st == SENSKLU_SINGULAR ? STAGEWALK_SINGULAR
                                                      : STAGEWALK_SOLVE_FAILED;
        }
        Pn = Xw + 2 * m;
        if (t->have_cols && nonfinite(t, Pn)) {
            /* (`advance`'s full-width step, the glue's) */
            t->j++;
            return t->status = STAGEWALK_NONFINITE;
        }
        advance(t, Pn);
    }
    return t->status = STAGEWALK_DONE;
}
