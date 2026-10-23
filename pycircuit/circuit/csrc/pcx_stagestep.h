/* pcx_stagestep.h -- the pieces of a Radau stage step both C stage drivers
 * make the same way: the shooting's plain stage walk (stagewalk.c) and the
 * plain run's stepping loop (loop.c's stage family) -- the step's nodes in
 * Python's stable-sort order, the three seeds by the predictor, the screen
 * of each stage, `_pred_note`'s merge of the step's nodes.  Static inline:
 * each driver compiles the one text.  Include after pcx_minmax.h, pred.h
 * and branch.h. */
#ifndef PCX_STAGESTEP_H
#define PCX_STAGESTEP_H
#include <string.h>
#include <math.h>

/* a call to a unit through its pointer, typed by the unit's own declaration */
#define PCX_SCALL(fn, entry, p) (((__typeof__(&entry)) (fn))(p))

/* the step's nodes -- its stage times `ts` and its end `tj` -- into `tn`,
   and the order Python's stable sort takes them in into `ord`: 1, or 0
   where the first is not after the history's newest (`ht[hk - 1]`: the
   general merge, the Python's) */
static inline int pcx_stage_order(const double *ts, double tj, long hk, const double *ht,
                                  double *tn, long *ord)
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
    return !(hk && !(tn[ord[0]] >= ht[hk - 1]));
}

/* the seeds, `_pred_or(np.array(x), t)` at each stage time `ts`, into
   `seed` (3 x n), the predictor's state `pr` reading `hk` nodes: the
   number made by its fast path (`pred:fast`), or -1 where the predictor's
   general path makes one */
static inline long pcx_stage_seeds(void *pred_fn, pred_t *pr, long hk, const double *ts,
                                   const double *x, double *seed, long n)
{
    long i, st, fast = 0;
    for (i = 0; i < 3; i++) {
        pr->hk = hk;
        pr->ttarget = ts[i];
        st = PCX_SCALL(pred_fn, pcx_pred, pr);
        if (st == PRED_DONE) {
            memcpy(seed + i * n, pr->out, (size_t) n * sizeof(double));
            fast++;
        } else if (st == PRED_NONE) {
            memcpy(seed + i * n, x, (size_t) n * sizeof(double));
        } else {
            return -1;
        }
    }
    return fast;
}

/* the screen of the three stages' C (`C3`: 3 x n x n), the running maximum
   carried in `*cmax` (written only where all three are QUIET): 1, or 0 at
   a stage not QUIET; nothing screened where `br` is NULL */
static inline int pcx_stage_screen(void *branch_fn, branch_t *br, const double *C3, long n,
                                   double *cmax)
{
    long i, nn = n * n;
    double cm = *cmax;
    if (!br) return 1;
    for (i = 0; i < 3; i++) {
        br->C = C3 + i * nn;
        br->cmax = cm;
        if (PCX_SCALL(branch_fn, pcx_branch, br) != BRANCH_QUIET) return 0;
        cm = br->ref;
    }
    *cmax = cm;
    return 1;
}

/* `_pred_note`'s merge of the step's four nodes (`tn`, the states `xv`, in
   the order `ord`), each after the history's newest: appended, one within
   `1e-14 * max(abs(t), 1.0)` of the newest replacing it, the newest
   PRED_HIST_MAX kept (an append past them drops the oldest) */
static inline void pcx_stage_note(long *hk, double *ht, double *hx, long n, const double *tn,
                                  const double *const *xv, const long *ord)
{
    long e;
    for (e = 0; e < 4; e++) {
        double te = tn[ord[e]];
        double mx = pcx_py_max(fabs(te), 1.0);
        if (*hk && fabs(te - ht[*hk - 1]) <= 1e-14 * mx) {
            (*hk)--;
        } else if (*hk == PRED_HIST_MAX) {
            memmove(ht, ht + 1, (PRED_HIST_MAX - 1) * sizeof(double));
            memmove(hx, hx + n, (size_t) ((PRED_HIST_MAX - 1) * n) * sizeof(double));
            *hk = PRED_HIST_MAX - 1;
        }
        ht[*hk] = te;
        memcpy(hx + *hk * n, xv[ord[e]], (size_t) n * sizeof(double));
        (*hk)++;
    }
}

#endif
