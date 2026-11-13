/* pcx_kluseq.h -- StructuralKLU.solve's own sequence (`linearsolver`),
 * written once for the units that repeat it (dcnewton.c, pcnrdc.c): a
 * refactor on the numeric it holds, a fresh factor on a zero pivot,
 * klu_solve, and after a refactor the residual `max |A x - b| / max(max |b|,
 * 1e-300)` against REFACTOR_RESIDUAL_TOL -- a fresh factor where it fails.
 * The struct and THE NUMERIC: pcx_klu.h, the verdicts pcx_verdict.h, THE
 * HELD VALUES pcx_held.h (all included before this).
 *
 * A step returns 0, or what stopped it: PCX_KLU_FACTOR (a fresh factor
 * failed: `_factor`'s LinAlgError('Singular matrix'), `klu_out` NULL),
 * PCX_KLU_SOLVE (klu_solve failed: its LinAlgError('klu_solve failed')),
 * PCX_KLU_FLAGS (a flag numpy reports in the residual's arithmetic).
 *
 * The residual's test needs only its verdict: `A x` by the CSC kernel on
 * KLU's own arrays, one order on every layout (`pcx_klu_resid_csc`,
 * `StructuralKLU._csc_residual`); where the caller asks (`bound`) a sparse
 * sum with its rounding bound first, which decides it where it does
 * whatever the order (`pcx_klu_resid_sparse`: its sums with subnormals
 * flushed, pcx_verdict.h). */
#include <float.h>
#include <math.h>
#include <fenv.h>

/* the flags numpy reports: underflow is not one (numpy ignores it by
   default; a unit's glue declines where it does not) */
#define PCX_NP_FLAGS (FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

enum pcx_klu_stop { PCX_KLU_FACTOR = 1, PCX_KLU_SOLVE = 2, PCX_KLU_FLAGS = 3 };



static void pcx_klu_free(pcx_klu_t *k, void *num)
{
    if (num) ((klu_free_numeric_t) k->klu_free)(&num, k->klu_common);
}

/* a fresh factor, as StructuralKLU's `_factor` after freeing its numeric
   (the entry's kept: `pcx_klu_finish`) */
static long pcx_klu_fresh(pcx_klu_t *k)
{
    if (k->klu_out != k->klu_num) pcx_klu_free(k, k->klu_out);
    k->klu_out = ((klu_factor_t) k->klu_factor)(k->kAp, k->kAi, k->kAx, k->klu_sym,
                                                k->klu_common);
    if (!k->klu_out) return PCX_KLU_FACTOR;
    k->factors++;
    return 0;
}

/* klu_solve of `b` into `x` */
static long pcx_klu_solve_b(pcx_klu_t *k, const double *b, double *x, long m)
{
    long i;
    for (i = 0; i < m; i++) x[i] = b[i];
    if (!((klu_solve_t) k->klu_solve)(k->klu_sym, k->klu_out, (int) m, 1, x, k->klu_common))
        return PCX_KLU_SOLVE;
    return 0;
}

/* the residual's test (`StructuralKLU._csc_residual`, every layout's): `A x`
   by the CSC kernel on KLU's own arrays -- each row from 0.0, each product
   rounded, then added in KLU's CSC order (np.bincount's) --, `- rhs`, |.|,
   the maxima -- numpy's max propagates a NaN: so does `nan_r` -- against the
   tolerance.  1 it fails (a fresh factor), 0 it holds, -1 a flag (the
   kernel's sums' too, which numpy's bincount does not report: a hand-back,
   the safe side) */
static int pcx_klu_resid_csc(pcx_klu_t *k, const double *b, const double *x, double *y,
                             long m)
{
    long i;
    double max_b = 0.0, max_r = 0.0;
    int nan_r = 0;
    pcx_csc_mv(k->kAp, k->kAi, k->kAx, x, y, m, 0);
    for (i = 0; i < m; i++) y[i] = y[i] - b[i];
    if (fetestexcept(PCX_NP_FLAGS)) return -1;
    for (i = 0; i < m; i++) {
        double a = fabs(b[i]), d = fabs(y[i]);
        if (a > max_b) max_b = a;
        if (isnan(d)) nan_r = 1;
        else if (d > max_r) max_r = d;
    }
    return nan_r || !(max_r / pcx_py_max(max_b, 1e-300) <= k->resid_tol);
}

/* the same verdict on KLU's CSC (A zero off its pattern): pcx_verdict.h's
   `pcx_resid_sparse`.  `t` 2 m scratch */
static int pcx_klu_resid_sparse(pcx_klu_t *k, const double *b, const double *x,
                                double *y, double *t, long m)
{
    return pcx_resid_sparse(k->kAp, k->kAi, k->kAx, k->resid_tol, b, x, y, t, t + m, m);
}

/* StructuralKLU.solve(A, b) on its own sequence, `kAx` A's values in its
   order (gathered by the caller); `bound` 1 the sparse sums' bound first,
   its undecided verdicts counted in `dense` (which may be NULL) -- the CSC
   kernel's test then, at once where `bound` is 0.  The solution in `x`, `y`
   (m) and `t` (2 m, read only where the bound is tried) scratch.  `kept`: `solve_kept`'s -- the numeric
   holds a refactor of these very values, none is made (nor counted).  `hv`, `hst`: StructuralKLU's held
   values (pcx_held.h) where `k` is its run numeric -- no refactor of the values held, counted there; NULL
   `hst` (the chord's numeric) the rule off */
static long pcx_klu_seq(pcx_klu_t *k, const double *b, double *x, double *y, double *t, long m,
                        int bound, long *dense, int kept, double *hv, int64_t *hst)
{
    long st, reused = k->klu_out != NULL, nnz = k->kAp[m];
    int fails = PCX_UNDECIDED;
    if (reused && !kept && pcx_held_same(hst, hv, k->kAx, nnz)) {
        /* (THE SAME VALUES: the numeric is their refactor already) */
    } else if (reused && !kept) {
        int ok = ((klu_refactor_t) k->klu_refactor)(k->kAp, k->kAi, k->kAx, k->klu_sym,
                                                  k->klu_out, k->klu_common);
        k->refactors++;
        pcx_held_refactored(hst, hv, k->kAx, nnz, ok);
        if (!ok) {
            k->pivot_fallbacks++;
            pcx_held_lost(hst, hv, nnz);
            if ((st = pcx_klu_fresh(k)) != 0) return st;
            reused = 0;
        }
    } else if (!reused) {
        pcx_held_lost(hst, hv, nnz);
        if ((st = pcx_klu_fresh(k)) != 0) return st;
    }
    if ((st = pcx_klu_solve_b(k, b, x, m)) != 0) return st;
    /* (what KLU raised inside: the Python's ctypes calls never see it) */
    feclearexcept(FE_ALL_EXCEPT);
    if (!reused) return 0;
    if (bound) {
        fails = pcx_klu_resid_sparse(k, b, x, y, t, m);
        if (fails == PCX_UNDECIDED) {
            /* (what the sparse sums raised is not the Python's) */
            feclearexcept(FE_ALL_EXCEPT);
            if (dense) (*dense)++;
        }
    }
    if (fails == PCX_UNDECIDED) fails = pcx_klu_resid_csc(k, b, x, y, m);
    if (fails < 0) return PCX_KLU_FLAGS;
    /* (THE REFACTOR KEPT ON ITS GROWTH: pcx_rgrowth_keeps) */
    if (fails && pcx_rgrowth_keeps(k->klu_rgrowth, k->rg_off, k->rg_keep, k->kAp, k->kAi, k->kAx,
                                   k->klu_sym, k->klu_out, k->klu_common, x, m)) {
        /* (counted where the Python counts it: the run numeric's held
           state, pcx_held.h, its checkpoint's to restore -- the chord's
           numeric, which holds none, in its struct) */
        if (hst) hst[PCX_HELD_KEPT]++;
        else k->residual_kept++;
        fails = 0;
    }
    if (fails) {
        k->residual_fallbacks++;
        pcx_held_lost(hst, hv, nnz);
        if ((st = pcx_klu_fresh(k)) != 0) return st;
        if ((st = pcx_klu_solve_b(k, b, x, m)) != 0) return st;
        feclearexcept(FE_ALL_EXCEPT);
    }
    return 0;
}

/* StructuralKLU.solve(A, b) (`pcx_klu_seq`) */
static long pcx_klu_step(pcx_klu_t *k, const double *b, double *x, double *y, double *t, long m,
                         int bound, long *dense, double *hv, int64_t *hst)
{
    return pcx_klu_seq(k, b, x, y, t, m, bound, dense, 0, hv, hst);
}

/* `_FrozenTransform.real_solve` on the StructuralKLU: `solve_kept` where
   `*kept` (the numeric holds the refactor of these values an earlier solve
   of them made), `solve` else; `*kept` then as it leaves it -- set where the
   solve ended on such a refactor (no fresh factor), cleared where not.  The
   residual's bound first */
static inline long pcx_klu_kept(pcx_klu_t *k, long *kept, const double *b, double *x, double *y,
                                double *t, long m, long *dense, double *hv, int64_t *hst)
{
    long f0 = k->factors, r0 = k->refactors;
    long st = pcx_klu_seq(k, b, x, y, t, m, 1, dense, *kept != 0, hv, hst);
    *kept = st == 0 && k->factors == f0 && (*kept || k->refactors == r0 + 1);
    return st;
}

/* the numeric as the outcome leaves it (THE NUMERIC, pcx_klu.h): a
   hand-back the entry's (any other freed), a served outcome the
   sequence's (the entry's freed where another replaced it) */
static void pcx_klu_finish(pcx_klu_t *k, int handed_back)
{
    if (handed_back) {
        if (k->klu_out != k->klu_num) pcx_klu_free(k, k->klu_out);
        k->klu_out = k->klu_num;
    } else if (k->klu_out != k->klu_num) {
        pcx_klu_free(k, k->klu_num);
    }
}
