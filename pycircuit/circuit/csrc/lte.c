
/* lte.c -- the adaptive error test in C: `StepController._max_error`'s
   chain on one attempt (`_charge_lte`, `_normalised`, `np.max`), bit for
   bit, or handed back.  Its glue is `_tran_lte_c`; the C loop (loop.c)
   calls it too.

   It tests floating-point flags under THE FLAG RULE (pcx_fenv.h, MUST
   KEEP): each operation a test covers stores its result through a buffer
   of the struct, or decides a branch taken before the test.  Its chain
   runs partly under `np.errstate(divide='ignore', invalid='ignore')`,
   so a flag there hands back where numpy would not warn.  Every return
   leaves the flags clear, as the chain's last numpy call leaves them.
   The early abort's estimate (`lu_mode`) tests none: its Python
   (`_StepNewton._early_estimate`) runs under `np.errstate(all='ignore')`,
   so its values are the IEEE ones whatever a flag says, and it warns
   about nothing. */
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <fenv.h>
#include "pcx_abi.h"
#include "pcx_minmax.h"
#include "pcx_fenv.h"
#include "pcx_dense.h"
#include "pcx_verdict.h"
#include "lte.h"

/* the flags numpy checks after its operations (inexact is not one) */

/* every bail: the flags cleared, `why` in `status` and returned */
static long lte_bail(lte_t *s, long why)
{
    feclearexcept(FE_ALL_EXCEPT);
    s->status = why;
    return why;
}

/* StructuralKLU.solve(the reduced J, the reduced Eg) on its refactor path
   (`_ours`, `klu_refactor` on the run's numeric, `klu_solve`, the
   residual's validation; `on_pattern`: the values gathered from J, the
   validation's verdict by the sparse sum's bound where that decides --
   lte.h): 0 with the solution in `b` (m), or why it is handed back -- LTE_BAIL_FOREIGN (a nonzero off the pattern: the
   Python's `_ours` sends the matrix to its fallback solver), LTE_BAIL_KLU
   (the refactor or the solve failed: the Python factors afresh),
   LTE_BAIL_FLAGS, LTE_BAIL_RESIDUAL (the validation failed or is NaN: the
   Python factors afresh).  Assumes the sparse block bound (lte.h) and
   `lte` holding Eg.  Writes `R`, `kAx`, `x`, `y`, `b`, and refactors
   `klu_num` in place: on a hand-back after the refactor the Python's own
   solve refactors it again, from the same values.  Flags: what KLU raised
   is cleared; on 0 the flags hold at most what the last division raised
   (pcx_lte's next test sees it), on a bail the caller clears them. */
static long sparse_solve(lte_t *s, long n, long r)
{
    long m = n - 1, k, i, nz_r = 0, nz_a = 0;
    double max_b = 0.0, max_r = 0.0, scale;
    /* (no finiteness test of J here: a non-finite entry off the pattern is a
       nonzero the count hands back, and one on it reaches the residual,
       which it makes non-finite -- handed back too) */
    if (s->on_pattern) {
        /* (J zero off the pattern: `_ours`'s count holds) */
        for (k = 0; k < s->knnz; k++) s->kAx[k] = s->J[s->kat[k]];
    } else {
        pcx_reduce_rm(s->J, s->R, n, r);
        /* `_ours`: the pattern's values, and as many nonzeros as the matrix
           (`v != 0.0` counts a NaN, as numpy's count_nonzero does) */
        nz_r = pcx_count_nonzero(s->R, m * m);
        for (k = 0; k < s->knnz; k++) {
            double v = s->R[s->ktake[k]];
            s->kAx[k] = v;
            nz_a += v != 0.0;
        }
        if (nz_r != nz_a) return LTE_BAIL_FOREIGN;
    }
    pcx_drop_ref(s->lte, s->b, n, r);
    for (i = 0; i < m; i++) s->x[i] = s->b[i];
    if (!((klu_refactor_t) s->klu_refactor)(s->kAp, s->kAi, s->kAx, s->klu_sym, s->klu_num,
                                            s->klu_common))
        return LTE_BAIL_KLU;
    if (!((klu_solve_t) s->klu_solve)(s->klu_sym, s->klu_num, (int) m, 1, s->x, s->klu_common))
        return LTE_BAIL_KLU;
    /* (what KLU raised inside: the Python's ctypes call never sees it) */
    feclearexcept(FE_ALL_EXCEPT);
    if (s->on_pattern) {
        /* the validation's verdict by the bound where it decides and its
           sums raised no flag; else the dense test on R as the Python holds
           it (zero off the pattern), what the sums raised cleared */
        int v = pcx_resid_sparse(s->kAp, s->kAi, s->kAx, s->resid_tol, s->b, s->x, s->y, s->vt,
                                 s->vt + n, m);
        if (v != PCX_UNDECIDED && !fetestexcept(PCX_FP_FLAGS)) {
            if (v == PCX_FAILS) return LTE_BAIL_RESIDUAL;
            for (i = 0; i < m; i++) s->b[i] = s->x[i];
            return 0;
        }
        feclearexcept(FE_ALL_EXCEPT);
        for (i = 0; i < m * m; i++) s->R[i] = 0.0;
        for (k = 0; k < s->knnz; k++) s->R[s->ktake[k]] = s->kAx[k];
        s->vdense++;
    }
    /* the refactor's validation, numpy's operations in their order: R x
       is `arr @ x` (numpy's matmul calls this dgemv), then |R x - b|, its
       maximum, the scale */
    ((dgemv64_t) s->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, m, m, 1.0, s->R, m, s->x, 1,
                           0.0, s->y, 1);
    /* (np.max of |b|: b is finite -- the estimate passed its flag test --
       and no |.| is -0.0, so a running `>` from +0.0 is numpy's maximum
       in any order) */
    for (i = 0; i < m; i++) {
        double a = fabs(s->b[i]);
        if (a > max_b) max_b = a;
    }
    /* (numpy's max propagates a NaN, which a running `>` would drop: a NaN
       is handed back -- the Python's NaN residual fails its test and it
       factors afresh) */
    for (i = 0; i < m; i++) {
        double d = fabs(s->y[i] - s->b[i]);
        if (isnan(d)) return LTE_BAIL_RESIDUAL;
        if (d > max_r) max_r = d;
    }
    if (fetestexcept(PCX_FP_FLAGS)) return LTE_BAIL_FLAGS;
    /* Python's `max(max_b, 1e-300)`: the first argument kept unless the
       second is strictly greater; `not (resid <= tol)` fails a NaN too */
    scale = pcx_py_max(max_b, 1e-300);
    if (!(max_r / scale <= s->resid_tol)) return LTE_BAIL_RESIDUAL;
    for (i = 0; i < m; i++) s->b[i] = s->x[i];
    return 0;
}

/* every input finite, and a running maximum not negative: anything else
   is the Python chain's.  1, or 0 */
static int inputs_finite(const lte_t *s, long n)
{
    long i;
    if (!(isfinite(s->h1) && isfinite(s->h2) && isfinite(s->h3) && isfinite(s->s12)
          && isfinite(s->s23) && isfinite(s->s123) && isfinite(s->coef)
          && isfinite(s->reltol) && isfinite(s->trtol) && isfinite(s->abstol)))
        return 0;
    for (i = 0; i < n; i++)
        if (!(isfinite(s->xc[i]) && isfinite(s->xl[i]) && isfinite(s->q[i])
              && isfinite(s->q1[i]) && isfinite(s->q2[i]) && isfinite(s->q3[i])))
            return 0;
    /* (a running maximum also not negative, -0.0 included: a -0.0 could
       reach sigglobal's np.max, whose pairing of +0.0 and -0.0 -- numpy's
       reduction is SIMD-paired -- a running `>` does not follow) */
    if (s->have_run)
        for (i = 0; i < n; i++)
            if (!isfinite(s->run[i]) || signbit(s->run[i])) return 0;
    if (s->abs_vec)
        for (i = 0; i < n; i++) if (!isfinite(s->abst[i])) return 0;
    if (s->formula == LTE_FORMULA_EULER)
        for (i = 0; i < n; i++) if (!isfinite(s->iq1[i])) return 0;
    return 1;
}

/* compute_lte into `lte`, under a cleared-and-tested region of its own.
   0, or LTE_BAIL_FLAGS */
static long charge_lte(lte_t *s, long n)
{
    long i;
    feclearexcept(FE_ALL_EXCEPT);
    /* compute_lte: coef * third_divided_difference, numpy's elementwise
       operations in their order (each result stored, so each operation
       precedes the flag test -- the flag rule, above); the step sums and
       the coefficient are the caller's (the glue's, or loop.c's), in the
       integrator's own scalar arithmetic.  The coefficient multiplies last: the Python makes the
       scalar `coef` first (`-h * (h + hl)`, `-(h**2)`, `-0.5 * scale`) and
       multiplies the finished vector by it */
    if (s->formula == LTE_FORMULA_EULER) {
        /* EulerIntegrator.compute_lte: gn = (q - q1) / h, gn - iq1, scaled */
        for (i = 0; i < n; i++) s->lte[i] = s->coef * ((s->q[i] - s->q1[i]) / s->h1 - s->iq1[i]);
    } else for (i = 0; i < n; i++) {
        /* third_divided_difference: d1..d3, dd_a, dd_b, then over the sum
           of the three steps -- each a numpy temporary, rounded as here */
        double d1 = (s->q[i] - s->q1[i]) / s->h1;
        double d2 = (s->q1[i] - s->q2[i]) / s->h2;
        double d3 = (s->q2[i] - s->q3[i]) / s->h3;
        double dda = (d1 - d2) / s->s12;
        double ddb = (d2 - d3) / s->s23;
        s->lte[i] = s->coef * ((dda - ddb) / s->s123);
    }
    if (fetestexcept(PCX_FP_FLAGS) && !s->lu_mode) return LTE_BAIL_FLAGS;
    return 0;
}

/* J^-1 Eg by the factors in hand (`lu_mode`, the early abort's): the
   reduced Eg back-substituted -- dgetrs on dgesv's LU is numpy's solve of
   that matrix, on dgetrf's SciPy's lu_solve, klu_solve the chord's
   `_KLUNumeric.solve`.  0, or why it is handed back (an info, a KLU failure:
   the flags cleared, as numpy's solve and the ctypes call leave them) */
static long factored_solve(lte_t *s, long n, long r)
{
    long m = n - 1;
    char tr = 'N';
    pcx_drop_ref(s->lte, s->b, n, r);
    if (s->lu_mode == LTE_LU_GETRS64) {
        int64_t m64 = m, one64 = 1, info64 = 0;
        ((dgetrs64_t) s->lu_dgetrs64)(&tr, &m64, &one64, s->luA, &m64, s->lu_ipiv64, s->b, &m64,
                                      &info64, 1);
        feclearexcept(FE_ALL_EXCEPT);
        if (info64 != 0) return LTE_BAIL_DGESV;
    } else if (s->lu_mode == LTE_LU_GETRS32) {
        int m32 = (int) m, one32 = 1, info32 = 0;
        ((dgetrs32_t) s->lu_dgetrs32)(&tr, &m32, &one32, s->luA, &m32, s->lu_ipiv32, s->b, &m32,
                                      &info32, 1);
        feclearexcept(FE_ALL_EXCEPT);
        if (info32 != 0) return LTE_BAIL_DGESV;
    } else {
        int ok = ((klu_solve_t) s->lu_klu_solve)(s->lu_sym, s->lu_num, (int) m, 1, s->b,
                                                 s->lu_common);
        feclearexcept(FE_ALL_EXCEPT);
        if (!ok) return LTE_BAIL_KLU;
    }
    return 0;
}

/* J^-1 Eg into `lte`, the reference row cut and restored as +0.0: by the
   factors in hand (`lu_mode`, `factored_solve`); dense, numpy's own dgesv
   on the reduced J; sparse, `sparse_solve`.  0, or why it is handed back */
static long solve_lte(lte_t *s, long n, long r)
{
    long m = n - 1, i;
    int64_t m64 = m, one64 = 1, info64 = 0;
    /* remove_row_col and numpy's solve: the reduced J column-major (numpy's
       copy for LAPACK), the reduced Eg, numpy's own dgesv */
    if (s->lu_mode) {
        long why = factored_solve(s, n, r);
        if (why) return why;
    } else if (s->klu_sym) {
        long why = sparse_solve(s, n, r);
        if (why) return why;
    } else {
        if (!pcx_reduce_cm(s->J, s->A, n, r)) return LTE_BAIL_NONFINITE;
        pcx_drop_ref(s->lte, s->b, n, r);
        ((dgesv64_t) s->dgesv)(&m64, &one64, s->A, &m64, s->ipiv, s->b, &m64, &info64);
        /* (numpy's solve clears what LAPACK raised; an info is its LinAlgError) */
        feclearexcept(FE_ALL_EXCEPT);
        if (info64 != 0) return LTE_BAIL_DGESV;
    }
    /* (a non-finite J^-1 Eg is the Python's: its normalised_error makes
       it +inf, a rejection) */
    for (i = 0; i < m; i++) if (!isfinite(s->b[i])) return LTE_BAIL_LTE;
    /* the full vector, the reference row restored as +0.0 (the Python's
       concatenate of the two halves around [0.0]) */
    pcx_put_ref(s->b, s->lte, n, r);
    return 0;
}

/* `_reference` by `relref` into `ref`, folded into the running maximum
   `out` */
static void reference(lte_t *s, long n)
{
    long i;
    /* _reference: numpy's maximum (the first where they compare equal) of
       |x_curr| and |x_last|, folded into the running maximum.  The ternary
       is numpy's maximum, `(a >= b || isnan(a)) ? a : b`, on finite values
       (fmax would not promise which of two equal zeros it keeps) */
    for (i = 0; i < n; i++) {
        double a = fabs(s->xc[i]), c = fabs(s->xl[i]);
        double loc = pcx_np_maximum(a, c);
        if (s->relref == LTE_RELREF_POINTLOCAL) {
            s->ref[i] = loc;
        } else if (s->have_run) {
            double rv = s->run[i];
            s->out[i] = pcx_np_maximum(rv, loc);
        } else {
            s->out[i] = loc;
        }
    }
    if (s->relref == LTE_RELREF_ALLLOCAL) {
        for (i = 0; i < n; i++) s->ref[i] = s->out[i];
    } else if (s->relref == LTE_RELREF_SIGGLOBAL) {
        /* sigglobal_reference: each unit group's maximum, broadcast back --
           np.max over values with no NaN and no -0.0 (both excluded
           above), so a running `>` is its answer in any order; one group
           where n_nodes is None, <= 0 or >= n, as the Python splits */
        long nn = s->nn, lo, hi, g;
        long cut = (nn <= 0 || nn >= n) ? n : nn;
        for (g = 0; g < 2; g++) {
            double mx;
            lo = g ? cut : 0;
            hi = g ? n : cut;
            if (lo >= hi) continue;
            mx = s->out[lo];
            for (i = lo + 1; i < hi; i++) if (s->out[i] > mx) mx = s->out[i];
            for (i = lo; i < hi; i++) s->ref[i] = mx;
        }
    }
}

/* `charge_floor`, `tolerance`, `normalised_error`, `_normalised`'s mask and
   `np.max`: `err`.  0, or why it is handed back */
static long normalised_error(lte_t *s, long n)
{
    long i;
    double err, ab, floor = 0.0;
    /* `charge_floor` (`lte_algebraic='charge'`): on the rows `charge` marks
       the absolute tolerance is np.maximum(abstol, chgtol / h) -- the
       Python's one division, h its `h_curr` -- elsewhere abstol.  Then the
       tolerance TRTOL (reltol ref + abstol) and |lte| / tol, in
       `tolerance`'s order.  normalised_error's two `where`s (0/0 -> 0, NaN
       -> inf) cannot fire past the tests: every tol positive and finite,
       every |lte| finite */
    if (s->has_charge) floor = s->chgtol / s->h1;
    for (i = 0; i < n; i++) {
        ab = s->abs_vec ? s->abst[i] : s->abstol;
        if (s->has_charge && s->charge[i] != 0.0) ab = pcx_np_maximum(ab, floor);
        s->tq[i] = s->trtol * (s->reltol * s->ref[i] + ab);
    }
    if (fetestexcept(PCX_FP_FLAGS) && !s->lu_mode) return LTE_BAIL_FLAGS;
    for (i = 0; i < n; i++) if (!(s->tq[i] > 0.0 && isfinite(s->tq[i]))) return LTE_BAIL_TOLERANCE;
    for (i = 0; i < n; i++) s->tq[i] = fabs(s->lte[i]) / s->tq[i];
    if (fetestexcept(PCX_FP_FLAGS) && !s->lu_mode) return LTE_BAIL_FLAGS;
    /* `_normalised`'s mask (`Transient.lte_algebraic`): an unknown not
       judged is +0.0, as `np.where(mask, e, 0.0)` makes it (a mask entry
       is true where nonzero: -0.0 false, NaN true, as `== 0.0` decides) */
    if (s->has_mask)
        for (i = 0; i < n; i++) if (s->mask[i] == 0.0) s->tq[i] = 0.0;
    /* np.max: no NaN, no negative zero (|lte| / a positive tolerance) */
    err = s->tq[0];
    for (i = 1; i < n; i++) if (s->tq[i] > err) err = s->tq[i];
    s->err = err;
    return 0;
}

/* `StepController._max_error` on one attempt: `compute_lte` by `formula`,
   `J^-1 Eg` with the reference row cut and restored (dense: numpy's
   `dgesv`; sparse: `sparse_solve`), `_reference` by `relref` folded into
   the running maximum, `tolerance`, `normalised_error`, `_normalised`'s
   mask and `np.max` -- numpy's operations in numpy's order.
   Assumes: the struct as lte.h says it (n >= 2, 0 <= iref < n, every
   pointer valid for its size; sparse, a bound StructuralKLU whose numeric
   exists).
   LTE_DONE: `err` the largest normalised error, `out` the new running
   maximum (relref not POINTLOCAL).
   A bail (enum lte_status): `err` is not written; `out` and the work
   buffers may hold part of the work -- no caller reads them on a bail --
   and the inputs are only read.  Every return leaves the flags clear. */
long pcx_lte(lte_t *s)
{
    long n = s->n, r = s->iref, why;
    if (!inputs_finite(s, n)) return lte_bail(s, LTE_BAIL_NONFINITE);
    if ((why = charge_lte(s, n)) != 0) return lte_bail(s, why);
    if ((why = solve_lte(s, n, r)) != 0) return lte_bail(s, why);
    reference(s, n);
    if ((why = normalised_error(s, n)) != 0) return lte_bail(s, why);
    feclearexcept(FE_ALL_EXCEPT);
    s->status = LTE_DONE;
    return LTE_DONE;
}
