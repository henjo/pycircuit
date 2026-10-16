
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include "pcx_abi.h"
#include "pcx_minmax.h"
#include "pcx_dense.h"
#include "newton.h"
#include "lte.h"

/* THE ITERATION: `StandardNewton.solve_system` (chord 0) or
   `ChordNewton.solve_system` up to its fallback (chord 1), on the reduced
   system `_newton` hands them (`refnode_removed`: each evaluation at
   insert_row(x), the reference row of F and of J's rows and columns
   dropped), and at convergence `jacobian_only` at the converged point.
   Every float operation is the Python's, in its order; where the Python
   would take a course this cannot follow (a non-finite value, a LAPACK or
   KLU failure, a limiter's decline, a chord break, maxiter) it returns a
   bail and the Python repeats the solve.
   Assumes what `pcx_newton` assumes (newton.h).  Returns the status, also left
   in `status`.  NEWTON_CONVERGED: `x` the reduced solution, `iters`, and
   where fold_j is nonzero xc, Cj, qj, iqj, Geqj (and Jj where it has
   PCX_PASS_G).  A bail writes none of xc .. Jj.  Every return: `*num` the
   sparse chord's numeric where it made one (the caller frees it).  Its
   phases come first below, then `iterate`, which runs them. */

/* the result, left in `status` too */
static long finish(newton_t *s, long code)
{
    s->status = code;
    return code;
}

/* `_C_lookup`'s equality test of the full seed it evaluates against the
   cached state: elementwise == (`alltrue(x_cached == x)`, a NaN never
   equal), the full seed being insert_row(x) -- +0.0 at the reference, so
   cx[r] is compared with 0.0.  The cached C where it holds (a hit), else
   NULL */
static const double *cache_hit(const newton_t *s, long m, long r)
{
    long i;
    if (!s->have_cache) return 0;
    for (i = 0; i < r; i++) if (!(s->cx[i] == s->x[i])) return 0;
    if (!(s->cx[r] == 0.0)) return 0;
    for (i = r; i < m; i++) if (!(s->cx[i + 1] == s->x[i])) return 0;
    return s->cC;
}

/* The sparse full Newton's solve (newton.h): `StructuralKLU.solve` on its
   refactor path, in `b` -- klu_refactor of the run's numeric, klu_solve,
   then the validation in numpy's operations and order: `arr @ x` (numpy's
   matmul calls this dgemv), |R x - rhs|, its maximum (a NaN handed back:
   numpy's max propagates one, and the Python's NaN residual fails its
   test), the scale `max(max |rhs|, 1e-300)` (rhs finite: a running `>`
   from +0.0 is numpy's maximum in any order), `not (resid <= tol)`.  Every
   other course is the Python's: no numeric, a non-finite matrix or
   right-hand side, a zero pivot, a failed solve or validation.  0, or the
   bail */
static long sparse_full(newton_t *s, long m)
{
    long i, k;
    double max_b = 0.0, max_r = 0.0, scale;
    if (!s->klu_num) return NEWTON_BAIL_KLU_NONUMERIC;
    for (k = 0; k < s->knnz; k++) if (!isfinite(s->kAx[k])) return NEWTON_BAIL_JNONFINITE;
    for (i = 0; i < m; i++) if (!isfinite(s->b[i])) return NEWTON_BAIL_BNONFINITE;
    for (i = 0; i < m; i++) s->rb[i] = s->b[i];
    if (!((klu_refactor_t) s->klu_refactor)(s->kAp, s->kAi, s->kAx, s->klu_sym, s->klu_num,
                                            s->klu_common))
        return NEWTON_BAIL_KLU_REFACTOR;
    if (!((klu_solve_t) s->klu_solve)(s->klu_sym, s->klu_num, (int) m, 1, s->b, s->klu_common))
        return NEWTON_BAIL_KLU_SOLVE;
    ((dgemv64_t) s->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, m, m, 1.0, s->R, m, s->b, 1,
                           0.0, s->ry, 1);
    for (i = 0; i < m; i++) {
        double a = fabs(s->rb[i]);
        if (a > max_b) max_b = a;
    }
    for (i = 0; i < m; i++) {
        double d = fabs(s->ry[i] - s->rb[i]);
        if (isnan(d)) return NEWTON_BAIL_KLU_RESIDUAL;
        if (d > max_r) max_r = d;
    }
    scale = pcx_py_max(max_b, 1e-300);
    if (!(max_r / scale <= s->resid_tol)) return NEWTON_BAIL_KLU_RESIDUAL;
    return 0;
}

/* The step, in `b` (-F reduced in, the step out), as the Python's solver
   makes it on iteration `it`: the sparse chord (`StructuralKLU.factor`'s
   numeric made on the first, `*num`; `_KLUNumeric.solve` every one), the
   dense chord (SciPy's lu_factor -- dgetrf -- once, lu_solve -- dgetrs, 'N'
   -- every one), the sparse full Newton (`sparse_full`), or the full
   Newton's numpy.linalg.solve (numpy's own dgesv on the column-major copy).
   A non-finite matrix or right-hand side, and any LAPACK info or KLU
   failure, hands back (the Python chord breaks to its fallback; numpy
   raises LinAlgError).  `finite`: the dense J's entries.  0, or the bail */
static long solve_step(newton_t *s, long it, int finite, void **num)
{
    long m = s->n - 1, i, k;
    int m32 = (int) m, one32 = 1, info32 = 0;
    int64_t m64 = m, one64 = 1, info64 = 0;
    char trans = 'N';
    if (s->chord && s->klu_sym) {
        if (it == 0) {
            for (k = 0; k < s->knnz; k++) if (!isfinite(s->kAx[k])) return NEWTON_BAIL_JNONFINITE;
            *num = ((klu_factor_t) s->klu_factor)(s->kAp, s->kAi, s->kAx, s->klu_sym,
                                                  s->klu_common);
            if (!*num) return NEWTON_BAIL_KLU_FACTOR;
        }
        for (i = 0; i < m; i++) if (!isfinite(s->b[i])) return NEWTON_BAIL_BNONFINITE;
        if (!((klu_solve_t) s->klu_solve)(s->klu_sym, *num, m32, 1, s->b, s->klu_common))
            return NEWTON_BAIL_KLU_SOLVE;
    } else if (s->chord) {
        if (it == 0) {
            if (!finite) return NEWTON_BAIL_JNONFINITE;
            ((dgetrf32_t) s->dgetrf)(&m32, &m32, s->A, &m32, s->ipiv32, &info32);
            if (info32 != 0) return NEWTON_BAIL_GETRF;
        }
        for (i = 0; i < m; i++) if (!isfinite(s->b[i])) return NEWTON_BAIL_BNONFINITE;
        ((dgetrs32_t) s->dgetrs)(&trans, &m32, &one32, s->A, &m32, s->ipiv32, s->b, &m32,
                               &info32, 1);
        if (info32 != 0) return NEWTON_BAIL_GETRS;
    } else if (s->klu_sym) {
        return sparse_full(s, m);
    } else {
        ((dgesv64_t) s->dgesv)(&m64, &one64, s->A, &m64, s->ipiv64, s->b, &m64, &info64);
        if (info64 != 0) return NEWTON_BAIL_DGESV;
    }
    return 0;
}

/* `xdiff`: `_newton_limiter`'s closure where an element limits -- `cir.limit`
   on insert_row of x_next and of x, the reference dropped after, then
   `xdiff = x_next - x` (x_next limited in place) -- else the step itself.
   0, or the bail (a limiter declined: the walk stopped) */
static long limit_step(newton_t *s, walk_fn_t walk, long n, long r)
{
    long m = n - 1, i, e;
    if (!walk) {
        for (i = 0; i < m; i++) s->xd[i] = s->b[i];
        return 0;
    }
    pcx_put_ref(s->xn, s->xnf, n, r);
    pcx_put_ref(s->x, s->x0f, n, r);
    e = walk(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK, s->xnf, s->x0f, 0, s->wn);
    if (e < s->wn) return NEWTON_BAIL_WALKSTOP;
    pcx_drop_ref(s->xnf, s->xn, n, r);
    for (i = 0; i < m; i++) s->xd[i] = s->xn[i] - s->x[i];
    return 0;
}

/* nrsolver's convergence test, entry by entry, in its operations and
   order: `I_scale = toolkit.dot(abs(J), ax_next) + aF` (numpy's own dgemv,
   row-major, as np.dot calls it on a C-ordered matrix), `abs(xdiff) <
   reltol * maximum(ax_next, abs(x)) + xtol`, `aF < reltol * I_scale +
   abstol` -- written `!(a < b)` so a NaN fails, as numpy's `<`.  Assumes x
   and x_next finite (checked) */
static int converged(newton_t *s, long m)
{
    long i;
    int conv_x = 1, conv_f = 1;
    for (i = 0; i < m; i++) s->ax[i] = fabs(s->xn[i]);
    ((dgemv64_t) s->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, m, m, 1.0, s->absJ, m,
                           s->ax, 1, 0.0, s->Is, 1);
    for (i = 0; i < m; i++) {
        double aF = fabs(s->Fr[i]);
        double Is = s->Is[i] + aF;
        double axi = s->ax[i], axo = fabs(s->x[i]);
        /* np.maximum(ax_next, abs(x)): both finite here, where this ternary
           is the same value; it would not propagate a NaN as np.maximum
           does */
        double mx = pcx_np_maximum(axi, axo);
        if (!(fabs(s->xd[i]) < s->reltol * mx + s->xtol[i])) conv_x = 0;
        if (!(aF < s->reltol * Is + s->abstol[i])) conv_f = 0;
    }
    return conv_x && conv_f;
}

/* `jacobian_only` at the converged point x: G, C, q -- C and q alone where
   J is unread (fold_j holds the passes).  Cj is both Cin and C: fold_j
   always has the C pass, so the companion reads the C just made; F is not
   written (no PCX_PASS_I) */
static void fold(newton_t *s, core_fn_t core, long n, long r)
{
    pcx_put_ref(s->x, s->xc, n, r);
    core(s->core, s->xc, s->T, s->fold_j | s->pattern, s->formula, s->a0, s->a1, s->a2, s->h,
         s->theta, s->q1, s->q2, s->iq1, s->u, s->Cj, s->Cj, s->qj, s->iqj, s->Geqj, s->F,
         s->Jj);
}

/* the chord's contraction test: `size = float(np.max(np.abs(xdiff)))` --
   np.max propagates a NaN, which a plain `>` scan would skip; from
   -INFINITY, below any |v| -- then `not np.isfinite(size) or not size <
   prev`, as the Python spells it.  1: the chord breaks; else `prev` is the
   size */
static int chord_breaks(const newton_t *s, long m, double *prev, int *have_prev)
{
    long i;
    double size = -INFINITY;
    for (i = 0; i < m; i++) {
        double v = fabs(s->xd[i]);
        if (isnan(v)) { size = v; break; }
        if (v > size) size = v;
    }
    if (!isfinite(size) || (*have_prev && !(size < *prev))) return 1;
    *prev = size;
    *have_prev = 1;
    return 0;
}

/* `_chord_settled`'s dominant component of the increment just taken:
   `int(np.argmax(np.abs(xdiff) / tol))`, tol the update
   test's (`reltol * np.maximum(ax_next, abs(x)) + xtol`), the first of
   equal maxima as the strict `>` scan finds it (every |xd| is finite here:
   the contraction test passed); `d` its increment, for the next iteration */
static long chord_dominant(const newton_t *s, long m, double *d)
{
    long i, k = 0;
    double best = -1.0;
    for (i = 0; i < m; i++) {
        double mx = pcx_np_maximum(fabs(s->xn[i]), fabs(s->x[i]));
        double w = fabs(s->xd[i]) / (s->reltol * mx + s->xtol[i]);
        if (w > best) { best = w; k = i; }
    }
    *d = s->xd[k];
    return k;
}

/* `_chord_settled`: 1 where an iterate the update test accepts is settled
   -- lambda (this increment over the last at the last's dominant component
   `k`, its increment `d`) under 1/2, or under 1 with max |xd|/tol x
   lambda/(1-lambda) under 1; the same operations in the same order */
static int chord_settled(const newton_t *s, long m, long k, double d)
{
    long i;
    double lam, w = 0.0;
    if (d == 0.0) return 1;
    lam = s->xd[k] / d;
    if (lam < 0.5) return 1;
    if (!(lam < 1.0)) return 0;
    for (i = 0; i < m; i++) {
        double mx = pcx_np_maximum(fabs(s->xn[i]), fabs(s->x[i]));
        double v = fabs(s->xd[i]) / (s->reltol * mx + s->xtol[i]);
        if (v > w) w = v;
    }
    return w * (lam / (1.0 - lam)) < 1.0;
}

/* One evaluation at insert_row(x) (`_tran_core.evaluate`): the
   iteration's buffers (the first's kept apart -- the chord's companion and
   the write-back read them), the full x, the core's passes -- every one
   where J is formed (`full`), else the chord's `residual_only` -- then the
   right-hand side `-F_s` and, formed, J reduced for the solve and |J| for
   the test (`finite`: every entry of J finite).  0, or the status to stop
   with */
/* THE EARLY ABORT's check (newton.h): the error test on the first iterate
   -- `xf`, its q -- by the factors the solve holds, the record's x, q and
   factorisation put back after the one call.  1 over `early_th` (`early_err`
   set), else 0 (a hand-back of the test is no abort: the Python's
   `_early_estimate` makes none there either) */
static int early(newton_t *s, const double *q, void *num)
{
    lte_t *le = (lte_t *) s->early_lte;
    const double *xc0 = le->xc, *q0 = le->q;
    long st;
    le->xc = s->xf;
    le->q = q;
    le->lu_mode = s->chord ? (s->klu_sym ? LTE_LU_KLU : LTE_LU_GETRS32) : LTE_LU_GETRS64;
    le->luA = s->A;
    le->lu_ipiv64 = s->ipiv64;
    le->lu_ipiv32 = s->ipiv32;
    le->lu_num = num;
    le->lu_sym = s->klu_sym;
    le->lu_common = s->klu_common;
    le->lu_klu_solve = s->klu_solve;
    le->lu_dgetrs64 = s->dgetrs64;
    le->lu_dgetrs32 = s->dgetrs;
    st = ((__typeof__(&pcx_lte)) s->early_lte_fn)(le);
    le->lu_mode = LTE_LU_NONE;
    le->xc = xc0;
    le->q = q0;
    if (st != LTE_DONE) return 0;
    s->early_checks = 1;
    if (!(le->err > s->early_th)) return 0;
    s->early_err = le->err;
    return 1;
}

static long evaluate(newton_t *s, long it, int full, long first_bits, const double *Cin0,
                     int *finite, void *num)
{
    long n = s->n, m = n - 1, r = s->iref, i, k;
    core_fn_t core = (core_fn_t) s->core_fn;
    const double *Cin;
    long bits;
    double *F, *J, *C, *q, *iq, *Geq;
    if (it == 0) { F = s->F0; J = s->J0; C = s->C0; q = s->q0; iq = s->iq0; Geq = s->Geq0; }
    else { F = s->F; J = s->J; C = s->C; q = s->q; iq = s->iq; Geq = s->Geq; }
    /* insert_row: the full x of this evaluation; a non-finite one hands
       back (`_tran_core.evaluate` declines one: the Python's course is
       not the core's) */
    pcx_put_ref(s->x, s->xf, n, r);
    for (i = 0; i < n; i++) if (!isfinite(s->xf[i])) return NEWTON_BAIL_X;
    if (full) {
        bits = (it == 0) ? first_bits : PCX_WANT_ALL;
        Cin = (it == 0 && Cin0) ? Cin0 : C;
    } else {
        /* the chord's `residual_only` (i, q), its companion from the
           held C: the cache's, else the first evaluation's -- what
           `tr._Cmat` holds after the first */
        bits = PCX_PASS_I;
        Cin = Cin0 ? Cin0 : s->C0;
    }
    core(s->core, s->xf, s->T, bits | s->pattern, s->formula, s->a0, s->a1, s->a2, s->h,
         s->theta,
         s->q1, s->q2, s->iq1, s->u, Cin, C, q, iq, Geq, F, J);
    /* the first iterate evaluated, the seed's factors still in place (A is
       written below, the numeric only by the next solve): the early abort's
       check where armed -- the chord on its factors, the full Newton dense */
    if (it == 1 && s->early_lte && (s->chord || !s->klu_sym) && early(s, q, num))
        return NEWTON_EARLY;
    /* `-F_s`: F reduced, negated (no scaler: `_tran_newton_c` declines
       one, and NoneScaler's leaves J and F as they are) */
    pcx_drop_ref(F, s->Fr, n, r);
    for (i = 0; i < m; i++) s->b[i] = -s->Fr[i];
    if (full && s->klu_sym) {
        /* the pattern's values, and |J| on it (zero off it); the full
           Newton's validation reads J reduced itself (`R`) */
        for (k = 0; k < s->knnz; k++) {
            double v = J[s->kat[k]];
            s->kAx[k] = v;
            s->absJ[s->kabs[k]] = fabs(v);
        }
        if (!s->chord && s->pattern) {
            for (k = 0; k < s->knnz; k++) s->R[s->kabs[k]] = s->kAx[k];
        } else if (!s->chord) {
            /* (no pattern: J reduced whole, and `_ours`'s count -- as many
               nonzeros as the pattern's values, `v != 0.0` counting a NaN
               as numpy's count_nonzero does) */
            long nz_a = 0;
            pcx_reduce_rm(J, s->R, n, r);
            for (k = 0; k < s->knnz; k++) nz_a += s->kAx[k] != 0.0;
            if (pcx_count_nonzero(s->R, m * m) != nz_a) return NEWTON_BAIL_KLU_FOREIGN;
        }
    } else if (full) {
        /* remove_row_col of J: A column-major, the copy numpy's solve
           (and SciPy's lu_factor) hands LAPACK, so LAPACK factors the
           same matrix in the same layout; absJ row-major, `abs(J)` as
           the convergence test's dgemv reads it */
        *finite = pcx_reduce_cm_abs(J, s->A, s->absJ, n, r);
    }
    return 0;
}

/* The Newton loop (`StandardNewton` / `ChordNewton.solve_system` on the
   reduced system `_newton` hands them): per iteration the evaluation at
   insert_row(x) (`evaluate`: the core's residual, and J where the
   iteration forms it), the step (`solve_step`), x_next and the limiter
   (`limit_step`), the convergence test (`converged`, a chord's also
   `chord_settled`, then `jacobian_only` there: `fold`), the chord's
   contraction (`chord_breaks`) and its dominant component (`chord_dominant`).
   Contract: THE ITERATION, at the top of this file. */
static long iterate(newton_t *s, void **num)
{
    long n = s->n, m = n - 1, r = s->iref, it, i, st;
    int finite = 1;   /* the last full iteration's J: every entry finite */
    core_fn_t core = (core_fn_t) s->core_fn;
    walk_fn_t walk = (walk_fn_t) s->walk_fn;
    double prev = 0.0;   /* the chord's last step size (`prev`) */
    int have_prev = 0;   /* (`prev is not None`) */
    long dom = -1;       /* the last increment's dominant component (`prev_dx`; none: -1) */
    double dom_d = 0.0;  /* its increment */
    const double *Cin0;  /* cC where the full seed equals cx (a cache hit), else NULL */
    long first_bits;     /* the first evaluation's want bits */
    s->iters = 0;
    s->early_checks = 0;
    s->status = NEWTON_RUNNING;
    /* the seed, reduced (remove_row_col); a hit evaluates without the C
       pass, its C the cached one: as `_tran_core.evaluate` drops its C bit
       on `_C_lookup`'s hit */
    pcx_drop_ref(s->xseed, s->x, n, r);
    Cin0 = cache_hit(s, m, r);
    first_bits = Cin0 ? (PCX_PASS_G | PCX_PASS_I) : PCX_WANT_ALL;
    for (it = 0; it < s->maxiter; it++) {
        /* J evaluated: every iteration of the full Newton, the chord's first */
        int full = (it == 0) || !s->chord;
        if ((st = evaluate(s, it, full, first_bits, Cin0, &finite, *num)) != 0)
            return finish(s, st);
        if ((st = solve_step(s, it, finite, num)) != 0) return finish(s, st);
        /* `x_next = x + xdiff` */
        for (i = 0; i < m; i++) s->xn[i] = s->x[i] + s->b[i];
        if ((st = limit_step(s, walk, n, r)) != 0) return finish(s, st);
        /* a non-finite x_next hands back: the next evaluation would not be
           the core's (and the test below reads every |x| and |x_next|
           finite) */
        for (i = 0; i < m; i++) if (!isfinite(s->xn[i])) return finish(s, NEWTON_BAIL_XNEXT);
        st = converged(s, m);
        /* the Python's count: `i + 1`, the chord's `done` */
        s->iters = it + 1;
        /* (and a chord iterate settled: `_chord_settled`) */
        if (st && s->chord && dom >= 0 && !chord_settled(s, m, dom, dom_d)) st = 0;
        if (st) {
            for (i = 0; i < m; i++) s->x[i] = s->xn[i];
            if (s->fold_j) fold(s, core, n, r);
            return finish(s, NEWTON_CONVERGED);
        }
        if (s->chord && chord_breaks(s, m, &prev, &have_prev))
            return finish(s, NEWTON_BAIL_CHORDBREAK);
        if (s->chord) dom = chord_dominant(s, m, &dom_d);
        /* `x = x_next` */
        for (i = 0; i < m; i++) s->x[i] = s->xn[i];
    }
    return finish(s, NEWTON_BAIL_MAXITER);
}

/* the entry (contract: newton.h); the sparse chord's numeric is the
   solve's own, freed on every return */
long pcx_newton(newton_t *s)
{
    void *num = 0;
    long status = iterate(s, &num);
    if (num) ((klu_free_numeric_t) s->klu_free)(&num, s->klu_common);
    return status;
}
