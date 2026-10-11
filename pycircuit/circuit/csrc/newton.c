
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include "pcx_abi.h"
#include "pcx_dense.h"
#include "newton.h"

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
   sparse chord's numeric where it made one (the caller frees it) */
static long iterate(newton_t *s, void **num)
{
    long n = s->n, m = n - 1, r = s->iref, it, i, k;
    int finite = 1;   /* the last full iteration's J: every entry finite */
    core_fn_t core = (core_fn_t) s->core_fn;
    walk_fn_t walk = (walk_fn_t) s->walk_fn;
    dgemv64_t gemv = (dgemv64_t) s->dgemv;
    double prev = 0.0;   /* the chord's last step size (`prev`) */
    int have_prev = 0;   /* (`prev is not None`) */
    int64_t m64 = m, one64 = 1, info64 = 0;
    int m32 = (int) m, one32 = 1, info32 = 0;
    char trans = 'N';
    s->iters = 0;
    s->status = NEWTON_RUNNING;
    s->last_eval = -1;
    s->hit = 0;
    /* the seed, reduced (remove_row_col), and _C_lookup's equality test of
       the full seed it evaluates against the cached state: elementwise ==
       (`alltrue(x_cached == x)`, a NaN never equal), the full seed being
       insert_row(x) -- +0.0 at the reference, so cx[r] is compared with
       0.0 */
    pcx_drop_ref(s->xseed, s->x, n, r);
    s->Cin0 = 0;
    if (s->have_cache) {
        int eq = 1;
        for (i = 0; i < r; i++) if (!(s->cx[i] == s->x[i])) { eq = 0; break; }
        if (eq && !(s->cx[r] == 0.0)) eq = 0;
        if (eq) for (i = r; i < m; i++) if (!(s->cx[i + 1] == s->x[i])) { eq = 0; break; }
        if (eq) { s->Cin0 = s->cC; s->hit = 1; }
    }
    /* a hit evaluates without the C pass, its C the cached one: as
       `_tran_core.evaluate` drops its C bit on `_C_lookup`'s hit */
    s->first_bits = s->Cin0 ? (PCX_PASS_G | PCX_PASS_I) : PCX_WANT_ALL;
    for (it = 0; it < s->maxiter; it++) {
        const double *Cin;
        long bits;
        double *F, *J, *C, *q, *iq, *Geq;
        /* J evaluated: every iteration of the full Newton, the chord's first */
        int full = (it == 0) || !s->chord;
        if (it == 0) { F = s->F0; J = s->J0; C = s->C0; q = s->q0; iq = s->iq0; Geq = s->Geq0; }
        else { F = s->F; J = s->J; C = s->C; q = s->q; iq = s->iq; Geq = s->Geq; }
        /* insert_row: the full x of this evaluation; a non-finite one hands
           back (`_tran_core.evaluate` declines one: the Python's course is
           not the core's) */
        pcx_put_ref(s->x, s->xf, n, r);
        for (i = 0; i < n; i++) if (!isfinite(s->xf[i])) { s->status = NEWTON_BAIL_X; return NEWTON_BAIL_X; }
        if (full) {
            bits = (it == 0) ? s->first_bits : PCX_WANT_ALL;
            Cin = (it == 0 && s->Cin0) ? s->Cin0 : C;
        } else {
            /* the chord's `residual_only` (i, q), its companion from the
               held C: the cache's, else the first evaluation's -- what
               `tr._Cmat` holds after the first */
            bits = PCX_PASS_I;
            Cin = s->Cin0 ? s->Cin0 : s->C0;
        }
        core(s->core, s->xf, s->T, bits | s->pattern, s->formula, s->a0, s->a1, s->a2, s->h,
             s->theta,
             s->q1, s->q2, s->iq1, s->u, Cin, C, q, iq, Geq, F, J);
        s->last_eval = it;
        /* `-F_s`: F reduced, negated (no scaler: `_tran_newton_c` declines
           one, and NoneScaler's leaves J and F as they are) */
        pcx_drop_ref(F, s->Fr, n, r);
        for (i = 0; i < m; i++) s->b[i] = -s->Fr[i];
        if (full && s->klu_sym) {
            /* the pattern's values, and |J| on it (zero off it) */
            for (k = 0; k < s->knnz; k++) {
                double v = J[s->kat[k]];
                s->kAx[k] = v;
                s->absJ[s->kabs[k]] = fabs(v);
            }
        } else if (full) {
            /* remove_row_col of J: A column-major, the copy numpy's solve
               (and SciPy's lu_factor) hands LAPACK, so LAPACK factors the
               same matrix in the same layout; absJ row-major, `abs(J)` as
               the convergence test's dgemv reads it */
            finite = pcx_reduce_cm_abs(J, s->A, s->absJ, n, r);
        }
        if (s->chord && s->klu_sym) {
            /* StructuralKLU.factor's numeric, then `_KLUNumeric.solve`; a
               non-finite matrix or right-hand side hands back */
            if (it == 0) {
                for (k = 0; k < s->knnz; k++) if (!isfinite(s->kAx[k])) { s->status = NEWTON_BAIL_JNONFINITE; return NEWTON_BAIL_JNONFINITE; }
                *num = ((klu_factor_t) s->klu_factor)(s->kAp, s->kAi, s->kAx, s->klu_sym,
                                                      s->klu_common);
                if (!*num) { s->status = NEWTON_BAIL_KLU_FACTOR; return NEWTON_BAIL_KLU_FACTOR; }
            }
            for (i = 0; i < m; i++) if (!isfinite(s->b[i])) { s->status = NEWTON_BAIL_BNONFINITE; return NEWTON_BAIL_BNONFINITE; }
            if (!((klu_solve_t) s->klu_solve)(s->klu_sym, *num, m32, 1, s->b, s->klu_common)) {
                s->status = NEWTON_BAIL_KLU_SOLVE;
                return NEWTON_BAIL_KLU_SOLVE;
            }
        } else if (s->chord) {
            /* the dense chord: SciPy's lu_factor (dgetrf) once, lu_solve
               (dgetrs, 'N') every iteration; both refuse a non-finite
               matrix or right-hand side (the Python chord then breaks to
               its fallback) -- handed back, as is any LAPACK info */
            if (it == 0) {
                if (!finite) { s->status = NEWTON_BAIL_JNONFINITE; return NEWTON_BAIL_JNONFINITE; }
                ((dgetrf32_t) s->dgetrf)(&m32, &m32, s->A, &m32, s->ipiv32, &info32);
                if (info32 != 0) { s->status = NEWTON_BAIL_GETRF; return NEWTON_BAIL_GETRF; }
            }
            for (i = 0; i < m; i++) if (!isfinite(s->b[i])) { s->status = NEWTON_BAIL_BNONFINITE; return NEWTON_BAIL_BNONFINITE; }
            ((dgetrs32_t) s->dgetrs)(&trans, &m32, &one32, s->A, &m32, s->ipiv32, s->b, &m32,
                                   &info32, 1);
            if (info32 != 0) { s->status = NEWTON_BAIL_GETRS; return NEWTON_BAIL_GETRS; }
        } else {
            /* numpy.linalg.solve: numpy's own dgesv on the column-major
               copy; info != 0 (singular: numpy's LinAlgError) hands back */
            ((dgesv64_t) s->dgesv)(&m64, &one64, s->A, &m64, s->ipiv64, s->b, &m64, &info64);
            if (info64 != 0) { s->status = NEWTON_BAIL_DGESV; return NEWTON_BAIL_DGESV; }
        }
        /* `x_next = x + xdiff` */
        for (i = 0; i < m; i++) s->xn[i] = s->x[i] + s->b[i];
        if (walk) {
            /* `_newton_limiter`'s closure: `cir.limit` on insert_row of
               x_next and of x, the reference dropped after; then
               `xdiff = x_next - x` */
            long e;
            pcx_put_ref(s->xn, s->xnf, n, r);
            pcx_put_ref(s->x, s->x0f, n, r);
            e = walk(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK, s->xnf, s->x0f, 0, s->wn);
            if (e < s->wn) { s->status = NEWTON_BAIL_WALKSTOP; return NEWTON_BAIL_WALKSTOP; }
            pcx_drop_ref(s->xnf, s->xn, n, r);
            for (i = 0; i < m; i++) s->xd[i] = s->xn[i] - s->x[i];
        } else {
            for (i = 0; i < m; i++) s->xd[i] = s->b[i];
        }
        /* a non-finite x_next hands back: the next evaluation would not be
           the core's (and below, every |x| and |x_next| is finite) */
        for (i = 0; i < m; i++) if (!isfinite(s->xn[i])) { s->status = NEWTON_BAIL_XNEXT; return NEWTON_BAIL_XNEXT; }
        /* `toolkit.dot(abs(J), ax_next)`: numpy's own dgemv, row-major, as
           np.dot calls it on a C-ordered matrix */
        for (i = 0; i < m; i++) s->ax[i] = fabs(s->xn[i]);
        gemv(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, m, m, 1.0, s->absJ, m, s->ax, 1, 0.0, s->Is, 1);
        {
            /* nrsolver's test, entry by entry, in its operations and order:
               `I_scale = dot(...) + aF`, `abs(xdiff) < reltol *
               maximum(ax_next, abs(x)) + xtol`, `aF < reltol * I_scale +
               abstol`; written `!(a < b)` so a NaN fails, as numpy's `<` */
            int conv_x = 1, conv_f = 1;
            for (i = 0; i < m; i++) {
                double aF = fabs(s->Fr[i]);
                double Is = s->Is[i] + aF;
                double axi = s->ax[i], axo = fabs(s->x[i]);
                /* np.maximum(ax_next, abs(x)): both finite here (x checked
                   at its evaluation, x_next above), where this ternary is
                   the same value; it would not propagate a NaN as
                   np.maximum does */
                double mx = (axi >= axo) ? axi : axo;
                if (!(fabs(s->xd[i]) < s->reltol * mx + s->xtol[i])) conv_x = 0;
                if (!(aF < s->reltol * Is + s->abstol[i])) conv_f = 0;
            }
            /* the Python's count: `i + 1`, the chord's `done` */
            s->iters = it + 1;
            if (conv_x && conv_f) {
                for (i = 0; i < m; i++) s->x[i] = s->xn[i];
                if (s->fold_j) {
                    /* jacobian_only at the converged point: G, C, q -- C and
                       q alone where J is unread (fold_j holds the passes).
                       Cj is both Cin and C: fold_j always has the C pass, so
                       the companion reads the C just made; F is not written
                       (no PCX_PASS_I) */
                    pcx_put_ref(s->x, s->xc, n, r);
                    core(s->core, s->xc, s->T, s->fold_j | s->pattern, s->formula, s->a0, s->a1,
                         s->a2, s->h,
                         s->theta, s->q1, s->q2, s->iq1, s->u, s->Cj, s->Cj, s->qj, s->iqj,
                         s->Geqj, s->F, s->Jj);
                }
                s->status = NEWTON_CONVERGED;
                return NEWTON_CONVERGED;
            }
        }
        if (s->chord) {
            /* the chord's contraction test: `size =
               float(np.max(np.abs(xdiff)))` -- np.max propagates a NaN,
               which a plain `>` scan would skip; from -INFINITY, below any
               |v| -- then `not np.isfinite(size) or not size < prev`, as
               the Python spells it */
            double size = -INFINITY;
            for (i = 0; i < m; i++) {
                double v = fabs(s->xd[i]);
                if (isnan(v)) { size = v; break; }
                if (v > size) size = v;
            }
            if (!isfinite(size) || (have_prev && !(size < prev))) { s->status = NEWTON_BAIL_CHORDBREAK; return NEWTON_BAIL_CHORDBREAK; }
            prev = size;
            have_prev = 1;
        }
        /* `x = x_next` */
        for (i = 0; i < m; i++) s->x[i] = s->xn[i];
    }
    s->status = NEWTON_BAIL_MAXITER;
    return NEWTON_BAIL_MAXITER;
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
