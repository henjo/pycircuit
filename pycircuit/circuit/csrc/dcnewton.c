/* dcnewton.c -- the DC operating point's base Newton solve in C (contract:
   dcnewton.h), `StandardNewton.solve_system` on the evaluation the DC's
   ladders compose, bit for bit, or handed back.  Its glue is `_dc_newton_c`.

   It tests floating-point flags under THE FLAG RULE (pcx_fenv.h, MUST
   KEEP): each operation a test covers stores its result through a buffer of
   the struct, or decides a branch taken before the test.  What the device
   kernels, KLU and the walk raise inside is cleared after them: the
   Python's calls of them (cffi, ctypes) never report it.  A flag numpy
   would report -- the residual's and the layers' arithmetic, the step,
   `arr @ x`, the convergence test's dot -- hands the solve back, and the
   Python's solve gives its warnings.  Every return leaves the flags clear. */
#include <stdint.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <fenv.h>
#include "pcx_abi.h"
#include "pcx_minmax.h"
#include "pcx_fenv.h"
#include "pcx_dense.h"
#include "core.h"
#include "pcx_klu.h"
#include "pcx_verdict.h"
#include "pcx_kluseq.h"
#include "dcnewton.h"

/* the flags numpy reports here: underflow is not one -- numpy ignores it by
   default, and the glue declines where it does not (`_dc_newton_c`) */
#define DC_FLAGS (FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

/* every return: the numeric as the outcome leaves it (`pcx_klu_finish`: a
   hand-back the entry's, CONVERGED and MAXITER the sequence's) -- the flags
   clear */
static long finish(dcnewton_t *s, long status)
{
    if (s->klu.klu_sym) pcx_klu_finish(&s->klu, status >= DCNEWTON_BAIL_FLAGS);
    feclearexcept(FE_ALL_EXCEPT);
    return s->status = status;
}

/* `x` (reduced) evaluated as the composition: the core's G and i at
   insert_row(x); F = cut(i + u) or cut(i + lam * u), J = cut(G) on the
   active entries; gshunt's `F + x g`, `J + 0.0` and g on the diagonal;
   junction gmin's terms, junction by junction.  0, or DCNEWTON_BAIL_FLAGS
   (a flag, or a value that is not finite) */
static long evaluate(dcnewton_t *s)
{
    const hdl_core_t *c = (const hdl_core_t *) s->core;
    long n = s->n, m = n - 1, r = s->iref, i, k;
    pcx_put_ref(s->x, s->xf, n, r);
    ((core_fn_t) s->core_fn)(s->core, s->xf, s->T, PCX_PASS_G | PCX_PASS_I | s->pattern,
                             PCX_FORMULA_PASSES, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL, NULL,
                             NULL, NULL, s->Cs, s->q, NULL, NULL, NULL, NULL);
    /* (the kernels' flags: the Python's cffi calls never report them) */
    feclearexcept(FE_ALL_EXCEPT);
    if (s->has_lam)
        for (i = 0; i < n; i++) s->Ff[i] = c->bini[i] + s->lam * s->u[i];
    else
        for (i = 0; i < n; i++) s->Ff[i] = c->bini[i] + s->u[i];
    pcx_drop_ref(s->Ff, s->Fr, n, r);
    for (k = 0; k < s->nact; k++) s->R[s->act_r[k]] = c->binG[s->act_f[k]];
    if (s->has_gshunt) {
        double g = s->gshunt;
        for (i = 0; i < m; i++) s->Fr[i] = s->Fr[i] + s->x[i] * g;
        for (k = 0; k < s->nact; k++) s->R[s->act_r[k]] = s->R[s->act_r[k]] + 0.0;
        for (i = 0; i < m; i++) s->R[i * m + i] += g;
    }
    if (s->has_jgmin) {
        double g = s->gjgmin;
        for (k = 0; k < s->nj; k++) {
            long a = s->ja[k], b = s->jb[k];
            double vj = s->x[a] - s->x[b];
            s->Fr[a] += g * vj;
            s->Fr[b] -= g * vj;
            s->R[a * m + a] += g;
            s->R[a * m + b] -= g;
            s->R[b * m + a] -= g;
            s->R[b * m + b] += g;
        }
    }
    if (fetestexcept(DC_FLAGS)) return DCNEWTON_BAIL_FLAGS;
    for (i = 0; i < m; i++) if (!isfinite(s->Fr[i])) return DCNEWTON_BAIL_FLAGS;
    for (k = 0; k < s->nact; k++) if (!isfinite(s->R[s->act_r[k]])) return DCNEWTON_BAIL_FLAGS;
    return 0;
}

/* StructuralKLU.solve(J, -F) on its own sequence (pcx_kluseq.h): J's active
   entries gathered in its order; the step in `kx`.  0, or why it is handed
   back (a failed factor or solve is the Python's raise) */
static long klu_step(dcnewton_t *s, long m)
{
    long i, st;
    for (i = 0; i < s->nact; i++) s->klu.kAx[i] = s->R[s->act_r[i]];
    st = pcx_klu_step(&s->klu, s->R, s->b, s->kx, s->ky, NULL, m, s->dgemv, 1, NULL);
    if (st == PCX_KLU_FLAGS) return DCNEWTON_BAIL_FLAGS;
    return st ? DCNEWTON_BAIL_SINGULAR : 0;
}

/* numpy.linalg.solve(J, -F): dgesv on J's column-major copy -- an info, or
   the invalid flag numpy's linalg turns into LinAlgError, is the Python's
   raise; what else it raises numpy's linalg ignores.  The step in `b` */
static long dense_step(dcnewton_t *s, long m)
{
    int64_t m64 = m, one64 = 1, info64 = 0;
    long i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < m; j++) s->A[j * m + i] = s->R[i * m + j];
    ((dgesv64_t) s->dgesv)(&m64, &one64, s->A, &m64, s->ipiv64, s->b, &m64, &info64);
    if (info64 != 0 || fetestexcept(FE_INVALID)) return DCNEWTON_BAIL_SINGULAR;
    feclearexcept(FE_ALL_EXCEPT);
    return 0;
}

/* `x_next = x + xdiff`, the DC's `limiter_func` -- insert_row of x_next
   and of x, `cir.limit` (the walk, x_next in place), the reference dropped
   -- and `xdiff = x_next - x`, always (the DC's limiter is always given).
   0, or why it is handed back */
static long step(dcnewton_t *s, const double *dx)
{
    long n = s->n, m = n - 1, r = s->iref, i;
    for (i = 0; i < m; i++) s->xn[i] = s->x[i] + dx[i];
    if (fetestexcept(DC_FLAGS)) return DCNEWTON_BAIL_FLAGS;
    for (i = 0; i < m; i++) if (!isfinite(s->xn[i])) return DCNEWTON_BAIL_FLAGS;
    if (s->walk_fn) {
        long e;
        pcx_put_ref(s->xn, s->xnf, n, r);
        pcx_put_ref(s->x, s->x0f, n, r);
        e = ((walk_fn_t) s->walk_fn)(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK,
                                     s->xnf, s->x0f, 0, s->wn);
        if (e < s->wn) return DCNEWTON_BAIL_WALKSTOP;
        /* (the limiters' flags: the walk's cffi call never reports them) */
        feclearexcept(FE_ALL_EXCEPT);
        pcx_drop_ref(s->xnf, s->xn, n, r);
    }
    for (i = 0; i < m; i++) s->xd[i] = s->xn[i] - s->x[i];
    if (fetestexcept(DC_FLAGS)) return DCNEWTON_BAIL_FLAGS;
    for (i = 0; i < m; i++) if (!isfinite(s->xn[i])) return DCNEWTON_BAIL_FLAGS;
    return 0;
}

/* nrsolver's convergence test, as newton.c's: `I_scale = dot(abs(J),
   ax_next) + aF` (numpy's own dgemv on the dense |J|: BLAS's summation
   order), `abs(xdiff) < reltol * maximum(ax_next, abs(x)) + xtol`, `aF <
   reltol * I_scale + abstol`.  1 converged, 0 not, -1 a flag (np.dot and
   the arithmetic report it: handed back); `*held` the residual's test */
static int converged(dcnewton_t *s, long m, int *held)
{
    long i, k;
    int conv_x = 1, conv_f = 1;
    for (k = 0; k < s->nact; k++) s->absJ[s->act_r[k]] = fabs(s->R[s->act_r[k]]);
    for (i = 0; i < m; i++) s->ax[i] = fabs(s->xn[i]);
    ((dgemv64_t) s->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, m, m, 1.0, s->absJ, m,
                           s->ax, 1, 0.0, s->Is, 1);
    for (i = 0; i < m; i++) {
        double aF = fabs(s->Fr[i]);
        double Is = s->Is[i] + aF;
        double mx = pcx_np_maximum(s->ax[i], fabs(s->x[i]));
        if (!(fabs(s->xd[i]) < s->reltol * mx + s->xtol[i])) conv_x = 0;
        if (!(aF < s->reltol * Is + s->abstol[i])) conv_f = 0;
    }
    if (fetestexcept(DC_FLAGS)) return -1;
    *held = conv_f;
    return conv_x && conv_f;
}

/* the update's largest entry, `float(np.max(np.abs(xdiff)))`: every |xd|
   finite here (`step`), where a scan from -INFINITY is numpy's maximum */
static double step_size(const dcnewton_t *s, long m)
{
    long i;
    double size = -INFINITY;
    for (i = 0; i < m; i++) if (fabs(s->xd[i]) > size) size = fabs(s->xd[i]);
    return size;
}

/* THE ROUNDING FLOOR (`nrsolver._within_rounding`, as newton.c's
   `floor_holds`), on an iteration whose residual's test held and whose
   update failed and stopped contracting: `I_scale` (Isf) the test's
   `Is + |F|`; the unknowns failing the update (`np.flatnonzero`'s order),
   `floor_rows` at a time: their rows of J^-1 -- `np.linalg.solve(J.T, E)`
   on the dense J (dgesv on R's row-major copy, J^T's column-major one,
   made afresh for each block as numpy's solve makes it; what it raises
   with info 0 numpy's linalg clears), klu_tsolve on the numeric this
   iteration's sequence holds (`klu.klu_out`) on KLU's (its flags, as KLU's
   every call's, never the Python's) -- each floor `eps * np.dot(|row|,
   I_scale)` (0.0 plus ddot) and the tolerance plus it, a block's all made
   before the flags are read (THE FLAG RULE: kept in Ef), then the update
   within them.  1 where every failing unknown is; 0 at the first block
   with one that is not, and where a solve fails (numpy's LinAlgError;
   KLU's); -1 a flag */
static int floor_holds(dcnewton_t *s, long m)
{
    long i, j, c, k = 0, lo, kk;
    for (i = 0; i < m; i++) s->Isf[i] = s->Is[i] + fabs(s->Fr[i]);
    for (i = 0; i < m; i++) {
        double mx = pcx_np_maximum(s->ax[i], fabs(s->x[i]));
        if (!(fabs(s->xd[i]) < s->reltol * mx + s->xtol[i])) s->ffail[k++] = i;
    }
    if (k == 0) return 0;
    for (lo = 0; lo < k; lo += s->floor_rows) {
        kk = k - lo < s->floor_rows ? k - lo : s->floor_rows;
        for (j = 0; j < kk * m; j++) s->Ef[j] = 0.0;
        for (c = 0; c < kk; c++) s->Ef[c * m + s->ffail[lo + c]] = 1.0;
        if (s->klu.klu_sym) {
            int ok = ((klu_solve_t) s->klu_tsolve)(s->klu.klu_sym, s->klu.klu_out, (int) m,
                                                   (int) kk, s->Ef, s->klu.klu_common);
            feclearexcept(FE_ALL_EXCEPT);
            if (!ok) return 0;
        } else {
            int64_t m64 = m, k64 = kk, info64 = 0;
            for (j = 0; j < m * m; j++) s->Tf[j] = s->R[j];
            ((dgesv64_t) s->dgesv)(&m64, &k64, s->Tf, &m64, s->ipivf, s->Ef, &m64, &info64);
            feclearexcept(FE_ALL_EXCEPT);
            if (info64 != 0) return 0;
        }
        for (c = 0; c < kk; c++) {
            double mx;
            i = s->ffail[lo + c];
            for (j = 0; j < m; j++) s->af[j] = fabs(s->Ef[c * m + j]);
            mx = pcx_np_maximum(s->ax[i], fabs(s->x[i]));
            s->Ef[c * m] = s->reltol * mx + s->xtol[i]
                           + DBL_EPSILON * (0.0 + ((ddot64_t) s->ddot)(m, s->af, 1, s->Isf, 1));
        }
        if (fetestexcept(DC_FLAGS)) return -1;
        for (c = 0; c < kk; c++)
            if (!(fabs(s->xd[s->ffail[lo + c]]) < s->Ef[c * m])) return 0;
    }
    return 1;
}

/* the Newton loop (StandardNewton.solve_system): evaluate, solve, step,
   test (the rounding floor's too); `x = x_next` after each.  CONVERGED (x the solution), MAXITER
   (the last iteration's Fr, xd, xn, and x = x_next as the Python's loop
   leaves it), or a hand-back */
static long iterate(dcnewton_t *s)
{
    long m = s->n - 1, it, i, st;
    int conv, held = 0;
    double size_prev = 0.0;  /* the last update's largest entry (`size_prev`) */
    int have_size = 0;       /* (`size_prev is not None`) */
    for (it = 0; it < s->maxiter; it++) {
        if ((st = evaluate(s)) != 0) return st;
        for (i = 0; i < m; i++) s->b[i] = -s->Fr[i];
        st = s->klu.klu_sym ? klu_step(s, m) : dense_step(s, m);
        if (st != 0) return st;
        if ((st = step(s, s->klu.klu_sym ? s->kx : s->b)) != 0) return st;
        if ((conv = converged(s, m, &held)) < 0) return DCNEWTON_BAIL_FLAGS;
        s->iters = it + 1;
        /* (or the update stalled -- not under half the last one's largest
           entry -- with the residual's test held: the rounding floor) */
        if (!conv) {
            double size = step_size(s, m);
            if (held && have_size && !(size < 0.5 * size_prev)
                    && (conv = floor_holds(s, m)) < 0)
                return DCNEWTON_BAIL_FLAGS;
            size_prev = size;
            have_size = 1;
        }
        for (i = 0; i < m; i++) s->x[i] = s->xn[i];
        if (conv) return DCNEWTON_CONVERGED;
    }
    return DCNEWTON_MAXITER;
}

/* the entry (contract: dcnewton.h) */
long pcx_dcnewton(dcnewton_t *s)
{
    s->iters = 0;
    s->klu.refactors = s->klu.factors = s->klu.pivot_fallbacks = 0;
    s->klu.residual_fallbacks = 0;
    s->klu.klu_out = s->klu.klu_num;
    feclearexcept(FE_ALL_EXCEPT);
    return finish(s, iterate(s));
}
