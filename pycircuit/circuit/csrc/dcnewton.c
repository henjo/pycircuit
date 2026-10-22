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
   the arithmetic report it: handed back) */
static int converged(dcnewton_t *s, long m)
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
    return conv_x && conv_f;
}

/* the Newton loop (StandardNewton.solve_system): evaluate, solve, step,
   test; `x = x_next` after each.  CONVERGED (x the solution), MAXITER
   (the last iteration's Fr, xd, xn, and x = x_next as the Python's loop
   leaves it), or a hand-back */
static long iterate(dcnewton_t *s)
{
    long m = s->n - 1, it, i, st;
    int conv;
    for (it = 0; it < s->maxiter; it++) {
        if ((st = evaluate(s)) != 0) return st;
        for (i = 0; i < m; i++) s->b[i] = -s->Fr[i];
        st = s->klu.klu_sym ? klu_step(s, m) : dense_step(s, m);
        if (st != 0) return st;
        if ((st = step(s, s->klu.klu_sym ? s->kx : s->b)) != 0) return st;
        if ((conv = converged(s, m)) < 0) return DCNEWTON_BAIL_FLAGS;
        s->iters = it + 1;
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
