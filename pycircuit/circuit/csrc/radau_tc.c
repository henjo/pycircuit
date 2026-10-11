/* radau_tc.c -- the Radau IIA(3) transform's simplified Newton in one call:
   `_RadauStages._transform_loop` (_tran_radau.py) against the step's frozen
   factors (`_FrozenTransform.solve`, `_transform_rhs`, `_transform_back`,
   `NumpyLU.solve`, `ComplexKLUSolver.solve_prepared`), bit for bit but for
   the residual check, which decides only outside a band (below).  The glue,
   `_tran_radau_tc.solve`, checks the case, fills radau_tf_t (radau_tc.h),
   writes the factorisations' state back on every return, and on a bail
   rolls back the source memo and runs the Python loop from the same seed. */
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
#include "pcx_abi.h"
#include "pcx_fenv.h"
#include "radau_tc.h"

/* the four flags numpy's errstate governs (inexact is not one).  One
   raised in the C's own arithmetic hands the loop back: the Python makes
   the same operation and numpy reports it as the caller's errstate says.
   The flags the core, the walk, LAPACK and KLU raise are cleared, not
   tested.
   MUST KEEP, at every fetestexcept below: each result it judges is stored
   to a buffer before it.  The compiler does not model the flags (no
   FENV_ACCESS), so an operation whose result stays in a register may move
   past the test; one whose result is stored to memory the call could read
   may not. */

/* every hand-back: the flags cleared, `status` set; returns `why` */
static long tf_bail(radau_tf_t *s, long why)
{
    feclearexcept(FE_ALL_EXCEPT);
    s->status = why;
    return why;
}

/* numpy's complex product (ar + i ai)(br + i bi), in the form numpy's loops
   take on this machine (`cmul_mode`): fused, or two products and a sum.
   The fused form as `cmul_mode` reads it from numpy: re rounds ai bi first
   and fuses ar br into it, im rounds ai br first and fuses ar bi -- the
   other pairing rounds apart */
static void cmul(long fused, double ar, double ai, double br, double bi, double *re, double *im)
{
    if (fused) {
        *re = fma(ar, br, -(ai * bi));
        *im = fma(ar, bi, ai * br);
    } else {
        *re = ar * br - ai * bi;
        *im = ar * bi + ai * br;
    }
}

/* the transform's Newton loop from `seed`.
   Assumes what `_tran_radau_tc.solve` checks and `_Ctx` sizes: n >= 2,
   0 <= iref < n, the buffers of radau_tc.h's sizes; lu_info -1 or 0; the
   KLU handles analysed and factored for the record's pattern; the core
   bound to this circuit and its readiness probed (T); walk_fn NULL or a
   walk ready; A, p, v, fused and bits set.
   Per iteration, as `_transform_loop`: the stages' q and i, the residual,
   `_FrozenTransform.solve` (the right-hand sides, the kept LU's solve, the
   complex factor's refactor -- at most once a record -- and solve with
   `solve_prepared`'s residual check, the update), per stage the step, the
   walk and the limited step's largest reduced entry; the test after the
   update.
   Returns RADAU_TC_CONVERGED: Y the converged stages, `iters`, `status`,
   and where end_bits the end passes (`end_done`, the e_ buffers).  Else a
   bail code, also in `status`: Y partly updated and not read, the flags
   cleared.  On every return lu_a, lu_ipiv and lu_info, the KLU numeric,
   `refactored` and `klu_fresh` hold the factorisations as made so far --
   what the Python's loop would have left: a refactor of the same values
   repeats its bits and dgetrs on dgesv's factors is dgesv's solve, so the
   Python loop run after a bail makes the C's iterates again.  Outside the
   struct's buffers and the factorisations only the core's and the walk's
   own writes are made (the core's bins; the temperature in the elements'
   parameter rows), as the Python's calls of the same C make them. */
long hdl_fn(radau_tf_t *s)
{
    long n = s->n, m = n - 1, r = s->iref, it, i, j, k, e;
    core_fn_t core = (core_fn_t) s->core_fn;
    walk_fn_t walk = (walk_fn_t) s->walk_fn;
    int64_t m64 = m, one64 = 1, info64 = 0;
    const double h = s->h;
    s->iters = 0;
    s->status = RADAU_TC_RUNNING;
    s->refactored = 0;
    /* every input finite: anything else is the Python loop's */
    if (!(isfinite(h) && isfinite(s->reltol) && isfinite(s->abstol))) return tf_bail(s, RADAU_TC_BAIL_NONFINITE);
    for (k = 0; k < 9; k++) if (!isfinite(s->A[k]) || !isfinite(s->v[k])) return tf_bail(s, RADAU_TC_BAIL_NONFINITE);
    for (k = 0; k < 12; k++) if (!isfinite(s->p[k])) return tf_bail(s, RADAU_TC_BAIL_NONFINITE);
    for (k = 0; k < n; k++) if (!isfinite(s->qn[k])) return tf_bail(s, RADAU_TC_BAIL_NONFINITE);
    for (k = 0; k < 3 * n; k++) if (!isfinite(s->u[k])) return tf_bail(s, RADAU_TC_BAIL_NONFINITE);
    memcpy(s->Y, s->seed, (size_t) (3 * n) * sizeof(double));
    for (it = 0; it < s->maxiter; it++) {
        double scale = 0.0, ynorm = 0.0;
        /* the stages' q and i (formula 4: i in the core's bini, q in place;
           a non-finite stage: `passes` declines, the Python evaluates) */
        for (j = 0; j < 3; j++) {
            for (k = 0; k < n; k++) if (!isfinite(s->Y[j * n + k])) return tf_bail(s, RADAU_TC_BAIL_NONFINITE);
            core(s->core, s->Y + j * n, s->T, s->bits, PCX_FORMULA_PASSES, 0.0, 0.0, 0.0, 1.0, 0.0,
                 s->dummy, s->dummy, s->dummy, s->dummy, s->dummy, s->dummy, s->q + j * n,
                 s->dummy, s->dummy, s->dummy, s->dummy);
            memcpy(s->iv + j * n, s->bini, (size_t) n * sizeof(double));
        }
        /* (the core's flags are not numpy's: cleared) */
        feclearexcept(FE_ALL_EXCEPT);
        /* `_coupled_stage_system`'s residual: K_j = -(i_j + u_j); Python's
           `sum` starts at the int 0 (a -0.0 becomes +0.0: `0.0 +`) */
        for (j = 0; j < 3; j++)
            for (k = 0; k < n; k++) s->K[j * n + k] = -(s->iv[j * n + k] + s->u[j * n + k]);
        for (i = 0; i < 3; i++) {
            const double a0 = s->A[3 * i], a1 = s->A[3 * i + 1], a2 = s->A[3 * i + 2];
            double *f = s->F + i * n;
            long rr;
            /* F_i = (q_i - q_n) - h (((0 + A_i0 K_0) + A_i1 K_1) + A_i2 K_2) */
            for (k = 0; k < n; k++) {
                double sm = 0.0 + a0 * s->K[k];
                sm = sm + a1 * s->K[n + k];
                sm = sm + a2 * s->K[2 * n + k];
                f[k] = (s->q[i * n + k] - s->qn[k]) - h * sm;
            }
            /* R_i = red(F_i): row iref dropped */
            for (k = 0, rr = 0; k < n; k++) {
                if (k == r) continue;
                s->R[i * m + rr] = f[k];
                rr++;
            }
        }
        if (fetestexcept(PCX_FP_FLAGS)) return tf_bail(s, RADAU_TC_BAIL_FLAGS);
        /* `_transform_rhs`: -((P_k0 F0 + P_k1 F1) + P_k2 F2), each F cast to
           complex (F + 0j) and multiplied as numpy multiplies; rhs0's
           imaginary part made as numpy makes it (its flags), then dropped.
           The sums left to right and the negation last, as the Python
           expression rounds them */
        for (k = 0; k < m; k++) {
            long kk;
            for (kk = 0; kk < 2; kk++) {
                const double *p = s->p + 6 * kk;
                double t0r, t0i, t1r, t1i, t2r, t2i, sr, si;
                cmul(s->fused, p[0], p[1], s->R[k], 0.0, &t0r, &t0i);
                cmul(s->fused, p[2], p[3], s->R[m + k], 0.0, &t1r, &t1i);
                sr = t0r + t1r;
                si = t0i + t1i;
                cmul(s->fused, p[4], p[5], s->R[2 * m + k], 0.0, &t2r, &t2i);
                sr = sr + t2r;
                si = si + t2i;
                if (kk == 0) {
                    s->b0[k] = -sr;
                    s->x0[k] = -si;
                } else {
                    s->r1[2 * k] = -sr;
                    s->r1[2 * k + 1] = -si;
                }
            }
        }
        if (fetestexcept(PCX_FP_FLAGS)) return tf_bail(s, RADAU_TC_BAIL_FLAGS);
        /* a non-finite real right-hand side is the analysis solver's */
        for (k = 0; k < m; k++) if (!isfinite(s->b0[k])) return tf_bail(s, RADAU_TC_BAIL_RHS);
        /* numpy's kept LU: `dgesv` at its first solve, `dgetrs` after
           (`NumpyLU.solve`; lu_info as radau_tc.h says) */
        memcpy(s->x0, s->b0, (size_t) m * sizeof(double));
        if (s->lu_info < 0) {
            ((dgesv64_t) s->dgesv)(&m64, &one64, s->lu_a, &m64, s->lu_ipiv, s->x0, &m64, &info64);
            s->lu_info = info64;
        } else if (s->lu_info == 0) {
            ((dgetrs64_t) s->dgetrs)("N", &m64, &one64, s->lu_a, &m64, s->lu_ipiv, s->x0, &m64,
                                   &info64, (size_t) 1);
        }
        /* (LAPACK's flags are not numpy's: cleared; a nonzero info is
           numpy's LinAlgError) */
        feclearexcept(FE_ALL_EXCEPT);
        if (s->lu_info != 0) return tf_bail(s, RADAU_TC_BAIL_LU);
        for (k = 0; k < m; k++) if (!isfinite(s->x0[k])) return tf_bail(s, RADAU_TC_BAIL_STEP);
        /* the complex factor: its refactor at most once a record, then the solve */
        if (!s->klu_fresh) {
            if (!((klu_refactor_t) s->klu_refactor)(s->Ap, s->Ai, s->Ax, s->Symbolic, s->Numeric,
                                               s->Common))
                return tf_bail(s, RADAU_TC_BAIL_REFACTOR);
            s->klu_fresh = 1;
            s->refactored++;
        }
        memcpy(s->x1, s->r1, (size_t) (2 * m) * sizeof(double));
        if (!((klu_solve_t) s->klu_solve)(s->Symbolic, s->Numeric, (int) m, 1, s->x1, s->Common))
            return tf_bail(s, RADAU_TC_BAIL_KLUSOLVE);
        /* (KLU's flags are not numpy's: cleared) */
        feclearexcept(FE_ALL_EXCEPT);
        /* the residual check, decided outside the band numpy's `abs` and
           SciPy's product can round apart from these: `solve_prepared`'s
           max|A x - b| / max(max|b|, 1e-300) <= rtol passes only where it
           passes by more than the band; anything else (not finite, a flag,
           failing or too near) hands back and the Python decides */
        {
            /* the arrays in locals, the inner bound read once: `hypot` is a
               libm call that may write memory as far as the compiler knows,
               so a struct field read after one is reloaded -- a local is not
               (the same values, the same operations in the same order) */
            const double *r1 = s->r1, *x1 = s->x1, *Ax = s->Ax;
            const int *Ap = s->Ap, *Ai = s->Ai;
            double *y1 = s->y1, *ab = s->ab;
            double sc = 0.0, res = 0.0, bnd = 0.0, band;
            for (k = 0; k < m; k++) {
                double a = hypot(r1[2 * k], r1[2 * k + 1]);
                if (!(a <= DBL_MAX)) return tf_bail(s, RADAU_TC_BAIL_RESIDUAL);
                if (a > sc) sc = a;
            }
            /* Python's max(sc, 1e-300); sc is finite here */
            if (!(sc >= 1e-300)) sc = 1e-300;
            memset(y1, 0, (size_t) (2 * m) * sizeof(double));
            memset(ab, 0, (size_t) m * sizeof(double));
            /* y1 = A x1 and ab_k = sum_j |A_kj| |x1_j|, column by column */
            for (j = 0; j < m; j++) {
                const double xr = x1[2 * j], xi = x1[2 * j + 1], xa = hypot(xr, xi);
                const long kend = Ap[j + 1];
                for (k = Ap[j]; k < kend; k++) {
                    const long row = Ai[k];
                    const double ar = Ax[2 * k], ai = Ax[2 * k + 1];
                    y1[2 * row] += ar * xr - ai * xi;
                    y1[2 * row + 1] += ar * xi + ai * xr;
                    ab[row] += hypot(ar, ai) * xa;
                }
            }
            for (k = 0; k < m; k++) {
                const double a = hypot(y1[2 * k] - r1[2 * k], y1[2 * k + 1] - r1[2 * k + 1]);
                if (!(a <= DBL_MAX)) return tf_bail(s, RADAU_TC_BAIL_RESIDUAL);
                if (a > res) res = a;
                if (ab[k] > bnd) bnd = ab[k];
            }
            /* an overflow, invalid or divide in the check's own arithmetic
               hands back (underflow is not tested here) */
            if (fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) return tf_bail(s, RADAU_TC_BAIL_RESIDUAL);
            res = res / sc;
            /* the band: 32 (m + 4) eps of max_k sum_j |A_kj| |x_j| over
               max|b|, plus 16 eps of the residual */
            band = (32.0 * (double) (m + 4)) * DBL_EPSILON * (bnd / sc) + 16.0 * DBL_EPSILON * res;
            if (!(res + band <= s->rtol)) return tf_bail(s, RADAU_TC_BAIL_RESIDUAL);
            feclearexcept(FE_ALL_EXCEPT);
        }
        /* `_transform_back`: dY_i = V_i0 w0 + 2 Re(V_i1 w1), the complex
           product whole (its flags), as numpy makes it; then a w0 and 2.0
           times the real part, each rounded, and their sum, a w0 left */
        for (i = 0; i < 3; i++) {
            const double a = s->v[3 * i], br = s->v[3 * i + 1], bi = s->v[3 * i + 2];
            for (k = 0; k < m; k++) {
                double tr, ti, t1, t2;
                cmul(s->fused, br, bi, s->x1[2 * k], s->x1[2 * k + 1], &tr, &ti);
                (void) ti;
                t1 = a * s->x0[k];
                t2 = 2.0 * tr;
                s->dY[i * m + k] = t1 + t2;
            }
        }
        if (fetestexcept(PCX_FP_FLAGS)) return tf_bail(s, RADAU_TC_BAIL_FLAGS);
        /* per stage: Y_prev + insert(dY_i, iref, 0.0), the walk against
           Y_prev, the limited step's largest reduced entry */
        memcpy(s->Yp, s->Y, (size_t) (3 * n) * sizeof(double));
        for (i = 0; i < 3; i++) {
            const double *yp = s->Yp + i * n, *d = s->dY + i * m;
            double *yn = s->Y + i * n, mx = 0.0;
            for (k = 0; k < r; k++) yn[k] = yp[k] + d[k];
            /* the reference row adds the inserted 0.0, not a copy: a -0.0
               there becomes +0.0, as numpy's elementwise sum makes it */
            yn[r] = yp[r] + 0.0;
            for (k = r + 1; k < n; k++) yn[k] = yp[k] + d[k - 1];
            if (fetestexcept(PCX_FP_FLAGS)) return tf_bail(s, RADAU_TC_BAIL_FLAGS);
            /* `cir.limit(Y_trial, Y_prev)`: the walk limits yn in place and
               returns wn, or the element whose limiter declined (the
               Python limits); its flags are not tested */
            if (walk) {
                e = walk(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK, yn, yp, 0, s->wn);
                if (e < s->wn) return tf_bail(s, RADAU_TC_BAIL_WALKSTOP);
                feclearexcept(FE_ALL_EXCEPT);
            }
            /* step_i = red(Y_new - Y_prev) (stored: its flags tested) */
            for (k = 0; k < n; k++) s->stp[k] = yn[k] - yp[k];
            if (fetestexcept(PCX_FP_FLAGS)) return tf_bail(s, RADAU_TC_BAIL_FLAGS);
            /* np.max(np.abs(step_i)), then Python's max(scale, it): each
               takes a later value only where it is greater.  Every value is
               finite and >= +0.0 here, so starting from 0.0 and any tie
               give the same bits */
            for (k = 0; k < n; k++) {
                double vv;
                if (k == r) continue;
                vv = fabs(s->stp[k]);
                if (!isfinite(vv)) return tf_bail(s, RADAU_TC_BAIL_STEP);
                if (vv > mx) mx = vv;
            }
            if (mx > scale) scale = mx;
        }
        /* `_stages_converged`: the step within reltol of the largest
           reduced stage plus abstol.  The reduced entries are finite here
           (each finite before the update, each reduced step finite), so
           the running maximum from 0.0 is numpy's; then two roundings (no
           fma), each
           stored before the test (numpy's scalar arithmetic reports their
           flags) */
        for (i = 0; i < 3; i++)
            for (k = 0; k < n; k++) {
                double vv;
                if (k == r) continue;
                vv = fabs(s->Y[i * n + k]);
                if (vv > ynorm) ynorm = vv;
            }
        feclearexcept(FE_ALL_EXCEPT);
        s->tolv[0] = s->reltol * ynorm;
        s->tolv[1] = s->tolv[0] + s->abstol;
        if (fetestexcept(PCX_FP_FLAGS)) return tf_bail(s, RADAU_TC_BAIL_FLAGS);
        s->iters = it + 1;
        if (scale <= s->tolv[1]) {
            s->status = RADAU_TC_CONVERGED;
            /* THE STEP'S END PASSES: the calls `_stage_end_passes` makes
               next and in its order -- all four passes at the last stage,
               then `C` and `G` at the first two where the shooting reads
               them -- each `passes`'s core call with its bits, outputs and
               dummies, into the caller's buffers (the core's arrays its
               scratch, each call writing what it reads).  A non-finite
               stage: none made, `end_done` 0 (`passes` declines such a
               state; the Python's calls then run).  Either way the result
               is RADAU_TC_CONVERGED */
            s->end_done = 0;
            if (s->end_bits) {
                for (k = 0; k < 3 * n; k++) if (!isfinite(s->Y[k])) return RADAU_TC_CONVERGED;
                core(s->core, s->Y + 2 * n, s->T, s->end_bits, PCX_FORMULA_PASSES, 0.0, 0.0, 0.0, 1.0, 0.0,
                     s->dummy, s->dummy, s->dummy, s->dummy, s->e_C, s->e_C, s->e_q,
                     s->dummy, s->dummy, s->dummy, s->dummy);
                memcpy(s->e_i, s->bini, (size_t) n * sizeof(double));
                memcpy(s->e_G, s->binG, (size_t) (n * n) * sizeof(double));
                if (s->end_cg) {
                    for (j = 0; j < 2; j++) {
                        core(s->core, s->Y + j * n, s->T, s->end_cg, PCX_FORMULA_PASSES, 0.0, 0.0, 0.0, 1.0, 0.0,
                             s->dummy, s->dummy, s->dummy, s->dummy, s->e_C01 + j * n * n,
                             s->e_C01 + j * n * n, s->e_q01, s->dummy, s->dummy, s->dummy,
                             s->dummy);
                        memcpy(s->e_G01 + j * n * n, s->binG, (size_t) (n * n) * sizeof(double));
                    }
                }
                /* (the core's flags are not numpy's: cleared) */
                feclearexcept(FE_ALL_EXCEPT);
                s->end_done = 1;
            }
            return RADAU_TC_CONVERGED;
        }
    }
    return tf_bail(s, RADAU_TC_BAIL_MAXITER);
}
