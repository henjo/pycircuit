/* radau.c -- the Radau IIA(3) coupled stage Newton in one call: `_stage_newton`
   (`_RadauStages._coupled_stage_solver`, _tran_radau.py) undamped, without a
   shunt or a stateful limiter, bit for bit -- its assembly
   (`_coupled_stage_system`), the solve -- numpy's, or StructuralKLU's own
   sequence where the transient solves the stage block on one -- the update
   and the limiter walk, and the test (`_stages_converged`).  The glue,
   `_tran_radau_c.solve`,
   checks the case, fills radau_t (radau.h), and on RADAU_CONVERGED records
   every assembly in the device memo; on any other status it rolls back the
   source memo and the Python Newton runs from the same seed. */
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#include "pcx_abi.h"
#include "pcx_minmax.h"
#include "pcx_fenv.h"
#include "pcx_klu.h"
#include "pcx_verdict.h"
#include "pcx_kluseq.h"
#include "radau.h"

/* the flags numpy reports here: underflow is not one -- numpy ignores it by
   default, and the glue declines where it does not (`_tran_radau_c`).  One
   raised in the C's own arithmetic hands the solve back: the Python makes
   the same operation and numpy reports it as the caller's errstate says.
   The flags the core, the walk and LAPACK raise are cleared, not tested:
   the Python calls the same core and walk C, and numpy's solve clears
   LAPACK's.
   MUST KEEP, at every fetestexcept below: each result it judges is stored
   to a buffer before it.  The compiler does not model the flags (no
   FENV_ACCESS), so an operation whose result stays in a register may move
   past the test; one whose result is stored to memory the call could read
   may not. */
#define RC_FLAGS (FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

/* every hand-back: the KLU's numeric the entry's (THE NUMERIC, pcx_klu.h),
   the flags cleared, `status` set; returns `why` */
static long rc_bail(radau_t *s, long why)
{
    if (s->klu.klu_sym) pcx_klu_finish(&s->klu, 1);
    feclearexcept(FE_ALL_EXCEPT);
    s->status = why;
    return why;
}

/* one assembly at the stages Ya[a]: each stage's passes through the core
   (formula 4: G in binG, i in bini, C and q in place), then the residual and
   the blocks in `_coupled_stage_system`'s operations over the full width
   (each result stored, so each operation precedes the flag test).
   Assumes: slot a exists (a < maxa) and Ya[a] holds three stages; u, qn, A,
   h set.  Writes slot a's qa, ia, Ca, Ga, and K, F, R, and the blocks: Jf
   (unfactored) and blk, or with the KLU bound the pattern's entries alone
   in klu.kAx -- the same operations on each, and the block zero off it.
   Returns 0; RADAU_BAIL_NONFINITE where a stage is not finite (checked
   before that stage's core call, the earlier stages' passes already made);
   RADAU_BAIL_FLAGS where the residual or the blocks raised a flag. */
static long assemble(radau_t *s, long a)
{
    long n = s->n, m = n - 1, r = s->iref, M = 3 * m, nn = n * n;
    long i, j, k, rr, cc, ri, ci;
    core_fn_t core = (core_fn_t) s->core_fn;
    double *Y = s->Ya + a * 3 * n, *q = s->qa + a * 3 * n, *iv = s->ia + a * 3 * n;
    double *C = s->Ca + a * 3 * nn, *G = s->Ga + a * 3 * nn;
    double h = s->h;
    for (j = 0; j < 3; j++) {
        /* (a non-finite stage: `passes` declines and the Python evaluates) */
        for (k = 0; k < n; k++) if (!isfinite(Y[j * n + k])) return RADAU_BAIL_NONFINITE;
        core(s->core, Y + j * n, s->T, PCX_WANT_ALL, PCX_FORMULA_PASSES, 0.0, 0.0, 0.0, 1.0, 0.0,
             s->dummy, s->dummy, s->dummy, s->dummy, s->dummy, C + j * nn, q + j * n,
             s->dummy, s->dummy, s->dummy, s->dummy);
        memcpy(G + j * nn, s->binG, (size_t) nn * sizeof(double));
        memcpy(iv + j * n, s->bini, (size_t) n * sizeof(double));
    }
    /* (the core's flags are not numpy's: cleared) */
    feclearexcept(FE_ALL_EXCEPT);
    /* K_j = -(i_j + u_j): the sum, then its negation, as the Python's
       `-(i_j + src(t_j))` */
    for (j = 0; j < 3; j++)
        for (k = 0; k < n; k++) s->K[j * n + k] = -(iv[j * n + k] + s->u[j * n + k]);
    for (i = 0; i < 3; i++) {
        const double a0 = s->A[3 * i], a1 = s->A[3 * i + 1], a2 = s->A[3 * i + 2];
        double *f = s->F + i * n;
        for (k = 0; k < n; k++) {
            /* F_i = (q_i - q_n) - h sum(A_ij K_j): Python's `sum` starts at
               the int 0, so its first addition is 0 + A_i0 K_0 -- a -0.0
               product becomes +0.0 (`0.0 +`, not the product alone); then
               left to right, as `sum` adds */
            double sm = 0.0 + a0 * s->K[k];
            sm = sm + a1 * s->K[n + k];
            sm = sm + a2 * s->K[2 * n + k];
            f[k] = (q[i * n + k] - s->qn[k]) - h * sm;
        }
        /* R_i = red(F_i): row iref dropped */
        for (k = 0, rr = 0; k < n; k++) {
            if (k == r) continue;
            s->R[i * m + rr] = f[k];
            rr++;
        }
    }
    if (s->klu.klu_sym) {
        /* the pattern's entries, in its order: entry p of block i j at the
           full C_i, G_j place act_f[p], as below; off the pattern each
           block is an exact zero -- a product with an exact zero, raising
           no flag */
        double hA[9];
        long p;
        for (k = 0; k < 9; k++) hA[k] = h * s->A[k];
        for (p = 0; p < s->nact; p++) {
            long ij = s->act_b[p], f = s->act_f[p], bi = ij / 3, bj = ij - 3 * bi;
            double v = hA[ij] * G[bj * nn + f];
            s->klu.kAx[p] = bi == bj ? C[bi * nn + f] + v : v;
        }
        /* (the reference row's and column's, for their flags) */
        for (p = 0; p < s->nref; p++) {
            long ij = s->ref_b[p], f = s->ref_f[p], bi = ij / 3, bj = ij - 3 * bi;
            double v = hA[ij] * G[bj * nn + f];
            s->refv[p] = bi == bj ? C[bi * nn + f] + v : v;
        }
        if (fetestexcept(RC_FLAGS)) return RADAU_BAIL_FLAGS;
        return 0;
    }
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            /* the block C_i + (h A_ij) G_j (diagonal) or (h A_ij) G_j: the
               scalar h A_ij rounded first, as Python's left-to-right
               `h * A[i, j] * Gi[j]` makes it -- h (A_ij G_j) rounds apart */
            const double hA = h * s->A[3 * i + j];
            const double *Gj = G + j * nn, *Ci = C + i * nn;
            double *blk = s->blk;
            if (i == j) for (k = 0; k < nn; k++) blk[k] = Ci[k] + hA * Gj[k];
            else for (k = 0; k < nn; k++) blk[k] = hA * Gj[k];
            /* `remove_row_col`, then the block at rows i m.., columns
               j m.. of Jbig, stored column-major (numpy's solve's copy) */
            for (ri = 0, rr = 0; ri < n; ri++) {
                if (ri == r) continue;
                for (ci = 0, cc = 0; ci < n; ci++) {
                    if (ci == r) continue;
                    s->Jf[(j * m + cc) * M + (i * m + rr)] = blk[ri * n + ci];
                    cc++;
                }
                rr++;
            }
        }
    }
    if (fetestexcept(RC_FLAGS)) return RADAU_BAIL_FLAGS;
    return 0;
}

/* StructuralKLU.solve(Jbig, -R) on its own sequence (pcx_kluseq.h), the
   block's pattern entries in klu.kAx (the assembly's) -- `_ours` takes it:
   the block is zero off the pattern -- the solution in kx.  Returns 0, or
   the bail: RADAU_BAIL_KLU (a factor or klu_solve failed: its
   LinAlgError), RADAU_BAIL_FLAGS (a flag numpy reports in the residual's
   arithmetic) */
static long stage_klu(radau_t *s, long M)
{
    long st = pcx_klu_step(&s->klu, s->Ad, s->b, s->kx, s->ky, s->kt, M, s->dgemv, 0, NULL);
    if (st == PCX_KLU_FLAGS) return RADAU_BAIL_FLAGS;
    if (st) return RADAU_BAIL_KLU;
    return 0;
}

/* every entry of klu.kAx finite: 1, else 0 */
static long klu_finite(const radau_t *s)
{
    long p;
    for (p = 0; p < s->nact; p++) if (!isfinite(s->klu.kAx[p])) return 0;
    return 1;
}

/* the undamped coupled Newton from `seed`.
   Assumes what `_tran_radau_c.solve` checks and `_Ctx` sizes: n >= 2,
   0 <= iref < n, the buffers of radau.h's sizes, the core bound to this
   circuit and its readiness probed (T), walk_fn NULL or a walk ready.
   Per iteration: dY = solve(Jbig, -R) by dgesv, or by StructuralKLU's
   sequence where klu.klu_sym is bound (`stage_klu`); per stage Y_prev +
   insert(dY_i, iref, 0.0), limited by the walk against Y_prev, and the
   limited step's largest reduced entry; the assembly at the new stages;
   then the test -- the Python assembles before it tests.
   Returns RADAU_CONVERGED with Ya..Ga slots 0 .. assemblies-1 written
   (radau.h), `assemblies`, `iters` and `status` set, and with the KLU bound
   fin_ok.  Else a bail code
   (radau_status), also in `status`: `assemblies` 0, the buffers partly
   written and not read by the glue, the flags cleared, the KLU's numeric
   the entry's.  Outside the
   struct's buffers only the core's and the walk's own writes are made (the
   core's bins, its scratch; the temperature in the elements' parameter
   rows), as the Python's calls of the same C make them. */
long pcx_radau(radau_t *s)
{
    long n = s->n, m = n - 1, r = s->iref, M = 3 * m, a = 0, it, i, k, e;
    int64_t M64 = M, one64 = 1, info64 = 0;
    walk_fn_t walk = (walk_fn_t) s->walk_fn;
    const double *d;
    s->klu.refactors = s->klu.factors = s->klu.pivot_fallbacks = 0;
    s->klu.residual_fallbacks = 0;
    s->klu.klu_out = s->klu.klu_num;
    s->assemblies = 0;
    s->iters = 0;
    s->status = RADAU_RUNNING;
    /* every input finite: anything else is the Python Newton's */
    if (!(isfinite(s->h) && isfinite(s->reltol) && isfinite(s->abstol))) return rc_bail(s, RADAU_BAIL_NONFINITE);
    for (k = 0; k < 9; k++) if (!isfinite(s->A[k])) return rc_bail(s, RADAU_BAIL_NONFINITE);
    for (k = 0; k < n; k++) if (!isfinite(s->qn[k])) return rc_bail(s, RADAU_BAIL_NONFINITE);
    for (k = 0; k < 3 * n; k++) if (!isfinite(s->u[k])) return rc_bail(s, RADAU_BAIL_NONFINITE);
    memcpy(s->Ya, s->seed, (size_t) (3 * n) * sizeof(double));
    e = assemble(s, 0);
    if (e) return rc_bail(s, e);
    for (it = 0; it < s->maxiter; it++) {
        double scale = 0.0, ynorm = 0.0;
        double *Yp = s->Ya + a * 3 * n, *Yn;
        /* (the new stages go to slot a + 1: past the buffers, hand back) */
        if (a + 1 >= s->maxa) return rc_bail(s, RADAU_BAIL_ASSEMBLIES);
        Yn = Yp + 3 * n;
        /* solve(Jbig, -R), then sum|R| (the undamped Newton reads it for its
           overflow alone: a bound, numpy's pairwise sum not mirrored).
           Every |R_k| below 2^1000, so no sum of 3m of them overflows; a
           larger entry or a NaN hands back.  A non-finite Jacobian entry
           hands back too (with the KLU bound, every entry off its pattern is
           an exact zero). */
        for (k = 0; k < M; k++) {
            if (!(fabs(s->R[k]) < 0x1p1000)) return rc_bail(s, RADAU_BAIL_RBOUND);
            s->b[k] = -s->R[k];
        }
        if (s->klu.klu_sym) {
            if (!klu_finite(s)) return rc_bail(s, RADAU_BAIL_RBOUND);
            e = stage_klu(s, M);
            if (e) return rc_bail(s, e);
            d = s->kx;
        } else {
            for (k = 0; k < M * M; k++) if (!isfinite(s->Jf[k])) return rc_bail(s, RADAU_BAIL_RBOUND);
            ((dgesv64_t) s->dgesv)(&M64, &one64, s->Jf, &M64, s->ipiv, s->b, &M64, &info64);
            /* (numpy's solve clears what LAPACK raised; an info is its LinAlgError) */
            feclearexcept(FE_ALL_EXCEPT);
            if (info64 != 0) return rc_bail(s, RADAU_BAIL_DGESV);
            d = s->b;
        }
        for (k = 0; k < M; k++) if (!isfinite(d[k])) return rc_bail(s, RADAU_BAIL_STEP);
        for (i = 0; i < 3; i++) {
            const double *yp = Yp + i * n, *di = d + i * m;
            double *yn = Yn + i * n, mx = 0.0;
            /* Y_prev + insert(dY_i, iref, 0.0), Y_prev the left operand */
            for (k = 0; k < r; k++) yn[k] = yp[k] + di[k];
            /* the reference row adds the inserted 0.0, not a copy: a -0.0
               there becomes +0.0, as numpy's elementwise sum makes it */
            yn[r] = yp[r] + 0.0;
            for (k = r + 1; k < n; k++) yn[k] = yp[k] + di[k - 1];
            if (fetestexcept(RC_FLAGS)) return rc_bail(s, RADAU_BAIL_FLAGS);
            /* `cir.limit(Y_trial, Y_prev)`: the walk limits yn in place and
               returns wn, or the element whose limiter declined (the
               Python limits); its flags are not tested */
            if (walk) {
                e = walk(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK, yn, yp, 0, s->wn);
                if (e < s->wn) return rc_bail(s, RADAU_BAIL_WALKSTOP);
                feclearexcept(FE_ALL_EXCEPT);
            }
            /* step_i = red(Y_new - Y_prev) (stored: its flags tested) */
            for (k = 0; k < n; k++) s->stp[k] = yn[k] - yp[k];
            if (fetestexcept(RC_FLAGS)) return rc_bail(s, RADAU_BAIL_FLAGS);
            /* np.max(np.abs(step_i)), then Python's max(scale, it): each
               takes a later value only where it is greater.  Every value is
               finite and >= +0.0 here, so starting from 0.0 and any tie
               give the same bits */
            for (k = 0; k < n; k++) {
                double v;
                if (k == r) continue;
                v = fabs(s->stp[k]);
                if (!isfinite(v)) return rc_bail(s, RADAU_BAIL_STEP);
                if (v > mx) mx = v;
            }
            if (mx > scale) scale = mx;
        }
        a++;
        e = assemble(s, a);
        if (e) return rc_bail(s, e);
        /* `_stages_converged`: ynorm = max over stages of np.max(np.abs(
           red(Y_i))) -- finite (`assemble` checked the stages) and >= +0.0,
           so the running maximum from 0.0 is numpy's */
        for (i = 0; i < 3; i++)
            for (k = 0; k < n; k++) {
                double v;
                if (k == r) continue;
                v = fabs(Yn[i * n + k]);
                if (v > ynorm) ynorm = v;
            }
        /* scale <= reltol * ynorm + abstol: two roundings (no fma), each
           stored before the test (numpy's scalar arithmetic reports their
           flags) */
        feclearexcept(FE_ALL_EXCEPT);
        s->tolv[0] = s->reltol * ynorm;
        s->tolv[1] = s->tolv[0] + s->abstol;
        if (fetestexcept(RC_FLAGS)) return rc_bail(s, RADAU_BAIL_FLAGS);
        s->iters = it + 1;
        if (scale <= s->tolv[1]) {
            if (s->klu.klu_sym) {
                pcx_klu_finish(&s->klu, 0);
                s->fin_ok = klu_finite(s);
            }
            s->assemblies = a + 1;
            s->status = RADAU_CONVERGED;
            return RADAU_CONVERGED;
        }
    }
    return rc_bail(s, RADAU_BAIL_MAXITER);
}
