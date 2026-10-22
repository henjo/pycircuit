/* pcnrdc.c -- PCNR's DC solve in C (contract: pcnrdc.h): `pcnr.solve_dc`,
   statement by statement, or handed back.  Its glue is `_pcnr_c`.

   THE FLAG RULE (pcx_fenv.h): what the core's kernels and the devices'
   kernels raise inside is cleared after them -- the Python's calls of them
   (cffi) never report it; a flag numpy would report (the assembly's, the
   Schur updates', the step's, the refine's, the convergence tests'
   arithmetic) hands the solve back, and the Python's solve gives its
   warnings.  A NaN is carried quietly, as numpy carries it, from the
   assembly through the solve (dense, or KLU's) to x_new, whose test the
   Python makes (PCNRDC_NONFINITE); an infinite KLU entry or rhs, and a
   value that is not finite in v_new or g_lim, hands back.  Every return
   leaves the flags clear. */
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <fenv.h>
#include "pcx_abi.h"
#include "pcx_minmax.h"
#include "pcx_fenv.h"
#include "core.h"
#include "pcx_klu.h"
#include "pcx_verdict.h"
#include "pcx_kluseq.h"
#include "pcnrdc.h"

/* the flags numpy reports here (pcx_kluseq.h) */
#define PC_FLAGS PCX_NP_FLAGS

/* one participant at its own unknowns (`PcnrDevice._kernel_stamp`): the local
   point by its walk and mirrors, its i and G kernels there, the affine
   remainder taken off, `block[:, j] = sum_r G_nl[:, r] T[r, j]`; its current
   into F (`np.add.at`'s order), its g_lim (`augmented_system`) */
static void device(pcnrdc_t *s, long d)
{
    long t = s->d_t[d], m = s->d_m[d], off = s->d_off[d], i, j, r, c;
    const int64_t *rows = s->rows + s->d_rows[d];
    const int64_t *wk = s->walk + 4 * s->d_walk[d], *mr = s->mir + 2 * s->d_mir[d];
    const double *Tm = s->Tm + s->d_tm[d];
    double *blk = s->blk + s->d_tm[d], *pot = s->pot, *cur = s->cur, *G = s->Gl;
    const double *v = s->v + off;
    for (r = 0; r < t; r++) pot[r] = 0.0;
    for (i = 0; i < s->n_walk[d]; i++)
        pot[wk[4 * i]] = pot[wk[4 * i + 1]] + (double) wk[4 * i + 3] * v[wk[4 * i + 2]];
    for (i = 0; i < s->n_mir[d]; i++) pot[mr[2 * i]] = pot[mr[2 * i + 1]];
    if (s->ti[d] >= 0) s->prow[d][s->ti[d]] = s->T;
    ((pass_fn_t) s->ifn[d])(pot, s->prow[d], cur);
    ((pass_fn_t) s->Gfn[d])(pot, s->prow[d], G);
    /* (the kernels' flags: the Python's cffi calls never report them) */
    feclearexcept(FE_ALL_EXCEPT);
    if (s->d_lin[d] >= 0) {
        const double *L = s->lin + s->d_lin[d];
        for (r = 0; r < t; r++) {
            double acc = 0.0;
            for (c = 0; c < t; c++) acc += L[r * t + c] * pot[c];
            cur[r] -= acc;
        }
        for (i = 0; i < t * t; i++) G[i] -= L[i];
    }
    for (i = 0; i < t * m; i++) blk[i] = 0.0;
    for (j = 0; j < m; j++)
        for (r = 0; r < t; r++) {
            double tr = Tm[r * m + j];
            if (tr != 0.0)
                for (i = 0; i < t; i++) blk[i * m + j] += G[i * t + r] * tr;
        }
    for (r = 0; r < t; r++) s->F[rows[r]] += cur[r];
    for (j = 0; j < m; j++) {
        long ra = s->prb[2 * (off + j)], rb = s->prb[2 * (off + j) + 1];
        s->gl[off + j] = v[j] - (s->x[ra] - s->x[rb]);
    }
}

/* `augmented_system` and the gshunt rows: the view's core at x (F = i + u,
   J = G), each participant, then `g[shunt] += gshunt x[shunt]`,
   `J[shunt, shunt] += gshunt`, then `gjgmin` across each junction pair
   (`F[a] += g vj`, `F[b] -= g vj`, J's four entries, in pair order).  F and J may carry a NaN on (the kernels'
   at a wild point); g_lim, from the finite x and v, may not.  0, or
   PCNRDC_BAIL_FLAGS */
static long assemble(pcnrdc_t *s)
{
    const hdl_core_t *c = (const hdl_core_t *) s->core;
    long n = s->n, i, d;
    ((core_fn_t) s->core_fn)(s->core, s->x, s->T, PCX_PASS_G | PCX_PASS_I | s->pattern,
                             PCX_FORMULA_PASSES, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL, NULL,
                             NULL, NULL, s->Cs, s->q, NULL, NULL, NULL, NULL);
    feclearexcept(FE_ALL_EXCEPT);
    for (i = 0; i < n; i++) s->F[i] = c->bini[i] + s->u[i];
    if (s->np)
        for (i = 0; i < s->np; i++) s->J[s->P[i]] = c->binG[s->P[i]];
    else
        for (i = 0; i < n * n; i++) s->J[i] = c->binG[i];
    for (d = 0; d < s->nd; d++) device(s, d);
    if (s->gshunt != 0.0) {
        for (i = 0; i < s->nn; i++) {
            if (i == s->iref) continue;
            s->F[i] = s->F[i] + s->gshunt * s->x[i];
            s->J[i * n + i] = s->J[i * n + i] + s->gshunt;
        }
    }
    if (s->gjgmin != 0.0) {
        double g = s->gjgmin;
        long j;
        for (j = 0; j < s->nj; j++) {
            long a = s->ja[j], b = s->jb[j];
            double vj = s->x[a] - s->x[b];
            s->F[a] = s->F[a] + g * vj;
            s->F[b] = s->F[b] - g * vj;
            s->J[a * n + a] = s->J[a * n + a] + g;
            s->J[a * n + b] = s->J[a * n + b] - g;
            s->J[b * n + a] = s->J[b * n + a] - g;
            s->J[b * n + b] = s->J[b * n + b] + g;
        }
    }
    if (fetestexcept(PC_FLAGS)) return PCNRDC_BAIL_FLAGS;
    for (i = 0; i < s->k; i++) if (!isfinite(s->gl[i])) return PCNRDC_BAIL_FLAGS;
    return 0;
}

/* `schur_reduce` in place: each unknown's column added at its probe's ra
   column, subtracted at rb (`np.add.at` / `np.subtract.at` on the device's
   rows), `col * g_lim` into the rhs corrections; f_eff = F - rc */
static long schur(pcnrdc_t *s)
{
    long n = s->n, i, d, j, r;
    for (i = 0; i < n; i++) s->rc[i] = 0.0;
    for (d = 0; d < s->nd; d++) {
        long t = s->d_t[d], m = s->d_m[d], off = s->d_off[d];
        const int64_t *rows = s->rows + s->d_rows[d];
        const double *blk = s->blk + s->d_tm[d];
        for (j = 0; j < m; j++) {
            long ra = s->prb[2 * (off + j)], rb = s->prb[2 * (off + j) + 1];
            double g = s->gl[off + j];
            for (r = 0; r < t; r++) s->J[rows[r] * n + ra] += blk[r * m + j];
            for (r = 0; r < t; r++) s->J[rows[r] * n + rb] -= blk[r * m + j];
            for (r = 0; r < t; r++) s->rc[rows[r]] += blk[r * m + j] * g;
        }
    }
    for (i = 0; i < n; i++) s->fe[i] = s->F[i] - s->rc[i];
    if (fetestexcept(PC_FLAGS)) return PCNRDC_BAIL_FLAGS;
    return 0;
}

/* `predict`'s solve, numpy.linalg.solve on the reduced Schur matrix:
   dgesv on its column-major copy, rhs `-f_eff` (in `b`).  An info is
   numpy's LinAlgError (PCNRDC_SINGULAR, served: the glue raises it); what
   dgesv raises with info 0 numpy's linalg clears (`set_fp_invalid_or_clear`),
   so it is cleared here -- a step that is not finite is `update`'s
   hand-back.  The solution in `b` */
static long solve_dense(pcnrdc_t *s)
{
    long n = s->n, m = s->m, i, j;
    const int64_t *red = s->red;
    int64_t m64 = m, one64 = 1, info64 = 0;
    for (i = 0; i < n; i++) {
        if (red[i] < 0) continue;
        for (j = 0; j < n; j++)
            if (red[j] >= 0) s->A[red[j] * m + red[i]] = s->J[i * n + j];
    }
    ((dgesv64_t) s->dgesv)(&m64, &one64, s->A, &m64, s->ipiv64, s->b, &m64, &info64);
    feclearexcept(FE_ALL_EXCEPT);
    return info64 != 0 ? PCNRDC_SINGULAR : 0;
}

/* StructuralKLU.solve(schur[keep, keep], -f_eff): the reduced matrix (R),
   `_ours`'s test (as many nonzeros on the pattern as in all of it, else
   its fallback solve: handed back), its entries in the pattern's order, its
   own sequence (pcx_kluseq.h).  A NaN entry or rhs goes through: KLU is the
   same library on the same values, and the sequence's residual guard
   propagates a NaN as numpy's max does.  An infinite one hands back (the
   residual's `arr @ x - rhs` may be numpy's warning).  A failed factor is
   its LinAlgError (PCNRDC_SINGULAR).  The solution in `kx`.
   On the pattern (`kfull`): the entries gathered from J where they stand,
   the foreign test over the pattern's entries off KLU's (J is zero
   elsewhere: the same count), the residual's test on the CSC -- R is not
   made */
static long solve_klu(pcnrdc_t *s)
{
    long n = s->n, m = s->m, i, j, nz = 0, nzp = 0, fin = 1, st;
    const int64_t *red = s->red;
    if (s->kfull) {
        for (i = 0; i < s->nact; i++) {
            s->klu.kAx[i] = s->J[s->kfull[i]];
            fin &= isinf(s->klu.kAx[i]) == 0;
        }
        for (i = 0; i < s->nfor; i++)
            if (s->J[s->ffull[i]] != 0.0) return PCNRDC_BAIL_FOREIGN;
    } else {
        for (i = 0; i < n; i++) {
            if (red[i] < 0) continue;
            for (j = 0; j < n; j++)
                if (red[j] >= 0) s->R[red[i] * m + red[j]] = s->J[i * n + j];
        }
        for (i = 0; i < m * m; i++) nz += s->R[i] != 0.0;
        for (i = 0; i < s->nact; i++) {
            s->klu.kAx[i] = s->R[s->act_r[i]];
            nzp += s->klu.kAx[i] != 0.0;
            fin &= isinf(s->klu.kAx[i]) == 0;
        }
        if (nz != nzp) return PCNRDC_BAIL_FOREIGN;
    }
    for (i = 0; i < m; i++) fin &= isinf(s->b[i]) == 0;
    if (!fin) return PCNRDC_BAIL_FLAGS;
    st = pcx_klu_step(&s->klu, s->R, s->b, s->kx, s->ky, s->kt, m, s->dgemv, s->kfull == NULL,
                      &s->dense_resid);
    if (st == PCX_KLU_FACTOR) return PCNRDC_SINGULAR;
    if (st == PCX_KLU_SOLVE) return PCNRDC_BAIL_KLU;
    if (st == PCX_KLU_FLAGS) return PCNRDC_BAIL_FLAGS;
    return 0;
}

/* the step dx_MNA in `dx` (the reference row and each held row 0.0):
   `-f_eff` reduced, the DC's solve.  0, or why not */
static long solve(pcnrdc_t *s)
{
    long n = s->n, i, st;
    const int64_t *red = s->red;
    const double *x;
    for (i = 0; i < n; i++) if (red[i] >= 0) s->b[red[i]] = -s->fe[i];
    st = s->klu.klu_sym ? solve_klu(s) : solve_dense(s);
    if (st != 0) return st;
    x = s->klu.klu_sym ? s->kx : s->b;
    for (i = 0; i < n; i++) s->dx[i] = red[i] < 0 ? 0.0 : x[red[i]];
    return 0;
}

/* x_new = x + dx (the reference 0.0); v_new = v + dx_lim, dx_lim =
   -(g_lim + (-dx[ra] + dx[rb])); x_new not finite is the Python's
   RuntimeError (PCNRDC_NONFINITE); then each device's refine (its class's
   kernel at its last accepted sub-vector).  0, or why not */
static long update(pcnrdc_t *s)
{
    long n = s->n, i, d, j, r;
    for (i = 0; i < n; i++) s->xn[i] = s->x[i] + s->dx[i];
    s->xn[s->iref] = 0.0;
    for (i = 0; i < s->k; i++) {
        long ra = s->prb[2 * i], rb = s->prb[2 * i + 1];
        s->vn[i] = s->v[i] + -(s->gl[i] + (-s->dx[ra] + s->dx[rb]));
    }
    if (fetestexcept(PC_FLAGS)) return PCNRDC_BAIL_FLAGS;
    for (i = 0; i < n; i++) if (!isfinite(s->xn[i])) return PCNRDC_NONFINITE;
    for (d = 0; d < s->nd; d++) {
        long t = s->d_t[d], m = s->d_m[d], off = s->d_off[d];
        const int64_t *rows = s->rows + s->d_rows[d];
        for (r = 0; r < t; r++) s->xs[r] = s->x[rows[r]];
        for (j = 0; j < m; j++) { s->rio[j] = s->vn[off + j]; s->rio[m + j] = s->v[off + j]; }
        if (s->ti[d] >= 0) s->prow[d][s->ti[d]] = s->T;
        if (((lim_fn_t) s->rfn[d])(s->xs, s->prow[d], s->rio)) return PCNRDC_BAIL_REFINE;
        if (fetestexcept(PC_FLAGS)) return PCNRDC_BAIL_FLAGS;
        for (j = 0; j < m; j++) s->vn[off + j] = s->rio[j];
    }
    for (i = 0; i < s->k; i++) if (!isfinite(s->vn[i])) return PCNRDC_BAIL_FLAGS;
    return 0;
}

/* conv_f's verdict on the pattern: each kept row's S_i = sum |J_ij|
   |x_new_j| over P (J is zero off it; every row holds its diagonal there),
   in P's order.  A sum rounds only where two operands holding a nonzero
   term meet -- a zero product joins exactly, fused or not -- so with k_i
   entries of P in the row numpy's dgemv (any summation tree of the same
   products, all >= 0) is within 2 gamma_(k_i + 1) S_i of it; with d_i =
   4 (k_i + 2) eps it lies in [S_i (1 - d_i), S_i (1 + d_i)] -- the
   endpoints rounded outward -- and the test, monotone in it, holds where
   it holds at the lower end and fails where it fails at the upper:
   PCX_HOLDS every row, PCX_FAILS one, or PCX_UNDECIDED -- also where a sum
   is not finite or reaches 1e300 (there numpy's flags could be raised).
   Where it decides, no flag numpy reports is raised here or in the dense
   test */
static int cf_sparse(pcnrdc_t *s)
{
    long n = s->n, e = 0, r, k;
    int open = 0;
    while (e < s->np) {
        int64_t base = (s->P[e] / n) * n;
        double S = 0.0, af, d;
        r = (long) (base / n);
        for (k = 0; e < s->np && s->P[e] < base + n; e++, k++)
            S += fabs(s->J[s->P[e]]) * fabs(s->xn[s->P[e] - base]);
        if (s->red[r] < 0) continue;
        if (!(S < 1e300)) return PCX_UNDECIDED;
        d = 4.0 * ((double) k + 2.0) * DBL_EPSILON;
        af = fabs(s->fe[r]);
        if (af < s->reltol * ((S - d * S) + af) + s->iabstol) continue;
        if (!(af < s->reltol * ((S + d * S) + af) + s->iabstol)) return PCX_FAILS;
        open = 1;
    }
    return open ? PCX_UNDECIDED : PCX_HOLDS;
}

/* the three tests (`solve_dc`): conv_x, conv_f with `I = |J_eff| |x_new|`
   (numpy's own dgemv on the full |J_eff|: BLAS's summation order; on the
   pattern its verdict, `cf_sparse`, where that decides) + |f_eff| -- each
   over the reduced rows (`red`) -- and `lim_converged`.  1 converged, 0
   not, -1 a flag */
static int converged(pcnrdc_t *s)
{
    long n = s->n, i;
    int cx = 1, cf = 1, cl = 1, verdict = PCX_UNDECIDED;
    fexcept_t entry;
    for (i = 0; i < n; i++) {
        if (s->red[i] < 0) continue;
        if (!(fabs(s->dx[i]) < s->reltol * pcx_np_maximum(fabs(s->xn[i]), fabs(s->x[i]))
                               + s->abstol)) cx = 0;
    }
    if (s->np) {
        fegetexceptflag(&entry, FE_ALL_EXCEPT);
        verdict = cf_sparse(s);
        if (verdict == PCX_UNDECIDED) {
            /* (the dense test: what the sparse sums raised is not the Python's) */
            fesetexceptflag(&entry, FE_ALL_EXCEPT);
            s->dense_conv++;
        }
        cf = verdict == PCX_HOLDS;
    }
    if (verdict == PCX_UNDECIDED) {
        cf = 1;
        for (i = 0; i < n * n; i++) s->Jabs[i] = fabs(s->J[i]);
        for (i = 0; i < n; i++) s->ax[i] = fabs(s->xn[i]);
        ((dgemv64_t) s->dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, n, n, 1.0, s->Jabs, n,
                               s->ax, 1, 0.0, s->Is, 1);
        for (i = 0; i < n; i++) {
            double af = fabs(s->fe[i]);
            if (s->red[i] < 0) continue;
            if (!(af < s->reltol * (s->Is[i] + af) + s->iabstol)) cf = 0;
        }
    }
    for (i = 0; i < s->k; i++)
        if (!(fabs(s->gl[i]) < s->reltol * pcx_np_maximum(fabs(s->vn[i]), 1.0) + s->abstol))
            cl = 0;
    if (fetestexcept(PC_FLAGS)) return -1;
    return cx && cf && cl;
}

/* the loop (`solve_dc`): assemble, the Schur step, the update, the tests;
   `x, v = x_new, v_new` after each.  CONVERGED, MAXITER, SINGULAR, or a
   hand-back */
static long iterate(pcnrdc_t *s)
{
    long it, i, st;
    int conv;
    for (it = 0; it < s->maxiter; it++) {
        if ((st = assemble(s)) != 0) return st;
        if ((st = schur(s)) != 0) return st;
        if (s->choose) return PCNRDC_CHOOSE;
        if ((st = solve(s)) != 0) return st;
        if ((st = update(s)) != 0) return st;
        if ((conv = converged(s)) < 0) return PCNRDC_BAIL_FLAGS;
        s->iters = it + 1;
        for (i = 0; i < s->n; i++) s->x[i] = s->xn[i];
        for (i = 0; i < s->k; i++) s->v[i] = s->vn[i];
        if (conv) return PCNRDC_CONVERGED;
    }
    return PCNRDC_MAXITER;
}

/* the entry (contract: pcnrdc.h) */
long pcx_pcnrdc(pcnrdc_t *s)
{
    long st;
    s->iters = 0;
    s->dense_conv = s->dense_resid = 0;
    s->klu.refactors = s->klu.factors = s->klu.pivot_fallbacks = 0;
    s->klu.residual_fallbacks = 0;
    s->klu.klu_out = s->klu.klu_num;
    feclearexcept(FE_ALL_EXCEPT);
    st = iterate(s);
    if (s->klu.klu_sym) pcx_klu_finish(&s->klu, st >= PCNRDC_CHOOSE);
    feclearexcept(FE_ALL_EXCEPT);
    return s->status = st;
}
