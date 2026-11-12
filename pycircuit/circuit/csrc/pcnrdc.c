/* pcnrdc.c -- PCNR's DC solve in C (contract: pcnrdc.h): `pcnr.solve_dc`,
   statement by statement, or handed back.  Its glue is `_pcnr_c`.

   THE FLAG RULE (pcx_fenv.h): what the core's kernels and the devices'
   kernels raise inside is cleared after them -- the Python's calls of them
   (cffi) never report it; a flag numpy would report (the assembly's, the
   step's, the refine's, the convergence tests' arithmetic) hands the solve
   back, and the Python's solve gives its warnings.  The Schur sum's are not
   reported (`schur_reduce`'s errstate ignores its overflow and invalid) and
   are cleared.  A NaN or an infinity is carried quietly, as numpy carries
   it, from the assembly through the solve (dense, or KLU's -- its residual's
   dense test handing back where numpy's would warn) to x_new, whose test
   the Python makes (PCNRDC_NONFINITE); a value that is not finite in v_new
   or g_lim hands back.  Every return leaves the flags clear. */
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
#include "pcx_held.h"
#include "pcx_kluseq.h"
#include "pcnrdc.h"

/* the flags numpy reports here (pcx_kluseq.h) */
#define PC_FLAGS PCX_NP_FLAGS

/* `feclearexcept(FE_ALL_EXCEPT)` where a device's kernel raised a flag this
   unit tests (`pc_clear_raised`): x87's flags
   (`fnclex`) and MXCSR's (FE_ALL_EXCEPT's bits), what glibc's clears --
   every flag `fetestexcept` and `fegetexceptflag` read the same after
   (`fnclex` drops x87's denormal bit too, which FE_ALL_EXCEPT excludes) --
   without its `fnstenv` / `fldenv` (microcoded, the most of its cost).
   Each asm clobbers
   memory, so it is ordered with the kernel's call before it and the
   device's arithmetic after it as the library call is (THE FLAG RULE) */
static inline void pc_clear_flags(void)
{
#if defined(__x86_64__) && defined(__GNUC__)
    unsigned int csr;
    __asm__ volatile ("fnclex" : : : "memory");
    __asm__ volatile ("stmxcsr %0" : "=m" (csr) : : "memory");
    csr &= ~(unsigned int) FE_ALL_EXCEPT;
    __asm__ volatile ("ldmxcsr %0" : : "m" (csr) : "memory");
#else
    feclearexcept(FE_ALL_EXCEPT);
#endif
}

/* a kernel's flags taken off where it raised one this unit tests
   (PC_FLAGS, x87's or MXCSR's: what `fetestexcept(PC_FLAGS)` reads), all of
   them then (`pc_clear_flags`); inexact, underflow and denormal, which
   nearly every kernel raises and no test here reads, left as they are --
   so every PC_FLAGS test after reads what it read with the clear made
   every time, without its serializing MXCSR write on every device */
static inline void pc_clear_raised(void)
{
#if defined(__x86_64__) && defined(__GNUC__)
    unsigned int csr;
    unsigned short sw;
    __asm__ volatile ("stmxcsr %0" : "=m" (csr) : : "memory");
    __asm__ volatile ("fnstsw %0" : "=m" (sw) : : "memory");
    if (((unsigned int) sw | csr) & (unsigned int) PC_FLAGS) pc_clear_flags();
#else
    if (fetestexcept(PC_FLAGS)) feclearexcept(FE_ALL_EXCEPT);
#endif
}

/* one participant at its own unknowns (`PcnrDevice._kernel_stamp`): the local
   point by its walk and mirrors, its class's kernel there (i and G in one
   call), the affine remainder taken off, `block[:, j] = sum_r G_nl[:, r]
   T[r, j]`; its current into F (`np.add.at`'s order), its g_lim
   (`augmented_system`) */
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
    ((fused_fn_t) s->kfn[d])(pot, s->prow[d], PCX_PASS_G | PCX_PASS_I, G, 0, cur, 0);
    /* (the kernels' flags: the Python's cffi calls never report them) */
    pc_clear_raised();
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
   J = G: packed, its pattern's entries), each participant, then `g[shunt] +=
   gshunt x[shunt]`, `J[shunt, shunt] += gshunt`, then `gjgmin` across each
   junction pair (`F[a] += g vj`, `F[b] -= g vj`, J's four entries, in pair
   order) -- J at its positions (pcnrdc.h's THE PATTERN).  F and J may carry
   a NaN on (the kernels' at a wild point); g_lim, from the finite x and v,
   may not.  0, or PCNRDC_BAIL_FLAGS */
static long assemble(pcnrdc_t *s)
{
    const hdl_core_t *c = (const hdl_core_t *) s->core;
    double *J = s->J;
    long n = s->n, i, d;
    ((core_fn_t) s->core_fn)(s->core, s->x, s->T, PCX_PASS_G | PCX_PASS_I | s->pattern,
                             PCX_FORMULA_PASSES, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL, NULL,
                             NULL, NULL, s->Cs, s->q, NULL, NULL, NULL, NULL);
    feclearexcept(FE_ALL_EXCEPT);
    for (i = 0; i < n; i++) s->F[i] = c->bini[i] + s->u[i];
    /* (with the pattern always THE PACKED G: `gmap` set where P is) */
    if (s->gmap)
        for (i = 0; i < s->np; i++) J[i] = s->gmap[i] >= 0 ? c->binGp[s->gmap[i]] : 0.0;
    else
        for (i = 0; i < n * n; i++) J[i] = c->binG[i];
    for (d = 0; d < s->nd; d++) device(s, d);
    if (s->gshunt != 0.0) {
        for (i = 0; i < s->nn; i++) {
            if (i == s->iref) continue;
            s->F[i] = s->F[i] + s->gshunt * s->x[i];
            J[s->dmap[i]] = J[s->dmap[i]] + s->gshunt;
        }
    }
    if (s->gjgmin != 0.0) {
        double g = s->gjgmin;
        long j;
        for (j = 0; j < s->nj; j++) {
            long a = s->ja[j], b = s->jb[j];
            const int64_t *q = s->jmap + 4 * j;
            double vj = s->x[a] - s->x[b];
            s->F[a] = s->F[a] + g * vj;
            s->F[b] = s->F[b] - g * vj;
            J[q[0]] = J[q[0]] + g;
            J[q[1]] = J[q[1]] - g;
            J[q[2]] = J[q[2]] - g;
            J[q[3]] = J[q[3]] + g;
        }
    }
    if (fetestexcept(PC_FLAGS)) return PCNRDC_BAIL_FLAGS;
    for (i = 0; i < s->k; i++) if (!isfinite(s->gl[i])) return PCNRDC_BAIL_FLAGS;
    return 0;
}

/* `schur_reduce` in place: each unknown's column added at its probe's ra
   column, subtracted at rb (`np.add.at` / `np.subtract.at` on the device's
   rows: J's positions `smap`), `col * g_lim` into the rhs corrections;
   f_eff = F - rc.  What the sum raises -- an overflow, an invalid: a
   derivative run away at a run-away iterate -- `schur_reduce` does not
   report (its errstate ignores both): the infinity or NaN is kept for the
   solve to fail on, as the Python's, and the flags are cleared so the steps
   after read their own.  Always 0 */
static long schur(pcnrdc_t *s)
{
    long n = s->n, i, d, j, r;
    const int64_t *sm = s->smap;
    double *J = s->J;
    for (i = 0; i < n; i++) s->rc[i] = 0.0;
    for (d = 0; d < s->nd; d++) {
        long t = s->d_t[d], m = s->d_m[d], off = s->d_off[d];
        const int64_t *rows = s->rows + s->d_rows[d];
        const double *blk = s->blk + s->d_tm[d];
        for (j = 0; j < m; j++, sm += 2 * t) {
            double g = s->gl[off + j];
            for (r = 0; r < t; r++) J[sm[r]] += blk[r * m + j];
            for (r = 0; r < t; r++) J[sm[t + r]] -= blk[r * m + j];
            for (r = 0; r < t; r++) s->rc[rows[r]] += blk[r * m + j] * g;
        }
    }
    for (i = 0; i < n; i++) s->fe[i] = s->F[i] - s->rc[i];
    feclearexcept(FE_ALL_EXCEPT);
    return 0;
}

/* `predict`'s solve, numpy.linalg.solve on the reduced Schur matrix:
   dgesv on its column-major copy (packed J: zero, then P's entries, as
   the n x n J is zero off P), rhs `-f_eff` (in `b`).  An info is
   numpy's LinAlgError (PCNRDC_SINGULAR, served: the glue raises it); what
   dgesv raises with info 0 numpy's linalg clears (`set_fp_invalid_or_clear`),
   so it is cleared here -- a step that is not finite is `update`'s
   hand-back.  The solution in `b` */
static long solve_dense(pcnrdc_t *s)
{
    long n = s->n, m = s->m, i, j;
    const int64_t *red = s->red;
    int64_t m64 = m, one64 = 1, info64 = 0;
    if (s->np) {
        for (i = 0; i < m * m; i++) s->A[i] = 0.0;
        for (i = 0; i < s->np; i++) {
            long r = (long) (s->P[i] / n), c = (long) (s->P[i] % n);
            if (red[r] >= 0 && red[c] >= 0) s->A[red[c] * m + red[r]] = s->J[i];
        }
    } else {
        for (i = 0; i < n; i++) {
            if (red[i] < 0) continue;
            for (j = 0; j < n; j++)
                if (red[j] >= 0) s->A[red[j] * m + red[i]] = s->J[i * n + j];
        }
    }
    ((dgesv64_t) s->dgesv)(&m64, &one64, s->A, &m64, s->ipiv64, s->b, &m64, &info64);
    feclearexcept(FE_ALL_EXCEPT);
    return info64 != 0 ? PCNRDC_SINGULAR : 0;
}

/* StructuralKLU.solve(schur[keep, keep], -f_eff): the reduced matrix (R),
   `_ours`'s test (as many nonzeros on the pattern as in all of it, else
   its fallback solve: handed back), its entries in the pattern's order, its
   own sequence (pcx_kluseq.h).  A NaN or an infinite entry or rhs goes
   through: KLU is the same library on the same values, the sequence's
   residual guard propagates a NaN as numpy's max does, and its dense test
   hands back where the residual's `arr @ x - rhs` raises numpy's warning
   (an infinity times a zero).  A failed factor is its LinAlgError
   (PCNRDC_SINGULAR).  The solution in `kx`.
   On the pattern (`kfull`): the entries gathered from J's positions (an
   entry off P the zero J[np]), the foreign test over the pattern's entries
   off KLU's (J is zero elsewhere: the same count), the residual's test on
   the CSC, its bound first -- R is not made (nor read) */
static long solve_klu(pcnrdc_t *s)
{
    long n = s->n, m = s->m, i, j, nz = 0, nzp = 0, st;
    const int64_t *red = s->red;
    if (s->kfull) {
        for (i = 0; i < s->nact; i++) s->klu.kAx[i] = s->J[s->kfull[i]];
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
        }
        if (nz != nzp) return PCNRDC_BAIL_FOREIGN;
    }
    st = pcx_klu_step(&s->klu, s->b, s->kx, s->ky, s->kt, m, !s->packed && s->kfull,
                      &s->dense_resid, s->held_v, s->held_st);
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
   |x_new_j| over P (packed J: zero off it; every row holds its diagonal
   there), in P's order.  A sum rounds only where two operands holding a nonzero
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
            S += fabs(s->J[e]) * fabs(s->xn[s->P[e] - base]);
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

/* conv_f on THE PACKED LAYOUT: each kept row's S_i = sum |J_ij| |x_new_j|
   over P in P's order -- row by row, columns ascending: the CSC kernel's
   order within a row (`_packed.PatternMatrix.matvec`), J zero off P, its
   terms adding +0.0 -- the test on S itself.  1 holds, 0 fails */
static int cf_packed(pcnrdc_t *s)
{
    long n = s->n, e = 0, r;
    int ok = 1;
    while (e < s->np) {
        int64_t base = (s->P[e] / n) * n;
        double S = 0.0, af;
        r = (long) (base / n);
        for (; e < s->np && s->P[e] < base + n; e++)
            S += fabs(s->J[e]) * fabs(s->xn[s->P[e] - base]);
        if (s->red[r] < 0) continue;
        af = fabs(s->fe[r]);
        if (!(af < s->reltol * (S + af) + s->iabstol)) ok = 0;
    }
    return ok;
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
    if (s->packed) {
        /* (THE PACKED LAYOUT: the CSC kernel's sums decide, every time) */
        cf = cf_packed(s);
        verdict = cf ? PCX_HOLDS : PCX_FAILS;
    } else if (s->np) {
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
        /* (|J| zero off the pattern as J is: written on it alone, from
           packed J) */
        if (s->np)
            for (i = 0; i < s->np; i++) s->Jabs[s->P[i]] = fabs(s->J[i]);
        else
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

/* `v_lim_init(limit=True)`, where the call gives no v (`seed`): each unknown
   at its branch's voltage at the start, `x[ra] - x[rb]`, then each device's
   refine kernel from v_old 0.0 at that x -- `refine` a device at a time
   through its class's kernel (`_pcnr_c.refine`), the temperature in its
   slot.  A kernel's decline or a flag numpy reports hands back, as
   `update`'s refine does.  0, or why not */
static long seed(pcnrdc_t *s)
{
    long d, j, r;
    for (d = 0; d < s->nd; d++) {
        long t = s->d_t[d], m = s->d_m[d], off = s->d_off[d];
        const int64_t *rows = s->rows + s->d_rows[d];
        for (j = 0; j < m; j++) {
            long ra = s->prb[2 * (off + j)], rb = s->prb[2 * (off + j) + 1];
            s->rio[j] = s->x[ra] - s->x[rb];
            s->rio[m + j] = 0.0;
        }
        for (r = 0; r < t; r++) s->xs[r] = s->x[rows[r]];
        if (s->ti[d] >= 0) s->prow[d][s->ti[d]] = s->T;
        if (((lim_fn_t) s->rfn[d])(s->xs, s->prow[d], s->rio)) return PCNRDC_BAIL_REFINE;
        if (fetestexcept(PC_FLAGS)) return PCNRDC_BAIL_FLAGS;
        for (j = 0; j < m; j++) s->v[off + j] = s->rio[j];
    }
    return 0;
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
    s->held_own = s->klu.klu_sym ? pcx_held_open(s->held_st, s->klu.klu_num) : 0;
    feclearexcept(FE_ALL_EXCEPT);
    st = s->seed ? seed(s) : 0;
    if (st == 0) st = iterate(s);
    if (s->klu.klu_sym) {
        const int back = st >= PCNRDC_CHOOSE;
        pcx_klu_finish(&s->klu, back);
        /* (the held values as the outcome leaves them: pcx_held.h) */
        if (back)
            pcx_held_back(s->held_st, s->held_v, s->klu.kAp[s->m], s->held_own,
                          s->klu.klu_refactor, s->klu.kAp, s->klu.kAi, s->klu.klu_sym,
                          s->klu.klu_common);
        else
            pcx_held_commit(s->held_st, s->held_own);
    }
    feclearexcept(FE_ALL_EXCEPT);
    return s->status = st;
}
