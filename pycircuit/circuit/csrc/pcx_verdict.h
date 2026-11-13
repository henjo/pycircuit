/* pcx_verdict.h -- a test's verdict from a sparse sum and its rounding
 * bound, where the Python makes the test on numpy's dense product (dgemv) of
 * a matrix zero off a pattern: written once for the units that make such a
 * test (pcx_kluseq.h's residual, newton.c, lte.c).
 *
 * A sum of terms rounds only where two operands holding a nonzero term meet
 * -- a zero product joins exactly, fused or not -- so a row with k_i entries
 * on the pattern rounds at most k_i times on any path of any summation tree:
 * its S_i = sum A_ij x_j, in the pattern's order, is within 2 gamma_(k_i + 1)
 * T_i of the dgemv's (T_i = sum |A_ij| |x_j|); e_i = 4 (k_i + 2) eps T_i
 * covers that and the rounding of S_i -+ e_i.  Gradual underflow adds an
 * absolute error -- at most 2^-1075 a product (a sum whose result is
 * subnormal is exact) -- to each side's sum and to e_i itself:
 * 4 (k_i + 2) (eps T_i + 2 DBL_TRUE_MIN) covers it too, so a verdict holds
 * whether or not a product underflowed.  A test monotone in the row's
 * product then holds where it holds at the bound's worse end and fails where
 * it fails at its better end; elsewhere -- and where a value is not finite
 * or a row's sum reaches 1e300 (there numpy's flags could be raised) -- the
 * whole computation decides (the residual's: the CSC kernel's test,
 * `pcx_klu_resid_csc`, `StructuralKLU._csc_residual` -- the bounds below
 * cover a row's sum in any order, its too).
 *
 * THE RESIDUAL'S SUMS WITH SUBNORMALS FLUSHED (`pcx_resid_sparse`): a solve's
 * update near convergence holds subnormal entries, and each operation on one
 * is a microcode assist -- made gradually, the verdict's sums are a sixth of
 * the C transform loop.  They are made with the SSE unit's flush-to-zero and
 * denormals-are-zero set (its control bits alone restored after, its flags as
 * they are), where each operation then lies within DBL_MIN of its gradual
 * result: a term within (|A_ij| + |x_j| + 1) DBL_MIN (an operand or the
 * product flushed), each addition within DBL_MIN -- so the bound grows by
 * (k_i + 2) (amax + xmax + 4) DBL_MIN + 8 DBL_MIN (`amax`, `xmax` the largest
 * |A| and |x|; the 8 the bound's own, |b|'s and the ends' flushes), which
 * covers the gradual underflow's 2^-1075 terms too.  A verdict it makes holds
 * whatever the order's rounding -- the CSC kernel's, as the gradual sums'
 * is.  Where it does not decide (a residual near its threshold, or one
 * near 1e-300 whose threshold the flush could hide), the gradual sums are
 * made (`pcx_resid_gradual`) and theirs is the verdict: the assists paid on
 * those calls alone.
 *
 * WHERE NEITHER DECIDES (`pcx_resid_last`; linearsolver's `_resid_sparse`
 * the same rule):
 *
 * (N) A VALUE NOT FINITE in x, b or A (`pcx_resid_nonfinite`): where A and x
 * hold no infinity, every NaN is quiet and each row's sum of |A_ij| |x_j|
 * over its finite terms, plus a finite |b_i|, is under 1e300, the CSC
 * kernel's sums and their `- b` raise no flag numpy reports (a quiet NaN's
 * arithmetic raises none; no partial sum of finite terms can overflow), and
 * the test FAILS: a NaN in x reaches its column's rows (each column holds an
 * entry), a NaN in A or b its own row, and without a NaN an infinite b_i
 * makes max |.| and max |b| infinite, inf / inf a NaN.  An infinity in A or
 * x stays the kernel's test (0 x inf raises invalid).
 *
 * (W) THE WIDE PASS (`pcx_resid_wide`, where long double holds 64 bits): the
 * row's sums in long double, where a product of doubles neither underflows
 * nor overflows and rounds by v = 2^-64 at most, so S_L lies within
 * gamma^L_k T (+ k 2^-1075: each row's sum is kept as two doubles, its
 * remainder rounded where subnormal) of the exact sum and any double sum of
 * the row's k terms (the CSC kernel's) within gamma_k T + k 2^-1075 of it.  E_i = k u (1 + 2^-6) T_L + 2^-60 |b_i| +
 * (2 k + 4) 2^-1074 (u = 2^-53, T_L the row's sum of |A_ij| |x_j| made the
 * same way) covers both, T's own error, the rounding of (S_L -+ E_i) - b_i
 * in long double and E_i's to double, for every k < 2^31; each end, rounded
 * to double, bounds `fl(sum - b_i)` as the kernel's test makes it (rounding
 * is monotone).  About ten times the narrower band at k ~ 7.
 */
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#if defined(__SSE2__)
#include <xmmintrin.h>
/* the SSE control register's flush-to-zero (bit 15) and denormals-are-zero
   (bit 6) */
#define PCX_FLUSH_BITS 0x8040u
#endif

/* a test's verdict where a rounding bound decides it: it holds, it fails, or
   the dense test decides */
enum pcx_verdict { PCX_HOLDS = 0, PCX_FAILS = 1, PCX_UNDECIDED = 2 };

/* the refactor's residual test, `max |A x - b| / max(max |b|, 1e-300) <=
   tol` (`arr @ x`: numpy's matmul, its dgemv), on A's CSC (`Ap`, `Ai`,
   `Ax`; A zero off the pattern): each row's `fl(dgemv - b)` lies between
   `(S_i - e_i) - b_i` and `(S_i + e_i) - b_i` (rounding is monotone), max
   |.| between the rows' lower and upper ends, and the test -- monotone in
   it, the scale exact -- holds where it holds at the upper end, fails where
   it fails at the lower.  Where it decides, no flag numpy reports is raised
   here or in the dense test.  `y`, `t`, `cnt` (m each) scratch */
__attribute__((noinline)) static int pcx_resid_flushed(const int *Ap, const int *Ai,
                                                        const double *Ax, double tol,
                                                        const double *b, const double *x,
                                                        double *y, double *t, double *cnt, long m);
static int pcx_resid_gradual(const int *Ap, const int *Ai, const double *Ax, double tol,
                             const double *b, const double *x, double *y, double *t,
                             double *cnt, long m);
static int pcx_resid_last(const int *Ap, const int *Ai, const double *Ax, double tol,
                          const double *b, const double *x, double *y, double *t, double *cnt,
                          long m);

static int pcx_resid_sparse(const int *Ap, const int *Ai, const double *Ax, double tol,
                            const double *b, const double *x, double *y, double *t,
                            double *cnt, long m)
{
    int v;
#if defined(__SSE2__)
    /* (the sums in a call of their own: nothing of theirs moved across the
       control register's writes) */
    unsigned int csr = _mm_getcsr();
    _mm_setcsr(csr | PCX_FLUSH_BITS);
    v = pcx_resid_flushed(Ap, Ai, Ax, tol, b, x, y, t, cnt, m);
    _mm_setcsr((_mm_getcsr() & ~PCX_FLUSH_BITS) | (csr & PCX_FLUSH_BITS));
    /* (its band the wider: where it decides, the gradual one decides alike;
       where it does not, the gradual one's own verdict) */
    if (v != PCX_UNDECIDED) return v;
#endif
    v = pcx_resid_gradual(Ap, Ai, Ax, tol, b, x, y, t, cnt, m);
    if (v != PCX_UNDECIDED) return v;
    return pcx_resid_last(Ap, Ai, Ax, tol, b, x, y, t, cnt, m);
}

/* a NaN whose quiet bit is clear: an operation on it raises invalid */
static inline int pcx_snan(double v)
{
    uint64_t u;
    memcpy(&u, &v, sizeof u);
    return isnan(v) && !(u & UINT64_C(0x0008000000000000));
}

/* (N): PCX_FAILS where provable, else PCX_UNDECIDED -- called with a value
   not finite in x, b or A.  `t` (m) scratch */
static int pcx_resid_nonfinite(const int *Ap, const int *Ai, const double *Ax, const double *b,
                               const double *x, double *t, long m)
{
    long i, j, p;
    for (i = 0; i < m; i++) {
        if (isinf(x[i]) || pcx_snan(x[i]) || pcx_snan(b[i])) return PCX_UNDECIDED;
        t[i] = 0.0;
    }
    for (j = 0; j < m; j++)
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            double a = Ax[p];
            if (isinf(a) || pcx_snan(a)) return PCX_UNDECIDED;
            if (!isnan(a) && !isnan(x[j])) t[Ai[p]] += fabs(a) * fabs(x[j]);
        }
    for (i = 0; i < m; i++)
        if (!(t[i] + (isfinite(b[i]) ? fabs(b[i]) : 0.0) < 1e300)) return PCX_UNDECIDED;
    return PCX_FAILS;
}

#if LDBL_MANT_DIG >= 64
/* a row's long double sum kept as two doubles, its value and remainder
   (the remainder exact but where subnormal: within 2^-1075), plus `v` */
static inline void pcx_wide_add(double *hi, double *lo, long double v)
{
    long double acc = ((long double) *hi + (long double) *lo) + v;
    double h = (double) acc;
    *hi = h;
    *lo = (double) (acc - (long double) h);
}

/* (W) on finite values.  `y`, `t` (m each): each row's sum, `cnt` (m) its
   entries, then its E */
static int pcx_resid_wide(const int *Ap, const int *Ai, const double *Ax, double tol,
                          const double *b, const double *x, double *y, double *t, double *cnt,
                          long m)
{
    long i, j, p;
    double max_b = 0.0, hi = 0.0, lo = 0.0, sc;
    for (i = 0; i < m; i++) y[i] = t[i] = cnt[i] = 0.0;
    for (j = 0; j < m; j++)
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            i = Ai[p];
            cnt[i] += 1.0;
            pcx_wide_add(&y[i], &t[i], (long double) fabs(Ax[p]) * (long double) fabs(x[j]));
        }
    for (i = 0; i < m; i++) {
        long double T = (long double) y[i] + (long double) t[i], k = cnt[i];
        double ab = fabs(b[i]);
        if (!(T + ab < (long double) 1e300)) return PCX_UNDECIDED;
        cnt[i] = (double) (k * (0.5L * DBL_EPSILON) * (1.0L + 0x1p-6L) * T + 0x1p-60L * ab
                           + (2.0L * k + 4.0L) * (long double) DBL_TRUE_MIN);
        y[i] = t[i] = 0.0;
    }
    for (j = 0; j < m; j++)
        for (p = Ap[j]; p < Ap[j + 1]; p++)
            pcx_wide_add(&y[Ai[p]], &t[Ai[p]], (long double) Ax[p] * (long double) x[j]);
    for (i = 0; i < m; i++) {
        long double S = (long double) y[i] + (long double) t[i], E = cnt[i];
        double ab = fabs(b[i]), u, v, au, av, h, l;
        u = (double) ((S + E) - (long double) b[i]);
        v = (double) ((S - E) - (long double) b[i]);
        au = fabs(u);
        av = fabs(v);
        h = au > av ? au : av;
        l = (u > 0.0 && v > 0.0) || (u < 0.0 && v < 0.0) ? (au < av ? au : av) : 0.0;
        if (h > hi) hi = h;
        if (l > lo) lo = l;
        if (ab > max_b) max_b = ab;
    }
    sc = pcx_py_max(max_b, 1e-300);
    if (hi / sc <= tol) return PCX_HOLDS;
    if (!(lo / sc <= tol)) return PCX_FAILS;
    return PCX_UNDECIDED;
}
#endif

/* where the passes above leave it undecided: (N) where a value is not
   finite, else (W) */
static int pcx_resid_last(const int *Ap, const int *Ai, const double *Ax, double tol,
                          const double *b, const double *x, double *y, double *t, double *cnt,
                          long m)
{
    long i, p;
    for (i = 0; i < m; i++)
        if (!isfinite(b[i]) || !isfinite(x[i])) return pcx_resid_nonfinite(Ap, Ai, Ax, b, x, t, m);
    for (p = 0; p < Ap[m]; p++)
        if (!isfinite(Ax[p])) return pcx_resid_nonfinite(Ap, Ai, Ax, b, x, t, m);
#if LDBL_MANT_DIG >= 64
    return pcx_resid_wide(Ap, Ai, Ax, tol, b, x, y, t, cnt, m);
#else
    (void) tol, (void) y, (void) cnt;
    return PCX_UNDECIDED;
#endif
}

/* the verdict's sums with gradual underflow, the bound 4 (k_i + 2) (eps T_i
   + 2 DBL_TRUE_MIN) */
static int pcx_resid_gradual(const int *Ap, const int *Ai, const double *Ax, double tol,
                             const double *b, const double *x, double *y, double *t,
                             double *cnt, long m)
{
    long i, j, p;
    double max_b = 0.0, hi = 0.0, lo = 0.0, sc;
    for (i = 0; i < m; i++) {
        if (!isfinite(b[i]) || !isfinite(x[i])) return PCX_UNDECIDED;
        y[i] = 0.0;
        t[i] = 0.0;
        cnt[i] = 0.0;
    }
    for (j = 0; j < m; j++)
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            double a = Ax[p];
            i = Ai[p];
            y[i] += a * x[j];
            t[i] += fabs(a) * fabs(x[j]);
            cnt[i] += 1.0;
        }
    for (i = 0; i < m; i++) {
        double ab = fabs(b[i]), e, u, v, au, av, h, l;
        if (!(t[i] + ab < 1e300)) return PCX_UNDECIDED;
        e = 4.0 * (cnt[i] + 2.0) * (DBL_EPSILON * t[i] + 2.0 * DBL_TRUE_MIN);
        u = (y[i] + e) - b[i];
        v = (y[i] - e) - b[i];
        au = fabs(u);
        av = fabs(v);
        h = au > av ? au : av;
        l = (u > 0.0 && v > 0.0) || (u < 0.0 && v < 0.0) ? (au < av ? au : av) : 0.0;
        if (h > hi) hi = h;
        if (l > lo) lo = l;
        if (ab > max_b) max_b = ab;
    }
    sc = pcx_py_max(max_b, 1e-300);
    if (hi / sc <= tol) return PCX_HOLDS;
    if (!(lo / sc <= tol)) return PCX_FAILS;
    return PCX_UNDECIDED;
}

__attribute__((noinline)) static int pcx_resid_flushed(const int *Ap, const int *Ai,
                                                        const double *Ax, double tol,
                                                        const double *b, const double *x,
                                                        double *y, double *t, double *cnt, long m)
{
    long i, j, p;
    double max_b = 0.0, hi = 0.0, lo = 0.0, sc, amax = 0.0, xmax = 0.0, w;
    for (i = 0; i < m; i++) {
        if (!isfinite(b[i]) || !isfinite(x[i])) return PCX_UNDECIDED;
        y[i] = 0.0;
        t[i] = 0.0;
        cnt[i] = 0.0;
        if (fabs(x[i]) > xmax) xmax = fabs(x[i]);
    }
    for (j = 0; j < m; j++)
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            double a = Ax[p];
            i = Ai[p];
            y[i] += a * x[j];
            t[i] += fabs(a) * fabs(x[j]);
            cnt[i] += 1.0;
            if (fabs(a) > amax) amax = fabs(a);
        }
    /* (a flushed operand's, product's or sum's part of the bound) */
    if (!(amax < 1e300)) return PCX_UNDECIDED;
    w = (amax + xmax + 4.0) * DBL_MIN;
    for (i = 0; i < m; i++) {
        double ab = fabs(b[i]), e, u, v, au, av, h, l;
        if (!(t[i] + ab < 1e300)) return PCX_UNDECIDED;
        e = 4.0 * (cnt[i] + 2.0) * (DBL_EPSILON * t[i]) + (cnt[i] + 2.0) * w + 8.0 * DBL_MIN;
        u = (y[i] + e) - b[i];
        v = (y[i] - e) - b[i];
        au = fabs(u);
        av = fabs(v);
        h = au > av ? au : av;
        l = (u > 0.0 && v > 0.0) || (u < 0.0 && v < 0.0) ? (au < av ? au : av) : 0.0;
        if (h > hi) hi = h;
        if (l > lo) lo = l;
        if (ab > max_b) max_b = ab;
    }
    sc = pcx_py_max(max_b, 1e-300);
    if (hi / sc <= tol) return PCX_HOLDS;
    if (!(lo / sc <= tol)) return PCX_FAILS;
    return PCX_UNDECIDED;
}

/* THE CSC KERNEL (THE PACKED LAYOUT, `_packed`: the reference's order at
   size): `y = A x` (`absv` 0) or `|A| x` (`absv` 1) on CSC arrays -- each
   row from 0.0, each product rounded, then added column by column
   (`_packed.PatternMatrix.matvec`'s np.bincount; on KLU's arrays a matrix
   zero at KLU's entries off its own pattern, which add +0.0) */
static inline void pcx_csc_mv(const int *Ap, const int *Ai, const double *Ax, const double *x,
                              double *y, long m, int absv)
{
    long i, j, p;
    for (i = 0; i < m; i++) y[i] = 0.0;
    if (absv)
        for (j = 0; j < m; j++)
            for (p = Ap[j]; p < Ap[j + 1]; p++) y[Ai[p]] += fabs(Ax[p]) * x[j];
    else
        for (j = 0; j < m; j++)
            for (p = Ap[j]; p < Ap[j + 1]; p++) y[Ai[p]] += Ax[p] * x[j];
}

/* nrsolver's residual test, `aF < reltol * (dot(abs(J), ax) + aF) +
   abstol` row by row (all rows: `alltrue`), on J's CSC (|J| zero off the
   pattern): every term non-negative, so T_i = S_i and the dgemv's row lies
   in [max(S_i - e_i, 0), S_i + e_i]; the row's test, monotone in it, holds
   where it holds at the lower end, fails where it fails at the upper (one
   row failing decides the whole).  `F` the residual, `ax` |x_next| (both
   finite: checked before), `t`, `cnt` (m each) scratch */
static int pcx_scale_sparse(const int *Ap, const int *Ai, const double *Ax, const double *ax,
                            const double *F, double reltol, const double *abstol, double *t,
                            double *cnt, long m)
{
    long i, j, p;
    int undecided = 0;
    for (i = 0; i < m; i++) {
        t[i] = 0.0;
        cnt[i] = 0.0;
    }
    for (j = 0; j < m; j++)
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            i = Ai[p];
            t[i] += fabs(Ax[p]) * ax[j];
            cnt[i] += 1.0;
        }
    for (i = 0; i < m; i++) {
        double aF = fabs(F[i]), e, lo, hi;
        if (!(t[i] + aF < 1e300)) return PCX_UNDECIDED;
        e = 4.0 * (cnt[i] + 2.0) * (DBL_EPSILON * t[i] + 2.0 * DBL_TRUE_MIN);
        lo = t[i] - e;
        if (lo < 0.0) lo = 0.0;
        hi = t[i] + e;
        if (!(aF < reltol * (hi + aF) + abstol[i])) return PCX_FAILS;
        if (!(aF < reltol * (lo + aF) + abstol[i])) undecided = 1;
    }
    return undecided ? PCX_UNDECIDED : PCX_HOLDS;
}

/* THE REFACTOR KEPT ON ITS GROWTH (`StructuralKLU._growth_keeps`): whether
   a refactored numeric `num` whose residual test failed is kept -- KLU's
   reciprocal pivot growth of it (`klu_rgrowth`, pcx_rgrowth_fn, read from the
   common at `rg_off`) at least `rg_keep` -- only where its values `Ax`
   (`Ap[m]` of them) and the solution `x` (m) are all finite.  No entry, no
   offset (-1), no numeric, a failed call or a NaN: not kept.  What the call
   raised is cleared, as the Python's ctypes call never sees it */
typedef int (*pcx_rgrowth_fn)(int *, int *, double *, void *, void *, void *);

static inline int pcx_rgrowth_keeps(void *rgrowth, long rg_off, double rg_keep, int *Ap, int *Ai,
                                    double *Ax, void *sym, void *num, void *common,
                                    const double *x, long m)
{
    double g;
    long i;
    int ok;
    if (!rgrowth || rg_off < 0 || !num) return 0;
    for (i = 0; i < Ap[m]; i++) if (!isfinite(Ax[i])) return 0;
    for (i = 0; i < m; i++) if (!isfinite(x[i])) return 0;
    ok = ((pcx_rgrowth_fn) rgrowth)(Ap, Ai, Ax, sym, num, common);
    feclearexcept(FE_ALL_EXCEPT);
    if (!ok) return 0;
    memcpy(&g, (char *) common + rg_off, sizeof g);
    return g >= rg_keep;
}
