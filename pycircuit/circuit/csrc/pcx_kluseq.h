/* pcx_kluseq.h -- StructuralKLU.solve's own sequence (`linearsolver`),
 * written once for the units that repeat it (dcnewton.c, pcnrdc.c): a
 * refactor on the numeric it holds, a fresh factor on a zero pivot,
 * klu_solve, and after a refactor the residual `max |A x - b| / max(max |b|,
 * 1e-300)` against REFACTOR_RESIDUAL_TOL -- a fresh factor where it fails.
 * The struct and THE NUMERIC: pcx_klu.h (included before this).
 *
 * A step returns 0, or what stopped it: PCX_KLU_FACTOR (a fresh factor
 * failed: `_factor`'s LinAlgError('Singular matrix'), `klu_out` NULL),
 * PCX_KLU_SOLVE (klu_solve failed: its LinAlgError('klu_solve failed')),
 * PCX_KLU_FLAGS (a flag numpy reports in the residual's arithmetic).
 *
 * The residual's test needs only its verdict: on KLU's own CSC (`a_built`
 * 0: the dense A not made) a sparse sum with its rounding bound decides it
 * where the bound does, the dense test where it does not
 * (`pcx_klu_resid_sparse`). */
#include <float.h>
#include <math.h>
#include <fenv.h>

/* the flags numpy reports: underflow is not one (numpy ignores it by
   default; a unit's glue declines where it does not) */
#define PCX_NP_FLAGS (FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

enum pcx_klu_stop { PCX_KLU_FACTOR = 1, PCX_KLU_SOLVE = 2, PCX_KLU_FLAGS = 3 };

/* a test's verdict where a rounding bound decides it (`pcx_klu_resid_sparse`,
   pcnrdc.c's `cf_sparse`): it holds, it fails, or the dense test decides */
enum pcx_verdict { PCX_HOLDS = 0, PCX_FAILS = 1, PCX_UNDECIDED = 2 };

static void pcx_klu_free(pcx_klu_t *k, void *num)
{
    if (num) ((klu_free_numeric_t) k->klu_free)(&num, k->klu_common);
}

/* a fresh factor, as StructuralKLU's `_factor` after freeing its numeric
   (the entry's kept: `pcx_klu_finish`) */
static long pcx_klu_fresh(pcx_klu_t *k)
{
    if (k->klu_out != k->klu_num) pcx_klu_free(k, k->klu_out);
    k->klu_out = ((klu_factor_t) k->klu_factor)(k->kAp, k->kAi, k->kAx, k->klu_sym,
                                                k->klu_common);
    if (!k->klu_out) return PCX_KLU_FACTOR;
    k->factors++;
    return 0;
}

/* klu_solve of `b` into `x` */
static long pcx_klu_solve_b(pcx_klu_t *k, const double *b, double *x, long m)
{
    long i;
    for (i = 0; i < m; i++) x[i] = b[i];
    if (!((klu_solve_t) k->klu_solve)(k->klu_sym, k->klu_out, (int) m, 1, x, k->klu_common))
        return PCX_KLU_SOLVE;
    return 0;
}

/* the residual's test as the Python makes it on the dense A: `arr @ x`
   (numpy's matmul: this dgemv), `- rhs`, |.|, the maxima -- numpy's max
   propagates a NaN: so does `nan_r` -- against the tolerance.  1 it fails (a
   fresh factor), 0 it holds, -1 a flag numpy reports */
static int pcx_klu_resid_dense(pcx_klu_t *k, const double *A, const double *b,
                               const double *x, double *y, long m, void *dgemv)
{
    long i;
    double max_b = 0.0, max_r = 0.0;
    int nan_r = 0;
    ((dgemv64_t) dgemv)(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, m, m, 1.0, A, m, x, 1,
                        0.0, y, 1);
    for (i = 0; i < m; i++) y[i] = y[i] - b[i];
    if (fetestexcept(PCX_NP_FLAGS)) return -1;
    for (i = 0; i < m; i++) {
        double a = fabs(b[i]), d = fabs(y[i]);
        if (a > max_b) max_b = a;
        if (isnan(d)) nan_r = 1;
        else if (d > max_r) max_r = d;
    }
    return nan_r || !(max_r / pcx_py_max(max_b, 1e-300) <= k->resid_tol);
}

/* the same verdict on KLU's CSC (A zero off its pattern).  A sum of terms
   rounds only where two operands holding a nonzero term meet -- a zero
   product joins exactly, fused or not -- so a row with k_i entries on the
   pattern rounds at most k_i times on any path of any summation tree: its
   S_i = sum A_ij x_j, in the pattern's order, is within 2 gamma_(k_i + 1)
   T_i of the dgemv's (T_i = sum |A_ij| |x_j|), and with e_i = 4 (k_i + 2)
   eps T_i the row's `fl(dgemv - b)` lies between `(S_i - e_i) - b_i` and
   `(S_i + e_i) - b_i` (rounding is monotone), max |.| between the rows'
   lower and upper ends, and the test
   -- monotone in it, the scale exact -- holds where it holds at the upper
   end, fails where it fails at the lower.  PCX_UNDECIDED also where a value
   is not finite or a row's T_i + |b_i| reaches 1e300 (there numpy's flags
   could be raised).  Where it decides, no flag numpy reports is raised here
   or in the dense test */
static int pcx_klu_resid_sparse(pcx_klu_t *k, const double *b, const double *x,
                                double *y, double *t, long m)
{
    long i, j, p;
    double max_b = 0.0, hi = 0.0, lo = 0.0, sc, *cnt = t + m;
    for (i = 0; i < m; i++) {
        if (!isfinite(b[i]) || !isfinite(x[i])) return PCX_UNDECIDED;
        y[i] = 0.0;
        t[i] = 0.0;
        cnt[i] = 0.0;
    }
    for (j = 0; j < m; j++)
        for (p = k->kAp[j]; p < k->kAp[j + 1]; p++) {
            double a = k->kAx[p];
            i = k->kAi[p];
            y[i] += a * x[j];
            t[i] += fabs(a) * fabs(x[j]);
            cnt[i] += 1.0;
        }
    for (i = 0; i < m; i++) {
        double ab = fabs(b[i]), e, u, v, au, av, h, l;
        if (!(t[i] + ab < 1e300)) return PCX_UNDECIDED;
        e = 4.0 * (cnt[i] + 2.0) * DBL_EPSILON * t[i];
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
    if (hi / sc <= k->resid_tol) return PCX_HOLDS;
    if (!(lo / sc <= k->resid_tol)) return PCX_FAILS;
    return PCX_UNDECIDED;
}

/* StructuralKLU.solve(A, b) on its own sequence, `kAx` A's values in its
   order (gathered by the caller) and `A` the same matrix dense (m x m,
   row-major: the residual's `arr @ x`, numpy's matmul -- its dgemv) where
   `a_built`, else a buffer the residual's dense test scatters A into from
   the CSC where the sparse one does not decide (counted in `dense`, which
   may be NULL).  The solution in `x`, `y` (m) and `t` (2 m, read only
   where not `a_built`) scratch */
static long pcx_klu_step(pcx_klu_t *k, double *A, const double *b, double *x, double *y,
                         double *t, long m, void *dgemv, int a_built, long *dense)
{
    long i, j, p, st, reused = k->klu_out != NULL;
    int fails = PCX_UNDECIDED;
    if (reused) {
        int ok = ((klu_refactor_t) k->klu_refactor)(k->kAp, k->kAi, k->kAx, k->klu_sym,
                                                  k->klu_out, k->klu_common);
        k->refactors++;
        if (!ok) {
            k->pivot_fallbacks++;
            if ((st = pcx_klu_fresh(k)) != 0) return st;
            reused = 0;
        }
    } else if ((st = pcx_klu_fresh(k)) != 0) {
        return st;
    }
    if ((st = pcx_klu_solve_b(k, b, x, m)) != 0) return st;
    /* (what KLU raised inside: the Python's ctypes calls never see it) */
    feclearexcept(FE_ALL_EXCEPT);
    if (!reused) return 0;
    if (!a_built) {
        fails = pcx_klu_resid_sparse(k, b, x, y, t, m);
        if (fails == PCX_UNDECIDED) {
            /* (the dense test on A as the Python holds it: zero off the
               pattern; what the sparse sums raised is not the Python's) */
            feclearexcept(FE_ALL_EXCEPT);
            for (i = 0; i < m * m; i++) A[i] = 0.0;
            for (j = 0; j < m; j++)
                for (p = k->kAp[j]; p < k->kAp[j + 1]; p++) A[k->kAi[p] * m + j] = k->kAx[p];
            if (dense) (*dense)++;
        }
    }
    if (fails == PCX_UNDECIDED) fails = pcx_klu_resid_dense(k, A, b, x, y, m, dgemv);
    if (fails < 0) return PCX_KLU_FLAGS;
    if (fails) {
        k->residual_fallbacks++;
        if ((st = pcx_klu_fresh(k)) != 0) return st;
        if ((st = pcx_klu_solve_b(k, b, x, m)) != 0) return st;
        feclearexcept(FE_ALL_EXCEPT);
    }
    return 0;
}

/* the numeric as the outcome leaves it (THE NUMERIC, pcx_klu.h): a
   hand-back the entry's (any other freed), a served outcome the
   sequence's (the entry's freed where another replaced it) */
static void pcx_klu_finish(pcx_klu_t *k, int handed_back)
{
    if (handed_back) {
        if (k->klu_out != k->klu_num) pcx_klu_free(k, k->klu_out);
        k->klu_out = k->klu_num;
    } else if (k->klu_out != k->klu_num) {
        pcx_klu_free(k, k->klu_num);
    }
}
