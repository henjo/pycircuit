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
 * covers that and the rounding of S_i -+ e_i.  A test monotone in the row's
 * product then holds where it holds at the bound's worse end and fails where
 * it fails at its better end; elsewhere -- and where a value is not finite
 * or a row's sum reaches 1e300 (there numpy's flags could be raised) -- the
 * dense test decides. */
#include <float.h>
#include <math.h>

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
static int pcx_resid_sparse(const int *Ap, const int *Ai, const double *Ax, double tol,
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
    if (hi / sc <= tol) return PCX_HOLDS;
    if (!(lo / sc <= tol)) return PCX_FAILS;
    return PCX_UNDECIDED;
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
        e = 4.0 * (cnt[i] + 2.0) * DBL_EPSILON * t[i];
        lo = t[i] - e;
        if (lo < 0.0) lo = 0.0;
        hi = t[i] + e;
        if (!(aF < reltol * (hi + aF) + abstol[i])) return PCX_FAILS;
        if (!(aF < reltol * (lo + aF) + abstol[i])) undecided = 1;
    }
    return undecided ? PCX_UNDECIDED : PCX_HOLDS;
}
