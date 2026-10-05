"""RADAU'S TRANSFORM NEWTON IN C (speed round 10, B3.3; 2026-10-05).

The cost transform (`_RadauStages._rk_step_transformed`: PSS's default on
PSP-class circuits under 'auto') is a SIMPLIFIED Newton -- `C` and `G`
frozen at ``x_n``, the coupled solve through ``A^{-1}``'s eigenbasis, its
factors made once a step (`_radau_frozen`, B3.1).  Its loop spent ~3.5 M
instructions a step on the PSP stage's radau PSS in Python and numpy
around ~0.4 M of device kernels: per iteration the three stages' passes
through `_tran_core.passes`, the residual (`_coupled_stage_system`), the
frozen solve (`_FrozenTransform.solve`: the right-hand sides, numpy's kept
LU, the complex factor's KLU solve and its residual check, the update),
and per stage the step, the limiter and the convergence test.  Here ONE C
call runs the loop: per iteration the stages' ``q`` and ``i`` through the
evaluate core's own C (formula 4, as `passes`), ``K_j = -(i_j + u_j)``,
``F_i = (q_i - q_n) - h (((0 + A_i0 K_0) + A_i1 K_1) + A_i2 K_2)`` reduced
(`_tran_radau_c`'s assembly, without the blocks); the right-hand sides
``-((P_k0 F0 + P_k1 F1) + P_k2 F2)`` -- a complex scalar times a real
vector, as numpy multiplies it: the vector cast to complex, numpy's
complex product; numpy's kept LU (`linearsolver.NumpyLU`'s protocol: its
OpenBLAS's `dgesv` at the first solve, `dgetrs` on the kept factors
after); the complex factor's KLU on the solver's own handles (the refactor
`solve_prepared` would make, at most once a record, then `klu_z_solve`);
``dY_i = V_i0 w0 + 2 Re(V_i1 w1)``; per stage ``Y_prev + insert(dY_i,
iref, 0.0)``, the limiter walk's own driver and the limited step's largest
reduced entry; the convergence test after the update, as the Python tests.
Bit for bit the Python loop (pinned in `tests/test_pinned_pairs.py`).

⚠ NUMPY'S COMPLEX PRODUCT IS READ, NOT ASSUMED.  On this box (numpy's
X86_V3 loops: AVX2 and FMA3) it FUSES -- ``re = fma(ar, br, -(ai*bi))``,
``im = fma(ar, bi, ai*br)``, at every position of every length -- and a
build or CPU without FMA multiplies separately.  `cmul_mode` asks numpy
itself once, on inputs where the two forms round apart (and on signed
zeros and underflow, where only the fused form's sign rules hold), and
the C takes the form numpy showed; neither: the path is off.

⚠ THE RESIDUAL CHECK IS DECIDED WITH A MARGIN.  `solve_prepared` falls
back to a fresh factor where ``max|A x - b| / max|b|`` exceeds its
tolerance (1e-8), read with numpy's complex `abs` (its own algorithm, not
`hypot`) of SciPy's product.  The C computes the same quantities in its own
rounding and decides only outside the band their difference can reach
(``32 (m + 4) eps`` of ``max_k sum_j |A_kj| |x_j|`` over ``max|b|``, plus
``16 eps`` of the residual): a residual nearer the tolerance than that, or
past it, hands the iteration back and the Python decides, as it does --
in practice the residual sits ~1e-15 below a 1e-8 tolerance.

TRANSACTIONAL.  The C writes nothing outside its buffers but the two
factorisations' state, which the Python would have left the same: numpy's
kept LU after its first `dgesv` (`NumpyLU._info`), and the KLU numeric
holding a refactor of the step's record (`ComplexKLUSolver._fresh`,
`refactors`).  A refactor of the same values with the same pivots repeats
its bits and `dgetrs` on `dgesv`'s factors is `dgesv`'s solve, so a Python
loop run after a hand-back makes the C's iterates again.  The step's
source memo holds the three stage times (their calls made here, before
the C) and its hit count for every later iteration's.  Anything unusual --
a non-finite value, a floating-point exception in its own arithmetic
(numpy would warn), a non-finite real right-hand side (the analysis
solver's), a singular real factor (numpy's error), the KLU refactor or
solve failing, the residual check not decided or failing (the fallback's
fresh factor), the walk stopping, maxiter -- hands the loop back: the
source memo's entries and counts rolled back, and the Python loop runs
from the same seed (the same answer, exception, warning).

DECLINES, before anything is touched (`_paths.no`): the switch (and the
evaluate core's two), no frozen factors or no kept LU or marshalled record
(`_radau_frozen`'s own), no build, numpy's complex product of neither
form, a toolkit that is not numeric, a subclass of the transient or a
patched piece of the loop, an instance shadow of the circuit's passes or
limiter or of the transient's pieces, a `provided_function`, a stateful
limiter or circuit-level limiting, a bypass, no source memo, a circuit the
core cannot serve, a limiting element that is not its C kernel, a batch
not ready, a KLU solver whose numeric is not of the record's pattern (the
Python analyses and factors), a singular kept LU, inputs of another type
or size.

`ENABLED` (env `PYCIRCUIT_RADAU_TC=0`) is read on every call; `STATUS` says
why the path is off where it is.  History: `doc/transient_history.md`,
`_tran_radau_tc.py`.
"""
import ctypes
import math
import os

import numpy as np

from pycircuit.circuit import _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_RADAU_TC', '1') != '0'
STATUS = 'not loaded'

RADAU_TC = r"""
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
typedef long (*core_fn_t)(const void *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J);
typedef long (*walk_fn_t)(void **F, double **PR, const int64_t *TI, double T,
            const int64_t *NM, const int64_t *OFF, const int64_t *K,
            double *xg, const double *x0g, long start, long n);
typedef void (*dgesv_t)(const int64_t *, const int64_t *, double *, const int64_t *,
                        int64_t *, double *, const int64_t *, int64_t *);
typedef void (*dgetrs_t)(const char *, const int64_t *, const int64_t *, const double *,
                         const int64_t *, const int64_t *, double *, const int64_t *,
                         int64_t *, size_t);
typedef int (*kluref_t)(int *, int *, double *, void *, void *, void *);
typedef int (*klusol_t)(void *, void *, int, int, double *, void *);
typedef struct {
    const void *core; void *core_fn; double T;
    const double *bini; double *dummy; long bits;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv; void *dgetrs; void *klu_refactor; void *klu_solve;
    void *Symbolic; void *Numeric; void *Common;
    int *Ap; int *Ai; double *Ax;
    long n, iref, maxiter, fused;
    double h, reltol, abstol, rtol;
    double A[9];
    double p[12];
    double v[9];
    const double *qn, *u, *seed;
    double *lu_a; int64_t *lu_ipiv; int64_t lu_info;
    long klu_fresh, refactored;
    double *Y, *Yp, *q, *iv, *K, *F, *R, *b0, *x0, *r1, *x1, *y1, *ab, *dY, *stp, *tolv;
    long iters, status;
} radau_tf_t;

#define TC_FLAGS (FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO)

static long tf_bail(radau_tf_t *s, long why)
{
    feclearexcept(FE_ALL_EXCEPT);
    s->status = why;
    return why;
}

/* numpy's complex product (ar + i ai)(br + i bi), in the form numpy's loops
   take on this machine (`cmul_mode`): fused, or two products and a sum */
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

long hdl_fn(radau_tf_t *s)
{
    long n = s->n, m = n - 1, r = s->iref, it, i, j, k, e;
    core_fn_t core = (core_fn_t) s->core_fn;
    walk_fn_t walk = (walk_fn_t) s->walk_fn;
    int64_t m64 = m, one64 = 1, info64 = 0;
    const double h = s->h;
    s->iters = 0;
    s->status = 0;
    s->refactored = 0;
    /* every input finite: anything else is the Python loop's */
    if (!(isfinite(h) && isfinite(s->reltol) && isfinite(s->abstol))) return tf_bail(s, 2);
    for (k = 0; k < 9; k++) if (!isfinite(s->A[k]) || !isfinite(s->v[k])) return tf_bail(s, 2);
    for (k = 0; k < 12; k++) if (!isfinite(s->p[k])) return tf_bail(s, 2);
    for (k = 0; k < n; k++) if (!isfinite(s->qn[k])) return tf_bail(s, 2);
    for (k = 0; k < 3 * n; k++) if (!isfinite(s->u[k])) return tf_bail(s, 2);
    memcpy(s->Y, s->seed, (size_t) (3 * n) * sizeof(double));
    for (it = 0; it < s->maxiter; it++) {
        double scale = 0.0, ynorm = 0.0;
        /* the stages' q and i (formula 4: i in the core's bini, q in place;
           a non-finite stage: `passes` declines, the Python evaluates) */
        for (j = 0; j < 3; j++) {
            for (k = 0; k < n; k++) if (!isfinite(s->Y[j * n + k])) return tf_bail(s, 2);
            core(s->core, s->Y + j * n, s->T, s->bits, 4, 0.0, 0.0, 0.0, 1.0, 0.0,
                 s->dummy, s->dummy, s->dummy, s->dummy, s->dummy, s->dummy, s->q + j * n,
                 s->dummy, s->dummy, s->dummy, s->dummy);
            memcpy(s->iv + j * n, s->bini, (size_t) n * sizeof(double));
        }
        feclearexcept(FE_ALL_EXCEPT);
        /* `_coupled_stage_system`'s residual: K_j = -(i_j + u_j); Python's
           `sum` starts at the int 0 (a -0.0 becomes +0.0: `0.0 +`) */
        for (j = 0; j < 3; j++)
            for (k = 0; k < n; k++) s->K[j * n + k] = -(s->iv[j * n + k] + s->u[j * n + k]);
        for (i = 0; i < 3; i++) {
            const double a0 = s->A[3 * i], a1 = s->A[3 * i + 1], a2 = s->A[3 * i + 2];
            double *f = s->F + i * n;
            long rr;
            for (k = 0; k < n; k++) {
                double sm = 0.0 + a0 * s->K[k];
                sm = sm + a1 * s->K[n + k];
                sm = sm + a2 * s->K[2 * n + k];
                f[k] = (s->q[i * n + k] - s->qn[k]) - h * sm;
            }
            for (k = 0, rr = 0; k < n; k++) {
                if (k == r) continue;
                s->R[i * m + rr] = f[k];
                rr++;
            }
        }
        if (fetestexcept(TC_FLAGS)) return tf_bail(s, 3);
        /* `_transform_rhs`: -((P_k0 F0 + P_k1 F1) + P_k2 F2), each F cast to
           complex (F + 0j) and multiplied as numpy multiplies; rhs0's
           imaginary part made as numpy makes it (its flags), then dropped */
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
        if (fetestexcept(TC_FLAGS)) return tf_bail(s, 3);
        /* a non-finite real right-hand side is the analysis solver's */
        for (k = 0; k < m; k++) if (!isfinite(s->b0[k])) return tf_bail(s, 4);
        /* numpy's kept LU: `dgesv` at its first solve, `dgetrs` after */
        memcpy(s->x0, s->b0, (size_t) m * sizeof(double));
        if (s->lu_info < 0) {
            ((dgesv_t) s->dgesv)(&m64, &one64, s->lu_a, &m64, s->lu_ipiv, s->x0, &m64, &info64);
            s->lu_info = info64;
        } else if (s->lu_info == 0) {
            ((dgetrs_t) s->dgetrs)("N", &m64, &one64, s->lu_a, &m64, s->lu_ipiv, s->x0, &m64,
                                   &info64, (size_t) 1);
        }
        feclearexcept(FE_ALL_EXCEPT);
        if (s->lu_info != 0) return tf_bail(s, 5);
        for (k = 0; k < m; k++) if (!isfinite(s->x0[k])) return tf_bail(s, 9);
        /* the complex factor: its refactor at most once a record, then the solve */
        if (!s->klu_fresh) {
            if (!((kluref_t) s->klu_refactor)(s->Ap, s->Ai, s->Ax, s->Symbolic, s->Numeric,
                                               s->Common))
                return tf_bail(s, 6);
            s->klu_fresh = 1;
            s->refactored++;
        }
        memcpy(s->x1, s->r1, (size_t) (2 * m) * sizeof(double));
        if (!((klusol_t) s->klu_solve)(s->Symbolic, s->Numeric, (int) m, 1, s->x1, s->Common))
            return tf_bail(s, 7);
        feclearexcept(FE_ALL_EXCEPT);
        /* the residual check, decided outside the band numpy's `abs` and
           SciPy's product can round apart from these */
        {
            double sc = 0.0, res = 0.0, bnd = 0.0, band;
            for (k = 0; k < m; k++) {
                double a = hypot(s->r1[2 * k], s->r1[2 * k + 1]);
                if (!(a <= DBL_MAX)) return tf_bail(s, 8);
                if (a > sc) sc = a;
            }
            if (!(sc >= 1e-300)) sc = 1e-300;
            memset(s->y1, 0, (size_t) (2 * m) * sizeof(double));
            memset(s->ab, 0, (size_t) m * sizeof(double));
            for (j = 0; j < m; j++) {
                const double xr = s->x1[2 * j], xi = s->x1[2 * j + 1], xa = hypot(xr, xi);
                for (k = s->Ap[j]; k < s->Ap[j + 1]; k++) {
                    const long row = s->Ai[k];
                    const double ar = s->Ax[2 * k], ai = s->Ax[2 * k + 1];
                    s->y1[2 * row] += ar * xr - ai * xi;
                    s->y1[2 * row + 1] += ar * xi + ai * xr;
                    s->ab[row] += hypot(ar, ai) * xa;
                }
            }
            for (k = 0; k < m; k++) {
                const double a = hypot(s->y1[2 * k] - s->r1[2 * k],
                                       s->y1[2 * k + 1] - s->r1[2 * k + 1]);
                if (!(a <= DBL_MAX)) return tf_bail(s, 8);
                if (a > res) res = a;
                if (s->ab[k] > bnd) bnd = s->ab[k];
            }
            if (fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)) return tf_bail(s, 8);
            res = res / sc;
            band = (32.0 * (double) (m + 4)) * DBL_EPSILON * (bnd / sc) + 16.0 * DBL_EPSILON * res;
            if (!(res + band <= s->rtol)) return tf_bail(s, 8);
            feclearexcept(FE_ALL_EXCEPT);
        }
        /* `_transform_back`: dY_i = V_i0 w0 + 2 Re(V_i1 w1), the complex
           product whole (its flags), as numpy makes it */
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
        if (fetestexcept(TC_FLAGS)) return tf_bail(s, 3);
        /* per stage: Y_prev + insert(dY_i, iref, 0.0), the walk against
           Y_prev, the limited step's largest reduced entry */
        memcpy(s->Yp, s->Y, (size_t) (3 * n) * sizeof(double));
        for (i = 0; i < 3; i++) {
            const double *yp = s->Yp + i * n, *d = s->dY + i * m;
            double *yn = s->Y + i * n, mx = 0.0;
            for (k = 0; k < r; k++) yn[k] = yp[k] + d[k];
            yn[r] = yp[r] + 0.0;
            for (k = r + 1; k < n; k++) yn[k] = yp[k] + d[k - 1];
            if (fetestexcept(TC_FLAGS)) return tf_bail(s, 3);
            if (walk) {
                e = walk(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK, yn, yp, 0, s->wn);
                if (e < s->wn) return tf_bail(s, 10);
                feclearexcept(FE_ALL_EXCEPT);
            }
            for (k = 0; k < n; k++) s->stp[k] = yn[k] - yp[k];
            if (fetestexcept(TC_FLAGS)) return tf_bail(s, 3);
            for (k = 0; k < n; k++) {
                double vv;
                if (k == r) continue;
                vv = fabs(s->stp[k]);
                if (!isfinite(vv)) return tf_bail(s, 9);
                if (vv > mx) mx = vv;
            }
            if (mx > scale) scale = mx;
        }
        /* `_stages_converged`: the step within reltol of the largest
           reduced stage plus abstol */
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
        if (fetestexcept(TC_FLAGS)) return tf_bail(s, 3);
        s->iters = it + 1;
        if (scale <= s->tolv[1]) {
            s->status = 1;
            return 1;
        }
    }
    return tf_bail(s, 11);
}
"""
RADAU_TC_CDEF = """
typedef struct {
    const void *core; void *core_fn; double T;
    const double *bini; double *dummy; long bits;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv; void *dgetrs; void *klu_refactor; void *klu_solve;
    void *Symbolic; void *Numeric; void *Common;
    int *Ap; int *Ai; double *Ax;
    long n, iref, maxiter, fused;
    double h, reltol, abstol, rtol;
    double A[9];
    double p[12];
    double v[9];
    const double *qn, *u, *seed;
    double *lu_a; int64_t *lu_ipiv; int64_t lu_info;
    long klu_fresh, refactored;
    double *Y, *Yp, *q, *iv, *K, *F, *R, *b0, *x0, *r1, *x1, *y1, *ab, *dY, *stp, *tolv;
    long iters, status;
} radau_tf_t;
long hdl_fn(radau_tf_t *s);
"""

#: the C's status codes past 1 (converged): why it handed the loop back
BAIL = {2: 'nonfinite', 3: 'flags', 4: 'rhs', 5: 'lu', 6: 'refactor', 7: 'klusolve',
        8: 'residual', 9: 'step', 10: 'walkstop', 11: 'maxiter'}

_driver = None
_MOD = {}


def _cmul_samples():
    """Pairs ``(a, b)`` of complex operands on which numpy's two possible
    complex products round apart -- drawn where the separate products
    differ from the fused ones -- plus signed zeros and underflow (the
    right-hand sides' ``P F``: a real vector cast to complex), as
    deterministic arrays."""
    rng = np.random.default_rng(20261005)
    ar, ai, br, bi = (rng.standard_normal(4000) * 10.0 ** rng.integers(-20, 20, size=4000)
                      for _ in range(4))
    keep = np.array([math.fma(a, c, -(d * e)) != a * c - d * e
                     for a, c, d, e in zip(ar, br, ai, bi)])
    a = (ar + 1j * ai)[keep][:64]
    b = (br + 1j * bi)[keep][:64]
    ## the real right-hand sides' corner: a tiny exact product rounding to
    ## a signed zero, against a zero imaginary product of either sign
    tiny = np.array([3e-200, -3e-200, 0.0, -0.0, 1.5, -2.5])
    za = np.array([complex(x, y) for x in (1e-200, -1e-200, 2.0, -0.0)
                   for y in (1.0, -1.0, 0.0, -0.0, 1e-200, -1e-200)])
    return a, b, za, tiny


def cmul_mode():
    """1 where numpy's complex product is the fused form, 0 where it is the
    separate one, None where it is neither -- read from numpy itself on
    `_cmul_samples`: a complex scalar times complex vectors (the update)
    and times real vectors (the right-hand sides), every position of
    lengths 1 to 17."""
    if 'cmul' in _MOD:
        return _MOD['cmul']
    a, b, za, tiny = _cmul_samples()

    def forms(x, y):
        return ((math.fma(x.real, y.real, -(x.imag * y.imag)),
                 math.fma(x.real, y.imag, x.imag * y.real)),
                (x.real * y.real - x.imag * y.imag, x.real * y.imag + x.imag * y.real))
    seen = set()
    for L in range(1, 18):
        for j in range(0, len(a), 4):
            x = np.complex128(a[j])
            ys = np.resize(b, L)
            got = x * ys
            for t in range(L):
                f, sp = forms(complex(x), complex(ys[t]))
                g = (float(got[t].real), float(got[t].imag))
                seen.add('f' if _same(g, f) else ('s' if _same(g, sp) else 'x'))
    for x in za:
        xc = np.complex128(x)
        got = xc * tiny
        for t in range(len(tiny)):
            f, sp = forms(complex(x), complex(float(tiny[t]), 0.0))
            g = (float(got[t].real), float(got[t].imag))
            ok_f, ok_s = _same(g, f), _same(g, sp)
            seen.add('f' if ok_f and not ok_s else ('s' if ok_s and not ok_f else
                                                       ('b' if ok_f else 'x')))
    seen.discard('b')
    _MOD['cmul'] = 1 if seen == {'f'} else (0 if seen == {'s'} else None)
    return _MOD['cmul']


def _same(g, f):
    """Two (re, im) pairs the same bits (the sign of a zero, a NaN's
    payload)."""
    return (np.float64(g[0]).tobytes() == np.float64(f[0]).tobytes()
            and np.float64(g[1]).tobytes() == np.float64(f[1]).tobytes())


def driver():
    """`(ffi, cfn, dgesv, dgetrs)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        from pycircuit.circuit import _hdl_cbackend as cb
        from pycircuit.circuit import _tran_radau_c, linearsolver
        if _tran_radau_c.driver() is None:
            STATUS = f'off (the radau C: {_tran_radau_c.STATUS})'
            return None
        lap = linearsolver._numpy_lapack()
        if lap is None:
            STATUS = "off (numpy's own OpenBLAS not found)"
            return None
        mode = cmul_mode()
        if mode is None:
            STATUS = "off (numpy's complex product is neither form)"
            return None
        try:
            ffi, cfn, _key, _cold, _secs = cb.load_kernel(RADAU_TC, RADAU_TC_CDEF)
        except (cb.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _ct, gesv, getrs = lap
        from pycircuit.circuit import _tran_newton_c, _tran_radau_c
        _MOD['NC'], _MOD['RC'] = _tran_newton_c, _tran_radau_c
        _driver = (ffi, cfn, ctypes.cast(gesv, ctypes.c_void_p).value,
                   ctypes.cast(getrs, ctypes.c_void_p).value)
        STATUS = 'c'
    return _driver or None


def _chain():
    """The loop's pieces the C stands in for, as defined (`_paths.genuine`),
    read once every piece is its module's own; None: not yet.  A call
    compares what the Python loop would call now (`_now`)."""
    if 'chain' not in _MOD:
        from pycircuit.circuit import _tran_radau, _tran_radau_c
        from pycircuit.circuit.linearsolver import ComplexKLUSolver, NumpyLU
        base = _tran_radau_c._MOD.get('chain') or _tran_radau_c._chain()
        if base is None:
            return None
        R, L = 'pycircuit.circuit._tran_radau', 'pycircuit.circuit.linearsolver'
        own = ((_tran_radau._RadauStages._transform_loop, '_RadauStages._transform_loop', R),
               (_tran_radau._FrozenTransform.solve, '_FrozenTransform.solve', R),
               (_tran_radau._transform_rhs, '_transform_rhs', R),
               (_tran_radau._transform_back, '_transform_back', R),
               (NumpyLU.solve, 'NumpyLU.solve', L),
               (ComplexKLUSolver.solve_prepared, 'ComplexKLUSolver.solve_prepared', L))
        if not all(_paths.genuine(*o) for o in own):
            return None
        _MOD['chain'] = base + tuple(o for o, _q, _m in own)
        _MOD['mods'] = (_tran_radau, _tran_radau_c, NumpyLU, ComplexKLUSolver)
    return _MOD['chain']


def _now(T, cir):
    """What the Python loop calls now, in `_chain`'s order."""
    tr_, rc, LU, ZS = _MOD['mods']
    return (rc._now(T, rc._MOD['mods'], cir)
            + (T._transform_loop, tr_._FrozenTransform.solve, tr_._transform_rhs,
               tr_._transform_back, LU.solve, ZS.solve_prepared))


class _Ctx:
    """One transient's struct and buffers for one core, size and reference
    row (the walk's fields named as `_tran_newton_c._Ctx`'s, whose
    `_walk_ready` sets them)."""

    __slots__ = (
        'A_bytes',
        'Ap_obj',
        'buf',
        'core',
        'ffi',
        'iref',
        'keep',
        'keep_ap',
        'klu_kk',
        'lu_addr',
        'n',
        'p_of',
        's',
        'walk',
        'walk_arms',
        'walk_dicts',
        'walk_els',
        'walk_hw_els',
        'walk_hwk',
        'walk_kern',
        'walk_ok',
        'walk_packs',
        'walk_stamp',
    )

    def __init__(self, core, n, iref, ffi, cfn_core, cffi_core, drv):
        m = n - 1
        self.core, self.n, self.iref, self.ffi = core, n, iref, ffi
        self.s = s = ffi.new('radau_tf_t *')
        self.buf = buf = {
            'Y': np.empty((3, n)), 'Yp': np.empty((3, n)), 'q': np.empty((3, n)),
            'iv': np.empty((3, n)), 'K': np.empty((3, n)), 'F': np.empty((3, n)),
            'R': np.empty(3 * m), 'b0': np.empty(m), 'x0': np.empty(m),
            'r1': np.empty(2 * m), 'x1': np.empty(2 * m), 'y1': np.empty(2 * m),
            'ab': np.empty(m), 'dY': np.empty(3 * m), 'stp': np.empty(n), 'tolv': np.empty(2),
            'qn': np.empty(n), 'u': np.empty((3, n)), 'seed': np.empty((3, n)),
            'dummy': np.zeros(max(n * n, 1)),
        }
        self.keep = []
        fb = ffi.from_buffer
        for k, arr in buf.items():
            h = fb('double *', arr)
            self.keep.append(h)
            setattr(s, k, h)
        _f, _c, dgesv, dgetrs = drv
        s.dgesv = ffi.cast('void *', dgesv)
        s.dgetrs = ffi.cast('void *', dgetrs)
        s.core = ffi.cast('void *', int(cffi_core.cast('uintptr_t', core.cs)))
        s.core_fn = ffi.cast('void *', int(cffi_core.cast('uintptr_t', cfn_core)))
        a = core.arrays['bini']
        self.keep.append(a)
        s.bini = ffi.cast('double *', a.ctypes.data)
        s.n, s.iref = n, iref
        s.walk_fn = ffi.NULL
        s.wn = 0
        ## (what does not change from step to step is set once: the product
        ## form, the passes' bits; the rest when it changes -- speed round 10,
        ## B3.5: a cast from `arr.ctypes.data` cost 13.5 k instructions, 21
        ## of `P`'s and `V`'s entries ~30 k a step)
        s.fused = _MOD['cmul']
        s.bits = _MOD['bits']
        self.A_bytes = self.p_of = self.Ap_obj = self.keep_ap = self.klu_kk = None
        self.lu_addr = None
        self.walk = self.walk_ok = None
        self.walk_els = self.walk_dicts = self.walk_packs = self.walk_kern = ()
        self.walk_hw_els = self.walk_hwk = ()
        self.walk_stamp, self.walk_arms = -1, 0


def solve(tr, ctx, fz, seed, src, provided_function, lims, nobypass, reltol, abstol, maxit):
    """The transform's Newton loop from the stages `seed` in C: the three
    converged stages as fresh full-width arrays, the source memo and the
    two factorisations' state as the Python loop leaves them -- or None:
    declined or handed back, nothing left behind that the Python loop
    would not make again, and the Python loop runs."""
    if not ENABLED:
        return _no('radau_tc:off')
    if fz is None or fz.lu is None or fz.prep is None:
        return _no('radau_tc:frozen')
    if provided_function is not None:
        return _no('radau_tc:pf')
    if lims:
        return _no('radau_tc:stateful')
    if not nobypass:
        return _no('radau_tc:bypass')
    drv = _driver or driver()
    if not drv:
        return _no('radau_tc:driver')
    ffi, cfn = drv[0], drv[1]
    chain = _MOD.get('chain') or _chain()
    if chain is None:
        return _no('radau_tc:patched')
    T = chain[0]
    if type(tr) is not T:
        return _no('radau_tc:class')
    cir = tr.cir
    if _now(T, cir) != chain[1:]:
        return _no('radau_tc:patched')
    NC = _MOD['NC']
    M = NC._MOD or NC._mods()
    if type(tr.toolkit) is not M['NumericToolkit']:
        return _no('radau_tc:toolkit')
    cd, td = cir.__dict__, tr.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd or 'limit' in cd:
        return _no('radau_tc:shadow_cir')
    if ('_memo_put' in td or '_source_at' in td or '_coupled_stage_system' in td
            or '_stages_converged' in td or '_stage_limiter_states' in td
            or '_limiters_at_rest' in td or '_newton_limiter' in td
            or '_transform_loop' in td):
        return _no('radau_tc:shadow_tr')
    if M['limiting'].CIRCUIT_LEVEL:
        return _no('radau_tc:stateful')
    memo = td.get('_u_memo')
    if memo is None or not M['companion'].U_MEMO:
        return _no('radau_tc:umemo')
    tc = M['core']
    if not (tc.CORE and tc.CORE_PASSES):
        return _no('radau_tc:core')
    core = tc.core_for(tr)
    if core is None:
        return _no('radau_tc:unservable')
    n, iref = core.n, ctx.iref
    m = n - 1
    if not (0 <= iref < n) or n < 2 or ctx.m != m or fz.m != m:
        return _no('radau_tc:iref')
    h = ctx.h
    if type(h) is not float and type(h) is not np.float64:
        return _no('radau_tc:h')
    if (type(reltol) is not float and type(reltol) is not int
            and type(reltol) is not np.float64) or type(reltol) is bool:
        return _no('radau_tc:tol')
    qn, Amat = ctx.qn, ctx.Amat
    if not (type(qn) is np.ndarray and qn.dtype == np.float64 and qn.shape == (n,)
            and type(Amat) is np.ndarray and Amat.dtype == np.float64 and Amat.shape == (3, 3)):
        return _no('radau_tc:ctx')
    if len(seed) != 3:
        return _no('radau_tc:seed')
    if not all(type(y) is np.ndarray and y.dtype == np.float64 for y in seed):
        seed = [np.array(y, dtype=float) for y in seed]
    if any(y.shape != (n,) for y in seed):
        return _no('radau_tc:seed')
    ## the two factorisations as `_FrozenTransform.solve` would use them
    lu, zs, prep = fz.lu, fz.zs, fz.prep
    _ZS, _LU = _MOD['mods'][3], _MOD['mods'][2]
    if (type(lu) is not _LU or type(zs) is not _ZS or lu.n != m or prep[1] != m
            or 'solve' in zs.__dict__ or 'solve_prepared' in zs.__dict__):
        return _no('radau_tc:frozen')
    if lu._info is not None and lu._info != 0:
        ## (singular: the Python raises numpy's error at its first solve)
        return _no('radau_tc:lu')
    if not (zs._symbolic and zs._numeric and prep[5] == zs._pattern):
        ## (the record's pattern not analysed and factored: the Python does it)
        return _no('radau_tc:klu')
    has_limiter = tr._newton_limiter() is not None
    if has_limiter:
        Tw = getattr(tr.epar, 'T', 300.0)
        if type(Tw) is not float and type(Tw) is not int and type(Tw) is not np.float64:
            return _no('radau_tc:T')
    rec = td.get('_radau_tc')
    if rec is None or rec.core is not core or rec.n != n or rec.iref != iref:
        _MOD['bits'] = tc._PASS_BITS['q'] | tc._PASS_BITS['i']
        cffi_core, cfn_core, _dg = tc.driver()
        rec = _Ctx(core, n, iref, ffi, cfn_core, cffi_core, drv)
        td['_radau_tc'] = rec
    s, buf = rec.s, rec.buf
    if has_limiter:
        if not NC._walk_ready(tr, rec):
            return _no('radau_tc:limit')
    elif s.wn:
        s.walk_fn = ffi.NULL
        s.wn = 0
        rec.walk = rec.walk_ok = None
    Tc = core.probe(tr.epar)
    if Tc.__class__ is str:
        return _no('radau_tc:ready')
    zk = _MOD.get('klu')
    if zk is None or zk[0] is not zs._lib:
        zk = _MOD['klu'] = (zs._lib, ctypes.cast(zs._lib.klu_z_refactor, ctypes.c_void_p).value,
                            ctypes.cast(zs._lib.klu_z_solve, ctypes.c_void_p).value)
    ## ---- from here the attempt leaves one trace: the step's source memo ----
    RC = _MOD['RC']
    analysis = tr.par.analysis
    keys = [(tt, analysis) for tt in ctx.tstage]
    try:
        had = [key in memo for key in keys]
    except TypeError:
        return _no('radau_tc:u')
    counts = {k: _PC.get(k, 0) for k in M['u_keys']}
    try:
        us = [src(tt) for tt in ctx.tstage]
    except Exception:                                          # noqa: BLE001
        ## (the call raises again in the Python loop, where it raised before)
        RC._undo_source(memo, keys, had, counts)
        return _no('radau_tc:u')
    if not all(type(u) is np.ndarray and u.dtype == np.float64 and u.shape == (n,)
               for u in us):
        RC._undo_source(memo, keys, had, counts)
        return _no('radau_tc:u')
    buf['u'][...] = us
    buf['seed'][...] = seed
    np.copyto(buf['qn'], qn)
    ## (the tableau once its values change; `P`'s and `V`'s entries once the
    ## step's cached tuple changes -- `_radau_frozen` keeps it per step size)
    ab = Amat.tobytes()
    if ab != rec.A_bytes:
        for k, v in enumerate(Amat.ravel().tolist()):
            s.A[k] = v
        rec.A_bytes = ab
    if fz.p is not rec.p_of:
        for k, z in enumerate(fz.p):
            s.p[2 * k] = float(z.real)
            s.p[2 * k + 1] = float(z.imag)
        for i, (a, b) in enumerate(fz.v):
            s.v[3 * i] = float(a)
            s.v[3 * i + 1] = float(b.real)
            s.v[3 * i + 2] = float(b.imag)
        rec.p_of = fz.p
    ## (the kept LU by its addresses -- a refilled one keeps them; the
    ## pattern's arrays by their identity -- a kept pattern keeps them)
    addr = lu.addr
    if addr != rec.lu_addr:
        s.lu_a = ffi.cast('double *', addr[0])
        s.lu_ipiv = ffi.cast('int64_t *', addr[1])
        rec.lu_addr = addr
    s.lu_info = -1 if lu._info is None else 0
    Ap, Ai = prep[2], prep[3]
    if Ap is not rec.Ap_obj or Ai is not rec.keep_ap[1]:
        rec.keep_ap = (Ap, Ai, ffi.from_buffer('int *', Ap), ffi.from_buffer('int *', Ai))
        s.Ap, s.Ai = rec.keep_ap[2], rec.keep_ap[3]
        rec.Ap_obj = Ap
    ax = ffi.from_buffer('double *', prep[4])
    s.Ax = ax
    kk = (zk[1], zk[2], zs._symbolic, zs._numeric, ctypes.addressof(zs._common))
    if kk != rec.klu_kk:
        s.klu_refactor = ffi.cast('void *', kk[0])
        s.klu_solve = ffi.cast('void *', kk[1])
        s.Symbolic = ffi.cast('void *', kk[2])
        s.Numeric = ffi.cast('void *', kk[3])
        s.Common = ffi.cast('void *', kk[4])
        rec.klu_kk = kk
    fresh = zs._fresh is prep
    s.klu_fresh = 1 if fresh else 0
    s.T = Tc
    s.h, s.reltol, s.abstol = float(h), float(reltol), float(abstol)
    s.rtol, s.maxiter = float(zs.REFACTOR_RESIDUAL_TOL), int(maxit)
    status = cfn(s)
    ## the factorisations' state, as the Python would have left it: the
    ## kept LU after its first `dgesv`, the numeric holding the record's
    ## refactor (a failed refactor is the Python's to make and count again)
    if s.lu_info >= 0:
        lu._info = int(s.lu_info)
    if s.refactored:
        zs.refactors += int(s.refactored)
        zs._fresh = prep
    if status != 1:
        RC._undo_source(memo, keys, had, counts)
        if status == 10:
            ## (the walk stopped: its tables are taken again next time)
            rec.walk_ok = None
        return _no('radau_tc:bail:' + BAIL.get(status, str(status)))
    it = int(s.iters)
    _PC['umemo:hit'] += 3 * (it - 1)
    _PC['radau_tc:served'] += 1
    _PC['radau_tc:iterations'] += it
    Y = buf['Y']
    return [Y[j].copy() for j in range(3)]
