"""THE NEWTON SOLVE IN C (speed round 8, stage 4; 2026-10-04).

A multistep step's Newton spent ~750 k instructions an iteration on a
20-MosLevel1 chain in Python and numpy around ~250 k of C: the evaluate
core's wrapper (readiness over the elements, sixteen buffer handles, the
lookup), the reference row cut out of `J` and `F` and put back into `x`,
numpy's solve wrapper, the limiter's wrapper and the walk's per-element
checks, the convergence test's dozen small-array calls -- with about one
iteration a step.  Here ONE C call does the whole solve: per iteration the
evaluate core's own C (`_tran_core.CORE_C`, through its struct pointer), the
reduced `J` written column-major as numpy's solve copies it and `-F`, the
LU (numpy's own `scipy_dgesv_64_` for the plain Newton; SciPy's
`scipy_dgetrf_`/`scipy_dgetrs_` for the chord, which factors with SciPy's
OpenBLAS), `x + dx`, the limiter walk's own driver (`_hdl_climit._WALK_C`)
and `dx = x_next - x` re-taken where a limiter exists, the convergence test
in `nrsolver`'s order (`|J||x|` through numpy's own `scipy_cblas_dgemv64_`)
and the chord's shrink test; then the converged point's evaluation
(`jacobian_only`: G, C, q) in the same call.  Bit for bit `StandardNewton`,
`ChordNewton` and `jacobian_only` (the stage-0 prototype, pss_log
2026-10-04; pinned in `tests/test_pinned_pairs.py`).

TRANSACTIONAL.  Nothing is written on the transient or the devices before
the C has converged; anything unusual -- a LAPACK info, a non-finite value,
the walk stopping, a chord break, maxiter -- rolls back the attempt's one
trace (the step's source-memo entry and its counts), and `_newton` repeats
the same arithmetic from the same seed: the same answer, exception, warning
and statistics.

DECLINES, before anything is touched (`_paths.no`): the switch (and the
evaluate core's, `_tran_core.CORE`), no build or no OpenBLAS, a toolkit
that is not numeric, an instance shadow of the
circuit's passes or of the transient's machinery, a `provided_function`, no
step source memo, a stateful limiter, the continuation rescue, a caller's
solver or scaler, a linear solver that is not `AutoSolver`'s dense choice, a
circuit the core cannot serve, an integrator outside the core's enum, a
limiting element that is not its C kernel (structural reasons kept with the
walk: such a circuit declines at its first check), the core's readiness
(`_Core.probe`: counted once, by the path that follows).  Between two
iterations only the C runs, so readiness checked once a solve holds for
every iteration.

THE BRANCH CHECK follows as in `_newton`; its screen finds `C` at the
converged point in `_C_cache`.  Where it fires, the confirmation's
speculative solves overwrite the step's state, and `jacobian_only` runs
again after it, as before.

`ENABLED` (env `PYCIRCUIT_NEWTON_C=0`) is read on every call; `STATUS` says
why the path is off where it is.  History: `doc/transient_history.md`,
`_tran_newton_c.py`.
"""
import ctypes
import glob
import operator
import os
from itertools import repeat

import numpy as np

from pycircuit.circuit import _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_NEWTON_C', '1') != '0'
STATUS = 'not loaded'

NEWTON_C = r"""
#include <stdint.h>
#include <stddef.h>
#include <math.h>
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
typedef void (*dgetrf_t)(const int *, const int *, double *, const int *, int *, int *);
typedef void (*dgetrs_t)(const char *, const int *, const int *, const double *, const int *,
                         const int *, double *, const int *, int *, size_t);
typedef void (*dgemv_t)(int, int, int64_t, int64_t, double, const double *, int64_t,
                        const double *, int64_t, double, double *, int64_t);
typedef struct {
    const void *core; void *core_fn;
    double T; long formula; double a0, a1, a2, h, theta;
    const double *q1, *q2, *iq1, *u;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv, *dgetrf, *dgetrs, *dgemv;
    long n, iref, maxiter, chord, first_bits;
    double reltol; const double *abstol, *xtol;
    const double *Cin0; const double *xseed; const double *cx, *cC; long have_cache, hit;
    double *C0, *q0, *iq0, *Geq0, *F0, *J0;
    double *C, *q, *iq, *Geq, *F, *J;
    double *xf, *xnf, *x0f;
    double *A, *absJ;
    double *b, *Fr, *x, *xn, *xd, *ax, *Is;
    int64_t *ipiv64; int *ipiv32;
    long iters, last_eval, status;
    long fold_j; double *xc, *Cj, *qj, *iqj, *Geqj, *Jj;
} newton_t;

long hdl_fn(newton_t *s)
{
    long n = s->n, m = n - 1, r = s->iref, it, i, j;
    core_fn_t core = (core_fn_t) s->core_fn;
    walk_fn_t walk = (walk_fn_t) s->walk_fn;
    dgemv_t gemv = (dgemv_t) s->dgemv;
    double prev = 0.0;
    int have_prev = 0;
    int64_t m64 = m, one64 = 1, info64 = 0;
    int m32 = (int) m, one32 = 1, info32 = 0;
    char trans = 'N';
    s->iters = 0;
    s->status = 0;
    s->last_eval = -1;
    s->hit = 0;
    /* the seed, reduced (remove_row_col), and _C_lookup's equality test of
       the full seed it evaluates against the cached state */
    for (i = 0, j = 0; i < n; i++) { if (i == r) continue; s->x[j++] = s->xseed[i]; }
    s->Cin0 = 0;
    if (s->have_cache) {
        int eq = 1;
        for (i = 0; i < r; i++) if (!(s->cx[i] == s->x[i])) { eq = 0; break; }
        if (eq && !(s->cx[r] == 0.0)) eq = 0;
        if (eq) for (i = r; i < m; i++) if (!(s->cx[i + 1] == s->x[i])) { eq = 0; break; }
        if (eq) { s->Cin0 = s->cC; s->hit = 1; }
    }
    s->first_bits = s->Cin0 ? 5 : 7;
    for (it = 0; it < s->maxiter; it++) {
        const double *Cin;
        long bits;
        double *F, *J, *C, *q, *iq, *Geq;
        int full = (it == 0) || !s->chord;
        if (it == 0) { F = s->F0; J = s->J0; C = s->C0; q = s->q0; iq = s->iq0; Geq = s->Geq0; }
        else { F = s->F; J = s->J; C = s->C; q = s->q; iq = s->iq; Geq = s->Geq; }
        /* insert_row: the full x of this evaluation */
        for (i = 0; i < r; i++) s->xf[i] = s->x[i];
        s->xf[r] = 0.0;
        for (i = r; i < m; i++) s->xf[i + 1] = s->x[i];
        for (i = 0; i < n; i++) if (!isfinite(s->xf[i])) { s->status = 2; return 2; }
        if (full) {
            bits = (it == 0) ? s->first_bits : 7;
            Cin = (it == 0 && s->Cin0) ? s->Cin0 : C;
        } else {
            bits = 4;
            Cin = s->Cin0 ? s->Cin0 : s->C0;
        }
        core(s->core, s->xf, s->T, bits, s->formula, s->a0, s->a1, s->a2, s->h, s->theta,
             s->q1, s->q2, s->iq1, s->u, Cin, C, q, iq, Geq, F, J);
        s->last_eval = it;
        for (i = 0, j = 0; i < n; i++) {
            if (i == r) continue;
            s->Fr[j] = F[i];
            s->b[j] = -F[i];
            j++;
        }
        if (full) {
            long ri, ci, rr, cc;
            for (ri = 0, rr = 0; ri < n; ri++) {
                if (ri == r) continue;
                for (ci = 0, cc = 0; ci < n; ci++) {
                    double v;
                    if (ci == r) continue;
                    v = J[ri * n + ci];
                    s->A[cc * m + rr] = v;
                    s->absJ[rr * m + cc] = fabs(v);
                    cc++;
                }
                rr++;
            }
        }
        if (s->chord) {
            if (it == 0) {
                for (i = 0; i < m * m; i++) if (!isfinite(s->A[i])) { s->status = 3; return 3; }
                ((dgetrf_t) s->dgetrf)(&m32, &m32, s->A, &m32, s->ipiv32, &info32);
                if (info32 != 0) { s->status = 4; return 4; }
            }
            for (i = 0; i < m; i++) if (!isfinite(s->b[i])) { s->status = 5; return 5; }
            ((dgetrs_t) s->dgetrs)(&trans, &m32, &one32, s->A, &m32, s->ipiv32, s->b, &m32,
                                   &info32, 1);
            if (info32 != 0) { s->status = 6; return 6; }
        } else {
            ((dgesv_t) s->dgesv)(&m64, &one64, s->A, &m64, s->ipiv64, s->b, &m64, &info64);
            if (info64 != 0) { s->status = 7; return 7; }
        }
        for (i = 0; i < m; i++) s->xn[i] = s->x[i] + s->b[i];
        if (walk) {
            long e;
            for (i = 0; i < r; i++) { s->xnf[i] = s->xn[i]; s->x0f[i] = s->x[i]; }
            s->xnf[r] = 0.0;
            s->x0f[r] = 0.0;
            for (i = r; i < m; i++) { s->xnf[i + 1] = s->xn[i]; s->x0f[i + 1] = s->x[i]; }
            e = walk(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK, s->xnf, s->x0f, 0, s->wn);
            if (e < s->wn) { s->status = 8; return 8; }
            for (i = 0, j = 0; i < n; i++) { if (i == r) continue; s->xn[j++] = s->xnf[i]; }
            for (i = 0; i < m; i++) s->xd[i] = s->xn[i] - s->x[i];
        } else {
            for (i = 0; i < m; i++) s->xd[i] = s->b[i];
        }
        for (i = 0; i < m; i++) if (!isfinite(s->xn[i])) { s->status = 9; return 9; }
        for (i = 0; i < m; i++) s->ax[i] = fabs(s->xn[i]);
        gemv(101, 111, m, m, 1.0, s->absJ, m, s->ax, 1, 0.0, s->Is, 1);
        {
            int conv_x = 1, conv_f = 1;
            for (i = 0; i < m; i++) {
                double aF = fabs(s->Fr[i]);
                double Is = s->Is[i] + aF;
                double axi = s->ax[i], axo = fabs(s->x[i]);
                double mx = (axi >= axo) ? axi : axo;
                if (!(fabs(s->xd[i]) < s->reltol * mx + s->xtol[i])) conv_x = 0;
                if (!(aF < s->reltol * Is + s->abstol[i])) conv_f = 0;
            }
            s->iters = it + 1;
            if (conv_x && conv_f) {
                for (i = 0; i < m; i++) s->x[i] = s->xn[i];
                if (s->fold_j) {
                    /* jacobian_only at the converged point: G, C, q */
                    for (i = 0; i < r; i++) s->xc[i] = s->x[i];
                    s->xc[r] = 0.0;
                    for (i = r; i < m; i++) s->xc[i + 1] = s->x[i];
                    core(s->core, s->xc, s->T, 3, s->formula, s->a0, s->a1, s->a2, s->h,
                         s->theta, s->q1, s->q2, s->iq1, s->u, s->Cj, s->Cj, s->qj, s->iqj,
                         s->Geqj, s->F, s->Jj);
                }
                s->status = 1;
                return 1;
            }
        }
        if (s->chord) {
            double size = -INFINITY;
            for (i = 0; i < m; i++) {
                double v = fabs(s->xd[i]);
                if (isnan(v)) { size = v; break; }
                if (v > size) size = v;
            }
            if (!isfinite(size) || (have_prev && !(size < prev))) { s->status = 10; return 10; }
            prev = size;
            have_prev = 1;
        }
        for (i = 0; i < m; i++) s->x[i] = s->xn[i];
    }
    s->status = 11;
    return 11;
}
"""
NEWTON_CDEF = """
typedef struct {
    const void *core; void *core_fn;
    double T; long formula; double a0, a1, a2, h, theta;
    const double *q1, *q2, *iq1, *u;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv, *dgetrf, *dgetrs, *dgemv;
    long n, iref, maxiter, chord, first_bits;
    double reltol; const double *abstol, *xtol;
    const double *Cin0; const double *xseed; const double *cx, *cC; long have_cache, hit;
    double *C0, *q0, *iq0, *Geq0, *F0, *J0;
    double *C, *q, *iq, *Geq, *F, *J;
    double *xf, *xnf, *x0f;
    double *A, *absJ;
    double *b, *Fr, *x, *xn, *xd, *ax, *Is;
    int64_t *ipiv64; int *ipiv32;
    long iters, last_eval, status;
    long fold_j; double *xc, *Cj, *qj, *iqj, *Geqj, *Jj;
} newton_t;
long hdl_fn(newton_t *s);
"""

#: the C's status codes past 1 (converged): why it handed the solve back
BAIL = {2: 'x', 3: 'jnonfinite', 4: 'getrf', 5: 'bnonfinite', 6: 'getrs', 7: 'dgesv',
        8: 'walkstop', 9: 'xnext', 10: 'chordbreak', 11: 'maxiter'}
#: the counts the step's source pass makes (`_source_at`, `cir.u`): the
#: attempt's one trace, rolled back with the memo entry when it bails
_U_KEYS = ('umemo:hit', 'umemo:miss', 'umemo:unhashable', 'u:called')

_driver = None
_MOD = {}
_GETDICT = operator.attrgetter('__dict__')
_GETPACK = operator.methodcaller('get', '_hdl_cp')
#: the declines that are a property of the circuit: kept with its stamp
#: plan, so the next step declines at its first check
_KEPT = ('unservable', 'stateful', 'limit')


def _openblas(pkg, pattern):
    """The OpenBLAS `pkg` itself loaded (`<pkg>.libs/`): `CDLL` on its real
    path returns the instance already mapped."""
    d = os.path.join(os.path.dirname(pkg.__file__), '..', pkg.__name__ + '.libs')
    paths = glob.glob(os.path.join(d, pattern))
    return ctypes.CDLL(os.path.realpath(paths[0])) if paths else None


def driver():
    """`(ffi, cfn, LAPACK addresses)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            import scipy
            import scipy.linalg  # (its OpenBLAS mapped: the chord's)
            npl = _openblas(np, 'libscipy_openblas64_*.so')
            spl = _openblas(scipy, 'libscipy_openblas-*.so')
            if npl is None or spl is None:
                STATUS = "off (numpy's or SciPy's own OpenBLAS not found)"
                return None

            def addr(f):
                return ctypes.cast(f, ctypes.c_void_p).value
            lapack = {'dgesv': addr(npl.scipy_dgesv_64_),
                      'dgemv': addr(npl.scipy_cblas_dgemv64_),
                      'dgetrf': addr(spl.scipy_dgetrf_), 'dgetrs': addr(spl.scipy_dgetrs_)}
        except (ImportError, OSError, AttributeError) as e:
            STATUS = f'off ({e})'
            return None
        from pycircuit.circuit import _hdl_cbackend as cb
        try:
            ffi, cfn, _key, _cold, _secs = cb.load_kernel(NEWTON_C, NEWTON_CDEF)
        except (cb.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, lapack)
        STATUS = 'c'
    return _driver or None


def _chain():
    """The machinery `solve` stands in for, as defined -- the transient's
    class and its Newton methods, `nrsolver`'s two loops, the reference-row
    helpers `_newton` calls, the evaluate core's Python entry -- read once a
    call finds every piece its module's own (`_paths.genuine`); None: not
    yet.  A call then compares what the Python path would call now."""
    if 'chain' not in _MOD:
        from pycircuit.circuit import _tran_core, _tran_newton, nrsolver
        from pycircuit.circuit.transient import Transient as T
        TN, AN = 'pycircuit.circuit._tran_newton', 'pycircuit.circuit.analysis'
        own = ((T, 'Transient', 'pycircuit.circuit.transient'),
               (T._newton, '_StepNewton._newton', TN),
               (T._newton_limiter, '_StepNewton._newton_limiter', TN),
               (T._residual_and_jacobian, '_CompanionModel._residual_and_jacobian',
                'pycircuit.circuit._tran_companion'),
               (T._get_nrsolver, 'Analysis._get_nrsolver', AN),
               (T._get_scaler, 'Analysis._get_scaler', AN),
               (_tran_newton.refnode_removed, 'refnode_removed', 'pycircuit.circuit.dcanalysis'),
               (_tran_newton.remove_row_col, 'remove_row_col', AN),
               (nrsolver.StandardNewton.solve_system, 'StandardNewton.solve_system',
                'pycircuit.circuit.nrsolver'),
               (nrsolver.ChordNewton.solve_system, 'ChordNewton.solve_system',
                'pycircuit.circuit.nrsolver'),
               (_tran_core.evaluate, 'evaluate', 'pycircuit.circuit._tran_core'))
        if not all(_paths.genuine(*o) for o in own):
            return None
        _MOD['chain'] = tuple(o for o, _q, _m in own)
        _MOD['chain_mods'] = (_tran_newton, nrsolver.StandardNewton, nrsolver.ChordNewton,
                              _tran_core)
    return _MOD['chain']


def _mods():
    """The modules and names `solve` reads, imported once (a function-level
    `from pycircuit.circuit import ...` runs importlib's `_handle_fromlist`
    on every call)."""
    if not _MOD:
        from pycircuit.circuit import (
            _hdl_climit,
            _limiting,
            _stamp_plan,
            _tran_companion,
            _tran_core,
        )
        from pycircuit.circuit._lte_kernels import bdf2_alphas
        from pycircuit.circuit.analysis import insert_row
        from pycircuit.circuit.dcanalysis import refnode_removed
        from pycircuit.circuit.linearsolver import AutoSolver
        from pycircuit.circuit.toolkit import NumericToolkit
        _MOD.update(climit=_hdl_climit, limiting=_limiting, companion=_tran_companion,
                    plan_for=_stamp_plan._plan_for,
                    core=_tran_core, bdf2_alphas=bdf2_alphas, insert_row=insert_row,
                    refnode_removed=refnode_removed, AutoSolver=AutoSolver,
                    NumericToolkit=NumericToolkit)
    return _MOD


class _Ctx:
    """One transient's struct and buffers for one core, size and reference
    row: every pointer the C holds is a buffer of this context's, set once
    (a solve copies its inputs in and its state out)."""

    __slots__ = (
        'buf',
        'core',
        'ffi',
        'ipiv32',
        'ipiv64',
        'iref',
        'keep',
        'n',
        's',
        'tol',
        'walk',
        'walk_dicts',
        'walk_els',
        'walk_hw_els',
        'walk_kern',
        'walk_ok',
        'walk_packs',
    )

    def __init__(self, core, n, iref, ffi, cfn_core, cffi_core, lapack):
        m = n - 1
        nn = (n, n)
        self.core, self.n, self.iref, self.ffi = core, n, iref, ffi
        self.s = s = ffi.new('newton_t *')
        self.buf = buf = {
            'xf': np.empty(n), 'xnf': np.empty(n), 'x0f': np.empty(n),
            'A': np.empty(m * m), 'absJ': np.empty(m * m),
            'b': np.empty(m), 'Fr': np.empty(m), 'x': np.empty(m), 'xn': np.empty(m),
            'xd': np.empty(m), 'ax': np.empty(m), 'Is': np.empty(m),
            'F0': np.empty(n), 'F': np.empty(n), 'J0': np.empty(nn), 'J': np.empty(nn),
            'C0': np.empty(nn), 'q0': np.empty(n), 'iq0': np.empty(n), 'Geq0': np.empty(nn),
            'C': np.empty(nn), 'q': np.empty(n), 'iq': np.empty(n), 'Geq': np.empty(nn),
            'q1': np.empty(n), 'q2': np.empty(n), 'iq1': np.empty(n), 'u': np.empty(n),
            'abstol': np.empty(m), 'xtol': np.empty(m), 'xseed': np.empty(n),
            'cx': np.empty(n), 'cC': np.empty(nn),
            'xc': np.empty(n), 'Cj': np.empty(nn), 'qj': np.empty(n), 'iqj': np.empty(n),
            'Geqj': np.empty(nn), 'Jj': np.empty(nn),
        }
        self.ipiv64 = np.empty(max(m, 1), dtype=np.int64)
        self.ipiv32 = np.empty(max(m, 1), dtype=np.int32)
        fb = ffi.from_buffer
        self.keep = []
        for k, arr in buf.items():
            h = fb('double *', arr)
            self.keep.append(h)
            setattr(s, k, h)
        self.keep.append(fb('int64_t *', self.ipiv64))
        s.ipiv64 = self.keep[-1]
        self.keep.append(fb('int *', self.ipiv32))
        s.ipiv32 = self.keep[-1]
        for k in ('dgesv', 'dgetrf', 'dgetrs', 'dgemv'):
            setattr(s, k, ffi.cast('void *', lapack[k]))
        s.core = ffi.cast('void *', int(cffi_core.cast('uintptr_t', core.cs)))
        s.core_fn = ffi.cast('void *', int(cffi_core.cast('uintptr_t', cfn_core)))
        s.n, s.iref = n, iref
        s.walk_fn = ffi.NULL
        s.wn = 0
        self.tol = None
        self.walk = self.walk_ok = None
        self.walk_els = self.walk_dicts = self.walk_packs = self.walk_kern = ()
        self.walk_hw_els = ()


def _walk_fast(rec):
    """The walk's tables still stand as the last full setup left them: each
    class's kernel the one taken, every element's `__dict__` the one seen
    and without a `limit` shadow (nor a `vlimit` one on a hand-written
    limiter's element: PSP's), every pack the one mirrored -- tuple compares
    and mapped tests, at C speed."""
    for info, kern in rec.walk_kern:
        if info.get('_c_limit') is not kern:
            return False
    dicts = tuple(map(_GETDICT, rec.walk_els))
    if (dicts != rec.walk_dicts or any(map(operator.contains, dicts, repeat('limit')))
            or any(map(operator.contains, map(_GETDICT, rec.walk_hw_els), repeat('vlimit')))):
        return False
    return all(map(operator.is_, map(_GETPACK, dicts), rec.walk_packs))


def _walk_full(w, ffi):
    """`limit_walk`'s per-element setup without its call and its counts:
    True where every limiting element runs its C kernel."""
    entries, capable, kerns, addr, F, PR, mirror = (
        w.entries, w.capable, w.kerns, w.addr, w.F, w.PR, w.mirror)
    hc = _MOD['climit']
    if w.cF is None:
        w.prepare(ffi)
    for e in w.cap_idx:
        el = entries[e][1]
        ck = capable[e].get('_c_limit')
        if ck is not kerns[e] and not w.retake(e, ck, ffi):
            return False
        kern = kerns[e]
        if kern is None:
            return False
        d = el.__dict__
        if 'limit' in d or (type(capable[e]) is hc._Handwritten and 'vlimit' in d):
            return False
        cp = d.get('_hdl_cp')
        if cp is None:
            try:
                cp = kern.pack(el)
            except (TypeError, ValueError):
                cp = False
            d['_hdl_cp'] = cp
        if cp is False:
            return False
        F[e] = addr[e]
        if cp is not mirror[e]:
            mirror[e] = cp
            PR[e] = cp[0].ctypes.data
    return True


def _walk_ready(tr, rec):
    """The struct's walk for this solve: True, or False where a limiting
    element is not its C kernel.  A structural reason (an element whose
    class has no C limiter at all: PSP's) is kept with the walk object, so
    such a circuit declines at its first check; the others are re-checked."""
    hc = _MOD['climit']
    if not (hc.WALK and hc.ENABLED):
        return False
    w = hc._walk_for(tr.cir)
    if w is rec.walk:
        if rec.walk_ok is False:
            return False
        if rec.walk_ok and _walk_fast(rec):
            return True
    if not w.any_capable or len(w.cap_idx) != len(w.entries):
        rec.walk, rec.walk_ok = w, False
        return False
    drv = hc._walk_drv()
    if drv is None:
        return False
    wffi, wfn = drv
    ok = _walk_full(w, wffi)
    rec.walk, rec.walk_ok = w, (True if ok else None)
    if ok:
        els = tuple(w.entries[e][1] for e in w.cap_idx)
        rec.walk_els = els
        rec.walk_hw_els = tuple(w.entries[e][1] for e in w.cap_idx
                                if type(w.capable[e]) is hc._Handwritten)
        rec.walk_dicts = tuple(map(_GETDICT, els))
        rec.walk_packs = tuple(w.mirror[e] for e in w.cap_idx)
        rec.walk_kern = tuple({id(w.capable[e]): (w.capable[e], w.kerns[e])
                               for e in w.cap_idx}.values())
        s, ffi = rec.s, rec.ffi
        s.walk_fn = ffi.cast('void *', int(wffi.cast('uintptr_t', wfn)))
        s.wF = ffi.cast('void **', w.F.ctypes.data)
        s.wPR = ffi.cast('double **', w.PR.ctypes.data)
        s.wTI = ffi.cast('int64_t *', w.TI.ctypes.data)
        s.wNM = ffi.cast('int64_t *', w.NM.ctypes.data)
        s.wOFF = ffi.cast('int64_t *', w.OFF.ctypes.data)
        s.wK = ffi.cast('int64_t *', w.K.ctypes.data)
        s.wn = len(w.entries)
    return ok


def _undo_source(memo, key, had, counts):
    """The attempt's one trace, removed: the source memo's entry it made and
    the counts the pass added -- `_newton` makes them again."""
    if not had and memo is not None:
        memo.pop(key, None)
    for k, v in counts.items():
        _PC[k] = v


def _keep(tr, td, M, why):
    """A decline that is a property of the circuit, counted and kept with
    its stamp plan (`_KEPT`): `(plan, its _paths key)`."""
    key = 'newton_c:' + why
    td['_newton_c_no'] = (M['plan_for'](tr.cir), key)
    return _no(key)


def solve(tr, func, t, provided_function, seed, residual):
    """The multistep step's Newton solve, and its converged-point evaluation,
    in C: `(x, fj)` -- `x` the full-width solution as `_newton` returns it,
    `fj` `jacobian_only`'s `(None, J)` at it, or None where it must run again
    (a branch confirmation re-solved) -- or None: declined or handed back,
    nothing left behind, `_newton` runs."""
    if not ENABLED:
        return _no('newton_c:off')
    td, cd = tr.__dict__, tr.cir.__dict__
    kept = td.get('_newton_c_no')
    if kept is not None:
        ## (the circuit's own reason, under the plan the circuit's dict
        ## holds -- not `_plan_for`'s check, ~1 us inside a step, on every
        ## step of every such circuit: a plan gone stale is replaced at the
        ## circuit's next evaluation, which follows, and until then the
        ## decline is today's Newton, never another answer)
        if cd.get('_stamp_plan') is kept[0]:
            return _no(kept[1])
        del td['_newton_c_no']
    drv = _driver or driver()
    if not drv:
        return _no('newton_c:driver')
    ffi, cfn, lapack = drv
    M = _MOD or _mods()
    ## (the machinery the C stands in for, as defined: a subclass, or a
    ## piece patched on its class or module, takes the Python Newton)
    chain = M.get('chain') or _chain()
    if chain is None:
        return _no('newton_c:patched')
    T = chain[0]
    if type(tr) is not T:
        return _no('newton_c:class')
    tn, SN, CN, tcm = M['chain_mods']
    if (T._newton, T._newton_limiter, T._residual_and_jacobian, T._get_nrsolver,
            T._get_scaler, tn.refnode_removed, tn.remove_row_col, SN.solve_system,
            CN.solve_system, tcm.evaluate) != chain[1:]:
        return _no('newton_c:patched')
    if type(tr.toolkit) is not M['NumericToolkit']:
        return _no('newton_c:toolkit')
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd or 'limit' in cd:
        return _no('newton_c:shadow_cir')
    if ('get_diff' in td or '_companion_at' in td or '_C_at_state' in td or '_C_lookup' in td
            or '_source_at' in td or '_newton' in td or '_newton_limiter' in td
            or '_branch_after_solve' in td):
        return _no('newton_c:shadow_tr')
    if provided_function is not None:
        return _no('newton_c:pf')
    memo = td.get('_u_memo')
    if memo is None or not M['companion'].U_MEMO:
        return _no('newton_c:umemo')
    if M['limiting'].CIRCUIT_LEVEL:
        return _no('newton_c:stateful')
    if td.get('_stateful_lims'):
        return _keep(tr, td, M, 'stateful')
    par = tr.par
    if (getattr(tr, '_continuation_rescue', False) or par.nrsolver is not None
            or getattr(par, 'scaler', None) is not None):
        return _no('newton_c:solver')
    ls = tr._get_linearsolver()
    if not (type(ls) is M['AutoSolver'] and ls._choice is not None
            and ls._choice is ls._dense):
        return _no('newton_c:linsolver')
    tc = M['core']
    if not tc.CORE:
        ## (built on the evaluate core: its switch is this path's too)
        return _no('newton_c:core')
    core = tc.core_for(tr)
    if core is None:
        return _keep(tr, td, M, 'unservable')
    n, iref = core.n, tr.irefnode
    if not (0 <= iref < n) or n < 2:
        return _no('newton_c:iref')
    if not (type(seed) is np.ndarray and seed.dtype == np.float64 and seed.ndim == 1
            and seed.shape[0] == n and np.isfinite(seed).all()):
        return _no('newton_c:x')
    h = tr._dt
    h_last = tr._dt_last if tr._dt_last is not None else h
    active = tr.base_integrator.check_order_drop(h, h_last, tr._is_first_step)
    formula = tc.FORMULA.get(type(active).__name__)
    if formula is None:
        return _no('newton_c:formula')
    q1 = tc._row(tr._qlast[0], n)
    q2 = tc._row(tr._qlast[1], n) if formula == 0 else q1
    iq1 = tc._row(tr._iqlast[0], n) if formula in (2, 3) else q1
    if q1 is None or q2 is None or iq1 is None:
        return _no('newton_c:history')
    a0 = a1 = a2 = theta = 0.0
    if formula == 0:
        a0, a1, a2 = M['bdf2_alphas'](h, h_last)
    elif formula == 3:
        theta = active.theta_at(h)
    has_limiter = tr._newton_limiter() is not None
    if has_limiter:
        ## (a limiting element whose class has no C limiter at all -- PSP's:
        ## a property of the circuit, before any context is built for it)
        w = M['climit']._walk_for(tr.cir)
        if not w.any_capable or len(w.cap_idx) != len(w.entries):
            return _keep(tr, td, M, 'limit')
    rec = td.get('_newton_c')
    if rec is None or rec.core is not core or rec.n != n or rec.iref != iref:
        cffi_core, cfn_core, _dg = tc.driver()
        rec = _Ctx(core, n, iref, ffi, cfn_core, cffi_core, lapack)
        td['_newton_c'] = rec
    s, buf = rec.s, rec.buf
    if has_limiter:
        if not _walk_ready(tr, rec):
            if rec.walk_ok is False:
                return _keep(tr, td, M, 'limit')
            return _no('newton_c:limit')
    elif s.wn:
        s.walk_fn = ffi.NULL
        s.wn = 0
        rec.walk = rec.walk_ok = None
    T = core.probe(tr.epar)
    if T.__class__ is str:
        return _no('newton_c:ready')
    chord = bool(residual is not None and par.nrsolver is None
                 and tr._newton_option(par.chord_jacobian, 'chord_jacobian'))
    ## ---- from here the attempt leaves one trace: the step's source memo ----
    key = (t, par.analysis)
    try:
        had = key in memo
    except TypeError:
        return _no('newton_c:u')
    counts = {k: _PC.get(k, 0) for k in _U_KEYS}
    try:
        u = tr._source_at(t, None)
    except Exception:                                          # noqa: BLE001
        ## (the pass raises again in `_newton`, where it raised before)
        _undo_source(memo, key, had, counts)
        return _no('newton_c:u')
    if not (type(u) is np.ndarray and u.dtype == np.float64 and u.ndim == 1
            and u.shape[0] == n):
        _undo_source(memo, key, had, counts)
        return _no('newton_c:u')
    tol = tr._newton_tolerances()
    if tol is not rec.tol:
        rec.tol = tol
        buf['abstol'][...] = tol[2]
        buf['xtol'][...] = tol[3]
    np.copyto(buf['q1'], q1)
    np.copyto(buf['q2'], q2)
    np.copyto(buf['iq1'], iq1)
    np.copyto(buf['u'], u)
    np.copyto(buf['xseed'], seed)
    ## `_C_lookup`'s cached state at the seed (the multistep memo is empty;
    ## no stateful limiter -- declined -- and no bypass: its other routes)
    cached = td.get('_C_cache')
    s.have_cache = 0
    if (cached is not None and float(getattr(tr.epar, 'bypasstol', -1.0) or -1.0) < 0.0):
        cxa, cCa = cached
        if (type(cxa) is np.ndarray and cxa.shape == (n,) and type(cCa) is np.ndarray
                and cCa.dtype == np.float64 and cCa.shape == (n, n)):
            np.copyto(buf['cx'], cxa)
            np.copyto(buf['cC'], cCa)
            s.have_cache = 1
    s.T, s.formula = T, formula
    s.a0, s.a1, s.a2, s.h, s.theta = a0, a1, a2, float(h), float(theta)
    s.maxiter, s.chord = int(par.maxiter), (1 if chord else 0)
    s.reltol = float(par.reltol)
    s.fold_j = 1
    status = cfn(s)
    if status != 1:
        _undo_source(memo, key, had, counts)
        if status == 8:
            ## (the walk stopped: its tables are taken again next time)
            rec.walk_ok = None
        return _no('newton_c:bail:' + BAIL.get(status, str(status)))
    ## ---- converged: the converged point's state, the counts, the check ----
    iters = int(s.iters)
    if chord:
        _PC['chord:chosen'] += 1
    _PC['newton_c:served'] += 1
    _PC['newton_c:iterations'] += iters
    _PC['newton_c:jfold'] += 1
    stats = getattr(tr, 'statistics', None)
    if stats is not None:
        stats.newton_iterations += iters
    x_res = buf['x'].copy()
    insert_row = M['insert_row']
    x = insert_row(x_res, iref, tr.toolkit)
    Cj = buf['Cj'].copy()
    tr._q_cache = (x, buf['qj'].copy())
    tr._C_cache = (x, Cj)
    tr._iq, tr._Geq = buf['iqj'].copy(), buf['Geqj'].copy()
    tr._Cmat = Cj
    tr.active_integrator = active
    tr._companion_coeffs = active.companion_coefficients(h, h_last)
    tr._effective_method = type(active).__name__
    fj = (None, buf['Jj'].copy())
    if tr._branch_on():
        target = stats if stats is not None else tr
        before = getattr(target, 'branch_screens', 0)
        tr._branch_after_solve(M['refnode_removed'](func, iref, tr.toolkit), x_res)
        if getattr(target, 'branch_screens', 0) != before:
            ## the confirmation re-solved: `jacobian_only` again, after it
            fj = None
    return x, fj
