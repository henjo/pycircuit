"""THE ADAPTIVE ERROR TEST IN C (speed round 8, stage 3; 2026-10-04).

The default transient (adaptive gear) judges every step attempt: the
charge LTE `Eg` of the active integrator (`compute_lte`: Gear-2's third
divided difference of the charges), mapped to solution units by `J^-1`
(the reference row cut and restored around `np.linalg.solve`), against the
tolerance `TRTOL (reltol ref + abstol)` with `ref` by `relref` (the
controller's running maximum, `sigglobal_reference`), normalised
(`normalised_error`) and its maximum taken -- ~370 k instructions an
attempt on a 20-MosLevel1 chain, some forty numpy calls around one LU of
82 k.  Here ONE C call: the third divided difference in numpy's
elementwise order, the reduced `J` column-major as numpy's solve copies it
and numpy's own `scipy_dgesv_64_`, the reference (the running maximum into
a fresh array), the tolerance and the normalised error's maximum.  Bit for
bit `StepController._charge_lte`, `_normalised` and `np.max` (pinned in
`tests/test_pinned_pairs.py`).

SERVED: `IntegralController` and `PIController` exactly (their laws stay
Python: `StepController._max_error` is the one call both make), Gear-2 and
the trapezoid on the charge form (`h_last2` known: from a run's third
step), the numeric toolkit, float64 vectors and `J` of the run's size.

EXACT OR NOT AT ALL.  Every input finite; the floating-point flags cleared
on entry and tested after each part whose numpy calls would check them --
an overflow, underflow, invalid operation or division by zero, a LAPACK
info, a non-finite solution, a tolerance that is not positive or a signed
zero where numpy's reductions may pair differently hands the attempt back
before anything is written, and the Python chain runs (numpy's warnings
and errors and the solve's fallback then happen as they did).  The
controller's running maximum is assigned only on success, and the flags
are left clear, as the chain's last numpy call leaves them.

DECLINES before anything is touched (`_paths.no`): the switch, no build,
another controller or integrator, a patched or shadowed piece of the chain
(the controller's methods, the module functions they call, the
integrator's `compute_lte`, the toolkit's solve), the g-form step and the
Euler order drop, a non-numeric toolkit, an input that is not a float64
array of the run's size, a step that is not a float.

`ENABLED` (env `PYCIRCUIT_LTE_C=0`) is read on every call; `STATUS` says
why the path is off where it is.  History: `doc/transient_history.md`,
`_tran_lte_c.py`.
"""
import math
import os

import numpy as np

from pycircuit.circuit import _cbuild, _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_LTE_C', '1') != '0'
STATUS = 'not loaded'

LTE_C = _cbuild.csrc('lte.c')
LTE_CDEF = """
typedef struct {
    long n, iref, nn, relref, have_run, abs_vec, has_mask;
    double h1, h2, h3, s12, s23, s123, coef, reltol, abstol, trtol;
    const double *xc, *xl, *q, *q1, *q2, *q3, *J, *run, *abst, *mask;
    double *A, *b, *lte, *ref, *out, *tq;
    int64_t *ipiv;
    void *dgesv;
    double err;
    long status;
} lte_t;
long hdl_fn(lte_t *s);
"""

#: the C's status codes past 1 (done): why it handed the attempt back
BAIL = {2: 'nonfinite', 3: 'flags', 4: 'dgesv', 5: 'lte', 6: 'tolerance'}

_driver = None
_G = {}
_F64 = np.dtype(np.float64)
_RELREF = {'pointlocal': 0, 'alllocal': 1}
_NPF = np.float64
#: numpy's pieces of the chain, as imported (before any caller could patch)
_NP = (np.linalg.solve, np.maximum, np.concatenate)


def driver():
    """`(ffi, cfn, dgesv address)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        import ctypes

        from pycircuit.circuit import _tran_newton_c
        try:
            npl = _tran_newton_c._openblas(np, 'libscipy_openblas64_*.so')
            if npl is None:
                STATUS = "off (numpy's own OpenBLAS not found)"
                return None
            dgesv = ctypes.cast(npl.scipy_dgesv_64_, ctypes.c_void_p).value
        except (OSError, AttributeError) as e:
            STATUS = f'off ({e})'
            return None
        from pycircuit.circuit import _hdl_cbackend as cb
        try:
            ffi, cfn, _key, _cold, _secs = cb.load_kernel(LTE_C, LTE_CDEF)
        except (cb.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, dgesv)
        STATUS = 'c'
    return _driver or None


def _globals():
    """The chain's pieces as defined, read once a call finds every one its
    module's own (one taken under a caller's patch would hold the patch:
    refused, and read again on a later call -- the call declines meanwhile);
    a call compares what it would call now against them.  None: not yet."""
    if not _G:
        from pycircuit.circuit import _numeric, analysis, integrator
        from pycircuit.circuit import stepcontroller as sc
        from pycircuit.circuit.toolkit import NumericToolkit
        SC, IG, AN, NU = (sc.StepController, 'pycircuit.circuit.integrator',
                          'pycircuit.circuit.analysis', 'pycircuit.circuit._numeric')
        SCM = 'pycircuit.circuit.stepcontroller'
        g2, tr = integrator.Gear2Integrator, integrator.TrapezoidalIntegrator
        own = ((SC._charge_lte, 'StepController._charge_lte'),
               (SC._normalised, 'StepController._normalised'),
               (SC._reference, 'StepController._reference'),
               (SC.tolerance, 'StepController.tolerance'),
               (sc.normalised_error, 'normalised_error'),
               (sc.sigglobal_reference, 'sigglobal_reference'),
               (sc.IntegralController, 'IntegralController'),
               (sc.PIController, 'PIController'),
               (integrator.third_divided_difference, 'third_divided_difference',
                'pycircuit.circuit._lte_kernels'),
               (g2, 'Gear2Integrator', IG), (g2.compute_lte, 'Gear2Integrator.compute_lte', IG),
               (tr, 'TrapezoidalIntegrator', IG),
               (tr.compute_lte, 'TrapezoidalIntegrator.compute_lte', IG),
               (analysis.remove_row_col, 'remove_row_col', AN),
               (analysis._reduce_ndarray, '_reduce_ndarray', AN),
               (_numeric.linearsolver, 'linearsolver', NU), (_numeric.array, 'array', NU),
               (NumericToolkit, 'NumericToolkit', 'pycircuit.circuit.toolkit'))
        if not all(_paths.genuine(o[0], o[1], o[2] if len(o) > 2 else SCM) for o in own):
            return None
        _G['mods'] = (sc, integrator, analysis)
        _G['controllers'] = (sc.IntegralController, sc.PIController)
        _G['tk'] = NumericToolkit
        _G['chain'] = (SC._charge_lte, SC._normalised, SC._reference, SC.tolerance,
                       sc.normalised_error, sc.sigglobal_reference,
                       integrator.third_divided_difference, analysis.remove_row_col,
                       analysis._reduce_ndarray, _NP[0], _numeric.linearsolver,
                       _NP[1], _NP[2], _numeric.array)
        _G['formula'] = {g2: (0, g2.compute_lte), tr: (1, tr.compute_lte)}
    return _G


class _Ctx:
    """One controller's struct and buffers for one size and reference row:
    every pointer the C holds is a buffer of this context's, set once."""

    __slots__ = ('buf', 'ipiv', 'iref', 'keep', 'n', 's')

    def __init__(self, ffi, dgesv, n, iref):
        m = n - 1
        self.n, self.iref = n, iref
        self.s = s = ffi.new('lte_t *')
        self.buf = buf = {
            'xc': np.empty(n), 'xl': np.empty(n), 'q': np.empty(n), 'ql': np.empty((3, n)),
            'J': np.empty((n, n)), 'run': np.empty(n),
            'abst': np.empty(n), 'A': np.empty(m * m), 'b': np.empty(m), 'lte': np.empty(n),
            'ref': np.empty(n), 'out': np.empty(n), 'tq': np.empty(n),
            'mask': np.empty(n),
        }
        self.ipiv = np.empty(m, dtype=np.int64)
        fb = ffi.from_buffer
        self.keep = []
        for k, arr in buf.items():
            if k == 'ql':
                ## (the three past charges, most recent first: one block)
                for row, name in enumerate(('q1', 'q2', 'q3')):
                    self.keep.append(fb('double *', arr[row]))
                    setattr(s, name, self.keep[-1])
                continue
            h = fb('double *', arr)
            self.keep.append(h)
            setattr(s, k, h)
        self.keep.append(fb('int64_t *', self.ipiv))
        s.ipiv = self.keep[-1]
        s.dgesv = ffi.cast('void *', dgesv)
        s.n, s.iref = n, iref


def _vec(a, n):
    return type(a) is np.ndarray and a.dtype is _F64 and a.shape == (n,)


def max_error(ctrl, s):
    """`StepController._max_error` in C: `(err, p)`, or None -- declined or
    handed back, nothing written, the Python chain runs."""
    if not ENABLED:
        return _no('lte_c:off')
    drv = _driver or driver()
    if not drv:
        return _no('lte_c:driver')
    G = _G or _globals()
    if G is None:
        return _no('lte_c:patched')
    cls = type(ctrl)
    if cls is not G['controllers'][0] and cls is not G['controllers'][1]:
        return _no('lte_c:controller')
    ai = s.active_integrator
    form = G['formula'].get(type(ai))
    if form is None:
        return _no('lte_c:integrator')
    tk = s.toolkit
    if type(tk) is not G['tk']:
        return _no('lte_c:toolkit')
    sc, ig, an = G['mods']
    cd = ctrl.__dict__
    if ('_charge_lte' in cd or '_normalised' in cd or '_reference' in cd or 'tolerance' in cd
            or 'compute_lte' in ai.__dict__ or type(ai).compute_lte is not form[1]
            or (cls._charge_lte, cls._normalised, cls._reference, cls.tolerance,
                sc.normalised_error, sc.sigglobal_reference, ig.third_divided_difference,
                an.remove_row_col, an._reduce_ndarray, np.linalg.solve, tk.linearsolver,
                tk.maximum, tk.concatenate, tk.array) != G['chain']):
        return _no('lte_c:patched')
    h, hl, hl2 = s.h_curr, s.h_last, s.h_last2
    if s.no_history or hl2 is None:
        return _no('lte_c:formula')
    th, tl, tl2 = type(h), type(hl), type(hl2)
    if not ((th is float or th is _NPF) and (tl is float or tl is _NPF)
            and (tl2 is float or tl2 is _NPF)):
        return _no('lte_c:h')
    ## (the integrator's scalars -- its coefficient, the divided
    ## difference's step sums -- in Python floats: the IEEE operations
    ## numpy's scalars make, without their warnings; one that overflows
    ## takes the Python chain, which warns as before)
    h, hl, hl2 = float(h), float(hl), float(hl2)
    try:
        coef = -h * (h + hl) if form[0] == 0 else -(h ** 2)
    except ArithmeticError:
        return _no('lte_c:h')
    s12 = h + hl
    s23 = hl + hl2
    s123 = s12 + hl2
    if not (math.isfinite(coef) and math.isfinite(s12) and math.isfinite(s23)
            and math.isfinite(s123)):
        return _no('lte_c:h')
    J, iref = s.J, s.irefnode
    if not (type(J) is np.ndarray and J.dtype is _F64 and J.ndim == 2):
        return _no('lte_c:x')
    n = J.shape[0]
    ql = s.q_last_hist
    xc, xl, q = s.x_curr, s.x_last, s.q_curr
    if not (n >= 2 and J.shape[1] == n and type(iref) is int and 0 <= iref < n
            and _vec(xc, n) and _vec(xl, n) and _vec(q, n)):
        return _no('lte_c:x')
    ## (the charge ring: a (k, n) block, most recent first, or a sequence)
    if type(ql) is np.ndarray:
        if not (ql.dtype is _F64 and ql.ndim == 2 and ql.shape[0] >= 3 and ql.shape[1] == n):
            return _no('lte_c:x')
        ql3 = ql[:3]
    elif type(ql) in (list, tuple) and len(ql) >= 3 and all(_vec(v, n) for v in ql[:3]):
        ql3 = ql[:3]
    else:
        return _no('lte_c:x')
    rr = ctrl.relref
    if type(rr) is not str:
        return _no('lte_c:relref')
    mode = _RELREF.get(rr, 2)
    run = getattr(ctrl, '_ref_running', None) if mode else None
    if run is not None and not _vec(run, n):
        return _no('lte_c:x')
    ab, reltol, trtol, nn = s.abstol, s.reltol, s.TRTOL, s.n_nodes
    abs_vec = type(ab) is np.ndarray
    if abs_vec:
        if not _vec(ab, n):
            return _no('lte_c:x')
    elif type(ab) not in (float, int, np.float64):
        return _no('lte_c:x')
    if type(reltol) not in (float, int) or type(trtol) not in (float, int):
        return _no('lte_c:x')
    if nn is not None and type(nn) is not int:
        return _no('lte_c:x')
    ctx = cd.get('_lte_c')
    if ctx is None or ctx.n != n or ctx.iref != iref:
        ffi, _cfn, dgesv = drv
        ctx = cd['_lte_c'] = _Ctx(ffi, dgesv, n, iref)
    buf, st = ctx.buf, ctx.s
    np.copyto(buf['xc'], xc)
    np.copyto(buf['xl'], xl)
    np.copyto(buf['q'], q)
    np.copyto(buf['ql'], ql3)
    np.copyto(buf['J'], J)
    if run is not None:
        np.copyto(buf['run'], run)
    if abs_vec:
        np.copyto(buf['abst'], ab)
        st.abstol = 0.0
    else:
        st.abstol = float(ab)
    st.abs_vec = 1 if abs_vec else 0
    mask = s.lte_mask
    if mask is None:
        st.has_mask = 0
    else:
        if not (type(mask) is np.ndarray and mask.shape == (n,)):
            return _no('lte_c:x')
        np.copyto(buf['mask'], mask)
        st.has_mask = 1
    st.have_run = 0 if run is None else 1
    st.relref = mode
    st.nn = -1 if nn is None else nn
    st.h1, st.h2, st.h3, st.coef = h, hl, hl2, coef
    st.s12, st.s23, st.s123 = s12, s23, s123
    st.reltol, st.trtol = float(reltol), float(trtol)
    status = drv[1](st)
    if status != 1:
        return _no('lte_c:bail:' + BAIL.get(status, str(status)))
    _PC['lte_c:served'] += 1
    if mode:
        ctrl._ref_running = buf['out'].copy()
    return float(st.err), 3.0
