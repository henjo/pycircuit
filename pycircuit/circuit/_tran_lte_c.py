"""THE ADAPTIVE ERROR TEST IN C: `StepController._max_error` --
`_charge_lte`, `_normalised` and their maximum, the one call the integral
and PI controllers make on every attempt -- as one call of `csrc/lte.c`.

The Python judges an attempt by the charge LTE `Eg` of the active
integrator (`compute_lte`: Gear-2's or the trapezoid's third divided
difference of the charges, Euler's first difference at an order drop),
mapped to solution units by `J^-1` (the reference row cut and restored
around `np.linalg.solve`), against the tolerance
`TRTOL (reltol ref + abstol)` with `ref` by `relref` (the controller's
running maximum, `sigglobal_reference`), normalised (`normalised_error`)
and its maximum taken.  Here ONE C call: the divided difference in numpy's
elementwise order, the reduced `J` column-major as numpy's solve copies it
and numpy's own `scipy_dgesv_64_` (or the run's `StructuralKLU`), the
reference (the running maximum into a fresh array), the tolerance and the
normalised error's maximum.  Bit for bit `StepController._charge_lte`,
`_normalised` and `np.max` (pinned in `tests/test_pinned_pairs.py`).

OWNS, on the step controller: `_lte_c` (a `_Ctx`: the struct, its buffers,
the KLU binding, the `served` mark -- the watch counter at its last serve,
None after any other outcome; the C loop enters only on a mark equal to
the counter now).  A serve writes the controller's `_ref_running` (any
`relref` but pointlocal) and adds the refactor `StructuralKLU.solve`
counts to the run's solver (`refactors`).  The C loop binds a `_Ctx` of
its own (`_tran_loop_c._bind`).  `_paths` prefix `lte_c:`:
`lte_c:served`, `lte_c:sparse_served`, `lte_c:euler`, the declines below
and the hand-backs (`lte_c:bail:<why>`); the C loop re-creates the first
three for the tests it makes.  Switch: `ENABLED` (env `PYCIRCUIT_LTE_C=0`),
read on every call and by the C loop on every adaptive entry.

SERVED: `IntegralController` and `PIController` exactly (their laws stay
Python: `StepController._max_error` is the one call both make), Gear-2 and
the trapezoid on the charge form (`h_last2` known: from a run's third
step), Euler's step at an order drop, the numeric toolkit, float64 vectors
and `J` of the run's size, a dense solve or the run's `StructuralKLU` once
its first factor is made.

EXACT OR NOT AT ALL.  Every input finite; the floating-point flags cleared
on entry and tested after each part whose numpy calls would check them --
an overflow, underflow, invalid operation or division by zero, a LAPACK
info, a non-finite solution, a tolerance that is not positive or a signed
zero where numpy's reductions may pair differently hands the attempt back
before anything is written, and the Python chain runs (numpy's warnings
and errors and the solve's fallback then happen as they did).  The
controller's running maximum is assigned only on success, and the flags
are left clear, as the chain's last numpy call leaves them.

DECLINES before anything is touched (`_paths.no`): the switch, no build (or
numpy's OpenBLAS not found), another controller or integrator, a patched or
shadowed piece of the chain (the controller's methods, the module functions
they call, the integrator's `compute_lte`, the toolkit's solve), a sparse
solver that is not the run's `AutoSolver` over `StructuralKLU` or has made
no factor yet, the g-form step (`h_last2` not yet known) and a step with no
history, a non-numeric toolkit, an input that is not a float64 array of the
run's size, a step that is not a float, a `relref` that is not a str.

`STATUS` says why the path is off where it is.  History:
`doc/transient_history.md`, `_tran_lte_c.py`.
"""
import math
import os

import numpy as np

from pycircuit.circuit import _cbuild, _paths, _watch

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_LTE_C', '1') != '0'
STATUS = 'not loaded'

LTE_C = _cbuild.csrc('lte.c')
LTE_CDEF = _cbuild.csrc('lte.h')

#: the codes lte.h names: its results, its formulas, its `relref` modes
_CODES = _cbuild.constants(LTE_CDEF)
#: `{code: why}` of the results past 1 (done): why the C handed the attempt
#: back
BAIL = _cbuild.bails(_CODES, 'LTE_BAIL_')

_driver = None
_G = {}
_F64 = np.dtype(np.float64)
_RELREF = {'pointlocal': _CODES['LTE_RELREF_POINTLOCAL'], 'alllocal': _CODES['LTE_RELREF_ALLLOCAL']}
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
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(LTE_C, LTE_CDEF, 'pcx_lte')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, dgesv)
        STATUS = 'c'
    return _driver or None


#: the companions by their pcx_abi.h codes: the estimate's coefficient is
#: Gear-2's, the trapezoid's or Euler's (the step after a breakpoint)
_ABI = _cbuild.constants(_cbuild.csrc('pcx_abi.h'))
_GEAR2, _TRAP, _EULER = (_ABI['PCX_FORMULA_GEAR2'], _ABI['PCX_FORMULA_TRAPEZOIDAL'],
                         _ABI['PCX_FORMULA_EULER'])


def _globals():
    """The chain's pieces as defined, read once a call finds every one its
    module's own (one taken under a caller's patch would hold the patch:
    refused, and read again on a later call -- the call declines meanwhile);
    a call compares what it would call now against them.  None: not yet."""
    if not _G:
        from pycircuit.circuit import _numeric, analysis, integrator
        from pycircuit.circuit import stepcontroller as sc
        from pycircuit.circuit.linearsolver import AutoSolver, StructuralKLU
        from pycircuit.circuit.toolkit import NumericToolkit
        SC, IG, AN, NU = (sc.StepController, 'pycircuit.circuit.integrator',
                          'pycircuit.circuit.analysis', 'pycircuit.circuit._numeric')
        SCM = 'pycircuit.circuit.stepcontroller'
        LS = 'pycircuit.circuit.linearsolver'
        g2, tr = integrator.Gear2Integrator, integrator.TrapezoidalIntegrator
        eu = integrator.EulerIntegrator
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
               (eu, 'EulerIntegrator', IG), (eu.compute_lte, 'EulerIntegrator.compute_lte', IG),
               (analysis.remove_row_col, 'remove_row_col', AN),
               (analysis._reduce_ndarray, '_reduce_ndarray', AN),
               (_numeric.linearsolver, 'linearsolver', NU), (_numeric.array, 'array', NU),
               (NumericToolkit, 'NumericToolkit', 'pycircuit.circuit.toolkit'),
               (AutoSolver.solve, 'AutoSolver.solve', LS),
               (StructuralKLU.solve, 'StructuralKLU.solve', LS),
               (StructuralKLU._ours, 'StructuralKLU._ours', LS))
        if not all(_paths.genuine(o[0], o[1], o[2] if len(o) > 2 else SCM) for o in own):
            return None
        _G['mods'] = (sc, integrator, analysis)
        _G['controllers'] = (sc.IntegralController, sc.PIController)
        _G['tk'] = NumericToolkit
        _G['sparse'] = (AutoSolver, StructuralKLU, AutoSolver.solve, StructuralKLU.solve,
                        StructuralKLU._ours)
        _G['chain'] = (SC._charge_lte, SC._normalised, SC._reference, SC.tolerance,
                       sc.normalised_error, sc.sigglobal_reference,
                       integrator.third_divided_difference, analysis.remove_row_col,
                       analysis._reduce_ndarray, _NP[0], _numeric.linearsolver,
                       _NP[1], _NP[2], _numeric.array)
        ## (the estimate's coefficient: gear's, the trapezoid's, Euler's --
        ## the step after a breakpoint, where the order drops)
        _G['formula'] = {g2: (_GEAR2, g2.compute_lte), tr: (_TRAP, tr.compute_lte),
                         eu: (_EULER, eu.compute_lte)}
    return _G


class _Ctx:
    """One controller's struct and buffers for one size and reference row:
    every pointer the C holds is a buffer of this context's, set once."""

    __slots__ = ('buf', 'ffi', 'ipiv', 'iref', 'j_inplace', 'jbuf', 'keep', 'klu', 'klu_keep',
                 'klu_sym', 'n', 's', 'served')

    def __init__(self, ffi, dgesv, n, iref):
        m = n - 1
        self.n, self.iref, self.ffi = n, iref, ffi
        self.s = s = ffi.new('lte_t *')
        self.buf = buf = {
            'xc': np.empty(n), 'xl': np.empty(n), 'q': np.empty(n), 'ql': np.empty((3, n)),
            'J': np.empty((n, n)), 'run': np.empty(n),
            'abst': np.empty(n), 'A': np.empty(m * m), 'b': np.empty(m), 'lte': np.empty(n),
            'ref': np.empty(n), 'out': np.empty(n), 'tq': np.empty(n),
            'mask': np.empty(n), 'iq1': np.empty(n),
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
        ## (J's own buffer, where it is copied: below
        ## `_cbuild.INPLACE_MIN_ENTRIES`, or not C-contiguous)
        self.jbuf = s.J
        self.j_inplace = n * n >= _cbuild.INPLACE_MIN_ENTRIES
        ## (the watch counter at this context's last serve, None after any
        ## other outcome: the adaptive loop in C enters only on a serve
        ## nothing watched has moved since -- `_tran_loop_c`)
        self.served = None
        s.n, s.iref = n, iref
        s.klu_sym = ffi.NULL
        self.klu = self.klu_sym = None
        self.klu_keep = ()


def _klu_bind(ctx, sk):
    """Point the C's sparse solve at `sk`, the run's `StructuralKLU`: its
    pattern in CSC, its places in the reduced row-major J (`take`), its
    analysis and its common, the library's entries; the reduced J and the
    residual's buffers.  False -- declined -- where its size is not this
    context's."""
    from pycircuit.circuit import _tran_newton_c
    m = ctx.n - 1
    nd = _tran_newton_c.driver()
    if sk.m != m or nd is None:
        return False
    ffi, s, fb = ctx.ffi, ctx.s, ctx.ffi.from_buffer
    ## (numpy's `cblas_dgemv`, as the C Newton's: `arr @ x` is that call)
    s.dgemv = ffi.cast('void *', nd[2]['dgemv'])
    ct, lib = sk._ct, sk._lib
    ktake = np.ascontiguousarray(sk.take, dtype=np.int64)
    kAx, R, x, y = np.empty(max(sk.nnz, 1)), np.empty(m * m), np.empty(m), np.empty(m)
    keep = (sk, ktake, kAx, R, x, y, fb('int *', sk.Ap), fb('int *', sk.Ai),
            fb('int64_t *', ktake), fb('double *', kAx), fb('double *', R), fb('double *', x),
            fb('double *', y))
    s.knnz = sk.nnz
    s.kAp, s.kAi, s.ktake, s.kAx, s.R, s.x, s.y = keep[6:]

    def addr(f):
        return ffi.cast('void *', ct.cast(f, ct.c_void_p).value)
    s.klu_refactor, s.klu_solve = addr(lib.klu_refactor), addr(lib.klu_solve)
    s.klu_common = ffi.cast('void *', ct.addressof(sk._common))
    s.klu_sym = ffi.cast('void *', sk._symbolic)
    ctx.klu, ctx.klu_sym, ctx.klu_keep = sk, sk._symbolic, keep
    return True


def _vec(a, n):
    return type(a) is np.ndarray and a.dtype is _F64 and a.shape == (n,)


def max_error(ctrl, s):
    """`StepController._max_error` in C for `ctrl` on the attempt `s` (a
    `StepLTEInputs`): `(err, p)`, or None -- declined or handed back,
    nothing of the run's written, the Python chain runs.

    Every call first clears the context's `served` mark.  A serve sets it
    to the watch counter, writes `ctrl._ref_running` (any `relref` but
    pointlocal) and, sparse, adds the refactor to the solver's `refactors`;
    the context (`ctrl._lte_c`) is made, or remade for another size or
    reference row, before the call.

    In order, its sections: 1. the checks that decline (module docstring:
    DECLINES); 2. the integrator's scalars, in Python floats; 3. the
    context and the attempt's arrays copied in (or read where they stand);
    4. the call and the serve."""
    last = ctrl.__dict__.get('_lte_c')
    if last is not None:
        last.served = None
    if not ENABLED:
        return _no('lte_c:off')
    drv = _driver or driver()
    if not drv:
        return _no('lte_c:driver')
    ## ---- 1. the checks that decline: the build, the pieces, the solver,
    ## the controller, the integrator, the steps, the arrays' types ----
    G = _G or _globals()
    if G is None:
        return _no('lte_c:patched')
    ## (a sparse run solves `J^-1 Eg` on its structural factorisation:
    ## `StructuralKLU.solve`, twinned -- its refactor needs the run's numeric,
    ## which only the Python's first factor makes)
    ls = s.linsolver
    sk = None
    if ls is not None:
        sk = ls._choice
        AS, SK, as_solve, sk_solve, sk_ours = G['sparse']
        if not (type(ls) is AS and type(sk) is SK and AS.solve is as_solve
                and SK.solve is sk_solve and SK._ours is sk_ours and 'solve' not in ls.__dict__
                and 'solve' not in sk.__dict__ and '_ours' not in sk.__dict__):
            return _no('lte_c:sparse')
        if not (sk._numeric and sk._symbolic):
            return _no('lte_c:sparse_first')
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
    euler = form[0] == _EULER
    if s.no_history or (hl2 is None and not euler):
        return _no('lte_c:formula')
    th, tl, tl2 = type(h), type(hl), type(hl2)
    if not ((th is float or th is _NPF) and (tl is float or tl is _NPF or (euler and hl is None))
            and (tl2 is float or tl2 is _NPF or euler)):
        return _no('lte_c:h')
    ## ---- 2. the integrator's scalars ----
    ## (the integrator's scalars -- its coefficient, the divided
    ## difference's step sums -- in Python floats: the IEEE operations
    ## numpy's scalars make, without their warnings; one that overflows
    ## takes the Python chain, which warns as before)
    h = float(h)
    try:
        if euler:
            ## (Euler's: `-0.5 * scale`, the variable-step correction
            ## `2 h / (h + h_last)` where there is a last step; no step sums)
            hl = None if hl is None else float(hl)
            coef = -0.5 * (2.0 * h / (h + hl) if hl else 1.0)
            hl = hl2 = s12 = s23 = s123 = 0.0
        else:
            hl, hl2 = float(hl), float(hl2)
            coef = -h * (h + hl) if form[0] == _GEAR2 else -(h ** 2)
            s12 = h + hl
            s23 = hl + hl2
            s123 = s12 + hl2
    except ArithmeticError:
        return _no('lte_c:h')
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
    mode = _RELREF.get(rr, _CODES['LTE_RELREF_SIGGLOBAL'])
    run = (getattr(ctrl, '_ref_running', None)
           if mode != _CODES['LTE_RELREF_POINTLOCAL'] else None)
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
    ## ---- 3. the context and the attempt's arrays ----
    ctx = cd.get('_lte_c')
    if ctx is None or ctx.n != n or ctx.iref != iref:
        ffi, _cfn, dgesv = drv
        ctx = cd['_lte_c'] = _Ctx(ffi, dgesv, n, iref)
    buf, st = ctx.buf, ctx.s
    if sk is not None:
        if (ctx.klu is not sk or ctx.klu_sym != sk._symbolic) and not _klu_bind(ctx, sk):
            return _no('lte_c:sparse')
        ## (the run's numeric as it stands, and the validation's bound)
        st.klu_num = ctx.ffi.cast('void *', sk._numeric)
        st.resid_tol = float(sk.REFACTOR_RESIDUAL_TOL)
    elif ctx.klu is not None:
        st.klu_sym = ctx.ffi.NULL
        ctx.klu = ctx.klu_sym = None
        ctx.klu_keep = ()
    np.copyto(buf['xc'], xc)
    np.copyto(buf['xl'], xl)
    np.copyto(buf['q'], q)
    np.copyto(buf['ql'], ql3)
    if euler:
        ## (the last companion current, `iq_last[0]`)
        iql = s.iq_last_hist
        iq1 = iql[0] if type(iql) is np.ndarray and iql.ndim == 2 else (
            iql[0] if type(iql) in (list, tuple) and len(iql) else None)
        if not _vec(iq1, n):
            return _no('lte_c:x')
        np.copyto(buf['iq1'], iq1)
    ## (J read where it stands -- the C only reads it: an n x n copy a call
    ## saved -- else copied, `INPLACE_MIN_ENTRIES`)
    if ctx.j_inplace and J.flags.c_contiguous:
        st.J = ctx.ffi.from_buffer('double *', J)
    else:
        np.copyto(buf['J'], J)
        if ctx.j_inplace:
            st.J = ctx.jbuf
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
    st.formula = _CODES['LTE_FORMULA_EULER'] if euler else _CODES['LTE_FORMULA_DD3']
    st.h1, st.h2, st.h3, st.coef = h, hl, hl2, coef
    st.s12, st.s23, st.s123 = s12, s23, s123
    st.reltol, st.trtol = float(reltol), float(trtol)
    ## ---- 4. the call and the serve ----
    status = drv[1](st)
    if status != _CODES['LTE_DONE']:
        return _no('lte_c:bail:' + BAIL.get(status, str(status)))
    _PC['lte_c:served'] += 1
    ctx.served = _watch.now()
    if sk is not None:
        ## (the refactor `StructuralKLU.solve` counts)
        sk.refactors += 1
        _PC['lte_c:sparse_served'] += 1
    if euler:
        _PC['lte_c:euler'] += 1
    if mode:
        ctrl._ref_running = buf['out'].copy()
    return float(st.err), (2.0 if euler else 3.0)
