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
(`jacobian_only`: G, C, q -- C and q where nothing reads the step's
Jacobian, a fixed multistep run: speed round 9, C0) in the same call.  Bit for bit `StandardNewton`,
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
every iteration.  Between two SOLVES the walk's readiness, the core's and
the parameter reads (the solver options, the linear solver, the
tolerances) hold as stamped (`_watch`, speed round 9, stage 4): unchecked
while no dict they read has changed, what no dict watcher sees compared
on every call.

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

from pycircuit.circuit import _cbuild, _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_NEWTON_C', '1') != '0'
STATUS = 'not loaded'

NEWTON_C = _cbuild.csrc('newton.c')
NEWTON_CDEF = _cbuild.csrc('newton.h')

#: the C's status codes past 1 (converged): why it handed the solve back
#: the C's results, named in newton.h: {code: why} of its bails
_CODES = _cbuild.constants(NEWTON_CDEF)
BAIL = _cbuild.bails(_CODES, 'NEWTON_BAIL_')
#: the counts the step's source pass makes (`_source_at`, `cir.u`; the
#: source plan's own, `_stamp_plan._SK`, are added at `_mods`): the
#: attempt's one trace, rolled back with the memo entry when it bails
_U_KEYS = ('umemo:hit', 'umemo:miss', 'umemo:unhashable', 'u:called')

_driver = None
_MOD = {}
#: a parameter not stored (`_ParRec.bp_obj`)
_ABSENT = object()
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
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(NEWTON_C, NEWTON_CDEF)
        except (_cbuild.CompileError, OSError) as e:
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
            _watch,
        )
        from pycircuit.circuit._lte_kernels import bdf2_alphas
        from pycircuit.circuit.analysis import insert_row
        from pycircuit.circuit.dcanalysis import refnode_removed
        from pycircuit.circuit.linearsolver import AutoSolver
        from pycircuit.circuit.toolkit import NumericToolkit
        _MOD.update(climit=_hdl_climit, limiting=_limiting, companion=_tran_companion, watch=_watch,
                    plan_for=_stamp_plan._plan_for,
                    core=_tran_core, bdf2_alphas=bdf2_alphas, insert_row=insert_row,
                    refnode_removed=refnode_removed, AutoSolver=AutoSolver,
                    NumericToolkit=NumericToolkit,
                    u_keys=_U_KEYS + tuple(_stamp_plan._SK['u'].values())
                    + tuple(_stamp_plan._SD.values()))
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
        'par_arms',
        'par_held',
        'par_rec',
        's',
        'tol',
        'walk',
        'walk_arms',
        'walk_dicts',
        'walk_els',
        'walk_held',
        'walk_hw_els',
        'walk_hwk',
        'walk_kern',
        'walk_ok',
        'walk_packs',
        'walk_stamp',
        'walk_watch',
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
        self.walk_hw_els = self.walk_hwk = ()
        self.walk_stamp, self.walk_arms, self.walk_held = -1, 0, 0
        self.walk_watch = ()
        self.par_rec, self.par_arms, self.par_held = None, 0, 0


class _ParRec:
    """`solve`'s reads of the analysis parameters and of the transient's
    own caches, stamped (`_watch`): `tr.par`'s dicts watched, and of
    `tr.epar`'s its instance dict and its parameter table -- its values
    are not: the analyses write `analysis_kind` there around every solve's
    operating point (`analysis.analysis_kind`), which moved the counter of
    every stamp twice a solve (2026-10-06) -- and the one value read there,
    `bypasstol`, compared by its object on every call (`bp_obj`); the
    methods that read them and the caches they keep in the transient's dict
    compared on every call (`_par_stamped`)."""

    __slots__ = ('analysis', 'bp_obj', 'bypass', 'chord', 'epar', 'f_ls', 'f_tol', 'iref',
                 'ls', 'ls_default', 'maxiter', 'nb', 'nn', 'nrsolver', 'par', 'reltol',
                 'scaler', 'stamp', 'tol', 'tolc')


def _par_stamped(tr, td, par, rec):
    """The stamped record of `solve`'s parameter reads, where it stands:
    its stamp the counter, the same `par` and `epar`, the methods that read
    them (the linear solver's, the tolerances') the class's and not
    shadowed, the transient's caches of them the objects seen, and the
    circuit's size and reference row the same.  (What the reads feed is
    still computed every call: the chord option from the recorded value,
    the bypass test from the recorded tolerance -- their errors and side
    effects where they were.)"""
    P = rec.par_rec if rec is not None else None
    ep = _MOD['watch'].EPOCH
    epar = tr.epar
    if (P is None or ep is None or P.stamp != ep.value or P.par is not par
            or P.epar is not epar or epar._values.get('bypasstol', _ABSENT) is not P.bp_obj):
        return None
    T = type(tr)
    cir = tr.cir
    if (T._get_linearsolver is not P.f_ls or T._newton_tolerances is not P.f_tol
            or '_get_linearsolver' in td or '_newton_tolerances' in td
            or td.get('_default_linearsolver') is not P.ls_default
            or td.get('_newton_tol_cache') is not P.tolc
            or len(cir.nodes) != P.nn or len(cir.branches) != P.nb or tr.irefnode != P.iref):
        return None
    rec.par_held += 1
    return P


def _par_stamp(tr, td, par, rec, before, nrsolver, scaler, ls, tol):
    """Record and stamp `solve`'s parameter reads (after a solve that read
    them in full and passed every check): the values as read, nothing
    computed from them -- a read that raises leaves no record."""
    rec.par_rec = None
    P = _ParRec()
    T, cir = type(tr), tr.cir
    P.par, P.epar = par, tr.epar
    P.f_ls, P.f_tol = T._get_linearsolver, T._newton_tolerances
    P.ls_default = td.get('_default_linearsolver')
    P.tolc = td.get('_newton_tol_cache')
    P.nn, P.nb, P.iref = len(cir.nodes), len(cir.branches), tr.irefnode
    P.nrsolver, P.scaler, P.ls, P.tol = nrsolver, scaler, ls, tol
    try:
        P.chord, P.analysis = par.chord_jacobian, par.analysis
        P.maxiter, P.reltol = par.maxiter, par.reltol
        P.bypass = getattr(tr.epar, 'bypasstol', -1.0)
        ep_values, ep_params = tr.epar._values, tr.epar._parameters
    except Exception:                                          # noqa: BLE001
        return
    if type(ep_values) is not dict or type(ep_params) is not dict:
        return
    ## (`bypasstol` as `getattr` resolved it: `epar`'s instance dict and its
    ## parameter table watched, the value stored for it -- or its absence --
    ## compared on every call)
    P.bp_obj = ep_values.get('bypasstol', _ABSENT)
    dicts = [par.__dict__, tr.epar.__dict__, ep_params]
    v = par.__dict__.get('_values')
    if isinstance(v, dict):
        dicts.append(v)
    W = _MOD['watch']
    rec.par_arms, rec.par_held = W.counted(rec.par_arms, rec.par_held), 0
    P.stamp = W.arm(dicts, before) if rec.par_arms <= W.MAX_ARMS else -1
    rec.par_rec = P


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
    such a circuit declines at its first check; the others are re-checked.
    As stamped (`_watch`): the walk's elements' dicts, their classes' info
    dicts, the circuit's dict (its plan, which holds the walk) and the
    limiter module's switches unchanged since it was ready, every element's
    `__dict__` the object seen and each hand-written limiter's twin its
    class's still (its `vlimit`, a class attribute) -- ready again.  Where
    the counter moved, the tables still standing as the last setup left
    them (`_walk_fast`) are ready and stamped again (2026-10-06: they were
    ready unstamped, the tuple checks on every call from the first move
    on)."""
    hc = _MOD['climit']
    if not (hc.WALK and hc.ENABLED):
        return False
    ep = _MOD['watch'].EPOCH
    if (ep is not None and rec.walk_stamp == ep.value and rec.walk_ok
            and tuple(map(_GETDICT, rec.walk_els)) == rec.walk_dicts
            and all(c.get('_c_limit') is k for c, k in rec.walk_hwk)):
        rec.walk_held += 1
        return True
    before = _MOD['watch'].now()
    w = hc._walk_for(tr.cir)
    if w is rec.walk:
        if rec.walk_ok is False:
            return False
        if rec.walk_ok and _walk_fast(rec):
            _walk_arm(rec, before)
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
        infos = {id(c): c for c in w.capable if isinstance(c, dict)}
        hw = [c.info for c in w.capable if type(c) is hc._Handwritten]
        rec.walk_hwk = tuple((c, k) for c, k in rec.walk_kern if type(c) is hc._Handwritten)
        rec.walk_watch = (*rec.walk_dicts, *infos.values(), *hw, tr.cir.__dict__, hc.__dict__)
        _walk_arm(rec, before)
    else:
        rec.walk_stamp = -1
    return ok


def _walk_arm(rec, before):
    """The walk's stamp, over the dicts its last setup read (`walk_watch`)."""
    W = _MOD['watch']
    rec.walk_arms, rec.walk_held = W.counted(rec.walk_arms, rec.walk_held), 0
    rec.walk_stamp = W.arm(rec.walk_watch, before) if rec.walk_arms <= W.MAX_ARMS else -1


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


def solve(tr, func, t, provided_function, seed, residual, want_j=True):
    """The multistep step's Newton solve, and its converged-point evaluation,
    in C: `(x, fj)` -- `x` the full-width solution as `_newton` returns it,
    `fj` `jacobian_only`'s `(None, J)` at it (`(None, None)`, and no G,
    where `want_j` is False: nothing reads it), or None where it must run
    again (a branch confirmation re-solved) -- or None: declined or handed
    back, nothing left behind, `_newton` runs."""
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
    ## (the parameter reads, as stamped where they stand: `_par_stamped`)
    P = _par_stamped(tr, td, par, td.get('_newton_c'))
    if P is None:
        before = M['watch'].now()
        nrsolver, scaler = par.nrsolver, getattr(par, 'scaler', None)
    else:
        nrsolver, scaler = P.nrsolver, P.scaler
    if getattr(tr, '_continuation_rescue', False) or nrsolver is not None or scaler is not None:
        return _no('newton_c:solver')
    ls = tr._get_linearsolver() if P is None else P.ls
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
    chord = bool(residual is not None and nrsolver is None
                 and tr._newton_option(par.chord_jacobian if P is None else P.chord,
                                       'chord_jacobian'))
    ## ---- from here the attempt leaves one trace: the step's source memo ----
    key = (t, par.analysis if P is None else P.analysis)
    try:
        had = key in memo
    except TypeError:
        return _no('newton_c:u')
    counts = {k: _PC.get(k, 0) for k in M['u_keys']}
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
    tol = tr._newton_tolerances() if P is None else P.tol
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
    if cached is not None and float((getattr(tr.epar, 'bypasstol', -1.0) if P is None
                                     else P.bypass) or -1.0) < 0.0:
        cxa, cCa = cached
        if (type(cxa) is np.ndarray and cxa.shape == (n,) and type(cCa) is np.ndarray
                and cCa.dtype == np.float64 and cCa.shape == (n, n)):
            np.copyto(buf['cx'], cxa)
            np.copyto(buf['cC'], cCa)
            s.have_cache = 1
    s.T, s.formula = T, formula
    s.a0, s.a1, s.a2, s.h, s.theta = a0, a1, a2, float(h), float(theta)
    if P is None:
        s.maxiter, s.reltol = int(par.maxiter), float(par.reltol)
    else:
        s.maxiter, s.reltol = int(P.maxiter), float(P.reltol)
    s.chord = 1 if chord else 0
    s.fold_j = 3 if want_j else 2
    if P is None:
        ## (every read above in full, every check passed: stamp them)
        _par_stamp(tr, td, par, rec, before, nrsolver, scaler, ls, tol)
    status = cfn(s)
    if status != _CODES['NEWTON_CONVERGED']:
        _undo_source(memo, key, had, counts)
        if status == _CODES['NEWTON_BAIL_WALKSTOP']:
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
    if not want_j:
        _PC['newton_c:j_unread'] += 1
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
    fj = (None, buf['Jj'].copy()) if want_j else (None, None)
    if tr._branch_on():
        target = stats if stats is not None else tr
        before = getattr(target, 'branch_screens', 0)
        tr._branch_after_solve(M['refnode_removed'](func, iref, tr.toolkit), x_res)
        if getattr(target, 'branch_screens', 0) != before:
            ## the confirmation re-solved: `jacobian_only` again, after it
            fj = None
    return x, fj
