"""THE NEWTON SOLVE IN C: a multistep step's Newton -- `_StepNewton._newton`
with `nrsolver.StandardNewton` or `ChordNewton` -- and the converged point's
`jacobian_only`, made in ONE C call (`csrc/newton.c`), bit for bit; pinned
against the Python in `tests/test_pinned_pairs.py`.  Called from
`Transient._solve_timestep` (`solve`); the C stepping loop (`_tran_loop_c`)
drives the same struct.

OWNS.  On the transient: `_newton_c` (its `_Ctx` record: the struct and
buffers, the LAPACK and BLAS addresses, the walk's and KLU's bindings, the
stamped parameter reads `par_rec`, the served mark `served`) and
`_newton_c_no` (a kept decline, `(stamp plan, _paths key)`: a property of
the circuit, `_KEPT`, answered at once while the circuit's dict holds that
plan).  On the circuit's limiter walk (`_hdl_climit._walk_for`): the tables
`limit_walk`'s setup fills (`F`, `PR`, `mirror`), and on each limiting
element its pack `_hdl_cp` where it has none.  On the `StructuralKLU` it
solves with: its `numerics` count.  `_paths`: `newton_c:<why>` a decline,
`newton_c:bail:<why>` a hand-back, and on a serve `newton_c:served`,
`newton_c:iterations`, `newton_c:jfold`, `newton_c:j_unread`,
`newton_c:sparse` -- with the Python's own `chord:chosen` where the Python
would count it.  Switch: `PYCIRCUIT_NEWTON_C=0` (`ENABLED`, read on every
call), and the evaluate core's (`PYCIRCUIT_TRAN_CORE`, `_tran_core.CORE`);
`STATUS` says why the path is off where it is.

ONE CALL.  Per iteration: the evaluate core's own C (`_tran_core.CORE_C`,
through its struct pointer); the reduced `J`, written column-major as
numpy's solve copies it, and `-F`; the LU -- numpy's own `scipy_dgesv_64_`
for the plain Newton, SciPy's `scipy_dgetrf_`/`scipy_dgetrs_` for the dense
chord (SciPy's OpenBLAS, as the Python's chord factors), KLU on the run's
own symbolic analysis for the sparse chord (`_klu_bind`); `x + dx`; the
limiter walk's own driver (`_hdl_climit._WALK_C`) and `dx = x_next - x`
re-taken where a limiter exists; the convergence test in `nrsolver`'s order
(`|J||x|` through numpy's own `scipy_cblas_dgemv64_`) and the chord's
shrink test.  Then, in the same call, the converged point's evaluation:
`jacobian_only`'s G, C and q -- C and q alone where nothing reads the
step's Jacobian (`want_j` False: a fixed multistep run).

TRANSACTIONAL.  Nothing is written on the transient or the devices before
the C has converged; anything unusual -- a LAPACK info, a non-finite value,
the walk stopping, a chord break, maxiter -- rolls back the attempt's one
trace (the step's source-memo entry and its counts), and `_newton` repeats
the same arithmetic from the same seed: the same answer, exception, warning
and statistics.

THE SERVED MARK.  `_newton_c.served` is the watch counter at a call this
path served (with no branch confirmation re-solving after it), None after
any other outcome: the C stepping loop enters only while it equals the
counter now.

DECLINES, before anything is touched (`_paths.no`, `newton_c:<why>`): the
switch (`off`) or the evaluate core's (`core`); no build or no OpenBLAS
(`driver`); a piece of the chain the C stands in for patched or subclassed
(`patched`, `class`: `_chain`); a toolkit that is not numeric (`toolkit`);
an instance shadow of the circuit's passes or of the transient's machinery
(`shadow_cir`, `shadow_tr`); a `provided_function` (`pf`); no step source
memo (`umemo`); a stateful limiter (`stateful`); the continuation rescue or
a caller's solver or scaler (`solver`); a linear solver that is not
`AutoSolver`'s dense or `StructuralKLU` choice (`linsolver`), a KLU not
analysed yet or whose pattern does not hold the core's (`klu`); a circuit the core cannot
serve (`unservable`); a reference row out of range (`iref`); a seed that is
not one finite float64 vector of the circuit's size (`x`); an integrator
outside the core's enum (`formula`); a history row of the wrong shape
(`history`); a limiting element that is not its C kernel (`limit`; kept
where its class has no C limiter at all); the core's readiness (`ready`:
`_Core.probe`, counted once, by the path that follows); a source pass that
raises or gives no float64 vector of the circuit's size (`u`: the
attempt's trace undone).  Between two iterations only the C runs, so
readiness checked once a solve holds for every iteration.  Between two
SOLVES the walk's readiness, the core's and the parameter reads (the solver
options, the linear solver, the tolerances) hold as stamped (`_watch`):
unchecked while no dict they read has changed, what no dict watcher sees
compared on every call.

THE BRANCH CHECK follows as in `_newton`; its screen finds `C` at the
converged point in `_C_cache`.  Where it fires, the confirmation's
speculative solves overwrite the step's state, and `jacobian_only` runs
again after it, as in the Python.

History: `doc/transient_history.md`, `_tran_newton_c.py`.
"""
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

#: newton.h's named codes: the C's results (`NEWTON_CONVERGED`, the bails
#: `NEWTON_BAIL_<WHY>`)
_CODES = _cbuild.constants(NEWTON_CDEF)
#: {code: why} of the C's bails: why it handed the solve back
BAIL = _cbuild.bails(_CODES, 'NEWTON_BAIL_')
#: the counts the step's source pass makes (`_source_at`, `cir.u`; the
#: source plan's own, `_stamp_plan._SK`, are added at `_mods`): the
#: attempt's one trace, rolled back with the memo entry when it bails
_U_KEYS = ('umemo:hit', 'umemo:miss', 'umemo:unhashable', 'u:called')

_driver = None
#: the modules and names `solve` reads (`_mods`; `_tran_radau_c` and
#: `_tran_radau_tc` read them too) and the chain it checks (`_chain`),
#: each filled once
_MOD = {}
#: a parameter not stored (`_ParRec.bp_obj`)
_ABSENT = object()
_GETDICT = operator.attrgetter('__dict__')
#: float64 (an exact array's dtype is this object: `is`, then `==`)
_F64 = np.dtype(np.float64)
#: the codes the core and the units share (pcx_abi.h)
_ABI = _cbuild.constants(_cbuild.csrc('pcx_abi.h'))
#: the core's want bit for "these matrices are zero off the pattern"
_WANT_PATTERN = _ABI['PCX_WANT_PATTERN']
#: the passes evaluated at the converged point (`fold_j`): `G` and `C` where an
#: attempt reads the step's Jacobian, `C` alone where none does
FOLD_J = _ABI['PCX_PASS_G'] | _ABI['PCX_PASS_C']
FOLD_C = _ABI['PCX_PASS_C']
_GETPACK = operator.methodcaller('get', '_hdl_cp')
#: the declines that are a property of the circuit: kept with its stamp
#: plan, so the next step declines at its first check
_KEPT = ('unservable', 'stateful', 'limit')


def driver():
    """`(ffi, cfn, LAPACK addresses)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            import scipy
            import scipy.linalg  # (its OpenBLAS mapped: the chord's)
            npl, spl = _cbuild.openblas(np), _cbuild.openblas(scipy)
            if npl is None or spl is None:
                STATUS = "off (numpy's or SciPy's own OpenBLAS not found)"
                return None
            addr = _cbuild.address
            lapack = {'dgesv': addr(npl.scipy_dgesv_64_),
                      'dgemv': addr(npl.scipy_cblas_dgemv64_),
                      ## (the early abort's back-substitution on dgesv's LU)
                      'dgetrs64': addr(npl.scipy_dgetrs_64_),
                      'dgetrf': addr(spl.scipy_dgetrf_), 'dgetrs': addr(spl.scipy_dgetrs_)}
        except (ImportError, OSError, AttributeError) as e:
            STATUS = f'off ({e})'
            return None
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(NEWTON_C, NEWTON_CDEF, 'pcx_newton')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, lapack)
        STATUS = 'c'
    return _driver or None


def _chain():
    """The machinery `solve` stands in for, as defined -- the transient's
    class and its Newton methods, `nrsolver`'s two loops and the chord's
    settled test, the reference-row helpers `_newton` calls, the evaluate
    core's Python entry -- read once a call finds every piece its module's
    own (`_paths.genuine`); None: not yet.  A call then compares what the
    Python path would call now."""
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
               (nrsolver._chord_settled, '_chord_settled', 'pycircuit.circuit.nrsolver'),
               (_tran_core.evaluate, 'evaluate', 'pycircuit.circuit._tran_core'))
        if not all(_paths.genuine(*o) for o in own):
            return None
        _MOD['chain'] = tuple(o for o, _q, _m in own)
        _MOD['chain_mods'] = (_tran_newton, nrsolver.StandardNewton, nrsolver.ChordNewton,
                              _tran_core, nrsolver)
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
            _tran_lte_c,
            _watch,
        )
        from pycircuit.circuit.integrator import EulerIntegrator
        from pycircuit.circuit._lte_kernels import bdf2_alphas
        from pycircuit.circuit.dcanalysis import refnode_removed
        from pycircuit.circuit.linearsolver import AutoSolver, StructuralKLU
        from pycircuit.circuit.toolkit import NumericToolkit
        _MOD.update(climit=_hdl_climit, limiting=_limiting, companion=_tran_companion, watch=_watch,
                    plan_for=_stamp_plan._plan_for,
                    core=_tran_core, bdf2_alphas=bdf2_alphas,
                    refnode_removed=refnode_removed, AutoSolver=AutoSolver,
                    StructuralKLU=StructuralKLU,
                    NumericToolkit=NumericToolkit,
                    lte_c=_tran_lte_c, Euler=EulerIntegrator,
                    u_keys=_U_KEYS + tuple(_stamp_plan._SK['u'].values())
                    + tuple(_stamp_plan._SD.values()))
    return _MOD


class _Ctx:
    """One transient's struct (`s`, a `newton_t`) and buffers (`buf`) for
    one core, size and reference row, set once: a solve copies its inputs
    in and its state out.  The struct's other pointers: the core's struct
    and entry, the LAPACK and BLAS addresses, the walk's tables
    (`_walk_ready`), KLU's arrays and objects (`_klu_bind`), and the cached
    `C` read where it stands from `_cbuild.INPLACE_MIN_ENTRIES` entries
    (`cC_inplace`; else copied into `cC_buf`).  Held as `tr._newton_c`; the
    C loop (`_tran_loop_c`) points its struct at `s` and reads `core`, `klu`
    and `served`; `_tran_radau_c` and `_tran_radau_tc` call `_walk_ready`
    with records of their own whose walk fields are named as here."""

    __slots__ = (
        'buf',
        'cC_buf',
        'cC_inplace',
        'core',
        'ffi',
        'ipiv32',
        'ipiv64',
        'iref',
        'keep',
        'klu',
        'klu_keep',
        'klu_sym',
        'n',
        'par_arms',
        'par_held',
        'par_rec',
        's',
        'served',
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
        ## (the n x n buffers start at zero: where the core writes only the
        ## circuit's pattern -- `s.pattern` -- the entries off it stay so)
        Z = np.zeros
        self.buf = buf = {
            'xf': np.empty(n), 'xnf': np.empty(n), 'x0f': np.empty(n),
            'A': np.empty(m * m), 'absJ': Z(m * m),
            'b': np.empty(m), 'Fr': np.empty(m), 'x': np.empty(m), 'xn': np.empty(m),
            'xd': np.empty(m), 'ax': np.empty(m), 'Is': np.empty(m),
            'F0': np.empty(n), 'F': np.empty(n), 'J0': Z(nn), 'J': Z(nn),
            'C0': Z(nn), 'q0': np.empty(n), 'iq0': np.empty(n), 'Geq0': Z(nn),
            'C': Z(nn), 'q': np.empty(n), 'iq': np.empty(n), 'Geq': Z(nn),
            'q1': np.empty(n), 'q2': np.empty(n), 'iq1': np.empty(n), 'u': np.empty(n),
            'abstol': np.empty(m), 'xtol': np.empty(m), 'xseed': np.empty(n),
            'cx': np.empty(n), 'cC': Z(nn),
            'xc': np.empty(n), 'Cj': Z(nn), 'qj': np.empty(n), 'iqj': np.empty(n),
            'Geqj': Z(nn), 'Jj': Z(nn),
        }
        self.ipiv64 = np.empty(max(m, 1), dtype=np.int64)
        self.ipiv32 = np.empty(max(m, 1), dtype=np.int32)
        fb = ffi.from_buffer
        self.keep = []
        for k, arr in buf.items():
            h = fb('double *', arr)
            self.keep.append(h)
            setattr(s, k, h)
        ## (the cached C's own buffer, where it is copied: below
        ## `_cbuild.INPLACE_MIN_ENTRIES`, or not C-contiguous)
        self.cC_buf = s.cC
        self.cC_inplace = n * n >= _cbuild.INPLACE_MIN_ENTRIES
        self.keep.append(fb('int64_t *', self.ipiv64))
        s.ipiv64 = self.keep[-1]
        self.keep.append(fb('int *', self.ipiv32))
        s.ipiv32 = self.keep[-1]
        for k in ('dgesv', 'dgetrf', 'dgetrs', 'dgemv', 'dgetrs64'):
            setattr(s, k, ffi.cast('void *', lapack[k]))
        s.core = ffi.cast('void *', int(cffi_core.cast('uintptr_t', core.cs)))
        s.core_fn = ffi.cast('void *', int(cffi_core.cast('uintptr_t', cfn_core)))
        s.n, s.iref = n, iref
        s.pattern = _WANT_PATTERN if core.npat else 0
        s.walk_fn = ffi.NULL
        s.wn = 0
        s.klu_sym = ffi.NULL
        self.klu = self.klu_sym = None
        self.klu_keep = ()
        self.tol = None
        self.walk = self.walk_ok = None
        self.walk_els = self.walk_dicts = self.walk_packs = self.walk_kern = ()
        self.walk_hw_els = self.walk_hwk = ()
        self.walk_stamp, self.walk_arms, self.walk_held = -1, 0, 0
        self.walk_watch = ()
        self.par_rec, self.par_arms, self.par_held = None, 0, 0
        ## (the watch counter at this record's last serve, None after any
        ## other outcome: the stepping loop in C enters only on a serve
        ## nothing watched has moved since -- `_tran_loop_c`)
        self.served = None


def _klu_bind(rec, sk):
    """Point the C's sparse solves at `sk`, the run's `StructuralKLU`: its
    pattern's places in the core's full J (`kat`) and in the reduced |J|
    (`kabs`), its analysis and its common, the library's entries, the full
    Newton's validation (its buffers, `REFACTOR_RESIDUAL_TOL`).  False -- declined --
    where it is not analysed yet (its first solve does that) or where the
    core's pattern does not lie within its own: there a nonzero off its
    pattern would go to its fallback in the Python."""
    n, iref = rec.n, rec.iref
    m = n - 1
    if sk.m != m or not sk._symbolic:
        return False
    r, c = np.divmod(sk.take, m)
    kat = np.ascontiguousarray((r + (r >= iref)) * n + (c + (c >= iref)), dtype=np.int64)
    pat = rec.core.arrays['pat']
    pr, pc = np.divmod(pat, n)
    if not np.isin(pat[(pr != iref) & (pc != iref)], kat).all():
        return False
    ffi, s = rec.ffi, rec.s
    lib, ct = sk._lib, sk._ct
    addr = _cbuild.address
    kabs = np.ascontiguousarray(sk.take, dtype=np.int64)
    kAx = np.empty(max(sk.nnz, 1))
    ## (the sparse full Newton's: J reduced, zero off the pattern, and its
    ## validation's two vectors -- newton.h)
    R, rb, ry = np.zeros(m * m), np.empty(m), np.empty(m)
    keep = (sk, kat, kabs, kAx, ffi.from_buffer('int *', sk.Ap), ffi.from_buffer('int *', sk.Ai),
            ffi.from_buffer('int64_t *', kat), ffi.from_buffer('int64_t *', kabs),
            ffi.from_buffer('double *', kAx), R, rb, ry, ffi.from_buffer('double *', R),
            ffi.from_buffer('double *', rb), ffi.from_buffer('double *', ry))
    s.knnz = sk.nnz
    s.kAp, s.kAi, s.kat, s.kabs, s.kAx = keep[4:9]
    s.R, s.rb, s.ry = keep[12:]
    s.klu_factor, s.klu_solve, s.klu_free, s.klu_refactor = (ffi.cast('void *', addr(f)) for f in (
        lib.klu_factor, lib.klu_solve, lib.klu_free_numeric, lib.klu_refactor))
    s.resid_tol = float(sk.REFACTOR_RESIDUAL_TOL)
    s.klu_common = ffi.cast('void *', ct.addressof(sk._common))
    s.klu_sym = ffi.cast('void *', sk._symbolic)
    rec.klu, rec.klu_sym, rec.klu_keep = sk, sk._symbolic, keep
    return True


class _ParRec:
    """`solve`'s reads of the analysis parameters and of the transient's
    own caches, stamped (`_watch`).  Watched: `tr.par`'s dicts, and
    `tr.epar`'s instance dict and parameter table -- not its values dict:
    the analyses write `analysis_kind` there around every solve's operating
    point (`analysis.analysis_kind`), which would move every stamp's
    counter twice a solve.  The one value read there, `bypasstol`, is
    compared by its object on every call (`bp_obj`); the methods that read
    the parameters and the caches they keep in the transient's dict are
    compared on every call too (`_par_stamped`)."""

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
    them (`_walk_fast`) are ready and stamped again.  `rec` is a `_Ctx`, or
    another unit's record with the same walk fields (`s`, `ffi`, `walk*`:
    `_tran_radau_c`, `_tran_radau_tc`)."""
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


def _undo_source(memo, keys, had, counts):
    """An attempt's one trace, removed (the Newton's, the Radau units'): the
    source memo's entries its calls made -- each key it had not before --
    and the counts they added; the Python makes them again."""
    if memo is not None:
        for key, h in zip(keys, had):
            if not h:
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
    back, nothing left behind, `_newton` runs.

    In order: the checks that decline (`_admit`; module docstring:
    DECLINES); the step's companion (`_companion`); the record and its
    bindings (`_record`); the source pass, the attempt's one trace
    (`_source`); the inputs copied in (`_load`), the parameter reads
    stamped, the call -- a bail undoes the trace (`_undo_source`); the
    serve (`_serve`).  A serve writes the step's state (`_q_cache`,
    `_C_cache`, `_C_pattern` where the core has a pattern, `_iq`, `_Geq`,
    `_Cmat`, `active_integrator`, `_companion_coeffs`, `_effective_method`),
    the statistics' Newton iterations and the counts, runs the branch
    check, and sets the record's `served`."""
    td = tr.__dict__
    last = td.get('_newton_c')
    if last is not None:
        last.served = None
    if not ENABLED:
        return _no('newton_c:off')
    adm = _admit(tr, td, provided_function, seed)
    if adm is None:
        return None    # paths: counted by _admit
    drv, M, P, before, nrsolver, scaler, ls, ch, core, memo = adm
    n, iref = core.n, tr.irefnode
    comp = _companion(tr, M, n)
    if comp is None:
        return None    # paths: counted by _companion
    h, h_last, active, formula, rows, coeffs = comp
    got = _record(tr, td, M, core, drv, P, nrsolver, residual, ls, ch)
    if got is None:
        return None    # paths: counted by _record
    rec, T, chord, sparse = got
    early = _arm(tr, td, M, chord, sparse, active)
    if early is False:
        return _no('newton_c:early')
    ## from here the attempt leaves one trace: the step's source memo
    src = _source(tr, M, P, memo, t, n)
    if src is None:
        return None    # paths: counted by _source
    u, trace = src
    tol = _load(tr, td, rec, P, rows, u, seed)
    s = rec.s
    s.T, s.formula = T, formula
    a0, a1, a2, theta = coeffs
    s.a0, s.a1, s.a2, s.h, s.theta = a0, a1, a2, float(h), float(theta)
    par = tr.par
    if P is None:
        s.maxiter, s.reltol = int(par.maxiter), float(par.reltol)
    else:
        s.maxiter, s.reltol = int(P.maxiter), float(P.reltol)
    s.chord = 1 if chord else 0
    s.fold_j = FOLD_J if want_j else FOLD_C
    if P is None:
        ## (every read above in full, every check passed: stamp them)
        _par_stamp(tr, td, par, rec, before, nrsolver, scaler, ls, tol)
    _arm_set(rec, M, early)
    status = drv[1](s)
    s.early_lte = rec.ffi.NULL
    if status == _CODES['NEWTON_EARLY']:
        _early(tr, rec, early, chord, sparse, ch)
    if early is not None and status == _CODES['NEWTON_CONVERGED']:
        _PC['early:checks'] += s.early_checks
    if status != _CODES['NEWTON_CONVERGED']:
        key, had, cvals = trace
        _undo_source(memo, (key,), (had,), dict(zip(M['u_keys'], cvals)))
        if status == _CODES['NEWTON_BAIL_WALKSTOP']:
            ## (the walk stopped: its tables are taken again next time)
            rec.walk_ok = None
        return _no('newton_c:bail:' + BAIL.get(status, str(status)))
    return _serve(tr, func, M, rec, core, iref, active, h, h_last, chord, sparse, ch, want_j)


def _admit(tr, td, provided_function, seed):
    """The checks that decline, in order (module docstring: DECLINES) --
    `(the driver, the modules, the stamped reads or None, the watch counter
    before them, nrsolver, scaler, the linear solver, its choice, the core,
    the source memo)`, or None: declined (counted)."""
    cd = tr.cir.__dict__
    kept = td.get('_newton_c_no')
    if kept is not None:
        ## (the circuit's own reason, under the plan its dict holds -- a stale
        ## plan is replaced at the next evaluation, which follows: until then
        ## the decline runs the Python Newton, never another answer)
        if cd.get('_stamp_plan') is kept[0]:
            return _no(kept[1])
        del td['_newton_c_no']
    drv = _driver or driver()
    if not drv:
        return _no('newton_c:driver')
    M = _MOD or _mods()
    ## (the machinery the C stands in for, as defined: a subclass, or a
    ## piece patched on its class or module, takes the Python Newton)
    chain = M.get('chain') or _chain()
    if chain is None:
        return _no('newton_c:patched')
    T = chain[0]
    if type(tr) is not T:
        return _no('newton_c:class')
    tn, SN, CN, tcm, nrs = M['chain_mods']
    if (T._newton, T._newton_limiter, T._residual_and_jacobian, T._get_nrsolver,
            T._get_scaler, tn.refnode_removed, tn.remove_row_col, SN.solve_system,
            CN.solve_system, nrs._chord_settled, tcm.evaluate) != chain[1:]:
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
    before = None
    if P is None:
        before = M['watch'].now()
        nrsolver, scaler = par.nrsolver, getattr(par, 'scaler', None)
    else:
        nrsolver, scaler = P.nrsolver, P.scaler
    if getattr(tr, '_continuation_rescue', False) or nrsolver is not None or scaler is not None:
        return _no('newton_c:solver')
    ls = tr._get_linearsolver() if P is None else P.ls
    ## (the dense choice -- LAPACK, as the Python's -- or KLU on the circuit's
    ## pattern: the sparse chord, `_klu_bind`)
    ch = ls._choice if type(ls) is M['AutoSolver'] else None
    if ch is None or not (ch is ls._dense or type(ch) is M['StructuralKLU']):
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
    if not (type(seed) is np.ndarray and (seed.dtype is _F64 or seed.dtype == np.float64) and seed.ndim == 1
            and seed.shape[0] == n and np.isfinite(seed).all()):
        return _no('newton_c:x')
    return drv, M, P, before, nrsolver, scaler, ls, ch, core, memo


def _companion(tr, M, n):
    """The step's companion, set up as `get_diff` sets it up -- the same
    lines as `_tran_core._companion`'s, kept twice on purpose (why:
    `doc/transient_history.md`, `_tran_core.evaluate`; a change here is a
    change there): `(h, h_last, the active integrator, its formula, (q1,
    q2, iq1), (a0, a1, a2, theta))`, or None: declined (counted)."""
    tc = M['core']
    h = tr._dt
    h_last = tr._dt_last if tr._dt_last is not None else h
    active = tr.base_integrator.check_order_drop(h, h_last, tr._is_first_step)
    formula = tc.FORMULA.get(type(active).__name__)
    if formula is None:
        return _no('newton_c:formula')
    q1 = tc._row(tr._qlast[0], n)
    q2 = tc._row(tr._qlast[1], n) if formula == tc._GEAR2 else q1
    iq1 = tc._row(tr._iqlast[0], n) if formula in tc._IQ_FORMS else q1
    if q1 is None or q2 is None or iq1 is None:
        return _no('newton_c:history')
    a0 = a1 = a2 = theta = 0.0
    if formula == tc._GEAR2:
        a0, a1, a2 = M['bdf2_alphas'](h, h_last)
    elif formula == tc._THETA:
        theta = active.theta_at(h)
    return h, h_last, active, formula, (q1, q2, iq1), (a0, a1, a2, theta)


def _record(tr, td, M, core, drv, P, nrsolver, residual, ls, ch):
    """The record (`_Ctx`, made again for another core, size or reference
    row) and its bindings: the walk where anything limits, the core's
    readiness, the chord, KLU -- `(the record, the temperature, chord,
    sparse)`, or None: declined (counted)."""
    ffi, _cfn, lapack = drv
    n, iref = core.n, tr.irefnode
    ## (whether anything limits: what the genuine `_newton_limiter` -- the
    ## chain `_admit` checks holds it so -- answers by its memo, read without
    ## building the closure it returns; the method only where the memo is
    ## unset.  A limiting element whose class has no C limiter at all --
    ## PSP's -- is a property of the circuit: `_walk_ready` finds it, kept)
    _has = getattr(tr, '_has_limiters', None)
    has_limiter = bool(_has) if _has is not None else tr._newton_limiter() is not None
    rec = td.get('_newton_c')
    if rec is None or rec.core is not core or rec.n != n or rec.iref != iref:
        cffi_core, cfn_core, _dg = M['core'].driver()
        rec = _Ctx(core, n, iref, ffi, cfn_core, cffi_core, lapack)
        td['_newton_c'] = rec
    s = rec.s
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
    par = tr.par
    chord = bool(residual is not None and nrsolver is None
                 and tr._newton_option(par.chord_jacobian if P is None else P.chord,
                                       'chord_jacobian'))
    sparse = ch is not ls._dense
    if sparse:
        ## (the full Newton's sparse solve refactors the run's own numeric,
        ## validated by its residual: newton.h's THE SPARSE FULL NEWTON)
        if (rec.klu is not ch or ch._symbolic != rec.klu_sym) and not _klu_bind(rec, ch):
            return _no('newton_c:klu')
        ## (the run's numeric as it stands: the Python's solves replace it)
        num = ch._numeric
        s.klu_num = ffi.cast('void *', num) if num else ffi.NULL
    elif rec.klu is not None:
        s.klu_sym = ffi.NULL
        rec.klu = rec.klu_sym = None
        rec.klu_keep = ()
    return rec, T, chord, sparse


def _source(tr, M, P, memo, t, n):
    """The source pass at the step's end, the attempt's one trace in the
    source memo `memo`: `(u, (its key, whether the memo held it, the counts
    before))` -- what `_undo_source` takes back -- or None: declined
    (counted), the trace undone."""
    key = (t, tr.par.analysis if P is None else P.analysis)
    try:
        had = key in memo
    except TypeError:
        return _no('newton_c:u')
    ## (the counts before the pass, as values -- their dict made only to
    ## roll back, which a served attempt never does)
    u_keys = M['u_keys']
    cvals = tuple(map(_PC.get, u_keys, repeat(0)))
    try:
        u = tr._source_at(t, None)
    except Exception:                                          # noqa: BLE001
        ## (the pass raises again in `_newton`, where it raised before)
        _undo_source(memo, (key,), (had,), dict(zip(u_keys, cvals)))
        return _no('newton_c:u')
    if not (type(u) is np.ndarray and (u.dtype is _F64 or u.dtype == np.float64) and u.ndim == 1
            and u.shape[0] == n):
        _undo_source(memo, (key,), (had,), dict(zip(u_keys, cvals)))
        return _no('newton_c:u')
    return u, (key, had, cvals)


def _load(tr, td, rec, P, rows, u, seed):
    """The record's inputs: the tolerances (copied when they move), the
    history rows, u, the seed, and `_C_lookup`'s cached state at the seed
    (the multistep memo is empty; no stateful limiter -- declined -- and no
    bypass: its other routes).  Returns the tolerances."""
    s, buf = rec.s, rec.buf
    tol = tr._newton_tolerances() if P is None else P.tol
    if tol is not rec.tol:
        rec.tol = tol
        buf['abstol'][...] = tol[2]
        buf['xtol'][...] = tol[3]
    q1, q2, iq1 = rows
    np.copyto(buf['q1'], q1)
    np.copyto(buf['q2'], q2)
    np.copyto(buf['iq1'], iq1)
    np.copyto(buf['u'], u)
    np.copyto(buf['xseed'], seed)
    cached = td.get('_C_cache')
    s.have_cache = 0
    if cached is not None and float((getattr(tr.epar, 'bypasstol', -1.0) if P is None
                                     else P.bypass) or -1.0) < 0.0:
        cxa, cCa = cached
        n = rec.n
        if (type(cxa) is np.ndarray and cxa.shape == (n,) and type(cCa) is np.ndarray
                and (cCa.dtype is _F64 or cCa.dtype == np.float64) and cCa.shape == (n, n)):
            np.copyto(buf['cx'], cxa)
            ## (the cached C read where it stands -- the C only reads it: an
            ## n x n copy a call saved -- else copied, `INPLACE_MIN_ENTRIES`)
            if rec.cC_inplace and cCa.flags.c_contiguous:
                s.cC = rec.ffi.from_buffer('double *', cCa)
            else:
                np.copyto(buf['cC'], cCa)
                if rec.cC_inplace:
                    s.cC = rec.cC_buf
            s.have_cache = 1
    return tol


def _arm(tr, td, M, chord, sparse, active):
    """THE EARLY ABORT (`Transient.early_abort`): an attempt the loop armed
    (`tr._early_arm`) is checked at its first iterate by newton.c, on an
    error-test context filled for it (`_tran_lte_c.early_context`) -- not on
    the full Newton of a sparse matrix or a 2nd-order step with no third
    point, where `_newton` makes no estimate either.  `(the context, the
    attempt's inputs, euler, the threshold)`, None: not armed, or False: the
    context declined -- the Python Newton makes the check."""
    arm = td.get('_early_arm')
    if (arm is None or not (chord or not sparse)
            or (tr._dt_last2 is None and type(active) is not M['Euler'])):
        return None
    got = M['lte_c'].early_context(tr, arm, active)
    return False if got is None else (*got, float(arm.th))


def _arm_set(rec, M, early):
    """The Newton record armed with `_arm`'s context, or not (`early_lte`
    NULL): the context's struct and lte.c's entry, by address -- another
    unit's cffi made them."""
    s = rec.s
    if early is None:
        s.early_lte = rec.ffi.NULL
        return
    lffi, lcfn, _dgesv = M['lte_c'].driver()
    s.early_lte = rec.ffi.cast('void *', int(lffi.cast('uintptr_t', early[0].s)))
    s.early_lte_fn = rec.ffi.cast('void *', int(lffi.cast('uintptr_t', lcfn)))
    s.early_th = early[3]


def _early(tr, rec, early, chord, sparse, ch):
    """The early abort's stop (NEWTON_EARLY): the counts the Python Newton
    makes up to its abort -- the chord chosen, the sparse chord's numeric,
    the check and the abort -- and the statistics' iteration (the trace is
    the attempt's: nothing repeats it), then the abort `_newton` raises:
    the estimate, the attempt's error-test inputs, its exponent's
    denominator."""
    from pycircuit.circuit._tran_newton import _EarlyAbort
    s = rec.s
    iters = int(s.iters)
    if chord:
        _PC['chord:chosen'] += 1
    if sparse:
        ch.numerics += 1
    _PC['newton_c:iterations'] += iters
    _PC['early:checks'] += 1
    _PC['early:aborts'] += 1
    stats = getattr(tr, 'statistics', None)
    if stats is not None:
        stats.newton_iterations += iters
    _ctx, inputs, euler, _th = early
    raise _EarlyAbort(float(s.early_err), inputs, 2.0 if euler else 3.0)


def _serve(tr, func, M, rec, core, iref, active, h, h_last, chord, sparse, ch, want_j):
    """Converged: the counts, the statistics' iterations, the converged
    point's state written as `_newton` and `get_diff` leave it, the branch
    check -- `(x, fj)`, `fj` None where the confirmation re-solved (the
    record's `served` set otherwise)."""
    s, buf = rec.s, rec.buf
    iters = int(s.iters)
    if chord:
        _PC['chord:chosen'] += 1
    if sparse:
        if chord:
            ## (the numeric `StructuralKLU.factor` would have made: counted)
            ch.numerics += 1
        else:
            ## (the full Newton's solves: a refactor of the run's numeric each)
            ch.refactors += iters
        _PC['newton_c:sparse'] += 1
    _PC['newton_c:served'] += 1
    _PC['newton_c:iterations'] += iters
    _PC['newton_c:jfold'] += 1
    if not want_j:
        _PC['newton_c:j_unread'] += 1
    stats = getattr(tr, 'statistics', None)
    if stats is not None:
        stats.newton_iterations += iters
    x_res = buf['x'].copy()
    ## (the full-width solution the C wrote at the converged point, the
    ## reference row's +0.0 in its place: `insert_row`'s array, byte for
    ## byte -- every serve evaluates the converged point, `fold_j` never 0)
    x = buf['xc'].copy()
    Cj = buf['Cj'].copy()
    tr._q_cache = (x, buf['qj'].copy())
    tr._C_cache = (x, Cj)
    if core.npat:
        ## (written on the core's pattern: the branch screen reads only it)
        tr._C_pattern = (Cj, core.arrays['pat'])
    ## (the step's companion left as `get_diff` leaves it -- the same writes
    ## as `_tran_core._leave`'s; `transient._RUN_STATE`'s 'step' block)
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
    if fj is not None:
        rec.served = M['watch'].now()
    return x, fj
