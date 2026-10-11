"""RADAU'S DENSE STAGE NEWTON IN C: a Radau IIA(3) step's coupled ``3m``
Newton (`_coupled_stage_solver`'s `_stage_newton`, undamped and without a
shunt) in ONE C call, bit for bit the Python Newton (pinned in
`tests/test_pinned_pairs.py`).

OWNS.  On the transient, `tr._radau_c` (`_Ctx`: the `radau_t` struct, its
buffers for `MAXA` assemblies, and the limiter walk's fields, which
`_tran_newton_c._walk_ready` sets; declared in `transient._RUN_STATE`).  On
a serve it writes what the Python Newton leaves: the device memo
(`tr._memo_put`), the step's source memo (`tr._u_memo`) and its count
`umemo:hit`.  `_paths`: `radau_c:` -- `radau_c:<why>` a decline,
`radau_c:bail:<why>` a hand-back (`BAIL`), `radau_c:served`,
`radau_c:assemblies`.  Switch: `PYCIRCUIT_RADAU_C=0` (`ENABLED`, read on
every call); `STATUS` says why the path is off where it is.

THE NEWTON.  Per iteration the Python makes the block system
(`_coupled_stage_system`: nine blocks, each through `remove_row_col`), the
three stage evaluations, the limiter per stage, numpy's solve and the
convergence test.  The C makes, per assembly, the three stages' passes
through the evaluate core's own C (`_tran_core.CORE_C`, formula 4, as
`_tran_core.passes`), ``K_j = -(i_j + u_j)``, ``F_i = (q_i - q_n) - h
(((0 + A_i0 K_0) + A_i1 K_1) + A_i2 K_2)`` and the blocks ``C_i + (h
A_ij) G_j`` (``(h A_ij) G_j`` off the diagonal) over the full width, then
reduced and written column-major as numpy's solve copies them; numpy's own
`scipy_dgesv_64_` on ``-R``; per stage ``Y_prev + insert(dY_i, iref,
0.0)`` (the reference row's ``+ 0.0``), the limiter walk's own driver
(`_hdl_climit._WALK_C`) against ``Y_prev`` and the limited step's largest
reduced entry; the convergence test (`_stages_converged`: the step within
``reltol`` of the largest reduced stage plus ``abstol``) after the
assembly at the new stages, as the Python assembles before it tests.

TRANSACTIONAL.  The C writes nothing outside its buffers.  Where it
converges, what the Python Newton leaves is made here: every assembly's
three stages recorded in the device memo (`_memo_put`, in the Python's
order: `_finish_stage_step`, the next step's ``q_n``, the transform path
and PSS read them), the step's source memo holding the three stage times
(their first assembly's calls are made here, before the C) and its hit
count for every later assembly's.  Anything unusual -- a non-finite value,
a floating-point exception in its own arithmetic (numpy would warn), a
residual too large for its ``sum |R|`` (the undamped Newton's only use of
it), a LAPACK info, the walk stopping, more assemblies than its buffers
hold, maxiter -- hands the solve back: the source memo's entries and counts
rolled back, and the Python Newton runs from the same seed (the same
answer, exception, warning).

DECLINES, before anything is touched (`_paths.no`): the switch (and the
evaluate core's two), no build, a toolkit that is not numeric, a subclass
of the transient or a patched piece of the Newton (`_chain`), an instance
shadow of the circuit's passes or limiter or of the transient's pieces, a
`provided_function`, a stateful limiter or circuit-level limiting, a bypass,
no source memo, a circuit the core cannot serve, a limiting element that
is not its C kernel, a batch not ready, inputs of another type or size.

History: `doc/transient_history.md`, `_tran_radau_c.py`.
"""
import os
from itertools import repeat

import numpy as np

from pycircuit.circuit import _cbuild, _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_RADAU_C', '1') != '0'
STATUS = 'not loaded'
#: the assemblies one call holds (a Newton of more iterations is handed back)
MAXA = 8

RADAU_C = _cbuild.csrc('radau.c')
RADAU_CDEF = _cbuild.csrc('radau.h')

#: the C's results, named in radau.h (`_CODES`); `BAIL`: {code: why} of
#: the codes past 1 (converged) -- why it handed the solve back
_CODES = _cbuild.constants(RADAU_CDEF)
BAIL = _cbuild.bails(_CODES, 'RADAU_BAIL_')

_driver = None
_MOD = {}
#: numpy's solve as imported: a caller's stand-in (a spy, an injected
#: failure) takes the Python Newton
_SOLVE = np.linalg.solve


def driver():
    """`(ffi, cfn, dgesv address)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        from pycircuit.circuit import _tran_newton_c
        nd = _tran_newton_c.driver()
        if nd is None:
            STATUS = f'off (the C Newton: {_tran_newton_c.STATUS})'
            return None
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(RADAU_C, RADAU_CDEF, 'pcx_radau')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, nd[2]['dgesv'])
        STATUS = 'c'
    return _driver or None


def _chain():
    """The machinery `solve` stands in for, as defined (`_paths.genuine`),
    read once every piece is its module's own; None: not yet.  A call then
    compares what the Python Newton would call now (`_now`)."""
    if 'chain' not in _MOD:
        from pycircuit.circuit import _hdl_climit, _tran_core, _tran_radau
        from pycircuit.circuit.circuit import SubCircuit
        from pycircuit.circuit.transient import Transient as T
        R, CM = 'pycircuit.circuit._tran_radau', 'pycircuit.circuit._tran_companion'
        own = ((T, 'Transient', 'pycircuit.circuit.transient'),
               (T._coupled_stage_context, '_RadauStages._coupled_stage_context', R),
               (T._coupled_stage_system, '_RadauStages._coupled_stage_system', R),
               (T._stage_limiter_states, '_RadauStages._stage_limiter_states', R),
               (T._stages_converged, '_RadauStages._stages_converged', R),
               (T._limiters_at_rest, '_RadauStages._limiters_at_rest', R),
               (T._memo_put, '_CompanionModel._memo_put', CM),
               (T._source_at, '_CompanionModel._source_at', CM),
               (T._stage_source, '_SequentialStages._stage_source',
                'pycircuit.circuit._tran_stages'),
               (T._newton_limiter, '_StepNewton._newton_limiter',
                'pycircuit.circuit._tran_newton'),
               (_tran_radau.remove_row_col, 'remove_row_col', 'pycircuit.circuit.analysis'),
               (_tran_core.passes, 'passes', 'pycircuit.circuit._tran_core'),
               (SubCircuit.limit, 'SubCircuit.limit', 'pycircuit.circuit.circuit'),
               (_hdl_climit.limit_walk, 'limit_walk', 'pycircuit.circuit._hdl_climit'))
        if not all(_paths.genuine(*o) for o in own):
            return None
        _MOD['chain'] = tuple(o for o, _q, _m in own)
        _MOD['mods'] = (_tran_radau, _tran_core, _hdl_climit)
    return _MOD['chain']


def _now(T, mods, cir):
    """What the Python Newton calls now, in `_chain`'s order."""
    tr_, tc, hc = mods
    return (T._coupled_stage_context, T._coupled_stage_system, T._stage_limiter_states,
            T._stages_converged, T._limiters_at_rest, T._memo_put, T._source_at,
            T._stage_source, T._newton_limiter, tr_.remove_row_col, tc.passes,
            getattr(type(cir), 'limit', None), hc.limit_walk)


class _Ctx:
    """One transient's struct and buffers for one core, size and reference
    row (the walk's fields named as `_tran_newton_c._Ctx`'s, whose
    `_walk_ready` sets them)."""

    __slots__ = (
        'buf',
        'core',
        'ffi',
        'iref',
        'keep',
        'n',
        's',
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

    def __init__(self, core, n, iref, ffi, cfn_core, cffi_core, dgesv):
        m = n - 1
        M = 3 * m
        A = MAXA
        self.core, self.n, self.iref, self.ffi = core, n, iref, ffi
        self.s = s = ffi.new('radau_t *')
        self.buf = buf = {
            'Ya': np.empty((A, 3, n)), 'qa': np.empty((A, 3, n)), 'ia': np.empty((A, 3, n)),
            'Ca': np.empty((A, 3, n, n)), 'Ga': np.empty((A, 3, n, n)),
            'qn': np.empty(n), 'u': np.empty((3, n)), 'seed': np.empty((3, n)),
            'K': np.empty((3, n)), 'F': np.empty((3, n)), 'R': np.empty(M),
            'Jf': np.empty(M * M), 'b': np.empty(M), 'blk': np.empty(n * n),
            'stp': np.empty(n), 'tolv': np.empty(2), 'dummy': np.zeros(max(n * n, 1)),
        }
        self.keep = []
        fb = ffi.from_buffer
        for k, arr in buf.items():
            h = fb('double *', arr)
            self.keep.append(h)
            setattr(s, k, h)
        ipiv = np.empty(max(M, 1), dtype=np.int64)
        buf['ipiv'] = ipiv
        self.keep.append(fb('int64_t *', ipiv))
        s.ipiv = self.keep[-1]
        s.dgesv = ffi.cast('void *', dgesv)
        s.core = ffi.cast('void *', int(cffi_core.cast('uintptr_t', core.cs)))
        s.core_fn = ffi.cast('void *', int(cffi_core.cast('uintptr_t', cfn_core)))
        ## the core's own bins: G and i of the pass just run
        for k in ('binG', 'bini'):
            a = core.arrays[k]
            self.keep.append(a)
            setattr(s, k, ffi.cast('double *', a.ctypes.data))
        s.n, s.iref, s.maxa = n, iref, A
        s.walk_fn = ffi.NULL
        s.wn = 0
        self.walk = self.walk_ok = None
        self.walk_els = self.walk_dicts = self.walk_packs = self.walk_kern = ()
        self.walk_hw_els = self.walk_hwk = ()
        self.walk_stamp, self.walk_arms, self.walk_held = -1, 0, 0
        self.walk_watch = ()


def solve(tr, ctx, seed, src, provided_function, lims, nobypass, reltol, abstol, maxit):
    """`_stage_newton(seed)` (undamped, no shunt) in C: the three converged
    stages as fresh full-width arrays, the device memo and the source memo
    as the Python Newton leaves them -- or None: declined (`radau_c:<why>`)
    or handed back (`radau_c:bail:<why>`), nothing left behind, the Python
    Newton runs.  `ctx` is the step's stage context (`h`, `iref`, `m`,
    `qn`, `Amat`, `tstage`), `src` the source at a stage time, `lims` the
    stateful limiters (any declines); the record `tr._radau_c` is made or
    remade where the core, size or reference row changed."""
    if not ENABLED:
        return _no('radau_c:off')
    if provided_function is not None:
        ## (the Python calls it at every assembly's stage times: a caller may
        ## count it)
        return _no('radau_c:pf')
    if lims:
        return _no('radau_c:stateful')
    if not nobypass:
        return _no('radau_c:bypass')
    drv = _driver or driver()
    if not drv:
        return _no('radau_c:driver')
    ffi, cfn, dgesv = drv
    chain = _MOD.get('chain') or _chain()
    if chain is None:
        return _no('radau_c:patched')
    T = chain[0]
    if type(tr) is not T:
        return _no('radau_c:class')
    cir = tr.cir
    if _now(T, _MOD['mods'], cir) != chain[1:] or np.linalg.solve is not _SOLVE:
        return _no('radau_c:patched')
    from pycircuit.circuit import _tran_newton_c as NC
    M = NC._MOD or NC._mods()
    if type(tr.toolkit) is not M['NumericToolkit']:
        return _no('radau_c:toolkit')
    cd, td = cir.__dict__, tr.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd or 'limit' in cd:
        return _no('radau_c:shadow_cir')
    if ('_memo_put' in td or '_source_at' in td or '_coupled_stage_system' in td
            or '_stages_converged' in td or '_stage_limiter_states' in td
            or '_limiters_at_rest' in td or '_newton_limiter' in td):
        return _no('radau_c:shadow_tr')
    if M['limiting'].CIRCUIT_LEVEL:
        return _no('radau_c:stateful')
    memo = td.get('_u_memo')
    if memo is None or not M['companion'].U_MEMO:
        return _no('radau_c:umemo')
    tc = M['core']
    if not (tc.CORE and tc.CORE_PASSES):
        return _no('radau_c:core')
    core = tc.core_for(tr)
    if core is None:
        return _no('radau_c:unservable')
    n, iref = core.n, ctx.iref
    if not (0 <= iref < n) or n < 2 or ctx.m != n - 1:
        return _no('radau_c:iref')
    h = ctx.h
    if type(h) is not float and type(h) is not np.float64:
        return _no('radau_c:h')
    if (type(reltol) is not float and type(reltol) is not int
            and type(reltol) is not np.float64) or type(reltol) is bool:
        return _no('radau_c:tol')
    qn, Amat = ctx.qn, ctx.Amat
    if not (type(qn) is np.ndarray and qn.dtype == np.float64 and qn.shape == (n,)
            and type(Amat) is np.ndarray and Amat.dtype == np.float64 and Amat.shape == (3, 3)):
        return _no('radau_c:ctx')
    Y0 = [np.array(y, dtype=float) for y in seed]
    if len(Y0) != 3 or any(y.shape != (n,) for y in Y0):
        return _no('radau_c:seed')
    has_limiter = tr._newton_limiter() is not None
    if has_limiter:
        Tw = getattr(tr.epar, 'T', 300.0)
        if type(Tw) is not float and type(Tw) is not int and type(Tw) is not np.float64:
            ## (`limit_walk` declines such a temperature: the Python loop limits)
            return _no('radau_c:T')
    rec = td.get('_radau_c')
    if rec is None or rec.core is not core or rec.n != n or rec.iref != iref:
        cffi_core, cfn_core, _dg = tc.driver()
        rec = _Ctx(core, n, iref, ffi, cfn_core, cffi_core, dgesv)
        td['_radau_c'] = rec
    s, buf = rec.s, rec.buf
    if has_limiter:
        if not NC._walk_ready(tr, rec):
            return _no('radau_c:limit')
    elif s.wn:
        s.walk_fn = ffi.NULL
        s.wn = 0
        rec.walk = rec.walk_ok = None
    Tc = core.probe(tr.epar)
    if Tc.__class__ is str:
        return _no('radau_c:ready')
    ## ---- from here the attempt leaves one trace: the step's source memo ----
    analysis = tr.par.analysis
    keys = [(tt, analysis) for tt in ctx.tstage]
    try:
        had = [key in memo for key in keys]
    except TypeError:
        return _no('radau_c:u')
    ## (the counts before the passes, as values -- their dict made only
    ## to roll back)
    u_keys = M['u_keys']
    cvals = tuple(map(_PC.get, u_keys, repeat(0)))
    try:
        us = [src(tt) for tt in ctx.tstage]
    except Exception:                                          # noqa: BLE001
        ## (the call raises again in the Python Newton, where it raised before)
        NC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
        return _no('radau_c:u')
    if not all(type(u) is np.ndarray and u.dtype == np.float64 and u.shape == (n,)
               for u in us):
        NC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
        return _no('radau_c:u')
    for j in range(3):
        buf['u'][j] = us[j]
        buf['seed'][j] = Y0[j]
    np.copyto(buf['qn'], qn)
    for k, v in enumerate(Amat.ravel().tolist()):
        s.A[k] = v
    s.T = Tc
    s.h, s.reltol, s.abstol, s.maxiter = float(h), float(reltol), float(abstol), int(maxit)
    status = cfn(s)
    if status != _CODES['RADAU_CONVERGED']:
        NC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
        if status == _CODES['RADAU_BAIL_WALKSTOP']:
            ## (the walk stopped: its tables are taken again next time)
            rec.walk_ok = None
        return _no('radau_c:bail:' + BAIL.get(status, str(status)))
    ## ---- converged: the memo as the Python Newton leaves it ----
    A = int(s.assemblies)
    Ya, qa, ia, Ca, Ga = buf['Ya'], buf['qa'], buf['ia'], buf['Ca'], buf['Ga']
    put = tr._memo_put
    for a in range(A):
        for j in range(3):
            put(Ya[a, j], {'q': qa[a, j].copy(), 'i': ia[a, j].copy(),
                           'C': Ca[a, j].copy(), 'G': Ga[a, j].copy()})
    _PC['umemo:hit'] += 3 * (A - 1)
    _PC['radau_c:served'] += 1
    _PC['radau_c:assemblies'] += A
    return [Ya[A - 1, j].copy() for j in range(3)]
