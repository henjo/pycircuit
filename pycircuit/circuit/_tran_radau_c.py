"""RADAU'S DENSE STAGE NEWTON IN C: a Radau IIA(3) step's coupled ``3m``
Newton (`_coupled_stage_solver`'s `_stage_newton`, undamped and without a
shunt) in ONE C call, bit for bit the Python Newton (pinned in
`tests/test_pinned_pairs.py`).

OWNS.  On the transient, `tr._radau_c` (`_Ctx`: the `radau_t` struct, its
buffers, the output sets for `MAXA` assemblies (`_out`), the stage KLU's
binding (`_klu_ready`), the last serve's converged block (`final_block`),
and the limiter walk's fields, which `_tran_newton_c._walk_ready` sets;
declared in `transient._RUN_STATE`).  On a serve it writes what the Python
Newton leaves: the device memo (`tr._memo_put`), the step's source memo
(`tr._u_memo`) and its count `umemo:hit`, and on the stage KLU its numeric
and counters.  `_paths`: `radau_c:` -- `radau_c:<why>` a decline,
`radau_c:bail:<why>` a hand-back (`BAIL`), `radau_c:served`,
`radau_c:assemblies`, `radau_c:on_klu` (a serve on the stage KLU),
`radau_c:final` (a converged block a walk's step took).  Switch:
`PYCIRCUIT_RADAU_C=0` (`ENABLED`, read on every call); `STATUS` says why
the path is off where it is.

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

ON A STAGE KLU.  A shooting walk's transient solves the block on the walk's
StructuralKLU (`tr._stage_klu`, `_InnerTransient._new_transient`), the
Python Newton by its `solve`: a refactor on the numeric it holds, a fresh
factor on a zero pivot, the residual's check.  The C makes that sequence
(pcx_kluseq.h, as the DC Newton's) on the block's pattern entries alone:
the core's pattern lies within the KLU's in every block (checked when it is
bound), so the block is zero off it -- what `_ours` takes -- and an
assembly makes the KLU's values and, for their flags, the core's entries on
the reference row and column (the Python makes the full blocks).  The
numeric is StructuralKLU's (pcx_klu.h): a hand-back leaves the entry's, a
serve the one the Python's sequence ends with.  The last assembly, at the
converged stages, is the walk's sensitivities' block at those stages: a
live stage step takes its values (`final_block`) instead of making them.

THE MEMO'S ARRAYS.  A serve records the output slots themselves, views of
a pooled set no earlier record still holds one of (`_out`: the arrays'
reference counts as made), or copies from the transient's own set where
none is free -- the same values, the Python's fresh arrays.

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
evaluate core's two), no build, numpy's error state not ignoring underflow
(the C does not test it: numpy's default ignores it), a toolkit that is
not numeric, a subclass
of the transient or a patched piece of the Newton (`_chain`), an instance
shadow of the circuit's passes or limiter or of the transient's pieces, a
`provided_function`, a stateful limiter or circuit-level limiting, a bypass,
no source memo, a circuit the core cannot serve, a limiting element that
is not its C kernel, a batch not ready, inputs of another type or size; a
stage KLU not StructuralKLU as defined, of another size, not analysed yet
or not holding the core's pattern.

History: `doc/transient_history.md`, `_tran_radau_c.py`.
"""
import os
import sys
from itertools import repeat

import numpy as np

from pycircuit.circuit import _cbuild, _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_RADAU_C', '1') != '0'
STATUS = 'not loaded'
#: the assemblies one call holds (a Newton of more iterations is handed back)
MAXA = 8
#: THE OUTPUT SETS: the slots the C writes its assemblies to (`Ya` .. `Ga`).
#: A served call's device memo takes them as they are -- views, no copies --
#: from a set of the pool (`_Ctx.pool`) no earlier record still holds a view
#: of (its arrays' reference counts as made); at most `OUT_SETS` sets and
#: `OUT_BYTES` of them.  Where none is free the call writes to the
#: transient's own set (`_Ctx.own`, never viewed) and the memo takes copies.
OUT_SETS = 3
OUT_BYTES = 48 * 2 ** 20
_OUTS = ('Ya', 'qa', 'ia', 'Ca', 'Ga')

RADAU_C = _cbuild.csrc('radau.c')
## (radau_t embeds pcx_klu.h's struct: its text first)
RADAU_CDEF = _cbuild.csrc('pcx_klu.h') + _cbuild.csrc('radau.h')

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
    """`(ffi, cfn, dgesv address, dgemv address)` once loaded, or None
    (`STATUS`)."""
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
        from pycircuit.circuit import func
        _MOD['func'] = func
        _driver = (ffi, cfn, nd[2]['dgesv'], nd[2]['dgemv'])
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
        from pycircuit.circuit import linearsolver as LS
        sk = LS.StructuralKLU
        L = 'pycircuit.circuit.linearsolver'
        klu_own = ((sk, 'StructuralKLU', L), (sk.solve, 'StructuralKLU.solve', L),
                   (sk._ours, 'StructuralKLU._ours', L), (sk._factor, 'StructuralKLU._factor', L),
                   (sk._analysed, 'StructuralKLU._analysed', L), (sk._free, 'StructuralKLU._free', L))
        if not all(_paths.genuine(*o) for o in klu_own):
            return None
        _MOD['klu'] = tuple(o for o, _q, _m in klu_own)
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
    `_walk_ready` sets them); `klu` the StructuralKLU the struct is bound to
    (`_klu_ready`), None: numpy's dense solve."""

    __slots__ = (
        'buf',
        'core',
        'ffi',
        'fin',
        'fin_buf',
        'iref',
        'keep',
        'klu',
        'klu_keep',
        'n',
        'own',
        'pool',
        'pre',
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

    def __init__(self, core, n, iref, ffi, cfn_core, cffi_core, dgesv, dgemv):
        m = n - 1
        M = 3 * m
        A = MAXA
        self.core, self.n, self.iref, self.ffi = core, n, iref, ffi
        self.s = s = ffi.new('radau_t *')
        ## (the output slots: `own`, and the pool's, `_out`)
        self.own, self.pool = self._new_out(), []
        self.buf = buf = {
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
        s.dgesv, s.dgemv = ffi.cast('void *', dgesv), ffi.cast('void *', dgemv)
        self.klu, self.klu_keep = None, ()
        ## (the last served call's converged block for the walk's
        ## sensitivities: `final_block`)
        self.fin = self.fin_buf = None
        ## (a converged result held for the next call with its inputs: `hold`)
        self.pre = None
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

    def _new_out(self):
        """An output set: `(its arrays, their pointers, their reference
        counts as made)` -- slots for `MAXA` assemblies."""
        A, n, fb = MAXA, self.n, self.ffi.from_buffer
        arrs = {'Ya': np.empty((A, 3, n)), 'qa': np.empty((A, 3, n)), 'ia': np.empty((A, 3, n)),
                'Ca': np.empty((A, 3, n, n)), 'Ga': np.empty((A, 3, n, n))}
        ptrs = {k: fb('double *', a) for k, a in arrs.items()}
        return arrs, ptrs, {k: sys.getrefcount(arrs[k]) for k in _OUTS}


def _out(rec):
    """The output set a call writes to, and whether the memo may take its
    slots as views: a pooled set no record holds a view of (each array's
    reference count as made), or a new one within `OUT_SETS` and
    `OUT_BYTES`; else the transient's own, copied from."""
    for o in rec.pool:
        arrs, rc = o[0], o[2]
        if all(sys.getrefcount(arrs[k]) == rc[k] for k in _OUTS):
            return o, True
    per = 2 * MAXA * 3 * (rec.n * rec.n + 2 * rec.n) * 8
    if len(rec.pool) < OUT_SETS and (len(rec.pool) + 1) * per <= OUT_BYTES:
        o = rec._new_out()
        rec.pool.append(o)
        return o, True
    return rec.own, False


def _klu_ready(rec, sk):
    """The record bound to the transient's stage-block solve -- numpy's
    dense solve where `sk` (`tr._stage_klu`) is None, else that
    StructuralKLU: its pattern (each entry's block and place in the full
    `C_i`, `G_j`; the core's pattern on the reference row and column, for
    the flags), arrays, analysis, common and the library's entries -- or why
    not: not StructuralKLU as defined, another size, not analysed yet (the
    Python's first solve analyses it), a core pattern it does not hold in
    every block (the block would not be zero off it)."""
    s, ffi = rec.s, rec.ffi
    if sk is None:
        if rec.klu is not None:
            rec.klu, rec.klu_keep = None, ()
            s.klu.klu_sym = ffi.NULL
        return None
    SK = _MOD['klu'][0]
    if type(sk) is not SK:
        return 'klu_kind'
    if (SK.solve, SK._ours, SK._factor, SK._analysed, SK._free) != _MOD['klu'][1:]:
        return 'patched'
    if any(k in sk.__dict__ for k in ('solve', '_ours', '_factor', '_analysed', '_free',
                                      '_other')):
        return 'shadow_klu'
    M = 3 * (rec.n - 1)
    if sk.m != M:
        return 'klu_size'
    if not sk._symbolic:
        return 'klu_unanalysed'
    sym = ffi.cast('void *', sk._symbolic)
    if rec.klu is sk and s.klu.klu_sym == sym:
        return None
    ## (the core's pattern -- every entry its G and C passes write, the
    ## block zero off it -- within the KLU's in each of the nine blocks: the
    ## C then makes the pattern's entries alone)
    n, iref, m = rec.n, rec.iref, rec.n - 1
    pat = np.asarray(rec.core.arrays['pat'], dtype=np.int64)
    pr, pc = np.divmod(pat, n)
    inner = (pr != iref) & (pc != iref)
    rr = pr[inner] - (pr[inner] > iref)
    cc = pc[inner] - (pc[inner] > iref)
    tiled = np.concatenate([(i * m + rr) * M + (j * m + cc) for i in range(3) for j in range(3)])
    take = np.asarray(sk.take, dtype=np.int64)
    if not np.isin(tiled, take).all():
        return 'klu_pattern'
    r, c = np.divmod(take, M)
    bi, rk = np.divmod(r, m)
    bj, ck = np.divmod(c, m)
    edge = pat[~inner]
    ints = (3 * bi + bj, (rk + (rk >= iref)) * n + (ck + (ck >= iref)),
            np.repeat(np.arange(9, dtype=np.int64), len(edge)), np.tile(edge, 9))
    ints = tuple(np.ascontiguousarray(a, dtype=np.int64) for a in ints)
    dbls = (np.empty(len(take)), np.empty(M), np.empty(M), np.empty(2 * M), np.empty(M * M),
            np.empty(max(len(ints[2]), 1)))
    fb = ffi.from_buffer
    ip = tuple(fb('int64_t *', a) for a in ints)
    dp = tuple(fb('double *', a) for a in dbls)
    keep = (sk, fb('int *', sk.Ap), fb('int *', sk.Ai)) + ints + dbls + ip + dp
    s.nact, s.nref = len(take), len(ints[2])
    s.act_b, s.act_f, s.ref_b, s.ref_f = ip
    s.klu.kAx, s.kx, s.ky, s.kt, s.Ad, s.refv = dp
    k = s.klu
    k.kAp, k.kAi = keep[1:3]
    lib = sk._lib
    k.klu_factor, k.klu_refactor, k.klu_solve, k.klu_free = (
        ffi.cast('void *', _cbuild.address(f))
        for f in (lib.klu_factor, lib.klu_refactor, lib.klu_solve, lib.klu_free_numeric))
    k.klu_common = ffi.cast('void *', sk._ct.addressof(sk._common))
    k.klu_sym = sym
    k.resid_tol = float(sk.REFACTOR_RESIDUAL_TOL)
    rec.klu, rec.klu_keep, rec.fin_buf = sk, keep, dbls[0]
    return None


def _ready(tr, ctx, seed, provided_function, lims, nobypass, reltol, step_memo=True):
    """`solve`'s checks before it touches anything, in its order, the record
    made or remade and bound (the limiter walk, the stage KLU): `(the
    record, the temperature, the stage KLU or None, the seeds copied, the
    driver, the Newton glue's modules)` -- or why not (the decline's name,
    `radau_c:<why>`).  `seed` None: no seeds checked, and `step_memo` False:
    no step's source memo asked for (the stage walk's, `_stagewalk_c`,
    which makes its seeds and its sources' table itself)."""
    if not ENABLED:
        return 'off'
    if provided_function is not None:
        ## (the Python calls it at every assembly's stage times: a caller may
        ## count it)
        return 'pf'
    if lims:
        return 'stateful'
    if not nobypass:
        return 'bypass'
    drv = _driver or driver()
    if not drv:
        return 'driver'
    if not _MOD['func']._underflow_ignored():
        return 'underflow'
    ffi, _cfn, dgesv, dgemv = drv
    chain = _MOD.get('chain') or _chain()
    if chain is None:
        return 'patched'
    T = chain[0]
    if type(tr) is not T:
        return 'class'
    cir = tr.cir
    if _now(T, _MOD['mods'], cir) != chain[1:] or np.linalg.solve is not _SOLVE:
        return 'patched'
    from pycircuit.circuit import _tran_newton_c as NC
    M = NC._MOD or NC._mods()
    if type(tr.toolkit) is not M['NumericToolkit']:
        return 'toolkit'
    cd, td = cir.__dict__, tr.__dict__
    if 'G' in cd or 'C' in cd or 'i' in cd or 'q' in cd or 'limit' in cd:
        return 'shadow_cir'
    if ('_memo_put' in td or '_source_at' in td or '_coupled_stage_system' in td
            or '_stages_converged' in td or '_stage_limiter_states' in td
            or '_limiters_at_rest' in td or '_newton_limiter' in td):
        return 'shadow_tr'
    if M['limiting'].CIRCUIT_LEVEL:
        return 'stateful'
    memo = td.get('_u_memo')
    if step_memo and (memo is None or not M['companion'].U_MEMO):
        return 'umemo'
    tc = M['core']
    if not (tc.CORE and tc.CORE_PASSES):
        return 'core'
    core = tc.core_for(tr)
    if core is None:
        return 'unservable'
    n, iref = core.n, ctx.iref
    if not (0 <= iref < n) or n < 2 or ctx.m != n - 1:
        return 'iref'
    h = ctx.h
    if type(h) is not float and type(h) is not np.float64:
        return 'h'
    if (type(reltol) is not float and type(reltol) is not int
            and type(reltol) is not np.float64) or type(reltol) is bool:
        return 'tol'
    qn, Amat = ctx.qn, ctx.Amat
    if not (type(qn) is np.ndarray and qn.dtype == np.float64 and qn.shape == (n,)
            and type(Amat) is np.ndarray and Amat.dtype == np.float64 and Amat.shape == (3, 3)):
        return 'ctx'
    Y0 = None
    if seed is not None:
        Y0 = [np.array(y, dtype=float) for y in seed]
        if len(Y0) != 3 or any(y.shape != (n,) for y in Y0):
            return 'seed'
    has_limiter = tr._newton_limiter() is not None
    if has_limiter:
        Tw = getattr(tr.epar, 'T', 300.0)
        if type(Tw) is not float and type(Tw) is not int and type(Tw) is not np.float64:
            ## (`limit_walk` declines such a temperature: the Python loop limits)
            return 'T'
    rec = td.get('_radau_c')
    if rec is None or rec.core is not core or rec.n != n or rec.iref != iref:
        cffi_core, cfn_core, _dg = tc.driver()
        rec = _Ctx(core, n, iref, ffi, cfn_core, cffi_core, dgesv, dgemv)
        td['_radau_c'] = rec
    s = rec.s
    if has_limiter:
        if not NC._walk_ready(tr, rec):
            return 'limit'
    elif s.wn:
        s.walk_fn = ffi.NULL
        s.wn = 0
        rec.walk = rec.walk_ok = None
    Tc = core.probe(tr.epar)
    if Tc.__class__ is str:
        return 'ready'
    sk = td.get('_stage_klu')
    why = _klu_ready(rec, sk)
    if why is not None:
        return why
    return rec, Tc, sk, Y0, drv, M


def solve(tr, ctx, seed, src, provided_function, lims, nobypass, reltol, abstol, maxit):
    """`_stage_newton(seed)` (undamped, no shunt) in C: the three converged
    stages as fresh full-width arrays, the device memo and the source memo
    as the Python Newton leaves them -- or None: declined (`radau_c:<why>`)
    or handed back (`radau_c:bail:<why>`), nothing left behind, the Python
    Newton runs.  `ctx` is the step's stage context (`h`, `iref`, `m`,
    `qn`, `Amat`, `tstage`), `src` the source at a stage time, `lims` the
    stateful limiters (any declines); the record `tr._radau_c` is made or
    remade where the core, size or reference row changed."""
    got = _ready(tr, ctx, seed, provided_function, lims, nobypass, reltol)
    if got.__class__ is str:
        return _no('radau_c:' + got)
    rec, Tc, sk, Y0, drv, M = got
    from pycircuit.circuit import _tran_newton_c as NC
    ffi, cfn, _dgesv, _dgemv = drv
    s, buf, h, qn, Amat, n = rec.s, rec.buf, ctx.h, ctx.qn, ctx.Amat, rec.n
    memo = tr.__dict__.get('_u_memo')
    rec.fin = None
    if sk is not None:
        num = sk._numeric
        s.klu.klu_num = ffi.cast('void *', num) if num else ffi.NULL
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
    out, viewed = _out(rec)
    for k in _OUTS:
        setattr(s, k, out[1][k])
    for j in range(3):
        buf['u'][j] = us[j]
        buf['seed'][j] = Y0[j]
    np.copyto(buf['qn'], qn)
    for k, v in enumerate(Amat.ravel().tolist()):
        s.A[k] = v
    s.T = Tc
    s.h, s.reltol, s.abstol, s.maxiter = float(h), float(reltol), float(abstol), int(maxit)
    pre, rec.pre = rec.pre, None
    if pre is not None and _held_for(pre, rec, sk, buf, Amat, Tc, h, reltol, abstol, maxit):
        status = _take_held(rec, pre, out[0])
    else:
        status = cfn(s)
    if status != _CODES['RADAU_CONVERGED']:
        NC._undo_source(memo, keys, had, dict(zip(u_keys, cvals)))
        if status == _CODES['RADAU_BAIL_WALKSTOP']:
            ## (the walk stopped: its tables are taken again next time)
            rec.walk_ok = None
        return _no('radau_c:bail:' + BAIL.get(status, str(status)))
    ## ---- converged: the memo as the Python Newton leaves it ----
    A = int(s.assemblies)
    o = out[0]
    Ya, qa, ia, Ca, Ga = o['Ya'], o['qa'], o['ia'], o['Ca'], o['Ga']
    put = tr._memo_put
    if viewed:
        ## (the slots themselves: `_out`'s set, which no call writes to
        ## again while a record holds one)
        for a in range(A):
            for j in range(3):
                put(Ya[a, j], {'q': qa[a, j], 'i': ia[a, j], 'C': Ca[a, j], 'G': Ga[a, j]})
    else:
        for a in range(A):
            for j in range(3):
                put(Ya[a, j], {'q': qa[a, j].copy(), 'i': ia[a, j].copy(),
                               'C': Ca[a, j].copy(), 'G': Ga[a, j].copy()})
    if sk is not None:
        ## (StructuralKLU as its own solve leaves it: the numeric the
        ## sequence ends with, its counters)
        out = int(ffi.cast('uintptr_t', s.klu.klu_out))
        sk._numeric = out or None
        k = s.klu
        sk.refactors += k.refactors
        sk.factors += k.factors
        sk.pivot_fallbacks += k.pivot_fallbacks
        sk.residual_fallbacks += k.residual_fallbacks
        _PC['radau_c:on_klu'] += 1
    _PC['umemo:hit'] += 3 * (A - 1)
    _PC['radau_c:served'] += 1
    _PC['radau_c:assemblies'] += A
    Y = [Ya[A - 1, j].copy() for j in range(3)]
    if sk is not None and s.fin_ok:
        rec.fin = (tuple(Y), float(h), Amat.tobytes(), rec.fin_buf.copy(), sk)
    return Y


def _held_key(rec, sk, qn, u, seed, Amat, Tc, h, reltol, abstol, maxit):
    """What a held result was made from, as `solve` compares it: the
    inputs' bytes, and the stage KLU's numeric and counts as the call finds
    them."""
    return (qn.tobytes(), u.tobytes(), seed.tobytes(), Amat.tobytes(), float(Tc), float(h),
            float(reltol), float(abstol), int(maxit), int(rec.s.wn),
            None if sk is None else (sk._numeric, sk.refactors, sk.factors, sk.pivot_fallbacks,
                                     sk.residual_fallbacks, sk.numerics))


def hold(rec, sk, qn, u, seed, Amat, Tc, h, reltol, abstol, maxit):
    """THE C STAGE NEWTON'S CONVERGED RESULT JUST MADE on the record `rec`
    -- its output slots (`own`'s), the final block, its counts, the stage
    KLU's numeric it ended on -- held for the next `solve`, which serves it
    as made, once, where that call's inputs are these, byte for byte, and
    the stage KLU (its numeric and counts) and the limiter walk are as they
    stand now.  The shooting's stage walk (`_stagewalk_c`) stops after a
    converged Newton it cannot finish the step of: the Python's repeat of
    the step takes the very Newton, its numeric moved once."""
    s = rec.s
    A = int(s.assemblies)
    k = s.klu
    rec.pre = (_held_key(rec, sk, qn, u, seed, Amat, Tc, h, reltol, abstol, maxit),
               (sk, rec.walk, rec.walk_ok), A,
               {nm: rec.own[0][nm][:A].copy() for nm in _OUTS}, int(s.fin_ok),
               None if rec.fin_buf is None else rec.fin_buf.copy(),
               int(rec.ffi.cast('uintptr_t', k.klu_out)),
               (int(k.refactors), int(k.factors), int(k.pivot_fallbacks),
                int(k.residual_fallbacks)), int(s.iters))


def _held_for(pre, rec, sk, buf, Amat, Tc, h, reltol, abstol, maxit):
    """The held result `pre` is this call's (`hold`)."""
    ids = pre[1]
    return (ids[0] is sk and ids[1] is rec.walk and ids[2] is rec.walk_ok
            and pre[0] == _held_key(rec, sk, buf['qn'], buf['u'], buf['seed'], Amat, Tc, h,
                                    reltol, abstol, maxit))


def _take_held(rec, pre, o):
    """The held result `pre` written where the C writes it -- the output set
    `o`, the final block, the struct's counts and numeric: the call's
    result as the C made it."""
    _key, _ids, A, slots, fin_ok, fin, klu_out, counts, iters = pre
    s = rec.s
    for nm in _OUTS:
        o[nm][:A] = slots[nm]
    s.assemblies, s.fin_ok, s.iters = A, fin_ok, iters
    if fin is not None:
        rec.fin_buf[:] = fin
    k = s.klu
    k.klu_out = rec.ffi.cast('void *', klu_out)
    k.refactors, k.factors, k.pivot_fallbacks, k.residual_fallbacks = counts
    _PC['radau_c:held'] += 1
    return _CODES['RADAU_CONVERGED']


def final_block(tr, Y, h, A):
    """The converged stages' block of the step just taken, in the stage
    StructuralKLU's order, where the transient's last call here served it on
    that KLU for `Y` -- the very arrays it returned -- at this `h` and
    tableau `A`: `(the values, the StructuralKLU)`, the block `_stage_block`
    makes from the same stages (the same operations, each entry the same
    bits) and `_ours` would take.  Else None (the walk makes it)."""
    rec = tr.__dict__.get('_radau_c')
    fin = None if rec is None else rec.fin
    if (fin is None or len(Y) != 3 or any(a is not b for a, b in zip(Y, fin[0]))
            or float(h) != fin[1] or np.asarray(A, dtype=float).tobytes() != fin[2]):
        return None
    _PC['radau_c:final'] += 1
    return fin[3], fin[4]
