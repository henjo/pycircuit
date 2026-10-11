"""THE STEPPING LOOP IN C: `_SteppingLoop`'s attempts with a multistep
family (`_LMMSteps`), fixed step or adaptive -- `where` (the grid or the
controller's step, the next breakpoint -- crossed by a grid, landed on by an
adaptive step -- the order drop a landing or a force asks for), `take` (the
stage predictor's seed, the companion, the C Newton, the branch screen),
`judge` (adaptive: the error test, the controller's law, a rejection
retried at the point within its budget) and `accept` (the rings, the
predictor's node, the outputs, the statistics) -- made by `csrc/loop.c`.
It drives the units the Python steps call -- the sources' table
(`_tran_sources`), the predictor (`_tran_pred_c`), the C Newton (its record,
as its last served call left it: `_tran_newton_c`), the branch screen
(`_tran_branch_c`), the error test (`_tran_lte_c`), the controller's law
(`_tran_ctrl_c`) -- and stops before any attempt it cannot make as the
Python would, which the Python loop then makes; the next attempt enters
the C again.

OWNS, on the `_SteppingLoop`: `c`, a `LoopState` (`attach`; None where the
loop never tries) -- `run` (the bound `_Run`: the struct, its buffers, the
units' structs it drives), `tab` (`(the direct source pass, its Table)`),
`kept` (the watch counter at the last kept decline), `skip` and `backoff`
(the backoff)
is `_SteppingLoop.__init__`'s.  `_paths` prefix `loop_c:`:
`loop_c:attempts`, `loop_c:steps`, the declines (`loop_c:<reason>`) and the
bails (`loop_c:bail:<why>`); for the attempts it makes it re-creates the
counts the Python attempts' glue would have made (`newton_c:`, `lte_c:`,
`branch_c:`, `chord:chosen`, `pred:fast`, `src.direct:served`, `u:called`,
`umemo:miss`, `memo.q:same`, `memo.C:equal`).  Switch: `ENABLED` (env
`PYCIRCUIT_LOOP_C=0`), read on every entry.

THE RUN CONTEXT.  Inside the C loop no Python runs between attempts: the
sources are the table, the devices the core's kernels, and a circuit with
an `accept_step` hook, a breakpoint the table does not make, a periodic
state, a provided function, a stateful limiter, an injected piece of the
chain declines.  So what the C Newton's glue checks on every call holds for
every attempt of an entry when its last call SERVED and nothing watched
moved since (`_tran_newton_c`'s `served` mark against `_watch.now()`), and
the same of the error test's glue (`_tran_lte_c`'s); the rest is checked
here, once an entry.  Watching off, the loop does not run.

THE STATE.  Copied in at entry (the rings, the predictor's nodes, the
cached C, the step record, the error test's running reference, the
controller's retry record, the loop's own) and written back after it, as
the Python attempts leave it: `X` and `timelist` a fresh array and float a
step, the statistics, every attribute the attempts write
(`transient._RUN_STATE`: the rings, `_q_cache`, `_C_cache`, `_C_pattern`,
`_iq`, `_Geq`, `_Cmat`, the active integrator and its coefficients, the
predictor's history, the step record, the screen's running maximum), the
controller's (`_ref_running`, `_retry`, `_q_carried`, `last_err`), the
loop's, the run's KLU counts (`numerics`, `refactors`) -- and the `_paths`
tallies the Python attempts make, under their keys.

ENTRY (`run`).  Not entered, and not counted: the switch, the backoff (an
entry that made no attempt doubles the Python attempts made before the
next, up to `BACKOFF_MAX`; one that made any resets it), watching off, a
kept decline while nothing watched has moved, the Newton's (adaptive: the
error test's) last call not served.  DECLINES (`loop_c:<reason>`, kept on
the loop until something watched moves where `_KEEPS` says the reason is
the run's, not the attempt's): no build, a patched or shadowed piece
(`_PIECES`, the table `_own`, `_pieces` and `_SHADOWS` are made from;
`_CIRCUIT_SHADOWS`), another statistics class, loop or family, a fixed-step
loop whose Jacobian something reads, an imposed step, events or a periodic
state, the continuation rescue or a stateful limiter, a circuit hook, an
integrator outside the three, the predictor's settings, no source memo, a
Newton record whose `J` is not formed where the attempts read it, the
sources (no direct pass, no table or a stale one, a breakpoint the table
does not make, the toolkit or numpy's error state), the state's types and
shapes, a device bypass (`bypasstol`), the screen, the predictor (off or
unbuilt), the adaptive judge (`_judge`), the bind.

BAILS (`loop_c:bail:<why>`, `csrc/loop.h`): the C stops before an attempt
a unit hands back or the Python must make (a source time, the predictor's
general path or a flag, the Newton, a screen that is not quiet, a
companion it does not make, the error test, a forced accept); the attempts
made before it are written back.  History: `doc/transient_history.md`.
"""
import collections
import operator
import os
import time

import numpy as np

from pycircuit.circuit import _cbuild, _paths, _watch

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_LOOP_C', '1') != '0'
STATUS = 'not loaded'

LOOP_C = _cbuild.csrc('loop.c')
LOOP_CDEF = _cbuild.csrc('loop.h')
_CODES = _cbuild.constants(LOOP_CDEF)
#: (the law's mark for a kept point that is none of the loop's: ctrl.h)
_POINT_OTHER = _cbuild.constants(_cbuild.csrc('ctrl.h'))['CTRL_POINT_OTHER']
#: (the predictor's history capacity: pred.h)
_PRED_CAP = _cbuild.constants(_cbuild.csrc('pred.h'))
_HIST_MAX, _MAX_DEG = _PRED_CAP['PRED_HIST_MAX'], _PRED_CAP['PRED_FIT_MAX'] - 1
BAIL = _cbuild.bails(_CODES, 'LOOP_BAIL_')
_ABI = _cbuild.constants(_cbuild.csrc('pcx_abi.h'))
#: the outputs' rows an entry holds (full: written back, entered again)
CAP_MAX = 4096
#: the most Python attempts made between two entries after entries that
#: made none (`run`'s backoff)
BACKOFF_MAX = 256
_F64 = np.dtype(np.float64)

_driver = None
_OWN = None


def driver():
    """`(ffi, cfn)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(LOOP_C, LOOP_CDEF, 'pcx_loop')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn)
        STATUS = 'c'
    return _driver or None


#: THE PIECES the attempts call -- the Python the C loop stands in for -- written
#: once: (the object it is read from, as `_own` names it; the attribute; where it
#: is defined: its qualified name and module).  `_own` checks each is its
#: module's own (`_paths.genuine`) once a process; an entry compares them all as
#: they stand (`_pieces`); one read from the transient must not be shadowed on the
#: instance either (`_SHADOWS`).
_PIECES = (
    ('T', '_next_breakpoint', '_StepEvents._next_breakpoint', '_tran_events'),
    ('T', '_attempt_step', 'Transient._attempt_step', 'transient'),
    ('T', '_pred_or', '_StagePredictor._pred_or', '_tran_predictor'),
    ('T', '_predict_state', '_StagePredictor._predict_state', '_tran_predictor'),
    ('T', '_predict_fast', '_StagePredictor._predict_fast', '_tran_predictor'),
    ('T', '_pred_fast_nodes', '_StagePredictor._pred_fast_nodes', '_tran_predictor'),
    ('T', '_fit_fast', '_StagePredictor._fit_fast', '_tran_predictor'),
    ('T', '_pred_note', '_StagePredictor._pred_note', '_tran_predictor'),
    ('T', '_pred_promote', '_StagePredictor._pred_promote', '_tran_predictor'),
    ('T', '_pred_degree', '_StagePredictor._pred_degree', '_tran_predictor'),
    ('T', '_push_history', '_RunHistory._push_history', '_tran_history'),
    ('T', '_roll_history', '_RunHistory._roll_history', '_tran_history'),
    ('history', '_ring_push', '_ring_push', '_tran_history'),
    ('T', '_q_at', '_CompanionModel._q_at', '_tran_companion'),
    ('T', '_branch_after_solve', '_BranchCheck._branch_after_solve', '_tran_branch'),
    ('T', '_branch_screen', '_BranchCheck._branch_screen', '_tran_branch'),
    ('T', '_branch_screen_C', '_BranchCheck._branch_screen_C', '_tran_branch'),
    ('T', '_branch_on', '_BranchCheck._branch_on', '_tran_branch'),
    ('stats', '_note_step', 'TransientStatistics._note_step', 'transient'),
    ('Sub', 'next_event', 'SubCircuit.next_event', 'circuit'),
    ('Sub', 'accept_step', 'SubCircuit.accept_step', 'circuit'),
    ('VS', 'next_event', 'VS.next_event', 'elements'),
    ('IS', 'next_event', 'IS.next_event', 'elements'),
    ('Behavioural', 'next_event', 'Behavioural.next_event', 'hdl'),
    ('T', '_source_at', '_CompanionModel._source_at', '_tran_companion'),
    ('Sub', 'u', 'SubCircuit.u', 'circuit'),
    ('T', '_branch_structural_rank', '_BranchCheck._branch_structural_rank', '_tran_branch'),
    ('family', 'next_breakpoint', '_StepFamily.next_breakpoint', 'transient'),
    ('family', 'after_accept', '_StepFamily.after_accept', 'transient'),
    ('newton_c', 'solve', 'solve', '_tran_newton_c'),
    ('Gear2', 'check_order_drop', 'Gear2Integrator.check_order_drop', 'integrator'),
    ('Euler', 'check_order_drop', 'EulerIntegrator.check_order_drop', 'integrator'),
    ('Trap', 'check_order_drop', 'TrapezoidalIntegrator.check_order_drop', 'integrator'),
    ('Loop', 'execute', '_SteppingLoop.execute', 'transient'),
    ('Loop', 'where', '_SteppingLoop.where', 'transient'),
    ('Loop', 'take', '_SteppingLoop.take', 'transient'),
    ('Loop', 'judge', '_SteppingLoop.judge', 'transient'),
    ('Loop', 'event', '_SteppingLoop.event', 'transient'),
    ('Loop', 'accept', '_SteppingLoop.accept', 'transient'),
    ('family', 'attempt', '_StepFamily.attempt', 'transient'),
    ('LMM', 'judge', '_LMMSteps.judge', 'transient'),
    ('LMM', 'h_after_force', '_LMMSteps.h_after_force', 'transient'),
    ('T', 'solve_timestep', 'Transient.solve_timestep', 'transient'),
    ('T', '_solve_timestep', 'Transient._solve_timestep', 'transient'),
    ('T', '_lte_inputs', '_RunHistory._lte_inputs', '_tran_history'),
    ('T', '_excursion_ratio', '_StepEvents._excursion_ratio', '_tran_events'),
    ('T', '_lte_mask', '_RunHistory._lte_mask', '_tran_history'),
    ('T', '_sparse_linsolver', '_RunHistory._sparse_linsolver', '_tran_history'),
    ('IC', '_evaluate', 'IntegralController._evaluate', 'stepcontroller'),
    ('SC', 'evaluate_step', 'StepController.evaluate_step', 'stepcontroller'),
    ('SC', '_max_error', 'StepController._max_error', 'stepcontroller'),
    ('SC', '_band_decision', 'StepController._band_decision', 'stepcontroller'),
    ('SC', '_band_target', 'StepController._band_target', 'stepcontroller'),
    ('SC', '_retry_exponent', 'StepController._retry_exponent', 'stepcontroller'),
    ('SC', '_judged', 'StepController._judged', 'stepcontroller'),
    ('SC', '_damp', 'StepController._damp', 'stepcontroller'),
    ('lte_c', 'max_error', 'max_error', '_tran_lte_c'),
)

#: (an instance attribute of one of these names is a piece the Python attempts
#: would call and the C would not)
_SHADOWS = tuple(attr for owner, attr, _q, _m in _PIECES if owner == 'T')


def _pieces(O):
    """The pieces as they stand now, a tuple per object they are read from
    (`O['getters']`, in `_PIECES`'s order): what an entry compares with
    `O['chain']`."""
    return tuple(get(obj) for obj, get in O['getters'])


def _own():
    """The pieces the loop stands in for, as defined -- read once every one
    is its module's own (`_paths.genuine`), else None -- with what an entry
    compares against them (`_pieces`)."""
    global _OWN
    if _OWN is None:
        from pycircuit.circuit import (
            _tran_companion,
            _tran_history,
            _tran_lte_c,
            _tran_newton_c,
            _tran_predictor,
            hdl,
        )
        from pycircuit.circuit import elements as E
        from pycircuit.circuit import integrator as I
        from pycircuit.circuit import stepcontroller as sc
        from pycircuit.circuit.circuit import Circuit, SubCircuit
        from pycircuit.circuit.transient import (
            Transient,
            TransientStatistics,
            _LMMSteps,
            _StepFamily,
            _SteppingLoop,
        )
        M = 'pycircuit.circuit.'
        T, F, L, LMM, SC = Transient, _StepFamily, _SteppingLoop, _LMMSteps, sc.StepController
        IC = sc.IntegralController
        owners = {'T': T, 'history': _tran_history, 'stats': TransientStatistics,
                  'Sub': SubCircuit, 'VS': E.VS, 'IS': E.IS, 'Behavioural': hdl.Behavioural,
                  'family': F, 'newton_c': _tran_newton_c, 'Gear2': I.Gear2Integrator,
                  'Euler': I.EulerIntegrator, 'Trap': I.TrapezoidalIntegrator, 'Loop': L,
                  'LMM': LMM, 'IC': IC, 'SC': SC, 'lte_c': _tran_lte_c}
        if not all(_paths.genuine(getattr(owners[o], a), q, M + m) for o, a, q, m in _PIECES):
            return None
        names = {}
        for o, a, _q, _m in _PIECES:
            names.setdefault(o, []).append(a)
        O = {
            'T': T, 'stats': TransientStatistics, 'Circuit': Circuit, 'Sub': SubCircuit,
            'history': _tran_history, 'E': E, 'hdl': hdl, 'I': I,
            'source_next': (E.VS.next_event, E.IS.next_event),
            'inf_next': hdl.Behavioural.next_event,
            'integrators': {I.Gear2Integrator: _ABI['PCX_FORMULA_GEAR2'],
                            I.EulerIntegrator: _ABI['PCX_FORMULA_EULER'],
                            I.TrapezoidalIntegrator: _ABI['PCX_FORMULA_TRAPEZOIDAL']},
            'Euler': I.EulerIntegrator,
            'zs_ratio': I.ZERO_STABILITY_RATIO,
            'predictor': _tran_predictor, 'companion': _tran_companion,
            'newton_c': _tran_newton_c, 'lte_c': _tran_lte_c, 'family': F, 'Loop': L,
            'LMM': LMM, 'SC': SC, 'IC': IC,
        }
        O['getters'] = tuple((owners[o], operator.attrgetter(*a)) for o, a in names.items())
        O['chain'] = _pieces(O)
        _OWN = O
    return _OWN


class _Run:
    """One stepping loop's C loop: the struct, its buffers, the units'
    structs it drives -- bound to the C Newton's record, the source table and
    the screen it was made for (and, adaptive, the error test's context)."""

    __slots__ = ('br', 'buf', 'ct', 'ffi', 'keep', 'lte', 'npat', 'nrec', 'pr', 's', 'sk', 'tab')


def _addr(ffi, obj):
    return int(ffi.cast('uintptr_t', obj))


def _bind(tr, rec, tab, n, hist_len, screen, judge):
    """The loop's `_Run` for the C Newton's record `rec`, the source table
    `tab`, the branch screen `screen` (`(its ffi, its entry, the core's
    pattern or None)`, or None: nothing screened) and, adaptive, `judge`
    (`_judge`'s `(the controller, the run's settings, the sparse solver or
    None, the error test's mask or None)`) -- or None: the sparse solver
    could not be bound."""
    ffi, _cfn = driver()
    run = _Run()
    run.ffi, run.nrec, run.tab = ffi, rec, tab
    s = run.s = ffi.new('loop_t *')
    _bind_buffers(run, n, hist_len)
    _bind_steps(run, tr, rec, tab, n)
    _bind_screen(run, tr, n, screen)
    run.lte = run.ct = run.sk = None
    s.lte = s.lte_fn = s.ctrl = s.ctrl_fn = ffi.NULL
    if judge is not None and not _bind_judge(run, tr, n, judge):
        return None
    s.n, s.hist_len, s.cap = n, hist_len, CAP_MAX
    return run


def _plug(s, ffi, unit, uffi, state, entry):
    """The loop pointed at a unit (`loop.h`'s `<unit>`, `<unit>_fn`): its
    struct and its entry, made by the unit's own `uffi`."""
    setattr(s, unit, ffi.cast('void *', _addr(uffi, state)))
    setattr(s, unit + '_fn', ffi.cast('void *', _addr(uffi, entry)))


def _bind_buffers(run, n, hist_len):
    """The loop's own buffers -- the rings, the predictor's nodes, the two
    cached C and Geq, the attempt's point, the outputs -- each pointed at by
    its field and kept (`run.keep`)."""
    ffi, s = run.ffi, run.s
    cap = CAP_MAX
    run.buf = buf = {
        'qring': np.empty(hist_len * n), 'iqring': np.empty(hist_len * n),
        'xlast': np.empty(n), 'ht': np.empty(_HIST_MAX), 'hx': np.empty(_HIST_MAX * n),
        ## (zero off the pattern: where the core writes only the pattern's
        ## entries, the rest stay as the Newton's own buffers keep them)
        'cC0': np.zeros(n * n), 'cC1': np.zeros(n * n),
        'Geq0': np.zeros(n * n), 'Geq1': np.zeros(n * n),
        'run': np.empty(n), 'xatt': np.empty(n), 'qatt': np.empty(n), 'iqatt': np.empty(n),
        'X': np.empty(cap * n), 'T': np.empty(cap),
    }
    keep = run.keep = []
    fb = ffi.from_buffer
    for k in ('qring', 'iqring', 'xlast', 'ht', 'hx', 'run', 'xatt', 'qatt', 'iqatt', 'X', 'T'):
        keep.append(fb('double *', buf[k]))
        setattr(s, k, keep[-1])
    for i in (0, 1):
        keep.append(fb('double *', buf[f'cC{i}']))
        s.cC[i] = keep[-1]
        keep.append(fb('double *', buf[f'Geq{i}']))
        s.Geq[i] = keep[-1]


def _bind_steps(run, tr, rec, tab, n):
    """The units every attempt calls: the C Newton (its record `rec`), the
    predictor (a state of its own, its nodes the loop's `ht` and `hx`) and
    the sources' table `tab`."""
    from pycircuit.circuit import _tran_newton_c, _tran_pred_c
    ffi, s = run.ffi, run.s
    nffi, ncfn, _lapack = _tran_newton_c.driver()
    _plug(s, ffi, 'newton', nffi, rec.s, ncfn)
    pffi, pcfn, dgesv, dgemv = _tran_pred_c.driver()
    ps, pbuf, pkeep = _tran_pred_c.new_state(pffi, dgesv, dgemv, n)
    ps.ht = pffi.cast('double *', _addr(ffi, s.ht))
    ps.hx = pffi.cast('double *', _addr(ffi, s.hx))
    ps.deg = int(tr._pred_degree())
    ps.clamp = float(tr.PRED_CLAMP)
    run.pr = (ps, pbuf, pkeep)
    _plug(s, ffi, 'pred', pffi, ps, pcfn)
    _plug(s, ffi, 'src', tab.ffi, tab.s, tab.cfn)


def _bind_screen(run, tr, n, screen):
    """The branch screen (`screen`: `(its ffi, its entry, the core's pattern
    or None)`), or none: `branch` NULL."""
    ffi, s = run.ffi, run.s
    run.br, run.npat = None, 0
    if screen is None:
        s.branch = s.branch_fn = ffi.NULL
        return
    bffi, bcfn, pat = screen
    bs = bffi.new('branch_t *')
    bs.n, bs.npat = n, 0
    bs.tol_rel = float(tr.BRANCH_SCREEN_TOL)
    bkeep = [bs]
    if pat is not None:
        bkeep += [pat, bffi.from_buffer('int64_t *', pat)]
        bs.pat, bs.npat = bkeep[-1], len(pat)
        run.npat = len(pat)
    run.br = bkeep
    _plug(s, ffi, 'branch', bffi, bs, bcfn)


def _bind_judge(run, tr, n, judge):
    """The adaptive loop's judge (`judge`, `_judge`'s): the error test -- a
    context of its own, set as the controller's `_tran_lte_c` sets one, KLU
    bound where the run solves sparse -- and the controller's law.  False:
    the sparse solver could not be bound."""
    from pycircuit.circuit import _tran_ctrl_c, _tran_lte_c
    ffi, s = run.ffi, run.s
    ctrl, rset, sk, mask = judge
    lffi, lcfn, ldgesv = _tran_lte_c.driver()
    ctx = _tran_lte_c._Ctx(lffi, ldgesv, n, tr.irefnode)
    if sk is not None and not _tran_lte_c._klu_bind(ctx, sk):
        return False
    le = ctx.s
    ab = rset.abstol
    if type(ab) is np.ndarray:
        ctx.buf['abst'][:] = ab
        le.abs_vec, le.abstol = 1, 0.0
    else:
        le.abs_vec, le.abstol = 0, float(ab)
    le.reltol, le.trtol = float(tr.par.reltol), float(tr.LTERATIO)
    le.relref = _tran_lte_c._RELREF.get(ctrl.relref, _tran_lte_c._CODES['LTE_RELREF_SIGGLOBAL'])
    le.nn = len(tr.cir.nodes)
    ## (the unknowns the test judges, `lte_algebraic`: the run's own)
    le.has_mask = int(mask is not None)
    if mask is not None:
        ctx.buf['mask'][:] = mask
    if sk is not None:
        le.klu_num = lffi.cast('void *', sk._numeric)
        le.resid_tol = float(sk.REFACTOR_RESIDUAL_TOL)
    run.lte, run.sk = ctx, sk
    _plug(s, ffi, 'lte', lffi, le, lcfn)
    cffi_, ccfn = _tran_ctrl_c.driver()
    ct = _tran_ctrl_c.new_state(cffi_, ctrl.retry_order)
    run.ct = ct
    _plug(s, ffi, 'ctrl', cffi_, ct, ccfn)
    s.keeps_run = int(le.relref != _tran_lte_c._CODES['LTE_RELREF_POINTLOCAL'])
    return True


#: THE DECLINES (`loop_c:<reason>`), each the run's -- kept on the loop
#: (`c.kept`) until something watched moves -- or the attempt's alone (the
#: next Python attempt enters again)
_KEEPS = {
    ## the run's: its pieces, classes, settings, circuit and sources
    'driver': True, 'patched': True, 'stats': True, 'loop': True, 'shadow': True,
    'events': True, 'state': True, 'circuit': True, 'integrator': True, 'predictor': True,
    'umemo': True, 'sources': True, 'table': True, 'breakpoints': True, 'minbreak': True,
    'screen': True, 'pred': True, 'bind': True,
    ## (the adaptive judge's, `_judge`)
    'controller': True, 'veto': True, 'mask': True, 'abstol': True, 'step_bounds': True,
    'tolerance': True, 'budget': True,
    ## the attempt's: what the last Python attempt left, or a state the next
    ## one may change
    'imposed': False, 'newton': False, 'table_stale': False, 'toolkit': False,
    'rings': False, 'history': False, 'cache': False, 'reference': False, 'sparse': False,
}


class LoopState:
    """What the C loop keeps on one `_SteppingLoop` (`loop.c`): the bound run
    (`run`, `_bound`), the source table it reads (`tab`: `(the direct pass,
    its Table)`, `_sources`), the watch counter at the last kept decline
    (`kept`, `_declined`), and the backoff (`skip`: the Python attempts left
    before the next entry; `backoff`: the spacing that set it)."""

    __slots__ = ('backoff', 'kept', 'run', 'skip', 'tab')

    def __init__(self):
        self.run = self.tab = self.kept = None
        self.skip = self.backoff = 0


def attach(loop):
    """The C loop's state for a `_SteppingLoop` that may enter it -- a fixed
    run whose attempts read no Jacobian (`loop.j_unread`), or an adaptive
    multistep one -- else None, and the loop never tries.  `run` checks the
    rest at every entry."""
    from pycircuit.circuit.transient import _LMMSteps
    if loop.j_unread or (not loop.fixed and type(loop.family) is _LMMSteps):
        return LoopState()
    return None


def _declined(loop, why):
    """The entry declined, counted -- and kept on the loop until something
    watched moves where the reason is the run's (`_KEEPS`)."""
    if _KEEPS[why]:
        loop.c.kept = _watch.now()
    _no('loop_c:' + why)
    return False


def _screen(tr, td, rec, n):
    """The branch screen the attempts make: `(None,)` -- nothing screened
    (off, or `C` structurally empty) -- `(ffi, entry, pattern or None)` in
    C, or None: not in C."""
    if not tr._branch_on():
        return (None,)
    r0, scale0 = tr._branch_structural_rank()
    if not (r0 > 0 and scale0 > 0.0):
        return (None,)
    from pycircuit.circuit import _tran_branch_c
    bd = _tran_branch_c.driver()
    if bd is None or not _tran_branch_c.ENABLED or type(tr.BRANCH_SCREEN_TOL) is not float \
            or type(td.get('_branch_cmax', 0.0)) is not float:
        return None
    core = rec.core
    pat = core.arrays['pat'] if core.npat else None
    if pat is None and n * n > _cbuild.INPLACE_MIN_ENTRIES:
        return None
    return (bd[0], bd[1], pat)


def _judge(loop, tr, O, n):
    """The adaptive attempts' judge, checked: `(the controller, the run's
    settings, the sparse solver or None, the error test's mask or None)`,
    or a decline reason (str)."""
    from pycircuit.circuit import _tran_ctrl_c, _tran_lte_c
    from pycircuit.circuit import stepcontroller as sc
    ctrl = tr.step_controller
    cd = ctrl.__dict__
    if type(ctrl) is not O['IC'] or any(k in cd for k in _CONTROLLER_SHADOWS):
        return 'controller'
    if not (ctrl.lte_gamma_min == 0.0 and ctrl.lte_gamma_max == 1.0 and ctrl.lte_eta is None
            and ctrl.retry_order in _tran_ctrl_c.MODES and type(ctrl.relref) is str):
        return 'controller'
    if (sc.MIN_SHRINK_RATIO, sc.MAX_GROWTH_RATIO, sc.RETRY_ORDER_MIN) != (0.2, 2.0, 1.0):
        return 'controller'
    if _tran_ctrl_c.driver() is None or not _tran_ctrl_c.ENABLED or not _tran_lte_c.ENABLED:
        return 'controller'
    if tr.par.max_dv_step is not None or tr.par.max_di_step is not None:
        return 'veto'
    mask = tr._lte_mask()
    if mask is not None and not (type(mask) is np.ndarray and mask.shape == (n,)):
        return 'mask'
    rset = loop.run
    ab = rset.abstol
    if not ((type(ab) is np.ndarray and ab.dtype is _F64 and ab.shape == (n,))
            or type(ab) is float):
        return 'abstol'
    if type(rset.max_step) is not float or type(loop.minstep) is not float:
        return 'step_bounds'
    if type(tr.par.reltol) not in (float, int) or type(tr.LTERATIO) is not float:
        return 'tolerance'
    fam = loop.family
    if type(fam.max_retries) is not int or type(fam.max_reject) is not int:
        return 'budget'
    run = cd.get('_ref_running')
    if run is not None and not (type(run) is np.ndarray and run.dtype is _F64
                                and run.shape == (n,)):
        return 'reference'
    q = cd.get('_q_carried')
    if q is not None and type(q) is not float:
        return 'controller'
    retry = cd.get('_retry')
    if retry is not None and not (type(retry) is tuple and len(retry) == 3
                                  and type(retry[1]) is float and type(retry[2]) is float):
        return 'controller'
    sk = tr._sparse_linsolver()
    if sk is not None:
        sk = sk._choice
        if not (sk._numeric and sk._symbolic):
            return 'sparse'
    return (ctrl, rset, sk, mask)


#: the controller's methods the judge calls, shadowed on the instance
_CONTROLLER_SHADOWS = ('evaluate_step', '_evaluate', '_max_error', '_band_decision',
                       '_band_target', '_retry_exponent', '_judged', '_damp', '_charge_lte',
                       '_normalised', '_reference', 'tolerance')


def _why_not_the_run(loop, tr, td, rec, O, adaptive):
    """Why the run itself is not one the C loop makes (a reason, `_KEEPS`),
    or None: the pieces patched or shadowed, another statistics class, loop
    or family, an imposed step, events or a periodic state, the rescue or a
    stateful limiter, a circuit hook, an integrator outside the three, the
    predictor's settings, no source memo, a Newton record whose J is not
    formed where the attempts read it."""
    if O is None or type(tr) is not O['T'] or _pieces(O) != O['chain']:
        return 'patched'
    if type(tr.statistics) is not O['stats']:
        return 'stats'
    if type(loop) is not O['Loop'] or type(loop.family) is not O['LMM'] \
            or loop.provided_function is not None or (not adaptive and not loop.j_unread):
        return 'loop'
    if loop.imposed:
        return 'imposed'
    for k in _SHADOWS:
        if k in td:
            return 'shadow'
    if getattr(tr, '_ev_rows', None) is not None or getattr(tr, '_periodic_rows', None):
        return 'events'
    if getattr(tr, '_continuation_rescue', False) or td.get('_stateful_lims'):
        return 'state'
    cir = tr.cir
    cd = cir.__dict__
    if (type(cir).next_event is not O['Sub'].next_event
            or type(cir).accept_step is not O['Sub'].accept_step
            or type(cir).u is not O['Sub'].u or any(k in cd for k in _CIRCUIT_SHADOWS)
            or cir._hook_elements('accept_step')):
        return 'circuit'
    base = tr.base_integrator
    code = O['integrators'].get(type(base))
    if code is None or 'check_order_drop' in base.__dict__ or (
            adaptive and 'compute_lte' in base.__dict__):
        return 'integrator'
    P = O['predictor']
    if (not P.PRED_FAST or tr.stage_predictor != 'on' or tr.PRED_HIST_MAX != _HIST_MAX
            or type(tr)._fit is not P._FIT or type(tr.PRED_CLAMP) is not float
            or not 1 <= tr._pred_degree() <= _MAX_DEG):
        return 'predictor'
    if not O['companion'].U_MEMO:
        return 'umemo'
    ## (the Newton's record as its last call left it: J formed where the
    ## attempts read it)
    NC = O['newton_c']
    if rec.s.fold_j != (NC.FOLD_J if adaptive else NC.FOLD_C):
        return 'newton'
    return None


def _sources(loop, tr, O):
    """The sources as the C reads them -- the direct pass's table (kept on
    the loop, `c.tab`, while it holds), every other element's breakpoint
    inf, the toolkit's and numpy's state -- and `minbreak`: `(table,
    minbreak)`, or why not (a reason, `_KEEPS`)."""
    cir = tr.cir
    cd = cir.__dict__
    ## the sources: the direct pass's table; every other element's
    ## breakpoint inf (`Behavioural.next_event`)
    from pycircuit.circuit import _tran_sources, func
    plan = cd.get('_stamp_plan')
    sp = plan.methods.get('src:u') if plan is not None else None
    sd = sp.direct if sp is not None else None
    ep = _watch.EPOCH
    if sd is None or ep is None or sd.stamp != ep.value or not sd.same(cir, O['Circuit']):
        return 'sources'
    ct_ = loop.c.tab
    if ct_ is not None and ct_[0] is sd:
        tab = ct_[1]
    else:
        tab = _tran_sources.table_for(sd)
        if tab is None:
            return 'table'
        loop.c.tab = (sd, tab)
    if not tab.holds():
        loop.c.tab = None
        return 'table_stale'
    tab_els = {id(e[0]) for e in tab.entries}
    src_next, inf_next = O['source_next'], O['inf_next']
    for _inst, el in cir._hook_elements('next_event'):
        f = type(el).next_event
        if 'next_event' in el.__dict__ or not (
                (f in src_next and id(el) in tab_els) or f is inf_next):
            return 'breakpoints'
    tk = tab.tk
    if not (tk.exp is np.exp and tk.sin is np.sin and tk.where is np.where) or (
            tab.uses_libm and ((tab.has_sin and not func.SIN_SCALAR)
                               or not func._underflow_ignored())):
        return 'toolkit'
    minbreak = tr.par.minbreak
    if type(minbreak) is not float:
        return 'minbreak'
    return tab, minbreak


#: (what an entry read of the Python loop's state: `_read`)
_Read = collections.namedtuple('_Read', 'n ql iql hist_len hist X xl cx cC dl dl2')


def _read(loop, tr, td, rec, adaptive):
    """The state the attempts read, as the Python loop leaves it (`_Read`):
    the rings, the predictor's history, the last point and the cached C,
    the last steps -- or why not (a reason, `_KEEPS`)."""
    ## the state the attempts read
    n = rec.n
    ql, iql = tr._qlast, tr._iqlast
    if not (type(ql) is np.ndarray and ql.dtype is _F64 and ql.ndim == 2 and ql.shape[1] == n
            and ql.shape[0] >= (3 if adaptive else 2) and type(iql) is np.ndarray
            and iql.dtype is _F64 and iql.shape == ql.shape):
        return 'rings'
    hist_len = ql.shape[0]
    hist = td.get('_pred_hist') or []
    if (type(hist) is not list or len(hist) > _HIST_MAX or td.get('_pred_pending') is not None
            or not all(type(e) is tuple and len(e) == 2 and type(e[0]) is float
                       and type(e[1]) is np.ndarray and e[1].dtype is _F64 and e[1].shape == (n,)
                       for e in hist)):
        return 'history'
    X = loop.X
    xl = X[-1] if X else None
    cached = td.get('_C_cache')
    if cached is None or float(getattr(tr.epar, 'bypasstol', -1.0) or -1.0) >= 0.0:
        return 'cache'
    cx, cC = cached
    if not (type(xl) is np.ndarray and xl.dtype is _F64 and xl.shape == (n,)
            and type(cx) is np.ndarray and cx.dtype is _F64 and cx.shape == (n,)
            and type(cC) is np.ndarray and cC.dtype is _F64 and cC.shape == (n, n)):
        return 'cache'
    dl, dl2 = tr._dt_last, tr._dt_last2
    stats = tr.statistics
    if (type(dl) is not float or not (dl2 is None or type(dl2) is float) or td.get('_no_history')
            or not (stats.min_step is None or type(stats.min_step) is float)
            or not (stats.max_step is None or type(stats.max_step) is float)
            or type(loop.t) is not float or type(loop.h) is not float):
        ## (declined as `history`, not `steps`: `loop_c:steps` counts the
        ## steps made in C)
        return 'history'
    return _Read(n, ql, iql, hist_len, hist, X, xl, cx, cC, dl, dl2)


def _bound(loop, tr, rec, tab, rd, screen, judge):
    """The loop's bound run (`c.run`), bound again where the Newton's record,
    the table, the screen's presence, the sizes or the judge's solver or
    mask moved since; None where it cannot be bound."""
    n, hist_len = rd.n, rd.hist_len
    ## bound once a loop (again where the Newton's record, the table, the
    ## screen's presence, the sizes or the judge's solver or mask moved)
    run_ = loop.c.run
    if (run_ is None or run_.nrec is not rec or run_.tab is not tab
            or (run_.br is None) != (screen is None) or len(run_.buf['xlast']) != n
            or len(run_.buf['qring']) != hist_len * n
            or (run_.lte is None) != (judge is None)
            or (judge is not None and (run_.sk is not judge[2] or run_.lte.s.has_mask
                                       != (judge[3] is not None)))):
        run_ = _bind(tr, rec, tab, n, hist_len, screen, judge)
        if run_ is not None:
            loop.c.run = run_
    return run_


def _load(loop, tr, td, run_, rd, minbreak, adaptive, judge, O):
    """The struct loaded with the state as the Python loop leaves it -- the
    rings, the history, the cached C, the step record, the run's constants
    and, adaptive, the controller's record.  Returns the controller's retry
    record as it was (`_write_back` reads it)."""
    n, ql, iql, _hist_len, hist, X, xl, cx, cC, dl, dl2 = rd
    base = tr.base_integrator
    code = O['integrators'][type(base)]
    stats = tr.statistics
    s, buf = run_.s, run_.buf
    ## in: the state as the Python loop leaves it
    buf['qring'][:] = ql.reshape(-1)
    buf['iqring'][:] = iql.reshape(-1)
    buf['xlast'][:] = xl
    buf['xatt'][:] = cx
    s.hk = len(hist)
    for i, (tt, xx) in enumerate(hist):
        buf['ht'][i] = tt
        buf['hx'][i * n:(i + 1) * n] = xx
    buf['cC0'][:] = cC.reshape(-1)
    s.cur = 0
    s.cmax = td.get('_branch_cmax', 0.0)
    s.t, s.h = loop.t, loop.h
    s.landing, s.is_first = int(bool(loop.landing)), int(bool(tr._is_first_step))
    s.order_drop = int(bool(loop.order_drop))
    s.rejects, s.point_retries = loop.rejects, loop.point_retries
    s.npoint = len(X)
    s.dt_last, s.have_dt_last = dl, 1
    s.dt_last2, s.have_dt_last2 = (dl2, 1) if dl2 is not None else (0.0, 0)
    s.tend, s.timestep, s.minbreak = float(loop.tend), float(loop.timestep), minbreak
    s.adaptive = int(adaptive)
    s.max_step = float(loop.run.max_step)
    s.minstep = float(loop.minstep)
    s.max_retries, s.max_reject = loop.family.max_retries, loop.family.max_reject
    s.base = code
    s.zs_ratio = float(O['zs_ratio'])
    s.shrink_guard = int(bool(getattr(base, 'shrink_guard', True)))
    s.have_hm = int(stats.min_step is not None)
    s.hmin = stats.min_step if stats.min_step is not None else 0.0
    s.hmax = stats.max_step if stats.max_step is not None else 0.0
    retry0 = None
    if adaptive:
        ctrl = judge[0]
        ref = ctrl.__dict__.get('_ref_running')
        s.have_run = int(ref is not None)
        if ref is not None:
            buf['run'][:] = ref
        ct = run_.ct
        retry0 = ctrl.__dict__.get('_retry')
        ct.have_retry = int(retry0 is not None)
        if retry0 is not None:
            ## (the point the controller kept, by `is`: the last accepted one
            ## or another -- never one of the C's)
            ct.retry_point = len(X) - 1 if retry0[0] is xl else _POINT_OTHER
            ct.retry_h, ct.retry_err = retry0[1], retry0[2]
        q = ctrl.__dict__.get('_q_carried')
        ct.have_q = int(q is not None)
        ct.q_carried = q if q is not None else 0.0
        if run_.sk is not None:
            run_.lte.s.klu_num = run_.lte.ffi.cast('void *', run_.sk._numeric)
    return retry0


def _call(drv, rec, s, stats):
    """The C's attempts: its status.  The Newton record's own pointers are
    saved before and restored after, whatever the status (the C points
    them at the loop's buffers while it runs); the time is the solve's."""
    ## (the Newton's own pointers: the C points them elsewhere while it runs)
    nt = rec.s
    saved = (nt.q1, nt.q2, nt.iq1, nt.u, nt.xseed, nt.cx, nt.cC, nt.Cj, nt.Geqj, nt.have_cache)
    t0 = time.perf_counter()
    try:
        status = drv[1](s)
    finally:
        (nt.q1, nt.q2, nt.iq1, nt.u, nt.xseed, nt.cx, nt.cC, nt.Cj, nt.Geqj,
         nt.have_cache) = saved
        stats.solve_seconds += time.perf_counter() - t0
    return status


def _eligible(loop, tr, td, rec, adaptive):
    """What an entry that passed the gate runs with -- `(the driver, the
    pieces (`_own`), the source table, minbreak, the state (`_read`), the
    screen or None, the judge or None)` -- or why not (a reason, `_KEEPS`):
    the run (`_why_not_the_run`), the sources (`_sources`), the state
    (`_read`), the screen, the predictor, the judge (`_judge`), in turn."""
    drv = _driver or driver()
    if not drv:
        return 'driver'
    O = _OWN or _own()
    why = _why_not_the_run(loop, tr, td, rec, O, adaptive)
    if why is not None:
        return why
    got = _sources(loop, tr, O)
    if type(got) is str:
        return got
    tab, minbreak = got
    rd = _read(loop, tr, td, rec, adaptive)
    if type(rd) is str:
        return rd
    screen = _screen(tr, td, rec, rd.n)
    if screen is None:
        return 'screen'
    if screen[0] is None:
        screen = None
    ## (the predictor: off by its switch, or unbuilt -- as the screen's and
    ## the law's, read on every entry)
    from pycircuit.circuit import _tran_pred_c
    if not _tran_pred_c.ENABLED or _tran_pred_c.driver() is None:
        return 'pred'
    judge = None
    if adaptive:
        judge = _judge(loop, tr, O, rd.n)
        if type(judge) is str:
            return judge
    return drv, O, tab, minbreak, rd, screen, judge


def run(loop):
    """Make as many of `loop`'s attempts in C as it can, from where the loop
    stands: True where it made any (the loop's, the transient's and the
    controller's state written back as the Python attempts leave it,
    `_write_back`), else False -- nothing of the run's written; the Python
    loop makes the next attempt, and the one after enters again.  The
    entry's gate first (the cheap returns, uncounted, made once a Python
    attempt: kept here, no call), then each job in turn: what the run
    reads, each check that declines (`_eligible`), the bind (`_bound`), the
    load (`_load`), the call (`_call`), the write-back.
    What keeps it out, what it declines and where the C bails: the module's
    ENTRY."""
    if not ENABLED:
        return False
    ## (entries that made no attempt -- a screen that fires every step, the
    ## SVD the Python's -- are spaced out: each such entry doubles the
    ## Python attempts made before the next one, up to `BACKOFF_MAX`, whatever
    ## else holds the loop back meanwhile; an entry that makes any starts
    ## again from none.  Checked first: made once a Python attempt while the
    ## entries are spaced, it is the cheapest return.)
    if loop.c.skip:
        loop.c.skip -= 1
        return False
    now = _watch.now()
    ## (watching off: nothing would tell the attempts' checks have moved --
    ## the counter stays -1 -- so the Python attempts make them, every one)
    if now == -1 or loop.c.kept == now:
        return False
    tr = loop.tr
    td = tr.__dict__
    rec = td.get('_newton_c')
    if rec is None or rec.served != now:
        ## (the C Newton's last call did not serve, or something watched
        ## moved since: the Python attempt, which serves or not by its checks)
        return False
    adaptive = not loop.fixed
    if adaptive:
        lc = tr.step_controller.__dict__.get('_lte_c') if tr.step_controller is not None else None
        if lc is None or lc.served != now:
            return False
    got = _eligible(loop, tr, td, rec, adaptive)
    if type(got) is str:
        return _declined(loop, got)
    drv, O, tab, minbreak, rd, screen, judge = got
    n = rd.n
    run_ = _bound(loop, tr, rec, tab, rd, screen, judge)
    if run_ is None:
        return _declined(loop, 'bind')
    retry0 = _load(loop, tr, td, run_, rd, minbreak, adaptive, judge, O)
    s = run_.s
    status = _call(drv, rec, s, tr.statistics)
    if status not in (_CODES['LOOP_DONE'], _CODES['LOOP_FULL']):
        _no('loop_c:bail:' + BAIL.get(status, str(status)))
    if not s.attempts:
        loop.c.backoff = min(2 * loop.c.backoff + 1, BACKOFF_MAX)
        loop.c.skip = loop.c.backoff
        return False
    loop.c.backoff = 0
    _write_back(loop, tr, rec, run_, s, run_.buf, n, rd.hist_len, O, retry0)
    return True


#: the circuit's methods the attempts call, shadowed on the instance (a
#: perturbed source, as the shooting tests make): the Python attempts call
#: the shadow, the C the table and the core
_CIRCUIT_SHADOWS = ('u', 'G', 'C', 'i', 'q', 'limit', 'next_event', 'accept_step')


def _write_back(loop, tr, rec, run_, s, buf, n, hist_len, O, retry0):
    """The state after the C's attempts, as the Python attempts leave it
    (the module's THE STATE), and the `_paths` counts their glue would have
    made."""
    k, adaptive = s.k, bool(s.adaptive)
    Xo, To = buf['X'], buf['T']
    X, tl = loop.X, loop.timelist
    for i in range(k):
        X.append(Xo[i * n:(i + 1) * n].copy())
        tl.append(float(To[i]))
    stats = tr.statistics
    stats.accepted_steps += k
    stats.rejected_steps += s.rejected
    stats.newton_iterations += s.iters
    stats.breakpoints_hit += s.breaks
    stats.order_drops += s.euler_steps
    if s.have_hm:
        stats.min_step, stats.max_step = float(s.hmin), float(s.hmax)
    ## the last attempt's write-back (`_tran_newton_c`'s), accepted or not
    x = buf['xatt'].copy()
    cur = s.cur
    Cj = buf[f'cC{cur}'].reshape(n, n).copy()
    tr._q_cache = (x, buf['qatt'].copy())
    tr._C_cache = (x, Cj)
    if rec.core.npat:
        tr._C_pattern = (Cj, rec.core.arrays['pat'])
    tr._iq = buf['iqatt'].copy()
    tr._Geq = buf[f'Geq{cur}'].reshape(n, n).copy()
    tr._Cmat = Cj
    h, h_last = float(s.att_h), float(s.att_h_last)
    base = tr.base_integrator
    active = (O['Euler']() if s.att_formula == _ABI['PCX_FORMULA_EULER']
              and type(base) is not O['Euler'] else base)
    tr.active_integrator = active
    tr._companion_coeffs = active.companion_coefficients(h, h_last)
    tr._effective_method = type(active).__name__
    tr._dt = h
    tr._pred_pending = (float(s.att_t), ()) if s.att_rejected else None
    if k:
        ## the accepted steps' rings, history and step record
        tr._qlast = buf['qring'].reshape(hist_len, n).copy()
        tr._iqlast = buf['iqring'].reshape(hist_len, n).copy()
        tr._pred_hist = [(float(buf['ht'][i]), buf['hx'][i * n:(i + 1) * n].copy())
                         for i in range(s.hk)]
        tr._dt_last = float(s.dt_last)
        tr._dt_last2 = float(s.dt_last2) if s.have_dt_last2 else None
        tr._no_history = False
    tr._is_first_step = bool(s.is_first)
    if run_.br is not None:
        tr._branch_cmax = float(s.cmax)
    if rec.klu is not None:
        rec.klu.numerics += s.attempts
    if adaptive:
        ctrl, ct = tr.step_controller, run_.ct
        if s.keeps_run and s.have_run:
            ctrl._ref_running = buf['run'].copy()
        ctrl.last_err = float(s.last_err)
        kept = (retry0 is not None, retry0[1] if retry0 else None, retry0[2] if retry0 else None)
        now_ = (bool(ct.have_retry), float(ct.retry_h) if ct.have_retry else None,
                float(ct.retry_err) if ct.have_retry else None)
        if not ct.have_retry:
            ctrl._retry = None
        elif now_ != kept or ct.retry_point != _POINT_OTHER:
            ## (a rejection the C made: its point the accepted one it left
            ## from -- the object in `X`, as `_judged` keeps it)
            ctrl._retry = (X[ct.retry_point], float(ct.retry_h), float(ct.retry_err))
        if ct.have_q:
            ctrl._q_carried = float(ct.q_carried)
        if run_.sk is not None:
            run_.sk.refactors += s.ltes
    loop.t, loop.h = float(s.t), float(s.h)
    loop.landing = bool(s.landing)
    loop.order_drop = bool(s.order_drop)
    loop.rejects, loop.point_retries = int(s.rejects), int(s.point_retries)
    loop.t_break = s.t_break
    loop.x_new, loop.J = x, None
    loop.h_next = loop.timestep if not adaptive else loop.h
    loop.clamped, loop.hold, loop.imposed = False, not adaptive, False
    _count(s, run_, rec, adaptive, k)


def _count(s, run_, rec, adaptive, k):
    """The `_paths` counts the Python attempts' calls make, made for the C's
    attempts -- each where the Python makes it.  An A/B of the loop compares
    them all (test_loop_c): a count a unit adds on its served path is added
    here too."""
    a = s.attempts
    _PC['loop_c:attempts'] += a
    _PC['loop_c:steps'] += k
    ## the C Newton's serve, an attempt each (`_tran_newton_c.solve`, with the
    ## Python's own `chord:chosen`)
    if rec.klu is not None:
        _PC['newton_c:sparse'] += a
    _PC['newton_c:served'] += a
    _PC['newton_c:iterations'] += s.iters
    _PC['newton_c:jfold'] += a
    if not adaptive:
        _PC['newton_c:j_unread'] += a
    if rec.s.chord:
        _PC['chord:chosen'] += a
    ## the seed (`_StagePredictor._predict_fast`), the sources
    ## (`_stamp_plan.source_direct`, its `u` calls, the source memo
    ## `_CompanionModel._source_at`), the charge memo (`_CompanionModel._q_at`:
    ## the accepted point's and, adaptive, the error test's reads)
    if s.pred_fast:
        _PC['pred:fast'] += s.pred_fast
    _PC['src.direct:served'] += a
    _PC['u:called'] += len(run_.tab.entries) * a
    _PC['umemo:miss'] += a
    _PC['memo.q:same'] += k + (s.ltes if adaptive else 0)
    ## the branch screen (`_tran_branch_c.screen`, its C from
    ## `_CompanionModel._C_lookup`)
    if s.screens:
        _PC['branch_c:served'] += s.screens
        _PC['memo.C:equal'] += s.screens
        if run_.npat:
            _PC['branch_c:pattern'] += s.screens
    ## the error test (`_tran_lte_c.max_error`)
    if adaptive and s.ltes:
        _PC['lte_c:served'] += s.ltes
        if s.lte_euler:
            _PC['lte_c:euler'] += s.lte_euler
        if run_.sk is not None:
            _PC['lte_c:sparse_served'] += s.ltes
