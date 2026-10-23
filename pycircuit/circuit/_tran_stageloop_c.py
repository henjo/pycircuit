"""THE PLAIN RUN'S RADAU STEPS IN C: a stepping loop's attempts under the
stage family (`_SteppingLoop` over `_StageSteps`: a Radau IIA(3) step, the
dense coupled Newton `_rk_step_coupled` makes, the cost transform off) in
C -- `csrc/loop.c`'s stage family: where (the breakpoints, the run's end),
take (the stage times from the step's end, the start's record, the seeds by
pred.c, the sources at the stage times, radau.c on the transient's record
and its stage KLU, each stage's screen, the sources at the end, the
embedded estimate with numpy's dgemv and dgesv), judge (lte.c's STAGE
judgement, `_StageSteps.judge`'s law, the loop's budget) and accept (the
rings, the stages and the point noted) -- bit for bit (pinned in
`tests/test_pinned_pairs.py`).

ENTRY: only where the C stage Newton's last call served at the watch
counter as it stands (`_tran_radau_c`'s `served`), its record's checks
made; then the C loop's checks of the run (`_tran_loop_c`'s: its pieces,
the statistics, the circuit, the sources' table and the breakpoints, the
screen) and the family's own: `_StageSteps` with a Radau IIA(3)
integrator, no transform (whose step alone reads the early abort's arm),
no PCNR, no excursion veto, no stateful limiter, the stage pieces genuine
and not shadowed on the transient, numpy's underflow ignored (not watched:
the C stage Newton's check made again), the core's passes as stamped
(`_tran_core.passes_held`), numpy's sum the pairwise sum the C makes
(checked once a process), the loop's own output set within the C stage
Newton's budget for its sets (`_sets_fit`).  A decline is counted
(`stageloop_c:<why>`) and, where the reason is the run's, kept until
something watched moves; a stop at an attempt leaves that attempt to the
Python before the next entry.

THE STATE.  Read before a call -- the loop's (the point, the step, its
flags and counts), the transient's (the rings, the predictor's nodes, the
step record, the screen's maximum), the family's (its running reference,
its retry record), the device memo where it holds the start (the first
attempt's record), the stage KLU's numeric -- and written back after it as
the Python attempts leave them: the loop's and the statistics' as the C
loop writes a multistep run's (`_tran_loop_c`), the last made attempt's
stage outputs (`_finish_stage_step`'s: `_rk_Y`, `_q_cache`, `_iq`, `_Cmat`,
`_Geq`, the pending nodes of a rejected attempt; `_rk_est`; the pattern
mark), the device memo's two generations (the last made attempt's set; the
one before rebuilt by the core's passes at its stages, the very calls its
Newton made -- or the memo the call began with), the stage KLU's numeric
and counts, and the counts the Python's calls make.

A STOP before an attempt's Newton hands it back, nothing of it made; after
its converged Newton (a screen not quiet, a source, the estimate, the
judgement, the budget spent) the Newton's result is held for the Python's
repeat of the attempt (`_tran_radau_c.hold`), its numeric moved once.

THE OUTPUT SETS.  The C stage Newton writes an attempt's stages and their
records into one of two sets, the last made attempt's kept in the other:
the radau record's own (each call that writes it copies from it, no record
views it) and one of the loop's own.  The memo written back takes copies.

OWNS: on the `_SteppingLoop`'s C state (`_tran_loop_c.LoopState`) the
bound run (`_StageRun`).  `_paths`: `stageloop_c:` (the declines, the bails
`stageloop_c:bail:<why>`, `stageloop_c:attempts`, `stageloop_c:steps`).
Switch: `PYCIRCUIT_STAGELOOP_C=0` (`ENABLED`, read on every entry).

History: `doc/transient_history.md`, `_tran_stageloop_c.py`.
"""
import math
import os
import time
from types import SimpleNamespace

import numpy as np

from pycircuit.circuit import _paths, _tran_loop_c, _watch

#: `PYCIRCUIT_STAGELOOP_C=0` makes the Python attempts every one
ENABLED = os.environ.get('PYCIRCUIT_STAGELOOP_C', '1') != '0'
_PC = _paths.COUNTS
_no = _paths.no
_F64 = np.dtype(np.float64)
_LC = _tran_loop_c
_CODES = _LC._CODES
BAIL = _LC.BAIL
#: the stops after an attempt's converged Newton (where `held` says so):
#: its result held for the Python's repeat
_AFTER = frozenset((_CODES['LOOP_BAIL_SCREEN'], _CODES['LOOP_BAIL_SOURCES'],
                    _CODES['LOOP_BAIL_ESTIMATE'], _CODES['LOOP_BAIL_LTE'],
                    _CODES['LOOP_BAIL_FORCE']))
_RETRY_ORDER = {'assumed': _CODES['LOOP_RETRY_ASSUMED'], 'measured': _CODES['LOOP_RETRY_MEASURED'],
                'carried': _CODES['LOOP_RETRY_CARRIED']}
_MOD = {}
#: the reasons kept until something watched moves: the run's
_KEEPS = {
    'driver': True, 'patched': True, 'stats': True, 'loop': True, 'shadow': True,
    'events': True, 'state': True, 'circuit': True, 'veto': True, 'family': True,
    'pcnr': True, 'transform': True, 'est': True, 'units': True, 'sum': True,
    'sources': True, 'table': True, 'breakpoints': True, 'minbreak': True, 'screen': True,
    'judge': True, 'bind': True,
    'passes': False, 'table_stale': False, 'toolkit': False, 'rings': False,
    'history': False, 'underflow': False, 'size': True,
}

_M = 'pycircuit.circuit.'
#: the pieces the stage family's attempts call, beyond the C loop's
#: (`_tran_loop_c._PIECES`): (where read -- 'T' the transient's class, 'F'
#: the family's, a module -- the attribute, its qualified name, its module)
_S_PIECES = (
    ('T', '_solve_timestep_rk', '_SequentialStages._solve_timestep_rk', '_tran_stages'),
    ('T', '_stage_source', '_SequentialStages._stage_source', '_tran_stages'),
    ('T', '_finish_stage_step', '_SequentialStages._finish_stage_step', '_tran_stages'),
    ('T', '_rk_step_coupled', '_RadauStages._rk_step_coupled', '_tran_radau'),
    ('T', '_coupled_stage_solver', '_RadauStages._coupled_stage_solver', '_tran_radau'),
    ('T', '_coupled_stage_context', '_RadauStages._coupled_stage_context', '_tran_radau'),
    ('T', '_finish_radau', '_RadauStages._finish_radau', '_tran_radau'),
    ('T', '_radau_error_estimate', '_RadauStages._radau_error_estimate', '_tran_radau'),
    ('T', '_branch_after_coupled', '_BranchCheck._branch_after_coupled', '_tran_branch'),
    ('T', '_pred_stages', '_StagePredictor._pred_stages', '_tran_predictor'),
    ('T', '_memo_step', '_CompanionModel._memo_step', '_tran_companion'),
    ('T', '_memo_get', '_CompanionModel._memo_get', '_tran_companion'),
    ('T', '_memo_put', '_CompanionModel._memo_put', '_tran_companion'),
    ('T', '_memo_ok', '_CompanionModel._memo_ok', '_tran_companion'),
    ('T', '_rk_use_pcnr', '_PCNRSteps._rk_use_pcnr', '_tran_pcnr'),
    ('T', '_newton_option', '_StepNewton._newton_option', '_tran_newton'),
    ('T', '_lte_mask', '_RunHistory._lte_mask', '_tran_history'),
    ('F', 'judge', '_StageSteps.judge', 'transient'),
    ('F', 'normalised', '_StageSteps.normalised', 'transient'),
    ('F', '_reference', '_StageSteps._reference', 'transient'),
    ('F', '_retry_after', '_StageSteps._retry_after', 'transient'),
    ('F', '_retry_order', '_StageSteps._retry_order', 'transient'),
    ('F', 'attempt', '_StepFamily.attempt', 'transient'),
    ('F', 'next_breakpoint', '_StepFamily.next_breakpoint', 'transient'),
    ('F', 'h_after_force', '_StepFamily.h_after_force', 'transient'),
    ('F', 'after_accept', '_StepFamily.after_accept', 'transient'),
    ('stepcontroller', 'charge_floor', 'charge_floor', 'stepcontroller'),
    ('stepcontroller', 'normalised_error', 'normalised_error', 'stepcontroller'),
    ('stepcontroller', 'sigglobal_reference', 'sigglobal_reference', 'stepcontroller'),
    ('_tran_radau_c', 'solve', 'solve', '_tran_radau_c'),
    ('_tran_radau_c', '_ready', '_ready', '_tran_radau_c'),
    ('_tran_core', 'passes', 'passes', '_tran_core'),
    ('_tran_pred_c', 'stages', 'stages', '_tran_pred_c'),
    ('_tran_branch_c', 'screen', 'screen', '_tran_branch_c'),
    ('_numeric', 'linearsolver', 'linearsolver', '_numeric'))
#: an instance attribute of these names on the family sends its attempt elsewhere
_F_NAMES = tuple(a for o, a, _q, _m in _S_PIECES if o == 'F') + ('SAFETY', 'K')
#: ... and of these on the transient: the pieces, and what the C stage
#: Newton's own checks refuse shadowed (`_tran_radau_c._ready`)
_T_NAMES = tuple(a for o, a, _q, _m in _S_PIECES if o == 'T') + (
    '_source_at', '_coupled_stage_system', '_stages_converged', '_stage_limiter_states',
    '_limiters_at_rest', '_newton_limiter')


def _owners():
    """The objects the pieces are read from, imported once."""
    if 'owners' not in _MOD:
        import importlib

        from pycircuit.circuit.transient import Transient, _StageSteps
        own = {'T': Transient, 'F': _StageSteps}
        for o, _a, _q, _m in _S_PIECES:
            if o not in own:
                own[o] = importlib.import_module(_M + o)
        _MOD['owners'] = own
    return _MOD['owners']


def _now():
    own = _owners()
    return tuple(getattr(own[o], a) for o, a, _q, _m in _S_PIECES)


def _pieces():
    """The pieces as their modules define them -- read once every one is
    (`_paths.genuine`); None: one is a stand-in, read again next entry."""
    g = _MOD.get('g')
    if g is None:
        own = _owners()
        objs = [(getattr(own[o], a), q, _M + m) for o, a, q, m in _S_PIECES]
        if not all(_paths.genuine(*o) for o in objs):
            return None
        g = _MOD['g'] = tuple(o[0] for o in objs)
    return g


def _declined(loop, why):
    """The entry declined, counted -- kept until something watched moves
    where the reason is the run's."""
    if _KEEPS[why]:
        loop.c.kept = _watch.now()
    _no('stageloop_c:' + why)
    return False


def run(loop):
    """Make as many of `loop`'s attempts in C as it can: True where it made
    any (the state written back as the Python attempts leave it), else False
    -- the Python loop makes the next attempt (`_tran_loop_c.run`'s entry
    for the stage family)."""
    if not ENABLED:
        return False
    c = loop.c
    if c.skip:
        c.skip -= 1
        return False
    now = _watch.now()
    if now == -1 or c.kept == now:
        return False
    tr = loop.tr
    td = tr.__dict__
    rec = td.get('_radau_c')
    ## (the C stage Newton's last call served, nothing watched moved since:
    ## its record's checks made)
    if rec is None or rec.served != now:
        return False
    e = _eligible(loop, tr, td, rec)
    if e.__class__ is str:
        return _declined(loop, e)
    run_ = _bound(loop, e)
    if run_ is None:
        return _declined(loop, 'bind')
    _load(loop, e, run_)
    status = _call(e, run_)
    s = run_.s
    if status not in (_CODES['LOOP_DONE'], _CODES['LOOP_FULL']):
        _no('stageloop_c:bail:' + BAIL.get(status, str(status)))
    _settle(e, run_, status)
    if not s.attempts:
        c.backoff = min(2 * c.backoff + 1, _LC.BACKOFF_MAX)
        c.skip = c.backoff
        return False
    c.backoff = 0
    ## (a stop at an attempt -- nothing of it made, or its Newton's result
    ## held -- leaves that attempt to the Python before the next entry: an
    ## entry now would make it again, on the stage KLU from the numeric its
    ## first Newton moved)
    c.skip = 0 if status in (_CODES['LOOP_DONE'], _CODES['LOOP_FULL']) else 1
    _write_back(loop, e, run_)
    return True


def _eligible(loop, tr, td, rec):
    """What an entry runs with (a namespace) -- or why not: the run
    (`_run_ready`), the family (`_family_ready`), the sources
    (`_tran_loop_c._sources`), the size (`_sets_fit`), the state
    (`_state_ready`), the screen (`_tran_loop_c._screen`), the judge's
    settings (`_judge_ready`)."""
    drv = _LC.driver()
    if drv is None:
        return 'driver'
    O = _LC._OWN or _LC._own()
    why = _run_ready(loop, tr, td, O) or _family_ready(loop, tr, td, rec)
    if why:
        return why
    got = _LC._sources(loop, tr, O)
    if got.__class__ is str:
        return got
    tab, minbreak = got
    n = rec.n
    if not _sets_fit(n):
        return 'size'
    st = _state_ready(loop, tr, td, n)
    if st.__class__ is str:
        return st
    screen = _LC._screen(tr, td, rec, n)
    if screen is None:
        return 'screen'
    cfg = _judge_ready(loop, tr, n) if not loop.fixed else None
    if cfg.__class__ is str:
        return cfg
    return SimpleNamespace(drv=drv, O=O, tr=tr, rec=rec, n=n, iref=rec.iref, tab=tab,
                           minbreak=minbreak, st=st, screen=None if screen[0] is None else screen,
                           cfg=cfg, sk=rec.klu, fam=loop.family, integ=tr.base_integrator)


def _sets_fit(n):
    """The loop's own output set within the C stage Newton's budget for its
    sets (`_tran_radau_c.OUT_BYTES`): above it an attempt's dense estimate
    outweighs what the loop saves."""
    from pycircuit.circuit import _tran_radau_c
    return 2 * _tran_radau_c.MAXA * 3 * (n * n + 2 * n) * 8 <= _tran_radau_c.OUT_BYTES


def _run_ready(loop, tr, td, O):
    """Why the run is not one the C makes (a reason), or '': the C loop's
    checks of the run but the multistep family's own (`_tran_loop_c.
    _why_not_the_run`), and the excursion veto."""
    if O is None or type(tr) is not O['T'] or _LC._pieces(O) != O['chain']:
        return 'patched'
    if type(tr.statistics) is not O['stats']:
        return 'stats'
    if type(loop) is not O['Loop'] or loop.provided_function is not None:
        return 'loop'
    if any(k in td for k in _LC._SHADOWS):
        return 'shadow'
    if getattr(tr, '_ev_rows', None) is not None or getattr(tr, '_periodic_rows', None):
        return 'events'
    if td.get('_stateful_lims'):
        return 'state'
    cir = tr.cir
    if (type(cir).next_event is not O['Sub'].next_event
            or type(cir).accept_step is not O['Sub'].accept_step
            or type(cir).u is not O['Sub'].u
            or any(k in cir.__dict__ for k in _LC._CIRCUIT_SHADOWS)
            or cir._hook_elements('accept_step')):
        return 'circuit'
    if tr.par.max_dv_step is not None or tr.par.max_di_step is not None:
        return 'veto'
    return ''


def _family_ready(loop, tr, td, rec):
    """Why the attempts are not the stage family's the C makes (a reason),
    or '': the family and integrator, the pieces shadowed on the transient,
    numpy's underflow, the step's options, the pieces, the units, the core's
    passes, numpy's sum."""
    from pycircuit.circuit import (
        _tran_core,
        _tran_lte_c,
        _tran_pred_c,
        _tran_radau_c,
        func,
    )
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    fam = loop.family
    if (type(fam) is not _owners()['F'] or type(tr.base_integrator) is not RadauIIA3Integrator
            or any(k in fam.__dict__ for k in _F_NAMES)):
        return 'family'
    if any(k in td for k in _T_NAMES):
        return 'shadow'
    ## (numpy's error state is not watched: the C stage Newton's check of it
    ## made again)
    if not func._underflow_ignored():
        return 'underflow'
    if tr._rk_use_pcnr():
        return 'pcnr'
    ## (the early abort's arm is read by the transform's step alone: the
    ## coupled step's is set and taken back unread, `_SteppingLoop.take`)
    if tr._newton_option(tr.par.radau_transform, 'radau_transform'):
        return 'transform'
    if bool(getattr(tr, '_rk_want_est', False)) != (not loop.fixed):
        return 'est'
    g = _MOD.get('g') or _pieces()
    if g is None or _now() != g:
        return 'patched'
    if not (_tran_radau_c.ENABLED and _tran_pred_c.ENABLED and _tran_pred_c.driver()
            and _tran_lte_c.driver() and _tran_radau_c.driver()):
        return 'units'
    if not _tran_core.passes_held(tr, rec.core):
        return 'passes'
    if not _sum_ok():
        return 'sum'
    return ''


def _state_ready(loop, tr, td, n):
    """`(the rings, the predictor's nodes, X, the last point, the step
    record)` as the Python leaves them -- or why not."""
    ql, iql = td.get('_qlast'), td.get('_iqlast')
    if not (type(ql) is np.ndarray and ql.dtype is _F64 and ql.ndim == 2 and ql.shape[1] == n
            and ql.shape[0] >= 1 and type(iql) is np.ndarray and iql.dtype is _F64
            and iql.shape == ql.shape):
        return 'rings'
    hist = td.get('_pred_hist') or []
    if (type(hist) is not list or len(hist) > _LC._HIST_MAX
            or not all(type(e) is tuple and len(e) == 2 and type(e[0]) is float
                       and type(e[1]) is np.ndarray and e[1].dtype is _F64 and e[1].shape == (n,)
                       for e in hist)):
        return 'history'
    X = loop.X
    xl = X[-1] if X else None
    dl, dl2 = tr._dt_last, tr._dt_last2
    if not (type(xl) is np.ndarray and xl.dtype is _F64 and xl.shape == (n,)
            and (dl is None or type(dl) is float) and (dl2 is None or type(dl2) is float)
            and type(loop.t) is float and type(loop.h) is float):
        return 'history'
    return ql, iql, hist, X, xl, dl, dl2


def _judge_ready(loop, tr, n):
    """The stage judgement's settings as `_StageSteps.normalised` reads
    them -- the reference, the absolute tolerance per unknown, the mask, the
    charge floor, the family's running reference and retry record -- or why
    not (`judge`)."""
    from pycircuit.circuit import _tran_lte_c
    lc = _tran_lte_c._CODES
    fam = loop.family
    ## (`_StageSteps._reference`: any other name the group's)
    relref = ({'pointlocal': lc['LTE_RELREF_POINTLOCAL'], 'alllocal': lc['LTE_RELREF_ALLLOCAL']}
              .get(fam.relref, lc['LTE_RELREF_SIGGLOBAL']) if type(fam.relref) is str else None)
    order = _RETRY_ORDER.get(fam.retry_order)
    nn = fam.n_nodes
    abstol = np.asarray(loop.run.abstol, dtype=float)
    mask, charge = tr._lte_mask(), tr._lte_charge()
    running, retry, qc = fam.running, fam._retry, fam._q_carried
    reltol, chgtol = tr.par.reltol, tr.par.lte_chgtol
    if (relref is None or order is None or not (nn is None or type(nn) is int)
            or abstol.shape not in ((), (n,)) or type(reltol) is not float
            or (mask is not None and np.shape(mask) != (n,))
            or (charge is not None and (np.shape(charge) != (n,) or type(chgtol) is not float))
            or not (running is None or (type(running) is np.ndarray and running.dtype is _F64
                                        and running.shape == (n,)))
            or not (qc is None or type(qc) is float)
            or not (retry is None or (type(retry[1]) is float and type(retry[2]) is float))):
        return 'judge'
    return SimpleNamespace(relref=relref, order=order, nn=-1 if nn is None else nn,
                           abstol=np.broadcast_to(abstol, (n,)).copy(), mask=mask, charge=charge,
                           reltol=reltol, chgtol=chgtol, running=running, retry=retry, qc=qc)


def _sum_ok():
    """numpy's `sum` of a float64 vector the pairwise sum lte.c makes (+0.0
    first), checked once a process on vectors of every blocking: the C's
    RMS through its own entry against numpy's."""
    got = _MOD.get('sum')
    if got is None:
        from pycircuit.circuit import _tran_lte_c
        ffi, cfn, dgesv = _tran_lte_c.driver()
        rng = np.random.default_rng(20261023)
        got = True
        for n in (2, 3, 8, 9, 17, 129, 130, 300, 1031):
            c = _tran_lte_c._Ctx(ffi, dgesv, n, 0)
            s, b = c.s, c.buf
            ek = rng.random(n - 1) * 10.0 ** rng.integers(-3, 3, n - 1)
            lc = _tran_lte_c._CODES
            s.formula, s.h1, s.reltol, s.trtol, s.abstol = lc['LTE_FORMULA_STAGE'], 1.0, 0.0, 1.0, 1.0
            s.relref, s.nn, s.abs_vec = lc['LTE_RELREF_POINTLOCAL'], -1, 0
            s.have_run, s.has_mask, s.has_charge = 0, 0, 0
            b['xc'][:] = 0.0
            b['xl'][:] = 0.0
            b['lte'][0], b['lte'][1:] = 0.0, ek
            want = float((np.sum(ek * ek) / (n - 1)) ** 0.5)
            if cfn(s) != _tran_lte_c._CODES['LTE_DONE'] or float(s.err) != want:
                got = False
                break
        _MOD['sum'] = got
    return got


class _StageRun(_LC._Run):
    """The stage family's bound run: the C loop's (`_tran_loop_c._Run`: the
    struct, its buffers, the predictor's state, the screen's) and its own --
    the two output sets, the stage buffers, the error test's context, the
    binding's key."""

    __slots__ = ('key', 'lctx', 'sbuf', 'sets', 'skeep')


def _bound(loop, e):
    """The loop's bound stage run, bound again where the record, the table,
    the screen's presence or certificate, the stage KLU or the sizes moved;
    None where it cannot be bound."""
    run_ = loop.c.run
    scr = e.screen
    key = (e.rec, e.tab, e.sk, None if scr is None else scr[3],
           None if scr is None else scr[4], e.n, e.st[0].shape[0])
    if not (type(run_) is _StageRun and len(run_.key) == len(key)
            and all(a is b or (type(a) is float and a == b) or (type(a) is int and a == b)
                    for a, b in zip(run_.key, key))):
        run_ = _bind(e)
        if run_ is not None:
            run_.key = key
            loop.c.run = run_
    return run_


def _bind(e):
    """The stage run for this entry's record, table and screen: the C loop's
    buffers, the predictor's state, the sources' table, the screen, the
    radau record, an error test of its own (lte.c), the stage buffers."""
    from pycircuit.circuit import _tran_lte_c, _tran_pred_c, _tran_radau_c
    tr, rec, n, iref = e.tr, e.rec, e.n, e.iref
    if int(rec.s.maxa) != _tran_radau_c.MAXA:
        ## (the record's slots and the sets made here: one count)
        return None
    hist_len = e.st[0].shape[0]
    ffi, _cfn = _LC.driver()
    run = _StageRun()
    run.ffi, run.nrec, run.tab = ffi, rec, e.tab
    s = run.s = ffi.new('loop_t *')
    _LC._bind_buffers(run, n, hist_len)
    pffi, pcfn, pdgesv, pdgemv = _tran_pred_c.driver()
    ps, pbuf, pkeep = _tran_pred_c.new_state(pffi, pdgesv, pdgemv, n)
    ps.ht = pffi.cast('double *', _LC._addr(ffi, s.ht))
    ps.hx = pffi.cast('double *', _LC._addr(ffi, s.hx))
    ps.deg = int(tr._pred_degree())
    ps.clamp = float(tr.PRED_CLAMP)
    run.pr = (ps, pbuf, pkeep)
    _LC._plug(s, ffi, 'pred', pffi, ps, pcfn)
    _LC._plug(s, ffi, 'src', e.tab.ffi, e.tab.s, e.tab.cfn)
    _LC._bind_screen(run, tr, n, e.screen)
    s.newton = s.newton_fn = s.ctrl = s.ctrl_fn = ffi.NULL
    rffi, rcfn, rdgesv, rdgemv = _tran_radau_c.driver()
    _LC._plug(s, ffi, 'radau', rffi, rec.s, rcfn)
    lffi, lcfn, ldgesv = _tran_lte_c.driver()
    run.lctx = _tran_lte_c._Ctx(lffi, ldgesv, n, iref)
    _LC._plug(s, ffi, 'lte', lffi, run.lctx.s, lcfn)
    run.lte = run.ct = run.sk = None
    _bind_stage(run, rec, n, rdgesv, rdgemv)
    s.n, s.hist_len, s.cap, s.iref = n, hist_len, _LC.CAP_MAX, iref
    s.family = _CODES['LOOP_FAMILY_STAGE']
    return run


def _bind_stage(run, rec, n, dgesv, dgemv):
    """The stage buffers -- the two output sets (`MAXA` slots: the
    record's own, which each call that writes it copies from and no record
    views, and one of the loop's own), the sources, seeds, the start's
    passes, the estimate's work, the older generation's stages, the held
    charge, the made attempt's `iq` and estimate -- each pointed at by its
    field and kept."""
    from pycircuit.circuit import _tran_radau_c
    ffi, s, m = run.ffi, run.s, n - 1
    run.sets = sets = [rec.own, rec._new_out()]
    run.sbuf = sb = {
        'ust': np.empty(5 * n), 'useed': np.empty(3 * n), 'xr_q': np.empty(n),
        'xr_C': np.empty(n * n), 'xr_G': np.empty(n * n), 'xr_i': np.empty(n), 'eF': np.empty(n),
        'eA': np.empty(max(m * m, 1)), 'eb': np.empty(max(m, 1)),
        'ygen1': np.empty(_tran_radau_c.MAXA * 3 * n), 'hq': np.empty(n), 'liq': np.empty(n),
        'lest': np.empty(n)}
    fb = ffi.from_buffer
    run.skeep = keep = []
    for k, a in sb.items():
        keep.append(fb('double *', a))
        setattr(s, k, keep[-1])
    ipiv = np.empty(max(m, 1), dtype=np.int64)
    sb['eipiv'] = ipiv
    keep.append(fb('int64_t *', ipiv))
    s.eipiv = keep[-1]
    for o in range(2):
        for k, nm in enumerate(_tran_radau_c._OUTS):
            keep.append(fb('double *', sets[o][0][nm]))
            s.oset[o][k] = keep[-1]
    s.dgesv, s.dgemv = ffi.cast('void *', dgesv), ffi.cast('void *', dgemv)


def _load(loop, e, run_):
    """The struct loaded with the state as the Python leaves it: the loop's,
    the transient's rings and nodes, the run's constants, the stage's, the
    judge's."""
    tr, n = e.tr, e.n
    s, buf = run_.s, run_.buf
    ql, iql, hist, X, xl, dl, dl2 = e.st
    buf['qring'][:] = ql.reshape(-1)
    buf['iqring'][:] = iql.reshape(-1)
    buf['xlast'][:] = xl
    s.hk = len(hist)
    for i, (tt, xx) in enumerate(hist):
        buf['ht'][i] = tt
        buf['hx'][i * n:(i + 1) * n] = xx
    s.cmax = tr.__dict__.get('_branch_cmax', 0.0)
    if run_.br is not None:
        run_.br[0].ncert = 0
    s.t, s.h = loop.t, loop.h
    s.landing, s.is_first = int(bool(loop.landing)), int(bool(tr._is_first_step))
    s.order_drop, s.imposed = int(bool(loop.order_drop)), int(bool(loop.imposed))
    s.rejects, s.point_retries, s.npoint = loop.rejects, loop.point_retries, len(X)
    s.dt_last, s.have_dt_last = (dl, 1) if dl is not None else (0.0, 0)
    s.dt_last2, s.have_dt_last2 = (dl2, 1) if dl2 is not None else (0.0, 0)
    s.tend, s.timestep, s.minbreak = float(loop.tend), float(loop.timestep), e.minbreak
    s.adaptive = int(not loop.fixed)
    s.max_step, s.minstep = float(loop.run.max_step), float(loop.minstep)
    s.max_retries, s.max_reject = e.fam.max_retries, e.fam.max_reject
    s.early_th, s.newton_reject = 0.0, 0
    stats = tr.statistics
    s.have_hm = int(stats.min_step is not None)
    s.hmin = stats.min_step if stats.min_step is not None else 0.0
    s.hmax = stats.max_step if stats.max_step is not None else 0.0
    _load_stage(e, run_, xl)
    if e.cfg is not None:
        _load_judge(e, run_, X)
    else:
        s.keeps_run = s.have_run = 0


def _load_stage(e, run_, xl):
    """The stage's: the tableau and the estimate's constants (as the Python
    computes them), the generations (none made yet), the start's record
    where the memo the Python left holds it, the passes' bits, the stage
    KLU's numeric, numpy's underflow state."""
    from pycircuit.circuit import func
    s, sb, ffi, rec = run_.s, run_.sbuf, run_.ffi, e.rec
    integ = e.integ
    cvec, Amat = np.array(integ.C, dtype=float), np.array(integ.A, dtype=float)
    for i in range(3):
        s.ct[i] = cvec[i]
    s.a33 = Amat[2, 2]
    s6 = math.sqrt(6.0)
    s.gamma_r = integ.GAMMA_REAL
    s.dd1, s.dd2, s.dd3 = -(13.0 + 7.0 * s6) / 3.0, (-13.0 + 7.0 * s6) / 3.0, -1.0 / 3.0
    s.gen0, s.A0, s.A1, s.made = -1, 0, 0, 0
    memo = e.tr.__dict__.get('_dev_memo') or ({}, {})
    e.memo0 = memo[0]
    got = memo[0].get(xl.tobytes()) if memo[0] else None
    s.first_rec = int(got is not None and all(k in got for k in ('q', 'i', 'C', 'G')))
    if s.first_rec:
        sb['xr_q'][:] = got['q']
        sb['xr_i'][:] = got['i']
        sb['xr_C'][:] = np.asarray(got['C']).reshape(-1)
        sb['xr_G'][:] = np.asarray(got['G']).reshape(-1)
    abi = _LC._ABI
    npat = bool(rec.core.npat)
    s.pass_bits = abi['PCX_WANT_ALL'] | (abi['PCX_WANT_PATTERN'] if npat else 0)
    s.pass_zero = int(npat)
    num = e.sk._numeric if e.sk is not None else None
    s.klu_cur = ffi.cast('void *', num) if num else ffi.NULL
    s.under_ok = int(func._underflow_ignored())
    run_.pr[0].under_ok = s.under_ok
    run_.lctx.s.under_ok = s.under_ok


def _load_judge(e, run_, X):
    """The judge's: the error test's settings (lte.c's STAGE run fields),
    the running reference, the retry record (its point an index into X, as
    the C loop keeps the controller's), the law's constants."""
    from pycircuit.circuit import _tran_lte_c
    from pycircuit.circuit.stepcontroller import MIN_SHRINK_RATIO, RETRY_ORDER_MIN
    s, buf, cfg = run_.s, run_.buf, e.cfg
    ls, lb = run_.lctx.s, run_.lctx.buf
    lc = _tran_lte_c._CODES
    ls.formula, ls.trtol, ls.abstol, ls.reltol = lc['LTE_FORMULA_STAGE'], 1.0, 0.0, cfg.reltol
    ls.relref, ls.nn, ls.abs_vec = cfg.relref, cfg.nn, 1
    lb['abst'][:] = cfg.abstol
    ls.has_mask = int(cfg.mask is not None)
    if cfg.mask is not None:
        lb['mask'][:] = np.asarray(cfg.mask, dtype=float)
    ls.has_charge = int(cfg.charge is not None)
    if cfg.charge is not None:
        lb['charge'][:] = np.asarray(cfg.charge, dtype=float)
        ls.chgtol = cfg.chgtol
    s.keeps_run = int(cfg.relref != lc['LTE_RELREF_POINTLOCAL'])
    s.have_run = int(cfg.running is not None)
    if cfg.running is not None:
        buf['run'][:] = cfg.running
    r = cfg.retry
    s.have_sretry = int(r is not None)
    if r is not None:
        s.sretry_point = len(X) - 1 if r[0] is X[-1] else _LC._POINT_OTHER
        s.sretry_h, s.sretry_worst = r[1], r[2]
    s.have_sq, s.sq = (1, cfg.qc) if cfg.qc is not None else (0, 0.0)
    s.retry_order = cfg.order
    s.safety, s.kclamp = float(e.fam.SAFETY), float(e.fam.K)
    s.pst = float(int(e.integ.EMBEDDED_ORDER) + 1)
    s.min_shrink, s.rq_min = float(MIN_SHRINK_RATIO), float(RETRY_ORDER_MIN)


def _call(e, run_):
    """The C's attempts: its status -- the radau record's own pointers saved
    before and restored after (the C points them at the loop's buffers),
    the time the solve's."""
    rs = e.rec.s
    saved = (rs.qn, rs.u, rs.seed, rs.Ya, rs.qa, rs.ia, rs.Ca, rs.Ga)
    t0 = time.perf_counter()
    try:
        return e.drv[1](run_.s)
    finally:
        rs.qn, rs.u, rs.seed, rs.Ya, rs.qa, rs.ia, rs.Ca, rs.Ga = saved
        e.tr.statistics.solve_seconds += time.perf_counter() - t0


def _settle(e, run_, status):
    """The stage KLU as the attempts leave it -- its numeric, the made
    attempts' counts -- and where the loop stopped after an attempt's
    converged Newton, its result held for the Python's repeat of it
    (`_tran_radau_c.hold`: the inputs it was made from, its set's slots)."""
    from pycircuit.circuit import _tran_radau_c
    s, sk, rec, n = run_.s, e.sk, e.rec, e.n
    if (status == _CODES['LOOP_BAIL_RADAU']
            and rec.s.status == _tran_radau_c._CODES['RADAU_BAIL_WALKSTOP']):
        ## (the limiter walk stopped: its tables are taken again, as `solve`
        ## drops them)
        rec.walk_ok = None
    if sk is not None:
        sk._numeric = int(run_.ffi.cast('uintptr_t', s.klu_cur)) or None
        sk.refactors += s.st_klu_refactors
        sk.factors += s.st_klu_factors
        sk.pivot_fallbacks += s.st_klu_pivot
        sk.residual_fallbacks += s.st_klu_resid
    if s.held and status in _AFTER:
        sb, rs = run_.sbuf, rec.s
        o = 1 if s.gen0 == 0 else 0
        _tran_radau_c.hold(rec, sk, sb['hq'], sb['ust'][:3 * n].reshape(3, n),
                           sb['useed'].reshape(3, n), np.array(e.integ.A, dtype=float),
                           float(rs.T), float(s.cur_h), float(rs.reltol), float(rs.abstol),
                           int(rs.maxiter), outs=run_.sets[o][0])


def _write_back(loop, e, run_):
    """The state after the C's attempts as the Python attempts leave it (the
    module's THE STATE), and the counts their calls make."""
    s, buf, tr, n = run_.s, run_.buf, e.tr, e.n
    k = int(s.k)
    Xo, To = buf['X'], buf['T']
    X, tl = loop.X, loop.timelist
    for i in range(k):
        X.append(Xo[i * n:(i + 1) * n].copy())
        tl.append(float(To[i]))
    _LC._stats_back(tr.statistics, s, k)
    hist_len = e.st[0].shape[0]
    if k:
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
    x = _stage_back(e, run_)
    _family_back(e, run_, X)
    adaptive = bool(s.adaptive)
    loop.t, loop.h = float(s.t), float(s.h)
    loop.landing = bool(s.landing)
    loop.order_drop = bool(s.order_drop)
    loop.rejects, loop.point_retries = int(s.rejects), int(s.point_retries)
    loop.t_break = s.t_break
    loop.x_new, loop.J = x, None
    loop.h_next = loop.timestep if not adaptive else loop.h
    loop.clamped, loop.hold, loop.imposed = False, not adaptive, bool(s.imposed)
    _count(e, run_, k)


def _stage_back(e, run_):
    """The last made attempt's stage outputs as `_finish_stage_step` and the
    estimate leave them, the device memo's two generations, the pattern mark
    the C stage Newton's serve leaves: its point."""
    s, sb, tr, rec = run_.s, run_.sbuf, e.tr, e.rec
    _memo_back(e, run_)
    o, A = run_.sets[int(s.gen0)][0], int(s.A0)
    Y = [o['Ya'][A - 1, j].copy() for j in range(3)]
    g0 = tr._dev_memo[0]
    recs = [g0[Y[j].tobytes()] for j in range(3)]
    h = float(s.att_h)
    tr._rk_Y = Y
    tr._q_cache = (Y[2], recs[2]['q'])
    tr._iq = sb['liq'].copy()
    tr._Cmat = recs[2]['C']
    tr._Geq = float(s.a33) * h * recs[2]['G']
    tr._effective_method = type(e.integ).__name__
    tr._companion_coeffs = None
    tr._dt = h
    if s.adaptive:
        tr._rk_est = sb['lest'].copy()
    ts = [np.float64(s.att_ts[i]) for i in range(3)]
    tr._pred_pending = ((float(s.att_tend), list(zip(ts, Y))) if s.att_rejected else None)
    if rec.core.npat:
        tr._C_pattern = (tuple(r['C'] for r in recs), rec.core.arrays['pat'])
    if run_.br is not None:
        ## (the screen's last stage, as `_branch_screen` keeps it)
        tr._C_cache = (Y[2], recs[2]['C'])
    return Y[2]


def _memo_back(e, run_):
    """The device memo the Python attempts leave (`_memo_step` rolling,
    `_memo_put` per assembly): the last made attempt's records -- each stage
    of each assembly, the latest put where several share bytes -- and the
    one before: the attempt before it rebuilt from its stages (`_rebuild`),
    or, after one made attempt, the memo the call began with."""
    s, n = run_.s, e.n
    o = run_.sets[int(s.gen0)][0]
    gen0 = {}
    for a in range(int(s.A0)):
        for j in range(3):
            gen0[o['Ya'][a, j].tobytes()] = {
                'q': o['qa'][a, j].copy(), 'i': o['ia'][a, j].copy(),
                'C': o['Ca'][a, j].copy(), 'G': o['Ga'][a, j].copy()}
    gen1 = e.memo0 if s.made == 1 else _rebuild(e, run_.sbuf['ygen1'], int(s.A1), n)
    e.tr._dev_memo = (gen0, gen1)
    e.tr._memo_rolling = True


def _rebuild(e, ys, A, n):
    """The records an attempt's Newton put at its stages `ys` (`A`
    assemblies of three), made again by the very core calls its assemblies
    made (radau.c's: formula 4, every output, C its own; G and i from the
    core's bins) -- nothing counted, nothing else touched."""
    from pycircuit.circuit import _tran_core
    ffi, cfn, _dgemv = _tran_core.driver()
    core = e.rec.core
    fb, dp = ffi.from_buffer, core.dptr
    dummy = np.zeros(max(n * n, 1))
    d = fb(dp, dummy)
    T = float(e.rec.s.T)
    bits = _LC._ABI['PCX_WANT_ALL']
    out = {}
    for k in range(3 * A):
        x = ys[k * n:(k + 1) * n].copy()
        q, C = np.empty(n), np.empty((n, n))
        cfn(core.cs, fb(dp, x), T, bits, 4, 0.0, 0.0, 0.0, 1.0, 0.0, d, d, d, d, d, fb(dp, C),
            fb(dp, q), d, d, d, d)
        out[x.tobytes()] = {'q': q, 'i': core.arrays['bini'].copy(), 'C': C,
                            'G': core.arrays['binG'].reshape(n, n).copy()}
    return out


def _family_back(e, run_, X):
    """The family's state the judgements leave: the running reference, the
    retry record (its point the object in X), the carried order."""
    s, fam, cfg = run_.s, e.fam, e.cfg
    if cfg is None:
        return
    if s.keeps_run and s.have_run:
        fam.running = run_.buf['run'].copy()
    r = cfg.retry
    if not s.have_sretry:
        fam._retry = None
    elif (r is None or s.sretry_point != _LC._POINT_OTHER
          or (float(s.sretry_h), float(s.sretry_worst)) != (r[1], r[2])):
        fam._retry = (X[int(s.sretry_point)], float(s.sretry_h), float(s.sretry_worst))
    fam._q_carried = float(s.sq) if s.have_sq else None


def _count(e, run_, k):
    """The `_paths` counts the Python attempts' calls make, made for the C's
    -- the C stage Newton's serves, the sources' passes and the step memo,
    the seeds, the screens, the estimate, the core's passes at the start,
    the accepted point's charge memo."""
    from pycircuit.circuit import _tran_core
    s, a = run_.s, int(run_.s.attempts)
    _PC['stageloop_c:attempts'] += a
    _PC['stageloop_c:steps'] += k
    _PC['radau_c:served'] += a
    _PC['radau_c:assemblies'] += s.st_assemblies
    if s.st_klu_on:
        _PC['radau_c:on_klu'] += s.st_klu_on
    _PC['umemo:hit'] += 3 * (s.st_assemblies - a) + s.st_uhit
    _PC['umemo:miss'] += s.st_umiss
    _PC['src.direct:served'] += s.st_umiss
    _PC['u:called'] += len(e.tab.entries) * s.st_umiss
    if s.st_fast:
        _PC['pred:fast'] += s.st_fast
    if run_.br is not None:
        _PC['branch_c:served'] += 3 * a
        _PC['memo.C:memo'] += 3 * a
        if run_.npat:
            _PC['branch_c:pattern'] += 3 * a
        if run_.br[0].ncert:
            _PC['branch_c:cert'] += run_.br[0].ncert
    if s.adaptive:
        _PC['radau.est:fresh'] += a
        if s.st_estmemo:
            _PC['radau.est:memo'] += s.st_estmemo
    for _ in range(int(s.st_passes)):
        _tran_core.passes_count(e.tr, e.rec.core)
    _PC['memo.q:same'] += k
