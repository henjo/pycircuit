"""THE PLAIN RUN'S RADAU STEPS ON THE COST TRANSFORM IN C: a stepping
loop's attempts under the stage family where the transform makes them
(`_SteppingLoop` over `_StageSteps`: `_rk_step_coupled` taking
`_rk_step_transformed`) in C -- `csrc/loop.c`'s stage family in its
transform mode (`tmode`): where (the breakpoints, the run's end), take
(x_n's charge, the device memo's roll, x_n's C and G frozen; the two
factors by radau_fold.c -- on the core's pattern alone where the core
writes only that, x_n's C and G read in place -- with P's rows by numpy's
division, the real one's
LU in the loop's buffer -- or its values at the analysis solver's
StructuralKLU pattern, each real solve that solver's sequence -- the
complex one's values for the transient's ComplexKLUSolver; the seeds by pred.c; the sources at the stage times;
radau_tc.c's Newton and its end passes on the loop's own struct; the
end's record; each stage's screen on C by a core pass where the memo lacks
it; the sources at the end and the start; the estimate on the frozen real
solve), judge (lte.c's STAGE judgement, `_StageSteps.judge`'s law) and
accept (the rings, the stages and the point noted) -- bit for bit (pinned
in `tests/test_pinned_pairs.py`).

ENTRY: only where the transform's C Newton's last call served at the
watch counter as it stands (`_tran_radau_tc`'s `served`), its checks made;
then the C stage loop's checks of the run and the family
(`_tran_stageloop_c`'s, the transform's in place of its 'transform'), and
the transform's own: its pieces genuine (the frozen factors' chain, the C
Newton's, the step's), its switches, the analysis solver the dense one
(`NumpyLU`) or an AutoSolver whose choice is the circuit's StructuralKLU,
analysed (`_tran_radau_tc._real_klu`; not yet: a later entry's), the
complex solver analysed on the last record's pattern, the
transform's matrices made, the early abort's threshold one the C arms as
`_SteppingLoop.radau_arm` does, the device memo's generations within the
loop's model (`LOOP_TG_MAX` records each).  A decline is
counted (`tcloop_c:<why>`) and, where the reason is the run's, kept until
something watched moves; a stop at an attempt leaves it to the Python
before the next entry.

THE STATE.  Read before a call as `_tran_stageloop_c` reads it, and the
transform's: the device memo's generations by state and passes, x_n's
record where the memo holds it, `_radau_Pc`'s step.  Written back as the
Python attempts leave it: the loop's, the statistics', the transient's
rings and nodes as `_tran_stageloop_c` writes them; the last made
attempt's stage outputs (`_rk_Y`, `_q_cache`, `_iq`, `_Cmat`, `_Geq`,
`_rk_est`, the pending nodes of a rejected attempt, the screen's
`_C_cache`); the device memo's two generations (each record made again by
the core's passes at its state, or the entry's own arrays); `_radau_Pc`
for the last made attempt's step; the complex solver's refactors; the
real factor's StructuralKLU, its numeric and counts as its solves leave
them; the family's; the counts the Python's calls make.

A STOP before an attempt is made hands it back, nothing of it made (its
memo roll taken back); none holds a result: the Python's repeat refactors
the same complex values, deterministic, and the LU is the loop's own --
on KLU the real factor's numeric the one the attempt found, any it made
freed.

OWNS: on the `_SteppingLoop`'s C state (`_tran_loop_c.LoopState`) the
bound run (`_TcRun`).  `_paths`: `tcloop_c:` (the declines, the bails
`tcloop_c:bail:<why>`, `tcloop_c:attempts`, `tcloop_c:steps`).  Switch:
`PYCIRCUIT_TCLOOP_C=0` (`ENABLED`, read on every entry).

History: `doc/transient_history.md`, `_tran_tcloop_c.py`.
"""
import os
import time
from types import SimpleNamespace

import numpy as np

from pycircuit.circuit import _paths, _tran_loop_c, _tran_stageloop_c, _watch

#: `PYCIRCUIT_TCLOOP_C=0` makes the Python attempts every one
ENABLED = os.environ.get('PYCIRCUIT_TCLOOP_C', '1') != '0'
_PC = _paths.COUNTS
_no = _paths.no
_F64 = np.dtype(np.float64)
_LC = _tran_loop_c
_SL = _tran_stageloop_c
_CODES = _LC._CODES
BAIL = _LC.BAIL
#: a record's passes in the loop's memo model (loop.h's LOOP_TG_*), by name
_TG = {k: _CODES['LOOP_TG_' + k.upper()] for k in 'qiCG'}
_MOD = {}
#: the reasons kept until something watched moves: the run's
_KEEPS = {
    'driver': True, 'patched': True, 'stats': True, 'loop': True, 'shadow': True,
    'events': True, 'state': True, 'circuit': True, 'veto': True, 'family': True,
    'pcnr': True, 'est': True, 'units': True, 'sum': True, 'size': True,
    'sources': True, 'table': True, 'breakpoints': True, 'minbreak': True, 'screen': True,
    'judge': True, 'bind': True, 'method': True, 'early': True, 'solver': True,
    'switch': True,
    'passes': False, 'table_stale': False, 'toolkit': False, 'rings': False,
    'history': False, 'underflow': False, 'complex': False, 'memo': False, 'mats': False,
    'unchosen': False,
}

_M = 'pycircuit.circuit.'
#: the pieces the transform's attempts call beyond the stage loop's
#: (`_tran_stageloop_c._S_PIECES`) and the two chains read in place
#: (`_tran_radau._frozen_chain`, `_tran_radau_tc._chain`): (a getter of the
#: owners, its qualified name, its module)
_TC_PIECES = (
    (lambda o: o['T']._rk_step_transformed, '_RadauStages._rk_step_transformed', '_tran_radau'),
    (lambda o: o['T']._radau_frozen, '_RadauStages._radau_frozen', '_tran_radau'),
    (lambda o: o['T']._radau_early_check, '_RadauStages._radau_early_check', '_tran_radau'),
    (lambda o: o['T']._stage_end_passes, '_RadauStages._stage_end_passes', '_tran_radau'),
    (lambda o: o['T']._C_at_state, '_CompanionModel._C_at_state', '_tran_companion'),
    (lambda o: o['T']._C_lookup, '_CompanionModel._C_lookup', '_tran_companion'),
    (lambda o: o['R']._FrozenTransform.real_solve, '_FrozenTransform.real_solve', '_tran_radau'),
    (lambda o: o['R']._FrozenTransform.solve, '_FrozenTransform.solve', '_tran_radau'),
    (lambda o: o['R']._lam_folds, '_lam_folds', '_tran_radau'),
    (lambda o: o['R']._same_matrix, '_same_matrix', '_tran_radau'),
    (lambda o: o['L'].NumpyLU.make.__func__, 'NumpyLU.make', 'linearsolver'),
    (lambda o: o['L']._csc_of_dense, '_csc_of_dense', 'linearsolver'),
    (lambda o: o['L']._csc_of_values, '_csc_of_values', 'linearsolver'),
    (lambda o: o['TC'].solve, 'solve', '_tran_radau_tc'),
    (lambda o: o['TC']._outcome, '_outcome', '_tran_radau_tc'),
    (lambda o: o['TC'].fold, 'fold', '_tran_radau_tc'),
    (lambda o: o['TC'].resume, 'resume', '_tran_radau_tc'),
    (lambda o: o['TC'].iterate, 'iterate', '_tran_radau_tc'),
    (lambda o: o['TC'].seed, 'seed', '_tran_radau_tc'),
)
#: ... and on the transient, these shadowed send its attempts elsewhere
_TC_SHADOWS = ('_rk_step_transformed', '_radau_frozen', '_radau_early_check',
               '_stage_end_passes', '_C_at_state', '_C_lookup', '_radau_transform_solve',
               '_radau_complex_solve', '_radau_transform_matrices', '_get_linearsolver',
               '_transform_loop')


def _owners():
    """The objects the pieces are read from, imported once."""
    if 'owners' not in _MOD:
        from pycircuit.circuit import _tran_radau, _tran_radau_tc, linearsolver
        from pycircuit.circuit.transient import Transient
        _MOD['owners'] = {'T': Transient, 'R': _tran_radau, 'L': linearsolver,
                          'TC': _tran_radau_tc}
    return _MOD['owners']


def _now():
    own = _owners()
    return tuple(get(own) for get, _q, _m in _TC_PIECES)


def _pieces():
    """The pieces as their modules define them -- read once every one is
    (`_paths.genuine`); None: one is a stand-in, read again next entry."""
    g = _MOD.get('g')
    if g is None:
        objs = _now()
        if not all(_paths.genuine(o, q, _M + m) for o, (_g, q, m) in zip(objs, _TC_PIECES)):
            return None
        g = _MOD['g'] = objs
    return g


def _declined(loop, why):
    """The entry declined, counted -- kept (the mode's own, `tkept`) until
    something watched moves where the reason is the run's."""
    if _KEEPS[why]:
        loop.c.tkept = _watch.now()
    _no('tcloop_c:' + why)
    return False


def serving(tr, now):
    """Whether the transform's C Newton's last call on `tr` served at the
    watch counter `now` -- its checks made: the entry this module's
    (`_tran_loop_c.run`)."""
    rec = tr.__dict__.get('_radau_tc')
    return rec is not None and now != -1 and rec.served == now


def run(loop):
    """Make as many of `loop`'s attempts in C as it can: True where it made
    any (the state written back as the Python attempts leave it), else False
    -- the Python loop makes the next attempt (`_tran_loop_c.run`'s entry
    for a stage family whose transform served)."""
    if not ENABLED:
        return False
    c = loop.c
    if c.skip:
        c.skip -= 1
        return False
    now = _watch.now()
    if now == -1 or c.tkept == now:
        return False
    tr = loop.tr
    td = tr.__dict__
    rec = td.get('_radau_tc')
    if rec is None or rec.served != now:
        return False
    e = _eligible(loop, tr, td, rec)
    if e.__class__ is str:
        return _declined(loop, e)
    run_ = _bound(loop, e)
    if run_ is None:
        return _declined(loop, 'bind')
    _SL._load(loop, e, run_, stage=_load_tc)
    ## (settled whatever happens after the call: the C may have freed the
    ## solver's handles for a new pattern's, moved the stage KLU's numeric)
    status = None
    try:
        status = _call(e, run_)
    finally:
        _settle(e, run_)
        if e.fb is not None:
            _fb_settle(e, run_, status)
    s = run_.s
    if status not in (_CODES['LOOP_DONE'], _CODES['LOOP_FULL']):
        _no('tcloop_c:bail:' + BAIL.get(status, str(status)))
    if not (s.attempts or s.early_aborts):
        c.backoff = min(2 * c.backoff + 1, _LC.BACKOFF_MAX)
        c.skip = c.backoff
        return False
    c.backoff = 0
    ## (a stop at an attempt leaves that attempt to the Python before the
    ## next entry: an entry now would make it again)
    c.skip = 0 if status in (_CODES['LOOP_DONE'], _CODES['LOOP_FULL']) else 1
    _SL._write_back(loop, e, run_, stage=_stage_back, count=_count)
    return True


def _eligible(loop, tr, td, rec):
    """What an entry runs with (a namespace) -- or why not: the run
    (`_tran_stageloop_c._run_ready`), the family (`_family_ready`), the
    transform's state (`_tc_ready`), the sources (`_tran_loop_c._sources`),
    the size (`_fits`), the state (`_tran_stageloop_c._state_ready`), the
    screen (`_tran_loop_c._screen`), the judge's settings, the memo
    (`_memo_ready`)."""
    drv = _LC.driver()
    if drv is None:
        return 'driver'
    O = _LC._OWN or _LC._own()
    why = (_SL._run_ready(loop, tr, td, O) or _family_ready(loop, tr, td, rec)
           or _tc_ready(tr, td))
    if why:
        return why
    got = _LC._sources(loop, tr, O)
    if got.__class__ is str:
        return got
    tab, minbreak = got
    n = rec.n
    if not _fits(n):
        return 'size'
    st = _SL._state_ready(loop, tr, td, n)
    if st.__class__ is str:
        return st
    screen = _LC._screen(tr, td, rec, n)
    if screen is None:
        return 'screen'
    cfg = _SL._judge_ready(loop, tr, n) if not loop.fixed else None
    if cfg.__class__ is str:
        return cfg
    mem = _memo_ready(td, n)
    if mem is None:
        return 'memo'
    rk = _real_solver(tr, n)
    if rk.__class__ is str:
        return rk
    ## (the fallback's dense estimate is numpy's solve; where the real factor
    ## is on KLU the Python's is the analysis solver's: not bound then)
    return SimpleNamespace(drv=drv, O=O, tr=tr, rec=rec, n=n, iref=rec.iref, tab=tab,
                           minbreak=minbreak, st=st, screen=None if screen[0] is None else screen,
                           cfg=cfg, sk=None, fam=loop.family, integ=tr.base_integrator, mem=mem,
                           zs=td['_radau_zsolver'], mats=tr._radau_transform_matrices(),
                           rk=rk, fb=None if rk is not None else _fb_ready(tr, td, rec, n))


def _fb_ready(tr, td, rec, n):
    """The fallback's C Newton as `_tran_radau_c.solve` would serve it now --
    `(the transient's record, its temperature, the stage KLU)` -- or None,
    the loop then stopping at a fallback (the Python's): the record the
    Python's first fallback made, on this core, `_ready`'s checks passing on
    it (nothing remade), its slots and the loop's set within the budget."""
    from pycircuit.circuit import _tran_core, _tran_radau_c
    rc = td.get('_radau_c')
    ## (`_ready` remakes a record whose core, size or reference row moved:
    ## here it never does)
    if (rc is None or rc.core is not rec.core or rc.n != n or rc.iref != rec.iref
            or int(rc.s.maxa) != _tran_radau_c.MAXA or not _SL._sets_fit(n)
            or _tran_core.core_for(tr) is not rc.core):
        return None
    ctx = SimpleNamespace(iref=rec.iref, m=n - 1, h=1.0, qn=np.zeros(n),
                          Amat=np.array(tr.base_integrator.A, dtype=float))
    nobypass = float(getattr(tr.epar, 'bypasstol', -1.0) or -1.0) < 0.0
    got = _tran_radau_c._ready(tr, ctx, None, None, td.get('_stateful_lims'), nobypass,
                               tr.par.reltol, step_memo=False)
    if got.__class__ is str or got[0] is not rc:
        return None
    return got[0], got[1], got[2]


def _fits(n):
    """The loop's buffers for the transform within the C stage Newton's
    budget for its output sets (`_tran_radau_c.OUT_BYTES`): ten n x n and
    the fold's three m x m, complex two of them."""
    from pycircuit.circuit import _tran_radau_c
    return (10 * n * n + 5 * (n - 1) * (n - 1)) * 8 <= _tran_radau_c.OUT_BYTES


def _family_ready(loop, tr, td, rec):
    """Why the attempts are not the transform's the C makes (a reason), or
    '': the family and integrator, the pieces shadowed on the transient,
    numpy's underflow, PCNR, the transform on, the estimate's flag, the
    stage loop's pieces and the transform's, the units, the core's passes,
    numpy's sum."""
    from pycircuit.circuit import (
        _tran_core,
        _tran_lte_c,
        _tran_pred_c,
        _tran_radau_tc,
        func,
    )
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    fam = loop.family
    if (type(fam) is not _SL._owners()['F'] or type(tr.base_integrator) is not RadauIIA3Integrator
            or any(k in fam.__dict__ for k in _SL._F_NAMES)):
        return 'family'
    if any(k in td for k in _SL._T_NAMES) or any(k in td for k in _TC_SHADOWS):
        return 'shadow'
    if not func._underflow_ignored():
        return 'underflow'
    if tr._rk_use_pcnr():
        return 'pcnr'
    if not tr._newton_option(tr.par.radau_transform, 'radau_transform'):
        return 'method'
    if bool(getattr(tr, '_rk_want_est', False)) != (not loop.fixed):
        return 'est'
    g, h = _SL._MOD.get('g') or _SL._pieces(), _MOD.get('g') or _pieces()
    if g is None or h is None or _SL._now() != g or _now() != h:
        return 'patched'
    if not (_tran_radau_tc.ENABLED and _tran_pred_c.ENABLED and _tran_pred_c.driver()
            and _tran_lte_c.driver() and _tran_radau_tc.driver()
            and _tran_radau_tc._fold_driver()):
        return 'units'
    if not _tran_core.passes_held(tr, rec.core):
        return 'passes'
    if not _SL._sum_ok():
        return 'sum'
    return ''


def _tc_ready(tr, td):
    """Why the transform's state is not the one the C makes its attempts
    from (a reason), or '': its switches (the frozen factors, their fold,
    the end passes, the fused readers); the early abort's threshold a number
    the C arms where `radau_arm` does (it arms on any true one, the C on one
    above zero: not a negative one, not NaN); the frozen factors' chain and
    the C Newton's as `_radau_frozen` and `_tran_radau_tc.solve` read them;
    the analysis solver the dense one; the complex solver analysed and
    factored on the last record's pattern; the transform's matrices
    made."""
    from pycircuit.circuit import _tran_radau, _tran_radau_tc
    if not (_tran_radau.RADAU_FROZEN and _tran_radau.STAGE_FUSE and _tran_radau_tc.END
            and _tran_radau_tc.FOLD):
        return 'switch'
    th = tr.par.radau_early_abort
    if type(th) not in (float, int) or not th >= 0:
        return 'early'
    T = type(tr)
    fc = _tran_radau._FROZEN.get('chain') or _tran_radau._frozen_chain()
    tc = _tran_radau_tc._MOD.get('chain') or _tran_radau_tc._chain()
    if (fc is None or tc is None or tc[0] is not T or _tran_radau_tc._now(T, tr.cir) != tc[1:]
            or (T._radau_transform_solve, T._radau_complex_solve, T._radau_transform_matrices,
                T._get_linearsolver) != fc[:4]):
        return 'patched'
    why = _solver_ready(tr, fc)
    if why:
        return why
    zs = td.get('_radau_zsolver')
    ZS = _tran_radau._FROZEN['types'][0]
    last = getattr(zs, '_csc_last', None) if type(zs) is ZS else None
    if (last is None or not (zs._symbolic and zs._numeric) or zs._pattern != last[3]
            or (ZS.solve, ZS.prepare, ZS.solve_prepared, ZS.prepare_values) != fc[8:12]
            or any(k in zs.__dict__ for k in ('solve', 'prepare', 'solve_prepared',
                                              'prepare_values'))):
        return 'complex'
    mats = td.get('_radau_Tmats')
    if mats is None or not _tran_radau._lam_folds(mats[0]):
        return 'mats'
    return ''


def _solver_ready(tr, fc):
    """Why the real factor's solve is not numpy's kept LU as
    `_radau_frozen` makes it (`solver`), or '': the analysis solver the dense
    one (an AutoSolver that chose it, or a DenseSolver), each piece as
    defined, on the numeric toolkit, numpy's LAPACK found -- an AutoSolver
    that chose another, '' here (`_real_solver` reads it)."""
    from pycircuit.circuit import _tran_radau, linearsolver
    _ZS, AS, DS, NTK, NUM = _tran_radau._FROZEN['types']
    ls = tr._get_linearsolver()
    if type(ls) is AS:
        dense = ls._choice
        if not (AS.solve is fc[4] and AS._select is fc[5] and 'solve' not in ls.__dict__
                and '_select' not in ls.__dict__):
            return 'solver'
        if type(dense) is not DS:
            return ''
    else:
        dense = ls
    if not (type(dense) is DS and DS.solve is fc[6] and 'solve' not in dense.__dict__
            and type(tr.toolkit) is NTK and tr.toolkit.linearsolver is fc[7]
            and NUM.linearsolver is fc[7] and np.linalg.solve is _tran_radau._NP_SOLVE
            and linearsolver._numpy_lapack() is not None):
        return 'solver'
    return ''


def _real_solver(tr, n):
    """The real factor's solver the C's real solves stand in for: None --
    numpy's kept LU (`_solver_ready`'s) -- or the analysis solver's
    StructuralKLU (`_tran_radau_tc._real_klu`: each piece as defined, of this
    size), or why not: `unchosen` (the AutoSolver's choice not made, or its
    StructuralKLU not analysed: a later entry's), `solver`."""
    from pycircuit.circuit import _tran_radau, _tran_radau_tc
    from pycircuit.circuit.linearsolver import StructuralKLU
    AS, DS = _tran_radau._FROZEN['types'][1:3]
    ls = tr._get_linearsolver()
    if type(ls) is not AS or type(ls._choice) is DS:
        return None
    sk = ls._choice
    if sk is None or (type(sk) is StructuralKLU and not sk._symbolic):
        return 'unchosen'
    sk = _tran_radau_tc._real_klu(ls, n - 1)
    return 'solver' if sk.__class__ is str else sk


def _memo_ready(td, n):
    """The device memo's two generations by state and passes -- `((states,
    passes), (states, passes))`, newest first, each at most LOOP_TG_MAX
    records -- or None: more records, or a record not of the run's width."""
    memo = td.get('_dev_memo') or ({}, {})
    mx = _CODES['LOOP_TG_MAX']
    out = []
    for gen in memo:
        if len(gen) > mx:
            return None
        states, passes = [], []
        for key, rec in gen.items():
            if type(key) is not bytes or len(key) != 8 * n or not rec.keys() <= _TG.keys():
                return None
            states.append(key)
            passes.append(sum(_TG[k] for k in rec if k in _TG))
        out.append((states, passes))
    return tuple(out)


class _TcRun(_SL._StageRun):
    """The transform's bound run: the stage loop's (`_tran_stageloop_c.
    _StageRun`) and its own -- the loop's radau_tf_t context (`tctx`), its
    fold struct (`fo`), the transform's buffers (`tbuf`, kept with their
    handles in `tkeep`)."""

    __slots__ = ('ccounts', 'fo', 'tbuf', 'tctx', 'tkeep')


def _bound(loop, e):
    """The loop's bound transform run, bound again where the record, the
    table, the screen's presence or certificate, the complex solver's
    handles or pattern, or the sizes moved; None where it cannot be."""
    run_ = loop.c.run
    scr, zs = e.screen, e.zs
    fb = e.fb or (None, None, None)
    key = (e.rec, e.tab, zs, zs._symbolic, zs._numeric, zs._csc_last[1],
           None if scr is None else scr[3], None if scr is None else scr[4], e.n,
           e.st[0].shape[0], fb[0], fb[2], e.rk)
    if not (type(run_) is _TcRun and len(run_.key) == len(key)
            and all(a is b or (type(a) is int and a == b) for a, b in zip(run_.key, key))):
        run_ = _bind(e)
        if run_ is not None:
            run_.key = key
            loop.c.run = run_
    return run_


def _bind(e):
    """The transform run for this entry: the stage loop's units
    (`_tran_stageloop_c.bind_common`), the C Newton on a radau_tf_t of the
    loop's own (`_tran_radau_tc._Ctx`, not the transient's), the fold on a
    radau_fold_t of its own, the transform's buffers."""
    from pycircuit.circuit import _tran_core, _tran_radau_tc
    run = _SL.bind_common(e, _TcRun())
    s, ffi, n, rec = run.s, run.ffi, e.n, e.rec
    tdrv = _tran_radau_tc.driver()
    cffi_core, cfn_core, _dg = _tran_core.driver()
    tctx = run.tctx = _tran_radau_tc._Ctx(rec.core, n, e.iref, tdrv[0], cfn_core, cffi_core, tdrv)
    _LC._plug(s, ffi, 'tc', tdrv[0], tctx.s, tdrv[1])
    _bind_newton(run, e)
    fdrv = _tran_radau_tc._fold_driver()
    run.fo = fdrv[0].new('radau_fold_t *')
    _LC._plug(s, ffi, 'fold', fdrv[0], run.fo, fdrv[1])
    _bind_buffers(run, e)
    s.fb_ok = 0
    if e.fb is not None:
        _bind_fb(run, e)
    s.rk_klu = run.tctx.s.rk_on = 0
    if e.rk is not None:
        _bind_rk(run, e)
    run.ccounts = _cir_c_counts(e)
    s.tmode = 1
    return run


def _bind_rk(run, e):
    """The real factor on KLU (radau_tc.h's `rk`, the loop's to hold:
    `rk_hold`): the analysis solver's pattern, library entries, common and
    analysis; its values' buffer; the residual's dense work; the pattern's
    flat indices for `tc_frozen`'s take; the solves' work (its numeric and
    tolerance each call's, `_load`)."""
    from pycircuit.circuit import _cbuild
    sk, tctx = e.rk, run.tctx
    tffi, ts, s, ffi = tctx.ffi, tctx.s, run.s, run.ffi
    tb, keep = run.tbuf, run.tkeep
    tb['rk_take'] = np.ascontiguousarray(sk.take, dtype=np.int64)
    tb['rk_ax'] = np.empty(max(int(sk.nnz), 1))
    tb['rkx'] = np.empty(max(e.n - 1, 1))
    k = ts.rk
    for name, arr, kind in (('kAp', sk.Ap, 'int *'), ('kAi', sk.Ai, 'int *'),
                            ('kAx', tb['rk_ax'], 'double *')):
        keep.append(tffi.from_buffer(kind, arr))
        setattr(k, name, keep[-1])
    tb['rk_A'] = np.empty(max((e.n - 1) ** 2, 1))
    keep.append(tffi.from_buffer('double *', tb['rk_A']))
    ts.rk_A = keep[-1]
    lib = sk._lib
    k.klu_factor, k.klu_refactor, k.klu_solve, k.klu_free = (
        tffi.cast('void *', _cbuild.address(f))
        for f in (lib.klu_factor, lib.klu_refactor, lib.klu_solve, lib.klu_free_numeric))
    k.klu_common = tffi.cast('void *', sk._ct.addressof(sk._common))
    k.klu_sym = tffi.cast('void *', sk._symbolic)
    keep.append(ffi.from_buffer('int64_t *', tb['rk_take']))
    s.rk_take = keep[-1]
    keep.append(ffi.from_buffer('double *', tb['rkx']))
    s.rkx = keep[-1]
    s.rk_nnz = int(sk.nnz)
    ts.rk_on = ts.rk_hold = 1
    s.rk_klu = 1


def _bind_fb(run, e):
    """The fallback's units: radau.c on the transient's record (plugged as
    `_tran_stageloop_c._bind` plugs it), the coupled family's two output sets
    (the record's own and one of the loop's), the dense estimate's work and
    numpy's dgesv, the held charge, the last fallback's stages."""
    from pycircuit.circuit import _tran_radau_c
    rec = e.fb[0]
    s, ffi, n, m = run.s, run.ffi, e.n, e.n - 1
    rffi, rcfn, rdgesv, _rdgemv = _tran_radau_c.driver()
    _LC._plug(s, ffi, 'radau', rffi, rec.s, rcfn)
    run.sets = sets = [rec.own, rec._new_out()]
    tb, keep, fb = run.tbuf, run.tkeep, ffi.from_buffer
    for k, size in (('eA', max(m * m, 1)), ('hq', n), ('fbY', 3 * n)):
        tb[k] = np.empty(size)
        keep.append(fb('double *', tb[k]))
        setattr(s, k, keep[-1])
    tb['eipiv'] = np.empty(max(m, 1), dtype=np.int64)
    keep.append(fb('int64_t *', tb['eipiv']))
    s.eipiv = keep[-1]
    for o in range(2):
        for k, nm in enumerate(_tran_radau_c._OUTS):
            keep.append(fb('double *', sets[o][0][nm]))
            s.oset[o][k] = keep[-1]
    s.dgesv = ffi.cast('void *', rdgesv)
    s.fb_ok = 1


def _cir_c_counts(e):
    """The counts one `cir.C` call makes (the screen's `_C_at_state` where
    the memo lacks the stage): a second call's at x_n, the first's
    taken back with it (a first may build) -- nothing left counted."""
    tr, x = e.tr, e.st[4]
    out = {}
    for _ in range(2):
        before = _paths.snapshot()
        tr.cir.C(x, tr.epar)
        out = _paths.since(before)
        for k, v in out.items():
            _PC[k] -= v
    return out


def _bind_newton(run, e):
    """The loop's radau_tf_t's run constants (`_tran_radau_tc.solve`'s for a
    record): the tableau, V's entries as `_FrozenTransform.v`, the end
    passes' bits, the complex solver's library entries, handles and the
    last record's pattern; the product M = diag(lam) Tinv's rows 0 and 1,
    numpy's, from which each step's P is made."""
    import ctypes

    from pycircuit.circuit import _cbuild, _tran_radau_tc
    ts, tffi, zs = run.tctx.s, run.tctx.ffi, e.zs
    Amat = np.array(e.integ.A, dtype=float)
    for k, v in enumerate(Amat.ravel().tolist()):
        ts.A[k] = v
    lam, V, Tinv = e.mats
    for i in range(3):
        ts.v[3 * i] = float(V[i, 0].real)
        ts.v[3 * i + 1] = float(V[i, 1].real)
        ts.v[3 * i + 2] = float(V[i, 1].imag)
    M = np.diag(lam) @ Tinv
    s = run.s
    for k, z in enumerate((M[0, 0], M[0, 1], M[0, 2], M[1, 0], M[1, 1], M[1, 2])):
        s.tcm[2 * k], s.tcm[2 * k + 1] = float(z.real), float(z.imag)
    s.lam0r, s.lam1r, s.lam1i = float(lam[0].real), float(lam[1].real), float(lam[1].imag)
    ts.end_bits, ts.end_cg, ts.fused = 7, 0, int(_tran_radau_tc._MOD['cmul'])
    ts.klu_refactor = tffi.cast('void *', _cbuild.address(zs._lib.klu_z_refactor))
    ts.klu_solve = tffi.cast('void *', _cbuild.address(zs._lib.klu_z_solve))
    ts.Symbolic = tffi.cast('void *', zs._symbolic)
    ts.Numeric = tffi.cast('void *', zs._numeric)
    ts.Common = tffi.cast('void *', ctypes.addressof(zs._common))
    _nz, Ap, Ai, _key = zs._csc_last
    run.tctx.keep_ap = (Ap, Ai, tffi.from_buffer('int *', Ap), tffi.from_buffer('int *', Ai))
    ts.Ap, ts.Ai = run.tctx.keep_ap[2], run.tctx.keep_ap[3]
    ## (a new pattern's analysis, fresh factor and frees, `tc_repattern`)
    ffi = run.ffi
    for field, name in (('klu_analyze', 'klu_analyze'), ('klu_zfactor', 'klu_z_factor'),
                        ('klu_free_sym', 'klu_free_symbolic'),
                        ('klu_free_znum', 'klu_z_free_numeric')):
        setattr(s, field, ffi.cast('void *', _cbuild.address(getattr(zs._lib, name))))


def _bind_buffers(run, e):
    """The transform's buffers, each pointed at by its field and kept: the
    stage loop's names (the sources, seeds, x_n's passes, the estimate's
    work, the made attempt's `iq` and estimate), the reduced C and G (the
    fold's inputs too), the LU and its pivots, the stages' C, the end
    passes' two sets, the memo's three slots, the last made attempt's
    stages; the fold's outputs and the last record's pattern."""
    from pycircuit.circuit import _tran_radau_c
    ffi, s, n, m = run.ffi, run.s, e.n, e.n - 1
    nn, mm, mx = n * n, max(m * m, 1), _CODES['LOOP_TG_MAX']
    tb = run.tbuf = {
        'ust': np.empty(5 * n), 'useed': np.empty(3 * n), 'xr_q': np.empty(n),
        'xr_C': np.zeros(nn), 'xr_G': np.empty(nn), 'xr_i': np.empty(n), 'eF': np.empty(n),
        'eb': np.empty(max(m, 1)), 'liq': np.empty(n), 'lest': np.empty(n),
        'crr': np.empty(mm), 'grr': np.empty(mm), 'lua': np.empty(mm), 'tcC': np.zeros(2 * nn),
        'tcq': np.empty(n), 'tcY': np.empty(3 * n),
        'tg': np.empty(_CODES['LOOP_TG_SLOTS'] * mx * n),
        'tcf0': np.empty(n)}
    run.tkeep = keep = []
    fb = ffi.from_buffer
    for k, a in tb.items():
        keep.append(fb('double *', a))
        setattr(s, k, keep[-1])
    tb['lupiv'] = np.zeros(max(m, 1), dtype=np.int64)
    keep.append(fb('int64_t *', tb['lupiv']))
    s.lupiv = keep[-1]
    for o in range(2):
        for k, size in enumerate((n, n, nn, nn)):
            a = tb[f'te{o}{k}'] = np.zeros(size)
            keep.append(fb('double *', a))
            s.te[o][k] = keep[-1]
        ## (a new pattern's buffers: Ap, Ai, its mask laid out as A.T)
        tb[f'pAp{o}'] = np.zeros(m + 1, dtype=np.int32)
        tb[f'pAi{o}'] = np.zeros(max(m * m, 1), dtype=np.int32)
        tb[f'pnz{o}'] = np.zeros(max(m * m, 1), dtype=np.uint8)
        keep.append(fb('int *', tb[f'pAp{o}']))
        s.pAp[o] = keep[-1]
        keep.append(fb('int *', tb[f'pAi{o}']))
        s.pAi[o] = keep[-1]
        keep.append(fb('uint8_t *', tb[f'pnz{o}']))
        s.pnz[o] = keep[-1]
    s.dgemv = ffi.cast('void *', _tran_radau_c.driver()[3])
    _bind_fold(run, e, m)


def _bind_fold(run, e, m):
    """The loop's fold struct: its inputs the loop's reduced C and G, its
    outputs its own, the complex solver's last pattern (`_csc_last`'s
    mask, laid out as `_tran_radau_tc._fold_nz` reads it); on the core's
    pattern, the pattern's positions (`_bind_fold_pattern`)."""
    fo, tb = run.fo, run.tbuf
    fffi = _owners()['TC']._fold_driver()[0]
    fo.m, fo.fused = m, int(_owners()['TC']._MOD['cmul'])
    nz, _Ap, Ai, _key = e.zs._csc_last
    at = 1 if nz.flags.c_contiguous else 0
    tb['fold_rf'] = np.empty(max(m * m, 1))
    tb['fold_cf'] = np.empty(max(2 * m * m, 2))
    tb['fold_ax'] = np.empty(max(2 * m * m, 2))
    tb['fold_nz'] = nz if at else nz.T
    keep = run.tkeep
    for k, name in (('Cr', 'crr'), ('Gr', 'grr')):
        keep.append(fffi.cast('double *', _LC._addr(run.ffi, getattr(run.s, name))))
        setattr(fo, k, keep[-1])
    for k in ('rf', 'cf', 'ax'):
        keep.append(fffi.from_buffer('double *', tb['fold_' + k]))
        setattr(fo, k, keep[-1])
    ## (the complex values the Newton's KLU refactors: the fold's `ax`)
    tffi = run.tctx.ffi
    keep.append(tffi.from_buffer('double *', tb['fold_ax']))
    run.tctx.s.Ax = keep[-1]
    keep.append(fffi.from_buffer('uint8_t *', tb['fold_nz']))
    fo.nz, fo.nz_at, fo.nnz = keep[-1], at, int(Ai.shape[0])
    fo.onpat = 0
    if e.rec.core.npat:
        _bind_fold_pattern(run, e, m, fffi)


def _bind_fold_pattern(run, e, m, fffi):
    """The fold on the core's pattern (radau_fold.h's `onpat`): where the
    core writes only its pattern into a C and G zero off it, the factors
    are made at its positions alone -- reduced (the reference row and
    column out), column by column and rows ascending: their flat indices
    into x_n's n x n matrices, into the factors and into their transpose;
    the factors zeroed off it by the first call."""
    fo, tb, keep = run.fo, run.tbuf, run.tkeep
    n, r0 = e.n, e.iref
    i, j = np.divmod(np.asarray(e.rec.core.arrays['pat'], dtype=np.int64), n)
    on = (i != r0) & (j != r0)
    i, j = i[on], j[on]
    ri, rj = i - (i > r0), j - (j > r0)
    o = np.lexsort((ri, rj))
    for name, at in (('pfull', i * n + j), ('pred', ri * m + rj), ('ptr', rj * m + ri)):
        tb['fold_' + name] = np.ascontiguousarray(at[o], dtype=np.int64)
        keep.append(fffi.from_buffer('int64_t *', tb['fold_' + name]))
        setattr(fo, name, keep[-1])
    fo.pn, fo.pclean, fo.onpat = int(o.shape[0]), 0, 1


def _load_tc(e, run_, xl):
    """The transform's state as the Python leaves it: x_n's passes where the
    memo holds them (`xr_*`, `xfields`); the memo's two generations into
    the model's slots 0 (the last attempt's) and 1; the step `_radau_Pc`
    holds P for (where its matrices are the transient's); the Newton's
    per-run settings; the core's passes' bits; numpy's underflow."""
    from pycircuit.circuit import func
    s, tr, ts = run_.s, e.tr, run_.tctx.s
    _load_x(e, run_, xl)
    _load_memo(e, run_)
    pc, mats = tr.__dict__.get('_radau_Pc'), e.mats
    s.tc_hp = int(pc is not None and type(pc[0]) is float and pc[3] is mats[0]
                  and pc[4] is mats[1] and pc[5] is mats[2])
    s.tc_h = pc[0] if s.tc_hp else 0.0
    if s.tc_hp:
        ## (P's rows for that step, as the Newton reads them)
        for k, z in enumerate(pc[1]):
            ts.p[2 * k], ts.p[2 * k + 1] = float(z.real), float(z.imag)
    s.have_ph = 0
    s.pcur = -1
    _load_early(e, run_)
    ts.T = e.rec.core.probe(tr.epar)
    ts.reltol, ts.abstol = float(tr.par.reltol), float(tr.par.vabstol)
    ts.rtol = float(e.zs.REFACTOR_RESIDUAL_TOL)
    s.tc_maxit = int(tr.par.maxiter)
    rs = e.rec.s
    for k in ('walk_fn', 'wF', 'wPR', 'wTI', 'wNM', 'wOFF', 'wK', 'wn'):
        setattr(ts, k, getattr(rs, k))
    abi = _LC._ABI
    npat = bool(e.rec.core.npat)
    s.pass_bits = abi['PCX_WANT_ALL'] | (abi['PCX_WANT_PATTERN'] if npat else 0)
    s.pass_zero = int(npat)
    s.under_ok = int(func._underflow_ignored())
    run_.pr[0].under_ok = s.under_ok
    run_.lctx.s.under_ok = s.under_ok
    s.a33 = float(np.array(e.integ.A, dtype=float)[2, 2])
    integ = e.integ
    cvec = np.array(integ.C, dtype=float)
    for i in range(3):
        s.ct[i] = cvec[i]
    s.dd1, s.dd2, s.dd3 = _owners()['R']._DD
    s.fbx = 0
    if e.fb is not None:
        _load_fb(e, run_)
    if e.rk is not None:
        ## (the analysis solver's numeric now: the attempts' sequences start
        ## from it, `tc_take`)
        tffi, k = run_.tctx.ffi, run_.tctx.s.rk
        num = e.rk._numeric
        k.klu_num = k.klu_out = tffi.cast('void *', num) if num else tffi.NULL
        k.resid_tol = float(e.rk.REFACTOR_RESIDUAL_TOL)


def _load_fb(e, run_):
    """The fallback's: its record's run settings as `_tran_radau_c.solve`
    sets them (the temperature, the tolerances, maxiter, the tableau), no
    coupled set made yet, the stage KLU's numeric, the estimate's gamma."""
    rec, Tc, sk = e.fb
    s, ffi, rs, tr = run_.s, run_.ffi, rec.s, e.tr
    rs.T = Tc
    rs.reltol, rs.abstol = float(tr.par.reltol), float(tr.par.vabstol)
    rs.maxiter = int(tr.par.maxiter)
    for k, v in enumerate(np.array(e.integ.A, dtype=float).ravel().tolist()):
        rs.A[k] = v
    s.gen0, s.A0 = -1, 0
    num = sk._numeric if sk is not None else None
    s.klu_cur = ffi.cast('void *', num) if num else ffi.NULL
    s.gamma_r = e.integ.GAMMA_REAL


def _load_early(e, run_):
    """The early abort's threshold (0 where the loop never arms it: a fixed
    grid) and the family's `early_f0` where it is this point's (`pt[0] is
    X[-1]`)."""
    s, tr = run_.s, e.tr
    s.tc_th = 0.0 if s.adaptive == 0 else float(tr.par.radau_early_abort)
    pt = getattr(e.fam, 'early_f0', None)
    X = e.st[3]
    s.have_f0 = int(pt is not None and pt[0] is X[-1])
    if s.have_f0:
        run_.tbuf['tcf0'][:] = pt[1]
        s.f0_point = len(X) - 1


def _load_x(e, run_, xl):
    """x_n's passes where the memo holds them -- either generation's: the
    values one state gives -- into `xr_*`, their passes in `xfields`;
    `te_x` -1 (the end sets hold no record yet)."""
    s, tb = run_.s, run_.tbuf
    memo = e.tr.__dict__.get('_dev_memo') or ({}, {})
    key = xl.tobytes()
    f = 0
    for gen in memo:
        rec = gen.get(key)
        if rec is None:
            continue
        for k, name in (('q', 'xr_q'), ('i', 'xr_i'), ('C', 'xr_C'), ('G', 'xr_G')):
            if k in rec and not f & _TG[k]:
                tb[name][:] = np.asarray(rec[k], dtype=float).reshape(-1)
                f |= _TG[k]
    s.xfields, s.te_x = f, -1


def _load_memo(e, run_):
    """The memo's two generations into the model's slots: 0 the last
    attempt's (`_dev_memo[0]`), 1 the one before; slot 2 the next
    attempt's."""
    s, n = run_.s, e.n
    tg = run_.tbuf['tg']
    mx = _CODES['LOOP_TG_MAX']
    for g, (states, passes) in enumerate(e.mem):
        for k, (key, f) in enumerate(zip(states, passes)):
            off = (g * mx + k) * n
            tg[off:off + n] = np.frombuffer(key, dtype=np.float64)
            s.tg_f[g][k] = f
        s.tg_n[g] = len(states)
    s.tg_n[2] = s.tg_n[3] = 0
    s.tg_gen, s.tg_prev = 0, 1


def _call(e, run_):
    """The C's attempts: its status, the time the solve's.  The fallback's
    record's inputs and outputs, pointed at the loop's buffers while the C
    runs, pointed back at its own after (`_tran_stageloop_c._call`): the
    Python's next `_tran_radau_c.solve` reads and writes those."""
    rs = None if e.fb is None else e.fb[0].s
    saved = None if rs is None else (rs.qn, rs.u, rs.seed, rs.Ya, rs.qa, rs.ia, rs.Ca, rs.Ga)
    t0 = time.perf_counter()
    try:
        return e.drv[1](run_.s)
    finally:
        if saved is not None:
            rs.qn, rs.u, rs.seed, rs.Ya, rs.qa, rs.ia, rs.Ca, rs.Ga = saved
        e.tr.statistics.solve_seconds += time.perf_counter() - t0


def _settle(e, run_):
    """The complex solver as the made attempts leave it: their refactors
    (each attempt's record refactored once; an attempt not made refactored
    values the Python's repeat refactors again), and a new pattern's
    analyses and fresh factors -- its record (`_csc_last`, `_pattern`, its
    `Ap`, `Ai` and mask as `_csc_of_dense` makes them) and its handles, the
    analysis and factor it replaced freed by the C."""
    s, zs = run_.s, e.zs
    if e.rk is not None:
        ## (the analysis solver as its own solves leave it: the numeric the
        ## made attempts' sequences end on, their counts)
        sk, kr = e.rk, run_.tctx.s.rk
        sk._numeric = int(run_.tctx.ffi.cast('uintptr_t', kr.klu_out)) or None
        sk.refactors += int(s.tc_rkref)
        sk.factors += int(s.tc_rkfac)
        sk.pivot_fallbacks += int(s.tc_rkpiv)
        sk.residual_fallbacks += int(s.tc_rkres)
    k = int(s.tc_refactors)
    if k:
        zs.refactors += k
    if s.tc_analyses:
        zs.analyses += int(s.tc_analyses)
        zs.factors += int(s.tc_factors)
        zs.analyses_freed += int(s.tc_symfreed)
        zs.factors_freed += int(s.tc_numfreed)
    if s.pcur < 0:
        return
    tb, ts, m = run_.tbuf, run_.tctx.s, e.n - 1
    o = int(s.pcur)
    Ap = tb[f'pAp{o}'][:m + 1].copy()
    Ai = tb[f'pAi{o}'][:int(Ap[m])].copy()
    nz = tb[f'pnz{o}'][:m * m].reshape(m, m).astype(bool)
    key = (m, Ap.tobytes(), Ai.tobytes())
    tffi = run_.tctx.ffi
    zs._symbolic = int(tffi.cast('uintptr_t', ts.Symbolic))
    zs._numeric = int(tffi.cast('uintptr_t', ts.Numeric))
    zs._pattern = key
    zs._csc_last = (nz, Ap, Ai, key)


def _fb_settle(e, run_, status):
    """The fallback's after the call: the stage KLU as its made Newtons left
    it (its numeric, their counts), and where the loop stopped after a
    fallback's converged Newton, its result held for the Python's repeat
    (`_tran_radau_c.hold`, as `_tran_stageloop_c._settle` holds it)."""
    from pycircuit.circuit import _tran_radau_c
    rec, _Tc, sk = e.fb
    s, n = run_.s, e.n
    if (status == _CODES['LOOP_BAIL_RADAU']
            and rec.s.status == _tran_radau_c._CODES['RADAU_BAIL_WALKSTOP']):
        rec.walk_ok = None
    if sk is not None:
        sk._numeric = int(run_.ffi.cast('uintptr_t', s.klu_cur)) or None
        sk.refactors += s.st_klu_refactors
        sk.factors += s.st_klu_factors
        sk.pivot_fallbacks += s.st_klu_pivot
        sk.residual_fallbacks += s.st_klu_resid
    if s.held and status in _SL._AFTER:
        tb, rs = run_.tbuf, rec.s
        o = 1 if s.gen0 == 0 else 0
        _tran_radau_c.hold(rec, sk, tb['hq'], tb['ust'][:3 * n].reshape(3, n),
                           tb['useed'].reshape(3, n), np.array(e.integ.A, dtype=float),
                           float(rs.T), float(s.cur_h), float(rs.reltol), float(rs.abstol),
                           int(rs.maxiter), outs=run_.sets[o][0])


def _stage_back(e, run_):
    """The last made attempt's stage outputs as `_finish_stage_step` and the
    estimate leave them, the device memo's two generations, `_radau_Pc`,
    the screen's `_C_cache`: its point."""
    s, tb, tr, n = run_.s, run_.tbuf, e.tr, e.n
    gens, have = _memo_back(e, run_)
    Y = [tb['tcY'][j * n:(j + 1) * n].copy() for j in range(3)]
    ## (its end record: the memo's -- or, rolled out by attempts abandoned
    ## since at the early abort, the same passes made again)
    key3 = Y[2].tobytes()
    rec3 = gens[0].get(key3) or gens[1].get(key3)
    if rec3 is None or len(rec3) < 4:
        rec3 = _record(e, Y[2], key3, sum(_TG.values()), have)
    h = float(s.att_h)
    tr._rk_Y = Y
    tr._q_cache = (Y[2], rec3['q'])
    tr._iq = tb['liq'].copy()
    tr._Cmat = rec3['C']
    tr._Geq = float(s.a33) * h * rec3['G']
    tr._effective_method = type(e.integ).__name__
    tr._companion_coeffs = None
    tr._dt = h
    if s.adaptive:
        tr._rk_est = tb['lest'].copy()
    ts = [np.float64(s.att_ts[i]) for i in range(3)]
    tr._pred_pending = ((float(s.att_tend), list(zip(ts, Y))) if s.att_rejected else None)
    if run_.br is not None:
        tr._C_cache = (Y[2], rec3['C'])
    if s.have_ph:
        _pc_back(e, float(s.ph_made))
    if s.tc_fb:
        _fb_back(e, run_, gens, have)
    return Y[2]


def _fb_back(e, run_, gens, have):
    """What the made fallbacks leave on the transient: their count
    (`_radau_transform_fallbacks`), and the last one's mark of its stages' C
    with the core's pattern (`_C_pattern`, its memo records' C where the memo
    still holds them, else made again)."""
    s, n, tr = run_.s, e.n, e.tr
    tr._radau_transform_fallbacks = getattr(tr, '_radau_transform_fallbacks', 0) + int(s.tc_fb)
    core = e.fb[0].core
    if not (s.have_fbY and core.npat):
        return
    recs = []
    for j in range(3):
        y = run_.tbuf['fbY'][j * n:(j + 1) * n].copy()
        key = y.tobytes()
        r = gens[0].get(key) or gens[1].get(key)
        recs.append(r if r is not None and 'C' in r else _record(e, y, key, _TG['C'], have))
    tr._C_pattern = (tuple(r['C'] for r in recs), core.arrays['pat'])


def _pc_back(e, h):
    """`_radau_Pc` as the last made attempt's `_radau_frozen` leaves it: P's
    entries for its step, by numpy as there."""
    lam, V, Tinv = e.mats
    P = (np.diag(lam) @ Tinv) / h
    e.tr._radau_Pc = (h, (P[0, 0], P[0, 1], P[0, 2], P[1, 0], P[1, 1], P[1, 2]),
                      tuple((V[i, 0].real, V[i, 1]) for i in range(3)), lam, V, Tinv)


def _memo_back(e, run_):
    """The device memo the Python attempts leave, from the model: each
    record of the last attempt's generation and the one before, its passes
    the entry's own arrays where the entry held them at that state, else
    made by one core pass there (the values each pass makes) -- and those
    arrays by state (`have`)."""
    s, tb, n = run_.s, run_.tbuf, e.n
    mx = _CODES['LOOP_TG_MAX']
    entry = e.tr.__dict__.get('_dev_memo') or ({}, {})
    have = {}
    for gen in (entry[1], entry[0]):
        for key, rec in gen.items():
            have.setdefault(key, {}).update(rec)
    out = []
    for g in (int(s.tg_gen), int(s.tg_prev)):
        d = {}
        for k in range(int(s.tg_n[g])):
            off = (g * mx + k) * n
            x = tb['tg'][off:off + n].copy()
            key = x.tobytes()
            d[key] = _record(e, x, key, int(s.tg_f[g][k]), have)
        out.append(d)
    e.tr._dev_memo = (out[0], out[1])
    e.tr._memo_rolling = True
    return out, have


def _record(e, x, key, f, have):
    """The record of the passes `f` at `x`: the entry's arrays where it held
    them, the rest by one core call (formula 4, every output, C its own; G
    and i from the core's bins), kept in `have` for the state's other
    records -- nothing counted."""
    from pycircuit.circuit import _tran_core
    got = have.setdefault(key, {})
    if any(_TG[k] & f and k not in got for k in _TG):
        ffi, cfn, _dgemv = _tran_core.driver()
        core, n = e.rec.core, e.n
        fb, dp = ffi.from_buffer, core.dptr
        dummy = np.zeros(max(n * n, 1))
        d = fb(dp, dummy)
        q, C = np.empty(n), np.empty((n, n))
        cfn(core.cs, fb(dp, x), float(core.probe(e.tr.epar)), _LC._ABI['PCX_WANT_ALL'], 4, 0.0,
            0.0, 0.0,
            1.0, 0.0, d, d, d, d, d, fb(dp, C), fb(dp, q), d, d, d, d)
        made = {'q': q, 'i': core.arrays['bini'].copy(), 'C': C,
                'G': core.arrays['binG'].reshape(n, n).copy()}
        for k in _TG:
            got.setdefault(k, made[k])
    return {k: got[k] for k in ('q', 'i', 'C', 'G') if _TG[k] & f}


def _count(e, run_, k):
    """The `_paths` counts the Python transform attempts' calls make, made
    for the C's -- every attempt its frozen factors reached, made or
    abandoned at the early abort: the C Newton's (serves, iterations,
    probes, end passes, the source memo's hits over its iterations), the
    frozen factors and their fold, the readers' fused passes, the early
    abort's checks and aborts, the core's passes at x_n and the end's, the
    screens, the sources' passes, the seeds, the estimates, the accepted
    points' charge memo."""
    from pycircuit.circuit import _tran_core
    s = run_.s
    _PC['tcloop_c:attempts'] += int(s.tc_frozen)
    _PC['tcloop_c:steps'] += k
    if s.tc_rkdense:
        ## (the real factor on KLU: the residual tests its sparse verdict
        ## left to the dense one)
        _PC['tcloop_c:klu_dense'] += int(s.tc_rkdense)
    for key, v in (('radau_tc:served', s.tc_served), ('radau_tc:iterations', s.tc_iters),
                   ('radau_tc:probe', s.tc_probes), ('radau_tc:end', s.tc_ends - s.tc_pyends),
                   ('radau_tc:klu', s.tc_klu),
                   ('radau.frozen:served', s.tc_frozen), ('radau.frozen:fold', s.tc_folds),
                   ('radau.frozen:fold.same', s.tc_folds - s.tc_foldoff),
                   ('radau.fuse:served', s.tc_ends),
                   ('radau.fuse:handed', s.tc_handed), ('radau.fuse:unhanded', s.tc_unhanded),
                   ('radau.early:checks', s.tc_checks), ('radau.early:aborts', s.tc_aborts),
                   ('radau.est:transform', s.tc_ests), ('radau.est:memo', s.st_estmemo),
                   ('pred:fast', s.st_fast), ('radau_tc:noconv', s.tc_noconv)):
        if v:
            _PC[key] += int(v)
    _PC['umemo:hit'] += int(s.tc_uhit + s.st_uhit)
    _PC['umemo:miss'] += int(s.st_umiss)
    _PC['src.direct:served'] += int(s.st_umiss)
    _PC['u:called'] += len(e.tab.entries) * int(s.st_umiss)
    _count_screen(e, run_)
    ## (the core's passes: at x_n, the end's, and the Python loop's -- each
    ## stage's 'qi' an iteration -- with its limiter walk where the
    ## circuit limits)
    for _ in range(int(s.tc_xpasses + s.tc_endpass + s.tc_unhanded + 3 * s.tc_pyiters)):
        _tran_core.passes_count(e.tr, e.rec.core)
    if s.tc_pyiters and int(run_.tctx.s.wn):
        _PC['walk:served'] += 3 * int(s.tc_pyiters)
    _PC['memo.q:same'] += k
    if s.tc_fb:
        _count_fb(e, run_)
    _f0_back(e, run_)


def _count_fb(e, run_):
    """The made fallbacks' counts, as the dense path's calls make them
    (`_tran_stageloop_c._count`'s for its attempts): the fallback itself, the
    C stage Newton's serves and assemblies on the stage KLU, the screens on
    the C's its serve marked, the fresh estimates."""
    s, fb = run_.s, int(run_.s.tc_fb)
    _PC['tcloop_c:fallbacks'] += fb
    _PC['radau.transform:fallback'] += fb
    _PC['radau_c:served'] += fb
    _PC['radau_c:assemblies'] += int(s.st_assemblies)
    if s.st_klu_on:
        _PC['radau_c:on_klu'] += int(s.st_klu_on)
    if run_.br is not None:
        _PC['branch_c:served'] += 3 * fb
        _PC['memo.C:memo'] += 3 * fb
        if e.fb[0].core.npat:
            _PC['branch_c:pattern'] += 3 * fb
    if s.adaptive:
        _PC['radau.est:fresh'] += fb


def _f0_back(e, run_):
    """The family's `early_f0` where the C made it: `(the point it is made
    for -- X's object --, f0)`."""
    s = run_.s
    if s.have_f0:
        e.fam.early_f0 = (e.st[3][int(s.f0_point)] if int(s.f0_point) < len(e.st[3])
                          else e.st[3][-1], run_.tbuf['tcf0'].copy())


def _count_screen(e, run_):
    """The screen's counts: its memo hits (`memo.C:memo`), the `cir.C`
    calls where the memo lacked a stage (each counting as one does), and
    the C screen's own -- three stages a screened attempt, none of their
    C's marked with the core's pattern (`_C_pattern` is the coupled
    Newton's): served whole up to `INPLACE_MIN_ENTRIES` entries, else
    declined `large` (the Python's screen makes it, the same verdict)."""
    from pycircuit.circuit import _cbuild
    s = run_.s
    if run_.br is None:
        return
    if s.tc_mcmemo:
        _PC['memo.C:memo'] += int(s.tc_mcmemo)
    for key, v in run_.ccounts.items():
        _PC[key] += v * int(s.tc_cpasses)
    a = 3 * int(s.tc_screened)
    if e.n * e.n > _cbuild.INPLACE_MIN_ENTRIES:
        if a:
            _PC['branch_c:large'] += a
        ## (the fallbacks' screens alone certified here: their C's marked)
        if s.tc_fbcert:
            _PC['branch_c:cert'] += int(s.tc_fbcert)
        return
    if a:
        _PC['branch_c:served'] += a
    if run_.br[0].ncert:
        _PC['branch_c:cert'] += run_.br[0].ncert
