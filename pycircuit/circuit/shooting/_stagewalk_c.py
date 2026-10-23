"""THE SHOOTING'S PLAIN STAGE WALK IN C: the steps of a shooting Newton's
walk under a fully implicit stage method on the walk's KLU
(`_PeriodWalks._walk_stage`, plain: dense, nothing kept, no period or event
column) in ONE C call (`csrc/stagewalk.c`), each step as the inner
transient's step and the walk's step map make it, bit for bit (pinned in
`tests/test_pinned_pairs.py`).

THE STEP.  The seeds (`_pred_stages`: pred.c on the walk's history of
nodes), the C stage Newton (radau.c, as `_tran_radau_c.solve` calls it),
the screen of each stage (branch.c, as `_branch_after_coupled` screens
them), the step's finish and its roll (`_finish_stage_step`,
`_roll_history`: the rings, `_pred_note`), the walk's step map on its
numeric (sensklu.c, `C_n P` numpy's own BLAS) and the live columns' drop
(`_LiveColumns.advance`).  The steps are made a chunk a call (`CHUNK_FIRST`
growing to `CHUNK_MAX`), the stage times and the sources at them and at
each step's end made for it first, the Python's way and in its order
(`_coupled_stage_context`'s `tn + c_i h`, `_stage_source`): the sources a
function of time, the shooting's own premise -- every walk re-enters the
period's times.

THE STATE.  Taken from the Python walk and the inner transient before the
call (the point, its charge and capacitance, the rings, the predictor's
nodes, the screen's running maximum, the live columns, both numerics) and
written back after it as the Python leaves them after the last step made:
the rings, `_pred_hist`, `_q_cache`, `_rk_Y`, `_iq`, `_Cmat`, `_Geq`,
`_C_cache` and `_branch_cmax` where it screens, `_dt` and its history, the
device memo's records of the last step's stages, the stage KLU's numeric
and counts, the walk numeric's numeric, generation and counts, the live
columns -- and each step made's records (`record`, `stage_record`).

A STOP hands back the step it stopped before, the Python walk going on from
it, its step made as the Python makes it alone: stopped before the step's
Newton, nothing of the step made; after the Newton converged, its result
held for the Python's repeat of the step (`_tran_radau_c.hold`: the stage
KLU's numeric moved once, its counts once); after the step map, the step
made and what the Python does next done here -- `advance` at the full
width on a map not finite, the Python's LinAlgError on a singular block or
a failed solve.

OWNS.  On the inner transient, `_stagewalk_c` (a `_Ctx`: the struct, its
buffers, the predictor state and the screen's struct it points at).
`_paths`: `pss.stagewalk:steps` (the steps made), `pss.stagewalk:<why>`
a decline or a stop (`_STOPS`).  Switch: `PYCIRCUIT_STAGEWALK_C=0`
(`STAGEWALK_C`, read on every call); `STATUS` says why the path is off.

DECLINES, nothing written that the Python step reads: the switch, or the
units' (`_tran_radau_c`, `_sensklu_c`, `_tran_pred_c`, the screen's); no
build; a piece it stands in for not its module's own, or shadowed on the
transient or the walk; the walk asking for more than the plain walk
(landed edges, event, period or error columns); PCNR, the cost transform or
the embedded estimate; the predictor off or a periodic row; a stateful
limiter, bypass; whatever `_tran_radau_c._ready` declines; another KLU
than the walk's; the screen's pieces patched; a history, rings or live columns not as the Python holds them.

History: `doc/shooting_history.md`, `_stagewalk_c.py`.
"""
import os
from types import SimpleNamespace

import numpy as np

from pycircuit.circuit import _cbuild, _paths, _tran_radau_c

#: `PYCIRCUIT_STAGEWALK_C=0` makes the Python walk every step
STAGEWALK_C = os.environ.get('PYCIRCUIT_STAGEWALK_C', '1') != '0'
STATUS = 'not loaded'
_driver = None
_MOD = {}
_PC = _paths.COUNTS
_no = _paths.no
_F64 = np.dtype(np.float64)

STAGEWALK_SRC = _cbuild.csrc('stagewalk.c')

STAGEWALK_CDEF = _cbuild.csrc('stagewalk.h')

#: the C's results, named in stagewalk.h (`_CODES`); `BAIL`: {code: why} of
#: the stops before a step; `_STOPS` names every result but done, and
#: `_SOURCE` the glue's own: a source that raised at a chunk's first step
_CODES = _cbuild.constants(STAGEWALK_CDEF)
BAIL = _cbuild.bails(_CODES, 'STAGEWALK_BAIL_')
_SOURCE = -1
_STOPS = {v: k[len('STAGEWALK_'):].lower() for k, v in _CODES.items()
          if k.startswith('STAGEWALK_') and k != 'STAGEWALK_DONE'}
_STOPS[_SOURCE] = 'source'
#: the stops after the step's converged Newton: its result held
_HELD = frozenset(_CODES['STAGEWALK_BAIL_' + w] for w in ('FINAL', 'REF', 'SCREEN', 'SENS'))
#: the stops whose last step made left the radau record's slots its own
_SLOTS = frozenset(_CODES['STAGEWALK_' + w] for w in (
    'DONE', 'BAIL_LIVE', 'BAIL_NODES', 'BAIL_PRED', 'NONFINITE', 'SINGULAR', 'SOLVE_FAILED'))
_SLOTS = _SLOTS | {_SOURCE}
_OUTS = ('Ya', 'qa', 'ia', 'Ca', 'Ga')
_DONE = _CODES['STAGEWALK_DONE']

#: the steps of a call: the first chunk's, doubling to the last
CHUNK_FIRST, CHUNK_MAX = 16, 256

#: the inner transient's pieces the walk stands in for: (attribute,
#: qualified name, module)
_M = 'pycircuit.circuit.'
_T_PIECES = (
    ('solve_timestep', 'Transient.solve_timestep', _M + 'transient'),
    ('_solve_timestep', 'Transient._solve_timestep', _M + 'transient'),
    ('_solve_timestep_rk', '_SequentialStages._solve_timestep_rk', _M + '_tran_stages'),
    ('_stage_source', '_SequentialStages._stage_source', _M + '_tran_stages'),
    ('_finish_stage_step', '_SequentialStages._finish_stage_step', _M + '_tran_stages'),
    ('_rk_step_coupled', '_RadauStages._rk_step_coupled', _M + '_tran_radau'),
    ('_coupled_stage_solver', '_RadauStages._coupled_stage_solver', _M + '_tran_radau'),
    ('_coupled_stage_context', '_RadauStages._coupled_stage_context', _M + '_tran_radau'),
    ('_finish_radau', '_RadauStages._finish_radau', _M + '_tran_radau'),
    ('_branch_after_coupled', '_BranchCheck._branch_after_coupled', _M + '_tran_branch'),
    ('_branch_screen', '_BranchCheck._branch_screen', _M + '_tran_branch'),
    ('_branch_on', '_BranchCheck._branch_on', _M + '_tran_branch'),
    ('_branch_structural_rank', '_BranchCheck._branch_structural_rank', _M + '_tran_branch'),
    ('_roll_history', '_RunHistory._roll_history', _M + '_tran_history'),
    ('_push_history', '_RunHistory._push_history', _M + '_tran_history'),
    ('_pred_promote', '_StagePredictor._pred_promote', _M + '_tran_predictor'),
    ('_pred_note', '_StagePredictor._pred_note', _M + '_tran_predictor'),
    ('_pred_stages', '_StagePredictor._pred_stages', _M + '_tran_predictor'),
    ('_pred_or', '_StagePredictor._pred_or', _M + '_tran_predictor'),
    ('_predict_state', '_StagePredictor._predict_state', _M + '_tran_predictor'),
    ('_predict_fast', '_StagePredictor._predict_fast', _M + '_tran_predictor'),
    ('_fit_fast', '_StagePredictor._fit_fast', _M + '_tran_predictor'),
    ('_pred_fast_nodes', '_StagePredictor._pred_fast_nodes', _M + '_tran_predictor'),
    ('_pred_degree', '_StagePredictor._pred_degree', _M + '_tran_predictor'),
    ('_source_at', '_CompanionModel._source_at', _M + '_tran_companion'),
    ('_memo_get', '_CompanionModel._memo_get', _M + '_tran_companion'),
    ('_memo_put', '_CompanionModel._memo_put', _M + '_tran_companion'),
    ('_q_at', '_CompanionModel._q_at', _M + '_tran_companion'),
    ('_C_at_state', '_CompanionModel._C_at_state', _M + '_tran_companion'),
    ('_rk_use_pcnr', '_PCNRSteps._rk_use_pcnr', _M + '_tran_pcnr'),
    ('_newton_option', '_StepNewton._newton_option', _M + '_tran_newton'))
#: the walk object's (`_InnerTransient`'s, `_PeriodWalks`')
_W_PIECES = (
    ('solve_timestep', '_InnerTransient.solve_timestep', _M + 'shooting._pss_inner'),
    ('_transient_step', '_InnerTransient._transient_step', _M + 'shooting._pss_inner'),
    ('_insert_refnode', '_InnerTransient._insert_refnode', _M + 'shooting._pss_inner'),
    ('_C_at', '_InnerTransient._C_at', _M + 'shooting._pss_inner'),
    ('_stage_step', '_PeriodWalks._stage_step', _M + 'shooting._pss_walks'))
#: the modules' functions: (module, name)
_F_PIECES = (('_tran_history', '_ring_push'), ('_tran_radau_c', 'solve'),
             ('_tran_radau_c', '_ready'), ('_tran_radau_c', 'final_block'),
             ('_tran_radau_c', 'hold'), ('_tran_branch_c', 'screen'),
             ('shooting._sensklu_c', 'step'), ('shooting._sensklu_c', 'block_step'),
             ('shooting._sensklu_c', 'full_step'))
#: names an instance shadow of would send the Python step elsewhere
_TR_NAMES = tuple(a for a, _q, _m in _T_PIECES) + ('_stage_klu',)
_W_NAMES = tuple(a for a, _q, _m in _W_PIECES) + ('_factorise', '_step_klu', '_transient')


def driver():
    """`(ffi, cfn, numpy's dgemm, numpy's dgemv)` once loaded, or None
    (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            npl = _cbuild.openblas(np)
            if npl is None:
                STATUS = "off (numpy's own OpenBLAS not found)"
                return None
            dgemm = _cbuild.address(npl.scipy_cblas_dgemm64_)
            dgemv = _cbuild.address(npl.scipy_cblas_dgemv64_)
        except (OSError, AttributeError) as e:
            STATUS = f'off ({e})'
            return None
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(STAGEWALK_SRC, STAGEWALK_CDEF,
                                                             'pcx_stagewalk')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, dgemm, dgemv)
        STATUS = 'c'
    return _driver or None


def _mods():
    """The modules the pieces are read from, imported once."""
    if 'mods' not in _MOD:
        import importlib

        from pycircuit.circuit import transient

        from . import _pss_walks
        _MOD['mods'] = (transient.Transient, _pss_walks._PeriodWalks,
                        tuple(importlib.import_module(_M + m) for m, _n in _F_PIECES))
    return _MOD['mods']


def _pieces():
    """The pieces as their modules define them -- read once every one is
    (`_paths.genuine`); None: one is a stand-in, read again next call."""
    g = _MOD.get('g')
    if g is None:
        from . import _pss_inner, _steps
        T, PW, fmods = _mods()
        IT = _pss_inner._InnerTransient
        objs = [(getattr(T, a), q, m) for a, q, m in _T_PIECES]
        objs += [(getattr(IT if q.startswith('_InnerTransient') else PW, a), q, m)
                 for a, q, m in _W_PIECES]
        objs += [(getattr(mod, nm), nm, _M + m) for mod, (m, nm) in zip(fmods, _F_PIECES)]
        objs += [(_steps._LiveColumns.advance, '_LiveColumns.advance', _M + 'shooting._steps')]
        if not all(_paths.genuine(*o) for o in objs):
            return None
        g = _MOD['g'] = tuple(o[0] for o in objs)
    return g


def _now():
    """What the Python walk calls now, in `_pieces`' order."""
    from . import _pss_inner, _steps
    T, PW, fmods = _mods()
    IT = _pss_inner._InnerTransient
    return (tuple(getattr(T, a) for a, _q, _m in _T_PIECES)
            + tuple(getattr(IT if q.startswith('_InnerTransient') else PW, a)
                    for a, q, _m in _W_PIECES)
            + tuple(getattr(mod, nm) for mod, (_m, nm) in zip(fmods, _F_PIECES))
            + (_steps._LiveColumns.advance,))


class _Ctx:
    """One inner transient's struct and buffers for one size -- the state,
    the scratch, a chunk's tables (`tab`) -- the predictor state its history
    feeds and the screen's struct."""

    __slots__ = ('bkeep', 'bs', 'buf', 'ffi', 'hist_len', 'i64', 'keep', 'm', 'n', 'pkeep',
                 'pkeep_pat', 'ps', 't', 'tab')

    def __init__(self, ffi, n, hist_len, pdrv, bdrv):
        from pycircuit.circuit import _tran_pred_c
        m = n - 1
        self.ffi, self.n, self.m, self.hist_len = ffi, n, m, hist_len
        t = self.t = ffi.new('stagewalk_t *')
        self.buf = buf = {
            'x': np.empty(n), 'xr': np.empty(m), 'qn': np.empty(n), 'Cn': np.empty(m * m),
            'ht': np.empty(8), 'hx': np.empty(8 * n), 'qring': np.empty(hist_len * n),
            'iqring': np.empty(hist_len * n), 'seed': np.empty(3 * n),
            'base': np.empty(m * m), 'Xs': np.empty(3 * m * m), 'Xs2': np.empty(3 * m * m),
            'Pc': np.empty(m * m), 'Pin': np.empty(m * m)}
        C = CHUNK_MAX
        self.tab = tab = {'hs': np.empty(C), 'tend': np.empty(C), 'tst': np.empty((C, 3)),
                          'U': np.empty((C, 4, n)), 'Xout': np.empty((C, m)),
                          'Yout': np.empty((C, 3, n))}
        self.i64 = np.empty(max(m, 1), dtype=np.int64)
        fb = ffi.from_buffer
        self.keep = {k: fb('double *', a) for k, a in list(buf.items()) + list(tab.items())}
        self.keep['cols'] = fb('int64_t *', self.i64)
        for k in ('x', 'xr', 'qn', 'Cn', 'ht', 'hx', 'qring', 'iqring', 'seed', 'base', 'Xs',
                  'Xs2', 'Pc', 'cols', 'hs', 'tend', 'tst', 'U', 'Xout', 'Yout'):
            setattr(t, k, self.keep[k])
        self.pkeep_pat = ()
        t.n, t.m, t.hist_len = n, m, hist_len
        pffi, pcfn, dgesv, dgemv = pdrv
        ps, _pbuf, pkeep = _tran_pred_c.new_state(pffi, dgesv, dgemv, n)
        ps.ht = pffi.cast('double *', int(ffi.cast('uintptr_t', t.ht)))
        ps.hx = pffi.cast('double *', int(ffi.cast('uintptr_t', t.hx)))
        self.ps, self.pkeep = ps, pkeep
        t.pred = ffi.cast('void *', int(pffi.cast('uintptr_t', ps)))
        t.pred_fn = ffi.cast('void *', int(pffi.cast('uintptr_t', pcfn)))
        self.bs, self.bkeep = None, ()
        if bdrv is not None:
            bffi = bdrv[0]
            self.bs = bffi.new('branch_t *')


def _plug(t, ffi, unit, uffi, state, entry):
    """The struct `t`'s pointers to a unit's struct and entry (another
    ffi's)."""
    setattr(t, unit, ffi.cast('void *', int(uffi.cast('uintptr_t', state))))
    setattr(t, unit + '_fn', ffi.cast('void *', int(uffi.cast('uintptr_t', entry))))


def _in(ffi, p, base, count):
    """The pointer `p` within the `count` doubles from `base`."""
    a, b = int(ffi.cast('uintptr_t', p)), int(ffi.cast('uintptr_t', base))
    return b <= a < b + 8 * count


def walk(walks, j0, x, times, hs, tab, live, walk_lu, record, stage_record):
    """The plain walk's steps from `j0` (its point `x`, reduced) in C, a
    chunk a call: `(the next step to make, the point there, the one before
    it)` with every step made recorded and the state written back -- or
    None (counted): no step made, the Python's next.  Raises the Python's
    LinAlgError where its step map raises one."""
    e = _eligible(walks, j0, x, times, hs, tab, live, walk_lu)
    if e.__class__ is str:
        return _no('pss.stagewalk:' + e)
    c = _bind(e, walks, x, live, walk_lu)
    run = _chunks(c, e, j0, x, times, hs, record is not None, stage_record is not None)
    _settle(c, e, run, walk_lu)
    if not run.made:
        return _no('pss.stagewalk:' + _STOPS[run.status])
    if run.status != _DONE:
        _no('pss.stagewalk:' + _STOPS[run.status])
    _PC['pss.stagewalk:steps'] += run.made
    _write_back(c, e, run, live)
    ## the records, each step made's, in its order
    for xr, tt in run.pts:
        record(xr, tt)
    if stage_record is not None:
        stage_record.extend(run.stg)
    _after(c, e, run, walks, tab, live, walk_lu)
    return run.j, run.prev, run.before


def _eligible(walks, j0, x, times, hs, tab, live, walk_lu):
    """What a call reads (a namespace), the point's charge made as
    `_coupled_stage_context` makes it -- or why it declines."""
    if not STAGEWALK_C:
        return 'off'
    drv = _driver or driver()
    if not drv:
        return 'driver'
    g = _MOD.get('g') or _pieces()
    if g is None or _now() != g:
        return 'patched'
    from pycircuit.circuit import _limiting, _tran_core, _tran_pred_c, _tran_predictor
    from pycircuit.circuit import linearsolver as LS

    from . import _pss_walks, _sensklu_c
    if not (_tran_radau_c.ENABLED and _sensklu_c.SENSKLU_C and _tran_pred_c.ENABLED
            and _tran_predictor.PRED_FAST):
        return 'units'
    tr = walks._transient()
    td, wd = tr.__dict__, walks.__dict__
    if any(nm in td for nm in _TR_NAMES[:-1]) or any(nm in wd for nm in _W_NAMES):
        return 'shadow'
    ## (the step map's own pieces, as `_sensklu_c.step` reads them)
    sg = _sensklu_c._MOD.get('g') or _sensklu_c._pieces()
    if (type(walk_lu) is not LS._KLUWalk or sg is None or _sensklu_c._now() != sg
            or not _pss_walks._own_factorise(walks)):
        return 'patched'
    if (wd.get('_landed_edges') is not None or wd.get('_want_event_cols')
            or wd.get('_want_dfdh') or wd.get('_want_lte')):
        return 'walk'
    if not (type(tab[0]) is np.ndarray and tab[0].shape == (3, 3)):
        return 'tableau'
    if (tr._rk_use_pcnr() or tr._newton_option(tr.par.radau_transform, 'radau_transform')
            or td.get('_rk_want_est')):
        return 'method'
    if tr.stage_predictor == 'off' or (getattr(tr, '_periodic_rows', None) or ()):
        return 'pred'
    nsteps = len(times) - 1
    if not 0 <= j0 < nsteps:
        return 'range'
    lims = td.get('_stateful_lims')
    if lims is None:
        lims = tr._stateful_lims = _limiting.stateful_limiters(tr.cir)
    nobypass = float(getattr(tr.epar, 'bypasstol', -1.0) or -1.0) < 0.0
    n, iref, integ = tr.cir.n, walks.irefnode, tr.base_integrator
    x_full = walks._insert_refnode(x)
    Amat = np.array(integ.A, dtype=float)
    ## (the charge at the point, as `_coupled_stage_context` reads it)
    rq = tr._memo_get(x_full)
    if rq is None or 'q' not in rq:
        P = _tran_core.passes(tr, x_full, 'q')
        qn = P['q'] if P is not None else np.asarray(tr.cir.q(x_full, tr.epar), dtype=float)
    else:
        qn = rq['q']
    ctx = SimpleNamespace(h=hs[min(j0, len(hs) - 1)], iref=iref, m=n - 1, qn=np.asarray(qn),
                          Amat=Amat)
    got = _tran_radau_c._ready(tr, ctx, None, None, lims, nobypass, tr.par.reltol,
                               step_memo=False)
    if got.__class__ is str:
        return 'radau_' + got
    rec, Tc, sk, _Y0, rdrv, _M2 = got
    if sk is None or sk is not walk_lu.solver or not sk._symbolic:
        return 'klu'
    sdrv, pdrv = _sensklu_c.driver(), _tran_pred_c.driver()
    if not sdrv or not pdrv:
        return 'units'
    scr = _screen_ready(tr, td, n, x_full)
    st = scr if scr.__class__ is str else _state_ready(td, n, live)
    if st.__class__ is str:
        return st
    return SimpleNamespace(
        drv=drv, tr=tr, n=n, m=n - 1, iref=iref, nsteps=nsteps, x_full=x_full, integ=integ,
        Amat=Amat, cvec=np.array(integ.C, dtype=float), qn=ctx.qn, rec=rec, Tc=Tc, sk=sk,
        rdrv=rdrv, sdrv=sdrv, pdrv=pdrv, scr=scr, hist=st[0], qlast=st[1], iqlast=st[2],
        P0=live.P)


def _screen_ready(tr, td, n, x_full):
    """The screen as `_branch_after_coupled` makes it: None (none made),
    `(its driver, the run's certificate data or None, the running maximum,
    the tolerance, the growth)`, or why it declines."""
    from pycircuit.circuit import _stamp_plan, _tran_branch_c
    from pycircuit.utilities.param import ParameterDict
    if not tr._branch_on():
        return None
    r0, scale0 = tr._branch_structural_rank()
    if not (r0 > 0 and scale0 > 0.0):
        return None
    bdrv = _tran_branch_c.driver() if _tran_branch_c.ENABLED else None
    own = _tran_branch_c._own() if bdrv else None
    if own is None or _tran_branch_c._patched(tr, own):
        return 'screen'
    ## (the varying `C`'s mark, as the Python's screen makes it: a constant
    ## one's verdict the Python keeps (`_branch_const`) is a screen of the
    ## same `C` at the same running maximum -- the C's)
    epoch = ParameterDict._epoch
    if (td.get('_branch_varying') is not epoch
            and _stamp_plan.constant_matrix(tr.cir, 'C', x_full) is None):
        tr._branch_varying = epoch
    known, data = tr._branch_cert_kept()
    cmax = td.get('_branch_cmax', 0.0)
    tol, growth = tr.BRANCH_SCREEN_TOL, tr.BRANCH_CERT_GROWTH
    if (type(cmax) is not float or type(tol) is not float or type(growth) is not float
            or (data is not None and (data.n != n or bdrv[7] is None))):
        return 'screen'
    return bdrv, (data if known else None), cmax, tol, growth


def _state_ready(td, n, live):
    """`(the predictor's nodes, the rings)` as the Python holds them, and
    the live columns -- or why they are not."""
    hist = td.get('_pred_hist') or []
    qlast, iqlast = td.get('_qlast'), td.get('_iqlast')
    if (type(hist) is not list or len(hist) > 8
            or any(type(e[1]) is not np.ndarray or e[1].shape != (n,) for e in hist)
            or type(qlast) is not np.ndarray or type(iqlast) is not np.ndarray
            or qlast.ndim != 2 or qlast.shape != iqlast.shape or qlast.shape[1] != n):
        return 'history'
    P0 = live.P
    if not (type(P0) is np.ndarray and P0.dtype == _F64 and P0.ndim == 2
            and P0.shape[0] == n - 1 and 1 <= P0.shape[1] <= n - 1):
        return 'live'
    return hist, qlast, iqlast


def _bind(e, walks, x, live, walk_lu):
    """The transient's context made or kept, the units plugged into its
    struct, the radau record set as `solve` sets it once a call, the state
    loaded from the Python's."""
    from pycircuit.circuit import _tran_branch_c

    from . import _sensklu_c
    ffi, n, m, tr = e.drv[0], e.n, e.m, e.tr
    td, hist_len = tr.__dict__, e.qlast.shape[0]
    c = td.get('_stagewalk_c')
    if (c is None or c.n != n or c.hist_len != hist_len
            or (e.scr is not None and c.bs is None)):
        c = td['_stagewalk_c'] = _Ctx(ffi, n, hist_len, e.pdrv, e.scr and e.scr[0])
    t, buf, rec, sk = c.t, c.buf, e.rec, e.sk
    rec.pre = None
    rs = rec.s
    _plug(t, ffi, 'radau', e.rdrv[0], rs, e.rdrv[1])
    sctx = sk.__dict__.get('_sensklu_c')
    if sctx is None or sctx.sym is not sk._symbolic:
        sctx = sk._sensklu_c = _sensklu_c._Ctx(e.sdrv[0], sk)
    _plug(t, ffi, 'sens', e.sdrv[0], sctx.t, e.sdrv[1])
    if e.scr is not None:
        bdrv, data, cmax, tol, growth = e.scr
        bffi, bs = bdrv[0], c.bs
        bs.n, bs.tol_rel, bs.ncert = n, tol, 0
        ## (the core's pattern where it has one, as the C loop's screen)
        pat = rec.core.arrays['pat'] if rec.core.npat else None
        c.pkeep_pat = () if pat is None else (pat, bffi.from_buffer('int64_t *', pat))
        if pat is None:
            bs.npat = 0
        else:
            bs.pat, bs.npat = c.pkeep_pat[1], len(pat)
        c.bkeep = _tran_branch_c.cert_bind(bffi, bs, data, growth, bdrv[7])
        _plug(t, ffi, 'branch', bffi, bs, bdrv[1])
        t.cmax = cmax
    else:
        t.branch = t.branch_fn = ffi.NULL
        t.cmax = 0.0
    t.dgemm, t.dgemv = ffi.cast('void *', e.drv[2]), ffi.cast('void *', e.drv[3])
    c.ps.deg = int(tr._pred_degree())
    c.ps.clamp = float(tr.PRED_CLAMP)
    ## the radau record: the transient's own output set (`own`: never
    ## viewed), the temperature, the tableau, the tolerances; its sources',
    ## seeds' and charge's pointers given back after the calls (`_chunks`)
    for k in _OUTS:
        setattr(rs, k, rec.own[1][k])
    rs.T = e.Tc
    for k, v in enumerate(e.Amat.ravel().tolist()):
        rs.A[k] = v
    e.tol = (tr.par.reltol, float(tr.par.vabstol), int(tr.par.maxiter))
    rs.reltol, rs.abstol, rs.maxiter = float(e.tol[0]), e.tol[1], e.tol[2]
    ## the state
    np.copyto(buf['x'], e.x_full)
    np.copyto(buf['xr'], x)
    np.copyto(buf['qn'], e.qn)
    np.copyto(buf['Cn'], np.asarray(walks._C_at(x)).ravel())
    t.hk = len(e.hist)
    for i, (tt, xx) in enumerate(e.hist):
        buf['ht'][i] = tt
        buf['hx'][i * n:(i + 1) * n] = xx
    np.copyto(buf['qring'], e.qlast.ravel())
    np.copyto(buf['iqring'], e.iqlast.ravel())
    k0 = e.P0.shape[1]
    buf['Pin'][:m * k0] = np.ascontiguousarray(e.P0).ravel()
    t.P = c.keep['Pin']
    t.k, t.ld, t.trans = k0, k0, 0
    t.have_cols = 0 if live.cols is None else 1
    if live.cols is not None:
        c.i64[:k0] = live.cols
    t.nadv = live.n
    t.live_on = 0 if live.off else 1
    t.walk_num = ffi.cast('void *', walk_lu._num or 0)
    t.klu_out = ffi.cast('void *', sk._numeric or 0)
    t.iref = e.iref
    return c


def _table(c, e, src, j, nst, times, hs):
    """A chunk's tables for steps `j ...`: the lengths, ends, stage times and
    the sources at them in the Python's order (`tn + c_i h`, then the
    end) -- the steps made before a source that raised (the Python's step
    raising where it raises)."""
    tb, cvec = c.tab, e.cvec
    for r in range(nst):
        jj = j + r
        h = hs[min(jj, len(hs) - 1)]
        te = times[jj + 1]
        tn = te - h
        tb['hs'][r], tb['tend'][r] = h, te
        try:
            for i in range(3):
                ti = tn + cvec[i] * h
                tb['tst'][r, i] = ti
                tb['U'][r, i] = src(ti)
            tb['U'][r, 3] = src(te)
        except Exception:                                      # noqa: BLE001
            return r
    return nst


def _chunks(c, e, j0, x, times, hs, want_pts, want_stg):
    """The calls, a chunk of steps each, its tables made first (the sources a
    function of time -- the shooting's own premise, every walk re-entering
    the period's times), the chunks growing: what they made (a namespace),
    each step's records kept for after the write-back."""
    t, tb, rs = c.t, c.tab, e.rec.s
    src = e.tr._stage_source(None)
    run = SimpleNamespace(status=_DONE, made=0, j=j0, prev=x, before=x, hl=[], pts=[], stg=[],
                          Y=None, cnt=[0] * 8)
    saved = (rs.qn, rs.u, rs.seed)
    size = CHUNK_FIRST
    try:
        while run.j < e.nsteps:
            plan = min(size, e.nsteps - run.j)
            nst = _table(c, e, src, run.j, plan, times, hs)
            if not nst:
                run.status = _SOURCE
                break
            t.steps, t.j = nst, 0
            run.status = e.drv[1](t)
            k = int(t.j)
            for i, v in enumerate((t.klu_refactors, t.klu_factors, t.klu_pivot_fallbacks,
                                   t.klu_residual_fallbacks, t.walk_refactors,
                                   t.walk_fallbacks, t.walk_factors, t.assemblies)):
                run.cnt[i] += v
            ## (the Python's own objects: a fresh point a step, the grid's
            ## times and lengths)
            for r in range(k):
                jj = run.j + r
                xr = tb['Xout'][r].copy()
                if want_pts:
                    run.pts.append((xr, times[jj + 1]))
                if want_stg:
                    run.stg.append((run.prev, hs[min(jj, len(hs) - 1)],
                                    [tb['Yout'][r, i].copy() for i in range(3)]))
                run.before, run.prev = run.prev, xr
                run.hl = run.hl[-1:] + [hs[min(jj, len(hs) - 1)]]
            if k:
                run.Y = [tb['Yout'][k - 1, i].copy() for i in range(3)]
            run.made += k
            run.j += k
            if run.status != _DONE or nst < plan:
                break
            size = min(2 * size, CHUNK_MAX)
    finally:
        rs.qn, rs.u, rs.seed = saved
        e.rec.fin = None
    return run


def _settle(c, e, run, walk_lu):
    """The numerics and their counts as the steps made leave them, whatever
    the stop; a converged Newton's result held for the Python's repeat of
    the step it stopped in (`_tran_radau_c.hold`); a limiter walk that
    stopped taken again, as `solve` takes it."""
    t, ffi, sk, cnt, rec = c.t, e.drv[0], e.sk, run.cnt, e.rec
    sk._numeric = int(ffi.cast('uintptr_t', t.klu_out)) or None
    sk.refactors += cnt[0]
    sk.factors += cnt[1]
    sk.pivot_fallbacks += cnt[2]
    sk.residual_fallbacks += cnt[3]
    walk_lu._num = int(ffi.cast('uintptr_t', t.walk_num)) or None
    walk_lu.gen += run.made
    walk_lu.refactors += cnt[4]
    walk_lu.fallbacks += cnt[5]
    walk_lu.factors += cnt[6]
    sk.numerics += cnt[6]
    _PC['radau_c:assemblies'] += cnt[7]
    if run.status in _HELD:
        r, tb = int(t.j), c.tab
        _tran_radau_c.hold(rec, sk, c.buf['qn'], tb['U'][r, :3], c.buf['seed'], e.Amat, e.Tc,
                           tb['hs'][r], *e.tol)
    elif (run.status == _CODES['STAGEWALK_BAIL_RADAU']
          and rec.s.status == _tran_radau_c._CODES['RADAU_BAIL_WALKSTOP']):
        rec.walk_ok = None


def _write_back(c, e, run, live):
    """The transient as the last step made leaves it, and the live
    columns."""
    from pycircuit.circuit import _tran_core
    t, buf, tr, n, m, ffi = c.t, c.buf, e.tr, e.n, e.m, e.drv[0]
    Y = run.Y
    xe, q3 = Y[2], buf['qring'][:n].copy()
    rbuf = e.rec.own[0]
    if run.status in _SLOTS:
        ## (no Newton since the last step made: its slots stand)
        a = int(t.last_assemblies) - 1
        C3, G3 = rbuf['Ca'][a, 2].copy(), rbuf['Ga'][a, 2].copy()
        if tr._memo_ok():
            for i in range(3):
                tr._memo_put(Y[i], {nm: rbuf[nm + 'a'][a, i].copy() for nm in ('q', 'i', 'C',
                                                                              'G')})
    else:
        ## (a later Newton wrote them: the point's own passes)
        vals = _tran_core.passes(tr, xe, 'CG')
        if vals is None:
            vals = {'C': np.asarray(tr.cir.C(xe, tr.epar), dtype=float),
                    'G': np.asarray(tr.cir.G(xe, tr.epar), dtype=float)}
        C3, G3 = vals['C'], vals['G']
    tr._rk_Y = Y
    tr._q_cache = (xe, q3)
    tr._iq = buf['iqring'][:n].copy()
    tr._Cmat = C3
    tr._Geq = e.Amat[2, 2] * run.hl[-1] * G3
    tr._effective_method = type(e.integ).__name__
    tr._companion_coeffs = None
    tr._pred_pending = None
    tr._qlast = buf['qring'].reshape(c.hist_len, n).copy()
    tr._iqlast = buf['iqring'].reshape(c.hist_len, n).copy()
    tr._pred_hist = [(float(buf['ht'][i]), buf['hx'][i * n:(i + 1) * n].copy())
                     for i in range(int(t.hk))]
    if e.scr is not None:
        tr._branch_cmax = float(t.cmax)
        tr._C_cache = (xe, C3)
    tr._dt_last2 = run.hl[-2] if len(run.hl) >= 2 else tr._dt_last
    tr._dt = tr._dt_last = run.hl[-1]
    tr._is_first_step = False
    tr._no_history = False
    k1 = int(t.k)
    if _in(ffi, t.P, c.keep['Pc'], m * m):
        live.P = buf['Pc'][:m * k1].reshape(k1, m).T.copy(order='F')
    elif int(t.trans):
        ## (the step map's last stage block, column-major in the buffer it
        ## was made in)
        Xp = buf['Xs'] if _in(ffi, t.P, c.keep['Xs'], 3 * m * m) else buf['Xs2']
        live.P = Xp.reshape(m, 3 * m)[:k1, 2 * m:].T.copy(order='F')
    live.cols = c.i64[:k1].copy() if t.have_cols else None
    live.n = int(t.nadv)


def _after(c, e, run, walks, tab, live, walk_lu):
    """What the Python does after the last step's map where the C stopped
    there: `advance` at the full width on a map not finite (as the Python's
    `_sensklu_c.step` hands it, the map in the buffer the carried P is not
    in), the LinAlgError of a singular block or a failed solve."""
    from . import _sensklu_c
    if run.status == _CODES['STAGEWALK_NONFINITE']:
        t, buf, m, ffi = c.t, c.buf, e.m, e.drv[0]
        Xw = buf['Xs2'] if _in(ffi, t.P, c.keep['Xs'], 3 * m * m) else buf['Xs']
        Pn = Xw.reshape(m, 3 * m)[:int(t.k), 2 * m:].T
        Cn = np.asarray(walks._C_at(run.before))
        live.advance(Pn, _sensklu_c.full_step(Cn, run.hl[-1], tab, walk_lu))
    elif run.status == _CODES['STAGEWALK_SINGULAR']:
        raise np.linalg.LinAlgError('Singular matrix')
    elif run.status == _CODES['STAGEWALK_SOLVE_FAILED']:
        raise np.linalg.LinAlgError('klu solve failed')
