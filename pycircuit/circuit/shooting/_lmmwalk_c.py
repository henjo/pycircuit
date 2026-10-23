"""THE SHOOTING'S PLAIN MULTISTEP WALK IN C: the steps of a shooting
Newton's walk under a multistep companion (`_PeriodWalks._walk_lmm`, plain:
dense, nothing kept, no period or event column) in C, a chunk a call --
each step the C stepping loop's own attempt on the walk's grid (`csrc/
loop.c`'s walk: `take` -- the companion with its order drop, the
predictor's seed, the C Newton with its J formed, the branch screen -- and
`accept` -- the rings, the predictor's node, the point, the step record),
then the walk's step map on the monodromy's live columns (lmmsens.c, as
`_lmmsens_c.step` makes it) -- bit for bit (pinned in
`tests/test_pinned_pairs.py`).

THE SOURCES ARE THE PYTHON'S: a chunk's table is made first, each step's
row its `_source_at` at the step's end as the Python's step makes it (its
step memo opened, the solve's source memo read and written, their counts
made: `_tran_newton_c._source`'s pass) -- the shooting's own premise, a
walk's sources a function of time.  The rows of steps no call made are
taken back (`_undo`, as `_tran_newton_c._undo_source` takes back an
attempt's): their counts and the entries they added to the solve's memo.

ENTRY, as the C loop enters: only where the C Newton's last call served at
the watch counter as it stands (`_newton_c.served`), its record bound and
its checks made -- the walk's first steps are the Python's.  Then the C
loop's checks of the run (its pieces, the circuit, the integrator: Gear-2,
Euler or the trapezoid; the predictor, the screen: `_tran_loop_c`'s) and
the walk's own: a plain dense walk (no landed edges, LTE, period or event
columns, nothing kept, no `_record_cj`), its pieces genuine (the walk
object's read off its class), the step map's toolkit and solver
(`_lmmsens_c`), the live monodromy as `_LiveMonodromy` holds it with every
array as the step map takes it (C-contiguous: the views a prune leaves are
the Python's to read).

THE STATE.  Read from the inner transient and the walk before a call --
the rings, the predictor's nodes, the cached C at the last point, the step
record, the screen's running maximum, the monodromy's live columns, the
capacitance ring -- and written back after it as the Python leaves them
after the last step made: the transient's as `_tran_loop_c._write_back`
writes a stepping loop's (the C Newton's serve, `_roll_history`; the step
lengths the grid's own objects), the walk's point, the one before, its ring
`Cs`, the monodromy, its step's `_Jf`, `_Geq`, `_C` and `_coeffs` -- each
step's point recorded (`record`) and the counts the Python's steps make
made.

A STOP hands back the step it stopped before, nothing of it made: a unit's
hand-back -- the predictor's general path, the Newton, the screen -- leaves
its state as the C loop's does; a source that raises is the Python's step's
to raise.  After a step whose map handed back (singular, not finite) the
step is made and its map is made as `_walk_lmm` makes it (`_lmmsens_c.step`,
then `_LiveMonodromy.step`); after a step whose map left a column dead, the
state one step on is the Python's `_LiveMonodromy.advance` (its prune, its
views).

OWNS: on the walk object `_lmmwalk_c` (a `_Ctx`: the loop's run bound for
the walk, the walk's buffers, its step map's struct).  `_paths`:
`pss.lmmwalk:steps`, the declines and stops (`pss.lmmwalk:<why>`).
Switch: `PYCIRCUIT_LMMWALK_C=0` (`LMMWALK_C`, read on every call).

History: `doc/shooting_history.md`, `_lmmwalk_c.py`.
"""
import os
from itertools import repeat
from types import SimpleNamespace

import numpy as np

from pycircuit.circuit import _paths, _tran_loop_c

#: `PYCIRCUIT_LMMWALK_C=0` makes the Python walk every step
LMMWALK_C = os.environ.get('PYCIRCUIT_LMMWALK_C', '1') != '0'
#: the steps of a call: the first chunk's, doubling to the last
CHUNK_FIRST, CHUNK_MAX = 16, 256
_MOD = {}
_PC = _paths.COUNTS
_no = _paths.no
_F64 = np.dtype(np.float64)
#: a step length as the walk's grid hands it (`hs[j]`) or a run's (float)
_FLOATS = (float, np.float64)
_CODES = _tran_loop_c._CODES
_DONE = _CODES['LOOP_DONE']
_MAP, _PRUNE = _CODES['LOOP_WALK_MAP'], _CODES['LOOP_WALK_PRUNE']
#: the glue's own stops: a step's source raised (or was not a vector), the
#: watch counter moved while a table was made
_SOURCE, _WATCH = -1, -2
_STOPS = dict(_tran_loop_c.BAIL)
_STOPS.update({_MAP: 'walk_map', _PRUNE: 'walk_prune', _SOURCE: 'source', _WATCH: 'watch'})
#: (where the inner transient keeps no statistics: what the call's timing
#: adds to, read by nothing)
_NOSTATS = SimpleNamespace(solve_seconds=0.0)

_M = 'pycircuit.circuit.'
#: the walk's pieces the C stands in for beyond the C loop's
#: (`_tran_loop_c._PIECES`) and the step map's (`_lmmsens_c`'s): (where it
#: is read from -- 'W' the walk object's class -- the attribute, its
#: qualified name, its module)
_W_PIECES = (
    ('W', 'solve_timestep', '_InnerTransient.solve_timestep', 'shooting._pss_inner'),
    ('W', '_transient_step', '_InnerTransient._transient_step', 'shooting._pss_inner'),
    ('W', '_insert_refnode', '_InnerTransient._insert_refnode', 'shooting._pss_inner'),
    ('LM', 'step', '_LiveMonodromy.step', 'shooting._steps'),
    ('LM', 'advance', '_LiveMonodromy.advance', 'shooting._steps'),
    ('LM', '_prune', '_LiveMonodromy._prune', 'shooting._steps'),
    ('Gear2', 'companion_coefficients', 'Gear2Integrator.companion_coefficients', 'integrator'),
    ('Euler', 'companion_coefficients', 'EulerIntegrator.companion_coefficients', 'integrator'),
    ('Trap', 'companion_coefficients', 'TrapezoidalIntegrator.companion_coefficients',
     'integrator'),
    ('integrator', 'bdf2_alphas', 'bdf2_alphas', '_lte_kernels'),
    ('lmmsens_c', 'step', 'step', 'shooting._lmmsens_c'))
#: names an instance attribute of on the walk object would send its step,
#: or its map, elsewhere
_W_SHADOWS = ('solve_timestep', '_transient_step', '_insert_refnode', '_transient',
              '_get_linearsolver')


def _owners():
    """The objects the walk's pieces are read from (the walk object's class
    aside), imported once."""
    if 'owners' not in _MOD:
        from pycircuit.circuit import integrator as I

        from . import _lmmsens_c, _pss_inner, _steps
        _MOD['owners'] = {'W': _pss_inner._InnerTransient, 'LM': _steps._LiveMonodromy,
                          'Gear2': I.Gear2Integrator, 'Euler': I.EulerIntegrator,
                          'Trap': I.TrapezoidalIntegrator, 'integrator': I,
                          'lmmsens_c': _lmmsens_c}
    return _MOD['owners']


def _now(walks):
    """What the Python walk calls now, in `_W_PIECES`' order -- the walk
    object's from its class."""
    own = _owners()
    W = type(walks)
    return tuple(getattr(W if o == 'W' else own[o], a) for o, a, _q, _m in _W_PIECES)


def _pieces():
    """The walk's pieces as their modules define them -- read once every one
    is (`_paths.genuine`); None: one is a stand-in, read again next call."""
    g = _MOD.get('g')
    if g is None:
        own = _owners()
        objs = [(getattr(own[o], a), q, _M + m) for o, a, q, m in _W_PIECES]
        if not all(_paths.genuine(*o) for o in objs):
            return None
        g = _MOD['g'] = tuple(o[0] for o in objs)
    return g


class _Ctx:
    """One walk object's C walk for one size: the loop's run bound for it
    (`_tran_loop_c._bind`, no sources unit: the table's rows), the walk's
    buffers and its step map's struct."""

    __slots__ = ('buf', 'keep', 'ldrv', 'lm', 'm', 'n', 'run', 'wcap')

    def __init__(self, run, n, m, wcap, ldrv):
        from . import _lmmsens_c
        self.run, self.n, self.m, self.wcap, self.ldrv = run, n, m, wcap, ldrv
        ffi, s = run.ffi, run.s
        C = CHUNK_MAX
        self.buf = buf = {'Jf': np.empty(m * m), 'wh': np.empty(C), 'wt': np.empty(C),
                          'wu': np.empty(C * n)}
        for i in range(3):
            buf[f'Cr{i}'] = np.empty(m * m)
            buf[f'Pb{i}'] = np.empty(m * wcap)
        for i in range(2):
            buf[f'Pqb{i}'] = np.empty(m * wcap)
        fb = ffi.from_buffer
        self.keep = keep = {k: fb('double *', a) for k, a in buf.items()}
        s.Jf, s.wh, s.wt, s.wu = keep['Jf'], keep['wh'], keep['wt'], keep['wu']
        for i in range(3):
            s.Cr[i], s.Pb[i] = keep[f'Cr{i}'], keep[f'Pb{i}']
        for i in range(2):
            s.Pqb[i] = keep[f'Pqb{i}']
        lffi, lcfn, gesv, gemm = ldrv
        self.lm = _lmmsens_c._Ctx(lffi, gesv, gemm, m, wcap)
        _tran_loop_c._plug(s, ffi, 'map', lffi, self.lm.t, lcfn)
        s.walk, s.m = 1, m


def walk(walks, j0, x, times, hs, mono, Cs, record):
    """The plain multistep walk's steps from `j0` (its point `x`, reduced;
    the monodromy `mono`, the capacitance ring `Cs`) in C, a chunk a call:
    `(the next step to make, the point there, the one before it, the ring)`
    with every step made recorded and the state written back -- or None
    (counted): no step made, the Python's next.  Raises where the Python's
    map of a step made raises (its LinAlgError)."""
    e = _eligible(walks, j0, x, times, hs, mono, Cs)
    if e.__class__ is str:
        return _no('pss.lmmwalk:' + e)
    c = _bind(walks, e)
    if c is None:
        return _no('pss.lmmwalk:bind')
    _load(c, e, mono, Cs)
    run = _chunks(c, e, j0, x, times, hs, record is not None)
    _settle(c, e, run)
    stop = 'pss.lmmwalk:' + _STOPS.get(run.status, str(run.status))
    if run.made == 0:
        return _no(stop)
    if run.status != _DONE:
        _no(stop)
    _PC['pss.lmmwalk:steps'] += run.made
    _write_back(c, e, run)
    ## the records, each step made's, in its order -- before a map the
    ## Python makes, which may raise as the Python's walk does after them
    for xr, tt in run.pts:
        record(xr, tt)
    Cs = _walk_back(walks, c, e, run, mono)
    return run.j, run.prev, run.before, Cs


def _eligible(walks, j0, x, times, hs, mono, Cs):
    """What a call reads (a namespace) -- or why it declines."""
    if not (LMMWALK_C and _tran_loop_c.ENABLED):
        return 'off'
    drv = _tran_loop_c.driver()
    if drv is None:
        return 'driver'
    from pycircuit.circuit import _tran_newton_c, _tran_pred_c, _watch

    from . import _lmmsens_c, _steps
    ldrv = _lmmsens_c.driver() if _lmmsens_c.LMMSENS_C else None
    if not ldrv:
        return 'map'
    if not _tran_pred_c.ENABLED or _tran_pred_c.driver() is None:
        return 'pred'
    g = _MOD.get('g') or _pieces()
    if g is None or _now(walks) != g:
        return 'patched'
    tr = walks._transient()
    td, wd = tr.__dict__, walks.__dict__
    now = _watch.now()
    rec = td.get('_newton_c')
    ## (the C Newton's last call served, nothing watched moved since: its
    ## record bound and its checks made, as the C loop enters)
    if now == -1 or rec is None or rec.served != now:
        return 'served'
    if rec.s.fold_j != _tran_newton_c.FOLD_J:
        return 'newton'
    if any(nm in wd for nm in _W_SHADOWS):
        return 'shadow'
    ## (what the Python's step and walk read off the walk object, as they
    ## read it)
    if (getattr(walks, '_landed_edges', None) is not None
            or getattr(walks, '_want_event_cols', False) or walks._want_dfdh
            or walks._want_lte or walks._record_cj or wd.get('_skip_step_mats')):
        return 'walk'
    why = _map_ready(walks, mono, Cs, _lmmsens_c, _steps)
    if why:
        return why
    O = _tran_loop_c._own()
    why = _run_ready(tr, td, O)
    if why:
        return why
    n, m, iref = tr.cir.n, tr.cir.n - 1, walks.irefnode
    nsteps = len(times) - 1
    if not (0 <= j0 < nsteps and type(x) is np.ndarray and x.dtype is _F64 and x.shape == (m,)):
        return 'range'
    screen = _tran_loop_c._screen(tr, td, rec, n)
    if screen is None:
        return 'screen'
    st = _state_ready(tr, td, n)
    if st.__class__ is str:
        return st
    return SimpleNamespace(
        drv=drv, ldrv=ldrv, tr=tr, rec=rec, O=O, n=n, m=m, iref=iref, now=now, nsteps=nsteps,
        x_full=walks._insert_refnode(x), screen=None if screen[0] is None else screen,
        wcap=mono.width, ql=st[0], iql=st[1], hist=st[2], cx=st[3], cC=st[4], dl=st[5], j0=j0,
        hs=hs, u_keys=_tran_newton_c._mods()['u_keys'])


def _map_ready(walks, mono, Cs, LS, ST):
    """Why the step map is not the C's -- `_lmmsens_c.step`'s checks, its
    inputs as the Python holds them (a prune's views decline: the Python's
    maps read them) -- or ''."""
    g = LS._MOD.get('g') or LS._pieces()
    if g is None or LS._now() != g:
        return 'patched'
    if LS._which(walks.toolkit, walks._get_linearsolver(), g):
        return 'map_solver'
    if type(mono) is not ST._LiveMonodromy or mono.off:
        return 'mono'
    m = walks.cir.n - 1
    P0, P1, Pq = mono.P[0], mono.P[1], mono.Pq
    w = P0.shape[1] if type(P0) is np.ndarray and P0.ndim == 2 else 0
    if not (2 <= m <= LS.MAXM and 2 <= w <= mono.width
            and all(LS._array(a, (m, w)) for a in (P0, P1, Pq))
            and len(Cs) == 2 and all(LS._array(C, (m, m)) for C in Cs)):
        return 'mono'
    return ''


def _run_ready(tr, td, O):
    """Why the run is not one the C loop's attempts make (a reason), or ''
    -- `_tran_loop_c._why_not_the_run`'s checks a walk's steps take (no
    breakpoint, no accept hook: the shooting calls neither)."""
    LC = _tran_loop_c
    if O is None or type(tr) is not O['T'] or LC._pieces(O) != O['chain']:
        return 'patched'
    stats = getattr(tr, 'statistics', None)
    if stats is not None and type(stats) is not O['stats']:
        return 'stats'
    if any(k in td for k in LC._SHADOWS):
        return 'shadow'
    if getattr(tr, '_ev_rows', None) is not None or getattr(tr, '_periodic_rows', None):
        return 'events'
    if getattr(tr, '_continuation_rescue', False) or td.get('_stateful_lims'):
        return 'state'
    cir = tr.cir
    if type(cir).u is not O['Sub'].u or any(k in cir.__dict__ for k in LC._CIRCUIT_SHADOWS):
        return 'circuit'
    base = tr.base_integrator
    if O['integrators'].get(type(base)) is None or 'check_order_drop' in base.__dict__:
        return 'integrator'
    P = O['predictor']
    if (not P.PRED_FAST or tr.stage_predictor != 'on' or tr.PRED_HIST_MAX != LC._HIST_MAX
            or type(tr)._fit is not P._FIT or type(tr.PRED_CLAMP) is not float
            or not 1 <= tr._pred_degree() <= LC._MAX_DEG):
        return 'predictor'
    if not O['companion'].U_MEMO:
        return 'umemo'
    return ''


def _state_ready(tr, td, n):
    """`(the rings, the predictor's nodes, the cached C's point and C, the
    last step's length)` as the Python leaves them -- `_tran_loop_c._read`'s
    checks -- or why not."""
    ql, iql = td.get('_qlast'), td.get('_iqlast')
    if not (type(ql) is np.ndarray and ql.dtype is _F64 and ql.ndim == 2 and ql.shape[1] == n
            and ql.shape[0] >= 2 and type(iql) is np.ndarray and iql.dtype is _F64
            and iql.shape == ql.shape):
        return 'rings'
    hist = td.get('_pred_hist') or []
    if (type(hist) is not list or len(hist) > _tran_loop_c._HIST_MAX
            or td.get('_pred_pending') is not None
            or not all(type(e) is tuple and len(e) == 2 and type(e[0]) is float
                       and type(e[1]) is np.ndarray and e[1].dtype is _F64 and e[1].shape == (n,)
                       for e in hist)):
        return 'history'
    cached = td.get('_C_cache')
    if cached is None or float(getattr(tr.epar, 'bypasstol', -1.0) or -1.0) >= 0.0:
        return 'cache'
    cx, cC = cached
    if not (type(cx) is np.ndarray and cx.dtype is _F64 and cx.shape == (n,)
            and type(cC) is np.ndarray and cC.dtype is _F64 and cC.shape == (n, n)):
        return 'cache'
    dl, dl2 = tr._dt_last, tr._dt_last2
    if (type(dl) not in _FLOATS or not (dl2 is None or type(dl2) in _FLOATS)
            or td.get('_no_history')):
        return 'history'
    return ql, iql, hist, cx, cC, dl


def _bind(walks, e):
    """The walk object's `_Ctx`, made again where the Newton's record, the
    screen's presence or certificate, the sizes, the monodromy's width or
    the step map's driver moved."""
    wd = walks.__dict__
    c = wd.get('_lmmwalk_c')
    hist_len = e.ql.shape[0]
    run = None if c is None else c.run
    if (c is None or run.nrec is not e.rec or c.n != e.n or c.wcap < e.wcap
            or c.ldrv is not e.ldrv or (run.br is None) != (e.screen is None)
            or (e.screen is not None and (run.cert[0] is not e.screen[3]
                                          or run.cert[1] != e.screen[4]))
            or len(run.buf['qring']) != hist_len * e.n):
        run = _tran_loop_c._bind(e.tr, e.rec, None, e.n, hist_len, e.screen, None)
        if run is None:
            return None
        c = wd['_lmmwalk_c'] = _Ctx(run, e.n, e.m, e.wcap, e.ldrv)
    return c


def _load(c, e, mono, Cs):
    """The struct loaded with the state as the Python leaves it -- the
    rings, the predictor's nodes, the cached C, the step record, the run's
    constants (`_tran_loop_c._load`'s, a fixed grid's, no judge) -- and the
    walk's: the monodromy's live columns, the capacitance ring."""
    tr, n, m = e.tr, e.n, e.m
    run = c.run
    s, buf = run.s, run.buf
    buf['qring'][:] = e.ql.reshape(-1)
    buf['iqring'][:] = e.iql.reshape(-1)
    buf['xlast'][:] = e.x_full
    buf['xatt'][:] = e.cx
    s.hk = len(e.hist)
    for i, (tt, xx) in enumerate(e.hist):
        buf['ht'][i] = tt
        buf['hx'][i * n:(i + 1) * n] = xx
    buf['cC0'][:] = e.cC.reshape(-1)
    s.cur = 0
    s.cmax = tr.__dict__.get('_branch_cmax', 0.0)
    if run.br is not None:
        run.br[0].ncert = 0
    s.landing = s.order_drop = s.imposed = 0
    s.is_first = int(bool(tr._is_first_step))
    s.rejects = s.point_retries = 0
    dl2 = tr._dt_last2
    s.dt_last, s.have_dt_last = e.dl, 1
    s.dt_last2, s.have_dt_last2 = (dl2, 1) if dl2 is not None else (0.0, 0)
    s.adaptive, s.early_th, s.newton_reject = 0, 0.0, 0
    s.tend, s.max_step, s.minstep, s.minbreak = float('inf'), float('inf'), 0.0, 0.0
    base = tr.base_integrator
    s.base = e.O['integrators'][type(base)]
    s.zs_ratio = float(e.O['zs_ratio'])
    s.shrink_guard = int(bool(getattr(base, 'shrink_guard', True)))
    s.have_hm = 0
    ## the walk's: the monodromy (P_n, P_{n-1}, Pq_n on their live columns,
    ## row-major), the capacitance ring (C_n, C_{n-1}, reduced)
    s.iref = e.iref
    w = mono.P[0].shape[1]
    c.buf['Pb0'][:m * w] = mono.P[0].reshape(-1)
    c.buf['Pb1'][:m * w] = mono.P[1].reshape(-1)
    c.buf['Pqb0'][:m * w] = mono.Pq.reshape(-1)
    s.p0, s.p1, s.pq, s.w, s.mono_n = 0, 1, 0, w, mono.n
    c.buf['Cr0'][:] = Cs[0].reshape(-1)
    c.buf['Cr1'][:] = Cs[1].reshape(-1)
    s.cn, s.cp = 0, 1


def _table(c, e, j, nst, times, hs):
    """A chunk's tables for steps `j ...`: the lengths, the ends and the
    sources at each end as the Python's step makes them (`_source_at` in
    its own step memo: `_tran_newton_c._source`'s pass) -- and what each
    row added, `(the counts before, its key in the solve's source memo,
    whether the memo held it)`, for `_undo`.  The rows before a source that
    raised or is not a vector: the Python's step raises where it raises."""
    from pycircuit.circuit import _stamp_plan
    tr, n, keys, buf = e.tr, e.n, e.u_keys, c.buf
    memo = _stamp_plan._MEMO.get(id(tr.cir)) if _stamp_plan._MEMO else None
    an = tr.par.analysis
    rows = []
    for r in range(nst):
        jj = j + r
        te = times[jj + 1]
        key = (te, an)
        try:
            had = memo is not None and key in memo[0]
        except TypeError:
            break
        cvals = tuple(map(_PC.get, keys, repeat(0)))
        u = None
        tr._u_memo = {}
        try:
            u = tr._source_at(te, None)
        except Exception:                                      # noqa: BLE001
            u = None
        finally:
            tr._u_memo = None
        if not (type(u) is np.ndarray and u.dtype == np.float64 and u.shape == (n,)):
            _undo(e, [(cvals, key, had)], 0)
            break
        buf['wh'][r] = hs[min(jj, len(hs) - 1)]
        buf['wt'][r] = te
        buf['wu'][r * n:(r + 1) * n] = u
        rows.append((cvals, key, had))
    return rows


def _undo(e, rows, k):
    """The rows from `k` on taken back (no step read them): the counts as
    they were before row `k`, the entries they added to the solve's source
    memo removed (one begun again among them stays begun: its vectors are
    the same)."""
    if k >= len(rows):
        return
    from pycircuit.circuit import _stamp_plan
    memo = _stamp_plan._MEMO.get(id(e.tr.cir)) if _stamp_plan._MEMO else None
    if memo is not None:
        for _cv, key, had in rows[k:]:
            if not had:
                memo[0].pop(key, None)
    for nm, v in zip(e.u_keys, rows[k][0]):
        _PC[nm] = v


def _chunks(c, e, j0, x, times, hs, want_pts):
    """The calls, a chunk of steps each, its tables made first (`_table`),
    the chunks growing: what they made (a namespace), each step's point and
    time kept for the records."""
    from pycircuit.circuit import _watch
    s, buf, n, iref, rec = c.run.s, c.run.buf, e.n, e.iref, e.rec
    stats = getattr(e.tr, 'statistics', None) or _NOSTATS
    run = SimpleNamespace(status=_DONE, made=0, j=j0, prev=x, before=x, pts=[], iters=0,
                          attempts=0, screens=0, pred_fast=0, maps=0)
    s.t = float(times[j0])
    size = CHUNK_FIRST
    while run.j < e.nsteps:
        plan = min(size, e.nsteps - run.j)
        rows = _table(c, e, run.j, plan, times, hs)
        if _watch.now() != e.now:
            ## (a source's check moved the watch: the Python's step re-checks)
            _undo(e, rows, 0)
            run.status = _WATCH
            break
        if not rows:
            run.status = _SOURCE
            break
        s.wsteps, s.wj = len(rows), 0
        if rec.klu is not None:
            ## (the run's numeric as it stands, for the sparse full Newton)
            num = rec.klu._numeric
            rec.s.klu_num = rec.ffi.cast('void *', num) if num else rec.ffi.NULL
        run.status = _tran_loop_c._call(e.drv, rec, s, stats)
        k = int(s.k)
        _undo(e, rows, k)
        for r in range(k):
            xf = buf['X'][r * n:(r + 1) * n]
            xr = np.concatenate((xf[:iref], xf[iref + 1:]))
            if want_pts:
                run.pts.append((xr, times[run.j + r + 1]))
            run.before, run.prev = run.prev, xr
        run.made += k
        run.j += k
        run.iters += int(s.iters)
        run.attempts += int(s.attempts)
        run.screens += int(s.screens)
        run.pred_fast += int(s.pred_fast)
        run.maps += int(s.maps)
        if run.status != _DONE:
            break
        if len(rows) < plan:
            run.status = _SOURCE
            break
        size = min(2 * size, CHUNK_MAX)
    return run


def _settle(c, e, run):
    """The counts and KLU tallies the Python's steps make, as the C loop
    re-creates them (`_tran_loop_c._count`) -- the C Newton's serves, the
    predictor's, the charge memo's, the screens', the step maps' (the
    sources' the tables' own calls made)."""
    rec, a = e.rec, run.attempts
    if rec.klu is not None:
        if rec.s.chord:
            rec.klu.numerics += a
        else:
            rec.klu.refactors += run.iters
        _PC['newton_c:sparse'] += a
    _PC['newton_c:served'] += a
    _PC['newton_c:iterations'] += run.iters
    _PC['newton_c:jfold'] += a
    if rec.s.chord:
        _PC['chord:chosen'] += a
    if run.pred_fast:
        _PC['pred:fast'] += run.pred_fast
    _PC['memo.q:same'] += run.made
    if run.screens:
        _PC['branch_c:served'] += run.screens
        _PC['memo.C:equal'] += run.screens
        if c.run.npat:
            _PC['branch_c:pattern'] += run.screens
        if c.run.br[0].ncert:
            _PC['branch_c:cert'] += c.run.br[0].ncert
    if run.maps:
        _PC['pss.lmm:served'] += run.maps
    stats = getattr(e.tr, 'statistics', None)
    if stats is not None:
        stats.newton_iterations += run.iters


def _write_back(c, e, run):
    """The inner transient as the last step made leaves it: the C Newton's
    serve (`_tran_newton_c._serve`'s writes) and `_roll_history`, as the C
    loop writes a stepping loop's back (`_tran_loop_c._write_back`) -- the
    step lengths the grid's own objects, as the Python's steps leave them."""
    tr, n, s, buf, rec = e.tr, e.n, c.run.s, c.run.buf, e.rec
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
    ## (the lengths the last step read: the grid's `hs` entries -- before
    ## the call's first step, the run's own)
    jl, nh = run.j - 1, len(e.hs)
    h = e.hs[min(jl, nh - 1)]
    h_last = e.hs[min(jl - 1, nh - 1)] if jl > e.j0 else e.dl
    base, O = tr.base_integrator, e.O
    active = (O['Euler']() if s.att_formula == _tran_loop_c._ABI['PCX_FORMULA_EULER']
              and type(base) is not O['Euler'] else base)
    tr.active_integrator = active
    tr._companion_coeffs = active.companion_coefficients(h, h_last)
    tr._effective_method = type(active).__name__
    tr._dt = h
    tr._pred_pending = None
    hist_len = e.ql.shape[0]
    tr._qlast = buf['qring'].reshape(hist_len, n).copy()
    tr._iqlast = buf['iqring'].reshape(hist_len, n).copy()
    tr._pred_hist = [(float(buf['ht'][i]), buf['hx'][i * n:(i + 1) * n].copy())
                     for i in range(s.hk)]
    tr._dt_last, tr._dt_last2 = h, h_last
    tr._no_history = False
    tr._is_first_step = bool(s.is_first)
    if c.run.br is not None:
        tr._branch_cmax = float(s.cmax)


def _walk_back(walks, c, e, run, mono):
    """The walk as the last step made leaves it: the step's reduced J, Geq
    and C and its coefficients on the walk object (`solve_timestep`'s), the
    monodromy one step on (its live columns, its step count) -- and the
    capacitance ring returned.  Where the last step's map left a column
    dead, the state one step on is `_LiveMonodromy.advance`'s (the prune's
    views); where its map handed back, the map is made as `_walk_lmm` makes
    it (`_lmmsens_c.step`, then `_LiveMonodromy.step`)."""
    from pycircuit.circuit.analysis import remove_row_col
    s, m, tr, tk = c.run.s, e.m, e.tr, walks.toolkit
    w = int(s.w)

    def mat(k):
        return c.buf[k].reshape(m, m).copy()

    def pick(k):
        return c.buf[k][:m * w].reshape(m, w).copy()
    J = e.rec.buf['Jj'].reshape(e.n, e.n)
    (walks._Jf, walks._Geq, walks._C) = remove_row_col((J, tr._Geq, tr._Cmat), e.iref, tk)
    walks._coeffs = tr._companion_coeffs
    p0, p1, pq, cn, cp = int(s.p0), int(s.p1), int(s.pq), int(s.cn), int(s.cp)
    Cs = [mat(f'Cr{cn}'), mat(f'Cr{cp}')]
    if run.status == _PRUNE:
        ## (the state the step ran from -- P_n, P_{n-1} in the rings' third
        ## place, Pq_n in its other -- moved on by the Python)
        mono.P = [pick(f'Pb{p1}'), pick(f'Pb{3 - p0 - p1}')]
        mono.Pq, mono.n = pick(f'Pqb{1 - pq}'), int(s.mono_n) - 1
        mono.advance(pick(f'Pb{p0}'), pick(f'Pqb{pq}'))
        return Cs
    mono.P = [pick(f'Pb{p0}'), pick(f'Pb{p1}')]
    mono.Pq, mono.n = pick(f'Pqb{pq}'), int(s.mono_n)
    if run.status != _MAP:
        return Cs
    from . import _lmmsens_c
    C_new = mat(f'Cr{3 - cn - cp}')
    alphas, b = walks._coeffs
    Jf = np.asarray(walks._Jf)
    solver = walks._get_linearsolver()

    def solve(S, _Jf=Jf):
        return solver.solve(_Jf, S, tk)
    got = _lmmsens_c.step(walks, mono, Cs, C_new, alphas, b, Jf)
    if got is None:
        mono.step(Cs, C_new, alphas, b, solve)
    else:
        mono.advance(*got)
    return [C_new, Cs[0]]
