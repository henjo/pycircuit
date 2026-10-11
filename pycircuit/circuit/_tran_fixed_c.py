"""THE FIXED-STEP MULTISTEP LOOP IN C (the C stepping loop plan's stage 4,
step b; 2026-10-11).  `_SteppingLoop` under `fixed_timestep` with a multistep
family is, step after step, `where` (the grid, the next breakpoint, the
order drop a landing asks for), `take` (the stage predictor's seed, the
companion, the C Newton, the branch screen) and `accept` (the rings, the
predictor's node, the outputs, the statistics) -- ~85 us of Python around
~20 us of C a step on the 20-MosLevel1 chain.  Here `csrc/fixed.c` runs those
steps itself, driving the units the Python steps call -- the sources' table
(stage 1), the predictor (4a), the C Newton (its record, as its last served
call left it), the branch screen (3d) -- and stops before any step it
cannot take as the Python would, which the Python loop then takes; the
next step enters the C again.

THE RUN CONTEXT.  Inside the C loop no Python runs between steps: the
sources are the table, the devices the core's kernels, and a circuit with
an `accept_step` hook, a breakpoint the table does not make, a periodic
state, a provided function, a stateful limiter, an injected piece of the
chain declines.  So what the C Newton's glue checks on every call holds
for every step of an entry when its last call SERVED and nothing watched
moved since (`_tran_newton_c`'s `served` mark against `_watch.now()`); the
rest is checked here, once an entry.

THE STATE.  Copied in at entry (the rings, the predictor's nodes, the
cached C, the step record, the loop's `t`, `h`, `landing`) and written back
after it, as the Python steps leave it: `X` and `timelist` a fresh array
and float a step, the statistics, and every attribute the steps write
(`transient._RUN_STATE`: the rings, `_q_cache`, `_C_cache`, `_C_pattern`,
`_iq`, `_Geq`, `_Cmat`, the active integrator and its coefficients, the
predictor's history, the step record, the screen's running maximum) --
and the `_paths` tallies the Python steps make, under their keys.
`ENABLED` (env `PYCIRCUIT_FIXED_C=0`) is read on every entry.
"""
import os
import time

import numpy as np

from pycircuit.circuit import _cbuild, _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_FIXED_C', '1') != '0'
STATUS = 'not loaded'

FIXED_C = _cbuild.csrc('fixed.c')
FIXED_CDEF = _cbuild.csrc('fixed.h')
_CODES = _cbuild.constants(FIXED_CDEF)
BAIL = _cbuild.bails(_CODES, 'FIXED_BAIL_')
_ABI = _cbuild.constants(_cbuild.csrc('pcx_abi.h'))
#: the outputs' rows an entry holds (full: written back, entered again)
CAP_MAX = 4096
_F64 = np.dtype(np.float64)

_driver = None
_OWN = None


def driver():
    """`(ffi, cfn)` once loaded, or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(FIXED_C, FIXED_CDEF)
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn)
        STATUS = 'c'
    return _driver or None


def _own():
    """The pieces the loop stands in for, as defined -- read once every one
    is its module's own (`_paths.genuine`), else None -- with what an entry
    compares against them (`_now`)."""
    global _OWN
    if _OWN is None:
        from pycircuit.circuit import (
            _tran_companion,
            _tran_history,
            _tran_newton_c,
            _tran_predictor,
            hdl,
        )
        from pycircuit.circuit import elements as E
        from pycircuit.circuit import integrator as I
        from pycircuit.circuit.circuit import Circuit, SubCircuit
        from pycircuit.circuit.transient import (
            Transient,
            TransientStatistics,
            _StepFamily,
        )
        M = 'pycircuit.circuit.'
        TR, PR, HI, BR, CO, EV, CI, EL = (M + m for m in (
            'transient', '_tran_predictor', '_tran_history', '_tran_branch', '_tran_companion',
            '_tran_events', 'circuit', 'elements'))
        T = Transient
        own = (
            (T._next_breakpoint, '_StepEvents._next_breakpoint', EV),
            (T._attempt_step, 'Transient._attempt_step', TR),
            (T._pred_or, '_StagePredictor._pred_or', PR),
            (T._predict_state, '_StagePredictor._predict_state', PR),
            (T._predict_fast, '_StagePredictor._predict_fast', PR),
            (T._pred_fast_nodes, '_StagePredictor._pred_fast_nodes', PR),
            (T._fit_fast, '_StagePredictor._fit_fast', PR),
            (T._pred_note, '_StagePredictor._pred_note', PR),
            (T._pred_promote, '_StagePredictor._pred_promote', PR),
            (T._pred_degree, '_StagePredictor._pred_degree', PR),
            (T._push_history, '_RunHistory._push_history', HI),
            (T._roll_history, '_RunHistory._roll_history', HI),
            (_tran_history._ring_push, '_ring_push', HI),
            (T._q_at, '_CompanionModel._q_at', CO),
            (T._branch_after_solve, '_BranchCheck._branch_after_solve', BR),
            (T._branch_screen, '_BranchCheck._branch_screen', BR),
            (T._branch_screen_C, '_BranchCheck._branch_screen_C', BR),
            (T._branch_on, '_BranchCheck._branch_on', BR),
            (TransientStatistics._note_step, 'TransientStatistics._note_step', TR),
            (SubCircuit.next_event, 'SubCircuit.next_event', CI),
            (SubCircuit.accept_step, 'SubCircuit.accept_step', CI),
            (E.VS.next_event, 'VS.next_event', EL),
            (E.IS.next_event, 'IS.next_event', EL),
            (hdl.Behavioural.next_event, 'Behavioural.next_event', M + 'hdl'),
            (T._source_at, '_CompanionModel._source_at', CO),
            (SubCircuit.u, 'SubCircuit.u', CI),
            (T._branch_structural_rank, '_BranchCheck._branch_structural_rank', BR),
            (_StepFamily.next_breakpoint, '_StepFamily.next_breakpoint', TR),
            (_StepFamily.after_accept, '_StepFamily.after_accept', TR),
            (_tran_newton_c.solve, 'solve', M + '_tran_newton_c'),
            (I.Gear2Integrator.check_order_drop, 'Gear2Integrator.check_order_drop', M + 'integrator'),
            (I.EulerIntegrator.check_order_drop, 'EulerIntegrator.check_order_drop', M + 'integrator'),
            (I.TrapezoidalIntegrator.check_order_drop, 'TrapezoidalIntegrator.check_order_drop',
             M + 'integrator'),
        )
        if not all(_paths.genuine(*o) for o in own):
            return None
        _OWN = {
            'chain': tuple(o for o, _q, _m in own),
            'T': T, 'stats': TransientStatistics, 'Circuit': Circuit, 'Sub': SubCircuit,
            'history': _tran_history, 'E': E, 'hdl': hdl,
            'source_next': (E.VS.next_event, E.IS.next_event),
            'inf_next': hdl.Behavioural.next_event,
            'integrators': {I.Gear2Integrator: _ABI['PCX_FORMULA_GEAR2'],
                            I.EulerIntegrator: _ABI['PCX_FORMULA_EULER'],
                            I.TrapezoidalIntegrator: _ABI['PCX_FORMULA_TRAPEZOIDAL']},
            'Euler': I.EulerIntegrator,
            'zs_ratio': I.ZERO_STABILITY_RATIO,
            'predictor': _tran_predictor, 'companion': _tran_companion,
            'newton_c': _tran_newton_c, 'family': _StepFamily, 'I': I,
        }
    return _OWN


def _now(O):
    """The pieces `_own` read, as they stand now."""
    T, E, Sub, F, I = O['T'], O['E'], O['Sub'], O['family'], O['I']
    return (T._next_breakpoint, T._attempt_step, T._pred_or, T._predict_state, T._predict_fast,
            T._pred_fast_nodes, T._fit_fast, T._pred_note, T._pred_promote, T._pred_degree,
            T._push_history, T._roll_history, O['history']._ring_push, T._q_at,
            T._branch_after_solve, T._branch_screen, T._branch_screen_C, T._branch_on,
            O['stats']._note_step, Sub.next_event, Sub.accept_step, E.VS.next_event,
            E.IS.next_event, O['hdl'].Behavioural.next_event, T._source_at, Sub.u,
            T._branch_structural_rank, F.next_breakpoint, F.after_accept, O['newton_c'].solve,
            I.Gear2Integrator.check_order_drop, I.EulerIntegrator.check_order_drop,
            I.TrapezoidalIntegrator.check_order_drop)


class _Run:
    """One stepping loop's C loop: the struct, its buffers, the units'
    structs it drives -- bound to the C Newton's record, the source table and
    the screen it was made for."""

    __slots__ = ('br', 'buf', 'ffi', 'keep', 'npat', 'nrec', 'pr', 's', 'tab')


def _addr(ffi, obj):
    return int(ffi.cast('uintptr_t', obj))


def _bind(tr, rec, tab, n, hist_len, screen):
    """The loop's `_Run` for the C Newton's record `rec`, the source table
    `tab` and the branch screen `screen` (`(its ffi, its entry, the
    core's pattern or None)`, or None: nothing screened)."""
    from pycircuit.circuit import _tran_newton_c, _tran_pred_c
    ffi, _cfn = driver()
    run = _Run()
    run.ffi, run.nrec, run.tab = ffi, rec, tab
    s = run.s = ffi.new('fixed_t *')
    cap = CAP_MAX
    run.buf = buf = {
        'qring': np.empty(hist_len * n), 'iqring': np.empty(hist_len * n),
        'xlast': np.empty(n), 'ht': np.empty(8), 'hx': np.empty(8 * n),
        ## (zero off the pattern: where the core writes only the pattern's
        ## entries, the rest stay as the Newton's own buffers keep them)
        'cC0': np.zeros(n * n), 'cC1': np.zeros(n * n),
        'Geq0': np.zeros(n * n), 'Geq1': np.zeros(n * n),
        'X': np.empty(cap * n), 'T': np.empty(cap),
    }
    keep = run.keep = []
    fb = ffi.from_buffer
    for k in ('qring', 'iqring', 'xlast', 'ht', 'hx', 'X', 'T'):
        keep.append(fb('double *', buf[k]))
        setattr(s, k, keep[-1])
    for i in (0, 1):
        keep.append(fb('double *', buf[f'cC{i}']))
        s.cC[i] = keep[-1]
        keep.append(fb('double *', buf[f'Geq{i}']))
        s.Geq[i] = keep[-1]
    nffi, ncfn, _lapack = _tran_newton_c.driver()
    s.newton = ffi.cast('void *', _addr(nffi, rec.s))
    s.newton_fn = ffi.cast('void *', _addr(nffi, ncfn))
    pffi, pcfn, dgesv, dgemv = _tran_pred_c.driver()
    ps, pbuf, pkeep = _tran_pred_c.new_state(pffi, dgesv, dgemv, n)
    ps.ht = pffi.cast('double *', _addr(ffi, s.ht))
    ps.hx = pffi.cast('double *', _addr(ffi, s.hx))
    ps.deg = int(tr._pred_degree())
    ps.clamp = float(tr.PRED_CLAMP)
    run.pr = (ps, pbuf, pkeep)
    s.pred = ffi.cast('void *', _addr(pffi, ps))
    s.pred_fn = ffi.cast('void *', _addr(pffi, pcfn))
    s.src = ffi.cast('void *', _addr(tab.ffi, tab.s))
    s.src_fn = ffi.cast('void *', _addr(tab.ffi, tab.cfn))
    run.br, run.npat = None, 0
    if screen is None:
        s.branch = s.branch_fn = ffi.NULL
    else:
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
        s.branch = ffi.cast('void *', _addr(bffi, bs))
        s.branch_fn = ffi.cast('void *', _addr(bffi, bcfn))
    s.n, s.hist_len, s.cap = n, hist_len, cap
    return run


def _declined(loop, why, keep=True):
    """The entry declined, counted -- kept on the loop (`keep`) until
    something watched moves, where the reason is the run's and not the
    step's."""
    from pycircuit.circuit import _watch
    if keep:
        loop.c_fixed = _watch.now()
    _no('fixed_c:' + why)
    return False


def _screen(tr, td, rec, n):
    """The branch screen the steps make: `(None,)` -- nothing screened (off,
    or `C` structurally empty) -- `(ffi, entry, pattern or None)` in C, or
    None: not in C."""
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


def run(loop):
    """Take as many of `loop`'s fixed steps in C as it can, from where the
    loop stands: True where it took any (the loop's and the transient's
    state written back as the Python steps leave it), else False -- the
    Python loop takes the next step, and the one after enters again."""
    if not ENABLED:
        return False
    from pycircuit.circuit import _watch
    now = _watch.now()
    ## (watching off: nothing would tell the steps' checks have moved --
    ## the counter stays -1 -- so the Python steps make them, every step)
    if now == -1 or loop.c_fixed == now:
        return False
    tr = loop.tr
    td = tr.__dict__
    rec = td.get('_newton_c')
    if rec is None or rec.served != now:
        ## (the C Newton's last call did not serve, or something watched
        ## moved since: the Python step, which serves or not by its checks)
        return False
    drv = _driver or driver()
    if not drv:
        return _declined(loop, 'driver')
    O = _OWN or _own()
    if O is None or type(tr) is not O['T'] or _now(O) != O['chain']:
        return _declined(loop, 'patched')
    if type(tr.statistics) is not O['stats']:
        return _declined(loop, 'stats')
    if not loop.j_unread or loop.provided_function is not None:
        return _declined(loop, 'loop')
    if loop.imposed or loop.order_drop:
        return _declined(loop, 'loop', keep=False)
    for k in _SHADOWS:
        if k in td:
            return _declined(loop, 'shadow')
    if getattr(tr, '_ev_rows', None) is not None or getattr(tr, '_periodic_rows', None):
        return _declined(loop, 'events')
    if getattr(tr, '_continuation_rescue', False) or td.get('_stateful_lims'):
        return _declined(loop, 'state')
    cir = tr.cir
    cd = cir.__dict__
    if (type(cir).next_event is not O['Sub'].next_event
            or type(cir).accept_step is not O['Sub'].accept_step
            or type(cir).u is not O['Sub'].u or any(k in cd for k in _CIRCUIT_SHADOWS)
            or cir._hook_elements('accept_step')):
        return _declined(loop, 'circuit')
    base = tr.base_integrator
    code = O['integrators'].get(type(base))
    if code is None or 'check_order_drop' in base.__dict__:
        return _declined(loop, 'integrator')
    P = O['predictor']
    if (not P.PRED_FAST or tr.stage_predictor != 'on' or tr.PRED_HIST_MAX != 8
            or type(tr)._fit is not P._FIT or type(tr.PRED_CLAMP) is not float
            or not 1 <= tr._pred_degree() <= 4):
        return _declined(loop, 'predictor')
    if not O['companion'].U_MEMO:
        return _declined(loop, 'umemo')
    ## the sources: the direct pass's table, every other breakpoint inf
    from pycircuit.circuit import _tran_sources, func
    plan = cir.__dict__.get('_stamp_plan')
    sp = plan.methods.get('src:u') if plan is not None else None
    sd = sp.direct if sp is not None else None
    ep = _watch.EPOCH
    if sd is None or ep is None or sd.stamp != ep.value or not sd.same(cir, O['Circuit']):
        return _declined(loop, 'sources')
    ct = loop.c_tab
    if ct is not None and ct[0] is sd:
        tab = ct[1]
    else:
        tab = _tran_sources.table_for(sd)
        if tab is None:
            return _declined(loop, 'table')
        loop.c_tab = (sd, tab)
    if not tab.holds():
        loop.c_tab = None
        return _declined(loop, 'table', keep=False)
    tab_els = {id(e[0]) for e in tab.entries}
    src_next, inf_next = O['source_next'], O['inf_next']
    for _inst, el in cir._hook_elements('next_event'):
        f = type(el).next_event
        if 'next_event' in el.__dict__ or not (
                (f in src_next and id(el) in tab_els) or f is inf_next):
            return _declined(loop, 'breakpoints')
    tk = tab.tk
    if not (tk.exp is np.exp and tk.sin is np.sin and tk.where is np.where) or (
            tab.uses_libm and ((tab.has_sin and not func.SIN_SCALAR)
                               or not func._underflow_ignored())):
        return _declined(loop, 'sources', keep=False)
    minbreak = tr.par.minbreak
    if type(minbreak) is not float:
        return _declined(loop, 'minbreak')
    ## the state the steps read
    n = rec.n
    ql, iql = tr._qlast, tr._iqlast
    if not (type(ql) is np.ndarray and ql.dtype is _F64 and ql.ndim == 2 and ql.shape[1] == n
            and ql.shape[0] >= 2 and type(iql) is np.ndarray and iql.dtype is _F64
            and iql.shape == ql.shape):
        return _declined(loop, 'rings', keep=False)
    hist_len = ql.shape[0]
    hist = td.get('_pred_hist') or []
    if type(hist) is not list or len(hist) > 8 or td.get('_pred_pending') is not None or not all(
            type(e) is tuple and len(e) == 2 and type(e[0]) is float and type(e[1]) is np.ndarray
            and e[1].dtype is _F64 and e[1].shape == (n,) for e in hist):
        return _declined(loop, 'history', keep=False)
    X = loop.X
    xl = X[-1] if X else None
    cached = td.get('_C_cache')
    if cached is None or float(getattr(tr.epar, 'bypasstol', -1.0) or -1.0) >= 0.0:
        return _declined(loop, 'cache', keep=False)
    cx, cC = cached
    if not (type(xl) is np.ndarray and xl.dtype is _F64 and xl.shape == (n,)
            and type(cx) is np.ndarray and cx.tobytes() == xl.tobytes()
            and type(cC) is np.ndarray and cC.dtype is _F64 and cC.shape == (n, n)):
        return _declined(loop, 'cache', keep=False)
    dl, dl2 = tr._dt_last, tr._dt_last2
    stats = tr.statistics
    if (type(dl) is not float or not (dl2 is None or type(dl2) is float)
            or not (stats.min_step is None or type(stats.min_step) is float)
            or not (stats.max_step is None or type(stats.max_step) is float)
            or type(loop.t) is not float or type(loop.h) is not float):
        return _declined(loop, 'steps', keep=False)
    screen = _screen(tr, td, rec, n)
    if screen is None:
        return _declined(loop, 'screen')
    if screen[0] is None:
        screen = None
    ## bound once a loop (again where the Newton's record or the table moved)
    run_ = loop.c_run
    if run_ is None or run_.nrec is not rec or run_.tab is not tab or \
            (run_.br is None) != (screen is None) or len(run_.buf['xlast']) != n or \
            len(run_.buf['qring']) != hist_len * n:
        run_ = loop.c_run = _bind(tr, rec, tab, n, hist_len, screen)
    s, buf = run_.s, run_.buf
    ## in: the state as the Python loop leaves it
    buf['qring'][:] = ql.reshape(-1)
    buf['iqring'][:] = iql.reshape(-1)
    buf['xlast'][:] = xl
    s.hk = len(hist)
    for i, (tt, xx) in enumerate(hist):
        buf['ht'][i] = tt
        buf['hx'][i * n:(i + 1) * n] = xx
    buf['cC0'][:] = cC.reshape(-1)
    s.cur = 0
    s.cmax = td.get('_branch_cmax', 0.0)
    s.t, s.h = loop.t, loop.h
    s.landing, s.is_first = int(bool(loop.landing)), int(bool(tr._is_first_step))
    s.dt_last, s.have_dt_last = dl, 1
    s.dt_last2, s.have_dt_last2 = (dl2, 1) if dl2 is not None else (0.0, 0)
    s.tend, s.timestep, s.minbreak = float(loop.tend), float(loop.timestep), minbreak
    s.base = code
    s.zs_ratio = float(O['zs_ratio'])
    s.shrink_guard = int(bool(getattr(base, 'shrink_guard', True)))
    s.have_hm = int(stats.min_step is not None)
    s.hmin = stats.min_step if stats.min_step is not None else 0.0
    s.hmax = stats.max_step if stats.max_step is not None else 0.0
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
    k = s.k
    if status not in (_CODES['FIXED_DONE'], _CODES['FIXED_FULL']):
        _no('fixed_c:bail:' + BAIL.get(status, str(status)))
    if not k:
        return False
    _write_back(loop, tr, rec, run_, s, buf, n, hist_len, k, O)
    return True


#: the circuit's methods the steps call, shadowed on the instance (a
#: perturbed source, as the shooting tests make): the Python steps call the
#: shadow, the C the table and the core
_CIRCUIT_SHADOWS = ('u', 'G', 'C', 'i', 'q', 'limit', 'next_event', 'accept_step')
#: the transient's methods the steps call, shadowed on the instance: the
#: Python steps call the shadow
_SHADOWS = ('_next_breakpoint', '_attempt_step', 'solve_timestep', '_solve_timestep', '_pred_or',
            '_predict_state', '_predict_fast', '_fit_fast', '_fit', '_pred_note', '_pred_promote',
            '_push_history', '_roll_history', '_q_at', '_branch_after_solve', '_branch_screen',
            '_branch_screen_C', '_branch_on', '_pred_degree', '_branch_structural_rank',
            '_source_at')


def _write_back(loop, tr, rec, run_, s, buf, n, hist_len, k, O):
    """The state after the C's `k` steps, as the Python steps leave it."""
    Xo, To = buf['X'], buf['T']
    X, tl = loop.X, loop.timelist
    for i in range(k):
        X.append(Xo[i * n:(i + 1) * n].copy())
        tl.append(float(To[i]))
    stats = tr.statistics
    stats.accepted_steps += k
    stats.newton_iterations += s.iters
    stats.breakpoints_hit += s.breaks
    stats.order_drops += s.euler_steps
    stats.min_step, stats.max_step = float(s.hmin), float(s.hmax)
    ## the rings and the converged point's state
    tr._qlast = buf['qring'].reshape(hist_len, n).copy()
    tr._iqlast = buf['iqring'].reshape(hist_len, n).copy()
    x = buf['xlast'].copy()
    cur = s.cur
    Cj = buf[f'cC{cur}'].reshape(n, n).copy()
    tr._q_cache = (x, buf['qring'][:n].copy())
    tr._C_cache = (x, Cj)
    if rec.core.npat:
        tr._C_pattern = (Cj, rec.core.arrays['pat'])
    tr._iq = buf['iqring'][:n].copy()
    tr._Geq = buf[f'Geq{cur}'].reshape(n, n).copy()
    tr._Cmat = Cj
    h, h_last = float(s.dt_last), float(s.last_h_last)
    base = tr.base_integrator
    active = (O['Euler']() if s.last_formula == _ABI['PCX_FORMULA_EULER']
              and type(base) is not O['Euler'] else base)
    tr.active_integrator = active
    tr._companion_coeffs = active.companion_coefficients(h, h_last)
    tr._effective_method = type(active).__name__
    tr._pred_hist = [(float(buf['ht'][i]), buf['hx'][i * n:(i + 1) * n].copy())
                     for i in range(s.hk)]
    tr._pred_pending = None
    tr._dt, tr._dt_last = h, h
    tr._dt_last2 = float(s.dt_last2) if s.have_dt_last2 else None
    tr._is_first_step = bool(s.is_first)
    tr._no_history = False
    if run_.br is not None:
        tr._branch_cmax = float(s.cmax)
    loop.t, loop.h = float(s.t), float(s.h)
    loop.landing = bool(s.landing)
    loop.t_break = s.t_break
    loop.x_new, loop.J = x, None
    loop.h_next = loop.timestep
    loop.clamped, loop.hold = False, True
    ## the tallies the Python steps make
    _PC['fixed_c:steps'] += k
    if rec.klu is not None:
        ## (each step's numeric `StructuralKLU.factor` would have made)
        rec.klu.numerics += k
        _PC['newton_c:sparse'] += k
    _PC['newton_c:served'] += k
    _PC['newton_c:iterations'] += s.iters
    _PC['newton_c:jfold'] += k
    _PC['newton_c:j_unread'] += k
    if rec.s.chord:
        _PC['chord:chosen'] += k
    if s.pred_fast:
        _PC['pred:fast'] += s.pred_fast
    _PC['src.direct:served'] += k
    _PC['u:called'] += 2 * k
    _PC['umemo:miss'] += k
    _PC['memo.q:same'] += k
    if s.screens:
        _PC['branch_c:served'] += s.screens
        _PC['memo.C:equal'] += s.screens
        if run_.npat:
            _PC['branch_c:pattern'] += s.screens
