"""The shooting <-> transient boundary, measured (2026-09-28).

Phase V of the interface plan: each difference between the transient's own
stepping loop and the PSS's use of `Transient` that no comment explains is
MEASURED here before it is called a defect or refactored away.  One function
per item, its prediction in its docstring, written before the run.

    python benchmarks/pss_transient_boundary.py [V6 V9 V3 ...]

Nothing here changes the code: the instruments wrap `Transient` methods for
the duration of one measurement and restore them.
"""
import sys
import warnings
from contextlib import contextmanager

import numpy as np

from pycircuit.circuit import circuit
from pycircuit.circuit.elements import C, R, VSin, SubCircuit, gnd
from pycircuit.circuit.shooting import PSS
from pycircuit.circuit.transient import Transient

circuit.default_toolkit = circuit.numeric


## ---------------------------------------------------------------------------
## fixtures
## ---------------------------------------------------------------------------

def cv_loop(per=1e-3):
    """The driven index-2 C-V loop of `test_glm.py`."""
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
    c['c1'] = C('a', 'b', c=1e-9)
    c['c2'] = C('b', gnd, c=1e-9)
    c['r'] = R('b', gnd, r=1e5)
    return c


def jittered_fracs(npts, amp=0.25):
    """Step fractions whose size varies ~3x along the period (`test_glm.py`
    `_grid`): the GLM rescales its Nordsieck vector at every step."""
    u = np.linspace(0.0, 1.0, npts + 1)
    u = u + amp * np.sin(2 * np.pi * u) / np.pi
    u = (u - u[0]) / (u[-1] - u[0])
    return list(np.diff(u))


def glm_cases():
    """(label, runner): GLM runs through every path that reaches the GLM
    step -- forward adaptive twice on one object, PSS on a uniform and a
    jittered grid, a staged oscillator (event-landed grid, restarts on
    growth), the factored period and the PAC read."""
    from pycircuit.circuit.integrator import GLM3Integrator
    from pycircuit.circuit.tests._shooting_fixtures import (
        _comparator_relaxation_oscillator, _relaxation_oscillator_seed)

    def forward_twice():
        c = cv_loop()
        tr = Transient(c, integrator=GLM3Integrator(), reltol=1e-6)
        tr.solve(tend=2e-3, timestep=2e-5, x0=np.zeros(c.n))
        tr.solve(tend=1.3e-3, timestep=2e-5, x0=np.zeros(c.n))

    def pss_uniform(method):
        def run():
            p = PSS(cv_loop(), method=method, reltol=1e-12)
            p.solve(period=1e-3, timestep=1e-3 / 40, maxiterations=40)
            assert p.converged
            p.factored_period()
        return run

    def pss_jittered(method):
        def run():
            p = PSS(cv_loop(), method=method, reltol=1e-12)
            p.solve(period=1e-3, timestep=1e-3 / 40, grid=jittered_fracs(40),
                    maxiterations=40)
            assert p.converged
            p.factored_period()
        return run

    def staged(method):
        def run():
            cir = _comparator_relaxation_oscillator()
            seed, Tl = _relaxation_oscillator_seed(cir)
            q = PSS(_comparator_relaxation_oscillator(), method=method,
                    reltol=1e-9)
            q.solve(period=Tl, timestep=Tl / 200, x0=seed, maxiterations=100,
                    state_events=True)
            assert q.converged
            q.ppv()
        return run

    return [('forward glm3, solved twice', forward_twice),
            ('PSS glm2 uniform', pss_uniform('glm2')),
            ('PSS glm3 uniform', pss_uniform('glm3')),
            ('PSS glm3 jittered', pss_jittered('glm3')),
            ('PSS glm4 jittered', pss_jittered('glm4')),
            ('PSS glm3 staged oscillator', staged('glm3')),
            ('PSS glm2 staged oscillator', staged('glm2'))]


@contextmanager
def patched(cls, name, make):
    """`cls.name` replaced by `make(original)` for the duration."""
    orig = getattr(cls, name)
    setattr(cls, name, make(orig))
    try:
        yield orig
    finally:
        setattr(cls, name, orig)


## ---------------------------------------------------------------------------
## V6 -- `_begin_run` does not reset `_glm_prev`
## ---------------------------------------------------------------------------

def V6():
    """`_begin_run` clears `_glm_Q` and `_glm_Q_at_entry` but not `_glm_prev`,
    the GLM stage predictor's record of the last step.  The predictor
    accepts a record only when its end time equals the step's start and its
    step equals the step's (`_glm_stage_predictor`).

    PREDICTED BENIGN: 0 accepts of a record written before the current run's
    `_begin_run`.  Every run starts at t = 0 (a PSS walk at `times[0] = 0`)
    and a stale record ends at t >= h > 0, so the first step declines it and
    overwrites it.  A count above 0 is a defect: a prediction carried across
    a run boundary."""
    stats = {'accepts': 0, 'stale_accepts': 0, 'steps': 0}

    def begin(orig):
        def f(self, *a, **k):
            self._v6_run = getattr(self, '_v6_run', 0) + 1
            return orig(self, *a, **k)
        return f

    def step(orig):
        def f(self, *a, **k):
            out = orig(self, *a, **k)
            stats['steps'] += 1
            self._v6_prev_run = getattr(self, '_v6_run', 0)
            return out
        return f

    def pred(orig):
        def f(self, h, tn, c):
            prev = getattr(self, '_glm_prev', None)
            stale = (prev is not None and getattr(self, '_v6_prev_run', None)
                     != getattr(self, '_v6_run', 0))
            out = orig(self, h, tn, c)
            if out is not None:
                stats['accepts'] += 1
                stats['stale_accepts'] += int(stale)
            return out
        return f

    rows = []
    with patched(Transient, '_begin_run', begin), \
            patched(Transient, '_solve_timestep_glm', step), \
            patched(Transient, '_glm_stage_predictor', pred):
        for label, run in glm_cases():
            before = dict(stats)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                run()
            rows.append((label, stats['steps'] - before['steps'],
                         stats['accepts'] - before['accepts'],
                         stats['stale_accepts'] - before['stale_accepts']))
    print('V6  GLM predictor accepts of a record from an EARLIER run '
          '(predicted 0):')
    for label, n, acc, stale in rows:
        print('    %-30s steps %5d  predictor accepts %5d  stale %d'
              % (label, n, acc, stale))
    return sum(r[3] for r in rows)


## ---------------------------------------------------------------------------
## V9 -- a time-key restart does not set `_glm_restarted`
## ---------------------------------------------------------------------------

def V9():
    """`_solve_timestep_glm` starts afresh (`_glm_startup`) when neither slot
    matches the step's start time, but flags `_glm_restarted` only for a
    GROWTH restart.  The shooting reads the flag to decide how a step's
    output enters the map (`_glm_period_blocks`), so an unflagged fresh
    start inside a period would be linearised as a continuation.

    PREDICTED BENIGN: inside a PSS walk every fresh start after the period's
    first step is flagged -- the grids hand `solve_timestep` times that match
    the stored slot's to rounding (`_period_grid`), so a time-key miss never
    happens.  Counted: steps with t_start > 0 that started afresh, and of
    those the unflagged."""
    stats = {'in_pss': False}
    rec = []

    def walk(orig):
        def f(self, *a, **k):
            stats['in_pss'] = True
            try:
                return orig(self, *a, **k)
            finally:
                stats['in_pss'] = False
        return f

    def startup(orig):
        def f(self, *a, **k):
            if getattr(self, '_v9_in_step', False):
                self._v9_started = True
            return orig(self, *a, **k)
        return f

    def step(orig):
        def f(self, x0, t, provided_function=None):
            self._v9_in_step, self._v9_started = True, False
            try:
                out = orig(self, x0, t, provided_function)
            finally:
                self._v9_in_step = False
            tn = float(t) - float(self._dt)
            rec.append((stats['in_pss'], tn, self._v9_started,
                        bool(getattr(self, '_glm_restarted', False))))
            return out
        return f

    from pycircuit.circuit.shooting._pss_walks import _PeriodWalks
    rows = []
    with patched(_PeriodWalks, '_glm_period_blocks', walk), \
            patched(Transient, '_glm_startup', startup), \
            patched(Transient, '_solve_timestep_glm', step):
        for label, run in glm_cases():
            n0 = len(rec)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                run()
            new = rec[n0:]
            pss = [r for r in new if r[0]]
            mid = [r for r in pss if r[1] > 0.0 and r[2]]
            rows.append((label, len(pss), len(mid),
                         sum(1 for r in mid if not r[3]),
                         sum(1 for r in pss if r[3])))
    print('V9  fresh GLM starts after a period\'s first step, and the '
          'unflagged ones (predicted 0 unflagged):')
    for label, n, mid, unflagged, flagged in rows:
        print('    %-30s PSS-walk steps %5d  fresh starts at t > 0 %4d  '
              'unflagged %d  (flagged restarts %d)'
              % (label, n, mid, unflagged, flagged))
    return sum(r[3] for r in rows)


## ---------------------------------------------------------------------------
## V3 -- `_history_is_solved` is never cleared per period
## ---------------------------------------------------------------------------

def V3():
    """`_install_history` (the pair map's opening) sets
    `pss._history_is_solved = True` and nothing clears it, so a PLAIN walk
    with the LTE measured, on a PSS object that once ran the pair map,
    computes `_lte_seam` with it stale.

    PREDICTED LATENT: (1) structurally, the pair map is exactly the methods
    with `companion_reach() >= 2` (Gear-2), and every Gear-2 walk opens
    through `_install_history`, so no production path runs a plain walk
    after a pair one on one object; (2) forced (the plain-walk pattern of
    `test_shooting_methods.py`), the stale object reports `_lte_seam` False
    where a fresh one reports True on the same steps."""
    from pycircuit.circuit.elements import Diode
    names = ['euler', 'trap', 'theta', 'gear', 'trbdf2', 'radau',
             'esdirk43', 'glm2', 'glm3', 'glm4']
    kinds = {}
    for name in names:
        p = PSS(cv_loop(), method=name)
        integ = p._integrator_for(name)
        reach = (None if integ.is_stage_method()
                 else int(integ.companion_reach()))
        kinds[name] = (p._map_kind(), reach)
    print('V3  map kind and companion reach per method:')
    for name, (kind, reach) in kinds.items():
        print('    %-9s %-6s reach %s' % (name, kind, reach))
    pair = sorted(n for n, (k, _r) in kinds.items() if k == 'pair')
    reach2 = sorted(n for n, (_k, r) in kinds.items() if r is not None and r >= 2)
    print('    pair map: %s   reach >= 2: %s   same: %s'
          % (pair, reach2, pair == reach2))

    def build():
        c = SubCircuit()
        c['vs'] = VSin(1, gnd, vac=1.0, va=0.5, freq=1e6, phase=20)
        c['R'] = R(1, 2, r=1e3)
        c['Rd'] = R(2, 3, r=1e3)
        c['D'] = Diode(3, gnd)
        c['C'] = C(2, gnd, c=1e-9)
        return c

    def plain_walk_seams(pss):
        T = 1e-6
        n = 20
        times = np.linspace(0.0, T, n + 1)
        hs = np.diff(times)
        m = pss.cir.n - 1
        tr_saved = getattr(pss, '_tran', None)
        pss._tran = pss._new_transient(pss._integrator_for('gear'))
        seams = []
        try:
            pss._want_dfdh = False
            pss._want_lte = True
            pss._want_event_cols = False
            x = np.zeros(m)
            pss._begin_period(x)
            for j, t in enumerate(times[1:]):
                x = np.asarray(pss.solve_timestep(x, t, hs[j]), float)
                seams.append(bool(getattr(pss, '_lte_seam', False)))
        finally:
            pss._tran = tr_saved
            pss._want_lte = False
        return seams

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        used = PSS(build(), method='gear', reltol=1e-9)
        used.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=40)
        flag_after = bool(getattr(used, '_history_is_solved', False))
        stale = plain_walk_seams(used)
        fresh = plain_walk_seams(PSS(build(), method='gear', reltol=1e-9))
    print('    after a Gear-2 solve `_history_is_solved` = %s' % flag_after)
    print('    forced plain walk, `_lte_seam` per step: stale %s' % stale[:5])
    print('                                            fresh %s' % fresh[:5])
    differs = stale != fresh
    print('    the stale flag %s the plain walk\'s seam report'
          % ('CHANGES' if differs else 'does not change'))
    return int(pair != reach2), int(differs)



## ---------------------------------------------------------------------------
## V12 -- the PSS's `epar` does not reach its inner transient
## ---------------------------------------------------------------------------

def diode_rc():
    """A driven diode into an RC: the diode's thermal voltage reads
    `epar.T`, so the temperature moves the waveform."""
    from pycircuit.circuit.elements import Diode
    c = SubCircuit()
    c['vs'] = VSin(1, gnd, va=1.0, freq=1e3)
    c['R'] = R(1, 2, r=1e3)
    c['D'] = Diode(2, gnd)
    c['C'] = C(2, gnd, c=1e-7)
    return c


def V12():
    """`_new_transient` passes no `epar`, and `Analysis.__init__` gives every
    analysis its own copy of `defaultepar`: a `PSS(epar=E)` would shoot at
    300 K whatever `E` says.

    PREDICTED DEFECT: the inner transient's `epar` is not the PSS's, and the
    PSS waveforms at T = 400 K and at the default are BIT-EQUAL, where a
    forward transient at 400 K differs visibly from one at 300 K."""
    from pycircuit.circuit.circuit import defaultepar
    E = defaultepar.copy()
    E.T = 400.0
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for label, kw in (('default', {}), ('T = 400 K', {'epar': E})):
            p = PSS(diode_rc(), method='radau', reltol=1e-9, **kw)
            p.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=40)
            assert p.converged
            out[label] = (np.asarray(p.waveform[1], float), p)
        fwd = {}
        for label, kw in (('default', {}), ('T = 400 K', {'epar': E})):
            c = diode_rc()
            tr = Transient(c, reltol=1e-9, **kw)
            res = tr.solve(tend=5e-3, timestep=1e-5, x0=np.zeros(c.n))
            fwd[label] = np.asarray(res.x, float)[:, -1]
    p_hot = out['T = 400 K'][1]
    tr = p_hot._transient()
    print('V12 the PSS at T = 400 K: pss.epar.T = %g, its transient\'s '
          'epar.T = %g, same object: %s'
          % (p_hot.epar.T, tr.epar.T, tr.epar is p_hot.epar))
    dW = float(np.max(np.abs(out['T = 400 K'][0] - out['default'][0])))
    dF = float(np.max(np.abs(fwd['T = 400 K'] - fwd['default'])))
    print('    PSS waveforms, 400 K against default: max |diff| %.3e '
          '(bit-equal: %s)' % (dW, dW == 0.0))
    print('    forward transients, 400 K against default, last point: '
          '%.3e' % dF)
    return int(dW == 0.0 and dF > 0.0)


## ---------------------------------------------------------------------------
## V13 -- the tests' hand-driven marches never clear `_is_first_step`
## ---------------------------------------------------------------------------

def V13():
    """`test_transient_branch_and_dae.py`'s `counts` (and `test_glm.py`'s
    `_march`, `benchmarks/defect_locality.py`) drive `solve_timestep` by
    hand: `_begin_run`, then per step `_dt_last`, `_dt`, `epar.t`, the step
    and `_push_history` -- but never `_is_first_step = False`.

    PREDICTED: under Gear-2 every step of such a march runs as Euler (the
    first-step order drop), so the march is not the multistep method its
    test names.  Counted: the integrator each step actually used
    (`_effective_method`), on the same march pattern."""
    from pycircuit.circuit.integrator import Gear2Integrator
    c = diode_rc()
    tr = Transient(c, integrator=Gear2Integrator(), reltol=1e-12)
    tr.irefnode = c.get_node_index(gnd)
    x = np.zeros(c.n)
    tr.epar.t = 0.0
    tr._begin_run(x, c.n)
    h = 1e-3 / 50
    used = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for j in range(1, 11):
            tr._dt_last = tr._dt if j > 1 else None
            tr._dt = h
            tr.epar.t = j * h
            x, _f, _J, _ = tr.solve_timestep(x, j * h)
            tr._push_history(x)
            used.append(getattr(tr, '_effective_method', None))
    print('V13 the march pattern under Gear2Integrator, integrator per step:')
    print('    %s' % used)
    return int(all(u == 'EulerIntegrator' for u in used))


## ---------------------------------------------------------------------------
## V7 -- GLM node startups run on the PSS's CURRENT transient
## ---------------------------------------------------------------------------

def V7():
    """`_factored_self_starting` walks the period in a transient of its own
    (`method=` may differ from the PSS's) and restores the PSS's; the node
    startups (`_glm_node_startups`, built lazily by
    `x_matvec_transposed(collect=True)`) then run on `self._transient()` --
    the PSS's own, not the walk's.

    PREDICTED DEFECT: on a glm2 PSS, `factored_period_glm(method='glm3')`
    collects node states that differ from a glm3 PSS's for the same period
    point (or fails on a shape mismatch, r = 3 against 4); on a PSS never
    solved, the lazily built transient has no `base_integrator` and the
    call raises."""
    per, npts = 1e-3, 40
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p2 = PSS(cv_loop(per), method='glm2', reltol=1e-12)
        p2.solve(period=per, timestep=per / npts, maxiterations=40)
        x0 = np.asarray(p2._period_state[1], float)
        p3 = PSS(cv_loop(per), method='glm3', reltol=1e-12)
        p3.solve(period=per, timestep=per / npts, maxiterations=40)
        m = p2.cir.n - 1
        v = np.linspace(1.0, 2.0, m)
        results = {}
        for label, pss, kw in (('glm2 PSS, method=glm3', p2, {'method': 'glm3'}),
                               ('glm3 PSS, own method', p3, {}),
                               ('never-solved glm3 PSS', PSS(cv_loop(per), method='glm3'), {})):
            try:
                fp = pss.factored_period_glm(x0, per, npts, **kw)
                out, _ts, states = fp.x_matvec_transposed(v, collect=True)
                results[label] = np.asarray([np.asarray(s_, float) for s_ in states])
            except Exception as e:                            # noqa: BLE001
                results[label] = '%s: %s' % (type(e).__name__, str(e)[:80])
    ref = results['glm3 PSS, own method']
    print('V7  GLM node startups, collected node states against a glm3 PSS\'s:')
    bad = 0
    for label, r in results.items():
        if isinstance(r, str):
            print('    %-26s %s' % (label, r))
            bad += 1
        elif isinstance(ref, str) or r.shape != ref.shape:
            print('    %-26s shape %s' % (label, getattr(r, 'shape', None)))
            bad += 1
        else:
            d = float(np.max(np.abs(r - ref)) / max(float(np.max(np.abs(ref))), 1e-300))
            print('    %-26s max rel diff %.3e' % (label, d))
            bad += int(d > 1e-9)
    return bad


## ---------------------------------------------------------------------------
## V8 -- the startup override skips recording the trace
## ---------------------------------------------------------------------------

def V8():
    """`_glm_startup` returns the override's vector BEFORE recording
    `_glm_startup_trace`; `_glm_period_blocks` reads the trace with no
    default.  The override is a test hook (`test_glm.py`) only.

    PREDICTED: with the override set, the walk reads a STALE trace (from
    the transient's previous startup) -- or raises AttributeError on a
    transient that never recorded one."""
    per, npts = 1e-3, 20
    times = np.linspace(0.0, per, npts + 1)
    hs = np.diff(times)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p = PSS(cv_loop(per), method='glm3', reltol=1e-12)
        p.solve(period=per, timestep=per / npts, maxiterations=40)
        x0 = np.asarray(p._period_state[1], float)
        other = p._new_transient(p._integrator_for('glm3'))
        other._begin_run(np.concatenate((x0[:p.irefnode], [0.0], x0[p.irefnode:])), p.cir.n)
        rows = []
        for label, fresh in (('fresh transient', True), ('used transient', False)):
            saved = p._tran
            if fresh:
                p._tran = p._new_transient(p._integrator_for('glm3'))
            tr = p._transient()
            before = getattr(tr, '_glm_startup_trace', None)
            tr._glm_startup_override = (
                lambda tn, x_, h: Transient._glm_startup(other, tn, x_, h))
            try:
                p._glm_period_blocks(x0 + 1e-3, times, hs)
                after = getattr(tr, '_glm_startup_trace', None)
                rows.append((label, 'ran; trace %s' % (
                    'unchanged (STALE)' if after is before else 'recorded')))
            except Exception as e:                            # noqa: BLE001
                rows.append((label, '%s: %s' % (type(e).__name__, str(e)[:70])))
            finally:
                try:
                    del tr._glm_startup_override
                except AttributeError:
                    pass
                p._tran = saved
    print('V8  a period walked with `_glm_startup_override` set:')
    for label, r in rows:
        print('    %-16s %s' % (label, r))
    return sum(1 for _l, r in rows if 'STALE' in r or 'Error' in r)


## ---------------------------------------------------------------------------
## V10 -- the periodic gauge shift skips the GLM's Nordsieck slots
## ---------------------------------------------------------------------------

def V10():
    """`_apply_periodic_shifts` wraps a periodic row (Idtmod) in the state,
    the `_qlast` ring and the predictor nodes, but not `_glm_Q` /
    `_glm_Q_at_entry`, whose first block is the charge vector
    (`Q_0 = q(x)` at a step's end, `_solve_timestep_glm`).

    PREDICTED DEFECT (forward GLM only): after K wraps the state row of
    `Q_0` sits K moduli away from `q(x)`; the waveform itself may still be
    right to the Newton tolerance.  Benign if the difference stays below
    1e-12."""
    from pycircuit.circuit.elements import VS, Idtmod
    from pycircuit.circuit.integrator import GLM2Integrator, RadauIIA3Integrator
    import pycircuit
    from pycircuit.circuit.toolkit import numeric
    pycircuit.circuit.circuit.default_toolkit = numeric

    def ramp():
        c = SubCircuit()
        nin, nout = c.add_node('in'), c.add_node('out')
        c['vin'] = VS(nin, gnd, v=1.0)
        c['R1'] = R(nout, gnd, r=1e3)
        c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0)
        return c
    res = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for label, integ in (('glm2', GLM2Integrator), ('radau', RadauIIA3Integrator)):
            c = ramp()
            tr = Transient(c, integrator=integ(), reltol=1e-9)
            r_ = tr.solve(tend=3.5, timestep=0.05, x0=np.zeros(c.n),
                          fixed_timestep=True)
            res[label] = (tr, np.asarray(r_.x, float))
    tr, X = res['glm2']
    rows = list(getattr(tr, '_periodic_rows', None) or [])
    x = X[:, -1]
    q = np.asarray(tr.cir.q(x, tr.epar), float)
    Q = getattr(tr, '_glm_Q', None)
    if Q is None or not rows:
        print('V10 no Nordsieck slot or no periodic row (rows %s)' % rows)
        return 0
    Q0 = np.asarray(Q[0][0], float)
    print('V10 forward glm2 on an Idtmod ramp (modulus 1, 3.5 s: 3 wraps), '
          'periodic rows %s:' % rows)
    worst = 0.0
    for row in rows:
        row = row if isinstance(row, (int, np.integer)) else row[0]
        d = float(Q0[row] - q[row])
        worst = max(worst, abs(d))
        print('    row %d: Q_0 - q(x) = %+.6f' % (row, d))
    dW = float(np.max(np.abs(res['glm2'][1] - res['radau'][1])))
    print('    waveform, glm2 against radau at the same fixed step: max |diff| %.3e' % dW)
    return int(worst > 1e-12)


## ---------------------------------------------------------------------------
## V2 -- `_push_history` without the caller's window
## ---------------------------------------------------------------------------

def V2():
    """The shooting rolls the history with `_push_history(x_full)`, no
    window, so when a periodic row (Idtmod) wraps on a step the states the
    WALK still holds (gear's `x_prev`, the pair's second block) keep the
    pre-wrap representative while `x_end` moves.

    PREDICTED BENIGN: the closure folds every block's periodic rows
    (`_fold_periodic`, ``F -= m round(F/m)``) and reads the wrap jump from
    fractional positions (`_wrap_jump`), so a whole-modulus difference in
    `x_prev` changes the closing residual by <= 1e-12 m.  Measured on the
    gear pair of `test_shooting_pss.py`'s Idtmod fixture (2000 V into a
    modulus-1 Idtmod: two wraps per period): the residual as is, and with
    `x_prev`'s periodic row moved by +-m."""
    from pycircuit.circuit.elements import VS, Idtmod
    per = 1e-3

    def build():
        c = SubCircuit()
        for nn in ('in', 'out', 'f'):
            c.add_node(nn)
        c['vin'] = VS('in', gnd, v=2000.0)
        c['I'] = Idtmod('in', gnd, 'out', gnd, modulus=1.0, ic=0.3)
        c['Rf'] = R('out', 'f', r=1e3)
        c['Cf'] = C('f', gnd, c=1e-7)
        c['Rl'] = R('out', gnd, r=1e5)
        return c
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p = PSS(build(), method='gear', reltol=1e-9)
        p.solve(period=per, timestep=per / 200, maxiterations=60)
        assert p.converged
        fp = p.factored_period()
    m_ = p.cir.n - 1
    solved, x0, xm1, times, hs, T, _x0u = p._period_state
    z0 = np.concatenate((np.asarray(x0, float)[:m_], np.asarray(xm1, float)[:m_]))
    xl = np.asarray(fp.x_last, float).ravel()[:m_]
    xp = np.asarray(fp.x_prev, float).ravel()[:m_]
    tms = np.asarray(times, float)
    F0 = p._close_periodic(z0, np.concatenate((xl, xp)), tms)
    rows = p._periodic_fold
    print('V2  gear pair on an Idtmod (periodic rows %s), closing residual '
          '|F| = %.3e:' % ([r for r, _m, _o in rows], float(np.max(np.abs(F0)))))
    worst = 0.0
    for r, mod, _o in rows:
        for sgn in (+1.0, -1.0):
            xp2 = xp.copy()
            xp2[r] += sgn * mod
            F1 = p._close_periodic(z0, np.concatenate((xl, xp2)), tms)
            d = float(np.max(np.abs(F1 - F0)))
            worst = max(worst, d / mod)
            print('    x_prev row %d moved by %+g: max |dF| = %.3e' % (r, sgn * mod, d))
    return int(worst > 1e-12)


## ---------------------------------------------------------------------------
## V5 -- the PSS never calls `cir.reset_state`
## ---------------------------------------------------------------------------

def V5():
    """`Transient.solve` resets element state (`cir.reset_state`) and feeds
    `cir.accept_step`; the PSS does neither, and `event_grid` walks
    `cir.next_event` before its solve.  An Idtmod's `next_event` is a
    prediction from its last ACCEPTED point (`_bp_cache`), "inf before a
    traversal has started" (`event_grid`'s docstring).

    PREDICTED: on a fresh circuit and after `lte_grid` (whose adaptive run
    ends at t >= T, so its predictions lie beyond the period) the event
    times agree; after a forward `Transient.solve(tend=0.3 T)` on the SAME
    circuit object a predicted wrap inside the period joins them -- the
    grid then depends on what ran before."""
    from pycircuit.circuit.elements import VPulse, Idtmod
    T = 1e-6

    def build():
        c = SubCircuit()
        for nn in ('in', 'out'):
            c.add_node(nn)
        c['vs'] = VPulse('in', gnd, v1=0.0, v2=1.0, td=0.05 * T, tr=0.01 * T,
                         tf=0.01 * T, pw=0.4 * T, per=T)
        c['I'] = Idtmod('in', gnd, 'out', gnd, modulus=0.1 * T)
        c['Rl'] = R('out', gnd, r=1e3)
        return c
    got = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        c = build()
        p = PSS(c, method='radau')
        p.event_grid(T, npts=100)
        got['fresh'] = list(p.event_times)
        c = build()
        Transient(c).solve(tend=0.3 * T, timestep=T / 200, x0=np.zeros(c.n))
        p = PSS(c, method='radau')
        p.event_grid(T, npts=100)
        got['after Transient(0.3 T)'] = list(p.event_times)
        c = build()
        q = PSS(c, method='radau')
        q.lte_grid(T, x0=np.zeros(c.n), tstab=3 * T)
        p = PSS(c, method='radau')
        p.event_grid(T, npts=100)
        got['after lte_grid'] = list(p.event_times)
    print('V5  event times (fractions of T) by what ran first on the circuit:')
    for k, v in got.items():
        print('    %-24s %s' % (k, ' '.join('%.4f' % e for e in v)))
    return int(any(v != got['fresh'] for v in got.values()))


## ---------------------------------------------------------------------------
## V1 -- no order drop after a landed edge on the shooting path
## ---------------------------------------------------------------------------

def V1():
    """The loop re-arms the first-step order drop after every breakpoint
    landing; the shooting never does, even at edges `break_events` placed
    on grid points.  `pss.py` states the cost as a CONSTANT, not an order
    ("a two-step formula takes an O(h^2 [x'']) hit at the ONE step after a
    corner").

    PREDICTED BENIGN: on the pulsed RC of `test_shooting_events.py` against
    its exact solution, gear and trap, the error without the drop over the
    error with it (the loop's rule, re-armed at the landed edges) stays in
    [0.3, 3] at every N, and both converge at slope ~2 (1.8 .. 2.2)."""
    import types
    from pycircuit.circuit.elements import VPulse
    T = 1e-6
    RR, CC = 1e3, 3e-10
    TAU = RR * CC
    TD, PW = 0.0125 * T, 0.4 * T
    trise = 0.02 * T

    def pulsed():
        c = SubCircuit()
        c['vs'] = VPulse(1, gnd, v1=0.0, v2=1.0, td=TD, tr=trise, tf=trise,
                         pw=PW, per=T)
        c['R'] = R(1, 2, r=RR)
        c['C'] = C(2, gnd, c=CC)
        return c

    def exact(ts):
        e = [0.0, TD, TD + trise, TD + trise + PW, TD + 2 * trise + PW, T]
        seg = [(e[0], e[1], 0.0, 0.0), (e[1], e[2], 0.0, 1.0 / trise),
               (e[2], e[3], 1.0, 0.0), (e[3], e[4], 1.0, -1.0 / trise),
               (e[4], e[5], 0.0, 0.0)]

        def prop(x0, t0, t1, a, b):
            one_e = -np.expm1(-(t1 - t0) / TAU)
            return x0 * (1.0 - one_e) + (a - b * TAU) * one_e + b * (t1 - t0)
        alpha, beta = 1.0, 0.0
        for (t0, t1, a, b) in seg:
            al = np.exp(-(t1 - t0) / TAU)
            alpha, beta = al * alpha, al * beta + prop(0.0, t0, t1, a, b)
        x0 = beta / (1.0 - alpha)
        out = np.empty(len(ts))
        for i, t in enumerate(ts):
            t = float(t) % T
            x = x0
            for (t0, t1, a, b) in seg:
                if t <= t1 + 1e-30:
                    out[i] = prop(x, t0, t, a, b)
                    break
                x = prop(x, t0, t1, a, b)
            else:
                out[i] = x
        return out
    edges = np.array([TD, TD + trise, TD + trise + PW, TD + 2 * trise + PW])

    def err(method, N, drop):
        c = pulsed()
        pss = PSS(c, method=method, reltol=1e-10)
        if drop:
            orig = pss.solve_timestep

            def st(self, x0, t, dt, *a, **k):
                if np.any(np.abs((float(t) - float(dt)) % T - edges) < 1e-9 * T):
                    self._transient()._is_first_step = True
                return orig(x0, t, dt, *a, **k)
            pss.solve_timestep = types.MethodType(st, pss)
        pss.solve(period=T, timestep=T / N, break_events=True, maxiterations=40)
        if not pss.converged:
            return float('nan')
        ts = np.asarray(pss.waveform[0], float).ravel()
        xs = np.asarray(pss.waveform[1], float)
        row = [str(n) for n in c.nodes].index('2')
        return float(np.abs(xs[row] - exact(ts)).max())
    print('V1  landed edges without / with the order drop (pulsed RC, '
          'tr = T/50), max error against the exact solution:')
    Ns = (100, 200, 400, 800)
    bad = 0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for method in ('gear', 'trap'):
            e0 = np.array([err(method, N, False) for N in Ns])
            e1 = np.array([err(method, N, True) for N in Ns])
            s0 = np.polyfit(np.log(Ns), np.log(e0), 1)[0]
            s1 = np.polyfit(np.log(Ns), np.log(e1), 1)[0]
            print('    %-5s no drop %s  slope %.2f' % (method, ' '.join('%.2e' % v for v in e0), -s0))
            print('    %-5s drop    %s  slope %.2f' % ('', ' '.join('%.2e' % v for v in e1), -s1))
            print('    %-5s ratio   %s' % ('', ' '.join('%.2f' % v for v in e0 / e1)))
            bad += int(np.any(e0 / e1 > 3.0) or -s0 < 1.5)
    return bad


## ---------------------------------------------------------------------------
## V4 -- the branch screen's scale spans the whole PSS
## ---------------------------------------------------------------------------

def V4():
    """The branch screen compares `C` against the LARGEST `C` seen
    (`_branch_cmax`); the loop resets it per run, the shooting never does,
    so on the shooting path the scale spans every period and iteration.

    PREDICTED BENIGN: the solution does not depend on it -- a fired screen
    only confirms by re-solving, and restores the limiter state -- so the
    PSS waveform is BIT-EQUAL with the scale reset at every `_begin_run`;
    only the screen counts may differ.  Needs a fixture on which the screen
    fires (a driven diode: its capacitance collapses under reverse bias)."""
    counts = {'calls': 0, 'fired': 0}

    def screen(orig):
        def f(self, x, *a, **k):
            out = orig(self, x, *a, **k)
            counts['calls'] += 1
            counts['fired'] += int(bool(out[0]))
            return out
        return f

    def begin(orig):
        def f(self, *a, **k):
            self._branch_cmax = 0.0
            self._branch_warned = False
            return orig(self, *a, **k)
        return f

    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from branch_selection import CubicCap
    from pycircuit.circuit.elements import G

    def build():
        ## (`C = c0 V^2` collapses at every zero crossing of the node)
        c = SubCircuit()
        c.add_node('in')
        c.add_node('a')
        c['vs'] = VSin('in', gnd, va=1.0, freq=0.2)
        c['R'] = R('in', 'a', r=1.0)
        c['cq'] = CubicCap('a', gnd, c0=1.0)
        c['g'] = G('a', gnd, g=0.1)
        return c
    ## a seed FAR off the orbit (V = 1e4, so C = 1e8): the largest C the
    ## analysis has seen stays there, and against it the orbit's C reads as
    ## collapsed -- the screen fires where a per-run scale would not
    c0_ = build()
    seed = np.zeros(c0_.n - 1)
    names = [str(n_) for n_ in c0_.nodes if str(n_) != 'gnd!']
    seed[names.index('a')] = 1e4
    seed[names.index('in')] = 0.0
    rows = {}
    for label, reset in (('as is', False), ('scale reset per run', True)):
        counts.update(calls=0, fired=0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            with patched(Transient, '_branch_screen', screen):
                if reset:
                    with patched(Transient, '_begin_run', begin):
                        p = PSS(build(), method='radau', reltol=1e-9)
                        p.solve(period=5.0, timestep=5.0 / 100, x0=seed,
                                maxiterations=60)
                else:
                    p = PSS(build(), method='radau', reltol=1e-9)
                    p.solve(period=5.0, timestep=5.0 / 100, x0=seed,
                            maxiterations=60)
        rows[label] = (np.asarray(p.waveform[1], float), dict(counts), p.converged)
    W0, c0, ok0 = rows['as is']
    W1, c1, ok1 = rows['scale reset per run']
    print('V4  branch screen on a driven cubic capacitor (C = V^2):')
    print('    as is:               screens %d, fired %d, converged %s' % (c0['calls'], c0['fired'], ok0))
    print('    scale reset per run: screens %d, fired %d, converged %s' % (c1['calls'], c1['fired'], ok1))
    same = W0.shape == W1.shape and np.array_equal(W0, W1)
    print('    waveforms bit-equal: %s' % same)
    return int(not same)

ITEMS = {'V6': V6, 'V9': V9, 'V3': V3, 'V12': V12, 'V13': V13,
         'V7': V7, 'V8': V8, 'V10': V10, 'V2': V2, 'V5': V5, 'V1': V1,
         'V4': V4}

if __name__ == '__main__':
    todo = sys.argv[1:] or list(ITEMS)
    for key in todo:
        ITEMS[key]()
