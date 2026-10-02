"""The sources are assembled once per step and time (speed round 4, stage
B; 2026-10-02).

`Transient._source_at` kept no memo: `u(t)` -- every source's `u`, and
the hdl elements' -- was assembled at every Newton iteration of a step,
at the same `t` each time: 2.03x per gear step on a PSP stage (31 us a
call in-run), 1.89x on a 20-PSP chain (126 us).  Now the first assembly
at a time serves every later request at that exact time within the step,
and the memo lives exactly as long as `solve_timestep`.  The same
function on the same inputs: the same bits (`U_MEMO = False` is the old
behaviour, for the comparison).
"""
import numpy as np
import pytest

from pycircuit.circuit import PSS, _tran_companion
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, VSin, gnd
from pycircuit.circuit.integrator import (
    Gear2Integrator,
    RadauIIA3Integrator,
    TrapezoidalIntegrator,
    TRBDF2Integrator,
)
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.transient import Transient

H = 2e-8


def _rect():
    """A half-wave rectifier: a sine through a SPICE diode into an RC."""
    cm.default_toolkit = numeric
    c = SubCircuit()
    for n in ('in', 'a', 'k'):
        c.add_node(n)
    c['vin'] = VSin('in', gnd, v=0.0, va=1.0, freq=1e6)
    c['rs'] = R('in', 'a', r=100.0)
    c['d'] = eh.DiodeSpiceHdl('a', 'k')
    c['rl'] = R('k', gnd, r=1e3)
    c['cl'] = C('k', gnd, c=1e-10)
    return c


def _pulsed():
    """A pulse train (breakpoints) and a sine into a diode-loaded RC."""
    cm.default_toolkit = numeric
    c = SubCircuit()
    for n in ('p', 's', 'a', 'k'):
        c.add_node(n)
    c['vp'] = VPulse('p', gnd, v1=0.0, v2=1.0, td=1e-7, tr=1e-8, tf=1e-8,
                     pw=2e-7, per=5e-7)
    c['vs'] = VSin('s', gnd, v=0.3, va=0.2, freq=3e6)
    c['r1'] = R('p', 'a', r=1e3)
    c['r2'] = R('s', 'a', r=2e3)
    c['d'] = eh.DiodeSpiceHdl('a', 'k')
    c['rl'] = R('k', gnd, r=1e3)
    c['cl'] = C('k', gnd, c=2e-11)
    return c


def _count(cir, name='u'):
    """`cir.<name>` wrapped on the instance to count its calls, as the
    counting tests of `test_transient_repairs` do."""
    real = getattr(cir, name)
    seen = [0]

    def counting(*a, **kw):
        seen[0] += 1
        return real(*a, **kw)
    setattr(cir, name, counting)
    return seen


def test_a_multistep_step_assembles_its_sources_once():
    """Gear, the full Newton (no chord), `uic` so no operating point adds
    calls of its own: `cir.u` is called once per accepted step although the
    Newton takes more than one iteration per step.  Before the memo it was
    called once per ITERATION (this fails on the parent)."""
    c = _rect()
    seen = _count(c)
    tr = Transient(c, toolkit=numeric, uic=True, chord_jacobian=False)
    tr.solve(tend=50 * H, timestep=H, fixed_timestep=True)
    steps = tr.statistics.accepted_steps
    assert tr.statistics.newton_iterations > 1.3 * steps, \
        'the circuit must need more than one iteration a step for this to test anything'
    assert seen[0] <= steps + 2, (seen[0], steps, tr.statistics.newton_iterations)


def test_the_chord_iterations_share_the_step_s_sources():
    """The chord (`chord_jacobian=True`) evaluates the residual alone on its
    later iterations (`residual_only`), at the same `t`: the memo serves
    those too."""
    c = _rect()
    seen = _count(c)
    tr = Transient(c, toolkit=numeric, uic=True, chord_jacobian=True)
    tr.solve(tend=50 * H, timestep=H, fixed_timestep=True)
    steps = tr.statistics.accepted_steps
    assert tr.statistics.newton_iterations > 1.3 * steps
    assert seen[0] <= steps + 2, (seen[0], steps)


def test_a_coupled_step_assembles_its_sources_once_per_stage_time():
    """Radau IIA(3): two stage times a step, every one of them read again
    at every iteration of the coupled Newton (`_stage_source`): with the
    memo, one assembly per stage time."""
    c = _rect()
    seen = _count(c)
    seen_i = _count(c, 'i')
    tr = Transient(c, toolkit=numeric, uic=True, integrator=RadauIIA3Integrator())
    tr.solve(tend=50 * H, timestep=H, fixed_timestep=True)
    steps = tr.statistics.accepted_steps
    ## three times a step: its start, the interior stage (c = 1/3) and its
    ## end (c = 1); measured 10.24 calls a step before the memo
    times = 3
    ## (the coupled Newton keeps its own count: `i` at every stage of every
    ## iteration says there was more than one iteration a step)
    assert seen_i[0] > 1.3 * 2 * steps
    assert seen[0] <= times * steps + 2, (seen[0], times, steps)


def test_provided_function_is_called_as_often_as_before(monkeypatch):
    """`provided_function(t)` is a caller's extra source and a caller may
    count it: it is called at every request, memo or not."""
    def run():
        c = _rect()
        calls = [0]

        def pf(t):
            calls[0] += 1
            return np.zeros(c.n)
        tr = Transient(c, toolkit=numeric, uic=True, chord_jacobian=False)
        tr.solve(tend=30 * H, timestep=H, fixed_timestep=True,
                 provided_function=pf)
        return calls[0], tr.statistics.newton_iterations
    on = run()
    monkeypatch.setattr(_tran_companion, 'U_MEMO', False)
    off = run()
    assert on == off
    assert on[0] >= on[1]                      # at least once per iteration


class _Stepping(VS):
    """A source whose value depends on a state only `accept_step` moves --
    what a transmission line is to the memo."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.k = 0

    def accept_step(self, t, x, epar):
        self.k += 1

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        return super().u(t, epar, analysis) * (1.0 + 0.1 * self.k)


def test_the_memo_lives_exactly_as_long_as_the_step():
    """Two steps solved from the same point at the same time, the source's
    state moved in between: the second sees the new state (the memo is
    opened and closed by `solve_timestep`), and outside a step there is no
    memo at all."""
    cm.default_toolkit = numeric
    c = SubCircuit()
    for n in ('a', 'k'):
        c.add_node(n)
    c['vs'] = _Stepping('a', gnd, v=0.8)
    c['rs'] = R('a', 'k', r=1e3)
    c['d'] = eh.DiodeSpiceHdl('k', gnd)
    c['cl'] = C('k', gnd, c=1e-11)
    tr = Transient(c, toolkit=numeric, uic=True)
    res = tr.solve(tend=10 * H, timestep=H, fixed_timestep=True)
    x0 = np.asarray(res.x, float)[:, -1].copy()
    t1 = float(res.sweep_values[-1]) + H
    tr._dt = H
    a = tr.solve_timestep(x0, t1)[0]
    assert tr._u_memo is None
    c['vs'].k += 5
    b = tr.solve_timestep(x0, t1)[0]
    assert not np.array_equal(a, b)
    ## the same step twice without a change: the same bits
    assert tr.solve_timestep(x0, t1)[0].tobytes() == b.tobytes()


def test_a_time_that_cannot_be_a_key_steps_aside():
    """A step asked at a time that is not hashable -- a 0-d numpy array
    here, a JAX array under that toolkit (the gate's first run failed
    `test_pss_and_pac_run_under_the_jax_toolkit` on exactly this) -- is
    assembled at every request, as before, and solves the same step."""
    c = _rect()
    seen = _count(c)
    tr = Transient(c, toolkit=numeric, uic=True, chord_jacobian=False)
    res = tr.solve(tend=10 * H, timestep=H, fixed_timestep=True)
    x0 = np.asarray(res.x, float)[:, -1].copy()
    t1 = float(res.sweep_values[-1]) + H
    tr._dt = H
    before = seen[0]
    a = tr.solve_timestep(x0, t1)[0]
    n_a = seen[0] - before
    before = seen[0]
    b = tr.solve_timestep(x0, np.array(t1))[0]
    n_b = seen[0] - before
    assert a.tobytes() == b.tobytes()
    assert n_a == 1 < n_b, (n_a, n_b)


def test_an_instance_patched_source_is_what_the_memo_serves():
    """`cir.u` replaced on the instance (as the shooting tests do with a
    perturbed source): the replacement is called once per step and its
    answer is what the step uses."""
    def run(memo):
        c = _rect()
        real = c.u
        calls = [0]

        def pert(t, epar=None, analysis=None, **kw):
            calls[0] += 1
            return np.asarray(real(t, epar, analysis=analysis, **kw), float) * 1.01
        c.u = pert
        _tran_companion.U_MEMO = memo
        try:
            tr = Transient(c, toolkit=numeric, uic=True, chord_jacobian=False)
            res = tr.solve(tend=40 * H, timestep=H, fixed_timestep=True)
        finally:
            _tran_companion.U_MEMO = True
        return np.asarray(res.x, float), calls[0], tr.statistics.accepted_steps
    x_on, n_on, steps = run(True)
    x_off, n_off, _ = run(False)
    assert x_on.tobytes() == x_off.tobytes()
    assert n_on <= steps + 2 < n_off


@pytest.mark.parametrize('integrator', [Gear2Integrator, TrapezoidalIntegrator,
                                        RadauIIA3Integrator, TRBDF2Integrator])
def test_transients_are_byte_identical_with_the_memo_off(integrator, monkeypatch):
    """Multistep and stage methods, a pulse train's breakpoints, the
    adaptive controller: the same bytes with the memo and without."""
    def run():
        tr = Transient(_pulsed(), toolkit=numeric, integrator=integrator())
        res = tr.solve(tend=1.2e-6, timestep=1e-8)
        return np.asarray(res.x, float), np.asarray(res.sweep_values, float)
    xa, ta = run()
    monkeypatch.setattr(_tran_companion, 'U_MEMO', False)
    xb, tb = run()
    assert ta.tobytes() == tb.tobytes()
    assert xa.tobytes() == xb.tobytes()


def test_a_shooting_pss_is_byte_identical_with_the_memo_off(monkeypatch):
    """The shooting re-enters every period at the same times with a new
    state: the memo, closed with every step, serves none of them a stale
    vector -- the same orbit to the bit."""
    def run():
        p = PSS(_rect(), method='gear', reltol=1e-6)
        ## (32 points do not resolve the rectifier's knee: the solve says so,
        ## and the bits are the question here)
        with quiet(AccuracyWarning):
            p.solve(period=1e-6, timestep=1e-6 / 32, maxiterations=40)
        return np.asarray(p.waveform[1], float)
    a = run()
    monkeypatch.setattr(_tran_companion, 'U_MEMO', False)
    b = run()
    assert a.tobytes() == b.tobytes()
