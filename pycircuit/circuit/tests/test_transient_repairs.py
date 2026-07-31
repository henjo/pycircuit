"""Regression tests for the transient repair work (`doc/transient_work_plan.md`).

Every test here corresponds to a **declared gate**. They exist because each defect
they cover produced a *confident wrong answer* rather than an error, which is the one
class of failure a test suite cannot be trusted to catch by accident.

The file was originally `test_transient_stage1.py` and outgrew that name as later
stages landed. It is named for the work rather than for a stage number, so it stays
true as stages 4-12 add to it; the plan's stage numbering is provenance, and it is
recorded per-section below rather than in the filename.

Sections, in file order:

    stage 1   the silent failures
              1-1  a failed operating point raises, and names the escape
              1-2  epar reaches the devices, in the inner DC and every step
              1-3  `bypass` is connected to something
              1-5  Newton's tolerance and the controller's are separate knobs
              0.1d the coupled path raises instead of livelocking

    stage 2+  the three improvements made after stage 2
              2+.1 `Toolkit.__getattr__` is memoised
              2+.2 the converged step skips the residual nobody reads
              2+.3 `relref` modes, and `sigglobal` on a quiet node

    stage 3   the first step
              3-1  `reltol` controls the global error
              3-3  a 2nd-order method beats backward Euler
                   plus the defect pinned deliberately, and `firststep` validation

    stage 4   step-control correctness
              4e   the order-drop guard watches growth, not shrink
                   plus the half of the gate that measurement refuted, pinned
              4b   the force-accept path warns, drops order, and stays inside
                   BDF-2's zero-stability bound

**A trap this file has fallen into three times** -- see `_pulsed_rc` and the note in
`test_sigglobal_does_not_collapse_on_a_quiet_node`: a transient whose `timestep` is
near the circuit's own time constant runs at `max_step` from end to end, so the
controller decides nothing and **no tolerance can be observed to act**. Two draft
gates "passed" that way while measuring nothing. Any new test of a step-control knob
must first establish that the controller is what binds -- start from `uic=True`, use a
`timestep` well above the step the controller would choose, and check that the step
count is well above `tend/timestep`.
"""

import numpy as np
import pytest

from pycircuit.circuit import circuit as circuit_mod
from pycircuit.circuit import numeric, gnd
from pycircuit.circuit.circuit import SubCircuit, defaultepar
from pycircuit.circuit.elements import VS, IS, R, C, Diode, VPulse
from pycircuit.circuit.analysis import NoConvergenceError, SingularMatrix
from pycircuit.circuit.transient import Transient


def _rc(toolkit=numeric):
    cir = SubCircuit(toolkit=toolkit)
    cir['VS'] = VS('n1', gnd, v=1.0)
    cir['R'] = R('n1', 'n2', r=1e3)
    cir['C'] = C('n2', gnd, c=1e-9)
    return cir


# ---------------------------------------------------------------------------
# Gate 1-1 -- the failure is visible
# ---------------------------------------------------------------------------

def test_gate_1_1_failed_operating_point_raises_and_names_uic():
    """A circuit with no DC solution must raise, not return a waveform.

    Before this gate the transient substituted a vector of zeros and ran to
    completion, so a circuit that was never biased produced a complete, plausible
    result indistinguishable from a correct one.
    """
    ## Two current sources in series into a floating node: no DC solution exists
    ## (KCL cannot be satisfied), and there is no DC path to ground from 'n2'.
    cir = SubCircuit(toolkit=numeric)
    cir['I1'] = IS(gnd, 'n1', i=1e-3)
    cir['I2'] = IS('n1', gnd, i=2e-3)
    cir['C1'] = C('n1', gnd, c=1e-9)

    tran = Transient(cir, toolkit=numeric)
    with pytest.raises((NoConvergenceError, SingularMatrix)) as excinfo:
        tran.solve(refnode=gnd, tend=1e-6, timestep=1e-7)

    ## The message must tell the caller how to proceed deliberately.
    assert 'uic' in str(excinfo.value), \
        'the error must name the uic escape, got: %s' % excinfo.value

    ## And no waveform may be left behind for a caller to pick up.
    assert getattr(tran, 'result', None) is None


def test_gate_1_1_uic_is_a_working_escape():
    """`uic=True` must still start from zeros, deliberately, without raising."""
    cir = SubCircuit(toolkit=numeric)
    cir['I1'] = IS(gnd, 'n1', i=1e-3)
    cir['I2'] = IS('n1', gnd, i=2e-3)
    cir['C1'] = C('n1', gnd, c=1e-9)

    tran = Transient(cir, toolkit=numeric)
    res = tran.solve(refnode=gnd, tend=1e-6, timestep=1e-7, uic=True) \
        if 'uic' in tran.solve.__code__.co_varnames else None
    if res is None:
        ## `uic` is a Parameter, not a solve() argument.
        tran = Transient(cir, toolkit=numeric, uic=True)
        res = tran.solve(refnode=gnd, tend=1e-6, timestep=1e-7)
    assert len(np.asarray(res.sweep_values)) > 1


# ---------------------------------------------------------------------------
# Gate 1-2 -- epar actually arrives
# ---------------------------------------------------------------------------

class _TemperatureProbe(R):
    """A resistor that records the temperature it was evaluated at."""

    def __init__(self, *args, **kwargs):
        super(_TemperatureProbe, self).__init__(*args, **kwargs)
        self.seen_T = []

    def i(self, x, epar=defaultepar):
        self.seen_T.append(getattr(epar, 'T', None))
        return super(_TemperatureProbe, self).i(x, epar)


def test_gate_1_2_epar_reaches_devices_in_dc_and_every_step():
    """A device must see the analysis's epar.T, not defaultepar's 300 K.

    The transient called `cir.i(x)` with no epar at all, so every device in every
    step saw the module-level default. Setting `T` had no effect on a transient,
    silently.
    """
    cir = SubCircuit(toolkit=numeric)
    cir['VS'] = VS('n1', gnd, v=1.0)
    probe = _TemperatureProbe('n1', 'n2', r=1e3)
    cir['R'] = probe
    cir['C'] = C('n2', gnd, c=1e-9)

    epar = defaultepar.copy()
    epar.T = 400.0

    tran = Transient(cir, toolkit=numeric, epar=epar)
    tran.solve(refnode=gnd, tend=1e-6, timestep=1e-7)

    seen = [t for t in probe.seen_T if t is not None]
    assert seen, 'the probe was never evaluated'
    assert all(t == 400.0 for t in seen), \
        'device saw temperatures %r; expected every evaluation at 400 K' \
        % sorted(set(seen))


# ---------------------------------------------------------------------------
# Gate 1-3 -- bypass is connected to something
# ---------------------------------------------------------------------------

def _diode_model_evals(bypass, bypasstol):
    """Count the model evaluations bypass exists to skip, and the step count.

    Counting calls to `Diode.i` does NOT work: the bypass test lives *inside*
    `i()` and skips the exponential, not the call. The observable is therefore
    `toolkit.exp`, which in this circuit only the diode uses.
    """
    cir = SubCircuit(toolkit=numeric)
    cir['VS'] = VS('n1', gnd, v=1.0)
    cir['R'] = R('n1', 'n2', r=1e3)
    cir['D'] = Diode('n2', gnd)
    cir['C'] = C('n2', gnd, c=1e-9)

    real_exp = numeric.exp
    count = {'n': 0}

    def counting_exp(*args, **kwargs):
        count['n'] += 1
        return real_exp(*args, **kwargs)

    numeric.exp = counting_exp
    try:
        kwargs = dict(bypass=bypass)
        if bypasstol is not None:
            kwargs['bypasstol'] = bypasstol
        tran = Transient(cir, toolkit=numeric, **kwargs)
        res = tran.solve(refnode=gnd, tend=2e-6, timestep=1e-7)
        return count['n'], len(np.asarray(res.sweep_values)), \
            float(np.asarray(res.v('n2'))[-1])
    finally:
        numeric.exp = real_exp


def test_gate_1_3_bypass_changes_evaluation_count():
    """`bypass=True` must do something measurable, and not change the answer.

    It previously did nothing at all: `Analysis.__init__` attaches `bypasstol` to
    the analysis's own epar, the transient never passed that epar to the devices,
    so every device took its `except AttributeError` branch and used a hard-coded
    1e-12 regardless of what the caller asked for.
    """
    off_evals, off_steps, off_v = _diode_model_evals(bypass=False, bypasstol=None)
    on_evals, on_steps, on_v = _diode_model_evals(bypass=True, bypasstol=1e-5)

    assert off_evals > 0 and on_evals > 0
    assert on_evals < off_evals, \
        'bypass=True/bypasstol=1e-5 evaluated the model %d times against %d ' \
        'without bypass -- the parameter is not connected' % (on_evals, off_evals)

    ## Bypassing IS an approximation -- it freezes the model within `bypasstol` --
    ## so the answer may move, but by less than the tolerance that licensed it.
    ## Measured: 1.7e-8 V of drift for a 1e-5 V bypass tolerance.
    assert on_steps == off_steps
    assert abs(on_v - off_v) < 1e-5, \
        'bypassing at 1e-5 V moved the answer by %g V, which is more than the ' \
        'tolerance that permitted it' % abs(on_v - off_v)


def test_gate_1_3_a_loose_bypasstol_is_not_silently_tolerated():
    """A bypass tolerance that breaks Newton must fail loudly.

    Found by this gate: with bypass genuinely connected, `bypasstol >= 1e-3`
    prevents the inner DC from converging on a plain diode circuit -- the model is
    frozen across steps far larger than the junction's own scale, so Newton is
    working with a stale Jacobian. That is the correct outcome for an absurd
    tolerance; what matters is that it raises rather than returning a waveform.
    """
    with pytest.raises((NoConvergenceError, SingularMatrix)):
        _diode_model_evals(bypass=True, bypasstol=1e-2)


def test_gate_1_3_defaultepar_carries_bypasstol():
    """Devices evaluated outside an analysis must find a defined bypasstol."""
    assert hasattr(defaultepar, 'bypasstol')
    ## Negative means "never bypass", which is what bypass=False should mean.
    assert defaultepar.bypasstol < 0


# ---------------------------------------------------------------------------
# Gate 1-5 -- the two tolerance roles are separate knobs
# ---------------------------------------------------------------------------

def test_gate_1_5_newton_and_lte_tolerances_are_separate():
    """`vabstol` is Newton's; `lte_vabstol` is the step controller's.

    **The two defaults are equal again** since gate D3-e returned `lte_vabstol` to
    1e-12, so the values no longer demonstrate anything about separateness -- they
    coincided before the split too, which is how one knob ended up serving both
    roles. What proves they are separate is behavioural and lives in the two tests
    below: `vabstol` must not move the step count, and `lte_vabstol` must (under
    `pointlocal`, where the absolute floor is load-bearing). Both are verified by
    injecting the defect, not by assertion alone.

    This test pins the *values*, which is still worth doing: a default that drifts
    silently is how the 10^6 asymmetry between the operating point and the steps
    after it arose in the first place.
    """
    tran = Transient(_rc(), toolkit=numeric)

    ## Newton's, and it agrees with DC's rather than being 10^6 looser.
    from pycircuit.circuit.dcanalysis import DC
    dc = DC(_rc(), toolkit=numeric)
    assert tran.par.vabstol == dc.par.vabstol, \
        'the operating point and the steps after it must be solved to the same ' \
        'tolerance; got transient %g vs DC %g' % (tran.par.vabstol, dc.par.vabstol)

    ## The controller's.  1e-12 since gate D3-e -- measured bit-identical to 1e-6
    ## under the shipped `sigglobal`, so the principled value costs nothing.
    assert tran.par.lte_vabstol == 1e-12
    assert tran.par.lte_iabstol == 1e-12
    assert tran.par.relref == 'sigglobal'  


def _pulsed_rc():
    """A circuit whose stepping is genuinely LTE-bound.

    A plain RC run at `timestep` >= its own time constant is bound by `max_step`
    from end to end, so the controller never decides anything and no tolerance can
    be observed to act. Pulse edges force the controller to take charge: this runs
    ~1150 steps against a `tend/timestep` of 30.
    """
    cir = SubCircuit(toolkit=numeric)
    cir['VS'] = VPulse('n1', gnd, v1=0, v2=5, td=1e-7, tr=1e-8, tf=1e-8,
                       pw=5e-7, per=1e-6)
    cir['R'] = R('n1', 'n2', r=1e3)
    cir['C'] = C('n2', gnd, c=1e-10)
    return cir


def _pulsed_run(**kwargs):
    tran = Transient(_pulsed_rc(), toolkit=numeric, **kwargs)
    res = tran.solve(refnode=gnd, tend=3e-6, timestep=1e-7)
    return (len(np.asarray(res.sweep_values)),
            float(np.asarray(res.v('n2'))[-1]))


def test_gate_1_5_lte_vabstol_moves_the_step_count_under_pointlocal():
    """The controller's tolerance must control the controller.

    **Pinned to `relref='pointlocal'` since decision D3.** The absolute floor is
    load-bearing only in that mode: `pointlocal` references the relative tolerance
    to each unknown's own present value, which goes to zero on a quiet node, so
    `abstol` is what stops the tolerance degenerating there. Under the shipped
    default `sigglobal` the reference is the largest signal in the circuit and
    cannot degenerate, so the floor is never reached on a healthy circuit and
    moving it changes nothing -- measured at gate D3-b as 403 steps at
    `lte_vabstol` 1e-3, 1e-6 and 1e-9 alike.

    That is the intended behaviour, not a regression, so this test keeps the
    `pointlocal` property (that mode still exists and is still selectable) and
    `test_gate_1_5_lte_vabstol_is_not_load_bearing_under_sigglobal` asserts the
    default's. Deleting this one would have dropped the only check that the
    controller's absolute tolerance reaches the controller at all.
    """
    ## Explicit values, not "the default" as the middle point: gate D3-e moved the
    ## default to 1e-12, which is at the tight END of any sweep rather than in the
    ## middle of one.  A test that assumes where the default sits breaks whenever
    ## the default moves, and says nothing useful when it does.
    loose, _ = _pulsed_run(lte_vabstol=1e-3, relref='pointlocal')
    mid, _ = _pulsed_run(lte_vabstol=1e-6, relref='pointlocal')
    tight, _ = _pulsed_run(lte_vabstol=1e-9, relref='pointlocal')

    assert loose < mid < tight, \
        'step count must fall monotonically as the LTE tolerance loosens; ' \
        'got %d (1e-3) / %d (1e-6) / %d (1e-9)' % (loose, mid, tight)


def test_gate_1_5_lte_vabstol_is_not_load_bearing_under_sigglobal():
    """And under the SHIPPED DEFAULT, the absolute floor stops mattering.

    This is the property that makes `sigglobal` worth having: an absolute tolerance
    exists to stop the relative one degenerating, and under `sigglobal` it cannot
    degenerate, so the floor is never reached. Asserted rather than assumed because
    it is what licenses `lte_vabstol` to hold any defensible value -- and because if
    it ever stops being true, the reference has started collapsing somewhere and
    that is worth knowing.
    """
    base, v_base = _pulsed_run()
    for value in (1e-3, 1e-9, 1e-12):
        steps, v = _pulsed_run(lte_vabstol=value)
        assert steps == base, \
            'lte_vabstol=%g moved the step count %d -> %d under sigglobal; the ' \
            'relative reference must be collapsing somewhere' % (value, base, steps)
        assert abs(v - v_base) < 1e-12 * max(abs(v_base), 1.0), \
            'lte_vabstol=%g moved the waveform under sigglobal (%g -> %g)' \
            % (value, v_base, v)


def test_gate_1_5_vabstol_does_not_move_the_step_count():
    """Newton's tolerance must NOT be a step-control knob.

    This is the property the shared parameter made impossible: there was no way to
    change one role without silently changing the other. A failure here means
    `vabstol` is feeding the step controller again.

    **RUN UNDER `pointlocal`, DELIBERATELY, AND THE REASON MATTERS.** This test was
    blind under the shipped default. Verified by injection: point the controller's
    `abstol` at `self.par.vabstol` -- i.e. reintroduce exactly the defect gate 1-5
    exists to catch -- and under `sigglobal` this test still **passes**, because the
    same property that makes the absolute floor non-load-bearing under `sigglobal`
    (the relative reference cannot degenerate, so `reltol*ref` dominates) also stops
    the injected `vabstol` from moving anything.

    A guard that cannot see the defect it guards is worse than no guard, because it
    reads as coverage. `pointlocal` is where an absolute floor is load-bearing, so
    that is where this has to run. Under injection it then fails as it should.
    """
    base_steps, base_v = _pulsed_run(relref='pointlocal')
    for value in (1e-9, 1e-14):
        steps, v = _pulsed_run(vabstol=value, relref='pointlocal')
        assert steps == base_steps, \
            'vabstol=%g moved the step count %d -> %d; it is still reaching the ' \
            'step controller' % (value, base_steps, steps)
        assert abs(v - base_v) < 1e-9 * max(abs(base_v), 1.0), \
            'vabstol=%g moved the waveform (%g -> %g)' % (value, base_v, v)


# ---------------------------------------------------------------------------
# 0.1d -- the coupled path must fail rather than grind
# ---------------------------------------------------------------------------

def test_coupled_non_convergence_raises_instead_of_livelocking():
    """`_solve_coupled` must not advance time on an unconverged step.

    Measured before the fix: each outer iteration cost 10 Newton solves and
    advanced simulated time by h*0.25^10, neither raising nor returning -- about
    2.1e7 further Newton attempts to finish a 5 us run.
    """
    cir = _rc()
    tran = Transient(cir, toolkit=numeric)
    real = tran.solve_timestep
    state = {'n': 0}

    def flaky(xlast, tnext, **kwargs):
        state['n'] += 1
        if state['n'] > 3:
            raise NoConvergenceError('forced, to exercise the exhaustion path')
        return real(xlast, tnext, **kwargs)

    tran.solve_timestep = flaky

    with pytest.raises(NoConvergenceError) as excinfo:
        tran.solve(refnode=gnd, tend=5e-6, timestep=1e-6, coupled_lte=True)

    ## Bounded work: 3 successes plus at most a few outer iterations' worth of
    ## retries, not the ~1e7 the livelock would take.
    assert state['n'] < 200, \
        'took %d solves before failing; the exhaustion path is not bounded' \
        % state['n']
    assert 'converge' in str(excinfo.value).lower()


# ---------------------------------------------------------------------------
# Post-stage-2 improvements (doc/transient_work_plan.md items 2+.1 .. 2+.3)
# ---------------------------------------------------------------------------

def test_toolkit_getattr_is_memoised():
    """2+.1 -- a resolved backend attribute lands in the instance __dict__.

    The point is that the SECOND access never enters ``__getattr__`` at all.
    Instance assignment must still win, because two test/benchmark files rely on
    monkeypatching a toolkit primitive.
    """
    ## A fresh instance, so the memo state is known.  `zeros` is used because it
    ## genuinely exists on the numeric backend -- `reshape` does not, which is
    ## itself why the hoisted probes in stage 2b mattered.
    fresh = type(numeric)(numeric._backend)
    assert 'zeros' not in fresh.__dict__
    first = fresh.zeros
    assert 'zeros' in fresh.__dict__, 'lookup was not memoised'
    assert fresh.zeros is first

    ## A missing attribute still raises, still names the toolkit, and is NOT
    ## cached -- so an attribute that appears later would be found.
    with pytest.raises(AttributeError) as excinfo:
        fresh.definitely_not_an_attribute
    assert 'toolkit' in str(excinfo.value)
    assert 'definitely_not_an_attribute' not in fresh.__dict__

    ## Assignment onto the instance overrides any memo.
    saved = fresh.zeros
    fresh.zeros = 'patched'
    assert fresh.zeros == 'patched'
    fresh.zeros = saved


def test_converged_step_skips_the_residual_nobody_reads():
    """2+.2 -- `i` and `u` are not assembled at the converged point.

    `solve_timestep` returns `(x, feval, J, f)` and `f` is never read on the
    standard path, so assembling it was pure waste. With `provided_function` set
    -- the one caller that does read `f` -- the full evaluation must come back.
    """
    def counts(provided_function):
        cir = _rc()
        seen = {}

        for name in ('i', 'u', 'G'):
            real = getattr(cir, name)

            def make(nm, fn):
                def counting(*a, **kw):
                    seen[nm] = seen.get(nm, 0) + 1
                    return fn(*a, **kw)
                return counting
            setattr(cir, name, make(name, real))

        tran = Transient(cir, toolkit=numeric)
        res = tran.solve(refnode=gnd, tend=20e-6, timestep=1e-6,
                         provided_function=provided_function)
        steps = len(np.asarray(res.sweep_values))
        return {k: v / steps for k, v in seen.items()}

    lean = counts(None)
    full = counts(lambda f, J, C: None)

    ## G is needed either way and is the control: if it moves, something other
    ## than the residual changed.
    assert abs(lean['G'] - full['G']) < 0.05
    assert lean['i'] < full['i'] - 0.5, \
        'i should be assembled once less per step without provided_function ' \
        '(%.2f vs %.2f)' % (lean['i'], full['i'])
    assert lean['u'] < full['u'] - 0.5


@pytest.mark.parametrize('relref', ['pointlocal', 'alllocal', 'sigglobal'])
def test_relref_modes_all_run_and_agree(relref):
    """2+.3 -- every relref mode produces the same answer, at its own cost.

    `relref` changes what the tolerance is measured against, not what is being
    solved, so the waveform must agree; only the step count may differ.
    """
    tran = Transient(_rc(), toolkit=numeric, relref=relref)
    res = tran.solve(refnode=gnd, tend=20e-6, timestep=1e-6)
    v = float(np.asarray(res.v('n2'))[-1])
    ## 20 us is 20 tau, so the RC has settled to the source voltage.
    assert abs(v - 1.0) < 1e-3, '%s gave v(n2)=%g' % (relref, v)


def test_relref_rejects_an_unknown_mode():
    tran = Transient(_rc(), toolkit=numeric, relref='nonsense')
    with pytest.raises(ValueError) as excinfo:
        tran.solve(refnode=gnd, tend=1e-6, timestep=1e-7)
    assert 'relref' in str(excinfo.value)


def test_sigglobal_does_not_collapse_on_a_quiet_node():
    """2+.3b, the point of the item, pinned as a test.

    A node carrying no signal makes `pointlocal`'s reference tend to zero, so its
    tolerance collapses to the absolute floor and the controller starts chasing
    numerical noise. `sigglobal` references it to the largest signal in the
    circuit instead. Measured on the leapfrog: 3.53x more steps under
    `pointlocal` at `lte_vabstol=1e-12`, against 1.49x under `sigglobal`.
    """
    def steps(relref):
        cir = _rc()
        ## The quiet node: driven only through 1 G-ohm, so |x| there is ~0.
        cir['Rq'] = R('n3', gnd, r=1e9)
        cir['Rq2'] = R('n3', 'n2', r=1e9)
        ## uic=True matters: started from the DC point the RC is ALREADY settled,
        ## so there is no truncation error, the run sits at max_step throughout and
        ## no tolerance can be observed to act.  Starting from zero forces the
        ## controller to actually control something.  (Third time this trap has
        ## been hit in this file -- see the note on `_pulsed_rc`.)
        tran = Transient(cir, toolkit=numeric, relref=relref,
                         lte_vabstol=1e-12, reltol=1e-6, uic=True)
        res = tran.solve(refnode=gnd, tend=20e-6, timestep=4e-6)
        return len(np.asarray(res.sweep_values))

    local = steps('pointlocal')
    glob = steps('sigglobal')
    assert glob < local, \
        'sigglobal (%d steps) should need no more steps than pointlocal (%d) ' \
        'when a quiet node is present' % (glob, local)


# ---------------------------------------------------------------------------
# Stage 3 -- the first step (doc/transient_work_plan.md)
# ---------------------------------------------------------------------------

def _rc_analytic(reltol, firststep='unset', integrator=None, tend_tau=10):
    """RC step response and its exact global error.

    Started from `uic=True` and with `timestep` deliberately far above the step
    the controller would pick, so the controller -- not `max_step` -- is what
    limits the run. Both are required for any tolerance to be observable.
    """
    from pycircuit.circuit.integrator import TrapezoidalIntegrator
    V, Rv, Cv = 1.0, 1e3, 1e-9
    tau = Rv * Cv
    cir = SubCircuit(toolkit=numeric)
    cir['VS'] = VS('n1', gnd, v=V)
    cir['R'] = R('n1', 'n2', r=Rv)
    cir['C'] = C('n2', gnd, c=Cv)

    opts = dict(reltol=reltol, vabstol=1e-12, lte_vabstol=1e-9, uic=True,
                integrator=integrator or TrapezoidalIntegrator())
    if firststep != 'unset':
        opts['firststep'] = firststep
    res = Transient(cir, toolkit=numeric, **opts).solve(
        refnode=gnd, tend=tend_tau * tau, timestep=tau)
    t = np.asarray(res.v('n2').x[0])
    v = np.asarray(res.v('n2').y)
    exact = V * (1.0 - np.exp(-t / tau))
    return len(t), float(np.max(np.abs(v - exact)))


def test_gate_3_1_reltol_controls_the_global_error():
    """Stage 3 -- the defect this stage exists to remove.

    Before the opening step was ramped, the global error was **1.3212e-01 at
    reltol 1e-3, 1e-4, 1e-5 and 1e-6 alike** -- identical to five digits -- while
    the step count went 24 -> 195. Eight times the work for the same answer,
    because the one step the controller cannot check was also the largest.
    """
    errors = [_rc_analytic(tol)[1] for tol in (1e-3, 1e-4, 1e-5, 1e-6)]

    assert all(errors[i] >= errors[i + 1] for i in range(len(errors) - 1)), \
        'global error must fall monotonically with reltol, got %r' % errors
    reduction = errors[0] / errors[-1]
    assert reduction >= 20, \
        'reltol 1e-3 -> 1e-6 reduced the error only %.1fx (%r)' % (reduction, errors)


def test_gate_3_3_second_order_beats_first_order():
    """A 2nd-order method must be more accurate than backward Euler.

    They used to agree to five digits, which is the clearest single symptom that
    the error was set by the opening step and not by the integration method.
    """
    from pycircuit.circuit.integrator import EulerIntegrator, TrapezoidalIntegrator

    _, err_euler = _rc_analytic(1e-5, integrator=EulerIntegrator())
    _, err_trap = _rc_analytic(1e-5, integrator=TrapezoidalIntegrator())
    assert err_trap < err_euler / 2.0, \
        'trapezoidal (%.4e) should be clearly better than euler (%.4e)' \
        % (err_trap, err_euler)


def test_opening_at_max_step_still_destroys_reltol():
    """The defect, pinned deliberately, so a regression is unmistakable.

    Asking for `firststep=timestep` restores the old behaviour exactly, and with
    it the symptom: reltol stops mattering. If this test ever starts failing,
    either the ramp has been made unconditional or something else now bounds the
    first step -- both worth knowing.
    """
    from pycircuit.circuit.integrator import TrapezoidalIntegrator
    tau = 1e-6
    hi = _rc_analytic(1e-3, firststep=tau)[1]
    lo = _rc_analytic(1e-6, firststep=tau)[1]
    assert abs(hi - lo) < 1e-9 * max(hi, 1.0), \
        'opening at max_step should make reltol irrelevant, got %.6e vs %.6e' \
        % (hi, lo)


def test_firststep_is_validated_and_capped():
    tran = Transient(_rc(), toolkit=numeric, firststep=-1.0)
    with pytest.raises(ValueError) as excinfo:
        tran.solve(refnode=gnd, tend=1e-6, timestep=1e-7)
    assert 'firststep' in str(excinfo.value)

    ## Larger than max_step is capped rather than honoured: the controller would
    ## clamp it on the next step anyway.
    tran = Transient(_rc(), toolkit=numeric)
    assert tran._opening_step(1e-6) == pytest.approx(1e-9)
    tran = Transient(_rc(), toolkit=numeric, firststep=1.0)
    assert tran._opening_step(1e-6) == 1e-6


# ---------------------------------------------------------------------------
# stage 4 -- step-control correctness
#
# 4e (the order-drop guard watched the wrong direction) and 4b (the force-accept
# path grew 10x, outside BDF-2's zero-stability bound) are the same defect from
# two sides: 4e's missing upper guard is why 4b's 10x had nothing to catch it.
# ---------------------------------------------------------------------------

def test_gate_4e_growth_past_the_stability_bound_drops_order():
    """The guard that was missing entirely.

    Variable-step BDF-2 is zero-stable only below `1 + sqrt(2)`. Before this,
    `Gear2.check_order_drop` tested the *shrink* direction and nothing else, so a
    ratio of 10 -- which `transient.py`'s force-accept path handed it -- passed
    through to a second-order method whose parasitic root is 4.76.
    """
    from pycircuit.circuit.integrator import (Gear2Integrator, EulerIntegrator,
                                              ZERO_STABILITY_RATIO)
    g = Gear2Integrator()

    just_over = ZERO_STABILITY_RATIO * 1.0001
    assert isinstance(g.check_order_drop(just_over, 1.0, False), EulerIntegrator), \
        'a growth ratio above the zero-stability bound must drop order'
    ## The 10x the old force-accept path took, which is what 4b removed.
    assert isinstance(g.check_order_drop(10.0, 1.0, False), EulerIntegrator)

    just_under = ZERO_STABILITY_RATIO * 0.9999
    assert g.check_order_drop(just_under, 1.0, False) is g, \
        'a growth ratio inside the bound must keep the 2nd-order method'
    ## The controller's own clamp: every normal accepted step lands here, so the
    ## guard above must be inert in a healthy run or it would cost an order of
    ## accuracy on every step that grows at all.
    from pycircuit.circuit.stepcontroller import MAX_GROWTH_RATIO
    assert MAX_GROWTH_RATIO < ZERO_STABILITY_RATIO
    assert g.check_order_drop(MAX_GROWTH_RATIO, 1.0, False) is g


def test_gate_4e_a_moderate_shrink_does_not_drop_order():
    """Half of gate 4e as declared: shrinking is not a stability problem."""
    from pycircuit.circuit.integrator import Gear2Integrator
    g = Gear2Integrator()
    for ratio in (0.9, 0.5, 0.2, 0.11):
        assert g.check_order_drop(ratio, 1.0, False) is g, \
            'ratio %g is a shrink and unconditionally zero-stable for BDF-2' % ratio


def test_gate_4e_deep_shrink_still_drops_order_deliberately():
    """The other half of gate 4e as declared was REFUTED, and this pins it.

    The gate said "a shrink does not" trigger the drop, and the plan said replace
    the shrink test with the growth test. Measurement disagreed: deleting the
    shrink branch took `Gear2('ywr')` and `Gear2('classic')` from 0 force-accepts
    to 1 each on the stiff RLC at reltol 1e-5, because that branch is not idle --
    it fires 3-6 times a run there and each firing is a step that then passes its
    Euler error test instead of being force-accepted over tolerance.

    So it is kept, and re-labelled: it is a stalled-estimate heuristic, not a
    stability guard. A step only shrinks 10x below the last accepted one after
    several consecutive rejections, and what rejects repeatedly is a 2nd-order
    estimate built on a third difference of something that is not three times
    differentiable. **If this test is ever deliberately inverted**, re-run the
    force-accept counts before believing the shrink branch is free to remove.
    """
    from pycircuit.circuit.integrator import Gear2Integrator, EulerIntegrator
    g = Gear2Integrator()
    assert isinstance(g.check_order_drop(0.05, 1.0, False), EulerIntegrator)


class _RejectionInjector:
    """A step controller that forces the rejection cap to be reached, once.

    Reaching the cap on a real circuit needs a configuration that is itself
    defective (see the end-to-end test below). This drives the mechanism
    directly instead, so the test survives the estimators being repaired.
    """

    def __init__(self, after=4, n_reject=4):
        from pycircuit.circuit.stepcontroller import IntegralController
        self.inner = IntegralController()
        self.after = after
        self.n_reject = n_reject
        self.accepted = 0
        self.injected = 0
        self.seen = []          # (h_curr, accepted, effective_method)
        self.force_idx = None   # index in `seen` of the force-accepted step
        self.transient = None

    def set_relref(self, relref):
        self.inner.set_relref(relref)
        return self

    def evaluate_step(self, *args, **kwargs):
        accept, h_next = self.inner.evaluate_step(*args, **kwargs)
        h = kwargs['h_curr']
        method = getattr(self.transient, '_effective_method', None)
        if self.accepted >= self.after and self.injected < self.n_reject:
            self.injected += 1
            ## The n_reject'th consecutive rejection is the one `transient.py`
            ## force-accepts.  Recording the index here rather than searching for
            ## it afterwards: a search for "the 4th rejection in the run" would
            ## silently pick the wrong step if the circuit rejects one of its own
            ## before the injection point.
            if self.injected == self.n_reject:
                self.force_idx = len(self.seen)
            self.seen.append((h, False, method))
            return False, h * 0.5
        if accept:
            self.accepted += 1
        self.seen.append((h, bool(accept), method))
        return accept, h_next


def _run_with_injected_rejections():
    import warnings as _w
    from pycircuit.circuit.integrator import Gear2Integrator
    tran = Transient(_pulsed_rc(), toolkit=numeric, integrator=Gear2Integrator())
    ctrl = _RejectionInjector()
    ctrl.transient = tran
    tran.step_controller = ctrl
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter('always')
        tran.solve(refnode=gnd, tend=3e-6, timestep=1e-7)
    forced = [c for c in caught
              if issubclass(c.category, RuntimeWarning)
              and 'truncation error' in str(c.message)]
    return ctrl, forced


def test_gate_4b_force_accept_warns_and_names_t_and_h():
    """An unbounded accepted truncation error must not be invisible.

    Same failure class as the silent operating point stage 1 removed: the result
    is still returned, and without the warning the caller has no way to learn
    that part of it was never error-controlled.
    """
    ctrl, forced = _run_with_injected_rejections()
    assert ctrl.injected == ctrl.n_reject, \
        'the injector did not reach the rejection cap -- MAX_REJECT may have moved'
    assert len(forced) == 1, \
        'expected exactly one force-accept warning, got %d' % len(forced)
    message = str(forced[0].message)
    assert 't=' in message and 'h=' in message, \
        'the warning must name the time and the step size, got: %s' % message


def test_gate_4b_force_accept_growth_stays_inside_the_stability_bound():
    """It grew 10x. BDF-2 is zero-stable only below 2.414214.

    Growing tenfold in answer to an error that was already too large is the wrong
    sign; what a stalled high-order estimate is asking for is a lower order. Both
    halves are asserted here: the growth is bounded, and the step that follows
    actually runs at order 1.
    """
    from pycircuit.circuit.stepcontroller import MAX_GROWTH_RATIO
    ctrl, forced = _run_with_injected_rejections()

    ## The force-accepted step is the last injected rejection; the accepted step
    ## after it is the one whose ratio used to be 10.
    idx = ctrl.force_idx
    assert idx is not None, 'the rejection cap was never reached'
    h_forced = ctrl.seen[idx][0]
    h_after, accepted_after, method_after = ctrl.seen[idx + 1]

    assert accepted_after, 'the step after a force-accept was itself rejected'
    assert h_after / h_forced <= MAX_GROWTH_RATIO * (1 + 1e-12), \
        'the step after a force-accept grew %.4gx, outside the %.4gx clamp every ' \
        'other accepted step obeys' % (h_after / h_forced, MAX_GROWTH_RATIO)
    assert h_after > h_forced, \
        'the escape hatch must still escape the collapsed regime'
    assert method_after == 'EulerIntegrator', \
        'the step after a force-accept must run at order 1, not %s' % method_after


def test_gate_4b_no_accepted_step_ratio_leaves_the_stability_bound():
    """End to end, in the SHIPPED DEFAULT configuration.

    `Gear2('ywr')` at reltol 1e-3 on this stiff RLC reaches the rejection cap
    once, and before the repair that one force-accept took a step ratio of
    **exactly 10.0** -- a single accepted step outside the zero-stability bound,
    in the configuration every caller gets by default. Measured before and after
    by re-running against the stashed source:

        before   179 accepted steps, 1 force-accept, 1 ratio above 2.414,
                 worst 10.0000
        after    181 accepted steps, 1 force-accept, 0 above the bound,
                 worst 2.0000

    **How this case was found is the point.** The first sweep for gate 4b covered
    reltol 1e-4/1e-5/1e-6, where `Gear2` never reaches the cap and only
    `Trapezoidal('ywr')` does -- which produced a written conclusion that the
    default no longer reached it at all. The warning added by 4b refuted that on
    its own first full-suite run, by firing here. It takes a *loose* tolerance,
    and 1e-3 was outside the swept range.

    `Trapezoidal('ywr')` reaches the cap far more often -- 25 times on this same
    circuit at reltol 1e-4, with 8 ratios above the bound, worst 6.4669 -- and
    would make a louder test. It is deliberately not used here: stage 4d deletes
    that branch as an alias of `'classic'`, so a test built on it would quietly
    stop testing anything the day 4d lands.
    """
    import warnings as _w
    from pycircuit.circuit.integrator import (Gear2Integrator,
                                              ZERO_STABILITY_RATIO)
    from pycircuit.circuit.stepcontroller import IntegralController
    from pycircuit.circuit.elements import L

    class _Spy(IntegralController):
        """Records the ACCEPTED step sequence, force-accepts included.

        A force-accepted step is an accepted step -- it enters the integrator
        history like any other, so it belongs in the ratio sequence. Filtering on
        `accept` alone drops exactly the steps this gate is about.
        """

        def __init__(self):
            super().__init__()
            self.h_accepted = []
            self.n_forced = 0
            self._consec = 0

        def evaluate_step(self, *args, **kwargs):
            accept, h_next = super().evaluate_step(*args, **kwargs)
            forced = False
            if not accept:
                self._consec += 1
                if self._consec > 3:        # MAX_REJECT
                    forced = True
                    self.n_forced += 1
                    self._consec = 0
            else:
                self._consec = 0
            if accept or forced:
                self.h_accepted.append(kwargs['h_curr'])
            return accept, h_next

    cir = SubCircuit(toolkit=numeric)
    cir['C1'] = C(1, gnd, c=1e-6)
    cir['R1'] = R(1, 2, r=1.0)
    cir['L1'] = L(2, gnd, L=1e-6)
    x0 = np.zeros(cir.n)
    x0[cir.get_node_index('1')] = 1.0
    x0[cir.get_node_index('2')] = 1.0

    tran = Transient(cir, toolkit=numeric,
                     integrator=Gear2Integrator('ywr'), reltol=1e-3)
    spy = _Spy()
    tran.step_controller = spy
    with _w.catch_warnings():
        _w.simplefilter('ignore', RuntimeWarning)
        tran.solve(refnode=gnd, tend=5e-3, timestep=2e-4, x0=x0)

    assert spy.n_forced > 0, \
        'this run no longer reaches the rejection cap, so it no longer tests 4b. ' \
        'That is worth knowing either way -- re-measure with ' \
        'benchmarks/transient_stage4.py --forceaccept and re-point this test at a ' \
        'configuration that does, rather than deleting it'

    ratios = [spy.h_accepted[i + 1] / spy.h_accepted[i]
              for i in range(len(spy.h_accepted) - 1)]
    worst = max(ratios)
    assert worst <= ZERO_STABILITY_RATIO, \
        '%d of %d accepted step ratios exceed the zero-stability bound %.6f, ' \
        'worst %.4f' % (sum(1 for r in ratios if r > ZERO_STABILITY_RATIO),
                        len(ratios), ZERO_STABILITY_RATIO, worst)


# ---------------------------------------------------------------------------
# 4g(b) -- the trapezoidal estimator differences a mode-free quantity
#
# The trapezoidal companion `iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}` has homogeneous
# solution `(-1)^n`, undamped.  Differencing `g = iq` therefore differences that
# mode.  `d_n = (q_n - q_{n-1})/h_n` does not: the trapezoidal relation gives
# `(iq_n + iq_{n-1})/2 = d_n`, and a component that flips sign each step cancels
# exactly in that sum.
# ---------------------------------------------------------------------------

## An analytic charge, so the local truncation error is exactly computable.
_4G_OMEGA = 2 * np.pi * 1e6


def _4g_q(t):
    return np.sin(_4G_OMEGA * t)


def _4g_dq(t):
    return _4G_OMEGA * np.cos(_4G_OMEGA * t)


def _4g_d3q(t):
    return -(_4G_OMEGA ** 3) * np.cos(_4G_OMEGA * t)


def _trap_est_over_true(h_curr, h_last, h_last2, t_n=0.371e-6, formula='classic'):
    """est / (local truncation error), for the trapezoidal estimator.

    The reference is the LOCAL truncation error `-(h^2/6) q'''`, not
    `iq_n - q'(t_n)`.  Gate 4g-b0 is why: the latter contains error propagated from
    earlier steps, which is O(h^3) here -- one order SMALLER than the local error --
    so dividing by it produced the "1/h law" that stage 4g was originally chasing.
    A per-step estimator cannot track the propagated part and must not try: the
    parasitic mode is undamped, so shrinking h does not reduce it.
    """
    from pycircuit.circuit.integrator import TrapezoidalIntegrator
    t_1 = t_n - h_curr
    t_2 = t_1 - h_last
    t_3 = t_2 - h_last2
    est, _p = TrapezoidalIntegrator(formula).compute_lte(
        q_curr=np.array([_4g_q(t_n)]), h_curr=h_curr,
        q_last=np.array([[_4g_q(t_1)], [_4g_q(t_2)], [_4g_q(t_3)]]),
        ## The estimator must never need these once h_last2 is available; passing
        ## deliberately wrong values is how this test proves it.
        iq_last=np.array([[np.nan], [np.nan]]),
        h_last=h_last, is_first_step=False, toolkit=numeric, h_last2=h_last2)
    ## The local truncation error is the EXACT residual when exact values are
    ## substituted for the history -- here, the trapezoidal companion built from
    ## exact charges and an exact `iq_{n-1}`, minus the exact derivative.
    ## `-(h^2/6) q3` is only its leading term and differs by 3.4% at h=1e-8, which
    ## is enough to fail a 10% gate for a reason that has nothing to do with the
    ## estimator.  The two are checked against each other in
    ## `test_gate_4g_b_local_error_matches_the_derivation`.
    true_local = (2.0 * (_4g_q(t_n) - _4g_q(t_1)) / h_curr - _4g_dq(t_1)) - _4g_dq(t_n)
    return float(est[0]) / true_local


def test_gate_4g_b_local_error_matches_the_derivation():
    """The derivation gives the trapezoidal LTE as -(h^2/6) times the third
    derivative of the charge.  Checked against the exact residual, which must
    approach it as h falls.

    This is the check that YWR eq (22) was evaluated correctly for TRAP: the
    coefficient is what the whole estimator is scaled by, and getting it wrong
    would bias every step of every trapezoidal run by a constant factor.
    """
    t_n = 0.371e-6
    prev = None
    for h in (1e-8, 1e-9, 1e-10):
        exact = (2.0 * (_4g_q(t_n) - _4g_q(t_n - h)) / h
                 - _4g_dq(t_n - h)) - _4g_dq(t_n)
        leading = -(h ** 2 / 6.0) * _4g_d3q(t_n)
        ratio = exact / leading
        assert abs(ratio - 1.0) < 0.05, \
            'exact LTE / leading term = %.6f at h=%g' % (ratio, h)
        if prev is not None:
            assert abs(ratio - 1.0) < abs(prev - 1.0), \
                'the leading term is not the limit: %s' % [prev, ratio]
        prev = ratio


def test_gate_4g_b_estimator_does_not_touch_the_companion_currents():
    """The whole point, expressed as an assertion.

    `iq_last` is passed as NaN. If the estimator reads it at all the result is NaN,
    so a finite answer is proof that the parasitic mode cannot reach the estimate.
    """
    r = _trap_est_over_true(1e-9, 1e-9, 1e-9)
    assert np.isfinite(r), \
        'the estimator still reads iq_last; the (-1)^n mode can reach the estimate'


def test_gate_4g_b1_estimator_is_asymptotically_exact():
    """Measured before: 0.8067 / 0.6780 / 0.6678 -- a ~33% underestimate that does
    not improve with h. After: 0.9273 / 0.9933 / 0.9993."""
    ratios = [_trap_est_over_true(h, h, h) for h in (1e-8, 1e-9, 1e-10)]
    for h, r in zip((1e-8, 1e-9, 1e-10), ratios):
        assert abs(r - 1.0) < 0.1, \
            'est/true = %.4f at h=%g; the estimator is not consistent' % (r, h)
    ## Consistency means it gets BETTER with h, which the g-based form did not.
    assert abs(ratios[2] - 1.0) < abs(ratios[0] - 1.0), \
        'est/true does not converge to 1 as h falls: %s' % ratios


def test_gate_4g_b2_estimator_is_flat_across_the_reachable_step_ratios():
    """Recorded spread: 1540x for the g-based 'ywr' form, 4.0x for 'classic'.

    The range asserted is the range the controller can actually produce.
    `MAX_GROWTH_RATIO` caps growth at 2.0 and 4b removed the only path that
    bypassed it, so ratios above 2.0 cannot occur; the shrink end is reachable
    (three consecutive rejections give 0.2**3) and is asserted more loosely,
    because the `(h^2/24) q'''` midpoint term does not cancel off a uniform grid --
    measured -12.2% at r=0.008.
    """
    h = 1e-9
    for ratio in (0.25, 0.5, 1.0, 2.0):
        r = _trap_est_over_true(h * ratio, h, h)
        assert abs(r - 1.0) < 0.10, \
            'est/true = %.4f at step ratio %g, outside 10%%' % (r, ratio)
    for ratio in (0.008, 0.05, 0.1):
        r = _trap_est_over_true(h * ratio, h, h)
        assert abs(r - 1.0) < 0.15, \
            'est/true = %.4f at deep-shrink ratio %g' % (r, ratio)


def test_gate_4g_b3_ywr_and_classic_now_coincide_for_trapezoidal():
    """4d resolved for trapezoidal by construction rather than by deletion.

    The mode-free branch does not read `lte_formula`, so the two selections return
    the same value bit for bit. The branches survive only on the one step of a run
    that has fewer than three past charges.
    """
    h = 1e-9
    for ratio in (0.25, 1.0, 2.0):
        a = _trap_est_over_true(h * ratio, h, h, formula='classic')
        b = _trap_est_over_true(h * ratio, h, h, formula='ywr')
        assert a == b, \
            "lte_formula still changes the trapezoidal estimate at ratio %g: %r vs %r" \
            % (ratio, a, b)


def test_gate_4g_b_the_fallback_is_used_for_exactly_one_step():
    """`h_last2 is None` must select the old two-past-point formula, not zeros.

    Returning zeros would make the second step of every run unchecked -- the defect
    stage 3 removed from the first step. And the fallback must be reachable: it is
    what runs before the ring buffer has three real charges.
    """
    from pycircuit.circuit.integrator import TrapezoidalIntegrator
    h = 1e-9
    t_n = 0.371e-6
    est, _p = TrapezoidalIntegrator('classic').compute_lte(
        q_curr=np.array([_4g_q(t_n)]), h_curr=h,
        q_last=np.array([[_4g_q(t_n - h)], [_4g_q(t_n - 2 * h)], [0.0]]),
        iq_last=np.array([[_4g_dq(t_n - h)], [_4g_dq(t_n - 2 * h)]]),
        h_last=h, is_first_step=False, toolkit=numeric, h_last2=None)
    assert float(est[0]) != 0.0, \
        'the h_last2=None fallback returns zero, making that step unchecked'
    assert np.isfinite(float(est[0]))


def test_both_second_order_methods_ask_for_three_past_charges():
    """The ring buffer must be sized for the ESTIMATOR, not just the method.

    Trapezoidal's `compute_derivatives` looks back one step and Gear-2's looks back
    two, but since 4i both estimate `q3` from a third divided difference of the
    charge, which needs `q_{n-3}`. If either of these drops back to the method's own
    depth the estimator silently reads a seeded initial charge as if it were
    history -- and produces a plausible wrong number rather than an error.

    Gear-2 was 2 until 4i. That is the churn 4i was expected to cause, and it is the
    whole reason the fix was available: `integrator.py` used to record "Gear-2 keeps
    just two past charges, so a third divided difference of q is not available at
    all" as the reason it had to difference the companion currents instead.
    """
    from pycircuit.circuit.integrator import (TrapezoidalIntegrator,
                                              Gear2Integrator, EulerIntegrator)
    assert TrapezoidalIntegrator().get_required_history() == 3
    assert Gear2Integrator().get_required_history() == 3
    assert EulerIntegrator().get_required_history() == 1


def test_gate_4g_b4_rejections_collapse_on_the_stiff_case():
    """End to end. Recorded before 4g(b): 757 rejections, 26 force-accepts, 923
    accepted steps. After: 23, 1 and 464.

    This is the test that would catch `_dt_last2` being rolled in the wrong order,
    which produces a plausible-looking but wrong grid rather than an error.
    """
    import warnings as _w
    from pycircuit.circuit.integrator import TrapezoidalIntegrator
    from pycircuit.circuit.stepcontroller import IntegralController
    from pycircuit.circuit.elements import L

    class _Counter(IntegralController):
        def __init__(self):
            super().__init__()
            self.rejections = 0

        def evaluate_step(self, *args, **kwargs):
            accept, h_next = super().evaluate_step(*args, **kwargs)
            if not accept:
                self.rejections += 1
            return accept, h_next

    cir = SubCircuit(toolkit=numeric)
    cir['C1'] = C(1, gnd, c=1e-6)
    cir['R1'] = R(1, 2, r=1.0)
    cir['L1'] = L(2, gnd, L=1e-6)
    x0 = np.zeros(cir.n)
    x0[cir.get_node_index('1')] = 1.0
    x0[cir.get_node_index('2')] = 1.0

    tran = Transient(cir, toolkit=numeric,
                     integrator=TrapezoidalIntegrator('ywr'), reltol=1e-5)
    ctrl = _Counter()
    tran.step_controller = ctrl
    with _w.catch_warnings():
        _w.simplefilter('ignore', RuntimeWarning)
        res = tran.solve(refnode=gnd, tend=5e-3, timestep=2e-4, x0=x0)

    ## 757 before; a 10x margin on the measured 23 leaves room for platform drift
    ## while still failing loudly if the mode comes back.
    assert ctrl.rejections < 230, \
        'trapezoidal rejected %d steps; it was 757 with the g-based estimator and ' \
        '23 with the mode-free one' % ctrl.rejections
    assert len(np.asarray(res.sweep_values)) < 700, \
        'step count %d; it was 923 before 4g(b) and 464 after' \
        % len(np.asarray(res.sweep_values))


def test_gate_4g_b_a_second_solve_is_identical_to_the_first():
    """The step history must be rebuilt per run, like the charge history.

    `_qlast` is re-seeded at the top of every `solve()`, so `_dt_last2` must be
    too. Otherwise a second run starts with a leftover step size from the first
    while `q_last[2]` is the freshly seeded initial charge, and the estimator
    differences a grid that never existed -- silently, and only on re-use.
    """
    import warnings as _w
    from pycircuit.circuit.integrator import TrapezoidalIntegrator
    tran = Transient(_pulsed_rc(), toolkit=numeric,
                     integrator=TrapezoidalIntegrator('classic'))
    with _w.catch_warnings():
        _w.simplefilter('ignore', RuntimeWarning)
        first = tran.solve(refnode=gnd, tend=3e-6, timestep=1e-7)
        v1 = np.asarray(first.v('n2'), dtype=float).copy()
        second = tran.solve(refnode=gnd, tend=3e-6, timestep=1e-7)
        v2 = np.asarray(second.v('n2'), dtype=float)

    assert len(v1) == len(v2), \
        'a re-solve took a different number of steps: %d then %d' % (len(v1), len(v2))
    assert np.array_equal(v1, v2), \
        're-solving the same circuit gave a different waveform; max delta %g' \
        % np.max(np.abs(v1 - v2))


# ---------------------------------------------------------------------------
# 4i -- one q3 estimator for both second-order methods
#
# Gear-2's LTE is -(1/6) h1 (h1+h2) q3 and trapezoidal's is -(1/6) h1^2 q3, so
# both need the same quantity and only the constant differs.  Estimating it from a
# third divided difference of the CHARGE is exact up to where the mean-value point
# sits, with no coefficient that depends on the step ratios -- which every earlier
# formulation in this file lacked.
# ---------------------------------------------------------------------------

def _gear2_true_lte(h_curr, h_last, t_n=0.371e-6):
    """The exact local truncation error: companion from exact charges, minus q'."""
    t1, t2 = t_n - h_curr, t_n - h_curr - h_last
    a0 = (2 * h_curr + h_last) / (h_curr * (h_curr + h_last))
    a1 = -(h_curr + h_last) / (h_curr * h_last)
    a2 = h_curr / (h_last * (h_curr + h_last))
    return (a0 * _4g_q(t_n) + a1 * _4g_q(t1) + a2 * _4g_q(t2)) - _4g_dq(t_n)


def _gear2_est_over_true(h_curr, h_last, h_last2, t_n=0.371e-6, formula='classic'):
    from pycircuit.circuit.integrator import Gear2Integrator
    t1, t2, t3 = t_n - h_curr, t_n - h_curr - h_last, t_n - h_curr - h_last - h_last2
    est, _p = Gear2Integrator(formula).compute_lte(
        q_curr=np.array([_4g_q(t_n)]), h_curr=h_curr,
        q_last=np.array([[_4g_q(t1)], [_4g_q(t2)], [_4g_q(t3)]]),
        ## Deliberately poisoned: the estimator must not read the companion
        ## currents at all once four charges exist.
        iq_last=np.array([[np.nan], [np.nan]]),
        h_last=h_last, is_first_step=False, toolkit=numeric, h_last2=h_last2)
    return float(est[0]) / _gear2_true_lte(h_curr, h_last, t_n)


## The ratios the controller can actually produce.  0.008 is three consecutive
## rejections (0.2**3); 2.414 is ZERO_STABILITY_RATIO, above which the integrator
## drops to order 1 and this estimator is not used.
_4I_RATIOS = (0.008, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 2.414)


def test_gate_4i1_gear2_is_flat_across_the_reachable_step_ratios():
    """Measured before 4i: 83.06 at ratio 0.008, falling to 0.695 at 4 -- a 119x
    spread on the SHIPPED DEFAULT integrator.

    The 83x is not hypothetical: ratio 0.008 is the step after three consecutive
    rejections, i.e. exactly when the controller has already collapsed the step and
    is deciding whether to collapse it further.
    """
    h = 1e-9
    for ratio in _4I_RATIOS:
        r = _gear2_est_over_true(h * ratio, h, h)
        assert abs(r - 1.0) < 0.05, \
            'gear2 est/true = %.4f at step ratio %g, outside 5%%' % (r, ratio)


def test_gate_4i3_trapezoidal_is_flat_across_the_reachable_step_ratios():
    """4g(b) got trapezoidal to a 1.26x spread by differencing `d`; 4i gets it to
    1.008x by differencing the charge.

    The residual 4g(b) had to declare -- `d_k` carries `(h_k^2/24) q3`, which
    cancels only on a uniform grid -- is gone, so the -12.2%/+12.9% that gate 4g-b2
    had to report as a partial failure is gone with it.
    """
    h = 1e-9
    for ratio in _4I_RATIOS:
        r = _trap_est_over_true(h * ratio, h, h)
        assert abs(r - 1.0) < 0.05, \
            'trap est/true = %.4f at step ratio %g, outside 5%%' % (r, ratio)


def test_gate_4i_the_estimator_never_reads_the_companion_currents():
    """Both second-order estimators, with `iq_last` poisoned to NaN.

    A finite result is proof that neither the parasitic mode (trapezoidal) nor the
    method's own truncation error (Gear-2) can reach the estimate -- which is the
    entire content of 4g(b) and 4i respectively.
    """
    h = 1e-9
    assert np.isfinite(_gear2_est_over_true(h, h, h)), \
        'the Gear2 estimator still reads iq_last'
    assert np.isfinite(_trap_est_over_true(h, h, h)), \
        'the trapezoidal estimator still reads iq_last'


def test_gate_4i2_gear2_is_consistent():
    """est/true must converge to 1 as h falls: 0.9442 / 0.9950 / 0.9995."""
    ratios = [_gear2_est_over_true(h, h, h) for h in (1e-8, 1e-9, 1e-10)]
    assert abs(ratios[2] - 1.0) < 0.02, \
        'gear2 est/true = %.4f at h=1e-10' % ratios[2]
    assert abs(ratios[2] - 1.0) < abs(ratios[0] - 1.0), \
        'gear2 est/true does not converge to 1 as h falls: %s' % ratios


@pytest.mark.parametrize('formula', ['classic', 'ywr'])
def test_gate_4i_lte_formula_no_longer_selects_anything(formula):
    """`lte_formula` is now vestigial for both second-order methods.

    The shared helper does not read it, so the two selections agree bit for bit
    once four charges exist. They still differ on the single fallback step, which
    is the only remaining effect the parameter has -- recorded rather than removed,
    because deleting a public parameter is the maintainer's call.
    """
    h = 1e-9
    for ratio in (0.05, 1.0, 2.0):
        a = _gear2_est_over_true(h * ratio, h, h, formula='classic')
        b = _gear2_est_over_true(h * ratio, h, h, formula='ywr')
        assert a == b, \
            'lte_formula still changes the Gear2 estimate at ratio %g: %r vs %r' \
            % (ratio, a, b)


def test_gate_4i5_the_fallback_is_finite_and_non_zero_for_both_methods():
    """`h_last2 is None` -- the one step of a run with fewer than four charges.

    Returning zeros there would make that step unchecked, which is the defect stage
    3 removed from the first step.
    """
    from pycircuit.circuit.integrator import (TrapezoidalIntegrator,
                                              Gear2Integrator)
    h = 1e-9
    t_n = 0.371e-6
    q_last = np.array([[_4g_q(t_n - h)], [_4g_q(t_n - 2 * h)], [0.0]])
    iq_last = np.array([[_4g_dq(t_n - h)], [_4g_dq(t_n - 2 * h)]])
    for integ in (TrapezoidalIntegrator('classic'), TrapezoidalIntegrator('ywr'),
                  Gear2Integrator('classic'), Gear2Integrator('ywr')):
        est, _p = integ.compute_lte(
            q_curr=np.array([_4g_q(t_n)]), h_curr=h, q_last=q_last,
            iq_last=iq_last, h_last=h, is_first_step=False, toolkit=numeric,
            h_last2=None)
        v = float(est[0])
        assert v != 0.0 and np.isfinite(v), \
            '%s fallback returned %r' % (type(integ).__name__, v)


def test_gate_4i_the_error_constants_come_from_the_derivation():
    """The two constants are the whole difference between the integrators.

    Gear-2: -(1/6) h1 (h1+h2) q3.  Trapezoidal: -(1/6) h1^2 q3.  On a uniform grid
    that makes Gear-2's estimate exactly twice trapezoidal's, and getting either
    constant wrong would bias every step of every run by a fixed factor while still
    looking perfectly flat in the ratio sweeps above.
    """
    h = 1e-9
    t_n = 0.371e-6
    from pycircuit.circuit.integrator import (TrapezoidalIntegrator,
                                              Gear2Integrator)
    kw = dict(q_curr=np.array([_4g_q(t_n)]), h_curr=h,
              q_last=np.array([[_4g_q(t_n - h)], [_4g_q(t_n - 2 * h)],
                               [_4g_q(t_n - 3 * h)]]),
              iq_last=np.array([[np.nan], [np.nan]]),
              h_last=h, is_first_step=False, toolkit=numeric, h_last2=h)
    g = float(Gear2Integrator('classic').compute_lte(**kw)[0][0])
    t = float(TrapezoidalIntegrator('classic').compute_lte(**kw)[0][0])
    assert g / t == pytest.approx(2.0, rel=1e-12), \
        'on a uniform grid Gear-2 must estimate exactly 2x trapezoidal, got %.9f' \
        % (g / t)


# ---------------------------------------------------------------------------
# 4h -- `fixed_timestep=True` must actually fix the timestep
#
# `dt` is loop-carried and breakpoint truncation overwrote it, while the restore
# at the bottom of the loop was guarded by `if not fixed_timestep`. Each later
# breakpoint then truncated the already-shrunken step again, so it collapsed
# geometrically.
# ---------------------------------------------------------------------------

def _fixed_run(tend=3e-6, timestep=1e-7, **kwargs):
    """A VPulse-driven fixed-step run; returns the time vector."""
    import warnings as _w
    tran = Transient(_pulsed_rc(), toolkit=numeric, **kwargs)
    with _w.catch_warnings():
        _w.simplefilter('ignore', RuntimeWarning)
        res = tran.solve(refnode=gnd, tend=tend, timestep=timestep,
                         fixed_timestep=True)
    return np.asarray(res.sweep_values, dtype=float)


def test_gate_4h1_fixed_timestep_takes_the_requested_number_of_steps():
    """Measured before: 292 steps against an expected 30, ending at dt=1.24e-19.

    The plan recorded 19 002 steps on a different circuit; the mechanism is the
    same and the magnitude depends only on how many breakpoints the drive has.
    """
    tend, timestep = 3e-6, 1e-7
    t = _fixed_run(tend, timestep)
    expected = round(tend / timestep)
    assert abs(len(t) - expected) <= 1, \
        'fixed_timestep took %d steps for tend/timestep = %d' % (len(t), expected)


def test_gate_4h2_the_fixed_grid_is_actually_uniform():
    """A step count can be right while the spacing is not.

    Declared separately for that reason: the `VSin` case that motivated this had a
    step count inside +-1 *and* a final step of 2.033e-20 s, because a uniform grid
    divides `tend` exactly and `dt = tend - t` turns the floating-point residue
    into a step. Every step must be `timestep`, and the last must be a real step
    rather than rounding noise.
    """
    timestep = 1e-7
    t = _fixed_run(3e-6, timestep)
    dt = np.diff(t)
    assert np.allclose(dt, timestep, rtol=1e-9, atol=1e-9 * timestep), \
        'fixed grid is not uniform: dt ranges %.6g .. %.6g against timestep %g' \
        % (dt.min(), dt.max(), timestep)
    assert dt.min() > 0.5 * timestep, \
        'a degenerate step of %.4g s survived in a %g s grid' % (dt.min(), timestep)


def test_gate_4h3_a_breakpoint_inside_a_fixed_step_still_drops_the_order():
    """Uniformity must not cost the discontinuity handling.

    The grid no longer moves for a breakpoint, but crossing one must still drop to
    order 1 for the next step so no 2nd-order polynomial is fitted across the
    corner. Instrumented rather than inspected: `_effective_method` is derived from
    the live integrator object, so this sees what actually ran.
    """
    import warnings as _w
    from pycircuit.circuit.integrator import Gear2Integrator
    tran = Transient(_pulsed_rc(), toolkit=numeric, integrator=Gear2Integrator())
    seen = []
    real = tran.get_diff

    def spy(q, C_, method=None):
        out = real(q, C_, method)
        seen.append(tran._effective_method)
        return out

    tran.get_diff = spy
    with _w.catch_warnings():
        _w.simplefilter('ignore', RuntimeWarning)
        tran.solve(refnode=gnd, tend=3e-6, timestep=1e-7, fixed_timestep=True)

    drops = sum(1 for m in seen if m == 'EulerIntegrator')
    assert drops > 0, \
        'no order drop in a fixed run over a VPulse: breakpoints have stopped ' \
        'being noticed entirely'
    assert drops < len(seen), \
        'every evaluation dropped order; the run is first-order throughout'


def test_gate_4h4_a_fixed_grid_that_cannot_be_honoured_warns():
    """Shrinking on non-convergence is the only way to make progress -- but under
    `fixed_timestep` it silently abandons the grid the caller asked for.

    Same failure class as 4b's force-accept, and the same treatment: the result is
    still returned, and without the warning the caller has no way to learn that
    part of their uniform grid is not uniform.
    """
    import warnings as _w
    tran = Transient(_pulsed_rc(), toolkit=numeric)
    real = tran.solve_timestep
    calls = {'n': 0}

    def flaky(*args, **kwargs):
        calls['n'] += 1
        if calls['n'] == 5:
            raise NoConvergenceError('injected')
        return real(*args, **kwargs)

    tran.solve_timestep = flaky
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter('always')
        tran.solve(refnode=gnd, tend=3e-6, timestep=1e-7, fixed_timestep=True)

    hits = [c for c in caught
            if issubclass(c.category, RuntimeWarning)
            and 'no longer uniform' in str(c.message)]
    assert len(hits) == 1, \
        'expected one fixed-grid fallback warning, got %d' % len(hits)
    msg = str(hits[0].message)
    assert 't=' in msg and 'falling back to' in msg, \
        'the warning must name the time and the step it fell back to: %s' % msg


def test_gate_4h_the_adaptive_path_is_untouched():
    """The degenerate-final-step guard is scoped to fixed-step, deliberately.

    Measured before the change: zero degenerate steps on the adaptive path for both
    circuits, because a controller-chosen `dt` does not land on `tend` to within
    1e-20. Scoping the guard keeps 4h attributable -- if the adaptive step count
    moves, something other than 4h did it.
    """
    import warnings as _w
    tran = Transient(_pulsed_rc(), toolkit=numeric)
    with _w.catch_warnings():
        _w.simplefilter('ignore', RuntimeWarning)
        res = tran.solve(refnode=gnd, tend=3e-6, timestep=1e-7)
    dt = np.diff(np.asarray(res.sweep_values, dtype=float))
    assert (dt < 1e-6 * np.median(dt)).sum() == 0, \
        'the adaptive path grew a degenerate step; 4h should not have touched it'
