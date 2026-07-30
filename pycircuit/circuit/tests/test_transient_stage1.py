"""Stage 1 of `doc/transient_work_plan.md`: the silent failures.

Every test here corresponds to a declared gate. They exist because each defect
they cover produced a **confident wrong answer** rather than an error, which is
the one class of failure a test suite cannot be trusted to catch by accident.

Gate 1-1  a failed operating point raises, and names the escape
Gate 1-2  epar reaches the devices, in the inner DC and in every transient step
Gate 1-3  `bypass` is connected to something
Gate 1-5  Newton's tolerance and the controller's are separate knobs
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
    """`vabstol` is Newton's; `lte_vabstol` is the step controller's."""
    tran = Transient(_rc(), toolkit=numeric)

    ## Newton's, and it agrees with DC's rather than being 10^6 looser.
    from pycircuit.circuit.dcanalysis import DC
    dc = DC(_rc(), toolkit=numeric)
    assert tran.par.vabstol == dc.par.vabstol, \
        'the operating point and the steps after it must be solved to the same ' \
        'tolerance; got transient %g vs DC %g' % (tran.par.vabstol, dc.par.vabstol)

    ## The controller's, which is the one the step-count result belongs to.
    assert tran.par.lte_vabstol == 1e-6
    assert tran.par.lte_iabstol == 1e-12


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


def test_gate_1_5_lte_vabstol_moves_the_step_count():
    """The controller's tolerance must control the controller."""
    loose, _ = _pulsed_run(lte_vabstol=1e-3)
    base, _ = _pulsed_run()
    tight, _ = _pulsed_run(lte_vabstol=1e-9)

    assert loose < base < tight, \
        'step count must fall monotonically as the LTE tolerance loosens; ' \
        'got %d (1e-3) / %d (default 1e-6) / %d (1e-9)' % (loose, base, tight)


def test_gate_1_5_vabstol_does_not_move_the_step_count():
    """Newton's tolerance must NOT be a step-control knob.

    This is the property the shared parameter made impossible: there was no way to
    change one role without silently changing the other. A failure here means
    `vabstol` is feeding the step controller again.
    """
    base_steps, base_v = _pulsed_run()
    for value in (1e-9, 1e-14):
        steps, v = _pulsed_run(vabstol=value)
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
