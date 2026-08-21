"""STAGE 12B -- `coupled_method`, selecting how the coupled path corrects `h`.

    'approx'   Fang sec. 3.4 -- the step comes from the error RATIO (eq 17) and
               the solution is corrected by eq (18). The default.
    'bordered' Fang eq (12)/(14) -- a linearised Newton step on the LTE equation,
               which additionally accounts for the pending solution update
               through `q^T dv`.

**Eq (12) is usable only because its denominator is computed analytically.** The
paper forms it as `q^T dxh + d`, and those two terms are how the solution moves
with the step size and how the extrapolation moves; both are approximately
`dv/dt`, so the difference is the derivative of the truncation error -- tiny by
construction and computed as a difference of two large numbers.

Measured at h = 3.48e-7 on a driven RC, against a ground truth obtained by
RE-SOLVING the circuit at perturbed `h` (`benchmarks/transient_review/`):

    ground truth   +4.678e6
    analytic       +4.392e6   ratio  0.939
    subtraction    -9.680e5   ratio -0.207   <- wrong SIGN

The subtracted terms are -1.310e8 and +1.300e8; the result is 0.74% of the
larger. That wrong sign is why an earlier attempt at eq (12) drove the step size
down four decades while the error sat far below its band.
"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VSin
from pycircuit.circuit.transient import Transient

TAU = 1e-4
W = 2 * np.pi * 1e3


def _rc():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c


def _analytic(t):
    A = 1.0 / np.sqrt(1.0 + (W * TAU) ** 2)
    phi = np.arctan(W * TAU)
    return A * (np.sin(W * t - phi) + np.sin(phi) * np.exp(-t / TAU))


def _run(method=None):
    tran = Transient(_rc(), toolkit=numeric, reltol=1e-5)
    if method is not None:
        tran.par.coupled_method = method
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=5e-4, timestep=1e-5, coupled_lte=True)
    t = np.asarray(res.v('b').x, dtype=float).ravel()
    v = np.asarray(res.v('b').y, dtype=float).ravel()
    return tran.statistics, float(np.max(np.abs(v - _analytic(t))[2:]))


def test_the_default_is_the_method_with_the_measured_record():
    """`approx` is sec. 3.4, and it is what every stage-12 number was taken on."""
    tran = Transient(_rc(), toolkit=numeric)
    assert tran.par.coupled_method == 'approx'


def test_both_methods_solve_the_circuit_and_take_no_rejections():
    """Figure 3 has no rejection branch, whichever correction is selected."""
    for method in ('approx', 'bordered'):
        st, err = _run(method)
        assert st.rejected_steps == 0, '%s took %d rejections' % (method, st.rejected_steps)
        assert st.accepted_steps > 10
        assert err < 5e-3, '%s: max error %g' % (method, err)


def test_the_two_methods_agree_on_a_smooth_circuit():
    """On a smooth drive the extra `q^T dv` term changes almost nothing.

    That is the expected result, not a disappointing one: the term accounts for
    the LTE moving as the solution converges, and on a smooth circuit the
    solution is already close when the LTE is evaluated. Pinned so that a future
    change making them diverge here is noticed.
    """
    st_a, err_a = _run('approx')
    st_b, err_b = _run('bordered')
    assert abs(st_b.accepted_steps - st_a.accepted_steps) <= 0.05 * st_a.accepted_steps
    assert err_b == pytest.approx(err_a, rel=0.05)


def test_an_unknown_method_is_refused():
    """A typo must not silently select the default."""
    tran = Transient(_rc(), toolkit=numeric)
    tran.par.coupled_method = 'schur'
    with pytest.raises(ValueError) as exc:
        tran.solve(tend=5e-5, timestep=1e-5, coupled_lte=True)
    assert 'schur' in str(exc.value)


def _pulsed_rc():
    from pycircuit.circuit.elements import VPulse
    c = SubCircuit()
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-5, tr=1e-6, tf=1e-6,
                     pw=2e-5, per=5e-5)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    return c


def _pulse_run(method):
    tran = Transient(_pulsed_rc(), toolkit=numeric, reltol=1e-5)
    tran.par.coupled_method = method
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tran.solve(tend=6e-5, timestep=1e-6, coupled_lte=True)
    return tran.statistics


def test_bordered_grows_the_step_back_where_the_error_is_zero():
    """THE DEFECT THAT MADE EQ (12) LOOK UNUSABLE, and it was not eq (12).

    `bordered`'s denominator is `err * w'(h)/w(h)`, which vanishes when the error
    does -- and on a pulsed circuit the error IS zero across the flat regions
    between edges, where the solution is constant and the extrapolation
    reproduces it exactly. Measured on this circuit: **76.1% of all step
    adjustments happen at err = 0**.

    Guarding that degenerate denominator by leaving `h` alone is the obvious
    reading and the wrong one. Zero error does not mean "leave the step alone",
    it means "the step is far too small". With the guard, once an edge forced `h`
    down it never grew back: 11831 of 12382 time points repeated the step before
    them, the median step came out 10x smaller than `approx`'s, and the run took
    5.6x the time points for the same waveform.

    Asserted against `approx` rather than against an absolute count, so the test
    survives step counts changing for unrelated reasons.
    """
    st_a = _pulse_run('approx')
    st_b = _pulse_run('bordered')
    assert st_b.accepted_steps < 2.0 * st_a.accepted_steps, (
        'bordered took %d time points against approx\'s %d -- the step is not '
        'growing back where the error is zero'
        % (st_b.accepted_steps, st_a.accepted_steps))
    ## And the Newton work must be comparable too: a step count held down by
    ## doing more iterations per point would pass the check above for the wrong
    ## reason.
    assert st_b.newton_iterations < 2.0 * st_a.newton_iterations


def test_neither_method_rejects_a_step_on_a_pulsed_circuit():
    """Figure 3's promise has to survive real discontinuities, not just smooth
    drives -- which is the whole reason a breakpoint circuit is in this file.

    RE-DERIVED at F13 (doc/transient_review_260820.md): the old `== 0`
    passed vacuously -- the coupled path never counted rejections at all.
    With the counter live, the honest statement of Figure 3 is narrower:
    the steps the method SOLVES for take no rejections (the smooth test
    above asserts exactly 0), while HELD steps -- breakpoint- or
    tend-truncated, whose size was never the method's to choose -- may
    retry when their imposed size fails the error test.  Measured at
    re-derivation: 16/985 and 14/1035 rejections, all at edges.  Bound
    them to a small fraction rather than pretending they are zero."""
    for method in ('approx', 'bordered'):
        st = _pulse_run(method)
        assert st.rejected_steps <= 0.05 * st.accepted_steps, \
            '%s: %d rejections against %d accepted -- held-step retries ' \
            'should be a few percent at worst' \
            % (method, st.rejected_steps, st.accepted_steps)
        assert st.breakpoints_hit > 0, '%s never hit an edge' % method
