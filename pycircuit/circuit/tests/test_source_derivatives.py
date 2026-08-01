"""STAGE 12B -- ``TimeFunction.dfdt``, needed for Fang's ``p = df_ckt/dh``.

The circuit residual is evaluated at ``t_{m-1} + h``, so once ``h`` is an unknown
(DAC 2013 eq 11) the source term contributes ``du/dt`` to ``p`` -- and on a driven
circuit that is often the largest part of it.  Every source therefore needs a time
derivative, and every one of them is checked here against a central difference of
its OWN ``f``, so a derivative can never quietly disagree with the waveform it is
supposed to differentiate.

Kinks are tested separately and deliberately: `Pulse` and `PWL` are piecewise
linear, so ``du/dt`` genuinely does not exist at their corners, and the choice
made there (right-hand limit) is a decision rather than a derivation.
"""
import numpy as np
import pytest

from pycircuit.circuit.func import Sin, Pulse, PWL, Exp, AM, SFFM, TimeFunction


SMOOTH_CASES = [
    ('Sin', Sin(offset=0.2, amplitude=1.5, freq=1e3, td=0.0, phase=30.0),
     [1e-5, 1e-4, 3.3e-4, 9e-4]),
    ('Sin damped', Sin(offset=0.0, amplitude=1.0, freq=2e3, td=1e-4, theta=3e3),
     [2e-4, 5e-4, 8e-4]),
    ('Exp', Exp(v1=0.0, v2=2.0, td1=1e-5, tau1=2e-5, td2=8e-5, tau2=3e-5),
     [3e-5, 5e-5, 9e-5, 2e-4]),
    ('AM', AM(vo=0.1, va=1.0, fc=1e4, fm=1e3, m=0.5),
     [1e-5, 5e-5, 2e-4]),
    ('SFFM', SFFM(vo=0.0, va=1.0, fc=1e4, mdi=2.0, fm=1e3),
     [1e-5, 5e-5, 2e-4]),
]


@pytest.mark.parametrize('name,src,times', SMOOTH_CASES,
                         ids=[c[0] for c in SMOOTH_CASES])
def test_dfdt_matches_a_central_difference_of_its_own_f(name, src, times):
    for t in times:
        eps = 1e-7 * max(t, 1e-6)
        fd = (src.f(t + eps) - src.f(t - eps)) / (2 * eps)
        got = float(src.dfdt(t))
        scale = max(abs(fd), abs(got), 1e-9)
        assert abs(got - fd) / scale < 1e-4, \
            '%s at t=%g: analytic %g, finite difference %g' % (name, t, got, fd)


def test_sin_derivative_is_zero_before_its_delay():
    """`f` holds the source at its offset before `td`, so `dfdt` must be flat.

    Applying the chain rule to the unclamped expression here would report a slope
    the source does not have -- the same class of defect stage 8(b) fixed in `f`
    itself, where the damping term grew backwards in time.
    """
    src = Sin(offset=0.7, amplitude=1.0, freq=1e3, td=1e-4, theta=5e3)
    for t in (0.0, 2e-5, 9.9e-5):
        assert float(src.dfdt(t)) == 0.0
    ## And nonzero just after, so the test is not passing on a dead function.
    assert abs(float(src.dfdt(1.1e-4))) > 1.0


def test_pulse_derivative_is_the_segment_slope():
    """0 / rise / 0 / fall / 0, with the slopes SPICE would give."""
    src = Pulse(v1=0.0, v2=1.0, td=1e-5, tr=1e-7, tf=2e-7, pw=2e-5, per=5e-5)
    assert float(src.dfdt(5e-6)) == 0.0                       # before td
    assert float(src.dfdt(1e-5 + 5e-8)) == pytest.approx(1.0 / 1e-7)   # rising
    assert float(src.dfdt(1.5e-5)) == 0.0                     # flat top
    assert float(src.dfdt(1e-5 + 1e-7 + 2e-5 + 1e-7)) == \
        pytest.approx(-1.0 / 2e-7)                            # falling
    assert float(src.dfdt(4.5e-5)) == 0.0                     # after the fall


def test_pulse_derivative_follows_the_period_fold():
    """The second period must look like the first, not like a flat tail."""
    src = Pulse(v1=0.0, v2=1.0, td=1e-5, tr=1e-7, tf=2e-7, pw=2e-5, per=5e-5)
    for k in (0, 1, 3):
        t = k * 5e-5 + 1e-5 + 5e-8
        assert float(src.dfdt(t)) == pytest.approx(1.0 / 1e-7), \
            'period %d did not fold' % k


def test_pwl_derivative_is_the_slope_of_the_segment_being_entered():
    src = PWL([0.0, 0.0, 1e-5, 2.0, 3e-5, 2.0, 4e-5, -1.0])
    assert float(src.dfdt(5e-6)) == pytest.approx(2.0 / 1e-5)
    assert float(src.dfdt(2e-5)) == 0.0
    assert float(src.dfdt(3.5e-5)) == pytest.approx(-3.0 / 1e-5)
    ## Flat outside the table, both ends.
    assert float(src.dfdt(-1e-6)) == 0.0
    assert float(src.dfdt(5e-5)) == 0.0


def test_at_a_corner_the_right_hand_limit_is_taken():
    """A DECISION, not a derivation: `du/dt` does not exist at a kink.

    `p` says how the residual moves when the step grows away from `t`, so the
    segment ahead is the one that governs. Pinned because the opposite choice is
    equally defensible in isolation and would silently change `p`'s sign at every
    pulse edge.
    """
    src = PWL([0.0, 0.0, 1e-5, 2.0, 3e-5, 2.0])
    ## Exactly on the corner at 1e-5: the segment ahead is the flat one.
    assert float(src.dfdt(1e-5)) == 0.0

    pulse = Pulse(v1=0.0, v2=1.0, td=1e-5, tr=1e-7, tf=1e-7, pw=2e-5, per=0)
    ## Exactly at td: the rise is ahead.
    assert float(pulse.dfdt(1e-5)) == pytest.approx(1.0 / 1e-7)


def test_zero_rise_time_reports_zero_rather_than_an_invented_slope():
    """`tr=0` is SPICE's default and makes the edge a true discontinuity.

    No finite derivative exists, so `p` must not be handed the enormous number a
    finite difference would produce. The edge is a breakpoint, and
    breakpoint-clamped steps are excluded from the coupled solve.
    """
    src = Pulse(v1=0.0, v2=1.0, td=1e-5, tr=0.0, tf=0.0, pw=2e-5, per=5e-5)
    assert np.isfinite(float(src.dfdt(1e-5)))
    assert float(src.dfdt(1e-5)) == 0.0


def test_base_class_falls_back_to_a_finite_difference():
    """A source that does not override `dfdt` must still work."""

    class Quadratic(TimeFunction):
        def f(self, t):
            return 3.0 + 2.0 * t + 5.0 * t ** 2

    src = Quadratic()
    for t in (0.0, 1e-3, 0.5):
        assert float(src.dfdt(t)) == pytest.approx(2.0 + 10.0 * t, rel=1e-5,
                                                   abs=1e-6)
