"""STAGE 12A -- Fang's LTE acceptance band (eq 15) and step damper (eq 16).

Two of these tests pin defects that the stage's own gates found in the first
implementation of the band, and both were invisible as exceptions -- they showed
up only as step counts:

  * aiming the step-size prediction AT `gamma_min` put the aim exactly on the
    rejection edge, so every undershoot rejected: 3172 rejections to accept 1187
    steps, to save 7.8%;
  * applying the damper to the REJECTION path capped how fast an over-tolerance
    step could retreat, so the stiff RLC ringdown exhausted its rejection budget,
    force-accepted, and crossed the whole ringing transient in 62 steps against
    the baseline's 490 -- while reporting an LTE of exactly zero, because by then
    it was integrating a signal that had already decayed.

The end-to-end test at the bottom is the one that would have caught the second
directly; the unit tests above it pin the arithmetic that caused both.
"""
import numpy as np
import pytest

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.elements import R, C, L
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.stepcontroller import IntegralController


SAFETY = 0.9      # must match the controller's own constant
P = 2             # (order+1) for the second-order integrators


def test_defaults_are_the_historical_one_sided_test():
    """Stage 12 is behind a flag until 12D, so the band must ship inert.

    Asserted on the aim point rather than on the attribute values alone: the aim
    is what actually decides step size, and `safety**p` is the fixed point the
    unbanded law converges to.
    """
    c = IntegralController()
    assert c.lte_gamma_min == 0.0
    assert c.lte_gamma_max == 1.0
    assert c.lte_eta is None
    assert c._band_target(SAFETY, P) == pytest.approx(SAFETY ** P, rel=1e-15)
    ## And the damper must be a no-op, not a wide clamp.
    assert c._damp(1e-3, 1e-9) == 1e-3


def test_a_band_containing_the_aim_point_does_not_move_it():
    """The paper's own 0.7/3.0 straddles this controller's 0.81 aim.

    That the band is then inert is the honest behaviour and the reason the
    stage-12 entry measurement predicted little gain -- it must not be disguised
    by re-aiming at the band centre unconditionally.
    """
    c = IntegralController().set_lte_band(0.7, 3.0)
    assert c._band_target(SAFETY, P) == pytest.approx(SAFETY ** P, rel=1e-15)


def test_aim_point_moves_inside_the_band_not_onto_its_edge():
    """The 3172-rejection defect, pinned.

    When the band excludes the natural aim the target has to move, but landing it
    ON `gamma_min` makes every undershoot a rejection.  A predict-then-test scheme
    must aim strictly inside; Fang does not need this because there `h` is solved
    to satisfy the band rather than tested against it.
    """
    lo, hi = 0.95, 1.0
    c = IntegralController().set_lte_band(lo, hi)
    target = c._band_target(SAFETY, P)
    assert lo < target < hi, \
        'aim %r must be strictly inside the band (%r, %r)' % (target, lo, hi)
    ## Specifically the geometric centre -- the point furthest from both edges in
    ## the ratio sense, which is the sense the step-size law works in.
    assert target == pytest.approx((lo * hi) ** 0.5, rel=1e-15)


def test_upper_bound_only_keeps_a_safety_margin():
    """With no lower edge there is nothing to be centred against."""
    c = IntegralController().set_lte_band(0.0, 3.0)
    ## 0.81 is inside [0, 3], so the aim is untouched.
    assert c._band_target(SAFETY, P) == pytest.approx(SAFETY ** P, rel=1e-15)
    ## But an aim above the ceiling is pulled under it with the usual margin.
    c2 = IntegralController().set_lte_band(0.0, 0.5)
    assert c2._band_target(SAFETY, P) == pytest.approx(0.5 * SAFETY, rel=1e-15)


def test_damper_limits_the_change_in_both_directions():
    c = IntegralController().set_lte_band(eta=0.15)
    h = 1e-6
    assert c._damp(10 * h, h) == pytest.approx(1.15 * h, rel=1e-15)
    assert c._damp(0.01 * h, h) == pytest.approx(0.85 * h, rel=1e-15)
    assert c._damp(1.05 * h, h) == pytest.approx(1.05 * h, rel=1e-15)


def test_band_validation_rejects_an_inverted_or_empty_band():
    with pytest.raises(ValueError):
        IntegralController().set_lte_band(1.0, 1.0)
    with pytest.raises(ValueError):
        IntegralController().set_lte_band(2.0, 1.0)
    with pytest.raises(ValueError):
        IntegralController().set_lte_band(-0.1, 1.0)
    with pytest.raises(ValueError):
        IntegralController().set_lte_band(eta=0.0)


def _stiff_rlc():
    """The ringdown the damper defect crossed in 62 steps."""
    c = SubCircuit()
    c['C1'] = C(1, gnd, c=1e-6)
    c['R1'] = R(1, 2, r=1.0)
    c['L1'] = L(2, gnd, L=1e-6)
    x0 = np.zeros(c.n)
    x0[c.get_node_index('1')] = 1.0
    x0[c.get_node_index('2')] = 1.0
    return c, dict(tend=5e-3, timestep=2e-4, x0=x0)


def _stiff_rlc_analytic(t, r=1.0, l=1e-6, c=1e-6):
    """``LC v'' + RC v' + v = 0`` with ``v(0)=1``, ``v'(0)=0``."""
    alpha = r / (2.0 * l)
    wd = np.sqrt(1.0 / (l * c) - alpha ** 2)
    return np.exp(-alpha * t) * (np.cos(wd * t) + (alpha / wd) * np.sin(wd * t))


def _run(band):
    cir, kw = _stiff_rlc()
    tran = Transient(cir, toolkit=numeric, reltol=1e-5, **band)
    res = tran.solve(coupled_lte=False, **kw)
    w = res.v('1')
    t = np.asarray(w.x, dtype=float).ravel()
    v = np.asarray(w.y, dtype=float).ravel()
    ## Skip the opening step: it is accepted unevaluated, so no tolerance governs
    ## it and it dominates the maximum on this circuit regardless of settings.
    err = float(np.max(np.abs(v - _stiff_rlc_analytic(t))[2:]))
    return tran.statistics.accepted_steps, err, float(t[-1]), kw['tend']


def test_damper_does_not_muzzle_the_rejection_path():
    """The 62-vs-490 defect, end to end.

    Eq (16) bounds how far one ACCEPTED step may sit from the one before it. It
    is not a limit on how fast a step that failed its error test may retreat, and
    applying it there leaves the error control unable to respond on a stiff
    circuit.  Checked on both step count and accuracy: the broken version got its
    step count DOWN, which is why a step-count-only check would have read as a
    win.
    """
    base_steps, base_err, base_t, tend = _run({})
    damped_steps, damped_err, damped_t, _ = _run(dict(lte_eta=0.15))

    ## Both runs must actually reach the end -- a truncated run reports a
    ## flattering step count for the worst possible reason.
    assert base_t >= tend * (1 - 1e-6)
    assert damped_t >= tend * (1 - 1e-6)

    ## The damper may cost steps (it limits growth); it must not save them by
    ## failing to resolve the transient.  The defect gave 0.13x here.
    assert damped_steps >= 0.8 * base_steps, \
        'damper cut steps %d -> %d: it is throttling error control, not ' \
        'smoothing step changes' % (base_steps, damped_steps)
    assert damped_err == pytest.approx(base_err, rel=0.05)


def test_band_defaults_leave_the_run_bit_identical():
    """The default path must not move while stage 12 is behind its flag."""
    steps_default, err_default, _, _ = _run({})
    steps_explicit, err_explicit, _, _ = _run(
        dict(lte_gamma_min=0.0, lte_gamma_max=1.0, lte_eta=None))
    assert steps_explicit == steps_default
    assert err_explicit == err_default
