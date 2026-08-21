"""Unit gates for stepcontroller.py's dark territory.

The Phase-0 coverage snapshot (doc/phase0_baseline_260821.md) flagged
this module at 72% with the PI controller and parts of the band logic
unexercised; re-measured after the campaign the dark regions were: both
controllers' LTE-solve fallback (the loud units-violation warning kept
by decision 0.3d), IntegralController's lower-band growth-retry with its
three suppression guards, PIController's gamma_min refusal (F10), and
the WHOLE of SolutionLTEController.evaluate_step -- a selectable
standard-path controller (``tran.step_controller = ...``) whose design
record cites measurements no test was pinning.  These gates drive each
region directly, stubbing the integrator, exactly so the behaviors the
comments claim stay true.
"""
import numpy as np
import pytest

from pycircuit.circuit import numeric
from pycircuit.circuit.stepcontroller import (IntegralController,
                                              PIController,
                                              SolutionLTEController,
                                              MAX_GROWTH_RATIO,
                                              MIN_SHRINK_RATIO)


class StubIntegrator:
    """compute_lte returns a fixed charge-error vector; ORDER for the
    solution controller's degree selection."""
    def __init__(self, Eg, p=3, order=1):
        self.Eg = np.asarray(Eg, float)
        self.p = p
        self.ORDER = order

    def compute_lte(self, **kw):
        return self.Eg, self.p


def _evaluate(ctrl, integ, J, x_curr=None, h_curr=1e-6, max_step=1e-3,
              h_clamped=False, x_hist=None, h_last=1e-6, h_last2=None,
              no_history=False):
    n = len(J)
    x_curr = np.ones(n) if x_curr is None else np.asarray(x_curr, float)
    return ctrl.evaluate_step(
        x_curr=x_curr, x_last=np.zeros(n), q_curr=np.zeros(n),
        q_last_hist=[np.zeros(n)], iq_last_hist=[np.zeros(n)],
        h_curr=h_curr, h_last=h_last, no_history=no_history,
        J=np.asarray(J, float), active_integrator=integ, irefnode=n - 1,
        reltol=1e-3, abstol=1e-12, toolkit=numeric, max_step=max_step,
        TRTOL=7.0, n_nodes=None, h_last2=h_last2, h_clamped=h_clamped,
        x_hist=x_hist)


@pytest.mark.parametrize('ctrl_cls', [IntegralController, PIController])
def test_lte_solve_fallback_warns_and_keeps_running(ctrl_cls):
    """Decision 0.3d's kept-but-loud fallback: a singular reduced
    Jacobian makes the J^-1*Eg solve fail, the controller WARNS that the
    charge-domain residual is standing in for the solution-domain error
    (a units violation -- a current where a voltage belongs), and the
    step verdict is still produced from that stand-in.  Gate 4-D
    measured this firing zero times in production; this gate keeps the
    mechanism honest for the day it does fire."""
    ctrl = ctrl_cls()
    if ctrl_cls is PIController:
        ctrl.last_err = 0.5          # PI history seed (normally no_history)
    integ = StubIntegrator(Eg=[1e-9, 0.0])
    with pytest.warns(RuntimeWarning, match='LTE solve failed'):
        accept, h_next = _evaluate(ctrl, integ, J=np.zeros((2, 2)))
    ## The tiny stand-in residual is far under tolerance: accepted, and
    ## the error the fallback computed is observable.
    assert accept
    assert ctrl.last_err is not None and ctrl.last_err < 1.0


def test_integral_lower_band_growth_retry_and_its_three_guards():
    """Eq (15)'s lower bound in IntegralController: a step far under
    gamma_min is REDONE larger at the same time point -- and the three
    suppressions the comments claim are each honoured: a clamped step
    (breakpoint/tend truncation) is never grown, a step already at
    max_step has nowhere to grow, and a damper that would prevent any
    ACTUAL growth turns the retry into a plain accept (the anti-livelock
    'never report a grow that does not grow')."""
    integ = StubIntegrator(Eg=[1e-30, 0.0])   # err ~ 0: far under any band

    ctrl = IntegralController().set_lte_band(gamma_min=0.5, gamma_max=3.0)
    accept, h_next = _evaluate(ctrl, integ, J=np.eye(2))
    assert not accept                     # too-accurate: growth retry
    assert h_next > 1e-6                  # ...and it actually grows
    assert h_next <= 1e-6 * MAX_GROWTH_RATIO * (1 + 1e-12)

    ## Guard 1: the step size was never the controller's to choose.
    accept, _ = _evaluate(ctrl, integ, J=np.eye(2), h_clamped=True)
    assert accept

    ## Guard 2: already at max_step -- nowhere to grow.
    accept, _ = _evaluate(ctrl, integ, J=np.eye(2), h_curr=1e-3,
                          max_step=1e-3)
    assert accept

    ## Guard 3: eta so tight the damped retry would not actually grow.
    ctrl = IntegralController().set_lte_band(gamma_min=0.5, gamma_max=3.0,
                                             eta=1e-12)
    accept, _ = _evaluate(ctrl, integ, J=np.eye(2))
    assert accept


def test_pi_refuses_the_lower_band_and_passes_the_sentinel():
    """F10: PIController honours gamma_max/eta but REFUSES gamma_min > 0
    (the growth-retry redo has undefined PI history semantics), after
    the base class resolves the 'auto' sentinel so shipped defaults
    pass."""
    with pytest.raises(NotImplementedError, match='lower LTE band'):
        PIController().set_lte_band(gamma_min=0.5, gamma_max=3.0)
    ## The sentinel and the default both resolve to gamma_min = 0: fine.
    PIController().set_lte_band('auto', 'auto', 'auto')
    PIController().set_lte_band()


def test_solution_lte_controller_standard_path_contract():
    """SolutionLTEController.evaluate_step -- the selectable standard-
    path form of Fang's solution-space estimator.  Gated: the
    accept-unevaluated openings (no history, no ORDER on the
    integrator), the controlling-index record on a normal accept, the
    upper-band rejection shrinking the step, and the lower-band
    growth-retry.  History convention: x_hist[0] is the most recent
    accepted point; degree = min(ORDER, len(x_hist)-1, len(hs))."""
    x0, x1 = np.array([0.0, 0.0]), np.array([1.0, 0.0])
    hist = [x1, x0]                       # linear extrapolation -> [2, 0]

    ## Openings: too little history, and an integrator without ORDER.
    ctrl = SolutionLTEController()
    accept, h_next = _evaluate(ctrl, StubIntegrator([0.0, 0.0]),
                               J=np.eye(2), x_hist=None)
    assert accept and h_next == 1e-6
    class NoOrder:                        # duck-typed, no ORDER attribute
        def compute_lte(self, **kw): raise AssertionError('not called')
    accept, h_next = _evaluate(ctrl, NoOrder(), J=np.eye(2), x_hist=hist)
    assert accept and h_next == 1e-6
    ## ...and an order-0 integrator leaves no degree to extrapolate with.
    accept, h_next = _evaluate(ctrl, StubIntegrator([0.0, 0.0], order=0),
                               J=np.eye(2), x_hist=hist)
    assert accept and h_next == 1e-6

    ## Normal accept: computed point ON the extrapolation, controlling
    ## index recorded, error observable and ~0.
    ctrl = SolutionLTEController()
    accept, _ = _evaluate(ctrl, StubIntegrator([0.0, 0.0], order=1),
                          J=np.eye(2), x_curr=[2.0, 0.0], x_hist=hist)
    assert accept
    assert ctrl.controlling_index == 0
    assert ctrl.last_err == pytest.approx(0.0, abs=1e-9)

    ## Upper band: a gross deviation rejects and shrinks.
    ctrl = SolutionLTEController()
    accept, h_next = _evaluate(ctrl, StubIntegrator([0.0, 0.0], order=1),
                               J=np.eye(2), x_curr=[10.0, 0.0], x_hist=hist)
    assert not accept
    assert h_next < 1e-6
    assert h_next >= 1e-6 * MIN_SHRINK_RATIO * (1 - 1e-12)
    assert ctrl.last_err > 1.0

    ## Lower band: far-under-tolerance deviation is redone larger.
    ctrl = SolutionLTEController().set_lte_band(gamma_min=0.5,
                                                gamma_max=3.0)
    accept, h_next = _evaluate(ctrl, StubIntegrator([0.0, 0.0], order=1),
                               J=np.eye(2), x_curr=[2.0 + 1e-9, 0.0],
                               x_hist=hist)
    assert not accept
    assert h_next > 1e-6
