# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

from pycircuit.circuit.elements import VSin, ISin, IS, R, L, C, SubCircuit, gnd
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.stepcontroller import IntegralController
from pycircuit.circuit import circuit #new
from math import floor
import numpy as np
import pytest
import unittest

from pycircuit.circuit import Circuit, defaultepar
from pycircuit.utilities.param import Parameter

class myC(Circuit):
    """Capacitor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(np.zeros(2))
    array([[ 0.,  0.],
           [ 0.,  0.]])
    >>> c.C(np.zeros(2))
    array([[  1.0000e-12,  -1.0000e-12],
           [ -1.0000e-12,   1.0000e-12]])

    """

    terminals = ('plus', 'minus')
    instparams = [Parameter(name='c0', desc='Capacitance', 
                            unit='F', default=1e-12),
                  Parameter(name='c1', desc='Nonlinear capacitance', 
                            unit='F', default=0.5e-12),
                  Parameter(name='v0', desc='Voltage for nominal capacitance', 
                            unit='V', default=1),
                  Parameter(name='v1', desc='Slope voltage ...?', 
                            unit='V', default=1)
                  ]

    def C(self, x, epar=defaultepar): 
        v=x[0]-x[1]
        c0 = self.ipar.c0
        c1 = self.ipar.c1
        v0 = self.ipar.v0
        v1 = self.ipar.v1
        c = c0+c1*self.toolkit.tanh((v-v0)/v1)
        return self.toolkit.array([[c, -c],
                                  [-c, c]])

    def q(self, x, epar=defaultepar):
        v=x[0]-x[1]
        c0 = self.ipar.c0
        c1 = self.ipar.c1
        v0 = self.ipar.v0
        v1 = self.ipar.v1
        q = c0*v+c1*v1*self.toolkit.ln(self.toolkit.cosh((v-v0)/v1))
        return self.toolkit.array([q, -q])
    


def test_transient_RC():
    """Test of the of transient simulation of RC-circuit
    """
    circuit.default_toolkit = circuit.numeric
    
    c = SubCircuit()

    n1 = c.add_node('net1')
    n2 = c.add_node('net2')
    c['Is'] = IS(gnd, n1, i=10)    
    c['R1'] = R(n1, gnd, r=1)
    c['R2'] = R(n1, n2, r=1e3)
    c['R3'] = R(n2, gnd, r=100e3)
    c['C'] = C(n2, gnd, c=1e-5)
    tran = Transient(c)
    x0_zeros = np.zeros(c.n)
    res = tran.solve(tend=10e-3,timestep=1e-4, x0=x0_zeros)
    expected = 6.3
    assert  abs(res.v(n2,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result.'

    
def test_transient_RLC():
    """Test of transient simulation of RLC-circuit
    """
    
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    
    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = C(2, gnd, c=1e-12)
    #c['L'] = L(2,gnd, L=1e-3)
    tran_imp = Transient(c)
    res_imp = tran_imp.solve(tend=40e-6,timestep=1e-6, fixed_timestep=True)
    expected = 2.58
    assert  abs(res_imp.v(2,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result.'

def test_transient_nonlinear_C():
    """Test of transient simulation of RLC-circuit,
    with nonlinear capacitor.
    """
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    
    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = myC(2, gnd)
    #c['L'] = L(2,gnd, L=1e-3)
    tran_imp = Transient(c)
    res_imp = tran_imp.solve(tend=40e-6,timestep=1e-6, fixed_timestep=True)
    expected = 3.4
    assert  abs(res_imp.v(2,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result:'

def test_transient_get_diff():
    """Test of differentiation method
    """
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = C(2, gnd, c=1e-12)
    tran = Transient(c)
    tran._dt=1e-6
    x0=np.ones(c.n)
    q=c.q(x0)
    Cmatrix=c.C(x0)
    print(tran.parameters)
    tran.base_integrator = tran._get_integrator()
    hist_len = max(2, tran.base_integrator.get_required_history())
    tran._qlast=np.zeros((hist_len,tran.cir.n))#initialize q-history vector
    tran._iqlast=np.zeros((hist_len,tran.cir.n))
    tran._is_first_step = False
    iq,geq = tran.get_diff(q,Cmatrix)
    print(iq,geq)


def test_transient_methods_step_response():
    """Targeted step-response integration test to explicitly verify solvers.
    """
    circuit.default_toolkit = circuit.numeric
    from pycircuit.circuit.elements import VS
    from pycircuit.circuit.integrator import EulerIntegrator, TrapezoidalIntegrator, Gear2Integrator

    integrators = {
        'euler': EulerIntegrator,
        'trapezoidal': TrapezoidalIntegrator,
        'trap': TrapezoidalIntegrator,
        'gear2': Gear2Integrator,
    }
    expected_results = {
        'euler': [0.5, 0.75, 0.875],
        'trapezoidal': [0.5, 5/6, 17/18],
        'trap': [0.5, 5/6, 17/18],
        'gear2': [0.5, 0.8, 0.94]
    }

    for method, expected in expected_results.items():
        c = SubCircuit()
        n1 = c.add_node('n1')
        n2 = c.add_node('n2')
        c['vin'] = VS(n1, gnd, v=1.0)
        c['R'] = R(n1, n2, r=1)
        c['C'] = C(n2, gnd, c=1.0)

        tran = Transient(c, integrator=integrators[method]())
        x0_zeros = np.zeros(c.n)
        result = tran.solve(tend=3.0, timestep=1.0, x0=x0_zeros, fixed_timestep=True)

        computed = result.v(n2, gnd).y

        np.testing.assert_allclose(computed, expected, rtol=1e-5, err_msg=f"Failed for method {method}")
def test_transient_adaptive_efficiency():
    """Test comparing adaptive time step efficiency versus fixed time step.
    Verifies adaptive takes fewer steps and both pass tolerance checks for all methods.
    """
    circuit.default_toolkit = circuit.numeric
    from pycircuit.circuit.elements import IS
    from pycircuit.circuit.integrator import EulerIntegrator, TrapezoidalIntegrator, Gear2Integrator

    integrators = {
        'euler': EulerIntegrator,
        'trap': TrapezoidalIntegrator,
        'trapezoidal': TrapezoidalIntegrator,
        'gear2': Gear2Integrator,
    }
    methods = ['euler', 'trap', 'trapezoidal', 'gear2']

    for method in methods:
        c = SubCircuit()
        n1 = c.add_node('net1')
        n2 = c.add_node('net2')
        c['Is'] = IS(gnd, n1, i=10)
        c['R1'] = R(n1, gnd, r=1)
        c['R2'] = R(n1, n2, r=1e3)
        c['R3'] = R(n2, gnd, r=100e3)
        c['C'] = C(n2, gnd, c=1e-5)

        tran = Transient(c, integrator=integrators[method]())
        x0_zeros = np.zeros(c.n)
        
        res_fixed = tran.solve(tend=10e-3, timestep=1e-4, x0=x0_zeros, fixed_timestep=True)
        res_adapt = tran.solve(tend=10e-3, timestep=1e-3, x0=x0_zeros, fixed_timestep=False)
        
        fixed_steps = len(res_fixed.sweep_values)
        adapt_steps = len(res_adapt.sweep_values)
        
        expected = 6.3
        
        assert abs(res_fixed.v(n2, gnd)[-1] - expected) < 2e-2*expected, f'Fixed step failed tolerance check for {method}.'
        assert abs(res_adapt.v(n2, gnd)[-1] - expected) < 2e-2*expected, f'Adaptive step failed tolerance check for {method}.'
        assert adapt_steps < fixed_steps, f'Adaptive step took {adapt_steps} steps, which is not less than fixed {fixed_steps} steps for {method}.'


if __name__ == '__main__':
    #test_transient_RC()
    test_transient_RLC()
    test_transient_nonlinear_C()
    #test_transient_get_diff()

def test_transient_pcnr_diode():
    """Test PCNR limiting in Transient simulation with a diode driven by a large voltage step.
    
    Reasoning:
    Without PCNR, driving a non-linear junction with such an extreme voltage step 
    (10V, 1ns rise time) usually causes the Newton-Raphson iteration's initial 
    forward-voltage prediction to wildly overshoot (e.g. hundreds of volts). 
    This would cause a Python OverflowError when evaluating I = Is * exp(V/Vt), 
    or completely diverge.
    
    With our PCNR pnjlim algorithm enabled in Diode, the iteration limits the 
    exponential argument gracefully. The solver successfully clamps and solves 
    the diode at a realistic forward bias voltage.
    """
    from pycircuit.circuit.elements import Diode, VPulse, R
    c = SubCircuit()
    
    # A 10V step source with very fast rise time (which typically causes NR overshoot)
    c['VPulse'] = VPulse(1, gnd, v1=0, v2=10, tr=1e-9, tf=1e-9, pw=1e-3, per=2e-3)
    c['R1'] = R(1, 2, r=10) # Small resistor to allow large current
    c['D1'] = Diode(2, gnd)
    
    tran = Transient(c)
    # The solver should converge without OverflowError due to PCNR pnjlim
    res = tran.solve(tend=2e-3, timestep=1e-5)
    
    # Check that diode forward voltage settles around 0.7-0.9V (typical for 1A current)
    # The current will be roughly (10V - 0.8V) / 10 = 0.92A
    v_diode_max = max(res.v(2, gnd))
    assert 0.5 < v_diode_max < 1.5, f"Diode voltage {v_diode_max} outside expected forward bias range"


def test_transient_coupled_lte():
    """Test Option A: Coupled Nonlinear System for Adaptive Timestepping."""
    from pycircuit.circuit.elements import VS, R, C
    c = SubCircuit()
    c['VS'] = VS(1, gnd, v=10)
    c['R1'] = R(1, 2, r=10)
    c['C1'] = C(2, gnd, c=1e-6)
    
    tran = Transient(c)
    # Using coupled_lte=True to trigger Schur Complement Option A solver
    res = tran.solve(tend=50e-6, timestep=1e-6, coupled_lte=True)
    
    # Simple check that simulation completed and final time is close to tend
    assert len(res.sweep_values) > 5, "Coupled simulation took too few steps"
    assert abs(res.sweep_values[-1] - 50e-6) < 1e-6, "Simulation did not reach tend"
    
    # Simple RC charging check at final time (tau = 10 us, tend = 50 us -> ~5 tau)
    # V_c(5tau) ≈ 10 * (1 - e^-5) ≈ 9.93V
    v_c_final = res.v(2, gnd)[-1]
    assert 9.0 < v_c_final < 10.1, f"RC voltage {v_c_final} is way off"


def test_transient_adaptive_vs_coupled():
    """Compare adaptive time stepping (Option B) vs coupled solver (Option A)."""
    from pycircuit.circuit.elements import VPulse, R, C, Diode
    c = SubCircuit()
    
    # Non-linear circuit to challenge the step size controllers
    c['VPulse'] = VPulse(1, gnd, v1=0, v2=5, tr=1e-6, tf=1e-6, pw=10e-6, per=20e-6)
    c['R1'] = R(1, 2, r=100)
    c['D1'] = Diode(2, 3)
    c['C1'] = C(3, gnd, c=1e-9)
    c['R2'] = R(3, gnd, r=1e3)
    
    tran = Transient(c)
    
    # 1. Option B (Standard Adaptive LTE)
    res_adapt = tran.solve(tend=40e-6, timestep=1e-7, coupled_lte=False)
    steps_adapt = len(res_adapt.sweep_values)
    
    # 2. Option A (Coupled Schur Complement Solver)
    res_coupled = tran.solve(tend=40e-6, timestep=1e-7, coupled_lte=True)
    steps_coupled = len(res_coupled.sweep_values)
    
    # Compare final state accuracy
    v_adapt_final = res_adapt.v(3, gnd)[-1]
    v_coupled_final = res_coupled.v(3, gnd)[-1]
    
    err = abs(v_adapt_final - v_coupled_final)
    assert err < 1e-2, f"Option A and Option B diverged! Error: {err}"
    
    # Compare step counts
    # The coupled solver should ideally require fewer or comparable steps 
    # since it never rejects steps, though it may take smaller steps on average.
    ratio = steps_coupled / max(1, steps_adapt)
    
    import warnings
    if ratio < 0.2 or ratio > 5.0:
        warnings.warn(f"Time steps differ widely! Adaptive: {steps_adapt}, Coupled: {steps_coupled}")
    
    # Just to ensure it's not going haywire, bounded check
    assert steps_coupled < steps_adapt * 10, "Coupled solver took way too many steps!"

def test_transient_pi_controller():
    """Test the newly extracted PIController for step size control."""
    from pycircuit.circuit.elements import VPulse, R, C
    from pycircuit.circuit.stepcontroller import PIController
    
    c = SubCircuit()
    c['V'] = VPulse(1, gnd, v1=0, v2=5, tr=1e-6, tf=1e-6, pw=10e-6, per=20e-6)
    c['R'] = R(1, 2, r=100)
    c['C'] = C(2, gnd, c=1e-9)
    
    tran = Transient(c)
    # Inject PIController
    tran.step_controller = PIController()
    res = tran.solve(tend=40e-6, timestep=1e-7, coupled_lte=False)
    
    assert len(res.sweep_values) > 10, "PIController did not take enough steps"
    assert abs(res.sweep_values[-1] - 40e-6) < 1e-6, "PIController did not reach tend"


## ---------------------------------------------------------------------------
## LTE-estimate helpers, shared by the tests below.
##
## `compute_lte` is required to estimate the error in the *companion current*,
## `iq - q'(t_n)` -- the residual the step controller then maps to solution
## units with `J^-1`.  These helpers evaluate that at unit level from an
## analytic q(t), with no circuit, which is the cheap check that was missing
## when Gear2's 'classic' branch estimated the wrong derivative for years.
## ---------------------------------------------------------------------------

_LTE_W = 2 * np.pi * 1e6


def _lte_q(t):
    """Analytic charge vector, two independent smooth components."""
    return np.array([np.sin(_LTE_W * t), 0.5 * np.cos(0.7 * _LTE_W * t)])


def _lte_dq(t):
    """Exact dq/dt of _lte_q."""
    return np.array([_LTE_W * np.cos(_LTE_W * t),
                     -0.35 * _LTE_W * np.sin(0.7 * _LTE_W * t)])


def _lte_vs_onestep_true(integ, h, t_n=3.0e-7):
    """(estimate, true) for one uniform step, against the one-step LTE.

    "One-step" is the textbook local truncation error: exact charge history
    *and* exact past derivatives, so the value measured is the error this single
    step commits and nothing accumulated.  Each formula then has an exactly
    derivable ratio to it -- 1/2 for Backward Euler, 5/6 for Trapezoidal, 2/3
    for Gear2 'classic', 1/2 for Gear2 'ywr' -- which makes these numbers pins
    rather than fitted constants.

    Deliberately *not* the alternative reference, in which the companion
    history is built by running the integrator's own recursion forward: that one
    is ill-posed for Trapezoidal, whose companion current depends on its own
    past value through a recursion with eigenvalue -1, so its error carries an
    undamped alternating component whose size depends on how the history was
    seeded.  See doc/transient_repair_plan.md, gate 1-1.
    """
    from pycircuit.circuit.toolkit import numeric
    nhist = max(2, integ.get_required_history())
    ts = [t_n - k * h for k in range(nhist + 2)]
    qs = [_lte_q(t) for t in ts]
    iq_hist = [_lte_dq(ts[1 + j]) for j in range(nhist)]
    q_last = [qs[1 + j] for j in range(nhist)]
    ident = np.eye(len(qs[0]))
    iq, _geq = integ.compute_derivatives(qs[0], ident, h, q_last, iq_hist, h,
                                         False, numeric)
    est, _p = integ.compute_lte(qs[0], h, q_last, iq_hist, h, False, numeric)
    return np.asarray(est, dtype=float), np.asarray(iq, dtype=float) - _lte_dq(t_n)


def _lte_ratio(integ, h=2.5e-10):
    """estimate/true at the largest-magnitude component."""
    est, true = _lte_vs_onestep_true(integ, h)
    i = int(np.argmax(np.abs(true)))
    return est[i] / true[i]


class _CountingController(IntegralController):
    """IntegralController that records rejections and unchecked accepts."""

    def __init__(self):
        self.rejections = 0
        self.unchecked = 0

    def evaluate_step(self, *args, **kwargs):
        ## The flag is named `no_history` from stage 3 of the transient repair
        ## on; before that it was `is_first_step`, re-armed at every breakpoint.
        if kwargs.get('no_history', kwargs.get('is_first_step', False)):
            self.unchecked += 1
        accept, h_next = super().evaluate_step(*args, **kwargs)
        if not accept:
            self.rejections += 1
        return accept, h_next


def test_lte_formula_ywr():
    """The two LTE formulas as step-control options, and what separates them.

    This test used to assert, for Gear2, that ``'ywr'`` takes *more* steps than
    ``'classic'`` and reaches a smaller error.  Both held -- but only because
    ``'classic'`` estimated q'' scaled by h**3 where BDF-2 needs q''' scaled by
    h**2, an estimate ~1e-15 of the true truncation error.  The controller then
    never rejected a step and took the fewest steps physically possible, so the
    old assertions were satisfied *by* the defect and had to fail once it was
    repaired.  Rewritten in stage 2 of doc/transient_repair_plan.md.

    What is asserted instead is true of a correct implementation and false of
    the old one: both formulas actually control the step (they reject steps and
    they respond to ``reltol``), both track the analytic solution, and the one
    thing that still separates them is a constant -- the YWR GEAR2 residual
    estimates (1/4) h^2 q''' against a true (1/3) h^2 q''', so it reports 3/4 of
    the truncation error where the corrected classic form is asymptotically
    exact.  That 4/3 ratio between the two formulas is the assertion that keeps
    this file distinguishing them.
    """
    from pycircuit.circuit.elements import VS
    from pycircuit.circuit.integrator import EulerIntegrator, Gear2Integrator

    def run(integrator_cls, lte, reltol=1e-4):
        c = SubCircuit()
        c['VS'] = VS(1, gnd, v=10)
        c['R1'] = R(1, 2, r=10)
        c['C1'] = C(2, gnd, c=1e-6)
        tran = Transient(c, integrator=integrator_cls(lte), uic=True,
                         reltol=reltol)
        tran.step_controller = _CountingController()
        res = tran.solve(tend=50e-6, timestep=5e-6, coupled_lte=False)
        t = np.asarray(res.sweep_values, dtype=float)
        v_analytic = 10 * (1 - np.exp(-t[-1] / 10e-6))
        return (len(t), abs(res.v(2, gnd)[-1] - v_analytic),
                tran.step_controller.rejections)

    # Backward Euler: the two formulas are mathematically identical.
    n_ec, e_ec, _r_ec = run(EulerIntegrator, 'classic')
    n_ey, e_ey, _r_ey = run(EulerIntegrator, 'ywr')
    assert n_ec == n_ey
    assert abs(e_ec - e_ey) < 1e-9

    # Gear2 under *both* formulas: the controller is alive and controlling.
    for lte in ('classic', 'ywr'):
        n4, e4, rej4 = run(Gear2Integrator, lte, reltol=1e-4)
        n6, _e6, _rej6 = run(Gear2Integrator, lte, reltol=1e-6)
        assert rej4 >= 1, \
            "gear2-%s rejected no step at all on a 10 us RC charge" % lte
        assert n6 > 1.2 * n4, \
            "gear2-%s barely responds to reltol: %d -> %d steps" % (lte, n4, n6)
        assert e4 < 2e-2, \
            "gear2-%s off the analytic RC solution by %.3g V of 10 V" % (lte, e4)

    # What still separates the formulas: a constant 4/3, from 2/3 against 1/2
    # relative to the one-step LTE.  Not a step count -- that was the mistake
    # this test used to make.
    r_gc = _lte_ratio(Gear2Integrator('classic'))
    r_gy = _lte_ratio(Gear2Integrator('ywr'))
    assert abs(r_gc - 2.0 / 3.0) < 0.02 * (2.0 / 3.0), \
        "gear2-classic estimates %.4g of the one-step LTE, expected 2/3" % r_gc
    assert abs(r_gy - 0.5) < 0.02 * 0.5, \
        "gear2-ywr estimates %.4g of the one-step LTE, expected 1/2" % r_gy
    assert abs(r_gc / r_gy - 4.0 / 3.0) < 0.03, \
        "classic/ywr = %.4g, expected 4/3" % (r_gc / r_gy)


## ---------------------------------------------------------------------------
## The three checks the suite lacked, for every integrator.
##
## Gear2Integrator.compute_lte's 'classic' branch estimated q'' scaled by h**3
## where BDF-2 needs q''' scaled by h**2 -- about 1e-15 of the true truncation
## error -- so the step controller never rejected a step and results were
## bit-identical across a 1e3 change in reltol.  Twelve tests touched Gear2 and
## none could see it: none called compute_lte at all, none asserted that a step
## is ever rejected, and none asserted that anything responds to a tolerance.
## test_transient_adaptive_efficiency's `adapt_steps < fixed_steps` is satisfied
## most emphatically of all by a *blind* controller, which takes the fewest
## steps physically possible, so it could not serve.
##
## Each of the three runs for every integrator and every LTE formula, because
## the defect lived in one of three near-identical branches.
## See doc/transient_repair_reasoning.md (E) and plan stage 4.
## ---------------------------------------------------------------------------

## (name, class name, formula, expected order in h, expected estimate/true).
## The ratios are against the one-step LTE and every one of them is derived on
## paper, not read off the code: expanding each companion formula with exact
## past derivatives gives 1/2 for Backward Euler, 5/6 for Trapezoidal (either
## formula), 2/3 for Gear2 'classic' and 1/2 for Gear2 'ywr'.  Four independent
## constants hitting at once is what makes this a check rather than a snapshot.
_LTE_CASES = [
    ('euler-classic', 'EulerIntegrator', 'classic', 1.0, 0.5),
    ('euler-ywr', 'EulerIntegrator', 'ywr', 1.0, 0.5),
    ('trap-classic', 'TrapezoidalIntegrator', 'classic', 2.0, 5.0 / 6.0),
    ('trap-ywr', 'TrapezoidalIntegrator', 'ywr', 2.0, 5.0 / 6.0),
    ('gear2-classic', 'Gear2Integrator', 'classic', 2.0, 2.0 / 3.0),
    ('gear2-ywr', 'Gear2Integrator', 'ywr', 2.0, 0.5),
]


@pytest.mark.parametrize('name,cls_name,formula,order,ratio', _LTE_CASES)
def test_compute_lte_order_and_scale(name, cls_name, formula, order, ratio):
    """compute_lte must estimate the truncation error, to the right power of h.

    This is the cheap one -- an analytic q(t), no circuit -- and it is the one
    that would have caught the original defect on its own: the pre-repair Gear2
    'classic' branch scaled as h**3.001 instead of h**2, and its ratio to the
    true error was ~1e-15 instead of 2/3.
    """
    import pycircuit.circuit.integrator as integrator_mod
    cls = getattr(integrator_mod, cls_name)

    hs = [4e-9 / 2 ** i for i in range(4)]
    mags, ratios = [], []
    for h in hs:
        est, true = _lte_vs_onestep_true(cls(formula), h)
        i = int(np.argmax(np.abs(true)))
        mags.append(abs(est[i]))
        ratios.append(est[i] / true[i])

    observed = np.log2(mags[-2] / mags[-1])
    assert abs(observed - order) < 0.1, \
        '%s: LTE scales as h**%.3f, expected h**%.1f' % (name, observed, order)
    assert abs(ratios[-1] - ratio) < 0.02 * ratio, \
        '%s: estimate/true = %.6g at h=%.3g, expected %.4g' \
        % (name, ratios[-1], hs[-1], ratio)


## Series RLC loop released from an initial condition:
##     C(1,gnd) -- R(1,2) -- L(2,gnd)
## with i_L and v1 as states,
##     di/dt = (v1 - R i)/L,   dv1/dt = -i/C
##     A = [[-R/L, 1/L], [-1/C, 0]],   y(t) = expm(A t) y0   exactly.
## R=1k, L=1uH, C=1uF puts the poles at -1.000e+03 and -1.000e+09: a stiffness
## ratio of 1e6, and the reference is a matrix exponential rather than another
## integration, so it cannot flatter the method under test.
_STIFF_R, _STIFF_L, _STIFF_C = 1e3, 1e-6, 1e-6
_STIFF_A = np.array([[-_STIFF_R / _STIFF_L, 1.0 / _STIFF_L],
                     [-1.0 / _STIFF_C, 0.0]])


def _make_integrator(name):
    import pycircuit.circuit.integrator as integrator_mod
    base, formula = name.split('-')
    cls = {'euler': 'EulerIntegrator', 'trap': 'TrapezoidalIntegrator',
           'gear2': 'Gear2Integrator'}[base]
    return getattr(integrator_mod, cls)(formula)


def _stiff_run(name, reltol, tend=5e-3, timestep=2e-4):
    """Run the stiff case; return (steps, rejections, error after start-up)."""
    from scipy.linalg import expm
    c = SubCircuit()
    c['C1'] = C(1, gnd, c=_STIFF_C)
    c['R1'] = R(1, 2, r=_STIFF_R)
    c['L1'] = L(2, gnd, L=_STIFF_L)
    tran = Transient(c, integrator=_make_integrator(name), reltol=reltol)
    tran.step_controller = _CountingController()
    x0 = np.zeros(c.n)
    x0[c.get_node_index('1')] = 1.0
    x0[c.get_node_index('2')] = 1.0
    res = tran.solve(tend=tend, timestep=timestep, x0=x0, coupled_lte=False)

    t = np.asarray(res.sweep_values, dtype=float)
    v1 = np.asarray(res.v(1, gnd), dtype=float)
    v2 = np.asarray(res.v(2, gnd), dtype=float)

    ## The genuine first step of a run has no history, so its LTE cannot be
    ## estimated and it is accepted unevaluated at max_step.  The O(h^2) Euler
    ## error it commits is identical for every method and no later step undoes
    ## it, so a plain max-error-over-the-run measures the start-up and nothing
    ## else.  Measure instead the error accumulated *after* start-up: propagate
    ## the exact solution from the simulated state at the third point.  i_L is
    ## recovered from the resistor as (v1 - v2)/R.
    k = 2
    y_k = np.array([(v1[k] - v2[k]) / _STIFF_R, v1[k]])
    err = max(abs(v1[j] - (expm(_STIFF_A * float(t[j] - t[k])) @ y_k)[1])
              for j in range(k, len(t)))
    return len(t), tran.step_controller.rejections, float(err)


_INTEGRATORS = [c[0] for c in _LTE_CASES]


@pytest.mark.parametrize('name', _INTEGRATORS)
def test_step_is_actually_rejected_on_stiff_case(name):
    """A step must actually be rejected somewhere on a stiff run.

    This is the assertion a blind controller cannot satisfy by being blind.
    """
    steps, rejections, _err = _stiff_run(name, 1e-4)
    assert rejections >= 1, \
        '%s: %d steps and not one rejected -- the LTE estimate is not ' \
        'controlling the step size' % (name, steps)


@pytest.mark.parametrize('name', _INTEGRATORS)
def test_step_count_and_error_respond_to_reltol(name):
    """Tightening reltol must cost steps and buy accuracy, monotonically.

    Monotonicity alone would not be enough: the pre-repair Gear2 'classic'
    branch produced bit-identical results across a 1e3 change in reltol, and a
    constant is monotone.  So the response is also required to be real.
    """
    reltols = [1e-3, 1e-4, 1e-5, 1e-6]
    steps, errs = [], []
    for reltol in reltols:
        n, _rej, err = _stiff_run(name, reltol)
        steps.append(n)
        errs.append(err)

    for i in range(len(reltols) - 1):
        assert steps[i + 1] >= steps[i], \
            '%s: step count fell when reltol tightened: %s' % (name, steps)
        ## 5% slack for platform float variation only; measured behaviour is a
        ## strict decrease at every step of the sweep.
        assert errs[i + 1] <= errs[i] * 1.05, \
            '%s: error grew when reltol tightened: %s' % (name, errs)
    assert steps[-1] > 1.2 * steps[0], \
        '%s: step count barely moves across a 1e3 tolerance change (%s) -- ' \
        'the controller is not reacting' % (name, steps)
    assert errs[-1] < 0.2 * errs[0], \
        '%s: less than 5x accuracy from a 1e3 tolerance change: %s' % (name, errs)
