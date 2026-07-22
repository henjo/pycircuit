# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

from pycircuit.circuit.elements import VSin, ISin, IS, R, L, C, SubCircuit, gnd
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import circuit #new
from math import floor
import numpy as np
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
    a,b,b_=tran._method[tran.par.method] 
    tran._qlast=np.zeros((len(a),tran.cir.n))#initialize q-history vector
    tran._iqlast=np.zeros((len(b),tran.cir.n))
    tran._is_first_step = False
    iq,geq = tran.get_diff(q,Cmatrix)
    print(iq,geq)


def test_transient_methods_step_response():
    """Targeted step-response integration test to explicitly verify solvers.
    """
    circuit.default_toolkit = circuit.numeric
    from pycircuit.circuit.elements import VS
    
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
        
        tran = Transient(c)
        tran.par.method = method
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
        
        tran = Transient(c)
        tran.par.method = method
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
    
    # Only report/assert if they differ by more than 3x
    if ratio > 3.0 or ratio < 0.33:
        import warnings
        warnings.warn(f"Time steps differ widely! Adaptive: {steps_adapt}, Coupled: {steps_coupled}")
    
    # Just to ensure it's not going haywire, bounded check
    assert steps_coupled < steps_adapt * 10, "Coupled solver took way too many steps!"

