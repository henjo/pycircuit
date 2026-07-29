import pytest
from pycircuit.circuit.elements import R, C, L, VS, VPulse, VSin, Diode, Transformer, Nullor, VCVS_limited, Idtmod
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.transient import Transient
import numpy as np
import warnings

## Every test in this module drives a full transient integration; together
## they are the largest single block of suite runtime.  Deselected by default
## via pytest.ini; run them with `-m slow` or `-m ""`.  See architecture P15.
pytestmark = pytest.mark.slow

def _compare_methods(c, tend, timestep, test_name):
    tran = Transient(c)
    
    # 1. Option B (Standard Adaptive LTE)
    res_adapt = tran.solve(tend=tend, timestep=timestep, coupled_lte=False)
    steps_adapt = len(res_adapt.sweep_values)
    
    # 2. Option A (Coupled Schur Complement Solver)
    res_coupled = tran.solve(tend=tend, timestep=timestep, coupled_lte=True)
    steps_coupled = len(res_coupled.sweep_values)
    
    ratio = steps_coupled / max(1, steps_adapt)
    
    if ratio > 3.0 or ratio < 0.33:
        warnings.warn(f"[{test_name}] Time steps differ widely! Adaptive: {steps_adapt}, Coupled: {steps_coupled}")
        
    return res_adapt, res_coupled

def test_stress_stiff_rlc_pulse():
    c = SubCircuit()
    c['VP'] = VPulse(1, gnd, v1=0, v2=5, tr=1e-6, tf=1e-6, pw=10e-6, per=20e-6)
    c['R1'] = R(1, 2, r=1)
    c['L1'] = L(2, 3, L=1e-3)
    c['C1'] = C(3, gnd, c=1e-9)
    
    res_adapt, res_coupled = _compare_methods(c, 25e-6, 1e-9, "Stiff RLC")
    # Assert final voltage is near 5V if we are in the high pulse state
    # at t=25us we are in the second pulse (20-30us is high)
    assert abs(res_adapt.v(3, gnd)[-1] - 5.0) < 0.1
    assert abs(res_coupled.v(3, gnd)[-1] - 5.0) < 0.1

def test_stress_lossless_lc_tank():
    c = SubCircuit()
    c['VP'] = VPulse(3, gnd, v1=0, v2=1, tr=1e-6, tf=1e-6, pw=1e-3, per=100e-3)
    c['R1'] = R(3, 1, r=100) # Inject through resistor
    c['L1'] = L(1, gnd, L=1e-6)
    c['C1'] = C(1, gnd, c=1e-9)
    
    res_adapt, res_coupled = _compare_methods(c, 10e-6, 1e-6, "Lossless LC")
    assert abs(res_adapt.v(1, gnd)[-1] - res_coupled.v(1, gnd)[-1]) < 0.2

def test_stress_half_wave_rectifier():
    c = SubCircuit()
    c['VS'] = VSin(1, gnd, va=10, freq=50)
    c['D1'] = Diode(1, 2)
    c['C1'] = C(2, gnd, c=100e-6)
    c['R1'] = R(2, gnd, r=1e3)
    
    res_adapt, res_coupled = _compare_methods(c, 0.05, 1e-4, "Half-Wave")
    assert res_adapt.v(2, gnd)[-1] > 0.0
    assert abs(res_adapt.v(2, gnd)[-1] - res_coupled.v(2, gnd)[-1]) < 0.5

def test_stress_full_wave_bridge():
    c = SubCircuit()
    c['VS'] = VSin(1, 2, va=10, freq=50)
    c['R_src'] = R(2, gnd, r=1e-3) # ground reference for source
    c['D1'] = Diode(1, 3)
    c['D2'] = Diode(gnd, 1)
    c['D3'] = Diode(2, 3)
    c['D4'] = Diode(gnd, 2)
    c['C1'] = C(3, gnd, c=10e-6)
    c['R1'] = R(3, gnd, r=100)
    
    res_adapt, res_coupled = _compare_methods(c, 0.04, 1e-4, "Full-Wave Bridge")
    assert abs(res_adapt.v(3, gnd)[-1] - res_coupled.v(3, gnd)[-1]) < 0.5

def test_stress_charge_pump():
    c = SubCircuit()
    c['VP'] = VPulse(1, gnd, v1=-5, v2=5, tr=1e-6, tf=1e-6, pw=5e-6, per=10e-6)
    c['C1'] = C(1, 2, c=1e-6)
    c['D1'] = Diode(gnd, 2)
    c['D2'] = Diode(2, 3)
    c['C2'] = C(3, gnd, c=1e-6)
    
    res_adapt, res_coupled = _compare_methods(c, 50e-6, 1e-7, "Charge Pump")
    # Output should pump up towards ~10V
    assert res_adapt.v(3, gnd)[-1] > 7.0
    assert abs(res_adapt.v(3, gnd)[-1] - res_coupled.v(3, gnd)[-1]) < 0.5

def test_stress_transformer_inrush():
    c = SubCircuit()
    c['VP'] = VPulse(1, gnd, v1=0, v2=10, tr=1e-6, tf=1e-6, pw=1e-3, per=2e-3)
    c['R1'] = R(1, 2, r=1)
    # L1 acts as primary magnetizing inductance
    c['Lm'] = L(2, gnd, L=1e-3)
    # Ideal transformer coupling primary to secondary
    c['TX'] = Transformer(2, gnd, 3, gnd, n=1)
    c['R2'] = R(3, gnd, r=10)
    
    res_adapt, res_coupled = _compare_methods(c, 2e-3, 1e-5, "Transformer Inrush")
    assert abs(res_adapt.v(3, gnd)[-1] - res_coupled.v(3, gnd)[-1]) < 0.5

def test_stress_nullor_bandpass():
    c = SubCircuit()
    c['VP'] = VPulse(1, gnd, v1=0, v2=1, tr=1e-6, tf=1e-6, pw=10e-6, per=20e-6)
    c['R1'] = R(1, 2, r=1e3)
    c['C1'] = C(2, 3, c=1e-9)
    c['C2'] = C(2, 4, c=1e-9)
    c['R2'] = R(3, gnd, r=10e3)
    c['R3'] = R(4, 5, r=10e3)
    c['N1'] = Nullor(4, gnd, 5, gnd) # Opamp: inputs 4, gnd; outputs 5, gnd
    
    res_adapt, res_coupled = _compare_methods(c, 40e-6, 1e-7, "Nullor Filter")
    assert abs(res_adapt.v(5, gnd)[-1] - res_coupled.v(5, gnd)[-1]) < 1.0

def test_stress_deep_forward_clipper():
    c = SubCircuit()
    c['VS'] = VSin(1, gnd, va=100, freq=1e3)
    c['R1'] = R(1, 2, r=10)
    c['D1'] = Diode(2, gnd)
    c['D2'] = Diode(gnd, 2)
    
    res_adapt, res_coupled = _compare_methods(c, 2e-3, 1e-5, "Deep Clipper")
    # Voltage should be clamped to ~0.7-1.0V
    assert abs(res_adapt.v(2, gnd)[-1]) < 1.5
    assert abs(res_adapt.v(2, gnd)[-1] - res_coupled.v(2, gnd)[-1]) < 0.1

def test_stress_vcvs_limit_cycle():
    c = SubCircuit()
    c['VP'] = VPulse(1, gnd, v1=0, v2=1, tr=1e-6, tf=1e-6, pw=1e-6, per=10e-3) # kickstart
    c['R_in'] = R(1, 2, r=1e3)
    # Phase shift network
    c['C1'] = C(2, gnd, c=1e-9)
    c['R1'] = R(2, 3, r=1e3)
    c['C2'] = C(3, gnd, c=1e-9)
    c['R2'] = R(3, 4, r=1e3)
    c['C3'] = C(4, gnd, c=1e-9)
    c['R3'] = R(4, 5, r=1e3)
    # Limiting amplifier A=-29
    c['Amp'] = VCVS_limited(5, gnd, 2, gnd, g=-29, level=5)
    
    res_adapt, res_coupled = _compare_methods(c, 50e-6, 1e-7, "VCVS Limit Cycle")
    assert abs(res_adapt.v(2, gnd)[-1] - res_coupled.v(2, gnd)[-1]) < 2.0

def test_stress_delayed_stiff_avalanche():
    c = SubCircuit()
    # Stays 0 for 5ms, then pulses
    c['VP'] = VPulse(1, gnd, v1=0, v2=1, td=5e-3, tr=1e-9, tf=1e-9, pw=1e-3, per=10e-3)
    c['R1'] = R(1, 2, r=1e3)
    # Ideal integrator
    c['Idt'] = Idtmod(2, gnd, 3, gnd) 
    
    res_adapt, res_coupled = _compare_methods(c, 6e-3, 1e-5, "Delayed Avalanche")
    assert abs(res_adapt.v(3, gnd)[-1] - res_coupled.v(3, gnd)[-1]) < 0.5
