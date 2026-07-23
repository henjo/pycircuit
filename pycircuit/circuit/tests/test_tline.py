import pytest
import numpy as np
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import VS, R, VSin, VSwitch, TLine
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient

def test_tline_dc():
    # Test that T-Line behaves correctly at DC
    # V1 = V2, I1 = -I2
    c = SubCircuit()
    c['V1'] = VS('in', gnd, v=5.0)
    c['R1'] = R('in', 'p1', r=50.0)
    c['T1'] = TLine('p1', gnd, 'p2', gnd, Z0=50.0, TD=1e-9)
    c['R2'] = R('p2', gnd, r=50.0)
    
    dc = DC(c)
    res = dc.solve()
    
    vp1 = res.v('p1')
    vp2 = res.v('p2')
    
    # Voltage divider: 50 ohm source, 50 ohm load.
    # DC behavior should be a wire (short). Vp1 = Vp2 = 2.5V
    assert abs(vp1 - 2.5) < 1e-6
    assert abs(vp2 - 2.5) < 1e-6

def test_tline_transient_pulse():
    # Test short pulse reflection
    c = SubCircuit()
    c['V1'] = VS('in', gnd, v=0.0)
    c['R1'] = R('in', 'p1', r=50.0)
    c['T1'] = TLine('p1', gnd, 'p2', gnd, Z0=50.0, TD=1.0)
    c['R2'] = R('p2', gnd, r=1e6)
    
    # We apply a short pulse by overriding V1 in the transient
    tr = Transient(c)
    
    # Custom stimulus: pulse of 1.0V between 0.1s and 0.3s
    def pulse_func(t):
        if 0.1 <= t <= 0.3:
            return 1.0
        return 0.0
        
    c['V1'].function.f = pulse_func
    
    res = tr.solve(tend=3.0, timestep=0.05, x0=np.zeros(c.n))
    
    t = res.sweep_values
    vp1 = res.v('p1')
    vp2 = res.v('p2')
    
    # Pulse hits p1 at t=0.1. vp1 should be 0.5V (voltage divider)
    idx_02 = int(np.argmin(np.abs(t - 0.2)))
    assert abs(vp1[idx_02] - 0.5) < 0.05
    assert abs(vp2[idx_02] - 0.0) < 0.05
    
    # Pulse is gone at t=0.4
    idx_04 = int(np.argmin(np.abs(t - 0.4)))
    assert abs(vp1[idx_04] - 0.0) < 0.05
    
    # Pulse arrives at p2 at t=1.1, ends at 1.3
    # At p2 it doubles due to open circuit: 0.5 * 2 = 1.0V
    idx_12 = int(np.argmin(np.abs(t - 1.2)))
    assert abs(vp2[idx_12] - 1.0) < 0.05
    
    # Reflected pulse arrives back at p1 at t=2.1, ends at 2.3
    # At p1 it is fully absorbed because R1 = 50 ohms matches T1
    # Thus vp1 = 0.5V
    idx_22 = int(np.argmin(np.abs(t - 2.2)))
    assert abs(vp1[idx_22] - 0.5) < 0.05

def test_tline_transient_reflection():
    c = SubCircuit()
    
    # We will use a DC source that is switched on, but since we want a step from 0 to 1,
    # we can use VSin or we can just start with 0 initial conditions and a DC source
    # Wait, if we use UIC, a DC source is 0 for t<0 and V for t>0!
    c['V1'] = VS('in', gnd, v=1.0)
    c['R1'] = R('in', 'p1', r=50.0) # Match
    c['T1'] = TLine('p1', gnd, 'p2', gnd, Z0=50.0, TD=1.0)
    c['R2'] = R('p2', gnd, r=1e6) # Open circuit
    
    # We expect a step of 0.5V at p1 at t=0.
    # The wave travels to p2, taking 1.0s. At t=1.0s, p2 jumps to 1.0V.
    # Reflection travels back, arriving at p1 at t=2.0s.
    # At t=2.0s, p1 jumps from 0.5V to 1.0V.
    
    tr = Transient(c)
    # Solve from 0 to 3s
    res = tr.solve(tend=3.0, timestep=0.1, x0=np.zeros(c.n)) # UIC
    
    t = res.sweep_values
    vp1 = res.v('p1')
    vp2 = res.v('p2')
    
    # At t=0.5, vp1 should be 0.5, vp2 should be 0.0
    idx_05 = int(np.argmin(np.abs(t - 0.5)))
    assert abs(vp1[idx_05] - 0.5) < 0.05
    assert abs(vp2[idx_05] - 0.0) < 0.05
    
    # At t=1.5, vp1 should be 0.5, vp2 should be 1.0
    idx_15 = int(np.argmin(np.abs(t - 1.5)))
    assert abs(vp1[idx_15] - 0.5) < 0.05
    assert abs(vp2[idx_15] - 1.0) < 0.05
    
    # At t=2.5, vp1 should be 1.0, vp2 should be 1.0
    idx_25 = int(np.argmin(np.abs(t - 2.5)))
    assert abs(vp1[idx_25] - 1.0) < 0.05
    assert abs(vp2[idx_25] - 1.0) < 0.05
