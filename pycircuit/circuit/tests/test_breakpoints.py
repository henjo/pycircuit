import numpy as np
import pytest
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import R, C, VSin, VPulse
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.toolkit import numeric

def test_transient_breakpoints():
    """Verify that Transient simulator correctly steps exactly on pulse breakpoints."""
    c = SubCircuit()
    c['VP'] = VPulse('in', gnd, v1=0, v2=1, td=0.1e-3, tr=1e-6, tf=1e-6, pw=0.1e-3, per=0.5e-3)
    c['R'] = R('in', 'out', r=1.0)
    c['C'] = C('out', gnd, c=1e-6)
    
    tran = Transient(c, method='gear2')
    res = tran.solve(tend=0.6e-3, timestep=10e-6)
    
    # Check if breakpoints were hit exactly using np.isclose
    t_points = res.sweep_values
    
    # Pulse edges should be strictly in the time points
    assert np.any(np.isclose(t_points, 0.1e-3))  # td
    assert np.any(np.isclose(t_points, 0.1e-3 + 1e-6))  # td + tr
    assert np.any(np.isclose(t_points, 0.2e-3 + 1e-6))  # td + tr + pw
    assert np.any(np.isclose(t_points, 0.2e-3 + 2e-6))  # td + tr + pw + tf
    
    # Second period
    assert np.any(np.isclose(t_points, 0.6e-3))  # per + td

def test_transient_minbreak():
    """Verify that Transient simulator respects minbreak parameter to skip exceedingly close events."""
    c = SubCircuit()
    # Create two pulses that have breakpoints ridiculously close to each other (e.g. 1e-15s apart)
    c['VP1'] = VPulse('in1', gnd, v1=0, v2=1, td=1e-6, tr=1e-6, tf=1e-6, pw=1e-6, per=5e-6)
    # Start VP2 just 1e-15 seconds after VP1
    c['VP2'] = VPulse('in2', gnd, v1=0, v2=1, td=1e-6 + 1e-15, tr=1e-6, tf=1e-6, pw=1e-6, per=5e-6)
    
    c['R1'] = R('in1', 'out', r=1.0)
    c['R2'] = R('in2', 'out', r=1.0)
    c['C'] = C('out', gnd, c=1e-6)
    
    # With minbreak = 1e-14, the simulator should skip the 1e-15s delta and combine them into a single step
    # or skip the second breakpoint because it's too close to the first.
    tran = Transient(c, method='euler', minbreak=1e-14)
    res = tran.solve(tend=2e-6, timestep=1e-6)
    
    t_points = res.sweep_values
    # Ensure there isn't a tiny dt = 1e-15 anywhere in the sweep (ignoring the last step which might clamp to tend)
    dts = np.diff(t_points)[:-1]
    assert np.min(dts) >= 1e-14, f"A dt smaller than minbreak was taken: {np.min(dts)}"

if __name__ == '__main__':
    test_transient_breakpoints()
    test_transient_minbreak()
