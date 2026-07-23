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

if __name__ == '__main__':
    test_transient_breakpoints()
