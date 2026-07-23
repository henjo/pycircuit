import pytest
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import R, C, VSin, VS
from pycircuit.circuit.transient import Transient

def test_transient_uic():
    """Test that uic=True skips DC operating point and starts from 0V."""
    c = SubCircuit()
    c['VS'] = VS('in', gnd, v=10) # 10V DC source
    c['R'] = R('in', 'out', r=1.0)
    c['C'] = C('out', gnd, c=1e-6)
    
    # Run with uic=False (default), should start at DC operating point (10V)
    tran_dc = Transient(c)
    res_dc = tran_dc.solve(tend=1e-3, timestep=1e-4)
    # The capacitor should immediately be charged to 10V because DC OP is run
    assert abs(res_dc.v('out', gnd)[0] - 10.0) < 1e-4
    
    # Run with uic=True, should skip DC operating point and start at 0V
    tran_uic = Transient(c, uic=True)
    res_uic = tran_uic.solve(tend=1e-3, timestep=1e-4)
    import numpy as np
    t_first = res_uic.sweep_values[0]
    expected_v = 10.0 * (1 - np.exp(-t_first / 1e-6))
    # It should be charging up from 0V, instead of starting fully charged
    assert abs(res_uic.v('out', gnd)[0] - expected_v) < 0.5

if __name__ == '__main__':
    test_transient_uic()
