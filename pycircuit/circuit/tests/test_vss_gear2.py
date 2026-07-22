import numpy as np
import pytest
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import Diode, IS, C
from pycircuit.circuit.transient import Transient

def test_vss_gear2_sudden_step_variation():
    """
    Tests that the VSS Gear-2 integration safely drops to Backward Euler 
    when a sudden, drastic time step reduction occurs, avoiding the massive 
    history polynomial distortions described by Bizzarri and Brambilla.
    """
    c = SubCircuit()
    c['is'] = IS(gnd, 1, i=1e-3) # 1mA source
    c['D1'] = Diode(1, gnd)
    c['C1'] = C(1, gnd, c=1e-9)  # 1nF capacitor

    # Initial DC state
    x0 = np.array([0.0, 0.0]) # Node 0 is gnd, Node 1 is 0V
    
    tran = Transient(c, method='gear2')
    
    # We will manually step it to simulate a sudden step reduction
    tran.irefnode = 0
    tran._dt = 1e-6
    tran._dt_last = 1e-6
    tran._is_first_step = False
    
    # Setup history
    tran._qlast = np.array([[1e-9], [0.0]]) # Past charges
    tran._iqlast = np.array([[1e-3], [1e-3]])
    
    # Normal step
    q = np.array([1.5e-9])
    cap = np.array([[1e-9]])
    
    # This should be standard gear2
    tran._dt = 1.1e-6
    iq, geq = tran.get_diff(q, cap)
    assert tran._effective_method == 'gear2', "Should remain Gear-2 for normal step variations"
    
    # Sudden drop (e.g. hitting a sharp discontinuity)
    tran._dt = 1e-8
    iq, geq = tran.get_diff(q, cap)
    assert tran._effective_method == 'euler', "Should drop to Euler when step size drops by >10x"
    
    print("VSS Gear-2 Order Drop Test Passed!")

if __name__ == '__main__':
    test_vss_gear2_sudden_step_variation()
