# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

from pycircuit.circuit.elements import IS, R, L, C, SubCircuit, gnd
from pycircuit.circuit.analysis import Tran_spec, Transient
from math import floor

def test_transient_RC():
    """Test of the of transient simulation of RC-circuit
    """
    
    c = SubCircuit()

    n1 = c.add_node('net1')
    n2 = c.add_node('net2')
    c['Is'] = IS(gnd, n1, i=10)    
    c['R1'] = R(n1, gnd, r=1)
    c['R2'] = R(n1, n2, r=1e3)
    c['R3'] = R(n2, gnd, r=100e3)
    c['C'] = C(n2, gnd, c=1e-5)
    tran = Transient(c)
    res = tran.solve(tend=10e-3,timestep=1e-4)
    expected = 6.3
    assert  abs(res[0][1] - expected) < 1e-2*expected,\
        'Does not math QUCS result.'

    
def test_transient_RLC():
    """Test of transient simulation of RLC-circuit
    """
    
    c = SubCircuit()

    n1 = c.add_node('net1')
    c['Is'] = IS(gnd, n1, i=10e-3)    
    c['R1'] = R(n1, gnd, r=1e3)
    c['C'] = C(n1, gnd, c=1e-5)
    c['L'] = L(n1,gnd, L=1e-3)
    tran_imp = Transient(c)
    tran_exp = Tran_spec(c)
    res_imp = tran_imp.solve(tend=150e-6,timestep=1e-6)
    expected = 0.099

    assert  abs(tran_imp.result[-1][0] - expected) < 1e-2*expected,\
        'Does not math QUCS result.'

if __name__ == '__main__':
    test_transient_RC()
    test_transient_RLC()
