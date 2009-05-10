# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

from pycircuit.circuit.elements import VSin, IS, R, L, C, SubCircuit, gnd
from pycircuit.circuit.analysis import Tran_spec, Transient
from math import floor
from myCap import myC
import pylab
from pycircuit.post import plotall

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
    assert  abs(res.v(n2,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result.'

    
def test_transient_RLC():
    """Test of transient simulation of RLC-circuit
    """
    
    c = SubCircuit()

    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = C(2, gnd, c=1e-12)
    #c['L'] = L(2,gnd, L=1e-3)
    tran_imp = Transient(c)
    res_imp = tran_imp.solve(tend=150e-6,timestep=1e-6)
    expected = 0.099
    #plotall(res_imp.v(1),res_imp.v(2))
    #pylab.show()
    assert  abs(res_imp.v(1,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result.'

def test_transient_nonlinear_C():
    """Test of transient simulation of RLC-circuit,
    with nonlinear capacitor.
    """
    c = SubCircuit()

    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = myC(2, gnd)
    #c['L'] = L(2,gnd, L=1e-3)
    tran_imp = Transient(c)
    res_imp = tran_imp.solve(tend=150e-6,timestep=1e-6)
    expected = 0.099
    plotall(res_imp.v(1),res_imp.v(2))
    pylab.show()
    assert  abs(res_imp.v(1,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result.'

if __name__ == '__main__':
    test_transient_RC()
    #test_transient_RLC()
    test_transient_nonlinear_C()
    
