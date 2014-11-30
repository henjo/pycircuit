# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

from pycircuit.circuit.elements import VSin, ISin, IS, R, L, C, SubCircuit, gnd, nC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import *#circuit #new
from math import floor
import numpy as np
import unittest

from pycircuit.circuit import Circuit, MNA, defaultepar, theanotk
from pycircuit.utilities.param import Parameter


#@unittest.skip("Skip failing test")
def test_transient_RC():
    """Test of the of transient simulation of RC-circuit
    """
    
    c = SubCircuit(toolkit=theanotk)

    n1 = c.add_node('net1')
    n2 = c.add_node('net2')
    c['ISin'] = ISin(gnd, n1, ia=10, freq=500,toolkit=theanotk)    
    c['R1'] = R(n1, gnd, r=1)
    c['R2'] = R(n1, n2, r=1e3)
    c['R3'] = R(n2, gnd, r=100e3)
    c['C'] = C(n2, gnd, c=1e-5)
    tran = Transient(c)
    res = tran.solve(tend=10e-3,timestep=1e-4)
#    from pylab import plot, show
#    plot(abs(res.v(2,gnd)))
#    show()
    expected = 6.3
    assert  abs(res.v(n2,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result.'

#@unittest.skip("Skip failing test")    
def test_transient_RLC():
    """Test of transient simulation of RLC-circuit
    """
    
    c = SubCircuit(toolkit=theanotk)
    
    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3, toolkit=theanotk)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = C(2, gnd, c=1e-12)
    c['L'] = L(2,gnd, L=1e-3)
    tran_imp = Transient(c, toolkit=theanotk)
    res_imp = tran_imp.solve(tend=40e-6,timestep=1e-7)
    expected = 2.58
    from pylab import plot, show
    plot(res_imp.v(2,gnd))
    show()
    assert  abs(res_imp.v(2,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result.'

#@unittest.skip("Skip failing test")
def test_transient_nonlinear_C():
    """Test of transient simulation of RLC-circuit, with nonlinear capacitor.
    """

    c = SubCircuit(toolkit=theanotk)
    
    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3, toolkit=theanotk)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = nC(2, gnd)
    c['L'] = L(2,gnd, L=1e-3)
    tran_imp = Transient(c)
    res_imp = tran_imp.solve(tend=40e-6,timestep=1e-6)
    expected = 3.4
    assert  abs(res_imp.v(2,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result:'

@unittest.skip("Skip failing test")
def test_transient_get_diff():
    """Test of differentiation method
    """
    c = SubCircuit(toolkit=theanotk)
    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3, toolkit=theanotk)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = C(2, gnd, c=1e-12)
    tran = Transient(c, toolkit=theanotk)
    tran._dt=1e-6

    na = MNA(c, toolkit=theanotk)
    x0=np.ones(c.n)
    na.update(x0, defaultepar)

    q       = na.q
    Cmatrix = na.C

    print tran.parameters
    a,b,b_=tran._method[tran.par.method] 
    tran._qlast=np.zeros((len(a),tran.cir.n))#initialize q-history vector
    iq,geq = tran.get_diff(q,Cmatrix)
    print iq,geq


if __name__ == '__main__':
    #test_transient_RC()
    #test_transient_RLC()
    test_transient_nonlinear_C()
    #test_transient_get_diff()
