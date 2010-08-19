# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

from pycircuit.circuit.elements import VSin, IS, R, L, C, SubCircuit, gnd
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import circuit, SubCircuit #new
from math import floor
import pylab
from pycircuit.post import plotall

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
        c = self.ipar.c0+self.ipar.c1*self.toolkit.tanh((v-self.ipar.v0)/self.ipar.v1)
        return self.toolkit.array([[c, -c],
                                  [-c, c]])

    def q(self, x, epar=defaultepar):
        v=x[0]-x[1]
        q = self.ipar.c0*v+c1*v1*ln(cosh((v-self.ipar.v0)/self.ipar.v1))
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
    res = tran.solve(tend=10e-3,timestep=1e-4)
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
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()

    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = myC(2, gnd)
    #c['L'] = L(2,gnd, L=1e-3)
    tran_imp = Transient(c)
    res_imp = tran_imp.solve(tend=150e-6,timestep=1e-6)
    expected = 0.099
    #plotall(res_imp.v(1),res_imp.v(2))
    #pylab.show()
    assert  abs(res_imp.v(1,gnd)[-1] - expected) < 1e-2*expected,\
        'Does not match QUCS result.'

if __name__ == '__main__':
    test_transient_RC()
    #test_transient_RLC()
    test_transient_nonlinear_C()
    
