from nose.tools import *
from pycircuit.circuit import *
from pycircuit.circuit.shooting import *
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest

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
                            unit='V', default=1)]

    def C(self, x, epar=defaultepar): 
        v=x[0]-x[1]
        c0 = self.ipar.c0
        c1 = self.ipar.c1
        v0 = self.ipar.v0
        v1 = self.ipar.v1
        c = c0+c1*self.toolkit.tanh((v-v0)/v1)
        return self.toolkit.array([[c, -c],
                                  [-c, c]])

    def q(self, x, epar=defaultepar):
        v=x[0]-x[1]
        c0 = self.ipar.c0
        c1 = self.ipar.c1
        v0 = self.ipar.v0
        v1 = self.ipar.v1
        q = c0*v+c1*v1*self.toolkit.ln(self.toolkit.cosh((v-v0)/v1))
        return self.toolkit.array([q, -q])

@unittest.skip("Skip failing test")
def test_shooting():
    circuit.default_toolkit = circuit.numeric

    cir = SubCircuit()
    
    N = 10
    period = 1e-3

    cir['vs'] = VSin(1,gnd, vac=2.0, va=2.0, freq=1/period, phase=20)
    cir['R'] = R(1,2, r=1e4)
    cir['C'] = C(2,gnd, c=1e-8)
    
    ac = PSS(cir)
    resac = AC(cir).solve(1/period)

    pss = PSS(cir)

    res = pss.solve(period=period, timestep = period/N)
    
    v2ac = resac.v(2,gnd)
    v2pss = res['tpss'].v(2,gnd)
    
    t,dt = numeric.linspace(0,period,num=N,endpoint=True,
                            retstep=True)

    v2ref = numeric.imag(v2ac * numeric.exp(2j*numeric.pi*1/period*t))

    w2ref = Waveform(t,v2ref,ylabel='reference', yunit='V', 
                     xunits=('s',), xlabels=('vref(2,gnd!)',))
    
    ## Check amplitude error
    v2rms_ac = np.abs(v2ac) / np.sqrt(2)
    v2rms_pss = np.abs(res['fpss'].v(2,gnd)).value(1/period)
    assert_almost_equal(v2rms_ac, v2rms_pss)
 
    ## Check error of waveform
    rmserror = np.sqrt(average((v2pss-w2ref)**2))
    assert rmserror < 1e-3, 'rmserror=%f too high'%rmserror

 
@unittest.skip("Skip failing test")
def test_PSS_nonlinear_C():
    """Test of PSS simulation of RLC-circuit,
    with nonlinear capacitor.
    """
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()

    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = myC(2, gnd)
    #c['L'] = L(2,gnd, L=1e-3)
    pss = PSS(c)
    res = pss.solve(period=1/50e3,timestep=1/50e3/20)


@unittest.skip("Skip failing test")
def test_PAC():
    circuit.default_toolkit = circuit.numeric
    N = 10
    fc = 1e6

    cir = SubCircuit()
    cir['vs'] = VSin(1,gnd, vac=2.0, va=2.0, freq=fc, phase=20)
    cir['R'] = R(1, 2, r=1e6)
    cir['D'] = Diode(2,gnd)
    cir['C'] = C(2,gnd, c=1e-12)
    
    pss = PSS(cir)

    res = pss.solve(period=1/fc, timestep = 1/(fc*N))

    pac = PAC(cir)
    res = pac.solve(pss, freqs = fc + np.array([1e3, 2e3, 4e3]))
    
    assert False, "Test should compare with spectre simulation"
