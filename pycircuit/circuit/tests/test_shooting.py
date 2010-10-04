from nose.tools import *
from pycircuit.circuit import *
from pycircuit.circuit.shooting import *
from pycircuit.post import astable, plotall, Waveform, average, db20
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from copy import copy
import pylab

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
    
#    print v2pss.astable
#    plotall(v2pss,w2ref)
#    pylab.grid()
#    pylab.show()

    ## Check amplitude error
    v2rms_ac = abs(v2ac) / np.sqrt(2)
    v2rms_pss = abs(res['fpss'].v(2,gnd)).value(1/period)
    assert_almost_equal(v2rms_ac, v2rms_pss)
 
    ## Check error of waveform
    rmserror = numeric.sqrt(average((v2pss-w2ref)**2))
    assert rmserror < 1e-3, 'rmserror=%f too high'%rmserror

 
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
    db20(res['fpss'].v(2)).stem()

#    pylab.show()
#   plotall(res['tpss'].v(1),res['tpss'].v(2))
#    pylab.show()


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

#    plotall(db20(res.v(2, gnd)), db20(res.v(1, gnd)), plotargs=('x',))
#    pylab.legend()
#    pylab.show()
    
    assert False, "Test should compare with spectre simulation"
