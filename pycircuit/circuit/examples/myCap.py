from __future__ import division
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

    def update(self, subject):
        c = self.ipar.c0+self.ipar.c1
        self._C =  self.toolkit.array([[c, -c],
                                       [-c, c]])

    def C(self, x, epar=defaultepar):
        v=x[0]-x[1]
        c = self.ipar.c0+self.ipar.c1*self.toolkit.tanh((v-self.ipar.v0)/self.ipar.v1)
        return self.toolkit.array([[c, -c],
                                  [-c, c]])

    def q(self, x, epar=defaultepar):
        v=x[0]-x[1]
        q = self.ipar.c0*v+c1*v1*ln(cosh((v-self.ipar.v0)/self.ipar.v1))
        return self.toolkit.array([q, -q])
    
