"""MOS transistor models
"""
from numpy import array
from circuit import SubCircuit, VCCS, R, C, Parameter, gnd

class MOSSmallSignal(SubCircuit):
    """Small-signal MOS model

    >>> from sympy import Symbol
    >>> from symbolic import SymbolicTwoPort
    >>> c = SubCircuit()
    >>> inp = c.addNode('inp')
    >>> out = c.addNode('outp')
    >>> c['q1'] = MOSSmallSignal(inp, out, gnd, gnd, gm = Symbol('gm'), gds = Symbol('gds'), \
                                 Cgs = Symbol('Cgs'), Cgd = 0*Symbol('Cgd'))
    >>> res = SymbolicTwoPort(c, inp, gnd, out, gnd).run(freqs = Symbol('s'), complexfreq=True)
    >>> print res['ABCD']
    
    """
    terminals = ('g', 'd', 's', 'b')
    instparams = [Parameter(name='gm', desc='Gate transconductance', unit='A/V', default=1e-4),
                  Parameter(name='gmb', desc='Bulk transconductance', unit='A/V', default=1e-5),
                  Parameter(name='gds', desc='Gate transconductance', unit='A/V', default=1e-6),
                  Parameter(name='Cgs', desc='Bulk transconductance', unit='A/V', default=0.),
                  Parameter(name='Cgd', desc='Bulk transconductance', unit='A/V', default=0.),
                  Parameter(name='Cdb', desc='Bulk transconductance', unit='A/V', default=0.)
                  ]
    def __init__(self, *args, **kvargs):
        super(MOSSmallSignal, self).__init__(*args, **kvargs)

        self['Igm'] = VCCS(self.nodenames['g'], self.nodenames['s'], self.nodenames['d'], self.nodenames['s'], gm=self.ipar.gm)

        self['Igmb'] = VCCS(self.nodenames['b'], self.nodenames['s'], self.nodenames['d'], self.nodenames['s'], gm=self.ipar.gmb)

        self['gds'] = R(self.nodenames['d'], self.nodenames['s'], r = 1/self.ipar.gds)

        self['Cgs'] = C(self.nodenames['g'], self.nodenames['s'], c = self.ipar.Cgs)
        self['Cgd'] = C(self.nodenames['g'], self.nodenames['d'], c = self.ipar.Cgd)
        self['Cdb'] = C(self.nodenames['d'], self.nodenames['b'], c = self.ipar.Cdb)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
