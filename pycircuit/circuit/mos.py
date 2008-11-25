"""MOS transistor models
"""
from numpy import array
from circuit import SubCircuit, VCCS, G, C, IS, Parameter, gnd
from sympy import Symbol

class MOS(SubCircuit):
    """Small-signal MOS model

    >>> from sympy import Symbol
    >>> from symbolicanalysis import SymbolicTwoPortAnalysis
    >>> c = SubCircuit()
    >>> inp = c.add_node('inp')
    >>> out = c.add_node('outp')
    >>> c['q1'] = MOS(inp, out, gnd, gnd, \
                  gm = Symbol('gm'), gds = Symbol('gds'), \
                  Cgs = Symbol('Cgs'), Cgd = 0*Symbol('Cgd'))
    >>> symtwoport = SymbolicTwoPortAnalysis(c, inp, gnd, out, gnd)
    >>> res = symtwoport.run(freqs = array([Symbol('s')]), complexfreq=True)
    >>> print res['twoport']
    
    """
    terminals = ('g', 'd', 's', 'b')
    instparams = [Parameter(name='gm', desc='Gate transconductance', 
                            unit='A/V', default=1e-4),
                  Parameter(name='gmb', desc='Bulk transconductance', 
                            unit='A/V', default=1e-5),
                  Parameter(name='gds', desc='Gate transconductance', 
                            unit='A/V', default=1e-6),
                  Parameter(name='Cgs', desc='Bulk transconductance', 
                            unit='A/V', default=0.),
                  Parameter(name='Cgd', desc='Bulk transconductance', 
                            unit='A/V', default=0.),
                  Parameter(name='Cdb', desc='Bulk transconductance', 
                            unit='A/V', default=0.),
                  Parameter(name='gamma', desc='Excessive noise factor', 
                            unit='', default=1)
                  ]
    def __init__(self, *args, **kvargs):
        super(MOS, self).__init__(*args, **kvargs)

        self['Igm'] = VCCS(self.nodenames['g'], self.nodenames['s'], 
                           self.nodenames['d'], self.nodenames['s'], 
                           gm=self.ipar.gm)

        self['Igmb'] = VCCS(self.nodenames['b'], self.nodenames['s'], 
                            self.nodenames['d'], self.nodenames['s'], 
                            gm=self.ipar.gmb)

        self['gds'] = G(self.nodenames['d'], self.nodenames['s'], 
                        g = self.ipar.gds, nonoise = True)

        self['Cgs'] = C(self.nodenames['g'], self.nodenames['s'], 
                        c = self.ipar.Cgs)
        self['Cgd'] = C(self.nodenames['g'], self.nodenames['d'], 
                        c = self.ipar.Cgd)
        self['Cdb'] = C(self.nodenames['d'], self.nodenames['b'], 
                        c = self.ipar.Cdb)

        inoisepsd = 4 * Symbol('kT') * self.ipar.gamma * self.ipar.gm
        self['idnoise'] = IS(self.nodenames['d'], self.nodenames['s'], 
                             i = 0, iac = 0, 
                             noisePSD = inoisepsd)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
