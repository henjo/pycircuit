# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""MOS transistor models
"""
from numpy import array
from pycircuit.circuit import SubCircuit, VCCS, G, C, IS, Parameter, gnd, \
    symbolic, TwoPortAnalysis
from sympy import Symbol

class MOS(SubCircuit):
    """Small-signal MOS model

    >>> from sympy import Symbol
    >>> c = SubCircuit()
    >>> inp = c.add_node('inp')
    >>> out = c.add_node('outp')
    >>> c['q1'] = MOS(inp, out, gnd, gnd, \
                  gm = Symbol('gm'), gds = Symbol('gds'), \
                  Cgs = Symbol('Cgs'), Cgd = 0*Symbol('Cgd'))
    >>> twoport = TwoPortAnalysis(c, inp, gnd, out, gnd, toolkit=symbolic)
    >>> res = twoport.solve(freqs = array([Symbol('s')]), complexfreq=True)
    
    """
    terminals = ('g', 'd', 's', 'b')
    instparams = [Parameter(name='gm', desc='Gate transconductance', 
                            unit='A/V', default=1e-4),
                  Parameter(name='gmb', desc='Bulk transconductance', 
                            unit='A/V', default=1e-5),
                  Parameter(name='gds', desc='Output transconductance', 
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

        self['Igm'] = VCCS('g', 's', 
                           'd', 's', 
                           gm=self.ipar.gm,
                           toolkit=self.toolkit)

        self['Igmb'] = VCCS('b', 's', 
                            'd', 's', 
                            gm=self.ipar.gmb,
                            toolkit=self.toolkit)

        self['gds'] = G('d', 's', 
                        g = self.ipar.gds, noisy = False)

        self['Cgs'] = C('g', 's', 
                        c = self.ipar.Cgs)
        self['Cgd'] = C('g', 'd', 
                        c = self.ipar.Cgd)
        self['Cdb'] = C('d', 'b', 
                        c = self.ipar.Cdb)

        toolkit = self.toolkit
        kt = toolkit.kboltzmann * Symbol('T')

        inoisepsd = 4 * kt * self.ipar.gamma * self.ipar.gm
        self['idnoise'] = IS('d', 's', 
                             i = 0, iac = 0, 
                             noisePSD = inoisepsd)

class MOS_ACM(SubCircuit):
    """Small-signal MOS model based on ACM model 

    >>> from sympy import Symbol
    >>> c = SubCircuit()
    >>> inp = c.add_node('inp')
    >>> out = c.add_node('outp')
    >>> c['q1'] = MOS_ACM(inp, out, gnd, gnd, \
                  gm = Symbol('gm'), gds = Symbol('gds'), \
                  Cgs = Symbol('Cgs'), Cgd = 0*Symbol('Cgd'))
    >>> twoport = TwoPortAnalysis(c, inp, gnd, out, gnd, toolkit=symbolic)
    >>> res = twoport.solve(freqs = array([Symbol('s')]), complexfreq=True)
    
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

        self['Igm'] = VCCS('g', 's', 
                           'd', 's', 
                           gm=self.ipar.gm,
                           toolkit=self.toolkit)

        self['Igmb'] = VCCS('b', 's', 
                            'd', 's', 
                            gm=self.ipar.gmb,
                            toolkit=self.toolkit)

        self['gds'] = G('d', 's', 
                        g = self.ipar.gds, noisy = False)

        self['Cgs'] = C('g', 's', 
                        c = self.ipar.Cgs)
        self['Cgd'] = C('g', 'd', 
                        c = self.ipar.Cgd)
        self['Cdb'] = C('d', 'b', 
                        c = self.ipar.Cdb)

        inoisepsd = 4 * Symbol('kT') * self.ipar.gamma * self.ipar.gm
        self['idnoise'] = IS('d', 's', 
                             i = 0, iac = 0, 
                             noisePSD = inoisepsd)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
