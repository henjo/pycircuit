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

## `MOS_ACM` WAS DELETED HERE (stage 5+.2), and the reason is worth keeping.
##
## It advertised an ACM model.  It was a verbatim copy of `MOS` with two
## differences, both copy-paste damage rather than a different model:
##
##   * `gds` described as "Gate transconductance" -- it is the output conductance;
##   * the noise PSD used `Symbol('kT')`, a free symbol nothing in this package
##     binds, where `MOS` uses `toolkit.kboltzmann * Symbol('T')`.
##
## And it could not be constructed at all: `__init__` called `super(MOS, self)`
## from a class whose MRO is [MOS_ACM, SubCircuit, Circuit, object], raising
## `TypeError`.  Zero references anywhere outside this file -- no test, no example,
## no doc page, not exported from `circuit/__init__.py`.  Deleting it turns a
## TypeError at construction into an ImportError at import, which is earlier and
## clearer.
##
## **The mechanism that let it survive matters more than the class.**  The only
## thing that would ever have caught it is the doctest in its own docstring, and
## `pytest.ini` collects no doctests -- so the test existed and never ran, which
## reads as coverage while providing none.  See `tests/test_doctests.py`, which is
## this project's answer to that, and gate 5+.2b in `doc/transient_work_plan.md`
## for why `mos.py` is not yet on its list.
##
## Reconsider if someone actually wants ACM -- in which case it is written from
## the paper, not recovered from this.

if __name__ == "__main__":
    import doctest
    doctest.testmod()
