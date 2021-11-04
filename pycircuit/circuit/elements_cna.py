from __future__ import division
# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from circuit import *

class L(SubCircuit):
    """Inductor

    Implemented with a gyrator connected to a grounded capacitor
   """
    terminals = ('plus', 'minus', 'capnode', 'refnode')
    instparams = [Parameter(name='L', desc='Inductance', 
                            unit='H', default=1e-9)]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        ## Gyrator
        self['gyrator'] = Gyrator('plus', 'minus', 'capnode', 'refnode', gm=1)
        ## Capacitor with a 
        self['cap']     = Cap('capnode', 'refnode', c = 'L')

class VS(SubCircuit):
    """Independent DC voltage source

    Implemented with a current source connected to an gyrator and a refnode
    """
    terminals = ('plus', 'minus','','refnode')
    instparams = [Parameter(name='v', desc='Source DC voltage', 
                            unit='V', default=0),
                  Parameter(name='vac', desc='AC analysis amplitude', 
                            unit='V', default=1),
                  Parameter(name='phase', desc='AC analysis phase', 
                            unit='deg', default=0),
                  Parameter(name='noisePSD', 
                            desc='Voltage noise power spectral density', 
                            unit='V^2/Hz', default=0)]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        ## Gyrator
        self['gyrator'] = Gyrator('plus', 'minus', 'capnode', 'refnode', gm=1)
        ## Capacitor with a 
        self['IS']     = IS('capnode', 'refnode', c = 'L')


class Nullor(Circuit):
    """Nullor

       For Nullor documentations se first the Nullor in elements.py

       This Nullor is for Compact Nodal Analysis (CNA) instead of
       Modified Nodal Analysis (MNA). The Nullor equations are 
       solved through matrix transformations before solving the 
       system of equations see [1]
       
       1. Esteban Tlelo Cuautle / Arturo Sarmiento Reyes
          A PURE NODAL-ANALYSIS METHOD SUITABLE FOR ANALOG CIRCUITS USING
          NULLORS
          Journal of Applied Research and Technology, october, año/vol. 1,
          número 003
       
    """

    terminals = ('inp', 'inn', 'outp', 'outn')

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        n = self.n
        G = self.toolkit.zeros((n,n))
        self._G = G

    def G(self, x, epar=defaultepar): return self._G

class VCCS(SubCircuit):
    """Voltage controlled current source


    """
    instparams = [Parameter(name='gm', desc='transconductance',unit='A/V', 
                            default=1)]
    terminals = ('inp', 'inn', 'outp', 'outn')
                               
    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        i_p, i_n = self.add_nodes('internal_p', 'internal_n')

        self['Nullor_p'] = Nullor('inp',i_p,'outp',i_p)
        self['Nullor_n'] = Nullor('inn',i_n,'outn',i_n)
        self['res_gain'] = R(i_p,i_n, R = gm)

class CCVS(SubCircuit):
    """Current controlled voltage source


    """
    instparams = [Parameter(name='g', desc='transresistance',unit='V/A', 
                            default=1)]
    terminals = ('inp', 'inn', 'outp', 'outn')
                               
    def __init__(self, *args, **kvargs):
        # Question: is this intenional or is it intended to call super of CCVS?
        # In the latter case use 'super().__init__' instead.
        # VCCS and CCVS both inherit from SubCircuit, so it shouldn't make a difference.
        super(VCCS, self).__init__(*args, **kvargs)

        self['Nullor'] = Nullor('inp', 'inn', 'outp', 'outn')
        self['res_gain'] = R('inp', 'outp', R = g)

class VCVS(SubCircuit):
    """Voltage controlled voltage source


    """
    instparams = [Parameter(name='g', desc='gain',unit='V/V', 
                            default=1)]
    terminals = ('inp', 'inn', 'outp', 'outn')
                               
    def __init__(self, *args, **kvargs):
        # Question: is this intenional or is it intended to call super of VCVS?
        # In the latter case use 'super().__init__' instead.
        # VCCS and VCVS both inherit from SubCircuit, so it shouldn't make a difference. 
        super(VCCS, self).__init__(*args, **kvargs)
        i_p, i_n = self.add_nodes('internal_p', 'internal_n')

        self['VCCS'] = VCCS('inp', 'inn', i_p, i_n, g = g)
        self['CCVS'] = CCVS(i_p, i_n, 'outp', 'outn', g = 1)

class CCCS(SubCircuit):
    """Voltage controlled voltage source


    """
    instparams = [Parameter(name='g', desc='gain',unit='V/V', 
                            default=1)]
    terminals = ('inp', 'inn', 'outp', 'outn')
                               
    def __init__(self, *args, **kvargs):
        # Question: is this intenional or is it intended to call super of CCCS?
        # In the latter case use 'super().__init__' instead.
        # VCCS and CCCS both inherit from SubCircuit, so it shouldn't make a difference. 
        super(VCCS, self).__init__(*args, **kvargs)
        i_p, i_n = self.add_nodes('internal_p', 'internal_n')

        self['CCVS'] = CCVS('inp', 'inn', i_p, i_n, g = g)
        self['VCCS'] = VCCS(i_p, i_n, 'outp', 'outn', g = 1)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
