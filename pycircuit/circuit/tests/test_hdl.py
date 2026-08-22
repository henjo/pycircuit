# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

""" Test high-level circuit definition
"""

from pycircuit.circuit.hdl import Behavioural, Parameter, Branch, Contribution, ddt
import sympy
import numpy as np

def test_resistor():
    """Verify simple resistor model"""
    
    class Resistor(Behavioural):
         instparams = [Parameter(name='r', desc='Resistance', unit='ohm')]
         @staticmethod
         def analog(plus, minus):
             b = Branch(plus, minus)
             return Contribution(b.I, 1/r * b.V),

    res = Resistor(r=1e3)
    
    v1,v2 = sympy.symbols(('v1', 'v2'))

    assert res.i([v1,v2]) == [1e-3*(v1-v2), -1e-3*(v1-v2)]

    assert np.all(res.G([v1,v2]) == 
                       np.array([[1e-3, -1e-3], [-1e-3, 1e-3]]))

    assert np.all(res.C([v1,v2]) == np.zeros((2,2)))

    assert np.all(res.CY([v1,v2]) == np.zeros((2,2)))

def test_capacitor():
    """Verify simple capacitance model"""
    
    class Capacitor(Behavioural):
         instparams = [Parameter(name='c', desc='Capacitance', unit='F')]
         @staticmethod
         def analog(plus, minus):
             b = Branch(plus, minus)
             return Contribution(b.I, ddt(c * b.V)),
         
    C = sympy.Symbol('C')

    cap = Capacitor(c=C)
    
    v1,v2 = sympy.symbols(('v1', 'v2'))

    assert cap.i([v1,v2]) == [0, 0]

    ## Compared MATHEMATICALLY, not structurally.  This used to assert
    ## against `sympy.expand(...)` because the term splitter expanded
    ## every right-hand side -- an implementation artifact that has since
    ## been removed, because expansion distributes products of sums and a
    ## real compact model is nothing but those: it multiplied the
    ## operation count by ~6 per nesting level and no device model could
    ## be compiled at all.  The charge is the same charge either way.
    got = cap.q([v1, v2])
    for g, want in zip(got, [C * (v1 - v2), -C * (v1 - v2)]):
        assert sympy.simplify(g - want) == 0

    assert np.all(cap.C([v1,v2]) == 
                       np.array([[C, -C], [-C, C]]))

    assert np.all(cap.G([v1,v2]) == np.zeros((2,2)))

    assert np.all(cap.CY([v1,v2]) == np.zeros((2,2)))


