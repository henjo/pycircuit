# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

""" Test high-level circuit definition
"""

from pycircuit.circuit.hdl import *
import sympy
import numpy as npy

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

    assert npy.alltrue(res.G([v1,v2]) == 
                       npy.array([[1e-3, -1e-3], [-1e-3, 1e-3]]))

    assert npy.alltrue(res.C([v1,v2]) == npy.zeros((2,2)))

    assert npy.alltrue(res.CY([v1,v2]) == npy.zeros((2,2)))


