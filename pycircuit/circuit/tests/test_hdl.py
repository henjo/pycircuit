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
    
    print "apa", res.i(npy.array([1, 0]))

