# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

from nose.tools import *
import pycircuit.circuit.circuit 
from pycircuit.circuit import *
from pycircuit.circuit.elements import *
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from numpy.testing.decorators import slow
from sympy import var, Symbol, simplify, symbols
import sympy

from pylab import plot, show
from pycircuit.circuit.transient import Transient


def test_Idtmod_sym():
    """Test modulus integrator element symbolically"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    c = SubCircuit()

    Vin = Symbol('Vin')
    R1 = Symbol('R1')
    
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, vac=Vin)
    c['R1'] = R(nout, gnd, r=R1)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd)
    
    result = AC(c, toolkit=symbolic).solve(Symbol('s'),complexfreq=True)
    
    vtr = simplify(result.v(nout)/result.v(nin))
    assert_equal(vtr, 1/Symbol('s'))



def test_Idtmod_tran():
    """Test modulo integrator element in transient"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    c = SubCircuit()
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, v=1.01-0.1)
    c['R1'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus = 1., offset = -0.5)
    
    tran = Transient(c)
    result = tran.solve(tend=1.0,timestep=1e-2)
    plot(result.v(nout))
    show()

if __name__ == "__main__":
    test_Idtmod_tran()
