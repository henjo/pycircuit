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


def no_test_Idt_sym():
    """Test integrator element symbolically"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    c = SubCircuit()

    Vin = Symbol('Vin')
    R1 = Symbol('R1')
    
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, vac=Vin)
    c['R1'] = R(nout, gnd, r=R1)
    c['Idt'] = Idt(nin, gnd, nout, gnd)
    
    result = AC(c, toolkit=symbolic).solve(Symbol('s'),complexfreq=True)
    
    vtr = simplify(result.v(nout)/result.v(nin))
    assert_equal(vtr, 1/Symbol('s'))


def test_Idt_tran():
    """Test integrator element in transient"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    c = SubCircuit()
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, v=70.)
    c['R1'] = R(nout, gnd, r=1e6)
    c['Idt'] = Idt(nin, gnd, nout, gnd)
    
    tran = Transient(c)
    result = tran.solve(tend=10e-3,timestep=1e-5)
    plot(result.v(nout))
    show()

# def test_Idtmod_tran():
#     """Test modulo integrator element in transient"""
#     pycircuit.circuit.circuit.default_toolkit = numeric

#     c = SubCircuit()
#     nin = c.add_node('in')
#     nout = c.add_node('out')
     
#     #c['vin'] = VSin(nin, gnd, va=1., freq = 1e3)
#     c['vin'] = VS(nin, gnd, v=20.)
#     c['R1'] = R(nout, gnd, r=1e3)
#     c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus = 1., offset = -0.5)
    
#     tran = Transient(c)
#     result = tran.solve(tend=10e-3,timestep=1e-5)
#     #plot(result.v(nin),result.v(nout))
#     plot(result.v(nin))
#     plot(result.v(nout))
#     show()

