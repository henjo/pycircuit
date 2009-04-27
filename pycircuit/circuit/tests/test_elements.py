# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

from nose.tools import *
from pycircuit.circuit import AC, VS, VSin, R, G, C, L, Nullor, SubCircuit, \
    gnd, symbolic
import numpy as np
from numpy.testing import assert_array_equal
from sympy import var, Symbol, simplify
import sympy

def test_nullor_vva():
    """Test nullor element by building a V-V amplifier"""
    
    c = SubCircuit()

    Vin = Symbol('Vin')
    R1 =Symbol('R1')
    R2 = Symbol('R2')
    
    nin = c.add_node('in')
    n1 = c.add_node('n1')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, vac=Vin)
    c['R1'] = R(n1, gnd, r=R1)
    c['R2'] = R(nout, n1, r=R2)
    c['nullor'] = Nullor(n1, nin, gnd, nout)
    
    result = AC(c, toolkit=symbolic).solve(Symbol('s'))
    
    vout = result.v(nout)

    assert simplify(vout - Vin * (R1 + R2) / R1) == 0, \
        'Did not get the expected result, %s != 0'% \
        str(simplify(vout - Vin * (R1 + R2) / R1))

def test_vsin():
    var('vo va freq td theta phase t')
    vsin = VSin(toolkit = symbolic,
                vo=vo, va=va, freq=freq, td=td, theta=theta, phase=phase)

    v = vo + va*sympy.exp(-theta*(t - td)) * \
        sympy.sin(2*sympy.pi*freq*(t-td)+phase*sympy.pi/180)
    assert_array_equal(vsin.u(t), np.array([0,0,-v]))
           

def gen_stamps():
    var('R1 C1 L1')

    yield(R(1,gnd, r=R1), 1/R1 * np.array([[1, -1], [-1, 1]]), np.zeros((2,2)))
    yield(G(1,gnd, g=1/R1), 1/R1 * np.array([[1, -1], [-1, 1]]), 
          np.zeros((2,2)))
    yield(C(1,gnd, c=C1), np.zeros((2,2)), C1 * np.array([[1, -1], [-1, 1]]))

    GL = np.array([[0,0,1], [0,0,-1], [1, -1, 0]])
    CL = np.zeros((3,3), dtype=object)
    CL[2,2] = -L1
    yield(L(1,gnd, L=L1), GL, CL)

def test_stamp():

    for cir, G, C in gen_stamps():
        assert_array_equal(cir.G(np.zeros(cir.n)), G)
        assert_array_equal(cir.C(np.zeros(cir.n)), C)

if __name__ == '__main__':
    test_nullor_vva()
