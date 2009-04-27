# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Circuit element tests
"""

from nose.tools import *
from pycircuit.circuit import *
from pycircuit.circuit.elements import *
import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing.decorators import slow
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

def test_VCVS_laplace_d1():
    """Test VCCS with a laplace defined transfer function, with on denominator coefficient"""
    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    a0,a1,Gdc = [sympy.Symbol(symname, real=True) for symname in 'a0,a1,Gdc'.split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = VCVS( n1, gnd, n2, gnd, g = Gdc, denominator = [a0, a1, 0])   

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.v(n2,gnd)),-Gdc/(-a0*s*s-a1*s))    

def test_VCVS_laplace_n2_d3():
    """Test VCCS with a laplace defined transfer function first order denominator and 
    second order numerator"""

    cir = SubCircuit()
                 

    n1,n2 = cir.add_nodes('1','2')

    b0,b1,a0,a1,a2,Gdc = [sympy.Symbol(symname, real=True) for symname in 'b0,b1,a0,a1,a2,Gdc'.split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = VCVS( n1, gnd, n2, gnd, g = Gdc, denominator = [a0, a1, a2], numerator = [b0, b1])   

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.v(n2,gnd)),(Gdc*b1+Gdc*b0*s)/(a0*s*s+a1*s+a2))

@slow
def test_VCVS_laplace_d4_n1_c():
    """Test VCCS with a laplace defined transfer function with first order numerator and fourth order denominator
    """

    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    b0,b1,a0,a1,a2,a3,a4,Gdc = [sympy.Symbol(symname, real=True) for 
                                symname in 'b0,b1,a0,a1,a2,a3,a4,Gdc'
                                .split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = VCVS( n1, gnd, n2, gnd, 
                        g = Gdc, denominator = [a0, a1, a2, a3, a4], 
                        numerator = [b0, b1])

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.v(n2,gnd)),(Gdc*b0*s+Gdc*b1)/(a0*s*s*s*s+a1*s*s*s+a2*s*s+a3*s+a4))

@slow
def test_VCVS_laplace_d5_n2_observable():
    """Test VCCS with a laplace defined transfer function with first order order numerator and fith order denominator"""
    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    b0,b1,b2,a0,a1,a2,a3,a4,a5,Gdc = [sympy.Symbol(symname, real=True) for 
                                      symname in 'b0,b1,b2,a0,a1,a2,a3,a4,a5,Gdc'.split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = VCVS( n1, gnd, n2, gnd, g = Gdc, denominator = [a0, a1, a2, a3, a4, a5], 
                        numerator = [b0, b1, b2], realisation = 'controlable')

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    print(res.v(n2,gnd))
    assert_equal(sympy.simplify(res.v(n2,gnd)),sympy.simplify((Gdc*b2+Gdc*b1*s+Gdc*b0*s*s)/(a0*s*s*s*s*s+a1*s*s*s*s+a2*s*s*s+a3*s*s+a4*s+a5)))   

if __name__ == '__main__':
    test_nullor_vva()
