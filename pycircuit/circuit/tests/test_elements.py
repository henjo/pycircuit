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

def test_vsin():
    var('vo va freq td theta phase t')
    vsin = VSin(toolkit = symbolic,
                vo=vo, va=va, freq=freq, td=td, theta=theta, phase=phase)

    v = vo + va*sympy.exp(-theta*(t - td)) * \
        sympy.sin(2*sympy.pi*freq*(t-td)+phase*sympy.pi/180)

    assert_array_equal(vsin.u(t, analysis='tran'), symbolic.array([0,0,-v]))

def test_vpulse():
    t = sympy.Symbol('t')

    v1 = 1.1
    v2 = -0.9

    td = 0.4
    tr = 0.1
    tf = 0.1
    pw = 0.5
    per = 2.0
    
    eps = 1e-6
    
    pulse = VPulse(toolkit = symbolic,
                   v1=v1, v2=v2, td=td, tr=tr, tf=tf, pw=pw, per=per)
    
    tpoints = np.array((0,td,td+tr,td+tr+pw,td+tr+pw+tf,10))
    vpoints = np.array((v1,v1,v2,v2,v1,v1))
    
    tref = np.arange(0,per, 0.005)
    
    for tstart in 0,per:
        for t in tref:
            uref = np.array([0,0,-np.interp(t,tpoints,vpoints)])
            u = np.array(pulse.u(t + tstart, analysis='tran')).astype(float).reshape(3,)
            assert_array_almost_equal(u, uref)
           
def gen_stamps(toolkit=symbolic):
    circuit.default_toolkit = toolkit
    if toolkit is symbolic:
        R1,C1,L1,gain,gm,N = symbols('R1 C1 L1 gain gm N')
    else:
        R1=1.1e3
        C1=1e-12
        L1=1e-5
        gain=2.4
        gm=1.e-3
        N=1.2

    yield(R(1,gnd, r=R1), 1/R1 * np.array([[1, -1], [-1, 1]]), np.zeros((2,2)))

    yield(G(1,gnd, g=1/R1), 1/R1 * np.array([[1, -1], [-1, 1]]), 
          np.zeros((2,2)))

    yield(C(1,gnd, c=C1), np.zeros((2,2)), C1 * np.array([[1, -1], [-1, 1]]))

    GL = np.array([[0,0,1], [0,0,-1], [1, -1, 0]])
    CL = np.zeros((3,3), dtype=object)
    CL[2,2] = -L1
    yield(L(1,gnd, L=L1), GL, CL)

    GVCVS = np.array([[0,       0, 0,0, 0],
                      [0,       0, 0,0, 0],
                      [0,       0, 0,0, 1], 
                      [0,       0, 0,0,-1], 
                      [gain,-gain,-1,1, 0]])
    CVCVS = np.zeros((5,5))
    yield(VCVS(1, gnd, 2, gnd, g=gain),GVCVS,CVCVS)

    GVCCS = np.zeros((4,4))
    GVCCS[2:4,0:2] =  np.array([[1, -1],[-1, 1]])
    yield(VCCS(1, gnd, 2, gnd, gm = gm), gm *GVCCS,
          np.zeros((4,4)))

    GNullor = np.array([[0,0,0, 0, 0],
                        [0,0,0, 0, 0],
                        [0,0,0, 0, 1], 
                        [0,0,0, 0,-1], 
                        [1,-1,0,0, 0]])
    yield(Nullor(1, gnd, 2, gnd), GNullor, np.zeros((5,5)))

    GTransformer = np.array([[0,0,0, 0,  N],
                             [0,0,0, 0, -N],
                             [0,0,0, 0, 1], 
                             [0,0,0, 0,-1], 
                             [-1,1,N,-N, 0]])
    yield(Transformer(1, gnd, 2, gnd, n = N), GTransformer, np.zeros((5,5)))

    GGyrator = array([[  0.,  0.,  -gm,  gm],
                      [  0.,  0.,   gm, -gm],
                      [  gm, -gm,   0.,  0.],
                      [ -gm,  gm,   0.,  0.]])
    yield(Gyrator(1, gnd, 2, gnd, gm = gm), GGyrator, np.zeros((4,4)))
    

def gen_stamps_sources(toolkit=symbolic):
    circuit.default_toolkit = toolkit
    if toolkit is symbolic:
        vac, phase = symbols('vac phase')
    else:
        vac = 1.2
        phase = 30

    v = vac * toolkit.exp(1j * toolkit.pi * phase / 180.)
    G = array([[0, 0, 1],
               [0, 0,-1],
               [1, -1, 0]])
    yield(VS(1,0,vac=vac, phase=phase), G, np.zeros((3,3)), toolkit.array([0,0,-v]))

    cir = SubCircuit()
    cir['vs'] = VS(1,0,vac=vac, phase=phase)
    yield(cir, G, np.zeros((3,3)), toolkit.array([0,0,-v]))

def test_stamp():
    circuit.default_toolkit = symbolic
    
    for toolkit in numeric, symbolic:
        for cir, G, C in gen_stamps(toolkit=toolkit):
            assert_array_equal(cir.G(np.zeros(cir.n)), G)
            assert_array_equal(cir.C(np.zeros(cir.n)), C)

        for cir, G, C, u in gen_stamps_sources(toolkit=toolkit):
            assert_array_equal(cir.G(np.zeros(cir.n)), G)
            assert_array_equal(cir.C(np.zeros(cir.n)), C)
            assert_array_equal(cir.u(np.zeros(cir.n), analysis='ac'), u)

def test_nullor_vva():
    """Test nullor element by building a V-V amplifier"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

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

def test_SVCVS_laplace_d1():
    """Test VCCS with a laplace defined transfer function, with on denominator coefficient"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    a0,a1,Gdc = [sympy.Symbol(symname, real=True) for symname in 'a0,a1,Gdc'.split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = SVCVS( n1, gnd, n2, gnd, g = Gdc, denominator = [a0, a1, 0])   

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.v(n2,gnd)),sympy.simplify(Gdc/(a0*s*s+a1*s)))

def test_SVCVS_laplace_n1_d2():
    """Test VCCS with a laplace defined transfer function first order denominator and 
    second order numerator"""

    pycircuit.circuit.circuit.default_toolkit = symbolic
    cir = SubCircuit()
                 

    n1,n2 = cir.add_nodes('1','2')

    b0,a0,a1,Gdc = [sympy.Symbol(symname, real=True) for symname in 'b0,a0,a1,Gdc'.split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = SVCVS( n1, gnd, n2, gnd, g = Gdc, denominator = [a0, a1], numerator = [b0])   

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.v(n2,gnd)),(Gdc*b0)/(a0*s+a1))

@slow
def test_SVCVS_laplace_d3_n1_c():
    """Test VCCS with a laplace defined transfer function with first order numerator and third order denominator
    """

    pycircuit.circuit.circuit.default_toolkit = symbolic
    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    b0,b1,a0,a1,a2,a3,Gdc = [sympy.Symbol(symname, real=True) for 
                                symname in 'b0,b1,a0,a1,a2,a3,Gdc'
                                .split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = SVCVS( n1, gnd, n2, gnd, 
                        g = Gdc, denominator = [a0, a1, a2, a3], 
                        numerator = [b0, b1])

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.v(n2,gnd)),sympy.simplify((-1.0*Gdc*b0*s-1.0*Gdc*b1)/(-a0*s*s*s-a1*s*s-a2*s-a3)))


def test_integer_component_values():
    """Test dc analysis with integer component values
    
       As python per default uses integer arithmetics integer
       component values can lead to problems. (Python < 3.0)
    """
    c = SubCircuit()

    n1,n2 = c.add_nodes('net1', 'net2')

    c['vs'] = VS(n1, gnd, v = 9)
    c['R1'] = R( n1,  n2, r = 50)
    c['R2'] = R( n2, gnd, r = 50)

    dc = DC(c)
    res = dc.solve()
    
    assert_equal(res.v('net1'), 9.0)

if __name__ == '__main__':
    test_nullor_vva()
