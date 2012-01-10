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
import unittest

from pylab import plot, show
from pycircuit.circuit.transient import Transient

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

    GVCCS = toolkit.zeros((4,4))
    GVCCS[2:4,0:2] =  np.array([[1, -1],[-1, 1]])
    yield(VCCS(1, gnd, 2, gnd, gm = gm), gm * GVCCS,
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

    GGyrator = np.array([[  0,  0,  -gm,  gm],
                      [  0,  0,   gm, -gm],
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
    G = np.array([[0, 0, 1],
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

def test_SVCVS_laplace_integrator():
    """Test SVCCS with a integrator transfer function

    """
    pycircuit.circuit.circuit.default_toolkit = symbolic

    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    a0,b0,Gdc = [sympy.Symbol(symname, real=True) for symname in
                 'a0,b0,Gdc'.split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = SVCVS( n1, gnd, n2, gnd,
                         denominator = (a0, 0),
                         numerator = (b0,))

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.v(n2,gnd)),sympy.simplify(b0/(a0*s)))

def test_SVCVS_laplace_n1_d2():
    """Test VCCS with a laplace defined transfer function first order numerator
    and second order denominator"""

    pycircuit.circuit.circuit.default_toolkit = symbolic
    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    b0,a0,a1,a2,Gdc = [sympy.Symbol(symname, real=True) for symname in
                       'b0,a0,a1,a2,Gdc'.split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = SVCVS( n1, gnd, n2, gnd,
                         denominator = (a0, a1, a2),
                         numerator   = (b0, 0))

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.v(n2,gnd)),b0*s/(a0*s*s+a1*s+a2))

@slow
def test_SVCVS_laplace_d3_n1():
    """Test VCCS with a laplace defined transfer function with second order
    numerator and third order denominator
    """

    pycircuit.circuit.circuit.default_toolkit = symbolic
    cir = SubCircuit()

    n1,n2 = cir.add_nodes('1','2')

    b0,a0,a1,a2,a3,Gdc = [sympy.Symbol(symname, real=True) for
                                symname in 'b0,a0,a1,a2,a3,Gdc'
                                .split(',')]

    s = sympy.Symbol('s', complex=True)

    cir['VS']   = VS( n1, gnd, vac=1)
    cir['VCVS'] = SVCVS( n1, gnd, n2, gnd,
                        denominator = [a0, a1, a2, a3],
                        numerator   = [b0, 0, 0])

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.v(n2,gnd)),
                 sympy.simplify((b0*s*s)/(a0*s*s*s+a1*s*s+a2*s+a3)))

def test_Idt_sym():
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

@unittest.skip("Skip failing test")
def test_Idt_tran():
    """Test integrator element in transient"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    c = SubCircuit()
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['Idt'] = Idt(nin, gnd, nout, gnd)
    
    tran = Transient(c, toolkit=numeric)
    result = tran.solve(tend=0.5,timestep=1e-2)
    y = result.v(nout).y
    x = result.v(nout).x[0]
    # vout = vin * t with constant input => test that v(nout) = t
    assert_array_equal(y[1:]/x[1:], np.ones(y[1:].size)) #avoid divide by t=0.0

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


@unittest.skip("Skip failing test")
def test_Idtmod_tran():
    """Test modulo integrator element in transient"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    c = SubCircuit()
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus = 1., offset = -0.)
    
    tran = Transient(c, toolkit=numeric)
    result = tran.solve(tend=0.5,timestep=1e-2)
    y = result.v(nout).y
    x = result.v(nout).x[0]
    # vout = vin * t with constant input => test that v(nout) = t
    assert_array_equal(y[1:]/x[1:], np.ones(y[1:].size)) #avoid divide by t=0.0

@unittest.skip("Skip failing test")
def test_Idtmod_modulo():
    """Test modulo integrator element in transient"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    c = SubCircuit()
    nin = c.add_node('in')
    nout = c.add_node('out')
     
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus = 1., offset = -0.)
    
    tran = Transient(c, toolkit=numeric)
    result = tran.solve(tend=2.0,timestep=1e-2)
    y = result.v(nout).y
    x = result.v(nout).x[0]
    # vout = vin * t with constant input => test that v(nout) = t
    assert_array_equal(y[1:]/(x[1:]%1.0), np.ones(y[1:].size)) #avoid divide by t=0.0

if __name__ == '__main__':
    test_nullor_vva()
