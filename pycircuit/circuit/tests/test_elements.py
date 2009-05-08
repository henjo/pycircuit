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
    assert_array_equal(vsin.u(t), symbolic.array([0,0,-v]))

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
            u = np.array(pulse.u(t + tstart)).astype(float).reshape(3,)
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

if __name__ == '__main__':
    test_nullor_vva()
