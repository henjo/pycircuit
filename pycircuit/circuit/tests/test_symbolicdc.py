from nose.tools import *

from pycircuit.circuit.elements import VS, IS, R, L, C, SubCircuit, gnd, Diode
from pycircuit.circuit.symbolicdc import SymbolicDC
from pycircuit.circuit import symbolic

import sympy
from sympy import Symbol, var, symbols, log

def test_nonlinear():
    var('k qelectron I0 Isat q')

    c = SubCircuit(toolkit=symbolic)
    c['I0'] = IS(gnd, 'net1', i=I0, toolkit=symbolic)
    c['D'] = Diode('net1', gnd, IS=Isat, toolkit=symbolic)

    dc = SymbolicDC(c)

    dc.epar.T = Symbol('T')

    res = dc.solve()

    assert_equal(res.v('net1'), k * T / q * log(I0/Isat+1))

def test_linear():
    var('R1 R2 V0')

    c = SubCircuit(toolkit=symbolic)
    c['V0'] = VS(1, gnd, v=V0, vac=1, toolkit=symbolic)
    c['L'] = L(1,2, L=1e-3)
    c['R1'] = R(2, 3, r = Symbol('R1'))
    c['R2'] = R(3, gnd, r = Symbol('R2'))

    dc = SymbolicDC(c)

    res = dc.solve()

    assert_equal(sympy.simplify(res.v(3, gnd) -  V0*R2/(R1+R2)), 0)

    
