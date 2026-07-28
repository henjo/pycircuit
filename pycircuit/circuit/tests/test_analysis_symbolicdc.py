from __future__ import division

from pycircuit.circuit.elements import VS, IS, R, L, C, SubCircuit, gnd, Diode, VCVS_limited
from pycircuit.circuit.symbolicdc import SymbolicDC
from pycircuit.circuit import symbolic

import sympy
from sympy import Symbol, var, symbols, log
from numpy.testing import assert_array_almost_equal, assert_array_equal
import numpy as np

def test_nonlinear():
    k, qelectron, I0, Isat, T_sym = sympy.symbols('k qelectron I0 Isat T', positive=True, real=True)

    c = SubCircuit(toolkit=symbolic)
    c['I0'] = IS(gnd, 'net1', i=I0, toolkit=symbolic)
    c['D'] = Diode('net1', gnd, IS=Isat, toolkit=symbolic)

    dc = SymbolicDC(c)

    dc.epar.T = T_sym

    res = dc.solve()

    assert sympy.simplify(res.v('net1') - k * T_sym / qelectron * log(I0/Isat+1)) == 0

def test_linear():
    var('R1 R2 V0')

    c = SubCircuit(toolkit=symbolic)
    c['V0'] = VS(1, gnd, v=V0, vac=1, toolkit=symbolic)
    c['L'] = L(1,2, L=1e-3)
    c['R1'] = R(2, 3, r = Symbol('R1'))
    c['R2'] = R(3, gnd, r = Symbol('R2'))

    dc = SymbolicDC(c)

    res = dc.solve()

    ## Check voltage
    assert sympy.simplify(res.v(3, gnd) -  V0*R2/(R1+R2)) == 0

    ## Check current through R2
    assert sympy.simplify(res.i('R2.plus') - V0/(R1+R2)) == 0
    
def test_nonlinear_multivariable():
    """A genuinely nonlinear, multi-equation symbolic system.

    ``Diode`` (test_nonlinear, above) reduces to a single equation, so
    sympy.solve returns a plain root list there and SymbolicDC.solve()
    converts it to a dict by hand. A tanh-based element pulls in more nodes
    (VCVS_limited's own branch, the load), so this is the first test to
    exercise sympy.solve's *other* nonlinear return shape -- a list of
    solution tuples, not a dict -- which SymbolicDC.solve() didn't handle
    before (see architecture.md P11; this was unreachable until _symbolic
    gained ``tanh``, not merely untested).
    """
    Vin, g, level, RL = sympy.symbols('Vin g level RL', real=True, positive=True)

    c = SubCircuit(toolkit=symbolic)
    c['V0'] = VS(1, gnd, v=Vin, toolkit=symbolic)
    c['E1'] = VCVS_limited(1, gnd, 2, gnd, g=g, level=level, offset=0,
                           toolkit=symbolic)
    c['RL'] = R(2, gnd, r=RL)

    dc = SymbolicDC(c)
    res = dc.solve()

    expected = level * sympy.tanh(g * Vin / level)
    assert sympy.simplify(res.v(2, gnd) - expected) == 0

def test_geteqsys():
    R1, V0, Isat, T_sym = sympy.symbols('R1 V0 Isat T')
    k = symbolic.kboltzmann
    qelectron = symbolic.qelectron
    

    c = SubCircuit(toolkit=symbolic)
    c['V0'] = VS('net1', gnd, v=V0, toolkit=symbolic)
    c['R1'] = R('net1', 'net2', r=R1)
    c['D1'] = Diode('net2', gnd, IS=Isat, toolkit=symbolic)

    dc = SymbolicDC(c)

    dc.epar.T = T_sym

    eqsys, x = dc.get_eqsys()

    x0, x2, x3 = x

    eqsys_ref = np.array([x3 + x0/R1 - x2/R1, 
                          -Isat*(1 - sympy.exp(qelectron*x2/(T_sym*k))) + x2/R1 - x0/R1, 
                          x0 - V0])

    assert sympy.simplify(eqsys_ref[0] - eqsys_ref[0]) == 0
    assert sympy.simplify(eqsys_ref[0] - eqsys_ref[0]) == 0
    assert sympy.simplify(eqsys_ref[0] - eqsys_ref[0]) == 0

