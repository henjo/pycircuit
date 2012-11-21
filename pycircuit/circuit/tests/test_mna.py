# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

""" Test mna module
"""

from nose.tools import *
from ..circuit import *
from ..elements import *
from ..na import MNA
from .. import symbolic

from sympy import var, Symbol, simplify, symbols
import sympy
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from numpy.testing.decorators import slow

def gen_circuit(tk):
    if tk is symbolic:
        L1,C2,RL,vs,v1,v2,i0,i1 = symbols('L1 C2 RL vs v1 v2 i0 i1')
    else:
        L1 = 1e-3
        C2 = 1e-12
        RL = 1e3
        v1 = 1e-2
        v2 = 2e-3
        i0 = 1e-5
        i1 = 2e-5

    cir = SubCircuit()
    cir['V1'] = VS(1, gnd, vac=1)
    cir['L1'] = L(1, 2, L=L1)
    cir['C2'] = C(2, gnd, c=C2)
    cir['RL'] = R(2, gnd, r=RL)

    x = np.array([v1,v2,i0,i1])

    iref = [i0 + i1, -i1 + v2/RL, -v1, v2 - v1]
    qref = [0, C2*v2, 0, L1*i1]
    uref = [0, 0, 1, 0]
    Gref = np.array([[0, 0, 1, 1],
                     [0, 1/RL, 0, -1],
                     [-1, 0, 0, 0],
                     [-1, 1, 0, 0]])
    Cref = np.array([[0, 0, 0, 0],
                     [0, C2, 0, 0],
                     [0, 0, 0, 0],
                     [0, 0, 0, L1]])

    return cir, x, iref, qref, uref, Gref, Cref

def do_mna(tk):
    cir, x, iref, qref, uref, Gref, Cref = gen_circuit(tk)
    
    mna = MNA(cir, toolkit=tk)

    ## Run update twice to assure that all vectors are cleared before update
    for count in range(2):
        mna.update(x, defaultepar)

        assert_array_equal(map(simplify, mna.i - iref), 0)
        assert_array_equal(map(simplify, mna.q - qref), 0)
        assert_array_equal(map(simplify, mna.u - uref), 0)

        assert_array_equal(tk.todense(mna.G), Gref)

        assert_array_equal(tk.todense(mna.C), Cref)

def test_mna_symbolic():
    do_mna(symbolic)

def test_mna_theano():
    from .. import symbolic, theanotk
    do_mna(theanotk)

def test_describe_state_vector():
    tk = symbolic
    cir, x, iref, qref, uref, Gref, Cref = gen_circuit(tk)
    
    mna = MNA(cir, toolkit=tk)
    mna.setup()

    assert_array_equal(mna.describe_x_vector(),
                       [Node(1).V, Node(2).V, cir['V1'].I, cir['L1'].I])

    assert_array_equal(mna.describe_i_vector(),
                       [Node(1).I, Node(2).I, cir['V1'].V, cir['L1'].V])


def test_parallel():
    pycircuit.circuit.circuit.default_toolkit = numeric

    cir=SubCircuit()

    res = 1e3
    cir['R1'] = R(1, 2, res)
    cir['R2'] = R(1, 2, res)

    G = cir.G(np.array([0,0]))

    assert_array_equal(G, np.array([[2/res, -2/res],
                                    [-2/res, 2/res]]))

def test_short_resistor():
    """Test shorting of instance terminals"""
    cir = SubCircuit()

    cir['R1'] = R(gnd, gnd)
    
    assert_equal(cir.G(np.zeros(1)), np.array([0]))
    
## Move to NA test!!
def test_VCCS_tied():
    """Test VCCS with some nodes tied together"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    cir = SubCircuit()

    n3,n2 = cir.add_nodes('3','2')

    gm1 = sympy.Symbol('gm1')

    cir['gm'] = VCCS(gnd, n3, n2, n3, gm = gm1)   
    
    assert_array_equal(cir.G(np.zeros(cir.n)),
                       np.array([[gm1, 0, -gm1],
                                 [-gm1, 0, gm1],
                                 [0, 0, 0]]))

