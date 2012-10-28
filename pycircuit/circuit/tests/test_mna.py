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

from sympy import var, Symbol, simplify
import sympy
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from numpy.testing.decorators import slow

def test_mna():
    var('L1 C2 RL vs v1 v2 i0 i1')

    cir = SubCircuit()
    cir['V1'] = VS(1, gnd, vac=1)
    cir['L1'] = L(1, 2, L=L1)
    cir['C2'] = C(2, gnd, c=C2)
    cir['RL'] = R(2, gnd, r=RL)

    mna = MNA(cir, toolkit=symbolic)

    x = np.array([v1,v2,i0,i1])
    mna.update(x, defaultepar)

    iref = [i0 + i1, -i1 + v2/RL, v1, v1 - v2]
    qref = [0, C2*v2, 0, -L1*i1]
    uref = [0, 0, -1, 0]

    assert_array_equal(map(simplify, mna.i - iref), 0)
    assert_array_equal(map(simplify, mna.q - qref), 0)
    assert_array_equal(map(simplify, mna.u - uref), 0)

    Gref = np.array([[0, 0, 1, 1],
                     [0, 1/RL, 0, -1],
                     [1, 0, 0, 0],
                     [1, -1, 0, 0]], dtype=object)
    assert_array_equal(mna.G, Gref)

    Cref = np.array([[0, 0, 0, 0],
                     [0, C2, 0, 0],
                     [0, 0, 0, 0],
                     [0, 0, 0, -L1]], dtype=object)
    assert_array_equal(mna.C, Cref)

