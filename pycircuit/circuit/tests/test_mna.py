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
    cir = SubCircuit()
    cir['R1'] = R(1, 2, r=1e3)
    cir['R2'] = R(2, gnd, r=1e3)

    mna = MNA(cir, toolkit=symbolic)

    x = np.zeros(mna.n)

    mna.update(x, 0)

    print mna.G

