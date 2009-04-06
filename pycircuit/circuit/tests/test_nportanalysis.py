# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""
Test n-port analysis module

"""

from nose.tools import *
from pycircuit.circuit import *
from pycircuit.circuit.nport import NPort, NPortY, NPortZ, NPortA, NPortS
from pycircuit.circuit.nportanalysis import TwoPortAnalysis
from pycircuit.circuit import symbolic

from math import sqrt
import numpy as np
from sympy import Matrix, var, simplify
from numpy.testing import assert_array_almost_equal, assert_array_equal

## Import test vehicle from test_nport
from test_nport import cir, Aref, CAref, nin, nout, NPortS, CSref, T

def test_twoportanalysis():
    result = TwoPortAnalysis(cir, nin, gnd, nout, gnd, method='aparam').solve(freqs = 0)

    assert isinstance(result['twoport'], NPortA)
    assert_array_almost_equal(result['twoport'].A.astype(float), Aref)

def test_twoportanalysis_sparam():

    ana = TwoPortAnalysis(cir, nin, gnd, nout, gnd, method = 'sparam')
    ana.epar.T = T

    result = ana.solve(freqs = 0)

    assert isinstance(result['twoport'], NPortS)
    assert_array_almost_equal(result['twoport'].A.astype(float), Aref)

    assert_array_almost_equal(result['twoport'].CA.astype(complex),
                              CAref, decimal=25)

def test_symbolic_twoport():
    cir = SubCircuit()

    var('R1 C1 w kT', real=True, positive=True)
    s = 1j*w

    cir['R1'] = R(1, 2, r=R1)
    cir['C1'] = C(2, gnd, c=C1)

    ## Run symbolic 2-port analysis
    twoport_ana = TwoPortAnalysis(cir, Node('1'), gnd, Node('2'), gnd, 
                                  noise = True, toolkit=symbolic)
    result = twoport_ana.solve(freqs=s, complexfreq=True)

    ABCD = Matrix(result['twoport'].A)
    ABCD.simplify()

    assert_array_equal(ABCD, array([[1 + R1*C1*s, R1],
                                    [C1*s,  1]]))

    ## Run symbolic 2-port analysis 
    twoport_ana = TwoPortAnalysis(cir, Node('1'), gnd, Node('2'), gnd, 
                                  noise = True, toolkit=symbolic,
                                  noise_outquantity='i')
    result = twoport_ana.solve(freqs=s, complexfreq=True)

    ABCD = Matrix(result['twoport'].A)
    ABCD.simplify()

    assert_array_equal(ABCD, array([[1 + R1*C1*s, R1],
                                    [C1*s,  1]]))

