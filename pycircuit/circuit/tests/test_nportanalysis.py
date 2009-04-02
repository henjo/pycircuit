# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""
Test n-port analysis module

"""

from pycircuit.circuit import *
from pycircuit.circuit.nport import NPort, NPortY, NPortZ, NPortA, NPortS
from pycircuit.circuit.nportanalysis import TwoPortAnalysis

from math import sqrt
import numpy as np
from sympy import Symbol
from numpy.testing import assert_array_almost_equal

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
    ## Create circuit
    cir = SubCircuit()
    ## n1,n2 = nodes('1','2')
    cir['R1'] = R(1, 2, r=Symbol('R1'))
    cir['C1'] = C(2, gnd, c=Symbol('C1'))

    ## Run symbolic 2-port analysis
    twoport_ana = SymbolicTwoPortAnalysis(cir, Node('1'), gnd, Node('2'), gnd)
    result = twoport_ana.solve(freqs=Symbol('s'), complexfreq=True)

    ## Print ABCD parameter matrix
    ABCD = Matrix(result['twoport'].A)
    ABCD.simplify()
    print ABCD
