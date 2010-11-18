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

def test_noise2():
    cir = SubCircuit(toolkit=symbolic)

    var('R1 R2 w k T', real=True, positive=True)

    cir['Rp'] = R(1, gnd, r=R1/2, toolkit=symbolic)
    cir['Rn'] = R(2, gnd, r=R1/2, toolkit=symbolic)
    
    twoport_ana = TwoPortAnalysis(cir, 1,2, 2, 1,
                                  noise = True, toolkit=symbolic,
                                  noise_outquantity = 'v')
    result = twoport_ana.solve(freqs=1j*w, complexfreq=True)

    assert_equal(result['Sin'], 4*k*T/R1)
    assert_equal(result['Svn'], 0)

def test_symbolic_twoport():
    circuit.default_toolkit = symbolic
    cir = SubCircuit()

    var('R1 R0 C1 w k T', real=True, positive=True)
    s = 1j*w

    cir['R0'] = R(1, gnd, r=R0)
    cir['R1'] = R(1, 2, r=R1)
#    cir['C1'] = C(2, gnd, c=C1)

    ## Add an AC source to verify that the source will not affect results
#    cir['IS'] = IS(1, gnd, iac=1) 

    ## Run symbolic 2-port analysis
    twoport_ana = TwoPortAnalysis(cir, Node('1'), gnd, Node('2'), gnd,
                                  noise = True, toolkit=symbolic,
                                  noise_outquantity = 'v')
    result = twoport_ana.solve(freqs=s, complexfreq=True)

    ABCD = Matrix(result['twoport'].A)
    ABCD.simplify()

    assert_array_equal(ABCD, array([[1 + 0*R1*C1*s, R1],
                                    [(1 + 0*R0*C1*s + 0*R1*C1*s) / R0,  (R0 + R1)/R0]]))

    assert_array_equal(simplify(result['Sin'] - (4*k*T/R0 + 4*R1*k*T/R0**2)), 0)
    assert_array_equal(simplify(result['Svn']), 4*k*T*R1)
