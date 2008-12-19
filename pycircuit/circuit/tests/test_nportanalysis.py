# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""
Test n-port analysis module

"""

from pycircuit.circuit.constants import *
from pycircuit.circuit.circuit import SubCircuit, R, gnd
from pycircuit.circuit.nport import NPort, NPortY, NPortZ, NPortA, NPortS
from pycircuit.circuit.nportanalysis import TwoPortAnalysis

from math import sqrt
import numpy as npy
from numpy.testing import assert_array_almost_equal

## Import test vehicle from test_nport
from test_nport import cir, Aref, CAref, nin, nout, NPortS

def test_twoportanalysis():
    result = TwoPortAnalysis(cir, nin, gnd, nout, gnd, method='aparam').solve(freqs = 0)

    assert isinstance(result['twoport'], NPortA)
    assert_array_almost_equal(result['twoport'].A.astype(float), Aref)

def test_twoportanalysis_sparam():
    result = TwoPortAnalysis(cir, nin, gnd, nout, gnd, method = 'sparam').solve(freqs = 0)

    assert isinstance(result['twoport'], NPortS)
    assert_array_almost_equal(result['twoport'].A.astype(float), Aref)

    assert_array_almost_equal(result['twoport'].CA.astype(float), CAref, decimal=25)


