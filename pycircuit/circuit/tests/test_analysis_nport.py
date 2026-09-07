# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""
Test n-port analysis module

"""

from pycircuit.circuit import *
from pycircuit.circuit.nport import NPort, NPortY, NPortZ, NPortA, NPortS
from pycircuit.circuit.nportanalysis import TwoPortAnalysis
from pycircuit.circuit import symbolic

from math import sqrt
import numpy as np
from sympy import Matrix, var, simplify
from numpy.testing import assert_array_almost_equal, assert_array_equal

## Import test vehicle from test_nport
from .test_nport import cir, Aref, CAref, nin, nout, NPortS, CSref, T

import unittest

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

    R1, R2, w, k, T_sym = sympy.symbols('R1 R2 w k T', real=True, positive=True)

    cir['Rp'] = R(1, gnd, r=R1/2, toolkit=symbolic)
    cir['Rn'] = R(2, gnd, r=R1/2, toolkit=symbolic)
    
    twoport_ana = TwoPortAnalysis(cir, 1,2, 2, 1,
                                  noise = True, toolkit=symbolic,
                                  noise_outquantity = 'v')
    result = twoport_ana.solve(freqs=1j*w, complexfreq=True)

    assert result['Sin'] == 4*k*T_sym/R1
    assert result['Svn'] == 0

def test_symbolic_twoport():
    circuit.default_toolkit = symbolic
    cir = SubCircuit()

    k = symbolic.kboltzmann
    R1, R0, C1, w, T_sym2 = sympy.symbols('R1 R0 C1 w T', real=True, positive=True)
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

    assert_array_equal(ABCD, np.array([[1 + 0*R1*C1*s, R1],
                                    [(1 + 0*R0*C1*s + 0*R1*C1*s) / R0,  (R0 + R1)/R0]]))

    assert_array_equal(simplify(result['Sin'] - (4*k*T_sym2/R0 + 4*R1*k*T_sym2/R0**2)), 0)
    assert_array_equal(simplify(result['Svn']), 4*k*T_sym2*R1)


def test_the_swept_noise_correlation_is_per_frequency_not_the_last_one():
    """`solve_s` computed `CS` ONCE for a whole frequency sweep.

    `S` comes from `AC` and sweeps correctly; the noise-wave correlation `CS`
    came from `TransimpedanceAnalysis`, which builds a single
    `Yreciprocal = G.T + s*C.T` that a vector `s` cannot enter. What that
    produced depended on the sweep length against the MNA size `m`:

        len(freqs) == 1   -> correct
        len(freqs) == m   -> NO ERROR; `CS` was the LAST frequency's matrix,
                             attached to a fully swept `S`
        otherwise         -> ValueError deep in the stamp

    ⚠ THE S-PARAMETERS WERE RIGHT THROUGHOUT, which is why this survived: a
    two-port sweep is usually read for `S`, and every check that reads `S`
    passes while the NOISE correlation is wrong. Measured on this fixture at
    1e9/1e11 Hz, the swept `CS` equalled the 1e11 matrix to every digit.

    So this test reads `CS`, at three frequencies (`m = 2` here, so the
    2-frequency case is the silent one and the 3-frequency case used to raise).
    """
    import pycircuit.circuit.circuit as _cc
    _cc.default_toolkit = numeric

    cir = SubCircuit()
    n1, n2 = cir.add_nodes('1', '2')
    cir['R1'] = R(n1, n2, r=50.0)
    cir['C1'] = C(n2, gnd, c=3.18e-12)

    an = TwoPortAnalysis(cir, n1, gnd, n2, gnd)
    freqs = [1e7, 1e9, 1e11]

    ## The per-frequency answer, which the scalar path always got right.
    want = [np.asarray(an.solve_s(freqs=f).CS, dtype=complex) for f in freqs]

    ## ⚠ Three frequencies: two would be the SILENT case and one would be
    ## trivially right, so neither alone would have caught this.
    swept = an.solve_s(freqs=np.array(freqs)).CS
    for i, f in enumerate(freqs):
        got = np.array([[swept[a, b][i] for b in range(2)] for a in range(2)])
        assert np.max(np.abs(got - want[i])) < 1e-30, \
            'CS at %g Hz is %r, per-frequency gives %r' % (f, got, want[i])

    ## And the last frequency's matrix must NOT be what the first slot holds --
    ## the exact shape of the defect.
    first = np.array([[swept[a, b][0] for b in range(2)] for a in range(2)])
    assert np.max(np.abs(first - want[-1])) > 1e-24, \
        'CS[0] equals the LAST frequency; the sweep collapsed again'

    ## The scalar call keeps its plain (nport, nport) matrix.
    assert np.shape(np.asarray(an.solve_s(freqs=1e9).CS)) == (2, 2)
