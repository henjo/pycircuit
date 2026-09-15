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


def test_nports_z_keeps_the_imaginary_part_and_works_on_a_sweep():
    """`NPortS.Z` cast S to float; `.Z` and `.Y` both failed on a sweep.

    Reported by a peer session, reproduced and widened 2026-09-15.
    `np.asarray(self.S).astype(float)` dropped Im S with only a ComplexWarning:
    on a synthetic complex S `.Z` was off by 1.027 at z0 = 1 and 50, and on the
    R = 50 ohm / C = 3.18 pF two-port from `TwoPortAnalysis` (z0 = 1) it was
    off by 5103x at 1e7 Hz, 51x at 1e9, 1.12x at 1e11 -- `E - S` is nearly
    singular near DC, so an Im S of 3.8e-4 dominates the inverse.  And on a
    swept `solve_s(freqs=array)` result, S holds Waveforms: `.Z` raised
    ("setting an array element with a sequence") and `.Y` raised too (its
    `_inv` sends object arrays to sympy).
    """
    import pycircuit.circuit.circuit as _cc
    from pycircuit.circuit.nport import NPortS
    _cc.default_toolkit = numeric
    E = np.eye(2)

    ## 1. a complex S, two references, against the textbook conversion
    S = np.array([[0.2 + 0.3j, 0.05 - 0.01j], [0.7 + 0.2j, -0.1 + 0.4j]])
    for z0 in (1.0, 50.0):
        Z = np.asarray(NPortS(S, z0=z0).Z, dtype=complex)
        want = z0 * (E + S) @ np.linalg.inv(E - S)
        assert np.max(np.abs(Z / want - 1.0)) < 1e-12, (z0, Z, want)
        ## and back: (Z - z0)(Z + z0)^-1 recovers S
        back = (Z - z0 * E) @ np.linalg.inv(Z + z0 * E)
        assert np.max(np.abs(back - S)) < 1e-12, (z0, back, S)

    ## 2. an ordinary circuit, scalar and swept
    cir = SubCircuit()
    n1, n2 = cir.add_nodes('1', '2')
    cir['R1'] = R(n1, n2, r=50.0)
    cir['C1'] = C(n2, gnd, c=3.18e-12)
    an = TwoPortAnalysis(cir, n1, gnd, n2, gnd)
    freqs = [1e7, 1e9, 1e11]
    swept = an.solve_s(freqs=np.array(freqs))
    Zs, Ys = swept.Z, swept.Y
    for i, f in enumerate(freqs):
        one = an.solve_s(freqs=f)
        Sf = np.asarray(one.S, dtype=complex)
        want = one.z0 * (E + Sf) @ np.linalg.inv(E - Sf)
        Zf = np.asarray(one.Z, dtype=complex)
        assert np.max(np.abs(Zf / want - 1.0)) < 1e-9, (f, Zf, want)
        zi = np.array([[complex(Zs[a, b][i]) for b in range(2)]
                       for a in range(2)])
        yi = np.array([[complex(Ys[a, b][i]) for b in range(2)]
                       for a in range(2)])
        assert np.max(np.abs(zi / want - 1.0)) < 1e-9, ('swept Z', f, zi, want)
        assert np.max(np.abs(zi @ yi - E)) < 1e-9, ('swept Z*Y', f)


def test_a_converted_nport_carries_its_reference_impedance():
    """`NPortY(nports).S` renormalised to 50 ohm whatever the source's z0 was.

    Reported by a peer session (2026-09-14).  `NPortY/Z/A.S` and `.CS` were
    properties with a `z0=50` ARGUMENT, which a property can never receive, and
    the conversion constructors dropped the source's `z0`.  `TwoPortAnalysis`
    returns z0 = 1, so `NPortY(result).S` silently described a different
    reference than `result.S` did.  Now: a conversion carries the source's z0,
    `.S`/`.CS` are at the stored z0, a directly constructed port keeps 50 ohm,
    and `to_s(z0)` renormalises explicitly on every class.
    """
    import pycircuit.circuit.circuit as _cc
    from pycircuit.circuit.constants import kboltzmann
    from .test_nport import Yref, CYref, Sref
    _cc.default_toolkit = numeric
    E = np.eye(2)

    def textbook_s(Z, z0):
        return (Z - z0 * E) @ np.linalg.inv(Z + z0 * E)

    ## 1. round trip through every class, complex S, at two references
    S = np.array([[0.2 + 0.3j, 0.05 - 0.01j], [0.7 + 0.2j, -0.1 + 0.4j]])
    for z0 in (1.0, 50.0):
        src = NPortS(S, z0=z0)
        for cls in (NPortY, NPortZ, NPortA):
            port = cls(src)
            assert port.z0 == z0, (cls.__name__, z0, port.z0)
            got = np.asarray(port.S, dtype=complex)
            assert np.max(np.abs(got - S)) < 1e-12, (cls.__name__, z0, got)
            back = NPortS(port)
            assert back.z0 == z0
            assert np.max(np.abs(np.asarray(back.S, dtype=complex) - S)) \
                < 1e-12, (cls.__name__, z0)

    ## 2. a directly constructed port is unchanged: 50 ohm
    direct = NPortY(Yref, CYref)
    assert direct.z0 == 50.0
    assert_array_almost_equal(np.asarray(direct.S, dtype=float), Sref)

    ## 3. to_s(z0) on every class against the textbook, S and (Bosma) CS
    Zref = np.linalg.inv(Yref)
    T = 290.0
    ports = (NPortY(Yref, CYref), NPortZ(NPortY(Yref, CYref)),
             NPortA(NPortY(Yref, CYref)), NPortS(NPortY(Yref, CYref)))
    for port in ports:
        for z0 in (1.0, 50.0, 75.0):
            s_port = port.to_s(z0)
            assert isinstance(s_port, NPortS) and s_port.z0 == z0
            want = textbook_s(Zref, z0)
            got = np.asarray(s_port.S, dtype=complex)
            assert np.max(np.abs(got - want)) < 1e-9, \
                (type(port).__name__, z0, got, want)
            cs_want = kboltzmann * T * (E - want @ want.conj().T)
            cs_got = np.asarray(s_port.CS, dtype=complex)
            assert np.max(np.abs(cs_got - cs_want)) < 1e-9 * np.max(
                np.abs(cs_want)), (type(port).__name__, z0, cs_got, cs_want)
        ## and to_s() with no argument is the stored reference
        assert port.to_s().z0 == port.z0

    ## 4. the analysis result (z0 = 1), scalar and swept, converted and read back
    cir = SubCircuit()
    n1, n2 = cir.add_nodes('1', '2')
    cir['R1'] = R(n1, n2, r=50.0)
    cir['C1'] = C(n2, gnd, c=3.18e-12)
    an = TwoPortAnalysis(cir, n1, gnd, n2, gnd)
    one = an.solve_s(freqs=1e9)
    for cls in (NPortY, NPortZ):
        got = np.asarray(cls(one).S, dtype=complex)
        assert np.max(np.abs(got - np.asarray(one.S, dtype=complex))) < 1e-9, \
            (cls.__name__, got)
    freqs = [1e7, 1e9, 1e11]
    swept = an.solve_s(freqs=np.array(freqs))
    for cls in (NPortY, NPortZ):
        conv = cls(swept)
        assert conv.z0 == swept.z0
        Sc = conv.S
        for i in range(len(freqs)):
            got = np.array([[complex(Sc[a, b][i]) for b in range(2)]
                            for a in range(2)])
            want = np.array([[complex(swept.S[a, b][i]) for b in range(2)]
                             for a in range(2)])
            assert np.max(np.abs(got - want)) < 1e-9, (cls.__name__, i, got)

    ## 5. a cascade keeps the operands' reference
    p = NPortS(NPortY(Yref, CYref).to_s(1.0))
    cas = p * p
    assert cas.z0 == 1.0
    Zc = np.asarray(NPortA(np.asarray(NPortA(p).A, dtype=complex)
                           @ np.asarray(NPortA(p).A, dtype=complex)).Z,
                    dtype=complex)
    assert np.max(np.abs(np.asarray(cas.S, dtype=complex)
                         - textbook_s(Zc, 1.0))) < 1e-9


def test_a_swept_nport_renormalises_s_and_cs_frequency_by_frequency():
    """`to_s(z0)` on a SWEPT `solve_s` result -- `S` as Waveforms, `CS` as
    per-frequency arrays -- must equal the scalar renormalisation at every
    frequency.

    ⚠ RECORDED AS AN OPEN GAP ON 2026-09-14 ("a swept CS is not renormalised
    per frequency") WITHOUT A MEASUREMENT, AND WRONG: written 2026-09-15 to pin
    a fix, this test already PASSED on the unmodified code (to 1e-12), so the
    existing Y/Z conversion path handles the sweep and no change was made.
    Kept as the pin that the claim now has."""
    import pycircuit.circuit.circuit as _cc
    _cc.default_toolkit = numeric
    cir = SubCircuit()
    n1, n2 = cir.add_nodes('1', '2')
    cir['R1'] = R(n1, n2, r=50.0)
    cir['C1'] = C(n2, gnd, c=3.18e-12)
    an = TwoPortAnalysis(cir, n1, gnd, n2, gnd)
    freqs = [1e7, 1e9, 1e11]
    swept = an.solve_s(freqs=np.array(freqs))
    assert swept.z0 != 50.0
    ren = swept.to_s(50.0)
    assert ren.z0 == 50.0
    for i, f in enumerate(freqs):
        one = an.solve_s(freqs=f).to_s(50.0)
        Si = np.array([[complex(ren.S[a, b][i]) for b in range(2)]
                       for a in range(2)])
        CSi = np.array([[complex(ren.CS[a, b][i]) for b in range(2)]
                        for a in range(2)])
        S1 = np.asarray(one.S, dtype=complex)
        CS1 = np.asarray(one.CS, dtype=complex)
        assert np.max(np.abs(Si - S1)) <= 1e-12 * np.max(np.abs(S1)), (f, Si, S1)
        assert np.max(np.abs(CSi - CS1)) <= 1e-12 * np.max(np.abs(CS1)), (f, CSi, CS1)
        assert np.max(np.abs(CS1)) > 0.0
