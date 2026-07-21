# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Tests for the experimental symbolic_poly toolkit and its features:

* fraction-free ``linearsolver`` / ``linearsolver_num_den`` (a drop-in for the
  stock symbolic solver, but returning an exact ``N/D`` form),
* the transfer-function facade (``res.tf``, poles/zeros/dcgain/bode),
* the ``use_toolkit`` context manager replacing the ``default_toolkit`` global.

``symbolic`` itself must be left pristine.
"""

import inspect

import numpy as np
import sympy
from sympy import Symbol

import pycircuit.circuit.circuit as ccirc
from pycircuit.circuit import (symbolic, symbolic_poly, numeric,
                               use_toolkit, set_default_toolkit)
from pycircuit.circuit import SubCircuit, R, C, VS, gnd
from pycircuit.circuit.analysis_ss import AC, CircuitResultAC, CircuitResultACPoly
from pycircuit.circuit.transferfunction import TransferFunction


# --------------------------------------------------------------------------
# Toolkit identity / capabilities
# --------------------------------------------------------------------------
def test_toolkit_identity():
    """symbolic_poly is a distinct symbolic toolkit; symbolic is left pristine."""
    assert symbolic_poly.symbolic is True
    assert symbolic_poly.poly is True
    assert not getattr(symbolic, 'poly', False)
    assert not getattr(numeric, 'poly', False)

    ## capability introspection
    assert symbolic_poly.supports('num_den') is True
    assert symbolic.supports('num_den') is False
    assert numeric.supports('num_den') is False

    ## pristine symbolic still uses LUsolve; poly uses the fraction-free solver
    assert 'LUsolve' in inspect.getsource(symbolic.linearsolver)
    assert 'solve_den' in inspect.getsource(symbolic_poly.linearsolver_num_den)


# --------------------------------------------------------------------------
# linearsolver / linearsolver_num_den
# --------------------------------------------------------------------------
def _rc_system():
    s, R0, R1, C0 = sympy.symbols('s R0 R1 C0', positive=True)
    A = np.array([[1/R0 + s*C0, -1/R0],
                  [-1/R0, 1/R0 + 1/R1]], dtype=object)
    b = np.array([1/R0, 0], dtype=object)
    return A, b


def test_linearsolver_matches_symbolic():
    """Fraction-free solve returns the same solution as LUsolve."""
    A, b = _rc_system()
    x_ref = symbolic.linearsolver(A, b)
    x_poly = symbolic_poly.linearsolver(A, b)
    assert len(x_poly) == len(x_ref)
    for a, c in zip(x_ref, x_poly):
        assert sympy.simplify(a - c) == 0


def test_linearsolver_num_den():
    """num/den form divides back to the same solution as the plain solver."""
    A, b = _rc_system()
    num, den = symbolic_poly.linearsolver_num_den(A, b)
    x = symbolic_poly.linearsolver(A, b)
    x_ref = symbolic.linearsolver(A, b)

    assert len(num) == len(x)
    for ni, xi, xref in zip(num, x, x_ref):
        assert sympy.simplify(ni/den - xi) == 0        # x_i == num_i/den
        assert sympy.simplify(ni/den - xref) == 0      # and == LUsolve


def test_linearsolver_transcendental():
    """Transcendental entries still solve correctly (via the EX domain or the
    LUsolve fallback -- either way the result must match LUsolve)."""
    w = Symbol('w')
    A = np.array([[sympy.exp(w), 1], [0, 2]], dtype=object)
    b = np.array([1, 1], dtype=object)

    num, den = symbolic_poly.linearsolver_num_den(A, b)
    x = symbolic_poly.linearsolver(A, b)
    x_ref = sympy.Matrix([[sympy.exp(w), 1], [0, 2]]).LUsolve(sympy.Matrix([[1], [1]]))
    for ni, xi, xref in zip(num, x, x_ref):
        assert sympy.simplify(ni/den - xi) == 0
        assert sympy.simplify(xi - xref) == 0


def test_linearsolver_fallback_path():
    """The LUsolve fallback is exercised when solve_den raises."""
    import unittest.mock as mock
    A, b = _rc_system()
    x_ref = symbolic.linearsolver(A, b)
    ## Force the polynomial-domain path to raise so the except branch runs.
    with mock.patch('pycircuit.circuit.toolkit.DomainMatrix.from_Matrix',
                    side_effect=TypeError('forced')):
        x = symbolic_poly.linearsolver(A, b)
    for a, c in zip(x, x_ref):
        assert sympy.simplify(a - c) == 0


# --------------------------------------------------------------------------
# Transfer function facade
# --------------------------------------------------------------------------
def _rc_lowpass(Rval, Cval):
    """First-order RC low-pass built under symbolic_poly, output at 'out'."""
    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        cir['VS'] = VS('in', gnd, vac=1)
        cir['R'] = R('in', 'out', r=Rval)
        cir['C'] = C('out', gnd, c=Cval)
    return cir


def test_ac_returns_poly_result():
    """AC with symbolic_poly yields a transfer-function-capable result."""
    Rv, Cv = sympy.symbols('R C', positive=True)
    cir = _rc_lowpass(Rv, Cv)
    s = Symbol('s', imaginary=True)
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)
    assert isinstance(res, CircuitResultACPoly)
    ## v() still works and equals num/den
    assert sympy.simplify(res.v('out', gnd) - 1/(1 + s*Rv*Cv)) == 0


def test_tf_poles_zeros_dcgain():
    Rv, Cv = sympy.symbols('R C', positive=True)
    cir = _rc_lowpass(Rv, Cv)
    s = Symbol('s', imaginary=True)
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)

    H = res.tf('out', gnd)
    assert isinstance(H, TransferFunction)
    assert sympy.simplify(H.canonical() - 1/(1 + s*Rv*Cv)) == 0
    assert H.poles() == {-1/(Rv*Cv): 1}
    assert H.zeros() == {}
    assert sympy.simplify(H.dcgain() - 1) == 0

    ## the shared-denominator poles on the result agree with the tf poles
    assert res.poles() == H.poles()


def test_tf_highpass_has_zero_at_origin():
    """Series-C / shunt-R high-pass has a zero at s=0 and one pole."""
    Rv, Cv = sympy.symbols('R C', positive=True)
    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        cir['VS'] = VS('in', gnd, vac=1)
        cir['C'] = C('in', 'out', c=Cv)
        cir['R'] = R('out', gnd, r=Rv)
    s = Symbol('s', imaginary=True)
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)

    H = res.tf('out', gnd)
    assert sympy.simplify(H.canonical() - s*Rv*Cv/(1 + s*Rv*Cv)) == 0
    assert H.zeros() == {0: 1}
    assert H.poles() == {-1/(Rv*Cv): 1}


def test_tf_bode_numeric_matches():
    """bode() lambdifies to a fast evaluator matching direct substitution."""
    cir = _rc_lowpass(1e3, 1e-9)          # numeric component values
    s = Symbol('s', imaginary=True)
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)

    H = res.tf('out', gnd)
    f = H.bode()                           # num/den are in s only
    jw = 2j*np.pi*np.array([1e3, 1e6, 1e9])
    got = f(jw)
    expected = 1/(1 + jw*1e3*1e-9)
    assert np.allclose(np.asarray(got, dtype=complex), expected)


def test_tf_i_current_transfer():
    """tf_i gives the current transfer function and matches res.i()."""
    Rv, Cv = sympy.symbols('R C', positive=True)
    cir = _rc_lowpass(Rv, Cv)              # VS -> R -> C to gnd
    s = Symbol('s', imaginary=True)
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)

    Hi = res.tf_i('R.plus')                # current out of the source through R
    assert sympy.simplify(Hi.canonical() - s*Cv/(1 + s*Rv*Cv)) == 0
    assert Hi.poles() == {-1/(Rv*Cv): 1}
    assert Hi.zeros() == {0: 1}            # zero at DC (no current at s=0)

    ## matches the existing current extraction
    i_ref = cir.extract_i(res.x, 'R.plus', xdot=res.xdot,
                          linearized=True, xdcop=res.xdcop)
    assert sympy.simplify(Hi.expr() - i_ref) == 0


def test_input_impedance_via_voltage_tf():
    """A unit current source makes the node-voltage tf the input impedance."""
    from pycircuit.circuit import IS
    Rv, Cv = sympy.symbols('R C', positive=True)
    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        cir['IS'] = IS(gnd, 'in', iac=1)   # 1 A into 'in'
        cir['R'] = R('in', gnd, r=Rv)
        cir['C'] = C('in', gnd, c=Cv)      # R || C
    s = Symbol('s', imaginary=True)
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)

    Zin = res.tf('in', gnd)                # v/1 = Zin
    assert sympy.simplify(Zin.canonical() - Rv/(1 + s*Rv*Cv)) == 0


def test_numeric_poles_and_zeros():
    """numeric=True finds poles/zeros with numpy.roots for numeric circuits."""
    R0, C0 = 1e3, 1e-9                       # numeric values -> pole at -1/RC
    cir = _rc_lowpass(R0, C0)
    s = Symbol('s', imaginary=True)
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)

    poles = res.poles(numeric=True)
    assert len(poles) == 1
    assert np.allclose(poles, [-1/(R0*C0)])

    ## tf-level numeric poles/zeros agree; low-pass has no finite zero
    H = res.tf('out', gnd)
    assert np.allclose(H.poles(numeric=True), [-1/(R0*C0)])
    assert len(H.zeros(numeric=True)) == 0


def test_numeric_poles_scale():
    """A larger numeric ladder yields the right number of poles quickly."""
    s = Symbol('s', imaginary=True)
    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        cir['VS'] = VS('in', gnd, vac=1)
        prev = 'in'
        for i in range(6):
            cir['R%d' % i] = R(prev, 'n%d' % i, r=1e3)
            cir['C%d' % i] = C('n%d' % i, gnd, c=1e-9)
            prev = 'n%d' % i
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)
    poles = res.poles(numeric=True)
    assert len(poles) == 6
    ## a passive RC network is stable: all poles on the negative real axis
    assert np.all(np.real(poles) < 0)


def test_numeric_roots_require_numeric_coefficients():
    """numeric=True on a symbolic-valued circuit raises a helpful error."""
    Rv, Cv = sympy.symbols('R C', positive=True)
    cir = _rc_lowpass(Rv, Cv)
    s = Symbol('s', imaginary=True)
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)
    import pytest
    with pytest.raises(ValueError):
        res.poles(numeric=True)
    with pytest.raises(ValueError):
        res.tf('out', gnd).poles(numeric=True)


def test_tf_frequencyresponse():
    """frequencyresponse() evaluates H over a frequency array in Hz."""
    R0, C0 = 1e3, 1e-9
    cir = _rc_lowpass(R0, C0)
    s = Symbol('s', imaginary=True)
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)

    freqs = np.array([1e3, 1e6, 1e9])
    got = res.tf('out', gnd).frequencyresponse(freqs)
    expected = 1/(1 + 2j*np.pi*freqs*R0*C0)
    assert np.allclose(got, expected)


def test_noise_shared_denominator():
    """symbolic_poly noise gives the same PSD as symbolic, but far more compact.

    The generic form builds zm^T CY conj(zm) as an O(n^2) sum of divided
    rationals; the polynomial toolkit keeps it as N^T CY conj(N) / (D conj(D)),
    a single rational with a compact intermediate.
    """
    from pycircuit.circuit.analysis_ss import Noise

    def rc_noise(tk):
        with use_toolkit(tk):
            c = SubCircuit()
            Rv, Cv = sympy.symbols('R C', positive=True)
            c['vs'] = VS(1, gnd, vac=1)
            c['R'] = R(1, 2, r=Rv)
            c['C1'] = C(2, gnd, c=Cv)
            noise = Noise(c, inputsrc='vs', outputnodes=('2', gnd), toolkit=tk)
            return noise.solve(Symbol('s', imaginary=True), complexfreq=True)['Svnout']

    S_sym = rc_noise(symbolic)
    S_poly = rc_noise(symbolic_poly)

    ## same value
    assert sympy.simplify(S_sym - S_poly) == 0
    ## compact: the shared-denominator intermediate is much smaller
    assert sympy.count_ops(S_poly) < sympy.count_ops(S_sym)
    ## imaginary s resolves conjugation -- no lingering conjugate() terms
    assert not S_poly.has(sympy.conjugate)


def test_toolkit_noise_psd_matches():
    """The noise_psd toolkit methods agree between symbolic and symbolic_poly."""
    s = Symbol('s', imaginary=True)
    Rv, Cv, k, T = sympy.symbols('R C k T', positive=True)
    ## 2-node reciprocal system with a diagonal noise matrix
    Y = np.array([[1/Rv + s*Cv, -1/Rv], [-1/Rv, 1/Rv]], dtype=object)
    u = np.array([1, 0], dtype=object)
    CY = np.array([[4*k*T/Rv, 0], [0, 4*k*T/Rv]], dtype=object)

    zm_s, psd_s = symbolic.noise_psd(Y, u, CY, s)
    zm_p, psd_p = symbolic_poly.noise_psd(Y, u, CY, s)
    assert sympy.simplify(psd_s - psd_p) == 0


def test_symbolic_result_has_no_tf():
    """The stock symbolic toolkit produces a plain result without tf()."""
    Rv, Cv = sympy.symbols('R C', positive=True)
    with use_toolkit(symbolic):
        cir = SubCircuit()
        cir['VS'] = VS('in', gnd, vac=1)
        cir['R'] = R('in', 'out', r=Rv)
        cir['C'] = C('out', gnd, c=Cv)
    s = Symbol('s')
    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)
    assert isinstance(res, CircuitResultAC)
    assert not isinstance(res, CircuitResultACPoly)
    assert not hasattr(res, 'tf')


# --------------------------------------------------------------------------
# use_toolkit / set_default_toolkit
# --------------------------------------------------------------------------
def test_use_toolkit_sets_and_restores():
    before = ccirc.default_toolkit
    with use_toolkit(symbolic_poly):
        assert ccirc.default_toolkit is symbolic_poly
        with use_toolkit(symbolic):                 # nesting
            assert ccirc.default_toolkit is symbolic
        assert ccirc.default_toolkit is symbolic_poly
    assert ccirc.default_toolkit is before


def test_use_toolkit_restores_on_exception():
    before = ccirc.default_toolkit
    try:
        with use_toolkit(symbolic_poly):
            raise ValueError('boom')
    except ValueError:
        pass
    assert ccirc.default_toolkit is before


def test_use_toolkit_affects_construction():
    """A circuit built inside the context uses the scoped toolkit."""
    before = ccirc.default_toolkit
    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        assert cir.toolkit is symbolic_poly
    assert ccirc.default_toolkit is before


def test_set_default_toolkit():
    before = ccirc.default_toolkit
    try:
        set_default_toolkit(symbolic)
        assert ccirc.default_toolkit is symbolic
    finally:
        set_default_toolkit(before)
    assert ccirc.default_toolkit is before
