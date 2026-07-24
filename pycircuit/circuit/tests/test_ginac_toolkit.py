"""Tests for the experimental GiNaC-backed symbolic toolkit.

Skipped unless the compiled ``_ginac_ext`` extension is built (see
``pycircuit/circuit/build_ginac_ext.sh``).
"""
import numpy as np
import pytest
import sympy

pytest.importorskip("pycircuit.circuit._ginac_ext")

from pycircuit.circuit import symbolic_poly, SubCircuit, R, C, VS, gnd
from pycircuit.circuit.toolkit import ginac_toolkit
from pycircuit.circuit.analysis_ss import AC


def _rc_poles(tk):
    s = sympy.Symbol('s', imaginary=True)
    cir = SubCircuit(toolkit=tk)
    cir['VS'] = VS('in', gnd, vac=1)
    cir['R'] = R('in', 'out', r=sympy.Symbol('R', positive=True))
    cir['C'] = C('out', gnd, c=sympy.Symbol('C', positive=True))
    res = AC(cir, toolkit=tk).solve(s, complexfreq=True)
    return res.tf('out', gnd).poles()


def test_ginac_toolkit_available():
    assert ginac_toolkit is not None


def test_ginac_rc_poles_match_symbolic_poly():
    R_, C_ = sympy.symbols('R C', positive=True)
    poles = _rc_poles(ginac_toolkit)
    assert poles == {-1 / (C_ * R_): 1}
    assert poles == _rc_poles(symbolic_poly)


def test_ginac_linearsolver_num_den_matches():
    s = sympy.Symbol('s')
    Rs, Cs = sympy.symbols('R C', positive=True)
    A = np.array([[1 / Rs + s * Cs, -1 / Rs], [-1 / Rs, 1 / Rs]], dtype=object)
    b = np.array([0, 1], dtype=object)

    n_poly, d_poly = symbolic_poly.linearsolver_num_den(A, b)
    n_g, d_g = ginac_toolkit.linearsolver_num_den(A, b)

    for i in range(2):
        assert sympy.simplify(n_poly[i] / d_poly - n_g[i] / d_g) == 0


def test_ginac_large_system_falls_back():
    # Above ginac_max_dim the toolkit must fall back to sympy (guards the CLN
    # coefficient-explosion hang) and still return a correct result.
    n = ginac_toolkit.ginac_max_dim + 3
    s = sympy.Symbol('s')
    A = np.array(sympy.eye(n).tolist(), dtype=object)
    A[0, 0] = 1 + s
    b = np.zeros(n, dtype=object)
    b[0] = 1
    num, den = ginac_toolkit.linearsolver_num_den(A, b)
    assert sympy.simplify(num[0] / den - 1 / (1 + s)) == 0
