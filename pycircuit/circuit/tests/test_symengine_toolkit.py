"""Tests for the symengine-backed symbolic toolkit (GiNaC-style compiled solve)."""
import pytest
import sympy

pytest.importorskip("symengine")

from pycircuit.circuit import symbolic_poly, SubCircuit, R, C, VS, gnd
from pycircuit.circuit.toolkit import symengine_toolkit
from pycircuit.circuit.analysis_ss import AC


def _rc_poles(tk):
    s = sympy.Symbol('s', imaginary=True)
    cir = SubCircuit(toolkit=tk)
    cir['VS'] = VS('in', gnd, vac=1)
    cir['R'] = R('in', 'out', r=sympy.Symbol('R', positive=True))
    cir['C'] = C('out', gnd, c=sympy.Symbol('C', positive=True))
    res = AC(cir, toolkit=tk).solve(s, complexfreq=True)
    return res.tf('out', gnd).poles()


def test_symengine_toolkit_available():
    assert symengine_toolkit is not None


def test_symengine_rc_poles_match_symbolic_poly():
    R_, C_ = sympy.symbols('R C', positive=True)
    poles = _rc_poles(symengine_toolkit)
    # Correct analytic pole of the RC low-pass, with symbol assumptions preserved
    assert poles == {-1 / (C_ * R_): 1}
    # And identical to the pure-sympy fraction-free toolkit.
    assert poles == _rc_poles(symbolic_poly)


def test_symengine_linearsolver_num_den_matches():
    import numpy as np
    s = sympy.Symbol('s')
    Rs, Cs = sympy.symbols('R C', positive=True)
    A = np.array([[1 / Rs + s * Cs, -1 / Rs], [-1 / Rs, 1 / Rs]], dtype=object)
    b = np.array([0, 1], dtype=object)

    n_poly, d_poly = symbolic_poly.linearsolver_num_den(A, b)
    n_se, d_se = symengine_toolkit.linearsolver_num_den(A, b)

    for i in range(2):
        r_poly = sympy.simplify(n_poly[i] / d_poly)
        r_se = sympy.simplify(n_se[i] / d_se)
        assert sympy.simplify(r_poly - r_se) == 0
