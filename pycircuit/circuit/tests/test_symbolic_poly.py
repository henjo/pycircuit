# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Tests for the experimental symbolic_poly toolkit.

The toolkit must be a drop-in for ``symbolic`` (same results) while using the
fraction-free polynomial-domain solver, and must not alter ``symbolic`` itself.
"""

import numpy as np
import sympy
from sympy import Symbol

import pycircuit.circuit.circuit
from pycircuit.circuit import symbolic, symbolic_poly, numeric
from pycircuit.circuit import SubCircuit, R, C, VS, gnd
from pycircuit.circuit.analysis_ss import AC


def test_toolkit_identity():
    """symbolic_poly is a distinct symbolic toolkit, symbolic left pristine."""
    assert symbolic_poly.symbolic is True
    assert getattr(symbolic_poly, 'poly', False) is True
    assert not getattr(symbolic, 'poly', False)
    ## pristine symbolic still uses LUsolve, not the fraction-free solver
    import inspect
    assert 'LUsolve' in inspect.getsource(symbolic.linearsolver)
    assert 'solve_den' in inspect.getsource(symbolic_poly.linearsolver)


def test_linearsolver_matches_symbolic():
    """Fraction-free solve returns the same solution as LUsolve."""
    s, R0, R1, C0 = sympy.symbols('s R0 R1 C0', positive=True)
    A = np.array([[1/R0 + s*C0, -1/R0],
                  [-1/R0, 1/R0 + 1/R1]], dtype=object)
    b = np.array([1/R0, 0], dtype=object)

    x_ref = symbolic.linearsolver(A, b)
    x_poly = symbolic_poly.linearsolver(A, b)

    for a, c in zip(x_ref, x_poly):
        assert sympy.simplify(a - c) == 0


def test_linearsolver_transcendental_fallback():
    """Entries outside a polynomial domain fall back to LUsolve."""
    w = Symbol('w')
    A = np.array([[sympy.exp(w), 1], [0, 2]], dtype=object)
    b = np.array([1, 1], dtype=object)

    x = symbolic_poly.linearsolver(A, b)
    x_ref = sympy.Matrix([[sympy.exp(w), 1], [0, 2]]).LUsolve(sympy.Matrix([[1], [1]]))
    for a, c in zip(x, x_ref):
        assert sympy.simplify(a - c) == 0


def test_ac_transfer_function():
    """End-to-end symbolic AC via symbolic_poly on a float-valued RC circuit."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    cir = SubCircuit()
    cir['VS'] = VS('in', gnd, vac=1)
    cir['R'] = R('in', 'out', r=1e3)
    cir['C'] = C('out', gnd, c=1e-9)

    s = Symbol('s')
    H = sympy.cancel(AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True).v('out', gnd))

    ## 1 / (1 + s R C) with R C = 1e-6
    assert sympy.simplify(H - 1/(1 + s*sympy.Rational(1, 10**6))) == 0
