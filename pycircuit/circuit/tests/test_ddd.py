# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Tests for the determinant decision diagram core.

Three kinds of check, in decreasing order of how much they would tell us if they
broke:

* the **structural identity** ``|DDD| == n * 2**(n-1)`` for a full matrix under
  row-wise expansion.  This is exact, comes from Shi (TCAS-II 2010), and is the
  one test that distinguishes "sharing works" from "sharing silently does
  nothing" -- a diagram whose minor cache never hits is still *correct*, just
  exponentially large, so correctness tests alone cannot catch it;
* **signs**, which the layered construction fixes as it goes rather than in a
  separate pass, so a sign error would be a construction bug with nowhere else
  to look.  Caught by evaluating numerically against ``numpy``/``sympy``
  determinants on random sparse matrices;
* **circuit-level agreement** with ``symbolic_poly``.
"""

import numpy as np
import pytest
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import (DDDSizeError, ONE, ZERO, ddd_cramer,
                                   ddd_of_matrix)


def _full_matrix(n):
    return sympy.Matrix(n, n, lambda i, j: sympy.Symbol('a%d_%d' % (i, j)))


## -- the structural identity ---------------------------------------------

@pytest.mark.parametrize('n', [1, 2, 3, 4, 5, 6, 7])
def test_full_matrix_has_the_optimal_ddd_size(n):
    """``|DDD| == n * 2**(n-1)`` exactly, for a full matrix expanded row-wise.

    Shi (TCAS-II 2010) proves this is the optimal size.  It holds here because
    minors are memoised on ``(rows, cols)``: after ``k`` eliminations the rows
    are a fixed suffix while the columns are any ``(n-k)``-subset, giving
    ``C(n,k)`` distinct minors each contributing ``n-k`` vertices, and
    ``sum_k C(n,k)*(n-k) == n*2**(n-1)``.
    """
    assert ddd_of_matrix(_full_matrix(n), order='row').size == n * 2 ** (n - 1)


def test_sharing_actually_happens():
    """Distinct minors must be far fewer than the product-term count."""
    n = 6
    D = ddd_of_matrix(_full_matrix(n), order='row')
    assert D.n_minors < 2 ** n            # vs 720 product terms at n=6
    assert D.size == n * 2 ** (n - 1)


## -- correctness and signs -----------------------------------------------

@pytest.mark.parametrize('order', ['row', 'min-degree'])
def test_small_determinants_are_exact(order):
    a, b, c, d = sympy.symbols('a b c d')
    A = sympy.Matrix([[a, b], [c, d]])
    assert sympy.simplify(ddd_of_matrix(A, order=order).eval() - (a * d - b * c)) == 0

    M = _full_matrix(3)
    assert sympy.simplify(ddd_of_matrix(M, order=order).eval() - M.det()) == 0


@pytest.mark.parametrize('seed', range(6))
def test_signs_against_numeric_determinant(seed):
    """Random sparse matrices, both orderings, evaluated against numpy.

    A sign error survives every structural check but shows up immediately here.
    """
    rng = np.random.default_rng(seed)
    n = int(rng.integers(2, 7))
    A = sympy.zeros(n, n)
    env = {}
    for i in range(n):
        for j in range(n):
            if rng.random() < 0.65:
                sym = sympy.Symbol('a%d_%d' % (i, j))
                A[i, j] = sym
                env[sym] = float(rng.uniform(-4, 4))

    An = np.array([[complex(A[i, j].subs(env)) if A[i, j] != 0 else 0j
                    for j in range(n)] for i in range(n)])
    ref = np.linalg.det(An)
    for order in ('row', 'min-degree'):
        got = ddd_of_matrix(A, order=order).eval(env)
        assert abs(got - ref) <= 1e-8 * max(1.0, abs(ref))


def test_orderings_agree():
    """The expansion order changes the graph, never the value."""
    M = _full_matrix(4)
    a = ddd_of_matrix(M, order='row').eval()
    b = ddd_of_matrix(M, order='min-degree').eval()
    assert sympy.simplify(a - b) == 0


def test_partial_substitution_stays_symbolic():
    a, b, c, d = sympy.symbols('a b c d')
    A = sympy.Matrix([[a, b], [c, d]])
    got = ddd_of_matrix(A).eval({b: 0, c: 0})
    assert sympy.simplify(got - a * d) == 0


## -- degenerate shapes ---------------------------------------------------

def test_structurally_singular_matrix_gives_zero():
    a, b = sympy.symbols('a b')
    A = sympy.Matrix([[a, b], [0, 0]])          # zero row
    D = ddd_of_matrix(A)
    assert D.root is ZERO
    assert D.size == 0
    assert D.eval() == 0


def test_scalar_matrix():
    a = sympy.Symbol('a')
    D = ddd_of_matrix(sympy.Matrix([[a]]))
    assert D.size == 1
    assert D.root.one_edge is ONE
    assert sympy.simplify(D.eval() - a) == 0


def test_non_square_is_rejected():
    with pytest.raises(ValueError, match='square'):
        ddd_of_matrix(sympy.Matrix([[1, 2, 3], [4, 5, 6]]))


def test_unknown_order_is_rejected():
    with pytest.raises(ValueError, match='order'):
        ddd_of_matrix(_full_matrix(2), order='nonsense')


## -- rendering -----------------------------------------------------------

def test_to_dot_is_wellformed():
    dot = ddd_of_matrix(_full_matrix(3), order='row').to_dot()
    assert dot.startswith('digraph')
    assert dot.rstrip().endswith('}')
    assert 'style=dashed' in dot                # 0-edges are distinguishable
    assert dot.count('->') == 2 * 12            # two edges per vertex; n=3 -> 12


def test_to_dot_refuses_a_large_graph():
    """The size guard is mandatory: a real circuit must never reach ``dot``."""
    D = ddd_of_matrix(_full_matrix(6), order='row')     # 192 vertices
    with pytest.raises(DDDSizeError, match='above the'):
        D.to_dot(max_vertices=50)
    assert D.to_dot(max_vertices=500)                   # fits when allowed


## -- Cramer solving ------------------------------------------------------

def test_cramer_solves_a_small_system():
    a, b, c, d = sympy.symbols('a b c d')
    A = sympy.Matrix([[a, b], [c, d]])
    den, num = ddd_cramer(A, [1, 0])
    assert sympy.simplify(num[0].eval() - d) == 0
    assert sympy.simplify(num[1].eval() + c) == 0
    assert sympy.simplify(den.eval() - (a * d - b * c)) == 0


@pytest.mark.parametrize('seed', range(4))
def test_cramer_matches_a_numeric_solve(seed):
    rng = np.random.default_rng(100 + seed)
    n = int(rng.integers(2, 6))
    A = sympy.zeros(n, n)
    env = {}
    for i in range(n):
        for j in range(n):
            sym = sympy.Symbol('a%d_%d' % (i, j))
            A[i, j] = sym
            env[sym] = float(rng.uniform(1, 9)) + (n if i == j else 0)
    b = sympy.Matrix([1] + [0] * (n - 1))

    den, num = ddd_cramer(A, b)
    An = np.array([[complex(A[i, j].subs(env)) for j in range(n)]
                   for i in range(n)])
    ref = np.linalg.solve(An, np.array([1.0] + [0.0] * (n - 1)))
    dv = complex(den.eval(env))
    for i in range(n):
        assert abs(complex(num[i].eval(env)) / dv - ref[i]) <= 1e-8


## -- real circuits -------------------------------------------------------

@pytest.mark.parametrize('N', [4, 8, 12])
def test_ladder_transfer_function_matches_a_numeric_solve(N):
    """Driven from a real pycircuit circuit, through the normal MNA path."""
    system = bc.rc_ladder(N)
    out, inp = system.out_index, system.in_index
    _, num = ddd_cramer(system.A, system.b, indices=[out, inp])

    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e6
    n = system.dim
    An = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                   for i in range(n)])
    bn = np.array([complex(system.b[i].subs(env)) for i in range(n)])
    ref = np.linalg.solve(An, bn)
    ref_tf = ref[out] / ref[inp]

    ## The shared determinant cancels in a ratio of two unknowns.
    got = complex(num[out].eval(env)) / complex(num[inp].eval(env))
    assert abs(got - ref_tf) <= 1e-9 * max(1.0, abs(ref_tf))


def test_mfb_with_a_nullor_agrees_with_symbolic_poly():
    """A nullor stamps unit entries and no diagonal -- a different sparsity
    pattern from the passive ladder, and a classic source of cancelling terms.
    """
    from pycircuit.circuit.toolkit import symbolic_poly

    system = bc.mfb_filter()
    num_ref, den_ref = symbolic_poly.linearsolver_num_den(system.A, system.b)
    den, num = ddd_cramer(system.A, system.b, indices=[system.out_index])

    ## symbolic_poly may carry an overall sign/scale convention of its own, so
    ## compare the transfer function rather than the raw parts.
    ratio_ddd = num[system.out_index].eval() / den.eval()
    ratio_ref = num_ref[system.out_index] / den_ref
    assert sympy.simplify(ratio_ddd - ratio_ref) == 0


@pytest.mark.parametrize('N', [4, 8, 12, 16])
def test_ladder_diagram_stays_small(N):
    """Sparsity, not magic: the ladder's diagram must grow linearly.

    Guards the property the whole approach rests on -- if a change makes the
    minor cache miss, this grows explosively while every correctness test still
    passes.
    """
    D = ddd_of_matrix(bc.rc_ladder(N).A)
    assert D.size < 8 * N
