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
                                   ddd_of_matrix, s_expand)


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


## -- s-expanded (multiroot) form -----------------------------------------

def test_s_expansion_coefficients_are_exact():
    """Coefficients must equal ``Poly(det, s)`` term for term."""
    s, g, c = sympy.symbols('s g c')
    A = sympy.Matrix([[g + s * c, -g], [-g, g + s * c]])
    E = s_expand(A, s)
    ref = list(reversed(sympy.Poly(sympy.expand(A.det()), s).all_coeffs()))
    got = E.eval_coeffs()
    assert len(got) == len(ref)
    for a, b in zip(got, ref):
        assert sympy.simplify(a - b) == 0


@pytest.mark.parametrize('N', [2, 3, 4])
def test_s_expansion_matches_sympy_on_a_real_circuit(N):
    system = bc.rc_ladder(N)
    E = s_expand(system.A, system.s)
    ref = list(reversed(sympy.Poly(sympy.expand(system.A.det()),
                                   system.s).all_coeffs()))
    got = E.eval_coeffs()
    assert len(got) == len(ref)
    for a, b in zip(got, ref):
        assert sympy.simplify(a - b) == 0


@pytest.mark.parametrize('N', [3, 5])
def test_s_expansion_agrees_with_symbolic_poly(N):
    """Same determinant as the fraction-free solver, up to an overall sign.

    ``symbolic_poly`` fixes its own sign convention during elimination, so the
    two agree up to a constant -- which is what a transfer function cancels.
    """
    from pycircuit.circuit.toolkit import symbolic_poly

    system = bc.rc_ladder(N)
    _, den = symbolic_poly.linearsolver_num_den(system.A, system.b)
    E = s_expand(system.A, system.s)
    det = sum(c * system.s ** k for k, c in enumerate(E.eval_coeffs()))
    ratio = sympy.simplify(sympy.cancel(det / den))
    assert not ratio.free_symbols            # a constant, not a function
    assert abs(abs(complex(ratio)) - 1.0) < 1e-12


@pytest.mark.parametrize('N', [4, 8, 12, 16])
def test_s_expansion_respects_the_theorem_1_bound(N):
    """Size must stay within ``q * d * |DDD|`` (Shi & Tan, TCAD 2001, Thm 1).

    ``q`` is the degree in ``s`` and ``d`` the maximum devices at a node, which
    is at most 2 for an MNA matrix in the compact-symbol form.  The point of the
    theorem is that s-expansion costs a *modest factor* over the plain
    determinant diagram, not another exponential.
    """
    system = bc.rc_ladder(N)
    D = ddd_of_matrix(system.A)
    E = s_expand(system.A, system.s)
    assert E.size <= E.degree * 2 * D.size


def test_s_expansion_shares_between_coefficients():
    """The union must be smaller than the coefficients taken separately."""
    system = bc.rc_ladder(8)
    E = s_expand(system.A, system.s)
    separate = sum(E.coefficient(k).size for k in range(E.degree + 1))
    assert E.size < separate


def test_s_expansion_degree_and_trimming():
    s, g, c = sympy.symbols('s g c')
    ## Purely resistive: no powers of s at all.
    A = sympy.Matrix([[g, -g], [-g, 2 * g]])
    assert s_expand(A, s).degree == 0
    ## One capacitor per node: degree 2.
    B = sympy.Matrix([[g + s * c, -g], [-g, g + s * c]])
    assert s_expand(B, s).degree == 2


def test_s_expansion_rejects_non_polynomial_entries():
    s, g, l = sympy.symbols('s g l')
    A = sympy.Matrix([[g + 1 / (s * l), -g], [-g, g]])
    with pytest.raises(ValueError, match='polynomial'):
        s_expand(A, s)


def test_roots_of_requires_numeric_coefficients():
    s, g, c = sympy.symbols('s g c')
    A = sympy.Matrix([[g + s * c, -g], [-g, g + s * c]])
    with pytest.raises(ValueError, match='symbolic'):
        s_expand(A, s).roots_of()


@pytest.mark.parametrize('N,tol', [(4, 1e-12), (8, 1e-11), (12, 1e-10)])
def test_poles_match_a_generalized_eigenvalue_solve(N, tol):
    """Pole *accuracy*, not merely coefficient agreement.

    Coefficients can be exactly right and the poles still wrong, because
    root-finding from an expanded polynomial is ill-conditioned.  The reference
    is the generalized eigenproblem ``G x = -s C x``, which never forms the
    polynomial and so is not subject to the same error.
    """
    import scipy.linalg as la

    system = bc.rc_ladder(N)
    env = dict(system.params)
    got = np.sort_complex(s_expand(system.A, system.s).roots_of(env))

    n = system.dim
    G = np.array([[complex(system.A[i, j].subs({**env, system.s: 0}))
                   for j in range(n)] for i in range(n)])
    C = np.array([[complex(sympy.diff(system.A[i, j], system.s).subs(env))
                   for j in range(n)] for i in range(n)])
    ref = np.sort_complex(np.array(
        [e for e in la.eig(G, -C, right=False) if np.isfinite(e)]))

    assert len(got) == len(ref)
    assert np.max(np.abs(got - ref) / np.abs(ref)) < tol


def test_mfb_s_expansion_is_second_order():
    """The multiple-feedback filter is a two-pole section."""
    system = bc.mfb_filter()
    E = s_expand(system.A, system.s)
    assert E.degree == 2
    poles = E.roots_of(dict(system.params))
    assert len(poles) == 2
    assert all(p.real < 0 for p in poles)          # stable


## -- dominant terms and approximation ------------------------------------

def _spread_env(system, spread):
    """Component values spread geometrically -- the regime pruning needs."""
    env = {}
    for sym in system.A.free_symbols - {system.s}:
        name = str(sym)
        i = int(''.join(ch for ch in name if ch.isdigit()) or 0)
        base = 100.0 if name.startswith('R') else 1e-9
        env[sym] = base * (spread ** i)
    return env


def test_terms_come_out_in_decreasing_magnitude():
    system = bc.rc_ladder(8)
    env = _spread_env(system, 3.0)
    env[system.s] = 1j * 2 * np.pi * 1e3
    mags = []
    for _, values in ddd_of_matrix(system.A).iter_terms(env):
        mags.append(abs(values[0]))
        if len(mags) == 12:
            break
    assert len(mags) == 12
    assert mags == sorted(mags, reverse=True)


def test_all_terms_sum_to_the_exact_value():
    """Enumerating every term must reconstruct the determinant."""
    system = bc.rc_ladder(6)
    env = _spread_env(system, 2.0)
    env[system.s] = 1j * 2 * np.pi * 1e4
    D = ddd_of_matrix(system.A)
    total = sum(v[0] for _, v in D.iter_terms(env))
    exact = complex(D.eval(env))
    assert abs(total - exact) <= 1e-9 * abs(exact)


def test_term_count_matches_enumeration():
    system = bc.rc_ladder(6)
    env = _spread_env(system, 2.0)
    env[system.s] = 1j * 2 * np.pi * 1e4
    D = ddd_of_matrix(system.A)
    assert sum(1 for _ in D.iter_terms(env)) == D.term_count()


@pytest.mark.parametrize('tol', [1e-1, 1e-2, 1e-3])
def test_approximation_meets_its_tolerance(tol):
    system = bc.rc_ladder(8)
    env = _spread_env(system, 5.0)
    env[system.s] = 1j * 2 * np.pi * 1e4
    _, n, err = ddd_of_matrix(system.A).approximate(env, tol=tol,
                                                    max_terms=100000)
    assert err <= tol
    assert n >= 1


def test_pruning_depends_on_parameter_spread():
    """With uniform components there is nothing to prune; with spread, plenty.

    This is the honest shape of the result -- dominant-term approximation is not
    a free win, it is a win exactly when the parameters differ in magnitude.
    """
    system = bc.rc_ladder(10)
    D = ddd_of_matrix(system.A)
    total = D.term_count()

    flat = _spread_env(system, 1.0)                # every component identical
    spread = _spread_env(system, 10.0)
    for e in (flat, spread):
        e[system.s] = 1j * 2 * np.pi * 1e3

    _, n_flat, _ = D.approximate(flat, tol=1e-2, max_terms=100000)
    _, n_spread, _ = D.approximate(spread, tol=1e-2, max_terms=100000)
    assert n_spread < n_flat
    assert n_spread < 0.1 * total                  # a real reduction


def test_corners_keep_terms_that_only_matter_somewhere():
    """Ranking at one nominal point can discard a term that dominates elsewhere.

    Yu & Sechen's argument for pruning over parameter *ranges*: satisfying both
    corners needs strictly more terms than either alone.
    """
    system = bc.rc_ladder(10)
    c = s_expand(system.A, system.s).coefficient(3)
    lo, hi = _spread_env(system, 3.0), _spread_env(system, 10.0)

    _, n_lo, _ = c.approximate(lo, tol=1e-2, max_terms=100000)
    _, n_hi, _ = c.approximate(hi, tol=1e-2, max_terms=100000)
    _, n_both, err = c.approximate([lo, hi], tol=1e-2, max_terms=100000)

    assert n_both >= max(n_lo, n_hi)
    assert err <= 1e-2                             # holds at *both* corners


def test_approximation_needs_numeric_values():
    system = bc.rc_ladder(4)
    with pytest.raises(ValueError, match='symbolic'):
        ddd_of_matrix(system.A).approximate({}, tol=0.1)


def test_s_expanded_approximation_is_per_coefficient():
    """Each power of s is approximated on its own -- the required order."""
    system = bc.rc_ladder(8)
    E = s_expand(system.A, system.s)
    env = _spread_env(system, 5.0)
    result = E.approximate(env, tol=1e-2, max_terms=100000)

    assert len(result) == E.degree + 1
    for k, (expr, n, err) in enumerate(result):
        assert err <= 1e-2, 'coefficient %d missed its tolerance' % k
        assert n <= E.coefficient(k).term_count()


def test_dominant_poles_approach_the_exact_roots_when_separated():
    """Ratios of consecutive coefficients estimate well-separated poles.

    Accuracy improves for the higher poles, which are the better separated --
    the expected behaviour of the estimate, and the reason it is stated as a
    dominant-pole method rather than an exact one.
    """
    system = bc.rc_ladder(10)
    env = _spread_env(system, 3.0)
    E = s_expand(system.A, system.s)

    est = sorted((complex(p) for p in E.dominant_poles(env)), key=abs)
    exact = sorted(E.roots_of(env), key=abs)

    errors = [abs(a - b) / abs(b) for a, b in zip(est, exact)]
    assert errors[0] < 0.2                         # dominant pole, within 20%
    assert errors[3] < errors[0]                   # better where more separated
