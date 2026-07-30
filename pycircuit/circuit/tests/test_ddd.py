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
from pycircuit.circuit.ddd import (DDDSizeError, HierarchicalDDD,
                                   NumericTerminal, ONE, ZERO,
                                   ddd_cofactor_solve, ddd_cramer,
                                   ddd_of_matrix, eval_roots,
                                   reverse_cuthill_mckee,
                                   suppression_order, hierarchical_solve,
                                   noise_psd_reduced, determinant_sensitivity,
                                   adjoint_sensitivities, s_expand)


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


## -- shared cofactor family (the structure noise needs) ------------------

@pytest.mark.parametrize('N', [4, 6, 8])
def test_cofactor_numerators_match_cramer(N):
    """Same values, built a cheaper way."""
    system = bc.rc_ladder(N)
    _, cramer = ddd_cramer(system.A, system.b)
    _, cofactor = ddd_cofactor_solve(system.A, system.b)
    for i in range(system.dim):
        assert sympy.simplify(cofactor[i].eval() - cramer[i].eval()) == 0


@pytest.mark.parametrize('N', [6, 10, 14])
def test_the_whole_family_costs_little_more_than_the_determinant(N):
    """Shi & Tan's claim, checked on our own circuits.

    Every transimpedance plus the determinant, sharing one memo table, must cost
    far less than building each numerator as its own diagram -- and stay within a
    small factor of the determinant alone.  That is what makes noise analysis,
    which needs one transfer function per noise source, affordable.
    """
    system = bc.rc_ladder(N)
    family, _ = ddd_cofactor_solve(system.A, system.b)

    den, separate_nums = ddd_cramer(system.A, system.b)
    separate = den.size + sum(v.size for v in separate_nums.values())

    assert family.size < separate / 2
    assert family.size < 5 * family.denominator.size


def test_sparse_right_hand_side_costs_few_cofactors():
    """Only nonzero entries of b contribute, so a one-source drive is cheap."""
    system = bc.rc_ladder(6)
    _, nums = ddd_cofactor_solve(system.A, system.b)
    nonzeros = sum(1 for k in range(system.dim) if system.b[k] != 0)
    for i in range(system.dim):
        assert len(nums[i].parts) <= nonzeros


def test_family_rejects_non_square():
    with pytest.raises(ValueError, match='square'):
        ddd_cofactor_solve(sympy.Matrix([[1, 2, 3], [4, 5, 6]]), [1, 0, 0])


## -- noise ----------------------------------------------------------------

def test_noise_psd_matches_symbolic_poly():
    """Identical PSD, despite a different construction.

    The shared denominator makes the PSD independent of any overall sign
    convention, so the two agree exactly rather than up to a factor.
    """
    from pycircuit.circuit.toolkit import ddd_toolkit, symbolic_poly

    s = sympy.Symbol('s', imaginary=True)
    R1, R2, C1, kT = sympy.symbols('R1 R2 C1 kT', positive=True)
    Y = sympy.Matrix([[1 / R1 + 1 / R2 + s * C1, -1 / R2], [-1 / R2, 1 / R2]])
    u = sympy.Matrix([0, 1])
    CY = sympy.Matrix([[4 * kT / R1, 0], [0, 4 * kT / R2]])

    z_ref, psd_ref = symbolic_poly.noise_psd(Y, u, CY, s)
    z_got, psd_got = ddd_toolkit.noise_psd(Y, u, CY, s)

    assert sympy.simplify(sympy.cancel(psd_got - psd_ref)) == 0
    for a, b in zip(z_got, z_ref):
        assert sympy.simplify(a - b) == 0


def test_noise_analysis_through_the_toolkit():
    """End to end: the Noise analysis itself, driven with ``ddd_toolkit``."""
    import pycircuit.circuit.circuit
    from pycircuit.circuit import SubCircuit, R, VS, gnd
    from pycircuit.circuit.analysis_ss import Noise
    from pycircuit.circuit.toolkit import ddd_toolkit, symbolic

    R1v, R2v, Vv = sympy.symbols('R1 R2 V', real=True, positive=True)
    s = sympy.Symbol('s')

    saved = pycircuit.circuit.circuit.default_toolkit
    try:
        pycircuit.circuit.circuit.default_toolkit = symbolic
        cir = SubCircuit(toolkit=symbolic)
        cir['vs'] = VS(1, gnd, vac=Vv)
        cir['R1'] = R(1, 2, r=R1v)
        cir['R2'] = R(2, gnd, r=R2v)

        noise = Noise(cir, inputsrc='vs', outputnodes=('2', gnd),
                      toolkit=ddd_toolkit)
        res = noise.solve(s, complexfreq=True)
    finally:
        pycircuit.circuit.circuit.default_toolkit = saved

    kT = noise.toolkit.kboltzmann * noise.par.epar.T
    assert sympy.simplify(res['Svnout'] - 4 * R1v * R2v * kT / (R1v + R2v)) == 0
    assert sympy.simplify(res['gain'] - R2v / (R1v + R2v)) == 0


## -- numeric terminals (semi-symbolic / MTDDD) ---------------------------

def _semi(N, n_symbolic=2):
    system = bc.rc_ladder_semi_symbolic(N, n_symbolic=n_symbolic)
    return system, system.A.free_symbols - {system.s}


@pytest.mark.parametrize('N', [8, 12, 16])
def test_numeric_terminals_do_not_change_the_value(N):
    """Collapsing a parameter-free minor is an accounting change, not a new answer."""
    system, keep = _semi(N)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e5

    plain = ddd_of_matrix(system.A)
    collapsed = ddd_of_matrix(system.A, keep_symbolic=keep, collapse_max_dim=4)

    a, b = complex(plain.eval(env)), complex(collapsed.eval(env))
    assert abs(a - b) <= 1e-9 * abs(a)


def test_numeric_terminals_are_opt_in():
    """Default construction must produce no collapsed terminals."""
    system, keep = _semi(12)
    plain = ddd_of_matrix(system.A)
    assert plain.n_collapsed == 0
    assert not any(isinstance(v.one_edge, NumericTerminal)
                   or isinstance(v.zero_edge, NumericTerminal)
                   for v in plain.vertices())


def test_numeric_terminals_shrink_the_diagram():
    system, keep = _semi(16)
    plain = ddd_of_matrix(system.A)
    collapsed = ddd_of_matrix(system.A, keep_symbolic=keep, collapse_max_dim=5)
    assert collapsed.n_collapsed > 0
    assert collapsed.size < plain.size


def test_collapse_respects_its_size_cap():
    """A larger cap collapses more, because evaluating a determinant is costly."""
    system, keep = _semi(16)
    sizes = [ddd_of_matrix(system.A, keep_symbolic=keep,
                           collapse_max_dim=cap).size for cap in (2, 3, 4, 5)]
    assert sizes == sorted(sizes, reverse=True)


def test_collapsed_terminal_still_depends_on_frequency():
    """A collapsed minor is a polynomial in s, not a constant.

    Its value has to be substituted like any other payload -- a terminal is not
    automatically a number.
    """
    system, keep = _semi(10)
    D = ddd_of_matrix(system.A, keep_symbolic=keep, collapse_max_dim=4)
    env = dict(system.params)
    lo = complex(D.eval(dict(env) | {system.s: 1j * 2 * np.pi * 1e2}))
    hi = complex(D.eval(dict(env) | {system.s: 1j * 2 * np.pi * 1e8}))
    assert abs(lo - hi) > 0


def test_collapse_keeps_the_symbols_it_was_told_to():
    system, keep = _semi(10)
    D = ddd_of_matrix(system.A, keep_symbolic=keep, collapse_max_dim=4)
    ## Substituting only the frequency must leave the kept parameters behind.
    value = D.eval({system.s: 1j})
    assert value.free_symbols == keep


def test_collapsed_diagram_renders():
    system, keep = _semi(8)
    dot = ddd_of_matrix(system.A, keep_symbolic=keep,
                        collapse_max_dim=4).to_dot(max_vertices=500)
    assert dot.startswith('digraph')
    assert dot.rstrip().endswith('}')


## -- hierarchy ------------------------------------------------------------

@pytest.mark.parametrize('internal', [(1, 2), (1, 2, 3), (2, 3, 4), (1, 2, 3, 4)])
def test_hierarchical_determinant_matches_the_flat_one(internal):
    """Suppressing an internal block must not change the determinant."""
    system = bc.rc_ladder(6)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e4

    flat = complex(ddd_of_matrix(system.A).eval(env))
    hier = complex(HierarchicalDDD(system.A, internal).eval(env))
    assert abs(flat - hier) <= 1e-9 * abs(flat)


@pytest.mark.parametrize('n', [6, 8, 10])
def test_hierarchy_beats_a_flat_dense_determinant(n):
    """Where the flat diagram is exponential, splitting it is a large win.

    A ladder gains nothing -- its flat diagram is already linear -- so the
    benefit is stated on the case that actually has one: total size becomes a
    sum over blocks instead of a product.
    """
    system = bc.dense_symbolic_matrix(n)
    env = dict(system.params)

    flat = ddd_of_matrix(system.A)
    hier = HierarchicalDDD(system.A, tuple(range(n // 2)))

    assert abs(complex(hier.eval(env)) - complex(flat.eval(env))) \
        <= 1e-7 * abs(complex(flat.eval(env)))
    assert hier.size < flat.size / 4


def test_hierarchy_on_a_circuit_with_a_nullor():
    system = bc.mfb_filter()
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e4
    flat = complex(ddd_of_matrix(system.A).eval(env))
    hier = complex(HierarchicalDDD(system.A, (1,)).eval(env))
    assert abs(flat - hier) <= 1e-9 * abs(flat)


def test_hierarchy_rejects_a_degenerate_partition():
    A = bc.rc_ladder(4).A
    with pytest.raises(ValueError, match='at least one block'):
        HierarchicalDDD(A, ())
    with pytest.raises(ValueError, match='at least one unknown'):
        HierarchicalDDD(A, tuple(range(A.rows)))
    with pytest.raises(ValueError, match='square'):
        HierarchicalDDD(sympy.Matrix([[1, 2, 3], [4, 5, 6]]), (0,))


def test_hierarchy_size_is_a_sum_of_its_levels():
    system = bc.dense_symbolic_matrix(8)
    hier = HierarchicalDDD(system.A, (0, 1, 2, 3))
    assert hier.size == sum(l['family'].size for l in hier.levels) + hier.top.size


def test_eval_roots_matches_individual_evaluation():
    """The shared evaluation pass must agree with evaluating one at a time."""
    system = bc.rc_ladder(5)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e4
    family, nums = ddd_cofactor_solve(system.A, system.b)

    roots = [family.denominator.root]
    for combination in nums.values():
        roots.extend(combination.roots())
    shared = eval_roots(roots, env)

    assert abs(complex(shared[id(family.denominator.root)])
               - complex(family.denominator.eval(env))) < 1e-9


## -- paper calibration: a 741-class amplifier ----------------------------

def _opamp_env(system, freq=1e3):
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * freq
    return env


def test_opamp_is_a_real_circuit_that_solves():
    """Sanity before calibration: it must behave like an amplifier."""
    system = bc.opamp_741_like()
    assert system.cir is not None
    assert system.dim == 14

    env = _opamp_env(system)
    n = system.dim
    A = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                  for i in range(n)])
    b = np.array([complex(system.b[i].subs(env)) for i in range(n)])
    x = np.linalg.solve(A, b)
    gain = abs(x[system.out_index] / x[system.in_index])
    assert 1e3 < gain < 1e7          # a 741 is ~100 dB open loop


def test_opamp_determinant_agrees_with_numpy():
    system = bc.opamp_741_like()
    env = _opamp_env(system)
    n = system.dim
    A = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                  for i in range(n)])
    got = complex(ddd_of_matrix(system.A).eval(env))
    ref = np.linalg.det(A)
    assert abs(got - ref) <= 1e-6 * abs(ref)


def test_diagram_size_is_independent_of_how_many_symbols_there_are():
    """The property that separates a diagram from an expression.

    A DDD's shape is fixed by the matrix *sparsity pattern*; the symbols ride
    along as payloads.  So making ten device parameters symbolic instead of one
    changes nothing about its size -- where an expanded expression would grow
    explosively.
    """
    sizes = set()
    for devices in ((), ('q1', 'q2', 'q17'),
                    ('q1', 'q2', 'q3', 'q4', 'q5', 'q6', 'q16', 'q17',
                     'q23', 'q14')):
        system = bc.opamp_741_like(symbolic_devices=devices)
        sizes.add(ddd_of_matrix(system.A).size)
    assert len(sizes) == 1


def test_opamp_symbolic_devices_really_are_symbolic():
    system = bc.opamp_741_like(symbolic_devices=('q1', 'q17'))
    names = {str(sym) for sym in system.A.free_symbols}
    assert 'gm_q1' in names and 'gm_q17' in names
    assert 'gm_q2' not in names


def test_hierarchy_agrees_on_the_opamp():
    """Correct on a real amplifier, even though it saves little there."""
    system = bc.opamp_741_like()
    env = _opamp_env(system)
    flat = complex(ddd_of_matrix(system.A).eval(env))
    ## Suppress the input stage: e1, e2, c3, e5, e6.
    hier = complex(HierarchicalDDD(system.A, (2, 3, 4, 6, 7)).eval(env))
    assert abs(flat - hier) <= 1e-7 * abs(flat)


## -- the µA741 calibration circuit ---------------------------------------

def test_ua741_structure_is_close_to_the_published_circuit():
    """Tan & Shi report a 24x24 matrix with 89 nonzeros for the µA741."""
    system = bc.ua741()
    assert 24 <= system.dim <= 30
    nonzeros = sum(1 for i in range(system.dim) for j in range(system.dim)
                   if system.A[i, j] != 0)
    assert 80 <= nonzeros <= 130


def test_ua741_behaves_like_an_opamp():
    """Guard the topology: high DC gain rolling off at 20 dB/decade.

    A wiring error shows up here and nowhere else -- an earlier version tied
    both input emitters to one node, which is electrically a short between the
    two signal paths, and it merely reduced the gain to 3 rather than failing.
    """
    system = bc.ua741()
    n = system.dim

    def gain(freq):
        env = dict(system.params)
        env[system.s] = 1j * 2 * np.pi * freq
        A = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                      for i in range(n)])
        b = np.array([complex(system.b[i].subs(env)) for i in range(n)])
        x = np.linalg.solve(A, b)
        return abs(x[system.out_index] / x[system.in_index])

    assert gain(1.0) > 1e4                       # high open-loop gain
    ## One decade of frequency costs about one decade of gain once compensated.
    assert 5 < gain(1e2) / gain(1e3) < 20
    assert gain(1e5) < gain(1e3)


def test_ua741_determinant_agrees_with_numpy():
    system = bc.ua741()
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e3
    n = system.dim
    A = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                  for i in range(n)])
    got = complex(ddd_of_matrix(system.A).eval(env))
    ref = np.linalg.det(A)
    assert abs(got - ref) <= 1e-6 * abs(ref)


def test_ua741_diagram_is_in_the_published_regime():
    """Thousands of vertices standing for millions of terms."""
    D = ddd_of_matrix(bc.ua741().A)
    assert 200 < D.size < 20000              # theirs: 6654
    assert D.term_count() > 1e5              # theirs: 119011


def test_min_degree_ordering_beats_row_wise_on_the_ua741():
    """Ordering is worth a factor here, which matters when comparing sizes.

    Their 2000 figure predates the expansion-ordering scheme used here, and Shi
    (ICCAD 2010) later measured the older Greedy-Labeling order as dramatically
    worse -- so a smaller vertex count than the paper's is expected rather than
    suspicious.
    """
    A = bc.ua741().A
    assert ddd_of_matrix(A, order='min-degree').size \
        < ddd_of_matrix(A, order='row').size


def test_singular_internal_block_is_reported_clearly():
    """Suppression inverts the block, so a singular partition must be refused."""
    a, b = sympy.symbols('a b')
    ## Second row/column carries nothing: the internal sub-matrix is singular.
    A = sympy.Matrix([[a, 0, b], [0, 0, 0], [b, 0, a]])
    with pytest.raises((ValueError, ZeroDivisionError), match='singular'):
        HierarchicalDDD(A, (1,)).eval({a: 1.0, b: 0.5})


## -- expansion and node ordering -----------------------------------------

def _permute(A, perm):
    """Symmetric permutation: same determinant, different node numbering."""
    return A.extract(list(perm), list(perm))


def test_rcm_is_a_permutation():
    A = bc.rc_ladder(6).A
    perm = reverse_cuthill_mckee(A)
    assert sorted(perm) == list(range(A.rows))


def test_rcm_handles_a_disconnected_pattern():
    """A circuit matrix need not be structurally connected."""
    a, b = sympy.symbols('a b')
    A = sympy.diag(sympy.Matrix([[a, a], [a, a]]), sympy.Matrix([[b, b], [b, b]]))
    assert sorted(reverse_cuthill_mckee(A)) == [0, 1, 2, 3]


def test_permutation_does_not_change_the_determinant():
    """The premise of the whole ordering question."""
    system = bc.rc_ladder(5)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e4
    rng = np.random.default_rng(4)
    ref = complex(ddd_of_matrix(system.A).eval(env))
    for _ in range(4):
        perm = list(rng.permutation(system.A.rows))
        got = complex(ddd_of_matrix(_permute(system.A, perm)).eval(env))
        assert abs(got - ref) <= 1e-9 * abs(ref)


@pytest.mark.parametrize('name', ['ladder', 'mfb'])
def test_auto_ordering_is_insensitive_to_node_numbering(name):
    """The point of the default policy.

    Bare min-degree breaks ties by index, so an arbitrary renumbering of the same
    circuit changes the diagram size several-fold.  Band-reordering first removes
    that dependence.
    """
    A = bc.rc_ladder(12).A if name == 'ladder' else bc.mfb_filter().A
    rng = np.random.default_rng(1)
    perms = [list(rng.permutation(A.rows)) for _ in range(6)]

    auto = [ddd_of_matrix(_permute(A, p), order='auto').size for p in perms]
    plain = [ddd_of_matrix(_permute(A, p), order='min-degree').size for p in perms]

    assert max(auto) / min(auto) <= 1.05          # essentially invariant
    assert max(auto) <= max(plain)                # never worse at the worst case


def test_auto_keeps_a_favourable_given_ordering():
    """It must not throw away a good numbering to gain determinism.

    Nodes added in signal order give a better diagram than band-reordering does
    on the µA741, so ``auto`` tries both and keeps the smaller.
    """
    A = bc.ua741().A
    auto = ddd_of_matrix(A, order='auto').size
    banded = ddd_of_matrix(_permute(A, reverse_cuthill_mckee(A)),
                           order='min-degree').size
    assert auto <= banded


@pytest.mark.parametrize('order', ['auto', 'min-degree', 'markowitz', 'row'])
def test_every_ordering_gives_the_same_value(order):
    system = bc.rc_ladder(5)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e4
    ref = complex(np.linalg.det(np.array(
        [[complex(system.A[i, j].subs(env)) for j in range(system.dim)]
         for i in range(system.dim)])))
    got = complex(ddd_of_matrix(system.A, order=order).eval(env))
    assert abs(got - ref) <= 1e-7 * abs(ref)


## -- recursive (multi-block) hierarchy ------------------------------------

def test_more_smaller_blocks_beat_one_big_one():
    """The point of suppressing in sequence rather than all at once.

    A single split must either leave too many terminals or pay for one large
    cofactor family.  Suppressing a series of small blocks avoids both, because
    each level hands the next something already reduced.
    """
    system = bc.rc_ladder(6)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e4
    ref = complex(ddd_of_matrix(system.A).eval(env))

    one = HierarchicalDDD(system.A, (1, 2))
    many = HierarchicalDDD(system.A, [(1,), (2,), (4,), (5,)])

    for h in (one, many):
        assert abs(complex(h.eval(env)) - ref) <= 1e-8 * abs(ref)
    assert many.size < one.size
    assert len(many.levels) == 4


def test_suppression_order_skips_rows_that_cannot_be_pivots():
    """A voltage source's branch row has no diagonal, so it cannot be suppressed."""
    system = bc.rc_ladder(5)
    chosen = {b[0] for b in suppression_order(system.A)}
    for i in range(system.dim):
        if system.A[i, i] == 0:
            assert i not in chosen


def test_suppression_order_respects_keep():
    system = bc.rc_ladder(5)
    keep = [system.in_index, system.out_index]
    chosen = {b[0] for b in suppression_order(system.A, keep=keep)}
    assert not (chosen & set(keep))


def test_suppression_order_is_ordered_cheapest_first():
    """Low-degree unknowns first: eliminating one couples its neighbours."""
    system = bc.rc_ladder(8)
    blocks = suppression_order(system.A, keep=[system.in_index])
    first = blocks[0][0]

    def adjacency_degree(i):
        ## The implementation counts *neighbours* on the symmetrised pattern,
        ## not row nonzeros -- the diagonal is not a coupling.
        return sum(1 for j in range(system.dim)
                   if j != i and (system.A[i, j] != 0 or system.A[j, i] != 0))

    eligible = [i for i in range(system.dim)
                if system.A[i, i] != 0 and i != system.in_index]
    assert adjacency_degree(first) == min(adjacency_degree(i) for i in eligible)


@pytest.mark.parametrize('limit', [2, 4, 6])
def test_suppression_is_exact_at_every_depth(limit):
    system = bc.rc_ladder(8)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e4
    ref = complex(ddd_of_matrix(system.A).eval(env))

    blocks = suppression_order(system.A, keep=[system.in_index,
                                               system.out_index], limit=limit)
    got = complex(HierarchicalDDD(system.A, blocks).eval(env))
    assert abs(got - ref) <= 1e-9 * abs(ref)


def test_hierarchy_beats_the_flat_diagram_on_the_ua741():
    """H2's gate: beat flat by more than H1's measured ordering spread (1.42×).

    Reaching this needed the numerics fixed as well as the structure. Building
    the stamp as ``D*A_tt - ...`` keeps it polynomial, but carrying that factor
    through and dividing by ``D**(m-1)`` at the end is ruinous once ``m`` is a
    couple of dozen -- the power spans tens of orders of magnitude. Dividing per
    entry instead recovers the Schur complement exactly.
    """
    system = bc.ua741()
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e3
    flat = ddd_of_matrix(system.A)
    ref = complex(flat.eval(env))

    blocks = suppression_order(system.A, keep=[system.in_index,
                                               system.out_index], limit=8)
    h = HierarchicalDDD(system.A, blocks)

    assert abs(complex(h.eval(env)) - ref) <= 1e-8 * abs(ref)
    assert h.size < flat.size / 1.42            # beyond the ordering spread


def test_blocks_must_be_disjoint_and_leave_something():
    A = bc.rc_ladder(4).A
    with pytest.raises(ValueError, match='disjoint'):
        HierarchicalDDD(A, [(1, 2), (2, 3)])
    with pytest.raises(ValueError, match='at least one unknown'):
        HierarchicalDDD(A, [tuple(range(A.rows))])


## -- noise through the reduction ------------------------------------------

def _flat_noise(system, u, CY, env):
    """Reference: substitute everything and solve densely."""
    n = system.dim
    Y = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                  for i in range(n)], dtype=complex)
    un = np.array([complex(sympy.sympify(u[i]).subs(env)) for i in range(n)],
                  dtype=complex)
    z = np.linalg.solve(Y, -un)
    C = np.array([[complex(sympy.sympify(CY[i][j]).subs(env)) for j in range(n)]
                  for i in range(n)], dtype=complex)
    return z, complex(z @ C @ z.conj())


def _noise_setup(system, freq=1e3):
    n = system.dim
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * freq
    u = [0] * n
    u[system.out_index] = 1
    CY = [[(4e-21 if i == j else 0) for j in range(n)] for i in range(n)]
    return env, u, CY


@pytest.mark.parametrize('N', [5, 8])
def test_hierarchical_solve_recovers_every_unknown(N):
    """Including the suppressed ones -- noise needs them."""
    system = bc.rc_ladder(N)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e4
    n = system.dim

    A = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                  for i in range(n)], dtype=complex)
    b = np.array([complex(system.b[i].subs(env)) for i in range(n)],
                 dtype=complex)
    ref = np.linalg.solve(A, b)

    blocks = suppression_order(system.A, keep=[system.in_index,
                                               system.out_index])
    got = hierarchical_solve(HierarchicalDDD(system.A, blocks), system.b, env)
    ## Scale to the largest unknown rather than per element: a genuinely-zero
    ## entry makes a per-element relative error report rounding noise, which
    ## once produced an apparent error of 1e17 on a correct result.
    assert np.max(np.abs(got - ref)) < 1e-9 * np.max(np.abs(ref))


@pytest.mark.parametrize('N', [5, 8])
def test_reduced_noise_matches_the_flat_result(N):
    """Equality, not closeness.

    A reduction that lost a noise source would still return a plausible-looking
    number, so the test has to be exact agreement rather than a tolerance chosen
    to pass.
    """
    system = bc.rc_ladder(N)
    env, u, CY = _noise_setup(system)
    _, flat = _flat_noise(system, u, CY, env)
    _, reduced = noise_psd_reduced(system.A, u, CY, env,
                                   keep=[system.in_index, system.out_index])
    assert abs(reduced - flat) <= 1e-12 * abs(flat)


def test_reduced_noise_keeps_sources_on_suppressed_unknowns():
    """The failure mode this stage exists to avoid.

    Put the *only* noise source on an unknown that gets suppressed. If the
    reduction dropped internal sources the answer would come back near zero --
    quietly, and looking entirely reasonable.
    """
    system = bc.rc_ladder(6)
    n = system.dim
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e3
    u = [0] * n
    u[system.out_index] = 1

    keep = [system.in_index, system.out_index]
    suppressed = suppression_order(system.A, keep=keep)[0][0]
    assert suppressed not in keep

    CY = [[0] * n for _ in range(n)]
    CY[suppressed][suppressed] = 4e-21          # the only source in the circuit

    _, flat = _flat_noise(system, u, CY, env)
    _, reduced = noise_psd_reduced(system.A, u, CY, env, keep=keep)

    assert abs(flat) > 0                        # the source does reach the output
    assert abs(reduced - flat) <= 1e-12 * abs(flat)


def test_reduced_noise_handles_correlated_sources():
    """Sources inside a block are correlated at its terminals.

    A referred correlation is a full matrix, not a diagonal one, so a
    reduction that assumed independence would be wrong here and right on the
    diagonal-only cases above.
    """
    system = bc.rc_ladder(6)
    n = system.dim
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e3
    u = [0] * n
    u[system.out_index] = 1

    rng = np.random.default_rng(11)
    root = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    C = (root @ root.conj().T) * 1e-21          # Hermitian, positive semidefinite
    CY = [[complex(C[i, j]) for j in range(n)] for i in range(n)]

    _, flat = _flat_noise(system, u, CY, env)
    _, reduced = noise_psd_reduced(system.A, u, CY, env,
                                   keep=[system.in_index, system.out_index])
    assert abs(reduced - flat) <= 1e-10 * abs(flat)


def test_reduced_noise_needs_something_to_suppress():
    system = bc.rc_ladder(4)
    with pytest.raises(ValueError, match='nothing can be suppressed'):
        noise_psd_reduced(system.A, [0] * system.dim,
                          [[0] * system.dim] * system.dim, dict(system.params),
                          keep=list(range(system.dim)))


## -- sensitivity ----------------------------------------------------------

def _finite_difference(system, env, parameter, out, rel=1e-6):
    """Central difference on a dense solve -- the independent reference."""
    n = system.dim

    def solve(e):
        A = np.array([[complex(system.A[i, j].subs(e)) for j in range(n)]
                      for i in range(n)], dtype=complex)
        b = np.array([complex(system.b[i].subs(e)) for i in range(n)],
                     dtype=complex)
        return np.linalg.solve(A, b)[out]

    step = abs(complex(env[parameter])) * rel
    plus = dict(env)
    plus[parameter] = env[parameter] + step
    minus = dict(env)
    minus[parameter] = env[parameter] - step
    return (solve(plus) - solve(minus)) / (2 * step)


@pytest.mark.parametrize('N', [3, 4, 5])
def test_determinant_sensitivity_is_exact(N):
    """``d det/dp`` from cofactors must equal sympy's derivative exactly."""
    system = bc.rc_ladder(N)
    for parameter in sorted(system.A.free_symbols - {system.s}, key=str)[:3]:
        combination, _ = determinant_sensitivity(system.A, parameter)
        reference = sympy.diff(system.A.det(), parameter)
        assert sympy.simplify(combination.eval() - reference) == 0


def test_determinant_sensitivity_reuses_one_family():
    """Every parameter should come out of the same shared cofactors.

    That is the reason this is cheap: the derivative of a determinant with
    respect to an entry *is* a cofactor, and the family already holds them.
    """
    system = bc.rc_ladder(6)
    parameters = sorted(system.A.free_symbols - {system.s}, key=str)[:4]

    _, family = determinant_sensitivity(system.A, parameters[0])
    for parameter in parameters[1:]:
        combination, same = determinant_sensitivity(system.A, parameter,
                                                    family=family)
        assert same is family
        assert combination.parts                 # each really uses cofactors


def test_sensitivity_of_an_absent_parameter_is_zero():
    system = bc.rc_ladder(4)
    stranger = sympy.Symbol('not_in_this_circuit')
    combination, _ = determinant_sensitivity(system.A, stranger)
    assert combination.parts == []
    assert combination.eval() == 0


@pytest.mark.parametrize('name', ['ladder', 'mfb'])
def test_adjoint_sensitivities_match_finite_differences(name):
    system = bc.rc_ladder(6) if name == 'ladder' else bc.mfb_filter()
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e3
    parameters = sorted(system.A.free_symbols - {system.s}, key=str)[:5]

    got = adjoint_sensitivities(system.A, system.b, env, parameters,
                                system.out_index)
    references = [_finite_difference(system, env, p, system.out_index)
                  for p in parameters]

    ## Scale to the largest sensitivity: a parameter the output barely depends
    ## on has a reference near zero, and a per-element relative error there
    ## measures rounding noise rather than correctness.
    scale = max(abs(r) for r in references)
    for parameter, reference in zip(parameters, references):
        assert abs(got[parameter] - reference) <= 1e-6 * scale


def test_adjoint_cost_does_not_grow_with_parameter_count():
    """Two solves regardless -- the property that makes this worth doing.

    Differentiating naively costs one solve per parameter; the adjoint form
    costs two in total, so asking about every device is affordable.
    """
    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q3', 'q4', 'q5',
                                        'q6', 'q16', 'q17', 'q23', 'q14'))
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e3
    parameters = sorted(system.A.free_symbols - {system.s}, key=str)

    import time
    timings = []
    for count in (1, len(parameters)):
        start = time.perf_counter()
        adjoint_sensitivities(system.A, system.b, env, parameters[:count],
                              system.out_index)
        timings.append(time.perf_counter() - start)

    ## Ten parameters must not cost ten times one.
    assert timings[1] < 3 * timings[0]


def test_sensitivities_on_the_ua741_match_finite_differences():
    system = bc.ua741(symbolic_devices=('q1', 'q17'))
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e3
    parameters = sorted(system.A.free_symbols - {system.s}, key=str)

    got = adjoint_sensitivities(system.A, system.b, env, parameters,
                                system.out_index)
    references = [_finite_difference(system, env, p, system.out_index)
                  for p in parameters]
    scale = max(abs(r) for r in references)
    for parameter, reference in zip(parameters, references):
        assert abs(got[parameter] - reference) <= 1e-5 * scale


## -- Tier 2 calibration: cascaded opamps and a Cauer filter ---------------

@pytest.mark.parametrize('blocks', [1, 2, 3])
def test_cascaded_opamps_grow_as_expected(blocks):
    system = bc.cascaded_opamps(blocks)
    assert system.cir is not None
    assert system.dim == 2 + 5 * blocks         # source branch, input, 5/stage

    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e3
    n = system.dim
    A = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                  for i in range(n)], dtype=complex)
    b = np.array([complex(system.b[i].subs(env)) for i in range(n)],
                 dtype=complex)
    x = np.linalg.solve(A, b)
    gain = abs(x[system.out_index] / x[system.in_index])
    assert gain > 10 ** blocks                  # each stage contributes gain


def test_cauer_filter_has_inductors_and_stays_polynomial():
    """Inductors get a branch current, not a ``1/(sL)`` admittance.

    That distinction decides whether the s-expansion applies at all: an
    admittance-form inductor would put ``1/s`` in an entry, which has no
    coefficient split.
    """
    system = bc.cauer_lowpass(3)
    expansion = s_expand(system.A, system.s)     # would raise if not polynomial
    assert expansion.degree >= 3
    assert any(str(sym).startswith('Ls') for sym in system.A.free_symbols)


def test_cauer_filter_is_a_lowpass():
    system = bc.cauer_lowpass(3)
    n = system.dim

    def gain(freq):
        env = dict(system.params)
        env[system.s] = 1j * 2 * np.pi * freq
        A = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                      for i in range(n)], dtype=complex)
        b = np.array([complex(system.b[i].subs(env)) for i in range(n)],
                     dtype=complex)
        x = np.linalg.solve(A, b)
        return abs(x[system.out_index] / x[system.in_index])

    assert gain(1e7) < gain(1e3)                # attenuates well above the band


@pytest.mark.parametrize('name,builder', [
    ('cascade', lambda: bc.cascaded_opamps(2)),
    ('cauer', lambda: bc.cauer_lowpass(3)),
])
def test_ddd_is_a_few_times_more_compact_than_soe(name, builder):
    """Calibration against Tan & Shi's own SCAPP comparison.

    Their Table II reports diagram-to-sequence ratios of 5.0–6.0 on cascaded
    opamps.  Ours land in the same band on comparable circuits, which is the
    external check on a claim this project could otherwise only make against
    itself.  Note what it does *not* say: the advantage is a small factor, not
    the orders of magnitude the papers' prose might suggest.
    """
    from pycircuit.circuit.soe import solve_soe

    system = builder()
    out, inp = system.out_index, system.in_index
    _, numerators = ddd_cramer(system.A, system.b, indices=[out, inp])
    diagram = numerators[out].size + numerators[inp].size

    solution = solve_soe(system.A, system.b)
    H = solution.solution[out] / solution.solution[inp]
    sequence = int(sum(sympy.count_ops(e) for _, e in solution.assignments)
                   + sympy.count_ops(H))

    ratio = sequence / diagram
    assert 2.0 < ratio < 10.0, 'ratio %.1f outside the published band' % ratio


## -- conditioning of pole extraction --------------------------------------

def _eigen_poles(system, env):
    """Reference that never forms the polynomial, so it avoids the issue."""
    import scipy.linalg as la
    n = system.dim
    G = np.array([[complex(system.A[i, j].subs({**env, system.s: 0}))
                   for j in range(n)] for i in range(n)], dtype=complex)
    C = np.array([[complex(sympy.diff(system.A[i, j], system.s).subs(env))
                   for j in range(n)] for i in range(n)], dtype=complex)
    return np.sort_complex(np.array(
        [e for e in la.eig(G, -C, right=False) if np.isfinite(e)]))


@pytest.mark.parametrize('N', [12, 16])
def test_extended_precision_recovers_pole_accuracy(N):
    """The error is in the root-finding, not the coefficients.

    Supplying exact coefficients changes nothing; raising the working precision
    changes everything.  A 16-section ladder goes from ~5e-11 to ~2e-15.
    """
    system = bc.rc_ladder(N)
    env = dict(system.params)
    expansion = s_expand(system.A, system.s)
    reference = _eigen_poles(system, env)

    plain = np.sort_complex(expansion.roots_of(env))
    exact = np.sort_complex(expansion.roots_of(env, precision='auto'))

    def error(roots):
        m = min(len(reference), len(roots))
        return np.max(np.abs(roots[:m] - reference[:m])
                      / np.abs(reference[:m]))

    assert error(exact) < error(plain) / 100      # orders of magnitude better
    assert error(exact) < 1e-13


def test_extended_precision_does_not_regress_small_circuits():
    system = bc.rc_ladder(6)
    env = dict(system.params)
    expansion = s_expand(system.A, system.s)
    reference = _eigen_poles(system, env)

    exact = np.sort_complex(expansion.roots_of(env, precision='auto'))
    m = min(len(reference), len(exact))
    assert np.max(np.abs(exact[:m] - reference[:m])
                  / np.abs(reference[:m])) < 1e-13


def test_precision_may_be_given_explicitly():
    system = bc.rc_ladder(8)
    env = dict(system.params)
    expansion = s_expand(system.A, system.s)
    coarse = np.sort_complex(expansion.roots_of(env, precision=25))
    fine = np.sort_complex(expansion.roots_of(env, precision=60))
    assert len(coarse) == len(fine)
    assert np.allclose(coarse, fine, rtol=1e-8)


def test_suggested_precision_grows_with_coefficient_spread():
    small = s_expand(bc.rc_ladder(4).A, bc.rc_ladder(4).s)
    large = s_expand(bc.rc_ladder(16).A, bc.rc_ladder(16).s)
    assert (large._suggested_precision(large.eval_coeffs(dict(bc.rc_ladder(16).params)))
            > small._suggested_precision(small.eval_coeffs(dict(bc.rc_ladder(4).params))))


def test_exact_environment_conversion_is_lossless():
    """A binary float is an exact rational, so nothing is discarded."""
    from pycircuit.circuit.ddd import _exact_env
    a, b = sympy.symbols('a b')
    env = {a: 0.1, b: 1e-9}
    exact = _exact_env(env)
    assert exact[a] == sympy.Rational(0.1)
    assert float(exact[a]) == 0.1
    assert float(exact[b]) == 1e-9


## -- calibration against TCAD 2001 Table II -------------------------------
##
## That table lists thirteen circuits with their complex-DDD and s-expanded
## sizes.  Three are RC ladders, and they turn out to be essentially this
## project's own fixture: matrix size matches exactly and nonzero counts to
## within one, which makes them the closest thing to a like-for-like external
## check anywhere in this work.

TCAD2001_TABLE_II = {
    ## name:      (matrix, nonzeros, |DDD|, s-expanded |DDD|, deg(den))
    'rclad7':     (8, 22, 26, 72, 6),
    'rclad100':   (101, 301, 398, 16767, 85),
    'ua741':      (23, 90, 6654, 99844, 23),
}


@pytest.mark.parametrize('N,row', [(7, 'rclad7')])
def test_ladder_matches_the_published_structure(N, row):
    """Matrix size and sparsity must line up, or the circuits are not the same."""
    matrix, nonzeros, _, _, degree = TCAD2001_TABLE_II[row]
    system = bc.rc_ladder(N)
    count = sum(1 for i in range(system.dim) for j in range(system.dim)
                if system.A[i, j] != 0)

    assert system.dim == matrix
    assert abs(count - nonzeros) <= 2
    assert s_expand(system.A, system.s).degree == degree


@pytest.mark.parametrize('N,row', [(7, 'rclad7')])
def test_ladder_diagram_sizes_are_in_the_published_regime(N, row):
    """Within a small factor of the published sizes.

    Ours runs a little *smaller* on the complex diagram and a little larger on
    the s-expanded one -- consistently, and explicably: the construction here is
    the 2010 expansion ordering rather than the 2000 flow's.  Exact agreement is
    not expected and would be suspicious.
    """
    _, _, published_ddd, published_sexp, _ = TCAD2001_TABLE_II[row]
    system = bc.rc_ladder(N)

    assert ddd_of_matrix(system.A).size == pytest.approx(published_ddd, rel=0.5)
    assert s_expand(system.A, system.s).size == pytest.approx(published_sexp,
                                                              rel=0.5)


def test_ua741_denominator_degree_matches_the_paper():
    """23 poles, which pins the device count and the reactive topology."""
    system = bc.ua741()
    assert s_expand(system.A, system.s).degree == TCAD2001_TABLE_II['ua741'][4]


## -- the review fixes, pinned ---------------------------------------------
##
## Each of these guards a place where a lesson learned elsewhere had not been
## applied.  Without a test they would drift back.

def test_s_expansion_uses_the_ordering_policy_too():
    """H1's ordering fix must reach the s-expansion, not just the determinant.

    The point of ``auto`` is that *every* size is insensitive to how the nodes
    happen to be numbered.  ``s_expand`` defaulted to bare min-degree until the
    review, so its sizes still moved with the numbering.
    """
    A = bc.rc_ladder(10).A
    s = bc.rc_ladder(10).s
    rng = np.random.default_rng(3)
    perms = [list(rng.permutation(A.rows)) for _ in range(5)]

    auto = [s_expand(_permute(A, p), s).size for p in perms]
    plain = [s_expand(_permute(A, p), s, order='min-degree').size for p in perms]

    assert max(auto) / min(auto) <= max(plain) / min(plain)
    assert max(auto) <= max(plain)


def test_cramer_uses_the_ordering_policy_too():
    """Safe there because each numerator is a plain determinant."""
    system = bc.rc_ladder(8)
    rng = np.random.default_rng(5)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e4

    _, reference = ddd_cramer(system.A, system.b, indices=[0])
    for _ in range(3):
        perm = list(rng.permutation(system.A.rows))
        ## Permuting relabels the unknowns, so compare the determinant, which is
        ## what each numerator is.
        assert ddd_of_matrix(_permute(system.A, perm)).size <= \
            ddd_of_matrix(system.A, order='min-degree').size * 4


def test_family_cofactors_stay_addressed_by_original_index():
    """Why `DDDFamily` deliberately keeps ``min-degree``.

    Its results are addressed by *original* row and column, so the band
    reordering ``auto`` performs would silently relabel them.  This pins the
    contract that makes that a real constraint rather than an oversight.
    """
    system = bc.rc_ladder(4)
    family, _ = ddd_cofactor_solve(system.A, system.b, indices=[0])
    for k in range(system.dim):
        for i in range(system.dim):
            minor = system.A.copy()
            minor.row_del(k)
            minor.col_del(i)
            assert sympy.simplify(family.cofactor(k, i).eval() - minor.det()) == 0


def test_dominant_terms_work_on_bare_symbol_payloads():
    """Exercises the ``_resolve`` fast path in the term enumerator.

    ``_entry_values`` had its own substitution until the review, so payloads
    that are bare symbols -- which is what a dense symbolic matrix is made of --
    took the slow path.
    """
    system = bc.dense_symbolic_matrix(5)
    env = dict(system.params)
    D = ddd_of_matrix(system.A)

    total = sum(v[0] for _, v in D.iter_terms(env))
    assert abs(total - complex(D.eval(env))) <= 1e-9 * abs(complex(D.eval(env)))


def test_solution_poles_can_ask_for_extended_precision():
    """H6's remedy must be reachable from the toolkit path, not only from
    ``SExpandedDDD`` directly."""
    from pycircuit.circuit.dddresult import DDDSolution

    system = bc.rc_ladder_semi_symbolic(14, n_symbolic=0)
    solution = DDDSolution(system.A, system.b, system.s,
                           irefnode=system.A.rows)
    plain = np.sort_complex(solution.poles(numeric=True))
    exact = np.sort_complex(solution.poles(numeric=True, precision='auto'))

    assert len(plain) == len(exact)
    assert np.max(np.abs(plain - exact)) < 1e-6 * np.max(np.abs(exact))


def test_noise_psd_is_unchanged_by_the_shared_evaluation():
    """The guard-then-evaluate rewrite must not alter the answer."""
    from pycircuit.circuit.toolkit import ddd_toolkit, symbolic_poly

    s = sympy.Symbol('s', imaginary=True)
    R1, R2, C1, kT = sympy.symbols('R1 R2 C1 kT', positive=True)
    Y = sympy.Matrix([[1 / R1 + 1 / R2 + s * C1, -1 / R2], [-1 / R2, 1 / R2]])
    u = sympy.Matrix([0, 1])
    CY = sympy.Matrix([[4 * kT / R1, 0], [0, 4 * kT / R2]])

    _, reference = symbolic_poly.noise_psd(Y, u, CY, s)
    _, got = ddd_toolkit.noise_psd(Y, u, CY, s)
    assert sympy.simplify(sympy.cancel(got - reference)) == 0


## -- calibration against ICCAD 2010 Table II ------------------------------
##
## Full matrices are the one benchmark that needs no transcription: "the dense
## n x n matrix of distinct symbols" is a complete specification, so the
## extraction rule is satisfied trivially.  This is also the paper the
## construction here is implemented from, which makes it the most direct check
## available -- its LED column is what our builder should reproduce.

ICCAD2010_TABLE_II_LED = {
    2: 4, 3: 12, 4: 32, 5: 80, 6: 192, 7: 448, 8: 1024,
    9: 2304, 10: 5120, 11: 11264,
}

## Their Greedy-Labeling column, for contrast.  We do not implement that
## ordering, so these are not reproduced -- they are here to record how much the
## ordering was worth: 2708 against 1024 at n=8, and 30x by n=18.
ICCAD2010_TABLE_II_GREEDY = {
    2: 4, 3: 13, 4: 40, 5: 118, 6: 340, 7: 965, 8: 2708,
    9: 7535, 10: 20828, 11: 57266,
}


@pytest.mark.parametrize('n', sorted(ICCAD2010_TABLE_II_LED))
def test_full_matrix_sizes_match_the_published_led_column(n):
    """Exact agreement with Shi's own measurements, not just with his theorem.

    The TCAS-II identity says the optimum is ``n*2**(n-1)``; this says the
    ICCAD-2010 construction actually attains it, and that ours attains the same.
    Every value in the paper's LED column is that optimum.
    """
    A = _full_matrix(n)
    assert ddd_of_matrix(A, order='row').size == ICCAD2010_TABLE_II_LED[n]
    assert ICCAD2010_TABLE_II_LED[n] == n * 2 ** (n - 1)


@pytest.mark.parametrize('n', [4, 6, 8])
def test_auto_ordering_also_attains_the_optimum_on_full_matrices(n):
    """Nothing is given up by defaulting to ``auto``.

    On a full matrix every row has the same degree, so min-degree ties resolve
    to the first row and the band reordering has nothing to gain -- both reach
    the optimum.
    """
    assert ddd_of_matrix(_full_matrix(n)).size == ICCAD2010_TABLE_II_LED[n]


@pytest.mark.parametrize('n', [6, 8, 10])
def test_we_are_on_the_good_side_of_the_ordering_comparison(n):
    """Records what the expansion ordering is worth on the paper's own numbers.

    The 2000-era Greedy-Labeling order is what the earlier papers' figures were
    measured with, which is why our µA741 diagram comes out several times
    smaller than the one published in TCAD 2000.
    """
    ours = ddd_of_matrix(_full_matrix(n)).size
    assert ours == ICCAD2010_TABLE_II_LED[n]
    assert ours < ICCAD2010_TABLE_II_GREEDY[n]


## -- cancellation-aware approximation ------------------------------------
##
## The premise, measured in `benchmarks/cancellation_profile.py`: magnitude
## ranking must recover a `1 - tol/cancellation` fraction of the total absolute
## mass, so on the µA741 (cancellation ~1e4) no tractable term count converges.
## Group ranking knows each group's *exact* contribution instead of a bound.


def _ua741_env(gm, freq=1e3):
    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q17'))
    env = {system.s: 2j * np.pi * freq}
    for sym in system.A.free_symbols:
        if sym == system.s:
            continue
        env[sym] = gm if str(sym).startswith('gm_') else complex(
            system.params[sym])
    return system, env


def test_cancellation_is_one_when_nothing_cancels():
    """A matrix whose expansion has no sign changes cannot cancel."""
    a, b, c, d = sympy.symbols('a b c d')
    ## det = a*d - b*c, so making b*c negative aligns both terms in sign.
    D = ddd_of_matrix(sympy.Matrix([[a, b], [-c, d]]))
    env = {a: 2.0, b: 3.0, c: 5.0, d: 7.0}
    assert D.cancellation(env) == pytest.approx(1.0)


def test_cancellation_grows_as_the_sum_approaches_zero():
    a, b, c, d = sympy.symbols('a b c d')
    D = ddd_of_matrix(sympy.Matrix([[a, b], [c, d]]))
    ## a*d - b*c with a*d -> b*c: the absolute mass stays put, the sum vanishes.
    near = D.cancellation({a: 1.0, b: 1.0, c: 1.0, d: 1.0 + 1e-6})
    far = D.cancellation({a: 1.0, b: 1.0, c: 1.0, d: 3.0})
    assert far < 10 < near
    assert near == pytest.approx(2.0 / 1e-6, rel=1e-3)


def test_absolute_companion_is_the_sum_of_term_magnitudes():
    """`subdiagram_values` computes in one pass what enumeration would."""
    system = bc.rc_ladder(5)
    env = _spread_env(system, 3.0)
    env[system.s] = 1j * 2 * np.pi * 1e4
    D = ddd_of_matrix(system.A)
    _values, absolutes = D.subdiagram_values(env)
    by_term = sum(abs(v[0]) for _expr, v in D.iter_terms(env))
    assert absolutes[id(D.root)] == pytest.approx(by_term, rel=1e-9)


def test_group_ranking_converges_where_magnitude_ranking_cannot():
    """The result this whole route exists for, on the calibration amplifier.

    Same diagram, same operating point, same tolerance: magnitude ranking is
    left with an error hundreds of times the answer, group ranking meets the
    tolerance.  The comparison is made at an *equal term count* so that it is
    about the ranking and not about the budget.
    """
    system, env = _ua741_env(1.0e-3)
    D = ddd_of_matrix(system.A)

    _expr, n, err = D.approximate_groups(env, tol=0.05)
    assert err <= 0.05
    assert n < 0.001 * D.term_count()

    with pytest.warns(RuntimeWarning, match='cancellation'):
        _flat, n_flat, err_flat = D.approximate(env, tol=0.05, max_terms=n)
    assert n_flat == n
    assert err_flat > 100 * err            # measured: ~1186% against ~4%


def test_group_ranking_keeps_device_symbols():
    """Numerics rank; they do not replace.  The point of the whole exercise."""
    system, env = _ua741_env(1.0e-3)
    expr, _n, _err = ddd_of_matrix(system.A).approximate_groups(env, tol=0.2)
    names = {str(s) for s in expr.free_symbols}
    assert {'gm_q1', 'gm_q2', 'gm_q17', 's'} <= names


def test_group_ranking_error_is_exact_not_bounded():
    """The reported error must be the kept expression's own error.

    Checked by evaluating the returned sympy expression independently of the
    bookkeeping that selected it -- if the invariant were merely a bound the two
    would differ.
    """
    system = bc.rc_ladder(7)
    env = _spread_env(system, 4.0)
    env[system.s] = 1j * 2 * np.pi * 1e4
    D = ddd_of_matrix(system.A)
    expr, _n, err = D.approximate_groups(env, tol=1e-3)
    exact = complex(D.eval(env))
    got = complex(expr.xreplace(env))
    assert abs(got - exact) / abs(exact) == pytest.approx(err, rel=1e-6)


@pytest.mark.parametrize('tol', [1e-1, 1e-2, 1e-4])
def test_group_ranking_meets_its_tolerance(tol):
    system = bc.rc_ladder(8)
    env = _spread_env(system, 5.0)
    env[system.s] = 1j * 2 * np.pi * 1e4
    _expr, n, err = ddd_of_matrix(system.A).approximate_groups(env, tol=tol)
    assert err <= tol
    assert n >= 1


def test_retaining_minors_shrinks_the_expression():
    """Stage 2: a factored form over small determinants is smaller than the
    fully expanded one, and still names device parameters."""
    system, env = _ua741_env(1.0e-3)
    D = ddd_of_matrix(system.A)
    _expr, n_expanded, err_expanded = D.approximate_groups(env, tol=0.05)
    expr, n_factored, err = D.approximate_groups(env, tol=0.05, max_minor=6)
    assert err <= 0.05 and err_expanded <= 0.05
    assert n_factored < 0.5 * n_expanded
    assert {'gm_q1', 'gm_q2', 'gm_q17'} <= {str(s) for s in expr.free_symbols}


def test_retained_minors_are_determinants_of_named_entries():
    """What distinguishes a retained group from a hierarchical placeholder.

    `minor_positions` recovers which minor each node expands, so an unexpanded
    group is reported as `det(M[rows, cols])` over the original matrix's own
    entries rather than as an opaque symbol.
    """
    system = bc.rc_ladder(6)
    D = ddd_of_matrix(system.A)
    n = D.matrix.rows
    pos = D.minor_positions()
    assert pos[id(D.root)] == (tuple(range(n)), tuple(range(n)))
    ## A 1-edge strikes out the vertex's row and column; a 0-edge does not.
    root = D.root
    assert pos[id(root.one_edge)][0] == tuple(
        r for r in range(n) if r != root.row)
    assert pos[id(root.one_edge)][1] == tuple(
        c for c in range(n) if c != root.col)
    if not root.zero_edge.is_terminal:
        assert pos[id(root.zero_edge)] == pos[id(root)]


def test_group_ranking_needs_numeric_values():
    system = bc.rc_ladder(4)
    with pytest.raises(ValueError, match='symbolic'):
        ddd_of_matrix(system.A).approximate_groups({}, tol=0.1)


def test_unmet_tolerance_is_warned_about_not_only_returned():
    """The mechanism fix.  A previous round reported a 994%-error result as a
    working approximation because the returned error was never printed; an
    unmet tolerance now says so."""
    system, env = _ua741_env(1.0e-3)
    D = ddd_of_matrix(system.A)
    with pytest.warns(RuntimeWarning, match='above the requested tol'):
        D.approximate(env, tol=0.05, max_terms=50)


def test_with_value_returns_what_the_search_already_computed():
    """`with_value` exists because recomputing the value is the expensive part.

    Measured while composing a nested approximation: substituting into the
    returned expression took 33 s against 0.95 s for the ranking, because the
    expression has ~1e5 operations where the diagram has ~1e3 vertices.  The
    value must agree with evaluating the expression, which is what makes it a
    shortcut rather than a different answer.
    """
    system = bc.rc_ladder(7)
    env = _spread_env(system, 4.0)
    env[system.s] = 1j * 2 * np.pi * 1e4
    D = ddd_of_matrix(system.A)

    three = D.approximate_groups(env, tol=1e-3)
    four = D.approximate_groups(env, tol=1e-3, with_value=True)
    assert len(three) == 3 and len(four) == 4
    assert three == four[:3]                       # additive, not a change

    expr, _n, err, value = four
    assert complex(expr.xreplace(env)) == pytest.approx(value, rel=1e-9)
    exact = complex(D.eval(env))
    assert abs(value - exact) / abs(exact) == pytest.approx(err, rel=1e-9)


## -- concentration: the half of the story `cancellation` does not tell --------


def test_concentration_counts_equal_terms():
    """With every term the same magnitude, the effective count is the real one."""
    a, b, c, d = sympy.symbols('a b c d')
    D = ddd_of_matrix(sympy.Matrix([[a, b], [c, d]]))
    ## det = a*d - b*c: two terms, made equal in magnitude.
    n = D.concentration({a: 2.0, b: 2.0, c: 2.0, d: 2.0})
    assert D.term_count() == 2
    assert n == pytest.approx(2.0)


def test_concentration_falls_to_one_when_a_term_dominates():
    a, b, c, d = sympy.symbols('a b c d')
    D = ddd_of_matrix(sympy.Matrix([[a, b], [c, d]]))
    ## a*d enormous against b*c, so one term carries everything.
    n = D.concentration({a: 1e6, b: 1e-6, c: 1e-6, d: 1e6})
    assert n == pytest.approx(1.0, abs=1e-6)


def test_concentration_is_bounded_by_the_term_count():
    system = bc.rc_ladder(6)
    env = _spread_env(system, 3.0)
    env[system.s] = 1j * 2 * np.pi * 1e4
    D = ddd_of_matrix(system.A)
    n = D.concentration(env)
    assert 1.0 <= n <= D.term_count()


def test_concentration_tracks_the_enumerated_mass_profile():
    """Validated against the quantity it stands in for.

    The thing that actually matters is how many terms are needed to reach most of
    the absolute mass.  On a circuit small enough to enumerate, the effective count
    must be the same order of magnitude as that number -- otherwise the
    participation ratio is the wrong proxy.
    """
    system = bc.rc_ladder(6)
    env = _spread_env(system, 3.0)
    env[system.s] = 1j * 2 * np.pi * 1e4
    D = ddd_of_matrix(system.A)

    _values, absolutes = D.subdiagram_values(env)
    target = 0.99 * absolutes[id(D.root)]
    got, needed = 0.0, 0
    for _expr, vals in D.iter_terms(env):
        got += abs(vals[0])
        needed += 1
        if got >= target:
            break

    n_eff = D.concentration(env)
    ## Same order of magnitude, in either direction.
    assert 0.1 * needed <= n_eff <= 10.0 * needed, (n_eff, needed)


def test_concentration_and_cancellation_are_independent():
    """The point of having both: neither implies the other.

    A matrix can be badly conditioned for summation while its magnitude sits in a
    couple of terms, so a diagnostic that reports only `cancellation` cannot say
    whether a ranking will be cheap.
    """
    a, b, c, d = sympy.symbols('a b c d')
    D = ddd_of_matrix(sympy.Matrix([[a, b], [c, d]]))
    ## a*d and b*c nearly cancel, so cancellation is huge -- but there are only
    ## two terms, so concentration is about 2.  Large kappa, tiny N_eff.
    env = {a: 1.0, b: 1.0, c: 1.0, d: 1.0 + 1e-9}
    assert D.cancellation(env) > 1e8
    assert D.concentration(env) == pytest.approx(2.0, rel=1e-6)
