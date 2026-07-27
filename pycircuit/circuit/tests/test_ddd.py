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
                                   suppression_order, s_expand)


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
    lo = complex(D.eval(dict(env, **{}) | {system.s: 1j * 2 * np.pi * 1e2}))
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
