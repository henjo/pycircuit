# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Tests for the DDD-backed AC solution and its toolkit integration.

The important property here is not that the diagram path produces the same
answers as ``symbolic_poly`` -- though it must, and does below -- but that it
keeps producing them at sizes where the expanded form has ceased to be a usable
object at all.  ``test_graph_path_works_where_expansion_cannot`` is the test that
states that, and it is the point of the whole exercise.
"""

import numpy as np
import pytest
import sympy

from pycircuit.circuit import SubCircuit, R, C, VS, gnd
from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.acsolution import ACSolution, NumDenSolution
from pycircuit.circuit.analysis_ss import AC
from pycircuit.circuit.ddd import DDDSizeError
from pycircuit.circuit.dddresult import DDDSolution
from pycircuit.circuit.toolkit import ddd_toolkit, symbolic, symbolic_poly


def _rc_circuit():
    """The textbook RC low-pass, fully symbolic."""
    Rv, Cv = sympy.symbols('R C', positive=True)
    cir = SubCircuit(toolkit=symbolic)
    nin, nout = cir.add_node('in'), cir.add_node('out')
    cir['vs'] = VS(nin, gnd, vac=1)
    cir['R'] = R(nin, nout, r=Rv)
    cir['C'] = C(nout, gnd, c=Cv)
    return cir, Rv, Cv


## -- the ACSolution contract ---------------------------------------------

def _solutions_for(system):
    """One solution of each flavour over the same system, for conformance."""
    num, den = symbolic_poly.linearsolver_num_den(system.A, system.b)
    num = np.concatenate((num, np.array([0], dtype=object)))
    return {
        'numden': NumDenSolution(num, den, system.s),
        'ddd': DDDSolution(system.A, system.b, system.s,
                           irefnode=system.A.rows),
    }


@pytest.mark.parametrize('flavour', ['numden', 'ddd'])
def test_solution_conformance(flavour):
    """Both representations must satisfy the same interface.

    A single parametrised conformance test is what keeps two implementations of
    a duck-typed contract from drifting apart -- cheap, and it fails loudly the
    moment one of them grows or loses a method.
    """
    system = bc.rc_ladder(3)
    solution = _solutions_for(system)[flavour]

    assert isinstance(solution, ACSolution)
    assert solution.s is system.s

    nums = solution.numerators()
    den = solution.denominator()
    assert len(nums) == system.dim + 1          # reference node re-inserted
    assert den != 0

    x = solution.node_voltages()
    assert len(x) == len(nums)

    tf = solution.transfer_function(nums[system.out_index])
    assert tf.poles() == solution.poles()


@pytest.mark.parametrize('flavour', ['numden', 'ddd'])
def test_solution_numerators_support_the_arithmetic_tf_i_needs(flavour):
    """``tf_i`` does ``0*N``, ``s*N`` and ``A + c0*(D-1)`` on the numerators.

    That contract is duck-typed rather than declared, so it is asserted here
    directly.
    """
    system = bc.rc_ladder(3)
    solution = _solutions_for(system)[flavour]
    N, D, s = solution.numerators(), solution.denominator(), solution.s
    assert len(0 * N) == len(N)
    assert len(s * N) == len(N)
    assert len(N + (0 * N) * (D - 1)) == len(N)


## -- agreement with the established backend ------------------------------

def test_ac_through_the_toolkit_matches_symbolic_poly():
    cir, Rv, Cv = _rc_circuit()
    s = sympy.Symbol('s')

    ref = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)
    got = AC(cir, toolkit=ddd_toolkit).solve(s, complexfreq=True)

    assert sympy.simplify(got.tf('out', gnd).expr()
                          - ref.tf('out', gnd).expr()) == 0
    assert sympy.simplify(got.tf_i('R.plus').expr()
                          - ref.tf_i('R.plus').expr()) == 0
    assert got.poles() == ref.poles()


def test_exact_values_stay_exact():
    """No stray floats: an exact circuit must give an exact transfer function."""
    cir, Rv, Cv = _rc_circuit()
    s = sympy.Symbol('s')
    H = sympy.simplify(AC(cir, toolkit=ddd_toolkit)
                       .solve(s, complexfreq=True).tf('out', gnd).expr())
    assert H == 1 / (Rv * Cv * s + 1)
    assert not H.atoms(sympy.Float)


def test_toolkit_advertises_ddd():
    assert ddd_toolkit.supports('ddd') is True
    assert ddd_toolkit.supports('num_den') is True
    assert symbolic_poly.supports('ddd') is False


@pytest.mark.parametrize('N', [3, 5])
def test_solution_numerators_match_symbolic_poly(N):
    system = bc.rc_ladder(N)
    ref_num, ref_den = symbolic_poly.linearsolver_num_den(system.A, system.b)
    sol = DDDSolution(system.A, system.b, system.s, irefnode=system.A.rows)

    ## Both are Cramer solutions of the same system, so the ratios agree even if
    ## an overall sign convention does not.
    got = sol.numerators()
    for i in (0, system.out_index):
        assert sympy.simplify(got[i] / sol.denominator()
                              - ref_num[i] / ref_den) == 0


## -- the point of the exercise -------------------------------------------

def test_graph_path_works_where_expansion_cannot():
    """A 32-section ladder: 94 vertices, 2.18 million product terms.

    The expanded determinant is not a useful object at this size, and the
    compatibility path says so rather than trying.  The diagram path answers
    anyway -- which is the entire justification for the representation.
    """
    system = bc.rc_ladder_semi_symbolic(32, n_symbolic=0)   # numeric components
    sol = DDDSolution(system.A, system.b, system.s, irefnode=system.A.rows)

    ## Compact: a few dozen vertices standing for millions of terms.
    diagram = sol.denominator_diagram()
    assert diagram.size < 200
    assert diagram.term_count() > 1_000_000

    ## Expanding is refused, with an error that says what to do instead.
    with pytest.raises(DDDSizeError, match='diagram API'):
        sol.denominator()

    ## The graph path is unbothered.
    poles = sol.poles(numeric=True)
    ## 31, not 32: the voltage source sits directly across the first capacitor,
    ## so that one is not an independent state.
    assert len(poles) == 31
    assert np.all(poles.real < 0)                          # a passive RC ladder


def test_numeric_poles_agree_with_a_reference_solve():
    """Poles from the diagram vs a generalized eigenvalue solve."""
    import scipy.linalg as la

    system = bc.rc_ladder_semi_symbolic(12, n_symbolic=0)
    sol = DDDSolution(system.A, system.b, system.s, irefnode=system.A.rows)
    got = np.sort_complex(sol.poles(numeric=True))

    n = system.dim
    G = np.array([[complex(system.A[i, j].subs({system.s: 0}))
                   for j in range(n)] for i in range(n)])
    Cm = np.array([[complex(sympy.diff(system.A[i, j], system.s))
                    for j in range(n)] for i in range(n)])
    ref = np.sort_complex(np.array(
        [e for e in la.eig(G, -Cm, right=False) if np.isfinite(e)]))

    assert len(got) == len(ref)
    assert np.max(np.abs(got - ref) / np.abs(ref)) < 1e-9


@pytest.mark.parametrize('N', [4, 8, 16])
def test_eval_tf_matches_a_numeric_solve(N):
    """Numeric transfer function straight off the diagrams."""
    system = bc.rc_ladder(N)
    sol = DDDSolution(system.A, system.b, system.s, irefnode=system.A.rows)
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * 1e6

    got = sol.eval_tf(system.out_index, system.in_index, env)

    n = system.dim
    An = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                   for i in range(n)])
    bn = np.array([complex(system.b[i].subs(env)) for i in range(n)])
    x = np.linalg.solve(An, bn)
    ref = x[system.out_index] / x[system.in_index]
    assert abs(got - ref) <= 1e-9 * max(1.0, abs(ref))


def test_reference_node_has_zero_numerator():
    system = bc.rc_ladder(3)
    sol = DDDSolution(system.A, system.b, system.s, irefnode=system.A.rows)
    assert sol.numerator_diagram(system.A.rows) is None
    assert sol.numerators()[system.A.rows] == 0
    assert sol.eval_node(system.A.rows, dict(system.params)) == 0


def test_symbolic_poles_are_guarded():
    """Symbolic poles must expand, so they are limited -- and say so."""
    system = bc.rc_ladder(24)
    sol = DDDSolution(system.A, system.b, system.s, irefnode=system.A.rows)
    with pytest.raises(DDDSizeError, match='numeric=True'):
        sol.poles(numeric=False)
