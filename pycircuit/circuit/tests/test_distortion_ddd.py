# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Stage A of ``doc/distortion_ddd_plan.md``.

The question is whether the cost recorded in ``distortion_mimo_plan.md``
section 8.3 -- ``sympy.cancel`` not finishing in 900 s at ``U^7`` -- is
intrinsic to the method or an artifact of asking sympy to *simplify*.

Every gate here compares against the **existing** implementation running the
identical recurrence in plain complex arithmetic.  That is the point of the
separate module: `distortion.py` is untouched, so the comparison is between
two representations of the same algebra rather than between a thing and
itself.
"""

import time

import numpy as np
import pytest
import sympy

from pycircuit.circuit.distortion import (graded_response_mimo, GradedSpectrum,
                                          GradedVector)
from pycircuit.circuit.distortion_ddd import (Expr, evaluate, evaluate_one,
                                              nodes_of)


def _cube(spectrum, power):
    return ((spectrum*spectrum).truncated(power) * spectrum).truncated(power)


## ------------------------------------------------------------- algebra --

def test_expressions_are_interned_and_commutative_forms_coincide():
    """Structural equality must be identity, or the memo cannot key on it.

    Sorting arguments is what makes ``a*b`` and ``b*a`` one node.  sympy
    hashes structurally too, but this is the property the evaluation memo
    depends on, so it is pinned rather than assumed.
    """
    a, b = Expr.atom(sympy.Symbol('a')), Expr.atom(sympy.Symbol('b'))
    assert a*b is b*a
    assert a + b is b + a
    assert (a*b)*a is a*(b*a)

    ## Identities must collapse, or every pass through the recurrence would
    ## add a layer of `*1` and `+0` wrappers.
    assert a * Expr.atom(1) is a
    assert a + Expr.atom(0) is a
    assert (a * Expr.atom(0))._is_zero()

    ## GradedSpectrum drops zero coefficients, so `!= 0` has to work.
    assert a != 0
    assert not (Expr.atom(0) != 0)


def test_evaluation_is_memoised_across_shared_structure():
    """A shared subgraph must be evaluated once, not once per reference."""
    x = Expr.atom(sympy.Symbol('x'))
    shared = (x + Expr.atom(1)) * (x + Expr.atom(2))
    roots = [shared * Expr.atom(k) for k in range(1, 6)]

    ## Five roots over one shared subgraph: the graph is far smaller than
    ## the five trees would be.
    assert nodes_of(*roots) < 20

    values = evaluate(roots, {sympy.Symbol('x'): 3.0})
    for k, root in enumerate(roots, start=1):
        assert abs(values[id(root)] - k*(4.0*5.0)) < 1e-12


## ------------------------------------------------- the stage A sweep --

def _two_node_symbols():
    names = 'a c1 c2 b d k3 A w0'
    return dict(zip(names.split(), sympy.symbols(names)))


def _run(power, solver, coefficient, drive, tone):
    source = GradedVector([GradedSpectrum.from_phasor(1, 1, drive),
                           GradedSpectrum()])

    def f(x):
        return GradedVector([GradedSpectrum(),
                             _cube(x[0], power).scaled(coefficient)])

    return graded_response_mimo(solver, source, f, (tone,), max_power=power)


def _expr_solver(S):
    def solve(s, rhs):
        sv = Expr.atom(s)
        det = (S['a'] + sv*S['c1'])*(sv*S['c2']) - S['b']*S['d']
        return [(sv*S['c2']*rhs[0] + S['b']*rhs[1])/det,
                (S['d']*rhs[0] + (S['a'] + sv*S['c1'])*rhs[1])/det]
    return solve


def _numeric_solver(env, sym):
    def solve(s, rhs):
        det = ((env[sym['a']] + s*env[sym['c1']])*(s*env[sym['c2']])
               - env[sym['b']]*env[sym['d']])
        return [(s*env[sym['c2']]*rhs[0] + env[sym['b']]*rhs[1])/det,
                (env[sym['d']]*rhs[0]
                 + (env[sym['a']] + s*env[sym['c1']])*rhs[1])/det]
    return solve


@pytest.mark.parametrize('power', [7, 9, 11, 13])
def test_stage_a_sweep_matches_the_numeric_path(power):
    """The stage A gate: ``U^7`` through ``U^13``, agreeing to 1e-12.

    ``sympy.cancel`` did not finish in 900 s at ``U^7`` on this same system.
    The representation here is a hash-consed graph that is never expanded and
    never simplified -- cost moves to evaluation, which is memoised.
    """
    sym = _two_node_symbols()
    atoms = {name: Expr.atom(s) for name, s in sym.items()}
    env = {sym['a']: 1e-3, sym['c1']: 1e-12, sym['c2']: 2e-12,
           sym['b']: 3e-4, sym['d']: 4e-4, sym['k3']: 5e-3,
           sym['A']: 1e-3, sym['w0']: 2*np.pi*1e6}

    symbolic = _run(power, _expr_solver(atoms), atoms['k3'], atoms['A'],
                    sym['w0'])
    got = evaluate_one(symbolic[0].phasor(3), env)

    numeric = _run(power, _numeric_solver(env, sym), complex(env[sym['k3']]),
                   complex(env[sym['A']]), complex(env[sym['w0']]))
    want = numeric[0].phasor(3)

    assert abs(got - want) <= 1e-12*abs(want)


def test_stage_a_cost_is_polynomial_in_truncation_order():
    """The sweep's real deliverable: *cost* scaling, not a single pass mark.

    Stage E of the previous plan established that expression **size** grows
    polynomially and said nothing about cost.  This asserts the graph grows
    sub-exponentially in the truncation order -- the property that makes
    ``U^13`` reachable at all.
    """
    sym = _two_node_symbols()
    atoms = {name: Expr.atom(s) for name, s in sym.items()}

    sizes = []
    for power in (7, 9, 11, 13):
        solution = _run(power, _expr_solver(atoms), atoms['k3'], atoms['A'],
                        sym['w0'])
        sizes.append(nodes_of(solution[0].phasor(3)))

    assert sizes == sorted(sizes), sizes
    ## Geometric growth would roughly square over three steps; polynomial
    ## growth does not.  Measured 89 -> 434, a factor of ~4.9.
    assert sizes[-1] < 8*sizes[0], sizes


def test_stage_a_cost_is_polynomial_in_circuit_size():
    """The axis that matters for bigger circuits, on the *unfavourable* case.

    A dense matrix, so Gaussian elimination fills in completely -- a ladder
    would be the easy case and would flatter the result.  Growth measured at
    about ``n**2``; this guards against it turning exponential.
    """
    power = 7
    w0, drive = sympy.Symbol('w0'), sympy.Symbol('A')

    def build(n):
        entries = [[sympy.Symbol('y%d_%d' % (i, j)) for j in range(n)]
                   for i in range(n)]
        caps = [sympy.Symbol('cc%d' % i) for i in range(n)]
        coeffs = [sympy.Symbol('kk%d' % i) for i in range(n)]
        Ee = [[Expr.atom(e) for e in row] for row in entries]
        Ec = [Expr.atom(c) for c in caps]
        Ek = [Expr.atom(k) for k in coeffs]

        def solve(s, rhs):
            sv = Expr.atom(s)
            M = [[Ee[i][j] + (sv*Ec[i] if i == j else Expr.atom(0))
                  for j in range(n)] for i in range(n)]
            b = [r if isinstance(r, Expr) else Expr.atom(r) for r in rhs]
            for i in range(n):
                for r in range(i+1, n):
                    factor = M[r][i]/M[i][i]
                    for j in range(i, n):
                        M[r][j] = M[r][j] - factor*M[i][j]
                    b[r] = b[r] - factor*b[i]
            x = [None]*n
            for i in range(n-1, -1, -1):
                acc = b[i]
                for j in range(i+1, n):
                    acc = acc - M[i][j]*x[j]
                x[i] = acc/M[i][i]
            return x

        def f(x):
            return GradedVector([_cube(x[i], power).scaled(Ek[i])
                                 for i in range(n)])

        source = GradedVector(
            [GradedSpectrum.from_phasor(1, 1, Expr.atom(drive))]
            + [GradedSpectrum() for _ in range(n-1)])
        solution = graded_response_mimo(solve, source, f, (w0,),
                                        max_power=power)
        return nodes_of(solution[0].phasor(3))

    small, large = build(3), build(6)
    ## Doubling the node count must not do worse than cubic.
    assert large < 8*small, (small, large)
    assert large > small, (small, large)
