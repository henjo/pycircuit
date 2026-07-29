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


## ------------------------------------------------ stage C2: real circuits --

def _bench_solver(A, s_sym, n):
    """Gaussian elimination over :class:`Expr`, from a circuit's own matrix.

    Sparsity is exploited by skipping zero pivots; a real MNA matrix is
    mostly zeros (the µA741 is 103 nonzeros out of 676) and eliminating
    against them would build nodes carrying no information.
    """
    rows = [[A[i, j] for j in range(n)] for i in range(n)]

    def solve(s, rhs):
        M = [[Expr.atom(e.subs(s_sym, s)
                        if getattr(e, 'free_symbols', set()) else e)
              if e != 0 else Expr.atom(0) for e in row] for row in rows]
        b = [r if isinstance(r, Expr) else Expr.atom(r) for r in rhs]
        for i in range(n):
            for r in range(i + 1, n):
                if M[r][i]._is_zero():
                    continue
                factor = M[r][i]/M[i][i]
                for j in range(i, n):
                    if not M[i][j]._is_zero():
                        M[r][j] = M[r][j] - factor*M[i][j]
                b[r] = b[r] - factor*b[i]
        x = [None]*n
        for i in range(n - 1, -1, -1):
            acc = b[i]
            for j in range(i + 1, n):
                if not M[i][j]._is_zero():
                    acc = acc - M[i][j]*x[j]
            x[i] = acc/M[i][i]
        return x
    return solve


def _numeric_bench_solver(A, s_sym, n, env):
    rows = sympy.Matrix(A)

    def solve(s, rhs):
        local = dict(env)
        local[s_sym] = s
        M = np.array([[complex(rows[i, j].subs(local)) for j in range(n)]
                      for i in range(n)], dtype=complex)
        return list(np.linalg.solve(M, np.asarray(rhs, dtype=complex)))
    return solve


def _bench_case(system, nonlinear_nodes, power, drive_index, out_index):
    n, s_sym = system.dim, system.s
    ks = [sympy.Symbol('kk%d' % i) for i in nonlinear_nodes]
    drive = sympy.Symbol('Adrv')

    env = {k: 1e-9 for k in ks}
    env[drive] = 1e-3
    for symbol in system.A.free_symbols:
        if str(symbol).startswith('gm_'):
            env[symbol] = 1e-3

    def make_f(coeffs):
        def f(x):
            out = [GradedSpectrum() for _ in range(n)]
            for idx, node in enumerate(nonlinear_nodes):
                out[node] = _cube(x[node], power).scaled(coeffs[idx])
            return GradedVector(out)
        return f

    def source(value):
        vec = GradedVector([GradedSpectrum() for _ in range(n)])
        vec.components[drive_index].terms.update(
            GradedSpectrum.from_phasor(1, 1, value).terms)
        return vec

    tone = 2*np.pi*1e3
    graph = graded_response_mimo(
        _bench_solver(system.A, s_sym, n), source(Expr.atom(drive)),
        make_f([Expr.atom(k) for k in ks]), (tone,), max_power=power)
    root = graph[out_index].phasor(3)

    numeric = graded_response_mimo(
        _numeric_bench_solver(system.A, s_sym, n, env),
        source(complex(env[drive])),
        make_f([complex(env[k]) for k in ks]), (tone,), max_power=power)

    return (evaluate_one(root, env), numeric[out_index].phasor(3),
            nodes_of(root))


def test_stage_c2_ua741_matches_the_numeric_path():
    """Stage C2's gate, on a real circuit rather than a synthetic matrix.

    The µA741 as transcribed for the DDD work: 26x26, 103 nonzeros, built
    through the same path AC analysis uses.  The gate asked for 1e-10.

    Driven at the circuit's *own* input and read at its *own* output.  An
    earlier run drove node 0 and read node 0 instead, which produced entirely
    plausible build times and node counts against a harmonic that was
    identically zero -- a reminder that a size measurement says nothing about
    whether the analysis did anything.
    """
    from pycircuit.circuit import benchmark_circuits as bc

    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q3', 'q4'))
    drive_index = [i for i in range(system.dim) if system.b[i] != 0][0]
    got, want, size = _bench_case(system, [17, 18, 22], 5,
                                  drive_index, system.out_index)

    assert abs(want) > 0, 'the oracle must produce a nonzero third harmonic'
    assert abs(got - want) <= 1e-10*abs(want)
    ## A few thousand nodes for a 26-node op-amp; guards a blow-up.
    assert size < 20000, size


def test_stage_c2_size_grows_linearly_with_circuit_size():
    """The axis that decides whether bigger circuits are reachable.

    Cascaded op-amps, with the nonlinearity in the **first** block and the
    output read at the last, so the representation has to carry the whole
    chain.  Placing it in the last block instead makes size *constant* in the
    number of blocks, which flatters the result and measures locality rather
    than scale -- so it is deliberately not done that way here.
    """
    from pycircuit.circuit import benchmark_circuits as bc

    sizes = []
    for blocks in (1, 3):
        system = bc.cascaded_opamps(blocks=blocks,
                                    symbolic_devices=('q1', 'q2'))
        drive_index = [i for i in range(system.dim) if system.b[i] != 0][0]
        got, want, size = _bench_case(system, [1, 2, 3], 5,
                                      drive_index, system.out_index)
        assert abs(want) > 0
        assert abs(got - want) <= 1e-10*abs(want)
        sizes.append(size)

    ## Linear, not exponential: tripling the blocks must not square the graph.
    assert sizes[1] < 4*sizes[0], sizes
    assert sizes[1] > sizes[0], sizes
