"""Stage 14 -- fully symbolic nonlinear analysis of the leapfrog, then evaluated.

The maintainer's question: nonlinear analysis of the 127-unknown leapfrog, fully
symbolic, with numerical evaluation afterwards.

**Two things to be clear about before any number appears.**

`leapfrog_5th_order` is built from `add_small_signal_bjt`, which stamps only ``R``,
``VCCS`` and ``C`` -- the fixture is **entirely linear and contains no nonlinearity**.
`graded_response_mimo` attaches nonlinearities at chosen ports, so this is not a
blocker, but *where* they attach is a modelling choice.  It is made here explicitly:
on one amplifier's **input differential pair**, because that is where an op-amp's
distortion originates.

"Fully symbolic" means the **expression-graph** sense.  `distortion_ddd.Expr` is a
straight-line program over the device symbols; it evaluates fast and, as the
distortion_ddd documentation already records, cannot rank terms or be read.  That is
exactly the right representation for "symbolic form, numeric evaluation" and exactly
the wrong one for the readability goal the rest of this plan chased -- the two should
not be conflated.

**The oracle is not optional.**  An independent fully-numeric graded solve is run
alongside, and the third harmonic must be non-zero.  The µA741 version of this test
records driving the wrong node and reading a harmonic that was identically zero, which
produced entirely plausible sizes and timings for an analysis that did nothing.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/nonlinear_leapfrog.py
"""

import sys
import time

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.distortion import (GradedSpectrum, GradedVector,
                                          graded_response_mimo)
from pycircuit.circuit.distortion_ddd import Expr, evaluate_one, nodes_of


def cube(spectrum, power):
    return ((spectrum * spectrum).truncated(power) * spectrum).truncated(power)


def expr_solver(A, s_sym, n):
    """Gaussian elimination over `Expr`, straight from the circuit's matrix.

    Sparsity is exploited by skipping zero pivots: a 127x127 MNA matrix is 536
    nonzeros out of 16 129, and eliminating against the zeros would build graph
    nodes carrying no information.  Same shape as the uA741 benchmark's solver,
    which is the one already verified against a numeric oracle.
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
                factor = M[r][i] / M[i][i]
                for j in range(i, n):
                    if not M[i][j]._is_zero():
                        M[r][j] = M[r][j] - factor * M[i][j]
                b[r] = b[r] - factor * b[i]
        x = [None] * n
        for i in range(n - 1, -1, -1):
            acc = b[i]
            for j in range(i + 1, n):
                if not M[i][j]._is_zero():
                    acc = acc - M[i][j] * x[j]
            x[i] = acc / M[i][i]
        return x
    return solve


def numeric_solver(A, s_sym, n, env):
    rows = sympy.Matrix(A)

    def solve(s, rhs):
        local = dict(env)
        local[s_sym] = s
        M = np.array([[complex(rows[i, j].subs(local)) for j in range(n)]
                      for i in range(n)], dtype=complex)
        return list(np.linalg.solve(M, np.asarray(rhs, dtype=complex)))
    return solve


def run(system, nonlinear_nodes, power, drive_index, out_index, harmonic=3):
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
                out[node] = cube(x[node], power).scaled(coeffs[idx])
            return GradedVector(out)
        return f

    def source(value):
        vec = GradedVector([GradedSpectrum() for _ in range(n)])
        vec.components[drive_index].terms.update(
            GradedSpectrum.from_phasor(1, 1, value).terms)
        return vec

    tone = 2 * np.pi * 1e3

    t0 = time.time()
    graph = graded_response_mimo(expr_solver(system.A, s_sym, n),
                                 source(Expr.atom(drive)),
                                 make_f([Expr.atom(k) for k in ks]),
                                 (tone,), max_power=power)
    root = graph[out_index].phasor(harmonic)
    build = time.time() - t0

    t0 = time.time()
    got = evaluate_one(root, env)
    ev = time.time() - t0

    t0 = time.time()
    numeric = graded_response_mimo(numeric_solver(system.A, s_sym, n, env),
                                   source(complex(env[drive])),
                                   make_f([complex(env[k]) for k in ks]),
                                   (tone,), max_power=power)
    want = numeric[out_index].phasor(harmonic)
    num = time.time() - t0

    return {'nodes': nodes_of(root), 'build': build, 'eval': ev,
            'numeric': num, 'got': got, 'want': want}


def main():
    t0 = time.time()
    sym = tuple('s%d_q%s' % (k, q) for k in range(5) for q in ('1', '2', '17'))
    system = bc.leapfrog_5th_order(symbolic_devices=sym)
    print('leapfrog: %d unknowns, %d symbols, built in %.1f s'
          % (system.dim, len(system.A.free_symbols), time.time() - t0))

    ## The modelling choice, stated: the nonlinearity sits on stage 0's input
    ## differential pair, where an op-amp's distortion originates.
    nodes = []
    for name in ('s0_e1', 's0_e2'):
        nodes.append(system.cir.get_node_index(name))
    nodes = [i for i in nodes if i < system.dim]
    drive_index = [i for i in range(system.dim) if system.b[i] != 0][0]
    print('nonlinearity on stage 0 input pair, rows %s; driven at row %d,'
          ' read at row %d' % (nodes, drive_index, system.out_index))
    print()

    for power in (3, 5):
        print('== U^%d ==' % power)
        try:
            r = run(system, nodes, power, drive_index, system.out_index)
        except Exception as exc:
            print('  FAILED: %s: %s' % (type(exc).__name__, exc))
            break
        rel = (abs(r['got'] - r['want']) / abs(r['want'])
               if abs(r['want']) > 0 else float('inf'))
        print('  graph %d nodes, symbolic build %.1f s, evaluate %.2f s,'
              ' numeric oracle %.1f s'
              % (r['nodes'], r['build'], r['eval'], r['numeric']))
        print('  third harmonic: symbolic %.6e, numeric %.6e'
              % (abs(r['got']), abs(r['want'])))
        print('  GATE 14-3 (harmonic is non-zero): %s'
              % ('PASS' if abs(r['want']) > 0 else 'FAIL -- measures nothing'))
        print('  GATE 14-2 (matches the oracle to 1e-10): %s  (rel %.2e)'
              % ('PASS' if rel <= 1e-10 else 'FAIL', rel))
    return 0


if __name__ == '__main__':
    sys.exit(main())
