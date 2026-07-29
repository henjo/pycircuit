"""Stage 0 of `doc/cancellation_ranking_plan.md` -- where the cancellation lives.

Term-ranked approximation ranks product terms by |product|.  On the uA741 that
leaves 994% error at 500 of 2.77M terms.  This script measures *why*, and
whether the second route of `doc/hierarchical_approximation_plan.md` section 5
-- ranking groups of terms that cancel together -- has a premise to stand on.

For every vertex of the diagram it computes

    V[v] = exact sum of the product terms below v      (signs kept)
    A[v] = sum of their absolute values                (the absolute companion)
    k[v] = A[v] / |V[v]|                               (condition number of the sum)

Both passes cost one traversal of the DAG, because absolute value distributes
over a product.  ``k[root]`` is the number that decides whether magnitude
ranking can converge at all: to reach relative error ``tol`` it must capture a
``1 - tol/k`` fraction of the total absolute mass.

Run:  PYTHONPATH=<repo> python3 benchmarks/cancellation_profile.py
"""

import math
import sys
import time

import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix, _resolve


FREQ = 1.0e3
## The previous round's operating points, kept identical across every stage of
## this plan so that its recorded magnitude-ranking failure is reproduced rather
## than quoted, and so that every kappa printed here belongs to the same circuit
## state as the term counts in `cancellation_groups.py`.
NOMINAL = {'gm_q1': 1.0e-3, 'gm_q2': 1.0e-3, 'gm_q17': 1.0e-3}
DEGRADED = {'gm_q1': 0.1e-3, 'gm_q2': 0.1e-3, 'gm_q17': 0.1e-3}


def environment(system, gms):
    """Numeric substitution: every symbol of the matrix, plus s at FREQ."""
    env = {system.s: 2j * math.pi * FREQ}
    for sym in sorted(system.A.free_symbols, key=str):
        if sym == system.s:
            continue
        name = str(sym)
        if name in gms:
            env[sym] = gms[name]
        else:
            env[sym] = complex(system.params[sym])
    return env


def profile(ddd, env):
    """Per-vertex V, A and depth, in one traversal each.

    Returns ``(values, absolutes, depths)`` keyed by ``id(node)``.  Depth is the
    shortest distance from the root, which is what "a cut of the diagram" in the
    plan's gate refers to.
    """
    entry = {}
    for v in ddd.vertices():
        key = id(v.entry)
        if key not in entry:
            entry[key] = complex(_resolve(v.entry, env))

    values, absolutes = {}, {}
    stack = [(ddd.root, False)]
    while stack:
        node, expanded = stack.pop()
        if node.is_terminal:
            val = complex(_resolve(node.value, env))
            values.setdefault(id(node), val)
            absolutes.setdefault(id(node), abs(val))
            continue
        if id(node) in values:
            continue
        if not expanded:
            stack.append((node, True))
            stack.append((node.one_edge, False))
            stack.append((node.zero_edge, False))
            continue
        e = entry[id(node.entry)]
        values[id(node)] = (node.sign * e * values[id(node.one_edge)]
                            + values[id(node.zero_edge)])
        absolutes[id(node)] = (abs(e) * absolutes[id(node.one_edge)]
                               + absolutes[id(node.zero_edge)])

    depths = {id(ddd.root): 0}
    frontier = [ddd.root]
    while frontier:
        nxt = []
        for node in frontier:
            if node.is_terminal:
                continue
            for child in (node.one_edge, node.zero_edge):
                d = depths[id(node)] + 1
                if id(child) not in depths or d < depths[id(child)]:
                    depths[id(child)] = d
                    nxt.append(child)
        frontier = nxt
    return values, absolutes, depths


def kappa_of(ddd, values, absolutes):
    """k[v] for every *non-terminal* vertex; terminals are single terms."""
    out = {}
    for v in ddd.vertices():
        val, av = values[id(v)], absolutes[id(v)]
        if av == 0.0:
            continue
        out[id(v)] = math.inf if val == 0 else av / abs(val)
    return out


def report(name, ddd, env):
    t0 = time.time()
    values, absolutes, depths = profile(ddd, env)
    kap = kappa_of(ddd, values, absolutes)
    secs = time.time() - t0

    root = ddd.root
    k_root = kap[id(root)]
    print('--- %s ---' % name)
    print('  V[root]           = %.6e%+.6ej' % (values[id(root)].real,
                                                values[id(root)].imag))
    print('  A[root]           = %.6e' % absolutes[id(root)])
    print('  k[root]           = %.3e' % k_root)
    print('  double-prec error = %.2e   (k * 2**-53)' % (k_root * 2.0 ** -53))
    for tol in (0.05, 1e-3):
        frac = 1.0 - tol / k_root
        print('  magnitude ranking at tol=%-6g must capture %.8f%% of A'
              % (tol, 100.0 * frac))
    print('  profiled %d vertices in %.3f s' % (len(kap), secs))

    ks = sorted(kap.values())
    n = len(ks)
    def q(p):
        return ks[min(n - 1, int(p * n))]
    print('  k over all vertices: median %.3g  p90 %.3g  p99 %.3g  max %.3g'
          % (q(0.5), q(0.9), q(0.99), ks[-1]))
    below = sum(1 for k in ks if k < 10.0)
    print('  vertices with k < 10: %d / %d  (%.1f%%)'
          % (below, n, 100.0 * below / n))

    print('  by depth from root:')
    print('    %-6s %-7s %-11s %-11s %-11s' % ('depth', 'count', 'median k',
                                               'max k', 'frac k<10'))
    by_depth = {}
    for v in ddd.vertices():
        if id(v) in kap:
            by_depth.setdefault(depths[id(v)], []).append(kap[id(v)])
    for d in sorted(by_depth):
        row = sorted(by_depth[d])
        lo = sum(1 for k in row if k < 10.0)
        print('    %-6d %-7d %-11.3g %-11.3g %.2f'
              % (d, len(row), row[len(row) // 2], row[-1], lo / len(row)))
    return k_root


def main():
    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q17'))
    print('uA741: dim %d, symbols %d' % (system.dim, len(system.A.free_symbols)))
    t0 = time.time()
    ddd = ddd_of_matrix(system.A)
    print('flat diagram: %d vertices, %d terms, built in %.2f s'
          % (ddd.size, ddd.term_count(), time.time() - t0))
    print('frequency: %g Hz' % FREQ)
    print()

    worst = 0.0
    for name, gms in (('nominal gm', NOMINAL), ('degraded gm', DEGRADED)):
        worst = max(worst, report(name, ddd, environment(system, gms)))
        print()

    print('PRECISION CHECK: k[root]*2**-53 = %.2e -- mpmath %s'
          % (worst * 2.0 ** -53,
             'REQUIRED (plan stage 0 item 3)' if worst * 2.0 ** -53 > 1e-6
             else 'not needed'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
