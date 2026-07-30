"""Stage 7g -- rank the de-cancelled expansion.  The payoff test.

Stages 7a-7f measured properties of the de-cancelled full-symbol expansion: it
preserves the determinant, costs a constant factor in memo states, and brings the
uA741's ``kappa`` from 6.66e3 to 99 at a fixed frequency and to 14-17 per
coefficient of ``s``.  None of that is an answer until something is *ranked*.

This ranks it, two ways, against a control measured at the same operating point:

* **group ranking** -- a prefix plus a state, contributing ``prod(prefix) *
  value(state)`` exactly, as `DDD.approximate_groups` does on the compact diagram;
* **magnitude ranking** -- a prefix plus an upper bound on what any completion can
  reach, as `DDD.approximate` does.  Tan & Shi's published pipeline is
  magnitude-ranked, so this is the comparison that speaks to their result.

**Scope cut, stated rather than skipped.** The ranking runs over the de-cancelled
*state graph* rather than a materialised `DDD`.  A per-coefficient diagram would be
several million vertices under first-row pivoting, and the question here is whether
the algorithm converges, not whether the object can be built.  The consequence is
that this duplicates logic `approximate_groups` already has; the mitigation is that
the kept sum is checked against the exact de-cancelled value, which is itself
checked against `numpy.linalg.det`.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/decancellation_ranking.py
"""

import heapq
import math
import sys
import time

import numpy as np

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix, _resolve

from decancellation_calibration import entry_addends, device_positions


def _is_complete(state):
    """A finished term, as opposed to a dead end.

    ``state[0]`` is the remaining rows.  Empty means every row has been assigned,
    so the path is one product term with value 1 below it.  **Non-empty with no
    children is a dead end** -- de-cancellation pruned every remaining option --
    and it contributes exactly zero.

    Conflating the two is not a small error: a dead end's prefix product is a real
    number, often a large one, so counting it as a term adds a contribution that
    the exact expansion does not contain.  It is what made magnitude ranking report
    a relative error of 6.7e+55, which is impossible when the total absolute mass
    is only 99x the value -- the impossibility is what exposed it.
    """
    return not state[0]


class Graph:
    """The de-cancelled expansion as an explicit DAG of states.

    A state is ``(rows, cols, live)`` where ``live`` is the set of still-reachable
    forbidden labels -- the keying that stages 7b-7e showed costs a constant factor
    rather than growing geometrically.  Each state carries its children, the exact
    value of the sub-expansion below it, and that sub-expansion's absolute mass, so
    both rankings have what they need without re-walking.
    """

    def __init__(self, A, devices, positions, env):
        self.addends = {}
        for i in range(A.rows):
            for j in range(A.cols):
                if A[i, j] == 0:
                    continue
                self.addends[(i, j)] = [
                    (complex(_resolve(co, env)), d, co)
                    for co, d in entry_addends(A[i, j], devices)]
        self.positions = positions
        self.n = A.rows
        self.children = {}          # state -> [(sign*coeff, label, child), ...]
        self.value = {}
        self.mass = {}
        self._build(tuple(range(self.n)), tuple(range(self.n)), frozenset())

    def _canon(self, rows, cols, forbidden):
        rowset, colset = frozenset(rows), frozenset(cols)
        return frozenset((i, j, d) for (i, j, d) in forbidden
                         if i in rowset and j in colset)

    def _build(self, rows, cols, forbidden):
        live = self._canon(rows, cols, forbidden)
        key = (rows, cols, live)
        if key in self.value:
            return key
        if not rows:
            self.children[key] = []
            self.value[key] = 1.0 + 0j
            self.mass[key] = 1.0
            return key
        ## Reserve the slot before recursing so a cycle would be visible; the
        ## expansion is acyclic (each step removes a row) but the guard is cheap.
        self.value[key] = None
        kids = []
        r = rows[0]
        for q, c in enumerate(cols):
            entry = self.addends.get((r, c))
            if entry is None:
                continue
            sign = -1 if q % 2 else 1
            for coeff, dev, sym in entry:
                if dev is not None and (r, c, dev) in live:
                    continue
                nxt = live
                if dev is not None:
                    nxt = live | frozenset(
                        (i, j, dev) for (i, j) in self.positions[dev]
                        if (i, j) != (r, c))
                child = self._build(rows[1:], cols[:q] + cols[q + 1:], nxt)
                kids.append((sign * coeff, sym, child))
        val, mass = 0.0 + 0j, 0.0
        for coeff, _sym, child in kids:
            val += coeff * self.value[child]
            mass += abs(coeff) * self.mass[child]
        self.children[key] = kids
        self.value[key] = val
        self.mass[key] = mass
        return key

    @property
    def root(self):
        return (tuple(range(self.n)), tuple(range(self.n)), frozenset())

    def size(self):
        return len(self.value)


def rank_groups(graph, tol, cap=400000):
    """Best-first over groups, ranked by exact contribution.

    ``kept + everything still on the frontier == the total``, exactly, so
    ``|total - kept| / |total|`` is the kept expression's own error at every step.
    Same invariant `DDD.approximate_groups` maintains.
    """
    total = graph.value[graph.root]
    heap = [(-abs(total), 0, graph.root, 1.0 + 0j)]
    counter = 0
    kept, kept_sum, splits = 0, 0.0 + 0j, 0
    while heap:
        err = abs(total - kept_sum) / abs(total)
        if err <= tol or splits >= cap:
            break
        neg, _, state, acc = heapq.heappop(heap)
        if -neg == 0.0:
            continue
        kids = graph.children[state]
        if not kids:
            if _is_complete(state):
                kept += 1
                kept_sum += acc
            continue                        # else a dead end: contributes zero
        splits += 1
        for coeff, _sym, child in kids:
            contribution = acc * coeff * graph.value[child]
            if contribution == 0:
                continue
            if not graph.children[child]:
                if _is_complete(child):
                    kept += 1
                    kept_sum += acc * coeff
                continue                    # else a dead end: contributes zero
            counter += 1
            heapq.heappush(heap, (-abs(contribution), counter, child,
                                  acc * coeff))
    return kept, abs(total - kept_sum) / abs(total), splits


def rank_magnitude(graph, tol, cap=400000):
    """Best-first over terms, ranked by the largest product still reachable.

    This is what `DDD.approximate` does: the bound is the mass-style optimum below
    a state, so terms come out in decreasing magnitude.  The stopping test is still
    the exact signed error -- criterion (8) of the literature.
    """
    ## A completed term contributes 1; a dead end contributes 0, and must not be
    ## given 1 or the bound sends the search straight into it.
    best = {state: (1.0 if _is_complete(state) else 0.0)
            for state in graph.value}
    ## Shortest remaining-row list first, so children are finished before parents.
    for state in sorted(graph.value, key=lambda k: len(k[0])):
        kids = graph.children[state]
        if kids:
            best[state] = max(abs(co) * best[ch] for co, _s, ch in kids)

    total = graph.value[graph.root]
    heap = [(-best[graph.root], 0, graph.root, 1.0 + 0j)]
    counter = 0
    kept, kept_sum, steps = 0, 0.0 + 0j, 0
    while heap:
        err = abs(total - kept_sum) / abs(total)
        if err <= tol or steps >= cap:
            break
        neg, _, state, acc = heapq.heappop(heap)
        steps += 1
        kids = graph.children[state]
        if not kids:
            if _is_complete(state):
                kept += 1
                kept_sum += acc
            continue                        # else a dead end: contributes zero
        for coeff, _sym, child in kids:
            nacc = acc * coeff
            bound = abs(nacc) * best[child]
            if bound == 0.0:
                continue
            counter += 1
            heapq.heappush(heap, (-bound, counter, child, nacc))
    return kept, abs(total - kept_sum) / abs(total), steps


def main():
    tol = 0.05
    system = bc.ua741(fully_symbolic=True)
    A = system.A
    devices = {s for s in A.free_symbols if s is not system.s}
    positions = device_positions(A, devices)
    env = {system.s: 2j * math.pi * 1e3}
    for sym, value in system.params.items():
        if sym is not system.s:
            env[sym] = complex(value)

    print('uA741, fully symbolic, 1 kHz, tol = %.2f' % tol)
    print()

    ## -- the control, measured HERE and not quoted -----------------------
    numeric = bc.ua741()
    cenv = {numeric.s: env[system.s]}
    D = ddd_of_matrix(A)
    k_compact = D.cancellation(env)
    print('control: compact-symbol diagram at the SAME operating point')
    print('  %d vertices, %d terms, kappa %.4e'
          % (D.size, D.term_count(), k_compact))
    t0 = time.time()
    _e, n_grp, err_grp = D.approximate_groups(env, tol=tol)
    print('  group ranking     : %6d terms, err %.3e, %.1f s'
          % (n_grp, err_grp, time.time() - t0))
    import warnings
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        _e, n_mag, err_mag = D.approximate(env, tol=tol, max_terms=200000)
    print('  magnitude ranking : %6d terms, err %.3e, %.1f s'
          % (n_mag, err_mag, time.time() - t0))
    print()

    ## -- the de-cancelled expansion -------------------------------------
    t0 = time.time()
    graph = Graph(A, devices, positions, env)
    build = time.time() - t0
    total = graph.value[graph.root]
    num = np.array([[complex(_resolve(A[i, j], env)) for j in range(A.cols)]
                    for i in range(A.rows)], dtype=complex)
    exact = np.linalg.det(num)
    rel = abs(total - exact) / abs(exact)
    k_dec = graph.mass[graph.root] / abs(total)
    print('de-cancelled expansion: %d states in %.1f s' % (graph.size(), build))
    print('  value vs numpy.linalg.det: rel %.2e  %s'
          % (rel, 'OK' if rel < 1e-9 else 'MISMATCH'))
    print('  kappa %.4e   (compact %.4e, %.1fx better)'
          % (k_dec, k_compact, k_compact / k_dec))

    t0 = time.time()
    d_grp, d_grp_err, splits = rank_groups(graph, tol)
    print('  group ranking     : %6d terms, err %.3e, %d splits, %.1f s'
          % (d_grp, d_grp_err, splits, time.time() - t0))
    t0 = time.time()
    d_mag, d_mag_err, steps = rank_magnitude(graph, tol)
    print('  magnitude ranking : %6d terms, err %.3e, %d steps, %.1f s'
          % (d_mag, d_mag_err, steps, time.time() - t0))
    print()

    print('== gates ==')
    print('  7g-4 (kept sum converges to the verified value): %s'
          % ('PASS' if rel < 1e-9 and d_grp_err <= tol else 'see errors above'))
    ratio = n_grp / d_grp if d_grp else float('inf')
    print('  7g-2 (group ranking >= 5x fewer terms de-cancelled): %s'
          '   %d -> %d, %.1fx'
          % ('PASS' if ratio >= 5.0 else 'FAIL', n_grp, d_grp, ratio))
    print('  7g-3 magnitude ranking, compact %d terms (err %.3e) vs'
          ' de-cancelled %d terms (err %.3e)'
          % (n_mag, err_mag, d_mag, d_mag_err))
    if err_mag > tol >= d_mag_err:
        print('       -> de-cancellation makes plain magnitude ranking converge')
        print('          where it could not before: Tan & Shi reproduced.')
    return 0


if __name__ == '__main__':
    sys.setrecursionlimit(10000)
    sys.exit(main())
