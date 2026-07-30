"""Stage 7a of `doc/cancellation_ranking_plan.md` -- calibrate de-cancellation.

Before attempting this on a real amplifier, verify it on a circuit whose
cancellation-free answer is *published and independently re-derivable*.  Song &
Shi (ASP-DAC 2012) §II give a four-conductance ladder whose determinant over
composite matrix entries is ``adg - aef - bcg`` and whose expansion over the
actual parameters is ``G1G2G3 + G1G2G4 + G1G3G4`` -- three terms, all positive,
so ``kappa = 1`` exactly.

Two things make this a real calibration rather than a restatement:

* the target was re-derived here from the matrix, and separately as the **spanning
  trees** of the circuit graph, which agree with the paper term for term;
* that agreement identifies the mechanism.  For a nodal admittance matrix the
  cancellation-free expansion is Kirchhoff's spanning-tree expansion, so the
  de-cancellation rule is simply **each device symbol at most once along a term**.
  A floating two-terminal device stamps ``+g`` at (i,i) and (j,j) and ``-g`` at
  (i,j) and (j,i); a term taking (i,i)(j,j) and one taking (i,j)(j,i) have equal
  magnitude and opposite permutation sign, so they cancel in pairs.

The second measurement is the hazard, and it is the reason this stage exists
separately from applying the idea.  Tan & Shi (TCAD 2004) state that
"cancellation-free s-expanded DDDs do not satisfy Theorem 1", i.e.
**de-cancellation destroys DDD canonicity**.  Canonicity is what lets a diagram
share a minor between paths; without it the memo key has to carry the set of
devices already used, and the sharing that makes a DDD compact is lost.  That cost
is measured here on RC ladders, before anything is built on top of it.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/decancellation_calibration.py
"""

import itertools
import sys
import time

import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix


def entry_addends(expr, devices):
    """Split one matrix entry into its per-device additive contributions.

    Returns ``[(coefficient, device_symbol_or_None), ...]``.  This is the
    full-symbol (one-symbol-per-device) view of an entry that the compact-symbol
    DDD keeps as a single opaque payload.
    """
    out = []
    for term in sympy.Add.make_args(sympy.expand(expr)):
        touched = term.free_symbols & devices
        if len(touched) > 1:
            raise ValueError('addend %s touches %d devices, so the entry is not '
                             'linear in the device parameters' % (term, len(touched)))
        out.append((term, next(iter(touched)) if touched else None))
    return out


def expand_det(A, devices, rows, cols, used, prune):
    """Laplace expansion over per-device addends, yielding signed monomials.

    ``prune`` applies the de-cancellation rule: a device may appear at most once
    in a term.  With it off, this is the plain full-symbol expansion, which is
    what the compact diagram would give if its entries were multiplied out.
    """
    if not rows:
        yield sympy.Integer(1)
        return
    r = rows[0]
    for q, c in enumerate(cols):
        if A[r, c] == 0:
            continue
        sign = -1 if q % 2 else 1
        for coeff, dev in entry_addends(A[r, c], devices):
            if prune and dev is not None and dev in used:
                continue
            nxt = used | {dev} if dev is not None else used
            for rest in expand_det(A, devices, rows[1:], cols[:q] + cols[q + 1:],
                                   nxt, prune):
                yield sign * coeff * rest


def terms_of(A, devices, prune):
    n = A.rows
    return list(expand_det(A, devices, tuple(range(n)), tuple(range(n)),
                           frozenset(), prune))


def kappa_of(terms, env):
    vals = [complex(t.xreplace(env)) for t in terms]
    total = sum(vals)
    mass = sum(abs(v) for v in vals)
    if mass == 0:
        return 0.0, total, 1.0
    return mass, total, (float('inf') if total == 0 else mass / abs(total))


def count_states(A, devices, prune):
    """Distinct memo states, i.e. how much sharing survives.

    Without pruning the state is the minor ``(rows, cols)`` -- DDD canonicity.
    With pruning it must also carry the devices already consumed, because what
    lies below a minor now depends on the path taken to reach it.  The ratio of
    the two counts is the price of de-cancellation.
    """
    memo = {}

    def rec(rows, cols, used):
        key = (rows, cols, used) if prune else (rows, cols)
        if key in memo:
            return memo[key]
        if not rows:
            memo[key] = 1
            return 1
        r = rows[0]
        total = 0
        for q, c in enumerate(cols):
            if A[r, c] == 0:
                continue
            for _coeff, dev in entry_addends(A[r, c], devices):
                if prune and dev is not None and dev in used:
                    continue
                nxt = used | {dev} if dev is not None else used
                total += rec(rows[1:], cols[:q] + cols[q + 1:], nxt)
        memo[key] = total
        return total

    n = A.rows
    count = rec(tuple(range(n)), tuple(range(n)), frozenset())
    return len(memo), count


def calibrate():
    G1, G2, G3, G4 = sympy.symbols('G1 G2 G3 G4', positive=True)
    devices = {G1, G2, G3, G4}
    ## Song & Shi Fig. 1 / Eq. (1): Iin into node 1, G1 from 1 to 2, G2 from 2 to
    ## ground, G3 from 2 to 3, G4 from 3 to ground, Vout = V3.
    A = sympy.Matrix([[G1, -G1, 0],
                      [-G1, G1 + G2 + G3, -G3],
                      [0, -G3, G3 + G4]])
    target = G1 * G2 * G3 + G1 * G2 * G4 + G1 * G3 * G4
    env = {G1: 1.0e-3, G2: 2.0e-4, G3: 5.0e-3, G4: 1.0e-4}

    print('== calibration: Song & Shi (ASP-DAC 2012) four-conductance ladder ==')
    print('  exact det (sympy)          : %s' % sympy.expand(A.det()))
    print('  published / re-derived     : %s' % target)
    print('  agree                      : %s'
          % (sympy.simplify(sympy.expand(A.det()) - target) == 0))

    ## The spanning trees, computed from the topology, as an independent check on
    ## *why* those three terms are the answer.
    edges = {'G1': (1, 2), 'G2': (2, 0), 'G3': (2, 3), 'G4': (3, 0)}
    trees = []
    for combo in itertools.combinations(edges, 3):
        parent = {k: k for k in (0, 1, 2, 3)}

        def find(x):
            while parent[x] != x:
                x = parent[x]
            return x
        ok = True
        for name in combo:
            u, v = edges[name]
            ru, rv = find(u), find(v)
            if ru == rv:
                ok = False
                break
            parent[ru] = rv
        if ok:
            trees.append('*'.join(sorted(combo)))
    print('  spanning trees of the graph: %s' % ', '.join(sorted(trees)))
    print()

    for label, prune in (('full symbols, NO de-cancellation', False),
                         ('full symbols, DE-CANCELLED', True)):
        terms = terms_of(A, devices, prune)
        mass, total, kap = kappa_of(terms, env)
        got = sympy.expand(sympy.Add(*terms))
        print('  %s' % label)
        print('    terms %d,  sum|term| %.6e,  |sum| %.6e,  kappa %.6f'
              % (len(terms), mass, abs(total), kap))
        print('    symbolic sum: %s' % got)
        print('    equals the exact determinant: %s'
              % (sympy.simplify(got - target) == 0))
        if prune:
            print('    GATE 7a-1 (determinant preserved): %s'
                  % ('PASS' if sympy.simplify(got - target) == 0 else 'FAIL'))
            print('    GATE 7a-2 (kappa == 1): %s  (%.10f)'
                  % ('PASS' if abs(kap - 1.0) < 1e-12 else 'FAIL', kap))
        print()

    ## And the composite-symbol diagram, for the contrast that started all this.
    D = ddd_of_matrix(A)
    print('  compact-symbol DDD of the same matrix: %d vertices, %d terms,'
          ' kappa %.6f' % (D.size, D.term_count(), D.cancellation(env)))
    print('  -> the composite form cancels; the de-cancelled full-symbol form'
          ' does not.')
    print()


def sharing_cost():
    print('== the price of de-cancellation: how much sharing survives ==')
    print('  A de-cancelled expansion cannot key its memo on the minor alone,')
    print('  because what lies below depends on which devices the path already')
    print('  used.  That is the canonicity loss Tan & Shi (TCAD 2004) record.')
    print()
    print('  %-4s %-9s %-13s %-13s %-9s %s'
          % ('N', 'dim', 'states plain', 'states pruned', 'ratio', 'terms kept'))
    for N in (3, 4, 5, 6, 7):
        system = bc.rc_ladder(N)
        A = system.A
        devices = {s for s in A.free_symbols if s is not system.s}
        t0 = time.time()
        plain_states, plain_terms = count_states(A, devices, False)
        pruned_states, pruned_terms = count_states(A, devices, True)
        print('  %-4d %-9d %-13d %-13d %-9.1fx %d/%d  (%.1f s)'
              % (N, A.rows, plain_states, pruned_states,
                 pruned_states / plain_states, pruned_terms, plain_terms,
                 time.time() - t0))
    print()
    print('  Read: the ratio column is the sharing lost.  If it grows with N the')
    print('  de-cancelled form is not a compact representation, and the route has')
    print('  to fold de-cancellation into construction differently -- which is')
    print('  what Tan & Shi do by de-cancelling during the s-expansion instead.')


def main():
    calibrate()
    sharing_cost()
    return 0


if __name__ == '__main__':
    sys.exit(main())
