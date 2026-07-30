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

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix, _resolve


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


def device_positions(A, devices):
    """``device -> frozenset of (row, col)`` it stamps into."""
    pos = {d: set() for d in devices}
    for i in range(A.rows):
        for j in range(A.cols):
            if A[i, j] == 0:
                continue
            for _coeff, dev in entry_addends(A[i, j], devices):
                if dev is not None:
                    pos[dev].add((i, j))
    return {d: frozenset(v) for d, v in pos.items()}


def count_states_canonical(A, devices, positions):
    """De-cancelled state count, keyed on STILL-REACHABLE forbidden positions.

    The naive predicate carries every device used so far.  Most of that is dead
    information: once a device has been consumed at ``(p,p)``, its partners
    ``(p,q)`` and ``(q,p)`` need row or column ``p``, which the minor no longer
    has, so at most one partner position stays reachable -- and it stops mattering
    as soon as the minor drops its row or its column.

    Keying on the reachable remainder is what Tan & Shi's ``REMAINDER`` set
    operation does implicitly, and it is the difference between a path-dependent
    state and a diagram one.  The surviving-term count must come out identical to
    the naive predicate, which is asserted by the caller.

    **The state is a set of LABELS, not of positions**, and the distinction is not
    cosmetic: a first attempt forbade matrix *positions* and silently lost terms,
    because one position can carry several devices' addends -- entry (2,2) of an RC
    ladder is ``C2*s + 1/R1``, so banning the position to keep ``R1`` out also
    banned ``C2``.  This is why the literature speaks of a canceling list per
    *label*, a label being one (device, position) pair.  Gate 7b-1 caught it.
    """
    memo = {}

    def rec(rows, cols, forbidden):
        rowset, colset = frozenset(rows), frozenset(cols)
        ## Drop any forbidden label the minor can no longer reach.  This is the
        ## whole trick: it makes the state depend on the minor, not the path.
        live = frozenset((i, j, d) for (i, j, d) in forbidden
                         if i in rowset and j in colset)
        key = (rows, cols, live)
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
                if dev is not None and (r, c, dev) in live:
                    continue              # this label is de-cancelled away
                nxt = live
                if dev is not None:
                    ## Forbid this device's other labels from here down.
                    nxt = live | frozenset((i, j, dev)
                                           for (i, j) in positions[dev]
                                           if (i, j) != (r, c))
                total += rec(rows[1:], cols[:q] + cols[q + 1:], nxt)
        memo[key] = total
        return total

    n = A.rows
    count = rec(tuple(range(n)), tuple(range(n)), frozenset())
    return len(memo), count


def polynomial_check():
    print('== stage 7b: is determinant-side de-cancellation polynomial? ==')
    print('  Same de-cancellation, but the memo key carries only the')
    print('  STILL-REACHABLE forbidden positions instead of every device used.')
    print()
    print('  %-4s %-6s %-9s %-11s %-11s %-8s %-8s %s'
          % ('N', 'dim', 'plain', 'naive key', 'canon key', 'naive x', 'canon x',
             'terms'))
    naive_ratios, canon_ratios = [], []
    for N in range(3, 11):
        system = bc.rc_ladder(N)
        A = system.A
        devices = {s for s in A.free_symbols if s is not system.s}
        positions = device_positions(A, devices)
        plain_states, plain_terms = count_states(A, devices, False)
        naive_states, naive_terms = count_states(A, devices, True)
        canon_states, canon_terms = count_states_canonical(A, devices, positions)
        ## GATE 7b-1: the canonicalisation must not change the answer.
        assert canon_terms == naive_terms, (
            'canonicalised key kept %d terms, naive kept %d -- the '
            'canonicalisation is wrong' % (canon_terms, naive_terms))
        naive_ratios.append(naive_states / plain_states)
        canon_ratios.append(canon_states / plain_states)
        print('  %-4d %-6d %-9d %-11d %-11d %-8.1f %-8.2f %d/%d'
              % (N, A.rows, plain_states, naive_states, canon_states,
                 naive_states / plain_states, canon_states / plain_states,
                 canon_terms, plain_terms))
    print()
    print('  GATE 7b-1 (canonicalisation preserves the term count): PASS'
          ' -- asserted at every N')
    flat = canon_ratios[-1] <= 3.0 and canon_ratios[-1] <= 1.5 * canon_ratios[0]
    print('  GATE 7b-2 (canon ratio <= 3x and flat in N): %s'
          '  -- %.2fx at N=3 rising to %.2fx at N=10'
          % ('PASS' if flat else 'FAIL', canon_ratios[0], canon_ratios[-1]))
    print('  for contrast, the naive key went %.1fx -> %.1fx'
          % (naive_ratios[0], naive_ratios[-1]))
    print()
    if flat:
        print('  Reading: determinant-side de-cancellation IS polynomial once the')
        print('  state is the reachable remainder rather than the path. Building')
        print('  the REMAINDER/CL diagram construction is justified.')
    else:
        print('  Reading: it is NOT rescued by canonicalisation, so Tan & Shi\'s')
        print('  REMAINDER construction cannot make it compact either, and the')
        print('  topological route is the answer rather than one option of three.')


def decancelled_value(A, devices, positions, env):
    """Signed value and absolute mass of the de-cancelled full-symbol expansion.

    Same recursion and same label-canonical state as `count_states_canonical`,
    but carrying numbers instead of counting.  ``kappa`` is the ratio of the two
    returns, and the signed value must equal ``det(A)`` -- which is the check that
    the de-cancellation removed only genuinely cancelling pairs and nothing else.
    """
    addends = {}
    for i in range(A.rows):
        for j in range(A.cols):
            if A[i, j] == 0:
                continue
            addends[(i, j)] = [(complex(_resolve(co, env)), d)
                               for co, d in entry_addends(A[i, j], devices)]
    memo = {}

    def rec(rows, cols, forbidden):
        rowset, colset = frozenset(rows), frozenset(cols)
        live = frozenset((i, j, d) for (i, j, d) in forbidden
                         if i in rowset and j in colset)
        key = (rows, cols, live)
        hit = memo.get(key)
        if hit is not None:
            return hit
        if not rows:
            memo[key] = (1.0 + 0j, 1.0)
            return memo[key]
        r = rows[0]
        val, mass = 0.0 + 0j, 0.0
        for q, c in enumerate(cols):
            entry = addends.get((r, c))
            if entry is None:
                continue
            sign = -1 if q % 2 else 1
            for coeff, dev in entry:
                if dev is not None and (r, c, dev) in live:
                    continue
                nxt = live
                if dev is not None:
                    nxt = live | frozenset((i, j, dev)
                                           for (i, j) in positions[dev]
                                           if (i, j) != (r, c))
                v, m = rec(rows[1:], cols[:q] + cols[q + 1:], nxt)
                val += sign * coeff * v
                mass += abs(coeff) * m
        memo[key] = (val, mass)
        return memo[key]

    n = A.rows
    val, mass = rec(tuple(range(n)), tuple(range(n)), frozenset())
    return val, mass, len(memo)


def kappa_on_real_circuits():
    print('== stage 7c: does de-cancellation move kappa on a real circuit? ==')
    print('  Verified against numpy.linalg.det at every N -- if the de-cancelled')
    print('  sum is not the determinant, the rule removed something it should not.')
    print()
    print('  %-4s %-5s %-11s %-11s %-9s %-10s %s'
          % ('N', 'dim', 'kappa cmpct', 'kappa decan', 'change', 'states',
             'det check'))
    for N in (4, 5, 6, 7, 8, 9):
        system = bc.rc_ladder(N)
        A = system.A
        devices = {s for s in A.free_symbols if s is not system.s}
        positions = device_positions(A, devices)
        env = {system.s: 2j * 3.14159265358979 * 1e3}
        for k, sym in enumerate(sorted(devices, key=str)):
            name = str(sym)
            env[sym] = ((100.0 * 2.0 ** k) if name.startswith('R')
                        else (1e-9 * 2.0 ** k))
        D = ddd_of_matrix(A)
        k_compact = D.cancellation(env)
        val, mass, states = decancelled_value(A, devices, positions, env)
        k_decan = float('inf') if val == 0 else mass / abs(val)
        num = np.array([[complex(_resolve(A[i, j], env)) for j in range(A.cols)]
                        for i in range(A.rows)], dtype=complex)
        exact = np.linalg.det(num)
        rel = abs(val - exact) / abs(exact)
        print('  %-4d %-5d %-11.4g %-11.4g %-9.2fx %-10d %.1e %s'
              % (N, A.rows, k_compact, k_decan, k_compact / k_decan, states,
                 rel, 'OK' if rel < 1e-9 else 'MISMATCH'))
    print()
    print('  A MISMATCH is not a bug in the arithmetic: it would mean the')
    print('  "each device once" rule is incomplete for this formulation. MNA adds')
    print('  a row per voltage source whose +-1 incidence entries carry the same')
    print('  crossing pattern a device does, but they are constants rather than')
    print('  symbols, so the rule cannot see them.')


def decancelled_by_degree(A, devices, positions, env, s):
    """De-cancelled expansion, split by power of ``s``: (value, mass) per degree.

    Stage 7c found the de-cancelled kappa sitting near 2 on an RC ladder even
    though every surviving term is a spanning tree.  That residue is not sign
    cancellation -- it is **phase**: at a fixed complex ``s`` a conductance term is
    real and a ``C*s`` term is imaginary, so the moduli cannot add up to the
    modulus of the sum.  Separating the powers of ``s`` puts them in different
    coefficients, and within one coefficient of a passive network every spanning
    tree product is positive.

    So the prediction is ``kappa == 1`` per coefficient, exactly -- and if it
    holds, it is why Tan & Shi s-expand *before* de-cancelling rather than after.
    """
    addends = {}
    for i in range(A.rows):
        for j in range(A.cols):
            if A[i, j] == 0:
                continue
            out = []
            for co, d in entry_addends(A[i, j], devices):
                deg = sympy.Poly(co, s).degree() if co.has(s) else 0
                coeff = co / s ** deg if deg else co
                out.append((complex(_resolve(coeff, env)), d, deg))
            addends[(i, j)] = out
    memo = {}

    def rec(rows, cols, forbidden):
        rowset, colset = frozenset(rows), frozenset(cols)
        live = frozenset((i, j, d) for (i, j, d) in forbidden
                         if i in rowset and j in colset)
        key = (rows, cols, live)
        hit = memo.get(key)
        if hit is not None:
            return hit
        if not rows:
            memo[key] = {0: (1.0 + 0j, 1.0)}
            return memo[key]
        r = rows[0]
        acc = {}
        for q, c in enumerate(cols):
            entry = addends.get((r, c))
            if entry is None:
                continue
            sign = -1 if q % 2 else 1
            for coeff, dev, deg in entry:
                if dev is not None and (r, c, dev) in live:
                    continue
                nxt = live
                if dev is not None:
                    nxt = live | frozenset((i, j, dev)
                                           for (i, j) in positions[dev]
                                           if (i, j) != (r, c))
                below = rec(rows[1:], cols[:q] + cols[q + 1:], nxt)
                for d2, (v, m) in below.items():
                    tot = d2 + deg
                    pv, pm = acc.get(tot, (0.0 + 0j, 0.0))
                    acc[tot] = (pv + sign * coeff * v, pm + abs(coeff) * m)
        memo[key] = acc
        return acc

    n = A.rows
    return rec(tuple(range(n)), tuple(range(n)), frozenset())


def per_coefficient_kappa():
    print('== stage 7d: de-cancellation TOGETHER WITH s-expansion ==')
    print('  Prediction from 7c: the residual kappa there was phase, not sign, so')
    print('  splitting the powers of s should give kappa == 1 in every coefficient.')
    print()
    for N in (4, 6, 8):
        system = bc.rc_ladder(N)
        A = system.A
        devices = {s for s in A.free_symbols if s is not system.s}
        positions = device_positions(A, devices)
        env = {}
        for k, sym in enumerate(sorted(devices, key=str)):
            name = str(sym)
            env[sym] = ((100.0 * 2.0 ** k) if name.startswith('R')
                        else (1e-9 * 2.0 ** k))
        by_deg = decancelled_by_degree(A, devices, positions, env, system.s)
        worst = 0.0
        print('  rc_ladder(%d), dim %d: %d coefficients' % (N, A.rows, len(by_deg)))
        for deg in sorted(by_deg):
            v, m = by_deg[deg]
            if m == 0:
                continue
            kap = float('inf') if v == 0 else m / abs(v)
            if kap != float('inf'):
                worst = max(worst, kap)
            print('    s^%-3d value %13.6e  mass %13.6e  kappa %.12f'
                  % (deg, abs(v), m, kap))
        print('    worst kappa over coefficients: %.12f  -> %s'
              % (worst, 'PASS (== 1)' if abs(worst - 1.0) < 1e-9 else 'FAIL'))
        print()


def main():
    calibrate()
    sharing_cost()
    polynomial_check()
    kappa_on_real_circuits()
    per_coefficient_kappa()
    return 0


if __name__ == '__main__':
    sys.exit(main())
