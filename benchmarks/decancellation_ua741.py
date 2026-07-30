"""Stage 7e -- de-cancellation on a fully symbolic uA741.

`decancellation_calibration.py` established, on RC ladders, that de-cancellation
keyed on the still-reachable **labels** costs a flat ~1.2x in memo states (against
170x for a path-keyed state), preserves the determinant, and -- combined with
s-expansion -- gives ``kappa = 1`` exactly.

The declared reconsider-if was that ladders are a benign topology: every minor
borders few devices, so few labels stay live.  A transistor amplifier is denser.
This script tests that on `ua741(fully_symbolic=True)`, which exists for this
purpose -- with the default fixture 168 of 168 additive contributions are pure
numbers, so there is nothing for a device-level rule to see.

Everything is capped, because the honest outcome of this measurement may be "it
does not fit", and a script that runs for hours to say so is worse than one that
says so in a minute.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/decancellation_ua741.py
"""

import math
import sys
import time

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix, s_expand, _resolve

from decancellation_calibration import entry_addends, device_positions


STATE_CAP = 2_000_000
TIME_CAP = 900.0


class Exceeded(Exception):
    pass


def count_states(A, devices, positions, prune, cap=STATE_CAP,
                 time_cap=TIME_CAP):
    """Memo-state count, with and without label-keyed de-cancellation.

    Raises `Exceeded` rather than running away: the point of the measurement is
    the ratio, and if it cannot be reached that is itself the answer.
    """
    memo = {}
    start = time.time()

    def rec(rows, cols, forbidden):
        if len(memo) > cap:
            raise Exceeded('state cap %d exceeded' % cap)
        if time.time() - start > time_cap:
            raise Exceeded('time cap %.0f s exceeded' % time_cap)
        if prune:
            rowset, colset = frozenset(rows), frozenset(cols)
            live = frozenset((i, j, d) for (i, j, d) in forbidden
                             if i in rowset and j in colset)
            key = (rows, cols, live)
        else:
            live = forbidden
            key = (rows, cols)
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
                if prune and dev is not None and (r, c, dev) in live:
                    continue
                nxt = live
                if prune and dev is not None:
                    nxt = live | frozenset((i, j, dev)
                                           for (i, j) in positions[dev]
                                           if (i, j) != (r, c))
                total += rec(rows[1:], cols[:q] + cols[q + 1:], nxt)
        memo[key] = total
        return total

    n = A.rows
    terms = rec(tuple(range(n)), tuple(range(n)), frozenset())
    return len(memo), terms, time.time() - start


def decancelled_value(A, devices, positions, env, cap=STATE_CAP,
                      time_cap=TIME_CAP):
    addends = {}
    for i in range(A.rows):
        for j in range(A.cols):
            if A[i, j] == 0:
                continue
            addends[(i, j)] = [(complex(_resolve(co, env)), d)
                               for co, d in entry_addends(A[i, j], devices)]
    memo = {}
    start = time.time()

    def rec(rows, cols, forbidden):
        if len(memo) > cap:
            raise Exceeded('state cap %d exceeded' % cap)
        if time.time() - start > time_cap:
            raise Exceeded('time cap %.0f s exceeded' % time_cap)
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
    return val, mass, len(memo), time.time() - start


def decancelled_by_degree(A, devices, positions, env, s, cap=STATE_CAP,
                          time_cap=TIME_CAP):
    """As `decancelled_value`, split by power of ``s``.

    Stage 7d showed ``kappa`` reaches 1 per coefficient on a *passive* network,
    because every surviving spanning-tree product is then positive.  An amplifier
    is not passive: its controlled sources make the surviving terms signed (common
    spanning trees of two different graphs), so 1 is not the expectation here.
    What the coefficient split should remove is the *phase* spread between real
    conductances and imaginary ``C*s`` terms.  The residue is the honest measure
    of how much sign cancellation survives de-cancellation on a real circuit.
    """
    addends = {}
    for i in range(A.rows):
        for j in range(A.cols):
            if A[i, j] == 0:
                continue
            out = []
            for co, d in entry_addends(A[i, j], devices):
                deg = sympy.Poly(co, s).degree() if co.has(s) else 0
                base = co / s ** deg if deg else co
                out.append((complex(_resolve(base, env)), d, deg))
            addends[(i, j)] = out
    memo = {}
    start = time.time()

    def rec(rows, cols, forbidden):
        if len(memo) > cap:
            raise Exceeded('state cap %d exceeded' % cap)
        if time.time() - start > time_cap:
            raise Exceeded('time cap %.0f s exceeded' % time_cap)
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
                for d2, (v, m) in rec(rows[1:], cols[:q] + cols[q + 1:],
                                      nxt).items():
                    tot = d2 + deg
                    pv, pm = acc.get(tot, (0.0 + 0j, 0.0))
                    acc[tot] = (pv + sign * coeff * v, pm + abs(coeff) * m)
        memo[key] = acc
        return acc

    n = A.rows
    return rec(tuple(range(n)), tuple(range(n)), frozenset()), len(memo), \
        time.time() - start


def per_coefficient(A, devices, positions, env, s, k_fixed):
    print()
    print('== stage 7f: de-cancelled AND s-expanded, per coefficient ==')
    try:
        by_deg, states, secs = decancelled_by_degree(A, devices, positions,
                                                     env, s)
    except Exceeded as exc:
        print('  EXCEEDED -- %s' % exc)
        return
    except RecursionError:
        print('  RecursionError')
        return
    kappas = []
    print('  %d coefficients, %d states, %.1f s' % (len(by_deg), states, secs))
    print('  %-6s %-15s %-15s %s' % ('s^k', '|coefficient|', 'mass', 'kappa'))
    for deg in sorted(by_deg):
        v, m = by_deg[deg]
        if m == 0:
            continue
        kap = float('inf') if v == 0 else m / abs(v)
        if math.isfinite(kap):
            kappas.append((kap, deg))
        print('  %-6d %-15.6e %-15.6e %.4e' % (deg, abs(v), m, kap))
    if kappas:
        worst, wdeg = max(kappas)
        best, bdeg = min(kappas)
        print('  worst kappa %.4e (s^%d), best %.4e (s^%d)'
              % (worst, wdeg, best, bdeg))
        print('  for comparison: compact-symbol whole determinant %.4e' % k_fixed)
        ## The control that makes the comparison fair: the SAME matrix and the
        ## same operating point, s-expanded but NOT de-cancelled.  Without this
        ## the de-cancelled per-coefficient numbers have nothing to be better or
        ## worse than.
        print()
        print('  control -- same matrix, same point, s-expanded, NOT'
              ' de-cancelled:')
        try:
            t0 = time.time()
            expanded = s_expand(A, s)
            cenv = {k: v for k, v in env.items() if k is not s}
            rows = []
            for k in range(expanded.degree + 1):
                coeff = expanded.coefficient(k)
                if coeff.root is None or coeff.term_count() == 0:
                    continue
                rows.append((k, coeff.cancellation(cenv)))
            print('    s_expand: %d vertices, %.1f s' % (expanded.size,
                                                         time.time() - t0))
            finite = [(kap, k) for k, kap in rows if math.isfinite(kap)]
            if finite:
                cw, cwd = max(finite)
                cb, cbd = min(finite)
                print('    compact per-coefficient kappa: worst %.4e (s^%d),'
                      ' best %.4e (s^%d)' % (cw, cwd, cb, cbd))
                pairs = {k: kap for k, kap in rows}
                print('    %-6s %-14s %-14s %s'
                      % ('s^k', 'compact', 'de-cancelled', 'improvement'))
                for kap, deg in sorted(kappas, key=lambda x: x[1]):
                    c = pairs.get(deg)
                    if c is None or not math.isfinite(c):
                        continue
                    print('    %-6d %-14.4e %-14.4e %.1fx'
                          % (deg, c, kap, c / kap))
        except Exception as exc:
            print('    control failed: %s: %s' % (type(exc).__name__, exc))
        print('  An amplifier is not passive, so kappa == 1 is NOT expected here')
        print('  -- controlled sources make the surviving terms signed. The')
        print('  question is only how much sign cancellation is left.')


def main():
    t0 = time.time()
    system = bc.ua741(fully_symbolic=True)
    A = system.A
    devices = {s for s in A.free_symbols if s is not system.s}
    positions = device_positions(A, devices)
    addend_count = sum(len(entry_addends(A[i, j], devices))
                       for i in range(A.rows) for j in range(A.cols)
                       if A[i, j] != 0)
    print('fully symbolic uA741: dim %d, %d device symbols, %d addends, %.1f s'
          % (A.rows, len(devices), addend_count, time.time() - t0))
    ## Labels per device: how many matrix positions each one occupies.  This is
    ## what decides how many stay live inside a minor, hence the overhead.
    sizes = sorted(len(v) for v in positions.values())
    print('  stamp positions per device: min %d, median %d, max %d'
          % (sizes[0], sizes[len(sizes) // 2], sizes[-1]))

    env = {system.s: 2j * math.pi * 1e3}
    for sym, value in system.params.items():
        if sym is not system.s:
            env[sym] = complex(value)

    t0 = time.time()
    D = ddd_of_matrix(A)
    k_compact = D.cancellation(env)
    print('  compact-symbol DDD: %d vertices, %d terms, kappa %.4e, %.1f s'
          % (D.size, D.term_count(), k_compact, time.time() - t0))
    print()

    print('== memo states, plain vs label-keyed de-cancellation ==')
    plain = pruned = None
    try:
        states, terms, secs = count_states(A, devices, positions, False)
        plain = states
        print('  plain (no de-cancellation): %d states, %d terms, %.1f s'
              % (states, terms, secs))
    except Exceeded as exc:
        print('  plain: EXCEEDED -- %s' % exc)
    except RecursionError:
        print('  plain: RecursionError')

    try:
        states, terms, secs = count_states(A, devices, positions, True)
        pruned = states
        print('  label-keyed de-cancelled : %d states, %d terms, %.1f s'
              % (states, terms, secs))
    except Exceeded as exc:
        print('  label-keyed: EXCEEDED -- %s' % exc)
    except RecursionError:
        print('  label-keyed: RecursionError')

    if plain and pruned:
        print('  overhead: %.2fx   (ladders gave a flat 1.2x)'
              % (pruned / plain))
        print('  GATE 7e-1 (overhead stays <= 3x on a dense circuit): %s'
              % ('PASS' if pruned / plain <= 3.0 else 'FAIL'))
    print()

    print('== kappa of the de-cancelled expansion ==')
    try:
        val, mass, states, secs = decancelled_value(A, devices, positions, env)
        k_dec = float('inf') if val == 0 else mass / abs(val)
        num = np.array([[complex(_resolve(A[i, j], env))
                         for j in range(A.cols)] for i in range(A.rows)],
                       dtype=complex)
        exact = np.linalg.det(num)
        rel = abs(val - exact) / abs(exact)
        print('  de-cancelled value %.6e%+.6ej  (%d states, %.1f s)'
              % (val.real, val.imag, states, secs))
        print('  vs numpy.linalg.det: rel %.2e  %s'
              % (rel, 'OK' if rel < 1e-9 else 'MISMATCH -- rule is unsound here'))
        print('  kappa: compact %.4e -> de-cancelled %.4e   (%.1fx)'
              % (k_compact, k_dec, k_compact / k_dec if k_dec else float('inf')))
        print('  GATE 7e-2 (determinant preserved): %s'
              % ('PASS' if rel < 1e-9 else 'FAIL'))
        print('  NOTE: at a fixed complex s a residual kappa is expected and is')
        print('  phase, not sign -- stage 7d showed kappa reaches 1 only per')
        print('  coefficient of s.  The number to compare against 9.4e3 is this')
        print('  one; the number that should be 1 needs the s-expansion.')
    except Exceeded as exc:
        print('  EXCEEDED -- %s' % exc)
        print('  Reading: de-cancellation does not fit this circuit in the form')
        print('  measured here, and the ladder overhead did not carry over.')
    except RecursionError:
        print('  RecursionError')

    per_coefficient(A, devices, positions, env, system.s, k_compact)
    return 0


if __name__ == '__main__':
    sys.setrecursionlimit(10000)
    sys.exit(main())
