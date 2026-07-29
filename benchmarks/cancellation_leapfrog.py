"""Stage 3 of `doc/cancellation_ranking_plan.md` -- group ranking on the target.

The 5th-order leapfrog of five µA741s: 127 unknowns, 121 symbols, 536 nonzeros.
A flat diagram does not reach it (killed at 15 min / 2.7 GB).  A hierarchical
one does, by suppressing the 111 internal nodes and keeping the 16 terminals.

The level walk below is `HierarchicalDDD.eval` stopped one step early: each
level is resolved into the environment, and then instead of evaluating the top
diagram, group ranking is run *on* it.  That is the recipe the previous round
found; it is copied rather than shared because `eval` has no hook for it.

Two things are measured, and the second is the honest limitation:

1. how many groups the top diagram needs to reach a tolerance, and
2. what those groups are *named over* -- the top diagram's entries are the
   level stamp symbols (`_lvl109_16_0`), not device parameters, which is
   exactly the interpretability failure recorded as stage A of
   `doc/hierarchical_approximation_plan.md`.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/cancellation_leapfrog.py
"""

import math
import sys
import time

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import (HierarchicalDDD, suppression_order,
                                   eval_roots, _resolve)

from cancellation_groups import group_rank, minor_positions
from cancellation_profile import FREQ, profile


SYMBOLIC = tuple('s%d_q%s' % (k, q) for k in range(5) for q in ('1', '2', '17'))


def environment(system, scale):
    """Every symbol numeric, with the symbolic gm's scaled by ``scale``."""
    env = {system.s: 2j * math.pi * FREQ}
    for sym in sorted(system.A.free_symbols, key=str):
        if sym == system.s:
            continue
        value = complex(system.params[sym])
        if str(sym).startswith('gm_'):
            value = value * scale
        env[sym] = value
    return env


def resolve_levels(hier, env):
    """`HierarchicalDDD.eval`'s level walk, stopped before the top diagram.

    Returns ``(env, factors)``: the environment extended with every level's
    stamp symbols resolved to numbers, and the per-level determinant factors
    that multiply the top diagram's value.
    """
    env = dict(env)
    factors = []
    for lvl in hier.levels:
        roots = [lvl['family'].denominator.root]
        for combination in lvl['blocks'].values():
            roots.extend(combination.roots())
        memo = eval_roots(roots, env)
        D = memo[id(lvl['family'].denominator.root)]
        if D == 0:
            raise ZeroDivisionError('an internal block is singular')
        env.update({sym: combination.resolve(env, memo) / D
                    for sym, combination in lvl['blocks'].items()})
        factors.append(D)
    return env, factors


def main():
    t0 = time.time()
    system = bc.leapfrog_5th_order(symbolic_devices=SYMBOLIC)
    n = system.dim
    nnz = sum(1 for i in range(n) for j in range(n) if system.A[i, j] != 0)
    print('leapfrog: %d unknowns, %d symbols, %d nonzeros, built in %.1f s'
          % (n, len(system.A.free_symbols), nnz, time.time() - t0))

    ## Terminals: each amplifier exposes inp, inn, out; plus the filter input.
    keep_names = ['in'] + ['s%d_%s' % (k, b) for k in range(5)
                           for b in ('inp', 'inn', 'out')]
    keep = [system.cir.get_node_index(name) for name in keep_names]
    keep = [i for i in keep if i < n]
    print('keeping %d terminals: %s' % (len(keep), ', '.join(keep_names)))

    t0 = time.time()
    blocks = suppression_order(system.A, keep=keep)
    hier = HierarchicalDDD(system.A, blocks)
    print('hierarchical diagram: %d levels, top %dx%d, %d vertices, %.1f s'
          % (len(hier.levels), hier.top.matrix.rows, hier.top.matrix.rows,
             hier.size, time.time() - t0))
    print('top diagram: %d vertices, %d terms'
          % (hier.top.size, hier.top.term_count()))
    print()

    pos = minor_positions(hier.top)
    for label, scale in (('nominal', 1.0), ('degraded', 0.1)):
        base = environment(system, scale)
        t0 = time.time()
        env, factors = resolve_levels(hier, base)
        walk = time.time() - t0

        values, absolutes, _d = profile(hier.top, env)
        root = hier.top.root
        kappa = absolutes[id(root)] / abs(values[id(root)])
        top = values[id(root)]

        ## det(A) itself is unrepresentable: 111 level factors multiply into
        ## something around 1e-400, which underflows to 0.0 -- and `hier.eval`
        ## duly returns 0 here.  That is a range problem in the *product*, not a
        ## ranking problem, so carry it as a log magnitude and rank on the top
        ## diagram, whose own value is O(1) because each level was divided by D.
        log10_det = math.log10(abs(top)) + sum(math.log10(abs(D))
                                              for D in factors)
        ## Independent oracle for the quantity actually being ranked.
        reduced = np.array([[complex(_resolve(hier.reduced[i, j], env))
                             for j in range(hier.reduced.cols)]
                            for i in range(hier.reduced.rows)])
        numpy_det = np.linalg.det(reduced)
        diagram_top = complex(hier.top.eval(env))

        print('  %s gm (x%g): level walk %.2f s' % (label, scale, walk))
        print('    log10|det(A)| = %.1f   (underflows double precision;'
              ' hier.eval returns %s)' % (log10_det, complex(hier.eval(base))))
        print('    top diagram value  = %.6e%+.6ej' % (top.real, top.imag))
        print('      vs DDD.eval      rel %.2e'
              % (abs(diagram_top - top) / abs(top)))
        print('      vs numpy.linalg.det of the reduced matrix  rel %.2e'
              % (abs(numpy_det - top) / abs(top)))
        print('    k[root] of the TOP diagram = %.3e' % kappa)

        for tol in (1e-3, 1e-8, 1e-16):
            t0 = time.time()
            res = group_rank(hier.top, env, tol=tol, values=values, pos=pos)
            expr = sympy.Add(*res['items'])
            free = sorted(str(s) for s in expr.free_symbols)
            device = [s for s in free if s.startswith('gm_') or s == 's']
            print('    tol=%-7g kept %-4d err %.3e  %-10s %.2f s'
                  % (tol, res['n_kept'], res['error'],
                     'CONVERGED' if res['converged'] else 'not converged',
                     time.time() - t0))
            print('               symbols: %d, all level stamps: %s,'
                  ' device symbols: %s'
                  % (len(free), all(s.startswith('_lvl') for s in free),
                     ', '.join(device) if device else 'NONE'))
        print()

    print('what a kept group looks like (nominal, tol=1e-8):')
    env, _factors = resolve_levels(hier, environment(system, 1.0))
    values = profile(hier.top, env)[0]
    res = group_rank(hier.top, env, tol=1e-8, values=values, pos=pos)
    for item in res['items'][:6]:
        print('   +', item)
    if res['n_kept'] > 6:
        print('   ... %d more' % (res['n_kept'] - 6))
    print()
    print('and one stamp symbol resolved back to what it stands for:')
    lvl = hier.levels[-1]
    sym = sorted(lvl['blocks'], key=str)[0]
    combination = lvl['blocks'][sym]
    print('   %s = a DDDCombination of %d parts over the block cofactors;'
          % (sym, len(combination.parts)))
    print('   its own diagram stands for %d product terms of device entries.'
          % combination.term_count())
    return 0


if __name__ == '__main__':
    sys.exit(main())
