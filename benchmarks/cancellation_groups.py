"""Stages 1-2 of `doc/cancellation_ranking_plan.md` -- ranking groups of terms.

A **group** is a pair ``(P, v)`` of a chosen prefix of matrix entries and a
diagram node.  It stands for every product term below ``v`` carrying that
prefix, and its contribution to the determinant is ``(prod P) * V[v]``
*exactly*, where ``V[v]`` is the value of the subdiagram -- one memoised pass.

That exactness is the whole point.  Magnitude ranking sees only
``|prod entries|``, an upper bound, so a million mutually cancelling terms must
all be kept or the answer is lost.  A group knows its true contribution, so the
same million are dropped in one step at their real -- tiny -- cost.

The invariant maintained throughout is

    kept_sum + (sum over frontier groups) = det, exactly

so ``|det - kept_sum| / |det|`` is the **exact error of the kept expression at
every step**, with no bound and no bookkeeping.  Termination is that number
reaching ``tol``; whatever is left in the frontier is simply not kept.

One parameter spans both stages.  ``max_minor``:

* ``0`` -- groups are split until they reach a terminal, so every kept item is
  a plain product term.  This is stage 1.
* ``k > 0`` -- a group whose node stands for a minor of dimension ``<= k`` is
  *retained whole*, contributing ``(prod P) * det(M[rows, cols])``.  The kept
  expression is then a factored form over determinants of named matrix entries.
  This is stage 2.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/cancellation_groups.py
"""

import heapq
import sys
import time

import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix, _resolve

from cancellation_profile import FREQ, environment, profile


## The operating points of the previous round (doc/hierarchical_approximation_
## plan.md), kept identical so that its recorded magnitude-ranking failure is
## reproduced here rather than quoted.
NOMINAL = {'gm_q1': 1.0e-3, 'gm_q2': 1.0e-3, 'gm_q17': 1.0e-3}
DEGRADED = {'gm_q1': 0.1e-3, 'gm_q2': 0.1e-3, 'gm_q17': 0.1e-3}

## The independent recheck lambdifies the kept expression, which costs minutes
## once it runs to thousands of terms -- and the case that has to be verified is
## the one the gate is stated at.  Above this many items it is skipped and said
## to be skipped.  (Measured: 734 terms recheck instantly, 7573 took ~4 min and
## agreed to every printed digit.)
RECHECK_CAP = 2000


def minor_positions(ddd):
    """``id(node) -> (rows, cols)``: the minor each node's subdiagram expands.

    Derived from the diagram alone, without the builder: the root is the whole
    matrix, a 0-edge stays inside the same minor (it is the next sibling of the
    same expansion), and a 1-edge strikes out the vertex's own row and column.
    Sharing is keyed on ``(rows, cols)`` by construction, so the assignment is
    consistent -- asserted rather than assumed.
    """
    n = ddd.matrix.rows
    pos = {id(ddd.root): (tuple(range(n)), tuple(range(n)))}
    stack = [ddd.root]
    while stack:
        node = stack.pop()
        if node.is_terminal:
            continue
        rows, cols = pos[id(node)]
        children = ((node.one_edge, (tuple(r for r in rows if r != node.row),
                                     tuple(c for c in cols if c != node.col))),
                    (node.zero_edge, (rows, cols)))
        for child, key in children:
            if child.is_terminal:
                continue
            seen = pos.get(id(child))
            if seen is None:
                pos[id(child)] = key
                stack.append(child)
            else:
                assert seen == key, 'node shared between distinct minors'
    return pos


def group_rank(ddd, env, tol=0.05, max_minor=0, max_splits=500000,
               values=None, pos=None):
    """Best-first expansion over groups.  Returns a result dict.

    Args:
        ddd: the diagram.
        env: numeric substitution for every symbol.
        tol: stop once the kept expression's *exact* relative error is here.
        max_minor: retain groups whose minor is at most this dimension instead
            of splitting them (0 = expand everything to product terms).
        max_splits: give up after this many splits, and say so.
        values, pos: precomputed subdiagram values and minor positions.
    """
    if values is None:
        values = profile(ddd, env)[0]
    if pos is None:
        pos = minor_positions(ddd)
    entry = {}
    for v in ddd.vertices():
        key = id(v.entry)
        if key not in entry:
            entry[key] = complex(_resolve(v.entry, env))

    exact = values[id(ddd.root)]
    if exact == 0:
        raise ValueError('determinant is zero at this operating point')

    ## Heap of non-terminal groups, largest true contribution first.  A counter
    ## breaks ties so that nodes are never compared.
    counter = 0
    heap = [(-abs(exact), counter, ddd.root, 1.0 + 0j, None)]
    kept, kept_sum, splits, dropped = [], 0.0 + 0j, 0, 0
    peak, retained_dims, minor_dets = 1, [], {}

    def term_of(chain, tail):
        factors = []
        while chain is not None:
            v, chain = chain
            factors.append(v.sign * v.entry)
        factors.append(tail)
        return sympy.Mul(*reversed(factors))

    while True:
        err = abs(exact - kept_sum) / abs(exact)
        if err <= tol or not heap or splits >= max_splits:
            break
        neg, _, node, acc, chain = heapq.heappop(heap)
        if -neg == 0.0:
            dropped += 1
            continue

        rows, cols = pos[id(node)]
        if max_minor and len(rows) <= max_minor:
            ## The same minor is retained under many different prefixes, and a
            ## symbolic 6x6 determinant is not cheap, so memoise on (rows, cols).
            key = (rows, cols)
            if key not in minor_dets:
                minor_dets[key] = ddd.matrix.extract(list(rows),
                                                     list(cols)).det()
            kept.append(term_of(chain, minor_dets[key]))
            retained_dims.append(len(rows))
            kept_sum += acc * values[id(node)]
            continue

        splits += 1
        take_acc = acc * node.sign * entry[id(node.entry)]
        for child, child_acc, child_chain in (
                (node.one_edge, take_acc, (node, chain)),
                (node.zero_edge, acc, chain)):
            contribution = child_acc * values[id(child)]
            if contribution == 0:
                dropped += 1
                continue
            if child.is_terminal:
                kept.append(term_of(child_chain, child.value))
                retained_dims.append(1)
                kept_sum += contribution
                continue
            counter += 1
            heapq.heappush(heap, (-abs(contribution), counter, child,
                                  child_acc, child_chain))
        peak = max(peak, len(heap))

    err = abs(exact - kept_sum) / abs(exact)
    return {'exact': exact, 'error': err, 'n_kept': len(kept),
            'splits': splits, 'peak_frontier': peak, 'zero_pruned': dropped,
            'converged': err <= tol, 'items': kept, 'dims': retained_dims,
            'exhausted': not heap, 'capped': splits >= max_splits}


def recheck(items, env, exact):
    """Independent value of the kept expression, via sympy rather than the
    bookkeeping that selected it.  ``None`` when the expression is too large."""
    if not items or len(items) > RECHECK_CAP:
        return None
    expr = sympy.Add(*items)
    free = sorted(expr.free_symbols, key=str)
    fn = sympy.lambdify(free, expr, modules='numpy')
    got = complex(fn(*[env[s] for s in free]))
    return abs(got - exact) / abs(exact)


def show(indent, res, env, secs, devices):
    expr = sympy.Add(*res['items'])
    free = sorted(str(s) for s in expr.free_symbols)
    dims = sorted(res['dims'])
    state = ('CONVERGED' if res['converged']
             else 'EXHAUSTED' if res['exhausted']
             else 'SPLIT CAP' if res['capped'] else 'NOT converged')
    print('%skept %-6d err %.3e  %-10s %.2f s  splits %d, peak frontier %d'
          % (indent, res['n_kept'], res['error'], state, secs, res['splits'],
             res['peak_frontier']))
    hist = {}
    for d in dims:
        hist[d] = hist.get(d, 0) + 1
    print('%s  retained minor dims: %s'
          % (indent, ', '.join('%dx%d:%d' % (d, d, n)
                               for d, n in sorted(hist.items()))))
    present = [s for s in free if s in devices]
    print('%s  device symbols: %s'
          % (indent, ', '.join(present) if present else 'NONE -- GATE FAILS'))
    check = recheck(res['items'], env, res['exact'])
    print('%s  independent sympy recheck: %s'
          % (indent, 'err %.3e' % check if check is not None else 'skipped'))
    return expr


def main():
    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q17'))
    devices = {'gm_q1', 'gm_q2', 'gm_q17', 's'}
    ddd = ddd_of_matrix(system.A)
    print('uA741: dim %d, %d vertices, %d terms, %g Hz'
          % (system.dim, ddd.size, ddd.term_count(), FREQ))
    print('operating points: nominal gm=1.0 mA/V, degraded gm=0.1 mA/V'
          '  (the previous round\'s points)')
    print()

    points = (('nominal', NOMINAL), ('degraded', DEGRADED))
    cache = {}
    for name, gms in points:
        env = environment(system, gms)
        values, absolutes, _d = profile(ddd, env)
        cache[name] = (env, values, minor_positions(ddd),
                       absolutes[id(ddd.root)] / abs(values[id(ddd.root)]))

    print('== STAGE 1: group ranking, expanded to product terms ==')
    group_counts = {}
    for name, _gms in points:
        env, values, pos, kappa = cache[name]
        print('  %s  (k[root] = %.3e)' % (name, kappa))
        for tol in (0.2, 0.05, 1e-3):
            t0 = time.time()
            res = group_rank(ddd, env, tol=tol, values=values, pos=pos)
            print('    tol=%-7g' % tol)
            show('      ', res, env, time.time() - t0, devices)
            if tol == 0.05:
                group_counts[name] = res['n_kept']
    print()

    print('== BASELINE: magnitude ranking (DDD.approximate), same diagram ==')
    print('   re-measured here, not quoted; the previous round recorded')
    print('   994% at nominal and 17 300% at degraded for 500 terms.')
    for name, _gms in points:
        env, values, _pos, kappa = cache[name]
        for max_terms in (500, group_counts[name]):
            t0 = time.time()
            _expr, n, err = ddd.approximate(env, tol=0.05, max_terms=max_terms)
            print('  %-9s max_terms=%-6d kept %-6d err %.3e  %.2f s'
                  '   err/k = %.3e'
                  % (name, max_terms, n, err, time.time() - t0, err / kappa))
    print('  err/k is the residual absolute-mass fraction.  If it agreed across')
    print('  operating points, k[root] alone would predict the error.')
    print()

    print('== STAGE 2: retain minors instead of expanding them ==')
    for max_minor in (2, 3, 6):
        print('  max_minor=%d' % max_minor)
        for name, _gms in points:
            env, values, pos, _k = cache[name]
            for tol in (0.05, 1e-3):
                t0 = time.time()
                res = group_rank(ddd, env, tol=tol, max_minor=max_minor,
                                 values=values, pos=pos)
                print('    %-9s tol=%-7g' % (name, tol))
                show('      ', res, env, time.time() - t0, devices)
    print()

    print('== the smallest expression that met tol=0.2, nominal ==')
    env, values, pos, _k = cache['nominal']
    res = group_rank(ddd, env, tol=0.2, max_minor=3, values=values, pos=pos)
    print('  %d items, err %.3e' % (res['n_kept'], res['error']))
    for item in res['items'][:12]:
        print('   +', item)
    if res['n_kept'] > 12:
        print('   ... %d more' % (res['n_kept'] - 12))
    return 0


if __name__ == '__main__':
    sys.exit(main())
