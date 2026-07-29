"""Stage 5, second attempt -- parallel block elimination, device symbols kept.

The first attempt (`cancellation_blocks.py`) failed gate 2:
`HierarchicalDDD._suppress` renames *every* nonzero of the reduced matrix, so
suppressing amplifier 0 launders amplifier 1's device entries into `_lvl0_*`
stamps even though the elimination never touches them.  Reasoning and the
refutation: `doc/cancellation_ranking_conclusions.md` sections 9-10.

The blocks are provably independent -- measured 0 entries coupling one
amplifier's internal nodes to another's -- so ``A_ii`` over the union of internal
nodes is block diagonal and the Schur complement factorises:

    det(A) = (prod_l D_l) * det( A_tt - sum_l A_ti,l * adj(A_ii,l)/D_l * A_it,l )

Every ``A_ii,l`` is taken from ``A`` itself, so **every cofactor is over device
entries by construction** and there is no sequence to launder them.  That is the
whole point of doing it this way rather than through `HierarchicalDDD`.

Three kinds of object then get group-ranked independently, each over device
parameters: the block determinants ``D_l``, the cofactors that make up each
reduced entry, and the top diagram.  The result is a *nested* expression --
declared in advance in the plan, because substituting the levels into each other
would multiply out to nothing readable.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/cancellation_parallel.py
"""

import math
import sys
import time

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import DDDFamily, ddd_of_matrix, _resolve

from cancellation_blocks import amplifier_blocks
from cancellation_leapfrog import SYMBOLIC, environment
from cancellation_profile import FREQ


class ParallelReduction:
    """Five independent Schur eliminations against the original matrix.

    Attributes:
        families: one `DDDFamily` per block, each over device entries.
        reduced: the 17x17 matrix of fresh ``T_a_b`` symbols.
        definition: ``T_a_b -> (A_tt entry, [(factors, block, cofactor DDD)])``,
            which is what makes a reduced entry *readable* -- it names device
            entries and cofactors of the original matrix, not a reduction.
    """

    def __init__(self, A, blocks, order='min-degree'):
        self.A = A
        self.blocks = [list(b) for b in blocks]
        internal = sorted(i for b in self.blocks for i in b)
        self.rest = [i for i in range(A.rows) if i not in set(internal)]

        ## The claim the whole construction rests on, checked here rather than
        ## trusted: no entry couples two different blocks' internal nodes.
        for l, bl in enumerate(self.blocks):
            for m, bm in enumerate(self.blocks):
                if l == m:
                    continue
                for i in bl:
                    for j in bm:
                        if A[i, j] != 0:
                            raise ValueError(
                                'blocks %d and %d couple at (%d,%d); A_ii is '
                                'not block diagonal and this factorisation '
                                'does not hold' % (l, m, i, j))

        self.families, self.couplings = [], []
        for block in self.blocks:
            Aii = A.extract(block, block)
            for r in range(Aii.rows):
                for c in range(Aii.cols):
                    if any(str(x).startswith('_lvl')
                           for x in Aii[r, c].free_symbols):
                        raise ValueError('A_ii carries a stamp symbol')
            self.families.append(
                DDDFamily(Aii, sympy.zeros(len(block), 1), order=order))
            self.couplings.append((A.extract(self.rest, block),
                                   A.extract(block, self.rest)))

        Att = A.extract(self.rest, self.rest)
        m = len(self.rest)
        self.reduced = sympy.zeros(m, m)
        self.definition = {}
        self._cof_cache = [{} for _ in self.blocks]

        for a in range(m):
            for b in range(m):
                parts = []
                for l, (Ati, Ait) in enumerate(self.couplings):
                    ni = len(self.blocks[l])
                    for p in range(ni):
                        if Ati[a, p] == 0:
                            continue
                        for k in range(ni):
                            if Ait[k, b] == 0:
                                continue
                            cache = self._cof_cache[l]
                            if (k, p) not in cache:
                                cache[(k, p)] = self.families[l].cofactor(k, p)
                            ## adj(A_ii)[p,k] = (-1)**(k+p) * minor(k,p), and the
                            ## Schur complement subtracts -- same convention as
                            ## HierarchicalDDD._suppress, deliberately copied so
                            ## the signs cannot drift apart.
                            sign = -1 if (k + p) % 2 else 1
                            ## Store the cache *key*, not the diagram: a
                            ## consumer needs to look the cofactor up by name to
                            ## pair it with an approximation of it.
                            parts.append(((-sign, Ati[a, p], Ait[k, b]),
                                          l, (k, p)))
                if Att[a, b] == 0 and not parts:
                    continue
                sym = sympy.Symbol('T_%d_%d' % (a, b))
                self.definition[sym] = (Att[a, b], parts)
                self.reduced[a, b] = sym

        self.top = ddd_of_matrix(self.reduced, order='min-degree')

    def block_determinants(self, env):
        """``D_l`` at ``env``, over device entries."""
        return [complex(f.denominator.eval(env)) for f in self.families]

    def reduced_env(self, env, dets=None):
        """Numeric value of every ``T_a_b``, from the device environment."""
        dets = dets or self.block_determinants(env)
        ## One memo pass over every cofactor a reduced entry can need, rather
        ## than one walk per entry -- they overlap almost completely.
        memo = {}
        for l, cache in enumerate(self._cof_cache):
            for diagram in cache.values():
                key = id(diagram.root)
                if key not in memo:
                    memo[key] = complex(diagram.eval(env))
        out = {}
        for sym, (att, parts) in self.definition.items():
            total = complex(_resolve(att, env)) if att != 0 else 0j
            for factors, l, key in parts:
                value = memo[id(self._cof_cache[l][key].root)] / dets[l]
                for f in factors:
                    value = complex(_resolve(f, env)) * value
                total += value
            out[sym] = total
        return out

    def log10_det(self, env):
        """``log10|det(A)|``, which underflows if formed as a product."""
        dets = self.block_determinants(env)
        top = complex(self.top.eval(self.reduced_env(env, dets)))
        return (math.log10(abs(top))
                + sum(math.log10(abs(D)) for D in dets)), top

    def size(self):
        total = sum(f.size for f in self.families) + self.top.size
        cofactors = sum(len(c) for c in self._cof_cache)
        return total, cofactors


def ops(expr):
    return 0 if expr == 0 else int(sympy.count_ops(expr))


def main():
    t0 = time.time()
    system = bc.leapfrog_5th_order(symbolic_devices=SYMBOLIC)
    print('leapfrog: %d unknowns, %d symbols, built in %.1f s'
          % (system.dim, len(system.A.free_symbols), time.time() - t0))

    blocks = amplifier_blocks(system)
    t0 = time.time()
    red = ParallelReduction(system.A, blocks)
    build = time.time() - t0
    size, ncof = red.size()
    print('parallel reduction: %d blocks x %d internal, %d terminals left'
          % (len(blocks), len(blocks[0]), len(red.rest)))
    print('  built in %.1f s: %d vertices total, %d distinct cofactors'
          % (build, size, ncof))
    print('  top diagram %dx%d: %d vertices, %d terms'
          % (red.top.matrix.rows, red.top.matrix.rows, red.top.size,
             red.top.term_count()))
    print('  GATE 1a (under 5 min): %s' % ('PASS' if build < 300 else 'FAIL'))
    print('  GATE 2 (no stamp symbol in any block, no cross-block coupling):'
          ' PASS -- both are constructor assertions')
    print()

    for label, scale in (('nominal', 1.0), ('degraded', 0.1)):
        env = environment(system, scale)
        t0 = time.time()
        log10, top = red.log10_det(env)
        secs = time.time() - t0

        ## Oracle: the same reduced matrix, resolved and handed to numpy.
        renv = red.reduced_env(env)
        M = np.array([[renv.get(red.reduced[i, j], 0j)
                       for j in range(red.reduced.cols)]
                      for i in range(red.reduced.rows)])
        rel = abs(np.linalg.det(M) - top) / abs(top)

        ## Independent oracle for the *whole* determinant, which the product
        ## cannot represent: numpy's slogdet on the full 127x127 matrix.
        full = np.array([[complex(_resolve(system.A[i, j], env))
                          for j in range(system.dim)]
                         for i in range(system.dim)])
        sign, logabs = np.linalg.slogdet(full)
        ref_log10 = float(logabs) / math.log(10.0)

        print('== %s gm (x%g) ==' % (label, scale))
        print('  log10|det(A)| = %.4f   vs numpy slogdet %.4f   diff %.2e'
              % (log10, ref_log10, abs(log10 - ref_log10)))
        print('  top value %.6e%+.6ej  vs numpy.linalg.det  rel %.2e  (%.2f s)'
              % (top.real, top.imag, rel, secs))
        print('  GATE 1b (top rel < 1e-10): %s'
              % ('PASS' if rel < 1e-10 else 'FAIL'))
        print('  GATE 1c (log10 det within 1e-6 of slogdet): %s'
              % ('PASS' if abs(log10 - ref_log10) < 1e-6 else 'FAIL'))

        print('  GATE 3: kappa per object, all over device entries')
        print('    %-26s %-9s %-12s %s'
              % ('object', 'vertices', 'terms', 'kappa'))
        k_top = red.top.cancellation(renv)
        print('    %-26s %-9d %-12d %.3e'
              % ('top diagram (T symbols)', red.top.size,
                 red.top.term_count(), k_top))
        for l, family in enumerate(red.families):
            den = family.denominator
            print('    %-26s %-9d %-12d %.3e'
                  % ('D_%d = det(A_ii,%d)' % (l, l), den.size,
                     den.term_count(), den.cancellation(env)))
        worst = 0.0
        for l, cache in enumerate(red._cof_cache):
            ks = [d.cancellation(env) for d in cache.values()]
            ks = [k for k in ks if math.isfinite(k)]
            if ks:
                worst = max(worst, max(ks))
                print('    %-26s %-9s %-12d median %.3e  max %.3e'
                      % ('block %d cofactors (%d)' % (l, len(cache)), '-',
                         sum(d.term_count() for d in cache.values()),
                         sorted(ks)[len(ks) // 2], max(ks)))
        print()

    print('gate items 4-6 (composition, device symbols, op count) follow in')
    print('cancellation_compose.py once these pass.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
