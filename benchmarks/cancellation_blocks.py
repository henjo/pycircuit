"""Stage 5 of `doc/cancellation_ranking_plan.md` -- device parameters back in.

Stage 3 ranked the leapfrog's *top* diagram and got compact expressions over
level stamps (`_lvl110_*`), naming no device.  The cause is the partition, not
the ranking: `suppression_order` chose 111 single-node levels, so level k's
entries are level k-1's stamps and reaching a device parameter means unwinding
all 111 -- the unfolding hierarchy exists to avoid.  Reasoning:
`doc/cancellation_ranking_conclusions.md` section 9.

Here each amplifier is suppressed as **one block** of 22 internal nodes.  Because
amplifier 1's internal nodes couple to amplifier 0's only through terminals,
block 0's elimination stamps only into terminal rows and columns, so **every
block's internal matrix stays over device entries**.  That claim is the stage's
foundation and is asserted, not assumed.

This script runs gate items 1-3 only: build, verify, and report kappa per block.
Composition is a separate step and is not attempted until these pass.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/cancellation_blocks.py
"""

import math
import sys
import time

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import HierarchicalDDD, eval_roots, _resolve

from cancellation_leapfrog import SYMBOLIC, environment, resolve_levels
from cancellation_profile import FREQ


## Each amplifier exposes inp, inn and out; everything else in
## `_UA741_NODE_NAMES` is internal to it.
TERMINALS = ('inp', 'inn', 'out')


def amplifier_blocks(system, stages=5):
    """One block of internal-node indices per amplifier, in stage order."""
    internal = [n for n in bc._UA741_NODE_NAMES if n not in TERMINALS]
    blocks = []
    for k in range(stages):
        names = ['s%d_%s' % (k, base) for base in internal]
        blocks.append([system.cir.get_node_index(n) for n in names])
    return blocks


def assert_blocks_are_device_level(hier):
    """Gate item 2: no block's ``A_ii`` may carry a stamp symbol.

    If a later block's internal matrix picked up a `_lvl` symbol, its cofactors
    would be over another block's reduction rather than over devices, and the
    whole point of this partition would be gone.  Checked rather than argued.
    """
    report = []
    for i, lvl in enumerate(hier.levels):
        symbols = set()
        Aii = lvl['Aii']
        for r in range(Aii.rows):
            for c in range(Aii.cols):
                symbols |= Aii[r, c].free_symbols
        stamps = sorted(str(x) for x in symbols if str(x).startswith('_lvl'))
        devices = sorted(str(x) for x in symbols if str(x).startswith('gm_'))
        report.append((i, Aii.rows, len(symbols), stamps, devices))
        assert not stamps, (
            'block %d picked up stamp symbols %s -- the partition does not give '
            'device-level cofactors and stage 5 cannot proceed' % (i, stamps[:5]))
    return report


def kappa_of(diagram, env):
    """kappa for a diagram, guarding the degenerate cases."""
    values, absolutes = diagram.subdiagram_values(env)
    total, mass = values[id(diagram.root)], absolutes[id(diagram.root)]
    if mass == 0.0:
        return 0.0, complex(total), 1.0
    if total == 0:
        return mass, complex(total), math.inf
    return mass, complex(total), mass / abs(total)


def main():
    t0 = time.time()
    system = bc.leapfrog_5th_order(symbolic_devices=SYMBOLIC)
    n = system.dim
    print('leapfrog: %d unknowns, %d symbols, built in %.1f s'
          % (n, len(system.A.free_symbols), time.time() - t0))

    blocks = amplifier_blocks(system)
    print('blocks: %d amplifiers x %d internal nodes = %d suppressed, %d left'
          % (len(blocks), len(blocks[0]),
             sum(len(b) for b in blocks), n - sum(len(b) for b in blocks)))

    t0 = time.time()
    hier = HierarchicalDDD(system.A, blocks)
    build = time.time() - t0
    print('hierarchy: %d levels, top %dx%d, %d vertices total, built in %.1f s'
          % (len(hier.levels), hier.top.matrix.rows, hier.top.matrix.rows,
             hier.size, build))
    print('top diagram: %d vertices, %d terms'
          % (hier.top.size, hier.top.term_count()))
    print('GATE 1a (build under 5 min): %s  (%.1f s)'
          % ('PASS' if build < 300 else 'FAIL', build))
    print()

    print('== GATE 2: every block internal matrix is over device entries ==')
    for i, dim, nsym, stamps, devices in assert_blocks_are_device_level(hier):
        print('  block %d: A_ii is %dx%d, %d free symbols, %d stamp symbols,'
              ' devices %s' % (i, dim, dim, nsym, len(stamps),
                               ', '.join(devices) if devices else 'none'))
    print('  GATE 2: PASS (assertion held for every block)')
    print()

    for label, scale in (('nominal', 1.0), ('degraded', 0.1)):
        base = environment(system, scale)
        t0 = time.time()
        env, factors = resolve_levels(hier, base)
        walk = time.time() - t0
        top = complex(hier.top.eval(env))

        reduced = np.array([[complex(_resolve(hier.reduced[i, j], env))
                             for j in range(hier.reduced.cols)]
                            for i in range(hier.reduced.rows)])
        numpy_det = np.linalg.det(reduced)
        rel = abs(numpy_det - top) / abs(top)
        log10_det = math.log10(abs(top)) + sum(math.log10(abs(D))
                                               for D in factors)

        print('== %s gm (x%g): level walk %.2f s ==' % (label, scale, walk))
        print('  log10|det(A)| = %.1f  (product underflows; carried as a log)'
              % log10_det)
        print('  top value %.6e%+.6ej  vs numpy.linalg.det  rel %.2e'
              % (top.real, top.imag, rel))
        print('  GATE 1b (rel < 1e-10): %s' % ('PASS' if rel < 1e-10 else 'FAIL'))

        print('  GATE 3: kappa per object')
        print('    %-22s %-9s %-12s %s' % ('object', 'vertices', 'terms',
                                           'kappa'))
        _m, _t, k_top = kappa_of(hier.top, env)
        print('    %-22s %-9d %-12d %.3e'
              % ('top diagram', hier.top.size, hier.top.term_count(), k_top))
        for i, lvl in enumerate(hier.levels):
            den = lvl['family'].denominator
            _m, _t, k = kappa_of(den, base)
            print('    %-22s %-9d %-12d %.3e'
                  % ('block %d determinant' % i, den.size, den.term_count(), k))
        ## One stamp per block, the first by name, as a representative cofactor
        ## combination.  These are what the top diagram's entries actually are.
        for i, lvl in enumerate(hier.levels):
            sym = sorted(lvl['blocks'], key=str)[0]
            parts = lvl['blocks'][sym].parts
            worst = 0.0
            terms = 0
            for _factors, diagram in parts:
                _m, _t, k = kappa_of(diagram, base)
                worst = max(worst, k if math.isfinite(k) else worst)
                terms += diagram.term_count()
            print('    %-22s %-9s %-12d %.3e'
                  % ('%s (%d parts)' % (sym, len(parts)), '-', terms, worst))
        print()

    print('Interpretation is in the plan; this script only measures.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
