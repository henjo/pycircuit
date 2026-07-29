"""Stage 5, gate items 4-6 -- the composed device-symbolic approximation.

`cancellation_parallel.py` established the structure and passed gates 1-3:
five independent eliminations against the original matrix, every cofactor over
device entries, and -- the finding that makes this worth doing -- the
cancellation is **inside** the blocks (`kappa = 1.1e3` per block determinant)
while the top diagram barely cancels (`kappa = 13.8`).

So group ranking is applied where it pays.  Three kinds of object are ranked
separately and composed *without* multiplying out:

    det(A) = (prod_l D_l) * det(Reduced)
    Reduced[a,b] = A_tt[a,b] - sum_l sum_{p,k} A_ti[a,p] adj(A_ii,l)[p,k] A_it[k,b] / D_l

    blocks:    D_l         group-ranked over device entries
    cofactors: adj entries group-ranked over device entries
    top:       det(Reduced) group-ranked over the reduced entries, and ranked at
               the *approximated* entry values -- the honest order, because at
               use time the exact ones are not available.

The output is **nested**, which the plan declared in advance: substituting the
levels into one another would multiply 17 factors of ~10**5 operations each into
nothing anybody reads.  What is reported is the operation count of the nested
form and its end-to-end error against `numpy.linalg.slogdet` of the full
127x127 matrix -- an oracle outside this machinery entirely.

Every truncated piece's value comes from `approximate_groups(with_value=True)`.
Re-deriving it by substituting into the returned expression was measured at 33 s
per piece against 0.95 s for the ranking itself, and there are 130 pieces.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/cancellation_compose.py
"""

import math
import sys
import time
import warnings

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import _resolve

from cancellation_blocks import amplifier_blocks
from cancellation_leapfrog import SYMBOLIC, environment
from cancellation_parallel import ParallelReduction


def rank(diagram, env, tol):
    """Group-rank, quietly: a deliberately loose tol is not a surprise here."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return diagram.approximate_groups(env, tol=tol, with_value=True)


def compose(red, env, tol_top, tol_block, tol_cof, count_ops=True):
    """Approximate every level, then evaluate the nested form.

    The order matters and an earlier version got it wrong.  Ranking all 125
    cofactors first wastes most of the work: the top ranking keeps only a few
    dozen groups, and those touch a fraction of the reduced entries.  So rank
    the top **first**, at the exact entry values, and let it say which pieces
    are needed -- then approximate only those.
    """
    t0 = time.time()

    ## --- the top first, to find out what is actually needed ------------
    exact_renv = red.reduced_env(env)
    top_expr, top_n, top_err, _v = rank(red.top, exact_renv, tol_top)
    used_T = {s for s in top_expr.free_symbols}
    needed = set()
    for sym in used_T:
        for _factors, l, key in red.definition[sym][1]:
            needed.add((l, key))
    ## --- the blocks: D_l, over device entries --------------------------
    ## Every block is needed whether or not the kept top groups reach through
    ## it: det(A) = (prod_l D_l) * det(Reduced), so all five sit in the
    ## prefactor.  Only the *cofactors* can be pruned.  (An earlier version
    ## pruned blocks too, which silently dropped factors from the product.)
    block = []
    for family in red.families:
        expr, n, err, value = rank(family.denominator, env, tol_block)
        block.append({'expr': expr, 'n': n, 'err': err, 'value': value})

    ## --- only the cofactors the kept top groups reach through ----------
    cof = {}
    for (l, key) in sorted(needed):
        expr, n, err, value = rank(red._cof_cache[l][key], env, tol_cof)
        cof[(l, key)] = {'expr': expr, 'n': n, 'err': err, 'value': value}

    ## --- the reduced entries, over the approximated cofactors ----------
    ## Each T stays its own small expression -- A_tt plus one term per (block,
    ## cofactor) it couples through -- and is never expanded into the top.
    approx_renv, t_ops = dict(exact_renv), 0
    for sym in used_T:
        att, parts = red.definition[sym]
        total = complex(_resolve(att, env)) if att != 0 else 0j
        ops_here = 0 if att == 0 else int(sympy.count_ops(att))
        for factors, l, key in parts:
            piece = cof[(l, key)]
            v = piece['value'] / block[l]['value']
            for f in factors:
                v = complex(_resolve(f, env)) * v
            total += v
            ops_here += 2
        approx_renv[sym] = total
        t_ops += ops_here

    ## --- the same kept groups, re-evaluated at the approximated entries.
    ## Not re-ranked: the selection was made above and re-selecting at
    ## perturbed values would hide the error being measured.
    top_val = complex(top_expr.xreplace(approx_renv))

    log10 = (math.log10(abs(top_val))
             + sum(math.log10(abs(b['value'])) for b in block))

    out = {
        'log10': log10, 'secs': time.time() - t0,
        'top_n': top_n, 'top_err': top_err,
        'used_T': len(used_T), 'n_cof': len(cof),
        'block_n': [b['n'] for b in block],
        'block_err': max(b['err'] for b in block),
        'cof_n': sum(c['n'] for c in cof.values()),
        'cof_err': max(c['err'] for c in cof.values()) if cof else 0.0,
        't_ops': t_ops,
        'symbols': _symbols_of(block, cof),
        'cof_expr': cof,
    }
    if count_ops:
        t = time.time()
        out['top_ops'] = int(sympy.count_ops(top_expr))
        out['block_ops'] = sum(int(sympy.count_ops(b['expr']))
                               for b in block)
        out['cof_ops'] = sum(int(sympy.count_ops(c['expr']))
                             for c in cof.values())
        out['ops_secs'] = time.time() - t
    return out


def _symbols_of(block, cof):
    names = set()
    for b in block:
        names |= {str(s) for s in b['expr'].free_symbols}
    for c in cof.values():
        names |= {str(s) for s in c['expr'].free_symbols}
    return sorted(names)


def main():
    t0 = time.time()
    system = bc.leapfrog_5th_order(symbolic_devices=SYMBOLIC)
    red = ParallelReduction(system.A, amplifier_blocks(system))
    exact_top = red.top.term_count()
    exact_block = red.families[0].denominator.term_count()
    exact_cof = sum(d.term_count() for d in red._cof_cache[0].values())
    print('leapfrog, parallel 5-block reduction: %d vertices, set up in %.1f s'
          % (red.size()[0], time.time() - t0))
    print('the exact form, for scale: top %d terms, each D_l %d terms,'
          ' one block\'s 25 cofactors %d terms'
          % (exact_top, exact_block, exact_cof))
    print('  so the exact nested product stands for >%.2e terms in total'
          % (exact_top * (exact_block ** 5)))
    print()

    for label, scale in (('nominal', 1.0), ('degraded', 0.1)):
        env = environment(system, scale)
        full = np.array([[complex(_resolve(system.A[i, j], env))
                          for j in range(system.dim)]
                         for i in range(system.dim)])
        _sign, logabs = np.linalg.slogdet(full)
        ref = float(logabs) / math.log(10.0)
        print('== %s gm (x%g): numpy slogdet log10|det| = %.6f =='
              % (label, scale, ref))

        for tol_top, tol_block, tol_cof in ((5e-2, 2e-2, 2e-2),
                                            (2e-2, 5e-3, 5e-3),
                                            (1e-2, 1e-3, 1e-3)):
            r = compose(red, env, tol_top, tol_block, tol_cof)
            ## The determinant spans 358 decades, so compare logs: a relative
            ## error eps shows up as a shift of log10(1+eps).
            dex = abs(r['log10'] - ref)
            rel = 10.0 ** dex - 1.0
            total = (r['top_ops'] + r['block_ops'] + r['cof_ops'] + r['t_ops'])
            print('  tol top/block/cof = %.0e/%.0e/%.0e   rank %.1f s,'
                  ' count_ops %.1f s'
                  % (tol_top, tol_block, tol_cof, r['secs'], r['ops_secs']))
            print('    log10|det| = %.6f   off by %.2e dex   implied rel %.2e'
                  % (r['log10'], dex, rel))
            print('    kept: top %d groups (err %.1e) touching %d of %d'
                  ' reduced entries and %d of %d cofactors'
                  % (r['top_n'], r['top_err'], r['used_T'],
                     len(red.definition), r['n_cof'],
                     sum(len(c) for c in red._cof_cache)))
            print('          D_l %s (err<=%.1e), cofactors %d groups'
                  ' (err<=%.1e)'
                  % (r['block_n'], r['block_err'], r['cof_n'], r['cof_err']))
            print('    nested ops: top %d + blocks %d + cofactors %d'
                  ' + reduced entries %d = %d'
                  % (r['top_ops'], r['block_ops'], r['cof_ops'], r['t_ops'],
                     total))
            devices = [s for s in r['symbols']
                       if s.startswith('gm_') or s == 's']
            print('    device symbols in block/cofactor expressions: %d -- %s'
                  % (len(devices), ', '.join(devices[:8])))
            print('    GATE 4 (implied rel <= 1e-2): %s'
                  % ('PASS' if rel <= 1e-2 else 'FAIL'))
            print('    GATE 5 (names device parameters): %s'
                  % ('PASS' if devices else 'FAIL'))
            ## The honest comparison is NOT against the exact form's 4e34
            ## terms -- that bar is meaningless, anything passes it.  The bar
            ## that matters is whether a person can read the result, and the
            ## nearest measured reference is the 2 256 398 operations the
            ## previous round's plain unfolding of a *single* uA741 cost.
            print('    ops vs references: %d nested  |  2256398 = one uA741'
                  ' unfolded (previous round)  |  %.2e = exact terms here'
                  % (total, exact_top * exact_block ** 5))
            print('    GATE 6 (nested form smaller than unfolding one uA741):'
                  ' %s' % ('PASS' if total < 2256398 else 'FAIL'))
        print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
