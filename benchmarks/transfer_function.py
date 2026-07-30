"""Stage 11 -- the transfer function, in operations, verified over a sweep.

Two parts, together because both are operation counts and stating either honestly
needs the other.

**Part 1** re-reads verdicts that were reached on *term counts* before §23 showed
those overstate the answer 4-8x: the whole determinant's ranking, and the leapfrog's
181 groups from stage 3.

**Part 2** is the deliverable.  Everything measured so far approximates ``det(A)``,
the *denominator*.  A designer wants ``H(s) = N(s)/D(s)``, and both are determinants
-- the numerator being ``A`` with the output column replaced by ``b`` -- so
``s_expand`` applies to each.

Two methodological requirements, both from lessons this plan recorded the hard way:

* **verify over a frequency sweep**, not at one point.  Every measurement here bar
  stage 6 has been at 1 kHz, and single-frequency error control guarantees nothing
  between samples.  Still sampled, and said to be.
* **choose the coefficients greedily with exact re-evaluation of H after each drop**,
  never by a per-coefficient budget.  That is the global-reference lesson of Guerra
  et al. and the non-monotonicity of stage 5; stage 6 showed exact re-evaluation
  removes it.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/transfer_function.py
"""

import math
import sys
import time
import warnings

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import (HierarchicalDDD, ddd_of_matrix, s_expand,
                                   suppression_order, _resolve)

from cancellation_profile import NOMINAL, DEGRADED, environment


SWEEP = np.logspace(2, 4, 13)          # two decades, 100 Hz .. 100 kHz
DEVICES = {'gm_q1', 'gm_q2', 'gm_q17'}


def ops(expr):
    return 0 if expr == 0 else int(sympy.count_ops(expr))


def rank(diagram, env, tol):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return diagram.approximate_groups(env, tol=tol)


## ---------------------------------------------------------------- part 1 --

def reread_old_verdicts():
    print('== PART 1: verdicts reached on term counts, re-read in operations ==')
    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q17'))
    D = ddd_of_matrix(system.A)
    env = environment(system, NOMINAL)
    print('  uA741 whole determinant (%d terms):' % D.term_count())
    for tol in (0.5, 0.2, 0.05):
        expr, n, err = rank(D, env, tol)
        print('    tol=%-6g %6d groups -> %5d operations, err %.3e'
              % (tol, n, ops(expr), err))

    ## The leapfrog, whose 181-group result was called uninterpretable -- it still
    ## is, but the size verdict was term-count based.
    t0 = time.time()
    sym = tuple('s%d_q%s' % (k, q) for k in range(5) for q in ('1', '2', '17'))
    filt = bc.leapfrog_5th_order(symbolic_devices=sym)
    keep_names = ['in'] + ['s%d_%s' % (k, b) for k in range(5)
                           for b in ('inp', 'inn', 'out')]
    keep = [filt.cir.get_node_index(x) for x in keep_names]
    keep = [i for i in keep if i < filt.dim]
    hier = HierarchicalDDD(filt.A, suppression_order(filt.A, keep=keep))
    base = {filt.s: 2j * math.pi * 1e3}
    for s in filt.A.free_symbols:
        if s is not filt.s:
            base[s] = complex(filt.params[s])
    ## The level walk, resolved into an environment for the top diagram.
    from cancellation_leapfrog import resolve_levels
    tenv, _factors = resolve_levels(hier, base)
    print('  leapfrog top diagram (%d terms), built in %.1f s:'
          % (hier.top.term_count(), time.time() - t0))
    for tol in (1e-2, 1e-3):
        expr, n, err = rank(hier.top, tenv, tol)
        names = sorted(str(x) for x in expr.free_symbols)
        print('    tol=%-6g %6d groups -> %5d operations, err %.3e,'
              ' all placeholders: %s'
              % (tol, n, ops(expr), err,
                 all(x.startswith('_lvl') for x in names)))
    print()


## ---------------------------------------------------------------- part 2 --

def numerator_matrix(system):
    """``A`` with the output column replaced by ``b``: Cramer's numerator."""
    A = sympy.Matrix(system.A)
    b = sympy.Matrix(system.b)
    vin = [s for s in b.free_symbols if str(s) == 'vin']
    if vin:
        b = b.subs({vin[0]: 1})
    for i in range(A.rows):
        A[i, system.out_index] = b[i]
    return A


def coefficients(expanded, env):
    """``{k: (diagram, exact numeric value)}`` for the non-empty coefficients."""
    out = {}
    for k in range(expanded.degree + 1):
        c = expanded.coefficient(k)
        if c.root is None or c.term_count() == 0:
            continue
        out[k] = (c, complex(c.eval(env)))
    return out


def h_from(num_vals, den_vals, freq):
    s = 2j * math.pi * freq
    n = sum(v * s ** k for k, v in num_vals.items())
    d = sum(v * s ** k for k, v in den_vals.items())
    return None if d == 0 else n / d


def sweep_error(num_vals, den_vals, reference):
    worst = 0.0
    for freq, ref in zip(SWEEP, reference):
        got = h_from(num_vals, den_vals, freq)
        if got is None:
            return math.inf
        worst = max(worst, abs(got - ref) / abs(ref))
    return worst


def choose_coefficients(num, den, reference, tol):
    """Drop coefficients greedily, re-evaluating H over the whole sweep each time.

    Never a per-coefficient budget -- that is what stage 5 measured going
    non-monotone.  A drop is accepted only if the *global* sweep error still holds.
    """
    nv = {k: v for k, (_c, v) in num.items()}
    dv = {k: v for k, (_c, v) in den.items()}
    changed = True
    while changed:
        changed = False
        for which, table in (('n', nv), ('d', dv)):
            for k in sorted(table, key=lambda x: abs(table[x])):
                if len(table) == 1:
                    break
                saved = table.pop(k)
                if sweep_error(nv, dv, reference) <= tol:
                    changed = True
                    break
                table[k] = saved
    return sorted(nv), sorted(dv)


def transfer_function():
    print('== PART 2: the transfer function H(s) = N(s)/D(s) ==')
    print('   sweep: %d points, %g Hz to %g Hz (SAMPLED)'
          % (len(SWEEP), SWEEP[0], SWEEP[-1]))
    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q17'))
    Anum = numerator_matrix(system)

    t0 = time.time()
    den_exp = s_expand(system.A, system.s)
    num_exp = s_expand(Anum, system.s)
    print('   s-expanded: denominator degree %d, numerator degree %d  (%.1f s)'
          % (den_exp.degree, num_exp.degree, time.time() - t0))
    print()

    for label, gms in (('nominal', NOMINAL), ('degraded', DEGRADED)):
        full = environment(system, gms)
        cenv = {k: v for k, v in full.items() if k != system.s}

        ## The oracle: solve the real system at each frequency.
        reference = []
        for freq in SWEEP:
            env = dict(cenv)
            env[system.s] = 2j * math.pi * freq
            M = np.array([[complex(_resolve(system.A[i, j], env))
                           for j in range(system.dim)]
                          for i in range(system.dim)], dtype=complex)
            rhs = np.array([complex(_resolve(system.b[i], {**env, **{
                s: 1 for s in system.b.free_symbols if str(s) == 'vin'}}))
                for i in range(system.dim)], dtype=complex)
            reference.append(np.linalg.solve(M, rhs)[system.out_index])
        reference = np.array(reference)

        num = coefficients(num_exp, cenv)
        den = coefficients(den_exp, cenv)
        exact_err = sweep_error({k: v for k, (_c, v) in num.items()},
                                {k: v for k, (_c, v) in den.items()}, reference)
        print('  == %s gm ==' % label)
        print('    exact N/D against the solve, over the sweep: %.2e' % exact_err)

        for tol in (0.2, 0.05):
            keep_n, keep_d = choose_coefficients(num, den, reference, tol)
            ## Now approximate each surviving coefficient symbolically, at the
            ## same tolerance, and re-check the sweep with the approximated values.
            nsym, nval, nops = {}, {}, 0
            for k in keep_n:
                expr, _n, _e = rank(num[k][0], cenv, tol)
                nsym[k] = expr
                nval[k] = complex(expr.xreplace(cenv)) if expr != 0 else 0j
                nops += ops(expr)
            dsym, dval, dops = {}, {}, 0
            for k in keep_d:
                expr, _n, _e = rank(den[k][0], cenv, tol)
                dsym[k] = expr
                dval[k] = complex(expr.xreplace(cenv)) if expr != 0 else 0j
                dops += ops(expr)
            err = sweep_error(nval, dval, reference)
            names = set()
            for e in list(nsym.values()) + list(dsym.values()):
                names |= {str(x) for x in e.free_symbols}
            present = sorted(names & DEVICES)
            total = nops + dops + 2 * (len(keep_n) + len(keep_d))
            print('    tol=%-6g N keeps %s, D keeps %s' % (tol, keep_n, keep_d))
            print('             operations: N %d + D %d + assembly = %d,'
                  ' sweep err %.3e' % (nops, dops, total, err))
            print('             device symbols: %s'
                  % (', '.join(present) if present else 'NONE'))
            print('             GATE 11-2 (<200 ops at <=20%% over the sweep): %s'
                  % ('PASS' if total < 200 and err <= 0.20 and present
                     else 'FAIL'))
            if tol == 0.2 and label == 'nominal':
                print()
                print('    H(s) at tol=0.2, so it can be read:')
                for k in sorted(dsym, reverse=True):
                    print('      D s^%-2d = %s' % (k, dsym[k]))
                for k in sorted(nsym, reverse=True):
                    print('      N s^%-2d = %s' % (k, nsym[k]))
                print()
    return 0


def main():
    reread_old_verdicts()
    transfer_function()
    return 0


if __name__ == '__main__':
    sys.exit(main())
