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


## Stage 11 used only the narrow band, which sits BELOW the uA741's compensation
## pole and unity-gain frequency, so the coefficients shaping the rolloff were never
## exercised.  Stage 12 adds the wide one; a symbolic transfer function that holds
## only below the dominant pole is not one a designer would use.
SWEEPS = (('narrow 100 Hz-10 kHz', np.logspace(2, 4, 13)),
          ('wide 10 Hz-10 MHz', np.logspace(1, 7, 25)))
SWEEP = SWEEPS[0][1]                   # rebound per band by `transfer_function`
DEVICES = {'gm_q1', 'gm_q2', 'gm_q17'}


def ops(expr):
    return 0 if expr == 0 else int(sympy.count_ops(expr))


## Coefficients that fail their tolerance, recorded by `rank` for the caller to
## report.  This exists because blanket-suppressing the warning cost a whole
## measurement: at a tight tolerance three coefficients hit `max_splits` and came
## back with 58-87% error instead of 1e-4, the composed error went to 145%, and the
## RuntimeWarning that says exactly that was being thrown away.
UNCONVERGED = []


def rank(diagram, env, tol):
    """Group-rank, recording rather than discarding a missed tolerance.

    `approximate_groups` warns when it returns above `tol`.  Suppressing that
    warning -- which an earlier version of this helper did with a bare
    ``simplefilter('ignore')`` -- defeats the mechanism added precisely so an
    unconverged approximation cannot be mistaken for a converged one.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', RuntimeWarning)
        expr, n, err = diagram.approximate_groups(env, tol=tol)
    if err > tol:
        UNCONVERGED.append((tol, n, err, len(caught)))
    return expr, n, err


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


def show_response(system, cenv, sweep):
    """Print |H| across the band, so the claimed band and its poles are visible."""
    print('    |H| across the band:')
    prev = None
    for freq in sweep[::4]:
        env = dict(cenv)
        env[system.s] = 2j * math.pi * freq
        M = np.array([[complex(_resolve(system.A[i, j], env))
                       for j in range(system.dim)] for i in range(system.dim)],
                     dtype=complex)
        rhs = np.array([complex(_resolve(system.b[i], {**env, **{
            s: 1 for s in system.b.free_symbols if str(s) == 'vin'}}))
            for i in range(system.dim)], dtype=complex)
        h = abs(np.linalg.solve(M, rhs)[system.out_index])
        slope = ''
        if prev is not None:
            f0, h0 = prev
            if h > 0 and h0 > 0:
                slope = '  %+.0f dB/dec' % (20 * math.log10(h / h0)
                                            / math.log10(freq / f0))
        print('      %10.3g Hz  |H| %11.4e  (%+6.1f dB)%s'
              % (freq, h, 20 * math.log10(h) if h > 0 else float('-inf'), slope))
        prev = (freq, h)


def h_from(num_vals, den_vals, freq):
    s = 2j * math.pi * freq
    n = sum(v * s ** k for k, v in num_vals.items())
    d = sum(v * s ** k for k, v in den_vals.items())
    return None if d == 0 else n / d


def sweep_error(num_vals, den_vals, reference, sweep=None):
    sweep = SWEEP if sweep is None else sweep
    worst = 0.0
    for freq, ref in zip(sweep, reference):
        got = h_from(num_vals, den_vals, freq)
        if got is None:
            return math.inf
        worst = max(worst, abs(got - ref) / abs(ref))
    return worst


def choose_coefficients(num, den, reference, tol, sweep=None):
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
                if sweep_error(nv, dv, reference, sweep) <= tol:
                    changed = True
                    break
                table[k] = saved
    return sorted(nv), sorted(dv)


def transfer_function():
    print('== PART 2: the transfer function H(s) = N(s)/D(s) ==')
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
        num = coefficients(num_exp, cenv)
        den = coefficients(den_exp, cenv)
        print('  ================ %s gm ================' % label)
        show_response(system, cenv, SWEEPS[1][1])

        for band, sweep in SWEEPS:
            reference = []
            for freq in sweep:
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
            exact = sweep_error({k: v for k, (_c, v) in num.items()},
                                {k: v for k, (_c, v) in den.items()},
                                reference, sweep)
            print('    -- %s (%d points): exact N/D vs the solve %.2e --'
                  % (band, len(sweep), exact))

            ## The subset is chosen against the GLOBAL sweep error, but each kept
            ## coefficient was then approximated at its own tolerance -- a
            ## per-piece budget, which is exactly what stage 5 measured as
            ## unsound and what blew the wide band to 100% error.  So the
            ## coefficient tolerance is now swept separately from the subset
            ## tolerance, and the composed error is measured for each.
            for tol in (0.2, 0.05):
                keep_n, keep_d = choose_coefficients(num, den, reference, tol,
                                                     sweep)
                for ctol in (tol, tol / 20.0, tol / 400.0):
                    nsym, nval, nops = {}, {}, 0
                    for k in keep_n:
                        expr, _n, _e = rank(num[k][0], cenv, ctol)
                        nsym[k] = expr
                        nval[k] = (complex(expr.xreplace(cenv))
                                   if expr != 0 else 0j)
                        nops += ops(expr)
                    dsym, dval, dops = {}, {}, 0
                    for k in keep_d:
                        expr, _n, _e = rank(den[k][0], cenv, ctol)
                        dsym[k] = expr
                        dval[k] = (complex(expr.xreplace(cenv))
                                   if expr != 0 else 0j)
                        dops += ops(expr)
                    e2 = sweep_error(nval, dval, reference, sweep)
                    tot = nops + dops + 2 * (len(keep_n) + len(keep_d))
                    print('       subset tol=%-5g coeff tol=%-8.3g -> %5d ops,'
                          ' err %.3e  %s'
                          % (tol, ctol, tot, e2,
                             'PASS' if tot < 200 and e2 <= 0.20 else
                             'holds' if e2 <= 0.20 else 'fail'))
        print()
    return 0


def main():
    reread_old_verdicts()
    transfer_function()
    return 0


if __name__ == '__main__':
    sys.exit(main())
