"""Stage 10 -- how few terms at a tolerance a designer would actually accept?

Every measurement in this plan has used ``tol = 0.05``, inherited from the first
round and never questioned.  Stage 9 got the uA741's dominant coefficient to 155
terms there.  But the hand analyses symbolic tools are compared against are
order-of-magnitude arguments: a designer reading one routinely accepts 20-30%.

So this measures the trade the existing machinery already offers -- terms against
tolerance, on the coefficient that carries the response -- before reaching for
another transformation.  It is nearly free and nobody had looked.

**The expression is printed.**  A small term count at a loose tolerance is only a
result if a person can read the result, and the gate for this stage says so: if the
output is twenty products of six factors each then it is met on its letter and
failed on its purpose, which is exactly what happened to stage A of the predecessor
plan.  Printing it is what makes that callable.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/tolerance_curve.py
"""

import sys
import time
import warnings

import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import s_expand

from cancellation_profile import FREQ, NOMINAL, DEGRADED, environment


TOLERANCES = (0.5, 0.3, 0.2, 0.1, 0.05, 0.01)
DEVICES = {'gm_q1', 'gm_q2', 'gm_q17'}


def dominant_coefficient(system, expanded, env):
    """The coefficient carrying most of the response at this frequency."""
    s_value = env[system.s]
    coeff_env = {k: v for k, v in env.items() if k != system.s}
    best, best_mag = None, 0.0
    for k in range(expanded.degree + 1):
        coeff = expanded.coefficient(k)
        if coeff.root is None or coeff.term_count() == 0:
            continue
        mag = abs(complex(coeff.eval(coeff_env)) * (s_value ** k))
        if mag > best_mag:
            best, best_mag = (k, coeff), mag
    return best, coeff_env


def main():
    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q17'))
    t0 = time.time()
    expanded = s_expand(system.A, system.s)
    print('uA741 at %g Hz; s-expanded in %.1f s' % (FREQ, time.time() - t0))
    print()

    shown = False
    for label, gms in (('nominal', NOMINAL), ('degraded', DEGRADED)):
        env = environment(system, gms)
        (k, coeff), coeff_env = dominant_coefficient(system, expanded, env)
        print('== %s gm: dominant coefficient is s^%d (%d terms, N_eff %.1f) =='
              % (label, k, coeff.term_count(), coeff.concentration(coeff_env)))
        ## Operation count *after* sympy collects the kept groups, which is the
        ## quantity a reader actually faces.  Term count is a poor proxy: five
        ## groups collect into one product of four factors, 11 operations.
        print('  %-8s %-8s %-8s %-12s %s'
              % ('tol', 'terms', 'ops', 'err', 'device symbols'))

        first_small = None
        smallest = None
        for tol in TOLERANCES:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                expr, n, err = coeff.approximate_groups(coeff_env, tol=tol)
            names = {str(x) for x in expr.free_symbols}
            present = sorted(names & DEVICES)
            print('  %-8g %-8d %-8d %-12.3e %s'
                  % (tol, n, int(sympy.count_ops(expr)), err,
                     ', '.join(present) if present else 'NONE'))
            if first_small is None and n <= 30 and err <= 0.20 and present:
                first_small = (tol, n, err, expr)
            if smallest is None or int(sympy.count_ops(expr)) < smallest[1]:
                smallest = (tol, int(sympy.count_ops(expr)), n, err, expr)

        if first_small:
            tol, n, err, expr = first_small
            print('  GATE 10-2 (<=30 terms at <=20%% error, symbols intact): PASS'
                  '  -- tol=%g, %d terms, err %.3e' % (tol, n, err))
            if not shown:
                shown = True
                print()
                print('  THE EXPRESSION, so "readable" can be judged not asserted:')
                terms = sympy.Add.make_args(expr)
                for i, term in enumerate(terms):
                    print('    %s %s' % ('  ' if i else '+ ', term))
                print()
                print('  %d terms, %d sympy operations in total'
                      % (len(terms), int(sympy.count_ops(expr))))
        else:
            print('  GATE 10-2: FAIL -- nothing reached 30 terms at 20%% error')
        ## Print the most compact result regardless of whether the gate passed:
        ## the gate is about a threshold, this is about whether anyone can read it.
        if smallest and not shown:
            shown = True
            tol, ops, n, err, expr = smallest
            print()
            print('  MOST COMPACT RESULT (tol=%g): %d groups, %d operations,'
                  ' err %.1f%%' % (tol, n, ops, 100 * err))
            print('  collected by sympy, so "readable" can be judged:')
            print('      %s' % expr)
        print()

    return 0


if __name__ == '__main__':
    sys.exit(main())
