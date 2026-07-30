"""Stage 9 -- what lowers ``N_eff``?  The readability question, made answerable.

Stage 8 gave the diagnostic that predicts a ranking's *term count*:
``N_eff = A**2 / S2``, the effective number of terms.  The compact uA741 has
``N_eff ~ 194`` and group ranking duly needs ~870 terms.  A *readable* expression
means tens of terms, so it means ``N_eff`` of order ten.

So the question stops being "approximate better" and becomes **"what transformation
lowers N_eff while keeping device symbols?"**  This measures the cheapest untested
lever: **per coefficient of s**.  `doc/cancellation_ranking_conclusions.md` §13
measured per-coefficient ``kappa`` and found it *worse* than the whole determinant's;
``N_eff`` was never looked at, and the two are independent.

The contribution magnitude of each coefficient is reported beside its ``N_eff``
deliberately: a coefficient that is concentrated but carries none of the response is
not a result, and printing them together makes that impossible to claim by accident.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/concentration_profile.py
"""

import math
import sys
import time
import warnings

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix, s_expand

from cancellation_profile import FREQ, NOMINAL, DEGRADED, environment


DEVICE_SYMBOLS = {'gm_q1', 'gm_q2', 'gm_q17', 's'}


def main():
    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q17'))
    print('uA741: dim %d, %g Hz' % (system.dim, FREQ))

    t0 = time.time()
    flat = ddd_of_matrix(system.A)
    expanded = s_expand(system.A, system.s)
    print('compact diagram %d vertices / %d terms; s-expanded degree %d,'
          ' %d vertices  (%.1f s)'
          % (flat.size, flat.term_count(), expanded.degree, expanded.size,
             time.time() - t0))
    print()

    for label, gms in (('nominal', NOMINAL), ('degraded', DEGRADED)):
        env = environment(system, gms)
        s_value = env[system.s]
        coeff_env = {k: v for k, v in env.items() if k != system.s}

        base_k = flat.cancellation(env)
        base_n = flat.concentration(env)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            _e, base_terms, base_err = flat.approximate_groups(env, tol=0.05)
        print('== %s gm ==' % label)
        print('  whole determinant: kappa %.3e, N_eff %.1f, ranking %d terms'
              ' (err %.2e)' % (base_k, base_n, base_terms, base_err))
        print()
        print('  per coefficient of s:')
        print('    %-4s %-12s %-11s %-9s %-13s %s'
              % ('s^k', 'terms', 'kappa', 'N_eff', '|coeff*s^k|', 'share'))

        rows = []
        total_mag = 0.0
        for k in range(expanded.degree + 1):
            coeff = expanded.coefficient(k)
            if coeff.root is None or coeff.term_count() == 0:
                continue
            value = complex(coeff.eval(coeff_env))
            contribution = abs(value * (s_value ** k))
            total_mag += contribution
            rows.append((k, coeff, contribution))

        for k, coeff, contribution in rows:
            share = contribution / total_mag if total_mag else 0.0
            print('    %-4d %-12d %-11.3e %-9.1f %-13.3e %.4f'
                  % (k, coeff.term_count(), coeff.cancellation(coeff_env),
                     coeff.concentration(coeff_env), contribution, share))

        ## Rank the coefficients that actually carry the response.  A small N_eff
        ## in a coefficient contributing 1e-6 of the magnitude is not a result.
        print()
        print('  ranking the coefficients that carry the response (share > 1%):')
        best = None
        for k, coeff, contribution in rows:
            if total_mag and contribution / total_mag <= 0.01:
                continue
            t0 = time.time()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                expr, n, err = coeff.approximate_groups(coeff_env, tol=0.05)
            names = sorted(str(x) for x in expr.free_symbols)
            devices = [x for x in names if x in DEVICE_SYMBOLS]
            print('    s^%-3d N_eff %-8.1f -> %5d terms, err %.2e, %.1f s,'
                  ' symbols: %s'
                  % (k, coeff.concentration(coeff_env), n, err,
                     time.time() - t0,
                     ', '.join(devices) if devices else 'NONE'))
            if best is None or n < best[1]:
                best = (k, n, coeff.concentration(coeff_env), bool(devices))

        if best:
            k, n, neff, has_dev = best
            print()
            print('  best coefficient: s^%d, N_eff %.1f, %d terms (whole'
                  ' determinant needed %d)' % (k, neff, n, base_terms))
            print('  GATE 9-2 (N_eff <= 20 and <= 30 terms with device symbols):'
                  ' %s' % ('PASS' if neff <= 20 and n <= 30 and has_dev
                           else 'FAIL'))
            print('  improvement over the whole determinant: N_eff %.1fx,'
                  ' terms %.1fx' % (base_n / neff if neff else float('inf'),
                                    base_terms / n if n else float('inf')))
        print()

    return 0


if __name__ == '__main__':
    sys.exit(main())
