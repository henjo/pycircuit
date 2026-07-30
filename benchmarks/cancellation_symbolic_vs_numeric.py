"""How much of the uA741's cancellation is *symbolic*, and so removable?

Tan & Shi (BMAS 2002), "Parametric Analog Behavioral Modeling Based on
Cancellation-Free DDDs", draw a distinction this project had not:

* **symbolic cancellation** -- terms that cancel as symbols, arising from the MNA
  formulation and from device matching.  Detectable from local 2x2 matrix
  patterns (their Fig. 5, e.g. ``[[p, -p], [-p, p]]`` from a floating resistor)
  and removable *during* the construction of the s-expanded DDD.  They report
  "up to 70-90 percent product terms can be canceling terms for a typical
  analog circuit".
* **numerical cancellation** -- terms that cancel only once numbers are put in.
  They defer this to the dominant-term generation step.

`DDD.cancellation` measures the **numerical** kind: kappa = sum|term|/|sum|.  So
the open question, and the one that decides whether the de-cancellation route is
worth building, is how much of the uA741's kappa = 9.4e3 survives once the
symbolic part is gone.

This script bounds that without implementing de-cancellation, by comparing three
representations of the *same* circuit:

1. **compact-symbol complex DDD** -- what pycircuit builds today: one vertex per
   whole matrix entry, ``g1 + g2 + s*c1``.  The kappa already recorded.
2. **s-expanded coefficient DDDs** -- ``s_expand`` splits each entry by power of
   s.  This removes the cancellation that comes from mixing powers of s at one
   frequency, and nothing else.  Reported per coefficient.
3. **the s-expanded determinant reassembled at the same frequency** -- sum_k
   coeff_k * s**k, as a check that 2 describes the same number as 1.

What this can and cannot settle, stated plainly: if the coefficient DDDs already
have kappa near 1, then s-expansion alone buys the property and the rest of the
de-cancellation machinery is unnecessary.  If they do not, the residue is
cancellation *within* one power of s -- which is what the MNA patterns cause, and
what a de-cancellation step would target.  Either way this is a bound, not the
answer a real implementation would give.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/cancellation_symbolic_vs_numeric.py
"""

import math
import sys
import time

import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix, s_expand

from cancellation_profile import FREQ, NOMINAL, DEGRADED, environment


def main():
    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q17'))
    print('uA741: %dx%d, %d symbols, %g Hz'
          % (system.dim, system.dim, len(system.A.free_symbols), FREQ))

    t0 = time.time()
    flat = ddd_of_matrix(system.A)
    print('compact-symbol complex DDD: %d vertices, %d terms, %.2f s'
          % (flat.size, flat.term_count(), time.time() - t0))

    t0 = time.time()
    expanded = s_expand(system.A, system.s)
    print('s-expanded DDD: degree %d, %d vertices total, %.2f s'
          % (expanded.degree, expanded.size, time.time() - t0))
    print()

    for label, gms in (('nominal', NOMINAL), ('degraded', DEGRADED)):
        env = environment(system, gms)
        s_value = env[system.s]
        ## The coefficient diagrams carry no s, so hand them an environment
        ## without it -- passing one would not be wrong, but it would invite the
        ## reader to think s mattered here.
        coeff_env = {k: v for k, v in env.items() if k != system.s}

        k_flat = flat.cancellation(env)
        print('== %s gm ==' % label)
        print('  1. compact-symbol complex DDD: kappa = %.3e' % k_flat)

        print('  2. s-expanded coefficients, one kappa each:')
        print('     %-4s %-9s %-12s %-11s %s'
              % ('k', 'vertices', 'terms', 'kappa', '|coeff * s**k|'))
        total = 0j
        mass = 0.0
        worst = 0.0
        for k in range(expanded.degree + 1):
            coeff = expanded.coefficient(k)
            if coeff.root is None:
                continue
            terms = coeff.term_count()
            if terms == 0:
                continue
            kap = coeff.cancellation(coeff_env)
            value = complex(coeff.eval(coeff_env))
            contribution = value * (s_value ** k)
            total += contribution
            mass += abs(contribution)
            if math.isfinite(kap):
                worst = max(worst, kap)
            print('     %-4d %-9d %-12d %-11.3e %.3e'
                  % (k, coeff.size, terms, kap, abs(contribution)))

        ## Reassembly check: does the s-expanded form describe the same number?
        exact = complex(flat.eval(env))
        rel = abs(total - exact) / abs(exact)
        print('  3. reassembled sum_k coeff_k*s**k vs complex DDD: rel %.2e'
              % rel)
        ## Cancellation *between* powers of s, which s-expansion does not remove
        ## but which per-coefficient approximation never has to face.
        print('     kappa between powers of s alone = %.3e'
              % (mass / abs(total)))
        print('  worst per-coefficient kappa = %.3e   (vs %.3e for the'
              ' compact form)' % (worst, k_flat))
        print('  => s-expansion alone changes kappa by %.2fx at best'
              % (k_flat / worst if worst else float('inf')))
        print()

    print('Reading: a per-coefficient kappa still far above 1 means the residue')
    print('is cancellation within one power of s -- the MNA/device-matching kind')
    print('that Tan & Shi remove during construction.  That is the quantity a')
    print('de-cancellation implementation would have to move.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
