"""What matrix-free shooting could buy, before anyone builds it.

RECORDED SCOPE ITEM 6 (Telichevesky, Kundert & White, DAC 1995).  The
standard argument is about the FACTORISATION: Kundert (ICCAD'97) puts
forming `J_phi` at O(N^2 S) and factoring it at O(N^3), "intractable when N
exceeds several hundred", and names matrix-implicit Krylov as the answer.

⚠ BUT IN THIS IMPLEMENTATION THE FACTORISATION IS NOT OBVIOUSLY THE COST.
`J_phi` is factored ONCE per shooting iteration.  The sensitivity
propagation that BUILDS it runs `linearsolver(Jf, S)` once per timestep
with an m-column right-hand side -- m back-substitutions per step, `N m`
over a period, against the factorisation's single O(m^3).  Matrix-free
replaces the propagation with ONE column per Krylov iteration, `k N`, so
the ceiling on the whole idea is `m / k`.

This measures where the time actually goes, by the only method that
settles it: run a traversal with the propagation and again without it, on
the same trajectory, at several circuit sizes.  Assume the propagation is
infinitely fast and see what the run does -- if the ceiling is not worth
the work, that is the answer in minutes and no design was needed.

Interleaved and min-of-N, because this box carries other work.
"""
import time
import warnings

import numpy as np

from pycircuit import circuit
from pycircuit.circuit import SubCircuit, gnd, R, C, L, VSin
from pycircuit.circuit.shooting import PSS

NPTS = 50


def rc_ladder(sections):
    """A driven RC ladder: `m` grows linearly with `sections`.

    Deliberately linear and benign -- this measures where TIME goes, not
    whether a circuit is hard, and a nonlinear circuit would add
    inner-Newton variance that has nothing to do with the question.
    """
    c = SubCircuit()
    c.add_node('n0')
    c['vs'] = VSin('n0', gnd, va=1.0, freq=1e3)
    for k in range(sections):
        a, b = 'n%d' % k, 'n%d' % (k + 1)
        c.add_node(b)
        c['R%d' % k] = R(a, b, r=1e3)
        c['C%d' % k] = C(b, gnd, c=1e-9)
    return c


def _prepared(cir, npts, method='gear'):
    """A PSS with its `Transient` built, ready to traverse."""
    circuit.default_toolkit = circuit.numeric
    pss = PSS(cir, method=method, reltol=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / npts, maxiterations=2)
    return pss


def solve_shares(pss, npts):
    """Fractions of a traversal spent in the two kinds of linear solve.

    ⚠ MEASURED FROM INSIDE ONE RUN, ON PURPOSE.  The obvious design --
    time a traversal with the propagation and again without it -- compares
    two SEPARATE runs, and on a loaded box that is worthless: it produced a
    NEGATIVE propagation cost (the "off" variant timing slower than "on"),
    and across three attempts at m>=110 between three and four pairs in
    thirteen came out below 1.0, which is physically impossible.  Load
    average was 17-27 throughout.

    Timing the solves from inside the run puts load in the numerator AND
    the denominator, so it cancels.  The structure being measured is exact
    and load-independent anyway: per step the traversal does TWO one-column
    solves (the inner Newton) and ONE solve with 2m columns (the
    sensitivity propagation), counted directly.
    """
    times, hs = pss._period_grid(1e-3, npts, None)
    m = pss.cir.n - 1
    tk = pss.toolkit
    acc = {'one': 0.0, 'wide': 0.0}
    orig = tk.linearsolver

    def spy(A, b, *a, **k):
        w = 1 if np.asarray(b).ndim == 1 else np.asarray(b).shape[1]
        t0 = time.perf_counter()
        r = orig(A, b, *a, **k)
        acc['one' if w == 1 else 'wide'] += time.perf_counter() - t0
        return r

    tk.linearsolver = spy
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            t0 = time.perf_counter()
            pss._traverse_solved_history(np.zeros(m), np.zeros(m), times, hs)
            total = time.perf_counter() - t0
    finally:
        tk.linearsolver = orig
    return total, acc['one'] / total, acc['wide'] / total


def main():
    print('one traversal, %d points/period, Gear-2 (2m-column sensitivity)'
          % NPTS)
    print('"propagation" is the wide solve matrix-free would replace\n')
    hdr = ('%5s | %10s | %8s %10s | %12s'
           % ('m', 'traversal', 'newton', 'propagation', 'ceiling k=20'))
    print(hdr)
    print('-' * len(hdr))
    for sections in (38, 108, 240):
        pss = _prepared(rc_ladder(sections), NPTS)
        m = pss.cir.n - 1
        total, fn, fp = solve_shares(pss, NPTS)
        best = (1.0 - fp) + fp * (20.0 / m if m > 20 else 1.0)
        print('%5d | %10.4f | %7.1f%% %9.1f%% | %11.2fx'
              % (m, total, 100 * fn, 100 * fp, 1.0 / best))

    print('\nRESULT (2026-09-01), and it is PARTIAL.')
    print('  m=40 is stable and small: the propagation is 2.2% of a')
    print('  traversal by this method and 5.8-6.1% by an independent one')
    print('  (timing a traversal with the propagation switched off), so the')
    print('  ceiling there is 1.01x to 1.03x either way.  Matrix-free buys')
    print('  nothing at that size.')
    print('')
    print('  ⚠ m >= 110 IS NOT MEASURABLE ON THIS MACHINE, and three')
    print('  methods failed to make it so.  Load average sat at 17-27')
    print('  throughout.  The with/without comparison produced a NEGATIVE')
    print('  propagation cost and three-to-four of thirteen paired ratios')
    print('  below 1.0 -- physically impossible.  This in-run method, which')
    print('  should cancel load, still swung 2.9% -> 19.1% at m=110 between')
    print('  runs and reported m=110 and m=242 taking the same time.  None')
    print('  of those numbers should be quoted.')
    print('')
    print('  SO THE ITEM IS NOT REFUSED AND NOT JUSTIFIED -- it is')
    print('  unmeasured above m=40.  Re-run this on a quiet box; the')
    print('  gate is whether the propagation share passes ~30%, and the')
    print('  structural counts below are load-independent and do hold:')
    print('  per step the traversal does TWO one-column solves (the inner')
    print('  Newton) and ONE solve with 2m columns (the propagation).')
    print('')
    print('  ⚠ And the shape to expect: at m=40 BOTH kinds of solve')
    print('  together are under 5% of a traversal, so most of the time is')
    print('  ASSEMBLY -- device evaluation and matrix building -- which')
    print('  matrix-free does not touch.  Whatever the crossover is, the')
    print('  O(N^3) argument alone will overstate it for this code.')


if __name__ == '__main__':
    main()
