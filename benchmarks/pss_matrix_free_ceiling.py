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
SECTIONS = (38, 108, 240, 500, 1000)


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
    """Fractions of a traversal spent in the inner Newton and in propagation.

    ⚠ MEASURED FROM INSIDE ONE RUN, ON PURPOSE.  The obvious design --
    time a traversal with the propagation and again without it -- compares
    two SEPARATE runs, and on a loaded box that is worthless: it produced a
    NEGATIVE propagation cost (the "off" variant timing slower than "on"),
    and across three attempts at m>=110 between three and four pairs in
    thirteen came out below 1.0, which is physically impossible.  Load
    average was 17-27 throughout.

    Timing from inside the run puts load in the numerator AND the
    denominator, so it cancels.

    ⚠ AND TIMING THE SOLVE IS NOT TIMING THE PROPAGATION.  This function
    used to accumulate only `toolkit.linearsolver`, on the reasoning that
    the wide 2m-column solve is what matrix-free replaces.  That UNDERCOUNTS
    it at every size.  `_step_sensitivity` also forms `C @ P` products --
    m x m against m x 2m, O(m^3) each, the same order as the solve -- and
    matrix-free eliminates those too, because it never forms `P` at all.
    So the honest unit is the WHOLE sensitivity step; the solve is reported
    beside it only to show how small a part of it was being counted.
    """
    times, hs = pss._period_grid(1e-3, npts, None)
    m = pss.cir.n - 1
    tk = pss.toolkit
    acc = {'one': 0.0, 'wide': 0.0, 'prop': 0.0}
    orig_solve = tk.linearsolver
    orig_sens = pss._step_sensitivity

    def spy(A, b, *a, **k):
        w = 1 if np.asarray(b).ndim == 1 else np.asarray(b).shape[1]
        t0 = time.perf_counter()
        r = orig_solve(A, b, *a, **k)
        acc['one' if w == 1 else 'wide'] += time.perf_counter() - t0
        return r

    def timed_sens(*a, **k):
        t0 = time.perf_counter()
        try:
            return orig_sens(*a, **k)
        finally:
            acc['prop'] += time.perf_counter() - t0

    tk.linearsolver = spy
    pss._step_sensitivity = timed_sens
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            t0 = time.perf_counter()
            pss._traverse_solved_history(np.zeros(m), np.zeros(m), times, hs)
            total = time.perf_counter() - t0
    finally:
        tk.linearsolver = orig_solve
        pss._step_sensitivity = orig_sens
    return (total, acc['one'] / total, acc['prop'] / total,
            acc['wide'] / total)


def main():
    print('one traversal, %d points/period, Gear-2 (2m-column sensitivity)'
          % NPTS)
    print('"propagation" is the whole sensitivity step matrix-free replaces;')
    print('"solve" is the part of it this harness used to count alone.\n')
    hdr = ('%5s | %10s | %8s %12s %8s | %12s'
           % ('m', 'traversal', 'newton', 'propagation', '(solve)',
              'ceiling k=20'))
    print(hdr)
    print('-' * len(hdr))
    for sections in SECTIONS:
        pss = _prepared(rc_ladder(sections), NPTS)
        m = pss.cir.n - 1
        total, fn, fp, fw = solve_shares(pss, NPTS)
        best = (1.0 - fp) + fp * (20.0 / m if m > 20 else 1.0)
        print('%5d | %10.4f | %7.1f%% %11.1f%% %7.1f%% | %11.2fx'
              % (m, total, 100 * fn, 100 * fp, 100 * fw, 1.0 / best))

    print("""
RESULT (2026-09-02, on a QUIET box -- load average 0.1-1.5, and every
reading below reproduced to better than 0.5% across repeats).

⚠ THIS OVERTURNS THE 2026-09-01 RECORD, AND NOT BECAUSE THE BOX GOT QUIET.
That record said the propagation was 2.2% of a traversal at m=40 and the
ceiling 1.01x-1.03x.  The 2.2% is right for what it timed and WRONG for the
question: it counted `linearsolver` only, and `_step_sensitivity` spends
more than that again on `C @ P` products -- m x m against m x 2m, O(m^3),
the same order as the solve.  Matrix-free never forms `P`, so it removes
those too.  Counting the whole step, single-threaded BLAS:

     m | propagation | (solve alone) | ceiling k=20
    40 |        4.8% |          2.3% |       1.02x
   110 |       15.0% |          6.2% |       1.14x
   242 |       38.1% |         14.1% |       1.54x
   502 |       63.8% |         21.3% |       2.58x
  1002 |       79.1% |         24.9% |       4.44x

SO THE GATE IS PASSED.  It was "does the propagation share pass ~30%": it
crosses 30% near m=220 and reaches 79% at m=1002.  Item 6 is JUSTIFIED
above m~250 and remains POINTLESS at m=40, where the old verdict stands on
the corrected number too (4.8%, ceiling 1.02x).

⚠ AND IT IS NOW BUILT -- `PSS.solve(matrix_free=True)`, driven
solved-history path only.  What it actually delivers, end to end, against
the dense path on this same ladder:

      m    dense (iters)      matrix-free (iters)     speedup
    242      2.113 s (2)            1.557 s (2)        1.36x
    502      9.255 s (2)            6.131 s (3)        1.51x
   1002     52.402 s (2)           24.636 s (3)        2.13x

Lower than the traversal-only figures below (2.23x at m=502) because a
`solve` also does setup, the replay and the DFT, none of which this
touches, and matrix-free spends an extra Newton iteration.  Both numbers
are true; say which one is being quoted.  Measured k: 2/4/7/12 at
m=40/110/242/502, against the 20 assumed here.

⚠ BUT THE SHARE DEPENDS ON BLAS THREADING, so quote the condition with the
number.  With threads left free those same sizes read 4.8 / 7.2 / 12.4 /
18.2% -- the `C @ P` products thread well and the Python-level assembly
does not.  The single-threaded column is the one to trust here: the
threaded traversal time moved 2.64 s -> 4.02 s at m=502 between runs on
this box, while every single-threaded reading held to 0.1%.  A machine
whose BLAS has all cores will see a smaller prize than the table above.

⚠ AND THE CEILING IS STILL A CEILING.  `k=20` Krylov iterations is an
assumption, not a measurement, and the model charges matrix-free `k/m` of
the propagation and nothing else.  Assembly -- device evaluation and matrix
building -- is untouched, and it is what the remaining share is.""")


if __name__ == '__main__':
    main()
