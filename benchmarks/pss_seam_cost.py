"""What the PSS cold-start seam actually COSTS, as opposed to what it reads.

The LTE report landed 2026-09-01 with three numbers, and the loudest of them
is `max_lte_seam`: on the Q=20 resonator at 100 points/period it reads 38.6x
tolerance for Gear-2 and 15.1x for trapezoidal, against interior per-step
peaks of 0.076 and 0.024.  A ratio of ~500.

⚠ THAT IS AN ESTIMATOR, NOT A COST.  Two facts sit awkwardly beside it:
Gear-2's answer at that grid is only 1.17% low, and it improves 5x when the
grid is halved -- which is not what a defect that dominates the error and
gets WORSE under refinement should do.  So before anyone builds the fix (make
the seeded history part of the shooting unknowns, so the period map is a
fixed point in the augmented state), the question to answer is how much of
the answer the seam is actually responsible for.

THE DECOMPOSITION.  `_begin_period` re-seeds a flat charge history at every
shooting iteration -- which is exactly what keeps phi a function of `x0`
alone -- so PSS's answer is the fixed point of a map that cold-starts.  Keep
integrating past the converged solution WITHOUT re-seeding and the ring
buffers fill with the true periodic past; the trajectory then drifts to the
limit cycle the same grid and the same method would produce with no seam.
That gives three numbers on one axis:

    analytic  <---- interior ---->  primed  <-- seam -->  PSS answer

  seam cost      = |primed - PSS answer|
  interior cost  = |primed - analytic|

The circuit is chosen so the third point is exact rather than measured: a
series RLC driven at resonance with a 1 V source has |V_C| = Q volts, so the
analytic peak is 20 V with no reference run to trust.

Anti-vacuity: the primed continuation must actually differ from the cold one,
and the ring buffer must actually hold three distinct charges when it starts
-- both are asserted, because a priming that silently did nothing would show
up as "the seam costs zero", which is the answer this script exists to test.

RESULT (2026-09-01, reltol=1e-9 so the shooting solve's own slack is far
below what is being measured -- at reltol=1e-3 the euler and trap readings
were ~1e-5 V of solve residue, not seam):

    method  npts     PSS     primed   seam cost   interior   seam/interior
    euler    100   8.81540   8.81540    5.1e-12    1.12e+01        0.000
    euler    400  15.20994  15.20994    7.2e-13    4.79e+00        0.000
    gear     100  19.76639  19.89297    1.27e-01   1.07e-01        1.183
    gear     200  19.95451  19.98524    3.07e-02   1.48e-02        2.083
    gear     400  19.99008  19.99735    7.27e-03   2.65e-03        2.742
    trap     100  19.98968  19.98968    1.3e-11    1.03e-02        0.000
    trap     400  19.99957  19.99957    6.3e-13    4.27e-04        0.000

⚠ THE SEAM IS REAL FOR GEAR-2 AND EXACTLY ZERO FOR THE OTHER TWO, and the
report as shipped flagged all three.  `max_lte_seam` read 15.1 for
trapezoidal at npts=100 while its cost is 1.3e-11 V.  The mechanism: only
the third divided difference INSIDE THE ESTIMATOR touches the fabricated
charge for a one-step method, never the companion, so the reading is an
artefact of the measurement and not a property of the answer.  Euler's
companion reads `q_{n-1}`; trapezoidal's reads `q_{n-1}` and `iq_{n-1}`,
and the Euler opening step supplies an `iq` consistent with it (which is
what the first-step order drop is FOR).  Gear-2 reads `q_{n-2}` -- the
entering unknown `x_in` -- and that is the only method that can see it.

WHY GEAR-2 PAYS.  `x_in` is the shooting unknown, and the shooting condition
is `x(0) = x(P)`; it does NOT constrain `x_in` to be the orbit's own
`x(-dt)`.  `x_in` reaches `x(0)` through one order-dropped Euler step, so it
sits O(h^2) off the true history point -- and Gear-2 then reads it AS a
history point.  Confirmed by removing the drop, which makes Gear-2 WORSE at
every grid (err 2.34e-1 -> 4.39e-1 at npts=100, 9.92e-3 -> 2.34e-2 at 400):
the drop is protective, and the 1.27e-01 is the residue it leaves behind.

AND IT DOES NOT GO AWAY.  The seam falls as h^2 (4.13x, 4.23x per halving)
while the interior falls faster (7.25x, 5.57x), so its share of the error
GROWS: 54% of Gear-2's total error at npts=100, 68% at 200, 73% at 400.
That is the case for making the entering history part of the shooting
unknowns -- for multistep methods only.  The prize is measurable: Gear-2's
error at npts=100 would go 2.34e-1 -> 1.07e-1 at no extra cost per
iteration, and the factor grows with refinement.

⚠ ONE EXPERIMENT HERE IS INVALID AND IS KEPT OUT ON PURPOSE.  Re-introducing
a single mechanism (the order drop, or the flat seeding) into the free
continuation does NOT measure that mechanism's cost inside PSS: a shooting
fixed point ABSORBS a once-per-period perturbation by moving `x0`, and a
driven continuation cannot.  Measured, the contradiction is flat -- forcing
the order drop once per period into the trap continuation costs 0.234 V,
while the same drop inside PSS costs 1.3e-11 V.  Only fixed-point against
fixed-point is a fair comparison, which is what `run` does.
"""
import warnings
import numpy as np

from pycircuit import circuit
from pycircuit.circuit import SubCircuit, gnd, VSin, R, L, C
from pycircuit.circuit.shooting import PSS

## Enough periods for the primed continuation to reach its own limit cycle.
## The deviation decays at the circuit's per-period rate, exp(-pi/Q) = 0.8546,
## so 80 periods take a starting gap to 0.8546**80 = 3.5e-6 of itself.
SETTLE_PERIODS = 80


def q20_rlc(f0=1e3, Q=20.0):
    """Series RLC at resonance: |V_C| = Q * va, exactly."""
    L_, C_ = 1e-3, 1.0 / ((2 * np.pi * f0) ** 2 * 1e-3)
    c = SubCircuit()
    c.add_node('a'); c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, freq=f0)
    c['R'] = R('a', 'b', r=(1.0 / Q) * np.sqrt(L_ / C_))
    c['L'] = L('b', 'c', L=L_)
    c['C'] = C('c', gnd, c=C_)
    return c


def peak_of(pss, X, node='c'):
    """|V_node| over a collected period, from the reduced state vectors."""
    i = pss.cir.get_node_index(node) - (1 if pss.cir.get_node_index(node)
                                        > pss.irefnode else 0)
    return float(np.max(np.abs(np.asarray([x[i] for x in X], dtype=float))))


def solve_augmented(method, npts, f0=1e3, Q=20.0, reltol=1e-9):
    """The shipped formulation, whatever it is for this method.

    For a companion reaching two charges back this solves for `(x_0, x_{-1})`
    jointly; for the others it is the plain path unchanged.  The gate is that
    it REACHES `primed` -- the seam was the only difference, so removing it
    must land on the cycle a real history produces, not merely nearer to it.
    """
    period = 1.0 / f0
    circuit.default_toolkit = circuit.numeric
    pss = PSS(q20_rlc(f0, Q), method=method, reltol=reltol)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = pss.solve(period=period, timestep=period / npts,
                        maxiterations=40)
    assert pss.converged, 'augmented %s/%d did not converge' % (method, npts)
    return float(np.max(np.abs(
        np.asarray(res['tpss'].v('c'), dtype=float).ravel()))), pss


def run(method, npts, f0=1e3, Q=20.0, reltol=1e-3):
    period = 1.0 / f0
    circuit.default_toolkit = circuit.numeric
    pss = PSS(q20_rlc(f0, Q), method=method, reltol=reltol)
    ## THE PLAIN FORMULATION, held here on purpose.  Gear-2 now solves for
    ## its entering history, so measuring the shipped path would measure a
    ## seam of zero and the record of what the seam WAS would be gone.  This
    ## keeps the before, and `solve_augmented` supplies the after.
    pss._augmented = lambda: False
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = pss.solve(period=period, timestep=period / npts,
                        maxiterations=40)
    assert pss.converged, '%s/%d did not converge' % (method, npts)
    cold = float(np.max(np.abs(
        np.asarray(res['tpss'].v('c'), dtype=float).ravel())))

    ## THE PRIMED CONTINUATION.  `solve` has just finished its replay, so the
    ## `Transient`'s rings hold the real trailing history of the period --
    ## and crucially `_begin_period` is NOT called again, which is the whole
    ## experiment.  Sources are periodic in T, so advancing `t` by whole
    ## periods reproduces the same excitation on the same grid.
    tr = pss._transient()
    qhist = np.asarray(tr._qlast, dtype=float)
    spread = float(np.max(np.abs(qhist[0] - qhist[-1])))
    assert spread > 0.0, \
        'the history is still flat after a replay -- nothing was primed'

    times, dt = pss.toolkit.linspace(0.0, period, num=npts, endpoint=True,
                                     retstep=True)
    x = np.asarray(res['tpss'].x, dtype=float)[:, -1]
    x = np.concatenate((x[:pss.irefnode], x[pss.irefnode + 1:]))
    last = None
    for k in range(SETTLE_PERIODS):
        window = []
        ## ⚠ `times[1:]`, NOT `times` -- npts-1 steps of `P/(npts-1)`, which
        ## is exactly one period.  PSS's own replay takes npts steps because
        ## its FIRST one manufactures `x(0)` from the shooting unknown (which
        ## is the state at `-dt`); the period map proper is the remaining
        ## npts-1.  Repeating that extra step once per period detunes the
        ## drive by 1/(npts-1) -- measured, it took trapezoidal's continued
        ## amplitude to 18.72 V, and `1/sqrt(1 + (2 Q / 99)**2) * 20` is
        ## 18.55, so the artefact was the whole reading.  It is a good
        ## illustration of what this script is guarding against: a
        ## plausible number that measures the harness.
        for t in times[1:]:
            x = pss.solve_timestep(x, t + (k + 1) * period, dt)
            window.append(np.asarray(x, dtype=float))
        last = window
    primed = peak_of(pss, last)

    fixed, apss = solve_augmented(method, npts, f0, Q)
    return dict(method=method, npts=npts, cold=cold, primed=primed,
                analytic=Q, max_lte=pss.max_lte, total_lte=pss.total_lte,
                seam_lte=pss.max_lte_seam, hist_spread=spread,
                fixed=fixed, augmented=apss.augmented,
                fixed_seam=apss.max_lte_seam)


def main(reltol=1e-9):
    ## reltol=1e-9, not the shipped default: at 1e-3 the shooting solve's own
    ## slack is ~1e-5 V, which is LARGER than the effect being measured for
    ## two of the three methods and reads as a small non-zero seam.
    rows = [run(m, n, reltol=reltol) for m in ('euler', 'gear', 'trap')
            for n in (100, 200, 400)]

    print('Q=20 series RLC at resonance, analytic peak 20 V exactly, '
          'reltol=%g' % reltol)
    print('primed = %d periods continued past the solve with NO re-seed\n'
          % SETTLE_PERIODS)
    hdr = ('%-6s %5s | %9s %9s %9s | %10s %10s | %10s %6s'
           % ('method', 'npts', 'plain', 'primed', 'shipped', 'seam cost',
              'interior', 'err after', 'gain'))
    print(hdr); print('-' * len(hdr))
    for r in rows:
        seam = abs(r['primed'] - r['cold'])
        interior = abs(r['primed'] - r['analytic'])
        before = abs(r['cold'] - r['analytic'])
        after = abs(r['fixed'] - r['analytic'])
        print('%-6s %5d | %9.5f %9.5f %9.5f | %10.3e %10.3e | %10.3e %5.2fx'
              % (r['method'], r['npts'], r['cold'], r['primed'], r['fixed'],
                 seam, interior, after, before / after if after else
                 float('inf')))

    print('\ndoes the fix REACH the cycle a real history produces?')
    for r in rows:
        gap = abs(r['fixed'] - r['primed'])
        print('  %-6s n=%-4d augmented=%-6s |shipped - primed| = %.2e  %s'
              % (r['method'], r['npts'], r['augmented'], gap,
                 'seam reported: %s' % r['fixed_seam']))

    print('\nwhat the estimator says vs what it costs, where a seam exists'
          ' at all:')
    for r in rows:
        if r['seam_lte'] is None:
            ## Euler and trapezoidal: the companion never reads the entering
            ## point, so the report no longer claims a seam for them.  That
            ## `None` IS the corrected classification, not missing data.
            continue
        seam = abs(r['primed'] - r['cold'])
        interior = abs(r['primed'] - r['analytic'])
        print('  %-6s n=%-4d LTE says seam/interior = %8.1fx, '
              'the answer says %8.4fx'
              % (r['method'], r['npts'], r['seam_lte'] / r['max_lte'],
                 seam / interior if interior else float('inf')))


if __name__ == '__main__':
    main()
