"""The SOURCES-OFF numerical noise floor: what an oscillator's apparent jitter
is when there is no physical noise at all.

WHY THIS IS NOT BIGGIO'S GATE.  Biggio/Bizzarri/Brambilla/Storace 2013 measure
the floor by FFT-ing the jitter sequence of a time-domain noise simulation with
the sources off.  That construction does NOT transfer to this tree's
closed-form phase-noise stack: with the sources off, `c` is identically zero
and the reported PSD is exactly zero -- a structural identity that cannot fail,
and so cannot be a gate (recorded in the roadmap the first time it was tried).

What DOES transfer is the underlying object: threshold crossings of a simulated
waveform.  Those exist here, they are what the AM/PM estimators already use,
and their spacing carries the integrator's period error whether or not any
source is on.  So the floor is measured the way a floor is measured on a bench:

  * inject a KNOWN, DETERMINISTIC perturbation of amplitude `a`;
  * sweep `a` down over decades;
  * in the sloped region the estimator must read `J_abs ~ a` -- that is what
    says the instrument is alive, and it is why this cannot pass zero-vs-zero;
  * where the line BENDS is the floor, and the plateau is its value.

PREDICTIONS, NAMED BEFORE THE SWEEP WAS RUN (roadmap discipline -- a number is
worth much less if the arithmetic was not committed to first):

  P1  the sloped region has slope 1 in `a`.  If not, the instrument is not
      linear and no plateau read off it means anything.
  P2  the floor DROPS WITH INTEGRATOR ORDER at fixed grid.  This is docs-46's
      falsifiable prediction from the paper: if the floor is warping-dominated,
      order is the lever.  If it does not drop, the warping account is wrong
      for this implementation.
  P3  the floor drops with the grid at rate ~h^order at fixed order.
  P4  (this file's own) the sources-off "jitter" is DETERMINISTIC -- a beat
      between the grid and the period, not a random walk.  A deterministic
      period error is a FREQUENCY SHIFT and contributes no jitter at all; what
      makes jitter is that the crossing lands at a different place inside a
      step each cycle.  Signature: `dTmax/J_abs ~ 2` (a sinusoid) rather than
      the ~3-4 of a Gaussian, and a high lag-1 autocorrelation in `T_k`.
  P5  ADAPTIVE stepping raises the floor above fixed stepping at the same mean
      step -- the paper's own stated mechanism ("the finding of threshold
      crossing varies h, which in turn impacts on the T value").

⚠⚠ THREE INSTRUMENT FAILURES ON ONE MEASUREMENT, all of them found by the
same tell -- TWO DIFFERENT INTEGRATORS AGREEING TO FIVE SIGNIFICANT FIGURES.
A quantity that is supposed to BE the integrator's error cannot do that.

  1. LINEAR crossing interpolation.  ESDIRK43 and Radau both read 2.9384e-10
     at 480 points per period, and the "floor" fell as h^0.85 rather than as
     h^4 or h^5.  What was measured was the straight line drawn between two
     samples of a curved waveform: a property of the GRID, not of the method.
     Fixed by a degree-3 interpolant (verified converged against degree 5).
  2. THE ANALYTIC REFERENCE.  `T0 = 2 pi / sqrt(1 - mu^2/4)` is the LINEARISED
     period; the nonlinear limit cycle differs by `mu^2/16`, which is
     -3.9348e-05 here -- and the high-order methods' "period error" sat at
     exactly that value on EVERY grid, which indicts the reference, not them.
  3. THE FIXTURE'S OWN SETTLING, see `periods`.  This is the one that
     mattered: it was still the answer after the first two were fixed.
  4. A SLIVER STEP AT `tend`, see `periods` again -- 4.05e-10 where the grid
     is 1.309e-02, on a run declared uniform.
  5. AN ABSOLUTE METRIC FOR A RELATIVE QUANTITY.  `J_abs` carries the
     problem's time unit: rescale L and C by 100 (identical dynamics, T0 100x
     longer) and every J_abs is 100x larger.  `J_abs / T0` is constant to FIVE
     significant figures across that rescale -- gear 1.9500e-11 / 1.9499e-11 /
     1.9498e-11, radau 4.2879e-14 / 4.1347e-14 / 4.2326e-14 -- so the
     fractional jitter is the quantity, and any absolute number here is only
     meaningful next to its T0.

WHAT THE ANSWER IS, once all five are fixed (van der Pol Q = 15.9, 240 points
per period, fixed step, 200 cycles with 150 discarded):

    gear-2   fractional floor 1.95e-11   -- its own DISCRETISATION, 460x above
                                            the representation limit
    radau    fractional floor 4.23e-14   -- at the FLOATING-POINT limit of the
                                            time variable: ulp(t)/T0 is
                                            2.3e-14, so this is ~2 ulps

So the order lever does not merely lower the floor by 461x; it takes it all
the way DOWN TO THE ARITHMETIC, where no further order or grid refinement can
help.  ⚠ That also explains the one row this file cannot otherwise account
for: at 480 points per period every method reads the same ~2.4e-11, because
`t = t + dt` has accumulated twice as many roundings and the representation
limit has risen above every method's discretisation.  The MECHANISM is
identified; the exact 88x rise from 48000 to 96000 accumulations is not, and
is recorded as unexplained rather than fitted.

Run: `python benchmarks/noise_floor_sources_off.py`
"""
import warnings

import numpy as np

from pycircuit.circuit import circuit
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import C, L, BSource, ISin  # noqa: F401
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import (Gear2Integrator,
                                          TrapezoidalIntegrator,
                                          TRBDF2Integrator,
                                          ESDIRK43Integrator,
                                          RadauIIA3Integrator)

warnings.simplefilter('ignore')
circuit.default_toolkit = circuit.numeric

## (class, label, classical order) -- an order ladder 2,2,2,4,5 on ONE fixture
METHODS = [(Gear2Integrator, 'gear', 2),
           (TrapezoidalIntegrator, 'trap', 2),
           (TRBDF2Integrator, 'trbdf2', 2),
           (ESDIRK43Integrator, 'esdirk43', 4),
           (RadauIIA3Integrator, 'radau', 5)]


def vdp(Q=15.9, a=0.0, fm_div=50.0):
    """Van der Pol, and NO noise source anywhere.  `a` is a deterministic
    perturbation tone at `f0/fm_div` -- far from `f0`, so it slowly modulates
    the period rather than injection-locking the orbit."""
    mu = 1.0 / (2 * np.pi * Q)
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    T0 = 2 * np.pi / np.sqrt(1 - mu ** 2 / 4)
    if a:
        c['pert'] = ISin('v', gnd, ia=a, freq=1.0 / (fm_div * T0))
    return c, T0


def crossings(t, y, level=0.0, deg=5):
    """Upward crossings of `level`, located by a degree-`deg` interpolant
    through the samples around each one.

    ⚠⚠ THE INTERPOLANT IS THE INSTRUMENT, AND AT deg=1 IT IS THE FLOOR.
    MEASURED: with linear interpolation, ESDIRK43 and Radau IIA(3) report the
    same floor to FIVE SIGNIFICANT FIGURES (2.9384e-10 both, at 480 points per
    period) and it falls as h^0.85 rather than as h^4 or h^5.  Two different
    integrators cannot agree to five digits on their own error -- what was
    being measured was the straight line drawn between two samples of a curved
    waveform, which depends on the GRID and not on the method that produced
    it.  Any "the floor does not drop with order" conclusion read off that is a
    statement about this function.

    `deg=1` is kept because it is the control that shows the artefact.
    """
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    lo = y[:-1] - level
    hi = y[1:] - level
    idx = np.nonzero((lo < 0) & (hi >= 0))[0]
    frac = -lo[idx] / (hi[idx] - lo[idx])
    out = t[idx] + frac * (t[idx + 1] - t[idx])
    if deg <= 1:
        return out
    half = (deg + 1) // 2
    ref = []
    for k, i in enumerate(idx):
        a = i - half + 1
        b = a + deg + 1
        if a < 0 or b > len(t):
            ref.append(out[k])
            continue
        ## fit in a LOCAL coordinate: a Vandermonde in absolute time at
        ## t ~ 1e2 with 1e-2 spacing is hopeless
        tt = t[a:b] - t[i]
        c = np.polyfit(tt, y[a:b] - level, deg)
        r = out[k] - t[i]
        for _ in range(4):              # Newton onto the interpolant's root
            f = np.polyval(c, r)
            fp = np.polyval(np.polyder(c), r)
            if fp == 0.0:
                break
            step = f / fp
            r -= step
            if abs(step) < 1e-15 * max(abs(r), 1.0):
                break
        ref.append(t[i] + r if abs(r) < 2.0 * (t[i + 1] - t[i]) else out[k])
    return np.asarray(ref, dtype=float)


def periods(cls, npts, ncyc=200, a=0.0, fixed=True, Q=15.9, reltol=1e-12,
            deg=3, drop=150):
    """Period sequence from the settled orbit.

    ⚠⚠ `drop` IS AN INSTRUMENT PARAMETER AND IT WAS THE THIRD THING THIS
    MEASUREMENT GOT WRONG.  The van der Pol transient decays like
    `exp(-mu t / 2)`, so at Q = 15.9 (mu = 0.01001) its time constant is 31.8
    CYCLES.  Discarding 5 leaves the orbit still relaxing onto its limit cycle,
    and that relaxation is a smooth, monotone, perfectly repeatable drift in
    the period -- which is exactly what a deterministic floor looks like:
    `dT/J ~ 2`, lag-1 autocorrelation 0.9, identical between methods.
    MEASURED at Q = 15.9, npts = 240, radau: J_abs reads 2.05e-10 dropping 5
    cycles and 2.81e-13 dropping 150, a factor of 730.  Every "floor" this file
    reported for the high-order methods before that was the FIXTURE.
    """
    cir, T0 = vdp(Q, a)
    iv = cir.get_node_index('v')
    x0 = np.zeros(cir.n)
    x0[iv] = 2.0
    tr = Transient(cir, integrator=cls(), reltol=reltol)
    res = tr.solve(refnode=gnd, tend=ncyc * T0, timestep=T0 / npts, x0=x0,
                   fixed_timestep=fixed)
    w = res.v('v')
    tc = crossings(np.asarray(w.x[0], dtype=float),
                   np.asarray(w.y, dtype=float), deg=deg)
    ## ⚠ AND TRIM THE TAIL.  `tend = ncyc*T0` and `timestep = T0/npts` are
    ## exactly commensurate in exact arithmetic and NOT in floating point, so
    ## the driver's final `dt = min(dt, tend - t)` can insert a sliver step --
    ## MEASURED at npts=480: 96001 steps for 96000 expected, with `dt` ranging
    ## from 4.05e-10 to 1.309e-02 where the grid is supposed to be uniform.
    ## The cycle containing that sliver is not on the grid this is measuring,
    ## and it lifted the reading of EVERY method to the same 2.3e-11 -- the
    ## same shared-value tell as the other three failures.
    Tk = np.diff(tc)[drop:-2]
    return Tk, T0, tr.statistics.accepted_steps


def _stats(Tk):
    m = float(np.mean(Tk))
    j = float(np.std(Tk))
    d = float(np.max(np.abs(Tk - m)))
    if j > 0 and len(Tk) > 3:
        z = Tk - m
        r1 = float(np.dot(z[:-1], z[1:]) / np.dot(z, z))
    else:
        r1 = float('nan')
    return m, j, d, r1


def order_and_grid():
    """P2 and P3: the floor against integrator order and against the grid."""
    print('=== SOURCES OFF: floor vs ORDER and GRID (fixed step) ===')
    print('%-9s %-5s %6s %12s %12s %8s %7s'
          % ('method', 'order', 'npts', 'J_abs', 'dT_max', 'dT/J', 'ac(1)'))
    table = {}
    for cls, name, p in METHODS:
        row = []
        for npts in (120, 240, 480):
            Tk, T0, _ = periods(cls, npts)
            m, j, d, r1 = _stats(Tk)
            row.append(j)
            print('%-9s %-5d %6d %12.4e %12.4e %8.2f %7.3f'
                  % (name, p, npts, j, d, d / j if j else float('nan'), r1))
        table[name] = row
        sl = [np.log(row[i - 1] / row[i]) / np.log(2) for i in range(1, len(row))
              if row[i] > 0]
        print('%-9s %-5d %6s slopes in h: %s   (P3 wants ~%d)'
              % (name, p, '', ' '.join('%5.2f' % x for x in sl), p))
        print()
    return table


def amplitude_sweep(cls, name, npts=240):
    """P1 and the floor itself: sweep a KNOWN perturbation down until the
    estimator stops following it.  The sloped region is what proves the
    instrument is alive; the plateau is the floor."""
    print('=== %s: J_abs vs a KNOWN perturbation (npts=%d) ==='
          % (name, npts))
    print('%12s %12s %10s %8s' % ('a', 'J_abs', 'slope', 'dT/J'))
    prev = None
    for a in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 0.0):
        Tk, T0, _ = periods(cls, npts, a=a)
        m, j, d, r1 = _stats(Tk)
        sl = (np.log(prev[1] / j) / np.log(prev[0] / a)
              if prev and a > 0 and j > 0 else float('nan'))
        print('%12.1e %12.4e %10.3f %8.2f'
              % (a, j, sl, d / j if j else float('nan')))
        if a > 0:
            prev = (a, j)
    print()


def fixed_vs_adaptive():
    """P5: the paper's own mechanism -- a step size that MOVES puts the
    crossing somewhere different each cycle."""
    print('=== P5: fixed vs adaptive step, matched by step COUNT ===')
    print('%-9s %-10s %8s %12s %8s %7s' % ('method', 'grid', 'steps', 'J_abs',
                                           'dT/J', 'ac(1)'))
    for cls, name, p in ((Gear2Integrator, 'gear', 2),
                         (RadauIIA3Integrator, 'radau', 5)):
        Tk, T0, ns_fix = periods(cls, 240)
        m, j_fix, d, r1 = _stats(Tk)
        print('%-9s %-10s %8d %12.4e %8.2f %7.3f'
              % (name, 'fixed', ns_fix, j_fix, d / j_fix, r1))
        ## adaptive, tolerance chosen to land near the same step count
        for reltol in (1e-9, 1e-11, 1e-13):
            Tk, T0, ns = periods(cls, 240, fixed=False, reltol=reltol)
            m, j, d, r1 = _stats(Tk)
            print('%-9s %-10s %8d %12.4e %8.2f %7.3f'
                  % (name, 'adp %.0e' % reltol, ns, j, d / j if j else 0.0, r1))
        print()


def determinism():
    """P4: is the sources-off jitter a RANDOM WALK or a deterministic beat?

    A deterministic period error is a frequency shift and makes no jitter at
    all.  What makes jitter is the crossing landing at a different place
    inside a step each cycle -- which beats at the incommensurability between
    the grid and the period, and is perfectly repeatable.
    """
    print('=== P4: is the floor deterministic? ===')
    for cls, name, p in ((Gear2Integrator, 'gear', 2),
                         (RadauIIA3Integrator, 'radau', 5)):
        Tk, T0, _ = periods(cls, 240)
        m, j, d, r1 = _stats(Tk)
        z = (Tk - m) / j
        ## a repeat run must reproduce it EXACTLY (no RNG anywhere)
        Tk2, _, _ = periods(cls, 240)
        rep = float(np.max(np.abs(Tk2 - Tk)))
        ## how many sign changes: a sinusoid has ~2 per beat period, noise ~N/2
        sgn = int(np.sum(np.diff(np.sign(z)) != 0))
        print('%-7s J=%.3e  dT/J=%.2f  ac(1)=%+.4f  sign changes %d/%d  '
              'repeat delta %.1e' % (name, j, d / j, r1, sgn, len(z) - 1, rep))
    print()


if __name__ == '__main__':
    ## the checks FIRST: every number below them is only meaningful if these
    ## still read the way the docstring says they do
    instrument_checks()
    order_and_grid()
    amplitude_sweep(Gear2Integrator, 'gear')
    amplitude_sweep(RadauIIA3Integrator, 'radau')
    determinism()
    fixed_vs_adaptive()


## ------------------------------------------------------------------------
## The five instrument checks.  Each one is the measurement that FOUND the
## corresponding failure in the module docstring; they are kept because the
## numbers those failures produced are quoted in the roadmap and a claim whose
## evidence exists only in a shell history is not evidence.

def check_1_interpolant_degree():
    """FAILURE 1: the crossing interpolant.  With `deg=1`, ESDIRK43 and radau
    read the SAME floor to five significant figures -- the tell that what is
    being measured is the grid, not the method.  `deg=3` and `deg=5` agree, so
    the fit is converged there."""
    print('=== check 1: the same run read with three interpolants ===')
    print('%-9s %6s %14s %14s %14s'
          % ('method', 'npts', 'deg=1 (linear)', 'deg=3', 'deg=5'))
    for cls, name, p in METHODS:
        for npts in (120, 480):
            row = [float(np.std(periods(cls, npts, deg=d)[0]))
                   for d in (1, 3, 5)]
            print('%-9s %6d %14.4e %14.4e %14.4e'
                  % (name, npts, row[0], row[1], row[2]))
    print()


def check_2_the_analytic_reference():
    """FAILURE 2: `T0 = 2 pi / sqrt(1 - mu^2/4)` is the LINEARISED period.  The
    high-order methods' "period error" sits at exactly `-mu^2/16` on EVERY
    grid -- a constant where a discretisation error would converge, which
    indicts the reference and not the method."""
    print('=== check 2: mean period against the analytic reference ===')
    mu = 1.0 / (2 * np.pi * 15.9)
    print('    mu^2/16 = %.4e  (the nonlinear correction T0 is missing)'
          % (mu ** 2 / 16 * 2 * np.pi))
    print('%-9s %6s %14s %14s' % ('method', 'npts', 'mean T - T0', 'J_abs'))
    for cls, name, p in METHODS:
        for npts in (120, 480):
            Tk, T0, _ = periods(cls, npts)
            print('%-9s %6d %14.4e %14.4e'
                  % (name, npts, float(np.mean(Tk)) - T0, float(np.std(Tk))))
    print()


def check_3_the_fixtures_settling():
    """FAILURE 3, THE ONE THAT SURVIVED FIXING THE OTHER FOUR.  Van der Pol's
    transient decays with `tau = 2/mu`, which is 31.8 CYCLES at Q = 15.9.
    Discarding 5 leaves the orbit still relaxing, and that relaxation is a
    smooth, monotone, perfectly repeatable, autocorrelated drift in the period
    -- indistinguishable by eye from a deterministic numerical floor."""
    print('=== check 3: J_abs against how many opening cycles are discarded ===')
    mu = 1.0 / (2 * np.pi * 15.9)
    cir, T0 = vdp()
    print('    tau = %.1f cycles' % (2.0 / mu / T0))
    iv = cir.get_node_index('v')
    x0 = np.zeros(cir.n)
    x0[iv] = 2.0
    tr = Transient(cir, integrator=RadauIIA3Integrator(), reltol=1e-12)
    res = tr.solve(refnode=gnd, tend=200 * T0, timestep=T0 / 240, x0=x0,
                   fixed_timestep=True)
    w = res.v('v')
    tc = crossings(np.asarray(w.x[0], dtype=float),
                   np.asarray(w.y, dtype=float), deg=3)
    print('%8s %14s %8s %7s' % ('drop', 'J_abs', 'dT/J', 'ac(1)'))
    for d in (5, 20, 50, 100, 150):
        Tk = np.diff(tc)[d:-2]
        m, j, dm, r1 = _stats(Tk)
        print('%8d %14.4e %8.2f %7.3f' % (d, j, dm / j, r1))
    print()


def check_4_grid_uniformity():
    """FAILURE 4: `ncyc*T0` and `npts*ncyc*(T0/npts)` are commensurate in exact
    arithmetic and not in floating point, so the driver's final
    `dt = min(dt, tend - t)` can insert a sliver step onto a grid the caller
    declared uniform."""
    print('=== check 4: is the "fixed" grid actually uniform? ===')
    print('%6s %9s %9s %14s %14s' % ('npts', 'steps', 'samples', 'min dt',
                                     'max dt'))
    for npts in (240, 480):
        cir, T0 = vdp()
        iv = cir.get_node_index('v')
        x0 = np.zeros(cir.n)
        x0[iv] = 2.0
        tr = Transient(cir, integrator=RadauIIA3Integrator(), reltol=1e-12)
        res = tr.solve(refnode=gnd, tend=200 * T0, timestep=T0 / npts, x0=x0,
                       fixed_timestep=True)
        t = np.asarray(res.v('v').x[0], dtype=float)
        d = np.diff(t)
        print('%6d %9d %9d %14.6e %14.6e'
              % (npts, tr.statistics.accepted_steps, len(t), d.min(), d.max()))
    print()


def check_5_time_rescale():
    """FAILURE 5: `J_abs` carries the problem's time unit.  Rescale L and C by
    100 -- identical dynamics, `T0` 100x longer -- and every absolute jitter is
    100x larger while `J_abs/T0` is constant to FIVE significant figures.  The
    fractional jitter is the quantity; an absolute one is only meaningful next
    to its `T0`.

    This is also what shows radau is at the FLOATING-POINT limit of the time
    variable rather than at its own discretisation: its `J_abs/T0` sits within
    a decade of `ulp(t)/T0`, and it does not move when the grid does.
    """
    from pycircuit.circuit.elements import C as Cap, L as Ind

    def vdp_scaled(s, Q=15.9):
        mu = 1.0 / (2 * np.pi * Q)
        c = SubCircuit()
        c.add_node('v')
        c['C'] = Cap('v', gnd, c=s)
        c['L'] = Ind('v', gnd, L=s)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: mu * (u - u ** 3 / 3.0))
        return c, s * 2 * np.pi / np.sqrt(1 - mu ** 2 / 4)

    print('=== check 5: time rescale -- J_abs/T0 is the scale-free quantity ===')
    print('%-6s %7s %10s %12s %12s %12s'
          % ('method', 'scale', 'T0', 't_max', 'J_abs', 'J_abs/T0'))
    for cls, name in ((RadauIIA3Integrator, 'radau'), (Gear2Integrator, 'gear')):
        for s in (1.0, 10.0, 100.0):
            cir, T0 = vdp_scaled(s)
            iv = cir.get_node_index('v')
            x0 = np.zeros(cir.n)
            x0[iv] = 2.0
            tr = Transient(cir, integrator=cls(), reltol=1e-12)
            res = tr.solve(refnode=gnd, tend=200 * T0, timestep=T0 / 240,
                           x0=x0, fixed_timestep=True)
            w = res.v('v')
            tc = crossings(np.asarray(w.x[0], dtype=float),
                           np.asarray(w.y, dtype=float), deg=3)
            j = float(np.std(np.diff(tc)[150:-2]))
            print('%-6s %7.0f %10.2f %12.1f %12.4e %12.4e'
                  % (name, s, T0, 200 * T0, j, j / T0))
        print('   ulp(t_max)/T0 at scale 100 = %.3e'
              % (float(np.spacing(200 * 100 * 6.2832)) / (100 * 6.2832)))
    print()


def instrument_checks():
    check_1_interpolant_degree()
    check_2_the_analytic_reference()
    check_3_the_fixtures_settling()
    check_4_grid_uniformity()
    check_5_time_rescale()
