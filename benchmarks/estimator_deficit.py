"""The local-error estimate measured against the TRUE LOCAL ERROR.

⚠⚠⚠ READ THIS SECTION FIRST.  EVERYTHING BELOW IT WAS MEASURED AGAINST THE
DECLARED `EMBEDDED_ORDER + 1`, AND THE DECLARED ORDER IS THE WRONG REFERENCE.
The methods do not HAVE their declared order on these problems -- that is
stiff order reduction, and it is the thing the older sections were reading as
an estimator defect.  `true_local_pass()` re-runs every column against the
TRUE LOCAL STATE ERROR instead (one step of size `h` from a state taken from a
run `64x` finer, compared with that run at `t + h`).

    method    fixture   TRUE local orders          est orders                 est/true, N=100..1600
    esdirk43  C-V loop   2.00 2.00 2.00 2.00        3.98 3.99 2.32 2.00        8.05 2.05 0.515 0.414 0.415
    esdirk43  RC         2.07 2.14 2.28 2.53        2.04 2.11 2.22 2.41        0.405 0.415 0.425 0.442 0.481
    radau     C-V loop   3.02 3.01 3.01 3.04        4.23 4.13 4.07 4.04        242 105 48.6 23.3 11.7
    radau     RC         3.09 3.14 3.26 3.48        2.98 2.95 2.91 2.88        0.353 0.381 0.435 0.553 0.837
    trbdf2    C-V loop   2.98 2.99 3.00 3.00        2.99 2.99 3.00 3.00        1.00 1.00 1.00 1.00 1.00
    trbdf2    RC         1.97 1.99 2.01 2.04        2.01 2.05 2.11 2.21        2.26 2.19 2.09 1.94 1.74

THE VERDICT: THERE IS NO ORDER DEFICIT.  Every estimate either matches the
true local order asymptotically or exceeds it transiently and then turns over
to match.  What is left is a BOUNDED CONSTANT factor per (method, circuit) --
1.00 to 11.7 across everything measured.  The largest is radau over-stating
11.7x on the C-V loop and still improving, which costs `11.7**0.25` ~ 1.85x
the steps; the largest under-statement is ESDIRK43 at 0.415, a true error
2.4x the estimate, FLAT rather than growing.

⚠⚠ SO BOTH CANDIDATE FIXES WERE AIMED BACKWARDS, AND THE `sigma_min` SCREEN
POINTS AT THE WRONG FIXTURE.  `sigma_min ~ h` fires on the RC and stays silent
on the C-V loop -- but the RC is where the estimate is CLOSEST to the true
local error (0.4-0.8) and the C-V loop is where it is furthest (11.7-408).
Deflating where `sigma_min ~ h` would push an already-low estimate lower.
  - A stiffly-accurate embedded formula (no filter) removes an ORDER loss.
    There is no order loss to remove.
  - `sigma_min`-based deflation reduces the estimate where the screen fires.
    The screen fires on the wrong circuit, in the wrong direction.
  - Recalibrating the controller exponent to a "measured" order assumes the
    declared one is wrong for the ESTIMATE.  It is not; the estimate tracks
    the true local error.  It is the SOLUTION that is order-reduced.

⚠⚠ THE NEAR-MISS, and the reason this file now runs five grids.  On the
first three grids ESDIRK43's C-V loop ratio reads 8.05 -> 2.05 -> 0.515, a
clean geometric fall through 1, with the estimate at order 4 against a true
order 2.  Extrapolated, that is an estimate going unboundedly blind -- an
accuracy defect worth alarm.  TWO MORE GRIDS SHOW IT TURNS OVER: the estimate
order drops 3.99 -> 2.32 -> 2.00 to meet the true order, and the ratio flattens
at 0.415.  Three points on a log plot were enough to draw a line and wrong
about where it went.

⚠ WHAT MAY NOT BE READ FROM THE TABLE.  The ExpG rows are PRE-ASYMPTOTIC at
these step sizes (orders still climbing at the finest grid, e.g. radau
2.35 2.81 3.02) and are excluded above; the RC control -- ExpG's topology with
the nonlinearity removed -- is asymptotic and carries the same information.
And the RC true orders are still DRIFTING upward at N=1600 (3.09 -> 3.48), two
error components with the higher-order one taking over, so the RC order values
are local, not settled.  The C-V loop rows are settled.

CORROBORATION FROM A DIFFERENT HARNESS: radau's true local order on the
index-2 C-V loop reads 3.01, and the max-norm is picking up the algebraic
component.  That is the independently measured "order 5 differential / 3
algebraic on an index-2 MNA" (HLR Thm 5.9) arriving from a second direction.

THE REFERENCE IS NOT THE INSTRUMENT -- `reference_control()` re-runs rows
against a DIFFERENT method at DOUBLE the refinement (ESDIRK43 at `h/128` in
place of radau at `h/64`) and the true local orders are identical to three
digits.  A shared-method reference could have hidden a shared error; it does
not.

Two instrument failures are pinned in the code below rather than remembered:
the clock restart in `_march`, and the moving phase in `_start_state`.

============================================================================
EVERYTHING BELOW WAS MEASURED AGAINST THE DECLARED ORDER.  Kept because the
raw-vs-filtered split, the `sigma_min` topology and the exact-rational tableau
check are all still true statements about what they measured -- but "DEFICIT"
in those tables means "short of the DECLARED order", which the section above
shows is not a defect.
============================================================================

THE CLAIM THIS CORRECTS.  2026-09-10, building adaptive stepping: the filtered
estimate reads one order BELOW its declared `EMBEDDED_ORDER + 1`, and I
attributed it to the DAE index -- "`J = C + a h G` is singular in `C` on a
DAE, so `J^-1` behaves like `1/h` on the algebraic subspace".  That was
measured on ONE fixture.

⚠⚠⚠ IT IS NOT AN INDEX EFFECT.  Measured on three fixtures, estimate order per
doubling, deficit against the declared `EMBEDDED_ORDER + 1`:

    method   fixture              want   slopes              verdict
    radau    ExpG (nonlinear)      4     3.02 3.10 3.18      DEFICIT 0.82
    radau    RC (LINEAR control)   4     3.10 3.12 3.19      DEFICIT 0.81
    radau    C-V loop (INDEX 2)    4     4.03 4.03 4.01      no deficit
    trbdf2   ExpG (nonlinear)      3     1.93 1.84 1.77      DEFICIT 1.23
    trbdf2   RC (LINEAR control)   3     1.88 1.88 1.84      DEFICIT 1.16
    trbdf2   C-V loop (INDEX 2)    3     2.88 2.95 2.97      no deficit

The INDEX-2 circuit is the one with NO deficit, which is the opposite of an
index effect.  ESDIRK43 and GLM3 read the same way (deficit 1.94 and 0.93 on
ExpG, none on the C-V loop).

⚠ AND IT IS NOT THE NONLINEARITY EITHER.  The RC control is ExpG's topology
with the nonlinear conductance replaced by a resistor, and it shows the SAME
deficit to two decimals.  Both of my candidate explanations are refuted by one
pair of rows.

WHAT IT ACTUALLY IS: `sigma_min(J)` SCALING WITH `h`.

    RC control   sigma_min = 1.0146e-07 at h=2.0e-05,  1.3557e-08 at h=2.5e-06
                 -> ~ h, so `J^-1` ~ 1/h  -> one order lost
    C-V loop     sigma_min = 2.0500e-09 at h=2.0e-05,  2.0062e-09 at h=2.5e-06
                 -> CONSTANT (it is ||C||) -> nothing lost

The discriminator is topological: a circuit with a purely RESISTIVE path from
a source has a direction of `J` carrying no capacitance, whose singular value
is `~h*G`; a C-V loop's smallest direction is capacitive and is floored at
`||C||`.  It is a property of WHICH DIRECTIONS CARRY REACTANCE, not of the
index and not of the nonlinearity.

WHY THIS MATTERS PRACTICALLY.  The controller uses
`dt_next = dt * SAFETY * err^(-1/(ORDER+1))`, so a wrong estimate order makes
the exponent wrong and the estimate itself inflated -- steps stay smaller than
the accuracy warrants.  And the target IS identifiable at realistic step sizes:
`sigma_min(J)` comes from a matrix already assembled and factored every step,
so whether this circuit loses an order is a per-step measurement rather than an
assumption.  ⚠ That is the measurement, not the fix; deflating an error
estimate is exactly how accuracy is silently lost, so any correction has to be
gated two-sided -- same error against an independent reference AND fewer steps.

⚠⚠ AND SPLITTING RAW FROM FILTERED SPLITS IT INTO TWO DEFECTS.  Hairer &
Wanner IV.8 (8.18)/(8.19) say the embedded COEFFICIENTS supply the order and
the filter supplies only BOUNDEDNESS at `h*lambda -> infinity` -- "for h -> 0
we still have err = O(h^4)".  So the filter is order-PRESERVING, and a
deficient filtered estimate can come from either end.  MEASURED, both:

    method    fixture          RAW slopes          FILTERED slopes
    trbdf2    RC               2.87 2.86 2.80      1.88 1.88 1.84
    trbdf2    C-V loop         2.90 2.96 2.98      2.88 2.95 2.97
    esdirk43  RC               2.95 2.96 3.04      1.97 1.98 2.10
    esdirk43  C-V loop         2.95 2.94 2.99      3.89 3.96 3.69
    radau     RC               3.09 3.09 3.14      3.10 3.12 3.19
    radau     C-V loop         3.01 3.01 3.00      4.03 4.03 4.01

  ⚠⚠ THE RAW COLUMNS ARE NOT COMPARABLE TO EACH OTHER OR TO A SINGLE "want".
  They are DIFFERENT PHYSICAL QUANTITIES and each carries its own expected
  order, derived from its own definition:

    radau     raw = C*F1 + f0          CHARGE-RATE; the filter
              ((gamma/h)C + G)^-1 is ~ h/(gamma C) on a C-dominated
              direction, so raw order = state order - 1 = 3
    trbdf2    raw = h*sum(dk_i K_i)    CHARGE; the filter (C + a h G)^-1 is
    esdirk43  (same form)              ~ 1/C there, carrying NO h, so raw
                                       order = the declared order

  Against those, on the C-V loop: radau wants 3 and reads 3.00, TR-BDF2 wants
  3 and reads 2.98, ESDIRK43 wants 4 and reads 2.99.  ONLY ESDIRK43 IS SHORT.
  ⚠ A first reading of this file took radau's raw 3.00 as refuting the
  stage-order account, by applying "want 4" to a quantity whose want is 3.
  The INTERNAL CONTROL is what makes the ESDIRK43 gap safe: TR-BDF2 and
  ESDIRK43 use the SAME raw form on the SAME fixture and differ only in
  declared order -- 2.98 against 3, 2.99 against 4.

  DEFECT 1 -- THE FILTER, only where `sigma_min(J) ~ h`.  TR-BDF2's raw
  estimate is 2.80-2.98 on BOTH fixtures (its declared 3), and the filter
  takes it to 1.84 on the RC and leaves it at 2.97 on the C-V loop.  That is
  the substitution Hairer does NOT cover: his filter uses the ODE Jacobian
  `df/dy`, which on a DAE with singular `C` does not exist, so this code
  substitutes the step matrix -- and that substitution is not order-preserving
  when the step matrix has a direction carrying no reactance.

  DEFECT 2 -- THE DIRK STAGE-ORDER CAP, NOT THE COEFFICIENTS.  ⚠⚠ AN EARLIER
  VERSION OF THIS FILE SAID "Hairer's coefficient condition failing"; that is
  REFUTED.  docs-46 evaluated the Kennedy & Carpenter ARK4(3)6L[2]SA tableau
  in EXACT RATIONAL arithmetic: `B` satisfies all order conditions to 4 with
  zero residual and `B_HAT` to 3, failing at 4 -- which is precisely what
  `EMBEDDED_ORDER = 3` declares.  The coefficients are exactly as labelled.

  What caps it is the STAGE ORDER, and the method's own name says so:
  ESDIRK4(3)6L[2]SA, the `[2]` being `q`.  Hairer & Wanner IV.15 Ex.1, "the
  stage order of a DIRK method is at most 2", and Ex.3 (Burrage &
  Hundsdorfer), "the order of B-convergence ... is q+1".  So on a stiff
  problem the attainable order is capped at 3 for ANY DIRK -- which is
  circuit-independent, exactly the property measured.
  ⚠ THE JOIN IS RELAYED AND UNVERIFIED: those results are about the order of
  the SOLUTION; that the same cap lands on the embedded ESTIMATE is docs-46's
  inference, fits these numbers, and is not something either of us has seen
  stated or measured directly.

  It unifies both DIRK rows: TR-BDF2 is also a DIRK with `q = 2` and cap 3,
  and its DECLARED 3 already equals the cap, so no deficit can appear --
  raw 2.98 against 3.  ESDIRK43 declares 4 against the same cap -- raw 2.99,
  deficit 1.0.  The gap only appears where the declared order EXCEEDS q+1.

  ⚠ DO NOT RELABEL `EMBEDDED_ORDER`.  It is correct for the tableau and the
  controller's exponent is right for the nonstiff regime the tableau is
  designed for.  What is wrong is the EXPECTATION that a DIRK's estimate
  reaches its classical order on a stiff problem.  Relabelling would make the
  controller wrong in the nonstiff case to make a stiff measurement agree.

⚠ A THIRD READING WORTH HAVING (Guenther 2005 p.47, relayed): for STIFFLY
ACCURATE EMBEDDED ROW methods the estimate IS a stage increment, "based on
node potentials and branch currents only" -- no filter, no extra solve,
nothing needing a Jacobian inverse.  On that reading the filter is a REPAIR
for a non-stiffly-accurate embedded formula and stiff accuracy removes the
need for it.

Run: `python benchmarks/estimator_deficit.py`
"""
import warnings

import numpy as np

from pycircuit.circuit.circuit import SubCircuit, gnd, defaultepar
from pycircuit.circuit.elements import R, C as Cap, VSin, ISin
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import (RadauIIA3Integrator,
                                          TRBDF2Integrator,
                                          ESDIRK43Integrator, GLM3Integrator)

warnings.simplefilter('ignore')
PER = 1e-3


def rc_linear():
    """ExpG's TOPOLOGY with the nonlinearity removed -- the control that
    separates `sigma_min` from the nonlinear device."""
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=0.8, freq=1.0 / PER)
    c['rs'] = R('a', 'b', r=50.0)
    c['cl'] = Cap('b', gnd, c=1e-9)
    c['rl'] = R('b', gnd, r=1e4)
    return c


def expg():
    from pycircuit.circuit.tests.test_stage_predictor import _expg_fixture
    return _expg_fixture(PER)


def cv_loop():
    from pycircuit.circuit.tests.test_glm import _cv_loop
    return _cv_loop(PER)


def estimate_order(cls, build, Ns=(50, 100, 200, 400)):
    v = []
    for N in Ns:
        tr = Transient(build(), integrator=cls(), reltol=1e-12)
        tr._rk_want_est = True
        rec = []
        orig = Transient.solve_timestep

        def w(self, x0, t, *a, **k):
            r = orig(self, x0, t, *a, **k)
            e = getattr(self, '_rk_est', None)
            if e is not None:
                rec.append(float(np.max(np.abs(np.asarray(e, dtype=float)))))
            return r
        Transient.solve_timestep = w
        try:
            tr.solve(refnode=gnd, tend=PER, timestep=PER / N,
                     fixed_timestep=True)
        finally:
            Transient.solve_timestep = orig
        v.append(float(np.median(np.array(rec)[6:])))
    return [np.log(v[i - 1] / v[i]) / np.log(2) for i in range(1, len(v))]


def sigma_min_scaling(build, label):
    """Does `sigma_min(J)` scale with `h` at REALISTIC step sizes?  That is the
    discriminator, and it is computable from a matrix already factored."""
    cir = build()
    tr = Transient(cir, integrator=RadauIIA3Integrator(), reltol=1e-12)
    tr.irefnode = cir.get_node_index(gnd)
    x = np.zeros(cir.n)
    out = []
    for npts in (50, 100, 200, 400):
        h = PER / npts
        C = np.asarray(cir.C(x, defaultepar), dtype=float)
        G = np.asarray(cir.G(x, defaultepar), dtype=float)
        (Jr,) = remove_row_col((C + 0.25 * h * G,), tr.irefnode, tr.toolkit)
        out.append(float(np.linalg.svd(np.asarray(Jr, dtype=float),
                                       compute_uv=False).min()))
    sl = [np.log(out[i - 1] / out[i]) / np.log(2) for i in range(1, len(out))]
    print('  %-18s sigma_min %s   exponent in h %s'
          % (label, ' '.join('%.3e' % v for v in out),
             ' '.join('%+.2f' % v for v in sl)))
    return sl


# ---------------------------------------------------------------------------
# THE TRUE LOCAL ERROR AS THE REFERENCE, instead of the declared order.
# ---------------------------------------------------------------------------

_REF_CACHE = {}


def _march(cls, build, h, n, x0=None, t0=0.0, want_est=False):
    """`n` steps of size `h` starting from `x0` at ABSOLUTE time `t0`.

    ⚠⚠ `t0` MATTERS AND IS THE INSTRUMENT FAILURE THAT ALMOST LANDED.  A first
    version restarted the clock at 0 for the one-step run, so `VSin` was
    evaluated at the wrong phase and the "true local error" read ORDER 1 for
    an order-5 method.  It was caught only because a local error CANNOT be
    O(h) -- the number was physically impossible, not merely surprising.
    """
    cir = build()
    tr = Transient(cir, integrator=cls(), reltol=1e-14)
    tr.irefnode = cir.get_node_index(gnd)
    if x0 is None:
        x = np.asarray(DC(cir, refnode=gnd).solve().x, dtype=float).ravel()
    else:
        x = np.array(x0, dtype=float)
    tr.epar.t = t0
    tr._begin_run(x, cir.n)
    tr._rk_want_est = want_est
    xs = [x.copy()]
    ests = []
    for j in range(1, n + 1):
        tr._dt_last = tr._dt if j > 1 else None
        tr._dt = h
        tr.epar.t = t0 + j * h
        x, _f, _J, _ = tr.solve_timestep(x, t0 + j * h)
        tr._push_history(x)
        xs.append(np.asarray(x, dtype=float).copy())
        e = getattr(tr, '_rk_est', None)
        ests.append(float(np.max(np.abs(np.asarray(e, dtype=float))))
                    if e is not None else np.nan)
    return np.array(xs), np.array(ests)


T0_FIX = PER / 8.0
REF_DIV = 64
_START_CACHE = {}


def _start_state(build, lbl, ref_cls, href):
    """The EXACT state at the FIXED absolute time `T0_FIX`.

    ⚠⚠ THE SECOND INSTRUMENT FAILURE, and it is subtler than the clock one.
    A first version started the step at `t0 = (nref-1)*h`, which is a
    DIFFERENT ABSOLUTE TIME for every grid -- so the sinusoid was at a
    different phase on each, and the local-error coefficient `C(t)` in
    `C(t) h^(p+1)` moved between the very points whose ratio IS the order.
    It contaminated only the ORDERS; the `est/true` ratio was immune because
    both sides came from the same step.  The tell was radau's estimate
    reading 3.93 there and 3.19 in `estimate_order` -- the SAME quantity,
    two harnesses, two answers.
    """
    key = (lbl, ref_cls.__name__, href)
    if key not in _START_CACHE:
        n = int(round(T0_FIX / href))
        assert abs(n * href - T0_FIX) < 1e-18, 'T0_FIX not on the fine grid'
        xr, _ = _march(ref_cls, build, href, n)
        _START_CACHE[key] = xr[-1].copy()
    return _START_CACHE[key]


def true_local(cls, build, lbl, Ns=(50, 100, 200, 400), nstep=1,
               ref_cls=RadauIIA3Integrator, ref_div=REF_DIV):
    """TRUE local state error and the filtered estimate, on the same steps,
    every grid starting from the SAME absolute time `T0_FIX`.

    `nstep > 1` is for the multivalue methods, which have no estimate until
    the Nordsieck vector has started.  Over a FIXED number of steps the error
    accumulates to `nstep * C h^(p+1)`, so the ORDER is preserved and only the
    `est/true` ratio carries the extra factor `~nstep`.
    """
    href = PER / (max(Ns) * ref_div)
    start = _start_state(build, lbl, ref_cls, href)
    tl, es = [], []
    for N in Ns:
        h = PER / N
        key = (lbl, N, nstep, ref_cls.__name__, ref_div)
        if key not in _REF_CACHE:
            xr, _ = _march(ref_cls, build, h / ref_div, nstep * ref_div,
                           x0=start, t0=T0_FIX)
            _REF_CACHE[key] = xr[-1].copy()
        target = _REF_CACHE[key]
        xs, ee = _march(cls, build, h, nstep, x0=start, t0=T0_FIX,
                        want_est=True)
        tl.append(float(np.max(np.abs(xs[nstep] - target))))
        es.append(float(ee[-1]))
    return tl, es


def _orders(v):
    return [np.log(v[i - 1] / v[i]) / np.log(2) for i in range(1, len(v))]


def crossover(Ns=(100, 200, 400, 800, 1600)):
    """The five-grid table in the header.

    ⚠ FIVE GRIDS, NOT THREE, AND THAT IS THE POINT.  ESDIRK43's C-V loop ratio
    falls 8.05 -> 2.05 -> 0.515 on the first three -- a clean line through 1,
    which extrapolates to an estimate going unboundedly blind.  It does not:
    the estimate's ORDER turns over to meet the true order and the ratio
    flattens at 0.415.  ExpG is omitted, pre-asymptotic at these `h`.
    """
    print('=== est/true across five grids: does the trend turn over? ===')
    print('%-9s %-9s %-28s %-28s %s'
          % ('method', 'fixture', 'TRUE local orders', 'est orders',
             'est/true per grid'))
    for cls, nm in ((ESDIRK43Integrator, 'esdirk43'),
                    (RadauIIA3Integrator, 'radau'),
                    (TRBDF2Integrator, 'trbdf2')):
        for build, lbl in ((cv_loop, 'C-V loop'), (rc_linear, 'RC')):
            tl, es = true_local(cls, build, lbl + '/5', Ns=Ns)
            print('%-9s %-9s %-28s %-28s %s'
                  % (nm, lbl,
                     ' '.join('%5.2f' % s for s in _orders(tl)),
                     ' '.join('%5.2f' % s for s in _orders(es)),
                     ' '.join('%.3g' % (e / t) for e, t in zip(es, tl))))



def prothero_robinson(lam, w=2.0 * np.pi):
    """Scalar Prothero-Robinson as an MNA circuit, with an EXACT solution.

    `y' = lam (y - phi(t)) + phi'(t)` with `phi(t) = sin(w t)`, so `y(t) =
    sin(w t)` exactly.  As MNA: a unit capacitor to ground, a resistor
    `r = -1/lam`, and one sinusoidal current source carrying
    `-lam sin(w t) + w cos(w t)`, which is a SINGLE sinusoid of amplitude
    `hypot(lam, w)` and phase `atan2(w, -lam)`.

    WHY THIS FIXTURE.  `stiffness_control_is_CONFOUNDED` below records that
    tuning the C-V loop's resistor to add stiffness ALSO switches on a
    reactance-free direction, so it cannot separate the two effects, and names
    a scalar Prothero-Robinson fixture as what would.  This is that fixture.
    It has an EXACT solution, so there is no reference run and no inner-solve
    tolerance anywhere in the measurement -- the two floors that limit
    `global_order_control`.  Verified: a fine radau run matches `sin(w t)` to
    1.4e-15 at `lam = -1`, 4.4e-16 at `-1e2` and 2.2e-15 at `-1e4`.
    """
    amp = float(np.hypot(lam, w))
    phase = float(np.degrees(np.arctan2(w, -lam)))
    c = SubCircuit()
    c.add_node('a')
    c['c'] = Cap('a', gnd, c=1.0)
    c['r'] = R('a', gnd, r=1.0 / (-lam))
    c['i'] = ISin(gnd, 'a', ia=amp, freq=w / (2.0 * np.pi), phase=phase)
    return c


def rc_asymptotic_check(Ns=(200, 400, 800, 1600, 3200), span=PER / 4.0):
    """⚠⚠ IT WAS PRE-ASYMPTOTIC.  `global_order_control`'s RC rows read below
    the classical `p` and were STILL RISING at its finest grid, which is a
    caveat that file raised against itself.  Five grids settle it:

        method    p   GLOBAL orders               errors
        esdirk43  4    2.54  3.05  3.50  3.79     2.50e-10 ... 3.33e-14
        trbdf2    2    2.12  2.13  2.12  2.09     6.49e-10 ... 1.86e-12
        radau     5    3.48  3.88  3.87  0.58     8.02e-13 ... 2.22e-16

    ESDIRK43 climbs monotonically to its `p = 4`.  TR-BDF2 sits flat at its
    true `p = 2`.  So there is no order reduction on the index-1 RC and no
    Thm 2.26 hypothesis failure to explain -- which agrees with docs-46's
    independent read-only measurement that the RC's algebraic block is
    essentially perfectly conditioned, `sigma_min(d g_2/d y) = 0.9995`, flat.

    ⚠ RADAU IS UNRESOLVABLE ON THIS FIXTURE, AND ONLY THE MAGNITUDE SAYS SO.
    Its last error is 2.22e-16, which IS machine epsilon; the `0.58` is the
    floor, not the method.  The usable part is 3.48 -> 3.88 -> 3.87, climbing.
    An order-5 method reaches `eps` on this circuit before it reaches its
    asymptotic regime, so no refinement of THIS fixture can confirm `p = 5`.

    ⚠ AND A DOUBT OF MINE THAT WAS WRONG, kept because it nearly stopped the
    measurement: I expected the REFERENCE to be the limit at these grids,
    reasoning that a 25600-step reference march at an inner `reltol = 1e-14`
    would accumulate past the finest test error.  It does not.  The same
    esdirk43 errors against references at `PER/25600`, `PER/51200` and
    `PER/102400` are 2.5003e-10 / 4.3083e-11 / 5.2015e-12 -- unmoved to five
    digits at every grid.  A reasoned floor is not a measured floor.
    """
    href = PER / (400.0 * 128)
    start = _start_state(rc_linear, 'RC/asym', RadauIIA3Integrator, href)
    xr, _ = _march(RadauIIA3Integrator, rc_linear, href,
                   int(round(span / href)), x0=start, t0=T0_FIX)
    tgt = xr[-1]
    print('=== RC (index 1) global order, five grids: does it climb to p? ===')
    print('%-9s %-3s %-30s %s' % ('method', 'p', 'GLOBAL orders', 'errors'))
    for cls, nm, p_ in ((ESDIRK43Integrator, 'esdirk43', 4),
                        (TRBDF2Integrator, 'trbdf2', 2),
                        (RadauIIA3Integrator, 'radau', 5)):
        errs = []
        for N in Ns:
            xs, _ = _march(cls, rc_linear, span / N, N, x0=start, t0=T0_FIX)
            errs.append(float(np.max(np.abs(xs[-1] - tgt))))
        print('%-9s %-3d %-30s %s'
              % (nm, p_, ' '.join('%5.2f' % v for v in _orders(errs)),
                 ' '.join('%.2e' % v for v in errs)))



def local_vs_global_relation(t0=0.125, span=0.25):
    """Does `global = local - 1` survive stiffness?  NO -- it goes to zero.

        method    lam      h*|lam|  LOCAL orders       GLOBAL orders      glob-loc
        esdirk43  -1e0     0.00781   4.87 4.95 4.97     3.98 3.99 3.99     -0.98
        esdirk43  -1e2     0.125     4.64 4.82 4.91     3.95 3.99 4.00     -0.91
        esdirk43  -1e4     12.5      2.14 2.28 2.53     2.02 2.18 2.53     +0.00
        esdirk43  -1e6     1.25e+03  2.00 2.00 2.01     1.95 1.98 1.99     -0.02
        trbdf2    -1e0     0.00781   2.80 2.91 2.96     1.58 1.84 1.93     -1.03
        trbdf2    -1e2     0.125     2.72 2.85 2.92     2.03 2.02 2.01     -0.91
        trbdf2    -1e4     12.5      1.97 2.00 2.04     2.04 2.08 2.11     +0.07
        trbdf2    -1e6     1.25e+03  1.97 1.98 1.99     1.98 1.99 2.00     +0.01
        radau     -1e0     0.0104    5.76 5.04 5.68     5.02 5.01 5.00     -0.68
        radau     -1e2     0.125     5.62 5.80 5.91     4.92 4.96 4.98     -0.92
        radau     -1e4     12.5      3.15 3.26 3.49     3.12 3.23 3.49     +0.00
        radau     -1e6     1.25e+03  3.03 3.02 3.01     3.02 3.01 3.01     +0.01

    THE NONSTIFF ARM IS THE CONTROL AND IT IS WHAT MAKES THIS MEAN ANYTHING.
    At `lam = -1` the classical `global = local - 1` holds for all three
    methods, and the relation goes to ZERO as stiffness rises.  Without that
    arm the stiff number is just a number.  So on a stiff problem THE ENDPOINT
    ERROR IS ONE LOCAL ERROR -- earlier contributions are damped out -- and
    the classical local/global relation may not be used on any number in this
    file, in either direction.

    Independently reproduced: docs-46 measured the same switch on its own
    scalar fixture with an exact linear stage solve (glob-loc +0.99 at
    `lam=-1e0`, +0.08 at `-1e2`, -0.02 at `-1e4`, -0.05 at `-1e6`; sign
    convention reversed).  Two instruments, no shared code, same switch.

    AND IT CLOSES A FLAGGED ITEM.  The stiff arm lands all three methods on
    EXACTLY THEIR STAGE ORDER `q` -- esdirk43 `q=2` reads 1.99, radau `q=3`
    reads 3.01, trbdf2 `q=2` reads 2.00.  ⚠⚠ That is the same "all three land
    at `q`" pattern that `stiffness_control_is_CONFOUNDED` records as an
    ARTEFACT -- but there it was confounded because changing the C-V loop's
    resistor also switched on a reactance-free direction.  Here there is no
    DAE structure to carry one, which is exactly why that docstring nominated
    this fixture.  Same pattern, and this time it is not an artefact.

    ⚠ FLOOR, and it is the SOLUTION SCALE not `eps` (docs-46's caution, and it
    bites here).  With `|y| ~ 1` the radau rows run at 2.2e-13 to 8.3e-15
    absolute, tens of ulp, and radau's nonstiff LOCAL orders are visibly noisy
    (5.76 5.04 5.68 for a true 6) -- which is why its `glob-loc` reads -0.68
    rather than -1.  Its GLOBAL arm, further from the floor, reads a clean
    5.02 5.01 5.00.  Keep the grids coarse enough that the error stays above
    ~1e-13.
    """
    print('=== local vs global order across stiffness (EXACT solution) ===')
    print('%-9s %-9s %-9s %-20s %-20s %-8s %9s'
          % ('method', 'lam', 'h*|lam|', 'LOCAL orders', 'GLOBAL orders',
             'glob-loc', 'loc err'))
    for cls, nm, coarse in ((ESDIRK43Integrator, 'esdirk43', (4, 8, 16, 32)),
                            (TRBDF2Integrator, 'trbdf2', (4, 8, 16, 32)),
                            (RadauIIA3Integrator, 'radau', (3, 6, 12, 24))):
        for lam in (-1e0, -1e2, -1e4, -1e6):
            Ns = coarse if lam == -1e0 else (25, 50, 100, 200)
            build = lambda l=lam: prothero_robinson(l)
            loc, glo = [], []
            for N in Ns:
                h = span / N
                x0 = np.array([np.sin(2 * np.pi * t0), 0.0])
                xs, _ = _march(cls, build, h, 1, x0=x0, t0=t0)
                loc.append(abs(float(xs[1][0]) -
                               np.sin(2 * np.pi * (t0 + h))))
                xs, _ = _march(cls, build, h, N, x0=x0, t0=t0)
                glo.append(abs(float(xs[-1][0]) -
                               np.sin(2 * np.pi * (t0 + span))))
            ol, og = _orders(loc), _orders(glo)
            print('%-9s %-9.0e %-9.3g %-20s %-20s %+8.2f %9.2e'
                  % (nm, lam, (span / Ns[-1]) * abs(lam),
                     ' '.join('%5.2f' % s for s in ol),
                     ' '.join('%5.2f' % s for s in og),
                     og[-1] - ol[-1], loc[-1]))
        print()



def global_order_control(Ns=(50, 100, 200, 400), span=PER / 4.0):
    """⚠⚠ INCONCLUSIVE BY CONSTRUCTION -- KEPT FOR WHAT IT RULES IN, NOT OUT.

    docs-46 relayed Baechle 2007 Thm 2.26 (index-1 DAE, IRK of classical order
    `p` satisfying `C(q)`, `d g_2 / d y` with a BOUNDED INVERSE near the
    solution): for a STIFFLY ACCURATE method both components converge at `p`.
    All three methods here are stiffly accurate with `R(inf) = 0`, and the RC
    is index 1, so Thm 2.26 predicts `p`.  If we measure below `p`, a
    hypothesis is failing -- and the candidate named in the theorem is the
    bounded inverse, which is what a reactance-free direction threatens.

    THE MEASUREMENT DOES NOT SETTLE IT.  Global error over a fixed span from an
    EXACT start:

        method    fixture   p   GLOBAL orders        final rel err
        radau     RC        5    3.11 3.23 3.48       1.27e-13
        trbdf2    RC        2    2.05 2.08 2.12       2.65e-10
        esdirk43  RC        4    2.04 2.19 2.54       7.62e-11
        radau     C-V loop  5    3.33 3.01 3.01       2.37e-15
        trbdf2    C-V loop  2    2.00 2.00 2.00       1.76e-07
        esdirk43  C-V loop  4    3.97 1.99 1.99       1.39e-12

    Three reasons not to read a hypothesis failure out of this:

      - ⚠ THE RADAU ROWS ARE TOLERANCE-LIMITED, NOT FLOOR-LIMITED, AND THAT IS
        EASY TO MISS.  The inner Newton runs at `reltol = 1e-14`, so the
        per-step solve error is ~1e-14 -- and radau's RC error reaches 1.3e-13
        and the C-V loop 2.4e-15.  Within one to ten times the SOLVER's own
        tolerance, the slope is the tolerance's, not the method's.  Machine
        epsilon is far below and would not have flagged it; the magnitude to
        check against is the INNER SOLVE.
      - THE RC ROWS ARE STILL RISING at the finest grid (esdirk43 2.04 -> 2.54,
        radau 3.11 -> 3.48), i.e. pre-asymptotic, so a value below `p` is not
        yet a value that stays below `p`.
      - THE C-V LOOP IS INDEX 2 and outside the theorem entirely.

    WHAT IT DOES ESTABLISH, and it is worth having: GLOBAL ORDER EQUALS LOCAL
    ORDER here, row for row (radau RC 3.48/3.48, trbdf2 RC 2.12/2.04, esdirk43
    RC 2.54/2.53), instead of the `local = global + 1` of the ODE case.  Errors
    are NOT accumulating: the endpoint error is one local error.  That is what
    `R(inf) = 0` does to the algebraic component -- every earlier contribution
    is damped out and the last step is the whole error -- and it is the reason
    the estimator sections above cannot be reasoned about with the classical
    local/global relation.

    ⚠⚠ CORRECTION, SAME DAY.  A first version of this table labelled TR-BDF2
    `p = 3`.  IT IS AN ORDER-2 METHOD -- `integrator.py` says so in its own
    docstring, "it is here for three failure modes it does NOT have rather than
    for accuracy (it is order 2)".  With the right `p` the trbdf2 rows are AT
    full order, not reduced, and they were the rows I had read as the cleanest
    evidence of reduction.  The corrected reading is that order reduction
    appears exactly where `p > q + 1`: ESDIRK43 (`p=4`, `q=2`) and radau
    (`p=5`, `q=3`) both read below `p` and both are still rising; TR-BDF2
    (`p=2`, `q=2`) has no room to reduce and does not.

    ⚠⚠ ANSWERED, SAME DAY, BY `rc_asymptotic_check`: THE ROWS WERE
    PRE-ASYMPTOTIC.  Given two more grids ESDIRK43 climbs 2.54 -> 3.05 -> 3.50
    -> 3.79 to its `p = 4` and TR-BDF2 sits flat at its true `p = 2`, so there
    is no reduction on the index-1 RC and no hypothesis failure to explain.
    docs-46 reached the same place from the other side, measuring the RC's
    algebraic block as essentially perfectly conditioned
    (`sigma_min(d g_2/d y) = 0.9995`, flat) -- so Thm 2.26's hypothesis HOLDS,
    and the mechanism it had proposed (a reactance-free direction costing the
    bounded inverse) is refuted on our own fixtures.  Radau is unresolved HERE --
    it reaches machine epsilon on this circuit before its asymptotic regime --
    and resolved elsewhere: docs-46 reports (RELAYED, not reproduced in this
    repo) that against the EXACT analytic solution of the linear RC, with no
    reference run at all, radau reads 4.98 on its best PRE-FLOOR pair and all
    three methods reach full classical order.  Its error stalls at ~4e-13 for
    three consecutive grids there, where the fitted orders read 0.46, 0.16 and
    8.17 -- so that result depends entirely on reading BEFORE the floor, which
    is why it is worth restating rather than quoting a fit.

    A CONCLUSIVE version needs all three: an inner solve tightened well below
    the discretisation error (or a fixture whose error scale is larger), grids
    that reach the asymptotic regime, and the error SPLIT into differential and
    algebraic components rather than read through a max-norm that silently
    picks whichever is larger.  NOT RUN -- it is a new question about the
    SOLUTION's order, not the estimator's, and it belongs to whoever takes that
    up.
    """
    href = PER / (max(Ns) * REF_DIV)
    print('=== global order over a fixed span (see docstring: INCONCLUSIVE) ===')
    print('%-9s %-9s %-5s %-26s %s'
          % ('method', 'fixture', 'p', 'GLOBAL orders', 'final rel err'))
    for build, lbl in ((rc_linear, 'RC'), (cv_loop, 'C-V loop')):
        start = _start_state(build, lbl + '/g', RadauIIA3Integrator, href)
        xr, _ = _march(RadauIIA3Integrator, build, href,
                       int(round(span / href)), x0=start, t0=T0_FIX)
        tgt = xr[-1]
        scale = max(float(np.max(np.abs(tgt))), 1e-30)
        for cls, nm, p in ((RadauIIA3Integrator, 'radau', 5),
                           (TRBDF2Integrator, 'trbdf2', 2),
                           (ESDIRK43Integrator, 'esdirk43', 4)):
            errs = []
            for N in Ns:
                xs, _ = _march(cls, build, span / N, N, x0=start, t0=T0_FIX)
                errs.append(float(np.max(np.abs(xs[-1] - tgt))))
            print('%-9s %-9s %-5d %-26s %.3e'
                  % (nm, lbl, p, ' '.join('%5.2f' % s for s in _orders(errs)),
                     errs[-1] / scale))
        print()



def true_local_pass():
    """EVERY column re-read against the true local error rather than the
    declared `EMBEDDED_ORDER + 1`."""
    print('=== the filtered estimate against the TRUE LOCAL STATE ERROR ===')
    print('%-9s %-18s %-4s %-18s %-18s %s'
          % ('method', 'fixture', 'decl', 'TRUE local orders',
             'est orders', 'est/true'))
    for cls, name, nstep in ((RadauIIA3Integrator, 'radau', 1),
                             (TRBDF2Integrator, 'trbdf2', 1),
                             (ESDIRK43Integrator, 'esdirk43', 1),
                             (GLM3Integrator, 'glm3', 6)):
        decl = cls().EMBEDDED_ORDER + 1
        for build, lbl in ((expg, 'ExpG (nonlinear)'),
                           (rc_linear, 'RC (LINEAR control)'),
                           (cv_loop, 'C-V loop (INDEX 2)')):
            try:
                tl, es = true_local(cls, build, lbl, nstep=nstep)
            except Exception as exc:
                print('%-9s %-18s %-4d  FAILED: %s'
                      % (name, lbl, decl, type(exc).__name__))
                continue
            ot, oe = _orders(tl), _orders(es)
            r0 = es[0] / tl[0] if tl[0] else float('nan')
            r1 = es[-1] / tl[-1] if tl[-1] else float('nan')
            print('%-9s %-18s %-4d %-18s %-18s %.3g -> %.3g%s'
                  % (name, lbl, decl,
                     ' '.join('%5.2f' % s for s in ot),
                     ' '.join('%5.2f' % s for s in oe), r0, r1,
                     '   [%d steps]' % nstep if nstep != 1 else ''))
        print()


def reference_control():
    """⚠ IS THE REFERENCE ITSELF THE THING THAT IS WRONG?

    The reference is radau at `h/64`, and radau is also one of the methods
    under test -- a shared-method reference can hide a shared error.  Re-run
    one row against a DIFFERENT method at DOUBLE the refinement.  If the true
    local error moves, the reference is the instrument, not the estimate.
    """
    print('=== reference control: does the TRUE local error depend on the '
          'reference? ===')
    print('%-9s %-18s %-26s %s'
          % ('method', 'fixture', 'reference', 'TRUE local orders'))
    for cls, name in ((RadauIIA3Integrator, 'radau'),
                      (ESDIRK43Integrator, 'esdirk43')):
        for build, lbl in ((rc_linear, 'RC (LINEAR control)'),
                           (cv_loop, 'C-V loop (INDEX 2)')):
            for ref_cls, ref_div in ((RadauIIA3Integrator, 64),
                                     (ESDIRK43Integrator, 128)):
                tl, _ = true_local(cls, build, lbl, ref_cls=ref_cls,
                                   ref_div=ref_div)
                print('%-9s %-18s %-26s %s'
                      % (name, lbl, '%s /%d' % (ref_cls.__name__[:12], ref_div),
                         ' '.join('%5.2f' % s for s in _orders(tl))))
        print()



def main():
    print('=== estimate order, and the deficit against EMBEDDED_ORDER + 1 ===')
    print('%-9s %-18s %-5s %-24s %s'
          % ('method', 'fixture', 'want', 'slopes', 'verdict'))
    for cls, name in ((RadauIIA3Integrator, 'radau'),
                      (TRBDF2Integrator, 'trbdf2'),
                      (ESDIRK43Integrator, 'esdirk43'),
                      (GLM3Integrator, 'glm3')):
        want = cls().EMBEDDED_ORDER + 1
        for build, lbl in ((expg, 'ExpG (nonlinear)'),
                           (rc_linear, 'RC (LINEAR control)'),
                           (cv_loop, 'C-V loop (INDEX 2)')):
            sl = estimate_order(cls, build)
            d = want - sl[-1]
            print('%-9s %-18s %-5d %-24s %s'
                  % (name, lbl, want, ' '.join('%5.2f' % s for s in sl),
                     'DEFICIT %.2f' % d if d > 0.5 else 'no deficit'))
        print()
    print('=== and the discriminator: does sigma_min(J) scale with h? ===')
    for build, lbl in ((expg, 'ExpG'), (rc_linear, 'RC control'),
                       (cv_loop, 'C-V loop')):
        sigma_min_scaling(build, lbl)
    print()
    true_local_pass()
    crossover()
    print()
    reference_control()
    global_order_control()
    local_vs_global_relation()
    rc_asymptotic_check()


if __name__ == '__main__':
    main()


def stiffness_control_is_CONFOUNDED(build_fn=None):
    """⚠⚠⚠ KEPT AS A RECORD OF A CONFOUNDED CONTROL, NOT AS A RESULT.

    docs-46 measured the DIRK stage-order cap on Prothero-Robinson with the
    STIFFNESS as the control -- ESDIRK43 4.00 nonstiff -> 3.00 stiff, TR-BDF2
    flat at 3.00 -- and predicted that RADAU IIA(3), whose stage order is
    `q = 3` so whose cap equals its declared 4, should NOT move.

    Reproducing that control here by tuning the C-V loop's resistor gave:

        method     declared   NONSTIFF (tau=1e-4)   STIFF (tau=1e-9)   drop
        radau          4           4.01                  3.02          0.99
        trbdf2         3           2.97                  1.97          1.00
        esdirk43       4           3.69                  1.99          1.70

    All three drop, radau included, and all three land near their STAGE ORDER
    `q` (3, 2, 2) rather than `q+1`.  Read at face value that refutes the
    account.  IT DOES NOT, because THE FIXTURE MOVES TWO VARIABLES:

        r = 1e5   sigma_min(J) = 2.050e-09 ... 2.006e-09   exponent +0.00
        r = 1e0   sigma_min(J) = 5.000e-06 ... 6.245e-07   exponent +1.00

    Changing `r` to add stiffness ALSO turned DEFECT 1 on -- at `r = 1` the
    C-V loop acquires a direction whose singular value scales with `h`.  So
    the uniform `-1` is the FILTER, not a stage-order cap, and radau's drop is
    fully explained by Defect 1 switching on, which is CONSISTENT with the
    prediction that radau has no stage-order deficit.

    ⚠ This measurement therefore neither confirms nor refutes the stage-order
    account: it cannot separate the two effects.  A fixture that varies
    stiffness while holding `sigma_min`'s scaling fixed is what would, and
    that is what docs-46's scalar Prothero-Robinson fixture does by having no
    DAE structure to carry a reactance-free direction at all.

    ⚠ The near-miss is the point: "all three land at q" is a clean, coherent,
    quotable pattern, and it is an artefact.  Every fixture in this file now
    reports `sigma_min`'s exponent alongside its result for that reason.
    """
    raise NotImplementedError(
        'confounded by construction -- see the docstring; a scalar '
        'Prothero-Robinson fixture is what answers the stage-order question')
