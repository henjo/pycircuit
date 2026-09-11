"""WHERE the filtered local-error estimate loses an order, and where it does not.

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
from pycircuit.circuit.elements import R, C as Cap, VSin
from pycircuit.circuit.analysis import remove_row_col
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
