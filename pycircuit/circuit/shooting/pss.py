"""`PSS`, the periodic steady state by shooting: the class, its parameters and
`solve` in its phases.  The class inherits its themes from the `_pss_*`
modules.
"""
from copy import copy
import numpy as np
import warnings
from pycircuit.circuit.analysis import Analysis
from pycircuit.circuit.analysis import Parameter
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.circuit import gnd
from pycircuit.post import InternalResultDict
import pycircuit.circuit.analysis as analysis
from ._numerics import freq_analysis
from ._numerics import periodic_spline_weights
from ._pss_accuracy import _AccuracyChecks
from ._pss_events import _StateEvents
from ._pss_grids import _PeriodGrids
from ._pss_inner import _InnerTransient
from ._pss_newton import _ShootingNewton
from ._pss_periodic import _PeriodicStates
from ._pss_ppv import _PPVFloquet
from ._pss_replays import _FactoredReplays
from ._pss_walks import _PeriodWalks
from .diagnostics import algebraic_conditioning
from .diagnostics import topological_index


## ⚠ THE SPECTRAL RADIUS CANNOT DECIDE THIS, and trying it first is the
## instructive part.  An autonomous orbit gives an eigenvalue at exactly 1
## -- but only AT its own period, and a run at any other period reads well
## below (measured 0.9615 on the quadrature phase element at the nominal
## period against 1.000226 at the corrected one).  Worse, a merely
## lightly-damped DRIVEN circuit sits near 1 too: a Q=1000 resonator has
## `exp(-pi/Q) = 0.99686`.  So no threshold separates the two -- one
## setting misses the autonomous case where users will actually run it, the
## other fires on every high-Q filter.
##
## The distinction is structural, not spectral, and it is exact: a circuit
## is autonomous when nothing in it depends on `t`.  `u(t)` is sampled
## across the period and compared; a phase accumulator driven by a DC
## source is autonomous however energetically it oscillates, which is
## precisely the case arc 5 asks about.
AUTONOMOUS_U_TOL = 1e-12


class _SolveRun(object):
    """The state one `PSS.solve` carries between its phases (see `solve`):
    set by `_solve_prepare`, extended by `_shoot`, read by the rest."""

    def __init__(self, **fields):
        self.__dict__.update(fields)


class PSS(_ShootingNewton, _PeriodGrids, _StateEvents,
          _InnerTransient, _PeriodWalks, _FactoredReplays, _PPVFloquet, _AccuracyChecks, _PeriodicStates,
          Analysis):
    """Periodic Steady-State using shooting Newton iterations

    The algorithm is described in [1] p65.

     1. Kenneth S. Kundert, Jacob K. White, Alberto Sangiovanni-Vincentelli
        (1990)
        Steady-State Methods for Simulating Analog and Microwave Circuits
        Kluwer Academic Publishers
        ISBN 0792390695

    **THREE CONVERGENCE CHECKS NEST HERE, AND THEY ARE NOT INTERCHANGEABLE.**

    1. the per-timestep Newton, which solves the discretised circuit
       equations at one time point;
    2. the local truncation error, which decides how far the DISCRETE
       trajectory is from the true one;
    3. the shooting Newton, which finds the periodic point of the discrete
       map.

    Two rules order them.  **(3) cannot be tighter than (1)**: the period
    map is only KNOWN to the accuracy of the per-timestep solves, so the
    shooting residual has a floor there whatever its Jacobian is.  And
    **LTE must not run per shooting iteration**: an
    adaptive grid makes the step sequence a function of `x0`, so the period
    map stops being smooth and (3) loses its quadratic rate.  Choose the
    grid once, freeze it, shoot on it.

    **(3) is a true Newton for both methods.**  It was neither, for years:
    the monodromy accumulation was missing a factor of `1/h`, and since `C`
    is singular the product collapsed to exactly zero -- so the Jacobian was
    the identity and this was successive substitution, `x0 <- phi(x0)`, on
    every circuit and without saying so.  It was found by its RATE: the
    residual fell 0.855 per iteration on a Q=20 resonator, and `exp(-pi/Q)
    = 0.8546` is that circuit's own per-period decay.  Non-convergence is
    reported now; it used to be discarded with `full_output=False`.

    Trapezoidal needed more than the factor, and Gear-2 more again.  Every
    method here writes its companion as

        iq_n = sum_k a_k q_{n-k}  +  b iq_{n-1}

    so ONE recursion differentiates all of them:

        S    = sum_{k>=1} a_k C_{n-k} Px_{n-k}  +  b Pq
        Px_n = -Jf_n^-1 S
        Pq_n = a_0 C_n Px_n + S

    Euler is `b = 0` reaching back one step, trapezoidal `b = -1` reaching
    back one, Gear-2 `b = 0` reaching back two.  The coefficients come from
    `Integrator.companion_coefficients` -- from the integrator that ACTUALLY
    ran, so an order-dropped opening step contributes its own -- rather than
    being transcribed here, which this tree has paid for three times.

    An x-only monodromy is not merely less accurate for the methods with
    memory: measured, applying the Euler form to trapezoidal converged
    SLOWER than no Jacobian at all (0.90 against 0.855).  A wrong Jacobian
    is worse than none.

    ⚠ **THE TWO FAILURE MODES ARE ORTHOGONAL, and level 3 cannot see level
    2.**  On the Q=20 resonator against a 20 V analytic peak, all three now
    converge -- and they do not agree:

        euler   5 iterations,  8.815 V   (56% low)
        gear2   7 iterations, 19.766 V   (1.2% low)
        trap    6 iterations, 19.990 V   (0.05% low)

    That is each method's own numerical damping, invisible to levels 1 and
    3, and it is why the LTE report in the recorded scope below is worth
    more than it looks: a converged shooting solve is not by itself
    evidence of a correct answer.

    (2) DOES NOT EXIST HERE, AND UNDER A SHOOTING METHOD IT CANNOT BE A
    CONTROLLER.  The grid is a fixed `linspace`; under `fixed_timestep` both
    transient backends skip the LTE verdict outright, keeping only the
    order drop that protects the integrator across a breakpoint.  Nor is
    that a wiring gap to be closed: if the step sequence adapts to `x0`
    then phi is a DIFFERENT discrete map for each `x0`, so it is not smooth
    in `x0` and the accumulated monodromy is the derivative of a
    neighbouring problem.  Freezing the grid is what makes (3) a Newton.

    So (2) changes kind here.  The estimator is still computable on a
    frozen grid, and what it measures is not convergence but ACCURACY:
    which of the three levels is limiting the answer.  That question is
    live rather than academic -- on a Q=20 resonator with `method='euler'`
    the shooting solve converges completely (5 iterations, residual
    3.9e-05) and lands at 8.815 V against a 20 V analytic peak.  The answer
    is 56% low for a reason (1) and (3) cannot see, and nothing currently
    says so.

      4. LTE AS A REPORT -- DONE (2026-09-01).  `max_lte`, `total_lte` and
         `max_lte_seam`, measured on the final replay through
         `Transient.step_lte`; see the block in `solve` for what each one
         means and why one number was not enough.  Two things it taught,
         neither of them anticipated by the paragraph above:

           - THE PER-STEP PEAK PASSES THE 56%-LOW ANSWER.  At reltol=1e-3
             euler's peak LTE on that resonator is 0.288, in tolerance,
             because the estimator bounds ONE step and the 56% is what 99
             of them do together.  A transient is right to control on that
             number; a periodic analysis cannot report only it.  The sum
             over the period reads 26.27, against gear2's 0.941 and trap's
             0.340 -- tracking amplitude errors of 55.9%, 1.17% and 0.05%.
           - THE COLD-START SEAM IS PART OF THE MAP.  `_begin_period`
             re-seeds a flat history every shooting iteration -- which is
             exactly what keeps phi a function of `x0` alone -- so the
             discrete period map opens off a past that never happened, and
             that defect is inside the map the solve converged on.

         ⚠ AND A THIRD THING, WHICH CORRECTED THE SECOND (same day, by
         measurement -- `benchmarks/pss_seam_cost.py`).  The report first
         flagged a seam for ALL THREE methods and called it the dominant
         term for the multistep pair.  Measured against the limit cycle the
         same grid and method reach with a real history, the seam costs
         **5.1e-12 V for euler and 1.3e-11 V for trapezoidal -- zero** --
         while the report was calling them 0.286 and 15.1 times tolerance.
         Only Gear-2 pays: 1.266e-01 V at 100 points/period, 54% of its
         total error, rising to 73% at 400 points because the seam falls as
         h^2 while the interior falls faster.

         The discriminator is HOW FAR THE COMPANION REACHES, not how far the
         estimator does.  Euler reads `q_{n-1}`; trapezoidal reads `q_{n-1}`
         and `iq_{n-1}`, which the order-dropped opening step supplies
         consistently.  Gear-2 reads `q_{n-2}`, which at that step is the
         entering unknown -- and shooting constrains `x(0) = x(P)`, NOT
         `x_in` to be the orbit's `x(-dt)`, so it is an O(h^2) stand-in
         being read as a history point.  Removing the order drop makes
         Gear-2 WORSE (2.34e-1 -> 4.39e-1 at 100 points): the drop is
         protective and 1.266e-01 is the residue it leaves.

         ⚠ AND ON EVERY SHIPPED METHOD IT IS `None`, so the treatment below
         describes a number the code as configured cannot produce.  Verified
         2026-09-02 across all three: the seam is collected only for a
         method whose companion reaches two charges back, which is Gear-2
         alone -- and Gear-2 always takes the solved-history path, which has
         no manufactured opening and clears the flag.  Euler and trapezoidal
         reach one.  So the figures below are reachable only through the
         superseded plain-Gear route (a test gets at them by monkeypatching
         it), and the `_limits` warning entry that quotes them is dead code
         on the shipped paths.  Kept because the reasoning is what the
         figure is FOR, and because a fourth method reaching two charges
         back on the plain path would revive both.

         So `max_lte_seam` is a FLAG, not a magnitude -- at 100 points the
         estimator's seam/interior ratio is 505x and the answer's is 1.18x.

    ⚠ THREE THINGS IN THIS CLASS ENLARGE SOMETHING, AND ALL THREE WERE
    ONCE CALLED "AUGMENTED".  One word, three referents, in one file --
    which is how a reader ends up applying a statement about one of them to
    another.  They are now named apart, and the names are worth learning:

      the `(x, iq)` MONODROMY   what trapezoidal's period map differentiates.
                                Its recursion carries a companion current, so
                                an x-only monodromy is structurally
                                incomplete.  About the DERIVATIVE, not the
                                unknowns.
      the FREE-PERIOD system    what an autonomous circuit solves: unknowns
                                `(x0, T)` with a phase condition, because the
                                period is not given.  `func_autonomous` in
                                `solve`, on every kind's `_pmap`.
      a SOLVED ENTERING HISTORY unknowns `(x0, x_{-1})`, 4b.  About where
                                the period map STARTS, not how long it runs.
                                DRIVEN circuits.  `self.solved_history`.
      the COMPOSED system       unknowns `(x0, x_{-1}, T)`, 4c -- the second
                                and third TOGETHER, for an autonomous
                                circuit under a two-step method.  ⚠ NOT a
                                synonym for either half: quoting 4b's
                                numbers as "the composed system's" is an
                                error this docstring has already made.

    ⚠ AND "THE SEAM IS REMOVED" IS NOT "THE ANSWER IS RIGHT".  Below,
    `exact` means the SEAM is gone -- the solve lands on the limit cycle the
    same grid and method reach from a real history.  It does NOT mean the
    analytic answer: 4b's Gear-2 at 100 points/period returns 19.89297 V
    against 20 V, still 1.070e-01 V out, and that residue is ordinary
    interior discretisation error which no history fix touches.  Read
    `exact` as "seam-free", never as "error-free".

      4b. A SOLVED ENTERING HISTORY FOR A TWO-STEP COMPANION (DRIVEN
         circuits; the autonomous composition is 4c) -- DONE (2026-09-01),
         and it is the remedy that measurement pointed at.  `(x_0, x_{-1})`
         are unknowns together and both must close; see
         `_traverse_solved_history`.  Applied where the COMPANION reaches two
         charges back (`_companion_reach`), which is Gear-2 alone -- euler
         and trapezoidal keep the plain path because their seam measured
         zero, and enlarging their system would double the unknowns to fix
         nothing.

         THE GATE WAS THE PREDICTION, NOT AN IMPROVEMENT.  If the seam is
         the only difference between PSS's answer and the cycle a real
         history produces, removing it must LAND on that cycle.  It does,
         to 2.5e-07 V:

              points   plain      primed     solved-hist error     gain
                 100   19.76639   19.89297   19.89297    2.34e-1 -> 1.07e-1  2.18x
                 200   19.95451   19.98524   19.98524    4.55e-2 -> 1.48e-2  3.08x
                 400   19.99008   19.99735   19.99735    9.92e-3 -> 2.65e-3  3.74x

         ⚠ Read the `error` column: 19.89297 is NOT 20 V.  The seam is gone;
         1.070e-01 V of interior discretisation error remains, untouched.

         ⚠ AND IT IS CHEAPER, WHICH WAS NOT THE EXPECTATION.  Two residual
         evaluations against twelve, 4.3x faster wall-clock on that circuit.
         The plain path seeds BOTH sensitivity rings with `I` -- which is
         the flat-history assumption written into the Jacobian -- so its
         Newton was inexact and nobody could see it, because Newton
         converges anyway from an approximate Jacobian.  The
         solved-history one is exact.

      4c. THE COMPOSED SYSTEM -- `(x_0, x_{-1}, T)`, AUTONOMOUS ONLY -- DONE
         (2026-09-01), the pair's free-period residual (`func_autonomous` on
         `_pmap` in `solve`).  It is 4b's enlargement AND the free-period one at once, because an autonomous
         circuit under a two-step method needs both: the period because it
         is not given, the history because the companion reads it.  4b alone
         does not cover such a circuit and 4b's numbers are not this one's.

         ⚠ REFUSED FIRST, ON THE MEASUREMENT BELOW, THEN BUILT BECAUSE IT
         WAS ASKED FOR.  The evidence did not change; the decision did, and
         it was the owner's to make.  The measurement still stands and is
         still the reason `trap` is the default -- see the value caveat at
         the end of this item.

         WHAT THE SEAM DOES TO AN OSCILLATOR IS NOT WHAT IT DOES TO A DRIVEN
         CIRCUIT, and the guess was backwards.  The expectation was that it
         would matter MORE, the period being an unknown a per-period kick
         could land in.  Measured on the quadrature phase element it is the
         opposite: the seam moved the solved period by 2.5 ppm of a 332 ppm
         error (0.75%, and 0.38% at 400 steps) against 54% for the driven
         resonator.  The reason is one item up -- a shooting fixed point
         ABSORBS a once-per-period perturbation, and a free period is one
         more degree of freedom to absorb it into -- so the kick landed in
         the orbit's SHAPE, as a radius wobble of 2.095e-04 on an orbit of
         radius 1, instead of its frequency.

         Composing fixes both, and the wobble is the visible half:

              points   plain       free-running   composed
                 200   +329.682    +332.184       +332.185   ppm
                 400    +82.342     +82.652        +82.651   ppm
              radius wobble at 200: 2.095e-04  ->  6.6e-12

         The free-phase eigenvalue survives the enlargement -- the composed
         run reads `spectral_radius` 1.000000 -- so the autonomous
         diagnostic still says what it said.  It now reads the FULL 2m x 2m
         map, whose spectrum carries the parasitic roots of the two-step
         discretisation alongside the physical multipliers.

         ⚠ THAT WAS RECORDED AS AN OPEN WORRY AND HAS BEEN MEASURED: the
         parasitic roots are not a problem, and the reason is quantitative.
         Gear-2's parasitic root is 1/3 per STEP -- the roots of
         `1.5z^2 - 2z + 0.5` are 1 and 1/3 -- so over a period it is
         `(1/3)^N`, about 1e-95 at 200 points.  The autonomous 16x16
         spectrum measures as [1.000, 3.5e-06, 5.2e-16, 4.7e-17, 0, ...]:
         one physical unit eigenvalue and nothing else above rounding, so
         `max |eig|` picks the physical one.  A method whose parasitic root
         sat nearer the unit circle would need them separated; Gear-2's does
         not, and a test pins the gap so a future method cannot inherit the
         assumption silently.

         ⚠ THAT CLEAN SEPARATION IS A PROPERTY OF THE CIRCUITS TESTED, NOT
         OF THE METHOD, and it is known to degrade in one specific place:
         HIGH-Q OSCILLATORS, whose PHYSICAL multipliers cluster near 1.
         FOUR INDEPENDENT WITNESSES, and the earliest is Demir (IJCTA 2000)
         on a Colpitts oscillator: "the four largest eigenvalues of the
         monodromy matrix ... all four eigenvalues are +1.  THIS IS USUALLY
         THE CASE FOR HIGH-Q OSCILLATORS.  In fact, several eigenvalues can
         become very close to 1 such that they are NOT NUMERICALLY
         DISTINGUISHABLE from the one that is supposed to be equal to 1
         theoretically."  Then Bizzarri et al. -- "[shooting] is not suited
         to simulate oscillators based on very high quality resonators since
         these lead to fundamental matrices with eigenvalues very close to
         1" -- and Demir & Roychowdhury on the PPV -- "the oscillatory-mode
         eigenvalue of 1 ... cannot be distinguished from other eigenvalues
         of the system that are close to 1.  This is particularly true for
         many LC oscillators." -- and this codebase's own measurement below.

         ⚠ AND THE HISTORY SAYS WHAT TO DO ABOUT IT, which is why it is
         recorded rather than just cited.  Demir's 2000 remedy was to SELECT
         the eigenvector with the largest inner product against
         `C(0) xdot(0)` -- measured 0.2 against 1e-5, 1e-7, 2e-5.  His own
         2003 paper rejects that heuristic: "no guarantee that any of the
         candidate eigenvectors will be appreciably more orthonormal than
         the others, leading to a potential breakdown."  The SAME VECTOR
         then changes role: sampled `C(t) u_1(t)` becomes the augmented row
         `q`, so no selection happens at all.  The quantity used to CHOOSE
         among candidates becomes the CONSTRAINT that makes the candidate
         unique.  That is the whole of the 2003 improvement, and it is why
         a PPV built here must go to the augmented solve and NOT to the
         eigenvectors this method returns.

         ⚠ ONE CAUSE, AND IT SURFACES IN THREE PLACES HERE, which is why it
         is written down once instead of three times as coincidences:
           - `max |eig|` above, and `_spectral_report`'s split, both assume
             the unit root is identifiable.
           - THE PHASE ROW below removes the singularity from the unit
             eigenvalue and only that one.  Bizzarri et al. again: "this is
             not enough if any other eigenvalue is close to 1 and in this
             case possibly ill conditioned matrices must be managed."  The
             bordered system stays formally nonsingular and gets badly
             conditioned, so the failure is a slow, loud Newton rather than
             a `LinAlgError` -- the shape that invites the wrong diagnosis.
           - PPV eigen-selection, if it is ever built, picks the same root
             and inherits the same limit.
         Measured here only that the unbordered composed null space is
         EXACTLY 1-D on the circuits tested (sigma_min/sigma_next = 1.2e-11)
         -- true, and NOT a claim about the high-Q case, which has not been
         measured.

         ⚠ WHAT THE SAME MEASUREMENT DID FIND, in the DRIVEN path rather
         than the autonomous one: `_traverse_solved_history` was handing
         back `Px[0][:, :m]`, the `d x_{N-1}/d x_0` CORNER of the
         sensitivity, as the monodromy.  A corner of a sensitivity is not a
         monodromy and its eigenvalues mean nothing -- it reported
         `spectral_radius` 1.279605 for the Q=20 resonator, ABOVE ONE and so
         reading as an unstable orbit, where the analytic per-period decay
         is exp(-pi/Q) = 0.854636 and every other path reports 0.855.  Fixed
         to the pair map, which now reads 0.854833 -- closer to analytic
         than the plain path's 0.853369 -- and pinned against the analytic
         decay rather than against itself.

         ⚠ THE VALUE CAVEAT, UNCHANGED BY BUILDING IT.  Gear-2's own phase
         error is +332 ppm against trapezoidal's +83 ppm at the same grid,
         both second order, so on THIS circuit composed Gear-2 is still the
         worse choice and `method='trap'` -- the default -- remains the
         right answer.

         ⚠⚠ AND THE STIFF-CIRCUIT JUSTIFICATION DOES NOT REPRODUCE.  This
         item used to say the composition earns its place because a stiff
         autonomous circuit needs a two-step method, citing
         `doc/transient_review.md` sec. 4.6 -- trapezoidal ringing at
         `|e_n/e_{n-1}| ~ 0.9960` at `h*lambda = -1000` where Gear-2 sits at
         0.0972.  ⚠ THOSE ARE RINGDOWN NUMBERS.  They measure a TRANSIENT,
         and a periodic steady state has no transient to ring; the citation
         was carried across contexts without checking that it transfers.

         Measured on two stiff autonomous circuits: the phase element plus a
         fast RC at exactly `h*lambda = -1000`, and a diode peak detector
         whose orbit has a fast edge every period.  Trapezoidal shows NO
         ringing in either -- the alternating signature is identical between
         the methods and falls at ~h^3 under refinement, so it is the sharp
         edge being resolved, not an undamped mode -- and trapezoidal is 4x
         BETTER on frequency at both grids (+83.084 against +332.180 ppm at
         200 points; +1.287 against +5.126 at 1600).

         So no circuit measured so far prefers Gear-2 for autonomous PSS.
         The composition still earns its place on its OWN evidence and does
         not need that story: without it an autonomous Gear-2 run is
         silently biased in the period by 2.5 ppm and its orbit does not
         close in radius, so `method='gear'` was quietly WRONG there rather
         than merely inferior.  Making an offered method correct is the
         justification; "and it is the better method for stiff oscillators"
         was mine and is unsupported.

         ⚠ THE UNTRIED CASE IS NOW TRIED AND THE CLAIM IS REFUTED, not
         merely unsupported.  Van der Pol at `mu = 100` -- the canonical
         stiff relaxation oscillator, fast mode IN the orbit, measured
         stiffness ratio 5443 (edge timescale 0.0299 against a period of
         162.842412):

              method   npts        h   outcome         period      err ppm
              trap     2000   0.0815   NoConvergence        -            -
              gear     2000   0.0815   NoConvergence        -            -
              trap     8000   0.0204   NOT converged   162.813755        -
              gear     8000   0.0204   NoConvergence        -            -
              trap    20000   0.0081   converged       162.832543    -60.6
              gear    20000   0.0081   converged       162.823215   -117.9

         Trapezoidal wins on both counts: it is the only method that
         produced a finite answer at 8000 points, and at 20000 it is TWICE
         as accurate.  Across three circuits no case has been found where
         Gear-2 is the better choice for autonomous PSS, so the default
         stands and this half of the justification is closed as refuted.

         ⚠ AND THE BINDING CONSTRAINT TURNED OUT TO BE THE GRID, NOT THE
         METHOD.  Neither method runs that circuit below 20000 points,
         because this analysis freezes a UNIFORM grid (which is what makes
         (3) a Newton) and the edge needs `h < 0.01` against a period of
         162.8.  The adaptive transient that produced the reference used
         ~1160 points per period, so the uniform grid costs about 17x the
         points on this circuit class.  That is a measured argument for
         RECORDED SCOPE ITEM 5, the LTE-chosen grid -- and the first one it
         has had.

    4d. THE CHEAP APPROXIMATE ALTERNATIVE, measured and NOT shipped.  The
        4b system is seam-free but doubles the unknowns, so it is worth
        knowing what an approximation buys.  (Method H was measured against
        4b, the DRIVEN system -- not against 4c.)  `q_{-1}` can be BUILT from
        `x_0` instead of solved for, keeping the system at m unknowns.

        ⚠ AND THE PLAIN PATH IS ALREADY THE FIRST-ORDER MEMBER OF THAT
        FAMILY, which is the fact that reframes the whole question.  Its
        entering charge is EXACTLY `q_0 - h qdot_0` -- checked against the
        converged iterate to 1.6e-38 relative -- because backward Euler on
        the opening step says exactly that.  So "add a pseudo-history" is
        not an alternative to what shipped before 4b; it IS what shipped
        before 4b, at first order, and it is what measured 1.266e-01 V.

        THE STRUCTURAL FACT UNDERNEATH: a converged step satisfies
        `i(x) + iq + u = 0`, so `qdot = -(i + u)` is EXACTLY available with
        no solve, while `qddot` needs `xdot = C^-1(...)` and C is singular
        in MNA.  The DAE gives away the first derivative of the charge and
        refuses the second.  ⚠ That also explains 4b's measurement rather
        than merely restating it: TRAPEZOIDAL needs only `iq_{-1}`, which is
        that free derivative, so it is exactly initialised and has no seam
        (1.3e-11 V).  Gear-2 needs a second CHARGE, which no residual
        equation supplies.  That asymmetry is the whole story.

        Carrying it one term further, using derivatives rather than fitting
        charges (better conditioned: a quadratic fit `q_{-1} = 3q_0 - 3q_1 +
        q_2` has coefficients summing to 7 in magnitude and amplifies
        inner-solve noise):

            q_{-1} = q_0 - (3h/2) qdot_0 + (h/2) qdot_1,  error (5/12) h^3

        `x_1` is one throwaway backward-Euler predictor; its O(h^2) error
        enters with coefficient h/2 and lands at O(h^3), so a first-order
        predictor suffices.  Measured on the Q=20 resonator
        (`benchmarks/pss_seam_cost.py`, `solve_back_extrapolated`):

              points   seam plain   seam H     share plain   share H
                 100   1.266e-01    1.985e-03      1.183      0.0185
                 200   3.074e-02    1.175e-04      2.083      0.0080
                 400   7.272e-03    8.769e-06      2.742      0.0033

        64x smaller at 100 points, and falling 16.9x then 13.4x per halving
        against the plain path's 4.1x -- at least the h^3 predicted, better
        than that on these grids.  The share is what matters: the plain
        seam GROWS as a fraction of the error and method H's VANISHES.

        ⚠ AND IT DOES NOT SUBSTITUTE FOR 4c.  Measured on the autonomous
        phase element at 200 steps/period, where the frequency is the
        unknown and so the thing to watch:

              formulation        its   period ppm   err vs ref   wobble
              plain (x_in, T)     18     329.682      2.50      2.095e-04
              method H (x0, T)     3     329.720      2.46      5.8e-06
              composed 4c         14     332.185      0.001     6.6e-12

        H removes the orbit's WOBBLE -- 5.8e-06 against 5.0e-04 for a flat
        seed, so the construction itself works -- and leaves 98% of the
        FREQUENCY error.  It fixes the shape of the orbit and not its
        period.

        ⚠ THAT SPLIT CORRECTS THIS DOCSTRING'S OWN EARLIER ATTRIBUTION.  4c
        above called the 2.5 ppm "what the seam does to the period".  It is
        not: the two are separate defects that one subtraction
        (`|plain - free-running|`) had lumped together, because that
        difference contains both and can be attributed to whichever one the
        reader has in mind.

              the SEAM (entering-history error)  ->  the radius WOBBLE.
                  H fixes it, 86x.
              the MISSING CLOSURE (one equation
                  where a two-step map needs two) ->  the FREQUENCY, 2.5 ppm.
                  H does not touch it; only 4c does.

        It also sharpens the literature note: "k conditions for a k-step
        method" is NOT a statement about initialising history.  History can
        be initialised perfectly -- H does -- and the answer is still wrong
        if k-1 conditions are missing.  The two are independent
        requirements.

        THE REASON IS A MISSING EQUATION, NOT AN INACCURATE HISTORY, and
        method H is what PROVES it rather than merely suggesting it, because
        it changes exactly one variable:

              formulation   |q_-1 - q_N-2|/|q|   order   freq err   order
              plain              4.984e-04        h^2     2.502      h^3
              method H           3.005e-05        h^3     2.464      h^3
              4c (solved)        0                 --     0.001       --

        H cuts the history mismatch 16.6x and lifts its order from h^2 to
        h^3 -- and the frequency moves 1.5%.  A 16.6x better history buys
        nothing, so the history is not what is wrong.

        What is wrong is the equation count.  For a two-step method the
        discrete state is the PAIR, so periodicity is a condition on the
        pair: 4c imposes BOTH closures (`x_0 = x_{N-1}` and `x_{-1} =
        x_{N-2}`), while H imposes one and CONSTRUCTS the other state.  The
        period that closes the first component is not the period that closes
        the pair, and that gap is O(h^3) however well the history is built.  This is the literature note's
        "k conditions for a k-step method" arriving as a number: H supplies
        one (plus the phase row), 4c supplies two.  On a DRIVEN circuit the
        period is given, so the missing condition has nowhere to go and H
        does fine; on an autonomous one it goes straight into the period.

        ⚠ AND H CONVERGES FASTEST, TO THE WRONG ANSWER -- three Newton
        iterations against 4c's fourteen and the plain path's eighteen, with
        `ier == 1` and the residual satisfied.  This class already records
        that a converged shooting solve is not evidence of a correct answer
        (see the successive-substitution defect above); here it is again, in
        a formulation built the same day.  Fast convergence on a smaller
        system is not a merit when the system is missing an equation.

        SO WHY IS IT NOT SHIPPED.  ⚠ NOT ON SPEED -- that argument was made
        here and was WRONG.  Method H's 15 residual evaluations are an
        artefact of the flat Jacobian it was handed for the ACCURACY
        measurement; given a good one (finite differences) it converges in
        3, against 4b's 2.  That is not a difference worth a decision, and
        comparing an exact-Jacobian formulation against a deliberately
        crippled one was not a fair comparison.

        What is actually left against it: it would be a SECOND formulation
        for a job 4b already does seam-free, in a tree that has paid more
        than once for duplicate paths -- and its exact Jacobian is
        unwritten.  `d q_{-1}/d x_0 = C_0 + (3h/2) G_0 - (h/2) G_1
        dx_1/dx_0` needs deriving, including the predictor's own
        sensitivity; until it exists H is either slow (flat Jacobian, 15) or
        expensive (finite differences, m+1 traversals per iteration).

        For the DRIVEN case that leaves a maintenance judgement -- a second
        formulation for a job 4b already does seam-free -- which belongs to
        whoever owns the trade.  For the AUTONOMOUS case it is no longer a
        judgement at all: H is measurably wrong there, by 98% of the
        frequency error, and 4c is not optional.

        WHAT WOULD REOPEN IT: a circuit large enough that the 2m x 2m dense
        `J_phi` factorisation dominates -- Kundert puts that above a few
        hundred unknowns.  ⚠ But note where that argument really points:
        with the matrix-free Krylov solve of recorded scope item 6, the
        enlargement costs 2x (vector length), not 8x (factorisation).  So
        the scaling case is an argument for item 6 first, and only then for
        approximating the history.

    WHAT THE LITERATURE SAYS ABOUT 4b (checked 2026-09-01, because the fix
    above looked like something that ought to be standard):

      THE GENERAL RESULT IS CLASSICAL, and it is exactly what 4b hit.  A
      k-step linear multistep method turns a first-order continuous problem
      into a k-th ORDER DISCRETE one, which introduces parasitic (spurious)
      solutions and needs **k conditions** to determine the discrete
      solution -- not one.  That is why a single periodicity condition on
      `x0` is under-determined for Gear-2.  It is the founding observation
      of BOUNDARY VALUE METHODS (Brugnano & Trigiante), which supply the k
      conditions as "one initial and k-1 final", chosen at both ends
      deliberately because it improves stability.  A periodic BVP hands
      them over for free: requiring the whole k-tuple to close is what
      gear's pair residual (`_pmap` in `solve`) does.  The parasitic roots
      are also the spurious
      eigenvalues that appear in the 2m x 2m monodromy, which is why the
      autonomous eigenvalue-at-1 diagnostic would need redefining if the two
      systems were ever composed.  Standard reference for the BVP side:
      Ascher, Mattheij & Russell, "Numerical Solution of Boundary Value
      Problems for ODEs" (SIAM, 1995).

      THE CIRCUIT LITERATURE ASSUMES THE PROBLEM AWAY, consistently and
      reasonably.  Kundert ("Simulation Methods for RF Integrated Circuits",
      ICCAD 1997) writes the shooting map as `phi_T(v0, 0)` with the state
      `v` alone and gives the sensitivity's component pieces as
      `Jf(v(ts)) = G(v(ts)) + C(v(ts))/hs` -- a ONE-STEP, backward-Euler
      shaped Jacobian, no `q_{n-2}` term, no history in the state.  Gourary,
      Rusakov, Ulyanov & Zharov (MES 2019,
      doi:10.31114/2078-7707-2019-1-25-30) likewise solve "with respect to
      the state vector at the beginning of one period".  Both are correct
      for the one-step methods they are written for; the plain path IS that
      formulation.  The gap opens only when a two-step companion is handed
      to it.  ⚠ Note the 2019 paper's title promises more than it delivers
      here: its case for single-step is A-STABILITY ("the common drawback of
      BDF methods is the lack of A-stability for order higher than 2"), not
      history at the period boundary.

      THE CLOSEST ANYONE COMES is Wambacq, Vandersteen, Phillips,
      Roychowdhury, Eberle, Yang, Long & Demir, "CAD for RF circuits", which
      argues for one-step Chebyshev-IRK discretisation and says of it:
      "Each step is independent of the ones before and after".  That
      independence is precisely the property whose ABSENCE is the seam --
      but they argue it from stability and step adaptivity, not from
      initialisation.

      ⚠ AND THE STANDARD CASE AGAINST GEAR-2 DOES NOT APPLY HERE, which is
      the part worth keeping.  The same paper's objections are all about an
      ADAPTIVE grid: BDF "not actually as numerically stable as is popularly
      believed"; "the second order Gear method is not A-stable for
      nonuniform steps, and in fact it is not stable for any timestep if the
      ratios between consecutive steps exceed about 2.4"; "a rapid change of
      timestep in a multistep code also necessarily comes with a loss of
      order".  PSS freezes a UNIFORM grid for the whole solve -- that is
      what makes (3) a Newton, see above -- so none of those bite inside a
      run.  Shooting is the one place Gear-2 is on its best behaviour, which
      is an argument for having repaired it rather than refused it.

      WHAT IS NOT IN ANY OF IT: a number.  No source found quantifies what a
      mis-initialised history costs on a circuit.  The figures in 4b are
      this tree's own; `benchmarks/pss_seam_cost.py` is the measurement.

    4e. WHAT THE SCALAR `Idtmod` FORMS DO -- idtmod arc 5, investigated
        2026-09-01, and it turned up a defect that is NOT about `Idtmod`.

        The scalar element runs: with `ic` set, DC solves and a transient
        integrates it correctly, the state advancing one modulus per output
        period and the gauge shift (`_apply_periodic_shifts`) wrapping it
        back.  ⚠ That wrap is why shooting does not simply refuse it -- the
        state DOES return to itself once wrapped, so `x0 - phi(x0)` is
        satisfiable, and `IdtmodQuadrature`'s docstring claim that the
        scalar form "structurally cannot" close is true of the RAW state
        and not of the wrapped one this tree actually integrates.

        But asking whether the answer was right exposed the general defect:

        ⚠ AN AUTONOMOUS PERIOD IS DETERMINED ONLY UP TO AN INTEGER MULTIPLE,
        AND THE SOLVE FOLLOWS ITS SEED.  `k*T` satisfies the periodicity
        condition whenever `T` does.  Measured on the quadrature element,
        true period 1.000e-03:

              seed      solved         converged?
              1e-3      1.000083e-03   yes
              2e-3      2.000665e-03   yes
              3e-3      3.002245e-03   yes

        Every one is a correct periodic solution and every one reports
        success; the FUNDAMENTAL FREQUENCY -- usually the thing a PSS user
        wanted -- is wrong by the factor.  This is a property of the
        free-period system (4c), not of `Idtmod`, and it was silent.
        Detected now: see the recurrence check in `solve`, which needs no
        extra solve and sets `fundamental_period`.  Driven runs are exempt,
        their period being the caller's.

        ⚠ AND THE SECOND DEFECT IS ALSO NOT ABOUT `Idtmod`.  What looked
        like the scalar form's seed-fragility -- Gear-2 collapsing to
        `T ~ 1e-18`, trapezoidal raising a bare `LinAlgError` -- reproduces
        on the QUADRATURE element too: from a 1e-4 seed against a 1e-3
        fundamental, Gear-2 returns -1.5e-20 there and trapezoidal dies from
        three seeds of five.  `T = 0` is a REGULAR ROOT of every autonomous
        shooting system, because `x0 - phi_T(x0)` vanishes identically
        there and the phase condition constrains `x0`, not the period.  So
        any seed below the fundamental is drawn to it.

        Neither outcome is silent -- the collapse reports
        `converged = False` (⚠ since 2026-09-06: it previously reported
        `converged = True`, because `T = 0` is a REGULAR root that `fsolve`
        reaches cleanly and reports success on; `_free_period_solve` now
        demotes `ier` when it detects the collapse) and the exception is
        loud -- but neither named
        its cause, and the generic non-convergence advice ("raise
        maxiterations") is actively wrong for it: no number of iterations
        reaches a fundamental from below.  `_free_period_solve` now names
        both, and says what to do instead (seed at or above the expected
        period; a short transient and the interval between two output
        recurrences gives one).

        ⚠ THREE TIMES IN THIS ITEM the defect looked like the element and
        was the formulation.  The scalar form's raw state not closing is
        answered by the gauge shift; the multiple-period ambiguity and the
        trivial root are properties of the free-period system that any
        autonomous circuit has.  Arc 5's deliverable turned out to be two
        diagnostics on 4c, not a change to `Idtmod` at all.

        STILL OPEN, and genuinely about the scalar form: seeded correctly
        it converges, but a pure phase accumulator has no amplitude to pin,
        so its orbit is a ramp and nothing distinguishes one starting phase
        from another beyond the phase condition itself.  Whether that is a
        limitation worth removing has not been established.

    RECORDED SCOPE, in order, neither of these planned work yet:

      5. LTE-CHOSEN grid -- THE MECHANISM IS BUILT (2026-09-01), the
         payoff case is not yet reached.  `PSS.solve(grid=...)` takes step
         FRACTIONS of the period and freezes them; see `_period_grid`.

         ⚠ THE RECORDED BLOCKER WAS STALE.  This item said "blocked on
         `Transient` accepting a non-uniform grid; `fixed_timestep` is
         uniform-only".  That loop IS uniform-only -- and PSS never uses
         it, driving `solve_timestep` one step at a time instead, where
         non-uniform steps worked unchanged.  Verified before anything was
         written.

         Verified on benign circuits: driven and autonomous, both methods,
         2:1 and smoothly-varying grids all converge, with the driven
         spectral radius still matching exp(-pi/Q).  The autonomous case is
         what makes FRACTIONS the contract rather than times -- the grid is
         rebuilt at the current `T` on every residual evaluation, so
         `dh/dT = h/T` keeps holding.

         VAN DER POL, THE CIRCUIT THAT MOTIVATED IT, NOW SOLVES THROUGH IT
         (2026-09-02).  On its own 1105-step LTE-chosen grid, one added
         opening step: trapezoidal converges at -73.8 ppm and Gear-2 at
         -100.6, against the 20000 UNIFORM points a uniform grid needs for
         -60.6.  18x fewer points, and Gear-2 solves it for the first time
         on any grid.  Pinned by
         `test_the_lte_chosen_grid_solves_van_der_pol_through_the_analysis`.

         ⚠ THE BLOCKER WAS THE MANUFACTURED OPENING STEP, AND THIS
         DOCSTRING PREVIOUSLY NAMED IT WRONG.  It read: "what is left is
         that the throwaway used FINITE DIFFERENCES ... a stiff relaxation
         oscillator is the first circuit that cannot tolerate the plain
         path's 30%".  MEASURED: an exact finite-difference Jacobian on the
         real analysis does NOT fix it -- it fails identically -- and once
         the opening step is subdivided the analytic and finite-difference
         Jacobians agree to SIX DIGITS (-73.8 ppm both).  The ~30% was
         never what stopped this case.

         What stopped it: `_traverse` manufactures `x(0)` with ONE
         order-dropped Euler step of `hs[0]`, and a grid taken from an
         adaptive transient opens wherever that transient's window happened
         to start.  Here `h[0] = 1.4845` against a MEDIAN step of 4.62e-04
         -- 3200x coarser -- and the INNER Newton fails on that one step
         (100 iterations, residual 1e7x over tolerance at `v`).  The
         throwaway never met this because its unknown was `x_0` itself, so
         it had no manufacturing step to take.  `_period_grid` now opens on
         the grid's own finest step; see the guard there.

         ⚠ THE OLD DIAGNOSIS HAD ALREADY SEEN THE CAUSE AND FILED IT AS AN
         EXONERATION -- "the opening step is large not small" was written
         down as one of the obvious causes ELIMINATED.  It was checked
         against the wrong worry (that the opening step might be too small
         to open the trajectory) and never against its own size.

         ⚠ THE FIX PROPOSED HERE WAS TRIED AND IS WRONG AS STATED.  It
         read: `iq_{-1} = -(i(x_0) + u(t_0))` is exactly available from the
         DAE (item 4d), so trapezoidal can take a solved-history
         formulation seeded that way -- an exact Jacobian, clear of the
         fold.  The seeding half is true and was built.

         The formulation half is not.  A one-step companion depends on
         `x_{-1}` ONLY through `iq_{-1}`, and `d(iq_{-1})/d x_{-1} = -G` is
         SINGULAR at every purely reactive node -- most of a resonator.  So
         admitting `x_{-1}` as m unknowns leaves the 2m x 2m system
         rank-deficient, and it fails exactly as it should:
         `LinAlgError: Singular matrix`, on 25 tests at once.

         ⚠ THE CORRECTED DESIGN WAS ALSO WRONG, and was also built and
         measured before that was known.  It read: for a `b != 0` method
         the second unknown is `iq_0` itself -- the `(x, iq)` state its
         monodromy already uses -- closed by `iq_0 = iq_{N-1}`.

         It works, and then it does not.  On the Q=20 resonator it cut the
         outer solve to TWO residual evaluations (against euler's ten) and
         returned a peak identical to the plain path's 19.98968, confirming
         both that the prize is real and that trapezoidal's fixed point was
         never in doubt.  But `spectral_radius` came back 1.000000 where
         the circuit decays by exp(-pi/Q) = 0.854636 -- and that is the
         tell:

              npts   steps           outcome
               200     199 (odd)     converged
               201     200 (even)    LinAlgError
               202     201 (odd)     converged
               203     202 (even)    LinAlgError

         Trapezoidal's `iq` recursion is `iq_n = ... - iq_{n-1}`, whose
         homogeneous mode is `(-1)^n` -- UNDAMPED, exactly on the unit
         circle.  Over an even number of steps it returns `+1`, so `I - M`
         is singular and the solve dies; over an odd number it returns `-1`
         and nothing shows.  A user choosing 201 points instead of 200 gets
         a crash.

         ⚠ BOTH FAILED ATTEMPTS SHARE ONE CAUSE: trapezoidal's `iq` is not
         a coordinate a periodicity condition can pin.  Solving for a
         previous STATE is rank-deficient (`d(iq)/dx = -G` is singular at
         reactive nodes); solving for `iq` itself is degenerate (its mode is
         marginally stable, so closing it is vacuous or singular).

         SO THE PLAIN PATH IS CORRECT FOR TRAPEZOIDAL, not a legacy wart.
         Its Euler manufacturing step SUPPLIES `iq_0` from the DAE
         deterministically -- `-(i(x_0) + u(t_0))`, item 4d -- instead of
         asking a marginal mode to close.  That is why trapezoidal has no
         seam (1.3e-11 V) and converges where the enlargements do not.  The
         inexact Jacobian (item 4b, ~30%) is the price, and well-posedness
         is what it buys.

         A THIRD DESIGN WAS DERIVED AND ALSO FAILS, IDENTICALLY.  Make
         `iq_0` DEPENDENT rather than free -- `iq_0 := -(i(x_0) + u(t_0))`,
         unknown `x_0` alone, seeds `P_x(0) = I` and `P_q(0) = -G(x_0)`.
         It looked strictly better: m x m, no manufacturing step, exact
         Jacobian, and it needs only the FORWARD derivative `-G`, never the
         inverse that killed the first attempt.

         Measured, it dies on the same line:

              npts   steps           outcome        evals      rho
               200     199 (odd)     converged          2   1.000000
               201     200 (even)    LinAlgError        1        --
               202     201 (odd)     converged          2   1.000000
               203     202 (even)    LinAlgError        1        --

         Making `iq_0` dependent does not avoid the alternating mode -- it
         EXCITES it.  A perturbation `dx_0` gives `diq_0 = -G dx_0`, which
         drives the `(-1)^n` mode, which returns undamped with multiplier
         `(-1)^N`.  So the unit eigenvalue lands in `d x_end/d x_0` itself:
         even the m x m monodromy reads 1.000000, against the circuit's
         true 0.854636.

         ⚠ THREE DESIGNS, ONE OBSTRUCTION, AND IT NAMES WHAT THE PLAIN PATH
         IS FOR.  Every reformulation couples `iq_0` to `x_0` through `-G`,
         and that coupling excites a mode that never decays.  The plain
         path couples them through the EULER COMPANION instead -- its seed
         is `a_0 C`, not `-G` -- and that coupling is not degenerate: it
         reports 0.855 where all three of these report 1.000000.  The
         manufacturing step is not scaffolding to be removed; it is what
         makes trapezoidal's shooting problem well-posed, and the ~30%
         Jacobian error is what that costs.

         ⚠ THE "~30% JACOBIAN ERROR" IS WRONG IN BOTH DIRECTIONS, and the
         real statement is stronger.  Measured 2026-09-02 against a
         finite-difference `dF/dx_in` on three circuits at 100 points:

               circuit   method   relerr(J_code, J_true)   rank(J_true)/m
               RC        euler/trap        1.56                 1/3
               Q=20 RLC  euler             6.6e-03              2/4
               Q=20 RLC  trap              3.7e-03              2/4
               nonlinear euler/trap        4.73                 1/3

         So it is 0.4% on one circuit and 470% on another -- "~30%" was a
         single-circuit number carried as if it were a property of the
         method.

         ⚠ AND THE FRAMING WAS WRONG TOO: `J_code` IS NOT AN APPROXIMATION
         TO `dF/dx_in`.  That derivative is SINGULAR on every circuit above
         -- `sigma_min` is exactly 0.0 and the rank is 1/3, 2/4, 1/3 -- so
         AN EXACT NEWTON FOR THIS FORMULATION DOES NOT EXIST, and a solver
         handed the true Jacobian fails on the first step rather than
         converging faster.  What `_traverse` returns is `I - dx_end/dx_0`
         with a flat-history seed: an approximation to a DIFFERENT and
         WELL-POSED derivative, taken with respect to `x_0` while the
         unknown is `x_in`.  That is why the plain path works at all, and
         why its convergence is LINEAR rather than quadratic -- measured on
         the autonomous phase circuit, trapezoidal's residual falls
         3.91e-03 -> 3.14e-04 -> 2.66e-05 at a constant ratio ~0.076 where
         Gear-2's solved-history route is quadratic.  It is a contraction,
         not a Newton.

         ⚠ AND THE LITERATURE SETTLES THE DESIGN QUESTION: `x_0` IS THE
        CANONICAL UNKNOWN, so this is a RETURN and not an invention.
        Aprille & Trick (Proc IEEE 60(1) 108-114, 1972), "The Problem":
        "determine the periodic state w(0) such that integrating (1) from
        the initial state w(0) over the interval [0, T] we obtain the
        periodic solution", with Step 1 "for the given initial state x_0^i
        compute the solution x^i(t; x_0^i), 0 <= t <= T".  The trajectory
        begins AT the unknown; there is no pre-image and nothing
        manufactured before the first step.  Same in the oscillator paper
        (IEEE TCT 19(4) 354-360) and in Trick, Colon & Fan (TCAS 22(5)
        391-396) eq.(2).

        ⚠ AND THE IDENTITY SEED IS THEIRS TOO, which is what makes the frame
        error precise.  [AT-P] gives the discrete sensitivity for backward
        Euler as `Phi(T,0;x_0) = PROD [I - h F(x)]^-1` -- k factors seeded by
        `z(0)` at `t = 0`, i.e. the identity AT `x_0`, which in their
        formulation IS the unknown.  This file kept A&T's JACOBIAN and
        changed A&T's UNKNOWN.  `Px = [I, I]` is not wrong in isolation; it
        is correct for the formulation it came from and wrong for the one it
        now sits in.

        ⚠ WHAT PRECEDENT DOES NOT SETTLE is the COST, because A&T never paid
        it: they used BACKWARD EULER, which is L-stable, so the `(-1)^n`
        obstruction that forces an L-stable opener here never arose for
        them.  Their Phi is derived at fixed step for one method; the k-step
        and non-uniform-grid questions are later work.

        ⚠ WHICH MAKES THE FIX NAMEABLE: make `x_0` the unknown.  The
         throwaway driver that solved van der Pol did exactly that (its
         unknown was `x_0` itself), and item 5's note already half-records
         it.  Not built here; written down so the next attempt starts from
         the right statement rather than from "the Jacobian is 30% off".

         NOT PROPOSING A FOURTH.  Bordering the `(x, iq)` system to remove
         the alternating mode remains formally available -- the analogy to
         the phase condition is exact -- but the null direction there is a
         property of the discretisation with no closed form to pin, unlike
         `xdot(0)`, and three derivations in this item have now looked sound
         and failed on contact.  Anyone taking it up should measure before
         building: the falsifier is cheap and is even/odd step counts.

         ⚠ AND THE CONCLUSION IS NOW A THEOREM, NOT A TALLY OF THREE FAILED
         ATTEMPTS.  This item blamed trapezoidal's `iq` RECURSION for the
         `(-1)^n` mode.  That is the wrong cause, and an external review
         (2026-09-02) found the right one: trapezoidal is A-stable but NOT
         L-STABLE, so it maps the null space of the singular MNA `C` by
         exactly `-1` per step.  For `C x' + G x + u = 0` the one-step
         amplification is

             A_trap  = (C/h + G/2)^-1 (C/h - G/2)   ->  -I  on null(C)
             A_euler = (C/h + G)^-1 (C/h)           ->   0  on null(C)

         so the count of `-1` modes is exactly `m - rank(C)` -- one per
         ALGEBRAIC variable -- on every MNA circuit, and it has nothing to
         do with `iq`.  Verified here to 1e-09: the Q=20 resonator has
         `m = 4`, `rank(C) = 2` and exactly two eigenvalues at `-1` under
         trapezoidal and two at `0` under Euler; a plain RC has `m = 3`,
         `rank(C) = 1` and two of each.  A review reproduced the same
         singularity in an X-ONLY formulation containing no `iq` variable at
         all (cond 6.9e+18 at even K, 8.2e+04 at odd), which is the
         falsifier for the old attribution.

         What that buys, beyond a correct cause: the conclusion generalises.
         ANY formulation whose period map is `A_trap^K` without an L-stable
         opening step is singular at even `K` on every MNA circuit -- so
         "the plain path is correct for trapezoidal" is provable rather than
         observed.  And the cure is not specific to Euler: any L-STABLE
         opening step works, and it rescues exactly `m - rank(C)` modes.
         Corroborated in the literature the review checked: Houben (2003,
         App. A) biases theta off 1/2 so that "the numerical oscillations
         due to the DAE character of the equations are damped out during
         the 'insensitive time'" (p. 22, verbatim; an earlier rendering here
         elided three words without an ellipsis, caught by the docs
         session's quotation audit 2026-09-08) -- the same mechanism,
         attributed there to the DAE.

         Kept from the attempt: `_install_history` now takes the entering
         step size `h_prev` separately.  `x_{-1}` sits one step BEFORE
         `x_0`, and on a periodic grid that step is `hs[-1]`, not `hs[0]`.
         Uniform grids hide this; a caller's grid with a 16438:1 spread
         does not.

      5b. LTE-CHOSEN grid, the original wording -- pick the step sequence from an adaptive run and
         freeze it, refining BETWEEN shooting solves.  The grid still never
         moves inside one, so (3) stays exact.  Blocked on `Transient`
         accepting a non-uniform grid; `fixed_timestep` is uniform-only.
         Note this is the better structure anyway: a transient adapts
         because it cannot see the future, while PSS re-solves the same
         interval repeatedly and can therefore choose its grid once, well.
      6. MATRIX-FREE variational shooting (Telichevesky, Kundert & White,
         DAC 1995) -- the only structure that permits per-iteration
         adaptivity, because `M v` is obtained by integrating the
         variational system along the stored trajectory on that iteration's
         own grid.  The outer solve becomes an INEXACT Newton, phi shifting
         slightly between iterations.  A rewrite, not an increment on the
         above, and it should not leak into one.

         ITS CASE GOT STRONGER THREE TIMES ON 2026-09-01, from work that was
         not about it:

           - 4b/4c DOUBLE THE UNKNOWNS for a two-step method, so whatever
             the dense `J_phi` costs, this analysis now pays it on a vector
             twice as long.
           - THE LITERATURE NOTE puts a number on when that bites: Kundert
             (ICCAD'97) has forming and factoring `J_phi` at O(N^2 S) and
             O(N^3), "intractable when N exceeds several hundred", and names
             matrix-implicit Krylov as the answer -- this item.
           - 4d's SCALING ARGUMENT for method H pointed here instead.  H's
             only advantage over 4b was keeping the system at m unknowns for
             a dense solve; with a matrix-free solve the enlargement costs
             2x (vector length) rather than 8x (factorisation), so item 6
             removes H's reason to exist rather than competing with it.

         ⚠ AND THE COST TO ATTACK IS NOT THE FACTORISATION.  The final
         `J_phi` factorisation is one O(m^3).  The SENSITIVITY PROPAGATION
         that builds it is `N` steps of `_step_sensitivity`, and each of
         those is O(m^3) TWICE OVER: `linearsolver(Jf, S)` with a 2m-column
         right-hand side, and the `C @ P` products that form `S` and `Pq`
         -- m x m against m x 2m.  Matrix-free replaces the whole step, not
         half of it: it never forms `P`, so both go, leaving one matvec per
         Krylov iteration.

         MEASURED (2026-09-02, `benchmarks/pss_matrix_free_ceiling.py`), on
         a quiet box, single-threaded BLAS, every reading reproducible to
         better than 0.5%:

               m   propagation   (solve alone)   ceiling k=20
              40          4.8%            2.3%          1.02x
             110         15.0%            6.2%          1.14x
             242         38.1%           14.1%          1.54x
             502         63.8%           21.3%          2.58x
            1002         79.1%           24.9%          4.44x

         THE GATE IS PASSED.  It was "does the propagation share pass
         ~30%": it crosses near m=220 and reaches 79% at m=1002.  Item 6 is
         JUSTIFIED above m~250 and stays POINTLESS at m=40.

         AND IT IS NOW BUILT (2026-09-02) for the DRIVEN solved-history
         path: `solve(matrix_free=True)`, on `_traverse_factored` +
         `_monodromy_matvec` + the matrix-free Newton.  The recursion is not
         duplicated -- `_step_sensitivity` took a `solve` argument and is
         still the one recursion, run at width 1 against a stored
         factorisation instead of width 2m.  End to end: 1.36x at m=242,
         1.51x at m=502, 2.13x at m=1002, and the answers agree with the
         dense path to 1.1e-16.

         ⚠ k WAS ALSO A GUESS AND IS NOW MEASURED: GMRES takes 2/4/7/12
         iterations at m=40/110/242/502, because `I - M` has almost every
         eigenvalue within 1% of 1.0 -- the fast modes decay to nothing
         over a period, so k tracks the number of SLOW MODES, not m.  That
         is the property the item actually rests on, and it was assumed
         rather than checked until now.

         ⚠ WHAT IT COSTS: `2 N m^2` doubles of stored factorisations
         (~800 MB at m=1002, 50 points) where the dense path holds O(m^2);
         one MORE Newton iteration at m>=502 (3 against 2) from a
         convergence test that cannot be identical; and no monodromy, so
         `spectral_radius` is None on this path.  The autonomous and plain
         systems are NOT converted and raise rather than silently going
         dense.

         ⚠ THIS OVERTURNS THIS DOCSTRING'S OWN 2026-09-01 RECORD, and the
         box being quiet is NOT why.  That record read "the propagation is
         2.2% of a traversal at m=40, ceiling 1.01x-1.03x", and it is what
         happens when a measurement is named for the thing it was meant to
         settle rather than the thing it timed: it accumulated
         `toolkit.linearsolver` alone, which is under a third of the step
         at every size measured.  The 2.2% was reproducible, stable, and
         answering a different question.  ⚠ It also produced the confident
         wrong sentence "most of the time is ASSEMBLY -- which matrix-free
         does not touch"; at m=502 the propagation is 64% and assembly is
         not what dominates.  The m=40 VERDICT survives on the corrected
         number (4.8%, 1.02x) -- being right about m=40 is what made the
         error cheap to keep.

         ⚠ QUOTE THE THREADING CONDITION WITH THE NUMBER.  With BLAS
         threads free the same sizes read 4.8 / 7.2 / 12.4 / 18.2%: the
         `C @ P` products thread and the Python-level assembly does not.
         Single-threaded is the trustworthy column here -- the threaded
         traversal moved 2.64 s -> 4.02 s at m=502 between runs on this
         box while every single-threaded reading held to 0.1%.  A machine
         giving BLAS all its cores sees a smaller prize.

         ⚠ AND `k=20` IS AN ASSUMPTION, NOT A MEASUREMENT.  The ceiling
         charges matrix-free `k/m` of the propagation and nothing else, so
         it is an upper bound on an upper bound.  What the Krylov solve
         actually costs on these systems is unmeasured.

    Driving `Transient` -- done -- buys one integrator definition, the
    limiting/PCNR machinery, breakpoints and the order drop.  It does NOT
    buy (2) as a CONTROLLER; saying otherwise was this docstring's own
    earlier error, and item 4 above is what (2) turned out to be able to be
    instead.

    `reltol`, `iabstol` and `vabstol` mean exactly what they mean on
    `Transient` -- the tolerances of the TRANSIENT solution, applied to the
    per-timestep Newton, with the two absolute floors applied PER UNKNOWN in
    both flavours by `analysis.newton_tolerance_vectors`, the single
    definition all three analyses read.  Nothing here rescales them.

    The shooting criterion is expressed against that one: `steadyratio`
    (>= 1, default 1) multiplies it, so by default the shooting solve is
    held to the SAME relative tolerance as the transient, and raising it
    buys fewer shooting iterations for a looser periodic steady state.

    CHOOSING `method`: WHAT EACH ALTERNATIVE GIVES UP (owner decision,
    2026-09-08: the default is chosen for ACCURACY, not cost -- "we do not
    want to fool the user; instead they should change integrator and know
    the impact").  Every entry below is MEASURED in this tree
    (doc/pss_roadmap_260902.md A10; the radau-default section and its
    monotonicity fixtures moved to doc/pss_log_260902.md in the 2026-09-11
    plan/log split, and the plan's index names it); nothing is quoted from a
    textbook order alone.

      method    order   period error   error        above its monotonicity limit
                        (ppm @ 400pt,  ESTIMABLE    (h_FE = 1/steepest slope)
                        Q=1e4 vdP)     by refinement
      radau     5 (6.1  5.7e-10        yes          OUTPUT stays in range at every
                on a                                step measured; only the STAGES
                smooth                              leave it, <= 0.3 % of the swing,
                orbit)                              above ~4 h_FE -- seen ONLY on a
                                                    scalar square-wave fixture (the
                                                    widest practical margin of the
                                                    four); on a driven detector
                                                    circuit the stages never leave
                                                    the hull at 20-320 pts/period,
                                                    with a diode OR a tanh device
      esdirk43  4       5.3e-05        yes          output in range; stages leave it
                                                    <= 3 % above ~2.5 h_FE
      trbdf2    2       10.1           yes          OUTPUT rings above 2.4 h_FE
      trap      2       20.8           yes -- its   OUTPUT rings above 2 h_FE
                                       PPV surfaces
                                       are its
                                       TR-BDF2
                                       twin's (the
                                       recorded
                                       "sign change"
                                       was a defect,
                                       fixed
                                       2026-09-14)
      gear      2       83.1           yes          not an RC; a small ring at
                                                    h >= 10 tau measured; does not
                                                    CERTIFY a free-period solve at
                                                    1e-14 below ~400 pts at Q=1e4.
                                                    FIRST-CLASS ON NON-UNIFORM
                                                    GRIDS (2026-09-20, all
                                                    measured): the PPV, the
                                                    diffusion constant, the period
                                                    (free or driven), the adjoint
                                                    modes and every noise fold are
                                                    second order on any grid whose
                                                    step ratios stay inside the
                                                    zero-stability bound 1+sqrt(2)
                                                    (a 2:1 alternating grid costs
                                                    NOTHING over uniform: +4.86 vs
                                                    +4.87 ppm); beyond it (3:1) the
                                                    period is FIRST order and
                                                    `_period_grid` warns.  Events:
                                                    `break_events` is on, second
                                                    order across a landed edge; the
                                                    one cost is a CONSTANT at the
                                                    step after a hard corner (30x
                                                    trap's there).  Index 2: order
                                                    2 in BOTH subspaces (no BDF
                                                    order reduction, measured on the
                                                    C-V loop 2.01/2.01 uniform,
                                                    2.00/2.00 smooth); the PPV and
                                                    c on a non-uniform index-2 grid
                                                    second order since 2026-09-20
                                                    (the index-2 fallback's samples
                                                    come from the continuous adjoint
                                                    there, and PAC's period weights
                                                    are a periodic trapezoid).  ⚠
                                                    Gear is the ONLY method whose
                                                    oscillator surfaces are on a
                                                    non-uniform grid at all: the
                                                    one-step kinds' factored
                                                    periods and the TR-BDF2 twin
                                                    replay on a uniform grid
                                                    whatever the solve used.
                                                    Adjoints: the exact transpose
                                                    for the PPV and the noise folds,
                                                    a second-order continuous
                                                    adjoint for the mode shapes on
                                                    non-uniform grids (dense to 8
                                                    unknowns, Ritz-certified
                                                    Arnoldi above).  Cost: one real
                                                    factorisation per step, 3.0-11.8x
                                                    cheaper than radau as n grows
                                                    (table below); TR-BDF2 matched
                                                    it 1.5-1.7x better at equal
                                                    accuracy on E2's fixture, so
                                                    price both with `grid_error`.
                                                    ADAPTIVE GRIDS (`lte_grid`,
                                                    2026-09-21, relaxation van
                                                    der Pol mu = 10, each method
                                                    on the grid its own run made,
                                                    seeded from `lte_period`):
                                                    gear 195 pts -1461 ppm (second
                                                    order: -361 / -92 under 2x /
                                                    4x splitting; a window that
                                                    read -22 ppm sat on a SIGN
                                                    CHANGE of the error, -22 ->
                                                    +7 -> +3 -- refine before
                                                    believing a good number),
                                                    trbdf2 286 pts -8.7 ppm,
                                                    radau 110 pts -0.32 ppm.
                                                    Windows differ for every
                                                    second-order method (trbdf2
                                                    on the same two gear grids
                                                    +10 / -180 ppm; radau -0.01 /
                                                    -0.04): the adaptive run's
                                                    growth-by-2 steps land
                                                    differently each period.
                                                    ⚠ THE FOLD'S "GEAR SWING"
                                                    WAS THE FOLD'S ALIGNMENT
                                                    (2026-09-21): a rescale
                                                    after overshooting the
                                                    period put the orbit's
                                                    edges on the boundary of
                                                    the fine groups -- six
                                                    folds split 5x for gear
                                                    (+1409..+1563 vs +273..+302
                                                    ppm) AND trbdf2 on the same
                                                    grids (+150..+173 vs +44..
                                                    +46).  Walk fixed: gear +16
                                                    ..+30 (-44 at 192 pts),
                                                    trbdf2 +4..+13.  Gear is a
                                                    steady ~8x trbdf2 on any
                                                    given grid, both second
                                                    order; choose by that
                                                    factor, not by a swing.
                                                    For NOISE on a coarse or
                                                    folded grid choose trbdf2
                                                    or radau (exact per-stage
                                                    injection); gear's
                                                    covariance is first order
                                                    in h/tau.  STATE EVENTS
                                                    (a comparator, a threshold
                                                    switch; `state_events`,
                                                    2026-09-22): radau's
                                                    staged solve 9-21x closer
                                                    at 50-200 points on a PWM
                                                    loop; trbdf2's own
                                                    second-order error on the
                                                    switch-off decay dominates
                                                    there (staged no better
                                                    than unstaged) -- use radau
                                                    on such circuits; gear
                                                    skips the stage.  The adaptive
                                                    run's reltol must be tighter
                                                    than the transient's habit:
                                                    1e-5 is too coarse above
                                                    mu ~ 10 for every method
                                                    (the mu = 100 benchmark used
                                                    1e-7).  Pass the period
                                                    `lte_grid` observed
                                                    (`pss.lte_period`) to
                                                    `solve`: the grid is
                                                    fractions of it.
      euler     1       --             yes          never rings; damps the orbit it
                                                    is asked to find (13 % of the
                                                    amplitude at 20 pts/period)

    Cost: radau is one real plus one complex factorisation per step (a
    coupled 3n system) against one per step for the others; at the point
    counts above it is cheaper on wall-clock anyway (60 points beat trap's
    480 on both axes, roadmap radau-default section), and `grid_error` /
    `warping_estimate` price the trade on YOUR circuit.  AS `n` GROWS
    (measured 2026-09-08, ladder oscillator, 200 points, same grid for
    every method): radau/gear on one transient period 3.0 / 3.5 / 3.9 /
    5.7 at n = 12 / 32 / 102 / 302, and on the whole PSS solve 3.0 / 3.9 /
    4.2 / 4.2 -- per-step assembly dominates to n ~ 100 and the 3n stage
    system only starts to show at 300; radau/trbdf2 1.5 -> 2.6; esdirk43
    is the MOST expensive method at every size (four sequential stage
    solves), never the cheap alternative.  Absolute: n = 302 at 200 points
    is 20 s per period and 271 s per PSS solve under radau, 64 s under
    gear.  At n = 1002 (transient period only, same grid): gear 26 s,
    trbdf2 55 s, radau 311 s -- radau/gear 11.8x, radau/trbdf2 5.7x -- so
    the 3n stage factorisation dominates from a few hundred unknowns and
    the ratio roughly doubles per 3.3x in n there; above ~300 unknowns
    the accuracy is bought at an order of magnitude in wall-clock, and
    trbdf2 (order 2, contractive, R = 1 + sqrt 2) is the alternative to
    price against `grid_error` on your circuit.  Index-2: the
    period keeps classical order under radau (6.1 measured on a smooth
    orbit); the algebraic unknowns converge at the stage order (3) -- and a
    relaxation oscillator whose PERIOD is timed by such a variable (a
    comparator on a current sense through a C-V loop) STILL keeps order 5
    (5.66 / 5.06 / 5.05 measured, identical errors to its voltage-sensed
    twin from 400 pts on), so no caveat applies.  ⚠ What that fixture DID
    show: its inner step Newton fails undamped at every grid (stage
    sensitivity tau/h, basin 0.94 h/(k tau)) -- hence the line search as
    the last resort -- and `warping_estimate` reads 0.09-0.65 of the true
    error on it: the septic interpolant does not resolve a comparator
    edge, and the estimate is only as good as the interpolant (Part I's
    "moderately smooth" boundary).  Roadmap, radau-default section.
    """

    parameters = Analysis.parameters + \
        [Parameter(name='analysis', desc='Analysis name',
                   ## Sources supply their time-domain waveform only for an
                   ## analysis name in timedomain_analyses (('dc','tran')); the
                   ## old default 'PSS' matched nothing, so cir.u(t) returned 0
                   ## and the whole shooting solve had no excitation.
                   default='tran'),
         Parameter(name='reltol', 
                   desc='Relative tolerance', unit='', 
                   default=1e-4),
         Parameter(name='iabstol', 
                   desc='Absolute current error tolerance', unit='A', 
                   default=1e-12),
         ## 1e-6 since 2026-09-19 -- DC, Transient, JAXTransient and PSS share one
         ## meaning and one default; the reason is at `Transient.vabstol`
         Parameter(name='vabstol', 
                   desc='Absolute voltage error tolerance', unit='V', 
                   default=1e-6),
         Parameter(name='maxiter',
                   desc='Maximum number of iterations', unit='',
                   default=100),
         ## Forwarded to the inner Transient so PCNR (the junction-continuation
         ## limiting) reaches the shooting per-step solve too -- it lives in
         ## `Transient.solve_timestep`, which PSS DOES call, so no per-accepted-
         ## step machinery is needed (unlike breakpoints / continuation rescue,
         ## which are armed in `Transient.solve` and stay out of reach).
         Parameter(name='pcnr',
                   desc='Use Predictor/Corrector Newton-Raphson instead of '
                        'limiting in the inner transient; off by default',
                   unit='', default=False),
         ## 2026-09-21 (B7c completed): which step lengths move with an
         ## unknown period.  'proportional' rescales every step with T;
         ## 'closing' keeps a caller's inner steps at their absolute lengths
         ## and lets the last step close the period.  'auto' is 'closing' on
         ## a caller's grid for an autonomous run and 'proportional' otherwise
         ## -- see `_period_grid` and the note at `_period_column`.
         Parameter(name='period_column',
                   desc="'auto' (= 'closing' + proportional polish on a caller's grid, 'proportional' on a uniform one), 'proportional' or 'closing': "
                        "which step lengths depend on an unknown period; see "
                        "the note at the policy in solve()",
                   unit='', default='auto'),
         ## ⚠⚠ THE ONE KNOB `method='theta'` HAS, AND IT WAS UNREACHABLE.
         ## `_integrator_for` builds `table[method]()`, so every shooting run
         ## took `ThetaIntegrator.DEFAULT_C` -- a RATE, calibrated on ONE
         ## fixture's period.  The transferable quantity is the DIMENSIONLESS
         ## `C T` (see `ThetaIntegrator.DEFAULT_CT` for why `h` cancels), so
         ## that is what this parameter is, and `_theta_biased` turns it into
         ## the rate THIS period needs.  `None` takes the measured knee.
         Parameter(name='theta_ct',
                   desc="method='theta' only: the null(C) damping one PERIOD "
                        'applies, as the dimensionless product C*T. None '
                        'takes ThetaIntegrator.DEFAULT_CT (the measured '
                        'knee). Ignored by every other method.',
                   unit='', default=None),
         ## `reltol` MEANS THE SAME THING IN EVERY ANALYSIS: the relative
         ## tolerance of the transient solution.  It is applied to the
         ## per-timestep Newton here exactly as `Transient` applies it, and
         ## nothing rescales it.
         ##
         ## `steadyratio` is how the SHOOTING criterion is expressed relative
         ## to it: shooting reltol = reltol * steadyratio, with 1 meaning the
         ## two are equal.  It is >= 1 because the period map is only KNOWN
         ## to the accuracy of the inner solves, so asking the shooting
         ## residual to beat that is asking it to resolve its own noise --
         ## refused rather than silently accepted.  Raise it to accept a
         ## looser periodic steady state for fewer shooting iterations.
         ## The LTE floors, separate from the Newton ones for the reason
         ## `Transient` records: one knob must not move both criteria.  Same
         ## names, same defaults, same meaning -- this analysis reports the
         ## number a transient would have controlled on.
         Parameter(name='lte_vabstol',
                   desc='Absolute voltage floor for the truncation-error '
                        'estimate', unit='V', default=1e-12),
         Parameter(name='lte_iabstol',
                   desc='Absolute current floor for the truncation-error '
                        'estimate', unit='A', default=1e-12),
         Parameter(name='TRTOL',
                   desc='Truncation error over-estimation factor (SPICE '
                        'TRTOL / lteratio in a commercial simulator)', unit='', default=7.0),
         Parameter(name='relref',
                   desc="What the relative LTE tolerance is measured "
                        "against: 'pointlocal', 'alllocal' or 'sigglobal'",
                   unit='', default='sigglobal'),
         Parameter(name='steadyratio',
                   desc='Shooting tolerance as a multiple of reltol (>= 1); '
                        '1 holds the shooting solve to the same relative '
                        'tolerance as the transient, larger relaxes it',
                   unit='', default=1.0),
         ## ⚠⚠ THE DEFAULT IS `radau` (owner decision, 2026-09-07), CHANGED
         ## FROM `trap`.  The floor of this stack is DISCRETISATION and it
         ## grows LINEARLY IN Q; at 240 points per period the relative error
         ## in the diffusion constant against the analytic high-Q reference is
         ##
         ##     Q      gear        trap        radau
         ##      100   1.79e-03    2.63e-05    6.97e-10
         ##      500   9.02e-03    1.32e-04    3.48e-09
         ##     1000   1.82e-02    2.63e-04    6.97e-09
         ##
         ## (240-against-960 grid difference.)  `trap`'s column is its TR-BDF2
         ## twin's `c` (the PPV comes from the twin), `~2.6e-07 Q` at order 3
         ## -- linear in Q like the others.
         ## ⚠⚠ The column recorded here until 2026-09-14 (1.49e-06 / 1.04e-04
         ## / 2.36e-04) was a DEFECT, not trapezoidal: `diffusion_constant`
         ## divided the twin's integral by trap's own period, adding an O(h^2)
         ## term of opposite sign.  That sum was the "error that CHANGES SIGN
         ## near Q = 100" which `grid_error` refused and which was cited for
         ## this default; fixed, trap is estimable.  The default still stands
         ## on accuracy, and more clearly than before (below).
         ## `radau` is order 5, self-starting (no manufactured opener, so no
         ## seam in the period map), L-stable, and carries its own monodromy,
         ## so an autonomous run takes NO TR-BDF2 twin and reads its own
         ## spectrum -- see `monodromy_twin` and `carries_own_monodromy`.
         ##
         ## ⚠ It IS more expensive per step at a fine grid (3-stage fully
         ## implicit, through the 1-real/1-complex transform): PURE SOLVE time
         ## on van der Pol at Q=100, 480 points, is 3.538 s against trap's
         ## 2.586 s.  At a coarse grid it is cheaper (1.088 vs 1.519 at 120),
         ## the shooting Newton needing fewer iterations without an opener
         ## seam in the period map.
         ##
         ## ⚠⚠ BUT FOR ANY OSCILLATOR SURFACE THE TWIN DOMINATES, AND THAT IS
         ## WHAT SETTLES THE COST QUESTION.  `trap` is not self-sufficient: an
         ## autonomous run must solve a SECOND, TR-BDF2 PSS for its monodromy
         ## (`monodromy_twin`), and `radau` carries its own.  Measured
         ## `diffusion_constant` cost, which pays for that twin:
         ##
         ##     npts   trap solve / c-eval   radau solve / c-eval
         ##      120     1.519 / 1.255         1.088 / 0.187   (no twin)
         ##      480     2.586 / 4.161         3.538 / 0.737   (no twin)
         ##
         ## ⚠ AND AT EQUAL ACCURACY IT IS NOT CLOSE.  Relative error in `c`
         ## against the analytic high-Q reference, with total wall-clock:
         ##
         ##     Q~100   trap  480 pts  3.320e-06   6.747 s
         ##             radau  60 pts  7.599e-07   0.908 s
         ##     Q~500   trap  480 pts  1.661e-05  12.434 s
         ##             radau  60 pts  3.802e-06   0.633 s
         ##
         ## Radau at SIXTY points beats trap at four hundred and eighty, on
         ## both axes at once.  (trap's errors re-measured 2026-09-14 after the
         ## period-normalisation fix; the old 4.061e-06 / 9.230e-06 were the
         ## defect's, and its "non-monotone" 2.913e-06 -> 4.061e-06 too.)
         ##
         ## `trap` remains one argument away for a cheap coarse answer.
         Parameter(name='method',
                   desc="Integration method for the inner transient: 'radau' "
                        "(default, order 5), 'esdirk43', 'trbdf2', 'theta', "
                        "'gear' (BDF-2), 'trap' or 'euler'. The default is "
                        "chosen for ACCURACY; the class docstring's 'CHOOSING "
                        "method' table states, from measurement, what each "
                        "alternative gives up in order, error estimability, "
                        "monotonicity and cost",
                   unit='',
                   default="radau")]        

    
    def __init__(self, cir, toolkit=None, irefnode=None, **kvargs):
        self.parameters = super(PSS, self).parameters + self.parameters
        super(PSS, self).__init__(cir, **kvargs)
        ## The reference row is fixed for the analysis, and both the shooting
        ## loop and the Transient this drives need it.  It was recomputed in
        ## every method from a `refnode` argument that no caller ever varied.
        self.irefnode = self.cir.get_node_index(
            gnd if irefnode is None else irefnode)
        self._tran = None
        ## Only assembled when the period is an unknown; an extra assembly
        ## per timestep is not worth paying on the fixed-period path.
        self._want_dfdh = False
        self._dfdT = None
        ## ⚠ WHICH PERIOD-COLUMN CONVENTION `want_dT` USES.  'proportional'
        ## is the shipped one: every step scales with `T`, `dh_i/dT = h_i/T`.
        ## 'closing' is the commercial one Andreas described -- the inner
        ## transient owns the steps and the LAST one is placed on the period
        ## boundary, so `dh_i/dT = 0` inside and `dh_N/dT = 1` at the close.
        ## MEASURED (roadmap B7c gates 1 and 4): the proportional column is
        ## `O(h)` wrong on a smooth uniform grid and 4.2% RELATIVE wrong on
        ## van der Pol at `mu = 100` (step ratio 16438x), where 'closing' is
        ## 46x closer.  Default unchanged pending the rest of B7c.
        self._period_column = 'proportional'
        self._dfdh = None
        self._want_lte = False
        ## the state-event fractions the last solve landed on (None: none),
        ## and `dtheta/dx_0` from the bordered solve (K x m)
        self._state_event_fracs = None
        self._event_sensitivity = None
        self._event_columns = None
        self._captured = {}
        ## The caller's step fractions, or None for the uniform grid.  Read
        ## by the autonomous closures, which rebuild the grid at the current
        ## `T` on every residual evaluation.
        self._grid_fracs = None
        ## the period `lte_grid` observed in its adaptive run (see there)
        self.lte_period = None
        self._lte = None
        self._lte_seam = False
        self._lte_valid = True
        self._history_is_solved = False
        ## Set by `solve`: whether the entering history joined the unknowns.
        self.solved_history = False
        ## Reported by `solve`: the peak normalised truncation error over the
        ## converged period, and where in the period it fell.  None until a
        ## solve has run, or when the grid was too short to difference.
        self.max_lte = None
        self.max_lte_time = None
        self.max_lte_seam = None
        self.total_lte = None
        ## Set by `solve` when an autonomous run lands on a multiple of the
        ## fundamental; None when it did not (or on a driven run).
        self.fundamental_period = None
        ## What `factored_period()` needs to replay the CONVERGED period:
        ## which seed, which grid, which opening.  Written at the end of
        ## `solve`, and cleared at its start so a failed or interrupted run
        ## cannot leave a previous solution's state readable as this one's.
        self._period_state = None
        self._factored_period_cache = None
        ## ⚠ WHICH MONODROMY THE OSCILLATOR SURFACES READ -- see
        ## `monodromy_twin`.  Selects the method that supplies the PPV,
        ## Floquet modes and factored period when a one-step LMM (trap/euler)
        ## solved an autonomous circuit, whose OWN monodromy is first-order
        ## on a limit cycle (the opener seam).  'trbdf2' (DEFAULT): a TR-BDF2
        ## twin on the same grid -- self-starting, no opener, measured 12-32x
        ## more accurate on lambda2 than the Gear-2 twin at practical step
        ## counts (see `monodromy_twin`).  'gear': the former default, a
        ## Gear-2 twin, kept selectable.  'native': the run's OWN plain
        ## factorisation.  ⚠ 'native' UNDER A ONE-STEP METHOD IS THE WORST OF
        ## THE THREE, not an "exact, no twin" escape: trapezoidal's own
        ## monodromy is FIRST order on a limit cycle and its `Q` DIVERGES
        ## under refinement (11.1 / 28.4 / 63.9 vs an exact 5.9083), euler's
        ## likewise -- it is for the gates that measure that defect, not for
        ## results.  gear and trbdf2 runs are self-sufficient (second-order
        ## native monodromy) and ignore this -- they never twin.
        self.monodromy = 'trbdf2'
        self._monodromy_twin = None
        self._twins = {}
        self._solve_kwargs = {}
        ## Set to None by `solve` on the tstab path only -- see there.
        self.tstab_state = None

    def _is_autonomous(self, times):
        """True when nothing in the circuit depends on `t`.

        Exact where a spectral test is not: see `AUTONOMOUS_U_TOL`.  The
        source vector is evaluated at EVERY point of the grid and compared
        with the first; a circuit driven only by DC -- a VCO macromodel, a
        phase accumulator, an LC or ring oscillator -- has a constant `u`
        and a one-parameter family of periodic solutions, which is what
        makes fixed-period shooting ill-posed for it.
        """
        u0 = np.asarray(self.cir.u(times[0], analysis=self.par.analysis),
                        dtype=float)
        scale = max(float(np.max(np.abs(u0))), 1.0)
        ## ⚠ EVERY POINT ON THE GRID, NOT A STRIDE THROUGH IT.  This used to
        ## sample `times[1::len(times)//8]` -- about nine points -- while
        ## calling itself "exact where a spectral test is not".  Nine points
        ## cannot see a narrow pulse: measured on an RC driven by a `VPulse`
        ## positioned BETWEEN two samples, 40% and 20% duty were read
        ## correctly and 5%, 1% and 0.5% all came back AUTONOMOUS.  A clock
        ## misread that way is routed to the free-period system, which
        ## solves for `T` and DISCARDS the period the caller asked for --
        ## and `DEGENERATE_PERIOD_FACTOR` cannot catch it, because it tests
        ## the magnitude of `T`, not whether the circuit was driven.  PWM,
        ## sampling clocks, S/H and mixer LOs are core PSS workload and are
        ## exactly the shapes a stride misses.
        ##
        ## The cost is `N` evaluations of `u` ONCE per solve, against `N`
        ## Newton solves in the traversal it decides -- and the loop exits
        ## at the first sample that differs, which is the common case for
        ## every driven circuit.
        for t in times[1:]:
            u = np.asarray(self.cir.u(t, analysis=self.par.analysis),
                           dtype=float)
            if float(np.max(np.abs(u - u0))) > AUTONOMOUS_U_TOL * scale:
                return False
        return True

    ## The one mapping from `method` to a class.  Read by `_transient` to
    ## build the integrator and by `_companion_reach` to ask how far it
    ## reaches, so the two can never disagree about which method is running.
    @classmethod
    def _integrator_for(cls, method):
        from pycircuit.circuit.integrator import (EulerIntegrator,
                                                  TrapezoidalIntegrator,
                                                  ThetaIntegrator,
                                                  Gear2Integrator)
        from pycircuit.circuit.integrator import (TRBDF2Integrator,
                                                  RadauIIA3Integrator,
                                                  ESDIRK43Integrator,
                                                  GLM2Integrator,
                                                  GLM3Integrator,
                                                  GLM4Integrator)
        ## THE single method -> integrator map, and the one place method names
        ## are validated: an unknown name raises the ValueError here rather than
        ## a KeyError three frames down.  The polymorphic predicates
        ## (`_companion_reach`, `needs_x0_unknown`, ...) call this before the
        ## solve()-level whitelist runs, so the validation must live where the
        ## lookup does.
        table = {'euler': EulerIntegrator,
                 'trap': TrapezoidalIntegrator,
                 'trapezoidal': TrapezoidalIntegrator,
                 ## `theta` is trapezoidal biased by `C h` -- see
                 ## `ThetaIntegrator`. It takes the L-stable opener OUT, which
                 ## is the one thing the other one-step LMMs cannot do.
                 'theta': ThetaIntegrator,
                 'gear': Gear2Integrator,
                 'gear2': Gear2Integrator,
                 'trbdf2': TRBDF2Integrator,
                 'radau': RadauIIA3Integrator,
                 'esdirk43': ESDIRK43Integrator,
                 ## Nordsieck GLMs: stage order = order, one factorisation per
                 ## step, no index-2 order split.  Driven PSS only (no free
                 ## period), and the period map is on the MULTIVALUE state --
                 ## see `factored_period_glm`.
                 'glm2': GLM2Integrator,
                 'glm3': GLM3Integrator,
                 'glm4': GLM4Integrator}
        try:
            return table[method]()
        except KeyError:
            raise ValueError(
                "method must be 'euler', 'trap', 'theta', 'gear', 'trbdf2', "
                "'radau', 'esdirk43', 'glm2', 'glm3' or 'glm4', not %r"
                % (method,))

    def _resolve_break_events(self, requested):
        """`break_events`, ON for every method when not given.

        ⚠⚠ THIS USED TO BE OFF FOR GEAR, ON A ONE-STEP-COUNT MEASUREMENT,
        AND THE LADDER OVERTURNED IT (2026-09-20).  The old table (gear
        uniform 8.23e-3, + events 1.29e-2, jittered 1.24e-2, "lost 7 of 9")
        was one N on the ramped pulse, and its jittered control did show
        that non-uniformity costs a multistep method something.  It does; it
        is just smaller than what the alignment buys.  Pulsed RC, analytic
        reference, max node error, edges ramped over 0.02 T::

            N        gear uniform   gear + events    trap uniform   trap + events
            50       1.39e-02       3.16e-03         7.78e-03       1.68e-04
            100      4.70e-03       5.76e-03         1.50e-03       7.14e-05
            400      3.72e-04       2.44e-04         1.30e-04       5.81e-06
            1600     2.41e-05       1.61e-05         8.13e-06       3.63e-07

        Gear + events wins at 5 of 6 N, by 1.2-4x, and is second order on
        both grids (3.9x per doubling from 400).  Its gain is small next to
        trap's 22x because a two-step formula takes an O(h^2 [x'']) hit at
        the ONE step after a corner, where its history straddles the jump
        in x'' -- measured 5e-7 -> 1.7e-4 across that step, 30x trap's --
        which is a constant, not an order, and is the honest cost of gear on
        hard corners.  At a TRUE jump (tr = 0) the gain is 130x (1.3e-3 ->
        1.0e-5 at 800) once `event_grid` keeps both ends of the clamped ramp
        (see there).  So the default is on for every method, and the old
        `companion_reach() == 1` predicate is gone.

        ⚠ An explicit `True`/`False` is honoured untouched; this only fills in
        `None`.
        """
        if requested is not None:
            return bool(requested)
        return True

    def _resolve_x0_unknown(self, requested):
        """`x0_unknown`, defaulted from the circuit's TOPOLOGY when not given.

        ⚠⚠ WHY THIS IS CONDITIONAL AND NOT A NEW GLOBAL DEFAULT. `x0_unknown`
        is NOT free: trapezoidal still needs an L-stable opener, so switching
        it on moves the Euler step INSIDE the period, where it degrades the
        ORBIT rather than just the opening. Measured on a `Q = 20` resonator
        against its analytic 20 V peak, `x0_unknown` is WORSE --
        20.01273 against 19.76939 at 100 points, 20.02208 against 19.96123 at
        200. Turning it on everywhere would trade a real defect on a few
        circuits for a real regression on most.

        ⚠ ON AN INDEX-2 NETLIST THE TRADE REVERSES, and not marginally. The
        manufactured opening step is INCONSISTENT there: the constraint fixes
        the algebraic variable at a value the step cannot produce, so
        trapezoidal returns EXACTLY 2x on an L-I cutset -- and on an even
        number of steps reports `converged` and a periodicity residual of
        1e-13 while doing it. See the roadmap's section 0k. `x0_unknown`
        removes it at every parity because `x(0)` becomes a genuine unknown.

        Three refusals, each deliberate:

          * **an explicit `True`/`False` is honoured untouched** -- this only
            fills in `None`;
          * **a two-step method is left alone**, because its solved-history
            formulation already solves for `x(0)` and was never affected;
          * ⚠⚠ **a PROVISIONAL verdict does not trigger it, and that is a
            REFUSAL ON THE THEORY RATHER THAN CAUTION.** With a controlled
            source in the loop or cutset the index is not bounded by 2 and
            need not be a function of the topology at all (see
            `topological_index`).  So the premise "the criterion PROVES index
            2, therefore switch" is unavailable — and so is the REMEDY's
            justification, because `x0_unknown` fixes an inconsistent opening
            step on an INDEX-2 algebraic row.  **If the true index is 3 the
            remedy is not known to apply, and switching it on would mask a
            worse problem while reporting a fix.**  Leaving a known defect
            visible is the better failure mode.  The same goes for a
            structurally singular netlist, which has no index at all.

            ⚠ MEASURED 2026-09-08 (docs session, two instruments sharing
            only C and G -- InitDAE eq. 8's 1-fullness of the derivative
            array, and the Kronecker index of the pencil -- 5/5 against
            `topological_index` where that is valid, exact on nilpotent
            pencils of degree 1..5): **index 3 IS reachable on this element
            set, it is VALUE-dependent, and on every provisional fixture
            `topological_index` returned 1 -- not a low-confidence 2.**  A
            VCVS of gain `g` inside a C-V loop (C1 v->b, C2 v->gnd, source
            b->gnd = g*v) is index 3 exactly on `g* = 1 + C2/C1` (six
            (C1, C2) pairs verified) -- the controlled source cancels the
            node's total capacitance, `[C2 + C1(1-g)] dv/dt`, and one more
            differentiation is needed -- and index 2 off it, with
            `|M^3| = 224 |g - 2|`, so a gain within 1e-3 of `g*` still
            carries a ~1 % nilpotent tail: the surface is measure-zero, the
            NEIGHBOURHOOD is what bites.  ⚠ A CCCS half of the same report
            was WITHDRAWN by its author the same day: its fixture put the
            CCCS's ammeter input ACROSS the inductor instead of in series
            (DC-singular at every gain), and rebuilt correctly the L-I
            cutset is genuine, `topological_index` says 2 (provisional), the
            index is 2 at every gain but `F = 1`, and at `F = 1` the row at
            the output node is all-zero so `_structural_singularity` fires
            and `GminAnchorNewton` refuses the rescue in as many words --
            the DC layer catches it loudly; nothing for `ill_posed` to do.
            So: on the VCVS netlists (P1, its self-controlled and
            CCVS-controlled variants) the topological number is 1, wrong by
            two, and `idx != 2` short-circuits before `provisional` is
            consulted -- two independent guards, both firing; and the
            emphasis is INVERTED from "the case that fails to solve": **the
            index-3 netlist DC-solves cleanly and silently at `g*`** (v =
            1e-3 A x 10 Ohm, b = g v, no warning), while the singular one
            is loud.  The case to guard against is the one that solves and
            looks healthy.  What is NOT established: that any of this is
            reachable through the analog blocks people write; the claim is
            about what the element set PERMITS, which is what a refusal has
            to be justified against.  The numerical test (one SVD of an n(k+1) matrix) is
            an OFFER, not built: exact for constant C and G only, so on a
            nonlinear netlist it is the index of the linearisation at one
            operating point, and controlled sources move the index with the
            operating point too.  Code and fixtures:
            `~/docs/.corpus/checks/numindex.py`, `dae_index_probe.py`
            (`--quick` asserts the surface and both pencil cases).

        Warns when it fires, because a silently different formulation is the
        kind of thing that makes a later measurement inexplicable.
        """
        if requested is not None:
            return bool(requested)
        ## ⚠ EVERYTHING BELOW IS BEST-EFFORT AND MUST NEVER RAISE.  This runs
        ## BEFORE `solve` validates its own arguments, so a bad `method` was
        ## reaching `_solves_history` and coming back as `KeyError: 'bogus'`
        ## instead of the `ValueError('method must be ...')` the caller is
        ## owed -- two tests caught exactly that.  A defaulting helper has no
        ## business changing which exception an invalid call raises.
        try:
            if self._solves_history():
                return False
            idx, info = topological_index(self.cir)
        except Exception:
            return False
        if idx != 2 or info['provisional'] or info['ill_posed']:
            self._warn_if_the_block_disagrees(idx, info)
            return False
        where = (('C-V loop: ' + ', '.join(info['loop'])) if info['loop']
                 else ('L-I cutset: ' + ', '.join(info['cutset'])))
        warnings.warn(
            'PSS: this netlist is index 2 (%s), where the manufactured '
            'opening step is INCONSISTENT -- it starts an algebraic variable '
            'at a value the constraint forbids, and trapezoidal carries that '
            'seed forever (exactly 2x on an L-I cutset, reported as CONVERGED '
            'on an even number of steps). Solving for x_0 directly instead; '
            'pass x0_unknown=False to override.' % where,
            RuntimeWarning, stacklevel=3)
        return True

    def _warn_if_the_block_disagrees(self, idx, info):
        """Say so when the TOPOLOGICAL index reads below 2 and the NUMERIC
        algebraic block says otherwise.  Never raises, never changes behaviour.

        ⚠⚠ THE CASE THIS EXISTS FOR IS THE ONE THAT SOLVES AND LOOKS HEALTHY.
        `_resolve_x0_unknown`'s own notes record it: the index-3 VCVS netlist
        "DC-solves cleanly and silently at `g*` ... no warning", while a
        structurally singular one is loud.  A netlist the classifier cannot
        read is declined HERE in silence, and the user is told nothing.

        MEASURED 2026-09-11 on that documented P1 fixture (a VCVS of gain `g`
        inside a C-V loop, index 2 off `g* = 1 + C2/C1` and index 3 on it):
        `topological_index` reads 1 (provisional -- the VCVS is outside its
        covered class), and :func:`algebraic_conditioning` reads the algebraic
        block as SINGULAR at every gain, which is index >= 2.  The two
        instruments disagree and the numeric one is right.  That disagreement
        is the whole signal, and it is free here: this path has already
        decided to decline.

        ⚠ NOT A BEHAVIOUR CHANGE, DELIBERATELY.  `x0_unknown` stays off.  The
        remedy is justified for index 2 and NOT known to apply at index 3, and
        switching it on would mask a worse problem while reporting a fix --
        the reasoning in `_resolve_x0_unknown` above, unchanged.  This says
        what was seen and names the explicit override; it does not take it.
        """
        try:
            if idx is None or idx >= 2:
                return
            _sigma, ac = algebraic_conditioning(self.cir)
            if ac.get('verdict') != 'singular':
                return
            extra = ''
            if info.get('unclassified'):
                extra = (' %d element(s) are outside the classifier\'s covered '
                         'class (%s), so the topological reading is partial.'
                         % (len(info['unclassified']),
                            ', '.join(info['unclassified'][:4])))
            warnings.warn(
                'PSS: the topological index reads %s for this netlist, but its '
                'ALGEBRAIC BLOCK is numerically SINGULAR '
                '(sigma_min(d g_2/d y) = 0), which means index >= 2.%s A '
                'circuit like this solves cleanly and reports nothing, so the '
                'disagreement is the only signal you get. If it has a C-V loop '
                'or an L-I cutset through an element the classifier does not '
                'recognise, the manufactured opening step may be INCONSISTENT; '
                'pass x0_unknown=True explicitly to apply the index-2 remedy, '
                'having checked that the index really is 2 and not 3 -- the '
                'remedy is not known to apply at index 3.'
                % (idx, extra), RuntimeWarning, stacklevel=4)
        except Exception:                                      # noqa: BLE001
            ## ⚠ Best-effort, like everything else on this path: a DIAGNOSTIC
            ## that raises inside a defaulting helper would change which
            ## exception an invalid call reports, which is the defect the
            ## comment above records.  The gates call the helper DIRECTLY so a
            ## bug in it cannot hide behind this.
            return

    def _companion_reach(self):
        """How many charges back the chosen integrator's companion reads.

        The mechanistic property that decides whether this analysis needs
        the entering history as an unknown: a method reaching one charge
        back can be started from a single unknown, one reaching two
        cannot.  Asked of the
        integrator rather than inferred from `method`, so a fourth method
        arrives with the right answer instead of the default one.
        """
        ## Polymorphic: the method answers.  A stage method returns 1 (self
        ## starting, reads only x_n); an LMM computes it from its companion
        ## coefficients.  No isinstance/name branch to extend per method.
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        return integ.companion_reach()

    def _solves_history(self):
        """Whether the period map needs the entering history as an unknown.

        MEASURED, NOT ASSUMED (`benchmarks/pss_seam_cost.py`).  A method
        whose companion reads only `q_{n-1}` cannot see the fabricated
        opening history at all -- euler's seam costs 5.1e-12 V and
        trapezoidal's 1.3e-11 V, both zero -- so enlarging their system
        would quadruple the shooting solve to fix nothing.  Gear-2 reads
        `q_{n-2}`, which in the plain formulation is the entering stand-in,
        and pays 1.266e-01 V at 100 points/period: 54% of its total error,
        rising to 73% at 400 as the seam falls one order slower than the
        interior.

        Autonomous runs take it too, through the composed system (the pair's
        free-period residual in `solve`): a free period does not remove
        the need for a history the companion can read.  See 4c in the class
        docstring for what the seam does to an oscillator, which is NOT what
        it does to a driven circuit -- 0.75% of the error against 54%,
        landing in the orbit's shape rather than its frequency.
        """
        ## ⚠ TRAPEZOIDAL CANNOT JOIN THIS FORMULATION BY SOLVING FOR
        ## `x_{-1}`, and trying it is how that was learned.  `iq_{-1}` is
        ## exactly derivable from `x_{-1}` -- `-(i(x_{-1}) + u)`, item 4d --
        ## but the derivative that matters runs the other way: a one-step
        ## companion reads ONLY `iq_{-1}`, so the trajectory depends on
        ## `x_{-1}` solely through it, and `d(iq_{-1})/d x_{-1} = -G` is
        ## SINGULAR wherever a node carries no conductance -- every purely
        ## reactive node, which is most of a resonator.  Adding `x_{-1}` as
        ## m unknowns then makes the 2m x 2m system rank-deficient, and it
        ## fails exactly as it should: `LinAlgError: Singular matrix`, on 25
        ## tests at once.
        ##
        ## The right second unknown for a `b != 0` method is `iq_{-1}`
        ## ITSELF -- the `(x, iq)` state its monodromy already uses -- with
        ## the closure `iq_{-1} = iq_{N-1}`.  See item 5's note; not built.
        return self._map_kind() == 'pair'

    def _map_kind(self):
        """Which period map the PSS's method shoots on -- decided HERE, once
        (2026-09-23: `solve`, `factored_period` and the state-event check
        each decided it in their own words):

        * 'plain': a linear-multistep companion reaching one charge back --
          one entering unknown (`_walk_lmm`);
        * 'pair':  one reaching two (Gear-2) -- the solved history ``(x_0,
          x_{-1})`` (`_walk_lmm`; see `_solves_history` for why);
        * 'stage': a self-starting Runge-Kutta method (`_walk_stage`);
        * 'glm':   a Nordsieck multivalue method (`_glm_period_blocks`).

        Asked of the integrator, never inferred from the method's name.
        `_solves_history` is this decision's answer for 'pair', and this is
        the one place to override it: the tests that run Gear-2 on the old
        plain formulation do it here (`_force_plain_map`)."""
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        if integ.is_stage_method():
            return ('glm' if getattr(integ, 'is_multivalue', lambda: False)()
                    else 'stage')
        return 'pair' if integ.companion_reach() >= 2 else 'plain'


    def find_initial_solution(self, period, x0=None, npts=60, method=None,
                              eps_rel_lin=1e-2, eps_abs_lin=1e-3, n_iter=7,
                              max_periods=200, zeta=1e-7):
        """A PROPER initial solution to start shooting from, by pre-integrating
        until the fixed-point iteration has entered its LINEAR region.

        De Luca, Bolcato & Schilders, *Proper Initial Solution to Start Periodic
        Steady-State-Based Methods*, IEEE TCAS-I 2019 -- their Algorithm 2.  The
        paper is in `~/docs/07-shooting-methods/`.

        Shooting-Newton needs a start inside its contraction region, and the
        usual remedy is to GUESS a number of pre-integration periods; if the
        guess is wrong the solve diverges and the guess is repeated with no clue
        for the next one.  This detects the handoff point instead, from
        quantities the integration already produces.

        ⚠ WHAT THE CRITERION COMPARES -- AND WHAT IT IS NOT.  It is NOT "the
        iterate stopped moving" and NOT "a carried probe settled": this project
        measured that guess and it moves the WRONG WAY (drift 1.4e-2 while the
        solve still fails, 1.2e-1 once it succeeds -- `benchmarks/
        pss_warm_start.py`), because a settled probe only says the Jacobian
        stopped changing, which is equally true at an equilibrium.  The paper
        compares TWO SEQUENCES: the LINEAR prediction of the next shooting error
        against the one the ACTUAL nonlinear integration produces.  They agree
        only where the fixed-point map really has become linear, which is
        exactly the region a Newton-type method needs::

            u_k         = x_k - phi(x_k)                        (eq. 4)
            u_{k+1}     = J_phi(x_khat) u_k                     (eq. 12)
            utilde_{k+1} = x_{k+1} - phi(x_{k+1})               (eq. 13)

        accepted, componentwise, when (eq. 16)::

            |u_{k+1,j} - utilde_{k+1,j}| <= eps_rel_lin |u_{khat,j}| + eps_abs_lin

        holds for ``n_iter`` CONSECUTIVE iterations (the paper's defaults, used
        for every experiment in it: ``eps_rel_lin=1e-2``, ``eps_abs_lin=1e-3``,
        ``n_iter=7``).  The scale on the right is the shooting error at the
        DETECTION index ``khat``, not at the current one.  A failed check resets
        the run AND re-freezes ``J_phi`` at the new index, which is why ``khat``
        can move.

        ⚠ NON-AUTONOMOUS ONLY, AND THAT IS THE PAPER'S SCOPE, NOT AN OVERSIGHT
        HERE.  Its title, abstract and index terms all say non-autonomous, and
        the reason bites: for a FORCED circuit the DC point is not a fixed point
        of the period map, so "the map became linear" can only mean the orbit.
        For an AUTONOMOUS oscillator the equilibrium IS a fixed point of the
        period map and the map is linear in a neighbourhood of it, so this
        criterion will happily certify the trivial root -- which is precisely the
        van der Pol failure recorded in `benchmarks/pss_warm_start.py`.  **That
        case is NOT solved by this method and must not be handed to it.**

        ``J_phi(x_khat) u`` is taken by the paper's own alternative, the
        directional derivative of eq. (15),
        ``[phi(x + zeta u) - phi(x)] / zeta``, rather than by its Alg. 1
        left-product.  Both are in the paper; this one costs one extra period
        integration per iteration and buys freedom from the opening-frame
        question (which state a stored factorisation is the derivative *about*),
        a distinction that has already cost this file one wrong answer.

        Returns ``(x, info)``: the reduced state to start shooting from -- the
        paper's line 21, "the last computed x_k" -- and a dict with ``khat``,
        ``periods``, ``found`` and the per-period ``history``.  ``found=False``
        means ``max_periods`` ran out with no linear region; the returned ``x``
        is then simply the last iterate and carries no promise.
        """
        T = float(period)
        n = self.cir.n
        iref = self.irefnode
        m = n - 1
        npts = int(npts)
        if npts < 1:
            raise ValueError('find_initial_solution: npts must be >= 1, got %r'
                             % (npts,))
        if n_iter < 1:
            raise ValueError('find_initial_solution: n_iter must be >= 1')
        times = np.linspace(0.0, T, npts + 1)
        hs = np.diff(times)
        integ = method if method is not None else getattr(self.par, 'method',
                                                          'euler')

        def phi(xr):
            """One period of the inner transient from `xr` -- the map the
            shooting residual is built on, run on a fixed uniform grid."""
            tr_saved = getattr(self, '_tran', None)
            self._tran = self._new_transient(self._integrator_for(integ))
            try:
                self._want_dfdh = False
                self._want_lte = False
                self._begin_period(np.asarray(xr, dtype=float))
                x = copy(np.asarray(xr, dtype=float))
                for j, t in enumerate(times[1:]):
                    x = copy(self.solve_timestep(x, t, hs[j]))
                return np.asarray(x, dtype=float).ravel()
            finally:
                self._tran = tr_saved

        if x0 is None:
            x = np.zeros(m)
        else:
            x = np.asarray(x0, dtype=float).ravel()
            if x.shape[0] == n:
                x = np.concatenate((x[:iref], x[iref + 1:]))

        phi_x = phi(x)
        k = 0
        i_iter = 0
        khat = 0
        x_khat = x.copy()
        phi_khat = phi_x.copy()
        u = x - phi_x
        u_khat = u.copy()
        history = []
        found = False
        while k < int(max_periods):
            if i_iter == 0:
                ## re-freeze the linear generator at the current index (Alg. 2
                ## lines 4-8): khat moves whenever the run of successes breaks.
                khat = k
                x_khat = x.copy()
                phi_khat = phi_x.copy()
                u = x - phi_x
                u_khat = u.copy()
            x_next = phi_x
            phi_next = phi(x_next)
            u_tilde_next = x_next - phi_next            ## eq. (13)
            ## eq. (12) via the paper's eq. (15) directional derivative, with
            ## the step scaled to the iterate so `zeta` is a RELATIVE size.
            un = float(np.linalg.norm(u))
            if un <= 0.0:
                u_next = np.zeros(m)
            else:
                z = zeta * max(float(np.linalg.norm(x_khat)), 1.0) / un
                u_next = (phi(x_khat + z * u) - phi_khat) / z
            gap = np.abs(u_next - u_tilde_next)
            ok = bool(np.all(gap <= eps_rel_lin * np.abs(u_khat) + eps_abs_lin))
            k += 1
            i_iter = i_iter + 1 if ok else 0
            history.append({'k': k, 'khat': khat, 'ok': ok, 'run': i_iter,
                            'gap': float(np.max(gap)),
                            'shooting_error': float(np.max(np.abs(u_tilde_next)))})
            x, phi_x, u = x_next, phi_next, u_next
            if i_iter >= int(n_iter):
                found = True
                break
        return x, {'khat': khat, 'periods': k, 'found': found,
                   'history': history}

    ## Set by `solve`; declared here so a caller may read them on a PSS that
    ## has not solved yet, and so `event_times` is never a stale leftover.
    break_events = False
    event_times = []

    def solve(self, refnode=gnd, period=1e-3, x0=None, timestep=1e-6,
              maxiterations=20, grid=None, matrix_free=False,
              x0_unknown=None, tstab=None, break_events=None,
              phase_rule='frozen', state_events=True):
        """Solve for the periodic steady state.

        ⚠ STATE EVENTS AS NEWTON UNKNOWNS (2026-09-21, the event half of
        B7).  `break_events` lands the SOURCE discontinuities, known before
        the solve.  A comparator's, a threshold switch's or a latch's edge
        happens where the SOLUTION crosses a condition -- `Circuit.
        state_events()`, `row . x = threshold` -- at a time no frozen grid
        can hold, and with the crossing inside a step every method was
        FIRST order (voltage-mode PWM loop, radau uniform 100 .. 800 points:
        mean error 2.9e-2 -> 3.3e-3 V, halving per doubling; trbdf2 and
        gear the same).  With `state_events` (the default) a DRIVEN solve
        under radau or trbdf2 runs a second stage: the crossings of the
        first stage's orbit become unknowns `theta_k` (fractions of the
        period) alongside `x_0`, the grid between consecutive events
        scales with its segment (a proportional column per event, the
        period column's own algebra), each event contributes the row
        `row . x(theta_k T) = threshold`, and the bordered Newton lands the
        grid on the crossing exactly.  The solved fractions become the
        grid every consumer replays on and the event nodes are breaks for
        the period quadrature.  `state_events=False` keeps the one-stage
        solve.  Gear runs the stage on its pair of unknowns `(x_0, x_{-1})`
        (2026-09-22): its step derivative carries the previous step's
        partial, assembled from `residual_dT` and `residual_dh`, and the
        source's motion once (`residual_dh` already holds `du/dt`); the
        plain one-step kinds skip the stage with a warning.  Measured on
        the PWM loop: gear staged 5.9e-3 / 4.0e-3 of the swing at 100 /
        200 points against 1.1e-2 / 4.4e-3 unstaged -- its own
        second-order error on the switch-off decay dominates, as trbdf2's
        does; radau is the method for these circuits.
        On an AUTONOMOUS circuit the period joins the unknowns -- `z =
        [x_0, theta, T]`, the period one more column of the event algebra
        (`d h_j / d T = fraction_j`), the phase row closing the system
        (2026-09-22, phase B); the polish convention is proportional.
        Measured on a comparator relaxation oscillator with a 0.2 mV window
        (a 0.04 ns ramp on the slope), against a staged radau-1600
        reference: unstaged uniform grids read a period error of +4.4e-3 /
        +2.0e-4 / -1.8e-3 / -7.8e-4 at 100 / 200 / 400 / 800 points -- set
        by where the crossing falls in its step, no convergence -- and the
        staged solve -6.4e-4 / -1.9e-4 / -6.8e-5 / +8.1e-7, a thousandfold
        at 800.  (With a 20 mV window the unstaged solve was already fifth
        order: a 4 ns ramp is a third of a step and three collocation
        points resolve it -- the stage is for crossings sharper than the
        grid.)

        `break_events` lands the circuit's source discontinuities on grid
        points (`event_grid`).  `None` -- the default -- decides from the
        METHOD: on for a one-step method, off for a multistep one, because it
        HELPS the first and HURTS the second.  See `_resolve_break_events` for
        the measurement and the control that pins the cause.  A circuit whose
        sources declare no discontinuity is untouched bit-for-bit either way.

        `grid` is RECORDED SCOPE ITEM 5: a sequence of step FRACTIONS of the
        period, summing to 1, used in place of the uniform `timestep` grid.
        Fractions rather than absolute times because an autonomous period is
        an unknown and every step has to scale with it.  See
        `_period_grid`; the grid is still frozen for the whole solve, so the
        shooting Newton stays exact.

        `tstab` runs a TRANSIENT for that many seconds before shooting and
        uses its final state as the seed -- the stabilisation time every
        commercial PSS offers, and the standard answer to a seed that is not
        close enough.  `None` (the default) shoots from `x0`, or from the
        operating point, exactly as before.

        ⚠ IT IS THE REMEDY FOR THE FAILURE THIS ANALYSIS FAILS MOST OFTEN.
        Seeded near the unstable DC point -- the trivial-root basin, which is
        where an unseeded autonomous run starts -- van der Pol does not
        solve at all, and one period of `tstab` fixes it:

              circuit                    without tstab      periods needed
              mu = 1  (strongly attracting)  LinAlgError            1
              mu = 0.05 (high-Q)             not converged         ~24

        The count is the `1/mu` amplitude-envelope constant, so it is a
        property of how strongly the limit cycle attracts and NOT of how bad
        the seed is: from 4x and even 20x the orbit amplitude, one period
        suffices at `mu = 1`.  A high-Q oscillator needs proportionally
        more, which is the usual guidance stated as a number.

        ⚠ AND THE STOPPING POINT IS THE CALLER'S, DELIBERATELY.  De Luca,
        Bolcato & Schilders (2019) give a criterion for detecting the handoff
        automatically, and the probe it rests on was measured here and does
        NOT identify it -- near the DC point the monodromy is nearly
        constant, so the probe settles into its own eigenvector and reports
        convergence while the state is stuck at the trivial root.  Every
        obvious substitute shares that defect, because the trivial root IS a
        fixed point of the period map and passes every periodicity test.  So
        the number is asked for rather than guessed.

        ⚠ AND THE COMMERCIAL AUTOMATIC CRITERION IS REPORTED NOT TO SURVIVE
        Q EITHER.  A commercial simulator does offer one; this tree's user reports from
        their own practice that it does not work properly on circuits of
        even moderate Q and does not work on high-Q ones.  That is field
        experience and not a measurement, and it is recorded as such --
        but it points the same way as the probe measured here, which
        reported "settled" while the state was still in the trivial-root
        basin, and it is the same axis: the harder the oscillator attracts,
        the longer the approach and the easier it is for a detector to stop
        early.  Measured above, `mu = 0.05` needs ~24 periods and FIVE is
        not enough.

        ⚠ AND IT CANNOT ESCAPE AN EQUILIBRIUM IT IS STARTED ON.  With
        `x0=None` the seed is the operating point, which on an autonomous
        circuit is an equilibrium -- a transient started exactly there never
        leaves, so no amount of `tstab` helps.  The pre-integration needs
        somewhere to go: pass an `x0` off the equilibrium, or an `ic` on a
        device.  That limit is why this is not a substitute for the probe
        technique, which pumps energy in precisely so the solve cannot fall
        to the DC point.

        ⚠ THE SAME STATEMENT EXISTS IN THE LITERATURE, in words, and the
        measurement above is its quantitative form.  Kundert
        (*Introduction to RF Simulation*, v2 2003; relayed, cited not
        verified here) on starting shooting from a plain transient: "this
        is sufficient to get convergence even on troublesome circuits
        EXCEPT WHEN THE TIME CONSTANTS IN THE CIRCUIT ARE MUCH LARGER THAN
        THE PERIOD OF THE SIGNAL."  That is exactly the mu = 1 -> 1 period
        against mu = 0.05 -> ~24 periods above, and the reason the count
        tracks the `1/mu` envelope rather than the seed: the envelope time
        constant IS the time constant he names.

        See `benchmarks/pss_warm_start.py` for the probe that failed and the
        counts that decided the interface.

        ⚠ `x0_unknown` DEFAULTS TO `None`, WHICH MEANS "DECIDE FROM THE
        TOPOLOGY": it is switched ON automatically for a netlist the index
        criterion proves is index 2, where the manufactured opening step is
        inconsistent, and left OFF otherwise -- because it is not free (see
        the measured regression below). An explicit `True` or `False` is
        honoured untouched. See `_resolve_x0_unknown`.

        `x0_unknown` solves for `x_0` itself instead of for `x_in`, the
        pre-image of a manufactured opening step.  The plain path's default
        is to manufacture `x(0)` with one order-dropped Euler step and hand
        `fsolve` that step's PRE-IMAGE, while the Jacobian it returns is
        taken with respect to `x_0` -- a frame error, and the reason the
        true `dF/dx_in` is SINGULAR (rank 1/3, 2/4, 1/3 on three circuits,
        `sigma_min` exactly 0) and the iteration is a CONTRACTION with a
        linear rate rather than a Newton.  With this set there is no
        manufacturing step, the unknown is the period's own start, and the
        Jacobian is exact.

        ⚠ IT IS NOT FREE, AND THE TRADE-OFF IS GRID-DEPENDENT.  Trapezoidal
        still needs an L-stable opener -- without one its period map is
        `A_trap^K`, singular at EVEN K on every MNA circuit (the
        `(-1)^n` obstruction, now a theorem; see item 4d) -- so the Euler
        step moves INSIDE the period, where it degrades the ORBIT and not
        just the opening.  Measured on a Q=20 resonator against its analytic
        20 V peak:

              npts    default (x_in)    x0_unknown
               100       20.01273        19.76939
               200       20.02208        19.96123

        18x worse at 100 points, 1.8x at 200; the order is preserved and the
        gap closes as the single first-order step is diluted, but the
        constant is real.

        ⚠ AND ON A HARD CASE IT IS THE OTHER WAY ROUND, which is why this is
        an option and not a repair.  Van der Pol on its 1105-step LTE grid
        reaches -47.3 ppm in four iterations against the default path's
        -73.8 -- exactly the figure the throwaway driver in
        `benchmarks/pss_lte_grid.py` has been recorded as "the target to
        hit" since item 5 was written, and for the same reason: its unknown
        was `x_0` too.  Convergence is quadratic at even and odd point
        counts alike (2.9e+00 -> 9.1e-08 -> 1.1e-14).

        So: uniform grid on a benign circuit, leave it off; non-uniform grid
        on a stiff one, turn it on, where the quadratic rate is what decides
        whether it converges at all.  `benchmarks/pss_x0_unknown.py` is the
        driver those numbers come from.  Euler is unchanged to the digit
        either way -- its manufacturing step IS an Euler step, so the two
        maps coincide, which is what makes it the control.

        `matrix_free` is RECORDED SCOPE ITEM 6: solve the outer system
        without ever forming the monodromy, propagating ONE vector per
        Krylov iteration instead of `2m` columns per step.  It is worth
        asking for on LARGE circuits only -- measured on the RC ladder it
        LOSES at m=40 (0.82x) and wins above roughly m=250: 1.40x at m=242,
        2.23x at m=502, 3.62x at m=1002.  ⚠ It buys that with memory, `2 N
        m^2` doubles of stored factorisations (~800 MB at m=1002, 50
        points), and it does NOT produce a monodromy, so `spectral_radius`
        is `None` after a matrix-free solve -- not forming that matrix is
        the entire point.

        ⚠ EVERY FIGURE BELOW HAS A DENSE BASELINE ON BOTH SIDES, and that
        is a property of how they were taken rather than a choice.  Until
        2026-09-02 this path reached past `linearsolver=` to
        `scipy.linalg.lu_factor`, and the dense propagation reached past it
        to `toolkit.linearsolver`; both go through the caller's solver now,
        so the comparison can be re-asked with a sparse one -- and it has
        NOT been, because the box was shared with another session
        benchmarking at load 31-43 when it was tried.  A circuit Jacobian at
        m=1002 is very sparse, so the m~250 gate and these ratios may move.
        `benchmarks/pss_matrix_free_sparse.py` is the harness; run it quiet
        before quoting any of this as the sparse answer.

        ALL FOUR SYSTEMS are converted.  Measured end to end,
        single-threaded, against a DENSE-solver dense path:

          driven solved-history (2m columns)   1.36x/1.51x/2.13x  m=242/502/1002
          driven plain          (m columns)    1.10x/1.34x/1.63x  m=242/502/1002
          autonomous plain      (m + border)   1.02x/1.11x/1.34x  m=128/308/608
          autonomous composed   (2m + border)  1.22x/1.65x        m=208/408

        The plain path's share is about HALF the solved-history path's --
        42.5% against 63.8% at m=502 -- which is what `m` columns instead of
        `2m` against the same assembly has to mean, and it is why the two
        `2m` systems win earlier and by more.  Every one of them agrees with
        the dense path to <= 2e-16, and both autonomous ones reproduce the
        period exactly.

        `phase_rule` chooses how an AUTONOMOUS solve removes the orbit's
        time-translation freedom.  `'frozen'` (the DEFAULT) picks `k` and
        its pinned value once at the seed.  `'reselect'` is Aprille &
        Trick's own Step 3: at every Newton iterate the pinned coordinate
        is `k = argmax |dphi/dT|` -- the vector field at the period's end --
        and it is held at the iterate's OWN value, so the Newton step moves
        every other coordinate and the period (`dz[k] = 0`).  A driven
        circuit ignores it.

        ⚠ RESELECT IS OPT-IN, AND THE REASON IS MEASURED, NOT CAUTION -- see
        the costs below.  Ask for it when a seed is far from the orbit and
        the default reports non-convergence; `tstab` is the other remedy and
        they compose.

        ⚠ WHAT IT BUYS, measured 2026-09-16 on the shipped
        residuals: a frozen pin is a VALUE the orbit must attain, and a far
        seed names one it never reaches.  Van der Pol (mu = 1), six seeds
        per amplitude on a circle, reaching the true period:

              seed amplitude     frozen     reselect    (trap and radau alike)
              on the orbit        6/6         6/6       same answer, <= iterations
              4x                  0/6         6/6
              10x                 0/6         2/6
              30x                 0/6         0/6

        and a Q = 8 tank with a series loss resistor (an algebraic node) at
        4x: 2/6 -> 6/6 under trap, gear (the solved-history system) and
        radau.
        Matrix-free, the same van der Pol seeds: 4x 0/6 -> 4/6 (trap and
        gear), 10x 2/6 (trap) and 0/6 (gear) -- the inexact inner solve
        keeps less of the gain than the dense one.  ⚠ Far seeds can still
        reach the `T = 0` root; `_free_period_solve` refuses it as before.

        ⚠⚠ AND WHAT IT COSTS -- THE FULL SUITE FOUND BOTH, which is why the
        default did not move.  "Every frozen success stays a success" was
        asserted here from two fixtures and is FALSE:

          * **A discontinuous period map stops converging.** `Idtmod` with
            the wrap landing exactly ON a grid point -- the case documented
            to converge anyway -- solves frozen and fails reselected (1200
            points, `ic = 0`; the off-grid wrap at 500 points is fine).  The
            seed choice is not the difference: both rules pick the same `k`
            there, so it is the moving pin that wanders on a map with no
            derivative.
          * **It lands on a different PHASE of the same orbit**, and a
            phase-sensitive surface then reads differently: on the slow-node
            van der Pol (`tau = 100 T`) the orbit matches to the digit
            (period 6.656833399, `v` in +-1.9985, `w` in +-3.238e-3) while
            `frequency_aware_ppv`'s mode content at 0.1 f0 reads 1.64e-6
            against the frozen 2.45e-6.  Anything whose value depends on
            where `t = 0` sits on the orbit moves with it.  ⚠ That surface
            is phase-specific under the FROZEN rule too -- six seeds around
            the same orbit give 2.4522 / 2.3161 / 1.9423 / 1.8306 / 2.2461 /
            2.4445 e-6 -- so the test pinning 2.45e-6 measures ITS seed's
            phase; re-selection simply lands just outside that spread.

        A bit-level consequence of the same mechanism: dense and matrix-free
        Newtons pick `k` from `dphi/dT` values that differ in the last ulp,
        so `lambda_2` agrees to ~1e-13 rather than to the bit.

        ⚠ WHERE THE CONVERGED PERIOD DIFFERS BETWEEN THE RULES, IT IS THE
        GRID, NOT THE RULE: trap with `x0_unknown=True` moves its Euler step
        inside the period, so its discrete orbit depends on the start phase
        (6.663745 / 6.662827 / 6.663525 from on-orbit seeds at three angles,
        identical under either rule), and a frozen solve re-seeded at a
        reselected answer reproduces its period to 1.8e-15.

        ⚠ IT IS NOT A NEW FORMULATION, and the substitution it is usually
        described with buys nothing by itself.  Aprille & Trick write the
        unknown as `[x_01, ..., T, ..., x_0n]` with `x_0k` a constant; with
        `k` FROZEN that solve took the bordered solve's iterates to <= 1e-9
        and failed on exactly the same seeds.  The gain is the re-selection.
        Holding `x_0k` at the iterate's own value is the bordered row with a
        zero residual, which is how it is built here -- inside the residual,
        a pure function of the iterate, so the damped Newton's carried trial
        evaluation stays consistent.

        ⚠ A 2026-09-02 attempt at this was measured and REVERTED (it failed an
        on-orbit seed), and the reason recorded then -- "the row carries no
        information, the orbit slides" -- was wrong.  Reproduced since: a
        pin taken from the iterate's UNKNOWN (`x_in`) but compared against
        the MANUFACTURED `x_0[k]` fails the on-orbit seed exactly so (40
        iterations, not converged), while the consistent row converges in 8;
        and a re-selecting callback fed the damping's cached trial row fails
        every 4x seed the consistent one solves.
        """
        ## ONE SHOOTING SOLVE, IN ITS PHASES (2026-09-24: this was one
        ## 1650-line body).  Each phase is a method with its own record; they
        ## share the run's state through `run`.
        run = self._solve_prepare(refnode, period, x0, timestep,
                                  maxiterations, grid, matrix_free,
                                  x0_unknown, tstab, break_events,
                                  phase_rule, state_events)
        self._shoot(run)
        self._report_convergence(run)
        X, walk, lte_seen = self._replay_orbit(run)
        self._report_lte(run, lte_seen)
        self._check_fundamental(X, walk, run.period)
        tpss, fpss = self._orbit_results(run, X)
        polished = self._closing_polish(run)
        if polished is not None:
            return polished
        return InternalResultDict({'tpss': tpss, 'fpss': fpss})

    def _solve_prepare(self, refnode, period, x0, timestep, maxiterations,
                       grid, matrix_free, x0_unknown, tstab, break_events,
                       phase_rule, state_events):
        """`solve`, phase 1: validate, build the grid, decide autonomy and
        the period column, seed (operating point, `tstab`), pin the phase,
        choose the formulation.  Returns the run's state."""
        self._solve_kwargs = dict(refnode=refnode, maxiterations=maxiterations,
                                  matrix_free=matrix_free, tstab=tstab,
                                  period_seed=float(period))
        self._monodromy_twin = None
        self._twins = {}
        ## ⚠ HIDDEN STATE IS REFUSED, NOT INTEGRATED AND HOPED OVER.
        ## `TLine.history` is filled by `cir.accept_step`, which the
        ## TRANSIENT calls at every accepted step and which this analysis
        ## never calls -- PSS drives `solve_timestep` directly.  With the
        ## buffer empty `TLine.G`/`u` stamp a DC SHORT, so the line is
        ## silently absent: measured on a quarter-wave open stub, PSS
        ## returned `converged = True`, `spectral_radius = 0.0`, NO warning,
        ## and an amplitude of 0.999969 where a transient gives 0.244201.
        ##
        ## ⚠ AND FILLING THE BUFFER IS NOT THE FIX.  Calling `accept_step`
        ## per step would populate it and make `phi` genuinely
        ## history-dependent -- so the monodromy would be the derivative of
        ## a neighbouring problem, which is the exact failure `_begin_period`
        ## exists to prevent.  `_begin_period` resets what is IN `x`; that is
        ## the right scope for the integrator rings and the wrong one for
        ## state living outside the vector, and no reset of the rings can
        ## fix a period map that is not a function of `x_0`.  The honest
        ## answers are to admit the delay state into the unknowns (a
        ## different analysis) or to refuse; this refuses.
        ##
        ## ⚠ AND THE SCOPE OF THAT REFUSAL IS THE ELEMENT, NOT THE CLASS.
        ## Two things look alike and are not.  KUNDERT'S HIDDEN STATE is a
        ## behavioural model carrying internal state the simulator does not
        ## know about -- genuinely broken, and a commercial RF simulator "outlaws [them]
        ## outright".  A DISTRIBUTED COMPONENT has a KNOWN
        ## infinite-dimensional structure described by frequency-dependent
        ## Y/Z/S parameters, and is tractable: "the convolution operation is
        ## diagonalized by the Fourier transform", so the component is
        ## applied spectrally while the STATE stays finite and lumped
        ## (Yang & Phillips, DAC 2002), and an autonomous time-domain
        ## steady-state solve with transmission lines and exact period
        ## derivatives exists independently.
        ##
        ## `TLine` here trips the first test because of HOW IT IS
        ## IMPLEMENTED -- a `history` buffer filled by `accept_step` -- not
        ## because a transmission line is unshootable.  The flag is opt-in
        ## per element (`Circuit.hidden_state` defaults False), so nothing
        ## refuses distributed components as a class, and an element that
        ## declared its state properly would pass.  Cited, not verified
        ## here.
        _hidden = self.cir.hidden_state_elements()
        if _hidden:
            raise NotImplementedError(
                'PSS: these elements carry HIDDEN STATE -- %s -- so THIS '
                'formulation cannot solve this circuit. The period map must '
                'be a function of x_0 alone, and they stamp from state that '
                'lives outside x and is filled by accept_step, which only a '
                'forward transient calls. Left alone the answer would be '
                'silently wrong rather than slow: an empty TLine history '
                'stamps the line as a DC short and the solve reports '
                'converged. Use Transient for this circuit. ⚠ This is a '
                'limit of the ELEMENT as implemented here, not of shooting '
                'or of distributed components as a class: a component with '
                'a KNOWN frequency-domain description (a transmission line, '
                'an S-parameter block) is tractable in a time-domain '
                'steady-state solve by either admitting the delay state '
                'into the unknowns, or applying the component spectrally -- '
                'the Fourier transform diagonalises the convolution, so its '
                'action becomes a multiply by Y/Z/S while the state stays '
                'finite. Both are different analyses than this one.'
                % ', '.join(sorted(_hidden)))

        self.period = period
        toolkit = self.toolkit

        ## ⚠ ONE REFERENCE NODE PER ANALYSIS, CHECKED.  `self.irefnode` is
        ## fixed in `__init__` from `irefnode=` and is what the TRAVERSAL
        ## uses -- `_transient.irefnode`, every `remove_row_col`, the
        ## monodromy's shape.  This local one comes from `solve`'s own
        ## `refnode=` and is what reinserts the zero row into the RESULT.
        ## They were never compared, so `PSS(cir).solve(refnode=b)` solved
        ## against ground and reported against `b`: each row sensible on its
        ## own, the set of them incoherent, with ground itself coming back
        ## non-zero.  Refused rather than silently rotated, because there is
        ## no answer to give -- the two choices disagree about which
        ## variable was eliminated before the solve began.
        irefnode = self.cir.get_node_index(refnode)
        if irefnode != self.irefnode:
            raise ValueError(
                'PSS: solve(refnode=...) names a different reference node '
                '(index %d) than the analysis was constructed with (index '
                '%d). The traversal eliminated one and the result would '
                'reinsert the other, so the waveform would be reported '
                'against a node the solve never used -- ground itself comes '
                'back non-zero. Pass the same node to both, or construct '
                'the analysis with PSS(cir, irefnode=...) and leave '
                "solve()'s refnode at its default."
                % (irefnode, self.irefnode))
        ## ⚠ CLEARED BEFORE THE RUN, not after it.  These describe the
        ## period this call is about to solve for; leaving the previous
        ## call's behind would let `factored_period()` hand back an operator
        ## for the LAST solve after this one failed, and `converged` alone
        ## would not catch it -- a second solve that fails leaves the first
        ## solve's `converged=True` nowhere in sight but its state very much
        ## in reach.
        self._period_state = None
        self._factored_period_cache = None
        self.waveform = None

        ## ⚠ `theta`'s bias is PER-PERIOD (`_theta_biased`), so it must be
        ## known before anything builds or reuses the inner transient.
        ## `_new_transient` only runs when there is no cache, so a SECOND
        ## `solve()` at a different period would otherwise silently keep the
        ## first one's bias -- re-bias the cached integrator in place rather
        ## than dropping the cache, which would rebuild a `Transient` per
        ## solve for every method that does not care.  A no-op for all of
        ## them (`_theta_biased` type-checks, and tolerates `None`).
        self._theta_period = float(period)
        _tr_cached = getattr(self, '_tran', None)
        if _tr_cached is not None:
            self._theta_biased(getattr(_tr_cached.par, 'integrator', None))

        ## Everything `grid_error` needs to repeat THIS solve on a finer grid.
        ## Recorded rather than re-derived so the refinement differs from the
        ## original in the timestep and in nothing else.
        self._solve_args = dict(refnode=refnode, period=period, x0=x0,
                                timestep=timestep,
                                maxiterations=maxiterations, grid=grid,
                                matrix_free=matrix_free,
                                x0_unknown=x0_unknown, tstab=tstab,
                                break_events=break_events,
                                phase_rule=phase_rule)

        n = self.cir.n
        dt = timestep
        if x0 is None:
            x = toolkit.zeros(n-1) #currently without reference node !
        else:
            x = x0 # reference node not included !


        #create vector with timepoints and a more fitting dt
        ## ⚠ the flag must be set BEFORE the grid is built, because
        ## `_period_grid` consults it to decide whether to subdivide a
        ## coarse opening step -- see the note there.
        ## ⚠ `None` means "decide from the topology" -- see
        ## `_resolve_x0_unknown`.  Resolved to a concrete bool HERE, before
        ## anything reads it, so every downstream use sees one value.
        x0_unknown = self._resolve_x0_unknown(x0_unknown)
        if self._integrator_for(getattr(self.par, 'method', 'euler')
                                ).needs_x0_unknown():
            ## Self-starting stage methods: `x_in` IS `x_0`, there is no
            ## manufacturing step to differentiate `x_0` back through, so the
            ## unknown is always `x_0` itself.  The method says so
            ## (`needs_x0_unknown`), which keeps the phase pin and every
            ## open-at-x0 branch consistent without a name check here.
            x0_unknown = True
        self._open_at_x0 = bool(x0_unknown)
        ## Break the traversal's steps at the source discontinuities, when the
        ## method is one whose accuracy that helps -- see `_resolve_break_events`
        ## for the measurement, and `event_grid` for the snap that keeps it from
        ## manufacturing slivers.
        self.break_events = self._resolve_break_events(break_events)
        if self.break_events:
            _ev = (self.event_grid(period, grid=grid) if grid is not None
                   else self.event_grid(period, npts=int(period / dt)))
            ## ⚠ ONLY replace the grid when there ARE events.  `event_grid`
            ## rebuilds a uniform grid from `linspace` even when it finds none,
            ## and that differs from `_period_grid`'s own in the last bit --
            ## enough to move every event-free solve in the suite for nothing.
            ## Touching the grid only when an event exists keeps every circuit
            ## without one BIT-IDENTICAL, the same guarantee `_fold_periodic`
            ## gives a circuit with no periodic state.
            if self.event_times:
                grid = _ev
        times, hs = self._period_grid(period, int(period / dt), grid)
        npts = len(times)
        self._grid_fracs = (None if grid is None
                            else np.asarray(grid, dtype=float))
        ## read by `_period_grid`, which is called from the residual
        ## closures and so cannot take it as an argument
        self._open_at_x0 = bool(x0_unknown)
        ## The fold gauge, collected once per solve (late-bound moduli are
        ## resolved by now).  See `_fold_periodic` for why the residual needs
        ## it and the Jacobian does not.
        self._periodic_fold = self._collect_periodic_fold()
        self._wrap_jump_warned = False
        alpha = 1

        ## AUTONOMY IS DECIDED BEFORE THE SOLVE, because it decides which
        ## system is solved.  Structural and exact -- see `_is_autonomous`.
        self.autonomous = self._is_autonomous(times)
        ## the state-event stage runs only under the stage kinds -- see the
        ## docstring; say so once when the circuit declares events
        self._state_event_fracs = None
        self._event_columns = None
        if state_events:
            _rows = self.cir.state_events() if hasattr(self.cir, 'state_events') else []
            _method_se = getattr(self.par, 'method', 'euler')
            if _rows and self._map_kind() not in ('stage', 'pair'):
                warnings.warn(
                    'PSS: this circuit declares %d state event(s) (a threshold '
                    'switch or comparator) but the state-event stage is built '
                    'for radau, trbdf2 and gear; under %r the crossing stays '
                    'inside a step and the solve is first order there. Use '
                    "method='radau' (best measured), or pass state_events=False "
                    'to silence this.' % (len(_rows), _method_se),
                    RuntimeWarning, stacklevel=3)
        ## the period-column convention for this solve (see the Parameter)
        _pc = str(getattr(self, '_force_period_column', None)
                  or getattr(self.par, 'period_column', 'auto'))
        if _pc not in ('auto', 'proportional', 'closing'):
            raise ValueError("period_column must be 'auto', 'proportional' "
                             "or 'closing', got %r" % (_pc,))
        ## ⚠ 'auto' IS 'closing' WITH THE PROPORTIONAL POLISH (2026-09-21,
        ## Andreas's call restored on the diagnosis of why it was reversed).
        ## The reversal's evidence -- "closing lands +690 ppm off on the
        ## fold, fails from 2 % off" -- was two things that were not closing:
        ## the polish fired only when the closing step left the
        ## zero-stability bound, so a stretch INSIDE the bound stayed in the
        ## answer (gear on its fold from seeds 0 / 0.5 / 2 % low: +1615 /
        ## +2470 / +3646 ppm at stretches 1.04 / 1.19 / 1.61, proportional
        ## +1431 at every seed; the raw window at 5 %: +3129 ppm at 2.28x,
        ## and its "-33.5 ppm" at 2 % was a cancellation), and the 2 %
        ## failure was the old mid-edge seam whose last step was too short
        ## to absorb the correction.  With the polish unconditional the
        ## closing answer is proportional's at the solved period, seed-
        ## independent to 0.4 ppm (+1408.5 / +1408.6 / +1408.9), and closing
        ## keeps its basin: 5 % low on the fold and 16 % on the window,
        ## where proportional's per-step Newton fails.  'proportional'
        ## stays selectable by name -- one solve, the seed's own grid.
        ## ⚠ ON A CALLER'S GRID FOR AN AUTONOMOUS RUN, as the Parameter
        ## says: on a UNIFORM grid 'auto' is proportional -- the closing
        ## column there makes one step of N absorb the period correction
        ## for no reason, and two uniform-grid gear solves that converge
        ## proportionally (the stall fixture at lambda_2 = 0.9, the
        ## periodic-state fold fixture) did not converge closing.
        self._period_column = ('closing' if (_pc == 'closing' or
                               (_pc == 'auto' and grid is not None
                                and getattr(self, 'autonomous', False)))
                               else 'proportional')
        self._closing_inner = None
        self._closing_warned = False
        phase_k, phase_pin = 0, 0.0
        if phase_rule not in ('reselect', 'frozen'):
            raise ValueError("phase_rule must be 'reselect' or 'frozen', not %r"
                             % (phase_rule,))
        self.phase_rule = phase_rule

        if self.autonomous:
            ## An unseeded autonomous run starts at the origin, which IS a
            ## periodic solution -- the trivial one -- and the free-period
            ## system would sit there just as contentedly as the fixed one
            ## did.  The operating point is the honest default: for a phase
            ## accumulator `ic` pins it on the orbit.
            ##
            ## ⚠ `tstab` RUNS AFTER THIS, NOT BEFORE, and the order is the
            ## whole of what makes it work.  A pre-integration seeded from
            ## `zeros` starts AT the equilibrium of an autonomous circuit and
            ## a transient from an exact equilibrium never leaves it, so the
            ## warm start would return the basin it was asked to escape.
            ## Starting it from the operating point is the honest version of
            ## the same statement.
            if x0 is None:
                from pycircuit.circuit.dcanalysis import DC
                xdc = np.asarray(DC(self.cir, toolkit=self.toolkit).solve().x,
                                 dtype=float).reshape(-1)
                x = np.concatenate((xdc[:irefnode], xdc[irefnode + 1:]))

        if tstab:
            ## ⚠ THE PRE-INTEGRATION IS A PLAIN TRANSIENT, and it has to be:
            ## its whole value is that it is NOT a shooting solve, so it
            ## cannot be captured by the basin that the shooting Newton is
            ## stuck in.  It runs on its own adaptive grid -- `timestep` is
            ## a first step, not an imposed one -- because nothing here
            ## needs `phi` to be a function of `x_0`; that requirement
            ## starts when the shooting does.
            from pycircuit.circuit.transient import Transient
            _xred = np.asarray(x, dtype=float).reshape(-1)
            _xfull = np.concatenate((_xred[:irefnode],
                                     np.zeros(1), _xred[irefnode:]))
            _pre = Transient(
                self.cir, toolkit=self.toolkit, reltol=self.par.reltol,
                iabstol=self.par.iabstol, vabstol=self.par.vabstol,
                nrsolver=self.par.nrsolver,
                linearsolver=self.par.linearsolver, scaler=self.par.scaler)
            _res = _pre.solve(refnode=refnode, tend=float(tstab),
                              timestep=dt, x0=_xfull)
            _last = np.asarray(_res.x, dtype=float)[:, -1]
            x = np.concatenate((_last[:irefnode], _last[irefnode + 1:]))
            self.tstab_state = x

        if self.autonomous:
            ## ⚠ AND IT IS SUFFICIENT IN PRINCIPLE, DEGRADED IN PRACTICE ON
            ## HIGH-Q: the row removes the singularity from the UNIT
            ## multiplier and only that one, so an oscillator whose other
            ## multipliers cluster near 1 gives a bordered system that is
            ## nonsingular and ill conditioned.  Same cause as the
            ## eigen-selection limit in the class docstring; read it there.
            ##
            ## THE PHASE CONDITION pins the coordinate moving FASTEST at the
            ## seed, so the orbit crosses the pinning hyperplane
            ## transversally.  Pin a slow one and the last row of the
            ## bordered Jacobian is nearly parallel to the null direction it
            ## exists to remove, which is a singular system wearing an extra
            ## equation.
            ##
            ## ⚠ THE `argmax` COMPARES VOLTS WITH AMPERES ON PURPOSE, and
            ## the obvious repair is WRONG.  The row can only remove the
            ## orbit's tangent in proportion to `|e_k . fhat|`; one step of
            ## `|dx_k|` IS `|f_k|` up to `h`, so this argmax is exactly
            ## `argmax |e_k . fhat|` -- it maximises the quantity the row
            ## needs.  Normalising each coordinate by its own swing, which
            ## is what the vector's mixed units invite, was MEASURED on a
            ## van der Pol carrying a VCVS-scaled copy of `v` and picks a
            ## row 704x WORSE aligned (1.4e-03 against 1.0000).  The
            ## scaling that lets a large coordinate win the argmax is the
            ## same scaling that makes it dominate `f`; the two cancel.
            ## Pinned by `test_the_phase_pin_compares_units_on_purpose`.
            ##
            ## ⚠ A REVIEW DISPUTED THESE CONDITION NUMBERS AND THEN
            ## RETRACTED, and the reason is worth keeping because it is a
            ## trap this file can fall into again.  The reviewer measured
            ## the pin and the orthogonality row as IDENTICAL (edge 1.000x)
            ## and argued a gap was arithmetically impossible, since
            ## `|fhat[k]| = 0.9999` makes `e_k` and `fhat` nearly parallel.
            ## Their harness took the CODE'S analytic `J` and swapped only
            ## the border row, so both readings were the conditioning of
            ## the SAME operator -- and on the plain path that operator is
            ## in the wrong frame (see item 4b's correction below), so its
            ## defect dominated `cond` in both cases.  Identical numbers
            ## were guaranteed by construction.
            ##
            ## ⚠ THE REASONING WAS ALSO WRONG, INDEPENDENTLY: `cond` is a
            ## function of the WHOLE row, not of its projection on `e_k`.
            ## Alignment 0.9999 makes two rows nearly PARALLEL, not equal --
            ## `e_k` is a unit vector and `fhat` is dense -- and it places
            ## no bound on the conditioning gap.  Switching the phase
            ## condition means solving a DIFFERENT system and living with
            ## ITS conditioning, so the comparison has to build both
            ## systems, which is what the numbers above do.
            ##
            ## ⚠ AND THE UPGRADE THIS INVITED WAS TESTED AND REJECTED.  An
            ## orthogonality (Poincare) row `<x0 - x_ref, f(x_ref)> = 0`
            ## looks strictly better -- it is the flow-aligned row by
            ## construction, and it cannot pin an unattainable VALUE.
            ## Measured against this rule on the case built to break it
            ## (seeded at `v`'s turning point, so the pin sits 1e-3 of the
            ## way into its coordinate's range): both converge, to the same
            ## answer, at every grid tried, and the bordered condition
            ## number is 1.2e3/3.0e2/8.2e1 for the pin against
            ## 2.0e2/6.0e1/3.7e1 for orthogonality at 200/800/3200 points --
            ## a 2-6x edge that never decides anything, and it SHRINKS as
            ## the grid refines.  The row's alignment at the solution stayed
            ## 1.0000 throughout: the tangency this was supposed to induce
            ## never materialised.
            ##
            ## What IS real: `phase_pin` is a VALUE the orbit must attain,
            ## so a seed far off the orbit can pin one outside its range and
            ## the system is then INCONSISTENT rather than merely hard --
            ## measured on van der Pol at mu=1, seeds at 4x/10x/30x the
            ## orbit amplitude pin `v` at -5.66/-9.50/-15.46 against an
            ## orbit range of [-2.01, 2.01], and all three report ordinary
            ## non-convergence.  ⚠ That is NOT an argument for the
            ## orthogonality row: its plane through the same far seed misses
            ## the orbit too (checked).  It is an argument about SEEDS, and
            ## the remaining gap is diagnostic, not formulational.
            self._begin_period(x)
            _x1 = self.solve_timestep(x, times[0], hs[0])
            _x2 = self.solve_timestep(_x1, times[1], hs[0], iq_last=self._iq)
            phase_k = int(np.argmax(np.abs(np.asarray(_x2) - np.asarray(_x1))))
            self.phase_k = phase_k                  # which coordinate is pinned, for diagnosis
            ## ⚠ THE RULE IS CANONICAL, AND IT DOES MIX UNITS -- KEPT ON
            ## MEASUREMENT.  Aprille & Trick's oscillator paper, Step 3:
            ## "select k by |f_k(x^i(T^i))| = max_k |f_k(x^i(T^i))|" -- argmax
            ## of the vector field, which is what an argmax over one step is up
            ## to `h`.  `h` is common to every coordinate and cancels; the
            ## UNITS do not (this note used to say they did): the components
            ## are V/s and A/s.  A peer session found the consequence and it
            ## reproduces: on an LC van der Pol with real units (1 nF, 1 uH)
            ## seeded at v's turning point this pins v, nearly along the flow,
            ## and the swing-scaled bordered Jacobian's condition number is
            ## 120-160 against 3.75-9.3 with the pin chosen by motion relative
            ## to each coordinate's swing (one Newton iteration saved).  ⚠ BUT
            ## THAT RULE WAS BUILT, MEASURED AND REVERTED (2026-09-23): on the
            ## comparator relaxation oscillator, seeded at ten phases of its
            ## exact orbit (unstaged, 40 iterations), the raw rule converged
            ## from 7/10 and the swing rule from 4/10 -- raw pinned coordinate
            ## 1 at nine phases, swing moved it to 2 or 3 at all ten, and the
            ## two fail at different phases (not diagnosed further; Andreas
            ## kept the raw rule on these counts).  Conditioning on a smooth
            ## circuit is not what decides convergence on a switching one.
            ## Both the rule and the decision not to replace it with a
            ## Poincare row have precedent as well as measurement behind them.
            ##
            ## ⚠ WHAT IS NOT CANONICAL IS FREEZING IT.  Their Step 3 sits
            ## INSIDE the iteration (Step 5 returns to Step 1), so `k` and
            ## the pinned value are re-chosen from the CURRENT trajectory
            ## every iterate -- "note that in this method, an initial k and
            ## w_0k are not required".  Pinning `w_0k = x_0k^i`, the
            ## iterate's OWN current value, is attainable by construction,
            ## which removes the failure mode measured below (a far seed
            ## pinning a value the orbit never reaches) structurally rather
            ## than by advice.  ✅ DONE 2026-09-16 as `phase_rule='reselect'`,
            ## OPT-IN (the default stays `'frozen'`) -- `_phase_row`, and the
            ## measurements, gains AND costs, in `solve`'s docstring.
            ##
            ## ⚠ THE PHASE ROW SITS OUTSIDE THE INTEGRATOR ON PURPOSE, and
            ## that placement is load-bearing rather than incidental.
            ## Brachtendorf et al. (TCAD 33(6) 867-878) warn of the
            ## alternative: "adding an algebraic equation transforms the
            ## system of (implicit) ODEs to a system of DAEs.  Transient
            ## methods may run into severe problems when the index of a
            ## system of DAEs is two or higher."  This augments the OUTER
            ## shooting system; the inner integration is unaugmented, so the
            ## phase condition cannot raise the index of the DAE actually
            ## being integrated -- which matters here in proportion to how
            ## badly index 2 already behaves (trap and euler both fail to
            ## converge on a V-source/C/C/R loop where gear does).
            ##
            ## ⚠ THE RULE IS CANONICAL; RE-SELECTING IT WAS TRIED AND
            ## REJECTED, WITH NUMBERS.  Aprille & Trick's oscillator paper
            ## picks `k` by `argmax |f_k(x^i(T^i))|` -- the same quantity an
            ## argmax over one step is, up to an `h` common to every
            ## coordinate -- and their Step 3 sits INSIDE the loop, pinning
            ## `w_0k = x_0k^i`, the iterate's OWN value, so that "an initial
            ## k and w_0k are not required".  That looks like a free repair
            ## for the far-seed failure recorded below.  It is not.
            ##
            ## Built and measured 2026-09-02: re-selecting `k` and the pin
            ## from the current trajectory between outer iterations REGRESSED
            ## the working case.  Van der Pol at mu=1 from an ON-ORBIT seed
            ## went from converged to NOT converged, and far seeds wandered
            ## to periods of -52, -1088 and +110 against a true 6.6633.
            ##
            ## ⚠⚠ WITHDRAWN 2026-09-16: THE REASON RECORDED THEN WAS WRONG, and
            ## so was "taking Step 3 means taking the substitution".  The
            ## zero residual does not make the row empty -- `dz[k] = 0` is
            ## exactly A&T's constraint, and their substitution with `k`
            ## frozen reproduces the bordered iterates to <= 1e-9 (the two
            ## are one linear system).  Re-selection built consistently
            ## converges on-orbit and widens the basin (4x seeds 0/6 -> 6/6)
            ## -- at a cost that keeps it OPT-IN: it fails the grid-aligned
            ## `Idtmod` wrap the frozen pin solves, and it lands on a
            ## different phase of the same orbit.
            ## The 09-02 failure REPRODUCES with a pin taken from the
            ## unknown `x_in` but compared against the manufactured `x_0[k]`
            ## -- the frame error the note below records for the frozen pin.
            ## Built as `phase_rule`; see `solve`'s docstring.
            ##
            ## ⚠ THE PIN MUST BE IN THE UNKNOWN'S OWN FRAME.  `_x1` is the
            ## state one step AFTER the seed, which is the right thing to
            ## pin when the unknown is `x_in` and `x_0` is manufactured from
            ## it -- and the wrong thing when the unknown IS `x_0`.  With a
            ## fine opening step the two are nearly equal and the mismatch
            ## hides; on van der Pol's own LTE grid, where the opening step
            ## is 1.4845 against a median of 4.6e-04, it pins a value the
            ## orbit need never attain and the solve dies with a bare
            ## non-convergence.
            phase_pin = float(np.asarray(x if x0_unknown else _x1)[phase_k])

        ## Resolved here as well as in `solve_timestep`, because the SHOOTING
        ## Jacobian depends on which integrator the inner steps used.
        ##
        ## ⚠ THE NAME IS VALIDATED BEFORE ANYTHING ASKS THE INTEGRATOR A
        ## QUESTION.  `_solves_history` resolves `method` to a class to ask how
        ## far its companion reaches; run first, it turned an unknown name
        ## into a `KeyError` from a dict several frames down, in place of the
        ## `ValueError` this raises.  Two tests caught it, both written for
        ## the class's earlier fall-through defects.
        method = getattr(self.par, 'method', 'euler')
        if method not in ('euler', 'trap', 'trapezoidal', 'theta', 'gear',
                          'gear2', 'trbdf2', 'radau', 'esdirk43',
                          'glm2', 'glm3', 'glm4'):
            raise ValueError(
                "method must be 'euler', 'trap', 'theta', 'gear', 'trbdf2', "
                "'radau', 'esdirk43', 'glm2', 'glm3' or 'glm4', not %r"
                % (method,))

        ## Whether the entering history joins the unknowns.  Decided once,
        ## here, because it chooses which system is solved -- like autonomy,
        ## and after it, since the two are not composed yet.
        solved_history = self._solves_history()
        self.solved_history = solved_history
        xm1_ss = None

        return _SolveRun(
            refnode=refnode, period=period, x=x, dt=dt,
            maxiterations=maxiterations, matrix_free=matrix_free,
            x0_unknown=x0_unknown, phase_rule=phase_rule,
            state_events=state_events, irefnode=irefnode, n=n, times=times,
            hs=hs, npts=npts, alpha=alpha, phase_k=phase_k,
            phase_pin=phase_pin, method=method,
            solved_history=solved_history, xm1_ss=xm1_ss)

    def _shoot(self, run):
        """`solve`, phase 2: THE SHOOTING NEWTON -- the fixed-period or
        free-period system on the method's period map (dense or matrix-free),
        then the state-event stage.  Leaves the solution in `run`."""
        toolkit = self.toolkit
        (n, x, period, times, hs, npts, alpha, irefnode) = (
            run.n, run.x, run.period, run.times, run.hs, run.npts, run.alpha,
            run.irefnode)
        (x0_unknown, solved_history, method, phase_rule, phase_k,
         phase_pin) = (run.x0_unknown, run.solved_history, run.method,
                       run.phase_rule, run.phase_k, run.phase_pin)
        (maxiterations, matrix_free, state_events, xm1_ss) = (
            run.maxiterations, run.matrix_free, run.state_events, run.xm1_ss)

        def _phase_row(x0_vec, tcol):
            """The autonomous phase row at THIS iterate: `(k, residual)`.

            `'reselect'` pins `k = argmax |dphi/dT|` over the `x_0` block at
            the iterate's own value, so the residual is zero and the row
            only fixes the step (`dz[k] = 0`); `'frozen'` compares the seed's
            `k` against the seed's value.  See `solve`'s docstring."""
            if phase_rule == 'reselect':
                tc = np.abs(np.asarray(tcol, dtype=float).ravel()[:n - 1])
                return int(np.argmax(tc)), 0.0
            return (phase_k,
                    float(np.asarray(x0_vec, dtype=float)[phase_k]) - phase_pin)

        def _closing(x0, x_end, M, tms_):
            """The fixed-period system at one iterate: ``F = x_0 - phi(x_0)``
            (folded on the idtmod rows, see `_close_periodic`), ``J = I -
            alpha M``."""
            F = self._close_periodic(x0, x_end, tms_)
            D = np.asarray(toolkit.eye(F.shape[0]))
            return F, D - alpha * M

        def _bordered(F, J, tcol, x0_vec, phase_col):
            """The FREE-PERIOD system: the fixed-period one `(F, J)` bordered
            by the period column and the phase row --

                F = [ F ,  x0[k] - pinned ]
                J = [[ J , -dphi/dT ],
                     [ e_k^T ,  0   ]]

            -- because without a phase condition the system is singular by
            construction (every point on the orbit is a solution, so `I - M`
            has a null direction along it).  `k` is chosen by `_phase_row`
            from `phase_col`; the row pins only the `x_0` block, whatever the
            width of `F` (gear's pair: one phase row still suffices)."""
            w = len(F)
            Jb = np.zeros((w + 1, w + 1))
            Jb[:w, :w] = J
            Jb[:w, w] = -np.asarray(tcol).ravel()
            _k, _r = _phase_row(x0_vec, phase_col)
            Jb[w, _k] = 1.0
            Fb = np.zeros(w + 1)
            Fb[:w] = F
            Fb[w] = _r
            return Fb, Jb
        ## THE SHOOTING JACOBIAN FOLLOWS THE INTEGRATOR'S OWN COEFFICIENTS.
        ##
        ## Backward Euler's per-step sensitivity is
        ##     dx_n/dx_{n-1} = Jf_n^-1 * C(x_{n-1})/h
        ## -- the COMPANION CONDUCTANCE at the previous point, not the raw
        ## capacitance matrix.  The `/h` was missing, and because C is
        ## singular the accumulated product collapsed to EXACTLY ZERO: the
        ## Jacobian handed to fsolve was `I - 0 = I`, so the "shooting
        ## Newton" was plain successive substitution `x0 <- phi(x0)`.  That
        ## converges at the circuit's own per-period decay -- measured on a
        ## Q=20 resonator as 0.855 per iteration against exp(-pi/Q) = 0.8546,
        ## which is how it was found -- and it never reached fsolve's
        ## tolerance, on any circuit, silently.  With the companion
        ## conductance the same resonator converges in FIVE iterations
        ## (2.64 -> 2.6e-2 -> 1.0e-2 -> 1.9e-4 -> 3.9e-5).
        ##
        ## TRAPEZOIDAL'S MONODROMY MUST CARRY `iq` AS WELL AS `x`.  Its
        ## recursion carries `iq` as well as `x`, so the period map is a
        ## function of (x, iq); an x-only monodromy is structurally
        ## incomplete, and measured, using the Euler form for trap converged
        ## SLOWER than no Jacobian at all (0.90 against 0.855 per iteration).
        ##
        ## Differentiating the two recursions together,
        ##
        ##     iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}
        ##     0    = i(x_n) + iq_n + u(t_n)
        ##
        ## gives a propagation of `d(x,iq)/dx0` that costs one extra matrix
        ## product over the Euler form:
        ##
        ##     rhs   = Geq_{n-1} Px + Pq
        ##     Px_n  = Jf_n^-1 rhs
        ##     Pq_n  = Geq_n Px_n - rhs
        ##
        ## Euler is the SAME recursion with `Pq == 0`: its companion carries
        ## no `iq_{n-1}` term, so dF/diq_{n-1} vanishes and the second row
        ## never enters.  One formula, two methods, which is why this is not
        ## a second code path.
        ## (`newton = True` lived here and fed a message branch that is
        ## gone: it was set unconditionally, so the 'successive
        ## substitution' alternative was dead and the claim it selected was
        ## false anyway -- see the non-convergence warning below.)

        ## THE PERIOD MAP, per kind -- the only thing about the Newton that
        ## depends on the method (2026-09-23: it was eight residual closures,
        ## one per kind and driven/free period).
        _kind = self._map_kind()

        def _pmap(z, T, tms_, hs_, want_dT):
            """One period from the unknown `z`: ``(z_0, z_end, M, Mt)`` with
            `M = dz_end/dz_0` and, on request, the period column `Mt`.

            * PLAIN (one-step LMM): `z` is the entering state; with a
              manufacturing step the period opens one step in, so `z_0` is
              the state it opened at.
            * gear's PAIR: `z = (x_0, x_{-1})`, and BOTH close --
              ``F = [x_0 - x_{N-1}, x_{-1} - x_{N-2}]``, the rows of
              ``M = [[A(N-1,0), A(N-1,-1)], [A(N-2,0), A(N-2,-1)]]`` from
              `_traverse_solved_history`.  A two-step companion needs two
              states to be continued, so periodicity of ONE is an
              under-determined statement about the orbit (1.266e-01 V at 100
              points under Gear-2 before this).  ⚠ THE HISTORY POINT MOVES
              WITH T: `x_{-1}` sits at `-T/(N-1)`, and `x_{-1}`, `x_{N-2}` are
              the same phase of the orbit at every `T`, so the residual is
              still right and its `T` column is the propagation to step N-2.
            * STAGE (Radau, TR-BDF2, ESDIRK): self-starting, `z` IS `x_0`,
              `M` the dense stage product (`_traverse_stage`).
            * GLM: the method's own map (the startup at the top of the
              period, then N multivalue steps); ⚠ `M` is APPROXIMATE -- the
              residual is exact, the Jacobian drops the startup's derivative
              (see `_walk_glm`); its period column carries the two
              explicit `T` dependences a multivalue method has.

            ⚠ THE PERIOD COLUMN IS TRACTABLE ONLY FOR AN AUTONOMOUS CIRCUIT:
            the grid is rebuilt at the current `T` (``dh/dT = h/T`` for every
            step, uniform or not) and the stage derivatives carry no time of
            their own."""
            w = self._walk(_kind, z, tms_, hs_, T=T, want_dT=want_dT,
                           open_at_x0=x0_unknown)
            M = w.monodromy()
            if _kind != 'glm':
                ## kept for the checks after the solve: the spectrum is the
                ## only place a free period announces itself.  (Never the
                ## GLM's: its `M` drops the startup's derivative, and its map
                ## acts on the Nordsieck state -- see `_walk_glm`.)
                self._monodromy = M
            return (w.z0(), w.end(), M,
                    w.period_column() if want_dT else None)

        def func(z):
            """The fixed-period system, ``x_0 - phi(x_0) = 0``."""
            z0_, z_end, M, _Mt = _pmap(z, period, times, hs, False)
            return _closing(z0_, z_end, M, times)

        def func_autonomous(zT):
            """The FREE-PERIOD system: unknowns `(z, T)`, the fixed-period
            equations bordered by the period column and the phase row (see
            `_bordered`) -- rebuilt at the CURRENT `T`, which is what keeps
            `dh/dT = h/T` true of every step."""
            z, T = zT[:-1], float(zT[-1])
            tms_, hs_T = self._period_grid(T, npts, self._grid_fracs)
            z0_, z_end, M, Mt = _pmap(z, T, tms_, hs_T, True)
            return _bordered(*_closing(z0_, z_end, M, tms_), Mt, z0_, Mt)

        ## THE SHOOTING RESIDUAL IS IN SOLUTION UNITS, NOT KCL UNITS.
        ## `x0 - phi(x0)` is a difference of SOLUTIONS -- volts on node rows,
        ## amps on branch rows -- so its absolute floor is the `xtol` flavour
        ## (vabstol on nodes, iabstol on branches), not the residual flavour
        ## the transient's Newton uses for `i(x)`.  Getting that backwards is
        ## F6(a)'s defect, and it is easy to walk into here because the
        ## quantity is called a residual.
        _tol = analysis.newton_tolerance_vectors(
            len(self.cir.nodes), len(self.cir.branches),
            self.par.iabstol, self.par.vabstol, self.toolkit)[1]
        (_tol,) = remove_row_col((_tol,), irefnode, self.toolkit)

        ## The shooting criterion, expressed against the transient one.
        _ratio = float(self.par.steadyratio)
        if _ratio < 1.0:
            raise ValueError(
                'steadyratio must be >= 1 (got %g): the period map is only '
                'known to the accuracy of the per-timestep solves, so a '
                'shooting tolerance tighter than reltol asks the outer '
                'residual to resolve its own noise.' % _ratio)
        _shoot_reltol = self.par.reltol * _ratio
        _tol = _tol * _ratio

        ## ⚠ REFUSED RATHER THAN SILENTLY IGNORED.  A caller asking for
        ## matrix-free on a path that has not got it wants the cost model it
        ## implies; quietly taking the dense route would be a performance
        ## surprise with no symptom, which is the shape of defect this tree
        ## has paid for before.
        ## ⚠ REFUSED RATHER THAN IGNORED, like `matrix_free` above.  A
        ## solved-history method already solves for `x_0` and `x_{-1}`
        ## directly and manufactures nothing, so there is no frame to
        ## correct and the flag would be a no-op -- and a no-op flag that
        ## the caller believes changed something is worse than an error.
        if x0_unknown and solved_history:
            raise NotImplementedError(
                'PSS: x0_unknown=True has nothing to change for a two-step '
                "method (method=%r). The solved-history formulation already "
                'solves for x_0 and x_{-1} as real trajectory states and '
                'manufactures no opening step, so its Jacobian is exact '
                'without this. Drop the flag, or use a one-step method '
                "(method='trap' or 'euler') where the manufactured opening "
                'is what this replaces.' % method)

        ## ⚠ THE OUTER NEWTON IS DAMPED, which it was not.  All three
        ## `fsolve` calls took the FULL step with `limiter=None` and no line
        ## search -- a departure from standard practice rather than a
        ## neutral choice: Brachtendorf et al. (TCAD 33(6) 867-878) describe
        ## "shooting, finite difference, or harmonic balance techniques IN
        ## CONJUNCTION WITH A DAMPED NEWTON METHOD" as what is "widely
        ## employed" for limit cycles.  The full step is still tried first
        ## and kept whenever it improves the residual, so a solve that was
        ## converging is unchanged; the halving only runs where the
        ## undamped iteration would have moved uphill.

        ## Find periodic steady state x-vector
        ## Find the periodic steady state: ONE Newton for every kind
        ## (2026-09-23; it was a five-way branch).  The unknown is the entering
        ## state -- gear's PAIR `(x_0, x_{-1})`, seeded `x_{-1} = x_0`, the
        ## old formulation's assumption written down, so a pair run starts
        ## where a plain one does and the comparison is about the SOLVE -- and
        ## on an autonomous circuit the period joins it.  Its row is the phase
        ## condition, in the units of the coordinate it pins (`_tol[phase_k]`),
        ## while the UNKNOWN it adds is a time whose own floor must be a time:
        ## mixing the two is flavour error F6(a) one row further out.
        m = n - 1
        _width = 2 if _kind == 'pair' else 1
        if matrix_free and _kind in ('stage', 'glm'):
            raise NotImplementedError(
                'PSS: matrix-free shooting is not built for %s; its '
                'monodromy is a dense stage product. Drop '
                'matrix_free, or use a one-step LMM.' % method)
        z0 = np.concatenate([np.asarray(x, dtype=float)] * _width)
        tol_z = np.concatenate([_tol] * _width)
        _mf = None
        if matrix_free:
            ## ⚠ THE MONODROMY IS NOT FORMED, so it must not be REPORTED
            ## either: `_monodromy` survives from any earlier traversal and
            ## `spectral_radius` reads it without knowing which run wrote it.
            ## (Only the pair's matrix-free paths cleared it; the plain ones
            ## reported the previous solve's radius.)
            self._monodromy = None

            def _mf_build(zz):
                """RECORDED SCOPE ITEM 6: the Newton's residual and its
                Jacobian as a MAT-VEC, from the factored period (`m` columns
                on the plain map, `2m` on the pair, never formed).  With the
                period an unknown, `dphi/dT` is ONE column independent of the
                Krylov direction, computed once per Newton iteration:

                    J [v; s] = [ (I - M) v - s dphi/dT ; v_k ]

                -- one phase row, pinning the `x_0` block only, as the dense
                system does."""
                if self.autonomous:
                    z, T_ = zz[:-1], float(zz[-1])
                    tms_, hs_ = self._period_grid(T_, npts, self._grid_fracs)
                else:
                    z, T_, tms_, hs_ = zz, period, times, hs
                w_ = self._walk(_kind, z, tms_, hs_, T=T_, dense=False,
                                keep=True, want_dT=self.autonomous,
                                open_at_x0=x0_unknown)
                fp_ = w_.factored(self)
                z0_ = w_.z0()
                Mt_ = w_.period_column() if self.autonomous else None
                F_ = self._close_periodic(z0_, w_.end(), tms_)
                if not self.autonomous:
                    return F_, (lambda v: v - alpha * fp_.matvec(v))
                Mt_ = np.asarray(Mt_, dtype=float).ravel()
                k_, r_ = _phase_row(z0_, Mt_)

                def mv_(w):
                    v_, s_ = w[:-1], float(w[-1])
                    top = (v_ - alpha * fp_.matvec(v_)) - s_ * Mt_
                    return np.concatenate((top, [v_[k_]]))
                return np.concatenate((F_, [r_])), mv_

            def _mf(z0_, ab_, xt_, rt_, mi_):
                return self._matrix_free_newton(_mf_build, z0_, ab_, xt_,
                                                rt_, mi_)
        if self.autonomous:
            zT0 = np.concatenate((z0, [period]))
            abstol_z = np.concatenate((tol_z, [_tol[phase_k]]))
            xtol_z = np.concatenate((tol_z, [1e-15 * period]))
            z_ss, _info, _ier, _mesg = self._free_period_solve(
                func_autonomous, zT0, abstol_z, xtol_z, _shoot_reltol,
                maxiterations, period, solver=_mf)
            self.period = period = float(z_ss[-1])
            z_ss = z_ss[:-1]
            ## the grid follows the solved period; everything downstream --
            ## the replay, the waveform, the DFT -- must use it, or the
            ## answer is reported on a period the solver rejected
            times, hs = self._period_grid(period, npts, self._grid_fracs)
        elif matrix_free:
            z_ss, _info, _ier, _mesg = _mf(z0, tol_z, tol_z, _shoot_reltol,
                                           maxiterations)
        else:
            z_ss, _info, _ier, _mesg = analysis.fsolve(
                func, z0, maxiter=maxiterations, reltol=_shoot_reltol,
                abstol=tol_z, xtol=tol_z, toolkit=self.toolkit,
                full_output=True, line_search=True, floor_detect=True)
        ## the state events as Newton unknowns: a second, bordered stage from
        ## the converged orbit (the stage methods, driven or free period, and
        ## gear's driven pair) -- see `_state_event_stage`
        if state_events and _ier == 1 and not matrix_free and (
                _kind == 'stage' or (_kind == 'pair' and not self.autonomous)):
            (z_ss, _info, _ier, _mesg, period, times,
             hs) = self._state_event_stage(
                _kind, z_ss, _info, _ier, _mesg, period, times, hs,
                maxiterations, _tol, _shoot_reltol, alpha,
                *((_phase_row, phase_k) if self.autonomous else ()))
            if self.autonomous:
                self.period = period
        x0_ss = z_ss[:m]
        if _kind == 'pair':
            xm1_ss = z_ss[m:]
        run.period, run.times, run.hs = period, times, hs
        run.x0_ss, run.xm1_ss = x0_ss, xm1_ss
        run.info, run.ier = _info, _ier

    def _report_convergence(self, run):
        """`solve`, phase 3: the convergence flag, the stalled-step and
        non-convergence diagnostics, the Floquet report."""
        _ier, _info = run.ier, run.info
        maxiterations, method = run.maxiterations, run.method
        self.converged = (_ier == 1)
        ## ⚠ WHY A SOLVE THAT HAD STOPPED MOVING STILL FAILED (see `fsolve`'s
        ## `floor_detect`: counted there, never acted on).  The generic
        ## non-convergence warning cannot tell a solve that is lost from one
        ## that is sitting ON its answer; this one can, and names both causes.
        self.step_floor = (_info.get('step_floor')
                           if isinstance(_info, dict) else None)
        if not self.converged and self.step_floor:
            _sf = self.step_floor
            warnings.warn(
                'PSS: the shooting solve STOPPED MOVING AND STILL FAILED ITS STEP '
                'TEST: the periodicity residual has met its tolerance since '
                'iteration %d, but the Newton step of unknown %d stays at %.1e '
                'against a tolerance of %.1e (%.0fx) and does not contract.  Two '
                'causes look like this. (1) THE ARITHMETIC FLOOR: the step is '
                'rounding -- typically a long traversal carrying an unknown of '
                'very different magnitude, tested on a node near a zero crossing '
                'where only the absolute tolerance is left -- and the waveform IS '
                'the solution; an absolute tolerance (vabstol / iabstol) at or '
                'above ~%.0e ends the solve as soon as it is there. (2) A '
                'SINGULAR I - M: a Floquet multiplier at 1 (an autonomous or '
                'marginally stable circuit solved at a fixed period), so the '
                'periodic solution is NOT UNIQUE and the step wanders along the '
                'null direction; check `spectral_radius`.'
                % (_sf['since'], _sf['index'], _sf['step'], _sf['tol'],
                   _sf['ratio'], _sf['step']),
                RuntimeWarning, stacklevel=3)
        self.shooting_iterations = maxiterations if not self.converged else None
        ## ⚠ AN AUTONOMOUS OSCILLATOR CANNOT BE SOLVED AT A FIXED PERIOD,
        ## and this is the only place it says so.
        ##
        ## A circuit whose oscillation is self-sustaining -- a VCO
        ## macromodel, an LC or ring oscillator, any phase accumulator
        ## driven by a DC source -- has a one-parameter family of periodic
        ## solutions, because rotating the starting point along the orbit
        ## gives another one.  Its monodromy therefore has an eigenvalue at
        ## exactly 1 and `I - M` is singular AT the true period.  Away from
        ## it the orbit does not close at all: measured on the quadrature
        ## phase element, the discretisation precesses by 2.1e-3 rad per
        ## cycle at 100 steps/period (falling as h^2), so the period map is
        ## a rotation by slightly less than 2*pi whose ONLY fixed point is
        ## the origin -- which is what an unseeded run returns, silently.
        ##
        ## Measured either side on that element: |eig(M)| = 0.968 and
        ## sigma_min(I-M) = 2.3e-02 at the nominal period, against
        ## |eig(M)| = 1.000226 and 1.6e-04 at the corrected one.  So there is
        ## no period at which this analysis both has a solution and an
        ## invertible Jacobian, and the answer is not a better seed: it is
        ## the autonomous formulation, which solves for the period jointly
        ## with a phase condition.  Not implemented -- but a run that
        ## returns the origin, or refuses to converge, deserves to be told
        ## why rather than left to look like a tolerance problem.
        rho, self.floquet_multipliers, self.parasitic_roots = \
            self._spectral_report(getattr(self, '_monodromy', None))
        self.spectral_radius = rho
        ## `self.autonomous` was decided before the solve and chose which
        ## system ran; nothing to re-derive here.  It used to WARN at this
        ## point that a self-oscillating circuit could not be solved at all,
        ## which was true of the fixed-period system and is no longer true
        ## of this one -- the period was an unknown and `self.period` holds
        ## what it came to.
        if not self.converged:
            ## ⚠ THIS USED TO BE SILENT.  `fsolve` builds the "No
            ## convergence" message and then discards it whenever
            ## `full_output=False`, which is how this call was written -- so
            ## a shooting solve that never converged returned a
            ## plausible-looking waveform with no diagnostic at all.  It was
            ## non-convergent on EVERY circuit, including a linear RLC whose
            ## answer was visibly close, which is why nobody noticed.
            ## ⚠ THIS MESSAGE USED TO CLAIM A "true Newton", AND THAT IS
            ## FALSE ON THE PLAIN PATH.  `newton` is set True
            ## unconditionally, so the alternative branch was dead and every
            ## non-convergence was reported as a Newton failing.  Measured
            ## 2026-09-02: the true `dF/dx_in` is SINGULAR on every circuit
            ## tried (rank 1/3, 2/4, 1/3; sigma_min exactly 0), so no method
            ## solves a true Newton in the frame the plain path's unknown
            ## lives in -- it is a contraction, and its residual falls
            ## LINEARLY at a constant ratio.  Advising `method='euler'` on
            ## the strength of a distinction that does not exist sent people
            ## sideways.
            ##
            ## The advice that IS backed: the solved-history route has an
            ## exact Jacobian (item 4b) and converges quadratically --
            ## measured 1.69e-01 -> 1.06e-02 -> 3.78e-06 -> 3.48e-12 against
            ## trapezoidal's linear 3.91e-03 -> 3.14e-04 -> 2.66e-05 on the
            ## same circuit.
            ## ⚠ THE GEAR ADVICE IS FOR A DRIVEN SOLVE ONLY (2026-09-16).  On an
            ## autonomous oscillator near a unit second multiplier gear is the
            ## method that STALLS (see `_diagnose_lmm_free_period_stall`), so
            ## recommending it there sent a gear user to gear.
            if getattr(self, 'autonomous', False):
                _advice = ("On an oscillator, method='radau' (the default) is "
                           'the robust choice; see any preceding diagnosis.')
            else:
                _advice = ("Raise maxiterations, or use method='gear', whose "
                           'solved-history formulation has an exact Jacobian '
                           'and converges quadratically where the plain path '
                           'is a contraction with a linear rate.')
            warnings.warn(
                'PSS: the shooting solve did not converge in %d iterations '
                '(method=%r). ⚠ The returned waveform IS STILL A FULL '
                'RESULT -- it is the last iterate, not a periodic steady '
                'state -- so a reader who does not check `converged` gets '
                'an array that looks like an answer and is not. %s'
                % (maxiterations, method, _advice),
                RuntimeWarning, stacklevel=3)
        

    def _replay_orbit(self, run):
        """`solve`, phase 4: replay the CONVERGED period the way the solve
        opened it, collecting the per-step LTE.  Returns ``(X, walk,
        lte_seen)`` -- the states, the `(t, h)` pairs, the LTE readings."""
        (solved_history, x0_ss, xm1_ss, times, hs, period, x0_unknown,
         method) = (run.solved_history, run.x0_ss, run.xm1_ss, run.times,
                    run.hs, run.period, run.x0_unknown, run.method)
        ## THE THIRD LEVEL, MEASURED ON THE WAY OUT.
        ##
        ## ⚠ WHY THE NESTING WORKS AT ALL, which this docstring stated the
        ## shape of and never the reason for.  Kundert (*Introduction to RF
        ## Simulation*, v2 2003; relayed from the docs session, cited not
        ## verified here): "the strong convergence properties of shooting
        ## methods result from its nature AS A MULTILEVEL NEWTON METHOD, and
        ## not from the fact it is a time-domain method.  Indeed, it is
        ## possible to formulate harmonic balance as a time-domain method
        ## yet its convergence properties do not fundamentally change."
        ##
        ## The mechanism is that `phi_T` "is a near linear function ... even
        ## when the underlying circuit is behaving in a strongly nonlinear
        ## fashion, because `phi_T` is evaluated over one period of the
        ## large periodic clock signal" -- the nonlinearity is absorbed by
        ## the INNER transient, which is "a natural continuation method,
        ## quite robust".  So the outer Newton sees a nearly linear map and
        ## the arrangement below is not an arbitrary ordering of three
        ## tests: each level exists because the level inside it has already
        ## made the level outside tractable.
        ##
        ## Three convergence criteria stand between a PSS run and its answer,
        ## and the first two are checked while it runs: the inner Newton
        ## (`i(x) + iq + u` under `reltol/iabstol/vabstol`) and the shooting
        ## Newton (`x0 - phi(x0)` under the same, times `steadyratio`).  Both
        ## ask whether an EQUATION was solved.  Neither asks whether the
        ## equation was the right one -- the discrete period map is not the
        ## continuous one, and the gap between them is truncation error.
        ##
        ## PSS imposes its grid (`h = T/(N-1)`, uniform, N from `timestep`),
        ## so this cannot be a CONTROL signal -- nothing here may shrink a
        ## step, and doing so would change the period map between shooting
        ## iterations and destroy the monodromy.  It is a MEASUREMENT, taken
        ## on the converged solution over the final replay and reported.
        ##
        ## ⚠ IT IS THE LEVEL THAT WAS SILENT, and the one that dominates.
        ## On the Q=20 resonator at 100 steps/period all three integrators
        ## report a converged shooting solve, and their amplitudes are
        ## 8.815 V (euler), 19.766 V (gear2) and 19.990 V (trap) against
        ## 20 V analytic -- a 56% disagreement between two "converged"
        ## answers.  Nothing in the two Newton criteria can see that, because
        ## each integrator solved ITS OWN equations to tolerance.  This
        ## number can: it is `|J^-1 Eg| / (TRTOL (reltol ref + lte_abstol))`,
        ## the quantity a transient would have rejected a step on.
        ## ⚠ THE REPLAY MUST OPEN THE WAY THE SOLVE DID, or the waveform is
        ## not the solution: a plain replay of a solved-history answer
        ## would reintroduce the very seam that formulation exists to
        ## remove,
        ## and the reported amplitude would not be the one the residual was
        ## driven to zero on.
        ## ⚠ AND IT MUST WALK THE SAME (t, h) PAIRS, which it did not.
        ## The two traversals pair them differently: `_traverse_solved_history`
        ## walks `times[1:]` with `hs[_j]`, while the plain `_traverse` takes
        ## the MANUFACTURING step at `(times[0], hs[0])` FIRST and only then
        ## walks `times[1:]` with `hs[_j]` -- so in the plain case the step
        ## AFTER the opening one uses `hs[0]` again, not `hs[1]`.  This
        ## replay set `walk = times` and indexed `hs[min(_j, ...)]`, which
        ## pairs `times[k]` with `hs[k]` from `k = 1` on and is off by one
        ## against the traversal.
        ##
        ## A UNIFORM GRID HIDES IT COMPLETELY -- every `hs` is the same
        ## number -- which is why it survived every uniform test in this
        ## file.  Measured on the Q=20 resonator at 200 points, closure
        ## `|x(T) - x(0)|` of the RETURNED waveform:
        ##
        ##       grid      trap            gear (control)
        ##       uniform   5.61e-13        4.62e-14
        ##       4:1       1.70e-02        5.33e-15
        ##       16:1      4.88e-03        1.78e-14
        ##
        ## with `converged = True` in every row.  Gear closes on every grid
        ## because it takes the solved-history branch, whose pairing was
        ## already right; the plain path returned a waveform that is not the
        ## solution its own residual was driven to zero on.
        ##
        ## Now built as explicit `(t, h)` PAIRS rather than two sequences
        ## indexed in parallel, because the parallel indexing is the bug and
        ## a pair cannot be misaligned by one.
        ## ⚠ WHAT A LATER FACTORED REPLAY NEEDS, and the reason it is kept
        ## HERE rather than inside the Newton.  `_traverse_factored*` runs
        ## inside the matrix-free Newton's `_mf_build` closure, whose `steps` go
        ## out of scope with the closure -- and the last `build` call is at
        ## the last TRIAL iterate, which is the converged one only by
        ## accident.  `PAC` wants the factors of the SOLUTION.
        ##
        ## So this keeps the four things a replay cannot re-derive from
        ## outside (which seed, which grid, which opening) and
        ## `factored_period()` runs the traversal on demand.  That costs one
        ## traversal for a caller who asks and nothing at all for one who
        ## does not -- where retaining `N` factorisations from every solve
        ## would cost `2 N m^2` doubles on every run, which is the memory
        ## trade `_traverse_factored` documents and most callers never want.
        self._period_state = (bool(solved_history), copy(x0_ss),
                              None if xm1_ss is None else copy(xm1_ss),
                              times, hs, float(period), bool(x0_unknown))
        self._factored_period_cache = None

        if solved_history:
            self._install_history(x0_ss, xm1_ss, hs[0], h_prev=hs[-1])
            tr = self._transient()
            X = [np.asarray(x0_ss, dtype=float)]
            walk = list(zip(times[1:], hs[:len(times) - 1]))
        else:
            X = [x0_ss]
            tr = self._begin_period(x0_ss)
            ## the manufacturing step, then the loop -- exactly `_traverse`
            ## ... unless there was no manufacturing step, in which case the
            ## replay opens AT `x_0` and walks the period alone.  Getting
            ## this wrong is the same class of defect as the grid shift
            ## above: a replay that does not reproduce its own traversal.
            walk = list(zip(times[1:], hs[:len(times) - 1]))
            if not x0_unknown:
                walk = [(times[0], hs[0])] + walk
        ## Fresh probe, so `relref='sigglobal'`'s running signal maximum is
        ## the period's, not something an earlier shooting iteration saw.
        tr._lte_probe = None
        ## A stage method has no LMM divided-difference LTE (compute_lte
        ## refuses), and the seam/interior split is a property of a manufactured
        ## opener it does not have -- so the replay collects no per-step LTE for
        ## it, and the three LTE figures below report None (honestly: the
        ## diagnostic does not apply to a self-starting stage method).  The
        ## method says which it is.
        self._want_lte = not self._integrator_for(method).is_stage_method()
        lte_seen = []
        for t, dt in walk:
            x = self.solve_timestep(X[-1], t, dt)
            if self._lte is not None:
                lte_seen.append((float(self._lte), float(t), self._lte_seam,
                                 self._lte_valid))
            X.append(copy(x))
        self._want_lte = False

        return X, walk, lte_seen

    def _report_lte(self, run, lte_seen):
        """`solve`, phase 5: the three truncation-error figures and their
        warning."""
        method, npts = run.method, run.npts
        ## THREE NUMBERS, BECAUSE THEY HAVE DIFFERENT REMEDIES.
        ##
        ## `max_lte` is the INTERIOR per-step peak -- steps whose estimator
        ## saw only real past charges.  It is exactly the quantity a
        ## transient controls its grid on, and it ranks the integrators the
        ## way their answers rank: on the Q=20 resonator at 100 points per
        ## period it reads euler 0.2876, gear2 0.0763, trap 0.0239 against
        ## amplitudes of 8.815 / 19.766 / 19.990 V (analytic 20 V).
        ##
        ## ⚠ BUT THE PEAK IS A PER-STEP NUMBER AND A LIMIT CYCLE IS WHAT A
        ## WHOLE PERIOD DOES.  At `reltol=1e-3` euler's peak is 0.288 -- in
        ## tolerance -- while its amplitude is 56% low, because a transient's
        ## criterion bounds each step and says nothing about the 99 of them
        ## that damp the orbit.  `total_lte`, the SUM over the interior, is
        ## the one that sees it: 26.27 for that run against 0.941 for gear2
        ## and 0.340 for trap, tracking the amplitude errors of 55.9%, 1.17%
        ## and 0.05%.  It is an upper bound -- it adds magnitudes, so it
        ## cannot see the cancellation that makes trapezoidal's real error
        ## far smaller than its summed one -- which is the right direction
        ## for a diagnostic to be wrong in.
        ##
        ## `max_lte_seam` is the peak over the opening steps of a method
        ## whose COMPANION reads the entering unknown -- Gear-2 here, and
        ## `None` for euler and trapezoidal, which cannot have a seam.  See
        ## `solve_timestep` for why that is the right condition and what the
        ## looser one reported.
        ##
        ## ⚠ IT IS A FLAG, NOT A MAGNITUDE.  The seam is real -- measured at
        ## 1.266e-01 V on the Q=20 resonator at 100 points/period, against
        ## an interior contribution of 1.070e-01, and it is the term that
        ## stops converging (54% of Gear-2's error there, 73% at 400 points)
        ## -- but the NUMBER printed overstates it by orders of magnitude,
        ## because the estimator differences a fabricated charge while the
        ## solution merely reads one.  At 100 points the estimator's
        ## seam/interior ratio is 505x and the answer's is 1.18x.  So use it
        ## to know the seam is there; use `benchmarks/pss_seam_cost.py` to
        ## know what it costs.
        ##
        ## The fix is not a smaller timestep -- refining makes its SHARE
        ## grow.  It is to make the entering history part of the shooting
        ## unknowns, so the map is a fixed point in the state a two-step
        ## method actually needs, rather than one that opens off a
        ## stand-in.  Not built; the prize
        ## measured on that resonator is Gear-2's error going 2.34e-1 ->
        ## 1.07e-1 at 100 points, at no extra cost per iteration.
        ## An unsound estimate is reported as neither: for trapezoidal that
        ## step is the ONLY one whose number was ever wrong, and dropping it
        ## is what keeps it out of both figures.
        interior = [p for p in lte_seen if not p[2] and p[3]]
        seam = [p for p in lte_seen if p[2]]
        if interior:
            self.max_lte, self.max_lte_time = max(interior)[:2]
            self.total_lte = float(sum(p[0] for p in interior))
        else:                                            # pragma: no cover
            self.max_lte = self.max_lte_time = self.total_lte = None
        self.max_lte_seam = max(seam)[0] if seam else None

        ## Named so the warning can lead with whichever is actually
        ## limiting: the three have three different answers.
        _limits = [
            (self.total_lte, 'accumulated over the period',
             'use a smaller timestep or a less damping method -- this is '
             'the figure that sets a limit cycle, and a per-step criterion '
             'can be in tolerance while it is not'),
            (self.max_lte, 'in one interior step',
             'use a smaller timestep or a higher-order method'),
            (self.max_lte_seam, 'over the opening steps',
             "this is the period map's own seam, where each shooting "
             'iteration cold-starts from a fabricated history; it does '
             'NOT improve with a smaller timestep'),
        ]
        _over = [(v, where, fix) for v, where, fix in _limits
                 if v is not None and v > 1.0]
        if _over:
            v, where, fix = max(_over, key=lambda r: r[0])
            ## ⚠ DO NOT ASSERT CONVERGENCE HERE.  This clause read "the
            ## shooting solve converged, but ..." unconditionally, and the
            ## LTE report is produced whether or not it did -- so a
            ## non-converged run emitted a warning whose first words said it
            ## had converged, directly beside the warning that said it had
            ## not.  Two reviewers read non-converged waveforms as answers
            ## in this file's history; contradictory warnings are not why,
            ## but they are not help either.
            warnings.warn(
                'PSS: the shooting solve %s, and the periodic '
                'solution is not resolved at this accuracy (method=%r, %d '
                'points per period). Local truncation error reaches %.3g '
                'times tolerance %s: %s. Neither Newton criterion can see '
                'this -- they ask whether the discrete equations were '
                'solved, not whether the discretisation is the right one. '
                '(peak interior %s at t=%.6g s, period total %s, opening '
                'steps %s; relax lte_vabstol/lte_iabstol/TRTOL if this '
                'accuracy is intended.)'
                % ('converged' if self.converged else 'did NOT converge',
                   method, npts, v, where, fix,
                   'n/a' if self.max_lte is None else '%.3g' % self.max_lte,
                   -1.0 if self.max_lte_time is None else self.max_lte_time,
                   'n/a' if self.total_lte is None
                   else '%.3g' % self.total_lte,
                   'n/a' if self.max_lte_seam is None
                   else '%.3g' % self.max_lte_seam),
                RuntimeWarning, stacklevel=3)

    def _check_fundamental(self, X, walk, period):
        """`solve`, phase 6: warn when an autonomous solve returned a
        MULTIPLE of the fundamental (sets `fundamental_period`)."""
        ## ⚠ AN AUTONOMOUS PERIOD IS ONLY DETERMINED UP TO AN INTEGER
        ## MULTIPLE, AND THE SOLVE FOLLOWS THE SEED.
        ##
        ## `k*T` is a period whenever `T` is, so `x0 - phi_{kT}(x0) = 0` has
        ## solutions at every multiple and the free-period system converges
        ## to whichever one the seed is nearest.  Measured on the quadrature
        ## phase element, whose true period is 1.000e-03: seeds of 1e-3,
        ## 2e-3 and 3e-3 return 1.000083e-03, 2.000665e-03 and 3.002245e-03
        ## and ALL report `converged`.  The reported waveform is a correct
        ## periodic solution in each case -- and its fundamental frequency
        ## is wrong by the factor, which is what a PSS user is usually
        ## after.  Nothing said so.
        ##
        ## The detector is cheap and needs no extra solve: an orbit
        ## traversed k times comes back near `x_0` partway through.  Grid
        ## points do not land on `T/k` in general (`T/2` at 199 steps is
        ## step 99.5), so this is a NEAREST-APPROACH test against the
        ## orbit's own diameter rather than an equality, and the endpoints
        ## are excluded because every orbit is near `x_0` there.
        ##
        ## Driven runs are exempt: their period is the caller's, and asking
        ## for two source periods is a legitimate request, not a mistake.
        if self.autonomous and len(X) > 8:
            ## ⚠ THE ORBIT IS ITS DIFFERENTIAL STATE, WITH THE PERIODIC ROWS
            ## FOLDED (2026-09-23).  Two points of an autonomous orbit
            ## coincide exactly when their differential states do -- the
            ## algebraic rows are functions of them -- and a phase coincides
            ## modulo its modulus.  Over the whole vector, an idtmod VCO
            ## could never recur (its unfolded phase state advances one
            ## modulus per fundamental), and its wrapped output jumps by a
            ## whole modulus mid-period.  Measured on the free-running
            ## `VcoHdl`: every one-fold solve warned "19.8 times", and true
            ## 2- and 3-fold orbits were named 19.9 / 8.65 / 11.96.
            _pts = np.array([np.asarray(v, dtype=float).ravel() for v in X])
            _dyn = np.any(np.asarray(self._C_at(_pts[0]), dtype=float) != 0.0,
                          axis=0)
            if not np.any(_dyn):
                _dyn[:] = True
            _dev = _pts - _pts[0]
            for _r, _m, _o in self._periodic_fold:
                _dev[:, _r] -= _m * np.round(_dev[:, _r] / _m)
            _d = np.max(np.abs(_dev[:, _dyn]), axis=1)
            _diam = float(np.max(_d))
            ## ⚠ THE THRESHOLD IS THE ORBIT'S OWN SPEED AT `x_0`.  A k-fold
            ## orbit passes `x_0` again at the SAME speed, between grid
            ## points, so its nearest point is within about half a step's
            ## displacement there -- and an orbit moving away from `x_0`
            ## is two steps' displacement off by the excluded edge.  The
            ## LARGEST step on the orbit is the wrong scale: it was an
            ## output wrap's jump (0.99 on the VCO) or a fast edge, and at
            ## three times that every point in the first quarter of the
            ## diameter passed.  `_h` is each step of the replay, so
            ## non-uniform grids scale locally.
            _h = np.array([h for _t, h in walk], dtype=float)
            _v0 = _d[1] / _h[0]
            _edge = max(2, len(_d) // 20)
            ## ⚠ THE EARLIEST RECURRENCE, NOT THE NEAREST.  A three-fold
            ## orbit passes close to `x_0` at both `T/3` and `2T/3`, and
            ## `argmin` picked whichever happened to be numerically nearer
            ## -- it reported `2T/3` as "the fundamental", which is wrong by
            ## a factor of two and would have sent the reader to a period
            ## that is itself a multiple.
            _near = [j for j in range(_edge, len(_d) - _edge)
                     if _d[j] < _v0 * max(_h[j - 1], _h[j])
                     and _d[j] < 0.25 * _diam]
            if _diam > 0.0 and _near:
                ## The closest approach WITHIN THE FIRST cluster: the first
                ## point over the threshold is up to a step early, which
                ## read 3% low and made the multiple look like 2.06 rather
                ## than 2.00.
                _run = [_near[0]]
                for _c in _near[1:]:
                    if _c != _run[-1] + 1:
                        break
                    _run.append(_c)
                _j = min(_run, key=lambda i: _d[i])
                ## the time from the reference point, summed over the
                ## replay's own steps: with a caller's grid the points are
                ## not evenly spaced, and on the plain path `X[0]` sits one
                ## step before `t = 0`, so `times[_j]` was a step off there
                self.fundamental_period = float(np.sum(_h[:_j]))
                warnings.warn(
                    'PSS: this autonomous solve returned a period that '
                    'is a MULTIPLE of the fundamental. The orbit comes '
                    'back within %.2g of its own diameter at t=%.6g s, '
                    'so the fundamental is about %.6g s and the '
                    'returned %.6g s traverses it about %.1f times. '
                    'k*T solves the periodicity condition whenever T '
                    'does, so the solve follows its seed -- re-run with '
                    'period=%.6g to get the fundamental. The waveform '
                    'is a correct periodic solution either way; its '
                    'FUNDAMENTAL FREQUENCY is what is off.'
                    % (_d[_j] / _diam, self.fundamental_period,
                       self.fundamental_period, period,
                       period / self.fundamental_period,
                       self.fundamental_period),
                    RuntimeWarning, stacklevel=3)

    def _orbit_results(self, run, X):
        """`solve`, phase 7: the reported waveform (`self.waveform`) and the
        two results, `tpss` in time and `fpss` in frequency."""
        toolkit = self.toolkit
        (times, irefnode, solved_history, x0_unknown) = (
            run.times, run.irefnode, run.solved_history, run.x0_unknown)
        ## ⚠ THE PLAIN PATH'S FIRST ENTRY IS A SEED, THE OTHER PATH'S IS A
        ## SOLUTION.  Plain takes N steps from `x0_ss` and reports their
        ## results; the other starts AT `x_0` and takes N-1, so dropping the
        ## first would drop a real point and shift the waveform by a step.
        ## ⚠ THE FIRST ENTRY IS DROPPED ONLY WHEN IT IS NOT PART OF THE
        ## PERIOD.  On the default plain path `X[0]` is `x_in`, the
        ## pre-image of the manufactured step, which sits one step BEFORE
        ## t=0 and is not a point of the orbit.  With `x0_unknown` -- and on
        ## the solved-history path -- `X[0]` IS `x(0)`, so dropping it both
        ## discards a real sample and leaves the waveform one column short
        ## of `times`.
        X = toolkit.array(X if (solved_history or x0_unknown) else X[1:]).T

        # Insert reference node voltage
        X = toolkit.concatenate((X[:irefnode], 
                                 toolkit.zeros((1,len(times))), 
                                 X[irefnode:]))

        tpss = analysis.CircuitResult(self.cir, x=X, xdot=None,
                                      sweep_values=times, sweep_label='time', 
                                      sweep_unit='s')

        ## ⚠ KEPT FOR THE CARRIER PHASOR, which AM/PM needs and which
        ## `fpss` below cannot supply: `freq_analysis` returns an RMS,
        ## energy-folded, positive-frequency spectrum -- right for
        ## reporting and wrong for a phasor, because folding destroys the
        ## phase relationship between a carrier and its sidebands, which is
        ## the entire content of an AM/PM decomposition.
        self.waveform = (np.asarray(times, dtype=float),
                         np.asarray(X, dtype=float))

        freqs, FX = freq_analysis(X[:,:-1], times[:-1])
        ## ⚠ ON A NON-UNIFORM GRID (`grid=`) THE INDEX DFT ABOVE IS NOT A
        ## FOURIER COEFFICIENT -- measured 7.5 % off in the carrier and not
        ## converging on a 3:1 grid.  Same layout and RMS fold, taken as the
        ## trapezoid-weighted sum at the true times (`_period_quadrature`'s
        ## weights; uniform grids never reach this branch).
        _h = np.diff(np.asarray(times, dtype=float))
        if len(_h) >= 2 and float(np.max(_h)) / float(np.min(_h)) - 1.0 \
                > self.UNIFORM_GRID_TOL:
            _tt = np.asarray(times, dtype=float)
            _Tp = float(_tt[-1] - _tt[0])
            ## the SAME rule `_period_quadrature` gives every consumer
            ## (2026-09-21): a periodic cubic spline on an event-free grid, a
            ## piecewise one breaking at the landed event nodes -- `fpss` and
            ## `carrier_phasor` are pinned equal to 1e-12, which a second
            ## rule here broke
            _wq = periodic_spline_weights(_tt[:-1], _Tp, self._event_nodes(_tt[:-1], _Tp)) / _Tp
            _ks = np.arange(len(freqs))
            freqs = _ks / _Tp
            _E = np.exp(-2j * np.pi * np.outer(_ks, (_tt[:-1] - _tt[0]) / _Tp)) \
                * _wq[None, :]
            FX = np.asarray(X[:, :-1], dtype=float) @ _E.T
            FX[:, 1:] *= np.sqrt(2)

        ## ⚠ `fpss` IS RMS, AND THE USUAL THING TO COMPARE IT AGAINST IS NOT.
        ## `freq_analysis` returns an RMS, energy-folded, positive-frequency
        ## spectrum (see the note above).  A commercial simulator's frequency
        ## -domain PSS output is conventionally PEAK, so a harmonic read from
        ## one and compared against the other differs by `sqrt(2)` -- 3.01 dB
        ## -- with nothing in either result announcing it.  Multiply `fpss`
        ## by `sqrt(2)` for a peak-convention comparison, or divide theirs.
        ## Recorded because this campaign has lost time to factor-of-two and
        ## factor-of-pi convention defects more than once, and a 3 dB offset
        ## is small enough to be mistaken for a modelling difference.
        fpss = analysis.CircuitResult(self.cir, x=FX, xdot=None,
                                      sweep_values=freqs, sweep_label='freq', 
                                      sweep_unit='Hz')
        
        return tpss, fpss

    def _closing_polish(self, run):
        """`solve`, phase 8: after a 'closing' free-period solve, solve once
        more proportionally on the caller's fractions from the converged
        state, and return THAT result; None when no second pass is due."""
        (refnode, period, hs, x0_ss, maxiterations, matrix_free, x0_unknown,
         phase_rule) = (run.refnode, run.period, run.hs, run.x0_ss,
                        run.maxiterations, run.matrix_free, run.x0_unknown,
                        run.phase_rule)
        ## ⚠ THE SECOND PASS (2026-09-21).  'closing' keeps a caller's inner
        ## steps where the transient validated them, so the free-period
        ## Newton converges from a seed period 16 % off where 'proportional'
        ## fails its per-step Newton (relaxation van der Pol, mu = 10, its
        ## own `lte_grid`).  But the closing step then absorbed the whole
        ## period correction -- 3 s against a 0.7 s neighbour, a growth far
        ## beyond a two-step method's zero-stability bound -- so the
        ## integrator dropped it to Euler and the transposed replay refused.
        ## Closing is the BASIN device; once converged, the grid is
        ## re-fractioned at the solved period and solved once more
        ## proportionally from the converged state, which is the sane grid
        ## the answer is reported on.  ⚠ ALWAYS, NOT ONLY BEYOND THE BOUND
        ## (2026-09-21).  This used to polish only when the closing step
        ## left the zero-stability bound, on the reasoning that a seed
        ## already close needs no second pass -- but ANY stretch of the last
        ## step is a change of discretisation and stays in the answer: gear
        ## on its fold, seeds 0 / 0.5 / 2 % low, +1615 / +2470 / +3646 ppm
        ## at stretches 1.04 / 1.19 / 1.61 where proportional reads +1431 at
        ## every seed and the polished answer +1408.5 / +1408.6 / +1408.9.
        ## The seed-exact case stretches too (the discretisation's own
        ## period error is absorbed by the last step).  Cost: one Newton
        ## from a converged state, one or two iterations.
        if (getattr(self, '_period_column', 'proportional') == 'closing'
                and self.converged and getattr(self, 'autonomous', False)
                and not getattr(self, '_closing_second_pass', False)
                and len(hs) > 2):
            from pycircuit.circuit.integrator import ZERO_STABILITY_RATIO
            _hs = np.asarray(hs, dtype=float)
            _r = float(_hs[-1] / _hs[-2])
            if _r > ZERO_STABILITY_RATIO or _r < 1.0 / ZERO_STABILITY_RATIO:
                warnings.warn(
                    'PSS: the closing step ended %.2fx its neighbour after the '
                    'free-period solve moved the period from %.6g to %.6g s; '
                    'solving once more on that grid re-fractioned at the '
                    'solved period (proportional), from the converged state.'
                    % (_r, float(self._solve_kwargs.get('period_seed', period)),
                       float(period)), RuntimeWarning, stacklevel=3)
            ## ⚠ ON THE CALLER'S FRACTIONS, not the closing-distorted grid:
            ## re-fractioning THAT grid keeps the giant last step (measured:
            ## trbdf2 -1.1 % in period, c -97 %, second pass or not)
            _fr_caller = (np.asarray(self._grid_fracs, dtype=float)
                          if self._grid_fracs is not None else None)
            self._closing_second_pass = True
            self._force_period_column = 'proportional'
            try:
                return self.solve(refnode=refnode, period=float(period),
                                  x0=copy(x0_ss),
                                  timestep=float(period) / len(_hs),
                                  maxiterations=maxiterations,
                                  grid=_fr_caller,
                                  matrix_free=matrix_free,
                                  x0_unknown=x0_unknown, tstab=None,
                                  break_events=self.break_events,
                                  phase_rule=phase_rule)
            finally:
                self._closing_second_pass = False
                self._force_period_column = None
