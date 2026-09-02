from pycircuit.post import InternalResultDict
from .circuit import gnd
from pycircuit.circuit.analysis import *
from copy import copy
import warnings
import pycircuit.circuit.analysis as analysis
import numpy as np

def freq_analysis(x, t, rms = True, axis=-1, freqoffset = 0):
    """Return dft of equidistant sampled signal x"""
    
    npoints = np.size(x, axis)

    dt = t[1] - t[0]

    if x.dtype in (np.cdouble, np.cdouble):
        X = np.fft.fftshift(np.fft.fft(x, axis=axis),axes=(axis,)) / npoints
        freqs = np.fft.fftshift(np.fft.fftfreq(npoints, d=dt))
    else:
        freqs = np.fft.fftfreq(npoints, d=dt)[:int(np.ceil(npoints / 2.))]
        slices = [slice(None)] * x.ndim
        slices[axis] = slice(0, len(freqs))
        X = np.fft.fft(x, axis=axis)[tuple(slices)] / npoints
        ## Fold energy from negative frequencies
        X[:,1:] *= np.sqrt(2)

    if not rms:
        X *= np.sqrt(2)

    return freqs, X

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


class PSS(Analysis):
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
                                period is not given.  `func_solved_history`'s
                                sibling `func_autonomous`.
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
         (2026-09-01), `func_autonomous_solved_history`.  It is 4b's
         enlargement AND the free-period one at once, because an autonomous
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
      `func_solved_history` does.  The parasitic roots are also the spurious
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

        Neither outcome was silent -- the collapse reports
        `converged = False` and the exception is loud -- but neither named
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

         NOT PROPOSING A FOURTH.  Bordering the `(x, iq)` system to remove
         the alternating mode remains formally available -- the analogy to
         the phase condition is exact -- but the null direction there is a
         property of the discretisation with no closed form to pin, unlike
         `xdot(0)`, and three derivations in this item have now looked sound
         and failed on contact.  Anyone taking it up should measure before
         building: the falsifier is cheap and is even/odd step counts.

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
         `_monodromy_matvec` + `_matrix_free_solve`.  The recursion is not
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
         Parameter(name='vabstol', 
                   desc='Absolute voltage error tolerance', unit='V', 
                   default=1e-12),
         Parameter(name='maxiter', 
                   desc='Maximum number of iterations', unit='', 
                   default=100),
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
                        'TRTOL / Spectre lteratio)', unit='', default=7.0),
         Parameter(name='relref',
                   desc="What the relative LTE tolerance is measured "
                        "against: 'pointlocal', 'alllocal' or 'sigglobal'",
                   unit='', default='sigglobal'),
         Parameter(name='steadyratio',
                   desc='Shooting tolerance as a multiple of reltol (>= 1); '
                        '1 holds the shooting solve to the same relative '
                        'tolerance as the transient, larger relaxes it',
                   unit='', default=1.0),
         Parameter(name='method',
                   desc="Integration method for the inner transient: 'trap' "
                        "(default), 'euler', or 'gear' (BDF-2)",
                   unit='',
                   default="trap")]        

    
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
        self._dfdh = None
        self._want_lte = False
        ## The caller's step fractions, or None for the uniform grid.  Read
        ## by the autonomous closures, which rebuild the grid at the current
        ## `T` on every residual evaluation.
        self._grid_fracs = None
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

    def _is_autonomous(self, times):
        """True when nothing in the circuit depends on `t`.

        Exact where a spectral test is not: see `AUTONOMOUS_U_TOL`.  The
        source vector is sampled across the period and compared with the
        first sample; a circuit driven only by DC -- a VCO macromodel, a
        phase accumulator, an LC or ring oscillator -- has a constant `u`
        and a one-parameter family of periodic solutions, which is what
        makes fixed-period shooting ill-posed for it.
        """
        u0 = np.asarray(self.cir.u(times[0], analysis=self.par.analysis),
                        dtype=float)
        scale = max(float(np.max(np.abs(u0))), 1.0)
        for t in times[1::max(1, len(times) // 8)]:
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
                                                  Gear2Integrator)
        return {'euler': EulerIntegrator,
                'trap': TrapezoidalIntegrator,
                'trapezoidal': TrapezoidalIntegrator,
                'gear': Gear2Integrator,
                'gear2': Gear2Integrator}[method]()

    ## Below this fraction of the seed, a solved period is the trivial
    ## root rather than an orbit.  Deliberately loose: a real fundamental
    ## reached from a seed one decade high is ~0.1 of it, and the trivial
    ## root lands 15 orders down, so nothing sits near this line.
    DEGENERATE_PERIOD_FACTOR = 1e-6

    def _free_period_solve(self, func, z0, abstol, xtol, reltol, maxiter,
                           seed_period):
        """Solve a free-period system, with its degenerate root named.

        ⚠ `T = 0` IS A REGULAR ROOT OF EVERY AUTONOMOUS SHOOTING SYSTEM.
        `x0 - phi_T(x0)` vanishes identically at `T = 0`, and the phase
        condition does not exclude it -- it constrains `x0`, not the
        period -- so Newton reaches it from any seed below the fundamental
        and the run returns a period of ~1e-18 with no orbit in it.

        Measured on BOTH autonomous elements in the tree, so it is a
        property of the formulation and not of any circuit: from a 1e-4
        seed against a 1e-3 fundamental, Gear-2 returned -1.5e-20 on the
        quadrature element and 3.9e-19 on the scalar `Idtmod`, and
        trapezoidal raised a bare `LinAlgError` from three seeds of five as
        its Jacobian went singular on the way down.

        Neither outcome was a silent wrong answer -- the collapse reports
        `converged = False` and the exception is loud -- but both told the
        user nothing about the cause, and the generic non-convergence
        advice ("raise maxiterations") is wrong for it: no number of
        iterations reaches a fundamental from below.
        """
        try:
            z, info, ier, mesg = analysis.fsolve(
                func, z0, maxiter=maxiter, reltol=reltol, abstol=abstol,
                xtol=xtol, toolkit=self.toolkit, full_output=True)
        except np.linalg.LinAlgError as exc:
            raise np.linalg.LinAlgError(
                'PSS: the free-period Jacobian went singular while solving '
                'for an autonomous period seeded at %.6g s (%s). The usual '
                'cause is a seed BELOW the fundamental: `T = 0` solves the '
                'periodicity condition identically, so the iteration is '
                'drawn to it and the Jacobian degenerates on the way. Seed '
                'at or above the expected period -- a short transient and '
                'the interval between two output recurrences is the usual '
                'way to get one.' % (seed_period, exc)) from exc

        T = float(z[-1])
        if not np.isfinite(T) or abs(T) < self.DEGENERATE_PERIOD_FACTOR * abs(
                seed_period):
            warnings.warn(
                'PSS: this autonomous solve collapsed onto the TRIVIAL root, '
                'returning a period of %.6g s from a seed of %.6g s. `T = 0` '
                'satisfies `x0 - phi_T(x0) = 0` identically and the phase '
                'condition does not exclude it (it constrains x0, not the '
                'period), so a seed below the fundamental is drawn there. '
                'The returned waveform is not a periodic steady state. '
                'Raising maxiterations will not help -- seed at or above the '
                'expected period instead; a short transient and the interval '
                'between two output recurrences gives one.'
                % (T, seed_period), RuntimeWarning, stacklevel=3)
        return z, info, ier, mesg

    def _period_grid(self, period, npts, grid):
        """`(times, hs)` for one period -- uniform, or a caller's own grid.

        RECORDED SCOPE ITEM 5.  A transient adapts because it cannot see
        the future; PSS re-solves the SAME interval over and over, so it can
        be handed a grid that was chosen well ONCE and then frozen.  The
        grid still never moves inside a solve, so the shooting Newton stays
        exact -- freezing is what makes it a Newton, and this changes only
        WHICH frozen grid.

        ⚠ THE RECORDED BLOCKER WAS STALE.  Item 5 said this was "blocked on
        `Transient` accepting a non-uniform grid; `fixed_timestep` is
        uniform-only".  `Transient.solve`'s loop is uniform-only and always
        was -- but PSS never uses that loop.  It drives `solve_timestep`
        directly, one step at a time, and non-uniform steps worked through
        that path unchanged.  Verified before any of this was written.

        `grid` is a sequence of step FRACTIONS of the period, summing to 1.
        Fractions rather than times because an autonomous period is an
        unknown: every step must scale with `T`, or `dh/dT = h/T` -- the
        identity the period column rests on -- stops holding.

        Measured on van der Pol at mu=100 (`benchmarks/pss_lte_grid.py`):
        1105 steps taken from an adaptive transient converge where 1105
        UNIFORM steps do not, and beat a 20000-point uniform grid on
        accuracy -- 18x fewer points and -47.3 ppm against -60.6.
        """
        if grid is None:
            times, dt = self.toolkit.linspace(0.0, period, num=npts,
                                              endpoint=True, retstep=True)
            return times, np.full(len(times), float(dt))
        fr = np.asarray(grid, dtype=float)
        if fr.ndim != 1 or len(fr) < 2:
            raise ValueError('grid must be a 1-D sequence of at least two '
                             'step fractions, got shape %r' % (fr.shape,))
        if not np.all(fr > 0.0):
            raise ValueError('every grid step fraction must be positive; '
                             'the smallest given is %g' % float(fr.min()))
        total = float(fr.sum())
        if abs(total - 1.0) > 1e-9:
            raise ValueError(
                'grid step fractions must sum to 1 (they are fractions of '
                'the period, so that every step scales with T when the '
                'period is an unknown); they sum to %.12g' % total)
        ## ⚠ THE OPENING STEP IS MANUFACTURED, SO IT MUST NOT BE THE GRID'S
        ## COARSE END.  `_traverse` builds `x(0)` from the unknown with ONE
        ## order-dropped step of `hs[0]`, and a grid taken from an adaptive
        ## transient opens wherever that transient's window happened to
        ## start -- on van der Pol at mu=100, `h[0] = 1.4845` against a
        ## MEDIAN of 4.62e-04, 3200x coarser.  That single Euler step moves
        ## the state 7.4%, and the shooting Newton then has to invert a map
        ## whose first act is that step.  Opening at the grid's own finest
        ## step instead costs ONE extra step in 1105 and is the difference
        ## between not converging and converging.
        ##
        ## ⚠ THE 8x IS MEASURED, NOT DERIVED, and it is a guard rather than
        ## a threshold: it exists so grids that already work are left
        ## exactly as the caller wrote them ('2:1' opens at 2x its finest,
        ## 'smooth' at 5x, and neither needs this).  The falsifier is in
        ## the record: on van der Pol's grid, an opening step of 1e-1 still
        ## fails for gear and 1e-2 converges, against a ratio here of
        ## 13939.  Anything between those bounds separates the two cases.
        if fr[0] > 8.0 * fr.min():
            d = float(fr.min())
            fr = np.concatenate(([d, fr[0] - d], fr[1:]))

        hs = fr * period
        times = np.concatenate(([0.0], np.cumsum(hs)))
        return times, hs

    def _companion_reach(self):
        """How many charges back the chosen integrator's companion reads.

        The mechanistic property that decides whether this analysis needs
        the entering history as an unknown: a method reaching one charge
        back can be started from a single unknown, one reaching two
        cannot.  Asked of the
        integrator rather than inferred from `method`, so a fourth method
        arrives with the right answer instead of the default one.
        """
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        alphas, _b = integ.companion_coefficients(1.0, 1.0)
        return len(alphas) - 1

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

        Autonomous runs take it too, through the composed system
        (`func_autonomous_solved_history`): a free period does not remove
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
        return self._companion_reach() >= 2

    def _step_sensitivity(self, Px, Cs, Pq, Jf, C_new, solve=None):
        """One step of the sensitivity recursion, for ANY seed width.

        ONE RECURSION FOR EVERY METHOD (and now for either formulation).
        Each writes its companion as `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`,
        so differentiating the step gives

            S    = sum_{k>=1} a_k C_{n-k} P_{n-k} + b Pq
            P_n  = -Jf_n^-1 S
            Pq_n = a_0 C_n P_n + S

        `P` is `d x_j / d(unknowns)`: one block wide in the plain
        formulation, two when the history is solved for.  Nothing here
        depends on that
        width, which is why the two systems share this and not a copy.

        ⚠ A SOLVE, NOT AN INVERSE (stage 11).  `inv(Jf) @ ...` formed a dense
        inverse per timestep per iteration and squared the condition number
        it then multiplied through.

        `solve` overrides how that solve is taken, and exists so the
        MATRIX-FREE path (item 6) can hand in a PRE-FACTORED `Jf` without
        this recursion being copied.  It is the same recursion either way,
        which is the point: a second copy would be a second thing to get
        wrong, and this one is already shared by every method and both
        formulations.  Without it, matrix-free would refactor `Jf` once per
        step PER KRYLOV ITERATION -- `k` times the factorisations the dense
        path takes, which is worse than the problem it set out to fix.
        """
        alphas, b = self._coeffs
        S = b * Pq if b else np.zeros_like(Px[0])
        for k in range(1, len(alphas)):
            S = S + alphas[k] * (Cs[k - 1] @ Px[k - 1])
        if solve is None:
            Px_new = -self.toolkit.linearsolver(Jf, S)
        else:
            Px_new = -solve(S)
        Pq_new = alphas[0] * (C_new @ Px_new) + S
        return Px_new, Pq_new

    def _traverse_solved_history(self, x0_in, xm1_in, times, hs,
                                 T=None, want_dT=False):
        """One period from a SOLVED history, for a two-step companion.

        The plain `_traverse` solves for a single entering state `x_in` and
        manufactures `x(0)` from it with one order-dropped Euler step.  That
        is sound for a companion reaching one charge back and measurably
        wrong for one reaching two: the shooting condition constrains
        `x(0) = x(P)`, it does NOT constrain `x_in` to be the orbit's own
        `x(-dt)`, so `x_in` is an O(h^2) stand-in -- and Gear-2 reads it as
        a history point.

        Here BOTH `x(0)` and `x(-dt)` are unknowns and both are required to
        close:

            F = [ x_{N-1} - x_0 ,  x_{N-2} - x_{-1} ]

        which is periodicity of the whole state a two-step method needs,
        rather than of one slice of it.  The trajectory then opens at full
        order off a history the solve is responsible for, so there is no
        opening step to drop and no fabricated charge to read.

        Returns `(x_last, x_prev, P_last, P_prev)`, the `P` being
        `d x / d(x_0, x_{-1})` as one `n x 2n` block -- and with `want_dT`,
        two more entries: `d x_{N-1}/dT` and `d x_{N-2}/dT`, for the
        autonomous system, where the period is an unknown too.  BOTH rows
        need a period column, which is the one thing the composed system
        needs that neither enlargement carried alone.
        """
        toolkit = self.toolkit
        self._want_dfdh = want_dT
        toolkit = self.toolkit
        m = self.cir.n - 1

        ## THE HISTORY IS INSTALLED, NOT SEEDED.  `_begin_run(x_{-1})` opens
        ## the rings on the earlier point and the push puts `x_0` in front of
        ## it, so the first real step reads `q(x_0)` and `q(x_{-1})` -- two
        ## genuine solved points.  The flags then say what is true of them:
        ## a step of `dt` has been taken (`_dt_last`), the run is no longer
        ## opening (`_is_first_step`, `_no_history`), and `_dt_last2` stays
        ## None because the THIRD charge is still `q(x_{-1})` repeated -- the
        ## LTE estimator differences three, so its opening reading remains
        ## unsound and the report goes on discarding it.
        tr = self._install_history(x0_in, xm1_in, hs[0], h_prev=hs[-1])

        ## `P_0 = [I 0]`, `P_{-1} = [0 I]` -- the two unknowns, exactly.  The
        ## plain path seeds BOTH rings with `I`, which is the flat-history
        ## assumption written into the Jacobian; here there is nothing to
        ## assume.
        eye = np.asarray(toolkit.eye(m))
        zero = np.zeros((m, m))
        Px = [np.hstack((eye, zero)), np.hstack((zero, eye))]
        Cs = [np.asarray(self._C_at(x0_in)), np.asarray(self._C_at(xm1_in))]
        ## `Pq` is `d(iq_{-1})/d(x_0, x_{-1})`.  For `b = 0` the recursion
        ## never reads it and zero is right.  For `b != 0` it is NOT zero,
        ## and it is exactly differentiable because `_install_history` seeds
        ## `iq_{-1} = -(i(x_{-1}) + u)` from the DAE:
        ##
        ##     d(iq_{-1})/d x_{-1} = -G(x_{-1}),   d(iq_{-1})/d x_0 = 0
        ##
        ## Leaving it at zero would make the Jacobian wrong for trapezoidal
        ## in the same way the plain path's flat seed is wrong -- which is
        ## the whole defect this formulation removes.
        Pq = np.zeros((m, 2 * m))

        ## Kept for PAC as the plain path does -- but opening EMPTY, because
        ## on this path `x_0` is an unknown rather than the result of a step,
        ## so there is no solved `(C, Jf)` pair at it to record.  The lists
        ## therefore align with `times[1:]`, one shorter than `X`.  PAC is
        ## withdrawn (`test_PAC_is_withdrawn`) and nothing else reads them;
        ## stated here so its rewrite does not read a stale alignment out of
        ## the plain path's shape.
        self.Cvec = []
        self.Jtvec = []
        self.times = times

        ## The period column, propagated by the SAME recursion with one
        ## extra source term.  Zero at the start: neither unknown depends on
        ## T -- they are states, and the solve owns them.  (In the plain
        ## autonomous system the single entering unknown is likewise
        ## T-independent, so this opens the same way.)
        Pt = [np.zeros(m), np.zeros(m)]
        Pqt = np.zeros(m)

        x, x_prev = copy(x0_in), copy(xm1_in)
        P_prev = Px[1]
        for _j, t in enumerate(times[1:]):
            dt = hs[_j]
            x_prev = x
            x = copy(self.solve_timestep(x, t, dt))
            self.Cvec.append(copy(self._C))
            self.Jtvec.append(copy(self._Jf))

            Jf = np.asarray(self._Jf)
            C_new = np.asarray(self._C)
            Px_new, Pq = self._step_sensitivity(Px, Cs, Pq, Jf, C_new)

            if want_dT:
                ## Every step scales together (`h = T/(N-1)`), so `dh/dT =
                ## h/T`, and `df/dh` at fixed solution is Fang's `p` --
                ## `residual_dh`, already shared.  For an AUTONOMOUS circuit
                ## its `du/dt` half vanishes, which is what makes solving
                ## for the period tractable at all.
                alphas, b = self._coeffs
                St = b * Pqt if b else np.zeros_like(Pt[0])
                for k in range(1, len(alphas)):
                    St = St + alphas[k] * (Cs[k - 1] @ Pt[k - 1])
                St = St + np.asarray(self._dfdh).ravel() * (dt / T)
                Pt_new = -toolkit.linearsolver(Jf, St)
                Pqt = alphas[0] * (C_new @ Pt_new) + St
                Pt = [Pt_new, Pt[0]]

            P_prev = Px[0]
            Px = [Px_new, Px[0]]
            Cs = [copy(C_new), Cs[0]]

        self._want_dfdh = False
        ## ⚠ THE FULL 2m x 2m MAP, not the `d x_{N-1}/d x_0` corner.  For a
        ## two-step method the one-period map acts on the PAIR, and its
        ## spectrum is the physical Floquet multipliers TOGETHER WITH the
        ## parasitic roots the k-step discretisation introduces (see the
        ## literature note in the class docstring -- controlling those roots
        ## is the whole subject).  The autonomous unit-circle eigenvalue is
        ## in there.  ⚠ `spectral_radius` NO LONGER TAKES A MAXIMUM OVER
        ## THIS MIXTURE (2026-09-02): `_spectral_report` separates the
        ## physical multipliers from the parasitic roots by the
        ## eigenvector's block structure and reports the maximum over the
        ## physical ones, so a method whose spurious root sits near the
        ## unit circle cannot have it read back as the orbit's stability.
        ## ⚠ ALWAYS THE FULL 2m x 2m MAP, NEVER THE `d x_{N-1}/d x_0`
        ## CORNER.  This used to hand the corner back on the driven path
        ## (`Px[0][:, :m]`), and a sub-block of a sensitivity is not a
        ## monodromy: its eigenvalues mean nothing.  Measured, it reported
        ## `spectral_radius` 1.279605 for the Q=20 resonator -- ABOVE ONE,
        ## reading as an unstable orbit -- where the analytic per-period
        ## decay is exp(-pi/Q) = 0.854636 and every other path reports
        ## 0.855.  For a two-step method the one-period map acts on the
        ## PAIR, so the monodromy is the pair's.
        ##
        ## Its spectrum carries the parasitic roots of the discretisation
        ## alongside the physical multipliers -- but they are not a problem
        ## here and the reason is quantitative: BDF-2's parasitic root is
        ## 1/3 per STEP (roots of `1.5z^2 - 2z + 0.5` are 1 and 1/3), so
        ## over a period it is (1/3)^N, which at 200 points is ~1e-95.
        ## Measured on the autonomous element the 16x16 spectrum is
        ## [1.000, 3.5e-06, 5.2e-16, 4.7e-17, 0, ...]: one physical unit
        ## eigenvalue and nothing else above rounding, so `max |eig|` picks
        ## the physical one.  A method whose parasitic root sat nearer the
        ## unit circle would need that separated; Gear-2's does not.
        self._monodromy = np.vstack((Px[0], Px[1]))
        if want_dT:
            return x, x_prev, Px[0], Px[1], Pt[0], Pt[1]
        return x, x_prev, Px[0], P_prev

    def _traverse_factored(self, x0_in, xm1_in, times, hs):
        """One period WITHOUT the sensitivities, keeping each step factored.

        RECORDED SCOPE ITEM 6, the trajectory half.  `_traverse_solved_history`
        does two things at once: it walks the period, and it propagates a
        `2m`-column sensitivity alongside.  Matrix-free needs the walk
        without the propagation, because the propagation is the thing it
        replaces -- so this returns what a MATVEC needs to replay the same
        steps: the two opening capacitances, and per step a FACTORED `Jf`
        with its `C`.

        ⚠ THE FACTORISATION IS THE POINT, not an optimisation on top.  The
        matvec runs one solve per step per Krylov iteration; refactoring
        `Jf` each time would cost `k` times the factorisations the dense
        path takes and lose before it started.  Factoring once here makes
        every later solve a back-substitution.

        ⚠ AND THE COST OF THAT IS MEMORY: `N` factorisations and `N`
        capacitances, `2 N m^2` doubles, where the dense path holds `O(m^2)`
        at a time.  At `m = 1002` and 50 points that is ~800 MB.  This is
        the trade matrix-free makes here and it is not free; a caller with a
        long period and a large circuit can run out of memory where the
        dense path merely ran slowly.
        """
        from scipy.linalg import lu_factor
        m = self.cir.n - 1
        self._install_history(x0_in, xm1_in, hs[0], h_prev=hs[-1])
        C0 = [np.asarray(self._C_at(x0_in)), np.asarray(self._C_at(xm1_in))]
        steps = []
        x, x_prev = copy(x0_in), copy(xm1_in)
        for _j, t in enumerate(times[1:]):
            x_prev = x
            x = copy(self.solve_timestep(x, t, hs[_j]))
            steps.append((lu_factor(np.asarray(self._Jf)),
                          np.asarray(self._C).copy()))
        return C0, steps, x, x_prev

    def _monodromy_matvec(self, C0, steps, v):
        """`M v` for the solved-history map, replaying stored steps at width 1.

        THE SAME RECURSION THE DENSE PATH USES -- `_step_sensitivity`, with
        its solve pointed at the stored factorisation.  The dense path seeds
        it with `P_0 = [I 0]`, `P_{-1} = [0 I]` and carries `2m` columns;
        here it is seeded with the two halves of `v` and carries one, which
        is the whole difference between the two paths.

        `v` is `(v_0, v_{-1})`; the result is `(P_last v, P_prev v)`, so the
        shooting Jacobian's action is `v - M v` -- see `_matrix_free_solve`.
        """
        from scipy.linalg import lu_solve
        m = self.cir.n - 1
        v = np.asarray(v, dtype=float)
        Px = [v[:m].copy(), v[m:].copy()]
        Cs = list(C0)
        Pq = np.zeros(m)
        for lu, C_new in steps:
            Px_new, Pq = self._step_sensitivity(
                Px, Cs, Pq, None, C_new, solve=lambda S: lu_solve(lu, S))
            Px = [Px_new, Px[0]]
            Cs = [C_new, Cs[0]]
        return np.concatenate((Px[0], Px[1]))

    ## How hard GMRES is asked to solve, relative to the shooting tolerance.
    ## An inexact Newton only needs the step accurate enough not to spoil
    ## the outer convergence; measured on the RC ladder, k is 2/4/7/12 at
    ## m=40/110/242/502 against this factor and the whole system's spectrum
    ## explains why -- `I - M` has almost every eigenvalue within 1% of 1.0
    ## because the fast modes decay to nothing over a period, leaving only
    ## the slow ones for GMRES to resolve.  So k tracks the number of SLOW
    ## MODES, not m, which is the property the item rests on.
    KRYLOV_TOLERANCE_FACTOR = 1e-2

    def _matrix_free_solve(self, z0, times, hs, abstol, xtol, reltol,
                           maxiter):
        """The outer Newton with the monodromy never formed (item 6).

        The dense path builds the `2m x 2m` Jacobian and factors it once per
        iteration; here the same iteration runs on a matvec, so the
        `2m`-column propagation never happens.  Measured against the dense
        path on the RC ladder, single-threaded, k=12:

              m     dense traversal   trajectory + 12 matvecs   speedup
             40             0.0843                    0.1025      0.82x
            110             0.2366                    0.2175      1.09x
            242             0.7503                    0.5378      1.40x
            502             3.4709                    1.5457      2.23x
           1002            20.1143                    5.5512      3.62x

        -- 82-87% of the predicted ceiling, and a LOSS at m=40, which the
        ceiling said too.

        ⚠ THOSE ARE TRAVERSAL FIGURES AND THE END-TO-END SOLVE GAINS LESS.
        A `solve` also does its setup, the replay that builds the waveform
        and the DFT, none of which this touches, and matrix-free spends an
        extra Newton iteration (below).  Measured end to end, same circuits:

              m    dense (iters)      matrix-free (iters)     speedup
            242      2.113 s (2)            1.557 s (2)        1.36x
            502      9.255 s (2)            6.131 s (3)        1.51x
           1002     52.402 s (2)           24.636 s (3)        2.13x

        Quote whichever answers the question being asked, and say which it
        is; 2.23x and 1.51x at m=502 are both true and are not the same
        measurement.

        ⚠ THE CONVERGENCE TEST IS NOT BIT-IDENTICAL TO `analysis.fsolve`'s,
        and it cannot be.  `fsolve` scales its residual test by
        `|J| . |x|`, an ELEMENTWISE absolute value of the Jacobian, which no
        matrix-free method has.  The substitute here is `|x| + |M x| + |F|`,
        one extra matvec per iteration.

        ⚠ AND IT IS NOT PROVABLY THE STRICT DIRECTION.  This docstring first
        claimed the substitute was a LOWER bound on `fsolve`'s scale, so
        that the test could only ever be stricter.  That is FALSE: at
        `M = I` the true scale `|I - M| . |x|` is zero while the substitute
        is `2|x|`, so the substitute is the LARGER one there, and at `M = 0`
        they are equal.  Neither dominates the other in general.

        What is measured, on the RC ladder at m=242/502/1002: the two paths
        agree on the converged waveform to 1.1e-16 and on the converged/not
        verdict, and matrix-free takes ONE MORE Newton iteration at m>=502
        (3 against 2) -- so it is stricter in practice here, and still wins
        on wall time while doing 50% more traversals.  One circuit is not a
        proof of direction, and this is the first thing to check if the two
        paths ever disagree on convergence.
        """
        import scipy.sparse.linalg as spla
        m = self.cir.n - 1
        z = np.asarray(z0, dtype=float).copy()
        ier, mesg, xdiff = 2, 'No convergence', None
        for _i in range(maxiter):
            C0, steps, x_last, x_prev = self._traverse_factored(
                z[:m], z[m:], times, hs)
            F = np.concatenate((z[:m] - np.asarray(x_last),
                                z[m:] - np.asarray(x_prev)))

            def _mv(v, _C0=C0, _st=steps):
                return v - self._monodromy_matvec(_C0, _st, v)

            J = spla.LinearOperator((2 * m, 2 * m), matvec=_mv, dtype=float)
            xdiff, _info = spla.gmres(
                J, -F, rtol=self.KRYLOV_TOLERANCE_FACTOR * reltol,
                restart=min(2 * m, 200), maxiter=min(2 * m, 400))
            z_new = z + xdiff

            ## `|J| . |x|` is not available without the matrix; see the
            ## docstring for why this substitute is the strict direction.
            I_scale = (np.abs(z_new)
                       + np.abs(self._monodromy_matvec(C0, steps, z_new))
                       + np.abs(F))
            conv_x = np.all(np.abs(xdiff)
                            < reltol * np.maximum(np.abs(z_new), np.abs(z))
                            + xtol)
            conv_f = np.all(np.abs(F) < reltol * I_scale + abstol)
            z = z_new
            if conv_x and conv_f:
                ier, mesg = 1, 'Success'
                break
        return z, {}, ier, mesg

    ## Below this, a multiplier says the mode decays by six decades in one
    ## period and no stability question turns on it -- so parasitic
    ## contamination at that level is not worth a warning.  Used only to
    ## keep the warning off circuits where the WHOLE spectrum is numerical
    ## noise; it never changes a reported number.
    SPECTRAL_NOISE_FLOOR = 1e-6

    def _spectral_report(self, M):
        """Split a composed spectrum into physical multipliers and parasitics.

        RECORDED SCOPE ITEM 3.  A k-step method turns an m-dimensional
        system into a k*m-dimensional discrete one, so the composed
        monodromy's spectrum carries `(k-1) m` PARASITIC roots beside the
        physical Floquet multipliers.  `max |eig|` over that mixture is only
        a stability verdict while the parasitic roots stay small -- which
        for Gear-2 they emphatically do (`(1/3)^N`, ~1e-95 at 200 points)
        and for a method whose spurious root sits nearer the unit circle
        they would not.  This separates them instead of hoping.

        THE DISCRIMINATOR IS THE EIGENVECTOR'S BLOCK STRUCTURE, not the
        eigenvalue.  The composed map acts on the PAIR `(x_0, x_{-1})`:

          - a PHYSICAL mode follows the linearised ODE, so its two halves
            are one timestep apart on a smooth trajectory and
            `v_{-1} = e^{-lambda h} v_0 -> v_0` as `h -> 0`;
          - a PARASITIC mode is `r^n u` for the method's spurious root `r`,
            so `v_{-1} = u / r` -- three times `v_0` for BDF-2, minus it for
            a trapezoidal-like root -- and the halves differ by O(1)
            whatever `h` is.

        So `||v_{-1} - v_0||` (against a unit-norm eigenvector) is O(h) for
        a physical mode and O(1) for a parasitic one.  MEASURED, and it is
        the h-scaling that makes it a prediction rather than a story: on the
        phase circuit the physical ratio falls 0.1281 -> 0.0316 when the
        grid goes from 50 to 200 points -- a factor of 4.05 for a factor of
        4 in `h` -- while the parasitic ratios sit at 1.0 to 10.  On the
        Q=20 RLC the parasitic ratio is 1.9997 against the 2.0 that BDF-2's
        `v_{-1} = 3 v_0` predicts exactly.

        ⚠ THE SPLIT IS BY RANK, NOT BY A THRESHOLD, and that was measured
        into the design rather than chosen.  A threshold of 0.25 was tried
        first and returned NO physical modes at all on a stiff RC ladder --
        `lambda h ~ 40` there, so every mode's halves differ by O(1) and the
        classifier called the entire spectrum parasitic, handing back a
        `spectral_radius` of `None` where the old code said 6e-15.  A
        `k`-step method on `m` states has EXACTLY `m` physical multipliers
        and `(k-1) m` spurious ones -- that is structural -- so the `m`
        smallest splits are the physical set by construction, and the
        question of where to put a cut never arises.

        ⚠ ON A STIFF CIRCUIT THE LABELS MAY STILL BE WRONG, and it does not
        matter: when the physical modes are themselves stiff, a parasitic
        root can have the smaller split and swap places with one.  Every
        mode involved then has `|mu|` at the noise floor, so the RADIUS is
        unaffected -- it is the labels, not the number, that degrade.  What
        this buys is the case that motivated the item: a method whose
        spurious root sits NEAR THE UNIT CIRCLE, where the physical modes
        are well resolved, the splits separate cleanly, and taking a
        maximum over the mixture would report the discretisation's own
        artefact as the orbit's stability.

        Returns `(rho, physical, parasitic)`: the spectral radius over the
        PHYSICAL multipliers only, and both sets sorted by magnitude.  An
        `m x m` monodromy (any one-step method, the plain path) has no pairs
        and no parasitic roots, so everything in it is physical.  `None`
        gives `(None, None, None)` -- the matrix-free path forms no
        monodromy at all.
        """
        if M is None:
            return None, None, None
        M = np.asarray(M)
        m = self.cir.n - 1
        try:
            ev, V = np.linalg.eig(M)
        except np.linalg.LinAlgError:                     # pragma: no cover
            return None, None, None

        if M.shape[0] != 2 * m:
            ## one-step method: the monodromy IS the physical map
            phys = np.sort(np.abs(ev))[::-1]
            return float(phys[0]), phys, np.array([])

        ## columns of `V` are unit-norm, so this needs no denominator and
        ## cannot divide by a vanishing block -- a mode living entirely in
        ## one half reads as O(1) here, which is what it is.
        split = np.linalg.norm(V[m:, :] - V[:m, :], axis=0)
        order = np.argsort(split)
        phys = np.sort(np.abs(ev[order[:m]]))[::-1]
        para = np.sort(np.abs(ev[order[m:]]))[::-1]
        rho = float(phys[0])
        ## ⚠ THE POINT AT WHICH THIS STOPS BEING BOOKKEEPING.  While the
        ## parasitic roots are 80 orders down, separating them changes no
        ## number and only documents why the maximum was safe.  Once one
        ## climbs to within a decade of the physical spectrum, the method's
        ## spurious roots are a real part of what the analysis reports and
        ## the user is entitled to know before reading a stability verdict.
        if len(para) and para[0] > 0.1 * rho and rho > self.SPECTRAL_NOISE_FLOOR:
            warnings.warn(
                'PSS: this method\'s PARASITIC roots are no longer '
                'negligible -- the largest is %.4g against a physical '
                'spectral radius of %.4g. A k-step method contributes '
                '(k-1)*m spurious roots to the composed monodromy, and '
                '`spectral_radius` now reports the maximum over the '
                'PHYSICAL multipliers only (separated by eigenvector block '
                'structure). Treat the separation as load-bearing here '
                'rather than cosmetic: check `floquet_multipliers` and '
                '`parasitic_roots` before drawing a stability conclusion.'
                % (para[0], rho), RuntimeWarning, stacklevel=3)
        return rho, phys, para

    def _install_history(self, x0_in, xm1_in, dt, h_prev=None):
        """Open a run ON a solved two-point history rather than a seed.

        `_begin_run(x_{-1})` opens the rings on the earlier point and the
        push puts `x_0` in front of it, so the first real step reads
        `q(x_0)` and `q(x_{-1})` -- two genuine solved points.  The flags
        then say what is true of them: a step of `dt` has been taken, and
        the run is no longer opening, so nothing drops order.

        `_dt_last2` stays None on purpose.  The THIRD charge in the ring is
        still `q(x_{-1})` repeated, and the LTE estimator differences three,
        so its opening reading remains unsound and the report goes on
        discarding it -- a solved history fixes what the SOLUTION reads,
        not what the estimator does.

        Shared by `_traverse_solved_history` and the final replay, because a
        replay that opened differently from the solve would report a
        waveform the residual was never driven to zero on.
        """
        tr = self._transient()
        _alphas, b = tr._get_integrator().companion_coefficients(dt, dt)
        tr._begin_run(self._insert_refnode(xm1_in), self.cir.n)
        tr._dt = dt

        ## ⚠ A `b != 0` COMPANION IS REFUSED HERE, AND THE FIRST REASON
        ## GIVEN FOR IT WAS WRONG.  It said a solved history carries CHARGES
        ## while such a method also reads `iq_{-1}`, "which no charge
        ## determines".  The DAE determines it exactly -- a converged point
        ## satisfies `i(x) + iq + u = 0`, so `iq_{-1} = -(i(x_{-1}) + u)`,
        ## the same identity item 4d rests on -- and seeding it was tried.
        ##
        ## It fails for the derivative running the OTHER way.  A one-step
        ## companion reads only `iq_{-1}`, so the trajectory depends on
        ## `x_{-1}` solely through it, and `d(iq_{-1})/d x_{-1} = -G` is
        ## SINGULAR wherever a node carries no conductance -- every purely
        ## reactive node, which is most of a resonator.  Admitting `x_{-1}`
        ## as m unknowns then leaves the 2m x 2m system rank-deficient:
        ## measured, `LinAlgError: Singular matrix` on 25 tests at once.
        ##
        ## The right second unknown for such a method is `iq_{-1}` ITSELF --
        ## the `(x, iq)` state its monodromy already uses -- closed by
        ## `iq_{-1} = iq_{N-1}`.  That is a different formulation, not a
        ## seeding fix, and it is not built.
        if b:
            raise NotImplementedError(
                'a solved entering history admits `x_{-1}` as the second '
                'unknown, and a companion with a b != 0 term depends on it '
                'only through `iq_{-1} = -(i(x_{-1}) + u)`, whose derivative '
                '-G is singular at every purely reactive node -- so the '
                'enlarged system would be rank-deficient. Such a method '
                'needs `iq_{-1}` itself as the unknown, which is a different '
                'formulation; this refuses rather than solving a singular '
                'one.')
        ## The charge half of `_push_history`, without its `_iq` roll: no
        ## step has been solved yet, so there is no companion current to
        ## push -- and a `b = 0` companion never reads one, which the guard
        ## above is what makes true.
        q0 = tr.cir.q(self._insert_refnode(x0_in), tr.epar)
        tr._qlast = self.toolkit.concatenate(
            (self.toolkit.array([q0]), tr._qlast))[:-1]
        tr._q_cache = None
        tr._is_first_step = False
        tr._no_history = False
        ## ⚠ THE STEP THAT PRODUCED `x_0` IS THE PERIOD'S LAST ONE, NOT ITS
        ## FIRST.  `x_{-1}` sits one step BEFORE `x_0`, and on a periodic
        ## grid that step is `hs[-1]`.  With a uniform grid the two are
        ## equal and this never showed; on a caller's grid (item 5) with a
        ## 16438:1 spread, handing `hs[0]` to a method that reads `h_last`
        ## states a step ratio that never happened.
        tr._dt_last = dt if h_prev is None else h_prev
        tr._dt_last2 = None
        self._history_is_solved = True
        return tr

    def _C_at(self, x_reduced):
        """The reduced capacitance at a point, without taking a step."""
        tr = self._transient()
        C = tr.cir.C(self._insert_refnode(x_reduced), tr.epar)
        (C,) = remove_row_col((C,), self.irefnode, self.toolkit)
        return C

    def _traverse(self, x_in, T, times, hs, want_dT):
        """One pass over the period, with the sensitivities accumulated.

        Returns ``(x0, x_end, dx_end/dx0, dx_end/dT)`` -- the last only when
        asked for.  Shared by the fixed-period and autonomous systems so the
        period map is written once; they differ only in what they build from
        it.

        EVERY SHOOTING ITERATION IS ITS OWN RUN.  phi must be a function of
        its arguments alone; if iteration k+1 inherited the ring buffers
        iteration k ended with, the period map would depend on which
        iteration it was and the monodromy would be the derivative of
        something else.  `_begin_run` also makes the first step
        `is_first_step`, so a multi-step method opens at order 1 -- the same
        restart the old `iq_last=None` produced, now for the integrator's
        own reason.
        """
        toolkit = self.toolkit
        n = self.cir.n
        self._want_dfdh = want_dT
        self._begin_period(x_in)
        x = self.solve_timestep(x_in, times[0], hs[0])
        iq_last = self._iq
        x0 = copy(x)

        ## `Px[k]` is d(x_{n-k})/d(x0), `Pq` is d(iq_n)/d(x0), `Cs[k]` the
        ## capacitance of step n-k.  Two of each -- as far back as any
        ## method here reaches.  Both rings open seeded with the entering
        ## step, mirroring how the transient seeds `_qlast` with `q0`
        ## repeated: at the start of a period there is no earlier point to
        ## differentiate against.
        eye = np.asarray(toolkit.eye(n - 1))
        Px = [eye, eye]
        Cs = [copy(np.asarray(self._C)), copy(np.asarray(self._C))]
        a_first, b_first = self._coeffs
        Pq = (a_first[0] * Cs[0] if b_first else np.zeros((n - 1, n - 1)))
        ## The period column, propagated the same way with one extra source
        ## term.  Zero at the start: the entering state does not depend on T.
        Pt = [np.zeros(n - 1), np.zeros(n - 1)]
        Pqt = np.zeros(n - 1)

        ## Kept for PAC.
        self.Cvec = [copy(self._C)]
        self.Jtvec = [copy(self._Jf)]
        self.times = times

        for _j, t in enumerate(times[1:]):
            dt = hs[min(_j, len(hs) - 1)]
            x = copy(self.solve_timestep(x, t, dt, iq_last=iq_last))
            iq_last = self._iq
            self.Cvec.append(copy(self._C))
            self.Jtvec.append(copy(self._Jf))

            ## ONE RECURSION FOR EVERY METHOD.  Each writes its companion as
            ## `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`, so differentiating
            ## the step gives
            ##
            ##     S    = sum_{k>=1} a_k C_{n-k} P_{n-k} + b Pq
            ##     P_n  = -Jf_n^-1 S
            ##     Pq_n = a_0 C_n P_n + S
            ##
            ## Euler is `b = 0` reaching back one step, trapezoidal `b = -1`
            ## reaching back one, Gear-2 `b = 0` reaching back two.  The
            ## coefficients come from the integrator that RAN, so an
            ## order-dropped step contributes its own.
            ##
            ## ⚠ A SOLVE, NOT AN INVERSE (stage 11).  `inv(Jf) @ ...` formed
            ## a dense inverse per timestep per iteration and squared the
            ## condition number it then multiplied through.
            alphas, b = self._coeffs
            Jf = np.asarray(self._Jf)
            C_new = np.asarray(self._C)

            S = b * Pq if b else np.zeros_like(Px[0])
            for k in range(1, len(alphas)):
                S = S + alphas[k] * (Cs[k - 1] @ Px[k - 1])
            Px_new = -self.toolkit.linearsolver(Jf, S)
            Pq = alphas[0] * (C_new @ Px_new) + S

            if want_dT:
                ## The SAME recursion, plus the step size's own dependence
                ## on T.  Every step scales together (`h = T/(N-1)`), so
                ## `dh/dT = h/T`, and `df/dh` at fixed solution is Fang's
                ## `p` -- already shared.  For an autonomous circuit its
                ## `du/dt` half vanishes, which is what makes this term the
                ## companion derivative alone.
                St = b * Pqt if b else np.zeros_like(Pt[0])
                for k in range(1, len(alphas)):
                    St = St + alphas[k] * (Cs[k - 1] @ Pt[k - 1])
                St = St + np.asarray(self._dfdh).ravel() * (dt / T)
                Pt_new = -self.toolkit.linearsolver(Jf, St)
                Pqt = alphas[0] * (C_new @ Pt_new) + St
                Pt = [Pt_new, Pt[0]]

            Px = [Px_new, Px[0]]
            Cs = [copy(C_new), Cs[0]]

        self._want_dfdh = False
        ## Kept for the autonomous check after the solve: its spectrum is
        ## the only place a free period announces itself.
        self._monodromy = Px[0]
        return x0, x, Px[0], (Pt[0] if want_dT else None)

    def _transient(self):
        """The `Transient` this analysis integrates with.

        PSS used to carry its OWN transcription of one integrator step --
        the third in the tree, after `Transient` and `JAXTransient` -- and
        it had already cost the two defects its docstring records: `method`
        declared and never read, and a companion current fed back from the
        iterate before the converged one.  Driving the real thing removes
        the copy and brings what came with it: the limiting machinery, PCNR,
        breakpoint order drops, and the continuation rescue.

        The tolerances are handed over unchanged, which is the point of
        `newton_tolerance_vectors`: `reltol`/`iabstol`/`vabstol` mean the
        same thing on both sides, so passing them through is a no-op in
        meaning.
        """
        if getattr(self, '_tran', None) is None:
            from pycircuit.circuit.transient import Transient
            from pycircuit.circuit.integrator import (EulerIntegrator,
                                                      TrapezoidalIntegrator,
                                                      Gear2Integrator)
            ## ⚠ A MAPPING, not an if/else on 'euler'.  Written as
            ## `EulerIntegrator() if method == 'euler' else Trapezoidal...`
            ## it silently ran trapezoidal for every other name -- caught
            ## while adding 'gear', which produced numbers identical to
            ## trap's to the last digit.  This class has already paid once
            ## for a `method` that selected nothing; a dict raises KeyError
            ## on a name nobody wired.
            integ = self._integrator_for(self.par.method)
            self._tran = Transient(
                self.cir, toolkit=self.toolkit, integrator=integ,
                reltol=self.par.reltol, iabstol=self.par.iabstol,
                vabstol=self.par.vabstol, maxiter=self.par.maxiter,
                analysis=self.par.analysis,
                lte_vabstol=self.par.lte_vabstol,
                lte_iabstol=self.par.lte_iabstol,
                TRTOL=self.par.TRTOL, relref=self.par.relref)
            self._tran.irefnode = self.irefnode
        return self._tran

    def _begin_period(self, x_reduced):
        """Start one traversal of the period from a clean integrator state.

        Every shooting iteration re-integrates the SAME interval from its own
        `x0`, so "begin a run" happens once per iteration here, not once per
        analysis.  Without the reset, iteration k+1 would inherit the ring
        buffers iteration k ended with and the period map would depend on
        which iteration it was -- phi must be a function of `x0` alone or the
        monodromy is the derivative of something else.
        """
        tr = self._transient()
        tr._begin_run(self._insert_refnode(x_reduced), self.cir.n)
        return tr

    def _insert_refnode(self, x):
        return self.toolkit.concatenate(
            (x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))

    def solve_timestep(self, x0, t, dt, refnode=gnd, iq_last=None):
        """One timestep of the inner transient, taken by `Transient`.

        This used to be a private transcription of one integrator step --
        the third in the tree -- and it had already cost two defects that
        its own comments recorded: `method` was declared and never read, so
        PSS was backward-Euler only, and the companion current fed back to
        the next step belonged to the iterate BEFORE the converged one.
        Both are structurally impossible now: the integrator is an
        `Integrator` object driven by `Transient.get_diff`, and the
        companion current is the one that class stores at its own converged
        point.

        What came with the change, none of which the copy had: the limiting
        machinery (measured -- a rectifier whose diode never turned on, so
        the non-conducting solution was returned as a converged periodic
        steady state), PCNR when the circuit and Parameters ask for it,
        breakpoint order drops, and the continuation rescue.

        `dt` is imposed by the caller: PSS owns the grid, which is what
        keeps the period map a smooth function of `x0`.  `iq_last` is
        retained in the signature for callers that pass it, but the
        companion history now lives in the `Transient`'s own ring buffers,
        rolled here through `_push_history`.

        The measured cost of backward Euler on a limit cycle is unchanged
        and still the reason `method` matters -- it damps exactly what PSS
        exists to find:

            steps/period    PSS peak    fraction of analytic
                      20      2.63 V       13.2%
                      50      5.61 V       28.1%
                     100      8.81 V       44.1%
                     200     12.20 V       61.0%
        """
        toolkit = self.toolkit
        irefnode = self.irefnode
        tr = self._transient()

        ## ONE INTEGRATOR STEP, taken by the class that owns the definition.
        ## `Transient.solve_timestep` applies the chosen integrator through
        ## `get_diff` (so `method` selects something because the integrator
        ## object does), the limiting machinery, PCNR when asked for, and the
        ## continuation rescue.  None of that existed on the copy this
        ## replaced.
        tr._dt = dt
        x_full, J_full = self._transient_step(tr, x0, t)

        ## d(residual)/dh at the converged point, for the period
        ## derivative -- BEFORE `_push_history`, because it reads `_qlast`
        ## as the PREVIOUS charges, which the push is about to overwrite.
        ## `Transient.residual_dh` is Fang's `p`, already shared: it is
        ## `d(iq)/dh + du/dt`, and for an AUTONOMOUS circuit the second term
        ## is identically zero -- which is exactly why solving for the
        ## period is tractable here and would not be on a driven circuit,
        ## where scaling T also moves every source evaluation.
        if self._want_dfdh:
            (self._dfdh,) = remove_row_col(
                (tr.residual_dh(x_full, t, dt),), irefnode, toolkit)
        ## Measured, not controlled: the grid is the caller's, so nothing can
        ## act on this.  Also before the push, for the same reason.
        if self._want_lte:
            self._lte = tr.step_lte(x_full, self._insert_refnode(x0), J_full)
            ## A SEAM STEP IS ONE WHOSE COMPANION READS THE ENTERING
            ## UNKNOWN, not merely one whose ESTIMATOR does.
            ##
            ## ⚠ THIS CONDITION WAS `h_last2 is None` ALONE, AND THAT FLAGGED
            ## A PHANTOM FOR TWO METHODS OF THREE.  `h_last2 is None` is the
            ## transient's statement that the third past charge is not real,
            ## which is the reach of the LTE estimator's third divided
            ## difference -- not the reach of the integrator.  Euler's
            ## companion reads `q_{n-1}`; trapezoidal's reads `q_{n-1}` and
            ## `iq_{n-1}`, and the order-dropped opening step supplies an
            ## `iq` consistent with it, which is what that drop is FOR.
            ## Neither can see the fabricated charge at all.  Measured
            ## (`benchmarks/pss_seam_cost.py`): trapezoidal's seam reading
            ## was 15.1 times tolerance while its cost is 1.3e-11 V, and
            ## euler's 0.286 against 5.1e-12 V.  Both are exactly zero; the
            ## reading was an artefact of the measurement.
            ##
            ## Gear-2 reads `q_{n-2}` -- which at that step IS the entering
            ## unknown -- and the shooting condition constrains `x(0)` to
            ## equal `x(P)`, NOT `x_in` to be the orbit's own `x(-dt)`.  So
            ## `x_in` sits O(h^2) off a real history point and Gear-2 reads
            ## it as one.  That one costs 1.266e-01 V at 100 points/period
            ## against an interior contribution of 1.070e-01 -- the seam is
            ## 54% of its total error -- and it is the term that STOPS
            ## converging: it falls as h^2 while the interior falls faster,
            ## so its share grows to 68% at 200 points and 73% at 400.
            ##
            ## TWO DIFFERENT THINGS ARE TRUE OF AN OPENING STEP, and the
            ## first version of this conflated them -- which showed up as
            ## trapezoidal's phantom simply MOVING from the seam into the
            ## interior total (0.340 -> 15.47) when the seam test was
            ## tightened.  Suppressing a bad number is not the same as
            ## classifying it.
            ##
            ##   the ESTIMATE is invalid when the ESTIMATOR differences a
            ##   charge that was never a real point.  Its divided difference
            ##   reaches `p = ORDER + 1` charges back -- 2 for Euler, 3 for
            ##   the second-order pair -- so at the step with only two real
            ##   past charges, euler's estimate is sound and trap's and
            ##   gear2's are not.  An unsound estimate is DISCARDED: it is
            ##   not an interior reading and, on its own, not a seam either.
            ##
            ##   a SEAM exists when the COMPANION reads the entering
            ##   unknown, i.e. when its charge reach `len(alphas) - 1` is
            ##   deep enough to touch it.  Only Gear-2's is.
            ##
            ## `_dt_last2 is None` says exactly "two real past charges" here
            ## (it is set from `_dt_last` one step later), so it is the step
            ## index in disguise; both tests are written against the count.
            _real_past = 2 if tr._dt_last2 is None else 3
            _p = getattr(tr.active_integrator, 'ORDER', 1) + 1
            _reach = len(tr._companion_coeffs[0]) - 1
            self._lte_valid = _real_past >= _p
            ## `_history_is_solved` is that formulation saying the
            ## deepest charge is an UNKNOWN the solve closed, not a stand-in
            ## -- so there is no seam to report even though the companion
            ## reaches that far.  Without this the fix would go on flagging
            ## the defect it removed.
            self._lte_seam = (_reach >= _real_past
                              and not self._history_is_solved)

        ## The history advance is the accept path's, called rather than
        ## copied -- and `_dt_last` must roll AFTER the step, because
        ## `get_diff` read it as `h_last` while solving.
        tr._push_history(x_full)
        tr._dt_last2 = tr._dt_last
        tr._dt_last = dt
        tr._is_first_step = False
        tr._no_history = False

        ## Reduced-system views for the shooting Jacobian.  `_Geq` is the
        ## companion conductance the step actually used, which is the factor
        ## the monodromy needs; `_iq` is kept for the caller's own bookkeeping
        ## as before.
        (self._Jf, self._Geq, self._C) = remove_row_col(
            (J_full, tr._Geq, tr._Cmat), irefnode, toolkit)
        ## The coefficients of the integrator that ACTUALLY ran this step --
        ## an order drop on the opening step reports Euler's, which is what
        ## the propagation must use for that step.
        self._coeffs = tr._companion_coeffs
        self._iq = tr._iq

        x = toolkit.concatenate((x_full[:irefnode], x_full[irefnode + 1:]))
        return x

    def _transient_step(self, tr, x0_reduced, t):
        """`Transient.solve_timestep` on the FULL vector, returning
        ``(x, J)``.  PSS works on the reduced system throughout; this is the
        one place the two conventions meet."""
        x_full = self._insert_refnode(x0_reduced)
        x, _feval, J, _f = tr.solve_timestep(x_full, t)
        return x, J


    def solve(self, refnode=gnd, period=1e-3, x0=None, timestep=1e-6,
              maxiterations=20, grid=None, matrix_free=False):
        """Solve for the periodic steady state.

        `grid` is RECORDED SCOPE ITEM 5: a sequence of step FRACTIONS of the
        period, summing to 1, used in place of the uniform `timestep` grid.
        Fractions rather than absolute times because an autonomous period is
        an unknown and every step has to scale with it.  See
        `_period_grid`; the grid is still frozen for the whole solve, so the
        shooting Newton stays exact.

        `matrix_free` is RECORDED SCOPE ITEM 6: solve the outer system
        without ever forming the monodromy, propagating ONE vector per
        Krylov iteration instead of `2m` columns per step.  It is worth
        asking for on LARGE circuits only -- measured on the RC ladder it
        LOSES at m=40 (0.82x) and wins above roughly m=250: 1.40x at m=242,
        2.23x at m=502, 3.62x at m=1002.  ⚠ It buys that with memory, `2 N
        m^2` doubles of stored factorisations (~800 MB at m=1002, 50
        points), and it does NOT produce a monodromy, so `spectral_radius`
        is `None` after a matrix-free solve -- not forming that matrix is
        the entire point.  Supported for the DRIVEN solved-history path
        only; anything else raises, rather than silently taking the dense
        route the caller asked to avoid.
        """
        self.period = period
        toolkit = self.toolkit

        irefnode=self.cir.get_node_index(refnode)
        n = self.cir.n
        dt = timestep
        if x0 is None:
            x = toolkit.zeros(n-1) #currently without reference node !
        else:
            x = x0 # reference node not included !

        #create vector with timepoints and a more fitting dt
        times, hs = self._period_grid(period, int(period / dt), grid)
        npts = len(times)
        self._grid_fracs = (None if grid is None
                            else np.asarray(grid, dtype=float))
        alpha = 1

        ## AUTONOMY IS DECIDED BEFORE THE SOLVE, because it decides which
        ## system is solved.  Structural and exact -- see `_is_autonomous`.
        self.autonomous = self._is_autonomous(times)
        phase_k, phase_pin = 0, 0.0
        if self.autonomous:
            ## An unseeded autonomous run starts at the origin, which IS a
            ## periodic solution -- the trivial one -- and the free-period
            ## system would sit there just as contentedly as the fixed one
            ## did.  The operating point is the honest default: for a phase
            ## accumulator `ic` pins it on the orbit.
            if x0 is None:
                from pycircuit.circuit.dcanalysis import DC
                xdc = np.asarray(DC(self.cir, toolkit=self.toolkit).solve().x,
                                 dtype=float).reshape(-1)
                x = np.concatenate((xdc[:irefnode], xdc[irefnode + 1:]))

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
            phase_pin = float(np.asarray(_x1)[phase_k])

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
        if method not in ('euler', 'trap', 'trapezoidal', 'gear', 'gear2'):
            raise ValueError(
                "method must be 'euler', 'trap' or 'gear', not %r" % (method,))

        ## Whether the entering history joins the unknowns.  Decided once,
        ## here, because it chooses which system is solved -- like autonomy,
        ## and after it, since the two are not composed yet.
        solved_history = self._solves_history()
        self.solved_history = solved_history
        xm1_ss = None

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
        newton = True

        def func(x):
            x0, x_end, Mx, _Mt = self._traverse(x, period, times, hs,
                                                want_dT=False)
            D = np.asarray(toolkit.eye(n - 1))
            return x0 - x_end, D - alpha * Mx

        def func_solved_history(z):
            """Residual and Jacobian when the entering history is an unknown.

            Unknowns are `(x_0, x_{-1})`; the equations are that BOTH close,

                F = [ x_0 - x_{N-1} ,  x_{-1} - x_{N-2} ]
                J = [[ I - A(N-1,0) , -A(N-1,-1) ],
                     [   -A(N-2,0)  , I - A(N-2,-1) ]]

            with `A(j,k) = d x_j / d x_k`, which is what `_traverse_solved_history`
            returns as one `n x 2n` block per row.  A two-step companion needs
            two states to be continued, so periodicity of ONE of them is an
            under-determined statement about the orbit -- which is the defect
            this replaces, measured at 1.266e-01 V for Gear-2 at 100 points
            per period.
            """
            m = n - 1
            x0_in, xm1_in = z[:m], z[m:]
            x_last, x_prev, P_last, P_prev = self._traverse_solved_history(
                x0_in, xm1_in, times, hs)

            D = np.asarray(toolkit.eye(m))
            J = np.zeros((2 * m, 2 * m))
            J[:m, :m] = D - alpha * P_last[:, :m]
            J[:m, m:] = -alpha * P_last[:, m:]
            J[m:, :m] = -alpha * P_prev[:, :m]
            J[m:, m:] = D - alpha * P_prev[:, m:]
            F = np.concatenate((np.asarray(x0_in) - np.asarray(x_last),
                                np.asarray(xm1_in) - np.asarray(x_prev)))
            return F, J

        def func_autonomous(z):
            """Residual and Jacobian of the FREE-PERIOD system.

            Unknowns are `(x0, T)`; equations are the period map's fixed
            point plus a phase condition, because without one the system is
            singular by construction -- every point on the orbit is a
            solution, so `I - M` has a null direction along it.

                F = [ x0 - phi_T(x0) ,  x0[k] - pinned ]
                J = [[ I - M , -dphi/dT ],
                     [ e_k^T ,     0    ]]

            `k` is the coordinate moving fastest at the seed, so the orbit
            crosses the pinning hyperplane transversally; pinning a slow
            coordinate makes the last row nearly parallel to the null
            direction it is there to remove.
            """
            x_in, T = z[:-1], float(z[-1])
            ## Rebuilt at the CURRENT T, which is what keeps `dh/dT = h/T`
            ## true of every step -- uniform or not.
            tms, hs_T = self._period_grid(T, npts, self._grid_fracs)
            x0, x_end, Mx, Mt = self._traverse(x_in, T, tms, hs_T,
                                               want_dT=True)

            D = np.asarray(toolkit.eye(n - 1))
            m = n - 1
            J = np.zeros((m + 1, m + 1))
            J[:m, :m] = D - alpha * Mx
            J[:m, m] = -np.asarray(Mt).ravel()
            J[m, phase_k] = 1.0
            F = np.zeros(m + 1)
            F[:m] = x0 - x_end
            F[m] = x0[phase_k] - phase_pin
            return F, J
        
        def func_autonomous_solved_history(z):
            """Both enlargements at once: `(x_0, x_{-1}, T)`.

            The free-period system and a solved entering history grow the
            SAME unknown vector in different directions, and a two-step
            method on an autonomous circuit needs both -- the period because
            it is not given, the history because the companion reads it.

                F = [ x_0  - x_{N-1} ,  x_{-1} - x_{N-2} ,  x_0[k] - pin ]

                J = [[ I - A(N-1,0) ,   -A(N-1,-1) , -dx_{N-1}/dT ],
                     [   -A(N-2,0)  , I - A(N-2,-1), -dx_{N-2}/dT ],
                     [     e_k^T    ,       0      ,       0      ]]

            ⚠ ONE PHASE ROW STILL SUFFICES, and it pins only the `x_0`
            block.  Time translation slides BOTH states along the orbit
            together -- the null vector is `(xdot(0) ds, xdot(-h) ds, dT)` --
            so the freedom stays one-dimensional and pinning `x_0[k]` kills
            it whenever `xdot(0)[k] != 0`.  That is the same reason `k` is
            chosen as the fastest-moving coordinate at the seed, and the
            same reason a slow one would leave the last row nearly parallel
            to the direction it exists to remove.

            ⚠ THE HISTORY POINT MOVES WITH T.  `x_{-1}` sits at `-T/(N-1)`,
            so its location tracks the unknown.  The residual is still the
            right statement -- `x_{-1}` and `x_{N-2}` are the same phase of
            the orbit at every T -- and the T column is just the propagation
            accumulated to step N-2, which the ring already carries.
            """
            m = n - 1
            x0_in, xm1_in, T = z[:m], z[m:2 * m], float(z[-1])
            tms, hs_T = self._period_grid(T, npts, self._grid_fracs)
            (x_last, x_prev, P_last, P_prev, Pt_last,
             Pt_prev) = self._traverse_solved_history(
                x0_in, xm1_in, tms, hs_T, T=T, want_dT=True)

            D = np.asarray(toolkit.eye(m))
            J = np.zeros((2 * m + 1, 2 * m + 1))
            J[:m, :m] = D - alpha * P_last[:, :m]
            J[:m, m:2 * m] = -alpha * P_last[:, m:]
            J[:m, 2 * m] = -np.asarray(Pt_last).ravel()
            J[m:2 * m, :m] = -alpha * P_prev[:, :m]
            J[m:2 * m, m:2 * m] = D - alpha * P_prev[:, m:]
            J[m:2 * m, 2 * m] = -np.asarray(Pt_prev).ravel()
            J[2 * m, phase_k] = 1.0

            F = np.zeros(2 * m + 1)
            F[:m] = np.asarray(x0_in) - np.asarray(x_last)
            F[m:2 * m] = np.asarray(xm1_in) - np.asarray(x_prev)
            F[2 * m] = np.asarray(x0_in)[phase_k] - phase_pin
            return F, J

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
        if matrix_free and not (solved_history and not self.autonomous):
            raise NotImplementedError(
                'PSS: matrix_free=True is built for the DRIVEN '
                'solved-history path only (a two-step method on a driven '
                'circuit), and this solve is %s. Item 6 replaces the '
                'sensitivity propagation, which is where a two-step '
                'method spends its 2m columns; the plain and autonomous '
                'systems have not been converted. Use matrix_free=False, '
                'or a method whose companion reaches two charges back '
                "(method='gear') on a driven circuit."
                % ('autonomous' if self.autonomous
                   else 'on the plain single-unknown path'))

        ## Find periodic steady state x-vector
        if self.autonomous and solved_history:
            ## BOTH unknowns and the period.  The floors follow the same
            ## rule as the plain autonomous system: the two state blocks
            ## take the solution-unit tolerance, and the row that adds a
            ## TIME as an unknown takes a time as its floor -- mixing them
            ## is flavour error F6(a) one row further out.
            m_ = n - 1
            xa = np.asarray(x, dtype=float)
            z0 = np.concatenate((xa, xa, [period]))
            abstol_z = np.concatenate((_tol, _tol, [_tol[phase_k]]))
            xtol_z = np.concatenate((_tol, _tol, [1e-15 * period]))
            z_ss, _info, _ier, _mesg = self._free_period_solve(
                func_autonomous_solved_history, z0, abstol_z, xtol_z,
                _shoot_reltol, maxiterations, period)
            x0_ss, xm1_ss = z_ss[:m_], z_ss[m_:2 * m_]
            self.period = period = float(z_ss[-1])
            times, hs = self._period_grid(period, npts, self._grid_fracs)
        elif self.autonomous:
            ## The period joins the unknowns.  Its residual row is the phase
            ## condition -- in the units of the coordinate it pins, hence
            ## `_tol[phase_k]` -- while the UNKNOWN it adds is a time, whose
            ## own floor has to be a time; mixing the two is the flavour
            ## error F6(a) names, one row further out.
            z0 = np.concatenate((np.asarray(x, dtype=float), [period]))
            abstol_z = np.concatenate((_tol, [_tol[phase_k]]))
            xtol_z = np.concatenate((_tol, [1e-15 * period]))
            z_ss, _info, _ier, _mesg = self._free_period_solve(
                func_autonomous, z0, abstol_z, xtol_z, _shoot_reltol,
                maxiterations, period)
            x0_ss = z_ss[:-1]
            self.period = period = float(z_ss[-1])
            ## The grid follows the solved period; everything downstream --
            ## the replay, the waveform, the DFT -- must use it or the
            ## answer is reported on a period the solver rejected.
            times, hs = self._period_grid(period, npts, self._grid_fracs)
        elif solved_history:
            ## THE SEED IS THE OLD FORMULATION'S ASSUMPTION, written down:
            ## `x_{-1} = x_0`.  It is what the plain path silently assumes
            ## (it seeds both charge rings with the entering state), so an
            ## solved-history run starts where a plain one starts and the
            ## comparison between them is about the SOLVE, not the seed.
            xa = np.asarray(x, dtype=float)
            z0 = np.concatenate((xa, xa))
            tol_z = np.concatenate((_tol, _tol))
            if matrix_free:
                ## ⚠ THE MONODROMY IS NOT FORMED, so it must not be
                ## REPORTED either.  `_monodromy` survives on the object
                ## from any earlier traversal, and `spectral_radius` reads
                ## it below without knowing which run wrote it -- so a
                ## matrix-free solve that left it alone would report the
                ## PREVIOUS solve's radius as this one's.  Cleared here,
                ## and `spectral_radius` is documented as None on this path.
                self._monodromy = None
                z_ss, _info, _ier, _mesg = self._matrix_free_solve(
                    z0, times, hs, tol_z, tol_z, _shoot_reltol,
                    maxiterations)
            else:
                z_ss, _info, _ier, _mesg = analysis.fsolve(
                    func_solved_history, z0, maxiter=maxiterations,
                    reltol=_shoot_reltol, abstol=tol_z, xtol=tol_z,
                    toolkit=self.toolkit, full_output=True)
            x0_ss, xm1_ss = z_ss[:n - 1], z_ss[n - 1:]
        else:
            x0_ss, _info, _ier, _mesg = analysis.fsolve(
                func, x, maxiter=maxiterations, reltol=_shoot_reltol,
                abstol=_tol, xtol=_tol, toolkit=self.toolkit, full_output=True)
        self.converged = (_ier == 1)
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
            warnings.warn(
                'PSS: the shooting solve did not converge in %d iterations '
                '(method=%r, %s). The returned waveform is the last iterate, '
                'not a periodic steady state -- raise maxiterations, or use '
                "method='euler', which solves a true Newton system while "
                'trapezoidal still uses successive substitution.'
                % (maxiterations, method,
                   'true Newton' if newton else 'successive substitution'),
                RuntimeWarning, stacklevel=2)
        
        ## THE THIRD LEVEL, MEASURED ON THE WAY OUT.
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
        if solved_history:
            self._install_history(x0_ss, xm1_ss, hs[0], h_prev=hs[-1])
            tr = self._transient()
            X = [np.asarray(x0_ss, dtype=float)]
            walk = times[1:]
        else:
            X = [x0_ss]
            tr = self._begin_period(x0_ss)
            walk = times
        ## Fresh probe, so `relref='sigglobal'`'s running signal maximum is
        ## the period's, not something an earlier shooting iteration saw.
        tr._lte_probe = None
        self._want_lte = True
        lte_seen = []
        for _j, t in enumerate(walk):
            x = self.solve_timestep(X[-1], t, hs[min(_j, len(hs) - 1)])
            if self._lte is not None:
                lte_seen.append((float(self._lte), float(t), self._lte_seam,
                                 self._lte_valid))
            X.append(copy(x))
        self._want_lte = False

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
            warnings.warn(
                'PSS: the shooting solve converged, but the periodic '
                'solution is not resolved at this accuracy (method=%r, %d '
                'points per period). Local truncation error reaches %.3g '
                'times tolerance %s: %s. Neither Newton criterion can see '
                'this -- they ask whether the discrete equations were '
                'solved, not whether the discretisation is the right one. '
                '(peak interior %s at t=%.6g s, period total %s, opening '
                'steps %s; relax lte_vabstol/lte_iabstol/TRTOL if this '
                'accuracy is intended.)'
                % (method, npts, v, where, fix,
                   'n/a' if self.max_lte is None else '%.3g' % self.max_lte,
                   -1.0 if self.max_lte_time is None else self.max_lte_time,
                   'n/a' if self.total_lte is None
                   else '%.3g' % self.total_lte,
                   'n/a' if self.max_lte_seam is None
                   else '%.3g' % self.max_lte_seam),
                RuntimeWarning, stacklevel=2)

        ## ⚠ THE PLAIN PATH'S FIRST ENTRY IS A SEED, THE OTHER PATH'S IS A
        ## SOLUTION.  Plain takes N steps from `x0_ss` and reports their
        ## results; the other starts AT `x_0` and takes N-1, so dropping the
        ## first would drop a real point and shift the waveform by a step.
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
            _pts = [np.asarray(v, dtype=float) for v in X]
            _d = np.array([float(np.max(np.abs(v - _pts[0]))) for v in _pts])
            _diam = float(np.max(_d))
            ## ⚠ THE THRESHOLD IS THE GRID, NOT A FIXED FRACTION.  A k-fold
            ## orbit returns BETWEEN grid points, so the nearest approach is
            ## bounded below by how far the solution moves in one step: at
            ## 200 points over two periods that is ~1.6% of the diameter,
            ## and a fixed 1% test therefore fired on nothing.  Comparing
            ## against the per-step displacement is scale-free and tightens
            ## automatically as the grid refines.
            _step = float(np.max([np.max(np.abs(_pts[i + 1] - _pts[i]))
                                  for i in range(len(_pts) - 1)]))
            _edge = max(2, len(_d) // 20)
            _inner = _d[_edge:-_edge]
            ## ⚠ THE EARLIEST RECURRENCE, NOT THE NEAREST.  A three-fold
            ## orbit passes close to `x_0` at both `T/3` and `2T/3`, and
            ## `argmin` picked whichever happened to be numerically nearer
            ## -- it reported `2T/3` as "the fundamental", which is wrong by
            ## a factor of two and would have sent the reader to a period
            ## that is itself a multiple.
            _near = [j for j in range(_edge, len(_d) - _edge)
                     if _d[j] < 3.0 * _step and _d[j] < 0.25 * _diam]
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
                if True:
                    ## `times[_j]`, not `period * _j/(N-1)`: with a caller's
                    ## grid the points are not evenly spaced.
                    self.fundamental_period = float(
                        times[min(_j, len(times) - 1)])
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
                        RuntimeWarning, stacklevel=2)

        X = toolkit.array(X if solved_history else X[1:]).T

        # Insert reference node voltage
        X = toolkit.concatenate((X[:irefnode], 
                                 toolkit.zeros((1,len(times))), 
                                 X[irefnode:]))

        tpss = analysis.CircuitResult(self.cir, x=X, xdot=None,
                                      sweep_values=times, sweep_label='time', 
                                      sweep_unit='s')

        freqs, FX = freq_analysis(X[:,:-1], times[:-1])
        
        fpss = analysis.CircuitResult(self.cir, x=FX, xdot=None,
                                      sweep_values=freqs, sweep_label='freq', 
                                      sweep_unit='Hz')
        
        return InternalResultDict({'tpss': tpss, 'fpss': fpss})

class PAC(Analysis):
    """Small-signal analysis over a time varying operating point"""

    parameters  = [Parameter(name='analysis', desc='Analysis name', 
                             default='ac')]

    def __init__(self, cir, toolkit=None, **kvargs):
        self.parameters = super(PAC, self).parameters + self.parameters            
        super(PAC, self).__init__(cir, toolkit=toolkit, **kvargs)
    
    def solve(self, pss, freqs, refnode=gnd, period=1e-3, x0=None, timestep=1e-6,
              maxiterations=20):
        raise NotImplementedError(
            "PAC is withdrawn as unimplemented (stage 11). It forms the whole "
            "(N*M)x(N*M) operator densely: at N=137 unknowns and M=1000 time "
            "points that is 419.5 GiB (279.7 GiB for L in complex128 plus 139.8 "
            "GiB for B), which is not a tuning problem but a consequence of the "
            "formulation. It has also never been validated -- its only test was "
            "@unittest.skip('Skip failing test'). A thin advertised feature is "
            "worse than an absent one, so it says so instead of allocating.\n\n"
            "The body below is kept, unreachable, because it is the starting "
            "point for a matrix-free rewrite (Telichevesky, Kundert & White, DAC "
            "1995): PAC needs only products with the operator, never the operator "
            "itself. Use PSS for periodic steady state in the meantime.")

        tk = self.toolkit
        analysis_name = self.par.analysis
        print('solve PAC analysis_name = ' + analysis_name)
        ## Create U vector which is the RHS evaluated at every time instant
        T = pss.period
        times = pss.times[:-1]
        hs = tk.diff(pss.times)

        N = self.cir.n - 1 ## ref node removed
        M = len(times)

        irefnode = self.cir.get_node_index(refnode)
        (u0,) = remove_row_col((self.cir.u(0, analysis_name),), irefnode, self.toolkit)

        ## Create LHS matrix using backward Euler discretization
        L = tk.zeros((N*M, N*M),dtype=tk.cdouble)
        ## 0.1c: no dtype, so B was float64 while L is complex -- it worked by
        ## promotion, not by intent.  Corrected even though unreachable, so the
        ## starting point for a rewrite is not itself wrong.
        B = tk.zeros(L.shape, dtype=tk.cdouble)
        for i, (t, h, J, C) in enumerate(zip(times, hs, pss.Jtvec, pss.Cvec)):
            L[i*N:(i+1)*N, i*N:(i+1)*N] = J
            if i > 0:
                L[i*N:(i+1)*N, (i-1)*N:i*N] = -C / h
        B[0:N,(M-1)*N:M*N] = -pss.Cvec[-1] / hs[0]

        outfreq = []
        outV = []
        for fs in freqs:
            phase_shift = tk.zeros(N * M, dtype=tk.cdouble)
            u = tk.zeros(N * M, dtype=tk.cdouble)
            for i,t in enumerate(times):
                phase_shift[i*N:(i+1)*N] = tk.exp(2j*tk.pi*fs*t)
                u[i*N:(i+1)*N] = u0
            
            u *= phase_shift
            
            alpha = tk.exp(-2j*tk.pi*fs*T)

            ## Solve discrete-time AC-voltage vector
            v = tk.linearsolver(L + alpha*B, -u)
            
            ## multiply v matrix by exp(-j*2*pi*fs) so the spectrum
            ## is evaluated at 2*pi*(fs + 1/T) instead of 2*pi/T
            ## this will also make v T-periodic
            v_shifted = (v / phase_shift)

            freqs, V = freq_analysis(v_shifted.reshape(M,N),
                                     times, axis=0)

            outfreq.extend((abs(freqs + fs)).tolist())
            outV.extend(V.tolist())
            
        ## Sort on frequency
        freqs, X = zip(*sorted(zip(outfreq, outV)))

        X = np.array(X)
        freqs = np.array(freqs)

        # Insert reference node voltage
        irefnode = self.cir.get_node_index(refnode)
        X = tk.concatenate((X[:,:irefnode], 
                            tk.zeros((len(freqs),1)), 
                            X[:,irefnode:]), axis=1)


        res = analysis.CircuitResult(self.cir, x = X.T, 
                                        xdot=None,
                                        sweep_values=freqs, 
                                        sweep_label='freq', 
                                        sweep_unit='Hz')

        
        return res

        
            
