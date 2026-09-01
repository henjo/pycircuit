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
         diagnostic still says what it said.  ⚠ But it now reads the FULL
         2m x 2m map, whose spectrum carries the parasitic roots of the
         two-step discretisation alongside the physical multipliers.  See
         the literature note: controlling those roots is the whole subject,
         and a future threshold on this number should know it is not reading
         a purely physical spectrum.

         ⚠ THE VALUE CAVEAT, UNCHANGED BY BUILDING IT.  Gear-2's own phase
         error is +332 ppm against trapezoidal's +83 ppm at the same grid,
         both second order, so on THIS circuit composed Gear-2 is still the
         worse choice and `method='trap'` -- the default -- remains the
         right answer.  What the composition buys is that a two-step method
         is CORRECT when it is the right tool, which is a stiff autonomous
         circuit: `doc/transient_review.md` sec. 4.6 measures trapezoidal
         ringing at `|e_n/e_{n-1}| ~ 0.9960` at `h*lambda = -1000` where
         Gear-2 sits at 0.0972.  Before this, such a circuit had no good
         option.

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

    RECORDED SCOPE, in order, neither of these planned work yet:

      5. LTE-CHOSEN grid -- pick the step sequence from an adaptive run and
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
        return self._companion_reach() >= 2

    def _step_sensitivity(self, Px, Cs, Pq, Jf, C_new):
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
        """
        alphas, b = self._coeffs
        S = b * Pq if b else np.zeros_like(Px[0])
        for k in range(1, len(alphas)):
            S = S + alphas[k] * (Cs[k - 1] @ Px[k - 1])
        Px_new = -self.toolkit.linearsolver(Jf, S)
        Pq_new = alphas[0] * (C_new @ Px_new) + S
        return Px_new, Pq_new

    def _traverse_solved_history(self, x0_in, xm1_in, times, dt,
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
        tr = self._install_history(x0_in, xm1_in, dt)

        ## `P_0 = [I 0]`, `P_{-1} = [0 I]` -- the two unknowns, exactly.  The
        ## plain path seeds BOTH rings with `I`, which is the flat-history
        ## assumption written into the Jacobian; here there is nothing to
        ## assume.
        eye = np.asarray(toolkit.eye(m))
        zero = np.zeros((m, m))
        Px = [np.hstack((eye, zero)), np.hstack((zero, eye))]
        Cs = [np.asarray(self._C_at(x0_in)), np.asarray(self._C_at(xm1_in))]
        ## Zero because `_install_history` has already refused any companion
        ## with a `b != 0` term: with `b = 0` the recursion never reads
        ## `d(iq_{n-1})/d(unknowns)`, so there is nothing to seed.  A future
        ## two-step method carrying `iq` would need it as a third unknown,
        ## which is what that refusal says.
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
        for t in times[1:]:
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
        ## in there; a threshold read off `spectral_radius` is reading a
        ## spectrum with spurious members and should know it.
        M = np.vstack((Px[0], Px[1]))
        self._monodromy = M if want_dT else Px[0][:, :m]
        if want_dT:
            return x, x_prev, Px[0], Px[1], Pt[0], Pt[1]
        return x, x_prev, Px[0], P_prev

    def _install_history(self, x0_in, xm1_in, dt):
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
        if b:
            raise NotImplementedError(
                'a solved entering history carries CHARGES, so a companion '
                'with a b != 0 term needs iq_{-1} as an unknown of its own. '
                'No such two-step method exists here (gear2 has b = 0); this '
                'refuses rather than seeding an iq nothing solved for.')
        tr._begin_run(self._insert_refnode(xm1_in), self.cir.n)
        tr._dt = dt
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
        tr._dt_last = dt
        tr._dt_last2 = None
        self._history_is_solved = True
        return tr

    def _C_at(self, x_reduced):
        """The reduced capacitance at a point, without taking a step."""
        tr = self._transient()
        C = tr.cir.C(self._insert_refnode(x_reduced), tr.epar)
        (C,) = remove_row_col((C,), self.irefnode, self.toolkit)
        return C

    def _traverse(self, x_in, T, times, dt, want_dT):
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
        x = self.solve_timestep(x_in, times[0], dt)
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

        for t in times[1:]:
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
              maxiterations=20):
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
        times,dt=toolkit.linspace(0,period,num=int(period/dt),endpoint=True,
                             retstep=True)
        npts = len(times)
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
            self._begin_period(x)
            _x1 = self.solve_timestep(x, times[0], dt)
            _x2 = self.solve_timestep(_x1, times[1], dt, iq_last=self._iq)
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
            x0, x_end, Mx, _Mt = self._traverse(x, period, times, dt,
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
                x0_in, xm1_in, times, dt)

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
            tms, h = toolkit.linspace(0.0, T, num=npts, endpoint=True,
                                      retstep=True)
            x0, x_end, Mx, Mt = self._traverse(x_in, T, tms, h, want_dT=True)

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
            tms, h = toolkit.linspace(0.0, T, num=npts, endpoint=True,
                                      retstep=True)
            (x_last, x_prev, P_last, P_prev, Pt_last,
             Pt_prev) = self._traverse_solved_history(
                x0_in, xm1_in, tms, h, T=T, want_dT=True)

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
            z_ss, _info, _ier, _mesg = analysis.fsolve(
                func_autonomous_solved_history, z0, maxiter=maxiterations,
                reltol=_shoot_reltol, abstol=abstol_z, xtol=xtol_z,
                toolkit=self.toolkit, full_output=True)
            x0_ss, xm1_ss = z_ss[:m_], z_ss[m_:2 * m_]
            self.period = period = float(z_ss[-1])
            times, dt = toolkit.linspace(0.0, period, num=npts,
                                         endpoint=True, retstep=True)
        elif self.autonomous:
            ## The period joins the unknowns.  Its residual row is the phase
            ## condition -- in the units of the coordinate it pins, hence
            ## `_tol[phase_k]` -- while the UNKNOWN it adds is a time, whose
            ## own floor has to be a time; mixing the two is the flavour
            ## error F6(a) names, one row further out.
            z0 = np.concatenate((np.asarray(x, dtype=float), [period]))
            abstol_z = np.concatenate((_tol, [_tol[phase_k]]))
            xtol_z = np.concatenate((_tol, [1e-15 * period]))
            z_ss, _info, _ier, _mesg = analysis.fsolve(
                func_autonomous, z0, maxiter=maxiterations,
                reltol=_shoot_reltol, abstol=abstol_z, xtol=xtol_z,
                toolkit=self.toolkit, full_output=True)
            x0_ss = z_ss[:-1]
            self.period = period = float(z_ss[-1])
            ## The grid follows the solved period; everything downstream --
            ## the replay, the waveform, the DFT -- must use it or the
            ## answer is reported on a period the solver rejected.
            times, dt = toolkit.linspace(0.0, period, num=npts,
                                         endpoint=True, retstep=True)
        elif solved_history:
            ## THE SEED IS THE OLD FORMULATION'S ASSUMPTION, written down:
            ## `x_{-1} = x_0`.  It is what the plain path silently assumes
            ## (it seeds both charge rings with the entering state), so an
            ## solved-history run starts where a plain one starts and the
            ## comparison between them is about the SOLVE, not the seed.
            xa = np.asarray(x, dtype=float)
            z0 = np.concatenate((xa, xa))
            tol_z = np.concatenate((_tol, _tol))
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
        rho = None
        if getattr(self, '_monodromy', None) is not None:
            try:
                rho = float(np.max(np.abs(np.linalg.eigvals(
                    np.asarray(self._monodromy)))))
            except np.linalg.LinAlgError:                 # pragma: no cover
                rho = None
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
            self._install_history(x0_ss, xm1_ss, dt)
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
        for t in walk:
            x = self.solve_timestep(X[-1], t, dt)
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

        
            
