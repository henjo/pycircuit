# The shooting analyses: the history kept out of the code

On 2026-09-24 the comments and docstrings of `pycircuit/circuit/shooting/`
were split in two.  The code keeps what it does and what someone calling or
changing it must know now, in the present tense.  This document keeps the
rest, VERBATIM: the dated narratives, the measurements behind each decision,
the alternatives that were tried and rejected, the defects and how they were
found, and the original wording of every paragraph that was condensed.

Each section is one module, and each subsection one member (`### \`name\``).
A docstring or comment whose text moved says so with a line
`History: \`doc/shooting_history.md\`, \`<Class.member>\``; search this file
for the member's name.  The text is the code's own, line for line, with the
comment markers and indentation removed -- read the "before" of a member as
its docstring or comment block at commit `cb42317` or earlier.

Until the same day the whole package was one module, `shooting.py`: see
`pycircuit/circuit/shooting/__init__.py` for where each part went.  The
chronological record of the work is `doc/pss_log_260902.md`.



## `pss.py` -- module level

### `AUTONOMOUS_U_TOL`

The comment above `AUTONOMOUS_U_TOL`, before the move:

⚠ THE SPECTRAL RADIUS CANNOT DECIDE THIS, and trying it first is the
instructive part.  An autonomous orbit gives an eigenvalue at exactly 1
-- but only AT its own period, and a run at any other period reads well
below (measured 0.9615 on the quadrature phase element at the nominal
period against 1.000226 at the corrected one).  Worse, a merely
lightly-damped DRIVEN circuit sits near 1 too: a Q=1000 resonator has
`exp(-pi/Q) = 0.99686`.  So no threshold separates the two -- one
setting misses the autonomous case where users will actually run it, the
other fires on every high-Q filter.

The distinction is structural, not spectral, and it is exact: a circuit
is autonomous when nothing in it depends on `t`.  `u(t)` is sampled
across the period and compared; a phase accumulator driven by a DC
source is autonomous however energetically it oscillates, which is
precisely the case arc 5 asks about.

## `pss.py` -- `PSS`

### class docstring

The class docstring before the move:

Periodic Steady-State using shooting Newton iterations

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

### `parameters`

The comment on `analysis`, before the move:

Sources supply their time-domain waveform only for an
analysis name in timedomain_analyses (('dc','tran')); the
old default 'PSS' matched nothing, so cir.u(t) returned 0
and the whole shooting solve had no excitation.

The comment on `vabstol`, before the move:

1e-6 since 2026-09-19 -- DC, Transient, JAXTransient and PSS share one
meaning and one default; the reason is at `Transient.vabstol`

The comment on `period_column`, before the move:

2026-09-21 (B7c completed): which step lengths move with an
unknown period.  'proportional' rescales every step with T;
'closing' keeps a caller's inner steps at their absolute lengths
and lets the last step close the period.  'auto' is 'closing' on
a caller's grid for an autonomous run and 'proportional' otherwise
-- see `_period_grid` and the note at `_period_column`.

The comment on `theta_ct`, before the move:

⚠⚠ THE ONE KNOB `method='theta'` HAS, AND IT WAS UNREACHABLE.
`_integrator_for` builds `table[method]()`, so every shooting run
took `ThetaIntegrator.DEFAULT_C` -- a RATE, calibrated on ONE
fixture's period.  The transferable quantity is the DIMENSIONLESS
`C T` (see `ThetaIntegrator.DEFAULT_CT` for why `h` cancels), so
that is what this parameter is, and `_theta_biased` turns it into
the rate THIS period needs.  `None` takes the measured knee.

The comment on `steadyratio` and the LTE floors, before the move:

`reltol` MEANS THE SAME THING IN EVERY ANALYSIS: the relative
tolerance of the transient solution.  It is applied to the
per-timestep Newton here exactly as `Transient` applies it, and
nothing rescales it.

`steadyratio` is how the SHOOTING criterion is expressed relative
to it: shooting reltol = reltol * steadyratio, with 1 meaning the
two are equal.  It is >= 1 because the period map is only KNOWN
to the accuracy of the inner solves, so asking the shooting
residual to beat that is asking it to resolve its own noise --
refused rather than silently accepted.  Raise it to accept a
looser periodic steady state for fewer shooting iterations.
The LTE floors, separate from the Newton ones for the reason
`Transient` records: one knob must not move both criteria.  Same
names, same defaults, same meaning -- this analysis reports the
number a transient would have controlled on.

The comment on `method`, before the move:

⚠⚠ THE DEFAULT IS `radau` (owner decision, 2026-09-07), CHANGED
FROM `trap`.  The floor of this stack is DISCRETISATION and it
grows LINEARLY IN Q; at 240 points per period the relative error
in the diffusion constant against the analytic high-Q reference is

    Q      gear        trap        radau
     100   1.79e-03    2.63e-05    6.97e-10
     500   9.02e-03    1.32e-04    3.48e-09
    1000   1.82e-02    2.63e-04    6.97e-09

(240-against-960 grid difference.)  `trap`'s column is its TR-BDF2
twin's `c` (the PPV comes from the twin), `~2.6e-07 Q` at order 3
-- linear in Q like the others.
⚠⚠ The column recorded here until 2026-09-14 (1.49e-06 / 1.04e-04
/ 2.36e-04) was a DEFECT, not trapezoidal: `diffusion_constant`
divided the twin's integral by trap's own period, adding an O(h^2)
term of opposite sign.  That sum was the "error that CHANGES SIGN
near Q = 100" which `grid_error` refused and which was cited for
this default; fixed, trap is estimable.  The default still stands
on accuracy, and more clearly than before (below).
`radau` is order 5, self-starting (no manufactured opener, so no
seam in the period map), L-stable, and carries its own monodromy,
so an autonomous run takes NO TR-BDF2 twin and reads its own
spectrum -- see `monodromy_twin` and `carries_own_monodromy`.

⚠ It IS more expensive per step at a fine grid (3-stage fully
implicit, through the 1-real/1-complex transform): PURE SOLVE time
on van der Pol at Q=100, 480 points, is 3.538 s against trap's
2.586 s.  At a coarse grid it is cheaper (1.088 vs 1.519 at 120),
the shooting Newton needing fewer iterations without an opener
seam in the period map.

⚠⚠ BUT FOR ANY OSCILLATOR SURFACE THE TWIN DOMINATES, AND THAT IS
WHAT SETTLES THE COST QUESTION.  `trap` is not self-sufficient: an
autonomous run must solve a SECOND, TR-BDF2 PSS for its monodromy
(`monodromy_twin`), and `radau` carries its own.  Measured
`diffusion_constant` cost, which pays for that twin:

    npts   trap solve / c-eval   radau solve / c-eval
     120     1.519 / 1.255         1.088 / 0.187   (no twin)
     480     2.586 / 4.161         3.538 / 0.737   (no twin)

⚠ AND AT EQUAL ACCURACY IT IS NOT CLOSE.  Relative error in `c`
against the analytic high-Q reference, with total wall-clock:

    Q~100   trap  480 pts  3.320e-06   6.747 s
            radau  60 pts  7.599e-07   0.908 s
    Q~500   trap  480 pts  1.661e-05  12.434 s
            radau  60 pts  3.802e-06   0.633 s

Radau at SIXTY points beats trap at four hundred and eighty, on
both axes at once.  (trap's errors re-measured 2026-09-14 after the
period-normalisation fix; the old 4.061e-06 / 9.230e-06 were the
defect's, and its "non-monotone" 2.913e-06 -> 4.061e-06 too.)

`trap` remains one argument away for a cheap coarse answer.

### `__init__`

The comment on `irefnode`, before the move:

The reference row is fixed for the analysis, and both the shooting
loop and the Transient this drives need it.  It was recomputed in
every method from a `refnode` argument that no caller ever varied.

The comment on `_period_column`, before the move:

⚠ WHICH PERIOD-COLUMN CONVENTION `want_dT` USES.  'proportional'
is the shipped one: every step scales with `T`, `dh_i/dT = h_i/T`.
'closing' is the commercial one Andreas described -- the inner
transient owns the steps and the LAST one is placed on the period
boundary, so `dh_i/dT = 0` inside and `dh_N/dT = 1` at the close.
MEASURED (roadmap B7c gates 1 and 4): the proportional column is
`O(h)` wrong on a smooth uniform grid and 4.2% RELATIVE wrong on
van der Pol at `mu = 100` (step ratio 16438x), where 'closing' is
46x closer.  Default unchanged pending the rest of B7c.

The comment on `monodromy`, before the move:

⚠ WHICH MONODROMY THE OSCILLATOR SURFACES READ -- see
`monodromy_twin`.  Selects the method that supplies the PPV,
Floquet modes and factored period when a one-step LMM (trap/euler)
solved an autonomous circuit, whose OWN monodromy is first-order
on a limit cycle (the opener seam).  'trbdf2' (DEFAULT): a TR-BDF2
twin on the same grid -- self-starting, no opener, measured 12-32x
more accurate on lambda2 than the Gear-2 twin at practical step
counts (see `monodromy_twin`).  'gear': the former default, a
Gear-2 twin, kept selectable.  'native': the run's OWN plain
factorisation.  ⚠ 'native' UNDER A ONE-STEP METHOD IS THE WORST OF
THE THREE, not an "exact, no twin" escape: trapezoidal's own
monodromy is FIRST order on a limit cycle and its `Q` DIVERGES
under refinement (11.1 / 28.4 / 63.9 vs an exact 5.9083), euler's
likewise -- it is for the gates that measure that defect, not for
results.  gear and trbdf2 runs are self-sufficient (second-order
native monodromy) and ignore this -- they never twin.

### `_is_autonomous`

The comment on the sampling loop, before the move:

⚠ EVERY POINT ON THE GRID, NOT A STRIDE THROUGH IT.  This used to
sample `times[1::len(times)//8]` -- about nine points -- while
calling itself "exact where a spectral test is not".  Nine points
cannot see a narrow pulse: measured on an RC driven by a `VPulse`
positioned BETWEEN two samples, 40% and 20% duty were read
correctly and 5%, 1% and 0.5% all came back AUTONOMOUS.  A clock
misread that way is routed to the free-period system, which
solves for `T` and DISCARDS the period the caller asked for --
and `DEGENERATE_PERIOD_FACTOR` cannot catch it, because it tests
the magnitude of `T`, not whether the circuit was driven.  PWM,
sampling clocks, S/H and mixer LOs are core PSS workload and are
exactly the shapes a stride misses.

The cost is `N` evaluations of `u` ONCE per solve, against `N`
Newton solves in the traversal it decides -- and the loop exits
at the first sample that differs, which is the common case for
every driven circuit.

### `_resolve_break_events`

The docstring before the move:

`break_events`, ON for every method when not given.

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

### `_resolve_x0_unknown`

The docstring before the move:

`x0_unknown`, defaulted from the circuit's TOPOLOGY when not given.

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

The comment on the best-effort block, before the move:

⚠ EVERYTHING BELOW IS BEST-EFFORT AND MUST NEVER RAISE.  This runs
BEFORE `solve` validates its own arguments, so a bad `method` was
reaching `_solves_history` and coming back as `KeyError: 'bogus'`
instead of the `ValueError('method must be ...')` the caller is
owed -- two tests caught exactly that.  A defaulting helper has no
business changing which exception an invalid call raises.

### `_warn_if_the_block_disagrees`

The docstring before the move:

Say so when the TOPOLOGICAL index reads below 2 and the NUMERIC
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

### `_solves_history`

The docstring before the move:

Whether the period map needs the entering history as an unknown.

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

The comment below the docstring, before the move:

⚠ TRAPEZOIDAL CANNOT JOIN THIS FORMULATION BY SOLVING FOR
`x_{-1}`, and trying it is how that was learned.  `iq_{-1}` is
exactly derivable from `x_{-1}` -- `-(i(x_{-1}) + u)`, item 4d --
but the derivative that matters runs the other way: a one-step
companion reads ONLY `iq_{-1}`, so the trajectory depends on
`x_{-1}` solely through it, and `d(iq_{-1})/d x_{-1} = -G` is
SINGULAR wherever a node carries no conductance -- every purely
reactive node, which is most of a resonator.  Adding `x_{-1}` as
m unknowns then makes the 2m x 2m system rank-deficient, and it
fails exactly as it should: `LinAlgError: Singular matrix`, on 25
tests at once.

The right second unknown for a `b != 0` method is `iq_{-1}`
ITSELF -- the `(x, iq)` state its monodromy already uses -- with
the closure `iq_{-1} = iq_{N-1}`.  See item 5's note; not built.

### `_map_kind`

The docstring before the move:

Which period map the PSS's method shoots on -- decided HERE, once
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
plain formulation do it here (`_force_plain_map`).

### `find_initial_solution`

The docstring before the move:

A PROPER initial solution to start shooting from, by pre-integrating
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

### `solve`

The docstring before the move:

Solve for the periodic steady state.

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

The comment at the top of the body, before the move:

ONE SHOOTING SOLVE, IN ITS PHASES (2026-09-24: this was one
1650-line body).  Each phase is a method with its own record; they
share the run's state through `run`.

### `_solve_prepare`

The comment on the hidden-state refusal, before the move:

⚠ HIDDEN STATE IS REFUSED, NOT INTEGRATED AND HOPED OVER.
`TLine.history` is filled by `cir.accept_step`, which the
TRANSIENT calls at every accepted step and which this analysis
never calls -- PSS drives `solve_timestep` directly.  With the
buffer empty `TLine.G`/`u` stamp a DC SHORT, so the line is
silently absent: measured on a quarter-wave open stub, PSS
returned `converged = True`, `spectral_radius = 0.0`, NO warning,
and an amplitude of 0.999969 where a transient gives 0.244201.

⚠ AND FILLING THE BUFFER IS NOT THE FIX.  Calling `accept_step`
per step would populate it and make `phi` genuinely
history-dependent -- so the monodromy would be the derivative of
a neighbouring problem, which is the exact failure `_begin_period`
exists to prevent.  `_begin_period` resets what is IN `x`; that is
the right scope for the integrator rings and the wrong one for
state living outside the vector, and no reset of the rings can
fix a period map that is not a function of `x_0`.  The honest
answers are to admit the delay state into the unknowns (a
different analysis) or to refuse; this refuses.

⚠ AND THE SCOPE OF THAT REFUSAL IS THE ELEMENT, NOT THE CLASS.
Two things look alike and are not.  KUNDERT'S HIDDEN STATE is a
behavioural model carrying internal state the simulator does not
know about -- genuinely broken, and a commercial RF simulator "outlaws [them]
outright".  A DISTRIBUTED COMPONENT has a KNOWN
infinite-dimensional structure described by frequency-dependent
Y/Z/S parameters, and is tractable: "the convolution operation is
diagonalized by the Fourier transform", so the component is
applied spectrally while the STATE stays finite and lumped
(Yang & Phillips, DAC 2002), and an autonomous time-domain
steady-state solve with transmission lines and exact period
derivatives exists independently.

`TLine` here trips the first test because of HOW IT IS
IMPLEMENTED -- a `history` buffer filled by `accept_step` -- not
because a transmission line is unshootable.  The flag is opt-in
per element (`Circuit.hidden_state` defaults False), so nothing
refuses distributed components as a class, and an element that
declared its state properly would pass.  Cited, not verified
here.

The comment on the reference-node check, before the move:

⚠ ONE REFERENCE NODE PER ANALYSIS, CHECKED.  `self.irefnode` is
fixed in `__init__` from `irefnode=` and is what the TRAVERSAL
uses -- `_transient.irefnode`, every `remove_row_col`, the
monodromy's shape.  This local one comes from `solve`'s own
`refnode=` and is what reinserts the zero row into the RESULT.
They were never compared, so `PSS(cir).solve(refnode=b)` solved
against ground and reported against `b`: each row sensible on its
own, the set of them incoherent, with ground itself coming back
non-zero.  Refused rather than silently rotated, because there is
no answer to give -- the two choices disagree about which
variable was eliminated before the solve began.

The comment on clearing the period state, before the move:

⚠ CLEARED BEFORE THE RUN, not after it.  These describe the
period this call is about to solve for; leaving the previous
call's behind would let `factored_period()` hand back an operator
for the LAST solve after this one failed, and `converged` alone
would not catch it -- a second solve that fails leaves the first
solve's `converged=True` nowhere in sight but its state very much
in reach.

The comment on `break_events`, before the move:

Break the traversal's steps at the source discontinuities, when the
method is one whose accuracy that helps -- see `_resolve_break_events`
for the measurement, and `event_grid` for the snap that keeps it from
manufacturing slivers.

The comment on the 'auto' period column, before the move:

⚠ 'auto' IS 'closing' WITH THE PROPORTIONAL POLISH (2026-09-21,
Andreas's call restored on the diagnosis of why it was reversed).
The reversal's evidence -- "closing lands +690 ppm off on the
fold, fails from 2 % off" -- was two things that were not closing:
the polish fired only when the closing step left the
zero-stability bound, so a stretch INSIDE the bound stayed in the
answer (gear on its fold from seeds 0 / 0.5 / 2 % low: +1615 /
+2470 / +3646 ppm at stretches 1.04 / 1.19 / 1.61, proportional
+1431 at every seed; the raw window at 5 %: +3129 ppm at 2.28x,
and its "-33.5 ppm" at 2 % was a cancellation), and the 2 %
failure was the old mid-edge seam whose last step was too short
to absorb the correction.  With the polish unconditional the
closing answer is proportional's at the solved period, seed-
independent to 0.4 ppm (+1408.5 / +1408.6 / +1408.9), and closing
keeps its basin: 5 % low on the fold and 16 % on the window,
where proportional's per-step Newton fails.  'proportional'
stays selectable by name -- one solve, the seed's own grid.
⚠ ON A CALLER'S GRID FOR AN AUTONOMOUS RUN, as the Parameter
says: on a UNIFORM grid 'auto' is proportional -- the closing
column there makes one step of N absorb the period correction
for no reason, and two uniform-grid gear solves that converge
proportionally (the stall fixture at lambda_2 = 0.9, the
periodic-state fold fixture) did not converge closing.

The comment on the phase condition, before the move:

⚠ AND IT IS SUFFICIENT IN PRINCIPLE, DEGRADED IN PRACTICE ON
HIGH-Q: the row removes the singularity from the UNIT
multiplier and only that one, so an oscillator whose other
multipliers cluster near 1 gives a bordered system that is
nonsingular and ill conditioned.  Same cause as the
eigen-selection limit in the class docstring; read it there.

THE PHASE CONDITION pins the coordinate moving FASTEST at the
seed, so the orbit crosses the pinning hyperplane
transversally.  Pin a slow one and the last row of the
bordered Jacobian is nearly parallel to the null direction it
exists to remove, which is a singular system wearing an extra
equation.

⚠ THE `argmax` COMPARES VOLTS WITH AMPERES ON PURPOSE, and
the obvious repair is WRONG.  The row can only remove the
orbit's tangent in proportion to `|e_k . fhat|`; one step of
`|dx_k|` IS `|f_k|` up to `h`, so this argmax is exactly
`argmax |e_k . fhat|` -- it maximises the quantity the row
needs.  Normalising each coordinate by its own swing, which
is what the vector's mixed units invite, was MEASURED on a
van der Pol carrying a VCVS-scaled copy of `v` and picks a
row 704x WORSE aligned (1.4e-03 against 1.0000).  The
scaling that lets a large coordinate win the argmax is the
same scaling that makes it dominate `f`; the two cancel.
Pinned by `test_the_phase_pin_compares_units_on_purpose`.

⚠ A REVIEW DISPUTED THESE CONDITION NUMBERS AND THEN
RETRACTED, and the reason is worth keeping because it is a
trap this file can fall into again.  The reviewer measured
the pin and the orthogonality row as IDENTICAL (edge 1.000x)
and argued a gap was arithmetically impossible, since
`|fhat[k]| = 0.9999` makes `e_k` and `fhat` nearly parallel.
Their harness took the CODE'S analytic `J` and swapped only
the border row, so both readings were the conditioning of
the SAME operator -- and on the plain path that operator is
in the wrong frame (see item 4b's correction below), so its
defect dominated `cond` in both cases.  Identical numbers
were guaranteed by construction.

⚠ THE REASONING WAS ALSO WRONG, INDEPENDENTLY: `cond` is a
function of the WHOLE row, not of its projection on `e_k`.
Alignment 0.9999 makes two rows nearly PARALLEL, not equal --
`e_k` is a unit vector and `fhat` is dense -- and it places
no bound on the conditioning gap.  Switching the phase
condition means solving a DIFFERENT system and living with
ITS conditioning, so the comparison has to build both
systems, which is what the numbers above do.

⚠ AND THE UPGRADE THIS INVITED WAS TESTED AND REJECTED.  An
orthogonality (Poincare) row `<x0 - x_ref, f(x_ref)> = 0`
looks strictly better -- it is the flow-aligned row by
construction, and it cannot pin an unattainable VALUE.
Measured against this rule on the case built to break it
(seeded at `v`'s turning point, so the pin sits 1e-3 of the
way into its coordinate's range): both converge, to the same
answer, at every grid tried, and the bordered condition
number is 1.2e3/3.0e2/8.2e1 for the pin against
2.0e2/6.0e1/3.7e1 for orthogonality at 200/800/3200 points --
a 2-6x edge that never decides anything, and it SHRINKS as
the grid refines.  The row's alignment at the solution stayed
1.0000 throughout: the tangency this was supposed to induce
never materialised.

What IS real: `phase_pin` is a VALUE the orbit must attain,
so a seed far off the orbit can pin one outside its range and
the system is then INCONSISTENT rather than merely hard --
measured on van der Pol at mu=1, seeds at 4x/10x/30x the
orbit amplitude pin `v` at -5.66/-9.50/-15.46 against an
orbit range of [-2.01, 2.01], and all three report ordinary
non-convergence.  ⚠ That is NOT an argument for the
orthogonality row: its plane through the same far seed misses
the orbit too (checked).  It is an argument about SEEDS, and
the remaining gap is diagnostic, not formulational.

The comment on the pinned value, before the move:

⚠ THE RULE IS CANONICAL, AND IT DOES MIX UNITS -- KEPT ON
MEASUREMENT.  Aprille & Trick's oscillator paper, Step 3:
"select k by |f_k(x^i(T^i))| = max_k |f_k(x^i(T^i))|" -- argmax
of the vector field, which is what an argmax over one step is up
to `h`.  `h` is common to every coordinate and cancels; the
UNITS do not (this note used to say they did): the components
are V/s and A/s.  A peer session found the consequence and it
reproduces: on an LC van der Pol with real units (1 nF, 1 uH)
seeded at v's turning point this pins v, nearly along the flow,
and the swing-scaled bordered Jacobian's condition number is
120-160 against 3.75-9.3 with the pin chosen by motion relative
to each coordinate's swing (one Newton iteration saved).  ⚠ BUT
THAT RULE WAS BUILT, MEASURED AND REVERTED (2026-09-23): on the
comparator relaxation oscillator, seeded at ten phases of its
exact orbit (unstaged, 40 iterations), the raw rule converged
from 7/10 and the swing rule from 4/10 -- raw pinned coordinate
1 at nine phases, swing moved it to 2 or 3 at all ten, and the
two fail at different phases (not diagnosed further; Andreas
kept the raw rule on these counts).  Conditioning on a smooth
circuit is not what decides convergence on a switching one.
Both the rule and the decision not to replace it with a
Poincare row have precedent as well as measurement behind them.

⚠ WHAT IS NOT CANONICAL IS FREEZING IT.  Their Step 3 sits
INSIDE the iteration (Step 5 returns to Step 1), so `k` and
the pinned value are re-chosen from the CURRENT trajectory
every iterate -- "note that in this method, an initial k and
w_0k are not required".  Pinning `w_0k = x_0k^i`, the
iterate's OWN current value, is attainable by construction,
which removes the failure mode measured below (a far seed
pinning a value the orbit never reaches) structurally rather
than by advice.  ✅ DONE 2026-09-16 as `phase_rule='reselect'`,
OPT-IN (the default stays `'frozen'`) -- `_phase_row`, and the
measurements, gains AND costs, in `solve`'s docstring.

⚠ THE PHASE ROW SITS OUTSIDE THE INTEGRATOR ON PURPOSE, and
that placement is load-bearing rather than incidental.
Brachtendorf et al. (TCAD 33(6) 867-878) warn of the
alternative: "adding an algebraic equation transforms the
system of (implicit) ODEs to a system of DAEs.  Transient
methods may run into severe problems when the index of a
system of DAEs is two or higher."  This augments the OUTER
shooting system; the inner integration is unaugmented, so the
phase condition cannot raise the index of the DAE actually
being integrated -- which matters here in proportion to how
badly index 2 already behaves (trap and euler both fail to
converge on a V-source/C/C/R loop where gear does).

⚠ THE RULE IS CANONICAL; RE-SELECTING IT WAS TRIED AND
REJECTED, WITH NUMBERS.  Aprille & Trick's oscillator paper
picks `k` by `argmax |f_k(x^i(T^i))|` -- the same quantity an
argmax over one step is, up to an `h` common to every
coordinate -- and their Step 3 sits INSIDE the loop, pinning
`w_0k = x_0k^i`, the iterate's OWN value, so that "an initial
k and w_0k are not required".  That looks like a free repair
for the far-seed failure recorded below.  It is not.

Built and measured 2026-09-02: re-selecting `k` and the pin
from the current trajectory between outer iterations REGRESSED
the working case.  Van der Pol at mu=1 from an ON-ORBIT seed
went from converged to NOT converged, and far seeds wandered
to periods of -52, -1088 and +110 against a true 6.6633.

⚠⚠ WITHDRAWN 2026-09-16: THE REASON RECORDED THEN WAS WRONG, and
so was "taking Step 3 means taking the substitution".  The
zero residual does not make the row empty -- `dz[k] = 0` is
exactly A&T's constraint, and their substitution with `k`
frozen reproduces the bordered iterates to <= 1e-9 (the two
are one linear system).  Re-selection built consistently
converges on-orbit and widens the basin (4x seeds 0/6 -> 6/6)
-- at a cost that keeps it OPT-IN: it fails the grid-aligned
`Idtmod` wrap the frozen pin solves, and it lands on a
different phase of the same orbit.
The 09-02 failure REPRODUCES with a pin taken from the
unknown `x_in` but compared against the manufactured `x_0[k]`
-- the frame error the note below records for the frozen pin.
Built as `phase_rule`; see `solve`'s docstring.

⚠ THE PIN MUST BE IN THE UNKNOWN'S OWN FRAME.  `_x1` is the
state one step AFTER the seed, which is the right thing to
pin when the unknown is `x_in` and `x_0` is manufactured from
it -- and the wrong thing when the unknown IS `x_0`.  With a
fine opening step the two are nearly equal and the mismatch
hides; on van der Pol's own LTE grid, where the opening step
is 1.4845 against a median of 4.6e-04, it pins a value the
orbit need never attain and the solve dies with a bare
non-convergence.

The comment on validating `method`, before the move:

Resolved here as well as in `solve_timestep`, because the SHOOTING
Jacobian depends on which integrator the inner steps used.

⚠ THE NAME IS VALIDATED BEFORE ANYTHING ASKS THE INTEGRATOR A
QUESTION.  `_solves_history` resolves `method` to a class to ask how
far its companion reaches; run first, it turned an unknown name
into a `KeyError` from a dict several frames down, in place of the
`ValueError` this raises.  Two tests caught it, both written for
the class's earlier fall-through defects.

The comment on `solved_history`, before the move:

Whether the entering history joins the unknowns.  Decided once,
here, because it chooses which system is solved -- like autonomy,
and after it, since the two are not composed yet.

### `_shoot`

The comment on the shooting Jacobian, before the move:

THE SHOOTING JACOBIAN FOLLOWS THE INTEGRATOR'S OWN COEFFICIENTS.

Backward Euler's per-step sensitivity is
    dx_n/dx_{n-1} = Jf_n^-1 * C(x_{n-1})/h
-- the COMPANION CONDUCTANCE at the previous point, not the raw
capacitance matrix.  The `/h` was missing, and because C is
singular the accumulated product collapsed to EXACTLY ZERO: the
Jacobian handed to fsolve was `I - 0 = I`, so the "shooting
Newton" was plain successive substitution `x0 <- phi(x0)`.  That
converges at the circuit's own per-period decay -- measured on a
Q=20 resonator as 0.855 per iteration against exp(-pi/Q) = 0.8546,
which is how it was found -- and it never reached fsolve's
tolerance, on any circuit, silently.  With the companion
conductance the same resonator converges in FIVE iterations
(2.64 -> 2.6e-2 -> 1.0e-2 -> 1.9e-4 -> 3.9e-5).

TRAPEZOIDAL'S MONODROMY MUST CARRY `iq` AS WELL AS `x`.  Its
recursion carries `iq` as well as `x`, so the period map is a
function of (x, iq); an x-only monodromy is structurally
incomplete, and measured, using the Euler form for trap converged
SLOWER than no Jacobian at all (0.90 against 0.855 per iteration).

Differentiating the two recursions together,

    iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}
    0    = i(x_n) + iq_n + u(t_n)

gives a propagation of `d(x,iq)/dx0` that costs one extra matrix
product over the Euler form:

    rhs   = Geq_{n-1} Px + Pq
    Px_n  = Jf_n^-1 rhs
    Pq_n  = Geq_n Px_n - rhs

Euler is the SAME recursion with `Pq == 0`: its companion carries
no `iq_{n-1}` term, so dF/diq_{n-1} vanishes and the second row
never enters.  One formula, two methods, which is why this is not
a second code path.
(`newton = True` lived here and fed a message branch that is
gone: it was set unconditionally, so the 'successive
substitution' alternative was dead and the claim it selected was
false anyway -- see the non-convergence warning below.)

The comment on the period map, before the move:

THE PERIOD MAP, per kind -- the only thing about the Newton that
depends on the method (2026-09-23: it was eight residual closures,
one per kind and driven/free period).

The docstring of the nested `_pmap`, before the move:

One period from the unknown `z`: ``(z_0, z_end, M, Mt)`` with
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
their own.

The two comments on the refusals, before the move:

⚠ REFUSED RATHER THAN SILENTLY IGNORED.  A caller asking for
matrix-free on a path that has not got it wants the cost model it
implies; quietly taking the dense route would be a performance
surprise with no symptom, which is the shape of defect this tree
has paid for before.
⚠ REFUSED RATHER THAN IGNORED, like `matrix_free` above.  A
solved-history method already solves for `x_0` and `x_{-1}`
directly and manufactures nothing, so there is no frame to
correct and the flag would be a no-op -- and a no-op flag that
the caller believes changed something is worse than an error.

The comment on damping, before the move:

⚠ THE OUTER NEWTON IS DAMPED, which it was not.  All three
`fsolve` calls took the FULL step with `limiter=None` and no line
search -- a departure from standard practice rather than a
neutral choice: Brachtendorf et al. (TCAD 33(6) 867-878) describe
"shooting, finite difference, or harmonic balance techniques IN
CONJUNCTION WITH A DAMPED NEWTON METHOD" as what is "widely
employed" for limit cycles.  The full step is still tried first
and kept whenever it improves the residual, so a solve that was
converging is unchanged; the halving only runs where the
undamped iteration would have moved uphill.

The comment on the one Newton, before the move:

Find periodic steady state x-vector
Find the periodic steady state: ONE Newton for every kind
(2026-09-23; it was a five-way branch).  The unknown is the entering
state -- gear's PAIR `(x_0, x_{-1})`, seeded `x_{-1} = x_0`, the
old formulation's assumption written down, so a pair run starts
where a plain one does and the comparison is about the SOLVE -- and
on an autonomous circuit the period joins it.  Its row is the phase
condition, in the units of the coordinate it pins (`_tol[phase_k]`),
while the UNKNOWN it adds is a time whose own floor must be a time:
mixing the two is flavour error F6(a) one row further out.

The comment on the unformed monodromy, before the move:

⚠ THE MONODROMY IS NOT FORMED, so it must not be REPORTED
either: `_monodromy` survives from any earlier traversal and
`spectral_radius` reads it without knowing which run wrote it.
(Only the pair's matrix-free paths cleared it; the plain ones
reported the previous solve's radius.)

Gear's autonomous pair and the state-event stage (2026-09-24, Andreas:
"Analyse, fix and create a test for the comparator relaxation oscillator
with gear").  Until then the stage's gate read `_kind == 'stage' or
(_kind == 'pair' and not self.autonomous)`: gear's free-period pair was the
one kind skipped, and silently -- the method check at the top of `solve`
counts gear as event-capable, so nothing warned.  On the comparator
relaxation oscillator gear returned the unstaged orbit: +5.8e-3 of the exact
period at 200 points (-2.9e-3 / +3.6e-3 at 100 / 300, the error set by where
the crossing falls in its step), and the fixed-grid map's dominant
multiplier 52.8 in place of 1.  The stage was already width-generic (gear's
driven pair used it), so opening the gate was the fix: +1.1e-4 at 200 points
(+8.7e-5 / -6.8e-4 / -1.9e-4 / +1.3e-5 at 150 / 100 / 300 / 800), a unit
multiplier.  The same day the stage began running from a first stage that
did not converge (unless it collapsed): unstaged, gear failed at 350 and 400
points -- the fixed-grid map's multiplier 100-290 -- and the stage converged
from its last iterate to 1.5e-4 / 1.2e-4.

### `_report_convergence`

The comment on autonomous oscillators, before the move:

⚠ AN AUTONOMOUS OSCILLATOR CANNOT BE SOLVED AT A FIXED PERIOD,
and this is the only place it says so.

A circuit whose oscillation is self-sustaining -- a VCO
macromodel, an LC or ring oscillator, any phase accumulator
driven by a DC source -- has a one-parameter family of periodic
solutions, because rotating the starting point along the orbit
gives another one.  Its monodromy therefore has an eigenvalue at
exactly 1 and `I - M` is singular AT the true period.  Away from
it the orbit does not close at all: measured on the quadrature
phase element, the discretisation precesses by 2.1e-3 rad per
cycle at 100 steps/period (falling as h^2), so the period map is
a rotation by slightly less than 2*pi whose ONLY fixed point is
the origin -- which is what an unseeded run returns, silently.

Measured either side on that element: |eig(M)| = 0.968 and
sigma_min(I-M) = 2.3e-02 at the nominal period, against
|eig(M)| = 1.000226 and 1.6e-04 at the corrected one.  So there is
no period at which this analysis both has a solution and an
invertible Jacobian, and the answer is not a better seed: it is
the autonomous formulation, which solves for the period jointly
with a phase condition.  Not implemented -- but a run that
returns the origin, or refuses to converge, deserves to be told
why rather than left to look like a tolerance problem.

The comment on `self.autonomous`, before the move:

`self.autonomous` was decided before the solve and chose which
system ran; nothing to re-derive here.  It used to WARN at this
point that a self-oscillating circuit could not be solved at all,
which was true of the fixed-period system and is no longer true
of this one -- the period was an unknown and `self.period` holds
what it came to.

The comment on the non-convergence warning, before the move:

⚠ THIS USED TO BE SILENT.  `fsolve` builds the "No
convergence" message and then discards it whenever
`full_output=False`, which is how this call was written -- so
a shooting solve that never converged returned a
plausible-looking waveform with no diagnostic at all.  It was
non-convergent on EVERY circuit, including a linear RLC whose
answer was visibly close, which is why nobody noticed.
⚠ THIS MESSAGE USED TO CLAIM A "true Newton", AND THAT IS
FALSE ON THE PLAIN PATH.  `newton` is set True
unconditionally, so the alternative branch was dead and every
non-convergence was reported as a Newton failing.  Measured
2026-09-02: the true `dF/dx_in` is SINGULAR on every circuit
tried (rank 1/3, 2/4, 1/3; sigma_min exactly 0), so no method
solves a true Newton in the frame the plain path's unknown
lives in -- it is a contraction, and its residual falls
LINEARLY at a constant ratio.  Advising `method='euler'` on
the strength of a distinction that does not exist sent people
sideways.

The advice that IS backed: the solved-history route has an
exact Jacobian (item 4b) and converges quadratically --
measured 1.69e-01 -> 1.06e-02 -> 3.78e-06 -> 3.48e-12 against
trapezoidal's linear 3.91e-03 -> 3.14e-04 -> 2.66e-05 on the
same circuit.
⚠ THE GEAR ADVICE IS FOR A DRIVEN SOLVE ONLY (2026-09-16).  On an
autonomous oscillator near a unit second multiplier gear is the
method that STALLS (see `_diagnose_lmm_free_period_stall`), so
recommending it there sent a gear user to gear.

### `_replay_orbit`

The comment on the replay, before the move:

THE THIRD LEVEL, MEASURED ON THE WAY OUT.

⚠ WHY THE NESTING WORKS AT ALL, which this docstring stated the
shape of and never the reason for.  Kundert (*Introduction to RF
Simulation*, v2 2003; relayed from the docs session, cited not
verified here): "the strong convergence properties of shooting
methods result from its nature AS A MULTILEVEL NEWTON METHOD, and
not from the fact it is a time-domain method.  Indeed, it is
possible to formulate harmonic balance as a time-domain method
yet its convergence properties do not fundamentally change."

The mechanism is that `phi_T` "is a near linear function ... even
when the underlying circuit is behaving in a strongly nonlinear
fashion, because `phi_T` is evaluated over one period of the
large periodic clock signal" -- the nonlinearity is absorbed by
the INNER transient, which is "a natural continuation method,
quite robust".  So the outer Newton sees a nearly linear map and
the arrangement below is not an arbitrary ordering of three
tests: each level exists because the level inside it has already
made the level outside tractable.

Three convergence criteria stand between a PSS run and its answer,
and the first two are checked while it runs: the inner Newton
(`i(x) + iq + u` under `reltol/iabstol/vabstol`) and the shooting
Newton (`x0 - phi(x0)` under the same, times `steadyratio`).  Both
ask whether an EQUATION was solved.  Neither asks whether the
equation was the right one -- the discrete period map is not the
continuous one, and the gap between them is truncation error.

PSS imposes its grid (`h = T/(N-1)`, uniform, N from `timestep`),
so this cannot be a CONTROL signal -- nothing here may shrink a
step, and doing so would change the period map between shooting
iterations and destroy the monodromy.  It is a MEASUREMENT, taken
on the converged solution over the final replay and reported.

⚠ IT IS THE LEVEL THAT WAS SILENT, and the one that dominates.
On the Q=20 resonator at 100 steps/period all three integrators
report a converged shooting solve, and their amplitudes are
8.815 V (euler), 19.766 V (gear2) and 19.990 V (trap) against
20 V analytic -- a 56% disagreement between two "converged"
answers.  Nothing in the two Newton criteria can see that, because
each integrator solved ITS OWN equations to tolerance.  This
number can: it is `|J^-1 Eg| / (TRTOL (reltol ref + lte_abstol))`,
the quantity a transient would have rejected a step on.
⚠ THE REPLAY MUST OPEN THE WAY THE SOLVE DID, or the waveform is
not the solution: a plain replay of a solved-history answer
would reintroduce the very seam that formulation exists to
remove,
and the reported amplitude would not be the one the residual was
driven to zero on.
⚠ AND IT MUST WALK THE SAME (t, h) PAIRS, which it did not.
The two traversals pair them differently: `_traverse_solved_history`
walks `times[1:]` with `hs[_j]`, while the plain `_traverse` takes
the MANUFACTURING step at `(times[0], hs[0])` FIRST and only then
walks `times[1:]` with `hs[_j]` -- so in the plain case the step
AFTER the opening one uses `hs[0]` again, not `hs[1]`.  This
replay set `walk = times` and indexed `hs[min(_j, ...)]`, which
pairs `times[k]` with `hs[k]` from `k = 1` on and is off by one
against the traversal.

A UNIFORM GRID HIDES IT COMPLETELY -- every `hs` is the same
number -- which is why it survived every uniform test in this
file.  Measured on the Q=20 resonator at 200 points, closure
`|x(T) - x(0)|` of the RETURNED waveform:

      grid      trap            gear (control)
      uniform   5.61e-13        4.62e-14
      4:1       1.70e-02        5.33e-15
      16:1      4.88e-03        1.78e-14

with `converged = True` in every row.  Gear closes on every grid
because it takes the solved-history branch, whose pairing was
already right; the plain path returned a waveform that is not the
solution its own residual was driven to zero on.

Now built as explicit `(t, h)` PAIRS rather than two sequences
indexed in parallel, because the parallel indexing is the bug and
a pair cannot be misaligned by one.
⚠ WHAT A LATER FACTORED REPLAY NEEDS, and the reason it is kept
HERE rather than inside the Newton.  `_traverse_factored*` runs
inside the matrix-free Newton's `_mf_build` closure, whose `steps` go
out of scope with the closure -- and the last `build` call is at
the last TRIAL iterate, which is the converged one only by
accident.  `PAC` wants the factors of the SOLUTION.

So this keeps the four things a replay cannot re-derive from
outside (which seed, which grid, which opening) and
`factored_period()` runs the traversal on demand.  That costs one
traversal for a caller who asks and nothing at all for one who
does not -- where retaining `N` factorisations from every solve
would cost `2 N m^2` doubles on every run, which is the memory
trade `_traverse_factored` documents and most callers never want.

### `_report_lte`

The comment on the three figures, before the move:

THREE NUMBERS, BECAUSE THEY HAVE DIFFERENT REMEDIES.

`max_lte` is the INTERIOR per-step peak -- steps whose estimator
saw only real past charges.  It is exactly the quantity a
transient controls its grid on, and it ranks the integrators the
way their answers rank: on the Q=20 resonator at 100 points per
period it reads euler 0.2876, gear2 0.0763, trap 0.0239 against
amplitudes of 8.815 / 19.766 / 19.990 V (analytic 20 V).

⚠ BUT THE PEAK IS A PER-STEP NUMBER AND A LIMIT CYCLE IS WHAT A
WHOLE PERIOD DOES.  At `reltol=1e-3` euler's peak is 0.288 -- in
tolerance -- while its amplitude is 56% low, because a transient's
criterion bounds each step and says nothing about the 99 of them
that damp the orbit.  `total_lte`, the SUM over the interior, is
the one that sees it: 26.27 for that run against 0.941 for gear2
and 0.340 for trap, tracking the amplitude errors of 55.9%, 1.17%
and 0.05%.  It is an upper bound -- it adds magnitudes, so it
cannot see the cancellation that makes trapezoidal's real error
far smaller than its summed one -- which is the right direction
for a diagnostic to be wrong in.

`max_lte_seam` is the peak over the opening steps of a method
whose COMPANION reads the entering unknown -- Gear-2 here, and
`None` for euler and trapezoidal, which cannot have a seam.  See
`solve_timestep` for why that is the right condition and what the
looser one reported.

⚠ IT IS A FLAG, NOT A MAGNITUDE.  The seam is real -- measured at
1.266e-01 V on the Q=20 resonator at 100 points/period, against
an interior contribution of 1.070e-01, and it is the term that
stops converging (54% of Gear-2's error there, 73% at 400 points)
-- but the NUMBER printed overstates it by orders of magnitude,
because the estimator differences a fabricated charge while the
solution merely reads one.  At 100 points the estimator's
seam/interior ratio is 505x and the answer's is 1.18x.  So use it
to know the seam is there; use `benchmarks/pss_seam_cost.py` to
know what it costs.

The fix is not a smaller timestep -- refining makes its SHARE
grow.  It is to make the entering history part of the shooting
unknowns, so the map is a fixed point in the state a two-step
method actually needs, rather than one that opens off a
stand-in.  Not built; the prize
measured on that resonator is Gear-2's error going 2.34e-1 ->
1.07e-1 at 100 points, at no extra cost per iteration.
An unsound estimate is reported as neither: for trapezoidal that
step is the ONLY one whose number was ever wrong, and dropping it
is what keeps it out of both figures.

The comment on the warning's wording, before the move:

⚠ DO NOT ASSERT CONVERGENCE HERE.  This clause read "the
shooting solve converged, but ..." unconditionally, and the
LTE report is produced whether or not it did -- so a
non-converged run emitted a warning whose first words said it
had converged, directly beside the warning that said it had
not.  Two reviewers read non-converged waveforms as answers
in this file's history; contradictory warnings are not why,
but they are not help either.

### `_check_fundamental`

The comment on the detector, before the move:

⚠ AN AUTONOMOUS PERIOD IS ONLY DETERMINED UP TO AN INTEGER
MULTIPLE, AND THE SOLVE FOLLOWS THE SEED.

`k*T` is a period whenever `T` is, so `x0 - phi_{kT}(x0) = 0` has
solutions at every multiple and the free-period system converges
to whichever one the seed is nearest.  Measured on the quadrature
phase element, whose true period is 1.000e-03: seeds of 1e-3,
2e-3 and 3e-3 return 1.000083e-03, 2.000665e-03 and 3.002245e-03
and ALL report `converged`.  The reported waveform is a correct
periodic solution in each case -- and its fundamental frequency
is wrong by the factor, which is what a PSS user is usually
after.  Nothing said so.

The detector is cheap and needs no extra solve: an orbit
traversed k times comes back near `x_0` partway through.  Grid
points do not land on `T/k` in general (`T/2` at 199 steps is
step 99.5), so this is a NEAREST-APPROACH test against the
orbit's own diameter rather than an equality, and the endpoints
are excluded because every orbit is near `x_0` there.

Driven runs are exempt: their period is the caller's, and asking
for two source periods is a legitimate request, not a mistake.

The comment on the differential state, before the move:

⚠ THE ORBIT IS ITS DIFFERENTIAL STATE, WITH THE PERIODIC ROWS
FOLDED (2026-09-23).  Two points of an autonomous orbit
coincide exactly when their differential states do -- the
algebraic rows are functions of them -- and a phase coincides
modulo its modulus.  Over the whole vector, an idtmod VCO
could never recur (its unfolded phase state advances one
modulus per fundamental), and its wrapped output jumps by a
whole modulus mid-period.  Measured on the free-running
`VcoHdl`: every one-fold solve warned "19.8 times", and true
2- and 3-fold orbits were named 19.9 / 8.65 / 11.96.

The comment on the threshold, before the move:

⚠ THE THRESHOLD IS THE ORBIT'S OWN SPEED AT `x_0`.  A k-fold
orbit passes `x_0` again at the SAME speed, between grid
points, so its nearest point is within about half a step's
displacement there -- and an orbit moving away from `x_0`
is two steps' displacement off by the excluded edge.  The
LARGEST step on the orbit is the wrong scale: it was an
output wrap's jump (0.99 on the VCO) or a fast edge, and at
three times that every point in the first quarter of the
diameter passed.  `_h` is each step of the replay, so
non-uniform grids scale locally.

The comment on the earliest recurrence, before the move:

⚠ THE EARLIEST RECURRENCE, NOT THE NEAREST.  A three-fold
orbit passes close to `x_0` at both `T/3` and `2T/3`, and
`argmin` picked whichever happened to be numerically nearer
-- it reported `2T/3` as "the fundamental", which is wrong by
a factor of two and would have sent the reader to a period
that is itself a multiple.

The comment on the first cluster, before the move:

The closest approach WITHIN THE FIRST cluster: the first
point over the threshold is up to a step early, which
read 3% low and made the multiple look like 2.06 rather
than 2.00.

The comment on the time sum, before the move:

the time from the reference point, summed over the
replay's own steps: with a caller's grid the points are
not evenly spaced, and on the plain path `X[0]` sits one
step before `t = 0`, so `times[_j]` was a step off there

The excluded edges were `max(2, N // 20)` POINTS until 2026-09-24.  Gear's
grid after a landed state event opens with ten doubling steps from 1e-5 T,
so six points excluded 3e-4 of the period, and the orbit still leaving `x_0`
-- its displacement about the speed times the step, the elapsed time being
about the step -- read as a 3182-fold recurrence (the comparator oscillator
at 100 and 150 points).  Now 5 % of the period in time: the same points on
a uniform grid.

### `_orbit_results`

The comment on the first entry, before the move:

⚠ THE PLAIN PATH'S FIRST ENTRY IS A SEED, THE OTHER PATH'S IS A
SOLUTION.  Plain takes N steps from `x0_ss` and reports their
results; the other starts AT `x_0` and takes N-1, so dropping the
first would drop a real point and shift the waveform by a step.
⚠ THE FIRST ENTRY IS DROPPED ONLY WHEN IT IS NOT PART OF THE
PERIOD.  On the default plain path `X[0]` is `x_in`, the
pre-image of the manufactured step, which sits one step BEFORE
t=0 and is not a point of the orbit.  With `x0_unknown` -- and on
the solved-history path -- `X[0]` IS `x(0)`, so dropping it both
discards a real sample and leaves the waveform one column short
of `times`.

The comment on the non-uniform DFT, before the move:

⚠ ON A NON-UNIFORM GRID (`grid=`) THE INDEX DFT ABOVE IS NOT A
FOURIER COEFFICIENT -- measured 7.5 % off in the carrier and not
converging on a 3:1 grid.  Same layout and RMS fold, taken as the
trapezoid-weighted sum at the true times (`_period_quadrature`'s
weights; uniform grids never reach this branch).

The comment on the quadrature rule, before the move:

the SAME rule `_period_quadrature` gives every consumer
(2026-09-21): a periodic cubic spline on an event-free grid, a
piecewise one breaking at the landed event nodes -- `fpss` and
`carrier_phasor` are pinned equal to 1e-12, which a second
rule here broke

The comment on `fpss` being RMS, before the move:

⚠ `fpss` IS RMS, AND THE USUAL THING TO COMPARE IT AGAINST IS NOT.
`freq_analysis` returns an RMS, energy-folded, positive-frequency
spectrum (see the note above).  A commercial simulator's frequency
-domain PSS output is conventionally PEAK, so a harmonic read from
one and compared against the other differs by `sqrt(2)` -- 3.01 dB
-- with nothing in either result announcing it.  Multiply `fpss`
by `sqrt(2)` for a peak-convention comparison, or divide theirs.
Recorded because this campaign has lost time to factor-of-two and
factor-of-pi convention defects more than once, and a 3 dB offset
is small enough to be mistaken for a modelling difference.

### `_closing_polish`

The comment on the second pass, before the move:

⚠ THE SECOND PASS (2026-09-21).  'closing' keeps a caller's inner
steps where the transient validated them, so the free-period
Newton converges from a seed period 16 % off where 'proportional'
fails its per-step Newton (relaxation van der Pol, mu = 10, its
own `lte_grid`).  But the closing step then absorbed the whole
period correction -- 3 s against a 0.7 s neighbour, a growth far
beyond a two-step method's zero-stability bound -- so the
integrator dropped it to Euler and the transposed replay refused.
Closing is the BASIN device; once converged, the grid is
re-fractioned at the solved period and solved once more
proportionally from the converged state, which is the sane grid
the answer is reported on.  ⚠ ALWAYS, NOT ONLY BEYOND THE BOUND
(2026-09-21).  This used to polish only when the closing step
left the zero-stability bound, on the reasoning that a seed
already close needs no second pass -- but ANY stretch of the last
step is a change of discretisation and stays in the answer: gear
on its fold, seeds 0 / 0.5 / 2 % low, +1615 / +2470 / +3646 ppm
at stretches 1.04 / 1.19 / 1.61 where proportional reads +1431 at
every seed and the polished answer +1408.5 / +1408.6 / +1408.9.
The seed-exact case stretches too (the discretisation's own
period error is absorbed by the last step).  Cost: one Newton
from a converged state, one or two iterations.

The comment on the caller's fractions, before the move:

⚠ ON THE CALLER'S FRACTIONS, not the closing-distorted grid:
re-fractioning THAT grid keeps the giant last step (measured:
trbdf2 -1.1 % in period, c -97 %, second pass or not)


## `_pss_newton.py` -- `_ShootingNewton`

### `TRIVIAL_ORBIT_FACTOR`

The comment on the constant:

The equilibrium test above: a returned state whose DC residual is
within this factor of `iabstol` is an equilibrium, not an orbit.  1e3
because the shooting Newton's own residual sits at `abstol`, and an
orbit's DC residual at t = 0 is a CIRCUIT-scale current (measured:
the van der Pol at 2 V has |i + u| ~ 1 A there against 1e-27 for the
collapsed state -- 27 orders apart, so the factor is not delicate).

### `_free_period_solve`

The docstring paragraphs before the move:

Measured on BOTH autonomous elements in the tree, so it is a
property of the formulation and not of any circuit: from a 1e-4
seed against a 1e-3 fundamental, Gear-2 returned -1.5e-20 on the
quadrature element and 3.9e-19 on the scalar `Idtmod`, and
trapezoidal raised a bare `LinAlgError` from three seeds of five as
its Jacobian went singular on the way down.

Neither outcome is a silent wrong answer -- the collapse reports
`converged = False` (⚠ ENFORCED HERE, by demoting `ier`; asserting it
in prose was not enough -- see the note at the demotion) and the
exception is loud -- but both told the
user nothing about the cause, and the generic non-convergence
advice ("raise maxiterations") is wrong for it: no number of
iterations reaches a fundamental from below.

The comment on the second trivial root:

⚠⚠ THE SECOND TRIVIAL ROOT, found 2026-09-08 while gating
`warping_estimate` on an index-2 oscillator: the guard below
catches `T -> 0`, and an autonomous solve has ANOTHER root that it
cannot see -- the EQUILIBRIUM, `x(t) = x_dc`, which is periodic at
EVERY `T`.  Radau seeded 10 % below the fundamental on the
index-1 van der Pol and on the index-2 fixture returned
amplitude 0.0000 (state 1e-27) at a period near the seed with
`converged = True`; trapezoidal on the same seed failed honestly.
`T` stays finite, so the period test passes, and the periodicity
residual is exactly zero because the equilibrium IS periodic --
at every `T`, so it is a whole LINE of roots in `(x0, T)`, not a
point (peer's sharpening): that is why the residual is exactly
zero rather than small, why no residual-based guard could have
caught it, and why the DC residual does in one evaluation.
The test that sees it is the DC residual of the returned state:
an orbit has `C x' != 0` somewhere at t = 0, so `i(x) + u` is far
from zero there; an equilibrium has it at solver tolerance.

The comment on the demotion:

⚠⚠ THE COLLAPSE MUST BE DEMOTED HERE, and for two turns of this
record it was not.  The docstrings above and on `solve` both
asserted "the collapse reports `converged = False`" -- and
NOTHING ENFORCED IT.  `self.converged` is `(_ier == 1)` and
nothing else, while `T = 0` is a REGULAR root: the solver
reaches it cleanly and reports success, so Gear-2 returned
`T = 5.42e-18` with `converged = True` on a circuit with no
orbit in it.  The warning fired correctly the whole time, which
is exactly what made this survive -- a reader who checks the
documented flag instead of catching warnings got `True`.

Demoting `ier` rather than assigning `self.converged` is
deliberate: all three autonomous call sites already feed this
return value into `self.converged`, so one demotion covers the
plain, solved-history and matrix-free paths, and any future
path inherits it by construction.  `ier = 5` is `fsolve`'s
"not making good progress" code -- the closest existing
meaning, and already handled everywhere `ier` is read.

### `_diagnose_lmm_free_period_stall`

The docstring paragraphs before the move:

⚠⚠ MEASURED 2026-09-14 (peer report, reproduced).  On a weakly limited
LC oscillator (second Floquet multiplier 0.99, set by the limiting),
Gear-2's autonomous shooting solve stalls at a residual FLOOR:
converged quadratically at 925..1600 points per period, never at 900
or 800.  Every property that would point at the solver was ruled
out -- the residual is a pure function of the unknowns; damped,
Armijo (20 halvings) and Levenberg-Marquardt steps from the stall all
stop at the floor; tighter inner tolerances change nothing; the grid
is exactly uniform; a scan of the Jacobian's two weakest directions
finds no lower point.  The discrete periodic solution of the
solved-history system `(x_0, x_{-1}, T)` CEASES TO EXIST below a grid
threshold: `sigma_min(J)` at the root falls linearly toward it (to
zero near 904 points), the discrete second multiplier stays at 0.990,
and the vanishing direction lies ~99 % in the `x_{-1}` block.  The
threshold grows as the multiplier approaches 1 (no convergence at
6400 points at 0.999).  Trapezoidal stalls too; `radau` -- no history,
no manufactured opening -- converged in 7-10 evaluations at 800
points at both multipliers.

So no number of iterations, damping or tolerance fixes THAT stall, and
the generic advice would send the caller the wrong way.  This measures
what it can at the last iterate (one residual evaluation, only on a
failed solve) and says so.  Stage methods are skipped: they did not
show this.

⚠⚠ THE "ITERATIONS DO NOT HELP" CLAIM WAS OVER-BROAD AND IS NARROWED
(2026-09-16, peer report + reproduced here).  It was measured for
GEAR-2's solved-history stall and then written as though it held for
every multistep free-period solve.  It does not: a trapezoidal solve at
multiplier 0.9 converges with a bigger budget on the reporter's
oscillator (200 iterations) and at the DEFAULT 25 on the van der Pol
fixture here.  What IS measured at 0.99, on the reporter's tank and
rebuilt independently: trapezoidal's residual RISES with budget --
6.765e-07 at 25 iterations, 9.011e-07 at 200 -- and the analytic step
is UPHILL against the true derivative, so the line search's halvings
cannot improve it.  `radau` converged on that same fixture in 25
iterations (period 6.2831863e-07, amplitude to five digits).

⚠ AND THE DISCRIMINATOR IS NOW REPORTED RATHER THAN GUESSED:
`analysis.fsolve` counts the iterations whose step it could not improve
(`infodict['ls_unimproved']`) and this message names the count.  Zero
means the solve was descending and a bigger budget is worth trying;
non-zero means the budget would repeat an uphill direction.

⚠ `x0_unknown=True` IS A DIAGNOSTIC HERE, NOT A FIX.  In that frame the
analytic Jacobian IS the derivative (worst entry 6e-07 against 1.000 in
the default frame, both measured by central differences on the
reporter's fixture) and the step becomes a descent direction -- but the
solve still does not converge there: `||F||` plateaus at 3.069e-05,
identical at 25 and 200 iterations.  The frame explains the behaviour;
it does not rescue the case.

The state-event clause (2026-09-24).  On the comparator oscillator the
diagnosis blamed weak damping ("second Floquet multiplier near 1") where the
second multiplier is 0.002: the cause there is the unstaged map across a
switch sharper than the grid, and the staged solve converges.

### `KRYLOV_TOLERANCE_FACTOR`

The comment on the constant:

How hard GMRES is asked to solve, relative to the shooting tolerance.
An inexact Newton only needs the step accurate enough not to spoil the
outer convergence; measured k is 2-12 on circuits whose `I - M`
clusters at 1 (the fast modes decay over a period, leaving the slow
ones), so k tracks the number of SLOW MODES, not m.

### `_matrix_free_newton`

The docstring paragraphs before the move:

pair, and the bordered autonomous versions of each (one builder
since 2026-09-23, `_mf_build` in `solve`).

MEASURED (moved here from the solved-history driver, 2026-09-23).
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

Continuing after the convergence-test paragraph:

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

The comment on the inner solve's verdict:

⚠ THE INNER SOLVE'S VERDICT IS NOT DISCARDED.  It used to be,
and a Krylov breakdown then surfaced as the generic outer
'No convergence' with nothing naming the cause -- in a file
whose whole standard is that a failure says what happened
(`T = 0`, the trivial root, the singular free-period
Jacobian).  An unconverged GMRES makes `xdiff` a direction
the Newton has no reason to trust, so the outer loop is told
to stop rather than iterate on it.


## `_pss_grids.py` -- `_PeriodGrids`

### `event_grid`

The docstring paragraphs before the move:

A6/B7b.  `Transient` breaks its steps at `cir.next_event`; the PSS
traversal does not, so a pulse edge inside a step is integrated straight
through.  This returns step FRACTIONS for `solve(grid=...)` with each
event in the period placed ON a grid point -- by SNAPPING the nearest
point onto it when one is close, and INSERTING otherwise, so no
arbitrarily small step is ever created (the B7c lesson).

Measured on an RC driven by a `VPulse`, against a 4000-point reference,
a 40-step uniform grid versus the same grid with its 3-4 event times
landed on::

    edge offset   uniform      + events        gain
    td = 0        6.787e-03    8.032e-04       8.5x
    td = 0.0125T  5.720e-03    2.099e-04        27x
    td = 0.0092T  2.483e-03    7.910e-04       3.1x

The rest of the time-driven-events paragraph, and the next:

⚠ THE REASON RECORDED HERE WAS WRONG, AND IS CORRECTED (2026-09-06).
This used to say the wrap time "has to become an unknown the Newton
solves for", citing a map "discontinuous by |dphi| ~ 8.2e-3 on a grid
point".  Measured: that jump is **grid-INDEPENDENT** -- 1.414214e+09
at seven different grids, with the wrap ON a node and OFF it alike,
and unchanged when the exact wrap times are added to the grid.  A
quantity that does not move when the grid moves is not a grid
artefact.  The real defect was the fold at the period ENDPOINT, in the
OUTPUT map, and it is fixed in the RESIDUAL by `_fold_periodic`.
The conclusion for THIS method is unchanged and now for the right
reason: a state-dependent reset is not a grid feature, so `event_grid`
does not help it -- but nor does it need the traversal surgery that
sentence implied.

⚠ AND THIS IS NOT SALTATION.  Saltation was measured and falsified twice
for this codebase (a switched conductance, a discontinuous injection and
an `Idtmod` wrap all give a monodromy-vs-FD gap falling at 2.00x per
doubling, i.e. O(h)); each step already uses its own `Jf`/`C`, which
describe whichever side of the switch that step is on.  The problem was
only ever that the grid could not BREAK at the event.

The comment on the snap:

⚠ A SNAP MAY MOVE A GRID POINT ONTO AN EVENT, NEVER AN EVENT
ONTO AN EVENT (2026-09-20).  A `tr = 0` pulse is clamped to a
`Pulse.MIN_EDGE` = 1e-18 ramp, so each edge is TWO events 1e-18
apart.  The snap used to land the first and then overwrite that
node with the second, collapsing the ramp onto ONE node on the
post-jump side -- and the step arriving there integrated its
whole length with the post-jump source.  Measured on the pulsed
RC (analytic reference), edges landed, N = 100 .. 1600: EVERY
method first order, gear 1.3e-3 / trap 9.5e-4 / radau 3.0e-4 at
1600 halving per doubling, the error injected AT the edge node
(radau's = its endpoint weight 0.111 x h dU/tau).  With both ramp
ends kept as nodes -- a 1e-18 step -- radau is exact (8e-10), trap
second order, and gear second order too (3.7e-4 -> 1.0e-5 over
100 -> 800): variable-step BDF2 with h_n/h_{n-1} -> inf over a
consistent tiny step degenerates to the trapezoidal rule, the
"parasitic" factor w/2 multiplying a difference that is itself
O(1/w).  So the tiny step is the ramp and it stays; the B7c
lesson ("no arbitrarily small step") is about grid points near an
event, not about two events.

### `LTE_GRID_RELTOL_MAX`

The comment on the constant:

: the coarsest reltol the adaptive run of `lte_grid`
: takes by default -- see the table at `lte_grid`

### `FOLD_MAX_STEP`

The comment's first line:

the fold's own resolution (2026-09-23): no step above this fraction

### `_fold_periods`

The comment on the driven boundaries:

⚠ A DRIVEN CIRCUIT'S PERIODS ARE THE DRIVE'S (2026-09-21): its
sources are functions of absolute time and its edges sit at
fixed phases of the drive, so the boundaries are k T exactly
and phase 0 is the drive's t = 0.  Folded on a CROSSING instead
(the oscillator rule), the switched-capacitor sampler's finest
steps landed at phases 0.26 .. 0.31 of the drive with the
switch edges at 0.75 -- the grid half a period off where it
mattered, and the default path.

The comment on the straddling step:

⚠ A STEP THAT STRADDLES THE BOUNDARY KEEPS ITS OWN WIDTH
(2026-09-21).  Clipped to its sliver inside the period, the
per-phase MIN over 24 periods whose natural steps drift
against the boundary drove the first and last bins to
0.0009 T on the switched-capacitor sampler (natural steps
there 0.034 T) -- a fine seam on a flat part of the drive.
The crossing seam of an oscillator sits on an edge that is
already fine, which is why the autonomous fold did not
show it.

The comment on the seam and the source corners:

the grid starts at the phase of COARSEST density -- the slow
branch.  At the crossing itself the density jumps ~4x (coarse steps
approach an edge, fine ones follow it: a 0.25 seam ratio).  At the
FINEST phase the seam sits mid-edge, where the PPV's sensitivity
is: measured on the same grid, gear's c +28 % against -3.9 % with
the seam on the slow branch, trbdf2's +8 % against +0.4 %, the
period the same either way (+23 / -28, +10 / +7 ppm).  The coarse
opener is `_period_grid`'s doubling ramp.
(a driven grid keeps phase 0 at the drive's t = 0: `solve` reads
its fractions from there, and the landed events with them)
⚠ THE FOLD OWNS ITS RESOLUTION AT A SOURCE CORNER (2026-09-23).
It used to inherit it by accident: a transient stepping OVER a
pulse corner left a cluster of rejection-driven tiny steps after
every edge, and the fold kept them.  Once the transient landed
its corners (the Runge-Kutta loop, item 2 after E7) that cluster
was gone -- radau, L-stable, damped the RC's decay after the edge
with a quiet estimate and 0.03 T steps, then doubled through the
flat hold phase up to 0.23 T -- and the PSS grid built from the
fold read its first harmonic 4.5e-2 off (was 5.3e-4) with the
period quadrature's spline weights at -0.24 / +0.36.  So: the
step density is funnelled towards every source corner (a step at
phase distance `s` from a corner is at most `s`, never below one
bin) and capped at FOLD_MAX_STEP of the period -- a period
integral needs points where a transient's tolerance does not.
Measured on that grid: H2 1.15 -> 3e-3 from the cap, H1 4.5e-2 ->
1.1e-3 from the funnel (a ramp from an eighth of the edge width;
from a sixteenth, 2.2e-4), the pin being 2e-3.

The comment that stood before the walk (it described the overshoot-and-rescale walk the next comment replaces):

the walk overshoots the period by at most one step and the whole
grid is then rescaled to sum 1: every ratio is preserved and no
tiny closing step is made (a closing remainder used to trip the
closing-step bound and fire a needless second pass)

The comment on the walk's end:

⚠ THE WALK ENDS EXACTLY AT THE PERIOD (2026-09-21, item 3
of the non-uniform-grid list).  It used to overshoot by
up to one step and RESCALE every fraction to sum 1 --
which shifts every phase by up to a coarse step times
its phase, so the grid's fine regions landed off the
orbit's edges by up to 2 % of T: six folds of vdP mu = 10
(205 points, the same edge groups to the point) split
into two clusters, gear +1409 .. +1563 ppm against +273
.. +302, trbdf2 on the same grids +150 .. +173 against
+44 .. +46 -- the "gear swing" of the record was this,
the grid's alignment, and it was 5x for every second-
order method.  The remainder is its own step, merged into
the last one when it would be a sliver and the merged
step keeps the controller's own 2x growth -- spread over
the trailing steps, each grown to at most 2x the one
before it, so no sliver and no ratio beyond the bound
(the tail is the seam, the coarsest phase: a step there
grown by a fraction of itself costs nothing measurable).

### `_crossing_chain`

The docstring paragraph before the move:

⚠ A STATE THAT CROSSES ITS MIDLINE MORE THAN ONCE PER PERIOD
(2026-09-21, item 2 of the non-uniform-grid list).  The detector
took EVERY rising crossing of the fastest-swinging state, and a
strong third harmonic (`sin + 1.2 sin 3`, three crossings) or two
pulses of different height per period (two crossings) spread the
spacings past the 1 % rule: `_observed_period` returned None, the
fold was refused with the "no consistent recurrence" warning and
the grid fell back to the single window.  Not silent, but a
capability lost on exactly the harmonic-rich oscillators whose
grids need folding.  The hint names the period the caller expects
and picks the branch; a single-crossing state gives the same chain
as before, bit for bit.

### `lte_grid`

The docstring paragraphs before the move:

B7a.  A transient adapts because it cannot see the future; PSS
re-solves the SAME interval over and over, so it can be handed a
grid that was chosen well ONCE and then frozen.  This is the
derivation side; `_period_grid` is the consumption side and has
been shipped since item 5.

Continuing after the usage example:

MEASURED on van der Pol at `mu = 100`
(`benchmarks/pss_lte_grid.py`, the gate this was promoted from):
1105 derived steps converge where 1105 UNIFORM steps do not, and
beat a 20000-point uniform grid -- 18x fewer points and -47.3 ppm
against -60.6.

⚠ ON A DRIVEN CIRCUIT (measured 2026-09-21, the switched-capacitor
sampler, sine and pulse clocks, each method's fold against uniform
grids of the same and double count, radau uniform-3200 reference):
the waveform is 6.5-40x better on the fold at the same count (radau
95 points 3.7e-5 of the swing vs 2.4e-4; trbdf2 261: 1.5e-4 vs
2.1e-3; gear 125: 1.2e-3 vs 4.5e-2) and the fold at N beats uniform
at 2N; a uniform grid puts ONE step across a clock's ramp whatever
N, the fold's finest steps sit inside it.  Held noise on the fold is
at kT/C for trbdf2 and radau.  ⚠ NOT FOR GEAR'S NOISE: gear's held
variance on its own fold reads 0.78 x kT/C against 0.98 uniform at
the same count.  The fold resolves the STATE -- the tracking phase,
where the state is flat, gets h ~ tau -- and gear's covariance
recursion carries its recorded O(h/tau) tracking floor through the
turn-off.  A noise quantity with its own time scale is resolved only
where the state's grid happens to be fine; use trbdf2 or radau for
noise on a folded grid, or a uniform grid for gear's.

⚠⚠ IT IS FOR STIFF SMOOTH PROBLEMS AND NOT FOR EVENTS, and that
boundary is measured rather than cautionary.  On a wrapping
`Idtmod` the derived grid is WORSE than a uniform grid of the same
count -- max LTE 2.64e+05 against 1.67e+05 times tolerance at
~1429 steps -- because the LTE peak sits at the RESET on every
grid, and no step size makes a discontinuity's local truncation
error small.  The event half of B7 is a different item: it needs
the event time to be an unknown the Newton solves for, because a
grid frozen from a PAST traversal cannot represent an event whose
time MOVES as the Newton iterates.

Continuing after the fractions and `tstab` paragraphs:

⚠ TWO THINGS ABOUT THE WINDOW, BOTH FOUND WHEN GEAR COULD NOT RUN
ITS PPV ON THE GRID GEAR HAD PRODUCED (2026-09-21, relaxation van der
Pol, mu = 10): the window ends at the last NATURAL step, not at the
`tend`-landing step `Transient.solve` truncates (a tenth of its
predecessor, hence a 10.5x growth across the period seam, beyond a
two-step method's zero-stability bound -- the integrator dropped two
steps to Euler and the transposed replay refused); and the cut is
ROTATED to the node whose seam ratio is closest to 1, the seed moving
with it.  The interior ratios are the adaptive run's own (growth
clamped at 2).  ⚠ AND THE PERIOD YOU PASS MUST BE CLOSE: the grid is
fractions of it, so with `period` 16 % off the fine regions sit 16 %
away from the edges and gear's and trbdf2's per-step Newton fail on
the coarse steps that land there; seeded within a per cent (a coarse
uniform solve, or the transient's own recurrence) both converge.
Measured on that fixture, 195 points: gear's period -519 ppm (second
order under 2x splitting; a uniform gear grid of the same count
-7214 ppm), trbdf2's -66 ppm; gear's diffusion constant 52 % high
with its unit multiplier 0.10 off the circle (see `ppv`'s warning),
trbdf2's 16 %.  On a relaxation oscillator's adaptive grid trbdf2
uses gear's own grid better than gear does -- by a steady ~8x, the
methods' constants; the 5x swing between folds of one orbit that
the record blamed on gear was the fold's own rescale (2026-09-21,
`_fold_periods`), and every second-order method paid it.

⚠ THE ADAPTIVE RUN NEEDS A TIGHTER TOLERANCE THAN THE PSS'S OWN
(2026-09-21, "Do 2").  Measured on the same fixture, each method on
the FOLDED grid its own run made, period error (ppm) and diffusion
constant error against radau N = 8N::

    reltol   gear (pts, ppm, c)         trbdf2                  radau
    1e-4     95    -900    +179 %      146   -2.0    +1.0 %     70   -3.9    -1.9 %
    1e-5     207   +20     +22 %       290   -7.7    +0.27 %    116  -0.17   -0.04 %
    1e-6     425   -3.6    +7.8 %      602   -1.0    +0.05 %    191  -0.004  -0.002 %
    1e-7     928   -1.6    +2.6 %      1294  -0.25   +0.01 %    331   0.000  -0.0002 %

The PSS default reltol 1e-4 makes a grid gear cannot use (-900 ppm,
c 2.8x); 1e-5 is the coarsest that serves every method at the 10-ppm
level, so the default is `min(reltol, LTE_GRID_RELTOL_MAX)`; gear's c
within a few per cent wants 1e-7.

⚠ TO REPAIR A COARSE SOLVE, CALL THIS FROM ITS CONVERGED STATE WITH A
SHORT `tstab` (2026-09-21; `PSS.refine_grid` was deleted on this
measurement).  The fold needs `LTE_FOLD_PERIODS` (24) settled
periods, not the default 200, and a converged state is settled::

    fracs, seed = pss.lte_grid(T, x0=x_converged,
                               tstab=(PSS.LTE_FOLD_PERIODS + 1) * T)

Measured on the same fixture from a converged uniform-200 solve
(radau +4.8 ppm, gear -6264): radau 114 points -0.3 ppm in 7 s, gear
204 points +28 ppm in 3 s -- the same grids a 200-period run makes.
`refine_grid`, one subdivision pass from the same solves, reached
the same accuracy at 1120 (radau) and 4617 (gear) points with step
ratios of 4 .. 6, because it subdivided every cell uniformly at the
finest step the transient wanted anywhere inside it and could
neither move nor remove a point; from a decimated fold likewise
(57 -> 1077, 102 -> 4507).  Ten to twenty times the points for the
same answer, and grids gear's replays refuse: it had no niche left.

The comment on the adaptive run's integrator:

⚠ THE PSS'S OWN METHOD SHAPES THE GRID (2026-09-21, "Do 2"): this
used to build `Transient(...)` with no integrator, i.e. the
Transient default `Gear2Integrator()` whatever `method` the PSS
was given -- a gear-shaped grid for a radau PSS.  The adaptive run
steps with the PSS's integrator, so the LTE profile the grid
freezes is the one the PSS will pay.

The comment on the measured period:

⚠ THE PERIOD IS MEASURED FROM THE RUN, NOT TAKEN FROM THE HINT
(2026-09-21).  The grid is FRACTIONS of the period, and with the
hint 16 % off (my relaxation estimate at mu = 10) the fine regions
sat 16 % away from the edges and every method's per-step Newton
failed; the free-period solve then needed the closing convention
and a second pass just to recover.  The adaptive run already
contains the period: the rising crossings of the fastest-swinging
state over the last periods, refined by linear interpolation, read
19.1003 against a reference 19.0986 (+9e-5, the transient's own
discretisation) from that 16 %-low hint, consistent across
crossings to 5e-5.  The window is cut at THAT period, it is left in
`self.lte_period` for `solve(period=...)`, and a hint more than
1 % off is WARNED.  An inconsistent detection (spacings spread
above 1 %, or fewer than two) keeps the hint, with a warning.

The comment on the fold:

⚠ THE GRID COMES FROM EVERY SETTLED PERIOD, NOT THE LAST WINDOW
(2026-09-21, Andreas: "Do 1").  A frozen window inherits whichever
rejection-and-growth pattern that one period got: two windows of
the same run parameters read gear -1461 / -22 ppm and trbdf2 -180 /
+10 on a relaxation van der Pol (mu = 10) -- an 8x lottery for
every second-order method (radau, order 5, does not care).  With
the period detected, each of the last `LTE_FOLD_PERIODS` periods'
accepted steps is folded onto a common phase (phase 0 at the
rising crossing the detector found), the per-phase MINIMUM of the
local step over periods is the density, and a grid is re-meshed
from it with growth capped at the controller's own 2x.  Measured on
the same run: gear +40 ppm, trbdf2 +10, radau -0.01 at 207 points
(the 25 % quantile +267 / +13, the median +957 / +16, the last
window re-meshed +150 / +6).  The seam sits at the COARSEST phase
(the slow branch; a mid-edge seam cost gear's c +28 %) and the
opener is `_period_grid`'s doubling ramp; the seed is the
interpolated state at that phase of the last full period.  A
driven circuit folds on the drive's own boundaries k T (phase 0
at the drive's t = 0), an oscillator on the rising-crossing chain
`_crossing_chain` picks, one per period of the hint.  `fold=False`
keeps the single-window cut; a run without a consistent
recurrence falls back to it, warned.

The comment on the window's end:

⚠ THE WINDOW ENDS AT THE LAST NATURAL STEP, NOT AT `tend`
(2026-09-21).  `Transient.solve` lands its final step exactly on
`tend`, truncating it -- on a relaxation van der Pol at mu = 10 to
a tenth of its predecessor -- and a periodic grid then carries a
10.5x GROWTH out of that step that no controller clamp ever saw:
beyond the two-step zero-stability bound, so the integrator
dropped two steps to Euler and gear's transposed replay refused
the grid gear itself had produced.  Rotating the seam only moved
the step into the interior.  The last accepted point before the
landing step ends the window instead.

The comment on the rotated cut:

⚠ THE CUT IS ROTATED TO THE SANEST SEAM (2026-09-21).  The window
is the last period of the adaptive run, so its seam -- the last
step against the first, which a PERIODIC grid joins -- fell
wherever `t[-1] - T` happened to land: on a relaxation van der Pol
at mu = 10 that was right after a rejected tiny step, a 10.5x
GROWTH across the period boundary.  That is beyond a two-step
method's zero-stability bound (1 + sqrt 2), so the integrator
dropped the first step (and the one after it, rebuilding history)
to Euler, and gear's transposed replay refused the two-alpha
steps -- gear could not run its PPV on the grid gear itself had
produced.  The interior ratios are the run's own (growth clamped
at 2 by the controller; shrinks unconditionally stable), so the
cut moves, when the seam is outside the bound, by the least to a
sane one, and the seed moves to that node.
⚠ ROTATED ONLY WHEN THE SEAM DEMANDS IT, and then by the smallest
rotation to a sane seam: an unconditional move to the flattest
seam put the seed at a phase from which the free-period Newton,
started 8 % off in period, ran away (mu = 4, the tree's own
`lte_grid` test) -- the seed's phase is part of the basin.

### `_period_grid`

The docstring paragraphs before the move:

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

Continuing after the fractions paragraph:

Measured on van der Pol at mu=100 (`benchmarks/pss_lte_grid.py`):
1105 steps taken from an adaptive transient converge where 1105
UNIFORM steps do not, and beat a 20000-point uniform grid on
accuracy -- 18x fewer points and -47.3 ppm against -60.6.

The comment on the opening step:

⚠ THE OPENING STEP IS MANUFACTURED, SO IT MUST NOT BE THE GRID'S
COARSE END.  `_traverse` builds `x(0)` from the unknown with ONE
order-dropped step of `hs[0]`, and a grid taken from an adaptive
transient opens wherever that transient's window happened to
start -- on van der Pol at mu=100, `h[0] = 1.4845` against a
MEDIAN of 4.62e-04, 3200x coarser.  That single Euler step moves
the state 7.4%, and the shooting Newton then has to invert a map
whose first act is that step.  Opening at the grid's own finest
step instead costs ONE extra step in 1105 and is the difference
between not converging and converging.

⚠ THE 8x IS MEASURED, NOT DERIVED, and it is a guard rather than
a threshold: it exists so grids that already work are left
exactly as the caller wrote them ('2:1' opens at 2x its finest,
'smooth' at 5x, and neither needs this).  The falsifier is in
the record: on van der Pol's grid, an opening step of 1e-1 still
fails for gear and 1e-2 converges, against a ratio here of
13939.  Anything between those bounds separates the two cases.
⚠ AND IT IS ONLY NEEDED WHEN THERE IS A MANUFACTURED STEP TO
PROTECT.  The subdivision exists because `_traverse` builds `x(0)`
with one order-dropped Euler step of `hs[0]` FROM `x_in` -- an
iterate that may be far from the orbit -- and a coarse `hs[0]`
there defeats the inner Newton.  With `x0_unknown` the first step
starts ON the orbit and the same coarse step is solvable:
measured on van der Pol's own LTE grid, the raw 1105-step grid
converges and reaches -47.3 ppm where the subdivided 1106-step
one reaches -73.8.  So the subdivision COSTS accuracy, and it is
skipped where it buys nothing.
⚠ A RAMP, NOT ONE TINY STEP (2026-09-21).  The subdivision used to
put a single step of `fr.min()` in front of the coarse first step:
on a relaxation van der Pol's own `lte_grid` (mu = 10) a 2.6e-4
step before a 3.5e-2 one, a 138x GROWTH beyond a two-step method's
zero-stability bound, so the integrator dropped two steps to Euler
and gear's transposed replay refused the grid gear had produced.
Skipping the subdivision for the solved-history kind was tried
first and the coarse first step then defeated the per-step Newton
(the entering history is flat at the seed), so the protection is
needed there too.  Doubling steps `d, 2d, 4d, ...` up to the first
step keep every ratio at 2, inside the bound, for ~log2 of the
span in extra points.

The comment on the ramp's construction:

⚠ TOP-DOWN (2026-09-21): built bottom-up as d, 2d, 4d, ..., rest,
the remainder could be a sliver before the next coarse step --
a growth beyond the bound, an Euler drop, gear +918 ppm on its
own folded grid.  Halving from the coarse step down to its
finest scale, with the smallest piece doubled up, keeps EVERY
ratio at 2 or less, the hand-off to the next step included.

The comment on the step-ratio warning:

⚠ A CALLER'S GRID CAN SILENTLY DEMOTE GEAR-2 TO FIRST ORDER.
`_period_grid` validated positivity and sum-to-1 and nothing
about the INTERIOR ratios.  A two-step method is zero-stable only
up to `h_n / h_{n-1} = 1 + sqrt(2)`, and past it the integrator's
own guard drops the step to Euler -- correct, and invisible.
Measured on a Q=20 resonator driven at resonance with an
alternating 3:1 grid, where half the steps are demoted:

      npts   uniform    3:1 grid
       100   19.91489    7.99821
       200   20.00960   11.42923
       400   20.02218   14.54985
       800   20.02443   16.85280

against an analytic peak of 20 V -- 60% low at 100 points,
crawling up at FIRST order, and `converged = True` every time.
Refining does not fix it, because refining a 3:1 grid keeps it
3:1.  So the warning names the ratio rather than suggesting a
smaller step, which is the advice that does not work here.

⚠ This is what item 5 removed the premise for: the literature
note in the class docstring argues Wambacq's objections to
non-uniform BDF "do not bite inside a run" because the grid is
UNIFORM and frozen.  A caller's grid is frozen but not uniform.
⚠ ONLY REPEATED UP-STEPS COMPOUND (2026-09-20).  An ISOLATED
up-step -- an event ramp, an inserted event -- is harmless at any
ratio: the factor w/2 the recursion applies acts on the difference
across the SMALL step, and their product is (h/2) x', the
trapezoidal predictor.  Measured: gear across a 1e-18 event ramp
(ratio 2.5e9) is second order, and gear on an event grid with
ratios up to 20 at N = 50 is 4x BETTER than uniform.  The 60 %-low
alternating 3:1 grid has a bad ratio every other step, and it is
that repetition this warns about, so a bad ratio counts only when
another lies within `RATIO_ISOLATION` steps of it.  ⚠ This
traversal never drops a step to Euler -- what the text used to say
is `Transient`'s `check_order_drop`, which does not run here.

The comment on the closing convention:

⚠ THE 'CLOSING' CONVENTION (2026-09-21, B7c completed): on a
caller's grid the inner steps keep the ABSOLUTE lengths they had
at the first call of this solve (the seed period) and only the
last step moves with `T`.  Under 'proportional' every step is
rescaled with each trial period, so an adaptive grid's fine
regions slide off the edges they were placed on while the outer
Newton hunts for T -- measured on a relaxation van der Pol
(mu = 10) on its own `lte_grid`: with the period seeded 16 % low,
gear's and trbdf2's per-step Newton FAIL on the coarse steps that
then land on the edges.  The closing step is bounded below so a
trial period that would swallow it falls back to proportional
scaling for that evaluation, with a one-time warning.

The ratio warning's repeat rule was "another bad ratio within
`RATIO_ISOLATION` steps" until 2026-09-24.  A landed switching window's
exit -- its ramp, then the partial step back to the base grid -- is two
adjacent bad ratios, so it counted as repeated: a false alarm under gear on
the comparator oscillator at 100 and 800 points.  Measured: smoothing those
pairs into doubling ramps never improved the period, and moved it by more
than the pairs could have cost -- -6.8e-4 to +3.3e-3 (100 points), +8.7e-5
to -7.4e-4 (150), +1.1e-4 to -5.8e-4 (200), +1.3e-5 to -2.9e-5 (800).  Now
two others must lie within the isolation distance; an alternating 3:1 grid
has them at every bad ratio.

### `_period_quadrature`

The docstring before the move (its first line, then the paragraph):

Trapezoid weights `W_n = (h_{n-1} + h_n)/(2T)` for samples at

Continuing from the middle of the paragraph:

non-uniform grid (`solve(grid=...)`, `lte_grid`) `1/N`
is not a quadrature at all: measured 2026-09-19, a flat 59 % error in
the first harmonic and pnoise converging 30 % off, rate 1.00.
Measured with these weights on the same 3:1 grid: 1.8e-04 / 4.4e-05 /
1.1e-05 at N = 200/400/800 -- second order, the rate of the gear
samples themselves (a rectangle `h_n/T` is first order and becomes
the bottleneck for every method).  ⚠ Still four orders behind a
uniform grid on that circuit, and second order caps radau: a derived
grid is now CORRECT for noise analysis, not competitive.  ⚠ The SAME

The comment's first line:

⚠ THE SECOND-ORDER CAP IS LIFTED (2026-09-21): a periodic cubic

### `_replay_grid`

The docstring paragraph before the move:

⚠ UNTIL 2026-09-20 THESE REPLAYS WERE ALWAYS UNIFORM.  `factored_period`
passed `len(times) - 1` and the builders built `linspace(0, T, npts+1)`
whatever grid the solve had run on, so the PPV, the Floquet modes and
every noise fold of a radau / trbdf2 / esdirk43 run -- and of trap and
euler through the TR-BDF2 twin -- on a non-uniform grid were computed
on a uniform replay from the converged `x0`.  Accurate (radau is
order 5 on most grids), and the reason radau read as "exact on the
3:1 grid": its adjoint never saw that grid.  Only the solved-history
(gear) kind replayed on the caller's grid.  Now `factored_period`
hands the solved fractions down and the replay is on the grid the
solve was on; a direct call with a bare `npts` is uniform as before.


## `_pss_events.py` -- `_StateEvents`

### `_state_event_rows` (the comment above it)

⚠ THE STATE-EVENT STAGE (2026-09-21) -- see `solve`'s docstring.

### `_land_fractions`

The comment on the window sub-grid, as it stood before the history moved
out (2026-09-24):

⚠ A WINDOW BETWEEN TWO EVENTS GETS ITS OWN SUB-GRID (2026-09-21).
A threshold switch declares both edges of its transition; the
segment between them holds the whole S-curve of the switch, and
as ONE step it is integrated by three collocation points across
the curve -- the stage solve stalled at 6e-4 of the swing
whatever the count.  A gap between two landed events narrower
than its neighbours is split into `EVENT_WINDOW_STEPS` steps,
which the remap then scales with the window.
(against the BASE grid's local step, not the immediate
neighbours: an inserted event leaves a sliver beside the window,
and measured against that the rule never fired on the PWM
loop's on-window)

### `_stage_one_crossings`

The docstring's last lines, as they stood before the history moved out:

the orbit crosses none.  (Refactor E9 item 2: the three stages each
carried this.)

### `_finish_state_events`

⚠ THE GRID EVERY CONSUMER REPLAYS ON IS THE ONE `_period_grid` MAKES
OF THESE FRACTIONS -- with its opener ramp, which the window
sub-grid's tiny steps trigger.  The stored pieces were first
computed on the unramped fractions, and `factored_period`'s node
indices were shifted by the ramp's pieces against them: the
bordered sideband solve read its event rows at the wrong nodes
(dtheta 12 % off, the response 51 %).  So the ramped grid IS the
grid from here on: `_grid_fracs` carries it (its first piece is the
smallest, so `_period_grid` will not ramp it again), and the
monodromy and the event columns are computed on it -- the identity
remap on it gives every step's sensitivity to the events (the
ramp's pieces included) and the event nodes.

⚠ THE TOTAL MONODROMY THROUGH A MOVING EVENT (2026-09-22, phase B).
The period map's derivative is not `M` (the grid frozen): a
perturbation of `x_0` moves the crossing, `dtheta/dx_0 = -Gt^-1 G`
from the event rows, and the state at the period moves with it
through the event columns -- the bordered system's Schur
complement, which is the saltation matrix derived rather than
guessed.  Measured on the PWM loop: `M` 109 % off the finite
difference of the staged period map, this 3.3e-8; the dominant
multiplier 0.691 where `M` read 0.632.

The one write of `_monodromy` (2026-09-24).  Until then the stage kind
wrote a PARTIAL map from inside the stage's Newton (`evmap`, a remnant of
the `_traverse_stage` view) and gear's pair wrote none, while this method
overwrote both with the total map when the event columns were built.  So
the normal path was consistent -- measured on the PWM loop, both kinds
reported the total map (radau 6e-11, gear 5e-13 from `M + P_end dtheta`)
-- and a comment in `evmap` saying "a partial map either way" was wrong.
The failure path was not: with the column assembly failing, radau reported
the grid-frozen map on the landed grid (0.698) and gear's pair STAGE 1's
map, from another grid and another orbit (0.905 against the landed map's
0.730).  Now this method writes it once, for every kind.

### `_state_event_stage`

The docstring's opening, as it stood before the history moved out:

The bordered second stage: the crossings of the first stage's
orbit become Newton unknowns (events phases A/B).  One stage for every
kind that has one since 2026-09-23 -- it was three, the driven and
the autonomous stage methods and gear's pair.

The comment above the `_finish_state_events` call:

(gear's pair builds its columns whether or not the stage converged,
as it did before the stages shared this finish)

### `_event_costate_injection`

The event nodes' costate injections for a reverse pass of the
TOTAL map (2026-09-22): `-zeta_k W_k` at node `nd_k`, `zeta = Gt^-T
P_theta^T v` -- the transpose of the saltation, carried by the pass
to every earlier node; `None` when the solve is not staged.  What
`ppv()` and `floquet_modes` sample a left vector along the orbit
with.


## `_pss_inner.py` -- `_InnerTransient`

### `_factorise`

2026-09-25: every stored step now factors through here -- Radau's coupled
block, the GLM startup and gear's continuous adjoint had called
`lu_factor` directly and ignored `linearsolver=`.  `_PerSolve` gained a
transposed solve (with none, a solver lacking `factor` crashed every
adjoint surface), and KLU/Auto gained `factor` (log, 2026-09-25).


The docstring's opening paragraphs before the move:

One step's `Jf`, factored by the CALLER'S linear solver.

⚠ THIS USED TO REACH FOR `scipy.linalg.lu_factor` DIRECTLY, which
made every matrix-free run dense-LAPACK whatever `linearsolver=`
said -- and a circuit Jacobian at m=1002 is very sparse, so a dense
LU is ~3e8 flops per step, `N` times over.  The recorded 2.13x at
m=1002 and the m~250 gate were therefore both measured against a
DENSE baseline; see `benchmarks/pss_matrix_free_ceiling.py` for what
they become when both sides get a sparse solver.

### `_k_at`

The docstring before the move:

The reduced STAGE DERIVATIVE `dq/dt = -(i(x) + u(t))` at a point --
what the DIRK and coupled-Radau period columns need, and with `t`
the stage time what a DRIVEN circuit's event columns need
(2026-09-21): at t = 0 the source is wrong on every other stage of
a driven circuit, and on the PWM fixture's ramp row that read as a
-3.9 in the event row's derivative where the FD said +4.2.

⚠ THIS USED TO BE `-i(x)` ALONE, on the reasoning that an autonomous
circuit has "no source term" (2026-09-20).  An autonomous circuit
has a CONSTANT source vector, not a zero one: a DC supply, a bias
current.  On a row a source pins, `i(x) = -u` at convergence, so
the true derivative is 0 and `-i(x)` is `u` -- and the period
column read `-u/T` there.  Measured on the tree's own phase
fixture (a 1 kV DC supply): radau's column `[-1e6, ~0, -1.8e-4,
...]` against the finite difference `[0, -3.18, 6.28e3, ...]`,
and on a van der Pol with a decoupled 1 kV node the oscillator rows
agreed with the FD to 4 digits while the supply node read
-157.9 = -1e3/T and its branch current +i_R/T.  The first Newton
step on that column threw `x0` to 2.8e6 and the stage matrix went
singular there, mislabelled "seed below the fundamental".  Both FD
checks in the build were on a source-free van der Pol, where the
two expressions coincide.  `u` is taken at t = 0 because this is
only built for autonomous circuits, where it is constant.

### `_install_history`

The comment on the `b != 0` refusal:

⚠ A `b != 0` COMPANION IS REFUSED HERE, AND THE FIRST REASON
GIVEN FOR IT WAS WRONG.  It said a solved history carries CHARGES
while such a method also reads `iq_{-1}`, "which no charge
determines".  The DAE determines it exactly -- a converged point
satisfies `i(x) + iq + u = 0`, so `iq_{-1} = -(i(x_{-1}) + u)`,
the same identity item 4d rests on -- and seeding it was tried.

It fails for the derivative running the OTHER way.  A one-step
companion reads only `iq_{-1}`, so the trajectory depends on
`x_{-1}` solely through it, and `d(iq_{-1})/d x_{-1} = -G` is
SINGULAR wherever a node carries no conductance -- every purely
reactive node, which is most of a resonator.  Admitting `x_{-1}`
as m unknowns then leaves the 2m x 2m system rank-deficient:
measured, `LinAlgError: Singular matrix` on 25 tests at once.

The right second unknown for such a method is `iq_{-1}` ITSELF --
the `(x, iq)` state its monodromy already uses -- closed by
`iq_{-1} = iq_{N-1}`.  That is a different formulation, not a
seeding fix, and it is not built.

The comment on `_dt_last`:

⚠ THE STEP THAT PRODUCED `x_0` IS THE PERIOD'S LAST ONE, NOT ITS
FIRST.  `x_{-1}` sits one step BEFORE `x_0`, and on a periodic
grid that step is `hs[-1]`.  With a uniform grid the two are
equal and this never showed; on a caller's grid (item 5) with a
16438:1 spread, handing `hs[0]` to a method that reads `h_last`
states a step ratio that never happened.

### `_sync_limit_at`

The docstring's middle paragraphs before the move:

(`x_n` and every stage).  Whatever the last step's solve happened to
leave behind (the LAST stage) was therefore used for all of them, so the
period map linearised the junction at the wrong voltage.

Measured against a finite-difference derivative of the discrete period
map (a reference this code cannot influence) on a diode loaded through a
series resistor: the analytic monodromy was off by a FIXED 1.65e-3
(Radau) / 9.1e-4 (TR-BDF2) relative -- flat across four decades of the
FD step, so a genuine error and not FD noise -- and the error grew with
how hard the junction was driven, vanishing when it was off.  With this
sync the same comparison lands at ~3e-9, the FD noise floor.

### `_C_at`

The docstring before the move:

The reduced capacitance at a point, without taking a step.

⚠ NO LIMITING SYNC HERE, AND THAT IS MEASURED, NOT ASSUMED.  `_G_at`
needs the device limiting state to be at the point it is evaluating,
because a junction's `i`/`G` are read at the stored `_vlim`.  CHARGE IS
NOT: surveyed across every limiter in the tree,

* `elements.Diode` is the only STATEFUL one (it keeps `_vlim`), and
    its `C`/`q` do not read it -- with the stored state moved far from
    the evaluation point, `dC = dq = 0` while the control `dG = 15.2`
    and `di = 3.9e-1` confirm the limiting was live;
* `Semiconductor` (BJT/JFET/ZenerDiode/Varactor) limits STATE-FREE by
    construction -- "Return a limited copy of `x` -- STATE-FREE, and
    that is the point";
* `compact.PspMosLongChannel` likewise returns a limited copy;
* the hdl devices keep no `_vlim` at all (it is a codegen local).

So there is no device whose capacitance a sync could correct.  A sync
was carried here for a while as "correct in principle" insurance and was
never exercised by any test -- this tree's own rule is that unexercised
machinery is a liability.  ⚠ If a stateful limiter whose CHARGE reads its
state is ever added, this is where the sync goes back; `_sync_limit_at`
is kept for that, and for `_G_at`'s no-junction path.

### `_G_at`

The docstring's middle paragraphs before the move:

⚠ THIS GOES THROUGH PCNR WHEN THERE ARE JUNCTIONS, and the reason is
STATELESSNESS, not accuracy.  `pcnr.augmented_system` + `schur_reduce`
build `G` from an explicitly-passed `v_lim` instead of from the device's
stored one, so the answer depends on the POINT ALONE.  The `limit(x, x)`
route `_C_at` still uses does not: `limit` clamps relative to the STORED
`_vlim`, so it lands on the true point only when the previous evaluation
was already nearby.  Measured, varying the prior `_vlim` before
evaluating at a fixed point: PCNR's `J_eff` moves by 0.0, the limit-sync
`G` by up to 15.15.  It was right in the traversal only BY LOCALITY
(steps are small, so the prior state is always close) -- the same
accident `_begin_period` warns about when it insists the period map be a
function of `x0` alone, applied to its linearisation.

Numerically this changes NOTHING today: against a finite difference of
the discrete period map both routes give the same monodromy to every
printed digit (3.025e-09 radau / 2.358e-09 trbdf2, identical either
way).  It removes a latent order-dependence, and it is what lets the
transient and the monodromy share ONE limiting.

### `_pq_seed_at_x0`

The docstring before the move:

``d(iq_0)/d(x_0)`` when the method SEEDS a consistent companion current.

⚠⚠ THE CHAIN RULE THE `open_at_x0` PATH ASSUMED AWAY.  Every branch
that opens at `x_0` seeds `Pq = 0` and says so in the same words --
"no companion current has been formed yet".  That was true of every
method in this tree until `theta`, which refuses the L-stable opener
and therefore READS `iq_{-1}` on its first step: `_begin_run` seeds it
at ``iq_0 = -(i(x_0) + u(t_0))``, the DAE's own `dq/dt`, and that is a
FUNCTION OF `x_0`.  Differentiating it gives `-G(x_0)`, and dropping
that term is not a small error -- it is the whole `null(C)` mode.

Measured on the B2 gate resonator at `K = 200`, against a
finite-difference of the shooting residual (delta-swept over six
decades, FLAT, so a real error and not FD noise): the analytic
monodromy mapped `null(C)` to ZERO -- exactly what an L-stable Euler
opener would do -- where the true map multiplies it by `-0.7778`,
which is `(-(1-theta)/theta)^K` from `ThetaIntegrator`'s own table.
Relative Jacobian error 6.344; with this seed, 1.4e-10.

The cost of that was NOT a wrong answer -- the residual is what it is,
so the solve still lands on the right orbit -- but the Newton lost its
quadratic step: on a LINEAR circuit an exact shooting Newton converges
in ONE iteration (`trap` with `x0_unknown` takes 3 evaluations at every
`K`), and `theta` was taking 9 / 64 / 99 at `K = 100 / 200 / 400`.

`None` -- the default for every method that does NOT declare
`needs_consistent_iq0` -- means the zero seed is exact, and those
methods stay bit-identical.

### `_transient`

The docstring before the move:

The `Transient` this analysis integrates with.

PSS used to carry its OWN transcription of one integrator step --
the third in the tree, after `Transient` and `JAXTransient` -- and
it had already cost the two defects its docstring records: `method`
declared and never read, and a companion current fed back from the
iterate before the converged one.  Driving the real thing removes
the copy and brings what came with it: the limiting machinery, PCNR,
breakpoint order drops, and the continuation rescue.

⚠ THAT LIST WAS 1-FOR-4 AS SHIPPED (external review, 2026-09-02); it is
now 2-FOR-4.  LIMITING reaches -- `cir.limit` is called on the inner
Newton and the rectifier measurably conducts.  PCNR now reaches too:
`PSS(cir, pcnr=True)` is a declared Parameter forwarded to the inner
`Transient` above, and PCNR lives in `Transient.solve_timestep` (the
LMM `_solve_timestep_pcnr` and, for stage methods, `_rk_stage_pcnr`),
which PSS DOES call -- so it needs no per-accepted-step machinery.  The
remaining two still do not: the continuation rescue (`_rescue_solver`)
and breakpoints (`cir.next_event`) are armed only inside
`Transient.solve`, which PSS never calls -- it drives `solve_timestep`
directly on its own frozen grid, so a breakpoint has nothing to move.
The same structural fact behind the TLine refusal above: what
`Transient.solve` does per accepted step, PSS does not do at all.

The tolerances are handed over unchanged, which is the point of
`newton_tolerance_vectors`: `reltol`/`iabstol`/`vabstol` mean the
same thing on both sides, so passing them through is a no-op in
meaning.

The comment on the integrator mapping:

⚠ A MAPPING, not an if/else on 'euler'.  Written as
`EulerIntegrator() if method == 'euler' else Trapezoidal...`
it silently ran trapezoidal for every other name -- caught
while adding 'gear', which produced numbers identical to
trap's to the last digit.  This class has already paid once
for a `method` that selected nothing; a dict raises KeyError
on a name nobody wired.

### `_theta_biased`

The docstring before the move:

Give a `ThetaIntegrator` the bias THIS period needs, not a fixture's.

⚠⚠ `theta - 1/2 = C h` makes `C` a RATE, but the quantity that decides
anything is the DIMENSIONLESS product `C T`: `null(C)` is damped over a
period by `((1-theta)/theta)^K ~ exp(-4 C h K) = exp(-4 C T)`, and `h`
cancels.  So `ThetaIntegrator.DEFAULT_C = 1e4` is not a recipe -- it is
the measured knee `C T = 0.0628` divided by ONE fixture's period
(6.283e-6 s).  On a circuit 159x slower it is 159x the calibrated bias,
and MEASURED on `_q20_rlc` (analytic 20 V) that is a peak of 15.91 at
K = 100 -- 20% low, `converged=True`, because it did converge: to its
own over-damped discretisation.  See `ThetaIntegrator.DEFAULT_CT`.

This is where a shooting run stops inheriting that.  `theta_ct` is the
dimensionless knob (`None` = the measured knee) and the period is this
solve's, so `method='theta'` is now correct on any circuit.

⚠ THE SEED PERIOD IS ENOUGH, and that is a measurement not a hope: the
gate's own table has `rcond(I - A^K)` at 4.4e-03 / 4.6e-03 / 4.3e-03
across `C T` = 0.0063 / 0.0628 / 0.628, i.e. FLAT over two decades.  An
autonomous solve moving `T` by a few percent moves the bias by the same
few percent, which the knee does not notice.

Every other method is returned untouched, so nothing else moves a bit.

### `_new_transient`

The docstring's last paragraph before the move:

⚠ THE SOLVER STRATEGIES GO THROUGH TOO, and they used not to.
`nrsolver`, `linearsolver` and `scaler` are declared on the base
`Analysis`, so `PSS(cir, linearsolver=...)` has always been ACCEPTED
-- and then dropped here, with the inner `Transient` resolving to
`DenseSolver`/`StandardNewton` whatever the caller asked for.  That
was the third time this class took a parameter it never read
(`method` declared and never read; `analysis='PSS'` matching
nothing), and the same shape each time: accepted at the constructor,
silently discarded at the boundary.

The comment on `_damped_last_resort`:

The line search as the last resort on the shooting path, which
never arms the transient's rescue ladder (owner decision
2026-09-08, "Do 2"; see `_rk_step_coupled` and `solve_timestep`).

### `solve_timestep`

The docstring's opening paragraphs before the move:

One timestep of the inner transient, taken by `Transient`.

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

The docstring's last paragraph before the move:

The measured cost of backward Euler on a limit cycle is unchanged
and still the reason `method` matters -- it damps exactly what PSS
exists to find:

    steps/period    PSS peak    fraction of analytic
              20      2.63 V       13.2%
              50      5.61 V       28.1%
             100      8.81 V       44.1%
             200     12.20 V       61.0%

The comment on the step:

ONE INTEGRATOR STEP, taken by the class that owns the definition.
`Transient.solve_timestep` applies the chosen integrator through
`get_diff` (so `method` selects something because the integrator
object does), the limiting machinery, PCNR when asked for, and the
continuation rescue.  None of that existed on the copy this
replaced.

The comment on the event columns:

the event columns of a two-step companion need BOTH partials
(2026-09-22): its own step's and the previous step's, the
latter assembled from the total under uniform scaling

The comment on `residual_dT`:

⚠ `residual_dT`, NOT `residual_dh`.  The grid is rebuilt at
the current `T` on every residual evaluation, so every step
scales and the partial `d/dh_n` is not the total -- for
Gear-2 it is 3/2 of it on a uniform grid, measured against
finite differences at 1.4859/1.4939/1.4972 for 100/200/400
points, converging on the exact 3/2.  Euler and trapezoidal
were never wrong: their coefficients depend on `h_n` alone,
so the partial IS the total, which is why only Gear-2 was
hit.  See `Integrator.companion_dT`.

The comment on seam steps, its opening paragraphs:

A SEAM STEP IS ONE WHOSE COMPANION READS THE ENTERING
UNKNOWN, not merely one whose ESTIMATOR does.

⚠ THIS CONDITION WAS `h_last2 is None` ALONE, AND THAT FLAGGED
A PHANTOM FOR TWO METHODS OF THREE.  `h_last2 is None` is the
transient's statement that the third past charge is not real,
which is the reach of the LTE estimator's third divided
difference -- not the reach of the integrator.  Euler's
companion reads `q_{n-1}`; trapezoidal's reads `q_{n-1}` and
`iq_{n-1}`, and the order-dropped opening step supplies an
`iq` consistent with it, which is what that drop is FOR.
Neither can see the fabricated charge at all.  Measured
(`benchmarks/pss_seam_cost.py`): trapezoidal's seam reading
was 15.1 times tolerance while its cost is 1.3e-11 V, and
euler's 0.286 against 5.1e-12 V.  Both are exactly zero; the
reading was an artefact of the measurement.

Gear-2 reads `q_{n-2}` -- which at that step IS the entering
unknown -- and the shooting condition constrains `x(0)` to
equal `x(P)`, NOT `x_in` to be the orbit's own `x(-dt)`.  So
`x_in` sits O(h^2) off a real history point and Gear-2 reads
it as one.  That one costs 1.266e-01 V at 100 points/period
against an interior contribution of 1.070e-01 -- the seam is
54% of its total error -- and it is the term that STOPS
converging: it falls as h^2 while the interior falls faster,
so its share grows to 68% at 200 points and 73% at 400.

TWO DIFFERENT THINGS ARE TRUE OF AN OPENING STEP, and the
first version of this conflated them -- which showed up as
trapezoidal's phantom simply MOVING from the seam into the
interior total (0.340 -> 15.47) when the seam test was
tightened.  Suppressing a bad number is not the same as
classifying it.

Its seam paragraph:

a SEAM exists when the COMPANION reads the entering
unknown, i.e. when its charge reach `len(alphas) - 1` is
deep enough to touch it.  Only Gear-2's is.

The comment on `_history_is_solved`:

`_history_is_solved` is that formulation saying the
deepest charge is an UNKNOWN the solve closed, not a stand-in
-- so there is no seam to report even though the companion
reaches that far.  Without this the fix would go on flagging
the defect it removed.


## `_pss_walks.py` -- `_PeriodWalks`

### `_step_sensitivity`

The docstring's lines before the move:

ONE RECURSION FOR EVERY METHOD (and now for either formulation).

Its paragraphs on the solve:

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

The comments on `coeffs` and `source`:

`coeffs` overrides the LIVE `_coeffs` for the same reason `solve`
overrides the solve: a matrix-free replay happens after the run
that produced the steps, when `_coeffs` no longer describes the
step being replayed.  See `_traverse_factored_plain`.
⚠ `source` ENTERS THE SOLVE AND NOT `Pq`, and the asymmetry is
the physics rather than a convenience.  A small-signal source
appears in the step's residual -- `Jf dx + S + du = 0` -- but NOT
in the companion, because `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`
is built from CHARGES, and an injected current is not one.
Adding it to `Pq` as well would feed a fictitious charge forward
into every later step, and the error would grow along the period
rather than announce itself.

This is what makes PAC share the recursion instead of copying it:
the homogeneous propagation (`source=None`) is the monodromy and
the driven one is the forced response, and they differ by this
one term.

The comment on the default solve:

⚠ THROUGH THE CALLER'S SOLVER, not `toolkit.linearsolver`.
This is the DENSE propagation -- the thing matrix-free is
measured against -- and it used to be hardcoded to the
toolkit, so `linearsolver=SuperLUSolver()` reached the inner
Newton (once forwarded) and never the propagation.  Comparing
a sparse matrix-free path against a dense baseline would have
flattered it; both sides go through the same strategy now.
`DenseSolver` IS `toolkit.linearsolver`, so the default is
unchanged.

### `_walk_lmm`

The docstring's opening paragraph before the move:

ONE WALK OF THE PERIOD UNDER A LINEAR-MULTISTEP COMPANION -- the
plain map and gear's solved-history pair, dense or factored
(2026-09-23: it was four walks, `_traverse`,
`_traverse_solved_history`, `_traverse_factored_plain` and
`_traverse_factored`, each with its own copy of the period column).

Its second paragraph:

EVERY SHOOTING ITERATION IS ITS OWN RUN.  phi must be a function of
its arguments alone; if iteration k+1 inherited the ring buffers
iteration k ended with, the period map would depend on which
iteration it was and the monodromy would be the derivative of
something else.  `_begin_period` / `_install_history` make the run
start fresh.

Its last paragraph:

⚠ ONE LINEAR SOLVE PER STEP, CHOSEN BY WHETHER THE STEP IS KEPT.  A
kept step is factored once (`_factorise`) and every column solves
against that; a dropped one solves through the caller's solver
directly.  Both are an LU with partial pivoting, but numpy's `solve`
and scipy's `lu_factor`/`lu_solve` differ in the last bits (measured
on this box: 140 of 360 random systems, up to 3.5e-13 relative), so
each walk keeps the arithmetic it had -- the merge moved no answer.

Returns a `_PeriodWalk`; the `_traverse*` names are its views.

The comment on the installed history:

THE HISTORY IS INSTALLED, NOT SEEDED.  `_begin_run(x_{-1})`
opens the rings on the earlier point and the push puts `x_0` in
front of it, so the first real step reads `q(x_0)` and
`q(x_{-1})` -- two genuine solved points.  The flags then say
what is true of them: a step of `dt` has been taken
(`_dt_last`), the run is no longer opening (`_is_first_step`,
`_no_history`), and `_dt_last2` stays None because the THIRD
charge is still `q(x_{-1})` repeated -- the LTE estimator
differences three, so its opening reading remains unsound and
the report goes on discarding it.

The comment on the opening's coefficients:

⚠ THE OPENING'S COEFFICIENTS ARE THE OPENING'S.  `_coeffs` is
live state and the manufacturing step is order-dropped, so it
reports Euler's `(alphas, b)` -- `b = 0` -- where the loop's
steps report the method's own (trapezoidal opens at
`((49000, -49000), 0.0)` and runs at `((98000, -98000),
-1.0)`).  Reading them once for the whole run put the period
column 40-50 % out for trap and gear; reading the loop's for
the opening made `Pq` non-zero where it is zero, a 100 %
error for trap.  Euler was exact both ways and would have
passed a one-method test.  ⚠ `_coeffs` DOES NOT EXIST YET when
opening AT `x_0` -- no step has run -- so it is not read:
`b_open = 0`, no companion current has been formed.

The comment on the event columns, its first line:

⚠ EVENT COLUMNS FOR A TWO-STEP COMPANION (2026-09-22, phase B).

The comment on the per-step coefficients, its end:

describes the step being replayed.  (Inside the loop they are
constant for every method in this tree, so storing them is
belt-and-braces today -- a mutation replacing them with a
post-run snapshot does NOT fail the tests -- but a
variable-order method would make that failure silent.)

The comment on `residual_dh`, its end:

shift of the step's START.  Counting `u_dot (tau + w)` on
top of `residual_dh` -- twice -- read the node after a
landed event 153 % off (2026-09-22).

### `_traverse_solved_history`

The comment on the returned map:

⚠ ALWAYS THE FULL 2m x 2m MAP, NEVER THE `d x_{N-1}/d x_0`
CORNER.  This used to hand the corner back on the driven path,
and a sub-block of a sensitivity is not a monodromy: it reported
`spectral_radius` 1.279605 for the Q=20 resonator -- ABOVE ONE --
where the analytic per-period decay is exp(-pi/Q) = 0.854636.
For a two-step method the one-period map acts on the PAIR, and
its spectrum carries the discretisation's parasitic roots beside
the physical multipliers; BDF-2's is 1/3 per STEP, (1/3)^N over a
period, and `_spectral_report` separates them by eigenvector
block structure anyway.

### `_walk_glm`

The docstring's last lines before the move:

what the approximation can cost is Newton iterations.  MEASURED on the
index-2 C-V loop: 3 iterations to 1e-12, the same count as radau's
exact monodromy on the same fixture (roadmap).

Until 2026-09-24 the walk seeded the recursion with the startup's first
component only.  The docstring said:

⚠ THE JACOBIAN IS APPROXIMATE BY CONSTRUCTION and the residual is not:
only ``dQ_0/dx_0 = C(x_0)`` is carried into the recursion, the startup's
dependence of the higher Nordsieck components on ``x_0`` (p Radau
substeps and an interpolant) is dropped.  So the converged fixed point is
the method's own, exactly; what the approximation can cost is Newton
iterations.

and the period column's comment:

The startup's own T-dependence enters twice: through the SCALING
`Q_k = h^k q^(k)` (kept -- `dQ_k/dT = (k/T) Q_k`) and through the Radau
substeps at `h/p` (dropped, as `Mx`'s is).

What the approximation cost, measured when `_GLMStartup` replaced it: the
map 3.8e-2 / 9.5e-2 from finite differences (glm2 / glm3, van der Pol),
the period column 6.9e-6 / 3.4e-10; a LINEAR driven RLC took 9 / 14
Newton evaluations; matrix-free, whose operator is the factored map,
DIVERGED on that RLC; and `solve` kept no GLM monodromy, so no spectrum.
With the startup linearised: 3e-9 / 2e-10, 8e-12 / 6e-11, one Newton
step, matrix-free equal to dense to 1e-15.

### `_walk_stage`

The docstring's opening lines before the move:

ONE WALK OF THE PERIOD UNDER A RUNGE-KUTTA STAGE METHOD (Radau
IIA, TR-BDF2, ESDIRK), dense or factored (2026-09-23: it was two
walks, `_traverse_stage` and `_traverse_factored_stage`).

Its paragraphs on the finite-difference check and the event columns:

Finite-difference checked before use (the dT column has been got
wrong in this file twice -- roadmap 0j).

⚠ EVENT COLUMNS (2026-09-21): `hsens[j, k] = d h_j / d theta_k` for
the state-event unknowns, propagated by the same stage algebra as the
period column with the per-step weight taken from the matrix instead
of `h/T`; `capture` names the nodes whose state and sensitivities the
bordered residual reads.  ⚠ A DRIVEN CIRCUIT'S SOURCES MOVE WITH THE
GRID: an event column shifts the TIMES the stages are evaluated at,
so `f = -(i + u(t))` changes by `-u_dot . dt_stage`, `dt_stage =
tau_n + c_i dh` with `tau_n` the shift of the step's start (the `U_i`
term).  Without it the FD check read the column 77 % off and the
event row's derivative with the WRONG SIGN on the PWM fixture.


## `_pss_replays.py` -- `_FactoredReplays`

### `_FactoredReplays` (the section comment above `_replay`)

The comment as it stood before the history moved out (2026-09-24):

THE REPLAYS -- one of each for every kind of factored period (plain
LMM, gear's pair, stage, GLM; 2026-09-23, they were a copy per kind).
A map is `extract . step_N ... step_1 . seed` (see `FactoredPeriod`);
the steps know their own algebra (`_LMMStep`, `_StageStep`,
`_GLMStep`).  The `_monodromy_matvec*` names below are the entry
points the matrix-free Newton builds on from raw step lists.

### `_replay`

⚠ COMPLEX `v` IS TWO REAL REPLAYS, NOT A COMPLEX FACTORISATION.  The
steps are real, so `M` is a REAL linear map and `M(a + ib) = Ma + i
Mb` exactly; PAC needs complex products (`I + alpha(f) H`), and
factoring in complex arithmetic would double the stored factors for
a map with no imaginary part.  ⚠ The float cast below is therefore a
GUARD, not a convenience: it used to swallow a complex `v` silently
by discarding the imaginary part, a wrong answer rather than an
error.

### `_replay_transposed`

⚠ THIS IS WHY IT COSTS NOTHING TO HAVE.  Demir & Roychowdhury (TCAD
22(2) 188-196) call reverse integration "often unavailable even in
existing time-domain simulators" -- true of a forward-only DENSE
implementation.  The factored period already stores every step's
factorisation, and every factorisation solves transposed, so the
reverse pass needs no new integrator and no second traversal
(measured against the dense `M^T` at 1.8e-15 when first built, 0.75x
the forward cost).  It is the shared dependency of the PPV, adjoint
noise and the sideband rows.

⚠ `collect` HANDS BACK THE PER-STEP COSTATES (`ts[j]`, the step's
transposed solve(s)) and the adjoint state after each step
(`states[j]`; for a PPV seed it IS the PPV there, `Phi(T,s_j)^T v(T)
= v(s_j)`), both in step order.  ⚠ `inject[j]` lands on the adjoint
state after `steps[j]` is applied backwards, so the output becomes a
functional over the whole period rather than a value at its end --
the difference between the response at `t = 0` and a sideband
coefficient.  ⚠ The collected lists may be NESTED (a stage method's
`ts` holds per-stage costates per step), so the complex split
recombines through `_cx_collect`: a flat `a + 1j*b` multiplied a
LIST by `1j` and `floquet_modes` under trbdf2 raised for as long as
that path existed.

### `_sideband_forced`

The forced (source-injected) part of sideband row `l`, and the
final costate `g` for the closure: the transposed replay with the
OUTPUT functional `d` injected at every node (weighted by
``exp(-j(l w0 + w) t_n)/N``, or the period quadrature's weight) and
the source coupling read at every step.  ⚠ The injection is added
AFTER the step's costate update, so the output at `t_n` couples to
the sources of steps `< n` (causality); the state at `t_n` is the
one step `n` enters from.  `extra` is a raw costate injection per
node (the bordered adjoint's event-row term, 2026-09-22), at `d`'s
position.  Returns `(forced, g)`.

### `_monodromy_matvec_transposed_plain`

`M^T v` for the PLAIN map from its raw steps (B8: derived, and
gated against the dense `M` built from the forward replay -- a
from-scratch adjoint in this file has come out sign-inverted before,
roadmap 0h).

### `factored_period`

The comment on the self-starting branch:

A self-starting method has its own factored map (no opener, no
pair), replayed on the SOLVED grid's fractions (None on a
uniform grid keeps the replay bit-identical to before) -- see
`_replay_grid` and `_factored_self_starting`

### `factored_period_stage`

Refactor E9 item 2 (2026-09-23): this was `factored_period_full` and
`factored_period_dirk`, routed by the caller.

### `_factored_self_starting`

The factored period of a SELF-STARTING method about `x0` -- 'stage'
(`_walk_stage`; a `FactoredPeriod` of kind 'full' or 'dirk') or 'glm'
(`_glm_period_blocks`; kind 'glm') -- integrated under `method` (or
the PSS's own) in a transient of its own, on `npts` steps: uniform,
or `grid`'s fractions (see `_replay_grid`).  One builder for both
(2026-09-23: `factored_period_stage` and `factored_period_glm` were a
copy each).


## `_pss_ppv.py` -- `_PPVFloquet`

### `PPV_DEFLATION_ITERS`

The comment before the move:

How many deflated power iterations estimate the second multiplier.
Convergence is at |lambda_3|/|lambda_2|, which is fast in the case
that matters -- a lone slow node leaves everything below it tiny.
⚠ RETIRED 2026-09-03, kept as a name so the history reads.  The
deflated power iteration this sized converged at
`|lambda_3|/|lambda_2|` and lost three digits at a ratio of 1.065;
Arnoldi replaced it at machine precision and fewer matvecs.

The constant itself was removed on 2026-09-24: it had no reader, and with
this document holding the history there was no longer a reason to keep the
name.

### `PPV_RITZ_BASIS`

The comment before the move:

Arnoldi basis size for the second-multiplier estimate.  Exact at
`k = n`; below that it is a truncation and `lam2` is a lower bound
ONLY for a normal `M` -- see `ppv`, where a circuit monodromy was
measured to break the bound in both directions.  This is now the
STARTING basis: `_ritz_second_multiplier` grows it until the pair's
own Ritz residual certifies it.

### `PPV_RITZ_RESIDUAL_TOL`

The comment before the move:

⚠ THE GATE ON A TRUNCATED `lam2`: the SELECTED PAIR's Ritz residual,
`|h_{k+1,k}| |y_i[last]| / ||y_i||`.  It needs no extra matvec -- both
factors are already in the `H` this class forms -- and it is the one
quantity that separates a converged pair from a leaked one.  The
GMRES-style residual cannot: `_arnoldi_gmres`'s own note says a
drifted basis "gives multipliers that are wrong in a way the residual
cannot see", which is true of the SOLVE residual and false of the
EIGENPAIR one.

MEASURED (`_osc_with_ladder`, k = 12): <= 3.1e-07 at every `nslow`
the truncated path gets right, and 2.8e-04 / 3.5e-04 / 2.1e-03 at the
three it gets wrong -- and 1.5e-16 once `k` is large enough to be
exact.  A peer session's independent sweep puts the two populations
thirteen decades apart at the median, ⚠ TOUCHING at ~1e-5 (right 90th
percentile 1.27e-05 against wrong 10th percentile 1.03e-05).  So the
robust band is BELOW ~1e-6, and that -- not a magic 1e-8 -- is what
this is set to.

### `PPV_RITZ_MAX_BASIS`

The comment before the move:

⚠ A COST CEILING, NOT A CORRECTNESS THRESHOLD.  `k` doubles until the
residual certifies or this is reached; overrunning it produces a
WARNING and an uncertified number, never a silently wrong one, which
is what makes an arbitrary-ish constant safe here.

The size is set by what `k` has to reach: measured, `k` tracks the
SLOW-MODE COUNT and not `n` (this tree's ladder cannot separate the
two -- it sets `nslow = nladder` -- but a peer's synthetic can, and
reports `k_min` flat under a doubling of `n` at fixed `nslow`).  The
largest published case is Lai's 64-gated-capacitor DCO, so 128 leaves
2x headroom on it while costing 1/6 of that circuit's `n = 813`.

### `_ritz_second_multiplier`

The docstring before the move:

`(lam2, residual)` from a `kk`-dimensional Arnoldi on `I - M`.

Garcia, Romero & Acha (IEEE Trans. Power Systems 37(1), 2022): the
Ritz values of `I - M` map back as `lambda = 1 - theta`.  Returns the
selected pair's own RITZ RESIDUAL alongside it,
`|h_{k+1,k}| |y_i[last]| / ||y_i||`, which costs nothing -- both
factors are already in `H` and its eigenvectors.

⚠ `eig`, NOT `eigvals`.  The eigenVECTOR's last component is half the
residual, so asking only for the values is what made this estimate
uncheckable for as long as it was.

⚠ EVERY LOCAL IS UNDERSCORED ON PURPOSE, inherited from when this was
inline in `ppv`: the first version used `q` for the Arnoldi start
vector, silently overwriting `C(0) xdot(0)` -- returned as
`info['q']` and consumed by two tests -- and the suite caught it as a
shape mismatch three frames away.

`residual` is `inf` when no Ritz value survives the deflation, so a
caller that gates on it cannot read "nothing found" as "certified".

### `_equation_row_ppv`

The docstring before the move:

`v_1` — the adjoint contracted with an EQUATION-ROW input.

⚠⚠ THERE ARE TWO ADJOINT VECTORS AND CONFLATING THEM WAS A DEFECT.
Demir 2000 puts both conventions on one page: a STATE initial
condition contracts as `v_1^T(0) C(0) x(0)` (eq 41, WITH `C`), while
an EQUATION-ROW input contracts as `v_1^T(s) b(s)` (eq 42, and the
phase equation 44, BARE).  `ppv()` returns `C^T v_1`, which is the
right object for a state perturbation and is documented as such.
`CY` is an equation-row covariance -- a current injected into a KCL
row -- so `diffusion_constant` and `colour_projection` need `v_1`,
and this produces it.

⚠ ONE DEFECT, TWO SYMPTOMS, both measured.  `(C^T v_1)_j` is COLUMN
`j` of `C` dotted with `v_1`, so on DIFFERENTIAL states `C^T`
multiplies by the capacitance -- `diffusion_constant` was wrong by
`C^2`, ratios 0.010003 / 1.000334 / 100.033536 over a 100x sweep --
and on ALGEBRAIC states `C^T` ANNIHILATES, so `c` came back EXACTLY
0.0 for an oscillator whose only noise was its series tank loss.
Neither was visible because every fixture here uses `C = 1 F`.

⚠⚠ AND `C` IS NEVER INVERTED.  The solve is on `C[D, NZ]^T`, the
block between the DIFFERENTIAL equations and the NON-algebraic
states, which is square and invertible by construction; the
singular `C` as a whole is not touched.  The algebraic entries come
from the constraint below instead.

### `_algebraic_adjoint_fill`

The docstring before the move:

Fill the adjoint's ALGEBRAIC entries, which are SLAVED, not free.

⚠ THE REPLAY LEAVES THEM AT ZERO AND ZERO IS NOT THEIR VALUE.  The
PPV entry for a row IS the phase sensitivity to a perturbation
entering that row, and an algebraic row's perturbation reaches the
dynamics through the CONSTRAINT rather than through its own row.
On a tank with series loss, eliminating `v_x = r (i_L + b)` puts
`-r b` into the inductor's equation, so the sensitivity to node
`x` is `r` times the branch row's -- and a noise current landing
there was being contracted against a structural zero, which made
`diffusion_constant` return EXACTLY 0.0 for an oscillator whose
only noise was its tank loss.  Measured against three independent
references; see the roadmap's section 0d.

The adjoint equation's ALGEBRAIC-STATE columns are what determines
them.  For a column `j` with no `C` entry the equation carries no
derivative, so it reads `sum_i G_ij v_i = 0`, and splitting `i`
into algebraic and differential rows gives

    v_A  =  (G[A, Z]^T)^-1 G[D, Z]^T v_D

⚠ THE MAGNITUDE IS STRUCTURAL AND THE SIGN IS MEASURED, and saying
which is which is the point.  On a RESISTIVE DIVIDER between the
inductor and ground the two algebraic nodes fold into the branch
row with coefficients `(r1 + r2)` and `r2`, so their entries must
stand in a ratio the topology fixes -- measured `10.000000` against
a chosen `10.000000`, on a fixture built so the single-resistor
degeneracy cannot hide a mistake.  ⚠ THE SINGLE SERIES RESISTOR
CANNOT SETTLE THIS: there `|integral v_0|` and `|r integral
v_branch|` agree to 1.5e-4, so BOTH SIGNS FIT and an agreement
there is no evidence at all.

The overall sign is then fixed by requiring algebraic and
differential rows to share ONE convention -- `dT/dA = +integral
v_j` -- and measured on the divider: with it, the DC-injection
probe reads +0.9999849 (node v, differential), +0.9999982 and
+0.9998744 (the two algebraic nodes); with the sign the naive
derivation gives, the last two come back NEGATIVE.

Returns `vblock` unchanged, with a warning, when the structure is
not index-1: `len(rows) != len(cols)` or a singular block.  That
case is section B4's, and guessing at it would be worse than
leaving a known zero.

The comment on the sign:

⚠ THE SIGN IS NOW THE DERIVED ONE, because this acts on `v_1`.
It was FLIPPED while this fill acted on `C^T v_1`: on that
fixture `C = diag(1, 0, -L)`, the INDUCTOR BRANCH ROW CARRIES
`-L`, and the term below is dominated by that branch -- so
reading `C^T v_1` as `v_1` negated it.  The derivation and the
measurement were describing DIFFERENT VECTORS and both were
right.  Pinned by the eq (24) constraint, which returns
0.0000e+00 exactly here and 1.9870e+00 for either alternative.

### `_ppv_propagate`

The docstring before the move:

The pair-consistent SECOND-ORDER propagation of an anchor vector
`v` over the period (the block `ppv()` applies to its null vector,
lifted 2026-09-08 so `frequency_aware_ppv` can run it on a COMPLEX
anchor).  Returns `(states, states_pair, ts, Xf)`: `states` the
per-step state-space samples (`C^T v`, second order, rescaled by
`v . xdot = 1` at the first sample), `states_pair` the raw
pair-space replay, `ts` the transposed per-step states, `Xf` the
orbit.  Only the `solved_history` (LMM) kind carries the
correction; the others return the replay as it is.

The comment on the pair's first block:

⚠⚠ THE PAIR'S FIRST BLOCK IS NOT THE PPV, AND THE ERROR IS FIRST
ORDER AND GROWS WITH Q.  For Gear-2 the adjoint state is the pair
`(w1, w2) = (dphi/dx_k, dphi/dx_{k-1})`, and `w1` alone is the
response to a perturbation of `x_k` WITH `x_{k-1}` HELD -- an
inconsistent history, which the two-step method resolves through
its parasitic root.  A physical state perturbation moves both:
`dx_{k-1} = Phi(t_{k-1}, t_k) dx_k`, so the phase functional is

    v(t_k) = w1 + Phi(t_{k-1}, t_k)^T w2,   Phi ~ I - h J + O(h^2)

and with `w2 = C_{k-1}^T z` (`z = -a2 t_k`, exact by the
recursion) that is `w1 + (C_{k-1} + h G)^T z` -- no inverse of
`C`, so it holds for a DAE.  Equivalently `w1` is orthogonal to
the amplitude eigenvector's first block, which is the true
amplitude direction ROTATED by `O(h)`; `v . xdot = 1` then
amplifies that rotation by `|v||xdot|`, the near-cancellation a
non-isochronous oscillator has (its PPV grows with `Q_lambda`).
MEASURED against the exact continuous adjoint (DOP853 at 1e-12,
no shooting code in the reference) on `vdp + 0.3 u^2`, whose
`c` is 100x van der Pol's: the first block gave `c` 16.6 / 8.0 /
3.9 / 1.9 / 1.0% high at 400..6400 points -- clean first order
-- and violated `v(t) . xdot(t) = 1` along the orbit by 12%
(std 2.7e-2).  This contraction holds the invariant to 8e-5 and
gives `c` to 1.8e-3 at 400 and 8e-5 at 1600, second order.  On
van der Pol both agree to 1e-4: the two rows are in quadrature
there, so the rotation averaged out of `<v^2>` -- the fixture
shared the claim's assumption (failure shape 0b), and `pnoise`,
which contracts in PAIR space, was right all along and 14%
below `c` on the fixture that could see it.
The seed's scale is `w1(0) . xdot = 1`; the consistent object
is renormalised ONCE by its own `v(0) . xdot`, which is why the
per-step invariant is the test and not the definition.

The comment on the differential rows:

⚠ DIFFERENTIAL ROWS ONLY.  `w2 = C^T z` does not see the
algebraic rows of `z` (their rows of `C` are zero), so
the decomposition is non-unique there, and those
multipliers are O(1/h): `h G^T z` would carry an O(1)
component along the constraint normal into the
differential entries.  The consistent propagation
`C_D dx_{k-1} = (C_D + h G_D) dx_k` involves only the
differential equations, which is the choice that makes
it unique.  Measured: with the algebraic rows in, a DC
injection probe on a series-loss tank flipped sign.

The comment on the Schur complement:

⚠ ON A DAE THE ALGEBRAIC STATE IS SLAVED, AND ITS
COUPLING INTO THE DIFFERENTIAL PROPAGATION IS O(h).
A consistent perturbation propagates as
`C_D dx_{k-1} = (C_D + h G_red) dx_k` on the
differential states, with `G_red` the Schur
complement `G[D,NZ] - G[D,Z] G[A,Z]^-1 G[A,NZ]`.
With the full `G` instead, the series-loss tank had
`c` 0.6 / 0.3 / 0.15% low at 240/480/960 points
(first order) and the invariant drifting at 1.1e-3;
with the complement `c` is 2.7e-4 / 7e-5 / 2e-5 from
the exact reduced-ODE value and the drift 4e-4 /
1.1e-4 / 2.7e-5 -- second order (found through the
review session's linear-DAE partition, 2026-09-05).

The comment on the index-1 condition:

⚠ AND `G[A,Z]` NONSINGULAR *IS* THE INDEX-1 CONDITION.
At index >= 2 (an L-I cutset, a C-V loop) it is
singular by definition and the algebraic variables
come from a differentiation, not a solve; the
complement does not exist.  Same shape as the fill
above: warn once with the reason and fall back to
the full-`G` propagation, which is then first order.
(Boundary named by the review session from the
pencil: `eig(-G_red, C[D,NZ])` equals the finite
generalised eigenvalues of `(C, G)` to 1e-12 on the
series-loss tank, and the reduction is undefined on
`li_plus_rc` and `cv_plus_rc`.)

The comment on the algebraic columns:

⚠ AND ZERO ON THE ALGEBRAIC COLUMNS, as `C^T v_1` is:
the state functional contracts a perturbation ON the
constraint manifold, whose algebraic components are
slaved, and `h G^T z` would otherwise leave 4e-3 there
(caught by the full suite's Demir-(24) gate).  The
equation-row conversion never reads these entries.

The comment on the DC content:

⚠ AND ITS DC CONTENT IS THE CONSISTENT OBJECT'S TOO -- taking
the mean from the raw block was TRIED AND MEASURED WRONG.
The raw block's orbit integral reproduces a same-grid
DC-injection probe to 1e-5 on the divider fixture (node row,
true mean 4e-6 |v|), where the consistent object's O(h^2)
pointwise errors leave an absolute floor of ~1e-5 |v| --
the wrong sign at 480 points.  But on the bias-sensitive
fixture's INDUCTOR row (a DC voltage in series with L, true
dT/dV = 16.20 by a second-order re-solve) the raw block
reads 17.49 / 16.83 / 16.51 at 400/800/1600 -- first order,
8% off -- while the consistent object holds `v . xdot = 1`
to 8e-5 along the orbit, which pins its mean in EVERY row to
~1e-5 |v|.  Stitching the raw mean in broke that invariant
by +-0.3.  So the raw block's DC exactness is row- or
fixture-specific (mechanism open, recorded in the roadmap),
and `samples` is one object, second order everywhere, with a
~1e-5 |v| absolute floor on its mean.  The raw pair is kept
as `samples_pair` for the structural gates that live on its
discrete identities.

### `ppv`

The docstring before the move:

The perturbation projection vector at `t = 0` (Demir & Roychowdhury).

Returns `(v, info)`.  `v` is the pair-space left null vector of
`I - M`, normalised so that `v . xdot(0) = 1` -- see the note on the
normalisation below, which was MEASURED rather than transcribed.
The phase shift caused by a state perturbation `delta` at `t = 0` is
then `v[:m] . delta`.  `info` carries both border residuals, the
null residual, `q` and the scaled tangent.

⚠ AN AUGMENTED SOLVE, NOT AN EIGENVECTOR -- AND IT IS THE FIX FOR A
NAMED FAILURE OUR OWN FIXTURES SIT INSIDE.  Demir &
Sangiovanni-Vincentelli, 1998 (the book, read firsthand by the docs
session 2026-09-04), report the eigenvector route BREAKING on a
high-Q circuit, with a table of the crowded eigenvalues (their
Table 6.4): "Since this circuit is a high-Q one, Phi(T,0) has
eigenvalues with magnitudes close to 1 other than the one which is
supposed to be equal to 1 ... Because of numerical errors, we can
not identify the eigenvalue that is supposed to be equal to 1 ...
so it is not feasible to identify the correct" one.  Not
ill-conditioned there -- INFEASIBLE.  That reported failure is the
stated motivation for the single-solve method two years later, so
the lineage is firsthand end to end: 1998 selection fails at high
Q; 2000 the single linear solve; 2001 "particularly useful for
high-Q oscillators"; 2003 the fuller procedure.  ⚠ It also changes
what the second-multiplier warning below MEANS: not "this result is
degrading" but "you are in the regime this method was invented to
escape".  And the same book's eq (6.72), `|exp(eta_i)| << 1`, is
the closed-form variance's validity condition -- the book says it
"is satisfied for 'most' oscillator circuits" and defers the rest
to the high-Q discussion above, so that condition and the crowding
are ONE condition seen from two sides.

⚠ ATTRIBUTION CORRECTED
2026-09-04: the idea ORIGINATES in Demir, Long & Roychowdhury,
ICCAD 2000 ("Computing Phase Noise Eigenfunctions Directly from
Steady-State Jacobian Matrices" -- "a single linear solution of the
oscillator's ... steady-state Jacobian matrix ... dispenses with the
need to select the correct one eigenfunction"), with the 2001
companion carrying it to HB/shooting matrices and noting the
selection heuristic is worst "for high-Q oscillators".  The 2003
paper is the fuller procedure and the source of the quote below,
not the origin.  Demir's IJCTA 2000 method SELECTED the
right eigenvector by its inner product against `C(0) xdot(0)` --
measured 0.2 against 1e-5, 1e-7, 2e-5 on a Colpitts.  His 2003
paper rejects that: "no guarantee that any of the candidate
eigenvectors will be appreciably more orthonormal than the others,
leading to a potential breakdown."  The same vector then changes
role -- it becomes the BORDER, so the candidate is unique and
nothing is selected:

    [ I - M^T   q ] [v]     [0]
    [   q^T     0 ] [y]  =  [1]

This matters here specifically.  `_spectral_report`'s eigenvector
split was measured labelling a parasitic root physical at ~2 points
per cycle, and multipliers cluster near 1 on exactly the high-Q
oscillators a PPV is wanted for -- four independent sources say so.
A bordered solve does not care how close the other multipliers are;
it never has to tell them apart.

⚠ AND THAT IS THE METHOD'S STATED DESIGN DRIVER, not a lucky
property of it.  Demir, Long & Roychowdhury (ICCAD 2000), who
introduced the single-solve route: it is "especially advantageous
for HIGH-Q OSCILLATORS, MONODROMY MATRICES OF WHICH OFTEN HAVE MANY
EIGENVALUES CLOSE TO 1 THAT ARE NUMERICALLY INDISTINGUISHABLE from
the oscillatory [unit eigenvalue]", and "a key advantage is that it
DISPENSES WITH THE NEED TO SELECT THE CORRECT ONE-EIGENFUNCTION
from amongst a potentially large set of choices".  So the hardest
case in this codebase is the case the method was aimed at.

⚠ WHICH DOES NOT RETRACT THE SLOW-NODE BOUNDARY BELOW, and keeping
the two apart is the point.  Near-degenerate multipliers make
EIGENANALYSIS ILL-POSED -- there is no fact of the matter about
which eigenfunction is the PPV -- while they make this bordered
solve merely ILL-CONDITIONED, which is measured above and warned
about.  Ill-conditioned beats ill-posed, and neither is the same as
the PHASE EQUATION's own limit, which is about the response being
treated as instantaneous and is not an extraction question at all.
Three separate things that a "high-Q oscillators are hard" summary
would blur into one.

⚠ THE QUADRATIC RUNG BELOW WOULD NOT FIX THE SLOW-NODE BOUNDARY,
which is the obvious hope and is wrong.  TWO INDEPENDENT
approximations are in play: the LINEAR ISOCHRON one is in the
perturbation's AMPLITUDE -- it treats isochrons as flat
hyperplanes, and the quadratic rung adds their curvature -- while
the INSTANTANEOUS-RESPONSE one is in the DYNAMICS, ignoring the
bandwidth between injection point and core.  Slow nodes are the
second.  Noise is small by construction, so the linear term
dominates it by definition; the quadratic rung would earn its cost
on LARGE perturbations -- injection locking, big supply or
substrate interferers -- not on phase noise.

⚠ AND THE PPV IS ONE RUNG ON A LADDER, worth knowing before it is
mistaken for exact.  Suvak & Demir (TCAD 2011) place it: an EXACT
phase equation exists and is "practically unusable"; the PPV
equation is its LINEAR isochron approximation; a QUADRATIC one is
more accurate.  Isochrons are the geometric form of asymptotic
phase, so an oscillator without asymptotic phase is one whose
isochrons do not exist -- the same fact as the Floquet condition,
seen in the geometry.  Computing exact isochrons is exponential, so
the only live question is which local approximation is affordable.

`c` -- the diffusion constant this vector feeds -- has the
designer-facing reading "JITTER PER SECOND".

⚠ THREE NAMES FOR THIS OBJECT, AND ONE NEAR-MISS THAT IS NOT IT.
The PPV, Kaertner's adjoint LPTV impulse response, the PRC of
mathematical biology, and Hajimiri's NUMERICAL ISF are the same
thing.  His CLOSED-FORM ISF is NOT: it is the normalised tangent,
and the difference is not a scale factor but a SIGN -- for a noise
impulse at one point in the cycle "the closed-form ISF predicts a
POSITIVE phase change, whereas in fact the correct phase change is
in the NEGATIVE direction and of a different magnitude".  It also
does not scale with the perturbation, where the PPV does.

⚠ AND `xdot` IS NOT A CHEAP SUBSTITUTE FOR IT -- "time-shifts and
amplitudes are both different ... the two waveforms scale in
OPPOSITE DIRECTIONS with respect to the RC time constant".  Nothing
here offers it as one: `xdot` appears only as the NORMALISATION
(`v . xdot = 1`, the PPV's defining property), as the border `q =
C(0) xdot(0)`, and in the record above of the SELECTION heuristic
that was rejected.  Stated because the substitution is a documented
point of common confusion, and the failure would be a sign error
rather than a visible one.

⚠ `y` COMING BACK ZERO IS A FREE CORRECTNESS CHECK, and it is not
decoration.  With a zero first block on the right-hand side,
`(I - M^T) v + y q = 0` forces `y q = 0`, so a nonzero `y` means the
border absorbed a residual the null space should have taken -- the
computed `v` is not in the null space.  Measured on van der Pol:
1.4e-11.

⚠ ITS VALIDITY BOUNDARY IS SLOW NODES, AND THE VAN DER POL GATE
CANNOT SEE IT.  The phase equation `alpha' = v_1^T(t+alpha) b(t)`
treats the oscillator's frequency response as INSTANTANEOUS; the
truth is a convolution, and the PPV form is what you get by
assuming the kernel is `v_1(t) delta(t - tau)`.  Real circuits have
finite bandwidth, so a slow node FILTERS the noise of devices near
it, the PPV cannot see the filtering, and phase noise is
OVER-ESTIMATED.  Lai (Cadence) is explicit that better extraction
does not help: "although the PPV can be extracted correctly, the
oscillator noise analysis is still inaccurate: the phase noise is
always over-estimated."

⚠ AND HE NAMES THIS TEST'S OWN REGIME AS THE BLIND SPOT: "the
phase equation was verified to be correct in many previous works
... because it was evaluated on SMALL, SIMPLE OSCILLATORS, and
perturbations were applied to OSCILLATOR CORES.  Since oscillator
cores have very wide bandwidth, ignoring the dynamics may not
compromise the macromodelling accuracy very much."  That is
`test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has`
exactly -- van der Pol, perturbed at its core.  It passes whether
or not this failure mode is present, so it establishes that the
extraction and normalisation are right and says NOTHING about the
model's range.  The fix is a frequency-aware PPV, which is this
same bordered system at nonzero `w_s` -- the classical PPV is its
DC point, and `PAC` already solves at nonzero frequency.  Not
built.  VERIFIED at the source (docs session, 2026-09-08): Lai
2008 eq. (24) at `w_s = 0` "is the augmented PPV extraction
equation (6) and (7)", verbatim; two scope limits: (24) is a
NEAR-DC approximation of (23) (the AC Toeplitz columns dropped,
"if we are only interested in ... w_s close to DC"), and (23)
"is very difficult to solve using iterative solvers (such as
GMRES)" because the border degrades the block preconditioner.
⚠ AND `_vdp_with_slow_node` CANNOT SHOW THE EFFECT AT ANY tau/T:
its slow node couples through `Rs = 1e6` against a tank impedance
of 1, so its PPV entry is 7.3e-6 of the core's, flat over
tau/T = 1e2..1e6 while lambda_2 moves four decades -- it tests
conditioning, not the filtering of noise that REACHES the phase.
The fixture owed is a slow node IN the phase path (small Rs, large
Cs), with tau/T and coupling as separate knobs.  BUILT AND MEASURED
the same day: with an asymmetric core AND tank loss (G_0 != 0 needs
both) the slow node's PPV entry is DC-dominated (|G0|/|G1| = 30) and
the Lorentzian over-states a source behind it by 1000x at 0.1 f0,
predicted to four digits from this PPV's harmonics and the RC
filter -- so `c` from this PPV is right and the SHAPE above
T/(2 pi tau) is what the frequency-aware PPV corrects; see
`oscillator_spectrum`.  A2's "gated at tau/T = 10" was a Monte
Carlo of `c`, which cannot see it.

⚠ `q` IS EXACT, NOT DIFFERENCED.  `q = C(0) xdot(0)` looks like it
needs the orbit's tangent, and differencing the waveform for it
would be O(h) at best.  The DAE gives it directly: `dq/dt + i(x) +
u(t) = 0` and `dq/dt = C xdot`, so `q = -(i(x_0) + u(0))` -- two
evaluations at the converged solution, no derivative anywhere.

⚠ ON A STAGED SOLVE (`state_events=True`) THE PPV IS BORDERED
(2026-09-22, events phase B): the null vector is that of the TOTAL
monodromy `M + P_theta dtheta/dx_0`, and the samples carry the
crossings' motion as costate injections `-zeta_k W_k` at the event
nodes, `zeta = Gt^-T P_theta^T v` -- one reverse pass, and it IS
the phase gradient at fixed time (the tail-restricted Newton plus
the node's own motion agrees with it to 1e-4).  VERIFIED against
the exact piecewise-linear PPV of a comparator relaxation
oscillator (linear flows joined by saltation matrices): every
sample within 0.4 % before, between and after the crossings,
where the fixed-grid PPV of the same solve is 140 % off before
the first crossing and its `M^T v - v` residual is 1.5.  ⚠ THE
TRANSIENT FD THAT WAS TO BE THE INSTRUMENT read `c` at +4.4e-8
for the exact -7.4e-9 s/V, converged in its own step, period to
4e-6, along the flow to 1e-4: its phase was read at the `c`
waveform's crossing of ``(max + min) / 2`` OF THE PERTURBED RECORD
-- a level the perturbation itself moved.  Read at the
comparator's own crossing it agrees with the exact value to 0.4 %.
An instrument's reference level must not come from the record it
measures.

The comment on the map kinds:

⚠ NO LONGER GEAR-ONLY (B8). The refusal that stood here said the
transposed replay was "implemented for the solved-history map
only"; since `_monodromy_matvec_transposed_plain` shipped that
sentence is false, and every call below goes through
`fp.matvec_transposed`/`fp.matvec` so the map's kind is the
dispatcher's business rather than this method's.

The comment on the staged solve:

⚠ ON A STAGED SOLVE THE MONODROMY IS THE TOTAL ONE (2026-09-22,
events phase B): a perturbation moves the crossings, `M_tot = M
+ P_theta dtheta/dx_0`, and the PPV is its left null vector.
Its samples along the orbit carry the same correction as costate
injections at the event nodes -- `-zeta_k w_k`, `zeta = Gt^-T
P_theta^T v` -- which the reverse pass carries to every earlier
node: the saltation matrix's transpose, derived from the bordered
system rather than guessed.  A comparator oscillator's PPV jumps
at its switching instants, and this is where the jump comes from.

The comment on the normalisation:

⚠ NORMALISED BY `v . xdot = 1`, WHICH IS NOT WHAT `v . q = 1`
GIVES, and the difference is not cosmetic: on `_vdp_ppv(400)`
`v . xdot = 1.0` against `v . q = -1.0696` -- a factor 2.07 AND
the opposite sign (an earlier version of this note said "7%",
which was the |v . q| - 1 residual and not the error; corrected
by the review session's audit, 2026-09-04).  The defining
property is that displacing the state
ALONG the orbit by `eps xdot` advances the phase by `eps`, so
`v . xdot = 1` is the normalisation a state perturbation sees.
Demir's Remark 3.1 reads `v_1^T C u_1 = 1`; the vector this
bordered solve returns behaves as `C^T v_1` -- it is contracted
with a state perturbation directly -- so the two statements agree
about different objects.

⚠ AND THAT IS NOT A QUIRK OF THIS FORMULATION, WHICH THIS
COMMENT USED TO IMPLY.  The conserved pairing propagates to
`M^T (C(0)^T v_1) = C(0)^T v_1`, so the left eigenvector of the
STATE-SPACE monodromy simply IS `C(0)^T v_1` -- for ANY `C`,
symmetric or not, and whatever the augmentation.  MEASURED on a
limit cycle with a constant NON-SYMMETRIC `C`: alignment with
`C(0)^T v_1` is 1.000000000000 against 0.9657 for `v_1` itself,
and bordering with `xdot(0)` reproduces Demir's normalisation
exactly while bordering with `C(0) xdot(0)` gives 0.805.
(Derived and measured by the docs session, 2026-09-04.)  ⚠ TREATING THEM AS THE SAME OBJECT WAS
MEASURED WRONG: predicting a state jump's phase shift as
`v^T C delta` gives residuals of 0.36/0.40/0.42 that GROW with
refinement and per-direction ratios scattering from -0.44 to
28.7, while `v . delta` converges at O(h).
⚠ THE ALGEBRAIC ENTRIES ARE FILLED *AFTER* THIS, AND THAT IS A
DECISION RATHER THAN AN ORDERING ACCIDENT.  Filling first was
tried and MEASURED WORSE: the DC-injection probe went from
0.9999849 to 0.9992364, because `v` here is also the REPLAY'S
SEED and the algebraic components are SLAVED -- propagating them
through the step map corrupts the differential ones.

⚠ AND THE NORMALISATION SHOULD NOT SEE THEM EITHER.  `v . xdot`
is about a STATE perturbation, and a state perturbation of a DAE
lies ON the constraint manifold: its algebraic components are
determined by its differential ones, not free.  The algebraic
entries of `v` answer a different question -- the sensitivity to
a perturbation of an EQUATION ROW, which is what a noise current
injected into an algebraic KCL row is.  So this line is
unchanged, and every PPV number on every circuit is
bit-for-bit what it was.

The comment on the unit multiplier:

⚠ THE UNIT MULTIPLIER OFF THE CIRCLE IS A SILENT ERROR IN `c`
(2026-09-21).  This solves the bordered system AT lambda = 1
whatever the discrete period map's own unit multiplier is, so
nothing here can see that multiplier sit at 1.10 -- and on a
relaxation van der Pol (mu = 10) on its own 195-point `lte_grid`,
gear's did: the PPV then gave a diffusion constant 52 % high (14.7 %
at 2x, 4.0 % under trbdf2 at |rho - 1| = 3.4e-3), tracking the
departure at 5-12x, with `converged = True` and no message.  The
solve already records `spectral_radius`; an autonomous run whose
unit multiplier is more than `PPV_UNIT_MULTIPLIER_WARN` off the
circle is told so here, once, with the size and the remedy.

The comment on the samples over the period:

⚠ THE PPV OVER THE PERIOD, not just at `t = 0`, because that is
what an oscillator noise calculation needs: Demir's diffusion
constant is `c = (1/T) integral v_1^T(t) B(t) B^T(t) v_1(t) dt`,
an integral over the orbit.  `Phi(T,s)^T v(T) = v(s)`, and the
reverse replay computes exactly that sequence on its way to the
answer -- it was being discarded.

The comment on index >= 2 on a non-uniform grid:

⚠ INDEX >= 2 ON A NON-UNIFORM SOLVED-HISTORY GRID (2026-09-20, item
3 of "gear as a first-class choice on non-uniform grids"): the
fallback above keeps only the differential block of the
pair-consistent correction, and the coupling it drops -- a
DERIVATIVE term at index 2, with no Schur complement to stand in
for it -- cancels between steps on a uniform grid and not on a
non-uniform one.  Measured on a van der Pol with a DC source inside
a capacitor loop, smooth grid, `c` against radau: -5.0e-3 / -2.4e-3
/ -1.2e-3 / -5.8e-4 at N = 100..800 (FIRST order) where the uniform
grid gives +1.2e-3 / +3.3e-4 / +8.7e-5 / +2.2e-5 and trap on the
same smooth grid is second order.  So on that grid the samples come
from the continuous adjoint's phase mode instead -- the object
`floquet_modes` already uses there, whose invariant quarters on
this fixture -- as the STATE-SPACE PPV `v_j = C_j^T q_j` (the
continuous `q` IS Demir's equation-row `v_1`, and `samples` carry
`C^T v_1` -- see `_equation_row_ppv`; `q` put there bare lands
16.7x off with the right shape), scaled so `q_0^T C_0 xdot_0 = 1`,
which is `ppv()`'s own `v . xdot = 1`.  Measured with a periodic
trapezoid over the actual grid: +4.0e-5 / +3.6e-5 / +1.3e-5 at
N = 200 / 400 / 800, the reference's floor.  ⚠ The second half of
that measurement was PAC's period weights (`_period_weights`):
a left-rectangle rule on a non-uniform grid is first order by
itself.  The anchor `v` and the pair block are untouched; the
uniform grid never comes here, and the one-step kinds cannot:
their `factored_period_full` / `_dirk` replays (now
`factored_period_stage`) were built on a uniform `linspace` grid
whatever the solve used, which is WHY
radau read as "exact on the 3:1 grid" -- its adjoint surfaces
never saw that grid.

The comment on the per-sample fill and the equation-row adjoint:

⚠ AND FILL EVERY SAMPLE TOO, at ITS OWN operating point, because
`G` is state-dependent and the algebraic entries are a pointwise
function of the differential ones.  Done here rather than by
seeding the replay: these components are SLAVED, so there is
nothing to propagate, and post-processing leaves the validated
step map untouched.  The pair's SECOND block is the history term
and is deliberately not filled -- `v(t)` is the first block.
⚠⚠ THE EQUATION-ROW ADJOINT IS A SECOND OBJECT, NOT A CORRECTION
TO THE FIRST.  `states` and `v` stay exactly what they were --
`C^T v_1`, the vector a STATE perturbation contracts with, which
is what this method documents and what every existing gate
measures.  Demir gives both conventions on one page (eq 41 with
`C`, eq 42 and the phase equation 44 bare), so naming both is
the fix; converting one into the other would have silently
changed what `ppv()` returns.
⚠ ONE WARNING PER CALL, NOT ONE PER SAMPLE: the fill warns when
`G[A,Z]` is singular, and at index 2 it is singular at every
sample -- 240 identical warnings for one call, which trains a
reader to filter this module's warnings and miss a real one.

The comment on the second multiplier (before `vu`):

⚠ A SECOND MULTIPLIER NEAR 1 BREAKS THIS SILENTLY, and none of
the residuals above can see it.  The border removes the PHASE
mode's singularity and does nothing about any OTHER root
approaching the unit circle -- which a slow node puts there.
MEASURED on van der Pol with one weakly coupled RC node:

    tau/T    |lambda_2|    sigma_min(bordered)   null residual
    none      0.000856          8.62e-01            4.1e-11
    1e2       0.990049          4.49e-03            4.6e-11
    1e4       0.999900          4.47e-05            4.6e-11
    1e6       0.999999          4.47e-07            4.4e-11

`sigma_min` tracks `T/tau` over six decades while the residual
does not move at all: GMRES converges, the answer looks clean,
and the conditioning has lost six digits.  So this estimates
`|lambda_2|` explicitly rather than trusting a small residual.

⚠ AND THE ACCURACY COST WAS GATED, WITH A NEGATIVE RESULT worth
recording so nobody re-derives a fix from the warning alone.
Monte Carlo on the FULL NONLINEAR circuit -- 200 realisations,
150 periods, phase read from zero-crossing timing, so no PPV
appears anywhere in the measurement:

    core injection (control)   c_ppv/c_meas = 0.9965
    slow node, tau/T = 10      c_ppv/c_meas = 0.8016

Within 20%, about 2 sigma at this sample count, and in the
UNDER-predicting direction.

⚠ BUT THAT IS NOT A FALSIFICATION, AND THIS DOCSTRING SAID IT
WAS.  `tau/T = 10` is OUTSIDE the regime the reported mechanism
needs: it bites through ill-conditioning, and by the table above
`sigma_min` at `tau/T = 10` is ~4.5e-02 -- healthy.  The PPV has
no large entries there and nothing is splitting into two nearly
cancelling components.  Lai's own case is a gated-capacitor
tuning bank (226 MOSFETs, 3.15 GHz) whose off-caps have RC
exceeding ~1 s, i.e. `tau/T ~ 3e9` -- eight orders from what was
tested.  A null result at 10 is what the mechanism PREDICTS, not
evidence against it.

⚠⚠ PROVENANCE, AND THE CHAIN IS NOW FULLY TRACED -- A UNIT WAS
MANUFACTURED IN TWO STEPS.  This used to render "larger than 1
second" AS A QUOTATION.  The primary source IS on disk, at
`~/docs/09-phase-macromodels-and-prc/Lai-2008-Frequency-Aware
PPV ... (Cadence).pdf`, and p.4 reads, verbatim:

    "Since the RC time constants of the "off" gated capacitors is
     very large (LARGER THAN 1), it is safe to assume that these
     gates have very small contribution to the total phase noise
     when offset frequency is reasonably large."

**NO UNIT.**  Our own reading of the paper
(`~/docs/pycircuit-frequency-aware-ppv.md`) paraphrased it as
"their RC constants exceed 1 s" -- ADDING the unit, and unmarked,
beside that file's properly marked quotations.  This comment then
promoted the paraphrase to a QUOTATION, carrying the added unit
with it.  Two steps, each small, and the result was a quoted unit
the source does not contain.

⚠ Seconds remains the natural reading (the `tau/T ~ 3e9` above
follows from it and nothing downstream moves), but it is OURS and
is marked as such.  ⚠⚠ AND THE FIRST VERSION OF THIS CORRECTION
SAID THE PDF WAS "NOT ON DISK AT ALL" -- it is, in a
SUBDIRECTORY, and the search that missed it looked only at the
top level of `~/docs`.  Search a library recursively before
reporting a source missing.

⚠ SO THE HONEST RECORD IS: not reproduced at `tau/T = 10`, which
is outside the regime where the mechanism predicts an effect;
UNTESTED at the `tau/T ~ 1e9` where it is reported.  And the
reason the fix is still not built is COST, not falsification:
the measurement needs ~15 time constants of settling, so at
`tau/T = 1e4` that is 150 000 periods per realisation.  That
argument stands on its own; the falsification framing does not,
and this codebase's ledger distinguishes them.

⚠ AND THE 0.80 IS IN THE OPPOSITE DIRECTION TO THE REPORTED
EFFECT.  If it survives the ~10% Monte Carlo uncertainty at 200
realisations it is a separate ~20% UNDER-prediction at a `tau/T`
where the conditioning is fine -- not a weak version of Lai's.
At ~2.5 sigma it is not established either way, and it is
recorded rather than resolved.

⚠ Larger `tau/T` is untested and the cost is why: the
measurement needs ~15 time constants of settling.
⚠ AND IT TOOK THREE ATTEMPTS.  A window of 2-4 time constants
read the slow mode's DECAY as diffusion; an impulse test could
not resolve a 1e-11 time shift; and one noise amplitude for both
circuits put a 2.5 V jump per step on an orbit of amplitude 2,
because the slow node's capacitance is 6.7e-5 F against the
core's 1.0.  Each time the number was read before the
MEASUREMENT was shown to be in the regime it assumes.

⚠ ARNOLDI RITZ VALUES, NOT A DEFLATED POWER ITERATION -- and the
replacement is BOTH more accurate and cheaper, which is rare
enough to state plainly.  Power iteration converges at
`|lambda_3|/|lambda_2|`, so it fails exactly where a parasitic
multiplier crowds the oscillatory one.  Arnoldi does not care
about that ratio.  MEASURED on van der Pol at `Q = 16` with one
parasitic RC swept through it:

    tau_p/T   lam2/lam3   POWER err   RITZ err
      1        2.554      3.53e-14    0
      4        1.206      1.52e-06    1.11e-15
      8        1.065      1.41e-03    0
     16        1.000      1.10e-05    2.22e-16
     32        1.032      4.61e-03    2.22e-16
    100        1.054      2.48e-03    2.22e-16

Machine precision everywhere INCLUDING at exact degeneracy,
against a power iteration losing three digits at a ratio of
1.065 -- and at `k` matvecs rather than
`PPV_DEFLATION_ITERS = 30`.

⚠ THE ROUTE IS GARCIA, ROMERO & ACHA (IEEE Trans. Power
Systems 37(1), 2022): Ritz values of `I - M` map back as
`lambda = 1 - theta`.  They take `H` from the GMRES that
already solved the Newton correction; here it is a small
dedicated Arnoldi, because this call has no GMRES of its own.

⚠ EXACT AT `k = n` AND A TRUNCATION OTHERWISE.  The cap keeps a
large circuit from paying `n` matvecs for a diagnostic.

⚠ A TRUNCATED `lam2` IS A LOWER BOUND, so the near-unit warning
below can only UNDER-fire -- and that is Cauchy interlacing,
not an accident: `theta_j >= lambda_j(A)` for a Rayleigh-Ritz
projection, hence `1 - theta_2 <= lam2`.  MEASURED on a
synthetic 40x40 with a verified-normal `M`
(`||M^H M - M M^H||/||M||^2 = 8.9e-16`): a lower bound in
100.0% of 200 draws at `k` = 3, 5, 8 and 12.

⚠⚠ AND THIS SURVIVED A ROUND TRIP THROUGH A FALSE REFUTATION,
which is why it is written out.  An intermediate version of
this comment said "over-estimate, not a lower bound", on a
measurement whose selection rule took the largest `|1 - theta|`
after discarding `|lam - 1| < 1e-8`.  On a truncated basis that
picks the UNCONVERGED UNIT-MODE Ritz value -- below 1 but above
`lam2` -- so it measured its own filter and attributed the
result to the phase mode contaminating `lam2`.  Selecting
`theta_2` as the second-smallest Ritz value of `I - M`, which
is what interlacing is about, restores the bound.

⚠ THE BOUND IS FOR A NORMAL `M`.  At forced eigenvector
conditioning `cond(V) = 1e6` it fails in 100% of draws at
`k = 20` -- but by exactly `1 - lam2`, i.e. the unit mode being
SELECTED as `lam2`, a selection failure rather than a Ritz one:
at that conditioning `|lam_1 - 1|` is 1.5e-8 to 4.3e-8 and no
value-based rule separates it.  A fixed absolute tolerance is
the weak point; deflating the phase mode explicitly with `q`,
which this method already has, would sidestep it.

⚠⚠⚠ AND THE ESCAPE CLAUSE IS LOAD-BEARING: A CIRCUIT MONODROMY IS
NOT NORMAL, AND THE BOUND FAILS ON ONE.  MEASURED 2026-09-07 on
THIS FILE'S OWN `_osc_with_ladder(Q, 14, nslow)` against the dense
spectrum of the SAME operator (`n` matvecs, the route the
`dirk`/`full` branch below already takes), scored in the GAP
because `Q ~ 1/(1 - lam2)`::

    nslow   dense lam2     k=12 Arnoldi   gap ratio   Q dense/Arnoldi
     <=11   0.995706203    0.995706197      1.000       232 / 232
       12   0.996324417    1.000114048     -0.031       271 / inf
       13   0.996818781    0.942674586     18.020       313 / 16.9
       14   0.997220139    0.999318472      0.245       359 / 1467

**THE ERROR IS NOT ONE-SIGNED**, so "can only UNDER-fire" does not
hold here: 13 under-estimates (19x low in `Q`), 12 and 14
OVER-estimate, and 12 returns `lam2 > 1` -- a spurious UNSTABLE
multiplier, which `Q` reports as `inf`.  The two failures are
different: at 13 the Arnoldi never resolves `0.99682` and selects
the next TRUE eigenvalue down (`0.9427`); at 14 it selects a
SPURIOUS Ritz value at `0.99932` that is no eigenvalue at all.

⚠⚠ SO "NOT LIVE ON A CIRCUIT MONODROMY" (below) WAS MEASURED ON THE
WRONG AXIS.  It was checked against EIGENVECTOR CONDITIONING; the
trigger here is the number of distinct near-unit CLUSTERS, which
`_osc_with_ladder` varies BY CONSTRUCTION and which the fixture's
own test already reports reaching ~29 at `nslow = 14`.  It is not
Q-specific either: the same `nslow` fails at Q = 8/16/256.

⚠⚠ IT IS A SIZING PROBLEM, AND RAISING THIS CONSTANT IS THE WRONG
FIX -- MEASURED.  `k = 16` is exact on the fixture above and FAILS
on a longer ladder, because the required `k` grows with `n`::

    ladder/nslow   n    k=12 ratio   k=16 ratio   k=20 ratio
       14 / 14     32      0.245        1.000        1.000
       20 / 20     44      0.277        0.279        1.000
       26 / 26     56      2.313        0.410        0.265

`k ~ n/2` and rising, against a DENSE route that costs `n` and needs
no threshold at all.

⚠⚠ BUT `k ~ n/2` IS AN ARTEFACT OF THIS FIXTURE, AND THE FIXTURE
CANNOT SEE IT.  `_osc_with_ladder` sets `nslow = nladder`, so `n`
and the slow-mode count move together here and no measurement on
it can separate "k tracks n" from "k tracks nslow".  A peer
session's synthetic CAN separate them and reports `k_min` rising
with `nslow` and FLAT under a doubling of `n` at fixed `nslow`
(24->24, 24->16, 48->48, 48->48).  If that transfers, the rule is
**cost tracks the SLOW-MODE COUNT, not the system size** -- a
large fast circuit is cheap and a small one with a big tuning
bank is not, which also says Lai's 813-equation oscillator is
expensive because of the BANK and not the 813.  Recorded with
that provenance: measured on a synthetic, consistent with
everything measured here, and NOT separable on this fixture.

⚠⚠ BUT THE DENSE ROUTE IS OUT ON THE CIRCUITS THAT MOTIVATE THIS.
`FLOQUET_DENSE_LIMIT = 400`, and the published cases are LARGER:
Lai's 64-gated-capacitor DCO is "about 200 transistors, and the
system size is more than 500 ... We have trouble to apply direct
harmonic balance in this case due to memory issue" (DAC 2006
p.1021, verified on disk), and [L08]'s tuning oscillator is 813.
So dense is the right default only in the `n <= 400` band this
class already draws, and the RITZ-RESIDUAL gate is what the large
end needs.  ⚠ And a 64-element bank with any realistic fraction
off is an order of magnitude past the >= 3 clusters that break
`k = 12` -- i.e. the extension AT ITS CURRENT BASIS SIZE would
fail on exactly the circuits it exists for.
THE DIAGNOSTIC THAT SEPARATES THEM CLEANLY IS THE PER-PAIR RITZ
RESIDUAL `|h_{k+1,k}| |y_i[last]|`, free from `H`: 1.0e-02 at
k=8, 2.1e-03 at k=12 (both wrong), 1.5e-16 at k=16 (right), and
<=3.1e-07 at every `nslow` the shipped path gets right.  Neither
is built -- see the roadmap; `lam2` and `Q` are REPORTED
DIAGNOSTICS with no non-test consumer, so nothing computes wrong,
but a caller reading `info['Q']` on a bias network with many long
time constants can be off by 4x to 19x, silently.

⚠ ONE FAILURE MODE CHECKED AND NOT LIVE HERE -- ⚠⚠ SUPERSEDED BY
THE MEASUREMENT ABOVE, KEPT BECAUSE IT RECORDS WHAT WAS TESTED.  At
`cond(V) >= 1e4` the `|lam - 1|` filter itself fails: the phase
mode stops being resolved to the tolerance, survives the
discard, and is selected as `lam2`, sending `Q` to infinity.
MEASURED on what was then this class's stiffest realistic fixture
-- a Q=60 oscillator with a 10-mode damped bulk, `m = 12` --
`cond(V) = 92` and `|lam_1 - 1| = 3.0e-13`, seven orders inside
the 1e-6 filter.  That axis is still clean; the CLUSTER-COUNT axis
is not.  Sorted by real part, not magnitude, because an
amplitude mode is real and positive while a complex pair of
larger modulus would be an oscillation about the orbit.

The comment on the dense route:

⚠⚠ DENSE WHENEVER IT IS AFFORDABLE, AND THAT IS NOW THE DEFAULT
RATHER THAN A STAGE-METHOD CARVE-OUT.  Forming `M` by `n` matvecs
and taking its exact spectrum has no threshold, no basis size and
no selection ambiguity; the truncated Arnoldi below has all three.

It used to run only for `dirk`/`full`, on the argument quoted
below -- and that argument was never stage-specific.  MEASURED
2026-09-07 on `_osc_with_ladder(16, 14, nslow)` (`gear`, so the
Arnoldi path), against this same dense spectrum, scored in the GAP
because `Q ~ 1/(1 - lam2)`::

    nslow   dense lam2     k=12 Arnoldi   gap ratio   Q dense/Arn
     <=11   0.995706203    0.995706197      1.000       232 / 232
       12   0.996324417    1.000114048     -0.031       271 / inf
       13   0.996818781    0.942674586     18.020       313 / 16.9
       14   0.997220139    0.999318472      0.245       359 / 1467

Not one-signed, so the Cauchy lower bound recorded above does not
hold on a circuit monodromy (it is stated for a NORMAL `M`, and
this is not one); and at `nslow = 12` it reports `lam2 > 1`, a
spurious UNSTABLE multiplier, which `Q` turns into `inf`.

⚠ RAISING `PPV_RITZ_BASIS` IS NOT THE FIX AND WAS MEASURED NOT TO
BE: `k = 16` is exact on that fixture and fails on a longer
ladder (20/20 -> 0.279, 26/26 -> 0.410), because the basis has to
grow with the problem.  A constant cannot.

⚠ `FLOQUET_DENSE_LIMIT` is the same cap `floquet_modes` applies to
the same assembly, so the two agree about what "affordable" means.
`dirk`/`full` keep the dense route ABOVE it as well: there it is
expensive, but the alternative is not slower, it is WRONG, and
those paths have never had the truncated one.

The comment on the stage map:

⚠ THE STAGE MAP IS DENSE AND WIDTH `m`, so its exact spectrum
is cheap -- and the Arnoldi below resolves it BADLY here.
`I - M` has `M`'s annihilated modes clustered at eigenvalue 1
and the physical unit root also at 1 after `1 - theta`;
measured, the Arnoldi left the unit root at `1 - 1.5e-6`, past
the `1e-6` deflation, so it reported the ORBIT TANGENT as the
second multiplier and `Q ~ 6e5`.  Forming `M` by `m` matvecs
and taking its eigenvalues directly gives the unit root to
machine precision (it deflates cleanly) and the true second
multiplier -- 8.59e-4 on van der Pol, matching Gear-2's
8.58e-4.

The comment on an empty window:

⚠ NO MULTIPLIER IN THE WINDOW: for a MULTISTEP or trapezoidal
solve the phase multiplier is 1 only as far as the
discretisation is time-translation invariant, and a
non-uniform `grid=` breaks that at O(h^2) (measured
1 - 5.1e-05 gear / 1 + 4.7e-05 trap at N = 400, 3:1, rate
4).  ⚠ NOT radau: its collocation solve keeps the
multiplier at 1 + 1.2e-11 on the same grid (2026-09-20).  Without
this the phase multiplier itself was reported as the
SECOND one -- silently, f_amp 600x too small.  Drop the
one nearest 1 instead; a uniform grid never gets here.

The comment on the truncated path:

⚠⚠ THE TRUNCATED PATH, NOW GATED ON THE PAIR'S OWN RITZ
RESIDUAL AND GROWN UNTIL IT CERTIFIES.  This branch runs only
where the dense spectrum is unaffordable -- which is exactly
where the truncation is least trustworthy, since a big circuit
is the one likely to carry the many slow nodes that break the
selection.  A fixed basis cannot work here: `k` has to track
the slow-mode count, so `PPV_RITZ_BASIS = 16` was measured to
be exact on one ladder and wrong on a longer one.  Doubling
until the residual certifies is the same rule at every size.

⚠ THE LOOP TERMINATES ON THREE THINGS and only one of them is
a threshold: the residual certifying, the basis reaching `n`
(where the Arnoldi IS the spectrum), or the cost ceiling
`PPV_RITZ_MAX_BASIS` -- which produces a WARNING and an
uncertified number, never a silently wrong one.

The comment on `Q`:

⚠ ONE NUMBER THAT SUBSUMES FOUR DIAGNOSTICS.  An amplitude
perturbation decays to `|lambda_2|` of its size each cycle, so
the cycles needed to fall below a threshold IS the oscillator's
Q: `Q = log(threshold)/log|lambda_2|` (Wang & Roychowdhury).
The usual definitions do not apply to an autonomous circuit --
`f_r/df` presumes a Bode plot of a BIBO-stable linear system,
and stored/dissipated presumes damping a self-sustaining
oscillator does not have.  Nor is it the resonator's Q.

⚠ AND IT IS THE SAME CONDITION AS EVERY FAILURE THIS CLASS
WARNS ABOUT.  "High Q", "a second multiplier near 1", "slow
amplitude restoration" and "a long time constant" are four
vocabularies for one thing -- which is why the same circuits
defeat the phase row, the eigen-split, the probe's continuation
and the PPV's instantaneous-response assumption.  Not four
coincidences.  It costs nothing here: the Arnoldi above already
produced `|lambda_2|`.

⚠⚠ AND ITS NAME IS ONLY RIGHT WHILE THE OSCILLATOR'S AMPLITUDE
MODE IS THE SLOWEST NON-UNIT MODE.  A parasitic with
`tau_p/T > Q_osc` simply IS the second multiplier -- by
definition, not by error -- and then this reports THE
PARASITIC'S DECAY TIME IN PERIODS under the name `Q`.  MEASURED:
at `tau_p/T` = 32 and 100 on a `Q = 16` oscillator, `lam2` is
the parasitic and `Q` returns 32 and 100.

⚠ THE NUMBER IS RIGHT AND ITS NAME IS WRONG, which is why
nothing misbehaves: every residual stays clean and the value is
well converged.  A DCO's gated capacitor sits at
`tau_p/T ~ 1e4`, i.e. permanently in that regime, so on exactly
the circuits a hierarchical DCO method exists for, a reported
`Q` would be the gated cap's RC in periods.  Read `Q` as
"cycles for the SLOWEST NON-UNIT MODE to decay by 1/e", which is
what it computes; it is the oscillator's Q only when that mode
is the oscillator's.

Reported for a `1/e` threshold, so `Q` is cycles-to-1/e.

⚠⚠ AND THIS LINE IS WANG & ROYCHOWDHURY'S IDENTITY
`Q = log(threshold)/log|lambda_2|`, which does double duty and
was shipped before either use was noticed.  It is what makes
"bounded by Q" and "bounded by lambda_2" the SAME SENTENCE --
the organising fact of this whole area, since a designer's
objective (raise Q) IS the numerics' failure mode (lambda_2 ->
1).  It is also an ERROR AMPLIFIER:

    (dQ/Q) / (dlambda_2/lambda_2)  =  -1/ln(lambda_2)  =  Q

⚠ SO THE RELATIVE ERROR IN `Q` IS `Q` TIMES THE RELATIVE ERROR
IN `lambda_2`, and a caller reading `Q` at high Q is reading a
quantity far less accurate than the multiplier behind it.
MEASURED end to end on van der Pol tuned by `mu = 1/(2 pi Q)`,
against the finest grid:

    Q      npts   rel err lam2   rel err Q   ratio
     3.18   120    1.04e-03      3.31e-03      3.2
    15.92   120    5.75e-04      9.23e-03     16.1
    63.66   120    4.89e-04      3.21e-02     65.7
    63.66   480    7.56e-06      4.81e-04     63.7

⚠ THAT IS A RESOLUTION REQUIREMENT SCALING WITH `Q`, NOT A
FIXED ACCURACY: 120 points/period gives `Q` to 0.3% at Q = 3
and only 3.2% at Q = 64.  Payable here because Gear-2's
`lambda_2` converges at better than second order (~8x per
doubling); a method that BIASES `lambda_2` at fixed order has
no such escape, and backward Euler's 5.6e-2 bias would become
85% in `Q` at Q = 100.

The comment on `info['null_residual_amplification']`:

⚠ MULTIPLY `null_residual` BY THIS TO GET THE RELATIVE
ERROR IN `v` THE RESIDUAL CANNOT EXCLUDE.  `null_residual`
is `||v - M^T v|| / ||v||`, so an error component along the
`lam2` left-eigendirection enters it scaled by `1 - lam2`
and is nearly INVISIBLE exactly when `lam2 -> 1`.
MEASURED on `_vdp_with_slow_node`, injecting a 1% error
into a converged `v` (floor 4.6e-11):

    lam2        r(random dir)   r(lam2 dir)   0.01*(1-lam2)
    0.000856      1.65e-02       1.003e-02      9.99e-03
    0.990049      1.65e-02       9.950e-05      9.95e-05
    0.999900      1.65e-02       1.000e-06      1.00e-06
    0.999999      1.65e-02       1.000e-08      1.00e-08

Exact to every digit printed.  A RANDOM error is caught
nine orders above the floor, so `null_residual` is a real
gate and this module's assertions on it can fail -- but it
loses sensitivity in the ONE direction that matters as the
circuit gets better, which is the opposite of the
reassurance a flat residual gives.

⚠⚠ THIS IS WHY A FLAT `null_residual` IS NOT EVIDENCE OF
ACCURACY.  A residual that does not move while `lam2`
sweeps toward 1 is not reporting that the answer stayed
good; the bordered system is well-conditioned BY
CONSTRUCTION, and the quantity it fails to see is
precisely the one that grows.  Read the two numbers
together or neither.

The comment on `info['period']`:

⚠ THE PERIOD OF THE ORBIT THESE SAMPLES LIVE ON, which is
not the caller's `period` when this came from a twin (trap
and euler read a TR-BDF2 twin whose period differs by
O(h^2)).  A quadrature over `times` divides by THIS; mixing
it with the host's period was the E3 "sign change".

### `frequency_aware_ppv`

The docstring before the move:

The PPV at a nonzero modulation frequency (Lai 2008, eq. 23).

⚠ EQ. 23, NOT 24 (docs session, 2026-09-09): Lai's eq. 24 drops the
AC columns of the Toeplitz block and is justified only "if we are
only interested in the transfer functions when w_s is close to
DC"; this shooting form has no such truncation -- `I - exp(-j w_s T)
M^T` is the exact sampled LPTV adjoint at ANY offset, the monodromy
already carrying the full time variation -- so it is eq. 23 for
what the object IS, exact for the discretised system, and eq. 24
only for the DC-reduction sentence pinned below.

The classical PPV is the left null vector of `I - M^T`, bordered by
`q = C(0) xdot(0)`; it is the phase response to a perturbation that
is SLOW against every other Floquet mode.  This is the SAME
bordered system at `alpha = exp(-j w_s T)`,

    [[I - alpha M^T,  q], [q^T, 0]] [v; y] = [0; 1],

whose solution `v(w_s)` is the phase sensitivity to a perturbation
modulated at `w_s`: at `w_s = 0` it IS `ppv()` (pinned), and away
from it the AMPLITUDE mode admixes with weight
`(1 - alpha)/(1 - alpha mu_2)` -- zero at DC, rising ten-fold per
decade, cornering where `2 pi f_s T = 1 - mu_2` and flat above
(docs session, 2026-09-08, on the slow-node fixture; the corner
tracks the slow multiplier over two decades).  Verified at the
source: Lai's eq. (24) at `w_s = 0` "is the augmented PPV
extraction equation", verbatim; his construction is harmonic
balance, this is the shooting basis, and the two agree on what the
object is.

Returns `(v, info)`: `v` the pair-space anchor vector (complex);
`info['samples_pair']` the T-periodic envelope over the period,
`lambda_k exp(+j w_s t_k)` with `lambda` the transposed replay of
`v` -- the sideband rows' own convention, so its Fourier
coefficient at harmonic `k` is the phase transfer of a source band
at `k f0 + f_s`; `info['samples']` the same in state space
(`C^T v`), SECOND order in the step through `ppv()`'s own
pair-consistent propagation run on the complex anchor
(`_ppv_propagate`, lifted 2026-09-08; at w_s = 0 it is `ppv()`'s
`samples`); `info['admixture']` the norm fraction of `v(w_s)`
orthogonal to the DC PPV; `info['corner']` the predicted corner
`|1 - mu_2| / (2 pi T)` in Hz; `info['alpha']`; `info['ppv']` the
DC object's info.

⚠ WHAT IT IS FOR.  A source that reaches the phase through a slow
path (an RC leg, tau >> T) is filtered at its own corner, and the
DC PPV cannot see that; the harmonic sum built from THIS object's
coefficients carries the filter inside `c_k(w_s)` with no explicit
model of the path (A2, roadmap).  MEASURED (2026-09-08): the ratio
`sum_k |c_k(w_s)|^2 / sum_k |c_k(0)|^2` reproduces `pnoise`'s
`S_pm/(4 S_v)` for a source behind the slow node within 0.6 % to
r = 1e-2 and 2 % at 5e-2 (a = 0.4, loss 0.2, tau/T = 100), where the
DC sum with the filter by hand was 4 % off; at r = 0.1 the two
differ by -6 %, unchanged at twice the grid -- a gap between PM by
sideband quadrature and phase-mode projection, both 1e-3 of DC
there -- and put to a nonlinear Monte Carlo that instantiates
neither construction (2026-09-09; 4 seeds x 10 000 periods per
point, phase read two ways).  ⚠ A first reading at one asymmetry
(a = 0.25) assigned each construction to ONE phase definition
crosswise (crossing phase to S_pm at 1.027, demodulated phase to
this sum at 0.977, each +-1.8 %); the asymmetry sweep a = 0.25 /
0.12 / 0.05 / 0 REPLACED it.  Double ratios (Monte Carlo slow/core
over predicted slow/core, each estimator calibrated on the core):

    a      crossing vs S_pm / this sum   demod vs S_pm / this sum
    0.25   1.032 / 1.066  (+-1.8 %)      0.954 / 0.985
    0.12   1.093 / 1.121  (+-3.5 %)      0.974 / 0.999
    0.05   1.122 / 1.149  (+-4 %)        1.004 / 1.027
    0.00   1.134 / 1.162  (+-4 %)        1.027 / 1.052

⚠⚠ AT 16 SEEDS PER POINT (Andreas, same day; +-1.5-1.7 %):

    a      crossing vs S_pm / this sum   demod vs S_pm / this sum
    0.25   1.038 / 1.072                 0.957 / 0.989
    0.12   1.085 / 1.113                 0.979 / 1.004
    0.05   1.117 / 1.143                 1.011 / 1.035
    0.00   1.129 / 1.156                 1.036 / 1.061

The statistic is the SLOPE in a, not any one point (the docs
session's framing): demod vs this sum -0.27 +- 0.09 per unit a
(3.2 sigma), vs S_pm 3.6 sigma; crossing 4.2 / 4.7 sigma.  So the
ratio is NOT constant in a at ~3 sigma for the demodulated phase
and above 4 for the crossing: NEITHER construction describes
EITHER measured phase across the range.  The two constructions
track each other to 1 % over the sweep while both estimators --
a point sample at a crossing and an average over a period -- drift
together, in the same direction, against both.  The a = 0.25
agreement of the demodulated phase with this sum is where its
curve crosses the sum, not a match.  What the crossing carries
beyond that: two thirds of its excess is waveform content beyond
0.5 f0 from the carrier (an instantaneous crossing aliases the
additive noise a one-period demodulation cannot see), the rest
the demod's own boxcar loss (1/sinc^2 = 1.045 at the band centre)
plus a common amplitude-to-crossing gain of ~0.6 (band-limited
coherence).  The common drift of BOTH estimators against BOTH
linear constructions as a -> 0 is the open object; candidate, a
second-order amplitude-to-phase conversion the linear theory
cannot contain -- REFUTED the same evening (PSD/4: the slow
fixture scales linearly, 0.993 +- 0.013; grid doubling moves the
constructions < 0.7 %; band conventions identical).  RESOLVED by a
forward tone-transient route on the MC's own discretised system
(no adjoint, no sideband assembly): pnoise's S_pm agrees with the
forward LPTV PM sidebands to 1 % at both asymmetries, and the
estimators' a-dependent double ratios are REPRODUCED by that
deterministic linear route (demod 0.938 -> 1.003 against the MC's
0.957 -> 1.036).  The drift is the ESTIMATORS: a one-period
fundamental demodulation leaks the other harmonics' sidebands
through its boxcar (sinc(pi(1-r)) ~ 0.1 for the second harmonic's,
which is ~a), zero crossings convert every harmonic's; both read a
given sideband PM with a fixture-dependent gain.  This object and
S_pm are PM by quadrature of the FUNDAMENTAL'S sidebands; compare
them with that, not with a demodulated or crossing phase.  ⚠ AND
COMPARE BAND WITH BAND: a source behind the slow node has an
in-band spectrum that is not 1/r^2 (its slow/core ratio swings
1.16 -> 0.87 across 0.08-0.15 f0 at a = 0.25), so a band mean and a
point value differ by ~4 % there; a 2 % constant between the noisy
Monte Carlo and the deterministic route remains after that, within
the excursion-amplitude bound, unresolved.  ⚠ The
premise "a -> 0 makes the definitions coincide" was wrong: the
asymmetry removes even harmonics only, van der Pol's third stays
at 9.7 % of the fundamental, and the construction gap GROWS as
a -> 0 (core 1.010 -> 1.055); the sinusoidal limit is mu -> 0.  The
instrument hypothesis (spectrum analyser <-> this sum, time-
interval analyser <-> S_pm) is refuted in its crossing half.  The slow multiplier's own coefficient
(`mode_content[0]`) corners at 1.6e-3 f0 for tau/T = 100 with a
plateau of 2.45e-6 (the docs session's 2.29e-6), scaling as
T/tau.  ⚠ DO NOT GATE ON `|v|`: with
`q^T v = 1` the `1/(1 - alpha)` pole cancels between numerator
and denominator and the norm is frequency-flat by construction; a
one-percent orthogonal admixture moves it by 5e-5.  The change is a
DIRECTION -- read `admixture`, or the per-harmonic coefficients.

The comment on the equation-row samples:

⚠ THE EQUATION-ROW SAMPLES, which is what `CY` contracts against
(see `ppv()`'s `samples_eq`).  `_equation_row_ppv` casts to float,
so the complex samples go through it as real and imaginary parts:
the map is linear (measured to round-off) and commutes with the
per-step phase factor.  At offset 0 these equal `ppv()`'s
`samples_eq` to 7e-14.  ⚠ `times` above is truncated to the sample
count, ONE ENTRY SHORT of the orbit's grid -- a quadrature over it
drops the last step (0.4 % of `c` on an asymmetric orbit, measured);
integrate over `info['ppv']['times']` and `['period']` instead.

### `floquet_modes`

The DENSE paragraph as it stood until 2026-09-25, when the dominant modes
above `FLOQUET_DENSE_LIMIT` were built (`_floquet_modes_ritz`):

⚠ DENSE, AND REFUSED ABOVE `FLOQUET_DENSE_LIMIT`: `n` matvecs, then
`eig`.  The extension is an Arnoldi that keeps its Ritz VECTORS, and ⚠ it
converges to the physical mode LAST, worse as Q rises (Garcia, Romero &
Acha 2022): on `A = I - M` the physical `lam_2 -> 1` maps to the SMALLEST
`theta`, while the fast parasitic modes (`theta ~ 1`) resolve first -- a
separation `1/theta_2 ~ Q_lambda`.  A truncated run also "cannot compute
ALL the Floquet multipliers", which eq (22) requires.  (That is an
eigenvalue question; the GMRES iterations of the bordered SOLVE are
Q-independent.)

(Built as an Arnoldi on `M`, not `I - M`, for exactly that reason; the
"ALL the multipliers" point is why `nmodes=None` still refuses above the
limit.)

The docstring before the move:

⚠ THE MODES' ACCURACY IS THE METHOD'S -- AND GEAR'S ADJOINT MODES
WERE FIRST ORDER UNTIL 2026-09-20, from a second-order method, on
every grid.  The reconstruction of `q` from the discrete adjoint took
the pair's first block through pinv(C^T), a fraction of a step off
the node (see the note at the reconstruction below).  The invariant
`q^T C p`, constant along the orbit for the true adjoint, measured on
the asymmetric van der Pol at N = 200 / 400 / 800:

    gear, uniform, before   9.9e-03  5.0e-03  2.5e-03   halving: FIRST order
    gear, uniform, now      9.3e-04  2.3e-04  5.6e-05   quartering: SECOND
    gear, 3:1 grid, now     2.2e-02  1.1e-02  5.9e-03   halving (was 3.4e-02):
                            a variable-step multistep adjoint is first
                            order and no rescaling makes it more

On a uniform grid at N = 400 gear's modal-spectrum PARTS now agree
with trap's (second order) to 4e-4 / 9e-4 / 5e-4 and the total closes
on pnoise at 1.0004 (it read 1.0145 before).  ⚠ A mode's k = 1
Fourier coefficient still differs from a uniform N = 3200 reference
by 3e-2 at N = 200 -- for `p` AND `q` alike, so it is the mode's
phase across N, not the adjoint; the invariant and the spectra are
phase-insensitive and are the gates.

Radau: on a 3:1 non-uniform grid its modes reproduce the uniform
grid's to 1e-10 and its phase multiplier stays at 1 + 1e-11; gear and
trap leave it at O(h^2) there and the modal spectra refuse.  For the
modes on a non-uniform grid use radau (the default).

The Floquet pairs `(λ_l, μ_l, p_l(t), q_l(t))` — A9's prerequisite.

Returns a list of dicts, one per mode, ordered by `|λ|` descending.
`nmodes=None` returns EVERY non-null mode, and that default is the
requirement rather than a convenience:

⚠⚠ ALL OF THEM ARE REQUIRED, BY THE SOURCE. Traversa & Bonani, IET
CDS 2011: "The calculation of orbital fluctuations and of the
phase-orbital correlation within Floquet-based noise analysis of
autonomous systems requires the availability of ALL the direct and
adjoint Floquet eigenvectors associated with the noiseless limit
cycle."  (Cited, not verified here; relayed from the paper.)  An
earlier cost estimate for A9 -- "a few more Floquet pairs" -- was
relayed without checking it against that sentence, and is wrong.

⚠ TRUNCATION IS LEGITIMATE ONLY WITH A BOUND. Traversa & Bonani,
TCAD 2013, compute a CHOSEN number of exponents and both
eigenvector sets for the linearisation of index-1 DAEs around a
limit cycle -- this formulation -- with the error "proved to tend
to zero along with the ratio between the norms of the NEGLECTED
AND RETAINED ROWS".  So passing `nmodes` is allowed, but a caller
who does owes that ratio as the gate; this routine does not
compute it.  The dense route below returns everything anyway, so
at the sizes it serves the question does not arise.


    lam    the Floquet MULTIPLIER, eigenvalue of the monodromy
    mu     the Floquet EXPONENT, `log(λ)/T` (complex)
    u0,v0  right and left eigenvectors at `t = 0`, biorthonormal
           (`v_k† u_l = δ_kl`)
    p      `p_l(t_j) = Φ(t_j,0) u_l(0) · exp(−μ_l t_j)` — the
           T-PERIODIC part, sampled on the PSS grid
    q      the adjoint counterpart from the reverse replay
    times  the grid `p` and `q` are sampled on

⚠⚠ **WHY THIS EXISTS: `S_yy` NEEDS THE EIGENVECTORS OVER THE
PERIOD, NOT JUST THE EXPONENTS.** Traversa & Bonani (TCAS-I 2011)
Lemma 3.5 makes the orbital spectrum a sum of Lorentzians centred
at `jω₀ + Im{μ_l}` with half-width `|Re{μ_l}| + ½h²ω₀²c`, weighted
by `C_lhj` (their eq 22) — and `C_lhj` is built from the FOURIER
COEFFICIENTS of `u_l(t)` and of `v_l(t)ᵀ B(t)`. Their own text is
explicit that the exponents alone do not order the result: *"a
major role in the C and D coefficients is also played by the
Floquet eigenvectors, which could determine large orbital
fluctuations contributions even when the Floquet exponents are not
near to zero."* So `|λ₂|` — the only mode information this class
used to expose — is not sufficient, by the source's own statement.

⚠ THE PERIODIC PART IS THE OUTPUT, NOT `Φ(t,0)u(0)`. Floquet's
theorem says the solution is `p_l(t)exp(μ_l t)` with `p_l`
T-periodic; the raw propagated vector is not periodic and its
Fourier series is not the one eq (22) wants. Dividing out
`exp(μ_l t)` is what makes `p_l(T) = p_l(0)` — which is also the
gate below, and the only check here that needs no reference.

⚠ DENSE, AND REFUSED ABOVE `FLOQUET_DENSE_LIMIT`. The monodromy is
assembled column by column (`n` matvecs) and diagonalised. That is
honest for the sizes this is useful at and wrong to hide at larger
ones: an Arnoldi route would return Ritz VECTORS rather than only
the Ritz values `ppv()` currently keeps, and is the extension.

⚠⚠ AND THAT EXTENSION CONVERGES TO THE PHYSICAL MODE *LAST*, WORSE
AS Q RISES -- a structural fact, not a measurement (Garcia, Romero
& Acha 2022, read firsthand by the docs session).  Arnoldi resolves
the LARGEST-magnitude eigenvalues first; the Ritz route works on
`A = I - M` and recovers `lam = 1 - theta`, so the physical
`lam_2 -> 1` maps to `theta_2 -> 0`, the SMALLEST, while the fast
parasitic modes (`lam ~ 0`) sit at `theta ~ 1` and are resolved
first.  The separation to resolve is `1/theta_2 ~ Q_lambda`: 3.7,
16.4, 64.5, 128.5 at `Q_lambda` = 3.18, 15.9, 64, 128.  The
difficulty scales with the very quantity being measured.  ⚠ This
is a DIFFERENT Krylov problem from B13's, which measured GMRES
iterations for the bordered SOLVE `(I - M) w = b` and found them
independent of Q -- solving a system and extracting its smallest
eigenvalue are not the same question, and B13 says nothing about
the second.  The paper is a sound source for the method and was
validated on power networks, not RF oscillators, so it reports no
evidence either way about the high-Q regime; and it states that a
truncated run "cannot compute ALL the Floquet multipliers" -- which
is the requirement eq (22) carries (IET CDS 2011, above).

The comment on `pss_unused`:

⚠⚠ `pss_unused` IS IGNORED, AS ITS NAME SAYS -- and it used to be
DEREFERENCED here, so the documented default call `floquet_modes()`
raised `AttributeError: 'NoneType' object has no attribute
'factored_period'` on every method (measured on radau, trap and a
GLM alike).  Every call site inside this file passes `self`, so
reading `self` changes no existing answer and makes the no-argument
call work; the parameter stays in the signature because callers
pass it positionally.

The comment on the staged solve:

on a staged solve the map is the TOTAL one (2026-09-22): the
crossings move with the state, `M + P_theta dtheta/dx_0`

The comment on null modes:

⚠ NULL MODES ARE DROPPED, NOT RETURNED WITH A BAD RESIDUAL. A
DAE's monodromy has exact zeros (the algebraic directions the
step map annihilates); their "eigenvectors" are arbitrary, the
exponent `log(0)` does not exist, and the periodic part comes
back as noise -- measured, residual 0.56 and periodicity 8.9
against 1e-15 for the physical pair. Returning them invites a
caller to average over a mode that means nothing.

The comment on truncation by magnitude:

⚠⚠ TRUNCATING BY MULTIPLIER MAGNITUDE IS REFUTED BY THE SOURCE'S
OWN WORKED EXAMPLE.  Traversa & Bonani TCAS-I 2011 Sec. V, on
their Colpitts: "six orders of magnitude separate mu_2 and mu_3,
while the corresponding contribution to orbital noise are not in
the same ratio.  Rather, far from the oscillator harmonics, the
contribution of mu_3 is dominant with respect to mu_2".  The
ordering INVERTS.  Four statements agree: eq (8)'s sum over
k = 2..n (structural), p.4 (asserted), Sec. V (measured on a
real circuit), this repo's concentration sweep (m/n = 0.97).
The caller who truncates owes the dropped weight as a gate;
this says so at the call rather than only in the docstring.

The comment on `nmodes=None`:

⚠ `None` means ALL non-null modes -- the default since the IET CDS
2011 correction -- and it used to fall into `int(None)` here because
the only test passed a number. A default nobody exercises is not a
default.

The comment on the adjoint replay:

adjoint: Phi(T,s_j)^T v_k(T) -- B8 made this available under
every integrator, not only the solved-history one
⚠⚠ THE ADJOINT IS THE PER-STEP TRANSPOSED SOLVE `t`, NOT THE
PAIR'S FIRST BLOCK -- AND IT BELONGS TO THE NEXT NODE
(2026-09-20, measured).  `collect` hands back both: `ts[k]`,
the solve `Jf_k^-T w1` made while replaying step k backwards,
and `states[k]`, the pair (w1; w2) it leaves behind.  This
took the pair's first block and mapped it through pinv(C^T);
since `w1 = Jf^T t = (a0 C + G)^T t` and the adjoint equation
`C^T dq/dt = G^T q` turns the `G^T t` part into a time
derivative, that `q` was `a0 * q(t + 2h/3)` -- staggered by a
fraction of a step, so the invariant `q^T C p` drifted along
the orbit by 1e-2 at N = 200 and HALVED per doubling: FIRST
order, from a second-order method, on every grid.  The solve
`t` itself obeys the BDF2-discretised adjoint recursion, and
the replay computes it for step k from the pair at node k+1,
so `ts[k]` is the adjoint at node k + 1.  Scored on the
invariant's spread at N = 200 / 400 / 800, uniform grid:

    pair block, pinv(C^T)  (this, before)   9.9e-03 5.0e-03 2.5e-03   x2 per doubling
    ts[k] at node k        (one node off)   1.5e-02 7.6e-03 3.8e-03   x2
    ts[k] at node k + 1    (this, now)      9.3e-04 2.3e-04 5.6e-05   x4  SECOND ORDER

⚠ THE SCALE.  `t = Jf^-T w1` carries the step through `a0 ~ 1/h`,
invisible on a uniform grid (absorbed by `c0` below) and a
factor-3 modulation on a 3:1 one, so `a0` of the node's own step
is put back.  ⚠ THAT IS FIRST ORDER ON A NON-UNIFORM GRID and
cannot be more: the discrete adjoint of a variable-step
multistep method draws `a1` and `a2` from LATER steps, so its
recursion is the continuous adjoint's only to O(h) once the step
changes (Sandu's inconsistency).  Four scalings were measured
on the 3:1 grid and all halve per doubling; this one is the
best of them at 2x the previous accuracy.  A non-uniform gear
grid is refused by the modal spectra anyway (its phase
multiplier leaves the unit circle); radau is exact there.
Node 0 is node N by periodicity of the periodic part.
⚠ GEAR ONLY (`solved_history`).  A one-step kind's `ts` is
NESTED -- per-stage solves per step -- and its state-block
adjoint through pinv(C^T) was measured exact (radau, 1e-10 on
a 3:1 grid) and second order (trap); those keep their path.

The comment on `C^T q`:

⚠⚠⚠ THE REPLAYED VECTOR IS `C^T q`, NOT `q`.  The conserved
bilinear form of the variational DAE is `w^T C delta`, so over
a period `M_a^T C M = C`, which makes the LEFT eigenvector of
the state monodromy `C^T w(0)` -- the adjoint mode in the
"left-eigenvector coordinates", one factor of `C^T` away from
the state-space adjoint `q` that eq (22) and every covariance
here need.  The transposed replay propagates that object, so
every sample of `q` above is `C(t)^T q_true(t)`.

⚠ INVISIBLE ON EVERY FIXTURE THIS REPO HAD, for a geometric
reason: van der Pol's reduced `C` is `diag(1, -1)`, and at
`t = 0` the orbit sits at `[2, 0]` where the adjoint is nearly
axis-aligned, so `C^T q` and `q` point the same way up to sign
(`|cos| = 0.9972`).  On an ASYMMETRIC orbit the seed is off-axis
and the two separate -- measured `|cos(v_k, q_true)| = 0.5738`
at `a = 0.30` on van der Pol + `a u^2` -- while the two adjoints
there are nearly PARALLEL (`|cos(q_2, q_1)| = 0.997`), so the
wrong vector is mostly phase adjoint.  Result: the orbital
covariance was 81x LOW against a Monte Carlo (0.0123 of the
truth), and `|cos(C^-T v_k, q_true)| = 1.0000` at both
asymmetries.  Applying `C^-T` here takes it to 1.06 of the
Monte Carlo at `a = 0.30` and 1.0004 at `a = 0`.

⚠ THIS ALSO DISSOLVES THE "q IS STORED IN REVERSE TIME"
finding recorded the same day: with the right vector,
`q(t)^T C p(t)` is conserved at the SAME index (4.2e-04) and
NOT the reversed one (2.0).  `diag(1, -1)` flips one
component, which on a half-wave symmetric orbit is exactly the
relation between `q(t)` and `q(T - t)` -- a sign flip read as a
time reversal.

⚠ Per sample, because `C` may depend on the state.  `pinv`
rather than `inv` so a singular reduced `C` (an index-2 MNA,
algebraic rows) does not raise; the algebraic components of `q`
are then the minimum-norm choice, which is a SCOPE LIMIT and
not a solution -- recorded, not hidden.
The pinv(C^T) map belongs to the STATE-BLOCK adjoint of the
one-step kinds; gear's transposed solve is already the adjoint
of the DAE variable (its invariant is `q^T C p`, see above).

The comment on the renormalisation:

⚠⚠ RENORMALISE ON THE STATE BLOCK. `v_k` was biorthonormalised
against `u_k` at the map's FULL width `n`; under a
solved-history map that is the pair `[x_n; x_{n-1}]`, and the
width-`m` state block then carries `q(0)^T p(0) = c0 != 1`.
MEASURED on van der Pol under gear: c0 = 1.324143, constant
around the cycle to four digits -- and the orbital covariance
assembled from these parts came out too large by EXACTLY
c0^2 = 1.7535 against two independent routes, because `q`
enters it quadratically. The periodicity gate p(T) = p(0)
cannot see this: periodicity is scale-free. On the plain path
n = m and c0 = 1, so this is a no-op there. The adjoint takes
the scale (the right vector is the physical direction).
⚠⚠ THE INNER PRODUCT IS `C`-WEIGHTED, AND THE UNWEIGHTED ONE
WAS WRONG BY A FACTOR OF `C` -- INVISIBLE ON EVERY FIXTURE
THIS REPO HAD.  The conserved bilinear form of the variational
DAE is `q(t)^T C(t) p(t)`, not `q(t)^T p(t)`: differentiating
`G p + d(C p)/dt = 0` against the adjoint gives
`d/dt [q^T C p] = 0`, so `q^T C p` is the invariant and the
biorthonormality that eq (22) assumes is `q_k^T C p_l = d_kl`.

⚠ ON A UNIT-REACTANCE FIXTURE THE TWO ARE THE SAME NUMBER,
which is exactly why this survived: van der Pol with
`c = L = 1` gives `q^T C p = 0.9992` against `q^T p = 1`.
Sweep the capacitance at fixed `w0` and the two separate --
MEASURED `q^T C p` = 0.2495 / 0.9992 / 3.9982 at
`C` = 0.25 / 1 / 4, i.e. exactly `C`, while `q^T p` stayed
pinned at 1.000000.

⚠⚠ AND `q` ENTERS THE COVARIANCE QUADRATICALLY, so the orbital
covariance came out too large by exactly `C^2`.  Measured
against the independent Lyapunov reference before the fix:
ratio 0.0624 / 1.0001 / 16.043 at those same `C` -- right ONLY
at `C = 1`, which is the only place A9's three-way gate ever
ran.  §D 0c, on the very circuit that produced that entry: a
unit reactance makes `C` the identity and the two inner
products indistinguishable.

The pair-slicing correction this block was written for is
subsumed: normalising on `q^T C p` fixes the slice scale and
the weighting in one step.

### `_continuous_adjoint`

The docstring before the move:

A mode's adjoint `q(t_j)` on a NON-UNIFORM gear grid, by integrating
the continuous adjoint equation SEPARATELY (2026-09-20).

The exact transpose of gear's two-step recursion draws `a1` and `a2`
from LATER steps, so on a grid whose step changes it is a consistent
scheme for the adjoint equation only to first order, and no per-node
rescaling lifts it (four measured, all halving per doubling).  Here
`C^T dq/dt = G^T q` is integrated backwards with BDF2 whose
coefficients belong to the REVERSE grid's own step pair --
`companion_coefficients(h_{n+1}, h_{n+2})`:

    (a0' C_n + G_n)^T q_n = -C_n^T (a1' q_{n+1} + a2' q_{n+2})

The state is the pair (q_{n+1}, q_{n+2}); the backward map over one
period is built from 2m basis propagations and the mode's eigenvector
matched to the forward multiplier `lam` (they agree to O(h^2):
1.2e-3 / 3.1e-4 / 7.7e-5 at N = 200 / 400 / 800).  MEASURED on the
asymmetric van der Pol, 3:1 grid, invariant `q^T C p` spread:

    transpose, a0-scaled   2.2e-02  1.1e-02  5.9e-03   halving
    this                   1.8e-03  4.6e-04  1.1e-04   QUARTERING

and `orbital_correlation`'s R against radau's (exact there) 2.3e-02 /
7.2e-03 / 2.0e-03, monotone at ~x3.5, onto the floor `p` sets.

⚠ TWO ADJOINTS IN THE TREE, ON PURPOSE.  This one is NOT the transpose
of the discrete map: the PPV, the adjoint noise folds and the sideband
rows keep the exact transpose (pinned at 1e-15) because their
identities need it.  This serves the MODE SHAPES only, where the
continuous adjoint is the object wanted.  On a uniform grid the two
coincide (`a0' = a0`) and the transpose is used; one-step kinds never
come here.  Dense 2m x 2m up to `CONTINUOUS_ADJOINT_DENSE_M`
unknowns; above that, Arnoldi on the backward map with Ritz-residual
certification (see the note at the branch).  ⚠ Biorthonormality against the OTHER modes is not
exact here (it is for the transpose) -- each mode is normalised on its
own `q(0)^T C p(0)`, as before.

The comment on the matrix-free branch:

⚠ MATRIX-FREE ABOVE THAT (2026-09-20, item 3 of "gear as a
first-class choice on non-uniform grids").  The wanted modes --
the phase mode and the slow orbital ones -- are the DOMINANT
eigenvalues of the backward map, so a plain Arnoldi on it with
the Ritz-residual certification `_ritz_second_multiplier` uses
(|h_{k+1,k}| |y_last| / ||y||) reaches them at a basis far below
2m: MEASURED equal to the dense eigenvector to cos 1.00000000 at
a basis of 3 / 16 / 24 for 2m = 4 / 32 / 124, residuals 1e-16 ..
1e-52, on the hostile fixture and the ladder oscillator
re-solved on a 3:1 grid.  Grown from `PPV_RITZ_BASIS` toward
`PPV_RITZ_MAX_BASIS` until the matched pair certifies; a pair
that never certifies is refused, not returned.


## `_pss_accuracy.py` -- `_AccuracyChecks`

### `monodromy_twin`

The docstring paragraphs before the move (the B16 paragraph stood twice):

⚠⚠ THE TWIN DEFAULTS TO TR-BDF2, MEASURED (2026-09-05), and Gear-2
is one setting away, not retired.  The twin exists because trap's
and euler's own monodromy is unusable on a limit cycle: trap's
diverges with refinement (B16, below) and euler's is first order.
Both TR-BDF2 and Gear-2 give a clean second-order twin, but TR-BDF2
is more accurate on `lambda2` -- measured against exact references
(`exp(A T)` on a linear oscillator; Abel's `exp(mu integral(1-v^2))`
on van der Pol) it is 12-32x better at practical step counts, and
the advantage GROWS with Q and with coarser grids -- the regime a
real oscillator PSS sits in.  At Q=100 and 50 points/period the
Gear-2 twin misreads Q by 22%, the TR-BDF2 twin by 0.27%.  The gap
is a coarse-grid/high-Q effect, not a fixed factor: refine the grid
or drop Q and both fall to the ordinary O(h^2) floor where the
difference is single digits and can even favour Gear-2.  So
`monodromy = 'trbdf2'` is the default, `'gear'` restores the former
twin, `'native'` reads the run's own.

⚠⚠ THE B16 DECISION, TAKEN ON A MEASUREMENT (2026-09-05).  On the
bias-sensitive oscillator (`vdp + 0.3 u^2`, exact `Q_lambda =
5.908`, `c_true = 5.3703e-06`) trapezoidal's monodromy is unusable
with EITHER opener: the default reads `Q_lambda` 11.1 / 28.4 / 63.9
at 400/800/1600 points and DIVERGES with refinement, and
`x0_unknown=True` reads 3086 / 12228 / 48699 -- a spurious
multiplier at 1 (the one-step companion's parasitic mode), while
the state and period are second order either way.  So "the most
accurate" is not a choice between openers: the STATE keeps the
method you asked for, and every monodromy-derived quantity -- `Q`,
the PPV and everything built on it, the Floquet modes, the
phase-noise surfaces -- comes from the twin on the same orbit.  The
twin's period differs from this one's by O(h^2); its orbit is
re-converged, not copied.

⚠⚠ THE B16 DECISION, TAKEN ON A MEASUREMENT (2026-09-05).  On the
bias-sensitive oscillator (`vdp + 0.3 u^2`, exact `Q_lambda =
5.908`, `c_true = 5.3703e-06`) trapezoidal's monodromy is unusable
with EITHER opener: the default reads `Q_lambda` 11.1 / 28.4 / 63.9
at 400/800/1600 points and DIVERGES with refinement, and
`x0_unknown=True` reads 3086 / 12228 / 48699 -- a spurious
multiplier at 1 (the one-step companion's parasitic mode), while
the state and period are second order either way.  Gear-2 reads
5.9094 / 5.9086 / 5.9084.  So "the most accurate" is not a choice
between openers: the STATE keeps the method you asked for, and
every monodromy-derived quantity -- `Q`, the PPV and everything
built on it, the Floquet modes, the phase-noise surfaces -- comes
from Gear-2 on the same orbit.  The twin's period differs from
this one's by O(h^2); its orbit is re-converged, not copied.

The comment on the self-sufficient methods:

⚠ SELF-SUFFICIENT METHODS TAKE NO TWIN, and the method says which
it is (`carries_own_monodromy`): Gear-2 (a second-order native
companion monodromy) and every stage method (TR-BDF2, Radau -- no
opener seam, verified against the pencil).  This twin exists only
because a one-step LMM's monodromy is first-order on a limit
cycle (its manufactured opener is dropped to Euler and that seam
sits in the period map); twinning a self-sufficient method would
replace its own map with another on a re-converged orbit -- pure
cost -- and hide the run's own spectrum.  They read `native`
regardless of `monodromy`; the knob governs which twin trap/euler
borrow.

The twin could be trbdf2 or gear until 2026-09-24.  Now any method whose
own map serves (`carries_own_monodromy`) may be the twin.  Measured on van
der Pol, a trap run at 200 points against radau at 800 (lambda2 / c /
oscillator d): gear 1e-3 / 1e-4 / 2e-3; trbdf2 1e-4 / 6e-6 / 8e-5; radau
1e-10 / 2e-11 / 6e-13 at about trbdf2's cost; esdirk43 1e-9 / 2e-9 / 3e-5.

### `TWIN_MAXITER`

The comment on the constant:

The iteration budget of a monodromy twin's re-solve (Andreas,
2026-09-23).  A twin is a POLISH, not a cold solve: it is seeded at
this run's converged state on the same grid, and from a good seed it
converges in 4-13 iterations (measured on B16's van der Pol).  It used
to inherit the caller's `maxiterations` -- 300 in B16 -- and from a
poor seed (Euler at 400 points, its orbit 55 % off) both twins refuse by
NOT converging, so each ran its whole budget times the solve's retry
ladder: 971 + 968 traversals, 790 s, to say "no".  Capped at this, and
a twin that hits the cap without converging WARNS before it refuses.

### `_solve_twin`

The comment on the re-solve:

⚠ THE STEP COUNT MUST SURVIVE `solve`'s `int(period / timestep)`.
`T / (T / N)` is not N in floating point: measured on B16's
fixture, T = 6.730731946457316 gives 399.99999999999994, so the
"same grid" twin ran on 399 points against the state's 400 (found
when `phase_rule='reselect'` moved T in its 15th digit).  Half a
step of slack makes the floor land on N for every T.

⚠ AND THE TWIN TAKES THE DEFAULT PHASE RULE RATHER THAN THIS RUN'S.
It is seeded at the converged state so that it lands on the SAME
point of the SAME orbit -- which is exactly what the agreement
check below tests -- and `phase_rule='reselect'` re-chooses the
pinned coordinate and lands on ANOTHER phase (measured; see
`solve`).  Inheriting it would move the twin off the orbit point
whose monodromy was asked for.

The comment on the agreement check:

⚠⚠ THE TWIN MUST HAVE CONVERGED TO THE SAME ORBIT IT WAS SEEDED
ON, and a MORE ROBUST twin makes this check load-bearing rather
than paranoid.  Measured: from a poor seed (euler at 400 pts, its
orbit 55% off) the Gear-2 twin fails to converge -- LOUD -- but
the TR-BDF2 twin, being more robust, CONVERGES to a SPURIOUS limit
cycle and reports `Q = 1.97` against the exact 5.91 with no error.
Improving the method degraded safety: the failure moved from a
refusal to a plausible wrong number.  So the twin's converged
orbit is checked against the seed it was handed: two convergent
methods on the SAME limit cycle agree on period and entering state
to O(h^p) -- measured 6e-6 / 1e-4 (trbdf2) and 3e-5 / 2e-2 (gear)
on a good seed -- while the spurious jump above sits at 2.17 /
0.91.  The 0.25 gate is ~12x above the worst good case and ~3.6x
below the spurious one; it is set from the (universal, tiny)
good-case agreement, not the (fixture-dependent) failure size, so
it transfers.  The definitive test is refinement (a spurious orbit
does not survive h/2); this cheap consistency check is the
conservative stand-in -- it REFUSES a too-poor seed rather than
risk trusting it, which is the safe direction.

### `METHOD_ORDER`

The comment on the table:

Nominal convergence order per `method`, for `grid_error`'s ceiling on
a plausible OBSERVED order.  Sourced from this file's own measured
records rather than from the literature: trap/gear/theta second order,
TR-BDF2 measured at 4.01x/4.01x/4.00x per halving (exact `O(h^2)`),
Radau IIA(3) at 31.50x/31.74x (`O(h^5)`, theoretical 32), euler first.

Replaced on 2026-09-24 by `_nominal_order`, which reads the order off the
method's integrator (its `ORDER`, or a GLM's `order`).  The table covered 6
of the 12 accepted names.  esdirk43, glm2-4 and the aliases `gear2` /
`trapezoidal` took the generic range, so the ceiling was not applied.

### `grid_error`

The docstring paragraphs before the move:

⚠⚠ WHY THIS EXISTS RATHER THAN A PER-METHOD FORMULA.  The floor of
this stack is DISCRETISATION, it grows LINEARLY IN Q, and it is a
METHOD property: measured on the analytic high-Q van der Pol
reference, the relative error in the diffusion constant at 240
points per period is

    Q      gear        trap        radau
     100   1.79e-03    2.63e-05    6.97e-10
     500   9.02e-03    1.32e-04    3.48e-09
    1000   1.82e-02    2.63e-04    6.97e-09

(trap's column corrected 2026-09-14: it read 1.49e-06 / 1.04e-04 /
2.36e-04, the period-normalisation defect's values.)

i.e. `~1.8e-05 Q` for gear against `~7.0e-12 Q` for radau -- SIX
ORDERS at the same cost per step.  Those constants are real but they
belong to THAT fixture at THAT grid: `gear` converges at `O(h^3)` on
an autonomous problem for `Q >= 5` and at `O(h^2)` at `mu = 1`, so a
shipped predictor built from them would extrapolate a fitted constant
across a regime change (roadmap D.0y).  Refining the actual circuit
measures the actual number instead, and needs no calibration.

⚠⚠ WHY THREE GRIDS AND NOT TWO.  With `f_h = f + C h^p`, two grids
give `|f_h - f_h/r| = |C| h^p (1 - r^-p)`, which over-states the fine
grid's own error `|C|(h/r)^p` by `r^p - 1` -- an upper bound, and a
tempting place to stop.  **IT IS NOT SAFE, AND THIS STACK CONTAINS A
COUNTEREXAMPLE.**  A quantity mixing TWO discretisations -- an
`O(h^3)` error plus an `O(h^2)` one of opposite sign -- changes sign:
the two terms cancel, the two-grid difference collapses, and the
estimate UNDER-STATES the true error by 3.6x (measured: change
1.15e-06 against a true 4.06e-06 at 240 points).  ⚠ That quantity
was `diffusion_constant` under `trap` until 2026-09-14 -- the twin's
integral over trap's own period, a DEFECT now fixed -- and the
validity check below is what refused it.  A bound that fails silently
where the error is interesting is worse than none.

So the third grid is not extra confidence, it is the VALIDITY CHECK.
From `d1 = |f_h - f_h/r|` and `d2 = |f_h/r - f_h/r^2|`,

    order = log(d1/d2) / log(r)

is the order the circuit ACTUALLY shows, and it is checked against the
single-power-law assumption before the error estimate built on it is
offered.  Measured orders on that fixture: `gear` 2.94 (its `O(h^3)`
autonomous rate), `radau` ~5, `trap` 3.02 (its twin's), and the
mixed quantity failing the check exactly where it cancels.  `error` is then `d2 / (r^order - 1)`, and
`power_law=False` means READ `d2` AS A RAW CHANGE AND NOTHING MORE.

⚠ AND IT IS AN ESTIMATE OF THE GRID ERROR ONLY.  It cannot see an
error both grids share -- a wrong stamp, a wrong tolerance
convention, a mis-specified circuit.  A small `rel_change` says the
grid is fine enough; it does NOT say the answer is right.  That is
the same trap `null_residual_amplification` documents one screen up,
and it is worth stating twice.

The comment on the order ceiling:

⚠⚠ THE CEILING IS THE POINT OF THIS CHECK, AND A GENERIC
RANGE IS NOT ENOUGH.  A method cannot converge faster than
its order; an observed order well above it means two error
terms nearly cancelled at this grid, which makes the
deltas shrink faster than the error and the estimate
UNDER-state.  MEASURED: a quantity mixing two
discretisations (a twin's `c` over the host's period --
what `diffusion_constant` under `trap` computed until the
2026-09-14 fix) shows an apparent order of 6.45 at 120
points -- monotone, same-signed deltas, nothing else
suspicious -- while its estimate under-states the true
error by 300x.  A plain `0.5 <= order <= 8` range
ACCEPTS that case; the ceiling below rejects it.
⚠ The `+ 1.5` allowance is not slack: on an AUTONOMOUS
problem the period is an unknown that absorbs the leading
frequency error, so `gear` (nominal 2) genuinely converges
at 3.01 here.  Without the allowance this would reject the
shipped default method on its own reference fixture.

### `IDEC_DEGREE`

The comment on `WARPING_CHECK_TOL` and `IDEC_DEGREE`:

B7: the interpolant degree the defect-correction estimate needs, per
method.  Two-part rule, MEASURED 2026-09-08 (doc/pss_roadmap_260902.md,
B7's gate): the interpolant's DEGREE must exceed the method's stage
count -- a cubic spline lies INSIDE Radau IIA(3)'s collocation
exactness class (degree s = 3), so the neighbouring problem is solved
EXACTLY and the estimate is 1e-10 ppm against a true 6e-3 -- and its
ORDER must exceed the method's effective order, or a constant bias
remains (a quintic against radau's measured 6.1 left 4.9 % at every
grid; a septic gave 1.0001 / 0.9999 / 0.9998).  Cubic reproduced trap
and TR-BDF2 to 0.9996 -> 1.0000.
⚠ THE STAGE-COUNT CLAUSE IS A COLLOCATION PROPERTY (peer, measured the
same day): ESDIRK43 has SIX stages and is not a collocation method, so
a quintic fails the clause literally -- and quintic and septic AGREE
through the stack (0.9998 / 1.0000 at 100 pts, 1.0000 / 1.0001 at 200).
For a non-collocation method only the ORDER clause is established;
written as "degree > stage count" the rule would over-constrain every
DIRK ever added.  The failure the clause guards against is SILENT (a
clean small number), which is why it was measured rather than argued.

Replaced on 2026-09-24 by `_idec_degree`: the smallest odd degree above
the method's order, at least 3.  That rule reproduces every measured entry
here (cubic for the second-order methods, quintic for esdirk43, septic for
radau).  The table had given glm3 and glm4 a cubic.  Measured on van der
Pol: glm4 at 120 points read -6.7e-7 against a true +2.8e-9, and +5.7e-10
quintic.  Even degrees trip the half-grid check.

### `warping_estimate`

The docstring paragraphs before the move:

B7's answer, MEASURED (2026-09-08).  A per-step local truncation
estimate cannot see an accumulating period error because warping is a
GLOBAL error; defect correction (Sickenberger, Weinmueller & Winkler,
"Local Error Estimates for Moderately Smooth ODEs and DAEs", Part I,
Sec. 1) estimates the global error directly:

Continuing after the numbered steps and the transient remark:

Measured on A10's van der Pol (Q = 1e4), estimate / true period error
(true = T_h - T_ref, radau at 3200 points), numpy prototype:

    trap   / cubic    0.9996  0.9999  1.0000  1.0000   (100..800 pts)
    trbdf2 / cubic    0.9999  1.0000  1.0000
    radau  / cubic    0.0000  0.0000  0.0000   (25..50 pts) -- INSIDE the
                                               exactness class: see IDEC_DEGREE
    radau  / quintic  1.0490  1.0494  1.0495   -- a constant bias where the
                                               orders tie (6 vs 6.1)
    radau  / septic   1.0001  0.9999  0.9998

Controls: with RADAU solving the neighbouring problem of TRAP's defect
the estimate is 0.0000 -- the drift is the METHOD's error, not the
defect's; a LINEAR interpolant gives 0.03 / 2.3 / 3.4 -- the
interpolant-order wall from below.

Through THIS method (the stack, same fixture, reference radau at 3200
points, 2026-09-08): trap 1.0002 at 400 and 200 pts; radau septic
1.0003 / 1.0002 at 50 / 35 pts and CUBIC 0.0001 (the exactness-class
zero, reproduced); esdirk43 quintic 0.9998 / 1.0000 and septic
1.0000 / 1.0001 at 100 / 200 pts -- so for a non-collocation method
the order clause alone is established, and the stage-count clause is
a collocation property -- and Part I sec 1.1 says why: "one of the
most attractive features of the IDeC procedure is, that its fixed
point is a certain superconvergent COLLOCATION solution", so the
exactness class the stage-count clause guards against is a
collocation object by construction (docs session, 2026-09-09; one
family at two degrees on one fixture, so a mechanism, not a proof
that the clause is harmless in general).  ⚠ The first stack gate's driven control came
back `autonomous=True`: `Circuit.u(t)` evaluates its time functions
only when told `analysis='tran'`, and without it every source
VANISHES (zeros, DC value included) -- so that control ran against a
circuit with NO source at all; fixed at both call sites; with the
flag the driven van der
Pol returns `autonomous=False`, `period_error=None`, and a bounded
lag series, as it must.

⚠⚠ THE INTERPOLANT IS THE LIMIT (measured 2026-09-08): on a
relaxation oscillator with a comparator edge a few points wide the
estimate reads 0.09 of the true period error at 200 points per
period and 0.65 at 400 -- uniformly in every component, so not a
collapse: the septic spline does not resolve the edge and the
defect is interpolation error, not the method's.  Part I's own
scope is "moderately smooth"; a relaxation orbit at PSS grids is
outside it, and the number returned is then wrong by a factor that
nothing in it announces.  Trust it on smooth orbits (1.000 to four
digits on the van der Pol, index 1 and 2); on an orbit with edges,
refine until the estimate converges in `periods` and grid, or use
`grid_error`.
⚠ PART I READ THROUGH (docs session, 2026-09-09): its motivation is
this domain -- "we are especially motivated by applications in
electrical circuit simulation, where the models often contain data
with poor smoothness" -- so "moderately smooth" is the case the
paper was built for, not a clause this is outside of.  Its Remark
2.9 names a failure with the SAME SIGN as the edge underestimate:
the local estimates assume the leading term `c_i h^(p+1) x^(p+1)`
does not vanish, and "at least in case of oscillatory solutions,
there always exist time points where the derivative x^(p+1)
vanishes ... our error estimates will tend to UNDERESTIMATE the
true size of the error" (footnote: the third derivative vanishes
where the curvature is extremal; remedy: assume C^(p+2) and match
the next coefficient with an auxiliary scheme).  That remark is
stated for the LOCAL-error route of their section 2; this method is
the GLOBAL route of section 1 (Zadunaisky), and whether the global
route inherits it is not established.  ⚠ Ruled out here by the
h-scaling: Remark 2.9's mechanism is keyed on isolated zeros of
x^(p+1), whose aggregate effect is roughly h-INDEPENDENT, while the
edge reading improved 7.2x for a 2x grid (0.09 -> 0.65) -- that is
interpolation error, as stated above.  Where Remark 2.9 would bite
is a step controller built on defect correction; the paper hands
the fix.  Cost lead, not worked out: section 1.1's cheap variant
runs the high-order method once and a cheap LOW-order method twice
(original and neighbouring problem) -- a different substitution
from the radau-on-trap's-defect control that zeroed the estimate.
⚠ THE LITERATURE'S ANSWER IS STRUCTURAL, NOT "REFINE" (docs session,
Part I p. 9, READING-LOG 2.165).  The gate this instrument should
test before returning a number is Part I's own "only if": the
estimate is asymptotically correct ONLY IF the interpolant's defect
error is o(h^{p+1}) -- asymptotically SMALLER than the truncation
error it is meant to reveal; on the comparator edge it is not, and
the number is wrong by a factor nothing announces.  And the fix for
a non-smooth orbit is to form the defect as a WEIGHTED SUM OF
f-VALUES with an auxiliary scheme sharing the base scheme's
left-hand side, so the solution terms cancel identically (their eq.
2.13, an extra factor h) -- not a higher-degree interpolant of the
solution, which is exactly the construction this one uses.  Scope:
their construction is the LOCAL error of an LMM; whether it
transfers to a period functional is unproven.  THE GATE IS BUILT
(2026-09-08, `check=True`): the same pass through every second
sample of the same solution, transient still at the solve's step;
`check_ratio` = half-grid slope / full-grid slope, `trusted` =
within `WARPING_CHECK_TOL` (5 %) of 1, else a warning and the number
still returned.  Measured: van der Pol 1.0000 (radau and trap, 50
and 100 points); the relaxation orbit 0.0056 / 0.72 at 200 / 400
points under radau and 0.41 at 200 under trap -- the cases that
read 0.09 / 0.65 of the truth are refused, the smooth case accepted
with four orders of margin.  ⚠ The prediction "below 0.5 at 400"
was wrong (0.72): the ratio approaches 1 as the edge resolves, so
the tolerance is the gate, not the ratio's distance from 0.  The
restructured (f-value) defect (Part I eq. 2.13; for trap the Milne
device) was GATED and REFUTED as the edge fix (2026-09-08): smooth
0.9990 / 0.9998, but on the edge orbit 0.06 / 0.61 / 1.29 at 200 /
400 / 800 points against the spline route's 0.04 / 0.29 / 0.63 --
faster with the grid and NOT monotone, so a reading near 1 is
indistinguishable from a wrong one; the paper's own remedy is mesh
adaptation.  Not built.  Cost of the check: it doubles the call (a
second `periods`-long transient).
⚠ Scope and limits.  The period reading needs an AUTONOMOUS solve;
on a driven circuit the lag is bounded (entrained) and `period_error`
is returned as None with the per-period lag series still filled.
The prototype ran on a 2-state ODE; on a DAE the differential and
algebraic components converge at different orders (H&W VI.7), so the
interpolant threshold binds per component and a component-wise
exactness collapse would be invisible in this scalar phase drift --
`component_rms` is returned so a caller can look.  Index-2 is outside
Part I's stated scope.  Cost: `periods` periods of transient at the
working grid -- no refinement sweep, no analytic reference.

The comment on `analysis='tran'`:

⚠ `analysis='tran'`, on BOTH calls.  `Circuit.u(t)` evaluates a
time function only when told which analysis is asking (`VS.u`:
`elif analysis in timedomain_analyses`); without it every source
VANISHES -- the else-branch returns zeros, and even the DC value
lives inside the gated branch (`timedomain_analyses = ('dc',
'tran')`).  The first gate's driven control -- an `ISin` on the
van der Pol -- came back `autonomous=True` for exactly that reason,
and the defect would have omitted the drive on a driven circuit.

The comment on the per-component threshold:

a RELATIVE threshold: a node pinned by a source has a
derivative of pure roundoff (measured 1e-32 rms), and
`den > 0` let it print a ratio of 96 where NaN was meant.

The self-diagnostic comment's first line:

THE SELF-DIAGNOSTIC (Part I's "only if", built 2026-09-08): the

### `_spectral_report`

The docstring paragraphs before the move:

RECORDED SCOPE ITEM 3.  A k-step method turns an m-dimensional
system into a k*m-dimensional discrete one, so the composed
monodromy's spectrum carries `(k-1) m` PARASITIC roots beside the
physical Floquet multipliers.  `max |eig|` over that mixture is only
a stability verdict while the parasitic roots stay small -- which
for Gear-2 they emphatically do (`(1/3)^N`, ~1e-95 at 200 points)
and for a method whose spurious root sits nearer the unit circle
they would not.  This separates them instead of hoping.

Continuing after the discriminator's two cases:

So `||v_{-1} - v_0||` (against a unit-norm eigenvector) is O(h) for
a physical mode and O(1) for a parasitic one.  MEASURED, and it is
the h-scaling that makes it a prediction rather than a story: on the
phase circuit the physical ratio falls 0.1281 -> 0.0316 when the
grid goes from 50 to 200 points -- a factor of 4.05 for a factor of
4 in `h` -- while the parasitic ratios sit at 1.0 to 10.  On the
Q=20 RLC the parasitic ratio is 1.9997 against the 2.0 that BDF-2's
`v_{-1} = 3 v_0` predicts exactly.

⚠ THE MODE COUNT HERE IS AN ODE COUNT AND THE OBJECT IS A DAE, and
the difference is structural rather than an off-by-`k`.  Demir
(IJCTA 28:163-185, 2000) gives the DAE monodromy as

    Phi(t,s) = U(t) D(t-s) V(s) C(s)

with `D = diag[exp(mu_1 (t-s)), ..., exp(mu_d (t-s)), 0, ..., 0]`
for `d = rank(C)`: "equation (19) has k = n - m Floquet multipliers
that are 0", and on a real circuit "there are also eigenvalues
exactly equal to 0 due to the ALGEBRAIC EQUATIONS in the MNA
formulation".  So the `m - rank(C)` structural zeros are the
theory's, not an artefact -- which is why `parasitic_roots` comes
back identically zero on every MNA circuit tried here.  ⚠ AT INDEX
1 ONLY (the docs session, checked against the paper 2026-09-09:
"We assume that the DAEs we are dealing with are index-1").
`rank(C)` is the differential dimension at index 1 and OVERCOUNTS
by one per index-2 constraint: measured on `floquet_modes`, an
index-1 tank and an index-1 tank + R node give modes = rank(C) = 2,
an L-I cutset and a C-V loop give rank(C) = 3 with 2 modes -- the
code returns the true count; it is the formula that stops where
Demir says it does.

⚠ AND THE FACTORISATION CARRIES A TRAILING `C(s)` WITH NO ODE
ANALOGUE (where `C = I` and it disappears).  A DAE monodromy is not
simply a product of state-transition blocks, so an ODE-shaped
count does not merely miscount -- it describes a different object.
Anyone revisiting this split should start there and not from the
eigenvector heuristic below.  Relayed from the docs session's read;
check it against the paper before building on it.

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

⚠ THE COUNT IS AN ODE COUNT, AND MNA CIRCUITS ARE DAEs.  This
splits `2m` eigenvalues as `m` physical and `m` parasitic, which is
right for an ODE.  An index-1 MNA system with `d = rank(C) < m` has
`d` physical multipliers, `d` parasitic ones and `2(m - d)`
STRUCTURAL ZEROS from the algebraic variables -- so on a real
circuit both arrays are mislabelled: measured on the Q=20 resonator
(`m = 4`, `rank(C) = 2`), `parasitic_roots` comes back identically
zero and `floquet_multipliers` carries two structural zeros beside
the two real multipliers.

⚠ `spectral_radius` IS UNAFFECTED, which is why this is recorded
rather than re-engineered.  The physical multipliers have the
SMALLEST block split by construction, so they are always inside the
first `m`, and the maximum over that set is the right number --
0.97531 on that circuit, against the analytic 0.9753.  What is
unreliable is the LABELLING of the diagnostic arrays.  And it cannot
be fixed by magnitude either: Gear-2's true parasitic roots are
`(1/3)^N`, about 1e-95, which is numerically indistinguishable from
a structural zero -- so on this method the two populations cannot be
told apart at all, by any test, and saying so is the honest
position.

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


## `_pss_periodic.py` -- `_PeriodicStates`

### `_fold_periodic`

The docstring as it stood before the history moved out (2026-09-24):

Fold a shooting residual's periodic rows into `[-m/2, m/2)`.

⚠ THIS IS WHAT MAKES THE PERIOD MAP'S FIXED-POINT PROBLEM WELL POSED
FOR A FOLDING STATE, and it belongs in the RESIDUAL rather than in the
traversal.  `Idtmod`'s state is defined only up to `n*modulus` -- it
says so itself through `periodic_states()`, and the transient engine
already uses that declaration to keep the state bounded by exact gauge
translations.  Shooting did not: it asked for `x_0 - phi(x_0) == 0`
literally, which on a folding row demands the SAME REPRESENTATIVE, not
the same state.  An orbit that closes after advancing exactly one
modulus -- the normal case for a phase -- then has NO root at all, and
near the fold the raw difference jumps by a whole modulus while the
state moves infinitesimally.

Measured (see `test_a_state_fold_breaks_the_period_map_at_the_ENDPOINT_
not_on_the_grid`): that jump is grid-INDEPENDENT -- 1.414214e+09 at
seven different grids, with the wrap on a node and off it alike -- so
it is not an event-localisation defect and no refinement of the time
grid can reach it.  It is the output map, and the output map is what
this folds.

⚠ THE JACOBIAN IS DELIBERATELY NOT TOUCHED.  `d/dx0` of
`wrap(x0 - phi(x0))` equals `d/dx0 (x0 - phi(x0))` almost everywhere --
the wrap has unit slope between its jumps -- so `D - alpha*Mx` is
already the right derivative of the folded residual.  The fold moves
the residual onto the branch the Jacobian was always describing; that
is the whole reason this is a local change and not surgery on the six
traversal loops.

The comment on the block loop:

every STATE of the unknown folds on its own rows -- gear's pair
`(x_0, x_{-1})` carries two (its autonomous and event residuals did
not fold at all before 2026-09-23; the driven one folded each half)

### `_wrap_jump`

⚠ THE STATE FOLD ALONE LEAVES THE OUTPUTS DISCONTINUOUS.  An
idtmod's OUTPUT is `wrap(state)`, and every algebraic quantity it
feeds (the phase node, the branch current into a load, a phase
detector's node) jumps with it.  When the orbit wraps exactly at
``t = 0`` -- a free-running VCO pinned at a zero crossing of
`sin(2 pi phase)`, or a reference phase that starts at 0 -- `z_0`
sits just after the wrap and `z_end` just before it, the state row
folds to zero, and those rows still differ by the whole jump.  At
rounding level the Newton iterates straddle the wrap and the jump
flips in and out of the residual: radau on the free-running
`VcoHdl` failed there (2026-09-23), with every other row at 1e-14.

⚠ WHICH SIDE A STATE IS ON IS DECIDED BY THE STATE, NEVER BY THE
RESIDUAL.  The declared window's edges are where the output wraps,
so the fractional window positions `f_0`, `f_end` of the two states
say whether the state fold's nearest-representative path crosses an
edge: ``k = -round(f_0 - f_end)``, which is -1, 0 or +1.  Folding an
output row by the modulus instead (round the residual) accepts a
start point one modulus off its own constraint -- measured: radau,
gear and trap all "converged" with `ph(0) = 1.1` against a phase of
0.1 -- and cannot reach a loaded output at all (1 mA on a load
resistor's current, 2 V on a gain-2 detector's node).


## `pac.py` -- `PAC`

### `(class docstring)`

The class docstring before the move:

THE OPERATOR IS THE MONODROMY, and the whole method is one line of
algebra on the withdrawn implementation's own system.  That system was

    (L + alpha B) v = -u,   alpha = exp(-2j pi f T)

with `L` the block lower bidiagonal discretisation over the period and
`B` the periodic wrap.  Telichevesky, Kundert & White (DAC 1996) reach
the iterative form by "reinterpreting the use of `L^-1` ... as a
preconditioner":

    (I + alpha L^-1 B) v = -L^-1 u

`L` is block lower bidiagonal, so applying `L^-1` is forward
substitution through the timesteps -- which is the recursion PSS already
runs against stored factors -- and `B` is confined to the first `m` rows
and last `m` columns, so `L^-1 B` acts only on the LAST block.  Both
claims are checked against our own matrices in
`test_the_pac_operator_is_the_monodromy_and_L_is_never_formed`.

What is left after that is `m x m`:

    (I - alpha M) y_0 = alpha w(f)

with `M` the monodromy and `w` the forced response over one period from
a zero initial state.  `y_0` is the small-signal state at `t = 0`; one
more driven replay gives the rest of the period.

⚠ WHY THE 419.5 GiB IS GONE, precisely.  It was never the operator: it
was the cost of FORMING `L` and `B`, `(N m)^2` complex entries, 279.7 +
139.8 GiB at `N = 137`, `m = 1000`.  Nothing here forms either.  The
stored per-step factors PSS already makes are the preconditioner, and
the only dense object is `m x m` and only if the caller asks for it.

⚠ AND THE OLD `L` WAS BACKWARD-EULER-SHAPED, which is the trap a rewrite
falls into.  It has two terms per row; a two-step method's variational
system has three.  Rebuilding it for `trap` or `gear` gives an operator
for a different recursion than the trajectory it came from -- measured,
spectral radius 0 against the analytic 0.8546
(`test_the_pac_L_is_backward_euler_only`).  Taking `M` from the
traversal cannot make that mistake, because every step carries its own
`(alphas, b)`.

### `solve`

The comment on the source vector:

⚠ `analysis=` BY KEYWORD.  `Circuit.u(t, epar, analysis, ...)`
takes `epar` second, and the withdrawn body wrote
`self.cir.u(0, analysis_name)` -- passing 'ac' as the element
parameter set and taking the TRANSIENT source vector, which is
zero at `t = 0` for every sinusoid.  The whole analysis would
have returned zeros, silently, with no source to speak of.

The comment on the manufacturing step:

⚠ THE MANUFACTURING STEP IS NOT IN `steps`, AND IT COSTS AN
ORDER.  On the plain path `_traverse_factored_plain` takes one
step OUTSIDE the loop to manufacture a history, and folds its
effect into the `opening` triple as a flat-history assumption.
For the HOMOGENEOUS map that is the documented approximation the
whole plain path is built on.  For the DRIVEN one it also means
the source is never applied at that step -- one step of `u` out
of `N`, i.e. a relative O(h).

⚠ MEASURED, on the Q=20 resonator against the AC analysis at
700 Hz, rel error per doubling of the grid:

    trap, plain            2.00x  (O(h))   4.13e-03 at 250 pts
    trap, x0_unknown=True  4.00x  (O(h^2)) 1.09e-04 at 250 pts
    euler, either          2.00x  (O(h))   1.40e-02, unchanged

The euler row is the control: `x0_unknown` does not move it at
all (identical to five digits), so the trapezoidal gain is the
manufacturing step and not something else the formulation does.
The trajectory is NOT the problem -- trap's waveform converges at
4.2x per doubling either way.

So this is a silent order loss for a caller who did nothing
wrong, which is the one thing worth a warning.  Gear-2 takes the
solved-history path and has no manufacturing step at all.

The comment on the deflated route:

⚠ ON AN OSCILLATOR THE OPERATOR HAS THE ANSWER'S OWN POLE at every
harmonic (see `_check_harmonic`), and a plain solve near one
carries relative error `eta / (2 pi df/f0)`, `eta = |lambda_1 - 1|`
the computed unit multiplier's displacement -- measured to four
digits over five decades (docs session, Gourary reading).  The
deflated route (`_deflated_solve`) borders the pole out and is
exact there; it was wired into `adjoint_sideband_row` only, and
this sweep solved plain outside HARMONIC_GUARD (2026-09-08).
Under the radau default eta ~ 1e-12 puts the unguarded band
inside the guard, so this is correctness hygiene, not a fix a
user would see; the subspace recycling across frequencies is
given up on the autonomous path (one bordered solve per point).

The comment on the staged oscillator:

⚠ ON A STAGED OSCILLATOR (2026-09-22, events phase B) the
bordered system collapses onto the total map: with `dtheta
= dtheta/dx_0 y_0 + dtheta_f`, `dtheta_f = -Gt^-1 W f_node`
the source's own motion of the crossings, `(I - a M_tot)
y_0 = a (w + P_theta dtheta_f)` -- the deflated solve with
the total operator and this source.  Exact against the
piecewise-linear forced response (see the test); the plain
deflated solve on such a solve was 0.3-400x off.

The comment on the bordered sideband response:

⚠ THE BORDERED SIDEBAND RESPONSE (2026-09-22, events phase B).
On a solve whose grid was landed on state events, a periodic
perturbation moves the crossings: `y_end = M y_0 + w + P_theta
dtheta`, and the event rows close it -- `w_k . y(node_k) = 0`
with `y(node) = P_node y_0 + f_node + Pk_node dtheta` (the
homogeneous map to the node, the forced response there, the
event column there).  Solved by block elimination: the m x m
solve for the source and for each event column, then the K x K
Schur complement for `dtheta`.  A per-step saltation instead of
this read the dominant multiplier 28 % short of the exact total
(see `_state_event_stage`); the bordered system IS the
linearisation of the solve that produced the orbit.

The comment on the fixed-time event columns:

the crossings' motion at every node, AT FIXED TIME --
`Pk_j - xdot_j tau_j^T` (2026-09-22): the response of "node
j" itself includes the node's motion along the orbit,
O(1) of the response on a staged oscillator; with the
motion removed the exact forced response is matched to
1e-3 at every node (see `_fixed_time_event_columns`)

The comment on the DFT of the response:

`v(t) = y(t) exp(-j w t)` is T-periodic; its DFT is the
sideband set, exactly as the withdrawn body intended

The comment on the two reporting defects:

⚠ TWO REPORTING DEFECTS, FOUND BY AN EXTERNAL REFERENCE CROSS-CHECK
(2026-09-05), neither in the solve.  (a) `fp.times` spans
`[0, T]` INCLUSIVE, so the last sample repeats the first on a
T-periodic `v` (|v[0] - v[-1]| / |v[0]| = 7e-18 measured) and
the DFT's `dt = T/(N-1)` put the sidebands at `f0 (N-1)/N`:
99 500 Hz for 100 000 at N = 200 -- and cost an ORDER, O(h)
for O(h^2), 68x at 800 points.  `PSS.solve` already drops
the endpoint one function away; this did not.  Guarded on
the window rather than sliced blind, since the plain path's
`[:len(y)]` need not be inclusive.  (b) `|sb + f|` folded a
NEGATIVE sideband frequency to positive and left the
coefficient alone; the physical response there is the
CONJUGATE.  Uncorrected, `l = -1` was 166% off and did not
converge under refinement; conjugated it lands on its
positive twin's error to three digits (4.873e-3 / 4.877e-3).
Both defects are invisible on a circuit whose `v(t)` is
constant over the period -- every earlier PAC gate.

### `adjoint_transfer_row`

The end of the docstring before the move:

replay, and the whole row falls out. MEASURED against `m` forward
solves on an RC ladder: agreement 9.6e-16.

⚠ WHAT THIS IS NOT, so the next reader does not over-read it. The
output here is the state at `t = 0`, a single linear functional.
A SIDEBAND coefficient `H_l` is a functional DISTRIBUTED over the
period -- `(1/N) sum_n exp(-j l w0 t_n) d^T y_n` -- and its adjoint
needs the reverse pass to take an injection at every step rather
than a seed at the end. That extension is the next piece of A3, and
it is not built.

⚠ WAS SOLVED-HISTORY ONLY until B8 gave the one-step companions
their own reverse recursion; it now runs under every method.

The comment at the top of the body:

⚠ NO LONGER GEAR-ONLY (B8): the transposed replay exists for the
one-step companions too, and every use below goes through
`fp.matvec_transposed`.

### `adjoint_sideband_row`

The end of the docstring before the move:

Dropping the second term would leave an answer that looks entirely
reasonable: MEASURED on an RC ladder the two terms are comparable
in size (303 against 498 at `l = 0`, 79 against 606 at `l = 1`), so
neither is a correction to the other.

Still ONE transposed solve per sideband whatever the number of
sources, which is the property pnoise needs.  Agreement with the
`m` forward driven solves: 9.2e-16 / 3.3e-16 / 1.3e-15 at
`l = 0 / 1 / -2`.

⚠ WAS SOLVED-HISTORY ONLY, like the reverse pass; B8 lifted both.

The comment at the top of the body:

⚠ NO LONGER GEAR-ONLY (B8) -- see the adjoint row.

The comment on the bordered row:

⚠ ON A STAGED SOLVE THE ROW IS BORDERED (events phase B,
2026-09-22): the transpose of `PAC.solve`'s bordered system,
the output read at FIXED times (`g_theta` over the fixed-time
columns), and the event rows' term -- `-zeta_k W_k` at node
k, the source coupling of `f_node_k` -- as a second reverse
pass.  On a driven solve `z, zeta` come from the block
elimination (`EventColumns.bordered_adjoint`); on an
oscillator the system collapses onto the TOTAL operator,
deflated, with `zeta` read off after it.  Radau/trbdf2 and
gear's pair map alike (an autonomous gear solve stays
unbordered).  Verified by dual consistency against the
bordered forward solve (the suite's own PAC test pattern);
unbordered, pnoise on a staged gear solve was 10-15 % off.

### `pnoise`

The comment above `pnoise` (its last sentence was cut off in the source):

⚠ THE FOLD BELOW IS FOR DRIVEN CIRCUITS.  It is a frequency-conversion
computation and is complete for one; for an AUTONOMOUS oscillator it is
structurally incomplete -- the near-carrier phase-noise skirt is not a
conversion effect (Rizzoli, Mastri & Masotti, MTT 42-807, 1994).  Free-
running phase noise goes through the Floquet/PPV stack instead; see
`oscillator_spectrum` for why the two cannot be unified and why the wrong

The docstring before the move (after its first line):

⚠ `cyclostationary=True` IS THE CONSTRUCTION FOR A BIAS-DEPENDENT
`CY` (2026-09-08, on the corrected Okumura reading).  A source whose
PSD follows the orbit is white noise `xi` MODULATED by
`B(t) = sqrt(CY(x(t)))`, a T-periodic matrix with Fourier
coefficients `B_k` -- read off the PSS samples by one DFT, no
window count `p` at all (Okumura's windows are a piecewise-constant
approximation of exactly this, and their boxcar coefficients its
crude version).  `xi`'s band at `g_p = f - p f0` reaches the output
at `f` through EVERY modulation harmonic `k` and the sideband row
`a_{p-k}` (source at `g_p + k f0`, output at `f`) -- the SAME rows
the stationary fold computes -- COHERENTLY over `k` (one white band,
one realisation) and incoherently over `p`.  Summing the bands
turns the square root into the PSD's OWN harmonics `P_j` (the DFT
of `CY(x(t))`, no matrix square root anywhere):

    S(f) = sum_{l,l'} a_l P_{l'-l} a_{l'}^H,

which is exact on the grid (⚠ the sqrt-modulation form, tried
first, left a 2.8e-5 residual tied to the modulation's zero
crossings; this form agrees with the stationary side to 9e-16).
Constant `CY` gives `P_0 = CY` and nothing else, and the sum
collapses to the stationary `sum_l a_l CY a_l^H` -- Okumura's
`p = 1` case, pinned to machine precision.  The cost is the
stationary fold's (the rows dominate; the double sum is free).
A coloured source is folded band by band (each white band the
modulation reaches carries its own `CY`; see `_cyclostationary_fold`).
Like the stationary fold this is a LOWER bound at a sideband cap.
⚠ Coherence is the whole content: `modulated=True` (the cycle-
averaged `CY`, Hull & Meyer's stationary equivalent) keeps the
power and drops the correlation between sidebands, and the two
differ wherever the modulation has harmonics -- measured on a
driven multiplier, and the identity against the STATIONARY fold of
the same physics written as a white source through a periodically
varying gain is the gate (`test_..._cyclostationary_...`).  ON A
MOS (2026-09-09): an EKV stage switched by a 1 MHz LO whose channel
noise passes through a second EKV switched by the same LO reads
0.376 of the cycle average (thermal, every offset) and 0.32-0.44
with flicker at ten times thermal; the switch's own channel noise
is largest when its channel shunts it.  ⚠ Noise that reaches the
output through a time-INVARIANT transfer (a single stage's drain
into an RC load) gives cyc = cycle average to four digits, since
only P_0 survives -- the construction shows only where the
modulated noise crosses a periodically varying transfer.  Cost:
white = the cycle average's; coloured (any frequency-dependent
CY, a negligible flicker coefficient included) ~6x.
⚠ FLICKER, AND WHAT OKUMURA'S EQ. 23 MEANS HERE (measured
2026-09-08): a coloured source is folded band by band, and against
the stationary fold of the same SEPARABLE physics (a stationary
flicker source through a periodically varying gain) it is exact --
1.000000 -- as long as the modulation is SIGN-DEFINITE.  When the
modulation changes sign the two are DIFFERENT physics (0.56 / 1.33,
grid-independent to six digits): for white noise `m xi` and `|m| xi`
are one process, for a coloured one whose correlation spans the
sign change they are not (`R(t,t') = m(t) m(t') R_c(t-t')` keeps
the sign product), and a PSD cannot carry the sign -- so this fold,
like the HDL model feeding it, is the `|m|` one.  That is Okumura's
"cannot be modeled as a cyclostationary process by using this
method, because it has very long time constants" in concrete form.
A flicker source with a bias-dependent coefficient gets the `|m|`
number, correct when its modulation does not change sign.

Returns `(S, sidebands_used)`.  `S` is the one-sided
**time-averaged** PSD at the output, in the same units as
`analysis_ss.Noise`'s `Svnout`.

⚠ "TIME-AVERAGED" IS NOT A HEDGE, IT IS THE SPECIFICATION, and
saying so is the whole of this paragraph's job.  TWO SEPARATE
MECHANISMS make an output noise cyclostationary, and only one of
them is about the sources:

  1. bias-dependent sources modulated by the time-varying operating
     point -- this is what `_cy_reduced` refuses, because the
     stationary sum would be the wrong model;
  2. the PERIODIC SOURCE-TO-OUTPUT TRANSFER FUNCTION -- which
     applies even when every source is stationary.  A circuit whose
     only noise is the thermal noise of constant resistors STILL
     has cyclostationary output noise.

The sideband sum here handles (2) correctly and returns its TIME
AVERAGE.  That is the right answer for most uses and it is
incomplete for two ordinary RF topologies, both named by Kundert
(*Introduction to RF Simulation*, v2 2003 -- relayed from the docs
session, cited not verified here): a NONLINEAR SUBSEQUENT STAGE
("an oscillator drives a limiter ... the same is true when an
oscillator drives a mixer"), and CASCADED STAGES OFF A SHARED
REFERENCE, where "the second mixer is synchronous with, and tracks
the variations in, the cyclostationary noise of the first."  The
test is whether anything downstream can track the PSD's variation:
if it cannot, the phase is unknown to it and the time average is
sufficient.

⚠ AND A SCALAR CANNOT CARRY WHAT IS MISSING.  Cyclostationary noise
is CORRELATED between frequencies separated by `k f0`, where
stationary noise has no correlation between different frequencies
at all.  This returns one number per output frequency, so it does
not represent that correlation -- deliberately, and stated here
rather than left for a caller to discover by getting a wrong answer
in one of the two topologies above.

    S(f) = sum_l  h_l CY h_l^H ,   h_l = H_l(f - l f0)

Noise entering at `f - l f0` leaves at `f` through sideband `l`, and
white sources in disjoint bands are uncorrelated, so the bands add
in POWER.  Each `h_l` is one adjoint row -- one transposed solve for
every source in the circuit, which is the whole reason this is
affordable.

⚠ A PRECONDITION FOR THE FIRST COLOURED SOURCE, recorded here
because it is unreachable today and will be silent when it is not.
A 1/f source is singular at DC, and folding puts a copy of that
singularity at EVERY harmonic.  A commercial RF simulator: "place a cluster of
frequencies near each harmonic ... but AVOID PUTTING FREQUENCY
POINTS PRECISELY ON THE HARMONICS ... you run the risk of
generating absurd noise totals because a very narrow noise peak
artificially has its apparent width greatly magnified by a large
frequency, and has its amplitude exaggerated by placing a point
precisely at the singularity."  Plausible nonsense, no error
raised.  Every source in the discrete library is white, so `freq`
landing on `k f0` is harmless now; it stops being harmless the day
one is not.

⚠ AND AN OSCILLATOR IS NOT THIS FUNCTION'S PROBLEM AT ALL.  A
driven circuit's output noise is cyclostationary; an AUTONOMOUS
one's is STATIONARY, and structurally so -- "cyclostationarity in
the oscillator's output would, by definition, imply a time
reference ... noisy autonomous systems cannot provide a perfect
time reference" (Demir 2002).  That is the physical counterpart of
`I - M kron M` being exactly singular for an oscillator
(`test_no_periodic_covariance_exists_for_an_oscillator`): there is
no cyclostationary object to compute, not a hard one.  Oscillator
phase noise is a different output shape entirely -- a closed form
in a few scalars with no frequency sweep -- and is not built.

⚠ WHY CYCLOSTATIONARY IS NOT BUILT -- AND THE REASON RECORDED HERE
FIRST WAS WRONG.  This said the cross terms need "the `R_{m,n}`
construction from section III-B", unread, as though the window
Fourier coefficients were an exotic object.  They are not.  The
model is `c(t) = sum_m n_m(t) w_m(t)` with `w_m` a T-periodic
RECTANGULAR window over interval `m`, non-overlapping -- so
`W_{m,k}` is the Fourier series of a BOXCAR, closed form, a `sinc`
times a phase.  The `n_m` are taken UNCORRELATED across intervals,
justified because `H(jw,t)` is time-invariant within each one, so
the sum is INCOHERENT over `m` and coherent only over `k` within a
single interval.  Nothing there is missing.
⚠ THE ACTUAL BARRIER IS COST, which is a different decision -- and
the cost as first recorded here was OVER-STATED (verified at the
source by the docs session, 2026-09-08, Okumura et al. 1993).  The
source count is p x (noisy devices) where p is the number of
intervals over which "H(jw,t) is time-invariant within each
interval" -- set by how fast the transfer varies (their Fig. 2 has
FIVE windows; a switching circuit moves fast only at transitions),
NOT by the integration grid: the earlier "500-point grid x 50
devices = 25 000 sources" tied p to the timestep and was high by
(timepoints)/p, an order or two.  The reported noise analysis ran at
14.1x the PSS per frequency point (1086 s vs 77 s), "because all
aliasing components need to be computed" -- and the NEXT sentence,
elided before: "it is expected that this problem can be greatly
alleviated using a vectorization technique, because most of the
computational power is used to solve linear problems" -- which is
the batched JAX path this tree already carries.  Whether that is
affordable is unmeasured here.
⚠ AND THE METHOD CANNOT MODEL FLICKER (p. 585, verbatim): "Flicker
noise generated under a periodic large signal excitation cannot be
modeled as a cyclostationary process by using this method, because
it has very long time constants and thus equation (23) does not
hold" -- eq. 23 being the uncorrelated-across-intervals assumption.
Their fallback is that flicker "may exist as independent noise
sources which are practically modeled as stationary random
processes".  So the construction covers cyclostationary thermal
and shot noise, and NOT one of the three mechanisms `_cy_reduced`
names it as the precondition for.
⚠ AND ITS AUTHORS LEFT THE PHYSICS OPEN: "it is further necessary
to discuss the correspondence between the actual physical phenomena
of noises and this modeling".  The windowed-stationary
decomposition is a numerical construct, and its fidelity to a real
device is not settled by its numerical validation.
(Verified at the source 2026-09-08.  The boxcar is the paper's own
closed form, R_{m,n} = (h_m/T) Sa(n w_s h_m/2) exp(-j n w_s (tau_{m-1}
+ h_m/2)), p. 585 -- an earlier line here claimed the observation as
ours.)

⚠ STATIONARY SOURCES ONLY, AND IT CHECKS -- mechanism (1) above.
Okumura's cyclostationary model windows each source to a single
timestep, and the windows'
Fourier coefficients then CORRELATE the sidebands -- they stop
adding in power, and the cross terms need the `R_{m,n}` construction
from his §III-B.  Every noise source in this element library is
bias-INdependent (a resistor's `4kT/R` does not read `x` at all), so
the stationary formula is exact for them; a compact device with a
bias-dependent `CY` is not covered, and this raises rather than
returning a number that is quietly the wrong model.

⚠ `maxsidebands` IS AN ACCURACY KNOB HERE AND A REPORTING KNOB IN
`PAC.solve`, WHICH IS THE OPPOSITE OF HOW IT READS.  A commercial RF simulator's
own documentation states the inversion (relayed, cited not verified
here): reducing sidebands "affects only the amount of information
generated, not its quality.  HOWEVER, NOISE SOURCES GENERATE
SIGNALS AT ALL FREQUENCIES, and therefore with PNoise, reducing the
number of sidebands acts to REDUCE THE NUMBER OF NOISE
CONTRIBUTIONS in the output and so REDUCES THE ACCURACY of the
result."  A driven signal lives at the frequencies it is driven at,
so dropping sidebands drops answers you did not ask for; noise
lives at all of them, so dropping sidebands drops power that
belonged in the total.  Capping it here always makes `S` a LOWER
bound, never a cheaper estimate of the same number.

⚠ AN OBSERVABLE SYMPTOM WORTH KNOWING BEFORE IT IS SEEN.  For an
oscillator `Phi(T) - I` is singular and its null vector IS THE PPV,
so a near-carrier noise computation is ill-conditioned by
construction.  Gourary et al. name what that looks like: "the
standard time domain noise analysis yields FLAT PSD CURVES OR
CURVES WITH UNEXPECTED SLOPE NEAR THE OSCILLATION FREQUENCY."  If
oscillator noise ever comes out flat near the carrier, that is the
singularity -- not the physics, the noise models or the source
definitions -- which points at the right layer immediately.  The
published removal (Gourary et al., eq. 27/28: replace the output
row of J^T by u^T; verified at the source by the docs session,
2026-09-08) IS built here as `_deflated_solve`, which borders with
BOTH null vectors and is the better conditioned of the two; it is
wired into `adjoint_sideband_row` (so into this method) and, since
2026-09-08, into `PAC.solve` and `adjoint_transfer_row` as well
(`PAC.deflated` says which route ran).  The plain solve's relative
error near a harmonic is `eta / (2 pi df/f0)` with `eta =
|lambda_1 - 1|` the computed unit multiplier's displacement --
measured 1.1e-12 (Q = 16) and 1.8e-13 (Q = 100) under radau, so
under the default integrator the unguarded band sat inside
`HARMONIC_GUARD` and the wiring is hygiene; under gear at
df/f0 = 1e-10 the plain solve refuses outright (GMRES residual
1.7e-6) where the deflated one carries the pole to 1 %.

⚠ TWO STOPPING RULES, AND THE BOUND IS NOT THE RATIO TEST.  The
accumulation stops when a sideband pair adds less than `ratio_tol`
of the running total -- and it can never pass `|l| <= N/2`, the
grid's own Nyquist, because nothing aliases down from above the
maximum frequency the grid represents (eq. 32).  An implementation
with only the ratio test terminates for the wrong reason and, on a
coarse grid, after summing harmonics its own grid cannot carry.

GATED against `analysis_ss.Noise` on a linear circuit, where the
sidebands vanish and this must reduce to the stationary answer --
Okumura's own `p = 1` case, "exactly the same as that derived for a
stationary noise".  Measured ratio 1.000000, with every `l != 0`
contributing ~1e-32 of the total.

⚠ AND THE TIME-AVERAGE CHOICE MATCHES THE REFERENCE IMPLEMENTATION,
which is worth recording because it was documented above as a
deliberate scope decision and could have been the wrong one.
a commercial RF simulator's theory notes on PNoise and QPnoise, both: "THE TIME-AVERAGE of
the noise at the output of the circuit is computed in the form of a
spectral density versus frequency."  Same quantity, same
limitation.  (Relayed from the docs session; cited, not verified
here.)

The comment on the colour model:

the stop rule and the harmonic probes below run on the
cycle-averaged power (the modulation's B_0 B_0^H); the fold
itself is the convolution after the rows are gathered.
The colour model (see `_cy_colour_model`) is fitted ONCE
here and serves both: the 34 orbit sweeps of the stop rule
were 1.1 s of a 4.1 s call after the fold itself was cut.
per-element components, so independent sources ADD in the
coloured fold (see `_cy_components_model`); its call is the
summed model the stop rule reads, as before

The comment on the fold to DC:

⚠⚠ ON A HARMONIC, A SIDEBAND FOLDS THE SOURCES TO DC -- AND
SOME DEVICE MODELS ARE NOT DEFINED THERE.  Sideband `l`
evaluates `CY` at `f - l f0`, so `f = k f0` evaluates it at
ZERO.  A `1/f` term is infinite there; and MEASURED on
`PspMosLongChannel`, a flicker term with its coefficient set to
ZERO is `0/0` and returns `nan`:

    fnt=1, nfa=0        CY(f=0) = nan   <- DISABLED flicker
    fnt=1, nfa=8e22     CY(f=0) = inf   <- the real singularity

⚠ THE FIRST IS THE NASTIER ONE: a caller who sets `nfa = 0`
believing flicker is off still gets `nan` out of `pnoise`, with
no exception anywhere.

⚠ AND IT IS NOT "HARMONICS ARE BAD".  A driven divider with
white sources returns 1.490351e-17 at exactly `f0`, and at
`2 f0` and `3 f0` -- the fold to DC is harmless when the
sources are defined there.  So this checks the SOURCES at the
frequency that will actually be used, rather than refusing a
harmonic on principle.

The comment on the finite probe:

⚠⚠ A FINITE PROBE IS NOT A SAFE ONE (peer report, 2026-09-15).
`1/T` rounds (99999.999999999985 Hz for T = 1e-5 s), so at
`f = f0` the folded band sits 1.5e-11 Hz from DC, not ON it: a
1/f source there is finite and enormous, and the cyclostationary
fold returned 6.3e-2 V^2/Hz against 9.2e-15 at 0.1 % either
side, with no warning.  So a frequency-DEPENDENT CY at the folded
band refuses too; a white one stays allowed (its fold to DC is
harmless, the divider value above).

The comment on the steep region:

⚠ AND THE STEEP REGION BESIDE IT IS A SWEEP HAZARD RATHER THAN
A WRONG NUMBER, so it warns instead of raising.  MEASURED with
a real flicker source: the plateau is 1.321483e-16, `f0 + 1` Hz
gives 1.321766e-16 and `f0 + 0.01` Hz gives 1.350047e-16 -- 2%
high, finite, entirely plausible.  The VALUE is right; a grid
that lands there by accident integrates a spike it never
resolved.  A commercial RF simulator: "you run the risk of generating absurd
noise totals because a very narrow noise peak artificially has
its apparent width greatly magnified".

The comment on the stopping rule:

⚠ WHICH RULE STOPPED IT IS PART OF THE ANSWER.  Ending on the
ratio test means the series converged; ending on the Nyquist
bound means the grid ran out before the series did, and the
number is a LOWER bound on the folded noise -- every sideband
above the grid's own maximum frequency is missing, not small.
A strongly switching circuit does this readily: measured on a
driven diode at 80 points per period, the accumulation reached
l = +-39 without the ratio test ever firing, while folding was
already contributing 62% of the total.

### `_cy_harmonics`

The docstring before the move:

for a bias-independent source.  No square root: the fold uses the
PSD's own harmonics (`a P a^H`), which is exact on the grid, where
a sqrt-modulation route (tried first) left a 2.8e-5 residual: the
square root of a PSD that crosses zero has a kink, its harmonic
tail decays slowly, and the convolution's window -- the sidebands
the ratio stop kept, 7 here -- truncated it (measured -1.6e-4 /
-2.8e-5 / -3e-7 at 5 / 7 / 17 sidebands).  The PSD's harmonics
decay fast, so this form is exact at any window.

### `_cy_colour_model`

The comment in `model`:

|w| -- see `_cy_components_model`: a negative band frequency
made this NaN and silently disabled pnoise's ratio stop

### `_warn_signed_unused`

The docstring before the move:

⚠ A FALLBACK THAT REPRODUCES THE OLD ANSWER LOOKS LIKE AGREEMENT
(peer, 2026-09-19: found only by a bit-for-bit A/B against an older
commit).  If an element STATED its signed amplitudes and the fold is
about to factor that component by sqrt(PSD) anyway, say so.

### `_cy_components_model`

The docstring paragraph before the move:

⚠⚠ WHY (2026-09-15, measured).  The coloured fold took ONE symmetric
square root of the SUMMED `CY(x(t), w)` per band.  Independent
sources whose modulations differ do not add under a joint root:
`sqrt(A(t) + B)` cross-couples them at `(t, t')`.  Switch white
noise `4kT g(t)` plus a constant 1/f source at the same node read
pnoise(both) = white + flicker + 7.3 % of the total at 0.013 f0
(+2.2 % at 0.137 f0), and the sampled variance in hold +9.3 % with
white and flicker inside ONE element.  One root per component
restores additivity: 1.7e-16 (separate elements) and 2.3e-11 (one
element, split) against the sum.

The comments in `model`:

⚠ |w|: the stop rule asks at NEGATIVE band frequencies
(f - l f0 < 0), and a non-integer power of a negative base is
NaN -- the ratio stop then never fired and every coloured call
ran to the Nyquist bound with a spurious warning (2026-09-15)
⚠ and numpy division: exactly ON a harmonic the folded band is
w = 0, where a Python float division raised ZeroDivisionError
from inside pnoise's harmonic guard instead of letting it see
the non-finite CY and refuse by name

The comment on the signed amplitudes:

⚠ THE SIGN (2026-09-19): where the element states its coloured
AMPLITUDES, the folds factor the component with them instead of
with sqrt(PSD) -- see `Element.noise_amplitudes` (hdl.py)

### `_uniform_exponent`

The comment before the move:

⚠ ONLY ENTRIES THAT CARRY WEIGHT VOTE (2026-09-19).  The exponent is
fitted from differences of `CY`, so an entry at 1e-12 of the
component's scale -- a MOS flicker source at the sample where Vds
crosses zero -- has its exponent in the rounding of the white part
beside it (measured 1 - 2.2e-09 on an entry of 1.4e-31 against
7.2e-20).  That one entry failed the whole component into the
per-band route at ONE clock amplitude of a sweep.
⚠⚠ AND A WEIGHT CUT-OFF WAS THE WRONG REPAIR (same day, peer, 1000
points): my first fix let entries above 1e-9 of the scale vote, and
a finer orbit landed a sample at 7.1e-09 with its exponent off by
8.1e-09 -- the fallback fired again at two amplitudes and the new
commit reproduced the OLD one to the last bit there.  The exponent's
noise goes as 1/weight, so no cut-off separates them.  What matters
is what a wrong exponent COSTS: giving entry i the exponent `ref`
misstates the component by `r_i |(w1/w)^d_i - 1| ~ r_i d_i |ln(w1/w)|`
of its scale (`r` relative weight, `d` deviation).  Bounded over 50
e-folds of band frequency and held to 1e-9: a genuinely different
exponent (d ~ 1) still fails from a weight of 2e-11 up, while the
two measured offenders cost 2e-19 and 3e-15.  The reference is the
LARGEST entry's exponent, not the first's.

### `_cyclostationary_fold`

The docstring before the move (after its first two sentences):

`CY(x(t))` itself -- no square root, exact on the grid (9e-16
against the stationary fold of the same physics).  COLOURED source
(`CY` depends on `w`; detected by comparing two bands): the white
band `p = l + k` shared by rows `l` and `l'` carries its OWN `CY`,
so `Q_{l,l'} = sum_k B_k^{(l+k)} B_{k+l-l'}^{(l+k) H}` with
`B^{(p)}` the sqrt-DFT at the band's frequency `|f - p f0|`, summed
over ALL `N` modulation harmonics `k` (which is what makes the
square root exact here: the 2.8e-5 of the first version came from
a window on `k`, not from the root).  ⚠ Measured by the docs
session on a flicker source: `||P_0||` differs 24x across the bands
the fold sums, so "the band of l" (the first version's shortcut)
was a 24x approximation on the case the feature exists for; the
band-resolved form is pinned against the stationary fold of a
stationary FLICKER source through the same multiplier.

COST (2026-09-09): the coloured branch was 6x the white one
because of the circuit's `CY` (231 bands x 230 samples), not the
algebra.  With the colour model (`_cy_colour_model`, fitted once
in `pnoise` and shared with the stop rule) and the pair sum
vectorised it is 2.2x the white call and below the cycle average
(1.7 s / 0.8 s / 2.0 s on the switched EKV fixture), exact to
1e-11 against the per-band evaluation, which remains the fallback
for a colour the model does not fit.

The comment on the band-resolved branch:

band-resolved: every white band the modulation harmonics reach.
⚠ THE COST WAS THE CIRCUIT'S CY, NOT THE ALGEBRA (profiled
2026-09-09 on a switched EKV stage: 231 bands x 230 samples =
53 000 CY evaluations, 7.8 s of a 10.4 s fold; the eigen-
decompositions 0.85 s).  Every colour in the library is thermal
plus flicker in 1/f^ef, so THREE evaluations per sample fix each
entry's shape (A + B (w1/w)^ef, ef by a root find on the ratio of
differences), a FOURTH frequency verifies the fit to 1e-8, and
all the bands come from the model with no further circuit
calls; a source whose colour is not of that shape fails the
check and gets the full evaluation as before.

The comment on the specification limit:

⚠ A SPECIFICATION LIMIT, NOT AN IMPLEMENTATION ONE (docs session,
2026-09-08): a coloured source under a modulation that CHANGES
SIGN is not representable by any fold built from a PSD -- the
correlation R(t,t') = m(t) m(t') R_c(t-t') keeps the sign product
and CY cannot carry it -- so this fold, like the HDL model
feeding it, computes the |m| process (measured 0.56 / 1.33 of
the signed one on a flicker source through a zero-crossing
gain, grid-independent; 1.000000000 for a sign-definite gain).
The sign is invisible here; its NECESSARY condition is a PSD
that touches zero along the orbit with a KINK in its square
root, so that is warned on.  ⚠⚠ THE 2026-09-09 SCOPE NOTE THAT
STOOD HERE WAS WRONG AND IS WITHDRAWN (2026-09-19): it said a
DEVICE's own flicker "has no sign to lose -- sqrt(PSD(x(t))) >= 0
IS the process".  A 1/f current is a slow conductance fluctuation
TIMES the current and follows its sign; on a PSP sampler whose Vds
crosses zero while it conducts that is +0.1 % at one clock
amplitude and 400x at another, against a commercial simulator.
Where the element states its signed amplitudes
(`Element.noise_amplitudes`) the folds use them and nothing below
applies; this limit is for sources that state none.  (White
sources are untouched: uncorrelated across the period, no sign
product survives.)  Okumura's eq. 23
objection to flicker is then the separate, physical question of
whether a trap process is "modulated coloured noise" at all.
The proxy's threshold: a zero crossing SAMPLED on an N-point grid
bottoms out near (pi/N)^2 of the maximum (6e-4 at 200 points on
the gate fixture), while a sign-definite PSD with a ten-fold
swing sits at 1e-2 -- so 1e-2 separates them here; a heuristic,
and it is a warning for that reason.

The comment on the order of the zero:

⚠ THE ORDER OF THE ZERO (peer): a LINEAR sign crossing m ~ a t
gives sqrt(PSD) ~ |a t|, a first-derivative KINK; a sign-definite
quadratic touch m ~ b t^2 gives sqrt(PSD) ~ b t^2, SMOOTH.  The
circular second difference of sqrt(PSD) is 2|a|h at a kink and
2b h^2 where smooth -- both shrink under refinement, the smooth
one faster -- so a RAW threshold encodes the grid (5e-3 was safe
at 240 points and a false positive below ~100; peer).  Divided
by h/T and by the maximum it is a DERIVATIVE JUMP, grid-
independent at a kink (2|a|T/s_max ~ 4 pi for a sinusoidal
slope, 12.6 here) and falling as h/T where smooth (~2 (2 pi)^2
h/T: 0.33 at 240 points, 1.0 at 80, 2.0 at 40), so 3 separates
them down to ~50 points per period and the separation grows
with refinement.  Clears the squared-gain case (k V_lo^2, exact
to nine digits) that the touch test alone flagged.  ⚠ STILL
NECESSARY, NOT SUFFICIENT, AND THE DETECTOR'S SENSITIVITY RUNS
INVERSE TO THE EFFECT (peer): the indicator is 12.57 for a
sinusoidal crossing, 0.24 for sign|sin|^1.5 and 0 for sign|sin|^2
-- all sign-changing -- while the discrepancy stays O(1):
MEASURED on the flicker identity with the LO shaped to
v |v|^(p-1), B/A = 0.187 / 1.895 (p = 1, warned), 0.204 / 1.779
(p = 1.5, silent), 0.217 / 1.699 (p = 2, silent) at 0.13 / 1.37
f0; steeper crossings (p = 0.5: 0.161 / 2.046, p = 0.8: 0.175 /
1.970) err MORE and are caught.  So detector and effect are
aligned for p <= 1 and the silent region is exactly p > 1 (the
crossing flatter than linear): the deviation stays O(1) there
while the indicator falls by orders.  A quiet warning is
therefore not evidence of a small discrepancy.

The comment on one root per component:

⚠⚠ ONE SQUARE ROOT PER INDEPENDENT COMPONENT, NOT OF THE SUM
(2026-09-15).  A joint `sqrt(CY)` made independent sources with
different modulations NON-ADDITIVE (+7.3 % of the total on a
switch + a 1/f source at one node) -- see `_cy_components_model`.

The comment on the stacked bands:

every band the sum reaches, stacked once: BB[pi, k] =
B_k^{(p)} with pi = p - pmin.  The (l, l') pair sum is then two
fancy indexings and one einsum instead of N small products in
Python (204 000 `_B` calls, 1.9 s of a 2.8 s fold, before).

### `_band_resolved_pairs`

The comment before the move:

(B B^H)_j = sum_k B_k B_{k-j}^H: the partner index is
k + l - l', NOT k + l' - l -- the mirror was invisible to
the constant-modulation reduction (only k = 0 there) and
read 0.49 / 0.17 on the smooth-modulation flicker identity.
⚠ NO CIRCULAR WRAP HERE: a partner beyond N/2 would be
paired with the wrong BAND (each band carries its own
weight), harmless in the white P-form and wrong here --
it read 0.56 / 1.33 on the kinked (zero-crossing)
modulation whose coefficients reach N/2.

### `_cy_cycle_averaged`

2026-09-26: the same full-width misread as `_phase_mode_split` (`xs[:m]`
of the full-width `waveform`, then the reference's zero a second time): a
source controlled past the reference node read the shifted state.  The
library fixtures place their nodes before the reference, so nothing
pinned moved (`pnoise(modulated=True)` bit-identical on `_nu_identity`).


The docstring before the move (after its first line):

⚠⚠ VALID FOR GENTLE MODULATION ONLY, AND IT FAILS AS A FACTOR, NOT A
PERCENTAGE.  Measured by an external reference-simulator cross-check (2026-09-05)
on a series switch + shunt capacitor, `pnoise` at 10 kHz against
a reference simulator, swept over the modulation depth `goff/gon`:

    goff/gon   1        1e-1     1e-2     1e-3     1e-6
    ratio      1.000    4.33     13.2     15.7     16.0

The 1.000 at the top is what makes the 16 readable: with no
modulation the cycle average IS the value.  The mechanism is not
subtle -- the averaged source injects `4kT <g>` (about half the
on-state current noise) for the WHOLE period, including the hold
phase, where the node it injects into is 1 Gohm in parallel with
100 pF; Hull & Meyer's own condition ("none of the large-signal
state variables may change significantly over the decay time of
the impulse response") fails there by six orders (100 ns closed,
0.1 s open).  So this route is for a mixer's `gm`, a bias-
dependent shot noise -- not for a switch.  ⚠ A switch's noise IS
reachable exactly: `covariance` and `oscillator_covariance`
evaluate `CY` at every step and need no averaging (the switched
capacitor's held variance reads `kT/C` to 1e-4 there).

⚠ THIS IS WHAT `_cy_reduced` REFUSES, DONE INSTEAD OF REFUSED, and
the literature's answer rather than ours.  Hull & Meyer (1993):
*"cyclostationary noise sources, such as shot noise, may be modeled
as MODULATED STATIONARY NOISE SOURCES.  The impulse response that
is calculated INCLUDES THE EFFECT OF THIS MODULATION.  In the case
of shot noise, the hypothetical stationary noise source has
spectral density `S_i = 2q Ibar_c`"* with `Ibar_c` the
cycle-averaged current.

So the modulation is carried by the RESPONSE, which `pnoise`
already computes as `H_l`, rather than by the SOURCES.  Okumura's
route puts one independent stationary source per timestep interval
per device -- `p` per device, ~25,000 sources on a real circuit.
This is ONE per device.  Same physics, `p` times cheaper.

⚠ AND ITS CONDITION IS CHECKABLE RATHER THAN A BLANKET REFUSAL:
*"valid when the impulse response duration is much less than the
time it takes for the mixer circuit to significantly change its
state ... NONE OF THE LARGE-SIGNAL STATE VARIABLES MAY CHANGE
SIGNIFICANTLY OVER THE DECAY TIME OF THE IMPULSE RESPONSE."*

⚠⚠ WHICH IS THE OPPOSITE OF HIGH-Q, AND THEY SAY SO: *"high-Q
filters should be avoided, since they cause the impulse response to
ring, and thus require a very large value of M."*  So this
construction degrades exactly where `lambda_2 -> 1` -- the same
boundary as everything else in this class, arriving from a fourth
direction.  That makes the two constructions COMPLEMENTARY rather
than competing: Hull & Meyer for fast-settling circuits, Okumura's
expensive one for the high-Q case that needs it.  `info` reports
`|lambda_2|` so the caller can see which regime they are in.

⚠ SAMPLED ON THE ORBIT, NOT AT THE OPERATING POINT.  The average
that matters is over the LARGE-SIGNAL waveform, so `CY` is
evaluated at every stored state and averaged with the step weights
-- the same quadrature `diffusion_constant` uses, so the two remain
comparable.

### `_cy_reduced`

The docstring before the move (after its first two paragraphs):

⚠⚠ AND ITS SCOPE IS WIDER THAN "AN UNUSUAL CASE" -- IT IS A
BLANKET REFUSAL OF MOS pnoise.  ⚠ AN EARLIER VERSION OF THIS NOTE
SAID IT WAS UNREACHABLE BECAUSE `PspMosLongChannel.CY` IS
IDENTICALLY ZERO.  That was measured on a DEFAULT-CONSTRUCTED
element: the model has channel thermal and flicker noise, and
`fnt = 0` by default because "an element built without a card is
noiseless".  With `fnt = 1` it is nonzero, white, and
bias-dependent -- so the refusal is REACHABLE from a real device
today, and `modulated=True` is the route past it.  There is no physically correct MOS noise model whose
`CY` is state-independent: thermal channel noise is
`4kT gamma g_d0` with `g_d0` bias-dependent, flicker goes as
`I_D^AF`, gate shot noise as `2qI_G`, and Mahmutoglu & Demir
(2015) are explicit that trap rates "depend on the voltages across
the MOSFET which can considerably vary with time during
large-signal operation".  So the answer to "will a real device
pass this check" is already determined, and it is no.

⚠ THE ROUTE OUT IS THE CYCLOSTATIONARY CONSTRUCTION, NOT A
DIFFERENT DEVICE MODEL, and that reorders the roadmap: the
cyclostationary path is not an enhancement for MOS pnoise, it is
the PRECONDITION -- for the thermal and shot mechanisms.  ⚠ NOT
for flicker: Okumura's own construction excludes it (p. 585,
"cannot be modeled as a cyclostationary process by using this
method, because it has very long time constants"; verified at the
source 2026-09-08), by the same long-time-constant physics as the
trap-rate caveat beside it, and falls back to a stationary flicker
source.  Hull & Meyer (1993) make it affordable -- one
stationary source per device at the cycle-averaged current, with
the modulation carried by the impulse response `H_l` that A1
already computes -- and their worked example IS shot noise
modulated by the collector current, i.e. exactly this case.  Their
condition is checkable rather than a blanket refusal, and it fails
in the familiar direction: a ringing impulse response breaks it,
so it degrades as `lambda_2 -> 1`.

⚠ SECOND-ORDER CONSEQUENCE, WORTH KNOWING BEFORE THE MODEL LANDS.
This same check is what keeps the Ito/Stratonovich choice out of
reach (`CY = GG^T`, so a state-dependent `CY` is a state-dependent
`G`).  Relaxing it for MOS makes the two interpretations diverge,
and Demir's escape -- "the noise signals are small compared with
the deterministic signals" -- may NOT carry for trap noise: a trap
occupancy is a two-state Markov chain rather than a small
perturbation of a large signal, and the same paper says the state
dependence "in fact makes the equation nonlinear".  The tell would
be a discrepancy in a MEAN but not in a variance.

The comment on the probe states:

⚠ THREE STATES ON THE ORBIT -- the first used to be the ZERO
VECTOR, which is on the orbit only by accident, and a linear
time-invariant RC held by a DC clock was refused as
cyclostationary because a switch model read `goff` at v(ck) = 0
(found by an external reference-simulator cross-check, 2026-09-05).  The third
probe is now the stored state half a period in.

### `HARMONIC_GUARD`

The comment (it sits above `_cy_at`) before the move:

How close to a harmonic of `f0` counts as "on" it, as a fraction of
`f0`.

⚠ TIGHTENED BY A7.  This used to be the conditioning floor being
accepted, because `sigma_min` of the plain operator falls LINEARLY
with the distance and everything nearer was unusable.  The deflated
solve removes that: its conditioning is FLAT (measured 2.04e-01 from
0.3 down to 1e-9 of `f0`), so the only remaining reason to refuse is
the physical one -- at an EXACT harmonic `1/(1 - alpha)` is a
division by zero and the response is genuinely unbounded.  So the
guard now excludes only what has no finite answer, not what was
merely hard to compute.

### `_check_harmonic`

The docstring paragraph before the move:

circuit's unit multiplier makes singular.  MEASURED on van der Pol,
`sigma_min(I - exp(-jwT) M)`:

    offset/f0    0     0.25    0.5    0.75    1      2      3
    sigma_min   2.8e-11 0.51   0.65   0.51  2.8e-11 2.8e-11 2.8e-11

and LINEAR in the distance to the nearest one -- 2.5e-1, 2.6e-2,
2.6e-3, 2.6e-4 at 0.9, 0.99, 0.999, 0.9999 of the way there.  The
same sweep on a DRIVEN ladder never drops below 0.78: no unit
multiplier, no singularity, harmonics included.

### `_gmres_checked`

The docstring and the comment at the top of the body before the move:

GMRES, judged by its RESIDUAL rather than by its status flag.

⚠ SCIPY REPORTS BREAKDOWN ON SYSTEMS IT HAS ALREADY SOLVED.  These
operators are `2m x 2m` and often tiny, so the Krylov space is
exhausted in a handful of steps; the next vector is then numerically
zero, which is a LUCKY breakdown -- the solution is exact -- and it
comes back as `info = 4` all the same.  Trusting the flag turns an
exact answer into a `RuntimeError`, which is what it did for AM/PM
at small offsets.

So the residual decides.  A genuine failure still fails, and it
fails with the residual quoted, because the real cause near a
harmonic is that the operator is nearly singular there and no
tolerance will fix it.

⚠ NOW OUR OWN ARNOLDI-GMRES, which returns the residual as its
verdict instead of a status flag that has to be overruled.
The workaround below survives as the TOLERANCE decision; what
has gone is the second opinion about whether the solve failed.

### `_refuse_coloured`

The comments before the move:

⚠ ONE state, two frequencies: the colour question is separable
from the bias question, and asking it through `_cy_reduced`
refused every MODULATED source before the covariance routes
(which evaluate `CY` per step and handle modulation exactly)
could reach it.
⚠⚠ PER ENTRY, AND AT MORE THAN ONE STATE (2026-09-15; reported by a
peer session, reproduced).  This compared the difference against
`1e-9 * max|CY|` over the WHOLE matrix, at `x_last` only.  Once the
PSP gate resistor carried its white 4kT/rg (`d0e7e4e`), a low-rg
device put `1.27e-20` in `CY` and a drain's flicker colour
(`2.5e-30` on its own `1.3e-27`) passed that global threshold --
`covariance` folded 1/f at w0 on a sample-and-hold and said
nothing.  And a source white at the period boundary (a switch OFF
at t = 0) but coloured elsewhere passed at any scale.  So each
entry is judged against ITS OWN magnitude, at `x_last` and at
states spread over the orbit.  Exact zeros and white entries give
identical values at both frequencies, so neither can fire.
`_cy_at` takes a REDUCED state: the reference row is removed here.

### `_lyapunov_pieces`

The end of the docstring before the move:

one implementation, so the pair cannot drift apart in the way that
`diffusion_constant` and `covariance` once did over exactly this
factor of two.

The comment on `CY` per step:

⚠ `CY` PER STEP, AT THE STEP'S OWN STATE -- not one `CY` for the
period.  The Lyapunov accumulation was already per step; hoisting
a single `CY` out of it put a MODULATED source (a switch's
`4kT g(t)`, a MOS channel's `4kT gamma gd0(t)`) outside the
formulation rather than outside the accuracy, and the
cyclostationarity refusal in `_cy_reduced` then closed the door
on exactly the circuits whose noise is the point (found by the
reference-simulator cross-check, 2026-09-05).  Evaluated at the state
the step's companion was factored at (the implicit step's own
solution); the colour refusal still applies -- colour is a
different axis.

The comment in `step_map`:

a one-step companion (gear's Euler backstop past the
zero-stability bound on an event grid, 2026-09-21) has no
third alpha: the pair map holds with it zero -- see
`_monodromy_matvec_transposed`

### `_lyapunov_pieces_plain`

The docstring paragraph before the move:

⚠ `oscillator_covariance` IS REFUSED FOR TRAP-PLAIN, with the
reason: it borders `I - M kron M` with `ppv()`'s null vectors,
which are width `m` on the plain path, and the pair map is
`2m x 2m`.  The pair's own null vectors would be needed, with a
normalisation this record has been burned on twice today
(`floquet_modes`' state-block scale, `ppv()`'s `v . xdot`).
Named rather than approximated; euler-plain and gear both work.

The comment on `CY` per step (a copy of the one in `_lyapunov_pieces`):

⚠ `CY` PER STEP, AT THE STEP'S OWN STATE -- not one `CY` for the
period.  The Lyapunov accumulation was already per step; hoisting
a single `CY` out of it put a MODULATED source (a switch's
`4kT g(t)`, a MOS channel's `4kT gamma gd0(t)`) outside the
formulation rather than outside the accuracy, and the
cyclostationarity refusal in `_cy_reduced` then closed the door
on exactly the circuits whose noise is the point (found by the
reference-simulator cross-check, 2026-09-05).  Evaluated at the state
the step's companion was factored at (the implicit step's own
solution); the colour refusal still applies -- colour is a
different axis.

The comment on the period map:

⚠ THE PERIOD MAP RE-SEEDS THE COMPANION, AND THAT IS
LOAD-BEARING.  The plain product of the A_k carries `iq`
across the boundary, and its `I - M kron M` is SINGULAR:
trapezoidal maps an algebraic row's companion by exactly -1
per step, so the un-reset pair carries a marginal mode --
the `(-1)^n` obstruction this file records for every
formulation that keeps `iq` across a period (measured here:
LinAlgError on a driven RLC).  The shooting solve is
well-posed because the manufactured opener re-seeds `iq` at
zero; the covariance's period map must do the same.  With
`iq` zeroed at the start, the x->x block of the product IS
`fp.matvec` (tied to 1e-12 above), and the map on the pair
is the product applied to `(x, 0)`.

### `_vanloan_step_injection`

The docstring before the move (two lines of its first paragraph, and its last paragraph):

(Roemisch & Winkler; confirmed by the naive two-stage scheme coming
out 27% biased on kT/C).  Van Loan evaluates it exactly: the

Verified against kT/C at second order on R||C (ODE) and VS-R-C
(DAE), matching an independent measurement to the digit; the
stationary error is the METHOD's O(h^2), NOT machine zero (a
machine-zero kT/C would mean a method-consistent `Q = P(1-A^2)`
fudge that corrupts the transient covariance).

### `_lyapunov_pieces_stage`

The docstring paragraph before the move:

Identical in structure to `_lyapunov_pieces_trbdf2`: the per-step
transition `A_n` is the stage step map (dense, `m x m`, via
`_monodromy_matvec_stage` one step at a time) and the per-step
injection `Q_n` is the stage injection (`_stage_injection`) -- the
source enters every STAGE; the end-of-step DAE-projected VAN LOAN
integral (`_vanloan_step_injection`) read a switch's held variance
0.876 kT/C at 400 points (O(h)) and stays the fallback.  State width
`m`, so `n = m`.

### `_lyapunov_pieces_glm`

Built 2026-09-25 (plan item 1, Andreas: "Native for all but covariance").
Until then `_lyapunov_pieces` refused a GLM's map ("a Nordsieck GLM's own
period map acts on its Nordsieck state ... Leave monodromy at a twin"),
and a GLM run's covariance read a radau twin -- which it still does by
default (`_lyapunov_host`).

The design question was the injection.  A white source over a GLM step
reaches the state through the effective stage weights ``w = l^T B`` (`l`
the left unit eigenvector of `V`): GLM2 1/6, 7/12, 1/4; GLM3 0.359,
-0.0167, 0.067, 0.591; GLM4 -26 .. +166.  Independent per-stage samples
with variances ``CY / (2 h w_i)`` (the stage methods' `_stage_injection`)
need every `w_i > 0`; GLM3 and GLM4 fail that.  So the injection is one
shared sample per step, entering every stage, output row and opening
startup substage -- first order, measured: the sampler's held variance
glm2 1.0153 / 1.0055 / 1.0022 kT/C at 200 / 400 / 800 points, glm3
1.096 / 1.045 / 1.022; the jitter sampler's sigma glm2 0.921 / 0.961 of
the analytic at 100 / 200; van der Pol's growth `c_from_growth / c`
glm3 1.00092 / 1.00025 at 60 / 120 points.  The recursion runs on the
map on the state's ``(x, P)`` and closes on `x` (step 0's startup reads
`x` alone); the product of the steps' `x` block is the map to 1e-15.

### `_stage_injection`

The docstring before the move (its first sentence, and the paragraphs after the formula):

the source entering EVERY stage (2026-09-15), or None when the

⚠⚠ WHY (measured on a switched capacitor, Ron 1 k, 100 pF, 100 kHz).
The Van Loan injection freezes `C`, `G`, `CY` at the END of the step,
so across the switch-off edge it integrates the injection with the
OFF conductance: held variance 1 - 0.876 / 0.934 / 0.966 kT/C at
400 / 800 / 1600 points under radau (first order), 0.87 under
trbdf2.  With the stage injection radau reads 1.3e-7 off at 400
points and trbdf2 2.5e-3 (second order, 4.1x per doubling); a
b-weighted Van Loan at the stage states read 5.5e-4 / 4.6e-3.  On a
constant-operating-point RC it is not an exactness fit: the error is
nonzero and falls with the grid (radau 1.4e-9 / 4.4e-11 / 1.4e-12,
trbdf2 2.6e-4 / 6.3e-5 / 1.6e-5 at 100 / 200 / 400 points).
⚠ Needs stiff accuracy (`x_{k+1} = Y_s`) and positive weights:
radau and trbdf2 qualify.  ⚠⚠ A TABLEAU WITH A NON-POSITIVE WEIGHT
(ESDIRK43: `b = 0.158, 0, 0.187, 0.681, -0.275, 0.25`) takes the Van
Loan injection at the STAGE states instead, averaged with POSITIVE
trapezoid weights over the stage abscissae in time.  Measured on the
same sampler, 1 - held at 100 / 200 / 400 points: end-of-step Van Loan
0.367 / 0.224 / 0.124; this 2.1e-2 / 3.4e-3 / 3.7e-4 (tracking and a
constant-operating-point RC identical to the end-of-step form -- the
weights sum to one); equal-variance stage samples, the other
tableau-independent candidate, 5.3e-2 / 3.6e-2 / 2.0e-2 and 6.8e-2 off
on the RC -- rejected.

### `_orbit_rate`

The docstring before the move:

`xdot` at every node of the solved orbit (reduced width) -- THE
DAE'S OWN DERIVATIVE at the node's state (2026-09-22): on the
differential rows ``C(x) xdot = -(i(x) + u(t))``, on the algebraic
rows (a zero row of `C`) the differentiated constraint ``G(x) xdot
= -du/dt``; one small solve per node, exact for the discrete state
and independent of the step.  The rate converts a node's motion in
time into a state change: see `_fixed_time_event_columns`.

⚠ IT WAS A THREE-NODE PARABOLA (one-sided at a landed event), and
that cost 5 % on the comparator oscillator's collapsed `c` node
inside its 10 ns ON phase, where a 7 ns step cannot fit a parabola
to a 10 ns exponential -- the one node of the staged-oscillator
sideband test that had to be excluded.  Against the exact
piecewise-linear rate the DAE form reads 1e-15 at every node
outside the windows (the stencil 10.6 % at worst, 0.45 % at the
node that had to be excluded); on the driven jitter sampler the
sawtooth's rate is its slope exactly and the held capacitor's its
leak.  The stencil is kept only as the fallback where the
assembled matrix is singular (an index above one), and says so.

The comment on node 0:

⚠ NODE 0 IS EVALUATED AS NODE N (2026-09-23).  A source's
derivative at exactly its start (`VSin` clamps `t - td` at 0:
SPICE's rule for a transient) is the LEFT one -- 0 -- where the
periodic steady state, t = 0 == t = T, has the right one: the
rate read 0 for an exact 6283 V/s at node 0 on a driven RC.
Harmless to the consumers (the node's own motion `tau_0` is 0)
and wrong as a function; node N is the same point.

### `_orbit_rate_stencil`

The docstring before the move:

The three-node parabola `_orbit_rate` used until 2026-09-22 --
one-sided AT a landed event and at the node after one, central
elsewhere, periodic at the ends.  Kept as the fallback.

### `_event_closure`

The docstring before the move (its first sentence and its second paragraph):

The BORDERED Lyapunov closure on a staged solve (2026-09-22,
events phase B), or `None` when the solve is not staged.

⚠ THE UNBORDERED CLOSURE ON A STAGED SOLVE IS NOT MERELY
INCOMPLETE, IT IS WRONG BY O(1): the landed window step's OWN
linearisation carries the threshold noise through the switch with
a gain the three collocation points invent -- a comparator-jitter
sampler (a noisy threshold, a hold capacitor on a second ramp)
read 3.8x its analytic held variance `(s_2/s_1)^2 kT/C_n`.  With
the crossing conditions pinned at both window edges the bordered
system cancels that internal sensitivity, which is why the
moving events must be unknowns of the noise problem too.

The comment on gear's pair form:

gear's PAIR form (item 2 of the list after E7, 2026-09-22): the
state is (x_j, x_{j-1}), the event row acts on the first block,
the per-node column of node j is the pair (Pk_j, Pk_{j-1}) and
the map to node j the pair of `P_nodes` rows; the samples come
out as pair covariances, as the plain gear path returns them

### `event_jitter`

The docstring before the move (after its first line):

driven solve (2026-09-22): ``sigma`` in seconds per crossing, and
the crossings' covariance in fractions of the period.

The bordered Lyapunov closure (`_event_closure`) already carries
it: the crossings move as ``dtheta = (dtheta/dx_0) dx_0 - Gt^-1
sum_j d_j w_j`` -- the stationary state at the period start
(covariance `K_0`, from the previous periods' noise) and this
period's per-step injections, independent of each other -- so
``Cov(dtheta) = dth K_0 dth^T + Gt^-1 D Gt^-T`` with ``D = sum_j
d_j Q_j d_j^T``.  Measured on the comparator-jitter sampler
(`_jitter_sampler`: a sawtooth of slope s_1 crossing a threshold
node with kT/C_n of noise): the turn-off crossing's sigma is
``sqrt(kT/C_n) / s_1`` -- 11.5873 ps measured against 11.5844 ps
analytic, 1.0002, and flat at 100 / 200 / 400 points -- the same
crossing motion that gives the held capacitor its
`(s_2/s_1)^2 kT/C_n`.  The reset edges of that fixture read
0.6437 ps, the threshold node's own faster slope there.

### `_coloured_covariance`

2026-09-25 (latest): the grid became ADAPTIVE (Simpson on the log axis,
`|S_2 - S_1|/15` per panel).  A fixed 40-per-decade grid read -10.5 % on a
driven Q = 20 tank and -21 % on a van der Pol's transverse covariance: the
response lines were narrower than the grid.  A first adaptive version
flagged on the plain rule's `|Q_h - Q_2h|` and refined smooth regions to
the point cap.  An oscillator's orbital lines are resolved first
(`_orbital_lines`), seeded only where the mode's adjoint couples.


2026-09-25 (later): the quadrature became the power law exact per interval
plus Richardson (`_power_law_weights`).  The trapezoid in ln nu was exact
for EF = 1 only (+2.8e-4 at EF = 2).  The product rule alone left the
response's bend at h^2, which telescopes only at EF = 1 (+1.7e-5 at
EF = 0.8).  Both measured before the fix; fourth order after.  A
stationary per-band colour (a Lorentzian) is integrated through its own
`CY(nu)` and unit sources; a modulated one is refused.


Built 2026-09-25 (log: "coloured noise in `covariance` and
`event_jitter`").  Before it, `covariance` and `event_jitter` refused any
coloured source through `_refuse_coloured`.

Why the frequency domain and not a shaping filter: the codebase's stance
(`coloured_diffusion`: "a SLOPE, NOT A STATE"), and against a fitted
Lorentzian ladder it has the exact power law, a hard band as
`sampled_variance` has, and the bordered event machinery for free
(`_forced_responses` is `PAC.solve`'s core).  No augmented Kronecker
closure.

Why the white part needs its own hook (`_lyap_cy`): the Lyapunov pieces
read `CY` at `w0`, which on a flicker circuit carries the flicker at the
clock frequency.  Folded as white, it counts the flicker a second time,
as a flat density.  Measured as a poison: +4.9 % on the RC gate.

Why the trapezoid in `ln(nu)`: exact on a pure 1/f density.  The linear
trapezoid on a log grid (`sampled_variance`'s) carries `(r - 1)^3 / 6`
per point: 6e-4 of a 1/f band at 40 per decade.  That was predicted and
then seen as the residual of the comparison with `sampled_variance`
(−5.6e-6 predicted at 400 per decade, −3e-6 .. −5e-6 measured).

### `covariance`

The docstring paragraphs before the move (the first, from the top of the docstring):

⚠ A GRID CHOSEN FOR `kT/C` IS NOT A GRID FOR THE PROFILE.  With
`CY` per step (2026-09-05) a switched capacitor's HELD variance
reads `kT/C` to 1e-4 at 1600 points and converges at better than
second order, while the TRACKING phase sits at this routine's
O(h/tau) floor -- 4% out at 800 points where the held value is
already 1.6e-4 -- and both agree with a reference simulator's sampled pnoise at
matched instants to 1e-3 (0.99878 track, 0.99915 edge, 0.99999
hold).  ⚠⚠ **AN EARLIER READING OF THAT TRACKING NUMBER IS
WITHDRAWN (2026-09-18).**  It said: "the tracked variance is 0.957
kT/C, NOT kT/C -- a sinusoidal clock holds the switch at full `gon`
only instantaneously, so the capacitor is never in equilibrium with
`Ron`; both tools agree on that independently."  That is a
DISCRETISATION FLOOR described as physics.  Within this fixture's
own definition the switch contributes `g V` and `white_noise(4 kT
g)` with the SAME `g`, so the variance obeys `dV/dt = -2(g/C)V +
2kT g/C^2`, for which `V(t) = kT/C` is an EXACT solution at every
instant and for ANY `g(t)` -- fluctuation-dissipation, and the
periodic solution is unique because `g > 0` contracts.  The
capacitor IS in equilibrium throughout.  Measured against that
exact profile (a stiff solve to periodic steady state, closing
error 0.0): the tracking value reads 0.7398 / 0.8467 / 0.9157 /
0.9556 kT/C at 200 / 400 / 800 / 1600 points -- errors 0.260,
0.153, 0.084, 0.044, HALVING per doubling, i.e. first order,
converging to kT/C.  ⚠ And the cross-tool agreement does not
rescue the claim: two tools discretising the same period at
comparable step counts agree about a SHARED artefact, which is
exactly why agreement between implementations cannot establish a
LIMIT -- the same lesson the edge-jitter work paid for when
cross-family agreement could not show a number was right.
✅ CONSISTENCY IS SEPARATELY CONFIRMED, and it is what this routine
should be judged on: against a closed-form time-varying reference
built by breaking the fluctuation-dissipation balance (an extra
white source not tied to `g`, giving a profile that spans 15x over
the period), `covariance` converges to the exact continuous answer
at FIRST order in both phases -- hold 1.07e-2 -> 1.09e-3 and track
2.60e-1 -> 4.44e-2 over 200..1600 points, ratios 2.07 and 1.90.
That is the O(h) its piecewise-constant injection predicts, and it
converges to the RIGHT limit.  ⚠ The "held converges at better
than second order" above is a property of the kT/C fixture, where
the exact profile is a CONSTANT and the leading terms cancel: with
the balance broken the held value converges at first order too.
⚠⚠ **AND THE kT/C ARGUMENT IS SPECIFIC TO A FIXTURE WHOSE NOISE IS
TIED TO ITS OWN CONDUCTANCE** -- the caveat is the reference suite's
and it is right.  It is a theorem about THIS circuit only because
the element's own definition contributes `g V` and `white_noise(4 kT
g)` with the SAME `g`.  A REAL DEVICE WOULD BREAK IT -- on a PSP
switch the measured `sid/(4kT g)` runs 1.09 during conduction to
3.17 through turn-off, so fluctuation-dissipation would not balance
and that circuit's tracking limit need not be kT/C.  ⚠ But that is
a caveat about OTHER fixtures, not about this one: the reference
side runs the SAME behavioural switch (`pcswitch.va`, the same `g`
and the same `white_noise(4 kT g)`), which is the point of the
fixture, so the theorem applies to both sides and there is no real
device near it.
⚠⚠ A REFERENCE-SIDE READING OF ~0.974 kT/C FOR THE TRACKING LIMIT
IS WITHDRAWN (2026-09-18, by the side that made it).  It came from
an Aitken extrapolation of 0.95489 / 0.96723 / 0.97128 / 0.97239,
whose error ratios DECELERATE -- 1.38, 1.14, 1.04.  A ratio heading
to ONE is a sequence that has stopped moving, and Aitken on a
stalled sequence returns approximately where it stalled rather than
a limit; a converging first-order ladder heads to TWO, as the one
above does (1.70, 1.82, 1.90).  So that number is a floor in the
reference's own sampled-noise integration, and the limit for this
fixture remains kT/C by the argument above.

⚠ AND THAT IS WHY IT REFUSES AN OSCILLATOR.  There `lambda_1 = 1`
gives `lambda_1^2 = 1`, so `I - M ⊗ M` is exactly as singular as
`I - M` — measured 3.1e-11 against 3.8e-11 — and the covariance
does not settle, it GROWS.  Variance linear in `t` is a random
walk, which is phase diffusion, which is the linewidth.  Demir 2002
gives the physical counterpart: an oscillator's output noise is
STATIONARY, not cyclostationary, because "noisy autonomous systems
cannot provide a perfect time reference".  There is no
cyclostationary object there to compute, and `oscillator_spectrum`
is the right route instead.

⚠ `CY/2` IS THE ONE-SIDED-TO-TWO-SIDED CONVERSION AND IT IS NOT
COSMETIC.  `CY` is a one-sided density (a resistor's `4kT/R`), so
the per-step injection is `Q_j = Jf_j^-1 (CY_j / 2h_j) Jf_j^-T`.
MEASURED against `kT/C` — exact, famously independent of `R` — on
an RC circuit, with and without the half:

    npts        100      200      400      800
    CY          1.861    1.928    1.963    1.981
    CY/2        0.931    0.964    0.982    0.991

The full-`CY` column converges to 2 and the halved one to 1, so the
factor is settled by the measurement rather than by argument.  The
residual halves per grid doubling — O(h), first order, which a
piecewise-constant approximation to white noise is.

⚠ AND THE GRID MUST RESOLVE THE NOISE BANDWIDTH, which is a real
precondition rather than an accuracy note.  The first attempt at
that gate read 0.517 because the RC pole at 159 kHz sat ABOVE the
grid's 100 kHz Nyquist: the discrete system genuinely does not
carry the noise the continuous one does.  A `kT/C` that comes back
low is the grid, not the code.

⚠ COST: the solve has `(2m)^2` unknowns and is dense here, so it is
`O(m^4)`.  Small circuits only until that is replaced.

⚠ ON A STAGED SOLVE (`state_events=True`) THE CLOSURE IS BORDERED
(2026-09-22): the noise moves the landed crossings, and the plain
closure on such a solve is wrong by O(1), not merely incomplete --
see `_event_closure` for the algebra and the comparator-jitter
sampler that measured 3.8x its analytic held variance unbordered
and 0.999 bordered.  Samples are the covariance at FIXED times.

### `sampled_noise`

The docstring before the move:

The one-sided PSD of the SAMPLE SERIES `y(t0 + kT)` -- DRIVEN
circuits, white AND coloured sources.  2026-09-15.

Returns `S` of shape `(len(times), len(freqs))` in `output`'s units
squared per Hz, for `0 < f <= f0/2`.  `sum` over the band of `S` is
the variance at `t0` that a sampler sees; see `sampled_variance`.
The instants actually used (the nearest period-grid points) are left
in `self.sampled_instants`.  ⚠ When comparing at a round instant,
pass the grid's own times (`pss.factored_period().times`): the grid
need not have the step count the `timestep` suggests (T/1000 can give
999 steps), and with a 1/f source the held variance moves measurably
between neighbouring grid points (peer report).

The docstring paragraphs before the move (from the end of the source-model paragraph):

A joint root of the summed `CY` made independent sources
non-additive (measured +9.3 % of the held variance, white + flicker
in one element).  ⚠ Mahmutoglu & Demir (TCAS-I 62(4), 2015) show
that a SWITCHED MOSFET's trap (1/f) noise is whitened below the
switching frequency and that a modulated-stationary 1/f model
over-predicts it; the physical fix needs trap states in the device
model, which is outside this analysis.  An agreement with another
tool on this quantity is agreement on the convention.  A flicker
spectrum is also singular at DC: keep `f` away from 0 (the band
integral takes an explicit `fmin`).

MEASURED (switched capacitor, Ron 1 k, Roff 1 G, 100 pF, 100 kHz,
gear, 400 points): the white held variance equals `covariance`'s
at the same instant to 1e-6 (0.998721 kT/C) and the tracking value
to 1.6e-4 (0.8466, covariance's own O(h/tau) floor); the seeded row
at `t0 = 0` equals `adjoint_transfer_row` to 1.5e-16; on an LTI
circuit the series sum equals the fold of `pnoise` over the same
sidebands to 1e-10, white and 1/f.  Each gate was checked to fail:
the source samples one step early read 1.30 kT/C, sidebands cut to
N/8 read 0.77 in track.

⚠ STAGE METHODS (radau, trbdf2) RUN NATIVELY, 2026-09-15: the source
enters every stage, so the sensitivities are read at the stage
abscissae and `CY` at the STAGE states (one re-traversal of the
orbit, cached).  Measured on the sampler: the LTI series sum equals
the pnoise fold to 1e-10 under both; the held variance reads
1.000000 (radau) / 0.999979 (trbdf2) kT/C at 400 points and stays
there at 800, the tracking value converges at first order (radau
0.949 -> 0.975, trbdf2 0.916 -> 0.957).  ⚠ `covariance` under these
methods WAS no reference at a switching edge until its injection
moved to the stages too (`_stage_injection`, 2026-09-15): the Van
Loan injection frozen at the END of each step read the held variance
0.876 / 0.934 (radau, 400 / 800 points); it now reads 1.3e-7 off at
400.  GLM period maps are refused.

⚠ TIME AVERAGE, PER FREQUENCY: the mean over `t0` of this PSD is
the fold of the time-averaged PSD, `sum_k pnoise(|f + k f0|)`, over
EVERY output band to the grid's Nyquist -- measured 4.4e-5 at 100
points (a fold cut at |k| <= 10 left 0.4 %).

⚠ TWO THINGS LIMIT A HELD VARIANCE, AND THEY ARE NOT THE SAME THING
(2026-09-21, from a peer's kT/C ladder on a half-on switch).  The
fold stops at `|n| <= maxsidebands`, so the source spectrum beyond
`F = (L + 1/2) f0` is dropped: for a Lorentzian that is the tail
`1 - (2/pi) atan(F/fc)` ~ `(2/pi) fc/F`, first order in the sideband
count.  `tail=True` adds, per instant and series frequency, the
`1/nu^2` extrapolation of the two OUTERMOST covered sidebands over
the uncovered half-lines -- `A = nu_L^2 dens_L`, `A / (f0 (F +
f0/2))` each side -- which is exact for a single pole and correct to
`(fc/F)^2` in general; nothing is fitted.  AND the covered sidebands
are computed on the grid: at `omega h ~ 3` per step (the top bands
of a fold to the grid's Nyquist) gear's discrete transfer is far
from the continuous one.  Measured on the LTI limit of the switched
capacitor (g = 0.5 mS, C = 100 pF, fc = 796 kHz, f0 = 100 kHz), the
deficit against the pure tail at a FIXED 100 sidebands: 3.03 /
1.93 / 1.33 / 1.12 / 1.05 at 204 / 400 / 800 / 1600 / 3200 points --
the peer's "coefficient 3.0" was this, at a fixed M/npts, not a
property of the kernel.  A run whose top sideband sits above
`omega h = 1` is WARNED (`SAMPLED_RESOLUTION_WARN`): the remedy is
more points; `tail=True` then closes what lies beyond the covered
edge, which needs that edge well above the spectrum's corner (edge
at 12 fc: radau 0.948 -> 0.998 x kT/C at 204 points, gear 0.947 ->
0.995 at 3200; edge at 0.8 fc: 0.71, the 1/nu^2 form is wrong
inside the corner).

### `sampled_variance`

2026-09-26: a per-band colour in the sample series was evaluated for EVERY
element at every point, per band frequency, sideband and instant (75 of
79 s for 38 frequencies).  Now it is the element alone, cached per
frequency, and classified once (stationary / separable / otherwise):
10 .. 41x, the values bit-identical or 2e-16.


2026-09-25: the linear trapezoid on the log grid was replaced by the log-log
rule (`_loglog_integral`), after measurement.  The linear trapezoid
overestimated a 1/f band by `(r - 1)^3 / 6` per point: +3.7e-4 .. +5.1e-4
at 40 per decade on the sampler.  The trapezoid in ln f, the first
proposal, overestimated a white band by +2.8e-4.  Log-log is exact for both
and for any power law: +6e-5 at 40 per decade, where the series PSD bends.
Table in the log, 2026-09-25.


The docstring paragraph before the move:

⚠ `fmin` AND `fmax` ARE REQUIRED.  With a 1/f source the integral
grows as `ln(fmax/fmin)` (measured 0.0010 kT/C per decade, flat from
1e-4 to 1e-1 f0, on a switched capacitor with a constant 1/f source)
and has no limit at `fmin -> 0`; with white sources only, the band
removes `fmin/(f0/2)` of the full variance because the series PSD is
flat.  Nothing is extrapolated into `[0, fmin]`.

### `jitter_metrics`

2026-09-25: `R_k` was a trapezoid on the LINEAR grid alone, whose first
interval spans the 1/f low end (R_0 +0.6 .. 0.8 % at `nfreq=601` on a 1/f
sampler, R_1..4 up to +1.8 %).  Split into `int S df` -- log-log on a log
grid joined with the linear one below their crossover, the plain
trapezoid above it -- plus `int S (cos - 1)` on the linear grid; `R_0 -
R_k` unchanged.  The log grid alone was tried first and moved a white
spectrum's R_0 by 3.5e-5 (its top intervals are 6 % wide).


The docstring paragraphs before the move:

MEASURED 2026-09-16 on a one-stage linear fixture built so the answer
is known exactly (`tau = RC = T`, LTI noise path, so `rho_k = e^-k`):
the transform reproduces `exp(-k)` at **1.0005 at every lag** k = 1..4
(0.368067 / 0.135404 / 0.049812 / 0.018325 against 0.367879 /
0.135335 / 0.049787 / 0.018316).  A Monte Carlo over 1176 noisy
crossings -- no PSS, no adjoint, no spectrum anywhere in it -- agrees
within 1 sigma at every lag it can resolve (0.346404 / 0.131894 at
k = 1, 2, i.e. 0.74 and 0.12 sigma), `sigma_t` matching at 0.9901.

⚠ `dc_rectangle` EXTRAPOLATES, WHICH THE REST OF THIS FAMILY REFUSES
TO DO.  `S` is known only on `[fmin, fmax]`; adding `S(fmin)*fmin`
assumes the series PSD is FLAT below `fmin`.  That is exact for white
sources and WRONG for `1/f`, where the integral has no limit as
`fmin -> 0` (see `sampled_variance`).  Hence off by default.
MEASURED: with it, `rho_k` is unchanged over a 100x range of `fmin`
(50 -> 0.5 Hz, identical to four digits); without it, `rho_k` drifts
with `fmin` exactly as truncation should and the drift GROWS with `k`
(0.9889 of analytic at k = 4, fmin = f0/2e4, recovering to 0.9994 at
f0/2e5).  ⚠ With a coloured source, LOWER `fmin` -- do not reach for
the rectangle.

⚠ `slew` IS A FINITE DIFFERENCE ON THE PSS GRID, central about the
instant actually used, and it converges at FIRST order (measured
0.9924 / 0.9962 / 0.9981 of the analytic slope at 200 / 400 / 800
points).  It is returned so a caller can check it rather than trust
it; every metric here is inversely proportional to it.

⚠ THE INSTANT IS THE CALLER'S, deliberately.  This does not hunt for a
threshold crossing: a threshold inferred from a simulated record can
be biased by startup, which moves the crossing off the steepest point
-- that cost half of an apparent deficit before it was caught (A8).
Pass the instant you mean; `instant` in the result is the grid point
actually used.

The comment on the slope guard:

⚠⚠ TWO EARLIER GUARDS HERE WERE WRONG, THE SECOND BECAUSE OF A
MISREAD NUMBER.  `slew == 0.0` is unreachable on a grid, so it was
no guard at all.  Replacing it with a threshold RELATIVE to the
steepest slope looked right only because I had read the peak's
central DIFFERENCE (-8.86e-07) as a SLOPE: over 2h = 5e-9 that is
-1.767e+02, which is 3.58e-04 of the steepest 4.94e+05 -- not the
1.8e-12 I inferred.  And on a 400-point grid 3.58e-04 is the
SMALLEST ratio any instant can have (a grid point never lands
exactly on the peak), so a relative threshold only ever describes
the grid and moves with N.

What does not move with the grid is the LINEARISATION this whole
family rests on: a crossing is displaced by delta_y/slew only while
that displacement stays small against the period.

### `_sampled_series`

The comment on where the source enters:

components, rolled to the INJECTION index: step j's source enters
at t_{j+1} (the reverse pass's own pairing) -- sampled one step
early the held variance read 1.30 kT/C instead of 0.9987
⚠ WHERE THE SOURCE ENTERS: an LMM step's source enters at
`t_{j+1}`, so `CY` is sampled at `x(t_{j+1})`; a stage method's
enters at every stage abscissa `t_j + c_k h`, so at the STAGE
states.  Measured on the sampler at 400 points: the end-of-step
state for every stage read the held variance 0.867 kT/C (radau)
against 1.000000 at the stage states, and it is the stage-state
value that holds under refinement.

The comment on the top sidebands:

⚠ THE TOP SIDEBANDS ARE COMPUTED ON THE GRID.  At `omega h > 1`
per step their discrete transfer is not the continuous one, and
the held variance is short by a factor that looked like a kernel
constant (3.0 x the pure tail on a half-on switch) and was this.

The comment on the bordered adjoint:

⚠ ON A STAGED SOLVE THE SAMPLE'S ADJOINT IS BORDERED (2026-09-22,
item 3 of the list after E7) -- the dual of the bordered forward
solve, as `adjoint_sideband_row`'s: the operator is the TOTAL
map's transpose, the sample is read at FIXED time (its costate
carries `dtheta/dx_0^T Pk_fixed[k0]^T d`), and the source's own
motion of the crossings enters as a third reverse pass carrying
`-zeta_k W_k` at the event nodes, `zeta = Gt^-T (g_theta + a
P_theta^T z)`.  Unbordered, the jitter sampler's held node had no
path to the threshold's noise at all.

### `oscillator_covariance`

The comment above `oscillator_covariance` before the move:

⚠⚠⚠ THE NOTE THAT WAS HERE ACCUSED THIS ROUTE AND WAS WRONG.  A
MONTE CARLO SETTLED IT THE OTHER WAY: this route is CORRECT and
`orbital_correlation` is the one that fails on an asymmetric orbit.
Direct SDE simulation of the variational system (trapezoidal, the
calibrated `Var(i) = CY/(2h)` injection, phase projected out every
step), sharing no Lyapunov solve and no modal sum:

    a      |R| MONTE CARLO  |R| modal    |Lyap| proj   MC/Lyap
    0.00   4.3629e-06       4.4515e-06   4.4437e-06    0.9818
    0.30   2.9992e-04       3.6913e-06   2.9986e-04    1.0002

`a = 0` is the CONTROL -- both routes agree there, so the MC had a
known answer to hit, and it did (2 %).  At `a = 0.30` it lands on this
route to 0.02 % and is 81x from the modal one.

⚠⚠ THE ARGUMENT THAT MISLED ME, RECORDED BECAUSE IT WAS PLAUSIBLE:
`|lam2|` FALLS with asymmetry, so relaxation gets FASTER, so the
transverse variance "should" shrink -- and the modal route did shrink
while this one grew 67x.  That reasoning is WRONG: asymmetry changes
the MODE SHAPES, so the noise projected onto the orbital direction
grows, and the variance rises DESPITE the faster relaxation.  A
physical argument is not a measurement.

The original (refuted) note follows for the record:
⚠ SUPERSEDED: THIS ROUTE DEPARTS FROM PHYSICS
ON AN ASYMMETRIC ORBIT, AND `orbital_correlation` DOES NOT.  van der
Pol + `a u^2`, sweeping `a` (orbit asymmetry 0 -> 0.41):

    a      |R| modal    |Lyap| proj  |K_orb raw|  |lam2|    amp
    0.00   4.4515e-06   4.4437e-06   5.0034e-06   0.882521  2.000
    0.05   4.4702e-06   4.7152e-06   4.7471e-06   0.881719  2.004
    0.15   4.9505e-06   1.7371e-05   2.7846e-05   0.874886  2.034
    0.30   3.6913e-06   2.9986e-04   9.5944e-04   0.844322  2.164

`|lam2|` FALLS (0.883 -> 0.844), so amplitude relaxation gets FASTER
and the transverse variance should get slightly SMALLER.  The modal
sum does exactly that (4.45e-06 -> 3.69e-06).  This route grows 67x
projected and 192x raw.  ⚠ The RAW covariance grows MORE than the
projected one, so it is not the oblique projection -- it is this
covariance.

⚠ FOUR EXPLANATIONS EXCLUDED BY MEASUREMENT, not by argument:
  * harmonic truncation -- the disagreement is FLAT at 9.877e-01 from
    `H = 4` to `H = 128`;
  * the projection's tangent proxy -- `|cos(u, tangent)| = 1.000000`
    at every asymmetry, against the Floquet phase mode;
  * the modal decomposition -- `|lam1| = 1.000000`, `lam2` real and
    well separated, one orbital mode, `p(T)-p(0) ~ 1e-14`,
    `q^T C p = 1.000000`;
  * a defect in `orbital_correlation` -- its two internal routes agree
    to 3.4e-04 independently of `a`.

⚠⚠ WHAT THIS DOES **NOT** INVALIDATE.  Every use of this function as a
reference in this file was on a SYMMETRIC orbit, where the two routes
agree to 0.3 % -- including A9 step 3's three-way gate and the C^2
biorthonormalisation defect it caught on 2026-09-07.  Those stand.
⚠ WHAT IS OPEN: which route is right is NOT settled.  The physical
argument favours the modal one, but that is an argument, and this file
does not close items on arguments.  A transient MONTE CARLO of the
orbital fluctuation is the decisive third route and has not been run.
Until then, treat this on a strongly asymmetric orbit as unvalidated.

The docstring before the move:

The state covariance of a FREE-RUNNING oscillator, split in two.

        Returns `(K_orb, d, info)`.  `K_orb` is the BOUNDED periodic
        (orbital) part of the covariance at `t = 0`; `d` is the growth per
        period along the orbit tangent, so

        ⚠ "BOUNDED" IS NOT "TRANSVERSE".  `K_orb` has the SECULAR growth
        removed and still contains the phase direction's bounded
        within-period variance.  Demir's orbital deviation `y` is the
        OBLIQUE projection `v_1^T y = 0`, so the transverse covariance is
        `Pi K_orb Pi^T` with `Pi = I - u v^T/(v^T u)` -- which is what
        `orbital_correlation`'s eq (23) sum equals (to 1e-4), and what
        `K_orb` itself does NOT equal (2-6 %, falling as 1/Q).  Read
        `K_orb` as the bounded part; project it if you want `R_yy(0)`.
        See `orbital_correlation`.

            K(t_0 + n T) = K_orb + n d u u^T

        exactly, for every integer `n`, with `u` the pair-space tangent
        scaled so its first block is `xdot(0)`.

        ⚠ WITH `samples=True` THE SPLIT MOVES WITH THE ORBIT, AND IT IS
        WORTH SAYING PRECISELY BECAUSE THE OBVIOUS READING IS WRONG.
        `info['orbital_samples'][j]` is `P(t_j)`, the solution started from
        `K_orb` at `t = 0`, and it satisfies

            K(t_j + n T) = P(t_j) + n d u_j u_j^T,   u_j = Phi(t_j, 0) u

        so `P` is periodic UP TO the growth -- `P(T) = P(0) + d u u^T`, not
        `P(T) = P(0)`.  The walk is along the orbit and the orbit turns, so
        the growth DIRECTION is the propagated tangent rather than a fixed
        `u`.  MEASURED: `P(T) - P(0)` matches `d u u^T` to 3.2e-15, and the
        full prediction holds to 2.5e-09 against a brute-force recursion
        run forty periods (9,600 steps) from `K = 0`.

        ⚠ THIS IS THE OBJECT `covariance` REFUSES TO RETURN, and the
        refusal was right: there is no periodic solution, so anything that
        returned one number would be hiding the physics.  `lambda_1 = 1`
        gives `lambda_1^2 = 1`, so `I - M kron M` is exactly singular --
        MEASURED here at `sigma_min` 2.3e-11 with the next singular value
        at 0.997, i.e. a null space that is cleanly ONE-DIMENSIONAL and
        spanned by `u kron u`, with left null `v kron v`.  So it borders
        exactly as the PPV and the deflated PAC solve do, and the border is
        the pair the rest of this class already computes.

        ⚠ THE SPLIT IS NOT A NUMERICAL DEVICE, IT IS THE ANSWER.  Demir
        2002: an oscillator's noise is STATIONARY, not cyclostationary,
        because "noisy autonomous systems cannot provide a perfect time
        reference".  `K_orb` is the part a designer can read as an
        amplitude/orbital noise -- it settles, it is periodic, it is
        finite.  `n d u u^T` is the random walk ALONG the orbit, which
        never settles and which no periodic object can hold.  Reporting
        only their sum at some finite time is what makes an oscillator
        covariance look divergent and useless; reporting the parts makes
        both usable.

            [ I - M kron M    u kron u ] [ vec(K_orb) ]   [ vec(K_1) ]
            [ (v kron v)^T        0    ] [     d      ] = [     0    ]

        ⚠ AND `d` HAS A CLOSED FORM THAT NEEDS NO KRONECKER AT ALL.
        Left-multiplying the first row by `(v kron v)^T` kills the singular
        block, leaving

            d = (v^T K_1 v) / (v . u)^2

        an `O(n^2)` contraction against the `n^4` solve.  Both are computed
        and `info['d_residual']` is their relative difference; they are the
        same quantity by construction, so a disagreement means the border
        pair is wrong rather than that one route is less accurate.

        ⚠ `(v . u)` IS NOT 1 AND ASSUMING IT IS COSTS A FACTOR OF 2.3.
        `ppv()` normalises on the FIRST BLOCK, `v[:m] . xdot = 1`, which is
        the normalisation a state perturbation entering the first block
        sees -- an injected current, and what every other shipped path
        does.  The FULL PAIR contraction is a different number: 0.663 on
        van der Pol, so `(v . u)^-2 = 2.28`.  That exact mistake produced a
        2.31x discrepancy that was chased as a code defect for a while; it
        is why `d` is written with the pair inner product spelled out.

        ⚠ `d` ALONE IS MEANINGLESS WITHOUT PINNING `u`'s SCALE.  Rescaling
        `u -> s u` sends `d -> d / s^2`, so only the PRODUCT `d u u^T` --
        returned as `info['growth']` -- is an invariant of the circuit.
        `u` is pinned here by `C u = q`, the same condition `ppv()` uses to
        scale the tangent, which makes its first block exactly `xdot(0)`
        and gives `d` its physical reading below.

        ⚠ WHICH MAKES `d / T` A COMPLETELY INDEPENDENT ROUTE TO THE
        DIFFUSION CONSTANT, and that is this method's real gate.  A phase
        deviation `alpha` displaces the state by `alpha u`, so the growing
        covariance is `Var(alpha) u u^T = c t u u^T`, giving `d = c T`.
        The two computations share only the `CY/2` convention: `c` is a
        quadratic form in the ADJOINT-replayed PPV, while `d` comes from a
        FORWARD Lyapunov recursion closed by a bordered Kronecker solve.
        Neither touches the other's machinery.

        ⚠ AND THE TWO ANCHORS BEHIND THEM ARE ALSO INDEPENDENT, which is
        the property that was missing when a 2x error survived a 0.9965
        agreement.  `covariance`'s injection is anchored to `kT/C`;
        `diffusion_constant` is anchored to a nonlinear Monte Carlo reading
        phase from zero crossings.  `info['c_from_growth']` against
        `diffusion_constant` therefore closes a loop between two separately
        anchored quantities rather than reproducing one of them.

        ⚠ COST: the bordered solve has `(2m)^2 + 1` unknowns and is dense,
        so it is `O(m^4)` like `covariance`.  Small circuits only.  The
        closed form for `d` is cheap; pass `samples=False` and read
        `info['c_from_growth']` if the orbital part is not wanted.

The comment on the staged oscillator:

a staged oscillator closes on the TOTAL map with the crossings'
noise-driven motion in the injection (2026-09-22) -- the same
`_event_closure` as `covariance`, whose `u`, `v` below are the
total map's already

### `oscillator_edge_jitter`

The docstring paragraphs before the move (from the `projection_share` sentence):

cost you here.  ⚠ IT IS A PROPERTY OF THE SOURCE MIX, NOT OF THIS
METHOD: measured 0.1611 / 0.0173 / 0.0017 / 0.0002 as the TANK's noise
falls 1e-6 -> 1e-9 against a fixed 1e-6 at the buffers, because the
phase direction is exactly what tank noise drives.  So a fixture whose
oscillator is quiet will show the projection doing nothing and teach
you the wrong lesson; one with a noisy tank shows 16 %.  It also falls
as 1/Q, so a high-Q fixture hides it too.

⚠ THE SLOPE IS A LOCAL QUADRATIC FIT, NOT A TWO-POINT DIFFERENCE, and
that is a correction rather than a preference.  A straddling
difference across a threshold crossing samples whichever pair of grid
points brackets it, and the crossing sits at a different fraction of a
step on every grid: measured -0.720 / -0.077 / +2.028 % under
refinement, changing SIGN, while the period and swing of the same runs
converged cleanly at first order.  The quadratic fit reads
1.527458 / 1.527793 / 1.528023 at 240 / 480 / 960 points -- flat.  The
straddling value on a 240-point grid was 2.8 % low, and since every
quantity here goes as `1/s^2` that inflated the variance by 5.6 % and
was briefly blamed on the other side's integrator.

MEASURED 2026-09-17 on a van der Pol tank driving three tanh buffers
(`tau_buf = RC = 0.5`), noise at every node: `2 sigma_t^2` reads
1.055708e-06 / 1.055112e-06 / 1.062936e-06 at 240 / 480 / 960 points,
i.e. FLAT to a few parts per thousand.  ⚠ An earlier record of this
said it CONVERGED (1.023870e-06 -> 1.055112e-06 -> 1.062936e-06, "the
increments shrinking ~4x").  That was the slope error above shrinking
with the grid, not the covariance converging -- with the slope taken
at the requested instant the grid dependence is essentially gone.  A
Monte Carlo over EIGHTY
seed-runs -- a noisy transient with no PSS, no adjoint and no Lyapunov
solve in it -- gives

    MC / analysis = 1.0066 +/- 0.0102   (0.64 sigma from 1.000)

over 124 seed-runs across three grids, ONCE the comparison is made
between the same orbit on both sides.

⚠⚠ AN EARLIER RECORD OF THIS CLAIMED A 4-SIGMA DEFECT IN THIS METHOD
AND IT WAS WRONG -- THE FAULT WAS IN THE COMPARISON, NOT HERE.  An
80-seed campaign reported `1.0571 +/- 0.0141`, "under-predicts by
5.7 %, grid-independent, therefore on the analysis side".  It was
neither.  `sigma_t = sqrt(var)/slew`, so the comparison goes as
`1/slew^2` -- and the Monte Carlo finds its crossings on an EULER
orbit while this method divides by the slope of the PSS's GEAR orbit.
At finite `h` those orbits differ: Euler's period runs 0.39 % short
and its slope at the crossing is low by 2.717 / 1.361 / 0.686 % at
240 / 480 / 960 points, HALVING as first-order convergence requires.
Squared, that predicts ratios of 1.0566 / 1.0278 / 1.0139.

MEASURED AGAINST THAT PREDICTION (the 960 value pinned before its
seeds ran): raw 1.0493 / 1.0513 / 1.0165, and dividing by each grid's
own independently measured slew correction collapses them onto
0.9931 +/- 0.0153, 1.0229 +/- 0.0160, 1.0025 +/- 0.0261 -- pooled
1.0066 +/- 0.0102, with a residual trend of +0.0047 per doubling
against a per-grid scatter of 0.019.

⚠ THE ERROR THAT PRODUCED THE FALSE ALARM IS WORTH MORE THAN THE
NUMBER: "grid-independent" was asserted from TWO points 0.99 sigma
apart.  Two noisy points cannot tell FLAT from HALVING, and halving is
what it was doing.  An absence of evidence was read as evidence of a
property, and a structural conclusion ("therefore the analysis's") was
built on it.

A factor of 2 in the variance remains excluded at more than 15 sigma.
⚠ When validating this against a transient, run the Monte Carlo on the
SAME integrator as the PSS, or divide by the slope of the orbit the
Monte Carlo actually runs on -- otherwise the mismatch enters squared.

The comment on the slope:

⚠ THE SLOPE IS TAKEN AT THE REQUESTED INSTANT, NOT AT THE SNAPPED
GRID POINT, and the difference is the whole correction.  `time` is
typically a threshold crossing, which sits at a different FRACTION
of a step on every grid; differentiating at the nearest sample
instead reproduces the straddling value (1.485472 against a
converged 1.5275 on a 240-point grid, 2.8 % low) and every quantity
here goes as 1/s^2.

The comment on the PPV samples:

⚠ `samples[j]` IS node j (2026-09-20, measured).  This list used to
be `[v0] + samples`, which paired node j's covariance with the
phase vector of node j - 1 -- a one-node shift, first order in the
step: the obliquely-projected transverse cycle mean read 8.0e-2 /
4.0e-2 / 1.7e-2 against its own N = 3200 value at N = 200 / 400 /
800 (halving), and 3.3e-4 / 5.1e-4 / 1.8e-4 with the list
unshifted.  The same prepend sat in three tests, and it was read
as "the Lyapunov route is first order" (5129485) -- it was this.

### `orbital_mode_weights`

The docstring paragraphs before the move (from the top, after its first line):

⚠⚠⚠ READ THIS FIRST: THE BASIS OMITS THE ANNIHILATED MODES, AND WHAT
THEY CARRY IS A FLOOR NOTHING BELOW CAN GO UNDER.  `floquet_modes`
returns the NON-NULL directions, so `sum cw[k,k'] u_k u_k'^H` reproduces
only the part of `K_orb` that lives on them.  How much that is depends
entirely on WHERE THE NOISE ENTERS -- measured on `_osc_with_ladder`'s
circuit at `nslow = 4`, moving one current source and changing nothing
else::

    injected at            ||K_orb||    reconstruction residual
    the oscillator node    2.70e-05     1.80e-03   (0.18%)
    a SLOW ladder node     3.94e-01     3.56e-01   (36%)
    a FAST ladder node     6.87e+02     9.996e-01  (99.96%)
    a faster one           3.33e+03     9.999e-01  (99.99%)

**When the injection lands in a fast branch the non-null modes capture
essentially NOTHING of the covariance.**  The annihilated modes are
killed by the period map, so they enter the stationary covariance only
through the `j = 0` term -- but that term is not small when the noise
is injected there, and THAT IS WHERE DEVICE NOISE ACTUALLY IS: every
resistor in a bias or tuning network.

⚠ SO A MODAL ORBITAL SPECTRUM BUILT ON THIS BASIS IS COMPLETE ONLY FOR
NOISE THAT ENTERS THE SLOW SUBSPACE, and the suite's own gate on this
(`rel < 1e-2`) holds because its fixture injects at the oscillator
node.  That is a property of the fixture, not of the method.

⚠ AND THE RESIDUAL IS A DETECTOR, NOT A TRUNCATION BOUND.  It catches a
DROPPED NON-NULL MODE well -- which is what the note below claims for
it -- but it SATURATES at the floor above, so it cannot certify a
truncation below whatever the null modes carry, however many modes are
kept.  Independently reproduced by a peer session on a different
oscillator with a different `K_orb` route (69% there, mechanism
identical, magnitude not transferable).

⚠ **NOT THE SPECTRUM.** `S_yy` additionally needs the Fourier
coefficients of the periodic parts (`floquet_modes` returns them
as `p`/`q`) and the resolvent `1/(i(j−j')ω₀ − μ_l' − μ_l*)` of
eq (22), and the OUTPUT spectrum then needs a layer that can carry
contributions asymmetric about the carrier. Those are not built.

### `ORBITAL_ASYMMETRY_LIMIT`

The comment before the move:

Half-wave asymmetry above which `orbital_correlation` is known to be
wrong.  MEASURED (below); 0.02 is a decade inside the smallest
asymmetry at which the error was already visible.

### `_warn_if_orbit_is_asymmetric`

The docstring paragraph before the move:

⚠⚠ THIS WARNING WAS WRITTEN FOR A DEFECT THAT IS NOW FIXED, and kept
for the residual.  On 2026-09-07 `orbital_correlation` read 81x LOW
against a Monte Carlo on van der Pol + `a u^2` at asymmetry 0.41.
The cause was in `floquet_modes`: the replayed adjoint is `C^T q`,
not `q`, and was used untransformed -- invisible on a unit-reactance
symmetric orbit, catastrophic off-axis where the two adjoints are
nearly parallel.  With the `C^-T` transform applied, against the
same Monte-Carlo-validated Lyapunov reference at `a = 0.30`:

    npts    eq22 / Lyapunov
     400       1.0595
     800       1.0300
    1600       1.0151

halving per doubling -- an `O(h)` DISCRETISATION residual of the
adjoint replay, converging to 1, not a defect.  At `a = 0` it is
1.0001.  So this warns that the residual is grid-limited on such an
orbit and says how to shrink it; it no longer says the answer is
wrong, because it is not.

### `orbital_correlation`

The docstring paragraphs before the move:

⚠ `CY/2`, NOT `CY`.  The library's `CY` is one-sided; eq (22)
integrates `B B^T` as a two-sided intensity.  Consistent with the
`kT/C`-calibrated Monte Carlo injection `Var(i) = CY/(2h)` in this
file's record, and confirmed here by three routes agreeing.

⚠⚠ GATED THREE WAYS, because a modal sum transcribed from an image
of an equation is exactly the object this record distrusts.  On
van der Pol under gear: (i) this sum against `R_yy(0)` evaluated
from its DEFINITION as a 1-D Lyapunov integral along the orbital
mode, no Fourier machinery — agree to 3.5e-4; (ii) both against
the CYCLE-MEAN transverse part of `oscillator_covariance`'s
samples (`P(t_j) - (t_j/T)·growth_samples[j]`), which shares no
machinery with either — magnitude to < 1e-3.  That third route is
what found the state-block scale defect in `floquet_modes`.

⚠ THE REFERENCE IS THE CYCLE MEAN, NOT `K_orb(0)`.  Lemma 3.5's
`R∞_yy` depends on `τ` only — the stationary part.  At `t = 0`
van der Pol's amplitude direction is pure-v while this is
isotropic, which is a rotating radial direction averaged over a
cycle, not a disagreement.

✅ THE 2-3 % SHAPE RESIDUAL WAS THE REFERENCE, NOT THIS SUM -- closed
2026-09-04.  Subtracting only the SECULAR growth `(t/T) d u u^T`
from the Lyapunov samples leaves the phase direction's BOUNDED
within-period variance, which eq (22)'s `l >= 2` sum correctly
excludes.  Demir's `y` is defined by the OBLIQUE projection
`v_1^T y = 0`; project the samples with `Pi = I - u v^T/(v^T u)`
and the three-way agreement is 5.9e-4 / 6.0e-4 / 3.1e-4 / 1.5e-4
at Q = 4 / 8 / 16 / 32 (euler-plain, n = m), improving with
refinement -- quadrature.  The old residual fell as 1/Q_lambda
(5.4 / 2.4 / 1.2 / 0.6 %), which is orbital variance ~ Q against a
constant phase-bounded part: the same fact, seen from the sweep.
⚠ A candidate recorded earlier -- the phase-orbital CORRELATION's
tau = 0 value -- was DISPROVED from eq (18a) before it was built:
at tau = 0 its brace is {1 - 1} = 0 identically, and (23) states
R_yy(0) = sum C_lhj alone.  Named so nobody rebuilds it.

The comment on the phase mode:

the phase mode by its tangent alignment (see `_phase_mode_split`),
and NEVER in the orbital sum: swept in with a near-zero exponent it
blew up as 1/|mu|^2 (146x .. 2449x radau's R on a 3:1 gear grid)

### `orbital_spectrum`

The docstring paragraphs before the move:

⚠⚠ MEASURED 2026-09-14 -- "negligible" is a property of THOSE
circuits, and the over-statement can be large.  Against pnoise (the
total linear sideband noise, on one absolute scale with the certified
Lorentzian), `(up+lo)/(4 (S_ph + S_orb))` at 10 f_amp on van der Pol
with C = 4, Q = 8 and an `a u^2` asymmetry:

    half-wave asymmetry   0      0.017   0.033   0.067   0.100
    R                     0.999  0.997   0.972   0.691   0.309

So on a SYMMETRIC orbit the amplitude is right to 0.1 % (also at
C = 1 and Q = 50 -- the first external check this spectrum had), and
on an asymmetric one the sum over-states by up to 3.2x (5 dB) with no
grid dependence.  A Monte Carlo of the SDE (64 oscillators x 4000
periods, trapezoidal, `Var(i) = PSD/(2h)`) settles which side is
right: at a = 0.30, MC/pnoise = 1.011 and MC/(S_ph + S_orb) = 0.313,
with the a = 0 control reading 1.009 / 1.008.  pnoise is the total.

⚠⚠ THE DROPPED CROSS TERM IS THE CAUSE -- BUT ONLY WITH EVERY
HARMONIC KEPT (2026-09-15; the reverse reading of 2026-09-14 held only
for the truncated form).  `S_corr` built from eq (92) on this fixture
is ~1e-8 of the total: that form keeps only the PPV's DC harmonic AT
THE NOISE SOURCE's row, and an ideal tank inductor shorts that node at
DC (`vbar` = [-2.3e-6, -0.114]).  The full-harmonic correlation is
-1.1 to -2.4x this spectrum at a = 0.30, and `modal_spectrum`'s three
terms sum to pnoise (1.006 at 10 f_amp, grid error that halves with
the grid).  The orbital mode's AM share at the output,
`sin^2 arg(U_{l,1}/U_{0,1})` = 0.307, is the flat factor below.
The true total is BELOW even the phase Lorentzian alone (pnoise/S_ph
= 0.61 at a = 0.30): the over-statement is in the decomposition's
frequency-independent terms above f_amp, not in a missing
correction.  LOCALISED: the PHASE half is the Lorentzian's
frequency-independent PPV -- with `c(f)` from `frequency_aware_ppv`
pnoise's PM content matches it to <= 2.3 % at 0.3-10 f_amp -- and
the ORBITAL half over-states by a factor FLAT in offset (AM content
0.317 of this spectrum at a = 0.30, at every offset), which is open.
Traversa & Bonani's own Figs 1-2 show the same limit on their
amplitude-phase-coupled test oscillator (theory above the exact
spectrum at high frequency, growing with the coupling).  A warning
fires above `ORBITAL_ASYMMETRY_LIMIT`.

The comment on the asymmetry warning:

⚠ NOT the grid residual `_warn_if_orbit_is_asymmetric` names:
measured against pnoise (2026-09-14), the sum this spectrum is
meant for over-states the TOTAL on an asymmetric orbit, and no
refinement changes it -- see the docstring.

The comment on the missing line:

⚠⚠ NO ORBITAL LINE AT THIS HARMONIC -- refuse rather than return the
tails of the others (2026-09-14).  The line weight at `j f0` is
`W_j = sum_{l,h} Re(row C_lhj row)`; where it is zero (a symmetric
orbit's even harmonics: the modes' own Fourier content vanishes)
what this would return is the neighbouring lines' Lorentzian tails,
which a Monte Carlo put 3.2x LOW at 2 f0 on van der Pol C=4 Q=8
(pnoise agreed with it to 1 %).  ⚠ NOT caught: DC, where a small
line can exist and the model read ~100x HIGH (the tank inductor
shorts the node, which Lorentzian tails do not know), and 2 f0 on
an asymmetric orbit (0.40) -- away from the fundamental use
`pnoise`.

### `modal_spectrum`

2026-09-26: a source whose LEVEL FOLLOWS THE ORBIT is taken instead of
refused (Andreas: "modal_spectrum with orbit-varying noise").  White parts
in the P-form `sum_{m,m'} T_m (P_{m'-m}/2) T_{m'}^H` (exact, no root: a
root of `(k V)^2` is `|k V|`, whose kink the sideband window cuts); each
coloured component as a unit process through its own columns `G(t)`,
rows the harmonics of `q_l^T G`, band `p` at `|w - p w0|`.  Gated against
the same physics as a stationary source times the modulating voltage:
every part ~1e-13 (white, signed flicker, a per-band Lorentzian); the
moving shape against `pnoise(cyclostationary=True)` 1.2e-13 on radau.
The stationary paths are bit-identical.  Log, 2026-09-26.


The docstring before the move (its first-line tail and its `WHY IT EXISTS` and near-carrier paragraphs):

transfer, which sum to the total.  E6, built 2026-09-15.

⚠⚠ WHY IT EXISTS.  `oscillator_spectrum(frequency_aware=False) +
orbital_spectrum` over-states an asymmetric orbit's total by up to
3.2x (Monte-Carlo-confirmed).  The missing piece IS the correlation --
but not in the form Traversa & Bonani keep: their eq (92) retains only
its DC harmonic, ~1e-8 of the total on van der Pol (the tank inductor
shorts the source node at DC).  With every harmonic kept it is -1.1 to
-2.4x the orbital term there.  It removes the DC-PPV phase excess above
f_amp AND the orbital mode's PM projection at the output -- the orbital
line's AM share `sin^2 arg(U_{l,1}/U_{0,1})` (0.307 at a = 0.30) is the
flat factor `orbital_spectrum` was measured to over-state by.
Measured on van der Pol C=4, Q=8, 400 points per period (H = 8, 16
sidebands), `total / (pnoise/2)`:

    a     harmonic   +1      +3      +10     -3      -10  f_amp
    0.00  1          1.0006  1.0006  1.0006  1.0006  1.0006
    0.00  2                  1.0010  1.0010  1.0010  1.0009
    0.30  1          1.054   1.014   1.006   1.008   1.006
    0.30  2                  1.014   1.005   1.008   1.005

⚠ THE a = 0.30 EXCESS IS GRID ERROR, NOT THE MODEL: on 800 points it
reads 1.026 / 1.006 / 1.002 at +1/+3/+10 f_amp (it halves -- the O(h)
of the modes), and the a = 0 control reads 1.00013; H 8 -> 12 and
sidebands 16 -> 24 move the fourth digit.  ⚠ Because the correlation
cancels most of the other two, a few percent of error in any part is
AMPLIFIED in the total -- which is why the three are computed together
here rather than the correlation being offered as an add-on to
`oscillator_spectrum + orbital_spectrum`: those line-shape spectra keep
only the resonant term of each line (2.5 % short at 10 f_amp even on a
symmetric orbit), which is harmless alone and not under cancellation.
The modal sum also reproduces pnoise's upper/lower sideband asymmetry,
which the two-term sum cannot.

Near the carrier `phase` IS the library Lorentzian (1.00022 of
`oscillator_spectrum(frequency_aware=False)` from 0 to 3 linewidths,
symmetric orbit) and `correlation` is ~1e-6 of it.  Above f_amp
`total` agrees with `pnoise`; within the linewidth pnoise has no
meaning and this is the route.

The comment on the phase mode:

the phase mode by its tangent alignment, on any grid -- see
`_phase_mode_split` (the 1e-6 window refused every gear solve on a
non-uniform grid)

The comment on the Fourier coefficients:

⚠ `_period_dft`, not an index DFT: measured 8-13 % off and NOT
converging on a 3:1 grid before (see `PSS._period_quadrature`)

### `_phase_mode_split`

2026-09-26: `waveform` is FULL width (the reference row is in it), and the
tangent was computed from it as if reduced, a second zero inserted, so
every unknown past the reference was read one slot late.  Van der Pol hid
it (the inductor current it misread is ~0 at the phase anchor); an idle
3 V source in the circuit read the phase mode's alignment as 0.57 and
refused.  Now `_orbit_states`.  No existing result moved (the same mode
was picked).


The docstring before the move (after its first line):

DEFINES IT, its eigenvector's alignment with the orbit tangent, not by
a window on `|lam| - 1` (2026-09-20, Andreas: gear as a first-class
choice on non-uniform grids).

On a uniform grid the phase multiplier is 1 to rounding.  On a grid
whose step varies, a multistep or trapezoidal solve loses time-
translation symmetry and the multiplier leaves the circle at O(h^2)
-- measured 1 - 5.1e-05 (gear) and 1 + 4.7e-05 (trap) at 400 points
on a 3:1 grid; radau keeps it at 1 + 1e-11 -- and a 1e-6 window
refused every gear solve there.  The right eigenvector of the phase
mode is the tangent `xdot(0)` (`C xdot = -i(x)` for the autonomous
circuit), which no orbital mode shares, so alignment picks it on any
grid; its exponent is then forced to 0 exactly, as the consumers
already do.  The departure is WARNED with its size when it exceeds
rounding, and the split is REFUSED when a second multiplier lies
within ten times that departure of the circle with any alignment --
the case a window ever protected against.

### `diffusion_constant`

2026-09-26: a MODULATED white source is taken: `CY` at each PPV sample's
state, Demir's `B(x(t))` form, where `_cy_reduced` refused it.  Against
the Lyapunov growth route (`CY` per step, no PPV): radau 4.6e-13 (van der
Pol, `(k V)^2`) and -4.4e-11 (`test_multiplicative_noise_...`'s fixture;
gear +1.0e-3, its own gap).  The Ito/Stratonovich reason the refusal also
carried concerns the DRIFT; `c` is the diffusion, where they agree.


The docstring paragraphs before the move:

⚠ ITS SCALE WAS WRONG BY 2x AND IS NOW FIXED -- kept because the
way it survived is the instructive part.  `diffusion_constant` used
the full `CY` while `covariance` used `CY/2`: two functions in this
class disagreeing about whether `CY` is one- or two-sided.  It was
validated against a Monte Carlo injecting `Var(i) = CY/h` per step
and agreed to 0.9965 -- because that Monte Carlo carried the SAME
hot convention.  A measurement built on the assumption under test
cannot test it.

SETTLED AGAINST `kT/C`, which is external to both: an injection of
`Var(i) = CY/h` reproduces 1.92x `kT/C` over ten independent runs
(1.75-2.04).  With `CY/2` throughout, `diffusion_constant` gives
7.9516e-08 against a correctly scaled Monte Carlo at 7.7083e-08 --
ratio 1.0316, inside that measurement's 4.1% uncertainty.

⚠ AND A SECOND DISCREPANCY WAS NOT A CODE DEFECT AT ALL.  Two Monte
Carlo routes disagreed by 2.31x, which looked like a third error.
It was in the DIAGNOSTIC: `ppv()` normalises on the FIRST BLOCK
(`v[:m] . xdot = 1`), which is right for a perturbation entering
the first block -- an injected current, and what every shipped path
does -- but wrong for contracting against a full PAIR deviation,
where the factor is `1/(v . u_pair) = 1.508`.  Correcting it turned
a 2.13 variance ratio into 1.07.  The sign difference alongside it
is a convention, not an error: a later zero crossing means DELAYED,
while projecting onto the tangent makes positive mean ADVANCED.

### `_period_weights`

The docstring paragraph before the move:

⚠ EVERY PERIOD INTEGRAL HERE USED `h = diff(times)` -- the LEFT
RECTANGLE rule (2026-09-20).  On a uniform periodic grid that is the
trapezoid rule, spectrally accurate; on a NON-UNIFORM grid it is
FIRST order by itself (its error is (1/2) integral h'(t) y(t) dt, not
zero).  Measured on a solved-history (gear) run of an index-2
oscillator on a smooth grid, `c` against radau: -1.9e-3 / -9.3e-4 /
-4.6e-4 at N = 200 / 400 / 800 with the rectangle weights and
+4.0e-5 / +3.6e-5 / +1.3e-5 with these, from the same samples.  On a
uniform grid `0.5 h + 0.5 h == h` exactly, so every uniform-grid
number is bit-identical to before.  The one-step kinds' replays are
uniform-grid replays (`factored_period_stage`), so only the
solved-history kind ever paid this.  ⚠ AND THE TRAPEZOID IS ITSELF A
SECOND-ORDER CAP on a smoothly varying grid (2026-09-21): with `pss`
given, a non-uniform grid and no landed events, these are the
periodic cubic-spline weights of `periodic_spline_weights` (radau's
`c` on a 1 + 0.5 sin grid 1.3e-5 -> 2.4e-10 at N = 200); under
landed events the spline breaks at their nodes (2026-09-21, it used
to keep the trapezoid), and a uniform grid is unchanged.

### `_white_diffusion_at`

The docstring line before the move:

value `diffusion_constant` used to return silently.  `phase_psd`

The comment on `samples_eq`:

⚠ `samples_eq`, NOT `samples`.  `CY` is an EQUATION-ROW
covariance and `samples` is `C^T v_1`; contracting it here made
`c` wrong by `C^2` on the differential rows and exactly zero on
the algebraic ones.  See `_equation_row_ppv`.

The comment on the index-2 constraint:

⚠ A NOISE SOURCE ON AN INDEX-2 CONSTRAINT GIVES c = 0, SILENTLY
(2026-09-21): a voltage noise in series with a DC source inside a
capacitor loop perturbs an algebraic constraint -- a DIFFERENTIATED
input, whose response is a charge jump the PPV projection cannot
represent -- and `c` came back exactly 0 for every method on an
index-2 van der Pol (the same circuit's current noise at the node
gives 3.2e-9).  Named here once: the PPV's algebraic fallback
fired (index >= 2) and `CY` has power on an algebraic row.
(the gate is the index-2 condition itself -- `G[A,Z]` singular at
the orbit point, the test the PPV's algebraic fallback makes --
computed here so every kind is covered, not only solved-history)

The comment on `cy/2`:

⚠ `cy/2`, THE SAME ONE-SIDED-TO-TWO-SIDED CONVERSION `covariance`
USES.  `CY` is a one-sided density (a resistor's `4kT/R`), and
these two functions disagreed about it until a Monte Carlo was
run against `kT/C`: an injection of `Var(i) = CY/h` per step
reproduces `1.92x kT/C` over ten independent runs (1.75-2.04),
so that convention carries TWICE the physical noise power.
`covariance` was already right; this was not, and its agreement
with a Monte Carlo built on the SAME hot convention is exactly
why the error survived.

### `colour_projection`

The docstring paragraphs before the move:

⚠ AND USING THE QUADRATIC ONE FOR A COLOURED SOURCE RETURNS A
PLAUSIBLE NUMBER, NOT AN ERROR.  It is never zero where the white
answer is not, so nothing downstream would look wrong.  The
measured separation is 22 ORDERS on van der Pol -- `c = 7.95e-08`
against `Gamma = 1.9e-29` -- so the two functionals are not close
approximations of each other and cannot be substituted.

⚠ TWO INDEPENDENT MECHANISMS FORCE `vbar` TO ZERO, AND ONLY ONE OF
THEM IS THE ONE DESIGNERS KNOW.  Measured on an LC oscillator,
sweeping an even term `a (u^2 - 2)` in the nonlinearity and a
series tank resistance `Rs`:

    a      Rs      Gamma/c
    0.00   0.00    2.4e-22
    0.00   0.20    9.7e-23
    0.25   0.00    4.9e-23
    0.25   0.05    2.1e-04
    0.25   0.20    4.1e-03

NEITHER ASYMMETRY ALONE NOR LOSS ALONE UPCONVERTS.  `c` is
7.9e-08 to 1.2e-07 in every row, so the quadratic functional
cannot produce that pattern.

### `coloured_diffusion`

The docstring paragraph before the move:

⚠ THE `CY/2` IS THE SAME ONE-SIDED-TO-TWO-SIDED CONVERSION THE
REST OF THIS CLASS USES, and it is shared rather than repeated so
the pair cannot drift the way `diffusion_constant` and
`covariance` once did over exactly that factor.

### `coloured_diffusion_resolved`

2026-09-26: a source whose level follows the orbit takes
`_coloured_diffusion_modulated`: the harmonics of the PRODUCT `v_1^T G`
per coloured component, the white parts as Demir's `c` (flat in `f`).  So
`phase_psd` takes it too.


The docstring paragraph before the move:

⚠ THIS IS THE OBJECT `c + Gamma(f)` STOOD IN FOR, and the stand-in
is wrong in two ways that the fixture could not show: `c` reads
`CY` at ONE frequency (`2 pi / T`) as if it held at every harmonic,
and `Gamma` is exactly the `l = 0` term of this sum, so `c + Gamma`
counts `l = 0` twice.  Neither was visible on van der Pol, whose
PPV at the tank node averages to zero (`|V_0|/|V_1| = 5e-13`: the
inductor shorts the node at DC, so no core can bias it) -- the
fixture shared the claim's assumption, failure shape 0b.

### `phase_psd`

The docstring paragraphs before the move:

a white source it is `c` exactly; the earlier `c + Gamma(f)` form
counted the `l = 0` term twice and is retired.

⚠ THE CONVENTION IS PINNED BY `oscillator_spectrum`, NOT ARGUED.
`lorentzian`'s far skirt is `i^2 f_0^2 c / f^2` exactly, and that
object was gated by power conservation to 1.000000.  So this
expression is the same quantity its tail already reports, with the
coloured term added -- no second convention is introduced, which
is the only reason a `S_phi` is shipped here at all after a
one-sided/two-sided error cost this class a factor of two.

The comment on the corner:

⚠ THE CORNER IS THE WHITE LORENTZIAN'S, read at the carrier as it
always was.  For a coloured source `f_h = pi i^2 f0^2 c` is not a
lineshape parameter at all -- there is no Lorentzian -- and
taking the folded value nearest the carrier instead put a 1/f
source's corner ABOVE the offsets, in front of the power bound
below, which is the floor that actually binds for colour.

The comment on the power bound:

⚠ POWER CONSERVATION AS A SECOND, INDEPENDENT FLOOR -- and for a
COLOURED source it is the binding one, by orders.  The
normalised lineshape integrates to 1, and the integral over one
box of width `df` on each side is a lower bound on it, so

    2 df S_phi(df) <= 1

is NECESSARY for the linearised skirt to be consistent with
unit power.  Vanassche, Gielen & Sansen (2003) derive the same
statement for a 1/f input and reduce it to
`df_c >= eps f0 sqrt(2 f_1f)`; the form here needs no
assumption about the source's colour, and REPRODUCES their
worked example exactly -- 100.000 Hz against their ">= 100 Hz"
at `eps^2 = 1e-19`, `f0 = 1 GHz`, `f_1f = 50 kHz`.

⚠ THE LORENTZIAN CORNER ABOVE DOES NOT CATCH THIS.  It is built
from `c` alone, so it knows nothing about a `Gamma(f)` that
grows as the offset falls.  MEASURED on this class's own
flicker fixture: the power bound bites at 2.5e-06 Hz while the
Lorentzian corner sits at 8.2e-09 Hz -- 306x too permissive,
and the swept spectrum was carrying 3.10x unit power at the
bottom of the range before this check existed.

⚠ IT IS A LOWER BOUND ON THE BREAKDOWN, NOT THE BREAKDOWN.
Passing it is not a guarantee: on Vanassche's own example the
observed flattening sits at ~300 Hz, 3x the bound.  So this
refuses what is definitely invalid and admits a band that is
already suspect -- deliberately, because refusing at 3x would
be fitting a threshold to one example.
⚠ AND THE DERIVATION HAS A PRECONDITION THE BOUND DOES NOT
STATE, so it is checked rather than assumed.  The box argument
is `2 df S(df) <= integral_{-df}^{+df} S <= 1`, and the FIRST
inequality needs `S(f) >= S(df)` for every `|f| <= df` -- the
spectrum must not dip below its edge value anywhere further in.
True of a monotone skirt; TRUE of the flattened near-carrier
shape; true even with a spur, which ADDS power inside rather
than creating a dip.

⚠ FALSE FOR A LOCKED PLL, whose phase-noise transfer function
is HIGH-PASS: the spectrum is SUPPRESSED at DC and rises to the
free-running level beyond the loop bandwidth, so it dips below
its edge value everywhere inside.  The bound is not thereby
shown to be violated there -- total power is still 1 -- it is
NO LONGER DERIVED, and a floor that is not derived cannot be
used as one.  Unreachable today because this method refuses a
driven circuit, and squarely in the way of the driven-oscillator
work, which is why it is a check and not a comment.

### `band_spread`

The docstring paragraph before the move:

⚠⚠ WHY THIS EXISTS.  Far above the AM corner both AM and PM fall as
`1/r²`, so `S·r²` is flat and a band mean IS a point value — that is
the case every gate in this tree was written on, and it makes the
distinction invisible.  It is NOT general: MEASURED 2026-09-09, a
source behind a slow RC node has an in-band spectrum that is not
`1/r²` at all (its `k = 0` term is filtered at the RC corner while the
`k >= 1` terms are not, and their mix moves across the band), and its
slow/core ratio swings 1.16 -> 0.87 across `0.08 … 0.15 f0` — so a
band mean and a point value differ by ~4 % there, which is larger
than most of the agreements this file asserts.  A comparison that
takes a band mean on one side and a point value on the other is then
measuring the convention, not the physics.

### `oscillator_spectrum`

The docstring before the move (after its first line):

⚠⚠ THIS DOES NOT GO THROUGH `pnoise`'s SIDEBAND FOLD, AND IT CANNOT.
The fold is a FREQUENCY-CONVERSION computation, and for a driven circuit
-- a mixer, the diode-mixer fold case -- that is complete.  For an
AUTONOMOUS oscillator it is structurally incomplete, and what it omits is
exactly the near-carrier phase-noise skirt this method returns.  Rizzoli,
Mastri & Masotti (IEEE MTT 42-807, 1994) state it directly: frequency
conversion alone is insufficient for autonomous circuits, because the
noise-induced FREQUENCY MODULATION OF THE CARRIER at low offsets is not a
frequency-conversion effect (verified at the source 2026-09-08: p. 807,
Introduction, verbatim "frequency-conversion techniques alone are not
sufficient to solve the noise analysis problem for general autonomous
circuits (oscillators). An important further aspect that must be taken
into account is the noise-induced frequency modulation of the carrier
taking place at low frequency offsets, which is not a
frequency-conversion effect").  Their Section III (p. 810) NAMES the
two stacks: CONVERSION noise, power exchanged among the sidebands of
the unperturbed steady state, "invariably raises as 1/f for f -> 0,
which is not consistent with the measured behavior"; MODULATION noise,
"a jitter of the oscillatory steady state", proportional to noise power
over f^2 so the PSD "raises as 1/f^3 for f -> 0 in agreement with the
measured performance"; and the two DECOUPLE exactly at the steady state
(M_BH = M_HB = 0).  They also say the two are "usually nearly equal" in
an INTERMEDIATE offset band, "so that (20) and (21) are
interchangeable" -- a cross-stack agreement test this tree does not yet
have (recorded in the roadmap, not built).  ⚠ Their construction is
harmonic balance; what transfers is the classification, the two slopes
and the interchangeability, none of which need HB.  Diagnostic value:
a FLAT PSD near the carrier is neither slope -- it is the Phi(T) - I
singularity, not the conversion model being the wrong physics.

So the two stacks -- the Floquet/PPV one (`ppv`, `diffusion_constant`,
this method) and the sideband fold (`pnoise`) -- ARE NOT TWO
IMPLEMENTATIONS OF ONE QUANTITY, and unifying them is not a
simplification waiting to be made.  ⚠ THE HAZARD IS THAT THE WRONG ONE
STILL RETURNS A NUMBER: deriving oscillator phase noise from the fold
alone would produce a spectrum -- the conversion terms are real and
non-zero -- just one missing the dominant contribution near the carrier.
A plausible wrong answer, which is the failure shape this whole area
keeps generating.  That is the completeness argument for the split; the
efficiency argument (Floquet is cheaper) is the weaker one and was for a
long time the only one written down.

Returns `(S_v, L_dBc)`.  ⚠ `S_v` is the Lorentzian lineshape scaled by
`|X_1|^2 = A^2/4`, the carrier PHASOR's square -- which is HALF the
carrier power `A^2/2` a one-sided PSD carries, so `S_v` is exactly
0.5000x a one-sided PSD of the output voltage (measured against
a reference simulator at every offset over four decades, 2026-09-05).  `L_dBc`
is unaffected, `|X_1|^2` dividing out of the ratio; the absolute
V^2/Hz matters to anyone integrating `S_v` to a power, and the
scale is kept rather than doubled because it is a return value
that callers may already divide by `|X_1|^2` themselves.  `S_v`
was documented as the one-sided PSD of the output
voltage; `L_dBc` is that normalised to the harmonic's own power,
in dBc/Hz.

⚠ NO SWEEP AND NO PER-FREQUENCY SOLVE.  Once the PSS waveform's
Fourier coefficients and the scalar `c` are known, "we have an
analytical expression that gives us the spectrum at any frequency.
The computation of the spectrum is not performed separately for
every frequency of interest."  Which also means it never meets the
near-carrier singularity that a swept small-signal computation
would, and never meets the 1/f sweep-grid trap — there is no sweep
to place a point on.

⚠⚠ SCOPE: A SOURCE BEHIND A SLOW NODE (A2, resolved 2026-09-08).
The Lorentzian uses the DC PPV, so for a noise source that reaches
the core through a slow path (RC leg, tau >> T) it holds only
BELOW the source's corner `T/(2 pi tau)`; above it the true skirt
is this one scaled by the PPV-harmonic-weighted filter
`sum_k |G_k|^2 F_k(f) / sum_k |G_k|^2 F_k(0)` (G_k the PPV entry's
Fourier coefficients at the source node, F_k the path's transfer at
k f0 + f), which is 1/1000 at 0.1 f0 on a one-RC-leg fixture with
an asymmetric core AND tank loss (both needed for G_0 != 0; an
ideal tank inductor shorts DC).  `c` is still right (the filter
removes only high-frequency content); `pnoise` computes the true
value at any offset, and a Monte Carlo of `c` cannot see it.
✅ SINCE 2026-09-14 THIS METHOD DOES TOO, by default:
`frequency_aware=True` replaces `c` by `c(f)` from the
frequency-aware PPV (`frequency_aware_diffusion`), which matches
pnoise's PM content to 0.4 % at 1e-3 and 1e-2 f0 on that fixture
(the DC-PPV Lorentzian: 0.73x and 0.027x), and to <= 2 % above
f_amp on an orbit with AM-to-PM coupling.  `frequency_aware=False`
is the closed form, one `c` for every offset, and costs no solve.  Measured against `pnoise` to four digits
through the corner (test ..._behind_a_slow_node_...).

⚠ AND IT IS THE ONLY ROUTE THAT IS VALID BELOW THE CORNER.  A
small-signal analysis cannot produce `L(f)` there however well
conditioned it is: the excess phase is a Wiener process, its
spectrum has a singularity at the origin and no physical meaning,
and the finite value `L` attains comes from the NONLINEAR
phase-to-voltage map — which is what this closed form carries.
Reporting `S_phi` near the carrier instead is the mistake that
object invites.

The comment on the missing carrier:

⚠⚠ NO CARRIER, NO LINE -- AND THE ANSWER WOULD BE A PLAUSIBLE ZERO.
This is a LINE-SHAPE model: it broadens the carrier's own harmonic.
Where the output has no component at `harmonic` (a half-wave
symmetric orbit's even harmonics, or `harmonic = 0`, which the
Lorentzian returns as zeros by construction) there is nothing to
broaden, and the true density is BROADBAND noise this method does
not represent.  Measured (2026-09-14, Monte Carlo, which agreed
with pnoise to 1 % at every harmonic): van der Pol C=4 Q=8, 2 f0 +
10 f_amp -- this returned ~0 against 5.6e-6 V^2/Hz.  Refused, as
`am_pm` refuses the same case.  ⚠ NOT caught, and recorded instead:
away from the fundamental the model also misses where a line DOES
exist (asymmetric orbit, 2 f0: 0.40 of the Monte Carlo) -- use
`pnoise` away from the fundamental.

### `frequency_aware_diffusion`

The docstring paragraph before the move:

⚠⚠ WHY IT EXISTS (2026-09-14).  The Lorentzian from `c` uses the DC
PPV at every offset: a noise current is assumed to move the phase
instantly.  Wherever part of that response goes THROUGH a slow mode --
the amplitude mode on an orbit with AM-to-PM coupling, or a slow node
in the source's path -- it is filtered above that mode's corner, and
the DC-PPV Lorentzian over-states.  Measured against pnoise's PM
content (itself Monte-Carlo-confirmed on the first case):

    van der Pol C=4 Q=8, a=0.30     0.3 / 1 / 3 / 10 f_amp
      c(f)/c                         0.939 0.649 0.371 0.308
      pm / 4 S_v(DC PPV)             0.942 0.647 0.366 0.302
      pm / 4 S_v(c(f))               1.003 0.998 0.986 0.980
    A2 slow node (tau/T=100)         1e-3 / 1e-2 f0
      pm / 4 S_v(DC PPV)             0.729 0.027
      pm / 4 S_v(c(f))               1.004 1.004
    symmetric control (a=0)          c(f)/c within 1e-3

### `_warn_above_amplitude_pole`

The docstring before the move (from its first paragraph's last sentence):

adds to the total, so this method UNDER-reports.  Relayed measurement
against a commercial simulator's total noise, as excess over the
phase-only prediction:

    offset     lam2 = 0.90            lam2 = 0.99
               (f_amp 26.7 kHz)       (f_amp 2.55 kHz)
    100 Hz     -0.00 dB               -0.01 dB
    1 kHz      -0.00 dB               -0.54 dB
    10 kHz     -0.50 dB               -2.90 dB
    100 kHz    -3.11 dB               -3.27 dB

⚠⚠ AND THE VALID REGION SHRINKS AS `1/Q`, which makes this section 0
again rather than a detail.  With `f_amp = -ln(lam2)/(2 pi T)` and
`Q = -1/ln(lam2)`,

    f_amp = f0 / (2 pi Q)

-- verified both ways at 26671.9 / 2544.2 / 253.3 Hz for
`lam2 = 0.90 / 0.99 / 0.999`.  So the better the oscillator, the
narrower the band in which its phase-only spectrum is the whole
answer; at `lam2 = 0.999` it has collapsed below ~253 Hz.

⚠ THIS IS THE OPPOSITE SIGN FROM THE ERROR `PSS.ppv` ALREADY WARNS
ABOUT.  That one says the instantaneous phase equation misses slow
nodes which FILTER device noise, so phase noise is OVER-estimated.
This one is a second, independent mechanism in which the phase-only
answer is UNDER-estimated.  Both are live and they are not the same
effect.

### `_output_waveform_row`

The docstring before the move (after its first line):

⚠ THE REST OF THIS CLASS TAKES A DIRECTION VECTOR AND THIS PAIR
TOOK AN INTEGER, which is not a style difference -- it meant
`am_pm` and `carrier_phasor` could not express a DIFFERENTIAL
output at all.  `pnoise`, `adjoint_transfer_row` and
`adjoint_sideband_row` all accept `d`; these did `int(output)`.
For an oscillator the output of interest is very often
differential, and for the coordinate-invariance property an AM/PM
split has to have (Kaertner 1990 section 3.2) a
reference-independent observable is the whole point.

An integer is still accepted, so callers that name a node keep
working; an array is contracted against the full waveform with the
reference row reinserted.

### `am_pm`

The docstring paragraph before the move:

⚠ THE ABSOLUTE MAGNITUDE ON AN OSCILLATOR IS SMALL FOR A REASON
(established 2026-09-08, the three-leg chain).  These are the p = 0
band of `am_pm_noise`: a source at BASEBAND `freq` reaching the
carrier sideband.  A baseband current moves the PHASE through the
PPV's DC coefficient (Hajimiri-Lee's c_0), and a half-wave-symmetric
orbit -- odd nonlinearity, `u(t + T/2) = -u(t)` -- has none, so on
such a fixture the rows measure a symmetry zero (6e-9 .. 1e-13,
proportional to 1/freq and to mu), the same zero the coloured
up-conversion gate records for Gamma.  Breaking the symmetry
(`_lc_osc(a)`) lifts |m_pm| at 1e-3 f0 from 1.2e-8 to 46.6 (a =
0.05) and 231 (a = 0.25) -- linear in `a`.  The DIRECT rows (source
at f0 + freq, sideband 0) agree with `pnoise` at every offset, and
the split lands on the externally certified Lorentzian.  So do not
read a small `am_pm` on a symmetric oscillator as a defect: it is
the 1/f^3 up-conversion coefficient, and it is zero there.

### `am_pm_noise`

The docstring paragraphs before the move:

⚠ THE AUTONOMOUS CAVEAT IS RETIRED (three-leg chain, 2026-09-08).  On a
free-running oscillator this split sits on the SAME absolute scale as
`pnoise` (identity to 1e-12 / 1e-16) and as the externally certified
`oscillator_spectrum` (`S_pm = 4 S_v` at every offset: the PM content
of the pair IS the Lorentzian, 2 S_v per sideband), with `S_am` rising
from ~0 below the AM corner `f0/(2 pi Q_lambda)` to `S_pm` above it
(⚠ this line said `f0/(4 pi Q)` until 2026-09-09 -- a factor of two
the docs session caught against this very function: the ratio is an
exact Lorentzian `u^2/(u_c^2 + u^2)` in `u = offset/f0` with
`u_c = 1/(2 pi Q_lambda)`, `Q_lambda = -1/ln|lambda_2|`, half-power
0.5007 there and 0.20 at the old corner, at Q = 8 and 16, 240 and
480 points; the old formula OVERSTATED the AM content at every
offset, 2.5x at its own corner) -- so the
pair total is 4 S_v there and 8 S_v far out.  The "~1e-12 rows" were
`am_pm`'s p = 0 band on a half-wave-symmetric fixture: a symmetry
zero, see `am_pm`.  Oscillator magnitudes from this are trustworthy.

⚠ WHAT A MEASUREMENT MUST BE TO BE COMPARED WITH `S_pm` (2026-09-09,
after a 16-seed Monte Carlo campaign and a deterministic forward
tone route; roadmap "Item 3, answered"): `S_pm` is PM BY QUADRATURE
OF THE FUNDAMENTAL'S SIDEBANDS.  A "phase" read by a one-period
demodulation of the fundamental leaks the other harmonics'
sidebands through its boxcar (sinc(pi(1 - r)) ~ 0.1 in amplitude
for the second harmonic's, which is ~ the orbit's asymmetry), and
a phase read from zero crossings converts EVERY harmonic's
sidebands; both drifted 6 % against this quantity as the asymmetry
of a harmonic-rich orbit was swept while the forward LPTV response
agreed with it to 1 %.  So compare `S_pm` with the fundamental's
sideband PM (a spectrum analyser's sidebands around f0, or the
forward-tone gate `test_pnoise_oscillator_pm_matches_a_forward_tone_
transient_with_no_adjoint`), never with a demodulated or
crossing-time phase, and BAND WITH BAND: a source behind a slow node
has an in-band spectrum that is not 1/r^2 (its slow/core ratio
swings 1.16 -> 0.87 across 0.08-0.15 f0), so a band mean and a
point value differ by ~4 % there.

The comment on the carrier's frame:

⚠⚠ THE SPLIT IS TAKEN IN THE CARRIER'S FRAME, NOT THE TIME ORIGIN'S
(defect reported by a peer session and reproduced 2026-09-14).  AM
is the envelope component ALONG the carrier phasor, so `a + conj(b)`
is right only for a cosine-phased carrier.  Until this rotation the
answer depended on where t = 0 sat: a driven diode gave am/pm =
0.305 / 3.28 / 0.651 at drive phase 0 / 90 / 37 degrees, and
3.316 at all three once rotated, `S_am + S_pm` unchanged to 1e-15
(the identity cannot see it -- `|a_r|`, `|b_r|` are `|a|`, `|b|`).
`am_pm` never had it: it divides by the COMPLEX carrier phasor.
With no carrier at this harmonic the phase is undefined and the
split is left unrotated, as `am_pm` refuses the same case.

### `_deflated_solve`

The docstring paragraph before the move:

MEASURED on van der Pol, offsets from 0.3 down to 1e-9 of `f0`:

    offset/f0    sigma_min(plain)   sigma_min(bordered)
    3e-01           5.68e-01            1.17e-01
    1e-03           2.61e-03            2.04e-01
    1e-06           2.61e-06            2.04e-01
    1e-09           2.61e-09            2.04e-01

The plain operator tracks the offset over nine decades; the
bordered one is FLAT.  The two solutions agree to 5.7e-12 where the
plain solve is still trustworthy, and their disagreement grows as
`1/df` -- that is the PLAIN solve losing digits, not this one.

The comment on the border vectors:

⚠ THE BORDER VECTORS ARE NORMALISED (2026-09-22): the recovery
`y = w + s col / (1 - alpha)` is scale-free in exact arithmetic,
but GMRES sees the bordered MATRIX, and with the tangent in V/s
(1e6 on a comparator oscillator) against the PPV in s/V (1e-6)
its condition number was 1.8e9 where `I - alpha M` alone was
2.4 -- the forward solve returned with a 5.5e-4 residual and
the transposed one did not converge at all.  Unit vectors put
the bordering at the operator's own conditioning; `s` absorbs
the scale.

The comment on the staged solve:

⚠ ON A STAGED SOLVE THE POLE IS THE TOTAL MAP'S (2026-09-22,
events phase B): `u`, `v` are the null vectors of `I - M_tot`,
`M_tot = M + P_theta dtheta/dx_0`, and the fixed-grid `M` has no
unit multiplier at all (|lambda - 1| = 0.99 on the comparator
oscillator) -- bordered with the total map's vectors it read
the sideband response 0.3-400x off the exact one.  The
operator here is the total map's.

The comment on the refinement:

⚠ REFINED ON THE PLAIN OPERATOR WHERE THAT IS WELL CONDITIONED
(2026-09-22): the recovery assumes `u`, `v` are EXACT null
vectors of `I - M`.  On a staged solve the discrete total map's
unit multiplier is displaced by O(h) (8e-4 at 200 points on the
comparator oscillator, 4e-4 at 400; an unstaged radau map sits
at 1e-11), and the recovered `y` then misses the true operator
by that much over `|1 - alpha|`: 5.5e-4 at 0.3 f0, 0.14 at
1.001 f0.  The plain operator is well conditioned there (2.4 at
0.3 f0, 310 at 1.001 f0, 2e3 at 1.00001 f0 -- its pole sits
where the DISCRETE multiplier is, not at alpha = 1), so the
deflated answer seeds a plain correction on its own residual,
kept only if it lowers the residual; below
`DEFLATION_REFINE_MIN` the deflated answer stands.  ⚠ This makes
the forward and adjoint solves the discrete operator's own,
dual-consistent -- and near a harmonic the discrete operator's
answer carries the multiplier's displacement: 15.7 % off the
exact forced response at 1.001 f0 on the 200-point staged
comparator oscillator (0.2 % at 0.3 f0 and 1.7 f0), where the
unrefined recovery, which carries the pole analytically at
alpha = 1, happened to read 0.2 % but is not dual-consistent
(4.5e-4 at 0.3 f0) and has its own O(h/|1 - alpha|) split
error.  The item is the staged map's multiplier (roadmap E8).
On an unstaged oscillator the residual is already at the
tolerance and nothing happens.


## `probe.py` -- `ProbeShooting`

### `ProbeShooting`

The class docstring's warnings before the move:

⚠ IT IS NOT A CONVERGENCE AID, AND THE PAPER SAYS SO ITSELF.  On its own
flagship high-Q Pierce example the authors report *"it is easy to assign a
tentative current to the Ls inductor ... and obtain convergence in a few
iterations (we did this with conventional SH)"*.  What this buys is the
SWEEP: unstable limit cycles, coexisting solutions, and a stability screen.

⚠⚠ ONE TONE GIVES A DESCRIBING-FUNCTION SOLVE, NOT THE ORBIT.  A single
tone forces a SINUSOID, so a non-sinusoidal orbit can only null the probe's
FUNDAMENTAL current.  Measured on van der Pol, the frequency error is
QUADRATIC in harmonic content -- ``df/f = 4.0 THD^2`` to 3% across two
decades -- and the probe sits at the LC resonance at every ``mu`` because
``mu(u - u^3/3)`` is odd and memoryless, so its describing function shifts
no phase.  ``harmonics=K`` forces K tones and nulls K harmonics, which is
harmonic balance with a shooting inner solve: at K=3 the van der Pol error
falls from 5.9e-02 to 3.8e-04.

⚠⚠ PROBE PLACEMENT IS CIRCUIT-SPECIFIC, AND ITS FAILURE IS NOT A SOLVER
FAILURE.  Forcing a node fixes every state the source reaches; any state it
does NOT reach whose DC level is then unconstrained makes the shooting
Jacobian SINGULAR, because a whole family satisfies periodicity.  Measured
across van der Pol's only node with no series resistance: periodicity error
**2.11e-15** -- already a periodic solution -- reported as
``converged = False``.  :meth:`degenerate_placement` names that pairing.

### `_build`

The comment on `vac=0`:

⚠⚠ `vac=0` EXPLICITLY.  `VS.vac` DEFAULTS TO 1, not 0, so every
probe source in the chain would be AC-excited at once -- and in
series they share one current, so the PAC response came back
exactly K times too large (measured ratios 1, 2, 3 at K = 1, 2,
3).  K=1 validated because there was nothing to contaminate it,
which is why the diagonal alone could not catch this.

### `_spectrum`

The comment on the warm start:

⚠ WARM START: the finite-difference columns perturb a parameter by
~1e-5, so the trajectory barely moves and solving each from cold is
waste.  Measured 2.09 s cold against 1.31 s warm -- 1.60x -- with
the answers agreeing to 3e-15.  A pure accelerator: it changes which
iterate the Newton starts from and nothing else, and carries no
assumption about the circuit.

### `even_harmonic_content`

The docstring's warning before the move:

⚠⚠ PRUNING THE EVEN HARMONICS DOES **NOT** HOLD IN GENERAL, and the
failure is silent: dropping a tone that is really there removes both an
unknown and the residual row constraining it, so the solve converges to
the wrong waveform.  Van der Pol is HALF-WAVE SYMMETRIC and its even
content is 7.1e-16; add an even term ``beta u^2`` to the same
nonlinearity and it is not::

    beta    H2/H1       H3/H1       H4/H1
    0.00    7.080e-16   1.168e-01   2.710e-16   <- symmetric
    0.05    2.649e-02   1.162e-01   9.239e-03
    0.20    1.058e-01   1.068e-01   3.600e-02   <- H2 EQUALS H3
    0.50    2.613e-01   5.994e-02   7.610e-02   <- H2 is 4x H3

At ``beta = 0.5`` pruning would discard the LARGEST correction after
the fundamental, and ``beta = 0.05`` already gives 2.6% -- there is no
margin to judge by eye.  So this measures rather than assumes, and a
caller passing `tones=[1, 3, 5]` should check it first.

### `solve_multitone`

The docstring's measurement and note before the move:

Measured on van der Pol (`mu = 1`, autonomous `f = 0.150229`)::

    K   f          df/f         note
    1   0.159134   +5.93e-02
    2   0.159134   +5.93e-02    A_2 = 2.3e-13 -- no change at all
    3   0.150172   -3.77e-04    157x better

⚠ K=2 buys NOTHING here because the even harmonics do not exist on a
half-wave-symmetric circuit -- NOT because two tones cannot help.  See
:meth:`even_harmonic_content` before concluding the same elsewhere.

The comment on the chain rule, its middle lines:

across positionally feeds the Newton a Jacobian for the
WRONG VARIABLES -- it diverged to f = 0.0348 against 0.1502
while `pac_jacobian` itself validated at 1e-04, which is how
a correct derivative and a broken solve coexisted.

### `_pac_response`

The docstring's opening warning and first bullet before the move:

⚠⚠ TWO CONVENTION FACTORS AND ONE INDEXING TRAP, all three pinned
against a circuit whose answer is analytic (a resistor across the
probe, where `dI/dV = 1/R` exactly) rather than against the finite
difference this is meant to replace.  Calibrating against FD would
make the agreement circular and would absorb a genuine sideband-index
error into the fitted constant.

* **The frequency list carries DUPLICATES.**  `0.159155` appears twice
in the returned sweep, one entry near zero and one carrying the
response; `argmin(|fs - target|)` picks whichever comes first and it
was the wrong one, reading 5.9e-21 where the answer is 1e-3.  So the
entry is chosen by LARGEST RESPONSE among those at the target
frequency, not by proximity alone.

The comment on the excitation offset, its first paragraph:

⚠⚠⚠ EXCITE SLIGHTLY OFF THE HARMONIC, WHICH REMOVES THE AMBIGUITY
INSTEAD OF GUESSING IT.  Exciting exactly at `j*f0` sends TWO
sidebands to the same absolute output frequency -- `k = m - j` and
`k = -m - j` -- and `PAC.solve` returns absolute frequencies with the
sideband index folded away, so the two arrive in an order that is not
stable.  Ordering them by magnitude worked AT THE SOLUTION and failed
away from it (validation 1.6e-05 at the solved amplitudes, 1.763 at
the Newton's starting point), which is the kind of heuristic that
passes a gate and then fails in use.

The comment on the direct and image terms, its end:

conjugates.  Taking only the direct one halves the answer
(measured: a uniform 0.5 at every K); taking both without being
able to tell them apart is what forced the magnitude-ordering
heuristic that passed at the solution and failed at the Newton's
start.  With the offset they are identified by FREQUENCY, so the
rule is derived rather than guessed.

### `_pac_response` -- a comment on the superseded magnitude-ordering method

The comment on dead code removed from `_pac_response`:

⚠⚠ THE DUPLICATES ARE A PAIR AND THEY SUBTRACT, NOT ADD.
PAC folds negative output frequencies onto |f| and conjugates
them, so two entries land at each harmonic.  Measured, at the
fundamental:

    resistor       e1 = -3.4e-21   e2 = +1.0e-03
    van der Pol    e1 = -9.900250e-01   e2 = +9.900507e-01

They are nearly EQUAL AND OPPOSITE on a circuit with harmonic
content, so SUMMING them cancels (2.6e-05 against a true
1.98) and taking the LARGEST halves the answer -- which is
exactly the factor 2 that appeared on van der Pol and not on
the resistor, where one member is ~0 and max == sum == diff.
The difference reproduces the finite difference on BOTH
fixtures at once, which is the falsifier a fitted constant
could not have passed.
⚠⚠⚠ ORDER THE PAIR BY MAGNITUDE, NOT BY ARRAY POSITION.
`PAC.solve` returns ABSOLUTE frequencies -- the sideband index
is folded away -- so the two entries at a harmonic arrive in an
order that is NOT stable across `m`.  Measured, exciting
harmonic 1 at K=3:

    m=1   [ +0.75j , 3.94e-04 - 1.041251j ]   larger is 2nd
    m=3   [ -0.75j , -0.25j             ]   larger is 1st

Taking `e[-1] - e[0]` therefore gave m=1 correctly and m=3 with
the RIGHT MAGNITUDE AND THE WRONG SIGN (ratio -0.999985).  The
larger entry is the DIRECT response and the smaller its image,
which is an ordering the array position does not carry.
⚠ Fragile where the two magnitudes are close; the `validate`
gate is what stands behind it.
exactly one entry now -- the offset separated the pair

### `pac_jacobian`

The docstring's status and validation paragraphs before the move:

⚠⚠⚠ **NOT SHIPPED-READY: THE NORMALISATION IS INCOMPLETE, AND THE
DEFAULT VALIDATION CORRECTLY REFUSES.**  The route is confirmed viable
-- PAC returns the right quantity, verified against a resistor where
`dI/dV = 1/R` analytically -- and two of the three discrepancies are
pinned and removed.  A third is not:

    fixture                 PAC / FD after normalisation
    resistor (linear)       -1.0        (sign only)
    van der Pol (K=1)       -0.499999   (sign AND a factor 2)

**A constant that differs between two circuits is not a convention,
it is a missing term**, so the remaining factor is NOT applied by
fitting it -- that would make the "independent" Jacobian a fit to the
finite difference it replaces, hide any sideband-index error inside the
fitted constant, and reproduce exactly the circular-verification
failure this file already records for the PPV.

Until it is derived, `validate=True` raises on real circuits and the
finite-difference Jacobian in :meth:`solve_multitone` remains the
shipped path.  What IS established and reusable:

* the response is present and correct (1/R recovered exactly);
* the 90 degree rotation is OURS -- `VS` builds its AC phasor as
`vac * exp(j*phase)` and the operating-point probe sets `phase = 90`
-- and dividing by the excitation phasor removes it, measured: the
residual ratio is real, not imaginary;
* PAC's frequency sweep contains DUPLICATE entries at the same
frequency, one near zero and one carrying the response, so
`argmin(|fs - target|)` reads 5.9e-21 where the answer is 1e-3.

⚠⚠ `validate=True` CHECKS ONE COLUMN AGAINST THE FINITE DIFFERENCE AND
RAISES ON DISAGREEMENT, and it is on by default deliberately.  A wrong
Jacobian does not announce itself: the Newton still converges, to the
wrong orbit -- the same silent failure that even-harmonic pruning
produces, which this file already has a falsifier for.  The check costs
ONE extra solve, amortised over the whole solve, and it exercises the
conventions in `_pac_response` on the circuit actually in hand rather
than on the resistor they were derived from.

The comment on the single PSS solve:

⚠⚠ ONE PSS SOLVE FOR EVERY COLUMN.  `vac` is read ONLY under
`analysis='ac'` -- it does not enter the transient residual, so it
cannot move the periodic operating point.  Building a fresh circuit
and re-solving the PSS per column therefore recomputed the SAME
orbit K times, which is the whole cost this method exists to avoid:
it made the "cheap" Jacobian K nonlinear solves plus K linear ones,
against the finite difference's 2K.  Solve once, then walk `vac`
across the probes and take K LINEAR PAC solves against that one
operating point.


## `events.py` -- `EventColumns`

### `EventColumns` (class docstring)

The docstring's opening line and last paragraph, as they stood before the
history moved out (2026-09-24):

A staged solve's EVENT COLUMNS and the linear algebra every bordered
consumer does with them (refactor E9 item 1, 2026-09-23).

Before this, the dict was read in 17 places and three derivations were
re-typed around it: the TOTAL map `M + P_end dth` at four sites, the
bordered adjoint's Schur elimination in two verbatim copies, the
`-zeta_k W_k` injection at four.  They are methods here; the numbers
are unchanged (the refactor's gate compares every consumer's output
before and after).


## `diagnostics.py` -- module level

### `noise_enters_constraints`

The docstring before the move:

`(bad, residual)` — Winkler's index-1 SDAE precondition, `im B ⊆ im C`.

Winkler (JCAM 163:435–463, 2004) Definition 2: an SDAE is **index 1**
when "the noise sources do not appear in the constraints", `im G ⊆ im A`
— in circuit terms, the noise input must lie in the image of the
capacitance matrix.  Otherwise it is an SDAE **WITH DIRECT NOISE**, which
is outside the class that theory covers.

⚠ `CY = B B^T`, so `im B = im CY` and no factorisation is needed: project
`CY`'s columns onto `im C` and look at the residual.

⚠⚠ WHAT IT MEANS WHEN IT FAILS, MEASURED (roadmap §0j).  White noise
applied to a variable fixed by a CONSTRAINT rather than by an integrator
is filtered by nothing, so **that node has no finite variance** — not a
missing term and not a discretisation artefact.  A series tank-loss
resistor with no capacitance at its node fails this; adding a parasitic
capacitor moves the circuit back into the class.

⚠ The PHASE is a different question and is NOT covered by this: `c` on
the failing fixture still agrees with an equivalent parallel-loss circuit
to 0.999973 and with the Lyapunov route.  So this bears on the
COVARIANCE, not on `diffusion_constant`.

### `topological_index`

The docstring before the move:

`(index, info)` — the DAE index from the netlist, WITHIN A STATED CLASS.

⚠ A second source with the same class boundary (docs session, 2026-09-09,
from the rendered pages of Lamour, März & Tischendorf): Lemma 3.45 gives
the same two criteria as RANK conditions on incidence matrices --
`[A_C A_R A_V]` full row rank iff no L-I cutset, `Q_C^T A_V` full column
rank iff no C-V loop -- under the same hypothesis, "let all current and
voltage sources be independent"; Theorem 3.47 adds the index-0 case
(a capacitive path from every node to datum AND no voltage sources) and
⚠⚠⚠ AND IT READS AN INDEX; IT DOES NOT CERTIFY THAT ONE EXISTS.  No
frozen-`t` / structural test can, and that is a theorem rather than a
caution.  Brenan, Campbell & Petzold's counterexample (quoted in Estevez
Schwarz, Lamour & März, "The common ground of DAE approaches"):

    E(t) = [[-t, t^2], [-1, t]],  F = I,  t in [-1, 1]

has `det(lam E(t) + F(t)) = 1` exactly -- EVERY LOCAL PENCIL IS REGULAR --
and yet `x(t) = gamma(t) [t, 1]^T` solves the homogeneous DAE for ARBITRARY
smooth `gamma`.  An infinite-dimensional solution family, so no index is
meaningful at all.  The failure is NON-LOCAL: the pair is pre-regular and
the REDUCED pair has `im[E_1 F_1] = {0}`, which lives in the reduction
sequence and not in any quantity computable at fixed `t`.  BCP's own text
says this regularity notion "does not imply solvability".

So: use this to READ an index on a circuit already known to be solvable,
never to ESTABLISH solvability.  A simulator that reports "index 1"
pointwise has not shown the problem is well posed.  (Relayed from a source
reading; the counterexample was reproduced by that session to 4.4e-16 for
three unrelated `gamma`, not by this one.)

⚠⚠ IT USED TO BE FLOORED AT 1 BY CONSTRUCTION, and that is fixed as of
2026-09-11.  Estevez Schwarz & Tischendorf's criterion is "index 2 IF AND
ONLY IF the network contains a C-V loop or an L-I cutset, OTHERWISE 1", so
it could not return 0 and answered 1 for an implicit ODE -- SILENTLY, while
the line below already recorded Theorem 3.47's index-0 case as something
the theory "adds".  Found because a numerical probe read index 0 for a van
der Pol and was CLAMPED to agree with this function, which was agreeing for
the wrong reason.  Theorem 3.47's case is now implemented: `index` can be
0, and `info['cap_path_to_datum']` says whether that test was reached and
what it found.  ⚠ The rank condition is the independent cross-check --
index 0 IFF the reduced `C` is NONSINGULAR, an implicit ODE -- and it is
gated as such.

closed-form projectors (3.61)/(3.62).  ⚠ THIS LINE SAID "not implemented"
UNTIL 2026-09-10 AND WAS STALE BY A DAY: the rank form IS implemented, as
`test_the_topological_index_agrees_with_an_incidence_RANK_criterion`,
which cross-checks this graph traversal on four topologies and makes each
condition FAIL on the one it names (an all-pass comparison would prove
nothing).

⚠ A THIRD SOURCE, corroborating the class boundary below rather than
extending it (docs session, 2026-09-10, Lamour, März & Tischendorf
Remark 3.49): Theorem 3.48's constant-projector structure "remains valid
also for CONTROLLED current and voltage sources IF they do not belong to
C-V loops or L-I cutsets and their controlling voltages and currents do
not belong to C-V loops or L-I cutsets".  That is the same exclusion the
Estevez Schwarz & Tischendorf note below already states, reached from a
different book -- it does not widen or narrow what this function claims.

⚠⚠ NOT "FROM THE NETLIST ALONE", WHICH AN EARLIER VERSION OF THIS LINE
CLAIMED AND WHICH IS FALSE FOR CONTROLLED SOURCES.  Estevez Schwarz &
Tischendorf close the paper by giving up BOTH halves of the criterion for
that case: "if arbitrary controlling elements for the controlled sources
are considered then THE INDEX OF THE NETWORK EQUATIONS MAY DEPEND ON THE
PARAMETERS", and "if controlled sources are allowed to form a part of L-I
cutsets or C-V loops then IT IS POSSIBLE TO BE CONFRONTED WITH HIGHER
INDEX (> 2) PROBLEMS".

So `provisional` is NOT a lower-confidence index-2 verdict — **it is not
an index-2 verdict at all.**  Two independent failures at once: the index
is no longer bounded by 2, and it is no longer a function of the topology
at all, because it can turn on element VALUES.

Estevez Schwarz & Tischendorf (IJCTA 28(2):131–162, 2000): for a
nonlinear time-independent network **without controlled sources**, and
assuming positive-definite element Jacobians,

    the index is 2 IF AND ONLY IF the network contains a C-V loop or an
    L-I cutset; otherwise it is 1.

⚠⚠ C-ONLY LOOPS ARE *NOT* COUNTED HERE, AGAINST THE RELAYED QUOTE, AND
THE MEASUREMENT IS WHY.  The paper is quoted as saying "C-only loops have
to be added to the class of C-V loops since the currents through C-only
loops belong to the network variables whereas these currents are excluded
in MNA formulations", and a first version of this function counted them.
**Measured against a direct computation on our own MNA matrices — `C`'s
null basis `N`, then the rank of `N^T G N` — all three C-only topologies
come out INDEX 1:**

    capacitor ring touching ground     topological 2, measured 1
    floating capacitor triangle        topological 2, measured 1
    the same with every node grounded  topological 2, measured 1

and the arithmetic is checkable by hand.  For the grounded ring
`C = [[c1+c3, -c1], [-c1, c1+c2]]` has determinant
`c1 c2 + c1 c3 + c2 c3 != 0` — NOT SINGULAR, so the system is not even a
DAE.  For the floating triangle `C` IS singular (it is the triangle's
Laplacian, null vector all-ones) but `N^T G N = (1/R)/3 != 0`, so the
constraint is uniquely solvable and the index is 1.

**A C-only loop makes `C` singular WITHOUT making the index 2; index 2
needs a VOLTAGE SOURCE fixing the loop.**  ✅ RESOLVED AT THE SOURCE
(docs session, 2026-09-08; Estévez Schwarz & Tischendorf, IJCTA 28(2)
2000, on disk): the quote is FAITHFUL and describes Chua & Lin's
variable set, not MNA.  Their MNA theorem, Thm 4.1 p.141, has no C-only
clause -- "the conventional MNA leads to an index-1 DAE if and only if
the network contains neither L-I cutsets nor C-V loops.  Otherwise ...
index-2" -- which is exactly what this function implements.  The
C-only sentence is Remark 4 p.143, comparing with Table 10-3-1 of Chua
& Lin (Reference [10], the normal-tree / state-variable formulation),
and states its own reason: "in this case, C-only loops have to be
added to the class of C-V loops SINCE THE CURRENTS THROUGH C-ONLY LOOPS
BELONG TO THE NETWORK VARIABLES WHEREAS THESE CURRENTS ARE EXCLUDED IN
MNA FORMULATIONS."  So the measurement above AGREES with the theorem;
the loop test requires a voltage source, and nothing is split.  Thm 4.2,
immediately below, extends the same conclusions to the CHARGE-ORIENTED
MNA -- this tree's formulation, `d/dt q(x) + i(x) + u(t) = 0` -- so the
tree is covered by name, under 4.1's hypotheses (positive-definite C, L
and conductance matrices, the controlled-source conditions of their
Tables I-VI), which this docstring already carries.

⚠⚠ THIS IS A DIAGNOSTIC, NOT A REFUSAL, AND THE DIFFERENCE IS MEASURED.
Roadmap C4 closed index-2 detect-and-refuse because `index > 1` is NOT
PREDICTIVE: all three integrators converge on an L-I cutset and Gear-2
fails on 2 of 4 index-2 topologies.  Deciding the index was never the
obstacle — knowing it exactly and still not knowing which method to use
is the actual state of things.  What this buys is a BETTER MESSAGE when
something does fail.

⚠ WHICH IS WHY IT LOCALISES.  The authors' stated design goal is
"topological criteria that can be checked very fast ... based on LOCAL
assumptions, i.e. we want to provide the opportunity to LOCALIZE
critical element modellings", for networks of ~1e7 elements where "it is
often difficult to find the circuit configurations that lead to
numerical difficulties".  So `info['loop']` and `info['cutset']` name the
ELEMENTS, not just the verdict.

⚠ **V-ONLY LOOPS AND I-ONLY CUTSETS ARE REPORTED SEPARATELY, AND THEY ARE
NOT INDEX 2.**  A loop of voltage sources over-determines KVL and a cutset
of current sources over-determines KCL, so the MNA system is
STRUCTURALLY SINGULAR — it has no solution at all, barring an exact
cancellation.  The index criterion presumes a well-posed network.  Calling
those "index 2" would send a reader hunting a solver problem instead of a
netlist error, so they come back as `info['v_loop']`, `info['i_cutset']`
and `info['ill_posed']` — and **`index` is then `None`**, because the DAE
index presumes a solvable system and there is no honest value to give.

`info['unclassified']` lists elements outside the covered class —
controlled sources above all, which the theorem excludes.  When it is
non-empty the verdict is PROVISIONAL and `info['provisional']` is True;
the criterion is reported rather than withheld, because a named
assumption beats a silent refusal.

The comment on the C-V loop search:

⚠ CAPACITORS FIRST, THEN SOURCES ONE AT A TIME, so the closing
element is ALWAYS a voltage source.  A first version searched C and V
together and discarded any loop that turned out to have no `V` in it
-- which is wrong on a netlist carrying BOTH a C-only loop and a C-V
loop, because union-find returns only the FIRST closing edge and the
C-only one can close first, hiding the real one.

The comment on V-only loops and I-only cutsets:

⚠ ANDREAS ASKED FOR THIS AND IT IS A DIFFERENT CATEGORY.  A loop of
voltage sources over-determines KVL and a cutset of current sources
over-determines KCL: the MNA system is STRUCTURALLY SINGULAR and has
no solution at all (barring an exact cancellation), rather than having
a higher index.  The index criterion presumes a well-posed network, so
these are reported separately -- calling them "index 2" would send a
reader looking for a solver problem instead of a netlist error.
⚠ the V loop is localised the SAME way as the C-V loop -- naming only
the closing source would point at one of three parallel sources and
leave the reader to find the rest, which is the opposite of the point.

The comment on the ill-posed return:

⚠⚠ AN ILL-POSED NETLIST HAS NO INDEX, AND REPORTING ONE IS WORSE THAN
REPORTING NOTHING.  The DAE index presumes a solvable system; a V loop
or an I cutset makes MNA structurally singular, so `index` comes back
`None`.  A first version returned 2 here and ALSO mislabelled the
offending set -- a loop of three voltage sources was reported as a
"C-V loop" containing no capacitor, because the C-V search unions
capacitors then sources and a pure-V loop closes on a source.  Both
symptoms point a reader at the solver when the netlist is the error,
which is precisely what this split exists to prevent.

The comment on the index-0 rung:

⚠⚠ THE INDEX-0 RUNG, Theorem 3.47 (added 2026-09-11).  Estevez Schwarz &
Tischendorf's criterion is "2 iff a C-V loop or an L-I cutset, otherwise
1" and is FLOORED AT 1 BY CONSTRUCTION, so this function used to answer
1 for an implicit ODE -- SILENTLY, and its own docstring had recorded
the omission all along.  MEASURED on a van der Pol (C, L and a BSource,
no voltage source): its reduced `C` is NONSINGULAR, rank 2 of 2 with
sigma_min 1.0, which is index 0.

The condition is "a capacitive path from every node to datum AND no
voltage sources": with no voltage sources the only branch-current
unknowns are inductive, and each carries `L di/dt` in `q`, so every row
of `C` has a reactive entry and `C` is nonsingular.

⚠ An INDUCTOR does not spoil it, which is the case a reading of the
theorem's wording alone might get wrong -- the flux term makes that row
differential, not algebraic.  The van der Pol is exactly that shape and
the rank cross-check agrees.

The comment on the absence test:

⚠⚠ A NEGATIVE CLAIM CANNOT BE MADE PROVISIONALLY, and that is why
`not unclassified` is a condition of this rung and not of the others.
Theorem 3.47's index-0 test asserts an ABSENCE -- a capacitive path from
every node to datum AND NO VOLTAGE SOURCES.  An element outside the
covered class is classified `'?'`, so `kinds[nm] == 'V'` is FALSE for it
and the absence test passes VACUOUSLY.  A `VCVS` is exactly that: a
voltage source this classifier does not recognise.
⚠ MEASURED 2026-09-11, and this rung shipped WITH the defect that
morning: on the documented P1 fixture (a VCVS of gain `g` inside a C-V
loop, index 2 off `g* = 1 + C2/C1` and index 3 on it) this returned
INDEX 0 at every gain, where the docstring above records that every
provisional fixture returned 1.  `algebraic_conditioning` says the
algebraic block is SINGULAR at every gain, i.e. index >= 2, and it is
right.
The other criteria survive `unclassified` because they assert
PRESENCE -- "I found a C-V loop" stands whatever else is in the netlist,
and `provisional` then says only that there may be MORE.  An absence
cannot be established from a partial reading at all.

⚠⚠ BUT "BLOCK ON ANY UNCLASSIFIED ELEMENT" IS TOO STRICT, and the van
der Pol fixture is the proof: its nonlinear conductance is unclassified
and it is GENUINELY INDEX 0 (recorded 2026-09-10, after a clamp-at-zero
was reverted for exactly that reason).  Blocking there trades a vacuous
TRUE for an avoidable FALSE.

The absence can be established from COMPLETE data instead of a partial
reading, which is the actual requirement.  A voltage source -- of any
kind, recognised or not -- contributes a BRANCH-CURRENT unknown to MNA;
a VCCS, a current source or a nonlinear conductance does not.  Measured:
VS 1, VCVS 1, L 1, VCCS 0, IS 0, R 0.  So `cir.n - len(cir.nodes)` is
the total branch-unknown count, `V` and `L` are the classified elements
that carry one, and any EXCESS is an unclassified element that could be
a voltage source.  No excess means no unrecognised voltage source can
exist -- established from the MNA dimension, which is complete.
⚠ An unclassified element carrying a branch unknown that is NOT a
voltage source (a transformer, an ammeter) blocks the rung too, so this
reports 1 where 0 is true.  Conservative, and `provisional` already
says the reading is partial; an absence asserted wrongly is the failure
that has no floor.

### `algebraic_conditioning`

The docstring before the move:

`(sigma, info)` — how well conditioned the circuit's ALGEBRAIC block is.

⚠ IT LINEARISES AT `x`, AND `G(x)` IS NOT ALWAYS A PURE FUNCTION OF `x`.
An element carrying Newton LIMITING state stamps from solver history as
well: measured on `Diode`, `|G(x)|` moves by 3.6e+02 after a `limit()`
call.  On a diode fixture the verdict and `sigma` here were UNMOVED to
2.2e-16, because the diode sat in the differential part -- but that is one
circuit, not a guarantee.  For a reading that cannot depend on solver
history, pass a converged `x` or reset the circuit's state first.

`sigma` is `sigma_min(d g_2 / d y)`, the smallest singular value of the
block that Bächle 2007 Thm 2.26 requires to have a BOUNDED INVERSE before
a stiffly accurate method with `R(inf) = 0` is entitled to its classical
order on an index-1 DAE.  `None` when there is no such block, or when the
question cannot be answered on this circuit — the info dict says which.

HOW, and why no index-1 splitting is needed.  With `N = ker C` (the
algebraic unknowns) and `Z = ker C^T` (the algebraic equations),
`d g_2/d y = Z^T G N`, and

    sigma_min(C + h G) / h  ->  sigma_min(Z^T G N)    as h -> 0

so two SVDs at different `h` answer the question from `C` and `G` alone,
with no basis extraction and no `(x, y)` form.  ⚠ THE CONVENTION IS OURS:
`J = C + a h G`, not the literature's `C/h + G`, and it flips every
exponent here.  (Identification relayed from the docs session, 2026-09-11;
verified here on six random systems with NON-symmetric `G` and mixed-rank
`C`, recovered to six significant figures, plus both negative controls.)

⚠⚠ THE VERDICT IS FLATNESS, NOT MAGNITUDE, AND THAT IS THE WHOLE POINT.
`sigma_min(C + hG) ~ h` fires whenever an algebraic block EXISTS, which is
essentially every circuit — a voltage source alone is enough.  So the
screen that fires on "`sigma_min ~ h`" is testing EXISTENCE, and existence
is not a defect.  The hypothesis fails only when `sigma_min(Z^T G N) -> 0`,
and THAT shows up as the ratio FAILING TO SETTLE.  Reading the magnitude
instead of the flatness is what made an earlier screen point at the
healthiest fixture we had.

⚠ AND THE RATIO DOES NOT SIMPLY FALL WHEN THE BLOCK IS SINGULAR.  On a
singular block it falls while the leading term is resolved and then GROWS,
because what is left is roundoff divided by `h`: measured 1.3e-07,
5.2e-09, 7.7e-07, 3.0e-05 over four decades.  Anything keyed on "is it
decreasing" reports a singular block as healthy at small enough `h`.
⚠⚠ BUT THOSE NUMBERS ARE FROM RANDOM DENSE SYSTEMS AND YOU WILL NOT SEE
THEM ON A NETLIST.  `sigma_min(C + hG)` is ORTHOGONALLY INVARIANT, so in
exact arithmetic no rotation of the same problem can change it -- the
turn-up is therefore PURE ROUNDOFF and depends on how `ker C` happens to
be represented.  A capacitor-free node gives MNA an EXACT structural zero,
and there the ratio falls cleanly with no turn-up at all (docs-46,
2026-09-11, measured on these same assembled matrices: axis-aligned zero
no turn-up, the same problem rotated turn-up, MNA as assembled no
turn-up).  So gate on FLATNESS, which is right in either basis, and never
on SEEING the turn-up -- on real circuits it is not there to see.

⚠ THERE IS A WINDOW AND IT CAN BE EMPTY.  Above `sigma_r(C)/sigma`, where
`sigma_r(C)` is the smallest NONZERO singular value of `C`, the ratio is
reading the differential directions instead; below `eps*||C||/sigma` it is
reading roundoff.  A badly conditioned `C` leaves no window at all, and
this returns `verdict='no-window'` rather than a number — the one answer
that must never be silently replaced by a plausible-looking value.

`info` carries `verdict`, one of:

``well-conditioned``   the ratio settled; `sigma` is the plateau, and
                         Thm 2.26's hypothesis holds.
``singular``           an algebraic block exists and its smallest
                         singular value is zero to working precision.  This
                         is `theta_0 > 0`, i.e. index >= 2 — cross-check it
                         against :func:`topological_index`.
``no-algebraic-block`` `C` is nonsingular, so there is nothing to condition
                         (the ratio grows like `1/h`; the ODE case).
                         ⚠ This is an ABSENCE claim, and absence claims
                         cannot be made from a PARTIAL reading -- that is
                         what broke `topological_index`'s index-0 rung,
                         whose "no voltage sources" test passed vacuously
                         on an element it could not classify.  This one is
                         safe for a reason worth stating rather than
                         assuming: it is established from the ASSEMBLED
                         `C` by SVD against a relative tolerance, so there
                         is no classifier and nothing can be outside its
                         covered class.  Complete data, not a partial
                         reading.
``no-window``          `C` is too ill-conditioned for a plateau to exist
                         between the turn and the roundoff floor.

plus `h`, `ratio` and `usable` (the probe ladder and which points were not
roundoff), `spread` (max/min over the usable points), and `window`.

⚠ `spread - 1` BOUNDS THE RELATIVE ERROR OF `sigma`, and it is the only
accuracy statement on offer -- so `flat_tol` is not a cosmetic threshold,
it is the worst-case accuracy the caller is agreeing to accept.  Measured
against an explicitly formed `Z^T G N`:

    fixture        sigma          explicit    rel err     spread-1
    RC             0.99005        0.99005     1.98e-12    1.96e-11
    ladder 1e16    0.001          0.001       2.04e-10    1.98e-07
    ladder 1e18    0.001          0.001       2.04e-11    1.98e-08
    e5 G22=2e-04   0.0002         0.0002      2.98e-08    2.95e-07
    e5 G22=2e-06   1.99999e-06    2e-06       3.09e-06    3.06e-05

⚠⚠ THERE IS A RESOLUTION LIMIT, AND A WIDTH WITHOUT ITS PARAMETERS IS
NOT A RESULT.  On the `e5` sweep (a VCCS cancelling a node's
self-conductance, so `G_22` passes through zero with `rank C` FIXED) the
last value still reported is `G_22 = 2e-07` against `||G|| = 2e-03` --
but that `1e-4 * ||G||` holds only AT `flat_tol = 1e-2` AND a 10%
acceptance criterion, and both move it:

    flat_tol   1e-1    1e-2    1e-3    1e-4       accepted within 10%
    last G_22  2e-07   2e-07   2e-06   2e-06
    accepted   50%     10%     1e-3    1e-6       at flat_tol = 1e-2
    last G_22  2e-07   2e-07   2e-07   2e-05

⚠ WHAT DOES *NOT* SET IT IS THE LADDER DEPTH.  Extending `decades` from 8
to 12, 16, 20, 24 and 30 leaves the limit at 2e-07, unmoved to every
digit, because THE FLOOR GUARD CAPS THE USABLE POINT COUNT: 7 points
whether `decades` is 8 or 30, every extra probe falling below
`floor_k * eps * ||C||` and never reaching the flatness test.  That is the
mechanism; the depth-insensitivity is its consequence.  (docs-46 predicted
one decade of limit per decade of ladder, offered the falsification
explicitly, and it FAILED here -- then read this code and identified the
guard, which is the arm that did hold.)

⚠⚠ BOTH KNOBS BIND, ABOUT A DECADE EACH, AND A SINGLE-CAUSE STORY IS
WRONG.  An earlier version of this docstring named `flat_tol` as the
cause; docs-46 then named the floor guard instead, on a transcription
where loosening `flat_tol` bought nothing.  Measured here on an
independent transcription with both exposed, at `decades = 30`:

    shipped (floor_k=1e3, flat_tol=1e-2)    2e-07
    flat_tol loosened 100x to 1e0           2e-08
    floor guard floor_k 1e3 -> 1            2e-08
    both                                    2e-09

so each is worth a decade and they compose.  `flat_tol` is what FIRES --
every run stops with the spread over threshold -- and `floor_k` decides
HOW MANY points the flatness test ever sees.  `floor_k = 1e3` is a
three-decade safety margin over the natural roundoff floor `eps*||C||`,
and it costs about a decade of resolution: a judgement call, exposed as a
parameter so it can be measured rather than argued.  The default stays
high deliberately -- the failure mode on the other side is a CONFIDENT
WRONG NUMBER built from roundoff, which is worse than `singular`.

⚠ One might expect that price to vary by circuit, since the guard is keyed
on `||C||` while the quantity it protects (`smin ~ h * sigma` near the
plateau) carries no `||C||` at all -- so the margin, expressed in units of
what is actually being guarded, moves with `||C||/||G||`.  It does not
bite: the capacitance sweep above varies `||C||/||G||` over twelve decades
at fixed `||G||` and the limit is 2e-07 throughout.  (Raised by docs-46 as
untested; it was already covered by that sweep.)

⚠⚠ AND THAT LIMIT IS INVARIANT TO THE CAPACITANCE UNIT, which is the whole
question.  An absolute rank test on these blocks smears the index-2
crossing into a false window of width `~ tol * ||G||/||C||`, so it widens
as `1/||C||` and at picofarads it is enormous.  Measured here across
twelve decades of `C` -- 1e6, 1e3, 1, 1e-3, 1e-6 times nominal -- the
limit sits at `G_22 = 2e-07` in EVERY case, exactly `1e-4` of `||G||`.
The window is set by `||G||` alone.  Two things make that so and both are
load-bearing: the null space of `C` is taken by a RELATIVE tolerance, and
the verdict is the FLATNESS of the ratio, which is scale-free, rather than
a magnitude compared against a fixed number.

The comment on the limiting-state re-sync:

⚠⚠ RE-SYNC THE LIMITING STATE TO `xv`, AND PUT IT BACK AFTERWARDS.
`G(x)` is NOT a pure function of `x` for a device with a Newton
limiter -- `Diode` linearises around a stored `_vlim` -- so without
this the answer depends on whatever solve ran last.  MEASURED: with a
poisoned `_vlim` this routine's `sigma` moved by a RELATIVE 1.0, while
DC, AC and transient were all unaffected to 0.0e+00 exactly, because a
converged solve leaves the state consistent by construction.  This was
the ONLY site in the tree that inherited it, and it was shipped the
same day the hazard was written down -- which is the argument for the
gate rather than for vigilance.
`limit(x, x)` at ZERO DELTA is the documented re-sync (the same one
`Transient._branch_restore_limits` uses, and the PCNR coupled step at
its own convergence).  ⚠ And the restore is not optional: "A DIAGNOSTIC
THAT CHANGES THE SIMULATION IS A DEFECT, AND THIS ONE DID" is recorded
on that method about `branch_check`, which left `_vlim` at a
speculative solve's value and moved the NEXT step's Jacobian.

The comment on the window choice:

⚠⚠ THE ERROR IS U-SHAPED, SO NEITHER END OF THE LADDER IS THE ANSWER.
Two error sources move in opposite directions: the identification is a
LIMIT, so its truncation falls like `h`; and `C`'s near-null singular
values (plus roundoff) contaminate more as `h` falls.  On a clean
fixture the error decreases all the way down and the LAST point is
best; on a circuit carrying a singular value just under the rank
tolerance the ladder passes THROUGH the true value and climbs again,
and the last point is the WORST.  Two earlier versions read the median
(1000x worse than attainable) and then the last point (wrong by 1e-4 on
the second shape).  Take the FLATTEST 3-POINT WINDOW instead: it finds
the turn wherever it is, and its own spread is the error bound.
⚠ Flatness over the WHOLE ladder is not the test either -- the coarse
end is simply unconverged, and judging on it reported a perfectly
healthy block (ratio converging to 2.0e-04 to ten digits) as SINGULAR
because the top of its ladder was 3% off.


## `_factored.py` -- `FactoredPeriod`

### `FactoredPeriod` (class docstring)

The two warnings as they stood before the history moved out (2026-09-24):

⚠ `matvec` IS THE WHOLE INTERFACE, and deliberately so.  The withdrawn
`PAC` read `pss.Jtvec` / `pss.Cvec` and rebuilt the `(NM)x(NM)` operator
from them, which is where its 419.5 GiB went.  Two things are wrong with
that route and only one of them is the memory: those two lists are
written by `_traverse` and `_traverse_solved_history` and by NEITHER
factored traversal, so after `solve(matrix_free=True)` they are stale or
absent -- an analysis reading them would rebuild an operator for a
DIFFERENT trajectory than the one that converged, silently.

⚠ AND THE REBUILT OPERATOR WAS EULER-SHAPED.  See
`test_the_pac_L_is_backward_euler_only`: the old `L` has two terms per
row, so for `trap` or `gear` it is not the discretisation the trajectory
was produced by and `L^-1 B` is not the monodromy at all -- measured,
spectral radius 0 against the analytic 0.8546.  A `matvec` cannot make
that mistake, because each step carries its own `(alphas, b)`.

### `is_plain`, `is_pair`, `is_stage`, `is_glm`

The comment on the kind flags:

ONE CLASS PER KIND (2026-09-24): `FactoredPeriod(kind, ...)` builds the
subclass that kind names (`_PERIOD_KINDS`), and what differs between
the kinds -- how a direction seeds the per-step state, what the map
reads out of it, where a costate injection lands -- is that class's
own methods, where it was a four-way branch in each of seven.
Consumers ask what a map IS through these flags, not its kind string.

### `__init__`

The comment on `self.times`:

the grid the steps were taken on -- a forced replay needs the
time of each step to evaluate `exp(j w t)` there, and reading it
off `pss.times` later is exactly the parallel-indexing trap the
final replay's `(t, h)` pairing was rewritten to remove

### `matvec_transposed`

⚠ B8: the PLAIN path now has one too, for the ONE-STEP companions.
It used to refuse outright, which made every adjoint surface --
`ppv`, PAC, `pnoise`, `covariance` -- Gear-2 only.

⚠ `collect` AND `inject` FORWARD TO BOTH, which is what lets the
callers stop naming a map. They used to reach past this method to
`_monodromy_matvec_transposed` DIRECTLY, so B8 shipped without
reaching any of them -- the machinery existed and every surface
still refused. A caller that goes through here gets whichever
recursion its `kind` calls for and cannot acquire a Gear-2
assumption by accident.

### `step_objects` (the section comment above it)

-- what differs between the kinds: one subclass each -----------------

(2026-09-23) Every kind's map is `extract . step_N ... step_1 . seed`:
the steps are objects with one algebra (`_StageStep`, `_LMMStep`,
`_GLMStep`: `solve`, `adjoint`, `sources`, `source_adjoint`), and what
is left per kind is how a direction seeds the per-step state, what the
map reads out of it, and where a costate injection lands.  The replays
(`PSS._replay`, `_replay_transposed`, `_forced_replay`,
`_forced_replay_transposed`, `_sideband_forced`) are one function each.
The interface, defined by `_PlainPeriod`, `_PairPeriod`,
`_StagePeriod` and `_GLMPeriod`:

## `_factored.py` -- `_PlainPeriod`

### `_PlainPeriod` (class docstring)

Its `seed` takes ``Pq`` from THE OPENING PAIR, not the loop's -- the
walk opens it right after the MANUFACTURING step, order-dropped to Euler
(``b = 0``), so in practice at zero; using the loop's made the map 100 %
wrong for `trap`.  Plus the consistent-``iq_0`` seed (`_pq_seed_at_x0`,
`theta`), `None` for every other method.


## `_steps.py` -- `_StageStep`

### `_StageStep` (class docstring)

The paragraph as it stood before the history moved out (2026-09-24):

⚠ ONE COPY OF EVERY REPLAY (refactor E9 item 2, 2026-09-23).  The two
structures had a copy each of the forward and transposed mat-vecs, the
forced replay and its transpose, the sideband fold, the stage pass, the
stage-source response and the Lyapunov pieces, dispatched on the map's class.
They differ only inside `solve` and `adjoint`; the rest reads the stage
costates those return and the step's OWN tableau (`A`, `b`, `c` -- the
coupled replays read `par.method`'s, wrong for a `factored_period_stage`
built with another `method=`).

## `_steps.py` -- `_LMMStep`

### `_LMMStep` (class docstring)

The paragraph as it stood before the history moved out (2026-09-24):

⚠ ONE TRANSPOSE FOR EVERY COMPANION (2026-09-23).  The plain map had a
reverse recursion derived for a ONE-STEP companion and refused anything
else; gear's pair had its own, for ``b = 0``, and refused anything else.
Both are this: with ``(w1, w2, w3)`` the adjoints of ``(P_n, P_{n-1},
Pq_n)``,

    Sbar  = w3 - Jf^-T (w1 + a_0 C_n^T w3)
    (P_{n-1}, P_{n-2}, Pq_{n-1})bar = (a_1 C_{n-1}^T Sbar + w2,
                                       a_2 C_{n-2}^T Sbar,  b Sbar)

-- trap's shared bracket (``b != 0``, one-step) and gear's ``(-a_1 C^T t
+ w2, -a_2 C^T t)`` (``b = 0``, ``Sbar = -t``) are its two special cases,
operation for operation.


## `_numerics.py` -- module level

### `_arnoldi_gmres`

The paragraph as it stood before the history moved out (2026-09-24):

SECOND, AND THE REASON THIS IS A CORRECTNESS CHANGE: SciPy REPORTS
BREAKDOWN ON SYSTEMS IT HAS ALREADY SOLVED.  When the Krylov space is
exhausted the next basis vector is numerically zero -- a HAPPY
breakdown, where the answer is EXACT -- and it comes back as
`info = 4`.  Trusting that flag turns an exact answer into a
`RuntimeError`, which is what it did for AM/PM at small offsets and
why `PAC._gmres_checked` exists to overrule it.  Here the breakdown is
detected where it happens and returned as the converged answer it is.

### `periodic_spline_weights`

The three paragraphs as they stood before the history moved out
(2026-09-24):

⚠ UNDER LANDED EVENTS THE PIECES BREAK AT THE EVENT NODES (2026-09-21,
item 2 of the non-uniform-grid list).  A landed edge puts a kink in the
integrand at its node; a spline that is C^2 across it rings, which is
why the callers used to keep the trapezoid whenever `event_times` was
non-empty -- and the trapezoid caps every method's period integrals at
second order on exactly the grids a clocked circuit gets.  A cubic
spline fitted PER SEGMENT (not-a-knot ends; a two-node segment is the
trapezoid, a three-node one the parabola `CubicSpline` builds) never
crosses a kink.  Node 0 is always a break when any event exists: an
edge at the drive's t = 0 is dropped by `event_grid` (the period
boundary is not its to move) and would otherwise sit inside the
periodic seam.  Measured on `exp(sin) + triangle` (kinks landed, a
1 + 0.15 sin grid): the periodic spline +4.0e-5 .. -2.2e-8 with no order
(ringing), the trapezoid +1.8e-5 .. -7.8e-9 at second order, the
piecewise rule -8.6e-7 .. -1.8e-15 over N = 52 .. 1602; on a smooth
integrand with only node 0 broken it is fourth order too (-8.9e-7 ..
-6.1e-13 against the periodic rule's +1.7e-8 .. +2.5e-13), so the seam
break costs a constant, not an order.

The higher-order period quadrature for a NON-UNIFORM, EVENT-FREE grid
(2026-09-21).  The periodic trapezoid rule is spectrally accurate on a
uniform grid and on an alternating one (two interleaved uniform sums) and
genuinely O(h^2) on a smoothly varying grid -- measured: radau's diffusion
constant on a 1 + 0.5 sin grid sat at +4.9e-5 / 1.3e-5 / 3.3e-6 / 8.4e-7
(N = 100 .. 800) while its period, multipliers and mode invariant were at
order 5 / 1e-11 / 1e-15, so every method's noise was capped at second
order by the quadrature on exactly the grids `lte_grid`
produces.  With these weights: -5.5e-9 / -2.4e-10 / -7.3e-12 / +2.6e-12,
the reference's floor from N = 400.

⚠ THE PERIODIC RULE IS FOR EVENT-FREE GRIDS.  A spline is C^2 across
every node; a landed source edge puts a KINK in the integrand at that
node, and a spline through a kink rings where the trapezoid is
exact-ish -- under events the callers pass `breaks` and get the
piecewise rule above (they used to keep the trapezoid).  ⚠ ON A UNIFORM
GRID these equal the trapezoid weights to 1.7e-18 (a periodic spline's
`sum M_j = 0`), and the callers keep their uniform path, which is
bit-identical to before.  The spline system is cyclic tridiagonal and
is solved sparse (O(n)); the weights are `wt - B^T A^{-T} d`, with
`wt` the trapezoid weights and `d_j = (h_j^3 + h_{j-1}^3) / 24` the
coefficient of the node's second derivative in the spline integral.

### `_lu_solve_split` (the comment above its module-level import)

Bound once: `_lu_solve_split` runs once per step of every coupled replay
(~20k calls per PAC + adjoint row on the PWM fixture), and a function-local
import plus `np.iscomplexobj` there was 6 % of that row's profile.
scipy.linalg is loaded by the package's import already.
