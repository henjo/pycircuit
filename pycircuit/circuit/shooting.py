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


def _complex_solve(lu, b):
    """`lu.solve(b)` for a complex `b` against a REAL factorisation.

    Two back-substitutions, not a complex refactorisation: the step
    Jacobians are real, so the solve is real and splits exactly.
    """
    b = np.asarray(b)
    if np.iscomplexobj(b):
        return lu.solve(b.real) + 1j * lu.solve(b.imag)
    return lu.solve(b)


class FactoredPeriod(object):
    """One converged period, kept FACTORED -- the hook PAC/PPV/pnoise share.

    `PSS.solve` throws its factorisations away: they live inside the Newton's
    `build` closure and go out of scope with it.  Every small-signal analysis
    over a periodic operating point needs the same thing back -- products
    with the monodromy, never the monodromy itself -- so this is what
    `PSS.factored_period()` hands out.

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
    """

    __slots__ = ('kind', 'opening', 'steps', 'x_last', 'x_prev', 'width',
                 'times', 'T', 'open_at_x0', '_pss')

    def __init__(self, kind, opening, steps, x_last, x_prev, pss,
                 times=None, T=None, open_at_x0=False):
        self.kind, self.opening, self.steps = kind, opening, steps
        self.x_last, self.x_prev, self._pss = x_last, x_prev, pss
        ## the grid the steps were taken on -- a forced replay needs the
        ## time of each step to evaluate `exp(j w t)` there, and reading it
        ## off `pss.times` later is exactly the parallel-indexing trap the
        ## final replay's `(t, h)` pairing was rewritten to remove
        self.times, self.T = times, T
        ## whether the period opened AT `x_0` -- no manufacturing step.  A
        ## driven replay needs to know: the manufacturing step is not in
        ## `steps`, so its share of the source is not applied, and that
        ## costs an order.  See `PAC.solve`.
        self.open_at_x0 = bool(open_at_x0)
        m = pss.cir.n - 1
        ## the solved-history map acts on the PAIR, the plain map on one
        ## state -- the same distinction `_monodromy` carries
        self.width = 2 * m if kind == 'solved_history' else m

    def matvec(self, v):
        """`M v`, real or complex, replaying the stored factors."""
        if self.kind == 'solved_history':
            return self._pss._monodromy_matvec(self.opening, self.steps, v)
        return self._pss._monodromy_matvec_plain(self.opening, self.steps, v)

    def matvec_transposed(self, v):
        """`M^T v` -- solved-history only; see `_monodromy_matvec_transposed`."""
        if self.kind != 'solved_history':
            raise NotImplementedError(
                'PSS: the transposed replay is implemented for the '
                'solved-history map only, and this period was traversed on '
                "the plain path (method=%r)." % self._pss.par.method)
        return self._pss._monodromy_matvec_transposed(
            self.opening, self.steps, v)


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
         App. A) biases theta off 1/2 so that "numerical oscillations due to
         the DAE character are damped out" -- the same mechanism, attributed
         there to the DAE.

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
        self._dfdT = None
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
        ## What `factored_period()` needs to replay the CONVERGED period:
        ## which seed, which grid, which opening.  Written at the end of
        ## `solve`, and cleared at its start so a failed or interrupted run
        ## cannot leave a previous solution's state readable as this one's.
        self._period_state = None
        self._factored_period_cache = None
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
                           seed_period, solver=None):
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
            if solver is None:
                z, info, ier, mesg = analysis.fsolve(
                    func, z0, maxiter=maxiter, reltol=reltol, abstol=abstol,
                    xtol=xtol, toolkit=self.toolkit, full_output=True,
                    line_search=True)
            else:
                ## ⚠ THE MATRIX-FREE ROUTE COMES THROUGH HERE TOO, so the
                ## trivial-root diagnosis below covers it.  Routing it around
                ## this wrapper would have lost the one message that makes an
                ## autonomous collapse readable.
                z, info, ier, mesg = solver(z0, abstol, xtol, reltol, maxiter)
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
                'between two output recurrences gives one. The literature '
                'remedy for this basin is the PROBE technique -- a periodic '
                'voltage source that feeds the oscillator until its own '
                'current reaches zero, widening the basin so a Newton is '
                'less likely to fall into the DC solution (Bizzarri et al., '
                '"Probe Based Shooting Method ..."); it is not implemented '
                'here, and it widens the basin rather than removing the '
                'seed dependence.'
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
        ## ⚠ AND IT IS ONLY NEEDED WHEN THERE IS A MANUFACTURED STEP TO
        ## PROTECT.  The subdivision exists because `_traverse` builds `x(0)`
        ## with one order-dropped Euler step of `hs[0]` FROM `x_in` -- an
        ## iterate that may be far from the orbit -- and a coarse `hs[0]`
        ## there defeats the inner Newton.  With `x0_unknown` the first step
        ## starts ON the orbit and the same coarse step is solvable:
        ## measured on van der Pol's own LTE grid, the raw 1105-step grid
        ## converges and reaches -47.3 ppm where the subdivided 1106-step
        ## one reaches -73.8.  So the subdivision COSTS accuracy, and it is
        ## skipped where it buys nothing.
        if fr[0] > 8.0 * fr.min() and not getattr(self, '_open_at_x0', False):
            d = float(fr.min())
            fr = np.concatenate(([d, fr[0] - d], fr[1:]))

        ## ⚠ A CALLER'S GRID CAN SILENTLY DEMOTE GEAR-2 TO FIRST ORDER.
        ## `_period_grid` validated positivity and sum-to-1 and nothing
        ## about the INTERIOR ratios.  A two-step method is zero-stable only
        ## up to `h_n / h_{n-1} = 1 + sqrt(2)`, and past it the integrator's
        ## own guard drops the step to Euler -- correct, and invisible.
        ## Measured on a Q=20 resonator driven at resonance with an
        ## alternating 3:1 grid, where half the steps are demoted:
        ##
        ##       npts   uniform    3:1 grid
        ##        100   19.91489    7.99821
        ##        200   20.00960   11.42923
        ##        400   20.02218   14.54985
        ##        800   20.02443   16.85280
        ##
        ## against an analytic peak of 20 V -- 60% low at 100 points,
        ## crawling up at FIRST order, and `converged = True` every time.
        ## Refining does not fix it, because refining a 3:1 grid keeps it
        ## 3:1.  So the warning names the ratio rather than suggesting a
        ## smaller step, which is the advice that does not work here.
        ##
        ## ⚠ This is what item 5 removed the premise for: the literature
        ## note in the class docstring argues Wambacq's objections to
        ## non-uniform BDF "do not bite inside a run" because the grid is
        ## UNIFORM and frozen.  A caller's grid is frozen but not uniform.
        if len(fr) > 1 and self._companion_reach() >= 2:
            from pycircuit.circuit.integrator import ZERO_STABILITY_RATIO
            ratios = fr[1:] / fr[:-1]
            worst = float(np.max(ratios))
            if worst > ZERO_STABILITY_RATIO:
                n_bad = int(np.sum(ratios > ZERO_STABILITY_RATIO))
                warnings.warn(
                    'PSS: this grid steps up by %.3fx where a two-step '
                    'method is zero-stable only to %.3fx, at %d of %d '
                    'interior ratios. Those steps are dropped to Euler by '
                    'the integrator, so the run is first-order there and '
                    'the answer can be far low while reporting converged -- '
                    'measured 60%% low on a Q=20 resonator with an '
                    'alternating 3:1 grid. Refining will NOT fix it: a '
                    'refined 3:1 grid is still 3:1. Smooth the grid so '
                    'adjacent steps stay within %.3fx, or use a one-step '
                    "method (method='trap')."
                    % (worst, ZERO_STABILITY_RATIO, n_bad, len(ratios),
                       ZERO_STABILITY_RATIO), RuntimeWarning, stacklevel=3)

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

    def _step_sensitivity(self, Px, Cs, Pq, Jf, C_new, solve=None,
                          coeffs=None, source=None):
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
        ## `coeffs` overrides the LIVE `_coeffs` for the same reason `solve`
        ## overrides the solve: a matrix-free replay happens after the run
        ## that produced the steps, when `_coeffs` no longer describes the
        ## step being replayed.  See `_traverse_factored_plain`.
        ## ⚠ `source` ENTERS THE SOLVE AND NOT `Pq`, and the asymmetry is
        ## the physics rather than a convenience.  A small-signal source
        ## appears in the step's residual -- `Jf dx + S + du = 0` -- but NOT
        ## in the companion, because `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`
        ## is built from CHARGES, and an injected current is not one.
        ## Adding it to `Pq` as well would feed a fictitious charge forward
        ## into every later step, and the error would grow along the period
        ## rather than announce itself.
        ##
        ## This is what makes PAC share the recursion instead of copying it:
        ## the homogeneous propagation (`source=None`) is the monodromy and
        ## the driven one is the forced response, and they differ by this
        ## one term.
        alphas, b = self._coeffs if coeffs is None else coeffs
        S = b * Pq if b else np.zeros_like(Px[0])
        for k in range(1, len(alphas)):
            S = S + alphas[k] * (Cs[k - 1] @ Px[k - 1])
        S_solve = S if source is None else S + source
        if solve is None:
            ## ⚠ THROUGH THE CALLER'S SOLVER, not `toolkit.linearsolver`.
            ## This is the DENSE propagation -- the thing matrix-free is
            ## measured against -- and it used to be hardcoded to the
            ## toolkit, so `linearsolver=SuperLUSolver()` reached the inner
            ## Newton (once forwarded) and never the propagation.  Comparing
            ## a sparse matrix-free path against a dense baseline would have
            ## flattered it; both sides go through the same strategy now.
            ## `DenseSolver` IS `toolkit.linearsolver`, so the default is
            ## unchanged.
            Px_new = -self._get_linearsolver().solve(Jf, S_solve, self.toolkit)
        else:
            Px_new = -solve(S_solve)
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
                St = St + np.asarray(self._dfdT).ravel() / T
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

    def _factorise(self, Jf):
        """One step's `Jf`, factored by the CALLER'S linear solver.

        ⚠ THIS USED TO REACH FOR `scipy.linalg.lu_factor` DIRECTLY, which
        made every matrix-free run dense-LAPACK whatever `linearsolver=`
        said -- and a circuit Jacobian at m=1002 is very sparse, so a dense
        LU is ~3e8 flops per step, `N` times over.  The recorded 2.13x at
        m=1002 and the m~250 gate were therefore both measured against a
        DENSE baseline; see `benchmarks/pss_matrix_free_ceiling.py` for what
        they become when both sides get a sparse solver.

        A solver whose `factor` returns `None` (the symbolic toolkits, and
        any solver that has not implemented one) falls back to `solve` per
        replayed step -- correct, and `k` times the factorisations, which is
        the cost matrix-free exists to avoid.  It is a fallback, not a mode
        to run in.
        """
        solver = self._get_linearsolver()
        fac = solver.factor(Jf, self.toolkit)
        if fac is not None:
            return fac
        toolkit, A = self.toolkit, Jf

        class _PerSolve(object):
            def solve(self, b):
                return solver.solve(A, b, toolkit)
        return _PerSolve()

    def _traverse_factored(self, x0_in, xm1_in, times, hs, T=None,
                           want_dT=False):
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
        m = self.cir.n - 1
        self._want_dfdh = want_dT
        self._install_history(x0_in, xm1_in, hs[0], h_prev=hs[-1])
        C0 = [np.asarray(self._C_at(x0_in)), np.asarray(self._C_at(xm1_in))]
        Cs = list(C0)
        Pt = [np.zeros(m), np.zeros(m)]
        Pqt = np.zeros(m)
        steps = []
        x, x_prev = copy(x0_in), copy(xm1_in)
        for _j, t in enumerate(times[1:]):
            x_prev = x
            x = copy(self.solve_timestep(x, t, hs[_j]))
            ## ⚠ THE COEFFICIENTS ARE STORED PER STEP, like `Jf` and `C`.
            ## `_coeffs` is live state; on THIS path the history is
            ## installed so no step is order-dropped and the pair is
            ## constant, which is why a single post-run snapshot was right
            ## here.  It was right by luck rather than by construction --
            ## the plain path, whose opening IS dropped to Euler, showed
            ## what that luck is worth.  See `_traverse_factored_plain`.
            alphas, b = self._coeffs
            C_new = np.asarray(self._C).copy()
            lu = self._factorise(np.asarray(self._Jf))
            steps.append((lu, C_new, alphas, b))
            if want_dT:
                ## BOTH ROWS NEED A PERIOD COLUMN -- `dx_{N-1}/dT` and
                ## `dx_{N-2}/dT` -- and the ring carries both, so this is the
                ## same one-column recursion the plain path runs, kept to the
                ## end.  Neither depends on the Krylov direction, so both are
                ## computed once per Newton iteration rather than per GMRES
                ## iteration.
                St = b * Pqt if b else np.zeros(m)
                for k in range(1, len(alphas)):
                    St = St + alphas[k] * (Cs[k - 1] @ Pt[k - 1])
                St = St + np.asarray(self._dfdT).ravel() / T
                Pt_new = -lu.solve(St)
                Pqt = alphas[0] * (C_new @ Pt_new) + St
                Pt = [Pt_new, Pt[0]]
            Cs = [C_new, Cs[0]]
        self._want_dfdh = False
        if want_dT:
            return C0, steps, x, x_prev, Pt[0], Pt[1]
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
        m = self.cir.n - 1
        ## ⚠ COMPLEX `v` IS TWO REAL REPLAYS, NOT A COMPLEX FACTORISATION.
        ## `Jf` and `C` are real, so `M` is a REAL linear map and
        ## `M(a + ib) = Ma + i Mb` exactly.  PAC needs complex products
        ## (`I + alpha(f) H` with `alpha = exp(-2j pi f T)`), and the
        ## alternative -- factoring `Jf` in complex arithmetic -- would
        ## double the stored factors and roughly quadruple the solve cost
        ## for a map that has no imaginary part to begin with.  Splitting
        ## costs two back-substitutions against the SAME factors.
        ##
        ## ⚠ The float cast below is therefore a GUARD, not a convenience:
        ## it used to swallow a complex `v` silently by discarding the
        ## imaginary part, which is a wrong answer rather than an error.
        v = np.asarray(v)
        if np.iscomplexobj(v):
            return (self._monodromy_matvec(C0, steps, v.real)
                    + 1j * self._monodromy_matvec(C0, steps, v.imag))
        v = v.astype(float)
        Px = [v[:m].copy(), v[m:].copy()]
        Cs = list(C0)
        Pq = np.zeros(m)
        for lu, C_new, alphas, b in steps:
            Px_new, Pq = self._step_sensitivity(
                Px, Cs, Pq, None, C_new,
                solve=lambda S, _l=lu: _l.solve(S), coeffs=(alphas, b))
            Px = [Px_new, Px[0]]
            Cs = [C_new, Cs[0]]
        return np.concatenate((Px[0], Px[1]))

    def _traverse_factored_plain(self, x_in, T, times, hs, want_dT=False,
                                 open_at_x0=False):
        """The PLAIN path's trajectory pass, factored -- item 6 for one-step
        methods and for the systems whose unknown is a single entering state.

        Mirrors `_traverse`'s opening exactly, which is the whole
        requirement: the manufacturing step first, then BOTH rings seeded
        with the same `C` and `Pq = a_0 C` for a `b != 0` method.  Getting
        that opening wrong would give a matvec for a DIFFERENT map than the
        dense path's, and the two would disagree only in the third figure --
        the kind of difference a converged answer absorbs.

        With `want_dT` it also propagates the PERIOD COLUMN.  That column
        does not depend on the Newton direction, so it is computed once here
        rather than once per Krylov iteration, and returned as a vector --
        which is exactly what the bordered autonomous matvec needs.
        """
        toolkit = self.toolkit
        m = self.cir.n - 1
        self._want_dfdh = want_dT
        self._begin_period(x_in)
        if open_at_x0:
            ## no manufacturing step -- see `_traverse`, which this mirrors
            x = copy(x_in)
            x0 = copy(x_in)
            C_open = np.asarray(self._C_at(x_in)).copy()
        else:
            x = self.solve_timestep(x_in, times[0], hs[0])
            x0 = copy(x)
            C_open = np.asarray(self._C).copy()

        ## ⚠ THE COEFFICIENTS BELONG TO THE STEP, NOT TO THE RUN.  `_coeffs`
        ## is live state and the MANUFACTURING step is order-dropped, so it
        ## reports Euler's `(alphas, b)` -- `b = 0` -- where the loop steps
        ## report the method's own: measured on this ladder, trapezoidal
        ## opens at `((49000, -49000), 0.0)` and then runs at
        ## `((98000, -98000), -1.0)`.
        ##
        ## Both halves of that were got wrong here first, in opposite
        ## directions, with the trajectory itself matching to ZERO both
        ## times: reading the coefficients ONCE BEFORE THE LOOP applied the
        ## OPENING's to every step and put the period column 40-50% out for
        ## `trap` and `gear`; reading them once in the matvec applied the
        ## LOOP's to the opening and made `Pq` non-zero where `_traverse`
        ## seeds it at zero, which is a 100% error for `trap`.  `euler` was
        ## exact under both and would have passed a one-method test.
        ##
        ## ⚠ WHAT IS LOAD-BEARING IS THE OPENING PAIR.  Inside the loop the
        ## coefficients are CONSTANT for every method in this tree, so
        ## storing them per step is currently belt-and-braces -- a mutation
        ## replacing them with a post-run snapshot does NOT fail the tests,
        ## and that is recorded rather than hidden.  They are stored anyway
        ## because a variable-order method would make the loop vary too, and
        ## that failure would be silent.
        ## ⚠ `b_open = 0` when opening AT `x_0`: no companion current has
        ## been formed, so the matvec must seed `Pq` at zero to match
        ## `_traverse`'s exact seed rather than the manufactured one.
        ## ⚠ `_coeffs` DOES NOT EXIST YET on the `open_at_x0` path -- no step
        ## has run to set it -- so it must not be read at all, not even for
        ## the half that is then discarded.  `b_open = 0` makes the matvec
        ## seed `Pq` at zero, which is what `_traverse` does there, and
        ## `a_open` is unused in consequence.
        a_open, b_open = ((None, 0.0) if open_at_x0 else self._coeffs)
        Pt = [np.zeros(m), np.zeros(m)]
        Pqt = np.zeros(m)
        Cs = [C_open, C_open]
        steps = []
        for _j, t in enumerate(times[1:]):
            dt = hs[_j]
            x = copy(self.solve_timestep(x, t, dt))
            alphas, b = self._coeffs
            Jf = np.asarray(self._Jf)
            C_new = np.asarray(self._C).copy()
            lu = self._factorise(Jf)
            steps.append((lu, C_new, alphas, b))
            if want_dT:
                ## the same recursion with the step size's own source term,
                ## one column wide -- see `_traverse`
                St = b * Pqt if b else np.zeros(m)
                for k in range(1, len(alphas)):
                    St = St + alphas[k] * (Cs[k - 1] @ Pt[k - 1])
                St = St + np.asarray(self._dfdT).ravel() / T
                Pt_new = -lu.solve(St)
                Pqt = alphas[0] * (C_new @ Pt_new) + St
                Pt = [Pt_new, Pt[0]]
            Cs = [C_new, Cs[0]]
        self._want_dfdh = False
        return ((C_open, a_open, b_open), steps, x0, x,
                (Pt[0] if want_dT else None))

    def _monodromy_matvec_plain(self, opening, steps, v):
        """`M v` for the plain path, one column through the stored steps."""
        m = self.cir.n - 1
        C_open, a_open, b_open = opening
        ## ⚠ COMPLEX `v` IS TWO REAL REPLAYS, NOT A COMPLEX FACTORISATION.
        ## `Jf` and `C` are real, so `M` is a REAL linear map and
        ## `M(a + ib) = Ma + i Mb` exactly.  PAC needs complex products
        ## (`I + alpha(f) H` with `alpha = exp(-2j pi f T)`), and the
        ## alternative -- factoring `Jf` in complex arithmetic -- would
        ## double the stored factors and roughly quadruple the solve cost
        ## for a map that has no imaginary part to begin with.  Splitting
        ## costs two back-substitutions against the SAME factors.
        ##
        ## ⚠ The float cast below is therefore a GUARD, not a convenience:
        ## it used to swallow a complex `v` silently by discarding the
        ## imaginary part, which is a wrong answer rather than an error.
        v = np.asarray(v)
        if np.iscomplexobj(v):
            return (self._monodromy_matvec_plain(opening, steps, v.real)
                    + 1j * self._monodromy_matvec_plain(opening, steps, v.imag))
        v = v.astype(float)
        Px = [v.copy(), v.copy()]
        Cs = [C_open, C_open]
        ## ⚠ THE OPENING PAIR, NOT THE LOOP'S.  `_traverse` opens `Pq` at
        ## `a_0 C` for a `b != 0` method reading `_coeffs` right after the
        ## MANUFACTURING step -- which is order-dropped to Euler, where
        ## `b = 0`, so in practice `Pq` opens at ZERO for every method here.
        ## Using the loop's pair instead makes it non-zero and the matvec
        ## 100% wrong for `trap`; this is the half that is load-bearing.
        Pq = a_open[0] * (C_open @ v) if b_open else np.zeros(m)
        for lu, C_new, alphas, b in steps:
            Px_new, Pq = self._step_sensitivity(
                Px, Cs, Pq, None, C_new,
                solve=lambda S, _l=lu: _l.solve(S), coeffs=(alphas, b))
            Px = [Px_new, Px[0]]
            Cs = [C_new, Cs[0]]
        return Px[0]

    def _monodromy_matvec_transposed(self, C0, steps, v, collect=False,
                                     inject=None):
        """`M^T v` -- the same stored steps, REPLAYED BACKWARDS.

        A shooting monodromy is a product of per-step solves, so its
        transpose is that product in reverse order with each solve
        transposed.  For Gear-2 (`b = 0`, three alphas) the per-step state
        is the PAIR `(Px_n, Px_{n-1})` and the step map is

            B_n = [ -a1 Jf^-1 C_{n-1}   -a2 Jf^-1 C_{n-2} ]
                  [        I                    0         ]

        so `B_n^T (v1; v2)` is
        `(-a1 C_{n-1}^T Jf^-T v1 + v2 ; -a2 C_{n-2}^T Jf^-T v1)`.

        ⚠ THIS IS WHY IT COSTS NOTHING TO HAVE.  Demir & Roychowdhury
        (TCAD 22(2) 188-196) call reverse integration "often unavailable
        even in existing time-domain simulators", requiring "significant
        changes to core simulation routines" -- true of a forward-only
        DENSE implementation.  `_traverse_factored` already stores every
        step's factorisation, and every factorisation here already knows how
        to solve transposed, so the reverse pass needs no new integrator, no
        refactorisation and no second traversal.

        MEASURED against the dense `M^T` built from the forward matvec:
        agreement 1.8e-15, and the reverse pass costs 0.75x the forward one
        -- CHEAPER, because it does two `C^T` products against one shared
        transposed solve where the forward does two `C` products and a
        solve.

        ⚠ NOT A CAPABILITY, A BUILDING BLOCK.  It is the shared dependency
        of a PPV (Demir & Roychowdhury's reverse Jacobian `J_r`) and of
        adjoint noise, neither of which is built.  It exists because the
        spike that established it is worth keeping, and it is pinned by a
        test rather than left in a scratch file.
        """
        m = self.cir.n - 1
        ## complex `v` is two real reverse replays -- see the note in
        ## `_monodromy_matvec`; `M^T` is real for exactly the same reason
        ## `M` is, and adjoint noise needs complex products.
        ## ⚠ `collect` HANDS BACK THE PER-STEP TRANSPOSED SOLVES, and it
        ## exists so ADJOINT NOISE does not get a second copy of this
        ## recursion.  The sensitivity of the final state to a source
        ## injected at step `j` is exactly `-Jf_j^-T` applied to the adjoint
        ## state there -- which is `t` below, already computed.  A driven
        ## reverse replay is therefore this pass plus a weighted sum, not a
        ## new traversal, and the two cannot drift apart because there is
        ## only one of them.  See `_forced_replay_transposed`.
        ## ⚠ `inject` MAKES THE OUTPUT A FUNCTIONAL OVER THE WHOLE PERIOD
        ## RATHER THAN A VALUE AT ITS END, and that is the difference
        ## between "the response at `t = 0`" and a SIDEBAND coefficient.
        ## `H_l` is `(1/N) sum_n exp(-j l w0 t_n) d^T y_n`: every state on
        ## the trajectory contributes, so its adjoint takes `c_n` INTO the
        ## reverse state at every step instead of seeding once at the end.
        ## Same recursion, one added term -- which is why this is a
        ## parameter and not a second function to keep in step.
        ##
        ## `inject[j]` lands on the state AFTER `steps[j]` has been applied
        ## backwards, i.e. on `P_j`, and `v` seeds `P_N`.
        v = np.asarray(v)
        inj = None if inject is None else np.asarray(inject)
        cx = np.iscomplexobj(v) or (inj is not None and np.iscomplexobj(inj))
        if cx:
            ir = None if inj is None else inj.real
            ii = None if inj is None else inj.imag
            re = self._monodromy_matvec_transposed(
                C0, steps, np.asarray(v).real, collect, ir)
            im = self._monodromy_matvec_transposed(
                C0, steps, np.asarray(v).imag, collect, ii)
            if collect:
                return (re[0] + 1j * im[0],
                        [a + 1j * b for a, b in zip(re[1], im[1])],
                        [a + 1j * b for a, b in zip(re[2], im[2])])
            return re + 1j * im
        v = v.astype(float)
        w1, w2 = v[:m].copy(), v[m:].copy()
        ts = []
        states = []

        ## The capacitance ring as the FORWARD pass saw it, so the reverse
        ## pass can consume it backwards.  Rebuilt rather than stored: it is
        ## `len(steps)` references to matrices that already exist.
        cs0, cs1, ring = [], [], list(C0)
        for _lu, C_new, _alphas, _b in steps:
            cs0.append(ring[0])
            cs1.append(ring[1])
            ring = [C_new, ring[0]]

        for j in range(len(steps) - 1, -1, -1):
            lu, _C_new, alphas, b = steps[j]
            if len(alphas) < 3:
                ## ⚠ A ONE-STEP COMPANION HAS NO `alphas[2]`, and this used
                ## to read it anyway -- an `IndexError` from inside a
                ## reverse loop, three frames from the cause, where the
                ## `b != 0` case one line below states its refusal
                ## plainly.  Unreachable through `solve`, which takes this
                ## path only for Gear-2, but PPV and adjoint noise both
                ## call the transposed replay directly and a one-step
                ## method is the obvious thing to try first.
                raise NotImplementedError(
                    'PSS: the transposed replay is derived for a two-step '
                    'companion (Gear-2) and this step has %d alpha '
                    'coefficients, so there is no second history term to '
                    'transpose. A one-step method needs its own reverse '
                    'recursion, not this one.' % len(alphas))
            if b:
                ## a `b != 0` companion carries `iq` in the state, so the
                ## step map is not the pair above.  Unreachable today --
                ## `_solves_history` is true only for Gear-2, whose `b` is
                ## zero -- and refused rather than silently transposing a
                ## different operator.
                raise NotImplementedError(
                    'PSS: the transposed replay is derived for a companion '
                    'with b = 0 (Gear-2), and this step reports b = %r. A '
                    'b != 0 method carries iq in the per-step state, so its '
                    'transpose is not the pair map this implements.' % (b,))
            t = lu.solve_transposed(w1)
            if t is None:
                raise NotImplementedError(
                    'PSS: this linear solver cannot solve transposed, so '
                    'the monodromy transpose cannot be replayed. Use '
                    'DenseSolver or SuperLUSolver.')
            if collect:
                ts.append(t)
            w1, w2 = (-alphas[1] * (cs0[j].T @ t) + w2,
                      -alphas[2] * (cs1[j].T @ t))
            if inj is not None:
                w1 = w1 + inj[j]
            if collect:
                ## ⚠ THE ADJOINT STATE AFTER THE STEP, which for a seed of
                ## the PPV IS the PPV at that time -- `Phi(T,s_j)^T v(T) =
                ## v(s_j)`.  Oscillator phase noise needs `v(t)` over the
                ## whole period, not `v(0)`, and it is already being
                ## computed here and thrown away.
                states.append(np.concatenate((w1.copy(), w2.copy())))
        if collect:
            ## reversed, so `ts[j]` lines up with `steps[j]` -- the reverse
            ## loop produced them last-first and a caller weighting them by
            ## a per-step time must not have to remember that
            ts.reverse()
            states.reverse()
            return np.concatenate((w1, w2)), ts, states
        return np.concatenate((w1, w2))

    ## How many deflated power iterations estimate the second multiplier.
    ## Convergence is at |lambda_3|/|lambda_2|, which is fast in the case
    ## that matters -- a lone slow node leaves everything below it tiny.
    PPV_DEFLATION_ITERS = 30
    ## Above this, the bordered extraction is losing digits AND the phase
    ## equation's instantaneous-response assumption is in doubt.
    PPV_SECOND_MULTIPLIER_WARN = 0.9

    def ppv(self, tol=None):
        """The perturbation projection vector at `t = 0` (Demir & Roychowdhury).

        Returns `(v, info)`.  `v` is the pair-space left null vector of
        `I - M`, normalised so that `v . xdot(0) = 1` -- see the note on the
        normalisation below, which was MEASURED rather than transcribed.
        The phase shift caused by a state perturbation `delta` at `t = 0` is
        then `v[:m] . delta`.  `info` carries both border residuals, the
        null residual, `q` and the scaled tangent.

        ⚠ AN AUGMENTED SOLVE, NOT AN EIGENVECTOR, and the difference is the
        whole content of the 2003 paper.  Demir's 2000 method SELECTED the
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
        built.  (Cited from the docs session, not verified here.)

        ⚠ `q` IS EXACT, NOT DIFFERENCED.  `q = C(0) xdot(0)` looks like it
        needs the orbit's tangent, and differencing the waveform for it
        would be O(h) at best.  The DAE gives it directly: `dq/dt + i(x) +
        u(t) = 0` and `dq/dt = C xdot`, so `q = -(i(x_0) + u(0))` -- two
        evaluations at the converged solution, no derivative anywhere.
        """
        import scipy.sparse.linalg as spla
        fp = self.factored_period()
        if fp.kind != 'solved_history':
            raise NotImplementedError(
                'PSS.ppv: the augmented solve needs the transposed replay, '
                "which is implemented for the solved-history map only. "
                "Re-solve with method='gear'. (method=%r.)" % self.par.method)
        if not self.autonomous:
            raise ValueError(
                'PSS.ppv: a perturbation projection vector describes the '
                'phase of a FREE-RUNNING oscillator. This circuit is '
                'driven, so its phase is the source\'s and there is no '
                'unit Floquet multiplier to project onto.')

        m = self.cir.n - 1
        n = fp.width
        irn = self.irefnode
        x0r = np.asarray(self._period_state[1], dtype=float).ravel()
        x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
        qf = -(np.asarray(self.cir.i(x0f)).ravel()
               + np.asarray(self.cir.u(0.0,
                                       analysis=self.par.analysis)).ravel())
        q = np.delete(np.asarray(qf, dtype=float), irn)
        qp = np.concatenate((q, np.zeros(n - m)))
        nq = float(np.linalg.norm(qp))
        if nq == 0.0:
            raise ValueError(
                'PSS.ppv: C(0) xdot(0) is zero, so the orbit has no tangent '
                'at t=0 and there is nothing to normalise against. That '
                'should not happen on a converged limit cycle.')

        def _mv(z):
            z = np.asarray(z)
            v_, y_ = z[:n], z[n]
            top = v_ - self._monodromy_matvec_transposed(
                fp.opening, fp.steps, v_) + y_ * qp
            return np.concatenate((top, [float(qp @ v_)]))

        rtol = max(self.par.reltol * 1e-2 if tol is None else tol, 1e-14)
        A = spla.LinearOperator((n + 1, n + 1), matvec=_mv, dtype=float)
        rhs = np.zeros(n + 1)
        rhs[n] = 1.0
        z, code = spla.gmres(A, rhs, rtol=rtol, restart=min(n + 1, 50),
                             maxiter=min(n + 1, 200))
        if code != 0:
            raise RuntimeError(
                'PSS.ppv: the augmented solve did not converge (info=%r). '
                'The border is `q` itself; if the orbit tangent is nearly '
                'orthogonal to the null direction the bordering is poor.'
                % code)
        v, y = z[:n], float(z[n])
        resid = float(np.linalg.norm(
            v - self._monodromy_matvec_transposed(fp.opening, fp.steps, v)))

        ## ⚠ THE SCALE NEEDS THE TANGENT, so the RIGHT null vector is solved
        ## for too -- by the same bordering, not by an eigendecomposition,
        ## for the same reason: on a high-Q oscillator the other multipliers
        ## crowd 1 and no selection among candidates is reliable.
        def _mvf(zz):
            zz = np.asarray(zz)
            u_, yy = zz[:n], zz[n]
            top = u_ - self._monodromy_matvec(
                fp.opening, fp.steps, u_) + yy * qp
            return np.concatenate((top, [float(qp @ u_)]))

        Af = spla.LinearOperator((n + 1, n + 1), matvec=_mvf, dtype=float)
        zf, codef = spla.gmres(Af, rhs, rtol=rtol, restart=min(n + 1, 50),
                               maxiter=min(n + 1, 200))
        if codef != 0:
            raise RuntimeError(
                'PSS.ppv: the tangent solve did not converge (info=%r).'
                % codef)
        u, yf = zf[:n], float(zf[n])

        ## `u` is the tangent's DIRECTION; its scale comes from `C u = q`,
        ## which is the definition of `q` read backwards.  Least squares
        ## because `C` is singular for a DAE and only its range is
        ## determined.
        x0red = np.asarray(self.cir.C(x0f))
        Cm = np.delete(np.delete(x0red, irn, 0), irn, 1)
        Cu = np.asarray(Cm, dtype=float) @ u[:m]
        denom = float(Cu @ Cu)
        if denom == 0.0:
            raise ValueError(
                'PSS.ppv: C(0) annihilates the orbit tangent, so its scale '
                'is not determined by C u = q.')
        xdot = u[:m] * (float(q @ Cu) / denom)

        ## ⚠ NORMALISED BY `v . xdot = 1`, WHICH IS NOT WHAT `v . q = 1`
        ## GIVES, and the difference is not cosmetic -- measured 0.93 apart
        ## on van der Pol, i.e. a 7% error in every phase-noise number
        ## downstream.  The defining property is that displacing the state
        ## ALONG the orbit by `eps xdot` advances the phase by `eps`, so
        ## `v . xdot = 1` is the normalisation a state perturbation sees.
        ## Demir's Remark 3.1 reads `v_1^T C u_1 = 1`; the vector this
        ## bordered solve returns behaves as `C^T v_1` -- it is contracted
        ## with a state perturbation directly -- so the two statements agree
        ## about different objects.  ⚠ TREATING THEM AS THE SAME OBJECT WAS
        ## MEASURED WRONG: predicting a state jump's phase shift as
        ## `v^T C delta` gives residuals of 0.36/0.40/0.42 that GROW with
        ## refinement and per-direction ratios scattering from -0.44 to
        ## 28.7, while `v . delta` converges at O(h).
        vx = float(v[:m] @ xdot)
        if vx == 0.0:
            raise ValueError(
                'PSS.ppv: the null vector is orthogonal to the orbit '
                'tangent, so no normalisation makes it a phase projector.')
        v = v / vx
        ## ⚠ THE PPV OVER THE PERIOD, not just at `t = 0`, because that is
        ## what an oscillator noise calculation needs: Demir's diffusion
        ## constant is `c = (1/T) integral v_1^T(t) B(t) B^T(t) v_1(t) dt`,
        ## an integral over the orbit.  `Phi(T,s)^T v(T) = v(s)`, and the
        ## reverse replay computes exactly that sequence on its way to the
        ## answer -- it was being discarded.
        _end, _ts, states = self._monodromy_matvec_transposed(
            fp.opening, fp.steps, v, collect=True)
        ## ⚠ A SECOND MULTIPLIER NEAR 1 BREAKS THIS SILENTLY, and none of
        ## the residuals above can see it.  The border removes the PHASE
        ## mode's singularity and does nothing about any OTHER root
        ## approaching the unit circle -- which a slow node puts there.
        ## MEASURED on van der Pol with one weakly coupled RC node:
        ##
        ##     tau/T    |lambda_2|    sigma_min(bordered)   null residual
        ##     none      0.000856          8.62e-01            4.1e-11
        ##     1e2       0.990049          4.49e-03            4.6e-11
        ##     1e4       0.999900          4.47e-05            4.6e-11
        ##     1e6       0.999999          4.47e-07            4.4e-11
        ##
        ## `sigma_min` tracks `T/tau` over six decades while the residual
        ## does not move at all: GMRES converges, the answer looks clean,
        ## and the conditioning has lost six digits.  So this estimates
        ## `|lambda_2|` explicitly rather than trusting a small residual.
        ##
        ## ⚠ AND THE ACCURACY COST WAS GATED, WITH A NEGATIVE RESULT worth
        ## recording so nobody re-derives a fix from the warning alone.
        ## Monte Carlo on the FULL NONLINEAR circuit -- 200 realisations,
        ## 150 periods, phase read from zero-crossing timing, so no PPV
        ## appears anywhere in the measurement:
        ##
        ##     core injection (control)   c_ppv/c_meas = 0.9965
        ##     slow node, tau/T = 10      c_ppv/c_meas = 0.8016
        ##
        ## Within 20%, about 2 sigma at this sample count, and in the
        ## UNDER-predicting direction.
        ##
        ## ⚠ BUT THAT IS NOT A FALSIFICATION, AND THIS DOCSTRING SAID IT
        ## WAS.  `tau/T = 10` is OUTSIDE the regime the reported mechanism
        ## needs: it bites through ill-conditioning, and by the table above
        ## `sigma_min` at `tau/T = 10` is ~4.5e-02 -- healthy.  The PPV has
        ## no large entries there and nothing is splitting into two nearly
        ## cancelling components.  Lai's own case is a gated-capacitor
        ## tuning bank whose off-caps have RC "larger than 1 second" at
        ## 3.15 GHz, i.e. `tau/T ~ 3e9` -- eight orders from what was
        ## tested.  A null result at 10 is what the mechanism PREDICTS, not
        ## evidence against it.
        ##
        ## ⚠ SO THE HONEST RECORD IS: not reproduced at `tau/T = 10`, which
        ## is outside the regime where the mechanism predicts an effect;
        ## UNTESTED at the `tau/T ~ 1e9` where it is reported.  And the
        ## reason the fix is still not built is COST, not falsification:
        ## the measurement needs ~15 time constants of settling, so at
        ## `tau/T = 1e4` that is 150 000 periods per realisation.  That
        ## argument stands on its own; the falsification framing does not,
        ## and this codebase's ledger distinguishes them.
        ##
        ## ⚠ AND THE 0.80 IS IN THE OPPOSITE DIRECTION TO THE REPORTED
        ## EFFECT.  If it survives the ~10% Monte Carlo uncertainty at 200
        ## realisations it is a separate ~20% UNDER-prediction at a `tau/T`
        ## where the conditioning is fine -- not a weak version of Lai's.
        ## At ~2.5 sigma it is not established either way, and it is
        ## recorded rather than resolved.
        ##
        ## ⚠ Larger `tau/T` is untested and the cost is why: the
        ## measurement needs ~15 time constants of settling.
        ## ⚠ AND IT TOOK THREE ATTEMPTS.  A window of 2-4 time constants
        ## read the slow mode's DECAY as diffusion; an impulse test could
        ## not resolve a 1e-11 time shift; and one noise amplitude for both
        ## circuits put a 2.5 V jump per step on an orbit of amplitude 2,
        ## because the slow node's capacitance is 6.7e-5 F against the
        ## core's 1.0.  Each time the number was read before the
        ## MEASUREMENT was shown to be in the regime it assumes.
        ##
        ## Deflated power iteration, using the eigenvectors already in
        ## hand: `u` spans the unit mode and `v` is its left partner, so
        ## `w -> M w - (v.Mw / v.u) u` removes it exactly.
        vu = float(v[:m] @ u[:m] + v[m:] @ u[m:])
        lam2 = 0.0
        if abs(vu) > 1e-300:
            rng = np.random.default_rng(12345)
            w = rng.standard_normal(n)
            for _k in range(self.PPV_DEFLATION_ITERS):
                w = w - (float(v @ w) / vu) * u
                nw = float(np.linalg.norm(w))
                if nw < 1e-300:
                    break
                w = w / nw
                Mw = np.asarray(self._monodromy_matvec(
                    fp.opening, fp.steps, w))
                lam2 = float(np.linalg.norm(Mw))
                w = Mw
        if lam2 > self.PPV_SECOND_MULTIPLIER_WARN:
            warnings.warn(
                'PSS.ppv: a SECOND Floquet multiplier sits at %.6f, near '
                'the unit circle. The bordered extraction removes only the '
                'phase mode, so its conditioning degrades as that root '
                'approaches 1 -- measured losing six digits over six '
                'decades of time constant while every residual stayed at '
                '1e-11. ⚠ AND THE PHASE EQUATION ITSELF IS THE DEEPER '
                'ISSUE: it treats the frequency response as instantaneous, '
                'so slow nodes that FILTER a device\'s noise are not seen '
                'and phase noise is OVER-ESTIMATED. Neither a smaller '
                'tolerance nor a better extraction fixes that; it needs a '
                'frequency-aware PPV. Treat this result as an upper bound.'
                % lam2, RuntimeWarning, stacklevel=2)
        ## ⚠ ONE NUMBER THAT SUBSUMES FOUR DIAGNOSTICS.  An amplitude
        ## perturbation decays to `|lambda_2|` of its size each cycle, so
        ## the cycles needed to fall below a threshold IS the oscillator's
        ## Q: `Q = log(threshold)/log|lambda_2|` (Wang & Roychowdhury).
        ## The usual definitions do not apply to an autonomous circuit --
        ## `f_r/df` presumes a Bode plot of a BIBO-stable linear system,
        ## and stored/dissipated presumes damping a self-sustaining
        ## oscillator does not have.  Nor is it the resonator's Q.
        ##
        ## ⚠ AND IT IS THE SAME CONDITION AS EVERY FAILURE THIS CLASS
        ## WARNS ABOUT.  "High Q", "a second multiplier near 1", "slow
        ## amplitude restoration" and "a long time constant" are four
        ## vocabularies for one thing -- which is why the same circuits
        ## defeat the phase row, the eigen-split, the probe's continuation
        ## and the PPV's instantaneous-response assumption.  Not four
        ## coincidences.  It costs nothing here: the deflated power
        ## iteration above already produced `|lambda_2|`.
        ##
        ## Reported for a `1/e` threshold, so `Q` is cycles-to-1/e.
        Q = (-1.0 / np.log(lam2) if 0.0 < lam2 < 1.0 else float('inf'))
        info = {'border_residual': y,
                'tangent_border_residual': yf,
                'Q': Q,
                'null_residual': resid / max(float(np.linalg.norm(v)), 1e-300),
                'second_multiplier': lam2,
                'q': q, 'xdot': xdot, 'tangent_pair': u,
                'samples': np.asarray(states),
                'times': np.asarray(fp.times, dtype=float)}
        return v, info

    def _forced_replay_transposed(self, fp, freq, xa):
        """`W^T xa` -- the transpose of the map `u -> w(freq)`.

        THE MANY-TO-ONE HALF, and the reason adjoint noise is affordable.
        The forward replay answers "what does THIS source do at the
        output"; one run per source.  This answers "what does the output
        owe to EVERY source", in one reverse pass -- which is the shape
        pnoise has, with hundreds of sources and one output.  Okumura et
        al. (1993) choose the adjoint for exactly this: "it is efficient to
        use the adjoint method ... because circuits have many noise
        sources."

        The identity is one line.  The forward replay makes each step
        `Px_j = -Jf_j^-1 (S_j + u e^{jw t_j})`, so the final state's
        sensitivity to `u` through step `j` is `-Jf_j^-T` applied to the
        adjoint state there, weighted by `e^{jw t_j}`.  That solve is
        already taken by the reverse pass -- it is `t` in
        `_monodromy_matvec_transposed` -- so this is that pass plus a
        weighted sum, with no second recursion to keep in step.

        ⚠ SOLVED-HISTORY ONLY, because the reverse pass is.  See there.
        """
        _end, ts, _states = self._monodromy_matvec_transposed(
            fp.opening, fp.steps, xa, collect=True)
        jw = 2j * np.pi * float(freq)
        acc = np.zeros(self.cir.n - 1, dtype=complex)
        for tvec, t in zip(ts, fp.times[1:]):
            acc = acc - np.exp(jw * float(t)) * np.asarray(tvec)
        return acc

    def factored_period(self):
        """The converged period's steps, kept factored -- see `FactoredPeriod`.

        Runs the factored traversal ONCE, at the solution, and caches it.
        Lazy on purpose: `_traverse_factored` stores `N` factorisations and
        `N` capacitances (`2 N m^2` doubles -- ~800 MB at m=1002 and 50
        points), which is a bad trade to impose on every `solve` for the
        callers who never ask.

        ⚠ IT RE-TRAVERSES RATHER THAN REUSING THE NEWTON'S FACTORS, and the
        difference is not efficiency.  The last `build` call inside the
        Newton is at the last TRIAL iterate; the converged answer is the one
        after it.  Reusing those factors would give an operator for a
        trajectory near the solution instead of at it -- a small error, in
        the third figure, of exactly the kind a converged answer absorbs
        without complaint.
        """
        if getattr(self, '_period_state', None) is None:
            raise RuntimeError(
                'PSS: no period to factor -- call solve() first. '
                '(factored_period() replays the CONVERGED trajectory, so '
                'there has to be one.)')
        if not self.converged:
            raise RuntimeError(
                'PSS: the shooting solve did not converge, so there is no '
                'periodic operating point to linearise about. A '
                'small-signal analysis over a non-solution is not a '
                'meaningful answer -- fix the PSS run first.')
        if self._factored_period_cache is not None:
            return self._factored_period_cache

        solved, x0, xm1, times, hs, T, x0_unknown = self._period_state
        if solved:
            C0, steps, x_last, x_prev = self._traverse_factored(
                x0, xm1, times, hs, T=T)
            fp = FactoredPeriod('solved_history', C0, steps, x_last, x_prev,
                                self, times=times, T=T)
        else:
            opening, steps, x0_out, x_last, _dT = \
                self._traverse_factored_plain(x0, T, times, hs,
                                              open_at_x0=x0_unknown)
            fp = FactoredPeriod('plain', opening, steps, x_last, x0_out,
                                self, times=times, T=T,
                                open_at_x0=x0_unknown)
        self._factored_period_cache = fp
        return fp

    def _forced_replay(self, fp, freq, u_ac, y0=None, collect=False):
        """One period of the LINEARISED circuit, driven at `freq`.

        The same recursion as `_monodromy_matvec`, with `source` switched
        on: `_step_sensitivity` is linear in `(y0, source)`, so this returns

            y_end = M y0 + w(freq)

        with `w` the particular response from a zero initial state.  That
        superposition is not an incidental property -- it is what lets PAC
        solve an `m x m` system instead of an `(N m) x (N m)` one, and it is
        asserted in the suite rather than assumed.

        ⚠ THE SOLVE IS REAL, THE REPLAY IS COMPLEX.  `Jf` is real and stays
        factored once; a complex right-hand side costs two back-substitutions
        against those same factors.  See the note in `_monodromy_matvec`.
        """
        m = self.cir.n - 1
        jw = 2j * np.pi * float(freq)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()

        if fp.kind == 'solved_history':
            if y0 is None:
                Px = [np.zeros(m, dtype=complex), np.zeros(m, dtype=complex)]
            else:
                y0 = np.asarray(y0, dtype=complex).ravel()
                Px = [y0[:m].copy(), y0[m:].copy()]
            Cs = list(fp.opening)
            Pq = np.zeros(m, dtype=complex)
        else:
            C_open, a_open, b_open = fp.opening
            v = (np.zeros(m, dtype=complex) if y0 is None
                 else np.asarray(y0, dtype=complex).ravel().copy())
            Px = [v.copy(), v.copy()]
            Cs = [C_open, C_open]
            ## the OPENING pair, exactly as `_monodromy_matvec_plain` -- see
            ## the note there; using the loop's makes it wrong for `trap`
            Pq = (a_open[0] * (C_open @ v) if b_open
                  else np.zeros(m, dtype=complex))

        ys = []
        for (lu, C_new, alphas, b), t in zip(fp.steps, fp.times[1:]):
            src = u_ac * np.exp(jw * float(t))
            Px_new, Pq = self._step_sensitivity(
                Px, Cs, Pq, None, C_new,
                solve=lambda S, _l=lu: _complex_solve(_l, S),
                coeffs=(alphas, b), source=src)
            Px = [Px_new, Px[0]]
            Cs = [C_new, Cs[0]]
            if collect:
                ys.append(Px_new.copy())

        end = (np.concatenate((Px[0], Px[1])) if fp.kind == 'solved_history'
               else Px[0])
        return end, ys

    ## How hard GMRES is asked to solve, relative to the shooting tolerance.
    ## An inexact Newton only needs the step accurate enough not to spoil
    ## the outer convergence; measured on the RC ladder, k is 2/4/7/12 at
    ## m=40/110/242/502 against this factor and the whole system's spectrum
    ## explains why -- `I - M` has almost every eigenvalue within 1% of 1.0
    ## because the fast modes decay to nothing over a period, leaving only
    ## the slow ones for GMRES to resolve.  So k tracks the number of SLOW
    ## MODES, not m, which is the property the item rests on.
    KRYLOV_TOLERANCE_FACTOR = 1e-2

    ## ⚠ THE BUDGET IS A CHOICE AND SCIPY'S UNITS ARE A TRAP.  `maxiter`
    ## counts RESTART CYCLES, not matvecs, so the pair multiplies: the
    ## earlier `restart=min(n, 200), maxiter=min(n, 400)` was a worst case
    ## near 80 000 matvecs, each a full N-step replay of the period, with no
    ## diagnostic when it was being spent.  Measured k is 2-12 on circuits
    ## whose `I - M` clusters at 1, which is the property the whole method
    ## rests on, so 200 x 20 is already four thousand times what a
    ## well-behaved system needs; a circuit that exceeds it is telling you it
    ## does not cluster, and the answer is the dense path, not a bigger
    ## budget.
    KRYLOV_RESTART = 200
    KRYLOV_MAX_CYCLES = 20

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
        m = self.cir.n - 1

        def build(z):
            C0, steps, x_last, x_prev = self._traverse_factored(
                z[:m], z[m:], times, hs)
            F = np.concatenate((z[:m] - np.asarray(x_last),
                                z[m:] - np.asarray(x_prev)))
            return F, (lambda v: v - self._monodromy_matvec(C0, steps, v))

        return self._matrix_free_newton(build, z0, abstol, xtol, reltol,
                                        maxiter)

    def _matrix_free_newton(self, build, z0, abstol, xtol, reltol, maxiter):
        """The Newton loop every matrix-free system shares.

        `build(z)` returns `(F, matvec)` for the current iterate -- one
        trajectory pass, then a linear operator that never forms its matrix.
        Written once because the four systems differ ONLY in those two
        things: the plain path's `I - M`, the solved-history path's `2m`
        pair, and the bordered autonomous versions of each.
        """
        import scipy.sparse.linalg as spla
        z = np.asarray(z0, dtype=float).copy()
        n = len(z)
        ier, mesg, xdiff = 2, 'No convergence', None
        for _i in range(maxiter):
            F, mv = build(z)

            def _mv(v, _f=mv):
                return _f(v)

            J = spla.LinearOperator((n, n), matvec=_mv, dtype=float)
            xdiff, info = spla.gmres(
                J, -F, rtol=self.KRYLOV_TOLERANCE_FACTOR * reltol,
                restart=min(n, self.KRYLOV_RESTART),
                maxiter=self.KRYLOV_MAX_CYCLES)
            ## ⚠ THE INNER SOLVE'S VERDICT IS NOT DISCARDED.  It used to be,
            ## and a Krylov breakdown then surfaced as the generic outer
            ## 'No convergence' with nothing naming the cause -- in a file
            ## whose whole standard is that a failure says what happened
            ## (`T = 0`, the trivial root, the singular free-period
            ## Jacobian).  An unconverged GMRES makes `xdiff` a direction
            ## the Newton has no reason to trust, so the outer loop is told
            ## to stop rather than iterate on it.
            if info != 0:
                warnings.warn(
                    'PSS: the matrix-free inner solve did not converge at '
                    'outer iteration %d -- GMRES returned info=%d (%s) on a '
                    '%d-unknown system, after at most %d matvecs, each a '
                    'full replay of the period. The Newton step it returned '
                    'is not a direction worth iterating on, so this solve '
                    'stops here and reports not-converged. Measured k on '
                    'well-behaved circuits is 2-12 because `I - M` clusters '
                    'at 1; needing more than %d means this system does not '
                    'cluster, and the dense path (matrix_free=False) is the '
                    'reliable answer for it.'
                    % (_i, info,
                       'breakdown' if info < 0 else 'iteration limit',
                       n, min(n, self.KRYLOV_RESTART) * self.KRYLOV_MAX_CYCLES,
                       min(n, self.KRYLOV_RESTART) * self.KRYLOV_MAX_CYCLES),
                    RuntimeWarning, stacklevel=3)
                return z, {}, 2, 'No convergence (inner Krylov solve failed)'
            z_new = z + xdiff

            ## `|J| . |x|` is not available without the matrix; see
            ## `_matrix_free_solve` for what this substitute is and is not.
            I_scale = np.abs(z_new) + np.abs(mv(z_new)) + np.abs(F)
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
        back identically zero on every MNA circuit tried here.

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

    def _traverse(self, x_in, T, times, hs, want_dT, open_at_x0=False):
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
        if open_at_x0:
            ## ⚠ NO MANUFACTURING STEP: the caller's unknown IS `x_0`.
            ## `_begin_period` has seeded both charge rings from `q(x_in)`
            ## and marked the next step `is_first_step`, so the first step
            ## INSIDE the period is order-dropped to Euler -- which is the
            ## L-stable opener this formulation cannot do without.  See
            ## `x0_unknown` on `solve` for why, and for what it costs.
            x = copy(x_in)
            x0 = copy(x_in)
            C_open = np.asarray(self._C_at(x_in))
        else:
            x = self.solve_timestep(x_in, times[0], hs[0])
            x0 = copy(x)
            C_open = np.asarray(self._C)
        ## ⚠ `None`, not `self._iq`, when nothing has been solved yet: no
        ## companion current exists at `x_0` because no step has formed one.
        ## `solve_timestep` reads `iq_last=None` as "open the run", which is
        ## the same restart `_begin_period` already asked for -- passing a
        ## stale `_iq` from a previous traversal would be the hidden-state
        ## defect this class refuses elsewhere.
        iq_last = None if open_at_x0 else self._iq

        ## `Px[k]` is d(x_{n-k})/d(x0), `Pq` is d(iq_n)/d(x0), `Cs[k]` the
        ## capacitance of step n-k.  Two of each -- as far back as any
        ## method here reaches.  Both rings open seeded with the entering
        ## step, mirroring how the transient seeds `_qlast` with `q0`
        ## repeated: at the start of a period there is no earlier point to
        ## differentiate against.
        eye = np.asarray(toolkit.eye(n - 1))
        Px = [eye, eye]
        Cs = [copy(C_open), copy(C_open)]
        if open_at_x0:
            ## ⚠ AND HERE THE SEED IS EXACT RATHER THAN ASSUMED.  With `x_0`
            ## the unknown and both rings holding `q(x_0)`, `dq_{-1}/dx_0`
            ## really IS `C` -- the same `I` the other branch writes down as
            ## an assumption about a history it did not solve for.  `Pq` is
            ## zero because no companion current has been formed yet: the
            ## opening step has not been taken.  That is the whole reason
            ## this path has an exact Jacobian and the other does not.
            Pq = np.zeros((n - 1, n - 1))
        else:
            a_first, b_first = self._coeffs
            Pq = (a_first[0] * Cs[0] if b_first else np.zeros((n - 1, n - 1)))
        ## The period column, propagated the same way with one extra source
        ## term.  Zero at the start: the entering state does not depend on T.
        Pt = [np.zeros(n - 1), np.zeros(n - 1)]
        Pqt = np.zeros(n - 1)

        ## Kept for PAC.
        ## ⚠ `C_open`, not `self._C`: on the `open_at_x0` path no step has
        ## run, so `_C` does not exist yet.  `C_open` is the same matrix the
        ## other branch would have found there, evaluated at `x_0` directly.
        ## (Kept for PAC, which is withdrawn; `Jtvec` opens empty for the
        ## same reason -- there is no solved `Jf` at `x_0`.)
        self.Cvec = [copy(C_open)]
        self.Jtvec = ([] if open_at_x0 else [copy(self._Jf)])
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

            ## ⚠ THIS USED TO BE AN INLINE COPY of `_step_sensitivity`,
            ## byte-for-byte, while that method's own docstring said the two
            ## systems "share this and not a copy".  They did not: only the
            ## solved-history path called it.  Found while measuring the
            ## plain path's propagation share -- a timer wrapped around
            ## `_step_sensitivity` reported 0.0%, which is not a small
            ## number but a wrong one.
            Px_new, Pq = self._step_sensitivity(Px, Cs, Pq, Jf, C_new)

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
                St = St + np.asarray(self._dfdT).ravel() / T
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

        ⚠ THAT LIST IS 1-FOR-4 AS SHIPPED, and an external review counted it
        (2026-09-02).  Verified here: LIMITING does reach -- `cir.limit` is
        called on the inner Newton and the rectifier measurably conducts,
        which is the defect that motivated driving `Transient` at all.  The
        other three do not.  `PSS(cir, pcnr=True)` raises `KeyError` because
        this class declares no `pcnr` Parameter; the continuation rescue
        (`_rescue_solver`) and breakpoints (`cir.next_event`) are armed only
        inside `Transient.solve`, which PSS never calls -- it drives
        `solve_timestep` directly, and imposes its own grid, so a breakpoint
        has nothing to move.  The same structural fact behind the TLine
        refusal above: what `Transient.solve` does per accepted step, PSS
        does not do at all.

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
                TRTOL=self.par.TRTOL, relref=self.par.relref,
                ## ⚠ THE SOLVER STRATEGIES GO THROUGH TOO, and they used not
                ## to.  `nrsolver`, `linearsolver` and `scaler` are declared
                ## on the base `Analysis`, so `PSS(cir, linearsolver=...)`
                ## has always been ACCEPTED -- and then dropped here, with
                ## the inner `Transient` resolving to `DenseSolver` /
                ## `StandardNewton` whatever the caller asked for.  That is
                ## the third time this class has taken a parameter it never
                ## read (`method` declared and never read; `analysis='PSS'`
                ## matching nothing), and the same shape each time: accepted
                ## at the constructor, silently discarded at the boundary.
                nrsolver=self.par.nrsolver,
                linearsolver=self.par.linearsolver,
                scaler=self.par.scaler)
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
            ## ⚠ `residual_dT`, NOT `residual_dh`.  The grid is rebuilt at
            ## the current `T` on every residual evaluation, so every step
            ## scales and the partial `d/dh_n` is not the total -- for
            ## Gear-2 it is 3/2 of it on a uniform grid, measured against
            ## finite differences at 1.4859/1.4939/1.4972 for 100/200/400
            ## points, converging on the exact 3/2.  Euler and trapezoidal
            ## were never wrong: their coefficients depend on `h_n` alone,
            ## so the partial IS the total, which is why only Gear-2 was
            ## hit.  See `Integrator.companion_dT`.
            (self._dfdT,) = remove_row_col(
                (tr.residual_dT(x_full, dt),), irefnode, toolkit)
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
              maxiterations=20, grid=None, matrix_free=False,
              x0_unknown=False, tstab=None):
        """Solve for the periodic steady state.

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
        Q EITHER.  Spectre does offer one; this tree's user reports from
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
        """
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
        ## know about -- genuinely broken, and SpectreRF "outlaws [them]
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
        self._open_at_x0 = bool(x0_unknown)
        times, hs = self._period_grid(period, int(period / dt), grid)
        npts = len(times)
        self._grid_fracs = (None if grid is None
                            else np.asarray(grid, dtype=float))
        ## read by `_period_grid`, which is called from the residual
        ## closures and so cannot take it as an argument
        self._open_at_x0 = bool(x0_unknown)
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
            ## ⚠ THE RULE ITSELF IS CANONICAL.  Aprille & Trick's
            ## oscillator paper, Step 3: "select k by
            ## |f_k(x^i(T^i))| = max_k |f_k(x^i(T^i))|" -- argmax of the
            ## vector field, which is what an argmax over one step is up to
            ## `h`, and `h` is common to every coordinate so the mixed units
            ## cancel exactly.  Both the rule and the decision not to
            ## replace it with a Poincare row have precedent as well as
            ## measurement behind them.
            ##
            ## ⚠ WHAT IS NOT CANONICAL IS FREEZING IT.  Their Step 3 sits
            ## INSIDE the iteration (Step 5 returns to Step 1), so `k` and
            ## the pinned value are re-chosen from the CURRENT trajectory
            ## every iterate -- "note that in this method, an initial k and
            ## w_0k are not required".  Pinning `w_0k = x_0k^i`, the
            ## iterate's OWN current value, is attainable by construction,
            ## which removes the failure mode measured below (a far seed
            ## pinning a value the orbit never reaches) structurally rather
            ## than by advice.  Not done here: `analysis.fsolve` exposes no
            ## per-iteration hook, so it needs the autonomous outer solve
            ## restructured, and re-selecting BETWEEN outer iterations does
            ## not threaten `phi` being a function of `x_0` alone.
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
            ## from the current trajectory between outer iterations REGRESSES
            ## the working case.  Van der Pol at mu=1 from an ON-ORBIT seed
            ## went from converged to NOT converged, and far seeds wandered
            ## to periods of -52, -1088 and +110 against a true 6.6633.
            ##
            ## ⚠ AND THE REASON IS STRUCTURAL, not a bug in the attempt.
            ## Pinning the iterate's own value makes the phase residual
            ## `x_0[k] - pin` IDENTICALLY ZERO at every iterate, so the row
            ## carries no information: it constrains the STEP (`dz[k] = 0`)
            ## and nothing else, and with `k` re-chosen each iteration a
            ## different coordinate is frozen each time, so the orbit slides
            ## along itself.
            ##
            ## A&T do not have this problem because THEY HAVE NO PHASE
            ## EQUATION.  Their unknown vector SUBSTITUTES the period for the
            ## pinned coordinate -- `v = [x_01, ..., x_0(k-1), T,
            ## x_0(k+1), ..., x_0n]` -- an n x n system in which `x_0k` is a
            ## CONSTANT rather than an unknown with a trivially satisfied
            ## equation.  The constraint is structural where ours is
            ## algebraic.  So Step 3 is not portable to a bordered
            ## formulation as a drop-in: taking it means taking the
            ## substitution with it.
            ##
            ## The frozen pin's failure mode below therefore STANDS, and its
            ## fix is that substitution, not a moving pin.
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
        ## (`newton = True` lived here and fed a message branch that is
        ## gone: it was set unconditionally, so the 'successive
        ## substitution' alternative was dead and the claim it selected was
        ## false anyway -- see the non-convergence warning below.)

        def func(x):
            x0, x_end, Mx, _Mt = self._traverse(x, period, times, hs,
                                                want_dT=False,
                                                open_at_x0=x0_unknown)
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
                                               want_dT=True,
                                               open_at_x0=x0_unknown)

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
            _mfc = None
            if matrix_free:
                ## BOTH ENLARGEMENTS AT ONCE, matrix-free.  With
                ## `w = (v_0, v_{-1}, s)` and `M v` the `2m` pair map,
                ##     J w = [ v - M v - s (Pt_last, Pt_prev) ; v[k] ]
                ## -- one phase row, pinning the `x_0` block only, exactly as
                ## the dense system does and for the same reason.
                self._monodromy = None

                def _build_comp(z):
                    x0_, xm1_, T_ = z[:m_], z[m_:2 * m_], float(z[-1])
                    tms_, hsT_ = self._period_grid(T_, npts, self._grid_fracs)
                    (C0_, st_, xl_, xp_, Ptl_,
                     Ptp_) = self._traverse_factored(
                        x0_, xm1_, tms_, hsT_, T=T_, want_dT=True)
                    Ptv_ = np.concatenate((np.asarray(Ptl_).ravel(),
                                           np.asarray(Ptp_).ravel()))
                    F_ = np.concatenate(
                        (np.asarray(x0_, dtype=float) - np.asarray(xl_, dtype=float),
                         np.asarray(xm1_, dtype=float) - np.asarray(xp_, dtype=float),
                         [float(np.asarray(x0_, dtype=float)[phase_k])
                          - phase_pin]))

                    def mv_(w):
                        v_, s_ = w[:2 * m_], float(w[2 * m_])
                        top = (v_ - alpha * self._monodromy_matvec(C0_, st_, v_)
                               - s_ * Ptv_)
                        return np.concatenate((top, [v_[phase_k]]))
                    return F_, mv_

                def _mfc(z0_, ab_, xt_, rt_, mi_):
                    return self._matrix_free_newton(_build_comp, z0_, ab_,
                                                    xt_, rt_, mi_)
            z_ss, _info, _ier, _mesg = self._free_period_solve(
                func_autonomous_solved_history, z0, abstol_z, xtol_z,
                _shoot_reltol, maxiterations, period, solver=_mfc)
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
            _mf = None
            if matrix_free:
                ## THE BORDERED SYSTEM, matrix-free.  `dphi/dT` is ONE column
                ## and does not depend on the Krylov direction, so the
                ## trajectory pass computes it once per Newton iteration and
                ## the matvec just uses it:
                ##     J [v; s] = [ (I - M) v - s dphi/dT ; v_k ]
                m_ = n - 1

                def _build_auto(z):
                    xin_, T_ = z[:m_], float(z[-1])
                    tms_, hsT_ = self._period_grid(T_, npts, self._grid_fracs)
                    op_, st_, x0_, xe_, Mt_ = self._traverse_factored_plain(
                        xin_, T_, tms_, hsT_, want_dT=True,
                        open_at_x0=x0_unknown)
                    x0_ = np.asarray(x0_, dtype=float)
                    Mt_ = np.asarray(Mt_, dtype=float).ravel()
                    F_ = np.concatenate((x0_ - np.asarray(xe_, dtype=float),
                                         [x0_[phase_k] - phase_pin]))

                    def mv_(w):
                        v_, s_ = w[:m_], float(w[m_])
                        top = (v_ - alpha * self._monodromy_matvec_plain(
                            op_, st_, v_)) - s_ * Mt_
                        return np.concatenate((top, [v_[phase_k]]))
                    return F_, mv_

                def _mf(z0_, ab_, xt_, rt_, mi_):
                    return self._matrix_free_newton(_build_auto, z0_, ab_,
                                                    xt_, rt_, mi_)
            z_ss, _info, _ier, _mesg = self._free_period_solve(
                func_autonomous, z0, abstol_z, xtol_z, _shoot_reltol,
                maxiterations, period, solver=_mf)
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
                    toolkit=self.toolkit, full_output=True, line_search=True)
            x0_ss, xm1_ss = z_ss[:n - 1], z_ss[n - 1:]
        elif matrix_free:
            ## RECORDED SCOPE ITEM 6 on the PLAIN path: `m` columns rather
            ## than `2m`, so a lower share and a lower ceiling than the
            ## solved-history route -- 42.5% and 1.71x at m=502, 60.7% and
            ## 2.50x at m=1002.
            def _build(z):
                opening, steps, x0_, xe_, _dT = self._traverse_factored_plain(
                    z, period, times, hs, want_dT=False,
                    open_at_x0=x0_unknown)
                return (np.asarray(x0_) - np.asarray(xe_),
                        lambda v: v - alpha * self._monodromy_matvec_plain(
                            opening, steps, v))
            x0_ss, _info, _ier, _mesg = self._matrix_free_newton(
                _build, np.asarray(x, dtype=float), _tol, _tol,
                _shoot_reltol, maxiterations)
        else:
            x0_ss, _info, _ier, _mesg = analysis.fsolve(
                func, x, maxiter=maxiterations, reltol=_shoot_reltol,
                abstol=_tol, xtol=_tol, toolkit=self.toolkit,
                full_output=True, line_search=True)
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
            warnings.warn(
                'PSS: the shooting solve did not converge in %d iterations '
                '(method=%r). ⚠ The returned waveform IS STILL A FULL '
                'RESULT -- it is the last iterate, not a periodic steady '
                'state -- so a reader who does not check `converged` gets '
                'an array that looks like an answer and is not. Raise '
                "maxiterations, or use method='gear', whose solved-history "
                'formulation has an exact Jacobian and converges '
                'quadratically where the plain path is a contraction with a '
                'linear rate.'
                % (maxiterations, method),
                RuntimeWarning, stacklevel=2)
        
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
        ## inside `_matrix_free_solve`'s `build` closure, whose `steps` go
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
        self._want_lte = True
        lte_seen = []
        for t, dt in walk:
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
        
        fpss = analysis.CircuitResult(self.cir, x=FX, xdot=None,
                                      sweep_values=freqs, sweep_label='freq', 
                                      sweep_unit='Hz')
        
        return InternalResultDict({'tpss': tpss, 'fpss': fpss})

class PAC(Analysis):
    """Small-signal analysis over a periodic operating point, matrix-free.

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
    """

    parameters  = [Parameter(name='analysis', desc='Analysis name',
                             default='ac')]

    ## How hard GMRES is asked to solve the `m x m` system, relative to the
    ## PSS reltol that produced the operating point.  Looser than the
    ## operating point itself would be answering a question the trajectory
    ## cannot support; much tighter buys nothing, because the linearisation
    ## is only as good as the trajectory.
    KRYLOV_FACTOR = 1e-2

    def __init__(self, cir, toolkit=None, **kvargs):
        self.parameters = super(PAC, self).parameters + self.parameters
        super(PAC, self).__init__(cir, toolkit=toolkit, **kvargs)

    def solve(self, pss, freqs, refnode=gnd, recycle=True):
        """Sideband response at each frequency in `freqs`.

        `pss` must be a CONVERGED `PSS` -- the periodic operating point is
        what this linearises about, and there is no meaningful small-signal
        answer about a non-solution.  `PSS.factored_period()` enforces it.

        `recycle` shares one Krylov subspace across the sweep, which is
        where the sweep's cost goes; see `_solve_subspace`.
        """
        toolkit = self.toolkit
        freqs = np.atleast_1d(np.asarray(freqs, dtype=float))
        fp = pss.factored_period()
        T = float(fp.T)
        m = self.cir.n - 1

        irefnode = self.cir.get_node_index(refnode)
        if irefnode != pss.irefnode:
            raise ValueError(
                'PAC: refnode (index %d) differs from the PSS the operating '
                'point came from (index %d). The monodromy eliminated one '
                'and this would report against the other.'
                % (irefnode, pss.irefnode))
        ## ⚠ `analysis=` BY KEYWORD.  `Circuit.u(t, epar, analysis, ...)`
        ## takes `epar` second, and the withdrawn body wrote
        ## `self.cir.u(0, analysis_name)` -- passing 'ac' as the element
        ## parameter set and taking the TRANSIENT source vector, which is
        ## zero at `t = 0` for every sinusoid.  The whole analysis would
        ## have returned zeros, silently, with no source to speak of.
        (u_ac,) = remove_row_col((self.cir.u(0, analysis=self.par.analysis),),
                                 irefnode, toolkit)
        if not np.any(np.asarray(u_ac)):
            raise ValueError(
                'PAC: the %r source vector is identically zero, so there is '
                'nothing to analyse. Independent sources take their '
                'small-signal amplitude from `vac`/`iac`, not from `va` -- '
                'a source with va= set and vac=0 drives the operating point '
                'and not this.' % self.par.analysis)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()

        ## the forced response at each frequency -- one period replay each,
        ## and unavoidable: the source is what changes across the sweep
        ## ⚠ THE MANUFACTURING STEP IS NOT IN `steps`, AND IT COSTS AN
        ## ORDER.  On the plain path `_traverse_factored_plain` takes one
        ## step OUTSIDE the loop to manufacture a history, and folds its
        ## effect into the `opening` triple as a flat-history assumption.
        ## For the HOMOGENEOUS map that is the documented approximation the
        ## whole plain path is built on.  For the DRIVEN one it also means
        ## the source is never applied at that step -- one step of `u` out
        ## of `N`, i.e. a relative O(h).
        ##
        ## ⚠ MEASURED, on the Q=20 resonator against the AC analysis at
        ## 700 Hz, rel error per doubling of the grid:
        ##
        ##     trap, plain            2.00x  (O(h))   4.13e-03 at 250 pts
        ##     trap, x0_unknown=True  4.00x  (O(h^2)) 1.09e-04 at 250 pts
        ##     euler, either          2.00x  (O(h))   1.40e-02, unchanged
        ##
        ## The euler row is the control: `x0_unknown` does not move it at
        ## all (identical to five digits), so the trapezoidal gain is the
        ## manufacturing step and not something else the formulation does.
        ## The trajectory is NOT the problem -- trap's waveform converges at
        ## 4.2x per doubling either way.
        ##
        ## So this is a silent order loss for a caller who did nothing
        ## wrong, which is the one thing worth a warning.  Gear-2 takes the
        ## solved-history path and has no manufacturing step at all.
        if (fp.kind == 'plain' and not fp.open_at_x0
                and pss.par.method != 'euler'):
            warnings.warn(
                'PAC: this operating point was solved on the PLAIN path '
                'with a manufacturing step (method=%r, x0_unknown=False). '
                'The manufacturing step carries no small-signal source, so '
                'the response is FIRST order in the timestep whatever the '
                "method's own order -- measured 2.00x per doubling against "
                '4.00x for the same run with x0_unknown=True. The answer is '
                'not wrong, it is one order less accurate than the '
                'trajectory it came from. Re-solve with x0_unknown=True, or '
                "with method='gear', to get the method's own order."
                % pss.par.method,
                RuntimeWarning, stacklevel=2)

        self._check_circuit(pss)
        for f in freqs:
            self._check_harmonic(pss, f, 'a sweep point')

        rhs = []
        for f in freqs:
            w, _ = pss._forced_replay(fp, f, u_ac)
            rhs.append(np.exp(-2j * np.pi * f * T) * np.asarray(w))

        alphas = [np.exp(-2j * np.pi * f * T) for f in freqs]
        tol = max(pss.par.reltol * self.KRYLOV_FACTOR, 1e-14)
        if recycle:
            ys, self.matvecs = self._solve_subspace(fp, alphas, rhs, tol)
        else:
            ys, self.matvecs = self._solve_each(fp, alphas, rhs, tol)

        ## one driven replay per frequency turns `y_0` into the period
        outfreq, outV = [], []
        for f, y0 in zip(freqs, ys):
            _end, ysteps = pss._forced_replay(fp, f, u_ac, y0=y0, collect=True)
            y = np.array([np.asarray(y0)[:m]] + [np.asarray(v)[:m]
                                                 for v in ysteps])
            ## `v(t) = y(t) exp(-j w t)` is T-periodic; its DFT is the
            ## sideband set, exactly as the withdrawn body intended
            tms = np.asarray(fp.times, dtype=float)[:len(y)]
            v = y * np.exp(-2j * np.pi * f * tms)[:, None]
            sb, V = freq_analysis(v, tms, axis=0)
            outfreq.extend((np.abs(sb + f)).tolist())
            outV.extend(V.tolist())

        order = np.argsort(np.asarray(outfreq))
        fout = np.asarray(outfreq)[order]
        X = np.asarray(outV)[order]
        X = np.concatenate((X[:, :irefnode],
                            np.zeros((len(fout), 1)),
                            X[:, irefnode:]), axis=1)
        self.result = analysis.CircuitResult(
            self.cir, x=X.T, xdot=None, sweep_values=fout,
            sweep_label='freq', sweep_unit='Hz')
        return self.result

    def adjoint_transfer_row(self, pss, freq, output, recycle_tol=None):
        """Every source to ONE output, in a single transposed solve.

        Returns a row `r` of length `m`: `r[i]` is the small-signal
        response at `output` (at `t = 0`) to a unit source injected at
        reduced coordinate `i` at `freq`. Forward, that is `m` separate
        solves; here it is one.

        ⚠ THIS IS THE ASYMMETRY pnoise IS SHAPED BY, and the reason
        Okumura et al. (1993) reach for the adjoint at all: "it is
        efficient to use the adjoint method ... BECAUSE CIRCUITS HAVE MANY
        NOISE SOURCES." Recycling does not help the forward route, because
        the right-hand side is what changes from source to source.

            output = d^T y_0 = alpha * d^T (I - alpha M)^-1 w(u)
                             = alpha * ((I - alpha M)^-T d)^T W u

        so one transposed solve for `x^a`, then `W^T x^a` from the reverse
        replay, and the whole row falls out. MEASURED against `m` forward
        solves on an RC ladder: agreement 9.6e-16.

        ⚠ WHAT THIS IS NOT, so the next reader does not over-read it. The
        output here is the state at `t = 0`, a single linear functional.
        A SIDEBAND coefficient `H_l` is a functional DISTRIBUTED over the
        period -- `(1/N) sum_n exp(-j l w0 t_n) d^T y_n` -- and its adjoint
        needs the reverse pass to take an injection at every step rather
        than a seed at the end. That extension is the next piece of A3, and
        it is not built.

        ⚠ SOLVED-HISTORY ONLY (`method='gear'`), because the reverse replay
        is; a one-step method needs its own reverse recursion.
        """
        import scipy.sparse.linalg as spla
        fp = pss.factored_period()
        if fp.kind != 'solved_history':
            raise NotImplementedError(
                'PAC: the adjoint row needs the transposed replay, which is '
                "implemented for the solved-history map only. Re-solve with "
                "method='gear'. (PSS used %r.)" % pss.par.method)
        self._check_circuit(pss)
        self._check_harmonic(pss, freq, 'the adjoint row')
        m = pss.cir.n - 1
        n = fp.width
        alpha = np.exp(-2j * np.pi * float(freq) * float(fp.T))

        d = np.zeros(n, dtype=complex)
        if np.isscalar(output):
            d[int(output)] = 1.0
        else:
            out = np.asarray(output, dtype=complex).ravel()
            d[:len(out)] = out

        count = [0]

        def _mv(v):
            count[0] += 1
            return np.asarray(v) - alpha * fp.matvec_transposed(v)

        tol = (self.KRYLOV_FACTOR * pss.par.reltol if recycle_tol is None
               else recycle_tol)
        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        xa = self._gmres_checked(A, d, max(tol, 1e-14), 'the adjoint solve')
        self.matvecs = count[0]
        return alpha * pss._forced_replay_transposed(fp, freq, xa)

    def adjoint_sideband_row(self, pss, freq, output, sidebands=0):
        """`H_l` rows: every source to ONE output's sideband `l`.

        Returns an array of shape `(len(sidebands), m)`. Entry `[li, i]` is
        the coefficient at sideband `l` of the output at `output`, for a
        unit source injected at reduced coordinate `i` at `freq`:

            H_l = (1/N) sum_n exp(-j l w0 t_n) d^T y_n

        ⚠ THIS IS THE ONE `adjoint_transfer_row` IS NOT.  That row is the
        response at a single instant -- one linear functional, adjointed by
        seeding the reverse pass at the end.  A SIDEBAND is a functional
        DISTRIBUTED over the period, so its adjoint takes an injection at
        EVERY step, and the answer comes in two pieces:

            dH/du  =  [the forced part, from the injected reverse pass]
                    + [alpha * W^T z, with z = (I - alpha M)^-T g]

        where `g` is the reverse pass's own final state -- the sensitivity
        of the functional to the initial state `y_0`, which is itself a
        function of the source through the periodic boundary condition.
        Dropping the second term would leave an answer that looks entirely
        reasonable: MEASURED on an RC ladder the two terms are comparable
        in size (303 against 498 at `l = 0`, 79 against 606 at `l = 1`), so
        neither is a correction to the other.

        Still ONE transposed solve per sideband whatever the number of
        sources, which is the property pnoise needs.  Agreement with the
        `m` forward driven solves: 9.2e-16 / 3.3e-16 / 1.3e-15 at
        `l = 0 / 1 / -2`.

        ⚠ SOLVED-HISTORY ONLY (`method='gear'`), like the reverse pass.
        """
        import scipy.sparse.linalg as spla
        fp = pss.factored_period()
        if fp.kind != 'solved_history':
            raise NotImplementedError(
                'PAC: the sideband row needs the transposed replay, which '
                "is implemented for the solved-history map only. Re-solve "
                "with method='gear'. (PSS used %r.)" % pss.par.method)

        self._check_circuit(pss)
        self._check_harmonic(pss, freq, 'the sideband row')
        m = pss.cir.n - 1
        n = fp.width
        T = float(fp.T)
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        w0 = 2.0 * np.pi / T
        alpha = np.exp(-2j * np.pi * float(freq) * T)
        ls = np.atleast_1d(np.asarray(sidebands, dtype=int))

        ## ⚠ A HARD BOUND, NOT A HEURISTIC.  Okumura et al. eq. (32): the
        ## maximum frequency the analysis can speak about is the grid's own
        ## `w_max`, so `|l| <= (w_max - w0)/ws`.  You cannot alias down from
        ## above what the grid can represent, and a ratio test on the
        ## accumulated power operates INSIDE this ceiling rather than
        ## instead of it -- an implementation carrying only the ratio test
        ## terminates for the wrong reason.  The grid's ceiling here is its
        ## Nyquist, `N/2` harmonics of the period.
        lmax = N // 2
        bad = ls[np.abs(ls) > lmax]
        if len(bad):
            raise ValueError(
                'PAC: sideband %s is above the grid\'s Nyquist (|l| <= %d '
                'at %d points per period). Nothing can alias down from '
                'above the maximum frequency the grid represents, so this '
                'is not a tolerance to relax -- use a finer period grid.'
                % (bad.tolist(), lmax, N))

        d = np.zeros(m, dtype=complex)
        if np.isscalar(output):
            d[int(output)] = 1.0
        else:
            out = np.asarray(output, dtype=complex).ravel()
            d[:len(out)] = out

        count = [0]

        def _mv(v):
            count[0] += 1
            return np.asarray(v) - alpha * fp.matvec_transposed(v)

        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        tol = max(self.KRYLOV_FACTOR * pss.par.reltol, 1e-14)
        phase = np.exp(2j * np.pi * float(freq) * tms[1:N + 1])

        rows = np.zeros((len(ls), m), dtype=complex)
        for li, l in enumerate(ls):
            ## ⚠ THE PHASE OF THE INPUT COMES OUT FIRST, and getting this
            ## wrong is self-consistent rather than loud.  What is
            ## T-PERIODIC is `v(t) = y(t) exp(-j w t)`, not `y` -- so the
            ## sideband set is the DFT of `v`, which is what `solve` takes.
            ## Decomposing `y` instead gives a Dirichlet kernel smeared
            ## across every `l` whenever `f` is not a multiple of `1/T`,
            ## and it AGREES with a forward reference written the same way,
            ## so only a check against a circuit whose answer is known
            ## independently catches it.
            inject = ((np.exp(-1j * (float(l) * w0 + 2.0 * np.pi
                                     * float(freq)) * tms[:N]) / N)[:, None]
                      * d[None, :])
            g, ts, _st = pss._monodromy_matvec_transposed(
                fp.opening, fp.steps, np.zeros(n, dtype=complex),
                collect=True, inject=inject)
            forced = -np.tensordot(phase, np.asarray(ts), axes=(0, 0))
            ## ⚠ ON AN OSCILLATOR THIS OPERATOR IS SINGULAR AT EVERY
            ## HARMONIC and near-singular around them, which is exactly
            ## where phase noise is measured.  The deflated route borders
            ## the pole out and carries `1/(1 - alpha)` analytically; on a
            ## driven circuit there is no pole and the plain solve is both
            ## correct and cheaper.
            if getattr(pss, 'autonomous', False):
                z = self._deflated_solve(pss, alpha, g, transposed=True,
                                         tol=tol)
            else:
                z = self._gmres_checked(
                    A, g, tol, 'the adjoint solve at sideband %d' % l)
            rows[li] = forced + alpha * pss._forced_replay_transposed(
                fp, freq, z)
        self.matvecs = count[0]
        return rows

    ## How small a sideband's contribution must be, relative to the running
    ## total, before the accumulation stops.  Okumura et al.: powers are
    ## "accumulated until their contributions become negligible".
    ALIAS_RATIO_TOL = 1e-9

    def pnoise(self, pss, freq, output, ratio_tol=None, maxsidebands=None):
        """TIME-AVERAGED output noise PSD at `freq`, sidebands folded in.

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
        singularity at EVERY harmonic.  SpectreRF: "place a cluster of
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
        ⚠ THE ACTUAL BARRIER IS COST, which is a different decision.  The
        source count is (timepoints per period) x (noisy devices) -- a
        500-point grid with 50 noisy devices is 25 000 stationary sources
        -- and the reported noise analysis ran at ~14x the PSS PER
        FREQUENCY POINT, "because all aliasing components need to be
        computed".  Whether that is affordable is unmeasured here.
        ⚠ AND ITS AUTHORS LEFT THE PHYSICS OPEN: "it is further necessary
        to discuss the correspondence between the actual physical phenomena
        of noises and this modeling".  The windowed-stationary
        decomposition is a numerical construct, and its fidelity to a real
        device is not settled by its numerical validation.
        (Cited, not verified here; the boxcar observation is ours.)

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
        `PAC.solve`, WHICH IS THE OPPOSITE OF HOW IT READS.  SpectreRF's
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
        definitions -- which points at the right layer immediately.  A
        published removal exists in a time-domain form written for
        shooting, using that same null vector; it is not built.  (Cited,
        not verified here.)

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
        SpectreRF Theory on PNoise and QPnoise, both: "THE TIME-AVERAGE of
        the noise at the output of the circuit is computed in the form of a
        spectral density versus frequency."  Same quantity, same
        limitation.  (Relayed from the docs session; cited, not verified
        here.)
        """
        self._check_circuit(pss)
        fp = pss.factored_period()
        m = pss.cir.n - 1
        N = len(fp.steps)
        T = float(fp.T)
        f0 = 1.0 / T
        tol = self.ALIAS_RATIO_TOL if ratio_tol is None else float(ratio_tol)
        lmax = N // 2 if maxsidebands is None else min(int(maxsidebands),
                                                       N // 2)

        w = 2.0 * np.pi * float(freq)
        cy = self._cy_reduced(pss, w)

        total = 0.0
        used = []
        quiet = 0
        self.alias_stop = 'bound'
        for l in range(0, lmax + 1):
            step = 0.0
            for sl in ((0,) if l == 0 else (l, -l)):
                fin = float(freq) - sl * f0
                h = self.adjoint_sideband_row(pss, fin, output, sl)[0]
                step += float(np.real(h @ self._cy_reduced(
                    pss, 2.0 * np.pi * fin) @ np.conj(h)))
                used.append(sl)
            total += step
            if total > 0 and abs(step) < tol * abs(total):
                ## ⚠ TWO QUIET PAIRS, NOT ONE.  A single sideband can come
                ## back near zero by symmetry while its neighbours do not,
                ## and stopping there would truncate a series that had not
                ## converged.
                quiet += 1
                if quiet >= 2:
                    self.alias_stop = 'ratio'
                    break
            else:
                quiet = 0
        self.sidebands_used = used
        ## ⚠ WHICH RULE STOPPED IT IS PART OF THE ANSWER.  Ending on the
        ## ratio test means the series converged; ending on the Nyquist
        ## bound means the grid ran out before the series did, and the
        ## number is a LOWER bound on the folded noise -- every sideband
        ## above the grid's own maximum frequency is missing, not small.
        ## A strongly switching circuit does this readily: measured on a
        ## driven diode at 80 points per period, the accumulation reached
        ## l = +-39 without the ratio test ever firing, while folding was
        ## already contributing 62% of the total.
        if self.alias_stop == 'bound' and lmax > 0:
            warnings.warn(
                'PAC.pnoise: the sideband accumulation stopped at the '
                "grid's Nyquist (|l| = %d at %d points per period), not "
                'because the contributions became negligible. Sidebands '
                'above the grid\'s maximum frequency are MISSING rather '
                'than small, so this is a lower bound on the folded noise. '
                'Re-solve the PSS on a finer period grid and compare.'
                % (lmax, N),
                RuntimeWarning, stacklevel=2)
        return total, used

    def _cy_reduced(self, pss, w):
        """`CY` with the reference node removed, refusing a moving one.

        ⚠ THE CHECK IS THE POINT.  A bias-dependent `CY` makes the sources
        CYCLOSTATIONARY, and then sidebands stop adding in power -- the
        stationary sum this class computes would be the wrong model, not
        merely an inaccurate one, and nothing downstream would say so.
        Sampled at three states on the converged orbit rather than argued
        from the element types, because a compact model's `CY` reads `x`
        and the discrete library's does not.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        states = [np.zeros(pss.cir.n - 1),
                  np.asarray(fp.x_last, dtype=float).ravel(),
                  np.asarray(fp.x_prev, dtype=float).ravel()[:pss.cir.n - 1]]
        mats = []
        for xr in states:
            xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cy = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cy,) = remove_row_col((cy,), irn, pss.toolkit)
            mats.append(np.asarray(cy, dtype=complex))
        scale = max(float(np.max(np.abs(mats[0]))), 1e-300)
        for other in mats[1:]:
            drift = float(np.max(np.abs(other - mats[0]))) / scale
            if drift > 1e-9:
                raise NotImplementedError(
                    'PAC.pnoise: this circuit has a BIAS-DEPENDENT CY '
                    '(varies by %.3g over the orbit), so its noise sources '
                    'are cyclostationary. The sidebands are then correlated '
                    'through the window Fourier coefficients and no longer '
                    'add in power, so the stationary sum here would be the '
                    "wrong model rather than an imprecise one. Okumura's "
                    'section III-B has the construction; it is not built.'
                    % drift)
        return mats[0]

    ## How close to a harmonic of `f0` counts as "on" it, as a fraction of
    ## `f0`.
    ##
    ## ⚠ TIGHTENED BY A7.  This used to be the conditioning floor being
    ## accepted, because `sigma_min` of the plain operator falls LINEARLY
    ## with the distance and everything nearer was unusable.  The deflated
    ## solve removes that: its conditioning is FLAT (measured 2.04e-01 from
    ## 0.3 down to 1e-9 of `f0`), so the only remaining reason to refuse is
    ## the physical one -- at an EXACT harmonic `1/(1 - alpha)` is a
    ## division by zero and the response is genuinely unbounded.  So the
    ## guard now excludes only what has no finite answer, not what was
    ## merely hard to compute.
    HARMONIC_GUARD = 1e-12

    def _check_circuit(self, pss):
        """The operating point must belong to THIS circuit, not a similar one.

        ⚠ A DRIVEN OSCILLATOR MAKES THIS A CORRECTNESS TRAP RATHER THAN A
        TYPO GUARD.  The natural way to model one is to solve the PSS of
        the bare oscillator and then treat the injection as a perturbation
        -- and it is wrong, because the injection DEVICE is present even
        when its SIGNAL is zero.  Buonomo & Lo Schiavo: "in absence of the
        injection signal, the injection circuit affects the basic LC
        oscillator by CHANGING THE NONLINEARITY OF THE FEEDBACK LOOP ...
        [it] can affect the start-up condition of the basic differential LC
        oscillator OR ITS OSCILLATION AMPLITUDE, or both."

        So the free-running orbit of the circuit-with-the-device is not the
        orbit of the circuit-without-it, and every Floquet quantity built
        on the wrong one inherits the error -- monodromy, PPV, phase noise.
        The analysis would converge and report a plausible number.

        The reference-node check below catches a mismatched `refnode`; it
        cannot catch this, because two circuits differing by one device
        have the same reference node and often the same node count.
        """
        if self.cir is not pss.cir:
            raise ValueError(
                'PAC: this analysis was built on a different circuit object '
                'than the PSS it was handed. If that is deliberate -- e.g. '
                'solving the PSS of a bare oscillator and perturbing a '
                'version with an injection device added -- it is a '
                'CORRECTNESS error, not a bookkeeping one: the injection '
                'device changes the free-running orbit even with its signal '
                'at zero, so the base solution is the wrong one to '
                'linearise about. Solve the PSS on the SAME circuit.')

    def _check_harmonic(self, pss, freq, what):
        """Refuse an autonomous small-signal solve sitting on a harmonic.

        ⚠ `I - exp(-j w T) M` IS SINGULAR AT EVERY HARMONIC OF `f0`, NOT
        JUST AT DC, AND ONLY FOR AN OSCILLATOR.  At `w = k w0` the factor
        `exp(-j w T)` is 1 and the operator is `I - M`, which an autonomous
        circuit's unit multiplier makes singular.  MEASURED on van der Pol,
        `sigma_min(I - exp(-jwT) M)`:

            offset/f0    0     0.25    0.5    0.75    1      2      3
            sigma_min   2.8e-11 0.51   0.65   0.51  2.8e-11 2.8e-11 2.8e-11

        and LINEAR in the distance to the nearest one -- 2.5e-1, 2.6e-2,
        2.6e-3, 2.6e-4 at 0.9, 0.99, 0.999, 0.9999 of the way there.  The
        same sweep on a DRIVEN ladder never drops below 0.78: no unit
        multiplier, no singularity, harmonics included.

        ⚠ AND IT IS PHYSICS, NOT CONDITIONING.  A perturbation at a
        harmonic is a perturbation along the orbit, and an oscillator's
        response to that is unbounded phase drift -- there is no bounded
        periodic answer to return.  So this refuses rather than tightening
        a tolerance, and says which quantity to ask for instead.
        """
        if not getattr(pss, 'autonomous', False):
            return
        f0 = 1.0 / float(pss.factored_period().T)
        r = abs(float(freq)) / f0
        d = abs(r - round(r))
        if d <= self.HARMONIC_GUARD:
            raise ValueError(
                'PAC: %s at %.6g Hz is on harmonic %d of this OSCILLATOR\'s '
                'own frequency (%.6g Hz), where I - exp(-j w T) M is '
                'singular -- the unit Floquet multiplier makes it exactly '
                'I - M there. That is physics, not conditioning: a '
                'perturbation along the orbit produces unbounded phase '
                'drift, so there is no bounded periodic response to '
                'return. Ask off-harmonic, or ask for the phase quantity '
                'instead (PSS.ppv()).' % (what, float(freq), round(r), f0))

    @staticmethod
    def _gmres_checked(A, b, rtol, what):
        """GMRES, judged by its RESIDUAL rather than by its status flag.

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
        """
        import scipy.sparse.linalg as spla
        n = b.shape[0]
        x, info = spla.gmres(A, b, rtol=rtol, restart=min(n, 50),
                             maxiter=min(n, 200))
        r = float(np.linalg.norm(b - A.matvec(x)))
        scale = max(float(np.linalg.norm(b)), 1e-300)
        if r / scale > max(1e3 * rtol, 1e-8):
            raise RuntimeError(
                'PAC: %s did not converge (info=%r, relative residual '
                '%.3e). Near a harmonic of the oscillator this operator is '
                'genuinely near-singular and a smaller tolerance will not '
                'help -- move the offset, or ask for the phase quantity.'
                % (what, info, r / scale))
        return x

    def covariance(self, pss, samples=False):
        """The periodic (cyclostationary) state covariance — DRIVEN circuits.

        Returns `K0`, the covariance at `t = 0`; with `samples=True`,
        `(K0, [K_j])`, the covariance at every step, which is the
        time-varying statistic this exists to produce.

        The noise covariance obeys a Lyapunov recursion alongside the
        trajectory, `K_{j+1} = A_j K_j A_jᵀ + Q_j`, so over one period
        `K_N = M K_0 Mᵀ + K_1`.  Periodicity closes it:

            (I - M ⊗ M) vec(K_0) = vec(K_1)

        ⚠ ONE LINEAR SOLVE, NO NEWTON.  The Lyapunov equation is LINEAR in
        `K`, so shooting on it is exact in a single step — unlike the
        trajectory it rides on.  The monodromy of the covariance system is
        the KRONECKER SQUARE of the circuit's, so its multipliers are the
        pairwise products `lambda_i lambda_j`.

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
        """
        self._check_circuit(pss)
        if getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.covariance: an OSCILLATOR has no periodic covariance. '
                'Its unit multiplier squares to one, so I - M kron M is '
                'singular and the covariance grows without bound rather '
                'than settling -- that growth IS the phase diffusion, and '
                'its output noise is stationary rather than '
                'cyclostationary. Use oscillator_spectrum().')
        fp = pss.factored_period()
        if fp.kind != 'solved_history':
            raise NotImplementedError(
                'PAC.covariance: the per-step maps are rebuilt from the '
                "solved-history factors. Re-solve with method='gear'.")
        m = pss.cir.n - 1
        n = fp.width
        hs = np.diff(np.asarray(fp.times, dtype=float))
        cy = np.real(self._cy_reduced(pss, 2.0 * np.pi / float(fp.T)))

        ## the C ring as the forward recursion sees it -- see the replays
        cs0, cs1, ring = [], [], list(fp.opening)
        for _lu, C_new, _a, _b in fp.steps:
            cs0.append(ring[0])
            cs1.append(ring[1])
            ring = [C_new, ring[0]]

        def step_map(k):
            lu, _Cn, alphas, b = fp.steps[k]
            if b:
                raise NotImplementedError(
                    'PAC.covariance: derived for a b = 0 companion (Gear-2).')
            A = np.zeros((n, n))
            for j in range(n):
                p0 = np.zeros(m)
                p1 = np.zeros(m)
                (p0 if j < m else p1)[j if j < m else j - m] = 1.0
                A[:m, j] = -lu.solve(alphas[1] * (cs0[k] @ p0)
                                     + alphas[2] * (cs1[k] @ p1))
                A[m:, j] = p0
            return A

        As, Qs = [], []
        for k, (lu, _Cn, _a, _b) in enumerate(fp.steps):
            ## Q = Jf^-1 (CY / 2h) Jf^-T, symmetrised against round-off
            half = cy / (2.0 * hs[k])
            left = np.column_stack([lu.solve(half[:, j]) for j in range(m)])
            Q1 = np.column_stack([lu.solve(left[j, :]) for j in range(m)]).T
            Q = np.zeros((n, n))
            Q[:m, :m] = 0.5 * (Q1 + Q1.T)
            Qs.append(Q)
            As.append(step_map(k))

        K = np.zeros((n, n))
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
        K1 = K

        M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        S = np.eye(n * n) - np.kron(M, M)
        K0 = np.linalg.solve(S, K1.reshape(-1)).reshape(n, n)
        K0 = 0.5 * (K0 + K0.T)
        if not samples:
            return K0
        seq, K = [K0], K0
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
            seq.append(0.5 * (K + K.T))
        return K0, seq

    def diffusion_constant(self, pss):
        """`c` — the phase diffusion constant, in seconds.

        `c = (1/T) ∫ v₁ᵀ(t) B(t) Bᵀ(t) v₁(t) dt` with `B Bᵀ = CY`, so this
        is the time-average of a QUADRATIC form in the PPV.  It is the one
        scalar the whole free-running phase-noise spectrum is built from,
        and it reads, for a designer, as JITTER PER SECOND.

        ⚠ QUADRATIC FOR WHITE SOURCES, LINEAR FOR COLOURED ONES, and the
        two are different functionals of the same vector: a coloured
        source contributes `V_0m = (1/T) ∫ v₁ᵀ B_cm dt`, with no square.
        Using this one for a coloured source returns a plausible non-zero
        number from the same PPV.  Only stationary white sources are
        supported here, which `_cy_reduced` enforces.

        ⚠ ALREADY VALIDATED AGAINST A PHYSICAL MEASUREMENT.  The same
        expression, evaluated for a single injected source, was compared
        against the phase diffusion of the full nonlinear circuit measured
        by Monte Carlo with zero-crossing timing: ratio 0.9965 on the
        control.  So the scalar this spectrum is built on is not a fresh
        claim.
        """
        self._check_circuit(pss)
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.diffusion_constant: phase diffusion is a property of a '
                "FREE-RUNNING oscillator. A driven circuit's phase is its "
                "source's, and its noise is pnoise's problem, not this one.")
        v, info = pss.ppv()
        m = pss.cir.n - 1
        S = np.asarray(info['samples'])[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        h = np.diff(tms)
        T = float(pss.period)
        cy = self._cy_reduced(pss, 2.0 * np.pi / T)
        quad = np.einsum('ij,jk,ik->i', S, np.real(cy), S)
        return float((quad * h).sum() / T)

    @staticmethod
    def lorentzian(offsets, c, f0, harmonic=1):
        """The `i`-th harmonic's normalised lineshape at `offsets` from it.

            S_i(f) = i² f₀² c / (π² i⁴ f₀⁴ c² + f²)

        ⚠ EXACT FOR WHITE SOURCES, not a limiting form.  With coloured
        sources the transform "does not have a simple closed form" and only
        two-regime approximations exist — which is one more reason this
        module supports white sources only.

        ⚠ AND ITS TOTAL POWER IS EXACTLY 1.  `∫ a/(b²+f²) df = aπ/b`, and
        here `a = i² f₀² c`, `b = π i⁴ f₀⁴ c² ^ ½`… concretely `b = π i²
        f₀² c`, so the integral is exactly one.  **The carrier's power is
        redistributed, never created or destroyed** — which is the
        invariant that separates this from LTV small-signal treatments,
        which "erroneously predict infinite noise power [at the carrier] as
        well as infinite total integrated power".  It is asserted in the
        suite.

        The half-width is `π i² f₀² c` and the peak `1/(π² i² f₀² c)`, so a
        higher harmonic has a skirt scaling as `i²` and a corner as `i⁴` —
        `20 log₁₀(i)` dB noisier far out.
        """
        i = int(harmonic)
        if i == 0:
            return np.zeros_like(np.asarray(offsets, dtype=float))
        f = np.asarray(offsets, dtype=float)
        a = (i * i) * f0 * f0 * c
        b = np.pi * (i * i) * f0 * f0 * c
        return a / (b * b + f * f)

    def oscillator_spectrum(self, pss, offsets, output, harmonic=1):
        """Free-running output spectrum at `offsets` from harmonic `harmonic`.

        Returns `(S_v, L_dBc)`.  `S_v` is the one-sided PSD of the output
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

        ⚠ AND IT IS THE ONLY ROUTE THAT IS VALID BELOW THE CORNER.  A
        small-signal analysis cannot produce `L(f)` there however well
        conditioned it is: the excess phase is a Wiener process, its
        spectrum has a singularity at the origin and no physical meaning,
        and the finite value `L` attains comes from the NONLINEAR
        phase-to-voltage map — which is what this closed form carries.
        Reporting `S_phi` near the carrier instead is the mistake that
        object invites.
        """
        c = self.diffusion_constant(pss)
        f0 = 1.0 / float(pss.period)
        X = self.carrier_phasor(pss, output, harmonic)
        Sv = abs(X) ** 2 * self.lorentzian(offsets, c, f0, harmonic)
        with np.errstate(divide='ignore'):
            L = 10.0 * np.log10(np.maximum(Sv / max(abs(X) ** 2, 1e-300),
                                           1e-300))
        return Sv, L

    @staticmethod
    def am_pm_indices(a, b):
        """Split a sideband pair into AM and PM modulation indices.

        `a` and `b` are the upper and lower sideband amplitudes, each
        already divided by the carrier phasor.  Returns `(m_am, m_pm)`.

        THE WHOLE THING IS ONE CONJUGATE.  Write the complex envelope's
        deviation as `a e^{j w_m t} + b e^{-j w_m t}`.  The two sidebands
        COUNTER-ROTATE about the carrier phasor, so the sum traces an
        ellipse; the component ALONG the carrier is amplitude modulation
        and the component PERPENDICULAR to it is phase modulation.

          - pure AM keeps the envelope on the carrier's axis, which forces
            `d = conj(d)` for all `t`, i.e. `a = conj(b)`;
          - pure PM keeps it perpendicular, `d = -conj(d)`, i.e.
            `a = -conj(b)`.

        so `m_am = a + conj(b)` and `m_pm = a - conj(b)` -- each vanishing
        exactly when the other case holds.  No new solve: this is a change
        of basis on transfer functions `adjoint_sideband_row` already
        returns.

        ⚠ `conj(b)`, NOT `b`.  Using `a +- b` looks equally plausible and
        is wrong for any modulation whose sidebands are not real relative
        to the carrier -- it would report a rotating ellipse as pure AM.
        The conjugate is what makes the lower sideband counter-rotate.
        """
        a = np.asarray(a, dtype=complex)
        b = np.asarray(b, dtype=complex)
        return a + np.conj(b), a - np.conj(b)

    def carrier_phasor(self, pss, output, carrier=1):
        """The `carrier`-th Fourier coefficient of the steady-state output.

        Computed here rather than taken from `fpss`, whose spectrum is RMS
        and energy-folded -- correct for reporting a magnitude and useless
        for a phasor, since folding discards the phase the AM/PM split is
        made of.
        """
        if getattr(pss, 'waveform', None) is None:
            raise RuntimeError(
                'PAC: the PSS has no stored waveform -- call solve() first.')
        times, X = pss.waveform
        irn = pss.irefnode
        k = int(output)
        row = np.asarray(X, dtype=float)[k if k < irn else k + 1]
        t = np.asarray(times, dtype=float)[:-1]
        v = row[:len(t)]
        w0 = 2.0 * np.pi / float(pss.period)
        return complex(np.sum(v * np.exp(-1j * carrier * w0 * t)) / len(t))

    def am_pm(self, pss, freq, output, carrier=1):
        """AM and PM modulation indices at `carrier`, per noise/signal source.

        Returns `(m_am, m_pm)`, each a row of length `m`: the modulation a
        unit source at reduced coordinate `i`, driven at `freq`, imposes on
        the `carrier`-th harmonic of the output.

        ⚠ TWO SOLVES AT ±freq, NOT ONE.  The upper sideband of harmonic `i`
        sits at `i f0 + freq` and the lower at `i f0 - freq`; with the
        convention that an input at `f` produces output at `f + l f0`,
        those are `H_i(freq)` and `H_i(-freq)`.  They are NOT conjugates of
        each other -- that would hold for an LTI circuit, and the whole
        point of an LPTV analysis is that it does not.  Taking one and
        conjugating it would silently force `m_pm = 0` or `m_am = 0`
        depending on which.

        ⚠ AND AN OSCILLATOR IS ALMOST PURE PM NEAR ITS CARRIER, which is
        the physical check: the phase response to a perturbation goes as
        `1/w_m` while the amplitude response stays bounded, so
        `|m_pm|/|m_am|` grows without bound as `freq -> 0`.  A
        decomposition that got the conjugate wrong gives a bounded ratio
        instead.
        """
        C = self.carrier_phasor(pss, output, carrier)
        ## ⚠ RELATIVE TO THE SIGNAL, NOT AGAINST ZERO.  A harmonic the
        ## circuit does not produce still has a phasor of ~1e-16 rather
        ## than exactly 0, and dividing by it turns "there is no carrier
        ## here" into an enormous, confident modulation index.
        times, X = pss.waveform
        irn = pss.irefnode
        k = int(output)
        row = np.asarray(X, dtype=float)[k if k < irn else k + 1]
        scale = float(np.max(np.abs(row)))
        if abs(C) <= 1e-9 * max(scale, 1e-300):
            raise ValueError(
                'PAC.am_pm: the output carries no component at harmonic %d '
                '(|C| = %.3e against a signal scale of %.3e), so there is '
                'no carrier to modulate and AM/PM are not defined. Dividing '
                'by it would report a huge modulation of nothing. Pick a '
                'harmonic the circuit actually produces.'
                % (carrier, abs(C), scale))
        upper = self.adjoint_sideband_row(pss, freq, output, carrier)[0]
        lower = self.adjoint_sideband_row(pss, -freq, output, carrier)[0]
        return self.am_pm_indices(upper / C, lower / C)

    def _deflated_solve(self, pss, alpha, b, transposed=False, tol=None):
        """`(I - alpha M) y = b` on an OSCILLATOR, with the pole taken out.

        ⚠ THE SINGULARITY IS THE ANSWER'S OWN POLE, NOT A NUMERICAL DEFECT,
        and that reframing is what makes the fix obvious.  At `alpha = 1`
        the operator is `I - M`, singular by the unit multiplier, and the
        solution really does diverge -- an oscillator's phase response to a
        perturbation goes as `1/df`.  What is wrong is COMPUTING a `1/eps`
        quantity through a system whose conditioning is also `1/eps`: the
        answer is genuinely large and the digits are genuinely gone.

        So the pole is carried ANALYTICALLY.  With `u`, `v` the right and
        left null vectors of `I - M` -- the orbit tangent and the PPV, both
        of which `ppv()` already returns -- border the system:

            [ I - alpha M   u ] [ w ]   [ b ]
            [     v^T       0 ] [ s ] = [ 0 ]

        Because `v^T (I - alpha M) = (1 - alpha) v^T` and `v^T w = 0`, the
        border variable comes out BOUNDED, `s = (v^T b)/(v^T u)`, with no
        `1/eps` in it.  And since `(I - alpha M) u = (1 - alpha) u`, the
        solution is recovered as

            y = w + s u / (1 - alpha)

        where `1 - alpha = 1 - exp(-j w T)` is the factor that vanishes at
        every harmonic, evaluated in closed form rather than inverted
        numerically.

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

        ⚠ IT STILL DIVERGES AT AN EXACT HARMONIC, and it should: `1/(1 -
        alpha)` is then a division by zero, and the physical response is
        unbounded.  What changes is that every offset NEAR a harmonic is
        now well conditioned, which is where phase noise is measured.

        `transposed` solves `(I - alpha M^T) x = b`, whose null space is
        spanned by `v` and whose left null space is spanned by `u`, so the
        borders swap.
        """
        import scipy.sparse.linalg as spla
        fp = pss.factored_period()
        n = fp.width
        _v, info = pss.ppv()
        v = np.asarray(_v, dtype=float)
        u = np.asarray(info['tangent_pair'], dtype=float)
        vu = float(v @ u)
        if abs(vu) < 1e-300:
            raise ValueError(
                'PAC: the PPV is orthogonal to the orbit tangent, so the '
                'bordering is singular and the pole cannot be removed.')
        col, row = (v, u) if transposed else (u, v)
        mv = (fp.matvec_transposed if transposed else fp.matvec)
        b = np.asarray(b, dtype=complex).ravel()

        def _mv(z):
            z = np.asarray(z)
            w_, s_ = z[:n], z[n]
            top = w_ - alpha * np.asarray(mv(w_)) + s_ * col
            return np.concatenate((top, [complex(row @ w_)]))

        A = spla.LinearOperator((n + 1, n + 1), matvec=_mv, dtype=complex)
        rhs = np.concatenate((b, [0.0 + 0.0j]))
        rt = max(self.KRYLOV_FACTOR * pss.par.reltol if tol is None else tol,
                 1e-14)
        z = self._gmres_checked(A, rhs, rt, 'the deflated solve')
        w, s = z[:n], z[n]
        denom = 1.0 - alpha
        if denom == 0:
            raise ValueError(
                'PAC: the deflated solve was asked for an EXACT harmonic, '
                'where 1/(1 - alpha) is a division by zero and the physical '
                'response is unbounded. The pole is removed from the '
                'CONDITIONING, not from the answer.')
        return w + s * col / denom

    def _op(self, fp, alpha):
        """`v -> (I - alpha M) v`, never forming `M`."""
        return lambda v: np.asarray(v) - alpha * fp.matvec(v)

    def _solve_each(self, fp, alphas, rhs, tol):
        """One GMRES per frequency -- the baseline the sweep is measured against."""
        import scipy.sparse.linalg as spla
        n = fp.width
        ys, count = [], [0]

        for alpha, b in zip(alphas, rhs):
            op = self._op(fp, alpha)

            def _mv(v, _op=op):
                count[0] += 1
                return _op(v)

            A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
            ys.append(self._gmres_checked(A, b, tol, 'the m x m solve'))
        return ys, count[0]

    def _solve_subspace(self, fp, alphas, rhs, tol):
        """ONE Krylov subspace for the whole sweep -- the recycling.

        ⚠ THE SUBSPACE IS FREQUENCY-INDEPENDENT AND THAT IS THE WHOLE POINT.
        `A(alpha) = I - alpha M`, so

            span{r, A r, A^2 r, ...} = span{r, M r, M^2 r, ...}

        for every `alpha` -- Telichevesky et al.'s Theorem 1.  A basis of
        `M`'s Krylov space therefore serves every frequency, and each one
        costs a small dense least-squares over it instead of its own run of
        full-period replays.

        ⚠ WHAT IS NOT FREE is the right-hand side: `w(f)` genuinely differs
        per frequency, so a basis grown from one frequency's residual is not
        the space GMRES would have chosen for another.  This does not
        guess -- it minimises the TRUE residual over the shared span, checks
        it, and extends the basis (one matvec, kept for every later
        frequency) until every frequency is inside tolerance.  So the answer
        is never worse than the per-frequency solve; only the matvec count
        varies.
        """
        n = fp.width
        V = np.zeros((n, 0), dtype=complex)
        MV = np.zeros((n, 0), dtype=complex)
        count = [0]

        def extend(seed):
            """One Arnoldi step of `M` from `seed`, orthogonal to `V`."""
            nonlocal V, MV
            v = np.asarray(seed, dtype=complex).ravel().copy()
            if V.shape[1]:
                v = v - V @ (V.conj().T @ v)
                v = v - V @ (V.conj().T @ v)   ## reorthogonalise once
            nv = np.linalg.norm(v)
            if nv < 1e-14:
                return False
            v = v / nv
            count[0] += 1
            Mv = np.asarray(fp.matvec(v), dtype=complex)
            V = np.hstack((V, v[:, None]))
            MV = np.hstack((MV, Mv[:, None]))
            return True

        extend(rhs[0])
        ys = [None] * len(alphas)
        pending = list(range(len(alphas)))
        for _round in range(min(n, 200)):
            still = []
            for i in pending:
                alpha, b = alphas[i], rhs[i]
                AV = V - alpha * MV
                y, *_ = np.linalg.lstsq(AV, b, rcond=None)
                x = V @ y
                r = b - (x - alpha * (MV @ y))
                ys[i] = x
                nb = np.linalg.norm(b)
                if np.linalg.norm(r) > tol * (nb if nb else 1.0):
                    still.append((i, r))
            if not still:
                return ys, count[0]
            if V.shape[1] >= n:
                break
            ## grow the shared basis on the worst residual -- one matvec,
            ## and every frequency gets to use it
            worst = max(still, key=lambda p: np.linalg.norm(p[1]))
            if not extend(worst[1]):
                break
            pending = [i for i, _r in still]
        return ys, count[0]
