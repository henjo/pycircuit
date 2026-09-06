from abc import ABC, abstractmethod

from pycircuit.circuit._lte_kernels import (bdf2_alphas, bdf2_companion,
                                            bdf2_derivative,
                                            euler_companion,
                                            trapezoidal_companion,
                                            second_divided_difference,
                                            third_divided_difference as _tdd,
                                            ## STAGE 12B -- d(iq)/dh, the
                                            ## integrator half of Fang's `p`.
                                            euler_companion_dh,
                                            trapezoidal_companion_dh,
                                            bdf2_companion_dh)
import math
import warnings

## STAGE 4e -- the zero-stability bound on the step-size ratio for variable-step
## BDF-2.  Grigorieff's result is that the homogeneous recursion's parasitic root
## stays inside the unit disc only while `h_n/h_{n-1} < 1 + sqrt(2)`:
##
##     ratio    parasitic root   growth over 20 steps
##     2.414214       1.000000                      1     <- the bound
##     2.5            1.041667                  2.262
##     3.0            1.285714                  152.4
##     10.0           4.761905               3.59e+13
##
## It bounds *growth* only.  Shrinking is unconditionally zero-stable, which is
## the whole content of 4e: the guard this constant now protects used to fire on
## the shrink and leave the growth unwatched.
ZERO_STABILITY_RATIO = 1.0 + math.sqrt(2.0)

## STAGE 9(a).  The definition moved to `_lte_kernels` so the JAX backend can
## reach the same one; re-exported here because this is where callers and tests
## have always imported it from, and moving a name is not this change's job.
third_divided_difference = _tdd



## WHY THERE IS NO `lte_formula` PARAMETER.  Removed in stage 9(f), 2026-07-31.
##
## It chose between the classic divided-difference estimates and the
## Yao-Wang-Roychowdhury Table I residuals.  Three changes removed its effect on
## this backend before it was removed as API:
##
##   4g(b)  the trapezoidal estimator stopped differencing `g` (the companion
##          current), which carries an undamped (-1)^n mode;
##   4i     both second-order estimators moved onto a shared third-derivative
##          estimate taken from a divided difference of the CHARGE, which reads
##          neither formula;
##   4d     the one-step fallback -- the last place either branch ran -- now takes
##          the divided-difference form unconditionally, because YWR's TRAP entry
##          is a uniform-grid formula and its GEAR2 residual is 3/4 of the true
##          truncation error.
##
## So `'ywr'` and `'classic'` produced bit-identical runs for every integrator, and
## the parameter was kept for a while, accepted and documented as inert.  What
## settled its removal was the OTHER backend: `jaxtransient.py` carried its own
## `lte_formula`, where `'classic'` selected a charge-domain estimator whose
## tolerance applied `lte_abs = 1e-6` -- a VOLTAGE floor -- to a CHARGE.  One
## microcoulomb, against node charges of pico- to femtocoulombs, so the normalized
## error could never reach 1 and no step was ever rejected: the controller ran
## open-loop.  One parameter name meant "selects nothing" here and "selects a
## broken estimator" there, which is worse than either alone.
##
## Both are gone.  Each backend now has exactly one estimator, and passing
## `lte_formula=` raises TypeError rather than being silently ignored -- a kwarg
## accepted and discarded is how the JAX defect stayed invisible.


class Integrator(ABC):
    """
    Abstract Base Class for Transient Numerical Integration Strategies.
    
    This interface isolates all mathematical discretization for time-stepping.
    By subclassing Integrator, a LINEAR-MULTISTEP method (Gear-3, say) can be
    added without modifying the core Transient solvers: it states its
    companion recursion and the loop drives it through `compute_derivatives`.
    ⚠ A STAGE method cannot -- TRBDF2 (2026-09-05) has two implicit stages,
    not a companion recursion, and needed a dedicated two-stage step in the
    Transient loop (the generic `_solve_timestep_rk`) and its own m x m shooting
    monodromy; its `compute_derivatives`/`companion_coefficients`/`companion_dT`
    raise.  The ABC still fits it structurally (order, history, order-drop),
    just not the single-companion time step.
    """

    ## --- CAPABILITY QUERIES (polymorphic dispatch) ---
    ## The shooting/transient stacks used to branch on `isinstance(...)` and
    ## method-name strings at ~35 sites; these let a caller ask the METHOD
    ## instead, so a new integrator arrives with the right answers rather than
    ## needing an edit at each site.  See doc/integrator_architecture_260906.md.

    def is_stage_method(self) -> bool:
        """True for a Runge-Kutta / stage method, False for a linear-multistep
        companion method.  The single predicate the shooting stack branches on
        to pick the stage machinery over the LMM machinery."""
        return False

    def companion_reach(self) -> int:
        """How many charges back this method's companion reads (Euler/trap 1,
        Gear-2 2) -- decides whether the shooting solve needs the entering
        history as an unknown.  Asked of the method, not inferred from a name."""
        alphas, _b = self.companion_coefficients(1.0, 1.0)
        return len(alphas) - 1

    def carries_own_monodromy(self) -> bool:
        """Whether this method's own period map is already second-order (or
        higher) on a limit cycle, so `monodromy_twin` need not borrow one:
        Gear-2 (reach 2) and every stage method (no opener seam).  The one-step
        LMMs (trap/euler) return False -- their native monodromy is first
        order and they take a twin."""
        return self.is_stage_method() or self.companion_reach() >= 2

    def needs_x0_unknown(self) -> bool:
        """Whether the shooting solve must carry `x_0` as the unknown -- the
        self-starting stage methods, which have no manufactured opening step to
        differentiate `x_0` back through."""
        return self.is_stage_method()

    @abstractmethod
    def get_required_history(self) -> int:
        """
        Returns the number of past time steps required by the method.
        For example: Euler=1, Gear-2=2.
        """
        pass
        
    @abstractmethod
    def check_order_drop(self, h_curr: float, h_last: float, is_first_step: bool) -> 'Integrator':
        """
        Evaluates sudden step variations. High-order methods can return 
        an instance of a lower-order method (like EulerIntegrator) to 
        safely cross discontinuities if the step size shrinks aggressively.
        
        Returns:
            The Integrator instance to use for the current step.
        """
        pass
        
    def companion_coefficients(self, h_curr, h_last):
        """``(alphas, b)`` of this method's companion recursion.

        Every method here writes the companion current as

            iq_n = sum_k alphas[k] q_{n-k}  +  b iq_{n-1}

        and those numbers are all a caller needs to differentiate ONE STEP
        with respect to the state it came from -- which is what a shooting
        method's monodromy is made of.  Returned rather than transcribed
        because the alternative is a third copy of the coefficients living
        in `shooting.py`, and this file already records what transcribing an
        integration constant costs (the 3/4 optimism, found three times).

        `alphas[0]` multiplies the CURRENT charge, so `alphas[0] * C` is the
        `geq` :meth:`compute_derivatives` returns; the rest are the past.
        """
        raise NotImplementedError(
            '%s does not state its companion coefficients, so it cannot be '
            'used where a per-step sensitivity is needed (shooting).'
            % type(self).__name__)

    def companion_dT(self, q_curr, q_last, h_curr, h_last):
        """``d(iq)/dT`` when EVERY step scales together -- the shooting half.

        ⚠ NOT `companion_dh`, AND THE DIFFERENCE IS A FACTOR OF 3/2 FOR
        GEAR-2.  `companion_dh` is a PARTIAL: `d/dh_n` with `h_{n-1}` held
        fixed, which is right for the coupled time-stepping method it was
        written for -- there `h_{n-1}` really is a past step that is not
        moving.  A shooting analysis solving for the PERIOD rebuilds its
        grid at the current `T`, so every step is `c_k T` and `h_{n-1}`
        moves too; the total derivative needs that route as well.

        Euler and trapezoidal do not notice, because their coefficients
        depend on `h_n` alone and the partial IS the total.  Gear-2's depend
        on both, and on a uniform grid the missing route is exactly half the
        one that is there -- measured against finite differences on the
        autonomous system, the ratio of the code's column to the true one
        was 1.4859 / 1.4939 / 1.4972 at 100 / 200 / 400 points per period,
        converging on 3/2.

        ⚠ ONE IMPLEMENTATION FOR EVERY METHOD, and it does not differentiate
        anything.  The `alphas` are homogeneous of degree -1 in the step
        sizes -- they must be, since `iq` approximates `dq/dt` -- so Euler's
        theorem gives

            sum_j h_j d(alpha_k)/d(h_j)  =  -alpha_k

        and therefore `T d(iq)/dT = -sum_k alpha_k q_{n-k}` exactly, with no
        per-method partial to write down or get wrong.  Verified numerically
        for the variable-step BDF-2 coefficients to 6.4e-08, which is the
        finite-difference floor.  Writing three more derivative routines
        instead is what this file already warns against; the coefficients
        are taken from `companion_coefficients` rather than transcribed.

        Returns `T d(iq)/dT`, i.e. WITHOUT the `1/T`, because the integrator
        does not know the period -- the caller divides.
        """
        alphas, _b = self.companion_coefficients(h_curr, h_last)
        acc = alphas[0] * q_curr
        for k in range(1, len(alphas)):
            acc = acc + alphas[k] * q_last[k - 1]
        return -acc

    def companion_dh(self, q_curr, q_last, h_curr, h_last):
        """``d(iq)/dh`` at fixed solution -- the integrator half of Fang's ``p``.

        STAGE 12B.  Not abstract: an integrator that does not implement it simply
        cannot be used with the coupled time-stepping method, and raising here
        says so at the point of use rather than silently contributing zero to
        ``p``, which would look like a converged solve of the wrong problem.
        """
        raise NotImplementedError(
            '%s does not provide d(iq)/dh, so it cannot be used with the '
            'coupled (Fang) time-stepping method' % type(self).__name__)

    @abstractmethod
    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        """
        Calculates the numerical derivative of the charge (iq) and the 
        equivalent conductance (geq).
        
        Returns:
            Tuple[iq, geq]
        """
        pass
        
    @abstractmethod
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit,
                    h_last2=None) -> tuple:
        """
        Computes the Local Truncation Error vector for the current step,
        along with the order 'p' of the LTE formula (used for step size prediction).

        ``h_last2`` is the step BEFORE ``h_last``, added by stage 4g(b) for the
        trapezoidal estimator.  It is ``None`` exactly when ``q_last[2]`` is not yet
        a real past point -- the two become available on the same step, so one
        optional argument carries both facts and an estimator that needs three past
        charges can test ``h_last2 is None`` and fall back.  Estimators that need
        only two past points ignore it.

        UNITS -- read this before comparing the return value to any tolerance.
        The vector is named for the local truncation error and written ``Eg``, but it
        is a **current**, not a charge and not a voltage.  Every second-order estimator
        here returns a multiple of ``h^2 q'''``, and with ``q'''`` in C/s^3 that is C/s.
        Backward Euler's ``(q_n - q_{n-1})/h - iq_{n-1}`` is a difference of companion
        currents, likewise amperes.

        The controller consumes it by mapping it through ``J^-1`` into the solution
        domain and comparing *that* against ``reltol``/``vabstol``/``iabstol``, which is
        dimensionally sound.  Comparing the raw return value against a charge tolerance
        is not: decision 0.3d's option (D) was designed around exactly that and was
        refuted on it, having reproduced the units defect gate 0.2b recorded on the JAX
        backend.  If a future estimator needs a charge, it is ``h`` times this.

        Returns:
            Tuple[lte_vector, p] -- lte_vector in AMPERES, p the order
        """
        pass

class EulerIntegrator(Integrator):
    """Backward Euler (1st order) Integration Method"""

    ## STAGE 12B.  The order of the method, needed by `SolutionLTEController`:
    ## Fang's eq (6) compares the solution against a polynomial extrapolation,
    ## and that difference is only a truncation-error estimate when the
    ## extrapolation degree EQUALS the method order (the Milne device).  Get it
    ## wrong and the two errors no longer cancel in a controlled way -- measured
    ## with a degree-2 predictor against this order-1 corrector as an error that
    ## moved 600x for a 2x step change, where the model says 4x.
    ORDER = 1

    def __init__(self):
        ## No `lte_formula`: see the module note above.  For
        ## Backward Euler the two formulas always coincided, so this class never
        ## had a choice to preserve across an order drop in the first place.
        pass

    def get_required_history(self) -> int:
        return 1

    def check_order_drop(self, h_curr: float, h_last: float, is_first_step: bool) -> Integrator:
        # Euler is 1st order, no lower order to drop to.
        return self

    def companion_coefficients(self, h_curr, _h_last):
        ## iq_n = (q_n - q_{n-1})/h, and no iq history.  `_h_last` is
        ## underscored because a one-step method genuinely does not read it
        ## -- the uniform signature exists so the caller can ask any
        ## integrator without knowing which, and the dead-knob scan is right
        ## to want the difference visible.
        return (1.0 / h_curr, -1.0 / h_curr), 0.0
        
    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        return euler_companion(q_curr, C_curr, q_last[0], h_curr)

    def companion_dh(self, q_curr, q_last, h_curr, h_last):
        return euler_companion_dh(q_curr, q_last[0], h_curr)
        
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit,
                    h_last2=None):
        ## `h_last2` is unused: Euler differences one past point.  Accepted so every
        ## implementation shares one signature -- see the ABC.
        if is_first_step:
            return toolkit.zeros(len(q_curr)), 1.0

        gn = (q_curr - q_last[0]) / h_curr
        gn_1 = iq_last[0]

        ## STAGE 4c -- THE VARIABLE-STEP CORRECTION.
        ##
        ## Backward Euler's companion current is `(q_n - q_{n-1})/h`, which is a
        ## centred approximation of `q'` at the MIDPOINT of the step, not at the
        ## node.  So `g_n - g_{n-1}` differences two midpoint derivatives separated
        ## by `(h_curr + h_last)/2`, not by `h_curr`, and it therefore estimates
        ## `((h1+h2)/2) q''` where the truncation error is `(h1/2) q''`.
        ##
        ## On a uniform grid the two coincide and the estimator is exact, which is
        ## why this went unnoticed.  Off it, the estimate is wrong by
        ## `(h1+h2)/(2 h1)` -- measured est/true, before this correction:
        ##
        ##     ratio  0.25    0.5     1.0     2.0     4.0
        ##     est/true  2.5246  1.5089  1.0040  0.7522  0.6265
        ##
        ## which is a 4.03x spread across the sweep, and it is the wrong direction
        ## twice over: on a shrinking step the error is OVERstated (so the
        ## controller shrinks further than it needs to) and on a growing step it is
        ## UNDERstated (so the controller grows past what the tolerance allows).
        ##
        ## Rescaling by `2 h1 / (h1 + h2)` converts the midpoint spacing back to the
        ## step.  After: 1.0098 / 1.0059 / 1.0040 / 1.0029 / 1.0024 -- flat to
        ## within 1% across the same sweep.
        scale = 2.0 * h_curr / (h_curr + h_last) if h_last else 1.0
        lte = -0.5 * scale * (gn - gn_1)
        return lte, 2.0  # p=2.0 for Euler LTE h-dependence in step controller

class TrapezoidalIntegrator(Integrator):
    """Trapezoidal (2nd order) Integration Method"""

    ORDER = 2

    def __init__(self):
        ## No `lte_formula`: removed in 9(f) -- see the module note above.
        pass

    def get_required_history(self) -> int:
        ## THREE, not one, since stage 4g(b).  The *method* still looks back one
        ## step -- `compute_derivatives` uses q_last[0] and iq_last[0] only -- but
        ## the ESTIMATOR differences `d_k = (q_k - q_{k-1})/h_k` at three points,
        ## and d_{n-2} needs q_{n-3}.  The two requirements are different and this
        ## returns the larger, because it is what sizes the ring buffer.
        return 3

    def check_order_drop(self, h_curr: float, h_last: float, is_first_step: bool) -> Integrator:
        # Trapezoidal rule only looks back 1 step, so its polynomial isn't
        # distorted by past step sizes. No drop needed.
        if is_first_step:
            return EulerIntegrator()
        return self
        
    def companion_coefficients(self, h_curr, _h_last):
        ## iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}.  The `-1` is the term that
        ## makes the period map carry `iq`, and the undamped (-1)^n mode the
        ## LTE estimator has to avoid differencing.  One step back, so
        ## `_h_last` is unread here -- see EulerIntegrator.
        return (2.0 / h_curr, -2.0 / h_curr), -1.0

    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        ## `2*C/h`, which is what `C/h/0.5` computed -- both are exact scalings by
        ## two in binary floating point, so this is bit-identical, not merely equal.
        return trapezoidal_companion(q_curr, C_curr, q_last[0], iq_last[0], h_curr)

    def companion_dh(self, q_curr, q_last, h_curr, h_last):
        return trapezoidal_companion_dh(q_curr, q_last[0], h_curr)
        
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit,
                    h_last2=None):
        if is_first_step:
            return toolkit.zeros(len(q_curr)), 1.0

        ## STAGE 4g(b) -- DIFFERENCE A MODE-FREE QUANTITY.
        ##
        ## What the estimator must produce, from YWR eq (22) with p=1, k=2 and the
        ## trapezoidal coefficients (alpha = [1/h, -1/h], beta = [1/2, 1/2]), is
        ##
        ##     Eg = -(h^2/6) q'''
        ##
        ## once the controller's own J^-1 has absorbed the (q_x + 0.5 h f_x)^-1
        ## factor.  Table I approximates q''' by a second difference of the
        ## companion current g; eq (22) does not require that, and the paper's own
        ## wording is "(22) AND FINITE DIFFERENCE APPROXIMATION" -- the choice of
        ## difference is free, and for TRAP the g-based choice is what goes wrong.
        ##
        ## g carries an undamped parasitic mode.  The trapezoidal companion
        ## `iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}` has homogeneous solution
        ## `iq_n = -iq_{n-1}`, i.e. `(-1)^n`, and nothing damps it.  Differencing g
        ## therefore differences that mode, and the estimate depends on the step
        ## history that preceded it: measured est/true_local at h=1e-9, ratio 1, on
        ## two different prefixes of the SAME problem, 1.3176 and 0.6780 -- a 1.9x
        ## swing from history alone.
        ##
        ## `d` has no such component.  The trapezoidal relation gives
        ## `(iq_n + iq_{n-1})/2 = (q_n - q_{n-1})/h = d_n`, and the mode flips sign
        ## every step, so it CANCELS EXACTLY in that sum.  Expanding about the
        ## interval midpoint `m_n = (t_n + t_{n-1})/2`,
        ##
        ##     d_n = q'(m_n) + (h_n^2/24) q'''(m_n) + O(h^4)
        ##
        ## so d samples q' at midpoints and a second divided difference of d OVER
        ## THE MIDPOINTS estimates q'''/2.  Measured against the local truncation
        ## error, est/true at ratio 1: 0.9273 / 0.9933 / 0.9993 as h falls
        ## 1e-8 -> 1e-10, against the g-based form's 0.8067 / 0.6780 / 0.6678 --
        ## asymptotically exact where the old one holds a 33% underestimate.
        ##
        ## The midpoint spacings are the part that must not be got wrong; using
        ## h_curr/h_last (the NODE spacings) is the same class of error as the
        ## backward-Euler defect stage 4c fixed.
        if h_last2 is not None:
            ## STAGE 4i.  Eg = -(h^2/6) q''' and the shared helper returns q'''/6,
            ## so the whole formula is -h^2 times it.
            ##
            ## 4g(b) differenced `d_k = (q_k - q_{k-1})/h_k` over the interval
            ## midpoints instead, which removed the parasitic mode and was
            ## asymptotically exact but kept a +-12% bias at the extremes of the
            ## reachable step-ratio range: `d_k` carries `(h_k^2/24) q'''`, and that
            ## term cancels only on a uniform grid.  The charge carries no method
            ## error at all, so differencing it directly removes the residual --
            ## measured spread over ratio 0.008..2.414 falls from 1.26x to 1.008x.
            return -(h_curr**2) * third_divided_difference(
                q_curr, q_last, h_curr, h_last, h_last2), 3.0

        ## `h_last2 is None` means q_last[2] is not yet a real past point, which is
        ## true for exactly ONE step of a run -- the second, where the ring buffer
        ## still holds the initial charge twice.  Falling back to the g-based form
        ## for that single step is better than returning zeros, which would make it
        ## an unchecked step of the kind stage 3 exists to remove.
        ## THE DIVIDED-DIFFERENCE FORM, AND `lte_formula` DOES NOT SELECT HERE.
        ## YWR's Table I TRAP entry -- `Eg = -(1/6)(g_n - 2 g_{n-1} + g_{n-2})` --
        ## carries a single `h` and an UNWEIGHTED second difference, i.e. it is a
        ## uniform-grid formula, where the same table's GEAR2 entry carries h1 and
        ## h2 explicitly.  Off a uniform grid it is wrong by O(1/h), and the grid is
        ## not uniform here: stage 3's opening ramp is still growing the step at the
        ## one point this fallback runs.  Taking the divided-difference form
        ## unconditionally is stage 4d's stated fix.
        gn = 2 * (q_curr - q_last[0]) / h_curr - iq_last[0]
        gn_1 = iq_last[0]
        gn_2 = iq_last[1] if len(iq_last) > 1 else iq_last[0]
        dd1 = (gn - gn_1) / h_curr
        dd2 = (gn_1 - gn_2) / h_last
        lte = -(1.0/3.0) * h_curr**2 * (dd1 - dd2) / (h_curr + h_last)
        return lte, 3.0  # p=3.0 for Trapezoidal


class Gear2Integrator(Integrator):
    """Gear-2 / BDF-2 (2nd order) Variable Step Size Integration Method"""

    ORDER = 2

    def __init__(self):
        ## No `lte_formula`: removed in 9(f) -- see the module note above.
        ## The history is kept because both entries are recorded results, and both
        ## are about a choice that no longer exists rather than about this class.
        ##
        ## THE DEFAULT USED TO BE 'ywr', chosen belt-and-braces when 'classic' was
        ## repaired, on the grounds that it had the longer track record.  The price
        ## was that the YWR GEAR2 residual estimates (1/4) h^2 q''' against a true
        ## (1/3) h^2 q''', so it reported 3/4 of the truncation error at every step
        ## where a corrected 'classic' is asymptotically exact.  Stage 4i moved both
        ## variants onto the divided-difference form, which is the 'classic' one, so
        ## that optimism is gone rather than merely defaulted around.
        ##
        ## THE "5/6" AN EARLIER COMMENT CLAIMED WAS AN ARTEFACT, worth keeping
        ## because the number is so clean.  5/6 = 0.8333 is what the trapezoidal
        ## estimator reads when it is handed EXACT derivatives as its `g` history
        ## instead of the companion currents a real run produces.  Measured against
        ## the local truncation error with the real history it reads 1.09 / 1.31 /
        ## 1.33 as h falls 1e-8 -> 1e-10 -- it does not converge at all, let alone
        ## to 5/6.  Decision 0.3b called the claim measurably wrong; this is the
        ## measurement, and the mechanism.
        pass

    def get_required_history(self) -> int:
        ## THREE since stage 4i, for the same reason trapezoidal needs three: the
        ## METHOD looks back two steps -- `compute_derivatives` uses q_last[0] and
        ## q_last[1] -- but the ESTIMATOR takes a third divided difference of the
        ## charge and so needs q_{n-3}.  Until 4g(b) built the `h_last2` plumbing
        ## this was not available, and the comment in `compute_lte` below recorded
        ## it as the reason the g-based form had to be used.
        return 3

    def check_order_drop(self, h_curr: float, h_last: float, is_first_step: bool) -> Integrator:
        if is_first_step:
            return EulerIntegrator()

        ## STAGE 4e -- THE GUARD USED TO WATCH THE WRONG DIRECTION.
        ##
        ## The only test here was `if h_curr / h_last < 0.1`, and it was labelled
        ## as protecting the validity of the high-order polynomial.  It does not:
        ## variable-step BDF-2's parasitic root leaves the unit disc only on
        ## *growth*, past `ZERO_STABILITY_RATIO`, and any ratio below 1 is
        ## unconditionally zero-stable.  So the one ratio that can actually
        ## destabilise the recursion was unwatched -- which is how the 10x growth
        ## on `transient.py`'s force-accept path (4b) survived: nothing downstream
        ## would have caught it.  (The shrink test itself is kept, for a different
        ## and measured reason; see below.)  Measured on the stiff RLC, this new
        ## branch fires 0 times, because the controller's own clamp
        ## (`MAX_GROWTH_RATIO` = 2.0) keeps every normal step inside the bound.
        ##
        ## **That makes this a backstop, and a backstop that never fires in a
        ## healthy run is the point of it** -- it is what turns "no accepted step
        ## ratio exceeds the bound" from an accident of two clamps agreeing into
        ## something the integrator enforces for itself.  Dropping to Euler is the
        ## right response rather than refusing the step: order 1 has no parasitic
        ## root to amplify, so the ratio becomes harmless instead of forbidden.
        if h_curr / h_last > ZERO_STABILITY_RATIO:
            return EulerIntegrator()

        ## THE SHRINK BRANCH IS KEPT, AND RE-LABELLED.  The plan said replace; the
        ## measurement said add, so it is added and the reason is written down.
        ##
        ## Removing it outright took `Gear2('ywr')` and `Gear2('classic')` from 0
        ## force-accepts to 1 each on the stiff RLC at reltol 1e-5, because it is
        ## not idle: it fires 3-6 times a run there.  What it is doing is nothing
        ## to do with zero-stability -- it is a STALLED-ESTIMATE heuristic.  A step
        ## only shrinks 10x below the last accepted one after several consecutive
        ## rejections, and what rejects repeatedly is a 2nd-order estimate built on
        ## a third difference of a solution that is not three times differentiable
        ## -- i.e. a discontinuity.  Dropping to order 1 there is what every
        ## simulator does across a corner, and it is the same medicine 4b now
        ## administers at the rejection cap, one retry later and after having
        ## accepted an over-tolerance step to get there.  Deleting it would have
        ## traded a controlled Euler step for a force-accepted 2nd-order one.
        ##
        ## So: the guard above is the stability bound, this one is economics, and
        ## the defect 4e names was never that this branch existed -- it was that
        ## this branch was ALL there was, and it was labelled as protecting a
        ## stability property it has nothing to do with.  **Reconsider if** the
        ## rejection cap ever becomes rejection-count-aware: `h_curr/h_last < 0.1`
        ## is a proxy for "we have rejected three times at this time point", and
        ## the thing it is a proxy for is known exactly one level up in
        ## `transient.py`, where it would not need a threshold at all.
        if h_curr / h_last < 0.1:
            return EulerIntegrator()

        return self
        
    def companion_coefficients(self, h_curr, h_last):
        ## The variable-step BDF-2 coefficients themselves -- iq_n is a
        ## combination of THREE charges and carries no iq history, so `b` is
        ## zero and the past reaches two steps back.  Taken from
        ## `bdf2_alphas` rather than restated: fixed-step formulas are wrong
        ## the moment the step changes, which is why they are computed at
        ## all.
        a0, a1, a2 = bdf2_alphas(h_curr, h_last)
        return (a0, a1, a2), 0.0

    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        # --- VARIABLE STEP-SIZE BDF-2 (GEAR-2) DERIVATION ---
        # Traditional SPICE2 uses fixed-step BDF formulas which fail when dt changes.
        # Here we calculate the true Variable Step-Size (VSS) coefficients for a 2nd-order 
        # polynomial fit through the current point (n) and two previous points (n-1, n-2).
        # These coefficients mathematically convert the continuous time derivative dq/dt 
        # into a discrete algebraic equivalent.
        ## STAGE 9(a) -- one definition, shared with jaxtransient.  These three
        ## coefficients had three copies in source and two in tests.
        return bdf2_companion(q_curr, C_curr, q_last[0], q_last[1],
                              h_curr, h_last)

    def companion_dh(self, q_curr, q_last, h_curr, h_last):
        return bdf2_companion_dh(q_curr, q_last[0], q_last[1], h_curr, h_last)
        
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit,
                    h_last2=None):
        if is_first_step:
            return toolkit.zeros(len(q_curr)), 1.0

        ## STAGE 4i -- THE ESTIMATOR USED TO DIFFERENCE THE METHOD'S OWN ERROR.
        ##
        ## Gear-2's local truncation error is `-(1/6) h1 (h1+h2) q'''`, so every
        ## companion current in the history carries an error of exactly that shape.
        ## Both branches below take a second divided difference of `g` at the
        ## nodes, which differences those errors along with the signal.  The
        ## damage was computed by hand before it was measured, and the two agree to
        ## 0.3% at every step ratio:
        ##
        ##     h1/h2      0.008    0.05     0.1    0.25       1       2       4
        ##     predicted  83.34  13.365   6.727   2.800  1.0000   0.778   0.700
        ##     measured   83.06   13.32    6.71    2.79   0.998   0.775   0.695
        ##
        ## It vanishes exactly at h1 = h2 = h3, which is why the estimator measured
        ## asymptotically exact (1.000282 against 2/9) on a uniform grid: that
        ## measurement was taken at the one ratio where the defect is zero.  A step
        ## ratio of 0.008 is reached after three consecutive rejections, so the
        ## worst case is not hypothetical -- it is the step where the controller
        ## has just collapsed the step size and is told the error is 83x worse than
        ## it is.
        ##
        ## The fix is to estimate q''' from the CHARGES, which carry no method
        ## error.  The obstacle used to be real and was recorded here: Gear-2 kept
        ## only two past charges, so a third divided difference was unavailable.
        ## 4g(b) lifted it.
        if h_last2 is not None:
            ## Eg = -(1/6) h1 (h1+h2) q''', and the helper returns q'''/6.
            return -h_curr * (h_curr + h_last) * third_divided_difference(
                q_curr, q_last, h_curr, h_last, h_last2), 3.0

        ## `h_last2 is None` for exactly one step of a run -- the second, before the
        ## ring buffer holds four real charges.  The g-based form below serves that
        ## step; returning zeros would make it unchecked, which is the defect stage
        ## 3 removed from the first step.
        ##
        ## IT IS THE DIVIDED-DIFFERENCE FORM, NOT YWR's, AND `lte_formula` DOES NOT
        ## SELECT HERE.  YWR's Table I GEAR2 residual estimates `(1/4) h^2 q'''`
        ## against a true `(1/3)`, so it reports 3/4 of the truncation error --
        ## measured on this exact fallback as -2.827659e+01 where the correct value
        ## is -3.770212e+01.  After 4i this was the ONLY step of a run where the
        ## choice still had any effect on the CPU path, so taking the accurate one
        ## unconditionally is what finishes 4d: "delete the branch and keep 'ywr' as
        ## an alias".
        ## (The YWR Table I GEAR2 residual that used to be selectable here,
        ##  `Eg = -(1/8)((h1+h2)/(h1 h2))(h2 g_n - (h1+h2) g_{n-1} + h1 g_{n-2})`,
        ##  is derived and compared against this one in doc/src/circuit/lte_dae.rst;
        ##  it is not kept as dead code.)

        # --- CLASSIC GEAR-2 LOCAL TRUNCATION ERROR ---
        # Taylor-expanding the VSS companion current above about t_n (the alpha
        # coefficients kill the q' and q'' terms by construction) leaves
        #
        #     iq - q'(t_n) = -(1/6) h1 (h1 + h2) q'''(t_n) + O(h^3)
        #
        # -- equal steps: -(1/3) h^2 q''', the textbook BDF-2 result.  So what
        # has to be estimated here is the THIRD derivative of the charge, scaled
        # by h^2.
        #
        # A second divided difference of q yields only q'', and Gear-2 keeps just
        # two past charges (get_required_history() == 2), so a third divided
        # difference of q is not available at all.  The third derivative is
        # therefore taken as the second divided difference of g = dq/dt, read off
        # the companion-current history -- the same information the YWR branch
        # above uses.  Estimating q'' here and multiplying by h^3 (as this branch
        # did until 2026-07) is dimensionally not a current: it undershoots the
        # truncation error by a factor of order h*omega, which on a 1 MHz signal
        # at nanosecond steps is ~1e-15.  The controller then never rejects a
        # step, saturates the growth limiter every step, pins h at max_step and
        # stops responding to reltol/abstol altogether.
        h1, h2 = h_curr, h_last
        g_n = bdf2_derivative(q_curr, q_last[0], q_last[1], h1, h2)
        g_nm1 = iq_last[0]
        g_nm2 = iq_last[1] if len(iq_last) > 1 else iq_last[0]

        # Second divided difference of g at t_n, t_{n-1}, t_{n-2}, which is
        # q'''/2, so the -(1/6) above becomes -(1/3) here.
        dd2_g = ((g_n - g_nm1) / h1 - (g_nm1 - g_nm2) / h2) / (h1 + h2)

        lte = -(1.0 / 3.0) * h1 * (h1 + h2) * dd2_g

        return lte, 3.0  # p=3.0 is order+1 (the estimate itself scales with h^2)


class RungeKuttaIntegrator(Integrator):
    """Base for one-step Runge-Kutta (stage) methods, carrying the Butcher
    tableau as the single source of truth.

    A concrete method supplies three class attributes -- ``A`` (the ``s x s``
    Butcher matrix), ``B`` (the ``s`` step weights) and ``C`` (the ``s`` node
    abscissae) -- and OPTIONALLY the embedded-estimator data.  Everything the
    transient loop and the shooting stack need for a stage method (the coupled
    stage step, the monodromy, the forward/adjoint source folds, the Lyapunov
    injection, the cost transform, the embedded error estimate) is a function of
    that tableau, so it is written ONCE against this base rather than re-derived
    per method.  See ``doc/integrator_architecture_260906.md``.

    The one-step, self-starting contract is answered here for every subclass:
    ``get_required_history() == 1``, ``companion_reach() == 1``,
    ``is_stage_method() == True``, and the linear-multistep companion methods
    (``companion_coefficients``/``companion_dT``/``compute_lte``/
    ``compute_derivatives``) inherit the base's refusal -- a stage method has no
    companion recursion.
    """

    #: structure tags returned by :meth:`stage_structure`
    FULL = 'full'
    DIRK = 'dirk'
    SDIRK = 'sdirk'
    ESDIRK = 'esdirk'

    def is_stage_method(self) -> bool:
        return True

    def get_required_history(self) -> int:
        ## self-starting: the only past state is x_n, already carried by the loop
        return 1

    def companion_reach(self) -> int:
        return 1

    def check_order_drop(self, h_curr, h_last, is_first_step):
        ## a one-step method carries no zero-stability step-ratio limit
        return self

    def butcher(self):
        """``(A, B, C)`` as float ndarrays, cached."""
        cache = getattr(self, '_butcher_cache', None)
        if cache is None:
            import numpy as np
            cache = (np.array(self.A, dtype=float),
                     np.array(self.B, dtype=float),
                     np.array(self.C, dtype=float))
            self._butcher_cache = cache
        return cache

    @property
    def stages(self) -> int:
        import numpy as np
        return int(np.array(self.C).shape[0])

    def is_stiffly_accurate(self) -> bool:
        """``b == last row of A`` and ``c[-1] == 1`` -- the step lands ON the
        constraint manifold, so ``x_{n+1}`` is the last stage."""
        import numpy as np
        A, B, C = self.butcher()
        return bool(np.allclose(B, A[-1]) and abs(C[-1] - 1.0) < 1e-14)

    def is_explicit_first_stage(self) -> bool:
        import numpy as np
        A, _B, _C = self.butcher()
        return bool(np.allclose(A[0], 0.0))

    def stage_structure(self):
        """Classify the tableau so the solver picks the cheapest correct path:
        ``ESDIRK`` (lower-triangular, explicit first stage, equal remaining
        diagonals -> one shared factorisation), ``SDIRK`` (lower-triangular,
        all diagonals equal), ``DIRK`` (lower-triangular), or ``FULL`` (fully
        implicit -> coupled solve or the eig(A^-1) cost transform)."""
        import numpy as np
        A, _B, _C = self.butcher()
        s = A.shape[0]
        lower = np.allclose(np.triu(A, 1), 0.0)
        if not lower:
            return self.FULL
        diag = np.diag(A)
        expl0 = abs(diag[0]) < 1e-14
        impl = diag[1:] if expl0 else diag
        equal = impl.size > 0 and np.allclose(impl, impl[0]) and impl[0] != 0.0
        if expl0 and equal:
            return self.ESDIRK
        if equal:
            return self.SDIRK
        return self.DIRK

    def is_fully_implicit(self) -> bool:
        return self.stage_structure() == self.FULL

    ## The two remaining abstract LMM methods: a stage method has no
    ## single-companion time step and no divided-difference LTE, so it refuses
    ## both here ONCE for every subclass -- a new RK method needs only its
    ## tableau, not a repeat of these.
    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last,
                            h_last, is_first_step, toolkit):
        raise NotImplementedError(
            '%s is a Runge-Kutta stage method; the Transient loop runs it via '
            'the generic _solve_timestep_rk (one coupled or sequential stage '
            'solve), not through the single-companion compute_derivatives.'
            % type(self).__name__)

    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last,
                    is_first_step, toolkit, h_last2=None):
        raise NotImplementedError(
            '%s states no linear-multistep companion, so the LMM '
            'divided-difference compute_lte does not apply; its embedded '
            'estimate is computed in the generic RK step and consumed by '
            '_run_rk_adaptive.' % type(self).__name__)


class TRBDF2Integrator(RungeKuttaIntegrator):
    """TR-BDF2: a trapezoid stage over ``gamma*h`` then a BDF2 stage over ``h``.

    A one-step, self-starting, L-stable, order-2 DIRK.  It is here for three
    failure modes it does NOT have rather than for accuracy (it is order 2,
    like trapezoidal and Gear-2): no manufactured opening step (so a shooting
    monodromy is not seeded at first order), no ``(-1)^n`` companion mode, and
    no zero-stability step-ratio limit (a frozen non-uniform grid carries no
    ``ZERO_STABILITY_RATIO`` penalty).

    ⚠ IT IS NOT A LINEAR-MULTISTEP METHOD, so it does not fit the
    companion-recursion shape the rest of this module and the ABC are built
    around: there is no ``iq_n = sum_k alphas[k] q_{n-k}`` -- there are STAGES.
    The Transient loop runs it through :meth:`Transient._solve_timestep_rk`
    (two Newton solves sharing ONE factorisation, because ``a22 == a33``), not
    through :meth:`compute_derivatives`, and shooting builds its ``m x m``
    monodromy directly from the two stage linearisations rather than from
    :meth:`companion_coefficients`.  Those three methods therefore raise here.

    ⚠ ADAPTIVE, but NOT through this class's LMM interface.  The embedded
    2(3) estimate (Hosea & Shampine 1996) is computed inside
    :meth:`Transient._solve_timestep_rk` from the three stage
    derivatives and filtered once through the stage matrix, and a dedicated
    driver :meth:`Transient._run_rk_adaptive` runs step control on it.
    :meth:`compute_lte` STILL raises: it is the LMM controller's
    divided-difference interface, which a two-stage DIRK does not fit -- the
    estimate is a stage combination, not a companion difference.  The
    estimator coefficients were DERIVED (Taylor matching) and verified
    against the analytic local error, not quoted.  :meth:`check_order_drop`
    returns ``self`` (a one-step method has no ratio limit to enforce).

    THE TABLEAU, derived not quoted (the corpus has nothing on TR-BDF2).
    ``gamma`` is fixed by the ONE-LU condition: the trapezoid stage diagonal
    ``gamma/2`` equals the BDF2 stage diagonal ``(1-gamma)/(2-gamma)`` iff
    ``gamma^2 - 4 gamma + 2 = 0``, i.e. ``gamma = 2 - sqrt(2)``.  Then the two
    stage Jacobians ``C + (gamma*h/2) G`` and ``C + a33*h G`` are the same
    matrix (``a22 == a33``): one factorisation per step, two solves.  Stiffly
    accurate by construction (``b`` is the last row of ``A``), so ``R(inf)=0``.

    ⚠ ONE QUADRATIC, REACHED FROM FOUR DIFFERENT REQUIREMENTS -- so
    ``gamma`` is not merely a cost constant.  Writing
    ``Q(gamma) = gamma^2 - 4 gamma + 2`` (whose root in ``(0,1)`` is
    ``2 - sqrt(2)``), each of the following is ``Q`` TIMES A FACTOR THAT
    DOES NOT VANISH ON ``(0,1)`` -- which is the load-bearing statement,
    because it makes ``Q = 0`` the ONLY solution of each, not merely one
    (verified symbolically here; the constants below are the factors'
    values at the root):

      * the ONE-LU condition ``a33 - gamma/2 = Q / (2(2-gamma))`` -- zero iff
        ``Q = 0`` (factor 0.354);
      * BANK et al. 1985 (eq. 38) give the principal truncation coefficient
        ``C(gamma) = (-3 gamma^2 + 4 gamma - 2)/(12(2-gamma))``, whose
        derivative is ``Q / (4 (2-gamma)^2)`` -- extremal iff ``Q = 0``
        (factor 0.125) -- and ``C(2-sqrt2) = 2/3 - sqrt(2)/2 ~ -0.0404``
        equals the embedded estimator's leading coefficient
        ``(4 - 3 sqrt2)/6`` (this file's own derivation, reached
        independently);
      * ROSENBROCK 1963 (Comput. J. 5(4), eq. 29) fixes the equal stage
        diagonal ``d = gamma/2`` by ``d^2 - 2d + 1/2 = 0``, which is
        ``Q/4 = 0``;
      * the RADIUS OF ABSOLUTE MONOTONICITY ``R(A,b) = 2(2-gamma) /
        (1 + (1-gamma)^2)`` (Bonaventura & Della Rocca 2015) has derivative
        ``2 Q / (1 + (1-gamma)^2)^2`` -- extremal iff ``Q = 0`` (factor
        1.457) -- and is MAXIMISED there at ``1 + sqrt(2) ~ 2.414``, against
        trapezoid's ``2`` and backward Euler's ``inf``.  Verified here from
        the Kraaijevanger 1991 definition directly (validated on CN = 2 and
        BE = inf first); it means TR-BDF2 tolerates a ~21% larger step than
        trapezoidal before monotonicity/positivity/TVD can fail.

    Four structurally unrelated requirements -- one LU, minimal truncation,
    L-stability of a linearly implicit process, maximal monotonicity radius
    -- all extremal at the one quadratic is a strong reason to trust the
    constant, stronger than any one alone.  (TR-BDF2 is L-stable for EVERY
    gamma -- ``R(inf) = 0`` identically -- so L-stability alone selects
    nothing; it is the equal-diagonal condition within Rosenbrock's family
    that gives ``Q``.)

    ⚠ CONVENTION: this is BANK / HOSEA & SHAMPINE's ``gamma`` (the stage
    ABSCISSA ``c2 = 2 - sqrt(2) ~ 0.586``), NOT Kennedy & Carpenter's or
    Rosenbrock's, whose ``gamma`` is this one halved (``= a22 = 0.293``, our
    ``STAGE_DIAG``).  Same method, symbol differing by exactly 2; do not
    "correct" the value against the other paper.
    """

    ORDER = 2
    #: order of the embedded estimate, so the adaptive controller uses the
    #: ``1/(EMBEDDED_ORDER+1)`` step exponent (2(3) pair -> 1/3).
    EMBEDDED_ORDER = 2
    #: embedded 2(3) weights on the three stage derivatives (Hosea & Shampine
    #: 1996): ``est_raw = h (c0 K0 + c1 K1 + c2 K2)`` is the leading LTE of the
    #: order-2 solution.  Read by the generic DIRK step's error estimate.
    EMBEDDED_DK = ((1.0 - math.sqrt(2.0)) / 3.0,
                   1.0 / 3.0,
                   -(2.0 - math.sqrt(2.0)) / 3.0)

    #: ``gamma = 2 - sqrt(2)``, from the one-LU condition ``g^2 - 4g + 2 = 0``.
    GAMMA = 2.0 - math.sqrt(2.0)
    #: The two stage diagonals, equal by the gamma condition: ``a22 == a33``.
    STAGE_DIAG = 1.0 - math.sqrt(2.0) / 2.0
    #: BDF2 stage weights on ``q(Y1)`` and ``q(xn)``; ``A1 + A0 == 1``.
    A1 = 0.5 + math.sqrt(2.0) / 2.0
    A0 = 0.5 - math.sqrt(2.0) / 2.0

    ## THE BUTCHER TABLEAU (3-stage ESDIRK form): an explicit first stage
    ## (``Y1 = x_n``, ``c1 = 0``) then the two implicit stages with equal
    ## diagonal ``d = STAGE_DIAG`` -- the trapezoid stage at ``c2 = 2d = gamma``
    ## and the stiffly-accurate BDF2 stage at ``c3 = 1``.  ``w = sqrt(2)/4`` are
    ## the outer weights (``2w + d == 1``).  Equal implicit diagonals ->
    ## ``stage_structure() == 'esdirk'`` -> one shared factorisation, which is
    ## exactly the one-LU property stated above.  This is the same method as the
    ## ``(GAMMA, A1, A0, STAGE_DIAG)`` companion form; both are kept until the
    ## generic RK step (which reads the tableau) replaces the bespoke step.
    _W = math.sqrt(2.0) / 4.0
    A = ((0.0, 0.0, 0.0),
         (STAGE_DIAG, STAGE_DIAG, 0.0),
         (_W, _W, STAGE_DIAG))
    B = (_W, _W, STAGE_DIAG)
    C = (0.0, GAMMA, 1.0)

    def __init__(self):
        pass

    def get_required_history(self) -> int:
        ## Self-starting and one-step: the only past state is ``x_n`` itself,
        ## which the loop already carries as the previous solution.
        return 1

    def check_order_drop(self, h_curr, h_last, is_first_step):
        ## No zero-stability ratio limit -- a one-step method cannot lose it
        ## across a step change -- so nothing to drop to.
        return self

    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last,
                            h_last, is_first_step, toolkit):
        raise NotImplementedError(
            'TR-BDF2 is a two-stage method; the Transient loop runs it via '
            'the generic RK stage step _solve_timestep_rk, not through the '
            'single-companion compute_derivatives.')

    def companion_coefficients(self, h_curr, h_last):
        raise NotImplementedError(
            'TR-BDF2 has stages, not a linear-multistep companion recursion, '
            'so it states no (alphas, b); shooting builds its m x m monodromy '
            'from the two stage linearisations directly.')

    def companion_dT(self, q_curr, q_last, h_curr, h_last):
        raise NotImplementedError(
            'TR-BDF2 states no companion coefficients, so the Euler-theorem '
            'd(iq)/dT shared by the LMMs does not apply; the autonomous '
            'shooting dT for a stage method is not yet built.')

    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last,
                    is_first_step, toolkit, h_last2=None):
        raise NotImplementedError(
            'TR-BDF2 states no linear-multistep companion, so the LMM '
            'controller\'s divided-difference compute_lte does not apply. Its '
            'embedded 2(3) estimate is computed in the generic RK step and '
            'consumed by Transient._run_rk_adaptive, the adaptive path for '
            'every stage method.')


class RadauIIA3Integrator(RungeKuttaIntegrator):
    """Radau IIA, 3 stages: the order-5, L-stable, stiffly-accurate collocation
    method on the two Radau points and the endpoint.

    Where TR-BDF2 earns its place by the failure modes it lacks, Radau IIA(3)
    earns its place by *coverage and accuracy*:

    * **Order 5, stage order 3.**  Five times the classical order of the LMMs
      and of TR-BDF2, and -- crucially for a DAE -- stage order 3, so the
      algebraic (index-1) components do not suffer the order reduction that
      collapses a DIRK to its stage order on the constraint manifold.
    * **The only candidate covered on a DAE.**  ``det A = 1/60 != 0`` makes the
      stage system invertible, which is exactly the hypothesis of Hairer &
      Wanner VI.2 Thm 2.3 (convergence on index-1 DAEs).  A method with a
      singular ``A`` (any explicit-first-stage DIRK, TR-BDF2 included) is not
      covered by that theorem on the algebraic block.
    * **L-stable, stiffly accurate.**  ``R(inf) = 1 - b^T A^{-1} 1 = 0`` and the
      last stage IS the step (``b == A[-1, :]``, ``c[-1] == 1``), so the
      numerical solution lands ON the constraint manifold each step.

    The price is that it is FULLY implicit: the three stages are coupled into
    one ``3n`` system, with no explicit first stage to unlock and no
    per-stage one-LU shortcut.  The Transient loop therefore runs it through a
    dedicated coupled solve (``_solve_timestep_radau``), not through the
    single-companion ``compute_derivatives`` path.

    COST TRANSFORM (documented; the coupled solve is the default).  The ``3n``
    Newton system ``(I3 (x) C/h + A (x) G) dY = -F`` block-diagonalises through
    the eigendecomposition of ``A^{-1}``:  its spectrum is one real eigenvalue
    ``GAMMA_REAL`` and one complex pair ``ALPHA +- i BETA`` (stored below).  In
    the eigenbasis the coupled solve becomes one real factorisation of
    ``(GAMMA_REAL/h) C + G`` plus one complex factorisation of
    ``((ALPHA + i BETA)/h) C + G`` -- 1 real + 1 complex LU instead of a dense
    ``3n`` solve.  Realising that win on the sparse backend needs the complex
    ``klu_z_*`` binding; until then the coupled real solve is correct (and, for
    the small circuits here, cheap), and the transform is an efficiency
    follow-up rather than a correctness requirement.

    ⚠ CONVENTION: ``c = [2/5 - sqrt(6)/10, 2/5 + sqrt(6)/10, 1]`` (the two
    interior Radau points and the endpoint) with ``b`` equal to the LAST row of
    ``A`` (stiff accuracy).  This is the IIA family (right Radau, endpoint
    included), NOT IA (left Radau, ``c1 = 0``); do not swap the abscissae.
    """

    #: Classical order (B(5) holds; C(3) gives stage order 3).
    ORDER = 5
    #: order of the embedded 5(3) estimate -> ``1/4`` adaptive step exponent.
    EMBEDDED_ORDER = 3
    #: Number of coupled stages.
    STAGES = 3

    _S6 = math.sqrt(6.0)
    #: Abscissae: two interior Radau points and the stiffly-accurate endpoint.
    C = (2.0 / 5.0 - _S6 / 10.0, 2.0 / 5.0 + _S6 / 10.0, 1.0)
    #: The 3x3 collocation matrix ``A`` (rows sum to ``C``; ``det A = 1/60``).
    A = (
        (11.0 / 45.0 - 7.0 * _S6 / 360.0,
         37.0 / 225.0 - 169.0 * _S6 / 1800.0,
         -2.0 / 225.0 + _S6 / 75.0),
        (37.0 / 225.0 + 169.0 * _S6 / 1800.0,
         11.0 / 45.0 + 7.0 * _S6 / 360.0,
         -2.0 / 225.0 - _S6 / 75.0),
        (4.0 / 9.0 - _S6 / 36.0,
         4.0 / 9.0 + _S6 / 36.0,
         1.0 / 9.0),
    )
    #: Step weights = last row of ``A`` (stiff accuracy).
    B = A[2]
    #: Eigenvalues of ``A^{-1}`` for the cost transform: one real, one pair.
    GAMMA_REAL = 3.6378342527444960
    ALPHA = 2.6810828736277523
    BETA = 3.0504301992474105

    def __init__(self):
        pass

    def get_required_history(self) -> int:
        ## Self-starting collocation: the only past state is ``x_n`` itself,
        ## already carried by the loop as the previous solution.
        return 1

    def check_order_drop(self, h_curr, h_last, is_first_step):
        ## A one-step collocation method carries no zero-stability step-ratio
        ## limit across a step change -- nothing to drop to.
        return self

    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last,
                            h_last, is_first_step, toolkit):
        raise NotImplementedError(
            'Radau IIA(3) is a three-stage fully-implicit method; the Transient '
            'loop runs it via _solve_timestep_radau (one coupled 3n Newton '
            'solve), not through the single-companion compute_derivatives.')

    def companion_coefficients(self, h_curr, h_last):
        raise NotImplementedError(
            'Radau IIA(3) states no linear-multistep companion; its stages are '
            'coupled through the collocation matrix A, not a scalar companion.')

    def companion_dT(self, q_curr, q_last, h_curr, h_last):
        raise NotImplementedError(
            'Radau IIA(3) states no companion coefficients, so the Euler-theorem '
            'd(iq)/dT shared by the LMMs does not apply.')

    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last,
                    is_first_step, toolkit, h_last2=None):
        raise NotImplementedError(
            'Radau IIA(3) states no linear-multistep companion, so the LMM '
            'divided-difference compute_lte does not apply. Its embedded 5(3) '
            'estimate is computed in Transient._solve_timestep_radau and '
            'consumed by _run_radau_adaptive.')


class ESDIRK43Integrator(RungeKuttaIntegrator):
    """ESDIRK4(3)6L[2]SA -- Kennedy & Carpenter's 6-stage, order-4, L-stable,
    stiffly-accurate ESDIRK (the "KenCarp4" scheme), with its embedded order-3
    estimate.

    ⚠ THIS IS A TEST VEHICLE, and its whole point is that adding it took ONLY
    the Butcher tableau below -- no new transient step, no new shooting
    monodromy, no new noise or forced-fold code.  It exercises the generic
    Runge-Kutta machinery at ``s = 6`` stages (against TR-BDF2's 3 and Radau's
    3), which is the check that the s-stage generalisation of the DIRK-sequential
    family is real and not tuned to a stage count.

    ⚠ LOWER-TRIANGULAR with an EXPLICIT first stage and a CONSTANT implicit
    diagonal ``gamma = 1/4`` (``stage_structure() == 'esdirk'``): it takes the
    sequential DIRK path, never the coupled/transform one -- ``A`` is singular
    (the explicit first stage), so it lacks the Butcher-Bickart transform and
    the H&W index-1 DAE theorem that Radau's ``det A != 0`` buys.  ESDIRK buys
    its cheapness from the triangular constant-diagonal structure instead (one
    shared factorisation per step), which is a different bargain.
    """

    ORDER = 4
    EMBEDDED_ORDER = 3
    STAGES = 6

    _g = 1.0 / 4.0
    #: Butcher matrix (6x6, lower-triangular, ESDIRK: explicit first stage,
    #: a_ii = 1/4 for i >= 1).  Kennedy & Carpenter 2003, ARK4(3)6L[2]SA.
    A = (
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.0 / 4.0, 1.0 / 4.0, 0.0, 0.0, 0.0, 0.0),
        (8611.0 / 62500.0, -1743.0 / 31250.0, 1.0 / 4.0, 0.0, 0.0, 0.0),
        (5012029.0 / 34652500.0, -654441.0 / 2922500.0, 174375.0 / 388108.0,
         1.0 / 4.0, 0.0, 0.0),
        (15267082809.0 / 155376265600.0, -71443401.0 / 120774400.0,
         730878875.0 / 902184768.0, 2285395.0 / 8070912.0, 1.0 / 4.0, 0.0),
        (82889.0 / 524892.0, 0.0, 15625.0 / 83664.0, 69875.0 / 102672.0,
         -2260.0 / 8211.0, 1.0 / 4.0),
    )
    #: step weights = last row of A (stiffly accurate).
    B = A[5]
    #: node abscissae (rows of A sum to these).
    C = (0.0, 1.0 / 2.0, 83.0 / 250.0, 31.0 / 50.0, 17.0 / 20.0, 1.0)
    #: order-3 embedded weights ``b_hat``.
    B_HAT = (4586570599.0 / 29645900160.0, 0.0, 178811875.0 / 945068544.0,
             814220225.0 / 1159782912.0, -3700637.0 / 11593932.0,
             61727.0 / 225920.0)
    #: embedded estimate weights on the stage derivatives: ``b_hat - b``, so
    #: ``yhat - y = h sum_i (b_hat_i - b_i) K_i`` (read by the generic DIRK
    #: step's error estimate).
    EMBEDDED_DK = tuple(bh - b for bh, b in zip(B_HAT, B))

    def __init__(self):
        pass
