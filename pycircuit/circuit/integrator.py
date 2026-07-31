from abc import ABC, abstractmethod
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

def third_divided_difference(q_curr, q_last, h_curr, h_last, h_last2):
    """Estimate `q'''/6` from four charges, exactly and without a ratio-dependent bias.

    STAGE 4i.  Both second-order methods here need the same quantity: their local
    truncation errors are `-(1/6) h1 (h1+h2) q3` for Gear-2 and `-(1/6) h1^2 q3` for
    trapezoidal, where `q3` is the third derivative of the charge.  Only the
    constant differs, so only the constant belongs in the integrators.

    THE POINT IS WHAT IS DIFFERENCED.  By the mean-value form of the polynomial
    interpolation error,

        q[t_n, t_{n-1}, t_{n-2}, t_{n-3}] = q3(zeta) / 6

    for some `zeta` in the span, **with no coefficient that depends on the step
    ratios**.  Every previous formulation in this file lacked that property and was
    wrong off a uniform grid in consequence:

      * differencing `g` (the companion currents) differences the method's OWN
        truncation error along with the signal, because each `g_k` carries
        `-(1/6) h_k (h_k + h_{k-1}) q3`.  Predicted and measured at 83x for Gear-2
        at a step ratio of 0.008, which three consecutive rejections reach.
      * differencing `d_k = (q_k - q_{k-1})/h_k` -- what trapezoidal used between
        4g(b) and 4i -- is far better, because `d` carries no parasitic mode, but
        `d_k = q'(m_k) + (h_k^2/24) q3(m_k)` still has a term that depends on
        `h_k` and so does not cancel off a uniform grid.  Worth +-12% at the
        extremes.

    The charge itself carries no method error at all, which is why this form
    measures flat to 1.004x across the whole reachable ratio range.
    """
    d1 = (q_curr - q_last[0]) / h_curr
    d2 = (q_last[0] - q_last[1]) / h_last
    d3 = (q_last[1] - q_last[2]) / h_last2
    dd_a = (d1 - d2) / (h_curr + h_last)
    dd_b = (d2 - d3) / (h_last + h_last2)
    return (dd_a - dd_b) / (h_curr + h_last + h_last2)


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
    By subclassing Integrator, new methods (like Gear-3 or TRBDF2) can be added 
    without modifying the core Transient solvers.
    """
    
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
        
    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        resultEuler = (q_curr - q_last[0]) / h_curr
        iq = resultEuler
        geq = C_curr / h_curr
        return iq, geq
        
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
        
    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        resultTrap = 2 * (q_curr - q_last[0]) / h_curr - iq_last[0]
        iq = resultTrap
        geq = C_curr / h_curr / 0.5
        return iq, geq
        
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
        
    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        # --- VARIABLE STEP-SIZE BDF-2 (GEAR-2) DERIVATION ---
        # Traditional SPICE2 uses fixed-step BDF formulas which fail when dt changes.
        # Here we calculate the true Variable Step-Size (VSS) coefficients for a 2nd-order 
        # polynomial fit through the current point (n) and two previous points (n-1, n-2).
        # These coefficients mathematically convert the continuous time derivative dq/dt 
        # into a discrete algebraic equivalent.
        alpha0 = (2 * h_curr + h_last) / (h_curr * (h_curr + h_last))
        alpha1 = -(h_curr + h_last) / (h_curr * h_last)
        alpha2 = h_curr / (h_last * (h_curr + h_last))
        
        # Equivalent conductance represents the 'algorithmic resistance' introduced by the timestep
        geq = C_curr * alpha0
        # The numerical derivative of charge
        iq = alpha0 * q_curr + alpha1 * q_last[0] + alpha2 * q_last[1]
        return iq, geq
        
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
        alpha0 = (2 * h1 + h2) / (h1 * (h1 + h2))
        alpha1 = -(h1 + h2) / (h1 * h2)
        alpha2 = h1 / (h2 * (h1 + h2))
        g_n = alpha0 * q_curr + alpha1 * q_last[0] + alpha2 * q_last[1]
        g_nm1 = iq_last[0]
        g_nm2 = iq_last[1] if len(iq_last) > 1 else iq_last[0]

        # Second divided difference of g at t_n, t_{n-1}, t_{n-2}, which is
        # q'''/2, so the -(1/6) above becomes -(1/3) here.
        dd2_g = ((g_n - g_nm1) / h1 - (g_nm1 - g_nm2) / h2) / (h1 + h2)

        lte = -(1.0 / 3.0) * h1 * (h1 + h2) * dd2_g

        return lte, 3.0  # p=3.0 is order+1 (the estimate itself scales with h^2)
