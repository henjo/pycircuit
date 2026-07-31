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
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit) -> tuple:
        """
        Computes the Local Truncation Error vector for the current step, 
        along with the order 'p' of the LTE formula (used for step size prediction).
        
        Returns:
            Tuple[lte_vector, p]
        """
        pass

class EulerIntegrator(Integrator):
    """Backward Euler (1st order) Integration Method"""

    def __init__(self, lte_formula='classic'):
        ## 'classic' or 'ywr' (Yao-Wang-Roychowdhury DAE LTE, ICECS 2014).  For
        ## Backward Euler the two coincide, so this is only carried through so an
        ## order-drop from Gear2/Trap preserves the chosen formula.
        self.lte_formula = lte_formula

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
        
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
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

    def __init__(self, lte_formula='classic'):
        self.lte_formula = lte_formula

    def get_required_history(self) -> int:
        return 1

    def check_order_drop(self, h_curr: float, h_last: float, is_first_step: bool) -> Integrator:
        # Trapezoidal rule only looks back 1 step, so its polynomial isn't
        # distorted by past step sizes. No drop needed.
        if is_first_step:
            return EulerIntegrator(self.lte_formula)
        return self
        
    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        resultTrap = 2 * (q_curr - q_last[0]) / h_curr - iq_last[0]
        iq = resultTrap
        geq = C_curr / h_curr / 0.5
        return iq, geq
        
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        if is_first_step:
            return toolkit.zeros(len(q_curr)), 1.0
            
        # g = dq/dt companion current at the last three points (g_n reconstructed
        # from the trapezoidal companion formula, g_{n-1}, g_{n-2} from history).
        gn = 2 * (q_curr - q_last[0]) / h_curr - iq_last[0]
        gn_1 = iq_last[0]
        gn_2 = iq_last[1] if len(iq_last) > 1 else iq_last[0]

        if self.lte_formula == 'ywr':
            # Yao-Wang-Roychowdhury DAE LTE (ICECS 2014, Table I, TRAP):
            #   eps = -(1/12)(q_x + 0.5 h f_x)^-1 (g_n - 2 g_{n-1} + g_{n-2}) h
            # With (q_x + 0.5 h f_x) = 0.5 h (G + 2C/h) = 0.5 h J, the h cancels
            # and, since the controller applies J^-1, Eg = -(1/6)(2nd diff of g).
            lte = -(1.0/6.0) * (gn - 2 * gn_1 + gn_2)
        else:
            # Classic divided-difference form.
            dd1 = (gn - gn_1) / h_curr
            dd2 = (gn_1 - gn_2) / h_last
            lte = -(1.0/3.0) * h_curr**2 * (dd1 - dd2) / (h_curr + h_last)
        return lte, 3.0  # p=3.0 for Trapezoidal


class Gear2Integrator(Integrator):
    """Gear-2 / BDF-2 (2nd order) Variable Step Size Integration Method"""

    def __init__(self, lte_formula='ywr'):
        ## Defaults to 'ywr' (Yao-Wang-Roychowdhury, ICECS 2014) rather than
        ## 'classic', belt and braces: 'classic' is now mathematically correct
        ## (see compute_lte below) but 'ywr' has the longer track record here.
        ## The price is that the YWR GEAR2 residual estimates (1/4) h^2 q'''
        ## against a true (1/3) h^2 q''', so the default reports 3/4 of the
        ## truncation error at every step where a corrected 'classic' is
        ## asymptotically exact -- mild optimism about the solver's own error,
        ## and TRTOL = 7.0 already absorbs more than that factor.
        ## Euler and Trapezoidal keep 'classic': for Euler the two formulas are
        ## identical, and for Trapezoidal they agree to the same 5/6 of the
        ## one-step LTE.
        self.lte_formula = lte_formula

    def get_required_history(self) -> int:
        return 2

    def check_order_drop(self, h_curr: float, h_last: float, is_first_step: bool) -> Integrator:
        if is_first_step:
            return EulerIntegrator(self.lte_formula)

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
            return EulerIntegrator(self.lte_formula)

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
            return EulerIntegrator(self.lte_formula)

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
        
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        if is_first_step:
            return toolkit.zeros(len(q_curr)), 1.0
            
        if self.lte_formula == 'ywr':
            # Yao-Wang-Roychowdhury DAE LTE (ICECS 2014, Table I, GEAR2).  Uses the
            # 2nd difference of g = dq/dt (not a q'' divided difference).  h1=h_curr,
            # h2=h_last; g_n is the Gear2 companion current reconstructed from the
            # charge history, g_{n-1}, g_{n-2} come from the iq history.  The paper's
            # (q_x + [h1(h1+h2)/(2h1+h2)] f_x)^-1 factor equals alpha0*J^-1, and since
            # the controller applies J^-1, the residual reduces to:
            #   Eg = -(1/8) ((h1+h2)/(h1 h2)) (h2 g_n - (h1+h2) g_{n-1} + h1 g_{n-2})
            h1, h2 = h_curr, h_last
            alpha0 = (2 * h1 + h2) / (h1 * (h1 + h2))
            alpha1 = -(h1 + h2) / (h1 * h2)
            alpha2 = h1 / (h2 * (h1 + h2))
            g_n = alpha0 * q_curr + alpha1 * q_last[0] + alpha2 * q_last[1]
            g_nm1 = iq_last[0]
            g_nm2 = iq_last[1] if len(iq_last) > 1 else iq_last[0]
            lte = -(1.0/8.0) * ((h1 + h2) / (h1 * h2)) * \
                (h2 * g_n - (h1 + h2) * g_nm1 + h1 * g_nm2)
            return lte, 3.0

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
