from abc import ABC, abstractmethod
import warnings

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
        lte = -0.5 * (gn - gn_1)
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

    def __init__(self, lte_formula='classic'):
        self.lte_formula = lte_formula

    def get_required_history(self) -> int:
        return 2

    def check_order_drop(self, h_curr: float, h_last: float, is_first_step: bool) -> Integrator:
        if is_first_step:
            return EulerIntegrator(self.lte_formula)

        # Bizzarri & Brambilla Order Drop Protection
        # If the step shrinks by more than 10x, high-order polynomials become invalid.
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
