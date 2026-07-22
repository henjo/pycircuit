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
    
    def get_required_history(self) -> int:
        return 1
        
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
        
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        if is_first_step:
            return toolkit.zeros(len(q_curr)), 1.0
            
        gn = 2 * (q_curr - q_last[0]) / h_curr - iq_last[0]
        gn_1 = iq_last[0]
        gn_2 = iq_last[1] if len(iq_last) > 1 else iq_last[0]
        
        # Trapezoidal LTE
        dd1 = (gn - gn_1) / h_curr
        dd2 = (gn_1 - gn_2) / h_last
        lte = -(1.0/3.0) * h_curr**2 * (dd1 - dd2) / (h_curr + h_last)
        return lte, 3.0  # p=3.0 for Trapezoidal


class Gear2Integrator(Integrator):
    """Gear-2 / BDF-2 (2nd order) Variable Step Size Integration Method"""
    
    def get_required_history(self) -> int:
        return 2
        
    def check_order_drop(self, h_curr: float, h_last: float, is_first_step: bool) -> Integrator:
        if is_first_step:
            return EulerIntegrator()
            
        # Bizzarri & Brambilla Order Drop Protection
        # If the step shrinks by more than 10x, high-order polynomials become invalid.
        if h_curr / h_last < 0.1:
            return EulerIntegrator()
            
        return self
        
    def compute_derivatives(self, q_curr, C_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        # Dynamically compute VSS BDF-2 coefficients
        alpha0 = (2 * h_curr + h_last) / (h_curr * (h_curr + h_last))
        alpha1 = -(h_curr + h_last) / (h_curr * h_last)
        alpha2 = h_curr / (h_last * (h_curr + h_last))
        
        geq = C_curr * alpha0
        iq = alpha0 * q_curr + alpha1 * q_last[0] + alpha2 * q_last[1]
        return iq, geq
        
    def compute_lte(self, q_curr, h_curr, q_last, iq_last, h_last, is_first_step, toolkit):
        if is_first_step:
            return toolkit.zeros(len(q_curr)), 1.0
            
        # Dynamically compute VSS Gear2 LTE
        # Note: the LTE formula currently used in PyCircuit computes difference of derivatives
        dd1_n = (q_curr - q_last[0]) / h_curr
        dd1_nm1 = (q_last[0] - q_last[1]) / h_last
        dd2_n = (dd1_n - dd1_nm1) / (h_curr + h_last)
        lte = (h_curr**2) * (h_curr + h_last) / 3.0 * dd2_n
        
        return lte, 2.0  # p=2.0 for gear2
