from abc import ABC, abstractmethod
import numpy as np

class StepController(ABC):
    """
    Abstract Strategy Interface for deciding and predicting time steps.
    """
    
    @abstractmethod
    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, is_first_step, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0):
        """
        Evaluates the Local Truncation Error (LTE) for the current step.
        Returns:
            accept_step: bool
            h_next: float
        """
        pass

class IntegralController(StepController):
    """
    Standard Integral Step Controller based on Yao et al. ICECS 2014.
    Rejects steps with LTE > 1.0, and predicts the next step size.
    """
    
    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, is_first_step, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0):
        if is_first_step:
            err = 0.5
            return True, h_curr
        
        # --- LOCAL TRUNCATION ERROR (LTE) CALCULATION ---
        # 1. Ask the active integrator (e.g. Gear2, Trapezoidal) to calculate the 
        #    raw unscaled truncation error vector (Eg) based on charge curvature.
        Eg, p = active_integrator.compute_lte(
            q_curr=q_curr,
            h_curr=h_curr,
            q_last=q_last_hist,
            iq_last=iq_last_hist,
            h_last=h_last,
            is_first_step=is_first_step,
            toolkit=toolkit
        )
        
        # 2. Convert charge error into voltage error by multiplying by the 
        #    inverse of the Jacobian matrix: lte = J^-1 * Eg
        from pycircuit.circuit.analysis import remove_row_col
        J_reduced, Eg_reduced = remove_row_col((J, Eg), irefnode, toolkit)
        
        try:
            lte_reduced = toolkit.linearsolver(J_reduced, Eg_reduced)
        except Exception:
            lte_reduced = Eg_reduced
            
        lte = toolkit.concatenate((lte_reduced[:irefnode], toolkit.array([0.0]), lte_reduced[irefnode:]))
        
        # 3. Calculate dynamic tolerance for each node based on its current voltage level
        etol = reltol * toolkit.maximum(abs(x_curr), abs(x_last)) + abstol
        
        # 4. Normalize the error
        err_array = abs(lte) / etol
        err = float(np.max(err_array))
        
        # --- STEP REJECTION / ACCEPTANCE ---
        if err > 1.0:
            # Step rejected: The curvature was too sharp for the current timestep.
            # We shrink h_next and return False to force the solver to recalculate.
            h_next = h_curr * max(0.2, (1.0 / err)**0.5)
            return False, h_next

            
        # Step accepted: predict next size
        h_next = h_curr * min(2.0, (TRTOL / max(err, 1e-12))**0.5)
        h_next = min(h_next, max_step)
        
        return True, h_next

class PIController(StepController):
    """
    Proportional-Integral Step Controller.
    Uses history of truncation error to provide smoother step size changes.
    """
    def __init__(self, k_i=0.3, k_p=0.4):
        self.k_i = k_i
        self.k_p = k_p
        self.last_err = None
        
    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, is_first_step, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0):
        if is_first_step:
            err = 0.5
            self.last_err = err
            return True, h_curr
        
        Eg, p = active_integrator.compute_lte(
            q_curr=q_curr,
            h_curr=h_curr,
            q_last=q_last_hist,
            iq_last=iq_last_hist,
            h_last=h_last,
            is_first_step=is_first_step,
            toolkit=toolkit
        )
        
        from pycircuit.circuit.analysis import remove_row_col
        J_reduced, Eg_reduced = remove_row_col((J, Eg), irefnode, toolkit)
        
        try:
            lte_reduced = toolkit.linearsolver(J_reduced, Eg_reduced)
        except Exception:
            lte_reduced = Eg_reduced
            
        lte = toolkit.concatenate((lte_reduced[:irefnode], toolkit.array([0.0]), lte_reduced[irefnode:]))
        
        etol = reltol * toolkit.maximum(abs(x_curr), abs(x_last)) + abstol
        err_array = abs(lte) / etol
        
        err = float(np.max(err_array))
        
        if err > 1.0:
            # Step rejected
            # When rejected, we use a standard backoff
            h_next = h_curr * max(0.2, (1.0 / err)**0.5)
            return False, h_next
            
        # Step accepted: use PI update
        if self.last_err is None:
            self.last_err = err
            
        # Standard PI formula for step size control
        err_norm = max(err, 1e-12) / TRTOL
        err_last_norm = max(self.last_err, 1e-12) / TRTOL
        
        factor = (err_norm ** (-self.k_i)) * ((err_last_norm / err_norm) ** self.k_p)
        
        # Limit the step size change
        factor = min(2.0, max(0.2, factor))
        
        h_next = h_curr * factor
        h_next = min(h_next, max_step)
        
        self.last_err = err
        
        return True, h_next
