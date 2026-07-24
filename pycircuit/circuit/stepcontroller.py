from abc import ABC, abstractmethod
import numpy as np

class StepController(ABC):
    """
    Abstract Strategy Interface for deciding and predicting time steps.
    """
    
    @abstractmethod
    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, is_first_step, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0):
        """Evaluate the Local Truncation Error (LTE) for the current step.

        Returns:
            tuple: ``(accept_step, h_next)`` -- whether the step is accepted and
            the predicted next step size.
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

        # 3. Dynamic per-node tolerance, relaxed by the transient error factor TRTOL.
        #    TRTOL is the SPICE "transient tolerance" (Spectre calls the equivalent
        #    `lteratio`): the LTE estimate is deliberately conservative, so the
        #    allowed truncation error is TRTOL times the Newton-solve tolerance.
        #    Folding TRTOL into etol makes the accept threshold (err<=1) and the
        #    step-size prediction aim at the same target instead of oscillating.
        etol = TRTOL * (reltol * toolkit.maximum(abs(x_curr), abs(x_last)) + abstol)

        # 4. Normalize the error
        err_array = abs(lte) / etol
        err = float(np.max(err_array))

        # Step-size prediction exponent is 1/(order+1); compute_lte returns that
        # (order+1) as ``p`` (2 for Euler, 3 for the 2nd-order methods).
        exponent = 1.0 / p
        safety = 0.9

        # --- STEP REJECTION / ACCEPTANCE ---
        if err > 1.0:
            # Step rejected: shrink and recalculate.
            h_next = h_curr * max(0.2, safety * (1.0 / err)**exponent)
            return False, h_next

        # Step accepted: predict next size with the same controller law.
        h_next = h_curr * min(2.0, safety * (1.0 / max(err, 1e-12))**exponent)
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

        # Relax the LTE tolerance by TRTOL (see IntegralController) so the accept
        # threshold matches the target the PI update drives toward.
        etol = TRTOL * (reltol * toolkit.maximum(abs(x_curr), abs(x_last)) + abstol)
        err_array = abs(lte) / etol

        err = float(np.max(err_array))
        exponent = 1.0 / p

        if err > 1.0:
            # Step rejected: standard backoff using the method order.
            h_next = h_curr * max(0.2, 0.9 * (1.0 / err)**exponent)
            return False, h_next
            
        # Step accepted: use PI update
        if self.last_err is None:
            self.last_err = err
            
        # Standard PI formula for step size control.  err is already normalized so
        # that err==1 is the target (TRTOL is folded into etol above).
        err_norm = max(err, 1e-12)
        err_last_norm = max(self.last_err, 1e-12)

        factor = (err_norm ** (-self.k_i)) * ((err_last_norm / err_norm) ** self.k_p)
        
        # Limit the step size change
        factor = min(2.0, max(0.2, factor))
        
        h_next = h_curr * factor
        h_next = min(h_next, max_step)
        
        self.last_err = err
        
        return True, h_next
