from abc import ABC, abstractmethod
import numpy as np

from pycircuit.circuit.analysis import NoConvergenceError
class NonLinearSolver(ABC):
    """
    Abstract Base Class for Non-Linear Algebraic Solvers.
    
    This interface isolates the Newton-Raphson iteration math and continuation 
    methods (like Gmin-stepping). It is entirely agnostic to whether it is 
    solving a DC operating point or a Transient time-step.
    """
    
    @abstractmethod
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None):
        pass


class StandardNewton(NonLinearSolver):
    """
    Standard Full-Matrix Newton-Raphson Solver.
    Used natively by DCAnalysis and standard adaptive Transient simulations.
    """
    
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None):
        x = x0
        for i in range(maxiter):
            F, J = eval_FJ(x)
            
            try:
                xdiff = toolkit.linearsolver(J, -F)
            except Exception as e:
                raise NoConvergenceError(f"Singular Jacobian: {str(e)}")
            
            x_next = x + xdiff
            if limiter is not None:
                x_next = limiter(x_next, x)
                # Recompute xdiff based on limited step
                xdiff = x_next - x
                
            # KCL Scale: Upper bound of absolute branch currents/voltages
            I_scale = toolkit.dot(abs(J), abs(x_next)) + abs(F)
            
            conv_x = toolkit.alltrue(abs(xdiff) < reltol * toolkit.maximum(abs(x_next), abs(x)) + xtol)
            conv_f = toolkit.alltrue(abs(F) < reltol * I_scale + abstol)
            
            if conv_x and conv_f:
                return x_next, i + 1
                
            x = x_next
            
        raise NoConvergenceError(f"StandardNewton failed to converge after {maxiter} iterations.")


class GminSteppingNewton(NonLinearSolver):
    """
    Continuation Method: Gmin-Stepping Decorator.
    
    If the base solver fails to converge, this wrapper iteratively injects
    a diagonal Gmin conductivity into the Jacobian to guide the solver 
    through highly non-linear regions.
    """
    
    def __init__(self, base_solver: NonLinearSolver):
        self.base_solver = base_solver
        
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None):
        try:
            # First, attempt to solve the pure system without Gmin injection
            return self.base_solver.solve_system(x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter)
        except NoConvergenceError:
            pass # Proceed to Gmin stepping
            
        x_curr = x0
        gmin_steps = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
        
        for gmin in gmin_steps:
            def eval_FJ_with_gmin(x):
                F, J = eval_FJ(x)
                # Inject Gmin into the diagonal of the Jacobian
                J_gmin = J + toolkit.eye(len(J)) * gmin
                # Add the leakage current to the RHS vector

                F_gmin = F + x * gmin
                return F_gmin, J_gmin
                
            try:
                x_curr, _ = self.base_solver.solve_system(x_curr, eval_FJ_with_gmin, toolkit, reltol, abstol, xtol, maxiter, limiter)
            except NoConvergenceError:
                raise NoConvergenceError(f"Gmin Stepping failed at gmin={gmin}")
                
        # Finally, solve the exact pure system using the guided initial guess
        return self.base_solver.solve_system(x_curr, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter)

class SchurCoupledNewton(NonLinearSolver):
    """
    Coupled Newton-Raphson Solver using the Schur Complement.
    
    Optimizes a dual-state (x, h) by partitioning the Jacobian into blocks.
    The callback `eval_FJ(x, h)` must return the full set of partitioned matrices:
    (F, J_x, J_h, E, E_x, E_h)
    """
    
    def solve_system(self, S0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None):
        x_curr, h_curr = S0
        
        for i in range(maxiter):
            F, J_x, J_h, E, E_x, E_h = eval_FJ(x_curr, h_curr)
            
            # Formulate the Schur Complement RHS
            import numpy as np
            # Note: We assume the caller handles reference node stripping within eval_FJ
            rhs = np.column_stack([-F, -J_h])
            try:
                dx_res = toolkit.linearsolver(J_x, rhs)
            except Exception:
                dx_res = np.zeros_like(rhs)
                
            dx_0 = dx_res[:, 0]
            dx_h = dx_res[:, 1]
            
            denom = toolkit.dot(E_x, dx_h) + E_h
            if abs(denom) < 1e-20:
                dh = 0.0
            else:
                dh = (-E - toolkit.dot(E_x, dx_0)) / denom
                
            dh = max(-0.5 * h_curr, min(2.0 * h_curr, dh))
            dx = dx_0 + dx_h * dh
            
            x_next = x_curr + dx
            h_next = h_curr + dh
            
            if limiter is not None:
                x_next = limiter(x_next, x_curr)
                # Re-compute dx after limiting
                dx = x_next - x_curr
                
            I_scale = toolkit.dot(abs(J_x), abs(x_next)) + abs(F)
            
            conv_x = toolkit.alltrue(abs(dx) < reltol * toolkit.maximum(abs(x_next), 1e-12) + xtol)
            conv_h = abs(dh) < 0.15 * h_curr
            conv_F = toolkit.alltrue(abs(F) < reltol * I_scale + abstol)
            
            if conv_x and conv_h and conv_F:
                return (x_next, h_next), i + 1
                
            x_curr = x_next
            h_curr = h_next
            
        raise NoConvergenceError(f"SchurCoupledNewton failed to converge after {maxiter} iterations.")
