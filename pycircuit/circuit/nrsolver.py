from abc import ABC, abstractmethod
import numpy as np

from pycircuit.circuit.analysis import NoConvergenceError, SingularMatrix
from pycircuit.circuit.scaler import NoneScaler

def _structural_singularity(J, row_names, toolkit):
    """Name a structurally singular row or column, or return None.

    STAGE 6(a).  A structurally singular matrix is not something continuation can
    repair: gmin-stepping and source-stepping both perturb a system that has a
    solution, and neither adds the missing equation.  Reporting it as
    "Source Stepping failed at lambda=0.0" -- which is what three layers of
    re-wrapping used to produce -- tells the user about the last thing that was
    tried rather than the first thing that was wrong.

    Two shapes, and they mean different things to the person reading the message:

      * an all-zero COLUMN means the unknown appears in no equation.  For a node
        that is "nothing constrains this voltage" -- the floating-node case.
      * an all-zero ROW means the equation constrains nothing.

    Detected on the assembled Jacobian rather than from an LU pivot because it is
    exact, needs no factorisation, and cannot be confused with a merely
    ill-conditioned matrix -- which continuation genuinely can help with and which
    must therefore keep its existing path.
    """
    try:
        A = abs(np.asarray(J, dtype=float))
    except (TypeError, ValueError):
        return None            # symbolic: no numeric structure to inspect
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        return None

    def _name(i):
        if row_names is not None and 0 <= i < len(row_names):
            return "'%s'" % row_names[i]
        return 'row %d' % i

    col_dead = np.flatnonzero(A.sum(axis=0) == 0.0)
    if len(col_dead):
        i = int(col_dead[0])
        return ('%s appears in no equation, so nothing determines it -- for a node '
                'that means no DC path to ground (add a resistor, or use uic=True '
                'to skip the operating point)' % _name(i))

    row_dead = np.flatnonzero(A.sum(axis=1) == 0.0)
    if len(row_dead):
        i = int(row_dead[0])
        return ('the equation at %s constrains nothing -- every entry in its row '
                'is zero' % _name(i))
    return None


def _worst_row_report(F, xdiff, x_next, x, reltol, abstol, xtol, row_names, toolkit):
    """Which unknown failed to converge, by how much, and against what tolerance.

    STAGE 6(b).  Everything here was already in scope at the point the old message
    was raised; it said "DampedNewton failed to converge after N iterations" and
    threw all of it away.
    """
    def _name(i):
        if row_names is not None and 0 <= i < len(row_names):
            return "'%s'" % row_names[i]
        return 'row %d' % i

    try:
        F = np.asarray(F, dtype=float)
        xd = abs(np.asarray(xdiff, dtype=float))
        xn = np.asarray(x_next, dtype=float)
        xo = np.asarray(x, dtype=float)
        ab = np.broadcast_to(np.asarray(abstol, dtype=float), F.shape)
        xt = np.broadcast_to(np.asarray(xtol, dtype=float), F.shape)
    except (TypeError, ValueError):
        return ''

    ## Normalised misses: >1 means that row failed its test.  Reporting the
    ## normalised value rather than the raw residual is what makes the number
    ## actionable -- "2.3x its tolerance" says how far off, "4.7e-9 A" does not.
    f_tol = reltol * abs(F) + ab
    x_tol = reltol * np.maximum(abs(xn), abs(xo)) + xt
    with np.errstate(divide='ignore', invalid='ignore'):
        f_miss = np.where(f_tol > 0, abs(F) / f_tol, 0.0)
        x_miss = np.where(x_tol > 0, xd / x_tol, 0.0)

    parts = []
    if f_miss.size and np.nanmax(f_miss) > 1.0:
        i = int(np.nanargmax(f_miss))
        parts.append('residual worst at %s: |f| = %.4g against a tolerance of '
                     '%.4g (%.3gx over)' % (_name(i), abs(F[i]), f_tol[i], f_miss[i]))
    if x_miss.size and np.nanmax(x_miss) > 1.0:
        i = int(np.nanargmax(x_miss))
        parts.append('update worst at %s: |dx| = %.4g against a tolerance of '
                     '%.4g (%.3gx over)' % (_name(i), xd[i], x_tol[i], x_miss[i]))
    return '; '.join(parts)


class NonLinearSolver(ABC):
    """
    Abstract Base Class for Non-Linear Algebraic Solvers.
    
    This interface isolates the Newton-Raphson iteration math and continuation 
    methods (like Gmin-stepping). It is entirely agnostic to whether it is 
    solving a DC operating point or a Transient time-step.
    """
    
    @abstractmethod
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, row_names=None, linsolver=None):
        pass


class StandardNewton(NonLinearSolver):
    """
    Standard Full-Matrix Newton-Raphson Solver.
    Used natively by DCAnalysis and standard adaptive Transient simulations.
    """
    
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, row_names=None, linsolver=None):
        x = x0
        if scaler is None:
            scaler = NoneScaler()
            
        for i in range(maxiter):
            F, J = eval_FJ(x)
            
            # --- EQUILIBRATION SCALING ---
            J_s, F_s, s_vec = scaler.scale(J, F, toolkit)
            
            try:
                ## STAGE 7b.  `linsolver is None` keeps the exact call this made
                ## before, so nothing moves unless a caller asks for a strategy.
                if linsolver is None:
                    xdiff = toolkit.linearsolver(J_s, -F_s)
                else:
                    xdiff = linsolver.solve(J_s, -F_s, toolkit)
                xdiff = scaler.unscale_solution(xdiff, s_vec, toolkit)
            except Exception as e:
                ## STAGE 6(a).  See the matching note in DampedNewton -- and note
                ## that THIS is the default solver, so this is the copy that
                ## actually runs unless a caller asks for another.
                why = _structural_singularity(J, row_names, toolkit)
                if why is not None:
                    raise SingularMatrix('singular Jacobian: %s' % why) from e
                raise NoConvergenceError(f"Singular Jacobian: {str(e)}")
            
            x_next = x + xdiff
            if limiter is not None:
                x_next = limiter(x_next, x)
                # Recompute xdiff based on limited step
                xdiff = x_next - x
                
            # --- CONVERGENCE CRITERIA ---
            # To declare victory, two conditions must be met:
            # 1. Voltage/Current Update (conv_x): The state vector (x) must stop changing.
            #    We check if the difference (xdiff) is smaller than a relative tolerance
            #    based on the current values, plus an absolute floor (xtol) to prevent 
            #    division-by-zero for zero-crossing nodes.
            
            # 2. KCL Residual (conv_f): The sum of currents at each node (F) must be near zero.
            #    However, 1mA error is huge for a micro-power circuit, but tiny for a 100A power 
            #    supply. Therefore, we scale the tolerance dynamically (I_scale) by looking at the 
            #    absolute magnitude of currents flowing into the node.
            I_scale = toolkit.dot(abs(J), abs(x_next)) + abs(F)
            
            conv_x = toolkit.alltrue(abs(xdiff) < reltol * toolkit.maximum(abs(x_next), abs(x)) + xtol)
            conv_f = toolkit.alltrue(abs(F) < reltol * I_scale + abstol)
            
            if conv_x and conv_f:
                return x_next, i + 1
                
            x = x_next
            
        ## STAGE 6(b).
        detail = _worst_row_report(F, xdiff, x_next, x, reltol, abstol, xtol,
                                   row_names, toolkit)
        raise NoConvergenceError(
            'StandardNewton failed to converge after %d iterations%s'
            % (maxiter, ('; ' + detail) if detail else ''))


class DampedNewton(NonLinearSolver):
    """
    Damped Newton-Raphson Solver (Backtracking Line Search).
    If a full step causes the residual to increase, the step size (alpha) is halved.
    """
    
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, row_names=None, linsolver=None):
        x = x0
        if scaler is None:
            scaler = NoneScaler()
            
        for i in range(maxiter):
            F, J = eval_FJ(x)
            F_norm = toolkit.sum(abs(F))
            
            J_s, F_s, s_vec = scaler.scale(J, F, toolkit)
            
            try:
                ## STAGE 7b.  `linsolver is None` keeps the exact call this made
                ## before, so nothing moves unless a caller asks for a strategy.
                if linsolver is None:
                    xdiff = toolkit.linearsolver(J_s, -F_s)
                else:
                    xdiff = linsolver.solve(J_s, -F_s, toolkit)
                xdiff = scaler.unscale_solution(xdiff, s_vec, toolkit)
            except Exception as e:
                ## STAGE 6(a) -- CLASSIFY BEFORE THE CONTINUATION LAYERS SEE IT.
                ## A structural singularity is not something gmin- or
                ## source-stepping can repair, so it is raised as SingularMatrix,
                ## which is NOT a NoConvergenceError and therefore passes straight
                ## through both decorators instead of being re-wrapped into
                ## "Source Stepping failed at lambda=0.0".
                why = _structural_singularity(J, row_names, toolkit)
                if why is not None:
                    raise SingularMatrix('singular Jacobian: %s' % why) from e
                raise NoConvergenceError(f"Singular Jacobian: {str(e)}")
            
            alpha = 1.0
            x_next = x + xdiff
            if limiter is not None:
                x_next = limiter(x_next, x)
                xdiff = x_next - x
                
            # Backtracking line search
            while alpha > 0.05:
                x_test = x + alpha * xdiff
                F_test, _ = eval_FJ(x_test)
                if toolkit.sum(abs(F_test)) <= F_norm * (1.0 - 1e-4 * alpha):
                    # Step accepted
                    x_next = x_test
                    F = F_test
                    break
                alpha *= 0.5
            
            I_scale = toolkit.dot(abs(J), abs(x_next)) + abs(F)
            
            conv_x = toolkit.alltrue(abs(alpha * xdiff) < reltol * toolkit.maximum(abs(x_next), abs(x)) + xtol)
            conv_f = toolkit.alltrue(abs(F) < reltol * I_scale + abstol)
            
            if conv_x and conv_f:
                return x_next, i + 1
                
            x = x_next
            
        ## STAGE 6(b) -- SAY WHICH UNKNOWN FAILED, BY HOW MUCH, AND AGAINST WHAT.
        ## All of this was already in scope when the old message threw it away.
        detail = _worst_row_report(F, xdiff, x_next, x, reltol, abstol, xtol,
                                   row_names, toolkit)
        raise NoConvergenceError(
            'DampedNewton failed to converge after %d iterations%s'
            % (maxiter, ('; ' + detail) if detail else ''))


class GminSteppingNewton(NonLinearSolver):
    """
    Continuation Method: Gmin-Stepping Decorator.
    
    If the base solver fails to converge, this wrapper iteratively injects
    a diagonal Gmin conductivity into the Jacobian to guide the solver 
    through highly non-linear regions.
    """
    
    def __init__(self, base_solver: NonLinearSolver):
        self.base_solver = base_solver
        
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, row_names=None, linsolver=None):
        try:
            # First, attempt to solve the pure system without Gmin injection
            return self.base_solver.solve_system(x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names)
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
                x_curr, _ = self.base_solver.solve_system(x_curr, eval_FJ_with_gmin, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names)
            except NoConvergenceError as e:
                ## STAGE 6 -- KEEP THE INNER DIAGNOSIS.  This used to discard it,
                ## so the caller was told which continuation rung was last tried
                ## and nothing about what actually failed.  The rung is useful
                ## context; it is not a substitute for the cause.
                raise NoConvergenceError(
                    'Gmin Stepping failed at gmin=%s: %s' % (gmin, e)) from e
                
        # Finally, solve the exact pure system using the guided initial guess
        return self.base_solver.solve_system(x_curr, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names)

class SourceSteppingNewton(NonLinearSolver):
    """
    Continuation Method: Source-Stepping Decorator.
    
    If the base solver fails to converge, this wrapper iteratively scales 
    the independent sources from 0 to 1 to guide the solver.
    """
    def __init__(self, base_solver: NonLinearSolver, source_callback):
        self.base_solver = base_solver
        self.source_callback = source_callback
        
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, row_names=None, linsolver=None):
        try:
            # Note: eval_FJ natively evaluates sources at 1.0
            return self.base_solver.solve_system(x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names)
        except NoConvergenceError:
            pass # Proceed to source stepping
            
        x_curr = x0
        lambdas = [0.0, 1e-2, 1e-1, 1.0]
        
        for lambda_ in lambdas:
            def eval_FJ_with_source(x):
                # The callback MUST scale the source term specifically.
                # In PyCircuit, F = i(x) + lambda_ * u(0)
                # But evaluating this requires access to the circuit elements.
                return self.source_callback(x, lambda_)
                
            try:
                x_curr, _ = self.base_solver.solve_system(x_curr, eval_FJ_with_source, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names)
            except NoConvergenceError as e:
                ## STAGE 6 -- keep the inner diagnosis; see the Gmin note above.
                raise NoConvergenceError(
                    'Source Stepping failed at lambda=%s: %s' % (lambda_, e)) from e
                
        return self.base_solver.solve_system(x_curr, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names)

class SchurCoupledNewton(NonLinearSolver):
    """
    Coupled Newton-Raphson Solver using the Schur Complement.
    
    Optimizes a dual-state (x, h) by partitioning the Jacobian into blocks.
    The callback `eval_FJ(x, h)` must return the full set of partitioned matrices:
    (F, J_x, J_h, E, E_x, E_h)
    """
    
    def solve_system(self, S0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, row_names=None, linsolver=None):
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

class JAXNewtonSolver(NonLinearSolver):
    """
    A Pure JAX Newton-Raphson Solver.
    Uses `jax.lax.while_loop` to compile the entire solution process into a single GPU kernel.
    Only supports Dense matrices and full-JAX compatible circuits (no Python fallback elements).
    """
    
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, row_names=None, linsolver=None):
        if not toolkit.supports('autodiff'):
            raise ValueError("JAXNewtonSolver requires the JAX toolkit (_jaxtoolkit.py).")
            
        import jax
        import jax.numpy as jnp
        
        # We need a fallback or check to ensure eval_FJ is fully jittable.
        # But we'll try to just JIT it and if it fails, it throws a Tracer error.
        
        def cond_fun(state):
            x, xdiff, F_norm, iters = state
            
            # Convergence conditions
            conv_f = F_norm < abstol
            conv_x = jnp.all(jnp.abs(xdiff) < reltol * jnp.abs(x) + xtol)
            
            converged = jnp.logical_and(conv_f, conv_x)
            not_converged = jnp.logical_not(converged)
            under_max = iters < maxiter
            
            return jnp.logical_and(not_converged, under_max)
            
        def body_fun(state):
            x, _, _, iters = state
            
            F, J = eval_FJ(x)
            F_norm = jnp.sum(jnp.abs(F))
            
            # Simple dense linear solve (no scaler implementation for now)
            xdiff = jnp.linalg.solve(J, -F)
            
            x_next = x + xdiff
            if limiter is not None:
                x_next = limiter(x_next, x)
                xdiff = x_next - x
                
            return (x_next, xdiff, F_norm, iters + 1)
            
        # Compile the while loop
        @jax.jit
        def compiled_nr(x_init):
            # Evaluate once to get initial F_norm and seed the loop state
            F0, J0 = eval_FJ(x_init)
            F_norm0 = jnp.sum(jnp.abs(F0))
            xdiff0 = jnp.linalg.solve(J0, -F0)
            
            x_next = x_init + xdiff0
            if limiter is not None:
                x_next = limiter(x_next, x_init)
                xdiff0 = x_next - x_init
                
            initial_state = (x_next, xdiff0, F_norm0, 1)
            
            # Run the compiled loop
            final_state = jax.lax.while_loop(cond_fun, body_fun, initial_state)
            return final_state
            
        final_x, final_xdiff, final_F_norm, final_iters = compiled_nr(x0)
        
        # In pure JAX, we can't raise Python exceptions inside the compiled loop easily.
        # We check convergence after the loop finishes.
        conv_f = final_F_norm < abstol
        conv_x = jnp.all(jnp.abs(final_xdiff) < reltol * jnp.abs(final_x) + xtol)
        if not (conv_f and conv_x):
            raise NoConvergenceError(f"JAXNewtonSolver failed to converge after {maxiter} iterations.")
            
        return final_x, int(final_iters)

