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
            x_full = x + xdiff
            if limiter is not None:
                x_full = limiter(x_full, x)
                xdiff = x_full - x

            ## Backtracking line search -- F15 (doc/transient_review_260820.md).
            ## On failure this used to take the FULL step (the candidate the
            ## search had just proved worst), keep the residual of the PREVIOUS
            ## iterate, and test convergence with |alpha*dx| at alpha ~ 0.031
            ## against a step of |dx| -- a 32x-lenient test of a point whose
            ## residual was never evaluated.  Now the step taken is always the
            ## last one TRIED (the least-bad damped step on failure), its own
            ## residual travels with it, and the convergence test measures the
            ## step actually taken.
            while True:
                x_next = x + alpha * xdiff
                F_next, _ = eval_FJ(x_next)
                if toolkit.sum(abs(F_next)) <= F_norm * (1.0 - 1e-4 * alpha):
                    break                      # Armijo satisfied
                if alpha * 0.5 <= 0.05:
                    break                      # floor: keep the smallest tried
                alpha *= 0.5
            F = F_next

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


## Both continuation wrappers below forward EVERY solver option, `linsolver`
## included.  All six inner `solve_system` calls used to drop it, so a
## caller-selected linear-solver strategy silently reverted to the default the
## moment continuation engaged -- the third option-dropped-in-transit defect in
## this file's lineage, after `row_names` was retrofitted the same way
## (doc/transient_review_260820.md, F9).  If `solve_system` grows another
## option, bundle them into one object rather than extending this list again.
def _adaptive_conductance_ladder(solve_rung, solve_pure, x0,
                                 e_start=-2.0, e_end=-12.0, e_max=0.0,
                                 min_step=0.25, max_rungs=60,
                                 label='continuation'):
    """The adaptive rung driver both ladders share (owner request).

    Works in EXPONENT space, SPICE3-style: march g = 10**e from e_start
    down to e_end in steps of `step` decades, then finish with the PURE
    solve.  On a failed rung the step HALVES and the march resumes from
    the last converged exponent -- the fixed-decade schedule was measured
    stranding a rung whose solution sat too far from its predecessor's
    (the P18 unit-mock needed a 55-iteration walk across one decade gap).
    If the FIRST rung fails, escalate upward (more conductance, toward
    e_max) before refining.  A step below `min_step` decades or more than
    `max_rungs` rungs gives up honestly with the last failure chained.
    """
    x_curr = x0
    step = 2.0
    e = e_start
    last_e = None          # last converged exponent
    rungs = 0
    while True:
        if rungs >= max_rungs:
            raise NoConvergenceError(
                '%s: gave up after %d adaptive rungs' % (label, max_rungs))
        rungs += 1
        try:
            x_curr = solve_rung(x_curr, 10.0 ** e)
            last_e = e
            if e <= e_end:
                break
            e = max(e - step, e_end)
        except NoConvergenceError as exc:
            if last_e is None:
                ## Nothing has converged yet: escalate to a LARGER
                ## conductance first (an easier problem exists upward).
                if e < e_max:
                    e = min(e + step, e_max)
                    continue
                step = step / 2.0
                if step < min_step:
                    raise NoConvergenceError(
                        '%s failed even at g=10**%g with step refinement '
                        'exhausted: %s' % (label, e_max, exc)) from exc
                e = e_start
                continue
            step = step / 2.0
            if step < min_step:
                raise NoConvergenceError(
                    '%s: refinement exhausted below 10**%g (step < %g '
                    'decades): %s' % (label, last_e, min_step, exc)) from exc
            e = max(last_e - step, e_end)
    ## The zero rung: ONLY a pure converged solution may be returned.
    return solve_pure(x_curr)


class JunctionGminSteppingNewton(NonLinearSolver):
    """Continuation: junction-parallel gmin stepping -- P18's PRIMARY ladder.

    `gmin` proper (industry vocabulary): a conductance in parallel with each
    pn junction, ramped in decades and finished with a PURE solve -- the
    deformed circuit is a physical leaky diode, the linear subnetwork is
    untouched, and the homotopy tracks the physical branch.  Junction row
    pairs come from `pcnr_junctions`; a circuit without declared junctions
    passes straight through to the base solver.
    """

    def __init__(self, base_solver: NonLinearSolver, junction_rows):
        self.base_solver = base_solver
        self.junction_rows = junction_rows   # [(ra, rb), ...] or empty

    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol,
                     maxiter, limiter=None, scaler=None, row_names=None,
                     linsolver=None):
        try:
            return self.base_solver.solve_system(
                x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter,
                limiter, scaler, row_names=row_names, linsolver=linsolver)
        except NoConvergenceError:
            if not self.junction_rows:
                raise

        import numpy as _np

        def solve_rung(x_seed, g):
            def eval_FJ_with_gmin(x):
                F, J = eval_FJ(x)
                F = _np.array(F, dtype=float, copy=True)
                J = _np.array(J, dtype=float, copy=True)
                for ra, rb in self.junction_rows:
                    vj = x[ra] - x[rb]
                    F[ra] += g * vj
                    F[rb] -= g * vj
                    J[ra, ra] += g; J[ra, rb] -= g
                    J[rb, ra] -= g; J[rb, rb] += g
                return F, J
            x, _ = self.base_solver.solve_system(
                x_seed, eval_FJ_with_gmin, toolkit, reltol, abstol,
                xtol, maxiter, limiter, scaler, row_names=row_names,
                linsolver=linsolver)
            return x

        def solve_pure(x_seed):
            return self.base_solver.solve_system(
                x_seed, eval_FJ, toolkit, reltol, abstol, xtol, maxiter,
                limiter, scaler, row_names=row_names, linsolver=linsolver)

        return _adaptive_conductance_ladder(
            solve_rung, solve_pure, x0, label='junction-gmin stepping')


class GminSteppingNewton(NonLinearSolver):
    """
    Continuation Method: diagonal conductance stepping.

    INDUSTRY VOCABULARY NOTE (owner correction, P18): what this class ramps
    is **gshunt** -- a node-to-ground (here: full-diagonal) conductance --
    not `gmin`, which properly names the junction-parallel conductance
    (see JunctionGminSteppingNewton).  The class name is kept for
    compatibility; it is the SINGULAR-MATRIX rescue of the chain, tried
    after the physical junction ladder.

    If the base solver fails to converge, this wrapper iteratively injects
    a diagonal conductivity into the Jacobian to guide the solver 
    through highly non-linear regions.
    """
    
    def __init__(self, base_solver: NonLinearSolver):
        self.base_solver = base_solver
        
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, row_names=None, linsolver=None):
        try:
            # First, attempt to solve the pure system without Gmin injection
            return self.base_solver.solve_system(x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names, linsolver=linsolver)
        except NoConvergenceError:
            pass # Proceed to Gmin stepping
            
        def solve_rung(x_seed, g):
            def eval_FJ_with_gshunt(x):
                F, J = eval_FJ(x)
                return F + x * g, J + toolkit.eye(len(J)) * g
            ## STAGE 6's keep-the-inner-diagnosis rule now lives in the
            ## adaptive driver's chained raises.
            x, _ = self.base_solver.solve_system(
                x_seed, eval_FJ_with_gshunt, toolkit, reltol, abstol,
                xtol, maxiter, limiter, scaler, row_names=row_names,
                linsolver=linsolver)
            return x

        def solve_pure(x_seed):
            return self.base_solver.solve_system(
                x_seed, eval_FJ, toolkit, reltol, abstol, xtol, maxiter,
                limiter, scaler, row_names=row_names, linsolver=linsolver)

        return _adaptive_conductance_ladder(
            solve_rung, solve_pure, x0, e_start=-3.0,
            label='gshunt stepping')

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
            return self.base_solver.solve_system(x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names, linsolver=linsolver)
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
                x_curr, _ = self.base_solver.solve_system(x_curr, eval_FJ_with_source, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names, linsolver=linsolver)
            except NoConvergenceError as e:
                ## STAGE 6 -- keep the inner diagnosis; see the Gmin note above.
                raise NoConvergenceError(
                    'Source Stepping failed at lambda=%s: %s' % (lambda_, e)) from e
                
        return self.base_solver.solve_system(x_curr, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler, row_names=row_names, linsolver=linsolver)


class PseudoTransientNewton(NonLinearSolver):
    """Continuation: pseudo-transient (Psi-tc) -- P25, the chain's LAST rung.

    Embeds F(x) = 0 in the pseudo-ODE dx/dtau = -F(x) and takes
    backward-Euler pseudo-time steps to steady state: each step solves
    F(x) + (1/delta)*(x - x_k) = 0 with Jacobian J + I/delta, anchored at
    the PREVIOUS iterate.  With g = 1/delta this is exactly a conductance
    rung with a MOVING anchor, so the adaptive exponent-space driver is
    reused verbatim: g marches 1 -> 1e-12 (delta grows 1 s -> 1e12 s, the
    classical exponential pseudo-timestep growth), halving on failure,
    escalating to heavier damping (e_max = +6) if even the first step
    fails, and finishing PURE -- only a g = 0 converged solution may be
    returned (the P22 rule).

    What the moving anchor buys over the zero-anchored gshunt ladder,
    both measured at landing (doc/backend_parity_260821.md, P25):
    the homotopy follows the pseudo-transient trajectory FROM THE SEED,
    so on a multi-stable system it lands in the caller's basin -- on
    F = x^3 - x from x0 = 0.9 the gshunt path (x^3 + (g-1)x has only the
    0 root for g > 1) commits to 0 while Psi-tc lands on +1 -- and on the
    classic Newton 2-cycle cubic x^3 - 2x + 2 it converges in 3.4x fewer
    base-solver calls (47 vs 158) because no rung has to walk back from
    the alien zero-anchored deformation.

    `rung_solver` (default: `base_solver`) solves the deformed pseudo
    steps and the final pure solve.  Pass it when `base_solver` is a full
    continuation chain: a chain's own ladders would otherwise re-engage
    inside every pseudo step -- and SourceSteppingNewton's rungs rebuild
    F from its callback WITHOUT the pseudo term (the option-dropped-in-
    transit class), so the chain must never solve the deformed system.
    """

    def __init__(self, base_solver: NonLinearSolver, rung_solver=None):
        self.base_solver = base_solver
        self.rung_solver = rung_solver if rung_solver is not None \
            else base_solver

    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol,
                     maxiter, limiter=None, scaler=None, row_names=None,
                     linsolver=None):
        try:
            return self.base_solver.solve_system(
                x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter,
                limiter, scaler, row_names=row_names, linsolver=linsolver)
        except NoConvergenceError:
            pass  # proceed to pseudo-transient stepping

        def solve_rung(x_seed, g):
            def eval_FJ_ptc(x):
                F, J = eval_FJ(x)
                return (F + g * (x - x_seed),
                        J + toolkit.eye(len(J)) * g)
            x, _ = self.rung_solver.solve_system(
                x_seed, eval_FJ_ptc, toolkit, reltol, abstol,
                xtol, maxiter, limiter, scaler, row_names=row_names,
                linsolver=linsolver)
            return x

        def solve_pure(x_seed):
            return self.rung_solver.solve_system(
                x_seed, eval_FJ, toolkit, reltol, abstol, xtol, maxiter,
                limiter, scaler, row_names=row_names, linsolver=linsolver)

        return _adaptive_conductance_ladder(
            solve_rung, solve_pure, x0, e_start=0.0, e_end=-12.0,
            e_max=6.0, label='pseudo-transient continuation')


class SchurCoupledNewton(NonLinearSolver):
    """
    Coupled Newton-Raphson Solver using the Schur Complement.
    
    Optimizes a dual-state (x, h) by partitioning the Jacobian into blocks.
    The callback `eval_FJ(x, h)` must return the full set of partitioned matrices:
    (F, J_x, J_h, E, E_x, E_h)
    """
    
    def solve_system(self, S0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, row_names=None, linsolver=None, hmin=0.0):
        """``hmin`` floors the solved step size.

        STAGE 12B-0.  Without it this iteration has no globalisation at all: `dh`
        is clamped per iteration to ``[-h/2, 2h]``, so a step that keeps failing
        its LTE test simply halves for every one of ``maxiter`` iterations and
        arrives at ``h * 2**-maxiter``.  Measured on the stress circuits, that is
        exactly what happens -- the solved step reached 1.6e-25, 4.5e-38 and
        2.5e-45 s inside a single time point, and the run then made no progress.

        The floor matters more than a step-size guess would suggest, because the
        LTE estimate itself stops being meaningful first: the divided differences
        that produce it carry a ``1/h`` factor, so as ``h`` shrinks they amplify
        rounding in the charges until the "error" being solved against is
        cancellation noise.  Newton is then chasing a quantity that grows as the
        step shrinks, which is a runaway rather than a convergence.
        """
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
            ## Floor the SOLVED step, and take `dh` from the floored result so the
            ## solution update `dx` stays consistent with the step it belongs to.
            ## Clamping `h_next` afterwards while leaving `dh` alone would apply a
            ## correction for a step that was not taken.
            h_next = max(h_curr + dh, hmin)
            dh = h_next - h_curr
            dx = dx_0 + dx_h * dh

            x_next = x_curr + dx
            
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

## JAXNewtonSolver WAS DELETED HERE (F15, doc/transient_review_260820.md).
## It was a second, unfixed copy of the stale-residual defect stage 9(e)
## repaired in jaxtransient's own loop -- its converged test read the residual
## of the iterate BEFORE the final update, so hitting maxiter raised even when
## the returned point was converged -- and nothing in the package or the tests
## ever constructed it.  The production JAX Newton is newton_inner_loop in
## jaxtransient.py, which carries the per-row criteria (F6(b)).
