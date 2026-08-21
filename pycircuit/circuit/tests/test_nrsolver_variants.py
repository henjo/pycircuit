import pytest
import numpy as np
from pycircuit.circuit.nrsolver import NonLinearSolver, StandardNewton, DampedNewton, GminSteppingNewton, SourceSteppingNewton, NoConvergenceError
from pycircuit.circuit.scaler import NoneScaler, RowMaxScaler, RowL2Scaler, SinkhornKnoppScaler
from pycircuit.circuit import numeric

class MockFailingSolver(NonLinearSolver):
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None):
        try:
            eval_FJ(x0)
        except Exception:
            pass
        raise NoConvergenceError("Mock failure")

class MockSourceSteppingBaseSolver(NonLinearSolver):
    ## `**kwargs` since stage 6, which added `row_names` to the `NonLinearSolver`
    ## contract so a diagnostic can say "'net3'" instead of "row 4".  Absorbing
    ## unknown keywords is what any out-of-tree solver should do; the decorators
    ## pass the new argument by NAME, so a subclass that does neither fails loudly
    ## with a TypeError rather than silently binding it to the wrong parameter.
    def __init__(self):
        self.call_count = 0
        
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, **kwargs):
        self.call_count += 1
        # The first call is the direct solve (no source stepping), we want it to fail
        if self.call_count == 1:
            raise NoConvergenceError("Mock failure on pure solve")
            
        try:
            eval_FJ(x0)
        except Exception:
            pass
        return np.array([1.0, 2.0]), 5

def test_source_stepping_trigger():
    lambdas_seen = []
    
    def mock_eval_FJ(x):
        return x, np.eye(len(x))
        
    def mock_source_callback(x, lam):
        lambdas_seen.append(lam)
        return x, np.eye(len(x))

    succeeding_base = MockSourceSteppingBaseSolver()
    source_stepper = SourceSteppingNewton(succeeding_base, mock_source_callback)
    
    # It will succeed, so no exception is raised
    source_stepper.solve_system(np.zeros(2), mock_eval_FJ, numeric, 1e-3, 1e-6, 1e-6, 10, None)
    
    assert lambdas_seen == [0.0, 1e-2, 1e-1, 1.0]

def test_gmin_stepping_trigger():
    gmins_seen = []
    
    def mock_eval_FJ(x):
        return x, np.eye(len(x))

    class MockFailingGmin(NonLinearSolver):
        def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None, **kwargs):
            try:
                F, J = eval_FJ(x0)
                gmin = J[0,0] - 1.0
                if gmin > 0:
                    gmins_seen.append(gmin)
            except Exception:
                pass
            raise NoConvergenceError("Mock failure")

    failing_base = MockFailingGmin()
    gmin_stepper = GminSteppingNewton(failing_base)
    
    with pytest.raises(NoConvergenceError):
        gmin_stepper.solve_system(np.zeros(2), mock_eval_FJ, numeric, 1e-3, 1e-6, 1e-6, 10, None)
        
    assert len(gmins_seen) > 0
    assert np.isclose(gmins_seen[0], 1e-3)

def test_schur_coupled_newton():
    from pycircuit.circuit.nrsolver import SchurCoupledNewton
    
    # System:
    # F(x, h) = x + h - 5 = 0
    # E(x, h) = 2x - 3h = 0
    # Solution: x = 3, h = 2
    
    def eval_FJ(x, h):
        F = np.array([x[0] + h - 5.0])
        J_x = np.array([[1.0]])
        J_h = np.array([1.0])
        E = 2.0 * x[0] - 3.0 * h
        E_x = np.array([2.0])
        E_h = -3.0
        return F, J_x, J_h, E, E_x, E_h
        
    solver = SchurCoupledNewton()
    x0 = np.array([0.0])
    h0 = 0.5
    
    # Needs a few iterations to scale dh if limited, but linear system should solve fast.
    (x_res, h_res), iters = solver.solve_system((x0, h0), eval_FJ, numeric, 1e-3, 1e-6, 1e-6, 20)
    
    assert np.isclose(x_res[0], 3.0)
    assert np.isclose(h_res, 2.0)

def test_damped_newton_backtracking():
    solver = DampedNewton()
    
    # x^2 - 4 = 0 (roots at -2, +2)
    # F = x^2 - 4
    # J = 2x
    # Standard Newton from x0=100 takes x1 = 100 - (10000-4)/200 = 50.02
    def eval_FJ(x):
        return np.array([x[0]**2 - 4.0]), np.array([[2.0 * x[0]]])
        
    x0 = np.array([100.0])
    x_res, iters = solver.solve_system(x0, eval_FJ, numeric, 1e-4, 1e-12, 1e-12, 100)
    
    assert np.isclose(x_res[0], 2.0)

@pytest.mark.parametrize("scaler_class", [NoneScaler, RowMaxScaler, RowL2Scaler, SinkhornKnoppScaler])
def test_scalers_linear_system(scaler_class):
    # Create a highly ill-conditioned linear system J * x = -F
    # Row 1 is huge, Row 2 is tiny
    J = np.array([[1e9, 2e9], [3e-9, 4e-9]])
    # True solution: x = [1, -1]
    # -F = J * x = [ -1e9, -1e-9 ] => F = [ 1e9, 1e-9 ]
    F = np.array([1e9, 1e-9])
    
    scaler = scaler_class()
    J_s, F_s, c_scale = scaler.scale(J, F, numeric)
    
    dx_s = numeric.linearsolver(J_s, -F_s)
    dx = scaler.unscale_solution(dx_s, c_scale, numeric)
    
    assert np.allclose(dx, [1.0, -1.0])

def test_sinkhorn_knopp_matrix_stochasticity():
    scaler = SinkhornKnoppScaler(max_iter=10)
    
    J = np.array([[1.0, 100.0, 0.5],
                  [10.0, 1.0, 50.0],
                  [0.1, 0.05, 1.0]])
    F = np.array([1.0, 1.0, 1.0])
    
    J_s, _, _ = scaler.scale(J, F, numeric)
    
    # Check if the matrix is doubly stochastic (row/col sums near 1)
    row_sums = np.sum(np.abs(J_s), axis=1)
    col_sums = np.sum(np.abs(J_s), axis=0)
    
    assert np.allclose(row_sums, 1.0, atol=1e-2)
    assert np.allclose(col_sums, 1.0, atol=1e-2)


def test_junction_gmin_stepping_is_the_primary_ladder():
    """P18: junction-parallel gmin (the proper `gmin`, owner vocabulary)
    rescues the classic junction slam -- an exponential whose pure Newton
    step from a zero seed is ~2.5e7 (divergence) -- by dominating the
    early rungs and walking the seed in; the final solve is PURE (the
    zero-rung rule: ladder residue in a committed point is the P22
    inconsistency-floor trap).  A circuit without junctions passes
    straight through untouched."""
    import numpy as np
    from pycircuit.circuit.nrsolver import JunctionGminSteppingNewton

    seen_J00 = []

    class MiniNewton(NonLinearSolver):
        def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol,
                         maxiter, limiter=None, scaler=None, row_names=None,
                         linsolver=None):
            x = np.asarray(x0, float).copy()
            ## 120 iterations: a Newton walk on exp(40 v) advances ~25 mV
            ## per iteration, and a rung overshoot walks back from the limexp
            ## region (~55 iterations) before descending.  The
            ## divergence guard at 1e3 catches the PURE first step (2.5e7)
            ## without tripping on the walk.
            for _ in range(120):
                F, J = eval_FJ(x)
                seen_J00.append(float(J[0, 0]))
                dx = np.linalg.solve(J, -F)
                x = x + dx
                if np.max(np.abs(x)) > 1e3:
                    raise NoConvergenceError('diverged')
                if np.max(np.abs(dx)) < 1e-9:
                    return x, 1
            raise NoConvergenceError('maxiter')

    def eval_FJ(x):
        ## The junction slam: i = IS*(limexp(40*v)-1) driven by 1 mA.
        ## limexp (linearized above arg=80) is what real device models do
        ## -- a RAW exp overflows to inf/NaN on a rung overshoot and
        ## Newton loops on NaN; the clamp keeps the walk-back finite.
        ## Pure Newton from v=0: dx = 1e-3/4e-11 = 2.5e7 -> diverges.
        a = 40.0 * x[0]
        if a > 80.0:
            e = np.exp(80.0) * (1.0 + (a - 80.0))
            de = np.exp(80.0) * 40.0
        else:
            e = np.exp(a)
            de = 40.0 * e
        F = np.array([1e-12 * (e - 1.0) - 1e-3, x[1]])
        J = np.array([[1e-12 * de, 0.0], [0.0, 1.0]])
        return F, J

    x_star = np.log(1e9 + 1.0) / 40.0
    solver = JunctionGminSteppingNewton(MiniNewton(), [(0, 1)])
    x, _ = solver.solve_system(np.zeros(2), eval_FJ, None, 1e-4,
                               np.full(2, 1e-12), np.full(2, 1e-12), 120)
    assert abs(x[0] - x_star) < 1e-6, (x[0], x_star)
    ## The FINAL solve saw the pure model Jacobian (no injected gmin).
    assert seen_J00[-1] == pytest.approx(4e-11 * np.exp(40.0 * x[0]),
                                         rel=1e-9)
    ## Early rungs were genuinely dominated by the injection.
    assert max(seen_J00) >= 1e-2

    ## ADAPTIVE REFINEMENT (owner request): at an iteration budget of 20,
    ## a fixed decade gap's ~55-iteration walk-back cannot fit in any
    ## single rung -- the halving driver inserts intermediate rungs and
    ## still lands at machine precision (measured err 1.1e-16, ~165 base
    ## calls of refinement work at landing).
    tight = []

    class TightNewton(NonLinearSolver):
        def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol,
                         maxiter, limiter=None, scaler=None, row_names=None,
                         linsolver=None):
            x = np.asarray(x0, float).copy()
            for _ in range(20):
                F, J = eval_FJ(x)
                tight.append(1)
                dx = np.linalg.solve(J, -F)
                x = x + dx
                if np.max(np.abs(x)) > 1e3:
                    raise NoConvergenceError('diverged')
                if np.max(np.abs(dx)) < 1e-9:
                    return x, 1
            raise NoConvergenceError('maxiter')

    xa, _ = JunctionGminSteppingNewton(TightNewton(), [(0, 1)]).solve_system(
        np.zeros(2), eval_FJ, None, 1e-4, np.full(2, 1e-12),
        np.full(2, 1e-12), 20)
    assert abs(xa[0] - x_star) < 1e-6
    ## More base work than any fixed 7-rung schedule could spend within
    ## the budget -- the signature of inserted rungs.
    assert len(tight) > 7 * 20 * 0.7

    ## No junctions -> pure passthrough: the base fails, the ladder must
    ## re-raise without attempting rungs.
    n_before = len(seen_J00)
    with pytest.raises(NoConvergenceError):
        JunctionGminSteppingNewton(MiniNewton(), []).solve_system(
            np.zeros(2), eval_FJ, None, 1e-4, np.full(2, 1e-12),
            np.full(2, 1e-12), 8)
    assert len(seen_J00) - n_before <= 120


def test_pseudo_transient_continuation_is_the_chains_last_rung():
    """P25: pseudo-transient (Psi-tc) continuation -- backward-Euler
    pseudo-time steps F + (1/delta)*(x - x_k) = 0 through the shared
    adaptive driver, with the anchor MOVING to each rung's seed.  Gated
    here, all measured at landing:

    (1) it rescues the classic Newton 2-cycle cubic x^3 - 2x + 2 from
        x0 = 0 (plain Newton cycles 0 -> 1 -> 0 forever), landing on the
        real root at machine precision with a PURE final solve;
    (2) it does so CHEAPER than the zero-anchored gshunt ladder measured
        on the same problem (47 vs 158 base calls at landing) -- no rung
        walks back from an alien deformation;
    (3) the moving anchor is basin-respecting: on tristable F = x^3 - x
        every seed lands on ITS OWN basin's root.  (The zero-anchor
        hazard was measured but is NOT gated cross-solver: gshunt's
        deformation x^3 + (g-1)x has only the 0 root for g > 1, and a
        ladder FORCED to start at g = 1 commits every seed to 0 -- but
        the shipped GminSteppingNewton starts at g = 1e-3, whose weak
        rung already converges from these seeds, so the shipped class
        does not exhibit the hazard on this problem);
    (4) `rung_solver` keeps the deformed solves OFF a failed chain (whose
        source-stepping rungs would rebuild F without the pseudo term)."""
    from pycircuit.circuit.nrsolver import PseudoTransientNewton

    calls = []

    class MiniNewton(NonLinearSolver):
        def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol,
                         maxiter, limiter=None, scaler=None, row_names=None,
                         linsolver=None):
            x = np.asarray(x0, float).copy()
            for _ in range(60):
                F, J = eval_FJ(x)
                calls.append(float(J[0, 0]))
                dx = np.linalg.solve(J, -F)
                x = x + dx
                if np.max(np.abs(x)) > 1e6:
                    raise NoConvergenceError('diverged')
                if np.max(np.abs(dx)) < 1e-12:
                    return x, 1
            raise NoConvergenceError('maxiter (cycle)')

    def eval_cubic(x):
        ## Newton from 0 cycles: F(0)=2, J=-2 -> x=1; F(1)=1, J=1 -> x=0.
        return (np.array([x[0]**3 - 2.0*x[0] + 2.0]),
                np.array([[3.0*x[0]**2 - 2.0]]))

    root = -1.7692923542386314
    tols = (numeric, 1e-6, np.full(1, 1e-12), np.full(1, 1e-12), 60)

    with pytest.raises(NoConvergenceError):
        MiniNewton().solve_system(np.zeros(1), eval_cubic, *tols)

    calls.clear()
    x, _ = PseudoTransientNewton(MiniNewton()).solve_system(
        np.zeros(1), eval_cubic, *tols)
    assert x[0] == pytest.approx(root, abs=1e-12)
    n_ptc = len(calls)
    ## The FINAL solve saw the pure Jacobian: no pseudo conductance left.
    assert calls[-1] == pytest.approx(3.0*x[0]**2 - 2.0, rel=1e-12)

    calls.clear()
    xg, _ = GminSteppingNewton(MiniNewton()).solve_system(
        np.zeros(1), eval_cubic, *tols)
    assert xg[0] == pytest.approx(root, abs=1e-12)
    ## (2): both rescue this one, but the moving anchor pays measurably
    ## less (47 vs 158 at landing; asserted loosely to survive tuning).
    assert n_ptc < len(calls)

    ## (3): the basin discriminator.  FailFirst makes the pure attempt
    ## fail so both ladders genuinely engage from x0 = 0.9.
    class FailFirst(NonLinearSolver):
        def __init__(self, inner):
            self.inner = inner
            self.n = 0
        def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol,
                         maxiter, limiter=None, scaler=None, row_names=None,
                         linsolver=None, **kw):
            self.n += 1
            if self.n == 1:
                raise NoConvergenceError('forced first failure')
            return self.inner.solve_system(
                x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter,
                limiter, scaler, row_names=row_names, linsolver=linsolver)

    def eval_tri(x):
        return (np.array([x[0]**3 - x[0]]),
                np.array([[3.0*x[0]**2 - 1.0]]))

    for seed, basin_root in ((0.9, 1.0), (-0.9, -1.0), (0.05, 0.0)):
        x_pt, _ = PseudoTransientNewton(FailFirst(MiniNewton())).solve_system(
            np.array([seed]), eval_tri, *tols)
        assert x_pt[0] == pytest.approx(basin_root, abs=1e-9), seed

    ## (4): a dead chain as base, a live plain solver for the rungs --
    ## the dcanalysis wiring.  The chain fails once (the pure attempt);
    ## every deformed solve goes to the rung solver.
    class DeadChain(NonLinearSolver):
        def solve_system(self, *a, **kw):
            raise NoConvergenceError('whole chain failed')

    x, _ = PseudoTransientNewton(DeadChain(),
                                 rung_solver=MiniNewton()).solve_system(
        np.zeros(1), eval_cubic, *tols)
    assert x[0] == pytest.approx(root, abs=1e-12)

    ## Passthrough: a base that succeeds is returned verbatim, no rungs.
    class Instant(NonLinearSolver):
        def solve_system(self, x0, *a, **kw):
            return x0 + 1.0, 1
    x, _ = PseudoTransientNewton(Instant()).solve_system(
        np.zeros(1), eval_cubic, *tols)
    assert x[0] == 1.0
