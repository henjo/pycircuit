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
    def __init__(self):
        self.call_count = 0
        
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None):
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
        def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None):
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
