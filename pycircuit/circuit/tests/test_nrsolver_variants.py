import pytest
import numpy as np
from pycircuit.circuit.nrsolver import NonLinearSolver, StandardNewton, GminSteppingNewton, SourceSteppingNewton, NoConvergenceError
from pycircuit.circuit import numeric

class MockFailingSolver(NonLinearSolver):
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None):
        try:
            eval_FJ(x0)
        except Exception:
            pass
        raise NoConvergenceError("Mock failure")

class MockSourceSteppingBaseSolver(NonLinearSolver):
    def __init__(self):
        self.call_count = 0
        
    def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None):
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
        def solve_system(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None):
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
