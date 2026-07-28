import pytest
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import R, C, Diode, VPulse
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import EulerIntegrator
from pycircuit.circuit.nrsolver import NoConvergenceError

def test_transient_minstep_convergence_failure():
    """Verify that Transient simulator correctly halts and raises RuntimeError if dt shrinks below minstep."""
    c = SubCircuit()
    c['VP'] = VPulse('in', gnd, v1=0, v2=1000, tr=1e-12, tf=1e-12, pw=1e-6, per=2e-6)
    c['R'] = R('in', 'out', r=1e-3)
    c['D'] = Diode('out', gnd) # High voltage into diode will cause convergence issues if maxiter is artificially constrained
    
    # We constrain maxiter to 2 so the Newton solver fails to converge on the massive voltage step
    # This will force the transient solver to cut the timestep.
    # We set minstep to 1e-10. Since it starts with dt=1e-9, it will cut 1e-9 -> 2.5e-10 -> 6.25e-11 < 1e-10
    # Thus, it should raise RuntimeError on the 2nd cut!
    tran = Transient(c, integrator=EulerIntegrator(), minstep=1e-10, maxiter=2, bypass=False)
    
    with pytest.raises(RuntimeError, match="timestep shrank below"):
        res = tran.solve(tend=2e-6, timestep=1e-9)
