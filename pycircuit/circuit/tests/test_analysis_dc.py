import numpy as np
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import Diode, IS
from pycircuit.circuit.dcanalysis import DC

def test_dc_false_convergence_kcl():
    """
    Tests that the KCL Verification properly prevents the False Convergence Phenomenon
    in highly non-linear, high-conductance scenarios.
    We set the temperature to 50K to increase the diode steepness (J/I ratio),
    guaranteeing that Delta-V can drop below xtol (1e-6) while the KCL error is still
    larger than reltol (1e-4).
    """
    from pycircuit.circuit.circuit import defaultepar, ParameterDict, Parameter
    from pycircuit.circuit.elements import Diode, IS
    from pycircuit.circuit.dcanalysis import DC

    cold_epar = ParameterDict(Parameter("T", "Temperature", unit="K", default=50))
    
    c = SubCircuit()
    c['is'] = IS(gnd, 1, i=1e6) 
    c['D1'] = Diode(1, gnd)

    dc = DC(c, epar=cold_epar)
    dc.par.maxiter = 1000
    res = dc.solve()
    
    # Verify KCL explicitly using PyCircuit's exact mathematical residual F
    F = dc.cir.i(res.x, cold_epar) + dc.cir.u(0.0, cold_epar, analysis='dc')
    
    # In a false convergence scenario without KCL verification, the solver stops early
    # because the Delta-V step is tiny, leaving a massive current error F.
    # With the new KCL verification, the solver strictly enforces that F is within
    # reltol * I_scale + abstol.
    
    # We expect the KCL error at node 1 to be completely resolved.
    F_node1 = abs(F[1])
    I_scale = 1e6 # The absolute scale of currents at node 1
    
    assert F_node1 <= 1e-4 * I_scale + 1e-12, f"False Convergence detected! KCL error: {F_node1}A"

def test_dc_continuation_methods():
    """
    Test the robust execution of GminSteppingNewton and SourceSteppingNewton.
    We create a highly non-linear circuit with a poor initial guess that normally
    fails StandardNewton, ensuring the decorators catch the failure and step it.
    """
    from pycircuit.circuit.circuit import ParameterDict, Parameter
    from pycircuit.circuit.elements import Diode, IS, R
    from pycircuit.circuit.dcanalysis import DC

    c = SubCircuit()
    c['IS1'] = IS(gnd, 1, i=1.0) # High current into a diode
    c['D1'] = Diode(1, gnd)
    
    # A 1A current into a standard generic diode is extremely high and usually fails
    # basic Newton-Raphson from an x0=0 guess due to exponential overflow or
    # severe oscillation. Gmin/Source stepping makes it converge safely.
    dc = DC(c)
    
    try:
        res = dc.solve()
        assert res.v(1) > 0.5, "Continuation methods failed to find the forward biased diode point"
    except Exception as e:
        assert False, f"Continuation method failed to catch non-convergence! {str(e)}"

