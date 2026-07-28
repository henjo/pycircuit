import pytest
import numpy as np

from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import Diode, IS, R, VS
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import EulerIntegrator

def test_bypass_dc_consistency():
    """Verify that DC analysis with bypass=True converges to the same result as bypass=False."""
    c = SubCircuit()
    c['is'] = IS(gnd, 1, i=1.0)
    c['D1'] = Diode(1, gnd)

    dc_off = DC(c, bypass=False)
    res_off = dc_off.solve()
    v_off = res_off.v(1)

    dc_on = DC(c, bypass=True)
    res_on = dc_on.solve()
    v_on = res_on.v(1)

    # They should be extremely close
    np.testing.assert_allclose(v_on, v_off, rtol=1e-5)

def test_bypass_transient_consistency():
    """Verify that Transient analysis with bypass=True yields consistent waveforms as bypass=False."""
    c = SubCircuit()
    c['vs'] = VS(1, gnd, v=5.0)
    c['r1'] = R(1, 2, r=100)
    c['d1'] = Diode(2, gnd)

    tran_off = Transient(c, bypass=False, integrator=EulerIntegrator())
    res_off = tran_off.solve(timestep=1e-6, tend=10e-6)
    v_off = res_off.v(2)

    tran_on = Transient(c, bypass=True, integrator=EulerIntegrator())
    res_on = tran_on.solve(timestep=1e-6, tend=10e-6)
    v_on = res_on.v(2)

    # Verify both simulated waveforms are essentially identical
    np.testing.assert_allclose(v_on, v_off, rtol=1e-4)
