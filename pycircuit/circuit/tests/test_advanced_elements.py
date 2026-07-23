import pytest
from numpy.testing import assert_almost_equal
import sympy
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import R, VPWL, VExp, CoupledInductors, VSwitch, VS
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.toolkit import numeric, symbolic, symbolic_poly

TOOLKITS = [numeric, symbolic]

@pytest.mark.parametrize('toolkit', [numeric])
def test_vpwl(toolkit):
    """
    Test Piecewise Linear Voltage Source (VPWL).
    
    Testbench Description:
    - V1 is a VPWL source connected to node '1'.
    - R1 is a 1-ohm load connected to node '1'.
    - We specify time-voltage pairs `tvpairs` and verify that the `function.f(t)` 
      method interpolates the correct voltages at t=0.5s and t=1.5s.
    - We also run a DC analysis which defaults to t=0, and verify node '1' is 0V.
    """
    c = SubCircuit(toolkit=toolkit)
    tvpairs = [0.0, 0.0, 1.0, 2.0, 2.0, 2.0]
    c['V1'] = VPWL('1', gnd, tvpairs=tvpairs)
    c['R1'] = R('1', gnd, r=1.0)
    
    assert c['V1'].function.f(0.5) == 1.0
    assert c['V1'].function.f(1.5) == 2.0
    
    dc = DC(c, toolkit=toolkit)
    res = dc.solve()
    assert abs(res.v('1') - 0.0) < 1e-6

@pytest.mark.parametrize('toolkit', [numeric])
def test_vexp(toolkit):
    """
    Test Exponential Voltage Source (VExp).
    
    Testbench Description:
    - V1 is a VExp source connected to node '1'.
    - R1 is a dummy 1-ohm load connected to node '1'.
    - We verify that `function.f(t)` returns the exact exponentially growing 
      voltage at t=2.0s based on the given time constants.
    """
    c = SubCircuit(toolkit=toolkit)
    c['V1'] = VExp('1', gnd, v1=0.0, v2=5.0, td1=1.0, tau1=1.0, td2=3.0, tau2=1.0)
    c['R1'] = R('1', gnd, r=1.0)
    
    v_at_2 = 0.0 + 5.0 * (1 - numeric.exp(-1.0))
    assert abs(c['V1'].function.f(2.0) - v_at_2) < 1e-6

@pytest.mark.parametrize('toolkit', TOOLKITS)
def test_coupled_inductors(toolkit):
    """
    Test Mutual Inductance (CoupledInductors).
    
    Testbench Description:
    - K1 couples two inductors between node '1'/gnd and node '2'/gnd.
    - L1 is 1uH, L2 is 4uH, and coupling factor K is 0.5.
    - We extract the Matrix C (capacitive/inductive term) from the element 
      and verify the off-diagonal coupling terms equal exactly -M 
      where M = K * sqrt(L1 * L2).
    - Tests numeric and symbolic mathematical representations.
    """
    c = SubCircuit(toolkit=toolkit)
    # Use symbolic symbols if toolkit is symbolic
    if toolkit is numeric:
        L1_val, L2_val, K_val = 1e-6, 4e-6, 0.5
    else:
        L1_val, L2_val, K_val = sympy.symbols('L1 L2 K', positive=True)
        
    c['K1'] = CoupledInductors('1', gnd, '2', gnd, L1=L1_val, L2=L2_val, K=K_val, toolkit=toolkit)
    c['R1'] = R('1', gnd, r=1.0)
    c['R2'] = R('2', gnd, r=1.0)
    
    C_mat = c['K1'].C(None)
    if toolkit is numeric:
        M = 1e-6
        assert abs(C_mat[4, 5] - (-M)) < 1e-15
        assert abs(C_mat[5, 4] - (-M)) < 1e-15
    else:
        M = K_val * sympy.sqrt(L1_val * L2_val)
        assert sympy.simplify(C_mat[4, 5] - (-M)) == 0
        assert sympy.simplify(C_mat[5, 4] - (-M)) == 0

def test_vswitch_numeric():
    """
    Test Voltage Controlled Switch (VSwitch) standalone evaluation.
    
    Testbench Description:
    - Directly tests the `i(x)` evaluation hook.
    - Evaluates the smoothed tanh mathematical step transition.
    """
    c = SubCircuit(toolkit=numeric)
    c['Vctrl'] = VS('ctrl', gnd, v=1.0) # Control voltage is 1.0V (Von)
    c['VSW'] = VSwitch('1', gnd, 'ctrl', gnd, Ron=1.0, Roff=1e6, Von=1.0, Voff=0.0)
    c['Iin'] = R('1', gnd, r=1.0) # Dummy
    
    i_func = c['VSW'].i
    x = numeric.array([1.0, 0.0, 1.0, 0.0]) # v=1, vc=1
    i_out = i_func(x)
    Gon = 1.0; Goff = 1e-6
    factor = (numeric.tanh(1.0) + 1.0) / 2.0
    g_expected = Goff + (Gon - Goff) * factor
    assert abs(i_out[0] - g_expected) < 1e-6

def test_vswitch_dc():
    """
    Test Voltage Controlled Switch (VSwitch) within a DC non-linear solver.
    
    Testbench Description:
    - Iin injects 1.0A into node '1'.
    - VSW is placed as the only path to ground for this current, acting as the load.
    - Vctrl applies 1.0V to the VSW control nodes, turning it ON.
    - With the switch ON, its resistance is ~1.135 ohms (based on the tanh transition).
    - Ohm's law: 1.0A flowing through ~1.135 ohms means node '1' must be ~1.135V.
    - The non-linear Newton-Raphson solver validates this analytically.
    """
    from pycircuit.circuit.elements import IS
    c = SubCircuit(toolkit=numeric)
    c['Iin'] = IS(gnd, '1', i=1.0) # Inject 1A
    c['Vctrl'] = VS('ctrl', gnd, v=1.0) # Vctrl = 1V => ON
    c['VSW'] = VSwitch('1', gnd, 'ctrl', gnd, Ron=1.0, Roff=1e6, Von=1.0, Voff=0.0)
    
    dc = DC(c)
    res = dc.solve()
    # Resistance is ~ 1.135 ohms, so 1A injected gives ~ 1.135V
    assert 1.13 < res.v('1') < 1.14

def test_cccs():
    """
    Test Current Controlled Current Source (CCCS).
    
    Testbench Description:
    - Iin injects 1.0A into node '1'. R1 provides a path to ground.
    - CCCS input branch measures the current from node '1' to '2'.
    - CCCS output branch draws F * I(1->2) from node '3'.
    - R3 provides a 1-ohm path to ground from node '3'.
    """
    from pycircuit.circuit.elements import IS, CCCS
    c = SubCircuit(toolkit=numeric)
    c['Iin'] = IS(gnd, '1', i=1.0)
    c['R1'] = R('1', '2', r=1.0)
    c['F1'] = CCCS('2', gnd, '3', gnd, F=5.0) # Gain of 5
    c['R3'] = R('3', gnd, r=1.0)
    
    dc = DC(c)
    res = dc.solve()
    # 1A flows through R1 and into F1 input branch
    # So F1 output branch draws 5A from node 3.
    # V(3) = -5A * 1ohm = -5V
    assert abs(res.v('3') - (-5.0)) < 1e-6

def test_iswitch():
    """
    Test Current Controlled Switch (ISwitch).
    
    Testbench Description:
    - Ictrl injects a control current of 2mA into 'cp'.
    - Iin injects 1.0A into node '1'.
    - ISwitch turns ON when control current exceeds 1mA, connecting '1' and '2'.
    - Computes resulting node voltages.
    """
    from pycircuit.circuit.elements import IS, ISwitch
    c = SubCircuit(toolkit=numeric)
    c['Ictrl'] = IS(gnd, 'cp', i=2e-3)
    c['SW'] = ISwitch('1', '2', 'cp', gnd, Ron=1.0, Roff=1e6, Ion=1e-3, Ioff=0.0)
    c['Iin'] = IS(gnd, '1', i=1.0)
    c['Rload'] = R('2', gnd, r=1.0)
    
    dc = DC(c)
    res = dc.solve()
    # Control current is 2mA. Thresholds are 1mA (on) and 0 (off).
    # Since 2mA > 1mA, the switch is fully ON.
    # Ron = 1.0 ohm. Rload = 1.0 ohm. Total R = 2.0 ohms.
    # So V(2) should be 1.0V and V(1) should be 2.0V.
    assert abs(res.v('2') - 1.0) < 0.05
    assert abs(res.v('1') - 2.0) < 0.05

def test_nonlinear_vccs():
    """
    Test Behavioral B-Source (NonLinearVCCS).
    
    Testbench Description:
    - V1 applies a voltage to the control nodes.
    - B1 evaluates i_out = v_ctrl^2.
    - R1 is a load resistor on the output.
    """
    from pycircuit.circuit.elements import NonLinearVCCS
    c = SubCircuit(toolkit=numeric)
    c['V1'] = VS('ctrl', gnd, v=3.0)
    
    # Quadratic VCCS: Iout = Vctrl^2
    c['B1'] = NonLinearVCCS('ctrl', gnd, 'out', gnd, i_func=lambda v: v**2)
    c['R1'] = R('out', gnd, r=2.0)
    
    dc = DC(c)
    res = dc.solve()
    # Vctrl = 3.0V. Iout = 9.0A drawn from 'out'.
    # Vout = -Iout * R1 = -9.0 * 2.0 = -18.0V.
    assert abs(res.v('out') - (-18.0)) < 1e-6

def test_bsource_capacitance():
    """
    Test BSource non-linear capacitance feature (Q-Source).
    
    Testbench Description:
    - V1 applies 2.0V AC to the control nodes.
    - BSource evaluates q_out = v_ctrl^3.
    - Capacitance should be dQ/dV = 3 * v_ctrl^2.
    - At Vctrl = 2.0V, Capacitance = 12.0 F.
    - An AC analysis is performed to verify the admittance Y = j*w*C.
    """
    from pycircuit.circuit.elements import BSource, IS, VS
    from pycircuit.circuit.analysis_ss import AC
    c = SubCircuit()
    c['V1'] = VS('ctrl', gnd, v=2.0, vac=1.0)
    
    # Non-linear transcapacitor: Q_out = Vctrl^3
    c['B1'] = BSource('ctrl', gnd, 'out', gnd, q_func=lambda v: v**3)
    
    # We load 'out' with a 1 ohm resistor to measure the current.
    # No AC current source needed on 'out'.
    c['R1'] = R('out', gnd, r=1.0)
    
    # AC analysis at 1 rad/s => w = 1.
    ac = AC(c)
    import numpy as np
    res = ac.solve(freqs=np.array([1.0/(2*3.141592653589793)])) # f = 1/(2*pi) Hz => w = 1 rad/s
    
    # At DC bias of 2.0V, C_trans = 3*(2.0)^2 = 12.0 F.
    # Vctrl has AC magnitude 1.0.
    # The BSource draws I_out = j * w * C_trans * Vctrl_ac = j * 1 * 12.0 * 1.0 = j12 A.
    # This current comes from 'out', which is connected to ground via R1 (1 ohm).
    # Therefore, Vout = -I_out * 1.0 = -j12 V.
    Vout = res.v('out')[0]
    expected_Vout = -12j
    
    # Verify magnitude and phase
    assert abs(Vout - expected_Vout) < 1e-6
