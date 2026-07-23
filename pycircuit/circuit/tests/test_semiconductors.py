import pytest
import numpy as np
from pycircuit.circuit.circuit import SubCircuit, Node, gnd
from pycircuit.circuit.elements import VS, IS, R
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.analysis_ss import AC
from pycircuit.circuit.semiconductors import BJT, JFET, ZenerDiode, Varactor

def test_bjt_forward_active():
    c = SubCircuit()
    c['VBE'] = VS('b', 'e', v=0.7)
    c['VCE'] = VS('c', 'e', v=5.0)
    
    # NPN BJT
    # VT = kT/q = 0.02585 V at 300K
    # i_c approx IS * exp(v_be / VT)
    c['Q1'] = BJT('c', 'b', 'e', IS=1e-14, BF=100.0)
    # Emitter tied to gnd for reference
    c['RE'] = R('e', gnd, r=1e-6)
    
    dc = DC(c)
    res = dc.solve()
    
    VT = dc.toolkit.kboltzmann * dc.epar.T / dc.toolkit.qelectron
    v_be = 0.7
    expected_ic = 1e-14 * (np.exp(v_be / VT) - 1.0)
    # With early effect VA=100 and VCE=5.0
    expected_ic *= (1.0 + 5.0 / 100.0)
    
    i_c = -res.i('VCE.plus') # current from VCE into C
    assert abs(i_c - expected_ic) / expected_ic < 0.01

def test_jfet_saturation():
    c = SubCircuit()
    # VTO = -2.0V
    # VGS = -1.0V (VGS > VTO => Not cutoff)
    # VDS = 5.0V
    # VDS > VGS - VTO => 5.0 > -1.0 - (-2.0) = 1.0 => Saturation
    c['VGS'] = VS('g', 's', v=-1.0)
    c['VDS'] = VS('d', 's', v=5.0)
    
    c['J1'] = JFET('d', 'g', 's', VTO=-2.0, BETA=1e-3, LAMBDA=0.01)
    c['RS'] = R('s', gnd, r=1e-6)
    
    dc = DC(c)
    res = dc.solve()
    
    expected_id = 1e-3 * (-1.0 - (-2.0))**2 * (1.0 + 0.01 * 5.0)
    
    i_d = -res.i('VDS.plus')
    assert abs(i_d - expected_id) / expected_id < 0.01

def test_zener_breakdown():
    c = SubCircuit()
    # Apply -6V to a 5V Zener
    c['V1'] = VS('cathode', 'anode', v=6.0)
    c['Z1'] = ZenerDiode('anode', 'cathode', BV=5.0, IBV=1e-10)
    c['R1'] = R('anode', gnd, r=1e-6)
    
    dc = DC(c)
    res = dc.solve()
    
    VT = dc.toolkit.kboltzmann * dc.epar.T / dc.toolkit.qelectron
    v_zener = 6.0 - 5.0 # 1.0V over breakdown
    expected_iz = 1e-10 * (np.exp(1.0 / VT) - 1.0)
    
    i_z = -res.i('V1.plus') # V1 cathode is 'plus'
    assert abs(i_z - expected_iz) / expected_iz < 0.01

def test_varactor():
    c = SubCircuit()
    # Varactor with minus at ground
    c['VAR'] = Varactor('plus', gnd, CJ0=1e-12, VJ=1.0, M=0.5)
    
    # Apply DC bias of -0.5V and AC of 1.0V to 'plus' using a VS
    # Wait, we want to measure current. If we use a VS on 'plus', the current
    # comes out of the VS. The VS branch current is what we want!
    c['V1'] = VS('plus', gnd, v=-0.5, vac=1.0)
    
    ac = AC(c)
    w = 1.0
    res = ac.solve(freqs=np.array([w / (2 * np.pi)]))
    
    expected_c = 1e-12 / np.sqrt(1.0 - (-0.5)/1.0)
    
    # AC V_plus = 1.0V.
    # Current leaving Varactor at 'plus' = I_var = j * w * C * V_plus
    # The current entering V1 is what we measure.
    # I_V1 + I_var = 0 => I_V1 = - I_var
    # We can measure current through V1 directly!
    I_v1 = res.i('V1.plus')[0] # Current from V1 plus terminal
    
    # V1 plus terminal current is the current entering the element from the plus node
    # I_v1_plus + I_var_plus = 0 => I_v1_plus = -I_var_plus = -j * w * C * 1.0
    expected_i = -1j * w * expected_c
    
    assert abs(I_v1 - expected_i) < 1e-15
