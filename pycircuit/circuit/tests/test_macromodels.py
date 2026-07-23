import pytest
import numpy as np
from pycircuit.circuit.circuit import SubCircuit, Node, gnd
from pycircuit.circuit.elements import VS, R
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.analysis_ss import AC
from pycircuit.circuit.macromodels import OpAmp

def test_opamp_dc_gain():
    c = SubCircuit()
    c['V1'] = VS('inp', gnd, v=1e-5) # Very small input to stay in linear region
    c['OP1'] = OpAmp('inp', gnd, 'outp', gnd, Aol=1e5, Vdd=15.0, Vss=-15.0)
    
    dc = DC(c)
    res = dc.solve()
    
    v_out = res.v('outp')
    expected_vout = 1e-5 * 1e5 # 1.0 V
    
    # We allow a small error due to numerical precision and tanh compression
    assert abs(v_out - expected_vout) / expected_vout < 2e-3

def test_opamp_saturation():
    c = SubCircuit()
    c['V1'] = VS('inp', gnd, v=1.0) # 1.0V input with 1e5 gain will saturate
    c['OP1'] = OpAmp('inp', gnd, 'outp', gnd, Aol=1e5, Vdd=15.0, Vss=-15.0)
    
    dc = DC(c)
    res = dc.solve()
    
    v_out = res.v('outp')
    
    # Should saturate near 15.0V
    assert abs(v_out - 15.0) < 0.1

def test_opamp_ac_bandwidth():
    c = SubCircuit()
    # Inverting amplifier with gain of -10
    c['V1'] = VS('vin', gnd, v=0.0, vac=1.0)
    c['R1'] = R('vin', 'inn', r=10e3)
    c['R2'] = R('inn', 'outp', r=100e3)
    
    # Non-inverting input grounded
    c['OP1'] = OpAmp(gnd, 'inn', 'outp', gnd, Aol=1e5, GBW=1e6)
    
    ac = AC(c)
    # The closed loop gain is -10 (20 dB).
    # The closed loop bandwidth is GBW / (1 + |Gain|) = 1e6 / 11 ~= 90.9 kHz
    # Let's test at DC (1 Hz) and near bandwidth
    
    f1 = 1.0
    f2 = 1e6 / 11.0 # -3dB point
    
    res = ac.solve(freqs=np.array([f1, f2]))
    
    vout_dc = res.v('outp')[0]
    vout_3db = res.v('outp')[1]
    
    # DC gain should be approximately -10
    assert abs(abs(vout_dc) - 10.0) < 0.1
    
    # At 3dB point, magnitude should be 10 / sqrt(2) ~= 7.07
    expected_3db = 10.0 / np.sqrt(2)
    assert abs(abs(vout_3db) - expected_3db) / expected_3db < 0.05
