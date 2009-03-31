# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

""" Test loopgain module
"""

from nose.tools import *
from pycircuit.circuit import analysis
from pycircuit.circuit import *
from pycircuit.circuit.feedback import FeedbackDeviceAnalysis, LoopProbe, FeedbackLoopAnalysis
import sympy

class MyFeedbackDeviceAnalysis(FeedbackDeviceAnalysis, SymbolicAnalysis):
    pass
class MyFeedbackLoopAnalysis(FeedbackLoopAnalysis, SymbolicAnalysis):
    pass

def test_deviceanalysis_sourcefollower():
    """Loopgain of a source follower"""

    gm,RL,CL,s = sympy.symbols('gm RL CL s')

    cir = SubCircuit()
    cir['M1'] = VCCS('g', 's', gnd, 's', gm = gm)
    cir['RL'] = R('s', gnd, r=RL)
    cir['CL'] = C('s', gnd, c=CL)
    cir['VS'] = VS('g', gnd)

    ana = MyFeedbackDeviceAnalysis(cir, 'M1')
    res = ana.solve(s, complexfreq=True)

    assert_equal(simplify(res['loopgain']), simplify(- gm / (1/RL + s*CL)))

def test_deviceanalysis_viiv():
    """Loopgain of a resistor V-I and a I-V amplifier with a vcvs as gain element"""

    sympy.var('R1 R2 CL A s')

    cir = SubCircuit()
    cir['A1'] = VCVS(gnd, 'int', 'out', gnd, g = A)
    cir['R1'] = R('in', 'int', r=R1)
    cir['R2'] = R('int', 'out', r=R2)
    cir['VS'] = VS('in', gnd)

    ana = MyFeedbackDeviceAnalysis(cir, 'A1')
    res = ana.solve(s, complexfreq=True)

    assert simplify(res['loopgain'] - (- A * R1 / (R1 + R2))) == 0

def test_loopanalysis_incorrect_circuit():
    cir = SubCircuit()
    assert_raises(ValueError, lambda: MyFeedbackLoopAnalysis(cir))

    cir['probe1'] = LoopProbe('out', gnd, 'out_R2', gnd)
    cir['probe2'] = LoopProbe('out', gnd, 'out_R2', gnd)
    assert_raises(ValueError, lambda: MyFeedbackLoopAnalysis(cir))

def test_loopanalysis_viiv():
    sympy.var('R1 R2 CL A s')

    cir = SubCircuit()
    cir['A1'] = VCVS(gnd, 'int', 'out', gnd, g = A)
    cir['R1'] = R('in', 'int', r=R1)
    cir['probe'] = LoopProbe('out', gnd, 'out_R2', gnd)
    cir['R2'] = R('int', 'out_R2', r=R2)
    cir['VS'] = VS('in', gnd)

    ana = MyFeedbackLoopAnalysis(cir)
    res = ana.solve(s, complexfreq=True)

    assert simplify(res['loopgain'] - (- A * R1 / (R1 + R2))) == 0
    

