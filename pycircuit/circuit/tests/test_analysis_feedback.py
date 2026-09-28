import pytest
# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

""" Test loopgain module
"""

from pycircuit.circuit import analysis, analysis_ss
from pycircuit.circuit import symbolic, SubCircuit, R, C, VS, VCCS, VCVS, gnd
from pycircuit.circuit.feedback import FeedbackDeviceAnalysis, LoopProbe, FeedbackLoopAnalysis
import sympy
from sympy import simplify
import numpy as np
from numpy.testing import *

def test_deviceanalysis_sourcefollower():
    """Loopgain of a source follower"""

    gm,RL,CL,s = sympy.symbols('gm RL CL s')

    cir = SubCircuit(toolkit=symbolic)
    cir['M1'] = VCCS('g', 's', gnd, 's', gm = gm,toolkit=symbolic)
    cir['RL'] = R('s', gnd, r=RL)
    cir['CL'] = C('s', gnd, c=CL)
    cir['VS'] = VS('g', gnd)

    ana = FeedbackDeviceAnalysis(cir, 'M1', toolkit=symbolic)
    res = ana.solve(s, complexfreq=True)

    assert simplify(res['loopgain']) == simplify(- gm / (1/RL + s*CL))

def test_deviceanalysis_viiv():
    """Loopgain of a resistor V-I and a I-V amplifier with a vcvs as gain element"""

    sympy.var('R1 R2 CL A s')

    cir = SubCircuit(toolkit=symbolic)
    cir['A1'] = VCVS(gnd, 'int', 'out', gnd, g = A,toolkit=symbolic)
    cir['R1'] = R('in', 'int', r=R1)
    cir['R2'] = R('int', 'out', r=R2)
    cir['VS'] = VS('in', gnd)

    ana = FeedbackDeviceAnalysis(cir, 'A1', toolkit=symbolic)
    res = ana.solve(s, complexfreq=True)

    assert simplify(res['loopgain'] - (- A * R1 / (R1 + R2))) == 0

def test_loopanalysis_incorrect_circuit():
    cir = SubCircuit()
    with pytest.raises(ValueError):
        FeedbackLoopAnalysis(cir)

    cir['probe1'] = LoopProbe('out', gnd, 'out_R2', gnd)
    cir['probe2'] = LoopProbe('out', gnd, 'out_R2', gnd)
    with pytest.raises(ValueError):
        FeedbackLoopAnalysis(cir)

def test_loopanalysis_numeric():
    cir = SubCircuit()
    cir['A1'] = VCVS(gnd, 'int', 'out', gnd, g = 10)
    cir['R1'] = R('in', 'int')
    cir['probe'] = LoopProbe('out', gnd, 'out_R2', gnd)
    cir['C2'] = C('int', 'out_R2')
    cir['VS'] = VS('in', gnd)

    ana = FeedbackLoopAnalysis(cir)
    res = ana.solve(np.logspace(4,6), complexfreq=True)
    print(abs(res['loopgain']))
    
def test_loopanalysis_viiv():
    sympy.var('R1 R2 CL A s')

    cir = SubCircuit(toolkit=symbolic)
    cir['A1'] = VCVS(gnd, 'int', 'out', gnd, g = A, toolkit=symbolic)
    cir['R1'] = R('in', 'int', r=R1)
    cir['probe'] = LoopProbe('out', gnd, 'out_R2', gnd, toolkit=symbolic)
    cir['R2'] = R('int', 'out_R2', r=R2)
    cir['VS'] = VS('in', gnd)

    ana = FeedbackLoopAnalysis(cir, toolkit=symbolic)

    res = ana.solve(s, complexfreq=True)

    assert sympy.simplify(res['loopgain'] - (- A * R1 / (R1 + R2))) == 0


def test_the_device_loop_gain_does_not_depend_on_an_earlier_analysis():
    """`FeedbackDeviceAnalysis` evaluates at its own fixed point, and reads a
    stateful limiter's device there (`limit_sync`).  ⚠ A `Diode` reads `G`
    as the tangent at its stored `_vlim`, so the loop gain depended on
    where an earlier analysis had left it: after a DC solve at a forward
    bias it read the diode's forward conductance (2026-09-28)."""
    from pycircuit.circuit import numeric
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.elements import Diode

    def build():
        cir = SubCircuit()
        cir['M1'] = VCCS('g', 's', gnd, 's', gm=20e-3)
        cir['RL'] = R('s', gnd, r=1e3)
        cir['D'] = Diode('s', gnd)
        cir['VS'] = VS('g', gnd, v=5.0)
        return cir

    def loopgain(cir):
        return complex(np.asarray(
            FeedbackDeviceAnalysis(cir, 'M1', toolkit=numeric).solve(1e3)['loopgain']).ravel()[0])
    cir = build()
    DC(cir, toolkit=numeric).solve()
    lg, lg_ref = loopgain(cir), loopgain(build())
    assert abs(lg / lg_ref - 1.0) < 1e-12, (lg, lg_ref)
