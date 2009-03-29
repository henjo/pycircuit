# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

""" Test loopgain module
"""

from nose.tools import *
from pycircuit.circuit import analysis
from pycircuit.circuit import *
from pycircuit.circuit.feedback import FeedbackDeviceAnalysis
import sympy

class MyFeedbackDeviceAnalysis(FeedbackDeviceAnalysis, SymbolicAnalysis):
    pass

def test_simple():
    """Loopgain of a source follower"""

    gm,RL,CL,s = sympy.symbols('gm RL CL s')

    cir = SubCircuit()
    cir['M1'] = VCCS('g', 's', gnd, 's', gm = gm)
    cir['RL'] = R('s', gnd, r=RL)
    cir['CL'] = C('s', gnd, c=CL)
    cir['VS'] = VS('g', gnd)

    ana = MyFeedbackDeviceAnalysis(cir, 'M1', 'inp', 'inn', 'outp', 'outn')
    res = ana.solve(s, complexfreq=True)

    assert_equal(simplify(res['loopgain']), simplify(- gm / (1/RL + s*CL)))


    

