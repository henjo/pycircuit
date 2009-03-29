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

    gm,RL = sympy.symbols('gm RL')

    cir = SubCircuit()
    cir['M1'] = VCCS('g', 's', gnd, 's', gm = gm)
    cir['RL'] = R('s', gnd, r=RL)
    cir['VS'] = VS('g', gnd)

    res = MyFeedbackDeviceAnalysis(cir, 'M1', 'inp', 'inn', 'outp', 'outn').solve([1e3])

    assert_equal(simplify(res['loopgain']), - gm * RL)


    

