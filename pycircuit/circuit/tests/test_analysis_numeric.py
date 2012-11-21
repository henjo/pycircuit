# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from nose.tools import *
import pycircuit.circuit.circuit 
from pycircuit.circuit import *
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from test_circuit import create_current_divider
import unittest
from .. import theanotk

def test_integer_component_values():
    """Test dc analysis with integer component values
    
       As python per default uses integer arithmetics integer
       component values can lead to problems. (Python < 3.0)
    """
    c = SubCircuit()

    c['vs'] = VS('net1', gnd, v = 9)
    c['R1'] = R( 'net1', 'net2', r = 50)
    c['R2'] = R( 'net2', gnd, r = 50)

    c = c.save_current('R1.plus')
    c = c.save_current('R2.plus')

    dc = DC(c, toolkit=theanotk)
    res = dc.solve()
    
    assert_equal(res.v('net2'), 4.5)

    assert_equal(res.i('R1.plus'), 0.09)
    assert_equal(res.i('R2.plus'), 0.09)

def TODOtest_noise_dc_steady_state():
    """Test that dc-steady state is accounted for in noise simulations
    """
    pass

@unittest.skip("Skip failing test")
def test_noise_with_frequency_vector():
    """Test that noise analysis support an array as input argument for frequency

    """
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)

    c['vs'] = VS(1, gnd, v = 9.)
    c['R1'] = R(1,  2,   r = 50.)
    c['R2'] = R(2,  gnd, r = 50.)
    
    noise = Noise(c, inputsrc='vs', outputnodes=(n2, gnd))
    should = np.array([noise.solve(0)['Svnout'],noise.solve(1)['Svnout']])
    res = noise.solve(np.array([0,1]))
    assert_array_equal(res['Svnout'], should)

