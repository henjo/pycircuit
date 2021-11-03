# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from nose.tools import *
import pycircuit.circuit.circuit 
from pycircuit.circuit import *
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from .test_circuit import create_current_divider
import unittest

def test_integer_component_values():
    """Test dc analysis with integer component values
    
       As python per default uses integer arithmetics integer
       component values can lead to problems. (Python < 3.0)
    """
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)

    c['vs'] = VS('net1', gnd, v = 9)
    c['R1'] = R( 'net1',  'net2', r = 50)
    c['R2'] = R( 'net2', gnd, r = 50)

    dc = DC(c)
    res = dc.solve()
    
    assert_equal(res.v('net2'), 4.5)

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

    n1,n2 = c.add_nodes('net1', 'net2')

    c['vs'] = VS(n1, gnd, v = 9.)
    c['R1'] = R( n1,  n2, r = 50.)
    c['R2'] = R( n2, gnd, r = 50.)
    
    noise = Noise(c, inputsrc='vs', outputnodes=(n2, gnd))
    should = np.array([noise.solve(0)['Svnout'],noise.solve(1)['Svnout']])
    res = noise.solve(np.array([0,1]))
    assert_array_equal(res['Svnout'], should)

