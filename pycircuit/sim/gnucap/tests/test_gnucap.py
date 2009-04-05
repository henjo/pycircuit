# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.sim
import pycircuit.sim.gnucap as gnucap
from pycircuit.post import Waveform
import numpy as np

import logging
from nose.tools import *

#logging.basicConfig(level = logging.DEBUG)

simple_netlist = """V1 1 0 1.0
R1 1 2 1k
R2 2 0 3k
"""

def test_start():
    sim = gnucap.Simulation(None)
    
    print sim.send('list')

def test_netlist():
    sim = gnucap.Simulation(None)

    sim.send_netlist(simple_netlist)

    netlist = sim.send('list')
    
    assert_equal(netlist, '\r\nV1 ( 1 0 )  DC  1.\r\nR1 ( 1 2 )  1.K\r\nR2 ( 2 0 )  3.K')

def test_raw_dcop():
    sim = gnucap.Simulation(None)

    sim.send_netlist(simple_netlist)
    
    sim.send('print op v(nodes)')

    res = sim.send_command('op 27')

    assert_equal(res['v(2)'].value(27), 0.75)

def test_dcop():
    cir = gnucap.Circuit(simple_netlist)
    sim = gnucap.Simulation(cir)

    res = sim.run_dcop()

    assert_equal(res['v(2)'].value(27), 0.75)

def test_ac():
    cir = gnucap.Circuit()
    cir['V1'] = gnucap.VS(1, 0, 1, ac=1)
    cir['R'] = gnucap.R(1,2,10e3)
    cir['C'] = gnucap.C(2,0,1e-12)

    sim = gnucap.Simulation(cir)
    
    res = sim.run_ac(1e3,'1G',decade=10)
    

def test_transient():
    cir = gnucap.Circuit(simple_netlist)
    sim = gnucap.Simulation(cir)

    res = sim.run_transient(0,1,0.2)

    refv2 = Waveform(np.array([0.,0.2,0.4,0.6,0.8,1.]), 
                     np.array([0.75, 0.75, 0.75, 0.75, 0.75,  0.75]))

    assert_equal(res['v(2)'], refv2)


def test_construct_circuit():
    cir = gnucap.Circuit()

    cir['V1'] = gnucap.VS(1,0,1.0)
    cir['R1'] = gnucap.R(1,2,'1k')
    cir['R2'] = gnucap.R(2,0,'3k')
    
    print gnucap.Circuit(simple_netlist)
    assert_equal(cir, gnucap.Circuit(simple_netlist))
