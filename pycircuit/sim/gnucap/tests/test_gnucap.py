# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.sim import *
import pycircuit.sim.gnucap as gnucap
from pycircuit.post import Waveform, db20
from pycircuit.post.testing import *

import numpy as np

import logging
import unittest
from nose.tools import *

#logging.basicConfig(level = logging.DEBUG)

simple_netlist = """V1 1 0 1.0
R1 1 2 1k
R2 2 0 3k
"""

class GnucapTest(unittest.TestCase):
    direct = False
    
    def tearDown(self):
        del gnucap.session.GnucapSessionDirect.instance
        gnucap.session.GnucapSessionDirect.instance = None
    
    def test_start(self):
        sim = gnucap.Simulation(None, direct=self.direct)

        print sim.command('list')

    def test_netlist(self):
        sim = gnucap.Simulation(None, direct=self.direct)

        sim.send_netlist(simple_netlist)

        netlist = sim.command('list')

        assert_equal(netlist, 'V1 ( 1 0 )  DC  1.\nR1 ( 1 2 )  1.K\nR2 ( 2 0 )  3.K')

    def test_raw_dcop(self):
        sim = gnucap.Simulation(None, direct=self.direct)

        sim.send_netlist(simple_netlist)

        sim.command('print op v(nodes)')

        res = sim.command('op 27', parse_result=True)

        assert_equal(res['v(2)'].value(27), 0.75)

    def test_dcop(self):
        cir = gnucap.Circuit(simple_netlist)
        sim = gnucap.Simulation(cir, direct=self.direct)

        res = sim.run_dcop()

        assert_equal(res['v(2)'].value(27), 0.75)

    def test_ac(self):
        cir = gnucap.Circuit()
        cir['V1'] = gnucap.VS(1, 0, 1, ac=1)
        cir['R'] = gnucap.R(1,2,10e3)
        cir['C'] = gnucap.C(2,0,1e-9)

        sim = gnucap.Simulation(cir, direct=self.direct)

        res = sim.run_ac(LogSweep(1e4, 1e5, decade=10))

        v2 = res.v(2)

        freq = v2.xval()
        
        v2ref = 1/(1 + 2j*np.pi * freq * 10e3 * 1e-9)

        assert_waveform_almost_equal(abs(v2), abs(v2ref), 5)
        
    def test_transient(self):
        cir = gnucap.Circuit(simple_netlist)
        sim = gnucap.Simulation(cir, direct=self.direct)

        res = sim.run_transient(LinSweep(0,1,n=6))

        refv2 = Waveform(np.array([0.,0.2,0.4,0.6,0.8,1.]), 
                         np.array([0.75, 0.75, 0.75, 0.75, 0.75,  0.75]))

        assert_waveform_equal(res.v(2), refv2)

        assert_equal(res.v(2).xlabels, ('Time',))


    def test_construct_circuit(self):
        cir = gnucap.Circuit()

        cir['V1'] = gnucap.VS(1,0,1.0)
        cir['R1'] = gnucap.R(1,2,'1k')
        cir['R2'] = gnucap.R(2,0,'3k')

        print gnucap.Circuit(simple_netlist)
        assert_equal(cir, gnucap.Circuit(simple_netlist))

