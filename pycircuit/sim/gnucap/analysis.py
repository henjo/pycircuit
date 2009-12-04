# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.sim
from pycircuit.sim import Simulation, LogSweep, LinSweep, Analysis
from pycircuit.utilities.param import Parameter

class SweptAnalysis(pycircuit.sim.Analysis):
    command = None
    sweep_param = None

    def run(self):
        self.sim.command('print %s v(nodes)'%self.command)

        cmdl = [self.command]

        paramvalue = getattr(self.par, self.sweep_param)
        if isinstance(paramvalue, LogSweep):
            cmdl += [paramvalue.start, paramvalue.stop, '*', paramvalue.factor]
        elif isinstance(paramvalue, LinSweep):
            cmdl += [paramvalue.start, paramvalue.stop, paramvalue.step]

        cmd = ' '.join([str(e) for e in cmdl])

        res = self.sim.command(cmd, parse_result=True)

        return res
    
class DCOP(Analysis):
    def __init__(self, sim, T=27):
        self.sim = sim
        self.T = T

    def run(self):
        self.sim.command('print op v(nodes)')
        return self.sim.command('op %s'%str(self.T), parse_result=True)

class Transient(SweptAnalysis):
    parameters = [Parameter('t', 'Time', 's')]
    sweep_param = 't'
    command = 'transient'

class AC(SweptAnalysis):
    parameters = [Parameter('freq', 'Frequency', 'Hz')]
    sweep_param = 'freq'
    command = 'ac'

## Register analyses in Simulation class
Simulation.supported_analyses = [DCOP, AC, Transient]
