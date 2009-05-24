# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.sim
from pycircuit.sim import Simulation, LogSweep, LinSweep, Analysis

class SweptAnalysis(pycircuit.sim.SweptAnalysis):
    command = None
    def __init__(self, circuit, sweep, **parvalues):
        super(SweptAnalysis, self).__init__(circuit, sweep, **parvalues)

        if not isinstance(sweep, (LogSweep, LinSweep)):
            raise ValueError('Unsupported sweep type')

    def run(self):
        self.sim.command('print %s v(nodes)'%self.command)

        cmdl = [self.command]

        if isinstance(self.sweep , LogSweep):
            cmdl += [self.sweep.start, self.sweep.stop, '*', self.sweep.factor]
        elif isinstance(self.sweep, LinSweep):
            cmdl += [self.sweep.start, self.sweep.stop, self.sweep.step]

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
    command = 'transient'

class AC(SweptAnalysis):
    command = 'ac'

## Register analyses in Simulation class
Simulation.supported_analyses = [DCOP, AC, Transient]
