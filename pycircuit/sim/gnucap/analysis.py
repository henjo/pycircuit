import pycircuit.sim
from simulation import Simulation

class Analysis(pycircuit.sim.Analysis):
    command = None
    xlabel = None
    def __init__(self, sim, *args, **kvargs):
        self.sim = sim
        self.args = args
        self.kvargs = kvargs

    def run(self):
        self.sim.send('print %s v(nodes)'%self.command)
        cmdl = (self.command,) + tuple(str(arg) for arg in self.args) + \
            tuple((k + ' ' + str(v) for k,v in self.kvargs.items()))
        cmd = ' '.join(cmdl)
        return self.sim.send_command(cmd, [self.xlabel])
    
class DCOP(Analysis):
    def __init__(self, sim, T=27):
        self.sim = sim
        self.T = T

    def run(self):
        self.sim.send('print op v(nodes)')
        return self.sim.send_command('op %s'%str(self.T))

class Transient(Analysis):
    command = 'transient'
    xlabel = 'Time'

class AC(Analysis):
    command = 'ac'
    xlabel = 'Freq'

## Register analyses in Simulation class
Simulation.supported_analyses = [DCOP, AC, Transient]
