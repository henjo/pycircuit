from pycircuit.sim.simulation import Sweep
from pycircuit.utilities.param import Parameter, ParameterDict

class Analysis(object):
    parameters = []
    def __init__(self, sim, **parvalues):
        self.sim = sim
        self.par = ParameterDict(*self.parameters, **parvalues)

