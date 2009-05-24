from simulation import Sweep
from pycircuit.utilities.param import Parameter, ParameterDict

class Analysis(object):
    parameters = []
    def __init__(self, sim, **parvalues):
        self.sim = sim
        self.par = ParameterDict(*self.parameters, **parvalues)

class SweptAnalysis(Analysis):
    def __init__(self, sim, sweep, **parvalues):
        super(SweptAnalysis, self).__init__(sim, **parvalues)

        if isinstance(sweep, Sweep):
            self.sweep = sweep
        else:
            self.sweep = Sweep(sweep)
        
