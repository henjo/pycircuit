# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.utilities.param import Parameter, ParameterDict
import types

class Simulation(object):
    """Base class for simulations
    
    A Simulation object is responsible for setting up and running 
    analyses and finally providing the results gathered from them.

    """
    supported_analyses = []

    def __init__(self, circuit):
        self.circuit = circuit

        ## Init environment parameters
        self.epar = ParameterDict(Parameter("T", "Temperature", "K"))
        
        ## Init design variables
        self.var = ParameterDict()

        ## Add run_XXXX methods for each analysis
        def method_factory(analysis_class):
             def func(self, *args, **kvargs):
                ana = analysis_class(self, *args, **kvargs)
                return self.run_analysis(ana)

             name = 'run_%s'%analysis_class.__name__.lower()
                
             return name, types.MethodType(func, self, self.__class__)

        for analysis_class in self.supported_analyses:
            name, method = method_factory(analysis_class)
            setattr(self, name, method)
        
        self.analyses = []

    def clear(self):
        """Clear analyses"""
        self.analyses = []

    def add_analysis(self, analysis):
        """Add an analysis to simulation"""
        self.analyses.append(analysis)

    def run_analysis(self, analysis):
        """Run an analysis by analysis object"""
        self.clear()
        self.add_analysis(analysis)
        return self.run()

    def set_sweep(*sweeps):
        pass
    
    def run(self):
        """Run all analyses"""

class Circuit(object):
    """Base class for circuits"""
    pass

class Analysis(object):
    parameters = []
    def __init__(self, circuit, **parvalues):
        self.cir = circuit
        self.par = ParameterDict(*self.parameters, **parvalues)

class IParam(Parameter):
    """Instance parameter"""

class Variable(Parameter):
    """Design variable"""

class Sweep(object):
    """Parametric sweep of parameters"""
    def __init__(self, pardict, parname, iter):
        self.pardict = pardict
        self.parname = parname
        self.iter = iter
    

