# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy as np
from pycircuit.utilities.param import Parameter, ParameterDict
import types

class Simulation(object):
    """Base class for simulations
    
    A Simulation object is responsible for setting up and running 
    analyses and finally providing the results gathered from them.

    """
    supported_analyses = []
    sim_options = []
    def __init__(self, circuit):
        self.circuit = circuit

        ## Init environment parameters
        self.epar = ParameterDict(Parameter("T", "Temperature", "K"))
        
        ## Init simulation option parameters
        self.options = ParameterDict(*self.sim_options)
        
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

class IParam(Parameter):
    """Instance parameter"""

class Variable(Parameter):
    """Design variable"""

class Sweep(object):
    """Parametric sweep of parameters"""
    def __init__(self, iter, pardict=None, parname=None):
        self.pardict = pardict
        self.parname = parname
        self.iter = iter.__iter__()

    def __iter__(self):
        return self

    def next(self):
        return self.iter.next()
    
class LinSweep(Sweep):
    """Linear sweep"""
    def __init__(self, start, stop, n, pardict=None, parname=None):
        self.start = start
        self.stop = stop
        self.n = n
        self.iter = np.linspace(start, stop, n).__iter__()

    @property
    def step(self):
        """Difference between adjacent steps"""
        return float(self.stop - self.start) / (self.n - 1.)

class LogSweep(Sweep):
    """Logarithmic sweep"""
    def __init__(self, start, stop, n=None, decade=None, 
                 pardict=None, parname=None):
        self.start = start
        self.stop = stop

        self.n = n
        self.decade = decade

        if n == None and decade == None or \
                n != None and decade != None:
            raise ValueError("Either n or decade must be != None")
        
        if decade != None:
            n = (np.log10(stop) - np.log10(start)) * decade + 1

        self.iter = np.logspace(np.log10(start), np.log10(stop), n).__iter__()

    @property
    def factor(self):
        """Ratio between adjacent steps"""
        if self.decade:
            n = (np.log10(self.stop) - np.log10(self.start)) * self.decade + 1
        else:
            n = self.n
        return (self.stop / self.start) ** (1. / (n-1))

def identify_sweep(iterator):
    """Identifies sweep from sweep values and return sweep instance"""
    values = list(iterator)
    
    
