# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import types

class Simulation(object):
    """Base class for simulations
    
    A Simulation object is responsible for setting up and running 
    analyses and finally providing the results gathered from them.

    """
    supported_analyses = []

    def __init__(self, circuit):
        self.circuit = circuit
        
        ## Add run_XXXX methods for each analysis
        def method_factory(analysis_class):
             def func(self, *args, **kvargs):
                ana = analysis_class(self, *args, **kvargs)
                return ana.run()

             name = 'run_%s'%analysis_class.__name__.lower()
                
             return name, types.MethodType(func, self, self.__class__)

        for analysis_class in self.supported_analyses:
            name, method = method_factory(analysis_class)
            setattr(self, name, method)
        
    def add_analysis(self, name, analysis):
        """Add an analysis"""
        
    def run_analysis(self, analysis):
        """Run an analysis either by name or by analysis object"""
        pass
    
    def run_all(self):
        """Run all stored analyses"""

class Circuit(object):
    """Base class for circuits"""
    pass

class Analysis(object):
    def __init__(self, circuit):
        self.cir = circuit


