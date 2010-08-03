# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import copy
import sympy
import misc

class Parameter(sympy.Symbol):
    def __init__(self, name, desc=None, unit=None, default=None):
        self.name = name
        self.desc = desc
        self.default = default
        self.unit = unit

    def __eq__(self, a): 
        return self.__class__ is a.__class__ and self.name == a.name 

    def __hash__(self):
        return self.name.__hash__()

    def copy(self):
        return copy.copy(self)

    def __str__(self):
        return self.name

    def __repr__(self):
        kvargs = ('desc', 'unit', 'default')
        args = [repr(self.name)] + \
            [k + '=' + repr(getattr(self,k)) for k in kvargs
             if getattr(self,k)]
        return self.__class__.__name__ + '(' + ', '.join(args) + ')'

class ParameterDict(misc.ObserverSubject):
    def __init__(self, *parameters, **kvargs):
        super(ParameterDict, self).__init__()
        self._paramnames = []
        self._parameters = {}
        self._values = {}
        self.append(*parameters)
        self.set(**kvargs)

    def __eq__(self, a):
        return self._parameters == a._parameters
        
    def append(self, *parameters):
        for param in parameters:
            if param.name not in self._parameters:
                self._parameters[param.name] = param
                self._paramnames.append(param.name)
                self._values[param.name] = param.default

        self.notify()
                
    def set(self, **kvargs):
        for k,v in kvargs.items():
            if k not in self._values:
                raise KeyError('parameter %s not in parameter dictionary'%k )
            self._values[k] = v
            
        self.notify()

    def get(self, param):
        """Get value by parameter object or parameter name"""
        if isinstance(param, Parameter):
            return self._values[param.name]
        else:
            return self._values[param]

    def copy(self, *parameters, **kvargs):
        newpd = ParameterDict()
        newpd._values = copy.copy(self._values)
        newpd._parameters = copy.copy(self._parameters)
        newpd._paramnames = copy.copy(self._paramnames)
        newpd.append(*parameters)
        newpd.set(**kvargs)

        return newpd

    def eval_expressions(self, values):
        """Evaluate expressions using parameter values from other ParameterDicts
        
        The function performs substitution of symbolic expressions using
        parameter values from other ParameterDict objects.

        The *values* argument is a ParameterDict or a sequence of 
        (Parameter class, ParameterDict) tuples. This allows for doing 
        substutions of different kind of Parameters using different 
        paramter dictionaries. The order sets the priority.

        """
        out = ParameterDict(*self.parameters)

        ## Change values argument form to dictionary if needed
        if isinstance(values, ParameterDict):
            values = ((Parameter, values),)
        
        ## Create a substition dictionary
        substdict = {}
        for paramclass, paramdict in values:
            if paramdict != None:
                for param in paramdict.parameters:
                    substdict[param] = getattr(paramdict, param.name)

        for param, expr in self.items():
            if expr != None:
                try:
                    value = expr.subs(substdict)
                except AttributeError:
                    value = expr
                setattr(out, param, value)

        return out
    
    def items(self):
        return [(param.name, getattr(self, param.name)) 
                for param in self.parameters]

    def update(self, d):
        if isinstance(d, ParameterDict):
            self.__dict__.update(d.__dict__)
        else:
            self.__dict__.update(d)            

    def __getitem__(self, key):
        return self._parameters[key]

    def __setitem__(self, key, parameter):
        if key not in self._parameters:
            self.append(parameter)
        else:
            self._parameters[key] = parameter
            self._values[key] = parameter.default

    def __getattr__(self, key):
        if key != '_parameters' and hasattr(self, '_parameters') and \
                key in self._parameters:
            return self._values[key]
        else:
            return self.__dict__[key]

    def __setattr__(self, key, value):
        if hasattr(self, '_parameters') and key in self._parameters:
            self._values[key] = value
            self.notify()
        else:
            self.__dict__[key] = value
    
    def __contains__(self, key):
        if isinstance(key, Parameter):
            return key.name in self._parameters and \
                self._parameters[key.name] == key 
        return key in self._parameters
    
    def __len__(self):
        return len(self._parameters)

    @property
    def parameters(self):
        return [self._parameters[name] for name in self._paramnames]
