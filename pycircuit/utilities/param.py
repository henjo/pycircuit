# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import copy
from .misc import ObserverSubject

class Parameter(object):
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

class EvalError(Exception): pass

class ParameterDict(ObserverSubject):
    def __init__(self, *parameters, **kvargs):
        super().__init__()
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
            else:
                self._values[param.name] = param.default

        self.notify([param.name for param in parameters])
                
    def set(self, **kvargs):
        for k,v in kvargs.items():
            if k not in self._values:
                raise KeyError('parameter %s not in parameter dictionary'%k )
            self._values[k] = v
            
        self.notify(kvargs.keys())

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

    def __copy__(self):
        return self.copy()

    def eval_expressions(self, values, parameters=None, ignore_errors=False):
        """Evaluate expressions using parameter values from other ParameterDicts
        
        The function performs substitution of symbolic expressions using
        parameter values from other ParameterDict objects and returns
        a new ParameterDict with the evaluated expressions

        The *values* argument is a sequence of ParameterDict objects
        This allows for doing substutions from different parameter dictionaries. 
        The order sets the priority. First element has lowest priority.

        The *parameters* argument is a sequence of parameter names to update.
        If the value is None, all parameters will be updated.

        """
        out = ParameterDict(*self.parameters)

        ## Create a substition dictionary
        substdict = {}
        for paramdict in values:
            if paramdict is not None:
                for param in paramdict.parameters:
                    substdict[param.name] = paramdict.get(param)

        if parameters is None:
            parameters = self.keys()

        updated_values = {}
        for paramname in parameters:
            expr = self.get(paramname)

            if expr is not None:
                try:
                    if isinstance(expr, str):
                        value = eval(expr, substdict)
                    else:
                        value = expr
                except:
                    if not ignore_errors:
                        msg = "Can't evaluate %s (%s)" % (paramname, expr)
                        raise EvalError(msg)
                else:
                    setattr(out, paramname, value)

        return out

    def keys(self):
        return [param.name for param in self.parameters]
    
    def items(self):
        return [(param.name, getattr(self, param.name)) 
                for param in self.parameters]

    def update_values(self, d):
        """Update values from another paramdict"""
        self.set(**d._values)

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
            if key=='_parameters' and '_parameters' not in self.__dict__:
                raise AttributeError
            return self.__dict__[key]

    def __setattr__(self, key, value):
        if hasattr(self, '_parameters') and key in self._parameters:
            self._values[key] = value
            self.notify(args=(key,))
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
