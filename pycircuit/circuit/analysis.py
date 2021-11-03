# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit import sim
from pycircuit.utilities import Parameter, ParameterDict, isiterable
from pycircuit.circuit import Circuit, SubCircuit, VS,IS,R,C,L,Diode, gnd, \
    defaultepar, instjoin, circuit
import pycircuit.circuit.circuit
import pycircuit.circuit.symbolic as symbolic
from pycircuit.post.waveform import Waveform
from pycircuit.post.result import IVResultDict
from pycircuit.post.internalresult import InternalResultDict
from copy import copy
import pycircuit.circuit.numeric as numeric
import types

class NoConvergenceError(Exception):
    pass

class SingularMatrix(Exception):
    pass

class CircuitResult(IVResultDict, InternalResultDict):
    """Result class for analyses that returns voltages and currents"""
    def __init__(self, circuit, x, xdot = None, 
                 sweep_values=[], sweep_label='', sweep_unit=''):
        super(CircuitResult, self).__init__()

        nodes = circuit.nodes

        self.circuit = circuit
        self.x = x
        self.xdot = xdot
        self.sweep_values = sweep_values
        self.sweep_label = sweep_label
        self.sweep_unit = sweep_unit

    def build_waveform(self, result, ylabel, yunit):
        if hasattr(result, '__iter__'):
            return Waveform(self.sweep_values, result,
                            ylabel = ylabel, yunit = 'V', 
                            xlabels = (self.sweep_label,), 
                            xunits=(self.sweep_unit,))
        else:
            return result

    def v(self, plus, minus=None):
        result = self.circuit.extract_v(self.x, plus, minus)

        if minus is not None:
            ylabel = 'v(%s,%s)'%(str(plus), str(minus))
        else:
            ylabel = 'v(%s)'%(str(plus))

        return self.build_waveform(result, ylabel, 'V')

    def i(self, term):
        """Return terminal current i(term)"""
        result = self.circuit.extract_i(self.x, term, xdot = self.xdot)    
        return self.build_waveform(result, 'i(%s)'%(str(term)), 'A')

def remove_row_col(matrices, n, toolkit):
    result = []
    for A in matrices:
        for axis in range(len(A.shape)):
            A=toolkit.delete(A, [n], axis=axis)
        result.append(A)
    return tuple(result)

class Analysis(sim.Analysis):
    parameters = [Parameter(name='analysis', desc='Analysis name', 
                            default=None),
                  Parameter(name='epar', desc='Environment parameters',
                            default=defaultepar)]

    def __init__(self, cir, toolkit=None, **kvargs):
        
        self.parameters = super(Analysis, self).parameters + self.parameters
        super(Analysis, self).__init__(cir, **kvargs)

        if toolkit is None:
            if cir.toolkit is None:
                toolkit = numeric
            else:
                toolkit = cir.toolkit

        self.toolkit = toolkit

        epar = self.par.epar

        if hasattr(toolkit, 'setup_analysis'):
            toolkit.setup_analysis(epar)

        self.cir = cir
        self.result = None
        self.epar = epar

def fsolve(f, x0, args=(), full_output=False, maxiter=200,
           xtol=1e-6, reltol=1e-4, abstol=1e-12, toolkit='Numeric'):
    """Solve a multidimensional non-linear equation with Newton-Raphson's method

    In each iteration the linear system

    M{J(x_n)(x_{n+1}-x_n) + F(xn) = 0

    is solved and a new value for x is obtained x_{n+1}
    
    """
    
    converged = False
    ier = 2
    for i in xrange(maxiter):
        F, J = f(x0, *args) # TODO: Make sure J is never 0, e.g. by gmin (stepping)
        xdiff = toolkit.linearsolver(J, -F)# TODO: Limit xdiff to improve convergence

        x = x0 + xdiff

        if toolkit.alltrue(abs(xdiff) < reltol * toolkit.maximum(x, x0) + xtol):
            ier = 1
            mesg = "Success"
            break
        if toolkit.alltrue(abs(F) < reltol * max(F) + abstol):
            ier = 1
            mesg = "Success"
            break
            
        x0 = x

    if ier == 2:
        mesg = "No convergence. xerror = "+str(xdiff)
    
    infodict = {}
    if full_output:
        return x, infodict, ier, mesg
    else:
        return x


if __name__ == "__main__":
    import doctest
    doctest.testmod()
