# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy as np
from pycircuit import sim
from pycircuit.utilities import Parameter, ParameterDict, isiterable
from numpy import array, delete, linalg, size, zeros, concatenate, pi, \
    zeros, alltrue, maximum, conj, dot, imag, eye
from scipy import optimize
from pycircuit.circuit import Circuit, SubCircuit, VS,IS,R,C,L,Diode, gnd, \
    defaultepar, instjoin, circuit
import pycircuit.circuit.circuit
import symbolic
from pycircuit.post.waveform import Waveform
from pycircuit.post.result import IVResultDict
from pycircuit.post.internalresult import InternalResultDict
from copy import copy
import numeric
import types

np.set_printoptions(precision=4)

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

        if minus != None:
            ylabel = 'v(%s,%s)'%(str(plus), str(minus))
        else:
            ylabel = 'v(%s)'%(str(plus))

        return self.build_waveform(result, ylabel, 'V')

    def i(self, term):
        """Return terminal current i(term)"""
        result = self.circuit.extract_i(self.x, term, xdot = self.xdot)    
        return self.build_waveform(result, 'i(%s)'%(str(term)), 'A')

class CircuitResultAC(CircuitResult):
    """Result class for analyses that returns voltages and currents"""
    def __init__(self, circuit, xdcop, x, xdot = None,
                 sweep_values=[], sweep_label='', sweep_unit=''):
        super(CircuitResultAC, self).__init__(circuit, x, xdot,
                                              sweep_values=sweep_values,
                                              sweep_label=sweep_label,
                                              sweep_unit=sweep_unit)
        self.xdcop = xdcop

    def i(self, term):
        """Return terminal current i(term)"""
        result = self.circuit.extract_i(self.x, term, xdot = self.xdot, 
                                        linearized=True, 
                                        xdcop = self.xdcop)
        return self.build_waveform(result, 'i(%s)'%(str(term)), 'A')

def remove_row_col(matrices, n):
    result = []
    for A in matrices:
        for axis in range(len(A.shape)):
            A=delete(A, [n], axis=axis)
        result.append(A)
    return tuple(result)

class Analysis(sim.Analysis):
    parameters = [Parameter(name='analysis', desc='Analysis name', 
                            default=None),
                  Parameter(name='epar', desc='Environment parameters',
                            default=defaultepar)]
    def __init__(self, cir, toolkit=None, **kvargs):
        super(Analysis, self).__init__(cir, **kvargs)

        if toolkit == None:
            if cir.toolkit == None:
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
           xtol=1e-6, reltol=1e-4, abstol=1e-12):
    """Solve a multidimensional non-linear equation with Newton-Raphson's method

    In each iteration the linear system

    M{J(x_n)(x_{n+1}-x_n) + F(xn) = 0

    is solved and a new value for x is obtained x_{n+1}
    
    """
    
    converged = False
    ier = 2
    for i in xrange(maxiter):
        F, J = f(x0, *args) # TODO: Make sure J is never 0, e.g. by gmin (stepping)
        xdiff = linalg.solve(J, -F)# TODO: Limit xdiff to improve convergence

        x = x0 + xdiff

        if alltrue(abs(xdiff) < reltol * maximum(x, x0) + xtol):
            ier = 1
            mesg = "Success"
            break
        if alltrue(abs(F) < reltol * max(F) + abstol):
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

class SSAnalysis(Analysis):
    """Super class for small-signal analyses"""
    parameters = [Parameter(name='analysis', desc='Analysis name', 
                   default='ac')] + Analysis.parameters
    def ac_map_function(self, func, ss, refnode):
        """Apply a function over a list of frequencies or a single frequency"""
        irefnode = self.cir.nodes.index(refnode)

        def myfunc(s):
            x = func(s)
            # Insert reference node voltage
            return concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
            
        if isiterable(ss):
            return self.toolkit.array([myfunc(s) for s in ss]).swapaxes(0,1)
        else:
            return myfunc(ss)

    def dc_steady_state(self, freqs, refnode, toolkit, complexfreq=False, u=None, epar=defaultepar):
        """Return G,C,u matrices at dc steady-state and complex frequencies"""
        return dc_steady_state(self.cir, freqs, refnode, toolkit, 
                               complexfreq = complexfreq, u = u, 
                               analysis=self.par.analysis,
                               epar=epar)

class AC(SSAnalysis):
    """
    AC analysis class

    Examples:
    
    >>> circuit.default_toolkit = symbolic
    >>> from sympy import Symbol, simplify
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> c['vs'] = VS(n1, gnd, vac=Symbol('V'))
    >>> c['R'] = R(n1, gnd, r=Symbol('R'))
    >>> res = AC(c, toolkit=symbolic).solve(Symbol('s'), complexfreq=True)
    >>> res.v('net1')
    V
    >>> res.i('vs.plus')
    -V/R

    >>> circuit.default_toolkit = numeric
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> c['vs'] = VS(n1, gnd, vac=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> ac = AC(c)
    >>> res = ac.solve(freqs=array([1e6, 2e6]))
    >>> ac.result.v('net1')
    Waveform(array([ 1000000.,  2000000.]), array([ 1.5+0.j,  1.5+0.j]))
    >>> res.v(n1, gnd)
    Waveform(array([ 1000000.,  2000000.]), array([ 1.5+0.j,  1.5+0.j]))
    >>> res.i('vs.minus')
    Waveform(array([ 1000000.,  2000000.]), array([ 0.0015 +9.4248e-06j,  0.0015 +1.8850e-05j]))
    
    """
    
    def solve(self, freqs, refnode=gnd, complexfreq = False, u = None):
        G, C, u, x, ss = self.dc_steady_state(freqs, refnode, self.toolkit,
                                              complexfreq = complexfreq, u = u,
                                              epar = self.epar)

        def acsolve(s):
            return self.toolkit.linearsolver(s*C + G, -u)

        xac = self.ac_map_function(acsolve, ss, refnode)

        self.result = CircuitResultAC(self.cir, x, xac, ss * xac, 
                                      sweep_values = freqs, 
                                      sweep_label='frequency',
                                      sweep_unit='Hz')

        return self.result

class TransimpedanceAnalysis(SSAnalysis):
    """Calculates transimpedance or current-gain vector

    This function calculates the transimpedances from a current injected
    in every node to a voltage or branch current in the circuit. If 
    current=*True* the current gain is calculated.

    The outbranches is a list of Branch objects and in the voltage mode
    the output voltage is taken between the positive and negative node of
    the branch. The result is a list of transimpedance or current-gain
    vectors.

    Note, the element corresponding the reference node is eliminated in the
    result

    """
    def solve(self, freqs, outbranches, currentoutput=False,
              complexfreq=False, refnode=gnd):
        toolkit = self.toolkit 

        n = self.cir.n
        x = zeros(n) # This should be the x-vector at the DC operating point

        ## Complex frequency variable
        if complexfreq:
            s = freqs
        else:
            s = 2j*pi*freqs

        epar = self.epar
        G = self.cir.G(x, epar)
        C = self.cir.C(x, epar)

        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.cir.get_node_index(refnode)
        G,C = remove_row_col((G,C), irefnode)

        # Calculate the reciprocal G and C matrices
        Yreciprocal = G.T + s*C.T

        Yreciprocal = toolkit.toMatrix(Yreciprocal)

        result = []
        for branch in outbranches:
            ## Stimuli
            if currentoutput:
                u = zeros(n, dtype=int)
                ibranch = self.cir.get_branch_index(branch)
                u[ibranch] = -1
            else:
                u = zeros(n, dtype=int)
                ## The signed is swapped because the u-vector appears in the lhs
                u[self.cir.get_node_index(branch.plus)] = -1
                u[self.cir.get_node_index(branch.minus)] = 1

            u, = remove_row_col((u,), irefnode)

            ## Calculate transimpedances from currents in each nodes to output
            result.append(toolkit.linearsolver(Yreciprocal, -u))

        return result


class Noise(SSAnalysis):
    """Noise analysis that calculates input and output referred noise.
    
    The analysis is using the adjoint admittance matrix method to calculate the 
    transfers from each noise source to the output.
    
    Example, calculate input referred noise of a voltage divider:

    >>> circuit.default_toolkit = numeric
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> n2 = c.add_node('net2')
    >>> c['vs'] = VS(n1, gnd, vac=1.0)
    >>> c['R1'] = R(n1, n2, r=9e3)
    >>> c['R2'] = R(n2, gnd, r=1e3)
    >>> res = Noise(c, inputsrc='vs', outputnodes=(n2, gnd)).solve(0)
    >>> print res['Svnout']
    (1.4904e-17+0j)
    >>> print res['Svninp']
    (1.4904e-15+0j)
    >>> print res['gain']
    (0.1+0j)
    
    Symbolic example:
    
    >>> circuit.default_toolkit = symbolic
    >>> from sympy import Symbol, simplify
    >>> c = SubCircuit()
    >>> kT = Symbol('kT')
    >>> R1=Symbol('R1', real=True)
    >>> R2=Symbol('R2', real=True)
    >>> n1,n2 = c.add_nodes('net1', 'net2')
    >>> c['vs'] = VS(n1, gnd, v=Symbol('V'))
    >>> c['R1'] = R(n1, n2, r=R1)
    >>> c['R2'] = R(n2, gnd, r=R2)
    >>> noise = Noise(c, inputsrc='vs', outputnodes=(n2, gnd), toolkit=symbolic)
    >>> res = noise.solve(Symbol('s'), complexfreq=True)
    >>> simplify(res['Svnout'])
    4*R1*R2*kT/(R1 + R2)
    >>> simplify(res['Svninp'])
    (4*R1*R2*kT + 4*kT*R1**2)/R2
    >>> simplify(res['gain'] - R2 / (R1 + R2))
    0

    """

    parameters =  SSAnalysis.parameters + \
                  [Parameter(name='inputsrc', desc='Input voltage source', 
                            unit='', 
                            default=None),
                  Parameter(name='outputnodes', 
                            desc='Output nodes (voltage output)', unit='', 
                            default=None),
                  Parameter(name='outputsrc', 
                            desc='Output voltage source (current output)',
                            unit='', 
                            default=None)
                  ]
    def __init__(self, circuit, toolkit=None, **parvalues):
        """
        Initiate a noise analysis.

        Parameters
        ----------
        circuit : Circuit instance
            The circuit to be analyzed
        inputsrc : VS or IS instance
            A voltage or current source in the circuit where the input noise
            should be referred to
        outputnodes : tuple
            A tuple with the output nodes (outputpos outputneg)
        outputsrc: VS instance
            The voltage source where the output current noise is measured
        """

        Analysis.__init__(self, circuit, toolkit=toolkit, 
                          **parvalues)
    
        if not (self.par.outputnodes != None or self.par.outputsrc != None):
            raise ValueError('Output is not specified')
        elif self.par.outputnodes != None and self.par.outputsrc != None:
            raise ValueError('Cannot measure both output current and voltage '
                             'noise')
        
        if not (type(self.par.inputsrc) is types.StringType and \
                self.par.outputsrc == None or \
                    type(self.par.outputsrc) is types.StringType):
            raise ValueError('Sources must be given as instance names')

        self.inputsrc_name = self.par.inputsrc
        self.inputsrc = self.cir[self.par.inputsrc]
        self.outputnodes = self.par.outputnodes
        self.outputsrc_name = self.par.outputsrc
        if self.outputsrc_name:
            self.outputsrc = self.cir[self.par.outputsrc]

    def solve(self, freqs, refnode=gnd, complexfreq=False):
        toolkit = self.toolkit

        n = self.cir.n
        x = zeros(n) # This should be the x-vector at the DC operating point

        ## Complex frequency variable
        if complexfreq:
            s = freqs
        else:
            s = 2j*pi*freqs

        epar = self.epar
        G = self.cir.G(x, epar)
        C = self.cir.C(x, epar)
        CY = self.cir.CY(x, imag(s),epar)

        # Calculate output voltage noise
        if self.outputnodes != None:
            u = zeros(n, dtype=int)
            ioutp, ioutn = (self.cir.get_node_index(node) 
                            for node in self.outputnodes)
            u[ioutp] = -1
            u[ioutn] = 1
        # Calculate output current noise
        else:
            u = zeros(n, dtype=int)
            plus_term = instjoin(self.outputsrc_name, 'plus')
            branch = self.cir.get_terminal_branch(plus_term)[0]
            ibranch = self.cir.get_branch_index(branch)
            u[ibranch] = -1

        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.cir.nodes.index(refnode)
        G,C,u,CY = remove_row_col((G,C,u,CY), irefnode)

        # Calculate the reciprocal G and C matrices
        Yreciprocal = G.T + s*C.T

        Yreciprocal, u = (toolkit.toMatrix(A) for A in (Yreciprocal, u))

        ## Calculate transimpedances from currents in each nodes to output
        zm = toolkit.linearsolver(Yreciprocal, -u)

        xn2out = dot(dot(zm.reshape(1,size(zm)), CY), conj(zm))

        xn2out = xn2out[0]

        # Store results
        result = InternalResultDict()

        if self.outputnodes != None:
            result['Svnout'] = xn2out
        elif self.outputsrc != None:
            result['Sinout'] = xn2out

        # Calculate the gain from the input voltage source by using the 
        # transimpedance vector to find the transfer from the branch voltage of
        # the input source to the output
        gain = None
        if isinstance(self.inputsrc, VS):
            gain = self.cir.extract_i(zm, 
                                      instjoin(self.inputsrc_name, 'plus'),
                                      refnode=refnode, 
                                      refnode_removed=True)
            result['gain'] = gain
            result['Svninp'] = xn2out / abs(gain)**2

        elif isinstance(self.inputsrc, IS):
            plus_node = instjoin(self.inputsrc_name, 'plus')
            minus_node = instjoin(self.inputsrc_name, 'minus')
            gain = self.cir.extract_v(zm, 
                                      self.cir.get_node(plus_node), 
                                      self.cir.get_node(minus_node), 
                                      refnode=refnode, refnode_removed=True)
            result['gain'] = gain
            result['Sininp'] = xn2out / abs(gain)**2

        return result


def dc_steady_state(cir, freqs, refnode, toolkit, complexfreq = False, 
                    analysis='ac', u = None, epar=defaultepar):
    """Return G,C,u matrices at dc steady-state and complex frequencies"""

    n = cir.n

    x = zeros(n) ## FIXME, this should be calculated from the dc analysis

    G = cir.G(x)
    C = cir.C(x)

    ## Allow for custom stimuli, mainly used by other analyses
    if u == None:
        u = cir.u(x, analysis=analysis, epar=epar)

    ## Refer the voltages to the reference node by removing
    ## the rows and columns that corresponds to this node
    irefnode = cir.get_node_index(refnode)
    G,C,u = remove_row_col((G,C,u), irefnode)

    if complexfreq:
        ss = freqs
    else:
        ss = 2j*pi*freqs

    return G, C, u, x, ss

if __name__ == "__main__":
    import doctest
    doctest.testmod()
