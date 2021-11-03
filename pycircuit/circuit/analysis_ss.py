# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.utilities import Parameter, ParameterDict, isiterable
from pycircuit.circuit.analysis import CircuitResult, Analysis, remove_row_col
from pycircuit.circuit import Circuit, SubCircuit, VS,IS,R,C,L,Diode, gnd, \
    defaultepar, instjoin
import pycircuit.circuit.circuit
import pycircuit.circuit.symbolic as symbolic
from pycircuit.post.waveform import Waveform
from pycircuit.post.result import IVResultDict
from pycircuit.post.internalresult import InternalResultDict
from pycircuit.circuit.dcanalysis import DC

import pycircuit.circuit.numeric as numeric
import types


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

class SSAnalysis(Analysis):
    """Super class for small-signal analyses"""

    parameters = [Parameter(name='analysis', desc='Analysis name', default='ss'),
                   Parameter(name='dcx', desc='Provided DC-solution vector', 
                             unit='', 
                             default=None)]

    def __init__(self, cir, toolkit=None, **kvargs):    
        self.parameters = super(SSAnalysis, self).parameters + self.parameters            
        super(SSAnalysis, self).__init__(cir, **kvargs)
        

    def ss_map_function(self, func, ss, refnode):
        """Apply a function over a list of frequencies or a single frequency"""
        irefnode = self.cir.nodes.index(refnode)

        def myfunc(s):
            x = func(s)
            # Insert reference node voltage
            return self.toolkit.concatenate((x[:irefnode], self.toolkit.array([0.0]), x[irefnode:]))
            
        if isiterable(ss):
            return self.toolkit.array([myfunc(s) for s in ss]).swapaxes(0,1)
        else:
            return myfunc(ss)

    def dc_steady_state(self, freqs, refnode, complexfreq=False, u=None):
        """Return G,C,u matrices at dc steady-state and complex frequencies"""
        return dc_steady_state(self.cir, freqs, refnode, self.toolkit, 
                               complexfreq = complexfreq, u = u, 
                               analysis=self.par.analysis,
                               epar=self.epar,x0=self.par.dcx)

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

    parameters  = [Parameter(name='analysis', desc='Analysis name', 
                             default='ac')]

    def __init__(self, cir, toolkit=None, **kvargs):
        self.parameters = super(AC, self).parameters + self.parameters            
        super(AC, self).__init__(cir, **kvargs)

    def solve(self, freqs, refnode=gnd, complexfreq = False, u = None):
        G, C, CY, u, x, ss = self.dc_steady_state(freqs, refnode,
                                              complexfreq = complexfreq, u = u)

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.cir.get_node_index(refnode)
        G,C,CY,u = remove_row_col((G,C,CY,u), irefnode, self.toolkit)

        def acsolve(s):
            return self.toolkit.linearsolver(s*C + G, -u)

        xac = self.ss_map_function(acsolve, ss, refnode)

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

    parameters = [Parameter(name='analysis', desc='Analysis name', 
                            default='TransimpedanceAnalysis')]

    def __init__(self, cir, toolkit=None, **kvargs):
        parameters = super(TransimpedanceAnalysis, self).parameters + \
            self.parameters
            
        super(TransimpedanceAnalysis, self).__init__(cir, **kvargs)

            
    def solve(self, freqs, outbranches, currentoutput=False,
              complexfreq=False, refnode=gnd):
        toolkit = self.toolkit 

        n = self.cir.n
        x = self.toolkit.zeros(n) # This should be the x-vector at the DC operating point

        ## Complex frequency variable
        if complexfreq:
            s = freqs
        else:
            s = 2j*self.toolkit.pi*freqs

        epar = self.epar
        G = self.cir.G(x, epar)
        C = self.cir.C(x, epar)

        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.cir.get_node_index(refnode)
        G,C = remove_row_col((G,C), irefnode, self.toolkit)

        # Calculate the reciprocal G and C matrices
        Yreciprocal = G.T + s*C.T

        Yreciprocal = self.toolkit.toMatrix(Yreciprocal)

        result = []
        for branch in outbranches:
            ## Stimuli
            if currentoutput:
                u = zeros(n, dtype=int)
                ibranch = self.cir.get_branch_index(branch)
                u[ibranch] = -1
            else:
                u = self.toolkit.zeros(n, dtype=int)
                ## The signed is swapped because the u-vector appears in the lhs
                u[self.cir.get_node_index(branch.plus)] = -1
                u[self.cir.get_node_index(branch.minus)] = 1

            u, = remove_row_col((u,), irefnode, self.toolkit)

            ## Calculate transimpedances from currents in each nodes to output
            result.append(self.toolkit.linearsolver(Yreciprocal, -u))

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
    >>> print(res['Svnout'])
    (1.4904e-17+0j)
    >>> print(res['Svninp'])
    (1.4904e-15+0j)
    >>> print(res['gain'])
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

    parameters = [Parameter(name='analysis', desc='Analysis name', 
                            default='Noise'),
                  Parameter(name='inputsrc', desc='Input voltage source', 
                            unit='', 
                            default=None),
                  Parameter(name='outputnodes', 
                            desc='Output nodes (voltage output)', unit='', 
                            default=None),
                  Parameter(name='outputsrc', 
                            desc='Output voltage source (current output)',
                            unit='', 
                            default=None)]

    def __init__(self, cir, toolkit=None, **kvargs):
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

        self.parameters = super(Noise, self).parameters + self.parameters            
        super(Noise, self).__init__(cir, **kvargs)

    
        if not (self.par.outputnodes is not None or self.par.outputsrc is not None):
            raise ValueError('Output is not specified')
        elif self.par.outputnodes is not None and self.par.outputsrc is not None:
            raise ValueError('Cannot measure both output current and voltage '
                             'noise')
        
        if not (type(self.par.inputsrc) is types.StringType and \
                self.par.outputsrc is None or \
                    type(self.par.outputsrc) is types.StringType):
            raise ValueError('Sources must be given as instance names')

        self.inputsrc_name = self.par.inputsrc
        self.inputsrc = self.cir[self.par.inputsrc]
        self.outputnodes = self.par.outputnodes
        self.outputsrc_name = self.par.outputsrc
        if self.outputsrc_name:
            self.outputsrc = self.cir[self.par.outputsrc]

    def noise_map_function(self, func, ss, refnode):
        """Apply a function over a list of frequencies or a single frequency"""
        irefnode = self.cir.nodes.index(refnode)

        def myfunc(s):
            x, g = func(s)
            return x,g 
            
        if isiterable(ss):
            xlist = []
            glist = []
            for s in ss:
                x,g = myfunc(s)
                xlist = xlist + [x]
                glist = glist + [g]
            return self.toolkit.array(xlist),self.toolkit.array(glist) 
        else:
            return myfunc(ss)

    def solve(self, freqs, refnode=gnd, complexfreq = False, u = None):
        G, C, CY, u, x, ss = self.dc_steady_state(freqs, refnode,
                                              complexfreq = complexfreq, u = u)

        tk = self.toolkit
        
        def noisesolve(s):

            # Calculate the reciprocal G and C matrices
            Yreciprocal = G.T + s*C.T
            
            Yreciprocal2, uu = (tk.toMatrix(A) for A in (Yreciprocal, u))
            
            ## Calculate transimpedances from currents in each nodes to output
            zm =  tk.linearsolver(Yreciprocal2, -uu)

            xn2out = tk.dot(tk.dot(zm.reshape(1, tk.size(zm)), CY), tk.conj(zm))

            ## Etract gain
            gain = None
            if isinstance(self.inputsrc, VS):
                gain = self.cir.extract_i(zm, 
                                          instjoin(self.inputsrc_name, 'plus'),
                                          refnode=refnode, 
                                          refnode_removed=True)
                
            elif isinstance(self.inputsrc, IS):
                plus_node = instjoin(self.inputsrc_name, 'plus')
                minus_node = instjoin(self.inputsrc_name, 'minus')
                gain = self.cir.extract_v(zm, 
                                          self.cir.get_node(plus_node), 
                                          self.cir.get_node(minus_node), 
                                          refnode=refnode, refnode_removed=True)
            return xn2out[0], gain

        # Calculate output voltage noise
        if self.outputnodes is not None:
            ioutp, ioutn = (self.cir.get_node_index(node) 
                            for node in self.outputnodes)
            u[ioutp] = -1
            u[ioutn] = 1
        # Calculate output current noise
        else:
            plus_term = instjoin(self.outputsrc_name, 'plus')
            branch = self.cir.get_terminal_branch(plus_term)[0]
            ibranch = self.cir.get_branch_index(branch)
            u[ibranch] = -1
            
        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.cir.nodes.index(refnode)
        G,C,CY,u = remove_row_col((G,C,CY,u), irefnode, tk)
        
        xn2out, gain = self.noise_map_function(noisesolve, ss, refnode)

        # Store results
        result = InternalResultDict()

        if self.outputnodes is not None:
            result['Svnout'] = xn2out
        elif self.outputsrc is not None:
            result['Sinout'] = xn2out

        # Calculate the gain from the input voltage source by using the 
        # transimpedance vector to find the transfer from the branch voltage of
        # the input source to the output
        if isinstance(self.inputsrc, VS):
            result['gain'] = gain
            result['Svninp'] = xn2out / abs(gain)**2

        elif isinstance(self.inputsrc, IS):
            result['gain'] = gain
            result['Sininp'] = xn2out / abs(gain)**2

        return result


def dc_steady_state(cir, freqs, refnode, toolkit, complexfreq = False, 
                    analysis='ac', u = None, epar=defaultepar, x0=None):
    """Return G,C,CY,u matrices at dc steady-state and complex frequencies"""

    n = cir.n

    if complexfreq:
        ss = freqs
    else:
        ss = 2j*toolkit.pi*freqs

    if x0 is None:
        if toolkit.symbolic:
            x=None
        else:
        #x = zeros(n) ## FIXME, this should be calculated from the dc analysis
            resdc=DC(cir).solve()
            x = resdc.x
    else:
        x = x0 #provide the DC steady-state FIXME: need to add parameter to AC

    G = cir.G(x, epar)
    C = cir.C(x, epar)
    CY = cir.CY(x, toolkit.imag(ss), epar)

    ## Allow for custom stimuli, mainly used by other analyses
    if u is None:
        u = cir.u(x, analysis=analysis, epar=epar)

    return G, C, CY, u, x, ss


if __name__ == "__main__":
    import doctest
    doctest.testmod()
