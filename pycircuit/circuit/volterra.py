"""Volterra analysis based on the method described in
"Distortion analysis of analog integrated circuits" by Piet Wambacq and Willy Sansen

p129 chapter 5.2 depicts a flowchart describing the basic algorithm used in this module

"""

from numpy import array, delete, linalg, size, zeros, concatenate, pi, zeros, alltrue, maximum
from analysis import Analysis, DC, removeRowCol, isiterable
from pycircuit.param import Parameter, ParameterDict
from pycircuit.internalresult import InternalResultSet, InternalResult
from pycircuit.result import Waveform
import circuit
import sympy
from symbolicanalysis import symbolic_linsolve
from sympy import Symbol, Matrix, symbols, simplify, together, factor, cancel, diff, Mul, factorial
from symbolicelements import *
from copy import copy

def K(cir, x, ordervec, epar = circuit.defaultepar):
    """Calculate the taylor series term of the given order of i(x)

    Example:

    cir.K([1,1,0] will return the vector 1/(1!*1!) * d2I(X)/dX_0dX_1

    >>> from circuit import Node
    >>> d = Diode(Node('plus'), Node('minus'))
    >>> d.mpar.IS=Symbol('IS')
    >>> epar = copy(defaultepar)
    >>> epar.T = Symbol('T')
    >>> K(d, [Symbol('v1'), Symbol('v2')], [1,0], epar=epar)
    array([IS*qelectron*exp(qelectron*(v1 - v2)/(T*kboltzmann))/(T*kboltzmann),
           -IS*qelectron*exp(qelectron*(v1 - v2)/(T*kboltzmann))/(T*kboltzmann)], dtype=object)

    """

    xsyms = array([Symbol('x[%d]'%i) for i in range(size(x,0))])

    didx = cir.i(xsyms, epar=epar)
    for xsym, order in zip(xsyms, ordervec):
        if order > 0:
            didx = diff(didx, xsym, order)

    K = 1 / Mul(*map(factorial, ordervec))

    return array([K * expr.subs(zip(xsyms, x)) for expr in didx])

class Volterra(Analysis):
    """
    Volterra analysis class
    
    """

    @staticmethod
    def linearsolver(*args):
        return linalg.solve(*args)

    @staticmethod
    def toMatrix(array): return array.astype('complex')
        
    def run(self, freqs, **kvargs):
        x = self.solve(freqs, **kvargs)
        
        result = InternalResult()

        nodes = self.c.nodes
        for xvalue, node in zip(x[:len(nodes)], nodes):
            if isiterable(freqs):
                wave = Waveform(freqs, xvalue)
            else:
                wave = xvalue[0]
            result.storeSignal(self.c.getNodeName(node), wave)
        for i, data in enumerate(zip(x[len(nodes):], self.c.branches)):
            xvalue, branch = data
            if isiterable(freqs):
                wave = Waveform(freqs, xvalue)
            else:
                wave = xvalue[0]
            result.storeSignal('i' + str(i),wave)

        self.result = result

        return result

    def v(self, node1, node2):
        """Return the voltage v(node1, node2) from last run"""

        if self.result != None:
            return self.result[self.c.getNodeName(node1)] - \
                   self.result[self.c.getNodeName(node2)]

    def i(self, term):
        """Return terminal current i(term) of a circuit element from last run"""
        if self.result != None:
            branch, sign = self.c.getTerminalBranch(term)
            ibranch = self.c.branches.index(branch)
            return sign * self.result['i%d'%ibranch]
        
    def solve(self, freqs, refnode=gnd, complexfreq = False):
        n=self.c.n

        x = array([Symbol('V'), 0.0]) ## FIXME, this should be calculated from the dc analysis
        
        G=K(self.c, x,[1,0])
        C=self.c.C(x)
        U=self.c.U(x, analysis='ac')

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.getNodeIndex(refnode)
        G,C,U = removeRowCol((G,C,U), irefnode)

        G,C,U = (self.toMatrix(A) for A in (G,C,U))

        out = []

        if complexfreq:
            ss = freqs
        else:
            ss = 2j*pi*freqs

        def solvecircuit(s):
            solver = self.linearsolver

            x = solver(s*C + G, -U)
            
            # Insert reference node voltage
            return concatenate((x[:irefnode], array([0.0]), x[irefnode:]))

        if isiterable(freqs):
            out = [solvecircuit(s) for s in ss]
            # Swap frequency and x-vector dimensions
            return array(out).swapaxes(0,1)
        else:
            return solvecircuit(ss)

class SymbolicVolterra(Volterra):
    """
    Symbolic volterra analysis
    
    >>> cir = SubCircuit()

    >>> n1 = cir.addNode('n1')

    >>> cir['is'] = IS(n1, gnd, i=1e-3, iac=1)
    >>> cir['d'] = Diode(n1, gnd)
    >>> cir['d'].mpar.IS = Symbol('IS')
    >>> cir['c'] = C(n1, gnd, c=Symbol('C'))

    >>> volterra = SymbolicVolterra(cir)
    >>> res = volterra.run(freqs=array([Symbol('s')]), complexfreq=True)
    >>> volterra.result.getSignalNames()
    ['i0', 'gnd', 'net1']
    >>> volterra.result['n1']


    """
    numeric = False

    @staticmethod
    def linearsolver(A,b): return symbolic_linsolve(A,b)
    
    @staticmethod
    def toMatrix(array):
        return Matrix(array.tolist())

if __name__ == "__main__":
    import doctest
    doctest.testmod()
