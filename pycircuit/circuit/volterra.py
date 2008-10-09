"""Volterra analysis based on the method described in
"Distortion analysis of analog integrated circuits" by Piet Wambacq and Willy Sansen

p129 chapter 5.2 depicts a flowchart describing the basic algorithm used in this module

"""

from numpy import array, delete, linalg, size, zeros, concatenate, pi, zeros, alltrue, maximum
from analysis import Analysis, DC, removeRowCol, isiterable
from pycircuit.param import Parameter, ParameterDict
from pycircuit.internalresult import InternalResultSet, InternalResult
from pycircuit.result import Waveform
import sympy
from symbolicanalysis import SymbolicAC
from sympy import Symbol, Matrix, symbols, simplify, together, factor, cancel, diff, Mul, factorial
from symbolicelements import R, defaultepar, gnd, Diode, SubCircuit, IS, C
from copy import copy

def K(cir, x, ordervec, epar = defaultepar):
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
    Symbolic volterra analysis class
    
    >>> cir = SubCircuit()

    >>> n1 = cir.addNode('n1')

    >>> cir['is'] = IS(n1, gnd, i=1e-3, iac=1)
    >>> cir['d'] = Diode(n1, gnd)
    >>> cir['d'].mpar.IS = Symbol('IS')
    >>> cir['c'] = C(n1, gnd, c=Symbol('C'))

    >>> volterra = Volterra(cir)
    >>> res = volterra.run()
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
        
    def run(self, **kvargs):
        x = self.solve(**kvargs)

        result = InternalResult()

        return result

    def solve(self, refnode=gnd):
        x = zeros(self.c.n)
        
        ac = SymbolicAC(self.c)
        xac = ac.solve(freqs = Symbol('s'), refnode = refnode, complexfreq = True)

        print K(self.c, x, [2, 0])

if __name__ == "__main__":
    import doctest
    doctest.testmod()
