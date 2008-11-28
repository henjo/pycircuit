# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

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
from symbolicelements import R, defaultepar, gnd, Diode, SubCircuit, IS, C, VCCS
from copy import copy

class NLVCCS(VCCS):
    linear = False
    def i(self, x, epar=defaultepar):
        """
        
        """
        v = x[0]-x[1]
        I = self.ipar.gm * v + 2 * v**2 + 3 * v**3
        return array([0,0, I, -I])

def product(self, factors):
    return Mul(*factors)

def K(cir, x, ordervec, epar = defaultepar):
    """Calculate the taylor series term of the given order of i(x)

    Example:

    cir.K([1,1,0] will return the vector 1/(1!*1!) * d2I(X)/dX_0dX_1

    >>> from circuit import Node
    >>> d = NLVCCS(Node('plus'), Node('minus'), Node('plus'), Node('minus'), gm=Symbol('gm'))
    >>> epar = copy(defaultepar)
    >>> epar.T = Symbol('T')
    >>> K(d, [0,0], [1,0], epar=epar)
    array([0, 0, gm, -gm], dtype=object)
    >>> K(d, [0,0], [2,0], epar=epar)
    array([0, 0, 2, -2], dtype=object)
    >>> K(d, [0,0], [3,0], epar=epar)
    array([0, 0, 3, -3], dtype=object)
    >>> K(d, [0,0], [0,2], epar=epar)
    array([0, 0, 2, -2], dtype=object)

    """

    ## Generate a list of symbols like [x[0], x[1] ..]
    xsyms = array([Symbol('x[%d]'%i) for i in range(size(x,0))])

    ## Calculate derivative
    didx = cir.i(xsyms, epar=epar)
    for xsym, order in zip(xsyms, ordervec):
        if order > 0:
            didx = diff(didx, xsym, order)

    ## Calculate taylor coefficient
    K = 1 / product(map(factorial, ordervec))

    ## Substitute x[k] with given x values
    return array([K * expr.subs(zip(xsyms, x)) for expr in didx])

class Volterra(Analysis):
    """
    Symbolic volterra analysis class
    
    >>> cir = SubCircuit()

    >>> n1 = cir.addNode('n1')

    >>> cir['is'] = IS(n1, gnd, i=1e-3, iac=1)
    >>> cir['x'] = NLVCCS(n1, gnd, n1, gnd, gm=Symbol('gm'))
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

        ## Find non-linear elements
        nlelements = [e for e in self.c.xflatelements if not e.linear]
        print self.c.n, nlelements
#        print K(self.c, x, [2, 0]), nlelements

if __name__ == "__main__":
    import doctest
    doctest.testmod()
