# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import analysis
from analysis import Analysis, remove_row_col
from nportanalysis import TwoPortAnalysis
from numpy import array, delete, linalg, size, zeros, concatenate, pi, asarray
import circuit
from circuit import *
from copy import copy
from pycircuit.post.internalresult import InternalResultDict
import sympy
from sympy import Symbol, Matrix, symbols, simplify, together, factor, \
    cancel, exp, diff, Mul, factorial
from types import TupleType
from pycircuit.utilities.param import Parameter, ParameterDict

def symbolic_linsolve(A, b):
    """Numpy compatible wrapper around sympy.solve_linear_system"""

    A = sympy.Matrix(A.tolist())
    b = sympy.Matrix(b.tolist())

#    A.simplify(); b.simplify()

    res = array(A.inverse_ADJ() * b)

    return res.reshape((size(res,0),) )

class SymbolicAnalysis(Analysis):
    @staticmethod
    def linearsolver(*args):
        return symbolic_linsolve(*args)

    @staticmethod
    def toMatrix(array):
        return Matrix(array.tolist())

    @staticmethod
    def det(x):
        return Matrix(x).det()
    

class SymbolicAC(analysis.AC, SymbolicAnalysis):
    """Circuit analysis that calculates symbolic expressions of the unknowns

    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> c['vs'] = VS(n1, gnd, vac=Symbol('V'))
    >>> c['R'] = R(n1, gnd, r=Symbol('R'))
    >>> res = SymbolicAC(c).solve(Symbol('s'), complexfreq=True)
    >>> res['net1']
    V
    >>> res['i0']
    -V/R
    """

class SymbolicTransimpedanceAnalysis(SymbolicAnalysis,
                                     analysis.TransimpedanceAnalysis): pass

class SymbolicNoise(analysis.Noise):
    """Symbolic noise analysis that calculates input and output referred noise.
    
    The analysis is using the adjoint admittance matrix method to calculate the 
    transfers from each noise source to the output.
    
    Example, calculate input referred noise of a voltage divider:

    >>> c = SubCircuit()
    >>> kT = Symbol('kT')
    >>> R1=Symbol('R1', real=True)
    >>> R2=Symbol('R2', real=True)
    >>> n1,n2 = c.add_nodes('net1', 'net2')
    >>> c['vs'] = VS(n1, gnd, v=Symbol('V'))
    >>> c['R1'] = R(n1, n2, r=R1)
    >>> c['R2'] = R(n2, gnd, r=R2)
    >>> symnoise = SymbolicNoise(c, inputsrc=c['vs'], outputnodes=(n2, gnd))
    >>> res = symnoise.solve(Symbol('s'), complexfreq=True)
    >>> simplify(res['Svnout'])
    4*R1*R2*kT/(R1 + R2)
    >>> simplify(res['Svninp'])
    (4*R1*R2*kT + 4*kT*R1**2)/R2
    >>> simplify(res['gain'] - R2 / (R1 + R2))
    0
    
    """
    @staticmethod
    def linearsolver(*args):
        return symbolic_linsolve(*args)

    @staticmethod
    def toMatrix(array):
        return Matrix(array.tolist())

    def __init__(self, *args, **kvargs):
        analysis.Noise.__init__(self, *args, **kvargs)
        self.epar.append(Parameter('kT', default=Symbol('kT')))

class SymbolicTwoPortAnalysis(TwoPortAnalysis):
    """Analysis to find the symbolic 2-ports parameters of a circuit

    The transmission parameters are found as:

    A = v(inp, inn)/v(outp, outn) | io = 0
    B = v(inp, inn)/i(outp, outn) | vo = 0
    C = i(inp, inn)/v(outp, outn) | io = 0
    D = i(inp, inn)/i(outp, outn) | vo = 0

    >>> c = SubCircuit()
    >>> n1, n2 = c.add_nodes('net1', 'net2')
    >>> c['R1'] = R(n1, n2, r=Symbol('R1',real=True))
    >>> c['R2'] = R(n2, gnd, r=Symbol('R2',real=True))
    >>> symnoise = SymbolicTwoPortAnalysis(c, n1, gnd, n2, gnd, noise=True)
    >>> res = symnoise.solve(freqs = array([Symbol('s')]), complexfreq=True)
    >>> simplify(res['mu'].y[0])
    R2/(R1 + R2)
    >>> simplify(res['gamma'].y[0])
    1/R1
    >>> simplify(res['zeta'].y[0])
    R2
    >>> simplify(res['beta'].y[0])
    1
    >>> simplify(res['Svn'])
    (4*R1*R2*kT + 4*kT*R1**2)/R2
    >>> simplify(res['Sin'])
    4*kT/R2

    """
    
    def __init__(self, *args, **kvargs):
        TwoPortAnalysis.__init__(self, *args, **kvargs)
        self.epar.append(Parameter('kT', default=Symbol('kT')))

    ACAnalysis = SymbolicAC
    NoiseAnalysis = SymbolicNoise
    TransimpedanceAnalysis = SymbolicTransimpedanceAnalysis

if __name__ == "__main__":
    import doctest
    doctest.testmod()
