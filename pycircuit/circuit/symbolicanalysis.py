import analysis
from analysis import Analysis, removeRowCol
from nport import TwoPortAnalysis
from numpy import array, delete, linalg, size, zeros, concatenate, pi, asarray
import circuit
from circuit import *
from copy import copy
from pycircuit.internalresult import InternalResultSet, InternalResult
import sympy
from sympy import Symbol, Matrix, symbols, simplify, together, factor, cancel, exp, diff, Mul, factorial
from types import TupleType
from pycircuit.param import Parameter, ParameterDict

class NoSolutionFound(Exception):
    pass

def symbolic_linsolve(A, b):
    """Numpy compatible wrapper around sympy.solve_linear_system"""

    A = sympy.Matrix(A.tolist())
    b = sympy.Matrix(b.tolist())

    res = array(A.inv("ADJ") * b)

    return res.reshape((size(res,0),) )

class SymbolicAnalysis(Analysis):
    @staticmethod
    def linearsolver(*args):
        return symbolic_linsolve(*args)

    @staticmethod
    def toMatrix(array):
        return Matrix(array.tolist())

class SymbolicAC(analysis.AC):
    """Circuit analysis that calculates symbolic expressions of the unknowns

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['vs'] = VS(n1, gnd, vac=Symbol('V'))
    >>> c['R'] = R(n1, gnd, r=Symbol('R'))
    >>> res = SymbolicAC(c).run(Symbol('s'), complexfreq=True)
    >>> res['net1']
    V
    >>> res['i0']
    -V/R
    """
    @staticmethod
    def linearsolver(*args):
        return symbolic_linsolve(*args)

    @staticmethod
    def toMatrix(array):
        return Matrix(array.tolist())

class SymbolicNoise(analysis.Noise):
    """Symbolic noise analysis that calculates input and output referred noise.
    
    The analysis is using the adjoint admittance matrix method to calculate the transfers from
    each noise source to the output.
    
    Example, calculate input referred noise of a voltage divider:

    >>> c = SubCircuit()
    >>> kT = Symbol('kT')
    >>> R1=Symbol('R1', real=True)
    >>> R2=Symbol('R2', real=True)
    >>> n1,n2 = c.addNodes('net1', 'net2')
    >>> c['vs'] = VS(n1, gnd, v=Symbol('V'))
    >>> c['R1'] = R(n1, n2, r=R1)
    >>> c['R2'] = R(n2, gnd, r=R2)
    >>> res = SymbolicNoise(c, inputsrc=c['vs'], outputnodes=(n2, gnd)).run(Symbol('s'), complexfreq=True)
    >>> res['Svnout']
    4*R1*R2*kT/(R1 + R2)
    >>> res['Svninp']
    4*R1*kT*(-R1 - R2)**2/(R2*(R1 + R2))
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
    >>> n1, n2 = c.addNodes('net1', 'net2')
    >>> c['R1'] = R(n1, n2, r=Symbol('R1',real=True))
    >>> c['R2'] = R(n2, gnd, r=Symbol('R2',real=True))
    >>> res = SymbolicTwoPortAnalysis(c, n1, gnd, n2, gnd, noise=True).run(freqs = array([Symbol('s')]), complexfreq=True)
    >>> simplify(res['mu'].y[0])
    -R2/(-R1 - R2)
    >>> res['gamma'].y[0]
    1/R1
    >>> res['zeta'].y[0]
    R2
    >>> res['beta'].y[0]
    1
    >>> res['Svn']
    4*kT*R1*(R2+R1)/R2
    >>> res['Sin']
    4*kT/R2
    """
    
    ACAnalysis = SymbolicAC
    NoiseAnalysis = SymbolicNoise

if __name__ == "__main__":
    import doctest
    doctest.testmod()
