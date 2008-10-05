import analysis
from analysis import Analysis, removeRowCol, TwoPort
from numpy import array, delete, linalg, size, zeros, concatenate, pi
import sympy
from circuit import Circuit, SubCircuit, VS,R,C, gnd, VS, IS
from pycircuit.internalresult import InternalResultSet, InternalResult
from sympy import Symbol, Matrix, symbols, simplify, together, factor, cancel
from types import TupleType
from pycircuit.param import Parameter, ParameterDict

class NoSolutionFound(Exception):
    pass

def symbolic_linsolve(A, b):
    """Numpy compatible wrapper around sympy.solve_linear_system"""

    A = sympy.Matrix(A)
    b = sympy.Matrix(b)

    outputvariables = map(Symbol, map(str, range(size(b))))
    resultdict =  sympy.solve_linear_system(A.row_join(b), *outputvariables)
    return array([[resultdict[var] for var in outputvariables]]).T


class SymbolicAC(analysis.AC):
    """Circuit analysis that calculates symbolic expressions of the unknowns

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['vs'] = VS(n1, gnd, v=Symbol('V'))
    >>> c['R'] = R(n1, gnd, r=Symbol('R'))
    >>> res = SymbolicAC(c).run(Symbol('s'), complexfreq=True)
    >>> res.getSignal('net1')
    V
    >>> res.getSignal('i0')
    -V/R
    """

    numeric = False

    @staticmethod
    def linearsolver(A,b): return symbolic_linsolve(A,b)


class SymbolicNoise(Analysis):
    """Symbolic noise analysis that calculates input and output referred noise.
    
    The analysis is using the adjoint admittance matrix method to calculate the transfers from
    each noise source to the output.
    
    Example, calculate input referred noise of a voltage divider:

    >>> c = SubCircuit()
    >>> kT = Symbol('kT')
    >>> R1=Symbol('R1')
    >>> R2=Symbol('R2')
    >>> n1 = c.addNode('net1')
    >>> n2 = c.addNode('net2')
    >>> c['vs'] = VS(n1, gnd, v=Symbol('V'))
    >>> c['R1'] = R(n1, n2, r=R1)
    >>> c['R2'] = R(n2, gnd, r=R2)
    >>> res = SymbolicNoise(c, inputsrc=c['vs'], outputnodes=(n2, gnd)).run()
    >>> res.o.vn2out
    4*R1*R2*kT/(R1 + R2)
    >>> res.o.vn2in
    4*R1*kT/R2*(R1 + R2)
    >>> simplify(res.o.gain - R2 / (R1 + R2))
    0
    
    """

    def __init__(self, circuit, inputsrc=None, outputnodes=None, outputsrc=None):
        """
        Initiate a symbolic noise analysis.

        Parameters
        ----------
        circuit : Circuit instance
            The circuit to be analyzed
        inputsrc : VS or IS instance
            A voltage or current source in the circuit where the input noise should be referred to
        outputnodes : tuple
            A tuple with the output nodes (outputpos outputneg)
        outputsrc: VS instance
            The voltage source where the output current noise is measured
        """

        Analysis.__init__(self, circuit)
    
        if not (outputnodes != None or outputsrc != None):
            raise ValueError('Output is not specified')
        elif outputnodes != None and outputsrc != None:
            raise ValueError('Cannot measure both output current and voltage noise')
        
        self.inputsrc = inputsrc
        self.outputnodes = outputnodes
        self.outputsrc = outputsrc

    def run(self, refnode=gnd):
        
        s = Symbol('s')
        kT = Symbol('kT')
        
        ## Set environment parameters
        epar = ParameterDict(Parameter('kT', default=kT))
        
        n = self.c.n
        x = zeros((n,1)) # This should be the x-vector at the DC operating point

        G = self.c.G(x, epar)
        C = self.c.C(x, epar)
        CY = self.c.CY(x, epar)

        # Calculate output voltage noise
        if self.outputnodes != None:
            U = zeros((n,1))
            ioutp, ioutn = (self.c.getNodeIndex(node) for node in self.outputnodes)
            U[ioutp] = -1.0
            U[ioutn] = 1.0
        # Calculate output current noise
        else:
            U - zeros((n,1))
            ibranch = self.c.getBranchIndex(self.outputsrc.branch)
            U[ibranch] = 1.0

        ## Convert to Sympy matrices
        G,C,U,CY = (sympy.Matrix(A) for A in (G, C, U, CY))

        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.nodes.index(refnode)
        for A in (G,C,U,CY):
            A.row_del(irefnode)
        for A in (G,C,CY):
            A.col_del(irefnode)

        # Calculate the reciprocal G and C matrices
        Yreciprocal = G.T + s*C.T

        ## Calculate transimpedances from currents in each nodes to output
        outputvariables = map(Symbol, map(str, range(size(G,0))))
        resultdict =  sympy.solve_linear_system(Yreciprocal.row_join(-U), *outputvariables)

        if resultdict == None:
            raise NoSolutionFound()            
        
        ## Collect transimpedance vector from result dictionary
        zm = sympy.Matrix([[resultdict[var] for var in outputvariables]]).T

        ## Simplify
        zm = zm.applyfunc(lambda x: cancel(together(x)))

        # Calculate output noise using correlation matrix
        xn2out = zm.T * CY * zm.applyfunc(sympy.conjugate)
        xn2out = cancel(together(xn2out[0]))

        # Store results
        result = InternalResult()

        if self.outputnodes != None:
            result.storeSignal('vn2out', xn2out)
        elif self.outputsrc != None:
            result.storeSignal('in2out', xn2out)

        # Calculate the gain from the input voltage source by using the transimpedance vector
        # to find the transfer from the branch voltage of the input source to the output
        gain = None
        if isinstance(self.inputsrc, VS):
            i = self.c.getBranchIndex(self.inputsrc.branches[0])
            if i > irefnode:
                i -= 1
            gain = zm[i]
        
            result.storeSignal('gain', gain)
            result.storeSignal('vn2in', xn2out/gain**2)

        return result

class SymbolicTwoPort(TwoPort):
    """Analysis to find the symbolic 2-ports parameters of a circuit

    The transmission parameters are found as:

    A = v(inp, inn)/v(outp, outn) | io = 0
    B = v(inp, inn)/i(outp, outn) | vo = 0
    C = i(inp, inn)/v(outp, outn) | io = 0
    D = i(inp, inn)/i(outp, outn) | vo = 0

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> n2 = c.addNode('net2')
    >>> c['R1'] = R(n1, n2, r=Symbol('R1'))
    >>> c['R2'] = R(n2, gnd, r=Symbol('R2'))
    >>> res = SymbolicTwoPort(c, n1, gnd, n2, gnd).run(freqs = array([Symbol('s')]), complexfreq=True)
    >>> simplify(res.getSignal('mu').y[0])
    -R2/(-R1 - R2)
    >>> res.getSignal('gamma').y[0]
    1/R1
    >>> res.getSignal('zeta').y[0]
    R2
    >>> res.getSignal('beta').y[0]
    1
    
    """
    
    ACAnalysis = SymbolicAC
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
