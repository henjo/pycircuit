from analysis import *
from numpy import array, delete, linalg, size, zeros, concatenate, pi
import sympy
from circuit import Circuit, SubCircuit, VS,R,C, gnd, VS, IS
from internalresult import InternalResultSet, InternalResult
from sympy import Symbol, Matrix, symbols, simplify
from types import TupleType

class SymbolicAC(Analysis):
    """Circuit analysis that calculates symbolic expressions of the unknowns

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['vs'] = VS(n1, gnd, v=Symbol('V'))
    >>> c['R'] = R(n1, gnd, r=Symbol('R'))
    >>> res = SymbolicAC(c).run()
    >>> res.getSignal('net1')
    V
    >>> res.getSignal('i0')
    -V/R
    """
    
    def solve(self, refnode=gnd):
        """Run a symbolic AC analysis with SymPy and store the results

        >>> c = SubCircuit()
        >>> n1 = c.addNode('net1')
        >>> c['vs'] = VS(n1, gnd, v=Symbol('V'))
        >>> c['R'] = R(n1, gnd, r=Symbol('R'))
        >>> SymbolicAC(c).solve()
        array([[V],
               [0.0],
               [-V/R]], dtype=object)


        """

        n=self.c.n()

        x = zeros((n,1)) # This should be the x-vector at the DC operating point

        G=self.c.G(x)
        C=self.c.C(x)
        U=self.c.U(x)

        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.nodes.index(refnode)
        G,C,U=removeRowCol((G,C,U), irefnode)

        G,C,U = (sympy.Matrix(A) for A in (G,C,U))

        outputvariables = map(Symbol, map(str, range(size(G,0))))
        resultdict =  sympy.solve_linear_system((Symbol('s')*C+G).row_join(-U), outputvariables)

        x = array([[resultdict[var] for var in outputvariables]]).T

        # Insert reference node voltage
        x = concatenate((x[:irefnode, :], array([[0.0]]), x[irefnode:,:]))
        return x

class SymbolicNoise(Analysis):
    """Circuit analysis that calculates symbolic expressions of the unknowns

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
    4*kT*R2**2/(R1+R2)
    >>> res.o.vn2in
    4*kT*(R1+R2)
    >>> simplify(res.o.gain-R2/(R1+R2))
    0
    
    """

    def __init__(self, circuit, inputsrc=None, outputnodes=None, outputsrc=None):
        """
        Initiate a symbolic noise analysis.

        @param     circuit: The circuit to be analyzed
        @type      circuit: Circuit instance
        @param    inputsrc: A voltage or current source in the circuit where the input noise should be referred to
        @type     inputsrc: VS or IS instance
        @param outputnodes: A tuple with the output nodes (outputpos outputneg)
        @param   outputsrc: A voltage source where the output current noise should be measured
        """
        Analysis.__init__(self, circuit)
    
        self.inputsrc = inputsrc
        self.outputnodes = outputnodes
        self.outputsrc = outputsrc

    def run(self, refnode=gnd):
        s = Symbol('s')
        kT = Symbol('kT')
        
        n=self.c.n()
        x = zeros((n,1)) # This should be the x-vector at the DC operating point
        G=self.c.G(x)
        C=self.c.C(x)
        CY=self.c.CY(x=x,kT=kT)

        # Calculate output voltage noise
        if self.outputnodes != None:
            U = zeros((n,1))
            ioutp, ioutn = (self.c.nodes.index(node) for node in self.outputnodes)
            U[ioutp] = -1.0
            U[ioutn] = 1.0

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
        resultdict =  sympy.solve_linear_system(Yreciprocal.row_join(-U), outputvariables)

        zm = sympy.Matrix([[resultdict[var] for var in outputvariables]]).T

        # Calculate output noise using correlation matrix
        xn2out = zm.T * CY * zm.applyfunc(sympy.conjugate)

        # Store results
        result = InternalResult()

        if self.outputnodes != None:
            result.storeSignal('vn2out', xn2out[0])
        elif self.outputsrc != None:
            result.storeSignal('in2out', xn2out[0])

        # Calculate the gain from the input voltage source by using the transimpedance vector
        # to find the transfer from the branch voltage of the input source to the output
        gain = None
        if isinstance(self.inputsrc, VS):
            i = self.c.getBranchIndex(self.inputsrc.branches[0])
            if i > irefnode:
                i -= 1
            gain = zm[i]
        
            result.storeSignal('gain', gain)
            result.storeSignal('vn2in', xn2out[0]/gain**2)

        return result
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
