from analysis import *
from numpy import array, delete, linalg, size, zeros, concatenate, pi, zeros
import sympy
from circuit import Circuit, SubCircuit, VS,R,C, gnd
from internalresult import InternalResultSet, InternalResult
from sympy import Symbol, Matrix, symbols, simplify

class Symbolic(Analysis):
    """Circuit analysis that calculates symbolic expressions of the unknowns

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['vs'] = VS(n1, gnd, v=Symbol('V'))
    >>> c['R'] = R(n1, gnd, r=Symbol('R'))
    >>> res = Symbolic(c).run()
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
        >>> Symbolic(c).solve()
        array([[V],
               [0.0],
               [-V/R]], dtype=object)


        """

        n=self.c.n()
        G=self.c.G(zeros((n,1)))
        C=self.c.C(zeros((n,1)))
        U=self.c.U()

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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
