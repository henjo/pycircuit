import logging

from pycircuit.circuit.analysis import *

import numpy as np
import sympy

import pycircuit.circuit.symbolic as symbolic

class SymbolicDC(Analysis):
    """Symbolic DC analyis class

    The SymbolicDC analysis finds the DC operating point symbolically. It works
    on linear or very simple non-linear circuits.
    """
    def __init__(self, cir, refnode=gnd, **kvargs):
        super(SymbolicDC, self).__init__(cir, toolkit=symbolic)
        
        self.irefnode = self.cir.get_node_index(refnode)

    def get_eqsys(self):
        """Return the equation system and variable list that gives the DC solution
        
        Returns eqsys, x
        """
        
        x = np.array([sympy.Symbol('x%d'%i) for i in range(self.cir.n)])

        eqsys = self.cir.i(x, epar=self.epar) + self.cir.u(0, analysis='dc', epar=self.epar)

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        return [eq.subs(x[self.irefnode], 0) for eq in np.delete(eqsys, self.irefnode)], \
            np.delete(x, self.irefnode)

    def solve(self):
        try:
            eqsys, x = self.get_eqsys()

            ## Currently sympy doesn't support non-linear eq. systems only single equation
            if len(eqsys) == 1:
                sol = sympy.solve(eqsys[0], x[0])

                if len(sol) > 1:
                    raise NotImplemented("Multiple solutions")

                sol = {x[0]: sol[0]}
            else:
                sol = sympy.solve(eqsys, *x)

        except NotImplementedError as last_e:
            logging.error('Solver for equation %s not implemented in Sympy'%str(eqsys))

            raise last_e
        else:
            ## Get solution vector
            x = [sol[x_n] for x_n in x]
        
            ## Insert reference node
            x.insert(self.irefnode, 0)

            self.result = CircuitResult(self.cir, np.array(x))

            return self.result
