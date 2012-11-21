import logging

from analysis import *
from .na import MNA

import numpy as np
import sympy

import symbolic

class SymbolicDC(Analysis):
    """Symbolic DC analyis class

    The SymbolicDC analysis finds the DC operating point symbolically. It works
    on linear or very simple non-linear circuits.
    """
    def __init__(self, cir, refnode=gnd, **kvargs):
        super(SymbolicDC, self).__init__(cir, toolkit=symbolic)
        
        self.na = MNA(cir, refnode=refnode, toolkit=symbolic)
        self.na.setup()
        
    def get_eqsys(self):
        """Return the equation system and variable list that gives the DC solution
        
        Returns eqsys, x
        """

        ## Update output vectors
        x = np.array([sympy.Symbol('x%d'%i) for i in range(self.na.n)])
        self.epar.analysis = 'dc'
        self.na.update(x, self.epar)
        return self.na.i + self.na.u, x

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

        except NotImplementedError, last_e:
            logging.error('Solver for equation %s not implemented in Sympy'%str(eqsys))

            raise last_e
        else:
            ## Get solution vector
            x = [sol[x_n] for x_n in x]
        
            self.result = CircuitResult(self.na, x)

            return self.result
