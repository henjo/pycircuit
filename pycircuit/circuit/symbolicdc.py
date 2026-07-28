import logging

from pycircuit.circuit.analysis import *

import numpy as np
import sympy

from .toolkit import symbolic

class SymbolicDC(Analysis):
    """Symbolic DC analyis class

    The SymbolicDC analysis finds the DC operating point symbolically. It works
    on linear or very simple non-linear circuits.
    """
    def __init__(self, cir, refnode=gnd, **kvargs):
        super().__init__(cir, toolkit=symbolic)
        
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

            if len(eqsys) == 1:
                sol = sympy.solve(eqsys[0], x[0])

                if len(sol) > 1:
                    raise NotImplementedError("Multiple solutions")

                sol = {x[0]: sol[0]}
            else:
                sol = sympy.solve(eqsys, *x)

                ## sympy.solve returns a dict for a linear system (a unique
                ## solution) but a list of solution tuples for a genuinely
                ## nonlinear one, even when there is only one solution. This
                ## branch was only ever exercised by linear circuits until
                ## _symbolic gained `tanh` -- nothing could build a
                ## multi-equation nonlinear system before that (see
                ## architecture.md P11).
                if isinstance(sol, list):
                    if len(sol) > 1:
                        raise NotImplementedError("Multiple solutions")
                    sol = dict(zip(x, sol[0]))

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
