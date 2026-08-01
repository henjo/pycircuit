import logging

import numpy as np
from numpy.linalg import LinAlgError

from pycircuit.circuit.analysis import *

class DC(Analysis):
    """DC analyis class
    
    Linear circuit example:
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> dc = DC(c)
    >>> res = dc.solve()
    >>> res.v('net1')
    1.5

    Non-linear example:

    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> c['is'] = IS(gnd, n1, i=57e-3)
    >>> c['D'] = Diode(n1, gnd)
    >>> dc = DC(c)
    >>> res = dc.solve()
    >>> print(np.around(res.v('net1'), 2))
    0.7

    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> n2 = c.add_node('net2')
    >>> c['is'] = IS(gnd, n1, i=57e-3)
    >>> c['R'] = R(n1, n2, r=1e1)
    >>> c['D'] = Diode(n2, gnd)
    >>> dc = DC(c)
    >>> res = dc.solve()
    >>> print(np.around(res.v('net2'), 2))
    0.7

    """
    parameters = [Parameter(name='reltol', desc='Relative tolerance', unit='', 
                            default=1e-4),
                  Parameter(name='iabstol', 
                            desc='Absolute current eror tolerance', unit='A', 
                            default=1e-12),
                  Parameter(name='vabstol', 
                            desc='Absolute voltage error tolerance', unit='V', 
                            default=1e-12),
                  Parameter(name='maxiter', 
                            desc='Maximum number of iterations', unit='', 
                            default=100),
                  Parameter(name='bypass',
                            desc='Enable device model bypassing', unit='',
                            default=False),
                  Parameter(name='bypasstol',
                            desc='Bypass tolerance for device models', unit='V',
                            default=None),
                  Parameter(name='epar', desc='Environment parameters',
                            default=defaultepar)
                  ]

    def __init__(self, cir, toolkit=None, refnode=gnd, **kvargs):
        self.parameters = super().parameters + self.parameters
        super().__init__(cir, toolkit=toolkit, **kvargs)
        
        self.irefnode = self.cir.get_node_index(refnode)
        
    def solve(self):
        ## STAGE 8(d) -- see Circuit.reset_state.  A DC solve must not inherit a
        ## previous transient's history: it selected the wrong stamp and returned
        ## v(b) = 0.0 where 0.5 is correct.
        if hasattr(self.cir, 'reset_state'):
            self.cir.reset_state(self.epar)
        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node

        x0 = self.toolkit.zeros(self.cir.n) # Would be good with a better initial guess

        def func(x):
            return self.cir.i(x, self.epar) + self.cir.u(0, analysis='dc', epar=self.epar), self.cir.G(x, self.epar)
            
        def source_callback(x, lambda_):
            f = self.cir.i(x, self.epar) + lambda_ * self.cir.u(0, analysis='dc', epar=self.epar)
            dFdx = self.cir.G(x, self.epar)
            return f, dFdx
            
        from pycircuit.circuit.nrsolver import GminSteppingNewton, SourceSteppingNewton
        
        base_solver = self._get_nrsolver()
        gmin_solver = GminSteppingNewton(base_solver)
        solver_chain = SourceSteppingNewton(gmin_solver, refnode_removed(source_callback, self.irefnode, self.toolkit))

        try:
            x = self._newton(func, x0, solver_chain)
        except (NoConvergenceError, SingularMatrix) as last_e:
            logging.warning('Problems encountered: ' + str(last_e))
            raise last_e

        self.result = CircuitResult(self.cir, x)
        return self.result

    def _newton(self, func, x0, solver):
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))

        abstol = self.toolkit.concatenate((self.par.iabstol * ones_nodes,
                                 self.par.vabstol * ones_branches))
        xtol = self.toolkit.concatenate((self.par.vabstol * ones_nodes,
                                 self.par.iabstol * ones_branches))

        (x0, abstol, xtol) = remove_row_col((x0, abstol, xtol), self.irefnode, self.toolkit)

        def limiter_func(xr, x0r):
            x = self.toolkit.insert(xr, self.irefnode, 0.0)
            x0_full = self.toolkit.insert(x0r, self.irefnode, 0.0)
            
            x = self.cir.limit(x, x0_full, self.epar)
            return self.toolkit.concatenate((x[:self.irefnode], x[self.irefnode+1:]))

        try:
            scaler = self._get_scaler()
            x_res, _ = solver.solve_system(
                x0,
                refnode_removed(func, self.irefnode, self.toolkit),
                self.toolkit,
                self.par.reltol,
                abstol,
                xtol,
                self.par.maxiter,
                limiter=limiter_func,
                scaler=scaler,
                ## Stage 6: lets the solver name a node instead of a row index.
                row_names=reduced_row_names(self.cir, self.irefnode),
            )
        ## NARROW, deliberately -- see the matching note in `transient.py:_newton`.
        ## `except Exception` here reported every device-model bug as a convergence
        ## failure, which is the wrong diagnosis and points the reader at the bias
        ## point instead of at the traceback.  Only the solvers' own exceptions and
        ## genuine linear-algebra failures are translated; the rest propagate intact.
        except SingularMatrix:
            raise
        except NoConvergenceError as e:
            if 'Singular' in str(e) or 'linearsolver' in str(e).lower():
                raise SingularMatrix(str(e)) from e
            raise
        except LinAlgError as e:
            raise SingularMatrix(str(e)) from e

        # Insert reference node voltage
        return self.toolkit.concatenate((x_res[:self.irefnode], self.toolkit.array([0.0]), x_res[self.irefnode:]))

def refnode_removed(func, irefnode,toolkit):
    def new(x, *args, **kvargs):
        newx = toolkit.concatenate((x[:irefnode], toolkit.array([0.0]), x[irefnode:]))
        f, J = func(newx, *args, **kvargs)
        return remove_row_col((f, J), irefnode, toolkit)
    return new

if __name__ == "__main__":
    import doctest
    doctest.testmod()
