import logging

import numpy as np

from analysis import *

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
    >>> print np.around(res.v('net1'), 2)
    0.7

    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> n2 = c.add_node('net2')
    >>> c['is'] = IS(gnd, n1, i=57e-3)
    >>> c['R'] = R(n1, n2, r=1e1)
    >>> c['D'] = Diode(n2, gnd)
    >>> dc = DC(c)
    >>> res = dc.solve()
    >>> print np.around(res.v('net2'), 2)
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
                  Parameter(name='epar', desc='Environment parameters',
                            default=defaultepar)
                  ]

    def __init__(self, cir, toolkit=None, refnode=gnd, **kvargs):
        self.parameters = super(DC, self).parameters + self.parameters
        super(DC, self).__init__(cir, toolkit=toolkit, **kvargs)
        
        self.irefnode = self.cir.get_node_index(refnode)
        
    def solve(self):
        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node

        convergence_helpers = [self._simple, self._homotopy_gmin, 
                               self._homotopy_source, 
                               None]

        x0 = self.toolkit.zeros(self.cir.n) # Would be good with a better initial guess

        for algorithm in convergence_helpers:
            if algorithm == None:
                raise last_e
            else:
                if algorithm.__doc__:
                    logging.info('Trying ' + algorithm.__doc__)
                try:
                    x = algorithm(x0)
                except (NoConvergenceError, SingularMatrix), last_e:
                    logging.warning('Problems encoutered: ' + str(last_e))
                else:
                    break

        self.result = CircuitResult(self.cir, x)

        return self.result

    def _simple(self, x0):
        """Simple Newton's method"""
        def func(x):
            return self.cir.i(x) + self.cir.u(0,analysis='dc'), self.cir.G(x)

        return self._newton(func, x0)

    def _homotopy_gmin(self, x0):
        """Newton's method with gmin stepping"""
        x = x0
        for gmin in (1, 1e-1, 1e-2, 0):
            n_nodes = len(self.cir.nodes)
            Ggmin = self.toolkit.zeros((self.cir.n, self.cir.n))
            Ggmin[0:n_nodes, 0:n_nodes] = gmin * self.toolkit.eye(n_nodes)

            def func(x):
                return self.cir.i(x) + self.cir.u(0,analysis='dc'), \
                       self.cir.G(x) + Ggmin

            x, x0 = self._newton(func, x0), x

        return x

    def _homotopy_source(self, x0):
        """Newton's method with source stepping"""
        x = x0
        for lambda_ in (0, 1e-2, 1e-1, 1):
            def func(x):
                f = self.cir.i(x) + lambda_ * self.cir.u(0,analysis='dc')
                dFdx = self.cir.G(x)
                return f, dFdx            
            x, x0 = self._newton(func, x0), x

        return x

    def _newton(self, func, x0):
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))

        abstol = self.toolkit.concatenate((self.par.iabstol * ones_nodes,
                                 self.par.vabstol * ones_branches))
        xtol = self.toolkit.concatenate((self.par.vabstol * ones_nodes,
                                 self.par.iabstol * ones_branches))

        (x0, abstol, xtol) = remove_row_col((x0, abstol, xtol), self.irefnode, self.toolkit)

        try:
            result = fsolve(refnode_removed(func, self.irefnode,self.toolkit), 
                            x0, 
                            full_output = True, 
                            reltol = self.par.reltol,
                            abstol = abstol, xtol=xtol,
                            maxiter = self.par.maxiter,
                            toolkit = self.toolkit)
        except self.toolkit.linearsolverError(), e:
            raise SingularMatrix(e.message)

        x, infodict, ier, mesg = result

        if ier != 1:
            raise NoConvergenceError(mesg)

        # Insert reference node voltage
        return self.toolkit.concatenate((x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))

def refnode_removed(func, irefnode,toolkit):
    def new(x, *args, **kvargs):
        newx = toolkit.concatenate((x[:irefnode], toolkit.array([0.0]), x[irefnode:]))
        f, J = func(newx, *args, **kvargs)
        return remove_row_col((f, J), irefnode, toolkit)
    return new

if __name__ == "__main__":
    import doctest
    doctest.testmod()
