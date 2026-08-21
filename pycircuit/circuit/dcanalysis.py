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
                  ## STAGE 13 -- solve with PCNR instead of limiting.
                  ##
                  ## Not the default, and that is a measured decision rather than
                  ## caution: see gates 13-4 and 13-5 in the transient work plan.
                  ## PCNR gives every limited quantity its own unknown, so devices
                  ## cannot interfere through a shared node; classic limiting is
                  ## what every existing test and circuit is tuned against.
                  Parameter(name='pcnr',
                            desc='Use Predictor/Corrector Newton-Raphson instead '
                                 'of limiting (Aadithya et al.); off by default',
                            unit='', default=False),
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
        
    def solve(self, x0=None):
        """Solve the DC operating point.

        ``x0`` is an optional starting guess.  STAGE 10.1 added it: a DC sweep
        needs to seed each point with the previous solution (continuation), and
        without a way in, every point of a sweep restarts from zeros and has to
        re-traverse the whole nonlinearity.  ``None`` keeps the historical
        behaviour exactly.
        """
        ## STAGE 8(d) -- see Circuit.reset_state.  A DC solve must not inherit a
        ## previous transient's history: it selected the wrong stamp and returned
        ## v(b) = 0.0 where 0.5 is correct.
        if hasattr(self.cir, 'reset_state'):
            self.cir.reset_state(self.epar)
        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node

        if x0 is None:
            x0 = self.toolkit.zeros(self.cir.n)
        else:
            x0 = self.toolkit.array(x0, dtype=float)
            if len(x0) != self.cir.n:
                raise ValueError(
                    'x0 has %d entries but the circuit has %d unknowns'
                    % (len(x0), self.cir.n))

        if self.par.pcnr:
            from pycircuit.circuit import pcnr as _pcnr
            if _pcnr.pcnr_junctions(self.cir):
                x, _v_lim, _its = _pcnr.solve_dc(
                    self.cir, self.cir.nodes[self.irefnode], x0=x0,
                    epar=self.epar, maxiter=self.par.maxiter,
                    reltol=self.par.reltol, abstol=self.par.vabstol)
                self.result = CircuitResult(self.cir, x)
                return self.result
            ## No participating device: PCNR has nothing to do, and falling
            ## through to the ordinary solver is the honest answer rather than
            ## raising -- the circuit simply has no limited quantities.

        def func(x):
            return self.cir.i(x, self.epar) + self.cir.u(0, analysis='dc', epar=self.epar), self.cir.G(x, self.epar)
            
        def source_callback(x, lambda_):
            f = self.cir.i(x, self.epar) + lambda_ * self.cir.u(0, analysis='dc', epar=self.epar)
            dFdx = self.cir.G(x, self.epar)
            return f, dFdx
            
        from pycircuit.circuit.nrsolver import (GminSteppingNewton,
                                                 JunctionGminSteppingNewton,
                                                 PseudoTransientNewton,
                                                 SourceSteppingNewton)
        from pycircuit.circuit.pcnr import pcnr_junctions

        ## P18 chain, physical-first: junction-gmin (the proper `gmin`,
        ## tracking the physical branch), then the diagonal/gshunt rescue,
        ## then source stepping.  Junction rows are reduced-system indices:
        ## the solve runs with the reference row removed, so rows above
        ## irefnode shift down by one (and the reference node itself cannot
        ## be a junction row it makes sense to perturb).
        _jrows = []
        for _i, _e, _ra, _rb in pcnr_junctions(self.cir):
            if self.irefnode in (_ra, _rb):
                continue
            _jrows.append((_ra - (_ra > self.irefnode),
                           _rb - (_rb > self.irefnode)))

        base_solver = self._get_nrsolver()
        jgmin_solver = JunctionGminSteppingNewton(base_solver, _jrows)
        gshunt_solver = GminSteppingNewton(jgmin_solver)
        source_chain = SourceSteppingNewton(gshunt_solver, refnode_removed(source_callback, self.irefnode, self.toolkit))
        ## P25: pseudo-transient continuation as the chain's LAST resort
        ## (industry order: gmin -> gshunt -> source stepping -> Psi-tc).
        ## Its pseudo steps are solved by the PLAIN base solver, never the
        ## chain: SourceSteppingNewton's rungs rebuild F from the callback
        ## WITHOUT the pseudo term, so handing the chain a deformed system
        ## would solve the wrong problem mid-ladder.
        solver_chain = PseudoTransientNewton(source_chain,
                                             rung_solver=base_solver)

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


class DCSweep(Analysis):
    """Sweep an instance parameter and solve the DC operating point at each value.

    STAGE 10.1.  This is SPICE's `.dc`, and it was the most conspicuous absence in
    the analysis inventory: there was no way to ask for a transfer curve, an I-V
    characteristic or a bias sweep without writing the loop by hand -- and a
    hand-written loop almost always restarts every point from zeros, because
    `DC.solve()` had no way to accept a starting guess until this item added one.

    CONTINUATION IS THE POINT, not a refinement.  Each solve is seeded with the
    previous point's solution, which is what makes a sweep across a nonlinearity
    converge at all: the step between adjacent points is small, so the previous
    answer is an excellent guess, whereas zeros is a cold start into the same
    exponential every time.  `continuation=False` is offered so the difference can
    be measured rather than asserted.

    >>> from pycircuit.circuit import numeric, gnd, SubCircuit
    >>> from pycircuit.circuit.elements import R, VS
    >>> import numpy as np
    >>> cir = SubCircuit(toolkit=numeric)
    >>> n = cir.add_node('a')
    >>> cir['V1'] = VS('a', gnd, v=0.0)
    >>> cir['R1'] = R('a', gnd, r=1e3)
    >>> res = DCSweep(cir, toolkit=numeric).solve('V1', 'v', np.linspace(0, 2, 3))
    >>> ['%.2f' % v for v in np.asarray(res.v('a', gnd), dtype=float)]
    ['0.00', '1.00', '2.00']
    """

    parameters = [Parameter(name='analysis', desc='Analysis name', default='dc')]

    def __init__(self, cir, toolkit=None, refnode=gnd, **kvargs):
        self.parameters = super(DCSweep, self).parameters + self.parameters
        super(DCSweep, self).__init__(cir, toolkit=toolkit, **kvargs)
        self.refnode = refnode
        self.irefnode = self.cir.get_node_index(refnode)

    def solve(self, instance, param, values, refnode=None, continuation=True):
        """Sweep ``cir[instance].ipar.<param>`` over ``values``.

        Returns a :class:`CircuitResult` whose sweep axis is ``values``, so
        ``res.v('out')`` is the swept curve.
        """
        import numpy

        if instance not in self.cir.elements:
            raise ValueError(
                '%r is not an instance in this circuit; have %s'
                % (instance, sorted(self.cir.elements)))
        element = self.cir[instance]
        if not hasattr(element.ipar, param):
            raise ValueError(
                '%r has no parameter %r; have %s'
                % (instance, param, sorted(p.name for p in element.instparams)))

        values = numpy.asarray(values, dtype=float)
        if values.size == 0:
            raise ValueError('values is empty; nothing to sweep')

        original = getattr(element.ipar, param)
        dc = DC(self.cir, toolkit=self.toolkit,
                refnode=self.refnode if refnode is None else refnode)

        columns = []
        x0 = None
        self.failures = []
        try:
            for value in values:
                setattr(element.ipar, param, float(value))
                self.cir.update_iparv()
                ## The previous solution, not zeros -- see the class docstring.
                res = dc.solve(x0=x0 if continuation else None)
                x = self.toolkit.array(res.x, dtype=float).reshape(-1)
                columns.append(x)
                x0 = x
        finally:
            ## Leave the circuit as it was found, whatever happened.  A sweep that
            ## silently leaves the last swept value behind would make every
            ## subsequent analysis on the same circuit depend on it -- exactly the
            ## defect stage 8(d) found in TLine.
            setattr(element.ipar, param, original)
            self.cir.update_iparv()

        X = numpy.array(columns).T
        self.result = CircuitResult(self.cir, x=X, xdot=None,
                                    sweep_values=values,
                                    sweep_label='%s.%s' % (instance, param),
                                    sweep_unit='')
        return self.result
