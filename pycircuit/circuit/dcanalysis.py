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
                  ## ⚠ `pcnr=True` is a REQUEST, not a guarantee, and the
                  ## outcome is reported on the analysis object as
                  ## `pcnr_status`, which is always one of:
                  ##
                  ##   'off'              the parameter was False
                  ##   'used'             PCNR solved this
                  ##   'no-participants'  asked for, but no device in the
                  ##                      circuit declares a probe
                  ##   'fell-back'        PCNR raised; the ordinary Newton
                  ##                      chain solved it instead
                  ##
                  ## `pcnr_fell_back` is the older boolean and is kept, but
                  ## it cannot separate the middle two: both read False.
                  ## Check `pcnr_status` when it matters which happened.
                  Parameter(name='pcnr',
                            desc='Use Predictor/Corrector Newton-Raphson instead '
                                 'of limiting (Aadithya et al.); off by default. '
                                 'The outcome is reported as `pcnr_status`',
                            unit='', default=False),
                  ## 1e-6 since 2026-09-19 -- DC, Transient, JAXTransient and PSS share one
                  ## meaning and one default; the reason is at `Transient.vabstol`
                  Parameter(name='vabstol', 
                            desc='Absolute voltage error tolerance', unit='V', 
                            default=1e-6),
                  Parameter(name='maxiter', 
                            desc='Maximum number of iterations', unit='', 
                            default=100),
                  ## ROADMAP 12.3 -- the simulator-level anchor.
                  ##
                  ## SPICE's `GMIN`, and the same number SPICE3, ngspice and
                  ## A commercial simulator all default to -- and the same number
                  ## `compact.PspMosLongChannel` already carries privately as
                  ## `GLEAK`, which is the precedent item 12.3 was written to
                  ## generalise.  At one volt it is one picoamp, which is
                  ## exactly the default `iabstol`, so an anchor on a 1 V node
                  ## sits at the solver's own current noise floor by
                  ## construction rather than by argument.
                  ##
                  ## It is NOT in the matrix of a solve that succeeds.  The
                  ## anchor is a rescue: it engages only after the chain has
                  ## already raised `SingularMatrix`, and even then the answer
                  ## returned is normally from a final `gmin = 0` solve.  That
                  ## is what lets a weak-inversion measurement keep its
                  ## picoamps -- see `GminAnchorNewton`.
                  ##
                  ## `gmin=0` disables it and restores the pre-12.3 behaviour
                  ## exactly.
                  ##
                  ## ⚠⚠ THIS IS NOT SPICE'S STANDING `gmin`, DESPITE THE NAME
                  ## AND THE VALUE.  A commercial simulator inserts
                  ## `gmin = 1e-12 S` ACROSS EVERY NONLINEAR JUNCTION, in every
                  ## solve, and leaves it there.  Ours is a RESCUE ANCHOR:
                  ##
                  ##   - to GROUND from each node row, not ACROSS a junction;
                  ##   - engaged only after the chain has raised
                  ##     `SingularMatrix`, not unconditionally;
                  ##   - and absent from the returned answer, which normally
                  ##     comes from a final `gmin = 0` solve.
                  ##
                  ## ⚠ The two carry the SAME NAME and the SAME 1e-12, so
                  ## "gmin is 1e-12 on both sides" is a false reconciliation.
                  ## The difference is worth ~0.1% on a 1 GOhm hold node, and
                  ## unlike the `sqrt(2)` and `k_B` conventions recorded
                  ## elsewhere it changes the ANSWER rather than the units --
                  ## so a cross-tool comparison on any high-impedance node
                  ## (a sampled `kT/C` hold, a switched-capacitor bucket, an
                  ## oscillator tank) is comparing two different circuits until
                  ## the standing conductance is disabled on the other side.
                  ##
                  ## A real standing-`gmin` option is a ROADMAP ITEM (A.g),
                  ## wanted for parasitic-realistic circuits and DECIDED to
                  ## default to 0/off when it lands.  It is not this parameter
                  ## and must not be built by widening this one.
                  Parameter(name='gmin',
                            desc='Conductance to ground added to every node row '
                                 'to rescue a numerically empty row; 0 disables. '
                                 'NOT SPICE standing gmin -- see the note above',
                            unit='S', default=1e-12),
                  Parameter(name='bypass',
                            desc='Enable device model bypassing', unit='',
                            default=False),
                  Parameter(name='bypasstol',
                            desc='Bypass tolerance for device models', unit='V',
                            default=None),
                  Parameter(name='epar', desc='Environment parameters',
                            default=defaultepar),
                  ## THE SPICE BENCHMARK PLAN'S STAGE 3 (2026-10-07) -- SPICE's
                  ## `.ic` for a transient's operating point and `.nodeset`.  A
                  ## held node is not an unknown of the solve: its KCL row and
                  ## its column leave the system and its value enters as
                  ## lam * volts, lam the source-stepping factor -- SPICE's
                  ## row replacement (a unit row, `x_k = srcFact * ic`) in its
                  ## eliminated form, so no gmin ladder touches a held row.
                  ## A held node a voltage source or an inductor holds at DC
                  ## is refused (SPICE lets the source win, through a 1e10 S
                  ## pin carrying a meaningless current).  `pcnr` with held
                  ## nodes is refused.  None (the default) is today's solve,
                  ## untouched.
                  Parameter(name='pin',
                            desc="Node voltages held while the operating point is "
                                 "solved, {node: volts} (SPICE's .ic for a "
                                 "transient's operating point)",
                            unit='V', default=None),
                  Parameter(name='nodeset',
                            desc="Node voltages to start from, {node: volts}: a "
                                 "solve with them held, then one without them from "
                                 "its answer (SPICE's .nodeset)",
                            unit='V', default=None),
                  ]

    ## ROADMAP 12.3.  A class-level default so it can be read after ANY
    ## solve, including the PCNR branch that returns before the chain runs
    ## -- a flag that exists only on one path is a flag callers cannot use.
    gmin_anchor_retained = False

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

        The solve runs with ``epar.analysis_kind == 'dc'`` so elements whose
        stamps differ at DC (the ``Idt``/``Idtmod`` ic pin, idtmod.md sec. 5.1)
        can see it from ``G``/``i``/``u`` alike -- ``analysis='dc'`` reaches
        only ``u``.  Scoping notes on the shared ``analysis_kind`` helper.
        """
        with analysis_kind(self.epar, 'dc'):
            return self._solve_dc(x0)

    def _solve_dc(self, x0=None):
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

        ## PCNR OUTCOME, always set -- see `pcnr_status` below.  It used
        ## to be assigned only inside `if self.par.pcnr`, so a caller
        ## that checked `dc.pcnr_fell_back` on an ordinary solve got an
        ## AttributeError rather than an answer.
        self.pcnr_status = 'off'
        self.pcnr_fell_back = False
        pin, nodeset = self._held(self.par.pin), self._held(self.par.nodeset)
        if (pin or nodeset) and self.par.pcnr:
            raise ValueError('DC(pcnr=True) with held nodes (pin, nodeset) is not '
                             'supported')

        if self.par.pcnr:
            from pycircuit.circuit import pcnr as _pcnr
            ## Gate PARTICIPATION on the device records, not on the
            ## pnj-only pair view: that view exists for the gmin ladders
            ## and is empty for a circuit of pure fetlim/limvds devices,
            ## so `pcnr=True` on a MOSFET differential pair used to fall
            ## through to the ordinary solver SILENTLY (vector PCNR
            ## Stage 2, 2026-08-26).
            if _pcnr.pcnr_devices(self.cir):
                try:
                    x, _v_lim, _its = _pcnr.solve_dc(
                        self.cir, self.cir.nodes[self.irefnode], x0=x0,
                        epar=self.epar, maxiter=self.par.maxiter,
                        reltol=self.par.reltol, abstol=self.par.vabstol)
                    self.pcnr_status = 'used'
                    self.result = CircuitResult(self.cir, x)
                    return self.result
                except Exception as exc:               # noqa: BLE001
                    ## PCNR HAS NO RESCUE LADDER OF ITS OWN, and building one
                    ## was tried and measured (2026-08-26): a gmin shunt,
                    ## source stepping and Jacobian damping, each as an
                    ## adaptive ladder around `solve_dc`, all fail on a BJT
                    ## mirror from a 20 V start that the ordinary chain solves
                    ## in 146 evaluations.  The failure is PCNR's undamped
                    ## first step from a wild start, and it is the same on
                    ## every rung.  So a PCNR failure falls through to the
                    ## ordinary chain -- losing order-independence for THIS
                    ## solve, never the answer -- and says so.
                    logging.warning(
                        'DC(pcnr=True): PCNR failed (%s: %s); falling back to '
                        'the ordinary Newton chain for this solve',
                        type(exc).__name__, str(exc)[:80])
                    self.pcnr_status = 'fell-back'
                    self.pcnr_fell_back = True
            else:
                ## No participating device: PCNR has nothing to do, and
                ## falling through to the ordinary solver is the honest
                ## answer rather than raising -- the circuit simply has no
                ## limited quantities.
                ##
                ## ⚠ But it must SAY SO.  `pcnr_fell_back` stayed False
                ## here, which is the same value it has when PCNR ran and
                ## succeeded, so a caller could not tell "PCNR solved this"
                ## from "PCNR was asked for and did nothing".  That is the
                ## silent-fallthrough the note above records fixing for the
                ## diff-pair case, in a second guise: the gate was fixed and
                ## the REPORTING was not.
                self.pcnr_status = 'no-participants'
                logging.warning(
                    'DC(pcnr=True): no device in this circuit declares a '
                    'PCNR probe, so the ordinary Newton chain is used. '
                    'Check `dc.pcnr_status` if this is unexpected.')

        def func(x):
            return self.cir.i(x, self.epar) + self.cir.u(0, analysis='dc', epar=self.epar), self.cir.G(x, self.epar)
            
        def source_callback(x, lambda_):
            f = self.cir.i(x, self.epar) + lambda_ * self.cir.u(0, analysis='dc', epar=self.epar)
            dFdx = self.cir.G(x, self.epar)
            return f, dFdx

        if nodeset:
            ## SPICE's `.nodeset`: a solve with the nodes held, then one
            ## without them (but for the pins) from its answer.  A held
            ## solve that fails leaves the hint unused, as SPICE does.
            both = dict(nodeset)
            both.update(pin)
            try:
                x0 = self._solve_held(func, source_callback, x0, both)
            except (NoConvergenceError, SingularMatrix) as exc:
                logging.getLogger(__name__).warning(
                    'DC: the solve with the nodesets held failed (%s); solving without '
                    'them', str(exc)[:120])
        if pin:
            x = self._solve_held(func, source_callback, x0, pin)
            self.result = CircuitResult(self.cir, x)
            return self.result

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

        ## Node rows of the REDUCED system, for the anchor (see `_chain`).
        _node_rows = [i - (i > self.irefnode)
                      for i in range(len(self.cir.nodes)) if i != self.irefnode]
        anchored_chain = self._chain(
            refnode_removed(source_callback, self.irefnode, self.toolkit), _jrows,
            _node_rows)

        try:
            x = self._newton(func, x0, anchored_chain)
        except (NoConvergenceError, SingularMatrix) as last_e:
            logging.warning('Problems encountered: ' + str(last_e))
            raise last_e
        ## Reported rather than hidden: a caller that needs to know whether the
        ## number it just got still has 1 pS of invented conductance in it can
        ## ask, and a `True` here is the signal to give the node a real path.
        self.gmin_anchor_retained = anchored_chain.anchor_retained
        if anchored_chain.anchor_retained:
            logging.warning(
                'the DC answer was anchored: a gmin of %g S to ground was '
                'RETAINED because the unanchored system would not solve from '
                'the anchored point.  It passed both of the anchor\'s gates -- '
                'no unknown is missing from the Jacobian, and moving gmin a '
                'decade does not move the answer -- so it is a solution; but '
                'this circuit has a subnetwork with no DC path of its own, and '
                'gmin chose among its solutions' % self.par.gmin)

        self.result = CircuitResult(self.cir, x)
        return self.result

    def _chain(self, source_callback, jrows, node_rows):
        """The solver chain on the reduced system: `source_callback` its
        source-stepping residual, `jrows` its junction row pairs, `node_rows`
        its node (KCL) rows."""
        from pycircuit.circuit.nrsolver import (GminAnchorNewton,
                                                 GminSteppingNewton,
                                                 JunctionGminSteppingNewton,
                                                 PseudoTransientNewton,
                                                 SourceSteppingNewton)
        base_solver = self._get_nrsolver()
        jgmin_solver = JunctionGminSteppingNewton(base_solver, jrows)
        gshunt_solver = GminSteppingNewton(jgmin_solver)
        source_chain = SourceSteppingNewton(gshunt_solver, source_callback)
        ## P25: pseudo-transient continuation as the chain's LAST resort
        ## (industry order: gmin -> gshunt -> source stepping -> Psi-tc).
        ## Its pseudo steps are solved by the PLAIN base solver, never the
        ## chain: SourceSteppingNewton's rungs rebuild F from the callback
        ## WITHOUT the pseudo term, so handing the chain a deformed system
        ## would solve the wrong problem mid-ladder.
        solver_chain = PseudoTransientNewton(source_chain,
                                             rung_solver=base_solver)

        ## ROADMAP 12.3 -- OUTERMOST, and on a DIFFERENT exception from every
        ## layer below it.  The four ladders above engage on
        ## `NoConvergenceError`; a `SingularMatrix` passes straight through all
        ## of them by design (stage 6(a)), which is precisely why an empty row
        ## had no rescue at all.  The anchor takes that exception, and only
        ## that one.  Node rows of the REDUCED system: branch rows must not be
        ## anchored, because a conductance in a voltage source's KVL equation
        ## is not a leaky element, it is a wrong equation.
        return GminAnchorNewton(solver_chain, node_rows, gmin=self.par.gmin,
                                rung_solver=base_solver)

    def unholdable(self, nodes):
        """The nodes of `nodes` ({node: volts}) that cannot be held -- a
        voltage source or an inductor holds them at DC, so a hold would
        leave that branch's current undetermined (SPICE lets the source
        win).  Their names."""
        held = self._held(nodes)
        if not held:
            return []
        with analysis_kind(self.epar, 'dc'):
            red = _Held(self.cir.n, self.irefnode, held)
            return self._source_held(held, red.full(red.reduce(np.zeros(self.cir.n))))

    def _source_held(self, held, x):
        """The held nodes (`held`: {index: volts}) whose branch current only
        their own rows saw at `x`: nothing is left to fix it."""
        J = np.asarray(self.cir.G(x, self.epar))
        keep = [i for i in range(self.cir.n) if i != self.irefnode and i not in held]
        Jr = J[np.ix_(keep, keep)]
        nodes, out = len(self.cir.nodes), []
        for k in np.flatnonzero(~Jr.any(axis=0)):
            col = keep[k]
            if col >= nodes:
                out += [str(self.cir.nodes[r].name) for r in np.flatnonzero(J[:, col])
                        if r in held and str(self.cir.nodes[r].name) not in out]
        return out

    def _held(self, nodes):
        """`nodes` ({node or name: volts}) as {node index: volts}."""
        if not nodes:
            return {}
        out = {}
        for node, volts in dict(nodes).items():
            try:
                idx = self.cir.get_node_index(node)
            except ValueError:
                raise ValueError(f'a held node {node!r} is not in the circuit')
            if idx == self.irefnode:
                raise ValueError('the reference node cannot be held: it is 0 V by '
                                 'construction')
            out[idx] = float(volts)
        return out

    def _solve_held(self, func, source_callback, x0, held):
        """The operating point with the nodes `held` ({index: volts}) at
        their values (see the `pin` parameter)."""
        if not isinstance(self.toolkit.zeros(1), np.ndarray):
            raise NotImplementedError('held nodes need the numeric toolkit')
        from pycircuit.circuit.pcnr import pcnr_junctions
        red = _Held(self.cir.n, self.irefnode, held)
        nodes = len(self.cir.nodes)
        x0 = np.asarray(x0, dtype=float)
        by = self._source_held(held, red.full(red.reduce(x0)))
        if by:
            raise ValueError(
                f'held node(s) {", ".join(by)}: a voltage source or an inductor holds '
                'them at DC, so the hold leaves its current undetermined (SPICE lets the '
                'source win, through a 1e10 S pin)')
        pos = {int(i): k for k, i in enumerate(red.keep)}
        jrows = [(pos[ra], pos[rb]) for _i, _e, ra, rb in pcnr_junctions(self.cir)
                 if ra in pos and rb in pos]
        node_rows = [pos[i] for i in range(nodes) if i in pos]

        def held_func(xr):
            red.lam = 1.0
            return red.system(*func(red.full(xr)))

        def held_source(xr, lam):
            red.lam = lam
            return red.system(*source_callback(red.full(xr), lam))

        solver = self._chain(held_source, jrows, node_rows)
        n_branches = len(self.cir.branches)
        abstol = red.reduce(np.concatenate((self.par.iabstol * np.ones(nodes),
                                            self.par.vabstol * np.ones(n_branches))))
        xtol = red.reduce(np.concatenate((self.par.vabstol * np.ones(nodes),
                                          self.par.iabstol * np.ones(n_branches))))

        def limiter_func(xr, x0r):
            return red.reduce(self.cir.limit(red.full(xr), red.full(x0r), self.epar))

        names = reduced_row_names(self.cir, -1)
        try:
            x_res, _ = solver.solve_system(
                red.reduce(x0), held_func, self.toolkit, self.par.reltol, abstol, xtol,
                self.par.maxiter, limiter=limiter_func, scaler=self._get_scaler(),
                row_names=None if names is None else [names[i] for i in red.keep])
        except SingularMatrix:
            raise
        except NoConvergenceError as e:
            if 'Singular' in str(e) or 'linearsolver' in str(e).lower():
                raise SingularMatrix(str(e)) from e
            raise
        except LinAlgError as e:
            raise SingularMatrix(str(e)) from e
        red.lam = 1.0
        self.gmin_anchor_retained = solver.anchor_retained
        return red.full(x_res)

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

class _Held:
    """A solve with nodes held: the unknowns it leaves out -- the reference
    node at 0, each held node at `lam` times its volts (`lam` the
    source-stepping factor, as SPICE scales a pinned row) -- and the maps
    between the full vector and the solved one."""

    def __init__(self, n, irefnode, held):
        out = sorted({irefnode} | set(held))
        self.keep = np.array([i for i in range(n) if i not in set(out)], dtype=np.intp)
        self.out = np.array(out, dtype=np.intp)
        self.values = np.array([0.0 if i == irefnode else held[i] for i in out])
        self.n = n
        self.lam = 1.0

    def full(self, xr):
        x = np.empty(self.n)
        x[self.keep] = xr
        x[self.out] = self.lam * self.values
        return x

    def reduce(self, x):
        return np.asarray(x)[self.keep]

    def system(self, F, J):
        return np.asarray(F)[self.keep], np.asarray(J)[np.ix_(self.keep, self.keep)]


def refnode_removed(func, irefnode,toolkit):
    def new(x, *args, **kvargs):
        newx = insert_row(x, irefnode, toolkit)
        f, J = func(newx, *args, **kvargs)
        return remove_row_col((f, J), irefnode, toolkit)
    return new

if __name__ == "__main__":
    import doctest
    doctest.testmod()


## module level: a class-body comprehension cannot see a class attribute
_DC_PASSTHROUGH = ('reltol', 'iabstol', 'vabstol', 'maxiter', 'gmin',
                   'pcnr', 'bypass', 'bypasstol')


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

    ## ⚠ THE INNER DC'S TOLERANCES ARE DECLARED HERE AND PASSED THROUGH
    ## (2026-09-20, found by the peer suite after the `vabstol` default moved
    ## to 1e-6 in 039a017).  Until then a sweep built its DC with the DEFAULTS
    ## and REJECTED `vabstol` -- "parameter vabstol not in parameter
    ## dictionary" -- so a caller had no way to name a tolerance for a sweep at
    ## all, while the rule for that change was "name the tolerance, do not
    ## re-pin".  The list is DC's own Newton parameters, taken from `DC` so the
    ## two cannot drift apart; `analysis` and `epar` come from the base.
    DC_PASSTHROUGH = _DC_PASSTHROUGH
    parameters = ([Parameter(name='analysis', desc='Analysis name', default='dc')]
                  + [p for p in DC.parameters if p.name in _DC_PASSTHROUGH])

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
                refnode=self.refnode if refnode is None else refnode,
                epar=self.par.epar,
                **{name: getattr(self.par, name) for name in self.DC_PASSTHROUGH})

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
