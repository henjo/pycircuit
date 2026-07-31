# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import contextlib
import warnings

from numpy.linalg import LinAlgError

from pycircuit.circuit.analysis import *
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.dcanalysis import refnode_removed
## The clamp the step controller applies to every accepted step.  The force-accept
## path in `solve()` is the one place that used to bypass it, and 4b's whole point
## is that it must not: one bound, named once.  `stepcontroller` imports nothing
## from this package, so this is import-safe at module level.
from pycircuit.circuit.stepcontroller import MAX_GROWTH_RATIO

## STAGE 2a -- BLAS thread control, discovered rather than required.
##
## Circuit matrices are small (n ~ 10^2), so a threaded LAPACK spends more time
## spawning and synchronising threads than doing the 1.7 MFLOP of work.  Measured
## on the 139-unknown leapfrog: the whole transient runs **1.72x faster** with BLAS
## limited to one thread, on a 4-core box.
##
## `threadpoolctl` is optional.  It is discovered the way `symengine` and
## `_ginac_ext` already are; when it is absent nothing changes and
## `blas_single_thread_available()` reports False, so the speedup can be had from
## the environment instead (OMP_NUM_THREADS=1 etc.).  Making it a hard dependency
## for a 1.72x is the maintainer's call, not this module's.
try:
    from threadpoolctl import threadpool_limits as _threadpool_limits
except ImportError:
    _threadpool_limits = None


def blas_single_thread_available():
    """True if the BLAS thread limit can be applied from inside the process."""
    return _threadpool_limits is not None


def _single_threaded_blas():
    """Limit BLAS to one thread for the duration, if that is possible here.

    NOTE ON REPRODUCIBILITY: a threaded LAPACK reduces in a different order, so
    results differ in the last bits between one and many threads.  Measured on the
    leapfrog: identical step count, `max|dv| = 2.16e-15 V`.  That is BLAS
    non-associativity, not a change of method -- but it does mean a run is only
    bit-reproducible against another run with the same thread count.
    """
    if _threadpool_limits is None:
        return contextlib.nullcontext()
    return _threadpool_limits(limits=1, user_api='blas')

class Transient(Analysis):
    """Simple transient analysis class.

    Time step is fixed.

    i(t) = c*dv/dt
    v(t) = L*di/dt

    The usual companion models are used.
    backward euler:
    i(n+1) = c/dt*(v(n+1) - v(n)) = geq*v(n+1) + Ieq
    v(n+1) = L/dt*(i(n+1) - i(n)) = req*i(n+1) + Veq

    def F(x): return i(x)+Geq(x)*x+u(x)+ueq(x0), G(x)+Geq(x)
    x0=x(n)
    x(n+1) = fsolve(F, x0, fprime=J)

    Linear circuit example:
    >>> circuit.default_toolkit = numeric
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> n2 = c.add_node('net2')
    >>> c['Is'] = IS(gnd, n1, i=10)    
    >>> c['R1'] = R(n1, gnd, r=1)
    >>> c['R2'] = R(n1, n2, r=1e3)
    >>> c['R3'] = R(n2, gnd, r=100e3)
    >>> c['C'] = C(n2, gnd, c=1e-5)
    >>> tran = Transient(c)
    >>> res = tran.solve(tend=10e-3,timestep=1e-4)
    >>> expected = 6.3
    >>> abs(res.v(n2, gnd)[-1] - expected) < 1e-2*expected #node 2 of last x
    True

    Linear circuit example:
    >>> from pycircuit.circuit.elements import ISin
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> c['Isin'] = ISin(gnd, n1, ia=1e-3, freq=16e3)
    >>> c['R'] = R(n1, gnd, r=200)
    >>> c['C'] = C(n1, gnd, c=1e-6)
    >>> c['L'] = L(n1, gnd, L=1e-4)
    >>> tran = Transient(c)
    >>> res = tran.solve(tend=260e-6,timestep=1e-6)
    >>> expected = 0.063
    >>> abs(res.v(n1,gnd)[-1]) < 1e-1*expected #node 2 of last x
    True
    
    """
    ## TODO:
    ## * Implement automatic timestep adjustment, using difference between
    ##   BE and trapezoidal as a measure of the error.
    ##   Reference: "Time Step Control in Transient Analysis", by SHUBHA VIJAYCHAND
    def _get_integrator(self):
        from pycircuit.circuit.integrator import Integrator, EulerIntegrator
        integrator = getattr(self.par, 'integrator', None)
        if integrator is None:
            return EulerIntegrator()
        if not isinstance(integrator, Integrator):
            raise TypeError(
                "integrator must be an Integrator instance (e.g. EulerIntegrator(), "
                "TrapezoidalIntegrator(), Gear2Integrator()), not %r" % (integrator,))
        return integrator
    

    parameters = Analysis.parameters + \
        [Parameter(name='analysis', desc='Analysis name', 
                   #default='transient'),
                   default='tran'),
         Parameter(name='reltol', 
                   desc='Relative tolerance', unit='', 
                   default=1e-4),
         Parameter(name='iabstol', 
                   desc='Absolute current error tolerance', unit='A', 
                   default=1e-12),
         ## NEWTON's x-tolerance on node rows, and nothing else.  Shared with `DC`,
         ## which uses the same 1e-12, so the operating point and the steps after it
         ## are solved to the same accuracy.
         ##
         ## This used to be 1e-6 and used for BOTH roles.  The 1e-12 -> 1e-6 change
         ## was reasoned about purely as a step-control knob (see `lte_vabstol`
         ## below) and its effect on Newton was never measured -- it loosened node
         ## convergence by 10^6 as a side effect, while `DC.vabstol` stayed at 1e-12,
         ## so every transient was seeded by an operating point solved a million
         ## times tighter than any step that followed it.  Decision 0.3a/0.3d in
         ## `doc/transient_work_plan.md` split the two roles; this is the Newton one.
         Parameter(name='vabstol',
                   desc='Absolute voltage error tolerance for the Newton solve',
                   unit='V',
                   default=1e-12),
         ## THE STEP CONTROLLER's tolerances, which are a different quantity: they
         ## apply to `lte = J^-1 Eg`, not to Newton's residual or its x-update.
         ##
         ## 1 uV is Spectre's `vabstol` default and SPICE's VNTOL.  The controller's
         ## value was 1e-12 V -- a million times tighter than either, and tighter
         ## than double precision can resolve against a 1 V signal.  On the
         ## 127-unknown leapfrog that alone collapsed the timestep to 5 ns against a
         ## 39 ns cap, because the controller accepts on max(|lte|/etol) over ALL
         ## unknowns and most of that circuit's nodes carry no signal, so etol
         ## degenerated to TRTOL*abstol on numerical noise.  Relaxing it to Spectre's
         ## value cut the step count 5.4x with the waveform unchanged to 0.5%.  That
         ## result belongs to THIS parameter, not to `vabstol`.
         Parameter(name='lte_vabstol',
                   desc='Absolute voltage tolerance for the local truncation error',
                   unit='V',
                   default=1e-6),
         Parameter(name='lte_iabstol',
                   desc='Absolute current tolerance for the local truncation error',
                   unit='A',
                   default=1e-12),
         ## What the RELATIVE part of the LTE tolerance is measured against --
         ## Spectre's parameter of the same name, and Spectre's default.
         ##
         ## `pointlocal` is what pycircuit did for its whole history: each unknown
         ## referenced to itself, at this instant.  On a node carrying no signal
         ## that reference tends to zero, so the tolerance collapses to the absolute
         ## floor and the controller chases numerical noise on an idle node -- which
         ## is the defect `lte_vabstol` was raised a millionfold to work around.
         ## Measured on the leapfrog: at `lte_vabstol = 1e-12`, `pointlocal` needs
         ## 3.53x the steps of the shipped configuration where `sigglobal` needs
         ## 1.49x, i.e. it removes 81% of the excess.
         ##
         ## DECISION D3 ADOPTED `sigglobal` AND GATE D3-a SENT IT BACK.  Under
         ## `sigglobal` the tolerance is referenced to the largest signal, so steps
         ## grow larger -- and on TRAPEZOIDAL that amplifies the `(-1)^n` mode 4g(b)
         ## has not yet removed, to the point where accuracy stops responding
         ## monotonically to `reltol`:
         ##
         ##   trap-classic  err 4.89e-4  2.23e-4  3.14e-4  7.14e-5   <- rises
         ##   trap-ywr      err 5.94e-4  3.82e-5  2.65e-4  3.18e-5   <- rises
         ##
         ## Euler and both Gear2 variants are monotone under `sigglobal` and need
         ## 1.7-2.5x FEWER steps, so the mode is right and only the timing is wrong.
         ## Default stays `pointlocal` until 4g(b) makes the trapezoidal estimator
         ## mode-free; flipping it is then a one-line change with the gate already
         ## written.  See the D3 gates in `doc/transient_work_plan.md`.
         Parameter(name='relref',
                   desc="Reference for the relative LTE tolerance: 'pointlocal' "
                        "(each unknown against itself, the default for now), "
                        "'alllocal' (against its own past maximum), or 'sigglobal' "
                        "(against the largest signal anywhere -- Spectre's default, "
                        "and pycircuit's once 4g(b) lands)",
                   unit='',
                   default='pointlocal'),
         Parameter(name='maxiter', 
                   desc='Maximum number of iterations', unit='', 
                   default=100),
         Parameter(name='integrator',
                   desc='Integration strategy (an Integrator instance, e.g. '
                        "EulerIntegrator(), TrapezoidalIntegrator(), "
                        "Gear2Integrator(lte_formula='ywr')); default EulerIntegrator()",
                   unit='',
                   default=None),
         Parameter(name='uic',
                   desc='Use initial conditions (skip DC OP computation)', unit='',
                   default=False),
         Parameter(name='minbreak',
                   desc='Minimum time difference for breakpoint events', unit='s',
                   default=1e-14),
         Parameter(name='bypass',
                   desc='Enable device model bypassing', unit='',
                   default=False),
         Parameter(name='minstep',
                   desc='Minimum timestep to prevent infinite loops', unit='s',
                   default=1e-18),
         ## STAGE 3.  The opening step, which the controller must accept
         ## unevaluated because there is no history to difference against.  `None`
         ## means `timestep * 1e-3`; see `Transient._opening_step` for why opening
         ## at `timestep` made `reltol` unable to influence the answer at all.
         Parameter(name='firststep',
                   desc='Size of the first timestep; None means timestep*1e-3. '
                        'The first step cannot be error-checked, so taking it at '
                        'max_step lets its error dominate the whole run.',
                   unit='s',
                   default=None),
         Parameter(name='bypasstol',
                   desc='Bypass tolerance for device models', unit='V',
                   default=None)]

    def __init__(self, cir, toolkit=None, irefnode=None, **kvargs):
        self.parameters = super(Transient, self).parameters + self.parameters
        ## `toolkit` was accepted and then DROPPED -- it was never forwarded, so
        ## `Transient(cir, toolkit=X)` silently ran on `cir.toolkit` instead.  It
        ## went unnoticed because callers pass the toolkit the circuit already has,
        ## which makes the two agree by coincidence rather than by construction.
        super(Transient, self).__init__(cir, toolkit=toolkit, **kvargs)
    
        self._qlast  = None #q history
        self._iqlast = None #dq/dt history
        
        self._dt = None
        self._dt_last = None
        self._is_first_step = True
        ## Distinct from _is_first_step, which is re-armed at every breakpoint to
        ## force an order drop.  This one is true only until the first step of a
        ## run has been accepted, and it is what the step controller is given.
        self._no_history = True
    ## import it from there instead.
    ## But it's an object method requiring a DC as self
    ## so using DC._newton doesn't work
    def _newton(self, func, x0): 
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

        from pycircuit.circuit.nrsolver import NoConvergenceError
        solver = self._get_nrsolver()
        scaler = self._get_scaler()
        try:
            x_res, _ = solver.solve_system(
                x0,
                refnode_removed(func, self.irefnode, self.toolkit),
                self.toolkit,
                self.par.reltol,
                abstol,
                xtol,
                self.par.maxiter,
                limiter=limiter_func,
                scaler=scaler
            )
        ## NARROW, deliberately.  This used to be `except Exception`, which turned
        ## every failure inside a device model into "the circuit did not converge" --
        ## a `ZeroDivisionError` from a source with `tr=0`, a `TypeError` from a
        ## mis-specified parameter, an `AttributeError` from a typo in a subclass.
        ## All of them were reported as a convergence failure, which sends the reader
        ## to look at the bias point of a circuit whose real problem is a bug three
        ## frames down.  The solvers already classify what they mean, so only their
        ## own exceptions and genuine linear-algebra failures are translated here;
        ## everything else propagates with its original type and traceback.
        except SingularMatrix:
            raise
        except NoConvergenceError as e:
            ## The solvers wrap a singular factorisation as NoConvergenceError
            ## ("Singular Jacobian: ..."); promote that to SingularMatrix so callers
            ## can tell "no solution here" from "could not get there".  Matching on
            ## the message is weak and stage 6 replaces it with a real classification
            ## off the zero pivot -- but it is what the previous code did, and
            ## changing the taxonomy is not this stage's job.
            if 'Singular' in str(e) or 'linalgerror' in str(e).lower():
                raise SingularMatrix(str(e)) from e
            raise
        except LinAlgError as e:
            raise SingularMatrix(str(e)) from e
        
        x = x_res
        
        # Insert reference node voltage
        return self.toolkit.concatenate((x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))
    

    
    def get_diff(self, q, C, method=None):
        """Method used to calculate time derivative for charge storing elements (i_eq and g_eq)."""
        # Determine the active integrator based on step size variations
        h_last = getattr(self, '_dt_last', self._dt)
        self.active_integrator = self.base_integrator.check_order_drop(
            self._dt, h_last, getattr(self, '_is_first_step', False)
        )
        
        iq, geq = self.active_integrator.compute_derivatives(
            q_curr=q, C_curr=C, h_curr=self._dt, 
            q_last=self._qlast, iq_last=self._iqlast, h_last=h_last,
            is_first_step=getattr(self, '_is_first_step', False),
            toolkit=self.toolkit
        )
        
        self._iq = iq
        ## The class actually running this step -- check_order_drop() may have
        ## dropped to a lower-order integrator than the one requested, so this
        ## is derived from the live object rather than compared against a
        ## removed user-supplied method name.
        self._effective_method = type(self.active_integrator).__name__
        return iq, geq

    
    
    def _opening_step(self, timestep):
        """The size of the first step of a run.

        STAGE 3.  A run used to open at `timestep`, which is also `max_step` -- the
        largest step the controller is ever allowed to take.  The step controller
        accepts the first step unevaluated, because with no history there is nothing
        to difference and no truncation error can be estimated, so that opening step
        was both the **largest** and the **only unchecked** step in the run.

        Its error then dominated everything after it.  Measured on an RC step
        response against the analytic solution, trapezoidal, before this method
        existed: the global error was **1.3212e-01 at reltol 1e-3, 1e-4, 1e-5 AND
        1e-6** -- identical to five digits -- while the step count went from 24 to
        195.  Eight times the work for the same answer.  Backward Euler and the
        two second-order methods also agreed to five digits, which is why the
        integrator choice appeared not to matter.

        Opening at `timestep * 1e-3` costs one cheap step and leaves the controller
        to grow the step from there, which it does geometrically, so the ramp is
        paid off within a handful of steps.

        The principled alternative is a Hairer-style estimate from `q'`/`q''` at the
        operating point.  The plan asks for the ramp first and a *measurement* of
        whether the estimate is worth the complexity, rather than an assumption --
        see the outcome recorded under gate 3-1.
        """
        firststep = getattr(self.par, 'firststep', None)
        if firststep is None:
            return timestep * 1e-3
        if firststep <= 0:
            raise ValueError(
                "firststep must be positive, not %r; pass None to use the default "
                "ramp of timestep*1e-3" % (firststep,))
        ## An opening step larger than max_step would be capped on the very next
        ## step anyway, and asking for one is more likely a mistake than an intent.
        return min(firststep, timestep)

    def _q_at(self, x):
        """``cir.q(x)``, reusing the value computed during the last assembly.

        STAGE 2c.  This is a *memoisation*, not an approximation: the cached value
        was produced by the same function at the same state, so it is bit-identical
        to recomputing, and the whole of stage 2 is defined as behaviour-preserving.
        The guard is deliberately strict -- identity first, then full equality --
        because serving a charge vector from the wrong state would corrupt the LTE
        estimate silently, which is precisely the failure class stage 1 removed.
        """
        cached = getattr(self, '_q_cache', None)
        if cached is not None:
            x_cached, q_cached = cached
            if x_cached is x:
                return q_cached
            if (x_cached is not None and x is not None
                    and getattr(x_cached, 'shape', None) == getattr(x, 'shape', None)
                    and bool(self.toolkit.alltrue(x_cached == x))):
                return q_cached
        return self.cir.q(x, self.epar)

    def _solve_operating_point(self, refnode):
        """Solve the DC operating point that seeds the transient.

        A failure here **raises**.  It used to substitute a vector of zeros, so a
        circuit that had no operating point at all -- or one whose solve hit a bug in
        a device model -- returned a complete, plausible-looking waveform computed
        from a bias point that was never found.  Nothing in the result distinguished
        that from a successful run, which makes it the most expensive class of defect
        this module had: it does not fail, it lies.

        The inner `DC` is constructed from *this* analysis's configuration rather
        than from `DC`'s defaults.  Before, `DC(self.cir)` inherited none of the
        transient's toolkit, environment parameters, tolerances, solver or scaler, so
        the operating point could be solved at a different temperature, to a
        different accuracy, and with a different Newton strategy than every step that
        followed it -- and the mismatch was invisible.
        """
        from pycircuit.circuit.dcanalysis import DC

        ## Only forward what DC actually declares, so a parameter that exists on
        ## Transient but not on DC (e.g. `integrator`) does not raise a KeyError,
        ## and so this keeps working if either parameter list changes.
        dc_par_names = {p.name for p in DC.parameters}
        shared = {}
        for name in ('reltol', 'iabstol', 'vabstol', 'maxiter',
                     'bypass', 'bypasstol', 'epar', 'nrsolver', 'scaler'):
            if name not in dc_par_names:
                continue
            try:
                value = getattr(self.par, name)
            except (AttributeError, KeyError):
                continue
            if value is not None:
                shared[name] = value

        dc = DC(self.cir, toolkit=self.toolkit, refnode=refnode, **shared)
        try:
            return dc.solve().x
        except (NoConvergenceError, SingularMatrix) as exc:
            raise NoConvergenceError(
                "Transient could not find a DC operating point to start from: %s\n"
                "The transient has NOT been run.  Either fix the bias condition, or "
                "start deliberately from a known state:\n"
                "  * uic=True            -- start from zeros (SPICE's 'use initial "
                "conditions'), or\n"
                "  * x0=<vector>         -- start from an operating point you supply.\n"
                "Both are explicit choices; substituting zeros silently is what this "
                "error replaced." % (exc,)) from exc

    def solve_timestep(self, x0, t, refnode=gnd, provided_function=None):
        n=self.cir.n
        dt = self._dt
        
        def func(x):
            ## `self.epar`, not the module-level `defaultepar`.  Omitting it meant
            ## every device in a transient was evaluated at defaultepar's T = 300 K
            ## whatever the caller asked for, and -- because `Analysis.__init__`
            ## attaches `bypasstol` to the analysis's own epar and nowhere else --
            ## every device took its `except AttributeError` branch and the `bypass`
            ## parameter did nothing at all.
            C = self.cir.C(x, self.epar)
            q = self.cir.q(x, self.epar)
            ## STAGE 2c.  Stash the charge vector alongside the state it belongs
            ## to.  `solve()` needs `q` at the converged point twice more -- once
            ## for the step controller and once for the history roll -- and was
            ## recomputing the whole assembly both times at an x it had already
            ## evaluated.  Measured 5.08 `q` assemblies per accepted step against
            ## 3.06 for every other stamp; the difference is exactly those two.
            ##
            ## Keyed by the state so a stale value can never be served: the check
            ## below is identity-then-equality on x, not a bare "did we cache".
            self._q_cache = (x, q)
            iq, Geq = self.get_diff(q, C)
            f = (self.cir.i(x, self.epar) + iq
                 + self.cir.u(t, self.epar, analysis=self.par.analysis))
            J = self.cir.G(x, self.epar) + Geq
            return self.toolkit.array(f, dtype=float), self.toolkit.array(J, dtype=float)

        def jacobian_only(x):
            """The converged-point evaluation, without the residual nobody reads.

            ITEM 2+.2.  Newton needs `f` on every iteration -- it is what drives the
            update.  The *final* evaluation at the converged point is different: its
            `f` is unpacked by both callers and then never referenced again (`solve`
            at the `x, feval, J, f = ...` site, and `_solve_coupled` likewise). Only
            `J` is consumed, by the step controller.

            So `cir.i(x)` and `cir.u(t)` were assembled once per accepted step and
            discarded.  `C`, `q` and `G` are all still needed -- `C` and `G` build
            `J`, and `q` feeds the charge cache and the history roll -- so this is
            not a cheaper approximation of the same work, it is the same work minus
            two vectors that had no consumer.

            `f` is returned as None rather than a zero vector, so any future caller
            that starts reading it fails loudly instead of silently using zeros --
            which is the whole lesson of stage 1.
            """
            C = self.cir.C(x, self.epar)
            q = self.cir.q(x, self.epar)
            self._q_cache = (x, q)
            iq, Geq = self.get_diff(q, C)
            J = self.cir.G(x, self.epar) + Geq
            return None, self.toolkit.array(J, dtype=float)

        x=self._newton(func,x0)
        ## `provided_function` is the one caller that consumes `f`, so it gets the
        ## full evaluation; everything else takes the reduced one.
        if provided_function is not None:
            f, J = func(x)
        else:
            f, J = jacobian_only(x)

        if provided_function is not None:
            result=x,provided_function(f,J,self.cir.C(x, self.epar)), J, f
        else:
            result=x,None, J, f
        return result
    
    def solve(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None, fixed_timestep=False, coupled_lte=False, analytical_eh=True):
        ## Stage 2a: hold BLAS to one thread for the whole run.  It wraps the whole
        ## transient rather than just the linear solve because the win is not in the
        ## solve -- that is ~2% of runtime, so even an infinite speedup there could
        ## not produce the measured 1.72x.  The cost is thread-pool overhead spread
        ## across the many small numpy operations in assembly, and that is only
        ## avoided by setting the limit once, outside the loop.
        with _single_threaded_blas():
            return self._solve(refnode, tend, x0, timestep, provided_function,
                               fixed_timestep, coupled_lte, analytical_eh)

    def _solve(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None, fixed_timestep=False, coupled_lte=False, analytical_eh=True):
        if coupled_lte:
            return self._solve_coupled(refnode, tend, x0, timestep, provided_function, analytical_eh)

        ## Respect a step controller injected by the caller (e.g. PIController);
        ## only fall back to the default IntegralController when none was set.
        if getattr(self, 'step_controller', None) is None:
            from pycircuit.circuit.stepcontroller import IntegralController
            self.step_controller = IntegralController()
        ## ITEM 2+.3.  Applied to whichever controller is in use, including one
        ## the caller injected, and re-applied every run so the running maximum
        ## a global mode keeps cannot leak from a previous solve.
        self.step_controller.set_relref(self.par.relref)

        X = []
        self.irefnode=self.cir.get_node_index(refnode)
        n = self.cir.n
        if x0 is None:
            if self.par.uic:
                # Use Initial Conditions = skip DC operating point, start at zero
                x0 = self.toolkit.zeros(n)
            else:
                x0 = self._solve_operating_point(refnode)
            x = x0
        else:
            x = x0
        
        self.base_integrator = self._get_integrator()
        hist_len = max(2, self.base_integrator.get_required_history())
        q0 = self.cir.q(x, self.epar)
        self._qlast = self.toolkit.array([q0 for _ in range(hist_len)])
        self._iqlast = self.toolkit.zeros((hist_len, n))
        
        X.append(copy(x))
        if hasattr(self.cir, 'accept_step'):
            self.cir.accept_step(0.0, X[-1], self.epar)
        
        timelist = []
        self._is_first_step = True
        self._no_history = True
        t = 0.0
        max_step = timestep
        ## The opening ramp exists to stop the ONE step the controller cannot check
        ## from dominating the run.  Under `fixed_timestep` there is no controller
        ## and `dt` is never updated, so ramping would not open small and grow -- it
        ## would run the ENTIRE simulation at `timestep*1e-3`, a thousand times more
        ## steps for a result the caller explicitly asked to be uniform.  Caught by
        ## `test_transient_RLC` and three others, which is what the suite gate is for.
        dt = timestep if fixed_timestep else self._opening_step(timestep)
        TRTOL = 7.0
        ## Bound the number of consecutive LTE rejections at a single time point.
        ## Near a source discontinuity (e.g. a VPulse corner) the truncation-error
        ## estimate can stay above tolerance for arbitrarily small steps while the
        ## stored history is frozen; without a cap the step size collapses and the
        ## solve grinds indefinitely.  After MAX_REJECT rejections we accept the
        ## already-converged Newton solution (only its LTE was too high) so time
        ## advances and the integrator history refreshes.
        reject_count = 0
        MAX_REJECT = 3
        ## Set by the force-accept path below and consumed at the top of the next
        ## iteration, exactly like `was_break_step`.  Both mean the same thing to
        ## the integrator -- "do not trust a 2nd-order polynomial through this
        ## point" -- and they arrive from opposite ends: a breakpoint knows the
        ## history is about to be discontinuous, a force-accept has just found out.
        force_order_drop = False

        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        ## SOLUTION-flavoured, not residual-flavoured.  This vector is used by the
        ## step controller as a tolerance on `lte = J^-1 * Eg`, which carries the
        ## units of the solution vector x -- volts on node rows, amps on branch rows.
        ## `_newton` needs the other flavour, because there the tolerance applies to
        ## the residual f (KCL currents at nodes), and it builds both separately as
        ## `abstol`/`xtol`.  This is the `xtol` one; using `_newton`'s `abstol` here
        ## silently applied iabstol (1 pA) as a *voltage* tolerance to every node,
        ## which is what made a larger `vabstol` have no effect on node rows at all.
        ##
        ## It reads `lte_vabstol`/`lte_iabstol`, NOT `vabstol`/`iabstol`.  Sharing
        ## the Newton parameters was decision 0.3a's defect: one knob could not be
        ## moved for the controller without silently moving Newton's convergence
        ## criterion with it.
        abstol = self.toolkit.concatenate((self.par.lte_vabstol * ones_nodes,
                                          self.par.lte_iabstol * ones_branches))

        was_break_step = False
        while t < tend:
            # --- BREAKPOINT HANDLING ---
            # A breakpoint generally signifies a mathematical discontinuity (such as 
            # the sharp corner of a VPulse/IPulse square wave).
            # 
            # When the solver approaches a breakpoint, it does three vital things to 
            # maintain mathematical stability:
            # 1. Truncates Step (dt): If normal step size overshoots the breakpoint, 
            #    it forces dt to land *exactly* on the breakpoint timestamp.
            # 2. Flags the Breakpoint: was_break_step is set to True.
            # 3. Drops the Integration Order: Immediately after crossing the
            #    breakpoint (was_break_step == True in the next iteration), it
            #    sets `self._is_first_step = True`.
            #
            # Why? Integrators (like Gear2 or Trapezoidal) use past state history
            # to fit a smooth mathematical polynomial. If they
            # tried to fit a polynomial across a sharp discontinuous edge, the
            # simulation would suffer from massive artificial ringing and overshoot.
            # So the method drops to a safer 1st-order one
            # (like Backward Euler) to gracefully navigate the corner and rebuild.
            #
            # `_is_first_step` does NOT mean "there is no history": the q and iq
            # ring buffers keep rolling across a breakpoint, and the step just
            # taken is perfectly good data.  It means "do not trust a 2nd-order
            # polynomial through this point".  Those are different claims, and
            # conflating them is what made `max_step` a correctness knob rather
            # than a cost knob: the step controller used this same flag to skip
            # the error check entirely, and `Sin.next_event` fires every quarter
            # period, so a VSin drive produced a periodic, drive-synchronous,
            # full-`max_step` step that was never checked at all.  The controller
            # is therefore handed `_no_history` instead, which is true only at
            # the genuine start of a run -- where the LTE really cannot be
            # estimated and accepting is the only option.
            if was_break_step or force_order_drop:
                self._is_first_step = True
            force_order_drop = False

            next_t_break = self.cir.next_event(t)
            
            # Ensure next_t_break strictly advances time to avoid infinite dt=0 loops
            if next_t_break <= t + self.par.minbreak * max(abs(t), 1.0):
                next_t_break = self.cir.next_event(t + (self.par.minbreak * 1e3) * max(abs(t), 1.0))
            
            if t + dt > next_t_break:
                dt = float(next_t_break - t)
                was_break_step = True
            else:
                was_break_step = False
                
            if t + dt > tend:
                dt = tend - t
            
            self._dt = dt
            next_t = t + dt
            
            try:
                x, feval, J, f = self.solve_timestep(X[-1], next_t, provided_function=provided_function)
            except NoConvergenceError:
                dt = dt * 0.25
                if dt < getattr(self.par, 'minstep', 1e-18):
                    raise RuntimeError(f"Transient solver failed to converge: timestep shrank below {getattr(self.par, 'minstep', 1e-18):g}s at t={t}")
                continue
                
            if not fixed_timestep:
                accept, dt_next = self.step_controller.evaluate_step(
                    x_curr=x,
                    x_last=X[-1],
                    q_curr=self._q_at(x),
                    q_last_hist=self._qlast,
                    iq_last_hist=self._iqlast,
                    h_curr=dt,
                    h_last=getattr(self, '_dt_last', dt),
                    no_history=self._no_history,
                    J=J,
                    active_integrator=self.active_integrator,
                    irefnode=self.irefnode,
                    reltol=self.par.reltol,
                    abstol=abstol,
                    toolkit=self.toolkit,
                    max_step=max_step,
                    TRTOL=TRTOL,
                    ## `relref`'s global modes must not mix volts with amps.
                    n_nodes=len(self.cir.nodes),
                )
                
                if not accept and reject_count < MAX_REJECT:
                    reject_count += 1
                    dt = dt_next
                    if dt < getattr(self.par, 'minstep', 1e-18):
                        raise RuntimeError(f"Transient solver integration error: timestep shrank below {getattr(self.par, 'minstep', 1e-18):g}s at t={t}")
                    continue
                elif not accept:
                    ## STAGE 4b -- THE ESCAPE HATCH USED TO GROW 10x, WHICH IS THE
                    ## WRONG SIGN AND OUTSIDE BDF-2'S STABILITY BOUND.
                    ##
                    ## Reaching here means the LTE estimate stayed over tolerance
                    ## for MAX_REJECT successively smaller steps.  The old response
                    ## was `next_dt = min(max_step, dt * 10.0)`: grow tenfold in
                    ## answer to an error that was already too large.  Variable-step
                    ## BDF-2 is zero-stable only below `ZERO_STABILITY_RATIO`
                    ## (2.414214); at 10x the parasitic root is 4.76, so the step
                    ## that follows a force-accept amplified the previous solution
                    ## instead of forgetting it -- and nothing warned.
                    ##
                    ## Measured before the change (stiff RLC, reltol 1e-5,
                    ## `Trapezoidal('ywr')`): 78 force-accepts in 873 accepted steps,
                    ## and every one of the 9 accepted-step ratios above 2.414
                    ## across the whole run sat immediately after one, the largest
                    ## being exactly 10.0.  **The shipped default is not exempt**:
                    ## the same circuit at reltol 1e-3 under `Gear2('ywr')` reached
                    ## here once and took the 10x once.  That case was found by
                    ## this warning, after a sweep of three tighter tolerances had
                    ## concluded the default no longer reached the path at all.
                    ##
                    ## What a stalled high-order estimate is actually asking for is
                    ## a LOWER ORDER, not a bigger step: order 1 differences one
                    ## past point instead of two, so it is far less sensitive to the
                    ## stale history that a discontinuity leaves behind, and it
                    ## still gets a real error estimate -- the controller is handed
                    ## `_no_history`, not `_is_first_step`, so an order-dropped step
                    ## is error-controlled rather than accepted blind.  Growth is
                    ## then bounded by the same clamp every other accepted step
                    ## obeys, which is what makes "no accepted ratio exceeds 2.414"
                    ## true of the run as a whole rather than of its quiet parts.
                    force_order_drop = True
                    next_dt = min(max_step, dt * MAX_GROWTH_RATIO)
                    ## An unbounded accepted truncation error must not be invisible.
                    ## This is the same failure class stage 1 exists to remove: the
                    ## result is still returned, and without this the caller has no
                    ## way to learn that part of it was not error-controlled.
                    warnings.warn(
                        'transient: local truncation error still above tolerance '
                        'after %d rejections at t=%.6g s; accepting the step at '
                        'h=%.6g s with an order drop. The accepted error is '
                        'unbounded -- treat the waveform near this time with '
                        'suspicion.' % (MAX_REJECT, t, dt),
                        ## 3, not 2: the loop lives in `_solve`, which `solve`
                        ## calls, so 2 attributes the warning to `solve`'s own body
                        ## and tells the caller nothing about which simulation
                        ## produced it.  Verified by reading the reported filename.
                        RuntimeWarning, stacklevel=3)
                else:
                    next_dt = dt_next
                reject_count = 0

            t = next_t
            timelist.append(t)
            X.append(copy(x))
            
            if hasattr(self.cir, 'accept_step'):
                self.cir.accept_step(t, X[-1], self.epar)
            
            # --- INTEGRATOR HISTORY RING BUFFERS ---
            # To support 2nd-order (and higher) integration methods, we must preserve the 
            # charge (q) and current/derivative (iq) of previous timesteps.
            # We push the newest values to index 0, and slice off the oldest `[:-1]` to 
            # maintain a constant buffer size (e.g. size 2 for Gear2).
            # This acts as a mathematical sliding window across the simulation time.
            self._iqlast = self.toolkit.concatenate((self.toolkit.array([self._iq]), self._iqlast))[:-1]
            self._qlast = self.toolkit.concatenate((self.toolkit.array([self._q_at(x)]), self._qlast))[:-1]
            self._dt_last = dt
            
            self._is_first_step = False
            self._no_history = False
            
            if not fixed_timestep:
                dt = next_dt
            
        X = self.toolkit.array(X[1:]).T
        timelist = self.toolkit.array(timelist)
        
        self.result = CircuitResult(self.cir, x=X, xdot=None,
                                    sweep_values=timelist, 
                                    sweep_label='time', 
                                    sweep_unit='s')
        
        return self.result


    def _solve_coupled(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None, analytical_eh=True):
        import numpy as np
        from copy import copy
        from pycircuit.circuit.analysis import CircuitResult
        
        X = []
        self.irefnode = self.cir.get_node_index(refnode)
        n = self.cir.n
        if x0 is None:
            ## Same fix as `solve()`: a failed operating point raises rather than
            ## silently becoming a vector of zeros.  This path had the defect too.
            ##
            ## And it must honour `uic`, which it previously ignored -- one of the
            ## four inputs 0.1d found this path silently dropping.  That went
            ## unnoticed because ignoring `uic` and failing the DC produced the same
            ## zeros the caller wanted; with the silent fallback gone, `uic=True`
            ## would have raised on a circuit that has no operating point by design.
            if self.par.uic:
                x0 = self.toolkit.zeros(n)
            else:
                x0 = self._solve_operating_point(refnode)

        x = x0
        self.base_integrator = self._get_integrator()
        hist_len = max(2, self.base_integrator.get_required_history())
        q0 = self.cir.q(x, self.epar)
        self._qlast = self.toolkit.array([q0 for _ in range(hist_len)])
        self._iqlast = self.toolkit.zeros((hist_len, n))
        
        X.append(copy(x))
        if hasattr(self.cir, 'accept_step'):
            self.cir.accept_step(0.0, X[-1], self.epar)
        timelist = []
        
        self._is_first_step = True
        self._no_history = True
        t = 0.0
        ## Same opening ramp as `solve()` -- this path had the same defect.
        h = self._opening_step(timestep)
        max_step = timestep
        TRTOL = 7.0
        minstep = getattr(self.par, 'minstep', 1e-18)

        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        ## Solution-flavoured, for the same reason as in `solve` above: the coupled
        ## controller also applies this to `lte`, not to the residual -- so it takes
        ## the `lte_*` tolerances, not Newton's.
        abstol = self.toolkit.concatenate((self.par.lte_vabstol * ones_nodes,
                                          self.par.lte_iabstol * ones_branches))
        reltol = self.par.reltol

        ## Coupled adaptive time-stepping, Fang, "A New Time-Stepping Method for
        ## Circuit Simulation" (DAC 2013).  The circuit solution and the step size
        ## are co-determined at each time point: converge the circuit at a fixed
        ## step (stage 1), then bring the local truncation error into the accept
        ## band by driving the step with Gear's formula and re-solving.  This is
        ## the paper's robust "approximate Newton" (sec. 3.4); it replaces the
        ## exact (N+1) Schur update, which is very sensitive to step changes and
        ## collapses the step size.  The LTE evaluation and Gear prediction are
        ## shared with the standard adaptive controller (IntegralController) so
        ## the coupled and adaptive paths stay consistent.
        from pycircuit.circuit.stepcontroller import IntegralController
        from pycircuit.circuit.nrsolver import NoConvergenceError
        controller = IntegralController()
        MAX_LTE_ITERS = 10

        while t < tend:
            if t + h > tend:
                h = tend - t

            ## Co-determine (x, h): converge the circuit at h_curr, evaluate the
            ## LTE, and while it is above the band shrink the step (Gear-predicted
            ## by the controller) and re-solve.
            h_curr = h
            x_curr = copy(X[-1])
            h_next = h
            ## Whether any solve at this time point actually succeeded.  Without
            ## this the loop could exhaust its retries and fall through to the
            ## accept block below, which advanced `t` by the collapsed `h_curr` and
            ## appended the PREVIOUS solution as if it were a new one -- while `h`
            ## was restored to full size for the next outer iteration.  The result
            ## was a livelock: 10 Newton solves bought `h*0.25^10` of simulated time,
            ## measured at 9.5367e-13 s per iteration against a predicted 9.5367e-13,
            ## i.e. ~2.1e7 further Newton attempts to finish a 5 us run.  It neither
            ## raised nor returned.  (On a FIRST-step failure it instead died with
            ## `AttributeError: _iq`, because `_iq` is set inside `solve_timestep`.)
            ## See `benchmarks/transient_review/stage0_1d_coupled_livelock.py`.
            converged = False
            for lte_iter in range(MAX_LTE_ITERS):
                self._dt = h_curr
                try:
                    x_new, feval, J, f = self.solve_timestep(
                        X[-1], t + h_curr, provided_function=provided_function)
                except NoConvergenceError:
                    ## Circuit did not converge at this step -> shrink and retry.
                    h_curr *= 0.25
                    if h_curr < minstep:
                        raise RuntimeError(f"Coupled transient: timestep shrank below {minstep:g}s at t={t}")
                    continue
                converged = True

                accept, h_next = controller.evaluate_step(
                    x_curr=x_new, x_last=X[-1],
                    q_curr=self._q_at(x_new),
                    q_last_hist=self._qlast, iq_last_hist=self._iqlast,
                    h_curr=h_curr, h_last=getattr(self, '_dt_last', h_curr),
                    no_history=self._no_history, J=J,
                    active_integrator=self.active_integrator,
                    irefnode=self.irefnode, reltol=reltol, abstol=abstol,
                    toolkit=self.toolkit, max_step=max_step, TRTOL=TRTOL,
                    n_nodes=len(self.cir.nodes))

                x_curr = x_new
                if accept or lte_iter == MAX_LTE_ITERS - 1:
                    break
                if h_next < minstep:
                    raise RuntimeError(f"Coupled transient: timestep shrank below {minstep:g}s at t={t}")
                h_curr = h_next

            if not converged:
                raise NoConvergenceError(
                    "Coupled transient: Newton failed to converge at t=%g after %d "
                    "step reductions (h fell from %g to %g s). The run is "
                    "abandoned rather than advanced -- continuing here would repeat "
                    "the same %d solves per outer iteration while advancing time by "
                    "h*0.25^%d, which neither converges nor terminates."
                    % (t, MAX_LTE_ITERS, h, h_curr, MAX_LTE_ITERS, MAX_LTE_ITERS))

            t += h_curr
            timelist.append(t)
            X.append(copy(x_curr))

            if hasattr(self.cir, 'accept_step'):
                self.cir.accept_step(t, X[-1], self.epar)

            self._dt = h_curr
            self._dt_last = h_curr
            self._is_first_step = False
            self._no_history = False
            self._iqlast = self.toolkit.concatenate((self.toolkit.array([self._iq]), self._iqlast))[:-1]
            self._qlast = self.toolkit.concatenate((self.toolkit.array([self._q_at(x_curr)]), self._qlast))[:-1]

            ## Next step: Gear-predicted size (already bounded by max_step).
            h = min(max_step, max(h_next, minstep))
            
        X = self.toolkit.array(X[1:]).T
        timelist = self.toolkit.array(timelist)
        self.result = CircuitResult(self.cir, x=X, xdot=None,
                                    sweep_values=timelist, 
                                    sweep_label='time', 
                                    sweep_unit='s')
        return self.result



if __name__ == "__main__":
    import doctest
    doctest.testmod()
