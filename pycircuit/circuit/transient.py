# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import contextlib
import time
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

def resample_uniform(t, x, npoints=None, step=None):
    """Interpolate a transient result onto a UNIFORM time grid.

    STAGE 10.2.  A transient returns the solver's own adaptive points, so their
    spacing is whatever the step controller chose -- measured at a **1000x**
    spread on an ordinary RC driven by a sine.  Anything that needs a uniform
    grid, above all an FFT, therefore has to resample, and every caller has been
    left to do that themselves: `benchmarks/nonlinear_leapfrog_sweep.py` reads

        spec = np.fft.rfft(np.interp(grid, t, v)) / npt

    -- `np.interp` is LINEAR, applied to a solution the integrator computed to
    second order.  That throws away an order of accuracy before the transform,
    which is a real hazard the caller was forced to own rather than a stylistic
    preference.

    QUADRATIC, to match the integrator.  Three-point Lagrange through the
    interval's neighbours, the same reasoning `TLine._interpolate_history` gives
    for its own history lookup: a first-order interpolant feeding a second-order
    method injects error the method never made.  Measured on the solver's real
    (non-uniform) grid, interpolating a known signal so only the interpolation
    error is present:

        signal                  linear      quadratic    ratio
        fundamental           1.12e-03      8.11e-05     13.8x
        5th harmonic          1.55e-02      1.20e-02      1.3x
        decaying exponential  2.00e-04      4.90e-06     40.8x

    **The 5th-harmonic row is the honest one**: where the adaptive grid is barely
    resolving the signal, no interpolant recovers what was not sampled, and the
    right fix there is a smaller `max_step`, not a cleverer resample.

    Give exactly one of ``npoints`` or ``step``.
    """
    import numpy

    t = numpy.asarray(t, dtype=float)
    x = numpy.asarray(x)
    if (npoints is None) == (step is None):
        raise ValueError('give exactly one of npoints or step')
    if t.size < 3:
        raise ValueError('need at least 3 points to interpolate quadratically, '
                         'got %d' % t.size)

    if npoints is not None:
        grid = numpy.linspace(t[0], t[-1], int(npoints))
    else:
        ## `t[-1] + step/2` so a grid that divides the interval exactly still
        ## includes the endpoint, without inventing a point beyond it.
        grid = numpy.arange(t[0], t[-1] + float(step) / 2.0, float(step))
        grid = grid[grid <= t[-1] * (1 + 1e-12)]

    ## The interval containing each output point, clamped so the three-point
    ## stencil stays inside the data at both ends.
    idx = numpy.clip(numpy.searchsorted(t, grid) - 1, 1, t.size - 2)
    t0, t1, t2 = t[idx - 1], t[idx], t[idx + 1]
    L0 = (grid - t1) * (grid - t2) / ((t0 - t1) * (t0 - t2))
    L1 = (grid - t0) * (grid - t2) / ((t1 - t0) * (t1 - t2))
    L2 = (grid - t0) * (grid - t1) / ((t2 - t0) * (t2 - t1))

    if x.ndim == 1:
        out = x[idx - 1] * L0 + x[idx] * L1 + x[idx + 1] * L2
    else:
        ## Rows are unknowns, columns are time -- the shape a CircuitResult holds.
        out = (x[:, idx - 1] * L0 + x[:, idx] * L1 + x[:, idx + 1] * L2)
    return grid, out


class TransientStatistics(object):
    """What a transient run actually did, as opposed to what it returned.

    STAGE 6(c).  Every number here was already being computed and thrown away --
    `solve_system` returns its iteration count and the call site bound it to `_`;
    the step controller knows what it rejected; the force-accept path from 4b
    counts nothing.  A run that takes 40x more steps than expected is currently
    indistinguishable, from the outside, from one that does not.

    The force-accept counter is the one to read first.  It counts steps accepted
    with an unbounded truncation error, and after 4d it should be zero on every
    circuit measured -- so a non-zero value is the run telling you that part of
    its own result is not error-controlled.
    """

    __slots__ = ('accepted_steps', 'rejected_steps', 'newton_iterations',
                 'force_accepts', 'order_drops', 'breakpoints_hit',
                 'min_step', 'max_step', 'solve_seconds', 'total_seconds')

    def __init__(self):
        self.accepted_steps = 0
        self.rejected_steps = 0
        self.newton_iterations = 0
        self.force_accepts = 0
        self.order_drops = 0
        self.breakpoints_hit = 0
        self.min_step = None
        self.max_step = None
        self.solve_seconds = 0.0
        self.total_seconds = 0.0

    def _note_step(self, dt):
        self.min_step = dt if self.min_step is None else min(self.min_step, dt)
        self.max_step = dt if self.max_step is None else max(self.max_step, dt)

    def as_dict(self):
        return {k: getattr(self, k) for k in self.__slots__}

    def __repr__(self):
        pct = (100.0 * self.solve_seconds / self.total_seconds
               if self.total_seconds else float('nan'))
        return (
            'accepted %d, rejected %d (%.1f%% of attempts), Newton iterations %d '
            '(%.1f per accepted step)\n'
            'force-accepts %d, order drops %d, breakpoints hit %d\n'
            'step %.4g .. %.4g s\n'
            'time %.3f s total, %.3f s in the Newton solve (%.1f%%)'
            % (self.accepted_steps, self.rejected_steps,
               100.0 * self.rejected_steps
               / max(1, self.accepted_steps + self.rejected_steps),
               self.newton_iterations,
               self.newton_iterations / max(1, self.accepted_steps),
               self.force_accepts, self.order_drops, self.breakpoints_hit,
               self.min_step if self.min_step is not None else float('nan'),
               self.max_step if self.max_step is not None else float('nan'),
               self.total_seconds, self.solve_seconds, pct))


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
         ## BACK TO 1e-12, matching `vabstol`, at gate D3-e.  The history is worth
         ## keeping because it is a clean example of a workaround outliving its
         ## defect.  This was 1e-12; it was raised to 1e-6 (Spectre's `vabstol`,
         ## SPICE's VNTOL) because on the 127-unknown leapfrog the timestep
         ## collapsed to 5 ns against a 39 ns cap -- the controller accepts on
         ## max(|lte|/etol) over ALL unknowns, most of that circuit's nodes carry no
         ## signal, and under `pointlocal` the relative reference on a quiet node
         ## tends to zero, so `etol` degenerated to TRTOL*abstol on numerical noise.
         ## Raising the floor cut the step count 5.4x.
         ##
         ## **That was treating the symptom.**  The cause was `pointlocal`, and
         ## `sigglobal` -- the default since decision D3 -- references every unknown
         ## to the largest signal in the circuit, so the reference cannot degenerate
         ## and the floor is never reached.  Measured at gate D3-e, `lte_vabstol` at
         ## 1e-6, 1e-9 and 1e-12 give **bit-identical** results under `sigglobal`:
         ## 403 steps on a pulsed RC and 601 on a circuit with a quiet node, at every
         ## one of the three values.  Under `pointlocal` the same change costs
         ## +8.5% and +9.2% -- which is what "load-bearing" looks like, and why the
         ## workaround was needed then and is not now.
         ##
         ## So the principled value returns at measured zero cost.  **If you select
         ## `relref='pointlocal'`, this floor becomes load-bearing again** and 1e-6
         ## may be the better choice for your circuit.
         Parameter(name='lte_vabstol',
                   desc='Absolute voltage tolerance for the local truncation error',
                   unit='V',
                   default=1e-12),
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
         ## DECISION D3, SECOND ATTEMPT -- `sigglobal` IS NOW THE DEFAULT, as in
         ## Spectre.  It was adopted, sent back by gate D3-a, and re-run.
         ##
         ## What sent it back: referencing the tolerance to the largest signal lets
         ## steps grow, and on an estimator carrying the trapezoidal `(-1)^n` mode
         ## that was enough to break the controller's response to `reltol` --
         ## accuracy stopped falling monotonically as the tolerance tightened.
         ## 4g(b) and 4i removed that contamination, and on re-run **all six
         ## integrator/formula combinations are monotone in both step count and
         ## error**.
         ##
         ## What it buys, measured at MATCHED ACCURACY rather than at matched
         ## `reltol` -- the two are not the same thing, and comparing step counts at
         ## equal `reltol` overstates the win because `sigglobal`'s error at a given
         ## `reltol` is ~1.5x larger:
         ##
         ##   euler 1.48-2.06x fewer steps, gear2 1.44-1.60x, trapezoidal 1.31-1.47x
         ##
         ## The figure previously recorded here was "1.7-2.5x fewer steps", taken at
         ## equal `reltol`.  That was a relabelling of the tolerance as much as a
         ## speedup; the honest worst case is 1.31x.
         Parameter(name='relref',
                   desc="Reference for the relative LTE tolerance: 'sigglobal' "
                        "(against the largest signal anywhere -- the default, as in "
                        "Spectre), 'pointlocal' (each unknown against itself, "
                        "pycircuit's historical behaviour), or 'alllocal' (against "
                        "its own past maximum)",
                   unit='',
                   default='sigglobal'),
         ## STAGE 12A -- Fang's acceptance band (DAC 2013 eq 15) and step-change
         ## damper (eq 16).  The defaults are the historical one-sided test, so
         ## nothing changes until a caller asks; see `StepController.set_lte_band`
         ## for why the paper's own 0.7/3.0 are not adopted as defaults.
         Parameter(name='lte_gamma_min',
                   desc="Lower edge of the LTE acceptance band, as a fraction of "
                        "the LTE tolerance. A step whose normalised error falls "
                        "below this is redone LARGER, rather than accepted as "
                        "wasted work. 0 disables the lower bound (default).",
                   unit='',
                   default=0.0),
         Parameter(name='lte_gamma_max',
                   desc="Upper edge of the LTE acceptance band, as a fraction of "
                        "the LTE tolerance. 1.0 is the historical rejection "
                        "threshold; above 1 fewer steps are redone.",
                   unit='',
                   default=1.0),
         Parameter(name='lte_eta',
                   desc="Relative limit on the change in step size between "
                        "consecutive steps, |dh| <= eta*h (Fang eq 16, ~0.15). "
                        "None leaves the change bounded only by zero stability.",
                   unit='',
                   default=None),
         Parameter(name='maxiter',
                   desc='Maximum number of iterations', unit='',
                   default=100),
         Parameter(name='integrator',
                   desc='Integration strategy (an Integrator instance, e.g. '
                        "EulerIntegrator(), TrapezoidalIntegrator(), "
                        "Gear2Integrator()); default EulerIntegrator()",
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
         ## STAGE 10.2 -- ask for a uniform output grid instead of resampling by
         ## hand.  `None` keeps the solver's own adaptive points, which is what
         ## every existing caller gets.
         Parameter(name='outputstep',
                   desc='Spacing of a UNIFORM output grid; None returns the '
                        "solver's own adaptive points. Results are interpolated "
                        'quadratically, to match the integrator -- see '
                        'resample_uniform',
                   unit='s',
                   default=None),
         Parameter(name='max_step',
                   desc='Largest accepted timestep; None means timestep. Raise it '
                        'to let the controller coast through quiescent intervals '
                        '(SPICE calls this TMAX and defaults it to tend/50)',
                   unit='s',
                   default=None),
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
        ## STAGE 4g(b).  The step before `_dt_last`.  It stays None until the run
        ## has taken two steps, which is exactly when `_qlast[2]` stops being the
        ## seeded initial charge and becomes a real past point -- so a single
        ## `None` tells an estimator both that the step and the charge are
        ## missing, and there is no second flag to keep in sync.
        self._dt_last2 = None
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
        linsolver = self._get_linearsolver()
        try:
            x_res, _iters = solver.solve_system(
                x0,
                refnode_removed(func, self.irefnode, self.toolkit),
                self.toolkit,
                self.par.reltol,
                abstol,
                xtol,
                self.par.maxiter,
                limiter=limiter_func,
                scaler=scaler,
                linsolver=linsolver,
                ## Stage 6: lets the solver name a node instead of a row index.
                row_names=reduced_row_names(self.cir, self.irefnode),
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
        
        ## Stage 6(c): this count was bound to `_` and discarded.
        stats = getattr(self, 'statistics', None)
        if stats is not None:
            stats.newton_iterations += int(_iters)

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
                "  * Transient(..., uic=True)  -- start from zeros (SPICE's 'use "
                "initial conditions'); note this is a Transient() argument, NOT a "
                "solve() one, or\n"
                "  * solve(x0=<vector>)        -- start from an operating point you "
                "supply.\n"
                "Both are explicit choices; substituting zeros silently is what this "
                "error replaced." % (exc,)) from exc

    ## STAGE 12B -- small helpers the coupled path needs, factored out of `_solve`
    ## rather than re-derived, so the two paths cannot drift apart on tolerances.

    ## The LTE tolerance multiplier.  `TRTOL` in this module, `lteratio` in
    ## Spectre: the LTE estimate is deliberately conservative, so the allowed
    ## truncation error is this many times the Newton-solve tolerance.
    LTERATIO = 7.0

    def _lte_abstol_vector(self):
        """The absolute floor of the LTE tolerance, per unknown.

        Note this is the `lte_*` pair, NOT the Newton `abstol` -- the two were
        split by stage 0.3d precisely because they are different quantities, and
        the coupled path must take the step-control flavour like every other
        controller here.
        """
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        return self.toolkit.concatenate(
            (self.par.lte_vabstol * ones_nodes,
             self.par.lte_iabstol * ones_branches))

    def _newton_abstol_vector(self):
        """The Newton RESIDUAL tolerance, per unknown.

        A node row of the residual is a current and a branch row is a voltage,
        hence `iabstol` on nodes and `vabstol` on branches.  This is the flavour
        for testing ``f``, not for testing an increment ``dx`` -- see
        :meth:`_newton_xtol_vector`, and `_newton` which builds both.
        """
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        return self.toolkit.concatenate(
            (self.par.iabstol * ones_nodes,
             self.par.vabstol * ones_branches))

    def _newton_xtol_vector(self):
        """The Newton SOLUTION tolerance, per unknown -- the other flavour.

        An increment on a node is a voltage and on a branch is a current, so the
        two vectors are transposed with respect to each other.  Getting this
        backwards is the same class of error stage 0.3d separated for the LTE
        tolerances: the numbers are dimensionally different quantities that
        happen to share default values (`iabstol` and `vabstol` are both 1e-12),
        so a swap is invisible until someone changes one of them.
        """
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        return self.toolkit.concatenate(
            (self.par.vabstol * ones_nodes,
             self.par.iabstol * ones_branches))

    def _residual_and_jacobian(self, x, t, provided_function=None):
        """``(f, J)`` at ``(x, t)`` using the current ``self._dt``.

        The same assembly `solve_timestep`'s inner function performs, reachable
        without running a Newton solve -- the coupled method needs one residual
        per iteration of its OWN loop.
        """
        C = self.cir.C(x, self.epar)
        q = self.cir.q(x, self.epar)
        self._q_cache = (x, q)
        iq, Geq = self.get_diff(q, C)
        u = self.cir.u(t, self.epar, analysis=self.par.analysis)
        if provided_function is not None:
            u = u + provided_function(t)
        f = self.cir.i(x, self.epar) + iq + u
        J = self.cir.G(x, self.epar) + Geq
        return (self.toolkit.array(f, dtype=float),
                self.toolkit.array(J, dtype=float))

    def fang_timestep(self, x_prev, t_prev, h, x_hist, refnode=gnd,
                      provided_function=None, gamma_min=0.7, gamma_max=3.0,
                      eta=0.15, maxiter=None, hmin=None, max_step=None,
                      hold_h=False, grid_locked=False):
        """One time point of Fang's coupled method: solve for ``(x, h)`` together.

        STAGE 12B.  Figure 4 of DAC 2013, and the structure is the substance:

          1. Solve the ordinary ``N`` circuit system at the current ``h`` and
             update the solution.  This is the existing Newton step, untouched.
          2. Estimate the LTE (eq 6) and find the controlling node.
          3. **If the LTE condition holds**, the step size needs no attention --
             check ordinary convergence and either finish or iterate again.
          4. **Only if it does not**, form the combined ``(N+1)`` system (eq 12)
             and solve for a solution update AND a step-size update at once.

        The (N+1) system is therefore NOT formed on every iteration, which is
        what makes the paper's overhead claim plausible.  There is no rejection
        path: Figure 3 has none, and the predicted ``h`` is only an initial
        guess.

        Eq (12) is solved by its Schur complement rather than by factorising an
        ``(N+1)`` matrix::

            dx0 = -J^-1 f_ckt          dxh = -J^-1 p
            dh  = -(f_lte + q^T dx0) / (q^T dxh + d)
            dx  = dx0 + dxh dh

        which needs two solves against the SAME ``J`` -- so with a factor/solve
        split the second is nearly free.  ``q^T`` has a single nonzero, so the
        two inner products are one multiply each.

        Returns ``(x, h, iterations, converged)``.
        """
        from pycircuit.circuit.stepcontroller import SolutionLTEController

        toolkit = self.toolkit
        n = self.cir.n
        irefnode = self.irefnode
        maxiter = self.par.maxiter if maxiter is None else maxiter
        hmin = getattr(self.par, 'minstep', 1e-18) if hmin is None else hmin
        if max_step is None:
            max_step = self.par.max_step
            if max_step is None or max_step <= 0:
                max_step = float('inf')

        ## GATE 12-4, last of the four inputs: honour a caller-injected step
        ## controller.
        ##
        ## A caller injects a step controller in order to control the steps.  On
        ## this path they cannot: the step-size law is Fang's, and the injected
        ## controller's own accept/predict logic is never consulted.  All the
        ## coupled path takes from it is `_reference`, the `relref` machinery,
        ## which every `StepController` has.
        ##
        ## So an injected controller is REFUSED unless it is one whose law this
        ## path actually implements.  Accepting it and using it only for
        ## `relref` would make `tran.step_controller = IntegralController()` look
        ## honoured while doing nothing -- the same class of defect as a
        ## documented feature that does not exist, which is what this path was
        ## until now (it silently built its own and ignored the caller's).
        ##
        ## NOT keyed on `lte_gradients`, which would be the obvious test and is
        ## wrong: `q^T` and `d` are implemented and gated but are NOT called on
        ## the shipped path, because sec. 3.4 replaced the eq (12) branch that
        ## used them.  Testing for a method nothing calls would pass controllers
        ## that cannot work and fail ones that can.
        injected = getattr(self, 'step_controller', None)
        if getattr(self, '_step_controller_is_auto', False):
            ## Auto-created by `_solve` on an earlier run of this object, not a
            ## caller's choice.  Nothing to honour and nothing to refuse.
            injected = None
        if injected is not None and not isinstance(injected, SolutionLTEController):
            raise ValueError(
                "the coupled (Fang) path cannot honour an injected %s: on this "
                "path the step size is solved from eq (6), so a controller's "
                "own accept/predict law is never used. Either drop the injected "
                "controller, pass a SolutionLTEController, or run with "
                "coupled_lte=False." % type(injected).__name__)

        ctrl = injected if injected is not None else \
            getattr(self, '_fang_controller', None)
        if ctrl is None:
            ctrl = self._fang_controller = SolutionLTEController()
        if getattr(ctrl, 'relref', None) != self.par.relref:
            ctrl.set_relref(self.par.relref)

        x = toolkit.array(x_prev, dtype=float).copy()
        ## The increment flavour: `dx0` is a solution update, not a residual.
        xtol = self._newton_xtol_vector()
        reltol = self.par.reltol

        ## Eq (16) bounds the step change BETWEEN ITERATIONS, and iterating it is
        ## how the step size collapses inside a single time point: 0.85 per
        ## iteration over `maxiter` iterations is seven decades.  Measured on the
        ## charge pump, `h` reached 8.75e-15 s at t=1.1e-5 before the solve gave
        ## up.  So the TOTAL excursion within one time point is bounded too, by
        ## the same window the standard controller allows for one step.
        from pycircuit.circuit.stepcontroller import (MIN_SHRINK_RATIO,
                                                      MAX_GROWTH_RATIO)
        h_entry = h
        h_floor = max(hmin, h_entry * MIN_SHRINK_RATIO)
        h_ceil = min(max_step, h_entry * MAX_GROWTH_RATIO)

        for it in range(maxiter):
            ## --- STAGE 1: the ordinary N system, at the current step size.
            self._dt = h
            t = t_prev + h
            f, J = self._residual_and_jacobian(x, t, provided_function)

            f_r = toolkit.delete(f, irefnode)
            J_r = toolkit.delete(toolkit.delete(J, irefnode, axis=0),
                                 irefnode, axis=1)
            dx0_r = toolkit.linearsolver(J_r, -f_r)
            dx0 = toolkit.insert(dx0_r, irefnode, 0.0)

            ## DEVICE LIMITING, the same the standard Newton applies.  Without
            ## it an undamped step across a diode's exponential is meaningless:
            ## the six nonlinear stress circuits all returned ~0 V where the
            ## standard path gives 8.9 V, and they did it silently -- the solve
            ## "converged", to the wrong thing.
            x_stage1 = self.cir.limit(x + dx0, x, self.epar)
            ## The limiter may shorten the step, so the convergence test must use
            ## what was actually taken, not what was asked for.
            dx0 = x_stage1 - x

            ## --- STAGE 2: has the step size earned any attention?
            h_hist = [hh for hh in (getattr(self, '_dt_last', None),
                                    getattr(self, '_dt_last2', None))
                      if hh is not None]
            etol = self._lte_tolerance(ctrl, x_stage1, x_prev, h_hist)

            eps_ok, err = self._lte_in_band(ctrl, x_stage1, x_hist, h_hist, h,
                                            etol, gamma_min, gamma_max)

            converged_x = bool(toolkit.alltrue(
                abs(dx0) < reltol * abs(x_stage1) + xtol))

            ## `hold_h` -- the step size is IMPOSED, not free.  A step truncated
            ## onto a breakpoint or onto `tend` has its size decided by where it
            ## must land, so there is nothing for the coupled system to SOLVE.
            ##
            ## Without it the truncation was pointless -- `fang_timestep` solved
            ## for its own `h` and walked straight off the edge again: 0 of 10
            ## pulse edges landed on, worst miss 1.24e-7 s, the whole rise time.
            ##
            ## BUT "DO NOT SOLVE FOR h" IS NOT "DO NOT CHECK THE ERROR", and
            ## conflating the two was a defect worth the same scrutiny as the one
            ## it replaced.  A held step was accepted blind, so its truncation
            ## error was governed by nothing: on the pulsed RC the maximum error
            ## sat at 1.465e-2 at BOTH reltol 1e-5 and 1e-6 -- identical across a
            ## decade of tolerance, the signature of a quantity no tolerance
            ## controls -- and the mean was 5.9x the standard path's at 1e-6,
            ## getting worse as the tolerance tightened.
            ##
            ## A held step whose error is over the band is reported to the caller,
            ## which shrinks and retries, exactly as the standard path does with a
            ## breakpoint-clamped step.  That is not the backup Figure 3 forbids:
            ## the paper's "no backup due to LTE" is about the step it SOLVES for,
            ## and this step's size was never its to choose.
            ## A held step whose error is over the band is reported so the
            ## caller can shrink and retry -- UNLESS the grid is locked.
            ##
            ## `fixed_timestep` is the caller stating that the output points are
            ## theirs, so shrinking is not an option available to us: the honest
            ## response to an over-tolerance step on a locked grid is to take it
            ## and let the run's accuracy be what the caller asked for, exactly
            ## as the standard path does. Conflating "truncated onto a
            ## breakpoint" with "grid imposed by the caller" broke
            ## `test_fixed_timestep_keeps_the_grid_on_the_coupled_path`: the
            ## retry shrank `h` and the uniform grid disappeared.
            if hold_h and not grid_locked and not eps_ok and err > gamma_max:
                return x_stage1, h, it + 1, False

            if hold_h or eps_ok or not h_hist or len(x_hist) < 2:
                ## The LTE condition holds (or cannot be evaluated yet, on the
                ## opening steps).  Nothing to solve for `h`; finish on the
                ## circuit equations alone.
                x = x_stage1
                if converged_x:
                    return x, h, it + 1, True
                continue

            ## --- The LTE condition failed, so the step size must move too.
            ##
            ## SEC. 3.4's APPROXIMATE NEWTON, NOT EQ (12), AND THE REASON IS
            ## MEASURED.  Eq (12) recovers `dh` from eq (14), whose denominator
            ## is `q^T dxh + d`.  Those two terms are the solution's sensitivity
            ## to the step size and the extrapolation's slope, and BOTH are
            ## approximately `dv/dt`: their difference is the truncation error's
            ## derivative, which is tiny by construction.  Measured on a driven
            ## RC at h = 1.6e-7: `q^T dxh = +1.818e9`, `d = -1.820e9`, denominator
            ## -2e6.  Three digits lost, and the SIGN of the denominator decided
            ## by the cancellation -- so `dh` saturated at the eta limit with an
            ## essentially arbitrary sign and the step drifted down four decades
            ## while `err` sat at 0.2, far BELOW the band that should have grown
            ## it.  Eq (12) computes a small quantity as the difference of two
            ## large ones; this is very likely what sec. 3.4 means by "the
            ## coupled nonlinear system sometimes is very sensitive to the change
            ## of step size".
            ##
            ## Eq (17) gets the new step from the error RATIO instead, which
            ## involves no cancellation at all.  `step_for_error_ratio` inverts
            ## the node polynomial rather than applying the (tau/eps)^(1/(n+1))
            ## power law, because that law only holds while h >> h_last -- see
            ## `extrapolation_error_weight`.
            from pycircuit.circuit._lte_kernels import step_for_error_ratio

            target = self._band_centre(ctrl, gamma_min, gamma_max)
            ratio = target / max(err, 1e-300)
            h_new = step_for_error_ratio(h, h_hist, ratio, 1.0 - eta, 1.0 + eta)

            ## WHAT THE STEP WANTS, ignoring every clamp.  Saturation has to be
            ## measured against this, not against the clamped result: once `h` is
            ## pinned at a bound the clamped `dh` is exactly 0.0, which is
            ## indistinguishable from "the step size has stopped moving" -- the
            ## definition of converged in eq (16) -- when in fact it stopped
            ## because it hit a wall.
            ##
            ## That hole is what let the first step after a pulse edge be
            ## accepted at 2.0e-7 s when it needed 3.55e-9 s, a factor of 56.
            ## The step came out of `fang_timestep` at exactly 0.2x its entry
            ## value -- MIN_SHRINK_RATIO, the within-time-point floor -- after 12
            ## iterations, reporting converged. It produced v = 0.033333 against
            ## an analytic 0.018731, a 78% single-step error, and the resulting
            ## 1.465e-2 was IDENTICAL at reltol 1e-5 and 1e-6 because nothing
            ## about it was tolerance-controlled.
            h_want = step_for_error_ratio(h, h_hist, ratio, 1e-6, 1e6)
            h_new = min(max(h_new, h_floor), h_ceil)
            dh = h_new - h

            ## DID THE CORRECTION SATURATE?  Eq (16), `|dh| <= eta*h`, is a
            ## convergence criterion meaning "the step size has stopped moving".
            ## A correction pinned AT the limiter has not stopped moving -- it
            ## was cut off -- and testing it with `<=` makes the two
            ## indistinguishable, because a clamped `dh` equals `eta*h` exactly.
            ##
            ## Measured on the pulsed RC: a step that needed to shrink tenfold
            ## just after a rising edge declared itself converged after a single
            ## 15% shrink, ran at h = 4.0e-8 s where the standard path used
            ## 4.4e-9 s, and left a maximum error of 1.465e-2 that was IDENTICAL
            ## at reltol 1e-5 and 1e-6 -- the signature of an error no tolerance
            ## governs.
            ## Only a clamped SHRINK counts.  A step that wants to grow and is
            ## held at the cap is in the normal state of every adaptive
            ## controller -- growth is bounded by zero stability, not by the
            ## error -- and treating that as unconverged made the run fail at
            ## t=1e-9 s with h driven to 9.5e-16, since the opening steps always
            ## want to grow faster than the cap allows.
            ## Only a thwarted SHRINK counts. A step that wants to grow and is
            ## held at the cap is the normal state of every adaptive controller
            ## -- growth is bounded by zero stability, not by the error -- and
            ## treating that as unconverged drove `h` to 9.5e-16 at t = 1e-9,
            ## because the opening steps always want to grow faster than allowed.
            saturated = h_want < h_floor * (1.0 - 1e-9)

            ## Eq (18): correct the solution already computed rather than
            ## re-solving at the new step size.  `dxh = -J^-1 p` reuses the
            ## factors from the stage-1 solve, which is the whole of sec. 3.4's
            ## "carries very little overhead".
            if dh != 0.0:
                p = self.residual_dh(x_stage1, t, h)
                p_r = toolkit.delete(p, irefnode)
                dxh_r = toolkit.linearsolver(J_r, -p_r)
                dxh = toolkit.insert(dxh_r, irefnode, 0.0)
                x = x_stage1 + dxh * dh
            else:
                x = x_stage1
            h = h_new

            ## `not saturated`, and a strict test: a step still moving at the
            ## limiter must keep iterating, and if the time point's own bound
            ## (`h_floor`) is what stopped it, the caller shrinks and retries.
            if converged_x and not saturated and abs(dh) < eta * h:
                return x, h, it + 1, True

        return x, h, maxiter, False

    def _band_centre(self, ctrl, gamma_min, gamma_max):
        """The normalised error eq (10) drives towards.

        Fang writes ``f_lte = eps_m - tau_m``, i.e. a target of exactly the
        tolerance.  With a BAND the sensible target is inside it rather than on
        either edge -- aiming at an edge makes every undershoot a violation, the
        defect gate 12A-1 measured as 3172 rejections against 1187 accepted
        steps.  The geometric centre is the point furthest from both edges in the
        ratio sense, which is the sense the step-size law works in.
        """
        return (gamma_min * gamma_max) ** 0.5

    def _lte_tolerance(self, ctrl, x_curr, x_last, h_hist):
        """``tau_m``, per unknown.  The paper does not specify it; this reuses
        the one every other controller here uses, so the coupled and standard
        paths are scored on the same scale."""
        ref = ctrl._reference(x_curr, x_last, not h_hist,
                              len(self.cir.nodes), self.toolkit)
        return self.LTERATIO * (self.par.reltol * ref + self._lte_abstol_vector())

    def _lte_in_band(self, ctrl, x_curr, x_hist, h_hist, h, etol,
                     gamma_min, gamma_max):
        """``(condition_holds, normalised_error)`` for eq (15)."""
        from pycircuit.circuit._lte_kernels import solution_lte

        if not h_hist or len(x_hist) < 2:
            return True, 0.0
        order = getattr(self.active_integrator, 'ORDER', 1)
        degree = min(order, len(x_hist) - 1, len(h_hist))
        if degree < 1:
            return True, 0.0
        lte = solution_lte(x_curr, list(x_hist[:degree + 1]),
                           list(h_hist[:degree]), h)
        import numpy as _np
        err = float(_np.max(abs(lte) / etol))
        return (gamma_min <= err <= gamma_max), err

    def residual_dh(self, x, t, h=None):
        """Fang's ``p = df_ckt/dh_m``, at fixed solution ``x``.

        STAGE 12B.  The residual assembled in :meth:`solve_timestep` is

            f = i(x) + iq(x, h) + u(t_{m-1} + h)

        so with ``x`` held fixed it depends on the step size through exactly two
        terms, and ``p`` is their sum:

          * ``d(iq)/dh`` from the integration coefficients, which eq (4) writes
            as explicit functions of ``h_m`` for this reason -- delegated to the
            ACTIVE integrator, so an order drop to Euler takes its derivative
            with it;
          * ``du/dt`` from the independent sources, since they are evaluated at
            ``t_{m-1} + h``.  On a driven circuit this is usually the larger of
            the two, and dropping it does not make the coupled system slightly
            wrong -- it makes it solve a different problem.

        ``i(x)`` is resistive and carries no ``h`` dependence at all.

        The solution's own dependence on ``h`` is deliberately absent: eq (12) is
        a block system of PARTIAL derivatives, and that coupling is what ``J``
        carries. Including it here would count it twice.
        """
        if h is None:
            h = self._dt
        q = self._q_at(x)
        h_last = getattr(self, '_dt_last', h)
        d_iq = self.active_integrator.companion_dh(q, self._qlast, h, h_last)
        d_u = self.cir.dudt(t, self.epar, analysis=self.par.analysis)
        return self.toolkit.array(d_iq, dtype=float) + \
            self.toolkit.array(d_u, dtype=float)

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
        ## STAGE 8(d) -- clear per-analysis element state BEFORE anything seeds it.
        ##
        ## Position matters and cost a test to learn: placed after the initial
        ## `accept_step(0.0, ...)` this wiped the very history that call had just
        ## seeded, so `TLine.G` saw an empty buffer and stamped the line as a DC
        ## SHORT -- v(p1) came out 1.0 where 0.5 is correct.  Elements that carry
        ## state must be reset before the run seeds them, not after.
        if hasattr(self.cir, 'reset_state'):
            self.cir.reset_state(self.epar)

        if coupled_lte:
            return self._solve_coupled(refnode, tend, x0, timestep,
                                       provided_function, analytical_eh,
                                       fixed_timestep=fixed_timestep)

        ## Respect a step controller injected by the caller (e.g. PIController);
        ## only fall back to the default IntegralController when none was set.
        if getattr(self, 'step_controller', None) is None:
            from pycircuit.circuit.stepcontroller import IntegralController
            self.step_controller = IntegralController()
            ## Marked so the coupled path can tell this apart from a controller
            ## the CALLER injected.  Without the distinction, any object that ran
            ## the standard path first presents this auto-created controller to
            ## the coupled path, which then refuses a controller nobody asked
            ## for -- 11 tests failed exactly that way.
            self._step_controller_is_auto = True
        ## ITEM 2+.3.  Applied to whichever controller is in use, including one
        ## the caller injected, and re-applied every run so the running maximum
        ## a global mode keeps cannot leak from a previous solve.
        self.step_controller.set_relref(self.par.relref)
        ## STAGE 12A -- same treatment, and for the same reason: applied to
        ## whichever controller is in use and re-applied every run, so a band set
        ## for one solve cannot leak into the next.
        self.step_controller.set_lte_band(self.par.lte_gamma_min,
                                          self.par.lte_gamma_max,
                                          self.par.lte_eta)

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
        ## The charge history is rebuilt per RUN, so the step history must be too.
        ## Without this a second `solve()` on the same object starts with a stale
        ## `_dt_last2` from the previous run while `_qlast[2]` is the freshly
        ## seeded initial charge -- breaking the one invariant 4g(b) relies on,
        ## that `h_last2 is not None` exactly when `q_last[2]` is a real past
        ## point.  `_dt_last` had the same staleness before 4g(b); nothing read it
        ## in a way that showed.
        self._dt_last = None
        self._dt_last2 = None
        
        X.append(copy(x))
        if hasattr(self.cir, 'accept_step'):
            self.cir.accept_step(0.0, X[-1], self.epar)
        
        timelist = []
        ## Stage 6(c).  Created per run, so a second `solve()` reports its own
        ## numbers rather than the sum of every run on this object.
        self.statistics = TransientStatistics()
        was_break_step = False
        force_order_drop = False
        _t_run_start = time.perf_counter()
        self._is_first_step = True
        self._no_history = True
        t = 0.0
        ## DECISION D2, 2026-08-01.  The clamp on how large an ACCEPTED step may
        ## grow.  It defaults to `timestep` -- the historical behaviour, and what
        ## `.tran tstep` means to most callers -- but it is now reachable.
        ##
        ## Measured before deciding.  On a run of ~5 tau the clamp is mostly
        ## irrelevant: above `timestep ~ 3e-4` the ERROR CONTROLLER becomes binding
        ## and `max dt` stalls at 2.97e-4 however much slack the clamp is given.
        ## But on a run of 100 tau, ~99% of it quiescent, it is clamped at every
        ## setting -- **1027 steps to traverse a dead-flat solution whose error is
        ## 4.4e-16**.  That cost is paid by exactly the circuits that idle, which is
        ## most mixed-signal ones.
        ##
        ## The DEFAULT is not changed: doing so would move every waveform in the
        ## package for a benefit only idle-heavy runs see.  SPICE's own default for
        ## the equivalent knob (`TMAX`) is `(tstop - tstart)/50`, which a caller can
        ## now ask for directly.
        max_step = self.par.max_step
        if max_step is None:
            max_step = timestep
        elif max_step < timestep:
            raise ValueError(
                "max_step (%g) is smaller than timestep (%g). The step can never "
                "exceed max_step, so this silently overrides the timestep you "
                "asked for; pass a smaller timestep instead, or max_step=None to "
                "use timestep as the clamp." % (max_step, timestep))

        ## STAGE 8(d) -- a delay element caps the step, per line.
        ##
        ## `TLine` interpolates its history at `t - TD`; with `dt` comparable to
        ## `TD` there is nothing to interpolate between, and the delay simply comes
        ## out wrong -- measured 2.00x TD at dt = 1e-9, 4.00x at 2e-9 and 8.00x at
        ## 5e-9 under `fixed_timestep`, with no warning of any kind.  The adaptive
        ## controller usually rescues it (1.01-1.05x), which is exactly why it went
        ## unnoticed: the defect only bites the configuration nobody checks.
        ##
        ## The cap is asked of the ELEMENTS rather than hard-coded here, so a future
        ## delay element gets it by implementing one method.
        element_cap = self.cir.max_timestep() if hasattr(self.cir, 'max_timestep') else None
        if element_cap is not None and element_cap < max_step:
            ## Under `fixed_timestep` the caller has taken the grid into their own
            ## hands, so silently substituting a finer one would be the wrong kind
            ## of help -- but running on regardless is worse, because the error is a
            ## WRONG DELAY and nothing else reports it.  Warn and obey, which is what
            ## the force-accept and non-convergence paths already do.
            if fixed_timestep:
                warnings.warn(
                    'transient: fixed_timestep=%g exceeds the %g s cap a delay '
                    'element needs (TD/2). The propagation delay will come out too '
                    'long -- measured 4x at twice the cap -- and nothing else will '
                    'report it. Use a timestep <= %g, or drop fixed_timestep.'
                    % (timestep, element_cap, element_cap), RuntimeWarning)
            else:
                max_step = element_cap
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
            
            ## STAGE 4h -- UNDER `fixed_timestep` THE GRID WINS.
            ##
            ## `dt` is loop-carried and truncation OVERWRITES it, while the restore
            ## at the bottom of this loop is guarded by `if not fixed_timestep`.  So
            ## a truncation used to be permanent, and because each later breakpoint
            ## truncated the already-shrunken `dt` again the step collapsed
            ## geometrically: a `VPulse` run that should take 30 steps took **292**,
            ## ending at `dt = 1.241e-19 s`.
            ##
            ## Restoring `dt` afterwards would fix the collapse but not the count --
            ## a `VPulse` at `per=1e-6` over `tend=3e-6` has ~12 edges, so the run
            ## would take `tend/timestep + 12` steps.  `fixed_timestep` exists so a
            ## caller can ask for exactly these output points, so the grid is what
            ## must be preserved: breakpoints no longer move it.
            ##
            ## What IS kept is the part that protects the integrator: crossing a
            ## breakpoint inside a step still drops the order for the next step, so
            ## no 2nd-order polynomial is fitted across a discontinuity.  It just
            ## costs no grid change.  A caller who needs edges resolved exactly
            ## wants the adaptive path with `max_step`, not a fixed grid.
            if fixed_timestep:
                ## `<=`, not `<`.  A breakpoint landing exactly ON a grid point is
                ## not a rounding curiosity here: with `td=1e-7` against a `1e-7`
                ## grid, 2 of the 9 edges in a 30-step VPulse run land exactly on
                ## one.  With a strict `<` those two produce no order drop at all --
                ## the step that ends on the edge does not see it, and the next
                ## iteration asks `next_event(t)` from the edge itself, whose
                ## fixed-point guard skips past it.
                was_break_step = next_t_break <= t + dt
            elif t + dt > next_t_break:
                dt = float(next_t_break - t)
                was_break_step = True
            else:
                was_break_step = False

            ## STAGE 12A -- was this step's size chosen by the controller, or
            ## imposed on it?  A truncated step is not LTE-limited, so its error
            ## says nothing about the integrator and Fang's lower bound must not
            ## try to grow it; see the guard in `IntegralController`.
            dt_clamped = was_break_step and not fixed_timestep

            if t + dt > tend:
                dt = tend - t
                dt_clamped = True
                ## A uniform grid divides `tend` exactly, so what is left at the end
                ## is floating-point residue rather than a step.  Turning it into
                ## one produced a final `dt` of 2.033e-20 s on a run whose other
                ## steps were 1e-6 -- 14 orders of magnitude down, and a step-size
                ## ratio no integrator should be asked to swallow.  Measured on the
                ## adaptive path for comparison: zero such steps, because a
                ## controller-chosen `dt` does not land on `tend` to within 1e-20.
                ## Hence the guard is scoped to fixed-step.
                if fixed_timestep and dt <= 1e-9 * timestep:
                    break
            
            self._dt = dt
            next_t = t + dt
            
            if was_break_step:
                self.statistics.breakpoints_hit += 1
            try:
                _t0 = time.perf_counter()
                x, feval, J, f = self.solve_timestep(X[-1], next_t, provided_function=provided_function)
                self.statistics.solve_seconds += time.perf_counter() - _t0
            except NoConvergenceError:
                ## STAGE 4h -- A FIXED GRID THAT CANNOT BE HONOURED MUST SAY SO.
                ## Shrinking is the only way to make progress, so it stays -- but
                ## under `fixed_timestep` the caller asked for exactly these output
                ## points and is no longer getting them.  Same failure class as 4b's
                ## force-accept: the result is still returned, and without this the
                ## caller has no way to learn that the grid they specified was
                ## abandoned partway through.
                if fixed_timestep:
                    warnings.warn(
                        'transient: Newton did not converge at t=%.6g s with the '
                        'requested fixed timestep %.6g s; falling back to %.6g s for '
                        'this step. The output grid is no longer uniform.'
                        % (t, timestep, dt * 0.25),
                        RuntimeWarning, stacklevel=3)
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
                    h_last2=getattr(self, '_dt_last2', None),
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
                    h_clamped=dt_clamped,
                    ## STAGE 12B -- accepted SOLUTION history, most recent first,
                    ## for Fang's eq (6) estimator.  The charge-based controllers
                    ## ignore it; `SolutionLTEController` extrapolates it to the
                    ## new time point and measures the deviation there.  Sliced
                    ## rather than passed whole so the controller cannot come to
                    ## depend on the full run being retained.
                    x_hist=X[-1:-4:-1],
                )

                if not accept and reject_count < MAX_REJECT:
                    self.statistics.rejected_steps += 1
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
                    self.statistics.force_accepts += 1
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
            self.statistics.accepted_steps += 1
            self.statistics._note_step(dt)
            if self._effective_method == 'EulerIntegrator' and \
                    type(self.base_integrator).__name__ != 'EulerIntegrator':
                self.statistics.order_drops += 1
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
            ## Roll before overwriting: _dt_last2 takes the value _dt_last is
            ## about to lose.  Reversing these two lines makes _dt_last2 equal
            ## _dt_last and the estimator silently differences the wrong grid.
            self._dt_last2 = self._dt_last
            self._dt_last = dt
            
            self._is_first_step = False
            self._no_history = False
            
            if not fixed_timestep:
                dt = next_dt
            
        X = self.toolkit.array(X[1:]).T
        timelist = self.toolkit.array(timelist)
        
        self.statistics.total_seconds = time.perf_counter() - _t_run_start

        self.result = CircuitResult(self.cir, x=X, xdot=None,
                                    sweep_values=timelist, 
                                    sweep_label='time', 
                                    sweep_unit='s')

        ## STAGE 10.2 -- resample onto a uniform grid if one was asked for.
        ## Done after the run rather than inside it, deliberately: the adaptive
        ## grid is what the error control is defined on, so the solver keeps
        ## choosing its own steps and only the REPORTED points change.
        outputstep = getattr(self.par, 'outputstep', None)
        if outputstep is not None:
            _grid, _Xg = resample_uniform(self.result.sweep_values,
                                          self.result.x, step=outputstep)
            self.result = CircuitResult(self.cir, x=_Xg, xdot=None,
                                        sweep_values=_grid,
                                        sweep_label='time', sweep_unit='s')
        ## Stage 6(c): reachable from the result, not only from the analysis, so a
        ## caller who kept only the waveform can still ask what produced it.
        self.result.statistics = self.statistics

        return self.result


    def _solve_coupled(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None, analytical_eh=True, fixed_timestep=False):
        ## STAGE 8(d) -- clear per-analysis element state BEFORE anything seeds it.
        ##
        ## Position matters and cost a test to learn: placed after the initial
        ## `accept_step(0.0, ...)` this wiped the very history that call had just
        ## seeded, so `TLine.G` saw an empty buffer and stamped the line as a DC
        ## SHORT -- v(p1) came out 1.0 where 0.5 is correct.  Elements that carry
        ## state must be reset before the run seeds them, not after.
        if hasattr(self.cir, 'reset_state'):
            self.cir.reset_state(self.epar)

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
        ## The charge history is rebuilt per RUN, so the step history must be too.
        ## Without this a second `solve()` on the same object starts with a stale
        ## `_dt_last2` from the previous run while `_qlast[2]` is the freshly
        ## seeded initial charge -- breaking the one invariant 4g(b) relies on,
        ## that `h_last2 is not None` exactly when `q_last[2]` is a real past
        ## point.  `_dt_last` had the same staleness before 4g(b); nothing read it
        ## in a way that showed.
        self._dt_last = None
        self._dt_last2 = None
        
        X.append(copy(x))
        if hasattr(self.cir, 'accept_step'):
            self.cir.accept_step(0.0, X[-1], self.epar)
        timelist = []
        
        self._is_first_step = True
        self._no_history = True
        t = 0.0
        ## Same opening ramp as `solve()` -- this path had the same defect.
        h = self._opening_step(timestep)
        ## DECISION D2 -- same clamp as `_solve`; see the note there.
        max_step = self.par.max_step
        if max_step is None:
            max_step = timestep
        elif max_step < timestep:
            raise ValueError(
                "max_step (%g) is smaller than timestep (%g); pass a smaller "
                "timestep instead, or max_step=None." % (max_step, timestep))

        ## STAGE 8(d) -- a delay element caps the step; see `_solve` for why.
        ## The coupled path is always adaptive, so the cap simply applies.
        element_cap = self.cir.max_timestep() if hasattr(self.cir, 'max_timestep') else None
        if element_cap is not None and element_cap < max_step:
            max_step = element_cap
        ## STAGE 12B -- the coupled path never created a statistics object, so
        ## `tran.statistics` raised AttributeError after any `coupled_lte=True`
        ## run and the two paths could not be compared on step counts at all.
        self.statistics = TransientStatistics()
        was_break_step = False
        force_order_drop = False
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

        ## Coupled time-stepping, Fang, "A New Time-Stepping Method for Circuit
        ## Simulation" (DAC 2013).  The solution and the step size are solved
        ## TOGETHER at each time point -- see `fang_timestep` for Figure 4's
        ## two-stage Newton and for why sec. 3.4's approximate step correction is
        ## used in place of eq (12).  The LTE is eq (6)'s solution-space estimate
        ## (`SolutionLTEController`), not the charge divided difference the rest
        ## of this module uses; `doc/fang_stage12_conclusions.md` sec. 3 has the
        ## reason, and it is the whole reason this path works at all.
        from pycircuit.circuit.nrsolver import NoConvergenceError
        MAX_LTE_ITERS = 10

        while t < tend:
            ## BREAKPOINTS, which this path did not have at all until now.
            ##
            ## The omission was worse here than it would be on the standard path,
            ## because Figure 3 has no rejection branch: a coupled step that runs
            ## past a pulse edge cannot back up from it.  It also quietly broke
            ## the invariant `p` depends on -- `TimeFunction.dfdt` takes the
            ## right-hand limit at a corner, which is exact ONLY because a step
            ## always starts on the corner rather than straddling it, and that is
            ## true only if the solver truncates to breakpoints.  It does now.
            if was_break_step or force_order_drop:
                ## Not "there is no history" -- the ring buffers keep rolling
                ## across a breakpoint.  It means "do not fit a 2nd-order
                ## polynomial through this point", which also drops the eq (6)
                ## predictor's degree, since that follows the active integrator.
                self._is_first_step = True
            force_order_drop = False

            next_t_break = self.cir.next_event(t)
            if next_t_break <= t + self.par.minbreak * max(abs(t), 1.0):
                next_t_break = self.cir.next_event(
                    t + (self.par.minbreak * 1e3) * max(abs(t), 1.0))

            ## GATE 12-4 -- `fixed_timestep` on the coupled path.
            ##
            ## The two are not in conflict so much as one simply wins: Fang's
            ## method exists to CHOOSE the step size, and `fixed_timestep` exists
            ## to say the caller has already chosen it.  So the grid is kept and
            ## the LTE equation is dropped on every step, exactly as it is for a
            ## breakpoint-truncated one -- the circuit is still solved coupled,
            ## it just has nothing to solve for.  Silently adapting anyway would
            ## return output points the caller did not ask for, which is what
            ## stage 4h fixed on the standard path.
            if fixed_timestep:
                h = timestep
                was_break_step = next_t_break <= t + h
            elif t + h > next_t_break:
                h = float(next_t_break - t)
                was_break_step = True
            else:
                was_break_step = False

            if t + h > tend:
                h = tend - t

            ## Co-determine (x, h): converge the circuit at h_curr, evaluate the
            ## LTE, and while it is above the band shrink the step (Gear-predicted
            ## by the controller) and re-solve.
            h_curr = h
            x_curr = copy(X[-1])
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
            ## STAGE 12B -- Fang's method, at last, replacing the rejection loop
            ## that used to live here.  `fang_timestep` solves for the solution
            ## and the step size together (Figures 3 and 4): there is no backup
            ## due to LTE, and `h_curr` is only the initial guess.
            ##
            ## The old loop is gone rather than kept behind a flag.  It re-solved
            ## the circuit from scratch whenever the LTE was over tolerance,
            ## which is a rejection loop with the retries hidden inside the time
            ## point -- the thing the citation in this module claimed not to be.
            ## A CONVERGENCE backup, which is NOT the thing Figure 3 forbids.
            ## The paper says "There is no backup due to LTE" -- the step size is
            ## solved rather than retried.  A Newton that fails to converge at
            ## all is a different failure, orthogonal to the LTE, and every
            ## production simulator retries it at a smaller step.  Removing this
            ## along with the LTE rejection loop is what broke the charge pump,
            ## the transformer inrush and the delayed avalanche.
            converged = False
            for _retry in range(MAX_LTE_ITERS):
                try:
                    x_curr, h_solved, _iters, converged = self.fang_timestep(
                        X[-1], t, h_curr, X[-1:-4:-1],
                        refnode=refnode, provided_function=provided_function,
                        hold_h=was_break_step or fixed_timestep,
                        grid_locked=fixed_timestep,
                        gamma_min=self.par.lte_gamma_min or 0.7,
                        gamma_max=(self.par.lte_gamma_max
                                   if self.par.lte_gamma_max != 1.0 else 3.0),
                        eta=self.par.lte_eta or 0.15,
                        hmin=minstep, max_step=max_step)
                except NoConvergenceError:
                    ## A device that cannot be evaluated at this step is the same
                    ## situation as a Newton that will not converge, and gets the
                    ## same response: shrink and try again, bounded.  Letting it
                    ## escape here would skip the backup entirely and abandon a
                    ## run that a smaller step would have completed.
                    converged = False
                ## STAGE 12-3.  Count the inner Newton iterations, including
                ## those of an attempt that then failed -- they are real work.
                ## Without this `newton_iterations` was flat zero on this path,
                ## so the only available cost comparison was per TIME POINT, and
                ## the coupled path runs several iterations per point (12 were
                ## measured at a pulse edge).  That made it look 14% cheaper per
                ## point than the standard path on rc-vsin, which was an artefact
                ## of the unit, not a property of the method.
                self.statistics.newton_iterations += int(_iters)
                if converged:
                    h_curr = h_solved
                    break
                h_curr *= 0.25
                if h_curr < minstep:
                    break

            if not converged:
                raise NoConvergenceError(
                    "Coupled transient: the coupled (x, h) Newton failed to "
                    "converge at t=%g s, with h reduced to %g s over %d "
                    "attempts. This is a convergence failure, not an LTE "
                    "rejection -- run with coupled_lte=False to use the "
                    "standard adaptive controller on this circuit."
                    % (t, h_curr, MAX_LTE_ITERS))

            t += h_curr
            self.statistics.accepted_steps += 1
            self.statistics._note_step(h_curr)
            ## The coupled path recorded only `accepted_steps`, so `order_drops`
            ## and `breakpoints_hit` read as zero on a circuit that was hitting
            ## ten pulse edges and dropping order at each -- a statistic that is
            ## silently always zero is worse than one that is absent.
            if was_break_step:
                self.statistics.breakpoints_hit += 1
            if self._effective_method == 'EulerIntegrator' and \
                    type(self.base_integrator).__name__ != 'EulerIntegrator':
                self.statistics.order_drops += 1
            timelist.append(t)
            X.append(copy(x_curr))

            if hasattr(self.cir, 'accept_step'):
                self.cir.accept_step(t, X[-1], self.epar)

            self._dt = h_curr
            self._dt_last2 = self._dt_last
            self._dt_last = h_curr
            self._is_first_step = False
            self._no_history = False
            self._iqlast = self.toolkit.concatenate((self.toolkit.array([self._iq]), self._iqlast))[:-1]
            self._qlast = self.toolkit.concatenate((self.toolkit.array([self._q_at(x_curr)]), self._qlast))[:-1]

            ## THE SOLVED STEP CARRIES FORWARD.  This is the whole point of the
            ## method -- `fang_timestep` returned the step size it solved for, and
            ## it becomes the initial guess for the next time point (Figure 3:
            ## "predict a step size", where the prediction is the previous
            ## answer).
            ##
            ## Writing anything else here is the defect gate 12B-0 found in the
            ## 2026-07 implementation, where `h = h_next` was followed two lines
            ## later by `h = h_curr` and the coupled system's answer was thrown
            ## away every step.  It came back the moment the old inner loop was
            ## deleted, because `h_next` was assigned inside it: the run then
            ## stayed at the opening step for its whole length, taking 151,176
            ## steps where the standard path takes 4,067, and -- the tell -- the
            ## step count did not move when `reltol` changed by two decades.
            h = min(max_step, max(h_curr, minstep))
            
        X = self.toolkit.array(X[1:]).T
        timelist = self.toolkit.array(timelist)
        self.result = CircuitResult(self.cir, x=X, xdot=None,
                                    sweep_values=timelist, 
                                    sweep_label='time', 
                                    sweep_unit='s')

        ## STAGE 10.2 -- resample onto a uniform grid if one was asked for.
        ## Done after the run rather than inside it, deliberately: the adaptive
        ## grid is what the error control is defined on, so the solver keeps
        ## choosing its own steps and only the REPORTED points change.
        outputstep = getattr(self.par, 'outputstep', None)
        if outputstep is not None:
            _grid, _Xg = resample_uniform(self.result.sweep_values,
                                          self.result.x, step=outputstep)
            self.result = CircuitResult(self.cir, x=_Xg, xdot=None,
                                        sweep_values=_grid,
                                        sweep_label='time', sweep_unit='s')
        return self.result



if __name__ == "__main__":
    import doctest
    doctest.testmod()
