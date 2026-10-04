# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import contextlib
import sys
import time

import numpy as np


from pycircuit.circuit.analysis import *
from pycircuit.circuit import _evalhint
from pycircuit.circuit import pcnr as _pcnr
from pycircuit.circuit._limiting import (limit_sync, stateful_limiters)
from pycircuit.circuit.simwarnings import (
    AccuracyWarning,
    ConvergenceWarning,
    UsageWarning,
    warn,
)
## The clamp the step controller applies to every accepted step, the force-accept
## path in `solve()` included: one bound, named once.  `stepcontroller` imports
## nothing from this package, so this is import-safe at module level.
## History: `doc/transient_history.md`, `transient.py`.
from pycircuit.circuit.stepcontroller import (MAX_GROWTH_RATIO,
                                              MIN_SHRINK_RATIO,
                                              normalised_error)
## THE THEMES `Transient` is assembled from (`_tran_*.py`, 2026-10-01; the
## review's O8) -- the stepping loop, its policy and the Parameters stay
## here.
from pycircuit.circuit._tran_newton import _StepNewton
from pycircuit.circuit._tran_branch import _BranchCheck
from pycircuit.circuit import _tran_companion, _tran_core, _tran_newton_c
from pycircuit.circuit._tran_companion import _CompanionModel
from pycircuit.circuit._tran_history import _RunHistory
from pycircuit.circuit._tran_predictor import _StagePredictor
from pycircuit.circuit._tran_initial import _InitialState
from pycircuit.circuit._tran_fang import _CoupledLTE
from pycircuit.circuit._tran_pcnr import _PCNRSteps
from pycircuit.circuit._tran_stages import _SequentialStages
from pycircuit.circuit._tran_glm import _NordsieckGLM
from pycircuit.circuit._tran_radau import _RadauStages
from pycircuit.circuit._tran_events import _StepEvents

## STAGE 2a -- BLAS thread control, discovered rather than required.
##
## Circuit matrices are small (n ~ 10^2), so a threaded LAPACK spends more time
## spawning and synchronising threads than doing the work.  The overhead scales
## with the core count the thread pool spans: on the leapfrog the transient runs
## 1.72x faster with BLAS limited to one thread on 4 cores, **14-20x** on 24.
##
## `threadpoolctl` is a dependency.  The import stays guarded anyway, and
## `blas_single_thread_available()` stays in the API, so a stripped environment
## degrades to threaded BLAS rather than failing to import a circuit simulator.
##
## ⚠ THIS LIMIT BELONGS TO THE TRANSIENT AND SHOULD NOT BE COPIED TO `DC` OR
## `AC`: the penalty is per-call thread-pool overhead, crossing over around
## 50-100 solves, NOT a property of problem size.  A transient runs thousands
## of small assemblies and solves in a Python loop; DC PREFERS threads at every
## size measured (limiting it would cost up to 2.3x), and AC never flips.
## History: `doc/transient_history.md`, `transient.py`.
try:
    from threadpoolctl import threadpool_limits as _threadpool_limits
except ImportError:
    _threadpool_limits = None


def blas_single_thread_available():
    """True if the BLAS thread limit can be applied from inside the process."""
    return _threadpool_limits is not None


## THE CONTROLLER, KEPT (speed round 8, stage 2c, 2026-10-04): building
## one rescans every shared library the process has loaded -- 3.4 M
## instructions, two steps of a 20-MosLevel1 transient, on every `solve`
## (`threadpool_limits` builds one per call) -- for a set that changes when
## a package is imported: a BLAS arrives with its package.  So one,
## rebuilt when `len(sys.modules)` moved; a library loaded through ctypes
## alone is not seen until then.  The limit itself is the same call on it.
## In a container: the leak detector reads module-level scalars.
_BLAS_CONTROLLER = {}


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
    n = len(sys.modules)
    rec = _BLAS_CONTROLLER.get('blas')
    if rec is None or rec[0] != n:
        from threadpoolctl import ThreadpoolController
        rec = _BLAS_CONTROLLER['blas'] = (n, ThreadpoolController())
    return rec[1].limit(limits=1, user_api='blas')

def resample_uniform(t, x, npoints=None, step=None, grid=None):
    """Interpolate a transient result onto a UNIFORM time grid.

    A transient returns the solver's own adaptive points, so their spacing is
    whatever the step controller chose.  Anything that needs a uniform grid,
    above all an FFT, resamples with this rather than with `np.interp`, which
    is LINEAR and throws away an order of accuracy before the transform.

    QUADRATIC, to match the integrator.  Three-point Lagrange through the
    interval's neighbours, the same reasoning `TLine._interpolate_history` gives
    for its own history lookup: a first-order interpolant feeding a second-order
    method injects error the method never made (measured 14-41x below linear
    where the grid resolves the signal, 1.3x where it barely does).  Where the
    adaptive grid is barely resolving the signal, no interpolant recovers what
    was not sampled, and the right fix there is a smaller `max_step`, not a
    cleverer resample.

    Give exactly one of ``npoints``, ``step`` or ``grid``.

    ``grid`` takes an explicit set of output times, for a caller whose window
    this function cannot otherwise express -- an FFT that needs a TRAILING
    window with ``endpoint=False``, because dropping the duplicated period
    boundary is what keeps the tones exactly on bins.  Neither ``npoints``
    (which spans ``t[0]..t[-1]`` inclusive) nor ``step`` (which starts at
    ``t[0]``) can say that.

    History: `doc/transient_history.md`, `resample_uniform`.
    """
    numpy = np

    t = numpy.asarray(t, dtype=float)
    x = numpy.asarray(x)
    given = sum(spec is not None for spec in (npoints, step, grid))
    if given != 1:
        raise ValueError('give exactly one of npoints, step or grid')
    if t.size < 3:
        raise ValueError('need at least 3 points to interpolate quadratically, '
                         'got %d' % t.size)

    if grid is not None:
        grid = numpy.asarray(grid, dtype=float)
    elif npoints is not None:
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

    Without these counts a run that takes 40x more steps than expected is
    indistinguishable, from the outside, from one that does not.

    On the COUPLED path a persistently failing point raises BY DESIGN (F13):
    its steps are solved, not rejected for error, so its only force-accept
    is a step the excursion veto (`max_dv_step`) still refused after
    `_CoupledSteps.max_reject` retries.  `rejected_steps` counts failed Newton
    attempts too, on every path.

    The force-accept counter is the one to read first.  It counts steps accepted
    with an unbounded truncation error, and it should be zero on every
    circuit measured -- so a non-zero value is the run telling you that part of
    its own result is not error-controlled.

    History: `doc/transient_history.md`, `TransientStatistics`.
    """

    __slots__ = ('accepted_steps', 'rejected_steps', 'newton_iterations',
                 'force_accepts', 'order_drops', 'breakpoints_hit',
                 'gmin_rescues', 'state_events_hit', 'state_event_cuts',
                 'branch_screens', 'branch_points', 'chord_fallbacks',
                 'min_step', 'max_step', 'solve_seconds', 'total_seconds')

    def __init__(self):
        self.accepted_steps = 0
        self.rejected_steps = 0
        self.newton_iterations = 0
        self.force_accepts = 0
        self.order_drops = 0
        self.breakpoints_hit = 0
        ## P18: failed time points recovered by the continuation ladder.
        self.gmin_rescues = 0
        ## E7: state events landed (a step cut to a declared crossing, the
        ## history restarted there) and the secant re-solves it took
        self.state_events_hit = 0
        self.state_event_cuts = 0
        ## `branch_check`: steps whose `rank C` fell below its structural
        ## value (screened), and those a re-solve found a second root of.
        ## ⚠ Both must be in `__slots__`: without them `_branch_count` raises
        ## on the first screen of every `solve`, and the check catches its
        ## own error and switches itself off.
        ## History: `doc/transient_history.md`, `TransientStatistics`.
        self.branch_screens = 0
        self.branch_points = 0
        ## `chord_jacobian`: steps whose chord iterations stopped contracting
        ## and were handed to the full Newton
        self.chord_fallbacks = 0
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
            'force-accepts %d, order drops %d, breakpoints hit %d, '
            'state events landed %d (%d cuts)\n'
            'step %.4g .. %.4g s\n'
            'time %.3f s total, %.3f s in the Newton solve (%.1f%%)'
            % (self.accepted_steps, self.rejected_steps,
               100.0 * self.rejected_steps
               / max(1, self.accepted_steps + self.rejected_steps),
               self.newton_iterations,
               self.newton_iterations / max(1, self.accepted_steps),
               self.force_accepts, self.order_drops, self.breakpoints_hit,
               self.state_events_hit, self.state_event_cuts,
               self.min_step if self.min_step is not None else float('nan'),
               self.max_step if self.max_step is not None else float('nan'),
               self.total_seconds, self.solve_seconds, pct))


#: What `Transient._glm_startup` did, for a caller that LINEARISES it (the
#: shooting's GLM period map, `_PeriodWalks._glm_startup_linearisation`): the
#: start time, the substep, the substep count, the substep states, each
#: substep's converged stages, and the Runge-Kutta method that took them.
from pycircuit.circuit._tran_glm import GLMStartupTrace  # noqa: F401  (re-exported)


class LastStep:
    """What the last `Transient.solve_timestep` left, for a caller that
    drives the steps itself (the shooting: `_InnerTransient.solve_timestep`,
    `_PeriodWalks`); `Transient.last_step` hands one out.

    ⚠ READ-ONLY AND LIVE, NOT A RECORD: each field is the transient's own
    object, not a copy.  The next step overwrites it, and the periodic
    gauge shift moves a stage method's `Y[-1]` -- which IS the step's `x`
    -- in place; so read it before the next step and copy what must
    outlive it.  (A snapshot would not be bit-identical to the reads it
    replaces, for exactly that reason.)

    Multistep: `C` (the charge Jacobian at the solution), `Geq` (the
    companion conductance), `coeffs` (the companion coefficients of the
    integrator that ACTUALLY ran -- Euler's on an order-dropped step).
    Stage methods: `Y` (stage states, full width), `K` (stage
    derivatives).  Nordsieck GLM: `nordsieck_in` (the vector the step
    entered with, after its rescale), `rho` (that rescale, `h / h_prev`; 1
    when none), `nordsieck` (the vector it left), `restarted` (it began
    afresh past the run's opening), `startup` (the last startup's
    `GLMStartupTrace`; None after `_glm_startup_override`)."""
    __slots__ = ('_tr',)

    def __init__(self, tr):
        self._tr = tr

    C = property(lambda self: self._tr._Cmat)
    Geq = property(lambda self: self._tr._Geq)
    coeffs = property(lambda self: self._tr._companion_coeffs)
    Y = property(lambda self: self._tr._rk_Y)
    K = property(lambda self: self._tr._rk_K)
    nordsieck_in = property(lambda self: self._tr._glm_Q_in)
    rho = property(lambda self: getattr(self._tr, '_glm_rho', 1.0))
    nordsieck = property(lambda self: self._tr._glm_Q[0])
    restarted = property(
        lambda self: bool(getattr(self._tr, '_glm_restarted', False)))
    startup = property(lambda self: self._tr._glm_startup_trace)


class TransientStepError(NoConvergenceError, RuntimeError):
    """A time point that could not be solved even at `minstep`, after the
    continuation rescue.  Both a `NoConvergenceError` and a `RuntimeError`,
    so a caller written against either catches it.

    History: `doc/transient_history.md`, `TransientStepError`."""


## ---------------------------------------------------------------------------
## THE STEP FAMILIES -- what differs between integrators in the one stepping
## loop of `Transient._solve` (`_SteppingLoop`).  The loop owns everything else: breakpoints, `tend`, the failure ladder and rescue, the
## excursion veto, the rejection budget and force-accept, state events, the
## bookkeeping of an accepted step.  A family supplies:
##
##   attempt(X, t, h, hold, provided_function) -> (x_new, h_taken, J)
##       one step from X[-1] (`hold`: its size is imposed), or raises
##       NoConvergenceError
##   judge(X, x_new, h, J, clamped)  ->  (ok, h_next)  its error test
##   h_after_force(h, h_next)  ->  the step after a force-accept
##   next_breakpoint(t), after_accept(landing), finish()
##
## and its retry budget at one time point: `max_retries` attempts in all, of
## which `max_reject` may be rejections for error -- then the step is
## force-accepted.
## History: `doc/transient_history.md`, `transient.py`.
## ---------------------------------------------------------------------------

class _StepFamily(object):
    """What the step families share unless they say otherwise: the
    breakpoints are the transient's own, a step is one `solve_timestep` from
    `X[-1]` at the size asked, a force-accept takes the step the error test
    proposed, and nothing is kept between accepts or after the run."""

    def next_breakpoint(self, t):
        return self.tr._next_breakpoint(t)

    def attempt(self, X, t, h, _hold, provided_function):
        x_new, _feval, J, _f = self.tr.solve_timestep(
            X[-1], t + h, provided_function=provided_function)
        return x_new, h, J

    def h_after_force(self, _h, h_next):
        return h_next

    def after_accept(self, landing):
        pass

    def finish(self):
        pass


class _LMMSteps(_StepFamily):
    """Linear multistep methods (euler, trap, theta, gear): a step is one
    companion-model Newton solve, judged by the step controller on the
    divided-difference charge LTE (`IntegralController` unless the caller
    injected one).

    ⚠ THREE REJECTIONS, THEN FORCE-ACCEPT.  Near a source corner the LTE
    estimate can stay above tolerance for arbitrarily small steps while the
    stored history is frozen; without a cap the step collapses and the solve
    grinds.  After `max_reject` the already-converged step is accepted (only
    its LTE was too high) with an order drop, so time advances and the history
    refreshes.  The controller's lower-band GROWTH retries (a too-accurate
    step redone larger, F14) are not rejections and do not count here."""

    max_reject = 3
    ## F14: the growth retries count here, bounding over/under alternation
    max_retries = 10

    def __init__(self, tr, run):
        self.tr, self.run = tr, run
        if getattr(tr, 'step_controller', None) is None:
            from pycircuit.circuit.stepcontroller import IntegralController
            tr.step_controller = IntegralController()
            ## Marked so the coupled family can tell this apart from a
            ## controller the CALLER injected.  Without the distinction, any
            ## object that ran an LMM run first presents this auto-created
            ## controller to a coupled run, which then refuses a controller
            ## nobody asked for.
            ## History: `doc/transient_history.md`, `_LMMSteps.__init__`.
            tr._step_controller_is_auto = True
        ## ITEM 2+.3 / STAGE 12A: applied to whichever controller is in use,
        ## including one the caller injected, and re-applied every run so a
        ## running maximum or a band cannot leak from a previous solve.
        tr.step_controller.set_relref(tr.par.relref)
        tr.step_controller.set_lte_band(tr.par.lte_gamma_min,
                                        tr.par.lte_gamma_max,
                                        tr.par.lte_eta)

    def judge(self, X, x_new, h, J, clamped):
        tr, run = self.tr, self.run
        return tr.step_controller.evaluate_step(**tr._lte_inputs(
            x_new, X[-1], J, h, run.abstol, run.max_step,
            ## STAGE 12A -- a truncated step is not LTE-limited, so Fang's
            ## lower bound must not try to grow it
            clamped=clamped, x_hist=X[-1:-4:-1]))

    def h_after_force(self, h, _h_next):
        ## the clamp every accepted step obeys -- 4b's point is that the
        ## force-accept must not bypass it
        return min(self.run.max_step, h * MAX_GROWTH_RATIO)


class _StageSteps(_StepFamily):
    """Runge-Kutta methods (radau, trbdf2, esdirk43) and Nordsieck GLMs:
    self-starting, so no divided-difference LTE -- a step is judged by the
    method's own embedded estimate, the filtered `_rk_est` the step leaves
    (a GLM delivers `_glm_error_estimate` through the same slot).  Accept at
    ``err <= 1``, the RMS over the unknowns of ``est / (reltol ref +
    abstol)``; the next step is ``h * clamp(0.9 err^(-1/(p+1)), 0.5, 2)``
    with ``p = EMBEDDED_ORDER`` (TR-BDF2 2(3) -> 1/3, Radau 5(3) -> 1/4), so
    one anomalous estimate cannot swing the step wildly.  No history to
    freeze at a corner, so more rejections are meaningful than for an LMM.

    ⚠ `ref` IS THE MULTISTEP FAMILY'S `relref`, default 'sigglobal': each
    unknown against the largest value of its unit group (node voltages,
    branch currents) so far.  A pointwise reference collapses to
    `lte_iabstol` wherever a source's branch current crosses zero, and that
    row then sets the step: on a hard-driven diode TR-BDF2 takes 1537
    attempts under 'pointlocal' and 324 under 'sigglobal', for the same
    error.  The price is where the pointwise reference had over-delivered:
    esdirk43 on a rectifier at reltol 1e-6 is 2e-6 V off under
    'pointlocal', 4e-5 under 'sigglobal'.  'pointlocal' is
    ``max(|x_new|, |x_n|)``, as in the multistep family and radau5.
    `TRTOL` is NOT applied: these are the embedded pairs' own estimates,
    not a divided difference to be over-estimated.

    The running maximum takes ACCEPTED points only (`X[-1]` at each
    judgement); a candidate counts for its own judgement, so a rejected
    step's overshoot never loosens the rest of the run.
    History: `doc/transient_history.md`, `_StageSteps`."""

    max_reject = max_retries = 12
    SAFETY, K = 0.9, 0.5

    def __init__(self, tr, run):
        self.tr, self.run = tr, run
        n = tr.cir.n
        self.keep = [i for i in range(n) if i != tr.irefnode]
        from pycircuit.circuit.stepcontroller import RELREF_MODES
        if tr.par.relref not in RELREF_MODES:
            raise ValueError("relref must be one of %r, not %r"
                             % (RELREF_MODES, tr.par.relref))
        self.relref = tr.par.relref
        self.n_nodes = len(tr.cir.nodes)
        self.running = None
        ## an adaptive run asks the step for its estimate; a fixed grid leaves
        ## the flag as the caller set it
        if not run.fixed:
            tr._rk_want_est = True

    def _reference(self, x_last, x_new):
        """The `ref` in ``reltol ref + abstol``, per `relref` (see the class
        note): the running maximum over the ACCEPTED points, `x_last` folded
        in, then the candidate `x_new` for this judgement only."""
        from pycircuit.circuit.stepcontroller import sigglobal_reference
        last = np.abs(np.asarray(x_last, dtype=float))
        new = np.abs(np.asarray(x_new, dtype=float))
        if self.relref == 'pointlocal':
            return np.maximum(last, new)
        self.running = last if self.running is None \
            else np.maximum(self.running, last)
        ref = np.maximum(self.running, new)
        if self.relref == 'alllocal':
            return ref
        return sigglobal_reference(ref, self.n_nodes)

    def judge(self, X, x_new, h, _J, _clamped):
        tk = self.tr.toolkit
        est = tk.array(self.tr._rk_est)
        wt = (self.tr.par.reltol * tk.array(self._reference(X[-1], x_new))
              + tk.array(self.run.abstol))
        ## (`normalised_error`: a 0/0 entry is 0, not a NaN that rejects the
        ## step and then grows it)
        ek = tk.array(normalised_error(np.asarray(est)[self.keep],
                                       np.asarray(wt)[self.keep]))
        err = float((tk.sum(ek * ek) / len(self.keep)) ** 0.5)
        order = int(self.tr.base_integrator.EMBEDDED_ORDER)
        grow = self.SAFETY * (err if err > 1e-16 else 1e-16) ** (
            -1.0 / (order + 1))
        return err <= 1.0, h * min(1.0 / self.K, max(self.K, grow))

    def finish(self):
        if not self.run.fixed:
            self.tr._rk_want_est = False


class _CoupledSteps(_StepFamily):
    """`coupled_lte=True`: Fang, "A New Time-Stepping Method for Circuit
    Simulation" (DAC 2013).  The solution and the step size are solved
    TOGETHER at each time point (`fang_timestep`, Figure 4's two-stage
    Newton), so there is no backup due to LTE -- the step size is solved, not
    retried, and `judge` always accepts; the solved step carries forward as
    the next guess (Figure 3).  The LTE is eq (6)'s solution-space estimate
    (`SolutionLTEController`), not the charge divided difference; see
    `doc/fang_stage12_conclusions.md` sec. 3.  A held step (landing, `tend`,
    a fixed grid, an event cut) has nothing to solve for; one whose LTE
    stays over the band returns unconverged and goes the way of a Newton
    failure -- smaller, and at `minstep` it raises.  So the only rejection
    the loop can make of a coupled step is the excursion veto.

    ⚠ TLINE WAVEFRONT ARRIVALS AND THE KINK DISCIPLINE ARE THIS FAMILY'S.
    A source corner reaching a delay line re-emerges at the far end TD later
    as a from-zero kink in an ALGEBRAIC variable that no element reports, so
    each corner schedules {corner + TD, corner + 2*TD} as breakpoints, and a
    landing empties the step ring (cold-start semantics for eq (6)).  History
    that straddles such a kink poisons the solution-space LTE (h cancels:
    err = 1/(TRTOL*reltol) = 1428.6 on the JAX probe; a ~10 fs crawl, 113
    points to cross one edge here).  The LMM family's integrator-side LTE
    decays with h and never livelocks, and applied there the reset MOVED the
    pulsed-RC comparison 9.8e-4 -> 5.4e-3 V (test_coupled_breakpoints) -- so
    it is gated on delay lines and on this family, by measurement."""

    max_reject = max_retries = 10

    def __init__(self, tr, run):
        import heapq
        self.tr, self.run = tr, run
        self._heapq = heapq
        ## R1 HARDENING: the sigglobal running maximum would otherwise
        ## survive on the cached controllers across runs of one object
        for ctrl in (getattr(tr, 'step_controller', None),
                     getattr(tr, '_fang_controller', None)):
            if ctrl is not None and hasattr(ctrl, 'set_relref'):
                ctrl.set_relref(tr.par.relref)
        ## the coupled step's once-per-run setup (`Transient._fang_setup`)
        tr._fang_run = tr._fang_setup(tr.par.coupled_method)
        ## F5: the 'auto' band sentinel resolved once, to Fang's values
        self.band = tr._coupled_band()
        self.tline_tds = sorted({float(e.iparv.TD)
                                 for _nm, e in tr.cir.elements.items()
                                 if type(e).__name__ == 'TLine'})
        self.arrivals = []
        self.seen_corners = set()

    def next_breakpoint(self, t):
        tr, run = self.tr, self.run
        t_break = tr._next_breakpoint(t)
        if self.tline_tds:
            if t_break < run.tend and t_break not in self.seen_corners:
                self.seen_corners.add(t_break)
                for td in self.tline_tds:
                    for k in (1, 2):
                        arrival = t_break + k * td
                        if arrival < run.tend:
                            self._heapq.heappush(self.arrivals, arrival)
            guard = t + tr.par.minbreak * max(abs(t), 1.0)
            while self.arrivals and self.arrivals[0] <= guard:
                self._heapq.heappop(self.arrivals)
            if self.arrivals:
                t_break = min(t_break, self.arrivals[0])
        self.t_break = t_break
        return t_break

    def attempt(self, X, t, h, hold, provided_function):
        tr, run = self.tr, self.run
        ## THE STEP MAY NOT GROW PAST THE BREAKPOINT: fang solves for its own
        ## h and could grow across the corner the entry h cleared (measured
        ## on the pulsed TLine: landing 8.9e-12 PAST it -- a straddled kink).
        ## Capping at the gap makes growth land exactly ON it.  A fixed grid
        ## keeps the caller's step uncapped.
        cap = run.max_step
        if not run.fixed and self.t_break < float('inf'):
            gap = self.t_break - t
            if gap > 0.0:
                cap = min(run.max_step, gap)
        gamma_min, gamma_max, eta = self.band
        x_new, h_solved, iters, converged = tr.fang_timestep(
            X[-1], t, h, X[-1:-4:-1],
            provided_function=provided_function,
            hold_h=hold, grid_locked=run.fixed,
            method=tr.par.coupled_method,
            gamma_min=gamma_min, gamma_max=gamma_max, eta=eta,
            hmin=tr.par.minstep, max_step=cap)
        ## STAGE 12-3: the inner Newton iterations are real work, counted
        ## for a failed attempt too
        tr.statistics.newton_iterations += int(iters)
        if not converged:
            raise NoConvergenceError(
                'coupled transient: the (x, h) Newton did not converge at '
                't=%g s, h=%g s (or a held step stayed over the LTE band)'
                % (t, h))
        return x_new, h_solved, None

    def judge(self, _X, _x_new, h, _J, _clamped):
        ## THE SOLVED STEP CARRIES FORWARD -- the whole point of the method
        ## (gate 12B-0: writing anything else here took 151,176 steps where
        ## the standard path takes 4,067, and the count did not move with
        ## `reltol`)
        return True, max(h, self.tr.par.minstep)

    def after_accept(self, landing):
        ## COUPLED KINK DISCIPLINE (see the class note): a landing empties
        ## the step ring, on delay-line circuits only
        if landing and self.tline_tds:
            self.tr._dt_last = None
            self.tr._dt_last2 = None

    def finish(self):
        self.tr._fang_run = None


class _SteppingLoop:
    """THE ONE STEPPING LOOP of `Transient._solve`, over a step family
    (`_LMMSteps`, `_StageSteps`, `_CoupledSteps`), the same for every
    method.  Each attempt is placed (`where`: a breakpoint, `tend`), taken
    (`take`: smaller on a Newton failure, the rescue at `minstep`), judged
    (`judge`: the family's error test, the excursion veto, the rejection
    budget and the force-accept), cut to a declared state event (`event`)
    and accepted (`accept`).  A phase returning False ends the attempt --
    `where` the run; the others retry from `where` at the size they set.

    What lives between attempts is here, not on the transient: `t`, `h`,
    `landing` (this attempt ends on a breakpoint, or the last one landed a
    state event) and `order_drop` (the last step was force-accepted) --
    both mean "do not fit a 2nd-order polynomial through this point" and
    take effect on the NEXT attempt, which, after a rejection, is the retry
    of the same point -- this time point's `rejects` / `point_retries`, the
    state-event secant cuts in flight (`ev_iter`, E7), and `imposed` (the
    next attempt's size was imposed: an event cut, the excursion veto); and
    the steps the run could not take as asked (`fixed_fallbacks`,
    `forced`), warned once after it.

    F14 (doc/transient_review_260820.md): a lower-band GROWTH retry -- the
    controller redoing a too-ACCURATE step larger -- is a voluntary redo,
    not a failure, and must not trip the force-accept (on a QUIESCENT
    circuit the opening ramp's growth retries would reach it and warn
    spuriously).  Over-tolerance rejections strictly shrink in every
    controller, growth retries return only behind a strict-growth guard, so
    `h_next > h` tells them apart; the family's `max_retries` bounds both
    against pathological alternation.

    Split out of `_solve` (2026-10-01), statement for statement.
    History: `doc/transient_history.md`, `Transient._solve`."""

    def __init__(self, tr, family, run, X, tend, timestep, provided_function):
        self.tr, self.family, self.run = tr, family, run
        self.X, self.timelist = X, []
        self.tend, self.timestep = tend, timestep
        self.fixed = run.fixed
        self.provided_function = provided_function
        self.minstep = tr.par.minstep
        ## The opening ramp exists to stop the ONE step the controller cannot
        ## check from dominating the run.  Under `fixed_timestep` there is no
        ## controller and the step is never adapted, so ramping would not
        ## open small and grow -- it would run the ENTIRE simulation at
        ## `timestep*1e-3`, a thousand times more steps for a result the
        ## caller explicitly asked to be uniform.
        ## History: `doc/transient_history.md`, `Transient._solve`.
        self.h = timestep if self.fixed else min(tr._opening_step(timestep),
                                                 run.max_step)
        self.t = 0.0
        self.landing = self.order_drop = False
        self.rejects = self.point_retries = 0
        self.fixed_fallbacks, self.forced = [], []
        self.ev_iter = 0
        self.imposed = False

    def execute(self):
        try:
            while self.t < self.tend:
                if not self.where():
                    break
                if not (self.take() and self.judge() and self.event()):
                    continue
                self.accept()
        finally:
            self.family.finish()

    def where(self):
        """1. Where this attempt may end.  False: the run is over.

        A breakpoint is a discontinuity (a VPulse corner): the step is
        truncated to land exactly on it, and the step after it drops the
        order -- `_is_first_step` means "do not trust a 2nd-order polynomial
        through this point", not "there is no history": the rings keep
        rolling, and the controller is handed `_no_history` (true only at
        the genuine start of a run) so a truncated step is still checked (a
        pulse train fires four per period; a sine none since stage 4g(a),
        `Sin.next_event`)."""
        tr, t, h = self.tr, self.t, self.h
        if self.landing or self.order_drop:
            tr._is_first_step = True
        self.order_drop = False
        t_break = self.t_break = self.family.next_breakpoint(t)
        if self.fixed:
            ## STAGE 4h -- UNDER `fixed_timestep` THE GRID WINS: a breakpoint
            ## does not move it (a truncation would be permanent and
            ## collapse the step geometrically), but crossing one still drops
            ## the order.  `<=`, TOLERANCED: `t` accumulates by `+= h`, so an
            ## edge exactly on a grid point is a float knife-edge (measured:
            ## the drop fired at the first edge and missed every later one).
            ## History: `doc/transient_history.md`, `Transient._solve`.
            self.landing = t_break <= t + h * (1.0 + 1e-9)
        elif t + h > t_break:
            h = float(t_break - t)
            self.landing = True
        else:
            self.landing = False
        ## STAGE 12A -- was this step's size chosen, or imposed?
        clamped = self.landing and not self.fixed
        if t + h > self.tend:
            h = self.tend - t
            clamped = True
            ## what a uniform grid leaves at `tend` is rounding residue, not
            ## a step (a final 2.033e-20 s step against 1e-6 ones)
            if self.fixed and h <= 1e-9 * self.timestep:
                self.h = h
                return False
        self.h, self.clamped = h, clamped
        ## a step whose size was decided by where it must land is HELD: the
        ## coupled family has nothing to solve for (F3: unheld, its final
        ## step grew past `tend` in 5 of 6 configurations)
        self.hold = clamped or self.fixed or self.imposed
        return True

    def take(self):
        """2. Take it: smaller on a Newton failure.  False: retry smaller."""
        tr, h = self.tr, self.h
        tr._dt = h
        try:
            self.x_new, self.h, self.J = tr._attempt_step(
                self.family, self.X, self.t, h, self.hold,
                self.provided_function)
        except NoConvergenceError:
            ## A failed attempt is retried smaller: a rejection in all but
            ## name, and counted as one (F13).
            tr.statistics.rejected_steps += 1
            if self.fixed:
                ## STAGE 4h -- a fixed grid that cannot be honoured must say
                ## so (counted; one warning a run, after the loop -- until
                ## 2026-10-01 one per step, each with its own `t`, which
                ## Python's once-per-location filter cannot collapse; the
                ## review's X8)
                self.fixed_fallbacks.append((self.t, h * 0.25))
            h = h * 0.25
            ## a step whose size was imposed is retried smaller STILL
            ## imposed: the coupled family then solves the circuit only
            ## (unheld, its (x, h) Newton cost +25-38 % Newton iterations on
            ## the pulsed RC for no accuracy)
            self.imposed = self.hold
            if h >= self.minstep:
                self.h = h
                return False
            ## P18 phase 3 (+P25): the LAST resort -- the point is re-solved
            ## at `minstep` through the junction-gmin -> gshunt ->
            ## pseudo-transient chain, and only a converged solution flows on
            h = self.minstep
            self.x_new, self.h, self.J = tr._rescue_step(
                self.family, self.X, self.t, h, self.hold,
                self.provided_function)
        return True

    def judge(self):
        """3. Judge it: the family's error test, then the excursion veto
        (`max_dv_step` / `max_di_step`).  False: retry at the size set."""
        if self.fixed:
            self.h_next = self.timestep
            return True
        tr, family, h = self.tr, self.family, self.h
        ok, h_next = family.judge(self.X, self.x_new, h, self.J, self.clamped)
        ratio = tr._excursion_ratio(self.x_new, self.X[-1]) if ok else None
        vetoed = ratio is not None and ratio > 1.0
        if vetoed:
            ok, h_next = False, h * max(MIN_SHRINK_RATIO, 0.9 / ratio)
        if not ok:
            growth = h_next > h
            if self.point_retries < family.max_retries and (
                    growth or (self.rejects < family.max_reject
                               and h > self.minstep)):
                tr.statistics.rejected_steps += 1
                self.point_retries += 1
                self.rejects += 0 if growth else 1
                self.h = max(h_next, self.minstep)
                ## the bound decided this size, so the coupled family must
                ## not solve it back up (unheld, its (x, h) Newton grew
                ## straight past the bound again)
                self.imposed = vetoed
                return False
            if not growth:
                ## FORCE-ACCEPT (4b): the error is still over tolerance and
                ## the budget is spent; the converged step is taken with an
                ## order drop, and the run says so -- the accepted error is
                ## unbounded.
                tr.statistics.force_accepts += 1
                self.order_drop = True
                h_next = family.h_after_force(h, h_next)
                ## (one warning a run, after the loop: X8)
                self.forced.append((self.t, h, self.rejects))
        self.rejects = self.point_retries = 0
        self.h_next = h_next
        return True

    def event(self):
        """4. A declared state event inside the step: cut to it (E7) by the
        secant on the fraction and re-solve, holding the cut step; the
        landed step restarts the history as a corner does.  Not from the
        initial condition (`len(X) > 1`).  False: retry the cut step."""
        tr = self.tr
        if tr._ev_rows is not None and not self.fixed and len(self.X) > 1:
            evs = tr._state_event_step(self.X[-1], self.x_new, self.h,
                                       self.t + self.h, self.ev_iter,
                                       self.minstep)
            if evs is not None and evs[0] == 'cut':
                self.ev_iter += 1
                self.h = evs[1] * self.h
                self.imposed = True
                return False
            self.ev_iter = 0
            if evs is not None and tr.EVENT_RESTART_HISTORY:
                self.landing = True
        return True

    def accept(self):
        """5. Accept it, and size the next attempt."""
        tr, h = self.tr, self.h
        self.imposed = False
        ## a step is judged a landing by where it ENDED: a coupled step that
        ## grew onto the corner is one in every way that matters
        on_break = self.t + h >= self.t_break * (1.0 - 1e-12)
        self.landing = self.landing or on_break
        t = self.t = self.t + h
        self.X.append(copy(self.x_new))
        self.timelist.append(t)
        tr.statistics.accepted_steps += 1
        tr.statistics._note_step(h)
        if on_break:
            tr.statistics.breakpoints_hit += 1
        if tr._effective_method == 'EulerIntegrator' and \
                type(tr.base_integrator).__name__ != 'EulerIntegrator':
            tr.statistics.order_drops += 1
        if hasattr(tr.cir, 'accept_step'):
            tr.cir.accept_step(t, self.X[-1], tr.epar)
        tr._roll_history(self.x_new, h, self.X)
        self.family.after_accept(self.landing)
        self.h = (min(self.h_next, self.run.max_step) if not self.fixed
                  else self.timestep)


class Transient(_StepNewton, _BranchCheck, _CompanionModel, _RunHistory, _StagePredictor, _InitialState, _CoupledLTE, _PCNRSteps, _SequentialStages, _NordsieckGLM, _RadauStages, _StepEvents, Analysis):
    """Simple transient analysis class.

    NOT REENTRANT: one Transient object runs one solve at a time.  The run
    threads state through the instance (_dt, _dt_last, _qlast, _iqlast,
    _is_first_step, _q_cache, the cached controllers), and fang_timestep
    communicates the trial step to get_diff via self._dt -- sharing an object
    across threads or interleaving solves corrupts all of it silently.
    (Review hygiene note; the fields are reset per run, so SEQUENTIAL solves
    on one object are fine.)

    The time step is adaptive (LTE-controlled, capped by `timestep_max`);
    `fixed_timestep=True` keeps the caller's uniform grid instead.

    **Backend parity** (doc/backend_parity_260821.md is the ledger): this
    class and `JAXTransient` share their parameter vocabulary and their
    defaults -- Gear-2 integration on the standard path, tend/50 as the
    default step cap, the same tolerances, band, `relref` modes, `uic`/`ic`
    machinery (the JAX class binds these methods verbatim), `outputstep`,
    and `provided_function`.  The coupled research path (`coupled_lte=True`,
    Fang DAC 2013) and PCNR run on both backends, sharing the Gear-2
    default since P22's state-row mask (eq (6) measured on state rows
    only; algebraic rows are slaved through the Jacobian).
    CPU-only, with cause: every STAGE and MULTIVALUE method (radau, trbdf2,
    esdirk43, the Nordsieck GLMs -- `JAXTransient` takes the three LMM
    companions, euler, trap and gear, and nothing else) and the
    `nrsolver`/`scaler`/`linearsolver` strategy objects -- per-iteration
    Python dispatch that a traced loop cannot host, so `JAXTransient`
    refuses them permanently (P17).  JAX-only by design: `solve_batched`
    (P20) -- one compiled kernel integrating every lane of a parameter sweep
    concurrently; this class gets no imitation of it, because a Python loop
    over `Transient` is already expressible and honest about its cost.

    Each step solves the method's companion system ``i(x) + iq + u(t) = 0``
    (a stage method: its stages) by Newton (`_newton`, or the strategy
    objects above), the integrator supplying the companion current `iq` --
    Gear-2 by default (the `integrator` Parameter).

    Example: an RC network charging from zero (`uic=True`; the default
    starts from the DC operating point, where it already sits), with a
    time constant of 9.9 ms and a final value of 9.9 V at `net2`:

    >>> circuit.default_toolkit = numeric
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> n2 = c.add_node('net2')
    >>> c['Is'] = IS(gnd, n1, i=10)
    >>> c['R1'] = R(n1, gnd, r=1)
    >>> c['R2'] = R(n1, n2, r=1e3)
    >>> c['R3'] = R(n2, gnd, r=100e3)
    >>> c['C'] = C(n2, gnd, c=1e-5)
    >>> tran = Transient(c, uic=True)
    >>> res = tran.solve(tend=10e-3, timestep=1e-4)
    >>> expected = 9.9 * (1 - np.exp(-10e-3 / 9.9e-3))
    >>> bool(abs(res.v(n2, gnd)[-1] - expected) < 1e-2 * expected)
    True
    """

    def _get_integrator(self):
        from pycircuit.circuit.integrator import (Integrator, Gear2Integrator)
        integrator = getattr(self.par, 'integrator', None)
        if integrator is None:
            ## Gear-2 is the shipped default (P6, owner decision), the same
            ## method as the JAX backend: half the steps and half the
            ## wall-clock of Euler on the same circuit at the same
            ## tolerance, and the conformance harness pins the pair.  The
            ## coupled path uses it too, with eq (6) measured on state rows
            ## only (`_state_row_mask`, P22).
            ## History: `doc/transient_history.md`, `Transient._get_integrator`.
            return Gear2Integrator()
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
         ## P2: same name and default as the JAX backend's Parameter.
         Parameter(name='TRTOL',
                   desc='LTE tolerance multiplier (lteratio in a commercial simulator): the '
                        'allowed truncation error is this many times the '
                        'Newton-solve tolerance',
                   unit='', default=7.0),
         ## NEWTON's x-tolerance on node rows, and nothing else.  Shared with `DC`,
         ## which uses the same value, so the operating point and the steps after it
         ## are solved to the same accuracy.
         ## ⚠ 1e-6 in DC, Transient, JAXTransient and PSS TOGETHER -- the four
         ## share one meaning and one default.  1e-12 is below what double
         ## precision can deliver on a node once ANY unknown in the circuit is
         ## large.  `lte_vabstol` is NOT this quantity.
         ## History: `doc/transient_history.md`, `Transient.parameters`.
         Parameter(name='vabstol',
                   desc='Absolute voltage error tolerance for the Newton solve',
                   unit='V',
                   default=1e-6),
         ## THE STEP CONTROLLER's tolerances, which are a different quantity: they
         ## apply to `lte = J^-1 Eg`, not to Newton's residual or its x-update.
         ##
         ## 1e-12.  `sigglobal` (the default `relref`) references every unknown to
         ## the largest signal in the circuit, so the reference cannot degenerate
         ## and the floor is never reached: `lte_vabstol` at 1e-6, 1e-9 and 1e-12
         ## gives **bit-identical** results under it.  **If you select
         ## `relref='pointlocal'`, this floor becomes load-bearing** (+8.5-9.2%
         ## steps at 1e-12) and 1e-6 may be the better choice for your circuit.
         ## History: `doc/transient_history.md`, `Transient.parameters`.
         Parameter(name='lte_vabstol',
                   desc='Absolute voltage tolerance for the local truncation error',
                   unit='V',
                   default=1e-12),
         Parameter(name='lte_iabstol',
                   desc='Absolute current tolerance for the local truncation error',
                   unit='A',
                   default=1e-12),
         ## What the RELATIVE part of the LTE tolerance is measured against --
         ## a commercial simulator's parameter of the same name, and its default.
         ##
         ## `pointlocal` references each unknown to itself, at this instant.  On a
         ## node carrying no signal that reference tends to zero, so the tolerance
         ## collapses to the absolute floor and the controller chases numerical
         ## noise on an idle node.
         ##
         ## `sigglobal` IS THE DEFAULT, as in a commercial simulator; all six
         ## integrator/formula combinations are monotone in both step count and
         ## error under it.  Measured at MATCHED ACCURACY (its error at a given
         ## `reltol` is ~1.5x larger, so equal `reltol` overstates the win):
         ##
         ##   euler 1.48-2.06x fewer steps, gear2 1.44-1.60x, trapezoidal 1.31-1.47x
         ## History: `doc/transient_history.md`, `Transient.parameters`.
         Parameter(name='relref',
                   desc="Reference for the relative LTE tolerance, for every "
                        "integrator: 'sigglobal' "
                        "(against the largest signal anywhere -- the default, as in "
                        "a commercial simulator), 'pointlocal' (each unknown against itself, "
                        "pycircuit's historical behaviour), or 'alllocal' (against "
                        "its own past maximum)",
                   unit='',
                   default='sigglobal'),
         ## STAGE 12A -- Fang's acceptance band (DAC 2013 eq 15) and step-change
         ## damper (eq 16).  The defaults are the historical one-sided test, so
         ## nothing changes until a caller asks; see `StepController.set_lte_band`
         ## for why the paper's own 0.7/3.0 are not adopted as defaults.
         ## THE DEFAULT IS THE STRING 'auto', NOT A NUMBER, AND NOT None --
         ## F5's lesson (doc/transient_review_260820.md).  Every documented
         ## value is meaningful (0.0 disables the lower bound, 1.0 is the
         ## historical threshold, None disables the damper), so a numeric or
         ## None default cannot be told apart from an explicit request for it.
         ## 'auto' resolves per path: the standard
         ## controller maps it to (0.0, 1.0, None) inside set_lte_band, the
         ## coupled path maps it to Fang's (0.7, 3.0, 0.15) in _coupled_band.
         ## Any explicit value is honoured verbatim on both paths.
         ## History: `doc/transient_history.md`, `Transient.parameters`.
         Parameter(name='lte_gamma_min',
                   desc="Lower edge of the LTE acceptance band, as a fraction of "
                        "the LTE tolerance. A step whose normalised error falls "
                        "below this is redone LARGER, rather than accepted as "
                        "wasted work. 0 disables the lower bound. 'auto' (the "
                        "default) means 0 on the standard path, 0.7 (Fang) on "
                        "the coupled path.",
                   unit='',
                   default='auto'),
         Parameter(name='lte_gamma_max',
                   desc="Upper edge of the LTE acceptance band, as a fraction of "
                        "the LTE tolerance. 1.0 is the historical rejection "
                        "threshold; above 1 fewer steps are redone. 'auto' (the "
                        "default) means 1.0 standard, 3.0 (Fang) coupled.",
                   unit='',
                   default='auto'),
         Parameter(name='lte_eta',
                   desc="Relative limit on the change in step size between "
                        "consecutive steps, |dh| <= eta*h (Fang eq 16, ~0.15). "
                        "None leaves the change bounded only by zero stability. "
                        "'auto' (the default) means None standard, 0.15 coupled.",
                   unit='',
                   default='auto'),
         ## STAGE 12B -- how the coupled path corrects the step size.
         ##
         ## 'approx'   Fang sec. 3.4: the new step comes from the error RATIO
         ##            (eq 17) and the solution is corrected by eq (18).  The
         ##            default, and the one with the measured record.
         ## 'bordered' Fang eq (12)/(14), RETIRED: asking for it raises.
         ## History: `doc/transient_history.md`, `Transient.parameters`.
         Parameter(name='coupled_method',
                   desc="Step-size correction for coupled_lte=True: 'approx' "
                        "(Fang sec 3.4, the only one; 'bordered' is retired "
                        "and raises)",
                   unit='',
                   default='approx'),
         ## THE RADAU COST TRANSFORM.  radau5's eig(A^-1) split: one real and one complex m x m solve per
         ## SIMPLIFIED-Newton iteration instead of the dense 3m system, and the
         ## dense full Newton when it stalls.  Measured on a diode-loaded RC
         ## ladder: fixed step 1.74x at m = 102 and 3.05x at m = 402, the answer
         ## the same to the Newton tolerance (1.2e-5 at the defaults; the
         ## dense-vs-transform test pins 1e-9 at tight ones); adaptive 1.5x
         ## (m = 52) to 3.1x (m = 402), at most one fallback per run.
         ## Simplified Newton can stall on a strongly nonlinear step, and the
         ## dense solve is the correctness reference, so it is on where the
         ## Jacobian is the expensive part: 'auto' (the default since
         ## 2026-10-01) where the circuit's compiled Jacobian reaches
         ## `AUTO_JACOBIAN_CODE` (the measurements are there).
         ## History: `doc/transient_history.md`, `Transient.parameters`.
         Parameter(name='radau_transform',
                   desc="Solve each Radau IIA step by the eig(A^-1) cost "
                        "transform: one real and one complex m x m solve per "
                        "iteration instead of the dense 3m one, falling back "
                        "to the dense solve if it stalls. 'auto' (default): "
                        "on where the circuit's compiled device Jacobians are "
                        "expensive (compact models: up to 68 % faster), off "
                        "otherwise; True / False force it. With pcnr=True on a "
                        "circuit PCNR applies to, PCNR's dense coupled solve is "
                        "used instead (warned where the transform was asked "
                        "for).",
                   unit='', default='auto'),
         ## THE CHORD JACOBIAN on the multistep step's Newton (`ChordNewton`):
         ## the Jacobian -- `G` and the companion's `C` -- evaluated and
         ## factored once at the step's seed, the residual alone (`i`, `q`,
         ## the source) every iteration; the full Newton from the same seed
         ## where it stops contracting.  The step ends at the converged point
         ## exactly as the full Newton's does (`jacobian_only`), so the step
         ## controller, the branch check and the shooting read the same
         ## state.  'auto' (the default): on where the circuit's compiled
         ## Jacobian reaches `AUTO_JACOBIAN_CODE`.  Measured (2026-10-01): a compact
         ## MOSFET (PSP) common-source stage, whose `G` and `C` are 92 % of
         ## the solve -- gear transient 5.5 -> 3.4 s, gear PSS 23.7 -> 15.5 s
         ## (`G` 3.18 -> 2.0 a step), the same iterations, the answers 1e-13
         ## apart; van der Pol -6 %, a comparator oscillator -1 %, a
         ## switching PWM loop +15 % (29 fallbacks), a diode mixer twice the
         ## Newton iterations for 39 % fewer `G`.
         ## History: `doc/transient_history.md`, `ChordNewton`.
         Parameter(name='chord_jacobian',
                   desc="Multistep methods (gear, trap, euler, theta): hold "
                        "each step's Newton Jacobian at its seed (the chord "
                        "method), iterating on the residual alone, and fall "
                        "back to the full Newton where it stops contracting. "
                        "'auto' (default): on where the circuit's compiled "
                        "device Jacobians are expensive (compact models: "
                        "5-32 % faster), off otherwise; True / False force "
                        "it.",
                   unit='', default='auto'),
         ## STAGE 13 -- PCNR instead of limiting, on the transient path too.
         ## Off by default for the same measured reason as on DC: gate 13-4 puts
         ## it at +60-80% per Newton iteration, for a consistency these circuits
         ## do not currently need.
         Parameter(name='pcnr',
                   desc='Use Predictor/Corrector Newton-Raphson instead of '
                        'limiting (Aadithya et al.); off by default',
                   unit='', default=False),
         Parameter(name='maxiter',
                   desc='Maximum number of iterations', unit='',
                   default=100),
         Parameter(name='integrator',
                   desc='Integration strategy (an Integrator instance, e.g. '
                        "EulerIntegrator(), TrapezoidalIntegrator(), "
                        "Gear2Integrator()); default Gear2Integrator() -- "
                        'the same method the JAX backend defaults to (P6)',
                   unit='',
                   default=None),
         Parameter(name='uic',
                   desc='Use initial conditions (skip DC OP computation)', unit='',
                   default=False),
         Parameter(name='state_events',
                   desc='Land the circuit\'s declared STATE events (a VSwitch '
                        'window edge, `Circuit.state_events()`): a crossing '
                        'inside an accepted step cuts the step to the '
                        'crossing and restarts the history there, as a '
                        'source corner does',
                   unit='', default=True),
         ## STAGE 10.3 -- SPICE's `.ic`, for `uic=True`.
         ##
         ## A start from a vector of zeros leaves a whole class of circuit
         ## unstartable: an LC tank at zero is AT an equilibrium and stays there,
         ## and a latch at zero sits on its metastable point.  Neither can be
         ## simulated at all without a way to say where it starts.
         ##
         ## Node voltages only: element initial conditions are the elements'
         ## own `ic` -- see `_initial_state`.
         ## History: `doc/transient_history.md`, `Transient.parameters`.
         Parameter(name='ic',
                   desc="Initial node voltages for uic=True, as {node: volts}. "
                        "Node may be a name or a Node instance.",
                   unit='V',
                   default=None),
         Parameter(name='minbreak',
                   desc='Minimum time difference for breakpoint events', unit='s',
                   default=1e-14),
         Parameter(name='bypass',
                   desc='Enable device model bypassing', unit='',
                   default=False),
         Parameter(name='minstep',
                   desc='Minimum timestep to prevent infinite loops', unit='s',
                   default=1e-18),
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
         ## The commercial-simulator-class VOLTAGE CHECK: on a purely resistive/
         ## algebraic network (a designer exploring an amplifier topology
         ## with Rs and controlled sources, no reactances yet) NO error
         ## estimator has anything to measure -- the charge-based LTE is
         ## identically zero and P22's mask excludes algebraic rows from the
         ## coupled band -- so the run samples at the step cap and nothing
         ## reports how coarse that is.  Bounding the PER-STEP node-voltage
         ## change controls output resolution directly, tracks the actual
         ## waveform (a slewing output gets dense points, a quiet one does
         ## not), and is robust where solution-LTE was not: |dv| ~ h*slew is
         ## h-proportional by construction, so it cannot h-cancel.
         Parameter(name='max_dv_step',
                   desc='Per-step node-voltage excursion bound (the commercial-'
                        'simulator-style voltage check), as a FACTOR: the bound '
                        'is max_dv_step * lte_vabstol (e.g. 2e11 at the '
                        "default 1e-12 bounds steps to 0.2 V); 'auto' "
                        'derives it from sampling theory (points_per_period), so it scales with the '
                        'tolerance family exactly as the LTE does. Factors '
                        'below 1 clamp to 1 (a bound below the Newton '
                        'accuracy measures noise). None disables. Controls '
                        'output resolution on circuits where no LTE exists '
                        '(resistive/algebraic networks).',
                   unit='',
                   default=None),
         Parameter(name='max_di_step',
                   desc='Per-step branch-current excursion bound, as a '
                        'FACTOR times lte_iabstol; the current-row sibling of '
                        "max_dv_step, same clamp-at-1 floor. 'auto' derives "
                        'the bound from sampling theory -- see '
                        'points_per_period. None disables.',
                   unit='',
                   default=None),
         ## THE SCIENTIFIC SETTING (owner request): N points per period of a
         ## sinusoid is per-step excursion <= 2*pi*swing/N, so 'auto' bounds
         ## the step to (2*pi/N) * max(static source swing, running signal
         ## maximum).  The source term (the signal_scale element hook)
         ## anchors the bound at signal BIRTH -- every running-reference
         ## scheme h-cancels there, the trap this file now documents three
         ## times -- and the running per-unit-group maximum grows the bound
         ## as an amplifier's output reveals gain the sources cannot know.
         Parameter(name='points_per_period',
                   desc="Sampling density behind max_dv_step/max_di_step = "
                        "'auto': at least this many points per period of a "
                        'full-swing sinusoid (per-step excursion '
                        '2*pi*swing/N). 64 resolves harmonics to ~20th '
                        'order with Nyquist margin.',
                   unit='',
                   default=64),
         ## The step cap (`timestep_max`) is DECOUPLED from `timestep` (owner
         ## decision): a shared value makes the step count on gentle circuits a
         ## property of the requested output density rather than of the error
         ## control.  `timestep` only sets the opening-step scale and the
         ## fixed_timestep grid.
         ## History: `doc/transient_history.md`, `Transient.parameters`.
         Parameter(name='timestep_max',
                   desc='Largest accepted timestep; None means tend/50, the '
                        'SPICE TMAX default. Decoupled from timestep, which '
                        'only sets the opening-step scale and the '
                        'fixed_timestep grid',
                   unit='s',
                   default=None),
         ## STAGE 3.  The opening step, which the controller must accept
         ## unevaluated because there is no history to difference against.  `None`
         ## means `timestep * 1e-3`; see `Transient._opening_step` for why opening
         ## at `timestep` made `reltol` unable to influence the answer at all.
         Parameter(name='firststep',
                   desc='Size of the first timestep; None means timestep*1e-3, '
                        'and a larger value is capped at timestep. The first '
                        'step cannot be error-checked, so taking it large lets '
                        'its error dominate the whole run.',
                   unit='s',
                   default=None),
         Parameter(name='bypasstol',
                   desc='Bypass tolerance for device models', unit='V',
                   default=None)]

    ## `irefnode` is not an argument: passing it fails loudly (F18).
    ## History: `doc/transient_history.md`, `Transient.__init__`.
    def __init__(self, cir, toolkit=None, **kvargs):
        ## `toolkit` is forwarded: `Transient(cir, toolkit=X)` runs on X, not
        ## on `cir.toolkit`.
        ## History: `doc/transient_history.md`, `Transient.__init__`.
        super(Transient, self).__init__(cir, toolkit=toolkit, **kvargs)

        ## ⚠ PCNR BOOKKEEPING IS INITIALISED HERE, NOT ONLY IN `_solve`.
        ## `_solve` resets these per analysis, which is right -- but SHOOTING
        ## never calls `_solve`: it drives `solve_timestep` directly on its own
        ## grid, and the PCNR paths read these attributes.  Defining them at
        ## construction makes every entry point safe, and `_solve`'s reset
        ## still gives each analysis a clean count.
        ## History: `doc/transient_history.md`, `Transient.__init__`.
        self._reset_pcnr_counts()

        ## THE REFERENCE NODE IS `self.irefnode`, SET HERE AND NOWHERE ELSE
        ## SPONTANEOUSLY -- `solve()` overwrites it from its `refnode`
        ## argument.  One fact, one home: a second copy (a `refnode=gnd`
        ## default somewhere) gives two reference nodes in one solve (F7).
        ## History: `doc/transient_history.md`, `Transient.__init__`.
        self.irefnode = self.cir.get_node_index(gnd)

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

    
    
    def _opening_step(self, timestep):
        """The size of the first step of a run.

        STAGE 3.  The step controller accepts the first step unevaluated, because
        with no history there is nothing to difference and no truncation error can
        be estimated, so the opening step is the **only unchecked** step in the
        run.  Opened large, its error dominates everything after it: at
        `timestep`, an RC step response's global error is the same to five digits
        from reltol 1e-3 to 1e-6, and for every integrator.

        Opening at `timestep * 1e-3` costs one cheap step and leaves the controller
        to grow the step from there, which it does geometrically, so the ramp is
        paid off within a handful of steps.

        History: `doc/transient_history.md`, `Transient._opening_step`.
        """
        firststep = self.par.firststep
        if firststep is None:
            return timestep * 1e-3
        if firststep <= 0:
            raise ValueError(
                "firststep must be positive, not %r; pass None to use the default "
                "ramp of timestep*1e-3" % (firststep,))
        ## capped at `timestep`, the opening-step scale: asking for more is more
        ## likely a mistake than an intent
        return min(firststep, timestep)

    @property
    def last_step(self):
        """What the last step left, read-only and live (`LastStep`)."""
        return LastStep(self)

    def solve_timestep(self, x0, t, provided_function=None):
        """One step of the run's method from `x0` to time `t`; the step
        length is `self._dt`, set by the caller (the stepping loop, or the
        shooting's inner transient).  Dispatched on the integrator: a
        Nordsieck GLM (`_solve_timestep_glm`), any Runge-Kutta tableau
        (`_solve_timestep_rk`), PCNR when asked for and a device takes part
        (`_solve_timestep_pcnr`), else the multistep companion Newton below.

        Returns ``(x, feval, J, f)``: the new state and the step's Jacobian at
        it (the step controller and the shooting walks read `J`).  `feval`
        and `f` are read by no caller -- None on every path but PCNR's, and
        `f` is None on the multistep path by design (`jacobian_only`).
        """
        ## THE STEP'S SOURCE MEMO lives exactly as long as this call
        ## (`_source_at`, speed round 4): `u(t)` assembled once per time
        ## within the step, and nothing of it left for the next one -- a
        ## later caller at the same `t` (the shooting re-entering a period,
        ## a source whose state `accept_step` moved) assembles afresh.
        self._u_memo = {} if _tran_companion.U_MEMO else None
        try:
            return self._solve_timestep(x0, t, provided_function)
        finally:
            self._u_memo = None

    def _solve_timestep(self, x0, t, provided_function=None):
        """`solve_timestep`'s body: the dispatch on the integrator, and the
        multistep companion Newton."""
        from pycircuit.circuit.integrator import RungeKuttaIntegrator
        ## ⚠ THE STEP STARTS FROM ITS ENTERING POINT, NOT FROM THE LAST
        ## ATTEMPT'S DEVICE STATE.  A stateful limiter (`Diode`) reads `i` /
        ## `G` as the tangent at its stored `_vlim`, and a REJECTED attempt
        ## leaves that at its own end: the retry's explicit first stage read
        ## `i(x_n)` there, and adaptive TR-BDF2 on a hard-driven diode took
        ## 1653 steps where the state-free twin takes 1117 (ESDIRK43 506 /
        ## 473).  PSS re-enters every period at a new `x0` the same way.
        ## After an accepted step the state already sits exactly on the
        ## entering point, so a run without rejections is unchanged.
        lims = getattr(self, '_stateful_lims', None)
        if lims is None:
            lims = self._stateful_lims = stateful_limiters(self.cir)
        if lims:
            limit_sync(self.cir, x0, self.epar, lims)
        if getattr(self.base_integrator, 'is_multivalue', lambda: False)():
            ## a Nordsieck general linear method: r = p + 1 values per unknown
            ## carried between steps, DIRK-like sequential stages, a starting
            ## vector computed at the first step -- see `_solve_timestep_glm`
            return self._solve_timestep_glm(x0, t, provided_function)
        if isinstance(self.base_integrator, RungeKuttaIntegrator):
            ## ANY Runge-Kutta method: the one tableau-driven stage step, which
            ## picks the DIRK-sequential or fully-implicit-coupled path from the
            ## tableau's structure.  PCNR (when `par.pcnr`) is the limiting on
            ## both paths -- see `_rk_stage_pcnr` and `_rk_step_coupled_pcnr`;
            ## it flows to shooting too, since the inner transient calls this
            ## same method.
            ## History: `doc/transient_history.md`, `Transient.solve_timestep`.
            return self._solve_timestep_rk(x0, t, provided_function)
        ## STAGE 13 -- the PCNR path, when asked for and when the circuit has a
        ## device that participates.  A circuit with no PCNR junction falls
        ## through: there is nothing for the method to do, and refusing would be
        ## a worse answer than solving it the ordinary way.
        if self.par.pcnr:
            ## Gate PARTICIPATION on the device records, not on the
            ## pnj-only pair view: that view exists for the gmin ladders
            ## and is empty for a circuit of pure fetlim/limvds devices,
            ## so gating on it would let `pcnr=True` on a MOSFET
            ## differential pair fall through to the ordinary solver
            ## SILENTLY.
            ## History: `doc/transient_history.md`, `Transient.solve_timestep`.
            if _pcnr.pcnr_devices(self.cir):
                ## Same fallback as DC(pcnr=True): a PCNR failure on one
                ## timestep falls through to the ordinary step solver rather
                ## than ending the transient.  See dcanalysis and
                ## `_pcnr_attempt`.
                out, ok = self._pcnr_attempt(
                    lambda: self._solve_timestep_pcnr(x0, t, provided_function),
                    lambda exc: self._pcnr_failed('PCNR', t, exc))
                if ok:
                    return out
            else:
                ## Asked for, and no device declares a probe.  Falling
                ## through is right -- refusing would be a worse answer --
                ## but it must SAY SO, or `pcnr=True` and `pcnr=False` are
                ## indistinguishable from outside.
                ## History: `doc/transient_history.md`, `Transient.solve_timestep`.
                self.pcnr_status = 'no-participants'

        n=self.cir.n
        dt = self._dt
        
        def func(x):
            return self._residual_and_jacobian(x, t, provided_function)

        def jacobian_only(x):
            """The converged-point evaluation, without the residual nobody reads.

            ITEM 2+.2.  Newton needs `f` on every iteration -- it is what drives the
            update.  The *final* evaluation at the converged point is different: its
            `f` is unpacked by both callers and then never referenced again (`solve`
            at the `x, feval, J, f = ...` site, and the coupled path likewise). Only
            `J` is consumed, by the step controller.

            So `cir.i(x)` and `cir.u(t)` are not assembled here.  `C`, `q` and `G`
            are all still needed -- `C` and `G` build `J`, and `q` feeds the charge
            cache and the history roll -- so this is not a cheaper approximation of
            the same work, it is the same work minus two vectors that have no
            consumer.

            `f` is returned as None rather than a zero vector, so any future caller
            that starts reading it fails loudly instead of silently using zeros --
            which is the whole lesson of stage 1.

            History: `doc/transient_history.md`, `Transient.solve_timestep`.
            """
            ## (the converged point's session, which the branch screen's `C`
            ## read at this state shared: `_conv_session`)
            with _evalhint.evaluating(session=conv):
                ## (the evaluate core, `_tran_core`: C, q and G in one call
                ## where it serves, the same state; the path below where not)
                r = _tran_core.evaluate(self, x, t, provided_function, 'j')
                if r is not None:
                    return r
                _iq, Geq = self._companion_at(x)
                J = self.cir.G(x, self.epar) + Geq
            return None, _tran_companion._as_float(J, self.toolkit)

        def residual_only(x):
            """The step residual alone, for the chord iterations
            (`chord_jacobian`): `i`, `q` and the source at `x`, the companion
            current from `q` (`get_diff`).  Every multistep companion's current
            is a function of the charges alone -- the conductance it also
            returns is the only use of `C`, and it is discarded here, so the
            seed's `C` (`_Cmat`, set by the evaluation at the seed) stands
            in.  The step's state is the full Newton's after it all the same:
            `jacobian_only` evaluates the converged point."""
            with _evalhint.evaluating('q', 'i'):
                ## (the evaluate core: i and q in one call, the conductance
                ## from the held `_Cmat`, where it serves)
                r = _tran_core.evaluate(self, x, t, provided_function, 'f')
                if r is not None:
                    return r
                q = self.cir.q(x, self.epar)
                iq, _geq = self.get_diff(q, self._Cmat)
                u = self._source_at(t, provided_function)
                return _tran_companion._as_float(
                    self.cir.i(x, self.epar) + iq + u, self.toolkit)

        ## STAGE PREDICTOR.  A multistep method has no stages, so its analogue
        ## is the classical one: extrapolate the accepted history to `t`.  The
        ## seed it replaces is `x_n`, a whole step behind.
        ## THE CONVERGED POINT IS ONE EVALUATION SESSION: the branch screen
        ## reads `C` there (inside `_newton`) and `jacobian_only` then `q`
        ## and `G` -- the session object both re-enter (`_evalhint`), so a
        ## compiled model computes the three in one fused pass
        conv = _evalhint.Session(('C', 'q', 'G'))
        self._conv_session = conv
        try:
            seed = self._pred_or(x0, t)
            ## THE NEWTON SOLVE IN C, the converged point with it, where it
            ## serves (`_tran_newton_c`); `_newton` where it does not
            fj = None
            r = _tran_newton_c.solve(self, func, t, provided_function, seed,
                                     residual_only)
            if r is None:
                x = self._newton(func, seed, residual=residual_only)
            else:
                x, fj = r
        finally:
            self._conv_session = None
        ## ⚠ AND IT MUST RECORD ITS OWN NODE.  The stage methods get theirs for
        ## free next to `_rk_Y`; this path has no stages, so without this line
        ## the history never reaches two entries and the predictor declines
        ## every step in silence -- measured, 200 calls and 0 predictions.
        self._pred_pending = (t, ())
        ## The source term does not enter `J`, and `jacobian_only` returns
        ## `f = None` by design, so the reduced evaluation stays correct with
        ## `provided_function` folded into `func` above (F4).
        f, J = fj if fj is not None else jacobian_only(x)
        return x, None, J, f
    
    ## `analytical_eh` is not an argument (F8): passing it raises TypeError.
    ## History: `doc/transient_history.md`, `Transient.solve`.
    def solve(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None, fixed_timestep=False, coupled_lte=False):
        """Integrate the circuit from 0 to `tend`; returns the run's
        `CircuitResult` over time (also `self.result`), the start point
        included and `statistics` attached.

        `x0` is the start state; None means the DC operating point, or with
        `uic=True` the `ic` vector (zeros where none is given).  `timestep`
        sets the opening-step scale (`firststep`) and, with
        `fixed_timestep=True`, the uniform grid the run keeps; otherwise the
        step is LTE-controlled and capped by `timestep_max`.
        `provided_function(t)` is an extra source term added to `u(t)`.
        `coupled_lte=True` runs Fang's coupled `(x, h)` step (`fang_timestep`)
        instead of the predict/accept loop.  The steps themselves are
        `_solve`'s.
        """
        ## (the caller may have changed the circuit since the last run)
        self._memo_clear()
        ## Stage 2a: hold BLAS to one thread for the whole run.  It wraps the whole
        ## transient rather than just the linear solve because the win is not in the
        ## solve -- that is ~2% of runtime, so even an infinite speedup there could
        ## not produce the measured 1.72x.  The cost is thread-pool overhead spread
        ## across the many small numpy operations in assembly, and that is only
        ## avoided by setting the limit once, outside the loop.
        with _single_threaded_blas():
            return self._solve(refnode, tend, x0, timestep, provided_function,
                               fixed_timestep, coupled_lte)

    def _finish_result(self, X, timelist, t_start):
        """The run's `CircuitResult`.

        ⚠ The t=0 point IS part of the result -- SPICE convention, and what
        the JAX backend does.  `X[0]` is the operating point (or the uic
        vector) the run worked to compute; dropping it makes every
        index-aligned backend comparison off by one and leaves
        resample_uniform unable to reproduce the initial value (F12(a)).

        STAGE 10.2 -- resampled onto a uniform grid if one was asked for,
        after the run rather than inside it, deliberately: the adaptive grid
        is what the error control is defined on, so the solver keeps
        choosing its own steps and only the REPORTED points change.  The
        statistics are reachable from the result, not only from the
        analysis (Stage 6(c), F13), so a caller who kept only the waveform
        can still ask what produced it.

        History: `doc/transient_history.md`, `Transient._finish_result`."""
        X = self.toolkit.array(X).T
        timelist = self.toolkit.array([0.0] + timelist)
        self.statistics.total_seconds = time.perf_counter() - t_start
        self.result = CircuitResult(self.cir, x=X, xdot=None,
                                    sweep_values=timelist,
                                    sweep_label='time', sweep_unit='s')
        outputstep = self.par.outputstep
        if outputstep is not None:
            _grid, _Xg = resample_uniform(self.result.sweep_values,
                                          self.result.x, step=outputstep)
            self.result = CircuitResult(self.cir, x=_Xg, xdot=None,
                                        sweep_values=_grid,
                                        sweep_label='time', sweep_unit='s')
        self.result.statistics = self.statistics
        return self.result

    @staticmethod
    def _is_stage_family(integ):
        """Whether `integ` steps as a stage family (`_StageSteps`): a
        Runge-Kutta method, or a Nordsieck GLM -- not a Runge-Kutta method,
        but self-starting, keeping no charge ring and delivering its own
        estimate through `_rk_est`."""
        from pycircuit.circuit.integrator import RungeKuttaIntegrator
        return (isinstance(integ, RungeKuttaIntegrator)
                or getattr(integ, 'is_multivalue', lambda: False)())

    def _run_max_step(self, tend, timestep, fixed_timestep):
        """The run's clamp on how large an ACCEPTED step may grow, and the
        delay elements' cap on it (with its warning under a fixed grid)."""
        ## The clamp on how large an ACCEPTED step may grow (`timestep_max`).
        ## Decoupled from `timestep` by owner decision (see the Parameter):
        ## None means SPICE's TMAX default, tend/50.  `timestep` promises no
        ## output density, so a cap below it just clamps the opening step like
        ## any other.
        ## History: `doc/transient_history.md`, `Transient._solve`.
        max_step = self.par.timestep_max
        if max_step is None or max_step <= 0:
            max_step = tend / 50.0

        ## STAGE 8(d) -- a delay element caps the step, per line.
        ##
        ## `TLine` interpolates its history at `t - TD`; with `dt` comparable to
        ## `TD` there is nothing to interpolate between, and the delay simply comes
        ## out wrong (4x TD at twice the cap under `fixed_timestep`).  The adaptive
        ## controller usually rescues it, so it only bites the configuration
        ## nobody checks.
        ##
        ## The cap is asked of the ELEMENTS rather than hard-coded here, so a future
        ## delay element gets it by implementing one method.
        ## History: `doc/transient_history.md`, `Transient._solve`.
        element_cap = self.cir.max_timestep() if hasattr(self.cir, 'max_timestep') else None
        if element_cap is not None:
            ## Under `fixed_timestep` the caller has taken the grid into their own
            ## hands, so silently substituting a finer one would be the wrong kind
            ## of help -- but running on regardless is worse, because the error is a
            ## WRONG DELAY and nothing else reports it.  Warn and obey, which is what
            ## the force-accept and non-convergence paths already do.  The
            ## comparison is against the GRID, not the cap: since timestep_max
            ## decoupled from timestep, the cap says nothing about how coarse
            ## the caller's fixed grid is.
            if fixed_timestep:
                if timestep > element_cap:
                    warn(
                        'transient: fixed_timestep=%g exceeds the %g s cap a delay '
                        'element needs (TD/2). The propagation delay will come out too '
                        'long -- measured 4x at twice the cap -- and nothing else will '
                        'report it. Use a timestep <= %g, or drop fixed_timestep.'
                        % (timestep, element_cap, element_cap), AccuracyWarning)
            elif element_cap < max_step:
                max_step = element_cap
        return max_step

    def _step_family(self, run, coupled_lte):
        """The run's step family (`_LMMSteps`, `_StageSteps`, or
        `_CoupledSteps` for `coupled_lte`) -- see the note at its use in
        `_solve`."""
        if coupled_lte:
            return _CoupledSteps(self, run)
        if self._is_stage_family(self.base_integrator):
            return _StageSteps(self, run)
        return _LMMSteps(self, run)

    def _solve(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None, fixed_timestep=False, coupled_lte=False):
        """`solve`'s run, inside its single-threaded BLAS: the start state
        (`_open_run`), then ONE stepping loop (`_SteppingLoop`) over the
        method's step family (`_LMMSteps`, `_StageSteps`, or `_CoupledSteps`
        for `coupled_lte`), the breakpoints and state events landed on the
        way, then `_finish_result`."""
        self._reset_pcnr_counts()
        if coupled_lte:
            self._refuse_coupled_on_stage()
        X = self._open_run(refnode, x0, provided_function, coupled_lte)

        ## Stage 6(c).  Created per run, so a second `solve()` reports its own
        ## numbers rather than the sum of every run on this object.
        self.statistics = TransientStatistics()
        ## the branch check reports (and may fail) once PER RUN, and its
        ## collapse scale is this run's (`_branch_screen`)
        self._branch_error = None
        self._branch_warned = False
        self._branch_cmax = 0.0
        _t_run_start = time.perf_counter()
        max_step = self._run_max_step(tend, timestep, fixed_timestep)

        ## SOLUTION-flavoured, not residual-flavoured.  This vector is used by the
        ## step controller as a tolerance on `lte = J^-1 * Eg`, which carries the
        ## units of the solution vector x -- volts on node rows, amps on branch rows.
        ## `_newton` needs the other flavour, because there the tolerance applies to
        ## the residual f (KCL currents at nodes), and it builds both separately as
        ## `abstol`/`xtol`.  This is the `xtol` one; `_newton`'s `abstol` here would
        ## apply iabstol (1 pA) as a *voltage* tolerance to every node.
        ##
        ## It reads `lte_vabstol`/`lte_iabstol`, NOT `vabstol`/`iabstol`, so the
        ## controller's knob moves without silently moving Newton's convergence
        ## criterion with it (decision 0.3a).
        ## History: `doc/transient_history.md`, `Transient._solve`.
        abstol = self._lte_abstol_vector()

        ## THE STEP FAMILY -- the only thing about stepping that depends on the
        ## integrator (see `_LMMSteps`, `_StageSteps`, `_CoupledSteps`): how
        ## one step is taken and judged.  Everything else is the same loop
        ## (`_SteppingLoop`) for every method.  A Nordsieck GLM is not a
        ## Runge-Kutta method, but it is self-starting, keeps no charge ring
        ## and delivers its own estimate through `_rk_est`, so it is a stage
        ## family.
        ## History: `doc/transient_history.md`, `Transient._solve`.
        from types import SimpleNamespace
        run = SimpleNamespace(tend=float(tend), max_step=max_step,
                              abstol=abstol, fixed=bool(fixed_timestep))
        family = self._step_family(run, coupled_lte)
        loop = _SteppingLoop(self, family, run, X, tend, timestep,
                             provided_function)
        loop.execute()
        self._warn_run_summary(loop.fixed_fallbacks, loop.forced, timestep)
        return self._finish_result(X, loop.timelist, _t_run_start)

    def _open_run(self, refnode, x0, provided_function, coupled_lte):
        """`_solve`'s start state: the elements' state cleared, the reference
        node and the state events set up, `x0` from the initial conditions
        (`uic`) or the operating point, the run's history begun on it and
        the elements told.  Returns the run's state list, `[x0]`."""
        ## STAGE 8(d) -- clear per-analysis element state BEFORE anything seeds it.
        ##
        ## Position matters: after the initial `accept_step(0.0, ...)` this would
        ## wipe the very history that call seeds, and `TLine.G` would see an empty
        ## buffer and stamp the line as a DC SHORT.  Elements that carry state
        ## must be reset before the run seeds them, not after.
        ## History: `doc/transient_history.md`, `Transient._solve`.
        if hasattr(self.cir, 'reset_state'):
            self.cir.reset_state(self.epar)

        ## The DC operating point that seeds a non-uic run knows nothing about
        ## `provided_function`, so with the extra source in the residual the
        ## t=0 state does not satisfy the t->0+ equations and the opening
        ## steps integrate a startup transient the fix itself manufactured
        ## (F4, second-order effect (i)).  Start from uic=True or an explicit
        ## x0 to avoid it -- or thread the term into the seeding DC, which is
        ## the fuller fix if this warning is ever load-bearing.  The source is
        ## also invisible to `next_event`, so a DISCONTINUOUS
        ## provided_function gets neither breakpoint truncation nor an order
        ## drop: it must be smooth.
        if provided_function is not None and x0 is None and not self.par.uic:
            warn(
                'transient: provided_function adds a source the DC operating '
                'point does not see, so the run opens from an inconsistent '
                'state and integrates a spurious startup transient. Pass '
                'uic=True or an explicit x0.', UsageWarning)

        X = []
        self.irefnode=self.cir.get_node_index(refnode)
        n = self.cir.n
        self._init_state_events(n)
        if x0 is None:
            if self.par.uic:
                ## Skip the operating point and start from the stated initial
                ## conditions -- zeros for anything `ic` does not name.
                x0 = self._initial_state(refnode)
            else:
                ## A failed operating point raises rather than silently
                ## becoming a vector of zeros.
                x0 = self._solve_operating_point(refnode)
        x = x0

        ## `ic` without `uic` is a request the operating point overwrites, so
        ## honouring it silently would be a lie in either direction: SPICE uses
        ## `.ic` to CONSTRAIN the operating point and then releases it, which is
        ## a different feature from the one implemented here. Raising says which
        ## one is missing rather than quietly doing neither.
        ## `include_state=False`: an Idt/Idtmod `ic` pins the DC operating
        ## point (LRM), so it is meaningful without uic and exempt here.
        if (self.par.ic or self._descendant_has_ic(self.cir,
                                                   include_state=False)) \
                and not self.par.uic:
            raise ValueError(
                "ic was given without uic=True. This implements SPICE's initial "
                "conditions for the uic case only -- starting values for the "
                "transient. Constraining the operating point with .ic and then "
                "releasing it is a separate feature and is not implemented. "
                "Pass uic=True, or drop ic.\n"
                "(This covers element initial conditions such as L(..., ic=...) "
                "as well as the analysis-level ic dict -- both are starting "
                "values, and both are ignored without uic.)")

        if coupled_lte:
            ## P22: eq (6)'s state-row mask, built once at the seed
            self._lte_state_mask = self._state_row_mask(x0)
        self._begin_run(x, n)

        X.append(copy(x))
        if hasattr(self.cir, 'accept_step'):
            self.cir.accept_step(0.0, X[-1], self.epar)
        return X

    def _warn_run_summary(self, fixed_fallbacks, forced, timestep):
        """One warning a run for each kind of step it could not take as
        asked (the fixed grid's Newton fallbacks, the force-accepts), and
        for its PCNR fallbacks (`_warn_pcnr_summary`)."""
        if fixed_fallbacks:
            t0, h0 = fixed_fallbacks[0]
            warn(
                f'transient: Newton did not converge at {len(fixed_fallbacks)} '
                f'step(s) with the requested fixed timestep {timestep:.6g} s '
                f'-- the first at t={t0:.6g} s, falling back to {h0:.6g} s for '
                'that step (each a quarter of the step it replaced). The '
                'output grid is no longer uniform.', ConvergenceWarning)
        if forced:
            t0, h0, r0 = forced[0]
            warn(
                'transient: local truncation error still above tolerance after '
                f'the rejection budget at {len(forced)} step(s) -- the first at '
                f't={t0:.6g} s (h={h0:.6g} s, {r0} rejections); those steps '
                'were accepted with an order drop. The accepted error is '
                'unbounded -- treat the waveform near those times with '
                'suspicion (statistics: force_accepts).', AccuracyWarning)
        self._warn_pcnr_summary('transient')

    def _attempt_step(self, family, X, t, h, hold, provided_function):
        """One attempt of the family's step, timed."""
        _t0 = time.perf_counter()
        try:
            return family.attempt(X, t, h, hold, provided_function)
        finally:
            self.statistics.solve_seconds += time.perf_counter() - _t0

    def _rescue_step(self, family, X, t, h, hold, provided_function):
        """The last resort at `minstep`: the step re-solved with the
        continuation chain armed (`_continuation_rescue`, read by `_newton`),
        or a `TransientStepError` naming the point.

        One loop, one ladder: every step family reaches it.

        History: `doc/transient_history.md`, `Transient._rescue_step`."""
        self._dt = h
        ## ⚠ WHETHER THIS PATH REACHES A LADDER AT ALL
        ## (`_honours_continuation_rescue`): the error must not claim a rescue
        ## that was never attempted, and a success is a rescue only where one
        ## could run
        ## History: `doc/transient_history.md`, `Transient._rescue_step`.
        honours = self._honours_continuation_rescue()
        self._continuation_rescue = True
        try:
            out = self._attempt_step(family, X, t, h, hold, provided_function)
            if honours:
                self.statistics.gmin_rescues += 1
            return out
        except NoConvergenceError as e:
            if honours:
                raise TransientStepError(
                    'Transient solver failed to converge: timestep shrank below '
                    'minstep=%gs at t=%s, and the gmin/gshunt/pseudo-transient '
                    'continuation could not rescue the point: %s'
                    % (self.par.minstep, t, e)) from e
            raise TransientStepError(
                'Transient solver failed to converge: timestep shrank below '
                'minstep=%gs at t=%s, and this step path carries no '
                'continuation rescue (neither a ladder nor a fallback to one), '
                'so none was attempted: %s' % (self.par.minstep, t, e)) from e
        finally:
            self._continuation_rescue = False



if __name__ == "__main__":
    import doctest
    doctest.testmod()
