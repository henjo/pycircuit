"""Fang's coupled (x, h) stepping (DAC 2013) and its LTE band.  A theme of
`Transient` (see `transient.py`).
"""

import numpy as np

from pycircuit.circuit import pcnr as _pcnr

## The clamp the step controller applies to every accepted step, the force-accept
## path in `solve()` included: one bound, named once.  `stepcontroller` imports
## nothing from this package, so this is import-safe at module level.
## History: `doc/transient_history.md`, `transient.py`.
from pycircuit.circuit.stepcontroller import normalised_error


class _CoupledLTE:
    """Fang's coupled (x, h) stepping (DAC 2013) and its LTE band.  A theme of
    `Transient` (see `transient.py`)."""

    def fang_timestep(self, x_prev, t_prev, h, x_hist,
                      provided_function=None, gamma_min=0.7, gamma_max=3.0,
                      eta=0.15, maxiter=None, hmin=None, max_step=None,
                      hold_h=False, grid_locked=False, method='approx'):
        """One time point of Fang's coupled method from `x_prev` at
        `t_prev`, solving for the state AND the step together from the trial
        step `h` (see `_fang_timestep_inner`, whose parameters these are).
        Returns ``(x, h, iterations, converged)``.  Around it, the
        `sigglobal` running reference is rebuilt from the RETURNED solution
        only, so an unconverged iterate never moves a later tolerance."""
        ## R2 HYGIENE (doc/transient_review_260820.md, refuted-but-latent):
        ## the Newton loop below folds every UNCONVERGED iterate into the
        ## sigglobal running maximum through _lte_tolerance -> _reference.
        ## Measured pollution was <= 7.5e-6 relative across six circuits --
        ## not a live bug -- but the JAX backend's accept-only update is the
        ## right hygiene, so the running max is snapshotted here and, on
        ## exit, rebuilt from the snapshot plus the solution actually
        ## returned: iterates influence tolerances only within this time
        ## point, never beyond it.
        ctrl_probe = getattr(self, 'step_controller', None)
        if getattr(self, '_step_controller_is_auto', False):
            ctrl_probe = None
        if ctrl_probe is None:
            ctrl_probe = getattr(self, '_fang_controller', None)
        snapshot = getattr(ctrl_probe, '_ref_running', None) \
            if ctrl_probe is not None else None
        result = self._fang_timestep_inner(
            x_prev, t_prev, h, x_hist, provided_function, gamma_min,
            gamma_max, eta, maxiter, hmin, max_step, hold_h, grid_locked,
            method)
        ## The inner call may have created the controller; re-resolve it.
        ctrl = getattr(self, '_fang_controller', None) or ctrl_probe
        if ctrl is not None and getattr(ctrl, 'relref', None) != 'pointlocal':
            local = np.abs(np.asarray(result[0], dtype=float))
            ctrl._ref_running = local if snapshot is None \
                else np.maximum(np.asarray(snapshot, dtype=float), local)
        return result

    def _fang_setup(self, method):
        """What `_fang_timestep_inner` needs that does not change within a
        run: the check of `method`, the refusal of an injected controller
        whose law the coupled path does not implement, the LTE controller
        (created on first use, `relref` applied), the PCNR devices and the
        Newton solution tolerance.  `_CoupledSteps` builds it once per run as
        `_fang_run`; a `fang_timestep` call outside a run builds its own."""
        from types import SimpleNamespace

        from pycircuit.circuit.stepcontroller import SolutionLTEController

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
        ## documented feature that does not exist.
        ##
        ## NOT keyed on `lte_gradients`, which would be the obvious test and is
        ## wrong: `q^T` and `d` are implemented and gated but are NOT called on
        ## the shipped path, because sec. 3.4 replaced the eq (12) branch that
        ## used them.  Testing for a method nothing calls would pass controllers
        ## that cannot work and fail ones that can.
        ## History: `doc/transient_history.md`, `Transient._fang_timestep_inner`.
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

        if method == 'bordered':
            raise ValueError(
                "coupled_method='bordered' (Fang eq 12/14) was retired on "
                "2026-09-27: once its double-counted q^T dv0 term was removed "
                "(2026-09-19) it took the same steps as 'approx' to every "
                "printed digit -- measured on a smooth and a pulsed RC, the "
                "same accepted steps, Newton iterations and error. Use "
                "'approx' (the default).")
        if method != 'approx':
            raise ValueError(
                "coupled_method must be 'approx' (Fang sec 3.4), not %r"
                % (method,))

        ## PCNR ON THE COUPLED PATH: `pcnr=True` is honoured here too.
        ## History: `doc/transient_history.md`, `Transient._fang_timestep_inner`.
        junctions = _pcnr.pcnr_devices(self.cir) if self.par.pcnr else []

        ## The increment flavour: `dx0` is a solution update, not a residual.
        xtol = self._newton_xtol_vector()
        return SimpleNamespace(method=method, ctrl=ctrl, junctions=junctions,
                               xtol=xtol)

    def _fang_timestep_inner(self, x_prev, t_prev, h, x_hist,
                             provided_function, gamma_min, gamma_max,
                             eta, maxiter, hmin, max_step,
                             hold_h, grid_locked, method):
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
        toolkit = self.toolkit
        n = self.cir.n
        irefnode = self.irefnode
        maxiter = self.par.maxiter if maxiter is None else maxiter
        hmin = self.par.minstep if hmin is None else hmin
        if max_step is None:
            max_step = self.par.timestep_max
            if max_step is None or max_step <= 0:
                max_step = float('inf')

        ## The once-per-run part (the controller, the PCNR devices, the Newton
        ## solution tolerance, the checks): `_CoupledSteps` builds it when the
        ## run starts; a direct call outside a run builds its own.
        run = getattr(self, '_fang_run', None)
        if run is None or run.method != method:
            run = self._fang_setup(method)
        ctrl, junctions, xtol = run.ctrl, run.junctions, run.xtol

        x = toolkit.array(x_prev, dtype=float).copy()

        ## `v_lim` is per-time-point state, seeded from the incoming solution and
        ## carried across the iterations below.
        v_lim = _pcnr.v_lim_init(junctions, x)

        reltol = self.par.reltol

        ## Eq (16) bounds the step change BETWEEN ITERATIONS, and iterating it is
        ## how the step size collapses inside a single time point: 0.85 per
        ## iteration over `maxiter` iterations is seven decades.  Measured on the
        ## charge pump, `h` reached 8.75e-15 s at t=1.1e-5 before the solve gave
        ## up.  So the TOTAL excursion within one time point is bounded too, by
        ## the same window the standard controller allows for one step.
        from pycircuit.circuit.stepcontroller import MAX_GROWTH_RATIO, MIN_SHRINK_RATIO

        h_entry = h
        h_floor = max(hmin, h_entry * MIN_SHRINK_RATIO)
        h_ceil = min(max_step, h_entry * MAX_GROWTH_RATIO)

        for it in range(maxiter):
            ## --- STAGE 1: the ordinary N system, at the current step size.
            self._dt = h
            t = t_prev + h
            if junctions:
                f, J, g_lim = self._residual_and_jacobian_pcnr(
                    x, v_lim, t, junctions, provided_function)
            else:
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
            if junctions:
                ## PCNR's CORRECT phase.  No `cir.limit` anywhere: each device
                ## limits ONLY the unknown it owns, so one device's limiter
                ## cannot disturb another's.  `dx0` is deliberately NOT shortened
                ## -- the MNA update is taken in full, which is the whole point.
                x_stage1 = x + dx0
                x_stage1[irefnode] = 0.0
                dx_lim = _pcnr.dx_lim_of(junctions, g_lim, dx0)
                v_stage1 = _pcnr.refine(junctions, v_lim, v_lim + dx_lim,
                                        self.epar, x_old=x)
            else:
                x_stage1 = self.cir.limit(x + dx0, x, self.epar)
                ## The limiter may shorten the step, so the convergence test must
                ## use what was actually taken, not what was asked for.
                dx0 = x_stage1 - x
                v_stage1 = v_lim

            ## --- STAGE 2: has the step size earned any attention?
            h_hist = [hh for hh in (self._dt_last,
                                    self._dt_last2)
                      if hh is not None]
            etol = self._lte_tolerance(ctrl, x_stage1, x_prev, h_hist)

            eps_ok, err = self._lte_in_band(ctrl, x_stage1, x_hist, h_hist, h,
                                            etol, gamma_min, gamma_max)

            converged_x = bool(toolkit.alltrue(
                abs(dx0) < reltol * abs(x_stage1) + xtol))
            if junctions:
                ## BOTH residuals, for the reason recorded in
                ## `_solve_timestep_pcnr`: converging on `dx_mna` alone can return
                ## with `v_lim != e_a - e_b`, i.e. the diode evaluated at a voltage
                ## that is not the node voltage, so the vector is not a solution of
                ## the circuit at all -- and the LTE built from it then reads low.
                converged_x = converged_x and _pcnr.lim_converged(
                    g_lim, v_stage1, reltol, self.par.vabstol)

            ## `hold_h` -- the step size is IMPOSED, not free.  A step truncated
            ## onto a breakpoint or onto `tend` has its size decided by where it
            ## must land, so there is nothing for the coupled system to SOLVE
            ## (solving for its own `h` walks straight off the edge again).
            ##
            ## BUT "DO NOT SOLVE FOR h" IS NOT "DO NOT CHECK THE ERROR": a held
            ## step accepted blind has a truncation error governed by nothing.
            ## A held step whose error is over the band is reported so the
            ## caller can shrink and retry -- UNLESS the grid is locked.
            ##
            ## `fixed_timestep` is the caller stating that the output points are
            ## theirs, so shrinking is not an option available to us: the honest
            ## response to an over-tolerance step on a locked grid is to take it
            ## and let the run's accuracy be what the caller asked for, exactly
            ## as the standard path does.
            ## History: `doc/transient_history.md`, `Transient._fang_timestep_inner`.
            if hold_h and not grid_locked and not eps_ok and err > gamma_max:
                return x_stage1, h, it + 1, False

            if hold_h or eps_ok or not h_hist or len(x_hist) < 2:
                ## The LTE condition holds (or cannot be evaluated yet, on the
                ## opening steps).  Nothing to solve for `h`; finish on the
                ## circuit equations alone.
                x, v_lim = x_stage1, v_stage1
                if converged_x:
                    return x, h, it + 1, True
                continue

            ## --- The LTE condition failed, so the step size must move too.
            ##
            ## SEC. 3.4's APPROXIMATE NEWTON, NOT EQ (12).  Eq (12) recovers
            ## `dh` from eq (14), whose denominator is `q^T dxh + d`.  Those two
            ## terms are the solution's sensitivity to the step size and the
            ## extrapolation's slope, and BOTH are approximately `dv/dt`: their
            ## difference is the truncation error's derivative, which is tiny
            ## by construction.  Eq (12) computes a small quantity as the
            ## difference of two large ones (on a driven RC, three digits lost
            ## and the SIGN of `dh` decided by the cancellation).
            ##
            ## Eq (17) gets the new step from the error RATIO instead, which
            ## involves no cancellation at all.  `step_for_error_ratio` inverts
            ## the node polynomial rather than applying the (tau/eps)^(1/(n+1))
            ## power law, because that law only holds while h >> h_last -- see
            ## `extrapolation_error_weight`.
            ## History: `doc/transient_history.md`, `Transient._fang_timestep_inner`.
            from pycircuit.circuit._lte_kernels import step_for_error_ratio

            target = self._band_centre(ctrl, gamma_min, gamma_max)

            ## (`coupled_method='bordered'` is retired: see the check above.)
            ## History: `doc/transient_history.md`, `Transient._fang_timestep_inner`.
            ratio = target / max(err, 1e-300)
            h_new = step_for_error_ratio(h, h_hist, ratio,
                                         1.0 - eta, 1.0 + eta)

            ## WHAT THE STEP WANTS, ignoring every clamp.  Saturation has to be
            ## measured against this, not against the clamped result: once `h` is
            ## pinned at a bound the clamped `dh` is exactly 0.0, which is
            ## indistinguishable from "the step size has stopped moving" -- the
            ## definition of converged in eq (16) -- when in fact it stopped
            ## because it hit a wall.
            ## History: `doc/transient_history.md`, `Transient._fang_timestep_inner`.
            h_want = step_for_error_ratio(h, h_hist, ratio, 1e-6, 1e6)
            h_new = min(max(h_new, h_floor), h_ceil)
            dh = h_new - h

            ## DID THE CORRECTION SATURATE?  Eq (16), `|dh| <= eta*h`, is a
            ## convergence criterion meaning "the step size has stopped moving".
            ## A correction pinned AT the limiter has not stopped moving -- it
            ## was cut off -- and testing it with `<=` makes the two
            ## indistinguishable, because a clamped `dh` equals `eta*h` exactly.
            ##
            ## Only a thwarted SHRINK counts. A step that wants to grow and is
            ## held at the cap is the normal state of every adaptive controller
            ## -- growth is bounded by zero stability, not by the error -- and
            ## treating that as unconverged drove `h` to 9.5e-16 at t = 1e-9,
            ## because the opening steps always want to grow faster than allowed.
            ## History: `doc/transient_history.md`, `Transient._fang_timestep_inner`.
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
                ## `v_lim` tracks the branch voltage, so the stage-4 correction
                ## has to move it too -- the loop can return immediately after
                ## this, and a `v_lim` left at its stage-1 value would be exactly
                ## the `v_lim != e_a - e_b` inconsistency the check above exists
                ## to catch.
                v_lim = np.array(
                    [v_stage1[k] + (dxh[ra] - dxh[rb]) * dh
                     for k, (ra, rb) in enumerate(
                         _pcnr.flat_probes(junctions))])
            else:
                x, v_lim = x_stage1, v_stage1
            h = h_new

            ## `not saturated`, and a strict test: a step still moving at the
            ## limiter must keep iterating, and if the time point's own bound
            ## (`h_floor`) is what stopped it, the caller shrinks and retries.
            if converged_x and not saturated and abs(dh) < eta * h:
                return x, h, it + 1, True

        return x, h, maxiter, False

    def _coupled_band(self):
        """The (gamma_min, gamma_max, eta) the coupled path runs with.

        Maps the 'auto' sentinel to Fang's sec. 4.1 values; any explicit
        number (including the documented 0.0 / 1.0 / None) passes through
        verbatim -- the property F5 exists to restore
        (doc/transient_review_260820.md).  The standard path's mapping lives
        in StepController.set_lte_band, whose own defaults ARE that path's
        'auto' resolution; the two sites partition cleanly because
        fang_timestep takes the band as kwargs and never calls set_lte_band.
        """
        gm, gx, eta = (self.par.lte_gamma_min, self.par.lte_gamma_max,
                       self.par.lte_eta)
        return (0.7 if gm == 'auto' else gm,
                3.0 if gx == 'auto' else gx,
                0.15 if eta == 'auto' else eta)

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
        etol = ctrl.tolerance(ref, self.par.reltol, self._lte_abstol_vector(),
                              self.LTERATIO)
        ## P22: eq (6) over the STATE rows only -- an infinite tolerance on
        ## algebraic rows removes them from the band test and the controlling-
        ## node argmax through this one mechanism.  See _state_row_mask for
        ## the derivation.
        ## History: `doc/transient_history.md`, `Transient._lte_tolerance`.
        ## 1e30, not inf: lte_gradients differentiates 1/etol terms, and an
        ## inf there turns a masked row's gradient into 0*inf = NaN.
        mask = getattr(self, '_lte_state_mask', None)
        if mask is not None:
            etol = np.where(mask, etol, 1e30)
        return etol

    def _state_row_mask(self, x_ref):
        """True where the unknown participates in ANY charge -- P22.

        Fang's eq (6) is a truncation estimate for the DIFFERENTIATED
        variables of the DAE; rows whose unknown appears in no charge (zero
        row AND zero column of C) are algebraic -- Lagrange-multiplier
        currents of voltage-defining branches, purely resistive node
        voltages -- and are slaved to the states through the Jacobian, not
        integrated.  Measuring eq (6) on them measures conventions, not
        truncation: the rectifier's source-current row carries the diode's
        dq/dt through KCL, its accepted value holds the OLD grid's
        derivative convention, the re-solve computes the NEW grid's, and
        the deviation floor (measured 2.5e-6 A against etol 3.6e-7) is
        h-independent -- the band can never be satisfied and the coupled
        run livelocks (the Gear-2 rectifier trace behind P22; the TLine
        campaign's from-zero kinks were the same class on port rows).

        Structural, from C at the seed point; a nonlinear charge whose C
        row vanishes AT the seed but not elsewhere would be misclassified
        -- accepted for now, recorded here, revisit if a circuit shows it.
        """
        C = self.toolkit.toMatrix(self.cir.C(x_ref, self.epar))
        ## `toMatrix` can hand back a complex matrix (imaginary part exactly
        ## zero); `abs` of it, not a float cast, which warned on every run
        Ca = np.abs(np.asarray(C))
        mask = (Ca.sum(axis=0) + Ca.sum(axis=1)) > 0.0
        return mask

    def _lte_in_band(self, ctrl, x_curr, x_hist, h_hist, h, etol,
                     gamma_min, gamma_max):
        """``(condition_holds, normalised_error)`` for eq (15): the solution
        controller's estimator (`SolutionLTEController.solution_deviation`;
        its degree cap, 2, never binds here: `h_hist` holds at most the two
        past steps)."""
        from pycircuit.circuit.stepcontroller import SolutionLTEController

        if not h_hist or len(x_hist) < 2:
            return True, 0.0
        order = getattr(self.active_integrator, 'ORDER', 1)
        lte, _degree = SolutionLTEController.solution_deviation(
            x_curr, x_hist, h_hist, h, order)
        if lte is None:
            return True, 0.0
        err = float(np.max(normalised_error(lte, etol)))
        return (gamma_min <= err <= gamma_max), err

    def _refuse_coupled_on_stage(self):
        """`coupled_lte=True` on a stage method or a GLM, refused by name
        before any work.

        ⚠ FANG'S COUPLED PATH IS BUILT ON A LINEAR MULTISTEP COMPANION: it
        solves the step from eq (6), a solution-space LTE over the step
        history.  A stage method or a GLM judges its step by its own
        embedded estimate.  Refused here, by name, before any work.
        ⚠ AND REFUSED ON MEASUREMENT: around a stage step the coupled loop
        does not keep Fang's no-rejection property, and on a stiff circuit
        its error is set by no tolerance, because eq (6) needs accepted
        history and so cannot judge the opening steps an embedded estimate
        judges from the first.  Script:
        `benchmarks/transient_review/stage12c_fang_stage_methods.py`.
        History: `doc/transient_history.md`, `Transient._solve`.
        """
        _integ = self._get_integrator()
        if self._is_stage_family(_integ):
            raise NotImplementedError(
                "coupled_lte=True (Fang's coupled (x, h) stepping) is "
                "built on a linear multistep companion (gear, trap, "
                "euler, theta), not on %s: it solves the step size from "
                "eq (6), a solution-space LTE over the step history, "
                "where a Runge-Kutta stage method or a GLM judges its "
                "step by its own embedded estimate. Measured on a "
                "prototype: around a stage step it re-solves more often "
                "than the standard run rejects, and on a stiff circuit "
                "its error does not follow the tolerance (eq (6) cannot "
                "judge the opening steps). Run coupled_lte=False (the "
                "embedded estimate drives the adaptive step), or choose "
                "a multistep integrator." % type(_integ).__name__)
