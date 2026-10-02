"""The run's history (begin, push, roll, freeze, the periodic shifts) and what
the shooting reads off a step (`step_lte`, `residual_dh`, `residual_dT`).  A
theme of `Transient` (see `transient.py`).
"""

import numpy as np

from pycircuit.circuit._limiting import limit_sync, stateful_limiters

## The clamp the step controller applies to every accepted step, the force-accept
## path in `solve()` included: one bound, named once.  `stepcontroller` imports
## nothing from this package, so this is import-safe at module level.
## History: `doc/transient_history.md`, `transient.py`.


class _RunHistory:
    """The run's history (begin, push, roll, freeze, the periodic shifts) and
    what the shooting reads off a step (`step_lte`, `residual_dh`,
    `residual_dT`).  A theme of `Transient` (see `transient.py`)."""

    def _collect_periodic_states(self):
        """Poll the circuit's periodic-state declarations, once per solve.

        Phase 2 of idtmod.md (sec. 5.2).  Polled AFTER parameters are
        resolved (`modulus` may be a late-bound expression) and verified
        against the one contract the gauge shift rests on: ``q[row] ==
        x[row]``, i.e. the declared row's C row is a unit diagonal, so
        shifting the charge ring by ``n*modulus`` IS shifting the state
        ring.  A non-conforming declaration fails loudly here rather than
        corrupting the LTE estimate silently mid-run.
        """
        if not hasattr(self.cir, 'periodic_states'):
            return []
        declared = self.cir.periodic_states()
        if not declared:
            return []
        C = self.cir.C(self.toolkit.zeros(self.cir.n), self.epar)
        rows = []
        for row, m, o in declared:
            c_row = np.asarray(C[int(row)], dtype=float)
            expected = np.zeros_like(c_row)
            expected[int(row)] = 1.0
            if not np.allclose(c_row, expected):
                raise ValueError(
                    "periodic_states row %d violates its contract: the gauge "
                    "shift requires q[row] == x[row] (a unit C diagonal on "
                    "that row), but the assembled C row is %r. The declaring "
                    "element cannot be wrapped this way." % (int(row), c_row))
            rows.append((int(row), float(m), float(o)))
        return rows

    ## -----------------------------------------------------------------------
    ## THE SHOOTING'S USE OF THIS CLASS.  `PSS` drives `solve_timestep` on a
    ## frozen grid and never calls `_solve`.  What it calls, in order (its
    ## side: `_InnerTransient`, `_PeriodWalks`):
    ##
    ##   built       `PSS._new_transient`: the PSS's settings (tolerances,
    ##               `epar`, `relref`, solvers, PCNR, the method), then
    ##               `_freeze_grid()` unless the grid is its own (`lte_grid`,
    ##               `tstab`)
    ##   per period  `_begin_run(x0)` -- or, for the pair map on a solved
    ##               two-point history, `_begin_run_on_history`
    ##   per step    `_dt = h`; `solve_timestep`; the reads that need the
    ##               previous step's history (`step_lte`, `residual_dh`,
    ##               `residual_dT`) BEFORE `_roll_history(x, h)` moves it;
    ##               what the step left through `last_step` (`LastStep`)
    ##
    ## Left out on purpose, so the period map is a function of `x0` alone:
    ## `cir.accept_step` (element state -- an Idtmod's wrap prediction; a
    ## TLine's history, which is why the PSS refuses hidden state), the
    ## statistics, the family's `after_accept`, the rescue ladder,
    ## breakpoints (the shooting lands source edges on its own grid,
    ## `PSS.event_grid`).  The order drop after a landing it TAKES, as the
    ## loop does (`_InnerTransient.solve_timestep`, keyed to the edge's
    ## node), since 2026-09-28 -- a TRADE: on a smooth state the shooting was
    ## more accurate without it (gear's error 1/3, trap's 1/50 .. 1/100;
    ## `benchmarks/pss_transient_boundary.py` V1), on a STIFF one it is what
    ## damps the corner -- without it trap's current was O(1) off after
    ## every landed edge (`benchmarks/landing_order_drop.py`).  The default
    ## method (radau) keeps no history across an edge.
    ## The gauge shift runs in `_roll_history` without the caller's window,
    ## and the closure folds the whole moduli between (V2).  A GLM is never
    ## handed its Nordsieck vector: `_begin_run` empties both slots, the
    ## first step starts afresh at the period's start, and the walk reads
    ## each step's records back (`last_step`).
    ## History: `doc/pss_log_260902.md`, 2026-09-28 (the interface review).
    ## -----------------------------------------------------------------------

    def _begin_run(self, x, n):
        """Reset every piece of PER-RUN integrator state and seed the rings.

        `_solve` calls it, and so does `PSS` -- it re-integrates one period
        from a fresh state on every shooting iteration, so "begin a run"
        happens many times per analysis there.

        The charge history is rebuilt per run, so the STEP history must be:
        without that, a second `solve()` on the same object starts with a
        stale `_dt_last2` while `_qlast[2]` is the freshly seeded initial
        charge, breaking the invariant 4g(b) relies on -- that
        `h_last2 is not None` exactly when `q_last[2]` is a real past point.

        History: `doc/transient_history.md`, `Transient._begin_run`.
        """
        self.base_integrator = self._get_integrator()
        hist_len = max(2, self.base_integrator.get_required_history())
        ## the stateful limiters AT the run's first point: what the devices
        ## are read at below, and a PSS period's start must not depend on
        ## where the previous period left them (see `solve_timestep`)
        self._stateful_lims = stateful_limiters(self.cir)
        if self._stateful_lims:
            limit_sync(self.cir, x, self.epar, self._stateful_lims)
        q0 = self.cir.q(x, self.epar)
        self._qlast = self.toolkit.array([q0 for _ in range(hist_len)])
        self._iqlast = self.toolkit.zeros((hist_len, n))
        self._pred_reset()
        ## (the step's source memo: none between steps -- `_source_at`)
        self._u_memo = None
        ## the Nordsieck slots are per-RUN too: a second solve() on the same
        ## object must not read the first run's vector -- nor its stage
        ## predictor the first run's last step.  (The entry slot is keyed by
        ## the time its step STARTED: a PSS period starts at the time the
        ## previous traversal's first step did, and a slot left behind would
        ## be picked up in preference to a fresh startup.)
        ## History: `doc/shooting_history.md`, `_PeriodWalks._glm_period_blocks`.
        self._glm_Q = None
        self._glm_Q_at_entry = None
        self._glm_prev = None
        ## ⚠ A ZERO `iq` RING IS ONLY SAFE FOR A METHOD OPENED BY EULER, which
        ## reads no past current.  A method that refuses that opener reads it on
        ## step one, and zero is wrong there -- measured at a full order of
        ## accuracy (see `Integrator.needs_consistent_iq0`).  Asked of the
        ## method, so every existing integrator keeps the zero ring exactly.
        if self.base_integrator.needs_consistent_iq0():
            t0 = float(getattr(self.epar, 't', 0.0) or 0.0)
            iq0 = -(np.asarray(self.cir.i(x, self.epar), dtype=float)
                    + np.asarray(self.cir.u(t0, analysis='tran'), dtype=float))
            self._iqlast = self.toolkit.array([iq0 for _ in range(hist_len)])
        self._dt_last = None
        self._dt_last2 = None
        ## Excursion-check running maxima are per-run state too.
        self._dv_run_v = 0.0
        self._dv_run_i = 0.0
        ## Phase-2 gauge shift (idtmod.md 5.2): polled after update_iparv so
        ## late-bound moduli are resolved; static for the run.
        self._periodic_rows = self._collect_periodic_states()
        self._is_first_step = True
        self._no_history = True
        self._transform_pcnr_warned = False
        ## The measurement probe's running reference is per-run state too:
        ## carrying one run's signal maximum into the next would make the
        ## relative floor depend on what ran before it.
        self._lte_probe = None

    def _push_history(self, x, X=None):
        """Push one ACCEPTED point onto the integrator's ring buffers.

        The one ring push for every accept site -- `_solve`, and `PSS`, which
        imposes its own grid and so cannot use the loop -- so the
        trailing-window gauge shift below cannot be applied to one and
        forgotten in another.

        `X` is the solution window the gauge shift also has to rewrap; a
        caller that keeps no waveform passes none, and then only the state
        and the `_qlast` ring are shifted.

        History: `doc/transient_history.md`, `Transient._push_history`.
        """
        self._iqlast = self.toolkit.concatenate(
            (self.toolkit.array([self._iq]), self._iqlast))[:-1]
        self._qlast = self.toolkit.concatenate(
            (self.toolkit.array([self._q_at(x)]), self._qlast))[:-1]
        self._pred_promote(x)
        ## AFTER the ring push, so the newest ring entry shares the old gauge
        ## with its elders when the increment lands on all of them.
        if self._periodic_rows:
            self._apply_periodic_shifts(x, [] if X is None else X)

    def _roll_history(self, x, h, X=None):
        """Roll the run past one ACCEPTED step of `h` ending at `x`: the ring
        push (`_push_history`, the gauge shift with it) and the record the
        next step reads -- `_dt`, the `_dt_last`/`_dt_last2` ring, and the
        flags that say the run is no longer opening.

        The one accept bookkeeping for the loop (`_solve`) and the shooting
        (`PSS.solve_timestep`).  What only the loop does at an accept --
        `cir.accept_step`, the statistics, the family's `after_accept`, the
        order drop re-armed after a landing -- stays in the loop, and the
        shooting leaves it out on purpose (the period map must be a function
        of `x0` alone)."""
        self._push_history(x, X)
        self._dt = h
        self._dt_last2 = self._dt_last
        self._dt_last = h
        self._is_first_step = False
        self._no_history = False

    def _begin_run_on_history(self, x0, x_m1, h, h_prev=None):
        """Open a run ON a solved two-point history, full-width `x0` and the
        point `x_m1` one step `h_prev` before it (default `h`), instead of on
        a seed -- the shooting's pair map (`PSS._install_history`).

        `_begin_run(x_m1)` opens the rings on the earlier point and the
        charge half of `_push_history` puts `q(x0)` in front of it (no `_iq`
        roll: no step has been solved, and only a `b = 0` companion, which
        never reads one, may open this way -- the caller refuses the rest).
        The flags then say what is true: a step of `h` has been taken and
        the run is not opening, so nothing drops order.  `_dt_last2` stays
        None on purpose: the THIRD charge in the ring is `q(x_m1)` repeated,
        so the LTE estimator's opening reading stays unsound and is
        discarded.  ⚠ `_dt_last` is the step that PRODUCED `x0`, on a
        periodic grid the period's LAST one (`h_prev`), not its first."""
        self._begin_run(x_m1, self.cir.n)
        self._dt = h
        q0 = self.cir.q(x0, self.epar)
        self._qlast = self.toolkit.concatenate(
            (self.toolkit.array([q0]), self._qlast))[:-1]
        self._q_cache = None
        self._is_first_step = False
        self._no_history = False
        self._dt_last = h if h_prev is None else h_prev
        self._dt_last2 = None

    def _freeze_grid(self):
        """Configure this transient for traversals of a FROZEN grid -- the
        shooting's (`PSS._new_transient`; its users: the period walks, the
        replay, `warping_estimate`'s fixed-step run and the transient PAC
        reaches through `pss._C_at`/`_G_at`/`_k_at`):

        - Gear-2's shrink drop to Euler off (`Gear2Integrator.shrink_guard`):
          a frozen grid has no stalled estimate, and the drop costs a
          two-alpha step the shooting adjoint cannot transpose; growth stays
          guarded;
        - the damped Newton as the last resort (`_damped_last_resort`), in
          place of the step reduction an adaptive run makes: a frozen grid
          cannot shrink (an owner decision; see `_rk_step_coupled`).

        A method, not a Parameter: nothing a caller of `Transient` sets."""
        try:
            self.par.integrator.shrink_guard = False
        except Exception:                                      # noqa: BLE001
            pass
        self._damped_last_resort = True

    def _apply_periodic_shifts(self, x, X):
        """Rewrap periodic states after an ACCEPTED step -- the gauge shift.

        Subtracts ``n*modulus`` from the accepted ``x`` and from every live
        history entry of that row: the trailing solution window (the last
        <=3 entries of ``X`` -- all any integrator or controller reads) and
        the whole ``_qlast`` ring.  ``_iqlast`` is derivative-domain and
        invariant under a constant shift.  A uniform translation of the
        entire read window is invisible to the multistep formulas and the
        divided-difference LTE estimates (idtmod.md sec. 3.2), so this is
        exact -- no event, no restart, no order drop.

        Never called inside a Newton solve: within-step values may exceed
        the window by up to the step's excursion, and the element's own
        output wrap covers that.

        ``_q_cache`` invalidation is load-bearing: `_q_at` memoises ``(x,
        q)`` identity-first, and the in-place shift of ``x`` would keep the
        identity while the cached charge goes stale -- serving it corrupts
        the LTE estimate silently (the exact failure class `_q_at`'s own
        docstring warns about).

        Older entries of ``X`` keep whatever gauge they were recorded in;
        only the element's branch OUTPUT is observable, and it is wrapped in
        ``i()`` at every point.  The private state row's recorded waveform
        is gauge-dependent and not user-meaningful (idtmod.md sec. 5.2).
        """
        shifted = False
        for row, m, o in self._periodic_rows:
            n_wraps = int(np.floor((float(x[row]) - o) / m))
            if n_wraps != 0:
                d = n_wraps * m
                x[row] -= d
                for k in range(1, min(3, len(X)) + 1):
                    X[-k][row] -= d
                self._qlast[:, row] -= d
                ## the predictor's nodes are a live history like the others;
                ## shifting them ALL by the same `d` keeps the trajectory it
                ## fits continuous, which is the only thing a polynomial
                ## through them reads
                for _k in range(len(getattr(self, '_pred_hist', ()) or ())):
                    self._pred_hist[_k][1][row] -= d
                ## ⚠ AND A GLM'S READ WINDOW IS ITS NORDSIECK VECTOR: its
                ## charge level (block 0) is `q(x)` and moves with the state;
                ## the scaled derivatives do not see a constant.  Left out,
                ## the next step continued the UNwrapped charge: measured, the
                ## state row sat 2 moduli further out after every wrap and
                ## `Q_0 - q(x)` was -3 after three (outputs right, as they
                ## read the state modulo the modulus;
                ## `benchmarks/pss_transient_boundary.py` V10).
                for slot in ('_glm_Q', '_glm_Q_at_entry'):
                    rec = getattr(self, slot, None)
                    if rec is not None:
                        Q = np.array(rec[0], dtype=float)
                        Q[0][row] -= d
                        setattr(self, slot, (Q,) + tuple(rec[1:]))
                shifted = True
        if shifted:
            self._q_cache = None

    def step_lte(self, x_curr, x_last, J):
        """Normalised local truncation error of the step just taken.

        The number the step controller tests against 1 -- `|J^-1 Eg| / etol`
        with `etol = TRTOL (reltol ref + lte_abstol)` -- computed by asking
        the controller itself rather than restating the chain, so the
        estimator, the `J^-1` map, `relref` and the TRTOL folding are the
        ones a transient would have used.  `last_err` is exposed on the
        controllers precisely so it can be read from outside.

        ⚠ Call BEFORE :meth:`_push_history`: it reads `_qlast`/`_iqlast` as
        the PREVIOUS charges, which the push overwrites.

        For a caller that IMPOSES the grid this is not a control signal --
        nothing can act on it -- it is a MEASUREMENT: how far the discrete
        trajectory is from the true one, which is the only thing a Newton
        residual cannot tell you.  Returns None on a step with no history
        to difference.
        """
        from pycircuit.circuit.stepcontroller import IntegralController
        if getattr(self, '_lte_probe', None) is None:
            self._lte_probe = IntegralController().set_relref(self.par.relref)
        ## ⚠ THE LTE FLOORS, NOT THE NEWTON ONES -- the shared helper, not a
        ## second transcription of it: the two `abstol` flavours were split
        ## by stage 0.3d precisely because they are different quantities.
        self._lte_probe.evaluate_step(**self._lte_inputs(
            x_curr, x_last, J, self._dt, self._lte_abstol_vector(), self._dt))
        return self._lte_probe.last_err

    def _lte_inputs(self, x_curr, x_last, J, h, abstol, max_step,
                    clamped=False, x_hist=None):
        """A step controller's inputs (`StepLTEInputs`, as keywords) for the
        step `x_last -> x_curr` of size `h` at this transient's state: its
        charges and past steps, active integrator, reference node, `reltol`,
        TRTOL and node count -- the stepping loop's judge and `step_lte`'s
        one set (written out in each until 2026-10-01)."""
        return {
            'x_curr': x_curr, 'x_last': x_last,
            'q_curr': self._q_at(x_curr),
            'q_last_hist': self._qlast, 'iq_last_hist': self._iqlast,
            'h_curr': h,
            'h_last': self._dt_last if self._dt_last is not None else h,
            'h_last2': self._dt_last2, 'no_history': self._no_history, 'J': J,
            'active_integrator': self.active_integrator,
            'irefnode': self.irefnode, 'reltol': self.par.reltol,
            'abstol': abstol, 'toolkit': self.toolkit, 'max_step': max_step,
            'TRTOL': self.LTERATIO, 'n_nodes': len(self.cir.nodes),
            'h_clamped': clamped, 'x_hist': x_hist}

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
        h_last = self._dt_last if self._dt_last is not None else h
        d_iq = self.active_integrator.companion_dh(q, self._qlast, h, h_last,
                                                   self._iqlast)
        d_u = self.cir.dudt(t, self.epar, analysis=self.par.analysis)
        return self.toolkit.array(d_iq, dtype=float) + \
            self.toolkit.array(d_u, dtype=float)

    def residual_dT(self, x, h=None):
        """``T df/dT`` at fixed solution, for a grid whose EVERY step scales.

        The shooting counterpart of :meth:`residual_dh`.  That one is a
        PARTIAL -- ``d/dh_m`` with the past steps held fixed, which is what
        the coupled time-stepping method needs because there they are held
        fixed.  A shooting analysis solving for the period rebuilds its grid
        at the current ``T``, so every step is ``c_k T`` and the total
        derivative needs every step's route; see
        :meth:`Integrator.companion_dT` for why that is one term rather than
        a second set of partials.

        ⚠ THE SOURCE HALF IS DELIBERATELY ABSENT, and this is only correct
        because of who calls it.  ``residual_dh`` carries ``du/dt`` because
        ``u`` is evaluated at ``t_{m-1} + h``; under a scaling grid the
        source term's route is ``du/dt * t_n``, not ``du/dt * h_n``.  The
        only caller is a shooting analysis's period column, which exists
        ONLY on the autonomous path -- and ``PSS._is_autonomous`` is
        structural: it is true exactly when no source varies with time, so
        ``du/dt`` is identically zero there.  A driven caller would need the
        ``t_n`` route added, and would be wrong without it.

        ⚠ AND THAT IS WHY THERE IS NO ``t`` ARGUMENT.  `residual_dh` takes
        one because it evaluates `du/dt` there; taking one here and not
        reading it is the F8/F18 defect class this tree has a standing guard
        against (`test_no_function_accepts_an_argument_it_never_reads`),
        which caught exactly that within one run of adding this method.  A
        future driven caller adds the argument back WITH the term that uses
        it.
        """
        if h is None:
            h = self._dt
        q = self._q_at(x)
        h_last = self._dt_last if self._dt_last is not None else h
        d_iq = self.active_integrator.companion_dT(q, self._qlast, h, h_last)
        return self.toolkit.array(d_iq, dtype=float)
