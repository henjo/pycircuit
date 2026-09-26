"""The inner transient the period walks step with, and the circuit linearised
at a point (C, G, the stage derivative).
"""
import numpy as np
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.circuit import gnd


class _InnerTransient(object):
    """The inner transient the period walks step with, and the circuit
    linearised at a point (C, G, the stage derivative).  A theme of `PSS`
    (see `pss.py`)."""

    def _factorise(self, Jf):
        """One step's `Jf`, factored by the CALLER'S linear solver -- never a
        dense LU directly, which would make every matrix-free run dense-LAPACK
        whatever `linearsolver=` said, on a very sparse circuit Jacobian
        (`benchmarks/pss_matrix_free_ceiling.py`).

        A solver whose `factor` returns `None` (the symbolic toolkits, and
        any solver that has not implemented one) falls back to `solve` per
        replayed step -- correct, and `k` times the factorisations, which is
        the cost matrix-free exists to avoid.  It is a fallback, not a mode
        to run in.

        History: `doc/shooting_history.md`, `_InnerTransient._factorise`.
        """
        solver = self._get_linearsolver()
        fac = solver.factor(Jf, self.toolkit)
        if fac is not None:
            return fac
        toolkit, A = self.toolkit, Jf

        class _PerSolve(object):
            def solve(self, b):
                return solver.solve(A, b, toolkit)

            ## ⚠ the adjoint surfaces transpose every stored step; without
            ## this a solver with no `factor` raised `AttributeError` there
            ## (gear's `ppv` under `KLUSolver`/`AutoSolver` until 2026-09-25)
            def solve_transposed(self, b):
                return solver.solve(A.T, b, toolkit)
        return _PerSolve()


    def _k_at(self, x_reduced, t=0.0):
        """The reduced STAGE DERIVATIVE `dq/dt = -(i(x) + u(t))` at a point --
        what the DIRK and coupled-Radau period columns need, and with `t`
        the stage time what a DRIVEN circuit's event columns need: at t = 0
        the source is wrong on every other stage of a driven circuit.

        ⚠ `u` IS NOT OPTIONAL ON AN AUTONOMOUS CIRCUIT.  Its source vector is
        CONSTANT, not zero (a DC supply, a bias current): on a row a source
        pins, `i(x) = -u` at convergence, so the true derivative is 0 and
        `-i(x)` alone would read `u`.  A source-free fixture cannot tell the
        two expressions apart.  The default `t = 0` is right for an
        autonomous circuit, where `u` is constant.

        History: `doc/shooting_history.md`, `_InnerTransient._k_at`."""
        tr = self._transient()
        xf = self._insert_refnode(x_reduced)
        k = -(np.asarray(tr.cir.i(xf, tr.epar), dtype=float)
              + np.asarray(tr.cir.u(float(t), tr.epar,
                                    analysis=self.par.analysis), dtype=float))
        iref = self.irefnode
        return self.toolkit.concatenate((k[:iref], k[iref + 1:]))

    def _install_history(self, x0_in, xm1_in, dt, h_prev=None):
        """Open a run ON a solved two-point history rather than a seed.

        `_begin_run(x_{-1})` opens the rings on the earlier point and the
        push puts `x_0` in front of it, so the first real step reads
        `q(x_0)` and `q(x_{-1})` -- two genuine solved points.  The flags
        then say what is true of them: a step of `dt` has been taken, and
        the run is no longer opening, so nothing drops order.

        `_dt_last2` stays None on purpose.  The THIRD charge in the ring is
        still `q(x_{-1})` repeated, and the LTE estimator differences three,
        so its opening reading remains unsound and the report goes on
        discarding it -- a solved history fixes what the SOLUTION reads,
        not what the estimator does.

        Shared by the pair walk (`_walk_lmm`) and the final replay, because a
        replay that opened differently from the solve would report a
        waveform the residual was never driven to zero on.
        """
        tr = self._transient()
        _alphas, b = tr._get_integrator().companion_coefficients(dt, dt)
        tr._begin_run(self._insert_refnode(xm1_in), self.cir.n)
        tr._dt = dt

        ## ⚠ A `b != 0` COMPANION IS REFUSED HERE.  Such a method depends on
        ## `x_{-1}` only through `iq_{-1} = -(i(x_{-1}) + u)` (exact at a
        ## converged point), whose derivative `-G` is SINGULAR at every purely
        ## reactive node -- so admitting `x_{-1}` as m unknowns leaves the
        ## 2m x 2m system rank-deficient.  The right second unknown for such a
        ## method is `iq_{-1}` ITSELF, closed by `iq_{-1} = iq_{N-1}` -- a
        ## different formulation, not a seeding fix, and not built.
        if b:
            raise NotImplementedError(
                'a solved entering history admits `x_{-1}` as the second '
                'unknown, and a companion with a b != 0 term depends on it '
                'only through `iq_{-1} = -(i(x_{-1}) + u)`, whose derivative '
                '-G is singular at every purely reactive node -- so the '
                'enlarged system would be rank-deficient. Such a method '
                'needs `iq_{-1}` itself as the unknown, which is a different '
                'formulation; this refuses rather than solving a singular '
                'one.')
        ## The charge half of `_push_history`, without its `_iq` roll: no
        ## step has been solved yet, so there is no companion current to
        ## push -- and a `b = 0` companion never reads one, which the guard
        ## above is what makes true.
        q0 = tr.cir.q(self._insert_refnode(x0_in), tr.epar)
        tr._qlast = self.toolkit.concatenate(
            (self.toolkit.array([q0]), tr._qlast))[:-1]
        tr._q_cache = None
        tr._is_first_step = False
        tr._no_history = False
        ## ⚠ THE STEP THAT PRODUCED `x_0` IS THE PERIOD'S LAST ONE, NOT ITS
        ## FIRST: `x_{-1}` sits one step BEFORE `x_0`, and on a periodic grid
        ## that step is `hs[-1]`.  On a non-uniform grid, handing `hs[0]` to a
        ## method that reads `h_last` states a step ratio that never happened.
        tr._dt_last = dt if h_prev is None else h_prev
        tr._dt_last2 = None
        self._history_is_solved = True
        return tr

    def _sync_limit_at(self, x_full):
        """Put every limiting device's internal state AT ``x_full``.

        ⚠ EVALUATING "AT A POINT" REQUIRES THE DEVICE LIMITING STATE TO BE AT
        THAT POINT.  A junction device's `i`/`G` are read at its stored `_vlim`,
        not at the vector handed in, and there is only ONE `_vlim` per device --
        while the monodromy evaluates `C`/`G` at SEVERAL distinct points per step
        (`x_n` and every stage).  Without the sync, whatever the last solve
        left behind (the LAST stage) is used for all of them, and the period
        map linearises the junction at the wrong voltage.

        `limit(x, x)` moves the state without perturbing the point -- but it
        clamps against the STORED state, so above a junction's critical
        voltage it lands on `x` only from a state already near it (measured
        on `algebraic_conditioning`, which resets first for that reason).
        `_G_at` calls this only where no junction limits.  Same defect and
        same remedy as the coupled stage solve in
        `Transient._rk_step_coupled`.

        History: `doc/shooting_history.md`, `_InnerTransient._sync_limit_at`.
        """
        tr = self._transient()
        tr.cir.limit(x_full, x_full, tr.epar)

    def _C_at(self, x_reduced):
        """The reduced capacitance at a point, without taking a step.

        ⚠ NO LIMITING SYNC HERE.  `_G_at` needs the device limiting state at
        the point it evaluates, because a junction's `i`/`G` are read at the
        stored `_vlim`; CHARGE is not.  `elements.Diode` is the only STATEFUL
        limiter and its `C`/`q` do not read `_vlim`; `Semiconductor` and
        `compact.PspMosLongChannel` limit state-free, and the hdl devices keep
        no `_vlim` at all.
        ⚠ If a stateful limiter whose CHARGE reads its state is ever added,
        this is where the sync goes back; `_sync_limit_at` is kept for that,
        and for `_G_at`'s no-junction path.

        History: `doc/shooting_history.md`, `_InnerTransient._C_at`.
        """
        tr = self._transient()
        xf = self._insert_refnode(x_reduced)
        C = tr.cir.C(xf, tr.epar)
        (C,) = remove_row_col((C,), self.irefnode, self.toolkit)
        return C

    def _G_at(self, x_reduced):
        """The reduced conductance `di/dx` at a point, without taking a step.

        The companion to `_C_at`.  The one-step-method monodromy needs the
        PHYSICAL `(C, G)` at each stage point -- not the companion
        `Geq = a h G` the accepted step happens to store -- because a
        two-stage DIRK linearises `q` and `i` at THREE distinct points per
        step (`x_n`, the internal stage, and `x_{n+1}`), each with its own
        stage coefficient.  Recovering a physical `G` from a single stored
        `Geq` would divide out only one of those coefficients and mislabel
        the other two.

        ⚠ THIS GOES THROUGH PCNR WHEN THERE ARE JUNCTIONS, and the reason is
        STATELESSNESS, not accuracy.  `pcnr.augmented_system` + `schur_reduce`
        build `G` from an explicitly-passed `v_lim` instead of from the device's
        stored one, so the answer depends on the POINT ALONE.  A `limit(x, x)`
        sync clamps relative to the STORED `_vlim`, so it lands on the true
        point only when the previous evaluation was already nearby -- the
        order-dependence `_begin_period` forbids for the period map, applied
        to its linearisation.  It also lets the transient and the monodromy
        share ONE limiting.

        ⚠ `_C_at` CANNOT JOIN: PCNR re-stamps `i`/`G` at `v_lim` but leaves `q`
        alone (`pcnr.py` treats the algebraic equations; diffusion charge is its
        stated caveat), so the capacitance is read directly -- with no limit
        sync either, since no charge in the library reads the limiting state
        (see `_C_at`).

        History: `doc/shooting_history.md`, `_InnerTransient._G_at`.
        """
        tr = self._transient()
        xf = self._insert_refnode(x_reduced)
        junctions = self._pcnr_junctions()
        if not junctions:
            ## no junction devices: nothing limits, so the plain read IS the
            ## physical G and PCNR would only add an assembly for no reason.
            self._sync_limit_at(xf)
            G = tr.cir.G(xf, tr.epar)
        else:
            from pycircuit.circuit import pcnr as _pcnr
            xfa = np.asarray(xf, dtype=float)
            v_lim = _pcnr.v_lim_init(junctions, xfa)
            g_mna, g_lim, J_mm, _J_ml, _J_lm, didv = _pcnr.augmented_system(
                tr.cir, xfa, v_lim, junctions, tr.epar,
                u_extra=0.0, dense_blocks=False, J_extra=0.0)
            _f_eff, G = _pcnr.schur_reduce(
                g_mna, g_lim, J_mm, junctions=junctions, didv=didv)
            G = np.asarray(G)
        (G,) = remove_row_col((G,), self.irefnode, self.toolkit)
        return G

    def _pq_seed_at_x0(self, x_reduced):
        """``d(iq_0)/d(x_0)`` when the method SEEDS a consistent companion current.

        ⚠ THE CHAIN RULE THE `open_at_x0` PATH ASSUMES AWAY.  Every branch
        that opens at `x_0` seeds `Pq = 0` ("no companion current has been
        formed yet").  A method that refuses the L-stable opener (`theta`)
        READS `iq_{-1}` on its first step instead: `_begin_run` seeds it at
        ``iq_0 = -(i(x_0) + u(t_0))``, a FUNCTION OF `x_0` whose derivative is
        `-G(x_0)`.  Dropping that term loses the whole `null(C)` mode from
        the monodromy -- the orbit is still right (the residual is what it
        is), but the Newton loses its quadratic step.

        `None` -- the default for every method that does NOT declare
        `needs_consistent_iq0` -- means the zero seed is exact, and those
        methods stay bit-identical.

        History: `doc/shooting_history.md`, `_InnerTransient._pq_seed_at_x0`.
        """
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        if not integ.needs_consistent_iq0():
            return None
        ## `u(t_0)` carries no `x`, so only `i` contributes: `d(-i)/dx = -G`.
        return -np.asarray(self._G_at(x_reduced), dtype=float)

    def _pcnr_junctions(self):
        """The circuit's PCNR-participating devices, found once and cached.

        `pcnr_devices` walks every element and rebuilds the node map, which is
        far too much to repeat inside `_G_at` -- the monodromy calls it once per
        stage per step per traversal.
        """
        junc = getattr(self, '_pcnr_junctions_cache', None)
        if junc is None:
            from pycircuit.circuit import pcnr as _pcnr
            junc = _pcnr.pcnr_devices(self.cir)
            self._pcnr_junctions_cache = junc
        return junc

    def _transient(self):
        """The `Transient` this analysis integrates with.

        PSS drives the real `Transient` rather than a transcription of one
        integrator step, so it inherits LIMITING (`cir.limit` is called on the
        inner Newton) and PCNR (`PSS(cir, pcnr=True)` is forwarded to the
        inner `Transient`, and PCNR lives in `Transient.solve_timestep`, which
        PSS calls).  ⚠ It does NOT inherit the continuation rescue
        (`_rescue_solver`) or breakpoints (`cir.next_event`): both are armed
        only inside `Transient.solve`, which PSS never calls -- it drives
        `solve_timestep` directly on its own frozen grid, so a breakpoint has
        nothing to move.  What `Transient.solve` does per accepted step, PSS
        does not do at all (the same structural fact behind the TLine
        refusal).

        The tolerances are handed over unchanged, which is the point of
        `newton_tolerance_vectors`: `reltol`/`iabstol`/`vabstol` mean the
        same thing on both sides, so passing them through is a no-op in
        meaning.

        History: `doc/shooting_history.md`, `_InnerTransient._transient`.
        """
        if getattr(self, '_tran', None) is None:
            ## ⚠ A MAPPING (`_integrator_for`), not an if/else on 'euler':
            ## an if/else silently runs its fallback for every other name,
            ## while a dict raises KeyError on a name nobody wired.
            self._tran = self._new_transient(
                self._integrator_for(self.par.method))
        return self._tran

    def _theta_biased(self, integ):
        """Give a `ThetaIntegrator` the bias THIS period needs, not a fixture's.

        ⚠ `theta - 1/2 = C h` makes `C` a RATE, but the quantity that decides
        anything is the DIMENSIONLESS product `C T`: `null(C)` is damped over a
        period by `((1-theta)/theta)^K ~ exp(-4 C h K) = exp(-4 C T)`, and `h`
        cancels.  So `ThetaIntegrator.DEFAULT_C` is not a recipe -- it is one
        fixture's knee divided by that fixture's period, and on a slower
        circuit it over-damps: the solve still converges, to its own
        over-damped discretisation.  See `ThetaIntegrator.DEFAULT_CT`.

        Here `theta_ct` is the dimensionless knob (`None` = the measured knee)
        and the period is this solve's, so `method='theta'` is correct on any
        circuit.  The SEED period is enough: the knee is flat over two decades
        of `C T`, so an autonomous solve moving `T` by a few percent does not
        matter.

        Every other method is returned untouched, so nothing else moves a bit.

        History: `doc/shooting_history.md`, `_InnerTransient._theta_biased`.
        """
        from pycircuit.circuit.integrator import ThetaIntegrator
        T = getattr(self, '_theta_period', None)
        if not isinstance(integ, ThetaIntegrator) or T is None:
            return integ
        ct = getattr(self.par, 'theta_ct', None)
        integ.cbias = float(ThetaIntegrator.DEFAULT_CT if ct is None
                            else ct) / float(T)
        return integ

    def _new_transient(self, integ):
        """A `Transient` on this circuit driven by `integ`, with every
        strategy and tolerance PSS was given handed through.

        Extracted so the TR-BDF2 monodromy can build a transient on a
        DIFFERENT integrator than `self.par.method` without duplicating the
        pass-through -- and so the pass-through cannot drift between the two
        call sites.

        ⚠ THE SOLVER STRATEGIES GO THROUGH TOO: `nrsolver`, `linearsolver`
        and `scaler` are declared on the base `Analysis`, so `PSS(cir,
        linearsolver=...)` is accepted at the constructor and must not be
        dropped at this boundary.

        History: `doc/shooting_history.md`, `_InnerTransient._new_transient`.
        """
        from pycircuit.circuit.transient import Transient
        ## a frozen grid has no stalled estimate: the shrink drop to Euler is
        ## off here (see `Gear2Integrator.shrink_guard`); growth stays guarded
        try:
            integ.shrink_guard = False
        except Exception:                                      # noqa: BLE001
            pass
        ## ⚠ ONE CHOKE POINT for the theta bias, because there are four call
        ## sites building a transient and a per-site fix would drift.  A no-op
        ## for every other integrator -- see `_theta_biased`.
        integ = self._theta_biased(integ)
        tr = Transient(
            self.cir, toolkit=self.toolkit, integrator=integ,
            reltol=self.par.reltol, iabstol=self.par.iabstol,
            vabstol=self.par.vabstol, maxiter=self.par.maxiter,
            analysis=self.par.analysis,
            lte_vabstol=self.par.lte_vabstol,
            lte_iabstol=self.par.lte_iabstol,
            TRTOL=self.par.TRTOL, relref=self.par.relref,
            nrsolver=self.par.nrsolver,
            linearsolver=self.par.linearsolver,
            scaler=self.par.scaler,
            pcnr=self.par.pcnr)
        ## The line search as the last resort on the shooting path, which
        ## never arms the transient's rescue ladder (an owner decision; see
        ## `_rk_step_coupled` and `solve_timestep`).
        tr._damped_last_resort = True
        tr.irefnode = self.irefnode
        return tr

    def _begin_period(self, x_reduced):
        """Start one traversal of the period from a clean integrator state.

        Every shooting iteration re-integrates the SAME interval from its own
        `x0`, so "begin a run" happens once per iteration here, not once per
        analysis.  Without the reset, iteration k+1 would inherit the ring
        buffers iteration k ended with and the period map would depend on
        which iteration it was -- phi must be a function of `x0` alone or the
        monodromy is the derivative of something else.
        """
        tr = self._transient()
        tr._begin_run(self._insert_refnode(x_reduced), self.cir.n)
        return tr

    def _insert_refnode(self, x):
        return self.toolkit.concatenate(
            (x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))

    def solve_timestep(self, x0, t, dt, refnode=gnd, iq_last=None):
        """One timestep of the inner transient, taken by `Transient`.

        The integrator is an `Integrator` object driven by
        `Transient.get_diff`, so `method` selects the step, and the companion
        current fed to the next step is the one that class stores at its own
        converged point.  What does and does not reach from `Transient` is
        listed at `_transient`.

        `dt` is imposed by the caller: PSS owns the grid, which is what
        keeps the period map a smooth function of `x0`.  `iq_last` is
        retained in the signature for callers that pass it, but the
        companion history now lives in the `Transient`'s own ring buffers,
        rolled here through `_push_history`.

        Backward Euler damps exactly what PSS exists to find (a limit
        cycle's amplitude), which is why `method` matters.

        History: `doc/shooting_history.md`, `_InnerTransient.solve_timestep`.
        """
        toolkit = self.toolkit
        irefnode = self.irefnode
        tr = self._transient()

        ## ONE INTEGRATOR STEP, taken by the class that owns the definition.
        ## `Transient.solve_timestep` applies the chosen integrator through
        ## `get_diff` (so `method` selects something because the integrator
        ## object does), the limiting machinery and PCNR when asked for.  Not
        ## the continuation rescue: only `Transient.solve`'s stepping loop
        ## arms it (see `_transient`); this path's last resort is the line
        ## search.
        tr._dt = dt
        x_full, J_full = self._transient_step(tr, x0, t)

        ## d(residual)/dh at the converged point, for the period
        ## derivative -- BEFORE `_push_history`, because it reads `_qlast`
        ## as the PREVIOUS charges, which the push is about to overwrite.
        ## `Transient.residual_dh` is Fang's `p`, already shared: it is
        ## `d(iq)/dh + du/dt`, and for an AUTONOMOUS circuit the second term
        ## is identically zero -- which is exactly why solving for the
        ## period is tractable here and would not be on a driven circuit,
        ## where scaling T also moves every source evaluation.
        if getattr(self, '_want_event_cols', False):
            ## the event columns of a two-step companion need BOTH partials:
            ## its own step's and the previous step's, the latter assembled
            ## from the total under uniform scaling
            (self._dfdh,) = remove_row_col(
                (tr.residual_dh(x_full, t, dt),), irefnode, toolkit)
            (self._dfdT,) = remove_row_col(
                (tr.residual_dT(x_full, dt),), irefnode, toolkit)
        elif self._want_dfdh:
            ## ⚠ `residual_dT`, NOT `residual_dh`.  The grid is rebuilt at
            ## the current `T` on every residual evaluation, so every step
            ## scales and the partial `d/dh_n` is not the total (for Gear-2
            ## the total is 3/2 of it on a uniform grid).  Euler and
            ## trapezoidal coefficients depend on `h_n` alone, so for them
            ## the partial IS the total.  See `Integrator.companion_dT`.
            if self._period_column == 'closing':
                ## ⚠ `residual_dh`, THE PARTIAL, NOT `residual_dT`.  The
                ## note above explains why the total is 3/2 of the partial
                ## for Gear-2: `residual_dT` accounts for EVERY step scaling
                ## with `T`.  Under the closing-step convention only ONE
                ## step's `h` moves, so the partial IS the derivative and
                ## the 3/2 would be exactly the error.
                (self._dfdh,) = remove_row_col(
                    (tr.residual_dh(x_full, t, dt),), irefnode, toolkit)
                ## and the total, from which a TWO-STEP method's opening step
                ## takes its previous-step partial (`_walk_lmm`)
                (self._dfdT,) = remove_row_col(
                    (tr.residual_dT(x_full, dt),), irefnode, toolkit)
            else:
                (self._dfdT,) = remove_row_col(
                    (tr.residual_dT(x_full, dt),), irefnode, toolkit)
        ## Measured, not controlled: the grid is the caller's, so nothing can
        ## act on this.  Also before the push, for the same reason.
        if self._want_lte:
            self._lte = tr.step_lte(x_full, self._insert_refnode(x0), J_full)
            ## A SEAM STEP IS ONE WHOSE COMPANION READS THE ENTERING
            ## UNKNOWN, not merely one whose ESTIMATOR does.  TWO DIFFERENT
            ## THINGS ARE TRUE OF AN OPENING STEP:
            ##
            ##   the ESTIMATE is invalid when the ESTIMATOR differences a
            ##   charge that was never a real point.  Its divided difference
            ##   reaches `p = ORDER + 1` charges back -- 2 for Euler, 3 for
            ##   the second-order pair -- so at the step with only two real
            ##   past charges, euler's estimate is sound and trap's and
            ##   gear2's are not.  An unsound estimate is DISCARDED: it is
            ##   not an interior reading and, on its own, not a seam either.
            ##
            ##   a SEAM exists when the COMPANION reads the entering
            ##   unknown, i.e. when its charge reach `len(alphas) - 1` is
            ##   deep enough to touch it.  Only Gear-2's is: it reads
            ##   `q_{n-2}`, the entering unknown, and the shooting condition
            ##   constrains `x(0) = x(P)`, NOT `x_in` to be the orbit's own
            ##   `x(-dt)` -- so `x_in` sits O(h^2) off a real history point,
            ##   an error that falls only as h^2 and so dominates as the grid
            ##   refines.
            ##
            ## `_dt_last2 is None` says exactly "two real past charges" here
            ## (it is set from `_dt_last` one step later), so it is the step
            ## index in disguise; both tests are written against the count.
            _real_past = 2 if tr._dt_last2 is None else 3
            _p = getattr(tr.active_integrator, 'ORDER', 1) + 1
            _reach = len(tr._companion_coeffs[0]) - 1
            self._lte_valid = _real_past >= _p
            ## `_history_is_solved` says the deepest charge is an UNKNOWN the
            ## solve closed (`_install_history`), not a stand-in -- so there
            ## is no seam to report even though the companion reaches that
            ## far.
            self._lte_seam = (_reach >= _real_past
                              and not self._history_is_solved)

        ## The history advance is the accept path's, called rather than
        ## copied -- and `_dt_last` must roll AFTER the step, because
        ## `get_diff` read it as `h_last` while solving.
        tr._push_history(x_full)
        tr._dt_last2 = tr._dt_last
        tr._dt_last = dt
        tr._is_first_step = False
        tr._no_history = False

        ## Reduced-system views for the shooting Jacobian.  `_Geq` is the
        ## companion conductance the step actually used, which is the factor
        ## the monodromy needs; `_iq` is kept for the caller's own bookkeeping
        ## as before.
        (self._Jf, self._Geq, self._C) = remove_row_col(
            (J_full, tr._Geq, tr._Cmat), irefnode, toolkit)
        ## The coefficients of the integrator that ACTUALLY ran this step --
        ## an order drop on the opening step reports Euler's, which is what
        ## the propagation must use for that step.
        self._coeffs = tr._companion_coeffs
        self._iq = tr._iq

        x = toolkit.concatenate((x_full[:irefnode], x_full[irefnode + 1:]))
        return x

    def _transient_step(self, tr, x0_reduced, t):
        """`Transient.solve_timestep` on the FULL vector, returning
        ``(x, J)``.  PSS works on the reduced system throughout; this is the
        one place the two conventions meet."""
        x_full = self._insert_refnode(x0_reduced)
        x, _feval, J, _f = tr.solve_timestep(x_full, t)
        return x, J
