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
        """One step's `Jf`, factored by the CALLER'S linear solver.

        ⚠ THIS USED TO REACH FOR `scipy.linalg.lu_factor` DIRECTLY, which
        made every matrix-free run dense-LAPACK whatever `linearsolver=`
        said -- and a circuit Jacobian at m=1002 is very sparse, so a dense
        LU is ~3e8 flops per step, `N` times over.  The recorded 2.13x at
        m=1002 and the m~250 gate were therefore both measured against a
        DENSE baseline; see `benchmarks/pss_matrix_free_ceiling.py` for what
        they become when both sides get a sparse solver.

        A solver whose `factor` returns `None` (the symbolic toolkits, and
        any solver that has not implemented one) falls back to `solve` per
        replayed step -- correct, and `k` times the factorisations, which is
        the cost matrix-free exists to avoid.  It is a fallback, not a mode
        to run in.
        """
        solver = self._get_linearsolver()
        fac = solver.factor(Jf, self.toolkit)
        if fac is not None:
            return fac
        toolkit, A = self.toolkit, Jf

        class _PerSolve(object):
            def solve(self, b):
                return solver.solve(A, b, toolkit)
        return _PerSolve()


    def _k_at(self, x_reduced, t=0.0):
        """The reduced STAGE DERIVATIVE `dq/dt = -(i(x) + u(t))` at a point --
        what the DIRK and coupled-Radau period columns need, and with `t`
        the stage time what a DRIVEN circuit's event columns need
        (2026-09-21): at t = 0 the source is wrong on every other stage of
        a driven circuit, and on the PWM fixture's ramp row that read as a
        -3.9 in the event row's derivative where the FD said +4.2.

        ⚠ THIS USED TO BE `-i(x)` ALONE, on the reasoning that an autonomous
        circuit has "no source term" (2026-09-20).  An autonomous circuit
        has a CONSTANT source vector, not a zero one: a DC supply, a bias
        current.  On a row a source pins, `i(x) = -u` at convergence, so
        the true derivative is 0 and `-i(x)` is `u` -- and the period
        column read `-u/T` there.  Measured on the tree's own phase
        fixture (a 1 kV DC supply): radau's column `[-1e6, ~0, -1.8e-4,
        ...]` against the finite difference `[0, -3.18, 6.28e3, ...]`,
        and on a van der Pol with a decoupled 1 kV node the oscillator rows
        agreed with the FD to 4 digits while the supply node read
        -157.9 = -1e3/T and its branch current +i_R/T.  The first Newton
        step on that column threw `x0` to 2.8e6 and the stage matrix went
        singular there, mislabelled "seed below the fundamental".  Both FD
        checks in the build were on a source-free van der Pol, where the
        two expressions coincide.  `u` is taken at t = 0 because this is
        only built for autonomous circuits, where it is constant."""
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

        Shared by `_traverse_solved_history` and the final replay, because a
        replay that opened differently from the solve would report a
        waveform the residual was never driven to zero on.
        """
        tr = self._transient()
        _alphas, b = tr._get_integrator().companion_coefficients(dt, dt)
        tr._begin_run(self._insert_refnode(xm1_in), self.cir.n)
        tr._dt = dt

        ## ⚠ A `b != 0` COMPANION IS REFUSED HERE, AND THE FIRST REASON
        ## GIVEN FOR IT WAS WRONG.  It said a solved history carries CHARGES
        ## while such a method also reads `iq_{-1}`, "which no charge
        ## determines".  The DAE determines it exactly -- a converged point
        ## satisfies `i(x) + iq + u = 0`, so `iq_{-1} = -(i(x_{-1}) + u)`,
        ## the same identity item 4d rests on -- and seeding it was tried.
        ##
        ## It fails for the derivative running the OTHER way.  A one-step
        ## companion reads only `iq_{-1}`, so the trajectory depends on
        ## `x_{-1}` solely through it, and `d(iq_{-1})/d x_{-1} = -G` is
        ## SINGULAR wherever a node carries no conductance -- every purely
        ## reactive node, which is most of a resonator.  Admitting `x_{-1}`
        ## as m unknowns then leaves the 2m x 2m system rank-deficient:
        ## measured, `LinAlgError: Singular matrix` on 25 tests at once.
        ##
        ## The right second unknown for such a method is `iq_{-1}` ITSELF --
        ## the `(x, iq)` state its monodromy already uses -- closed by
        ## `iq_{-1} = iq_{N-1}`.  That is a different formulation, not a
        ## seeding fix, and it is not built.
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
        ## FIRST.  `x_{-1}` sits one step BEFORE `x_0`, and on a periodic
        ## grid that step is `hs[-1]`.  With a uniform grid the two are
        ## equal and this never showed; on a caller's grid (item 5) with a
        ## 16438:1 spread, handing `hs[0]` to a method that reads `h_last`
        ## states a step ratio that never happened.
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
        (`x_n` and every stage).  Whatever the last step's solve happened to
        leave behind (the LAST stage) was therefore used for all of them, so the
        period map linearised the junction at the wrong voltage.

        Measured against a finite-difference derivative of the discrete period
        map (a reference this code cannot influence) on a diode loaded through a
        series resistor: the analytic monodromy was off by a FIXED 1.65e-3
        (Radau) / 9.1e-4 (TR-BDF2) relative -- flat across four decades of the
        FD step, so a genuine error and not FD noise -- and the error grew with
        how hard the junction was driven, vanishing when it was off.  With this
        sync the same comparison lands at ~3e-9, the FD noise floor.

        `limit(x, x)` sets `_vlim` to `x`'s own branch voltage at zero delta, so
        it moves the state without perturbing the point.  Same defect and same
        remedy as the coupled stage solve in `Transient._rk_step_coupled`.
        """
        tr = self._transient()
        tr.cir.limit(x_full, x_full, tr.epar)

    def _C_at(self, x_reduced):
        """The reduced capacitance at a point, without taking a step.

        ⚠ NO LIMITING SYNC HERE, AND THAT IS MEASURED, NOT ASSUMED.  `_G_at`
        needs the device limiting state to be at the point it is evaluating,
        because a junction's `i`/`G` are read at the stored `_vlim`.  CHARGE IS
        NOT: surveyed across every limiter in the tree,

          * `elements.Diode` is the only STATEFUL one (it keeps `_vlim`), and
            its `C`/`q` do not read it -- with the stored state moved far from
            the evaluation point, `dC = dq = 0` while the control `dG = 15.2`
            and `di = 3.9e-1` confirm the limiting was live;
          * `Semiconductor` (BJT/JFET/ZenerDiode/Varactor) limits STATE-FREE by
            construction -- "Return a limited copy of `x` -- STATE-FREE, and
            that is the point";
          * `compact.PspMosLongChannel` likewise returns a limited copy;
          * the hdl devices keep no `_vlim` at all (it is a codegen local).

        So there is no device whose capacitance a sync could correct.  A sync
        was carried here for a while as "correct in principle" insurance and was
        never exercised by any test -- this tree's own rule is that unexercised
        machinery is a liability.  ⚠ If a stateful limiter whose CHARGE reads its
        state is ever added, this is where the sync goes back; `_sync_limit_at`
        is kept for that, and for `_G_at`'s no-junction path.
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
        stored one, so the answer depends on the POINT ALONE.  The `limit(x, x)`
        route `_C_at` still uses does not: `limit` clamps relative to the STORED
        `_vlim`, so it lands on the true point only when the previous evaluation
        was already nearby.  Measured, varying the prior `_vlim` before
        evaluating at a fixed point: PCNR's `J_eff` moves by 0.0, the limit-sync
        `G` by up to 15.15.  It was right in the traversal only BY LOCALITY
        (steps are small, so the prior state is always close) -- the same
        accident `_begin_period` warns about when it insists the period map be a
        function of `x0` alone, applied to its linearisation.

        Numerically this changes NOTHING today: against a finite difference of
        the discrete period map both routes give the same monodromy to every
        printed digit (3.025e-09 radau / 2.358e-09 trbdf2, identical either
        way).  It removes a latent order-dependence, and it is what lets the
        transient and the monodromy share ONE limiting.

        ⚠ `_C_at` CANNOT JOIN: PCNR re-stamps `i`/`G` at `v_lim` but leaves `q`
        alone (`pcnr.py` treats the algebraic equations; diffusion charge is its
        stated caveat), so the capacitance keeps the limit-sync.
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

        ⚠⚠ THE CHAIN RULE THE `open_at_x0` PATH ASSUMED AWAY.  Every branch
        that opens at `x_0` seeds `Pq = 0` and says so in the same words --
        "no companion current has been formed yet".  That was true of every
        method in this tree until `theta`, which refuses the L-stable opener
        and therefore READS `iq_{-1}` on its first step: `_begin_run` seeds it
        at ``iq_0 = -(i(x_0) + u(t_0))``, the DAE's own `dq/dt`, and that is a
        FUNCTION OF `x_0`.  Differentiating it gives `-G(x_0)`, and dropping
        that term is not a small error -- it is the whole `null(C)` mode.

        Measured on the B2 gate resonator at `K = 200`, against a
        finite-difference of the shooting residual (delta-swept over six
        decades, FLAT, so a real error and not FD noise): the analytic
        monodromy mapped `null(C)` to ZERO -- exactly what an L-stable Euler
        opener would do -- where the true map multiplies it by `-0.7778`,
        which is `(-(1-theta)/theta)^K` from `ThetaIntegrator`'s own table.
        Relative Jacobian error 6.344; with this seed, 1.4e-10.

        The cost of that was NOT a wrong answer -- the residual is what it is,
        so the solve still lands on the right orbit -- but the Newton lost its
        quadratic step: on a LINEAR circuit an exact shooting Newton converges
        in ONE iteration (`trap` with `x0_unknown` takes 3 evaluations at every
        `K`), and `theta` was taking 9 / 64 / 99 at `K = 100 / 200 / 400`.

        `None` -- the default for every method that does NOT declare
        `needs_consistent_iq0` -- means the zero seed is exact, and those
        methods stay bit-identical.
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

        PSS used to carry its OWN transcription of one integrator step --
        the third in the tree, after `Transient` and `JAXTransient` -- and
        it had already cost the two defects its docstring records: `method`
        declared and never read, and a companion current fed back from the
        iterate before the converged one.  Driving the real thing removes
        the copy and brings what came with it: the limiting machinery, PCNR,
        breakpoint order drops, and the continuation rescue.

        ⚠ THAT LIST WAS 1-FOR-4 AS SHIPPED (external review, 2026-09-02); it is
        now 2-FOR-4.  LIMITING reaches -- `cir.limit` is called on the inner
        Newton and the rectifier measurably conducts.  PCNR now reaches too:
        `PSS(cir, pcnr=True)` is a declared Parameter forwarded to the inner
        `Transient` above, and PCNR lives in `Transient.solve_timestep` (the
        LMM `_solve_timestep_pcnr` and, for stage methods, `_rk_stage_pcnr`),
        which PSS DOES call -- so it needs no per-accepted-step machinery.  The
        remaining two still do not: the continuation rescue (`_rescue_solver`)
        and breakpoints (`cir.next_event`) are armed only inside
        `Transient.solve`, which PSS never calls -- it drives `solve_timestep`
        directly on its own frozen grid, so a breakpoint has nothing to move.
        The same structural fact behind the TLine refusal above: what
        `Transient.solve` does per accepted step, PSS does not do at all.

        The tolerances are handed over unchanged, which is the point of
        `newton_tolerance_vectors`: `reltol`/`iabstol`/`vabstol` mean the
        same thing on both sides, so passing them through is a no-op in
        meaning.
        """
        if getattr(self, '_tran', None) is None:
            ## ⚠ A MAPPING, not an if/else on 'euler'.  Written as
            ## `EulerIntegrator() if method == 'euler' else Trapezoidal...`
            ## it silently ran trapezoidal for every other name -- caught
            ## while adding 'gear', which produced numbers identical to
            ## trap's to the last digit.  This class has already paid once
            ## for a `method` that selected nothing; a dict raises KeyError
            ## on a name nobody wired.
            self._tran = self._new_transient(
                self._integrator_for(self.par.method))
        return self._tran

    def _theta_biased(self, integ):
        """Give a `ThetaIntegrator` the bias THIS period needs, not a fixture's.

        ⚠⚠ `theta - 1/2 = C h` makes `C` a RATE, but the quantity that decides
        anything is the DIMENSIONLESS product `C T`: `null(C)` is damped over a
        period by `((1-theta)/theta)^K ~ exp(-4 C h K) = exp(-4 C T)`, and `h`
        cancels.  So `ThetaIntegrator.DEFAULT_C = 1e4` is not a recipe -- it is
        the measured knee `C T = 0.0628` divided by ONE fixture's period
        (6.283e-6 s).  On a circuit 159x slower it is 159x the calibrated bias,
        and MEASURED on `_q20_rlc` (analytic 20 V) that is a peak of 15.91 at
        K = 100 -- 20% low, `converged=True`, because it did converge: to its
        own over-damped discretisation.  See `ThetaIntegrator.DEFAULT_CT`.

        This is where a shooting run stops inheriting that.  `theta_ct` is the
        dimensionless knob (`None` = the measured knee) and the period is this
        solve's, so `method='theta'` is now correct on any circuit.

        ⚠ THE SEED PERIOD IS ENOUGH, and that is a measurement not a hope: the
        gate's own table has `rcond(I - A^K)` at 4.4e-03 / 4.6e-03 / 4.3e-03
        across `C T` = 0.0063 / 0.0628 / 0.628, i.e. FLAT over two decades.  An
        autonomous solve moving `T` by a few percent moves the bias by the same
        few percent, which the knee does not notice.

        Every other method is returned untouched, so nothing else moves a bit.
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

        ⚠ THE SOLVER STRATEGIES GO THROUGH TOO, and they used not to.
        `nrsolver`, `linearsolver` and `scaler` are declared on the base
        `Analysis`, so `PSS(cir, linearsolver=...)` has always been ACCEPTED
        -- and then dropped here, with the inner `Transient` resolving to
        `DenseSolver`/`StandardNewton` whatever the caller asked for.  That
        was the third time this class took a parameter it never read
        (`method` declared and never read; `analysis='PSS'` matching
        nothing), and the same shape each time: accepted at the constructor,
        silently discarded at the boundary.
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
        ## never arms the transient's rescue ladder (owner decision
        ## 2026-09-08, "Do 2"; see `_rk_step_coupled` and `solve_timestep`).
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

        This used to be a private transcription of one integrator step --
        the third in the tree -- and it had already cost two defects that
        its own comments recorded: `method` was declared and never read, so
        PSS was backward-Euler only, and the companion current fed back to
        the next step belonged to the iterate BEFORE the converged one.
        Both are structurally impossible now: the integrator is an
        `Integrator` object driven by `Transient.get_diff`, and the
        companion current is the one that class stores at its own converged
        point.

        What came with the change, none of which the copy had: the limiting
        machinery (measured -- a rectifier whose diode never turned on, so
        the non-conducting solution was returned as a converged periodic
        steady state), PCNR when the circuit and Parameters ask for it,
        breakpoint order drops, and the continuation rescue.

        `dt` is imposed by the caller: PSS owns the grid, which is what
        keeps the period map a smooth function of `x0`.  `iq_last` is
        retained in the signature for callers that pass it, but the
        companion history now lives in the `Transient`'s own ring buffers,
        rolled here through `_push_history`.

        The measured cost of backward Euler on a limit cycle is unchanged
        and still the reason `method` matters -- it damps exactly what PSS
        exists to find:

            steps/period    PSS peak    fraction of analytic
                      20      2.63 V       13.2%
                      50      5.61 V       28.1%
                     100      8.81 V       44.1%
                     200     12.20 V       61.0%
        """
        toolkit = self.toolkit
        irefnode = self.irefnode
        tr = self._transient()

        ## ONE INTEGRATOR STEP, taken by the class that owns the definition.
        ## `Transient.solve_timestep` applies the chosen integrator through
        ## `get_diff` (so `method` selects something because the integrator
        ## object does), the limiting machinery, PCNR when asked for, and the
        ## continuation rescue.  None of that existed on the copy this
        ## replaced.
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
            ## the event columns of a two-step companion need BOTH partials
            ## (2026-09-22): its own step's and the previous step's, the
            ## latter assembled from the total under uniform scaling
            (self._dfdh,) = remove_row_col(
                (tr.residual_dh(x_full, t, dt),), irefnode, toolkit)
            (self._dfdT,) = remove_row_col(
                (tr.residual_dT(x_full, dt),), irefnode, toolkit)
        elif self._want_dfdh:
            ## ⚠ `residual_dT`, NOT `residual_dh`.  The grid is rebuilt at
            ## the current `T` on every residual evaluation, so every step
            ## scales and the partial `d/dh_n` is not the total -- for
            ## Gear-2 it is 3/2 of it on a uniform grid, measured against
            ## finite differences at 1.4859/1.4939/1.4972 for 100/200/400
            ## points, converging on the exact 3/2.  Euler and trapezoidal
            ## were never wrong: their coefficients depend on `h_n` alone,
            ## so the partial IS the total, which is why only Gear-2 was
            ## hit.  See `Integrator.companion_dT`.
            if self._period_column == 'closing':
                ## ⚠ `residual_dh`, THE PARTIAL, NOT `residual_dT`.  The
                ## note above explains why the total is 3/2 of the partial
                ## for Gear-2: `residual_dT` accounts for EVERY step scaling
                ## with `T`.  Under the closing-step convention only ONE
                ## step's `h` moves, so the partial IS the derivative and
                ## the 3/2 would be exactly the error.
                (self._dfdh,) = remove_row_col(
                    (tr.residual_dh(x_full, t, dt),), irefnode, toolkit)
            else:
                (self._dfdT,) = remove_row_col(
                    (tr.residual_dT(x_full, dt),), irefnode, toolkit)
        ## Measured, not controlled: the grid is the caller's, so nothing can
        ## act on this.  Also before the push, for the same reason.
        if self._want_lte:
            self._lte = tr.step_lte(x_full, self._insert_refnode(x0), J_full)
            ## A SEAM STEP IS ONE WHOSE COMPANION READS THE ENTERING
            ## UNKNOWN, not merely one whose ESTIMATOR does.
            ##
            ## ⚠ THIS CONDITION WAS `h_last2 is None` ALONE, AND THAT FLAGGED
            ## A PHANTOM FOR TWO METHODS OF THREE.  `h_last2 is None` is the
            ## transient's statement that the third past charge is not real,
            ## which is the reach of the LTE estimator's third divided
            ## difference -- not the reach of the integrator.  Euler's
            ## companion reads `q_{n-1}`; trapezoidal's reads `q_{n-1}` and
            ## `iq_{n-1}`, and the order-dropped opening step supplies an
            ## `iq` consistent with it, which is what that drop is FOR.
            ## Neither can see the fabricated charge at all.  Measured
            ## (`benchmarks/pss_seam_cost.py`): trapezoidal's seam reading
            ## was 15.1 times tolerance while its cost is 1.3e-11 V, and
            ## euler's 0.286 against 5.1e-12 V.  Both are exactly zero; the
            ## reading was an artefact of the measurement.
            ##
            ## Gear-2 reads `q_{n-2}` -- which at that step IS the entering
            ## unknown -- and the shooting condition constrains `x(0)` to
            ## equal `x(P)`, NOT `x_in` to be the orbit's own `x(-dt)`.  So
            ## `x_in` sits O(h^2) off a real history point and Gear-2 reads
            ## it as one.  That one costs 1.266e-01 V at 100 points/period
            ## against an interior contribution of 1.070e-01 -- the seam is
            ## 54% of its total error -- and it is the term that STOPS
            ## converging: it falls as h^2 while the interior falls faster,
            ## so its share grows to 68% at 200 points and 73% at 400.
            ##
            ## TWO DIFFERENT THINGS ARE TRUE OF AN OPENING STEP, and the
            ## first version of this conflated them -- which showed up as
            ## trapezoidal's phantom simply MOVING from the seam into the
            ## interior total (0.340 -> 15.47) when the seam test was
            ## tightened.  Suppressing a bad number is not the same as
            ## classifying it.
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
            ##   deep enough to touch it.  Only Gear-2's is.
            ##
            ## `_dt_last2 is None` says exactly "two real past charges" here
            ## (it is set from `_dt_last` one step later), so it is the step
            ## index in disguise; both tests are written against the count.
            _real_past = 2 if tr._dt_last2 is None else 3
            _p = getattr(tr.active_integrator, 'ORDER', 1) + 1
            _reach = len(tr._companion_coeffs[0]) - 1
            self._lte_valid = _real_past >= _p
            ## `_history_is_solved` is that formulation saying the
            ## deepest charge is an UNKNOWN the solve closed, not a stand-in
            ## -- so there is no seam to report even though the companion
            ## reaches that far.  Without this the fix would go on flagging
            ## the defect it removed.
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
