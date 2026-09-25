"""The period walks: one per step family (`_walk_lmm`, `_walk_stage`,
`_walk_glm`), dense and/or factored, and their views.
"""
from copy import copy
import numpy as np
from ._factored import _PeriodWalk
from ._steps import _GLMStartup
from ._steps import _StageStep
from ._steps import _butcher
from ._steps import _GLMStep
from ._steps import _lmm_recursion


class _PeriodWalks(object):
    """The period walks: one per step family (`_walk_lmm`, `_walk_stage`,
    `_walk_glm`), dense and/or factored, and their views.  A theme of `PSS`
    (see `pss.py`)."""

    def _step_sensitivity(self, Px, Cs, Pq, Jf, C_new, solve=None,
                          coeffs=None, source=None):
        """One step of the sensitivity recursion, for ANY seed width.

        ONE RECURSION FOR EVERY METHOD AND BOTH FORMULATIONS.
        Each writes its companion as `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`,
        so differentiating the step gives

            S    = sum_{k>=1} a_k C_{n-k} P_{n-k} + b Pq
            P_n  = -Jf_n^-1 S
            Pq_n = a_0 C_n P_n + S

        `P` is `d x_j / d(unknowns)`: one block wide in the plain
        formulation, two when the history is solved for.  Nothing here
        depends on that
        width, which is why the two systems share this and not a copy.

        ⚠ A SOLVE, NOT AN INVERSE: a dense `inv(Jf)` per step would square
        the condition number it then multiplies through.

        `solve` overrides how that solve is taken, so the MATRIX-FREE path
        can hand in a PRE-FACTORED `Jf` without a second copy of this
        recursion (which every method and both formulations share).  Without
        it, matrix-free would refactor `Jf` once per step PER KRYLOV
        ITERATION -- `k` times the factorisations the dense path takes.

        History: `doc/shooting_history.md`, `_PeriodWalks._step_sensitivity`.
        """
        ## `coeffs` overrides the LIVE `_coeffs` for the same reason `solve`
        ## overrides the solve: a matrix-free replay happens after the run,
        ## when `_coeffs` no longer describes the step being replayed.  See
        ## `_walk_lmm`.
        ## ⚠ `source` ENTERS THE SOLVE AND NOT `Pq`: a small-signal source
        ## appears in the step's residual -- `Jf dx + S + du = 0` -- but NOT
        ## in the companion, which is built from CHARGES.  Adding it to `Pq`
        ## would feed a fictitious charge forward into every later step, an
        ## error that grows along the period rather than announcing itself.
        ## That one term is all that separates the monodromy (`source=None`)
        ## from PAC's forced response, which is why they share this recursion.
        alphas, b = self._coeffs if coeffs is None else coeffs
        if solve is None:
            ## ⚠ THROUGH THE CALLER'S SOLVER, not `toolkit.linearsolver`:
            ## this is the DENSE propagation -- the thing matrix-free is
            ## measured against -- so both sides go through the same
            ## strategy.  `DenseSolver` IS `toolkit.linearsolver`, so the
            ## default is unchanged.
            def solve(S_solve):
                return self._get_linearsolver().solve(Jf, S_solve, self.toolkit)
        return _lmm_recursion(Px, Cs, Pq, C_new, alphas, b, solve, source)

    def _walk(self, kind, z, times, hs, T=None, dense=True, keep=False,
              want_dT=False, hsens=None, capture=None, open_at_x0=False):
        """ONE WALK OF THE PERIOD for every map, from the unknown `z` (gear's
        pair stacked): `_walk_stage` for 'stage', `_walk_glm` for 'glm',
        `_walk_lmm` for 'plain' and 'pair'.  Returns the `_PeriodWalk`."""
        if kind == 'glm':
            return self._walk_glm(z, T, times, hs, dense=dense, keep=keep,
                                  want_dT=want_dT, hsens=hsens,
                                  capture=capture)
        if kind == 'stage':
            return self._walk_stage(z, T, times, hs, dense=dense, keep=keep,
                                    want_dT=want_dT, hsens=hsens,
                                    capture=capture)
        if kind == 'pair':
            m = self.cir.n - 1
            w = self._walk_lmm(('pair', z[:m], z[m:]), times, hs, T=T,
                               dense=dense, keep=keep, want_dT=want_dT,
                               hsens=hsens, capture=capture)
            w.z = z
            return w
        return self._walk_lmm(('plain', z, open_at_x0), times, hs, T=T,
                              dense=dense, keep=keep, want_dT=want_dT,
                              hsens=hsens, capture=capture)

    def _walk_lmm(self, opening, times, hs, T=None, dense=True, keep=False,
                  want_dT=False, hsens=None, capture=None):
        """ONE WALK OF THE PERIOD UNDER A LINEAR-MULTISTEP COMPANION -- the
        plain map and gear's solved-history pair, dense or factored.

        EVERY SHOOTING ITERATION IS ITS OWN RUN (`_begin_period` /
        `_install_history` start it fresh): inherited ring buffers would make
        phi depend on the iteration, and the monodromy the derivative of
        something else.

        `opening` says where the period starts:

        * ``('plain', x_in, open_at_x0)``: ONE entering unknown.  By default
          `x(0)` is manufactured from `x_in` with one order-dropped step;
          with `open_at_x0` the unknown IS `x(0)` and the first step inside
          the period is the order-dropped opener.
        * ``('pair', x0_in, xm1_in)``: gear's SOLVED HISTORY.  The plain map
          manufactures `x(0)` from a single state, which is sound for a
          companion reaching one charge back and measurably wrong for one
          reaching two: the shooting condition constrains `x(0) = x(P)`, it
          does NOT constrain `x_in` to be the orbit's own `x(-dt)`, so `x_in`
          is an O(h^2) stand-in -- and Gear-2 reads it as a history point.
          Here BOTH `x(0)` and `x(-dt)` are unknowns, required to close as
          ``F = [x_{N-1} - x_0, x_{N-2} - x_{-1}]``, so the trajectory opens
          at full order off a history the solve is responsible for.

        What rides along the walk:

        * `dense`: the monodromy columns (`m` for the plain map, `2m` for the
          pair), propagated step by step and DROPPED -- the Newton's
          Jacobian.  The walk also records `Cvec`/`Jtvec`/`times`, as it
          always has on this path.
        * `keep`: each step's FACTORED `Jf` with its `C` and coefficients,
          for a later replay (`FactoredPeriod`) -- the matrix-free Newton's
          matvec and every small-signal consumer.  ⚠ THE FACTORISATION IS
          THE POINT: a replay runs one solve per step per Krylov iteration,
          and refactoring `Jf` each time would cost `k` times the
          factorisations.  ⚠ AND THE COST IS MEMORY: `N` factorisations and
          `N` capacitances, `2 N m^2` doubles (~800 MB at `m = 1002` and 50
          points) -- a caller with a long period and a large circuit can run
          out of memory where the dense walk merely runs slowly.
        * `want_dT`: the PERIOD column, the same recursion with the step
          size's own source term (`dh/dT = h/T` on every step, or `1` on
          the closing step alone).  It does not depend on the Newton
          direction, so a factored walk computes it once here rather than
          once per Krylov iteration.
        * `hsens` / `capture`: the state-event columns of the pair (phase B)
          and the nodes whose state and sensitivities the bordered residual
          reads.

        ⚠ ONE LINEAR SOLVE PER STEP, CHOSEN BY WHETHER THE STEP IS KEPT.  A
        kept step is factored once (`_factorise`) and every column solves
        against that; a dropped one solves through the caller's solver
        directly.  Both are an LU with partial pivoting, but numpy's `solve`
        and scipy's `lu_factor`/`lu_solve` differ in the last bits (up to
        ~3.5e-13 relative), so a dense and a factored walk of the same period
        agree to rounding, not bit for bit.

        Returns a `_PeriodWalk`; `_walk` is the one dispatch over the kinds.

        History: `doc/shooting_history.md`, `_PeriodWalks._walk_lmm`.
        """
        pair = opening[0] == 'pair'
        toolkit = self.toolkit
        m = self.cir.n - 1
        solver = self._get_linearsolver()
        eye = np.asarray(toolkit.eye(m))
        self._want_dfdh = want_dT
        Px = Pq = None
        if pair:
            _kind, x0_in, xm1_in = opening
            ## THE HISTORY IS INSTALLED, NOT SEEDED (`_install_history`): the
            ## first real step reads `q(x_0)` and `q(x_{-1})`, two genuine
            ## solved points; `_dt_last2` stays None, so the LTE estimator's
            ## opening reading is still discarded.
            self._install_history(x0_in, xm1_in, hs[0], h_prev=hs[-1])
            opened = [np.asarray(self._C_at(x0_in)),
                      np.asarray(self._C_at(xm1_in))]
            Cs = list(opened)
            x, x_prev, x0 = copy(x0_in), copy(xm1_in), x0_in
            if dense:
                ## `P_0 = [I 0]`, `P_{-1} = [0 I]` -- the two unknowns,
                ## exactly.  The plain map seeds BOTH rings with `I`, which
                ## is the flat-history assumption written into the Jacobian;
                ## here there is nothing to assume.
                zero = np.zeros((m, m))
                Px = [np.hstack((eye, zero)), np.hstack((zero, eye))]
                ## `Pq` is `d(iq_{-1})/d(x_0, x_{-1})`.  For `b = 0` the
                ## recursion never reads it and zero is right.  (A `b != 0`
                ## companion would need ``-G(x_{-1})`` in the second block --
                ## `_install_history` seeds `iq_{-1} = -(i(x_{-1}) + u)` --
                ## but only companions reaching two charges back take the
                ## pair, `_solves_history`, and Gear-2's `b` is 0.)
                Pq = np.zeros((m, 2 * m))
                ## Kept as the plain map does -- but opening EMPTY: `x_0` is
                ## an unknown here rather than the result of a step, so
                ## there is no solved `(C, Jf)` pair at it to record.
                self.Cvec, self.Jtvec = [], []
        else:
            _kind, x_in, open_at_x0 = opening
            self._begin_period(x_in)
            if open_at_x0:
                ## ⚠ NO MANUFACTURING STEP: the caller's unknown IS `x_0`.
                ## `_begin_period` has seeded both charge rings from
                ## `q(x_in)` and marked the next step `is_first_step`, so
                ## the first step INSIDE the period is order-dropped to
                ## Euler -- the L-stable opener this formulation cannot do
                ## without.  See `x0_unknown` on `solve`.
                x, x0 = copy(x_in), copy(x_in)
                C_open = np.asarray(self._C_at(x_in))
            else:
                x = self.solve_timestep(x_in, times[0], hs[0])
                x0 = copy(x)
                C_open = np.asarray(self._C)
            x_prev = None
            ## `Px[k]` is d(x_{n-k})/d(x0), `Pq` is d(iq_n)/d(x0), `Cs[k]`
            ## the capacitance of step n-k.  Two of each -- as far back as
            ## any method here reaches.  Both rings open seeded with the
            ## entering step, mirroring how the transient seeds `_qlast`
            ## with `q0` repeated.
            Cs = [copy(C_open), copy(C_open)]
            ## ⚠ THE OPENING'S COEFFICIENTS ARE THE OPENING'S.  `_coeffs` is
            ## live state and the manufacturing step is order-dropped, so it
            ## reports Euler's `(alphas, b)` -- `b = 0` -- where the loop's
            ## steps report the method's own; neither may stand in for the
            ## other (Euler is exact both ways, so a one-method test cannot
            ## tell).  ⚠ `_coeffs` DOES NOT EXIST YET when opening AT `x_0`
            ## -- no step has run -- so it is not read: `b_open = 0`, no
            ## companion current has been formed.
            a_open, b_open = ((None, 0.0) if open_at_x0 else self._coeffs)
            ## ⚠ UNLESS THE METHOD SEEDS ONE.  `theta` refuses the opener
            ## and so reads `iq_{-1}` on step one, where `_begin_run` puts
            ## `-(i(x_0) + u(t_0))` -- a function of the unknown.  See
            ## `_pq_seed_at_x0`; `None` restores the zero seed exactly.
            pq_open = self._pq_seed_at_x0(x_in) if open_at_x0 else None
            opened = (copy(C_open), a_open, b_open, pq_open)
            if dense:
                Px = [eye, eye]
                if open_at_x0:
                    ## the seed is EXACT here rather than assumed: with `x_0`
                    ## the unknown and both rings holding `q(x_0)`,
                    ## `dq_{-1}/dx_0` really IS `C`
                    Pq = pq_open if pq_open is not None else np.zeros((m, m))
                else:
                    Pq = a_open[0] * Cs[0] if b_open else np.zeros((m, m))
                ## ⚠ `C_open`, not `self._C`: on the `open_at_x0` path no
                ## step has run, so `_C` does not exist yet.
                self.Cvec = [copy(C_open)]
                self.Jtvec = [] if open_at_x0 else [copy(self._Jf)]
        if dense:
            self.times = times
        if dense or capture is not None:
            ## the nodes a staged residual reads (`capture`); on a factored
            ## walk (the matrix-free event stage) without the dense map
            self._captured = {}
        ## the period column, zero at the start: the unknowns are states,
        ## and the solve owns them
        Pt = [np.zeros(m), np.zeros(m)]
        Pqt = np.zeros(m)
        ## ⚠ EVENT COLUMNS FOR A TWO-STEP COMPANION.
        ## `hsens[j, k] = d h_j / d theta_k`; a BDF-2 step's residual depends
        ## on ITS step and on the PREVIOUS one (the 3/2 of `residual_dT`),
        ## so each column carries `dr/dh_n hsens[j] + dr/dh_{n-1} hsens[j-1]`
        ## with `dr/dh_{n-1} = (residual_dT T - residual_dh h_n) / h_{n-1}`
        ## (Euler's theorem on the homogeneous coefficients gives the total
        ## under uniform scaling; the difference is the previous step's
        ## partial), plus the source at the NEW time, `u_dot(t_{n+1}) (tau +
        ## hsens[j])`, for a driven circuit.
        self._want_event_cols = hsens is not None
        Pk = ([[np.zeros(m), np.zeros(m)] for _ in range(hsens.shape[1])]
              if hsens is not None else None)
        Pqk = ([np.zeros(m) for _ in range(hsens.shape[1])]
               if hsens is not None else None)
        _tau = np.zeros(hsens.shape[1]) if hsens is not None else None
        steps = [] if keep else None
        last = len(times) - 2
        for _j, t in enumerate(times[1:]):
            dt = hs[min(_j, len(hs) - 1)]
            x_prev = x
            x = copy(self.solve_timestep(x, t, dt))
            ## ⚠ THE COEFFICIENTS BELONG TO THE STEP, NOT TO THE RUN -- read
            ## per step and, on a kept step, stored with it: a matrix-free
            ## replay happens after the run, when `_coeffs` no longer
            ## describes the step being replayed.  (Inside the loop they are
            ## constant for every method in this tree, so no test catches a
            ## post-run snapshot; a variable-order method would make that
            ## failure silent.)
            alphas, b = self._coeffs
            Jf = np.asarray(self._Jf)
            C_new = np.asarray(self._C).copy()
            if dense:
                self.Cvec.append(copy(self._C))
                self.Jtvec.append(copy(self._Jf))
            if keep:
                lu = self._factorise(Jf)
                steps.append((lu, C_new, alphas, b))
                solve = lu.solve
            else:
                ## ⚠ THROUGH THE CALLER'S SOLVER, not `toolkit.linearsolver`
                ## (`DenseSolver` IS that, so the default is unchanged): a
                ## sparse matrix-free path must not be measured against a
                ## hard-wired dense baseline.
                def solve(S, _Jf=Jf):
                    return solver.solve(_Jf, S, toolkit)
            ## ONE RECURSION FOR EVERY METHOD AND EVERY COLUMN
            ## (`_lmm_recursion`): the monodromy, the period column and the
            ## event columns differ only in the forcing each step adds.
            if dense:
                Px_new, Pq = _lmm_recursion(Px, Cs, Pq, C_new, alphas, b,
                                            solve)
            if Pk is not None:
                ## ⚠ `residual_dh` is Fang's `p`: the coefficients' partial
                ## PLUS `du/dt`, because the residual sits at `t_n + h_n` and
                ## that time moves with `h_n`; `residual_dT` is the
                ## coefficients' total under uniform scaling.  So the
                ## previous step's coefficient partial is `(residual_dT -
                ## (residual_dh - u_dot) h_n) / h_{n-1}`, and the source's
                ## motion enters once: `residual_dh`'s `u_dot h_n` part
                ## through this step's weight, and `u_dot tau_n` for the
                ## shift of the step's START -- never `u_dot (tau + w)` on
                ## top of `residual_dh`, which counts it twice.
                dr_dhn_full = np.asarray(self._dfdh, dtype=float).ravel()
                Ud = np.delete(np.asarray(self.cir.dudt(
                    float(t), analysis=self.par.analysis), dtype=float),
                    self.irefnode)
                dr_dhn_iq = dr_dhn_full - Ud
                dt_prev = float(hs[_j - 1]) if _j > 0 else float(hs[-1])
                dr_dhprev = (np.asarray(self._dfdT, dtype=float).ravel()
                             - dr_dhn_iq * float(dt)) / dt_prev
                hprev = hsens[_j - 1] if _j > 0 else hsens[-1]
                for k in range(len(Pk)):
                    if not pair:
                        ## ⚠ ON THE PLAIN MAP THE NEXT STEP MAY READ `Pq` --
                        ## trapezoidal's `b != 0` companion carries `iq`, and
                        ## its first step is an order-dropped Euler whose `Pq`
                        ## the second reads -- and the source's motion is a
                        ## CURRENT, not a charge: it enters the solve
                        ## (`source`) and not the carried companion
                        ## sensitivity.  (Gear's pair has `b = 0` on every
                        ## step and never reads `Pq`, which is why the lumped
                        ## forcing below is exact there; on trap's columns on
                        ## a driven PWM loop it was 0.3-380x off, 2026-09-24.)
                        ## ⚠ AND NO PREVIOUS-STEP PARTIAL: a one-step
                        ## companion's coefficients depend on `h_n` alone --
                        ## the previous step enters through the CARRIED state
                        ## `(x, iq)` -- so `dr/dh_{n-1}` is zero by structure.
                        ## The Euler-theorem estimate above assumes
                        ## coefficients homogeneous in `h`, which theta's
                        ## `1/2 + c h` is not (its columns read up to 358x off).
                        Pk_new, Pqk[k] = _lmm_recursion(
                            Pk[k], Cs, Pqk[k], C_new, alphas, b, solve,
                            source=Ud * float(hsens[_j, k] + _tau[k]),
                            forcing=(dr_dhn_iq * float(hsens[_j, k]),))
                    else:
                        Pk_new, Pqk[k] = _lmm_recursion(
                            Pk[k], Cs, Pqk[k], C_new, alphas, b, solve,
                            forcing=(dr_dhn_full * float(hsens[_j, k]),
                                     dr_dhprev * float(hprev[k]),
                                     Ud * float(_tau[k])))
                    Pk[k] = [Pk_new, Pk[k][0]]
                _tau = _tau + hsens[_j]
                if capture is not None and (_j + 1) in capture:
                    self._captured[_j + 1] = (copy(x),
                                              (np.asarray(Px_new).copy()
                                               if dense else None),
                                              [pk[0].copy() for pk in Pk])
            if want_dT:
                ## Every step scales together, so `dh/dT = h/T`, and `df/dh`
                ## at fixed solution is Fang's `p` (`residual_dT` for the
                ## total under uniform scaling -- for Gear-2 the partial is
                ## 3/2 of it).  For an AUTONOMOUS circuit its `du/dt` half
                ## vanishes, which is what makes solving for the period
                ## tractable at all.  'closing': only the closing step's
                ## length depends on `T`, with `dh/dT = 1`, so the term
                ## appears once, undivided, on the last step.
                if self._period_column == 'closing':
                    fT = ((np.asarray(self._dfdh).ravel(),) if _j == last
                          else ())
                else:
                    fT = (np.asarray(self._dfdT).ravel() / T,)
                Pt_new, Pqt = _lmm_recursion(Pt, Cs, Pqt, C_new, alphas, b,
                                             solve, forcing=fT)
                Pt = [Pt_new, Pt[0]]
            if dense:
                Px = [Px_new, Px[0]]
            Cs = [C_new, Cs[0]]
        self._want_dfdh = False
        self._want_event_cols = False
        return _PeriodWalk(kind=opening[0],
                           fp_kind='solved_history' if pair else 'plain',
                           x0=x0, x_end=x, x_prev=x_prev, P=Px,
                           Pt=Pt if want_dT else None, Pk=Pk, steps=steps,
                           opening=opened,
                           open_at_x0=(not pair) and bool(opening[2]))

    ## ------------------------------------------------------------------
    ## The Runge-Kutta STAGE family (Radau IIA, TR-BDF2, ESDIRK) --
    ## self-starting, tableau-generic over any number of stages.  One
    ## `_StageStep` per step carries the factored stage system in the
    ## structure its tableau calls for (the coupled block for a fully
    ## implicit one, per-stage factors for a lower-triangular one: an
    ## explicit first stage makes the coupled block singular on a DAE, and
    ## the sequential recursion never forms it).  See
    ## doc/integrator_architecture_260906.md.
    ## ------------------------------------------------------------------

    def _stage_step(self, xn, h, tab, coupled):
        """The factored stage system of the step just taken from `xn` -- its
        stage states read off the inner transient -- as a `_StageStep`, and
        those stage states (reduced).  ``C(Y_i)`` and ``G(Y_j)`` are at the
        DISTINCT stage points, which is why `_C_at`/`_G_at` exist rather
        than a stored `Geq`."""
        import scipy.linalg as sla
        A, b, c = tab
        s = A.shape[0]
        iref = self.irefnode
        Ys = [self.toolkit.concatenate((yf[:iref], yf[iref + 1:]))
              for yf in self._transient()._rk_Y]
        Cn = np.asarray(self._C_at(xn))
        m = Cn.shape[0]
        if coupled:
            Cs = [np.asarray(self._C_at(y)) for y in Ys]
        Gs = [np.asarray(self._G_at(y)) for y in Ys]
        if coupled:
            Jb = np.zeros((s * m, s * m))
            for i in range(s):
                for jj in range(s):
                    blk = h * A[i, jj] * Gs[jj]
                    if i == jj:
                        blk = Cs[i] + blk
                    Jb[i * m:(i + 1) * m, jj * m:(jj + 1) * m] = blk
            return _StageStep(Cn, Gs, h, A, b, c, lu=sla.lu_factor(Jb)), Ys
        Kf = []
        for i in range(s):
            if abs(A[i, i]) < 1e-14:
                if i != 0:
                    raise NotImplementedError(
                        'PSS: an explicit stage other than the first is not '
                        'supported (its solve would need C^-1, singular on a '
                        'DAE).')
                Kf.append(None)
            else:
                Ci = np.asarray(self._C_at(Ys[i]))
                Kf.append(self._factorise(Ci + h * A[i, i] * Gs[i]))
        return _StageStep(Cn, Gs, h, A, b, c, Kf=Kf), Ys

    def _glm_period_blocks(self, x_in, times, hs):
        """One period under a Nordsieck GLM, collecting per-step factors.

        Drives `Transient._solve_timestep_glm` with the Nordsieck vector fed in
        explicitly, so the startup runs ONCE at `t = 0` and every later step
        continues the multivalue state (no seam inside the period).  Returns
        ``(steps, xs, Q0, Qend, x_end, trace0)``, each step a `_GLMStep`
        record (its fields are documented there): the stage factors and
        conductances, the tableau, the stage derivatives, abscissae and
        states, the rescale ``rho = h / h_prev`` the step applied and the
        vector it entered with, the state it started from, the next node's
        startup where the next step RESTARTS on growth, and whether this
        step entered through such a restart.  `trace0` is the period's own
        startup (`Transient._glm_startup_trace` after the first step: a
        restart later in the period overwrites the transient's).
        """
        tr = self._transient()
        iref = self.irefnode
        self._want_dfdh = False
        self._want_lte = False
        self._begin_period(x_in)          # sets up the transient's integrator
        integ = tr.base_integrator
        A, U, B, V, c, p = integ.tableau()
        s = A.shape[0]
        x = copy(x_in)
        tr._glm_Q = None                      # force the startup at t = 0
        ## ⚠ AND THE ENTRY SLOT.  It is keyed by the time the step STARTED
        ## at, so the previous shooting iteration's first step left one at
        ## t = 0 -- exactly the time this period's first step asks for.  Left
        ## behind it is picked up in preference to a fresh startup and the
        ## period map silently reads the LAST iteration's Nordsieck vector.
        tr._glm_Q_at_entry = None
        tr._glm_prev = None                   # no stage predictor across the seam
        tr._pred_reset()                      # nor its node history: the period
                                              # seam is a DISCONTINUITY in the
                                              # trajectory a predictor fits
        steps, xs = [], []
        Q0 = None
        trace0 = None
        for _j, t in enumerate(times[1:]):
            h = hs[min(_j, len(hs) - 1)]
            xn = x
            x = copy(self.solve_timestep(xn, t, h))
            restarted = bool(getattr(tr, '_glm_restarted', False)) and _j > 0
            if _j == 0:
                trace0 = tr._glm_startup_trace
            elif restarted:
                ## the step before ends on this node's startup
                ## (`_GLMStep.forward`)
                steps[-1].restart_out = self._glm_startup_linearisation(
                    tr._glm_startup_trace)
            Qn = np.asarray(tr._glm_Q[0], dtype=float)
            if Q0 is None:
                ## the vector the first step actually entered with, i.e. what
                ## the startup produced from `x_in`.  ⚠ REDUCED: the transient
                ## carries the Nordsieck vector at FULL width (its rows are
                ## charges, reference row included); every sensitivity here is
                ## on the reduced state, so the seed must be too.
                Q0 = np.asarray(tr._glm_Q_in, dtype=float)[:, [i for i in
                                                               range(self.cir.n)
                                                               if i != iref]]
            Ys = [self.toolkit.concatenate((yf[:iref], yf[iref + 1:]))
                  for yf in tr._rk_Y]
            Gs = [np.asarray(self._G_at(Ys[i])) for i in range(s)]
            Kfacs = [self._factorise(np.asarray(self._C_at(Ys[i]))
                                     + h * A[i, i] * Gs[i]) for i in range(s)]
            Ks = [np.asarray(self.toolkit.concatenate((kf[:iref], kf[iref + 1:])),
                             dtype=float) for kf in tr._rk_K]
            ## ⚠ THE RESCALE IS PART OF THE MAP: where the step changes, the
            ## transient scales the Nordsieck vector by `rho^k` before the
            ## step, and the sensitivities must be scaled with it.  Missing,
            ## the map on a 3:1 grid was 0.7 % (glm2) / 2 % (glm3) off its
            ## finite difference, exact on a uniform one.
            Qin = np.delete(np.asarray(tr._glm_Q_in, dtype=float), iref, axis=1)
            steps.append(_GLMStep(Kfacs, Gs, float(h), A, U, B, V, Ks,
                                  rho=float(getattr(tr, '_glm_rho', 1.0)),
                                  Qin=Qin, x_in=np.asarray(xn, dtype=float),
                                  restart_out=None, restarted=restarted,
                                  c=np.asarray(c, dtype=float),
                                  Ys=[np.array(yf, dtype=float)
                                      for yf in tr._rk_Y]))
            xs.append(np.asarray(x, dtype=float))
        return (steps, xs, Q0, np.asarray(tr._glm_Q[0], dtype=float), x,
                trace0)

    @staticmethod
    def _glm_propagate(steps, P, T=None, closing=False):
        """The Nordsieck sensitivity recursion, one period.

        ``P`` is a list of `r` blocks ``dQ_k/d(unknown)``; each step maps it by

            D_i  = K_i^{-1} ( sum_j U_ij P_j - h sum_{j<i} A_ij G_j D_j )
            P'_k = sum_j V_kj P_j - h sum_i B_ki G_i D_i

        -- the DIRK recursion with the entering ``C_n`` term replaced by the
        Nordsieck combination ``U P``.  Returns ``(P_out, D_last)``; ``D_last``
        is ``d x_N / d(unknown)`` because the method is stiffly accurate
        (``x_N`` is the last stage).

        With ``T`` given the recursion carries the PERIOD column instead: the
        grid is ``h_j = frac_j T`` so ``dh/dT = h/T``, and on an autonomous
        circuit the stage derivatives carry no time of their own, which adds
        ``(h/T) sum_j A_ij K_j`` to each stage's right-hand side and
        ``(h/T) sum_i B_ki K_i`` to each output component.
        """
        D = None
        _nsteps = len(steps)
        def _dhdT(i):
            return ((1.0 if i == _nsteps - 1 else 0.0) if closing
                    else steps[i].h / T)
        for _si, rec in enumerate(steps):
            ## 'closing': dh/dT = 1 on the last step, 0 elsewhere.  The
            ## explicit `h = frac T` in the stage: for an AUTONOMOUS circuit
            ## `K_j = -i(Y_j)` carries no time of its own, so the only new
            ## term is the stage sum (see `_GLMStep.forward`)
            fT = None if T is None else _dhdT(_si)
            ## ⚠ AND UNDER 'closing' THE LAST STEP'S RESCALE MOVES WITH `T`
            ## (`rho = h_N / h_{N-1}`, `d rho / dT = rho / h_N`), even on a
            ## uniform grid where `rho` is 1; proportionally scaled steps keep
            ## every `rho` fixed.  (A step entered by a restart has none.)
            drho = (rec.rho / rec.h
                    if (T is not None and closing and _si == _nsteps - 1
                        and not rec.restarted)
                    else None)
            ## a step ending on a restart: the startup moves with the next
            ## step's length
            nxt = ((_dhdT(_si + 1), 0.0)
                   if (T is not None and rec.restart_out is not None
                       and _si + 1 < _nsteps) else None)
            P, D = rec.forward(P, fT, drho, None, nxt)
        return P, (D[-1] if D is not None else None)

    def _walk_glm(self, x_in, T, times, hs, dense=True, keep=False,
                  want_dT=False, hsens=None, capture=None):
        """ONE WALK OF THE PERIOD UNDER A NORDSIECK GLM: `_glm_period_blocks`,
        then, when `dense`, the sensitivity recursion (`_glm_propagate`) for
        the map and, with `want_dT`, its period column.  The blocks are always
        collected -- the dense map is propagated FROM them -- so `keep` only
        decides whether they are handed back (a `FactoredPeriod` of kind
        'glm', on the Nordsieck state, `width` ``r*m``).

        The map shot on is ``x_0 -> x_N``: the Nordsieck vector is built from
        ``x_0`` by `Transient._glm_startup` at the top of the period and
        propagated by the method to the end.  The recursion is seeded with
        the startup LINEARISED (`_GLMStartup`: its p Radau substeps and its
        interpolant, ``dQ_k/dx_0``), so the map on ``x`` is exact, and the
        period column carries the substeps' own motion with `T`.  (Until
        2026-09-24 only ``dQ_0/dx_0 = C(x_0)`` was carried: the Jacobian was
        approximate, and the factored map not the Newton's.)

        ⚠ EVENT COLUMNS (`hsens`, ``d h_j / d theta_k``; `capture`, the nodes
        the bordered residual reads) run through the same step with the
        three ways a step's length enters it: its own `h` in the stage and
        output rows (``dh sum A K``, ``dh sum B K``), the entering rescale
        ``rho = h_j / h_{j-1}`` (``d rho = (dh_j - rho dh_{j-1}) / h_{j-1}``)
        and, on a driven circuit, the source moving with the stage times
        (``u_dot (tau_n + c_i dh)``, `tau_n` the shift of the step's start);
        at node 0 the startup's own substeps move (`_GLMStartup.dh`).

        History: `doc/shooting_history.md`, `_PeriodWalks._walk_glm`.
        """
        steps, _xs, Q0, _Qend, x_end, trace0 = self._glm_period_blocks(
            x_in, times, hs)
        m = self.cir.n - 1
        r = len(Q0)
        su = self._glm_startup_linearisation(trace0)
        Mx = Mt = None
        Pk_end = None
        if hsens is not None:
            Mx, Pk_end = self._glm_event_columns(steps, _xs, su, times,
                                                 hsens, capture, dense=dense)
        elif dense:
            _Pout, Mx = self._glm_propagate(steps, su.matrix())
        if want_dT:
            ## the period column: the startup's substeps scale with `T`
            ## (`_GLMStartup.dT`) -- unless the convention is 'closing',
            ## where only the LAST step's length moves and the startup, at
            ## the first, does not
            closing = self._period_column == 'closing'
            Pt = ([np.zeros(m) for _ in range(r)] if closing
                  else su.dT(float(T)))
            _Ptout, Mt = self._glm_propagate(steps, Pt, T=float(T),
                                             closing=closing)
            Mt = np.asarray(Mt).ravel()
        return _PeriodWalk(kind='glm', fp_kind='glm',
                           x0=np.asarray(x_in, dtype=float),
                           x_end=np.asarray(x_end, dtype=float), P=Mx, Pt=Mt,
                           Pk=Pk_end, steps=steps if keep else None,
                           width=r * m, startup=su)

    def _glm_event_columns(self, steps, xs, su, times, hsens, capture,
                           dense=True):
        """The monodromy on the state and the event columns of one GLM
        period, step by step (`_walk_glm`), capturing ``(x_j, dx_j/dx_0,
        [dx_j/dtheta_k])`` at the nodes in `capture` into `_captured`.
        Returns ``(M, [dx_N/dtheta_k])``; without `dense`, the columns alone
        (`M` and the captured maps None: the matrix-free event stage)."""
        iref = self.irefnode
        integ = self._transient().base_integrator
        c = np.asarray(integ.tableau()[4], dtype=float)
        K = hsens.shape[1]
        P = su.matrix() if dense else None
        Pk = [su.dh(float(hsens[0, k])) for k in range(K)]
        tau = np.zeros(K)
        self._captured = {}
        D = None
        Dk = [None] * K
        for j, rec in enumerate(steps):
            if dense:
                P, D = rec.forward(P)
            h = rec.h
            t0 = float(times[j])
            Ud = [np.delete(np.asarray(self.cir.dudt(t0 + float(ci) * h,
                                                      analysis=self.par.analysis),
                                       dtype=float), iref) for ci in c]
            for k in range(K):
                w = float(hsens[j, k])
                drho = ((w - rec.rho * float(hsens[j - 1, k])) / steps[j - 1].h
                        if (j > 0 and not rec.restarted) else None)
                src = [Ud[i] * (tau[k] + float(c[i]) * w) for i in range(len(c))]
                nxt = ((float(hsens[j + 1, k]), float(tau[k]) + w)
                       if (rec.restart_out is not None and j + 1 < len(steps)) else None)
                Pk[k], Dk[k] = rec.forward(Pk[k], fT=w, drho=drho, src=src,
                                           nxt=nxt)
            tau = tau + hsens[j]
            if capture is not None and (j + 1) in capture:
                self._captured[j + 1] = (
                    np.asarray(xs[j], dtype=float),
                    np.asarray(D[-1]).copy() if dense else None,
                    [np.asarray(d[-1]).copy() for d in Dk])
        return (np.asarray(D[-1], dtype=float) if dense else None,
                [np.asarray(d[-1], dtype=float).ravel() for d in Dk])

    def _glm_node_startups(self, fp):
        """The startup AT EACH NODE of a GLM period, linearised: ``S_j``, how
        the Nordsieck vector a startup would build at ``t_j`` from ``x_j``
        (with that node's step, ``h_j``) moves with ``x_j``.  Node 0's is the
        period's own (`fp.startup`).  What turns a costate on the Nordsieck
        vector into one on the STATE: ``v_j = S_j^T lambda_j`` (`ppv`'s
        samples) -- the first block alone, `dphi/dQ_0` with the higher
        components HELD, is the inconsistent object gear's pair also had.
        Costs one startup (p Radau substeps) per node."""
        tr = self._transient()
        iref = self.irefnode
        times = np.asarray(fp.times, dtype=float)
        out = [fp.startup]
        saved = getattr(tr, '_glm_startup_trace', None)
        try:
            for j in range(1, len(fp.steps)):
                rec = fp.steps[j]
                xf = np.insert(np.asarray(rec.x_in, dtype=float), iref, 0.0)
                tr._glm_startup(float(times[j]), xf, float(rec.h))
                out.append(self._glm_startup_linearisation(
                    tr._glm_startup_trace))
        finally:
            tr._glm_startup_trace = saved
        return out

    def _glm_startup_linearisation(self, trace=None):
        """The startup of the period just walked (or the one `trace`
        records), linearised (`_GLMStartup`), from what
        `Transient._glm_startup` recorded -- the substep states and stages --
        with the point evaluations the rest of the map uses (`_C_at`,
        `_G_at`, `_k_at`)."""
        from math import factorial
        from scipy.linalg import lu_factor
        from pycircuit.circuit.integrator import RadauIIA3Integrator
        tn, hs, p, xs, Ys = (self._transient()._glm_startup_trace
                             if trace is None else trace)
        A = np.array(RadauIIA3Integrator.A, dtype=float)
        c = np.array(RadauIIA3Integrator.C, dtype=float)
        iref = self.irefnode
        m = self.cir.n - 1

        def red(v):
            return np.delete(np.asarray(v, dtype=float), iref)
        Cx = [np.asarray(self._C_at(red(x)), dtype=float) for x in xs]
        lus, fT, Ud = [], [], []
        for j in range(1, p + 1):
            Yj = [red(y) for y in Ys[j - 1]]
            Ci = [np.asarray(self._C_at(y), dtype=float) for y in Yj]
            Gi = [np.asarray(self._G_at(y), dtype=float) for y in Yj]
            J = np.zeros((3 * m, 3 * m))
            for i in range(3):
                for l_ in range(3):
                    blk = hs * A[i, l_] * Gi[l_]
                    if i == l_:
                        blk = blk + Ci[i]
                    J[i * m:(i + 1) * m, l_ * m:(l_ + 1) * m] = blk
            lus.append(lu_factor(J))
            ts = [tn + (j - 1) * hs + c[l_] * hs for l_ in range(3)]
            K = [np.asarray(self._k_at(Yj[l_], ts[l_]), dtype=float)
                 for l_ in range(3)]
            fT.append(np.concatenate([sum(A[i, l_] * K[l_] for l_ in range(3))
                                      for i in range(3)]))
            Ud.append([red(self.cir.dudt(ts[l_], analysis=self.par.analysis))
                       for l_ in range(3)])
        Vd = np.array([[float(k) ** jj for jj in range(p + 1)]
                       for k in range(p + 1)])
        Vi = np.linalg.solve(Vd, np.eye(p + 1))
        W = np.array([[p ** k * factorial(k) * Vi[k, jj] for jj in range(p + 1)]
                      for k in range(p + 1)])
        W[0] = 0.0
        W[0, 0] = 1.0
        return _GLMStartup(Cx, lus, fT, W, hs, A=A, c=c, Ud=Ud, tn=float(tn),
                           Ys=[[np.array(y, dtype=float) for y in Yj]
                               for Yj in Ys])

    def _walk_stage(self, x_in, T, times, hs, dense=True, keep=False,
                    want_dT=False, hsens=None, capture=None):
        """ONE WALK OF THE PERIOD UNDER A RUNGE-KUTTA STAGE METHOD (Radau
        IIA, TR-BDF2, ESDIRK), dense or factored.
        Self-starting: `x_in` IS `x_0`, so there is no opener seam and the
        map keeps the method's order round the whole period.

        Every step is one `_StageStep` -- the stage system ``J Z = B`` with
        ``J[i][j] = delta_ij C(Y_i) + h A_ij G(Y_j)``, coupled for a fully
        implicit tableau, stage by stage for a lower-triangular one.  What
        rides along, as for `_walk_lmm`: `dense` propagates the monodromy
        `P = dx/dx0` through it and drops it; `keep` keeps it for a replay
        (`FactoredPeriod`); `want_dT` and `hsens` propagate the period and
        event columns through the SAME solve, their right-hand sides
        differing only in their forcing:

            P:   [C_n P]_i
            Pt:  [C_n Pt + (dh/dT) S_i]_i,       S_i = sum_j A_ij K_j
            Pk:  [C_n Pk + w S_i - h U_i]_i

        and the last stage is the step's result (stiff accuracy).

        ⚠ THE PERIOD COLUMN IS TRACTABLE ONLY BECAUSE THE CIRCUIT IS
        AUTONOMOUS.  With `T` unknown the grid rebuilds as ``h_j = frac_j T``,
        so ``dh_j/dT = h_j/T`` and the stage derivative ``K_j = -i(Y_j)`` has
        no explicit time dependence; differentiating the stage residuals
        ``F_i = q(Y_i) - q(x_n) - h sum_j A_ij K_j`` w.r.t. `T` gives the `Pt`
        row above.  ('closing': only the last step's length depends on `T`.)
        Finite-difference checked before use (the dT column is easy to get
        wrong).

        ⚠ EVENT COLUMNS: `hsens[j, k] = d h_j / d theta_k` for the
        state-event unknowns, propagated by the same stage algebra as the
        period column with the per-step weight taken from the matrix instead
        of `h/T`; `capture` names the nodes whose state and sensitivities the
        bordered residual reads.  ⚠ A DRIVEN CIRCUIT'S SOURCES MOVE WITH THE
        GRID: an event column shifts the TIMES the stages are evaluated at,
        so `f = -(i + u(t))` changes by `-u_dot . dt_stage`, `dt_stage =
        tau_n + c_i dh` with `tau_n` the shift of the step's start (the `U_i`
        term).  Without it the column is badly wrong, down to the sign of the
        event row's derivative.

        Returns a `_PeriodWalk`; `_factored_self_starting` keeps its
        steps.

        History: `doc/shooting_history.md`, `_PeriodWalks._walk_stage`.
        """
        from pycircuit.circuit.integrator import RungeKuttaIntegrator
        toolkit = self.toolkit
        m = self.cir.n - 1
        self._want_dfdh = False
        self._want_lte = False
        self._begin_period(x_in)
        integ = self._transient().base_integrator
        if not isinstance(integ, RungeKuttaIntegrator):
            raise ValueError('_walk_stage needs a Runge-Kutta stage inner '
                             'integrator, got %r' % (integ,))
        tab = _butcher(integ)
        Amat = tab[0]
        s = Amat.shape[0]
        coupled = integ.is_fully_implicit()
        iref = self.irefnode
        x = copy(x_in)
        x0 = copy(x_in)
        x_prev = copy(x_in)
        P = np.asarray(toolkit.eye(m), dtype=float) if dense else None
        Pt = np.zeros(m)
        Pk = ([np.zeros(m) for _ in range(hsens.shape[1])]
              if hsens is not None else None)
        if dense or capture is not None:
            self._captured = {}
        _cabs = tab[2] if Pk is not None else None
        _tau = np.zeros(hsens.shape[1]) if Pk is not None else None
        steps = [] if keep else None
        for _j, t in enumerate(times[1:]):
            h = hs[min(_j, len(hs) - 1)]
            xn = x
            x = copy(self.solve_timestep(xn, t, h))
            x_prev = xn
            st, Ys = self._stage_step(xn, h, tab, coupled)
            if keep:
                steps.append(st)
            if dense:
                P = st.solve(P)
            if Pk is not None:
                _t0 = float(times[_j])
                Ks = [np.asarray(self._k_at(Ys[jj], _t0 + float(_cabs[jj]) * h))
                      for jj in range(s)]
                Ss = [sum(Amat[i, jj] * Ks[jj] for jj in range(s)) for i in range(s)]
                Ud = [np.delete(np.asarray(self.cir.dudt(_t0 + float(_cabs[jj]) * h,
                                                          analysis=self.par.analysis),
                                           dtype=float), iref)
                      for jj in range(s)]
                for k in range(len(Pk)):
                    w = float(hsens[_j, k])
                    forcing = []
                    for i in range(s):
                        Ui = sum(Amat[i, jj] * Ud[jj] * (_tau[k] + float(_cabs[jj]) * w)
                                 for jj in range(s))
                        forcing.append((w * Ss[i], -(h * Ui)))
                    Pk[k] = st.solve(Pk[k], forcing)
                _tau = _tau + hsens[_j]
                if capture is not None and (_j + 1) in capture:
                    self._captured[_j + 1] = (copy(x),
                                              None if P is None else P.copy(),
                                              [pk.copy() for pk in Pk])
            if want_dT:
                Ks = [np.asarray(self._k_at(y)) for y in Ys]
                ## 'closing': only the last step's length depends on T
                _dhdT = ((1.0 if _j == len(times) - 2 else 0.0)
                         if self._period_column == 'closing' else h / float(T))
                forcing = [(_dhdT * sum(Amat[i, jj] * Ks[jj] for jj in range(s)),)
                           for i in range(s)]
                Pt = st.solve(Pt, forcing)
        self._want_dfdh = False
        return _PeriodWalk(kind='stage', fp_kind='full' if coupled else 'dirk',
                           x0=x0, x_end=x, x_prev=x_prev, P=P,
                           Pt=Pt if want_dT else None, Pk=Pk, steps=steps)
