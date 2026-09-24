"""The period walks: one per step family (`_walk_lmm`, `_walk_stage`,
`_walk_glm`), dense and/or factored, and their views.
"""
from copy import copy
import numpy as np
from ._factored import _PeriodWalk
from ._steps import _StageStep
from ._steps import _butcher
from ._steps import _glm_step
from ._steps import _lmm_recursion


class _PeriodWalks(object):
    """The period walks: one per step family (`_walk_lmm`, `_walk_stage`,
    `_walk_glm`), dense and/or factored, and their views.  A theme of `PSS`
    (see `pss.py`)."""

    def _step_sensitivity(self, Px, Cs, Pq, Jf, C_new, solve=None,
                          coeffs=None, source=None):
        """One step of the sensitivity recursion, for ANY seed width.

        ONE RECURSION FOR EVERY METHOD (and now for either formulation).
        Each writes its companion as `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`,
        so differentiating the step gives

            S    = sum_{k>=1} a_k C_{n-k} P_{n-k} + b Pq
            P_n  = -Jf_n^-1 S
            Pq_n = a_0 C_n P_n + S

        `P` is `d x_j / d(unknowns)`: one block wide in the plain
        formulation, two when the history is solved for.  Nothing here
        depends on that
        width, which is why the two systems share this and not a copy.

        ⚠ A SOLVE, NOT AN INVERSE (stage 11).  `inv(Jf) @ ...` formed a dense
        inverse per timestep per iteration and squared the condition number
        it then multiplied through.

        `solve` overrides how that solve is taken, and exists so the
        MATRIX-FREE path (item 6) can hand in a PRE-FACTORED `Jf` without
        this recursion being copied.  It is the same recursion either way,
        which is the point: a second copy would be a second thing to get
        wrong, and this one is already shared by every method and both
        formulations.  Without it, matrix-free would refactor `Jf` once per
        step PER KRYLOV ITERATION -- `k` times the factorisations the dense
        path takes, which is worse than the problem it set out to fix.
        """
        ## `coeffs` overrides the LIVE `_coeffs` for the same reason `solve`
        ## overrides the solve: a matrix-free replay happens after the run
        ## that produced the steps, when `_coeffs` no longer describes the
        ## step being replayed.  See `_traverse_factored_plain`.
        ## ⚠ `source` ENTERS THE SOLVE AND NOT `Pq`, and the asymmetry is
        ## the physics rather than a convenience.  A small-signal source
        ## appears in the step's residual -- `Jf dx + S + du = 0` -- but NOT
        ## in the companion, because `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`
        ## is built from CHARGES, and an injected current is not one.
        ## Adding it to `Pq` as well would feed a fictitious charge forward
        ## into every later step, and the error would grow along the period
        ## rather than announce itself.
        ##
        ## This is what makes PAC share the recursion instead of copying it:
        ## the homogeneous propagation (`source=None`) is the monodromy and
        ## the driven one is the forced response, and they differ by this
        ## one term.
        alphas, b = self._coeffs if coeffs is None else coeffs
        if solve is None:
            ## ⚠ THROUGH THE CALLER'S SOLVER, not `toolkit.linearsolver`.
            ## This is the DENSE propagation -- the thing matrix-free is
            ## measured against -- and it used to be hardcoded to the
            ## toolkit, so `linearsolver=SuperLUSolver()` reached the inner
            ## Newton (once forwarded) and never the propagation.  Comparing
            ## a sparse matrix-free path against a dense baseline would have
            ## flattered it; both sides go through the same strategy now.
            ## `DenseSolver` IS `toolkit.linearsolver`, so the default is
            ## unchanged.
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
                                  want_dT=want_dT)
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
        plain map and gear's solved-history pair, dense or factored
        (2026-09-23: it was four walks, `_traverse`,
        `_traverse_solved_history`, `_traverse_factored_plain` and
        `_traverse_factored`, each with its own copy of the period column).

        EVERY SHOOTING ITERATION IS ITS OWN RUN.  phi must be a function of
        its arguments alone; if iteration k+1 inherited the ring buffers
        iteration k ended with, the period map would depend on which
        iteration it was and the monodromy would be the derivative of
        something else.  `_begin_period` / `_install_history` make the run
        start fresh.

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
        and scipy's `lu_factor`/`lu_solve` differ in the last bits (measured
        on this box: 140 of 360 random systems, up to 3.5e-13 relative), so
        each walk keeps the arithmetic it had -- the merge moved no answer.

        Returns a `_PeriodWalk`; the `_traverse*` names are its views.
        """
        pair = opening[0] == 'pair'
        if hsens is not None and not pair:
            raise NotImplementedError(
                'PSS: state-event columns are built for the stage methods '
                "and gear's solved-history pair, not the plain map")
        toolkit = self.toolkit
        m = self.cir.n - 1
        solver = self._get_linearsolver()
        eye = np.asarray(toolkit.eye(m))
        self._want_dfdh = want_dT
        Px = Pq = None
        if pair:
            _kind, x0_in, xm1_in = opening
            ## THE HISTORY IS INSTALLED, NOT SEEDED.  `_begin_run(x_{-1})`
            ## opens the rings on the earlier point and the push puts `x_0` in
            ## front of it, so the first real step reads `q(x_0)` and
            ## `q(x_{-1})` -- two genuine solved points.  The flags then say
            ## what is true of them: a step of `dt` has been taken
            ## (`_dt_last`), the run is no longer opening (`_is_first_step`,
            ## `_no_history`), and `_dt_last2` stays None because the THIRD
            ## charge is still `q(x_{-1})` repeated -- the LTE estimator
            ## differences three, so its opening reading remains unsound and
            ## the report goes on discarding it.
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
            ## steps report the method's own (trapezoidal opens at
            ## `((49000, -49000), 0.0)` and runs at `((98000, -98000),
            ## -1.0)`).  Reading them once for the whole run put the period
            ## column 40-50 % out for trap and gear; reading the loop's for
            ## the opening made `Pq` non-zero where it is zero, a 100 %
            ## error for trap.  Euler was exact both ways and would have
            ## passed a one-method test.  ⚠ `_coeffs` DOES NOT EXIST YET when
            ## opening AT `x_0` -- no step has run -- so it is not read:
            ## `b_open = 0`, no companion current has been formed.
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
            self._captured = {}
        ## the period column, zero at the start: the unknowns are states,
        ## and the solve owns them
        Pt = [np.zeros(m), np.zeros(m)]
        Pqt = np.zeros(m)
        ## ⚠ EVENT COLUMNS FOR A TWO-STEP COMPANION (2026-09-22, phase B).
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
            ## constant for every method in this tree, so storing them is
            ## belt-and-braces today -- a mutation replacing them with a
            ## post-run snapshot does NOT fail the tests -- but a
            ## variable-order method would make that failure silent.)
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
                ## shift of the step's START.  Counting `u_dot (tau + w)` on
                ## top of `residual_dh` -- twice -- read the node after a
                ## landed event 153 % off (2026-09-22).
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
                    Pk_new, Pqk[k] = _lmm_recursion(
                        Pk[k], Cs, Pqk[k], C_new, alphas, b, solve,
                        forcing=(dr_dhn_full * float(hsens[_j, k]),
                                 dr_dhprev * float(hprev[k]),
                                 Ud * float(_tau[k])))
                    Pk[k] = [Pk_new, Pk[k][0]]
                _tau = _tau + hsens[_j]
                if capture is not None and (_j + 1) in capture:
                    self._captured[_j + 1] = (copy(x),
                                              np.asarray(Px_new).copy(),
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

    def _traverse(self, x_in, T, times, hs, want_dT, open_at_x0=False):
        """The PLAIN map, dense: ``(x0, x_end, dx_end/dx0, dx_end/dT)`` --
        the last only when asked for.  A view of `_walk_lmm`."""
        w = self._walk_lmm(('plain', x_in, open_at_x0), times, hs, T=T,
                           want_dT=want_dT)
        ## Kept for the autonomous check after the solve: its spectrum is
        ## the only place a free period announces itself.
        self._monodromy = w.P[0]
        return w.x0, w.x_end, w.P[0], (w.Pt[0] if want_dT else None)

    def _traverse_solved_history(self, x0_in, xm1_in, times, hs,
                                 T=None, want_dT=False, hsens=None,
                                 capture=None):
        """Gear's PAIR map, dense: ``(x_last, x_prev, P_last, P_prev)``, the
        `P` being ``d x / d(x_0, x_{-1})`` as `n x 2n` blocks; with `want_dT`
        two more entries, ``d x_{N-1}/dT`` and ``d x_{N-2}/dT`` -- BOTH rows
        need a period column -- and with `hsens` the two blocks of each
        event column instead.  A view of `_walk_lmm`."""
        w = self._walk_lmm(('pair', x0_in, xm1_in), times, hs, T=T,
                           want_dT=want_dT, hsens=hsens, capture=capture)
        if w.Pk is not None:
            return (w.x_end, w.x_prev, w.P[0], w.P[1],
                    [pk[0] for pk in w.Pk], [pk[1] for pk in w.Pk])
        ## ⚠ ALWAYS THE FULL 2m x 2m MAP, NEVER THE `d x_{N-1}/d x_0`
        ## CORNER.  This used to hand the corner back on the driven path,
        ## and a sub-block of a sensitivity is not a monodromy: it reported
        ## `spectral_radius` 1.279605 for the Q=20 resonator -- ABOVE ONE --
        ## where the analytic per-period decay is exp(-pi/Q) = 0.854636.
        ## For a two-step method the one-period map acts on the PAIR, and
        ## its spectrum carries the discretisation's parasitic roots beside
        ## the physical multipliers; BDF-2's is 1/3 per STEP, (1/3)^N over a
        ## period, and `_spectral_report` separates them by eigenvector
        ## block structure anyway.
        self._monodromy = np.vstack((w.P[0], w.P[1]))
        if want_dT:
            return w.x_end, w.x_prev, w.P[0], w.P[1], w.Pt[0], w.Pt[1]
        return w.x_end, w.x_prev, w.P[0], w.P[1]

    def _traverse_factored(self, x0_in, xm1_in, times, hs, T=None,
                           want_dT=False):
        """Gear's PAIR map, factored: ``(C0, steps, x_last, x_prev[, Pt_last,
        Pt_prev])`` -- the two opening capacitances and per step
        ``(lu, C, alphas, b)``.  A view of `_walk_lmm`."""
        w = self._walk_lmm(('pair', x0_in, xm1_in), times, hs, T=T,
                           dense=False, keep=True, want_dT=want_dT)
        if want_dT:
            return w.opening, w.steps, w.x_end, w.x_prev, w.Pt[0], w.Pt[1]
        return w.opening, w.steps, w.x_end, w.x_prev

    def _traverse_factored_plain(self, x_in, T, times, hs, want_dT=False,
                                 open_at_x0=False):
        """The PLAIN map, factored: ``(opening, steps, x0, x_end, Pt)`` with
        ``opening = (C_open, a_open, b_open, pq_open)``.  A view of
        `_walk_lmm`."""
        w = self._walk_lmm(('plain', x_in, open_at_x0), times, hs, T=T,
                           dense=False, keep=True, want_dT=want_dT)
        return (w.opening, w.steps, w.x0, w.x_end,
                (w.Pt[0] if want_dT else None))

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
        ``(steps, xs, Q0, Qend, x_end)`` where each step record is
        ``(Kfacs, Gs, h, A, U, B, V)``: the stage factors
        ``K_i = LU(C(Y_i) + h lambda G(Y_i))``, the stage conductances, and the
        tableau blocks the sensitivity recursion needs.
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
        for _j, t in enumerate(times[1:]):
            h = hs[min(_j, len(hs) - 1)]
            xn = x
            x = copy(self.solve_timestep(xn, t, h))
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
            steps.append((Kfacs, Gs, float(h), A, U, B, V, Ks))
            xs.append(np.asarray(x, dtype=float))
        return steps, xs, Q0, np.asarray(tr._glm_Q[0], dtype=float), x

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
        for _si, rec in enumerate(steps):
            ## 'closing': dh/dT = 1 on the last step, 0 elsewhere.  The
            ## explicit `h = frac T` in the stage: for an AUTONOMOUS circuit
            ## `K_j = -i(Y_j)` carries no time of its own, so the only new
            ## term is the stage sum (see `_glm_step`)
            fT = (None if T is None else
                  ((1.0 if _si == _nsteps - 1 else 0.0) if closing
                   else rec[2] / T))
            P, D = _glm_step(rec, P, fT)
        return P, (D[-1] if D is not None else None)

    def _walk_glm(self, x_in, T, times, hs, dense=True, keep=False,
                  want_dT=False):
        """ONE WALK OF THE PERIOD UNDER A NORDSIECK GLM: `_glm_period_blocks`,
        then, when `dense`, the sensitivity recursion (`_glm_propagate`) for
        the map and, with `want_dT`, its period column.  The blocks are always
        collected -- the dense map is propagated FROM them -- so `keep` only
        decides whether they are handed back (a `FactoredPeriod` of kind
        'glm', on the Nordsieck state, `width` ``r*m``).

        The map shot on is ``x_0 -> x_N``: the Nordsieck vector is built from
        ``x_0`` by `Transient._glm_startup` at the top of the period and
        propagated by the method to the end.  ⚠ THE JACOBIAN IS APPROXIMATE
        BY CONSTRUCTION and the residual is not: only ``dQ_0/dx_0 = C(x_0)``
        is carried into the recursion, the startup's dependence of the higher
        Nordsieck components on ``x_0`` (p Radau substeps and an interpolant)
        is dropped.  So the converged fixed point is the method's own, exactly;
        what the approximation can cost is Newton iterations.  MEASURED on the
        index-2 C-V loop: 3 iterations to 1e-12, the same count as radau's
        exact monodromy on the same fixture (roadmap).
        """
        steps, _xs, Q0, _Qend, x_end = self._glm_period_blocks(x_in, times, hs)
        m = self.cir.n - 1
        r = len(Q0)
        Mx = Mt = None
        if dense:
            P = [np.asarray(self._C_at(np.asarray(x_in, dtype=float)),
                            dtype=float)] + [np.zeros((m, m)) for _ in range(r - 1)]
            _Pout, Mx = self._glm_propagate(steps, P)
        if want_dT:
            ## the period column.  The startup's own T-dependence enters
            ## twice: through the SCALING `Q_k = h^k q^(k)` (kept --
            ## `dQ_k/dT = (k/T) Q_k`) and through the Radau substeps at `h/p`
            ## (dropped, as `Mx`'s is).
            Pt = [(k / float(T)) * np.asarray(Q0[k], dtype=float)
                  for k in range(r)]
            _Ptout, Mt = self._glm_propagate(
                steps, Pt, T=float(T),
                closing=(self._period_column == 'closing'))
            Mt = np.asarray(Mt).ravel()
        return _PeriodWalk(kind='glm', fp_kind='glm',
                           x0=np.asarray(x_in, dtype=float),
                           x_end=np.asarray(x_end, dtype=float), P=Mx, Pt=Mt,
                           steps=steps if keep else None, width=r * m)

    def _walk_stage(self, x_in, T, times, hs, dense=True, keep=False,
                    want_dT=False, hsens=None, capture=None):
        """ONE WALK OF THE PERIOD UNDER A RUNGE-KUTTA STAGE METHOD (Radau
        IIA, TR-BDF2, ESDIRK), dense or factored (2026-09-23: it was two
        walks, `_traverse_stage` and `_traverse_factored_stage`).
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
        Finite-difference checked before use (the dT column has been got
        wrong in this file twice -- roadmap 0j).

        ⚠ EVENT COLUMNS (2026-09-21): `hsens[j, k] = d h_j / d theta_k` for
        the state-event unknowns, propagated by the same stage algebra as the
        period column with the per-step weight taken from the matrix instead
        of `h/T`; `capture` names the nodes whose state and sensitivities the
        bordered residual reads.  ⚠ A DRIVEN CIRCUIT'S SOURCES MOVE WITH THE
        GRID: an event column shifts the TIMES the stages are evaluated at,
        so `f = -(i + u(t))` changes by `-u_dot . dt_stage`, `dt_stage =
        tau_n + c_i dh` with `tau_n` the shift of the step's start (the `U_i`
        term).  Without it the FD check read the column 77 % off and the
        event row's derivative with the WRONG SIGN on the PWM fixture.

        Returns a `_PeriodWalk`; `_traverse_stage` is its dense view, and
        `_factored_self_starting` keeps its steps."""
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
        if dense:
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
                    self._captured[_j + 1] = (copy(x), P.copy(),
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

    def _traverse_stage(self, x_in, T, times, hs, want_dT=False, hsens=None,
                        capture=None):
        """The stage map, dense: ``(x0, x_end, M, Mt)`` -- or the event
        columns in place of `Mt` when `hsens` is given.  A view of
        `_walk_stage`."""
        w = self._walk_stage(x_in, T, times, hs, want_dT=want_dT,
                             hsens=hsens, capture=capture)
        self._monodromy = w.P
        if w.Pk is not None:
            return w.x0, w.x_end, w.P, w.Pk
        return w.x0, w.x_end, w.P, (w.Pt if want_dT else None)
