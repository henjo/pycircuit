"""State events as Newton unknowns: the crossings, the remap of the grid onto
them, and the bordered event stage.
"""
import numpy as np
import warnings
import pycircuit.circuit.analysis as analysis
from .events import EventColumns


class _StateEvents(object):
    """State events as Newton unknowns: the crossings, the remap of the grid
    onto them, and the bordered event stage.  A theme of `PSS` (see
    `pss.py`)."""

    ## ⚠ THE STATE-EVENT STAGE -- see `solve`'s docstring.
    ## (History: `doc/shooting_history.md`, `_state_event_rows`.)
    def _state_event_rows(self):
        """`(W, c)`: the circuit's state-event rows on the REDUCED state and
        their thresholds, or `(None, None)` when the circuit declares none."""
        rows = self.cir.state_events() if hasattr(self.cir, 'state_events') else []
        if not rows:
            return None, None
        W = np.array([np.delete(np.asarray(r, dtype=float).ravel(), self.irefnode)
                      for r, _t in rows])
        c = np.array([float(t) for _r, t in rows])
        return W, c

    @classmethod
    def _land_fractions(cls, base, theta, min_sep=0.25, window_steps=None):
        """`(fractions, theta)`: the base grid with each fraction of `theta`
        landed on a node -- a node within `min_sep` of the local step is
        MOVED onto it (no sliver), otherwise one is inserted; `event_grid`'s
        rule.  Endpoints never move.  A switching window is cut into
        `window_steps` steps (the solve's `event_window_steps`; the class's
        `EVENT_WINDOW_STEPS` when not given)."""
        pts = np.concatenate(([0.0], np.cumsum(np.asarray(base, dtype=float))))
        pts[-1] = 1.0
        landed = []
        for f in np.asarray(theta, dtype=float):
            j = int(np.argmin(np.abs(pts - f)))
            if 0 < j < len(pts) - 1:
                h = min(pts[j] - pts[j - 1], pts[j + 1] - pts[j])
                if abs(pts[j] - f) < min_sep * h and pts[j] not in landed:
                    pts[j] = f
                    landed.append(f)
                    continue
            pts = np.sort(np.append(pts, f))
            landed.append(f)
        pts = np.unique(pts)
        ## ⚠ A WINDOW BETWEEN TWO EVENTS GETS ITS OWN SUB-GRID.  A
        ## threshold switch declares both edges of its transition; the
        ## segment between them holds the whole S-curve of the switch, and
        ## in a few steps the stage solve cannot resolve it, whatever the
        ## step count.  A gap between two landed events holding fewer than
        ## `EVENT_WINDOW_STEPS` steps -- and narrower than that many BASE
        ## steps -- is re-cut into `EVENT_WINDOW_STEPS` equal steps (its
        ## base nodes dropped), which the remap then scales with the window.
        ## (Not the immediate neighbours': an inserted event leaves a sliver
        ## beside the window, and against that the rule does not fire.)
        ## ⚠ Until 2026-09-25 it fired only on a window with NO base node
        ## inside (narrower than one base step): a window spanning 2-7 base
        ## steps kept them, fewer where a snapped edge stretched one, and the
        ## staged PWM loop's error was non-monotone in N for every method --
        ## radau 6.5e-5 / 1.0e-5 / 7.0e-5 / 4.2e-6 / 7.5e-5 / 3.9e-6 of the
        ## swing at N = 780..805 against radau at 3200 (5e-6 to 1.1e-5 with
        ## this rule, 5e-7 to 7.6e-6 at 16 steps).
        ## History: `doc/shooting_history.md`, `_land_fractions`.
        base_pts = np.concatenate(([0.0], np.cumsum(np.asarray(base, dtype=float))))
        base_h = np.diff(base_pts)
        landed_sorted = sorted(landed)
        K = int(cls.EVENT_WINDOW_STEPS if window_steps is None else window_steps)
        for a, b in zip(landed_sorted[:-1], landed_sorted[1:]):
            jb = min(int(np.searchsorted(base_pts, 0.5 * (a + b), side='right')) - 1,
                     len(base_h) - 1)
            inside = pts[(pts > a) & (pts < b)]
            if (b - a) < K * base_h[max(jb, 0)] and len(inside) + 1 < K:
                sub = a + (b - a) * np.arange(1, K) / K
                pts = np.sort(np.concatenate((pts[(pts <= a) | (pts >= b)], sub)))
        pts = np.unique(pts)
        return np.diff(pts), np.asarray(landed_sorted, dtype=float)

    #: steps the segment between a switching window's two edges is cut into
    #: -- the default of the solve's `event_window_steps` Parameter.  16 since
    #: 2026-09-25 (Andreas: "Cut by 16"; 8 before): on the PWM loop the
    #: staged radau error against a fine reference was 5e-6 to 1.1e-5 of the
    #: swing at 8, 5e-7 to 7.6e-6 at 16 (N = 780..805)
    EVENT_WINDOW_STEPS = 16

    @staticmethod
    def _event_remap(base, theta0, theta, T):
        """`(fractions, hsens, nodes)`: the base grid (with `theta0` on
        nodes) mapped piecewise-linearly so the anchors `theta0` land on
        `theta`, every step between two anchors scaling with its segment;
        `hsens[j, k] = d h_j / d theta_k` in seconds (`h_j = fraction_j T`)
        and `nodes[k]` the node each event sits on -- the grid's topology
        is frozen, which is what keeps the Newton exact."""
        base = np.asarray(base, dtype=float)
        p = np.concatenate(([0.0], np.cumsum(base)))
        p[-1] = 1.0
        theta0 = np.asarray(theta0, dtype=float)
        theta = np.asarray(theta, dtype=float)
        nodes = [int(np.argmin(np.abs(p - t))) for t in theta0]
        a0 = np.concatenate(([0.0], theta0, [1.0]))
        a1 = np.concatenate(([0.0], theta, [1.0]))
        seg = np.clip(np.searchsorted(a0, p, side='right') - 1, 0, len(a0) - 2)
        L0 = a0[1:] - a0[:-1]
        L1 = a1[1:] - a1[:-1]
        pn = a1[seg] + (p - a0[seg]) * (L1[seg] / L0[seg])
        pn[0] = 0.0
        pn[-1] = 1.0
        for k, j in enumerate(nodes):
            pn[j] = theta[k]
        fr = np.diff(pn)
        K = len(theta0)
        N = len(fr)
        hsens = np.zeros((N, K))
        blen = np.diff(p)
        for j in range(N):
            sj = int(seg[j])
            g = blen[j] / L0[sj] * float(T)
            if sj >= 1:
                hsens[j, sj - 1] -= g
            if sj + 1 <= K:
                hsens[j, sj] += g
        return fr, hsens, nodes

    def _stage_one_crossings(self, x0_ss, times, period, W, c):
        """The crossings of a stage-1 orbit, read off the states the
        capturing traversal just left in `_captured`: every sign change of
        `W_r . x - c_r` between two nodes, interpolated, kept inside (1e-6,
        1 - 1e-6) of the period.  Returns `(base2, th0, Wk, ck)` -- the steps
        the traversal took with the crossings landed on nodes, the
        crossings as fractions, their rows and thresholds -- or `None` when
        the orbit crosses none.

        History: `doc/shooting_history.md`, `_stage_one_crossings`."""
        N = len(times) - 1
        X = np.array([np.asarray(x0_ss, dtype=float)]
                     + [np.asarray(self._captured[j][0], dtype=float)
                        for j in range(1, N + 1)])
        g = X @ W.T - c[None, :]
        found = []
        for r in range(W.shape[0]):
            gr = g[:, r]
            for j in range(N):
                if gr[j] == 0.0 or gr[j] * gr[j + 1] < 0.0:
                    f = (times[j] if gr[j] == gr[j + 1] else
                         times[j] + (times[j + 1] - times[j]) * gr[j] / (gr[j] - gr[j + 1]))
                    f = float(f) / float(period)
                    if 1e-6 < f < 1.0 - 1e-6:
                        found.append((f, r))
        if not found:
            return None
        found.sort()
        ## the steps the traversal takes: `times[1:]` (the autonomous
        ## rebuild hands N times for N fractions, the driven one N + 1)
        base = np.diff(np.asarray(times, dtype=float))[:N] / float(period)
        base2, th0 = self._land_fractions(
            base, [f for f, _r in found],
            window_steps=getattr(self.par, 'event_window_steps', None))
        rows = [r for _f, r in found]
        return base2, th0, W[rows], c[rows]

    def _event_rows_into(self, F, J, off, nodes, Wk, ck, width, ncols):
        """The event rows of a bordered Newton system, rows `off .. off + K`:
        ``F = W_k . x_{nd_k} - c_k``, ``J[:, :width] = W_k P_{nd_k}`` (the
        state's map to the event node: `m` wide, gear's `2m`) and ``J[:,
        off + l] = W_k Pk_{nd_k, l}`` for the `ncols` columns the traversal
        carried (the events', and the period's on an oscillator)."""
        for k, jn in enumerate(nodes):
            xj, Pj, Pkj = self._captured[jn]
            F[off + k] = float(Wk[k] @ np.asarray(xj)) - ck[k]
            J[off + k, :width] = Wk[k] @ np.asarray(Pj)
            for l in range(ncols):
                J[off + k, off + l] = float(Wk[k] @ np.asarray(Pkj[l]).ravel())

    def _finish_state_events(self, th, base2, th0, T, Wk, ck, attempt, columns):
        """After a stage's Newton: the solved crossings `th` become the grid,
        and -- when `attempt` -- the bordered consumers' columns are built on
        it.  `columns(tms, hs, hsens)` traverses that grid capturing every
        node and returns `(M, P_end, P0, fp, P_nodes)`: the fixed-grid map,
        the event columns at the period, the map to node 0 (`I`, or gear's
        `[I, 0]`), and -- on a MATRIX-FREE stage, whose `M` is None -- the
        factored period, whose reverse replays give the event rows'
        derivatives (`EventColumns.from_capture_rows`), and None; the dense
        maps to every node (`columns(..., dense=True)`'s last item) are then
        built on first read.  Returns the grid `(tms, hs)`.

        ⚠ THE GRID EVERY CONSUMER REPLAYS ON IS THE ONE `_period_grid` MAKES
        OF THESE FRACTIONS -- with its opener ramp, which the window
        sub-grid's tiny steps trigger.  So the ramped grid IS the grid from
        here on: `_grid_fracs` carries it (its first piece is the
        smallest, so `_period_grid` will not ramp it again), and the
        monodromy and the event columns are computed on it -- the identity
        remap on it gives every step's sensitivity to the events (the
        ramp's pieces included) and the event nodes.  Pieces computed on
        the unramped fractions would not match `factored_period`'s node
        indices (shifted by the ramp's pieces against them), and a bordered
        consumer would read its event rows at the wrong nodes.

        ⚠ THE TOTAL MONODROMY THROUGH A MOVING EVENT.
        The period map's derivative is not `M` (the grid frozen): a
        perturbation of `x_0` moves the crossing, `dtheta/dx_0 = -Gt^-1 G`
        from the event rows, and the state at the period moves with it
        through the event columns -- the bordered system's Schur
        complement, which is the saltation matrix derived rather than
        guessed.

        The landed crossings join `event_times` (fractions of the period)
        for every kind, so `covariance`, `sampled_noise` and a reader of
        the solve see the same list.

        History: `doc/shooting_history.md`, `_finish_state_events`."""
        m = self.cir.n - 1
        fr, _hsens, _nodes = self._event_remap(base2, th0, th, T)
        _tms_r, _hs_r = self._period_grid(float(T), len(fr), np.asarray(fr, dtype=float))
        fr_r = np.asarray(_hs_r, dtype=float) / float(T)
        self._grid_fracs = fr_r
        self._state_event_fracs = th
        ## ⚠ THE ONE WRITE OF `_monodromy` AFTER A STAGE, for every kind: the
        ## map at the stage's final state on the landed grid -- the TOTAL map
        ## through the events when their columns are built, the grid-frozen
        ## map when not (the one the unbordered consumers then use), None
        ## when neither can be formed.  History:
        ## `doc/shooting_history.md`, `_finish_state_events`.
        M = None
        try:
            _fr_id, hsens_r, nodes_r = self._event_remap(fr_r, th, th, T)
            hs_r = np.asarray(_hs_r, dtype=float)
            tms_r = np.asarray(_tms_r, dtype=float)
            M, P_end, P0, fp, _Pn = columns(tms_r, hs_r, hsens_r)
            if attempt and M is None:
                ev = EventColumns.from_capture_rows(
                    self._captured, len(hs_r), m, nodes_r, Wk, ck, P_end,
                    self._event_rows_T(fp, nodes_r, Wk),
                    lambda: columns(tms_r, hs_r, hsens_r, dense=True)[4])
                self._event_sensitivity = ev.dth
                self._event_columns = ev
            elif attempt:
                ev = EventColumns.from_capture(self._captured, len(hs_r), m, P0,
                                               nodes_r, Wk, ck, P_end)
                self._event_sensitivity = ev.dth
                self._event_columns = ev
                M = M + P_end @ ev.dth
        except (np.linalg.LinAlgError, ValueError, KeyError) as _exc:
            if attempt:
                warnings.warn('PSS: the state-event stage could not assemble its event '
                              'columns (%s); the bordered consumers run unbordered on '
                              'this solve.' % _exc, RuntimeWarning, stacklevel=3)
        self._monodromy = M
        self.event_times = sorted(set([float(e) for e in self.event_times]
                                      + [float(t) for t in th]))
        return np.asarray(_tms_r, dtype=float), np.asarray(_hs_r, dtype=float)

    @staticmethod
    def _event_rows_T(fp, nodes, Wk):
        """``G[k] = W_k P_{nd_k}`` without the dense map to the node: one
        reverse replay each, the row injected at its node (the costate on
        the node's STATE, `FactoredPeriod.inject`; a GLM's map on the state
        takes it as its last stage's)."""
        sm = fp.state_map() if fp.is_glm else fp
        N = len(sm.steps)
        m = np.asarray(Wk).shape[1]
        G = []
        for k, nd in enumerate(nodes):
            inj = np.zeros((N, m))
            inj[nd] = np.asarray(Wk[k], dtype=float)
            G.append(np.asarray(sm.matvec_transposed(np.zeros(sm.width),
                                                     inject=inj), dtype=float))
        return np.array(G)

    @staticmethod
    def _stack_columns(Pk):
        return np.column_stack([np.asarray(pk, dtype=float).ravel() for pk in Pk])

    def _state_event_stage(self, kind, z_ss, info, ier, mesg, period, times,
                           hs, maxiterations, tol, shoot_reltol, alpha,
                           phase_row=None, phase_k=None, matrix_free=False):
        """The bordered second stage: the crossings of the first stage's
        orbit become Newton unknowns.  One stage for every kind that has
        one: the stage methods and gear's pair, driven or free period.

        Unknowns ``(z, theta[, T])``: `z` the entering state (gear's PAIR
        `(x_0, x_{-1})`), `theta` the `K` crossing fractions, and the period
        when the circuit is autonomous (`phase_row` given).  Equations: the
        orbit closes (folded on the idtmod rows), each crossing sits on its
        threshold (`_event_rows_into`), and the phase row.  The period enters
        as one more column of the step-size sensitivities (``d h_j / d T =
        frac_j``).  Returns the stage-1 tuple unchanged when the circuit
        declares no state event or its orbit crosses none; else ``(z, info,
        ier, mesg, period, times, hs)`` on the landed grid (see
        `_finish_state_events`).

        ⚠ `matrix_free`: the same bordered system as a MAT-VEC, never forming
        the period map.  The walk is factored (`keep`) and carries only the
        event columns (and the period's); a Krylov direction ``[v; s]`` costs
        ONE forward replay, which gives ``M v`` and the state at every event
        node together:

            J [v; s] = [ v - a M v - a P_theta s ;
                         W_k (P_{nd_k} v + Pk_{nd_k} s) ;  v[k_phase] ]

        (`_matrix_free_newton`; an oscillator through `_free_period_solve`).
        The stage then leaves `_monodromy` None, as the unstaged matrix-free
        solve does, and the event rows' derivatives `G` come from reverse
        replays (`_finish_state_events`).  Until 2026-09-25 a matrix-free
        solve warned and skipped the stage.

        History: `doc/shooting_history.md`, `_state_event_stage`."""
        autonomous = phase_row is not None
        W, c = self._state_event_rows()
        if W is None:
            return z_ss, info, ier, mesg, period, times, hs
        N = len(times) - 1
        m = self.cir.n - 1
        wm = len(z_ss)

        def evmap(z, T, tms_, hs_, hsens, capture, dense=True):
            """`(z_end, M, Pk)`: the period map with its event columns --
            without `dense`, the walk itself in `M`'s place (factored, its
            steps kept: the matrix-free stage replays them)."""
            w = self._walk(kind, z, tms_, hs_, T=T, hsens=hsens,
                           capture=capture, dense=dense, keep=not dense,
                           open_at_x0=(kind == 'plain'
                                       and getattr(self, '_open_at_x0', False)))
            if kind == 'plain':
                ## the plain walk carries each column as a two-entry ring,
                ## `[P_n, P_{n-1}]`: the map's column is the current one
                cols = (self._stack_columns([pk[0] for pk in w.Pk])
                        if len(w.Pk) else None)
            elif kind == 'pair':
                cols = (np.vstack((self._stack_columns([pk[0] for pk in w.Pk]),
                                   self._stack_columns([pk[1] for pk in w.Pk])))
                        if len(w.Pk) else None)
            else:
                ## (`_monodromy` is written once, after the stage: see
                ## `_finish_state_events`)
                cols = self._stack_columns(w.Pk) if len(w.Pk) else None
            return (np.asarray(w.end(), dtype=float),
                    (np.asarray(w.monodromy(), dtype=float) if dense else w),
                    cols)

        ## the first stage's orbit, captured at every node, for the crossings
        evmap(np.asarray(z_ss, dtype=float), period, times, hs,
              np.zeros((N, 0)), set(range(1, N + 1)))
        found = self._stage_one_crossings(z_ss[:m], times, period, W, c)
        if found is None:
            return z_ss, info, ier, mesg, period, times, hs
        base2, th0, Wk, ck = found
        K = len(th0)
        ncol = K + (1 if autonomous else 0)

        def func_ev(zz):
            z = np.asarray(zz[:wm], dtype=float)
            th = np.asarray(zz[wm:wm + K], dtype=float)
            T = float(zz[-1]) if autonomous else period
            fr, hsens, nodes = self._event_remap(base2, th0, th, T)
            hs_ = fr * float(T)
            tms_ = np.concatenate(([0.0], np.cumsum(hs_)))
            if autonomous:
                hsens = np.column_stack((hsens, fr))   # d h_j / d T = fraction_j
            z_end, M, Pkm = evmap(z, T, tms_, hs_, hsens, set(nodes))
            F = np.zeros(wm + ncol)
            J = np.zeros((wm + ncol, wm + ncol))
            F[:wm] = self._close_periodic(z, z_end, tms_)
            J[:wm, :wm] = np.eye(wm) - alpha * M
            J[:wm, wm:] = -alpha * Pkm
            self._event_rows_into(F, J, wm, nodes, Wk, ck, wm, ncol)
            if autonomous:
                _k, _r = phase_row(z, Pkm[:, K])
                J[wm + K, _k] = 1.0
                F[wm + K] = _r
            return F, J

        def build_mf(zz):
            """`func_ev`'s residual, and its Jacobian as a mat-vec (see the
            docstring)."""
            z = np.asarray(zz[:wm], dtype=float)
            th = np.asarray(zz[wm:wm + K], dtype=float)
            T = float(zz[-1]) if autonomous else period
            fr, hsens, nodes = self._event_remap(base2, th0, th, T)
            hs_ = fr * float(T)
            tms_ = np.concatenate(([0.0], np.cumsum(hs_)))
            if autonomous:
                hsens = np.column_stack((hsens, fr))   # d h_j / d T = fraction_j
            z_end, w, Pkm = evmap(z, T, tms_, hs_, hsens, set(nodes),
                                  dense=False)
            fp = w.factored(self, times=tms_, T=T)
            sm = fp.state_map() if fp.is_glm else fp
            steps = sm.step_objects()
            F = np.zeros(wm + ncol)
            F[:wm] = self._close_periodic(z, z_end, tms_)
            Gt_ = np.zeros((K, ncol))
            for k, jn in enumerate(nodes):
                xj, _Pj, Pkj = self._captured[jn]
                F[wm + k] = float(Wk[k] @ np.asarray(xj)) - ck[k]
                for l in range(ncol):
                    Gt_[k, l] = float(Wk[k] @ np.asarray(Pkj[l]).ravel())
            _k = None
            if autonomous:
                _k, _r = phase_row(z, Pkm[:, K])
                F[wm + K] = _r
            where = {jn: k for k, jn in enumerate(nodes)}

            def mv(ww):
                v, sv = np.asarray(ww[:wm], dtype=float), np.asarray(ww[wm:], dtype=float)
                ## one replay: the map and the state at every event node
                c = sm.seed(v)
                at = {}
                for j, st in enumerate(steps):
                    c = st.solve(c)
                    if (j + 1) in where:
                        at[j + 1] = np.asarray(sm.node(c), dtype=float)
                out = np.zeros(wm + ncol)
                out[:wm] = v - alpha * np.asarray(sm.extract(c), dtype=float) - alpha * (Pkm @ sv)
                for k, jn in enumerate(nodes):
                    out[wm + k] = float(Wk[k] @ at[jn]) + float(Gt_[k] @ sv)
                if autonomous:
                    out[wm + K] = v[_k]
                return out
            ## the columns scaled to unit norm: the crossing and period
            ## columns are exact here, the state's are `I - M`'s, near one
            _cs = np.ones(wm + ncol)
            for l in range(ncol):
                _n = float(np.sqrt(np.sum((alpha * Pkm[:, l]) ** 2)
                                   + np.sum(Gt_[:, l] ** 2)))
                if _n > 0.0:
                    _cs[wm + l] = 1.0 / _n
            return F, mv, _cs

        tol = np.asarray(tol, dtype=float)
        tail = ([float(period)],) if autonomous else ()
        z0 = np.concatenate((np.asarray(z_ss, dtype=float), th0) + tail)
        abst = np.concatenate(tuple([tol] * (wm // m))
                              + (np.full(K, float(np.max(tol))),)
                              + (([tol[phase_k]],) if autonomous else ()))
        xt = np.concatenate(tuple([tol] * (wm // m)) + (np.full(K, 1e-12),)
                            + (([1e-15 * float(period)],) if autonomous
                               else ()))
        ## ⚠ A STAGE THAT FAILS HANDS BACK STAGE 1, NOT AN EXCEPTION: an
        ## iterate can leave the orbit far enough that an inner step does
        ## not converge (measured: trap's first columns on a driven PWM loop,
        ## before they were right, sent `vin` to 1e9 and the whole solve
        ## raised).  The stage is an improvement on stage 1, never a
        ## condition of having a result.
        _solver = None
        if matrix_free:
            def _solver(z0_, ab_, xt_, rt_, mi_):
                return self._matrix_free_newton(build_mf, z0_, ab_, xt_, rt_,
                                                mi_)
        try:
            if autonomous:
                z_new, info, ier, mesg = self._free_period_solve(
                    func_ev, z0, abst, xt, shoot_reltol, maxiterations,
                    float(period), solver=_solver)
            elif matrix_free:
                z_new, info, ier, mesg = _solver(z0, abst, xt, shoot_reltol,
                                                 maxiterations)
            else:
                z_new, info, ier, mesg = analysis.fsolve(
                    func_ev, z0, maxiter=maxiterations, reltol=shoot_reltol,
                    abstol=abst, xtol=xt, toolkit=self.toolkit,
                    full_output=True, line_search=True, floor_detect=True)
        except (analysis.NoConvergenceError, np.linalg.LinAlgError) as _exc:
            warnings.warn(
                'PSS: the state-event stage failed (%s); the solve returns '
                'the first stage, the crossings inside their steps.'
                % str(_exc)[:160], RuntimeWarning, stacklevel=3)
            return z_ss, info, ier, mesg, period, times, hs
        zn = np.asarray(z_new[:wm], dtype=float)
        Tn = float(z_new[-1]) if autonomous else period

        P0 = np.hstack((np.eye(m), np.zeros((m, wm - m))))

        def columns(tms_r, hs_r, hsens_r, dense=None):
            ## the monodromy at FIXED period: the orbit's own map -- or, on a
            ## matrix-free stage, the factored walk (M None); `dense=True`
            ## there is `EventColumns`' deferred `P_nodes` (the dense map to
            ## every node), the stage's captures left as they were
            dense = (not matrix_free) if dense is None else dense
            saved = getattr(self, '_captured', None)
            _ze, M, Pk = evmap(zn, Tn, tms_r, hs_r, hsens_r,
                               set(range(1, len(hs_r) + 1)), dense=dense)
            if not dense:
                return (None, Pk, P0, M.factored(self, times=tms_r, T=Tn),
                        None)
            Pn = None
            if matrix_free:
                Pn = np.array([P0] + [np.asarray(self._captured[j][1], dtype=float)
                                      for j in range(1, len(hs_r) + 1)])
                self._captured = saved
            return M, Pk, P0, None, Pn
        ## (gear's pair builds its columns whether or not the stage
        ## converged; history: `doc/shooting_history.md`, `_state_event_stage`)
        tms_r, hs_r = self._finish_state_events(
            np.asarray(z_new[wm:wm + K], dtype=float), base2, th0, Tn, Wk, ck,
            ier == 1 or kind == 'pair', columns)
        if matrix_free:
            self._monodromy = None
        return zn, info, ier, mesg, Tn, tms_r, hs_r

    def _event_costate_injection(self, fp, v, n):
        """The event nodes' costate injections for a reverse pass of the
        TOTAL map: `-zeta_k W_k` at node `nd_k`, `zeta = Gt^-T
        P_theta^T v` -- the transpose of the saltation, carried by the pass
        to every earlier node; `None` when the solve is not staged.  What
        `ppv()` and `floquet_modes` sample a left vector along the orbit
        with.

        History: `doc/shooting_history.md`, `_event_costate_injection`."""
        _ev = EventColumns.of(self, n)
        if _ev is None:
            return None
        ## ⚠ AT THE CIRCUIT'S WIDTH, NOT THE MAP'S: the event row reads the
        ## node's circuit state, so its costate enters the circuit block --
        ## gear's pair map is `2m` wide and its reverse step adds the
        ## injection to that `m` block (a `2m` row broke `floquet_modes` on
        ## the first staged gear oscillator, 2026-09-24)
        return _ev.costate_injection(v, len(fp.steps), self.cir.n - 1)
