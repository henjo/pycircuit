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

    ## ⚠ THE STATE-EVENT STAGE (2026-09-21) -- see `solve`'s docstring.
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

    @staticmethod
    def _land_fractions(base, theta, min_sep=0.25):
        """`(fractions, theta)`: the base grid with each fraction of `theta`
        landed on a node -- a node within `min_sep` of the local step is
        MOVED onto it (no sliver), otherwise one is inserted; `event_grid`'s
        rule.  Endpoints never move."""
        from .pss import PSS    # imported when called: pss.py imports this module
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
        ## ⚠ A WINDOW BETWEEN TWO EVENTS GETS ITS OWN SUB-GRID (2026-09-21).
        ## A threshold switch declares both edges of its transition; the
        ## segment between them holds the whole S-curve of the switch, and
        ## as ONE step it is integrated by three collocation points across
        ## the curve -- the stage solve stalled at 6e-4 of the swing
        ## whatever the count.  A gap between two landed events narrower
        ## than its neighbours is split into `EVENT_WINDOW_STEPS` steps,
        ## which the remap then scales with the window.
        ## (against the BASE grid's local step, not the immediate
        ## neighbours: an inserted event leaves a sliver beside the window,
        ## and measured against that the rule never fired on the PWM
        ## loop's on-window)
        base_pts = np.concatenate(([0.0], np.cumsum(np.asarray(base, dtype=float))))
        base_h = np.diff(base_pts)
        landed_sorted = sorted(landed)
        for a, b in zip(landed_sorted[:-1], landed_sorted[1:]):
            ia = int(np.argmin(np.abs(pts - a)))
            ib = int(np.argmin(np.abs(pts - b)))
            if ib == ia + 1:
                gap = pts[ib] - pts[ia]
                jb = min(int(np.searchsorted(base_pts, 0.5 * (a + b), side='right')) - 1,
                         len(base_h) - 1)
                if gap < base_h[max(jb, 0)]:
                    sub = pts[ia] + gap * np.arange(1, PSS.EVENT_WINDOW_STEPS) / PSS.EVENT_WINDOW_STEPS
                    pts = np.sort(np.concatenate((pts, sub)))
        pts = np.unique(pts)
        return np.diff(pts), np.asarray(landed_sorted, dtype=float)

    #: steps the segment between a switching window's two edges is cut into
    EVENT_WINDOW_STEPS = 8

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
        the orbit crosses none.  (Refactor E9 item 2: the three stages each
        carried this.)"""
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
        base2, th0 = self._land_fractions(base, [f for f, _r in found])
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
        node and returns `(M, P_end, P0)`: the fixed-grid map, the event
        columns at the period and the map to node 0 (`I`, or gear's `[I,
        0]`).  Returns the grid `(tms, hs)`.

        ⚠ THE GRID EVERY CONSUMER REPLAYS ON IS THE ONE `_period_grid` MAKES
        OF THESE FRACTIONS -- with its opener ramp, which the window
        sub-grid's tiny steps trigger.  The stored pieces were first
        computed on the unramped fractions, and `factored_period`'s node
        indices were shifted by the ramp's pieces against them: the
        bordered sideband solve read its event rows at the wrong nodes
        (dtheta 12 % off, the response 51 %).  So the ramped grid IS the
        grid from here on: `_grid_fracs` carries it (its first piece is the
        smallest, so `_period_grid` will not ramp it again), and the
        monodromy and the event columns are computed on it -- the identity
        remap on it gives every step's sensitivity to the events (the
        ramp's pieces included) and the event nodes.

        ⚠ THE TOTAL MONODROMY THROUGH A MOVING EVENT (2026-09-22, phase B).
        The period map's derivative is not `M` (the grid frozen): a
        perturbation of `x_0` moves the crossing, `dtheta/dx_0 = -Gt^-1 G`
        from the event rows, and the state at the period moves with it
        through the event columns -- the bordered system's Schur
        complement, which is the saltation matrix derived rather than
        guessed.  Measured on the PWM loop: `M` 109 % off the finite
        difference of the staged period map, this 3.3e-8; the dominant
        multiplier 0.691 where `M` read 0.632.

        The landed crossings join `event_times` (fractions of the period)
        for every kind, so `covariance`, `sampled_noise` and a reader of
        the solve see the same list."""
        m = self.cir.n - 1
        fr, _hsens, _nodes = self._event_remap(base2, th0, th, T)
        _tms_r, _hs_r = self._period_grid(float(T), len(fr), np.asarray(fr, dtype=float))
        fr_r = np.asarray(_hs_r, dtype=float) / float(T)
        self._grid_fracs = fr_r
        self._state_event_fracs = th
        if attempt:
            try:
                _fr_id, hsens_r, nodes_r = self._event_remap(fr_r, th, th, T)
                hs_r = np.asarray(_hs_r, dtype=float)
                tms_r = np.asarray(_tms_r, dtype=float)
                M, P_end, P0 = columns(tms_r, hs_r, hsens_r)
                ev = EventColumns.from_capture(self._captured, len(hs_r), m, P0,
                                               nodes_r, Wk, ck, P_end)
                self._event_sensitivity = ev.dth
                self._monodromy = M + P_end @ ev.dth
                self._event_columns = ev
            except (np.linalg.LinAlgError, ValueError, KeyError) as _exc:
                warnings.warn('PSS: the state-event stage could not assemble its event '
                              'columns (%s); the bordered consumers run unbordered on '
                              'this solve.' % _exc, RuntimeWarning, stacklevel=3)
        self.event_times = sorted(set([float(e) for e in self.event_times]
                                      + [float(t) for t in th]))
        return np.asarray(_tms_r, dtype=float), np.asarray(_hs_r, dtype=float)

    @staticmethod
    def _stack_columns(Pk):
        return np.column_stack([np.asarray(pk, dtype=float).ravel() for pk in Pk])

    def _state_event_stage(self, kind, z_ss, info, ier, mesg, period, times,
                           hs, maxiterations, tol, shoot_reltol, alpha,
                           phase_row=None, phase_k=None):
        """The bordered second stage: the crossings of the first stage's
        orbit become Newton unknowns (events phases A/B).  One stage for every
        kind that has one since 2026-09-23 -- it was three, the driven and
        the autonomous stage methods and gear's pair.

        Unknowns ``(z, theta[, T])``: `z` the entering state (gear's PAIR
        `(x_0, x_{-1})`), `theta` the `K` crossing fractions, and the period
        when the circuit is autonomous (`phase_row` given).  Equations: the
        orbit closes (folded on the idtmod rows), each crossing sits on its
        threshold (`_event_rows_into`), and the phase row.  The period enters
        as one more column of the step-size sensitivities (``d h_j / d T =
        frac_j``).  Returns the stage-1 tuple unchanged when the circuit
        declares no state event or its orbit crosses none; else ``(z, info,
        ier, mesg, period, times, hs)`` on the landed grid (see
        `_finish_state_events`)."""
        autonomous = phase_row is not None
        W, c = self._state_event_rows()
        if W is None:
            return z_ss, info, ier, mesg, period, times, hs
        N = len(times) - 1
        m = self.cir.n - 1
        wm = len(z_ss)

        def evmap(z, T, tms_, hs_, hsens, capture):
            """`(z_end, M, Pk)`: the period map with its event columns."""
            if kind == 'pair':
                (x_last, x_prev, P_last, P_prev, Pk_last,
                 Pk_prev) = self._traverse_solved_history(
                    z[:m], z[m:], tms_, hs_, T=T, hsens=hsens, capture=capture)
                return (np.concatenate((np.asarray(x_last, dtype=float),
                                        np.asarray(x_prev, dtype=float))),
                        np.vstack((np.asarray(P_last, dtype=float),
                                   np.asarray(P_prev, dtype=float))),
                        (np.vstack((self._stack_columns(Pk_last),
                                    self._stack_columns(Pk_prev)))
                         if len(Pk_last) else None))
            _x0, x_end, Mx, Pk = self._traverse_stage(
                z, T, tms_, hs_, hsens=hsens, capture=capture)
            return (np.asarray(x_end, dtype=float), np.asarray(Mx, dtype=float),
                    self._stack_columns(Pk) if len(Pk) else None)

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

        tol = np.asarray(tol, dtype=float)
        tail = ([float(period)],) if autonomous else ()
        z0 = np.concatenate((np.asarray(z_ss, dtype=float), th0) + tail)
        abst = np.concatenate(tuple([tol] * (wm // m))
                              + (np.full(K, float(np.max(tol))),)
                              + (([tol[phase_k]],) if autonomous else ()))
        xt = np.concatenate(tuple([tol] * (wm // m)) + (np.full(K, 1e-12),)
                            + (([1e-15 * float(period)],) if autonomous
                               else ()))
        if autonomous:
            z_new, info, ier, mesg = self._free_period_solve(
                func_ev, z0, abst, xt, shoot_reltol, maxiterations,
                float(period))
        else:
            z_new, info, ier, mesg = analysis.fsolve(
                func_ev, z0, maxiter=maxiterations, reltol=shoot_reltol,
                abstol=abst, xtol=xt, toolkit=self.toolkit, full_output=True,
                line_search=True, floor_detect=True)
        zn = np.asarray(z_new[:wm], dtype=float)
        Tn = float(z_new[-1]) if autonomous else period

        def columns(tms_r, hs_r, hsens_r):
            ## the monodromy at FIXED period: the orbit's own map
            _ze, M, Pk = evmap(zn, Tn, tms_r, hs_r, hsens_r,
                               set(range(1, len(hs_r) + 1)))
            return M, Pk, np.hstack((np.eye(m), np.zeros((m, wm - m))))
        ## (gear's pair builds its columns whether or not the stage converged,
        ## as it did before the stages shared this finish)
        tms_r, hs_r = self._finish_state_events(
            np.asarray(z_new[wm:wm + K], dtype=float), base2, th0, Tn, Wk, ck,
            ier == 1 or kind == 'pair', columns)
        return zn, info, ier, mesg, Tn, tms_r, hs_r

    def _event_costate_injection(self, fp, v, n):
        """The event nodes' costate injections for a reverse pass of the
        TOTAL map (2026-09-22): `-zeta_k W_k` at node `nd_k`, `zeta = Gt^-T
        P_theta^T v` -- the transpose of the saltation, carried by the pass
        to every earlier node; `None` when the solve is not staged.  What
        `ppv()` and `floquet_modes` sample a left vector along the orbit
        with."""
        _ev = EventColumns.of(self, n)
        if _ev is None:
            return None
        return _ev.costate_injection(v, len(fp.steps), n)
