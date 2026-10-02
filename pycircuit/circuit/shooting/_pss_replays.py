"""The factored period and its replays: `M v`, `M^T v`, the forced replays and
the sideband fold.
"""
import numpy as np
from ._factored import FactoredPeriod
from ._numerics import _cx_collect
from ._steps import _at_point, _butcher


class _FactoredReplays(object):
    """The factored period and its replays: `M v`, `M^T v`, the forced
    replays and the sideband fold.  A theme of `PSS` (see `pss.py`)."""

    ## ------------------------------------------------------------------
    ## THE REPLAYS -- one of each for every kind of factored period (plain
    ## LMM, gear's pair, stage, GLM).
    ## A map is `extract . step_N ... step_1 . seed` (see `FactoredPeriod`);
    ## the steps know their own algebra (`_LMMStep`, `_StageStep`,
    ## `_GLMStep`).  The `_monodromy_matvec*` names below are the entry
    ## points the matrix-free Newton builds on from raw step lists.
    ## History: `doc/shooting_history.md`, `_FactoredReplays`.
    ## ------------------------------------------------------------------

    def _replay(self, fp, v):
        """`M v` for any factored period: the seed, the steps, the output.

        ⚠ COMPLEX `v` IS TWO REAL REPLAYS, NOT A COMPLEX FACTORISATION.  The
        steps are real, so `M` is a REAL linear map and `M(a + ib) = Ma + i
        Mb` exactly; PAC needs complex products (`I + alpha(f) H`), and
        factoring in complex arithmetic would double the stored factors for
        a map with no imaginary part.  ⚠ The complex split must stay ahead
        of the float cast below: the cast alone discards an imaginary part
        silently, a wrong answer rather than an error.

        History: `doc/shooting_history.md`, `_replay`."""
        v = np.asarray(v)
        if np.iscomplexobj(v):
            return self._replay(fp, v.real) + 1j * self._replay(fp, v.imag)
        c = fp.seed(v.astype(float))
        for st in fp.step_objects():
            c = st.solve(c)
        return fp.extract(c)

    def _replay_transposed(self, fp, v, collect=False, inject=None,
                           with_seed=False):
        """`M^T v` -- the same stored steps REPLAYED BACKWARDS, each solve
        transposed (``M = M_{N-1} ... M_0``, so ``M^T = M_0^T ...
        M_{N-1}^T``).

        ⚠ THIS IS WHY IT COSTS NOTHING TO HAVE.  Demir & Roychowdhury (TCAD
        22(2) 188-196) call reverse integration "often unavailable even in
        existing time-domain simulators" -- true of a forward-only DENSE
        implementation.  The factored period already stores every step's
        factorisation, and every factorisation solves transposed, so the
        reverse pass needs no new integrator and no second traversal.  It
        is the shared dependency of the PPV, adjoint noise and the sideband
        rows.

        ⚠ `collect` HANDS BACK THE PER-STEP COSTATES (`ts[j]`, the step's
        transposed solve(s)) and the adjoint state after each step
        (`states[j]`; for a PPV seed it IS the PPV there, `Phi(T,s_j)^T v(T)
        = v(s_j)`), both in step order.  ⚠ `inject[j]` lands on the adjoint
        state after `steps[j]` is applied backwards, so the output becomes a
        functional over the whole period rather than a value at its end --
        the difference between the response at `t = 0` and a sideband
        coefficient.  ⚠ The collected lists may be NESTED (a stage method's
        `ts` holds per-stage costates per step), so the complex split
        recombines through `_cx_collect`: a flat `a + 1j*b` would multiply
        a LIST by `1j`.

        History: `doc/shooting_history.md`, `_replay_transposed`."""
        v = np.asarray(v)
        inj = None if inject is None else [np.asarray(z) for z in inject]
        if np.iscomplexobj(v) or (
                inj is not None and any(np.iscomplexobj(z) for z in inj)):
            ii = (None, None) if inj is None else (
                [np.real(z) for z in inj], [np.imag(z) for z in inj])
            re = self._replay_transposed(fp, v.real, collect, ii[0], with_seed)
            im = self._replay_transposed(fp, v.imag, collect, ii[1], with_seed)
            if collect:
                out = (re[0] + 1j * im[0], _cx_collect(re[1], im[1]),
                       _cx_collect(re[2], im[2]))
                if with_seed:
                    out = out + (None if re[3] is None
                                 else re[3] + 1j * im[3],)
                return out
            return re + 1j * im
        v = v.astype(float)
        steps = fp.step_objects()
        if not steps:
            if collect:
                return (v.copy(), [], []) + ((None,) if with_seed else ())
            return v.copy()
        w = fp.extract_T(v)
        ts, states = [], []
        for j in range(len(steps) - 1, -1, -1):
            w, r = steps[j].adjoint(w)
            if inj is not None:
                w = fp.inject(w, inj[j])
            if collect:
                ts.append(r)
                states.append(fp.collected(w))
        out = fp.seed_T(w)
        if collect:
            ts.reverse()
            states.reverse()
            if with_seed:
                wq = fp.seed_source_T(w)
                return out, ts, states, (None if wq is None
                                         else np.array(wq, dtype=float))
            return out, ts, states
        return out

    def _forced_replay(self, fp, freq, u_ac, y0=None, collect=False,
                       u_points=None):
        """One period of the LINEARISED circuit, driven at `freq`: the same
        steps as `_replay`, each with its source switched on (`sources`),
        so ``y_end = M y0 + w(freq)`` with `w` the particular response from
        a zero state.  That superposition is not incidental -- it is what
        lets PAC solve an `m x m` system instead of an `(N m) x (N m)` one --
        and it holds because the source enters the SAME steps the monodromy
        replays (and, on the plain map, the same consistent-``iq_0`` seed).
        With `collect`, the circuit state at every node.

        ⚠ THE SOLVE IS REAL, THE REPLAY IS COMPLEX: a complex right-hand
        side costs two back-substitutions against the same factors.

        `u_points` (per step, the source at each of its injection points,
        `injection_times` order: `(points, m)`) replaces the constant `u_ac`
        with a MODULATED source ``W(x(t)) e^{jw t}`` -- a coloured component
        of a noise element (`PAC._coloured_covariance`)."""
        jw = 2j * np.pi * float(freq)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()
        tms = np.asarray(fp.times, dtype=float)
        v = (np.zeros(fp.width, dtype=complex) if y0 is None
             else np.asarray(y0, dtype=complex).ravel().copy())
        c = fp.seed(v)
        ## the source at the period's start, where the opening reads it
        ## (theta's consistent `iq_0`: `seed_source`); periodic, so a
        ## modulated one's value there is the last step's last point's
        u0 = u_ac if u_points is None else _at_point(u_points[-1], -1)
        c = fp.seed_source(c, u0 * np.exp(jw * tms[0]))
        ys = []
        for j, st in enumerate(fp.step_objects()):
            u_j = u_ac if u_points is None else u_points[j]
            c = st.solve(c, st.sources(u_j, jw, tms[j], tms[j + 1]))
            if collect:
                ys.append(fp.node(c).copy())
        return fp.extract(c), ys

    def _forced_replay_cols(self, fp, freqs, u_ac, Y0=None, collect=False,
                            u_points=None):
        """`_forced_replay` at every frequency of `freqs` at once: one
        state COLUMN per frequency, so each step is one block solve instead
        of one solve per frequency (2026-10-02; `PAC.solve` was 94 % its
        forced replays after the P4 work).  `Y0` an `(width, F)` block of
        starting states.  Returns `(Y_end, ys)`, `Y_end` `(width, F)` and
        `ys` (with `collect`) the node states, each `(m, F)`: column `k` is
        `_forced_replay(fp, freqs[k], u_ac, y0=Y0[:, k], ...)` to rounding
        -- a block solve is not the column solves bit for bit.  The stage
        and multistep maps (`sources_cols`); a GLM's replays stay per
        frequency (`PAC._forced_responses`)."""
        freqs = np.asarray(freqs, dtype=float).ravel()
        F = len(freqs)
        jws = 2j * np.pi * freqs
        u_ac = np.asarray(u_ac, dtype=complex).ravel()
        tms = np.asarray(fp.times, dtype=float)
        V = (np.zeros((fp.width, F), dtype=complex) if Y0 is None
             else np.array(Y0, dtype=complex).reshape(fp.width, F))
        c = fp.seed(V)
        u0 = u_ac if u_points is None else _at_point(u_points[-1], -1)
        c = fp.seed_source(c, np.asarray(u0)[:, None]
                           * np.exp(jws * tms[0])[None, :])
        ys = []
        for j, st in enumerate(fp.step_objects()):
            u_j = u_ac if u_points is None else u_points[j]
            c = st.solve(c, st.sources_cols(u_j, jws, tms[j], tms[j + 1]))
            if collect:
                ys.append(fp.node(c).copy())
        return fp.extract(c), ys

    def _forced_replay_transposed(self, fp, freq, xa):
        """`W^T xa` -- the transpose of the map `u -> w(freq)`, the
        many-to-one half and the reason adjoint noise is affordable: the
        forward replay answers what THIS source does at the output, one run
        per source; this answers what the output owes to EVERY source in one
        reverse pass (Okumura et al. 1993 choose the adjoint for exactly
        this).  It is the transposed replay plus each step's source coupling
        (`source_adjoint`), with no second recursion to keep in step."""
        jw = 2j * np.pi * float(freq)
        tms = np.asarray(fp.times, dtype=float)
        w = fp.extract_T(np.asarray(xa, dtype=complex).ravel().copy())
        acc = np.zeros(self.cir.n - 1, dtype=complex)
        steps = fp.step_objects()
        for j in range(len(steps) - 1, -1, -1):
            w, r = steps[j].adjoint(w)
            acc = steps[j].source_adjoint(acc, r, jw, tms[j], tms[j + 1])
        wq = fp.seed_source_T(w)
        if wq is not None:
            acc = acc - np.exp(jw * tms[0]) * np.asarray(wq)
        return acc

    def _reverse_points(self, fp, f_out, d=None, extra=None, lam0=None):
        """One reverse pass that COLLECTS the source coupling rather than
        applying it at one frequency: `(times, cps, g)` with the forced part
        of a row at ANY input frequency `fin` equal to ``-sum_p cps[p]
        exp(2j pi fin times[p])`` (each step's `source_points`).  `d`, if
        given, is the output functional injected at every node with phase
        ``exp(-2j pi f_out t_n)`` times `1/N` or the quadrature weight
        (`l w0 + w` IS `2 pi f_out` for every row sharing that output
        frequency), AFTER the step's costate update -- so the output at
        `t_n` couples to the sources of steps `< n` (causality), and the
        state at `t_n` is the one step `n` enters from; `extra` a raw
        costate injection per node; `lam0` a starting costate (the replay
        from `z`, `_forced_replay_transposed`'s).  `g` is the final costate.

        History: `doc/pss_log_260902.md`, 2026-09-30 (review batch 4)."""
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        if lam0 is None:
            lam = fp.extract_T(np.zeros(fp.width, dtype=complex))
        else:
            lam = fp.extract_T(np.asarray(lam0, dtype=complex).ravel().copy())
        if d is not None:
            d = np.asarray(d, dtype=complex).ravel()
            _wq = self._period_quadrature(fp)
        times, cps = [], []
        steps = fp.step_objects()
        for j in range(N - 1, -1, -1):
            st = steps[j]
            ts = tms[j]
            lam, r = st.adjoint(lam)
            for t, cp in st.source_points(r, ts, tms[j + 1]):
                times.append(t)
                cps.append(cp)
            if d is not None:
                _e = np.exp(-2j * np.pi * float(f_out) * ts)
                lam = fp.inject(lam, (_e / N if _wq is None else _e * _wq[j]) * d)
            if extra is not None and j in extra:
                lam = fp.inject(lam, np.asarray(extra[j], dtype=complex))
        ## (the source the opening reads at the period's start: `seed_source`)
        wq = fp.seed_source_T(lam)
        if wq is not None:
            times.append(float(tms[0]))
            cps.append(np.asarray(wq))
        m = self.cir.n - 1
        return (np.asarray(times, dtype=float),
                np.asarray(cps, dtype=complex).reshape(len(times), m),
                fp.seed_T(lam))

    def _reverse_points_cols(self, fp, f_outs, d=None, lam0=None):
        """`_reverse_points` for several OUTPUT frequencies at once: one
        reverse pass carrying one costate COLUMN per frequency of `f_outs`,
        so each step is one block solve instead of one solve per frequency
        (the speed plan's P4, 2026-10-02).  `lam0` an `(width, F)` block of
        starting costates; no `extra` (an event map's rows stay per
        frequency).  Returns `(times, cps, g)`, `cps` of shape `(points, m,
        F)` and `g` `(width, F)`: column `k` is `_reverse_points(fp,
        f_outs[k], d, lam0=lam0[:, k])` to rounding -- a block solve is not
        the column solves bit for bit."""
        f_outs = np.asarray(f_outs, dtype=float).ravel()
        F = len(f_outs)
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        if lam0 is None:
            lam = fp.extract_T(np.zeros((fp.width, F), dtype=complex))
        else:
            lam = fp.extract_T(
                np.array(lam0, dtype=complex).reshape(fp.width, F))
        if d is not None:
            d = np.asarray(d, dtype=complex).ravel()
            _wq = self._period_quadrature(fp)
        times, cps = [], []
        steps = fp.step_objects()
        for j in range(N - 1, -1, -1):
            st = steps[j]
            ts = tms[j]
            lam, r = st.adjoint(lam)
            for t, cp in st.source_points(r, ts, tms[j + 1]):
                times.append(t)
                cps.append(cp)
            if d is not None:
                _e = np.exp(-2j * np.pi * f_outs * ts)
                _e = _e / N if _wq is None else _e * _wq[j]
                lam = fp.inject(lam, d[:, None] * _e[None, :])
        wq = fp.seed_source_T(lam)
        if wq is not None:
            times.append(float(tms[0]))
            cps.append(np.asarray(wq))
        m = self.cir.n - 1
        return (np.asarray(times, dtype=float),
                np.asarray(cps, dtype=complex).reshape(len(times), m, F),
                fp.seed_T(lam))

    def _monodromy_matvec(self, C0, steps, v):
        """`M v` for gear's solved-history PAIR map from its raw steps --
        ``(v_0, v_{-1}) -> (P_last v, P_prev v)``."""
        return self._replay(FactoredPeriod('solved_history', C0, steps,
                                           None, None, self), v)

    def _monodromy_matvec_transposed(self, C0, steps, v, collect=False,
                                     inject=None):
        """`M^T v` for gear's PAIR map from its raw steps."""
        return self._replay_transposed(
            FactoredPeriod('solved_history', C0, steps, None, None, self), v,
            collect=collect, inject=inject)

    def _monodromy_matvec_plain(self, opening, steps, v):
        """`M v` for the PLAIN map (one entering state) from its raw
        steps and its opening."""
        return self._replay(FactoredPeriod('plain', opening, steps, None,
                                           None, self), v)

    def _monodromy_matvec_transposed_plain(self, opening, steps, v,
                                           collect=False, inject=None):
        """`M^T v` for the PLAIN map from its raw steps (derived, and gated
        against the dense `M` built from the forward replay -- gate any
        from-scratch adjoint here that way: sign inversion is the known
        failure).

        History: `doc/shooting_history.md`,
        `_monodromy_matvec_transposed_plain`."""
        return self._replay_transposed(
            FactoredPeriod('plain', opening, steps, None, None, self), v,
            collect=collect, inject=inject)

    def _monodromy_matvec_stage(self, steps, v):
        """`M v` for a stage method's steps (`_StageStep`)."""
        return self._replay(FactoredPeriod('full', None, steps, None, None,
                                           self), v)

    def factored_period_glm(self, x0, T, npts, method=None, grid=None):
        """The factored period map of a Nordsieck GLM about a periodic point.

        ⚠ THE MAP IS ON THE NORDSIECK STATE, width ``r*m = (p+1)*m``, not on
        ``x``: a multivalue method's period map carries the scaled derivatives
        between steps, so the object that returns to itself is the whole
        vector and its multipliers are ``r*m`` in number -- ``m`` of them the
        circuit's, the rest the method's own, clustered at zero because ``V``'s
        lower block is nilpotent (that is what `verify()` asserts).  A caller
        reading Floquet data off this map must expect the extra zeros, exactly
        as it expects the DAE's structural zeros.
        """
        return self._factored_self_starting('glm', x0, T, npts, method, grid)

    def factored_period(self):
        """The converged period's steps, kept factored -- see `FactoredPeriod`.

        Runs the factored traversal ONCE, at the solution, and caches it --
        on a multistep map that fits `REPLAY_FACTOR_BUDGET`, `solve`'s own
        converged replay WAS that traversal (`_replay_keeps_factors`), and
        this returns it.

        ⚠ ON A TRAP/EULER OSCILLATOR IT IS THE TWIN'S (`monodromy_twin`):
        its grid, orbit and period, not this run's -- read `fp.times` and
        `fp.T`, never this PSS's, alongside it.
        Lazy past the budget: a factored walk stores `N` factorisations and
        `N` capacitances (`2 N m^2` doubles -- ~800 MB at m=1002 and 50
        points), which is a bad trade to impose on every `solve` for the
        callers who never ask.

        ⚠ IT RE-TRAVERSES RATHER THAN REUSING THE NEWTON'S FACTORS, and the
        difference is not efficiency.  The last `build` call inside the
        Newton is at the last TRIAL iterate; the converged answer is the one
        after it.  Reusing those factors would give an operator for a
        trajectory near the solution instead of at it -- a small error, in
        the third figure, of exactly the kind a converged answer absorbs
        without complaint.
        """
        _tw = self.monodromy_twin()
        if _tw is not self:
            return _tw.factored_period()
        if getattr(self, '_period_state', None) is None:
            raise RuntimeError(
                'PSS: no period to factor -- call solve() first. '
                '(factored_period() replays the CONVERGED trajectory, so '
                'there has to be one.)')
        if not self.converged:
            raise RuntimeError(
                'PSS: the shooting solve did not converge, so there is no '
                'periodic operating point to linearise about. A '
                'small-signal analysis over a non-solution is not a '
                'meaningful answer -- fix the PSS run first.')
        if self._factored_period_cache is not None:
            return self._factored_period_cache

        solved, x0, xm1, times, hs, T, x0_unknown = self._period_state
        kind = self._map_kind()
        if (kind == 'stage'
                and getattr(self, '_stage_replay', None) is not None):
            ## the converged replay WAS this walk: factor its record
            ## (`_replay_orbit`, `_stage_period_from`)
            fp = self._stage_period_from(self._stage_replay)
            self._stage_replay = None
        elif kind in ('stage', 'glm'):
            ## A self-starting method has its own factored map (no opener, no
            ## pair), walked on the SOLVED grid itself -- see
            ## `_factored_self_starting` (history: `doc/shooting_history.md`,
            ## `factored_period`)
            fp = self._factored_self_starting(kind, x0, T, len(times) - 1,
                                              solved=(times, hs))
        else:
            w = self._walk('pair' if solved else 'plain',
                           np.concatenate((x0, xm1)) if solved else x0,
                           times, hs, T=T, dense=False, keep=True,
                           open_at_x0=x0_unknown)
            fp = w.factored(self, times=times, T=T)
        self._factored_period_cache = fp
        return fp

    def factored_period_stage(self, x0, T, npts, method=None, grid=None):
        """The factored period map of ANY Runge-Kutta stage method about a
        periodic point `x0` -- a `FactoredPeriod` of kind 'full' (a fully
        implicit tableau: Radau IIA, and any Gauss/Lobatto/higher-Radau one)
        or 'dirk' (lower triangular: TR-BDF2, ESDIRK), one `_StageStep` per
        step in the structure its tableau calls for.

        Integrates the orbit under the stage method (`method`, or the PSS's
        own) on `npts` steps -- uniform, or `grid`'s fractions (see
        `_replay_grid`) -- and returns the `m x m` monodromy, kept factored.

        ⚠ SELF-STARTING, SO NO TWIN.  No order-dropped opener, so this does
        not consult `monodromy_twin`.

        History: `doc/shooting_history.md`, `factored_period_stage`.
        """
        return self._factored_self_starting('stage', x0, T, npts, method, grid)

    def _factored_self_starting(self, kind, x0, T, npts, method=None,
                                grid=None, solved=None):
        """The factored period of a SELF-STARTING method about `x0` -- 'stage'
        (`_walk_stage`; a `FactoredPeriod` of kind 'full' or 'dirk') or 'glm'
        (`_glm_period_blocks`; kind 'glm') -- integrated under `method` (or
        the PSS's own) in a transient of its own, on `npts` steps: uniform,
        or `grid`'s fractions (see `_replay_grid`), or EXACTLY the solve's own
        `solved = (times, hs)` (`factored_period`).  One builder for both.

        ⚠ THE SOLVED GRID ITSELF, NOT ITS FRACTIONS RE-SUMMED (2026-10-02).
        Rebuilt from fractions, a node can come out an ulp off the solved one
        -- and where a landed edge with no ramp (`tr = 0`) sits on that node,
        the step reads the source on the OTHER side of the jump: the walk
        then does not close (`x_N - x_0` 5.7e-11 and 6.9e-9 on a pulsed RC
        at 100 and 200 points, 1.4e-20 where the node happened to land), and
        the map linearised is not the one solved.
        A stage period is built in two halves, the walk recording its stage
        states (`_stage_walk_recorded`) and the factors from that record
        (`_stage_period_from`) -- the converged replay takes the first half
        and leaves the second to `factored_period()`.

        History: `doc/shooting_history.md`, `_factored_self_starting`."""
        if kind == 'stage':
            return self._stage_period_from(self._stage_walk_recorded(
                x0, T, npts, method, grid, solved=solved))
        if method is None:
            method = getattr(self.par, 'method', 'euler')
        x0 = np.asarray(x0, dtype=float)
        if x0.shape[0] == self.cir.n:
            x0 = np.concatenate((x0[:self.irefnode], x0[self.irefnode + 1:]))
        times, hs = (self._replay_grid(T, npts, grid) if solved is None
                     else solved)
        tr_saved = getattr(self, '_tran', None)
        self._tran = self._new_transient(self._integrator_for(method))
        try:
            w = self._walk(kind, x0, times, hs, dense=False, keep=True)
            walk_tr = self._tran
        finally:
            self._tran = tr_saved
        fp = w.factored(self, times=times, T=float(T))
        if kind == 'glm':
            ## ⚠ ITS NODE STARTUPS RUN ON THE TRANSIENT THAT WALKED IT
            ## (`_glm_node_startups`), not on the PSS's own: `method` may
            ## differ from the PSS's, and a PSS never solved has none.
            ## Measured before (2026-09-28): a glm2 PSS asked for a glm3
            ## period collected node states 2.9e-5 off a glm3 PSS's
            ## (`benchmarks/pss_transient_boundary.py` V7).
            fp._walk_tr = walk_tr
        return fp

    def _stage_walk_recorded(self, x0, T, npts, method=None, grid=None,
                             record=None, solved=None):
        """The first half of a stage method's factored period (see
        `_factored_self_starting`, whose grid arguments these are): the walk,
        in a transient of its own, recording each step's stage states and
        factoring none.  `record` is `_walk_stage`'s.  Returns what
        `_stage_period_from` takes -- a few states per step, not the `(1 + s
        + s^2) m^2` of the factors."""
        if method is None:
            method = getattr(self.par, 'method', 'euler')
        x0 = np.asarray(x0, dtype=float)
        if x0.shape[0] == self.cir.n:
            x0 = np.concatenate((x0[:self.irefnode], x0[self.irefnode + 1:]))
        times, hs = (self._replay_grid(T, npts, grid) if solved is None
                     else solved)
        steps = []
        tr_saved = getattr(self, '_tran', None)
        self._tran = self._new_transient(self._integrator_for(method))
        try:
            w = self._walk('stage', x0, times, hs, dense=False, keep=False,
                           record=record, stage_record=steps)
            walk_tr = self._tran
        finally:
            self._tran = tr_saved
        return (walk_tr, w, steps, times, float(T))

    def _stage_period_from(self, walked):
        """The second half (see `_factored_self_starting`): the recorded
        steps factored -- under the transient that walked them, as during
        the walk -- and returned as the `FactoredPeriod`."""
        walk_tr, w, steps, times, T = walked
        integ = walk_tr.base_integrator
        tab = _butcher(integ)
        coupled = integ.is_fully_implicit()
        tr_saved = getattr(self, '_tran', None)
        self._tran = walk_tr
        try:
            w.steps = [self._stage_step(xn, h, tab, coupled, Y=Y)[0]
                       for xn, h, Y in steps]
        finally:
            self._tran = tr_saved
        return w.factored(self, times=times, T=T)
