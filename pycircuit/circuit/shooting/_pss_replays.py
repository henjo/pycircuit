"""The factored period and its replays: `M v`, `M^T v`, the forced replays and
the sideband fold.
"""
import numpy as np
from ._factored import FactoredPeriod
from ._numerics import _cx_collect


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

    def _replay_transposed(self, fp, v, collect=False, inject=None):
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
            re = self._replay_transposed(fp, v.real, collect, ii[0])
            im = self._replay_transposed(fp, v.imag, collect, ii[1])
            if collect:
                return (re[0] + 1j * im[0], _cx_collect(re[1], im[1]),
                        _cx_collect(re[2], im[2]))
            return re + 1j * im
        v = v.astype(float)
        steps = fp.step_objects()
        if not steps:
            return (v.copy(), [], []) if collect else v.copy()
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
            return out, ts, states
        return out

    def _forced_replay(self, fp, freq, u_ac, y0=None, collect=False):
        """One period of the LINEARISED circuit, driven at `freq`: the same
        steps as `_replay`, each with its source switched on (`sources`),
        so ``y_end = M y0 + w(freq)`` with `w` the particular response from
        a zero state.  That superposition is not incidental -- it is what
        lets PAC solve an `m x m` system instead of an `(N m) x (N m)` one --
        and it holds because the source enters the SAME steps the monodromy
        replays (and, on the plain map, the same consistent-``iq_0`` seed).
        With `collect`, the circuit state at every node.

        ⚠ THE SOLVE IS REAL, THE REPLAY IS COMPLEX: a complex right-hand
        side costs two back-substitutions against the same factors."""
        jw = 2j * np.pi * float(freq)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()
        tms = np.asarray(fp.times, dtype=float)
        v = (np.zeros(fp.width, dtype=complex) if y0 is None
             else np.asarray(y0, dtype=complex).ravel().copy())
        c = fp.seed(v)
        ys = []
        for j, st in enumerate(fp.step_objects()):
            c = st.solve(c, st.sources(u_ac, jw, tms[j], tms[j + 1]))
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
        return acc

    def _sideband_forced(self, fp, freq, l, d, extra=None):
        """The forced (source-injected) part of sideband row `l`, and the
        final costate `g` for the closure: the transposed replay with the
        OUTPUT functional `d` injected at every node (weighted by
        ``exp(-j(l w0 + w) t_n)/N``, or the period quadrature's weight) and
        the source coupling read at every step.  ⚠ The injection is added
        AFTER the step's costate update, so the output at `t_n` couples to
        the sources of steps `< n` (causality); the state at `t_n` is the
        one step `n` enters from.  `extra` is a raw costate injection per
        node (the bordered adjoint's event-row term), at `d`'s position.
        Returns `(forced, g)`.

        History: `doc/shooting_history.md`, `_sideband_forced`."""
        jw = 2j * np.pi * float(freq)
        T = float(fp.T)
        w0 = 2.0 * np.pi / T
        N = len(fp.steps)
        tms = np.asarray(fp.times, dtype=float)
        d = np.asarray(d, dtype=complex).ravel()
        lam = fp.extract_T(np.zeros(fp.width, dtype=complex))
        forced = np.zeros(self.cir.n - 1, dtype=complex)
        _wq = self._period_quadrature(fp)
        steps = fp.step_objects()
        for j in range(N - 1, -1, -1):
            st = steps[j]
            ts = tms[j]
            lam, r = st.adjoint(lam)
            forced = st.source_adjoint(forced, r, jw, ts, tms[j + 1])
            _e = np.exp(-1j * (float(l) * w0 + 2.0 * np.pi * float(freq)) * ts)
            lam = fp.inject(lam, (_e / N if _wq is None else _e * _wq[j]) * d)
            if extra is not None and j in extra:
                lam = fp.inject(lam, np.asarray(extra[j], dtype=complex))
        return forced, fp.seed_T(lam)

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

        Runs the factored traversal ONCE, at the solution, and caches it.
        Lazy on purpose: a factored walk stores `N` factorisations and
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
        if kind in ('stage', 'glm'):
            ## A self-starting method has its own factored map (no opener, no
            ## pair), replayed on the SOLVED grid's fractions (None on a
            ## uniform grid: the uniform replay) -- see `_replay_grid` and
            ## `_factored_self_starting` (history: `doc/shooting_history.md`,
            ## `factored_period`)
            _hs = np.asarray(hs, dtype=float).ravel()
            _uniform = (len(_hs) < 2 or float(np.max(_hs)) / float(np.min(_hs))
                        - 1.0 <= self.UNIFORM_GRID_TOL)
            _fr = None if _uniform else _hs / float(_hs.sum())
            fp = self._factored_self_starting(kind, x0, T, len(times) - 1,
                                              grid=_fr)
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
                                grid=None):
        """The factored period of a SELF-STARTING method about `x0` -- 'stage'
        (`_walk_stage`; a `FactoredPeriod` of kind 'full' or 'dirk') or 'glm'
        (`_glm_period_blocks`; kind 'glm') -- integrated under `method` (or
        the PSS's own) in a transient of its own, on `npts` steps: uniform,
        or `grid`'s fractions (see `_replay_grid`).  One builder for both.

        History: `doc/shooting_history.md`, `_factored_self_starting`."""
        if method is None:
            method = getattr(self.par, 'method', 'euler')
        x0 = np.asarray(x0, dtype=float)
        if x0.shape[0] == self.cir.n:
            x0 = np.concatenate((x0[:self.irefnode], x0[self.irefnode + 1:]))
        times, hs = self._replay_grid(T, npts, grid)
        tr_saved = getattr(self, '_tran', None)
        self._tran = self._new_transient(self._integrator_for(method))
        try:
            w = self._walk(kind, x0, times, hs, dense=False, keep=True)
        finally:
            self._tran = tr_saved
        return w.factored(self, times=times, T=float(T))
