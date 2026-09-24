"""The period grids: the event grid, the LTE grid and its fold, the solve grid,
the replay grid, and the period quadrature.
"""
import numpy as np
import warnings
from pycircuit.circuit.circuit import gnd
from ._numerics import periodic_spline_weights


class _PeriodGrids(object):
    """The period grids: the event grid, the LTE grid and its fold, the
    solve grid, the replay grid, and the period quadrature.  A theme of
    `PSS` (see `pss.py`)."""

    def event_grid(self, period, npts=None, grid=None, min_sep=0.05):
        """A step grid with the circuit's EVENT TIMES landed on exactly.

        `Transient` breaks its steps at `cir.next_event`; the PSS traversal
        does not, so a pulse edge inside a step is integrated straight
        through.  This returns step FRACTIONS for `solve(grid=...)` with each
        event in the period placed ON a grid point -- by SNAPPING the nearest
        point onto it when one is close, and INSERTING otherwise, so no
        arbitrarily small step is ever created next to an event.

        ⚠ TIME-DRIVEN EVENTS ONLY, AND THE LIMIT IS STRUCTURAL.  `next_event(t)`
        is parameterised by time, so a source's edges can be walked out once and
        placed.  A STATE-DEPENDENT reset -- `Idtmod`'s wrap -- cannot: its
        `next_event` is a linear prediction from the last accepted point and
        returns `inf` before a traversal has started, so there is nothing to walk.
        The period-map jump at a wrap lies in the OUTPUT map and is fixed in
        the residual (`_close_periodic`), not by the grid; a declared state
        event (`state_events`) is landed by the solve itself, as a Newton
        unknown.

        ⚠ AND FOR A TIME-DRIVEN EVENT THIS IS NOT SALTATION: each step
        already uses its own `Jf`/`C`, which describe whichever side of the
        switch that step is on; the problem is only that the grid cannot
        BREAK at the event.  (A state event's time moves with the state,
        and that motion is what the state-event columns carry.)

        History: `doc/shooting_history.md`, `_PeriodGrids.event_grid`.
        """
        T = float(period)
        if not T > 0.0:
            raise ValueError('event_grid: period must be positive, got %g' % T)
        if grid is not None:
            fr = np.asarray(grid, dtype=float).ravel()
            tot = float(np.sum(fr))
            if not np.isclose(tot, 1.0, rtol=0, atol=1e-9):
                raise ValueError('event_grid: `grid` fractions must sum to 1, '
                                 'they sum to %.12g' % tot)
            pts = np.concatenate(([0.0], np.cumsum(fr)))
            pts[-1] = 1.0
        else:
            if npts is None:
                raise ValueError('event_grid: give either `npts` or `grid`')
            pts = np.linspace(0.0, 1.0, int(npts) + 1)

        ## walk the events across one period
        ev = []
        t = 0.0
        for _ in range(10 * len(pts) + 100):
            e = float(self.cir.next_event(t))
            if not np.isfinite(e) or e >= T * (1.0 - 1e-15):
                break
            if e > T * 1e-15:
                ev.append(e / T)
            if e <= t:
                break
            t = e
        self.event_times = list(ev)
        if not ev:
            return list(np.diff(pts))

        ## ⚠ A SNAP MAY MOVE A GRID POINT ONTO AN EVENT, NEVER AN EVENT
        ## ONTO AN EVENT.  A `tr = 0` pulse edge is TWO events 1e-18 apart
        ## (`Pulse.MIN_EDGE`); collapsed onto one node, the step arriving
        ## there integrates its whole length with the post-jump source and
        ## every method drops to first order.  Both ramp ends stay nodes --
        ## the 1e-18 step IS the ramp, and BDF2 across it degenerates to the
        ## trapezoidal rule.
        landed = set()
        for f in ev:
            j = int(np.argmin(np.abs(pts - f)))
            if j == 0 or j == len(pts) - 1:
                ## never move an endpoint: the period boundary is not ours
                k = 1 if j == 0 else len(pts) - 2
                h = abs(pts[k] - pts[j])
                if abs(pts[k] - f) < min_sep * h and pts[k] not in landed:
                    pts[k] = f
                    landed.add(f)
                    continue
            else:
                h = min(pts[j] - pts[j - 1], pts[j + 1] - pts[j])
                if abs(pts[j] - f) < min_sep * h and pts[j] not in landed:
                    pts[j] = f          ## SNAP -- no tiny step created
                    landed.add(f)
                    continue
            pts = np.append(pts, f)
            pts = np.sort(pts)
            landed.add(f)
        pts = np.unique(pts)
        return list(np.diff(pts))

    #: how many settled periods `lte_grid` folds onto one phase for its grid
    LTE_FOLD_PERIODS = 24

    #: the coarsest reltol the adaptive run of `lte_grid`
    #: takes by default -- see `lte_grid`
    LTE_GRID_RELTOL_MAX = 1e-5

    ## the fold's own resolution: no step above this fraction
    ## of the period, and a funnel of steps no larger than their distance
    ## from a source corner within this reach of it
    FOLD_MAX_STEP = 1.0 / 16.0
    FOLD_EDGE_REACH = 1.0 / 16.0

    def _fold_periods(self, t, xs, iref_probe, T_obs, nbins=4000,
                      boundaries=None, rotate=True):
        """`(fracs, seed)` from the accepted steps of the last settled
        periods folded onto a common phase: the per-phase minimum of the
        local step, re-meshed with growth capped at 2, the seed the state at
        the last rising crossing.  `(None, None)` when fewer than 4 periods
        fold."""
        ## History: `doc/shooting_history.md`, `_PeriodGrids._fold_periods`.
        t = np.asarray(t, dtype=float).ravel()
        xs = np.asarray(xs, dtype=float)
        keep = [i for i in range(xs.shape[0]) if i != iref_probe]
        nper = int(self.LTE_FOLD_PERIODS) + 1
        if boundaries is not None:
            ## ⚠ A DRIVEN CIRCUIT'S PERIODS ARE THE DRIVE'S: its sources are
            ## functions of absolute time and its edges sit at fixed phases of
            ## the drive, so the boundaries are k T exactly and phase 0 is the
            ## drive's t = 0.  Folded on a CROSSING instead (the oscillator
            ## rule), the finest steps can land half a period from the edges.
            tc = np.asarray(boundaries, dtype=float).ravel()
            if len(tc) < 5:
                return None, None
        else:
            j0 = int(np.searchsorted(t, t[-1] - nper * T_obs))
            if j0 >= len(t) - 3:
                return None, None
            W = xs[keep][:, j0:]
            tw = t[j0:]
            k = int(np.argmax(W.max(axis=1) - W.min(axis=1)))
            y = W[k] - 0.5 * (W[k].max() + W[k].min())
            ## the same chain the period came from (one crossing per period)
            tc = self._crossing_chain(tw, y, T_obs)
            if len(tc) < 5:
                return None, None
        phases = (np.arange(nbins) + 0.5) / nbins
        dens = None
        for a, b in zip(tc[:-1], tc[1:]):
            Tk = float(b - a)
            i0 = int(np.searchsorted(t, a))
            i1 = int(np.searchsorted(t, b))
            ts_k = np.concatenate(([a], t[i0:i1], [b]))
            hs_k = np.diff(ts_k)
            if len(hs_k) < 2 or np.any(hs_k <= 0.0):
                continue
            ## ⚠ A STEP THAT STRADDLES THE BOUNDARY KEEPS ITS OWN WIDTH:
            ## clipped to its sliver inside the period, the per-phase MIN over
            ## periods whose natural steps drift against the boundary drives
            ## the first and last bins far below the natural step -- a fine
            ## seam on a flat part of a drive.
            if i0 > 0:
                hs_k[0] = t[i0] - t[i0 - 1]
            if i1 < len(t):
                hs_k[-1] = t[i1] - t[i1 - 1]
            edges = (ts_k - a) / Tk
            idx = np.clip(np.searchsorted(edges, phases, side='right') - 1,
                          0, len(hs_k) - 1)
            d = hs_k[idx] / Tk
            dens = d if dens is None else np.minimum(dens, d)
        if dens is None:
            return None, None
        ## the grid starts at the phase of COARSEST density -- the slow
        ## branch; at the FINEST phase the seam would sit mid-edge, where the
        ## PPV's sensitivity is.  The coarse opener is `_period_grid`'s
        ## doubling ramp.  (a driven grid keeps phase 0 at the drive's t = 0:
        ## `solve` reads its fractions from there, and the landed events
        ## with them)
        ## ⚠ THE FOLD OWNS ITS RESOLUTION AT A SOURCE CORNER: a transient
        ## that lands its corners may step coarsely through the decay after
        ## an edge and a flat hold, too coarse for a period integral.  So the
        ## density is funnelled towards every source corner (a step at phase
        ## distance `s` from a corner is at most `s`, never below one bin)
        ## and capped at FOLD_MAX_STEP of the period.
        if not rotate:
            corners = []
            tt = 0.0
            for _ in range(64):
                tt = float(self.cir.next_event(tt))
                if not np.isfinite(tt) or tt >= T_obs * (1.0 - 1e-9):
                    break
                corners.append(tt / T_obs)
            if corners:
                phase_axis = (np.arange(nbins) + 0.5) / nbins
                for c in corners:
                    sdist = np.abs(phase_axis - c)
                    sdist = np.minimum(sdist, 1.0 - sdist)
                    funnel = np.where(sdist < self.FOLD_EDGE_REACH,
                                      np.maximum(sdist, 1.0 / nbins), np.inf)
                    dens = np.minimum(dens, funnel)
        dens = np.minimum(dens, self.FOLD_MAX_STEP)
        k0 = int(np.argmax(dens)) if rotate else 0
        dens = np.roll(dens, -k0)
        ph0 = float(phases[k0]) if rotate else 0.0
        fr = []
        ph = 0.0
        last = None
        while ph < 1.0:
            h = float(dens[min(int(ph * nbins), nbins - 1)])
            if last is not None:
                h = min(h, 2.0 * last)
            if ph + h >= 1.0:
                ## ⚠ THE WALK ENDS EXACTLY AT THE PERIOD: rescaling an
                ## overshooting walk to sum 1 would shift every phase and land
                ## the fine regions off the orbit's edges.  The remainder is
                ## its own step or, when it would be a sliver, spread over the
                ## trailing steps, each grown to at most 2x the one before it
                ## (the tail is the seam, the coarsest phase).
                rem = 1.0 - ph
                if last is not None and rem < 0.5 * last:
                    j = 1
                    while rem > 0.0 and j < len(fr):
                        room = 2.0 * fr[-j - 1] - fr[-j]
                        take = min(rem, room) if room > 0.0 else 0.0
                        fr[-j] += take
                        rem -= take
                        j += 1
                    if rem > 0.0:
                        fr.append(rem)
                else:
                    fr.append(rem)
                break
            fr.append(h)
            ph += h
            last = h
        fr = np.asarray(fr, dtype=float)
        fr = fr / float(fr.sum())          # 1 up to rounding
        ## the seed: the state at that phase of the last full period
        t_seed = float(tc[-2] + ph0 * (tc[-1] - tc[-2]))
        i = int(np.searchsorted(t, t_seed)) - 1
        w = (t_seed - t[i]) / (t[i + 1] - t[i])
        x = (1.0 - w) * xs[:, i] + w * xs[:, i + 1]
        seed = np.concatenate((x[:iref_probe], x[iref_probe + 1:]))
        return fr, seed

    @staticmethod
    def _crossing_chain(tw, y, T_hint):
        """The rising crossings of `y` (times `tw`, linear interpolation),
        ONE PER PERIOD OF THE HINT: every rising crossing when they are
        already one per period, else the chain walked back from the last
        crossing choosing, per step of `T_hint`, the crossing nearest the
        expected time within a quarter period.

        ⚠ A STATE CAN CROSS ITS MIDLINE MORE THAN ONCE PER PERIOD (a strong
        third harmonic, two pulses of different height); taking EVERY rising
        crossing would spread the spacings past `_observed_period`'s 1 %
        rule and refuse the fold, on exactly the harmonic-rich oscillators
        whose grids need folding.  The hint names the period the caller
        expects and picks the branch; a single-crossing state gives every
        crossing.

        History: `doc/shooting_history.md`, `_PeriodGrids._crossing_chain`."""
        up = np.flatnonzero((y[:-1] < 0.0) & (y[1:] >= 0.0))
        if len(up) == 0:
            return np.zeros(0)
        tc = tw[up] - y[up] * (tw[up + 1] - tw[up]) / (y[up + 1] - y[up])
        d = np.diff(tc)
        if len(d) < 2 or float(np.std(d)) <= 1e-2 * float(np.mean(d)):
            return tc
        chain = [float(tc[-1])]
        while True:
            target = chain[-1] - float(T_hint)
            j = int(np.argmin(np.abs(tc - target)))
            if abs(tc[j] - target) > 0.25 * float(T_hint) or tc[j] >= chain[-1]:
                break
            chain.append(float(tc[j]))
        return np.asarray(chain[::-1], dtype=float)

    @classmethod
    def _observed_period(cls, t, xs, iref_probe, T_hint, nper=8):
        """The period an adaptive run actually shows, from the rising
        crossings of the fastest-swinging state over the last `nper` hints,
        each crossing refined by linear interpolation; None when fewer than
        two spacings exist or they spread by more than 1 % (a run that has
        not settled, or a hint off by more than the window can hold)."""
        t = np.asarray(t, dtype=float).ravel()
        xs = np.asarray(xs, dtype=float)
        keep = [i for i in range(xs.shape[0]) if i != iref_probe]
        j0 = int(np.searchsorted(t, t[-1] - nper * float(T_hint)))
        if j0 >= len(t) - 3:
            return None
        W = xs[keep][:, j0:]
        tw = t[j0:]
        swing = W.max(axis=1) - W.min(axis=1)
        if not np.any(swing > 0.0):
            return None
        k = int(np.argmax(swing))
        y = W[k] - 0.5 * (W[k].max() + W[k].min())
        tc = cls._crossing_chain(tw, y, T_hint)
        if len(tc) < 3:
            return None
        d = np.diff(tc)
        if len(d) < 2 or float(np.std(d)) > 1e-2 * float(np.mean(d)):
            return None
        return float(np.mean(d))

    def lte_grid(self, period, x0=None, refnode=gnd, tstab=None,
                 reltol=None, timestep=None, fold=True):
        """Step FRACTIONS for `solve(grid=...)`, derived from an adaptive run.

        A transient adapts because it cannot see the future; PSS re-solves
        the SAME interval over and over, so it can be handed a grid that was
        chosen well ONCE and then frozen.  This is the derivation side;
        `_period_grid` is the consumption side.

        Returns `(fracs, seed)` -- step fractions of one period, summing to
        1, folded from the last `LTE_FOLD_PERIODS` settled periods of the
        run (`fold=False`, or a run the fold cannot use: the accepted steps
        of the last period alone), and the state at the grid's start --
        ready to pass straight back in, with the period the run showed:

            fracs, seed = pss.lte_grid(period=T)
            pss.solve(period=pss.lte_period, grid=fracs, x0=seed)

        Gate: `benchmarks/pss_lte_grid.py` (van der Pol at `mu = 100`).

        ⚠ ON A DRIVEN CIRCUIT the fold beats a uniform grid of twice the
        count (a uniform grid puts ONE step across a clock's ramp).  ⚠ NOT
        FOR GEAR'S NOISE: the fold resolves the STATE, so the flat tracking
        phase gets h ~ tau, and gear's covariance recursion carries its
        O(h/tau) tracking floor through the turn-off.  Use trbdf2 or radau
        for noise on a folded grid, or a uniform grid for gear's.

        ⚠⚠ IT IS FOR STIFF SMOOTH PROBLEMS AND NOT FOR EVENTS such as a
        wrapping `Idtmod`: the LTE peak sits at the RESET on every grid, so
        the derived grid is worse than uniform, and a grid frozen from a
        PAST traversal cannot follow an event whose time MOVES as the
        Newton iterates.

        ⚠ FRACTIONS, NOT TIMES, and that is load-bearing rather than a
        convenience.  An autonomous period is an unknown, so every step
        must scale with `T` or `dh/dT = h/T` -- the identity the period
        column rests on -- stops holding.  See `_period_grid`.

        `tstab` is how long to run before the window is taken; it defaults
        to 200 periods, which is a settling heuristic and not a
        convergence criterion.  ⚠ A RUN THAT HAS NOT SETTLED YIELDS A GRID
        FOR THE WRONG TRAJECTORY, silently -- the fractions will still sum
        to 1 and `solve` will still accept them.  Pass a longer `tstab`,
        or seed `x0` on the orbit, when the answer matters.

        ⚠ The grid is cut at the period the run itself shows, left in
        `self.lte_period`: pass THAT to `solve`, the fractions are of it.

        ⚠ THE ADAPTIVE RUN NEEDS A TIGHTER TOLERANCE THAN THE PSS'S OWN: the
        PSS default 1e-4 makes a grid gear cannot use, 1e-5 serves every
        method at the 10-ppm level, hence `min(reltol, LTE_GRID_RELTOL_MAX)`;
        gear's diffusion constant within a few per cent wants 1e-7.

        ⚠ TO REPAIR A COARSE SOLVE, CALL THIS FROM ITS CONVERGED STATE WITH A
        SHORT `tstab`.  The fold needs `LTE_FOLD_PERIODS` (24) settled
        periods, not the default 200, and a converged state is settled::

            fracs, seed = pss.lte_grid(T, x0=x_converged,
                                       tstab=(PSS.LTE_FOLD_PERIODS + 1) * T)

        (Subdividing a coarse grid's cells instead -- the deleted
        `PSS.refine_grid` -- takes ten to twenty times the points.)

        History: `doc/shooting_history.md`, `_PeriodGrids.lte_grid`.
        """
        import warnings as _warnings
        from pycircuit.circuit.transient import Transient
        T = float(period)
        if not T > 0.0:
            raise ValueError('lte_grid: period must be positive, got %g' % T)
        tstab = 200.0 * T if tstab is None else float(tstab)
        if tstab < 0.0:
            raise ValueError('lte_grid: tstab must not be negative, got %g'
                             % tstab)
        rt = (min(float(self.par.reltol), self.LTE_GRID_RELTOL_MAX)
              if reltol is None else float(reltol))
        h0 = (T / 200.0) if timestep is None else float(timestep)

        ## ⚠ THE PSS'S OWN METHOD SHAPES THE GRID: the adaptive run steps
        ## with the PSS's integrator, not the `Transient` default, so the LTE
        ## profile the grid freezes is the one the PSS will pay.
        tr = Transient(self.cir, toolkit=self.toolkit, reltol=rt,
                       integrator=self._integrator_for(self.par.method))
        with _warnings.catch_warnings():
            _warnings.simplefilter('ignore')
            res = tr.solve(refnode=refnode, tend=tstab + T, timestep=h0,
                           x0=x0)
        t = np.asarray(res.sweep_values, dtype=float).ravel()
        xs = np.asarray(res.x, dtype=float)
        ## ⚠ THE PERIOD IS MEASURED FROM THE RUN, NOT TAKEN FROM THE HINT:
        ## the grid is FRACTIONS of it, and a period off puts the fine regions
        ## off the edges.  `_observed_period` reads it from the rising
        ## crossings; the window is cut at THAT period, left in
        ## `self.lte_period`, and a hint more than 1 % off is WARNED.  An
        ## inconsistent detection keeps the hint, with a warning.
        ## (drivenness read on the run's own last period of natural steps:
        ## break_events landed them on every edge, so a narrow pulse that 64
        ## even points would straddle is seen -- `_is_autonomous`'s own rule)
        _driven = not self._is_autonomous(np.asarray(t, float)[np.asarray(t, float) >= t[-1] - T])
        T_obs = (T if _driven else
                 self._observed_period(t, xs, self.cir.get_node_index(refnode), T))
        if T_obs is not None:
            if abs(T_obs / T - 1.0) > 1e-2:
                _warnings.warn(
                    'lte_grid: the adaptive run recurs every %.6g s but the '
                    'period passed was %.6g s (%.1f %% off); the grid is cut '
                    'at the observed period -- pass period=pss.lte_period to '
                    'solve(), the grid\'s fractions are of it.'
                    % (T_obs, T, 100.0 * abs(T_obs / T - 1.0)),
                    RuntimeWarning, stacklevel=2)
            T = float(T_obs)
        else:
            _warnings.warn(
                'lte_grid: no consistent recurrence was found in the last '
                'periods of the adaptive run (unsettled, or the period hint '
                'is far off); the grid is cut at the period passed, %.6g s.'
                % T, RuntimeWarning, stacklevel=2)
        self.lte_period = float(T)
        ## ⚠ THE GRID COMES FROM EVERY SETTLED PERIOD, NOT THE LAST WINDOW:
        ## one window inherits whatever rejection-and-growth pattern its
        ## period got.  The last `LTE_FOLD_PERIODS` periods are folded onto a
        ## common phase (`_fold_periods`) -- a driven circuit on the drive's
        ## boundaries k T, an oscillator on the crossing chain
        ## `_crossing_chain` picks.  `fold=False` keeps the single-window cut;
        ## a run without a consistent recurrence falls back to it, warned.
        if fold and T_obs is not None:
            if _driven:
                ## the drive's own period boundaries, the last ones the run holds
                _kmax = int(np.floor(t[-2] / float(T)))
                _bnd = np.arange(max(_kmax - int(self.LTE_FOLD_PERIODS), 1), _kmax + 1) * float(T)
                _fr, _seed = self._fold_periods(t, xs, self.cir.get_node_index(refnode), float(T),
                                                boundaries=_bnd, rotate=False)
            else:
                _fr, _seed = self._fold_periods(t, xs, self.cir.get_node_index(refnode), float(T))
            if _fr is not None:
                return _fr, _seed
        ## a settled window of exactly one period, taken from the END
        ## ⚠ THE WINDOW ENDS AT THE LAST NATURAL STEP, NOT AT `tend`: the
        ## truncated `tend`-landing step would put a growth no controller
        ## clamp saw across the seam, beyond the two-step bound.
        j1 = len(t) - 2 if len(t) > 3 else len(t) - 1
        j0 = int(np.searchsorted(t, t[j1] - T))
        win_t, win_x = t[j0:j1 + 1], xs[:, j0:j1 + 1]
        if len(win_t) < 3:
            raise RuntimeError(
                'lte_grid: the adaptive run put only %d points in the last '
                'period, which is not a grid. Either the transient took '
                'steps larger than the period (raise tstab or lower '
                'timestep) or the period given is wrong.' % len(win_t))
        hs = np.diff(win_t)
        total = float(hs.sum())
        if total <= 0.0:
            raise RuntimeError('lte_grid: the derived window has zero span.')
        ## ⚠ THE CUT IS ROTATED TO A SANE SEAM, ONLY WHEN THE SEAM DEMANDS
        ## IT.  The seam (last step against first, joined by a PERIODIC grid)
        ## falls wherever `t[-1] - T` lands; beyond a two-step method's
        ## zero-stability bound (1 + sqrt 2) it drops steps to Euler and
        ## gear's transposed replay refuses.  The interior ratios are the
        ## run's own (growth clamped at 2), so the cut moves by the least
        ## rotation to a sane seam, the seed with it -- never further: the
        ## seed's phase is part of the free-period Newton's basin (the tree's
        ## own `lte_grid` test).
        from pycircuit.circuit.integrator import ZERO_STABILITY_RATIO
        k = 0
        if len(hs) > 2:
            seam = float(hs[0] / hs[-1])
            if seam > ZERO_STABILITY_RATIO or seam < 1.0 / ZERO_STABILITY_RATIO:
                ratios = hs[1:] / hs[:-1]
                sane = np.flatnonzero((ratios < ZERO_STABILITY_RATIO)
                                      & (ratios > 1.0 / ZERO_STABILITY_RATIO)) + 1
                if len(sane):
                    n_ = len(hs)
                    dist = np.minimum(sane, n_ - sane)
                    k = int(sane[int(np.argmin(dist))])
        fr = np.roll(hs, -k) / total
        ## the seed is the state at the (rotated) window start, with the
        ## reference row removed -- the shape `solve(x0=...)` takes
        iref = self.cir.get_node_index(refnode)
        seed = np.concatenate((win_x[:iref, k], win_x[iref + 1:, k]))
        return fr, seed

    def _period_grid(self, period, npts, grid):
        """`(times, hs)` for one period -- uniform, or a caller's own grid.

        A caller's grid (`lte_grid`) is chosen once and then frozen: it
        never moves inside a solve, so the shooting Newton stays exact --
        freezing is what makes it a Newton, and this changes only WHICH
        frozen grid.  Non-uniform steps work because the traversal drives
        `solve_timestep` directly, one step at a time, not
        `Transient.solve`'s uniform-only loop.

        `grid` is a sequence of step FRACTIONS of the period, summing to 1.
        Fractions rather than times because an autonomous period is an
        unknown: every step must scale with `T`, or `dh/dT = h/T` -- the
        identity the period column rests on -- stops holding.

        History: `doc/shooting_history.md`, `_PeriodGrids._period_grid`.
        """
        if grid is None:
            times, dt = self.toolkit.linspace(0.0, period, num=npts,
                                              endpoint=True, retstep=True)
            return times, np.full(len(times), float(dt))
        fr = np.asarray(grid, dtype=float)
        if fr.ndim != 1 or len(fr) < 2:
            raise ValueError('grid must be a 1-D sequence of at least two '
                             'step fractions, got shape %r' % (fr.shape,))
        if not np.all(fr > 0.0):
            raise ValueError('every grid step fraction must be positive; '
                             'the smallest given is %g' % float(fr.min()))
        total = float(fr.sum())
        if abs(total - 1.0) > 1e-9:
            raise ValueError(
                'grid step fractions must sum to 1 (they are fractions of '
                'the period, so that every step scales with T when the '
                'period is an unknown); they sum to %.12g' % total)
        ## ⚠ THE OPENING STEP IS MANUFACTURED, SO IT MUST NOT BE THE GRID'S
        ## COARSE END.  The plain walk builds `x(0)` with ONE order-dropped Euler
        ## step of `hs[0]` FROM `x_in`, an iterate that may be far from the
        ## orbit, and an adaptive grid can open thousands of times coarser
        ## than its median step -- which defeats the inner Newton.  So a
        ## coarse first step is split into a doubling ramp (every ratio 2,
        ## inside the two-step bound; one tiny step would not be).  The 8x
        ## is a guard, not a threshold: grids that already work ('2:1',
        ## 'smooth') are left as written.  Skipped with `x0_unknown`
        ## (`_open_at_x0`), where the first step starts ON the orbit and the
        ## ramp would only cost accuracy; the solved-history kind needs it.
        if fr[0] > 8.0 * fr.min() and not getattr(self, '_open_at_x0', False):
            ## ⚠ TOP-DOWN: halving from the coarse step down, the smallest
            ## piece doubled up, keeps EVERY ratio at 2 or less, the hand-off
            ## included (bottom-up, the remainder could be a sliver).
            f0 = float(fr[0])
            kk = int(np.ceil(np.log2(f0 / float(fr.min()))))
            pieces = [f0 / 2.0 ** j for j in range(kk, 0, -1)]      # f0/2^k .. f0/2
            pieces = [f0 / 2.0 ** kk] + pieces                       # the remainder, first
            fr = np.concatenate((np.asarray(pieces, dtype=float), fr[1:]))

        ## ⚠ A CALLER'S GRID CAN SILENTLY DEMOTE A TWO-STEP METHOD TO FIRST
        ## ORDER: zero-stable only up to `h_n / h_{n-1} = 1 + sqrt(2)`, it
        ## converges at first order on a grid whose up-steps past that REPEAT
        ## (an alternating 3:1 grid: 60 % low on a Q=20 resonator, reported
        ## converged).  Refining keeps a 3:1 grid 3:1, so the warning names
        ## the ratio, not a smaller step.  An ISOLATED up-step (an event ramp)
        ## is harmless at any ratio -- the recursion's factor w/2 acts on the
        ## difference across the SMALL step, giving the trapezoidal
        ## predictor -- and so is a PAIR, a landed switching window's exit
        ## (its ramp, then the partial step back to the base grid): one
        ## amplification of the parasitic root, bounded, not compounding.
        ## So a bad ratio counts only when TWO others lie within
        ## `RATIO_ISOLATION` steps.  (Measured on the comparator oscillator
        ## under gear: smoothing those pairs into doubling ramps never
        ## improved the period, at 100, 150, 200 or 800 points.)  This
        ## traversal never drops a step to Euler (`Transient`'s
        ## `check_order_drop` does not run here).
        if len(fr) > 1 and self._companion_reach() >= 2:
            from pycircuit.circuit.integrator import ZERO_STABILITY_RATIO
            ratios = fr[1:] / fr[:-1]
            bad = np.flatnonzero(ratios > ZERO_STABILITY_RATIO)
            rep = [i for i in bad
                   if np.sum((bad != i) & (np.abs(bad - i) <= self.RATIO_ISOLATION)) >= 2]
            if rep:
                worst = float(np.max(ratios[rep]))
                warnings.warn(
                    'PSS: this grid steps up by %.3fx where a two-step '
                    'method is zero-stable only to %.3fx, at %d of %d '
                    'interior ratios, and those up-steps REPEAT within %d '
                    'steps of each other, which is what compounds: the '
                    'answer can be far low while reporting converged -- '
                    'measured 60%% low on a Q=20 resonator with an '
                    'alternating 3:1 grid. Refining will NOT fix it: a '
                    'refined 3:1 grid is still 3:1. Smooth the grid so '
                    'adjacent steps stay within %.3fx, or use a one-step '
                    "method (method='trap')."
                    % (worst, ZERO_STABILITY_RATIO, len(rep), len(ratios),
                       self.RATIO_ISOLATION, ZERO_STABILITY_RATIO),
                    RuntimeWarning, stacklevel=3)

        hs = fr * period
        ## ⚠ THE 'CLOSING' CONVENTION: on a caller's grid the inner steps keep
        ## the ABSOLUTE lengths of the seed period and only the last step
        ## moves with `T`; 'proportional' rescaling would slide an adaptive
        ## grid's fine regions off their edges while the Newton hunts for T.
        ## A trial period that would swallow the closing step falls back to
        ## proportional scaling for that evaluation, warned once.
        if getattr(self, '_period_column', 'proportional') == 'closing':
            inner = getattr(self, '_closing_inner', None)
            if inner is None or len(inner) != len(hs) - 1:
                self._closing_inner = np.array(hs[:-1], dtype=float, copy=True)
            else:
                last = float(period) - float(np.sum(inner))
                if last > 0.05 * float(np.min(inner)):
                    hs = np.concatenate((inner, [last]))
                elif not getattr(self, '_closing_warned', False):
                    self._closing_warned = True
                    warnings.warn(
                        'PSS: the closing step would collapse at this trial '
                        'period (%.6g s against inner steps summing to %.6g '
                        's); this evaluation scales the grid proportionally '
                        'instead. Seed the period closer, or pass '
                        "period_column='proportional'."
                        % (float(period), float(np.sum(inner))),
                        RuntimeWarning, stacklevel=3)
        times = np.concatenate(([0.0], np.cumsum(hs)))
        return times, hs

    #: Relative step spread below which the period grid counts as uniform.
    #: Not zero: `event_grid`'s uniform grid differs from `_period_grid`'s in
    #: the last bit.
    UNIFORM_GRID_TOL = 1e-9

    def _event_nodes(self, tms, T):
        """Indices of the grid nodes the landed events sit on (fractions of
        `T` from `tms[0]`, within 1e-9 of a node), for the piecewise
        quadrature's breaks.  Empty when nothing is landed."""
        ev = getattr(self, 'event_times', None)
        if not ev:
            return []
        fr = (np.asarray(tms, dtype=float) - float(tms[0])) / float(T)
        out = []
        for e in ev:
            j = int(np.argmin(np.abs(fr - float(e))))
            if abs(fr[j] - float(e)) < 1e-9:
                out.append(j)
        return out

    def _period_quadrature(self, fp):
        """Period quadrature weights for samples at
        `fp.times[:N]`, or **None on a uniform grid**.

        Every Fourier coefficient the noise folds take over the period is an
        INTEGRAL, `(1/T) int c(t) e^{-j k w0 t} dt`.  On a uniform grid the
        index DFT (`fft/N`, the inject's `/ N`) IS the periodic trapezoid rule
        and nothing changes -- callers keep their original expression when
        this returns None, so uniform results are bit-identical.  On a
        non-uniform grid (`solve(grid=...)`, `lte_grid`) `1/N` is not a
        quadrature at all (measured: a flat 59 % error in the first
        harmonic).  ⚠ The SAME
        vector must be used by the adjoint inject and by every harmonic of
        `CY`: the dual-consistency identity closes for ANY weights, so it
        cannot catch a mismatch."""
        ## History: `doc/shooting_history.md`, `_PeriodGrids._period_quadrature`.
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        h = np.diff(tms[:N + 1])
        if len(h) < 2 or float(np.max(h)) / float(np.min(h)) - 1.0 <= self.UNIFORM_GRID_TOL:
            return None
        T = float(tms[N] - tms[0])
        ## ⚠ PAST THE TRAPEZOID'S SECOND-ORDER CAP: a periodic cubic
        ## spline rule on an event-free grid, a PIECEWISE one breaking at
        ## the landed event nodes otherwise -- see `periodic_spline_weights`.
        ## Same vector for every consumer.
        return periodic_spline_weights(tms[:N], T, self._event_nodes(tms[:N], T)) / T

    @staticmethod
    def _replay_grid(T, npts, grid):
        """The grid a one-step factored period replays on: uniform for a
        bare point count, the CALLER'S fractions when `grid` is given.

        `factored_period` hands the solved fractions down, so the replay is
        on the grid the solve was on; a direct call with a bare `npts` is
        uniform.
        The fractions are scaled to `T` (an autonomous period is solved, and
        the fractions are of the period, not of the seed)."""
        ## History: `doc/shooting_history.md`, `_PeriodGrids._replay_grid`.
        if grid is None:
            times = np.linspace(0.0, float(T), int(npts) + 1)
            return times, np.diff(times)
        fr = np.asarray(grid, dtype=float).ravel()
        if len(fr) != int(npts):
            raise ValueError('factored period: %d step fractions for %d steps'
                             % (len(fr), int(npts)))
        hs = fr / float(fr.sum()) * float(T)
        times = np.concatenate(([0.0], np.cumsum(hs)))
        times[-1] = float(T)
        return times, hs

    #: a step-ratio above `ZERO_STABILITY_RATIO` is reported by `_period_grid`
    #: only when two others lie within this many steps -- an isolated up-step
    #: (an event ramp) or a pair (a window's exit) does not compound; see the
    #: note there
    RATIO_ISOLATION = 4
