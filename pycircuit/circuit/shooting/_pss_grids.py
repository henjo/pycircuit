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

        A6/B7b.  `Transient` breaks its steps at `cir.next_event`; the PSS
        traversal does not, so a pulse edge inside a step is integrated straight
        through.  This returns step FRACTIONS for `solve(grid=...)` with each
        event in the period placed ON a grid point -- by SNAPPING the nearest
        point onto it when one is close, and INSERTING otherwise, so no
        arbitrarily small step is ever created (the B7c lesson).

        Measured on an RC driven by a `VPulse`, against a 4000-point reference,
        a 40-step uniform grid versus the same grid with its 3-4 event times
        landed on::

            edge offset   uniform      + events        gain
            td = 0        6.787e-03    8.032e-04       8.5x
            td = 0.0125T  5.720e-03    2.099e-04        27x
            td = 0.0092T  2.483e-03    7.910e-04       3.1x

        ⚠ TIME-DRIVEN EVENTS ONLY, AND THE LIMIT IS STRUCTURAL.  `next_event(t)`
        is parameterised by time, so a source's edges can be walked out once and
        placed.  A STATE-DEPENDENT reset -- `Idtmod`'s wrap -- cannot: its
        `next_event` is a linear prediction from the last accepted point and
        returns `inf` before a traversal has started, so there is nothing to walk.
        ⚠ THE REASON RECORDED HERE WAS WRONG, AND IS CORRECTED (2026-09-06).
        This used to say the wrap time "has to become an unknown the Newton
        solves for", citing a map "discontinuous by |dphi| ~ 8.2e-3 on a grid
        point".  Measured: that jump is **grid-INDEPENDENT** -- 1.414214e+09
        at seven different grids, with the wrap ON a node and OFF it alike,
        and unchanged when the exact wrap times are added to the grid.  A
        quantity that does not move when the grid moves is not a grid
        artefact.  The real defect was the fold at the period ENDPOINT, in the
        OUTPUT map, and it is fixed in the RESIDUAL by `_fold_periodic`.
        The conclusion for THIS method is unchanged and now for the right
        reason: a state-dependent reset is not a grid feature, so `event_grid`
        does not help it -- but nor does it need the traversal surgery that
        sentence implied.

        ⚠ AND THIS IS NOT SALTATION.  Saltation was measured and falsified twice
        for this codebase (a switched conductance, a discontinuous injection and
        an `Idtmod` wrap all give a monodromy-vs-FD gap falling at 2.00x per
        doubling, i.e. O(h)); each step already uses its own `Jf`/`C`, which
        describe whichever side of the switch that step is on.  The problem was
        only ever that the grid could not BREAK at the event.
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
        ## ONTO AN EVENT (2026-09-20).  A `tr = 0` pulse is clamped to a
        ## `Pulse.MIN_EDGE` = 1e-18 ramp, so each edge is TWO events 1e-18
        ## apart.  The snap used to land the first and then overwrite that
        ## node with the second, collapsing the ramp onto ONE node on the
        ## post-jump side -- and the step arriving there integrated its
        ## whole length with the post-jump source.  Measured on the pulsed
        ## RC (analytic reference), edges landed, N = 100 .. 1600: EVERY
        ## method first order, gear 1.3e-3 / trap 9.5e-4 / radau 3.0e-4 at
        ## 1600 halving per doubling, the error injected AT the edge node
        ## (radau's = its endpoint weight 0.111 x h dU/tau).  With both ramp
        ## ends kept as nodes -- a 1e-18 step -- radau is exact (8e-10), trap
        ## second order, and gear second order too (3.7e-4 -> 1.0e-5 over
        ## 100 -> 800): variable-step BDF2 with h_n/h_{n-1} -> inf over a
        ## consistent tiny step degenerates to the trapezoidal rule, the
        ## "parasitic" factor w/2 multiplying a difference that is itself
        ## O(1/w).  So the tiny step is the ramp and it stays; the B7c
        ## lesson ("no arbitrarily small step") is about grid points near an
        ## event, not about two events.
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
    #: takes by default -- see the table at `lte_grid`
    LTE_GRID_RELTOL_MAX = 1e-5

    ## the fold's own resolution (2026-09-23): no step above this fraction
    ## of the period, and a funnel of steps no larger than their distance
    ## from a source corner within this reach of it
    FOLD_MAX_STEP = 1.0 / 16.0
    FOLD_EDGE_REACH = 1.0 / 16.0

    def _fold_periods(self, t, xs, iref_probe, T_obs, nbins=4000,
                      boundaries=None, rotate=True):
        """`(fracs, seed)` from the accepted steps of the last settled
        periods folded onto a common phase: the per-phase minimum of the
        local step, re-meshed with growth capped at 2, the seed the state at
        the last rising crossing.  None when fewer than 4 periods fold."""
        t = np.asarray(t, dtype=float).ravel()
        xs = np.asarray(xs, dtype=float)
        keep = [i for i in range(xs.shape[0]) if i != iref_probe]
        nper = int(self.LTE_FOLD_PERIODS) + 1
        if boundaries is not None:
            ## ⚠ A DRIVEN CIRCUIT'S PERIODS ARE THE DRIVE'S (2026-09-21): its
            ## sources are functions of absolute time and its edges sit at
            ## fixed phases of the drive, so the boundaries are k T exactly
            ## and phase 0 is the drive's t = 0.  Folded on a CROSSING instead
            ## (the oscillator rule), the switched-capacitor sampler's finest
            ## steps landed at phases 0.26 .. 0.31 of the drive with the
            ## switch edges at 0.75 -- the grid half a period off where it
            ## mattered, and the default path.
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
            ## ⚠ A STEP THAT STRADDLES THE BOUNDARY KEEPS ITS OWN WIDTH
            ## (2026-09-21).  Clipped to its sliver inside the period, the
            ## per-phase MIN over 24 periods whose natural steps drift
            ## against the boundary drove the first and last bins to
            ## 0.0009 T on the switched-capacitor sampler (natural steps
            ## there 0.034 T) -- a fine seam on a flat part of the drive.
            ## The crossing seam of an oscillator sits on an edge that is
            ## already fine, which is why the autonomous fold did not
            ## show it.
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
        ## branch.  At the crossing itself the density jumps ~4x (coarse steps
        ## approach an edge, fine ones follow it: a 0.25 seam ratio).  At the
        ## FINEST phase the seam sits mid-edge, where the PPV's sensitivity
        ## is: measured on the same grid, gear's c +28 % against -3.9 % with
        ## the seam on the slow branch, trbdf2's +8 % against +0.4 %, the
        ## period the same either way (+23 / -28, +10 / +7 ppm).  The coarse
        ## opener is `_period_grid`'s doubling ramp.
        ## (a driven grid keeps phase 0 at the drive's t = 0: `solve` reads
        ## its fractions from there, and the landed events with them)
        ## ⚠ THE FOLD OWNS ITS RESOLUTION AT A SOURCE CORNER (2026-09-23).
        ## It used to inherit it by accident: a transient stepping OVER a
        ## pulse corner left a cluster of rejection-driven tiny steps after
        ## every edge, and the fold kept them.  Once the transient landed
        ## its corners (the Runge-Kutta loop, item 2 after E7) that cluster
        ## was gone -- radau, L-stable, damped the RC's decay after the edge
        ## with a quiet estimate and 0.03 T steps, then doubled through the
        ## flat hold phase up to 0.23 T -- and the PSS grid built from the
        ## fold read its first harmonic 4.5e-2 off (was 5.3e-4) with the
        ## period quadrature's spline weights at -0.24 / +0.36.  So: the
        ## step density is funnelled towards every source corner (a step at
        ## phase distance `s` from a corner is at most `s`, never below one
        ## bin) and capped at FOLD_MAX_STEP of the period -- a period
        ## integral needs points where a transient's tolerance does not.
        ## Measured on that grid: H2 1.15 -> 3e-3 from the cap, H1 4.5e-2 ->
        ## 1.1e-3 from the funnel (a ramp from an eighth of the edge width;
        ## from a sixteenth, 2.2e-4), the pin being 2e-3.
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
        ## the walk overshoots the period by at most one step and the whole
        ## grid is then rescaled to sum 1: every ratio is preserved and no
        ## tiny closing step is made (a closing remainder used to trip the
        ## closing-step bound and fire a needless second pass)
        while ph < 1.0:
            h = float(dens[min(int(ph * nbins), nbins - 1)])
            if last is not None:
                h = min(h, 2.0 * last)
            if ph + h >= 1.0:
                ## ⚠ THE WALK ENDS EXACTLY AT THE PERIOD (2026-09-21, item 3
                ## of the non-uniform-grid list).  It used to overshoot by
                ## up to one step and RESCALE every fraction to sum 1 --
                ## which shifts every phase by up to a coarse step times
                ## its phase, so the grid's fine regions landed off the
                ## orbit's edges by up to 2 % of T: six folds of vdP mu = 10
                ## (205 points, the same edge groups to the point) split
                ## into two clusters, gear +1409 .. +1563 ppm against +273
                ## .. +302, trbdf2 on the same grids +150 .. +173 against
                ## +44 .. +46 -- the "gear swing" of the record was this,
                ## the grid's alignment, and it was 5x for every second-
                ## order method.  The remainder is its own step, merged into
                ## the last one when it would be a sliver and the merged
                ## step keeps the controller's own 2x growth -- spread over
                ## the trailing steps, each grown to at most 2x the one
                ## before it, so no sliver and no ratio beyond the bound
                ## (the tail is the seam, the coarsest phase: a step there
                ## grown by a fraction of itself costs nothing measurable).
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

        ⚠ A STATE THAT CROSSES ITS MIDLINE MORE THAN ONCE PER PERIOD
        (2026-09-21, item 2 of the non-uniform-grid list).  The detector
        took EVERY rising crossing of the fastest-swinging state, and a
        strong third harmonic (`sin + 1.2 sin 3`, three crossings) or two
        pulses of different height per period (two crossings) spread the
        spacings past the 1 % rule: `_observed_period` returned None, the
        fold was refused with the "no consistent recurrence" warning and
        the grid fell back to the single window.  Not silent, but a
        capability lost on exactly the harmonic-rich oscillators whose
        grids need folding.  The hint names the period the caller expects
        and picks the branch; a single-crossing state gives the same chain
        as before, bit for bit."""
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

    @staticmethod
    def _observed_period(t, xs, iref_probe, T_hint, nper=8):
        """The period an adaptive run actually shows, from the rising
        crossings of the fastest-swinging state over the last `nper` hints,
        each crossing refined by linear interpolation; None when fewer than
        two spacings exist or they spread by more than 1 % (a run that has
        not settled, or a hint off by more than the window can hold)."""
        from .pss import PSS    # imported when called: pss.py imports this module
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
        tc = PSS._crossing_chain(tw, y, T_hint)
        if len(tc) < 3:
            return None
        d = np.diff(tc)
        if len(d) < 2 or float(np.std(d)) > 1e-2 * float(np.mean(d)):
            return None
        return float(np.mean(d))

    def lte_grid(self, period, x0=None, refnode=gnd, tstab=None,
                 reltol=None, timestep=None, fold=True):
        """Step FRACTIONS for `solve(grid=...)`, derived from an adaptive run.

        B7a.  A transient adapts because it cannot see the future; PSS
        re-solves the SAME interval over and over, so it can be handed a
        grid that was chosen well ONCE and then frozen.  This is the
        derivation side; `_period_grid` is the consumption side and has
        been shipped since item 5.

        Returns `(fracs, seed)` -- the accepted steps of one settled
        period as fractions summing to 1, and the state at the start of
        that window, ready to pass straight back in:

            fracs, seed = pss.lte_grid(period=T)
            pss.solve(period=T, grid=fracs, x0=seed)

        MEASURED on van der Pol at `mu = 100`
        (`benchmarks/pss_lte_grid.py`, the gate this was promoted from):
        1105 derived steps converge where 1105 UNIFORM steps do not, and
        beat a 20000-point uniform grid -- 18x fewer points and -47.3 ppm
        against -60.6.

        ⚠ ON A DRIVEN CIRCUIT (measured 2026-09-21, the switched-capacitor
        sampler, sine and pulse clocks, each method's fold against uniform
        grids of the same and double count, radau uniform-3200 reference):
        the waveform is 6.5-40x better on the fold at the same count (radau
        95 points 3.7e-5 of the swing vs 2.4e-4; trbdf2 261: 1.5e-4 vs
        2.1e-3; gear 125: 1.2e-3 vs 4.5e-2) and the fold at N beats uniform
        at 2N; a uniform grid puts ONE step across a clock's ramp whatever
        N, the fold's finest steps sit inside it.  Held noise on the fold is
        at kT/C for trbdf2 and radau.  ⚠ NOT FOR GEAR'S NOISE: gear's held
        variance on its own fold reads 0.78 x kT/C against 0.98 uniform at
        the same count.  The fold resolves the STATE -- the tracking phase,
        where the state is flat, gets h ~ tau -- and gear's covariance
        recursion carries its recorded O(h/tau) tracking floor through the
        turn-off.  A noise quantity with its own time scale is resolved only
        where the state's grid happens to be fine; use trbdf2 or radau for
        noise on a folded grid, or a uniform grid for gear's.

        ⚠⚠ IT IS FOR STIFF SMOOTH PROBLEMS AND NOT FOR EVENTS, and that
        boundary is measured rather than cautionary.  On a wrapping
        `Idtmod` the derived grid is WORSE than a uniform grid of the same
        count -- max LTE 2.64e+05 against 1.67e+05 times tolerance at
        ~1429 steps -- because the LTE peak sits at the RESET on every
        grid, and no step size makes a discontinuity's local truncation
        error small.  The event half of B7 is a different item: it needs
        the event time to be an unknown the Newton solves for, because a
        grid frozen from a PAST traversal cannot represent an event whose
        time MOVES as the Newton iterates.

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

        ⚠ TWO THINGS ABOUT THE WINDOW, BOTH FOUND WHEN GEAR COULD NOT RUN
        ITS PPV ON THE GRID GEAR HAD PRODUCED (2026-09-21, relaxation van der
        Pol, mu = 10): the window ends at the last NATURAL step, not at the
        `tend`-landing step `Transient.solve` truncates (a tenth of its
        predecessor, hence a 10.5x growth across the period seam, beyond a
        two-step method's zero-stability bound -- the integrator dropped two
        steps to Euler and the transposed replay refused); and the cut is
        ROTATED to the node whose seam ratio is closest to 1, the seed moving
        with it.  The interior ratios are the adaptive run's own (growth
        clamped at 2).  ⚠ AND THE PERIOD YOU PASS MUST BE CLOSE: the grid is
        fractions of it, so with `period` 16 % off the fine regions sit 16 %
        away from the edges and gear's and trbdf2's per-step Newton fail on
        the coarse steps that land there; seeded within a per cent (a coarse
        uniform solve, or the transient's own recurrence) both converge.
        Measured on that fixture, 195 points: gear's period -519 ppm (second
        order under 2x splitting; a uniform gear grid of the same count
        -7214 ppm), trbdf2's -66 ppm; gear's diffusion constant 52 % high
        with its unit multiplier 0.10 off the circle (see `ppv`'s warning),
        trbdf2's 16 %.  On a relaxation oscillator's adaptive grid trbdf2
        uses gear's own grid better than gear does -- by a steady ~8x, the
        methods' constants; the 5x swing between folds of one orbit that
        the record blamed on gear was the fold's own rescale (2026-09-21,
        `_fold_periods`), and every second-order method paid it.

        ⚠ THE ADAPTIVE RUN NEEDS A TIGHTER TOLERANCE THAN THE PSS'S OWN
        (2026-09-21, "Do 2").  Measured on the same fixture, each method on
        the FOLDED grid its own run made, period error (ppm) and diffusion
        constant error against radau N = 8N::

            reltol   gear (pts, ppm, c)         trbdf2                  radau
            1e-4     95    -900    +179 %      146   -2.0    +1.0 %     70   -3.9    -1.9 %
            1e-5     207   +20     +22 %       290   -7.7    +0.27 %    116  -0.17   -0.04 %
            1e-6     425   -3.6    +7.8 %      602   -1.0    +0.05 %    191  -0.004  -0.002 %
            1e-7     928   -1.6    +2.6 %      1294  -0.25   +0.01 %    331   0.000  -0.0002 %

        The PSS default reltol 1e-4 makes a grid gear cannot use (-900 ppm,
        c 2.8x); 1e-5 is the coarsest that serves every method at the 10-ppm
        level, so the default is `min(reltol, LTE_GRID_RELTOL_MAX)`; gear's c
        within a few per cent wants 1e-7.

        ⚠ TO REPAIR A COARSE SOLVE, CALL THIS FROM ITS CONVERGED STATE WITH A
        SHORT `tstab` (2026-09-21; `PSS.refine_grid` was deleted on this
        measurement).  The fold needs `LTE_FOLD_PERIODS` (24) settled
        periods, not the default 200, and a converged state is settled::

            fracs, seed = pss.lte_grid(T, x0=x_converged,
                                       tstab=(PSS.LTE_FOLD_PERIODS + 1) * T)

        Measured on the same fixture from a converged uniform-200 solve
        (radau +4.8 ppm, gear -6264): radau 114 points -0.3 ppm in 7 s, gear
        204 points +28 ppm in 3 s -- the same grids a 200-period run makes.
        `refine_grid`, one subdivision pass from the same solves, reached
        the same accuracy at 1120 (radau) and 4617 (gear) points with step
        ratios of 4 .. 6, because it subdivided every cell uniformly at the
        finest step the transient wanted anywhere inside it and could
        neither move nor remove a point; from a decimated fold likewise
        (57 -> 1077, 102 -> 4507).  Ten to twenty times the points for the
        same answer, and grids gear's replays refuse: it had no niche left.
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

        ## ⚠ THE PSS'S OWN METHOD SHAPES THE GRID (2026-09-21, "Do 2"): this
        ## used to build `Transient(...)` with no integrator, i.e. the
        ## Transient default `Gear2Integrator()` whatever `method` the PSS
        ## was given -- a gear-shaped grid for a radau PSS.  The adaptive run
        ## steps with the PSS's integrator, so the LTE profile the grid
        ## freezes is the one the PSS will pay.
        tr = Transient(self.cir, toolkit=self.toolkit, reltol=rt,
                       integrator=self._integrator_for(self.par.method))
        with _warnings.catch_warnings():
            _warnings.simplefilter('ignore')
            res = tr.solve(refnode=refnode, tend=tstab + T, timestep=h0,
                           x0=x0)
        t = np.asarray(res.sweep_values, dtype=float).ravel()
        xs = np.asarray(res.x, dtype=float)
        ## ⚠ THE PERIOD IS MEASURED FROM THE RUN, NOT TAKEN FROM THE HINT
        ## (2026-09-21).  The grid is FRACTIONS of the period, and with the
        ## hint 16 % off (my relaxation estimate at mu = 10) the fine regions
        ## sat 16 % away from the edges and every method's per-step Newton
        ## failed; the free-period solve then needed the closing convention
        ## and a second pass just to recover.  The adaptive run already
        ## contains the period: the rising crossings of the fastest-swinging
        ## state over the last periods, refined by linear interpolation, read
        ## 19.1003 against a reference 19.0986 (+9e-5, the transient's own
        ## discretisation) from that 16 %-low hint, consistent across
        ## crossings to 5e-5.  The window is cut at THAT period, it is left in
        ## `self.lte_period` for `solve(period=...)`, and a hint more than
        ## 1 % off is WARNED.  An inconsistent detection (spacings spread
        ## above 1 %, or fewer than two) keeps the hint, with a warning.
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
        ## ⚠ THE GRID COMES FROM EVERY SETTLED PERIOD, NOT THE LAST WINDOW
        ## (2026-09-21, Andreas: "Do 1").  A frozen window inherits whichever
        ## rejection-and-growth pattern that one period got: two windows of
        ## the same run parameters read gear -1461 / -22 ppm and trbdf2 -180 /
        ## +10 on a relaxation van der Pol (mu = 10) -- an 8x lottery for
        ## every second-order method (radau, order 5, does not care).  With
        ## the period detected, each of the last `LTE_FOLD_PERIODS` periods'
        ## accepted steps is folded onto a common phase (phase 0 at the
        ## rising crossing the detector found), the per-phase MINIMUM of the
        ## local step over periods is the density, and a grid is re-meshed
        ## from it with growth capped at the controller's own 2x.  Measured on
        ## the same run: gear +40 ppm, trbdf2 +10, radau -0.01 at 207 points
        ## (the 25 % quantile +267 / +13, the median +957 / +16, the last
        ## window re-meshed +150 / +6).  The seam sits at the COARSEST phase
        ## (the slow branch; a mid-edge seam cost gear's c +28 %) and the
        ## opener is `_period_grid`'s doubling ramp; the seed is the
        ## interpolated state at that phase of the last full period.  A
        ## driven circuit folds on the drive's own boundaries k T (phase 0
        ## at the drive's t = 0), an oscillator on the rising-crossing chain
        ## `_crossing_chain` picks, one per period of the hint.  `fold=False`
        ## keeps the single-window cut; a run without a consistent
        ## recurrence falls back to it, warned.
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
        ## ⚠ THE WINDOW ENDS AT THE LAST NATURAL STEP, NOT AT `tend`
        ## (2026-09-21).  `Transient.solve` lands its final step exactly on
        ## `tend`, truncating it -- on a relaxation van der Pol at mu = 10 to
        ## a tenth of its predecessor -- and a periodic grid then carries a
        ## 10.5x GROWTH out of that step that no controller clamp ever saw:
        ## beyond the two-step zero-stability bound, so the integrator
        ## dropped two steps to Euler and gear's transposed replay refused
        ## the grid gear itself had produced.  Rotating the seam only moved
        ## the step into the interior.  The last accepted point before the
        ## landing step ends the window instead.
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
        ## ⚠ THE CUT IS ROTATED TO THE SANEST SEAM (2026-09-21).  The window
        ## is the last period of the adaptive run, so its seam -- the last
        ## step against the first, which a PERIODIC grid joins -- fell
        ## wherever `t[-1] - T` happened to land: on a relaxation van der Pol
        ## at mu = 10 that was right after a rejected tiny step, a 10.5x
        ## GROWTH across the period boundary.  That is beyond a two-step
        ## method's zero-stability bound (1 + sqrt 2), so the integrator
        ## dropped the first step (and the one after it, rebuilding history)
        ## to Euler, and gear's transposed replay refused the two-alpha
        ## steps -- gear could not run its PPV on the grid gear itself had
        ## produced.  The interior ratios are the run's own (growth clamped
        ## at 2 by the controller; shrinks unconditionally stable), so the
        ## cut moves, when the seam is outside the bound, by the least to a
        ## sane one, and the seed moves to that node.
        ## ⚠ ROTATED ONLY WHEN THE SEAM DEMANDS IT, and then by the smallest
        ## rotation to a sane seam: an unconditional move to the flattest
        ## seam put the seed at a phase from which the free-period Newton,
        ## started 8 % off in period, ran away (mu = 4, the tree's own
        ## `lte_grid` test) -- the seed's phase is part of the basin.
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

        RECORDED SCOPE ITEM 5.  A transient adapts because it cannot see
        the future; PSS re-solves the SAME interval over and over, so it can
        be handed a grid that was chosen well ONCE and then frozen.  The
        grid still never moves inside a solve, so the shooting Newton stays
        exact -- freezing is what makes it a Newton, and this changes only
        WHICH frozen grid.

        ⚠ THE RECORDED BLOCKER WAS STALE.  Item 5 said this was "blocked on
        `Transient` accepting a non-uniform grid; `fixed_timestep` is
        uniform-only".  `Transient.solve`'s loop is uniform-only and always
        was -- but PSS never uses that loop.  It drives `solve_timestep`
        directly, one step at a time, and non-uniform steps worked through
        that path unchanged.  Verified before any of this was written.

        `grid` is a sequence of step FRACTIONS of the period, summing to 1.
        Fractions rather than times because an autonomous period is an
        unknown: every step must scale with `T`, or `dh/dT = h/T` -- the
        identity the period column rests on -- stops holding.

        Measured on van der Pol at mu=100 (`benchmarks/pss_lte_grid.py`):
        1105 steps taken from an adaptive transient converge where 1105
        UNIFORM steps do not, and beat a 20000-point uniform grid on
        accuracy -- 18x fewer points and -47.3 ppm against -60.6.
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
        ## COARSE END.  `_traverse` builds `x(0)` from the unknown with ONE
        ## order-dropped step of `hs[0]`, and a grid taken from an adaptive
        ## transient opens wherever that transient's window happened to
        ## start -- on van der Pol at mu=100, `h[0] = 1.4845` against a
        ## MEDIAN of 4.62e-04, 3200x coarser.  That single Euler step moves
        ## the state 7.4%, and the shooting Newton then has to invert a map
        ## whose first act is that step.  Opening at the grid's own finest
        ## step instead costs ONE extra step in 1105 and is the difference
        ## between not converging and converging.
        ##
        ## ⚠ THE 8x IS MEASURED, NOT DERIVED, and it is a guard rather than
        ## a threshold: it exists so grids that already work are left
        ## exactly as the caller wrote them ('2:1' opens at 2x its finest,
        ## 'smooth' at 5x, and neither needs this).  The falsifier is in
        ## the record: on van der Pol's grid, an opening step of 1e-1 still
        ## fails for gear and 1e-2 converges, against a ratio here of
        ## 13939.  Anything between those bounds separates the two cases.
        ## ⚠ AND IT IS ONLY NEEDED WHEN THERE IS A MANUFACTURED STEP TO
        ## PROTECT.  The subdivision exists because `_traverse` builds `x(0)`
        ## with one order-dropped Euler step of `hs[0]` FROM `x_in` -- an
        ## iterate that may be far from the orbit -- and a coarse `hs[0]`
        ## there defeats the inner Newton.  With `x0_unknown` the first step
        ## starts ON the orbit and the same coarse step is solvable:
        ## measured on van der Pol's own LTE grid, the raw 1105-step grid
        ## converges and reaches -47.3 ppm where the subdivided 1106-step
        ## one reaches -73.8.  So the subdivision COSTS accuracy, and it is
        ## skipped where it buys nothing.
        ## ⚠ A RAMP, NOT ONE TINY STEP (2026-09-21).  The subdivision used to
        ## put a single step of `fr.min()` in front of the coarse first step:
        ## on a relaxation van der Pol's own `lte_grid` (mu = 10) a 2.6e-4
        ## step before a 3.5e-2 one, a 138x GROWTH beyond a two-step method's
        ## zero-stability bound, so the integrator dropped two steps to Euler
        ## and gear's transposed replay refused the grid gear had produced.
        ## Skipping the subdivision for the solved-history kind was tried
        ## first and the coarse first step then defeated the per-step Newton
        ## (the entering history is flat at the seed), so the protection is
        ## needed there too.  Doubling steps `d, 2d, 4d, ...` up to the first
        ## step keep every ratio at 2, inside the bound, for ~log2 of the
        ## span in extra points.
        if fr[0] > 8.0 * fr.min() and not getattr(self, '_open_at_x0', False):
            ## ⚠ TOP-DOWN (2026-09-21): built bottom-up as d, 2d, 4d, ..., rest,
            ## the remainder could be a sliver before the next coarse step --
            ## a growth beyond the bound, an Euler drop, gear +918 ppm on its
            ## own folded grid.  Halving from the coarse step down to its
            ## finest scale, with the smallest piece doubled up, keeps EVERY
            ## ratio at 2 or less, the hand-off to the next step included.
            f0 = float(fr[0])
            kk = int(np.ceil(np.log2(f0 / float(fr.min()))))
            pieces = [f0 / 2.0 ** j for j in range(kk, 0, -1)]      # f0/2^k .. f0/2
            pieces = [f0 / 2.0 ** kk] + pieces                       # the remainder, first
            fr = np.concatenate((np.asarray(pieces, dtype=float), fr[1:]))

        ## ⚠ A CALLER'S GRID CAN SILENTLY DEMOTE GEAR-2 TO FIRST ORDER.
        ## `_period_grid` validated positivity and sum-to-1 and nothing
        ## about the INTERIOR ratios.  A two-step method is zero-stable only
        ## up to `h_n / h_{n-1} = 1 + sqrt(2)`, and past it the integrator's
        ## own guard drops the step to Euler -- correct, and invisible.
        ## Measured on a Q=20 resonator driven at resonance with an
        ## alternating 3:1 grid, where half the steps are demoted:
        ##
        ##       npts   uniform    3:1 grid
        ##        100   19.91489    7.99821
        ##        200   20.00960   11.42923
        ##        400   20.02218   14.54985
        ##        800   20.02443   16.85280
        ##
        ## against an analytic peak of 20 V -- 60% low at 100 points,
        ## crawling up at FIRST order, and `converged = True` every time.
        ## Refining does not fix it, because refining a 3:1 grid keeps it
        ## 3:1.  So the warning names the ratio rather than suggesting a
        ## smaller step, which is the advice that does not work here.
        ##
        ## ⚠ This is what item 5 removed the premise for: the literature
        ## note in the class docstring argues Wambacq's objections to
        ## non-uniform BDF "do not bite inside a run" because the grid is
        ## UNIFORM and frozen.  A caller's grid is frozen but not uniform.
        ## ⚠ ONLY REPEATED UP-STEPS COMPOUND (2026-09-20).  An ISOLATED
        ## up-step -- an event ramp, an inserted event -- is harmless at any
        ## ratio: the factor w/2 the recursion applies acts on the difference
        ## across the SMALL step, and their product is (h/2) x', the
        ## trapezoidal predictor.  Measured: gear across a 1e-18 event ramp
        ## (ratio 2.5e9) is second order, and gear on an event grid with
        ## ratios up to 20 at N = 50 is 4x BETTER than uniform.  The 60 %-low
        ## alternating 3:1 grid has a bad ratio every other step, and it is
        ## that repetition this warns about, so a bad ratio counts only when
        ## another lies within `RATIO_ISOLATION` steps of it.  ⚠ This
        ## traversal never drops a step to Euler -- what the text used to say
        ## is `Transient`'s `check_order_drop`, which does not run here.
        if len(fr) > 1 and self._companion_reach() >= 2:
            from pycircuit.circuit.integrator import ZERO_STABILITY_RATIO
            ratios = fr[1:] / fr[:-1]
            bad = np.flatnonzero(ratios > ZERO_STABILITY_RATIO)
            rep = [i for i in bad
                   if np.any((bad != i) & (np.abs(bad - i) <= self.RATIO_ISOLATION))]
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
        ## ⚠ THE 'CLOSING' CONVENTION (2026-09-21, B7c completed): on a
        ## caller's grid the inner steps keep the ABSOLUTE lengths they had
        ## at the first call of this solve (the seed period) and only the
        ## last step moves with `T`.  Under 'proportional' every step is
        ## rescaled with each trial period, so an adaptive grid's fine
        ## regions slide off the edges they were placed on while the outer
        ## Newton hunts for T -- measured on a relaxation van der Pol
        ## (mu = 10) on its own `lte_grid`: with the period seeded 16 % low,
        ## gear's and trbdf2's per-step Newton FAIL on the coarse steps that
        ## then land on the edges.  The closing step is bounded below so a
        ## trial period that would swallow it falls back to proportional
        ## scaling for that evaluation, with a one-time warning.
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
        """Trapezoid weights `W_n = (h_{n-1} + h_n)/(2T)` for samples at
        `fp.times[:N]`, or **None on a uniform grid**.

        Every Fourier coefficient the noise folds take over the period is an
        INTEGRAL, `(1/T) int c(t) e^{-j k w0 t} dt`.  On a uniform grid the
        index DFT (`fft/N`, the inject's `/ N`) IS the periodic trapezoid rule
        and nothing changes -- callers keep their original expression when
        this returns None, so uniform results are bit-identical.  On a
        non-uniform grid (`solve(grid=...)`, `lte_grid`) `1/N`
        is not a quadrature at all: measured 2026-09-19, a flat 59 % error in
        the first harmonic and pnoise converging 30 % off, rate 1.00.
        Measured with these weights on the same 3:1 grid: 1.8e-04 / 4.4e-05 /
        1.1e-05 at N = 200/400/800 -- second order, the rate of the gear
        samples themselves (a rectangle `h_n/T` is first order and becomes
        the bottleneck for every method).  ⚠ Still four orders behind a
        uniform grid on that circuit, and second order caps radau: a derived
        grid is now CORRECT for noise analysis, not competitive.  ⚠ The SAME
        vector must be used by the adjoint inject and by every harmonic of
        `CY`: the dual-consistency identity closes for ANY weights, so it
        cannot catch a mismatch."""
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        h = np.diff(tms[:N + 1])
        if len(h) < 2 or float(np.max(h)) / float(np.min(h)) - 1.0 <= self.UNIFORM_GRID_TOL:
            return None
        T = float(tms[N] - tms[0])
        ## ⚠ THE SECOND-ORDER CAP IS LIFTED (2026-09-21): a periodic cubic
        ## spline rule on an event-free grid, a PIECEWISE one breaking at
        ## the landed event nodes otherwise -- see `periodic_spline_weights`.
        ## Same vector for every consumer.
        return periodic_spline_weights(tms[:N], T, self._event_nodes(tms[:N], T)) / T

    @staticmethod
    def _replay_grid(T, npts, grid):
        """The grid a one-step factored period replays on: uniform for a
        bare point count, the CALLER'S fractions when `grid` is given.

        ⚠ UNTIL 2026-09-20 THESE REPLAYS WERE ALWAYS UNIFORM.  `factored_period`
        passed `len(times) - 1` and the builders built `linspace(0, T, npts+1)`
        whatever grid the solve had run on, so the PPV, the Floquet modes and
        every noise fold of a radau / trbdf2 / esdirk43 run -- and of trap and
        euler through the TR-BDF2 twin -- on a non-uniform grid were computed
        on a uniform replay from the converged `x0`.  Accurate (radau is
        order 5 on most grids), and the reason radau read as "exact on the
        3:1 grid": its adjoint never saw that grid.  Only the solved-history
        (gear) kind replayed on the caller's grid.  Now `factored_period`
        hands the solved fractions down and the replay is on the grid the
        solve was on; a direct call with a bare `npts` is uniform as before.
        The fractions are scaled to `T` (an autonomous period is solved, and
        the fractions are of the period, not of the seed)."""
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
    #: only when another lies within this many steps -- isolated up-steps
    #: (event ramps) do not compound; see the note there
    RATIO_ISOLATION = 4
