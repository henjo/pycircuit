"""How much to trust an answer: the monodromy twin, `grid_error`,
`warping_estimate` and the Floquet spectral report.
"""
import numpy as np
import warnings
from pycircuit.circuit.analysis import NoConvergenceError
from pycircuit.circuit.circuit import gnd


class _AccuracyChecks(object):
    """How much to trust an answer: the monodromy twin, `grid_error`,
    `warping_estimate` and the Floquet spectral report.  A theme of `PSS`
    (see `pss.py`)."""

    def monodromy_twin(self):
        """The `PSS` whose monodromy the oscillator surfaces read.

        `self` for a driven circuit, when `self.monodromy == 'native'`,
        when this run's own method already carries a second-order monodromy
        (gear or a stage method), or for a Nordsieck GLM, whose own map on the
        state serves `ppv` and `floquet_modes` (`_GLMPeriod.state_map`; the
        consumers built only on a map on `x` take its twin through
        `_state_twin`).  Otherwise (an autonomous circuit solved with a
        one-step LMM, trap or euler) a twin `PSS` of the same circuit under
        the method `self.monodromy` names, solved once on the SAME grid from
        this orbit's converged state and cached.

        The twin defaults to Radau IIA (`monodromy='radau'`, the most
        accurate by far, at about trbdf2's cost; Andreas, 2026-09-24:
        "Set radau as default twin"); `'trbdf2'` (the default until then),
        `'esdirk43'` and `'gear'` select those twins, and `'native'` reads
        the run's own.  Trap's and
        euler's own monodromy is unusable on a limit cycle (trap's diverges
        with refinement, with either opener; euler's is first order), so
        the STATE keeps the method you asked for and every
        monodromy-derived quantity -- `Q`, the PPV and everything built on
        it, the Floquet modes, the phase-noise surfaces -- comes from the
        twin on the same orbit, re-converged rather than copied (its period
        differs from this one's by O(h^2)).

        History: `doc/shooting_history.md`, `_AccuracyChecks.monodromy_twin`.
        """
        mono = getattr(self, 'monodromy', 'radau')
        ## ⚠ ANY METHOD WHOSE OWN MAP SERVES may be the twin -- the ones that
        ## carry their own monodromy (`carries_own_monodromy`: gear and the
        ## stage methods; a Nordsieck GLM's map is on its Nordsieck state).
        ## Measured on van der Pol, a trap run at 200 points against radau
        ## at 800: lambda2 / c / oscillator d off by 1e-3 / 1e-4 / 2e-3
        ## under a gear twin, 1e-4 / 6e-6 / 8e-5 under trbdf2, 1e-10 /
        ## 2e-11 / 6e-13 under radau (9.5 s against 8.6 s).
        _ok = mono == 'native'
        if not _ok:
            try:
                _ok = self._integrator_for(mono).carries_own_monodromy()
            except ValueError:
                _ok = False
        if not _ok:
            raise ValueError(
                "PSS.monodromy must be 'native' or a method whose own period "
                "map serves as the twin -- 'radau' (the default), 'trbdf2', "
                "'esdirk43' or 'gear' -- not %r" % (mono,))
        _integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        ## ⚠ A GLM USES ITS OWN MAP WHERE IT IS BUILT (Andreas, 2026-09-24:
        ## "For twin use the own map when possible"): measured on van der
        ## Pol, its `ppv` / `floquet_modes` on the state are second / third
        ## order, and glm3's spectrum is 4.6e-8 of radau's where the trbdf2
        ## twin read 9.7e-7.  `carries_own_monodromy` stays False: a GLM may
        ## not BE the twin -- its map's unit multiplier sits O(h^p) off 1
        ## (the startup seam) and its covariance injection is first order;
        ## the surfaces that feel those take a twin (`_state_twin`).
        if (mono == 'native'
                or _integ.carries_own_monodromy()
                or self._map_kind() == 'glm'
                or not getattr(self, 'autonomous', False)):
            ## ⚠ SELF-SUFFICIENT METHODS TAKE NO TWIN
            ## (`carries_own_monodromy`: Gear-2 and every stage method).  Only
            ## a one-step LMM needs one -- its opener seam makes its monodromy
            ## first-order on a limit cycle; twinning any other method is pure
            ## cost and hides the run's own spectrum.  `monodromy` governs only
            ## the twin trap/euler borrow.
            return self
        if getattr(self, '_period_state', None) is None or not self.converged:
            return self
        twin = self._solve_twin(mono)
        self._monodromy_twin = twin
        return twin

    ## The iteration budget of a monodromy twin's re-solve.  A twin is a
    ## POLISH seeded at the converged orbit: a good seed converges in a
    ## handful of iterations, a poor one refuses by NOT converging, so the
    ## caller's `maxiterations` (times the retry ladder) would only make the
    ## refusal slow.  Hitting the cap WARNS before the twin refuses.
    ## History: `doc/shooting_history.md`, `_AccuracyChecks.TWIN_MAXITER`.
    TWIN_MAXITER = 40

    def _solve_twin(self, method):
        """A converged twin of this analysis (the same class) under
        `method`, re-solved once on the SAME grid from this orbit's
        converged state, cached per method.  The machinery behind
        `monodromy_twin`, which picks the method by `self.monodromy`; the
        noise surfaces' hosts (`_lyapunov_host`, `_adjoint_host`) are that
        twin.  Raises if the re-solve does not converge.

        History: `doc/shooting_history.md`, `_AccuracyChecks._solve_twin`.
        """
        cache = self._twins
        if method in cache:
            return cache[method]
        solved, x0, xm1, times, hs, T, x0_unknown = self._period_state
        kw = dict(self._solve_kwargs)
        twin = type(self)(self.cir, toolkit=self.toolkit, irefnode=None,
                          method=method, reltol=self.par.reltol,
                          iabstol=self.par.iabstol, vabstol=self.par.vabstol,
                          event_window_steps=self.par.event_window_steps)
        hs = np.asarray(hs, dtype=float)
        ## the same grid: its fractions when it is not uniform, else the
        ## uniform step (a one-step plain state can carry an `hs` whose
        ## sum is a trial period, so the fractions are the safe object)
        nonuniform = float(hs.max() / hs.min()) > 1.0 + 1e-9
        grid = (hs / float(hs.sum())) if nonuniform else None
        x0r = np.asarray(x0, dtype=float)[:self.cir.n - 1]
        ## ⚠ THE STEP COUNT MUST SURVIVE `solve`'s `int(period / timestep)`:
        ## `T / (T / N)` is not N in floating point; half a step of slack
        ## makes the floor land on N for every T.
        ## ⚠ THE TWIN TAKES THE DEFAULT PHASE RULE, NOT THIS RUN'S: it must
        ## land on the SAME point of the SAME orbit (what the agreement check
        ## below tests), and `phase_rule='reselect'` lands on another phase.
        _asked = max(int(kw.get('maxiterations', 20)), 20)
        _budget = min(_asked, int(self.TWIN_MAXITER))
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                twin.solve(refnode=kw.get('refnode', gnd), period=float(T),
                           x0=x0r, timestep=float(T) / (len(hs) + 0.5),
                           grid=grid, maxiterations=_budget,
                           matrix_free=bool(kw.get('matrix_free', False)))
            _ok = bool(twin.converged)
        except NoConvergenceError:
            _ok = False
        if not _ok and _budget < _asked:
            warnings.warn(
                'PSS.monodromy_twin: the %s re-solve did not converge within '
                'its capped iteration budget (PSS.TWIN_MAXITER = %d; the solve '
                'itself allowed %d).  A twin is seeded at the converged orbit '
                'and a good seed converges in a handful of iterations, so this '
                'usually means the %s orbit is too poor to seed it; if the seed '
                'is good but slow, raise PSS.TWIN_MAXITER.'
                % (method, _budget, _asked, getattr(self.par, 'method', '?')),
                RuntimeWarning, stacklevel=3)
        if not _ok:
            raise RuntimeError(
                'PSS.monodromy_twin: the %s re-solve from the converged '
                '%s orbit did not converge, so no second-order monodromy is '
                'available; try pss.monodromy = "gear" (the other twin) or '
                '"native" to read the one-step method\'s own (first-order, '
                'and on an oscillator its second multiplier is not the '
                'physical one).'
                % (method, getattr(self.par, 'method', '?')))

        ## ⚠⚠ THE TWIN MUST HAVE CONVERGED TO THE ORBIT IT WAS SEEDED ON: from
        ## a poor seed a robust twin (TR-BDF2) can converge to a SPURIOUS
        ## limit cycle and report a plausible wrong `Q` with no error.  Two
        ## convergent methods on one limit cycle agree on period and entering
        ## state to O(h^p), far inside the 0.25 gate (set from that good-case
        ## agreement, not the fixture-dependent failure size, so it
        ## transfers); a too-poor seed is REFUSED, the safe direction.
        Th = float(T)
        dT = abs(float(twin.period) - Th) / max(abs(Th), 1e-30)
        x0t = np.asarray(twin._period_state[1],
                         dtype=float)[:self.cir.n - 1]
        dx = (float(np.linalg.norm(x0t - x0r))
              / (float(np.linalg.norm(x0r)) + 1e-30))
        if dT > 0.25 or dx > 0.25:
            raise RuntimeError(
                'PSS.monodromy_twin: the %s twin converged to a DIFFERENT '
                'orbit than the seeding %s orbit (period Delta = %.2e, state '
                'Delta = %.2e, gate 0.25): that orbit is too poor to seed a '
                'monodromy twin -- the free-period Newton reached a spurious '
                'limit cycle.  Refine the grid so the %s state is a good '
                'seed, or set pss.monodromy = "native" to read the run\'s own '
                '(defective) monodromy rather than a wrong number.'
                % (method, getattr(self.par, 'method', '?'), dT, dx,
                   getattr(self.par, 'method', '?')))
        twin.monodromy = 'native'
        cache[method] = twin
        return twin

    def _state_twin(self):
        """`monodromy_twin`, except that a Nordsieck GLM takes its twin
        DRIVEN OR NOT -- for the Lyapunov surfaces (`_lyapunov_host`) and,
        on an OSCILLATOR, the small-signal ones (`_small_signal_host`).
        Everything else reads the GLM's own map on the state
        (`PSS._state_map`; Andreas, 2026-09-25: "Native for all but
        covariance").  Both exceptions are MEASURED:

        * the covariance: a white source over a GLM step reaches the state
          through its effective weights ``w = l^T B`` (GLM3 0.359, -0.0167,
          0.067, 0.591; GLM4 -26 .. +166), which no per-stage sampling can
          carry with positive variances, so a native injection is one shared
          sample per step -- first order;
        * an oscillator's small-signal response NEAR A HARMONIC: see
          `_small_signal_host`.

        `monodromy='native'` keeps the GLM's own map for both.
        """
        host = self.monodromy_twin()
        if (host is self and self._map_kind() == 'glm'
                and getattr(self, 'monodromy', 'radau') != 'native'
                and getattr(self, '_period_state', None) is not None
                and self.converged):
            return self._solve_twin(self.monodromy)
        return host

    def _small_signal_host(self):
        """The `PSS` whose map on the state `PAC.solve`, the adjoint rows and
        `pnoise` read: this one -- a Nordsieck GLM included, on its own map
        (`PSS._state_map`) -- except an OSCILLATOR solved with a GLM, which
        hands them to its twin (`_state_twin`: radau unless
        `monodromy='native'`).

        ⚠ MEASURED (2026-09-25), NOT ASSUMED.  A GLM's map on the state
        opens with its startup, which breaks the discrete phase symmetry:
        its unit multiplier sits ``eta = O(h^p)`` off 1 (van der Pol:
        glm3 2.4e-5 / 1.1e-6 at 60 / 120 points, glm2 7e-6 / 1.1e-6; radau
        3e-11).  The deflated solve returns the DISCRETE operator's answer
        (refined on it, dual-consistent: `PAC._deflated_solve`), which
        near a harmonic misses the physical one by ``eta / (2 pi r)``, `r`
        the offset in units of f0 -- glm3 at 60 points 0.4 % at 1e-3, 35 %
        at 1e-5, and below ``r ~ eta / 2 pi`` the pole is gone (the answer
        bounded).  A driven circuit has no such pole: there the GLM's own
        map is exact at its order (PAC on an LTI RLC against AC: 2.0 /
        3.0 / 4.1 per halving for glm2 / 3 / 4).  (The deflated answer
        UNREFINED is O(h^p) at every offset -- glm3 3.3e-5 at 60 points
        from 1e-3 to 1e-7 -- but forward and adjoint then agree only to
        O(eta); not built.)"""
        if self._map_kind() == 'glm' and getattr(self, 'autonomous', False):
            return self._state_twin()
        return self

    def _lyapunov_host(self):
        """The `PSS` the Lyapunov noise surfaces (`covariance`,
        `oscillator_covariance`) read.

        The covariance is propagated on the SAME orbit the Floquet quantities
        use, so this is just the monodromy twin: a trap/euler autonomous run
        hands its covariance to the twin (`TR-BDF2` by default, `gear` if
        `monodromy='gear'`) so both come from one consistent orbit; a
        gear/trbdf2 host is its own host.  TR-BDF2's per-step injection is
        built (`_lyapunov_pieces_trbdf2`, DAE-projected Van Loan), so there
        is no Gear-2 fallback -- the injection follows the chosen twin.
        A Nordsieck GLM hands them to its twin driven or not (`_state_twin`:
        radau unless `monodromy='native'`).
        """
        return self._state_twin()

    def _adjoint_host(self):
        """Host for the ADJOINT SIDEBAND noise surface (`pnoise`, via
        `adjoint_sideband_row`).

        The two-stage sideband fold IS built for TR-BDF2
        (`_sideband_forced_trbdf2`, the two-vector injected reverse pass that
        carries the source coupling through both stages), so this is the
        monodromy twin -- the same orbit the Floquet and Lyapunov surfaces
        use, no Gear-2 fallback.  `monodromy='gear'` still routes to the
        Gear-2 twin if asked.  A Nordsieck GLM reads its own map on the
        state (`PSS._state_map`; until 2026-09-25 a radau twin), except on
        an oscillator (`_small_signal_host`).
        """
        return self._small_signal_host().monodromy_twin()

    ## `grid_error`'s ceiling on a plausible OBSERVED order is the method's
    ## nominal order (`_nominal_order`), read off its integrator for every
    ## accepted name.  The table it replaced covered 6 of 12 -- esdirk43,
    ## the GLMs and the aliases `gear2` / `trapezoidal` silently took the
    ## generic range, the ceiling not applied.
    ## History: `doc/shooting_history.md`, `_AccuracyChecks.METHOD_ORDER`.

    def grid_error(self, evaluate, refine=2, levels=3, label=None):
        """How much of a scalar is DISCRETISATION rather than answer.

        Re-solves this circuit on a `refine`x finer grid, with every other
        argument identical to the original `solve()`, and reports how far
        `evaluate` moves.  `evaluate` takes a solved `PSS` and returns a
        float -- `lambda p: PAC(cir).diffusion_constant(p)`, a Floquet
        multiplier, a harmonic amplitude, the period.

        ⚠⚠ WHY THIS EXISTS RATHER THAN A PER-METHOD FORMULA: the floor is
        DISCRETISATION, linear in Q and a METHOD property (on the high-Q van
        der Pol reference at 240 points, `~1.8e-05 Q` for gear against
        `~7.0e-12 Q` for radau), and such constants belong to one fixture
        and grid -- a fitted predictor would extrapolate across a regime
        change.  Refining the actual circuit measures the actual number.

        ⚠⚠ WHY THREE GRIDS AND NOT TWO.  With `f_h = f + C h^p`, the
        two-grid difference `|C| h^p (1 - r^-p)` bounds the fine grid's
        error only under a single power law: a quantity mixing two
        discretisations of opposite sign changes sign, the difference
        collapses, and the estimate UNDER-STATES the error.  The third grid
        is the VALIDITY CHECK: from `d1 = |f_h - f_h/r|` and
        `d2 = |f_h/r - f_h/r^2|`,

            order = log(d1/d2) / log(r)

        is the order the circuit ACTUALLY shows, checked against the
        single-power-law assumption before `error = d2 / (r^order - 1)` is
        offered; `power_law=False` means READ `d2` AS A RAW CHANGE AND
        NOTHING MORE.

        ⚠ IT ESTIMATES THE GRID ERROR ONLY: it cannot see an error both
        grids share (a wrong stamp, tolerance convention or circuit).  A
        small `rel_error` says the grid is fine enough, NOT that the answer
        is right.

        `levels=3` (the default) costs two extra solves, at `r` and `r^2`
        times the points.  `levels=2` is the cheap two-grid difference with
        no order and no validity check -- use it only where the method's
        order on this problem is already known.

        Returns a dict with `values` (coarse to finest), `deltas`, `order`,
        `error` (of the FINEST value), `rel_error`, `power_law`, `refine`
        and `npts`.

        History: `doc/shooting_history.md`, `_AccuracyChecks.grid_error`.
        """
        args = getattr(self, '_solve_args', None)
        if args is None:
            raise RuntimeError(
                'PSS.grid_error: call solve() before grid_error() -- the '
                'refinement repeats THIS solve and there is nothing to '
                'repeat yet.')
        refine = int(refine)
        if refine < 2:
            raise ValueError(
                'PSS.grid_error: refine must be >= 2, got %r; a refinement '
                'that does not refine reports 0.0 and means nothing.'
                % (refine,))

        levels = int(levels)
        if levels not in (2, 3):
            raise ValueError('PSS.grid_error: levels must be 2 or 3, got %r'
                             % (levels,))
        ## ⚠ A CALLER-SUPPLIED GRID CANNOT BE REFINED BY A TIMESTEP.  `grid`
        ## fixes the sample fractions outright, so dividing `timestep` would
        ## change nothing and this would report a confident 0.0.
        if args.get('grid') is not None:
            raise ValueError(
                'PSS.grid_error: this solve used an explicit `grid`, whose '
                'fractions fix the samples regardless of `timestep`. The '
                'refinement would return the same grid and report 0.0 -- '
                'pass a refined `grid` and compare directly instead.')

        ## Same class, same parameters, same solve arguments -- only the
        ## timestep changes.  Fresh instances rather than re-solving `self`,
        ## so the caller's solved state survives the call.
        kv = {}
        for _p in self.parameters:
            try:
                kv[_p.name] = getattr(self.par, _p.name)
            except AttributeError:
                pass

        values = [float(evaluate(self))]
        for _k in range(1, levels):
            ## `self.irefnode` is an INDEX; `__init__` wants the node, and
            ## feeding the index back through `get_node_index` would not
            ## round-trip.
            twin = type(self)(self.cir, toolkit=self.toolkit,
                              irefnode=self.cir.nodes[self.irefnode], **kv)
            ## ⚠ AND THE SAME MONODROMY TWIN: `monodromy` is an attribute,
            ## not a Parameter, and the refined runs took the default --
            ## measured, a trap run under `monodromy='gear'` read the gear
            ## twin at level 0 and the default twin at levels 1 and 2
            twin.monodromy = getattr(self, 'monodromy', 'radau')
            sub = dict(args)
            sub['timestep'] = args['timestep'] / float(refine ** _k)
            twin.solve(**sub)
            if not getattr(twin, 'converged', False):
                raise RuntimeError(
                    'PSS.grid_error: the %dx refined solve did not converge, '
                    'so there is no comparison to report.'
                    % (refine ** _k,))
            values.append(float(evaluate(twin)))

        deltas = [abs(values[i + 1] - values[i])
                  for i in range(len(values) - 1)]
        _tiny = np.finfo(float).tiny
        order, power_law = None, None
        if levels == 3:
            d1, d2 = deltas
            ## ⚠ `d2 >= d1` means the sequence is NOT settling: either the
            ## error is not a single power law (two terms of opposite sign) or the
            ## finest grid has reached a roundoff floor.  Either way the
            ## Richardson step below would be arithmetic on noise.
            _sgn = ((values[1] - values[0]) * (values[2] - values[1]) > 0.0)
            if d2 > _tiny and d1 > d2:
                order = float(np.log(d1 / d2) / np.log(float(refine)))
                ## ⚠⚠ THE CEILING IS THE POINT OF THIS CHECK: a method cannot
                ## converge faster than its order, and an observed order well
                ## above it means two error terms nearly cancelled, so the
                ## estimate UNDER-states -- with same-signed deltas that a
                ## generic `0.5 <= order <= 8` range would accept.  The `+ 1.5`
                ## is not slack: on an AUTONOMOUS problem the period absorbs
                ## the leading frequency error, and `gear` (nominal 2)
                ## genuinely converges at ~3.  ⚠ AND IT IS THE HIGHER OF THE
                ## RUN'S METHOD AND ITS TWIN'S: an oscillator quantity of a
                ## trap run (its `c`) comes from the twin, and under the radau
                ## twin it converges at radau's 5.03 -- a trap-only ceiling
                ## (3.5) withheld a correct estimate.
                _hi = max(self._nominal_order(),
                          self._twin_nominal_order()) + 1.5
                power_law = bool(0.5 <= order <= _hi and _sgn)
            else:
                order, power_law = None, False
            if not power_law:
                warnings.warn(
                    'PSS.grid_error: the refinement does not follow a single '
                    'power law (changes %.3e then %.3e over %dx refinements, '
                    'implied order %s). The error estimate is WITHHELD -- '
                    'read the raw change and nothing more. This happens when '
                    'two error terms of opposite sign cancel at some grid, '
                    'which makes a two-grid difference UNDER-state the true '
                    'error, and when the finest grid has hit a roundoff '
                    'floor.'
                    % (d1, d2, refine,
                       ('%.2f' % order) if order is not None else 'none'),
                    RuntimeWarning, stacklevel=2)

        if power_law:
            err = deltas[-1] / (float(refine) ** order - 1.0)
        else:
            ## No validated order: report the raw change, which is what a
            ## two-grid call gets and is honest about being unbounded.
            err = deltas[-1]
        rel = err / max(abs(values[-1]), _tiny)
        _npts_c = int(round(args['period'] / args['timestep']))
        return {'values': values, 'deltas': deltas, 'order': order,
                'power_law': power_law, 'error': err, 'rel_error': rel,
                'refine': refine,
                'npts': [_npts_c * refine ** k for k in range(levels)],
                'label': label}

    ## The interpolant degree `warping_estimate` needs, per method: its ORDER
    ## must exceed the method's effective order (else a constant bias), and
    ## for a COLLOCATION method its DEGREE must exceed the stage count (else
    ## the spline lies inside the exactness class -- a cubic for Radau
    ## IIA(3) -- and the estimate reads a SILENT ~0).  For a non-collocation
    ## method (ESDIRK43) only the order clause is established.
    ## History: `doc/shooting_history.md`, `_AccuracyChecks.IDEC_DEGREE`.
    ## ⚠ THE SMALLEST ODD DEGREE ABOVE THE METHOD'S ORDER, at least 3: the
    ## rule every measured entry of the table this replaced followed (3 for
    ## the second-order methods, quintic for esdirk43, septic for radau),
    ## and it keeps the stage clause for every method here.  The table
    ## missed the GLMs: glm3 and glm4 fell back to a cubic, and glm4's
    ## estimate read -6.7e-7 against a true +2.8e-9 (van der Pol, 120
    ## points; +5.7e-10 quintic).  Even degrees trip the check.
    WARPING_CHECK_TOL = 0.05     # |half-grid / full-grid - 1| above this: the interpolant sets the reading

    def _idec_degree(self, method=None):
        """`warping_estimate`'s interpolant degree for `method` (this
        run's by default): the smallest odd degree above its order, at
        least 3."""
        p = self._nominal_order(method)
        k = (p + 1) if p % 2 == 0 else (p + 2)
        return max(3, k)

    def _twin_nominal_order(self):
        """The nominal order of the twin this run hands surfaces to -- a
        one-step LMM's oscillator surfaces (`monodromy_twin`), a GLM's
        covariance surfaces (`_state_twin`) -- or 0 when it hands none."""
        mono = getattr(self, 'monodromy', 'radau')
        if mono == 'native':
            return 0
        try:
            _own = self._integrator_for(getattr(self.par, 'method', 'euler'))
        except ValueError:
            return 0
        uses = (self._map_kind() == 'glm'
                or (not _own.carries_own_monodromy()
                    and getattr(self, 'autonomous', False)))
        return self._nominal_order(mono) if uses else 0

    def _nominal_order(self, method=None):
        """The method's nominal convergence order, as its integrator
        states it (`ORDER`, or a GLM's `order`) -- every accepted name,
        aliases included."""
        integ = self._integrator_for(
            method if method is not None
            else getattr(self.par, 'method', 'euler'))
        p = getattr(integ, 'ORDER', None)
        if p is None:
            p = getattr(integ, 'order')
        return int(p)

    def warping_estimate(self, periods=20, degree=None, check=True):
        """Estimate THIS solve's period (warping) error at ITS OWN grid, with
        no reference solution and no refinement -- by defect correction.

        A per-step local truncation estimate cannot see warping, a GLOBAL
        error; defect correction (Sickenberger, Weinmueller & Winkler,
        "Local Error Estimates for Moderately Smooth ODEs and DAEs", Part I,
        Sec. 1) estimates the global error directly:

          1. p(t)  -- a periodic spline through this solve's own grid values
                      (`_idec_degree`, or `degree`);
          2. r(t)  = C(p) p' + i(p) + u(t)   -- the DEFECT of the interpolant
                      against the circuit's own equations, T-periodic;
          3. the NEIGHBOURING problem  d/dt q(y) + i(y) + u(t) - r(t) = 0,
                      whose exact solution is p by construction, integrated as
                      a TRANSIENT with the SAME method at the SAME step for
                      `periods` periods (`Transient.solve` takes `-r` through
                      `provided_function`, an extra source term on every path);
          4. the phase lag of y against p, per period, by projection onto p':
                      tau_k = <p'.(y - p)> / <p'.p'>;  its slope against the
                      period count is the estimated period error.

        A transient and not a periodic solve, deliberately: a forcing at T
        fixes the period, and warping cannot present as a period change.

        ⚠⚠ THE INTERPOLANT IS THE LIMIT: the estimate is right only if the
        interpolant's defect error is o(h^{p+1}) (Part I's "only if").  On
        an orbit with an edge a few points wide the spline does not resolve
        it, and the estimate under-reads by a factor nothing announces.
        Trust it on smooth orbits; on edges refine until it converges in
        `periods` and grid, or use `grid_error`.  (A restructured f-value
        defect, Part I eq. 2.13, is not the fix: it is not monotone in the
        grid.)

        `check=True` (the gate) repeats the pass through every second sample
        of the same solution, the transient still at the solve's step:
        `check_ratio` = half-grid slope / full-grid slope, `trusted` =
        within `WARPING_CHECK_TOL` of 1, else a warning and the number still
        returned.  It doubles the cost.

        ⚠ The period reading needs an AUTONOMOUS solve; on a driven circuit
        `period_error` is None and the (bounded) lag series is still filled.
        On a DAE the components converge at different orders (H&W VI.7), so
        a component-wise exactness collapse is invisible in the scalar drift
        -- see `component_rms` and `lag_components`.  Index-2 is outside
        Part I's stated scope.  Cost: `periods` periods of transient at the
        working grid -- no refinement sweep, no analytic reference.

        Returns a dict: `period_error` (s, signed: positive = this solve's
        period is LONG), `ppm`, `lag` (per-period phase lag, s), `degree`,
        `periods`, `autonomous`, `component_rms` (RMS of y - p per unknown
        over the last period, the global error estimate in state space),
        `check_ratio` and `trusted` (the half-grid self-check above; both
        None with `check=False` or when a slope is not finite),
        `lag_components` (periods x unknowns: the same lag per component --
        every row's slope should equal `period_error`; a row at ~0 while
        the others agree is the per-component exactness-class collapse the
        DAE caveat names; a row of NaN is a component with no motion, e.g.
        a node pinned by a source).

        History: `doc/shooting_history.md`, `_AccuracyChecks.warping_estimate`.
        """
        import numpy as _np
        from scipy.interpolate import make_interp_spline
        if self.waveform is None:
            raise ValueError('warping_estimate needs a solved PSS -- call solve() first')
        T = float(self.period)
        times = _np.asarray(self.waveform[0], dtype=float)
        X = _np.asarray(self.waveform[1], dtype=float)         # (n, m)
        if abs(times[-1] - T) > 1e-12 * T:
            times = _np.r_[times, T]; X = _np.column_stack([X, X[:, 0]])
        X = X.copy(); X[:, -1] = X[:, 0]                        # close the orbit exactly
        method = str(self.par.method)
        k = int(self._idec_degree(method) if degree is None else degree)
        if X.shape[1] <= k + 1:
            raise ValueError('warping_estimate: %d points per period cannot carry a degree-%d '
                             'periodic spline' % (X.shape[1] - 1, k))
        cir, epar = self.cir, self.epar
        ## ⚠ `analysis='tran'`, on BOTH calls.  `Circuit.u(t)` evaluates a
        ## time function only when told which analysis is asking; without it
        ## every source VANISHES (zeros, DC value included), and the defect
        ## would omit the drive on a driven circuit.
        _u = lambda t: _np.asarray(cir.u(t, epar, analysis='tran'), dtype=float)
        u0 = _u(0.0)
        autonomous = all(_np.allclose(u0, _u(f * T)) for f in (0.37, 0.71))
        n_per = X.shape[1] - 1

        def _run(times_i, X_i):
            """One defect-correction pass: the periodic spline through
            (times_i, X_i), its defect, the neighbouring transient at the
            SOLVE's step T/n_per -- the method's error at ITS grid is what is
            read, and only the interpolant's grid differs between the main
            pass and the self-check below -- and the per-period lag."""
            p = make_interp_spline(times_i, X_i.T, k=k, bc_type='periodic')
            dp = p.derivative()
            def _defect_source(t):
                tt = t % T
                x = _np.asarray(p(tt), dtype=float); xd = _np.asarray(dp(tt), dtype=float)
                r = (_np.asarray(cir.C(x, epar), dtype=float) @ xd
                     + _np.asarray(cir.i(x, epar), dtype=float) + _u(t))
                return -r
            tr = self._new_transient(self._integrator_for(self.par.method))
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = tr.solve(tend=periods * T, x0=X_i[:, 0].copy(), timestep=T / n_per,
                               provided_function=_defect_source, fixed_timestep=True)
            ty = _np.asarray(res.sweep_values, dtype=float)
            Y = _np.asarray(res.x, dtype=float)                     # (n, steps+1)
            if Y.shape[0] != X_i.shape[0]:
                Y = Y.T
            P = _np.asarray(p(ty % T), dtype=float).T; dP = _np.asarray(dp(ty % T), dtype=float).T
            E = Y - P
            lag = []; lag_c = []
            for j in range(periods):
                sl = (ty >= j * T - 1e-12 * T) & (ty < (j + 1) * T - 1e-12 * T)
                num = float(_np.sum(dP[:, sl] * E[:, sl])); den = float(_np.sum(dP[:, sl] ** 2))
                lag.append(num / den if den > 0 else _np.nan)
            ## per component: the same projection restricted to one unknown.
                ## A phase shift moves every component by the same lag, so on a
                ## healthy estimate every row's slope equals the period error;
                ## a row reading ~0 while the others read the period error is
                ## the exactness-class collapse on THAT component (the DAE
                ## caveat), invisible in the scalar `lag` above.
                num_c = _np.sum(dP[:, sl] * E[:, sl], axis=1); den_c = _np.sum(dP[:, sl] ** 2, axis=1)
                ## a RELATIVE threshold: a node pinned by a source has a
                ## derivative of pure roundoff, and `den > 0` would print a
                ## ratio where NaN is meant.
                with _np.errstate(divide='ignore', invalid='ignore'):
                    lag_c.append(_np.where(den_c > 1e-20 * den_c.max(), num_c / den_c, _np.nan))
            lag = _np.asarray(lag); lag_c = _np.asarray(lag_c)
            ok = _np.isfinite(lag)
            slope = float(_np.polyfit(_np.arange(periods)[ok], lag[ok], 1)[0]) if ok.sum() >= 2 else _np.nan
            return slope, lag, lag_c, E, ty

        slope, lag, lag_c, E, ty = _run(times, X)
        ## SIGN: a positive lag means y is AHEAD of p; a LONG period makes y
        ## fall BEHIND, so the period error is minus the slope.
        period_error = -slope if autonomous else None
        last = ty >= (periods - 1) * T - 1e-12 * T
        component_rms = _np.sqrt(_np.mean(E[:, last] ** 2, axis=1))
        ## THE SELF-DIAGNOSTIC (Part I's "only if"): the
        ## estimate is the method's error only while the interpolant's own
        ## defect is asymptotically smaller than it, and then it does NOT
        ## depend on the interpolant: the same pass through EVERY SECOND
        ## sample of the same solution (the transient still at the solve's
        ## step) must read the same slope.  Where the interpolation error
        ## dominates -- an edge a few points wide -- the two passes disagree,
        ## and the reading is refused (`trusted=False`, a warning) instead
        ## of returned as a number wrong by a factor nothing announces.
        check_ratio = None; trusted = None
        if check:
            sub = _np.arange(0, n_per + 1, 2)
            if sub[-1] != n_per:
                sub = _np.r_[sub, n_per]
            if len(sub) > k + 1:
                slope2 = _run(times[sub], X[:, sub])[0]
                if _np.isfinite(slope) and _np.isfinite(slope2) and slope != 0.0:
                    check_ratio = float(slope2 / slope)
                    trusted = bool(abs(check_ratio - 1.0) <= self.WARPING_CHECK_TOL)
                    if not trusted:
                        warnings.warn('warping_estimate: the interpolant, not the method, sets '
                                      'this reading (half-grid pass / full-grid pass = %.3f); '
                                      'refine the grid until the two agree' % check_ratio)
        return dict(period_error=period_error,
                    ppm=(period_error / T * 1e6) if period_error is not None else None,
                    lag=lag, lag_components=lag_c, degree=k, periods=periods,
                    autonomous=autonomous, component_rms=component_rms,
                    check_ratio=check_ratio, trusted=trusted)

    ## Below this, a multiplier says the mode decays by six decades in one
    ## period and no stability question turns on it -- so parasitic
    ## contamination at that level is not worth a warning.  Used only to
    ## keep the warning off circuits where the WHOLE spectrum is numerical
    ## noise; it never changes a reported number.
    SPECTRAL_NOISE_FLOOR = 1e-6

    def _spectral_report(self, M):
        """Split a composed spectrum into physical multipliers and parasitics.

        A k-step method's composed monodromy carries `(k-1) m` PARASITIC
        roots beside the `m` physical Floquet multipliers, and `max |eig|`
        over that mixture is a stability verdict only while they stay small
        -- true for Gear-2 (`(1/3)^N`), not for a method whose spurious root
        sits nearer the unit circle.

        THE DISCRIMINATOR IS THE EIGENVECTOR'S BLOCK STRUCTURE, not the
        eigenvalue.  The composed map acts on the PAIR `(x_0, x_{-1})`:

          - a PHYSICAL mode follows the linearised ODE, so its two halves
            are one timestep apart on a smooth trajectory and
            `v_{-1} = e^{-lambda h} v_0 -> v_0` as `h -> 0`;
          - a PARASITIC mode is `r^n u` for the method's spurious root `r`,
            so `v_{-1} = u / r` -- three times `v_0` for BDF-2, minus it for
            a trapezoidal-like root -- and the halves differ by O(1)
            whatever `h` is.

        So `||v_{-1} - v_0||` (unit-norm eigenvector) is O(h) for a physical
        mode and O(1) for a parasitic one, and the split is by RANK, not by
        a threshold: the `m` smallest splits are the physical set by
        construction (a threshold finds no physical mode at all on a stiff
        circuit, where `lambda h ~ 40`).

        ⚠ THE COUNT IS AN ODE COUNT, AND MNA CIRCUITS ARE DAEs: an index-1
        system with `d = rank(C) < m` has `d` physical multipliers, `d`
        parasitic and `2(m - d)` STRUCTURAL ZEROS (Demir, IJCTA 28:163-185,
        2000: `Phi(t,s) = U(t) D(t-s) V(s) C(s)`), so on a real circuit
        `parasitic_roots` comes back identically zero and
        `floquet_multipliers` carries structural zeros; `rank(C)` also
        overcounts by one per index-2 constraint.  The factorisation's
        trailing `C(s)` has no ODE analogue -- revisit the split from there,
        not from the eigenvector heuristic (check it against the paper
        first).  Magnitude cannot separate the populations either: Gear-2's
        parasitic roots (~1e-95) are numerically structural zeros.

        ⚠ `spectral_radius` IS UNAFFECTED: the physical multipliers have the
        smallest splits, so the maximum over the first `m` is right; only
        the LABELS of the diagnostic arrays degrade -- also on a stiff
        circuit, where a parasitic root can swap with a stiff physical mode
        at the noise floor.  What the split buys is a method whose spurious
        root sits NEAR THE UNIT CIRCLE, where a maximum over the mixture
        would report the discretisation's artefact as the orbit's stability.

        Returns `(rho, physical, parasitic)`: the spectral radius over the
        PHYSICAL multipliers only, and both sets sorted by magnitude.  An
        `m x m` monodromy (any one-step method, the plain path) has no pairs
        and no parasitic roots, so everything in it is physical.  `None`
        gives `(None, None, None)` -- the matrix-free path forms no
        monodromy at all.

        History: `doc/shooting_history.md`, `_AccuracyChecks._spectral_report`.
        """
        if M is None:
            return None, None, None
        M = np.asarray(M)
        m = self.cir.n - 1
        try:
            ev, V = np.linalg.eig(M)
        except np.linalg.LinAlgError:                     # pragma: no cover
            return None, None, None

        if M.shape[0] != 2 * m:
            ## one-step method: the monodromy IS the physical map
            phys = np.sort(np.abs(ev))[::-1]
            return float(phys[0]), phys, np.array([])

        ## columns of `V` are unit-norm, so this needs no denominator and
        ## cannot divide by a vanishing block -- a mode living entirely in
        ## one half reads as O(1) here, which is what it is.
        split = np.linalg.norm(V[m:, :] - V[:m, :], axis=0)
        order = np.argsort(split)
        phys = np.sort(np.abs(ev[order[:m]]))[::-1]
        para = np.sort(np.abs(ev[order[m:]]))[::-1]
        rho = float(phys[0])
        ## ⚠ THE POINT AT WHICH THIS STOPS BEING BOOKKEEPING.  While the
        ## parasitic roots are 80 orders down, separating them changes no
        ## number and only documents why the maximum was safe.  Once one
        ## climbs to within a decade of the physical spectrum, the method's
        ## spurious roots are a real part of what the analysis reports and
        ## the user is entitled to know before reading a stability verdict.
        if len(para) and para[0] > 0.1 * rho and rho > self.SPECTRAL_NOISE_FLOOR:
            warnings.warn(
                'PSS: this method\'s PARASITIC roots are no longer '
                'negligible -- the largest is %.4g against a physical '
                'spectral radius of %.4g. A k-step method contributes '
                '(k-1)*m spurious roots to the composed monodromy, and '
                '`spectral_radius` now reports the maximum over the '
                'PHYSICAL multipliers only (separated by eigenvector block '
                'structure). Treat the separation as load-bearing here '
                'rather than cosmetic: check `floquet_multipliers` and '
                '`parasitic_roots` before drawing a stability conclusion.'
                % (para[0], rho), RuntimeWarning, stacklevel=3)
        return rho, phys, para
