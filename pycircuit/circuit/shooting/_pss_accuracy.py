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

        `self` for a driven circuit, when `self.monodromy == 'native'`, or
        when this run's own method already carries a second-order monodromy
        (`gear`/`trbdf2`).  Otherwise (an autonomous circuit solved with a
        one-step LMM, trap or euler) a twin `PSS` of the same circuit under
        the method `self.monodromy` names, solved once on the SAME grid from
        this orbit's converged state and cached.

        ⚠⚠ THE TWIN DEFAULTS TO TR-BDF2, MEASURED (2026-09-05), and Gear-2
        is one setting away, not retired.  The twin exists because trap's
        and euler's own monodromy is unusable on a limit cycle: trap's
        diverges with refinement (B16, below) and euler's is first order.
        Both TR-BDF2 and Gear-2 give a clean second-order twin, but TR-BDF2
        is more accurate on `lambda2` -- measured against exact references
        (`exp(A T)` on a linear oscillator; Abel's `exp(mu integral(1-v^2))`
        on van der Pol) it is 12-32x better at practical step counts, and
        the advantage GROWS with Q and with coarser grids -- the regime a
        real oscillator PSS sits in.  At Q=100 and 50 points/period the
        Gear-2 twin misreads Q by 22%, the TR-BDF2 twin by 0.27%.  The gap
        is a coarse-grid/high-Q effect, not a fixed factor: refine the grid
        or drop Q and both fall to the ordinary O(h^2) floor where the
        difference is single digits and can even favour Gear-2.  So
        `monodromy = 'trbdf2'` is the default, `'gear'` restores the former
        twin, `'native'` reads the run's own.

        ⚠⚠ THE B16 DECISION, TAKEN ON A MEASUREMENT (2026-09-05).  On the
        bias-sensitive oscillator (`vdp + 0.3 u^2`, exact `Q_lambda =
        5.908`, `c_true = 5.3703e-06`) trapezoidal's monodromy is unusable
        with EITHER opener: the default reads `Q_lambda` 11.1 / 28.4 / 63.9
        at 400/800/1600 points and DIVERGES with refinement, and
        `x0_unknown=True` reads 3086 / 12228 / 48699 -- a spurious
        multiplier at 1 (the one-step companion's parasitic mode), while
        the state and period are second order either way.  So "the most
        accurate" is not a choice between openers: the STATE keeps the
        method you asked for, and every monodromy-derived quantity -- `Q`,
        the PPV and everything built on it, the Floquet modes, the
        phase-noise surfaces -- comes from the twin on the same orbit.  The
        twin's period differs from this one's by O(h^2); its orbit is
        re-converged, not copied.

        ⚠⚠ THE B16 DECISION, TAKEN ON A MEASUREMENT (2026-09-05).  On the
        bias-sensitive oscillator (`vdp + 0.3 u^2`, exact `Q_lambda =
        5.908`, `c_true = 5.3703e-06`) trapezoidal's monodromy is unusable
        with EITHER opener: the default reads `Q_lambda` 11.1 / 28.4 / 63.9
        at 400/800/1600 points and DIVERGES with refinement, and
        `x0_unknown=True` reads 3086 / 12228 / 48699 -- a spurious
        multiplier at 1 (the one-step companion's parasitic mode), while
        the state and period are second order either way.  Gear-2 reads
        5.9094 / 5.9086 / 5.9084.  So "the most accurate" is not a choice
        between openers: the STATE keeps the method you asked for, and
        every monodromy-derived quantity -- `Q`, the PPV and everything
        built on it, the Floquet modes, the phase-noise surfaces -- comes
        from Gear-2 on the same orbit.  The twin's period differs from
        this one's by O(h^2); its orbit is re-converged, not copied.
        """
        mono = getattr(self, 'monodromy', 'trbdf2')
        if mono not in ('trbdf2', 'gear', 'native'):
            raise ValueError(
                "PSS.monodromy must be 'trbdf2', 'gear' or 'native', not %r"
                % (mono,))
        _integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        if (mono == 'native'
                or _integ.carries_own_monodromy()
                or not getattr(self, 'autonomous', False)):
            ## ⚠ SELF-SUFFICIENT METHODS TAKE NO TWIN, and the method says which
            ## it is (`carries_own_monodromy`): Gear-2 (a second-order native
            ## companion monodromy) and every stage method (TR-BDF2, Radau -- no
            ## opener seam, verified against the pencil).  This twin exists only
            ## because a one-step LMM's monodromy is first-order on a limit
            ## cycle (its manufactured opener is dropped to Euler and that seam
            ## sits in the period map); twinning a self-sufficient method would
            ## replace its own map with another on a re-converged orbit -- pure
            ## cost -- and hide the run's own spectrum.  They read `native`
            ## regardless of `monodromy`; the knob governs which twin trap/euler
            ## borrow.
            return self
        if getattr(self, '_period_state', None) is None or not self.converged:
            return self
        twin = self._solve_twin(mono)
        self._monodromy_twin = twin
        return twin

    ## The iteration budget of a monodromy twin's re-solve (Andreas,
    ## 2026-09-23).  A twin is a POLISH, not a cold solve: it is seeded at
    ## this run's converged state on the same grid, and from a good seed it
    ## converges in 4-13 iterations (measured on B16's van der Pol).  It used
    ## to inherit the caller's `maxiterations` -- 300 in B16 -- and from a
    ## poor seed (Euler at 400 points, its orbit 55 % off) both twins refuse by
    ## NOT converging, so each ran its whole budget times the solve's retry
    ## ladder: 971 + 968 traversals, 790 s, to say "no".  Capped at this, and
    ## a twin that hits the cap without converging WARNS before it refuses.
    TWIN_MAXITER = 40

    def _solve_twin(self, method):
        """A converged twin `PSS` of this circuit under `method`, re-solved
        once on the SAME grid from this orbit's converged state, cached per
        method.  The shared machinery behind `monodromy_twin` (which picks
        the method by `self.monodromy`) and `_lyapunov_host` (which forces
        `gear`, because the noise-injection surfaces cannot use a TR-BDF2
        twin yet).  Raises if the re-solve does not converge.
        """
        from .pss import PSS    # imported when called: pss.py imports this module
        cache = self._twins
        if method in cache:
            return cache[method]
        solved, x0, xm1, times, hs, T, x0_unknown = self._period_state
        kw = dict(self._solve_kwargs)
        twin = PSS(self.cir, toolkit=self.toolkit, irefnode=None,
                   method=method, reltol=self.par.reltol,
                   iabstol=self.par.iabstol, vabstol=self.par.vabstol)
        hs = np.asarray(hs, dtype=float)
        ## the same grid: its fractions when it is not uniform, else the
        ## uniform step (a one-step plain state can carry an `hs` whose
        ## sum is a trial period, so the fractions are the safe object)
        nonuniform = float(hs.max() / hs.min()) > 1.0 + 1e-9
        grid = (hs / float(hs.sum())) if nonuniform else None
        x0r = np.asarray(x0, dtype=float)[:self.cir.n - 1]
        ## ⚠ THE STEP COUNT MUST SURVIVE `solve`'s `int(period / timestep)`.
        ## `T / (T / N)` is not N in floating point: measured on B16's
        ## fixture, T = 6.730731946457316 gives 399.99999999999994, so the
        ## "same grid" twin ran on 399 points against the state's 400 (found
        ## when `phase_rule='reselect'` moved T in its 15th digit).  Half a
        ## step of slack makes the floor land on N for every T.
        ##
        ## ⚠ AND THE TWIN TAKES THE DEFAULT PHASE RULE RATHER THAN THIS RUN'S.
        ## It is seeded at the converged state so that it lands on the SAME
        ## point of the SAME orbit -- which is exactly what the agreement
        ## check below tests -- and `phase_rule='reselect'` re-chooses the
        ## pinned coordinate and lands on ANOTHER phase (measured; see
        ## `solve`).  Inheriting it would move the twin off the orbit point
        ## whose monodromy was asked for.
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

        ## ⚠⚠ THE TWIN MUST HAVE CONVERGED TO THE SAME ORBIT IT WAS SEEDED
        ## ON, and a MORE ROBUST twin makes this check load-bearing rather
        ## than paranoid.  Measured: from a poor seed (euler at 400 pts, its
        ## orbit 55% off) the Gear-2 twin fails to converge -- LOUD -- but
        ## the TR-BDF2 twin, being more robust, CONVERGES to a SPURIOUS limit
        ## cycle and reports `Q = 1.97` against the exact 5.91 with no error.
        ## Improving the method degraded safety: the failure moved from a
        ## refusal to a plausible wrong number.  So the twin's converged
        ## orbit is checked against the seed it was handed: two convergent
        ## methods on the SAME limit cycle agree on period and entering state
        ## to O(h^p) -- measured 6e-6 / 1e-4 (trbdf2) and 3e-5 / 2e-2 (gear)
        ## on a good seed -- while the spurious jump above sits at 2.17 /
        ## 0.91.  The 0.25 gate is ~12x above the worst good case and ~3.6x
        ## below the spurious one; it is set from the (universal, tiny)
        ## good-case agreement, not the (fixture-dependent) failure size, so
        ## it transfers.  The definitive test is refinement (a spurious orbit
        ## does not survive h/2); this cheap consistency check is the
        ## conservative stand-in -- it REFUSES a too-poor seed rather than
        ## risk trusting it, which is the safe direction.
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
        """
        return self.monodromy_twin()

    def _adjoint_host(self):
        """Host for the ADJOINT SIDEBAND noise surface (`pnoise`, via
        `adjoint_sideband_row`).

        The two-stage sideband fold IS built for TR-BDF2
        (`_sideband_forced_trbdf2`, the two-vector injected reverse pass that
        carries the source coupling through both stages), so this is the
        monodromy twin -- the same orbit the Floquet and Lyapunov surfaces
        use, no Gear-2 fallback.  `monodromy='gear'` still routes to the
        Gear-2 twin if asked.
        """
        return self.monodromy_twin()

    ## Nominal convergence order per `method`, for `grid_error`'s ceiling on
    ## a plausible OBSERVED order.  Sourced from this file's own measured
    ## records rather than from the literature: trap/gear/theta second order,
    ## TR-BDF2 measured at 4.01x/4.01x/4.00x per halving (exact `O(h^2)`),
    ## Radau IIA(3) at 31.50x/31.74x (`O(h^5)`, theoretical 32), euler first.
    ## ⚠ An unlisted method falls back to a generic range and the ceiling is
    ## not applied -- add it here rather than letting it default silently.
    METHOD_ORDER = {'euler': 1, 'trap': 2, 'gear': 2, 'theta': 2,
                    'trbdf2': 2, 'radau': 5}

    def grid_error(self, evaluate, refine=2, levels=3, label=None):
        """How much of a scalar is DISCRETISATION rather than answer.

        Re-solves this circuit on a `refine`x finer grid, with every other
        argument identical to the original `solve()`, and reports how far
        `evaluate` moves.  `evaluate` takes a solved `PSS` and returns a
        float -- `lambda p: PAC(cir).diffusion_constant(p)`, a Floquet
        multiplier, a harmonic amplitude, the period.

        ⚠⚠ WHY THIS EXISTS RATHER THAN A PER-METHOD FORMULA.  The floor of
        this stack is DISCRETISATION, it grows LINEARLY IN Q, and it is a
        METHOD property: measured on the analytic high-Q van der Pol
        reference, the relative error in the diffusion constant at 240
        points per period is

            Q      gear        trap        radau
             100   1.79e-03    2.63e-05    6.97e-10
             500   9.02e-03    1.32e-04    3.48e-09
            1000   1.82e-02    2.63e-04    6.97e-09

        (trap's column corrected 2026-09-14: it read 1.49e-06 / 1.04e-04 /
        2.36e-04, the period-normalisation defect's values.)

        i.e. `~1.8e-05 Q` for gear against `~7.0e-12 Q` for radau -- SIX
        ORDERS at the same cost per step.  Those constants are real but they
        belong to THAT fixture at THAT grid: `gear` converges at `O(h^3)` on
        an autonomous problem for `Q >= 5` and at `O(h^2)` at `mu = 1`, so a
        shipped predictor built from them would extrapolate a fitted constant
        across a regime change (roadmap D.0y).  Refining the actual circuit
        measures the actual number instead, and needs no calibration.

        ⚠⚠ WHY THREE GRIDS AND NOT TWO.  With `f_h = f + C h^p`, two grids
        give `|f_h - f_h/r| = |C| h^p (1 - r^-p)`, which over-states the fine
        grid's own error `|C|(h/r)^p` by `r^p - 1` -- an upper bound, and a
        tempting place to stop.  **IT IS NOT SAFE, AND THIS STACK CONTAINS A
        COUNTEREXAMPLE.**  A quantity mixing TWO discretisations -- an
        `O(h^3)` error plus an `O(h^2)` one of opposite sign -- changes sign:
        the two terms cancel, the two-grid difference collapses, and the
        estimate UNDER-STATES the true error by 3.6x (measured: change
        1.15e-06 against a true 4.06e-06 at 240 points).  ⚠ That quantity
        was `diffusion_constant` under `trap` until 2026-09-14 -- the twin's
        integral over trap's own period, a DEFECT now fixed -- and the
        validity check below is what refused it.  A bound that fails silently
        where the error is interesting is worse than none.

        So the third grid is not extra confidence, it is the VALIDITY CHECK.
        From `d1 = |f_h - f_h/r|` and `d2 = |f_h/r - f_h/r^2|`,

            order = log(d1/d2) / log(r)

        is the order the circuit ACTUALLY shows, and it is checked against the
        single-power-law assumption before the error estimate built on it is
        offered.  Measured orders on that fixture: `gear` 2.94 (its `O(h^3)`
        autonomous rate), `radau` ~5, `trap` 3.02 (its twin's), and the
        mixed quantity failing the check exactly where it cancels.  `error` is then `d2 / (r^order - 1)`, and
        `power_law=False` means READ `d2` AS A RAW CHANGE AND NOTHING MORE.

        ⚠ AND IT IS AN ESTIMATE OF THE GRID ERROR ONLY.  It cannot see an
        error both grids share -- a wrong stamp, a wrong tolerance
        convention, a mis-specified circuit.  A small `rel_change` says the
        grid is fine enough; it does NOT say the answer is right.  That is
        the same trap `null_residual_amplification` documents one screen up,
        and it is worth stating twice.

        `levels=3` (the default) costs two extra solves, at `r` and `r^2`
        times the points.  `levels=2` is the cheap two-grid difference with
        no order and no validity check -- use it only where the method's
        order on this problem is already known.

        Returns a dict with `values` (coarse to finest), `deltas`, `order`,
        `error` (of the FINEST value), `rel_error`, `power_law`, `refine`
        and `npts`.
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
                ## ⚠⚠ THE CEILING IS THE POINT OF THIS CHECK, AND A GENERIC
                ## RANGE IS NOT ENOUGH.  A method cannot converge faster than
                ## its order; an observed order well above it means two error
                ## terms nearly cancelled at this grid, which makes the
                ## deltas shrink faster than the error and the estimate
                ## UNDER-state.  MEASURED: a quantity mixing two
                ## discretisations (a twin's `c` over the host's period --
                ## what `diffusion_constant` under `trap` computed until the
                ## 2026-09-14 fix) shows an apparent order of 6.45 at 120
                ## points -- monotone, same-signed deltas, nothing else
                ## suspicious -- while its estimate under-states the true
                ## error by 300x.  A plain `0.5 <= order <= 8` range
                ## ACCEPTS that case; the ceiling below rejects it.
                ## ⚠ The `+ 1.5` allowance is not slack: on an AUTONOMOUS
                ## problem the period is an unknown that absorbs the leading
                ## frequency error, so `gear` (nominal 2) genuinely converges
                ## at 3.01 here.  Without the allowance this would reject the
                ## shipped default method on its own reference fixture.
                _nom = self.METHOD_ORDER.get(
                    str(getattr(self.par, 'method', '')).lower())
                _hi = 8.0 if _nom is None else (_nom + 1.5)
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

    ## B7: the interpolant degree the defect-correction estimate needs, per
    ## method.  Two-part rule, MEASURED 2026-09-08 (doc/pss_roadmap_260902.md,
    ## B7's gate): the interpolant's DEGREE must exceed the method's stage
    ## count -- a cubic spline lies INSIDE Radau IIA(3)'s collocation
    ## exactness class (degree s = 3), so the neighbouring problem is solved
    ## EXACTLY and the estimate is 1e-10 ppm against a true 6e-3 -- and its
    ## ORDER must exceed the method's effective order, or a constant bias
    ## remains (a quintic against radau's measured 6.1 left 4.9 % at every
    ## grid; a septic gave 1.0001 / 0.9999 / 0.9998).  Cubic reproduced trap
    ## and TR-BDF2 to 0.9996 -> 1.0000.
    ## ⚠ THE STAGE-COUNT CLAUSE IS A COLLOCATION PROPERTY (peer, measured the
    ## same day): ESDIRK43 has SIX stages and is not a collocation method, so
    ## a quintic fails the clause literally -- and quintic and septic AGREE
    ## through the stack (0.9998 / 1.0000 at 100 pts, 1.0000 / 1.0001 at 200).
    ## For a non-collocation method only the ORDER clause is established;
    ## written as "degree > stage count" the rule would over-constrain every
    ## DIRK ever added.  The failure the clause guards against is SILENT (a
    ## clean small number), which is why it was measured rather than argued.
    WARPING_CHECK_TOL = 0.05     # |half-grid / full-grid - 1| above this: the interpolant sets the reading
    IDEC_DEGREE = {'euler': 3, 'trap': 3, 'gear': 3, 'theta': 3,
                   'trbdf2': 3, 'esdirk43': 5, 'radau': 7}

    def warping_estimate(self, periods=20, degree=None, check=True):
        """Estimate THIS solve's period (warping) error at ITS OWN grid, with
        no reference solution and no refinement -- by defect correction.

        B7's answer, MEASURED (2026-09-08).  A per-step local truncation
        estimate cannot see an accumulating period error because warping is a
        GLOBAL error; defect correction (Sickenberger, Weinmueller & Winkler,
        "Local Error Estimates for Moderately Smooth ODEs and DAEs", Part I,
        Sec. 1) estimates the global error directly:

          1. p(t)  -- a periodic spline through this solve's own grid values
                      (`IDEC_DEGREE[method]`, or `degree`);
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

        Measured on A10's van der Pol (Q = 1e4), estimate / true period error
        (true = T_h - T_ref, radau at 3200 points), numpy prototype:

            trap   / cubic    0.9996  0.9999  1.0000  1.0000   (100..800 pts)
            trbdf2 / cubic    0.9999  1.0000  1.0000
            radau  / cubic    0.0000  0.0000  0.0000   (25..50 pts) -- INSIDE the
                                                       exactness class: see IDEC_DEGREE
            radau  / quintic  1.0490  1.0494  1.0495   -- a constant bias where the
                                                       orders tie (6 vs 6.1)
            radau  / septic   1.0001  0.9999  0.9998

        Controls: with RADAU solving the neighbouring problem of TRAP's defect
        the estimate is 0.0000 -- the drift is the METHOD's error, not the
        defect's; a LINEAR interpolant gives 0.03 / 2.3 / 3.4 -- the
        interpolant-order wall from below.

        Through THIS method (the stack, same fixture, reference radau at 3200
        points, 2026-09-08): trap 1.0002 at 400 and 200 pts; radau septic
        1.0003 / 1.0002 at 50 / 35 pts and CUBIC 0.0001 (the exactness-class
        zero, reproduced); esdirk43 quintic 0.9998 / 1.0000 and septic
        1.0000 / 1.0001 at 100 / 200 pts -- so for a non-collocation method
        the order clause alone is established, and the stage-count clause is
        a collocation property -- and Part I sec 1.1 says why: "one of the
        most attractive features of the IDeC procedure is, that its fixed
        point is a certain superconvergent COLLOCATION solution", so the
        exactness class the stage-count clause guards against is a
        collocation object by construction (docs session, 2026-09-09; one
        family at two degrees on one fixture, so a mechanism, not a proof
        that the clause is harmless in general).  ⚠ The first stack gate's driven control came
        back `autonomous=True`: `Circuit.u(t)` evaluates its time functions
        only when told `analysis='tran'`, and without it every source
        VANISHES (zeros, DC value included) -- so that control ran against a
        circuit with NO source at all; fixed at both call sites; with the
        flag the driven van der
        Pol returns `autonomous=False`, `period_error=None`, and a bounded
        lag series, as it must.

        ⚠⚠ THE INTERPOLANT IS THE LIMIT (measured 2026-09-08): on a
        relaxation oscillator with a comparator edge a few points wide the
        estimate reads 0.09 of the true period error at 200 points per
        period and 0.65 at 400 -- uniformly in every component, so not a
        collapse: the septic spline does not resolve the edge and the
        defect is interpolation error, not the method's.  Part I's own
        scope is "moderately smooth"; a relaxation orbit at PSS grids is
        outside it, and the number returned is then wrong by a factor that
        nothing in it announces.  Trust it on smooth orbits (1.000 to four
        digits on the van der Pol, index 1 and 2); on an orbit with edges,
        refine until the estimate converges in `periods` and grid, or use
        `grid_error`.
        ⚠ PART I READ THROUGH (docs session, 2026-09-09): its motivation is
        this domain -- "we are especially motivated by applications in
        electrical circuit simulation, where the models often contain data
        with poor smoothness" -- so "moderately smooth" is the case the
        paper was built for, not a clause this is outside of.  Its Remark
        2.9 names a failure with the SAME SIGN as the edge underestimate:
        the local estimates assume the leading term `c_i h^(p+1) x^(p+1)`
        does not vanish, and "at least in case of oscillatory solutions,
        there always exist time points where the derivative x^(p+1)
        vanishes ... our error estimates will tend to UNDERESTIMATE the
        true size of the error" (footnote: the third derivative vanishes
        where the curvature is extremal; remedy: assume C^(p+2) and match
        the next coefficient with an auxiliary scheme).  That remark is
        stated for the LOCAL-error route of their section 2; this method is
        the GLOBAL route of section 1 (Zadunaisky), and whether the global
        route inherits it is not established.  ⚠ Ruled out here by the
        h-scaling: Remark 2.9's mechanism is keyed on isolated zeros of
        x^(p+1), whose aggregate effect is roughly h-INDEPENDENT, while the
        edge reading improved 7.2x for a 2x grid (0.09 -> 0.65) -- that is
        interpolation error, as stated above.  Where Remark 2.9 would bite
        is a step controller built on defect correction; the paper hands
        the fix.  Cost lead, not worked out: section 1.1's cheap variant
        runs the high-order method once and a cheap LOW-order method twice
        (original and neighbouring problem) -- a different substitution
        from the radau-on-trap's-defect control that zeroed the estimate.
        ⚠ THE LITERATURE'S ANSWER IS STRUCTURAL, NOT "REFINE" (docs session,
        Part I p. 9, READING-LOG 2.165).  The gate this instrument should
        test before returning a number is Part I's own "only if": the
        estimate is asymptotically correct ONLY IF the interpolant's defect
        error is o(h^{p+1}) -- asymptotically SMALLER than the truncation
        error it is meant to reveal; on the comparator edge it is not, and
        the number is wrong by a factor nothing announces.  And the fix for
        a non-smooth orbit is to form the defect as a WEIGHTED SUM OF
        f-VALUES with an auxiliary scheme sharing the base scheme's
        left-hand side, so the solution terms cancel identically (their eq.
        2.13, an extra factor h) -- not a higher-degree interpolant of the
        solution, which is exactly the construction this one uses.  Scope:
        their construction is the LOCAL error of an LMM; whether it
        transfers to a period functional is unproven.  THE GATE IS BUILT
        (2026-09-08, `check=True`): the same pass through every second
        sample of the same solution, transient still at the solve's step;
        `check_ratio` = half-grid slope / full-grid slope, `trusted` =
        within `WARPING_CHECK_TOL` (5 %) of 1, else a warning and the number
        still returned.  Measured: van der Pol 1.0000 (radau and trap, 50
        and 100 points); the relaxation orbit 0.0056 / 0.72 at 200 / 400
        points under radau and 0.41 at 200 under trap -- the cases that
        read 0.09 / 0.65 of the truth are refused, the smooth case accepted
        with four orders of margin.  ⚠ The prediction "below 0.5 at 400"
        was wrong (0.72): the ratio approaches 1 as the edge resolves, so
        the tolerance is the gate, not the ratio's distance from 0.  The
        restructured (f-value) defect (Part I eq. 2.13; for trap the Milne
        device) was GATED and REFUTED as the edge fix (2026-09-08): smooth
        0.9990 / 0.9998, but on the edge orbit 0.06 / 0.61 / 1.29 at 200 /
        400 / 800 points against the spline route's 0.04 / 0.29 / 0.63 --
        faster with the grid and NOT monotone, so a reading near 1 is
        indistinguishable from a wrong one; the paper's own remedy is mesh
        adaptation.  Not built.  Cost of the check: it doubles the call (a
        second `periods`-long transient).
        ⚠ Scope and limits.  The period reading needs an AUTONOMOUS solve;
        on a driven circuit the lag is bounded (entrained) and `period_error`
        is returned as None with the per-period lag series still filled.
        The prototype ran on a 2-state ODE; on a DAE the differential and
        algebraic components converge at different orders (H&W VI.7), so the
        interpolant threshold binds per component and a component-wise
        exactness collapse would be invisible in this scalar phase drift --
        `component_rms` is returned so a caller can look.  Index-2 is outside
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
        k = int(self.IDEC_DEGREE.get(method, 3) if degree is None else degree)
        if X.shape[1] <= k + 1:
            raise ValueError('warping_estimate: %d points per period cannot carry a degree-%d '
                             'periodic spline' % (X.shape[1] - 1, k))
        cir, epar = self.cir, self.epar
        ## ⚠ `analysis='tran'`, on BOTH calls.  `Circuit.u(t)` evaluates a
        ## time function only when told which analysis is asking (`VS.u`:
        ## `elif analysis in timedomain_analyses`); without it every source
        ## VANISHES -- the else-branch returns zeros, and even the DC value
        ## lives inside the gated branch (`timedomain_analyses = ('dc',
        ## 'tran')`).  The first gate's driven control -- an `ISin` on the
        ## van der Pol -- came back `autonomous=True` for exactly that reason,
        ## and the defect would have omitted the drive on a driven circuit.
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
                ## derivative of pure roundoff (measured 1e-32 rms), and
                ## `den > 0` let it print a ratio of 96 where NaN was meant.
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
        ## THE SELF-DIAGNOSTIC (Part I's "only if", built 2026-09-08): the
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

        RECORDED SCOPE ITEM 3.  A k-step method turns an m-dimensional
        system into a k*m-dimensional discrete one, so the composed
        monodromy's spectrum carries `(k-1) m` PARASITIC roots beside the
        physical Floquet multipliers.  `max |eig|` over that mixture is only
        a stability verdict while the parasitic roots stay small -- which
        for Gear-2 they emphatically do (`(1/3)^N`, ~1e-95 at 200 points)
        and for a method whose spurious root sits nearer the unit circle
        they would not.  This separates them instead of hoping.

        THE DISCRIMINATOR IS THE EIGENVECTOR'S BLOCK STRUCTURE, not the
        eigenvalue.  The composed map acts on the PAIR `(x_0, x_{-1})`:

          - a PHYSICAL mode follows the linearised ODE, so its two halves
            are one timestep apart on a smooth trajectory and
            `v_{-1} = e^{-lambda h} v_0 -> v_0` as `h -> 0`;
          - a PARASITIC mode is `r^n u` for the method's spurious root `r`,
            so `v_{-1} = u / r` -- three times `v_0` for BDF-2, minus it for
            a trapezoidal-like root -- and the halves differ by O(1)
            whatever `h` is.

        So `||v_{-1} - v_0||` (against a unit-norm eigenvector) is O(h) for
        a physical mode and O(1) for a parasitic one.  MEASURED, and it is
        the h-scaling that makes it a prediction rather than a story: on the
        phase circuit the physical ratio falls 0.1281 -> 0.0316 when the
        grid goes from 50 to 200 points -- a factor of 4.05 for a factor of
        4 in `h` -- while the parasitic ratios sit at 1.0 to 10.  On the
        Q=20 RLC the parasitic ratio is 1.9997 against the 2.0 that BDF-2's
        `v_{-1} = 3 v_0` predicts exactly.

        ⚠ THE MODE COUNT HERE IS AN ODE COUNT AND THE OBJECT IS A DAE, and
        the difference is structural rather than an off-by-`k`.  Demir
        (IJCTA 28:163-185, 2000) gives the DAE monodromy as

            Phi(t,s) = U(t) D(t-s) V(s) C(s)

        with `D = diag[exp(mu_1 (t-s)), ..., exp(mu_d (t-s)), 0, ..., 0]`
        for `d = rank(C)`: "equation (19) has k = n - m Floquet multipliers
        that are 0", and on a real circuit "there are also eigenvalues
        exactly equal to 0 due to the ALGEBRAIC EQUATIONS in the MNA
        formulation".  So the `m - rank(C)` structural zeros are the
        theory's, not an artefact -- which is why `parasitic_roots` comes
        back identically zero on every MNA circuit tried here.  ⚠ AT INDEX
        1 ONLY (the docs session, checked against the paper 2026-09-09:
        "We assume that the DAEs we are dealing with are index-1").
        `rank(C)` is the differential dimension at index 1 and OVERCOUNTS
        by one per index-2 constraint: measured on `floquet_modes`, an
        index-1 tank and an index-1 tank + R node give modes = rank(C) = 2,
        an L-I cutset and a C-V loop give rank(C) = 3 with 2 modes -- the
        code returns the true count; it is the formula that stops where
        Demir says it does.

        ⚠ AND THE FACTORISATION CARRIES A TRAILING `C(s)` WITH NO ODE
        ANALOGUE (where `C = I` and it disappears).  A DAE monodromy is not
        simply a product of state-transition blocks, so an ODE-shaped
        count does not merely miscount -- it describes a different object.
        Anyone revisiting this split should start there and not from the
        eigenvector heuristic below.  Relayed from the docs session's read;
        check it against the paper before building on it.

        ⚠ THE SPLIT IS BY RANK, NOT BY A THRESHOLD, and that was measured
        into the design rather than chosen.  A threshold of 0.25 was tried
        first and returned NO physical modes at all on a stiff RC ladder --
        `lambda h ~ 40` there, so every mode's halves differ by O(1) and the
        classifier called the entire spectrum parasitic, handing back a
        `spectral_radius` of `None` where the old code said 6e-15.  A
        `k`-step method on `m` states has EXACTLY `m` physical multipliers
        and `(k-1) m` spurious ones -- that is structural -- so the `m`
        smallest splits are the physical set by construction, and the
        question of where to put a cut never arises.

        ⚠ THE COUNT IS AN ODE COUNT, AND MNA CIRCUITS ARE DAEs.  This
        splits `2m` eigenvalues as `m` physical and `m` parasitic, which is
        right for an ODE.  An index-1 MNA system with `d = rank(C) < m` has
        `d` physical multipliers, `d` parasitic ones and `2(m - d)`
        STRUCTURAL ZEROS from the algebraic variables -- so on a real
        circuit both arrays are mislabelled: measured on the Q=20 resonator
        (`m = 4`, `rank(C) = 2`), `parasitic_roots` comes back identically
        zero and `floquet_multipliers` carries two structural zeros beside
        the two real multipliers.

        ⚠ `spectral_radius` IS UNAFFECTED, which is why this is recorded
        rather than re-engineered.  The physical multipliers have the
        SMALLEST block split by construction, so they are always inside the
        first `m`, and the maximum over that set is the right number --
        0.97531 on that circuit, against the analytic 0.9753.  What is
        unreliable is the LABELLING of the diagnostic arrays.  And it cannot
        be fixed by magnitude either: Gear-2's true parasitic roots are
        `(1/3)^N`, about 1e-95, which is numerically indistinguishable from
        a structural zero -- so on this method the two populations cannot be
        told apart at all, by any test, and saying so is the honest
        position.

        ⚠ ON A STIFF CIRCUIT THE LABELS MAY STILL BE WRONG, and it does not
        matter: when the physical modes are themselves stiff, a parasitic
        root can have the smaller split and swap places with one.  Every
        mode involved then has `|mu|` at the noise floor, so the RADIUS is
        unaffected -- it is the labels, not the number, that degrade.  What
        this buys is the case that motivated the item: a method whose
        spurious root sits NEAR THE UNIT CIRCLE, where the physical modes
        are well resolved, the splits separate cleanly, and taking a
        maximum over the mixture would report the discretisation's own
        artefact as the orbit's stability.

        Returns `(rho, physical, parasitic)`: the spectral radius over the
        PHYSICAL multipliers only, and both sets sorted by magnitude.  An
        `m x m` monodromy (any one-step method, the plain path) has no pairs
        and no parasitic roots, so everything in it is physical.  `None`
        gives `(None, None, None)` -- the matrix-free path forms no
        monodromy at all.
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
