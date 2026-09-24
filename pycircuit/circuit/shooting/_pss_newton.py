"""The outer solves: the free-period Newton (with its stall diagnosis) and the
matrix-free Newton.
"""
import numpy as np
import warnings
import pycircuit.circuit.analysis as analysis


class _ShootingNewton(object):
    """The outer solves: the free-period Newton (with its stall diagnosis)
    and the matrix-free Newton.  A theme of `PSS` (see `pss.py`)."""

    ## Below this fraction of the seed, a solved period is the trivial
    ## root rather than an orbit.  Deliberately loose: a real fundamental
    ## reached from a seed one decade high is ~0.1 of it, and the trivial
    ## root lands 15 orders down, so nothing sits near this line.
    DEGENERATE_PERIOD_FACTOR = 1e-6
    ## The equilibrium test above: a returned state whose DC residual is
    ## within this factor of `iabstol` is an equilibrium, not an orbit.  1e3
    ## because the shooting Newton's own residual sits at `abstol`, and an
    ## orbit's DC residual at t = 0 is a CIRCUIT-scale current (measured:
    ## the van der Pol at 2 V has |i + u| ~ 1 A there against 1e-27 for the
    ## collapsed state -- 27 orders apart, so the factor is not delicate).
    TRIVIAL_ORBIT_FACTOR = 1e3

    def _free_period_solve(self, func, z0, abstol, xtol, reltol, maxiter,
                           seed_period, solver=None):
        """Solve a free-period system, with its degenerate root named.

        ⚠ `T = 0` IS A REGULAR ROOT OF EVERY AUTONOMOUS SHOOTING SYSTEM.
        `x0 - phi_T(x0)` vanishes identically at `T = 0`, and the phase
        condition does not exclude it -- it constrains `x0`, not the
        period -- so Newton reaches it from any seed below the fundamental
        and the run returns a period of ~1e-18 with no orbit in it.

        Measured on BOTH autonomous elements in the tree, so it is a
        property of the formulation and not of any circuit: from a 1e-4
        seed against a 1e-3 fundamental, Gear-2 returned -1.5e-20 on the
        quadrature element and 3.9e-19 on the scalar `Idtmod`, and
        trapezoidal raised a bare `LinAlgError` from three seeds of five as
        its Jacobian went singular on the way down.

        Neither outcome is a silent wrong answer -- the collapse reports
        `converged = False` (⚠ ENFORCED HERE, by demoting `ier`; asserting it
        in prose was not enough -- see the note at the demotion) and the
        exception is loud -- but both told the
        user nothing about the cause, and the generic non-convergence
        advice ("raise maxiterations") is wrong for it: no number of
        iterations reaches a fundamental from below.
        """
        try:
            if solver is None:
                z, info, ier, mesg = analysis.fsolve(
                    func, z0, maxiter=maxiter, reltol=reltol, abstol=abstol,
                    xtol=xtol, toolkit=self.toolkit, full_output=True,
                    line_search=True, floor_detect=True)
            else:
                ## ⚠ THE MATRIX-FREE ROUTE COMES THROUGH HERE TOO, so the
                ## trivial-root diagnosis below covers it.  Routing it around
                ## this wrapper would have lost the one message that makes an
                ## autonomous collapse readable.
                z, info, ier, mesg = solver(z0, abstol, xtol, reltol, maxiter)
        except np.linalg.LinAlgError as exc:
            raise np.linalg.LinAlgError(
                'PSS: the free-period Jacobian went singular while solving '
                'for an autonomous period seeded at %.6g s (%s). The usual '
                'cause is a seed BELOW the fundamental: `T = 0` solves the '
                'periodicity condition identically, so the iteration is '
                'drawn to it and the Jacobian degenerates on the way. Seed '
                'at or above the expected period -- a short transient and '
                'the interval between two output recurrences is the usual '
                'way to get one. A LinAlgError raised INSIDE the traversal '
                'is wrapped here too: a stage matrix going singular after a '
                'Newton step left the orbit reads the same -- check the '
                'period column first (2026-09-20: it omitted a constant '
                'source vector and threw x0 to 2.8e6 on the first step).'
                % (seed_period, exc)) from exc

        T = float(z[-1])
        ## ⚠⚠ THE SECOND TRIVIAL ROOT, found 2026-09-08 while gating
        ## `warping_estimate` on an index-2 oscillator: the guard below
        ## catches `T -> 0`, and an autonomous solve has ANOTHER root that it
        ## cannot see -- the EQUILIBRIUM, `x(t) = x_dc`, which is periodic at
        ## EVERY `T`.  Radau seeded 10 % below the fundamental on the
        ## index-1 van der Pol and on the index-2 fixture returned
        ## amplitude 0.0000 (state 1e-27) at a period near the seed with
        ## `converged = True`; trapezoidal on the same seed failed honestly.
        ## `T` stays finite, so the period test passes, and the periodicity
        ## residual is exactly zero because the equilibrium IS periodic --
        ## at every `T`, so it is a whole LINE of roots in `(x0, T)`, not a
        ## point (peer's sharpening): that is why the residual is exactly
        ## zero rather than small, why no residual-based guard could have
        ## caught it, and why the DC residual does in one evaluation.
        ## The test that sees it is the DC residual of the returned state:
        ## an orbit has `C x' != 0` somewhere at t = 0, so `i(x) + u` is far
        ## from zero there; an equilibrium has it at solver tolerance.
        trivial_orbit = False
        try:
            xr = np.asarray(z[:-1], dtype=float)
            irn = self.irefnode
            xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            r_dc = (np.asarray(self.cir.i(xf, self.epar), dtype=float).ravel()
                    + np.asarray(self.cir.u(0.0, self.epar, analysis='dc'),
                                 dtype=float).ravel())
            r_dc = np.delete(r_dc, irn)
            tol = float(getattr(self.par, 'iabstol', 1e-12))
            trivial_orbit = bool(np.abs(r_dc).max() <= self.TRIVIAL_ORBIT_FACTOR * tol)
        except Exception:
            trivial_orbit = False
        if trivial_orbit:
            ier = 5
            mesg = ('collapsed onto the EQUILIBRIUM (a trivial orbit, periodic '
                    'at every T) at T = %.6g s from a seed of %.6g s' % (T, seed_period))
            warnings.warn(
                'PSS: this autonomous solve returned an EQUILIBRIUM, not an '
                'orbit: the state at t = 0 satisfies the DC equations to '
                '%.1e (max |i(x) + u|), so the periodicity residual is zero '
                'at ANY period and the solver reported success at T = %.6g s '
                'from a seed of %.6g s. The basin of the equilibrium is '
                'entered from a seed below the fundamental (measured: 10 %% '
                'low under radau); seed at or above the expected period, or '
                'from a transient that is already on the orbit. '
                '`converged` is False.' % (np.abs(r_dc).max(), T, seed_period),
                RuntimeWarning, stacklevel=3)
        elif not np.isfinite(T) or abs(T) < self.DEGENERATE_PERIOD_FACTOR * abs(
                seed_period):
            ## ⚠⚠ THE COLLAPSE MUST BE DEMOTED HERE, and for two turns of this
            ## record it was not.  The docstrings above and on `solve` both
            ## asserted "the collapse reports `converged = False`" -- and
            ## NOTHING ENFORCED IT.  `self.converged` is `(_ier == 1)` and
            ## nothing else, while `T = 0` is a REGULAR root: the solver
            ## reaches it cleanly and reports success, so Gear-2 returned
            ## `T = 5.42e-18` with `converged = True` on a circuit with no
            ## orbit in it.  The warning fired correctly the whole time, which
            ## is exactly what made this survive -- a reader who checks the
            ## documented flag instead of catching warnings got `True`.
            ##
            ## Demoting `ier` rather than assigning `self.converged` is
            ## deliberate: all three autonomous call sites already feed this
            ## return value into `self.converged`, so one demotion covers the
            ## plain, solved-history and matrix-free paths, and any future
            ## path inherits it by construction.  `ier = 5` is `fsolve`'s
            ## "not making good progress" code -- the closest existing
            ## meaning, and already handled everywhere `ier` is read.
            ier = 5
            mesg = ('collapsed onto the trivial root T = %.6g s from a seed '
                    'of %.6g s' % (T, seed_period))
            warnings.warn(
                'PSS: this autonomous solve collapsed onto the TRIVIAL root, '
                'returning a period of %.6g s from a seed of %.6g s. `T = 0` '
                'satisfies `x0 - phi_T(x0) = 0` identically and the phase '
                'condition does not exclude it (it constrains x0, not the '
                'period), so a seed below the fundamental is drawn there. '
                'The returned waveform is not a periodic steady state. '
                'Raising maxiterations will not help -- seed at or above the '
                'expected period instead; a short transient and the interval '
                'between two output recurrences gives one. The literature '
                'remedy for this basin is the PROBE technique -- a periodic '
                'voltage source that feeds the oscillator until its own '
                'current reaches zero, widening the basin so a Newton is '
                'less likely to fall into the DC solution (Bizzarri et al., '
                '"Probe Based Shooting Method ..."); it is not implemented '
                'here, and it widens the basin rather than removing the '
                'seed dependence.'
                % (T, seed_period), RuntimeWarning, stacklevel=3)
        if (ier != 1 and not trivial_orbit and solver is None
                and np.isfinite(T)
                and abs(T) >= self.DEGENERATE_PERIOD_FACTOR * abs(seed_period)):
            self._diagnose_lmm_free_period_stall(func, z, info)
        return z, info, ier, mesg

    def _diagnose_lmm_free_period_stall(self, func, z, info=None):
        """Say WHY a multistep free-period solve stalled, when it is not the
        iteration count.

        ⚠⚠ MEASURED 2026-09-14 (peer report, reproduced).  On a weakly limited
        LC oscillator (second Floquet multiplier 0.99, set by the limiting),
        Gear-2's autonomous shooting solve stalls at a residual FLOOR:
        converged quadratically at 925..1600 points per period, never at 900
        or 800.  Every property that would point at the solver was ruled
        out -- the residual is a pure function of the unknowns; damped,
        Armijo (20 halvings) and Levenberg-Marquardt steps from the stall all
        stop at the floor; tighter inner tolerances change nothing; the grid
        is exactly uniform; a scan of the Jacobian's two weakest directions
        finds no lower point.  The discrete periodic solution of the
        solved-history system `(x_0, x_{-1}, T)` CEASES TO EXIST below a grid
        threshold: `sigma_min(J)` at the root falls linearly toward it (to
        zero near 904 points), the discrete second multiplier stays at 0.990,
        and the vanishing direction lies ~99 % in the `x_{-1}` block.  The
        threshold grows as the multiplier approaches 1 (no convergence at
        6400 points at 0.999).  Trapezoidal stalls too; `radau` -- no history,
        no manufactured opening -- converged in 7-10 evaluations at 800
        points at both multipliers.

        So no number of iterations, damping or tolerance fixes THAT stall, and
        the generic advice would send the caller the wrong way.  This measures
        what it can at the last iterate (one residual evaluation, only on a
        failed solve) and says so.  Stage methods are skipped: they did not
        show this.

        ⚠⚠ THE "ITERATIONS DO NOT HELP" CLAIM WAS OVER-BROAD AND IS NARROWED
        (2026-09-16, peer report + reproduced here).  It was measured for
        GEAR-2's solved-history stall and then written as though it held for
        every multistep free-period solve.  It does not: a trapezoidal solve at
        multiplier 0.9 converges with a bigger budget on the reporter's
        oscillator (200 iterations) and at the DEFAULT 25 on the van der Pol
        fixture here.  What IS measured at 0.99, on the reporter's tank and
        rebuilt independently: trapezoidal's residual RISES with budget --
        6.765e-07 at 25 iterations, 9.011e-07 at 200 -- and the analytic step
        is UPHILL against the true derivative, so the line search's halvings
        cannot improve it.  `radau` converged on that same fixture in 25
        iterations (period 6.2831863e-07, amplitude to five digits).

        ⚠ AND THE DISCRIMINATOR IS NOW REPORTED RATHER THAN GUESSED:
        `analysis.fsolve` counts the iterations whose step it could not improve
        (`infodict['ls_unimproved']`) and this message names the count.  Zero
        means the solve was descending and a bigger budget is worth trying;
        non-zero means the budget would repeat an uphill direction.

        ⚠ `x0_unknown=True` IS A DIAGNOSTIC HERE, NOT A FIX.  In that frame the
        analytic Jacobian IS the derivative (worst entry 6e-07 against 1.000 in
        the default frame, both measured by central differences on the
        reporter's fixture) and the step becomes a descent direction -- but the
        solve still does not converge there: `||F||` plateaus at 3.069e-05,
        identical at 25 and 200 iterations.  The frame explains the behaviour;
        it does not rescue the case.
        """
        try:
            integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
            if integ.is_stage_method():
                return
            _F, J = func(z)
            J = np.asarray(J, dtype=float)
            _U, s, Vt = np.linalg.svd(J)
            v = Vt[-1]
            m = self.cir.n - 1
            hist = (float(np.linalg.norm(v[m:2 * m]))
                    if J.shape[0] == 2 * m + 1 else None)
            cond = float(s[-1] / max(s[0], 1e-300))
            method = getattr(self.par, 'method', '?')
        except Exception:
            return
        ## the line search's own verdict on this solve -- see the docstring
        ls = int((info or {}).get('ls_unimproved', 0))
        if ls:
            ls_note = (
                'The line search could not improve %d of the steps taken here: '
                'the full step and four halvings each left the residual no '
                'better than the point they started from, and the step was '
                'taken anyway because there is no other candidate. That is a '
                'direction UPHILL against the true derivative, which halving '
                'cannot cure, so a larger budget would repeat it -- measured '
                'on a weakly limited LC tank at multiplier 0.99, trapezoidal '
                'goes 6.765e-07 at 25 iterations to 9.011e-07 at 200. ' % ls)
        else:
            ls_note = (
                'Every step here was one the line search could improve, so '
                'the solve was descending when the budget ran out: raising '
                'maxiterations is worth trying before anything else. ')
        common = ls_note + (
            "method='radau' (the default) converges where these do not -- 7-10 "
            'evaluations at 800 points on the oscillator above, and 25 '
            'iterations on a weakly limited LC tank at multiplier 0.99 '
            '(2026-09-16) where trapezoidal never did -- so use it, or refine '
            'the grid. ⚠ Whether MORE ITERATIONS help depends on the circuit, '
            'and this message over-claimed until 2026-09-16: gear-2 stalling '
            'with the solved-history signature below is a residual floor that '
            'iterations, damping and tolerances do not move, but a trapezoidal '
            'solve at multiplier 0.9 converges with a 200-iteration budget on '
            'one circuit and at the default 25 on another. The count above is '
            'the discriminator. ⚠ `x0_unknown=True` is a DIAGNOSTIC here, not '
            'a fix: it makes the Jacobian the actual derivative (worst entry '
            '6e-07 against 1.000) and the step a descent direction, but on '
            'that 0.99 fixture the residual plateaus at 3.069e-05 and the '
            'solve still does not converge.')
        if hist is not None and hist > 0.5:
            msg = ('PSS: this autonomous solve (method=%r) stopped at a '
                   'residual floor, and the Jacobian at the last iterate '
                   'shows the signature of a solved-history stall: its '
                   'weakest direction (sigma_min/sigma_max = %.1e) lies %.0f %% '
                   'in the entering-history block x_{-1}. On a weakly damped '
                   'orbit (second Floquet multiplier near 1) the two-step '
                   'discrete periodic solution exists only above a grid '
                   'threshold that grows as the multiplier approaches 1 '
                   '(measured: about 905 points per period at 0.99, above 6400 '
                   'at 0.999). ' % (method, cond, 100.0 * hist)) + common
        else:
            msg = ('PSS: this autonomous solve (method=%r) did not converge. '
                   'On a weakly damped oscillator (second Floquet multiplier '
                   'near 1) the multistep free-period solves stall at a '
                   'residual floor that iterations and tolerances do not move '
                   '(measured for gear and trapezoidal at 0.99 on 800 points; '
                   'sigma_min/sigma_max at the last iterate here: %.1e). '
                   % (method, cond)) + common
        warnings.warn(msg, RuntimeWarning, stacklevel=4)

    ## How hard GMRES is asked to solve, relative to the shooting tolerance.
    ## An inexact Newton only needs the step accurate enough not to spoil the
    ## outer convergence; measured k is 2-12 on circuits whose `I - M`
    ## clusters at 1 (the fast modes decay over a period, leaving the slow
    ## ones), so k tracks the number of SLOW MODES, not m.
    KRYLOV_TOLERANCE_FACTOR = 1e-2
    ## ⚠ THE BUDGET IS A CHOICE AND SCIPY'S UNITS ARE A TRAP: `maxiter` counts
    ## RESTART CYCLES, not matvecs, so the pair multiplies. 200 x 20 is far
    ## more than a clustered system needs; a circuit that exceeds it does not
    ## cluster, and the answer is the dense path, not a bigger budget.
    KRYLOV_RESTART = 200
    KRYLOV_MAX_CYCLES = 20

    def _matrix_free_newton(self, build, z0, abstol, xtol, reltol, maxiter):
        """The Newton loop every matrix-free system shares.

        `build(z)` returns `(F, matvec)` for the current iterate -- one
        trajectory pass, then a linear operator that never forms its matrix.
        Written once because the four systems differ ONLY in those two
        things: the plain path's `I - M`, the solved-history path's `2m`
        pair, and the bordered autonomous versions of each (one builder
        since 2026-09-23, `_mf_build` in `solve`).

        MEASURED (moved here from the solved-history driver, 2026-09-23).
        The dense path builds the `2m x 2m` Jacobian and factors it once per
        iteration; here the same iteration runs on a matvec, so the
        `2m`-column propagation never happens.  Measured against the dense
        path on the RC ladder, single-threaded, k=12:

              m     dense traversal   trajectory + 12 matvecs   speedup
             40             0.0843                    0.1025      0.82x
            110             0.2366                    0.2175      1.09x
            242             0.7503                    0.5378      1.40x
            502             3.4709                    1.5457      2.23x
           1002            20.1143                    5.5512      3.62x

        -- 82-87% of the predicted ceiling, and a LOSS at m=40, which the
        ceiling said too.

        ⚠ THOSE ARE TRAVERSAL FIGURES AND THE END-TO-END SOLVE GAINS LESS.
        A `solve` also does its setup, the replay that builds the waveform
        and the DFT, none of which this touches, and matrix-free spends an
        extra Newton iteration (below).  Measured end to end, same circuits:

              m    dense (iters)      matrix-free (iters)     speedup
            242      2.113 s (2)            1.557 s (2)        1.36x
            502      9.255 s (2)            6.131 s (3)        1.51x
           1002     52.402 s (2)           24.636 s (3)        2.13x

        Quote whichever answers the question being asked, and say which it
        is; 2.23x and 1.51x at m=502 are both true and are not the same
        measurement.

        ⚠ THE CONVERGENCE TEST IS NOT BIT-IDENTICAL TO `analysis.fsolve`'s,
        and it cannot be.  `fsolve` scales its residual test by
        `|J| . |x|`, an ELEMENTWISE absolute value of the Jacobian, which no
        matrix-free method has.  The substitute here is `|x| + |M x| + |F|`,
        one extra matvec per iteration.

        ⚠ AND IT IS NOT PROVABLY THE STRICT DIRECTION.  This docstring first
        claimed the substitute was a LOWER bound on `fsolve`'s scale, so
        that the test could only ever be stricter.  That is FALSE: at
        `M = I` the true scale `|I - M| . |x|` is zero while the substitute
        is `2|x|`, so the substitute is the LARGER one there, and at `M = 0`
        they are equal.  Neither dominates the other in general.

        What is measured, on the RC ladder at m=242/502/1002: the two paths
        agree on the converged waveform to 1.1e-16 and on the converged/not
        verdict, and matrix-free takes ONE MORE Newton iteration at m>=502
        (3 against 2) -- so it is stricter in practice here, and still wins
        on wall time while doing 50% more traversals.  One circuit is not a
        proof of direction, and this is the first thing to check if the two
        paths ever disagree on convergence.
        """
        import scipy.sparse.linalg as spla
        z = np.asarray(z0, dtype=float).copy()
        n = len(z)
        ier, mesg, xdiff = 2, 'No convergence', None
        ## the "stopped moving and still failing its step test" signature, as
        ## `analysis.fsolve(floor_detect=True)` records it: COUNTED, never
        ## acted on -- the iteration below is what it was
        _floor_q, step_floor = [], None
        for _i in range(maxiter):
            F, mv = build(z)

            def _mv(v, _f=mv):
                return _f(v)

            J = spla.LinearOperator((n, n), matvec=_mv, dtype=float)
            xdiff, info = spla.gmres(
                J, -F, rtol=self.KRYLOV_TOLERANCE_FACTOR * reltol,
                restart=min(n, self.KRYLOV_RESTART),
                maxiter=self.KRYLOV_MAX_CYCLES)
            ## ⚠ THE INNER SOLVE'S VERDICT IS NOT DISCARDED.  It used to be,
            ## and a Krylov breakdown then surfaced as the generic outer
            ## 'No convergence' with nothing naming the cause -- in a file
            ## whose whole standard is that a failure says what happened
            ## (`T = 0`, the trivial root, the singular free-period
            ## Jacobian).  An unconverged GMRES makes `xdiff` a direction
            ## the Newton has no reason to trust, so the outer loop is told
            ## to stop rather than iterate on it.
            if info != 0:
                warnings.warn(
                    'PSS: the matrix-free inner solve did not converge at '
                    'outer iteration %d -- GMRES returned info=%d (%s) on a '
                    '%d-unknown system, after at most %d matvecs, each a '
                    'full replay of the period. The Newton step it returned '
                    'is not a direction worth iterating on, so this solve '
                    'stops here and reports not-converged. Measured k on '
                    'well-behaved circuits is 2-12 because `I - M` clusters '
                    'at 1; needing more than %d means this system does not '
                    'cluster, and the dense path (matrix_free=False) is the '
                    'reliable answer for it.'
                    % (_i, info,
                       'breakdown' if info < 0 else 'iteration limit',
                       n, min(n, self.KRYLOV_RESTART) * self.KRYLOV_MAX_CYCLES,
                       min(n, self.KRYLOV_RESTART) * self.KRYLOV_MAX_CYCLES),
                    RuntimeWarning, stacklevel=3)
                return z, {}, 2, 'No convergence (inner Krylov solve failed)'
            z_new = z + xdiff

            ## `|J| . |x|` is not available without the matrix; see
            ## this docstring for what this substitute is and is not.
            I_scale = np.abs(z_new) + np.abs(mv(z_new)) + np.abs(F)
            conv_x = np.all(np.abs(xdiff)
                            < reltol * np.maximum(np.abs(z_new), np.abs(z))
                            + xtol)
            conv_f = np.all(np.abs(F) < reltol * I_scale + abstol)
            z = z_new
            if conv_x and conv_f:
                ier, mesg = 1, 'Success'
                break
            if conv_f:
                _tolx = reltol * np.maximum(np.abs(z_new), np.abs(z)) + xtol
                _dx = np.abs(xdiff)
                _k = int(np.argmax(_dx / _tolx))
                _floor_q.append((float(_dx[_k] / _tolx[_k]), _k,
                                 float(_dx[_k]), float(_tolx[_k])))
            else:
                _floor_q = []
            if len(_floor_q) >= 3 and any(
                    _floor_q[-j_][0] >= _floor_q[-j_ - 1][0] for j_ in (1, 2)):
                _w = max(_floor_q[-3:])
                if step_floor is None:
                    step_floor = dict(since=_i + 1 - len(_floor_q))
                step_floor.update(index=_w[1], step=_w[2], tol=_w[3],
                                  ratio=_w[0])
            elif not _floor_q:
                step_floor = None
        return z, {'step_floor': step_floor}, ier, mesg
