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
    ## The equilibrium test of `_free_period_solve`: a returned state whose
    ## DC residual is within this factor of `iabstol` is an equilibrium, not
    ## an orbit.  Not delicate: the shooting Newton's own residual sits at
    ## `abstol`, while an orbit's DC residual at t = 0 is a CIRCUIT-scale
    ## current, many orders above it.
    ## History: `doc/shooting_history.md`, `_ShootingNewton.TRIVIAL_ORBIT_FACTOR`.
    TRIVIAL_ORBIT_FACTOR = 1e3

    def _free_period_solve(self, func, z0, abstol, xtol, reltol, maxiter,
                           seed_period, solver=None, defer_diagnosis=False):
        """Solve a free-period system, with its degenerate root named.

        A collapse onto a trivial root sets ``info['collapsed']``.  With
        `defer_diagnosis`, a stall's diagnosis is not emitted but left as a
        callable in ``info['stall_diagnosis']``, for a caller whose next
        stage may still converge (`solve`'s state-event stage).

        ⚠ `T = 0` IS A REGULAR ROOT OF EVERY AUTONOMOUS SHOOTING SYSTEM.
        `x0 - phi_T(x0)` vanishes identically at `T = 0`, and the phase
        condition does not exclude it -- it constrains `x0`, not the
        period -- so Newton reaches it from any seed below the fundamental
        and the run returns a period of ~1e-18 with no orbit in it.

        It is a property of the formulation, not of any circuit.  The
        collapse is demoted here to `converged = False`, with a warning that
        names the cause, and a Jacobian going singular on the way down is
        re-raised with the same diagnosis: the generic advice ("raise
        maxiterations") is wrong for it -- no number of iterations reaches a
        fundamental from below.

        History: `doc/shooting_history.md`, `_ShootingNewton._free_period_solve`.
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
        ## ⚠⚠ THE SECOND TRIVIAL ROOT: the EQUILIBRIUM, `x(t) = x_dc`, is
        ## periodic at EVERY `T` -- a whole LINE of roots in `(x0, T)` -- so
        ## `T` stays finite, the period test passes, and the periodicity
        ## residual is exactly zero: no residual-based guard can catch it
        ## (radau seeded 10 % low returns it with `converged = True`).  The
        ## DC residual of the returned state does, in one evaluation: an
        ## orbit has `C x' != 0` somewhere at t = 0, so `i(x) + u` is far
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
            if isinstance(info, dict):
                info['collapsed'] = True
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
            ## ⚠⚠ THE COLLAPSE MUST BE DEMOTED HERE: `self.converged` is
            ## `(_ier == 1)` alone and `T = 0` is a REGULAR root the solver
            ## reports as success, so a warning alone leaves `converged =
            ## True`.  Demoting `ier` (not assigning `self.converged`) covers
            ## all three autonomous call sites -- plain, solved-history,
            ## matrix-free -- and any future path by construction; `ier = 5`
            ## is `fsolve`'s "not making good progress", handled everywhere.
            ier = 5
            if isinstance(info, dict):
                info['collapsed'] = True
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
            if defer_diagnosis and isinstance(info, dict):
                info['stall_diagnosis'] = (
                    lambda: self._diagnose_lmm_free_period_stall(func, z, info))
            else:
                self._diagnose_lmm_free_period_stall(func, z, info)
        return z, info, ier, mesg

    def _diagnose_lmm_free_period_stall(self, func, z, info=None):
        """Say WHY a multistep free-period solve stalled, when it is not the
        iteration count.

        On a weakly damped oscillator (second Floquet multiplier near 1) a
        multistep free-period solve can stall at a residual FLOOR.  Gear-2's
        solved-history system `(x_0, x_{-1}, T)` has no discrete periodic
        solution below a grid threshold that grows as the multiplier nears
        1 (`sigma_min(J)` -> 0, the vanishing direction almost wholly in
        `x_{-1}`): no iterations, damping or tolerance fix that; `radau`
        converges there.  This measures what it can at the last iterate (one
        residual evaluation, failed solves only; stage methods skipped).

        ⚠ Whether MORE ITERATIONS help depends on the circuit (a trapezoidal
        residual can fall or RISE with the budget, its step then UPHILL), so
        the discriminator is reported: `analysis.fsolve` counts the steps
        the line search could not improve (`infodict['ls_unimproved']`) and
        the message names it -- zero means a bigger budget is worth trying,
        non-zero that it would repeat an uphill direction.
        ⚠ `x0_unknown=True` is a diagnostic here, not a fix: the Jacobian
        becomes the true derivative, and the solve still plateaus.

        History: `doc/shooting_history.md`, `_ShootingNewton._diagnose_lmm_free_period_stall`.
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
        ## ⚠ ON A CIRCUIT WITH STATE EVENTS THE LIKELIER CAUSE IS THE SWITCH,
        ## NOT THE DAMPING: across a crossing sharper than the grid the
        ## unstaged map moves with where the crossing falls in its step.
        ## Measured: gear on the comparator oscillator stalls at 350 and 400
        ## points with a second multiplier of 0.002 -- and converges staged.
        try:
            _ev = (self.cir.state_events()
                   if hasattr(self.cir, 'state_events') else [])
        except Exception:                                      # noqa: BLE001
            _ev = []
        if _ev:
            msg += (' ⚠ This circuit declares %d state event(s): across a '
                    'switch sharper than the grid the UNSTAGED map moves with '
                    'where the crossing falls in its step, and its Newton can '
                    'fail on that whatever the damping (measured: gear on a '
                    'comparator oscillator, second multiplier 0.002). With '
                    '`state_events=True` (the default) the crossings become '
                    'Newton unknowns, and the stage runs even from a first '
                    'stage that did not converge.' % len(_ev))
        warnings.warn(msg, RuntimeWarning, stacklevel=4)

    ## How hard GMRES is asked to solve, relative to the shooting tolerance.
    ## An inexact Newton only needs the step accurate enough not to spoil the
    ## outer convergence.  `I - M` clusters at 1 (the fast modes decay over a
    ## period), so k tracks the number of SLOW MODES, not m (2-12 measured).
    ## History: `doc/shooting_history.md`, `_ShootingNewton.KRYLOV_TOLERANCE_FACTOR`.
    KRYLOV_TOLERANCE_FACTOR = 1e-2
    ## ⚠ THE BUDGET IS A CHOICE AND SCIPY'S UNITS ARE A TRAP: `maxiter` counts
    ## RESTART CYCLES, not matvecs, so the pair multiplies. 200 x 20 is far
    ## more than a clustered system needs; a circuit that exceeds it does not
    ## cluster, and the answer is the dense path, not a bigger budget.
    KRYLOV_RESTART = 200
    ## the relative residual below which an inner solve that missed its
    ## Krylov target still gives the Newton its step (`_matrix_free_newton`)
    KRYLOV_ACCEPT_RESIDUAL = 1e-6
    ## the matrix-free Newton's line search, as `analysis.fsolve`'s: four
    ## halvings reach 1/16 of the step
    LINE_SEARCH_MAX_HALVINGS = 4
    KRYLOV_MAX_CYCLES = 20

    def _matrix_free_newton(self, build, z0, abstol, xtol, reltol, maxiter,
                            line_search=True):
        """The Newton loop every matrix-free system shares.

        `build(z)` returns `(F, matvec)` for the current iterate -- one
        trajectory pass, then a linear operator that never forms its matrix
        -- or `(F, matvec, d)` with `d` per-unknown COLUMN SCALES: the Krylov
        solve is then on ``J diag(d)`` and the step ``diag(d) y``.  (The
        event stage's bordered system needs them: its crossing and period
        columns are ``dx/dtheta`` and ``dx/dT``, up to 1e6 against the unit
        columns of ``I - M``, and GMRES could not reach its tolerance on the
        raw operator -- condition 3.7e6 on the comparator oscillator, 16
        scaled.)

        `line_search`: `analysis.fsolve`'s -- the full step tried first and
        kept when it lowers ``||F||``, else halved up to four times, the
        accepted trial's `build` carried into the next iteration (so a
        converging solve costs what the undamped loop costs).  A trial the
        builder cannot evaluate (a `ValueError`, a failed inner step) counts
        as uphill.  ON, as the dense path's `fsolve(line_search=True)` is
        in every shooting stage (until 2026-09-25 the matrix-free Newton
        took full steps: the unstaged solve of the PWM loop failed under
        every method, radau included, and a first staged step from stage
        1's crossings moved two of them by a whole period).
        Written once because the four systems differ ONLY in those two
        things: the plain path's `I - M`, the solved-history path's `2m`
        pair, and the bordered autonomous versions of each (one builder,
        `_mf_build` in `solve`).

        The dense path builds the `2m x 2m` Jacobian and factors it once per
        iteration; here the same iteration runs on a matvec, so the
        `2m`-column propagation never happens.  It LOSES at small m and wins
        increasingly with m (on the RC ladder at m=1002, 3.6x per traversal
        and 2.1x end to end: setup, replay and DFT are untouched and
        matrix-free can spend an extra Newton iteration).  Traversal and
        end-to-end figures are different measurements; say which is quoted.

        ⚠ THE CONVERGENCE TEST IS NOT BIT-IDENTICAL TO `analysis.fsolve`'s,
        and it cannot be.  `fsolve` scales its residual test by
        `|J| . |x|`, an ELEMENTWISE absolute value of the Jacobian, which no
        matrix-free method has.  The substitute here is `|x| + |M x| + |F|`,
        one extra matvec per iteration.

        ⚠ AND IT IS NOT PROVABLY THE STRICT DIRECTION: at `M = I` the true
        scale `|I - M| . |x|` is zero while the substitute is `2|x|`, and at
        `M = 0` they are equal, so neither dominates.  On the RC ladder the
        two paths agree on the converged waveform to 1.1e-16 and on the
        verdict, matrix-free taking one more Newton iteration -- one circuit
        is not a proof of direction, and this is the first thing to check if
        the two paths ever disagree on convergence.

        History: `doc/shooting_history.md`, `_ShootingNewton._matrix_free_newton`.
        """
        import scipy.sparse.linalg as spla
        z = np.asarray(z0, dtype=float).copy()
        n = len(z)
        ier, mesg, xdiff = 2, 'No convergence', None
        ## the "stopped moving and still failing its step test" signature, as
        ## `analysis.fsolve(floor_detect=True)` records it: COUNTED, never
        ## acted on -- the iteration below is what it was
        _floor_q, step_floor = [], None
        _cached = None
        for _i in range(maxiter):
            _b = build(z) if _cached is None else _cached
            _cached = None
            F, mv = _b[0], _b[1]
            _d = None if len(_b) < 3 or _b[2] is None else np.asarray(_b[2], dtype=float)

            def _mv(v, _f=mv, _d=_d):
                return _f(v) if _d is None else _f(_d * np.asarray(v))

            J = spla.LinearOperator((n, n), matvec=_mv, dtype=float)
            xdiff, info = spla.gmres(
                J, -F, rtol=self.KRYLOV_TOLERANCE_FACTOR * reltol,
                restart=min(n, self.KRYLOV_RESTART),
                maxiter=self.KRYLOV_MAX_CYCLES)
            ## ⚠ THE INNER SOLVE'S VERDICT IS NOT DISCARDED: an unconverged
            ## GMRES makes `xdiff` a direction the Newton has no reason to
            ## trust, so the outer loop stops, and the warning names the cause
            ## rather than a generic 'No convergence'.
            ## ⚠ BUT THE VERDICT IS THE RESIDUAL IT REACHED, NOT THE FLAG: the
            ## Krylov target (`KRYLOV_TOLERANCE_FACTOR * reltol`, 1e-11 at
            ## reltol 1e-9) can sit below what the operator's conditioning
            ## allows, and a Newton step needs far less (inexact Newton).
            ## Measured on the comparator oscillator's staged glm2 system:
            ## the scaled operator's condition 1.6e6, GMRES stalled at 4.4e-11
            ## -- flagged a failure, and an excellent step.  A step whose true
            ## residual is below `KRYLOV_ACCEPT_RESIDUAL` is taken.
            if info != 0:
                _res = float(np.linalg.norm(_mv(xdiff) + F))
                _nF = float(np.linalg.norm(F))
                if _nF > 0.0 and _res <= self.KRYLOV_ACCEPT_RESIDUAL * _nF:
                    info = 0
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
            if _d is not None:
                xdiff = _d * xdiff
            z_new = z + xdiff
            if line_search:
                F0 = float(np.linalg.norm(F))

                def _trial(zt):
                    try:
                        bt = build(zt)
                    except (ValueError, np.linalg.LinAlgError,
                            analysis.NoConvergenceError):
                        return None, np.inf
                    return bt, float(np.linalg.norm(bt[0]))
                step = 1.0
                bt, nt = _trial(z_new)
                for _k in range(self.LINE_SEARCH_MAX_HALVINGS):
                    if nt < F0:
                        break
                    step *= 0.5
                    z_new = z + step * xdiff
                    bt, nt = _trial(z_new)
                if bt is None:
                    ## not one trial could be evaluated: stop at the last
                    ## iterate that could, rather than raise from it
                    return (z, {'step_floor': step_floor}, 2,
                            'No convergence (no trial step could be evaluated)')
                xdiff = z_new - z
                _cached = bt

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
