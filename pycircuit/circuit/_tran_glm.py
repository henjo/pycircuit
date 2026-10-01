"""The Nordsieck general linear methods (glm2, glm3, glm4).  A theme of
`Transient` (see `transient.py`).
"""
from collections import namedtuple

import numpy as np

from pycircuit.circuit.analysis import (
    remove_row_col,
)


GLMStartupTrace = namedtuple('GLMStartupTrace', 'tn hs p xs Ys starter')


class _NordsieckGLM:
    """The Nordsieck general linear methods (glm2, glm3, glm4).  A theme of
    `Transient` (see `transient.py`)."""

    ## ⚠ A NORDSIECK GLM RESTARTS WHERE ITS STEP GROWS BY MORE THAN THIS,
    ## rather than rescaling (`Q_k <- rho^k Q_k`).  The vector leaving a
    ## short step carries the derivatives of a fast stretch (a switch's
    ## transition), and the rescale extrapolates them over a step `rho`
    ## times longer.  Measured on a PWM loop's event-landed grid (jumps up
    ## to 1125x): glm3 1.2e-1 of the swing off a radau reference, glm2
    ## 1.5e-3; restarting, 7.8e-5 and 1.8e-4 (the same at a bound of 1.5).
    ## The adaptive controller grows a step at most 2x, so it never reaches
    ## this.
    GLM_RESTART_GROWTH = 4.0

    def _glm_startup(self, tn, x0, h, provided_function=None):
        """The Nordsieck starting vector ``Q_k = h^k q^(k)(t_n)``, k = 0..p, to
        O(h^p) -- Voigtmann's Theorem 9.5 hypothesis (c), "computed by
        generalised Runge-Kutta methods taking only the initial value".

        p substeps of Radau IIA(3) (order 5 >= p) from ``x0`` at ``h_s = h/p``
        through this Transient's own coupled RK step (so the sources are read
        at the right absolute time), a degree-p interpolant through the p + 1
        charge vectors ``q(x(t_n + k h_s))``, and its scaled derivatives at
        ``t_n``: the interpolant's k-th derivative is O(h_s^{p+1-k}) accurate,
        times ``h^k`` that is O(h^{p+1}) in every component.  Algebraic rows
        (``q == 0``) come out identically zero.  ``_glm_startup_override``, if
        set, replaces this (the gate feeds the EXACT vector through it to
        show the computed one is not what limits the order).
        """
        from math import factorial

        from pycircuit.circuit.integrator import RadauIIA3Integrator
        integ = self.base_integrator
        A, U, B, V, c, p = integ.tableau()
        starter = RadauIIA3Integrator()
        override = getattr(self, '_glm_startup_override', None)
        if override is not None:
            ## (nothing computed, nothing to linearise: no trace, not the
            ## previous startup's -- measured, a stale one was read silently)
            self._glm_startup_trace = None
            return np.asarray(override(tn, x0, h), dtype=float)
        tk = self.toolkit
        epar = self.epar
        hs = h / float(p)
        xs = [np.asarray(x0, dtype=float)]
        Ys = []
        saved = (self.base_integrator, self._dt)
        try:
            self.base_integrator = starter
            self._dt = hs
            for k in range(1, p + 1):
                xk = self._rk_step_coupled(xs[-1], tn + k * hs, provided_function)[0]
                xs.append(np.asarray(xk, dtype=float))
                Ys.append([np.asarray(y, dtype=float) for y in self._rk_Y])
        finally:
            self.base_integrator, self._dt = saved
        ## what the startup did, for a caller that LINEARISES it (the
        ## shooting analysis's GLM period map: `_GLMStartup`) -- the substep
        ## states and each substep's converged stages
        self._glm_startup_trace = GLMStartupTrace(float(tn), float(hs), int(p),
                                                  xs, Ys, starter)
        qs = np.array([np.asarray(self.cir.q(x, epar), dtype=float) for x in xs])   # (p+1, n)
        ## Vandermonde in the scaled variable tau = (t - t_n)/h_s = k: q(tau) = sum_j a_j tau^j
        Vd = np.array([[float(k) ** j for j in range(p + 1)] for k in range(p + 1)])
        a = np.linalg.solve(Vd, qs)                       # (p+1, n): a_j in tau
        ## d^k q / dt^k at t_n = k! a_k / h_s^k ;  Q_k = h^k d^k q/dt^k = (h/h_s)^k k! a_k
        Q = np.array([(h / hs) ** k * factorial(k) * a[k] for k in range(p + 1)])
        Q[0] = qs[0]
        return Q

    def _glm_stage_predictor(self, h, tn, c):
        """The GLM's view of the shared stage predictor: a callable
        ``pred(i, Y)`` giving stage ``i``'s starting guess, or ``None``.

        The nodes are :meth:`_predict_state`'s -- accepted states and stages
        at absolute times -- plus this step's already converged stages, which
        being the nearest ones are what continues the polynomial by ONE STAGE
        GAP instead of a whole step.  ``stage_predictor='off'`` restores the
        old guess (the previous stage's converged value), which is the control
        the gate measures against, not a knob to tune.

        ⚠⚠ THE SIMPLER PREDICTORS LOSE (``Y_i^prev + dx``, and a polynomial
        continued a whole step): a gate on the MEAN passes them, but their
        worst seed and worst stage are several times the old guess's, where
        this one is -12% to -29% device evaluations with the old guess's worst
        case (``benchmarks/stage_predictor.py``).  The mechanism is
        structural: the old guess is always a value the circuit ACTUALLY
        ATTAINED, so it can never sit in a device's overflow region, while a
        polynomial continued a whole step can, and on an exponential a 3x
        overshoot is ``exp(3 dV / VT)``.  ⚠ A ratio test against the step's
        own motion does NOT screen the bad case: the bad prediction's
        displacement and the TRUE stage spread are the same size.

        ⚠ The Nordsieck state itself needs a constant step and an unbroken
        predecessor, so this declines exactly where `_glm_Q` would be rescaled
        or restarted; the shared predictor is happy with a variable step, the
        method is not.

        History: `doc/transient_history.md`, `Transient._glm_stage_predictor`.
        """
        if self.stage_predictor == 'off':
            return None
        prev = getattr(self, '_glm_prev', None)
        if prev is None:
            return None
        Yp, xp, hp, tp = prev
        ## (the step compared RELATIVE to itself: `1e-14 max(h, 1)` read every
        ## step below 1 s as the same to 1e-14 absolute -- a 1e-5 change at h
        ## = 1 ns -- until 2026-10-01, the review's X5)
        if len(Yp) != len(c) or abs(hp - h) > 1e-14 * h \
                or abs(tp - tn) > 1e-12 * max(abs(tn), h):
            return None

        def _pred(i, Y, _c=c, _h=h, _tn=tn):
            return self._predict_state(
                _tn + _c[i] * _h,
                extra=[(_tn + _c[j] * _h, Y[j]) for j in range(i)
                       if Y[j] is not None])
        return _pred

    def _solve_timestep_glm(self, x0, t, provided_function=None):
        """One step of a Nordsieck general linear method (:class:`integrator.
        NordsieckGLMIntegrator`).  In the charge formulation, stage ``i``:

            q(Y_i) = sum_j U_ij Q^[n-1]_j + h sum_{j<=i} A_ij K_j,
            K_j = -(i(Y_j) + u(t_n + c_j h)),

        an ``m x m`` Newton per stage through :meth:`_newton` (A is lower
        triangular with one diagonal, so the Jacobian ``C + h a K G`` is the
        same operator at every stage), then the output

            Q^[n]_k = sum_j V_kj Q^[n-1]_j + h sum_j B_kj K_j,

        with ``x_{n+1} = Y_s`` (stiffly accurate; ``Q^[n]_0 = q(Y_s)`` to
        rounding).  The Nordsieck state lives in ``self._glm_Q`` with the time
        and step it was formed at; the first step (or a step that does not
        continue from the stored time) computes it by :meth:`_glm_startup`,
        a changed step rescales it (``Q_k <- (h/h_old)^k Q_k``).  Algebraic
        rows have ``q == 0`` and every stage satisfies their constraint
        exactly (by induction down the triangular A), so their Nordsieck
        components stay zero.
        """
        integ = self.base_integrator
        A, U, B, V, c, p = integ.tableau()
        s = A.shape[0]
        r = p + 1
        h = self._dt
        tn = t - h
        epar = self.epar
        tk = self.toolkit
        arr = lambda v: tk.array(v, dtype=float)
        src = self._stage_source(provided_function)

        ## ⚠⚠ TWO SLOTS, AND THE SECOND ONE IS WHAT MAKES ADAPTIVE STEPPING
        ## POSSIBLE AT ALL.  `_glm_Q` holds the vector this method last
        ## PRODUCED, valid at the end of that step; `_glm_Q_at_entry` holds the
        ## one it last CONSUMED, valid at the start.  A step that the
        ## controller REJECTS has already overwritten the first, and the retry
        ## -- which starts from the same `x` at the same `tn`, only with a
        ## smaller `h` -- would find no vector valid at `tn` without the second
        ## and run the full startup: a vicious cycle rather than a slow path
        ## (every step rejected once, 42x radau's device evaluations), because
        ## a fresh startup's top component makes the next estimate spurious too.
        ## History: `doc/transient_history.md`, `Transient._solve_timestep_glm`.
        state = getattr(self, '_glm_Q', None)
        entry = getattr(self, '_glm_Q_at_entry', None)
        Q = None
        ## the rescale this step applied, recorded for the shooting
        ## analysis's linearisation (`_GLMStep.forward` scales its sensitivities by
        ## the same `rho^k`)
        self._glm_rho = 1.0
        ## and whether it restarted on growth (`GLM_RESTART_GROWTH`), which
        ## the shooting analysis's linearisation follows
        self._glm_restarted = False
        for cand in (state, entry):
            if cand is not None and abs(cand[1] - tn) <= 1e-12 * max(abs(tn), h):
                Q, _t_old, h_old = cand
                if h > self.GLM_RESTART_GROWTH * h_old:
                    Q = None
                    self._glm_restarted = True
                elif abs(h_old - h) > 1e-14 * h:
                    rho = h / h_old
                    Q = np.array([rho ** k * Q[k] for k in range(r)])
                    self._glm_rho = float(rho)
                break
        started = Q is None
        if started and (state is not None or entry is not None):
            ## a slot existed and did not continue: a fresh start past the
            ## run's opening is a restart, on growth or on a time-key miss
            ## alike (the shooting reads the flag, `_glm_period_blocks`)
            self._glm_restarted = True
        if started:
            Q = self._glm_startup(tn, x0, h, provided_function)
            self.statistics_glm_startups = getattr(self, 'statistics_glm_startups', 0) + 1
        xn = np.asarray(x0, dtype=float)
        ## the vector this step ENTERED with -- the shooting traversal reads it
        ## after the first step to learn what the startup produced
        self._glm_Q_in = np.asarray(Q, dtype=float)
        self._glm_Q_at_entry = (np.array(Q, dtype=float), tn, h)
        tstage = [tn + c[i] * h for i in range(s)]
        pred = self._glm_stage_predictor(h, tn, c)
        Y = [None] * s
        K = [None] * s
        for i in range(s):
            target = sum(U[i, j] * Q[j] for j in range(r)) \
                + h * sum(A[i, j] * K[j] for j in range(i))
            aii = A[i, i]
            ti = tstage[i]
            guess = pred(i, Y) if pred is not None else None
            if guess is None:
                guess = Y[i - 1] if i > 0 else xn
            ## A GLM stage IS the DC-flow form `_rk_stage_pcnr` solves, and
            ## every shipped tableau is DIRK-like with ONE nonzero diagonal
            ## (0.25 / 0.25 / 0.258 for GLM2/3/4), so no explicit stage to
            ## except -- see `_solve_implicit_stage`
            Y[i] = self._solve_implicit_stage(target, aii, h, ti, guess, src,
                                              provided_function, 'GLM stage')
            K[i] = -(arr(self.cir.i(Y[i], epar)) + src(ti))
        Qn = np.array([sum(V[k, j] * Q[j] for j in range(r))
                       + h * sum(B[k, j] * K[j] for j in range(s)) for k in range(r)])
        self._glm_Q = (Qn, t, h)
        ## what the NEXT step's stage predictor extrapolates from
        self._glm_prev = ([np.asarray(y, dtype=float) for y in Y],
                          np.asarray(xn, dtype=float), float(h), float(t))
        xnp1, J = self._finish_stage_step(t, tstage, Y, A[s - 1, s - 1], h,
                                          src, K=K)
        if getattr(self, '_rk_want_est', False):
            if started:
                ## ⚠ NO ESTIMATE ACROSS A RESTART.  The entering vector came
                ## from the startup's interpolant, not from a step of this
                ## method, so `Q_p` is not the same quantity `Qn_p` is and
                ## their difference is not a local error -- it reads large and
                ## the controller rejects a step that was fine, which then
                ## restarts again.  An LMM's first step has no
                ## divided-difference LTE for the same reason and is treated
                ## the same way.
                self._rk_est = np.zeros_like(np.asarray(xnp1, dtype=float))
            else:
                self._rk_est = self._glm_error_estimate(Q[p], Qn[p], J)
        return xnp1, None, J, None

    def _glm_error_estimate(self, Qp_in, Qp_out, J):
        """The local error of one Nordsieck GLM step, in STATE units.

        The top Nordsieck component is ``Q_p = h^p q^(p)`` (this file's
        convention carries no ``1/k!``), so the change across a step is

            Q_p^[n] - Q_p^[n-1] = h^p (q^(p)(t_n) - q^(p)(t_{n-1}))
                                ~ h^(p+1) q^(p+1),

        which is the order the local error of a ``p``-th order method has.
        Both vectors are at THIS step's scale: the entering one was rescaled by
        ``rho^k`` when the step changed, before the step used it, so no second
        rescaling belongs here.

        ⚠ FILTERED through the step's own last-stage operator, exactly as
        :meth:`_rk_dirk_estimate` filters a DIRK's embedded estimate and for
        the same reason: the raw difference is in CHARGE units and an
        unfiltered charge residual GROWS like ``|a h|`` on a stiff mode where
        the true error is L-damped to zero.  One back-substitution, no new
        factorisation.
        """
        iref = self.irefnode
        tk = self.toolkit
        raw = np.asarray(Qp_out, dtype=float) - np.asarray(Qp_in, dtype=float)
        (Jr,) = remove_row_col((J,), iref, tk)
        raw_r = tk.concatenate((raw[:iref], raw[iref + 1:]))
        Est_r = tk.linearsolver(Jr, raw_r)
        return tk.insert(Est_r, iref, 0.0)
