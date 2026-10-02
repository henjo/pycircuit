"""The sequential stage methods (DIRK, ESDIRK, TR-BDF2) and what every stage
step shares.  A theme of `Transient` (see `transient.py`).
"""

import numpy as np

from pycircuit.circuit import _evalhint
from pycircuit.circuit.analysis import (
    remove_row_col,
)
from pycircuit.circuit.dcanalysis import refnode_removed


class _SequentialStages:
    """The sequential stage methods (DIRK, ESDIRK, TR-BDF2) and what every
    stage step shares.  A theme of `Transient` (see `transient.py`)."""

    def _solve_timestep_rk(self, x0, t, provided_function=None):
        """One step of ANY Runge-Kutta method, driven by its Butcher tableau --
        the single entry point that replaces the per-method stage steps.

        Structure-aware (``integ.stage_structure()``):

        * a lower-triangular tableau (DIRK / SDIRK / ESDIRK -- TR-BDF2) is solved
          stage by stage, each an ``m x m`` Newton through :meth:`_newton` (so
          limiting, the nrsolver and the continuation rescue come for free) --
          :meth:`_rk_step_dirk`;
        * a fully-implicit tableau (FULL -- Radau) is solved coupled, dense
          ``3m`` or the ``eig(A^-1)`` cost transform when opted in --
          :meth:`_rk_step_coupled`.

        Stiffly-accurate only for now (``x_{n+1}`` == last stage); both shipped
        stage methods are.  Every path leaves the stage record in ``_rk_Y`` /
        ``_rk_K`` for the shooting monodromy and, when ``_rk_want_est`` is set,
        the filtered embedded estimate in ``_rk_est``.
        """
        integ = self.base_integrator
        if not integ.is_stiffly_accurate():
            raise NotImplementedError(
                'RK step: only stiffly-accurate tableaux are wired (x_{n+1} == '
                'last stage); a non-stiffly-accurate method needs the extra '
                'final-weight solve, not built.')
        if integ.stage_structure() == integ.FULL:
            return self._rk_step_coupled(x0, t, provided_function)
        return self._rk_step_dirk(x0, t, provided_function)

    ## -- what every stage step shares: the source closure, the implicit stage
    ## solve of the sequential methods, and the epilogue ---------------------
    ## History: `doc/transient_history.md`, `Transient._stage_source`.

    def _stage_source(self, provided_function):
        """`u(t)` for a stage step: the circuit's sources at `t`, plus the
        caller's `provided_function`."""
        epar, ana, tk = self.epar, self.par.analysis, self.toolkit

        def src(tt):
            u = tk.array(self.cir.u(tt, epar, analysis=ana), dtype=float)
            if provided_function is not None:
                u = u + provided_function(tt)
            return u
        return src

    def _solve_implicit_stage(self, target, aii, h, ti, guess, src,
                              provided_function, what):
        """One implicit stage of a sequential stage method (DIRK, ESDIRK, a
        Nordsieck GLM): ``q(Y) - target - h a_ii K(Y) = 0`` with ``K =
        -(i(Y) + u(t_i))``, an ``m x m`` Newton via :meth:`_newton`
        (limiting and the continuation rescue included).

        With `pcnr` the stage is first the DC-flow solve `_rk_stage_pcnr`
        takes (divided by ``h a_ii`` it is ``i(Y) + iq_eff + u = 0``) -- the
        same augmented junction continuation the LMM step uses -- and ⚠ IT
        FALLS BACK PER STAGE, as the LMM and coupled steps do: PCNR has no
        continuation ladder, `_newton` does, so a stage PCNR cannot solve is
        handed to the limiting solve rather than ending the transient.

        History: `doc/transient_history.md`, `Transient._solve_implicit_stage`."""
        epar = self.epar
        arr = lambda v: self.toolkit.array(v, dtype=float)

        def func_i(x):
            ## (one evaluation session: `_evalhint`)
            with _evalhint.evaluating('i', 'q', 'C', 'G'):
                Ki = -(arr(self.cir.i(x, epar)) + src(ti))
                f = arr(self.cir.q(x, epar)) - target - h * aii * Ki
                J = (arr(self.cir.C(x, epar))
                     + h * aii * arr(self.cir.G(x, epar)))
            return f, J
        if self._rk_use_pcnr():
            Y, ok = self._pcnr_attempt(
                lambda: self._rk_stage_pcnr(target, aii, h, ti, guess,
                                            provided_function),
                lambda exc: self._pcnr_failed(f'{what} PCNR (a stage)',
                                              ti, exc))
            if ok:
                ## THE BRANCH CHECK, CONFIRMED on the stage equation `_newton`
                ## would have solved
                ## History: `doc/transient_history.md`, `Transient._solve_implicit_stage`.
                if self._branch_on():
                    iref = self.irefnode
                    self._branch_after_solve(
                        refnode_removed(func_i, iref, self.toolkit),
                        self.toolkit.concatenate((Y[:iref], Y[iref + 1:])))
                return Y

        return self._newton(func_i, guess)

    def _finish_stage_step(self, t, tstage, Y, a_last, h, src, K=None):
        """What a stage step leaves for the machinery that reads it, from its
        stages `Y` (stiffly accurate: ``x_{n+1} = Y_{s-1}``): the charge
        cache, `_iq` (the charge derivative at the new point), the per-step
        `C` and companion conductance ``a_ss h G`` shooting's monodromy
        reads, the stages (`_rk_Y`, `_rk_K`) and the predictor nodes for the
        NEXT step (promoted only on accept).  `screen` names the step for the
        branch screen of the paths that confirm no branch themselves.
        Returns ``(x_{n+1}, J)``, ``J = C + a_ss h G``."""
        epar = self.epar
        arr = lambda v: self.toolkit.array(v, dtype=float)
        xnp1 = Y[-1]
        ## what the device memo holds at this very state (the stage Newton's
        ## own, the branch screen's `C`) is read; the rest is evaluated in
        ## one evaluation session (`_evalhint`) and recorded where the memo
        ## rolls -- the next step opens here, and the shooting's stage step
        ## reads `C` and `G` here again
        rec = self._memo_get(xnp1) or {}
        need = [k for k in ('q', 'i', 'C', 'G') if k not in rec]
        if need:
            vals = {}
            with _evalhint.evaluating(*need):
                if 'q' in need:
                    vals['q'] = self.cir.q(xnp1, epar)
                if 'i' in need:
                    vals['i'] = arr(self.cir.i(xnp1, epar))
                if 'C' in need:
                    vals['C'] = arr(self.cir.C(xnp1, epar))
                if 'G' in need:
                    vals['G'] = arr(self.cir.G(xnp1, epar))
            if self._memo_ok():
                self._memo_put(xnp1, vals)
            rec = dict(rec, **vals)
        qY, i_n, Cm, Gm = rec['q'], rec['i'], rec['C'], rec['G']
        self._q_cache = (xnp1, qY)
        self._iq = -(i_n + src(t))
        self._Cmat = Cm
        self._Geq = a_last * h * Gm
        self._effective_method = type(self.base_integrator).__name__
        self._companion_coeffs = None
        self._rk_Y = list(Y)
        self._pred_pending = (t, list(zip(tstage, list(Y))))
        if K is not None:
            self._rk_K = list(K)
        return xnp1, Cm + a_last * h * Gm

    def _rk_step_dirk(self, x0, t, provided_function=None):
        """One step of a lower-triangular (DIRK/SDIRK/ESDIRK) RK method, solved
        stage by stage from the Butcher tableau.

        For stage ``i`` the residual in the charge formulation is
        ``q(Y_i) - [q(x_n) + h sum_{j<i} A_ij K_j] - h A_ii K_i(Y_i) = 0`` with
        ``K_j = -(i(Y_j) + u(t_n + c_j h))`` -- one ``m x m`` Newton via
        :meth:`_newton` (limiting included).  An explicit first stage
        (``A[0]==0``, ``c0==0``) is just ``Y_0 = x_n``.  Stiffly accurate, so
        ``x_{n+1} = Y_{s-1}``.

        For TR-BDF2's ESDIRK tableau the two implicit stages share the
        diagonal ``d = STAGE_DIAG`` (the one-LU property), and the stage form
        is algebraically identical to the TR + BDF2-companion writing.  Any
        SDIRK/ESDIRK to come reuses this untouched.

        History: `doc/transient_history.md`, `Transient._rk_step_dirk`.
        """
        integ = self.base_integrator
        A, B, C = integ.butcher()
        s = A.shape[0]
        h = self._dt
        tn = t - h
        epar = self.epar
        tk = self.toolkit
        arr = lambda v: tk.array(v, dtype=float)
        src = self._stage_source(provided_function)

        xn = x0
        qn = arr(self.cir.q(xn, epar))
        tstage = [tn + C[i] * h for i in range(s)]
        Y = [None] * s
        K = [None] * s
        for i in range(s):
            target = qn + h * sum(A[i, j] * K[j] for j in range(i)) \
                if i > 0 else qn
            aii = A[i, i]
            if abs(aii) < 1e-14:
                ## explicit stage; for the usual explicit FIRST stage (c0==0)
                ## this is x_n exactly, else solve q(Y_i) = target.
                if i == 0:
                    Y[i] = np.array(xn, dtype=float)
                else:
                    def func_e(x, _tgt=target):
                        f = arr(self.cir.q(x, epar)) - _tgt
                        J = arr(self.cir.C(x, epar))
                        return f, J
                    Y[i] = self._newton(func_e, xn)
            else:
                ti = tstage[i]
                ## STAGE PREDICTOR.  The already converged stages of THIS step
                ## are the nearest nodes there are, so the polynomial is
                ## continued by one stage gap rather than a whole step -- the
                ## property that decides its worst case, not its mean.
                guess = self._predict_state(
                    ti, extra=[(tstage[j], Y[j]) for j in range(i)
                               if Y[j] is not None])
                if guess is None:
                    guess = Y[i - 1] if i > 0 else xn
                Y[i] = self._solve_implicit_stage(target, aii, h, ti, guess,
                                                  src, provided_function,
                                                  'stage')
            K[i] = -(arr(self.cir.i(Y[i], epar)) + src(tstage[i]))

        xnp1, J = self._finish_stage_step(t, tstage, Y, A[s - 1, s - 1], h,
                                          src, K=K)
        if getattr(self, '_rk_want_est', False):
            self._rk_est = self._rk_dirk_estimate(K, h, J)
        return xnp1, None, J, None

    def _rk_dirk_estimate(self, K, h, J):
        """The filtered embedded error estimate for a DIRK step, in STATE units.

        ``est_raw = h sum_i dk_i K_i`` (the method's embedded weights
        ``integ.EMBEDDED_DK`` on the stage derivatives -- for TR-BDF2 the
        Hosea & Shampine 2(3) coefficients), then filtered through the last
        stage operator ``(C + a_last h G)^{-1}`` -- the same factor the step's
        final stage used, so the estimate is one back-substitution and stays
        bounded on a stiff mode (an unfiltered ``est_raw`` grows like ``|a h|``
        while the true error is L-damped to zero)."""
        integ = self.base_integrator
        dk = integ.EMBEDDED_DK
        iref = self.irefnode
        tk = self.toolkit
        est_raw = h * sum(dk[i] * np.asarray(K[i]) for i in range(len(dk)))
        (Jr,) = remove_row_col((J,), iref, tk)
        er_r = tk.concatenate((est_raw[:iref], est_raw[iref + 1:]))
        Est_r = tk.linearsolver(Jr, er_r)
        return tk.insert(Est_r, iref, 0.0)
