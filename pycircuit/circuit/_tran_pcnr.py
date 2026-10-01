"""Predictor/corrector Newton-Raphson (PCNR) steps and their counts.  A theme
of `Transient` (see `transient.py`).
"""

import numpy as np

from pycircuit.circuit import pcnr as _pcnr
from pycircuit.circuit._limiting import limit_sync
from pycircuit.circuit.analysis import (
    NoConvergenceError,
)
from pycircuit.circuit.dcanalysis import refnode_removed
from pycircuit.circuit.simwarnings import (
    ConvergenceWarning,
    warn,
)


class _PCNRSteps:
    """Predictor/corrector Newton-Raphson (PCNR) steps and their counts.  A
    theme of `Transient` (see `transient.py`)."""

    def _pcnr_augmented(self, x, v_lim, t, junctions, provided_function=None):
        """The PCNR augmented system at `(x, v_lim)` with the step's companion
        and sources folded in as `u_extra` / `J_extra` -- ``(g_mna, g_lim,
        J_mm, J_ml, J_lm, didv)``, see `pcnr.augmented_system`."""
        iq, Geq = self._companion_at(x)
        u = self._source_at(t, provided_function)
        return _pcnr.augmented_system(
            self.cir, x, v_lim, junctions, self.epar,
            u_extra=np.asarray(iq, dtype=float) + np.asarray(u, dtype=float),
            dense_blocks=False, J_extra=Geq)

    def _residual_and_jacobian_pcnr(self, x, v_lim, t, junctions,
                                    provided_function=None):
        """``(f_eff, J_eff, g_lim)`` -- the PCNR system reduced to MNA size.

        The coupled path writes its Newton out by hand rather than delegating to
        `_newton`, so PCNR cannot attach to it the way `_solve_timestep_pcnr`
        attaches to `solve_timestep`.  It does not need to: the Schur-reduced
        system IS an n-sized system whose Newton step equals `predict`'s
        ``dx_mna``, so handing `fang_timestep` ``(f_eff, J_eff)`` in place of
        ``(f, J)`` makes its existing solve work unchanged.
        """
        g_mna, g_lim, J_mm, J_ml, J_lm, didv = self._pcnr_augmented(
            x, v_lim, t, junctions, provided_function)
        f_eff, J_eff = _pcnr.schur_reduce(g_mna, g_lim, J_mm, J_ml, J_lm,
                                          junctions, didv)
        return (self.toolkit.array(f_eff, dtype=float),
                self.toolkit.array(J_eff, dtype=float), g_lim)

    def _pcnr_attempt(self, solve, on_fail, catch=Exception):
        """Run `solve()` by PCNR and keep the run's PCNR bookkeeping:
        `(result, True)` and a solve counted, or -- `catch` raised --
        `on_fail(exc)` (the path's own warning), a fallback counted and
        `(None, False)`, for the caller to solve by device limiting.
        `pcnr_status` is 'used' / 'partial' / 'fell-back' from the two
        counts.

        ⚠ ANY EXCEPTION FALLS BACK, on every path (Andreas), as in DC: a
        PCNR failure on one point -- a singular Jacobian included -- must not
        end the run.  Every path's warning names the exception's type, so a
        genuine bug on the PCNR path is still visible in the log.

        History: `doc/transient_history.md`, `Transient._pcnr_attempt`."""
        try:
            out = solve()
        except catch as exc:
            on_fail(exc)
            self.pcnr_fallbacks += 1
            self.pcnr_status = ('partial' if self.pcnr_solves
                                else 'fell-back')
            return None, False
        self.pcnr_solves += 1
        self.pcnr_status = ('used' if not self.pcnr_fallbacks else 'partial')
        return out, True

    def _pcnr_newton(self, x, v_lim, junctions, assemble, msg, t):
        """The PCNR iteration both step families run: the LMM step
        (`_solve_timestep_pcnr`) and an implicit Runge-Kutta stage
        (`_rk_stage_pcnr`) differ only in the augmented system
        `assemble(x, v_lim)` returns -- the multistep companion, or the
        stage's effective one.  Returns the converged `(x, v_lim, feval)`;
        raises `NoConvergenceError` (`msg % (t, maxiter)`) otherwise.

        History: `doc/transient_history.md`, `Transient._pcnr_newton`."""
        irefnode = self.irefnode
        xtol = self._newton_xtol_vector()
        reltol = self.par.reltol
        feval = 0
        for _it in range(int(self.par.maxiter)):
            g_mna, g_lim, J_mm, J_ml, J_lm, didv = assemble(x, v_lim)
            feval += 1

            dx_mna, dx_lim = _pcnr.predict(g_mna, g_lim, J_mm, J_ml, J_lm,
                                           irefnode, junctions=junctions,
                                           didv=didv)
            x_new = x + dx_mna
            x_new[irefnode] = 0.0
            v_new = _pcnr.refine(junctions, v_lim, v_lim + dx_lim, self.epar,
                                 x_old=x)

            ## BOTH residuals, not just the MNA one, as `solve_dc` does.
            ##
            ## Converging on the MNA one alone can return with
            ## `v_lim != e_a - e_b` -- the diode evaluated at a voltage that is
            ## not the node voltage, so the returned vector is not a solution
            ## of the circuit at all.  Everything downstream then inherits it:
            ## the charge history is wrong, and the LTE estimate built from
            ## that history reads low, so the step controller takes large steps
            ## believing they are accurate.
            ## History: `doc/transient_history.md`, `Transient._pcnr_newton`.
            lim_ok = _pcnr.lim_converged(g_lim, v_new, reltol,
                                         self.par.vabstol)
            done = lim_ok and bool(self.toolkit.alltrue(
                abs(dx_mna) < reltol * abs(x_new) + xtol))
            x, v_lim = x_new, v_new
            if done:
                return x, v_lim, feval
        raise NoConvergenceError(msg % (t, self.par.maxiter))

    def _solve_timestep_pcnr(self, x0, t, provided_function=None):
        """One time point by PCNR rather than by limiting -- STAGE 13.

        The transient residual is ``f = i(x) + iq + u(t)`` with ``J = G + Geq``,
        so the companion terms enter the coupled system as the extra blocks
        `augmented_system` takes; everything else is the DC flow unchanged.

        Returns the same 4-tuple `solve_timestep` does, so the caller cannot tell
        which produced it -- the step controller and history roll are downstream
        of both and must stay so.
        """

        junctions = _pcnr.pcnr_devices(self.cir)
        irefnode = self.irefnode
        ## STAGE PREDICTOR -- and it has to be here, not only on the limiting
        ## path, or the two stop agreeing: with the predictor on one path
        ## only, the two converge to values that differ in the last digits,
        ## that moves the LTE estimate, and the step sequences part company
        ## (`test_gate_13_6_pcnr_and_limiting_take_the_same_steps`).  It also
        ## seeds `v_lim`.
        ## History: `doc/transient_history.md`, `Transient._solve_timestep_pcnr`.
        x = self.toolkit.array(self._pred_or(x0, t), dtype=float).copy()
        v_lim = _pcnr.v_lim_init(junctions, x)

        x, v_lim, feval = self._pcnr_newton(
            x, v_lim, junctions,
            lambda x_, v_: self._pcnr_augmented(x_, v_, t, junctions,
                                                 provided_function),
            'PCNR did not converge at t=%g after %d iterations', t)
        ## SYNC the stateful limiters to the solution (`limit_sync`, as the
        ## stage path does): PCNR never moved them, and what reads the
        ## devices after this step -- the branch check's confirming
        ## re-solve, a limiting fallback's first iteration on the next
        ## step -- would otherwise start from wherever the last limiting
        ## solve left them.  (`augmented_system` excludes the PCNR devices
        ## from the ordinary assembly, so nothing PCNR solves reads them.)
        limit_sync(self.cir, x, self.epar)
        iq, Geq = self._companion_at(x)

        ## THE JACOBIAN HANDED TO THE STEP CONTROLLER MUST BE THE ONE
        ## THIS PATH ACTUALLY SOLVED, and `cir.G(x) + Geq` is not it.
        ##
        ## `cir.G(x)` is the ordinary assembly, and PCNR solved with each
        ## junction's current taken at its OWN unknown `v_lim` instead;
        ## the controller computes `lte = J^-1 Eg`, so it must see the
        ## matrix that was solved.
        ##
        ## The right matrix is the one `predict` factorises: the non-PCNR
        ## part plus each probe's `didv` column as a rank-one update.  At
        ## convergence `v_lim == e_a - e_b`, so it is exactly the
        ## Jacobian of the residual with respect to `x` -- and it is
        ## `schur_reduce`'s matrix, taken from there rather than
        ## written out a second time.
        ## History: `doc/transient_history.md`, `Transient._solve_timestep_pcnr`.
        _g2, _gl2, J_mm2, _Jml2, _Jlm2, didv2 = _pcnr.augmented_system(
            self.cir, x, v_lim, junctions, self.epar,
            u_extra=np.asarray(iq, dtype=float),
            dense_blocks=False, J_extra=Geq)
        _f2, J = _pcnr.schur_reduce(_g2, _gl2, J_mm2, junctions=junctions,
                                    didv=didv2)
        J = self.toolkit.array(J, dtype=float)
        f = self.toolkit.array(
            self.cir.i(x, self.epar) + iq
            + self.cir.u(t, self.epar, analysis=self.par.analysis),
            dtype=float)
        ## and its own predictor node, exactly as the limiting
        ## path records one -- symmetry is what gate 13-6 asks for
        self._pred_pending = (t, ())
        ## THE BRANCH CHECK, CONFIRMED: the step equation is the one
        ## the limiting path solves, so the speculative re-solve is
        ## that Newton's (`_branch_after_solve`)
        if self._branch_on():
            self._branch_after_solve(
                refnode_removed(
                    lambda xx: self._residual_and_jacobian(
                        xx, t, provided_function),
                    irefnode, self.toolkit),
                self.toolkit.concatenate((x[:irefnode],
                                          x[irefnode + 1:])))
        return x, feval, J, f

    def _rk_use_pcnr(self):
        """Whether the per-stage solve should use PCNR (the first-class
        junction-continuation limiting) rather than device `limit()` -- true
        when ``pcnr`` is asked for and the circuit has a participating device."""
        if not self.par.pcnr:
            return False
        return bool(_pcnr.pcnr_devices(self.cir))

    def _rk_stage_pcnr(self, target, aii, h, ti, guess, provided_function=None):
        """Solve ONE implicit RK stage by PCNR instead of device limiting.

        The stage residual ``q(Y) - target - h a_ii K(Y) = 0`` (``K = -(i+u)``)
        divided by ``h a_ii`` is the DC-flow form PCNR solves,
        ``i(Y) + iq_eff + u = 0`` with the EFFECTIVE companion

            iq_eff = (q(Y) - target) / (h a_ii),   Geq_eff = C(Y) / (h a_ii)

        so the same augmented junction-continuation the LMM step uses
        (`pcnr.augmented_system`/`predict`/`refine`) applies unchanged, with
        the stage's ``iq_eff``/``Geq_eff`` in place of the multistep companion.
        Returns the converged full-size stage value ``Y``; raises
        `NoConvergenceError` if PCNR does not converge (the caller does not fall
        back -- PCNR is the chosen limiting)."""
        junctions = _pcnr.pcnr_devices(self.cir)
        tk = self.toolkit
        epar = self.epar
        ana = self.par.analysis
        scale = h * aii
        target = np.asarray(target, dtype=float)
        x = tk.array(guess, dtype=float).copy()
        v_lim = _pcnr.v_lim_init(junctions, x)
        def assemble(x_, v_):
            ## the stage's EFFECTIVE companion in place of the multistep one
            q = np.asarray(self.cir.q(x_, epar), dtype=float)
            C = self.cir.C(x_, epar)
            iq_eff = (q - target) / scale
            Geq_eff = np.asarray(C, dtype=float) / scale
            u = np.asarray(self.cir.u(ti, epar, analysis=ana), dtype=float)
            if provided_function is not None:
                u = u + np.asarray(provided_function(ti), dtype=float)
            return _pcnr.augmented_system(
                self.cir, x_, v_, junctions, epar,
                u_extra=iq_eff + u, dense_blocks=False, J_extra=Geq_eff)
        x, _v_lim, _feval = self._pcnr_newton(
            x, v_lim, junctions, assemble,
            'Radau/DIRK stage PCNR did not converge at t=%g after %d iterations',
            ti)
        ## SYNC the devices' internal limiting voltage to the converged
        ## solution.  PCNR never calls `cir.limit`, so each junction's
        ## `_vlim` is left stale -- and the caller's downstream
        ## `i(Y)`/`G(Y)`/`C(Y)` (the stage derivative K, the returned J,
        ## the estimate) linearise there.  At convergence the junction
        ## voltage IS the node voltage, and `limit_sync` puts `_vlim` on it
        ## without altering the solution.  ⚠ Not `limit(x, x)` alone, which
        ## clamps against the stale state: 50 mV short on a hard-driven
        ## diode, and the next stage's K read there put the TR-BDF2
        ## waveform 6.7e-3 V off the exact one.
        limit_sync(self.cir, x, epar)
        return x

    def _reset_pcnr_counts(self):
        """Zero the run's PCNR counters (`pcnr_solves`, `pcnr_fallbacks`,
        `pcnr_status`) and its first recorded failure -- at construction and
        at every run.

        PCNR OUTCOME, per run (roadmap sec. 47).  `pcnr=True` is a
        request: PCNR can decline for the whole run (no device
        declares a probe) or fail on individual timesteps and fall
        through to the ordinary step solver.  DC reports a single
        `pcnr_status`; a transient cannot, because the answer differs
        per step -- so it COUNTS.

        These count SOLVER INVOCATIONS, not accepted steps: a rejected
        step is solved and then thrown away, and it is still a step
        PCNR did or did not carry.

        SETTLED 2026-08-31 by the branch author, asked directly: this is
        the honest number.  The consequence is intended -- `pcnr_solves +
        pcnr_fallbacks` will generally EXCEED `statistics.accepted_steps`,
        and that is not a bug to reconcile.  Do not "fix" these to track
        accepted steps; that would hide work that actually happened.
        """
        self.pcnr_solves = 0
        self.pcnr_fallbacks = 0
        self.pcnr_status = 'off'
        self._pcnr_first = None

    def _pcnr_failed(self, where, t, exc, note=''):
        """A PCNR fallback, kept for the run's one warning
        (`_warn_pcnr_summary`): until 2026-10-01 a `logging.warning` per
        failed step (or stage), outside the warnings machinery altogether
        (the review's X8)."""
        if getattr(self, '_pcnr_first', None) is None:
            self._pcnr_first = (where, float(t), type(exc).__name__,
                                str(exc)[:80], note)

    def _warn_pcnr_summary(self, who):
        """The run's PCNR fallbacks, one warning: how many, and the first."""
        first = getattr(self, '_pcnr_first', None)
        if first is None or not self.pcnr_fallbacks:
            return
        where, t0, ename, emsg, note = first
        total = self.pcnr_fallbacks + self.pcnr_solves
        warn(
            f'{who} pcnr=True: PCNR failed on {self.pcnr_fallbacks} solve(s) '
            f'of {total} -- the first {where} at t={t0:g} ({ename}: {emsg}); '
            'each fell back to device limiting or the ordinary solver '
            f'(pcnr_status {self.pcnr_status!r}){note}', ConvergenceWarning)
        self._pcnr_first = None
