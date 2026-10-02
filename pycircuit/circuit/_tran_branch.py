"""Branch detection: did the step equation have more than one root.  A theme of
`Transient` (see `transient.py`).
"""

import contextlib

import numpy as np

from pycircuit.circuit import _evalhint
from pycircuit.circuit._limiting import state_restore, state_snapshot
from pycircuit.circuit.simwarnings import (
    ModelWarning,
    warn,
)


class _BranchCheck:
    """Branch detection: did the step equation have more than one root.  A
    theme of `Transient` (see `transient.py`)."""

    ## ------------------------------------------------------------------
    ## BRANCH DETECTION: did the step equation have more than one root?
    ##
    ## At a point where `rank C(x)` DROPS and the surviving dynamics REPEL,
    ## the DAE's solution is genuinely non-unique (Lamour, März & Tischendorf
    ## Thm 3.53) -- and numerically that becomes MULTIPLICITY OF ROOTS OF THE
    ## STEP EQUATION, which the Newton resolves by taking whichever one its
    ## seed is nearest, silently.  MEASURED (`benchmarks/branch_selection.py`,
    ## 2026-09-10): one netlist, one grid, one tolerance returns -1.414744, 0
    ## or +1.414744 depending on the seed alone, with a spread that does NOT
    ## shrink under refinement.  Both ingredients are needed -- with a passive
    ## conductance the equilibrium attracts and uniqueness is safe however
    ## badly `C` degenerates.  A negative conductance is not exotic: it is what
    ## an oscillator's active device is.
    ##
    ## ⚠⚠ THE OBVIOUS SCREEN DOES NOT WORK HERE: firing on the STEP MATRIX
    ## ACQUIRING A NEGATIVE EIGENVALUE fires on EVERYTHING on real MNA.  Two
    ## independent reasons, both structural: MNA WITH A VOLTAGE SOURCE IS A
    ## SADDLE-POINT SYSTEM and is indefinite by construction, and an
    ## OSCILLATOR'S `G` HAS A NEGATIVE EIGENVALUE BY DESIGN with no rank drop
    ## anywhere.  4 false fires out of 4 ordinary circuits, so it cannot gate
    ## anything.
    ##
    ## What is screened instead is the condition itself: `rank C(x)` below the
    ## STRUCTURAL rank -- what `C` has at a generic operating point.  A
    ## structurally zero row (a resistive node, a source branch) is not a drop;
    ## a capacitance that VANISHES is.  That is exactly "im D(t) is not
    ## time-invariant", and it is quiet on all four circuits above.
    ##
    ## ⚠ THE ASYMMETRY IS THE POINT: the confirmation perturbs the seed by a
    ## HEURISTIC amount, so it can MISS a second root -- but when it fires it
    ## has an actual second solution in hand and converged to it.  A warning is
    ## evidence; silence is not.
    ## History: `doc/transient_history.md`, `Transient.BRANCH_SCREEN_TOL`.
    BRANCH_SCREEN_TOL = 1e-9

    def _branch_structural_rank(self):
        """`rank C` at a GENERIC operating point, computed once per run.

        The comparison has to be against this and not against a running
        maximum: on the degenerate branch `C` can be identically zero for the
        whole run, so its rank never "drops" -- it was never up.  Measured:
        screening against a running maximum reads QUIET on the very fixture
        that motivated this.
        """
        cached = getattr(self, '_branch_rank0', None)
        if cached is not None:
            return cached
        rng = np.random.RandomState(20260910)
        best = 0
        scale = 0.0
        for _ in range(3):
            xr = rng.uniform(-1.0, 1.0, self.cir.n)
            try:
                C = np.asarray(self.cir.C(xr, self.epar), dtype=float)
            except Exception as exc:                           # noqa: BLE001
                ## (the check switches itself off for this analysis: said
                ## once, not silently -- the review's F11, 2026-10-01)
                warn(
                    'transient: the branch check could not read the '
                    f'structural rank of C ({type(exc).__name__}: '
                    f'{str(exc)[:80]}) and is OFF for this analysis.', ModelWarning)
                self._branch_rank0 = (0, 0.0)
                return self._branch_rank0
            if C.size == 0:
                continue
            nrm = float(np.max(np.abs(C)))
            scale = max(scale, nrm)
            if nrm > 0:
                sv = np.linalg.svd(C, compute_uv=False)
                best = max(best, int((sv > self.BRANCH_SCREEN_TOL * nrm).sum()))
        self._branch_rank0 = (best, scale)
        return self._branch_rank0

    def _branch_screen(self, x):
        """`(fired, null_direction)` -- has `rank C(x)` fallen below the
        structural rank?

        ⚠ The cheap proxy runs first and is what keeps this affordable: an SVD
        every step would be the same order as the factorisation.  `C`'s
        sparsity pattern is fixed by topology, so a vanishing reactance shows
        up as a structurally nonzero DIAGONAL entry collapsing, which is
        `O(m)` to check.  The SVD only runs when that fires.
        """
        r0, scale0 = self._branch_structural_rank()
        if r0 <= 0 or scale0 <= 0.0:
            return False, None
        ## (inside the converged point's evaluation session, where the step
        ## sets one: its `jacobian_only` reads `q` and `G` here next --
        ## `Transient.solve_timestep`, `_evalhint`)
        conv = getattr(self, '_conv_session', None)
        with (_evalhint.evaluating(session=conv) if conv is not None
              else contextlib.nullcontext()):
            _C = self._C_at_state(x)
        ## (kept for the step's own assembly at this state, which follows)
        self._C_cache = (x, _C)
        C = np.asarray(_C, dtype=float)
        if C.size == 0:
            return False, None
        ## ⚠ THE SCALE IS THE LARGEST `C` THIS RUN HAS SEEN, not the random
        ## states' (`scale0`).  Those sit anywhere in +-1 V, where a forward
        ## junction's diffusion capacitance is e^38 above any state the
        ## circuit visits, so a real junction reads as ZERO against it and
        ## the screen fires (and the confirmation re-solves) on most steps.
        ## The RANK stays structural -- a running maximum of the rank reads
        ## quiet where `C` is zero for the whole run -- and for a linear `C`
        ## the two scales are one number.
        ## History: `doc/transient_history.md`, `Transient._branch_screen`.
        ref = max(getattr(self, '_branch_cmax', 0.0),
                  float(np.max(np.abs(C))))
        self._branch_cmax = ref
        d = np.abs(np.diag(C))
        if float(np.max(np.abs(C))) > self.BRANCH_SCREEN_TOL * ref \
                and float(np.max(d)) > self.BRANCH_SCREEN_TOL * ref:
            ## nothing has collapsed at the cheap level
            if float(np.min(d[d > 0.0]) if np.any(d > 0.0) else 0.0) \
                    > 1e-6 * ref:
                return False, None
        nrm = float(np.max(np.abs(C)))
        if nrm <= self.BRANCH_SCREEN_TOL * ref:
            return True, None            # C has collapsed entirely
        U, sv, _Vt = np.linalg.svd(C)
        r = int((sv > self.BRANCH_SCREEN_TOL * ref).sum())
        if r >= r0:
            return False, None
        return True, U[:, -1]

    def _branch_confirm(self, func, x_res, direction):
        """Re-solve the SAME step from a perturbed seed; return the other root
        if the Newton lands somewhere materially different.

        ⚠ The perturbation magnitude is a HEURISTIC.  The other roots of the
        step equation sit an `O(sqrt(h))` distance away -- measured on the
        reference fixture, they are at 0.061 at 800 points per period and a
        seed of 0.1 reached them while 0.01 did not -- and that distance is not
        knowable in general, so several magnitudes are tried.  A miss is
        therefore possible and silence proves nothing; a fire has an actual
        second solution in hand.
        """
        xr = np.asarray(x_res, dtype=float)
        n = len(xr)
        if direction is None or len(direction) != n:
            direction = np.ones(n) / np.sqrt(n)
        scale = max(float(np.max(np.abs(xr))), 1.0)
        tol = self.par.reltol * scale
        for mag in (1.0, 0.1):
            for sign in (1.0, -1.0):
                seed = xr + sign * mag * scale * np.asarray(direction,
                                                            dtype=float)
                ## ⚠ WITH THE STEP'S OWN LIMITER.  Without it no
                ## `cir.limit` runs, a stateful device (`Diode`) stays
                ## linearised at the ACCEPTED point, and the solve and the
                ## root test below both see that linearisation: measured, a
                ## diode on the rank-dropping node at 0.85 V gave an
                ## "alternative" with residual 0 there and 489 A in the
                ## step equation itself.
                try:
                    alt, _it = self._get_nrsolver().solve_system(
                        seed, func, self.toolkit, self.par.reltol,
                        self._newton_abstol_vector_reduced(),
                        self._newton_xtol_vector_reduced(), self.par.maxiter,
                        limiter=self._newton_limiter())
                except Exception:                              # noqa: BLE001
                    continue
                alt = np.asarray(alt, dtype=float)
                if float(np.max(np.abs(alt - xr))) <= 1e3 * tol:
                    continue
                ## ⚠⚠ AND IT MUST ACTUALLY BE A ROOT.  A solver that hands the
                ## SEED BACK -- converged at iteration zero, or bailed -- looks
                ## exactly like a second solution to a pure distance test, and
                ## that is a FALSE ALARM on a default-on diagnostic.
                ## History: `doc/transient_history.md`, `Transient._branch_confirm`.
                if not self._branch_is_root(func, alt):
                    continue
                return alt
        return None

    def _branch_is_root(self, func, x):
        """Does `x` satisfy the step residual?  The alternative has to be a
        SOLUTION, not merely a different vector."""
        try:
            F, J = func(x)
        except Exception:                                      # noqa: BLE001
            return False
        F = np.asarray(F, dtype=float)
        J = np.asarray(J, dtype=float)
        scale = float(np.max(np.abs(J) @ np.abs(np.asarray(x, dtype=float)))) \
            if J.size else 0.0
        return bool(np.all(np.abs(F)
                           <= self.par.reltol * max(scale, 1.0)
                           + float(np.max(self._newton_abstol_vector()))))

    def _branch_after_solve(self, func, x_res):
        """Screen the converged point, and confirm only if the screen fires.

        ⚠ `x_res` is the REDUCED vector the solver returned; the circuit's
        `C` wants the full one, so the reference row goes back in before the
        screen and comes out again for the re-solve.
        """
        ## ⚠ THE FINDING IS WARNED AFTER THE `try`, NOT IN IT: raised inside,
        ## a warning made an error (`-W error`, the suite's own policy) was
        ## swallowed and reported as the check failing (the review's W3)
        found = None
        try:
            xf = self.toolkit.concatenate(
                (x_res[:self.irefnode], self.toolkit.array([0.0]),
                 x_res[self.irefnode:]))
            fired, direction = self._branch_screen(xf)
            if not fired:
                return
            self._branch_count('branch_screens')
            if direction is not None:
                direction = self.toolkit.concatenate(
                    (direction[:self.irefnode], direction[self.irefnode + 1:]))
            ## ⚠⚠ A DIAGNOSTIC THAT CHANGES THE SIMULATION IS A DEFECT.  Every
            ## speculative solve writes the devices' limiting state (`Diode`'s
            ## `_vlim`, which its `i` and `G` read), so the NEXT step would
            ## start from the alternative's linearisation.  It goes back
            ## EXACTLY, from a snapshot: a `limit(x, x)` re-sync clamps against
            ## the stored state and lands short above a junction's critical
            ## voltage (`state_snapshot`).
            ## History: `doc/transient_history.md`, `Transient._branch_after_solve`.
            _snap = state_snapshot(self.cir)
            try:
                alt = self._branch_confirm(func, x_res, direction)
            finally:
                state_restore(_snap)
            if alt is None:
                return
            self._branch_count('branch_points')
            t = float(getattr(self.epar, 't', 0.0) or 0.0)
            if not getattr(self, '_branch_warned', False):
                self._branch_warned = True
                found = (t, float(np.max(np.abs(
                    alt - np.asarray(x_res, dtype=float)))))
        except Exception as exc:                               # noqa: BLE001
            ## ⚠⚠ A DIAGNOSTIC MUST NOT BE ABLE TO FAIL A SOLVE -- it runs
            ## after the answer is in hand and only reports.  But a bare
            ## `pass` here would hide the check's own bugs (it would silently
            ## do nothing while looking healthy), so the failure is recorded
            ## and announced ONCE.
            ## History: `doc/transient_history.md`, `Transient._branch_after_solve`.
            self._warn_branch_error(exc)
        if found is not None:
            warn(
                f'transient: THE STEP EQUATION HAD MORE THAN ONE ROOT at '
                f't={found[0]:.6g} s. rank C fell below its structural value '
                'there, and re-solving the same step from a different seed '
                'converged to a DIFFERENT solution (largest component '
                f'differs by {found[1]:.3g}). The answer returned is one of '
                'several valid ones and the choice was made by the Newton '
                'seed. See branch_points in the statistics; set '
                'branch_check="off" to skip this test.', ModelWarning)

    def _warn_branch_error(self, exc):
        """The branch check's own failure, warned ONCE a run."""
        if not getattr(self, '_branch_error', None):
            self._branch_error = repr(exc)
            warn(f'transient: the branch check itself failed '
                 f'({self._branch_error}); it is disabled for this run and the '
                 'solve is unaffected', ModelWarning)

    def _branch_after_coupled(self, stage_newton, seed0, Y, residual,
                              build=None):
        """The coupled-path branch check: screen every stage, and if any is at
        a rank drop, re-solve the WHOLE BLOCK from a perturbed seed.

        ⚠ Perturbing one stage is not enough and would be a different question:
        the three stages are unknowns of ONE Newton, so a second root of the
        block is what "the step equation had several roots" means here.
        """
        if not self._branch_on():
            return
        ## (the finding warned after the `try`, as `_branch_after_solve`'s)
        found = None
        try:
            ## ⚠ "FIRED" AND "GAVE A DIRECTION" ARE DIFFERENT ANSWERS: when
            ## `C` collapses ENTIRELY the screen returns `(True, None)` --
            ## there is no null direction to name because every direction is
            ## one -- and a `fired = direction` idiom reads it as "did not
            ## fire".
            ## History: `doc/transient_history.md`, `Transient._branch_after_coupled`.
            fired = False
            direction = None
            for Yi in Y:
                hit, dvec = self._branch_screen(Yi)
                if hit:
                    fired, direction = True, dvec
                    break
            if not fired:
                return
            self._branch_count('branch_screens')
            if stage_newton is None:
                ## built only now the screen has fired: a path that solved
                ## with its own Newton (PCNR, the transform) confirms with
                ## the dense coupled one, on the same step equation
                _ctx, stage_newton, residual, seed0 = build()
            ## ⚠⚠ THE PERTURBATION MUST NOT BE A GAUGE SHIFT.  These are
            ## FULL-WIDTH stage vectors, so a direction of `ones` moves the
            ## REFERENCE NODE too -- a common-mode shift the circuit cannot
            ## see.  The solve leaves the pinned row alone and hands it back
            ## unchanged, the "alternative" differs from the base only in that
            ## row, and its residual is EXACTLY ZERO because it is the same
            ## physical solution -- which the residual check cannot catch.
            ## History: `doc/transient_history.md`, `Transient._branch_after_coupled`.
            n = len(np.asarray(Y[0], dtype=float))
            if direction is not None:
                d = np.asarray(direction, dtype=float).copy()
            else:
                d = np.zeros(n)
                d[(self.irefnode + 1) % n] = 1.0
            if 0 <= self.irefnode < n:
                d[self.irefnode] = 0.0
            dn = float(np.linalg.norm(d))
            if dn <= 0.0:
                return
            d = d / dn
            scale = max(float(np.max([np.max(np.abs(np.asarray(Yi,
                                                               dtype=float)))
                                      for Yi in Y])), 1.0)
            base = np.array([np.asarray(Yi, dtype=float) for Yi in Y])
            _snap = state_snapshot(self.cir)
            try:
                gap = self._branch_coupled_scan(stage_newton, seed0, base, d,
                                                scale, residual)
            finally:
                ## ⚠ ALWAYS, however the scan leaves: every speculative solve
                ## has written the devices' limiting state
                state_restore(_snap)
            if gap is None:
                return
            self._branch_count('branch_points')
            if not getattr(self, '_branch_warned', False):
                self._branch_warned = True
                found = (float(getattr(self.epar, 't', 0.0) or 0.0), gap)
        except Exception as exc:                               # noqa: BLE001
            self._warn_branch_error(exc)
        if found is not None:
            warn(
                f'transient: THE COUPLED STAGE SYSTEM HAD MORE THAN ONE ROOT '
                f'at t={found[0]:.6g} s. rank C fell below its structural '
                'value at one of the stages, and re-solving the block from a '
                'different seed converged to a DIFFERENT solution (largest '
                f'component differs by {found[1]:.3g}). The answer returned is '
                'one of several valid ones. Set branch_check="off" to skip '
                'this test.', ModelWarning)

    def _branch_coupled_scan(self, stage_newton, seed0, base, d, scale,
                             residual):
        """Try the perturbations; return the gap to a genuine second root, or
        `None`.  Split out so the caller can restore limiting state in a
        `finally` however this leaves."""
        for mag in (1.0, 0.1):
            for sign in (1.0, -1.0):
                seed = [np.asarray(seed0[i], dtype=float)
                        + sign * mag * scale * d for i in range(len(seed0))]
                try:
                    alt = stage_newton(seed)
                except Exception:                              # noqa: BLE001
                    continue
                altm = np.array([np.asarray(a, dtype=float) for a in alt])
                gap = float(np.max(np.abs(altm - base)))
                if gap <= 1e3 * self.par.reltol * scale:
                    continue
                ## ⚠⚠ AND IT MUST ACTUALLY BE A ROOT.  `_stage_newton` can
                ## hand the SEED BACK, which a pure distance test reads as a
                ## second solution.  ⚠ A FIXED-POINT test does not catch that:
                ## re-solving from a seed the solve hands back returns it again
                ## and the test passes vacuously.  So the BLOCK RESIDUAL is
                ## assembled and measured -- the alternative has to solve the
                ## equations, not merely be a different vector.
                try:
                    r_alt = residual([a.copy() for a in altm])
                    r_base = residual([b.copy() for b in base])
                except Exception:                              # noqa: BLE001
                    continue
                if r_alt > max(1e3 * r_base, 1e-9 * scale):
                    continue
                return gap
        return None

    def _branch_on(self):
        """`branch_check` is 'on' and has not failed in this run -- a
        failure is announced once and disables the check until the next
        `solve` (`_solve` clears `_branch_error`; PSS, which never calls it,
        keeps a failure for the analysis' life)."""
        return (getattr(self, 'branch_check', 'on') == 'on'
                and not getattr(self, '_branch_error', None))

    def _branch_count(self, name):
        """Count on `statistics` when there is one, on the instance otherwise
        -- a hand-driven march has no `statistics` object."""
        stats = getattr(self, 'statistics', None)
        target = stats if stats is not None else self
        setattr(target, name, getattr(target, name, 0) + 1)
    ## `'on'` (the default) or `'off'`: after each Newton, test whether the
    ## step equation had more than one root -- see `_branch_after_solve`
    branch_check = 'on'
