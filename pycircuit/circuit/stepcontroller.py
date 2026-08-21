from abc import ABC, abstractmethod
import warnings

import numpy as np

## Valid values for `relref`, matching Spectre's parameter of the same name.
RELREF_MODES = ('pointlocal', 'alllocal', 'sigglobal')

## STAGE 4b -- the largest factor by which one ACCEPTED step may exceed the one
## before it.  This is a stability bound, not a taste parameter.  Variable-step
## BDF-2 is zero-stable only for `h_n/h_{n-1} < 1 + sqrt(2) = 2.414214`
## (Grigorieff); above it the parasitic root of the homogeneous recursion leaves
## the unit disc and the previous solution is amplified rather than forgotten.
##
## The value below sits *inside* that bound rather than on it: at exactly
## 2.414214 the parasitic root is 1.000000, i.e. marginal, and a method that runs
## permanently at its own stability boundary has no margin for the rounding that
## put it there.  2.0 is what both controllers already used as an unexplained
## literal; naming it is what lets `transient.py`'s force-accept path honour the
## same bound instead of quietly taking 10x -- see `ZERO_STABILITY_RATIO` in
## `integrator.py`, which is the backstop for anything that still gets past this.
MAX_GROWTH_RATIO = 2.0

## The matching floor on shrink, applied on the rejection path.  Unlike the
## growth bound this one is pure economics -- shrinking is unconditionally
## zero-stable for BDF-2 -- and it exists so a single over-tolerance estimate
## cannot collapse the step by orders of magnitude in one move.
MIN_SHRINK_RATIO = 0.2


class StepController(ABC):
    """
    Abstract Strategy Interface for deciding and predicting time steps.
    """

    ## ITEM 2+.3 -- what the RELATIVE part of the LTE tolerance is measured against.
    ##
    ## The tolerance is `lteratio * (reltol*ref + abstol)`.  Until now `ref` was
    ## hard-coded to `max(|x_curr|, |x_last|)` -- each unknown against itself, at
    ## this instant.  That is Spectre's `pointlocal`, and it has a failure mode that
    ## is easy to miss: on a node carrying no signal `ref -> 0`, so the tolerance
    ## collapses to `abstol` and the controller starts chasing numerical noise on a
    ## quiet node.  On the leapfrog that alone cut the step size 5.4x, and the fix
    ## applied at the time was to raise the absolute floor a millionfold
    ## (`lte_vabstol` 1e-12 -> 1e-6), which treats the symptom.
    ##
    ## Spectre's answer is `relref`, and its default is `sigglobal`: measure each
    ## signal against the largest signal anywhere in the circuit, over all past
    ## time, so a quiet node inherits a sane reference instead of degenerating.
    ##
    ##   pointlocal  each unknown against itself, now.  (pycircuit's historical
    ##               behaviour, and still selectable.)
    ##   alllocal    each unknown against its OWN largest value so far.
    ##   sigglobal   each unknown against the largest value of ANY unknown so far.
    ##
    ## THE WORKAROUND ABOVE IS NOW GONE.  With `sigglobal` shipped, `lte_vabstol` is
    ## back to 1e-12: measured at gate D3-e, 1e-6 / 1e-9 / 1e-12 give bit-identical
    ## runs under `sigglobal` (403 steps on a pulsed RC, 601 with a quiet node, at
    ## every value), where under `pointlocal` the same change costs 8.5-9.2%.  That
    ## difference IS the symptom, and it is what the floor was raised to hide.
    ##
    ## DEFAULT IS `sigglobal` SINCE DECISION D3's SECOND ATTEMPT, matching Spectre.
    ## It was adopted, sent back by its own gate, and re-run once the reason for the
    ## failure was removed -- see the D3 gates in `doc/transient_work_plan.md`.
    relref = 'sigglobal'

    def set_relref(self, relref):
        if relref not in RELREF_MODES:
            raise ValueError(
                "relref must be one of %r, not %r" % (RELREF_MODES, relref))
        self.relref = relref
        self._ref_running = None
        return self

    ## STAGE 12A -- Fang's two-sided acceptance band, eq (15), and the step-change
    ## damper, eq (16).  (DAC 2013; see doc/transient_work_plan.md STAGE 12.)
    ##
    ## The classical test is one-sided: accept whatever falls under tolerance.
    ## Fang accepts only inside a BAND, `gamma_min*tau <= eps <= gamma_max*tau`,
    ## and sec. 4.1 attributes the paper's headline result to the LOWER bound:
    ## "the lower bound gamma_min prevents step sizes from being unnecessarily
    ## small."  A step far under tolerance is not free -- it is a step that could
    ## have covered more time for the same work.
    ##
    ## In this module's normalisation `err = eps/tau` (the tolerance already folds
    ## in TRTOL), so the band is simply `gamma_min <= err <= gamma_max`.
    ##
    ## THE DEFAULTS BELOW REPRODUCE THE PREVIOUS BEHAVIOUR EXACTLY: `gamma_min=0`
    ## makes the lower test vacuous and `gamma_max=1` is the historical `err > 1`
    ## rejection.  That is deliberate -- stage 12 is behind a flag until 12D, and
    ## a band that changed the default path would make its own gate unreadable.
    ##
    ## The paper's own values (`ltemin=0.7`, `ltemax=3.0` in sec. 4.1) are NOT
    ## adopted as defaults, and not because they are unattractive: they are quoted
    ## against a comparison method that redid a step above a normalised LTE of
    ## 4.63, so "1.0" there is not "1.0" here.  Copying them would be using an
    ## external number whose normalisation is not established -- so they are
    ## available to a caller and not baked in.
    lte_gamma_min = 0.0
    lte_gamma_max = 1.0
    lte_eta = None

    def set_lte_band(self, gamma_min=0.0, gamma_max=1.0, eta=None):
        """Select Fang's acceptance band and step-change damper.

        ``eta`` is eq (16)'s relative limit on how far one step may move from the
        one before it, ``|dh| <= eta*h`` (the paper suggests ~15%); ``None``
        leaves the step change limited only by the zero-stability bound.

        ``'auto'`` -- the Transient parameters' unset sentinel (F5,
        doc/transient_review_260820.md) -- resolves HERE to this method's own
        defaults, which are the standard path's historical one-sided test.
        The mapping lives inside set_lte_band so the stored band attributes
        are always numeric (eta: float or None) and no controller ever sees
        the sentinel; the coupled path resolves the same sentinel to Fang's
        values in Transient._coupled_band, which bypasses this method.
        """
        if gamma_min == 'auto':
            gamma_min = 0.0
        if gamma_max == 'auto':
            gamma_max = 1.0
        if eta == 'auto':
            eta = None
        if not (0.0 <= gamma_min < gamma_max):
            raise ValueError(
                "LTE band requires 0 <= gamma_min < gamma_max, got %r, %r"
                % (gamma_min, gamma_max))
        if eta is not None and not eta > 0.0:
            raise ValueError("LTE band damper eta must be positive, got %r" % (eta,))
        self.lte_gamma_min = float(gamma_min)
        self.lte_gamma_max = float(gamma_max)
        self.lte_eta = None if eta is None else float(eta)
        return self

    def _band_target(self, safety, p):
        """The normalised error the step prediction aims at.

        Without a band this is ``safety**p`` -- which is not a new choice, it is
        what the existing law ``h*safety*(1/err)**(1/p)`` already converges to,
        and the stage-12 entry measurement found the accepted steps sitting on it
        to within half a percent.

        With a band, the aim only moves if the band actually excludes where the
        controller was already going -- a band containing the natural aim point is
        inert by construction, which is the honest behaviour and must not be
        disguised by silently re-aiming at the band centre.

        WHEN the aim does have to move it goes to the band's geometric centre,
        NOT to the edge that excluded it.  Clipping to the edge was the first
        implementation and gate 12A-1 measured what it costs: aiming exactly at
        `gamma_min` makes every undershoot a rejection, so `gamma_min=0.95` took
        3172 rejections to accept 1187 steps -- more than two redos per step, to
        save 7.8%.  The geometric mean is the point furthest from both edges in
        the ratio sense, which is the sense the step-size law works in.

        Fang does not need this because there `h` is an unknown solved to satisfy
        the band, not a prediction tested against it; a predict-then-test scheme
        must aim strictly inside or it rejects on its own rounding.
        """
        lo, hi = self.lte_gamma_min, self.lte_gamma_max
        target = safety ** p
        if lo > 0.0 or hi != 1.0:
            if not (lo <= target <= hi):
                ## `lo == 0` is an upper bound only; there is no lower edge to be
                ## centred against, so keep the usual safety margin under `hi`.
                target = (lo * hi) ** 0.5 if lo > 0.0 else hi * safety
        return target

    def _damp(self, h_next, h_curr):
        """Eq (16): limit the step change to ``eta`` of the current step."""
        if self.lte_eta is None:
            return h_next
        return min(max(h_next, h_curr * (1.0 - self.lte_eta)),
                   h_curr * (1.0 + self.lte_eta))

    def _reference(self, x_curr, x_last, no_history, n_nodes, toolkit):
        """The `ref` in `reltol*ref + abstol`, per the selected `relref`.

        `n_nodes` splits node voltages from branch currents.  The global modes
        must NOT mix them: a circuit with amperes of branch current and millivolts
        of node signal would otherwise reference every node to a current, which is
        dimensional nonsense and would silently disable node error control.  When
        `n_nodes` is None the vector is treated as one group, which is correct only
        if every entry shares a unit -- callers that know better should say so.
        """
        local = toolkit.maximum(abs(x_curr), abs(x_last))
        if self.relref == 'pointlocal':
            return local

        if no_history or getattr(self, '_ref_running', None) is None:
            ## First step of a run: nothing is remembered yet, so the running
            ## reference starts from what is in front of us.
            self._ref_running = local
        else:
            self._ref_running = toolkit.maximum(self._ref_running, local)

        if self.relref == 'alllocal':
            return self._ref_running

        ## sigglobal: collapse each unit group to its maximum and broadcast back.
        running = self._ref_running
        out = np.array(running, dtype=float, copy=True)
        if n_nodes is None or n_nodes <= 0 or n_nodes >= len(out):
            out[:] = np.max(out) if len(out) else 0.0
            return out
        out[:n_nodes] = np.max(running[:n_nodes])
        out[n_nodes:] = np.max(running[n_nodes:])
        return out

    @abstractmethod
    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, no_history, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0, n_nodes=None, h_last2=None, h_clamped=False, x_hist=None):
        """Evaluate the Local Truncation Error (LTE) for the current step.

        ``no_history`` means there is genuinely no past point to difference
        against -- the first step of a run -- so the LTE cannot be estimated and
        the step has to be accepted unevaluated.  It is deliberately *not* the
        transient's ``_is_first_step``, which is re-armed at every breakpoint to
        force an order drop: history still exists there, and a step that can be
        bounded should be.

        Returns:
            tuple: ``(accept_step, h_next)`` -- whether the step is accepted and
            the predicted next step size.
        """
        pass

class IntegralController(StepController):
    """
    Standard Integral Step Controller based on Yao et al. ICECS 2014.
    Rejects steps with LTE > 1.0, and predicts the next step size.
    """
    
    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, no_history, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0, n_nodes=None, h_last2=None, h_clamped=False, x_hist=None):
        ## Cleared on entry so `last_err` is None wherever no error was computed,
        ## rather than silently holding the previous step's value -- a stale
        ## reading is worse than a missing one for anything measuring the
        ## distribution of accepted-step errors.
        self.last_err = None

        ## No past point exists yet, so there is nothing to difference and the
        ## step is accepted unevaluated.  This is the only place in a run where
        ## that is correct, and it costs one uncontrolled step of O(h^2) Euler
        ## error at max_step -- which is why it used to dominate every accuracy
        ## measurement when breakpoints re-armed it periodically.  It still does
        ## on the stiff ringdown: gate 12A-2 measured the total error there
        ## saturating at 1.3589e-02 from this one step, unchanged across four
        ## decades of `reltol`.
        if no_history:
            return True, h_curr

        # --- LOCAL TRUNCATION ERROR (LTE) CALCULATION ---
        # 1. Ask the active integrator (e.g. Gear2, Trapezoidal) to calculate the 
        #    raw unscaled truncation error vector (Eg) based on charge curvature.
        Eg, p = active_integrator.compute_lte(
            q_curr=q_curr,
            h_curr=h_curr,
            q_last=q_last_hist,
            iq_last=iq_last_hist,
            h_last=h_last,
            is_first_step=no_history,
            toolkit=toolkit,
            ## None until the run has three real past charges; see the ABC.
            h_last2=h_last2,
        )
        
        # 2. Convert charge error into voltage error by multiplying by the 
        #    inverse of the Jacobian matrix: lte = J^-1 * Eg
        from pycircuit.circuit.analysis import remove_row_col
        J_reduced, Eg_reduced = remove_row_col((J, Eg), irefnode, toolkit)
        
        try:
            lte_reduced = toolkit.linearsolver(J_reduced, Eg_reduced)
        except Exception as exc:
            ## DECISION 0.3d called this "the unlogged half-(B) fallback" and said
            ## to delete it.  What made it a defect is not the fallback but the
            ## SILENCE: `Eg` is a current and `J^-1 Eg` is a voltage (see the
            ## units note in the plan), so this substitutes one flavour for the
            ## other and then compares the result against a voltage tolerance.
            ## An error that is wrong by a factor of `h` either never fires or
            ## always does, and either way looks like it is working.
            ##
            ## Kept rather than removed, because removing it would turn a
            ## degenerate Jacobian into an exception from inside step control
            ## rather than from the operating point, where the diagnostic is much
            ## better.  Made loud instead: gate 4-D measured this firing ZERO
            ## times on a circuit with cond(J) = 1.0e12, so if it ever does fire
            ## that is news.
            warnings.warn(
                'transient step control: the LTE solve failed (%s), so the '
                'charge-domain residual is being used in place of the '
                'solution-domain error. These are different quantities -- the '
                'first is a current, the second a voltage -- so the step size '
                'from this point is not error-controlled in the usual sense.'
                % exc, RuntimeWarning, stacklevel=2)
            lte_reduced = Eg_reduced
            
        lte = toolkit.concatenate((lte_reduced[:irefnode], toolkit.array([0.0]), lte_reduced[irefnode:]))

        # 3. Dynamic per-node tolerance, relaxed by the transient error factor TRTOL.
        #    TRTOL is the SPICE "transient tolerance" (Spectre calls the equivalent
        #    `lteratio`): the LTE estimate is deliberately conservative, so the
        #    allowed truncation error is TRTOL times the Newton-solve tolerance.
        #    Folding TRTOL into etol makes the accept threshold (err<=1) and the
        #    step-size prediction aim at the same target instead of oscillating.
        ref = self._reference(x_curr, x_last, no_history, n_nodes, toolkit)
        etol = TRTOL * (reltol * ref + abstol)

        # 4. Normalize the error
        err_array = abs(lte) / etol
        err = float(np.max(err_array))
        ## Exposed under the same name `PIController` uses, so the normalised
        ## error of whichever controller is running can be read from outside.
        ## Not used by this controller's own law -- it is pure integral -- but a
        ## step-size band (Fang's gamma_min, stage 12) is a statement about this
        ## number, and it was otherwise a local that nothing could observe.
        self.last_err = err

        # Step-size prediction exponent is 1/(order+1); compute_lte returns that
        # (order+1) as ``p`` (2 for Euler, 3 for the 2nd-order methods).
        exponent = 1.0 / p
        safety = 0.9

        ## STAGE 12A.  `target` is the normalised error the prediction aims at.
        ## With no band it is `safety**p`, so `(target/err)**(1/p)` is identically
        ## the old `safety*(1/err)**(1/p)` -- the rewrite below changes no
        ## arithmetic on the default path, it only names the aim point so a band
        ## can move it.
        target = self._band_target(safety, p)

        # --- STEP REJECTION / ACCEPTANCE ---
        if err > self.lte_gamma_max:
            # Step rejected: shrink and recalculate.
            ##
            ## THE DAMPER IS DELIBERATELY NOT APPLIED HERE.  Eq (16) bounds how
            ## far one accepted step may sit from the one before it; it is not a
            ## limit on how fast a step that failed its error test may retreat.
            ## Applying it to the rejection path was measured on the stiff RLC
            ## ringdown: with eta=0.15 the step could only shrink 15% per retry,
            ## so it exhausted MAX_REJECT, force-accepted, and crossed the whole
            ## ringing transient in 62 steps against the baseline's 490 -- with a
            ## reported LTE of exactly zero, because by then it was integrating a
            ## signal that had already decayed.  A limiter that makes the error
            ## control unable to respond is not a damper, it is a muzzle.
            h_next = h_curr * max(MIN_SHRINK_RATIO, (target / err) ** exponent)
            return False, h_next

        ## Eq (15)'s LOWER bound: the step landed so far under tolerance that it
        ## was wasted work, so it is redone LARGER at the same time point.  This
        ## is the one place in the controller that rejects a step for being too
        ## ACCURATE, and it is the mechanism sec. 4.1 credits for the paper's 39%.
        ##
        ## `h_clamped` suppresses it, and that suppression is not a detail: a step
        ## truncated onto a breakpoint or onto `tend` is not LTE-limited, so its
        ## small error says nothing about the integrator and growing it is either
        ## impossible or wrong.  The stage-12 entry measurement found exactly this
        ## population -- at loose tolerance the steps sitting far below target were
        ## breakpoint-clamped, not controller-chosen.  Without this guard the band
        ## would spend its retries re-solving steps whose size was never its to
        ## choose.  Likewise a step already at `max_step` has nowhere to grow.
        if (err < self.lte_gamma_min and not h_clamped
                and h_curr < max_step * (1.0 - 1e-12)):
            h_next = h_curr * min(MAX_GROWTH_RATIO,
                                  (target / max(err, 1e-12)) ** exponent)
            h_next = min(self._damp(h_next, h_curr), max_step)
            ## Never report a "grow" that does not actually grow: with the damper
            ## or `max_step` binding, the retry would re-solve the same step and
            ## reject it again, which is a livelock with extra steps.
            if h_next > h_curr * (1.0 + 1e-9):
                return False, h_next

        # Step accepted: predict next size with the same controller law.
        h_next = h_curr * min(MAX_GROWTH_RATIO,
                              (target / max(err, 1e-12)) ** exponent)
        h_next = min(self._damp(h_next, h_curr), max_step)

        return True, h_next

class PIController(StepController):
    """
    Proportional-Integral Step Controller.
    Uses history of truncation error to provide smoother step size changes.
    """
    ## STAGE 4a -- THE GAINS ARE PER UNIT ORDER, AND THAT DIVISION WAS MISSING.
    ##
    ## Gustafsson's classic values are `k_I = 0.3/k` and `k_P = 0.4/k` where `k` is
    ## the order the error estimate follows (`err ~ h^k`).  They were used here
    ## undivided.  Linearising the update about the fixed point with
    ## `x_n = ln(h_n/h*)` and `err_n = (h_n/h*)^p` gives
    ##
    ##     x_{n+1} = x_n (1 - p(k_I + k_P)) + x_{n-1} (p k_P)
    ##     char:    z^2 - (1 - p(k_I + k_P)) z - p k_P = 0
    ##
    ##     gains              p=2 roots            p=3 roots           radius
    ##     0.3, 0.4     -1.1165, +0.7165     -1.7758, +0.6758     1.117 / 1.776
    ##     0.3/p, 0.4/p +0.8000, -0.5000     +0.8000, -0.5000     0.800 / 0.800
    ##
    ## So the loop was unstable at both orders.  It never looked like a divergence
    ## because `min(2, max(0.2, .))` clamps the factor, which converts the growing
    ## oscillation into a PERMANENT period-2 limit cycle: measured h alternating
    ## 0.8572 / 0.4286, i.e. running against the growth clamp every other step, for
    ## as long as the run lasts.  The only test asserted `len(steps) > 10`.
    ##
    ## Note the corrected radius is the SAME at both orders.  That is the point of
    ## dividing by `k` rather than choosing two smaller constants: it makes the
    ## closed loop order-independent, so the controller behaves the same behind
    ## Euler as behind Gear-2.  Gains that were merely smaller would be tuned to
    ## whichever order they were measured on.
    def __init__(self, k_i=0.3, k_p=0.4):
        self.k_i = k_i
        self.k_p = k_p
        self.last_err = None
        
    def pi_factor(self, err, last_err, p):
        """The step-size factor for one accepted step, clamped.

        Public and separate so the gate-4a harness can drive the REAL update law
        instead of transcribing it.  `benchmarks/transient_stage4.py --pi` used to
        keep its own copy, which meant the gate could pass against a formula the
        simulator no longer used -- the same two-transcriptions failure this plan
        has already paid for twice elsewhere.
        """
        err_norm = max(err, 1e-12)
        err_last_norm = max(last_err if last_err is not None else err, 1e-12)
        ## `/p` on both gains -- see the derivation on __init__.
        factor = ((err_norm ** (-self.k_i / p))
                  * ((err_last_norm / err_norm) ** (self.k_p / p)))
        return min(MAX_GROWTH_RATIO, max(MIN_SHRINK_RATIO, factor))

    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, no_history, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0, n_nodes=None, h_last2=None, h_clamped=False, x_hist=None):
        ## As in IntegralController: nothing to difference on the first step of a
        ## run.  Unlike there, the 0.5 is not dead -- it seeds the PI history so
        ## the first real update has a previous error to work from.
        if no_history:
            self.last_err = 0.5
            return True, h_curr
        
        Eg, p = active_integrator.compute_lte(
            q_curr=q_curr,
            h_curr=h_curr,
            q_last=q_last_hist,
            iq_last=iq_last_hist,
            h_last=h_last,
            is_first_step=no_history,
            toolkit=toolkit,
            ## None until the run has three real past charges; see the ABC.
            h_last2=h_last2,
        )
        
        from pycircuit.circuit.analysis import remove_row_col
        J_reduced, Eg_reduced = remove_row_col((J, Eg), irefnode, toolkit)
        
        try:
            lte_reduced = toolkit.linearsolver(J_reduced, Eg_reduced)
        except Exception as exc:
            ## DECISION 0.3d called this "the unlogged half-(B) fallback" and said
            ## to delete it.  What made it a defect is not the fallback but the
            ## SILENCE: `Eg` is a current and `J^-1 Eg` is a voltage (see the
            ## units note in the plan), so this substitutes one flavour for the
            ## other and then compares the result against a voltage tolerance.
            ## An error that is wrong by a factor of `h` either never fires or
            ## always does, and either way looks like it is working.
            ##
            ## Kept rather than removed, because removing it would turn a
            ## degenerate Jacobian into an exception from inside step control
            ## rather than from the operating point, where the diagnostic is much
            ## better.  Made loud instead: gate 4-D measured this firing ZERO
            ## times on a circuit with cond(J) = 1.0e12, so if it ever does fire
            ## that is news.
            warnings.warn(
                'transient step control: the LTE solve failed (%s), so the '
                'charge-domain residual is being used in place of the '
                'solution-domain error. These are different quantities -- the '
                'first is a current, the second a voltage -- so the step size '
                'from this point is not error-controlled in the usual sense.'
                % exc, RuntimeWarning, stacklevel=2)
            lte_reduced = Eg_reduced
            
        lte = toolkit.concatenate((lte_reduced[:irefnode], toolkit.array([0.0]), lte_reduced[irefnode:]))

        # Relax the LTE tolerance by TRTOL (see IntegralController) so the accept
        # threshold matches the target the PI update drives toward.
        ref = self._reference(x_curr, x_last, no_history, n_nodes, toolkit)
        etol = TRTOL * (reltol * ref + abstol)
        err_array = abs(lte) / etol

        err = float(np.max(err_array))
        exponent = 1.0 / p

        if err > 1.0:
            # Step rejected: standard backoff using the method order.
            ## STAGE 4a -- AND THE HISTORY IS DROPPED, DELIBERATELY.
            ##
            ## This used to return without touching `last_err`, so the next accepted
            ## step differenced against an error two steps stale -- and one measured
            ## at a DIFFERENT step size at the same time point, which is not a
            ## sequence the P term is meaningful over.
            ##
            ## Setting it to None makes the next accepted step take the elementary
            ## (pure-I) update, because `if self.last_err is None: self.last_err =
            ## err` below drives the P factor `(err_last/err)^k_P` to exactly 1.
            ## That is the textbook response to a rejection (Hairer & Wanner II.4),
            ## and it needs no new mode -- the machinery was already here.
            h_next = h_curr * max(MIN_SHRINK_RATIO, 0.9 * (1.0 / err)**exponent)
            self.last_err = None
            return False, h_next
            
        # Step accepted: use PI update
        if self.last_err is None:
            self.last_err = err
            
        # Standard PI formula for step size control.  err is already normalized so
        # that err==1 is the target (TRTOL is folded into etol above).
        factor = self.pi_factor(err, self.last_err, p)
        
        h_next = h_curr * factor
        h_next = min(h_next, max_step)
        
        self.last_err = err
        
        return True, h_next


class SolutionLTEController(StepController):
    """Fang's LTE (DAC 2013 eq 6): computed solution minus polynomial extrapolation.

    STAGE 12B.  Every other controller here estimates the truncation error from
    divided differences of the CHARGE vector and converts to solution units by
    solving against ``J``.  Fang's is a different estimator: extrapolate the
    accepted solution history to the new time point with a polynomial, and take
    the largest deviation of the computed solution from it.

        eps_m = | v_i(t_m) - v_{i,extrapolated} |

    ``i`` is the *controlling LTE node*, and the paper is explicit that it "may
    vary from time point to time point" -- so it is recorded on the controller
    after every evaluation rather than assumed fixed.

    **Why the estimator was changed rather than reformulated.** The charge form
    divides by ``h`` repeatedly, so as the step shrinks it amplifies rounding in
    ``q`` and the estimate GROWS.  That is harmless when the error is only tested
    against a threshold, and disqualifying when it is the equation being solved
    for ``h``: Newton walks the step size down until it underflows, measured at
    gate 12B-0.  This form differences two solution values and falls as the step
    falls.  See ``doc/fang_dac2013_math.md`` sec. 3 and
    ``test_solution_lte.py::test_it_does_not_blow_up_as_the_step_shrinks``.

    **Two things here are ours, not the paper's** (it does not specify either;
    see the extraction doc sec. 6):

      * ``tau_m``.  Reused from the rest of this module -- ``TRTOL * (reltol*ref
        + abstol)`` with ``relref`` choosing ``ref`` -- so this controller is
        comparable with `IntegralController` at the same settings rather than
        being scored on a different tolerance.
      * the extrapolation degree, taken as ``len(x_hist)-1`` capped at 2, which
        leaves an ``O(h^3)`` deviation to match the second-order integrators.
    """

    ## The controlling LTE node of the most recent evaluation, or None.  Public
    ## because Fang's `q^T = df_lte/dv` is a signed unit vector on exactly this
    ## index -- which is what makes the coupled system cheap to form.
    controlling_index = None

    ## Degree cap: a degree-2 polynomial through three accepted points leaves an
    ## O(h^3) deviation, the order of the trapezoidal and Gear-2 LTE.  Raising it
    ## without also raising the integrator order would measure a truncation the
    ## method does not commit.
    MAX_DEGREE = 2

    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, no_history, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0, n_nodes=None, h_last2=None, h_clamped=False, x_hist=None):
        from pycircuit.circuit._lte_kernels import (solution_lte,
                                                    step_for_error_ratio)

        self.last_err = None
        self.controlling_index = None

        ## Needs at least two past points to extrapolate with, i.e. a degree-1
        ## line.  With fewer there is nothing to compare against and the step is
        ## accepted unevaluated, exactly as the charge-based controllers do on
        ## their opening step.
        hs = [h for h in (h_last, h_last2) if h is not None]
        if no_history or not x_hist or len(x_hist) < 2 or not hs:
            return True, h_curr

        ## THE DEGREE MUST EQUAL THE METHOD ORDER.  Fang's eq (6) is the Milne
        ## device: the difference between the corrector's solution and an
        ## explicit predictor of the SAME order is a truncation-error estimate
        ## with a known constant.  A predictor that is more accurate than the
        ## corrector does not give a smaller estimate, it gives an erratic one --
        ## the deviation becomes the corrector's own error minus a much smaller
        ## number, and the two stop cancelling in any controlled way.  Measured
        ## with a fixed degree of 2 against the default order-1 backward Euler:
        ## the normalised error moved 600x for a 2x step change (0.005 -> 2.86)
        ## where the node-polynomial model says 4x, and the controller then
        ## oscillated accept/reject on alternate steps -- 2911 rejections against
        ## 3004 accepted on `rc-vsin` at reltol 1e-6.
        ##
        ## `active_integrator` is the one actually running this step, so an
        ## order drop to Euler drops the predictor degree with it.
        order = getattr(active_integrator, 'ORDER', None)
        if order is None:
            return True, h_curr
        degree = min(order, len(x_hist) - 1, len(hs), self.MAX_DEGREE)
        if degree < 1:
            return True, h_curr
        v_hist = list(x_hist[:degree + 1])
        h_hist = list(hs[:degree])

        lte = solution_lte(x_curr, v_hist, h_hist, h_curr)

        ref = self._reference(x_curr, x_last, no_history, n_nodes, toolkit)
        etol = TRTOL * (reltol * ref + abstol)

        err_array = abs(lte) / etol
        ## The reference node is held at zero by construction, so its deviation
        ## is identically zero and cannot be the controlling one; taking the
        ## argmax over the full vector is safe and keeps the index in the
        ## caller's numbering rather than a reduced one.
        self.controlling_index = int(np.argmax(err_array))
        err = float(err_array[self.controlling_index])
        self.last_err = err

        safety = 0.9
        target = self._band_target(safety, degree + 1)

        ## THE PREDICTION IS AN INVERSION, NOT A POWER LAW.  The deviation scales
        ## with h(h+h1)(h+h1+h2), which is ~h^(k+1) only while h >> h1 and goes
        ## LINEAR in h once h << h1.  Using (target/err)**(1/(k+1)) everywhere --
        ## the law the charge estimators use -- mis-predicts by orders of
        ## magnitude on the shrinking side: measured at 2049 rejections against
        ## 2143 accepted steps on rc-vsin at reltol 1e-6, the controller
        ## oscillating between a step far under tolerance and one far over.
        def predict(ratio):
            return step_for_error_ratio(h_curr, h_hist, ratio,
                                        MIN_SHRINK_RATIO, MAX_GROWTH_RATIO)

        if err > self.lte_gamma_max:
            return False, predict(target / err)

        if (err < self.lte_gamma_min and not h_clamped
                and h_curr < max_step * (1.0 - 1e-12)):
            h_next = min(self._damp(predict(target / max(err, 1e-12)), h_curr),
                         max_step)
            if h_next > h_curr * (1.0 + 1e-9):
                return False, h_next

        h_next = min(self._damp(predict(target / max(err, 1e-12)), h_curr),
                     max_step)
        return True, h_next

    def lte_gradients(self, x_curr, x_hist, h_hist, h_curr, etol):
        """Fang's ``q^T`` and ``d`` for the LTE equation, both in closed form.

        Returns ``(index, q_value, d)`` for

            f_lte(v_m, h_m) = |v_i - P_i(h_m)| / etol_i - target

        where ``i`` is the controlling LTE node and ``P`` is the extrapolation
        from the accepted history.

        **``q^T`` is a single nonzero.**  ``P`` is built from PAST solutions
        only, so it does not depend on ``v_m`` at all, and

            q^T = d f_lte / d v_m = sign(v_i - P_i) / etol_i  on column ``i``,
                                    zero everywhere else.

        Returned as ``(index, value)`` rather than a dense row because a dense
        ``R^{1xN}`` of zeros is a poor way to say "one entry", and the bordered
        solve wants the index anyway.

        **``d`` is the derivative of the extrapolation.**  Holding ``v_m`` fixed,

            d = d f_lte / d h_m = -sign(v_i - P_i) P'_i(h_m) / etol_i.

        Together these are what section 3.2's "p, q^T and d can be computed
        explicitly with negligible computational costs" actually refers to; the
        claim is only true for an LTE in solution space, which is why it was a
        useful signal that the charge-based estimator was the wrong quantity.

        TWO APPROXIMATIONS, both deliberate and both ours rather than the
        paper's:

          * ``etol`` is treated as constant.  Under `relref='pointlocal'` it
            depends on ``|x_curr|`` and so contributes a second term; under the
            default `sigglobal` it is a running maximum over the whole run and
            is genuinely constant with respect to this time point's unknowns.
          * the coupling of ``v_m`` to ``h_m`` through the circuit equations is
            NOT in ``d``.  That is correct: eq (12) is a partial-derivative
            block system, and that coupling is exactly what ``p`` and ``J``
            carry.
        """
        from pycircuit.circuit._lte_kernels import extrapolate_with_derivative

        P, dP = extrapolate_with_derivative(x_hist, h_hist, h_curr)
        dev = x_curr - P
        scaled = abs(dev) / etol
        i = int(np.argmax(scaled))

        ## `sign` of exactly zero would make both gradients vanish and the
        ## bordered row degenerate; +1 is the correct one-sided choice because
        ## the quantity being differentiated is an absolute value at its minimum.
        s = 1.0 if dev[i] >= 0.0 else -1.0
        return i, s / etol[i], -s * dP[i] / etol[i]
