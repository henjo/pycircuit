from abc import ABC, abstractmethod
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
    ##   pointlocal  each unknown against itself, now.  (Previous behaviour.)
    ##   alllocal    each unknown against its OWN largest value so far.
    ##   sigglobal   each unknown against the largest value of ANY unknown so far.
    ##
    ## Default stays `pointlocal` here so this change is inert until asked for;
    ## whether to adopt Spectre's default belongs with stage 4's `lteratio` work.
    relref = 'pointlocal'

    def set_relref(self, relref):
        if relref not in RELREF_MODES:
            raise ValueError(
                "relref must be one of %r, not %r" % (RELREF_MODES, relref))
        self.relref = relref
        self._ref_running = None
        return self

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
    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, no_history, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0, n_nodes=None, h_last2=None):
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
    
    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, no_history, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0, n_nodes=None, h_last2=None):
        ## No past point exists yet, so there is nothing to difference and the
        ## step is accepted unevaluated.  This is the only place in a run where
        ## that is correct, and it costs one uncontrolled step of O(h^2) Euler
        ## error at max_step -- which is why it used to dominate every accuracy
        ## measurement when breakpoints re-armed it periodically.
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
        except Exception:
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

        # Step-size prediction exponent is 1/(order+1); compute_lte returns that
        # (order+1) as ``p`` (2 for Euler, 3 for the 2nd-order methods).
        exponent = 1.0 / p
        safety = 0.9

        # --- STEP REJECTION / ACCEPTANCE ---
        if err > 1.0:
            # Step rejected: shrink and recalculate.
            h_next = h_curr * max(MIN_SHRINK_RATIO, safety * (1.0 / err)**exponent)
            return False, h_next

        # Step accepted: predict next size with the same controller law.
        h_next = h_curr * min(MAX_GROWTH_RATIO,
                              safety * (1.0 / max(err, 1e-12))**exponent)
        h_next = min(h_next, max_step)

        return True, h_next

class PIController(StepController):
    """
    Proportional-Integral Step Controller.
    Uses history of truncation error to provide smoother step size changes.
    """
    def __init__(self, k_i=0.3, k_p=0.4):
        self.k_i = k_i
        self.k_p = k_p
        self.last_err = None
        
    def evaluate_step(self, x_curr, x_last, q_curr, q_last_hist, iq_last_hist, h_curr, h_last, no_history, J, active_integrator, irefnode, reltol, abstol, toolkit, max_step, TRTOL=7.0, n_nodes=None, h_last2=None):
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
        except Exception:
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
            h_next = h_curr * max(MIN_SHRINK_RATIO, 0.9 * (1.0 / err)**exponent)
            return False, h_next
            
        # Step accepted: use PI update
        if self.last_err is None:
            self.last_err = err
            
        # Standard PI formula for step size control.  err is already normalized so
        # that err==1 is the target (TRTOL is folded into etol above).
        err_norm = max(err, 1e-12)
        err_last_norm = max(self.last_err, 1e-12)

        factor = (err_norm ** (-self.k_i)) * ((err_last_norm / err_norm) ** self.k_p)
        
        # Limit the step size change
        factor = min(MAX_GROWTH_RATIO, max(MIN_SHRINK_RATIO, factor))
        
        h_next = h_curr * factor
        h_next = min(h_next, max_step)
        
        self.last_err = err
        
        return True, h_next
