from typing import NamedTuple, Any
import numpy as np
import jax
import warnings

import jax.numpy as jnp

from pycircuit.circuit._lte_kernels import (bdf2_companion, euler_companion,
                                            trapezoidal_companion,
                                            second_divided_difference,
                                            euler_lte, gear2_lte)
from pycircuit.circuit.analysis import Analysis
from pycircuit.circuit.nrsolver import NoConvergenceError
from pycircuit.utilities.param import Parameter
from pycircuit.circuit.analysis import CircuitResult as Result

## Depth of the per-TLine delay-line ring buffer, in accepted steps.  It was the
## bare literal `10000` in seven places, three of them inside modulo arithmetic,
## which made it impossible to see that the buffer is a RING: past this many
## accepted steps `tline_head` wraps and `interp_tlines` starts interpolating
## against entries from a previous lap.  `cond_fun` stops at the end of the buffer
## and returns whatever is there, so the result is a plausible wrong waveform with
## no error and no warning.  `JAXTransient.solve` now checks for the wrap and
## raises; sizing the buffer from the run instead of fixing it belongs to stage 9.
TLINE_HISTORY_DEPTH = 10000


class NewtonState(NamedTuple):
    x: Any
    xdiff: Any
    F_norm: Any
    iters: Any

    ## STAGE 9(e).  The loop exits on `F_norm <= conv_tol` OR `iters >= maxiter`
    ## and the caller took `x` either way, so an unconverged iterate was committed
    ## as the step's solution and its LTE computed from it.  Measured on a LINEAR
    ## RC at maxiter=1: 4.97e-2 max error against 4.28e-3, reported as nothing at
    ## all.  A traced loop cannot raise, so the fact travels out as a flag.
    converged: Any = True

class TransientState(NamedTuple):
    t: Any
    dt: Any
    step_idx: Any
    
    x_history: Any     
    q_history: Any     
    iq_history: Any    
    h_history: Any     
    
    results_buffer: Any
    time_buffer: Any   
    
    tline_history: Any  # (N_tlines, buf_size, 5) -> t, v1, v2, i1, i2
    tline_head: Any     # int32 scalar

    ## STAGE 9(c) -- the running maximum of |x| over ALL past steps and ALL
    ## unknowns.  This is Spectre's `sigglobal`, which `Transient` already ships as
    ## its default `relref`, and it is what makes an absolute LTE floor of 1e-12
    ## safe: under `pointlocal` -- each unknown against itself, now -- a node
    ## carrying no signal drives `ref -> 0`, the tolerance collapses to the floor,
    ## and the controller chases numerical noise on a quiet node.  The CPU hid that
    ## for a while by raising the floor a millionfold; `sigglobal` is the fix, and
    ## this field is what carries it through a traced loop.
    sig_max: Any = 0.0

    ## STAGE 9, gate 9-1(b) -- how many steps the controller REJECTED.  Nothing
    ## reported this, so "a step is actually rejected" -- one of the three CPU
    ## step-control gates -- was not expressible on this backend, which is the
    ## asymmetry that let a copied LTE defect survive being fixed twice.  A
    ## rejected step advances neither `t` nor `step_idx`, so it leaves no trace in
    ## the output buffers; it has to be counted where it happens.
    n_rejected: Any = 0

    ## Steps rejected because the Newton solve did not converge, and steps
    ## ACCEPTED anyway because `dt` was already at the floor.  The second is the
    ## one to read: it counts places where the returned waveform is not a solution
    ## of the circuit equations, which is exactly what 9(e) found happening
    ## silently.
    n_nonconverged: Any = 0
    n_forced_nonconverged: Any = 0

    ## F19(c) (doc/transient_review_260820.md): a step at the dt floor whose
    ## Newton CONVERGED but whose LTE still failed used to be accepted with no
    ## trace at all -- the CPU counts these as force_accepts and warns that
    ## the accepted error is unbounded.  Counted here for the same reason.
    n_forced_lte: Any = 0

    ## F11/F19: "do not trust a 2nd-order polynomial through this point" --
    ## set when the PREVIOUS accepted step landed on a breakpoint, consumed by
    ## effective_first_order.  An explicit flag rather than an h_history
    ## sentinel, per the review's revised recommendation: a sentinel in the
    ## step history would also falsify the LTE-exemption and predictor guards
    ## that key on "no history yet", which are different claims -- the same
    ## conflation the CPU path documents having paid for.
    force_first_order: Any = False

# ---------------------------------------------------------------------------
# Phase 1: Pure Functional Integrators
# ---------------------------------------------------------------------------

## STAGE 9(a) -- these three were transcriptions of `integrator.py`'s.  They now
## call the one definition in `_lte_kernels`, which is plain arithmetic and so
## traces under `jax.jit` unchanged.
##
## Euler and Gear-2 were already bit-for-bit the same expression, so those runs do
## not move.  TRAPEZOIDAL DID DIFFER: this file computed `(2/dt) * (q_n - q_{n-1})`
## where `integrator.py` computes `2 * (q_n - q_{n-1}) / dt` -- algebraically equal,
## but two roundings against one, so the last bits differed between the backends.
## Unifying necessarily moves one of them; it moves this one, onto the form with the
## fewer roundings.  `eval_method` is hard-coded to 'gear' at both call sites, so no
## production path reaches `trap_step` today.

def backward_euler_step(q_curr, C_curr, q_prev, dt):
    return euler_companion(q_curr, C_curr, q_prev, dt)


def trap_step(q_curr, C_curr, q_prev, iq_prev, dt):
    return trapezoidal_companion(q_curr, C_curr, q_prev, iq_prev, dt)


def gear2_step(q_curr, C_curr, q_prev1, q_prev2, dt_curr, dt_prev):
    return bdf2_companion(q_curr, C_curr, q_prev1, q_prev2, dt_curr, dt_prev)

def effective_first_order(state: TransientState):
    """Is this step effectively integrated at order 1?

    True while the run has fewer than two completed steps recorded
    (``h_history[0] == 0`` or ``h_history[1] == 0`` -- RUN-global facts, the
    buffers carry across chunk boundaries) and when the previous accepted
    step landed on a breakpoint (``force_first_order``, F11).  ONE
    definition, consumed by the integration dispatch, the LTE estimator, and
    the step-size exponent alike -- F19's point
    (doc/transient_review_260820.md): the integration used to fall back to
    Euler dynamically while the estimator stayed statically Gear-2, so an
    order-1 step was scored by the order-2 formula on a zero-seeded history
    with the wrong exponent handed to the step controller.

    Deliberately NOT ``step_idx < 2``: step_idx is CHUNK-local and resets at
    every chunk boundary, so the old predicate re-dropped the order (and,
    once the estimator followed it, re-scored steps) at each boundary --
    measured as chunking changing a 59-step run to 61 steps.  The history
    buffers are the run-global truth.
    """
    no_history = jnp.logical_or(state.h_history[0] == 0.0,
                                state.h_history[1] == 0.0)
    return jnp.logical_or(no_history, state.force_first_order)


def compute_integration(q_curr, C_curr, state: TransientState, method='trap',
                        first_order=None):

    q_prev = state.q_history[0]
    iq_prev = state.iq_history[0]
    dt = state.dt

    def do_euler():
        return backward_euler_step(q_curr, C_curr, q_prev, dt)

    def do_trap():
        return trap_step(q_curr, C_curr, q_prev, iq_prev, dt)

    def do_gear2():
        q_prev2 = state.q_history[1]
        dt_prev = state.h_history[0]
        fallback = (effective_first_order(state) if first_order is None
                    else first_order)

        def _euler():
            return do_euler()

        def _gear():
            return gear2_step(q_curr, C_curr, q_prev, q_prev2, dt, dt_prev)

        return jax.lax.cond(fallback, _euler, _gear)

    if method == 'euler':
        return do_euler()
    elif method in ('gear', 'gear2'):
        return do_gear2()
    else:
        return do_trap()

# ---------------------------------------------------------------------------
# Phase 2: Inner Newton Loop
# ---------------------------------------------------------------------------

def extrapolate_predictor(state: TransientState):
    x1 = state.x_history[0]
    x2 = state.x_history[1]
    dt = state.dt
    dt_old = state.h_history[0]
    
    dt_old_safe = jnp.where(dt_old == 0.0, 1.0, dt_old)
    x_pred = x1 + (dt / dt_old_safe) * (x1 - x2)
    x_pred = jnp.where(dt_old == 0.0, x1, x_pred)
    return x_pred

def newton_inner_loop(state: TransientState, circuit, irefnode, tline_params, tline_indices, eval_method='euler', reltol=1e-3, abstol=1e-6, xtol=1e-12, maxiter=50, params_tree=None, max_dv=jnp.inf):
    x_init = extrapolate_predictor(state)
    def interp_tlines(t_target, params, history, head):
        # t_target is scalar, history is (10000, 5)
        def find_idx(idx, _):
            curr_t = history[(head - idx) % TLINE_HISTORY_DEPTH, 0]
            return jax.lax.cond(curr_t <= t_target, lambda: idx, lambda: idx + 1)
    
        # Simple while loop to find the interval
        def cond_fun(val):
            idx = val
            curr_t = history[(head - idx) % TLINE_HISTORY_DEPTH, 0]
            # stop if curr_t <= t_target or we hit max history (10000)
            return jnp.logical_and(curr_t > t_target, idx < TLINE_HISTORY_DEPTH - 1)
        
        def body_fun(val):
            return val + 1
        
        idx1 = jax.lax.while_loop(cond_fun, body_fun, 0)
        idx0 = jnp.maximum(0, idx1 - 1)
    
        idx1_mapped = (head - idx1) % TLINE_HISTORY_DEPTH
        idx0_mapped = (head - idx0) % TLINE_HISTORY_DEPTH
    
        val1 = history[idx1_mapped]
        val0 = history[idx0_mapped]
    
        t1 = val1[0]
        t0 = val0[0]
    
        # linear interpolation
        frac = jnp.where(t0 != t1, (t_target - t1) / (t0 - t1), 0.0)
        interp_val = val1 + frac * (val0 - val1)
    
        v1_past = interp_val[1]
        v2_past = interp_val[2]
        i1_past = interp_val[3]
        i2_past = interp_val[4]
    
        Z0 = params[1]
        e1 = v2_past + Z0 * i2_past
        e2 = v1_past + Z0 * i1_past
    
        return e1, e2

    def apply_tlines(I_u, t_curr):
        n_tlines = tline_params.shape[0]
        if n_tlines == 0:
            return I_u
        t_targets = t_curr - tline_params[:, 0]
        e1s, e2s = jax.vmap(interp_tlines, in_axes=(0, 0, 0, None))(t_targets, tline_params, state.tline_history, state.tline_head)
    
        I_u = I_u.at[tline_indices[:, 4]].add(-e1s)
        I_u = I_u.at[tline_indices[:, 5]].add(-e2s)
        return I_u


    ## F6(b) (doc/transient_review_260820.md): THE CPU'S PER-ROW CRITERIA,
    ## replacing `conv_tol = abstol + reltol * F_norm0` -- which was wrong
    ## three ways at once.  Flavour: the scalar floor threaded here was
    ## vabstol, a voltage, against a residual of KCL currents (F6(a) picked
    ## the majority flavour; the per-row vectors finish the job -- iabstol on
    ## node rows, vabstol on branch rows, and the transposed pair for the
    ## update test).  Reference: relative to the INITIAL residual, so a bad
    ## predictor loosened the target -- false convergence -- while a good
    ## one collapsed it toward the absolute floor -- spurious
    ## non-convergence.  Norm: a summed L1, so one badly-failed row diluted
    ## with circuit size.  The body below mirrors StandardNewton: per-row
    ## conv_f against reltol*I_scale + abstol, per-row conv_x against the
    ## step actually taken, both computed on the consistent (F(x), dx) pair
    ## -- which also retires stage 9(e)'s trailing re-evaluation, since the
    ## converged flag now travels in state and refers to the pair that
    ## produced the returned x, exactly as on the CPU.
    def cond_fun(nr_state: NewtonState):
        return jnp.logical_and(jnp.logical_not(nr_state.converged),
                               nr_state.iters < maxiter)

    def body_fun(nr_state: NewtonState):
        x = nr_state.x

        I_G = circuit.i(x, params_tree=params_tree)
        G_G = circuit.G(x, params_tree=params_tree)
        q_C = circuit.q(x, params_tree=params_tree)
        C_C = circuit.C(x, params_tree=params_tree)
        I_u = circuit.u(state.t + state.dt, analysis='tran', params_tree=params_tree)
        I_u = apply_tlines(I_u, state.t + state.dt)

        i_C, G_C_eq = compute_integration(q_C, C_C, state, method=eval_method)

        F = I_G + i_C + I_u
        J = G_G + G_C_eq

        J_sub = jnp.delete(jnp.delete(J, irefnode, axis=0), irefnode, axis=1)
        F_sub = jnp.delete(F, irefnode)
        xdiff_sub = jnp.linalg.solve(J_sub, -F_sub)
        xdiff = jnp.insert(xdiff_sub, irefnode, 0.0)

        ## F16 (doc/transient_review_260820.md): the update clamp is a
        ## PARAMETER, default disabled (jnp.inf -> alpha 1.0).  The hardcoded
        ## 0.5 V made any swing beyond ~maxiter*0.5 V non-convergent by
        ## construction with a false "failed to converge" diagnosis -- a 48 V
        ## rail was unreachable in 100 iterations -- and under the per-row
        ## criteria (F6(b)) it even made a LINEAR 1 V step need multiple
        ## clamped iterations.  A flat voltage clamp punishes the linear part
        ## of the circuit for the nonlinear part's sins; per-junction
        ## limiting (pnjlim's shape) is the eventual replacement if the JAX
        ## element set grows junctions.
        max_diff = jnp.max(jnp.abs(xdiff))
        alpha = jnp.where(max_diff > max_dv, max_dv / max_diff, 1.0)
        step = alpha * xdiff
        x_next = x + step

        I_scale = jnp.abs(J) @ jnp.abs(x_next) + jnp.abs(F)
        conv_f = jnp.all(jnp.abs(F) < reltol * I_scale + abstol)
        conv_x = jnp.all(jnp.abs(step)
                         < reltol * jnp.maximum(jnp.abs(x_next), jnp.abs(x))
                         + xtol)

        return NewtonState(x=x_next, xdiff=step, F_norm=jnp.sum(jnp.abs(F)),
                           iters=nr_state.iters + 1,
                           converged=jnp.logical_and(conv_f, conv_x))

    initial_nr_state = NewtonState(x=x_init, xdiff=jnp.zeros_like(x_init),
                                   F_norm=jnp.asarray(jnp.inf), iters=0,
                                   converged=jnp.asarray(False))
    return jax.lax.while_loop(cond_fun, body_fun, initial_nr_state)

# ---------------------------------------------------------------------------
# Phase 3: Outer Time Loop & Adaptive Control
# ---------------------------------------------------------------------------

def ywr_error_ratio(i_curr, x_curr, J, state: TransientState, irefnode,
                    method='trap', trtol=7.0, lte_rel=1e-4, lte_abstol=1e-12,
                    first_order=None):
    """Yao-Wang-Roychowdhury DAE LTE, returned as a normalized error ratio.

    Mirrors the CPU transient's estimator: forms the residual as a
    difference of the charge derivative ``g = dq/dt`` (the companion current),
    maps it to the solution space with ``J^-1`` (the DAE Jacobian factor), and
    normalizes by a voltage-domain tolerance.  ``g_n`` is the just-computed
    companion current, ``g_{n-1}``/``g_{n-2}`` come from the iq history.

    Returns ``(error_ratio, order_plus_one)``.
    """
    dt = state.dt
    dt_prev = jnp.where(state.h_history[0] == 0.0, dt, state.h_history[0])
    g_n = i_curr
    g_nm1 = state.iq_history[0]
    g_nm2 = state.iq_history[1]

    if method == 'euler':
        Eg = euler_lte(g_n, g_nm1)
        order_p1 = 2.0
    elif method in ('gear', 'gear2'):
        ## F19: THE ESTIMATOR FOLLOWS THE EFFECTIVE ORDER.  `method` is a
        ## static string, but the integration falls back to Euler dynamically
        ## (opening steps; the F11 breakpoint sentinel), and an order-2
        ## formula scoring an order-1 step differences a zero-seeded or
        ## corner-straddling g-history with the wrong exponent downstream --
        ## the CPU's estimator follows the dropped order, and now this one
        ## does too.
        first = (effective_first_order(state) if first_order is None
                 else first_order)

        def _euler_branch(_):
            return euler_lte(g_n, g_nm1), jnp.asarray(2.0)

        def _gear_branch(_):
        ## STAGE 9, gate 9-1(a) -- THE 3/4 OPTIMISM, FOUND A THIRD TIME.
        ##
        ## This was YWR's Table I GEAR2 residual,
        ##     -(1/8) ((h1+h2)/(h1 h2)) (h2 g_n - (h1+h2) g_{n-1} + h1 g_{n-2}),
        ## which on a uniform grid reduces to -(1/4) h^2 q''' against a true BDF-2
        ## local truncation error of -(1/3) h^2 q'''.  So it reported 3/4 of the
        ## error at every step -- the solver was 25% optimistic about its own
        ## accuracy, on the ONLY eval_method either entry point uses.
        ##
        ## The CPU found and fixed this in stage 4i and the fix never crossed to
        ## this file, which is precisely the divergence stage 9 exists to close:
        ## the same defect has now been found three times in two transcriptions.
        ## Measured here at 2.5e5 against the CPU's 3.333e5 for q''' = 1e6.
        ##
        ## The form below is 4i's: estimate q''' from the second divided difference
        ## of the companion current and multiply by the method's own error
        ## constant, so the coefficient is derived rather than transcribed.
            h1, h2 = dt, dt_prev
            ## `second_divided_difference` returns q'''/2 and `gear2_lte`
            ## wants q'''/6, hence the /3.  Both live in `_lte_kernels` so the
            ## normalisation is stated once -- pairing an estimate with a
            ## differently-normalised constant is exactly how the 3/4 optimism
            ## survived in two places.
            q3_over_6 = second_divided_difference(g_n, g_nm1, g_nm2, h1, h2) / 3.0
            return gear2_lte(h1, h2, q3_over_6), jnp.asarray(3.0)

        Eg, order_p1 = jax.lax.cond(first, _euler_branch, _gear_branch, None)
    else:  # trapezoidal
        Eg = -(1.0 / 6.0) * (g_n - 2.0 * g_nm1 + g_nm2)
        order_p1 = 3.0

    # lte = J^-1 Eg in the reduced (reference-node-removed) space.
    J_r = jnp.delete(jnp.delete(J, irefnode, axis=0), irefnode, axis=1)
    Eg_r = jnp.delete(Eg, irefnode)
    lte_r = jnp.linalg.solve(J_r, Eg_r)
    lte = jnp.insert(lte_r, irefnode, 0.0)

    ## `sigglobal`: the reference is the largest signal ANY unknown has reached so
    ## far, not this unknown right now.  `state.sig_max` carries it; the
    ## `maximum(..., |x_curr|)` keeps the very first step -- where sig_max is still
    ## the seed -- from measuring against zero.
    ref = jnp.maximum(state.sig_max, jnp.max(jnp.abs(x_curr)))
    ## `lte_abstol` is a per-row VECTOR (volts on node rows, amps on branch rows),
    ## built by `JAXTransient._lte_abstol`.  A scalar here applies one physical
    ## kind of tolerance to rows of another -- 0.3a's residual-vs-solution defect,
    ## which 9(c) exists to avoid repeating.  A scalar still broadcasts, so a
    ## direct caller that passes one gets the old behaviour rather than an error.
    etol = trtol * (lte_rel * ref + lte_abstol)
    return jnp.max(jnp.abs(lte) / etol), order_p1


def collect_breakpoints(cir, tend):
    """Every source discontinuity in ``(0, tend]``, plus ``tend`` itself.

    STAGE 9(d).  This replaces two copies that disagreed.  ``solve`` iterated
    ``for elem in cir.elements`` -- a **dict**, so it yielded string keys and
    ``hasattr('V1', 'next_event')`` was False, giving **0 breakpoints, always**.
    ``solve_batched`` iterated ``.items()`` and was correct, and therefore hit the
    second bug instead: ``Pulse.next_event`` returned ``t`` itself at ``t = 0``, so
    the enumeration never advanced and the call **hung**.  Two bugs cancelling, one
    per copy, which is the argument for having one copy.

    ``Pulse.next_event`` is fixed at the source, but the progress guard below stays:
    it is the difference between a wrong breakpoint list and a wall-clock hang, and
    only one of those is diagnosable from a stack trace.
    """
    ## Bounded per element (review hygiene item, Phase 3): a periodic source
    ## enumerates 4*f*tend events -- a 1 GHz sine over 1 ms is 4 million list
    ## entries, materialised before the first step runs.  Past the cap the
    ## element's later events are simply not breakpoints (the adaptive
    ## controller still resolves them, at rejection cost), and the run says
    ## so instead of silently thrashing.
    MAX_EVENTS_PER_ELEMENT = 10000
    breakpoints = []
    for _inst_name, elem in cir.elements.items():
        if not hasattr(elem, 'next_event'):
            continue
        n_events = 0
        t_event = 0.0
        while t_event < tend:
            if n_events >= MAX_EVENTS_PER_ELEMENT:
                warnings.warn(
                    'transient: %r produced %d breakpoints by t=%g of tend=%g; '
                    'later events from this element are not breakpoint-'
                    'truncated (the step controller still resolves them, at '
                    'extra rejections). A longer timestep or shorter tend '
                    'avoids this.' % (type(elem).__name__,
                                     MAX_EVENTS_PER_ELEMENT, t_event, tend),
                    RuntimeWarning)
                break
            nxt = elem.next_event(t_event)
            if nxt is None:
                break
            nxt = float(nxt)
            ## The guard.  `next_event` is contracted to return strictly more than
            ## its argument (see `pycircuit.circuit.func`); an implementation that
            ## breaks that contract stops this element's enumeration instead of
            ## spinning.  Warn rather than raise: one bad source should not take
            ## down a run whose other sources are fine.
            if not nxt > t_event:
                warnings.warn(
                    'transient: %r.next_event(%g) returned %g, which does not '
                    'advance; breakpoints for this element are incomplete. Its '
                    'next_event violates the strictly-increasing contract.'
                    % (type(elem).__name__, t_event, nxt),
                    RuntimeWarning)
                break
            t_event = nxt
            if t_event <= tend:
                breakpoints.append(t_event)
                n_events += 1
            else:
                break
    breakpoints.append(tend)
    return sorted(set(breakpoints))


def calculate_next_dt(dt, error_ratio, dt_min, dt_max, t_breaks_array, current_t, order_p1=2.0):
    ## THE CPU'S LAW, F17 (doc/transient_review_260820.md): aim at
    ## target = safety**p, not at err = 1.0 -- the rejection threshold
    ## itself.  Aiming at the edge meant every successor step that landed a
    ## hair over 1.0 was a rejected re-solve; measured on rc-vsin at
    ## reltol 1e-4/1e-6 as a 44%/40% rejection rate, cut to a few percent by
    ## the safety margin.  `(0.9**p / err)**(1/p) = 0.9 * err**(-1/p)` --
    ## identical arithmetic to a bare safety multiplier, written in the CPU's
    ## vocabulary (stepcontroller._band_target) so the backends read as one
    ## law.  The clamps are the CPU's named constants rather than re-derived
    ## literals: the shrink floor was 0.1 here against MIN_SHRINK_RATIO=0.2.
    from pycircuit.circuit.stepcontroller import (MIN_SHRINK_RATIO,
                                                  MAX_GROWTH_RATIO)
    safety = 0.9
    target = safety ** order_p1
    factor = jnp.where(error_ratio <= 0.0, MAX_GROWTH_RATIO,
                       (target / jnp.maximum(error_ratio, 1e-12))
                       ** (1.0 / order_p1))
    factor = jnp.clip(factor, MIN_SHRINK_RATIO, MAX_GROWTH_RATIO)
    dt_new = dt * factor
    dt_new = jnp.clip(dt_new, dt_min, dt_max)

    future_breaks = jnp.where(t_breaks_array > current_t + 1e-12, t_breaks_array, jnp.inf)
    next_break = jnp.min(future_breaks)

    dt_new = jnp.where(current_t + dt_new > next_break, next_break - current_t, dt_new)
    dt_new = jnp.maximum(dt_new, dt_min)
    return dt_new

def outer_time_loop(initial_state: TransientState, circuit, tend, chunk_size, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method='euler', params_tree=None, reltol=1e-4, abstol=1e-12, xtol=1e-12, maxiter=100, trtol=7.0, lte_reltol=1e-4, lte_abstol=1e-12, max_dv=jnp.inf):

    ## The same epsilon `calculate_next_dt` uses to decide a breakpoint is "already
    ## reached".  They disagreed: after 500 steps of 1e-5 the accumulated `t` sits
    ## ~1e-18 short of `tend`, which the breakpoint filter treats as arrived (so it
    ## drops the clamp) while `t < tend` treats as not arrived (so it takes another
    ## FULL step).  That is the residual overshoot -- exactly one timestep, measured
    ## at t[-1] = 5.010e-3 for a requested 5e-3.
    t_eps = 1e-12 * max(abs(float(tend)), 1.0)

    def time_cond(state: TransientState):
        under_time = state.t < tend - t_eps
        under_chunk = state.step_idx < chunk_size
        return jnp.logical_and(under_time, under_chunk)

    def time_body(state: TransientState):
        nr_state = newton_inner_loop(state, circuit, irefnode, tline_params, tline_indices,
                                     eval_method=eval_method, reltol=reltol, abstol=abstol,
                                     xtol=xtol, maxiter=maxiter,
                                     params_tree=params_tree, max_dv=max_dv)
        x_curr = nr_state.x
        q_curr = circuit.q(x_curr, params_tree=params_tree)
        C_curr = circuit.C(x_curr, params_tree=params_tree)
        ## F19: one effective-order flag for this step, consumed by the
        ## integration AND the estimator so they cannot disagree.
        first_order = effective_first_order(state)
        i_curr, Geq = compute_integration(q_curr, C_curr, state,
                                          method=eval_method,
                                          first_order=first_order)

        ## STAGE 9(f) -- ONE ESTIMATOR.  There used to be a charge-domain branch
        ## here, selected by `lte_formula='classic'`.  It was deleted rather than
        ## repaired: its tolerance applied `lte_abs = 1e-6`, a VOLTAGE floor, to a
        ## CHARGE -- one microcoulomb, against node charges of pico- to
        ## femtocoulombs -- so the normalized error could not reach 1 and no step
        ## was ever rejected.  The controller ran open-loop.  Under `solve()` that
        ## degenerates to a fixed-step run (`dt_max = timestep`), which is why it
        ## looked plausible; under `solve_batched` (`dt_max = tend/10`) it also
        ## costs accuracy.  0.2b measured the `J^-1` mapping this branch avoided at
        ## 1-3% of a step, well under its own 10% keep-it threshold.
        J = circuit.G(x_curr, params_tree=params_tree) + Geq
        error_ratio, order_p1 = ywr_error_ratio(
            i_curr, x_curr, J, state, irefnode,
            method=eval_method, trtol=trtol, lte_rel=lte_reltol,
            lte_abstol=lte_abstol, first_order=first_order)

        ## Accept when the LTE is within tolerance.  Always accept the first step
        ## (no history for the estimate) and when dt has reached the floor -- the
        ## latter guarantees forward progress so a rejection loop cannot deadlock.
        at_floor = state.dt <= dt_min * (1.0 + 1e-9)
        ## Run-global "genuinely nothing to difference against", NOT step_idx
        ## (chunk-local) and NOT force_first_order (an order-dropped step is
        ## error-controlled -- the CPU documents conflating those two as what
        ## made max_step a correctness knob).
        first = state.h_history[0] == 0.0

        ## STAGE 9(e) -- A NON-CONVERGED NEWTON IS NOT A SOLUTION.
        ##
        ## `error_ratio` is computed from `x_curr`; if Newton did not converge,
        ## `x_curr` does not satisfy the circuit equations and its truncation error
        ## is meaningless, so accepting on a small ratio accepts a number about a
        ## point the circuit never occupied.  SPICE's answer to a failed Newton is
        ## a SMALLER STEP, not a committed non-solution.
        ##
        ## `at_floor` still forces progress -- otherwise a circuit that cannot
        ## converge at any step size hangs instead of finishing -- but that path is
        ## counted and warned about after the run, exactly as the CPU's force-accept
        ## is.  The `first` step is no longer exempt: a first step that did not
        ## converge is the worst one to commit, since every later step extrapolates
        ## from it.
        lte_ok = jnp.logical_or(error_ratio <= 1.0, first)
        accept = jnp.logical_or(
            jnp.logical_and(nr_state.converged, lte_ok), at_floor)
        forced = jnp.logical_and(at_floor, jnp.logical_not(nr_state.converged))
        ## F19(c): converged Newton, failing LTE, at the floor -- accepted
        ## with an unbounded truncation error, which the CPU force-accept
        ## warns about and this path used to swallow silently.
        forced_lte = jnp.logical_and(
            at_floor, jnp.logical_and(nr_state.converged,
                                      jnp.logical_not(lte_ok)))

        def do_accept(_):
            ## `state.t + state.dt`, NOT `state.t`.  This step has been accepted, so
            ## the next one starts where this one ENDS; passing the old time sized the
            ## next step to reach the breakpoint from the PREVIOUS position and so
            ## overshot it by about one step.  Measured before the fix: t[-1] of
            ## 5.0559e-3 against a requested tend of 5e-3, where `Transient` lands on
            ## it exactly.  `tend` is itself in `t_breaks_array`, which is why the
            ## overshoot showed up at the end of every run.
            next_dt = calculate_next_dt(state.dt, error_ratio, dt_min, dt_max,
                                        t_breaks_array, state.t + state.dt, order_p1)

            x_hist_new = jnp.roll(state.x_history, shift=1, axis=0).at[0].set(x_curr)
            q_hist_new = jnp.roll(state.q_history, shift=1, axis=0).at[0].set(q_curr)
            iq_hist_new = jnp.roll(state.iq_history, shift=1, axis=0).at[0].set(i_curr)
            h_hist_new = jnp.roll(state.h_history, shift=1, axis=0).at[0].set(state.dt)

            res_buf = state.results_buffer.at[state.step_idx].set(x_curr)
            time_buf = state.time_buffer.at[state.step_idx].set(state.t + state.dt)

            n_tlines = tline_params.shape[0]
            def update_tlines():
                v1 = x_curr[tline_indices[:, 0]] - x_curr[tline_indices[:, 1]]
                v2 = x_curr[tline_indices[:, 2]] - x_curr[tline_indices[:, 3]]
                i1 = x_curr[tline_indices[:, 4]]
                i2 = x_curr[tline_indices[:, 5]]
                tline_data = jnp.stack([jnp.full(n_tlines, state.t + state.dt), v1, v2, i1, i2], axis=1)
                new_head = (state.tline_head + 1) % TLINE_HISTORY_DEPTH
                new_history = state.tline_history.at[:, new_head, :].set(tline_data)
                return new_history, new_head

            tline_history_new, tline_head_new = jax.lax.cond(
                n_tlines > 0, update_tlines,
                lambda: (state.tline_history, state.tline_head))

            return TransientState(
                t=state.t + state.dt, dt=next_dt, step_idx=state.step_idx + 1,
                x_history=x_hist_new, q_history=q_hist_new,
                iq_history=iq_hist_new, h_history=h_hist_new,
                results_buffer=res_buf, time_buffer=time_buf,
                tline_history=tline_history_new, tline_head=tline_head_new,
                ## Updated only on ACCEPT: a rejected step's x is not a signal the
                ## circuit ever had, and letting it raise the reference would loosen
                ## every later tolerance on the strength of a discarded iterate.
                sig_max=jnp.maximum(state.sig_max, jnp.max(jnp.abs(x_curr))),
                ## CARRIED, not defaulted.  `do_accept` builds a fresh
                ## TransientState rather than `_replace`-ing, so every field it
                ## omits silently reverts to the NamedTuple default -- and these are
                ## cumulative counters.  Omitting them reset the rejection count to
                ## zero on each accepted step, which made "16 configurations, zero
                ## rejections" a measurement of the bug rather than of the solver.
                n_rejected=state.n_rejected,
                n_nonconverged=state.n_nonconverged,
                n_forced_nonconverged=state.n_forced_nonconverged
                + jnp.where(forced, 1, 0),
                n_forced_lte=state.n_forced_lte
                + jnp.where(forced_lte, 1, 0),
                ## F11: the CPU's breakpoint discipline, ported.  A step that
                ## LANDS on a breakpoint (calculate_next_dt truncates onto
                ## them, so landing is exact to t_eps) marks the NEXT step
                ## "do not trust a 2nd-order polynomial through this point":
                ## effective_first_order then drops both the integration and
                ## the error estimate to order 1 for exactly one step, instead
                ## of differencing a g-history that straddles the corner --
                ## which cost a rejection burst at every edge (measured on a
                ## VPulse RC before this line: 38 rejections in 183 accepted
                ## steps, edge-synchronous).
                force_first_order=jnp.any(
                    jnp.abs((state.t + state.dt) - t_breaks_array) <= t_eps))

        def do_reject(_):
            ## LTE above tolerance: shrink the step (bounded below by dt_min) and
            ## retry the same time point without committing anything.
            ## Two different reasons to shrink, and they need different factors.
            ## The predictor's `ratio^(-1/p)` is only meaningful when `x_curr` IS a
            ## solution; after a failed Newton the ratio describes a point the
            ## circuit never occupied, so a fixed cut is used instead.
            lte_shrink = jnp.clip(error_ratio ** (-1.0 / order_p1), 0.1, 0.9)
            shrink = jnp.where(nr_state.converged, lte_shrink, 0.25)
            retry_dt = jnp.maximum(state.dt * shrink, dt_min)
            return state._replace(
                dt=retry_dt,
                n_rejected=state.n_rejected + 1,
                n_nonconverged=state.n_nonconverged
                + jnp.where(nr_state.converged, 0, 1))

        return jax.lax.cond(accept, do_accept, do_reject, None)

    return jax.lax.while_loop(time_cond, time_body, initial_state)

# ---------------------------------------------------------------------------
# Phase 4: Class Wrapper & Chunking Engine
# ---------------------------------------------------------------------------

class JAXTransientStatistics(object):
    """What a JAX transient run did, as opposed to what it returned.

    The subset of `transient.TransientStatistics` a traced `while_loop` can count.
    `rejected_steps` is the load-bearing one: a rejected step advances neither `t`
    nor `step_idx`, so it leaves no trace in the output buffers, and without this
    counter the CPU gate "a step is actually rejected" could not be stated on this
    backend at all.  That is the asymmetry stage 9 exists to remove -- it is why a
    copied LTE defect had to be found and fixed twice.
    """

    __slots__ = ('accepted_steps', 'rejected_steps', 'signal_max',
                 'nonconverged_steps', 'forced_nonconverged_steps',
                 'forced_lte_steps')

    def __init__(self):
        self.accepted_steps = 0
        self.rejected_steps = 0
        self.signal_max = 0.0
        self.nonconverged_steps = 0
        self.forced_nonconverged_steps = 0
        ## F19(c): the JAX analogue of the CPU's force_accepts -- converged
        ## Newton, failing LTE, accepted at the dt floor.
        self.forced_lte_steps = 0

    def __repr__(self):
        return ('<JAXTransientStatistics accepted=%d rejected=%d '
                'nonconverged=%d forced=%d signal_max=%.4g>'
                % (self.accepted_steps, self.rejected_steps,
                   self.nonconverged_steps, self.forced_nonconverged_steps,
                   self.signal_max))


class JAXTransient(Analysis):
    """
    Fully compiled GPU-native transient analysis.
    Implements adaptive LTE and Predictor-Corrector implicitly in the JAX kernel.
    """

    ## STAGE 9(b)/(c) -- TOLERANCES ARE SETTABLE, AND THEY ARE THE CPU'S.
    ##
    ## This class declared no tolerances at all, so `JAXTransient(cir, reltol=1e-6)`
    ## raised `KeyError: 'parameter reltol not in parameter dictionary'` and there
    ## was no supported way to ask for a tighter run.  The values were hard-coded
    ## at the kernel's own defaults (`reltol=1e-3, abstol=1e-6` for Newton;
    ## `trtol=7.0, lte_rel=1e-3, lte_abs=1e-6` for the LTE, never threaded from
    ## the caller at all).  Two of those disagreed with `Transient`'s shipped
    ## defaults, so the two backends were solving to different accuracies while
    ## presenting as the same analysis.
    ##
    ## THE NAMES AND FLAVOURS ARE `Transient`'s DELIBERATELY.  0.3a's residual-vs-
    ## solution defect on the CPU came from applying one scalar absolute tolerance
    ## to rows of different physical kinds; the plan's 9(c) warns in as many words
    ## that threading a scalar here would re-create it.  So the absolute
    ## tolerances are VECTORS built the same way `Transient` builds them --
    ## `lte_vabstol` on node rows, `lte_iabstol` on branch rows, because the LTE
    ## lives in the SOLUTION domain where nodes carry volts and branches carry
    ## amps.  `relref` is not offered yet: the CPU's `sigglobal` needs a reduction
    ## over the whole vector inside the traced loop, and that is 9(c)'s remaining
    ## work rather than something to fake with a scalar.
    parameters = Analysis.parameters + [
        Parameter(name='reltol', desc='Relative tolerance', unit='',
                  default=1e-4),
        Parameter(name='iabstol',
                  desc='Absolute current error tolerance', unit='A',
                  default=1e-12),
        Parameter(name='vabstol',
                  desc='Absolute voltage error tolerance for the Newton solve',
                  unit='V', default=1e-12),
        Parameter(name='lte_vabstol',
                  desc='Absolute voltage tolerance for the local truncation error',
                  unit='V', default=1e-12),
        Parameter(name='lte_iabstol',
                  desc='Absolute current tolerance for the local truncation error',
                  unit='A', default=1e-12),
        Parameter(name='TRTOL',
                  desc='Ratio used to compute LTE tolerances from the Newton '
                       'tolerance (Spectre calls this lteratio)',
                  unit='', default=7.0),
        Parameter(name='maxiter', desc='Maximum number of iterations', unit='',
                  default=100),
        ## F16: opt-in, OFF by default -- the old hardcoded 0.5 V made any
        ## swing beyond ~maxiter*0.5 V non-convergent by construction.
        Parameter(name='max_dv',
                  desc='Largest Newton update per iteration, as a voltage; a '
                       'convergence aid for stiff nonlinear circuits. None '
                       '(default) disables the clamp.',
                  unit='V', default=None),
        ## STAGE 9(g).  Ported from `Transient`, same name, default and validation.
        ## DECISION D2 -- the same knob as `Transient.max_step`, same name and same
        ## default, because 9(d) and 9(g) both exist because these two backends had
        ## drifted apart on exactly this kind of detail.
        Parameter(name='max_step',
                  desc='Largest accepted timestep; None means timestep. Raise it to '
                       'let the controller coast through quiescent intervals (SPICE '
                       'calls this TMAX and defaults it to tend/50)',
                  unit='s', default=None),
        Parameter(name='firststep',
                  desc='Size of the first timestep; None means timestep*1e-3. '
                       'The first step cannot be error-checked, so taking it at '
                       'dt_max lets its error dominate the whole run.',
                  unit='s', default=None),
    ]

    def __init__(self, cir, **kwargs):
        super().__init__(cir, **kwargs)
        self.toolkit = getattr(self.cir, 'toolkit', None)
        if not (self.toolkit and self.toolkit.supports('autodiff')):
            raise ValueError("JAXTransient requires the circuit to use _jaxtoolkit.py.")

        ## `nrsolver` and `scaler` are inherited from `Analysis` and this backend
        ## honours NEITHER -- the Newton loop is a `jax.lax.while_loop` with its own
        ## fixed algorithm, and there is no scaling step.  They used to be accepted
        ## in silence, which is the "thin advertised feature" 0.1c warns about and
        ## the same shape of defect as the `lte_formula` knob removed in 9(f).
        ## Rejecting them loudly is the honest option until the loop can dispatch.
        for unsupported in ('nrsolver', 'scaler'):
            if kwargs.get(unsupported) is not None:
                raise NotImplementedError(
                    "JAXTransient does not support %r: its Newton loop is a "
                    "traced jax.lax.while_loop with a fixed algorithm and no "
                    "scaling step. It was previously accepted and ignored. Use "
                    "Transient for a circuit that needs %r."
                    % (unsupported, unsupported))

    def _max_step(self, timestep):
        """The clamp on how large an accepted step may grow -- decision D2.

        `Transient` has the identical method; a `max_step` below `timestep` is a
        caller error rather than a silent override, because the step can never
        exceed it and accepting one quietly would discard the timestep asked for.
        """
        max_step = getattr(self.par, 'max_step', None)
        if max_step is None:
            return timestep
        if max_step < timestep:
            raise ValueError(
                "max_step (%g) is smaller than timestep (%g); pass a smaller "
                "timestep instead, or max_step=None." % (max_step, timestep))
        return float(max_step)

    def _opening_step(self, timestep):
        """The size of the first step of a run.

        STAGE 9(g), ported from `Transient._opening_step`, which records the
        measurement: a run opening at `timestep` -- which is also `dt_max`, the
        largest step the controller may ever take -- makes the first step both the
        **largest** and the **only unchecked** one, because with no history there is
        nothing to difference and no truncation error can be estimated.  Its error
        then dominates everything after it.

        On the CPU that showed up as a global error of 1.3212e-01 at reltol 1e-3,
        1e-4, 1e-5 AND 1e-6 -- identical to five digits -- while the step count went
        from 24 to 195.  **Gate 9-1(c) measured the identical signature here**:
        4.2535e-3 across the same four decades, always at index 1, step count 53 to
        85.  Eight times the work, or 1.6x, for the same answer.

        Opening at `timestep * 1e-3` costs one cheap step and leaves the controller
        to grow from there, which it does geometrically, so the ramp is paid off
        within a handful of steps.
        """
        firststep = getattr(self.par, 'firststep', None)
        if firststep is None:
            return timestep * 1e-3
        if firststep <= 0:
            raise ValueError(
                "firststep must be positive, not %r; pass None to use the default "
                "ramp of timestep*1e-3" % (firststep,))
        ## Larger than dt_max would be capped on the very next step anyway, and
        ## asking for one is more likely a mistake than an intent.
        return min(float(firststep), float(timestep))

    def _newton_abstol(self, n):
        """Per-row Newton RESIDUAL floor: currents on node rows, volts on
        branch rows -- the same split Transient._newton builds (F6(b))."""
        import jax.numpy as jnp
        n_nodes = len(self.cir.nodes)
        return jnp.concatenate((
            self.par.iabstol * jnp.ones(n_nodes),
            self.par.vabstol * jnp.ones(n - n_nodes)))

    def _newton_xtol(self, n):
        """Per-row Newton UPDATE floor -- the transposed flavour: volts on
        node rows, amps on branch rows."""
        import jax.numpy as jnp
        n_nodes = len(self.cir.nodes)
        return jnp.concatenate((
            self.par.vabstol * jnp.ones(n_nodes),
            self.par.iabstol * jnp.ones(n - n_nodes)))

    def _lte_abstol(self, n):
        """Per-row absolute LTE tolerance, in the SOLUTION domain.

        Volts on node rows, amps on branch rows -- the same split
        `Transient._solve` makes.  A single scalar here is the residual-vs-solution
        flavour bug 0.3a fixed on the CPU, and 9(c) exists to not repeat it.
        """
        import jax.numpy as jnp
        n_nodes = len(self.cir.nodes)
        return jnp.concatenate((
            self.par.lte_vabstol * jnp.ones(n_nodes),
            self.par.lte_iabstol * jnp.ones(n - n_nodes)))

    def solve_batched(self, irefnode, override_params_tree, tend, timestep=1e-12, CHUNK_SIZE=5000, dt_min=1e-15, dt_max=None, uic=False):
        import jax

        self.cir.update_iparv()
        import jax.numpy as jnp
        import numpy as np

        ## An override for a class that is not in the vmap evaluation groups was
        ## SILENTLY IGNORED: `params_tree` is only consumed for classes with
        ## `eval_i_pure`/`eval_q_pure`, so `{'R': {'r': ...}}` produced N
        ## bit-identical lanes presented as N samples -- a parameter sweep that
        ## sweeps nothing, with no symptom (doc/transient_review_260820.md, F2;
        ## measured with r = 1 ohm vs 1 kohm).  Refuse loudly until the class is
        ## made batchable.
        batchable = {cls.__name__
                     for cls in self.toolkit.evaluation_groups(self.cir)}
        unbatchable = set(override_params_tree) - batchable
        if unbatchable:
            raise NotImplementedError(
                "override_params_tree names %s, which cannot be batched: only "
                "element classes providing eval_i_pure/eval_q_pure participate "
                "in vmapped evaluation (currently: %s). An override for any "
                "other class would be silently ignored and every lane would "
                "return the same waveform. Make the class batchable, or drop "
                "the override."
                % (', '.join(sorted(repr(n) for n in unbatchable)),
                   ', '.join(sorted(batchable)) or 'none'))

        if hasattr(self.cir, 'get_node_index'):
            irefnode = self.cir.get_node_index(irefnode)
            
        # Determine batch size from override_params_tree
        batch_size = 1
        for cls_name, params in override_params_tree.items():
            for p_name, p_val in params.items():
                batch_size = p_val.shape[0]
                break
            break
            
        n = self.cir.n
        if not uic:
            ## This used to be `pass`, under a comment saying a DC solve belonged
            ## here -- so every batched run silently started from zeros, whatever
            ## the circuit's bias point was.  That is the same defect `Transient`
            ## had, and worse: `Transient` at least attempted the solve before
            ## falling back.  Raising is the honest option until a batched DC
            ## exists, because a batched run seeded from zeros is not an
            ## approximation of the right answer, it is a different problem.
            raise NotImplementedError(
                "solve_batched has no DC operating-point solve, so it cannot start "
                "from a bias point. Pass uic=True to start from zeros deliberately. "
                "(The unbatched JAXTransient.solve() does solve a DC operating "
                "point; only the batched path lacks one.)")


        if dt_max is None:
            ## STAGE 9(g) -- `timestep`, not `tend/10`.  `solve` and `solve_batched`
            ## are two entry points of one class and disagreed about this by ~50x,
            ## so the same circuit at the same requested timestep was error-
            ## controlled to two different standards depending on which was called.
            ## `timestep` matches `solve` and matches `Transient.max_step`.
            dt_max = self._max_step(timestep)
            
        t_breaks_array = jnp.array(collect_breakpoints(self.cir, tend))
        
        # Setup batched initial state
        x0_batch = jnp.zeros((batch_size, n))
        x_hist = jnp.stack([x0_batch, x0_batch, x0_batch], axis=1) # (batch_size, 3, n)
        
        # We need a batched q0
        batched_q = jax.vmap(lambda x, p: self.cir.q(x, params_tree=p))(x0_batch, override_params_tree)
        q_hist = jnp.stack([batched_q, batched_q, batched_q], axis=1) # (batch_size, 3, n)
        
        iq_hist = jnp.zeros((batch_size, 3, n))
        h_hist = jnp.zeros((batch_size, 3))
        
        from pycircuit.circuit.elements import TLine
        tlines = []
        for instance_name, elem in self.cir.elements.items():
            if isinstance(elem, TLine):
                tlines.append((instance_name, elem))
        n_tlines = len(tlines)
        tline_params = jnp.zeros((n_tlines, 2))
        tline_indices = jnp.zeros((n_tlines, 6), dtype=jnp.int32)
        if n_tlines > 0:
            tline_history = jnp.zeros((batch_size, n_tlines, TLINE_HISTORY_DEPTH, 5))
            for i, (name, tline) in enumerate(tlines):
                nodemap = self.cir.elementnodemap[name]
                tline_indices = tline_indices.at[i].set(nodemap)
                tline_params = tline_params.at[i, 0].set(tline.iparv.TD)
                tline_params = tline_params.at[i, 1].set(tline.iparv.Z0)
                # Init with zero for now
        else:
            tline_history = jnp.zeros((batch_size, 0, TLINE_HISTORY_DEPTH, 5))
            
        tline_head = jnp.zeros(batch_size, dtype=jnp.int32)
        
        def run_chunk(s, p_tree):
            return outer_time_loop(s, self.cir, tend, CHUNK_SIZE, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method='gear', params_tree=p_tree, reltol=self.par.reltol, abstol=self._newton_abstol(n),
                                   xtol=self._newton_xtol(n),
                                   maxiter=self.par.maxiter, trtol=self.par.TRTOL,
                                   lte_reltol=self.par.reltol,
                                   lte_abstol=self._lte_abstol(n),
                                   max_dv=(jnp.inf if self.par.max_dv is None
                                           else float(self.par.max_dv)))
            
        # JIT the vmapped run_chunk
        batched_run_chunk = jax.jit(jax.vmap(run_chunk))
        
        ## Per-lane accumulators, each seeded with the lane's t=0 point (F12's
        ## guarantee on this path; see the collection note below for why the
        ## lanes are not kept rectangular).
        x0_np = np.array(x0_batch)
        results_list = [[x0_np[b:b+1, :]] for b in range(batch_size)]
        times_list = [[np.zeros(1)] for b in range(batch_size)]
        
        current_t = jnp.zeros(batch_size)
        ## Same ramp as `solve` -- see `_opening_step`.
        current_dt = jnp.full(batch_size, self._opening_step(timestep))
        b_sig_max = jnp.zeros(batch_size)
        ## dtype follows `jnp.array(0)` exactly, as the unbatched path uses:
        ## `lax.cond` requires both branches to agree, and a mismatched
        ## int32/int64 here fails at trace time rather than at the call.
        b_n_rejected = jnp.full(batch_size, 0)
        b_n_nonconverged = jnp.full(batch_size, 0)
        b_n_forced_nc = jnp.full(batch_size, 0)
        b_n_forced_lte = jnp.full(batch_size, 0)
        b_force_first = jnp.full(batch_size, False)
        
        while np.any(current_t < tend):
            res_buf = jnp.zeros((batch_size, CHUNK_SIZE, n))
            time_buf = jnp.zeros((batch_size, CHUNK_SIZE))
            
            ## Every field of the batched state must carry a leading batch axis:
            ## `batched_run_chunk` is a `jax.vmap` over axis 0, and a field left at
            ## its scalar default fails with "rank should be at least 1".  The
            ## counters added by 9(c)/9(e) -- sig_max, n_rejected, n_nonconverged,
            ## n_forced_nonconverged -- default to scalars for the unbatched path
            ## and so must be widened here explicitly.  They are also RUNNING
            ## TOTALS, carried across chunks for the same reason `solve` carries
            ## them: rebuilding them per chunk resets the sigglobal reference.
            state = TransientState(
                t=current_t, dt=current_dt, step_idx=jnp.zeros(batch_size, dtype=jnp.int32),
                x_history=x_hist, q_history=q_hist, iq_history=iq_hist, h_history=h_hist,
                results_buffer=res_buf, time_buffer=time_buf,
                tline_history=tline_history, tline_head=tline_head,
                sig_max=b_sig_max, n_rejected=b_n_rejected,
                n_nonconverged=b_n_nonconverged,
                n_forced_nonconverged=b_n_forced_nc,
                n_forced_lte=b_n_forced_lte,
                force_first_order=b_force_first
            )
            
            final_state = batched_run_chunk(state, override_params_tree)

            ## PER-LANE COLLECTION, NO PADDING.  The old code trimmed every lane
            ## to the batch's rectangular [0, max_steps) window and FILLED
            ## SHORTER LANES FORWARD -- duplicating both the state AND the
            ## timestamp, so shorter lanes' results carried repeated abscissae
            ## (which break interpolation downstream), and a lane that had
            ## already reached tend in a previous chunk had b_len == 0, making
            ## the fill source `x_chunk[b, -1:0, :]` an EMPTY slice: numpy
            ## raised "could not broadcast (0,n) into (max_steps,n)" and any
            ## heterogeneous sweep spanning more than one chunk crashed.
            ##
            ## Nothing needs the lanes to be rectangular: chunk-to-chunk state
            ## flows through `final_state`, and the method returns a separate
            ## Result per lane, each with its own time base.  So each lane
            ## simply keeps its own valid slice and the padding -- the only
            ## code that could crash or corrupt -- is deleted rather than
            ## guarded (doc/transient_review_260820.md, F1(c); measured with
            ## c = 1e-9 vs 1e-6 F, CHUNK_SIZE=10).
            valid_steps = np.array(final_state.step_idx)
            if int(np.max(valid_steps)) == 0:
                break

            res_buf = np.array(final_state.results_buffer)
            time_buf = np.array(final_state.time_buffer)
            for b in range(batch_size):
                b_len = int(valid_steps[b])
                if b_len > 0:
                    results_list[b].append(res_buf[b, :b_len, :])
                    times_list[b].append(time_buf[b, :b_len])

            current_t = final_state.t
            current_dt = final_state.dt
            b_sig_max = final_state.sig_max
            b_n_rejected = final_state.n_rejected
            b_n_nonconverged = final_state.n_nonconverged
            b_n_forced_nc = final_state.n_forced_nonconverged
            b_n_forced_lte = final_state.n_forced_lte
            b_force_first = final_state.force_first_order
            x_hist = final_state.x_history
            q_hist = final_state.q_history
            iq_hist = final_state.iq_history
            h_hist = final_state.h_history
            tline_history = final_state.tline_history
            tline_head = final_state.tline_head
            
        ## Each lane concatenates its own slices -- lengths may differ, which
        ## is the point.  The t=0 seed row (results_list[b][0]) is what keeps
        ## F12's include-the-initial-point guarantee on this path.
        res = []
        for b in range(batch_size):
            xb = np.concatenate(results_list[b], axis=0)
            tb = np.concatenate(times_list[b])
            res.append(Result(self.cir, x=xb.T, xdot=None,
                              sweep_values=tb,
                              sweep_label='time',
                              sweep_unit='s'))
        return res

    def solve(self, refnode=0, tend=1e-3, x0=None, timestep=1e-6, CHUNK_SIZE=5000, **kwargs):
        n = self.cir.n
        irefnode = self.cir.get_node_index(refnode)
    
        if x0 is None:
            if kwargs.get('uic', False):
                x0 = np.zeros(n)
                # Load node initial conditions if they exist (cir.nodes is a list).
                for node in self.cir.nodes:
                    if hasattr(node, 'ic'):
                        idx = self.cir.get_node_index(node)
                        x0[idx] = node.ic
            else:
                # Run a fast DC operating point if none provided
                from pycircuit.circuit.dcanalysis import DC
                x0 = DC(self.cir).solve().x
            
        dt_min = kwargs.get('minstep', 1e-18)
        dt_max = self._max_step(timestep)
    
        t_breaks_array = jnp.array(collect_breakpoints(self.cir, tend))
    
        # Extract TLines
        from pycircuit.circuit.elements import TLine
        tlines = []
        for instance_name, elem in self.cir.elements.items():
            if isinstance(elem, TLine):
                tlines.append((instance_name, elem))
    
        tline_params_list = []
        tline_indices_list = []
        for name, tline in tlines:
            tline_params_list.append([tline.iparv.TD, tline.iparv.Z0])
            tline_indices_list.append(self.cir.elementnodemap[name])
        
        tline_params = jnp.array(tline_params_list) if tline_params_list else jnp.zeros((0, 2))
        tline_indices = jnp.array(tline_indices_list, dtype=jnp.int32) if tline_indices_list else jnp.zeros((0, 6), dtype=jnp.int32)
    
        # We will use a fixed buffer size of 10000 steps for TLines
        # This covers a 10ns simulation with 1ps steps.
        # Shape: (N_tlines, 10000, 5) => (t, v1, v2, i1, i2)
        n_tlines = len(tlines)
        if n_tlines > 0:
            tline_history = jnp.zeros((n_tlines, TLINE_HISTORY_DEPTH, 5))
            # Initialize with DC values for all history!
            for i, (name, tline) in enumerate(tlines):
                # evaluate tline DC state from x0
                nodemap = self.cir.elementnodemap[name]
                subx = x0[nodemap] if len(nodemap) > 0 else jnp.array([])
                v1 = subx[0] - subx[1]
                v2 = subx[2] - subx[3]
                i1 = subx[4]
                i2 = subx[5]
                init_val = jnp.array([-1.0, v1, v2, i1, i2]) # t=-1.0 to avoid t=0 collisions
                tline_history = tline_history.at[i].set(init_val)
                # Ensure t=0 is at head 0
                tline_history = tline_history.at[i, 0].set(jnp.array([0.0, v1, v2, i1, i2]))
        else:
            tline_history = jnp.zeros((0, TLINE_HISTORY_DEPTH, 5))
        
        tline_head = jnp.array(0, dtype=jnp.int32)
    
        # Setup initial state
        x_hist = jnp.array([x0, x0, x0])
        q0 = self.cir.q(x0)
        q_hist = jnp.array([q0, q0, q0])
        iq_hist = jnp.zeros((3, n))
        h_hist = jnp.zeros(3)
    
        # We need a JIT-able wrapper that burns the chunk_size and tend into static parameters if needed,
        # but tend is a runtime parameter.
        @jax.jit
        def run_chunk(s):
            return outer_time_loop(s, self.cir, tend, CHUNK_SIZE, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method='gear', reltol=self.par.reltol, abstol=self._newton_abstol(n),
                                   xtol=self._newton_xtol(n),
                                   maxiter=self.par.maxiter, trtol=self.par.TRTOL,
                                   lte_reltol=self.par.reltol,
                                   lte_abstol=self._lte_abstol(n),
                                   max_dv=(jnp.inf if self.par.max_dv is None
                                           else float(self.par.max_dv)))
        
        results_list = [np.array([x0])]
        times_list = [np.array([0.0])]
    
        current_t = 0.0
        current_dt = self._opening_step(timestep)
        ## Accepted steps so far, tracked only to detect the TLine ring wrapping.
        ## The check cannot live inside the compiled loop -- it cannot raise there --
        ## so it is done here, per chunk, which is soon enough: the wrap corrupts the
        ## interpolation from the step after it happens, and a chunk boundary is at
        ## most CHUNK_SIZE steps later.
        total_steps = 0
        sig_max = jnp.array(0.0)
        n_rejected = jnp.array(0)
        n_nonconverged = jnp.array(0)
        n_forced_nc = jnp.array(0)
        n_forced_lte = jnp.array(0)
        force_first = jnp.array(False)

        while current_t < tend:
            res_buf = jnp.zeros((CHUNK_SIZE, n))
            time_buf = jnp.zeros(CHUNK_SIZE)
        
            ## `sig_max` and `n_rejected` are RUNNING TOTALS and must cross the chunk
            ## boundary.  Rebuilding the state without them reset the `sigglobal`
            ## reference to zero every CHUNK_SIZE steps -- so a long run silently
            ## reverted to a `pointlocal`-like tolerance at each boundary -- and threw
            ## the rejection count away.  Same shape as the `_dt_last2` reset the CPU
            ## side had: a per-run quantity re-seeded by a per-call constructor.
            state = TransientState(
                t=jnp.array(current_t), dt=jnp.array(current_dt), step_idx=jnp.array(0),
                x_history=x_hist, q_history=q_hist, iq_history=iq_hist, h_history=h_hist,
                results_buffer=res_buf, time_buffer=time_buf,
                tline_history=tline_history, tline_head=tline_head,
                sig_max=sig_max, n_rejected=n_rejected,
                n_nonconverged=n_nonconverged,
                n_forced_nonconverged=n_forced_nc,
                n_forced_lte=n_forced_lte,
                force_first_order=force_first
            )
        
            final_state = run_chunk(state)
        
            # Extract valid steps
            valid_idx = int(final_state.step_idx)
        
            # If no steps were taken (e.g. tend was already reached or error), break
            if valid_idx == 0:
                break
            
            total_steps += valid_idx
            if n_tlines > 0 and total_steps >= TLINE_HISTORY_DEPTH:
                raise RuntimeError(
                    "TLine delay-line history overflowed: %d accepted steps against a "
                    "ring buffer of %d (TLINE_HISTORY_DEPTH). Past this point the "
                    "buffer wraps and the delay interpolation silently reads entries "
                    "from a previous lap, so the waveform would be wrong with no "
                    "other symptom. Shorten tend, use a larger timestep, or raise "
                    "TLINE_HISTORY_DEPTH."
                    % (total_steps, TLINE_HISTORY_DEPTH))

            x_chunk = np.array(final_state.results_buffer[:valid_idx])
            t_chunk = np.array(final_state.time_buffer[:valid_idx])

            results_list.append(x_chunk)
            times_list.append(t_chunk)
        
            # Update loop variables for next chunk
            current_t = float(final_state.t)
            current_dt = float(final_state.dt)
            x_hist = final_state.x_history
            q_hist = final_state.q_history
            iq_hist = final_state.iq_history
            h_hist = final_state.h_history
            tline_history = final_state.tline_history
            tline_head = final_state.tline_head
            sig_max = final_state.sig_max
            n_rejected = final_state.n_rejected
            n_nonconverged = final_state.n_nonconverged
            n_forced_nc = final_state.n_forced_nonconverged
            n_forced_lte = final_state.n_forced_lte
            force_first = final_state.force_first_order

            ## STAGE 9(e) -- CHECKED PER CHUNK, NOT AT THE END, BECAUSE THE END MAY
            ## NEVER ARRIVE.  Rejecting a non-converged step shrinks `dt`; a circuit
            ## that cannot converge at any step size drives `dt` to `dt_min` and then
            ## advances by `dt_min` forever.  Measured while writing this: a linear RC
            ## at maxiter=1 needs ~5e12 steps to reach tend that way, so warning
            ## "after the run" is a warning that never prints.  That is the trade the
            ## plan flags for 9(d) -- turning a silent wrong answer into a hang -- and
            ## it is avoided by raising here, where the loop is bounded by CHUNK_SIZE.
            if int(n_forced_nc) > 0:
                raise NoConvergenceError(
                    "jaxtransient: the Newton solve failed to converge at t=%g s "
                    "even at the minimum step dt_min=%g, so %d step(s) would have "
                    "been committed without solving the circuit equations. The "
                    "transient is NOT complete. Raise maxiter (currently %d), "
                    "loosen reltol (currently %g), or lower dt_min. Previously "
                    "these points were committed silently -- measured at 11.6x the "
                    "error on a linear RC, reported as nothing."
                    % (float(final_state.t), dt_min, int(n_forced_nc),
                       int(self.par.maxiter), float(self.par.reltol)))
        
        ## STAGE 9 -- what the run actually did, as opposed to what it returned.
        ## The CPU's `TransientStatistics` exists for the same reason; this is the
        ## subset a traced loop can count.  `rejected_steps` is what makes gate
        ## 9-1(b) -- "a step is actually rejected" -- expressible on this backend.
        self.statistics = JAXTransientStatistics()
        self.statistics.accepted_steps = int(total_steps)
        self.statistics.rejected_steps = int(n_rejected)
        self.statistics.signal_max = float(sig_max)
        self.statistics.nonconverged_steps = int(n_nonconverged)
        self.statistics.forced_nonconverged_steps = int(n_forced_nc)
        self.statistics.forced_lte_steps = int(n_forced_lte)

        ## F19(c): mirror the CPU force-accept warning -- an unbounded
        ## accepted truncation error must not be invisible.
        if int(n_forced_lte) > 0:
            warnings.warn(
                'jaxtransient: %d step(s) were accepted at dt_min with the '
                'local truncation error still above tolerance. The accepted '
                'error there is unbounded -- treat the waveform near those '
                'times with suspicion (the CPU path calls these '
                'force_accepts).' % int(n_forced_lte), RuntimeWarning)

        ## Steps REJECTED for non-convergence are recoverable -- the controller
        ## shrank and retried -- but a run full of them is one to look at, so it is
        ## surfaced rather than left in the statistics object alone.  The
        ## unrecoverable case raises inside the chunk loop above.
        if int(n_nonconverged) > 0:
            warnings.warn(
                'jaxtransient: the Newton solve failed to converge on %d step '
                'attempt(s); each was rejected and retried at a smaller dt, so the '
                'result is still error-controlled, but the circuit is hard for the '
                'solver at this tolerance.' % int(n_nonconverged),
                RuntimeWarning)

        # Concatenate
        all_results = np.vstack(results_list)
        all_times = np.concatenate(times_list)
    
        # Build Result
        from pycircuit.circuit.analysis import CircuitResult
        res = CircuitResult(self.cir, x=all_results.T, xdot=None,
                            sweep_values=all_times,
                            sweep_label='time',
                            sweep_unit='s')
        ## Reachable from the result, matching the CPU paths (F13): a caller
        ## who kept only the waveform can still ask what produced it.
        res.statistics = self.statistics
        return res
