from typing import NamedTuple, Any
import numpy as np
import jax
import warnings

import jax.numpy as jnp

from pycircuit.circuit._lte_kernels import (bdf2_companion, euler_companion,
                                            second_divided_difference,
                                            euler_lte, gear2_lte)
from pycircuit.circuit.analysis import Analysis
from pycircuit.circuit.nrsolver import NoConvergenceError
from pycircuit.utilities.param import Parameter
from pycircuit.circuit.analysis import CircuitResult as Result
from pycircuit.circuit.analysis import gnd

## Depth of the per-TLine delay-line ring buffer, in accepted steps.  It was the
## bare literal `10000` in seven places, three of them inside modulo arithmetic,
## which made it impossible to see that the buffer is a RING: past this many
## accepted steps `tline_head` wraps and `interp_tlines` starts interpolating
## against entries from a previous lap.  `cond_fun` stops at the end of the buffer
## and returns whatever is there, so the result is a plausible wrong waveform with
## no error and no warning.  `JAXTransient.solve` now checks for the wrap and
## raises; sizing the buffer from the run instead of fixing it belongs to stage 9.
TLINE_HISTORY_DEPTH = 10000

## The CPU's within-time-point excursion clamps (stepcontroller
## MIN_SHRINK_RATIO / MAX_GROWTH_RATIO), re-exported as plain floats for the
## traced Fang loop.
from pycircuit.circuit.stepcontroller import (MIN_SHRINK_RATIO as MIN_SHRINK,
                                              MAX_GROWTH_RATIO as MAX_GROWTH)


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
    ## P7: the per-row running signal maximum, relref's 'alllocal' and the
    ## unit-split 'sigglobal' reference.  Shape (n,), updated on ACCEPT only
    ## (the same hygiene as sig_max -- a rejected iterate must not loosen
    ## later tolerances).
    ref_running: Any = 0.0

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
## not move.
##
## THE TRAPEZOIDAL BRANCH WAS DELETED (review hygiene): no production path
## could reach it (at the time, `eval_method` was hardcoded 'gear' at both
## call sites; P6 has since exposed the gear/euler choice as the
## `integrator` Parameter), and its LTE formula was the uniform-grid one -- wrong the
## moment the step changed, had anyone ever wired it up.  Dead-but-plausible
## solver branches are exactly how the 3/4-optimism defect survived twice.
## The CPU's trapezoidal integrator, with the correct variable-step
## estimator, lives in integrator.py.

def backward_euler_step(q_curr, C_curr, q_prev, dt):
    return euler_companion(q_curr, C_curr, q_prev, dt)


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


def compute_integration(q_curr, C_curr, state: TransientState, method='gear',
                        first_order=None):

    q_prev = state.q_history[0]
    dt = state.dt

    def do_euler():
        return backward_euler_step(q_curr, C_curr, q_prev, dt)

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
    raise ValueError("eval_method must be 'euler' or 'gear', not %r" % (method,))

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

def _tline_emfs(t_target, params, history, head):
    """One line's reflected-wave EMFs (e1, e2) AND their time derivatives.

    Extracted from newton_inner_loop's closures (TLine-under-coupled/PCNR
    port): the delay-line history lives in TransientState and is updated by
    the shared accept machinery, so the only thing that kept the other
    assemblies TLine-blind was this code being trapped in one function.  The
    derivatives are the linear interpolant's segment slope -- free from the
    same lookup -- and are what Fang's ``p = df/dh`` needs from a source
    whose value depends on the step size through ``t - TD``.
    """
    def cond_fun(idx):
        curr_t = history[(head - idx) % TLINE_HISTORY_DEPTH, 0]
        return jnp.logical_and(curr_t > t_target, idx < TLINE_HISTORY_DEPTH - 1)

    idx1 = jax.lax.while_loop(cond_fun, lambda idx: idx + 1, 0)
    idx0 = jnp.maximum(0, idx1 - 1)
    ## P16: the third stencil point, one OLDER than the bracket -- the CPU's
    ## quadratic (its `i-1` point), so a first-order interpolant no longer
    ## feeds a second-order method error the method never made.
    idx2 = jnp.minimum(idx1 + 1, TLINE_HISTORY_DEPTH - 1)

    val1 = history[(head - idx1) % TLINE_HISTORY_DEPTH]
    val0 = history[(head - idx0) % TLINE_HISTORY_DEPTH]
    val2 = history[(head - idx2) % TLINE_HISTORY_DEPTH]
    t1, t0, t2 = val1[0], val0[0], val2[0]

    denom_ok = t0 != t1
    inv_dt = jnp.where(denom_ok, 1.0 / jnp.where(denom_ok, t0 - t1, 1.0), 0.0)
    frac = (t_target - t1) * inv_dt
    lin_val = val1 + frac * (val0 - val1)
    lin_slope = (val0 - val1) * inv_dt

    ## Quadratic Lagrange over (t2, t1, t0), value and d/dt, guarded so a
    ## degenerate stencil (ring not yet filled: t2 == t1) falls back to the
    ## linear form.  MONOTONE-LIMITED per channel exactly like the CPU's
    ## `_interpolate_history`: a quadratic value outside the bracket's own
    ## [min, max] is evidence the stencil spans a kink -- the phantom-EMF
    ## defect the CPU limiter was measured against (reflected EMF read
    ## 1.009 against samples bounded by 1.000) -- and that channel keeps
    ## the linear value AND the linear slope, so du/dt stays d/dt of u.
    quad_ok = jnp.logical_and(denom_ok, t2 != t1)
    d20 = jnp.where(quad_ok, t2 - t0, 1.0)
    d21 = jnp.where(quad_ok, t2 - t1, 1.0)
    d10 = jnp.where(denom_ok, t1 - t0, 1.0)
    L2 = (t_target - t1) * (t_target - t0) / (d21 * d20)
    L1 = (t_target - t2) * (t_target - t0) / (-d21 * d10)
    L0 = (t_target - t2) * (t_target - t1) / (d20 * d10)
    q_val = val2 * L2 + val1 * L1 + val0 * L0
    dL2 = (2.0 * t_target - t1 - t0) / (d21 * d20)
    dL1 = (2.0 * t_target - t2 - t0) / (-d21 * d10)
    dL0 = (2.0 * t_target - t2 - t1) / (d20 * d10)
    q_slope = val2 * dL2 + val1 * dL1 + val0 * dL0

    lo = jnp.minimum(val1, val0)
    hi = jnp.maximum(val1, val0)
    in_env = jnp.logical_and(q_val >= lo, q_val <= hi)
    use_quad = jnp.logical_and(quad_ok, in_env)
    interp_val = jnp.where(use_quad, q_val, lin_val)
    slope = jnp.where(use_quad, q_slope, lin_slope)

    Z0 = params[1]
    e1 = interp_val[2] + Z0 * interp_val[4]
    e2 = interp_val[1] + Z0 * interp_val[3]
    de1 = slope[2] + Z0 * slope[4]
    de2 = slope[1] + Z0 * slope[3]
    return e1, e2, de1, de2


def tline_stamp_correction(n, tlines_meta):
    """The constant (n, n) correction (transient stamp - DC stamp), per line.

    THE STANDARD-PATH DEFECT THIS FIXES (found while porting TLine under the
    coupled/PCNR paths, present since the JAX TLine landed): ``TLine.G``
    selects its stamp on ``len(self.history) == 0`` -- Python state that only
    the CPU's accept_step writes.  On this backend the traced ring buffer
    replaced accept_step, the Python history stays empty forever, and every
    assembly therefore stamped the DC form (the line as an ideal short:
    v1 - v2 = 0, i1 + i2 = 0) while the injected EMFs assumed the transient
    form (v - Z0*i = e).  Measured on a matched 1 V line: |v(far end)| =
    24.5 V, delay 0 -- and NO e2e JAX TLine test existed to see it.

    Both stamps are constants, so their difference is a constant matrix built
    once at trace time; adding it to the assembled J (and dG @ x to the
    residual) converts the DC-stamped rows to the transient form everywhere.

    ``tlines_meta`` is [(nodemap6, Z0), ...].
    """
    dG = np.zeros((n, n))
    for rows, Z0 in tlines_meta:
        r4, r5 = rows[4], rows[5]
        p1, m1, p2, m2 = rows[0], rows[1], rows[2], rows[3]
        ## DC row 4: v1 - v2 = 0 ; transient row 4: v1 - Z0*i1 = e1
        dG[r4, p2] += 1.0
        dG[r4, m2] -= 1.0
        dG[r4, r4] += -Z0
        ## DC row 5: i1 + i2 = 0 ; transient row 5: v2 - Z0*i2 = e2
        dG[r5, p2] += 1.0
        dG[r5, m2] -= 1.0
        dG[r5, r4] += -1.0
        dG[r5, r5] += -Z0 - 1.0
    return jnp.asarray(dG)


def apply_tline_sources(I_u, t_curr, tline_params, tline_indices,
                        tline_history, tline_head):
    """Add every line's reflected-wave EMFs into the source vector."""
    if tline_params.shape[0] == 0:
        return I_u
    t_targets = t_curr - tline_params[:, 0]
    e1s, e2s, _d1, _d2 = jax.vmap(_tline_emfs, in_axes=(0, 0, 0, None))(
        t_targets, tline_params, tline_history, tline_head)
    I_u = I_u.at[tline_indices[:, 4]].add(-e1s)
    I_u = I_u.at[tline_indices[:, 5]].add(-e2s)
    return I_u


def tline_source_dudt(n, t_curr, tline_params, tline_indices,
                      tline_history, tline_head):
    """d/dt of the TLine source contribution -- the piece of Fang's ``p``
    a delay line owns."""
    out = jnp.zeros(n)
    if tline_params.shape[0] == 0:
        return out
    t_targets = t_curr - tline_params[:, 0]
    _e1, _e2, d1s, d2s = jax.vmap(_tline_emfs, in_axes=(0, 0, 0, None))(
        t_targets, tline_params, tline_history, tline_head)
    out = out.at[tline_indices[:, 4]].add(-d1s)
    out = out.at[tline_indices[:, 5]].add(-d2s)
    return out


def newton_inner_loop(state: TransientState, circuit, irefnode, tline_params, tline_indices, eval_method='euler', reltol=1e-3, abstol=1e-6, xtol=1e-12, maxiter=50, params_tree=None, max_dv=jnp.inf, tline_dG=None, analysis='tran', provided_function=None):
    x_init = extrapolate_predictor(state)

    def apply_tlines(I_u, t_curr):
        ## Extracted to module level so the fang and pcnr assemblies share it.
        return apply_tline_sources(I_u, t_curr, tline_params, tline_indices,
                                   state.tline_history, state.tline_head)


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
        I_u = circuit.u(state.t + state.dt, analysis=analysis, params_tree=params_tree)
        if provided_function is not None:
            ## P11: the F4 contract -- an extra source, folded where u is.
            I_u = I_u + provided_function(state.t + state.dt)
        I_u = apply_tlines(I_u, state.t + state.dt)

        i_C, G_C_eq = compute_integration(q_C, C_C, state, method=eval_method)

        F = I_G + i_C + I_u
        J = G_G + G_C_eq
        if tline_dG is not None:
            ## Stamp correction: see tline_stamp_correction.
            F = F + tline_dG @ x
            J = J + tline_dG

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
        ## THE REDUCED ROWS, not all of them (P11's port found this): the
        ## reference row is not an equation of the solved system -- its
        ## residual is whatever KCL imbalance the source terms carry (an
        ## unbalanced provided_function put the full injection there), and
        ## scoring it can NEVER be satisfied by any x, which livelocked the
        ## run at maxiter on every step, silently force-accepting at the dt
        ## floor.  The CPU's Newton has always tested the reduced system.
        conv_f = jnp.all(jnp.delete(
            jnp.abs(F) < reltol * I_scale
            + jnp.broadcast_to(abstol, F.shape), irefnode))
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
                    method='gear', trtol=7.0, lte_rel=1e-4, lte_abstol=1e-12,
                    first_order=None, relref='sigglobal', n_nodes=None):
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
        ## accuracy, on the default eval_method of both entry points.
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
    else:
        raise ValueError(
            "method must be 'euler' or 'gear', not %r (the trapezoidal "
            "branch was deleted -- its formula was uniform-grid-only and no "
            "production path could reach it)" % (method,))

    # lte = J^-1 Eg in the reduced (reference-node-removed) space.
    J_r = jnp.delete(jnp.delete(J, irefnode, axis=0), irefnode, axis=1)
    Eg_r = jnp.delete(Eg, irefnode)
    lte_r = jnp.linalg.solve(J_r, Eg_r)
    lte = jnp.insert(lte_r, irefnode, 0.0)

    ## P7: the CPU's relref modes, selected by a trace-static string.  The
    ## reference is a per-row VECTOR now; 'sigglobal' collapses each UNIT
    ## GROUP (node voltages vs branch currents) to its own maximum -- the
    ## old scalar max over everything mixed volts and amps, which the CPU
    ## comment warns silently disables node error control on circuits with
    ## large branch currents.  `local` includes the current iterate, exactly
    ## as the CPU's `_reference` does.
    local = jnp.maximum(jnp.abs(x_curr), jnp.abs(state.x_history[0]))
    if relref == 'pointlocal':
        ref = local
    else:
        run = jnp.maximum(state.ref_running, local)
        if relref == 'alllocal':
            ref = run
        else:  ## sigglobal
            n_all = x_curr.shape[0]
            if n_nodes is None or n_nodes <= 0 or n_nodes >= n_all:
                ref = jnp.full(n_all, jnp.max(run))
            else:
                ref = jnp.concatenate([
                    jnp.full(n_nodes, jnp.max(run[:n_nodes])),
                    jnp.full(n_all - n_nodes, jnp.max(run[n_nodes:]))])
    ## `lte_abstol` is a per-row VECTOR (volts on node rows, amps on branch rows),
    ## built by `JAXTransient._lte_abstol`.  A scalar here applies one physical
    ## kind of tolerance to rows of another -- 0.3a's residual-vs-solution defect,
    ## which 9(c) exists to avoid repeating.  A scalar still broadcasts, so a
    ## direct caller that passes one gets the old behaviour rather than an error.
    etol = trtol * (lte_rel * ref + lte_abstol)
    return jnp.max(jnp.abs(lte) / etol), order_p1


def collect_breakpoints(cir, tend, minbreak=1e-14):
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
    ## TLINE WAVEFRONT ARRIVALS: a source corner reaching a delay line
    ## re-emerges at the far end TD later -- a from-zero kink in an ALGEBRAIC
    ## variable that no element reports.  Registering {corner + TD, corner +
    ## 2*TD} makes those steps breakpoint-landings, which the coupled path's
    ## post-breakpoint grace (see fang_inner_loop) can then absorb -- the two
    ## mechanisms only work TOGETHER: arrivals-as-breakpoints alone was tried
    ## first and falsified, because a registered kink without the grace still
    ## livelocks on the h-independent relative LTE.  Deeper bounce ancestry
    ## is truncated, as every SPICE truncates it.
    tline_tds = sorted({float(elem.iparv.TD)
                        for _n, elem in cir.elements.items()
                        if type(elem).__name__ == 'TLine'})
    if tline_tds and breakpoints:
        base = sorted(set(breakpoints))
        arrivals = []
        for td_ in tline_tds:
            arrivals += [b + td_ for b in base if b + td_ < tend]
            arrivals += [b + 2.0 * td_ for b in base if b + 2.0 * td_ < tend]
        breakpoints += arrivals
    ## P14: MERGE breakpoints closer than `minbreak`, the CPU's spacing rule.
    ## Without it two sources with edges 1e-16 apart produce two breakpoints
    ## and a degenerate truncated step between them.  Clusters keep their
    ## FIRST member -- except at `tend`, which must survive verbatim so the
    ## run lands on it exactly (F3): a breakpoint within `minbreak` of `tend`
    ## is dropped in its favour.
    merged = []
    for b in sorted(set(breakpoints)):
        tol = minbreak * max(abs(merged[-1]), 1.0) if merged else 0.0
        if merged and b - merged[-1] <= tol:
            continue
        merged.append(b)
    if merged and merged[-1] != tend \
            and tend - merged[-1] <= minbreak * max(abs(tend), 1.0):
        merged.pop()
    if not merged or merged[-1] != tend:
        merged.append(tend)
    return merged


def calculate_next_dt(dt, error_ratio, dt_min, dt_max, t_breaks_array, current_t, order_p1=2.0, eta=jnp.inf):
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
    ## P8: eq (16)'s damper, the CPU's `_damp` -- inert at the default
    ## eta=inf, a clamp on the step CHANGE when the band sets one.
    factor = jnp.clip(factor, jnp.maximum(MIN_SHRINK_RATIO, 1.0 - eta),
                      jnp.minimum(MAX_GROWTH_RATIO, 1.0 + eta))
    dt_new = dt * factor
    dt_new = jnp.clip(dt_new, dt_min, dt_max)

    future_breaks = jnp.where(t_breaks_array > current_t + 1e-12, t_breaks_array, jnp.inf)
    next_break = jnp.min(future_breaks)

    dt_new = jnp.where(current_t + dt_new > next_break, next_break - current_t, dt_new)
    dt_new = jnp.maximum(dt_new, dt_min)
    return dt_new


## ---------------------------------------------------------------------------
## Fang's coupled (x, h) time step -- the 'approx' branch, traced (P19).
##
## A faithful port of Transient._fang_timestep_inner's DEFAULT path, scoped the
## way doc/backend_parity_260821.md P19 records: sec. 3.4's error-ratio step
## correction with the eq (18) solution update, hold_h for imposed step sizes,
## the within-point excursion clamps and the thwarted-shrink saturation test.
## NOT ported (CPU-only until separately justified): the 'bordered' eq (12)
## branch, PCNR-inside-Fang, grid_locked (no fixed_timestep on this backend
## yet), and cir.limit (this element set has no limiter; nonlinear robustness
## on this backend is PCNR's job, next stage).  The LTE degree follows the
## effective order (F19): both degrees are computed and selected, which keeps
## every shape static under the trace.
## ---------------------------------------------------------------------------

class FangState(NamedTuple):
    x: Any
    h: Any
    iters: Any
    done: Any
    converged: Any
    ## PCNR-inside-Fang: the per-junction limited unknowns.  Zero-length when
    ## PCNR is off, so the state shape stays static either way.
    v_lim: Any = None


def _fang_lte(x_new, state, h, first_order):
    """Fang eq (6) deviation for degree 2 and 1, selected by effective order.

    ``x_history``/``h_history`` are the accepted rings, most recent first;
    zero step entries (run start) are guarded the way ``ywr_error_ratio``
    guards ``dt_prev``.
    """
    from pycircuit.circuit._lte_kernels import solution_lte
    h0 = jnp.where(state.h_history[0] == 0.0, h, state.h_history[0])
    h1 = jnp.where(state.h_history[1] == 0.0, h0, state.h_history[1])
    lte2 = solution_lte(x_new, [state.x_history[0], state.x_history[1],
                                state.x_history[2]], [h0, h1], h)
    lte1 = solution_lte(x_new, [state.x_history[0], state.x_history[1]],
                        [h0], h)
    return jnp.where(first_order, lte1, lte2), h0, h1


def fang_inner_loop(state: TransientState, circuit, irefnode,
                    hold_h, tline_params, tline_indices,
                    eval_method='gear', reltol=1e-4,
                    xtol=1e-12, lte_abstol=1e-12, trtol=7.0, maxiter=100,
                    params_tree=None, gamma_min=0.7, gamma_max=3.0, eta=0.15,
                    dt_min=1e-18, dt_max=jnp.inf, pcnr_meta=None,
                    pcnr_VT=0.025, tline_dG=None, analysis='tran',
                    provided_function=None):
    from pycircuit.circuit._lte_kernels import (step_for_error_ratio,
                                                euler_companion_dh,
                                                bdf2_companion_dh)

    t_prev = state.t
    h_entry = state.dt
    h_floor = jnp.maximum(dt_min, h_entry * MIN_SHRINK)
    h_ceil = jnp.minimum(dt_max, h_entry * MAX_GROWTH)
    target = (gamma_min * gamma_max) ** 0.5
    first_order = effective_first_order(state)
    ## The band cannot be evaluated before two accepted points exist.  The
    ## post-breakpoint grace on delay-line circuits arrives through this same
    ## test: the accept path EMPTIES h_history on a breakpoint landing when
    ## the circuit has TLines (see do_accept), so the next step reads as
    ## history-free here.  An explicit `or force_first_order` term was tried
    ## and reverted -- it duplicated the reset's effect on TLine circuits and
    ## extended a band-blind step to every breakpoint on circuits that never
    ## needed it, the exact accuracy cost the CPU measured at 9.8e-4 ->
    ## 5.4e-3 median on the pulsed RC when its reset ran ungated.
    no_hist = state.h_history[0] == 0.0

    def assemble(x, h, v_lim):
        st_h = state._replace(dt=h)
        I_G = circuit.i(x, params_tree=params_tree)
        G_G = circuit.G(x, params_tree=params_tree)
        q_C = circuit.q(x, params_tree=params_tree)
        C_C = circuit.C(x, params_tree=params_tree)
        I_u = circuit.u(t_prev + h, analysis=analysis, params_tree=params_tree)
        if provided_function is not None:
            I_u = I_u + provided_function(t_prev + h)
        I_u = apply_tline_sources(I_u, t_prev + h, tline_params,
                                  tline_indices, state.tline_history,
                                  state.tline_head)
        i_C, G_eq = compute_integration(q_C, C_C, st_h, method=eval_method,
                                        first_order=first_order)
        F, J = I_G + i_C + I_u, G_G + G_eq
        if tline_dG is not None:
            F = F + tline_dG @ x
            J = J + tline_dG
        if pcnr_meta is None:
            return F, J, q_C, None
        ## PCNR-INSIDE-FANG: the CPU's design note is the whole argument --
        ## the Schur-reduced (f_eff, J_eff) IS an n-sized system whose Newton
        ## step equals predict's dx_mna, so fang's machinery works on it
        ## unchanged, eq (18)'s second solve included (same factors).
        j_ra, j_rb, j_IS = pcnr_meta
        v_node = x[j_ra] - x[j_rb]
        i_n, g_n = _junction_terms(v_node, j_IS, pcnr_VT)
        i_l, g_l = _junction_terms(v_lim, j_IS, pcnr_VT)
        F = F.at[j_ra].add(i_l - i_n)
        F = F.at[j_rb].add(-(i_l - i_n))
        J = _scatter_junction_G(J, j_ra, j_rb, g_n, -1.0)
        g_lim = v_lim - v_node
        J = _scatter_junction_G(J, j_ra, j_rb, g_l, +1.0)
        F = F.at[j_ra].add(-g_l * g_lim)
        F = F.at[j_rb].add(g_l * g_lim)
        return F, J, q_C, g_lim

    def body(fs: FangState):
        from pycircuit.circuit._limiting import _pnjlim_branchless
        x, h, v_lim = fs.x, fs.h, fs.v_lim
        f, J, _q_x, g_lim = assemble(x, h, v_lim)
        J_sub = jnp.delete(jnp.delete(J, irefnode, axis=0), irefnode, axis=1)
        dx0_sub = jnp.linalg.solve(J_sub, -jnp.delete(f, irefnode))
        dx0 = jnp.insert(dx0_sub, irefnode, 0.0)
        x1 = x + dx0
        if pcnr_meta is not None:
            ## PCNR's CORRECT phase: the MNA update is taken in FULL (that is
            ## the method's point); only each device's own unknown is limited.
            j_ra, j_rb, j_IS = pcnr_meta
            dx_lim = -(g_lim + (-dx0[j_ra] + dx0[j_rb]))
            v1 = _pnjlim_branchless(v_lim + dx_lim, v_lim, pcnr_VT, j_IS,
                                    jnp.log)
        else:
            v1 = v_lim

        ## Stage 2: has the step size earned any attention?
        lte, h0, _h1 = _fang_lte(x1, state, h, first_order)
        ## P7 deliberately does NOT reach this line: the coupled path keeps
        ## the scalar sigglobal reference its whole measured record -- the
        ## TLine kink analysis included (err = 1/(TRTOL*reltol) is derived
        ## against exactly this ref) -- was built on.  Reconsider with the
        ## coupled estimator, not before.
        ref = jnp.maximum(state.sig_max, jnp.max(jnp.abs(x1)))
        etol = trtol * (reltol * ref + lte_abstol)
        err = jnp.max(jnp.abs(lte) / etol)
        eps_ok = jnp.logical_and(gamma_min <= err, err <= gamma_max)

        conv_x = jnp.all(jnp.abs(dx0)
                         < reltol * jnp.maximum(jnp.abs(x1), jnp.abs(x))
                         + xtol)
        if pcnr_meta is not None:
            ## BOTH residuals: converging on dx alone can return
            ## v_lim != e_a - e_b -- a vector that is not a solution of the
            ## circuit, whose LTE then reads low (the CPU's measured
            ## half-wave lesson, applied to this path from the start).
            conv_x = jnp.logical_and(conv_x, jnp.max(jnp.abs(g_lim)) < (
                reltol * jnp.maximum(jnp.max(jnp.abs(v1)), 1.0) + 1e-12))

        ## A held step whose error is over the band is reported unconverged so
        ## the caller shrinks and retries -- exactly the CPU's hold_h return.
        held_over = jnp.logical_and(hold_h,
                                    jnp.logical_and(~eps_ok, err > gamma_max))

        def plain(_):
            ## hold_h, in-band, or no history: nothing to solve for h.
            return x1, h, v1, conv_x, jnp.asarray(False)

        def solve_h(_):
            ## Sec. 3.4: the new step from the error RATIO, eq (17) inverted
            ## on the node polynomial; eq (18) corrects the solution with the
            ## factors already in hand.
            ratio = target / jnp.maximum(err, 1e-300)
            args1, args2 = [h0], [h0, _h1]
            h_new = jnp.where(
                first_order,
                step_for_error_ratio(h, args1, ratio, 1.0 - eta, 1.0 + eta),
                step_for_error_ratio(h, args2, ratio, 1.0 - eta, 1.0 + eta))
            h_want = jnp.where(
                first_order,
                step_for_error_ratio(h, args1, ratio, 1e-6, 1e6),
                step_for_error_ratio(h, args2, ratio, 1e-6, 1e6))
            h_new = jnp.clip(h_new, h_floor, h_ceil)
            dh = h_new - h
            ## Only a thwarted SHRINK counts (the CPU's saturation lesson).
            saturated = h_want < h_floor * (1.0 - 1e-9)

            ## eq (18): p = d(iq)/dh at fixed solution + du/dt.
            q1 = circuit.q(x1, params_tree=params_tree)
            q_prev = state.q_history[0]
            q_prev2 = state.q_history[1]
            p_e = euler_companion_dh(q1, q_prev, h)
            p_g = bdf2_companion_dh(q1, q_prev, q_prev2, h, h0)
            ## Ordinary sources: analytic du/dt where the circuit provides
            ## it.  A circuit containing a TLine cannot -- TLine.dudt raises
            ## (on the CPU its u IS the delayed history, and the CPU coupled
            ## path is therefore TLine-refused to this day) -- so fall back
            ## to a central difference of u, which is exact off the kinks and
            ## the kinks are breakpoint-held steps where p is never solved.
            ## The TLine part of df/dh is added analytically below either
            ## way: on this backend the delayed EMF is injected separately,
            ## and its derivative is the interpolant's segment slope.
            try:
                du = circuit.dudt(t_prev + h, analysis=analysis,
                                  params_tree=params_tree)
            except NotImplementedError:
                t_c = t_prev + h
                ## eps scales with the STEP, not with absolute time: an
                ## absolute-scaled eps (1e-9 at t << 1) straddled the first
                ## pulse edge from inside the opening ramp, so the derivative
                ## saw the future corner and eq (18) corrupted the solve --
                ## measured as the coupled path dying at t ~ 1e-12 on any
                ## driven TLine circuit.
                eps = 1e-3 * h
                u_hi = circuit.u(t_c + eps, analysis=analysis,
                                 params_tree=params_tree)
                u_lo = circuit.u(t_c - eps, analysis=analysis,
                                 params_tree=params_tree)
                if provided_function is not None:
                    ## The FD pair carries pf' for free; the analytic-dudt
                    ## branch above omits it, exactly as the CPU's
                    ## residual_dh does (recorded F4 gap, both backends).
                    u_hi = u_hi + provided_function(t_c + eps)
                    u_lo = u_lo + provided_function(t_c - eps)
                du = (u_hi - u_lo) / (2.0 * eps)
            p = (jnp.where(first_order, p_e, p_g) + du
                 + tline_source_dudt(x.shape[0], t_prev + h, tline_params,
                                     tline_indices, state.tline_history,
                                     state.tline_head))
            dxh_sub = jnp.linalg.solve(J_sub, -jnp.delete(p, irefnode))
            dxh = jnp.insert(dxh_sub, irefnode, 0.0)
            x_corr = x1 + dxh * dh
            if pcnr_meta is not None:
                ## v_lim tracks the branch voltage: the eq (18) correction
                ## must move it too, or the loop can return with the exact
                ## v_lim != e_a - e_b inconsistency the g_lim test guards
                ## (the CPU learned this as a stage-4 defect).
                j_ra2, j_rb2, _ = pcnr_meta
                v_corr = v1 + (dxh[j_ra2] - dxh[j_rb2]) * dh
            else:
                v_corr = v1

            done = jnp.logical_and(conv_x,
                                   jnp.logical_and(~saturated,
                                                   jnp.abs(dh) < eta * h_new))
            return x_corr, h_new, v_corr, done, jnp.asarray(False)

        take_plain = jnp.logical_or(hold_h,
                                    jnp.logical_or(eps_ok, no_hist))
        x_out, h_out, v_out, done, _ = jax.lax.cond(take_plain, plain,
                                                    solve_h, None)
        done = jnp.logical_or(done, held_over)
        conv = jnp.logical_and(done, ~held_over)
        return FangState(x=x_out, h=h_out, iters=fs.iters + 1,
                         done=done, converged=conv, v_lim=v_out)

    def cond(fs: FangState):
        return jnp.logical_and(~fs.done, fs.iters < maxiter)

    x0 = extrapolate_predictor(state)
    if pcnr_meta is not None:
        _ra0, _rb0, _ = pcnr_meta
        v0 = x0[_ra0] - x0[_rb0]
    else:
        v0 = jnp.zeros(0)
    init = FangState(x=x0, h=h_entry,
                     iters=jnp.asarray(0), done=jnp.asarray(False),
                     converged=jnp.asarray(False), v_lim=v0)
    final = jax.lax.while_loop(cond, body, init)
    ## maxiter exhaustion without done: unconverged.
    return final



## ---------------------------------------------------------------------------
## PCNR on the traced loop -- P19 stage 2.
##
## Fang/Aadithya's replacement for limiting, and the ONLY junction-robustness
## mechanism this backend has: classic limiting cannot run here at all, because
## Diode.limit keeps per-device Python state (_vlim) that a lax.while_loop
## cannot host.  PCNR's limited quantities are just more state vector -- a
## per-junction v_lim array -- which is why the parity review reversed its
## "one-sided by design" verdict.
##
## The junction set is static at trace time (row indices and IS gathered from
## pcnr_junctions); everything per-iteration is pure array arithmetic: the
## MNA assembly minus each junction's stamp at the NODE voltage plus its stamp
## at its OWN unknown, the rank-k Schur update as four scatter-adds per
## junction, and the CORRECT phase as the branchless pnjlim over the v_lim
## vector.  Junction devices with charge storage are refused at setup, exactly
## as the CPU's augmented_system refuses them.
## ---------------------------------------------------------------------------

class PcnrState(NamedTuple):
    x: Any
    v_lim: Any
    iters: Any
    converged: Any


def _junction_arrays(circuit):
    """Static junction metadata: (ra[], rb[], IS[]) -- or None if none."""
    from pycircuit.circuit.pcnr import pcnr_junctions
    junctions = pcnr_junctions(circuit)
    if not junctions:
        return None
    for _inst, element, _ra, _rb in junctions:
        if hasattr(element, 'eval_q_pure') or hasattr(type(element), 'q'):
            ## Mirror of the CPU augmented_system refusal: a charge-storing
            ## junction's iq would be evaluated at the node voltage while the
            ## current moves to v_lim -- the exact inconsistency PCNR removes,
            ## reintroduced through the charge.
            if hasattr(element, 'eval_q_pure'):
                raise NotImplementedError(
                    '%r has charge storage and a PCNR junction; not '
                    'implemented (same refusal as the CPU path).' % _inst)
    ra = jnp.array([j[2] for j in junctions], dtype=jnp.int32)
    rb = jnp.array([j[3] for j in junctions], dtype=jnp.int32)
    IS = jnp.array([float(getattr(j[1].iparv, 'IS', 0.0)) for j in junctions])
    return ra, rb, IS


def _junction_terms(v, j_IS, VT):
    """Per-junction diode current and conductance at voltage vector ``v``."""
    i = j_IS * (jnp.exp(v / VT) - 1.0)
    g = j_IS * jnp.exp(v / VT) / VT
    return i, g


def _scatter_junction_G(J, ra, rb, g, sign):
    """The four-entry-per-junction conductance pattern, scattered."""
    J = J.at[ra, ra].add(sign * g)
    J = J.at[ra, rb].add(-sign * g)
    J = J.at[rb, ra].add(-sign * g)
    J = J.at[rb, rb].add(sign * g)
    return J


def pcnr_inner_loop(state: TransientState, circuit, irefnode, j_ra, j_rb,
                    j_IS, VT, tline_params, tline_indices,
                    eval_method='gear', reltol=1e-4, abstol=1e-12,
                    xtol=1e-12, maxiter=100, params_tree=None, tline_dG=None,
                    analysis='tran', provided_function=None):
    """One time point by PCNR: predict on the Schur-reduced system, correct
    each junction's own unknown with pnjlim.  Returns PcnrState."""
    from pycircuit.circuit._limiting import _pnjlim_branchless

    t_new = state.t + state.dt
    first_order = effective_first_order(state)

    def assemble_base(x):
        I_G = circuit.i(x, params_tree=params_tree)
        G_G = circuit.G(x, params_tree=params_tree)
        q_C = circuit.q(x, params_tree=params_tree)
        C_C = circuit.C(x, params_tree=params_tree)
        I_u = circuit.u(t_new, analysis=analysis, params_tree=params_tree)
        if provided_function is not None:
            I_u = I_u + provided_function(t_new)
        I_u = apply_tline_sources(I_u, t_new, tline_params, tline_indices,
                                  state.tline_history, state.tline_head)
        i_C, G_eq = compute_integration(q_C, C_C, state, method=eval_method,
                                        first_order=first_order)
        F, J = I_G + i_C + I_u, G_G + G_eq
        if tline_dG is not None:
            F = F + tline_dG @ x
            J = J + tline_dG
        return F, J

    def body(ps: PcnrState):
        x, v_lim = ps.x, ps.v_lim
        F, J = assemble_base(x)

        ## Remove each junction's stamp at the NODE voltage (it is inside the
        ## base assembly) and re-enter its current through its OWN unknown.
        v_node = x[j_ra] - x[j_rb]
        i_n, g_n = _junction_terms(v_node, j_IS, VT)
        i_l, g_l = _junction_terms(v_lim, j_IS, VT)
        F = F.at[j_ra].add(i_l - i_n)
        F = F.at[j_rb].add(-(i_l - i_n))
        J = _scatter_junction_G(J, j_ra, j_rb, g_n, -1.0)

        g_lim = v_lim - v_node

        ## Schur: J_eff = J + rank-k(didv at v_lim); f_eff = F - J_ml g_lim.
        J_eff = _scatter_junction_G(J, j_ra, j_rb, g_l, +1.0)
        F_eff = F.at[j_ra].add(-g_l * g_lim)
        F_eff = F_eff.at[j_rb].add(g_l * g_lim)

        J_sub = jnp.delete(jnp.delete(J_eff, irefnode, axis=0),
                           irefnode, axis=1)
        dx_sub = jnp.linalg.solve(J_sub, -jnp.delete(F_eff, irefnode))
        dx = jnp.insert(dx_sub, irefnode, 0.0)
        x_new = x + dx

        dx_lim = -(g_lim + (-dx[j_ra] + dx[j_rb]))
        ## CORRECT: each device limits only the unknown it owns.
        v_new = _pnjlim_branchless(v_lim + dx_lim, v_lim, VT, j_IS, jnp.log)

        ## BOTH residuals, per the CPU's measured lesson: converging on dx
        ## alone can return v_lim != e_a - e_b, a vector that is not a
        ## solution of the circuit, whose LTE then reads low.
        I_scale = jnp.abs(J_eff) @ jnp.abs(x_new) + jnp.abs(F_eff)
        conv_f = jnp.all(jnp.abs(F_eff) < reltol * I_scale + abstol)
        conv_x = jnp.all(jnp.abs(dx)
                         < reltol * jnp.maximum(jnp.abs(x_new), jnp.abs(x))
                         + xtol)
        lim_ok = jnp.max(jnp.abs(g_lim)) < (
            reltol * jnp.maximum(jnp.max(jnp.abs(v_new)), 1.0) + 1e-12)
        conv = jnp.logical_and(jnp.logical_and(conv_x, conv_f), lim_ok)
        return PcnrState(x=x_new, v_lim=v_new, iters=ps.iters + 1,
                         converged=conv)

    def cond(ps: PcnrState):
        return jnp.logical_and(~ps.converged, ps.iters < maxiter)

    x0 = extrapolate_predictor(state)
    init = PcnrState(x=x0, v_lim=x0[j_ra] - x0[j_rb],
                     iters=jnp.asarray(0), converged=jnp.asarray(False))
    return jax.lax.while_loop(cond, body, init)


def pcnr_controller_jacobian(circuit, state, x, v_lim, j_ra, j_rb, j_IS, VT,
                             eval_method, first_order, params_tree=None,
                             tline_dG=None):
    """The Jacobian the step controller must see: the one PCNR SOLVED.

    ``G(x) + Geq`` is not it -- the base G carries the junction stamp at the
    node voltage; the solved matrix replaces it with didv at v_lim (the CPU
    path re-learned this at a 6.6x step-count cost).
    """
    G_G = circuit.G(x, params_tree=params_tree)
    q_C = circuit.q(x, params_tree=params_tree)
    C_C = circuit.C(x, params_tree=params_tree)
    _i_C, G_eq = compute_integration(q_C, C_C, state, method=eval_method,
                                     first_order=first_order)
    J = G_G + G_eq
    if tline_dG is not None:
        J = J + tline_dG
    v_node = x[j_ra] - x[j_rb]
    _i_n, g_n = _junction_terms(v_node, j_IS, VT)
    _i_l, g_l = _junction_terms(v_lim, j_IS, VT)
    J = _scatter_junction_G(J, j_ra, j_rb, g_n, -1.0)
    return _scatter_junction_G(J, j_ra, j_rb, g_l, +1.0)


def outer_time_loop(initial_state: TransientState, circuit, tend, chunk_size, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method='euler', params_tree=None, reltol=1e-4, abstol=1e-12, xtol=1e-12, maxiter=100, trtol=7.0, lte_reltol=1e-4, lte_abstol=1e-12, max_dv=jnp.inf, coupled=False, gamma_min=0.7, gamma_max=3.0, eta=0.15, pcnr_meta=None, pcnr_VT=0.025, tline_dG=None, analysis='tran', s_gamma_min=0.0, s_gamma_max=1.0, s_eta=jnp.inf, fixed_timestep=False, grid_dt=None, relref='sigglobal', n_nodes=None, provided_function=None):

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
        if coupled:
            ## FANG'S COUPLED PATH (P19).  A step that must LAND somewhere --
            ## the next breakpoint or tend (tend is in t_breaks_array) -- has
            ## its size imposed, so it is HELD: fang solves x only, and an
            ## over-band held step is reported for the shrink-and-retry below,
            ## exactly the CPU's hold_h contract (F3 included: tend-truncated
            ## steps cannot overshoot).
            future = jnp.where(t_breaks_array > state.t + t_eps,
                               t_breaks_array, jnp.inf)
            next_break = jnp.min(future)
            overshoots = state.t + state.dt > next_break
            h_entry = jnp.where(overshoots, next_break - state.t, state.dt)
            hold_h = jnp.logical_or(
                overshoots,
                jnp.abs((state.t + h_entry) - next_break) <= t_eps)
            fstate = fang_inner_loop(
                state._replace(dt=h_entry), circuit, irefnode,
                hold_h, tline_params, tline_indices,
                eval_method=eval_method, reltol=reltol,
                xtol=xtol, lte_abstol=lte_abstol, trtol=trtol,
                maxiter=maxiter, params_tree=params_tree,
                gamma_min=gamma_min, gamma_max=gamma_max, eta=eta,
                dt_min=dt_min, dt_max=dt_max, pcnr_meta=pcnr_meta,
                pcnr_VT=pcnr_VT, tline_dG=tline_dG, analysis=analysis,
                provided_function=provided_function)
            ## Re-enter the shared accept machinery with fang's verdict: the
            ## solved h is both the step taken and (clipped) the next guess --
            ## "the solved step carries forward" is the method's whole point.
            state = state._replace(dt=fstate.h)
            nr_state = NewtonState(x=fstate.x, xdiff=jnp.zeros_like(fstate.x),
                                   F_norm=jnp.asarray(0.0),
                                   iters=fstate.iters,
                                   converged=fstate.converged)
        elif pcnr_meta is not None:
            ## PCNR (P19 stage 2): the standard path's Newton, with each
            ## junction's limited quantity as its own unknown.  This backend's
            ## ONLY junction-robustness mechanism -- classic limiting cannot
            ## run in a traced loop.
            j_ra, j_rb, j_IS = pcnr_meta
            pstate = pcnr_inner_loop(
                state, circuit, irefnode, j_ra, j_rb, j_IS, pcnr_VT,
                tline_params, tline_indices,
                eval_method=eval_method, reltol=reltol, abstol=abstol,
                xtol=xtol, maxiter=maxiter, params_tree=params_tree,
                tline_dG=tline_dG, analysis=analysis,
                provided_function=provided_function)
            nr_state = NewtonState(x=pstate.x, xdiff=jnp.zeros_like(pstate.x),
                                   F_norm=jnp.asarray(0.0),
                                   iters=pstate.iters,
                                   converged=pstate.converged)
            pcnr_v_lim = pstate.v_lim
        else:
            nr_state = newton_inner_loop(state, circuit, irefnode, tline_params, tline_indices,
                                         eval_method=eval_method, reltol=reltol, abstol=abstol,
                                         xtol=xtol, maxiter=maxiter,
                                         params_tree=params_tree, max_dv=max_dv,
                                         tline_dG=tline_dG, analysis=analysis,
                                         provided_function=provided_function)
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
        if pcnr_meta is not None and not coupled:
            j_ra, j_rb, j_IS = pcnr_meta
            J = pcnr_controller_jacobian(circuit, state, x_curr, pcnr_v_lim,
                                         j_ra, j_rb, j_IS, pcnr_VT,
                                         eval_method, first_order,
                                         params_tree=params_tree,
                                         tline_dG=tline_dG)
        else:
            J = circuit.G(x_curr, params_tree=params_tree) + Geq
            if tline_dG is not None:
                J = J + tline_dG
        if coupled:
            ## Fang solved (x, h) against its own solution-space LTE; a second
            ## verdict from the charge estimator would re-litigate it.
            error_ratio = jnp.where(nr_state.converged, 0.0, jnp.inf)
            order_p1 = jnp.asarray(2.0)
        else:
            error_ratio, order_p1 = ywr_error_ratio(
                i_curr, x_curr, J, state, irefnode,
                method=eval_method, trtol=trtol, lte_rel=lte_reltol,
                lte_abstol=lte_abstol, first_order=first_order,
                relref=relref, n_nodes=n_nodes)

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
        ## P8: the CPU standard band, same semantics (stepcontroller
        ## set_lte_band / IntegralController.evaluate_step).  At the shipped
        ## defaults (s_gamma_min=0, s_gamma_max=1, s_eta=inf) every line
        ## below reduces to the historical err<=1.0 accept -- default runs
        ## are bit-identical, which the P8 gate pins.
        landed = jnp.any(
            jnp.abs((state.t + state.dt) - t_breaks_array) <= t_eps)
        lte_ok = jnp.logical_or(error_ratio <= s_gamma_max, first)
        ## Eq (15)'s LOWER bound, the CPU's too-accurate growth-reject: a
        ## step far under tolerance is wasted work, redone LARGER at the
        ## same time point -- UNLESS its size was never the controller's to
        ## choose (breakpoint landing), it has nowhere to grow (at dt_max),
        ## or the grow would not actually grow (damper/cap binding would
        ## livelock the retry).
        from pycircuit.circuit.stepcontroller import MAX_GROWTH_RATIO as _MG
        _safety_target = 0.9 ** order_p1
        _aim = jnp.where(
            jnp.logical_and(s_gamma_min <= _safety_target,
                            _safety_target <= s_gamma_max),
            _safety_target, jnp.sqrt(jnp.maximum(s_gamma_min * s_gamma_max,
                                                 1e-30)))
        _gfac = jnp.minimum(_MG, (_aim / jnp.maximum(error_ratio, 1e-12))
                            ** (1.0 / order_p1))
        _gfac = jnp.minimum(_gfac, 1.0 + s_eta)
        grow_dt = jnp.minimum(state.dt * _gfac, dt_max)
        too_accurate = jnp.logical_and(
            jnp.logical_and(nr_state.converged, error_ratio < s_gamma_min),
            jnp.logical_and(
                jnp.logical_not(jnp.logical_or(first, landed)),
                jnp.logical_and(state.dt < dt_max * (1.0 - 1e-12),
                                grow_dt > state.dt * (1.0 + 1e-9))))
        if coupled:
            ## Fang solved (x, h) against its own band; the sentinel
            ## error_ratio=0.0 above must not trip the standard path's
            ## too-accurate reject when the user set a nonzero gamma_min.
            too_accurate = jnp.asarray(False)
        if fixed_timestep:
            ## P9, the CPU's stage 4h: under a fixed grid there is no
            ## controller -- the LTE verdict is skipped entirely, and only a
            ## failed Newton rejects (shrink-retry below; the next accepted
            ## step returns to the grid).
            too_accurate = jnp.asarray(False)
            lte_ok = jnp.asarray(True)
        accept = jnp.logical_or(
            jnp.logical_and(
                jnp.logical_and(nr_state.converged, lte_ok),
                jnp.logical_not(too_accurate)),
            at_floor)
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
            if coupled:
                ## The solved step carries forward; breakpoint truncation for
                ## the NEXT point happens at the next fang entry.
                next_dt = jnp.clip(state.dt, dt_min, dt_max)
            elif fixed_timestep:
                ## P9: THE GRID WINS.  Breakpoints do not move it (the CPU's
                ## stage 4h measured a VPulse run collapsing 30 -> 292 steps
                ## when they did); only `tend` truncates, and a failed step's
                ## shrunken retry returns to the grid here on its next accept
                ## -- the CPU's own warn text ("for this step") states that
                ## intent, though its loop today keeps the shrunken step;
                ## recorded as a finding in the parity doc.
                next_dt = jnp.maximum(
                    jnp.minimum(grid_dt, tend - (state.t + state.dt)), dt_min)
            else:
                next_dt = calculate_next_dt(state.dt, error_ratio, dt_min, dt_max,
                                            t_breaks_array, state.t + state.dt, order_p1,
                                            eta=s_eta)

            x_hist_new = jnp.roll(state.x_history, shift=1, axis=0).at[0].set(x_curr)
            q_hist_new = jnp.roll(state.q_history, shift=1, axis=0).at[0].set(q_curr)
            iq_hist_new = jnp.roll(state.iq_history, shift=1, axis=0).at[0].set(i_curr)
            h_hist_new = jnp.roll(state.h_history, shift=1, axis=0).at[0].set(state.dt)
            if coupled and tline_params.shape[0] > 0:
                ## COUPLED ONLY: a landing zeroes the step ring, restarting the
                ## history as a cold start does.  WHY: history that STRADDLES a
                ## kink poisons the coupled LTE two distinct ways, both traced
                ## on the TLine pulse probe -- (1) the step AFTER the landing
                ## needs its band skipped, because at a from-zero kink dev is
                ## proportional to ref and h cancels (err = 1/(TRTOL*reltol) =
                ## 1428.6 measured, h-independent); (2) the step after THAT
                ## extrapolates degree-2 through a pre-kink point, and a
                ## quadratic through two flat points and one ramp point misses
                ## a line by err = 139.2/1e-11 * h -- in-band only at h ~
                ## 2e-13, unreachable above the within-point floor 0.2*h_entry.
                ## Zeroed history makes no_hist skip the band for one step and
                ## effective_first_order hold degree 1 for two (the F19
                ## run-global facts), after which every history point is
                ## post-kink.  The standard path keeps its measured one-step
                ## force_first_order (F11) -- its integrator-side LTE decays
                ## with h, so it never livelocks; a two-step reset would
                ## change its recorded step counts for no defect.  GATED ON
                ## DELAY LINES for the same measured reason as the CPU site:
                ## unconditional, the reset moved the pulsed-RC coupled-vs-
                ## standard median from 9.8e-4 to 5.4e-3 V at reltol 1e-5 --
                ## capacitive rows never needed the grace, and their kink-
                ## spanning estimate is conservative rather than trapped.
                h_hist_new = jnp.where(landed,
                                       jnp.zeros_like(h_hist_new),
                                       h_hist_new)

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
                ## P7: accepted-only, like sig_max; `local` matches the
                ## estimator's (current point vs the one before it).
                ref_running=jnp.maximum(
                    state.ref_running,
                    jnp.maximum(jnp.abs(x_curr),
                                jnp.abs(state.x_history[0]))),
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
                force_first_order=(jnp.any(jnp.logical_and(
                    t_breaks_array > state.t + t_eps,
                    t_breaks_array <= state.t + state.dt + t_eps))
                    if fixed_timestep else landed))

        def do_reject(_):
            ## LTE above tolerance: shrink the step (bounded below by dt_min) and
            ## retry the same time point without committing anything.
            ## Two different reasons to shrink, and they need different factors.
            ## The predictor's `ratio^(-1/p)` is only meaningful when `x_curr` IS a
            ## solution; after a failed Newton the ratio describes a point the
            ## circuit never occupied, so a fixed cut is used instead.
            lte_shrink = jnp.clip(error_ratio ** (-1.0 / order_p1), 0.1, 0.9)
            shrink = jnp.where(nr_state.converged, lte_shrink, 0.25)
            if coupled:
                ## The CPU coupled path's outer retry: h *= 0.25 whatever the
                ## failure flavour (Newton or held-step LTE).
                shrink = jnp.asarray(0.25)
            retry_dt = jnp.maximum(state.dt * shrink, dt_min)
            ## P8: a too-accurate rejection retries LARGER (the CPU's eq (15)
            ## lower bound); grow_dt's guards above ensure it genuinely grows.
            retry_dt = jnp.where(too_accurate, grow_dt, retry_dt)
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
                 'force_accepts')

    def __init__(self):
        self.accepted_steps = 0
        self.rejected_steps = 0
        self.signal_max = 0.0
        self.nonconverged_steps = 0
        self.forced_nonconverged_steps = 0
        ## F19(c): the JAX analogue of the CPU's force_accepts -- converged
        ## Newton, failing LTE, accepted at the dt floor.
        self.force_accepts = 0

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
        ## P19: Fang's coupled (x, h) stepping, 'approx' branch -- the same
        ## contract as Transient's coupled_lte solve argument, as a Parameter
        ## here because the traced loop is built per run.
        ## P19 stage 2: PCNR instead of limiting -- and on this backend,
        ## instead of NOTHING, since Python-stateful limiting cannot run in a
        ## traced loop.  Off by default, matching the CPU.
        Parameter(name='pcnr',
                  desc='Use Predictor/Corrector Newton-Raphson for junction '
                       'devices (Aadithya et al.); the only junction-'
                       'robustness mechanism available on this backend.',
                  unit='', default=False),
        Parameter(name='coupled_lte',
                  desc="Solve solution and step size together (Fang DAC 2013, "
                       "sec 3.4 'approx' correction). The 'bordered' branch "
                       "and PCNR-inside-Fang remain CPU-only.",
                  unit='', default=False),
        ## The LTE acceptance band, with the CPU's 'auto' sentinel semantics
        ## (F5): 'auto' resolves to Fang's (0.7, 3.0, 0.15) on the coupled
        ## path; explicit values pass through verbatim.  The standard JAX
        ## path does not read the band yet (parity item P8).
        Parameter(name='lte_gamma_min',
                  desc="Lower edge of the LTE acceptance band; 'auto' means "
                       "0.7 (Fang) on the coupled path and 0.0 (no lower "
                       "test) on the standard path -- the CPU's F5 split.",
                  unit='', default='auto'),
        Parameter(name='lte_gamma_max',
                  desc="Upper edge of the LTE acceptance band (coupled path); "
                       "'auto' means 3.0 (Fang).",
                  unit='', default='auto'),
        Parameter(name='lte_eta',
                  desc="Relative per-iteration step-change limit, eq (16); "
                       "'auto' means 0.15 (Fang).",
                  unit='', default='auto'),
        ## F16: opt-in, OFF by default -- the old hardcoded 0.5 V made any
        ## swing beyond ~maxiter*0.5 V non-convergent by construction.
        Parameter(name='max_dv',
                  desc='Largest Newton update per iteration, as a voltage; a '
                       'convergence aid for stiff nonlinear circuits. None '
                       '(default) disables the clamp.',
                  unit='V', default=None),
        ## STAGE 9(g).  Ported from `Transient`, same name, default and validation.
        ## DECISION D2 -- the same knob as the CPU's, same name and same
        ## default, because 9(d) and 9(g) both exist because these two backends
        ## had drifted apart on exactly this kind of detail.  OWNER DECISION
        ## (2026-08-21): decoupled from `timestep` and renamed -- `timestep`
        ## doubling as the cap made gentle circuits step-cap-limited, where no
        ## tolerance knob could move the run (measured: identical 209-step
        ## rc-vsin runs at reltol 1e-4 and 1e-6).
        Parameter(name='timestep_max',
                  desc='Largest accepted timestep; None means tend/50, the '
                       'SPICE TMAX default. Decoupled from timestep, which '
                       'only sets the opening-step scale',
                  unit='s', default=None),
        ## P1: `uic`/`minstep` were **kwargs reads -- `solve(uicc=True)` ran
        ## silently with defaults, the dead-knob defect class at the call
        ## boundary.  Declared as Parameters (CPU names, CPU defaults); the
        ## solve()/solve_batched() arguments remain as explicit per-call
        ## overrides, None meaning "use the Parameter".
        ## P5's other half: the inherited `analysis` default is not 'tran',
        ## and a threaded self.par.analysis of anything else makes every
        ## source's u() return zeros -- caught by the rc-charging gate as an
        ## all-zero waveform the moment the hardcode was removed.  Same
        ## re-declaration the CPU Transient makes.
        Parameter(name='analysis', desc='Analysis name', default='tran'),
        ## P6: the integrator choice, reachable at last -- the traced loop
        ## implemented both methods all along, but eval_method was hardcoded
        ## at both call sites.  String-valued on this backend (the traced
        ## kernels select by name; there is no Integrator instance to hold
        ## state).  Trapezoidal stays CPU-only until someone ports a
        ## VARIABLE-STEP trap estimator: the uniform-grid trap branch was
        ## deleted for cause (its LTE formula assumed equal steps).
        Parameter(name='integrator',
                  desc="Integration method: 'gear' (default, order 2, the "
                       "CPU's default too) or 'euler' (order 1). "
                       'Trapezoidal is CPU-only -- see the source note',
                  unit='', default='gear'),
        Parameter(name='uic',
                  desc='Use initial conditions (skip DC OP computation)',
                  unit='', default=False),
        ## P12: SPICE's .ic for the uic case, shared with the CPU -- the
        ## initial-state machinery (ic dict -> element ICs -> spanning-tree
        ## capacitor solve) is pure pre-loop Python on names and indices, so
        ## the CPU's methods are bound below unchanged.
        ## P7: the CPU's relref, same values, same default; 'sigglobal' on
        ## this backend now splits unit groups (node voltages vs branch
        ## currents) exactly as the CPU does -- the pre-P7 scalar reference
        ## mixed volts and amps.  The COUPLED path keeps its scalar
        ## reference (its measured record depends on it; see fang's note).
        Parameter(name='relref',
                  desc="Reference for the relative LTE tolerance: 'sigglobal' "
                       "(against the largest signal in the unknown's unit "
                       "group -- the default, as in Spectre), 'pointlocal' "
                       "(each unknown against itself), or 'alllocal' "
                       "(against its own past maximum)",
                  unit='', default='sigglobal'),
        Parameter(name='ic',
                  desc="Initial node voltages for uic=True, as {node: volts}. "
                       "Node may be a name or a Node instance.",
                  unit='V', default=None),
        Parameter(name='minstep',
                  desc='Minimum timestep to prevent infinite loops', unit='s',
                  default=1e-18),
        ## P14: same spacing rule as the CPU's breakpoint handling.
        Parameter(name='minbreak',
                  desc='Minimum time difference for breakpoint events',
                  unit='s', default=1e-14),
        ## P10: same uniform-output-grid contract as the CPU path;
        ## resample_uniform is a free function on arrays and the result is
        ## numpy by the time it exists.
        Parameter(name='outputstep',
                  desc='Spacing of a UNIFORM output grid; None returns the '
                       "solver's own adaptive points. Results are "
                       'interpolated quadratically, to match the integrator '
                       '-- see resample_uniform',
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

    ## P12: the CPU's initial-state machinery, SHARED rather than ported --
    ## `_initial_state` and its helpers are pre-loop Python over names and
    ## indices (they build a numpy vector; the chunk loop converts), so the
    ## bound functions run here unchanged, spanning-tree capacitor solve
    ## included.  The old `node.ic` attribute walk this replaces was dead
    ## code posing as a feature: nothing in the package or tests ever SET
    ## node.ic, and a misspelled node name in it could not even be detected.
    from pycircuit.circuit.transient import Transient as _CPU
    _initial_state = _CPU._initial_state
    _apply_element_ics = _CPU._apply_element_ics
    _apply_voltage_ics = _CPU._apply_voltage_ics
    _descendant_has_ic = _CPU._descendant_has_ic
    del _CPU

    def _relref(self):
        """Validated relref mode -- P7."""
        mode = self.par.relref
        if mode not in ('sigglobal', 'pointlocal', 'alllocal'):
            raise ValueError(
                "relref must be 'sigglobal', 'pointlocal' or 'alllocal', "
                "not %r" % (mode,))
        return mode

    def _eval_method(self):
        """The traced loop's method name, validated -- P6."""
        method = self.par.integrator
        if method not in ('gear', 'euler'):
            raise ValueError(
                "integrator must be 'gear' or 'euler' on this backend, not %r. "
                "Trapezoidal is CPU-only: the uniform-grid trap branch was "
                "deleted for cause and a variable-step trap estimator has not "
                "been written." % (method,))
        return method

    def _timestep_max(self, tend):
        """The clamp on how large an accepted step may grow -- decision D2.

        `Transient` resolves identically: None means tend/50, SPICE's TMAX
        default.  Decoupled from `timestep` by owner decision -- the old
        None-means-timestep coupling made gentle circuits step-cap-limited,
        with no tolerance knob able to move the run.
        """
        timestep_max = self.par.timestep_max
        if timestep_max is None or timestep_max <= 0:
            return float(tend) / 50.0
        return float(timestep_max)

    def _element_cap(self, dt_max):
        """Clamp ``dt_max`` to the tightest element step cap (stage 8(d)
        parity).  A TLine cannot be resolved by steps as long as its delay --
        the CPU measured the observed delay at 2.00x TD when dt = TD, with no
        warning -- and this backend never applied the cap at all: standard-path
        runs were correct only because the adaptive controller happened to
        keep steps small, and the coupled path solves h freely and would walk
        straight past TD/2 on any quiet stretch."""
        cap = self.cir.max_timestep() \
            if hasattr(self.cir, 'max_timestep') else None
        if cap is not None and cap < dt_max:
            return float(cap)
        return dt_max

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
        firststep = self.par.firststep
        if firststep is None:
            return timestep * 1e-3
        if firststep <= 0:
            raise ValueError(
                "firststep must be positive, not %r; pass None to use the default "
                "ramp of timestep*1e-3" % (firststep,))
        ## Opening above the suggested scale is more likely a mistake than
        ## an intent; the controller grows geometrically from here anyway.
        return min(float(firststep), float(timestep))

    def _pcnr_setup(self):
        """(ra, rb, IS) arrays + VT, or (None, VT): static trace-time junction
        metadata.  Charge-storing junction devices are refused inside."""
        from pycircuit.circuit.circuit import defaultepar
        epar = getattr(self, 'epar', defaultepar) or defaultepar
        T = getattr(epar, 'T', 300.0)
        VT = float(self.toolkit.kboltzmann) * float(T) \
            / float(self.toolkit.qelectron)
        if not self.par.pcnr:
            return None, VT
        meta = _junction_arrays(self.cir)
        if meta is None:
            ## No participating device: fall through to the plain Newton,
            ## as the CPU's solve_timestep does.
            return None, VT
        return meta, VT

    def _standard_band(self):
        """The STANDARD path's (gamma_min, gamma_max, eta) -- P8.

        Mirrors the CPU's set_lte_band resolution exactly (F5 semantics):
        'auto' means (0.0, 1.0, None) here -- the historical err<=1.0 accept
        -- while the coupled path's 'auto' means Fang's (0.7, 3.0, 0.15);
        the same three Parameters serve both, resolved per path.  eta=None
        becomes inf so the traced damper is branchless and inert."""
        gm, gx, eta = (self.par.lte_gamma_min, self.par.lte_gamma_max,
                       self.par.lte_eta)
        gm = 0.0 if gm == 'auto' else float(gm)
        gx = 1.0 if gx == 'auto' else float(gx)
        eta = None if eta == 'auto' else eta
        if not (0.0 <= gm < gx):
            raise ValueError(
                'LTE band requires 0 <= gamma_min < gamma_max, got %r, %r'
                % (gm, gx))
        if eta is not None and not eta > 0.0:
            raise ValueError(
                'LTE band damper eta must be positive, got %r' % (eta,))
        return gm, gx, float('inf') if eta is None else float(eta)

    def _coupled_band(self):
        """The (gamma_min, gamma_max, eta) the coupled path runs with --
        Transient._coupled_band's exact mapping (F5): 'auto' resolves to
        Fang's values, explicit numbers (0.0 / 1.0 / None included) pass
        through verbatim."""
        gm, gx, eta = (self.par.lte_gamma_min, self.par.lte_gamma_max,
                       self.par.lte_eta)
        return (0.7 if gm == 'auto' else gm,
                3.0 if gx == 'auto' else gx,
                0.15 if eta == 'auto' else eta)

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

    ## P3: `dt_min=1e-15` was a second, differently-named, three-decades-
    ## looser floor for the same physical quantity `solve` calls `minstep`
    ## -- both entry points now read the one `minstep` Parameter (per-call
    ## override under the same name).  The first argument accepted a node
    ## object and called get_node_index on it, so its old name `irefnode`
    ## lied; renamed `refnode`, defaulting to the same object as everywhere
    ## else (P4).
    def solve_batched(self, refnode=gnd, override_params_tree=None, tend=1e-3,
                      timestep=1e-12, CHUNK_SIZE=5000, dt_max=None, uic=None,
                      minstep=None):
        if override_params_tree is None:
            override_params_tree = {}
        dt_min = float(self.par.minstep) if minstep is None else minstep
        uic = bool(self.par.uic) if uic is None else uic
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
            irefnode = self.cir.get_node_index(refnode)
        else:
            irefnode = refnode
            
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
            ## `tend/50` matches `solve` and matches the CPU's resolution.
            dt_max = self._timestep_max(tend)
        dt_max = self._element_cap(dt_max)
            
        t_breaks_array = jnp.array(collect_breakpoints(
            self.cir, tend, float(self.par.minbreak)))
        
        # Setup batched initial state -- the SAME shared machinery as solve
        # (P12): every lane starts from the ic-resolved vector, since the
        # overrides sweep element parameters, not starting states.
        x0_batch = jnp.tile(jnp.asarray(self._initial_state(refnode),
                                        dtype=jnp.float64), (batch_size, 1))
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
        _tline_dG = (tline_stamp_correction(
            n, [(self.cir.elementnodemap[name], float(t.iparv.Z0))
                for name, t in tlines]) if n_tlines > 0 else None)
        
        _gm, _gx, _eta = self._coupled_band()
        _sgm, _sgx, _seta = self._standard_band()
        _pcnr_meta, _pcnr_VT = self._pcnr_setup()

        def run_chunk(s, p_tree):
            return outer_time_loop(s, self.cir, tend, CHUNK_SIZE, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method=self._eval_method(), params_tree=p_tree, reltol=self.par.reltol, abstol=self._newton_abstol(n),
                                   xtol=self._newton_xtol(n),
                                   maxiter=self.par.maxiter, trtol=self.par.TRTOL, analysis=self.par.analysis, s_gamma_min=_sgm, s_gamma_max=_sgx, s_eta=_seta, relref=self._relref(), n_nodes=len(self.cir.nodes),
                                   lte_reltol=self.par.reltol,
                                   lte_abstol=self._lte_abstol(n),
                                   max_dv=(jnp.inf if self.par.max_dv is None
                                           else float(self.par.max_dv)),
                                   coupled=bool(self.par.coupled_lte),
                                   gamma_min=_gm, gamma_max=_gx, eta=_eta,
                                   pcnr_meta=_pcnr_meta, pcnr_VT=_pcnr_VT,
                                   tline_dG=_tline_dG)
            
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
        b_ref_running = jnp.zeros((batch_size, n))
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
                sig_max=b_sig_max, ref_running=b_ref_running,
                n_rejected=b_n_rejected,
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
            b_ref_running = final_state.ref_running
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

    ## P1: no **kwargs sink -- a misspelled argument raises TypeError instead
    ## of running silently with defaults.  P4: `refnode=gnd`, the CPU's
    ## default object; `0` resolved gnd's index only by accident of ordering.
    def solve(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6,
              CHUNK_SIZE=5000, uic=None, minstep=None, fixed_timestep=False,
              provided_function=None):
        n = self.cir.n
        irefnode = self.cir.get_node_index(refnode)
        uic = bool(self.par.uic) if uic is None else uic
        ## P11: the CPU's F4 caveat, verbatim -- the seeding DC knows
        ## nothing about the extra source, and a discontinuous pf gets
        ## neither breakpoint truncation nor an order drop.  On THIS backend
        ## pf must additionally be jax-traceable (pure, jnp-composable); it
        ## is baked into the compiled chunk at trace time.
        if provided_function is not None and x0 is None and not uic:
            warnings.warn(
                'jaxtransient: provided_function adds a source the DC '
                'operating point does not see, so the run opens from an '
                'inconsistent state and integrates a spurious startup '
                'transient. Pass uic=True or an explicit x0.',
                RuntimeWarning, stacklevel=2)
        if fixed_timestep and self.par.coupled_lte:
            ## The CPU's coupled path locks the grid via grid_locked; wiring
            ## the same through fang's hold machinery is real work not yet
            ## done on this backend.  Refused rather than approximated.
            raise NotImplementedError(
                'fixed_timestep with coupled_lte is not implemented on this '
                'backend yet; the CPU path supports it (grid_locked). Use '
                'coupled_lte=False for a fixed grid here.')

        ## Same contract as the CPU (P12): ic without uic is a different
        ## feature (constraining the operating point) and is refused, not
        ## silently ignored.
        if (self.par.ic or self._descendant_has_ic(self.cir)) and not uic:
            raise ValueError(
                "ic was given without uic=True. This implements SPICE's "
                "initial conditions for the uic case only -- starting values "
                "for the transient. Pass uic=True, or drop ic.")
        if x0 is None:
            if uic:
                x0 = self._initial_state(refnode)
            else:
                # Run a fast DC operating point if none provided
                from pycircuit.circuit.dcanalysis import DC
                x0 = DC(self.cir).solve().x
            
        dt_min = float(self.par.minstep) if minstep is None else minstep
        dt_max = self._element_cap(self._timestep_max(tend))
        if fixed_timestep and timestep > dt_max * (1.0 + 1e-12):
            ## The CPU's stage-8(d) warning, same trade: the caller owns the
            ## grid, so obey it and say what it costs.
            warnings.warn(
                'jaxtransient: fixed_timestep=%g exceeds the %g s cap the '
                'circuit needs (element TD/2 or timestep_max). The result on '
                'this grid is degraded and nothing else will report it.'
                % (timestep, float(dt_max)), RuntimeWarning)
    
        t_breaks_array = jnp.array(collect_breakpoints(
            self.cir, tend, float(self.par.minbreak)))
    
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
        _tline_dG = (tline_stamp_correction(
            n, [(self.cir.elementnodemap[name], float(t.iparv.Z0))
                for name, t in tlines]) if n_tlines > 0 else None)
    
        # Setup initial state
        x_hist = jnp.array([x0, x0, x0])
        q0 = self.cir.q(x0)
        q_hist = jnp.array([q0, q0, q0])
        iq_hist = jnp.zeros((3, n))
        h_hist = jnp.zeros(3)
    
        # We need a JIT-able wrapper that burns the chunk_size and tend into static parameters if needed,
        # but tend is a runtime parameter.
        _gm, _gx, _eta = self._coupled_band()
        _sgm, _sgx, _seta = self._standard_band()
        _pcnr_meta, _pcnr_VT = self._pcnr_setup()

        @jax.jit
        def run_chunk(s):
            return outer_time_loop(s, self.cir, tend, CHUNK_SIZE, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method=self._eval_method(), reltol=self.par.reltol, abstol=self._newton_abstol(n),
                                   fixed_timestep=fixed_timestep, grid_dt=timestep,
                                   provided_function=provided_function,
                                   xtol=self._newton_xtol(n),
                                   maxiter=self.par.maxiter, trtol=self.par.TRTOL, analysis=self.par.analysis, s_gamma_min=_sgm, s_gamma_max=_sgx, s_eta=_seta, relref=self._relref(), n_nodes=len(self.cir.nodes),
                                   lte_reltol=self.par.reltol,
                                   lte_abstol=self._lte_abstol(n),
                                   max_dv=(jnp.inf if self.par.max_dv is None
                                           else float(self.par.max_dv)),
                                   coupled=bool(self.par.coupled_lte),
                                   gamma_min=_gm, gamma_max=_gx, eta=_eta,
                                   pcnr_meta=_pcnr_meta, pcnr_VT=_pcnr_VT,
                                   tline_dG=_tline_dG)
        
        results_list = [np.array([x0])]
        times_list = [np.array([0.0])]
    
        current_t = 0.0
        current_dt = timestep if fixed_timestep else self._opening_step(timestep)
        ## Accepted steps so far, tracked only to detect the TLine ring wrapping.
        ## The check cannot live inside the compiled loop -- it cannot raise there --
        ## so it is done here, per chunk, which is soon enough: the wrap corrupts the
        ## interpolation from the step after it happens, and a chunk boundary is at
        ## most CHUNK_SIZE steps later.
        total_steps = 0
        sig_max = jnp.array(0.0)
        ref_running = jnp.zeros(n)
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
                sig_max=sig_max, ref_running=ref_running,
                n_rejected=n_rejected,
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
            ref_running = final_state.ref_running
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
        self.statistics.force_accepts = int(n_forced_lte)

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
            if fixed_timestep:
                ## The CPU warns per event; a traced loop cannot, so the count
                ## arrives once at the end -- same information, same class.
                warnings.warn(
                    'jaxtransient: Newton did not converge on %d step attempt(s) '
                    'at the requested fixed timestep %g s; each fell back to a '
                    'smaller step. The output grid is no longer uniform.'
                    % (int(n_nonconverged), timestep), RuntimeWarning)
            else:
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
        ## P10: the CPU's uniform-output-grid contract, same Parameter, same
        ## free function -- the arrays are numpy by now, nothing about the
        ## resample is backend-specific.
        if self.par.outputstep is not None:
            from pycircuit.circuit.transient import resample_uniform
            _grid, _Xg = resample_uniform(res.sweep_values, res.x,
                                          step=self.par.outputstep)
            res = CircuitResult(self.cir, x=_Xg, xdot=None,
                                sweep_values=_grid,
                                sweep_label='time', sweep_unit='s')
        ## Reachable from the result, matching the CPU paths (F13): a caller
        ## who kept only the waveform can still ask what produced it.
        res.statistics = self.statistics
        return res
