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


def _gmin_junction_rows(circuit):
    """(ra[], rb[]) of every declared junction -- the P18 gmin ladder's
    scatter targets.  No charge refusal here: that constraint is PCNR's
    (its limited unknown must not split a charge evaluation); a parallel
    conductance across a charge-storing junction is unconditionally fine,
    and at DC the charge is inert anyway."""
    from pycircuit.circuit.pcnr import pcnr_junctions
    junctions = pcnr_junctions(circuit)
    if not junctions:
        return None
    ra = jnp.asarray([r for _i, _e, r, _ in junctions], dtype=jnp.int32)
    rb = jnp.asarray([r for _i, _e, _, r in junctions], dtype=jnp.int32)
    return ra, rb


## The continuation ladders are ADAPTIVE (owner request): exponent-space
## marching with halving-on-failure, mirroring the CPU driver -- the fixed
## decade schedule was measured stranding rungs whose solutions sat too
## far apart.  They still end at EXACTLY zero: only the pure-rung solution
## may be accepted (residue in a committed point is the P22
## inconsistency-floor trap by its measured signature).


def _adaptive_ladder_traced(rung_solve, x0, e_start, e_end, e_max=0.0,
                            min_step=0.25, max_rungs=60):
    """Traced twin of nrsolver._adaptive_conductance_ladder: one compiled
    rung solver (g traced) inside a lax.while_loop carrying (seed,
    exponent, step, last-converged exponent, phase flags).  Returns
    (x, converged) where converged means the g = 0 rung landed."""
    def cond(st):
        (_x, _e, step, _le, _hl, _pure, done, fail, rungs) = st
        return jnp.logical_and(
            jnp.logical_not(jnp.logical_or(done, fail)),
            rungs < max_rungs)

    def body(st):
        (x_seed, e, step, last_e, has_last, pure, done, fail, rungs) = st
        g = jnp.where(pure, 0.0, 10.0 ** e)
        x_try, conv = rung_solve(x_seed, g)

        ## Success bookkeeping.
        x_next = jnp.where(conv, x_try, x_seed)
        done_n = jnp.logical_and(conv, pure)
        pure_n = jnp.where(conv,
                           jnp.logical_or(pure, e <= e_end + 1e-12),
                           jnp.asarray(False))
        last_e_n = jnp.where(jnp.logical_and(conv, jnp.logical_not(pure)),
                             e, last_e)
        has_last_n = jnp.logical_or(has_last,
                                    jnp.logical_and(conv,
                                                    jnp.logical_not(pure)))
        e_march = jnp.maximum(e - step, e_end)

        ## Failure bookkeeping: escalate before any rung has landed,
        ## refine (halve the step from the last landing) afterwards.
        step_f = step / 2.0
        e_escalate = jnp.minimum(e + step, e_max)
        can_escalate = jnp.logical_and(jnp.logical_not(has_last),
                                       e < e_max - 1e-12)
        e_refine = jnp.maximum(
            jnp.where(has_last, last_e, e_start) - step_f, e_end)
        e_fail = jnp.where(can_escalate, e_escalate, e_refine)
        step_n_fail = jnp.where(can_escalate, step, step_f)
        fail_n = jnp.logical_and(jnp.logical_not(can_escalate),
                                 step_f < min_step)

        e_n = jnp.where(conv, e_march, e_fail)
        step_n = jnp.where(conv, step, step_n_fail)
        return (x_next, e_n, step_n, last_e_n, has_last_n, pure_n,
                done_n, jnp.where(conv, jnp.asarray(False), fail_n),
                rungs + 1)

    init = (x0, jnp.asarray(e_start), jnp.asarray(2.0),
            jnp.asarray(e_start), jnp.asarray(False), jnp.asarray(False),
            jnp.asarray(False), jnp.asarray(False), jnp.asarray(0))
    final = jax.lax.while_loop(cond, body, init)
    return final[0], final[6]


def dc_operating_point(circuit, irefnode, n, params_tree=None,
                       reltol=1e-4, abstol=1e-12, xtol=1e-12, maxiter=100,
                       analysis='dc', x0=None, gmin_rows=None, gmin=0.0,
                       gshunt_nodes=0, gshunt=0.0, ptc_g=0.0,
                       ptc_anchor=None):
    """The batched DC operating point -- P21, the roadmap's last refusal.

    ``F = i(x) + u(0, 'dc')``, ``J = G(x)``: exactly the CPU DC class's
    assembly, as a traced Newton on the reduced system so ``jax.vmap`` can
    run every lane of a parameter sweep concurrently.  Convergence is
    scored on the REDUCED rows (the P11 lesson: the reference row is not
    an equation of the solved system).  Returns ``(x, converged)``; the
    caller decides what a non-converged lane means.  No limiting and no
    PCNR here yet -- a junction circuit that needs help at DC should start
    from ``uic=True`` (or wait for gmin continuation, P18); the raise at
    the call site says so.
    """
    if x0 is None:
        x0 = jnp.zeros(n)

    def cond_fun(st):
        x, converged, iters = st
        return jnp.logical_and(jnp.logical_not(converged), iters < maxiter)

    def body_fun(st):
        x, _converged, iters = st
        F = circuit.i(x, params_tree=params_tree) \
            + circuit.u(0.0, analysis=analysis, params_tree=params_tree)
        J = circuit.G(x, params_tree=params_tree)
        ## P18: junction-parallel gmin (the physical homotopy -- a leaky
        ## diode; the linear subnetwork is untouched) and, separately, the
        ## node-to-ground gshunt (SPICE3's diagonal stepping, the
        ## singular-G rescue).  Both are static per rung; gmin/gshunt = 0
        ## compiles them away.
        ## g operands are TRACED (adaptive rungs reuse ONE compiled rung
        ## solver); the structural skip stays static, and g = 0.0 adds
        ## exact zeros.
        if gmin_rows is not None:
            ra, rb = gmin_rows
            vj = x[ra] - x[rb]
            F = F.at[ra].add(gmin * vj)
            F = F.at[rb].add(-gmin * vj)
            J = _scatter_junction_G(J, ra, rb,
                                    jnp.full(ra.shape, 1.0) * gmin, 1.0)
        if gshunt_nodes > 0:
            idx = jnp.arange(gshunt_nodes)
            F = F.at[idx].add(gshunt * x[idx])
            J = J.at[idx, idx].add(gshunt)
        ## P25: the pseudo-transient (Psi-tc) term -- a backward-Euler step
        ## of dx/dtau = -F(x) anchored at the previous pseudo iterate,
        ## g = 1/delta on the full diagonal (the refnode row's share falls
        ## out of the reduced solve).  g = 0 makes the term exact zeros,
        ## so the pure rung IS the pure system.
        if ptc_anchor is not None:
            F = F + ptc_g * (x - ptc_anchor)
            J = J + ptc_g * jnp.eye(n)
        J_sub = jnp.delete(jnp.delete(J, irefnode, axis=0),
                           irefnode, axis=1)
        F_sub = jnp.delete(F, irefnode)
        dx = jnp.insert(jnp.linalg.solve(J_sub, -F_sub), irefnode, 0.0)
        x_next = x + dx
        I_scale = jnp.abs(J) @ jnp.abs(x_next) + jnp.abs(F)
        conv_f = jnp.all(jnp.delete(
            jnp.abs(F) < reltol * I_scale
            + jnp.broadcast_to(abstol, F.shape), irefnode))
        conv_x = jnp.all(jnp.abs(dx)
                         < reltol * jnp.maximum(jnp.abs(x_next), jnp.abs(x))
                         + xtol)
        return (x_next, jnp.logical_and(conv_f, conv_x), iters + 1)

    x, converged, _iters = jax.lax.while_loop(
        cond_fun, body_fun, (x0, jnp.asarray(False), 0))
    return x, converged


def dc_with_continuation(circuit, irefnode, n, n_nodes, params_tree=None,
                         reltol=1e-4, abstol=1e-12, xtol=1e-12, maxiter=100,
                         gmin_rows=None):
    """P18 + P25: plain Newton, then the junction-gmin ladder, then the
    gshunt ladder, then pseudo-transient (Psi-tc) continuation -- each
    rung seeded from the last, each ladder ending at EXACTLY zero, and
    only a zero-rung converged solution reported as converged.  The
    Psi-tc ladder anchors each rung at ITS OWN seed (the moving anchor:
    a backward-Euler pseudo-time step), which is what rescues seeds the
    zero-anchored deformations abandon.  All rungs are compiled
    unconditionally (small DC graphs); `lax.cond` skips executing the
    ladders when an earlier stage lands."""
    def plain(x_seed):
        return dc_operating_point(circuit, irefnode, n,
                                  params_tree=params_tree, reltol=reltol,
                                  abstol=abstol, xtol=xtol, maxiter=maxiter,
                                  x0=x_seed)

    x, conv = plain(None)

    def run_ladder(x_conv):
        x_in, conv_in = x_conv

        def gmin_ladder(x_seed):
            def rung(xs, g):
                return dc_operating_point(
                    circuit, irefnode, n, params_tree=params_tree,
                    reltol=reltol, abstol=abstol, xtol=xtol,
                    maxiter=maxiter, x0=xs, gmin_rows=gmin_rows, gmin=g)
            return _adaptive_ladder_traced(rung, x_seed,
                                           e_start=-2.0, e_end=-12.0)

        def gshunt_ladder(x_seed):
            def rung(xs, g):
                return dc_operating_point(
                    circuit, irefnode, n, params_tree=params_tree,
                    reltol=reltol, abstol=abstol, xtol=xtol,
                    maxiter=maxiter, x0=xs, gshunt_nodes=n_nodes, gshunt=g)
            return _adaptive_ladder_traced(rung, x_seed,
                                           e_start=-3.0, e_end=-12.0)

        def ptc_ladder(x_seed):
            ## P25: each rung is anchored at ITS seed -- the traced driver
            ## already reseeds every rung from the last converged iterate,
            ## so passing the seed as the anchor IS the pseudo-transient
            ## march.  delta = 1/g grows 1 s -> 1e12 s; heavier damping
            ## (e_max = +6) is the escalation if even the first step fails.
            def rung(xs, g):
                return dc_operating_point(
                    circuit, irefnode, n, params_tree=params_tree,
                    reltol=reltol, abstol=abstol, xtol=xtol,
                    maxiter=maxiter, x0=xs, ptc_g=g, ptc_anchor=xs)
            return _adaptive_ladder_traced(rung, x_seed,
                                           e_start=0.0, e_end=-12.0,
                                           e_max=6.0)

        if gmin_rows is not None:
            xg, cg = gmin_ladder(jnp.zeros(n))
        else:
            xg, cg = x_in, conv_in
        ## gshunt only when the physical ladder failed -- the crude
        ## rescue -- and Psi-tc only when that failed too (P25, the
        ## chain's last rung, mirroring the CPU order).
        xs_, cs_ = jax.lax.cond(cg, lambda _: (xg, cg),
                                lambda _: gshunt_ladder(jnp.zeros(n)), None)
        return jax.lax.cond(cs_, lambda _: (xs_, cs_),
                            lambda _: ptc_ladder(jnp.zeros(n)), None)

    return jax.lax.cond(conv, lambda xc: xc, run_ladder, (x, conv))


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
                    provided_function=None, state_mask=None):
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
        j_ra, j_rb, j_IS, j_VT, j_fns = pcnr_meta
        v_node = x[j_ra] - x[j_rb]
        if j_fns is not None:
            i_n, g_n = _junction_terms_generic(v_node, j_fns)
            i_l, g_l = _junction_terms_generic(v_lim, j_fns)
        else:
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
            j_ra, j_rb, j_IS, j_VT, _fns = pcnr_meta
            dx_lim = -(g_lim + (-dx0[j_ra] + dx0[j_rb]))
            v1 = _pnjlim_branchless(v_lim + dx_lim, v_lim, j_VT, j_IS,
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
        ## P22: eq (6) over the STATE rows only, same mechanism as the CPU's
        ## _lte_tolerance -- algebraic rows (zero row AND column of C) carry
        ## grid-convention deviations with an h-independent floor, measured
        ## as the Gear-2 rectifier livelock on the source-current row.  1e30
        ## rather than inf keeps every downstream gradient finite.
        if state_mask is not None:
            etol = jnp.where(state_mask, etol, 1e30)
        err = jnp.max(jnp.abs(lte) / etol)
        ## NaN -> inf, found by the P22 mask port: a diverging cold-start
        ## Newton makes lte NaN on the state rows, and with the algebraic
        ## rows masked to ~0 the max is NaN -- every band predicate then
        ## reads False, fang's h-update propagates the NaN, and the outer
        ## loop can neither accept nor reach the dt floor: an in-trace
        ## livelock.  inf turns it into the unmasked path's deterministic
        ## hard shrink toward the floor, where the forced-accept escape and
        ## the post-chunk non-convergence raise take over.
        err = jnp.where(jnp.isfinite(err), err, jnp.inf)
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
                j_ra2, j_rb2 = pcnr_meta[0], pcnr_meta[1]
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
        _ra0, _rb0 = pcnr_meta[0], pcnr_meta[1]
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

## ---------------------------------------------------------------------------
## VECTOR PCNR (roadmap sec. 49): the DEVICE-shaped path.
##
## `_junction_arrays` below is the pnj-only PAIR view, and it keeps its own
## consumer -- the gmin ladders put a conductance across each junction pair on
## the ordinary `pcnr=False` path, where a gmin across a FET's `vgs` would be a
## gate leak.  Vector PCNR's unit is the DEVICE: one unknown per limited
## quantity, `m` of them, of mixed kind.  Measured on a MosLevel1Hdl,
## `pcnr_devices` reports 1 device with m=3 where `pcnr_junction_pairs`
## reports 2 pairs -- different objects, different laws -- which is the
## distinction Stage 2 turned on when `pcnr=True` on a MOSFET differential
## pair fell through to the ordinary solver in silence.
## ---------------------------------------------------------------------------


class PcnrVectorDevice(NamedTuple):
    """One participating device, in the shapes a trace can use."""
    instance: Any
    element: Any
    rows: Any            # int32 (t,)   -- circuit rows, elementnodemap order
    probe_a: Any         # int32 (m,)   -- circuit row of each probe's +
    probe_b: Any         # int32 (m,)   -- and its -
    off: int             # offset into the flat v_lim
    m: int
    t: int
    params: Any
    lin_G: Any           # sec. 40's affine remainder, or None


def _device_arrays(circuit, epar):
    """Static per-device metadata for the traced vector path, or None.

    Returns None when no device declares `pcnr_probes`, so a circuit of
    scalar-protocol junctions falls through to `_junction_arrays` and P19's
    ported path unchanged -- Stage 3's G4.

    Refusals happen HERE, at setup, with the device named, rather than
    somewhere inside a jitted loop:

    * a REDUNDANT probe (`MosLevel1Hdl`) -- `pcnr.limit_block`'s rule for an
      over-determined set is a data-dependent sort, which a traced loop
      cannot host.  `pcnr_limit_branchless` refuses it too; refusing at setup
      is what makes the message arrive before the compile.
    * a device with no `pcnr_limit_branchless` at all -- it would have no law
      on this backend.
    """
    from pycircuit.circuit.pcnr import pcnr_devices
    devs = [d for d in pcnr_devices(circuit) if not d.scalar]
    if not devs:
        return None
    out = []
    for dev in devs:
        el = dev.element
        if dev.redundant:
            raise NotImplementedError(
                "vector PCNR on the JAX backend cannot take %r (%s): it "
                "declares a REDUNDANT probe, and `pcnr.limit_block`'s rule "
                "for an over-determined set -- decreasing correction, drop "
                "what closes a cycle -- is a data-dependent sort a traced "
                "loop cannot host. Use Transient for this circuit."
                % (dev.instance, type(el).__name__))
        if getattr(el, 'pcnr_limit_branchless', None) is None:
            raise NotImplementedError(
                "vector PCNR on the JAX backend cannot take %r (%s): it "
                "declares `pcnr_probes` but no traced limiter "
                "(`pcnr_limit_branchless`), so it would have no law here."
                % (dev.instance, type(el).__name__))
        lin = getattr(el, 'pcnr_lin_G', None)
        lin_G = None if lin is None else jnp.asarray(
            np.asarray(lin(dev.params(), epar, el.toolkit), dtype=float))
        out.append(PcnrVectorDevice(
            instance=dev.instance, element=el,
            rows=jnp.asarray(dev.rows, dtype=jnp.int32),
            probe_a=jnp.asarray([a for a, _b in dev.probes], dtype=jnp.int32),
            probe_b=jnp.asarray([b for _a, b in dev.probes], dtype=jnp.int32),
            off=dev.off, m=dev.m, t=dev.t, params=dev.params(), lin_G=lin_G))
    return out


def _shadow_participants(devices):
    """Zero (or affine-shadow) each participant for the ordinary assembly.

    The CPU does this per Newton call because it must; here the assembly is
    traced ONCE, so the shadow is applied while the body is traced and never
    again.  Same rule, a fraction of the work.

    sec. 40's affine remainder: a participant whose current carries a CONSTANT
    conductance -- a series resistance to an internal node -- has that part
    shadowed IN rather than zeroed, because it reaches the solution through
    node voltages, not through any probe, and Newton solves a linear term
    exactly.
    """
    saved = []
    for d in devices:
        el = d.element
        saved.append((el, el.__dict__.get('i', None), el.__dict__.get('G', None)))
        t, lin = d.t, d.lin_G
        if lin is None:
            el.i = lambda xx, epar=None, _t=t: jnp.zeros(_t)
            el.G = lambda xx, epar=None, _t=t: jnp.zeros((_t, _t))
        else:
            el.i = lambda xx, epar=None, _g=lin: _g @ jnp.asarray(xx)
            el.G = lambda xx, epar=None, _g=lin: _g
    return saved


def _restore_participants(saved):
    for el, old_i, old_G in saved:
        for name, old in (('i', old_i), ('G', old_G)):
            if old is None:
                el.__dict__.pop(name, None)
            else:
                el.__dict__[name] = old


def pcnr_vector_blocks(circuit, devices, x, v_lim, epar, n,
                       u_extra=None, J_extra=None):
    """The augmented system's blocks, traced -- `pcnr.augmented_system`'s twin.

    Returns ``(g_mna, g_lim, J_mm, J_ml, J_lm)``.  ``J_ll`` is the identity by
    construction and is not formed, exactly as on the CPU.

    The per-device loop is PYTHON, run at trace time: the device set is static
    (it comes from the circuit, not from the state), so the traced graph gets
    one block per device and the shapes never have to be padded to a common
    `m`.  That costs graph size on a circuit with very many participants,
    which is a recorded trade rather than a surprise.
    """
    k = sum(d.m for d in devices)
    saved = _shadow_participants(devices)
    try:
        g_mna = jnp.asarray(circuit.i(x, epar))
        J_mm = jnp.asarray(circuit.G(x, epar))
    finally:
        _restore_participants(saved)
    if u_extra is not None:
        g_mna = g_mna + jnp.asarray(u_extra)
    if J_extra is not None:
        J_mm = J_mm + jnp.asarray(J_extra)

    J_ml = jnp.zeros((n, k))
    J_lm = jnp.zeros((k, n))
    g_lim = jnp.zeros(k)
    for d in devices:
        v = jax.lax.dynamic_slice(v_lim, (d.off,), (d.m,))
        el = d.element
        i_dev = jnp.asarray(el.pcnr_i(v, d.params, epar, el.toolkit))
        block = jnp.asarray(el.pcnr_didv(v, d.params, epar, el.toolkit))
        ## The device's current now enters through its OWN unknowns, so its
        ## contribution to J_MNA/MNA is zero and all of it lands in J_MNA/lim.
        g_mna = g_mna.at[d.rows].add(i_dev)
        cols = jnp.arange(d.off, d.off + d.m)
        J_ml = J_ml.at[jnp.ix_(d.rows, cols)].add(block)
        ## g_lim = v - (e_a - e_b)
        g_lim = jax.lax.dynamic_update_slice(
            g_lim, v - (x[d.probe_a] - x[d.probe_b]), (d.off,))
        ## ACCUMULATED, not assigned: a diode-connected transistor's probe can
        ## span one row (ra == rb), whose incidence is -1 + 1 = 0.
        J_lm = J_lm.at[cols, d.probe_a].add(-1.0)
        J_lm = J_lm.at[cols, d.probe_b].add(+1.0)
    return g_mna, g_lim, J_mm, J_ml, J_lm


class PcnrVectorState(NamedTuple):
    x: Any
    v_lim: Any
    iters: Any
    converged: Any


def pcnr_vector_inner_loop(circuit, devices, x_init, epar, n, irefnode,
                           reltol=1e-4, abstol=1e-12, xtol=1e-12,
                           maxiter=100, u_extra=None, J_extra=None):
    """Newton on ``[x_MNA ; v_lim]``, device-shaped, inside a traced loop.

    The CPU's `pcnr.solve_dc` structure exactly -- assemble the augmented
    blocks, collapse them by the Schur identity `pcnr.schur_reduce` writes
    once, solve the reduced MNA system, recover the limited unknowns from
    ``J_ll = I``, then CORRECT each device's own block with its own law.
    """
    def v_from_x(x):
        if not devices:
            return jnp.zeros(0)
        return jnp.concatenate([x[d.probe_a] - x[d.probe_b] for d in devices])

    def body(ps: PcnrVectorState):
        x, v_lim = ps.x, ps.v_lim
        g_mna, g_lim, J_mm, J_ml, J_lm = pcnr_vector_blocks(
            circuit, devices, x, v_lim, epar, n, u_extra, J_extra)

        ## `pcnr.schur_reduce`'s formula, and the only place it is written
        ## on this backend: J_ll is the identity by construction, so the
        ## border collapses to a rank-k update.
        J_eff = J_mm - J_ml @ J_lm
        F_eff = g_mna - J_ml @ g_lim

        J_sub = jnp.delete(jnp.delete(J_eff, irefnode, axis=0),
                           irefnode, axis=1)
        dx_sub = jnp.linalg.solve(J_sub, -jnp.delete(F_eff, irefnode))
        dx = jnp.insert(dx_sub, irefnode, 0.0)
        x_new = x + dx

        ## J_ll = I, so the limited unknowns follow from the MNA step.
        v_raw = v_lim - (g_lim + J_lm @ dx)

        ## CORRECT: each device limits only the variables it owns, with its
        ## OWN law, handed its whole block at once and its own slice of the
        ## last accepted iterate (a parameter that follows the bias --
        ## SPICE's `von` -- is evaluated there).  Nothing is shared, so one
        ## device's limiter cannot disturb another's; that is the property
        ## PCNR buys and the reason there is no ordering here.
        if devices:
            blocks = []
            for d in devices:
                el = d.element
                blocks.append(jnp.asarray(el.pcnr_limit_branchless(
                    jax.lax.dynamic_slice(v_raw, (d.off,), (d.m,)),
                    jax.lax.dynamic_slice(v_lim, (d.off,), (d.m,)),
                    d.params, epar, el.toolkit, x[d.rows])))
            v_new = jnp.concatenate(blocks)
        else:                                        # pragma: no cover
            v_new = v_raw

        I_scale = jnp.abs(J_eff) @ jnp.abs(x_new) + jnp.abs(F_eff)
        ## THE REDUCED ROWS: the reference row is not an equation of the
        ## solved system, and scoring it can never be satisfied.
        conv_f = jnp.all(jnp.delete(
            jnp.abs(F_eff) < reltol * I_scale + abstol, irefnode))
        conv_x = jnp.all(jnp.abs(dx)
                         < reltol * jnp.maximum(jnp.abs(x_new), jnp.abs(x))
                         + xtol)
        ## ⚠ PER COMPONENT, which is `pcnr.lim_converged`'s rule and NOT
        ## what the scalar traced loop above does.  Against `max|v_new|`
        ## over the whole vector, a 40 V `vds` sharing a device with a
        ## `vbe` would loosen the `vbe` component forty-fold.  For a
        ## circuit of scalar diodes the two criteria coincide, which is
        ## why the older loop is right for its own case and wrong here.
        lim_ok = jnp.all(jnp.abs(g_lim)
                         < reltol * jnp.maximum(jnp.abs(v_new), 1.0) + abstol)
        conv = jnp.logical_and(jnp.logical_and(conv_x, conv_f), lim_ok)
        return PcnrVectorState(x=x_new, v_lim=v_new, iters=ps.iters + 1,
                               converged=conv)

    def cond(ps: PcnrVectorState):
        return jnp.logical_and(jnp.logical_not(ps.converged),
                               ps.iters < maxiter)

    init = PcnrVectorState(x=x_init, v_lim=v_from_x(x_init),
                           iters=jnp.asarray(0),
                           converged=jnp.asarray(False))
    return jax.lax.while_loop(cond, body, init)


class PcnrState(NamedTuple):
    x: Any
    v_lim: Any
    iters: Any
    converged: Any


def _junction_arrays(circuit):
    """Static junction metadata for the traced PCNR path.

    Returns ``(ra, rb, IS, VT, funcs)``.  ``funcs`` is a list of
    ``(i_fn, g_fn)`` pairs -- one per junction, each a traced callable of
    the junction's own voltage -- or ``None`` when some device cannot
    supply one, in which case the caller falls back to rebuilding the
    textbook diode from ``IS`` and refuses anything that is not one.

    The traced loop used to rebuild EVERY junction itself as
    ``IS*(exp(v/VT)-1)`` with a single global VT, reading ``IS`` by name.
    That silently gave a device whose saturation current is not called
    ``IS`` a junction carrying no current, and it made multi-junction and
    charge-storing participants impossible here while the CPU path
    accepted them.  Asking the device instead removes all three limits at
    once: any shape, any number of junctions, and charge simply stays in
    the MNA block exactly as it does on the CPU (this assembly subtracts
    the junction term at the node voltage and adds it at ``v_lim``; it
    never touched the charge).
    """
    from pycircuit.circuit.pcnr import pcnr_junctions
    from pycircuit.circuit.circuit import defaultepar as _dp
    junctions = pcnr_junctions(circuit)
    if not junctions:
        return None

    VT_global = float(_dp.kboltzmann if hasattr(_dp, 'kboltzmann')
                      else 1.38e-23) * float(getattr(_dp, 'T', 300.0)) \
        / 1.602e-19

    ra = jnp.array([j[2] for j in junctions], dtype=jnp.int32)
    rb = jnp.array([j[3] for j in junctions], dtype=jnp.int32)

    IS_list, VT_list, funcs, seen = [], [], [], {}
    generic = True
    for inst, element, _ra, _rb in junctions:
        jn = seen.get(inst, 0)
        seen[inst] = jn + 1
        params = {q.name: getattr(element.iparv, q.name)
                  for q in element.instparams}
        ## Per-junction (VT, IS): what pnjlim needs.  A device that knows
        ## its own scales says so; otherwise fall back to the `IS`
        ## parameter and the global thermal voltage.
        scales = getattr(element, 'pcnr_scales', None)
        if scales is not None:
            VT_j, IS_j = scales(params, _dp, jn=jn)
        else:
            VT_j, IS_j = VT_global, float(getattr(element.iparv, 'IS', 0.0))
        VT_list.append(VT_j)
        IS_list.append(IS_j)

        pi = getattr(element, 'pcnr_i', None)
        pg = getattr(element, 'pcnr_didv', None)
        if pi is None or pg is None:
            generic = False
            continue
        tk = element.toolkit

        def _mk(pi=pi, pg=pg, params=params, tk=tk, jn=jn):
            def i_fn(v):
                return pi(v, params, _dp, tk, jn=jn)[0]

            def g_fn(v):
                return pg(v, params, _dp, tk, jn=jn)[0]
            return i_fn, g_fn
        funcs.append(_mk())

    if generic and funcs:
        ## Probe once, outside the trace, so a device that cannot be
        ## evaluated at all fails here with its own name rather than
        ## somewhere inside a jitted loop.
        for k, (i_fn, _g) in enumerate(funcs):
            try:
                float(np.asarray(i_fn(0.3), dtype=float))
            except Exception as exc:
                raise NotImplementedError(
                    'junction %d of %r could not be evaluated for the '
                    'traced backend (%s); use the CPU backend or '
                    'pcnr=False.' % (k, junctions[k][0], exc)) from exc
        return (ra, rb, jnp.array(IS_list), jnp.array(VT_list), funcs)

    ## FALLBACK: no device-supplied evaluation, so the traced loop must
    ## rebuild the textbook diode -- and may only do so for devices that
    ## really are one.  Everything else is refused rather than quietly
    ## computed wrong.
    for k, (inst, element, _ra, _rb) in enumerate(junctions):
        _has_charge = False
        _probes = (np.linspace(0.1, 0.4, element.n),
                   np.linspace(-0.2, 0.5, element.n))
        try:
            for _xp in _probes:
                if np.any(np.abs(np.asarray(element.C(_xp),
                                            dtype=float)) > 0.0):
                    _has_charge = True
                    break
        except Exception:
            _has_charge = True
        if _has_charge:
            raise NotImplementedError(
                '%r stores charge and declares a PCNR junction, and does '
                'not supply pcnr_i/pcnr_didv for the traced backend, '
                'which would have to rebuild the junction itself.' % inst)
        pcnr_i = getattr(element, 'pcnr_i', None)
        if pcnr_i is None:
            continue
        params = {q.name: getattr(element.iparv, q.name)
                  for q in element.instparams}
        for v_probe in (0.3, 0.6):
            got = float(np.asarray(pcnr_i(v_probe, params, _dp,
                                          element.toolkit), dtype=float)[0])
            want = IS_list[k] * (np.exp(v_probe / VT_global) - 1.0)
            if not np.isclose(got, want, rtol=1e-6,
                              atol=1e-30 + 1e-6 * abs(want)):
                raise NotImplementedError(
                    "%r declares a PCNR junction the traced backend cannot "
                    "reproduce (it gives %.6g A at %.2f V where "
                    "IS*(exp(v/VT)-1) gives %.6g A)." % (inst, got,
                                                         v_probe, want))
    return (ra, rb, jnp.array(IS_list), jnp.array(VT_list), None)


def _junction_terms_generic(v, funcs):
    """Ask each device for its own junction current and conductance.

    A Python loop over a STATIC list, so `jit` unrolls it; the junction
    set is fixed at trace time."""
    i = jnp.stack([f_i(v[k]) for k, (f_i, _g) in enumerate(funcs)])
    g = jnp.stack([f_g(v[k]) for k, (_i, f_g) in enumerate(funcs)])
    return i, g


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
                    j_IS, VT, tline_params, tline_indices, j_VT=None,
                    j_fns=None,
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
        if j_fns is not None:
            i_n, g_n = _junction_terms_generic(v_node, j_fns)
            i_l, g_l = _junction_terms_generic(v_lim, j_fns)
        else:
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
        v_new = _pnjlim_branchless(v_lim + dx_lim, v_lim,
                                   VT if j_VT is None else j_VT, j_IS,
                                   jnp.log)

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
                             tline_dG=None, j_fns=None):
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
    if j_fns is not None:
        _i_n, g_n = _junction_terms_generic(v_node, j_fns)
        _i_l, g_l = _junction_terms_generic(v_lim, j_fns)
    else:
        _i_n, g_n = _junction_terms(v_node, j_IS, VT)
        _i_l, g_l = _junction_terms(v_lim, j_IS, VT)
    J = _scatter_junction_G(J, j_ra, j_rb, g_n, -1.0)
    return _scatter_junction_G(J, j_ra, j_rb, g_l, +1.0)


def outer_time_loop(initial_state: TransientState, circuit, tend, chunk_size, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method='euler', params_tree=None, reltol=1e-4, abstol=1e-12, xtol=1e-12, maxiter=100, trtol=7.0, lte_reltol=1e-4, lte_abstol=1e-12, max_dv=jnp.inf, coupled=False, gamma_min=0.7, gamma_max=3.0, eta=0.15, pcnr_meta=None, pcnr_VT=0.025, tline_dG=None, analysis='tran', s_gamma_min=0.0, s_gamma_max=1.0, s_eta=jnp.inf, fixed_timestep=False, grid_dt=None, relref='sigglobal', n_nodes=None, provided_function=None, state_mask=None, dv_bounds=((jnp.inf, 0.0), (jnp.inf, 0.0)), periodic_states=None):

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
        ## A forced NON-converged accept is already an unconditional raise
        ## after the chunk (stage 9(e)) -- no run continues past one -- so
        ## the chunk exits the moment it happens instead of marching to
        ## chunk_size at the dt floor first.  P22's mask port measured why
        ## this matters: with the algebraic rows masked, a cold-start
        ## circuit whose Newton fails at every h reached the floor honestly
        ## and then ground out 500 forced steps at ~1 s each on GPU (the
        ## coupled point is ~100 serial kernel launches per attempt) before
        ## the raise -- an effective hang.  Pre-mask, the algebraic-row
        ## error was ACCIDENTALLY load-bearing as a step governor that
        ## rescued the Newton; the mask removed the accident, this exit
        ## replaces it with the deliberate escape.
        alive = state.n_forced_nonconverged == 0
        return jnp.logical_and(jnp.logical_and(under_time, under_chunk),
                               alive)

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
                provided_function=provided_function, state_mask=state_mask)
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
            j_ra, j_rb, j_IS, j_VT, j_fns = pcnr_meta
            pstate = pcnr_inner_loop(
                state, circuit, irefnode, j_ra, j_rb, j_IS, pcnr_VT,
                tline_params, tline_indices, j_VT=j_VT, j_fns=j_fns,
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
            j_ra, j_rb, j_IS, j_VT, j_fns = pcnr_meta
            J = pcnr_controller_jacobian(circuit, state, x_curr, pcnr_v_lim,
                                         j_ra, j_rb, j_IS, pcnr_VT,
                                         eval_method, first_order,
                                         params_tree=params_tree,
                                         tline_dG=tline_dG, j_fns=j_fns)
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
        ## VOLTAGE CHECK (max_dv_step, inf when disabled): node rows only,
        ## proportional retry -- the industry semantics; see the Parameter.
        ## Skipped at the dt floor (progress guarantee) and under a fixed
        ## grid (the caller owns the resolution).
        _nn_dv = n_nodes if n_nodes is not None else x_curr.shape[0]
        _dstep = jnp.abs(x_curr - state.x_history[0])
        ## 'auto' (sampling theory): bound = max(static source anchor,
        ## (2*pi/N) * running unit-group max).  ref_running is accepted-only
        ## and so anchored BEFORE this candidate -- the relative term cannot
        ## h-cancel at a signal birth.  Manual factors carry rel = 0.
        (_bvs, _cvr), (_bis, _cir) = dv_bounds
        _bv = jnp.maximum(_bvs, _cvr * jnp.max(state.ref_running[:_nn_dv])) \
            if _cvr else _bvs
        dv_ratio = jnp.max(_dstep[:_nn_dv]) / _bv
        if _nn_dv < _dstep.shape[0]:
            _bi = jnp.maximum(
                _bis, _cir * jnp.max(state.ref_running[_nn_dv:])) \
                if _cir else _bis
            dv_ratio = jnp.maximum(dv_ratio,
                                   jnp.max(_dstep[_nn_dv:]) / _bi)
        dv_ok = jnp.logical_or(dv_ratio <= 1.0,
                               jnp.asarray(fixed_timestep))
        accept = jnp.logical_or(
            jnp.logical_and(
                jnp.logical_and(nr_state.converged,
                                jnp.logical_and(lte_ok, dv_ok)),
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
                ## The SOLVED step carries forward; breakpoint truncation for
                ## the NEXT point happens at the next fang entry.  Solved is
                ## the operative word (P22's mask port measured the hole): a
                ## NON-converged fang's h is not a solved h -- with the band
                ## below tolerance it GROWS within the point, so a forced
                ## floor-accept came back at ~2x dt_min, the at_floor escape
                ## went un-sticky, and the run crawled at ~100 fang
                ## iterations per 1e-18 s -- an effective hang where the
                ## unmasked path reached its post-chunk non-convergence
                ## raise in seconds.  A forced accept keeps the floor.
                next_dt = jnp.where(nr_state.converged,
                                    jnp.clip(state.dt, dt_min, dt_max),
                                    jnp.asarray(dt_min))
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

            ## PHASE-2 GAUGE SHIFT (idtmod.md 5.2), ACCEPT-ONLY AND BRANCHLESS.
            ## Each declared periodic row is rewrapped by an exact n*modulus
            ## translation applied to the accepted point AND the whole live
            ## ring -- a uniform shift the multistep formulas and the divided-
            ## difference LTE cannot see -- so the state stays bounded (the
            ## floating-point property idtmod exists for).  `iq` is
            ## derivative-domain and invariant; sig_max/ref_running use the
            ## post-shift (bounded) values, which is itself part of the
            ## payoff: an unbounded state inflates the sigglobal reference
            ## and loosens every later tolerance in proportion to run time.
            ## The candidate step's error_ratio was computed pre-shift, in
            ## one consistent gauge.  `periodic_states=None` compiles to the
            ## identical program as before -- the branch is trace-time.
            if periodic_states is not None:
                p_rows, p_mods, p_offs = periodic_states
                _n_wraps = jnp.floor((x_curr[p_rows] - p_offs) / p_mods)
                _d = _n_wraps * p_mods
                x_acc = x_curr.at[p_rows].add(-_d)
                q_acc = q_curr.at[p_rows].add(-_d)
                x_hist_base = state.x_history.at[:, p_rows].add(-_d)
                q_hist_base = state.q_history.at[:, p_rows].add(-_d)

                ## ARC 7 (idtmod.md sec. 5.3): the in-trace wrap-crossing dt
                ## cap, which that section listed as future work because
                ## `t_breaks_array` is static and a state-dependent event
                ## cannot enter it.  It enters HERE instead: not as a
                ## breakpoint, but as a cap on the next step, computed from
                ## state the trace already holds.  Branchless -- no host
                ## control flow, so it compiles inside `lax.while_loop`.
                ##
                ## ⚠ WHAT THIS DOES AND DOES NOT BUY, measured before it was
                ## built.  A wrap is a discontinuity in the OUTPUT MAP, not in
                ## the ODE: the gauge shift above keeps the state continuous
                ## and the sawtooth is an exact `floored_wrap` of it, so there
                ## is no kink for the integrator to resolve.  Accordingly it
                ## buys NOTHING in accuracy or cost -- against a tight
                ## reference on a varying integrand, JAX and the CPU agreed to
                ## within 1-5% of each other's error (6.171e-03 vs 6.112e-03 at
                ## reltol 1e-4; 2.925e-04 vs 2.775e-04 at 1e-6) at the same
                ## step counts, and sec. 5.3's predicted step-collapse at wraps
                ## DID NOT REPRODUCE (70 points against the CPU's 74, no
                ## rejection storm).
                ##
                ## What it buys is SAMPLE PLACEMENT.  Measured on the exact
                ## ramp, distance from the true corners: CPU (which has
                ## breakpoints) 3.55e-15, JAX 3.29e-02 -- up to three output
                ## timesteps past the corner.  That matters to a consumer that
                ## RESAMPLES or edge-detects the output, because interpolating
                ## a sawtooth across an unmarked corner returns values the
                ## signal never takes (sec. 3.4's consumer discontinuity).
                ##
                ## THE PREDICTION is linear, from the step just accepted --
                ## the CPU's own rule in sec. 5.3 ("the linearly predicted
                ## crossing time").  It only has to BRACKET the corner; the
                ## step is capped, not landed exactly, and the next step
                ## re-predicts from there.
                if not fixed_timestep:
                    ## The previous ACCEPTED point is the head of the ring;
                    ## `x_curr` is pre-shift and the ring is not yet shifted,
                    ## so both are in one gauge and the slope is the step's.
                    _rate = ((x_curr[p_rows] - state.x_history[0, p_rows])
                             / state.dt)
                    ## Post-shift the state lies in [offs, offs+mods), so the
                    ## boundary being approached is the top going up and the
                    ## bottom going down.
                    _x_s = x_acc[p_rows]
                    _dist = jnp.where(_rate > 0.0,
                                      (p_offs + p_mods) - _x_s,
                                      p_offs - _x_s)
                    ## A point sitting exactly ON the boundary it is leaving
                    ## must travel a whole modulus, not zero -- otherwise the
                    ## cap collapses the next step to dt_min at every wrap,
                    ## turning an accuracy-neutral change into a cost one.
                    ##
                    ## ⚠ DEFENSIVE AND UNEXERCISED, and said plainly rather
                    ## than left to look tested.  The case needs a DESCENDING
                    ## state landing exactly on `offs`: ascending, a landing
                    ## on the top is shifted to `offs` and the next boundary
                    ## going up is a full modulus away, so no zero arises.
                    ## Descending was built and driven (v = -1 into modulus 1,
                    ## three wraps) and the branch still never fires, because
                    ## the cap lands just SHORT of the boundary in floating
                    ## point -- the gauge shift has then already wrapped the
                    ## state to the top and the distance is again a modulus.
                    ## Removing this branch changed that run not at all
                    ## (64 points, corners at 8e-15, both ways).  It is kept
                    ## because "unreachable in the arithmetic I could
                    ## construct" is not "unreachable", and the failure it
                    ## prevents is a step collapse at every wrap.
                    _eps = 1e-12 * p_mods
                    _dist = jnp.where(jnp.abs(_dist) < _eps,
                                      jnp.where(_rate > 0.0, p_mods, -p_mods),
                                      _dist)
                    ## Rate ~ 0 means no crossing is predicted at all.
                    _moving = jnp.abs(_rate) > _eps
                    _dt_cross = jnp.where(
                        _moving, _dist / jnp.where(_moving, _rate, 1.0),
                        jnp.inf)
                    _wrap_dt = jnp.min(_dt_cross, initial=jnp.inf)
                    next_dt = jnp.maximum(jnp.minimum(next_dt, _wrap_dt),
                                          dt_min)
            else:
                x_acc, q_acc = x_curr, q_curr
                x_hist_base = state.x_history
                q_hist_base = state.q_history

            x_hist_new = jnp.roll(x_hist_base, shift=1, axis=0).at[0].set(x_acc)
            q_hist_new = jnp.roll(q_hist_base, shift=1, axis=0).at[0].set(q_acc)
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

            res_buf = state.results_buffer.at[state.step_idx].set(x_acc)
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
                sig_max=jnp.maximum(state.sig_max, jnp.max(jnp.abs(x_acc))),
                ## P7: accepted-only, like sig_max; `local` matches the
                ## estimator's (current point vs the one before it).
                ## Post-shift values on both: bounded periodic states must
                ## not inflate the running references (see the gauge-shift
                ## note above).
                ref_running=jnp.maximum(
                    state.ref_running,
                    jnp.maximum(jnp.abs(x_acc),
                                jnp.abs(x_hist_base[0]))),
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
            ## A dv-vetoed step retries at 0.9*bound/dv (industry semantics);
            ## a converged-but-LTE-rejected one keeps the predictor's cut.
            dv_shrink = jnp.clip(0.9 / jnp.maximum(dv_ratio, 1e-300),
                                 0.1, 0.9)
            lte_shrink = jnp.where(jnp.logical_not(dv_ok), dv_shrink,
                                   lte_shrink)
            shrink = jnp.where(nr_state.converged, lte_shrink, 0.25)
            retry_dt = jnp.maximum(state.dt * shrink, dt_min)
            if coupled:
                ## The CPU coupled path's outer retry: h *= 0.25 whatever the
                ## failure flavour (Newton or held-step LTE) -- applied to
                ## the ENTRY h, not to fang's returned h.  fang GROWS h
                ## within the point when the (masked) band reads the error
                ## as tiny, by up to MAX_GROWTH=5x; shrinking the grown h by
                ## 4x nets a 1.25x GROWTH per reject cycle, so a cold-start
                ## circuit whose Newton fails at every h never reached the
                ## dt floor where the forced-accept escape lives -- an
                ## in-trace livelock the CPU cannot express (its retry loop
                ## is Python-bounded at 10 attempts).  Measured behind
                ## P22's mask port; the mask only made it reachable.
                retry_dt = jnp.maximum(
                    jnp.minimum(state.dt, h_entry) * 0.25, dt_min)
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
    """Fully compiled transient analysis: the whole time loop is one traced
    ``jax.lax.while_loop``, jitted per chunk, with adaptive LTE control and
    the same parameter vocabulary as `Transient`.

    **Backend parity** (doc/backend_parity_260821.md is the ledger; the
    conformance harness in test_backend_conformance.py pins the agreement):

    * **Shared and aligned**: tolerances, `TRTOL`, `timestep_max` (None means
      tend/50), `minstep`, `minbreak`, `outputstep`, `uic`/`ic` (the CPU's
      initial-state machinery runs here verbatim, spanning-tree capacitor
      solve included), `relref` with the unit-group split, the standard LTE
      acceptance band (`lte_gamma_*`/`lte_eta`), `fixed_timestep` (bit-equal
      fixed-grid waveforms across backends), `provided_function` (must be
      jax-traceable; baked in at jit time), and the `integrator` choice --
      'gear' (default, matching the CPU's Gear2 default) or 'euler'.
      Research paths run on both backends too: `coupled_lte` (Fang),
      `pcnr`, PCNR-inside-Fang, and TLine on every path.

    * **P13, `bypass`**: a non-concept here, not a missing feature -- device
      bypassing skips evaluating quiescent elements, and a vmapped
      evaluation group computes all lanes of all instances in one kernel;
      there is nothing to skip.

    * **P17, strategy objects**: `nrsolver`/`scaler`/`linearsolver` are
      refused in __init__, PERMANENTLY -- a traced while_loop cannot
      dispatch into per-iteration Python objects.  Use `Transient` for a
      circuit that needs them.

    * **P20, `solve_batched`**: this backend's purpose -- one compiled
      kernel integrating every parameter lane of a sweep concurrently, each
      lane with its own adaptive (or coupled (x, h)) step sequence.  The
      CPU deliberately has no imitation of it.

    * **CPU-only, with cause**: trapezoidal integration (the uniform-grid
      trap branch was deleted; a variable-step trap estimator has not been
      written), the coupled 'bordered' branch, and coupled+fixed_timestep
      (grid_locked is not wired here yet -- refused, not approximated).
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
        ## The Spectre/Mica-class VOLTAGE CHECK -- see the CPU's Parameter
        ## note: on algebraic networks no LTE exists (the charge estimator
        ## is identically zero; P22's mask excludes algebraic rows from the
        ## coupled band), and this bounds the per-step node-voltage change
        ## instead.  |dv| ~ h*slew is h-proportional, so it cannot h-cancel.
        Parameter(name='max_dv_step',
                  desc='Per-step node-voltage excursion bound (the Spectre/'
                       'Mica-style voltage check), as a FACTOR: the bound '
                       'is max_dv_step * lte_vabstol (e.g. 2e11 at the '
                       'default 1e-12 bounds steps to 0.2 V), scaling with the '
                       'tolerance family as the LTE does; factors below 1 '
                       'clamp to 1. None disables.',
                  unit='', default=None),
        Parameter(name='max_di_step',
                  desc='Per-step branch-current excursion bound, as a '
                       'FACTOR times lte_iabstol; same clamp-at-1 floor. '
                       "'auto' derives it from sampling theory "
                       '(points_per_period). None disables.',
                  unit='', default=None),
        Parameter(name='points_per_period',
                  desc="Sampling density behind max_dv_step/max_di_step = "
                       "'auto': at least this many points per period of a "
                       'full-swing sinusoid (per-step excursion '
                       '2*pi*swing/N).',
                  unit='', default=64),
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

        ## P17 -- and this refusal is the PERMANENT contract, not a "not
        ## yet": `nrsolver`/`scaler`/`linearsolver` are Python strategy
        ## objects dispatched per iteration, and a traced `jax.lax.while_loop`
        ## cannot call into them -- ever; that is what tracing means.  They
        ## used to be accepted in silence (the "thin advertised feature" 0.1c
        ## warns about, the same defect shape as the `lte_formula` knob
        ## removed in 9(f)); `linearsolver` joined the refusal at Phase C --
        ## the traced loop solves with jnp.linalg.solve and a passed solver
        ## object was still being swallowed.
        for unsupported in ('nrsolver', 'scaler', 'linearsolver'):
            if kwargs.get(unsupported) is not None:
                raise NotImplementedError(
                    "JAXTransient does not support %r, permanently: its Newton "
                    "loop is a traced jax.lax.while_loop, and a traced loop "
                    "cannot dispatch into per-iteration Python strategy "
                    "objects. Use Transient for a circuit that needs %r."
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

    ## P23 sharing: the CPU's excursion-bound resolver and source-scale
    ## walker run here verbatim (pre-loop Python; the loop takes floats).
    from pycircuit.circuit.transient import Transient as _CPU2
    _source_signal_scales = _CPU2._source_signal_scales
    _dv_step_bounds = _CPU2._dv_step_bounds
    del _CPU2

    def _jax_state_row_mask(self, x_ref):
        """P22: the CPU's _state_row_mask, computed numpy-side at trace
        prep -- True where the unknown participates in any charge; algebraic
        rows get a 1e30 tolerance in fang's band so eq (6) measures states
        only.  Structural at the seed point, same caveat as the CPU's."""
        C = np.abs(np.asarray(self.cir.C(jnp.asarray(np.asarray(x_ref,
                                                                dtype=float)))))
        return jnp.asarray((C.sum(axis=0) + C.sum(axis=1)) > 0.0)

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

    def _periodic_state_arrays(self):
        """Static (rows, moduli, offsets) arrays for the Phase-2 gauge shift,
        or ``None`` when nothing declares (idtmod.md 5.2)."""
        if not hasattr(self.cir, 'periodic_states'):
            return None
        entries = list(self.cir.periodic_states())
        if not entries:
            return None
        return (jnp.array([int(e[0]) for e in entries], dtype=jnp.int32),
                jnp.array([float(e[1]) for e in entries]),
                jnp.array([float(e[2]) for e in entries]))

    def _periodic_state_arrays_batched(self, override_params_tree,
                                       batch_size):
        """Per-lane (rows, moduli, offsets) for the Phase-2 shift under
        ``solve_batched`` -- roadmap item 3 (idtmod.md sec. 8).

        A swept ``modulus``/``offset`` used to DROP the element from the
        shift (correct Phase-1 fallback, unbounded state); now the
        declaration is re-evaluated per lane with the lane's parameter
        column applied, so swept lanes KEEP the bounded-state property.
        Pre-trace numpy work: swept top-level instances get their
        parameters set per lane, ``periodic_states()`` is re-collected,
        and the originals are restored in ``finally``.  Column order
        within a class override is the group (elements-iteration) order,
        same as ``batched_contributions`` consumes.

        Lanes must agree on WHICH rows are periodic -- a modulus swept
        across the idt-degradation boundary (positive to non-positive)
        changes the declaration set -- and a disagreeing element's rows
        are dropped to the Phase-1 fallback with a warning rather than
        mis-shifted.

        Returns ``(rows int32 (n,), mods (batch, n), offs (batch, n))``
        or ``None`` when nothing declares.
        """
        if not hasattr(self.cir, 'periodic_states'):
            return None
        elements = getattr(self.cir, 'elements', {})

        ## Which top-level declaring instances are swept, and their column
        ## in the class's override arrays.
        swept = []      # (inst_name, element, {param: (batch,) values})
        for inst_name, element in elements.items():
            if not hasattr(element, 'periodic_states'):
                continue
            cls_name = type(element).__name__
            if cls_name not in (override_params_tree or {}):
                continue
            siblings = [nm for nm, el in elements.items()
                        if type(el).__name__ == cls_name]
            col = siblings.index(inst_name)
            lane_params = {}
            for key, arr in override_params_tree[cls_name].items():
                if key not in getattr(element, 'iparv', ()):
                    continue
                a = np.atleast_2d(np.asarray(arr))
                lane_params[key] = a[:, col]
            if lane_params:
                swept.append((inst_name, element, lane_params))

        base = list(self.cir.periodic_states())
        if not base and not swept:
            return None

        if not swept:
            if not base:
                return None
            rows = [int(e[0]) for e in base]
            mods = np.tile([float(e[1]) for e in base], (batch_size, 1))
            offs = np.tile([float(e[2]) for e in base], (batch_size, 1))
            return (jnp.array(rows, dtype=jnp.int32), jnp.array(mods),
                    jnp.array(offs))

        ## Re-evaluate the swept instances' declarations per lane.
        swept_rows_of = {}   # inst_name -> [global rows] (lane-0 reference)
        lane_decls = []      # per lane: {global_row: (m, o)}
        saved = {name: dict(el.iparv._values)
                 for name, el, _ in swept}
        try:
            for lane in range(batch_size):
                decls = {}
                for inst_name, element, lane_params in swept:
                    element.iparv.set(**{k: float(v[lane])
                                         for k, v in lane_params.items()})
                    rows_map = self.cir.elementnodemap[inst_name]
                    rows_here = []
                    for local_row, m, o in element.periodic_states():
                        g = int(rows_map[local_row])
                        decls[g] = (float(m), float(o))
                        rows_here.append(g)
                    if lane == 0:
                        swept_rows_of[inst_name] = rows_here
                    elif rows_here != swept_rows_of[inst_name]:
                        swept_rows_of[inst_name] = None   # inconsistent
                lane_decls.append(decls)
        finally:
            for name, el, _ in swept:
                el.iparv.set(**saved[name])

        dropped = [nm for nm, rows in swept_rows_of.items() if rows is None]
        if dropped:
            warnings.warn(
                'solve_batched: %s declare periodic states on some lanes '
                'but not others (a modulus swept across the degradation '
                'boundary?); the gauge shift is disabled for them and their '
                'state runs unbounded (Phase-1 behaviour, still correct).'
                % ', '.join(sorted(repr(n) for n in dropped)),
                RuntimeWarning)
        swept_rows = set()
        for nm, rows in swept_rows_of.items():
            if rows is not None:
                swept_rows.update(rows)

        ## Constant (unswept) rows come from the base declaration; the rows
        ## of DROPPED instances must go too -- their trace-time modulus is
        ## wrong for the swept lanes.
        base_rows = {int(e[0]): (float(e[1]), float(e[2])) for e in base}
        all_dropped_rows = set()
        for nm in dropped:
            rows_map = self.cir.elementnodemap[nm]
            el = elements[nm]
            for local_row, _m, _o in el.periodic_states():
                all_dropped_rows.add(int(rows_map[local_row]))
        rows = sorted((set(base_rows) - all_dropped_rows) | swept_rows)
        if not rows:
            return None

        mods = np.empty((batch_size, len(rows)))
        offs = np.empty((batch_size, len(rows)))
        for j, row in enumerate(rows):
            for lane in range(batch_size):
                if row in lane_decls[lane]:
                    mods[lane, j], offs[lane, j] = lane_decls[lane][row]
                else:
                    mods[lane, j], offs[lane, j] = base_rows[row]
        return (jnp.array(rows, dtype=jnp.int32), jnp.array(mods),
                jnp.array(offs))

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
        ## THE TREE IS KEYED BY CLASS NAME -- the vmapped groups are
        ## per-class stacks, and `batched_contributions` consumes
        ## `params_tree[cls.__name__]`.  Every test happened to name its
        ## instance after its class (`c['R'] = R(...)`), which hid this from
        ## users until an instance named 'R1' was refused with a message
        ## about classes (roadmap item 1, idtmod.md sec. 8).  An INSTANCE
        ## key is remapped when unambiguous -- its class has exactly one
        ## instance -- and refused with the correct spelling otherwise,
        ## because with several instances the override's per-instance column
        ## order is the group order, which the caller must address by class.
        elements = getattr(self.cir, 'elements', {})
        remapped = {}
        for key, value in override_params_tree.items():
            if key in batchable:
                remapped[key] = value
                continue
            element = elements.get(key)
            cls_name = type(element).__name__ if element is not None else None
            if cls_name in batchable:
                siblings = [nm for nm, el in elements.items()
                            if type(el).__name__ == cls_name]
                if len(siblings) == 1:
                    remapped[cls_name] = value
                    continue
                raise NotImplementedError(
                    "override_params_tree names instance %r, but the tree is "
                    "keyed by CLASS name and class %r has %d instances (%s): "
                    "use {%r: ...} with one column per instance, in that "
                    "order." % (key, cls_name, len(siblings),
                                ', '.join(siblings), cls_name))
            remapped[key] = value
        override_params_tree = remapped
        unbatchable = set(override_params_tree) - batchable
        if unbatchable:
            raise NotImplementedError(
                "override_params_tree names %s, which cannot be batched: only "
                "element classes providing eval_i_pure/eval_q_pure participate "
                "in vmapped evaluation (currently: %s), and the tree is keyed "
                "by CLASS name. An override for any other class would be "
                "silently ignored and every lane would return the same "
                "waveform. Make the class batchable, or drop the override."
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
        ## P21: the batched DC operating point, retiring the roadmap's last
        ## refusal.  Each lane's bias is solved with ITS OWN swept
        ## parameters (that is the whole point -- an R sweep moves the bias
        ## per lane); a lane whose plain Newton cannot converge at DC
        ## raises with the honest alternatives, because a batched run
        ## seeded from a non-solution is not an approximation of the right
        ## answer, it is a different problem.
        _dc_x0 = None
        if not uic:
            _n = self.cir.n
            _iref = self.cir.get_node_index(refnode) \
                if hasattr(self.cir, 'get_node_index') else refnode

            _gmin_rows = _gmin_junction_rows(self.cir)

            def _lane_dc(p_tree):
                return dc_with_continuation(
                    self.cir, _iref, _n, len(self.cir.nodes),
                    params_tree=p_tree,
                    reltol=self.par.reltol,
                    abstol=self._newton_abstol(_n),
                    xtol=self._newton_xtol(_n),
                    maxiter=int(self.par.maxiter),
                    gmin_rows=_gmin_rows)

            ## Roadmap item 2 (idtmod.md sec. 8): the ic-pin flag for the
            ## traced DC.  The traced assemblies call circuit.i/G without an
            ## epar, so elements read `defaultepar` -- the flag goes there,
            ## and it is baked in at TRACE time (the jitted closure below is
            ## fresh per call, so nothing stale is reused; the batched-eval
            ## cache keys on epar values and separates the two states).
            from pycircuit.circuit.circuit import defaultepar
            from pycircuit.circuit.analysis import analysis_kind
            with analysis_kind(defaultepar, 'dc'):
                _dc_x0, _dc_conv = jax.jit(jax.vmap(_lane_dc))(
                    override_params_tree)
            _bad = np.where(~np.asarray(_dc_conv))[0]
            if _bad.size:
                raise NoConvergenceError(
                    'solve_batched: the DC operating point did not converge '
                    'for lane(s) %s, even through the gmin, gshunt and '
                    'pseudo-transient continuation ladders (P18/P25). The '
                    'circuit likely has no DC solution on these lanes (a '
                    'floating node with no DC path stays singular at g=0). '
                    'Start from uic=True with explicit ic if the start '
                    'state is known.'
                    % _bad.tolist())


        if dt_max is None:
            ## STAGE 9(g) -- `timestep`, not `tend/10`.  `solve` and `solve_batched`
            ## are two entry points of one class and disagreed about this by ~50x,
            ## so the same circuit at the same requested timestep was error-
            ## controlled to two different standards depending on which was called.
            ## `tend/50` matches `solve` and matches the CPU's resolution.
            dt_max = self._timestep_max(tend)
        dt_max = self._element_cap(dt_max)
            
        ## Stage 8(d)'s contract ("called at the start of every analysis")
        ## applied here too: a previous CPU run leaves per-element state --
        ## e.g. Idtmod's Phase-3 crossing cache -- that would otherwise leak
        ## into the STATIC breakpoint enumeration below.
        if hasattr(self.cir, 'reset_state'):
            self.cir.reset_state(getattr(self, 'epar', None))
        t_breaks_array = jnp.array(collect_breakpoints(
            self.cir, tend, float(self.par.minbreak)))
        
        # Setup batched initial state -- the SAME shared machinery as solve
        # (P12): every lane starts from the ic-resolved vector, since the
        # overrides sweep element parameters, not starting states.
        x0_batch = (_dc_x0 if _dc_x0 is not None else
                    jnp.tile(jnp.asarray(self._initial_state(refnode),
                                         dtype=jnp.float64), (batch_size, 1)))
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
        _state_mask = (self._jax_state_row_mask(self._initial_state(refnode))
                       if self.par.coupled_lte else None)
        _pcnr_meta, _pcnr_VT = self._pcnr_setup()
        ## Phase-2 gauge shift, PER LANE (roadmap item 3): swept moduli/
        ## offsets are re-declared per lane, so swept lanes keep the
        ## bounded-state property.  The lane arrays travel as vmapped
        ## arguments; the row set is static and shared.  Zero-width arrays
        ## when nothing declares -- the shift compiles to a no-op then.
        _pb = self._periodic_state_arrays_batched(override_params_tree,
                                                  batch_size)
        if _pb is None:
            _p_rows = jnp.zeros((0,), dtype=jnp.int32)
            _p_mods_lanes = jnp.ones((batch_size, 0))
            _p_offs_lanes = jnp.zeros((batch_size, 0))
        else:
            _p_rows, _p_mods_lanes, _p_offs_lanes = _pb

        def run_chunk(s, p_tree, p_mods_lane, p_offs_lane):
            return outer_time_loop(s, self.cir, tend, CHUNK_SIZE, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method=self._eval_method(), params_tree=p_tree, reltol=self.par.reltol, abstol=self._newton_abstol(n), periodic_states=(_p_rows, p_mods_lane, p_offs_lane),
                                   xtol=self._newton_xtol(n),
                                   maxiter=self.par.maxiter, trtol=self.par.TRTOL, analysis=self.par.analysis, s_gamma_min=_sgm, s_gamma_max=_sgx, s_eta=_seta, relref=self._relref(), n_nodes=len(self.cir.nodes), state_mask=_state_mask, dv_bounds=self._dv_step_bounds(),
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
            
            final_state = batched_run_chunk(state, override_params_tree,
                                            _p_mods_lanes, _p_offs_lanes)

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
        ## silently ignored.  `include_state=False`, as on the CPU: an
        ## Idt/Idtmod ic PINS the DC operating point (LRM semantics), so it
        ## is meaningful without uic -- and the DC below goes through
        ## dcanalysis.DC, which sets the pin flag.
        if (self.par.ic or self._descendant_has_ic(self.cir,
                                                   include_state=False)) \
                and not uic:
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
    
        ## Stage 8(d)'s contract ("called at the start of every analysis")
        ## applied here too: a previous CPU run leaves per-element state --
        ## e.g. Idtmod's Phase-3 crossing cache -- that would otherwise leak
        ## into the STATIC breakpoint enumeration below.
        if hasattr(self.cir, 'reset_state'):
            self.cir.reset_state(getattr(self, 'epar', None))
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
        _state_mask = self._jax_state_row_mask(x0) if self.par.coupled_lte else None
        _pcnr_meta, _pcnr_VT = self._pcnr_setup()
        _periodic = self._periodic_state_arrays()

        @jax.jit
        def run_chunk(s):
            return outer_time_loop(s, self.cir, tend, CHUNK_SIZE, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method=self._eval_method(), reltol=self.par.reltol, abstol=self._newton_abstol(n), periodic_states=_periodic,
                                   fixed_timestep=fixed_timestep, grid_dt=timestep,
                                   provided_function=provided_function,
                                   xtol=self._newton_xtol(n),
                                   maxiter=self.par.maxiter, trtol=self.par.TRTOL, analysis=self.par.analysis, s_gamma_min=_sgm, s_gamma_max=_sgx, s_eta=_seta, relref=self._relref(), n_nodes=len(self.cir.nodes), state_mask=_state_mask, dv_bounds=self._dv_step_bounds(),
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
