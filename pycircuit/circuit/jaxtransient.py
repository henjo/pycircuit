from typing import NamedTuple, Any
import numpy as np
import jax
import jax.numpy as jnp

from pycircuit.circuit.analysis import Analysis
from pycircuit.circuit.analysis import CircuitResult as Result

class NewtonState(NamedTuple):
    x: Any
    xdiff: Any
    F_norm: Any
    iters: Any

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

# ---------------------------------------------------------------------------
# Phase 1: Pure Functional Integrators
# ---------------------------------------------------------------------------

def backward_euler_step(q_curr, C_curr, q_prev, dt):
    i_curr = (q_curr - q_prev) / dt
    g_curr = C_curr / dt
    return i_curr, g_curr

def trap_step(q_curr, C_curr, q_prev, iq_prev, dt):
    i_curr = (2.0 / dt) * (q_curr - q_prev) - iq_prev
    g_curr = (2.0 / dt) * C_curr
    return i_curr, g_curr


def gear2_step(q_curr, C_curr, q_prev1, q_prev2, dt_curr, dt_prev):
    alpha0 = (2.0 * dt_curr + dt_prev) / (dt_curr * (dt_curr + dt_prev))
    alpha1 = -(dt_curr + dt_prev) / (dt_curr * dt_prev)
    alpha2 = dt_curr / (dt_prev * (dt_curr + dt_prev))
    
    i_curr = alpha0 * q_curr + alpha1 * q_prev1 + alpha2 * q_prev2
    g_curr = alpha0 * C_curr
    return i_curr, g_curr

def compute_integration(q_curr, C_curr, state: TransientState, method='trap'):

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
        fallback = jnp.logical_or(state.step_idx < 2, dt_prev == 0.0)
        
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

def newton_inner_loop(state: TransientState, circuit, irefnode, tline_params, tline_indices, eval_method='euler', reltol=1e-3, abstol=1e-6, maxiter=50, params_tree=None):
    x_init = extrapolate_predictor(state)
    def interp_tlines(t_target, params, history, head):
        # t_target is scalar, history is (10000, 5)
        def find_idx(idx, _):
            curr_t = history[(head - idx) % 10000, 0]
            return jax.lax.cond(curr_t <= t_target, lambda: idx, lambda: idx + 1)
    
        # Simple while loop to find the interval
        def cond_fun(val):
            idx = val
            curr_t = history[(head - idx) % 10000, 0]
            # stop if curr_t <= t_target or we hit max history (10000)
            return jnp.logical_and(curr_t > t_target, idx < 9999)
        
        def body_fun(val):
            return val + 1
        
        idx1 = jax.lax.while_loop(cond_fun, body_fun, 0)
        idx0 = jnp.maximum(0, idx1 - 1)
    
        idx1_mapped = (head - idx1) % 10000
        idx0_mapped = (head - idx0) % 10000
    
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

    
    ## Convergence uses the passed tolerances (previously hardcoded to
    ## F_norm > 1e-6 and iters < 20, ignoring the reltol/abstol/maxiter args).
    ## conv_tol is set below once the initial residual F_norm0 is known.
    def make_cond(conv_tol):
        def cond_fun(nr_state: NewtonState):
            return jnp.logical_and(nr_state.F_norm > conv_tol, nr_state.iters < maxiter)
        return cond_fun
    
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
    
        F_norm = jnp.sum(jnp.abs(F))
    
        J_sub = jnp.delete(jnp.delete(J, irefnode, axis=0), irefnode, axis=1)
        F_sub = jnp.delete(F, irefnode)
        xdiff_sub = jnp.linalg.solve(J_sub, -F_sub)
        xdiff = jnp.insert(xdiff_sub, irefnode, 0.0)
    
        max_dv = 0.5
        max_diff = jnp.max(jnp.abs(xdiff))
        alpha = jnp.where(max_diff > max_dv, max_dv / max_diff, 1.0)
        x_next = x + alpha * xdiff
    
        return NewtonState(x=x_next, xdiff=xdiff, F_norm=F_norm, iters=nr_state.iters + 1)
    
    I_G0 = circuit.i(x_init, params_tree=params_tree)
    G_G0 = circuit.G(x_init, params_tree=params_tree)
    q_C0 = circuit.q(x_init, params_tree=params_tree)
    C_C0 = circuit.C(x_init, params_tree=params_tree)
    I_u0 = circuit.u(state.t + state.dt, analysis='tran', params_tree=params_tree)
    I_u0 = apply_tlines(I_u0, state.t + state.dt)
    i_C0, G_C_eq0 = compute_integration(q_C0, C_C0, state, method=eval_method)

    F0 = I_G0 + i_C0 + I_u0
    J0 = G_G0 + G_C_eq0
    F_norm0 = jnp.sum(jnp.abs(F0))

    J_sub0 = jnp.delete(jnp.delete(J0, irefnode, axis=0), irefnode, axis=1)
    F_sub0 = jnp.delete(F0, irefnode)
    xdiff_sub0 = jnp.linalg.solve(J_sub0, -F_sub0)
    xdiff0 = jnp.insert(xdiff_sub0, irefnode, 0.0)

    max_dv = 0.5
    max_diff = jnp.max(jnp.abs(xdiff0))
    alpha = jnp.where(max_diff > max_dv, max_dv / max_diff, 1.0)
    x_next0 = x_init + alpha * xdiff0

    initial_nr_state = NewtonState(x=x_next0, xdiff=xdiff0, F_norm=F_norm0, iters=1)

    ## Relative + absolute residual test on the L1 residual norm.
    conv_tol = abstol + reltol * F_norm0
    return jax.lax.while_loop(make_cond(conv_tol), body_fun, initial_nr_state)

# ---------------------------------------------------------------------------
# Phase 3: Outer Time Loop & Adaptive Control
# ---------------------------------------------------------------------------

def estimate_lte(q_curr, state: TransientState, method='trap'):
    """Charge-domain LTE estimate and the method order+1.

    Returns ``(lte_vector, order_plus_one)``.  ``order_plus_one`` is the step
    predictor exponent denominator (2 for 1st-order Euler, 3 for the 2nd-order
    methods).  Unlike the CPU transient this estimator is charge-domain (it does
    not apply ``J^-1``), matching the existing JAX design.
    """
    dt = state.dt
    q_prev = state.q_history[0]
    dt_prev = jnp.where(state.h_history[0] == 0.0, dt, state.h_history[0])

    q_prev2 = state.q_history[1]

    if method == 'euler':
        # Backward Euler: LTE ~ 0.5 h^2 q''.  A second divided difference of the
        # charge IS q'', so this branch was always right, and it sets the
        # convention the two below must follow: return the absolute LTE in charge
        # units, with the method's own error constant.
        dd1_n = (q_curr - q_prev) / dt
        dd1_nm1 = (q_prev - q_prev2) / dt_prev
        q_double_prime = (dd1_n - dd1_nm1) / dt
        return 0.5 * (dt**2) * jnp.abs(q_double_prime), 2.0

    # --- BOTH 2ND-ORDER METHODS NEED q''', NOT q'' ---
    #
    # Until 2026-07 both branches below estimated a second divided difference of
    # the charge -- which is q'' -- and scaled it by h^3.  That is the same defect
    # the CPU `Gear2Integrator`'s 'classic' branch had (fixed in `integrator.py`,
    # see its comment): dimensionally it is not a charge-error at all, and it
    # undershoots the truncation error by a factor of order h*omega.  The old
    # comment here even asserted "h^3 q''' via a second divided difference of the
    # charge", which is self-contradictory -- that difference is q''.
    #
    # A third divided difference of q is unavailable (only two past charges are
    # kept), so q''' is taken as twice the SECOND divided difference of the
    # companion current g = dq/dt, read off `iq_history` -- exactly what the CPU
    # implementation and the YWR estimator both do.
    if method in ('gear', 'gear2'):
        # VSS BDF-2 companion current at t_n; the alpha coefficients annihilate
        # the q' and q'' terms by construction.
        h1, h2 = dt, dt_prev
        alpha0 = (2 * h1 + h2) / (h1 * (h1 + h2))
        alpha1 = -(h1 + h2) / (h1 * h2)
        alpha2 = h1 / (h2 * (h1 + h2))
        g_n = alpha0 * q_curr + alpha1 * q_prev + alpha2 * q_prev2
        error_constant = 2.0 / 9.0          # BDF-2: LTE = (2/9) h^3 q'''
    else:
        # Trapezoidal companion current at t_n.
        g_n = 2.0 * (q_curr - q_prev) / dt - state.iq_history[0]
        error_constant = 1.0 / 12.0         # TRAP: LTE = (1/12) h^3 q'''

    g_nm1 = state.iq_history[0]
    g_nm2 = state.iq_history[1]
    # Second divided difference of g over t_n, t_{n-1}, t_{n-2} is q'''/2.
    dd2_g = ((g_n - g_nm1) / dt - (g_nm1 - g_nm2) / dt_prev) / (dt + dt_prev)
    q_triple_prime = 2.0 * dd2_g
    return error_constant * (dt**3) * jnp.abs(q_triple_prime), 3.0

def ywr_error_ratio(i_curr, x_curr, x_last, J, state: TransientState, irefnode,
                    method='trap', trtol=7.0, lte_rel=1e-3, lte_abs=1e-6):
    """Yao-Wang-Roychowdhury DAE LTE, returned as a normalized error ratio.

    Mirrors the CPU transient's ``lte_formula='ywr'``: forms the residual as a
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
        Eg = -0.5 * (g_n - g_nm1)
        order_p1 = 2.0
    elif method in ('gear', 'gear2'):
        h1, h2 = dt, dt_prev
        Eg = -(1.0 / 8.0) * ((h1 + h2) / (h1 * h2)) * \
            (h2 * g_n - (h1 + h2) * g_nm1 + h1 * g_nm2)
        order_p1 = 3.0
    else:  # trapezoidal
        Eg = -(1.0 / 6.0) * (g_n - 2.0 * g_nm1 + g_nm2)
        order_p1 = 3.0

    # lte = J^-1 Eg in the reduced (reference-node-removed) space.
    J_r = jnp.delete(jnp.delete(J, irefnode, axis=0), irefnode, axis=1)
    Eg_r = jnp.delete(Eg, irefnode)
    lte_r = jnp.linalg.solve(J_r, Eg_r)
    lte = jnp.insert(lte_r, irefnode, 0.0)

    etol = trtol * (lte_rel * jnp.maximum(jnp.abs(x_curr), jnp.abs(x_last)) + lte_abs)
    return jnp.max(jnp.abs(lte) / etol), order_p1


def lte_error_ratio(lte, q_curr, trtol=7.0, lte_rel=1e-3, lte_abs=1e-6):
    """Normalized LTE: <=1 means the step is within tolerance.

    ``trtol`` is folded into the tolerance so the accept test and the step-size
    predictor share the same target (the CPU transient had these disagree,
    causing an accept/reject oscillation)."""
    tol = trtol * (lte_rel * jnp.maximum(jnp.abs(q_curr), 1e-12) + lte_abs)
    return jnp.max(lte / tol)


def calculate_next_dt(dt, error_ratio, dt_min, dt_max, t_breaks_array, current_t, order_p1=2.0):
    # Predictor exponent 1/(order+1): sqrt for 1st order, cube root for 2nd.
    factor = jnp.where(error_ratio <= 0.0, 2.0, error_ratio ** (-1.0 / order_p1))
    factor = jnp.clip(factor, 0.1, 2.0)
    dt_new = dt * factor
    dt_new = jnp.clip(dt_new, dt_min, dt_max)

    future_breaks = jnp.where(t_breaks_array > current_t + 1e-12, t_breaks_array, jnp.inf)
    next_break = jnp.min(future_breaks)

    dt_new = jnp.where(current_t + dt_new > next_break, next_break - current_t, dt_new)
    dt_new = jnp.maximum(dt_new, dt_min)
    return dt_new

def outer_time_loop(initial_state: TransientState, circuit, tend, chunk_size, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method='euler', params_tree=None, reltol=1e-3, abstol=1e-6, maxiter=50, lte_formula='ywr'):

    def time_cond(state: TransientState):
        under_time = state.t < tend
        under_chunk = state.step_idx < chunk_size
        return jnp.logical_and(under_time, under_chunk)

    def time_body(state: TransientState):
        nr_state = newton_inner_loop(state, circuit, irefnode, tline_params, tline_indices,
                                     eval_method=eval_method, reltol=reltol, abstol=abstol,
                                     maxiter=maxiter, params_tree=params_tree)
        x_curr = nr_state.x
        q_curr = circuit.q(x_curr, params_tree=params_tree)
        C_curr = circuit.C(x_curr, params_tree=params_tree)
        i_curr, Geq = compute_integration(q_curr, C_curr, state, method=eval_method)

        ## lte_formula is a Python-static string, so this branch is resolved at
        ## trace time (no jax.lax.cond needed).
        if lte_formula == 'ywr':
            # DAE LTE: g-difference mapped through J^-1 (one extra G eval + solve).
            J = circuit.G(x_curr, params_tree=params_tree) + Geq
            error_ratio, order_p1 = ywr_error_ratio(
                i_curr, x_curr, state.x_history[0], J, state, irefnode,
                method=eval_method)
        else:
            # Charge-domain estimate (no J^-1).
            lte, order_p1 = estimate_lte(q_curr, state, method=eval_method)
            error_ratio = lte_error_ratio(lte, q_curr)

        ## Accept when the LTE is within tolerance.  Always accept the first step
        ## (no history for the estimate) and when dt has reached the floor -- the
        ## latter guarantees forward progress so a rejection loop cannot deadlock.
        at_floor = state.dt <= dt_min * (1.0 + 1e-9)
        first = state.step_idx < 1
        accept = jnp.logical_or(jnp.logical_or(error_ratio <= 1.0, at_floor), first)

        def do_accept(_):
            next_dt = calculate_next_dt(state.dt, error_ratio, dt_min, dt_max,
                                        t_breaks_array, state.t, order_p1)

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
                new_head = (state.tline_head + 1) % 10000
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
                tline_history=tline_history_new, tline_head=tline_head_new)

        def do_reject(_):
            ## LTE above tolerance: shrink the step (bounded below by dt_min) and
            ## retry the same time point without committing anything.
            shrink = jnp.clip(error_ratio ** (-1.0 / order_p1), 0.1, 0.9)
            retry_dt = jnp.maximum(state.dt * shrink, dt_min)
            return state._replace(dt=retry_dt)

        return jax.lax.cond(accept, do_accept, do_reject, None)

    return jax.lax.while_loop(time_cond, time_body, initial_state)

# ---------------------------------------------------------------------------
# Phase 4: Class Wrapper & Chunking Engine
# ---------------------------------------------------------------------------

class JAXTransient(Analysis):
    """
    Fully compiled GPU-native transient analysis.
    Implements adaptive LTE and Predictor-Corrector implicitly in the JAX kernel.
    """

    def __init__(self, cir, **kwargs):
        super().__init__(cir, **kwargs)
        self.toolkit = getattr(self.cir, 'toolkit', None)
        if not (self.toolkit and self.toolkit.supports('autodiff')):
            raise ValueError("JAXTransient requires the circuit to use _jaxtoolkit.py.")
        

    def solve_batched(self, irefnode, override_params_tree, tend, timestep=1e-12, CHUNK_SIZE=5000, dt_min=1e-15, dt_max=None, uic=False, lte_formula='ywr'):
        import jax
        
        self.cir.update_iparv()
        import jax.numpy as jnp
        import numpy as np
        
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
        x0 = jnp.zeros(n)
        if not uic:
            # We would normally solve DC here. For now, use 0 or user-provided x0.
            pass
            
        if dt_max is None:
            dt_max = tend / 10.0
            
        breakpoints = []
        for inst_name, elem in self.cir.elements.items():
            if hasattr(elem, 'next_event'):
                t_event = 0.0
                while t_event < tend:
                    t_event = elem.next_event(t_event)
                    if t_event is not None and t_event <= tend:
                        breakpoints.append(t_event)
                    else:
                        break
        breakpoints.append(tend)
        breakpoints = sorted(list(set(breakpoints)))
        t_breaks_array = jnp.array(breakpoints)
        
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
            tline_history = jnp.zeros((batch_size, n_tlines, 10000, 5))
            for i, (name, tline) in enumerate(tlines):
                nodemap = self.cir.elementnodemap[name]
                tline_indices = tline_indices.at[i].set(nodemap)
                tline_params = tline_params.at[i, 0].set(tline.iparv.TD)
                tline_params = tline_params.at[i, 1].set(tline.iparv.Z0)
                # Init with zero for now
        else:
            tline_history = jnp.zeros((batch_size, 0, 10000, 5))
            
        tline_head = jnp.zeros(batch_size, dtype=jnp.int32)
        
        def run_chunk(s, p_tree):
            return outer_time_loop(s, self.cir, tend, CHUNK_SIZE, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method='gear', params_tree=p_tree, lte_formula=lte_formula)
            
        # JIT the vmapped run_chunk
        batched_run_chunk = jax.jit(jax.vmap(run_chunk))
        
        results_list = [np.expand_dims(np.array(x0_batch), axis=1)]
        times_list = [np.zeros((batch_size, 1))]
        
        current_t = jnp.zeros(batch_size)
        current_dt = jnp.full(batch_size, timestep)
        
        while np.any(current_t < tend):
            res_buf = jnp.zeros((batch_size, CHUNK_SIZE, n))
            time_buf = jnp.zeros((batch_size, CHUNK_SIZE))
            
            state = TransientState(
                t=current_t, dt=current_dt, step_idx=jnp.zeros(batch_size, dtype=jnp.int32),
                x_history=x_hist, q_history=q_hist, iq_history=iq_hist, h_history=h_hist,
                results_buffer=res_buf, time_buffer=time_buf,
                tline_history=tline_history, tline_head=tline_head
            )
            
            final_state = batched_run_chunk(state, override_params_tree)
            
            # Since batching can lead to different step sizes, we just append the raw buffers for now
            valid_steps = np.array(final_state.step_idx)
            max_steps = int(np.max(valid_steps))
            
            if max_steps == 0:
                break
                
            x_chunk = np.array(final_state.results_buffer)[:, :max_steps, :]
            t_chunk = np.array(final_state.time_buffer)[:, :max_steps]
            
            # For sequences that ended early, fill their remaining chunk with their last valid value
            for b in range(batch_size):
                b_len = valid_steps[b]
                if b_len < max_steps:
                    x_chunk[b, b_len:max_steps, :] = x_chunk[b, b_len-1:b_len, :]
                    t_chunk[b, b_len:max_steps] = t_chunk[b, b_len-1]
                    
            results_list.append(x_chunk)
            times_list.append(t_chunk)
            
            current_t = final_state.t
            current_dt = final_state.dt
            x_hist = final_state.x_history
            q_hist = final_state.q_history
            iq_hist = final_state.iq_history
            h_hist = final_state.h_history
            tline_history = final_state.tline_history
            tline_head = final_state.tline_head
            
        # Concat along time axis (axis=1)
        all_results = np.concatenate(results_list, axis=1)
        all_times = np.concatenate(times_list, axis=1)
        
        res = []
        for i in range(batch_size):
            r = Result(self.cir, x=all_results[i].T, xdot=None,
                            sweep_values=all_times[i],
                            sweep_label='time',
                            sweep_unit='s')
            res.append(r)
        return res

    def solve(self, refnode=0, tend=1e-3, x0=None, timestep=1e-6, CHUNK_SIZE=5000, lte_formula='ywr', **kwargs):
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
        dt_max = timestep # Legacy transient limits dt_max to tstep by default
    
        breakpoints = []
        for elem in self.cir.elements:
            if hasattr(elem, 'next_event'):
                t_event = 0.0
                while t_event < tend:
                    t_event = elem.next_event(t_event)
                    if t_event is not None and t_event <= tend:
                        breakpoints.append(t_event)
                    else:
                        break
        breakpoints.append(tend)
        breakpoints = sorted(list(set(breakpoints)))
        t_breaks_array = jnp.array(breakpoints)
    
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
            tline_history = jnp.zeros((n_tlines, 10000, 5))
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
            tline_history = jnp.zeros((0, 10000, 5))
        
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
            return outer_time_loop(s, self.cir, tend, CHUNK_SIZE, irefnode, dt_min, dt_max, t_breaks_array, tline_params, tline_indices, eval_method='gear', lte_formula=lte_formula)
        
        results_list = [np.array([x0])]
        times_list = [np.array([0.0])]
    
        current_t = 0.0
        current_dt = timestep
    
        while current_t < tend:
            res_buf = jnp.zeros((CHUNK_SIZE, n))
            time_buf = jnp.zeros(CHUNK_SIZE)
        
            state = TransientState(
                t=jnp.array(current_t), dt=jnp.array(current_dt), step_idx=jnp.array(0),
                x_history=x_hist, q_history=q_hist, iq_history=iq_hist, h_history=h_hist,
                results_buffer=res_buf, time_buffer=time_buf,
                tline_history=tline_history, tline_head=tline_head
            )
        
            final_state = run_chunk(state)
        
            # Extract valid steps
            valid_idx = int(final_state.step_idx)
        
            # If no steps were taken (e.g. tend was already reached or error), break
            if valid_idx == 0:
                break
            
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
        
        # Concatenate
        all_results = np.vstack(results_list)
        all_times = np.concatenate(times_list)
    
        # Build Result
        from pycircuit.circuit.analysis import CircuitResult
        res = CircuitResult(self.cir, x=all_results.T, xdot=None,
                            sweep_values=all_times,
                            sweep_label='time',
                            sweep_unit='s')
        return res
