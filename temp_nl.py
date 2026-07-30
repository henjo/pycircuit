def newton_inner_loop(state: TransientState, circuit, irefnode, tline_params, tline_indices, eval_method='euler', reltol=1e-3, abstol=1e-6, maxiter=50):
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

    
    def cond_fun(nr_state: NewtonState):
        return jnp.logical_and(nr_state.F_norm > 1e-6, nr_state.iters < 20)
        
    def body_fun(nr_state: NewtonState):
        x = nr_state.x
        
        I_G = circuit.i(x)
        G_G = circuit.G(x)
        q_C = circuit.q(x)
        C_C = circuit.C(x)
        I_u = circuit.u(state.t + state.dt, analysis='tran')
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
        
    I_G0 = circuit.i(x_init)
    G_G0 = circuit.G(x_init)
    q_C0 = circuit.q(x_init)
    C_C0 = circuit.C(x_init)
    I_u0 = circuit.u(state.t + state.dt, analysis='tran')
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
    
    return jax.lax.while_loop(cond_fun, body_fun, initial_nr_state)
