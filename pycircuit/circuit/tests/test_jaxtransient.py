"""Tests for the JAX transient solver (pycircuit.circuit.jaxtransient)."""
import numpy as np
import pytest

jax = pytest.importorskip("jax")
import jax.numpy as jnp

from pycircuit.circuit.jaxtransient import (TransientState, compute_integration,
                                            ywr_error_ratio)


def _state(dt, q_hist, iq_hist, h_hist, step_idx):
    return TransientState(
        t=0.0, dt=dt, step_idx=step_idx, x_history=None,
        q_history=jnp.array(q_hist), iq_history=jnp.array(iq_hist),
        h_history=jnp.array(h_hist), results_buffer=None, time_buffer=None,
        tline_history=None, tline_head=None)


def test_jax_integration_dispatch_gear2():
    """method='gear' must run Gear2, not silently fall through to Euler.

    Regression: the dispatch used ``method in ('euler', 'gear')`` with an
    unreachable ``elif method == 'gear'``, so 'gear' ran Backward Euler and
    gear2_step was dead code.
    """
    q_curr = jnp.array([3.0])
    C_curr = jnp.array([[2.0]])
    # step_idx >= 2 and a non-zero previous step, so Gear2 does not fall back.
    st = _state(1.0, [[2.0], [0.0], [0.0]], [[0.5], [0.0], [0.0]], [1.0, 1.0, 0.0], 5)

    i_euler, _ = compute_integration(q_curr, C_curr, st, method='euler')
    i_gear, _ = compute_integration(q_curr, C_curr, st, method='gear')

    # Euler: (3-2)/1 = 1.0 ; Gear2: 1.5*3 - 2.0*2 + 0.5*0 = 0.5
    assert float(i_euler[0]) == pytest.approx(1.0)
    assert float(i_gear[0]) == pytest.approx(0.5)


def test_jax_lte_order_per_method():
    """The estimator returns the method's (order+1); Gear2 is 2nd order, not 1st.

    This used to test `estimate_lte`, the charge-domain path deleted in 9(f).  The
    order it returns drives `calculate_next_dt`'s exponent, so the coverage is kept
    and re-pointed at the surviving estimator rather than dropped with the old one.
    """
    n = 2
    x_curr = jnp.array([1.0, 0.5])
    x_last = jnp.array([0.9, 0.4])
    i_curr = jnp.array([1.0, 0.0])
    J = jnp.eye(n)
    st = TransientState(
        t=0.0, dt=1.0, step_idx=5, x_history=None,
        q_history=jnp.zeros((3, n)),
        iq_history=jnp.array([[0.5, 0.0], [0.25, 0.0], [0.0, 0.0]]),
        h_history=jnp.array([1.0, 1.0, 0.0]),
        results_buffer=None, time_buffer=None,
        tline_history=None, tline_head=None)

    orders = {}
    for method in ('euler', 'gear', 'trap'):
        _, p = ywr_error_ratio(i_curr, x_curr, x_last, J, st, irefnode=1,
                               method=method)
        orders[method] = float(p)

    assert orders['euler'] == 2.0
    assert orders['gear'] == 3.0
    assert orders['trap'] == 3.0


def test_jaxtransient_rc_charging():
    """End-to-end JAXTransient: RC charging from 0 matches the analytic curve.

    One estimator since 9(f): the charge-domain path was deleted, not defaulted
    away, so there is nothing left to parametrize over.
    """
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.elements import SubCircuit, R, C, VS
    from pycircuit.circuit import gnd
    from pycircuit.circuit.jaxtransient import JAXTransient

    ## The toolkit is the JAXToolkit *instance*, not the backend module.
    ## Elements reach differentiation through toolkit.jacobian(), which is a
    ## method on the class -- a bare module has no such thing, and passing one
    ## used to work only because every element carried its own `import jax`.
    saved_toolkit = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in')
        cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=1.0)
        cir['R1'] = R('in', 'out', r=1e3)      # tau = R*C = 1e-3
        cir['C1'] = C('out', gnd, c=1e-6)
        res = JAXTransient(cir).solve(gnd, tend=5e-3, timestep=1e-4, uic=True)
        out_idx = cir.get_node_index('out')
    finally:
        circuit_mod.default_toolkit = saved_toolkit

    t = np.array(res.sweep_values)
    vout = np.array(res.x[out_idx])
    tau = 1e-3
    v_analytic = 1.0 * (1.0 - np.exp(-t[-1] / tau))

    assert abs(t[-1] - 5e-3) < 1e-4               # reached tend
    assert len(t) > 10                            # actually took steps
    assert abs(vout[-1] - v_analytic) < 5e-3      # correct within tolerance
