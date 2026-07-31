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


def _pulse_circuit():
    from pycircuit.circuit.elements import SubCircuit, R, C, VPulse
    from pycircuit.circuit import gnd
    cir = SubCircuit()
    cir.add_node('in'); cir.add_node('out')
    cir['V1'] = VPulse('in', gnd, v1=0, v2=1, td=1e-4, tr=1e-5, tf=1e-5,
                       pw=5e-4, per=1e-3)
    cir['R1'] = R('in', 'out', r=1e3)
    cir['C1'] = C('out', gnd, c=1e-6)
    return cir


def test_jax_breakpoints_are_found_and_exact():
    """Stage 9(d): the scan used to find ZERO breakpoints, always.

    `JAXTransient.solve` iterated `for elem in self.cir.elements` -- a dict, so it
    yielded string keys and `hasattr('V1', 'next_event')` was False.  Nothing in
    the suite noticed because a transient with no breakpoints still runs; it just
    steps straight through every source discontinuity.

    Asserted against the analytic edge times rather than a recorded count, so the
    test says what a breakpoint IS rather than how many there were on the day.
    """
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import collect_breakpoints

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        bps = collect_breakpoints(_pulse_circuit(), 3e-3)
    finally:
        circuit_mod.default_toolkit = saved

    td, tr, pw, tf, per = 1e-4, 1e-5, 5e-4, 1e-5, 1e-3
    expected = sorted({k * per + e
                       for k in range(4)
                       for e in (td, td + tr, td + tr + pw, td + tr + pw + tf, per)
                       if 0 < k * per + e <= 3e-3} | {3e-3})

    assert len(bps) == len(expected), \
        'got %d breakpoints, expected %d: %r' % (len(bps), len(expected), bps)
    for got, want in zip(bps, expected):
        assert got == pytest.approx(want, rel=1e-12, abs=0.0), \
            'breakpoint %r != analytic edge %r' % (got, want)


def test_jax_breakpoint_scan_cannot_hang():
    """Stage 9(d): a non-advancing next_event must not spin forever.

    `Pulse.next_event` returned `t` itself at t=0, so the enumeration never
    advanced and `solve_batched` HUNG -- a wall-clock timeout rather than a stack
    trace.  That is fixed at the source, so this test breaks the contract on
    purpose to prove the second line of defence still holds.
    """
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import collect_breakpoints

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = _pulse_circuit()
        cir.elements['V1'].next_event = lambda t: t      # never advances
        with pytest.warns(RuntimeWarning, match='does not advance'):
            bps = collect_breakpoints(cir, 3e-3)
    finally:
        circuit_mod.default_toolkit = saved

    ## It returns, and returns something usable: tend is always present.
    assert bps == [3e-3]
