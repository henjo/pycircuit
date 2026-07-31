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


def _rc_circuit():
    from pycircuit.circuit.elements import SubCircuit, R, C, VS
    from pycircuit.circuit import gnd
    cir = SubCircuit()
    cir.add_node('in'); cir.add_node('out')
    cir['V1'] = VS('in', gnd, v=1.0)
    cir['R1'] = R('in', 'out', r=1e3)
    cir['C1'] = C('out', gnd, c=1e-6)
    return cir


def _with_jax_toolkit(fn):
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        return fn()
    finally:
        circuit_mod.default_toolkit = saved


def test_jax_tolerances_are_settable():
    """Stage 9(b): `JAXTransient(cir, reltol=...)` raised KeyError.

    There was no supported way to ask for a tighter run, and the kernel's
    hard-coded reltol=1e-3/abstol=1e-6 disagreed with Transient's 1e-4/1e-12, so
    the two backends solved to different accuracies while presenting as one
    analysis.  The names and defaults here are Transient's deliberately.
    """
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.transient import Transient

    tran = _with_jax_toolkit(lambda: JAXTransient(_rc_circuit(), reltol=1e-6))
    assert float(tran.par.reltol) == 1e-6

    dflt = _with_jax_toolkit(lambda: JAXTransient(_rc_circuit()))
    cpu = {p.name: p.default for p in Transient.parameters}
    for name in ('reltol', 'iabstol', 'vabstol', 'lte_vabstol', 'lte_iabstol'):
        assert float(getattr(dflt.par, name)) == float(cpu[name]), \
            '%s defaults differ between backends: JAX %r vs CPU %r' \
            % (name, getattr(dflt.par, name), cpu[name])


@pytest.mark.parametrize('unsupported', ['nrsolver', 'scaler'])
def test_jax_rejects_strategies_it_cannot_honour(unsupported):
    """Stage 9(b): both were accepted and silently ignored.

    The Newton loop is a traced while_loop with a fixed algorithm and no scaling
    step, so honouring them is not a matter of wiring.  Silent acceptance is the
    same shape of defect as the `lte_formula` knob removed in 9(f).
    """
    from pycircuit.circuit.jaxtransient import JAXTransient
    with pytest.raises(NotImplementedError, match=unsupported):
        _with_jax_toolkit(
            lambda: JAXTransient(_rc_circuit(), **{unsupported: object()}))


def test_jax_lte_abstol_is_per_row_not_scalar():
    """Stage 9(c): volts on node rows, amps on branch rows.

    A single scalar applies one physical kind of tolerance to rows of another --
    0.3a's residual-vs-solution defect, which 9(c) exists to not repeat.  Asserted
    by setting the two to different values so a scalar cannot pass.
    """
    from pycircuit.circuit.jaxtransient import JAXTransient
    tran = _with_jax_toolkit(
        lambda: JAXTransient(_rc_circuit(), lte_vabstol=1e-9, lte_iabstol=1e-15))
    v = np.asarray(tran._lte_abstol(tran.cir.n))
    n_nodes = len(tran.cir.nodes)

    assert v.shape == (tran.cir.n,)
    assert np.all(v[:n_nodes] == 1e-9), 'node rows must carry lte_vabstol'
    assert np.all(v[n_nodes:] == 1e-15), 'branch rows must carry lte_iabstol'


def test_jax_error_responds_to_reltol():
    """Stage 9, gate 9-1(c): the CPU gate that 'is not currently expressible'.

    It is expressible now that tolerances are settable.  Measured with a timestep
    fine enough that `dt_max = timestep` is not the binding constraint -- with the
    clamp binding, the whole-waveform error is pinned by the FIRST step and does
    not move with reltol at all, which is a property of the clamp rather than of
    the controller.  See 9(b)/(c) in the plan.
    """
    import warnings
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit import gnd
    tau, tend = 1e-3, 5e-3

    def run(reltol):
        def go():
            cir = _rc_circuit()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = JAXTransient(cir, reltol=reltol).solve(
                    gnd, tend=tend, timestep=1e-5, uic=True)
            return res, cir.get_node_index('out')
        res, idx = _with_jax_toolkit(go)
        t = np.asarray(res.sweep_values, dtype=float).reshape(-1)
        v = np.asarray(res.x[idx], dtype=float).reshape(-1)
        return float(np.abs(v - (1.0 - np.exp(-t / tau))).max())

    loose, tight = run(1e-3), run(1e-5)
    assert tight < loose, \
        'tightening reltol did not reduce the error: %.4e -> %.4e' % (loose, tight)
    assert loose / tight > 1.5, \
        'reltol barely moves the error: %.4e -> %.4e (%.2fx)' \
        % (loose, tight, loose / tight)


## ---------------------------------------------------------------------------
## Stage 9, gates 9-1 .. 9-3: the CPU's step-control gates, ported.
##
## The plan records these as "None of these is currently expressible, which is
## the asymmetry that let the copied LTE survive."  They are expressible now that
## tolerances are settable (9(b)/(c)), rejections are counted, and breakpoints
## work (9(d)).  Porting them immediately found the 3/4 Gear2 optimism for the
## THIRD time -- on the backend that never received stage 4i's fix.
## ---------------------------------------------------------------------------

def _lte_ratio(h, method, order, G=1.0e6):
    """estimate/h^order for a synthetic companion-current history.

    `g` is given degree `order` so the derivative THIS method differences is
    constant and non-zero.  One quadratic for all three makes g'(0) = 0, which
    kills Euler's leading term and makes it read as second order -- the same
    degeneracy the CPU units test hit, walked into again here before it was seen.
    """
    from pycircuit.circuit.jaxtransient import TransientState, ywr_error_ratio
    n = 2
    g = (lambda t: G * t) if order == 1 else (lambda t: G * t ** 2 / 2.0)
    st = TransientState(
        t=0.0, dt=h, step_idx=5, x_history=None,
        q_history=jnp.zeros((3, n)),
        iq_history=jnp.array([[g(-h), 0.0], [g(-2 * h), 0.0], [0.0, 0.0]]),
        h_history=jnp.array([h, h, 0.0]),
        results_buffer=None, time_buffer=None,
        tline_history=None, tline_head=None,
        sig_max=jnp.array(1.0), n_rejected=0)
    r, _p = ywr_error_ratio(
        jnp.array([g(0.0), 0.0]), jnp.array([1.0, 0.0]), jnp.array([1.0, 0.0]),
        jnp.eye(n), st, irefnode=1, method=method,
        trtol=1.0, lte_rel=0.0, lte_abstol=1.0)
    return float(r) / h ** order


@pytest.mark.parametrize('method,order', [('euler', 1), ('trap', 2), ('gear', 2)])
def test_gate_9_1a_jax_lte_scales_with_the_right_power_of_h(method, order):
    """Gate 9-1(a): the estimate must follow h^order, not merely be non-zero."""
    hs = (1e-4, 5e-5, 2.5e-5, 1.25e-5)
    scaled = [_lte_ratio(h, method, order) for h in hs]
    for got in scaled[1:]:
        assert got == pytest.approx(scaled[0], rel=1e-9, abs=0.0), \
            '%s: estimate/h^%d is not constant in h: %r' % (method, order, scaled)


def test_gate_9_1a_jax_gear2_error_constant_matches_the_cpu():
    """The 3/4 optimism, pinned so it cannot come back a fourth time.

    YWR's Table I GEAR2 residual gives (1/4) h^2 q''' against a true (1/3), so the
    JAX solver reported 3/4 of its own truncation error at every step -- on the
    only eval_method either entry point uses.  The CPU found and fixed this in
    stage 4i; the fix never crossed to jaxtransient.py.  Asserted against the
    DERIVED constant rather than a recorded number.
    """
    G = 1.0e6
    assert _lte_ratio(1e-4, 'gear', 2, G) == pytest.approx(G / 3.0, rel=1e-9), \
        'gear2 error constant is not q\'\'\'/3 -- the YWR 3/4 optimism is back'
    assert _lte_ratio(1e-4, 'trap', 2, G) == pytest.approx(G / 6.0, rel=1e-9)
    assert _lte_ratio(1e-4, 'euler', 1, G) == pytest.approx(G / 2.0, rel=1e-9)


def test_gate_9_3_jax_transient_lands_on_the_pulse_edges():
    """Gate 9-3: a VPulse under JAX hits the pulse edges.

    Blocked until 9(d): the scan found zero breakpoints, so the solver stepped
    straight through every discontinuity.  Also asserts the run lands exactly on
    `tend`, which it used to overshoot by up to one timestep.
    """
    import warnings
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import SubCircuit, R, C, VPulse
    from pycircuit.circuit import gnd

    td, tr_, pw, tf, per = 1e-4, 1e-6, 5e-4, 1e-6, 1e-3
    tend = 3e-3

    def go():
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VPulse('in', gnd, v1=0, v2=1, td=td, tr=tr_, tf=tf,
                           pw=pw, per=per)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return JAXTransient(cir).solve(gnd, tend=tend, timestep=1e-4, uic=True)

    res = _with_jax_toolkit(go)
    t = np.asarray(res.sweep_values, dtype=float).reshape(-1)

    edges = sorted({k * per + e
                    for k in range(4)
                    for e in (td, td + tr_, td + tr_ + pw, td + tr_ + pw + tf)
                    if 0 < k * per + e <= tend})
    for edge in edges:
        assert np.abs(t - edge).min() < 1e-12, \
            'no time point lands on the pulse edge at %g (closest %g)' \
            % (edge, t[int(np.abs(t - edge).argmin())])

    assert t[-1] == pytest.approx(tend, rel=0.0, abs=1e-12), \
        'run ended at %.9e, not the requested tend %.9e' % (t[-1], tend)


def test_gate_9_2_cpu_and_jax_agree_on_an_rc_transient():
    """Gate 9-2: a CPU/JAX agreement test IN THE SUITE.

    The stage-5 cross-check was run by hand and written into prose, so the next
    divergence between the two backends was invisible -- and there was one: the
    Gear2 error constant, which had disagreed by 4/3 since stage 4i.  Compared on
    a common time grid by interpolation, since the two choose their own steps.
    """
    import warnings
    from pycircuit.circuit import numeric, gnd
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import Gear2Integrator
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import SubCircuit, R, C, VS

    tend, tau = 5e-3, 1e-3

    cpu = SubCircuit(toolkit=numeric)
    cpu['V1'] = VS(1, gnd, v=1.0)
    cpu['R1'] = R(1, 2, r=1e3)
    cpu['C1'] = C(2, gnd, c=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rc_cpu = Transient(cpu, toolkit=numeric, integrator=Gear2Integrator(),
                           uic=True).solve(refnode=gnd, tend=tend, timestep=1e-5)
    t_cpu = np.asarray(rc_cpu.sweep_values, dtype=float)
    v_cpu = np.asarray(rc_cpu.v(2, gnd), dtype=float)

    def go():
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=1.0)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = JAXTransient(cir).solve(gnd, tend=tend, timestep=1e-5, uic=True)
        return res, cir.get_node_index('out')

    rj, idx = _with_jax_toolkit(go)
    t_jax = np.asarray(rj.sweep_values, dtype=float).reshape(-1)
    v_jax = np.asarray(rj.x[idx], dtype=float).reshape(-1)

    ## Both against the analytic solution, and against each other on a shared
    ## grid.  Agreement alone would be satisfied by two backends wrong the same
    ## way, which is exactly what a copied transcription produces.
    grid = np.linspace(5e-4, tend, 200)
    exact = 1.0 - np.exp(-grid / tau)
    e_cpu = np.abs(np.interp(grid, t_cpu, v_cpu) - exact).max()
    e_jax = np.abs(np.interp(grid, t_jax, v_jax) - exact).max()
    between = np.abs(np.interp(grid, t_cpu, v_cpu)
                     - np.interp(grid, t_jax, v_jax)).max()

    assert e_cpu < 5e-4, 'CPU is not tracking the analytic solution: %.3e' % e_cpu
    assert e_jax < 5e-4, 'JAX is not tracking the analytic solution: %.3e' % e_jax
    assert between < 5e-4, \
        'the two backends disagree by %.3e on a common grid' % between


def test_gate_9_1b_statistics_report_what_the_run_did():
    """Gate 9-1(b): 'a step is actually rejected' -- now countable, still 0.

    A rejected step advances neither `t` nor `step_idx`, so it leaves no trace in
    the output buffers; without this counter the gate could not be stated on this
    backend at all.  With it, the answer across 16 measured configurations --
    RC, pulse, stiff, oversized initial step, sine at two frequencies, tolerances
    to 1e-8 -- is that the reject branch is NEVER reached: the controller shrinks
    predictively on the accept path and `dt_max = timestep` caps growth, so the
    error ratio does not exceed 1 in practice.

    This asserts the bookkeeping, NOT that rejections are zero.  Pinning zero
    would freeze in a property of the current clamp, and if `dt_max` is ever
    changed (an open decision -- see 9(b)/(c) in the plan) rejections should start
    happening and that must not read as a regression.
    """
    import warnings
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit import gnd

    def go():
        cir = _rc_circuit()
        tran = JAXTransient(cir, reltol=1e-6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=5e-3, timestep=1e-4, uic=True)
        return tran, res

    tran, res = _with_jax_toolkit(go)
    st = tran.statistics
    n_returned = len(np.asarray(res.sweep_values).reshape(-1))

    assert st.accepted_steps > 0
    assert st.rejected_steps >= 0
    ## The buffer holds one row per accepted step, plus the seeded initial point.
    assert n_returned == st.accepted_steps + 1, \
        'returned %d points for %d accepted steps' % (n_returned, st.accepted_steps)
    ## sigglobal's running reference must have seen the 1 V drive.
    assert st.signal_max == pytest.approx(1.0, rel=1e-6), \
        'signal_max = %r; the sigglobal reference is not tracking the solution' \
        % st.signal_max


def test_gate_9_1b_sigglobal_reference_survives_a_chunk_boundary():
    """`sig_max` and `n_rejected` are running totals and must cross chunks.

    Rebuilding TransientState per chunk without them reset the sigglobal
    reference to zero every CHUNK_SIZE steps, so a long run silently reverted to
    a pointlocal-like tolerance at each boundary.  Same shape as the CPU's
    `_dt_last2` reset: a per-run quantity re-seeded by a per-call constructor.
    Forced here with a CHUNK_SIZE far smaller than the run.
    """
    import warnings
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit import gnd

    def go(chunk):
        cir = _rc_circuit()
        tran = JAXTransient(cir)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tran.solve(gnd, tend=5e-3, timestep=1e-4, uic=True, CHUNK_SIZE=chunk)
        return tran.statistics

    one_chunk = _with_jax_toolkit(lambda: go(5000))
    many_chunks = _with_jax_toolkit(lambda: go(7))

    assert many_chunks.signal_max == pytest.approx(one_chunk.signal_max, rel=1e-9), \
        'sigglobal reference differs across chunking: %r vs %r' \
        % (many_chunks.signal_max, one_chunk.signal_max)
    assert many_chunks.accepted_steps == one_chunk.accepted_steps, \
        'chunking changed the step count: %d vs %d' \
        % (many_chunks.accepted_steps, one_chunk.accepted_steps)
