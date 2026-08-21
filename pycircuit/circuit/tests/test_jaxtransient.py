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
        _, p = ywr_error_ratio(i_curr, x_curr, J, st, irefnode=1,
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

    It is expressible now that tolerances are settable.

    THE TIMESTEP HERE IS COARSE, AND THAT INVERTED WITH 9(g).  Before the opening
    step was ramped, the whole-waveform error was pinned by the FIRST step -- the
    largest and only unchecked one -- so the response only showed at a timestep
    fine enough to make that first step small.  With the ramp, the first step is
    `timestep*1e-3` and no longer dominates; what binds at a FINE timestep is now
    `dt_max = timestep` itself, which leaves the tolerance no room (measured:
    1.2148e-5 at reltol 1e-3 and 1.2141e-5 at 1e-5, a ratio of 1.00). At a coarse
    timestep the controller has room and `reltol` is what decides.
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
                    gnd, tend=tend, timestep=1e-4, uic=True)
            return res, cir.get_node_index('out')
        res, idx = _with_jax_toolkit(go)
        t = np.asarray(res.sweep_values, dtype=float).reshape(-1)
        v = np.asarray(res.x[idx], dtype=float).reshape(-1)
        return float(np.abs(v - (1.0 - np.exp(-t / tau))).max())

    ## THREE decades, not two.  Over 1e-3 -> 1e-5 the response on this circuit is
    ## 1.36x, which is real but too close to the bar to be a useful check; over
    ## 1e-3 -> 1e-6 it is 5.96x.  The span is part of the assertion, so it is stated
    ## rather than left to whichever pair happened to be picked.
    loose, tight = run(1e-3), run(1e-6)
    assert tight < loose, \
        'tightening reltol did not reduce the error: %.4e -> %.4e' % (loose, tight)
    assert loose / tight > 2.0, \
        'reltol barely moves the error: %.4e -> %.4e (%.2fx); measured 5.96x when ' \
        'this was written' % (loose, tight, loose / tight)


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
        jnp.array([g(0.0), 0.0]), jnp.array([1.0, 0.0]),
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


def test_gate_9e_nonconverged_newton_is_not_committed():
    """Stage 9(e): an unconverged iterate used to be committed, silently.

    `newton_inner_loop`'s while_loop exits on `F_norm <= conv_tol` OR
    `iters >= maxiter`, and the caller took `x` either way -- so the step's
    "solution" need not solve the circuit equations, and its LTE was computed from
    it.  Measured on a LINEAR RC, where the exact answer is known: maxiter=1 gave
    4.97e-2 max error against 4.28e-3, and reported nothing at all.

    A traced loop cannot raise, so the check is per chunk in Python -- deliberately
    not at the end of the run, because rejecting non-converged steps drives dt to
    dt_min and the end may never arrive (~5e12 steps on this circuit).
    """
    import warnings
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit import gnd

    def go(maxiter):
        cir = _rc_circuit()
        tran = JAXTransient(cir, maxiter=maxiter)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=5e-3, timestep=1e-4, uic=True)
        return tran, res, cir.get_node_index('out')

    with pytest.raises(NoConvergenceError, match='failed to converge'):
        _with_jax_toolkit(lambda: go(1))

    ## And a solve that DOES converge is untouched, including at a maxiter low
    ## enough that the old stale-by-one convergence flag would have failed it.
    tau = 1e-3
    errs = []
    for maxiter in (2, 3, 100):
        tran, res, idx = _with_jax_toolkit(lambda m=maxiter: go(m))
        t = np.asarray(res.sweep_values, dtype=float).reshape(-1)
        v = np.asarray(res.x[idx], dtype=float).reshape(-1)
        errs.append(float(np.abs(v - (1.0 - np.exp(-t / tau))).max()))
        assert tran.statistics.forced_nonconverged_steps == 0
    ## Close, not bit-identical.  Both solves converge by the flag's own test, but
    ## a lower cap can stop an iteration earlier within the tolerance band, and with
    ## 9(g)'s ramp the step sequences then diverge slightly.  A 1e-12 equality was
    ## asserting that `maxiter` cannot influence a converged answer AT ALL, which was
    ## true only while the opening step dominated everything.
    assert errs[0] == pytest.approx(errs[-1], rel=0.05), \
        'maxiter materially changed a converged answer: %r' % errs


def test_gate_9e_converged_flag_is_evaluated_at_the_returned_point():
    """The flag must describe the x that is returned, not the one before it.

    `NewtonState.F_norm` is the residual at the iterate BEFORE the update that
    produced `final.x`, so the loop always does one update beyond the one its test
    approved.  Reading `F_norm <= conv_tol` as "converged" reports failure whenever
    the iteration cap is hit even when the returned point is good -- caught here at
    maxiter=2, which had been giving exactly the maxiter=100 answer and would have
    started raising.
    """
    import warnings
    from pycircuit.circuit.jaxtransient import newton_inner_loop, TransientState
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit import circuit as circuit_mod

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = _rc_circuit()
        cir.update_iparv()
        n = cir.n
        x0 = jnp.zeros(n)
        q0 = cir.q(x0)
        st = TransientState(
            t=jnp.array(0.0), dt=jnp.array(1e-5), step_idx=jnp.array(1),
            x_history=jnp.array([x0, x0, x0]), q_history=jnp.array([q0, q0, q0]),
            iq_history=jnp.zeros((3, n)), h_history=jnp.array([1e-5, 1e-5, 0.0]),
            results_buffer=None, time_buffer=None,
            tline_history=jnp.zeros((0, 10000, 5)),
            tline_head=jnp.array(0, dtype=jnp.int32))
        irefnode = cir.get_node_index(gnd_node())
        out = newton_inner_loop(
            st, cir, irefnode, jnp.zeros((0, 4)),
            jnp.zeros((0, 6), dtype=jnp.int32), eval_method='gear',
            reltol=1e-4, abstol=1e-12, maxiter=40)
    finally:
        circuit_mod.default_toolkit = saved

    ## A linear RC converges in a couple of iterations, so a flag computed at the
    ## returned point must say so.
    assert bool(out.converged), \
        'a converged linear solve reported non-convergence (F_norm=%r)' % out.F_norm


def gnd_node():
    from pycircuit.circuit import gnd
    return gnd


def test_gate_9_1b_the_reject_path_actually_fires():
    """Gate 9-1(b): a step is actually rejected -- a NON-ZERO reading.

    The counter added to state this gate carried a defect that produced exactly
    the answer needing no explanation: `do_accept` rebuilt TransientState without
    naming `n_rejected`, so the NamedTuple default zeroed it on every accepted
    step and the end-of-run value was structurally 0.  "Zero rejections across 16
    configurations" was reported as a finding.  It was a measurement of the bug.

    This test therefore demands a NON-ZERO count on a case chosen to produce one.
    A test asserting `>= 0` -- which is what was written first -- passes against
    the defect and is worth nothing here.
    """
    import warnings
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit import gnd

    ## RE-DERIVED at F17 (doc/transient_review_260820.md): the safety
    ## factor made the plain RC at reltol=1e-8 take ZERO rejections -- the
    ## old case stopped producing the nonzero reading this gate exists for.
    ## A sinusoidal drive still rejects (the error genuinely oscillates
    ## against the tolerance), so the gate keeps its teeth there.
    def go():
        from pycircuit.circuit.elements import SubCircuit, R, C, VSin
        cir = SubCircuit()
        cir['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
        cir['R'] = R('a', 'b', r=1e3)
        cir['C'] = C('b', gnd, c=1e-7)
        tran = JAXTransient(cir, reltol=1e-6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tran.solve(gnd, tend=5e-3, timestep=1e-4, uic=True)
        return tran.statistics

    st = _with_jax_toolkit(go)
    assert st.rejected_steps > 0, \
        'the reject branch never fired on rc-vsin at reltol=1e-6 (%r) -- either ' \
        'step control is not working or the counter is being reset again' % st
    assert st.accepted_steps > 0


def test_gate_9g_opening_step_is_ramped():
    """Stage 9(g): a run opening at `timestep` -- which is also dt_max -- makes the
    first step both the largest and the only unchecked one, and its error then
    dominates.  Ported from Transient._opening_step (stage 3).

    Pinned by the consequence rather than the value: with the ramp, `reltol` must
    move the answer.  Before it, the error sat at 4.2535e-3 across FOUR DECADES of
    reltol, identical to five figures, always at the first step.
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
                    gnd, tend=tend, timestep=1e-4, uic=True)
            return res, cir.get_node_index('out')
        res, idx = _with_jax_toolkit(go)
        t = np.asarray(res.sweep_values, dtype=float).reshape(-1)
        v = np.asarray(res.x[idx], dtype=float).reshape(-1)
        return float(np.abs(v - (1.0 - np.exp(-t / tau))).max())

    loose, tight = run(1e-3), run(1e-6)
    assert tight < loose / 2.0, \
        'reltol still barely moves the answer: %.4e -> %.4e; the opening step is ' \
        'probably dominating again' % (loose, tight)

    ## The ramp itself, and its validation, follow Transient's.
    tran = _with_jax_toolkit(lambda: JAXTransient(_rc_circuit()))
    assert tran._opening_step(1e-4) == pytest.approx(1e-7, rel=1e-12)
    tran2 = _with_jax_toolkit(lambda: JAXTransient(_rc_circuit(), firststep=1e-9))
    assert tran2._opening_step(1e-4) == pytest.approx(1e-9, rel=1e-12)
    with pytest.raises(ValueError, match='must be positive'):
        _with_jax_toolkit(
            lambda: JAXTransient(_rc_circuit(), firststep=-1.0))._opening_step(1e-4)


def test_solve_batched_runs_and_honours_timestep():
    """`solve_batched` had NO test coverage at all, and it was broken.

    The counters added by 9(c)/9(e) default to scalars; every field of the batched
    state must carry a leading batch axis for `jax.vmap`, so `solve_batched` raised
    "rank should be at least 1" while the whole suite stayed green.  It also used
    `dt_max = tend/10` where `solve` uses `timestep`, so the requested timestep was
    very nearly ignored -- 29 vs 32 steps and the same 3.67e-3 error whether 1e-4
    or 1e-5 was asked for.
    """
    import warnings
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit import gnd
    tau, tend = 1e-3, 5e-3

    def run(timestep):
        def go():
            cir = _rc_circuit()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                ## 'C', not 'R': R has no eval_*_pure, so an R override is not
                ## batchable and now raises (doc/transient_review_260820.md,
                ## F2(a)).  The old {'R': ...} here was a silent no-op -- equal
                ## values in both lanes, ignored either way.  C is batchable and
                ## the value matches the circuit's own, so this test keeps
                ## measuring what it always measured: timestep honouring.
                res = JAXTransient(cir).solve_batched(
                    gnd, override_params_tree={'C': {'c': jnp.array([[1e-6], [1e-6]])}},
                    tend=tend, timestep=timestep, uic=True)
            return res
        res = _with_jax_toolkit(go)
        r0 = res[0] if isinstance(res, (list, tuple)) else res
        t = np.asarray(r0.sweep_values, dtype=float).reshape(-1)
        v = np.asarray(r0.v('out'), dtype=float).reshape(-1)
        return len(t), float(np.diff(t).max()), \
            float(np.abs(v - (1.0 - np.exp(-t / tau))).max())

    n_coarse, dt_coarse, e_coarse = run(1e-4)
    n_fine, dt_fine, e_fine = run(1e-5)

    assert dt_coarse <= 1e-4 * 1.0001, 'dt_max is not the requested timestep'
    assert dt_fine <= 1e-5 * 1.0001, 'dt_max is not the requested timestep'
    assert n_fine > 3 * n_coarse, \
        'a 10x finer timestep barely changed the run: %d -> %d steps' \
        % (n_coarse, n_fine)
    assert e_fine < e_coarse, \
        'a finer timestep did not improve the answer: %.3e -> %.3e' \
        % (e_coarse, e_fine)
