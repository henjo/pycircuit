"""The JAX transient's failed-time-point continuation chain (sec. 50).

The CPU has rescued a failed time point since P18/P25: below `minstep` it
swaps in `Transient._rescue_solver` and walks junction-gmin -> gshunt ->
pseudo-transient.  That is Python control flow -- an exception per failed
point and a runtime choice of solver object -- and none of it survives
being traced into a `lax.while_loop`.  This is the ported twin.

What the tests below actually protect, in order of how easy it is to
break: the chain's ORDER and its pure-only acceptance; the fact that an
inactive ladder runs ZERO rungs (the whole cost argument, and the one
property a `lax.cond` port would silently lose under `vmap`); and
non-interference -- a run that never reaches the dt floor must produce the
same numbers bit for bit whether the chain is armed or not.
"""
import warnings

import numpy as np
import pytest

jax = pytest.importorskip('jax')
import jax.numpy as jnp

from pycircuit.circuit import circuit as circuit_mod, gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS, VSin, Diode
from pycircuit.circuit.jaxtransient import (
    JAXTransient, NewtonState, TransientState, _adaptive_ladder_traced,
    _gmin_junction_rows, newton_inner_loop, transient_point_rescue)


def _with_jax(fn):
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        return fn()
    finally:
        circuit_mod.default_toolkit = saved


def _nr(x, converged):
    return NewtonState(x=x, xdiff=jnp.zeros_like(x), F_norm=jnp.asarray(0.0),
                       iters=0, converged=jnp.asarray(converged))


def _synthetic_rung(lands):
    """A rung solver that converges only for the named ladder.

    Every converged rung adds 1.0 to the seed, so the returned vector
    counts the rungs the driver actually ran -- which is how the tests
    below distinguish "skipped" from "ran and failed".
    """
    def rung_newton(xs, kind, g):
        conv = (kind == lands)          # `kind` is static, chosen at trace
        return _nr(xs + (1.0 if conv else 0.0), conv)
    return rung_newton


# ---------------------------------------------------------------------------
# The driver's `active` gate -- the cost argument
# ---------------------------------------------------------------------------

def test_inactive_ladder_runs_zero_rungs():
    """`active=False` must return the seed untouched, having run nothing.

    This is why `active` lands in the loop CONDITION rather than in a
    wrapping `lax.cond`: under `vmap` a batched predicate cannot branch,
    so `cond` degenerates to a `select` that evaluates both sides and the
    ladder runs in full on every lane.  Measured at 59x against no ladder
    with ZERO lanes failing; gating the condition measured 0.9x.
    """
    def rung(xs, g):
        return xs + 1.0, jnp.asarray(True)

    seed = jnp.zeros(3)
    x_off, c_off = _adaptive_ladder_traced(rung, seed, e_start=-2.0,
                                           e_end=-12.0,
                                           active=jnp.asarray(False))
    assert not bool(c_off)
    ## zero rungs ran: the counter never moved off the seed
    assert np.all(np.asarray(x_off) == 0.0)

    x_on, c_on = _adaptive_ladder_traced(rung, seed, e_start=-2.0,
                                         e_end=-12.0,
                                         active=jnp.asarray(True))
    assert bool(c_on)
    assert float(np.asarray(x_on)[0]) > 0.0


def test_inactive_ladder_is_free_under_vmap_across_lanes():
    """The gate must survive batching: a lane with active=False runs no
    rungs even while a sibling lane is marching."""
    def rung(xs, g):
        return xs + 1.0, jnp.asarray(True)

    def one(active):
        return _adaptive_ladder_traced(rung, jnp.zeros(2), e_start=-2.0,
                                       e_end=-12.0, active=active)[0]

    out = np.asarray(jax.vmap(one)(jnp.asarray([True, False, True])))
    assert out[1, 0] == 0.0, 'an inactive lane ran rungs'
    assert out[0, 0] > 0.0 and out[2, 0] > 0.0


# ---------------------------------------------------------------------------
# The chain: order, pure-only acceptance, gating
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('lands', ['gmin', 'gshunt', 'ptc'])
def test_chain_lands_on_each_rung_in_turn(lands):
    """Whichever ladder can solve the point is the one that reports it."""
    rows = (jnp.asarray([0], jnp.int32), jnp.asarray([1], jnp.int32))
    st, rescued = transient_point_rescue(
        _synthetic_rung(lands), jnp.zeros(2), _nr(jnp.full(2, 9.0), False),
        jnp.asarray(True), gmin_rows=rows, n_nodes=2)
    assert bool(rescued)
    assert bool(st.converged)
    ## the rescued point is the ladder's solution, not the failed iterate
    assert float(np.asarray(st.x)[0]) != 9.0


def test_chain_gives_up_when_no_ladder_lands():
    """Nothing converges: the failed iterate and its False flag stand, so
    the caller still force-accepts and still counts it."""
    rows = (jnp.asarray([0], jnp.int32), jnp.asarray([1], jnp.int32))
    st, rescued = transient_point_rescue(
        _synthetic_rung('nothing'), jnp.zeros(2),
        _nr(jnp.full(2, 9.0), False), jnp.asarray(True),
        gmin_rows=rows, n_nodes=2)
    assert not bool(rescued)
    assert not bool(st.converged)
    assert np.all(np.asarray(st.x) == 9.0)


def test_chain_does_nothing_when_not_asked():
    """`needs_rescue=False` -- the ordinary step, above the dt floor.

    Even though this rung solver WOULD land, no ladder may run: this is
    the property that keeps a healthy run at its pre-chain cost.
    """
    rows = (jnp.asarray([0], jnp.int32), jnp.asarray([1], jnp.int32))
    st, rescued = transient_point_rescue(
        _synthetic_rung('gmin'), jnp.zeros(2), _nr(jnp.full(2, 9.0), False),
        jnp.asarray(False), gmin_rows=rows, n_nodes=2)
    assert not bool(rescued)
    assert np.all(np.asarray(st.x) == 9.0)


def test_chain_never_overrides_a_converged_newton():
    """A point the plain Newton solved is not re-solved by the chain."""
    rows = (jnp.asarray([0], jnp.int32), jnp.asarray([1], jnp.int32))
    st, rescued = transient_point_rescue(
        _synthetic_rung('gmin'), jnp.zeros(2), _nr(jnp.full(2, 4.0), True),
        jnp.asarray(True), gmin_rows=rows, n_nodes=2)
    assert not bool(rescued)
    assert bool(st.converged)
    assert np.all(np.asarray(st.x) == 4.0)


def test_junctionless_circuit_skips_the_gmin_ladder_structurally():
    """No declared junction -> no scatter targets -> the physical homotopy
    is not in the graph at all, exactly as at DC."""
    st, rescued = transient_point_rescue(
        _synthetic_rung('gmin'), jnp.zeros(2), _nr(jnp.full(2, 9.0), False),
        jnp.asarray(True), gmin_rows=None, n_nodes=2)
    ## gmin was the only ladder that could land, and it was skipped
    assert not bool(rescued)


# ---------------------------------------------------------------------------
# The deformation terms
# ---------------------------------------------------------------------------

def _rc_diode():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e4)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    c['D'] = Diode('b', gnd)
    return c


def test_zero_conductance_rung_is_the_undeformed_system():
    """The P22 rule: an accepted rescued point carries NO residue.

    Every ladder ends at exactly g = 0, and that final rung must be the
    circuit's own equations -- not the circuit plus a small conductance.
    Bit-for-bit, against a Newton called with no deformation at all.
    """
    def go():
        cir = _rc_diode()
        n = cir.n
        irefnode = cir.get_node_index(gnd)
        st = TransientState(
            t=0.0, dt=1e-7, step_idx=0,
            x_history=jnp.zeros((3, n)), q_history=jnp.zeros((3, n)),
            iq_history=jnp.zeros((3, n)), h_history=jnp.zeros(3),
            results_buffer=None, time_buffer=None,
            tline_history=None, tline_head=None)
        tp = jnp.zeros((0, 2))
        ti = jnp.zeros((0, 6), dtype=jnp.int32)
        rows = _gmin_junction_rows(cir)
        common = dict(eval_method='euler', reltol=1e-5, abstol=1e-12,
                      xtol=1e-12, maxiter=100, x0=jnp.zeros(n))
        plain = newton_inner_loop(st, cir, irefnode, tp, ti, **common)
        deformed = newton_inner_loop(
            st, cir, irefnode, tp, ti, gmin_rows=rows, gmin=0.0,
            gshunt_nodes=len(cir.nodes), gshunt=0.0, ptc_g=0.0,
            ptc_anchor=jnp.zeros(n), **common)
        return (np.asarray(plain.x, float), np.asarray(deformed.x, float),
                bool(plain.converged), bool(deformed.converged))

    xp, xd, cp, cd = _with_jax(go)
    assert cp and cd
    assert np.array_equal(xp, xd), 'g = 0 did not reproduce the pure system'


def test_a_deformed_rung_actually_moves_the_answer():
    """The companion of the test above: a NON-zero g must change the
    solution, or the zero-g agreement above would be vacuous."""
    def go():
        cir = _rc_diode()
        n = cir.n
        irefnode = cir.get_node_index(gnd)
        st = TransientState(
            t=0.0, dt=1e-7, step_idx=0,
            x_history=jnp.zeros((3, n)), q_history=jnp.zeros((3, n)),
            iq_history=jnp.zeros((3, n)), h_history=jnp.zeros(3),
            results_buffer=None, time_buffer=None,
            tline_history=None, tline_head=None)
        tp = jnp.zeros((0, 2))
        ti = jnp.zeros((0, 6), dtype=jnp.int32)
        common = dict(eval_method='euler', reltol=1e-5, abstol=1e-12,
                      xtol=1e-12, maxiter=100, x0=jnp.zeros(n))
        plain = newton_inner_loop(st, cir, irefnode, tp, ti, **common)
        shunted = newton_inner_loop(st, cir, irefnode, tp, ti,
                                    gshunt_nodes=len(cir.nodes),
                                    gshunt=1e-2, **common)
        return np.asarray(plain.x, float), np.asarray(shunted.x, float)

    xp, xs = _with_jax(go)
    assert not np.allclose(xp, xs), 'gshunt = 1e-2 changed nothing'


# ---------------------------------------------------------------------------
# Non-interference at the analysis level
# ---------------------------------------------------------------------------

def test_continuation_is_off_by_default():
    """It compiles into every step, so it is opt-in -- see the parameter's
    recorded measurement.  A flipped default would tax every run."""
    par = _with_jax(lambda: JAXTransient(_rc_diode()).par)
    assert bool(par.continuation) is False


def test_armed_chain_does_not_change_a_healthy_run():
    """The run never reaches the dt floor, so no rung may execute -- and
    the waveform must be identical bit for bit, not merely close."""
    def run(cont):
        tran = JAXTransient(_rc_diode(), continuation=cont, pcnr=False)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=2e-4, timestep=1e-6, uic=True)
        return (np.asarray(res.v('b'), float).reshape(-1),
                tran.statistics.gmin_rescues,
                tran.statistics.forced_nonconverged_steps)

    v_off, r_off, f_off = _with_jax(lambda: run(False))
    v_on, r_on, f_on = _with_jax(lambda: run(True))
    assert r_off == 0 and r_on == 0
    assert f_off == 0 and f_on == 0
    assert np.array_equal(v_off, v_on), 'the armed chain perturbed the run'


def test_rescue_counter_is_reported():
    """`gmin_rescues` is named for the CPU's counter so the two backends'
    run reports can be compared field by field."""
    from pycircuit.circuit.jaxtransient import JAXTransientStatistics
    s = JAXTransientStatistics()
    assert s.gmin_rescues == 0
    assert 'rescued' in repr(s)


# ---------------------------------------------------------------------------
# The value demonstration
# ---------------------------------------------------------------------------

def _overshoot():
    """A coarse step on a fast signal into a junction.

    The plain Newton starts from the EXTRAPOLATED PREDICTOR, which at ten
    points per period overshoots the diode into its exponential; the
    chain's ladders start from the last ACCEPTED state, which is still
    good.  `minstep` is pinned at the step size on purpose: the step-size
    controller would otherwise absorb this by shrinking dt, which is
    exactly why a triggering circuit is so hard to fabricate (P18's scope
    finding, and the mechanism behind it).
    """
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=5.0, freq=2e6)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-12)
    c['D'] = Diode('b', gnd)
    return c


_TS = 2e-7


def _run_overshoot(cont, minstep=_TS):
    def go():
        tran = JAXTransient(_overshoot(), pcnr=False, continuation=cont,
                            reltol=1e-6, minstep=minstep, firststep=_TS)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=20 * _TS, timestep=_TS, uic=True)
        return (np.asarray(res.sweep_values, float),
                np.asarray(res.v('b'), float).reshape(-1), tran.statistics)
    return _with_jax(go)


def test_the_chain_completes_a_run_the_plain_newton_cannot():
    """THE point of the feature: without it the run dies at the dt floor.

    This is the reopening trigger `doc/backend_parity_260821.md` stated in
    the abstract -- a forced-non-converged exit the chain can complete --
    made concrete.
    """
    from pycircuit.circuit.nrsolver import NoConvergenceError

    with pytest.raises(NoConvergenceError):
        _run_overshoot(False)

    _t, v, stats = _run_overshoot(True)
    assert stats.gmin_rescues > 0, 'the run completed without the chain'
    assert stats.forced_nonconverged_steps == 0
    assert np.all(np.isfinite(v))


def test_the_rescued_waveform_is_right_not_merely_finite():
    """A rescued point that is not a solution would be worse than a raise.

    Reference: the same circuit with the floor lowered, so the controller
    solves every point itself and no rescue happens.  Agreement is checked
    against the SIGNAL SWING, because a 5 V circuit resolved at ten points
    per period carries real truncation error -- the claim being protected
    is "the rescued point is the circuit's solution", not "the coarse run
    is accurate".
    """
    t_res, v_res, s_res = _run_overshoot(True)
    t_ref, v_ref, s_ref = _run_overshoot(False, minstep=1e-15)
    assert s_res.gmin_rescues > 0 and s_ref.gmin_rescues == 0

    swing = float(v_ref.max() - v_ref.min())
    dev = float(np.max(np.abs(v_res - np.interp(t_res, t_ref, v_ref))))
    assert dev < 0.02 * swing, 'rescued waveform off by %.3e V of %.3f V' % (
        dev, swing)


# ---------------------------------------------------------------------------
# The chain layered on PCNR -- where the CPU has always had it
# ---------------------------------------------------------------------------

def _diode_chain(nd=2, va=5.0, freq=2e6, rs=1e3):
    """A series diode string driven hard and fast.

    Two junctions in series with tiny node capacitances: the limiter alone
    walks into points it cannot leave, because limiting each junction says
    nothing about the STRING.  This is what is left for a conductance
    homotopy once limiting is already in place.
    """
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=va, freq=freq)
    c['R'] = R('a', 'n0', r=rs)
    for i in range(nd):
        c['D%d' % i] = Diode('n%d' % i, 'n%d' % (i + 1))
        c['C%d' % i] = C('n%d' % (i + 1), gnd, c=1e-13)
    c['RL'] = R('n%d' % nd, gnd, r=1e5)
    return c


_PTS = 2e-7


def _run_chain(cont, minstep=_PTS):
    def go():
        tran = JAXTransient(_diode_chain(), pcnr=True, continuation=cont,
                            reltol=1e-6, minstep=minstep, firststep=_PTS)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=20 * _PTS, timestep=_PTS, uic=True)
        return (np.asarray(res.sweep_values, float),
                np.asarray(res.v('n0'), float).reshape(-1), tran.statistics)
    return _with_jax(go)


def test_the_chain_completes_a_pcnr_run_pcnr_alone_cannot():
    """The point of LAYERING it on PCNR rather than the plain Newton.

    On the plain Newton the chain is nearly useless for a junction
    circuit: that solver does not limit, so the failure is the diode
    exponential overflowing, and no conductance in parallel with a
    junction repairs a term that is already infinite.  PCNR limits, so it
    is the analogue of the CPU's base solver -- and gmin stepping on top
    of a limited Newton is the combination SPICE has always shipped.
    """
    from pycircuit.circuit.nrsolver import NoConvergenceError

    with pytest.raises(NoConvergenceError):
        _run_chain(False)

    _t, v, stats = _run_chain(True)
    assert stats.gmin_rescues > 0
    assert stats.forced_nonconverged_steps == 0
    assert np.all(np.isfinite(v))


def test_pcnr_rescued_waveform_is_right_not_merely_finite():
    """Same gate as the plain path: a rescued point that is not a solution
    would be worse than the raise it replaced."""
    t_res, v_res, s_res = _run_chain(True)
    t_ref, v_ref, s_ref = _run_chain(False, minstep=1e-15)
    assert s_res.gmin_rescues > 0 and s_ref.gmin_rescues == 0

    swing = float(v_ref.max() - v_ref.min())
    dev = float(np.max(np.abs(v_res - np.interp(t_res, t_ref, v_ref))))
    assert dev < 0.02 * swing, 'rescued waveform off by %.3e V of %.3f V' % (
        dev, swing)


def test_vector_pcnr_traces_with_the_chain_armed():
    """The chain must also compile into the DEVICE-view PCNR path, and
    must not perturb it.

    Routed separately from the pair view (sec. 49), so an armed chain that
    only ever traced through the scalar path would be untested here --
    and the vector path is the one every real compact model takes.
    """
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit import elements_hdl as _eh
    from pycircuit.circuit.elements import VS
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit import numeric

    def build():
        c = SubCircuit(toolkit=jaxtoolkit)
        nd, ng, nv = c.add_node('d'), c.add_node('g'), c.add_node('vdd')
        c['vg'] = VSin(ng, gnd, va=0.4, vo=1.1, freq=1e6, toolkit=jaxtoolkit)
        c['vd'] = VS(nv, gnd, v=1.8, toolkit=jaxtoolkit)
        c['Rd'] = R(nv, nd, r=1e4, toolkit=jaxtoolkit)
        c['Cd'] = C(nd, gnd, c=1e-13, toolkit=jaxtoolkit)
        c['M'] = _eh.EkvNmosHdl(nd, ng, gnd, gnd, toolkit=jaxtoolkit)
        return c

    saved = _cm.default_toolkit
    _cm.default_toolkit = jaxtoolkit
    try:
        out = {}
        for cont in (False, True):
            tran = JAXTransient(build(), pcnr=True, reltol=1e-8,
                                continuation=cont)
            meta, _vt = tran._pcnr_setup()
            assert meta[0] == 'vector', 'test no longer covers the vector path'
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = tran.solve(refnode=gnd, tend=1e-6, timestep=2e-8,
                                 uic=True)
            out[cont] = (np.asarray(res.v('d'), float).reshape(-1),
                         tran.statistics.gmin_rescues)
    finally:
        _cm.default_toolkit = saved

    assert out[True][1] == 0 and out[False][1] == 0
    assert np.array_equal(out[True][0], out[False][0]), \
        'the armed chain perturbed the vector PCNR path'


# ---------------------------------------------------------------------------
# The ladder must actually search the range it declares
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('k', [2, 0, -1, -2, -4, -8, -12])
def test_the_ladder_reaches_every_decade_it_declares(k):
    """`e_end` is a contract, and the driver used to break it.

    With no rung yet landed the failure path ESCALATED toward `e_max` and,
    once that was spent, REFINED just below `e_start` -- halving from
    `step/2` and giving up at `min_step`.  So the exponents it could ever
    try were

        {e_start, e_start+step, ... e_max}  and
        {e_start-1, e_start-0.5, e_start-0.25}

    and nothing below, however low `e_end` was set.  Measured with the rung
    here: the ladder landed for k >= -1 and failed for every k <= -2,
    against an `e_end` of -12 advertising a search down to 1e-12.

    That mattered because the ladder's whole job is to find SOME rung that
    converges and then anneal from it; a driver that cannot look below one
    decade under its start is not searching, it is sampling.
    """
    def rung(xs, g):
        ## converges only once the deformation is weak enough, plus the
        ## pure rung (g == 0), which is the solution the ladder must land
        conv = g <= 10.0 ** k
        return xs + jnp.where(conv, 1.0, 0.0), conv

    _x, converged = _adaptive_ladder_traced(rung, jnp.zeros(1),
                                            e_start=0.0, e_end=-12.0,
                                            e_max=6.0)
    assert bool(converged), \
        'the ladder could not reach a rung at g = 1e%d, inside its own ' \
        'declared range' % k


def test_the_descent_does_not_run_when_a_rung_lands():
    """The new phase must be inert on the ordinary path.

    Descending exists for "nothing has landed anywhere".  A ladder whose
    first rung converges must still march down and finish on the pure rung
    exactly as before -- the counter proves the rung count did not blow up.
    """
    calls = []

    def rung(xs, g):
        calls.append(1)
        return xs + 1.0, jnp.asarray(True)

    x, converged = _adaptive_ladder_traced(rung, jnp.zeros(1), e_start=-2.0,
                                           e_end=-12.0)
    assert bool(converged)
    ## -2 to -12 in steps of 2, then the pure rung: 7 traced rung builds.
    ## The assertion is on the ORDER, not the exact number -- a descent
    ## wrongly triggered here would multiply it.
    assert len(calls) <= 10, 'the ordinary march grew to %d rungs' % len(calls)
