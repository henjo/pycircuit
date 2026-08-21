"""Regression tests for the Phase-0 fixes of doc/transient_review_260820.md.

One test (or small group) per finding, appended in the order the fixes landed
so each commit carries the test that fails before it and passes after it.
"""

import numpy as np
import pytest

from pycircuit.circuit import gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, IS
from pycircuit.circuit.transient import Transient


def _rc():
    c = SubCircuit()
    c['is'] = IS(gnd, 'a', i=1e-3)
    c['R'] = R('a', gnd, r=1e3)
    c['C'] = C('a', gnd, c=1e-9)
    return c


## F8 / F18 -- dead arguments deleted; passing them must now fail loudly
## rather than being silently discarded.

def test_f8_analytical_eh_is_gone():
    with pytest.raises(TypeError):
        Transient(_rc()).solve(tend=1e-6, timestep=1e-8, analytical_eh=True)


def test_f18_init_irefnode_is_gone():
    ## Not TypeError: the deleted argument now falls into **kvargs, where the
    ## parameter dictionary rejects it by name -- louder and more specific.
    with pytest.raises(KeyError, match='irefnode'):
        Transient(_rc(), irefnode=0)


def test_f18_get_diff_method_is_gone():
    tran = Transient(_rc())
    with pytest.raises(TypeError):
        tran.get_diff(np.zeros(2), np.zeros((2, 2)), method='euler')


## F9 -- the continuation wrappers must forward a caller-selected linear
## solver; all six inner call sites used to drop it silently.

class _CountingLinsolver:
    def __init__(self):
        self.calls = 0

    def solve(self, A, b, toolkit):
        self.calls += 1
        return toolkit.linearsolver(A, b)


def test_f9_continuation_forwards_linsolver():
    from pycircuit.circuit.nrsolver import GminSteppingNewton, StandardNewton
    from pycircuit.circuit import numeric

    ls = _CountingLinsolver()
    ## A trivially convergent 1-unknown linear system: F(x) = x - 1.
    def eval_FJ(x):
        return np.asarray(x) - 1.0, np.eye(1)

    x, _iters = GminSteppingNewton(StandardNewton()).solve_system(
        np.zeros(1), eval_FJ, numeric, reltol=1e-6,
        abstol=np.array([1e-12]), xtol=np.array([1e-12]), maxiter=20,
        linsolver=ls)
    assert abs(float(x[0]) - 1.0) < 1e-9
    assert ls.calls > 0          # was 0: linsolver dropped at every call site


## F2 -- an override for a class outside the vmap evaluation groups must
## refuse loudly (F2(a)); and for batchable classes the override must
## actually differentiate the lanes (F2(b) -- R was the headline offender:
## its override was silently ignored and every lane returned one waveform).

def test_f2a_batched_override_of_unbatchable_class_raises():
    jax = pytest.importorskip('jax')
    import jax.numpy as jnp
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=1.0)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        tran = JAXTransient(cir)
        ## VS carries a branch unknown and no eval_*_pure -- still honestly
        ## unbatchable.  (R moved to the behavioral test below when F2(b)
        ## made it batchable.)
        with pytest.raises(NotImplementedError, match="'VS'"):
            tran.solve_batched(
                gnd, override_params_tree={'VS': {'v': jnp.array([[1.0], [2.0]])}},
                tend=1e-5, timestep=1e-6, CHUNK_SIZE=50, uic=True)
    finally:
        circuit_mod.default_toolkit = saved


def test_f2b_batched_r_override_differentiates_lanes():
    jax = pytest.importorskip('jax')
    import warnings
    import jax.numpy as jnp
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=1.0)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = JAXTransient(cir, reltol=1e-6).solve_batched(
                gnd, override_params_tree={'R': {'r': jnp.array([[1e2], [1e4]])}},
                tend=5e-3, timestep=1e-4, uic=True)
        v_fast = np.asarray(res[0].v('out'), float).reshape(-1)
        v_slow = np.asarray(res[1].v('out'), float).reshape(-1)
        ## tau = 1e-4 settles inside tend; tau = 1e-2 cannot: the lanes must
        ## end far apart, at their analytic levels.
        assert abs(v_fast[-1] - 1.0) < 1e-3
        assert abs(v_slow[-1] - (1.0 - np.exp(-0.5))) < 1e-2
    finally:
        circuit_mod.default_toolkit = saved


## F6(a) -- the JAX Newton residual is a vector of KCL currents, so its
## absolute floor is iabstol.  vabstol was threaded instead, invisible while
## both defaulted to 1e-12.  Pin the flavour by its observable asymmetry: a
## huge vabstol must leave the waveform untouched, a huge iabstol must
## loosen Newton visibly.  (Node rows only -- F6(b)'s per-row criterion is
## the complete fix.)

def _jax_rc_run(**tol):
    pytest.importorskip('jax')
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=1.0)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        res = JAXTransient(cir, **tol).solve(gnd, tend=2e-4, timestep=1e-5,
                                             uic=True)
        return np.asarray(res.x)
    finally:
        circuit_mod.default_toolkit = saved


def test_f6_newton_tolerances_are_per_row_and_both_live():
    ## RE-DERIVED at F6(b): with per-row criteria BOTH flavours are live --
    ## iabstol floors the node (KCL) rows, vabstol the branch (KVL) rows --
    ## so the interim F6(a) assertion (vabstol bit-identical) no longer
    ## holds, correctly.  Loosening either flavour to absurdity must now
    ## visibly degrade the run; the defaults must stay accurate.
    ## Behavioral loosening probes are UNFALSIFIABLE here: with the clamp
    ## off (F16), one full Newton step from any start lands exactly on a
    ## linear circuit's solution, so even all-floors-huge commits the right
    ## answer.  Enforcement is proven behaviorally elsewhere -- gate 9(e)
    ## raises at maxiter=1 because the per-row pair genuinely fails, and the
    ## F16 test reaches a 400 V rail to 1e-6 -- so what this test pins is
    ## the flavour CONSTRUCTION, the thing F6 was actually about: currents
    ## floor the node rows of the residual test, volts the branch rows, and
    ## the update test is the transpose.
    x_default = _jax_rc_run()
    assert np.all(np.isfinite(x_default))

    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=1.0)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        tran = JAXTransient(cir, iabstol=3.0, vabstol=7.0)
        n = cir.n
        n_nodes = len(cir.nodes)
        ab = np.asarray(tran._newton_abstol(n))
        xt = np.asarray(tran._newton_xtol(n))
        assert np.all(ab[:n_nodes] == 3.0) and np.all(ab[n_nodes:] == 7.0), \
            'residual floor: iabstol on node rows, vabstol on branch rows'
        assert np.all(xt[:n_nodes] == 7.0) and np.all(xt[n_nodes:] == 3.0), \
            'update floor is the transposed flavour'
        assert n > n_nodes, 'circuit must actually have a branch row'
    finally:
        circuit_mod.default_toolkit = saved


## F1(c) -- lanes finishing in different chunks used to crash the fill-forward
## with an empty-slice broadcast; and even without the crash, padding
## duplicated timestamps.  Per-lane collection must survive heterogeneous
## lanes with strictly increasing time and a t=0 first point per lane.

def test_f1_batched_lanes_finishing_in_different_chunks():
    pytest.importorskip('jax')
    import jax.numpy as jnp
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=1.0)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        res = JAXTransient(cir, reltol=1e-6).solve_batched(
            gnd, override_params_tree={'C': {'c': jnp.array([[1e-9], [1e-6]])}},
            tend=5e-3, timestep=1e-3, CHUNK_SIZE=10, uic=True)
        assert len(res) == 2
        lengths = []
        for r in res:
            t = np.asarray(r.sweep_values, dtype=float).reshape(-1)
            assert t[0] == 0.0                    # per-lane t=0 seed survives
            assert np.all(np.diff(t) > 0)         # no duplicated timestamps
            assert t[-1] <= 5e-3 * (1 + 1e-9)
            lengths.append(len(t))
        ## the 1nF and 1uF lanes must genuinely differ -- the heterogeneity
        ## that used to trigger the crash
        assert lengths[0] != lengths[1]
    finally:
        circuit_mod.default_toolkit = saved


## F12(a) -- the CPU result must include the t=0 point (SPICE convention,
## JAX parity), on both the standard and the coupled path.

def test_f12_cpu_result_includes_t0_standard_and_coupled():
    import warnings
    for coupled in (False, True):
        tran = Transient(_rc())
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(tend=1e-6, timestep=1e-8, coupled_lte=coupled)
        t = np.asarray(res.sweep_values, dtype=float).reshape(-1)
        assert t[0] == 0.0, 'coupled=%s missing t=0' % coupled
        assert np.all(np.diff(t) > 0)
        ## and the t=0 state is the operating point, not garbage:
        x0 = np.asarray(res.x)[:, 0]
        assert np.all(np.isfinite(x0))


## F4 -- provided_function has ONE contract: an extra source term pf(t), on
## every path.  The standard path used to call it post-solve as pf(f, J, C)
## and discard the result.

def test_f4_provided_function_is_an_extra_source_on_every_path():
    import warnings

    def run(coupled, pf):
        c = SubCircuit()
        c['is'] = IS(gnd, 'a', i=1e-3)
        c['R'] = R('a', gnd, r=1e3)
        c['C'] = C('a', gnd, c=1e-9)
        tran = Transient(c, uic=True)      # uic: the seed knows nothing of pf
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(tend=2e-6, timestep=1e-8, coupled_lte=coupled,
                             provided_function=pf)
        return float(np.asarray(res.v('a'), dtype=float).reshape(-1)[-1])

    n = 2  # one node + one branch? (sized below from the circuit)
    ## pf injects a constant extra current into node 'a'; called with t only.
    seen_arity = []

    def make_pf(ncirc, idx):
        def pf(*args):
            seen_arity.append(len(args))
            u = np.zeros(ncirc)
            u[idx] = 1e-3          # doubles the source current
            return u
        return pf

    probe = SubCircuit()
    probe['is'] = IS(gnd, 'a', i=1e-3)
    probe['R'] = R('a', gnd, r=1e3)
    probe['C'] = C('a', gnd, c=1e-9)
    idx_a = probe.get_node_index('a')
    ncirc = probe.n

    v_std = run(False, make_pf(ncirc, idx_a))
    v_cpl = run(True, make_pf(ncirc, idx_a))
    v_ref = run(False, None)

    assert set(seen_arity) == {1}, 'pf must be called as pf(t) everywhere'
    ## the extra source moves the settled level, identically on both paths
    assert abs(v_std - v_ref) > 0.1
    assert v_std == pytest.approx(v_cpl, rel=1e-3)


def test_f4_inconsistent_seed_warns():
    import warnings as w
    c = _rc()
    pf = lambda t: np.zeros(c.n)
    with w.catch_warnings(record=True) as caught:
        w.simplefilter('always')
        Transient(c).solve(tend=1e-7, timestep=1e-9, provided_function=pf)
    assert any('does not see' in str(x.message) for x in caught)


## F7 -- the PCNR path must pin the CALLER'S reference node, not gnd.
## Before: solve() never forwarded refnode to solve_timestep, whose gnd
## default fed _solve_timestep_pcnr while step control stripped the caller's
## row -- two reference nodes in one solve.

def test_f7_pcnr_honours_refnode():
    import warnings
    from pycircuit.circuit.elements import Diode, VSin

    def rectifier():
        c = SubCircuit()
        c['vs'] = VSin('a', gnd, va=5.0, freq=1e3)
        c['D'] = Diode('a', 'b')
        c['R'] = R('b', gnd, r=1e3)
        c['C'] = C('b', gnd, c=1e-7)
        return c

    for pcnr in (False, True):
        c = rectifier()
        ib = c.get_node_index('b')
        tran = Transient(c, pcnr=pcnr, reltol=1e-5)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(c.get_node('b'), tend=2e-4, timestep=1e-5)
        xb = np.asarray(res.x)[ib]
        assert float(np.max(np.abs(xb))) == 0.0, \
            'pcnr=%s: requested reference node is not pinned' % pcnr


## F3 -- the coupled path must land on tend exactly; the LTE equation used to
## grow the tend-truncated final step past it (13-24% of the last step on
## quiet tails, 5 of 6 measured configurations).

@pytest.mark.parametrize('tend,ts', [(1e-3, 1e-5), (5e-4, 1e-5), (1e-3, 1e-6)])
def test_f3_coupled_lands_on_tend(tend, ts):
    import warnings
    c = SubCircuit()
    c['is'] = IS(gnd, 'a', i=1e-3)
    c['R'] = R('a', gnd, r=1e3)
    c['C'] = C('a', gnd, c=1e-8)          # tau = 10us: quiet tail long before tend
    tran = Transient(c)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=tend, timestep=ts, coupled_lte=True)
    t = np.asarray(res.sweep_values, dtype=float).reshape(-1)
    assert t[-1] <= tend * (1 + 1e-12)
    assert t[-1] >= tend * (1 - 1e-9)     # and it actually reaches tend


## F5 -- the band parameters' 'auto' sentinel resolves per path, and every
## DOCUMENTED explicit value (0.0 / 1.0 / None included) is honoured
## verbatim.  Before: `par.lte_gamma_min or 0.7` silently replaced an
## explicit 0.0 with 0.7 on the coupled path.

def test_f5_coupled_band_honours_explicit_values():
    tran = Transient(_rc(), lte_gamma_min=0.0, lte_gamma_max=1.0, lte_eta=None)
    assert tran._coupled_band() == (0.0, 1.0, None)


def test_f5_coupled_band_defaults_when_unset():
    assert Transient(_rc())._coupled_band() == (0.7, 3.0, 0.15)


def test_f5_standard_band_defaults_when_unset():
    import warnings
    tran = Transient(_rc())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tran.solve(tend=1e-7, timestep=1e-9)
    ctrl = tran.step_controller
    assert (ctrl.lte_gamma_min, ctrl.lte_gamma_max, ctrl.lte_eta) \
        == (0.0, 1.0, None)


## F14 -- lower-band growth retries (a step redone LARGER because it was too
## accurate) must not trip the force-accept escape hatch: a quiescent circuit
## with a band used to emit "error still above tolerance" warnings during its
## opening ramp, about steps whose error was below the band.

def test_f14_growth_retries_do_not_force_accept():
    import warnings as w
    c = SubCircuit()
    c['is'] = IS(gnd, 'a', i=1e-3)
    c['R'] = R('a', gnd, r=1e3)
    c['C'] = C('a', gnd, c=1e-8)
    tran = Transient(c, lte_gamma_min=0.5, lte_gamma_max=3.0)
    with w.catch_warnings(record=True) as caught:
        w.simplefilter('always')
        res = tran.solve(tend=1e-3, timestep=1e-5)
    spurious = [x for x in caught if 'still above tolerance' in str(x.message)]
    assert not spurious, [str(x.message) for x in spurious]
    assert res.statistics.force_accepts == 0
    ## the band's growth retries are real re-solves and stay in the accounting
    assert res.statistics.rejected_steps > 0


## F13 -- the coupled path's statistics must be complete and attached to the
## result; the JAX path attaches too (parity).

def test_f13_coupled_statistics_complete_and_attached():
    import warnings
    from pycircuit.circuit.elements import VPulse
    ## Same drive/timestep proportions as test_coupled_method's _pulse_run --
    ## a pulse the coupled path is known to complete.
    c = SubCircuit()
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-5, tr=1e-6, tf=1e-6,
                     pw=2e-5, per=4e-5)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    tran = Transient(c, reltol=1e-5)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=6e-5, timestep=1e-6, coupled_lte=True)
    st = res.statistics                      # attached to the RESULT
    assert st.total_seconds > 0
    assert st.solve_seconds > 0
    assert st.accepted_steps > 0
    assert st.breakpoints_hit > 0
    assert st.force_accepts == 0             # by design on this path


def test_f13_jax_statistics_attached_to_result():
    pytest.importorskip('jax')
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=1.0)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        res = JAXTransient(cir).solve(gnd, tend=1e-4, timestep=1e-5, uic=True)
        assert res.statistics.accepted_steps > 0
    finally:
        circuit_mod.default_toolkit = saved


## F19 -- the JAX error estimator follows the EFFECTIVE integration order:
## a step integrated at order 1 (opening steps, breakpoint sentinel) must be
## scored by the order-1 formula with the order-1 exponent, not Gear-2's.

def test_f19_estimator_follows_effective_order():
    pytest.importorskip('jax')
    import jax.numpy as jnp
    from pycircuit.circuit.jaxtransient import (TransientState,
                                                ywr_error_ratio,
                                                effective_first_order)

    def state(h0, h1, forced=False):
        n = 2
        return TransientState(
            t=0.0, dt=1.0, step_idx=5, x_history=None,
            q_history=jnp.zeros((3, n)),
            iq_history=jnp.array([[0.5, 0.0], [0.25, 0.0], [0.0, 0.0]]),
            h_history=jnp.array([h0, h1, 0.0]),
            results_buffer=None, time_buffer=None,
            tline_history=None, tline_head=None,
            force_first_order=jnp.asarray(forced))

    i_curr = jnp.array([1.0, 0.0])
    x_curr = jnp.array([1.0, 0.5])
    J = jnp.eye(2)

    ## Two completed steps recorded, no drop flag: gear scoring, order_p1 = 3.
    _, p_gear = ywr_error_ratio(i_curr, x_curr, J, state(1.0, 1.0),
                                irefnode=1, method='gear')
    assert float(p_gear) == 3.0

    ## Order-1 cases, all through the SAME 'gear' method string: fewer than
    ## two completed steps (run-global h_history facts, chunk-independent),
    ## and the F11 forced-drop flag.
    for st in (state(0.0, 0.0), state(1.0, 0.0), state(1.0, 1.0, forced=True)):
        assert bool(effective_first_order(st))
        _, p = ywr_error_ratio(i_curr, x_curr, J, st, irefnode=1,
                               method='gear')
        assert float(p) == 2.0


def test_f19_forced_lte_counter_exists_and_is_zero_on_clean_runs():
    pytest.importorskip('jax')
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=1.0)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        tran = JAXTransient(cir)
        tran.solve(gnd, tend=1e-4, timestep=1e-5, uic=True)
        assert tran.statistics.forced_lte_steps == 0
    finally:
        circuit_mod.default_toolkit = saved


## F11 -- a step landing on a breakpoint marks the next step first-order, so
## no 2nd-order polynomial is differenced across the corner.  Before, every
## edge cost a rejection burst: 38 rejections / 183 accepted on this circuit;
## after, 16 / 190.  The ratio bound distinguishes the two regimes.

def test_f11_breakpoint_order_drop_tames_edge_rejections():
    pytest.importorskip('jax')
    import warnings
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VPulse
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VPulse('in', gnd, v1=0.0, v2=1.0, td=1e-4, tr=1e-5,
                           tf=1e-5, pw=5e-4, per=1e-3)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        tran = JAXTransient(cir, reltol=1e-4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=3e-3, timestep=2e-5, uic=True)
        st = tran.statistics
        assert st.accepted_steps > 100
        assert st.rejected_steps <= 0.15 * st.accepted_steps, \
            'edge rejection burst is back: %d rejected / %d accepted' \
            % (st.rejected_steps, st.accepted_steps)
        ## and the waveform still tracks the CPU Gear2 reference across the
        ## pulse train (measured 1.2 mV max deviation at landing; bound 10x)
        t = np.asarray(res.sweep_values, float)
        v = np.asarray(res.v('out'), float).reshape(-1)
    finally:
        circuit_mod.default_toolkit = saved

    from pycircuit.circuit import numeric
    from pycircuit.circuit.integrator import Gear2Integrator
    from pycircuit.circuit.elements import VPulse
    c2 = SubCircuit()
    c2.add_node('in'); c2.add_node('out')
    c2['V1'] = VPulse('in', gnd, v1=0.0, v2=1.0, td=1e-4, tr=1e-5,
                      tf=1e-5, pw=5e-4, per=1e-3)
    c2['R1'] = R('in', 'out', r=1e3)
    c2['C1'] = C('out', gnd, c=1e-6)
    tran_c = Transient(c2, toolkit=numeric, integrator=Gear2Integrator(),
                       reltol=1e-4, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res_c = tran_c.solve(tend=3e-3, timestep=2e-5)
    tc = np.asarray(res_c.sweep_values, float)
    vc = np.asarray(res_c.v('out'), float).reshape(-1)
    dev = float(np.max(np.abs(vc - np.interp(tc, t, v))))
    assert dev < 0.012, 'JAX pulse waveform drifted from CPU: %.3e' % dev


## F17 -- the JAX step prediction aims at safety**p, not at the rejection
## edge.  Edge-aiming cost a 44%/40% rejection rate on rc-vsin at reltol
## 1e-4/1e-6; the safety margin cuts it to 20%/4%.  Bounds sit between the
## two regimes.

def test_f17_safety_factor_tames_rejection_rate():
    pytest.importorskip('jax')
    import warnings
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VSin
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        for reltol, bound in ((1e-4, 0.30), (1e-6, 0.10)):
            c = SubCircuit()
            c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
            c['R'] = R('a', 'b', r=1e3)
            c['C'] = C('b', gnd, c=1e-7)
            tran = JAXTransient(c, reltol=reltol)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tran.solve(gnd, tend=5e-3, timestep=1e-4, uic=True)
            st = tran.statistics
            ratio = st.rejected_steps / max(1, st.accepted_steps)
            assert ratio < bound, \
                'reltol %g: rejection ratio %.3f (edge-aiming regime was ~0.4)' \
                % (reltol, ratio)
    finally:
        circuit_mod.default_toolkit = saved


## F16 -- the Newton update clamp is a parameter, default OFF.  The hardcoded
## 0.5 V made any swing beyond ~maxiter*0.5 V non-convergent by construction:
## a 400 V uic start needed 800 clamped iterations against maxiter=100.

def test_f16_high_voltage_converges_with_clamp_off():
    pytest.importorskip('jax')
    import warnings
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        cir = SubCircuit()
        cir.add_node('in'); cir.add_node('out')
        cir['V1'] = VS('in', gnd, v=400.0)
        cir['R1'] = R('in', 'out', r=1e3)
        cir['C1'] = C('out', gnd, c=1e-6)
        tran = JAXTransient(cir)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=1e-4, timestep=1e-5, uic=True)
        assert tran.statistics.nonconverged_steps == 0
        assert tran.statistics.forced_nonconverged_steps == 0
        v = np.asarray(res.v('in'), float).reshape(-1)
        assert abs(v[-1] - 400.0) < 1e-6      # the rail is actually reached
    finally:
        circuit_mod.default_toolkit = saved


## Hygiene (Phase 3 item 22) -- collect_breakpoints is bounded per element:
## a periodic source over a long tend used to materialise millions of list
## entries before the first step ran.

def test_breakpoint_scan_is_bounded_per_element():
    pytest.importorskip('jax')
    import warnings as w
    from pycircuit.circuit import circuit as circuit_mod
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import collect_breakpoints
    from pycircuit.circuit.elements import VPulse
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        c = SubCircuit()
        ## ~5 events/us over 1 s = ~5M events without the cap.  (VSin's
        ## next_event returns inf on this path -- pulses are the event
        ## generators here.)
        c['vs'] = VPulse('a', gnd, v1=0, v2=1, td=1e-7, tr=1e-7, tf=1e-7,
                         pw=4e-7, per=1e-6)
        c['R'] = R('a', gnd, r=1e3)
        with w.catch_warnings(record=True) as caught:
            w.simplefilter('always')
            bps = collect_breakpoints(c, 1.0)
        assert len(bps) <= 10001                     # cap + tend
        assert bps[-1] == 1.0                        # tend always present
        assert any('breakpoints by t=' in str(x.message) for x in caught)
    finally:
        circuit_mod.default_toolkit = saved


## F15 -- a failed backtracking search takes the LAST-TRIED (smallest) step
## with its own residual, never the full step with a stale one.  Probe: a
## residual landscape where every step from the minimum increases |F|, so the
## search always fails; the solver must creep (least-bad damped steps), not
## leap (full steps), and must report non-convergence honestly.

def test_f15_damped_newton_failed_search_takes_smallest_step():
    from pycircuit.circuit.nrsolver import DampedNewton
    from pycircuit.circuit import numeric
    from pycircuit.circuit.analysis import NoConvergenceError

    evals = []

    def eval_FJ(x):
        evals.append(float(x[0]))
        return np.array([1.0 + 100.0 * abs(x[0])]), np.array([[1.0]])

    with pytest.raises(NoConvergenceError):
        DampedNewton().solve_system(
            np.zeros(1), eval_FJ, numeric, reltol=1e-3,
            abstol=np.array([1e-12]), xtol=np.array([1e-6]), maxiter=5)

    ## The Jacobian lies about the slope, so even damped iterates grow --
    ## but the growth rate is the discriminator: full steps multiply |x| by
    ## ~101 per outer iteration (|x| ~ 1e8 after 5), the least-bad 0.0625
    ## step by ~6.3 (measured max ~2.8e3).  1e5 separates the regimes.
    assert max(abs(v) for v in evals) < 1e5, \
        'iterates ran away -- the failed search took full steps again: %r' \
        % evals[-3:]


## F10 -- PIController honours the band it accepts: gamma_max moves the
## rejection threshold, eta damps, and the unimplemented lower band is
## refused loudly instead of ignored.

def _pi_evaluate(pi, err_value):
    """Drive evaluate_step with a stub integrator so err comes out exactly
    err_value (J = I, TRTOL = 1, reltol = 0, abstol = 1)."""
    from pycircuit.circuit import numeric

    class StubIntegrator:
        ORDER = 1

        def compute_lte(self, **kw):
            return np.array([err_value, 0.0]), 2.0

    return pi.evaluate_step(
        x_curr=np.array([1.0, 0.0]), x_last=np.array([1.0, 0.0]),
        q_curr=None, q_last_hist=None, iq_last_hist=None,
        h_curr=1e-6, h_last=1e-6, no_history=False,
        J=np.eye(2), active_integrator=StubIntegrator(), irefnode=1,
        reltol=0.0, abstol=np.ones(2), toolkit=numeric,
        max_step=1e-3, TRTOL=1.0)


def test_f10_pi_honours_gamma_max():
    from pycircuit.circuit.stepcontroller import PIController
    ## err = 2.0: rejected at the historical threshold, accepted with the
    ## band's upper edge at 3.0 -- which was silently ignored before.
    accept, _ = _pi_evaluate(PIController(), 2.0)
    assert not accept
    accept, _ = _pi_evaluate(PIController().set_lte_band(gamma_max=3.0), 2.0)
    assert accept


def test_f10_pi_honours_eta_damper():
    from pycircuit.circuit.stepcontroller import PIController
    pi = PIController().set_lte_band(gamma_max=3.0, eta=0.1)
    accept, h_next = _pi_evaluate(pi, 1e-6)   # tiny error wants huge growth
    assert accept
    assert h_next <= 1e-6 * 1.1 * (1 + 1e-12)   # damped to +10%


def test_f10_pi_refuses_the_lower_band():
    from pycircuit.circuit.stepcontroller import PIController
    with pytest.raises(NotImplementedError, match='lower LTE band'):
        PIController().set_lte_band(gamma_min=0.5, gamma_max=3.0)
    ## the 'auto' sentinel (the shipped default) must pass untouched
    PIController().set_lte_band('auto', 'auto', 'auto')
