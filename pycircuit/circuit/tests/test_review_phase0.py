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


## F2(a) -- an override for a class outside the vmap evaluation groups was
## silently ignored (N bit-identical lanes presented as N samples); it must
## refuse loudly until the class is batchable.

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
        with pytest.raises(NotImplementedError, match="'R'"):
            tran.solve_batched(
                gnd, override_params_tree={'R': {'r': jnp.array([[1.0], [1e3]])}},
                tend=1e-5, timestep=1e-6, CHUNK_SIZE=50, uic=True)
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


def test_f6a_newton_residual_floor_is_current_flavoured():
    x_default = _jax_rc_run()
    ## vabstol no longer reaches the Newton criterion: bit-identical run.
    assert np.array_equal(x_default, _jax_rc_run(vabstol=1e6))
    ## iabstol does: a huge floor accepts the predictor and the waveform moves
    ## (measured 4.48 V max deviation on this circuit when the test was written).
    x_loose = _jax_rc_run(iabstol=1e6)
    assert x_default.shape != x_loose.shape or \
        float(np.max(np.abs(x_default - x_loose))) > 1e-3


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
