"""Phase A of doc/backend_parity_260821.md -- the eight alignment items.

Each test here fails against the pre-Phase-A code by construction:
`solve(uicc=True)` ran silently with defaults (P1), `Transient(TRTOL=...)`
raised KeyError (P2), `solve_batched` took a second differently-named floor
three decades looser (P3), the JAX `refnode` default was the integer 0 (P4),
`outputstep` was a KeyError on JAX (P10), sub-`minbreak` breakpoints were
kept distinct (P14), and `force_accepts` was an AttributeError on the JAX
statistics object (P15).  P5's gate is the dead-knob scan itself: its
allowlist entry for ('JAXTransient', 'analysis') is deleted, so the scan
fails unless the Parameter is actually read.
"""

import warnings

import numpy as np
import pytest

jax = pytest.importorskip('jax')
import jax.numpy as jnp

from pycircuit.circuit import circuit as circuit_mod, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS, VSin, VPulse


def _with_jax(fn):
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        return fn()
    finally:
        circuit_mod.default_toolkit = saved


def _rc():
    c = SubCircuit()
    c['V1'] = VSin('in', gnd, va=1.0, freq=1e3)
    c['R1'] = R('in', 'out', r=1e3)
    c['C1'] = C('out', gnd, c=1e-6)
    return c


def test_p1_typo_raises_instead_of_running_with_defaults():
    """The **kwargs sink is gone: a misspelled argument is a TypeError."""
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        tran = JAXTransient(_rc(), reltol=1e-4)
        with pytest.raises(TypeError):
            tran.solve(gnd, tend=1e-4, timestep=1e-5, uicc=True)
        with pytest.raises(TypeError):
            tran.solve(gnd, tend=1e-4, timestep=1e-5, min_step=1e-12)
    _with_jax(go)


def test_p1_uic_works_as_parameter_and_as_argument():
    """Constructor Parameter and per-call override give the same run."""
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        outs = []
        for kind in ('param', 'arg'):
            tran = (JAXTransient(_rc(), reltol=1e-4, uic=True) if kind == 'param'
                    else JAXTransient(_rc(), reltol=1e-4))
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = (tran.solve(gnd, tend=1e-4, timestep=1e-5)
                       if kind == 'param' else
                       tran.solve(gnd, tend=1e-4, timestep=1e-5, uic=True))
            outs.append(np.asarray(res.v('out'), float).reshape(-1))
        assert np.array_equal(outs[0], outs[1])
    _with_jax(go)


def test_p2_trtol_settable_and_live_on_both_backends():
    """One name for the LTE ratio, on both sides.

    CPU: behavioral -- `Transient(TRTOL=...)` used to raise KeyError; now the
    LTERATIO property follows the Parameter and TRTOL=1 (7x tighter) takes
    measurably more steps (77 -> 154 measured on the step-response RC at
    landing).  JAX: kernel-level -- on the probed configurations the gear-2
    charge LTE sits so far under tolerance that the run is step-cap-limited
    and NO tolerance knob moves the step count (reltol 1e-4 vs 1e-6:
    identical 209-step runs, an order/headroom artifact, not a dead knob), so
    the behavioral assertion would measure headroom rather than plumbing.
    Instead the kernel contract is pinned directly: `ywr_error_ratio` scales
    exactly as 1/TRTOL, which is what makes the Parameter live wherever the
    estimate DOES bind."""
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.elements import VS

    def rc_step():
        c = SubCircuit()
        c['V1'] = VS('in', gnd, v=1.0)
        c['R1'] = R('in', 'out', r=1e3)
        c['C1'] = C('out', gnd, c=1e-6)
        return c

    counts = {}
    for trtol in (7.0, 1.0):
        tran = Transient(rc_step(), toolkit=numeric, reltol=1e-4, uic=True,
                         TRTOL=trtol)
        assert tran.LTERATIO == trtol
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tran.solve(gnd, tend=5e-3, timestep=1e-4)
        counts[trtol] = tran.statistics.accepted_steps
    assert counts[1.0] > counts[7.0], counts

    def go():
        from pycircuit.circuit.jaxtransient import (JAXTransient,
                                                    ywr_error_ratio,
                                                    TransientState)
        ## The Parameter must exist and be threaded (KeyError before P2's
        ## JAX-side declaration; the two outer_time_loop call sites pass it).
        JAXTransient(_rc(), reltol=1e-4, TRTOL=1.0)

        ## Kernel contract: etol = TRTOL * (...), so ratio_1 == 7 * ratio_7.
        ## The history rings are depth 3, as solve builds them.
        n, depth = 3, 3
        rng = np.random.default_rng(42)
        st = TransientState(
            t=jnp.asarray(1e-4), dt=jnp.asarray(1e-5), step_idx=5,
            x_history=jnp.asarray(rng.normal(size=(depth, n))),
            q_history=jnp.asarray(rng.normal(size=(depth, n)) * 1e-9),
            iq_history=jnp.asarray(rng.normal(size=(depth, n)) * 1e-3),
            h_history=jnp.full((depth,), 1e-5),
            results_buffer=jnp.zeros((1, n)), time_buffer=jnp.zeros(1),
            tline_history=jnp.zeros((0, 1, 5)),
            tline_head=jnp.asarray(0, dtype=jnp.int32),
            sig_max=jnp.asarray(1.0))
        i_curr = jnp.asarray(rng.normal(size=n) * 1e-3)
        x_curr = jnp.asarray(rng.normal(size=n))
        J = jnp.eye(n)
        r7, _ = ywr_error_ratio(i_curr, x_curr, J, st, 0, method='gear',
                                trtol=7.0, lte_rel=1e-4, lte_abstol=1e-12,
                                first_order=jnp.asarray(False))
        r1, _ = ywr_error_ratio(i_curr, x_curr, J, st, 0, method='gear',
                                trtol=1.0, lte_rel=1e-4, lte_abstol=1e-12,
                                first_order=jnp.asarray(False))
        assert float(r7) > 0.0
        np.testing.assert_allclose(float(r1), 7.0 * float(r7), rtol=1e-12)
    _with_jax(go)


def test_p3_one_floor_one_vocabulary_for_solve_batched():
    """`solve_batched` has no `dt_min`, its first argument is `refnode`, and
    both entry points read the same `minstep` Parameter."""
    import inspect
    from pycircuit.circuit.jaxtransient import JAXTransient

    sig = inspect.signature(JAXTransient.solve_batched)
    assert 'dt_min' not in sig.parameters
    assert 'irefnode' not in sig.parameters
    assert 'minstep' in sig.parameters
    assert list(sig.parameters)[1] == 'refnode'
    assert 'minstep' in inspect.signature(JAXTransient.solve).parameters

    def go():
        tran = JAXTransient(_rc(), reltol=1e-4, minstep=1e-16)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve_batched(
                refnode=gnd,
                override_params_tree={'R': {'r': jnp.array([[1e2], [1e4]])}},
                tend=1e-4, timestep=1e-5, uic=True)
        assert len(res) == 2
    _with_jax(go)


def test_p4_refnode_defaults_to_gnd_object():
    """A bare solve() must pin gnd, not whichever node holds index 0."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    import inspect
    assert inspect.signature(JAXTransient.solve).parameters['refnode'].default is gnd

    def go():
        outs = []
        for explicit in (False, True):
            tran = JAXTransient(_rc(), reltol=1e-4)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = (tran.solve(gnd, tend=1e-4, timestep=1e-5, uic=True)
                       if explicit else
                       tran.solve(tend=1e-4, timestep=1e-5, uic=True))
            outs.append(np.asarray(res.v('out'), float).reshape(-1))
        assert np.array_equal(outs[0], outs[1])
    _with_jax(go)


def test_p10_outputstep_resamples_the_jax_result():
    """Same contract as the CPU: uniform grid out, quadratic values within
    tolerance of the adaptive run."""
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r_adaptive = JAXTransient(_rc(), reltol=1e-4).solve(
                gnd, tend=2e-3, timestep=1e-5, uic=True)
            r_uniform = JAXTransient(_rc(), reltol=1e-4,
                                     outputstep=5e-5).solve(
                gnd, tend=2e-3, timestep=1e-5, uic=True)
        tu = np.asarray(r_uniform.sweep_values, float)
        vu = np.asarray(r_uniform.v('out'), float).reshape(-1)
        ta = np.asarray(r_adaptive.sweep_values, float)
        va = np.asarray(r_adaptive.v('out'), float).reshape(-1)
        dsteps = np.diff(tu)
        assert np.allclose(dsteps, 5e-5, rtol=1e-9)
        dev = np.max(np.abs(np.interp(tu, ta, va) - vu))
        assert dev < 1e-3, 'resampled values drifted: %.3e' % dev
    _with_jax(go)


def test_p14_breakpoints_closer_than_minbreak_merge():
    """Two pulses offset by less than minbreak yield merged breakpoints, and
    tend always survives the merge verbatim."""
    from pycircuit.circuit.jaxtransient import collect_breakpoints

    def go():
        c = SubCircuit()
        c.add_node('b')
        c['V1'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-6, tr=1e-7,
                         tf=1e-7, pw=1e-5, per=1e-4)
        c['V2'] = VPulse('b', gnd, v1=0.0, v2=1.0, td=1e-6 + 1e-16, tr=1e-7,
                         tf=1e-7, pw=1e-5, per=1e-4)
        c['R1'] = R('a', 'b', r=1e3)
        breaks = collect_breakpoints(c, tend=5e-6, minbreak=1e-14)
        arr = np.asarray(breaks, float)
        assert np.min(np.diff(arr)) > 1e-14, 'sub-minbreak spacing survived'
        assert arr[-1] == 5e-6
        ## The un-merged scan would carry both 1e-6 and 1e-6 + 1e-16.
        near_td = arr[np.abs(arr - 1e-6) < 1e-9]
        assert len(near_td) == 1
    _with_jax(go)


def test_p15_statistics_speak_the_cpu_name():
    """`forced_lte_steps` was the CPU's `force_accepts` under another name --
    same definition, accepted with LTE above tolerance."""
    from pycircuit.circuit.jaxtransient import JAXTransientStatistics
    st = JAXTransientStatistics()
    assert hasattr(st, 'force_accepts')
    assert not hasattr(st, 'forced_lte_steps')
