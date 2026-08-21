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
            ## timestep_max pins the adaptive grid density the 1e-3 bound
            ## was derived at (cap decoupled from timestep 2026-08-21) --
            ## this test measures the resample, not the default cap.
            r_adaptive = JAXTransient(_rc(), reltol=1e-4,
                                      timestep_max=1e-5).solve(
                gnd, tend=2e-3, timestep=1e-5, uic=True)
            r_uniform = JAXTransient(_rc(), reltol=1e-4, timestep_max=1e-5,
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


def test_p6_integrator_defaults_agree_and_the_choice_is_live():
    """P6 (Phase B, recorded here with the other parity gates): the shipped
    standard-path defaults agree at last -- Gear-2 on both -- and the JAX
    `integrator` Parameter actually reaches the traced loop: euler and gear
    produce measurably different runs (order-1 vs order-2 error at the same
    config: measured 8.2e-3 vs 1.1e-3 on the step-response RC), and
    trapezoidal is refused with the reason (the uniform-grid trap branch was
    deleted for cause)."""
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import Gear2Integrator, EulerIntegrator
    from pycircuit.circuit.elements import VS

    tran = Transient(_rc(), toolkit=numeric)
    ## P22: one default for both paths since the state-row mask -- the
    ## coupled= parameter itself is gone (dead-knob scan).
    assert isinstance(tran._get_integrator(), Gear2Integrator)

    def rc_step():
        c = SubCircuit()
        c['V1'] = VS('in', gnd, v=1.0)
        c['R1'] = R('in', 'out', r=1e3)
        c['C1'] = C('out', gnd, c=1e-6)
        return c

    def go():
        from pycircuit.circuit.jaxtransient import JAXTransient
        assert JAXTransient(rc_step()).par.integrator == 'gear'
        errs = {}
        for m in ('gear', 'euler'):
            tran = JAXTransient(rc_step(), reltol=1e-4, integrator=m)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = tran.solve(gnd, tend=5e-3, timestep=1e-4, uic=True)
            t = np.asarray(res.sweep_values, float)
            v = np.asarray(res.v('out'), float).reshape(-1)
            errs[m] = float(np.max(np.abs(v - (1 - np.exp(-t / 1e-3)))))
        assert errs['euler'] > 3.0 * errs['gear'], errs
        with pytest.raises(ValueError, match='CPU-only'):
            JAXTransient(rc_step(), integrator='trap')._eval_method()
    _with_jax(go)


def test_p8_standard_band_on_jax():
    """P8 (Phase B): the CPU standard band's semantics on the traced loop.

    Four facts, each measured at landing: (1) the shipped defaults reduce
    every band line to the historical err<=1.0 accept -- a default run and
    an explicit (0, 1) run are bit-identical; (2) gamma_min fires the
    too-accurate growth-reject (10 rejections, 63 -> 58 accepted on the
    step-response RC at reltol 1e-5); (3) eta clamps the step ratio at
    exactly 1 + eta across the run; (4) an invalid band is refused with the
    CPU's message."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS

    def rc_step():
        c = SubCircuit()
        c['V1'] = VS('in', gnd, v=1.0)
        c['R1'] = R('in', 'out', r=1e3)
        c['C1'] = C('out', gnd, c=1e-6)
        return c

    def run(**kw):
        tran = JAXTransient(rc_step(), reltol=1e-5, **kw)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=5e-3, timestep=1e-3, uic=True)
        return (np.asarray(res.sweep_values, float),
                np.asarray(res.v('out'), float).reshape(-1),
                tran.statistics)

    def go():
        t0, v0, st0 = run()
        _t, v0b, _s = run(lte_gamma_min=0.0, lte_gamma_max=1.0)
        assert np.array_equal(v0, v0b), 'the default band is not inert'

        _t, _v, st1 = run(lte_gamma_min=0.3)
        assert st1.rejected_steps > 0, \
            'gamma_min never fired the too-accurate reject'
        assert st1.accepted_steps < st0.accepted_steps, \
            'growth-rejects did not buy larger steps'

        t3, _v, _s = run(lte_eta=0.10)
        d = np.diff(t3)
        ratio = float(np.max(d[1:] / d[:-1]))
        assert ratio <= 1.10 * (1.0 + 1e-9), \
            'eta damper violated: max step ratio %.4f' % ratio

        with pytest.raises(ValueError, match='gamma_min < gamma_max'):
            run(lte_gamma_min=2.0, lte_gamma_max=1.0)
    _with_jax(go)


def test_p9_fixed_timestep_on_jax():
    """P9 (Phase B): the CPU's stage-4h fixed-grid semantics on the traced
    loop -- the grid wins, breakpoints do not move it, crossing an edge
    still drops the order for the next step, tend is landed on exactly, and
    coupled+fixed is refused rather than approximated.

    The cross-backend agreement here is the strongest gate this file has:
    fixing a float knife-edge in the CPU's crossing test (found by this
    port -- the order drop fired at the first pulse edge and silently
    missed every later one once accumulated t flipped a <=), the two
    backends produce the SAME fixed-grid waveform to 5.6e-17."""
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.integrator import Gear2Integrator
    from pycircuit.circuit.elements import VPulse

    def pc():
        c = SubCircuit()
        c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-6, tr=1e-7,
                         tf=1e-7, pw=1e-6, per=3e-6)
        c['R'] = R('a', 'b', r=1e3)
        c['C'] = C('b', gnd, c=1e-10)
        return c

    cpu = Transient(pc(), toolkit=numeric, uic=True,
                    integrator=Gear2Integrator())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rc_ = cpu.solve(gnd, tend=6e-6, timestep=1e-7, fixed_timestep=True)
    vc = np.asarray(rc_.v('b'), float).reshape(-1)

    def go():
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            rj = JAXTransient(pc()).solve(gnd, tend=6e-6, timestep=1e-7,
                                          uic=True, fixed_timestep=True)
        t = np.asarray(rj.sweep_values, float)
        v = np.asarray(rj.v('b'), float).reshape(-1)
        assert np.allclose(np.diff(t), 1e-7, rtol=1e-9), 'grid not uniform'
        assert t[-1] == pytest.approx(6e-6, rel=1e-12)
        m = min(len(vc), len(v))
        dev = float(np.max(np.abs(vc[:m] - v[:m])))
        assert dev < 1e-12, 'backends disagree on the fixed grid: %.3e' % dev

        with pytest.raises(NotImplementedError, match='grid_locked'):
            JAXTransient(pc(), coupled_lte=True).solve(
                gnd, tend=1e-6, timestep=1e-7, uic=True, fixed_timestep=True)
    _with_jax(go)


def test_p7_relref_modes_on_jax():
    """P7 (Phase B): the CPU's relref modes on the traced estimator.

    Measured at landing on the step-response RC at reltol 1e-5: sigglobal
    63 steps / 6.9e-4, alllocal 72 / 4.9e-4, pointlocal 95 / 3.4e-4 -- the
    tighter the reference, the more steps and the better the answer, in the
    CPU's ordering.  The kernel half pins the unit-group split: with a
    branch current a thousand times the node signals, the mixed (pre-P7)
    scalar reference deadens node error control; the split reference keeps
    the node rows measured against node signals."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VS

    def rc_step():
        c = SubCircuit()
        c['V1'] = VS('in', gnd, v=1.0)
        c['R1'] = R('in', 'out', r=1e3)
        c['C1'] = C('out', gnd, c=1e-6)
        return c

    def go():
        from pycircuit.circuit.jaxtransient import (ywr_error_ratio,
                                                    TransientState)
        import jax.numpy as jnp
        stats = {}
        for mode in ('sigglobal', 'alllocal', 'pointlocal'):
            tran = JAXTransient(rc_step(), reltol=1e-5, relref=mode)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = tran.solve(gnd, tend=5e-3, timestep=1e-3, uic=True)
            t = np.asarray(res.sweep_values, float)
            v = np.asarray(res.v('out'), float).reshape(-1)
            stats[mode] = (len(t), float(np.max(np.abs(v - (1 - np.exp(-t / 1e-3))))))
        assert stats['pointlocal'][0] > stats['alllocal'][0] > stats['sigglobal'][0], stats
        assert stats['pointlocal'][1] < stats['sigglobal'][1], stats

        with pytest.raises(ValueError, match='relref'):
            JAXTransient(rc_step(), relref='bogus').solve(
                gnd, tend=1e-4, timestep=1e-5, uic=True)

        ## Unit-group split, kernel level: rows 0-1 are node voltages at
        ## ~1e-3 V, row 2 a branch current at 1 A.
        n, depth = 3, 3
        rng = np.random.default_rng(7)
        x = jnp.asarray([1e-3, -8e-4, 1.0])
        st = TransientState(
            t=jnp.asarray(1e-4), dt=jnp.asarray(1e-5), step_idx=5,
            x_history=jnp.tile(x * 0.9, (depth, 1)),
            q_history=jnp.asarray(rng.normal(size=(depth, n)) * 1e-9),
            iq_history=jnp.asarray(rng.normal(size=(depth, n)) * 1e-6),
            h_history=jnp.full((depth,), 1e-5),
            results_buffer=jnp.zeros((1, n)), time_buffer=jnp.zeros(1),
            tline_history=jnp.zeros((0, 1, 5)),
            tline_head=jnp.asarray(0, dtype=jnp.int32),
            sig_max=jnp.asarray(1.0), ref_running=jnp.abs(x))
        kw = dict(method='gear', trtol=7.0, lte_rel=1e-4,
                  lte_abstol=jnp.asarray([1e-12, 1e-12, 1e-12]),
                  first_order=jnp.asarray(False), relref='sigglobal')
        r_split, _ = ywr_error_ratio(jnp.zeros(n) + 1e-6, x, jnp.eye(n), st, 0,
                                     n_nodes=2, **kw)
        r_mixed, _ = ywr_error_ratio(jnp.zeros(n) + 1e-6, x, jnp.eye(n), st, 0,
                                     n_nodes=None, **kw)
        ## Mixed: node rows referenced to the 1 A branch -> etol inflated a
        ## thousandfold -> the ratio collapses.  Split restores node control.
        assert float(r_split) > 100.0 * float(r_mixed), (float(r_split),
                                                        float(r_mixed))
    _with_jax(go)


def test_p11_provided_function_on_jax():
    """P11 (Phase B): the F4 contract -- an extra source pf(t) -- on the
    traced loop.  The same sinusoidal injection produces the same waveform
    on both backends; the DC-seed caveat warns with the CPU's wording; and
    a pf must actually change the answer (a gate against silently dropping
    the term, the failure mode F4 exists to prevent)."""
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.integrator import Gear2Integrator

    def rc():
        c = SubCircuit()
        c['R1'] = R('out', gnd, r=1e3)
        c['C1'] = C('out', gnd, c=1e-6)
        return c

    import numpy as _np
    def pf_cpu(t):
        ## 1 mA * sin injected into the 'out' row; the vector is full-length
        ## (n = 2 with gnd) with zero on the reference row.  DELIBERATELY
        ## unbalanced (no return entry on gnd): the reference-row residual
        ## then carries the full injection, which pinned a JAX defect -- its
        ## Newton scored the UNSOLVED reference row and livelocked at
        ## maxiter on every step.  The reduced-rows convergence test is what
        ## this gate regresses.
        return _np.array([-1e-3 * _np.sin(2 * _np.pi * 1e3 * t), 0.0])

    cpu = Transient(rc(), toolkit=numeric, reltol=1e-5, uic=True,
                    integrator=Gear2Integrator(), timestep_max=1e-5)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rc_ = cpu.solve(gnd, tend=1e-3, timestep=1e-5,
                        provided_function=pf_cpu)
    tc = np.asarray(rc_.sweep_values, float)
    vc = np.asarray(rc_.v('out'), float).reshape(-1)
    assert np.max(np.abs(vc)) > 1e-2, 'the CPU pf did nothing'

    def go():
        import jax.numpy as jnp

        def pf_jax(t):
            return jnp.array([-1e-3 * jnp.sin(2 * jnp.pi * 1e3 * t), 0.0])

        tran = JAXTransient(rc(), reltol=1e-5, timestep_max=1e-5)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            rj = tran.solve(gnd, tend=1e-3, timestep=1e-5, uic=True,
                            provided_function=pf_jax)
        tj = np.asarray(rj.sweep_values, float)
        vj = np.asarray(rj.v('out'), float).reshape(-1)
        dev = float(np.max(np.abs(np.interp(tc, tj, vj) - vc)))
        scale = float(np.max(np.abs(vc)))
        assert dev < 0.02 * scale, \
            'pf waveforms disagree across backends: %.3e (scale %.3e)' % (
                dev, scale)

        ## The DC-seed caveat, CPU wording.
        tran2 = JAXTransient(rc(), reltol=1e-4, timestep_max=1e-5)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            tran2.solve(gnd, tend=1e-4, timestep=1e-5,
                        provided_function=pf_jax)
        assert any('DC operating point does not see' in str(x.message)
                   for x in w), 'the DC-seed warning did not fire'
    _with_jax(go)


def test_p16_tline_ring_interpolation_is_limited_quadratic():
    """P16 (Phase B): the delay lookup is the CPU's 3-point limited
    quadratic now.  Kernel-level: on a smoothly curved history the lookup
    reproduces a parabola exactly (linear cannot); at a recorded kink the
    limiter clamps to the bracket instead of overshooting -- the phantom-EMF
    defect the CPU limiter was measured against.  Cross-backend note,
    recorded per house rules: matching the order left the matched-line
    coupled agreement bit-close (2.8e-16) and did NOT move the rc-load
    corner gap (2.66e-2, before and after) -- that gap is controller step
    placement, not interpolation."""
    def go():
        import jax.numpy as jnp
        from pycircuit.circuit.jaxtransient import (_tline_emfs,
                                                    TLINE_HISTORY_DEPTH)

        params = jnp.asarray([1e-9, 50.0])  # TD, Z0

        def hist(rows):
            h = jnp.zeros((TLINE_HISTORY_DEPTH, 5))
            for i, r in enumerate(rows):
                h = h.at[i].set(jnp.asarray(r))
            return h, len(rows) - 1

        ## Parabolic v2(t) = t^2 (scaled), i2 = 0: quadratic recovery exact.
        ts = [0.0, 1.0, 2.0, 3.0]
        rows = [(t, 0.0, t * t, 0.0, 0.0) for t in ts]
        h, head = hist(rows)
        e1, _e2, de1, _de2 = _tline_emfs(2.5, params, h, head)
        assert abs(float(e1) - 2.5 ** 2) < 1e-12, float(e1)
        assert abs(float(de1) - 2 * 2.5) < 1e-12, float(de1)

        ## Kink: flat 0 then ramp -- a quadratic through (ramp, corner,
        ## flat) overshoots below zero at 0.5; the limiter must clamp into
        ## the bracket [0, 0].
        rows = [(0.0, 0.0, 2.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0, 0.0),
                (2.0, 0.0, 0.0, 0.0, 0.0)]
        h, head = hist(rows)
        e1, _e2, _d, _d2 = _tline_emfs(1.5, params, h, head)
        assert abs(float(e1)) < 1e-12, \
            'kink overshoot survived the limiter: %g' % float(e1)
    _with_jax(go)


def test_p17_strategy_objects_are_refused_permanently():
    """Phase C / P17: the traced loop cannot dispatch into per-iteration
    Python strategy objects -- the __init__ refusal is the permanent
    contract, and `linearsolver` (silently swallowed until Phase C) is
    covered along with `nrsolver` and `scaler`."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.linearsolver import DenseSolver

    def go():
        for knob, value in (('nrsolver', object()), ('scaler', object()),
                            ('linearsolver', DenseSolver())):
            with pytest.raises(NotImplementedError, match='permanently'):
                JAXTransient(_rc(), **{knob: value})
    _with_jax(go)


def test_p22_state_row_mask_and_shared_coupled_default():
    """P22: eq (6) measured on STATE rows only, and the coupled Euler
    carve-out retired.

    The mask's proof is the rectifier that livelocked under coupled+Gear-2
    (h collapsed to 6.3e-12 at t=1.25e-4; the source-current row's
    deviation floor of 2.5e-6 A against etol 3.6e-7 was h-independent --
    an ALGEBRAIC row carrying the diode's dq/dt through KCL measures grid
    conventions, not truncation).  With the mask it completes in ~259
    points against Euler's ~769 at the same accuracy (9.7e-3 vs 9.9e-3
    against a fine reference), and the coupled default is Gear-2 on both
    backends.  The cold-start contract survives: plain coupled still
    RAISES (in seconds, not a dt-floor death march -- the forced-accept
    early exit), which pins three in-trace livelock fixes at once."""
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.elements import VSin, Diode
    from pycircuit.circuit.integrator import Gear2Integrator, EulerIntegrator

    def rect():
        c = SubCircuit()
        c['vs'] = VSin('a', gnd, va=5.0, freq=1e3)
        c['D'] = Diode('a', 'b')
        c['R'] = R('b', gnd, r=1e3)
        c['C'] = C('b', gnd, c=1e-7)
        return c

    counts = {}
    for name, integ in (('gear2', Gear2Integrator()),
                        ('euler', EulerIntegrator())):
        tran = Transient(rect(), toolkit=numeric, pcnr=True, reltol=1e-5,
                         uic=True, timestep_max=2e-5, integrator=integ)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=2e-3, timestep=2e-5, coupled_lte=True)
        t = np.asarray(res.sweep_values, float)
        assert t[-1] == pytest.approx(2e-3, rel=1e-9), name
        counts[name] = len(t)
    ## The order-2 payoff the flip was made for: measured 259 vs 769.
    assert counts['gear2'] < 0.6 * counts['euler'], counts

    ## The mask itself, structurally: source current and source node are
    ## algebraic on the rectifier; the capacitor node is a state.
    tran = Transient(rect(), toolkit=numeric)
    mask = tran._state_row_mask(np.zeros(rect().n))
    c = rect()
    assert bool(mask[c.get_node_index('b')])
    assert not bool(mask[c.get_node_index('a')])


def test_max_dv_step_voltage_check_on_algebraic_networks():
    """The Spectre/Mica-style voltage check (owner request, follow-on to
    P22): on a purely resistive/algebraic amplifier network -- the topology-
    exploration scenario, Rs + VCCS driven by a sine, no reactances -- NO
    error estimator has anything to measure, and the default run samples at
    the step cap: measured max per-step output change 2.256 V on a 10 V
    swing.  max_dv_step bounds the per-step node-voltage change instead,
    with a proportional retry; it is h-proportional by construction (|dv| ~
    h*slew), so it cannot h-cancel like the solution-LTE did.  Measured at
    landing: 0.2 V bound held at 0.199 on every path -- CPU 413 points, JAX
    412, JAX-coupled 777, from a blind 61.  Off by default (None): default
    runs are untouched, which the suite as a whole pins."""
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.elements import VCCS, VSin

    def amp():
        c = SubCircuit()
        c['vs'] = VSin('in', gnd, va=1.0, freq=1e6)
        c['Ri'] = R('in', 'x', r=1e3)
        c['Rg'] = R('x', gnd, r=9e3)
        c['gm'] = VCCS('x', gnd, 'out', gnd, gm=5e-3)
        c['RL'] = R('out', gnd, r=2e3)
        return c

    def max_step_dv(res):
        vo = np.asarray(res.v('out'), float).reshape(-1)
        return float(np.max(np.abs(np.diff(vo))))

    ## CPU: blind by default, bounded with the check.
    tran = Transient(amp(), toolkit=numeric, reltol=1e-4, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r0 = tran.solve(gnd, tend=2e-6, timestep=2e-8)
    assert max_step_dv(r0) > 1.0, 'premise gone: the default is not blind'

    ## FACTOR semantics: 2e11 * lte_vabstol(1e-12) = the 0.2 V bound.
    tran = Transient(amp(), toolkit=numeric, reltol=1e-4, uic=True,
                     max_dv_step=2e11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r1 = tran.solve(gnd, tend=2e-6, timestep=2e-8)
    assert max_step_dv(r1) <= 0.2 * (1.0 + 1e-9)
    assert len(np.asarray(r1.sweep_values)) > 5 * len(np.asarray(r0.sweep_values))

    ## The COUPLED veto is a genuinely separate code path only on the CPU
    ## (its retry lives in the Python retry loop); the JAX coupled accept
    ## shares the traced dv_ok expression with the standard path, and its
    ## chunk costs ~110 s of XLA compile per bound constant -- so coupled
    ## coverage runs here, on the CPU, in milliseconds.
    tran = Transient(amp(), toolkit=numeric, reltol=1e-4, uic=True,
                     max_dv_step=2e11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rc2 = tran.solve(gnd, tend=2e-6, timestep=2e-8, coupled_lte=True)
    assert max_step_dv(rc2) <= 0.2 * (1.0 + 1e-9)

    def go():
        tran = JAXTransient(amp(), reltol=1e-4, max_dv_step=2e11)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=2e-6, timestep=2e-8, uic=True)
        assert max_step_dv(res) <= 0.2 * (1.0 + 1e-9)
        t = np.asarray(res.sweep_values, float)
        assert t[-1] >= 2e-6 * (1.0 - 1e-9)
    _with_jax(go)

    ## FACTOR semantics with the clamp-at-1 floor (owner decision: the
    ## bound is factor*abstol, coupled to the tolerance family exactly as
    ## the LTE is, and a factor below 1 clamps to the solver-noise floor).
    tran = Transient(amp(), toolkit=numeric, max_dv_step=0.01,
                     max_di_step=0.5)
    (bv, cv), (bi, ci) = tran._dv_step_bounds()
    assert bv == float(tran.par.lte_vabstol) and cv == 0.0
    assert bi == float(tran.par.lte_iabstol) and ci == 0.0
    tran = Transient(amp(), toolkit=numeric, max_dv_step=2e11,
                     max_di_step=1e6, lte_vabstol=2e-12)
    (bv, _), (bi, _) = tran._dv_step_bounds()
    assert bv == pytest.approx(0.4)          # scales WITH lte_vabstol
    assert bi == pytest.approx(1e-6)

    ## 'auto' (owner request, sampling theory): static anchor 2*pi/N times
    ## the declared source swing, growing with the running unit-group max
    ## -- measured 0.0982 V anchor (va=1, N=64), per-step max 0.847 <=
    ## 2*pi/64 * 10 V once the output reveals the gain, 130/129 points
    ## CPU/JAX.
    import math
    tran = Transient(amp(), toolkit=numeric, reltol=1e-4, uic=True,
                     max_dv_step='auto')
    (bvs, cvr), (_bis, _cir) = tran._dv_step_bounds()
    assert cvr == pytest.approx(2 * math.pi / 64)
    assert bvs == pytest.approx(2 * math.pi / 64 * 1.0)   # va = 1 V
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ra = tran.solve(gnd, tend=2e-6, timestep=2e-8)
    assert max_step_dv(ra) <= 2 * math.pi / 64 * 10.0 * (1.0 + 1e-9)
    assert len(np.asarray(ra.sweep_values)) > 100

    def go2():
        tran = JAXTransient(amp(), max_dv_step=0.01, max_di_step=0.5)
        (bv, _), (bi, _) = tran._dv_step_bounds()
        assert bv == float(tran.par.lte_vabstol)
        assert bi == float(tran.par.lte_iabstol)
        import math
        tran = JAXTransient(amp(), reltol=1e-4, max_dv_step='auto')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ra = tran.solve(gnd, tend=2e-6, timestep=2e-8, uic=True)
        assert max_step_dv(ra) <= 2 * math.pi / 64 * 10.0 * (1.0 + 1e-9)

        ## max_di_step is live: the circuit's one branch unknown is the
        ## source current (0.1 mA peak -- the 1k/9k divider leaves 0.1*vin
        ## across Ri; VCCS is G-only, no branch row), whose per-step change
        ## at the dv-controlled resolution is ~3 uA -- a 1 uA bound
        ## (factor 1e6 * lte_iabstol) must
        ## force a visibly denser run than the voltage check alone.
        import warnings as _w
        def run(**kw):
            tran = JAXTransient(amp(), reltol=1e-4, **kw)
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                res = tran.solve(gnd, tend=2e-6, timestep=2e-8, uic=True)
            return len(np.asarray(res.sweep_values))
        n_v = run(max_dv_step=2e11)
        n_vi = run(max_dv_step=2e11, max_di_step=1e6)
        assert n_vi > 1.5 * n_v, (n_v, n_vi)
    _with_jax(go2)


def test_p24_stamps_build_identically_on_every_toolkit():
    """P24 (owner request): the stamp-construction hygiene sweep.

    16 zeros-then-mutate sites across 12 element classes were migrated to
    Toolkit.matrix_from_entries by a mechanical transformer.  The roadmap's
    open question resolved decisively first: NONE of these classes had ever
    been constructed under the JAX toolkit -- every one crashed on
    immutable-array assignment at the probe (the fail-first record).  The
    migration gates: numeric stamps stayed bit-identical (asserted 16/16 at
    landing against pre-sweep captures), and here, enduringly: each class
    builds the SAME stamp under the numeric and JAX toolkits, and
    constructs under the symbolic one."""
    from pycircuit.circuit import elements as E
    from pycircuit.circuit.toolkit import jaxtoolkit

    probes = [
        ('VCVS', lambda: E.VCVS(1, gnd, 2, gnd, g=2.0)),
        ('CCVS', lambda: E.CCVS(1, gnd, 2, gnd, r=10.0)),
        ('Nullor', lambda: E.Nullor(1, gnd, 2, gnd)),
        ('Transformer', lambda: E.Transformer(1, gnd, 2, gnd, n=2.0)),
        ('Gyrator', lambda: E.Gyrator(1, gnd, 2, gnd, gm=1e-3)),
        ('CCCS', lambda: E.CCCS(1, gnd, 2, gnd, F=2.0)),
        ('Idt', lambda: E.Idt(1, gnd, 2, gnd)),
        ('Idtmod', lambda: E.Idtmod(1, gnd, 2, gnd)),
        ('CoupledInductors', lambda: E.CoupledInductors(
            1, gnd, 2, gnd, L1=1e-6, L2=2e-6, K=0.9)),
        ('VCVS_limited', lambda: E.VCVS_limited(1, gnd, 2, gnd, g=3.0,
                                                level=2.0, offset=0.1)),
        ('ISwitch', lambda: E.ISwitch(1, gnd, 2, gnd, Ron=2.0, Roff=1e5,
                                      Ion=1e-3)),
        ('BSource', lambda: E.BSource(1, gnd, 2, gnd,
                                      i_func=lambda v: 1e-3 * v ** 3)),
    ]

    def stamps(toolkit_):
        saved = circuit_mod.default_toolkit
        circuit_mod.default_toolkit = toolkit_
        try:
            out = {}
            for name, factory in probes:
                c = SubCircuit()
                c.add_node('1'); c.add_node('2')
                c['el'] = factory()
                el = c['el']
                x = np.linspace(0.1, 0.5, el.n)
                out[name] = np.asarray(el.G(
                    x if toolkit_ is numeric else
                    __import__('jax.numpy', fromlist=['jnp']).asarray(x)),
                    dtype=float)
            return out
        finally:
            circuit_mod.default_toolkit = saved

    ref = stamps(numeric)
    jx = stamps(jaxtoolkit)
    for name in ref:
        if name in ('VCVS_limited', 'BSource'):
            ## Owner decision (option 2): compared separately at fp
            ## tolerance -- these two compute G by DIFFERENT methods per
            ## toolkit (VCVS_limited: autodiff of eval_i_pure vs the
            ## hand-derived formula; BSource: autodiff vs the numeric
            ## toolkit's finite-difference `derivative`).  Same
            ## derivative, different rounding.  Everything else is the
            ## SAME construction on both toolkits and must be bit-exact.
            np.testing.assert_allclose(jx[name], ref[name], rtol=1e-9,
                                       atol=1e-15, err_msg=name)
        else:
            np.testing.assert_array_equal(jx[name], ref[name],
                                          err_msg=name)

    ## Symbolic: constructs with exact symbols intact (the failure mode the
    ## first VCCS fix hit).
    from pycircuit.circuit.toolkit import SymbolicToolkit
    import pycircuit.circuit.toolkit as tk
    sym = getattr(tk, 'symbolic', None)
    if sym is not None:
        saved = circuit_mod.default_toolkit
        circuit_mod.default_toolkit = sym
        try:
            c = SubCircuit()
            c.add_node('1'); c.add_node('2')
            import sympy
            c['e'] = E.VCVS(1, gnd, 2, gnd, g=sympy.Symbol('A'))
            assert sympy.Symbol('A') in sympy.Matrix(np.asarray(
                c['e']._G, dtype=object)).free_symbols
        finally:
            circuit_mod.default_toolkit = saved


def test_p21_batched_dc_operating_point():
    """P21: solve_batched starts from bias -- the roadmap's last refusal
    retired.  Each lane's DC is solved with ITS OWN swept parameters (the
    point of the feature: an R sweep moves the bias per lane).  Measured at
    landing: a 2 V divider with R1 swept {1k, 3k} against Rg=1k starts its
    lanes at exactly 1.0 V and 0.5 V (analytic) and holds them flat.  The
    honest failure contract is gated too: a junction slammed at DC has no
    limiting/PCNR on this path, so the per-lane Newton reports the failing
    lanes and raises with the alternatives."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.elements import VS, Diode

    def go():
        import jax.numpy as jnp
        c = SubCircuit()
        c['V1'] = VS('in', gnd, v=2.0)
        c['R1'] = R('in', 'out', r=1e3)
        c['Rg'] = R('out', gnd, r=1e3)
        c['C1'] = C('out', gnd, c=1e-9)
        r_lanes = jnp.asarray([[1e3, 1e3], [3e3, 1e3]])  # (lane, instance)
        tran = JAXTransient(c, reltol=1e-5)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve_batched(
                gnd, override_params_tree={'R': {'r': r_lanes}},
                tend=5e-6, timestep=1e-7)          # uic defaults False now
        for lane, r1 in ((0, 1e3), (1, 3e3)):
            v = np.asarray(res[lane].v('out'), float).reshape(-1)
            exact = 2.0 * 1e3 / (r1 + 1e3)
            assert v[0] == pytest.approx(exact, rel=1e-6), (lane, v[0])
            assert v[-1] == pytest.approx(exact, rel=1e-4), (lane, v[-1])

        ## The failure contract: 5 V across a junction at DC, from zeros,
        ## with no limiting -- must raise naming the lanes, not hang or
        ## seed from a non-solution.
        c2 = SubCircuit()
        c2['vs'] = VS('a', gnd, v=5.0)
        c2['D'] = Diode('a', 'b')
        c2['R'] = R('b', gnd, r=1e3)
        c2['C'] = C('b', gnd, c=1e-9)
        tran2 = JAXTransient(c2, reltol=1e-5)
        with pytest.raises(NoConvergenceError, match='lane'):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tran2.solve_batched(
                    gnd,
                    override_params_tree={'R': {'r': jnp.asarray([[1e3], [2e3]])}},
                    tend=1e-6, timestep=1e-7)
    _with_jax(go)
