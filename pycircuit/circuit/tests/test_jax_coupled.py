"""Fang's coupled (x, h) stepping on the JAX backend (parity item P19).

The 'approx' branch of the CPU's fang_timestep, traced.  Measured at landing:
rc step counts 127/127 CPU/JAX, cross-backend deviation 6.8e-3 (rc) and
4.3e-3 (vsin) at reltol 1e-4, JAX-coupled MORE accurate than CPU-coupled at
every tolerance probed (1.0e-3/1.0e-3/2.9e-4 vs 1.45e-2/5.8e-3/1.9e-3 at
reltol 1e-3/1e-4/1e-5 on the uic RC), and per-lane batched coupled errors of
1.3e-3 / 1.1e-5 against analytic references.

The vsin STEP-COUNT comparison is deliberately absent: CPU takes ~3x the
steps there because Sin.next_event fires quarter-period breakpoints on the
CPU and returns inf on this backend -- a breakpoint-feature divergence
(parity doc P14 note), not a property of the coupled method.
"""

import warnings

import numpy as np
import pytest

jax = pytest.importorskip('jax')
import jax.numpy as jnp

from pycircuit.circuit import circuit as circuit_mod, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS, VSin


def _build(kind='rc'):
    c = SubCircuit()
    c.add_node('in'); c.add_node('out')
    c['V1'] = VS('in', gnd, v=1.0) if kind == 'rc' else \
        VSin('in', gnd, va=1.0, freq=1e3)
    c['R1'] = R('in', 'out', r=1e3)
    c['C1'] = C('out', gnd, c=1e-6)
    return c


def _jax_coupled(kind, uic, reltol):
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        tran = JAXTransient(_build(kind), reltol=reltol, coupled_lte=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=5e-3, timestep=1e-4, uic=uic)
        return (np.asarray(res.sweep_values, float),
                np.asarray(res.v('out'), float).reshape(-1))
    finally:
        circuit_mod.default_toolkit = saved


def _cpu_coupled(kind, uic, reltol):
    from pycircuit.circuit.transient import Transient
    tran = Transient(_build(kind), toolkit=numeric, reltol=reltol, uic=uic)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=5e-3, timestep=1e-4, coupled_lte=True)
    return (np.asarray(res.sweep_values, float),
            np.asarray(res.v('out'), float).reshape(-1))


@pytest.mark.parametrize('kind,uic,bound', [('rc', True, 2e-2),
                                            ('vsin', False, 2e-2)])
def test_jax_coupled_matches_cpu_coupled(kind, uic, bound):
    tc, vc = _cpu_coupled(kind, uic, 1e-4)
    tj, vj = _jax_coupled(kind, uic, 1e-4)
    for t in (tc, tj):
        assert t[0] == 0.0
        assert t[-1] == pytest.approx(5e-3, rel=1e-12)   # F3's guarantee holds
        assert np.all(np.diff(t) > 0)
    dev = float(np.max(np.abs(vc - np.interp(tc, tj, vj))))
    assert dev < bound, 'coupled backends diverge: %.3e' % dev
    if kind == 'rc':
        ## no breakpoint divergence on this drive: step parity is tight
        ## (127/127 measured at landing)
        assert 0.8 < len(tc) / len(tj) < 1.25


def test_jax_coupled_reltol_governs():
    errs = []
    for rt in (1e-3, 1e-5):
        tj, vj = _jax_coupled('rc', True, rt)
        errs.append(float(np.max(np.abs(vj - (1 - np.exp(-tj / 1e-3))))))
    assert errs[1] < errs[0], \
        'two decades of reltol did not move the error: %r' % errs


def test_jax_coupled_batched_lanes():
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        tran = JAXTransient(_build('rc'), reltol=1e-4, coupled_lte=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve_batched(
                gnd, override_params_tree={'R': {'r': jnp.array([[1e2], [1e4]])}},
                tend=5e-3, timestep=1e-4, uic=True)
        for b, tau, bound in ((0, 1e-4, 5e-3), (1, 1e-2, 1e-3)):
            t = np.asarray(res[b].sweep_values, float)
            v = np.asarray(res[b].v('out'), float).reshape(-1)
            assert t[0] == 0.0 and np.all(np.diff(t) > 0)
            assert t[-1] == pytest.approx(5e-3, rel=1e-12)
            err = float(np.max(np.abs(v - (1 - np.exp(-t / tau)))))
            assert err < bound, 'lane %d err %.2e' % (b, err)
    finally:
        circuit_mod.default_toolkit = saved


def test_jax_coupled_band_sentinel():
    """The 'auto' band resolves to Fang's values; explicit values survive --
    the F5 contract, on this backend."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        assert JAXTransient(_build())._coupled_band() == (0.7, 3.0, 0.15)
        assert JAXTransient(_build(), lte_gamma_min=0.0, lte_gamma_max=1.0,
                            lte_eta=None)._coupled_band() == (0.0, 1.0, None)
    finally:
        circuit_mod.default_toolkit = saved


def test_jax_coupled_tline_matches_cpu_standard():
    """RE-DERIVED from a refusal test (this used to expect
    NotImplementedError).  The refusal's cause was a livelock, and the fix
    took three composed pieces, each earned by a traced measurement:

    - a from-zero kink makes every signal scale together, so the coupled
      LTE's dev is proportional to its ref and h CANCELS -- err =
      1/(TRTOL*reltol) = 1428.6 measured, constant from h=4.4e-11 down to
      the excursion floor.  No step size can satisfy the band; the outer
      retry collapses dt to dt_min.
    - one graced step is not enough: the NEXT step extrapolates degree-2
      through a PRE-kink history point, and a quadratic through two flat
      points and one ramp point misses a line by err = 139.2e11*h --
      in-band only below the within-point floor.  Hence the coupled accept
      path ZEROES h_history on a breakpoint landing (band skipped one step,
      degree held at 1 for two, cold-start semantics).
    - the far-end wavefront onset is the same kink NOT at any element-
      reported breakpoint, so collect_breakpoints registers source corners
      + k*TD arrivals.  (Arrivals alone were tried first and falsified --
      a registered kink without the history reset still livelocks.)

    The gate: the pulsed matched line that livelocked at every tend past
    the first arrival now lands on tend and matches the CPU standard path
    to 5e-16 (measured 4.996e-16 at tend=8e-9, 112 points)."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.integrator import Gear2Integrator
    from pycircuit.circuit.elements import TLine, VPulse

    def line():
        c = SubCircuit()
        c.add_node('a'); c.add_node('b')
        c['V1'] = VPulse('s', gnd, v1=0.0, v2=1.0, td=1e-9, tr=2e-10,
                         tf=2e-10, pw=1e-8, per=1e-7)
        c['Rs'] = R('s', 'a', r=50.0)
        c['T1'] = TLine('a', gnd, 'b', gnd, Z0=50.0, TD=1e-9)
        c['Rl'] = R('b', gnd, r=50.0)
        return c

    cpu = Transient(line(), toolkit=numeric, reltol=1e-4,
                    integrator=Gear2Integrator(), uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res_c = cpu.solve(gnd, tend=8e-9, timestep=2e-10)
    tc = np.asarray(res_c.sweep_values, float)
    vc = np.asarray(res_c.v('b'), float).reshape(-1)

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        tran = JAXTransient(line(), reltol=1e-4, coupled_lte=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=8e-9, timestep=2e-10, uic=True)
        tj = np.asarray(res.sweep_values, float)
        vj = np.asarray(res.v('b'), float).reshape(-1)
    finally:
        circuit_mod.default_toolkit = saved

    ## The run must reach tend (the livelock died mid-ramp at ~1.05e-9)...
    assert tj[-1] > 8e-9 * (1.0 - 1e-6)
    ## ...must carry the delayed wavefront (the DC-stamp defect gave 24.5 V
    ## garbage; a dropped history gives vb == 0 on a matched line)...
    vb_delay = tj[np.argmax(vj > 0.25)]
    assert 1.0e-9 < vb_delay < 2.5e-9
    ## ...and must match the CPU standard path.  5e-13 leaves three orders
    ## over the measured 5e-16 without letting a controller regression hide.
    dev = float(np.max(np.abs(vc - np.interp(tc, tj, vj))))
    assert dev < 5e-13, 'coupled+TLine drifted from CPU: %.3e' % dev
