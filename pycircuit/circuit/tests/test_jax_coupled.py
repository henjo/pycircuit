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


def test_jax_coupled_refuses_tline():
    """The coupled assembly does not apply the delay-line history; running a
    TLine circuit would silently drop the reflections -- refused instead.
    (The dead-knob scan caught the original port accepting the TLine buffers
    and reading neither: this is that finding, resolved honestly.)"""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.elements import TLine
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        c = SubCircuit()
        c.add_node('a'); c.add_node('b')
        c['V1'] = VS('a', gnd, v=1.0)
        c['T1'] = TLine('a', gnd, 'b', gnd, Z0=50.0, TD=1e-9)
        c['R1'] = R('b', gnd, r=50.0)
        tran = JAXTransient(c, coupled_lte=True)
        with pytest.raises(NotImplementedError, match='TLine'):
            tran.solve(gnd, tend=1e-8, timestep=1e-10, uic=True)
    finally:
        circuit_mod.default_toolkit = saved
