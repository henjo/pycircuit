"""Cross-backend conformance harness (doc/transient_review_260820.md, Phase 1).

The divergence-class counterpart of the dead-knob scan: the same circuits run
through `Transient` and `JAXTransient` at matched reltol AND matched
integration order -- the CPU side pins Gear2Integrator, because the shipped
defaults differ (CPU order 1, JAX order 2; F19(b)) and a comparison across
orders measures the method gap, not backend drift.

Contract points asserted on every pair, plus a waveform-agreement bound.

Bound history: at landing (pre-Phase-3) the measured deviations were
rc 7.4e-4 / vsin 3.5e-3 at reltol 1e-4, bounds set ~5x.  At the Phase-3
exit gate the cluster (F19, F11, F17, F6(b), F16) had IMPROVED them to
rc 4.3e-6 / vsin 2.9e-3 -- the per-row Newton criterion accounts for most
of the rc gain -- and step counts matched 120/120.  Bounds re-tightened
per the gate's rule to ~5x the new measurements.
"""

import warnings

import numpy as np
import pytest

jax = pytest.importorskip('jax')

from pycircuit.circuit import circuit as circuit_mod, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS, VSin
from pycircuit.circuit.integrator import Gear2Integrator, TrapezoidalIntegrator
from pycircuit.circuit.toolkit import jaxtoolkit
from pycircuit.circuit.transient import Transient

TEND, TIMESTEP, RELTOL = 5e-3, 1e-4, 1e-4


def _build(kind):
    c = SubCircuit()
    c.add_node('in'); c.add_node('out')
    if kind == 'rc':
        c['V1'] = VS('in', gnd, v=1.0)
    else:
        c['V1'] = VSin('in', gnd, va=1.0, freq=1e3)
    c['R1'] = R('in', 'out', r=1e3)
    c['C1'] = C('out', gnd, c=1e-6)
    return c


def _cpu(kind, uic, integrator=None):
    tran = Transient(_build(kind), toolkit=numeric,
                     integrator=integrator or Gear2Integrator(),
                     reltol=RELTOL, uic=uic)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=TEND, timestep=TIMESTEP)
    return (np.asarray(res.sweep_values, float),
            np.asarray(res.v('out'), float).reshape(-1))


def _jax(kind, uic, integrator='gear'):
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        from pycircuit.circuit.jaxtransient import JAXTransient
        tran = JAXTransient(_build(kind), reltol=RELTOL, integrator=integrator)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=TEND, timestep=TIMESTEP, uic=uic)
        return (np.asarray(res.sweep_values, float),
                np.asarray(res.v('out'), float).reshape(-1))
    finally:
        circuit_mod.default_toolkit = saved


@pytest.mark.parametrize('kind,uic,bound', [('rc', True, 2e-5),
                                            ('vsin', False, 1.5e-2)])
def test_backends_agree_at_matched_order(kind, uic, bound):
    tc, vc = _cpu(kind, uic)
    tj, vj = _jax(kind, uic)

    ## Contract points (the F12/F3 class):
    for t in (tc, tj):
        assert t[0] == 0.0                       # t=0 included
        assert t[-1] == pytest.approx(TEND, rel=1e-12)   # lands on tend
        assert np.all(np.diff(t) > 0)            # strictly increasing

    ## Waveform agreement, JAX interpolated onto the CPU grid:
    dev = float(np.max(np.abs(vc - np.interp(tc, tj, vj))))
    assert dev < bound, 'backend deviation %.3e exceeds %.0e' % (dev, bound)

    ## Step counts within 2x of each other -- same method, same tolerance,
    ## grossly different work means one controller is not doing its job:
    assert 0.5 < len(tc) / len(tj) < 2.0


@pytest.mark.parametrize('kind,uic,bound', [('rc', True, 2e-5),
                                            ('vsin', False, 2e-2)])
def test_backends_agree_on_trapezoidal(kind, uic, bound):
    """The second integrator pairing this harness has ever carried.

    Trapezoidal was CPU-only until 2026-09-01 -- the JAX branch had been
    deleted for cause and the record said a variable-step estimator "has not
    been written", which was wrong: the kernels existed, the wiring did not.

    Bounds follow the harness rule, ~5x the measurement at landing:
    rc 4.2779e-06, vsin 4.1643e-03, step counts 61/60 and 92/89 CPU/JAX.
    For reference the Gear-2 pairing measures rc 4.3e-6 / vsin 2.9e-3, so
    trapezoidal agrees to the same order on rc and is somewhat looser on
    vsin -- expected, and NOT to be tuned away: the two backends drop order
    on different rules (JAX's `effective_first_order` opens with two Euler
    steps and re-drops after every breakpoint, where the CPU's
    `TrapezoidalIntegrator.check_order_drop` drops on the first step only).

    ⚠ Trapezoidal's ringing on stiff modes is UNMITIGATED on both backends
    (doc/transient_review.md 4.6: |e_n/e_(n-1)| ~ 0.9960 at h*lambda = -1000,
    against Gear-2's 0.0972).  This test pins AGREEMENT, not quality; a JAX
    trap run that hid the ringing would be a defect, not an improvement.
    """
    tc, vc = _cpu(kind, uic, integrator=TrapezoidalIntegrator())
    tj, vj = _jax(kind, uic, integrator='trap')

    for t in (tc, tj):
        assert t[0] == 0.0
        assert t[-1] <= TEND * (1 + 1e-9)
        assert np.all(np.diff(t) > 0)

    dev = float(np.max(np.abs(vj - np.interp(tj, tc, vc))))
    assert dev < bound, '%s: trapezoidal backends disagree by %.3e' % (kind, dev)


def test_trapezoidal_rings_on_both_backends_identically():
    """The gate that proves the JAX branch is really trapezoidal.

    A trap estimator that silently fell back to Euler would still agree with
    the CPU on smooth waveforms and would still pass the pairing above -- and
    it would be wrong.  Trapezoidal's signature is its UNDAMPED behaviour on a
    stiff mode: at `h*lambda = -1000` its amplification factor is
    `|R_TR(-1000)| = 0.996008` (doc/transient_review.md 4.6), so the error
    alternates sign and decays 0.4% per step, where Gear-2 decays ~0.07 and
    Euler ~0.001.  A fixed grid puts both backends on identical steps, so what
    is left is the method.

    Measured 2026-09-01: cpu 0.9960, jax 0.9960 -- the ringing is REPRODUCED,
    which is the correct outcome.  It is not mitigated on either backend and
    this test must not be "fixed" by damping the JAX side; that would be a
    parity break dressed as an improvement.
    """
    from pycircuit.circuit.integrator import TrapezoidalIntegrator

    dt = 1e-3                      # tau = 1e-6, so h/tau = 1000

    def stiff():
        c = SubCircuit()
        c.add_node('in'); c.add_node('out')
        c['V1'] = VS('in', gnd, v=1.0)
        c['R1'] = R('in', 'out', r=1e3)
        c['C1'] = C('out', gnd, c=1e-9)
        return c

    def decay(v):
        e = np.abs(np.asarray(v, float).reshape(-1)[3:] - 1.0)
        e = e[e > 0]
        return float(np.median(e[1:] / e[:-1]))

    tran = Transient(stiff(), toolkit=numeric,
                     integrator=TrapezoidalIntegrator(), reltol=1e-4, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rc = tran.solve(tend=40 * dt, timestep=dt, fixed_timestep=True)
    d_cpu = decay(rc.v('out'))

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        from pycircuit.circuit.jaxtransient import JAXTransient
        j = JAXTransient(stiff(), integrator='trap', reltol=1e-4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            rj = j.solve(gnd, tend=40 * dt, timestep=dt, uic=True,
                         fixed_timestep=True)
        d_jax = decay(rj.v('out'))
    finally:
        circuit_mod.default_toolkit = saved

    ## the textbook value, on both sides
    assert d_cpu == pytest.approx(0.996008, abs=5e-3), d_cpu
    assert d_jax == pytest.approx(0.996008, abs=5e-3), d_jax
    ## and, decisively, NOT an Euler fallback in disguise
    assert d_jax > 0.9, 'the JAX trap branch damped like a lower order: %.4f' \
        % d_jax
