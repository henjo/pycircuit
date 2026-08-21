"""Cross-backend conformance harness (doc/transient_review_260820.md, Phase 1).

The divergence-class counterpart of the dead-knob scan: the same circuits run
through `Transient` and `JAXTransient` at matched reltol AND matched
integration order -- the CPU side pins Gear2Integrator, because the shipped
defaults differ (CPU order 1, JAX order 2; F19(b)) and a comparison across
orders measures the method gap, not backend drift.

Contract points asserted on every pair, plus a waveform-agreement bound.
Bounds are ~5x the deviation measured when the harness landed (rc 7.4e-4,
vsin 3.5e-3 at reltol 1e-4), so genuine drift trips them while Phase-3
behavior changes have headroom; tighten at the Phase-3 exit gate if the
measured deltas allow.  This harness is that gate's acceptance instrument.
"""

import warnings

import numpy as np
import pytest

jax = pytest.importorskip('jax')

from pycircuit.circuit import circuit as circuit_mod, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS, VSin
from pycircuit.circuit.integrator import Gear2Integrator
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


def _cpu(kind, uic):
    tran = Transient(_build(kind), toolkit=numeric,
                     integrator=Gear2Integrator(), reltol=RELTOL, uic=uic)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=TEND, timestep=TIMESTEP)
    return (np.asarray(res.sweep_values, float),
            np.asarray(res.v('out'), float).reshape(-1))


def _jax(kind, uic):
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        from pycircuit.circuit.jaxtransient import JAXTransient
        tran = JAXTransient(_build(kind), reltol=RELTOL)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=TEND, timestep=TIMESTEP, uic=uic)
        return (np.asarray(res.sweep_values, float),
                np.asarray(res.v('out'), float).reshape(-1))
    finally:
        circuit_mod.default_toolkit = saved


@pytest.mark.parametrize('kind,uic,bound', [('rc', True, 5e-3),
                                            ('vsin', False, 2e-2)])
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
