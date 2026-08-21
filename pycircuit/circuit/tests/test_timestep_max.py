"""timestep_max -- the step cap, decoupled from `timestep` (owner decision).

`timestep` used to double as the largest accepted step on both backends.
That silently made the step count on gentle circuits a property of the
requested output density rather than of the error control: measured before
the change, the JAX rc-vsin run took IDENTICAL 209-step runs at reltol 1e-4
and 1e-6 -- no tolerance knob could move it.  The cap is now its own
Parameter, `timestep_max`, defaulting to tend/50 (SPICE's TMAX); `timestep`
only sets the opening-step scale and the fixed_timestep grid.
"""

import warnings

import numpy as np
import pytest

jax = pytest.importorskip('jax')

from pycircuit.circuit import circuit as circuit_mod, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VSin


def _with_jax(fn):
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        return fn()
    finally:
        circuit_mod.default_toolkit = saved


def _rc_vsin():
    c = SubCircuit()
    c['V1'] = VSin('in', gnd, va=1.0, freq=1e3)
    c['R1'] = R('in', 'out', r=1e3)
    c['C1'] = C('out', gnd, c=1e-6)
    return c


def _max_h(res):
    return float(np.max(np.diff(np.asarray(res.sweep_values, float))))


def test_the_old_name_is_gone_on_both_backends():
    from pycircuit.circuit.transient import Transient
    with pytest.raises(KeyError):
        Transient(_rc_vsin(), toolkit=numeric, max_step=1e-4)

    def go():
        from pycircuit.circuit.jaxtransient import JAXTransient
        with pytest.raises(KeyError):
            JAXTransient(_rc_vsin(), max_step=1e-4)
    _with_jax(go)


def test_steps_grow_past_timestep_up_to_tend_over_50():
    """The decoupling itself: with timestep=1e-5 and tend=2e-3, both
    backends take steps larger than timestep, bounded by tend/50 = 4e-5.
    Measured at landing: CPU max_h 3.88e-5, JAX max_h 4.00e-5."""
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.jaxtransient import JAXTransient

    tran = Transient(_rc_vsin(), toolkit=numeric, reltol=1e-4, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(gnd, tend=2e-3, timestep=1e-5)
    assert _max_h(res) > 1e-5
    assert _max_h(res) <= 4e-5 * (1.0 + 1e-9)

    def go():
        tran = JAXTransient(_rc_vsin(), reltol=1e-4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=2e-3, timestep=1e-5, uic=True)
        assert _max_h(res) > 1e-5
        assert _max_h(res) <= 4e-5 * (1.0 + 1e-9)
    _with_jax(go)


def test_an_explicit_cap_binds_exactly_on_both_backends():
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.jaxtransient import JAXTransient

    tran = Transient(_rc_vsin(), toolkit=numeric, reltol=1e-4, uic=True,
                     timestep_max=2e-5)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(gnd, tend=2e-3, timestep=1e-5)
    assert _max_h(res) <= 2e-5 * (1.0 + 1e-9)

    def go():
        tran = JAXTransient(_rc_vsin(), reltol=1e-4, timestep_max=2e-5)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=2e-3, timestep=1e-5, uic=True)
        assert _max_h(res) <= 2e-5 * (1.0 + 1e-9)
    _with_jax(go)


def test_tolerance_knobs_are_live_on_jax_again():
    """THE POINT of the decision: this exact configuration measured
    identical 209-step runs at reltol 1e-4 and 1e-6 while the cap bound
    every step.  Decoupled, the error control takes over: measured 63
    steps at 1e-4 against 216 at 1e-6 at landing."""
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        counts = {}
        for rel in (1e-4, 1e-6):
            tran = JAXTransient(_rc_vsin(), reltol=rel)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = tran.solve(gnd, tend=2e-3, timestep=1e-5, uic=True)
            counts[rel] = len(np.asarray(res.sweep_values))
        assert counts[1e-6] > 2 * counts[1e-4], counts
    _with_jax(go)
