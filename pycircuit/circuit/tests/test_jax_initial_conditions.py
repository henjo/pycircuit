"""P12: real initial conditions on the JAX backend, via the CPU's machinery.

The `uic=True` path used to walk a `node.ic` attribute that nothing in the
package or tests ever set -- dead code posing as a feature.  The CPU's
`_initial_state` chain (ic dict -> element ICs -> spanning-tree capacitor
solve) is pre-loop Python over names and indices, so `JAXTransient` now binds
those methods unchanged; these tests run the CPU suite's key cases through
the JAX class, spanning-tree cases included.
"""

import warnings

import numpy as np
import pytest

jax = pytest.importorskip('jax')

from pycircuit.circuit import circuit as circuit_mod, gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, L, VS


def _with_jax(fn):
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        return fn()
    finally:
        circuit_mod.default_toolkit = saved


def test_node_ic_reaches_the_waveform():
    """End-to-end: an RC released from ic={'out': 0.5} decays from 0.5."""
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        c = SubCircuit()
        c['R'] = R('out', gnd, r=1e3)
        c['C'] = C('out', gnd, c=1e-6)
        tran = JAXTransient(c, reltol=1e-4, ic={'out': 0.5})
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=2e-3, timestep=1e-5, uic=True)
        t = np.asarray(res.sweep_values, float)
        v = np.asarray(res.v('out'), float).reshape(-1)
        assert v[0] == pytest.approx(0.5)
        exact = 0.5 * np.exp(-t / 1e-3)
        assert np.max(np.abs(v - exact)) < 5e-3
    _with_jax(go)


def test_inductor_ic_lands_on_its_own_branch_row():
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        c = SubCircuit()
        c['L'] = L(1, gnd, L=1e-6, ic=0.5)
        c['C'] = C(1, gnd, c=1e-9)
        tran = JAXTransient(c, uic=True)
        x0 = tran._initial_state(gnd)
        row = c.instance_branch_indices('L')[0]
        assert x0[row] == pytest.approx(0.5)
        assert np.count_nonzero(np.asarray(x0)) == 1
    _with_jax(go)


def test_capacitor_ics_solve_as_a_spanning_tree():
    """The chain vC1 + vC2 stacked from ground: node voltages come out as
    cumulative sums, exactly the CPU case."""
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        c = SubCircuit()
        c['C1'] = C(1, gnd, c=1e-9, ic=1.5)
        c['C2'] = C(2, 1, c=1e-9, ic=0.25)
        c['R'] = R(2, gnd, r=1e6)
        tran = JAXTransient(c, uic=True)
        x0 = tran._initial_state(gnd)
        assert x0[c.get_node_index('1')] == pytest.approx(1.5)
        assert x0[c.get_node_index('2')] == pytest.approx(1.75)
    _with_jax(go)


def test_a_floating_ic_group_raises():
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        c = SubCircuit()
        c['C1'] = C(1, 2, c=1e-9, ic=1.0)
        c['R'] = R(1, 2, r=1e6)
        ## Ground the circuit somewhere OUTSIDE the constrained pair.
        c['R2'] = R(3, gnd, r=1e3)
        c['R3'] = R(2, 3, r=1e33)
        tran = JAXTransient(c, uic=True)
        with pytest.raises(ValueError, match='up to a constant'):
            tran._initial_state(gnd)
    _with_jax(go)


def test_ic_naming_the_reference_node_raises():
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        c = SubCircuit()
        c['R'] = R(1, gnd, r=1e3)
        c['C'] = C(1, gnd, c=1e-9)
        tran = JAXTransient(c, ic={gnd: 1.0}, uic=True)
        with pytest.raises(ValueError, match='reference node'):
            tran._initial_state(gnd)
    _with_jax(go)


def test_ic_without_uic_raises():
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        c = SubCircuit()
        c['V'] = VS(1, gnd, v=1.0)
        c['R'] = R(1, 2, r=1e3)
        c['C'] = C(2, gnd, c=1e-9)
        tran = JAXTransient(c, ic={2: 0.5})
        with pytest.raises(ValueError, match='uic=True'):
            tran.solve(gnd, tend=1e-6, timestep=1e-7)
    _with_jax(go)


def test_the_dead_node_ic_walk_is_gone():
    """A `node.ic` ATTRIBUTE (the old dead feature) must have no effect --
    only the `ic` Parameter reaches the starting vector."""
    from pycircuit.circuit.jaxtransient import JAXTransient

    def go():
        c = SubCircuit()
        c['R'] = R('out', gnd, r=1e3)
        c['C'] = C('out', gnd, c=1e-6)
        for node in c.nodes:
            if str(node) == 'out':
                node.ic = 123.0   ## the old walk would have read this
        tran = JAXTransient(c, uic=True)
        x0 = np.asarray(tran._initial_state(gnd))
        assert np.count_nonzero(x0) == 0
    _with_jax(go)
