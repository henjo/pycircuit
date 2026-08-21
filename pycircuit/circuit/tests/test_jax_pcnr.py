"""PCNR on the traced loop (parity item P19, stage 2).

On this backend PCNR is not an alternative to limiting -- it is the ONLY
junction-robustness mechanism possible: Diode.limit keeps per-device Python
state (_vlim), which a lax.while_loop cannot host, while PCNR's limited
quantities are a per-junction state vector.

Measured at landing: the cold-start probe below FAILS outright on the plain
per-row Newton (NoConvergenceError with dt driven to dt_min) and completes in
35 steps under PCNR; on a converging rectifier the two paths agree to 1e-15
(both solve the same equations to the same tolerance -- identical results on
easy circuits are PCNR working, not PCNR inert) and JAX-PCNR tracks CPU-PCNR
to 1.1e-2; wall-clock cost +29% on the rectifier, against the CPU's measured
+60-80% per iteration -- the CPU figure did not transfer, as the parity doc
predicted it must be re-measured.
"""

import warnings

import numpy as np
import pytest

jax = pytest.importorskip('jax')

from pycircuit.circuit import circuit as circuit_mod, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS, VSin, Diode


def _with_jax(fn):
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        return fn()
    finally:
        circuit_mod.default_toolkit = saved


def _cold_start():
    c = SubCircuit()
    c['vs'] = VS('a', gnd, v=5.0)
    c['D'] = Diode('a', 'b')
    c['R'] = R('b', gnd, r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    return c


def _rectifier():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=5.0, freq=1e3)
    c['D'] = Diode('a', 'b')
    c['R'] = R('b', gnd, r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c


def test_pcnr_solves_the_cold_start_plain_newton_cannot():
    """The value demonstration: 5 V slammed across a junction in one
    full-size step (firststep=timestep kills the opening ramp)."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.nrsolver import NoConvergenceError

    def run(pcnr):
        tran = JAXTransient(_cold_start(), reltol=1e-5, firststep=1e-6,
                            pcnr=pcnr)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=2e-5, timestep=1e-6, uic=True)
        return np.asarray(res.v('b'), float).reshape(-1)

    with pytest.raises(NoConvergenceError):
        _with_jax(lambda: run(False))

    v = _with_jax(lambda: run(True))
    assert np.all(np.isfinite(v))
    ## settled at source minus one diode drop
    assert v[-1] == pytest.approx(4.367, abs=0.01)


def test_pcnr_matches_cpu_pcnr_on_the_rectifier():
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.transient import Transient

    tran_c = Transient(_rectifier(), toolkit=numeric, pcnr=True,
                       reltol=1e-5, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res_c = tran_c.solve(tend=2e-3, timestep=2e-5)
    tc = np.asarray(res_c.sweep_values, float)
    vc = np.asarray(res_c.v('b'), float).reshape(-1)

    def run():
        tran = JAXTransient(_rectifier(), reltol=1e-5, pcnr=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=2e-3, timestep=2e-5, uic=True)
        assert tran.statistics.nonconverged_steps == 0
        return (np.asarray(res.sweep_values, float),
                np.asarray(res.v('b'), float).reshape(-1))

    tj, vj = _with_jax(run)
    dev = float(np.max(np.abs(vc - np.interp(tc, tj, vj))))
    assert dev < 3e-2, 'JAX-PCNR drifted from CPU-PCNR: %.3e' % dev


def test_pcnr_without_junctions_falls_through():
    """A circuit with no participating device runs the plain Newton --
    bit-identically, as the CPU's solve_timestep dispatch does."""
    from pycircuit.circuit.jaxtransient import JAXTransient

    def run(pcnr):
        c = SubCircuit()
        c.add_node('in'); c.add_node('out')
        c['V1'] = VS('in', gnd, v=1.0)
        c['R1'] = R('in', 'out', r=1e3)
        c['C1'] = C('out', gnd, c=1e-6)
        tran = JAXTransient(c, pcnr=pcnr)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=1e-4, timestep=1e-5, uic=True)
        return np.asarray(res.x)

    a = _with_jax(lambda: run(False))
    b = _with_jax(lambda: run(True))
    assert np.array_equal(a, b)


def test_pcnr_inside_coupled_is_refused():
    from pycircuit.circuit.jaxtransient import JAXTransient

    def run():
        tran = JAXTransient(_rectifier(), pcnr=True, coupled_lte=True)
        with pytest.raises(NotImplementedError, match='coupled'):
            tran.solve(gnd, tend=1e-5, timestep=1e-6, uic=True)

    _with_jax(run)


def test_pcnr_refuses_charge_storing_junctions():
    """Mirror of the CPU augmented_system refusal."""
    from pycircuit.circuit.jaxtransient import JAXTransient

    class ChargedDiode(Diode):
        @staticmethod
        def eval_q_pure(x, params, epar, toolkit):
            q = 1e-12 * (x[0] - x[1])
            return toolkit.array([q, -q])

    def run():
        c = SubCircuit()
        c['vs'] = VS('a', gnd, v=1.0)
        c['D'] = ChargedDiode('a', 'b')
        c['R'] = R('b', gnd, r=1e3)
        tran = JAXTransient(c, pcnr=True)
        with pytest.raises(NotImplementedError, match='charge'):
            tran.solve(gnd, tend=1e-6, timestep=1e-7, uic=True)

    _with_jax(run)
