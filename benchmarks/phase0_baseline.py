"""Phase-0 performance baseline (doc/transient_review_260820.md, Phase 0 item 6).

Records wall-clock and step counts for representative transients on both
backends, so the Phase-3 JAX-behavior fixes are judged against a number rather
than an impression.  Bounded on purpose: two circuits, three configurations,
single repetition after a warm-up -- a yardstick, not a benchmark suite.
"""
import time
import warnings
import numpy as np


def rc_vsin_cpu(integrator=None):
    from pycircuit.circuit import gnd, numeric
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import R, C, VSin
    from pycircuit.circuit.transient import Transient
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    kw = {} if integrator is None else {'integrator': integrator}
    return Transient(c, toolkit=numeric, **kw)


def run_cpu(label, integrator=None):
    tran = rc_vsin_cpu(integrator)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tran.solve(tend=2e-3, timestep=1e-5)        # warm-up
        t0 = time.perf_counter()
        res = tran.solve(tend=2e-3, timestep=1e-5)
        dt = time.perf_counter() - t0
    st = res.statistics
    print(f"{label:34s} wall={dt*1e3:8.1f} ms  accepted={st.accepted_steps:5d} "
          f"rejected={st.rejected_steps:4d}")


def run_jax(label):
    from pycircuit.circuit import circuit as circuit_mod, gnd
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        from pycircuit.circuit.circuit import SubCircuit
        from pycircuit.circuit.elements import R, C, VSin
        from pycircuit.circuit.jaxtransient import JAXTransient
        c = SubCircuit()
        c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
        c['R'] = R('a', 'b', r=1e3)
        c['C'] = C('b', gnd, c=1e-7)
        tran = JAXTransient(c)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tran.solve(gnd, tend=2e-3, timestep=1e-5)   # warm-up incl. jit
            t0 = time.perf_counter()
            tran.solve(gnd, tend=2e-3, timestep=1e-5)
            dt = time.perf_counter() - t0
        st = tran.statistics
        print(f"{label:34s} wall={dt*1e3:8.1f} ms  accepted={st.accepted_steps:5d} "
              f"rejected={st.rejected_steps:4d}")
    finally:
        circuit_mod.default_toolkit = saved


if __name__ == '__main__':
    from pycircuit.circuit.integrator import Gear2Integrator
    run_cpu('CPU rc-vsin Euler (default)')
    run_cpu('CPU rc-vsin Gear2', Gear2Integrator())
    run_jax('JAX rc-vsin gear (default)')
