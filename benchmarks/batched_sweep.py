"""The batched-sweep payoff benchmark (doc/batched_sweep_260821.md).

The branch's reason to exist, measured: an N-lane resistor corner sweep of
the half-wave rectifier (VSin -> Diode -> R || C, r log-spaced 300 ohm ..
30 kohm) through `JAXTransient.solve_batched` -- per-lane DC bias via the
P18/P25 continuation ladders included -- against the honest CPU
alternative, a Python loop over `Transient` rebuilding the circuit per
lane.  Every number is measured, none extrapolated: the CPU loop really
runs all N lanes at every N.  Bounded on purpose, like phase0_baseline:
one circuit, one tolerance, single timed repetition after a warm-up --
a yardstick, not a benchmark suite.

Run:  XLA_PYTHON_CLIENT_PREALLOCATE=false python benchmarks/batched_sweep.py
"""
import time
import warnings
import numpy as np

TEND = 2e-3          # two periods of the 1 kHz drive
TIMESTEP = 1e-5
R_LO, R_HI = 3e2, 3e4


def lane_values(n):
    return np.logspace(np.log10(R_LO), np.log10(R_HI), n)


def build(R, C, VSin, Diode, SubCircuit, gnd, r):
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=5.0, freq=1e3)
    c['D'] = Diode('a', 'b')
    c['R'] = R('b', gnd, r=float(r))
    c['C'] = C('b', gnd, c=1e-7)
    return c


def run_cpu(n):
    """The honest baseline: one Transient per lane, circuit rebuilt per
    lane (what a user sweeping corners today actually writes)."""
    from pycircuit.circuit import gnd, numeric
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import R, C, VSin, Diode
    from pycircuit.circuit.transient import Transient
    finals = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t0 = time.perf_counter()
        for r in lane_values(n):
            c = build(R, C, VSin, Diode, SubCircuit, gnd, r)
            res = Transient(c, toolkit=numeric).solve(
                tend=TEND, timestep=TIMESTEP)
            finals.append(float(np.asarray(res.v('b'))[-1]))
        wall = time.perf_counter() - t0
    return wall, np.asarray(finals)


def run_jax(n):
    """One solve_batched call, all lanes vmapped; cold run includes the
    jit compile for this lane count, warm run is the steady state."""
    from pycircuit.circuit import circuit as circuit_mod, gnd
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        import jax.numpy as jnp
        from pycircuit.circuit.circuit import SubCircuit
        from pycircuit.circuit.elements import R, C, VSin, Diode
        from pycircuit.circuit.jaxtransient import JAXTransient
        c = build(R, C, VSin, Diode, SubCircuit, gnd, 1e3)
        tree = {'R': {'r': jnp.asarray(lane_values(n)).reshape(-1, 1)}}
        tran = JAXTransient(c)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            t0 = time.perf_counter()
            tran.solve_batched(gnd, override_params_tree=tree,
                               tend=TEND, timestep=TIMESTEP)
            cold = time.perf_counter() - t0
            t0 = time.perf_counter()
            res = tran.solve_batched(gnd, override_params_tree=tree,
                                     tend=TEND, timestep=TIMESTEP)
            warm = time.perf_counter() - t0
        finals = np.asarray([float(np.asarray(r.v('b'))[-1]) for r in res])
        return cold, warm, finals
    finally:
        circuit_mod.default_toolkit = saved


def main():
    import jax
    print('devices:', jax.devices())
    print('%6s %12s %14s %14s %9s %12s' % (
        'lanes', 'cpu loop', 'jax cold', 'jax warm', 'speedup',
        'max |dv|'))
    for n in (8, 32, 128, 512):
        cpu_wall, cpu_f = run_cpu(n)
        cold, warm, jax_f = run_jax(n)
        ## Cross-backend agreement on the final sample, per lane.  Both
        ## backends default to Gear-2 (P22), so this is a matched-order
        ## comparison; the bound is loose because the final points land
        ## on different adaptive grids.
        dv = np.max(np.abs(cpu_f - jax_f))
        print('%6d %10.1f ms %12.1f ms %12.1f ms %8.1fx %11.2e' % (
            n, cpu_wall * 1e3, cold * 1e3, warm * 1e3,
            cpu_wall / warm, dv))


if __name__ == '__main__':
    main()
