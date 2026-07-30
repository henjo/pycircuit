"""Fully symbolic nonlinear analysis of the leapfrog, 127 unknowns:
amplitude swept, orders to U^13, accuracy against **transient simulation**.

Stage 14 attached a nonlinearity *inside* the perturbation solve, as a callable.
Nothing in the circuit was nonlinear, so a transient run would have been perfectly
linear and produced no harmonics -- there would have been nothing to compare against.

Here the nonlinearity is a real element, a `BSource` cubic ``i = g3*v**3`` across one
amplifier's input node, so **the transient engine and the perturbation analysis see
the same circuit**.  Two independent paths, no shared code on the measurement side --
a different solver, a different representation, a different failure mode -- which is
the property the existing distortion-vs-transient gate was built around.

One simplification makes this exact rather than approximate: ``d/dv (g3 v^3) = 0`` at
``v = 0``, so the cubic contributes nothing to the linearised matrix.  The matrix the
perturbation solves against is therefore *exactly* the linear leapfrog's, and the
nonlinearity enters only through ``f(x)``.  No separate linearisation is needed and
none can drift.

The transient is the expensive half; it is cached per amplitude and its cost is
printed as it goes, so a run that is going to take too long says so early.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/nonlinear_leapfrog_sweep.py
"""

import sys
import time

import numpy as np
import sympy

from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit import circuit as circuit_module
from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.benchmark_circuits import (_UA741_NODE_NAMES,
                                                 _UA741_ROLES,
                                                 _stamp_ua741_devices,
                                                 add_small_signal_bjt)
from pycircuit.circuit.elements import BSource, C, R, VSin
from pycircuit.circuit.distortion import (GradedSpectrum, GradedVector,
                                          graded_response_mimo)
from pycircuit.circuit.transient import Transient


F0 = 1.0e3                      # in the leapfrog's passband (cutoff ~16 kHz)
G3 = 5.0e-3                     # cubic coefficient of the added nonlinearity
STAGES = 5
AMPLITUDES = (1e-3, 3e-3, 1e-2, 3e-2, 1e-1)
ORDERS = (3, 5, 7, 9, 11, 13)


def build_transient_leapfrog(amplitude):
    """The same topology as `leapfrog_5th_order`, numeric, driven in time.

    Replicated rather than reused because the fixture builds with the *symbolic*
    toolkit and forces a symbolic `VS`; a transient needs numeric elements and a
    time-domain source.  The device stamping itself is shared, via
    `_stamp_ua741_devices`, so the transistor network cannot drift from the
    fixture's.
    """
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = numeric
    try:
        cir = SubCircuit(toolkit=numeric)
        node = {'in': cir.add_node('in')}
        cir['vs'] = VSin(node['in'], gnd, va=amplitude, freq=F0)

        amp = []
        for k in range(STAGES):
            local = {}
            for base in _UA741_NODE_NAMES:
                local[base] = cir.add_node('s%d_%s' % (k, base))
            amp.append(local)

        def bjt_for(stage):
            def bjt(dev, role, b, c, e):
                gm, rpi, ro, cpi, cmu = _UA741_ROLES[role]
                add_small_signal_bjt(cir, 's%d_%s' % (stage, dev), b, c, e,
                                     gm, rpi, ro, cpi, cmu)
            return bjt

        for k in range(STAGES):
            _stamp_ua741_devices(cir, amp[k], bjt_for(k), True, pfx='s%d_' % k)
            cir['sg%d' % k] = R(amp[k]['inp'], gnd, r=1.0)
        for k in range(STAGES):
            cir['ci%d' % k] = C(amp[k]['out'], amp[k]['inn'], c=1e-9)
        cir['rd0'] = R(amp[0]['out'], amp[0]['inn'], r=10e3)
        cir['rd%d' % (STAGES - 1)] = R(amp[STAGES - 1]['out'],
                                       amp[STAGES - 1]['inn'], r=10e3)
        cir['rin'] = R(node['in'], amp[0]['inn'], r=10e3)
        for k in range(1, STAGES):
            cir['rf%d' % k] = R(amp[k - 1]['out'], amp[k]['inn'], r=10e3)
        for k in range(STAGES - 1):
            cir['rb%d' % k] = R(amp[k + 1]['out'], amp[k]['inn'], r=10e3)

        ## THE nonlinearity, and the only nonlinear element in the circuit: a
        ## cubic conductance to ground on stage 0's input-pair emitter.
        nl = amp[0]['e1']
        cir['nl'] = BSource(nl, gnd, nl, gnd, i_func=lambda v: G3 * v ** 3)
        return cir, amp
    finally:
        circuit_module.default_toolkit = saved


def transient_hd3(amplitude, periods=20, points=256, keep=8):
    """HD3 at the output from a transient run and an FFT of the steady state."""
    cir, amp = build_transient_leapfrog(amplitude)
    out = amp[STAGES - 1]['out']
    t0 = time.time()
    res = Transient(cir, toolkit=numeric).solve(
        refnode=gnd, tend=periods / F0, timestep=1.0 / (F0 * points))
    secs = time.time() - t0
    wave = res.v(out)
    t = np.asarray(wave.x[0])
    v = np.asarray(wave.y)
    ## Resample the last whole periods onto an exact grid so FFT bins land on
    ## the harmonics; the solver's own steps are not uniform.
    grid = np.linspace(t[-1] - keep / F0, t[-1], keep * points, endpoint=False)
    spec = np.fft.rfft(np.interp(grid, t, v)) / (keep * points)
    k = keep
    fund, third = abs(spec[k]), abs(spec[3 * k])
    return (third / fund if fund else float('nan')), fund, third, secs, len(t)


def perturbation_hd3(system, node_row, out_row, drive_row, amplitude, order):
    """HD3 from the graded perturbation series truncated at ``U^order``."""
    n, s_sym = system.dim, system.s
    env = {}
    for symbol in system.A.free_symbols:
        if symbol is not s_sym:
            env[symbol] = complex(system.params[symbol])
    rows = sympy.Matrix(system.A)

    def solve(s, rhs):
        local = dict(env)
        local[s_sym] = s
        M = np.array([[complex(rows[i, j].subs(local)) for j in range(n)]
                      for i in range(n)], dtype=complex)
        return list(np.linalg.solve(M, np.asarray(rhs, dtype=complex)))

    def f(x):
        outv = [GradedSpectrum() for _ in range(n)]
        cube = ((x[node_row] * x[node_row]).truncated(order)
                * x[node_row]).truncated(order)
        outv[node_row] = cube.scaled(G3)
        return GradedVector(outv)

    def source(value):
        vec = GradedVector([GradedSpectrum() for _ in range(n)])
        vec.components[drive_row].terms.update(
            GradedSpectrum.from_phasor(1, 1, value).terms)
        return vec

    t0 = time.time()
    g = graded_response_mimo(solve, source(complex(amplitude)), f,
                             (2 * np.pi * F0,), max_power=order)
    fund = abs(g[out_row].phasor(1))
    third = abs(g[out_row].phasor(3))
    return (third / fund if fund else float('nan')), fund, third, time.time() - t0


def main():
    t0 = time.time()
    system = bc.leapfrog_5th_order()
    ## Row indices: the fixture's node order is 'in' then s<k>_<name>.
    names = ['in'] + ['s%d_%s' % (k, b) for k in range(STAGES)
                      for b in _UA741_NODE_NAMES]
    node_row = names.index('s0_e1')
    out_row = system.out_index
    drive_row = [i for i in range(system.dim) if system.b[i] != 0][0]
    print('leapfrog: %d unknowns, built in %.1f s' % (system.dim,
                                                      time.time() - t0))
    print('nonlinearity: BSource i = %g*v^3 at s0_e1 (row %d); drive row %d,'
          ' output row %d; f0 = %g Hz' % (G3, node_row, drive_row, out_row, F0))
    print('the cubic has zero slope at v=0, so the linearised matrix IS the'
          ' linear leapfrog\'s')
    print()

    header = ('%-9s %-11s %-11s' % ('amplitude', 'HD3 trans', 'transient s')
              + ''.join('%-13s' % ('U^%d' % o) for o in ORDERS))
    print(header)
    print('-' * len(header))

    for amp in AMPLITUDES:
        try:
            hd3_t, fund_t, third_t, secs, nsteps = transient_hd3(amp)
        except Exception as exc:
            print('%-9.3g TRANSIENT FAILED: %s: %s'
                  % (amp, type(exc).__name__, exc))
            break
        row = '%-9.3g %-11.4e %-11.1f' % (amp, hd3_t, secs)
        for order in ORDERS:
            try:
                hd3_p, _f, _th, psecs = perturbation_hd3(
                    system, node_row, out_row, drive_row, amp, order)
                rel = abs(hd3_p - hd3_t) / hd3_t if hd3_t else float('nan')
                row += '%-13s' % ('%.2e' % rel)
            except Exception as exc:
                row += '%-13s' % ('ERR', )
        print(row)
        sys.stdout.flush()

    print()
    print('cells are |HD3(perturbation) - HD3(transient)| / HD3(transient).')
    print('The transient is an independent solver; nothing on that path shares')
    print('code with the analysis under test.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
