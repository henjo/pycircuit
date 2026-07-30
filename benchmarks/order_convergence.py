"""At which perturbation order does the leapfrog's distortion stop changing?

Symbolic throughout; numbers substituted only at the end.

The transient comparison was abandoned on a measured cost wall: the leapfrog is stiff
-- time constants from ~0.5 ns (the uA741's 27 ohm output degeneration into 20 pF) to
10 us (10 kohm x 1 nF integrators), with 1 ohm input tie-downs -- so 2 ms of simulated
time did not complete in 6 minutes and the table needed hundreds of times more.  So
the question becomes self-convergence, which needs no external oracle: at which
truncation order does the third harmonic stop moving?

**The analysis stays symbolic and the amplitude is a symbol too.**  That is not a
constraint here, it is the whole economy of the thing: with the drive amplitude
carried as ``Adrv`` and the nonlinear coefficient as ``kk``, **one symbolic build per
order serves the entire amplitude sweep**.  The expensive part happens once; every
amplitude is then a substitution into the same expression graph.  Substituting numbers
into the circuit first would give up exactly that.

The representation is `distortion_ddd.Expr`, a hash-consed straight-line program over
the device symbols -- verified against an independent fully-numeric graded solve to
4.6e-13 on this circuit in stage 14, so what is being measured here is the *series*,
not the representation.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/order_convergence.py
"""

import sys
import time

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.benchmark_circuits import _UA741_NODE_NAMES
from pycircuit.circuit.distortion import (GradedSpectrum, GradedVector,
                                          graded_response_mimo)
from pycircuit.circuit.ddd import _resolve
from pycircuit.circuit.distortion_ddd import Expr, evaluate_one, nodes_of

from nonlinear_leapfrog import expr_solver


## THE FREQUENCY MATTERS AS MUCH AS THE AMPLITUDE.  At 1 kHz the amplifier's loop
## gain is high, the input pair is held at a virtual ground (2.5 uV at a 10 mV drive),
## and distortion is suppressed -- every amplitude converged at U^3 with HD3 ~ 1e-14,
## which is physically right and tells you nothing about the series.  Distortion rises
## where loop gain falls.  The uA741 here measures +87.9 dB at 10 Hz falling at
## -20 dB/decade, so at 100 kHz its open-loop gain is down to roughly 8 dB and the
## virtual ground no longer holds.  That is where the nonlinearity acts and where the
## higher orders are actually exercised.
F0 = 1.0e5
STAGES = 5
ORDERS = (3, 5, 7, 9, 11, 13, 15, 17)
AMPLITUDES = (1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0)
## Back on the input differential pair, which is where an op-amp's distortion
## physically originates.  It was useless at 1 kHz only because the loop gain held it
## at a virtual ground; with the loop gain reduced it sees real signal.
NODE_NAME = 's0_e1'
## kk is a SYMBOL, so sweeping it costs nothing after the build.  That is the whole
## point of substituting numbers last, and it is what lets one build answer both
## "how does order matter with amplitude" and "how does it matter with strength".
KK_VALUES = (5e-3, 5e-1, 5e+1)
G3_VALUE = 5.0e-2


def turning_point(system, node_row, kk, freq):
    """``v_turn = 1/sqrt(3|a|)`` with ``a = kk/g``, the cubic's validity limit.

    A truncated cubic ``i = g(v + a v^3)`` has negative differential conductance beyond
    this and **is not a physical device** (see `doc/distortion_mimo_plan.md` 8.3).  An
    order-convergence study is only meaningful inside it: without this number a
    divergence caused by the *model* looks exactly like one caused by the *series*,
    which cost a withdrawn table row.
    """
    env = {}
    for sym in system.A.free_symbols:
        if sym is not system.s:
            env[sym] = complex(system.params[sym])
    env[system.s] = 2j * np.pi * freq
    n = system.dim
    M = np.array([[complex(_resolve(system.A[i, j], env)) for j in range(n)]
                  for i in range(n)], dtype=complex)
    g = 1.0 / abs(np.linalg.inv(M)[node_row, node_row])
    a = kk / g
    return 1.0 / np.sqrt(3 * abs(a)), g


def build_symbolic(system, node_row, out_row, drive_row, order):
    """One symbolic build for a whole amplitude sweep.

    ``Adrv`` and ``kk`` stay symbols, so the returned graphs are functions of the
    drive amplitude and the nonlinearity strength.  Substituting either is then a
    walk over the graph, not a rebuild.
    """
    n, s_sym = system.dim, system.s
    drive = sympy.Symbol('Adrv')
    kk = sympy.Symbol('kk')

    def f(x):
        out = [GradedSpectrum() for _ in range(n)]
        cube = ((x[node_row] * x[node_row]).truncated(order)
                * x[node_row]).truncated(order)
        out[node_row] = cube.scaled(Expr.atom(kk))
        return GradedVector(out)

    vec = GradedVector([GradedSpectrum() for _ in range(n)])
    vec.components[drive_row].terms.update(
        GradedSpectrum.from_phasor(1, 1, Expr.atom(drive)).terms)

    t0 = time.time()
    g = graded_response_mimo(expr_solver(system.A, s_sym, n), vec, f,
                             (2 * np.pi * F0,), max_power=order)
    fund = g[out_row].phasor(1)
    third = g[out_row].phasor(3)
    nodev = g[node_row].phasor(1)
    return {'fund': fund, 'third': third, 'nodev': nodev,
            'build': time.time() - t0,
            'nodes': nodes_of(fund, third, nodev),
            'drive': drive, 'kk': kk}


def main():
    t0 = time.time()
    system = bc.leapfrog_5th_order()
    names = ['in'] + ['s%d_%s' % (k, b) for k in range(STAGES)
                      for b in _UA741_NODE_NAMES]
    node_row = names.index(NODE_NAME)
    out_row = system.out_index
    drive_row = [i for i in range(system.dim) if system.b[i] != 0][0]
    print('leapfrog: %d unknowns, built in %.1f s' % (system.dim,
                                                      time.time() - t0))
    print('nonlinearity i = kk*v^3 at %s (row %d), drive Adrv at row %d,'
          ' output row %d, f0 = %g Hz'
          % (NODE_NAME, node_row, drive_row, out_row, F0))
    print('Adrv and kk stay SYMBOLIC: one build per order, then substitution.')
    ## Evidence that the loop gain really has reduced: the node the loop holds at a
    ## virtual ground when the gain is high should now carry real signal.
    print()

    base = {}
    for symbol in system.A.free_symbols:
        if symbol is not system.s:
            base[symbol] = complex(system.params[symbol])

    built = {}
    print('%-6s %-12s %-11s %s' % ('order', 'graph nodes', 'build s',
                                   'evaluate s (6 amplitudes)'))
    print('-' * 62)
    for order in ORDERS:
        try:
            b = build_symbolic(system, node_row, out_row, drive_row, order)
        except Exception as exc:
            print('U^%-4d FAILED: %s: %s' % (order, type(exc).__name__, exc))
            break
        t1 = time.time()
        vals = {}
        for amp in AMPLITUDES:
            env = dict(base)
            env[b['drive']] = complex(amp)
            env[b['kk']] = complex(G3_VALUE)
            vals[amp] = (evaluate_one(b['fund'], env),
                         evaluate_one(b['third'], env),
                         evaluate_one(b['nodev'], env))
        ev = time.time() - t1
        b['vals'] = vals
        built[order] = b
        print('U^%-4d %-12d %-11.1f %.2f' % (order, b['nodes'], b['build'], ev))
        sys.stdout.flush()

    if len(built) < 2:
        print('\nnot enough orders completed to compare')
        return 1

    done = sorted(built)
    print()
    ## system.dim, NOT a literal.  This header read "127 unknowns" while printing a
    ## 136-unknown run, because the fixture gained nine nodes from GBW compensation and
    ## a pasted number is correct exactly once -- in the script whose output is pasted
    ## into doc/distortion_ddd_conclusions.md, which is the worst place for it.
    print('== Fully symbolic nonlinear analysis of the leapfrog, %d unknowns =='
          % system.dim)
    print('Third harmonic at the output; one symbolic build per order, evaluated')
    print('at each amplitude.  d(U^k) = |H3(U^k) - H3(U^k-2)| / |H3(U^k)|.')
    print()
    vturn, gnode = turning_point(system, node_row, G3_VALUE, F0)
    print('cubic validity limit: g at %s = %.4e S, a = kk/g = %.4e,'
          ' v_turn = %.4e V' % (NODE_NAME, gnode, G3_VALUE / gnode, vturn))
    print('rows beyond 100% of v_turn are NOT a physical device and are marked.')
    print()
    head = ('%-9s %-12s %-9s %-12s' % ('amp (V)', '|v| %s' % NODE_NAME,
                                       '%v_turn', 'HD3')
            + ''.join('%-11s' % ('d(U^%d)' % o) for o in done[1:])
            + ' agrees from')
    print(head)
    print('-' * len(head))
    for amp in AMPLITUDES:
        top = built[done[-1]]['vals'][amp]
        hd3 = abs(top[1]) / abs(top[0]) if abs(top[0]) else float('nan')
        frac = 100.0 * abs(top[2]) / vturn
        row = '%-9.3g %-12.4e %-9s %-12.4e' % (
            amp, abs(top[2]),
            ('%.0f%%' % frac) + ('!' if frac > 100 else ''), hd3)
        agrees = None
        for prev, order in zip(done, done[1:]):
            a = built[prev]['vals'][amp][1]
            c = built[order]['vals'][amp][1]
            rel = abs(c - a) / abs(c) if abs(c) else float('nan')
            row += '%-11s' % ('%.1e' % rel)
            if agrees is None and rel < 1e-6:
                agrees = prev
        row += ' %s' % ('U^%d' % agrees if agrees else 'not by U^%d' % done[-1])
        print(row)
    print()
    print('One symbolic build per order serves every amplitude: the build column')
    print('is paid once, the evaluate column covers all six substitutions.')

    ## The free second sweep: kk is a symbol, so this needs no rebuild at all.
    print()
    print('== and because kk is symbolic, sweeping the nonlinearity is free ==')
    print('Convergence order at a 1 V drive, as the cubic coefficient grows:')
    print('%-12s %-11s %-9s %s' % ('kk', 'HD3', '%v_turn', 'agrees from'))
    print('-' * 40)
    for kk in KK_VALUES:
        env = dict(base)
        env[built[done[-1]]['drive']] = 1.0 + 0j
        vals = {}
        for order in done:
            e = dict(base)
            e[built[order]['drive']] = 1.0 + 0j
            e[built[order]['kk']] = complex(kk)
            vals[order] = (evaluate_one(built[order]['fund'], e),
                           evaluate_one(built[order]['third'], e))
        top = vals[done[-1]]
        hd3 = abs(top[1]) / abs(top[0]) if abs(top[0]) else float('nan')
        vt, _g = turning_point(system, node_row, kk, F0)
        agrees = None
        for prev, order in zip(done, done[1:]):
            a, c = vals[prev][1], vals[order][1]
            rel = abs(c - a) / abs(c) if abs(c) else float('nan')
            if agrees is None and rel < 1e-6:
                agrees = prev
        ## The node voltage barely moves with kk at fixed drive, so the 1 V
        ## amplitude row's value is a fair stand-in for the fraction.
        nodev = abs(built[done[-1]]['vals'][1.0][2])
        frac = 100.0 * nodev / vt
        print('%-12.3g %-11.4e %-9s %s'
              % (kk, hd3, ('%.0f%%' % frac) + ('!' if frac > 100 else ''),
                 ('UNPHYSICAL -- past v_turn' if frac > 100 else
                  ('U^%d' % agrees if agrees else 'not by U^%d' % done[-1]))))
    return 0


if __name__ == '__main__':
    sys.exit(main())
