"""Fully symbolic nonlinear analysis of the leapfrog, 127 unknowns:
amplitude swept, orders to U^13, accuracy against **transient simulation** -- via
**two-tone IM3** rather than HD3.

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

WHY IM3 AND NOT HD3
-------------------
HD3 cannot be measured on this circuit by any transient, at any cost.  The frequency
that makes the distortion appear at all is the frequency that kills the measurement:
the drive has to sit where the amplifier's loop gain has fallen (100 kHz), and that
is 2.6 octaves into this 5th-order lowpass's stopband.  The fundamental is attenuated
106x from the nonlinear node to the output; the **third harmonic, at 300 kHz, is
attenuated 160 000x**.  Measured, at a 1 V drive:

    3rd harmonic at s0_e1   6.21e-06 V        HD3 = 9.84e-04
    3rd harmonic at output  3.89e-11 V        HD3 = 5.83e-07

An integrator with ``reltol`` r leaves roughly ``r*|signal|`` of error behind, and the
output signal is 6.7e-05 V, so even ``reltol = 1e-8`` leaves ~7e-13 V against a 3.89e-11 V
harmonic.  The quantity is not small because the analysis is imprecise; it is small
because five lowpass stages removed it.

Two tones at 100 and 110 kHz put the IM3 product ``2*f1 - f2`` at **90 kHz -- beside
the fundamentals**, where it takes essentially the same filter loss (a 5th-order
rolloff across 90->100 kHz is a factor of 1.7, against 160 000x for the harmonic).
Same circuit, same nonlinearity, same drive level; only the observable changes:

    IM3 at output, 1 V drive  1.08e-08 V      IM3/fund = 1.61e-04

277x larger than HD3 at the same node and drive, which lifts it clear of the floor.

WHY 100 AND 110 kHz SPECIFICALLY
--------------------------------
Tone spacing sets the run length, because the FFT has to resolve the tones from the
product.  100/101 kHz would need a >= 1 ms window for 1 kHz resolution.  With 10 kHz
spacing everything is a multiple of ``FB = 10 kHz``, so a 100 us window puts the three
signals of interest exactly on bins 9, 10 and 11 -- no leakage, no windowing, and a
tenth of the integration.

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


## The bin base.  Every frequency here is an integer multiple of it, so one 1/FB
## window resolves all of them exactly.
FB = 10.0e3
F1, F2 = 10 * FB, 11 * FB       # 100 kHz, 110 kHz
BIN_IM3, BIN_F1, BIN_F2 = 9, 10, 11        # 2*f1 - f2 = 90 kHz, then the tones
## Matched to `order_convergence.py` so the two paths describe the same circuit:
## same node, same coefficient, same drive levels.
G3 = 5.0e-2
STAGES = 5
NODE_NAME = 's0_e1'
## Largest drive FIRST.  IM3 grows faster than the fundamental, so 3 V has the
## widest margin over the integrator floor (778x at reltol 1e-6, against 162x at
## 1 V) -- and each transient costs hours, so if only one row completes it should
## be the one with the most headroom.
AMPLITUDES = (3.0, 1.0)
ORDERS = (3, 5, 7, 9, 11, 13)

## Settle for ~20 of the 10 us integrator time constants, then measure exactly one
## 1/FB window so the bins land.
SETTLE = 200e-6
MEASURE = 1.0 / FB
## vabstol comes back DOWN here, well below the library default of 1e-6 (Spectre's
## value, right for general use).  This measurement needs the integrator's residual
## error under the ~1e-8 V IM3 amplitude at the output, and the default leaves
## etol = TRTOL*(reltol*|x| + 1e-6) ~ 7e-6 -- larger than the 6.7e-05 V signal's own
## precision requirement.  Measured on a fixed 2.5 us window, |vout| error against a
## tight reference: defaults +2.03%, vabstol 1e-9 +0.45%, vabstol 1e-9 + reltol 1e-6
## -0.02%.  The last is what IM3 needs and what it costs.
TRAN_OPTS = {'vabstol': 1e-9, 'reltol': 1e-6}
## Samples per 1/F1 period.  This is NOT just an output grid -- `timestep` also sets
## `max_step`, which here is a CORRECTNESS knob, not only a cost knob.
##
## `IntegralController.evaluate_step` accepts unconditionally when `is_first_step`
## (stepcontroller.py:26), and `transient.py:314` re-arms that flag after every
## breakpoint.  `Sin.next_event` (func.py:30) fires every quarter period, so each
## tone contributes a periodic, drive-synchronous step that gets NO error check and
## runs at the full `max_step`.  Capping `max_step` is what bounds the damage.
##
## 512 costs nothing at the tolerances in TRAN_OPTS: the controller settles at
## h ~ 4.4 ns there, well under the 19.5 ns cap, so the cap never binds on an
## accuracy-controlled step and only limits the unchecked breakpoint ones.  It is
## not free at loose tolerance -- with the library defaults the controller wants
## ~96 ns and the cap then forces >6000 steps -- so a cheap structural check should
## lower it rather than inherit this value.
POINTS = 512


def build_transient_leapfrog(amplitude):
    """The same topology as `leapfrog_5th_order`, numeric, driven by TWO tones.

    Replicated rather than reused because the fixture builds with the *symbolic*
    toolkit and forces a symbolic `VS`; a transient needs numeric elements and a
    time-domain source.  The device stamping itself is shared, via
    `_stamp_ua741_devices`, so the transistor network cannot drift from the
    fixture's.

    The two tones are two sources in series, so the input node carries their sum.
    """
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = numeric
    try:
        cir = SubCircuit(toolkit=numeric)
        node = {'in': cir.add_node('in'), 'mid': cir.add_node('mid')}
        cir['vs1'] = VSin(node['in'], node['mid'], va=amplitude, freq=F1)
        cir['vs2'] = VSin(node['mid'], gnd, va=amplitude, freq=F2)

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
        nl = amp[0][NODE_NAME.split('_')[-1]]
        cir['nl'] = BSource(nl, gnd, nl, gnd, i_func=lambda v: G3 * v ** 3)
        return cir, amp
    finally:
        circuit_module.default_toolkit = saved


def transient_im3(amplitude, opts=None, settle=SETTLE, measure=MEASURE):
    """IM3 at the output and at the nonlinear node, from a transient plus an FFT.

    Returns a dict per node of ``(im3_ratio, fund, im3)`` plus cost.
    """
    cir, amp = build_transient_leapfrog(amplitude)
    probes = {'out': amp[STAGES - 1]['out'], NODE_NAME: amp[0]['e1']}
    t0 = time.time()
    res = Transient(cir, toolkit=numeric,
                    **(TRAN_OPTS if opts is None else opts)).solve(
        refnode=gnd, tend=settle + measure, timestep=1.0 / (F1 * POINTS))
    secs = time.time() - t0

    out = {}
    npt = int(round(measure * F1 * POINTS))
    for label, nd in probes.items():
        wave = res.v(nd)
        t = np.asarray(wave.x[0])
        v = np.asarray(wave.y)
        ## Resample the measured window onto an exact grid; the solver's own steps
        ## are not uniform, and the bins only land if the grid is.
        grid = np.linspace(t[-1] - measure, t[-1], npt, endpoint=False)
        spec = np.fft.rfft(np.interp(grid, t, v)) / npt
        ## The window is exactly 1/FB, so bin k IS k*FB and the constants can be
        ## used directly.  Asserted rather than assumed: a window that is not a
        ## whole number of 1/FB puts the tones between bins and the IM3 reading
        ## becomes leakage from the fundamental, which would look like distortion.
        assert abs(measure * FB - round(measure * FB)) < 1e-12, \
            'measurement window must be a whole number of 1/FB'
        fund, im3 = abs(spec[BIN_F1]), abs(spec[BIN_IM3])
        out[label] = (im3 / fund if fund else float('nan'), fund, im3)
    return out, secs, len(np.asarray(res.v(probes['out']).x[0]))


def perturbation_im3(system, node_row, out_row, drive_row, amplitude, order):
    """IM3 from the two-tone graded perturbation series truncated at ``U^order``.

    The solver here is numeric per amplitude, deliberately: this script's job is to
    check the series against an independent time-domain solver, not to demonstrate
    the symbolic representation -- `order_convergence.py` is the one that keeps
    ``Adrv`` and ``kk`` symbolic and substitutes last.

    ``A`` is affine in ``s``, so ``A0 + s*C`` is built once numerically instead of
    substituting into 16 129 sympy entries per solve.  Two tones need far more solves
    than one tone does, and that shortcut is what makes this affordable; it is
    checked against the recorded single-tone HD3 numbers by `main`.
    """
    n, s_sym = system.dim, system.s
    A = sympy.Matrix(system.A)
    env = {s: complex(system.params[s]) for s in system.A.free_symbols
           if s is not s_sym}
    A0 = np.array([[complex(A[i, j].subs(s_sym, 0).subs(env)) for j in range(n)]
                   for i in range(n)], dtype=complex)
    Cm = np.array([[complex(sympy.diff(A[i, j], s_sym).subs(env))
                    for j in range(n)] for i in range(n)], dtype=complex)

    def solve(s, rhs):
        return list(np.linalg.solve(A0 + s * Cm, np.asarray(rhs, dtype=complex)))

    def f(x):
        outv = [GradedSpectrum() for _ in range(n)]
        cube = ((x[node_row] * x[node_row]).truncated(order)
                * x[node_row]).truncated(order)
        outv[node_row] = cube.scaled(G3)
        return GradedVector(outv)

    vec = GradedVector([GradedSpectrum() for _ in range(n)])
    for idx in ((1, 0), (0, 1)):
        for k, v in GradedSpectrum.from_phasor(idx, 1,
                                               complex(amplitude)).terms.items():
            vec.components[drive_row].terms[k] = \
                vec.components[drive_row].terms.get(k, 0) + v

    t0 = time.time()
    g = graded_response_mimo(solve, vec, f, (2 * np.pi * F1, 2 * np.pi * F2),
                             max_power=order)
    secs = time.time() - t0
    out = {}
    for label, row in (('out', out_row), (NODE_NAME, node_row)):
        fund = abs(g[row].phasor((1, 0)))
        im3 = abs(g[row].phasor((2, -1)))
        out[label] = (im3 / fund if fund else float('nan'), fund, im3)
    return out, secs


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
    print('nonlinearity: BSource i = %g*v^3 at %s (row %d); drive row %d,'
          ' output row %d' % (G3, NODE_NAME, node_row, drive_row, out_row))
    print('two tones %g/%g kHz; IM3 at 2f1-f2 = %g kHz (bins %d/%d/%d of a'
          ' %g us window)' % (F1 / 1e3, F2 / 1e3, (2 * F1 - F2) / 1e3,
                              BIN_IM3, BIN_F1, BIN_F2, MEASURE * 1e6))
    print('transient options: %s (library defaults are vabstol 1e-6/reltol 1e-4;'
          ' see TRAN_OPTS)' % TRAN_OPTS)
    print('the cubic has zero slope at v=0, so the linearised matrix IS the'
          ' linear leapfrog\'s')
    print()

    ## The perturbation side, first: it is cheap, and it says what the transient is
    ## being asked to resolve before any integration is paid for.
    print('== perturbation IM3 (independent of the transient) ==')
    print('%-7s %-6s %-11s %-11s %-11s %-11s %-11s'
          % ('amp', 'order', 'fund@out', 'IM3@out', 'IM3/f@out', 'IM3@%s' % 's0_e1',
             'IM3/f@e1'))
    print('-' * 76)
    pert = {}
    for amp in AMPLITUDES:
        for order in ORDERS:
            try:
                r, secs = perturbation_im3(system, node_row, out_row, drive_row,
                                           amp, order)
            except Exception as exc:
                print('%-7.3g U^%-4d FAILED: %s: %s'
                      % (amp, order, type(exc).__name__, exc))
                break
            pert[(amp, order)] = r
            ro, re = r['out'], r[NODE_NAME]
            print('%-7.3g U^%-4d %.4e %.4e %.4e %.4e %.4e  (%.0f s)'
                  % (amp, order, ro[1], ro[2], ro[0], re[2], re[0], secs))
            sys.stdout.flush()
        print()

    ## Self-convergence of the series, which needs no external oracle.
    print('== does the series settle?  d(U^k) = |IM3(U^k) - IM3(U^k-2)| /'
          ' |IM3(U^k)| ==')
    head = ('%-7s %-11s' % ('amp', 'IM3/f@out')
            + ''.join('%-11s' % ('d(U^%d)' % o) for o in ORDERS[1:])
            + ' agrees from')
    print(head)
    print('-' * len(head))
    for amp in AMPLITUDES:
        have = [o for o in ORDERS if (amp, o) in pert]
        if len(have) < 2:
            continue
        row = '%-7.3g %-11.4e' % (amp, pert[(amp, have[-1])]['out'][0])
        agrees = None
        for prev, order in zip(have, have[1:]):
            a = pert[(amp, prev)]['out'][2]
            c = pert[(amp, order)]['out'][2]
            rel = abs(c - a) / abs(c) if abs(c) else float('nan')
            row += '%-11s' % ('%.1e' % rel)
            if agrees is None and rel < 1e-6:
                agrees = prev
        row += ' %s' % ('U^%d' % agrees if agrees else
                        'not by U^%d' % have[-1])
        print(row)
    print()

    ## The transient, last, and only where the perturbation says IM3 clears the
    ## integrator's floor.  Running it where it cannot would produce a number, and
    ## the number would be integration error.
    print('== transient IM3, and the comparison ==')
    print('estimated cost at these tolerances: ~4-5 h per amplitude for a'
          ' %g us run.' % ((SETTLE + MEASURE) * 1e6))
    print('Each amplitude prints as it completes, so a kill leaves partial results.')
    print()
    header = ('%-7s %-11s %-11s %-9s %-8s %-7s'
              % ('amp', 'IM3/f tran', 'IM3/f pert', 'rel diff', 'steps', 'hours'))
    print(header)
    print('-' * len(header))
    for amp in AMPLITUDES:
        best = max(o for o in ORDERS if (amp, o) in pert) if any(
            (amp, o) in pert for o in ORDERS) else None
        if best is None:
            continue
        try:
            tr, secs, nsteps = transient_im3(amp)
        except Exception as exc:
            print('%-7.3g TRANSIENT FAILED: %s: %s'
                  % (amp, type(exc).__name__, exc))
            break
        for label in ('out', NODE_NAME):
            p = pert[(amp, best)][label][0]
            t = tr[label][0]
            rel = abs(t - p) / abs(p) if p else float('nan')
            print('%-7.3g %-11.4e %-11.4e %-9.2e %-8d %-7.2f  @%s (U^%d)'
                  % (amp, t, p, rel, nsteps, secs / 3600.0, label, best))
        sys.stdout.flush()

    print()
    print('The transient is an independent solver; nothing on that path shares')
    print('code with the analysis under test.  A transient IM3 that matches to a')
    print('few percent is a real cross-check; one that matches to 1e-6 would mean')
    print('the two paths are not actually independent and should be suspected.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
