"""STAGE 0.2a -- is the IM3 harness's integrator choice contaminated by breakpoints?

`doc/transient_work_plan.md` 0.2a.  `benchmarks/nonlinear_leapfrog_sweep.py` sets
`TrapezoidalIntegrator()` on the strength of a step-count comparison:

    euler        2896 steps      trapezoidal   288 steps  <- 10x, and the reason for the choice
    gear2 ywr     347 steps      gear2 classic 376 steps

That comparison was taken on a two-tone `VSin` drive, which is the exact configuration
that triggers the seeding mechanism in stage 4g: `Sin.next_event` plants a breakpoint
every quarter period on a C-infinity waveform, and `transient.py` re-arms the order drop
after every breakpoint.  A trapezoidal method restarted periodically is precisely how the
`(-1)^n` parasitic mode is fed.  So the step count that chose the integrator could be an
artefact of the defect the integrator choice is supposed to be independent of.

THE QUESTION, narrowly.  Not "is Trapezoidal the best integrator" -- that is gate 4f.
Only: **does removing the breakpoints change the ranking?**  If it does not, the recorded
10x stands and stage 4 may proceed on it.  If it does, the 10x is void and must be
re-measured after 4g.

WHAT IS COMPARED.  Five configurations over one fixed window, from the same initial
state, same tolerances, same `max_step`:

    A  Trapezoidal, breakpoints as shipped   (what the harness runs today)
    B  Gear2 'ywr', breakpoints as shipped   (the shipped default integrator)
    C  Gear2 'classic', breakpoints as shipped
    D  Trapezoidal, breakpoints SUPPRESSED   (stage 4g(a) applied)
    E  Gear2 'ywr', breakpoints SUPPRESSED

Suppression is exactly what 4g(a) proposes -- `Sin.next_event` returns `inf` instead of
the next quarter period -- applied by monkeypatch here so that no shipped behaviour
changes while a stage-0 question is being answered.

WHY A SHORT WINDOW.  The harness's own run is `SETTLE = 5*208 us` plus a `1/FB = 100 us`
measurement window, which at these tolerances is hours per configuration and five
configurations is not affordable for a stage-0 gate.  The comparison that chose the
integrator was itself taken "over a fixed 2.5 us window" (the harness comment), so this
reproduces that window and reports steps-per-microsecond, which is the quantity the 10x
claim is actually about.  The waveform check is therefore a *transient-accuracy*
comparison and NOT an IM3 comparison; that limitation is reported with the result rather
than papered over.  See the LIMITATION note printed at the end.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks MPLBACKEND=Agg python3 -u \
        benchmarks/transient_review/stage0_2a_integrator_choice.py
"""

import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
for _p in (_REPO, os.path.join(_REPO, 'benchmarks')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import TrapezoidalIntegrator, Gear2Integrator
from pycircuit.circuit import func as func_module

import nonlinear_leapfrog_sweep as H


## The window the harness's own integrator comparison used.  2.5 us is a QUARTER of one
## 100 kHz period, which is short, and that is what the original comparison used: the
## cost here is set by resolving the Q ~ 16.8 resonance rather than by the drive period,
## and the controller reaches its steady step size well inside it.  Note what this
## implies for the breakpoint question -- `Sin.next_event` fires every quarter period,
## so this window spans about ONE breakpoint per tone.  It is therefore the shortest
## window in which a breakpoint effect can appear at all, and it UNDERSTATES a
## cumulative one.  If the effect is already visible here it is real; if it is not, the
## honest reading is "not visible in one breakpoint interval", not "absent".
WINDOW = 2.5e-6
AMPLITUDE = 1.0


def _steps(res, node):
    return len(np.asarray(res.v(node).x[0]))


def run_one(label, integrator_factory, suppress_breakpoints):
    """One configuration.  Returns (label, steps, secs, t, v_out, v_nl)."""
    saved_next_event = func_module.Sin.next_event
    if suppress_breakpoints:
        ## Stage 4g(a), in the minimal form: a C-infinity source has no
        ## discontinuity to land on, so it contributes no breakpoint.  The shipped
        ## fix also returns `td` once, to catch the start of a delayed source; this
        ## harness drives from t=0 with td=0, so `inf` alone is the same thing.
        func_module.Sin.next_event = lambda self, t: float('inf')
    try:
        cir, amp = H.build_transient_leapfrog(AMPLITUDE)
        probes = {'out': amp[H.STAGES - 1]['out'],
                  'nl': amp[0]['e1']}
        opts = dict(H.TRAN_OPTS)
        opts['integrator'] = integrator_factory()
        t0 = time.time()
        res = Transient(cir, toolkit=numeric, **opts).solve(
            refnode=gnd, tend=WINDOW, timestep=1.0 / (H.F1 * H.POINTS))
        secs = time.time() - t0
        wave = res.v(probes['out'])
        t = np.asarray(wave.x[0])
        v_out = np.asarray(wave.y)
        v_nl = np.asarray(res.v(probes['nl']).y)
        return dict(label=label, steps=len(t), secs=secs,
                    t=t, v_out=v_out, v_nl=v_nl)
    finally:
        func_module.Sin.next_event = saved_next_event


def compare(ref, other, npt=2001):
    """Max abs difference of the two waveforms on a common uniform grid.

    Both runs choose their own non-uniform steps, so a pointwise comparison is not
    available; interpolating onto a shared grid is the same thing the harness does
    before its FFT, and the same hazard applies -- this bounds the difference, it does
    not resolve it below the interpolation error.
    """
    lo = max(ref['t'][0], other['t'][0])
    hi = min(ref['t'][-1], other['t'][-1])
    grid = np.linspace(lo, hi, npt)
    out = {}
    for key in ('v_out', 'v_nl'):
        a = np.interp(grid, ref['t'], ref[key])
        b = np.interp(grid, other['t'], other[key])
        scale = max(np.max(np.abs(a)), 1e-30)
        out[key] = (np.max(np.abs(a - b)), np.max(np.abs(a - b)) / scale)
    return out


def main():
    configs = [
        ('A  trapezoidal      bp on ', lambda: TrapezoidalIntegrator(), False),
        ('B  gear2 ywr        bp on ', lambda: Gear2Integrator(), False),
        ('C  gear2 classic    bp on ',
         lambda: Gear2Integrator(lte_formula='classic'), False),
        ('D  trapezoidal      bp OFF', lambda: TrapezoidalIntegrator(), True),
        ('E  gear2 ywr        bp OFF', lambda: Gear2Integrator(), True),
    ]

    print('STAGE 0.2a -- does breakpoint seeding decide the integrator choice?')
    print('window %.3g s, amplitude %g, %s' % (WINDOW, AMPLITUDE, H.TRAN_OPTS))
    print('max_step = timestep = %.4g s (POINTS=%d)'
          % (1.0 / (H.F1 * H.POINTS), H.POINTS))
    print()

    results = []
    for label, factory, suppress in configs:
        try:
            r = run_one(label, factory, suppress)
        except Exception as exc:                       # noqa: BLE001
            print('%s  FAILED: %s: %s' % (label, type(exc).__name__, exc))
            results.append(dict(label=label, steps=None, secs=None, exc=exc))
            continue
        results.append(r)
        print('%s  %6d steps  %8.1f s  %8.1f steps/us'
              % (label, r['steps'], r['secs'], r['steps'] / (WINDOW * 1e6)))
        sys.stdout.flush()

    ok = [r for r in results if r.get('steps') is not None]
    if len(ok) < 2:
        print('\nToo few configurations completed to compare.')
        return 1

    print()
    print('STEP-COUNT RATIOS (the quantity the 10x claim is about)')
    base = ok[0]
    for r in ok[1:]:
        print('  %s / %s : %.3f'
              % (r['label'].strip(), base['label'].strip(),
                 r['steps'] / base['steps']))

    print()
    print('WAVEFORM DIFFERENCE against A (max abs, and relative to signal peak)')
    for r in ok[1:]:
        d = compare(base, r)
        print('  %s  out %.3e (%.2e rel)   nl %.3e (%.2e rel)'
              % (r['label'].strip(), d['v_out'][0], d['v_out'][1],
                 d['v_nl'][0], d['v_nl'][1]))

    ## THE GATE, evaluated in the script rather than by eye.
    print()
    print('GATE 0.2a')
    by = {r['label'][0]: r for r in ok}
    verdict = []
    if 'A' in by and 'D' in by:
        ratio = by['D']['steps'] / by['A']['steps']
        print('  breakpoint effect on trapezoidal: D/A = %.3f  (within 20%%: %s)'
              % (ratio, 'YES' if 0.8 <= ratio <= 1.25 else 'NO'))
        verdict.append(0.8 <= ratio <= 1.25)
    if 'B' in by and 'E' in by:
        ratio = by['E']['steps'] / by['B']['steps']
        print('  breakpoint effect on gear2 ywr:   E/B = %.3f  (within 20%%: %s)'
              % (ratio, 'YES' if 0.8 <= ratio <= 1.25 else 'NO'))
        verdict.append(0.8 <= ratio <= 1.25)
    if 'A' in by and 'B' in by and 'D' in by and 'E' in by:
        before = by['B']['steps'] / by['A']['steps']
        after = by['E']['steps'] / by['D']['steps']
        print('  RANKING trapezoidal-vs-gear2: %.3f with breakpoints, %.3f without'
              % (before, after))
        print('  -> the choice %s'
              % ('SURVIVES breakpoint removal' if (before > 1) == (after > 1)
                 else 'REVERSES when breakpoints are removed'))

    print()
    print('LIMITATION, stated rather than buried: this measures step count and')
    print('transient waveform over %.3g s, NOT the IM3 figure itself, which needs the'
          % WINDOW)
    print('harness\'s full %.3g s settle.  A ranking that survives here is necessary'
          % H.SETTLE)
    print('for the 10x to stand, not sufficient -- the IM3 number itself is only')
    print('re-validated by rerunning the harness, which gate 4g-3 requires anyway.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
