"""The measurements that CHANGED A DECISION or REORDERED THE PLAN.

`benchmarks/transient_review/README.md` records that the numerical lens's probe
scripts "were not recovered before the scratchpad was reused", so its results are
reproducible only by reconstruction. This file exists so that does not happen again
to the 2026-07-31 work: every measurement here is one whose *outcome* is quoted in
`doc/transient_work_plan.md` as the reason something was reversed, reordered or
declined. A quoted number whose script is gone is a number nobody can argue with.

Not included: the routine gate measurements, which already live in
`transient_stage2.py`, `transient_stage3.py` and `transient_stage4.py`.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks MPLBACKEND=Agg python3 -u \\
          benchmarks/transient_decisions.py [--parity|--relref|--vabstol|--all]
"""

import argparse
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
for _p in (_REPO, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def parity():
    """WHY 4d IS BLOCKED ON 4g, and the stage was reordered.

    Gate 4d divides by the trapezoidal rule's true truncation error. If that error is
    dominated by the alternating `(-1)^n` mode, the ratio measures 4g and not 4d.

    The test holds the step size and ratio fixed and varies ONLY the number of steps
    taken beforehand. A method whose error depends on step-count *parity* has an
    alternating component; one whose error does not is clean.

    Recorded outcome: euler and gear2 identical to five digits at every prefix,
    trapezoidal swinging 48x between -28.0 and -0.575.
    """
    from pycircuit.circuit.integrator import (EulerIntegrator,
                                              TrapezoidalIntegrator, Gear2Integrator)
    import transient_stage4 as S

    h = 1e-9
    print('TRUE error at a fixed step and ratio, vs the number of PRECEDING steps')
    print('%-16s %s' % ('', ''.join('%14s' % ('prefix=%d' % n) for n in (3, 4, 5, 6, 7))))
    for name, integ in (('euler', EulerIntegrator()),
                        ('trap-classic', TrapezoidalIntegrator(lte_formula='classic')),
                        ('gear2-classic', Gear2Integrator(lte_formula='classic'))):
        row = []
        for n in (3, 4, 5, 6, 7):
            times = [S.T_N - 2 * h - (n - i) * h for i in range(n)]
            times += [S.T_N - 2 * h, S.T_N - h, S.T_N]
            iq = [S.dq(times[0])]
            for i in range(1, len(times)):
                q_c = np.array([S.q(times[i])])
                q_l = np.array([[S.q(times[i - 1])],
                                [S.q(times[i - 2] if i >= 2 else times[i - 1])]])
                iq_l = np.array([[iq[-1]], [iq[-2] if len(iq) >= 2 else iq[-1]]])
                v, _ = integ.compute_derivatives(
                    q_curr=q_c, C_curr=np.array([[1.0]]), h_curr=h, q_last=q_l,
                    iq_last=iq_l, h_last=h, is_first_step=False, toolkit=None)
                iq.append(float(v[0]))
            row.append(iq[-1] - S.dq(times[-1]))
        print('%-16s %s' % (name, ''.join('%14.4e' % v for v in row)))
    print()
    print('A value that flips with prefix parity IS the (-1)^n mode.  Trapezoidal has')
    print('one and it is ~48x larger than the smooth truncation error, so no est/true')
    print('figure for trapezoidal means anything until 4g(b) lands.')


def relref():
    """WHY DECISION D3 WAS REVERTED after being adopted.

    `relref='sigglobal'` was decided, implemented, and sent back by its own gate.
    Referencing the tolerance to the largest signal lets steps grow, and on an
    estimator still carrying the `(-1)^n` mode that is enough to break the
    controller's response to `reltol`: accuracy stops falling monotonically.

    Recorded outcome: euler and both gear2 variants monotone under `sigglobal` and
    needing 1.7-2.5x FEWER steps; both trapezoidal variants non-monotone in ERROR.
    """
    sys.path.insert(0, os.path.join(_REPO, 'pycircuit', 'circuit', 'tests'))
    import test_analysis_transient as T
    from pycircuit.circuit.transient import Transient

    for mode in ('pointlocal', 'sigglobal'):
        for p in Transient.parameters:
            if getattr(p, 'name', '') == 'relref':
                p.default = mode
        print('=== relref=%s' % mode)
        for name in T._INTEGRATORS:
            rows = [T._stiff_run(name, r) for r in (1e-3, 1e-4, 1e-5, 1e-6)]
            print('  %-14s steps %-24s err %s'
                  % (name, [x[0] for x in rows], ['%.2e' % x[2] for x in rows]))
    print()
    print('Look at the ERROR column for trap-*: under sigglobal it RISES when reltol')
    print('tightens.  That is why the default stayed pointlocal until 4g(b).')
    print()
    print('RESOLVED 2026-07-31: re-run after 4g(b) and 4i, every configuration is')
    print('monotone in BOTH columns and `sigglobal` is now the default.  This script')
    print('is kept as the record of why it was reverted once -- and note that the')
    print('step counts here are at equal reltol, which is NOT equal accuracy;')
    print('the matched-accuracy gain is 1.31-2.06x, not the 1.7-2.5x once claimed.')


def vabstol():
    """WHY THE 0.3a TOLERANCE SPLIT COST NOTHING -- the measurement 0.3a asked for.

    Decision 0.3a was raised because `vabstol` 1e-12 -> 1e-6 loosened Newton's node
    convergence by 10^6 as an unmeasured side effect.  Reverting it for Newton could
    have cost iterations; it costs nothing at all.

    Recorded outcome when 0.3a was decided: 4094 steps both ways, endpoint delta
    exactly 0.000e+00.  The absolute step count has since moved (stage 3's opening
    ramp changed the trajectory) -- what matters, and what this script re-measures
    live, is that the delta stays exactly 0.  Newton converges quadratically, so the
    node x-tolerance was never the binding criterion.
    """
    import time
    from pycircuit.circuit import gnd, numeric
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import R, C, L, VPulse
    from pycircuit.circuit.transient import Transient

    def build():
        c = SubCircuit(toolkit=numeric)
        c['VP'] = VPulse(1, gnd, v1=0, v2=5, tr=1e-6, tf=1e-6, pw=10e-6, per=20e-6)
        c['R1'] = R(1, 2, r=1)
        c['L1'] = L(2, 3, L=1e-3)
        c['C1'] = C(3, gnd, c=1e-9)
        return c

    print('%-22s %8s %8s %16s' % ('vabstol (Newton)', 'secs', 'steps', 'v(3) end'))
    out = {}
    for vab in (1e-6, 1e-12):
        best = None
        for _ in range(3):
            t0 = time.perf_counter()
            r = Transient(build(), toolkit=numeric, vabstol=vab).solve(
                refnode=gnd, tend=4e-6, timestep=1e-9)
            dt = time.perf_counter() - t0
            best = dt if best is None else min(best, dt)
        n = len(np.asarray(r.sweep_values))
        v = float(np.asarray(r.v(3, gnd))[-1])
        out[vab] = (n, v)
        print('%-22g %8.2f %8d %16.9f' % (vab, best, n, v))
    print()
    print('steps %d -> %d, endpoint delta %.3e'
          % (out[1e-6][0], out[1e-12][0], out[1e-12][1] - out[1e-6][1]))
    print('A delta of exactly 0 means the 10^6 loosening bought nothing, so')
    print('reverting it for Newton costs nothing.  That is decision 0.3a, measured.')


def main():
    p = argparse.ArgumentParser()
    for name in ('parity', 'relref', 'vabstol', 'all'):
        p.add_argument('--' + name, action='store_true')
    a = p.parse_args()
    chosen = [n for n in ('parity', 'relref', 'vabstol') if getattr(a, n)]
    if a.all or not chosen:
        chosen = ['parity', 'relref', 'vabstol']
    for n in chosen:
        print('=' * 72)
        globals()[n]()
        print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
