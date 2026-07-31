"""STAGE 3 gates: the first step.

`doc/transient_work_plan.md` stage 3.  The claim under test is that `reltol` does
not control transient accuracy at all, because the run opens at `max_step` and that
opening step is accepted unevaluated -- so its error is the maximum of the whole run
and no tolerance can move it.

The reference circuit is an RC step response, `v(t) = V(1 - exp(-t/tau))`, because it
has an **exact** solution: the error reported here is the true global error, not a
comparison against another approximation that might be wrong in the same direction.

TWO TRAPS THIS HARNESS IS BUILT TO AVOID, both of which produced a vacuous "pass"
elsewhere in this work:

  * **max_step binding.**  If `timestep` is at or below the step the controller would
    choose, the run sits at `max_step` from end to end, the controller decides
    nothing, and every tolerance gives the same answer.  `assert_controller_binds()`
    checks that the step count is well above `tend/timestep` before any result is
    reported.
  * **starting at the answer.**  Seeded from the DC operating point an RC is already
    settled, so there is no transient and no truncation error to control.  Every run
    here starts from `uic=True`.

Usage::

    python -u benchmarks/transient_stage3.py            # all gates
    python -u benchmarks/transient_stage3.py --gate 1   # one of them
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

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import (EulerIntegrator, TrapezoidalIntegrator,
                                          Gear2Integrator)

V, RES, CAP = 1.0, 1e3, 1e-9
TAU = RES * CAP                 # 1 us
TEND = 10 * TAU
## Deliberately far ABOVE the step the controller would pick, so `max_step` is not
## what is limiting the run.  See the note on max_step binding above.
TIMESTEP = TAU


def build():
    cir = SubCircuit(toolkit=numeric)
    cir['VS'] = VS('n1', gnd, v=V)
    cir['R'] = R('n1', 'n2', r=RES)
    cir['C'] = C('n2', gnd, c=CAP)
    return cir


def run(reltol=1e-4, integrator=None, firststep='unset', **kwargs):
    opts = dict(reltol=reltol, vabstol=1e-12, lte_vabstol=1e-9, uic=True)
    opts['integrator'] = integrator if integrator is not None else TrapezoidalIntegrator()
    if firststep != 'unset':
        opts['firststep'] = firststep
    opts.update(kwargs)
    try:
        res = Transient(build(), toolkit=numeric, **opts).solve(
            refnode=gnd, tend=TEND, timestep=TIMESTEP)
    except KeyError:
        ## `firststep` does not exist yet -- this is the pre-fix measurement.
        opts.pop('firststep', None)
        res = Transient(build(), toolkit=numeric, **opts).solve(
            refnode=gnd, tend=TEND, timestep=TIMESTEP)
    t = np.asarray(res.v('n2').x[0])
    v = np.asarray(res.v('n2').y)
    exact = V * (1.0 - np.exp(-t / TAU))
    return len(t), float(np.max(np.abs(v - exact))), t, v


def assert_controller_binds(steps, label):
    floor = int(TEND / TIMESTEP)
    ok = steps > 2 * floor
    print('    [%s] %d steps against tend/timestep = %d -- controller %s'
          % (label, steps, floor, 'BINDS' if ok else 'IS NOT BINDING'))
    if not ok:
        print('    !! the run is max_step-limited, so no tolerance can be observed')
        print('    !! to act.  Any "pass" from this configuration is vacuous.')
    return ok


def gate_1():
    """reltol must actually control the global error."""
    print('GATE 3-1 -- does reltol control accuracy?')
    print('  RC step response against the analytic solution, trapezoidal.')
    tols = (1e-3, 1e-4, 1e-5, 1e-6)
    rows = []
    for tol in tols:
        n, err, _, _ = run(reltol=tol)
        rows.append((tol, n, err))
    print('  %-10s %8s %14s' % ('reltol', 'steps', 'max err'))
    for tol, n, err in rows:
        print('  %-10g %8d %14.4e' % (tol, n, err))
    assert_controller_binds(rows[-1][1], 'reltol=1e-6')

    first, last = rows[0][2], rows[-1][2]
    ratio = first / last if last > 0 else float('inf')
    monotonic = all(rows[i][2] >= rows[i + 1][2] for i in range(len(rows) - 1))
    print('  error 1e-3 -> 1e-6: %.4e -> %.4e  (%.1fx reduction)'
          % (first, last, ratio))
    print('  monotonic: %s' % monotonic)
    print('  GATE 3-1 (>= 20x and monotonic): %s'
          % ('PASS' if (ratio >= 20 and monotonic) else 'FAIL'))
    return ratio >= 20 and monotonic


def gate_2():
    """The ramp must not cost more than 20% in steps."""
    print('GATE 3-2 -- cost of the ramp, in steps')
    n_off, _, _, _ = run(firststep=TIMESTEP)      # open at max_step, i.e. old behaviour
    n_on, _, _, _ = run()                          # default (ramped once implemented)
    print('  opening at max_step : %d steps' % n_off)
    print('  opening ramped      : %d steps' % n_on)
    growth = (n_on - n_off) / float(n_off) * 100.0
    print('  growth: %+.1f%%' % growth)
    print('  GATE 3-2 (<= +20%%): %s' % ('PASS' if growth <= 20.0 else 'FAIL'))
    return growth <= 20.0


def gate_3():
    """A 2nd-order method must beat backward Euler once the opening step is fixed."""
    print('GATE 3-3 -- does the order of the method show up in the error?')
    print('  %-14s %8s %14s' % ('integrator', 'steps', 'max err'))
    out = {}
    for name, integ in (('euler', EulerIntegrator()),
                        ('trapezoidal', TrapezoidalIntegrator()),
                        ('gear2', Gear2Integrator())):
        n, err, _, _ = run(reltol=1e-5, integrator=integ)
        out[name] = err
        print('  %-14s %8d %14.4e' % (name, n, err))
    better = out['trapezoidal'] < out['euler']
    print('  trapezoidal / euler = %.4f' % (out['trapezoidal'] / out['euler']))
    print('  GATE 3-3 (2nd order beats 1st): %s' % ('PASS' if better else 'FAIL'))
    return better


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gate', type=int, default=0)
    args = p.parse_args()
    gates = {1: gate_1, 2: gate_2, 3: gate_3}
    if args.gate:
        return 0 if gates[args.gate]() else 1
    ok = True
    for i in (1, 2, 3):
        ok = gates[i]() and ok
        print()
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
