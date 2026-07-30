"""STAGE 0.1d -- what `_solve_coupled` does when Newton fails to converge.

`doc/transient_work_plan.md` 0.1d, the numerical-analyst lens on the third
transcription of the time loop (`transient.py:429-551`).

Reading the code, `:499-526` looks like a bounded retry: on `NoConvergenceError`
shrink `h_curr` by 0.25 and try again, at most `MAX_LTE_ITERS = 10` times, with a
`RuntimeError` guard if the step falls below `minstep`.  It is not bounded, and the
guard does not fire.  Two things combine:

  * after 10 failures the `for` loop simply EXITS -- there is no `else` clause and no
    raise -- and `:528-530` then advance time by the collapsed `h_curr` and append the
    PREVIOUS solution as though it were a new one;
  * `h` for the next outer iteration is restored from `h_next` at `:543`, which the
    failure path never updates, so the next iteration starts again at FULL step size.

So each outer iteration costs 10 Newton solves and buys `h * 0.25^10` of simulated
time, forever.  `minstep` is 1e-18 and `0.25^10` of a microsecond is 9.5e-13, so the
`RuntimeError` at `:508` is never reached.

This script measures the progress rate and checks it against that predicted mechanism,
so the diagnosis is confirmed rather than inferred.  It also shows the different --
and not better -- behaviour when the failure happens on the FIRST step, where the run
dies with an `AttributeError` on the `_iq` side channel instead.

NOTE ON THE ARITHMETIC.  The rate must be measured INCREMENTALLY, between two marks
well inside the stalled region.  An average taken from t=0 still contains the
successful steps and is meaningless; the first version of this probe made exactly that
mistake and reported a finishing time five orders of magnitude too optimistic.

Run:  PYTHONPATH=<repo> MPLBACKEND=Agg python3 -u \
        benchmarks/transient_review/stage0_1d_coupled_livelock.py
"""

import sys
import time

import numpy as np

from pycircuit.circuit import SubCircuit, R, C, VS, gnd, numeric
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.nrsolver import NoConvergenceError

TEND, TSTEP = 5e-6, 1e-6
MAX_LTE_ITERS = 10          # transient.py:487, mirrored here to derive the rate
MARK_A, MARK_B, STOP = 500, 2500, 2501


class Stop(Exception):
    """Raised to end the grind once enough of it has been measured."""


def build():
    cir = SubCircuit(toolkit=numeric)
    cir['VS'] = VS('n1', gnd, v=1.0)
    cir['R'] = R('n1', 'n2', r=1e3)
    cir['C'] = C('n2', gnd, c=1e-9)
    return cir


def midrun_failure():
    """Fail from the 4th call onward, i.e. once history exists."""
    tran = Transient(build(), toolkit=numeric)
    real = tran.solve_timestep
    st = {'n': 0, 'marks': {}}

    def flaky(xlast, tnext, **kw):
        st['n'] += 1
        if st['n'] > 3:
            if st['n'] in (MARK_A, MARK_B):
                st['marks'][st['n']] = tnext
            if st['n'] >= STOP:
                st['t_stop'] = tnext
                raise Stop
            raise NoConvergenceError('probe: forced non-convergence')
        return real(xlast, tnext, **kw)

    tran.solve_timestep = flaky

    print('MID-RUN non-convergence   tend=%.1e  timestep=%.1e' % (TEND, TSTEP))
    try:
        tran.solve(refnode=gnd, tend=TEND, timestep=TSTEP, coupled_lte=True)
        print('  completed normally -- the defect did NOT reproduce')
        return False
    except Stop:
        ta, tb = st['marks'][MARK_A], st['marks'][MARK_B]
        douter = (MARK_B - MARK_A) / float(MAX_LTE_ITERS)
        per_outer = (tb - ta) / douter
        predicted = TSTEP * 0.25 ** MAX_LTE_ITERS
        print('  stalled at simulated t = %.6e s of %.1e  (%.1f%% of the run)'
              % (tb, TEND, 100 * tb / TEND))
        print('  between call %d and %d, simulated time advanced %.4e s'
              % (MARK_A, MARK_B, tb - ta))
        print('  = %.0f outer iterations -> %.4e s each' % (douter, per_outer))
        print()
        print('  MECHANISM CHECK  h*0.25^%d = %.4e s predicted'
              % (MAX_LTE_ITERS, predicted))
        print('                   %.4e s measured -- ratio %.4f'
              % (per_outer, per_outer / predicted))
        print()
        print('  outer iterations still needed : %.3e'
              % ((TEND - tb) / per_outer))
        print('  further Newton attempts       : %.3e'
              % (MAX_LTE_ITERS * (TEND - tb) / per_outer))
        print('  -> it neither raises nor returns.')
        return True
    except Exception as exc:                                   # noqa: BLE001
        print('  raised %s: %s' % (type(exc).__name__, exc))
        return False


def firststep_failure():
    """Fail from the very first call, before any history exists."""
    tran = Transient(build(), toolkit=numeric)
    st = {'n': 0}

    def always_fail(*a, **kw):
        st['n'] += 1
        raise NoConvergenceError('probe: forced non-convergence')

    tran.solve_timestep = always_fail

    print()
    print('FIRST-STEP non-convergence')
    try:
        tran.solve(refnode=gnd, tend=TEND, timestep=TSTEP, coupled_lte=True)
        print('  returned a result after %d failures -- silently' % st['n'])
    except Exception as exc:                                   # noqa: BLE001
        print('  raised after %d attempts: %s: %s'
              % (st['n'], type(exc).__name__, exc))
        print('  -> not the RuntimeError at transient.py:508, and not a diagnostic:')
        print('     `_iq` is a side channel set inside solve_timestep, which never ran.')


def main():
    ok = midrun_failure()
    firststep_failure()
    print()
    print('Both paths reached from `coupled_lte=True`; neither exists in `solve()`.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
