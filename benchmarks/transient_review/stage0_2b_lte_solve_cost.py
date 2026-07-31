"""STAGE 0.2b -- what does the LTE's J^-1 solve actually cost per step?

`doc/transient_work_plan.md` 0.2b asks for "the per-step cost of the extra `G` +
solve as a fraction of the whole step, on a circuit of realistic size", with the
decision rule: under 10%, delete the JAX `estimate_lte` charge path and make
`lte_formula` select an algebra only; over 10%, fix its tolerance and keep both.

FIRST, A CORRECTION TO THE QUESTION.  There is no extra `G` evaluation to measure.
Both estimators receive the Jacobian as an argument from the Newton solve that just
converged -- CPU `IntegralController.evaluate_step(..., J, ...)` and JAX
`ywr_error_ratio(i_curr, x_curr, x_last, J, ...)`.  Neither calls `cir.G`.  So the
charge-domain path's advantage was never "a `G` evaluation and a linear solve"; it is
a linear solve and two `np.delete` copies, and half the claimed advantage does not
exist.  What follows measures what is actually there.

WHAT IS TIMED, on the 136-unknown leapfrog:

  1. the whole transient
  2. `IntegralController.evaluate_step` in total  -- the entire step-control cost
  3. inside it, `remove_row_col` + `toolkit.linearsolver` -- the part a charge-domain
     estimator would avoid

(3)/(1) is the number the decision rule wants.  (2)-(3) is what a charge-domain path
would still pay, and it is reported so the saving is not overstated: switching
estimators removes (3), not (2).

Timing is by wrapping the methods, so it includes Python call overhead on the
measured side -- which biases the result *towards* "the solve is expensive" and
therefore towards keeping both paths.  A result under 10% despite that bias is a
safe result; one just over 10% would not be.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks MPLBACKEND=Agg python3 -u \
        benchmarks/transient_review/stage0_2b_lte_solve_cost.py
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
from pycircuit.circuit import stepcontroller as sc_module
from pycircuit.circuit import analysis as analysis_module

import nonlinear_leapfrog_sweep as H

WINDOW = 2.5e-6
AMPLITUDE = 1.0


class Timers(object):
    def __init__(self):
        self.eval_step = 0.0
        self.eval_calls = 0
        self.solve = 0.0
        self.solve_calls = 0
        self.removerc = 0.0
        self.removerc_calls = 0


def instrument(timers):
    """Wrap evaluate_step, linearsolver and remove_row_col.  Returns an undo fn."""
    real_eval = sc_module.IntegralController.evaluate_step
    real_rrc = analysis_module.remove_row_col
    real_solve = numeric.linearsolver

    def eval_step(self, *a, **kw):
        t0 = time.perf_counter()
        try:
            return real_eval(self, *a, **kw)
        finally:
            timers.eval_step += time.perf_counter() - t0
            timers.eval_calls += 1

    def rrc(*a, **kw):
        t0 = time.perf_counter()
        try:
            return real_rrc(*a, **kw)
        finally:
            timers.removerc += time.perf_counter() - t0
            timers.removerc_calls += 1

    def solve(*a, **kw):
        t0 = time.perf_counter()
        try:
            return real_solve(*a, **kw)
        finally:
            timers.solve += time.perf_counter() - t0
            timers.solve_calls += 1

    sc_module.IntegralController.evaluate_step = eval_step
    analysis_module.remove_row_col = rrc
    numeric.linearsolver = solve

    def undo():
        sc_module.IntegralController.evaluate_step = real_eval
        analysis_module.remove_row_col = real_rrc
        numeric.linearsolver = real_solve
    return undo


def measure(label, integrator_factory):
    cir, amp = H.build_transient_leapfrog(AMPLITUDE)
    out_node = amp[H.STAGES - 1]['out']
    opts = dict(H.TRAN_OPTS)
    opts['integrator'] = integrator_factory()

    ## Baseline: the same run, uninstrumented, so the wrapper overhead is visible
    ## rather than assumed negligible.
    t0 = time.perf_counter()
    res = Transient(cir, toolkit=numeric, **opts).solve(
        refnode=gnd, tend=WINDOW, timestep=1.0 / (H.F1 * H.POINTS))
    plain = time.perf_counter() - t0
    steps = len(np.asarray(res.v(out_node).x[0]))

    cir, amp = H.build_transient_leapfrog(AMPLITUDE)
    opts = dict(H.TRAN_OPTS)
    opts['integrator'] = integrator_factory()
    timers = Timers()
    undo = instrument(timers)
    try:
        t0 = time.perf_counter()
        Transient(cir, toolkit=numeric, **opts).solve(
            refnode=gnd, tend=WINDOW, timestep=1.0 / (H.F1 * H.POINTS))
        total = time.perf_counter() - t0
    finally:
        undo()

    ## The controller calls remove_row_col and linearsolver exactly once each per
    ## evaluate_step, but the Newton solve calls linearsolver many times more.  The
    ## LTE share is therefore the per-call solve cost times the controller's call
    ## count, not the whole linearsolver total.
    per_solve = timers.solve / max(timers.solve_calls, 1)
    per_rrc = timers.removerc / max(timers.removerc_calls, 1)
    lte_solve_share = per_solve * timers.eval_calls
    lte_rrc_share = per_rrc * timers.eval_calls
    lte_linear = lte_solve_share + lte_rrc_share

    print('%s' % label)
    print('  steps                       %6d' % steps)
    print('  transient, uninstrumented   %8.3f s' % plain)
    print('  transient, instrumented     %8.3f s   (wrapper overhead %+.1f%%)'
          % (total, 100.0 * (total - plain) / plain))
    print('  evaluate_step total         %8.3f s   %5.2f%% of the run   (%d calls)'
          % (timers.eval_step, 100.0 * timers.eval_step / total, timers.eval_calls))
    print('  all linearsolver calls      %8.3f s   %5.2f%% of the run   (%d calls,'
          ' %.1f per LTE evaluation)'
          % (timers.solve, 100.0 * timers.solve / total, timers.solve_calls,
             timers.solve_calls / max(timers.eval_calls, 1)))
    print('  ---- the part a charge-domain estimator would avoid ----')
    print('  LTE J^-1 solve              %8.3f s   %5.2f%% of the run'
          % (lte_solve_share, 100.0 * lte_solve_share / total))
    print('  LTE remove_row_col copies   %8.3f s   %5.2f%% of the run'
          % (lte_rrc_share, 100.0 * lte_rrc_share / total))
    print('  TOTAL AVOIDABLE             %8.3f s   %5.2f%% of the run'
          % (lte_linear, 100.0 * lte_linear / total))
    print('  still paid either way       %8.3f s   %5.2f%% of the run'
          % (timers.eval_step - lte_linear,
             100.0 * (timers.eval_step - lte_linear) / total))
    print()
    return 100.0 * lte_linear / total


def main():
    print('STAGE 0.2b -- cost of the LTE J^-1 solve, %d-unknown leapfrog, %.3g s window'
          % (H.build_transient_leapfrog(AMPLITUDE)[0].n, WINDOW))
    print('decision rule: under 10% -> the charge-domain path buys too little to keep')
    print()

    shares = []
    for label, factory in (
            ('TRAPEZOIDAL (the harness default)', lambda: TrapezoidalIntegrator()),
            ('GEAR2 (the shipped default)', lambda: Gear2Integrator()),
    ):
        shares.append(measure(label, factory))

    worst = max(shares)
    print('GATE 0.2b')
    print('  largest avoidable share across the configurations: %.2f%%' % worst)
    print('  decision rule threshold: 10%%')
    print('  -> %s' % ('OVER: keep both paths and fix the charge tolerance'
                       if worst >= 10.0 else
                       'UNDER: the charge-domain shortcut does not pay for itself'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
