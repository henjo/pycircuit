"""Adaptive stepping for the Nordsieck GLM: cost, restarts, and the order on a
grid the controller did not choose.

Two things worth running here:

* `cost()` -- device evaluations, accepted/rejected steps and STARTUP count
  against the Runge-Kutta methods.  The startup count is the interesting
  column: a GLM carries a Nordsieck vector between steps, and a REJECTED step
  has already overwritten it, so without a second slot holding the vector the
  step CONSUMED, every retry pays a full startup and the run livelocks into
  rejecting almost every step.  Measured before that slot existed: GLM2 at
  2732 accepted / 2732 rejected / 2733 startups, 102034 device evaluations
  against radau's 2419.

* `order_on_a_jittered_grid()` -- Thm 9.5 is stated at CONSTANT stepsize, so
  the variable-step Nordsieck rescale `Q_k <- (h_new/h_old)^k Q_k` is outside
  the result the method rests on.  This refines a smoothly non-uniform grid
  (h varying ~3x) and reports the observed order, with the uniform grid as its
  own control.
"""
import warnings

import numpy as np

from pycircuit.circuit.circuit import gnd
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.integrator import (GLM2Integrator, GLM3Integrator,
                                          GLM4Integrator, RadauIIA3Integrator,
                                          TRBDF2Integrator)

warnings.simplefilter('ignore')
PER = 1e-3


def _expg():
    from pycircuit.circuit.tests.test_stage_predictor import _expg_fixture
    return _expg_fixture(PER)


def _cv():
    from pycircuit.circuit.tests.test_glm import _cv_loop
    return _cv_loop(PER)


def cost():
    tr = Transient(_expg(), integrator=RadauIIA3Integrator(), reltol=1e-13)
    rw = tr.solve(refnode=gnd, tend=PER, timestep=PER / 8000,
                  fixed_timestep=True).v('b')
    t_ref = np.asarray(rw.x[0], dtype=float)
    y_ref = np.asarray(rw.y, dtype=float)
    print('%-7s %-8s %8s %7s %7s %9s %11s'
          % ('method', 'reltol', 'i-calls', 'steps', 'reject', 'startups',
             'err vs ref'))
    for cls, name in ((GLM2Integrator, 'glm2'), (GLM3Integrator, 'glm3'),
                      (GLM4Integrator, 'glm4'),
                      (RadauIIA3Integrator, 'radau'),
                      (TRBDF2Integrator, 'trbdf2')):
        for reltol in (1e-6, 1e-9):
            n = {'i': 0}
            cir = _expg()
            fi = cir.i
            cir.i = lambda *a, **k: (n.__setitem__('i', n['i'] + 1),
                                     fi(*a, **k))[1]
            tr = Transient(cir, integrator=cls(), reltol=reltol)
            res = tr.solve(refnode=gnd, tend=PER, timestep=PER / 200,
                           fixed_timestep=False)
            w = res.v('b')
            t = np.asarray(w.x[0], dtype=float)
            y = np.asarray(w.y, dtype=float)
            err = float(np.max(np.abs(y - np.interp(t, t_ref, y_ref))))
            print('%-7s %-8.0e %8d %7d %7d %9d %11.3e'
                  % (name, reltol, n['i'], tr.statistics.accepted_steps,
                     tr.statistics.rejected_steps,
                     getattr(tr, 'statistics_glm_startups', 0), err))
    print()


def _march(cls, times, build):
    cir = build()
    tr = Transient(cir, integrator=cls(), reltol=1e-13)
    tr.irefnode = cir.get_node_index(gnd)
    x = np.asarray(DC(cir, refnode=gnd).solve().x, dtype=float).ravel()
    tr.epar.t = 0.0
    tr._begin_run(x, cir.n)
    for j in range(1, len(times)):
        tr._dt_last = tr._dt if j > 1 else None
        tr._dt = times[j] - times[j - 1]
        tr.epar.t = times[j]
        x, _f, _J, _ = tr.solve_timestep(x, times[j])
        tr._push_history(x)
    return np.asarray(x, dtype=float).ravel()


def _grid(npts, jitter):
    u = np.linspace(0.0, 1.0, npts + 1)
    if jitter:
        u = u + 0.25 * np.sin(2 * np.pi * u) / np.pi
        u = (u - u[0]) / (u[-1] - u[0])
    return PER * u


def order_on_a_jittered_grid():
    ref = _march(RadauIIA3Integrator, _grid(12000, False), _cv)
    print('%-6s %-9s %-38s %s' % ('method', 'grid', 'endpoint error', 'slopes'))
    for cls, name, p in ((GLM2Integrator, 'glm2', 2),
                         (GLM3Integrator, 'glm3', 3),
                         (GLM4Integrator, 'glm4', 4),
                         (RadauIIA3Integrator, 'radau', 5)):
        for jit in (False, True):
            errs = [float(np.max(np.abs(_march(cls, _grid(n, jit), _cv) - ref)))
                    for n in (100, 200, 400, 800)]
            sl = [np.log(errs[i - 1] / errs[i]) / np.log(2)
                  for i in range(1, len(errs))]
            print('%-6s %-9s %-38s %s  (uniform order %d)'
                  % (name, 'JITTERED' if jit else 'uniform',
                     ' '.join('%.2e' % x for x in errs),
                     ' '.join('%5.2f' % x for x in sl), p))
    print()


def estimator_order():
    """⚠ The FILTERED estimate reads one order below `EMBEDDED_ORDER + 1` on a
    DAE -- for every method here, not just the GLM.  Run this before believing
    any calibration claim about the controller."""
    from pycircuit.circuit.integrator import ESDIRK43Integrator
    print('%-9s %-8s %-34s %s'
          % ('method', 'EMB_ORD', 'estimate (median)', 'slopes'))
    for cls, name in ((RadauIIA3Integrator, 'radau'),
                      (TRBDF2Integrator, 'trbdf2'),
                      (ESDIRK43Integrator, 'esdirk43'),
                      (GLM2Integrator, 'glm2'), (GLM3Integrator, 'glm3'),
                      (GLM4Integrator, 'glm4')):
        v = []
        for N in (50, 100, 200, 400):
            tr = Transient(_expg(), integrator=cls(), reltol=1e-12)
            tr._rk_want_est = True
            rec = []
            orig = Transient.solve_timestep

            def w(self, x0, t, *a, **k):
                r = orig(self, x0, t, *a, **k)
                e = getattr(self, '_rk_est', None)
                if e is not None:
                    rec.append(float(np.max(np.abs(np.asarray(e, float)))))
                return r
            Transient.solve_timestep = w
            try:
                tr.solve(refnode=gnd, tend=PER, timestep=PER / N,
                         fixed_timestep=True)
            finally:
                Transient.solve_timestep = orig
            v.append(float(np.median(np.array(rec)[6:])))
        sl = [np.log(v[i - 1] / v[i]) / np.log(2) for i in range(1, len(v))]
        print('%-9s %-8d %-34s %s'
              % (name, cls().EMBEDDED_ORDER, ' '.join('%.2e' % x for x in v),
                 ' '.join('%5.2f' % x for x in sl)))
    print()


if __name__ == '__main__':
    cost()
    order_on_a_jittered_grid()
    estimator_order()
