"""STAGE 12B -- Fang's eq (6) LTE against pycircuit's charge-based one.

The stage-12 work was built against a charge divided-difference LTE; the paper's
eq (6) is a different estimator,

    eps_m = | v_i(t_m) - v_{i,extrapolated} |

measured in SOLUTION space.  Before anything is coupled to it, the two are
compared as drop-in step controllers on the same circuits at the same
tolerances, against CLOSED FORMS rather than a fine-mesh reference -- a
fine-mesh run would compare the integrator with itself and could not detect the
failure mode that matters here, a step count bought by integrating less
accurately.

What this is looking for, declared before running:

  1. eq (6) must produce a comparable step count at comparable accuracy.  If it
     needs far more steps for the same error it is not a viable replacement,
     however well conditioned it is.
  2. the controlling LTE node must actually move between time points -- the
     paper says it "may vary from time point to time point", and if it never
     moves here then `q^T` is effectively constant and the two-stage Newton of
     Fig. 4 is solving a problem we do not have.
  3. the conditioning claim, end to end: the normalised error must not run away
     as the step shrinks.

Run:
    MPLBACKEND=Agg PYTHONPATH=.:benchmarks python benchmarks/transient_review/stage12b_solution_lte.py
"""
import warnings

import numpy as np

from pycircuit.circuit import numeric
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import stepcontroller as sc

import transient_stage4 as S

RELTOLS = (1e-4, 1e-5, 1e-6)


def _rc_vsin_analytic(t, r=1e3, c=1e-7, va=1.0, freq=1e3):
    tau = r * c
    w = 2.0 * np.pi * freq
    A = va / np.sqrt(1.0 + (w * tau) ** 2)
    phi = np.arctan(w * tau)
    return A * (np.sin(w * t - phi) + np.sin(phi) * np.exp(-t / tau))


def _stiff_rlc_analytic(t, r=1.0, l=1e-6, c=1e-6):
    alpha = r / (2.0 * l)
    wd = np.sqrt(1.0 / (l * c) - alpha ** 2)
    return np.exp(-alpha * t) * (np.cos(wd * t) + (alpha / wd) * np.sin(wd * t))


REFERENCE = {'rc-vsin': ('b', _rc_vsin_analytic),
             'stiff-rlc': ('1', _stiff_rlc_analytic)}


def _run(builder, reltol, controller_cls, reference):
    cir, kw = builder()
    tran = Transient(cir, toolkit=numeric, reltol=reltol)
    if controller_cls is not None:
        tran.step_controller = controller_cls()

    errs, nodes = [], []
    orig = controller_cls.evaluate_step if controller_cls \
        else sc.IntegralController.evaluate_step
    target = controller_cls or sc.IntegralController

    def cap(self, *a, **k):
        accept, h_next = orig(self, *a, **k)
        if accept and self.last_err is not None:
            errs.append(float(self.last_err))
            ci = getattr(self, 'controlling_index', None)
            if ci is not None:
                nodes.append(int(ci))
        return accept, h_next

    target.evaluate_step = cap
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(coupled_lte=False, **kw)
    finally:
        target.evaluate_step = orig

    w = res.v(reference[0])
    t = np.asarray(w.x, dtype=float).ravel()
    v = np.asarray(w.y, dtype=float).ravel()
    d = np.abs(v - reference[1](t))
    ## Past the opening step, which is accepted unevaluated and which no
    ## tolerance governs -- on the stiff circuit it otherwise dominates.
    post = float(np.max(d[2:])) if len(d) > 2 else float(np.max(d))

    return dict(steps=tran.statistics.accepted_steps,
                rejected=tran.statistics.rejected_steps,
                err=post, errs=np.array(errs), nodes=nodes,
                done=float(t[-1]) >= kw['tend'] * (1 - 1e-6))


def main():
    circuits = dict(S._fa_circuits())

    print('eq (6) solution-space LTE  vs  charge divided difference')
    print('(accuracy against closed forms; "err" excludes the opening step)')
    print()
    for cname, ref in REFERENCE.items():
        print('--- %s' % cname)
        print('%-8s  %-22s %7s %7s %12s %9s  %s'
              % ('reltol', 'estimator', 'steps', 'reject', 'max|err|',
                 'med err', 'done'))
        for reltol in RELTOLS:
            for label, cls in (('charge (current)', None),
                               ('eq (6) solution', sc.SolutionLTEController)):
                try:
                    r = _run(circuits[cname], reltol, cls, ref)
                except Exception as e:
                    print('%-8.0e  %-22s  raised %s: %s'
                          % (reltol, label, type(e).__name__, str(e)[:40]))
                    continue
                med = np.median(r['errs']) if len(r['errs']) else float('nan')
                print('%-8.0e  %-22s %7d %7d %12.3e %9.4f  %s'
                      % (reltol, label, r['steps'], r['rejected'], r['err'],
                         med, 'yes' if r['done'] else 'NO'))
        print()

    print('--- does the controlling LTE node move?  (Fang: it "may vary from')
    print('    time point to time point"; if it never moves, q^T is constant')
    print('    and Fig. 4\'s two-stage split is solving a problem we do not have)')
    print()
    print('%-12s %-8s %10s %10s  %s'
          % ('circuit', 'reltol', 'distinct', 'switches', 'share of steps'))
    for cname, ref in REFERENCE.items():
        for reltol in RELTOLS:
            r = _run(circuits[cname], reltol, sc.SolutionLTEController, ref)
            nodes = r['nodes']
            if not nodes:
                print('%-12s %-8.0e   (no controlling index recorded)'
                      % (cname, reltol))
                continue
            switches = sum(1 for a, b in zip(nodes, nodes[1:]) if a != b)
            print('%-12s %-8.0e %10d %10d  %.1f%%'
                  % (cname, reltol, len(set(nodes)), switches,
                     100.0 * switches / max(1, len(nodes) - 1)))


if __name__ == '__main__':
    main()
