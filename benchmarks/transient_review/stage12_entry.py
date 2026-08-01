"""STAGE 12 entry measurement -- does Fang's method have anything left to buy?

The stage-12 stub set an entry condition on ONE of Fang's two benefit mechanisms
and picked the wrong one.  It said:

    "Fang's benefit is the elimination of *rejected* steps [...] if rejections
     are still material -- say above 10% of accepted steps on the stress
     circuits -- then the bordered system is buying something a threshold
     cannot."

But the paper's own results (DAC 2013 sec. 4.1) attribute the headline number --
39% fewer time points, 17% wall clock -- to the *lower* bound of the two-sided
band eq (15), not to rejection elimination:

    "The introduction of the lower bound gamma_min prevents step sizes from
     being unnecessarily small."

So there are two questions, and this script measures both:

  1. REJECTIONS.  What fraction of steps are discarded and re-solved?  This is
     the stub's own entry condition, threshold 10%.

  2. WASTEFULLY SMALL STEPS.  Fang's gamma_min exists to stop the integrator
     taking steps whose LTE lands far below tolerance.  A step at err = 0.01
     could have been ~10x larger (err ~ h^p, p = 2 here).  If pycircuit's
     accepted steps already sit close to tolerance, gamma_min has nothing to
     tighten and the paper's mechanism is already present under another name.

Question 2 turns out to be the decisive one, and the answer is structural rather
than circuit-dependent: `IntegralController` predicts

    h_next = h * safety * (1/err)**(1/p)

which is a deadbeat law -- it drives the NEXT step's error to `safety**p`
exactly.  With safety = 0.9 and p = 2 that is 0.81.  The measurement below
recovers `p` from the realised step ratios rather than assuming it, so the
target is derived from the running code and not from this docstring.

Run:
    MPLBACKEND=Agg PYTHONPATH=.:benchmarks python benchmarks/transient_review/stage12_entry.py
"""
import warnings

import numpy as np

from pycircuit.circuit import numeric
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import stepcontroller as sc

import transient_stage4 as S

RELTOLS = (1e-4, 1e-6)
SAFETY = 0.9          # must match IntegralController's own constant


def _run(builder, reltol):
    """Solve once, returning (statistics, [(err, h_curr, h_next)] for accepts)."""
    cir, kw = builder()
    tran = Transient(cir, toolkit=numeric, reltol=reltol)

    rec = []
    orig = sc.IntegralController.evaluate_step

    def cap(self, *a, **k):
        h_curr = k['h_curr'] if 'h_curr' in k else a[5]
        accept, h_next = orig(self, *a, **k)
        ## `last_err` is None on the paths that accept without evaluating (the
        ## opening step); it is cleared on entry precisely so those cannot be
        ## mistaken for a real reading of zero error.
        if accept and self.last_err is not None:
            rec.append((float(self.last_err), float(h_curr), float(h_next)))
        return accept, h_next

    sc.IntegralController.evaluate_step = cap
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tran.solve(coupled_lte=False, **kw)
    finally:
        sc.IntegralController.evaluate_step = orig

    return tran.statistics, rec


def _recover_p(err, ratio):
    """Recover the controller's order exponent from realised step ratios.

    From h_next = h * safety * err**(-1/p), taking logs gives
    p = -log(err) / log(ratio/safety).  Restricted to steps that neither hit the
    growth cap nor sat at err == 1, where the law is not the binding constraint.
    """
    m = (ratio < 1.99) & (err > 1e-6)
    m &= np.abs(np.log(np.where(m, err, 1.0))) > 1e-3
    if not m.any():
        return float('nan')
    return float(np.median(-np.log(err[m]) / np.log(ratio[m] / SAFETY)))


def main():
    print(__doc__.split('Run:')[0].strip().splitlines()[0])
    print()
    print('1. REJECTIONS  (stub entry condition: material means >10%)')
    print()
    print('%-12s %8s %9s %9s %9s' % ('circuit', 'reltol', 'accepted',
                                     'rejected', 'reject %'))
    rej_max = 0.0
    results = {}
    for name, builder in S._fa_circuits():
        for reltol in RELTOLS:
            st, rec = _run(builder, reltol)
            results[(name, reltol)] = (st, rec)
            acc = st.accepted_steps
            ## Read directly, NOT via getattr with a default: a misremembered
            ## attribute name then reports a comfortable 0% instead of failing,
            ## and this measurement's whole job is to decide a stage on that
            ## number.  (It reported 0.0% everywhere once, from `n_rejected`.)
            rej = st.rejected_steps
            pct = 100.0 * rej / max(acc + rej, 1)
            rej_max = max(rej_max, pct)
            print('%-12s %8.0e %9d %9d %8.1f%%' % (name, reltol, acc, rej, pct))

    print()
    print('   worst rejection rate: %.1f%%  -> entry condition %s'
          % (rej_max, 'MET' if rej_max > 10.0 else 'NOT met'))

    print()
    print('2. WASTEFULLY SMALL STEPS  (what Fang\'s gamma_min actually targets)')
    print()
    print('%-12s %8s %6s %8s %8s %8s %9s' % ('circuit', 'reltol', 'p',
                                             'target', 'median', 'IQR',
                                             'within5%'))
    for name, builder in S._fa_circuits():
        for reltol in RELTOLS:
            _st, rec = results[(name, reltol)]
            err = np.array([r[0] for r in rec])
            ratio = np.array([r[2] / r[1] for r in rec])
            p = _recover_p(err, ratio)
            target = SAFETY ** round(p) if p == p else float('nan')
            q1, q3 = np.percentile(err, [25, 75])
            near = 100.0 * np.mean(np.abs(err - target) < 0.05 * target)
            print('%-12s %8.0e %6.2f %8.4f %8.4f %8.4f %8.1f%%'
                  % (name, reltol, p, target, np.median(err), q3 - q1, near))

    print()
    print('   "target" is safety**p, the fixed point of the controller\'s own')
    print('   prediction law -- steps clustered there are already as large as')
    print('   the LTE permits, which is the state gamma_min exists to reach.')


if __name__ == '__main__':
    main()
