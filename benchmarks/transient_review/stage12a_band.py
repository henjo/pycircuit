"""STAGE 12A gates -- Fang's acceptance band (eq 15) and step damper (eq 16).

Gate 12A-1: does the band change step placement at all?
Gate 12A-2: does any step it saves come out of accuracy rather than out of waste?

The second gate is the one that matters, and it needs an accuracy reference that
is understood rather than approximated.  `rc-vsin` is a linear RC driven by a
sine, so it has a closed form and no reference run is needed:

    v_b(t) = A [ sin(wt - phi) + sin(phi) exp(-t/tau) ],
    tau = RC,  A = 1/sqrt(1 + (w tau)^2),  phi = atan(w tau)

(the particular solution plus the homogeneous term that puts v_b(0) at zero).
Comparing against a fine-mesh run instead would compare the integrator with
itself, which cannot detect the failure mode this gate exists to catch -- a
step count bought by quietly integrating less accurately.

THE HYPOTHESIS UNDER TEST, stated before running (stage-12 entry measurement):
`IntegralController` is deadbeat, driving the accepted error to the fixed point
`safety**p = 0.81`.  For a controller already sitting on its aim point, raising
`gamma_min` above 0.81 does not recover waste -- there is none -- it just moves
the aim point up, which is arithmetically a tolerance change.  So the `gmin=0.95`
configuration below is paired with a `reltol` scaled by 0.95/0.81, and if the two
agree in BOTH step count and error, the lower bound is a re-parameterisation of
the tolerance on this controller and not an independent mechanism.

Run:
    MPLBACKEND=Agg PYTHONPATH=.:benchmarks python benchmarks/transient_review/stage12a_band.py
"""
import warnings

import numpy as np

from pycircuit.circuit import numeric, gnd
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import stepcontroller as sc

import transient_stage4 as S

## The aim point of the unbanded controller: safety**p with safety=0.9, p=2.
## Measured, not assumed -- benchmarks/transient_review/stage12_entry.py recovers
## p from the realised step ratios and gets 2.00.
DEADBEAT_TARGET = 0.81
GMIN_HIGH = 0.95

RELTOL = 1e-5

## (label, reltol scale, band kwargs)
CONFIGS = [
    ('baseline',            1.0,                        {}),
    ('paper 0.7/3.0',       1.0,                        dict(gamma_min=0.7, gamma_max=3.0)),
    ('gmax=3 only',         1.0,                        dict(gamma_min=0.0, gamma_max=3.0)),
    ('gmin=0.95',           1.0,                        dict(gamma_min=GMIN_HIGH, gamma_max=1.0)),
    ('reltol x %.3f' % (GMIN_HIGH / DEADBEAT_TARGET),
                            GMIN_HIGH / DEADBEAT_TARGET, {}),
    ('eta=0.15 damper',     1.0,                        dict(eta=0.15)),
]


def _rc_vsin_analytic(t, r=1e3, c=1e-7, va=1.0, freq=1e3):
    """Closed form for the `rc-vsin` output node, zero initial charge."""
    tau = r * c
    w = 2.0 * np.pi * freq
    A = va / np.sqrt(1.0 + (w * tau) ** 2)
    phi = np.arctan(w * tau)
    return A * (np.sin(w * t - phi) + np.sin(phi) * np.exp(-t / tau))


def _stiff_rlc_analytic(t, r=1.0, l=1e-6, c=1e-6):
    """Closed form for `stiff-rlc`'s node 1.

    ADDED AFTER the eta=0.15 damper crossed this circuit's whole ringing
    transient in 62 steps against the baseline's 490 and was visible only as a
    step count, because accuracy was being measured on `rc-vsin` alone.  A
    configuration that integrates a stiff ringdown wrongly must fail an ACCURACY
    gate, not merely look odd in a step column -- so the circuit where the defect
    appeared gets its own reference.

    C1 sits from node 1 to ground, with R1 in series with L1 as the other branch
    across it, so with `i` the branch current:

        C dv/dt = -i,        L di/dt = v - R i     =>    LC v'' + RC v' + v = 0

    `x0` sets both node voltages to 1 and leaves the inductor current at zero, so
    v(0) = 1 and v'(0) = -i(0)/C = 0 (consistent with v2(0) = v1 - R i = 1).
    """
    alpha = r / (2.0 * l)
    w0_sq = 1.0 / (l * c)
    wd = np.sqrt(w0_sq - alpha ** 2)          # underdamped for these values
    return np.exp(-alpha * t) * (np.cos(wd * t) + (alpha / wd) * np.sin(wd * t))


## Which node carries the closed form, per circuit.
REFERENCE = {'rc-vsin': ('b', _rc_vsin_analytic),
             'stiff-rlc': ('1', _stiff_rlc_analytic)}


def _run(builder, band, reltol, reference=None):
    cir, kw = builder()
    tran = Transient(cir, toolkit=numeric, reltol=reltol,
                     lte_gamma_min=band.get('gamma_min', 0.0),
                     lte_gamma_max=band.get('gamma_max', 1.0),
                     lte_eta=band.get('eta', None))

    errs = []
    orig = sc.IntegralController.evaluate_step

    def cap(self, *a, **k):
        accept, h_next = orig(self, *a, **k)
        if accept and self.last_err is not None:
            errs.append(float(self.last_err))
        return accept, h_next

    sc.IntegralController.evaluate_step = cap
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(coupled_lte=False, **kw)
    finally:
        sc.IntegralController.evaluate_step = orig

    ## Did the run actually finish?  A configuration that stops early reports a
    ## flattering step count for the worst possible reason, and "fewer steps" is
    ## precisely the number these gates turn on -- so the end time is checked
    ## rather than assumed.  (`eta=0.15` on the stiff ringdown reported -87.3%
    ## steps in the first run of this benchmark, which is what prompted this.)
    probe = reference[0] if reference is not None \
        else next(nd for nd in cir.nodes if nd.name != 'gnd')
    w = res.v(probe)
    t = np.asarray(w.x, dtype=float).ravel()
    t_reached = float(t[-1])
    complete = t_reached >= kw['tend'] * (1.0 - 1e-6)

    max_err = post_err = float('nan')
    if reference is not None:
        v = np.asarray(w.y, dtype=float).ravel()
        d = np.abs(v - reference[1](t))
        max_err = float(np.max(d))
        ## The opening step is accepted UNEVALUATED -- there is no past point to
        ## difference, so no tolerance governs it (see `IntegralController`'s
        ## `no_history` branch).  On the stiff ringdown that single step dominates
        ## the maximum: the error saturates at 1.3589e-02 at t=2e-7 for every
        ## reltol from 1e-5 to 1e-8, i.e. four decades of tolerance change move it
        ## not at all, which makes the all-points maximum blind to exactly what
        ## this gate is trying to compare.  Both are reported: the total is the
        ## honest accuracy, the second column is the part step control owns.
        post = d[2:] if len(d) > 2 else d
        post_err = float(np.max(post))

    return tran.statistics, np.array(errs), max_err, post_err, complete, t_reached


def main():
    circuits = dict(S._fa_circuits())

    print('GATE 12A-1 -- does the band move the accepted-error distribution?')
    print('(reltol=%g; "median err" is the normalised LTE of accepted steps;' % RELTOL)
    print(' the unbanded controller aims at %.2f)' % DEADBEAT_TARGET)
    print()
    for cname in ('rc-vsin', 'stiff-rlc', 'rc-pulse'):
        print('--- %s' % cname)
        print('%-18s %8s %8s %10s %10s %10s  %s'
              % ('config', 'accept', 'reject', 'median', 'IQR', 'vs base', 'done'))
        base = None
        for label, scale, band in CONFIGS:
            st, errs, _, _, complete, t_end = _run(circuits[cname], band,
                                                   RELTOL * scale)
            q1, q3 = np.percentile(errs, [25, 75])
            if base is None:
                base = st.accepted_steps
            rel = 100.0 * (st.accepted_steps - base) / base
            print('%-18s %8d %8d %10.4f %10.4f %+9.1f%%  %s'
                  % (label, st.accepted_steps, st.rejected_steps,
                     np.median(errs), q3 - q1, rel,
                     'yes' if complete else 'NO (t=%.3g)' % t_end))
        print()

    print('GATE 12A-2 -- is a saved step paid for out of accuracy?')
    print('(against closed forms -- no reference run, so the integrator is not')
    print(' being compared with itself)')
    for cname, ref in REFERENCE.items():
        print()
        print('--- %s, node %s' % (cname, ref[0]))
        print('%-18s %8s %8s %12s %12s %9s  %s'
              % ('config', 'accept', 'reject', 'max|err|', 'after open',
                 'vs base', 'done'))
        base_err = None
        for label, scale, band in CONFIGS:
            st, _errs, max_err, post_err, complete, t_end = _run(
                circuits[cname], band, RELTOL * scale, reference=ref)
            if base_err is None:
                base_err = post_err
            print('%-18s %8d %8d %12.3e %12.3e %8.2fx  %s'
                  % (label, st.accepted_steps, st.rejected_steps, max_err,
                     post_err, post_err / base_err,
                     'yes' if complete else 'NO (t=%.3g)' % t_end))

    print()
    print('Read `gmin=0.95` against the scaled-reltol row: if they agree in both')
    print('step count and error, the lower bound is a re-parameterisation of the')
    print('tolerance on a deadbeat controller, not an independent mechanism --')
    print('and the rejection columns then say which parameterisation is cheaper.')


if __name__ == '__main__':
    main()
