"""The measurements that CHANGED A DECISION or REORDERED THE PLAN.

`benchmarks/transient_review/README.md` records that the numerical lens's probe
scripts "were not recovered before the scratchpad was reused", so its results are
reproducible only by reconstruction. This file exists so that does not happen again
to the 2026-07-31 work: every measurement here is one whose *outcome* is quoted in
`doc/transient_work_plan.md` as the reason something was reversed, reordered or
declined. A quoted number whose script is gone is a number nobody can argue with.

Not included: the routine gate measurements, which already live in
`transient_stage2.py`, `transient_stage3.py` and `transient_stage4.py`.

    --parity            why 4d was thought to be blocked on 4g  (SEE THE NOTE BELOW)
    --relref            why decision D3 was reverted the first time
    --vabstol           why 0.3a's tolerance split cost nothing
    --matched-accuracy  why "1.7-2.5x fewer steps" was struck to 1.31-2.06x  (D3-c)
    --lte-floor         why `lte_vabstol` went back to 1e-12                 (D3-e)
    --gear2-accuracy    what 4i's extra steps bought, after a gate clause FAILED
    --pi-gains          what fixing the PI gains bought end to end            (4a-4)

**`--parity` is kept as a record of a MISTAKE, not of a result.** Its table is
measured against the *propagated* companion-current error, which is O(h^3) where the
local truncation error is O(h^2) -- so it divides by a quantity an order smaller that
is passing through a zero crossing, and the "48x parity swing" it reports is an
artefact. Against the local error trapezoidal is flat. See gate 4g-b0. The script
stays because the plan quotes its numbers as the reason 4d was reordered, and a
reader chasing that needs to be able to see what went wrong.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks MPLBACKEND=Agg python3 -u \\
          benchmarks/transient_decisions.py [--<one of the above>|--all]
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
    print()
    print('*** THE CONCLUSION ABOVE IS WRONG, AND THIS MODE IS KEPT TO SHOW WHY. ***')
    print()
    print('The quantity differenced here is the PROPAGATED companion-current error,')
    print('which is O(h^3) where the local truncation error -- the part a step')
    print('controller must estimate -- is O(h^2).  So the table divides by something an')
    print('order SMALLER that is passing through a zero crossing, and the 48x swing is')
    print('an artefact of the denominator, not a property of the method.  The tell is')
    print('that it changes sign at h=1e-11; a diverging quantity does not.')
    print()
    print('Measured against the LOCAL truncation error, trapezoidal reads -2.8395e+01 at')
    print('every prefix length from 3 to 7 -- no parity dependence at all.  See gate')
    print('4g-b0.  So "4d is blocked on 4g" was never true: 4d was blocked on a')
    print('contaminated reference.  Trap 6 in the plan is the general lesson --')
    print('before believing a ratio, look at numerator and denominator separately.')


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


def matched_accuracy():
    """WHY THE "1.7-2.5x FEWER STEPS" CLAIM FOR `sigglobal` WAS STRUCK (gate D3-c).

    The recorded justification for `relref='sigglobal'` compared step counts at equal
    `reltol`.  The two modes do not deliver equal ACCURACY at equal `reltol` --
    `sigglobal`'s error is ~1.5x larger there -- so a step count taken at a looser
    effective tolerance is a relabelling of the tolerance, not a speedup.

    This sweeps `reltol` finely in each mode and interpolates the steps each needs to
    reach the SAME global error, inside the error range both modes cover so neither is
    extrapolated.  The gain is real, but it is 1.31-2.06x and not 1.7-2.5x.

    Recorded outcome: worst-case ratio per integrator 1.48x (euler), 1.31x
    (trap-classic), 1.44x (both gear2), against a declared threshold of 1.2x.
    """
    sys.path.insert(0, os.path.join(_REPO, 'pycircuit', 'circuit', 'tests'))
    import warnings
    import test_analysis_transient as T
    from pycircuit.circuit.transient import Transient

    reltols = [10 ** (-e) for e in (2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5)]

    def _set(pname, value):
        for par in Transient.parameters:
            if getattr(par, 'name', '') == pname:
                par.default = value

    def sweep(name, mode, floor):
        ## `lte_vabstol` MUST be pinned, and this is not a detail.  Gate D3-c ran
        ## while the floor was 1e-6; decision D3-e then moved the shipped default to
        ## 1e-12.  The floor is load-bearing under `pointlocal` and inert under
        ## `sigglobal` (that is D3-e's whole finding), so leaving it at the default
        ## penalises the `pointlocal` arm by 8.5-9.2% and inflates the ratio this
        ## function exists to report -- from 1.31-2.06x to 1.70-2.94x.  Isolating
        ## `relref` means holding the floor fixed across both arms.
        _set('relref', mode)
        _set('lte_vabstol', floor)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            return [T._stiff_run(name, rt)[::2] for rt in reltols]

    def steps_at(rows, target):
        """log-log interpolation of step count against global error."""
        xs = np.log([r[1] for r in rows])
        ys = np.log([r[0] for r in rows])
        o = np.argsort(xs)
        if not (xs[o][0] <= np.log(target) <= xs[o][-1]):
            return None
        return float(np.exp(np.interp(np.log(target), xs[o], ys[o])))

    for floor, label in ((1e-6, "the floor gate D3-c ran under"),
                         (1e-12, "the floor shipped since D3-e")):
        print('--- lte_vabstol = %g   (%s)' % (floor, label))
        print('  %-16s %26s' % ('integrator',
                                'matched-accuracy ratio, over 5 error targets'))
        lo_all, hi_all = 1e9, 0.0
        for name in ('euler-classic', 'trap-classic', 'gear2-classic', 'gear2-ywr'):
            pl = sweep(name, 'pointlocal', floor)
            sg = sweep(name, 'sigglobal', floor)
            lo = max(min(r[1] for r in pl), min(r[1] for r in sg))
            hi = min(max(r[1] for r in pl), max(r[1] for r in sg))
            ratios = []
            for target in np.exp(np.linspace(np.log(lo * 1.05), np.log(hi * 0.95), 5)):
                a, b = steps_at(pl, target), steps_at(sg, target)
                if a and b:
                    ratios.append(a / b)
            lo_all, hi_all = min(lo_all, min(ratios)), max(hi_all, max(ratios))
            print('  %-16s %8.2fx .. %.2fx   %s'
                  % (name, min(ratios), max(ratios),
                     ' '.join('%.2f' % v for v in ratios)))
        print('  RANGE OVER ALL INTEGRATORS AND TARGETS: %.2fx .. %.2fx'
              % (lo_all, hi_all))
        print()
    ## Leave the defaults as this file found them.
    _set('relref', 'sigglobal')
    _set('lte_vabstol', 1e-12)
    print('At the gate\'s own floor the range is 1.31-2.06x, which is what the plan')
    print('records and what struck the "1.7-2.5x" taken at equal reltol.  At the')
    print('shipped floor it reads larger -- but that number is `relref` AND `lte_vabstol`')
    print('together, since the floor is load-bearing under `pointlocal` only.  Either')
    print('way the gain clears the 1.2x the gate declared.')


def lte_floor():
    """WHY `lte_vabstol` WENT BACK TO 1e-12 (gate D3-e), and what it used to buy.

    1e-6 was a workaround: under `pointlocal` the relative reference degenerates to
    zero on a node carrying no signal, so the absolute floor is what stops the
    tolerance collapsing there.  `sigglobal` removes the cause, and the floor stops
    being load-bearing -- so the principled value (matching `vabstol`) can come back.

    The `pointlocal` rows are the control: they are what "load-bearing" looks like.

    Recorded outcome: bit-identical under `sigglobal` at 1e-6 / 1e-9 / 1e-12 on both
    circuits (403 and 601 steps); +8.5% and +9.2% under `pointlocal`.
    """
    import warnings
    from pycircuit.circuit import gnd, numeric
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import R, C, VS, VPulse
    from pycircuit.circuit.transient import Transient

    def healthy(floor, mode):
        c = SubCircuit(toolkit=numeric)
        c['VS'] = VPulse('n1', gnd, v1=0, v2=5, td=1e-7, tr=1e-8, tf=1e-8,
                         pw=5e-7, per=1e-6)
        c['R'] = R('n1', 'n2', r=1e3)
        c['C'] = C('n2', gnd, c=1e-10)
        t = Transient(c, toolkit=numeric, lte_vabstol=floor, relref=mode)
        return t.solve(refnode=gnd, tend=3e-6, timestep=1e-7), 'n2'

    def quiet(floor, mode):
        """A node driven only through 1 Gohm -- |x| there is ~0, which is what makes
        `pointlocal`'s reference degenerate."""
        c = SubCircuit(toolkit=numeric)
        c['VS'] = VS('n1', gnd, v=1.0)
        c['R'] = R('n1', 'n2', r=1e3)
        c['C'] = C('n2', gnd, c=1e-9)
        c['Rq'] = R('n3', gnd, r=1e9)
        c['Rq2'] = R('n3', 'n2', r=1e9)
        t = Transient(c, toolkit=numeric, relref=mode, lte_vabstol=floor,
                      reltol=1e-6, uic=True)
        return t.solve(refnode=gnd, tend=20e-6, timestep=4e-6), 'n2'

    print('%-10s %-12s %13s %8s %14s %9s'
          % ('circuit', 'relref', 'lte_vabstol', 'steps', 'v_end', 'vs 1e-6'))
    for label, fn in (('healthy', healthy), ('quiet', quiet)):
        for mode in ('sigglobal', 'pointlocal'):
            base = None
            for floor in (1e-6, 1e-9, 1e-12):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    res, node = fn(floor, mode)
                n = len(np.asarray(res.sweep_values))
                v = float(np.asarray(res.v(node))[-1])
                base = n if base is None else base
                print('%-10s %-12s %13.0e %8d %14.6g %+8.1f%%'
                      % (label, mode, floor, n, v, 100.0 * (n - base) / base))
        print()
    print('Zero under `sigglobal` on both circuits, 8.5-9.2% under `pointlocal`.')
    print('The floor is only load-bearing in the mode whose reference can collapse.')


def gear2_accuracy():
    """WHAT 4i's EXTRA STEPS BOUGHT -- the measurement a FAILED gate clause needed.

    Gate 4i-4 declared that Gear2's rejection counts must not rise.  They rose, on 5
    of 9 configurations, because the old estimator did not merely over-report: it
    over-reported on shrink (83x at ratio 0.008) and UNDER-reported on growth (0.775
    at ratio 2), and a controller spends most of a run growing.  It rarely rejected
    because it was under-reporting, not because its steps were good.

    That claim is not self-evident, so this is the measurement that supports it:
    global error against the analytic solution, at matched `reltol`.

    Recorded outcome: `gear2-ywr` (the shipped default) 1.21-1.22x better accuracy
    for ~10% more steps; `gear2-classic` neutral; euler bit-identical.  The two Gear2
    rows now AGREE because 4i made `lte_formula` select nothing.
    """
    sys.path.insert(0, os.path.join(_REPO, 'pycircuit', 'circuit', 'tests'))
    import warnings
    import test_analysis_transient as T

    print('%-16s %8s %8s %9s %14s' % ('config', 'reltol', 'steps', 'rejects', 'error'))
    for name in ('gear2-classic', 'gear2-ywr', 'trap-classic', 'euler-classic'):
        for reltol in (1e-3, 1e-4, 1e-5, 1e-6):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                n, rej, err = T._stiff_run(name, reltol)
            print('%-16s %8.0e %8d %9d %14.6e' % (name, reltol, n, rej, err))
        print()
    print('Run this against the pre-4i source to see the other half: `gear2-ywr` read')
    print('2.4325e-03 / 5.4326e-04 / 1.1917e-04 at reltol 1e-4 / 1e-5 / 1e-6, in 63 /')
    print('102 / 185 steps.  More steps, and they buy accuracy -- but at MATCHED')
    print('accuracy the new estimator is marginally more expensive.  4i is `reltol`')
    print('meaning what it says, not a speedup, and the plan says so.')


def pi_gains():
    """WHAT FIXING THE PI GAINS BOUGHT END TO END (gate 4a-4).

    Gate 4a-4 compared `PIController` against the DEFAULT `IntegralController` and
    failed its rejections clause on a 0 -> 4 change in a 238-step run.  That clause
    was answering the wrong question: what judges 4a is PI before vs after its own
    fix, which is what this measures.

    Gustafsson's gains were used undivided, putting the closed loop outside the unit
    circle (radius 1.117 at p=2, 1.776 at p=3).  The growth clamp turned that into a
    permanent period-2 limit cycle rather than a blow-up, so it never looked like
    divergence -- but a controller alternating against its clamp overshoots the
    tolerance half the time, and each overshoot is a rejected step.

    Recorded outcome, before -> after: 34 -> 4 rejections (rc-vsin 1e-5), 13 -> 2 and
    29 -> 3 (stiff), 132 -> 61 (rc-pulse 1e-5), at essentially unchanged step counts.
    Up to 9.7x.
    """
    import warnings
    sys.path.insert(0, _HERE)
    import transient_stage4 as S
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import Gear2Integrator, ZERO_STABILITY_RATIO
    from pycircuit.circuit.stepcontroller import PIController

    class Spy(PIController):
        def __init__(self):
            super().__init__()
            self.h = []
            self.rej = 0
            self.forced = 0
            self._c = 0

        def evaluate_step(self, *a, **kw):
            acc, hn = super().evaluate_step(*a, **kw)
            forced = False
            if not acc:
                self.rej += 1
                self._c += 1
                if self._c > 3:            # transient.py's MAX_REJECT
                    forced, self._c = True, 0
                    self.forced += 1
            else:
                self._c = 0
            if acc or forced:
                self.h.append(kw['h_curr'])
            return acc, hn

    print('%-12s %8s %8s %8s %8s %10s %8s'
          % ('circuit', 'reltol', 'steps', 'reject', 'forced', 'worst', 'outside'))
    for cname, builder in S._fa_circuits():
        for reltol in (1e-3, 1e-4, 1e-5):
            cir, kw = builder()
            tran = Transient(cir, integrator=Gear2Integrator(), reltol=reltol)
            spy = Spy()
            tran.step_controller = spy
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                tran.solve(coupled_lte=False, **kw)
            r = [spy.h[i + 1] / spy.h[i] for i in range(len(spy.h) - 1)]
            print('%-12s %8.0e %8d %8d %8d %10.4f %8d'
                  % (cname, reltol, len(spy.h), spy.rej, spy.forced, max(r),
                     sum(1 for v in r if v > ZERO_STABILITY_RATIO)))
        print()
    print('Pre-4a the same runs gave 34 / 13 / 29 / 35 / 132 rejections where they now')
    print('give 4 / 2 / 3 / 24 / 61.  Step counts move by a few percent either way, so')
    print('the rejections were pure waste: the limit cycle overshooting its own clamp.')
    print()
    print('NOT a reason to make PI the default -- measured against IntegralController')
    print('it is comparable, at 0-13% MORE steps.  4a fixes an opt-in path.')


def chgtol_guard():
    """Decision 0.3d: the charge-domain guard, option (D), and why it is not built.

    Two independent refutations, both reproduced here:

      1. `Eg` is a CURRENT, so (D)'s `err_q = |Eg|/etol_q` divides amperes by
         coulombs.  Settled by scaling: `Eg/h^2` is constant, `Eg/h^3` is not.
      2. Repairing the units -- comparing `h*Eg`, a genuine charge, against
         `etol_q` -- makes the guard fire on most REJECTED steps, i.e. it
         disables step control.  That is gate 0.2b's JAX failure mode.
    """
    import warnings
    from pycircuit.circuit import numeric
    from pycircuit.circuit.integrator import (TrapezoidalIntegrator,
                                              Gear2Integrator)
    from pycircuit.circuit.stepcontroller import IntegralController
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.analysis import remove_row_col
    import transient_stage4 as S

    CHGTOL = 1e-14                      # SPICE's default

    print('0.3d -- the charge guard, option (D): REFUTED, resolved to (A)')
    print()
    print('(1) Eg is a current, not a charge.  A cubic charge history q = q3 t^3/6,')
    print('    uniform grid, halving h.  Whichever column is CONSTANT gives the units.')
    print()
    q3 = 1.0e6
    print('    %-22s %10s %14s %14s' % ('integrator', 'h', 'Eg/h^2', 'Eg/h^3'))
    for cls in (TrapezoidalIntegrator, Gear2Integrator):
        integ = cls()
        for h in (1e-6, 5e-7, 2.5e-7):
            ts = [0.0, -h, -2 * h, -3 * h]
            qs = [q3 * t ** 3 / 6.0 for t in ts]
            Eg, _ = integ.compute_lte(
                q_curr=np.array([qs[0]]),
                q_last=[np.array([qs[1]]), np.array([qs[2]]), np.array([qs[3]])],
                iq_last=[np.array([0.0]), np.array([0.0])],
                h_curr=h, h_last=h, h_last2=h,
                is_first_step=False, toolkit=numeric)
            e = abs(float(np.asarray(Eg)[0]))
            print('    %-22s %10.3g %14.6g %14.6g'
                  % (cls.__name__, h, e / h ** 2, e / h ** 3))
    print()
    print('    Eg/h^2 constant => Eg ~ h^2 q\'\'\' => C/s^3 * s^2 = C/s = AMPERES.')
    print('    So (D)\'s err_q = |Eg| / (TRTOL*(reltol*|q| + chgtol)) has units 1/s.')
    print()

    class _Probe(IntegralController):
        """Computes both criteria alongside the real one; changes no behaviour."""

        def __init__(self):
            super().__init__()
            self.rejects = self.fire_D = self.fire_repaired = 0

        def evaluate_step(self, *a, **kw):
            accept, h_next = super().evaluate_step(*a, **kw)
            if kw['no_history']:
                return accept, h_next
            tk, integ = kw['toolkit'], kw['active_integrator']
            Eg, _ = integ.compute_lte(
                q_curr=kw['q_curr'], h_curr=kw['h_curr'], q_last=kw['q_last_hist'],
                iq_last=kw['iq_last_hist'], h_last=kw['h_last'],
                is_first_step=False, toolkit=tk, h_last2=kw.get('h_last2'))
            Jr, Egr = remove_row_col((kw['J'], Eg), kw['irefnode'], tk)
            try:
                lte = np.abs(np.asarray(tk.linearsolver(Jr, Egr), dtype=float))
            except Exception:
                return accept, h_next
            Egv = np.abs(np.delete(np.asarray(Eg, dtype=float), kw['irefnode']))
            q = np.abs(np.asarray(kw['q_curr'], dtype=float))
            ab = np.delete(np.asarray(kw['abstol'], dtype=float), kw['irefnode'])
            x = np.abs(np.asarray(kw['x_curr'], dtype=float))
            etol_v = kw['TRTOL'] * (kw['reltol'] * np.delete(x, kw['irefnode']) + ab)
            etol_q = kw['TRTOL'] * (kw['reltol'] * q.max() + CHGTOL)
            err_v = float(np.max(lte / etol_v)) if etol_v.size else 0.0
            if err_v > 1.0:
                self.rejects += 1
                if float(np.max(Egv)) / etol_q <= 1.0:
                    self.fire_D += 1
                if float(kw['h_curr'] * np.max(Egv)) / etol_q <= 1.0:
                    self.fire_repaired += 1
            return accept, h_next

    print('(2) What each version of the guard would do, over real runs.')
    print()
    print('    %-12s %8s %9s %10s %20s'
          % ('circuit', 'reltol', 'rejects', '(D) fires', 'repaired fires'))
    for name, builder in S._fa_circuits():
        for reltol in (1e-4,):
            cir, kwargs = builder()
            tran = Transient(cir, toolkit=numeric, reltol=reltol)
            pr = _Probe()
            tran.step_controller = pr
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tran.solve(coupled_lte=False, **kwargs)
            frac = (100.0 * pr.fire_repaired / pr.rejects) if pr.rejects else 0.0
            print('    %-12s %8.0e %9d %10d %13d (%4.1f%%)'
                  % (name, reltol, pr.rejects, pr.fire_D, pr.fire_repaired, frac))
    print()
    print('    (D) as specified is inert by accident -- it fires only where BOTH')
    print('    tolerances sit on their absolute floors, so its threshold is the ratio')
    print('    chgtol/vabstol rather than anything about the circuit.  The repaired')
    print('    version overrides the controller on most rejections, which IS gate')
    print('    0.2b\'s "never rejects a step".  Neither is a guard.  Resolution: (A).')


def main():
    p = argparse.ArgumentParser()
    names = ('parity', 'relref', 'vabstol', 'matched_accuracy', 'lte_floor',
             'gear2_accuracy', 'pi_gains', 'chgtol_guard')
    for name in names + ('all',):
        p.add_argument('--' + name.replace('_', '-'), dest=name,
                       action='store_true')
    a = p.parse_args()
    chosen = [n for n in names if getattr(a, n)]
    if a.all or not chosen:
        chosen = list(names)
    for n in chosen:
        print('=' * 72)
        globals()[n]()
        print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
