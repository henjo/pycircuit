"""STAGE 12B -- Fang's coupled method against the standard adaptive controller.

`coupled_lte=True` now solves for the solution and the step size together
(DAC 2013 Figures 3 and 4): the LTE is eq (6)'s solution-space estimate, the
step correction is sec. 3.4's approximate Newton (eqs 17 and 18), and there is
no rejection path.

FOUR ERROR STATISTICS ARE REPORTED, because they answer different questions and
the usual single maximum answers the least interesting one.

  min     the quietest part of the run.  A controller that is working keeps this
          off the floor: an error orders below the rest means steps that were
          far smaller than they needed to be, which is the waste Fang's lower
          bound gamma_min exists to remove.
  mean    sensitive to the tail, so `mean` well above `median` says the error is
          concentrated in a few points rather than spread.
  median  what the error looks like across the run as a whole -- the quantity a
          step controller actually shapes.
  max     dominated by wherever the waveform is hardest, usually the opening
          step (accepted unevaluated, governed by no tolerance) or the fastest
          transient.

Sec. 4.1's claim is that the coupled method distributes the LTE more evenly
along the non-uniform grid.  That is a statement about the SHAPE of this
distribution, and a maximum on its own cannot show it in either direction.

Accuracy is measured against CLOSED FORMS, never a fine-mesh reference run --
comparing the integrator with itself cannot detect a step count bought by
integrating less accurately.

Run:
    MPLBACKEND=Agg PYTHONPATH=.:benchmarks python benchmarks/transient_review/stage12b_coupled.py
"""
import time
import warnings

import numpy as np

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, L, VSin
from pycircuit.circuit.transient import Transient

RELTOLS = (1e-4, 1e-5, 1e-6)


def rc_vsin():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c, dict(tend=2e-3, timestep=1e-5), 'b'


def rc_vsin_analytic(t, r=1e3, c=1e-7, va=1.0, freq=1e3):
    tau = r * c
    w = 2.0 * np.pi * freq
    A = va / np.sqrt(1.0 + (w * tau) ** 2)
    phi = np.arctan(w * tau)
    return A * (np.sin(w * t - phi) + np.sin(phi) * np.exp(-t / tau))


def stiff_rlc():
    c = SubCircuit()
    c['C1'] = C(1, gnd, c=1e-6)
    c['R1'] = R(1, 2, r=1.0)
    c['L1'] = L(2, gnd, L=1e-6)
    x0 = np.zeros(c.n)
    x0[c.get_node_index('1')] = 1.0
    x0[c.get_node_index('2')] = 1.0
    return c, dict(tend=5e-3, timestep=2e-4, x0=x0), '1'


def stiff_rlc_analytic(t, r=1.0, l=1e-6, c=1e-6):
    alpha = r / (2.0 * l)
    wd = np.sqrt(1.0 / (l * c) - alpha ** 2)
    return np.exp(-alpha * t) * (np.cos(wd * t) + (alpha / wd) * np.sin(wd * t))


## A circuit with REAL breakpoints, and the reason it is here: the coupled path
## had no breakpoint handling at all until gate 12-4, and the omission was
## invisible because neither circuit above has an edge.  An RC driven by a
## trapezoidal pulse has a closed form -- each segment is either an exponential
## relaxation toward a constant or toward a ramp -- so it can be measured the
## same way as the others rather than against a fine-mesh proxy.
PULSE = dict(td=1e-5, tr=1e-6, pw=2e-5, tf=1e-6, per=5e-5, v1=0.0, v2=1.0)
PULSE_TAU = 1e3 * 1e-9          # R*C below


def rc_pulse():
    from pycircuit.circuit.elements import VPulse
    c = SubCircuit()
    c['vs'] = VPulse('a', gnd, **PULSE)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    return c, dict(tend=1.2e-4, timestep=1e-6), 'b'


def rc_pulse_analytic(t):
    """Exact response of an RC low-pass to the trapezoidal drive.

    Integrated segment by segment.  On a segment where the source is the ramp
    ``u = a + b(t-t0)`` the solution of ``tau v' + v = u`` is

        v(t) = a + b(s - tau) + (v0 - a + b tau) exp(-s/tau),   s = t - t0

    which covers the flat segments too (b = 0).  Marching the segments rather
    than writing one formula keeps each piece checkable against that identity.
    """
    tau = PULSE_TAU
    td, tr, pw, tf, per = (PULSE['td'], PULSE['tr'], PULSE['pw'],
                           PULSE['tf'], PULSE['per'])
    v1, v2 = PULSE['v1'], PULSE['v2']

    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.empty_like(t)

    ## Segment table for one period: (duration, start value, slope).
    segs = [(td, v1, 0.0),
            (tr, v1, (v2 - v1) / tr),
            (pw, v2, 0.0),
            (tf, v2, (v1 - v2) / tf),
            (per - td - tr - pw - tf, v1, 0.0)]

    for i, ti in enumerate(t):
        v = v1                      # the source starts at v1 and so does the RC
        base = 0.0
        k = 0
        done = False
        while not done:
            for j, (dur, a, b) in enumerate(segs):
                ## EVERY period repeats the `td` lead-in.  `Pulse.f` folds time
                ## with `t % per` and then applies td/tr/pw/tf inside the folded
                ## period, so the shape is identical each time round -- skipping
                ## `td` after the first period put every later edge 10 us early
                ## and gave a "closed form" that disagreed with a reltol=1e-8 run
                ## by 1.0 V, i.e. full scale.
                if dur <= 0.0:
                    continue
                seg_end = base + dur
                s = min(ti, seg_end) - base
                if s > 0.0:
                    v = a + b * (s - tau) + (v - a + b * tau) * np.exp(-s / tau)
                if ti <= seg_end:
                    done = True
                    break
                base = seg_end
            k += 1
            if k > 100:
                done = True
        out[i] = v
    return out


CIRCUITS = [('rc-vsin', rc_vsin, rc_vsin_analytic),
            ('stiff-rlc', stiff_rlc, stiff_rlc_analytic),
            ('rc-pulse', rc_pulse, rc_pulse_analytic)]


def run(builder, analytic, reltol, coupled, repeats=3):
    ## GATE 12-3, the overhead claim.  Wall clock is the only thing that can
    ## answer it: sec. 3.4 says the approximate Newton "carries very little
    ## overhead", and a solve count cannot see the extra `J^-1 p` solve, the
    ## extrapolation, or the gradient assembly.
    ##
    ## Best-of-N, not mean: the minimum is the run least disturbed by whatever
    ## else the machine was doing, and scheduling noise only ever adds time.
    best = None
    for _ in range(repeats):
        cir, kw, node = builder()
        tran = Transient(cir, toolkit=numeric, reltol=reltol)
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(coupled_lte=coupled, **kw)
        wall = time.perf_counter() - t0
        if best is None or wall < best:
            best, best_res, best_tran = wall, res, tran
    res, tran, wall = best_res, best_tran, best

    w = res.v(node)
    t = np.asarray(w.x, dtype=float).ravel()
    v = np.asarray(w.y, dtype=float).ravel()
    d = np.abs(v - analytic(t))

    ## Past the opening step, which is accepted unevaluated -- no tolerance
    ## governs it, and on the stiff circuit it otherwise dominates the maximum
    ## at 1.359e-2 regardless of settings.
    d = d[2:] if len(d) > 2 else d

    st = tran.statistics
    return dict(wall=wall, steps=st.accepted_steps, rejected=st.rejected_steps,
                solves=st.accepted_steps + st.rejected_steps,
                min_err=float(np.min(d)), max_err=float(np.max(d)),
                med_err=float(np.median(d)), avg_err=float(np.mean(d)),
                t_end=float(t[-1]), tend=kw['tend'])


def main():
    print('Fang coupled method vs the standard adaptive controller')
    print('(errors against closed forms; opening step excluded)')
    print()
    for name, builder, analytic in CIRCUITS:
        print('--- %s' % name)
        print('%-8s %-9s %7s %7s %7s %9s %9s %11s %11s %11s'
              % ('reltol', 'path', 'steps', 'reject', 'solves',
                 'wall (s)', 'us/solve', 'mean|err|', 'median|err|', 'max|err|'))
        for reltol in RELTOLS:
            base = None
            for label, coupled in (('standard', False), ('coupled', True)):
                try:
                    r = run(builder, analytic, reltol, coupled)
                except Exception as e:
                    print('%-8.0e %-9s  raised %s: %s'
                          % (reltol, label, type(e).__name__, str(e)[:44]))
                    continue
                if base is None:
                    base = r
                done = '' if r['t_end'] >= r['tend'] * (1 - 1e-6) else '  RUN INCOMPLETE'
                print('%-8.0e %-9s %7d %7d %7d %9.3f %9.1f %11.3e %11.3e %11.3e%s'
                      % (reltol, label, r['steps'], r['rejected'], r['solves'],
                         r['wall'], 1e6 * r['wall'] / max(r['solves'], 1),
                         r['avg_err'], r['med_err'], r['max_err'], done))
        print()

    print('`us/solve` is what gate 12-3 turns on: the coupled path does more')
    print('work per Newton solve (an extra J^-1 p, the extrapolation, the')
    print('gradients), and the gate declared that per-step cost must stay within')
    print('15% of the standard path. Total wall clock also carries the step')
    print('count, so the two columns say different things.')


if __name__ == '__main__':
    main()
