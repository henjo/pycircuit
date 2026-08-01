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


CIRCUITS = [('rc-vsin', rc_vsin, rc_vsin_analytic),
            ('stiff-rlc', stiff_rlc, stiff_rlc_analytic)]


def run(builder, analytic, reltol, coupled):
    cir, kw, node = builder()
    tran = Transient(cir, toolkit=numeric, reltol=reltol)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(coupled_lte=coupled, **kw)

    w = res.v(node)
    t = np.asarray(w.x, dtype=float).ravel()
    v = np.asarray(w.y, dtype=float).ravel()
    d = np.abs(v - analytic(t))

    ## Past the opening step, which is accepted unevaluated -- no tolerance
    ## governs it, and on the stiff circuit it otherwise dominates the maximum
    ## at 1.359e-2 regardless of settings.
    d = d[2:] if len(d) > 2 else d

    st = tran.statistics
    return dict(steps=st.accepted_steps, rejected=st.rejected_steps,
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
        print('%-8s %-9s %7s %7s %7s %11s %11s %11s %11s %7s'
              % ('reltol', 'path', 'steps', 'reject', 'solves',
                 'min|err|', 'mean|err|', 'median|err|', 'max|err|',
                 'med/max'))
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
                print('%-8.0e %-9s %7d %7d %7d %11.3e %11.3e %11.3e %11.3e %7.3f%s'
                      % (reltol, label, r['steps'], r['rejected'], r['solves'],
                         r['min_err'], r['avg_err'], r['med_err'], r['max_err'],
                         r['med_err'] / r['max_err'], done))
        print()

    print('`med/max` is the shape of the error distribution, not its size.')
    print('Sec. 4.1 claims the coupled method spreads the LTE more evenly along')
    print('the non-uniform grid; a HIGHER ratio is a flatter distribution, and a')
    print('maximum on its own cannot show that in either direction.')


if __name__ == '__main__':
    main()
