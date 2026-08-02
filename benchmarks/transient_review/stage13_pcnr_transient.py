"""STAGE 13 -- does PCNR's step reduction survive beyond one circuit?

On a half-wave rectifier, PCNR reached the same accuracy as classic limiting
using 4.1x to 6.8x FEWER time steps.  One circuit is not evidence about a method,
and gate 13-4's "off by default" was decided on DC measurements, where there is
no step controller and the effect could not appear at all.

This runs the diode-bearing stress circuits both ways and asks two questions
separately:

  1. Does the step reduction hold, at matched accuracy?  Accuracy is measured
     against a TIGHT run of the classic path, so a method cannot win by being
     fast and wrong.
  2. Is the mechanism the one claimed?  The hypothesis is that limiting leaves
     `G` linearised around `_vlim` while `i` is evaluated at the node voltage, so
     the Jacobian is not the derivative of the residual and the LTE estimate
     built from that solve is corrupted into forcing small steps.  If that is
     right, the ACCEPTED-STEP ERROR distribution should differ: limiting should
     sit well below its target, having been scared into small steps, where PCNR
     should sit on it.

Run:
    MPLBACKEND=Agg PYTHONPATH=. python benchmarks/transient_review/stage13_pcnr_transient.py
"""
import warnings

import numpy as np

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VSin, VPulse, Diode
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import stepcontroller as sc


def half_wave():
    c = SubCircuit()
    c['VS'] = VSin(1, gnd, va=10, freq=50)
    c['D1'] = Diode(1, 2)
    c['C1'] = C(2, gnd, c=100e-6)
    c['R1'] = R(2, gnd, r=1e3)
    return c, dict(tend=0.05, timestep=1e-4), (2, gnd)


def full_wave():
    c = SubCircuit()
    c['VS'] = VSin(1, 2, va=10, freq=50)
    c['R_src'] = R(2, gnd, r=1e-3)
    c['D1'] = Diode(1, 3); c['D2'] = Diode(gnd, 1)
    c['D3'] = Diode(2, 3); c['D4'] = Diode(gnd, 2)
    c['C1'] = C(3, gnd, c=10e-6)
    c['R1'] = R(3, gnd, r=100)
    return c, dict(tend=0.04, timestep=1e-4), (3, gnd)


def charge_pump():
    c = SubCircuit()
    c['VP'] = VPulse(1, gnd, v1=-5, v2=5, tr=1e-6, tf=1e-6, pw=5e-6, per=10e-6)
    c['C1'] = C(1, 2, c=1e-6)
    c['D1'] = Diode(gnd, 2); c['D2'] = Diode(2, 3)
    c['C2'] = C(3, gnd, c=1e-6)
    return c, dict(tend=50e-6, timestep=1e-7, uic=True), (3, gnd)


def clipper():
    c = SubCircuit()
    c['VS'] = VSin(1, gnd, va=20, freq=1e3)
    c['R'] = R(1, 2, r=100.0)
    c['D1'] = Diode(2, gnd); c['D2'] = Diode(gnd, 2)
    c['C'] = C(2, gnd, c=1e-9)
    return c, dict(tend=2e-3, timestep=1e-5), (2, gnd)


CASES = [('half wave', half_wave), ('full wave', full_wave),
         ('charge pump', charge_pump), ('clipper', clipper)]


def run(build, pcnr, reltol, collect_err=False):
    c, kw, probe = build()
    tran = Transient(c, toolkit=numeric, reltol=reltol, pcnr=pcnr)
    errs = []
    if collect_err:
        real = sc.IntegralController.evaluate_step

        def cap(self, *a, **k):
            accept, h = real(self, *a, **k)
            if accept and self.last_err is not None:
                errs.append(float(self.last_err))
            return accept, h
        sc.IntegralController.evaluate_step = cap
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(**kw)
    except Exception as e:
        return None, None, -1, errs, '%s' % type(e).__name__
    finally:
        if collect_err:
            sc.IntegralController.evaluate_step = real
    w = res.v(*probe)
    return (np.asarray(w.x, dtype=float).ravel(),
            np.asarray(w.y, dtype=float).ravel(),
            tran.statistics.accepted_steps, errs, None)


def main():
    print('QUESTION 1 -- does the step reduction hold, at matched accuracy?')
    print('(error against the classic path at reltol 1e-7)')
    print()
    print('%-13s %-9s %8s %13s %9s' %
          ('circuit', 'path', 'steps', 'max|err|', 'steps x'))
    for name, build in CASES:
        tr, vr, sr, _e, err = run(build, False, 1e-7)
        if err:
            print('%-13s reference failed: %s' % (name, err))
            continue
        base = None
        for label, use in (('limiting', False), ('pcnr', True)):
            t, v, st, _e, err = run(build, use, 1e-5)
            if err:
                print('%-13s %-9s %8s   %s' % (name, label, 'FAIL', err))
                continue
            e = float(np.max(np.abs(np.interp(tr, t, v) - vr)))
            if base is None:
                base = st
                print('%-13s %-9s %8d %13.4e %9s' % (name, label, st, e, '-'))
            else:
                print('%-13s %-9s %8d %13.4e %8.2fx' %
                      (name, label, st, e, base / max(st, 1)))
        print()

    print('QUESTION 2 -- is it the LTE estimate?')
    print('(accepted-step normalised error; the controller aims at 0.81)')
    print()
    print('%-13s %-9s %9s %9s %9s' % ('circuit', 'path', 'median', 'p10', 'p90'))
    for name, build in CASES:
        for label, use in (('limiting', False), ('pcnr', True)):
            _t, _v, _st, errs, err = run(build, use, 1e-5, collect_err=True)
            if err or not errs:
                print('%-13s %-9s   (no data)' % (name, label))
                continue
            e = np.array(errs)
            print('%-13s %-9s %9.4f %9.4f %9.4f'
                  % (name, label, np.median(e), np.percentile(e, 10),
                     np.percentile(e, 90)))
        print()


if __name__ == '__main__':
    main()
