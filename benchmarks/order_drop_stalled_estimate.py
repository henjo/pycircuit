"""The stalled-estimate order drop, measured on JAX -- and refused.

The CPU drops a Gear-2 step to Euler in two places, both aimed at a stalled
2nd-order LTE estimate (one built on a THIRD divided difference, which stops
improving across a discontinuity):

  * `Gear2Integrator.check_order_drop`'s `h_curr/h_last < 0.1` -- a PROXY,
    whose own comment says it stands for "we have rejected three times at this
    time point", something an Integrator object cannot see;
  * stage 4b's `MAX_REJECT = 3` cap, which force-accepts WITH an order drop
    rather than growing the step 10x.

Measured value on the CPU: the proxy fires 3-6 times a run on the stiff RLC at
reltol 1e-5, and removing it took `Gear2('ywr')`/`Gear2('classic')` from 0
force-accepts to 1 each.  A force-accept commits an unbounded truncation
error, so trading one for a controlled Euler step is the point.

BOTH WERE IMPLEMENTED ON THIS BACKEND (2026-09-01) AND BOTH WERE REVERTED.
This script reproduces the two order-rule variants against the shipped code by
wrapping `effective_first_order`; run it before re-proposing either.

  1. Ported LITERALLY, neither rule ever fires where it matters, because the
     backends give up in different places.  The CPU force-accepts after three
     rejections; this loop keeps shrinking to the `dt_min` floor (~130
     halvings at the default 1e-18) and force-accepts there.  At the floor
     `dt` equals the previous accepted step, so the ratio never trips, and the
     point has been rejected zero times, so the count never reaches three.
     Measured: 1900 force-accepts, unchanged by either rule; waveforms
     identical to ~1e-7.

  2. Ported FAITHFULLY -- sending a converged floor step whose LTE failed back
     for one retry at order 1, which is where this backend actually gives up
     -- it works as advertised and is still wrong.  It cut force-accepts hard
     (stiff RLC 1900 -> 516/688/864 at reltol 1e-6/-7/-8) and made the ANSWER
     WORSE every time, against an error-controlled reference run with the
     floor lowered until nothing was force-accepted:

         reltol   circuit      force-accepts        max|run - ref|
                               off -> ON            off        ON        ratio
         1e-6     stiff RLC    1900 -> 516          5.142e-01  6.949e-01  1.35x
         1e-6     pulsed RC      12 -> 12           7.156e-04  1.114e-03  1.56x
         1e-7     stiff RLC    1900 -> 688          5.196e-01  6.946e-01  1.34x
         1e-7     pulsed RC      55 -> 55           6.935e-04  4.505e-03  6.50x
         1e-8     stiff RLC    1900 -> 864          5.208e-01  6.946e-01  1.33x
         1e-8     pulsed RC     600 -> 480          7.035e-04  4.498e-03  6.39x

     (That variant needed a change to the ACCEPT path, so it cannot be
     reproduced by patching alone; the numbers are recorded here instead.)

THE FINDING, which is the transferable part: **fewer force-accepts is not a
better answer.**  Dropping to order 1 makes the ESTIMATE pass while making the
actual error larger -- the order-1 estimator is agreeing with a less accurate
method about a smaller claim.  The force-accept count is a label on a step,
not a measurement of its error, and optimising the label made the waveform
worse in 6 of 6 configurations.

Why the CPU benefits and this backend does not: its `MAX_REJECT` cap turns a
stalled estimate into a force-accept within three retries, and the order drop
prevents that.  Here the deep `dt_min` floor already prevents it -- there is
nothing left for the order drop to save, only accuracy for it to spend.
"""
import argparse, warnings

import numpy as np
import jax.numpy as jnp

import pycircuit.circuit.circuit as _cm
from pycircuit.circuit import gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, L, VPulse
from pycircuit.circuit.toolkit import jaxtoolkit
import pycircuit.circuit.jaxtransient as jt

MAX_POINT_REJECTS = 3
_BASE = jt.effective_first_order


def _with_rules(shrink, rejects):
    """`effective_first_order` plus the two CPU order-drop rules.

    `point_rejects` does not exist on the shipped state, so the count rule is
    approximated here by the run-global `n_rejected` being non-zero at this
    point -- enough to show the rule firing, which is all this needs to do.
    """
    def f(state):
        out = _BASE(state)
        if shrink:
            h_prev = state.h_history[0]
            out = jnp.logical_or(out, jnp.logical_and(
                h_prev > 0.0, state.dt < 0.1 * h_prev))
        if rejects:
            out = jnp.logical_or(out, state.n_rejected >= MAX_POINT_REJECTS)
        return out
    return f


def stiff_rlc():
    c = SubCircuit(toolkit=jaxtoolkit)
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-6, tr=1e-9, tf=1e-9,
                     pw=5e-6, per=1e-5, toolkit=jaxtoolkit)
    c['R'] = R('a', 'b', r=1.0, toolkit=jaxtoolkit)
    c['L'] = L('b', 'c', L=1e-6, toolkit=jaxtoolkit)
    c['C'] = C('c', gnd, c=1e-9, toolkit=jaxtoolkit)
    c['R2'] = R('c', gnd, r=1e4, toolkit=jaxtoolkit)
    return c


def pulsed_rc():
    c = SubCircuit(toolkit=jaxtoolkit)
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-5, tr=1e-7, tf=1e-7,
                     pw=2e-5, per=5e-5, toolkit=jaxtoolkit)
    c['R'] = R('a', 'b', r=1e3, toolkit=jaxtoolkit)
    c['C'] = C('b', gnd, c=1e-9, toolkit=jaxtoolkit)
    return c


CASES = [('stiff RLC', stiff_rlc, 2e-5, 1e-7, 'c'),
         ('pulsed RC', pulsed_rc, 1.2e-4, 1e-6, 'b')]
MODES = {'neither': (False, False), 'shrink': (True, False),
         'rejects': (False, True), 'both': (True, True)}


def run(build, tend, ts, node, mode, minstep, reltol):
    jt.effective_first_order = _with_rules(*MODES[mode])
    saved = _cm.default_toolkit
    _cm.default_toolkit = jaxtoolkit
    try:
        tran = jt.JAXTransient(build(), reltol=reltol, minstep=minstep)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=tend, timestep=ts, uic=True)
        return (np.asarray(res.sweep_values, float).reshape(-1),
                np.asarray(res.v(node), float).reshape(-1), tran.statistics)
    finally:
        jt.effective_first_order = _BASE
        _cm.default_toolkit = saved


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reltol', type=float, default=1e-7)
    ap.add_argument('--minstep', type=float, default=1e-8,
                    help='raise it to make the dt floor reachable')
    a = ap.parse_args()
    print('reltol %g, minstep %g\n' % (a.reltol, a.minstep))
    print('%-10s %-8s %8s %8s %13s %12s' %
          ('circuit', 'mode', 'accept', 'reject', 'force_accept', 'vs ref'))
    for name, build, tend, ts, node in CASES:
        tr, vr, sr = run(build, tend, ts, node, 'neither', 1e-18, a.reltol)
        for mode in MODES:
            t, v, st = run(build, tend, ts, node, mode, a.minstep, a.reltol)
            dev = float(np.max(np.abs(v - np.interp(t, tr, vr))))
            print('%-10s %-8s %8d %8d %13d %12.3e' %
                  (name, mode, st.accepted_steps, st.rejected_steps,
                   st.force_accepts, dev))
        print('   (reference: %d points, %d force-accepts)\n'
              % (len(tr), sr.force_accepts))


if __name__ == '__main__':
    main()
