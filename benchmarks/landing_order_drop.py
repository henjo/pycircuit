"""The order drop after a LANDED breakpoint: what it buys and what it costs,
in the stepping loop and in the shooting, against exact solutions.

THE MECHANISM.  `Transient._solve` re-arms `_is_first_step` after every step
that lands on a breakpoint (and after a force-accept), so a two-step method
(gear, trap) takes ONE backward-Euler step there.  The shooting did not
until 2026-09-28 -- it landed edges on its own frozen grid (`PSS.event_grid`)
and stepped on -- and `pss_transient_boundary.py` V1, on a smooth RC, found
that MORE accurate; this measured the other side, and the shooting now drops
too (`_InnerTransient.ORDER_DROP_AT_EDGES`, keyed to the edge's node).

'nodrop' (forward): the flag is cleared before each step unless the run is
genuinely opening (`_no_history`) -- which would also remove the drop after
a force-accept; the counts are printed (none on these fixtures).
'drop' (shooting): `ORDER_DROP_AT_EDGES`, the shooting's DEFAULT since
2026-09-28 (Andreas: "add the drop"); 'nodrop' switches it off.

MEASURED 2026-09-28 (predictions were written first; B's did not bind):

 A  smooth RC (tau = 0.3 T, tr = T/50), the node voltage, fixed grid: no drop
    is better, gear 3.0x and trap 86x at N = 400-800 (V1's numbers, forward).
    Adaptive: mixed -- trap 7.8x better without at reltol 1e-3, gear 1.6x
    WORSE; equal at 1e-5.
 B  a SINE drive: nothing to measure.  `Sin.next_event` returns only its
    start delay (stage 4g(a): a sine has no discontinuity), so no step ever
    lands and the two variants are bit-equal -- the comment in `_solve` that
    said a VSin fires every quarter period was stale.
 C  STIFF RC (tau = 1e-4 T << h), the resistor CURRENT k points after a
    corner (relative to its peak), N = 400:
      forward trap   drop 0.039 0.033 0.028 0.021   nodrop 0.85 0.73 0.62 0.45
      forward gear   drop 0.039 0.016 0.002 5e-5    nodrop 0.42 0.050 0.004 2e-5
      PSS trap       drop 0.039 ...                 nodrop 0.85 0.73 0.62 0.45
      PSS gear       drop 0.039 ...                 nodrop 0.42 0.050 0.004 2e-5
      PSS radau      0.060 0.0036 2e-4 8e-7 ; trbdf2 0.13 0.018 0.0025 5e-5
    Without the drop trap RINGS -- O(1), decaying by its stiff factor
    (1 - h/2tau)/(1 + h/2tau) = -0.85 per step here -- and gear is 10x off for
    one step (L-stable, it then damps).  The shooting, which never drops,
    shows exactly the forward no-drop numbers: its max |v| error is 22x
    (trap) and 11x (gear) what re-arming the drop gives.  The stage methods
    keep no history across the edge and need no drop.  Adaptive forward
    (reltol 1e-4): trap without the drop 0.70 at k = 1 against 0.18; gear
    about equal at k = 1 and BETTER without it after (0.013 against 0.086
    at k = 5; not explained).

So the drop is a TRADE: it costs a smooth state accuracy and it is what keeps
a stiff state (an algebraic current, a fast parasitic pole) from ringing.
Which way to go for the shooting's trap/gear is the owner's decision.
"""
import sys
import warnings

import numpy as np

from pycircuit.circuit import circuit
from pycircuit.circuit import transient as _T
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import C, R, VPulse
from pycircuit.circuit.integrator import Gear2Integrator, TrapezoidalIntegrator

circuit.default_toolkit = circuit.numeric
T = 1e-6
RR = 1e3
TD, PW = 0.0125 * T, 0.4 * T
MODE = ['drop']
_orig_step = _T.Transient.solve_timestep


def _step(self, *a, **k):
    if MODE[0] == 'nodrop' and not self._no_history:
        self._is_first_step = False
    return _orig_step(self, *a, **k)


_T.Transient.solve_timestep = _step


def pulsed(tau, tr):
    """`(build, exact, u_at, corners)` for the pulsed RC: `v' = (u - v)/tau`,
    `u` piecewise linear; `exact` from `v(0) = 0`."""
    def build():
        c = SubCircuit()
        c['vs'] = VPulse(1, gnd, v1=0.0, v2=1.0, td=TD, tr=tr, tf=tr, pw=PW,
                         per=T)
        c['R'] = R(1, 2, r=RR)
        c['C'] = C(2, gnd, c=tau / RR)
        return c

    def segments(tend):
        segs, k = [], 0
        while k * T < tend:
            o = k * T
            e = [o, o + TD, o + TD + tr, o + TD + tr + PW, o + TD + 2 * tr + PW,
                 o + T]
            segs += [(e[0], e[1], 0.0, 0.0), (e[1], e[2], 0.0, 1.0 / tr),
                     (e[2], e[3], 1.0, 0.0), (e[3], e[4], 1.0, -1.0 / tr),
                     (e[4], e[5], 0.0, 0.0)]
            k += 1
        return segs

    def prop(v0, t0, t, a, b):
        s = t - t0
        return a + b * s - b * tau + (v0 + b * tau - a) * np.exp(-s / tau)

    def exact(ts):
        ts = np.asarray(ts, float)
        out = np.empty_like(ts)
        v0, i = 0.0, 0
        for (t0, t1, a, b) in segments(ts[-1] + T):
            while i < len(ts) and ts[i] <= t1:
                out[i] = prop(v0, t0, ts[i], a, b)
                i += 1
            v0 = prop(v0, t0, t1, a, b)
            if i >= len(ts):
                break
        return out

    def u_at(ts):
        ts = np.asarray(ts, float)
        out = np.zeros_like(ts)
        for (t0, t1, a, b) in segments(ts[-1] + T):
            m = (ts >= t0) & (ts <= t1)
            out[m] = a + b * (ts[m] - t0)
        return out
    corners = np.array([TD, TD + tr, TD + tr + PW, TD + 2 * tr + PW])
    return build, exact, u_at, corners


def forward(fix, integ, mode, npts=None, reltol=None):
    build, exact, u_at, _c = fix
    MODE[0] = mode
    c = build()
    kw = {} if reltol is None else {'reltol': reltol}
    tr = _T.Transient(c, integrator=integ(), **kw)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(refnode=gnd, tend=2 * T, timestep=T / (npts or 100),
                       x0=np.zeros(c.n), fixed_timestep=npts is not None)
    ts = np.asarray(res.sweep_values, float)
    v = np.asarray(res.x, float)[[str(n) for n in c.nodes].index('2')]
    return ts, v, res.statistics


def after_corners(ts, err, corners, period, ks=(1, 2, 3, 5)):
    """max over the corners (in every period) of `err` k points after one"""
    cs = np.concatenate([corners + k * period
                         for k in range(int(np.ceil(ts[-1] / period)) + 1)])
    out = []
    for k in ks:
        m = 0.0
        for tc in cs:
            j = np.searchsorted(ts, tc * (1 + 1e-9), side='right') + k - 1
            if j < len(ts):
                m = max(m, err[j])
        out.append(m)
    return out


def current_err(fix, ts, v, shift=0.0):
    _b, exact, u_at, _c = fix
    ie = (u_at(ts) - exact(ts + shift)) / RR
    return np.abs((u_at(ts) - v) / RR - ie) / np.max(np.abs(ie))


def table_A():
    fix = pulsed(0.3 * T, 0.02 * T)
    print('A  smooth pulsed RC (tau = 0.3 T): max |v - exact|')
    print('   fixed grid   method   N    drop       nodrop     ratio')
    for name, integ in (('gear', Gear2Integrator), ('trap', TrapezoidalIntegrator)):
        for n in (100, 200, 400, 800):
            e = []
            for mode in ('drop', 'nodrop'):
                ts, v, _st = forward(fix, integ, mode, npts=n)
                e.append(np.max(np.abs(v - fix[1](ts))))
            print('                %-5s %4d  %.3e  %.3e  %6.2f'
                  % (name, n, e[0], e[1], e[0] / e[1]))
    print('   adaptive     method  reltol   drop (err, steps)  nodrop (err, steps)')
    for name, integ in (('gear', Gear2Integrator), ('trap', TrapezoidalIntegrator)):
        for rt in (1e-3, 1e-4, 1e-5):
            row = []
            for mode in ('drop', 'nodrop'):
                ts, v, st = forward(fix, integ, mode, reltol=rt)
                row += [np.max(np.abs(v - fix[1](ts))), st.accepted_steps,
                        st.force_accepts]
            print('                %-5s  %.0e   %.3e %5d      %.3e %5d   fa %d/%d'
                  % (name, rt, row[0], row[1], row[3], row[4], row[2], row[5]))


def table_C():
    fix = pulsed(1e-4 * T, 0.02 * T)
    print('C  STIFF pulsed RC (tau = 1e-4 T): the current, k = 1 2 3 5 points '
          'after a corner, relative')
    for name, integ in (('gear', Gear2Integrator), ('trap', TrapezoidalIntegrator)):
        for lab, kw in (('N=400', {'npts': 400}), ('N=800', {'npts': 800}),
                        ('rt=1e-4', {'reltol': 1e-4})):
            for mode in ('drop', 'nodrop'):
                ts, v, _st = forward(fix, integ, mode, **kw)
                p = after_corners(ts, current_err(fix, ts, v), fix[3], T)
                print('   forward %-5s %-8s %-6s ' % (name, lab, mode)
                      + ' '.join('%.2e' % x for x in p))
    from pycircuit.circuit.shooting import PSS
    MODE[0] = 'drop'                     # (the forward hook is inert here)
    for method in ('trap', 'gear', 'radau', 'trbdf2'):
        for N in (400, 800):
            for drop in ((False, True) if method in ('trap', 'gear') else (False,)):
                c = fix[0]()
                pss = PSS(c, method=method, reltol=1e-10)
                ## (the default since 2026-09-28; False is the old walk)
                pss.ORDER_DROP_AT_EDGES = bool(drop)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    pss.solve(period=T, timestep=T / N, break_events=True,
                              maxiterations=40)
                if not pss.converged:
                    print('   PSS %-6s N=%d not converged' % (method, N))
                    continue
                ts = np.asarray(pss.waveform[0], float).ravel()
                v = np.asarray(pss.waveform[1], float)[
                    [str(n) for n in c.nodes].index('2')]
                ## periodic exact: the transient one period on (settled)
                p = after_corners(ts, current_err(fix, ts, v, shift=T),
                                  fix[3], 2 * T)
                dv = np.max(np.abs(v - fix[1](ts + T)))
                print('   PSS     %-6s N=%-5d %-6s ' % (
                    method, N, 'drop' if drop else 'nodrop')
                    + ' '.join('%.2e' % x for x in p) + '   max|dv| %.2e' % dv)


if __name__ == '__main__':
    which = sys.argv[1:] or ['A', 'C']
    if 'A' in which:
        table_A()
    if 'C' in which:
        table_C()
