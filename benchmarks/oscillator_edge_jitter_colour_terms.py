"""Step 3 of the edge-jitter plan: do the DEFERRED coloured terms matter?

`PAC.oscillator_edge_jitter` with a coloured source takes its phase as the
colour fold's instant-specific increment and its transverse part as `2 A_col`
at every k -- deferring the transverse part's own correlation across periods
and its cross term with the coloured phase.  `pair` measures the omission:
the element (a Lorentzian `IS(noiseTau=...)`) against its exact realisation
(white noise through an RC, read by the exact white law, which contains every
term), `k_cycle^2` element / realisation - 1 at k = 1..8, on van der Pol
(Q = 8, radau) with the source on the tank or behind a slow RC node (the
coloured current into an RC node of time constant `slow` feeding the tank
through a unit transconductance: a transverse state with memory across
periods; the orbit is unchanged).

MEASURED 2026-09-29 (400 points; worst over k, always at k = 1):
    on the tank, tau = 0.3 / 1 / 3 / 10 T:        1.3e-3 / 1.3e-3 / 1.9e-3 / 4.0e-3
    tau = 0.3 T behind a slow node, 1 T / 3 T:    1.9e-1 / 5.9e-2
    tau = 0.01 T (nearly white) behind 1 T:       8.3e-2
  unchanged at 800 points and at 40 points per decade.  The gap is the slow
  node's memory in the transverse part -- the deferred terms -- whatever the
  source's spectrum.

`bank` measures what state augmentation would cost a 1/f source: the ripple
of a Lorentzian bank (the midpoint rule on 1/f = (2/pi) int (1/fc) /
(1 + (f/fc)^2) dln fc) -- set by where the bank stops beyond the band, not by
its density: 3-6 % with one decade beyond each edge, 0.3-0.8 % with two.

Usage:  python oscillator_edge_jitter_colour_terms.py pair TAU_T [SLOW_T] [NPTS]
        python oscillator_edge_jitter_colour_terms.py bank
"""
import sys
import time
import warnings

import numpy as np

from pycircuit.circuit import IS, PSS, BSource, C, L, R, SubCircuit, circuit, gnd
from pycircuit.circuit.shooting import PAC


def vdp(kind, tau_T, slow_T=None, npts=400):
    """`(cir, pss, reduced index of the tank, its rising mid-level crossing)`:
    the element ('coloured') or its realisation ('filtered')."""
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    Rf, Cf, g, Pw = 1.0, tau_T * T, 1e-2, 1e-4
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: mu * (u - u ** 3 / 3.0))
    tgt = 'v'
    if slow_T is not None:
        c.add_node('s')
        c['rs'] = R('s', gnd, r=1.0, noisy=False)
        c['cs'] = C('s', gnd, c=slow_T * T)
        c['gs'] = BSource('s', gnd, gnd, 'v', i_func=lambda u: 1.0 * u)
        tgt = 's'
    if kind == 'coloured':
        c['n'] = IS(tgt, gnd, i=0.0, noisePSD=g * g * Pw * Rf * Rf, noiseTau=Rf * Cf)
    else:
        c.add_node('f')
        c['nw'] = IS('f', gnd, i=0.0, noisePSD=Pw)
        c['rf'] = R('f', gnd, r=Rf, noisy=False)
        c['cf'] = C('f', gnd, c=Cf)
        c['gm'] = BSource('f', gnd, gnd, tgt, i_func=lambda u, _g=g: _g * u)
    pss = PSS(c, method='radau', reltol=1e-12)
    full = c.get_node_index('v')
    red = full if full < c.get_node_index(gnd) else full - 1
    x0 = np.zeros(c.n - 1)
    x0[red] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=300)
    assert pss.converged
    Xw = np.asarray(pss.waveform[1], float)
    grid = np.asarray(pss.factored_period().times, float)[:Xw.shape[1]]
    v = Xw[full]
    mid = 0.5 * (v.max() + v.min())
    j = next(k for k in range(2, len(v) - 2) if (v[k - 1] - mid) < 0 <= (v[k] - mid))
    tc = grid[j - 1] + (mid - v[j - 1]) / (v[j] - v[j - 1]) * (grid[j] - grid[j - 1])
    return c, pss, red, float(tc)


def pair(tau_T, slow_T=None, npts=400, ppd=20):
    t0 = time.time()
    cw, pw, rw, tw = vdp('filtered', tau_T, slow_T, npts)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        w = PAC(cw, toolkit=circuit.numeric).oscillator_edge_jitter(pw, rw, tw)
    t1 = time.time()
    ce, pe, re_, te = vdp('coloured', tau_T, slow_T, npts)
    f0 = 1.0 / float(pe.period)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        e = PAC(ce, toolkit=circuit.numeric).oscillator_edge_jitter(
            pe, re_, te, colour_fmin=1e-6 * f0, points_per_decade=ppd)
    t2 = time.time()
    r = e['k_cycle'] ** 2 / w['k_cycle'] ** 2 - 1.0
    print(f'tau {tau_T} T, slow {slow_T} T, {npts} points: element / realisation - 1 '
          f'at k = 1..8: {" ".join(f"{x:+.1e}" for x in r)}; transverse share '
          f'{2 * e["A"] / e["k_cycle"][0] ** 2:.1e}; {t1 - t0:.0f} s + {t2 - t1:.0f} s')
    return r


def bank():
    for decades in (4, 10):
        f = np.logspace(0, decades, 4000)
        for ext in (1, 2):
            for ppd in (1, 2, 3):
                npole = round((decades + 2 * ext) * ppd) + 1
                x = np.linspace(-ext * np.log(10), (decades + ext) * np.log(10), npole)
                fc = np.exp(x)
                w = (2.0 / np.pi) * (x[1] - x[0]) / fc
                S = (w[None, :] / (1.0 + (f[:, None] / fc[None, :]) ** 2)).sum(axis=1)
                print(f'{decades:2d} decades, {ext} beyond each edge, {ppd} per decade '
                      f'({npole} states): ripple {np.max(np.abs(S * f - 1.0)):.1e}')


if __name__ == '__main__':
    if sys.argv[1] == 'bank':
        bank()
    else:
        slow = float(sys.argv[3]) if len(sys.argv) > 3 and float(sys.argv[3]) > 0 else None
        pair(float(sys.argv[2]), slow, int(sys.argv[4]) if len(sys.argv) > 4 else 400)
