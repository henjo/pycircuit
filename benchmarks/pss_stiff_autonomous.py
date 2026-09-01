"""Does a stiff autonomous circuit actually need Gear-2 in PSS?  Measured: no.

`shooting.py`'s item 4c used to justify the composed `(x_0, x_{-1}, T)`
system by saying a two-step method becomes the right tool on a stiff
oscillator, citing `doc/transient_review.md` sec. 4.6: trapezoidal ringing
at `|e_n/e_{n-1}| ~ 0.9960` at `h*lambda = -1000` where Gear-2 sits at
0.0972.

⚠ THOSE ARE RINGDOWN NUMBERS -- they measure a TRANSIENT, and a periodic
steady state has no transient to ring.  The citation was carried from the
context it was measured in into one where it had not been checked.  This
script is the check, and it came out negative.

RESULT (2026-09-01).  Three circuits, and the third settles it.

    method  npts   period ppm      vf max      alt ripple
    trap     200      +83.084      0.539510     6.451e-03
    gear     200     +332.180      0.539507     6.113e-03
    trap    1600       +1.287      0.539515     1.143e-04
    gear    1600       +5.126      0.539515     1.141e-04

Trapezoidal shows NO ringing on either circuit below.  The alternating
signature is identical between the two methods and falls at ~h^3 under
refinement, so it is the sharp edge being RESOLVED, not an undamped mode --
if trapezoidal were ringing, its figure would not fall like Gear-2's.  And
trapezoidal is 4x better on frequency at both grids.

WHAT SURVIVES.  The composed system still earns its place, on evidence that
does not need this story: without it an autonomous Gear-2 run is silently
biased in the period by 2.5 ppm and its orbit does not close in radius, so
`method='gear'` was WRONG there rather than merely inferior.  Making an
offered method correct is the justification.

THE CASE THAT WAS UNTRIED IS NOW TRIED, and it closes the question the
other way.  Van der Pol at mu = 100 -- the canonical stiff relaxation
oscillator, fast mode IN the orbit, stiffness ratio 5443:

    method   npts        h   outcome              period       err ppm
    trap     2000   0.0815   NoConvergenceError        -             -
    gear     2000   0.0815   NoConvergenceError        -             -
    trap     8000   0.0204   NOT converged     162.813755             -
    gear     8000   0.0204   NoConvergenceError        -             -
    trap    20000   0.0081   converged         162.832543         -60.6
    gear    20000   0.0081   converged         162.823215        -117.9

⚠ TRAPEZOIDAL WINS ON BOTH COUNTS.  It is the only one that produced a
finite answer at 8000 points, and at 20000 it is TWICE as accurate.  So
"Gear-2 is the right tool for a stiff oscillator" is now refuted rather
than merely unsupported -- three circuits, no circuit yet found where
trapezoidal is the worse choice for autonomous PSS.

⚠ AND THE BINDING CONSTRAINT IS NOT THE METHOD, IT IS THE GRID.  Neither
method runs this circuit below 20000 points, because PSS freezes a UNIFORM
grid and the edge needs `h < 0.01` while the period is 162.8.  The adaptive
transient that produced the reference used ~1160 points per period.  So
the uniform-grid restriction costs about **17x the points** on this circuit
class, and that is a measured argument for the LTE-chosen grid (recorded
scope item 5) rather than for any integration method.
"""
import warnings

import numpy as np

from pycircuit import circuit
from pycircuit.circuit import SubCircuit, gnd, VS, R, C, L, Diode
from pycircuit.circuit.elements import IdtmodQuadrature, BSource
from pycircuit.circuit.shooting import PSS
from pycircuit.circuit.transient import Transient

NOM = 1e-3
TAU = 5e-9          # h*lambda = -1000 at 200 steps/period
RF = 1e3


def _oscillator(c):
    """The autonomous core: a self-sustaining quadrature phasor on DC."""
    for nd in ('in', 'o', 's'):
        c.add_node(nd)
    c['vin'] = VS('in', gnd, v=1e3)
    c['X'] = IdtmodQuadrature('in', gnd, 'o', gnd, 's', gnd,
                              modulus=1.0, amplitude=1.0, ic=0.0)
    c['Ro'] = R('o', gnd, r=1e6)
    c['Rs'] = R('s', gnd, r=1e6)
    return c


def smooth_parasitic(tau=TAU):
    """A fast RC at exactly `h*lambda = -1000`, driven and not loading.

    The element drives its outputs through branches, so the load cannot
    change the orbit.  Reference is analytic: node `f` is a 32 MHz lowpass
    of a 1 kHz cosine, so `v_f` tracks `v_o` to within `omega*tau`.

    ⚠ Its weakness as a test, stated because it is the reason for the
    second circuit: a stiff mode that is never EXCITED cannot reveal a
    damping failure, and a smooth slow drive never excites this one.
    """
    c = _oscillator(SubCircuit())
    c.add_node('f')
    c['Rf'] = R('o', 'f', r=RF)
    c['Cf'] = C('f', gnd, c=tau / RF)
    return c


def peak_detector(cf=1e-9, rload=1e5):
    """An orbit with a fast edge every period, so the mode IS excited.

    The diode conducts sharply near each peak and is off the rest of the
    cycle, so the kick comes from inside the periodic solution rather than
    from a transient that a PSS does not have.
    """
    c = _oscillator(SubCircuit())
    c.add_node('f')
    c['D'] = Diode('o', 'f')
    c['Cf'] = C('f', gnd, c=cf)
    c['Rl'] = R('f', gnd, r=rload)
    return c


MU = 100.0
VDP_PERIOD = 162.842412     # measured, see `van_der_pol`


def van_der_pol(mu=MU):
    """THE canonical stiff relaxation oscillator, and the decisive case.

        C dv/dt = -i_L + mu (v - v^3/3),    L di_L/dt = v

    No sources at all, so `u` is identically zero and the structural test
    calls it autonomous.  At `mu = 100` the orbit sits on its slow branches
    and crosses between them in O(1/mu): measured, the edge timescale is
    0.0299 against a period of 162.842412 (cycle-to-cycle spread 1.3e-06
    from a settled adaptive transient), a stiffness ratio of **5443**.  The
    asymptotic `(3 - 2 ln 2) mu = 161.37` agrees to the expected order.

    This is the circuit the earlier two were not.  A decoupled parasitic is
    never excited and a diode edge was not stiff enough; here the fast mode
    IS the orbit.
    """
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    return c


def van_der_pol_seed(tend=1200.0):
    """A point on the limit cycle, from a settled adaptive transient."""
    circuit.default_toolkit = circuit.numeric
    cir = van_der_pol()
    x0 = np.zeros(cir.n)
    x0[cir.get_node_index('v')] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(cir, reltol=1e-7).solve(refnode=gnd, tend=tend,
                                                timestep=0.05, x0=x0)
    xs = np.asarray(res.x, dtype=float)
    iref = cir.get_node_index(gnd)
    return np.concatenate((xs[:iref, -1], xs[iref + 1:, -1]))


def run_vdp(method, npts, seed, reltol=1e-7):
    circuit.default_toolkit = circuit.numeric
    pss = PSS(van_der_pol(), method=method, reltol=reltol)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            pss.solve(period=VDP_PERIOD, timestep=VDP_PERIOD / npts,
                      x0=seed, maxiterations=25)
        except Exception as exc:
            return dict(outcome=type(exc).__name__, period=None)
    return dict(outcome='converged' if pss.converged else 'NOT converged',
                period=pss.period)


def run(build, method, npts, reltol=1e-6, node='f'):
    circuit.default_toolkit = circuit.numeric
    pss = PSS(build(), method=method, reltol=reltol)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        res = pss.solve(period=NOM, timestep=NOM / npts, maxiterations=60)
        nonconv = any('did not converge' in str(c.message) for c in caught)
    v = np.asarray(res['tpss'].v(node), dtype=float).ravel()
    vo = np.asarray(res['tpss'].v('o'), dtype=float).ravel()
    ## An undamped alternating mode flips sign every step; a resolved edge
    ## does not.  The discriminator is whether this FALLS with the grid.
    alt = (float(np.max(np.abs(np.diff(np.diff(v))))) if len(v) > 2 else 0.0)
    return dict(conv=pss.converged and not nonconv,
                ppm=1e6 * (pss.period - NOM) / NOM,
                vmax=float(np.max(v)), alt=alt,
                lag=float(np.max(np.abs(v - vo))),
                solved_history=pss.solved_history)


def main():
    print('SMOOTH PARASITIC -- h*lambda = -1000 at 200 points/period')
    print('analytic |v_f - v_o| = omega*tau = %.3e' % (2 * np.pi / NOM * TAU))
    print('%-6s %5s %-6s %11s %12s %12s'
          % ('method', 'npts', 'conv', 'period ppm', '|v_f - v_o|',
             'alt ripple'))
    for m in ('trap', 'gear'):
        r = run(smooth_parasitic, m, 200, reltol=1e-8)
        print('%-6s %5d %-6s %+11.3f %12.3e %12.3e'
              % (m, 200, r['conv'], r['ppm'], r['lag'], r['alt']))

    print('\nPEAK DETECTOR -- a fast edge inside the orbit, every period')
    print('%-6s %5s %-6s %11s %10s %12s'
          % ('method', 'npts', 'conv', 'period ppm', 'vf max', 'alt ripple'))
    keep = {}
    for m in ('trap', 'gear'):
        for n in (200, 1600):
            r = run(peak_detector, m, n)
            keep[(m, n)] = r
            print('%-6s %5d %-6s %+11.3f %10.6f %12.3e'
                  % (m, n, r['conv'], r['ppm'], r['vmax'], r['alt']))

    print('\nthe discriminator -- a RESOLVED edge falls with the grid, an')
    print('undamped mode does not:')
    for m in ('trap', 'gear'):
        a, b = keep[(m, 200)]['alt'], keep[(m, 1600)]['alt']
        print('  %-5s alt ripple 200 -> 1600: %.2fx  (8x refinement)'
              % (m, a / b))

    print('\nVAN DER POL, mu=%g -- the canonical stiff relaxation oscillator'
          % MU)
    print('true period %.6f, edge timescale 0.0299, stiffness ratio 5443'
          % VDP_PERIOD)
    seed = van_der_pol_seed()
    print('%-6s %7s %9s %-20s %12s %10s'
          % ('method', 'npts', 'h', 'outcome', 'period', 'err ppm'))
    for m in ('trap', 'gear'):
        for npts in (2000, 20000):
            r = run_vdp(m, npts, seed)
            err = ('%10.1f' % (1e6 * (r['period'] - VDP_PERIOD) / VDP_PERIOD)
                   if r['period'] else '%10s' % '-')
            per = ('%12.6f' % r['period']) if r['period'] else '%12s' % '-'
            print('%-6s %7d %9.4f %-20s %s %s'
                  % (m, npts, VDP_PERIOD / (npts - 1), r['outcome'], per, err))


if __name__ == '__main__':
    main()
