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

RESULT (2026-09-01):

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

WHAT WOULD REOPEN IT: a circuit where trapezoidal actually fails inside a
PSS run.  Two attempts here did not find one, and two attempts are not
proof of absence -- a decoupled parasitic is never excited, and the diode
below may not be stiff enough in conduction.  The case still untried is a
genuinely stiff RELAXATION oscillator, where the fast mode lives in the
orbit itself rather than on a side node.
"""
import warnings

import numpy as np

from pycircuit import circuit
from pycircuit.circuit import SubCircuit, gnd, VS, R, C, Diode
from pycircuit.circuit.elements import IdtmodQuadrature
from pycircuit.circuit.shooting import PSS

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


if __name__ == '__main__':
    main()
