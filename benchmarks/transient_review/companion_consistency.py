"""Is `Diode.G` the exact derivative of `Diode.i`?

`doc/src/circuit/pcnr_limiting.rst` asserted for some time that it is NOT -- that
`G` linearises around the limited voltage while `i` evaluates at the node voltage,
"so the Jacobian and the residual are taken at different points". That claim is
the main argument for PCNR being a *correctness* fix rather than a robustness one,
so it is worth more than a reading of the source.

It is false. `Diode.i` ends with

    I_eff = I - g * (VD - v_nodes)

i.e. the companion form `i(x) = I(VD) + g(VD)(v(x) - VD)`, which is affine in `x`
with slope exactly `g(VD)` -- and `G` returns `g(VD)`.

This script re-measures that by central difference rather than asserting it, and
prints the true conductance at the node voltage alongside, to show how far the
companion can sit from the device while remaining a perfectly valid tangent.

Run:
    MPLBACKEND=Agg PYTHONPATH=. python benchmarks/transient_review/companion_consistency.py
"""
import numpy as np

from pycircuit.circuit import numeric
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import Diode

VLIM = 0.55
IS = 1e-14


def main():
    d = Diode(1, 2, IS=IS)
    d.toolkit = numeric
    d._vlim = VLIM

    VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
    print('Diode companion linearised at VD = %.3f V   (VT = %.6f V)' % (VLIM, VT))
    print()
    print('%8s %14s %14s %8s %16s' %
          ('node v', 'G[0,0]', 'd(i)/dx[0]', 'equal?', 'g at node v'))

    worst = 0.0
    for v in (0.50, VLIM, 0.70, -2.0):
        x = np.array([v, 0.0])
        G = np.asarray(d.G(x, defaultepar), dtype=float)

        ## Central difference, and the step matters: too small and the exponential's
        ## own round-off dominates, too large and the affine claim is not what is
        ## being tested.  The companion is exactly affine, so any h in a wide band
        ## gives the same answer -- which is itself the property under test.
        h = 1e-7
        up = np.asarray(d.i(x + np.array([h, 0.0]), defaultepar), dtype=float)
        dn = np.asarray(d.i(x - np.array([h, 0.0]), defaultepar), dtype=float)
        fd = (up - dn) / (2.0 * h)

        rel = abs(G[0, 0] - fd[0]) / max(abs(fd[0]), 1e-300)
        worst = max(worst, rel)
        g_true = IS * np.exp(v / VT) / VT
        print('%8.2f %14.6e %14.6e %8s %16.6e' %
              (v, G[0, 0], fd[0], 'yes' if rel < 1e-6 else 'NO', g_true))

    print()
    print('worst relative disagreement: %.3e' % worst)
    print()
    print('So the residual and the Jacobian ARE mutually consistent, and PCNR is')
    print('not a correctness fix here. What is questionable is the CHOICE of')
    print('linearisation point: VD comes from a limiter that depends on the')
    print('previous iterate and on the order devices are visited.')
    print()
    print('Note the last row. At -2 V the companion presents 6.77e-04 S where the')
    print('device has 9.49e-47 S. That is a valid tangent, not a bug -- but it is')
    print('why a STALE VD is silently a plausible matrix rather than an obviously')
    print('broken one, which is how the gate 13-6 defect survived three gates.')


if __name__ == '__main__':
    main()
