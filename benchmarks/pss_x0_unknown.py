"""`x_0` as the unknown instead of `x_in` -- the formulation change three
findings point at, prototyped and measured before anything is shipped.

WHY IT IS PROPOSED.  The plain path solves for `x_in`, the pre-image of a
manufactured opening step, while its Jacobian `I - dx_end/dx_0` is taken
with respect to `x_0`.  Three independent findings converge on that frame
error:

  1. The true `dF/dx_in` is SINGULAR on every circuit measured (rank 1/3,
     2/4, 1/3; sigma_min exactly 0), so no exact Newton exists for the
     present formulation and the path is a CONTRACTION -- its residual falls
     linearly at a constant ratio (~0.076 on the autonomous phase circuit)
     where a solved-history route is quadratic.
  2. Trapezoidal's period column is O(h) because `Pt` opens at ZERO although
     `x_0` is manufactured from a step of size `c_0 T` and does depend on T.
  3. Index-2 consistent initialisation at the period boundary, which the
     manufactured flat history destroys.

⚠ THE OBVIOUS VERSION IS THE FOURTH DESIGN TO DIE ON THE (-1)^n MODE, and
this time it was PREDICTED rather than discovered.  Dropping the
manufacturing step removes the L-STABLE opener, so trapezoidal's period map
becomes `A_trap^K` -- and the L-stability theorem says that is singular at
even K on every MNA circuit.  Measured on the model problem before any
shooting code was written, `sigma_min(I - A^K)`:

      circuit  variant       K=100(even)   K=101(odd)   K=200(even)
      Q=20     trap  bare     0.0000e+00   6.18e-03     0.0000e+00
      Q=20     euler bare     6.31e-03     6.32e-03     7.33e-03
      Q=20     trap  E-open   6.04e-03     6.10e-03     9.84e-03
      RC       trap  bare     0.0000e+00   9.99e-01     0.0000e+00
      RC       trap  E-open   9.99e-01     9.99e-01     9.99e-01

THE DESIGN THAT SURVIVES keeps ONE L-stable step and moves it INSIDE the
period: the unknown is `x_0`, there is no manufacturing step, and the first
step of the period is order-dropped to Euler -- which `_begin_period`
already arranges, so the period map is `A_trap^{K-1} A_euler` and is
non-singular at every K.

WHAT IT BUYS, MEASURED (this file):

  * QUADRATIC CONVERGENCE.  Q=20 resonator, 2 iterations:
    2.9e+00 -> 9.1e-08 -> 1.1e-14, at EVEN and odd point counts alike.
  * VAN DER POL, item 5's payoff case, on its 1105-step LTE grid:
    -47.3 ppm in 4 iterations (1.5e-02 -> 6.5e-04 -> 1.9e-06 -> 2.6e-10)
    against the shipped path's -73.8 ppm.  That -47.3 is exactly what the
    throwaway driver in `pss_lte_grid.py` reached, and for the same reason:
    its unknown was `x_0` and it manufactured nothing.
  * EULER IS UNCHANGED to the digit (8.83917 / 12.28401 at 100 / 200
    points, both formulations), which is the control -- euler's
    manufacturing step IS an Euler step, so the two maps coincide.

⚠ AND WHAT IT COSTS, WHICH IS WHY THIS IS NOT SHIPPED ON A HUNCH.  Moving
the Euler step inside the period degrades the ORBIT, not just the opening.
Q=20 resonator against the analytic 20 V peak:

      npts   shipped (x_in)   x_0 formulation
       100      20.01273         19.76939
       200      20.02208         19.96123

18x worse at 100 points, 1.8x at 200 -- the gap closes as the single
first-order step is diluted, and the ORDER is preserved (one Euler step
contributes O(h^2) globally, the same order as trapezoidal's own error), but
the constant is real.

SO THE TRADE-OFF IS GRID-DEPENDENT, and both directions are measured:
uniform grid on a benign circuit favours the shipped path (the Euler step
stays outside the period); non-uniform grid on a stiff one favours `x_0`
(no flat-history seam, and quadratic convergence is what decides whether a
relaxation oscillator converges at all).  That is a decision about defaults,
not a bug, which is why this driver exists and the change does not.
"""
import sys
import os
import warnings
from copy import copy

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pycircuit import circuit
from pycircuit.circuit import SubCircuit, gnd, R, C as Cap, L, VSin
from pycircuit.circuit.shooting import PSS


def q20():
    """Q=20 series resonator DRIVEN AT RESONANCE; analytic peak 20 V."""
    Lv, Cv, Q = 1e-3, 1e-9, 20.0
    f0 = 1.0 / (2 * np.pi * np.sqrt(Lv * Cv))
    c = SubCircuit()
    c.add_node('n1')
    c.add_node('n2')
    c['vs'] = VSin(gnd, 'n1', va=1.0, freq=f0)
    c['L'] = L('n1', 'n2', L=Lv)
    c['C'] = Cap('n2', gnd, c=Cv)
    c['R'] = R('n1', 'n2', r=Q * np.sqrt(Lv / Cv))
    return c, 1.0 / f0


def shoot_x0(cir_f, method, npts, per, maxit=25):
    """`x_0` is the unknown and the traversal opens AT it.

    `_begin_period` seeds the rings from `q(x_0)` and marks the next step
    `is_first_step`, so the first step inside the period is the Euler opener
    the theorem requires.  Finite-difference Jacobian on purpose: the point
    is the FORMULATION, and an analytic one would mix the two questions.
    """
    pss = PSS(cir_f()[0], method=method, reltol=1e-12)
    times, hs = pss._period_grid(per, npts, None)
    m = pss.cir.n - 1

    def phi(x0):
        pss._begin_period(np.asarray(x0, dtype=float))
        x = copy(np.asarray(x0, dtype=float))
        for j, t in enumerate(times[1:]):
            x = copy(pss.solve_timestep(x, t, hs[j]))
        return np.asarray(x, dtype=float)

    def F(z):
        return np.asarray(z, dtype=float) - phi(z)

    z = np.zeros(m)
    norms = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(maxit):
            f0 = F(z)
            norms.append(float(np.linalg.norm(f0)))
            J = np.zeros((m, m))
            for j in range(m):
                d = 1e-7 * max(abs(z[j]), 1.0)
                zp = z.copy()
                zp[j] += d
                J[:, j] = (F(zp) - f0) / d
            try:
                dz = np.linalg.solve(J, -f0)
            except np.linalg.LinAlgError:
                return norms, None
            z = z + dz
            if np.linalg.norm(dz) < 1e-13 * max(np.linalg.norm(z), 1.0):
                break
        x = copy(z)
        vals = [x]
        pss._begin_period(z)
        for j, t in enumerate(times[1:]):
            x = copy(pss.solve_timestep(x, t, hs[j]))
            vals.append(x)
    V = np.asarray(vals, dtype=float)
    return norms, 0.5 * (V[:, 1].max() - V[:, 1].min())


def main():
    circuit.default_toolkit = circuit.numeric
    _cir, per = q20()
    print('Q=20 resonator at resonance, analytic peak 20 V.\n')
    print('%-6s %-5s %-6s | %-10s | %s'
          % ('method', 'npts', 'parity', 'peak', 'residual sequence'))
    print('-' * 88)
    for method in ('trap', 'euler'):
        for npts in (100, 101, 200):
            norms, peak = shoot_x0(q20, method, npts, per)
            print('%-6s %-5d %-6s | %-10s | %s'
                  % (method, npts, 'even' if npts % 2 == 0 else 'odd',
                     ('%.5f' % peak) if peak is not None else 'SINGULAR',
                     '  '.join('%.1e' % v for v in norms[:6])))

    print('\nshipped plain path (x_in the unknown), same circuit:')
    for method in ('trap', 'euler'):
        for npts in (100, 200):
            pss = PSS(q20()[0], method=method, reltol=1e-12)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = pss.solve(period=per, timestep=per / npts,
                                maxiterations=40)
            v = np.asarray(res['tpss'].v('n2', gnd), dtype=float).ravel()
            print('   %-6s %-5d peak %.5f  converged=%s'
                  % (method, npts, 0.5 * (v.max() - v.min()), pss.converged))

    print('\nSee this module\'s docstring for the van der Pol figures and '
          'the trade-off; that case needs benchmarks/pss_lte_grid.py\'s grid.')


if __name__ == '__main__':
    main()
