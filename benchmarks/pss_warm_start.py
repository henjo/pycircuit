"""A4 PREMISE CHECK: is a warm start worth building, and when does it stop?

De Luca, Bolcato & Schilders (2019, TCAS-I) propose pre-integrating a
transient until the iterate is inside the shooting Newton's contraction
region, then handing off.  Their Algorithm 1 is our `_monodromy_matvec`, so
the machinery exists; what has to be decided before building is whether
pre-integration actually rescues the cases that fail, and how many periods
it needs.

⚠ THE VALUE IS CONFIRMED AND IT IS LARGE.  Seeded near the unstable DC point
-- the trivial-root basin, which is the case the literature describes --
van der Pol fails cold and converges after pre-integration:

      circuit                 cold solve        periods needed
      mu = 1    (strong)      LinAlgError             1
      mu = 0.05 (high-Q)      not converged          ~24

The ~24 is the `1/mu` amplitude-envelope time constant, so the count is a
property of how strongly the cycle attracts, not of the seed.  A large-
amplitude seed is far easier: from 4x and even 20x the orbit amplitude ONE
period suffices at mu = 1.

⚠ AND THE STOPPING CRITERION IS THE OPEN PART.  Carrying a probe vector
through the variational system each period -- `u^{k+1} = M_k u^k` -- and
watching its direction settle DOES NOT identify the handoff point.  Measured
on the high-Q case, drift against whether shooting then converges:

      k       carried probe      shoot from x_k
      0-20    1.4e-02 .. 1.8e-02   NOT converged
      24      5.7e-02              converged
      32      1.2e-01              converged
      59      2.1e-02              converged

The drift is SMALL while the solve still fails and LARGE once it succeeds --
it moves the wrong way.  The reason is diagnostic rather than numerical:
near the DC point the monodromy is nearly constant, so the probe settles
into ITS dominant eigenvector and reports "converged" while the state is
stuck at the trivial root.  A settled probe means "the Jacobian is not
changing", which is true both on the orbit and at the equilibrium.

⚠ THE OBVIOUS ALTERNATIVES HAVE THE SAME DEFECT, which is why this is not
worth guessing at.  "The state stopped changing period to period" and "the
shooting residual is small" are both ALSO true at the trivial root -- it is
a fixed point of the period map, so it satisfies every periodicity test.
Distinguishing it needs something amplitude-like, which is what
`DEGENERATE_PERIOD_FACTOR` does for the period.

⚠ SO THE BLOCKER IS THE PAPER, NOT THE CODE.  The criterion is stated there
as `||u^k - u~^k|| <= eps` and what `u~^k` is has been GUESSED at here as
the previous iterate.  That guess produces the behaviour above.  Get the
definition before building; the machinery is ready and the value is
established, so this is the only thing standing in the way.

Run: `python benchmarks/pss_warm_start.py`
"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pycircuit import circuit
from pycircuit.circuit import SubCircuit, gnd, C as Cap, L, BSource
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.shooting import PSS


def van_der_pol(mu):
    """Small `mu` is the HIGH-Q case: Q ~ 1/mu and the amplitude envelope
    settles over ~1/mu periods, which is the slow approach a warm start
    exists for."""
    c = SubCircuit()
    c.add_node('v')
    c['C'] = Cap('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    return c


def study(name, mu, period, seed_scale, nper, stride=1, npts=200):
    cir = van_der_pol(mu)
    iref = cir.get_node_index(gnd)
    m = cir.n - 1

    def reduce_(f):
        return np.concatenate((f[:iref], f[iref + 1:]))

    x0 = np.zeros(cir.n)
    x0[cir.get_node_index('v')] = 2.0 * seed_scale
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(cir, reltol=1e-9).solve(
            refnode=gnd, tend=nper * period, timestep=period / (2 * npts),
            x0=x0)
    t = np.asarray(res.sweep_values, dtype=float).ravel()
    X = np.asarray(res.x, dtype=float)
    states = [reduce_(X[:, min(int(np.searchsorted(t, k * period)),
                               X.shape[1] - 1)])
              for k in range(nper)]

    probe_pss = PSS(van_der_pol(mu), method='trap', reltol=1e-9)
    times, hs = probe_pss._period_grid(period, npts, None)
    u = np.ones(m) / np.sqrt(m)
    prev = None

    print('\n%s -- seed at %gx the orbit amplitude' % (name, seed_scale))
    print('%-5s %-13s %-14s %s'
          % ('k', '|x_k - x_end|', 'carried probe', 'shoot from x_k'))
    print('-' * 60)
    for k, xk in enumerate(states):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _a, _e, M, _c = probe_pss._traverse(xk, period, times, hs,
                                                want_dT=False)
        u_new = np.asarray(M, dtype=float) @ u
        nrm = np.linalg.norm(u_new)
        if nrm > 0:
            u_new = u_new / nrm
        drift = (np.linalg.norm(u_new - prev) if prev is not None
                 else float('nan'))
        prev, u = u_new, u_new
        if k % stride and k < len(states) - 1:
            continue
        pss = PSS(van_der_pol(mu), method='trap', reltol=1e-9)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pss.solve(period=period, timestep=period / npts, x0=xk,
                          maxiterations=40)
            out = 'converged' if pss.converged else 'NOT converged'
        except Exception as exc:                          # noqa: BLE001
            out = type(exc).__name__
        print('%-5d %-13.4e %-14s %s'
              % (k, np.linalg.norm(xk - states[-1]),
                 ('%.4e' % drift) if k else '-', out))


def main():
    circuit.default_toolkit = circuit.numeric
    study('mu=1 strongly attracting, seed near DC', 1.0, 6.663293, 0.02, 6)
    study('mu=0.05 high-Q, seed near DC', 0.05, 6.2832, 0.02, 60, stride=4)


if __name__ == '__main__':
    main()
