"""GEAR'S COVARIANCE FLOOR (plan item 3, 2026-09-25): gear's (BDF2) periodic
covariance on a SCALAR switched RC, exact 3x3 algebra, for candidate noise
injections -- the model that located the floor and decided that no fix is
buildable for the switches this simulator has.  See `doc/pss_log_260902.md`,
2026-09-25 (two entries).

    C x' = -G(t) x + w,   E[w w] = (CY/2) delta,   CY = 4 kT G(t)

so the exact variance is kT/C at every instant (fluctuation-dissipation).

Two switch shapes:
  abrupt  G jumps between two steps, the jump landed on a node (duty 0.5);
  smooth  the kT/C sampler's own `_SwitchHdl`: a tanh window of 50 mV on a
          1 V cosine clock.  This model reproduces the sampler's RECORDED
          numbers to 3-4 digits (rectangle hold 0.963 / 0.992 / 0.9987 /
          0.9998; stochastic-BDF2 hold 0.569 / 0.751 / 0.868 / 0.932) --
          it is the instrument, not an idealisation of it.
  (VSwitch is a C2 smoothstep whose window the event stage re-cuts into
  `event_window_steps` = 16 steps: smooth at the step level too.)

State z_n = (x_n, x_{n-1}, e_n), e_n the carried noise sample.  One BDF2 step,
constant h, (a0, a1, a2) = (3/2, -2, 1/2), Jf = a0 C/h + G_{n+1}.
Candidates:
  A    the rectangle (today's): one sample per step, variance CY/2h;
  E1   stochastic BDF2 (BDF2-Maruyama): (3/2) xi_{n+1} - (1/2) xi_n;
  E3   the exact solution's BDF2 residual, its noise part:
       Jf eta_{n+1} + (Jf Phi_{n+1,n} + a1 C/h) eta_n, eta the exact
       (Van Loan) increment of each step, Phi its propagator;
  E1R, E3R   E1 / E3 with the step after a switch taken as Euler in the
       step maps (a jump of more than e in G between steps);
  E3T  E3 with Euler wherever |d ln G| > THR between steps (the smooth
       transition's rule; `THR` below).
Run:  python benchmarks/pss_gear_covariance_floor.py [abrupt|smooth] [candidates...]
"""
import sys
import numpy as np

kT = 1.380649e-23 * 300.0
C = 100e-12
Gon, Goff = 1e-3, 1e-12
T = 1e-5
tau = C / Gon
a0, a1, a2 = 1.5, -2.0, 0.5
THR = 0.1


def G_at(t, shape):
    ph = (t / T) % 1.0
    if shape == 'abrupt':
        return Gon if ph < 0.5 else Goff
    v = np.cos(2.0 * np.pi * ph)
    return Goff + (Gon - Goff) * 0.5 * (1.0 + np.tanh(v / 0.05))


def step(cand, g, gp, h):
    """`(A, b, var)` of one step: `z_{n+1} = A z_n + b e_{n+1}`."""
    Phi = np.exp(-g * h / C)
    exact = cand.startswith('E3')
    var = (kT / C) * (1.0 - np.exp(-2.0 * g * h / C)) if exact else 2.0 * kT * g / h
    jump = abs(np.log(g / gp))
    euler = ((cand in ('E1R', 'E3R') and jump > 1.0)
             or (cand == 'E3T' and jump > THR))
    if euler:
        Jf = C / h + g
        if exact:
            A = [[(C / h) / Jf, 0.0, Phi - (C / h) / Jf], [1, 0, 0], [0, 0, 0]]
            b = [1.0, 0.0, 1.0]
        else:
            A = [[(C / h) / Jf, 0.0, 0.0], [1, 0, 0], [0, 0, 0]]
            b = [1.0 / Jf, 0.0, 1.0]
    else:
        Jf = a0 * C / h + g
        det = [-(a1 * C / h) / Jf, -(a2 * C / h) / Jf]
        if cand == 'A':
            A, b = [det + [0.0], [1, 0, 0], [0, 0, 0]], [1.0 / Jf, 0.0, 0.0]
        elif cand in ('E1', 'E1R'):
            A, b = [det + [-0.5 / Jf], [1, 0, 0], [0, 0, 0]], [1.5 / Jf, 0.0, 1.0]
        else:
            A = [det + [Phi + (a1 * C / h) / Jf], [1, 0, 0], [0, 0, 0]]
            b = [1.0, 0.0, 1.0]
    return np.array(A, dtype=float), np.array(b, dtype=float), var


def periodic_var(N, cand, shape):
    """The periodic variance of x at every node, in units of kT/C."""
    h = T / N
    G = np.array([G_at((j + 1) * h, shape) for j in range(N)])
    As, Qs = [], []
    for j in range(N):
        A, b, var = step(cand, G[j], G[j - 1], h)
        As.append(A)
        Qs.append(np.outer(b, b) * var)
    M, K1 = np.eye(3), np.zeros((3, 3))
    for A, Q in zip(As, Qs):
        M = A @ M
        K1 = A @ K1 @ A.T + Q
    K = np.linalg.solve(np.eye(9) - np.kron(M, M), K1.reshape(-1)).reshape(3, 3)
    xs = [K[0, 0]]
    for A, Q in zip(As, Qs):
        K = A @ K @ A.T + Q
        xs.append(K[0, 0])
    return np.array(xs) / (kT / C)


if __name__ == '__main__':
    args = sys.argv[1:]
    shapes = [a for a in args if a in ('abrupt', 'smooth')] or ['abrupt', 'smooth']
    cands = [a for a in args if a not in ('abrupt', 'smooth')] or \
        ['A', 'E1', 'E3', 'E1R', 'E3R', 'E3T']
    Ns = (100, 200, 400, 800, 1600)
    print('tau = %.1e s, T = %.1e s, h/tau = %s' % (
        tau, T, ' / '.join('%.3g' % (T / N / tau) for N in Ns)))
    for shape in shapes:
        ## mid-tracking and mid-hold phases of each shape
        tr_ph, ho_ph = (0.25, 0.75) if shape == 'abrupt' else (0.85, 0.35)
        print('%s switch' % shape.upper())
        for cand in cands:
            rows = [periodic_var(N, cand, shape) for N in Ns]
            print('  %-4s tracking %s | hold %s' % (
                cand,
                ' '.join('%.5f' % r[int(tr_ph * N)] for r, N in zip(rows, Ns)),
                ' '.join('%.5f' % r[int(ho_ph * N)] for r, N in zip(rows, Ns))))
