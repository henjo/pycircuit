"""GEAR'S COVARIANCE FLOOR (plan item 3, 2026-09-25): gear's (BDF2) periodic
covariance on a SCALAR switched RC, exact 3x3 algebra, for candidate noise
injections -- the model that located the floor and verified a fix.  The
rectangle (A) reproduces the circuit's recorded tracking floor (0.600 /
0.741 / 0.847 / 0.916 kT/C at 100 .. 800 points, measured on the sampler
0.598 / 0.740 / 0.847 / 0.916).  See `doc/pss_log_260902.md`, 2026-09-25.
  C x' = -G(t) x + w, E[w w] = (CY/2) delta,
CY = 4kT G(t) -> the exact variance is kT/C at every instant (FDT).

State z_n = (x_n, x_{n-1}, e_n), e_n the previous step's noise sample
(variance CY_n / 2h).  One BDF2 step, constant h:
  Jf x_{n+1} = -(a1 C/h) x_n - (a2 C/h) x_{n-1} + c1 e_{n+1} + c2 e_n
with (a0, a1, a2) = (3/2, -2, 1/2), Jf = a0 C/h + G_{n+1}.
Candidates:
  A   rectangle:              (c1, c2) = (1, 0)
  E1  stochastic BDF2:        (3/2, -1/2) everywhere
  E2  E1 except the step after a landed switch: (1, 0)  [history across a
      stiffness jump is not consistent]
  E1R E1, and the step after a landed switch taken as Euler (order
      dropped) in the step maps
  AR  the rectangle, with the same Euler step
The switch lands on a node (event breaking), duty 0.5.
Run:  python benchmarks/pss_gear_covariance_floor.py [A E1 E2 E1R AR]
"""
import sys
import numpy as np

kT = 1.380649e-23 * 300.0
C = 100e-12
Gon, Goff = 1e-3, 1e-12
T = 1e-5
tau = C / Gon


def periodic_cov(N, cand, duty=0.5):
    h = T / N
    a0, a1, a2 = 1.5, -2.0, 0.5
    non = int(round(duty * N))
    G = np.array([Gon if j < non else Goff for j in range(N)])   # G of step j (-> node j+1)
    As, Qs = [], []
    for j in range(N):
        g = G[j]
        Jf = a0 * C / h + g
        switched = (j == 0 and G[-1] != g) or (j > 0 and G[j - 1] != g)
        if cand == 'A':
            c1, c2 = 1.0, 0.0
        elif cand == 'E1':
            c1, c2 = 1.5, -0.5
        elif cand == 'E2':
            c1, c2 = (1.0, 0.0) if switched else (1.5, -0.5)
        elif cand in ('E1R', 'AR'):
            ## an order-dropped (Euler) step on the step after a landed
            ## switch: the deterministic history is not carried across it
            if switched:
                Jf = C / h + g
                c1, c2 = 1.0, 0.0
                A = np.array([[(C / h) / Jf, 0.0, 0.0],
                              [1.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])
                var_e = 4 * kT * g / (2 * h)
                b = np.array([c1 / Jf, 0.0, 1.0])
                As.append(A)
                Qs.append(np.outer(b, b) * var_e)
                continue
            c1, c2 = (1.5, -0.5) if cand == 'E1R' else (1.0, 0.0)
        A = np.array([[-(a1 * C / h) / Jf, -(a2 * C / h) / Jf, c2 / Jf],
                      [1.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0]])
        var_e = 4 * kT * g / (2 * h)
        b = np.array([c1 / Jf, 0.0, 1.0])
        As.append(A)
        Qs.append(np.outer(b, b) * var_e)
    M = np.eye(3)
    K1 = np.zeros((3, 3))
    for A, Q in zip(As, Qs):
        M = A @ M
        K1 = A @ K1 @ A.T + Q
    n = 3
    K0 = np.linalg.solve(np.eye(n * n) - np.kron(M, M), K1.reshape(-1)).reshape(n, n)
    seq = [K0]
    K = K0
    for A, Q in zip(As, Qs):
        K = A @ K @ A.T + Q
        seq.append(K)
    x = np.array([s[0, 0] for s in seq]) / (kT / C)
    return x, non


print('tau = %.1e, T = %.1e' % (tau, T))
for cand in sys.argv[1:] or ('A', 'E1', 'E2', 'E1R', 'AR'):
    rows = []
    for N in (100, 200, 400, 800, 1600):
        x, non = periodic_cov(N, cand)
        track = x[int(0.25 * N)]          # mid tracking
        hold = x[int(0.75 * N)]           # mid hold
        rows.append((N, track, hold, x[non], x[non + 1]))
    print(cand)
    for N, tr, ho, sw, sw1 in rows:
        print('   N %5d  h/tau %.3f  tracking %.6f  hold %.6f  (at switch %.6f, next %.6f)' % (N, T / N / tau, tr, ho, sw, sw1))
