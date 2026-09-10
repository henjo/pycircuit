"""Wright's route WITH the stiffly accurate sub-class (3.9.2): beta_p = 0 (so strict stiff accuracy falls out of
the IRKS condition), B row 1 = A row s, B row 2 = e_s, on top of A~ strictly lower triangular.  Free: lambda in
Table 3.3's band, c distinct with c_s = 1, beta_1..beta_{p-1}, T = perm x unit lower triangular.  argv: p"""
import os, sys, itertools, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wright_full import method
from glm_tools import order_by_polynomials, stability_report
from scipy.optimize import least_squares
np.set_printoptions(precision=15, linewidth=230)
p = int(sys.argv[1]); r = p + 1
band = {2: (0.1804, 2.1856), 3: (0.2236, 0.5728), 4: (0.2480, 0.6760)}[p]
nlow = p * (p - 1) // 2
def Tof(x, perm):
    T = np.eye(p); k = 0
    for i in range(p):
        for j in range(i):
            T[i, j] = x[k]; k += 1
    P = np.zeros((p, p))
    for i, pi in enumerate(perm):
        P[i, pi] = 1.0
    TT = np.eye(r); TT[1:, 1:] = P @ T
    return TT
def unpack(x, lam, perm):
    beta = list(x[:p - 1]) + [0.0, 0.0]          # beta_p = 0 (Wright 3.9.2), beta_{p+1} irrelevant
    c = list(x[p - 1:2 * p - 1]) + [1.0]
    return method(p, lam, c, beta, 0.0, Tmat=Tof(x[2 * p - 1:], perm)), c
best = None
rng = np.random.default_rng(0)
for perm in itertools.permutations(range(p)):
    for lam in np.linspace(band[0] + 0.01, band[1] - 0.01, 7):
        for trial in range(3):
            def resid(x):
                try:
                    m, c = unpack(x, float(lam), perm)
                except Exception:
                    return np.ones(r * r + 2 * r + r * (r - 1) // 2) * 1e3
                res = list(np.triu(m['Atil']).ravel())
                res += list(m['B'][0] - m['A'][-1])
                res += list(m['B'][1] - np.eye(r)[r - 1])
                cc = np.array(c)
                res += [max(0.02 - abs(cc[i] - cc[j]), 0.0) * 10 for i in range(r) for j in range(i)]
                return np.array(res)
            x0 = np.concatenate((np.full(p - 1, 0.2) + rng.normal(0, 0.1, p - 1),
                                 np.array([(k + 1.0) / r for k in range(p)]) + rng.normal(0, 0.03, p),
                                 rng.normal(0, 0.2, nlow)))
            try:
                sol = least_squares(resid, x0, max_nfev=4000, xtol=1e-15, ftol=1e-15)
            except Exception:
                continue
            cost = float(np.sum(sol.fun ** 2))
            if best is None or cost < best[0]:
                best = (cost, float(lam), perm, sol.x)
                if cost < 1e-22:
                    break
        if best and best[0] < 1e-22: break
    if best and best[0] < 1e-22: break
cost, lam, perm, x = best
try:
    m, c = unpack(x, lam, perm)
except Exception as e:
    print('p=%d: the best point (cost %.3e, lambda %.4f, perm %s) does not reconstruct: %s'
          % (p, cost, lam, perm, type(e).__name__)); raise SystemExit(1)
print('p=%d cost %.3e  lambda %.6f  perm %s' % (p, cost, lam, perm))
print('A =', repr(m['A'])); print('c =', repr(np.array(c))); print('B =', repr(m['B']))
print('stiff acc: |B0-As| %.2e  |B1-e| %.2e   triu(A~) %.2e' % (np.max(np.abs(m['B'][0]-m['A'][-1])), np.max(np.abs(m['B'][1]-np.eye(r)[r-1])), np.max(np.abs(np.triu(m['Atil'])))))
print('poly:', order_by_polynomials(m['A'], m['U'], m['B'], m['V'], c, p))
print('stab:', stability_report(m['A'], m['U'], m['B'], m['V']))
