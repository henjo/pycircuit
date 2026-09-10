"""p = 3 and 4 by Wright's route: the only condition left on the free parameters (lambda in Table 3.3's band,
c distinct with c_s = 1, eps = 0, beta, T = permutation x unit lower triangular) is that
A~ = B~^-1 J B~ come out STRICTLY LOWER TRIANGULAR -- i.e. that A be diagonally implicit with diagonal lambda.
Least squares on that, per lambda and per permutation.  argv: p"""
import os, sys, itertools, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wright_full import method
from glm_tools import order_by_polynomials, stability_report
from scipy.optimize import least_squares
np.set_printoptions(precision=12, suppress=True, linewidth=220)
p = int(sys.argv[1]); r = p + 1
c = [(k + 1.0) / r for k in range(p)] + [1.0]
band = {2: (0.1804, 2.1856), 3: (0.2236, 0.5728), 4: (0.2480, 0.6760)}[p]
nlow = p * (p - 1) // 2
def Tof(x, perm):
    T = np.eye(p)
    k = 0
    for i in range(p):
        for j in range(i):
            T[i, j] = x[k]; k += 1
    P = np.zeros((p, p))
    for i, pi in enumerate(perm):
        P[i, pi] = 1.0
    TT = np.eye(r); TT[1:, 1:] = P @ T
    return TT
best = None
for perm in itertools.permutations(range(p)):
    for lam in np.linspace(band[0] + 0.01, band[1] - 0.01, 9):
        def resid(x):
            beta = list(x[:p]) + [0.0]
            try:
                m = method(p, float(lam), c, beta, 0.0, Tmat=Tof(x[p:], perm))
            except Exception:
                return np.ones(r * r) * 10.0
            return np.triu(m['Atil']).ravel()
        x0 = np.concatenate((np.full(p, 0.2), np.zeros(nlow)))
        try:
            sol = least_squares(resid, x0, max_nfev=3000, xtol=1e-15, ftol=1e-15)
        except Exception:
            continue
        cost = float(np.sum(sol.fun ** 2))
        if best is None or cost < best[0]:
            best = (cost, float(lam), perm, sol.x)
            if cost < 1e-24:
                break
    if best and best[0] < 1e-24:
        break
cost, lam, perm, x = best
print('p=%d  best |triu(A~)|^2 = %.3e  lambda %.4f  perm %s' % (p, cost, lam, perm))
m = method(p, lam, c, list(x[:p]) + [0.0], 0.0, Tmat=Tof(x[p:], perm))
print('A =\n', m['A']); print('c =', c); print('B =\n', m['B'])
print('poly:', order_by_polynomials(m['A'], m['U'], m['B'], m['V'], c, p))
print('stab:', stability_report(m['A'], m['U'], m['B'], m['V']))
