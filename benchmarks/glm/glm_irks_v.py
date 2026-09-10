"""IRKS construction proper (Wright): impose BA = XB and BU = XV - VX with X doubly companion, so M(z) has the
single nonzero eigenvalue R(z) by theorem; plus V spectrum {1, 0..}, M_inf nilpotent (epsilon = 0 -> L-stable),
stiff accuracy (B row 1 = A row s), strict (B row 2 = e_s), c_s = 1, lambda pinned in Table 3.3's band.
argv: p lambda [seed]"""
import os, sys, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from glm_tools import uv_from, M, stability_report, order_by_polynomials
from scipy.optimize import least_squares
p = int(sys.argv[1]); lam = float(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
s = r = p + 1
nA = s * (s - 1) // 2; nc = s - 1; nB = (r - 2) * s
def unpack(x):
    A = lam * np.eye(s); k = 0
    for i in range(s):
        for j in range(i):
            A[i, j] = x[k]; k += 1
    c = np.concatenate((x[k:k + nc], [1.0])); k += nc
    B = np.zeros((r, s)); B[0] = A[-1]; B[1, s - 1] = 1.0
    B[2:] = x[k:k + nB].reshape(r - 2, s); k += nB
    alpha = x[k:k + r]; k += r
    beta = x[k:k + r]; k += r
    ## doubly companion in the layout MEASURED on Wright's p = 2 tableau (X = B A B^-1):
    ## first row alpha, ones on the SUBdiagonal, last column beta; all eigenvalues = lambda
    X = np.zeros((r, r))
    for i in range(1, r):
        X[i, i - 1] = 1.0
    X[0, :] += alpha; X[:, -1] += beta
    return A, c, B, X
def resid(x):
    A, c, B, X = unpack(x); U, V = uv_from(A, c, B, p)
    res = list((B @ A - X @ B)[1:].ravel())   # Wright: BA == XB is ALSO modulo the first row (book p.62)
    res += list((B @ U - X @ V + V @ X)[1:].ravel())          # holds modulo the FIRST row (measured)
    cp = np.poly(X); target = np.poly(lam * np.ones(r))
    res += list((cp - target) * 3.0)                            # sigma(X) = {lambda}: denominator (1 - lam z)^(p+1)
    ## peer (Wright p.83): V = [[1, v~^T], [0, Vdot]] with Vdot NILPOTENT is a STRUCTURE -- impose it as a
    ## polynomial condition (Vdot^p = 0, and V[1:, 0] = 0, V[0, 0] = 1) rather than an eigenvalue condition on a
    ## nearly defective block; same for M_inf (M_inf^r = 0)
    res += [V[0, 0] - 1.0] + list(V[1:, 0]) + list(np.linalg.matrix_power(V[1:, 1:], p).ravel() * 3.0)
    Mi = V - B @ np.linalg.solve(A, U)
    res += list(np.linalg.matrix_power(Mi, r).ravel() * 3.0)      # epsilon = 0: L-stable
    for y in (0.3, 1.0, 3.0, 10.0, 30.0, 100.0):
        rho = np.max(np.abs(np.linalg.eigvals(M(complex(0, y), A, U, B, V)))); res.append(max(rho - 1.0, 0.0) * 30.0)
    res += list(0.003 * (c[:-1] - np.linspace(1.0 / s, (s - 1.0) / s, s - 1)))
    ## non-confluence (Theorem 3.24): distinct abscissae -- penalise any pair closer than 0.05
    for i in range(s):
        for j in range(i):
            res.append(max(0.05 - abs(c[i] - c[j]), 0.0) * 30.0)
    return np.array(res)
rng = np.random.default_rng(seed); best = None
for t in range(80):
    x0 = np.concatenate((rng.normal(0, 0.3, nA), np.linspace(1.0 / s, (s - 1.0) / s, s - 1) + rng.normal(0, 0.05, s - 1), rng.normal(0, 0.8, nB), rng.normal(0, 0.5, r), rng.normal(0, 0.5, r)))
    try:
        sol = least_squares(resid, x0, max_nfev=8000, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    except Exception:
        continue
    cost = float(np.sum(sol.fun[:-(s - 1)] ** 2))
    if best is None or cost < best[0]:
        best = (cost, sol.x)
        if cost < 1e-20: break
A, c, B, X = unpack(best[1]); U, V = uv_from(A, c, B, p)
np.set_printoptions(precision=15, linewidth=240)
print('IRKS p=%d lambda=%.4f cost %.3e' % (p, lam, best[0])); print('A=', repr(A)); print('c=', repr(c)); print('B=', repr(B)); print('X=', repr(X))
print('stability', stability_report(A, U, B, V)); print('poly', order_by_polynomials(A, U, B, V, c, p))
