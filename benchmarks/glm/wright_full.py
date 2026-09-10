"""The whole method from the free parameters: B~ by (3.7.8), A~ = B~^-1 J B~ by the transformed IRKS
condition (3.5.9), back transform B = Psi B~, A = A~ + lam I (W = I), U/V from the order conditions."""
import os, sys, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wright_ctor import construct, shifts, Cmat, Emat
np.set_printoptions(precision=12, suppress=True, linewidth=210)

def method(p, lam, c, beta, eps=0.0, Tmat=None):
    out = construct(p, lam, c, [0.0] * (p + 1), beta, eps, Tmat=Tmat)
    Btil, Psi = out['Btil'], out['Psi']
    J, K = shifts(p)
    Atil = np.linalg.solve(Btil, J @ Btil)          # (3.5.9): B~ A~ = J B~
    A = Atil + lam * np.eye(p + 1)
    B = Psi @ Btil
    C = Cmat(c, p); E = Emat(p)
    U = C - A @ C @ K
    V = E - B @ C @ K
    return dict(A=A, B=B, U=U, V=V, Atil=Atil, Btil=Btil, Psi=Psi,
                upper=float(np.max(np.abs(np.triu(Atil)))))

if __name__ == '__main__':
    p, lam = 2, 0.25
    m = method(p, lam, [0.25, 0.5, 1.0], [0.25, 0.0, 0.0])
    print('A =\n', m['A'])
    print("Wright's printed A = [[1/4,0,0],[1/6,1/4,0],[1/6,1/2,1/4]]")
    Ap = np.array([[0.25, 0, 0], [1/6, 0.25, 0], [1/6, 0.5, 0.25]])
    print('max|A - A_printed| = %.3e   strictly-upper part of A~ = %.3e' % (np.max(np.abs(m['A'] - Ap)), m['upper']))
    from glm_tools import order_by_polynomials, stability_report
    print('poly exactness:', order_by_polynomials(m['A'], m['U'], m['B'], m['V'], [0.25, 0.5, 1.0], p))
    print('stability:', stability_report(m['A'], m['U'], m['B'], m['V']))
