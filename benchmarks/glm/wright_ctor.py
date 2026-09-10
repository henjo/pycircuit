"""Wright's approach one (Thm 3.24, eq 3.7.8) implemented from the docs session's transcription, validated in
the prescribed order: (1) delta, (2) N_k from JF-FJ against the closed form, (3) B~ against Wright's printed B."""
import numpy as np
from math import factorial, comb
np.set_printoptions(precision=12, suppress=True, linewidth=200)

def shifts(p):
    r = p + 1
    J = np.zeros((r, r)); K = np.zeros((r, r))
    for i in range(r - 1):
        J[i + 1, i] = 1.0      # ones on the SUBdiagonal
        K[i, i + 1] = 1.0      # ones on the SUPERdiagonal
    return J, K

def Emat(p):
    r = p + 1
    E = np.zeros((r, r))
    for i in range(r):
        for j in range(i, r):
            E[i, j] = 1.0 / factorial(j - i)
    return E

def Cmat(c, p):
    c = np.asarray(c, float)
    return np.array([[ck ** k / factorial(k) for k in range(p + 1)] for ck in c])

def Fmat(lam, p):
    """F = exp( K (I + lam K)^-1 )  -- exp OF THE WHOLE PRODUCT (the trap)."""
    _J, K = shifts(p)
    r = p + 1
    Mprod = K @ np.linalg.inv(np.eye(r) + lam * K)
    F = np.zeros((r, r))
    term = np.eye(r)
    for k in range(r + 1):
        F = F + term / factorial(k)
        term = term @ Mprod
    return F

def N_closed(n, lam):
    return sum(comb(n - 1, i) * (-lam) ** i / factorial(n - i) for i in range(n))

def M_closed(n, lam):
    return sum(comb(n, i) * (-lam) ** i / factorial(n - i) for i in range(n + 1))

def polyK(coeffs, K):
    """1 + c_1 K + ... + c_{p+1} K^{p+1}"""
    r = K.shape[0]
    out = np.eye(r); term = np.eye(r)
    for c in coeffs:
        term = term @ K
        out = out + c * term
    return out

def LUsplit(M):
    """M = L(M) U(M): unit lower triangular L, upper triangular U (Doolittle, no pivoting)."""
    n = M.shape[0]
    L = np.eye(n); U = np.array(M, dtype=float)
    for k in range(n):
        if abs(U[k, k]) < 1e-300:
            raise np.linalg.LinAlgError('LU without pivoting hit a zero pivot at %d' % k)
        for i in range(k + 1, n):
            f = U[i, k] / U[k, k]
            L[i, k] = f
            U[i, k:] = U[i, k:] - f * U[k, k:]
        U[k + 1:, k] = 0.0
    return L, U

def DeltaLower(M):
    return np.tril(M)

def Psi_mats(alpha, beta, lam, p):
    """Psi from phi/chi recurrences at x = lam (book p.63/65)."""
    r = p + 1
    B = [1.0]; A = [1.0]
    for k in range(1, r):
        B.append(lam * B[-1] + beta[k - 1])
        A.append(lam * A[-1] + alpha[k - 1])
    ## phi(x) = [B_p ... B_1 1]^T ; derivatives w.r.t. x via sympy-free finite difference is unsafe --
    ## use exact polynomial derivatives: B_k(x) is a polynomial in x, build coefficient lists.
    import numpy.polynomial.polynomial as P
    Bp = [np.array([1.0])]
    for k in range(1, r):
        Bp.append(np.polyadd(np.polymul(Bp[-1], [1.0, 0.0]), [beta[k - 1]]))
    Ap = [np.array([1.0])]
    for k in range(1, r):
        Ap.append(np.polyadd(np.polymul(Ap[-1], [1.0, 0.0]), [alpha[k - 1]]))
    def dpoly(coef, order):
        c = coef
        for _ in range(order):
            c = np.polyder(c)
        return c
    Psi = np.zeros((r, r))
    for j in range(r):                       # column j: (1/(p-j)!) phi^(p-j)(lam)
        order = p - j
        col = []
        for i in range(r):                   # phi = [B_p, B_{p-1}, ..., B_1, 1]
            k = p - i
            coef = Bp[k] if k >= 1 else np.array([1.0])
            d = dpoly(coef, order)
            col.append(np.polyval(d, lam) / factorial(order) if len(d) else 0.0)
        Psi[:, j] = col
    Psinv = np.zeros((r, r))
    for i in range(r):                       # row i: (1/i!) chi^(i)(lam)
        row = []
        for j in range(r):                   # chi^T = [1, A_1, ..., A_p]
            coef = Ap[j] if j >= 1 else np.array([1.0])
            d = dpoly(coef, i)
            row.append(np.polyval(d, lam) / factorial(i) if len(d) else 0.0)
        Psinv[i, :] = row
    return Psi, Psinv

def construct(p, lam, c, alpha, beta, eps, Tmat=None):
    r = p + 1
    J, K = shifts(p)
    C = Cmat(c, p); F = Fmat(lam, p)
    bK = polyK(beta, K)
    e1 = np.eye(r)[:, 0:1]; ep1 = np.eye(r)[:, r - 1:r]
    Psi, _Psinv = Psi_mats(alpha, beta, lam, p)
    Omega = C @ (bK @ ep1 @ e1.T + K @ Psi)              # W = I
    OI = np.zeros((p, r));  OsI = np.zeros((r, p))
    for i in range(p):
        OI[i, i + 1] = 1.0; OsI[i + 1, i] = 1.0
    delta = np.zeros((r, 1))
    delta[0, 0] = eps + lam * M_closed(p, lam)
    for i in range(1, r):
        delta[i, 0] = N_closed(p + 1 - i, lam)
    Gamma = OsI @ OI @ F @ OsI @ OI + delta @ e1.T
    TT = np.eye(r) if Tmat is None else Tmat
    LT, UT = LUsplit(TT)
    OT = Omega @ TT
    LO, UO = LUsplit(OT)
    inner = UT @ np.linalg.inv(TT) @ Gamma @ TT @ np.linalg.inv(UO)
    Btil = LT @ DeltaLower(inner) @ np.linalg.inv(LO)
    return dict(delta=delta.ravel(), Omega=Omega, Gamma=Gamma, F=F, Btil=Btil, Psi=Psi)

if __name__ == '__main__':
    p, lam, eps = 2, 0.25, 0.0
    c = [0.25, 0.5, 1.0]
    beta = [0.0, 0.25, 1.0]              # peer: printed beta in DESCENDING powers -> b1=1/4? test both
    ## (2) N_k from JF - FJ against the closed form
    J, K = shifts(p); F = Fmat(lam, p)
    comm = J @ F - F @ J
    print('JF-FJ =\n', comm)
    print('closed N_1..N_p:', [N_closed(k, lam) for k in range(1, p + 1)])
    print('identity N_n = M_nn + lam M_(n-1)(n-1): max dev',
          max(abs(N_closed(n, lam) - (M_closed(n, lam) + lam * M_closed(n - 1, lam))) for n in range(1, 6)))
    ## (1) delta
    for bb in ([0.0, 0.25, 1.0], [1.0, 0.25, 0.0]):
        out = construct(p, lam, c, [0.0, 0.0, 0.0], bb, eps)
        print('beta=%s -> delta = %s   (target [1/64, 1/4, 1] = [0.015625, 0.25, 1.0])' % (bb, out['delta']))

    ## (3) end to end: Wright's printed p = 2, lambda = 1/4 method.
    ## beta(K) e_{p+1} = [b_p ... b_1 1] = [0, 1/4, 1] -- the PRINTED vector IS
    ## that product (K^{p+1} = 0, so b_{p+1} never enters), giving b_1 = 1/4,
    ## b_2 = 0; alpha enters only through Psi^-1, which (3.7.6) does not use.
    print()
    B_printed = np.array([[1/6, 1/2, 1/4], [0.0, 0.0, 1.0], [0.0, -2.0, 2.0]])
    for lab, bb in (('b=[1/4,0,0]', [0.25, 0.0, 0.0]), ('b=[0,1/4,0]', [0.0, 0.25, 0.0])):
        out = construct(p, lam, c, [0.0, 0.0, 0.0], bb, eps)
        Bt = out['Btil']
        print('%s  B~ =\n%s' % (lab, Bt))
        print('   max|B~ - B_printed| = %.3e' % np.max(np.abs(Bt - B_printed)))

    ## (3) CLOSED: B~ is the TRANSFORMED method (3.5.10), B = Psi B~ W^-1, W = I.
    ## Wright p.73: "A~ is strictly lower triangular and B~ is LOWER TRIANGULAR";
    ## Psi is unit UPPER triangular, so B = Psi B~ is full -- which is exactly
    ## why B~ from (3.7.8) is triangular and the printed B is not.
    print()
    out = construct(p, lam, c, [0.0, 0.0, 0.0], [0.25, 0.0, 0.0], eps)
    Psi = out['Psi']; Btil = out['Btil']
    print('Psi =\n', Psi)
    print('B = Psi B~ =\n', Psi @ Btil)
    print('printed B =\n', B_printed)
    print('max|Psi B~ - B_printed| = %.3e' % np.max(np.abs(Psi @ Btil - B_printed)))
