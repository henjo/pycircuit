"""Nordsieck-form GLM tools: build U, V for stage order q = p from (A, c, B) (Butcher: U = C - A C K,
V = E - B C K with C = [c^k/k!], K the shift, E = exp(K)); check the order by exactness on polynomials
(the method applied to y' = f(t) with y = t^k must be exact for k <= p, stages included); the stability
function M(z) = V + z B (I - zA)^-1 U; M_inf = V - B A^-1 U (nilpotent for Voigtmann's Theorem 9.5);
and a numerical constructor for a diagonally implicit, stiffly accurate, L-stable method with s = r = p + 1."""
import numpy as np
from math import factorial


def nordsieck_matrices(p):
    r = p + 1
    K = np.zeros((r, r))
    for i in range(r - 1):
        K[i, i + 1] = 1.0
    E = np.zeros((r, r))              # exp(K): E[i, j] = 1/(j-i)! for j >= i
    for i in range(r):
        for j in range(i, r):
            E[i, j] = 1.0 / factorial(j - i)
    return K, E


def cmat(c, p):
    c = np.asarray(c, float)
    return np.array([[ck ** k / factorial(k) for k in range(p + 1)] for ck in c])   # s x (p+1)


def uv_from(A, c, B, p):
    K, E = nordsieck_matrices(p)
    C = cmat(c, p)
    U = C - A @ C @ K
    V = E - B @ C @ K
    return U, V


def step(A, U, B, V, c, y_in, h, f, t):
    """One GLM step on y' = f(t, y): returns (Y stages, y_out).  Stages solved by simple fixed-point for
    the polynomial exactness check (f independent of y there)."""
    s = A.shape[0]
    F = np.zeros(s); Y = np.zeros(s)
    for i in range(s):
        Y[i] = U[i] @ y_in + h * sum(A[i, j] * F[j] for j in range(i))
        # implicit diagonal for f(t) only: F_i = f(t + c_i h) does not depend on Y
        F[i] = f(t + c[i] * h)
        Y[i] += h * A[i, i] * F[i]
    y_out = V @ y_in + h * (B @ F)
    return Y, y_out


def order_by_polynomials(A, U, B, V, c, p, h=0.37, t0=0.21):
    """Max error of one step on y = t^k, k = 0..p+1: exact (~1e-15) for k <= p means order/stage order p;
    the k = p+1 error is the leading term (must NOT vanish or the check saw nothing)."""
    out = []
    for k in range(p + 2):
        f = lambda t: k * t ** (k - 1) if k > 0 else 0.0 * t
        ## Butcher's Nordsieck convention: y_j = h^j y^(j) WITHOUT the 1/j! (so E = exp(K) is the Taylor shift)
        y_in = np.array([h ** j * (factorial(k) / factorial(k - j) * t0 ** (k - j) if j <= k else 0.0) for j in range(p + 1)])
        Y, y_out = step(A, U, B, V, c, y_in, h, f, t0)
        t1 = t0 + h
        y_ex = np.array([h ** j * (factorial(k) / factorial(k - j) * t1 ** (k - j) if j <= k else 0.0) for j in range(p + 1)])
        Y_ex = np.array([(t0 + ci * h) ** k for ci in c])
        out.append((k, float(np.max(np.abs(y_out - y_ex))), float(np.max(np.abs(Y - Y_ex)))))
    return out


def M(z, A, U, B, V):
    s = A.shape[0]
    return V + z * B @ np.linalg.solve(np.eye(s) - z * A, U)


def M_inf(A, U, B, V):
    return V - B @ np.linalg.solve(A, U)


def stability_report(A, U, B, V, grid=None):
    if grid is None:
        grid = [complex(x, y) for x in (-1e-3, -0.1, -1, -10, -100, -1e4) for y in (0, 0.5, 2, 10, 100)]
    rho = max(np.max(np.abs(np.linalg.eigvals(M(z, A, U, B, V)))) for z in grid)
    Mi = M_inf(A, U, B, V)
    r = Mi.shape[0]
    nil = float(np.max(np.abs(np.linalg.matrix_power(Mi, r))))
    rhoV = np.sort(np.abs(np.linalg.eigvals(V)))[::-1]
    return dict(rho_lhp=float(rho), nilpotency_residual=nil, eig_V=rhoV, rho_Minf=float(np.max(np.abs(np.linalg.eigvals(Mi)))))


def construct(p, lam=None, seed=0, tries=40):
    """Search: A lower triangular with diagonal lam (s = p+1 stages), c free with c[-1] = 1, B free (r x s),
    U,V from the order conditions; constraints: stiff accuracy (B[0] = A[-1], V[0] = U[-1] -> the latter is
    automatic when B[0] = A[-1] and c[-1] = 1? no -- imposed), eig(V) = {1, 0..}, M_inf nilpotent, rho on the
    LHP grid <= 1."""
    from scipy.optimize import least_squares
    s = r = p + 1
    rng = np.random.default_rng(seed)
    K, E = nordsieck_matrices(p)
    nA = s * (s - 1) // 2
    def unpack(x):
        A = np.zeros((s, s)); k = 0
        for i in range(s):
            for j in range(i):
                A[i, j] = x[k]; k += 1
        l = x[k] if lam is None else lam
        if lam is None: k += 1
        A += l * np.eye(s)
        c = np.concatenate((x[k:k + s - 1], [1.0])); k += s - 1
        B = x[k:k + r * s].reshape(r, s); k += r * s
        return A, c, B
    def resid(x):
        A, c, B = unpack(x)
        U, V = uv_from(A, c, B, p)
        res = []
        res += list(B[0] - A[-1])                       # stiff accuracy: y_out[0] = last stage (with c_s = 1, U[-1] = V[0] follows from order conditions? impose too)
        res += list(V[0] - U[-1])
        ev = np.sort(np.abs(np.linalg.eigvals(V)))[::-1]
        res += [ev[0] - 1.0] + list(ev[1:] * 3.0)         # V spectrum {1, 0, ..}
        Mi = V - B @ np.linalg.solve(A, U)
        res += list(np.linalg.matrix_power(Mi, r).ravel() * 3.0)   # nilpotent
        res += list(np.abs(np.linalg.eigvals(Mi)) * 3.0)
        for z in (-1.0, -10.0, -100.0, complex(-1, 3), complex(-0.1, 1), complex(-10, 30)):
            rho = np.max(np.abs(np.linalg.eigvals(M(z, A, U, B, V))))
            res.append(max(rho - 0.999, 0.0) * 10.0)
        res += list(0.02 * (c[:-1] - np.linspace(0.2, 0.8, s - 1)))   # mild preference for spread abscissae in (0, 1)
        return np.array(res)
    best = None
    for t in range(tries):
        x0 = rng.normal(0, 0.5, nA + (0 if lam is not None else 1) + (s - 1) + r * s)
        if lam is None: x0[nA] = 0.3
        x0[nA + (0 if lam is not None else 1):nA + (0 if lam is not None else 1) + s - 1] = np.linspace(0.2, 0.8, s - 1)
        try:
            sol = least_squares(resid, x0, max_nfev=4000)
        except Exception:
            continue
        cost = float(np.sum(sol.fun[:-(s - 1)] ** 2))
        if best is None or cost < best[0]:
            best = (cost, sol.x)
    A, c, B = unpack(best[1]); U, V = uv_from(A, c, B, p)
    return A, U, B, V, c, best[0]


if __name__ == '__main__':
    import sys
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    A, U, B, V, c, cost = construct(p, lam=float(sys.argv[2]) if len(sys.argv) > 2 else None)
    np.set_printoptions(precision=6, suppress=True, linewidth=150)
    print('p = %d, construction cost %.2e' % (p, cost)); print('A =\n', A); print('c =', c); print('B =\n', B); print('U =\n', U); print('V =\n', V)
    print('polynomial exactness (k, |y_out err|, |stage err|):'); [print('  ', o) for o in order_by_polynomials(A, U, B, V, c, p)]
    print('stability:', stability_report(A, U, B, V))
