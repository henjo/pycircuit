"""B7's gate: does a DEFECT-CORRECTION global error estimate see the WARPING that a per-step LTE
cannot?  (Sickenberger, Weinmueller & Winkler Part I: IDeC estimates LOCAL AND GLOBAL errors; the
neighbouring problem's exact solution is known by construction.)

Fixture: A10's van der Pol, C = L = 1, mu = 1/(2 pi Q), Q = 1e4 (nearly-harmonic, where trap's
period error at 400 pts is ~21 ppm -- the warping the LTE was blind to).  Method under test: trap
(order 2), ODE form  x' = f(x)  with x = (v, iL):  v' = -(iL + i_B(v)),  iL' = v.

Steps:
 1. x_h, T_h  = trap's PERIODIC solution at N points per period (numpy shooting in (v0, T)).
 2. p(t)      = periodic cubic-spline interpolant of x_h  (degree 3 >= method order 2).
 3. d(t)      = p'(t) - f(p(t))                       -- the defect, T_h-periodic, O(h^2).
 4. y' = f(y) + d(t), y(0) = p(0): the NEIGHBOURING problem, whose exact solution IS p(t).
    Integrate it with the SAME trap at the SAME h for K periods (a transient, B7's setting).
 5. Phase drift of y_h against p:  tau(t) = <p'(t) . (y_h(t) - p(t))> / <|p'|^2>  per period;
    the slope of tau against period count is the ESTIMATED period error  delta' (in time units).
 6. Compare with the TRUE period error  delta = T_h - T_ref  (T_ref from radau at 3200 points).

PREDICTION NAMED BEFORE RUNNING: if IDeC sees warping, delta'/delta is within ~20 % of 1 at each N
and both scale as h^2.  If delta' is ~0 or off by 10x, this construction does not see it.
Control: the same drift measured on the ORIGINAL problem's transient against the exact-period
reference orbit must give delta itself (that is what T_h - T_ref means).
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import solve

Q = 1e4; mu = 1.0 / (2 * np.pi * Q)

def f(x):
    v, iL = x
    return np.array([-(iL + mu * (v - v ** 3 / 3.0)), v])
def J(x):
    v, iL = x
    return np.array([[-mu * (1 - v * v), -1.0], [1.0, 0.0]])

## --- Runge-Kutta stepper from a tableau (trap = CN; radau = coded Radau IIA(3)) on x' = f(x) + g(t)
def rk_step(A, b, c, x, t, h, g=None):
    s = len(b); Y = np.tile(x, (s, 1))
    for _ in range(50):                      # Newton on the stage system
        F = np.zeros((s, 2)); Jb = np.zeros((2 * s, 2 * s))
        for i in range(s):
            acc = np.zeros(2)
            for j in range(s):
                gj = g(t + c[j] * h) if g is not None else 0.0
                acc += A[i, j] * (f(Y[j]) + gj)
            F[i] = Y[i] - x - h * acc
            for j in range(s):
                Jb[2*i:2*i+2, 2*j:2*j+2] = -h * A[i, j] * J(Y[j])
            Jb[2*i:2*i+2, 2*i:2*i+2] += np.eye(2)
        dY = solve(Jb, -F.ravel()).reshape(s, 2); Y += dY
        if np.abs(dY).max() < 1e-13: break
    acc = np.zeros(2)
    for i in range(s):
        gi = g(t + c[i] * h) if g is not None else 0.0
        acc += b[i] * (f(Y[i]) + gi)
    return x + h * acc

TRAP = (np.array([[0, 0], [0.5, 0.5]]), np.array([0.5, 0.5]), np.array([0.0, 1.0]))
S6 = np.sqrt(6.0)
RADAU = (np.array([[11/45-7*S6/360, 37/225-169*S6/1800, -2/225+S6/75],
                   [37/225+169*S6/1800, 11/45+7*S6/360, -2/225-S6/75],
                   [4/9-S6/36, 4/9+S6/36, 1/9]]),
         np.array([4/9-S6/36, 4/9+S6/36, 1/9]), np.array([2/5-S6/10, 2/5+S6/10, 1.0]))

def integrate(tab, x0, T, N, K=1, g=None):
    """K periods of length T at N steps per period; returns times and states (K*N+1)."""
    A, b, c = tab; h = T / N; xs = [x0.copy()]; x = x0.copy()
    for n in range(K * N):
        x = rk_step(A, b, c, x, n * h, h, g); xs.append(x.copy())
    return np.arange(K * N + 1) * h, np.array(xs)

def periodic(tab, N, T0, v0=2.0):
    """Shooting for (v0, T): x(T; (v0, 0)) - (v0, 0) = 0, phase condition iL(0) = 0."""
    u = np.array([v0, T0])
    for it in range(30):
        def Fres(u):
            _, xs = integrate(tab, np.array([u[0], 0.0]), u[1], N); return xs[-1] - np.array([u[0], 0.0])
        F0 = Fres(u); Jm = np.zeros((2, 2)); eps = 1e-7
        for k in range(2):
            d = np.zeros(2); d[k] = eps * max(1.0, abs(u[k])); Jm[:, k] = (Fres(u + d) - F0) / d[k]
        du = solve(Jm, -F0); u = u + du
        if np.abs(du).max() < 1e-12: break
    _, xs = integrate(tab, np.array([u[0], 0.0]), u[1], N)
    return u[0], u[1], xs

T_ref = 6.283185307279     # radau IIA(3) at 3200 pts, computed by this script 2026-09-08 (first run)
print('T_ref = %.12f (cached from the radau@3200 run)' % T_ref, flush=True)
from scipy.interpolate import make_interp_spline

def estimate(N, tab_basic, tab_neigh, degree=3, K=20):
    v0, T_h, xs = periodic(tab_basic, N, T_ref)
    delta = T_h - T_ref
    xs = xs.copy(); xs[-1] = xs[0]                     # close the orbit (shooting residual ~1e-12)
    tg = np.arange(N + 1) * (T_h / N)
    if degree == 3:
        p = CubicSpline(tg, xs, bc_type='periodic')
    else:                                               # linear: periodic by construction of the data
        p = make_interp_spline(tg, xs, k=1)
    dp = p.derivative()
    per = lambda t: t % T_h
    d = lambda t: dp(per(t)) - f(p(per(t)))
    tt, ys = integrate(tab_neigh, xs[0].copy(), T_h, N, K=K, g=d)
    pt = p(per(tt)); dpt = dp(per(tt)); e = ys - pt
    tau = np.array([np.sum(dpt[k*N:(k+1)*N] * e[k*N:(k+1)*N]) / np.sum(dpt[k*N:(k+1)*N] ** 2) for k in range(K)])
    slope = np.polyfit(np.arange(K), tau, 1)[0]
    ## SIGN: tau > 0 means y is AHEAD of p; a LONGER period makes y fall BEHIND, so the estimated
    ## period error is delta' = -slope.  (The first run printed ratio -1.000: the convention, not the estimate.)
    return delta, -slope

print('== trap basic / trap neighbouring / cubic interpolant  (the gate)')
print('%-6s %12s %12s %8s' % ('N', 'delta(ppm)', "delta'(ppm)", 'ratio'))
for N in (100, 200, 400, 800):
    delta, dprime = estimate(N, TRAP, TRAP, 3)
    print('%-6d %12.4f %12.4f %8.4f' % (N, delta / T_ref * 1e6, dprime / T_ref * 1e6, dprime / delta), flush=True)
print('== CONTROL (a): neighbouring problem solved by RADAU -- the drift must vanish (ratio ~ 0) if it is the METHOD\'s error')
for N in (100, 200):
    delta, dprime = estimate(N, TRAP, RADAU, 3)
    print('%-6d %12.4f %12.4f %8.4f' % (N, delta / T_ref * 1e6, dprime / T_ref * 1e6, dprime / delta), flush=True)
print('== CONTROL (b): LINEAR interpolant (degree 1 < method order 2) -- the peer\'s interpolant-order limit, measured')
for N in (100, 200, 400):
    delta, dprime = estimate(N, TRAP, TRAP, 1)
    print('%-6d %12.4f %12.4f %8.4f' % (N, delta / T_ref * 1e6, dprime / T_ref * 1e6, dprime / delta), flush=True)
