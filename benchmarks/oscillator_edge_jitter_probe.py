"""B0 of the #17 plan: does `PAC.oscillator_edge_jitter`'s k-cycle law hold?

The shipped law is `k_cycle_bound^2 = c k T + 2A`, `A = e' Pi P_j Pi' e / s^2`
(the bounded covariance at the edge, obliquely projected), called "the
large-k form with rho_k -> 0".  The EXACT linear k-lag variance of the edge
time, from the bounded split `K(t_j + nT) = P_j + n G_j` (`G_j = d u_j u_j'`,
`M_j u_j = u_j`), is

    Var_k = e' (2 P_j + k G_j - M_j^k P_j - P_j M_j^k') e / s^2,

whose large-k intercept is `2 (A + X)`, `X = e' Pi P_j (I - Pi)' e / s^2` --
the phase-transverse cross term the shipped law drops.

PREDICTIONS (a design pass, 2026-09-29, read-only, radau 240 on the A11
fixture), pinned BEFORE this probe ran:
    Var_1 exact = 1.5274e-6,  k_cycle_bound^2(k=1) = 1.6969e-6,
    X / A = -0.150,  2(A + X) = 0.9865e-6,  2A = 1.1608e-6.

`mc` measures `Var(t_{n+k} - t_n)` from NOISY radau transients at the PSS's
own step (the recorded lesson: an Euler Monte Carlo against a gear slope cost
5.7 %), noise piecewise-constant per step with variance `S/(2h)` (the
convention pinned against kT/C in `_a8_mc_edge_scatter`).

MEASURED 2026-09-29 (this probe, committed):
  * `exact radau 240` reproduces every prediction to its stated digits:
    Var_1 = 1.52740e-6, bound^2 = 1.69690e-6, X/A = -0.1502,
    2(A + X) = 9.86462e-7, 2A = 1.16079e-6.
  * `mc`, 8 seeds x 720 periods (5760 crossings), MC against the exact law
    and the shipped bound:
        k   MC                  exact            bound
        1   1.5001e-6 +- 2.2e-8 1.5274e-6 -1.2s  1.6969e-6 -8.9s
        2   2.0227e-6 +- 3.0e-8 2.0683e-6 -1.5s  2.2330e-6 -7.1s
        4   3.058e-6  +- 5.5e-8 3.150e-6  -1.7s  3.305e-6  -4.5s
        8   5.05e-6   +- 2.0e-7 5.31e-6   -1.3s  5.45e-6   -2.0s
    At k = 1, MC - c k T = 0.964e-6 +- 0.022e-6: 2(A + X) = 0.986e-6 at -1s,
    the shipped 2A = 1.161e-6 at -9s.  So the exact linear law holds and the
    shipped intercept 2A does not; the additive variance is
    `e' Pi P_j e / s^2 = A + X` (the one-sided projection), not
    `e' Pi P_j Pi' e / s^2 = A`.  The exact law sits ~2 % above the MC at
    every k (each 1.2-1.7 sigma, one set of seeds, so correlated): noted,
    not significant.  ~9 min per seed at 720 periods on a loaded box.

RE-RUN 2026-09-29 (step 2 of the edge-jitter plan), after the law became
each node's own blended to the instant (1.5182e-6 at k = 1, grid-converged
to 8e-5 between 240 and 480 points; the nearest-node law above moved 1.2e-2):
the Monte Carlo on the PSS's OWN grid (h = T/239 -- the first ran T/240)
with a QUADRATIC crossing estimator beside the linear one, seeds 11-18 x
720 periods, 14-20 min per seed on a loaded box:
        k   MC (quadratic)        per-node law       nearest-node law
        1   1.5049e-6 +- 1.8e-8   1.5182e-6 +0.75s   1.5274e-6 +1.27s
        2   2.1007e-6 +- 4.7e-8   2.0543e-6 -0.99s   2.0683e-6 -0.69s
        4   3.2000e-6 +- 1.1e-7   3.1265e-6 -0.64s   3.1501e-6 -0.44s
        8   5.2888e-6 +- 3.4e-7   5.2710e-6 -0.05s   5.3137e-6 +0.07s
  (the linear estimator within 0.3 % of the quadratic).  Both laws sit
  within 1.3 sigma: the Monte Carlo cannot separate them (0.6 % apart, a
  1.2 % standard error) -- the grid convergence does.

Usage:  python oscillator_edge_jitter_probe.py exact [radau|gear] [npts]
        python oscillator_edge_jitter_probe.py mc SEED NPER [npts] [save.npy]
"""
import sys
import time
import warnings

import numpy as np

from pycircuit.circuit import IS, PSS, BSource, C, L, R, SubCircuit, circuit, gnd
from pycircuit.circuit.shooting import PAC

T_SEED = 6.664052486


def a11_chain(psd_tank=1e-6, psd_buf=1e-6, nstage=3, cb=0.5, noise_fns=None):
    """`_a11_osc_chain` of test_shooting_oscillator_covariance.py; with
    `noise_fns` the noise sources are ideal TIME-DOMAIN currents instead
    (for the Monte Carlo), noiseless otherwise."""
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
    if noise_fns is None and psd_tank > 0:
        c['n'] = IS('v', gnd, i=0.0, noisePSD=psd_tank)
    prev = 'v'
    for k in range(nstage):
        nd = f'o{k}'
        c.add_node(nd)
        c[f'B{k}'] = BSource(prev, gnd, nd, gnd,
                               i_func=lambda u: 1.0 * np.tanh(2.0 * u))
        c[f'R{k}'] = R(nd, gnd, r=1.0, noisy=False)
        c[f'C{k}'] = C(nd, gnd, c=cb)
        if noise_fns is None and psd_buf > 0:
            c[f'nb{k}'] = IS(nd, gnd, i=0.0, noisePSD=psd_buf)
        prev = nd
    if noise_fns is not None:
        for name, (node, f) in noise_fns.items():
            c[name] = IS(node, gnd, i=0.0)
            c[name].function.f = f
    return c


def solved(method='radau', npts=240):
    circuit.default_toolkit = circuit.numeric
    cir = a11_chain()
    m = cir.n - 1
    pss = PSS(cir, method=method, reltol=1e-12)
    x0 = np.zeros(m)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T_SEED, timestep=T_SEED / npts, x0=x0, maxiterations=80)
    assert pss.converged
    full = cir.get_node_index('o2')
    red = full if full < cir.get_node_index(gnd) else full - 1
    Xw = np.asarray(pss.waveform[1], dtype=float)
    grid = np.asarray(pss.factored_period().times, dtype=float)[:Xw.shape[1]]
    v = Xw[full]
    mid = 0.5 * (v.max() + v.min())
    j = next(k for k in range(2, len(v) - 2)
             if (v[k - 1] - mid) < 0 <= (v[k] - mid))
    tc = grid[j - 1] + (mid - v[j - 1]) / (v[j] - v[j - 1]) * (grid[j] - grid[j - 1])
    return cir, pss, red, float(tc), float(mid)


def exact(method='radau', npts=240, ks=(1, 2, 4, 8, 16, 64)):
    cir, pss, red, tc, _mid = solved(method, npts)
    pac = PAC(cir, toolkit=circuit.numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = pac.oscillator_edge_jitter(pss, red, tc)
        host = pss._lyapunov_host()
        _K, info = pac.oscillator_covariance(host, samples=True)
        As, _Qs, _K1, _M, m, n = pac._lyapunov_pieces(host, 'probe')
        _v0, pinfo = host.ppv()
    assert n == m, 'the exact law is written for a one-step host (n == m)'
    from pycircuit.circuit.shooting._numerics import edge_slope
    fp = host.factored_period()
    times = np.asarray(fp.times, dtype=float)
    row = np.asarray(pac._output_waveform_row(host, red), dtype=float)
    nt = min(len(times), len(row))
    s, j = edge_slope(times[:nt], row[:nt], tc)
    N = len(As)
    jj = min(j, len(info['samples']) - 1)
    P = np.asarray(info['samples'][jj], dtype=float)
    G = np.asarray(info['growth_samples'][jj], dtype=float)
    ## the one-period map from node jj: A_jj first, around the period
    Mj = np.eye(m)
    for i in list(range(jj, N)) + list(range(jj)):
        Mj = np.asarray(As[i], dtype=float) @ Mj
    e = np.zeros(m)
    e[red] = 1.0
    w, U = np.linalg.eigh(G)
    u = U[:, int(np.argmax(w))] * np.sqrt(float(w.max()))
    vj = np.asarray(pinfo['samples'][jj], dtype=float)[:m]
    Pi = np.eye(m) - np.outer(u, vj) / float(vj @ u)
    A = float(e @ Pi @ P @ Pi.T @ e) / s ** 2
    X = float(e @ Pi @ P @ (np.eye(m) - Pi).T @ e) / s ** 2
    T = float(host.period)
    print(f"{method} npts={npts}  node j={jj}  s={s:.6f}  shipped: "
          f"A={r['A']:.5e} c={r['c']:.5e}  probe A={A:.5e}")
    print(f'X/A = {X / A:.4f}   2A = {2 * A:.5e}   2(A+X) = {2 * (A + X):.5e}')
    Mk = np.eye(m)
    for k in range(1, max(ks) + 1):
        Mk = Mj @ Mk
        if k in ks:
            V = float(e @ (2 * P + k * G - Mk @ P - P @ Mk.T) @ e) / s ** 2
            bound = r['c'] * k * T + 2 * r['A']
            print(f"k={k:3d}  exact Var_k={V:.5e}  bound^2={bound:.5e}  "
                  f"exact-ckT={V - r['c'] * k * T:.5e}  ratio={V / bound:.5f}")
    return r


def crossings(v, h, mid, t_start, quadratic=True):
    """The rising crossings of `mid` after `t_start` in `v` sampled at `k h`:
    linear between the two samples, or the root of the quadratic through
    the three samples nearest the crossing."""
    out = []
    for k in range(2, len(v) - 1):
        if (k - 1) * h < t_start or not ((v[k - 1] - mid) < 0 <= (v[k] - mid)):
            continue
        frac = (mid - v[k - 1]) / (v[k] - v[k - 1])
        if quadratic:
            idx = [k - 2, k - 1, k] if frac < 0.5 else [k - 1, k, k + 1]
            q = np.polyfit(np.asarray(idx, dtype=float) - (k - 1), v[idx] - mid, 2)
            r = np.roots(q)
            r = r[np.abs(r.imag) < 1e-12].real
            if len(r):
                frac = float(r[np.argmin(np.abs(r - frac))])
        out.append((k - 1 + frac) * h)
    return np.asarray(out)


def mc(seed, nper, npts=240, burn=20, save=None):
    """`Var(t_{n+k} - t_n)` from one noisy radau transient ON THE PSS'S OWN
    GRID (`h = T/N`, N its steps: the PSS's discrete orbit, crossing at the
    same fraction of a step every period up to the noise), by a linear and
    a quadratic crossing estimator.  Until 2026-09-29 it ran `h = T/npts`
    (240 against the PSS's 239) and the linear estimator alone."""
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    from pycircuit.circuit.transient import Transient
    _cir, pss, _red, _tc, mid = solved('radau', npts)
    T = float(pss.period)
    grid = np.asarray(pss.factored_period().times, dtype=float)
    N = len(grid) - 1
    h = T / N
    assert np.allclose(np.diff(grid), h, rtol=1e-9, atol=0.0), 'a non-uniform PSS grid'
    nst = int((nper + burn) * N) + 16
    rng = np.random.default_rng(seed)
    S = 1e-6
    draws = {name: rng.normal(0.0, np.sqrt(S / (2.0 * h)), size=nst)
             for name in ('n', 'nb0', 'nb1', 'nb2')}
    ## piecewise-constant over (t_i, t_i + h]: every stage of step i reads draw i
    def f_of(d):
        return lambda t, d=d: float(d[max(int(np.ceil(t / h - 1e-9)) - 1, 0)])
    fns = {'n': ('v', f_of(draws['n']))}
    for k in range(3):
        fns[f'nb{k}'] = (f'o{k}', f_of(draws[f'nb{k}']))
    cir = a11_chain(noise_fns=fns)
    x0 = np.asarray(pss.waveform[1], dtype=float)[:, 0]
    tran = Transient(cir, integrator=RadauIIA3Integrator())
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sol = tran.solve(tend=(nper + burn) * T, timestep=h, x0=x0,
                         fixed_timestep=True)
    v = np.asarray(sol.v('o2', gnd), dtype=float)
    if save:
        np.save(save, v)
    out = {'N': N, 'seconds': time.time() - t0}
    for est in ('lin', 'quad'):
        cross = crossings(v, h, mid, burn * T, quadratic=(est == 'quad'))
        out[est] = {'n': len(cross)}
        for k in (1, 2, 4, 8):
            dk = cross[k:] - cross[:-k]
            out[est][k] = float(np.var(dk - np.mean(dk), ddof=1))
    print(seed, out)
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'exact':
        exact(sys.argv[2] if len(sys.argv) > 2 else 'radau',
              int(sys.argv[3]) if len(sys.argv) > 3 else 240)
    else:
        mc(int(sys.argv[2]), int(sys.argv[3]),
           int(sys.argv[4]) if len(sys.argv) > 4 else 240,
           save=sys.argv[5] if len(sys.argv) > 5 else None)
