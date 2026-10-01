"""The normalised oscillator lineshape when the phase noise is COLOURED.

For a Gaussian excess phase the `i`-th harmonic's normalised spectrum is

    S(f) = 2 int_0^inf exp(-D(tau)/2) cos(2 pi f tau) dtau,
    D(tau) = 4 int_0^inf S_phi(nu) (1 - cos 2 pi nu tau) dnu,

with `S_phi(nu) = i^2 f0^2 c(nu) / nu^2` the phase PSD of `PAC.phase_psd`
(the lineshape's own far skirt).  The WHITE part of `c` gives `D = 2 a tau`,
`a = 2 pi^2 i^2 f0^2 c_w`, and the Lorentzian in closed form; the
COLOURED part `c_c(nu) = c(nu) - c_w` is integrated over `[fmin, fmax]`:

    S(f) = e^{-D_inf/2} L_w(f) + 2 int_0^inf e^{-a tau}
           (e^{-D_c(tau)/2} - e^{-D_inf/2}) cos(2 pi f tau) dtau

where `D_inf = D_c(inf)` is finite because of the band.  A white-only
circuit has `D_c = 0` and gets `L_w` exactly.

Numerics (gated against closed forms and an independent mpmath
reference, 2026-09-26, in `doc/pss_log_260902.md`):
  * `c_c` is a piecewise POWER LAW through nodes refined until the
    interpolant meets `c_c` at every geometric midpoint (`NODE_TOL`), so
    a power-law colour is exact and its moments are closed forms;
  * `D_c(tau)`: below `2 pi nu tau = 2` the cosine series on the pieces
    (exact, no cancellation), above it the closed moment minus QUADPACK's
    oscillatory rule (QAWO), per decade;
  * `D_c` on a log tau grid (`TAU_PER_DECADE`), a cubic spline of
    `log D_c` in `log tau`, the quadratic asymptote below and the
    band edge's ringing above;
  * the transform by QUADPACK's cosine rules (QAWO per decade, QAWF on
    the tail).

A SIGNED correction to `c` (the frequency-aware PPV's `c_fa - c_dc`,
`PAC._fa_lineshape`) cannot ride on the power-law pieces, which
cannot pass a zero.  It is tabulated (`SignedTable`) and its structure
function taken by QUADPACK (`correction_structure`, `2 sin^2` below
`beta nu = 2`), splined in its VALUE beside `D_c`, and added.  `D_inf` may
then be NEGATIVE (white noise taken away over a band holding the core), and
past `SPLIT_GAIN` the transform integrates `exp(-D/2)` itself, with the
Lorentzian a closed-form tail.  With a correction the TOTAL `D` is
splined (the white part and a correction that removes it cancel), and past
the band the white part is held at its corrected level (`ConstantTail`).
Gated against a closed form: <= 4.8e-7.
Far from the carrier the transform is a cancellation, and each offset
takes it or the SECOND-order skirt (`second_order_skirt`), whichever
`handover` estimates the smaller error.  Against the mpmath reference
(`benchmarks/lineshape_reference.py`) the worst is 3.4e-7 on a real
oscillator's line.

The correction's `rho = c_fa/c_dc - 1` comes from bordered solves.
`RationalRho` fits it as one rational from ~25 of them; `LogChebyshev` is
the fallback, taking ~65 to ~139.  Both hold the tolerance relative to
`1 + rho`.

History: `doc/shooting_history.md`, `_lineshape`.
"""
import math
import warnings

import numpy as np
from scipy import integrate, interpolate, special

#: the power-law interpolant meets `c_c` to this, relative, at every
#: geometric midpoint
NODE_TOL = 1e-6
#: initial nodes per decade of `nu`
NODES_PER_DECADE = 20
#: `D_c` samples per decade of `tau` (10 left the far skirt at 1.8e-4)
TAU_PER_DECADE = 40
#: the node refinement stops here, with a warning
MAX_NODES = 4000
#: a SIGNED correction's `rho` is resolved until its Chebyshev tail is below
#: this (`LogChebyshev`)
CHEB_TOL = 1e-6
#: and its degree doubles (nested points) up to this
CHEB_MAX = 128
#: `RationalRho` gives up (the caller falls back to `LogChebyshev`) past
#: this many samples
RATIONAL_BUDGET = 40
#: `RationalRho`'s candidate grid, points per decade of `nu`
RATIONAL_PER_DECADE = 20
#: a signed correction is tabulated at this many points per decade of `nu`
#: for the quadrature (`SignedTable`)
TABLE_PER_DECADE = 2000
#: past this gain `exp(-D_inf/2)` the transform integrates `exp(-D/2)`
#: itself (`ColouredLineshape.__call__`): the subtracted form amplifies
#: QUADPACK's 1.5e-8 by the gain, and below it is the (slightly) better one
#: -- 2.1e-7 against 3.5e-7 at gain 19 on the closed form
SPLIT_GAIN = 50.0
#: the colour band's lower edge `nu0` rings in `D(tau)`; past `beta nu0 =
#: EDGE_SWITCH` that ringing is carried exactly (`ColouredLineshape`)
EDGE_SWITCH = 2.0
#: and past the lag where its amplitude `|gamma|/2` falls below this, the
#: transform takes it by its harmonics of `nu0` ...
EDGE_HARMONIC_R = 0.02
#: ... at most this many
EDGE_MAX_HARMONICS = 40
#: with no correction the integrand is `e^{-a tau}` small; past `EDGE_DECAY /
#: a` it is taken as zero
EDGE_DECAY = 80.0


def _ratio(E, lr):
    """`(exp(E lr) - 1) / E`, stable at `E -> 0`."""
    el = E * lr
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        return np.where(np.abs(el) < 1e-8, lr * (1.0 + 0.5 * el),
                        np.expm1(el) / np.where(E == 0.0, 1.0, E))


class PowerPieces:
    """`c(nu)`: a piecewise power law through `(nus, cs)`, `cs > 0`."""

    def __init__(self, nus, cs):
        self.nu = np.asarray(nus, dtype=float)
        self.c = np.maximum(np.asarray(cs, dtype=float), 1e-300)
        self.p = np.diff(np.log(self.c)) / np.diff(np.log(self.nu))

    def __call__(self, v):
        v = np.asarray(v, dtype=float)
        k = np.clip(np.searchsorted(self.nu, v, side='right') - 1,
                    0, len(self.p) - 1)
        return self.c[k] * (v / self.nu[k]) ** self.p[k]

    def _pieces(self, a, b):
        lo = np.clip(a, self.nu[:-1], self.nu[1:])
        hi = np.clip(b, self.nu[:-1], self.nu[1:])
        ok = hi > lo
        return (lo[ok], hi[ok], self.c[:-1][ok], self.nu[:-1][ok],
                self.p[ok])

    def moment(self, e, a, b):
        """`int_a^b c(nu) nu^e dnu`, exact on the pieces."""
        if b <= a:
            return 0.0
        lo, hi, c, nk, p = self._pieces(a, b)
        return float(np.sum(c * (lo / nk) ** p * lo ** (e + 1.0)
                            * _ratio(p + e + 1.0, np.log(hi / lo))))

    def cosine_series(self, beta, a, b, mmax=40):
        """`int_a^b c(nu) (1 - cos beta nu) / nu^2 dnu` by the cosine
        series, for `beta b <~ 2`; scaled by `(beta lo)^{2m}` so no power
        overflows."""
        if b <= a:
            return 0.0
        lo, hi, c, nk, p = self._pieces(a, b)
        lr = np.log(hi / lo)
        base = c * (lo / nk) ** p / lo
        tot, first = 0.0, None
        for m in range(1, mmax):
            t = ((-1.0) ** (m + 1) / math.factorial(2 * m)
                 * float(np.sum(base * (beta * lo) ** (2 * m)
                                * _ratio(p + 2.0 * m - 1.0, lr))))
            tot += t
            first = abs(t) if first is None else first
            if abs(t) <= 1e-18 * max(first, 1e-300):
                break
        return tot


def refine(cfun, nu0, nuN, per_decade=NODES_PER_DECADE, tol=NODE_TOL,
           maxpts=MAX_NODES):
    """`PowerPieces` for `cfun` (vectorised, positive) on `[nu0, nuN]`,
    each interval halved (geometrically) until the interpolant meets
    `cfun` at its midpoint to `tol`.  Returns `(pieces, converged)`."""
    n = max(int(np.ceil(per_decade * np.log10(nuN / nu0))), 2)
    nus = np.geomspace(nu0, nuN, n + 1)
    cs = np.asarray(cfun(nus), dtype=float)
    converged = False
    while True:
        mid = np.sqrt(nus[:-1] * nus[1:])
        cm = np.asarray(cfun(mid), dtype=float)
        ci = np.sqrt(np.maximum(cs[:-1], 1e-300) * np.maximum(cs[1:], 1e-300))
        bad = np.abs(cm / np.maximum(ci, 1e-300) - 1.0) > tol
        if not np.any(bad):
            converged = True
            break
        if len(nus) + int(bad.sum()) > maxpts:
            break
        order = np.argsort(np.concatenate((nus, mid[bad])))
        nus = np.concatenate((nus, mid[bad]))[order]
        cs = np.concatenate((cs, cm[bad]))[order]
    return PowerPieces(nus, cs), converged


class _Quad:
    """QUADPACK with its warnings silenced and its error estimates kept: a
    caller sums them per result and `note`s the ratio to that result.
    `last_flag`: whether QUADPACK reported the last call as FAILED (a
    message, `ier != 0`) -- its bound is then worth trusting; on a normal
    exit it is an upper bound far above the error (see
    `ColouredLineshape.__call__`)."""

    def __init__(self):
        self.worst = 0.0
        self.last_flag = False

    def __call__(self, f, a, b, **kw):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', integrate.IntegrationWarning)
            r = integrate.quad(f, a, b, full_output=1, **kw)
        self.last_flag = len(r) > 3
        return r[0], abs(r[1])

    def note(self, err, value):
        self.worst = max(self.worst, err / max(abs(value), 1e-300))


def _aaa(z, F, tol, mmax, W=None):
    """Set-valued AAA (Nakatsukasa, Sete & Trefethen 2018): the barycentric
    rational through support points picked greedily at the largest residual,
    SHARED by every column of `F` (M, k); the residual is weighted by `W`
    (M, k) when given.  Returns the history of `(zj, fj, wj)` fits, the last
    the one that met `tol` on the samples (or used `mmax` support points)."""
    M, k = F.shape
    W = np.ones_like(F) if W is None else W
    rest = np.ones(M, bool)
    R = np.tile(F.mean(axis=0), (M, 1))
    J, hist = [], []
    for m in range(min(mmax, M - 1)):
        err = np.max(np.abs(F - R) * W, axis=1)
        err[~rest] = 0.0
        if m > 0 and err.max() <= tol:
            break
        j = int(np.argmax(err))
        J.append(j)
        rest[j] = False
        C = 1.0 / (z[rest][:, None] - z[J][None, :])
        A = np.vstack([(F[rest, i][:, None] - F[J, i][None, :]) * C
                       for i in range(k)])
        _u, _s, Vh = np.linalg.svd(A, full_matrices=False)
        w = Vh[-1].conj()
        R = F.copy()
        R[rest] = (C @ (w[:, None] * F[J])) / (C @ w)[:, None]
        hist.append((z[J].copy(), F[J].copy(), w.copy()))
    return hist


def _barycentric(z, zj, fj, wj):
    """The rational `(zj, fj, wj)` at `z`, `(len(z), k)`; exact at a support
    point."""
    z = np.asarray(z, dtype=float).ravel()
    d = z[:, None] - zj[None, :]
    hit = d == 0.0
    d[hit] = 1.0
    C = 1.0 / d
    out = (C @ (wj[:, None] * fj)) / (C @ wj)[:, None]
    for r, c in zip(*np.nonzero(hit)):
        out[r] = fj[c]
    return out


def _barycentric_poles(zj, wj):
    """The poles of the barycentric rational: the finite eigenvalues of the
    arrowhead pencil."""
    from scipy.linalg import eig
    m = len(zj)
    B = np.eye(m + 1)
    B[0, 0] = 0.0
    E = np.zeros((m + 1, m + 1), dtype=complex)
    E[0, 1:] = wj
    E[1:, 0] = 1.0
    E[1:, 1:] = np.diag(zj)
    ev = eig(E, B, right=False)
    return ev[np.isfinite(ev)]


class RationalRho:
    """`rho(nu)` on `[nu_a, nu_b]` as ONE barycentric rational in `nu` (set-
    valued AAA, the components sharing support), from ADAPTIVELY chosen
    solves -- the all-orders lineshape's correction in ~25 bordered solves
    where `LogChebyshev` takes ~65.  `fun` takes one `nu` and returns a
    number or a 1-D array; the call returns what `LogChebyshev`'s does.

    Why rational: a white source's `c_fa(nu)` IS a rational function of the
    frequency (its poles the Floquet exponents), and the coloured parts are
    close to one.  Measured offline on four fixtures and a closed form
    (175 captured solves each): 22 .. 26 solves, true error 1e-13 .. 5.9e-7.

    ⚠ THE TOLERANCE IS RELATIVE TO `1 + rho` (floored at `REL_FLOOR`), per
    component: what the lineshape needs is `c_fa = c_dc (1 + rho)`, and
    behind a slow node `1 + rho` falls to ~1e-3, where an absolute 1e-6 on
    `rho` would be 1e-3 of `c_fa`.

    The samples lie on a grid of `RATIONAL_PER_DECADE` per decade down from
    `nu_b`, its decades formed by repeated division as the lineshape's probes
    are, so a caller's cache of those serves them.  Starting from the
    decades, each round fits (tolerance `tol / 10` on the samples) and
    solves where the fit and its one-support-poorer predecessor disagree
    most.  It stops after three consecutive such solves are predicted to
    `tol`.  ⚠⚠ THAT STOP CAN BE FOOLED: a fit in `nu^2` stopped on the slow-
    node 1/f fixture with a true error of 8e-4, and "the fit no longer moves"
    did not catch it either.  So three GUARDS:
      * VERIFY: every half-decade point not yet sampled is solved and must
        be predicted to `tol`; a miss goes back into the samples and the
        search resumes (this caught the 8e-4);
      * POLES: a pole of the fit on the real band (the spurious doublets
        rational fits are known for; one sat at 0.23 f0 in a 5-point fit on
        the coloured slow-node LC, where `rho` has none) has the grid points
        around it solved and the fit redone, up to `POLE_ROUNDS` times;
        then the fit is refused;
      * BUDGET: past `budget` samples it gives up.
    A refusal leaves `converged = False` and `reason`, and the caller falls
    back.  `err` is the worst verified miss, relative; `calls` the solves.

    ⚠ NOT SMOOTH IN ITS INPUT below its accuracy.  The samples are chosen
    greedily, so inputs 1e-13 apart can take different sample paths.  Two
    physically identical noise setups gave lineshapes 7.8e-9 apart, each
    correct to its fit; `LogChebyshev`'s fixed nodes kept them 1.9e-13
    apart.  For finite-difference sensitivities through the lineshape, set
    `PAC.FA_RHO_FIT = 'chebyshev'`.

    History: `doc/shooting_history.md`, `_lineshape.RationalRho`."""

    #: the relative tolerance's floor on `1 + rho`
    REL_FLOOR = 1e-4
    #: rounds of sampling around a pole on the band before the fit is refused
    POLE_ROUNDS = 3
    #: the two guards, switchable only to prove in a test that they bind
    VERIFY = True
    POLE_GUARD = True

    def __init__(self, fun, nu_a, nu_b, tol=CHEB_TOL, budget=RATIONAL_BUDGET,
                 per_decade=RATIONAL_PER_DECADE):
        pd = int(per_decade)
        tops, v = [], float(nu_b)
        while v > float(nu_a) * (1.0 + 1e-12):
            tops.append(v)
            v /= 10.0
        nus = [t * 10.0 ** (-i / pd) for t in tops for i in range(pd)]
        nus = [x for x in nus if x > float(nu_a) * (1.0 + 1e-12)] + [float(nu_a)]
        nus = np.asarray(nus)[::-1]                        # ascending
        self.nu_a, self.nu_b = float(nu_a), float(nu_b)
        self.scale = float(nu_b)
        z = nus / self.scale
        M = len(z)
        vals = {}

        def at(i):
            if i not in vals:
                vals[i] = np.atleast_1d(np.asarray(fun(float(nus[i])), dtype=float))
            return vals[i]

        def weight(F):
            ## 1 / max(|1 + rho|, floor): the error allowed scales with c_fa
            return 1.0 / np.maximum(np.abs(1.0 + np.asarray(F)), self.REL_FLOOR)

        def miss(r_i, i):
            return float(np.max(np.abs(r_i - at(i)) * weight(at(i))))
        ## the decades (the caller's probes) and both ends
        dec = [M - 1 - j * pd for j in range(len(tops)) if M - 1 - j * pd >= 0]
        S = sorted(set(dec) | {0, M - 1})
        half = [i for i in range(M) if (M - 1 - i) % pd == pd // 2]
        self.converged, self.reason, self.err = False, None, np.inf

        def fit(S):
            F = np.asarray([at(i) for i in S])
            return _aaa(z[S], F, 0.1 * tol, len(S) - 1, W=weight(F))

        oos = []                   # out-of-sample misses of the last search

        def search(S):
            conf = 0
            del oos[:]
            while len(S) <= budget:
                hist = fit(S)
                r = _barycentric(z, *hist[-1])
                if len(hist) > 1:
                    d = np.max(np.abs(r - _barycentric(z, *hist[-2])), axis=1)
                else:
                    d = np.ones(M)
                d[S] = -1.0
                if d.max() < 0.0:
                    return S, True
                j = int(np.argmax(d))
                oos.append(miss(r[j], j))
                ok = oos[-1] <= tol
                S = sorted(set(S) | {j})
                conf = conf + 1 if ok else 0
                if conf >= 3:
                    return S, True
            return S, False
        pole_rounds = 0
        while True:
            S, ok = search(S)
            if not ok:
                self.reason = 'budget'
                break
            zj, fj, wj = fit(S)[-1]
            r = _barycentric(z, zj, fj, wj)
            todo = [i for i in half if i not in S] if self.VERIFY else []
            misses = [miss(r[i], i) for i in todo]
            self.err = max(misses + oos[-3:], default=0.0)
            S = sorted(set(S) | set(todo))
            if self.err > tol:
                if len(S) > budget:
                    self.reason = 'budget'
                    break
                continue
            zj, fj, wj = fit(S)[-1]
            lo, hi = z.min(), z.max()
            bad = [p for p in _barycentric_poles(zj, wj)
                   if lo <= p.real <= hi and abs(p.imag) <= 1e-3 * abs(p.real)
                   ] if self.POLE_GUARD else []
            if not bad:
                self.converged = True
                self.zj, self.fj, self.wj = zj, fj, wj
                break
            pole_rounds += 1
            if pole_rounds > self.POLE_ROUNDS:
                self.reason = 'pole'
                break
            ## sample the grid points around each pole, then search again
            for p in bad:
                k = int(np.searchsorted(z, p.real))
                S = sorted(set(S) | {max(k - 1, 0), min(k, M - 1)})
            if len(S) > budget:
                self.reason = 'budget'
                break
        self.calls = len(vals)
        self.vector = next(iter(vals.values())).shape[0] > 1

    def __call__(self, nu):
        nu = np.asarray(nu, dtype=float)
        out = _barycentric(np.clip(nu.ravel(), self.nu_a, self.nu_b)
                           / self.scale, self.zj, self.fj, self.wj)
        if self.vector:
            return out.T.reshape((out.shape[1],) + nu.shape)
        return out[:, 0].reshape(nu.shape)


class LogChebyshev:
    """`rho(nu)` on `[nu_a, nu_b]` as a Chebyshev series in `ln nu`: the
    degree doubles from 16 on nested points (every value reused) until the
    last three coefficients fall below `tol`.  `fun` takes one `nu` and
    returns a number or a 1-D array (several functions on the same points,
    each held to `tol`; the call then returns `(k, len(nu))`).  `err` is
    the tail, `converged` whether it met `tol`, `calls` the number of `fun`
    evaluations."""

    def __init__(self, fun, nu_a, nu_b, tol=CHEB_TOL, nmax=CHEB_MAX, n0=16):
        self.lo, self.hi = math.log(nu_a), math.log(nu_b)
        vals = {}
        top = 1 << int(math.ceil(math.log2(max(int(nmax), int(n0)))))
        N = int(n0)
        while True:
            k = np.arange(N + 1)
            u = 0.5 * (self.lo + self.hi) + 0.5 * (self.hi - self.lo) * np.cos(np.pi * k / N)
            ## the nested points: index j at degree N is index 2j at 2N
            rows = []
            for j in range(N + 1):
                key = j * (top // N)
                if key not in vals:
                    ## the ends at the band edges themselves, not
                    ## exp(log(.)) of them: a caller's cache of those
                    ## (the lineshape's probes) then serves them
                    nu = (nu_b if j == 0 else nu_a if j == N
                          else math.exp(float(u[j])))
                    vals[key] = np.atleast_1d(np.asarray(fun(nu), dtype=float))
                rows.append(vals[key])
            f = np.asarray(rows)                                   # (N+1, k)
            ## the coefficients by the DCT-I
            w = np.ones(N + 1)
            w[0] = w[-1] = 0.5
            c = (2.0 / N) * (np.cos(np.pi * np.outer(k, k) / N) @ (w[:, None] * f))
            c[0] *= 0.5
            c[-1] *= 0.5
            tail = float(np.max(np.abs(c[-3:])))
            if tail <= tol or 2 * N > nmax:
                break
            N *= 2
        self.vector = f.shape[1] > 1
        self.c = c if self.vector else c[:, 0]
        self.err, self.converged = tail, tail <= tol
        self.calls = len(vals)

    def __call__(self, nu):
        x = (2.0 * np.log(np.asarray(nu, dtype=float)) - self.lo - self.hi) / (self.hi - self.lo)
        return np.polynomial.chebyshev.chebval(np.clip(x, -1.0, 1.0), self.c)


class SignedTable:
    """A smooth SIGNED `delta(nu)` on `[nu_a, nu_b]`, tabulated on a uniform
    `ln nu` grid (`TABLE_PER_DECADE`) and read back by four-point Lagrange
    interpolation in pure Python -- cheap enough to sit inside QUADPACK's
    integrand.  Zero outside the band."""

    def __init__(self, fun, nu_a, nu_b, per_decade=TABLE_PER_DECADE):
        self.nu_a, self.nu_b = float(nu_a), float(nu_b)
        n = max(int(math.ceil(per_decade * math.log10(nu_b / nu_a))), 4)
        self.u0 = math.log(nu_a)
        self.du = (math.log(nu_b) - self.u0) / n
        u = self.u0 + self.du * np.arange(-1, n + 3)
        self.h = [float(v) for v in np.asarray(fun(np.exp(u)), dtype=float)]
        self.n = n

    def __call__(self, v):
        if v < self.nu_a or v > self.nu_b:
            return 0.0
        t = (math.log(v) - self.u0) / self.du
        k = min(int(t), self.n - 1)
        x = t - k
        h = self.h
        ## h[k + 1] is the node at u0 + k du (one guard point below)
        xm, x1, x2 = x + 1.0, x - 1.0, x - 2.0
        return (-x * x1 * x2 * h[k] / 6.0 + xm * x1 * x2 * h[k + 1] / 2.0
                - xm * x * x2 * h[k + 2] / 2.0 + xm * x * x1 * h[k + 3] / 6.0)


class ConstantTail:
    """A constant change `level` of `c` on `[nu_b, inf)`: the frequency-aware
    white part held past the band at its level at `fmax`.  Without it the
    white returned UNFILTERED above `fmax`, a jump of ~`c_w` whose ringing
    `sin(2 pi fmax tau)/tau` was as large as `D` itself behind a slow node
    (non-monotone `D`, the tau densities 4.7e-3 apart at 100 linewidths).
    Its structure function is exact:
    `int_{nu_b}^inf (1 - cos beta nu)/nu^2 = (1 - cos beta nu_b)/nu_b +
    beta (pi/2 - Si(beta nu_b))`, `beta pi/2` for small lags and `1/nu_b` at
    infinity."""

    def __init__(self, level, nu_b):
        self.level, self.nu_b = float(level), float(nu_b)

    def structure(self, tau, pref):
        from scipy.special import sici
        beta = 2.0 * np.pi * float(tau)
        x = beta * self.nu_b
        ## neither term cancels: `2 sin^2` and `pi/2 - Si` are both formed
        ## directly
        val = (2.0 * math.sin(0.5 * x) ** 2 / self.nu_b
               + beta * (0.5 * np.pi - float(sici(x)[0])))
        return pref * self.level * val

    def inf(self, pref):
        return pref * self.level / self.nu_b


class ClampedTable(SignedTable):
    """A smooth `c(nu)` tabulated as `SignedTable` is, but HELD at its end
    values outside `[nu_a, nu_b]` (the white level beyond the band) -- the
    second-order skirt's integrand, evaluated tens of thousands of times per
    offset."""

    def __call__(self, v):
        return SignedTable.__call__(self, min(max(float(v), self.nu_a), self.nu_b))


def correction_structure(tau, tab, pref, quad):
    """`pref int delta(nu) (1 - cos 2 pi nu tau) / nu^2 dnu` over the table's
    band, SIGNED: `2 sin^2(beta nu / 2)` below `beta nu = 2` (no
    cancellation), the plain moment minus QAWO above, per decade.  ⚠ `pref`
    is INSIDE the integrand, so QUADPACK's absolute tolerance is one on `D`
    itself -- the quantity `exp(-D/2)` needs absolutely -- whatever the
    frequency scale (`pref = 4 i^2 f0^2` is 4e18 at 1 GHz).  Measured: a
    purely relative 1e-10 changed nothing against the closed form and a
    160-per-decade reference, at 15x the time."""
    beta = 2.0 * np.pi * tau
    na, nb = tab.nu_a, tab.nu_b
    nus = min(max(2.0 / beta, na), nb)
    tot, err = 0.0, 0.0
    if nus > na:
        nd = max(int(math.ceil(math.log10(nus / na))), 1)
        edges = np.geomspace(na, nus, nd + 1)
        f = lambda v: pref * tab(v) * 2.0 * math.sin(0.5 * beta * v) ** 2 / (v * v)
        for lo, hi in zip(edges[:-1], edges[1:]):
            val, e = quad(f, lo, hi, limit=200)
            tot += val
            err += e
    if nus < nb:
        nd = max(int(math.ceil(math.log10(nb / nus))), 1)
        edges = np.geomspace(nus, nb, nd + 1)
        f = lambda v: pref * tab(v) / (v * v)
        for lo, hi in zip(edges[:-1], edges[1:]):
            m, e1 = quad(f, lo, hi, limit=200)
            c, e2 = quad(f, lo, hi, weight='cos', wvar=beta, limit=200)
            tot += m - c
            err += e1 + e2
    return tot


#: Gauss-Legendre points per panel of `edge_increment`'s integral
INCREMENT_GL = 12
#: log panels per decade below its first linear panel
INCREMENT_LOG_PER_DECADE = 8


def increment_nodes(fmin, fmax, T, kmax, lines=()):
    """Nodes and weights on `[fmin, fmax]` for `edge_increment`: panels of
    width `1/(2 kmax T)` from `1/(kmax T)` up -- the kernel's features are
    `1/(k T)` wide and centred on the harmonics, which fall on panel edges --
    and log panels below, down to `fmin` (a 1/f source's weight sits
    there); `INCREMENT_GL` Gauss-Legendre points each.

    `lines`: `(centre, half-width)` pairs in the same frequency -- an
    integrand's own narrow features (`PAC._edge_coloured_law` passes the
    orbital lines, folded).  Each gets panel edges graded geometrically from
    its centre, a quarter half-width wide there and doubling to the base
    panel.  ⚠ A line narrower than a panel falls between its points: a
    Q = 1000 resonator's line (half-width 1.2e-3 f0) read the edge jitter
    29 % low, and 4x the points per panel still 5.7e-4 off.  With no
    `lines` the nodes are the same as before."""
    x, w = np.polynomial.legendre.leggauss(INCREMENT_GL)
    d = 1.0 / (2.0 * float(kmax) * float(T))
    fmin, fmax = float(fmin), float(fmax)
    f_lin = min(max(fmin, 2.0 * d), fmax)
    edges = []
    if f_lin > fmin:
        nd = max(1, int(np.ceil(np.log10(f_lin / fmin)
                                * INCREMENT_LOG_PER_DECADE)))
        edges.extend(np.geomspace(fmin, f_lin, nd + 1)[:-1])
    lin = list(np.arange(f_lin, fmax, d))
    if not lin or fmax - lin[-1] > 1e-9 * d:
        lin.append(fmax)
    else:
        lin[-1] = fmax
    edges.extend(lin)
    edges = np.asarray(edges, dtype=float)
    extra = []
    for c, hw in lines:
        if not float(hw) > 0.0:
            continue
        step = 0.25 * float(hw)
        extra.append(float(c))
        while step < d:
            extra.extend((float(c) - step, float(c) + step))
            step *= 2.0
    if extra:
        extra = np.asarray(extra)
        edges = np.unique(np.concatenate(
            (edges, extra[(extra > fmin) & (extra < fmax)])))
    a, b = edges[:-1], edges[1:]
    nodes = (0.5 * (b - a))[:, None] * x[None, :] + (0.5 * (a + b))[:, None]
    weights = (0.5 * (b - a))[:, None] * w[None, :]
    return nodes.ravel(), weights.ravel()


def edge_increment(coef, ls, s, T, ks, fmin, fmax):
    """`Var_j(k)`, k in `ks`: the variance of the phase increment over `k`
    periods from an instant `t_j` -- `int_{t_j}^{t_j + kT} v^T G xi dt` for
    unit processes `xi` of one-sided power `s` through columns whose PPV
    projection is `v^T G = sum_l V_l e^{j l w0 t}`:

        Var_j(k) = sum_r int_{fmin <= |f| <= fmax} (s(2 pi |f|) / 2)
                   | kT sum_l (-1)^{k l} coef[l, r] sinc(kT (f + l / T)) |^2 df

    with `coef[l] = V_l e^{j l w0 t_j}` and `f` the SOURCE frequency.  The
    phase integral's poles at the harmonics are removable and gone here:
    `sin(pi f kT) / (pi (f + l/T)) = (-1)^{kl} kT sinc(kT (f + l/T))`.

    White (`s` flat) over every frequency it is `c k T` at ANY `t_j`
    (Parseval over whole periods), and its average over `t_j` is the
    stationary structure function `D(kT) / (2 pi f0)^2`; a coloured source
    has memory, and the increment depends on where in the period `t_j`
    sits.  `coef`: an `(nl, r)` array, or a callable of the frequencies
    returning `(nf, nl, r)` (a per-band root, whose power is in its columns:
    then `s` is None).  Both signs of `f` are integrated (complex columns
    need not pair `V_{-l}` with `conj(V_l)`)."""
    ks = np.atleast_1d(np.asarray(ks, dtype=int))
    ls = np.asarray(ls, dtype=float)
    nodes, weights = increment_nodes(fmin, fmax, T, int(ks.max()))
    out = np.zeros(len(ks))
    for sign in (1.0, -1.0):
        f = sign * nodes
        a = coef(f) if callable(coef) else np.asarray(coef)
        sw = (np.full(f.shape, 0.5) if s is None
              else 0.5 * np.asarray(s(2.0 * np.pi * np.abs(f)), dtype=float))
        for i, k in enumerate(ks):
            K = (k * T) * np.sinc((k * T) * (f[:, None] + ls[None, :] / T))
            K = K * ((-1.0) ** ((k * ls).astype(int) % 2))[None, :]
            amp = (K @ a) if a.ndim == 2 else np.einsum('fl,flr->fr', K, a)
            out[i] += float(np.sum(weights * sw
                                   * np.sum(np.abs(amp) ** 2, axis=1)))
    return out


def phase_structure(tau, pc, pref, quad):
    """`D_c(tau) = pref int c(nu) (1 - cos 2 pi nu tau) / nu^2 dnu` over the
    pieces' span."""
    beta = 2.0 * np.pi * tau
    nu0, nuN = pc.nu[0], pc.nu[-1]
    nus = min(max(2.0 / beta, nu0), nuN)
    low = pc.cosine_series(beta, nu0, nus) if nus > nu0 else 0.0
    high, err = 0.0, 0.0
    if nus < nuN:
        high = pc.moment(-2.0, nus, nuN)
        nd = max(int(np.ceil(np.log10(nuN / nus))), 1)
        edges = np.geomspace(nus, nuN, nd + 1)
        f = lambda v: pc(v) / (v * v)
        for lo, hi in zip(edges[:-1], edges[1:]):
            val, e = quad(f, lo, hi, weight='cos', wvar=beta, limit=200)
            high -= val
            err += e
    quad.note(err, low + high)
    return pref * (low + high)


def band_envelope(beta, fun, nu_a, nu_b, nu0, scale, quad):
    """``scale int_{nu_a}^{nu_b} fun(nu)/nu^2 e^{i beta (nu - nu0)} dnu`` --
    a band's oscillatory integral DEMODULATED at `nu0`, by QAWO per decade
    (its cosine and sine parts).  Smooth in `beta` where the band starts at
    `nu0` (`ColouredLineshape`'s edge envelope).  `scale` is INSIDE the
    integrand, so QUADPACK's absolute tolerance is one on `D` (as
    `correction_structure`'s)."""
    nd = max(math.ceil(math.log10(nu_b / nu_a)), 1)
    edges = np.geomspace(nu_a, nu_b, nd + 1)
    f = lambda v: scale * float(fun(v)) / (v * v)
    cs, sn = 0.0, 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        cs += quad(f, lo, hi, weight='cos', wvar=beta, limit=200)[0]
        sn += quad(f, lo, hi, weight='sin', wvar=beta, limit=200)[0]
    return complex(cs, sn) * complex(math.cos(beta * nu0), -math.sin(beta * nu0))


def _i0m1(r):
    """`I_0(r) - 1`, without the cancellation at small `r`."""
    if r < 0.1:
        q = 0.25 * r * r
        term, tot = 1.0, 0.0
        for m in range(1, 8):
            term *= q / (m * m)
            tot += term
        return tot
    return float(special.iv(0, r)) - 1.0


class _ScaledTable:
    """A `SignedTable` read in the lineshape's own units (`ColouredLineshape`:
    frequencies over `F`, the density times `F`)."""

    def __init__(self, tab, F):
        self.tab, self.F = tab, float(F)
        self.nu_a, self.nu_b = tab.nu_a / F, tab.nu_b / F

    def __call__(self, v):
        return self.F * self.tab(v * self.F)


class ColouredLineshape:
    """`S(f)` for a white rate `a` (the Lorentzian's `2 pi^2 i^2 f0^2 c_w`)
    and a coloured `c_c` (`PowerPieces`, or None), `D_c = pref int ...`,
    `pref = 4 i^2 f0^2` -- plus, optionally, a SIGNED correction `corr` to
    `c`: a `SignedTable`, or a list of them on their own bands and
    `ConstantTail`s, summed.  Its `D` (`correction_structure`, the tails
    exact) enters the spline of the TOTAL `D`.

    ⚠ IN UNITS OF `F = i f0` INSIDE (frequencies over `F`, lags times `F`,
    `c` times `F`, `pref` 4): QUADPACK's tolerances are absolute, and on
    the physical scale they did not bind -- the same physics at `f0 = 1e9`
    read the lineshape 1.9 % off at 100 times the band's lower edge, 3.4e-3
    at 5 times, its own QUADPACK estimate 16 % (until 2026-10-01; every
    reference test runs at `f0 = 1`, where nothing moves).

    ⚠ THE BAND'S LOWER EDGE RINGS.  The colour stops at `nu0 = fmin`, so
    `D_c` oscillates at `nu0` in the lag, ``D_c = D_inf - Re[e^{i 2 pi nu0
    tau} gamma(tau)]`` with the envelope ``gamma = pref int h e^{i beta (nu -
    nu0)}``, ``h = c/nu^2``, smooth.  A spline of `log D` in `log tau`
    ALIASES that ringing: around the edge the lineshape read 29 % high at
    0.55 `fmin` on an LC oscillator (its two grid phases 46 % apart), and
    moved with `offset_fmax` (the grid's phase) -- the "1 % far skirt" of
    the review (until 2026-10-01).  So past ``beta nu0 = EDGE_SWITCH`` the
    ringing is carried EXACTLY: `gamma` splined (smooth), its asymptotic
    series past the grid; the transform there on half-cycle panels until the
    ringing's amplitude falls below `EDGE_HARMONIC_R`, and beyond by the
    harmonics of `nu0` (``e^{r cos t} = I_0(r) + 2 sum I_k(r) cos k t``), each
    by QUADPACK's Fourier rule.  A correction table that starts at the same
    edge joins the envelope; the rest of `D` (the white part, the other
    tables, the tails) is splined beside it, as before.  Against an mpmath
    reference around the edge (`benchmarks/lineshape_reference.py`, regime
    'edge', 0.3 .. 20 `fmin`): <= 8.2e-9 to 5 `fmin`, 1.6e-7 at 20, the
    two grid phases within ~10x (were 1.3e-3 .. 4.3e-3).  Where the
    integrand is dead before the switch -- the white part's `e^{-a tau}`, or
    a strong colour's `D` itself -- nothing of this runs."""

    def __init__(self, a, pc, pref, per_decade=TAU_PER_DECADE, corr=None,
                 shift=False):
        ## the lineshape's own units (see the class note); `F = 1` at `pref
        ## = 4` changes nothing
        F = 0.5 * math.sqrt(float(pref))
        self.F = F
        corr = ([] if corr is None else
                list(corr) if isinstance(corr, (list, tuple)) else [corr])
        if F != 1.0:
            if pc is not None:
                pc = PowerPieces(pc.nu / F, pc.c * F)
            corr = [ConstantTail(t.level * F, t.nu_b / F)
                    if isinstance(t, ConstantTail) else _ScaledTable(t, F)
                    for t in corr]
            pref = 4.0
        self.a = float(a) / F
        self.pc, self.pref = pc, float(pref)
        a, pref = self.a, self.pref
        self.corr = corr or None
        self.quad = _Quad()
        self.last_err = 0.0
        self.tails = [t for t in corr if isinstance(t, ConstantTail)]
        tables = [t for t in corr if not isinstance(t, ConstantTail)]
        bands = ([(pc.nu[0], pc.nu[-1])] if pc is not None else []) + \
            [(t.nu_a, t.nu_b) for t in tables]
        nu0, nuN = min(b[0] for b in bands), max(b[1] for b in bands)
        self.tlo, self.thi = 1e-3 / nuN, 1e3 / nu0
        n = int(np.ceil(per_decade * np.log10(self.thi / self.tlo)))
        self.taus = np.geomspace(self.tlo, self.thi, n + 1)
        if shift:
            ## the interior nodes half a step along: the same density at the
            ## other phase (the error estimate's second grid)
            lt = np.log(self.taus)
            self.taus = np.exp(np.concatenate(
                ([lt[0]], 0.5 * (lt[:-1] + lt[1:]), [lt[-1]])))
        ## the colour's edge: below `ta` the splines, past it the ringing
        ## carried exactly (`_edge_init`); the splines on the nodes to `2 ta`.
        ## (Not where the white part has killed the integrand before `ta`:
        ## there the splines run to `thi`, as before.)
        ## (Nor where `D` itself has: a strong colour, `D(ta)/2` past
        ## `EDGE_DECAY` -- `D` rings about `D_inf` from there, so the whole
        ## integrand beyond is `e^{-EDGE_DECAY}` small.)
        self.ta = np.inf
        self.edge = False
        if pc is not None:
            ta = EDGE_SWITCH / (2.0 * np.pi * pc.nu[0])
            dead = corr == [] and self.a > 0.0 and EDGE_DECAY / self.a <= ta
            if not dead and corr == []:
                Dta = phase_structure(ta, pc, pref, _Quad())
                dead = min(Dta, pref * pc.moment(-2.0, pc.nu[0], pc.nu[-1])) \
                    > 4.0 * EDGE_DECAY
            self.edge = not dead
            self.ta = ta if self.edge else np.inf
        low = self.taus[self.taus <= 2.0 * self.ta]
        self.Dinf = 0.0
        if pc is not None:
            p0, pN = pc.nu[0], pc.nu[-1]
            self.Dinf = pref * pc.moment(-2.0, p0, pN)
            ## D_c ~ q2 tau^2 below the grid; the band edge rings above it
            self.q2 = pref * 2.0 * np.pi ** 2 * pc.moment(0.0, p0, pN)
            self.h0 = float(pc(p0)) / p0 ** 2
            Ds = np.array([phase_structure(t, pc, pref, self.quad)
                           for t in low])
            self.spline = interpolate.CubicSpline(
                np.log(low), np.log(np.maximum(Ds, 1e-300)))
        self.Dpow_inf = self.Dinf
        self._tab_inf = {}
        if corr:
            ## its moments by QUADPACK on each table, per decade
            self.cinf = sum(t.inf(pref) for t in self.tails)
            self.cq2 = 0.0
            for tab in tables:
                ed = np.geomspace(tab.nu_a, tab.nu_b,
                                  max(int(np.ceil(np.log10(tab.nu_b / tab.nu_a))), 1) + 1)
                f2 = lambda v, tab=tab: pref * tab(v) / (v * v)
                tinf = sum(self.quad(f2, lo, hi, limit=200)[0]
                           for lo, hi in zip(ed[:-1], ed[1:]))
                self._tab_inf[id(tab)] = tinf
                self.cinf += tinf
                f0_ = lambda v, tab=tab: pref * tab(v)
                self.cq2 += 2.0 * np.pi ** 2 * sum(
                    self.quad(f0_, lo, hi, limit=200)[0]
                    for lo, hi in zip(ed[:-1], ed[1:]))
            ## ⚠ THE TOTAL `D` IS SPLINED, IN LOGS.  A correction that takes
            ## the white part away up to `fmax` (a slow node) is `~ -2 a tau`
            ## there, and the total is ~1e-3 of either term: a spline of the
            ## correction alone, good to 3e-8 of ITSELF, left the total
            ## ~3e-5 off.  The node values are sums of accurately computed
            ## parts, and the spline's error then scales with `D`.
            Dt = np.array([2.0 * self.a * t
                           + (float(np.exp(self.spline(np.log(t))))
                              if pc is not None else 0.0)
                           + sum(correction_structure(t, tab, pref, self.quad)
                                 for tab in tables)
                           + sum(tl.structure(t, pref) for tl in self.tails)
                           for t in low])
            self.tspline = interpolate.CubicSpline(
                np.log(low), np.log(np.maximum(Dt, 1e-300)))
            ## `Dinf` is the non-white TOTAL at infinite lag, which `g` and
            ## the head subtract
            self.Dinf = self.Dinf + self.cinf
        if self.edge:
            self._edge_init(tables)

    def _edge_init(self, tables):
        """The edge's envelope `gamma` (and, with a correction, the rest of
        `D`, `Q`) past `ta`, splined; where the transform hands over to the
        harmonics (`tb`)."""
        pc, pref, a = self.pc, self.pref, self.a
        nu0 = float(pc.nu[0])
        self.nu_e, self.Om = nu0, 2.0 * np.pi * nu0
        self._edge_tabs = [t for t in tables
                           if abs(t.nu_a - nu0) <= 1e-12 * nu0]
        rest = [t for t in tables if not any(t is e for e in self._edge_tabs)]
        ## (no correction: the integrand decays as `e^{-a tau}`; past
        ## `EDGE_DECAY / a` it is below `e^{-EDGE_DECAY}`, and nothing is
        ## computed there -- but the envelope's grid reaches `30 / nu0`,
        ## where its endpoint series holds to ~1e-5 of the ringing)
        self.tend = self.thi
        if self.corr is None and a > 0.0:
            self.tend = min(self.thi, max(EDGE_DECAY / a, 30.0 / nu0))
        ## (strictly inside, by more than rounding: the grid's own end is
        ## `thi` to within an ulp)
        tg = self.taus[(self.taus > self.ta * (1.0 + 1e-9))
                       & (self.taus < self.tend * (1.0 - 1e-9))]
        tg = np.concatenate(([self.ta], tg, [self.tend]))
        if tg.size < 4:
            tg = np.geomspace(self.ta, max(self.tend, 2.0 * self.ta), 4)
        gam = []
        for t in tg:
            b = 2.0 * np.pi * t
            g = band_envelope(b, pc, nu0, pc.nu[-1], nu0, pref, self.quad)
            for tab in self._edge_tabs:
                g += band_envelope(b, tab, tab.nu_a, tab.nu_b, nu0, pref,
                                   self.quad)
            gam.append(g)
        gam = np.asarray(gam)
        lb = np.log(2.0 * np.pi * tg)
        self._gre = interpolate.CubicSpline(lb, 2.0 * np.pi * tg * gam.real)
        self._gim = interpolate.CubicSpline(lb, 2.0 * np.pi * tg * gam.imag)
        self._tg_end = float(tg[-1])
        ## past the grid the endpoint series, ``G ~ i h0/b - h1/b^2 - i
        ## h2/b^3`` (`h = c/nu^2` at the edge, the first piece's power law;
        ## a correction table's leading two terms)
        p = float(pc.p[0])
        h0 = float(pc.c[0]) / nu0 ** 2
        hs = [h0, h0 * (p - 2.0) / nu0, h0 * (p - 2.0) * (p - 3.0) / nu0 ** 2]
        for tab in self._edge_tabs:
            e = 1e-4 * nu0
            t0 = float(tab(nu0)) / nu0 ** 2
            t1 = (float(tab(nu0 + e)) / (nu0 + e) ** 2 - t0) / e
            hs[0] += t0
            hs[1] += t1
        self._hs = hs
        ## the rest of the TOTAL `D` past `ta`, `Q = D + Re[e^{i Om t} gamma]`:
        ## with no correction `2 a t + D_inf` exactly; with one, its parts at
        ## the nodes (the edge's at infinite lag), splined in logs
        self._Qspl = None
        if self.corr is not None:
            base = self.Dpow_inf + sum(self._tab_inf[id(t)]
                                       for t in self._edge_tabs)
            Q = np.array([2.0 * a * t + base
                          + sum(correction_structure(t, tab, pref, self.quad)
                                for tab in rest)
                          + sum(tl.structure(t, pref) for tl in self.tails)
                          for t in tg])
            self._Qlog = bool(np.all(Q > 0.0))
            self._Qspl = interpolate.CubicSpline(
                np.log(tg), np.log(Q) if self._Qlog else Q)
        ## the hand-over to the harmonics: the first node past which the
        ## ringing's amplitude `r = |gamma|/2` stays below EDGE_HARMONIC_R
        r = np.abs(gam) / 2.0
        above = np.nonzero(r >= EDGE_HARMONIC_R)[0]
        k = int(above[-1]) + 1 if above.size else 0
        self.tb = float(tg[min(k, tg.size - 1)])
        self.tb = min(max(self.tb, 2.0 * self.ta), self._tg_end)

    def _gamma(self, t):
        """The edge's envelope at lag `t` (own units)."""
        b = 2.0 * np.pi * t
        if t <= self._tg_end:
            lb = math.log(b)
            return complex(float(self._gre(lb)), float(self._gim(lb))) / b
        h0, h1, h2 = self._hs
        return self.pref * complex(-h1 / b ** 2, h0 / b - h2 / b ** 3)

    def _ring(self, t):
        g = self._gamma(t)
        ph = self.Om * t
        return g.real * math.cos(ph) - g.imag * math.sin(ph)

    def _Q(self, t):
        """`D` without the edge's ringing past `ta` (own units)."""
        if self._Qspl is None or t > self._tg_end:
            return 2.0 * self.a * t + self.Dinf
        v = float(self._Qspl(math.log(t)))
        return math.exp(v) if self._Qlog else v

    def _past_decay(self, t):
        ## (no correction, past `EDGE_DECAY / a`: the integrand is zero)
        return self.corr is None and self.a > 0.0 and t > self.tend

    def Dtotal(self, tau):
        """`D` with its white part `2 a tau`: with a correction, from the
        spline of the total (see `__init__`), never a difference."""
        return self._Dtotal(float(tau) * self.F)

    def D(self, tau):
        """`D` without its white part (the power-law pieces', and the
        correction's past the grid; on the grid a correction's `D` is
        `Dtotal - 2 a tau`)."""
        return self._D(float(tau) * self.F)

    def _Dtotal(self, tau):
        if self.edge and tau >= self.ta:
            return self._Q(tau) - self._ring(tau)
        if self.corr is None:
            return 2.0 * self.a * tau + self._D(tau)
        if tau <= 0.0:
            return 0.0
        if tau < self.tlo:
            q2 = self.q2 if self.pc is not None else 0.0
            return (2.0 * self.a * tau + (q2 + self.cq2) * tau * tau
                    + sum(tl.structure(tau, self.pref) for tl in self.tails))
        if tau > self.thi:
            return 2.0 * self.a * tau + self._D(tau)
        return float(np.exp(self.tspline(np.log(tau))))

    def _D(self, tau):
        if tau <= 0.0:
            return 0.0
        if self.edge and tau >= self.ta:
            return self._Dtotal(tau) - 2.0 * self.a * tau
        if self.corr is not None and self.tlo <= tau <= self.thi:
            return self._Dtotal(tau) - 2.0 * self.a * tau
        out = 0.0
        if self.pc is not None:
            if tau < self.tlo:
                out = self.q2 * tau * tau
            elif tau > self.thi:
                beta = 2.0 * np.pi * tau
                out = (self.Dpow_inf + self.pref * self.h0
                       * np.sin(beta * self.pc.nu[0]) / beta)
            else:
                out = float(np.exp(self.spline(np.log(tau))))
        if self.corr is not None:
            if tau < self.tlo:
                out += self.cq2 * tau * tau + sum(
                    tl.structure(tau, self.pref) for tl in self.tails)
            else:
                out += self.cinf
        return out

    def g(self, tau):
        """`e^{-a tau} (e^{-D/2} - e^{-D_inf/2})`, in logs (no 0 * inf);
        with a correction as `e^{-Dtotal/2} - e^{-a tau - D_inf/2}`, so the
        white part is never subtracted from the total.  (Own units.)"""
        if self.edge and self._past_decay(tau):
            return 0.0
        if self.corr is None:
            x, y = -0.5 * self._D(tau), -0.5 * self.Dinf
            m = max(x, y)
            diff = np.exp(m - self.a * tau) * (-np.expm1(min(x, y) - m))
            return diff if x >= y else -diff
        x, y = -0.5 * self._Dtotal(tau), -self.a * tau - 0.5 * self.Dinf
        m = max(x, y)
        diff = np.exp(m) * (-np.expm1(min(x, y) - m))
        return diff if x >= y else -diff

    @property
    def line_weight(self):
        """The carrier line's weight when no white part broadens it."""
        return float(np.exp(-0.5 * self.Dinf))

    def whole(self, tau):
        """`exp(-D/2)` itself (with the white part), at most 1.  (Own units.)"""
        return float(np.exp(-0.5 * self._Dtotal(tau)))

    def __call__(self, f):
        return self._transform(abs(float(f)) / self.F) / self.F

    def _transform(self, f):
        w = 2.0 * np.pi * f
        a = self.a
        ## ⚠ A NEGATIVE `D_inf` (a signed correction that takes white noise
        ## away over a band holding the core): `exp(-D_inf/2) L_w` and the
        ## integral of `g` are then each `exp(|D_inf|/2)` large and cancel
        ## (e^33 against e^33 at a corner 0.03 linewidths out, every digit
        ## lost; 4.8e-7 split, against the closed form).  So the body
        ## integrates `exp(-D/2) <= 1` itself over the grid, and the
        ## Lorentzian enters as its closed-form tail past the grid's end,
        ## `2 Re[e^{(-a + i w) t_end} / (a - i w)]`, weighted `exp(-D_inf/2)`.
        ## A white part is then present (`D` grows as `2 a tau`, so `D_inf`
        ## < 0 needs `a > 0`).  Below `SPLIT_GAIN` the subtracted form is
        ## the more accurate one and stays.
        t_end = self.tb if self.edge else self.thi
        split = a > 0.0 and -0.5 * self.Dinf > math.log(SPLIT_GAIN)
        if split:
            z = complex(-a, w) * t_end
            head = (np.exp(-0.5 * self.Dinf + z.real) * 2.0
                    * (complex(np.cos(z.imag), np.sin(z.imag)) / complex(a, -w)).real)
            body = self.whole
        else:
            lw = (2.0 * a / (a * a + w * w)) if a > 0.0 else 0.0
            head = np.exp(-0.5 * self.Dinf) * lw
            body = self.g
        ## per decade of the whole grid (to `ta` on an edge: the same panels
        ## below it as without one)
        nd = int(np.ceil(np.log10(self.thi / self.tlo)))
        edges = np.concatenate(([0.0], np.geomspace(self.tlo, self.thi, nd + 1)))
        if self.edge:
            edges = np.concatenate((edges[edges < self.ta * (1.0 - 1e-9)],
                                    [self.ta]))
            ## past `ta` the ringing, exactly, half a cycle a panel
            npan = max(int(np.ceil((self.tb - self.ta) * 2.0 * self.nu_e)), 1)
            edges = np.concatenate((edges, np.linspace(self.ta, self.tb,
                                                       npan + 1)[1:]))
        tot, err, err_failed = 0.0, 0.0, 0.0
        for lo, hi in zip(edges[:-1], edges[1:]):
            kw = ({'limit': 400} if w == 0.0 else
                  {'weight': 'cos', 'wvar': w, 'limit': 400})
            val, e = self.quad(body, lo, hi, **kw)
            tot += val
            err += e
            err_failed += e if self.quad.last_flag else 0.0
        if self.edge:
            val, e, ef = self._edge_tail(w)
            tot += val
            err += e
            err_failed += ef
        else:
            kw = ({'limit': 400} if w == 0.0 else
                  {'weight': 'cos', 'wvar': w, 'limlst': 200})
            val, e = self.quad(self.g, self.thi, np.inf, **kw)
            tot += val
            err += e
            err_failed += e if self.quad.last_flag else 0.0
        out = head + 2.0 * tot
        self.quad.note(2.0 * err, out)
        ## ⚠ QUADPACK'S BOUND IS NOT THIS VALUE'S ERROR.  It sums absolute
        ## bounds on O(1) pieces whose tiny difference is the far skirt: 100
        ## .. 10^4 times the true error against an mpmath reference.  The
        ## grid-phase difference (`handover`) is the honest part; QUADPACK's
        ## bound counts only where QUADPACK itself reported failure.
        ## History: `doc/shooting_history.md`, `_lineshape.ColouredLineshape`.
        self.last_bound = 2.0 * err / max(abs(out), 1e-300)
        self.last_err = 2.0 * err_failed / max(abs(out), 1e-300)
        return out

    def _edge_tail(self, w):
        """``int_tb^inf g(tau) cos(w tau)`` by the harmonics of the edge:
        ``g = e^{-Q/2} e^{r cos(Om tau + psi)} - e^{-a tau - D_inf/2}``, ``r =
        |gamma|/2``, ``psi = arg gamma``, and ``e^{r cos t} = I_0(r) + 2 sum_k
        I_k(r) cos k t`` -- a smooth part and, per harmonic, two smooth
        envelopes against ``cos/sin((k Om +- w) tau)``, QUADPACK's Fourier
        rule each.  Returns `(value, bound, bound where QUADPACK failed)`."""
        a, Om = self.a, self.Om
        rb = abs(self._gamma(self.tb)) / 2.0
        K = 1
        i1 = max(float(special.iv(1, rb)), 1e-300)
        while K < EDGE_MAX_HARMONICS and float(special.iv(K + 1, rb)) > 1e-17 * i1:
            K += 1
        self.last_harmonics = K

        def parts(t):
            ## `(e^{-Q/2}, r, psi, the smooth part)`, None past the decay
            if self._past_decay(t):
                return None
            g = self._gamma(t)
            r, psi = abs(g) / 2.0, math.atan2(g.imag, g.real)
            Q = self._Q(t)
            eq = math.exp(-0.5 * Q)
            ref = -a * t - 0.5 * self.Dinf
            ## the smooth part, both differences formed without cancelling
            sm = eq * _i0m1(r) + math.exp(ref) * math.expm1(
                -0.5 * (Q - 2.0 * a * t - self.Dinf))
            return eq, r, psi, sm

        def smooth(t):
            p = parts(t)
            return 0.0 if p is None else p[3]

        def harmonic(k, sine):
            def fn(t):
                p = parts(t)
                if p is None:
                    return 0.0
                eq, r, psi = p[0], p[1], p[2]
                amp = 2.0 * eq * float(special.iv(k, r))
                return -amp * math.sin(k * psi) if sine else amp * math.cos(k * psi)
            return fn
        tot, err, err_failed = 0.0, 0.0, 0.0
        ## per decade to the envelope's grid end, then QUADPACK's Fourier
        ## rule (a term at `W = 0` -- an offset ON a harmonic of `nu0` --
        ## integrates its slowly decaying envelope a decade at a time)
        t1 = max(self._tg_end, self.tb)
        nd = max(math.ceil(math.log10(t1 / self.tb)), 1) if t1 > self.tb else 0
        dec = np.geomspace(self.tb, t1, nd + 1) if nd else np.array([self.tb])
        beyond = not (self.corr is None and a > 0.0)

        def add(fn, W, sine, fac):
            nonlocal tot, err, err_failed
            if W == 0.0 and sine:
                return
            sg = -1.0 if (sine and W < 0.0) else 1.0
            kw = ({} if W == 0.0 else
                  {'weight': 'sin' if sine else 'cos', 'wvar': abs(W)})
            pieces = [(lo, hi, dict(kw, limit=400))
                      for lo, hi in zip(dec[:-1], dec[1:])]
            if beyond:
                pieces.append((t1, np.inf, dict(kw, limlst=200) if kw
                               else {'limit': 400}))
            for lo, hi, kwp in pieces:
                val, e = self.quad(fn, lo, hi, **kwp)
                tot += fac * sg * val
                err += fac * abs(e)
                err_failed += fac * abs(e) if self.quad.last_flag else 0.0
        add(smooth, w, False, 1.0)
        for k in range(1, K + 1):
            for sine in (False, True):
                fn = harmonic(k, sine)
                for W in (k * Om + w, k * Om - w):
                    add(fn, W, sine, 0.5)
        return tot, err, err_failed


#: the second-order skirt splits the phase at `f / SKIRT_SPLIT`, and its
#: error estimate is the move to `f / SKIRT_SPLIT_ALT`
SKIRT_SPLIT = 10.0
SKIRT_SPLIT_ALT = 4.0


def _log_quad(fn, a, b):
    """`int_a^b fn(nu) dnu`, `0 < a < b`, by QUADPACK per decade in `ln nu`,
    held relative (every integrand here is positive)."""
    if not b > a:
        return 0.0
    nd = max(int(math.ceil(math.log10(b / a))), 1)
    edges = np.geomspace(a, b, nd + 1)
    tot = 0.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', integrate.IntegrationWarning)
        for lo, hi in zip(edges[:-1], edges[1:]):
            tot += integrate.quad(lambda u: fn(math.exp(u)) * math.exp(u),
                                  math.log(lo), math.log(hi), limit=200,
                                  epsabs=0.0, epsrel=1e-11)[0]
    return tot


def second_order_skirt(f, cfun, q, nu_low, nu_high, c_at_f=None,
                       split=SKIRT_SPLIT):
    """The lineshape's far skirt to SECOND order in the phase.

    The linear skirt is `P(f) = q c(f) / f^2` (`q = i^2 f0^2`, `c` the whole
    `c(nu)` from `cfun`; `c_at_f` its exact value at `f` where the caller
    has one).  Split the phase at `nu_x = f / split` into a core `L` and the
    rest `H`; then exactly `S = S_L * S_H`, and to second order

        S2 = exp(-sH2) [ P(f) + (M2/2) P''(f) + (1/2) (P_H * P_H)(f) ],
        sH2 = 2 int_{nu_x}^inf P,   M2 = 2 int_0^{nu_x} nu^2 P,

    `sH2` the phase variance above the split, `M2` the core's spread in
    frequency, the convolution over `|g|, |f - g| > nu_x`.  These are the
    two terms `linear_error` bounds, now taken as corrections, so the error
    is third order.  Against an mpmath reference on a real oscillator's
    line (linewidth 5.5e-8 f0, 1/f from 1e-7 f0): 7.8e-4 / 6.5e-5 / 6.9e-6
    first order -> 3.0e-6 / 9.0e-8 / 8.1e-9 at 3e-3 / 1e-2 / 3e-2 f0.  On a
    very broad line (3e-4 f0) the expansion parameter is ~1e-2 and it gains
    only ~10x.  `c` is taken constant below `nu_low` and above ~`10
    nu_high` (the white level).

    History: `doc/shooting_history.md`, `_lineshape.second_order_skirt`."""
    f = abs(float(f))
    x = f / float(split)
    c = lambda v: float(cfun(abs(v)))
    P = lambda v: q * c(v) / (v * v)
    G = max(1e3 * f, 10.0 * float(nu_high))
    cG = c(G)
    lo = min(float(nu_low), 0.5 * x)
    sH2 = 2.0 * q * (_log_quad(lambda v: c(v) / (v * v), x, G) + cG / G)
    M2 = 2.0 * q * (_log_quad(c, lo, x) + c(lo) * lo)
    tail = q * q * cG * cG / (3.0 * G ** 3)
    conv = (_log_quad(lambda g: P(g) * P(f + g), x, G) + tail
            + _log_quad(lambda g: P(g) * P(g - f), f + x, G) + tail)
    if f / 2.0 > x:
        conv += 2.0 * _log_quad(lambda g: P(g) * P(f - g), x, f / 2.0)
    ## P'' by central differences, Richardson over two steps
    d = lambda h: (P(f + h) - 2.0 * P(f) + P(f - h)) / (h * h)
    P2 = (4.0 * d(0.02 * f) - d(0.04 * f)) / 3.0
    Pf = q * (float(c_at_f) if c_at_f is not None else c(f)) / (f * f)
    return float(np.exp(-sH2) * (Pf + 0.5 * M2 * P2 + 0.5 * conv))


def handover(shapes, off, skirt):
    """The lineshape per offset: the TRANSFORM (`shapes[1]`) or the second-
    order SKIRT (`skirt(o)` -> its value at the two splits and the linear
    skirt `P`), whichever carries the smaller ESTIMATED error.  Returns
    `(vals, errs)`.

    The transform's estimate is its move to the other grid phase
    (`shapes[0]`) plus QUADPACK's bound only where QUADPACK reported
    failure.  The skirt's is its move between the two splits PLUS the square
    of its own second-order correction, `((S2 - P)/P)^2`, the size of the
    next term.  ⚠ The split move alone reads ZERO where the expansion has
    collapsed (`exp(-sH2)` underflowing to 0 at both splits).  A
    non-positive skirt is never taken.  Against the mpmath reference both
    estimates track the true error within ~2.4x and ~5x.

    History: `doc/shooting_history.md`, `_lineshape.handover`."""
    vals, errs = [], []
    for o in np.atleast_1d(off).ravel():
        s1, s2 = shapes[0](o), shapes[1](o)
        err = abs(s2 - s1) / max(abs(s2), 1e-300) + shapes[1].last_err
        if o != 0.0:
            k1, k2, p1 = skirt(o)
            if k1 > 0.0 and p1 > 0.0:
                el = (abs(k1 - k2) / k1) + ((k1 - p1) / p1) ** 2
            else:
                el = np.inf
            if el < err:
                s2, err = k1, el
        vals.append(s2)
        errs.append(err)
    return vals, errs


def linear_error(f, i, f0, c_w, pc):
    """How far the lineshape at offset `f` may sit from its LINEAR skirt
    `S_phi(f)`: the core's spread, `6 sigma^2 / f^2` (a skirt of slope 2 .. 3
    convolved with a core of variance `sigma^2`), plus the phase variance
    above `nu_x`, `D_>(nu_x) / 2`, split at `nu_x = f / sqrt(6)` where the two
    balance for a white source.  An ESTIMATE (the second order of the
    expansion), not a bound.  `pc` None: a white line."""
    f = abs(float(f))
    nx = f / np.sqrt(6.0)
    if pc is None:
        ## (a white line: no coloured pieces)
        sig2 = 2.0 * i * i * f0 * f0 * c_w * nx
        above = 2.0 * i * i * f0 * f0 * c_w / nx
        return 6.0 * sig2 / (f * f) + above
    lo = min(max(nx, pc.nu[0]), pc.nu[-1])
    sig2 = 2.0 * i * i * f0 * f0 * (c_w * nx + pc.moment(0.0, pc.nu[0], lo))
    above = 2.0 * i * i * f0 * f0 * (c_w / nx + pc.moment(-2.0, lo, pc.nu[-1]))
    return 6.0 * sig2 / (f * f) + above
