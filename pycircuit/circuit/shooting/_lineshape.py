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
2026-09-26; `PAC._fa_lineshape`) cannot ride on the power-law pieces, which
cannot pass a zero.  It is tabulated (`SignedTable`) and its structure
function taken by QUADPACK (`correction_structure`, `2 sin^2` below
`beta nu = 2`), splined in its VALUE beside `D_c`, and added.  `D_inf` may
then be NEGATIVE (white noise taken away over a band holding the core), and
past `SPLIT_GAIN` the transform integrates `exp(-D/2)` itself, with the
Lorentzian a closed-form tail.  With a correction the TOTAL `D` is
splined (the white part and a correction that removes it cancel), and past
the band the white part is held at its corrected level (`ConstantTail`).
Gated against a closed form: <= 4.8e-7.
`LogChebyshev` represents the correction's `rho = c_fa/c_dc - 1` from its
solves.
"""
import math
import warnings

import numpy as np
from scipy import integrate, interpolate

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
#: a signed correction is tabulated at this many points per decade of `nu`
#: for the quadrature (`SignedTable`)
TABLE_PER_DECADE = 2000
#: past this gain `exp(-D_inf/2)` the transform integrates `exp(-D/2)`
#: itself (`ColouredLineshape.__call__`): the subtracted form amplifies
#: QUADPACK's 1.5e-8 by the gain, and below it is the (slightly) better one
#: -- 2.1e-7 against 3.5e-7 at gain 19 on the closed form
SPLIT_GAIN = 50.0


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
    caller sums them per result and `note`s the ratio to that result."""

    def __init__(self):
        self.worst = 0.0

    def __call__(self, f, a, b, **kw):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', integrate.IntegrationWarning)
            val, err = integrate.quad(f, a, b, **kw)[:2]
        return val, abs(err)

    def note(self, err, value):
        self.worst = max(self.worst, err / max(abs(value), 1e-300))


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
                    vals[key] = np.atleast_1d(np.asarray(
                        fun(math.exp(float(u[j]))), dtype=float))
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


class ColouredLineshape:
    """`S(f)` for a white rate `a` (the Lorentzian's `2 pi^2 i^2 f0^2 c_w`)
    and a coloured `c_c` (`PowerPieces`, or None), `D_c = pref int ...`,
    `pref = 4 i^2 f0^2` -- plus, optionally, a SIGNED correction `corr` to
    `c`: a `SignedTable`, or a list of them on their own bands and
    `ConstantTail`s, summed.  Its `D` (`correction_structure`, the tails
    exact) enters the spline of the TOTAL `D`."""

    def __init__(self, a, pc, pref, per_decade=TAU_PER_DECADE, corr=None):
        self.a = float(a)
        self.pc, self.pref = pc, float(pref)
        corr = ([] if corr is None else
                list(corr) if isinstance(corr, (list, tuple)) else [corr])
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
        self.Dinf = 0.0
        if pc is not None:
            p0, pN = pc.nu[0], pc.nu[-1]
            self.Dinf = pref * pc.moment(-2.0, p0, pN)
            ## D_c ~ q2 tau^2 below the grid; the band edge rings above it
            self.q2 = pref * 2.0 * np.pi ** 2 * pc.moment(0.0, p0, pN)
            self.h0 = float(pc(p0)) / p0 ** 2
            Ds = np.array([phase_structure(t, pc, pref, self.quad)
                           for t in self.taus])
            self.spline = interpolate.CubicSpline(
                np.log(self.taus), np.log(np.maximum(Ds, 1e-300)))
        self.Dpow_inf = self.Dinf
        if corr:
            ## its moments by QUADPACK on each table, per decade
            self.cinf = sum(t.inf(pref) for t in self.tails)
            self.cq2 = 0.0
            for tab in tables:
                ed = np.geomspace(tab.nu_a, tab.nu_b,
                                  max(int(np.ceil(np.log10(tab.nu_b / tab.nu_a))), 1) + 1)
                f2 = lambda v, tab=tab: pref * tab(v) / (v * v)
                self.cinf += sum(self.quad(f2, lo, hi, limit=200)[0]
                                 for lo, hi in zip(ed[:-1], ed[1:]))
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
                           for t in self.taus])
            self.tspline = interpolate.CubicSpline(
                np.log(self.taus), np.log(np.maximum(Dt, 1e-300)))
            ## `Dinf` is the non-white TOTAL at infinite lag, which `g` and
            ## the head subtract
            self.Dinf = self.Dinf + self.cinf

    def Dtotal(self, tau):
        """`D` with its white part `2 a tau`: with a correction, from the
        spline of the total (see `__init__`), never a difference."""
        tau = float(tau)
        if self.corr is None:
            return 2.0 * self.a * tau + self.D(tau)
        if tau <= 0.0:
            return 0.0
        if tau < self.tlo:
            q2 = self.q2 if self.pc is not None else 0.0
            return (2.0 * self.a * tau + (q2 + self.cq2) * tau * tau
                    + sum(tl.structure(tau, self.pref) for tl in self.tails))
        if tau > self.thi:
            return 2.0 * self.a * tau + self.D(tau)
        return float(np.exp(self.tspline(np.log(tau))))

    def D(self, tau):
        """`D` without its white part (the power-law pieces', and the
        correction's past the grid; on the grid a correction's `D` is
        `Dtotal - 2 a tau`)."""
        tau = float(tau)
        if tau <= 0.0:
            return 0.0
        if self.corr is not None and self.tlo <= tau <= self.thi:
            return self.Dtotal(tau) - 2.0 * self.a * tau
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
        white part is never subtracted from the total."""
        if self.corr is None:
            x, y = -0.5 * self.D(tau), -0.5 * self.Dinf
            m = max(x, y)
            diff = np.exp(m - self.a * tau) * (-np.expm1(min(x, y) - m))
            return diff if x >= y else -diff
        x, y = -0.5 * self.Dtotal(tau), -self.a * tau - 0.5 * self.Dinf
        m = max(x, y)
        diff = np.exp(m) * (-np.expm1(min(x, y) - m))
        return diff if x >= y else -diff

    @property
    def line_weight(self):
        """The carrier line's weight when no white part broadens it."""
        return float(np.exp(-0.5 * self.Dinf))

    def whole(self, tau):
        """`exp(-D/2)` itself (with the white part), at most 1."""
        return float(np.exp(-0.5 * self.Dtotal(tau)))

    def __call__(self, f):
        w = 2.0 * np.pi * abs(float(f))
        a = self.a
        ## ⚠ A NEGATIVE `D_inf` (a signed correction that takes white noise
        ## away over a band holding the core): `exp(-D_inf/2) L_w` and the
        ## integral of `g` are then each `exp(|D_inf|/2)` large and cancel
        ## (e^33 against e^33 at a corner 0.03 linewidths out, every digit
        ## lost; 4.8e-7 split, against the closed form).  So the body
        ## integrates `exp(-D/2) <= 1` itself over the grid, and the
        ## Lorentzian enters as its closed-form tail past `thi`,
        ## `2 Re[e^{(-a + i w) thi} / (a - i w)]`, weighted `exp(-D_inf/2)`.
        ## A white part is then present (`D` grows as `2 a tau`, so `D_inf`
        ## < 0 needs `a > 0`).  Below `SPLIT_GAIN` the subtracted form is
        ## the more accurate one and stays.
        split = a > 0.0 and -0.5 * self.Dinf > math.log(SPLIT_GAIN)
        if split:
            z = complex(-a, w) * self.thi
            head = (np.exp(-0.5 * self.Dinf + z.real) * 2.0
                    * (complex(np.cos(z.imag), np.sin(z.imag)) / complex(a, -w)).real)
            body = self.whole
        else:
            lw = (2.0 * a / (a * a + w * w)) if a > 0.0 else 0.0
            head = np.exp(-0.5 * self.Dinf) * lw
            body = self.g
        nd = int(np.ceil(np.log10(self.thi / self.tlo)))
        edges = np.concatenate(([0.0], np.geomspace(self.tlo, self.thi,
                                                    nd + 1)))
        tot, err = 0.0, 0.0
        for lo, hi in zip(edges[:-1], edges[1:]):
            kw = ({'limit': 400} if w == 0.0 else
                  {'weight': 'cos', 'wvar': w, 'limit': 400})
            val, e = self.quad(body, lo, hi, **kw)
            tot += val
            err += e
        kw = ({'limit': 400} if w == 0.0 else
              {'weight': 'cos', 'wvar': w, 'limlst': 200})
        val, e = self.quad(self.g, self.thi, np.inf, **kw)
        tot += val
        err += e
        out = head + 2.0 * tot
        self.quad.note(2.0 * err, out)
        ## this value's own QUADPACK estimate, relative
        self.last_err = 2.0 * err / max(abs(out), 1e-300)
        return out


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
