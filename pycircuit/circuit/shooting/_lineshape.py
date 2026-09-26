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
    and a coloured `c_c` (`PowerPieces`), `D_c = pref int ...`,
    `pref = 4 i^2 f0^2`."""

    def __init__(self, a, pc, pref, per_decade=TAU_PER_DECADE):
        self.a = float(a)
        self.pc, self.pref = pc, float(pref)
        self.quad = _Quad()
        self.last_err = 0.0
        nu0, nuN = pc.nu[0], pc.nu[-1]
        self.Dinf = pref * pc.moment(-2.0, nu0, nuN)
        ## D_c ~ q2 tau^2 below the grid; the band edge rings above it
        self.q2 = pref * 2.0 * np.pi ** 2 * pc.moment(0.0, nu0, nuN)
        self.h0 = float(pc(nu0)) / nu0 ** 2
        self.tlo, self.thi = 1e-3 / nuN, 1e3 / nu0
        n = int(np.ceil(per_decade * np.log10(self.thi / self.tlo)))
        self.taus = np.geomspace(self.tlo, self.thi, n + 1)
        Ds = np.array([phase_structure(t, pc, pref, self.quad)
                       for t in self.taus])
        self.spline = interpolate.CubicSpline(np.log(self.taus),
                                              np.log(np.maximum(Ds, 1e-300)))

    def D(self, tau):
        tau = float(tau)
        if tau <= 0.0:
            return 0.0
        if tau < self.tlo:
            return self.q2 * tau * tau
        if tau > self.thi:
            beta = 2.0 * np.pi * tau
            return (self.Dinf + self.pref * self.h0
                    * np.sin(beta * self.pc.nu[0]) / beta)
        return float(np.exp(self.spline(np.log(tau))))

    def g(self, tau):
        """`e^{-a tau} (e^{-D/2} - e^{-D_inf/2})`, in logs (no 0 * inf)."""
        x, y = -0.5 * self.D(tau), -0.5 * self.Dinf
        m = max(x, y)
        diff = np.exp(m - self.a * tau) * (-np.expm1(min(x, y) - m))
        return diff if x >= y else -diff

    @property
    def line_weight(self):
        """The carrier line's weight when no white part broadens it."""
        return float(np.exp(-0.5 * self.Dinf))

    def __call__(self, f):
        w = 2.0 * np.pi * abs(float(f))
        a = self.a
        lw = (2.0 * a / (a * a + w * w)) if a > 0.0 else 0.0
        head = np.exp(-0.5 * self.Dinf) * lw
        nd = int(np.ceil(np.log10(self.thi / self.tlo)))
        edges = np.concatenate(([0.0], np.geomspace(self.tlo, self.thi,
                                                    nd + 1)))
        tot, err = 0.0, 0.0
        for lo, hi in zip(edges[:-1], edges[1:]):
            kw = ({'limit': 400} if w == 0.0 else
                  {'weight': 'cos', 'wvar': w, 'limit': 400})
            val, e = self.quad(self.g, lo, hi, **kw)
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
    expansion), not a bound."""
    f = abs(float(f))
    nx = f / np.sqrt(6.0)
    lo = min(max(nx, pc.nu[0]), pc.nu[-1])
    sig2 = 2.0 * i * i * f0 * f0 * (c_w * nx + pc.moment(0.0, pc.nu[0], lo))
    above = 2.0 * i * i * f0 * f0 * (c_w / nx + pc.moment(-2.0, lo, pc.nu[-1]))
    return 6.0 * sig2 / (f * f) + above
