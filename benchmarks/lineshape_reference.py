"""THE COLOURED LINESHAPE'S INDEPENDENT REFERENCE (2026-09-27): white c_w
plus a pure 1/f colour c_c(nu) = k / nu on [fmin, fmax], i = 1, f0 = 1
(normalised; a circuit's c scales as f0 c(f0 nu)):

    S(f) = 2 int_0^inf exp(-D(tau)/2) cos(2 pi f tau) dtau,
    D = 2 a tau + D_c,   a = 2 pi^2 c_w,
    D_c(tau) = 4 k int_fmin^fmax (1 - cos beta nu) / nu^3 dnu,  beta = 2 pi tau,

the colour's integral EXACT through Ci (by parts):
    F(nu) = -(1 - cos b nu)/(2 nu^2) - b sin(b nu)/(2 nu) + (b^2/2) Ci(b nu),
the transform by mpmath's quadosc at DPS digits.

Validated: the Ci form against direct quadrature, 1e-25 .. 1e-38; quadosc
against the Lorentzian-colour closed form, 1e-15.  30 .. 130 s per offset,
so the gate test (`test_the_coloured_lineshape_meets_an_independent_
reference_across_the_handover`) holds these values as constants; this
script regenerates them:

    python benchmarks/lineshape_reference.py

Regimes: a real oscillator's line (`_lc_osc`'s levels: c_w 1.7496e-8,
k 7.5414e-11, linewidth 5.5e-8 f0) and a broad one (c_w 1e-4) with the
flicker 1x and 100x the white at 1e-3 f0; fmin 1e-7, fmax 0.5.
"""
import mpmath as mp

DPS = 40


def Dc(tau, k, fmin, fmax):
    b = 2 * mp.pi * mp.mpf(tau)

    def F(nu):
        nu = mp.mpf(nu)
        return (-(1 - mp.cos(b * nu)) / (2 * nu ** 2)
                - b * mp.sin(b * nu) / (2 * nu) + b ** 2 / 2 * mp.ci(b * nu))
    return 4 * mp.mpf(k) * (F(fmax) - F(fmin))


def S(f, c_w, k, fmin, fmax, dps=DPS):
    with mp.workdps(dps):
        a = 2 * mp.pi ** 2 * mp.mpf(c_w)
        g = lambda t: mp.exp(-(2 * a * t + Dc(t, k, fmin, fmax)) / 2) if t > 0 else mp.mpf(1)
        w = 2 * mp.pi * mp.mpf(f)
        if f == 0:
            val = mp.quad(g, [0, 1 / mp.mpf(fmax), 1 / mp.mpf(fmin), mp.inf])
        else:
            val = mp.quadosc(lambda t: g(t) * mp.cos(w * t), [0, mp.inf], omega=w)
        return float(2 * val)


def S_panels(f, c_w, k, fmin, fmax, dps=30, panel=0.125):
    """`S` AROUND THE BAND'S LOWER EDGE (regime 'edge', 2026-10-01): `D`
    rings at `fmin` there, which `quadosc` (cut at the zeros of the offset's
    cosine) does not follow; the transform is taken on panels of `panel /
    fmin` (and at most a quarter of the offset's period) out to where
    `exp(-a tau)` is below 1e-20."""
    with mp.workdps(dps):
        a = 2 * mp.pi ** 2 * mp.mpf(c_w)
        w = 2 * mp.pi * mp.mpf(f)
        g = lambda t: (mp.exp(-(2 * a * t + Dc(t, k, fmin, fmax)) / 2) * mp.cos(w * t)
                       if t > 0 else mp.mpf(1))
        tend = 46 / a
        step = min(panel / mp.mpf(fmin), 0.25 / mp.mpf(f)) if f else panel / mp.mpf(fmin)
        pts = [mp.mpf(0)] + [mp.mpf(10) ** e / fmax for e in range(-2, 4)]
        pts = [p for p in pts if p < step] + [step]
        t = step
        while t < tend:
            t = t + step
            pts.append(t)
        return float(2 * mp.quad(g, pts))


#: the edge regime: a 1/f colour with `D_inf = 1.5` on [1e-5, 0.5], offsets
#: in units of fmin (`test_the_coloured_lineshape_around_the_band_edge_...`)
EDGE = (1e-7, 0.75e-10, 1e-5, (0.3, 0.55, 0.8, 1.0, 1.2, 2.0, 5.0, 20.0))

REGIMES = (('narrow', 1.7496e-8, 7.5414e-11, (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3)),
           ('broad r1', 1e-4, 1e-7, (3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3)),
           ('broad r100', 1e-4, 1e-5, (3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3)))


if __name__ == '__main__':
    for name, c_w, k, offs in REGIMES:
        print(name, [(o, S(o, c_w, k, 1e-7, 0.5)) for o in offs], flush=True)
    c_w, k, fmin, offs = EDGE
    print('edge', [(o, S_panels(o * fmin, c_w, k, fmin, 0.5)) for o in offs], flush=True)
