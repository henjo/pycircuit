"""Shooting tests: shooting lineshape.  Split out of test_analysis_shooting.py on
2026-09-27 (in its original order); shared helpers are in
`_shooting_fixtures.py`, the HDL elements in `_shooting_elements.py`.
"""
from pycircuit.circuit import *
from pycircuit.circuit.shooting import (PAC, algebraic_conditioning,
                                        topological_index)
import warnings
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter as _HdlParameter, white_noise)
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest
import functools as _functools
from pycircuit.circuit.tests._shooting_fixtures import (_Flicker,
    _a10_vdp,
    _a9_vdp,
    _adjoint_ladder,
    _coloured_vdp,
    _lc_osc,
    _orbit_modulated_vdp,
    _solve_vdp_noise,
    _vdp_with_noise)


def test_the_lorentzian_conserves_the_carrier_power_exactly():
    """⚠ THE INVARIANT THAT SEPARATES THIS FROM AN LTV TREATMENT.

    Noise spreads the carrier's power into a line of finite width; it does
    not create any. `∫ S_i df = 1` exactly, for every harmonic and every
    `c`. Analyses "based on linear time-invariant or linear time-varying
    concepts erroneously predict infinite noise power [at the carrier] as
    well as infinite total integrated power" — so this is the property
    that says the closed form is doing the nonlinear thing.

    Integrated numerically over the IMPLEMENTED function, not re-derived:
    a re-derivation would only be checking the algebra against itself.
    """
    from scipy.integrate import quad
    f0 = 150.0
    for c in (1e-9, 1e-7, 1e-5):
        for i in (1, 2, 5):
            tot, err = quad(lambda f: float(PAC.lorentzian(f, c, f0, i)),
                            -np.inf, np.inf, limit=400)
            assert abs(tot - 1.0) < 1e-6, \
                'harmonic %d at c=%g integrates to %.9f, not 1 — the ' \
                'carrier power is not conserved' % (i, c, tot)


def test_the_lineshape_is_lorentzian_where_it_should_be():
    """Finite at the carrier, `1/f²` far out, and the corner where predicted.

    Three properties, each of which a wrong constant would break
    differently: the peak is `1/(π² i² f₀² c)`, the far skirt falls as
    `1/f²` (a RATE, asserted as one), and the half-width is `π i² f₀² c`.
    """
    f0, c, i = 150.0, 1e-7, 1
    peak = float(PAC.lorentzian(0.0, c, f0, i))
    want = 1.0 / (np.pi ** 2 * i * i * f0 * f0 * c)
    assert abs(peak / want - 1.0) < 1e-12, \
        'peak %.6e against the analytic %.6e' % (peak, want)

    ## far out: doubling the offset must quarter the density
    far = [float(PAC.lorentzian(f, c, f0, i)) for f in (1e5, 2e5, 4e5)]
    for a, b in zip(far, far[1:]):
        assert abs(a / b - 4.0) < 0.02, \
            'the skirt falls %.3fx per doubling, not the 4x of 1/f²' % (a / b)

    ## half-width
    hw = np.pi * i * i * f0 * f0 * c
    assert abs(float(PAC.lorentzian(hw, c, f0, i)) / peak - 0.5) < 1e-12


def test_higher_harmonics_are_noisier_by_20log10i():
    """The skirt scales as `i²` and the corner as `i⁴`.

    So harmonic `i` sits `20 log₁₀(i)` dB above the fundamental far from
    the carrier — 6.02 dB for the second, 9.54 for the third. A designer
    reads that as "the divider makes it better, the multiplier worse", and
    it is a free consequence of the closed form rather than a separate
    calculation.
    """
    f0, c = 150.0, 1e-7
    far = 1e6
    base = float(PAC.lorentzian(far, c, f0, 1))
    for i in (2, 3, 5):
        db = 10.0 * np.log10(float(PAC.lorentzian(far, c, f0, i)) / base)
        want = 20.0 * np.log10(i)
        assert abs(db - want) < 0.05, \
            'harmonic %d is %.3f dB above the fundamental, not %.3f' \
            % (i, db, want)
        ## and its corner is i^4 wider
        hw_i = np.pi * i * i * f0 * f0 * c
        assert abs(float(PAC.lorentzian(hw_i, c, f0, i))
                   / float(PAC.lorentzian(0.0, c, f0, i)) - 0.5) < 1e-12


def test_the_oscillator_spectrum_is_built_and_refuses_a_driven_circuit():
    """End to end, and the one circuit class it does not describe."""
    import warnings
    _cir, pss, pac = _solve_vdp_noise()
    offs = np.array([1e-4, 1e-3, 1e-2, 1e-1])
    Sv, L = pac.oscillator_spectrum(pss, offs, 0, harmonic=1)
    assert np.all(np.isfinite(Sv)) and np.all(np.diff(Sv) < 0), \
        'the spectrum should be finite and falling with offset: %s' % Sv
    assert L[0] > L[-1], 'L(f) should fall with offset'
    ## far out it is 1/f^2, which in dB is -20 dB/decade
    slope = (L[-1] - L[-2]) / np.log10(offs[-1] / offs[-2])
    assert abs(slope + 20.0) < 1.0, \
        'the far skirt is %.2f dB/decade, not the -20 of 1/f²' % slope

    circuit.default_toolkit = circuit.numeric
    driven = _adjoint_ladder(3)
    p2 = PSS(driven, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p2.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    with pytest.raises(ValueError, match='FREE-RUNNING'):
        PAC(driven, toolkit=circuit.numeric).diffusion_constant(p2)


def test_the_phase_psd_convention_is_the_lorentzians_own():
    """⚠ NO SECOND CONVENTION IS INTRODUCED, and that is the whole point.

    A one-sided/two-sided slip already cost this class a factor of two, so
    `phase_psd` is not given an independent normalisation: it is checked
    against `lorentzian`'s far skirt, an object already gated by power
    conservation to 1.000000.

        S_phi,i(f) = i² f₀² (c + Γ(f)) / f²   and   lorentzian → i² f₀² c / f²

    ⚠ AND THE RESIDUAL IS FULLY ACCOUNTED FOR, which is stronger than it
    being small. The Lorentzian carries an `f_h²` term its skirt drops, so
    the disagreement must be exactly `(i²·corner/f)²` — quartic in the
    harmonic. Measured at offsets starting 1e3 corners out:

        harmonic     1        2        3
        max |1-r|  1.0e-06  1.6e-05  8.1e-05
        ratio        1        16       81      = 1 : 2⁴ : 3⁴

    A convention error would not reproduce that ratio.
    """
    _cir, pss, pac = _lc_osc()
    c = pac.diffusion_constant(pss)
    f0 = 1.0 / float(pss.period)
    corner = np.pi * f0 ** 2 * c
    offs = np.logspace(np.log10(corner * 1e3), np.log10(corner * 1e9), 7)
    errs = []
    for i in (1, 2, 3):
        S = pac.phase_psd(pss, offs, harmonic=i,
                          frequency_aware=False)  # the DC closed form's convention
        Lz = PAC.lorentzian(offs, c, f0, harmonic=i)
        errs.append(float(np.max(np.abs(S / Lz - 1.0))))
    assert errs[0] < 2e-6, 'harmonic 1 off by %.3e' % errs[0]
    for i, e in zip((2, 3), errs[1:]):
        assert abs(e / errs[0] / i ** 4 - 1.0) < 0.05, \
            'harmonic %d residual is %.3e, %.2fx the fundamental rather ' \
            'than the %d predicted by the Lorentzian curvature — the ' \
            'agreement is not the analytic one' % (i, e, e / errs[0], i ** 4)


class _StateDependentNoise(IS):
    """A source whose `CY` reads `x` — multiplicative noise, `G = G(x)`.

    ⚠ NOTHING IN THE TREE DOES THIS, which is why it exists here. Every
    shipped source has a constant `noisePSD`, and both compact MOS models
    have `CY` identically ZERO (no noise model at all).
    """

    def CY(self, x, w, epar=None):
        p = self.iparv.noisePSD * (1.0 + 0.5 * float(np.asarray(x).ravel()[0]))
        return self.toolkit.array([[p, -p], [-p, p]])


def test_multiplicative_noise_is_refused_only_where_the_sum_is_stationary():
    """⚠⚠ ONE GUARD COVERS TWO UNRELATED THEORETICAL HAZARDS.

    `_cy_reduced` samples `CY` at three states on the orbit and refuses a
    bias-dependent one. It was built for CYCLOSTATIONARITY: a
    bias-dependent `CY` correlates the sidebands through the window
    Fourier coefficients, so they stop adding in power and the stationary
    sum would be the wrong model.

    ⚠ IT ALSO CLOSES THE ITÔ/STRATONOVICH AMBIGUITY, WHICH IS A DIFFERENT
    QUESTION ENTIRELY. Demir ch.2: the Itô SDE `dX = f dt + G dW` and the
    Stratonovich one agree *"as long as `G(t,x) = G(t)` is independent of
    `x`"*; otherwise they are **two distinct Markov processes** differing
    *"in the systematic (drift) behavior but not in the fluctuational
    (diffusion) behavior"*. `CY = GGᵀ`, so a state-dependent `CY` is
    exactly a state-dependent `G` — and the guard refuses it.

    ⚠ SO THE SHIPPED CODE NEVER FACES THE INTERPRETATION CHOICE. Demir's
    own resolution is that the drift shift is *"on the order of the noise
    source intensity"* and *"for most practical physical systems the noise
    signals are small compared with the deterministic signals"* — a
    small-noise assumption. We do not need to lean on it, because the case
    where it matters raises instead.

    ⚠ AND THE FIXTURE POINT IS SHARPER THAN IT LOOKS: for ADDITIVE noise
    Itô and Stratonovich are IDENTICAL, so a suite whose sources are all
    state-independent could not detect an interpretation error even in
    principle. Ours are all additive. **The reason that is safe here is
    not the fixtures — it is that the code path does not exist**, which is
    a stronger position than an untested one and worth distinguishing.

    ⚠ THE TELL, IF IT EVER ARRIVES: the drift shifts by `½G∂ₓG` and the
    diffusion does not. A discrepancy in a MEAN but not in a VARIANCE is
    where to look.

    ⚠ AND THAT IS WHY THE DIFFUSION SURFACES NOW TAKE IT.  The routes that
    evaluate `CY` along the orbit report second-order statistics of the
    linearised response -- DIFFUSION, where the two interpretations agree
    -- and no MEAN: `oscillator_covariance` (2026-09-05),
    `pnoise(cyclostationary=True)`, and since 2026-09-26
    `diffusion_constant` (Demir's `B(x(t))`, `CY` at each PPV sample's
    state), `frequency_aware_diffusion` / `oscillator_spectrum`,
    `coloured_diffusion` and `modal_spectrum`.  On this fixture `c` matches the Lyapunov
    growth route (`CY` per step, no PPV) to +1.0e-3 on gear at 240 points
    (gear's own gap between the two; radau -4.4e-11).  The paths that sum
    ONE stationary `CY` still refuse.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
    cir['n'] = _StateDependentNoise('v', gnd, i=0.0, noisePSD=1e-6)
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=6.6634, timestep=6.6634 / 240,
                  x0=np.array([2.0, 0.0]), maxiterations=60)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    ## the paths that sum one stationary `CY` funnel through `_cy_reduced`
    for name, call in (
            ('orbital_correlation', lambda: pac.orbital_correlation(pss)),
    ):
        with pytest.raises(NotImplementedError, match='BIAS-DEPENDENT CY'):
            call()
    ## ⚠ NOT `oscillator_covariance` ANY MORE (2026-09-05): the covariance
    ## routes evaluate `CY` at every step, so a modulated source is inside
    ## their formulation and they RUN -- see `_lyapunov_pieces` and the
    ## switched-capacitor gate below.
    K, _d, _info = pac.oscillator_covariance(pss)
    assert np.all(np.isfinite(np.asarray(K, dtype=float)))
    ## nor `diffusion_constant` (2026-09-26): the same diffusion by the PPV
    c = float(pac.diffusion_constant(pss))
    assert abs(c / float(_info['c_from_growth']) - 1.0) < 3e-3, \
        (c, _info['c_from_growth'])
    ## nor `oscillator_spectrum` (its `frequency_aware_diffusion`, 2026-09-26)
    Sv, _L = pac.oscillator_spectrum(pss, [1e-3], 0)
    assert np.all(np.isfinite(Sv)) and np.all(np.asarray(Sv) > 0.0)
    ## nor `coloured_diffusion` (the l = 0 term of the modulated fold)
    assert np.all(np.isfinite(pac.coloured_diffusion(pss, [1e-3])))


def _carrier_power(pss, ov):
    X = np.asarray(pss.waveform[1], dtype=float)[ov][:-1]
    A1 = 2.0 * abs(np.fft.rfft(X)[1]) / len(X)
    return 0.5 * A1 * A1


def test_the_modal_spectrum_takes_a_white_source_that_follows_the_orbit():
    """`modal_spectrum` and `diffusion_constant` with a WHITE source whose
    level follows the orbit (2026-09-26; Andreas: "modal_spectrum with
    orbit-varying noise").  Both refused it (`_cy_reduced`).  The sum is
    now the P-form

        sum_{m,m'} T_m (P_{m'-m}/2) T_{m'}^H,   P_k = harmonics of CY(x(t))

    exact and without a square root (a root of `(k V)^2` is `|k V|`, whose
    kink spreads it over harmonics the sideband window cuts), and `c` is
    Demir's own `B(x(t))` form, `CY` at each PPV sample's state.

    Gated against the same physics built as a STATIONARY source times the
    tank voltage (`_orbit_modulated_vdp`), the verified stationary path,
    on an asymmetric orbit (a = 0.3, where the correlation cancels the
    parts to 1/40 of each), gear 400 points, at 0.3 / 3 / 10 / -10 f_amp:
    every part to <= 1.7e-13 (4.4e-13 at a = 0), `c` to 2.2e-16 -- gear's
    algebraic adjoint IS `k V_v q_v`, so the two discretisations coincide.
    And `c` against the LYAPUNOV route (`oscillator_covariance`, `CY` per
    step, no PPV) on radau: 4.6e-13.  Poisons (a = 0.3): the P-form cut to
    `P_0` (the stationary assumption) 0.42; `P_{m-m'}` for `P_{m'-m}` 0.58;
    `CY` read at one state 1.04 (`c` 0.39).  A byproduct: `pnoise(
    cyclostationary=True)` on the modulated circuit reproduces the
    stationary pnoise of its realisation digit for digit -- its first
    measurement on an oscillator.

    `frequency_aware_diffusion` too (2026-09-26), and with it the default
    `oscillator_spectrum`: `CY` at each sample of the frequency-aware PPV
    (the solve's own orbit, no twin).  Same pair, a = 0.3, where `c(f)`
    falls 40x over 0.3 .. 10 f_amp: `c(f)` 1.7e-14, `S_v` 1.2e-14 (radau
    3.4e-10 / 7.3e-10, its GMRES tolerance); `CY` at one state 0.39 ..
    0.51."""
    import warnings as _w
    res = {}
    for kind in ('white_ref', 'white'):
        _c, pss, pac, ov = _orbit_modulated_vdp(kind, a=0.3)
        f0 = 1.0 / float(pss.period)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            _v, info = pss.ppv()
            f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
            offs = np.array([0.3, 3.0, 10.0, -10.0]) * f_amp
            res[kind] = (pac.modal_spectrum(pss, offs, ov, H=8, sidebands=16),
                         float(pac.diffusion_constant(pss)),
                         np.array([pac.frequency_aware_diffusion(pss, o)
                                   for o in offs]),
                         np.asarray(pac.oscillator_spectrum(pss, offs, ov)[0]))
    A, B = res['white_ref'], res['white']
    for k in ('phase', 'orbital', 'correlation', 'total'):
        err = np.max(np.abs(B[0][k] / A[0][k] - 1.0))
        assert err < 1e-9, (k, B[0][k] / A[0][k] - 1.0)
    assert abs(B[1] / A[1] - 1.0) < 1e-12, (B[1], A[1])
    for i, name in ((2, 'frequency_aware_diffusion'), (3, 'oscillator_spectrum')):
        err = np.max(np.abs(B[i] / A[i] - 1.0))
        assert err < 1e-9, (name, B[i] / A[i] - 1.0)
    _c, pss, pac, ov = _orbit_modulated_vdp('white', method='radau')
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        _K, _d, oi = pac.oscillator_covariance(pss)
        c = float(pac.diffusion_constant(pss))
    assert abs(c / float(oi['c_from_growth']) - 1.0) < 1e-9, \
        (c, oi['c_from_growth'])

def _banded_lineshape_reference(f, a, eps, nu1, nu2):
    """The exact lineshape for `c(nu) = c_w [1 - eps (k2 - k1)]`, `k_j = 1 /
    (1 + (nu/nu_j)^2)`: a SIGNED change of `c` that vanishes at both ends.
    `D = 2 a tau - A1 (1 - e^{-b1 tau}) + A2 (1 - e^{-b2 tau})`, `A_j = a eps
    / (pi nu_j)`, `b_j = 2 pi nu_j`, so `exp(-D/2)` is a double series of
    exponentials and `S` one of Lorentzians (mpmath, 40 digits)."""
    import mpmath as mp
    mp.mp.dps = 40
    a, eps, nu1, nu2 = (mp.mpf(v) for v in (a, eps, nu1, nu2))
    w = 2 * mp.pi * abs(mp.mpf(f))
    A1, A2 = a * eps / (mp.pi * nu1), a * eps / (mp.pi * nu2)
    b1, b2 = 2 * mp.pi * nu1, 2 * mp.pi * nu2
    tot, ck = mp.mpf(0), mp.mpf(1)
    for k in range(400):
        if k:
            ck *= -A1 / 2 / k
        inner, cl = mp.mpf(0), mp.mpf(1)
        for l in range(400):
            if l:
                cl *= A2 / 2 / l
            r = a + k * b1 + l * b2
            t = cl * r / (r * r + w * w)
            inner += t
            if l > 5 and abs(t) < mp.mpf(10) ** -35 * abs(inner):
                break
        tot += ck * inner
        if k > 5 and abs(ck * inner) < mp.mpf(10) ** -35 * abs(tot):
            break
    return float(2 * mp.e ** ((A1 - A2) / 2) * tot)


def test_the_lineshape_takes_a_signed_correction_to_c_against_a_closed_form():
    """A SIGNED change of `c(nu)` inside the lineshape (2026-09-26, for the
    frequency-aware PPV to all orders): `_lineshape.SignedTable` +
    `correction_structure` put its structure function into `D` beside the
    power-law pieces, and `LogChebyshev` represents `rho = c_fa/c_dc - 1`
    from a few solves.  Against `_banded_lineshape_reference` (white `c_w`
    with a band taken away, or added), 0 .. 1000 linewidths:

        eps    corners (linewidths)   worst       first order (the old path)
        0.9      0.3 / 30             2.1e-7      -0.57 at the carrier, -4.55
                                                  at 0.1 (a NEGATIVE line)
        -0.5     0.3 / 30             2.5e-8
        0.99     0.03 / 3             4.8e-7      (D_inf = -65: the transform
                                                  integrates exp(-D/2) itself)

    ⚠ The third row is the cancellation the old transform could not take:
    `exp(-D_inf/2) L_w` and the integral of `g` are each e^33 there.
    Before `SPLIT_GAIN` it returned -75x the true value, flagged by its own
    QUADPACK estimate (1e6).  The Chebyshev fit of the analytic `rho`:
    129 points, 1.1e-7."""
    from pycircuit.circuit.shooting import _lineshape
    cw = 1e-3
    a = 2 * np.pi ** 2 * cw
    fc = a / (2 * np.pi)
    for eps, r1, r2, tol in ((0.9, 0.3, 30.0, 1e-6), (-0.5, 0.3, 30.0, 1e-6),
                             (0.99, 0.03, 3.0, 2e-6)):
        nu1, nu2 = r1 * fc, r2 * fc
        delta = lambda v, eps=eps, nu1=nu1, nu2=nu2: -cw * eps * (
            1.0 / (1.0 + (v / nu2) ** 2) - 1.0 / (1.0 + (v / nu1) ** 2))
        tab = _lineshape.SignedTable(delta, 1e-9 * nu1, 1e3 * nu2)
        shape = _lineshape.ColouredLineshape(a, None, 4.0, corr=tab)
        for x in (0.0, 0.1, 1.0, 10.0, 100.0, 1000.0):
            ref = _banded_lineshape_reference(x * fc, a, eps, nu1, nu2)
            assert abs(shape(x * fc) / ref - 1.0) < tol, (eps, x, shape(x * fc) / ref - 1.0)
        if eps == 0.9:
            first = 2 * a / (a * a + (2 * np.pi * 0.1 * fc) ** 2) + delta(0.1 * fc) / (0.1 * fc) ** 2
            assert first < 0.0, 'the first-order path was the case this replaces'
    nu1, nu2 = 0.3 * fc, 30.0 * fc
    rho = lambda v: -0.9 * (1.0 / (1.0 + (v / nu2) ** 2) - 1.0 / (1.0 + (v / nu1) ** 2))
    ch = _lineshape.LogChebyshev(rho, 1e-6 * nu1, 1e3 * nu2)
    vs = np.geomspace(1e-6 * nu1, 1e3 * nu2, 333)
    assert ch.converged and np.max(np.abs(ch(vs) - rho(vs))) < 1e-6


def test_the_rational_rho_fit_and_its_guards():
    """`_lineshape.RationalRho` (2026-09-27): the all-orders lineshape's
    correction `rho(nu)` as one barycentric rational (set-valued AAA) from
    adaptively chosen samples, ~25 where the Chebyshev series takes ~65 --
    the tolerance RELATIVE to `1 + rho`.  On analytic functions:

      * the banded Lorentzian's `rho` (rational): 20 samples, 7.8e-13;
      * a WEAK resonance (1e-5 of it, width 0.03, on a half-decade point):
        the adaptive search alone stops at 13 samples 7.5e-5 off -- the
        stop is FOOLED -- and the half-decade VERIFICATION catches it,
        2.1e-10 at 27 (fooled in 11 of 18 such cases, caught in all);
      * a GENUINE pole in the band: refused ('pole'); with the guard off
        the fit keeps it.
    The guards' switches exist for this test alone."""
    from pycircuit.circuit.shooting import _lineshape
    nu_a, nu_b, fc = 0.5e-8, 0.5, 1e-4
    base = lambda v: -0.9 * (1.0 / (1.0 + (v / (30 * fc)) ** 2)
                             - 1.0 / (1.0 + (v / (0.3 * fc)) ** 2))
    dense = np.geomspace(nu_a, nu_b, 20000)

    def rel(fit, f):
        return float(np.max(np.abs(fit(dense) - f(dense))
                            / np.maximum(np.abs(1.0 + f(dense)), 1e-4)))
    vec = lambda v: np.array([base(v), 0.0])
    fit = _lineshape.RationalRho(vec, nu_a, nu_b)
    assert fit.converged and fit.calls <= 25, (fit.converged, fit.calls)
    got = fit(dense)
    assert got.shape == (2, dense.size)
    assert np.max(np.abs(got[0] - base(dense))) < 1e-9
    vb = nu_b * 10 ** -3.5
    g = 0.03 * vb
    bump = lambda v: base(v) + 1e-5 * g * g / ((v - vb) ** 2 + g * g)
    fit = _lineshape.RationalRho(bump, nu_a, nu_b)
    assert fit.converged and rel(fit, bump) < 1e-6, rel(fit, bump)
    _lineshape.RationalRho.VERIFY = False
    try:
        blind = _lineshape.RationalRho(bump, nu_a, nu_b)
    finally:
        _lineshape.RationalRho.VERIFY = True
    assert rel(blind, bump) > 1e-5, 'the verification no longer binds'
    pole = lambda v: base(v) + 1e-6 / (v / 0.7e-3 - 1.0)
    fit = _lineshape.RationalRho(pole, nu_a, nu_b)
    assert not fit.converged and fit.reason == 'pole', (fit.converged, fit.reason)
    _lineshape.RationalRho.POLE_GUARD = False
    try:
        kept = _lineshape.RationalRho(pole, nu_a, nu_b)
    finally:
        _lineshape.RationalRho.POLE_GUARD = True
    assert kept.converged, 'the pole guard no longer binds'

def _fa_core_oscillator(psd, flicker_rel=1e-8):
    """The slow-node LC (tau = 100 T) with white `psd` AND a 1/f source at
    the slow node, `flicker_rel` of it at f0.  That is just above the 1e-9
    `_coloured_present` asks of a colour; with fmin = 1e-5 f0 its
    `D_c(inf)` is ~0.01, so the line is the white one.  psd = 0.7 puts the
    slow corner 10 linewidths from the core; the math is linear in the
    level."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0 = 6.6634
    c = SubCircuit()
    c.add_node('v'); c.add_node('w'); c.add_node('x')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', 'x', L=1.0)
    c['Rl'] = R('x', gnd, r=0.2)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: (u - u ** 3 / 3.0) + 0.25 * (u ** 2 - 2.0))
    c['Rs'] = R('v', 'w', r=1e2, noisy=False)
    c['Cs'] = C('w', gnd, c=100.0 * T0 / 1e2)
    c['nw'] = IS('w', gnd, i=0.0, noisePSD=psd)
    c['nf'] = _Flicker('w', gnd, i=0.0, noisePSD=psd * flicker_rel, fref=1.0 / T0)
    pss = PSS(c, method='gear', reltol=1e-11)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T0, timestep=T0 / 240, x0=x0, maxiterations=200)
    assert pss.converged
    return pss, PAC(c, toolkit=circuit.numeric), [str(n) for n in c.nodes].index('v')


def test_the_frequency_aware_lineshape_goes_to_all_orders_where_the_first_does_not_hold():
    """The frequency-aware coloured lineshape (2026-09-26): FIRST ORDER in
    the change while its estimated error is below `FA_FIRST_ORDER_TOL`,
    to ALL ORDERS above it (`_fa_lineshape`).

      * an LC with no slow node (`_lc_osc`, flicker + white): first order,
        estimate 2.4e-6, the probes' 7 solves (all orders would take 42
        solves and 21.6 s against 5.9 s for a ~1e-6 change);
      * slow corner 10 linewidths out (psd 0.7): all orders.  The first order read
        8.4 % low at the carrier and 13 % at 10 linewidths.
        - The probe estimate of `D_corr(inf)` is within 6 % of the full one.
        - The Chebyshev `rho` takes 65 solves.
        - The skirt meets the frequency-aware linear one where the handover
          takes it.
    """
    import warnings as _w
    pss, pac, ov = _fa_core_oscillator(0.7)
    f0 = 1.0 / float(pss.period)
    ## the carrier's ONE-SIDED power 2|X|^2: `S_v` over it is `L(f)`
    X2 = 2.0 * abs(pac.carrier_phasor(pss, ov, 1)) ** 2
    fcore = np.pi * f0 * f0 * pac._colour_fold(pss, 1e-5 * f0, None, 'x').c_white
    offs = np.array([0.0, 1.0, 10.0, 100.0, 300.0, 1000.0]) * fcore
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        S_all = pac.oscillator_spectrum(pss, offs, ov, fmin=1e-5 * f0)[0] / X2
        info = dict(pac.lineshape_info)
    ## ⚠ the far skirt resolved: past fmax the white part is held at its
    ## corrected level (`ConstantTail`).  Returned to `c_w` there, its edge
    ## rang through `D` and the two tau densities parted 4.7e-3 at 100
    ## linewidths (warned); held, 6.8e-5, and 1000 linewidths is the
    ## frequency-aware `S_phi` to 3.9e-6.  ⚠ And at 300 linewidths the
    ## ESTIMATE: against half the density it read 2e-3 (warned) for a true
    ## 7.8e-6; one density at two phases reads 1.5e-5 (2026-09-27)
    assert not [r for r in rec if 'estimated relative error' in str(r.message)], \
        [str(r.message) for r in rec]
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        sphi = float(pac.phase_psd(pss, [offs[-1]])[0])
    assert abs(S_all[-1] / sphi - 1.0) < 1e-4, S_all[-1] / sphi - 1.0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        S_first = pac.oscillator_spectrum(pss, offs, ov, fmin=1e-5 * f0,
                                          all_orders=False)[0] / X2
    assert info['frequency_aware'] == 'all orders', info
    assert abs(info['cinf_estimate'] / info['cinf'] - 1.0) < 0.1, info
    assert info['rho_fit'] == 'rational' and info['rho_err'] < 1e-6, info
    assert info['solves'] <= 40, info
    ## the rational fit against the Chebyshev series (both relative to 1 +
    ## rho; the series takes ~139 solves): measured 5.9e-8
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pac.FA_RHO_FIT = 'chebyshev'
        try:
            S_ch = pac.oscillator_spectrum(pss, offs, ov, fmin=1e-5 * f0)[0] / X2
        finally:
            del pac.FA_RHO_FIT
    ## ⚠ vacuous if the knob is not read on the instance
    assert pac.lineshape_info['rho_fit'] == 'chebyshev', pac.lineshape_info
    assert np.max(np.abs(S_all / S_ch - 1.0)) < 1e-6, S_all / S_ch - 1.0
    first_err = (S_first / S_all - 1.0)[:3]
    assert np.all(first_err < -0.05), first_err
    ## the first-order path's carrier gap is its core weight, exp(-D_corr/2)
    assert abs((1.0 + first_err[0]) * np.exp(-0.5 * info['cinf']) - 1.0) < 0.03, \
        (first_err[0], info['cinf'])
    ## and a line with no slow path near its core stays first order
    _c, pss, pac = _lc_osc(a=0.25, rs=0.2, flicker=True, psd=1e-6,
                           fref=1.0 / 6.66, white=1e-6)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pac.oscillator_spectrum(pss, [0.0, 1e-3 / float(pss.period)], 0,
                                fmin=1e-7 / float(pss.period))
    assert pac.lineshape_info['frequency_aware'] == 'first order', pac.lineshape_info
    assert pac.lineshape_info['estimate'] < pac.FA_FIRST_ORDER_TOL, pac.lineshape_info


def test_the_white_lineshape_takes_the_frequency_aware_ppv_to_all_orders_when_asked():
    """`oscillator_spectrum(all_orders=True)` on a WHITE source (2026-09-27;
    Andreas: "Fix the oscillator_spectrum but put it off by default").  The
    default stays the Lorentzian with `c(f)` per offset, which is first
    order in the frequency-aware change: the core keeps its DC weight.  On
    the slow-node LC with the corner 10 linewidths out (psd 0.7, white
    only), against all orders:

        linewidths      0        1        10       100      1000
        default      -8.4 %   -7.7 %   -11 %    -0.95 %  -8.8e-5

    All orders takes 77 bordered solves, 10.8 s, against 0.7 s.  The
    cross-check is the COLOURED path (`_fa_lineshape`, `rho` from the fold's
    parts, not `frequency_aware_diffusion`) on the same line with a
    vanishing flicker.  After the flicker's own effect, taken at DC, is
    removed, the two agree to 2.3e-6 at every offset."""
    import warnings as _w
    pss, pac, ov = _fa_core_oscillator(0.7, flicker_rel=0.0)
    f0 = 1.0 / float(pss.period)
    ## the carrier's ONE-SIDED power 2|X|^2: `S_v` over it is `L(f)`
    X2 = 2.0 * abs(pac.carrier_phasor(pss, ov, 1)) ** 2
    fcore = np.pi * f0 * f0 * pac.diffusion_constant(pss)
    offs = np.array([0.0, 1.0, 10.0, 100.0, 1000.0]) * fcore
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        S_def = pac.oscillator_spectrum(pss, offs, ov)[0] / X2
        assert pac.lineshape_info['frequency_aware'] == 'first order'
        S_all = pac.oscillator_spectrum(pss, offs, ov, all_orders=True)[0] / X2
        info = dict(pac.lineshape_info)
        S_dc = pac.oscillator_spectrum(pss, offs, ov, frequency_aware=False)[0] / X2
        sphi = float(pac.phase_psd(pss, [offs[-1]])[0])
    assert info['frequency_aware'] == 'all orders', info
    assert not [r for r in rec if 'estimated relative error' in str(r.message)], \
        [str(r.message) for r in rec]
    err = S_def / S_all - 1.0
    assert np.all(err[:3] < -0.05), err
    ## the default's carrier gap is its core weight, exp(-D_corr(inf)/2)
    assert abs((1.0 + err[0]) * np.exp(-0.5 * info['cinf']) - 1.0) < 0.03, \
        (err[0], info['cinf'])
    assert abs(S_all[-1] / sphi - 1.0) < 2e-4, S_all[-1] / sphi - 1.0
    ## the coloured path, with the flicker's own (DC) effect taken out
    pss2, pac2, ov2 = _fa_core_oscillator(0.7, flicker_rel=1e-8)
    X22 = 2.0 * abs(pac2.carrier_phasor(pss2, ov2, 1)) ** 2
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        S_col = pac2.oscillator_spectrum(pss2, offs, ov2, fmin=1e-4 * f0,
                                         all_orders=True)[0] / X22
        S_col_dc = pac2.oscillator_spectrum(pss2, offs, ov2, fmin=1e-4 * f0,
                                            frequency_aware=False)[0] / X22
        ## and `all_orders=False` forces the coloured first order
        pac2.oscillator_spectrum(pss2, offs[:2], ov2, fmin=1e-4 * f0,
                                 all_orders=False)
    assert pac2.lineshape_info['frequency_aware'] == 'first order', pac2.lineshape_info
    gap = (S_col / S_all - 1.0) - (S_col_dc / S_dc - 1.0)
    assert np.max(np.abs(gap)) < 1e-5, gap

def test_the_oscillator_spectrum_takes_a_coloured_source():
    """`oscillator_spectrum` with a 1/f source (2026-09-26; Andreas: "Do as
    you suggest").  The phase is then no Wiener process and the line no
    Lorentzian: the lineshape is the transform of `exp(-D(tau)/2)`, `D`
    built from `c(f)` -- the white part in closed form, the coloured part
    over `[fmin, fmax]` (`fmin` required: a 1/f^3 phase has no stationary
    lineshape without a low cutoff).  On the asymmetric, lossy LC
    (`_lc_osc`, flicker = white at the fixture's reference frequency):

      * the core WIDENS as `fmin` falls -- S(0) 1.9e7 / 1.2e5 / 9.0e4 /
        7.6e4 at fmin 1e-5 .. 1e-8 f0 (the white Lorentzian's is 4.0e7);
      * the skirt meets the LINEAR one (`phase_psd`) as the transform's
        own estimate says it should: S / S_phi - 1 = 7.5e-2 / 6.5e-3 /
        6.9e-4 / 7.6e-5 at 3e-4 / 1e-3 / 3e-3 / 1e-2 f0 against estimates
        7.6e-2 / 8.1e-3 / 1.1e-3 / 1.3e-4; further out, where the transform
        is cancellation-limited (+-2e-3 at 0.1 f0), the linear skirt is
        returned;
      * a flicker 1e-8 / 1e-7 of the white moves the line centre by
        -6.4e-5 / -6.4e-4, LINEAR in that level (9.996x), and the skirt by
        < 1e-6 (the colour enters additively);
      * no `fmin`, `frequency_aware=True`, `fmax > f0/2`: refused."""
    import warnings as _w
    _c, pss, pac = _lc_osc(a=0.25, rs=0.2, flicker=True, psd=1e-6,
                           fref=1.0 / 6.66, white=1e-6)
    f0 = 1.0 / float(pss.period)
    ## the carrier's ONE-SIDED power 2|X|^2: `S_v` over it is `L(f)`
    X2 = 2.0 * abs(pac.carrier_phasor(pss, 0, 1)) ** 2
    offs = np.array([0.0, 1e-3, 0.1]) * f0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        S = {fr: pac.oscillator_spectrum(pss, offs, 0, fmin=fr * f0)[0] / X2
             for fr in (1e-5, 1e-7)}
        sphi = pac.phase_psd(pss, offs[1:])
    assert S[1e-7][0] < 0.5 * S[1e-5][0], 'the core did not widen'
    dev = S[1e-7][1] / sphi[0] - 1.0
    assert 3e-3 < dev < 1e-2, dev
    assert abs(S[1e-7][2] / sphi[1] - 1.0) < 1e-6, S[1e-7][2] / sphi[1]
    for bad, exc in (({}, NotImplementedError),
                     ({'fmin': 1e-7 * f0, 'fmax': f0}, ValueError)):
        with pytest.raises(exc):
            pac.oscillator_spectrum(pss, offs, 0, **bad)
    devs = []
    for psd in (1e-14, 1e-13):
        _c, pss, pac = _lc_osc(a=0.25, rs=0.2, flicker=True, psd=psd,
                               fref=1.0 / 6.66, white=1e-6)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            ## the colour's additivity onto the DC Lorentzian (the
            ## frequency-aware skirt moves off it on this asymmetric orbit)
            Sv = pac.oscillator_spectrum(pss, offs, 0, fmin=1e-7 * f0,
                                         frequency_aware=False)[0]
            cw = pac._colour_fold(pss, 1e-7 * f0, None, 'x').c_white
        Lw = pac.lorentzian(offs, cw, f0, 1)
        ## (over the carrier's one-sided power 2|X|^2)
        devs.append(Sv / (2.0 * abs(pac.carrier_phasor(pss, 0, 1)) ** 2) / Lw - 1.0)
    ## the line centre moves LINEARLY with the flicker level (9.996x for
    ## 10x; -6.4e-4 at 1e-7 of the white -- the 1/f^3 phase wanders at
    ## the lags that set the core); the skirt stays at the Lorentzian
    assert abs(devs[1][0] / devs[0][0] - 10.0) < 0.1, devs
    assert np.max(np.abs(devs[1][1:])) < 1e-6, devs[1]


def test_the_coloured_lineshape_takes_a_source_that_follows_the_orbit():
    """The coloured lineshape reads `c(f)` from the same fold as `phase_psd`,
    modulated sources included: a SIGNED flicker `k V_v flicker_noise(1)` on
    van der Pol against the same physics as a stationary 1/f source times
    V_v (`_orbit_modulated_vdp`), at 0 .. 1e-2 f0: 1.9e-13.

    ⚠ ON THE CHEBYSHEV `rho` (`FA_RHO_FIT = 'chebyshev'`, 2026-09-27).  Its
    nodes are fixed, so two inputs 1e-13 apart stay 1e-13 apart.  The
    default rational fit picks its samples adaptively, and a 1e-13
    difference can change a greedy choice: the two setups then agree to
    7.8e-9, each correct to its fit (2.3e-7 on `c_fa`, verified).  This
    test is about the FOLD's equivalence, so it takes the fixed nodes."""
    import warnings as _w
    res = {}
    for kind in ('flicker_ref', 'flicker'):
        _c, pss, pac, ov = _orbit_modulated_vdp(kind)
        f0 = 1.0 / float(pss.period)
        pac.FA_RHO_FIT = 'chebyshev'
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res[kind] = pac.oscillator_spectrum(
                pss, np.array([0.0, 1e-4, 1e-2]) * f0, ov, fmin=1e-7 * f0)[1]
        ## ⚠ vacuous if the knob is not read on the instance
        assert pac.lineshape_info['rho_fit'] == 'chebyshev', pac.lineshape_info
    err = np.max(np.abs(10 ** ((res['flicker'] - res['flicker_ref']) / 10)
                        - 1.0))
    assert err < 1e-9, err


def _slow_node_oscillator(kind, tau_over_T=100.0):
    """The A2 fixture: an asymmetric, lossy LC whose noise source sits
    behind a slow RC node `w` (tau = 100 T): 'white', 'lorentz' or
    'flicker'."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0 = 6.6634
    c = SubCircuit()
    c.add_node('v'); c.add_node('w'); c.add_node('x')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', 'x', L=1.0)
    c['Rl'] = R('x', gnd, r=0.2)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: (u - u ** 3 / 3.0) + 0.25 * (u ** 2 - 2.0))
    c['Rs'] = R('v', 'w', r=1e2, noisy=False)
    c['Cs'] = C('w', gnd, c=tau_over_T * T0 / 1e2)
    c['n'] = {'white': lambda: IS('w', gnd, i=0.0, noisePSD=1e-6),
              'lorentz': lambda: IS('w', gnd, i=0.0, noisePSD=1e-6,
                                    noiseTau=0.3 * T0),
              'flicker': lambda: _Flicker('w', gnd, i=0.0, noisePSD=1e-6,
                                          fref=1.0 / T0)}[kind]()
    pss = PSS(c, method='gear', reltol=1e-11)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T0, timestep=T0 / 240, x0=x0, maxiterations=200)
    assert pss.converged
    return pss, PAC(c, toolkit=circuit.numeric), [str(n) for n in c.nodes].index('v')


def test_phase_psd_is_frequency_aware_for_a_coloured_source_behind_a_slow_node():
    """⚠ The coloured folds took the DC PPV until 2026-09-26 (Andreas:
    "Continue as you suggest"; default ON, his call): a source behind a
    slow path moves the phase through that path's filter, which the DC PPV
    -- the response to a perturbation slow against every mode -- cannot
    see.  `coloured_diffusion_resolved` / `phase_psd` now read `V_l` from
    `frequency_aware_ppv(f)` per offset, the bands unchanged.  Against
    pnoise's PM content, `pm / (4 |X|^2 S_phi)`, behind a tau = 100 T node
    at 1e-3 / 1e-2 f0:

        source      frequency-aware      DC PPV
        white       1.0001  1.0002       0.7291  0.0268
        Lorentzian  1.0001  1.0003       0.7290  0.0264
        1/f         1.0001  1.0003       0.7290  0.0263

    For a white source the fold IS `frequency_aware_diffusion(f)` (4.7e-10
    here).  The AM-to-PM control (van der Pol C = 4, a = 0.3, a Lorentzian
    at the tank): 0.9954 / 0.9588 at 1 / 10 f_amp against DC 0.6510 /
    0.3077 (white: 0.998 / 0.980 against 0.647 / 0.302)."""
    import warnings as _w
    ## (pnoise's PM costs 13 s an offset: the flicker source at 1e-2 f0,
    ## where the two separate most; white is gated by the identity, and by
    ## `test_oscillator_spectrum_is_frequency_aware_above_the_slow_corner`)
    pss, pac, ov = _slow_node_oscillator('white')
    f0 = 1.0 / float(pss.period)
    offs = np.array([1e-3, 1e-2]) * f0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        cfa = pac.coloured_diffusion_resolved(pss, offs)
        fad = np.array([pac.frequency_aware_diffusion(pss, o) for o in offs])
    assert np.max(np.abs(cfa / fad - 1.0)) < 1e-8, cfa / fad
    pss, pac, ov = _slow_node_oscillator('flicker')
    ## the carrier's ONE-SIDED power 2|X|^2: `S_v` over it is `L(f)`
    X2 = 2.0 * abs(pac.carrier_phasor(pss, ov, 1)) ** 2
    o = 1e-2 * f0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        sfa = float(pac.phase_psd(pss, [o])[0])
        sdc = float(pac.phase_psd(pss, [o], frequency_aware=False)[0])
        pm = pac.am_pm_noise(pss, o, ov, carrier=1, maxsidebands=16)[1]
    assert abs(pm / (X2 * sfa) - 1.0) < 5e-3, pm / (X2 * sfa)
    assert pm / (X2 * sdc) < 0.03, pm / (X2 * sdc)
    ## and the coloured LINESHAPE (`oscillator_spectrum`), frequency-aware
    ## by default: its skirt there is the frequency-aware `S_phi`, to
    ## second order
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        L_fa = pac.oscillator_spectrum(pss, [0.0, o], ov, fmin=1e-7 * f0)[0] / X2
        info = dict(pac.lineshape_info)
        L_dc = pac.oscillator_spectrum(pss, [0.0, o], ov, fmin=1e-7 * f0,
                                       frequency_aware=False)[0] / X2
    ## the skirt is the frequency-aware S_phi PLUS its second-order
    ## correction (2026-09-27, `_lineshape.second_order_skirt`): 2.1e-4 here,
    ## the size the linear skirt was measured to miss by at this offset on a
    ## real oscillator's line (6.5e-5 at 1e-2 f0, against mpmath).  Until then
    ## the linear skirt was returned, identical to S_phi to 1e-6
    assert 1e-5 < abs(L_fa[1] / sfa - 1.0) < 1e-3, L_fa[1] / sfa
    ## the correction reaches this 1/f core (|D_corr(inf)|/2 ~ 5e-4 > the
    ## first order's tolerance), so it is taken to all orders and the
    ## carrier moves by about that (2.2e-4, measured); the first order left
    ## it untouched
    assert info['frequency_aware'] == 'all orders', info
    assert 1e-5 < abs(L_fa[0] / L_dc[0] - 1.0) < 1e-3, L_fa[0] / L_dc[0] - 1.0
    assert L_dc[1] / L_fa[1] > 30.0, L_dc[1] / L_fa[1]
    ## AM-to-PM above f_amp, a coloured source at the core
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=4.0)
    cir['L'] = L('v', gnd, L=0.25)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.3 * u * u)
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6, noiseTau=0.3 * T)
    pss = PSS(cir, method='gear', reltol=1e-12)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 400, x0=np.array([2.0, 0.0]),
                  maxiterations=300)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        f0 = 1.0 / float(pss.period)
        _v, info = pss.ppv()
        f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
        o = f_amp
        ## the carrier's ONE-SIDED power 2|X|^2: `S_v` over it is `L(f)`
        X2 = 2.0 * abs(pac.carrier_phasor(pss, 0, 1)) ** 2
        pm = pac.am_pm_noise(pss, o, 0, carrier=1, maxsidebands=16)[1]
        rfa = pm / (X2 * float(pac.phase_psd(pss, [o])[0]))
        rdc = pm / (X2 * float(pac.phase_psd(pss, [o],
                                                 frequency_aware=False)[0]))
    assert abs(rfa - 1.0) < 0.02, rfa
    assert rdc < 0.7, rdc


def test_a_coloured_source_folds_per_harmonic_and_agrees_with_pnoise():
    """`phase_psd` on a coloured source is `pnoise/P_carrier`, per harmonic.

    `pnoise` folds `CY` at the source-side frequency `f - l f_0` for each
    sideband (A3's design), so it is the reference the resolved fold must
    meet -- and it does, to 2e-3 at `df/f0 = 1e-3` and 3e-4 at `1e-4`,
    the same agreement the WHITE source shows (the control row), which
    says the residual is the sideband discretisation and not the colour.

    The white-only routines now REFUSE the coloured source with the
    reason, instead of folding `CY` at one frequency and returning a
    plausible number: `diffusion_constant`, and `oscillator_spectrum`
    through it (its closed form is exact for white only, by its own
    docstring).
    """
    for kind, tol in (('white', 3e-3), ('coloured', 3e-3)):
        _c, pss, pac, ov = _coloured_vdp(kind)
        f0 = 1.0 / float(pss.period)
        Pc = _carrier_power(pss, ov)
        ks = (1e-3, 1e-4)
        offs = np.array(ks) * f0
        sphi = pac.phase_psd(pss, offs)
        for k, o, s in zip(ks, offs, sphi):
            pn = float(np.real(pac.pnoise(pss, f0 * (1.0 + k), ov)[0])) / Pc
            assert abs(s / pn - 1.0) < tol, \
                '%s df/f0=%g: phase_psd %.6e vs pnoise/Pc %.6e' \
                % (kind, k, s, pn)
        if kind == 'coloured':
            with pytest.raises(NotImplementedError, match='COLOURED'):
                pac.diffusion_constant(pss)
            with pytest.raises(NotImplementedError, match='COLOURED'):
                pac.oscillator_spectrum(pss, offs, ov)
        else:
            c = pac.diffusion_constant(pss)
            assert np.allclose(sphi, (f0 ** 2) * c / offs ** 2, rtol=1e-12), \
                'for white the spectrum is the Lorentzian skirt in c exactly'


def test_phase_noise_stack_works_over_trbdf2():
    """ppv, the diffusion constant, and the oscillator spectrum all run over
    the TR-BDF2 monodromy and agree with Gear-2.

    The autonomous phase-noise surfaces ride on the PPV, which rides on the
    monodromy transpose -- so a correct two-stage adjoint makes the whole
    stack available without a Gear-2 twin. On van der Pol with a white
    source the diffusion constant `c` and the lineshape must match the
    Gear-2 numbers to O(h^2).
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def solve(method):
        cir = _vdp_with_noise(1e-6)
        m = cir.n - 1
        pss = PSS(cir, method=method, reltol=1e-12)
        x0 = np.zeros(m); x0[0] = 2.0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=6.6634, timestep=6.6634 / 240, x0=x0,
                      maxiterations=60)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        c = pac.diffusion_constant(pss)
        _Sv, L = pac.oscillator_spectrum(pss, [1e-2, 1e-1, 1.0], 0,
                                         harmonic=1)
        return c, np.asarray(L, dtype=float)

    c_g, L_g = solve('gear')
    c_t, L_t = solve('trbdf2')
    assert abs(c_t - c_g) < 1e-2 * c_g, (c_t, c_g)
    assert np.max(np.abs(L_t - L_g)) < 0.05, (L_t, L_g)


def test_phase_noise_stack_works_over_radau():
    """ppv, the diffusion constant, and the oscillator spectrum all run over
    the Radau IIA(3) monodromy and agree with Gear-2.

    The autonomous phase-noise surfaces ride on the PPV, which rides on the
    monodromy transpose and the coupled forced adjoint
    (`_forced_replay_transposed_radau`) -- so a correct three-stage adjoint
    makes the whole stack available without a Gear-2 twin.  On van der Pol
    with a white source the diffusion constant `c` and the lineshape must
    match the Gear-2 numbers (Radau is order 5, gear order 2, so they agree
    to the coarser of the two).
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def solve(method):
        cir = _vdp_with_noise(1e-6)
        m = cir.n - 1
        pss = PSS(cir, method=method, reltol=1e-12)
        x0 = np.zeros(m); x0[0] = 2.0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=6.6634, timestep=6.6634 / 240, x0=x0,
                      maxiterations=60)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        c = pac.diffusion_constant(pss)
        _Sv, L = pac.oscillator_spectrum(pss, [1e-2, 1e-1, 1.0], 0,
                                         harmonic=1)
        return c, np.asarray(L, dtype=float)

    c_g, L_g = solve('gear')
    c_r, L_r = solve('radau')
    assert abs(c_r - c_g) < 1e-2 * c_g, (c_r, c_g)
    assert np.max(np.abs(L_r - L_g)) < 0.05, (L_r, L_g)


def test_the_phase_only_spectrum_warns_above_the_amplitude_pole():
    """`oscillator_spectrum` is a LOWER BOUND above `f_amp`, and now says so.

    The method returns the PHASE contribution only. A real oscillator also has
    AMPLITUDE noise, suppressed near the carrier because the limit cycle
    restores the amplitude — but only at the amplitude-relaxation rate. Above
    that pole the amplitude noise no longer decays within a period and adds to
    the total, so this method under-reports. Relayed measurement of a
    commercial simulator's excess over the phase-only prediction: −0.54 dB at
    1 kHz and −2.90 dB at 10 kHz for `λ₂ = 0.99`.

    ⚠⚠ AND THE VALID BAND SHRINKS AS `1/Q`, which makes this §0 again rather
    than a detail. With `f_amp = -ln(λ₂)/(2πT)` and `Q = -1/ln(λ₂)`,

        f_amp = f0 / (2 pi Q)

    verified both ways at 26671.9 / 2544.2 / 253.3 Hz for
    `λ₂ = 0.90 / 0.99 / 0.999`. At `λ₂ = 0.999` the phase-only window has
    collapsed below ~253 Hz — the better the oscillator, the narrower the band
    in which this answer is the whole answer.

    ⚠ This is the OPPOSITE SIGN from the error `PSS.ppv` already warns about
    (slow nodes filtering device noise, phase noise OVER-estimated). Two
    independent mechanisms; only one was in the code before this.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0, mu, psd = 2.0 * np.pi, 0.005, 1e-6

    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=psd)

    pss = PSS(cir, method='radau', reltol=1e-12)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T0, timestep=T0 / 240, x0=np.array([2.0, 0.0]),
                  maxiterations=80)
    pac = PAC(cir, toolkit=circuit.numeric)

    f0 = 1.0 / float(pss.period)
    ## van der Pol at small `mu`: `|λ₂| = exp(-2 pi mu)`, checked by the
    ## sibling high-Q test rather than assumed here.
    lam2 = float(np.exp(-2.0 * np.pi * mu))
    f_amp = -np.log(lam2) * f0 / (2.0 * np.pi)

    ## The identity the warning is built on, both ways.
    Q = -1.0 / np.log(lam2)
    assert abs(f_amp - f0 / (2.0 * np.pi * Q)) < 1e-12 * f_amp, \
        'f_amp = f0/(2 pi Q) does not hold: %.6g vs %.6g' % (
            f_amp, f0 / (2.0 * np.pi * Q))

    def fires(offset):
        with _w.catch_warnings(record=True) as caught:
            _w.simplefilter('always')
            pac.oscillator_spectrum(pss, np.array([offset]), 0)
        return any('amplitude-relaxation' in str(c.message) for c in caught)

    assert not fires(0.1 * f_amp), \
        'warned below f_amp (%.4g Hz), where the phase-only answer is the ' \
        'whole answer' % f_amp
    assert fires(10.0 * f_amp), \
        'did NOT warn a decade above f_amp (%.4g Hz), where a commercial ' \
        'simulator measures several dB of excess' % f_amp


def test_the_orbital_spectrum_is_a_lorentzian_of_half_width_f_amp():
    """A9 step 4: `orbital_spectrum`, and the ONE gate that survived.

    Lemma 3.5 makes the orbital spectrum Lorentzians centred at
    `j*w0 + Im(mu_l)` with half-width `|Re(mu_l)| + (1/2) h^2 w0^2 c`. The
    `h = 0` term therefore has half-width exactly

        |Re(mu_2)| / (2 pi) = |ln(lam2)| f0 / (2 pi) = f_amp

    the SAME amplitude-relaxation pole `oscillator_spectrum` warns above, and
    that pole was derived independently — from a commercial simulator's
    measured excess over our phase-only answer. Two routes, one from a parity
    table and one from this paper's modal sum, landing on one quantity.

    ⚠⚠ A HEADLINE RETRACTED AND THEN RESTORED, WHICH IS THE WHOLE STORY.
    `S_orb/S_ph` reads 0.500 at `f_amp` on the default fixture. That looks
    like a law — *the orbital term reaches half the phase term at f_amp*. It
    was RETRACTED when sweeping `C` at fixed `w0` gave 0.500 / 8.002 / 0.031
    at `C` = 1 / 4 / 0.25, a 256x swing, which read as §D 0c fixture
    blindness.

    ✅ **The retraction was right on the evidence and wrong about the cause,
    and chasing the cause found a REAL DEFECT.** `floquet_modes`
    biorthonormalised on `q^T p = 1`; the variational DAE's conserved form is
    `q^T C p`, and `q` enters the covariance quadratically, so the orbital
    covariance was too large by exactly `C^2`. Against the independent
    Lyapunov reference, before: 0.0624 / 1.0001 / 16.043 — right ONLY at
    `C = 1`, the one place A9's three-way gate ever ran. After the fix:
    1.0027 / 1.0018 / 1.0036, and `w` scales as `C^-1` as predicted.

    So the law is real; the defect was hiding it. Measured at each circuit's
    OWN `f_amp` after the fix: **0.5052 / 0.5010 / 0.5005** at
    `C` = 0.25 / 1 / 4. ⚠ The first version of this test compared both
    circuits at the `C = 1` circuit's `f_amp` — but `f_amp ~ 1/C`, so it
    sampled the `C = 4` orbit at four times its own pole and read a
    difference that was its own bug.
    """
    cir, pss = _a9_vdp()
    pac = PAC(cir)
    f0 = 1.0 / float(pss.period)
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        _v, info = pss.ppv()
    lam2 = float(info['second_multiplier'])
    f_amp = -np.log(lam2) * f0 / (2.0 * np.pi)

    ## 1. THE SHAPE, tied to lam2. A Lorentzian of half-width `f_amp` is at
    ##    half its peak exactly `f_amp` away.
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        s_peak = float(pac.orbital_spectrum(
            pss, np.array([1e-6 * f_amp]), 0, H=4)[0])
        s_half = float(pac.orbital_spectrum(
            pss, np.array([f_amp]), 0, H=4)[0])
    assert abs(s_half / s_peak - 0.5) < 5e-3, \
        'the orbital line is at %.6f of its peak one f_amp out, not 0.5 — ' \
        'its half-width is not |Re(mu_2)|/(2 pi)' % (s_half / s_peak)

    ## 2. INVARIANCE UNDER THE NOISE SCALE. Both spectra are linear in the
    ##    source PSD, so their RATIO must not move; if it does, the one-sided
    ##    / two-sided conventions have drifted apart between them — which is
    ##    the defect this file has caught more than once.
    ratios = []
    for psd in (1e-6, 1e-4):
        c2, p2 = _a9_vdp(psd=psd)
        pac2 = PAC(c2)
        offs = np.array([f_amp, 100.0 * f_amp])
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            sp, _ = pac2.oscillator_spectrum(p2, offs, 0, harmonic=1)
            so = pac2.orbital_spectrum(p2, offs, 0, harmonic=1, H=4)
        ratios.append(np.asarray(so) / np.asarray(sp))
    drift = float(np.max(np.abs(ratios[0] / ratios[1] - 1.0)))
    assert drift < 5e-3, \
        'the orbital/phase ratio moved by %.3e over a 100x change in source ' \
        'PSD; both are linear in it, so a scaling convention disagrees' % drift

    ## 3. ⚠⚠ THE RATIO AT `f_amp` IS 0.5 AND IS FIXTURE-INDEPENDENT — the
    ##    assertion that would have caught the `C^2` biorthonormalisation
    ##    defect, and the one the three-way gate could not because it only
    ##    ever ran at `C = 1`.  Each circuit is evaluated at ITS OWN `f_amp`;
    ##    `f_amp ~ 1/C`, so a shared one silently samples the wrong offset.
    seen = []
    for cval, lval in ((0.25, 4.0), (1.0, 1.0), (4.0, 0.25)):
        ck, pk = _a9_vdp(cval=cval, lval=lval)
        pk_pac = PAC(ck)
        f0k = 1.0 / float(pk.period)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            _vk, ik = pk.ppv()
            fak = -np.log(float(ik['second_multiplier'])) * f0k / (2.0 * np.pi)
            spk, _ = pk_pac.oscillator_spectrum(
                pk, np.array([fak]), 0, harmonic=1)
            sok = pk_pac.orbital_spectrum(
                pk, np.array([fak]), 0, harmonic=1, H=4)
        seen.append(float(sok[0] / spk[0]))
    for cval, r in zip((0.25, 1.0, 4.0), seen):
        assert abs(r - 0.5) < 0.02, \
            'at C=%g the orbital/phase ratio at its own f_amp is %.4f, not ' \
            '0.5. A C-dependent value here is the signature of the '  \
            'biorthonormalisation defect: floquet_modes must normalise on ' \
            'q^T C p, not q^T p' % (cval, r)
    assert max(seen) / min(seen) < 1.05, \
        'the ratio at f_amp moved by %.3f across a 16x sweep in C (%r); it ' \
        'must be fixture-independent' % (max(seen) / min(seen), seen)


def test_the_orbital_spectrum_amplitude_matches_pnoise_on_a_symmetric_orbit():
    """`orbital_spectrum`'s ABSOLUTE amplitude, against a reference outside it.

    Until 2026-09-14 the amplitude was validated against NOTHING external:
    its SHAPE was tied to `lambda_2` (half-width `f_amp`) and its INTEGRAL to
    the Lyapunov covariance, but the V^2/Hz at an offset rested on the modal
    sum alone.  pnoise is the linear LPTV sideband noise -- adjoint sideband
    fold, no Floquet modes, no eq (22) -- and the three-leg chain puts it on
    one absolute scale with the externally certified Lorentzian.  Above the
    phase linewidth it is the TOTAL linear noise, so

        R = (up + lo) / (4 (S_ph + S_orb))

    must be 1.  Measured (H = 8, 16 sidebands), `R` at 0.1 / 1 / 3 / 10 f_amp:

        van der Pol Q=8, C=1       0.9997  0.9993  0.9981  0.9855 (0.2 f0)
        C=4, L=1/4 (same w0)       0.9997  0.9996  0.9995  0.9986
        Q=50                       0.9997  0.9996  0.9995  0.9992

    and the phase term alone reads 1.50 at f_amp and 1.99 above it, so the
    orbital term is carrying HALF the noise there and is right to 0.1 %.
    ⚠ The C = 4 row is the one that matters: the orbital covariance was once
    too large by exactly `C^2`, and a unit-reactance fixture cannot see that.
    The drift toward f0 (0.9855 at 0.2 f0) is shared by pnoise's PM and AM
    parts equally, i.e. the Lorentzian approximation, not the orbital term.

    ⚠⚠ SYMMETRIC ORBITS ONLY.  On `_hostile_oscillator` the sum over-states
    pnoise by ~3.3x above f_amp, and a Monte Carlo sides with pnoise.  The
    cross term in Traversa & Bonani's DC-harmonic form is ~1e-8 there; with
    EVERY harmonic kept it closes the gap (2026-09-15) -- see
    `test_the_modal_spectrum_with_the_full_correlation_closes_on_pnoise`.
    See `test_the_orbital_spectrum_sum_over_states_on_an_asymmetric_orbit_and_says_so`.
    """
    import warnings as _w
    cir, pss = _a9_vdp(cval=4.0, lval=0.25)
    pac = PAC(cir, toolkit=circuit.numeric)
    f0 = 1.0 / float(pss.period)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        _v, info = pss.ppv()
        f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
        offs = np.array([0.3, 3.0, 10.0]) * f_amp
        Sph = np.asarray(pac.oscillator_spectrum(pss, offs, 0)[0])
        Sorb = pac.orbital_spectrum(pss, offs, 0, harmonic=1, H=8)
        for f, sph, sorb in zip(offs, np.asarray(Sph), np.asarray(Sorb)):
            up, _ = pac.pnoise(pss, f0 + f, 0, maxsidebands=16)
            lo, _ = pac.pnoise(pss, f0 - f, 0, maxsidebands=16)
            ## (both one-sided: the mean of the two sidebands)
            tot = float(np.real(up) + np.real(lo)) / 2.0
            ## the orbital term must be carrying real weight, or `R` is a
            ## statement about the phase term alone
            if f > f_amp:
                assert sorb / sph > 0.5, (f / f_amp, sorb / sph)
            assert abs(tot / (sph + sorb) - 1.0) < 5e-3, \
                'at %.1f f_amp pnoise/(S_ph+S_orb) = %.4f; the orbital ' \
                'spectrum amplitude is off (phase alone would read %.4f)' \
                % (f / f_amp, tot / (sph + sorb), tot / sph)


def test_the_line_shape_spectra_refuse_a_harmonic_that_has_no_line():
    """No carrier, no line: refuse instead of returning a plausible zero.

    `oscillator_spectrum` and `orbital_spectrum` are LINE-SHAPE models -- they
    broaden the orbit's own harmonic lines.  A Monte Carlo of the SDE (64 x 4000
    periods), which agreed with pnoise to 1 % at every harmonic, measured them
    AWAY from the fundamental on van der Pol C=4 Q=8 (2026-09-14):

        symmetric orbit    k=0     k=1      k=2     k=3
        MC / model         0.0074  1.0015   3.21    3.70

    At an even harmonic of a symmetric orbit the output has no component and
    the modes have no Fourier content there, so `oscillator_spectrum` returned
    ~0 and `orbital_spectrum` returned other lines' tails, against a true
    5.6e-6 V^2/Hz.  Both now REFUSE that case (as `am_pm` refuses a harmonic
    with no carrier), and `harmonic=0`, which the Lorentzian returned as zeros
    by construction.  ⚠ The asymmetric orbit's k = 2 has a real line and is
    NOT refused -- its mismatch (0.40) is recorded, not guarded.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    offs = None
    for a_, expect_refusal in ((0.0, True), (0.30, False)):
        cir, pss = _a9_vdp(cval=4.0, lval=0.25, a=a_)
        pac = PAC(cir, toolkit=circuit.numeric)
        f0 = 1.0 / float(pss.period)
        offs = np.array([0.01 * f0])
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            ## the fundamental always works
            assert np.all(np.isfinite(np.asarray(
                pac.oscillator_spectrum(pss, offs, 0, frequency_aware=False)[0])))
            for call in (
                    lambda: pac.oscillator_spectrum(pss, offs, 0, harmonic=2,
                                                    frequency_aware=False),
                    lambda: pac.orbital_spectrum(pss, offs, 0, harmonic=2, H=8)):
                if expect_refusal:
                    with pytest.raises(ValueError, match='pnoise'):
                        call()
                else:
                    got = call()
                    got = got[0] if isinstance(got, tuple) else got
                    assert np.all(np.isfinite(np.asarray(got)))
            if expect_refusal:
                with pytest.raises(ValueError, match='harmonic 0'):
                    pac.oscillator_spectrum(pss, offs, 0, harmonic=0,
                                            frequency_aware=False)


def test_oscillator_spectrum_is_frequency_aware_above_the_slow_corner():
    """The phase Lorentzian with `c(f)` instead of `c` -- built 2026-09-14.

    The DC-PPV Lorentzian assumes a noise current moves the phase instantly.
    Where the response goes through a slow mode it is filtered above that
    mode's corner, and the closed form over-states.  `frequency_aware=True`
    (the default) uses `c(f)` from `PSS.frequency_aware_ppv`, integrated by
    `diffusion_constant`'s own quadrature.  Gated against pnoise's PM content
    (Monte-Carlo-confirmed on the first fixture), measured:

        van der Pol C=4 Q=8 a=0.30    1 f_amp   10 f_amp
          pm / S_v, closed form        0.647     0.302
          pm / S_v, frequency-aware    0.998     0.980
        A2 slow node, tau/T = 100     1e-3 f0   1e-2 f0
          closed form                  0.729     0.027
          frequency-aware              1.004     1.004
        symmetric control (a=0): unchanged to 1e-3

    ⚠ The prototype first read `c(0)/c = 0.9963` on the asymmetric fixtures and
    1.0000 on the symmetric one with IDENTICAL samples -- the quadrature used
    `frequency_aware_ppv`'s `times`, which is one entry short and drops the last
    step.  Over the orbit's full grid it is 1.000000000000; asserted below.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric

    def pm_ratio(pac, pss, off, ov, fa, sb=16):
        S = float(np.asarray(pac.oscillator_spectrum(
            pss, np.array([off]), ov, frequency_aware=fa)[0])[0])
        _am, pm, _ = pac.am_pm_noise(pss, off, ov, carrier=1, maxsidebands=sb)
        return pm / S

    ## 1. AM-to-PM coupling above f_amp
    cir, pss = _a9_vdp(cval=4.0, lval=0.25, a=0.30)
    pac = PAC(cir, toolkit=circuit.numeric)
    f0 = 1.0 / float(pss.period)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        _v, info = pss.ppv()
        f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
        c = pac.diffusion_constant(pss)
        assert abs(pac.frequency_aware_diffusion(pss, 1e-9 * f0) / c - 1.0) < 1e-6
        for mult, closed in ((1.0, 0.647), (10.0, 0.302)):
            r_fa = pm_ratio(pac, pss, mult * f_amp, 0, True)
            r_cf = pm_ratio(pac, pss, mult * f_amp, 0, False)
            assert abs(r_fa - 1.0) < 0.03, \
                'at %g f_amp the frequency-aware phase spectrum is off pnoise ' \
                'PM by %.4f' % (mult, r_fa)
            assert abs(r_cf - closed) < 0.03, \
                'the closed form must still over-state here (%.4f, measured ' \
                '%.3f) or this gate no longer separates the two' % (r_cf, closed)

    ## 2. symmetric control: nothing to correct
    cir_s, pss_s = _a9_vdp(cval=4.0, lval=0.25, a=0.0)
    pac_s = PAC(cir_s, toolkit=circuit.numeric)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        f0s = 1.0 / float(pss_s.period)
        _v, info_s = pss_s.ppv()
        fa_s = -np.log(float(info_s['second_multiplier'])) * f0s / (2 * np.pi)
        o = np.array([10.0 * fa_s])
        s_fa = float(np.asarray(pac_s.oscillator_spectrum(pss_s, o, 0)[0])[0])
        s_cf = float(np.asarray(pac_s.oscillator_spectrum(
            pss_s, o, 0, frequency_aware=False)[0])[0])
    assert abs(s_fa / s_cf - 1.0) < 2e-3, s_fa / s_cf

    ## 3. a source behind a slow node (A2)
    T0 = 6.6634
    c2 = SubCircuit()
    c2.add_node('v'); c2.add_node('w'); c2.add_node('x')
    c2['C'] = C('v', gnd, c=1.0)
    c2['L'] = L('v', 'x', L=1.0); c2['Rl'] = R('x', gnd, r=0.2)
    c2['B'] = BSource('v', gnd, gnd, 'v',
                      i_func=lambda u: (u - u ** 3 / 3.0) + 0.25 * (u ** 2 - 2.0))
    c2['Rs'] = R('v', 'w', r=1e2); c2['Cs'] = C('w', gnd, c=100.0 * T0 / 1e2)
    c2['n'] = IS('w', gnd, i=0.0, noisePSD=1e-6)
    p2 = PSS(c2, method='gear', reltol=1e-11)
    x0 = np.zeros(c2.n - 1); x0[0] = 2.0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p2.solve(period=T0, timestep=T0 / 240, x0=x0, maxiterations=200)
        assert p2.converged
        pac2 = PAC(c2, toolkit=circuit.numeric)
        ov = [str(n_) for n_ in c2.nodes].index('v')
        f02 = 1.0 / float(p2.period)
        for r_, closed in ((1e-3, 0.729), (1e-2, 0.027)):
            r_fa = pm_ratio(pac2, p2, r_ * f02, ov, True, sb=32)
            r_cf = pm_ratio(pac2, p2, r_ * f02, ov, False, sb=32)
            assert abs(r_fa - 1.0) < 0.02, ('A2 frequency-aware', r_, r_fa)
            assert abs(r_cf / closed - 1.0) < 0.05, ('A2 closed form', r_, r_cf)


def test_the_orbital_spectrum_sum_over_states_on_an_asymmetric_orbit_and_says_so():
    """The other half of the amplitude check: where `S_ph + S_orb` is WRONG.

    Same fixture family as the symmetric test (van der Pol, C = 4, L = 1/4,
    Q = 8) with an `a u^2` asymmetry, against pnoise at 10 f_amp:

        a      half-wave asymmetry   R = (up+lo)/4(S_ph+S_orb)
        0.00   0.000                 0.9986
        0.05   0.017                 0.9969
        0.10   0.033                 0.9717
        0.20   0.067                 0.6914
        0.30   0.100                 0.3090

    Smooth, monotone, 1 at `a = 0` -- the prediction named before the sweep.

    ⚠⚠ WHICH SIDE IS RIGHT WAS SETTLED BY MONTE CARLO, AND IT OVERTURNED THE
    FIRST ATTRIBUTION.  A vectorised trapezoidal SDE (64 x 4000 periods,
    `Var(i) = PSD/(2h)`), one-sided PSD in 8-12 f_amp on both sidebands:
    a = 0 control MC/pnoise 1.009, MC/modal 1.008; a = 0.30 MC/pnoise
    **1.011**, MC/modal **0.313** -- and MC reproduces pnoise's sideband
    ASYMMETRY (6.07e-4 / 7.64e-4 against 6.01e-4 / 7.55e-4), which the
    modal sum does not have.  `S_corr` in eq (92)'s DC-harmonic form is ~1e-8
    of the total on this fixture (it needs the PPV's DC at the noise source's
    row, which the tank inductor shorts) -- ⚠⚠ but the phase-orbital
    correlation with EVERY harmonic kept is -1.1 to -2.4x the orbital term
    and closes the gap (2026-09-15, `PAC.modal_spectrum`): the "cross term is
    not the cause" reading held only for the truncated form.

    ⚠ The warning that fired here before blamed an O(h) grid residual of the
    adjoint replay and said to refine -- true of `orbital_correlation`'s
    covariance, and wrong for this sum, which over-states by 1.45x at this
    asymmetry whatever the grid.  So this test pins BOTH the gap (a presence
    claim: it must stay large) and that the warning names it.
    """
    import warnings as _w
    cir, pss = _a9_vdp(cval=4.0, lval=0.25, a=0.20)
    pac = PAC(cir, toolkit=circuit.numeric)
    f0 = 1.0 / float(pss.period)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        _v, info = pss.ppv()
    f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
    off = np.array([10.0 * f_amp])
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter('always')
        Sorb = float(pac.orbital_spectrum(pss, off, 0, harmonic=1, H=8)[0])
    assert any('over-states' in str(w.message)
               and 'pnoise' in str(w.message) for w in caught), \
        'orbital_spectrum must warn that the sum over-states on an asymmetric ' \
        'orbit; got %r' % [str(w.message) for w in caught]
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        ## the closed form: this test's claim is about the modal SUM as
        ## `orbital_spectrum`'s docstring states it (DC-PPV phase term)
        Sph = float(np.asarray(pac.oscillator_spectrum(
            pss, off, 0, frequency_aware=False)[0])[0])
        up, _ = pac.pnoise(pss, f0 + off[0], 0, maxsidebands=16)
        lo, _ = pac.pnoise(pss, f0 - off[0], 0, maxsidebands=16)
    R = float(np.real(up) + np.real(lo)) / (2.0 * (Sph + Sorb))
    assert 0.6 < R < 0.8, \
        'at a = 0.20 pnoise/(S_ph+S_orb) = %.4f at 10 f_amp (measured 0.6914). ' \
        'Near 1 means the cross term became negligible or the sum changed; ' \
        'far below means something else moved' % R


def test_the_three_leg_chain_puts_pnoise_the_am_pm_split_and_the_lorentzian_on_one_absolute_scale():
    """2026-09-08 (Andreas: "run the three-leg experiment").  The chain

        reference simulator --(S_v = 0.5000 x one-sided PSD, four decades,
        2026-09-05)--> lorentzian --(Rizzoli overlay)--> pnoise --(AM/PM
        identity)--> S_am, S_pm

    has an EXTERNAL anchor at one end and no two links share machinery, so
    it is the one structure that can localise a common scale factor -- the
    signature of the "oscillator sideband rows ~1e-12, cause unestablished"
    caveat.  Measured on the single-cluster van der Pol at Q = 100 (mu =
    1/(2 pi Q), white 1e-6 A^2/Hz current noise) and on the LC oscillator
    the modulation stack was certified on:

        f/f0     (up+lo)/(2 S_v)   Q=100    LC mu=1
        1e-4          1.003                 0.9993
        1e-3          1.284                 0.9993
        1e-2          1.974                 0.9992
        1e-1          1.993                 0.9940
        identity residual: 1e-16 (Q=100), 2.5e-12 at 64 sidebands and
        2.8e-9 at 8 (LC) -- truncation, converging away as in the driven case.

    The ratio is 1 below the AM corner f0/(2 pi Q_lambda) = 1.6e-3 f0 (⚠ was
    written f0/(4 pi Q) = 8e-4 until 2026-09-09 -- the measured transition
    1e-3 .. 3e-3 f0 already said 2 pi; pinned by the corner test below) (the amplitude
    mode restores AM, only PM survives, and PM IS the Lorentzian) and 2
    above it where an LTI tank splits additive noise equally between AM and
    PM -- the prediction named before running, to the corner's decade.  So
    pnoise, the split and the externally certified closed form sit on ONE
    absolute scale and the "~1e-12" symptom is not in the noise split.
    Pinned: identity < 1e-9; ratio within 2 % of 1 at 1e-4 f0 and within
    3 % of 2 at 1e-2 and 1e-1 f0; and the LC fixture the closed form was
    certified on, overlay 1 within 1 % over three decades (below); and S_pm alone within 1 % of S_v at
    every offset -- S_pm is the PM content of ONE sideband (the pair's until
    2026-09-28, when `S_v` became one-sided and `S_pm` per sideband) (the PM part is the Lorentzian everywhere, the AM part is what
    the ratio adds).  ⚠ First written as 2 S_v and failed at 0.997 off:
    the pair, not one sideband.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir, mu = _a10_vdp(100.0)
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=2 * np.pi, timestep=2 * np.pi / 240,
                  x0=np.array([2.0, 0.0]), maxiterations=100)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    f0 = 1.0 / float(pss.period)
    ov = [str(n) for n in cir.nodes].index('v')
    offs = f0 * np.array([1e-4, 1e-2, 1e-1])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        Sv, _ = pac.oscillator_spectrum(pss, offs, ov)
        rows = []
        for f in offs:
            up, _ = pac.pnoise(pss, f0 + f, ov, maxsidebands=16)
            lo, _ = pac.pnoise(pss, f0 - f, ov, maxsidebands=16)
            am, pm, _ = pac.am_pm_noise(pss, f, ov, carrier=1, maxsidebands=16)
            rows.append((float(np.real(up)), float(np.real(lo)),
                         float(np.real(am)), float(np.real(pm))))
    for (up, lo, am, pm), sv, f, want in zip(rows, Sv, offs, (1.0, 2.0, 2.0)):
        assert abs(2.0 * (am + pm) - (up + lo)) / (up + lo) < 1e-9, (f / f0, am + pm, up + lo)
        ratio = (up + lo) / (2.0 * sv)
        assert abs(ratio / want - 1.0) < 0.03, \
            'overlay (up+lo)/(2 S_v) = %.4f at %.0e f0, expected %.0f' % (ratio, f / f0, want)
        assert abs(pm / sv - 1.0) < 0.01, (f / f0, pm, sv)

    ## AND ON THE FIXTURE THE MODULATION STACK WAS CERTIFIED ON (Andreas,
    ## 2026-09-08 evening): the LC oscillator at mu = 1 -- the external
    ## anchor's own circuit, so the chain's first link is tied to it directly.
    ## Q = 0.5 puts the AM corner (~0.5 f0) outside the band, so the overlay
    ## reads 1 throughout: measured 0.9993 / 0.9993 / 0.9992 at 1e-4 / 1e-3 /
    ## 1e-2 f0, identity 2.5e-12 at 64 sidebands (2.8e-9 at 8).
    cir, pss, pac = _lc_osc(psd=1e-6, npts=240)
    f0 = 1.0 / float(pss.period)
    ov = [str(n) for n in cir.nodes].index('v')
    offs = f0 * np.array([1e-4, 1e-3, 1e-2])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        Sv, _ = pac.oscillator_spectrum(pss, offs, ov)
        for f, sv in zip(offs, Sv):
            up, _ = pac.pnoise(pss, f0 + f, ov, maxsidebands=32)
            lo, _ = pac.pnoise(pss, f0 - f, ov, maxsidebands=32)
            am, pm, _ = pac.am_pm_noise(pss, f, ov, carrier=1, maxsidebands=32)
            up, lo, am, pm = (float(np.real(x)) for x in (up, lo, am, pm))
            assert abs(2.0 * (am + pm) - (up + lo)) / (up + lo) < 1e-8, (f / f0, am + pm, up + lo)
            assert abs((up + lo) / (2.0 * sv) - 1.0) < 0.01, ('LC overlay', f / f0, (up + lo) / (2.0 * sv))
            assert abs(pm / sv - 1.0) < 0.01, ('LC S_pm', f / f0, pm, sv)


def test_a_source_behind_a_slow_node_rolls_off_the_lorentzian_as_the_ppv_harmonics_say():
    """A2 RESOLVED (Andreas, 2026-09-08).  A white current source behind a
    slow RC node reaches the tank through `1/(1 + j w tau)`.  The DC PPV
    uses the DC transfer, so the Lorentzian from `c` holds only below the
    source's corner `T/(2 pi tau)`; above it the true skirt (`pnoise`, the
    conversion computation) is the DC-PPV one scaled by the PPV-harmonic-
    weighted filter, with NO free constant:

        ratio(r) = (|c_0|^2 F_0(f) + 2 sum_{k>=1} |c_k|^2) / sum_k |c_k|^2 (two-sided),
        F_0(f) = 1 / (1 + (2 pi f tau)^2),

    `c_k` the Fourier coefficients of the PPV entry at the source node,
    which ALREADY carry the path's transfer at k f0 (so the filter enters
    only as F_0(f)/F_0(0) on k = 0; the k >= 1 terms are the floor, and
    |c_1|^2/|c_0|^2 at w carries the (T/tau)^2 of the path at f0).
    Only k = 0 varies over the band, so the effect needs `c_0 != 0`, which
    needs BOTH an asymmetric core AND tank loss (an ideal tank inductor
    shorts DC; a half-wave-symmetric orbit has no DC PPV) -- the coloured-
    upconversion gate's pattern, seen from the transfer side.  Measured
    (a = 0.25, loss 0.2, tau/T = 100, Rs = 1e2): |c0|/|c1| at w = 59.2,
    floor 5.9e-4; ratio 0.9964 / 0.7291 / 0.2124 / 0.0268 / 0.00335 at
    r = 1e-4 / 1e-3 / 3.2e-3 / 1e-2 / 3.2e-2 against 0.9963 / 0.7298 /
    0.2125 / 0.02685 / 0.00328 predicted (within 2 %); at 0.1 f0 measured
    0.00096 vs 0.00086 -- on S_pm with the core control subtracted a
    residual of 7 % that grows with the asymmetry (0.9 / 7.1 / 11.9 % at
    a = 0.1 / 0.25 / 0.4), confined to r >= 5e-2 where the k = 0 term has
    fallen to the k >= 1 floor; candidate AM-to-PM conversion through the
    asymmetric core, not modelled; the DC-PPV Lorentzian over-states the slow-node
    source by 1000x there (Lai's sign).  At tau/T = 10 the same sum holds
    to 2 % up to 3.2e-2 f0 with the corner and floor shifted 10x and 100x.
    Controls: the source at the core is flat at 0.999; the odd core
    with loss has G_0 = 2.6e-10 and is flat.  `c` itself is unaffected (the
    filter removes only high-frequency content), which is why a Monte Carlo
    of `c` -- the 2026-09-03 gate -- read a null at tau/T = 10.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    T0 = 6.6634
    tau_over_T = 100.0

    def build(src, asym, loss):
        c = SubCircuit()
        c.add_node('v'); c.add_node('w'); c.add_node('x')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', 'x', L=1.0); c['Rl'] = R('x', gnd, r=loss)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: (u - u ** 3 / 3.0) + asym * (u ** 2 - 2.0))
        c['Rs'] = R('v', 'w', r=1e2); c['Cs'] = C('w', gnd, c=tau_over_T * T0 / 1e2)
        c['n'] = IS(src, gnd, i=0.0, noisePSD=1e-6)
        return c

    def measure(src, asym, loss, rs):
        cir = build(src, asym, loss)
        pss = PSS(cir, method='gear', reltol=1e-11)
        x0 = np.zeros(cir.n - 1); x0[0] = 2.0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T0, timestep=T0 / 240, x0=x0, maxiterations=200)
            assert pss.converged
            pac = PAC(cir, toolkit=circuit.numeric)
            f0 = 1.0 / float(pss.period)
            names = [str(n) for n in cir.nodes]; ov = names.index('v'); isrc = names.index(src)
            _v1, info = pss.ppv()
            V = np.array([np.asarray(s_, float) for s_ in info['samples']])
            if abs(V[0] - V[-1]).max() < 1e-9 * abs(V).max():
                V = V[:-1]
            ## two-sided POWER per harmonic, |c_-k|^2 + |c_k|^2.  ⚠ The PPV
            ## entry at the source node ALREADY carries the path's transfer
            ## at each harmonic, so the RC filter enters ONLY as the ratio
            ## F_0(f)/F_0(0) on the k = 0 term; the first version applied
            ## F_k(k f0) to the k >= 1 terms too (filtered twice) and doubled
            ## already-doubled rfft amplitudes, putting the floor at 6e-9
            ## instead of 5.9e-4 -- and the top of the sweep then read as
            ## "1e5 above the floor" (peer's arithmetic caught the wording).
            G = np.fft.rfft(V[:, isrc]) / V.shape[0]
            G2 = np.abs(G) ** 2; G2[1:] *= 2.0
            tau = tau_over_T * T0
            offs = f0 * np.array(rs)
            ## the CLOSED FORM, deliberately: this test pins how the DC-PPV
            ## Lorentzian departs from pnoise; the default now corrects it
            ## (`test_oscillator_spectrum_is_frequency_aware_above_the_slow_corner`)
            Sv, _ = pac.oscillator_spectrum(pss, offs, ov, frequency_aware=False)
            out = []
            for f, sv in zip(offs, Sv):
                up, _ = pac.pnoise(pss, f0 + f, ov, maxsidebands=32)
                lo, _ = pac.pnoise(pss, f0 - f, ov, maxsidebands=32)
                F0 = 1.0 / (1.0 + (2 * np.pi * f * tau) ** 2)
                pred = (G2[0] * F0 + float(np.sum(G2[1:]))) / float(np.sum(G2)) if src == 'w' else 1.0
                out.append(((float(np.real(up)) + float(np.real(lo))) / (2.0 * sv), pred))
        return abs(G[0]) / abs(G[1]), out

    rs = (1e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2)
    g01, rows = measure('w', 0.25, 0.2, rs)
    assert g01 > 1.0, 'the slow node\'s PPV entry should be DC-dominated with asymmetry AND loss; |c0|/|c1| = %.3e' % g01
    for r, (ratio, pred) in zip(rs, rows):
        assert abs(ratio / pred - 1.0) < 0.03, (r, ratio, pred)
    assert rows[-1][0] < 0.05, 'the Lorentzian should over-state the slow-node source >= 20x at 1e-2 f0; ratio %.4f' % rows[-1][0]
    g01_v, rows_v = measure('v', 0.25, 0.2, rs)
    for r, (ratio, _p) in zip(rs, rows_v):
        assert abs(ratio - 1.0) < 0.01, ('core-injection control', r, ratio)
    g01_s, rows_s = measure('w', 0.0, 0.2, rs)
    assert g01_s < 1e-6, g01_s
    for r, (ratio, _p) in zip(rs, rows_s):
        assert abs(ratio - 1.0) < 0.01, ('odd-core control', r, ratio)


def test_the_frequency_aware_ppv_is_the_ppv_at_dc_and_corners_at_the_slow_multiplier():
    """2026-09-08 night (Andreas: "do the frequency-aware PPV first").
    `PSS.frequency_aware_ppv(offset)` is `ppv()`'s bordered system at
    `alpha = exp(-j w_s T)` (Lai 2008 eq. 24, verified: its w_s = 0 point
    "is the augmented PPV extraction equation").  Three gates, named
    before running and scored:
      1. at offset 0 it IS `ppv()`: anchor and pair samples to 1e-12
         (measured 3e-16) -- after normalising by `v . xdot = 1` as `ppv`
         does; the border alone left a factor -0.94;
      2. the SLOW multiplier's own coefficient (dense eigenbasis, the
         docs session's construction) rises ten-fold per decade below the
         corner `(1 - mu_2)/(2 pi T)` and plateaus above; corner between
         1e-3 and 3e-3 at tau/T = 100 (predicted 1.58e-3) and ten times
         higher at tau/T = 10, plateau 2.45e-6 vs their 2.29e-6 and
         scaling as T/tau (2.5e-5 at tau/T = 10).  ⚠ The NORM admixture
         is dominated by the core's fast amplitude mode (linear in r, no
         corner in band, tau-independent) -- the wrong metric, as they
         warned; and the first version of this gate printed the THIRD
         multiplier's coefficient by an index slip;
      3. application (A2): the ratio of its harmonic powers at w_s to the
         DC ones, with NO explicit model of the RC path, reproduces
         pnoise's S_pm/(4 S_v) for a source behind the slow node within
         0.6 % to r = 1e-2 and 2 % at 5e-2 on the lossy asymmetric fixture
         (the DC harmonic sum needed the path's filter by hand and was
         4 % off at 5e-2).  With the PAIR-space envelope the two differed
         by -6 to -7 % at r = 0.1, unchanged at 480 points, which I read
         as "not the envelope"; the SECOND-order state-space envelope
         (`samples`, ppv's own propagation on the complex anchor) moves
         that to -2 to -3 % and 5e-2 to within 0.7 % -- the pair-space
         object differs by a consistency correction, not an h-error, so
         the grid test was blind.  What remains at 0.1 f0 is recorded, not
         resolved.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    T0 = 6.6634

    def build(tau_over_T, rs, asym, loss, src=None):
        c = SubCircuit(); c.add_node('v'); c.add_node('w')
        c['C'] = C('v', gnd, c=1.0)
        if loss > 0:
            c.add_node('x'); c['L'] = L('v', 'x', L=1.0); c['Rl'] = R('x', gnd, r=loss)
        else:
            c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: (u - u ** 3 / 3.0) + asym * (u ** 2 - 2.0))
        c['Rs'] = R('v', 'w', r=rs); c['Cs'] = C('w', gnd, c=tau_over_T * T0 / rs)
        if src:
            c['n'] = IS(src, gnd, i=0.0, noisePSD=1e-6)
        return c

    def solve(c):
        pss = PSS(c, method='gear', reltol=1e-11)
        x0 = np.zeros(c.n - 1); x0[0] = 2.0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T0, timestep=T0 / 240, x0=x0, maxiterations=200)
        assert pss.converged
        return pss

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ## 1. identity at DC, and 2. the corner, on the odd lossless core
        for tau_over_T, lo, hi in ((100.0, 1e-3, 3e-3), (10.0, 1e-2, 3e-2)):
            pss = solve(build(tau_over_T, 1e2, 0.0, 0.0))
            f0 = 1.0 / float(pss.period)
            v0, i0 = pss.ppv()
            vf, inf = pss.frequency_aware_ppv(0.0)
            assert np.linalg.norm(vf - v0) < 1e-12 * np.linalg.norm(v0)
            assert np.linalg.norm(inf['samples_pair'] - i0['samples_pair']) < 1e-12 * np.linalg.norm(i0['samples_pair'])
            ## the SECOND-order state-space envelope (ppv's own propagation
            ## on the complex anchor, lifted into _ppv_propagate 2026-09-08)
            assert np.linalg.norm(inf['samples'] - i0['samples']) < 1e-12 * np.linalg.norm(i0['samples'])
            rc = inf['corner'] / f0
            assert lo < rc < hi, (tau_over_T, rc)
            c = {r: pss.frequency_aware_ppv(r * f0)[1]['mode_content'][0] for r in (1e-5, 1e-4, 1e-1)}
            assert abs(c[1e-4] / c[1e-5] / 10.0 - 1.0) < 0.02, c        # ten-fold per decade below the corner
            ## plateau: far below the linear extrapolation (1000x); measured 16x
            ## at tau/T = 100 and 157x at tau/T = 10 (the corner ten times higher)
            assert c[1e-1] / c[1e-4] < 300.0, c
            plateau = c[1e-1]
            if tau_over_T == 100.0:
                ## ⚠⚠ THE ABSOLUTE PLATEAU IS PHASE-SPECIFIC -- WIDE ON PURPOSE.
                ## `mode_content` is read from the PPV AT t = 0, so it depends on
                ## which phase of the orbit the shooting solve lands on, and the
                ## solve's phase is set by its seed.  Three independent samplings
                ## of the SAME orbit (period and |mu_2| identical to five digits,
                ## orbit min/max to four):
                ##   six seeds, frozen rule : 2.45 / 2.32 / 1.94 / 1.83 / 2.25 / 2.44 e-6
                ##   this file's seed sweep : 2.256 .. 2.641 e-6 (seeds 1.5 .. 2.5)
                ##   the B3 pair            : 2.2461 / 2.4445 e-6
                ## and `phase_rule='reselect'`, which lands on a different phase of
                ## that same orbit, reads 1.6437e-6.  The old pin was `2.45e-6
                ## +-10 %` = [2.205, 2.695]e-6, which HALF of the recorded seeds
                ## fail; it passed only because this test's seed (x0[0] = 2.0) is
                ## fixed.
                ## ⚠ AN INTERVAL, NOT centre+-percent: a first rewrite of this line
                ## used `2.2e-6 +-25 %` = [1.65, 2.695]e-6, whose LOWER EDGE sits
                ## 0.08 % above the known-legitimate reselect reading 1.6437e-6 --
                ## a knife edge that would flip on rounding, i.e. the same defect
                ## in a new costume.  [1.3e-6, 3.0e-6] clears every phase on record
                ## by >= 1.13x at BOTH ends and still rejects a wrong POWER of
                ## T/tau: the tau/T = 10 plateau is 2.3e-5, 7.7x outside it.
                assert 1.3e-6 < plateau < 3.0e-6, plateau
                p100 = plateau
            else:
                ## ⚠ THE RATIO IS THE PHASE-INVARIANT QUANTITY, and it is what this
                ## gate is really about: across the seed sweep above the absolute
                ## moved 17 % while `plateau_10 / plateau_100` read 10.349, 10.350,
                ## 10.351, 10.352, 10.353 -- a 0.0 % spread, because both plateaus
                ## are read at the same phase and the phase factor cancels.
                ## ⚠ Tolerance stays 10 %, NOT tightened: the measured ratio is
                ## 10.35, not 10.0, and whether that 3.5 % is physical or a grid
                ## effect is NOT established -- tightening onto 10.35 would pin an
                ## unexplained number at 240 points.
                assert abs(plateau / p100 / 10.0 - 1.0) < 0.1, (plateau, p100)   # scales as T/tau
        ## 3. the application: a source behind the slow node
        cir = build(100.0, 1e2, 0.4, 0.2, src='w')
        pss = solve(cir); f0 = 1.0 / float(pss.period)
        pac = PAC(cir, toolkit=circuit.numeric)
        names = [str(n) for n in cir.nodes]; ov = names.index('v'); iw = names.index('w')
        S0 = pss.frequency_aware_ppv(0.0)[1]['samples'][:, iw]
        P0 = float(np.sum(np.abs(np.fft.fft(S0) / S0.shape[0]) ** 2))
        for r, tol in ((1e-3, 0.01), (1e-2, 0.01), (5e-2, 0.02), (1e-1, 0.05)):
            f = r * f0
            ## the CLOSED FORM, scaled by the frequency-aware PPV's own ratio
            ## `Pf/P0` -- which is what the default now does internally
            Sv, _ = pac.oscillator_spectrum(pss, np.array([f]), ov,
                                            frequency_aware=False)
            _am, pm, _ = pac.am_pm_noise(pss, f, ov, carrier=1, maxsidebands=32)
            Sf = pss.frequency_aware_ppv(f)[1]['samples'][:, iw]
            Pf = float(np.sum(np.abs(np.fft.fft(Sf) / Sf.shape[0]) ** 2))
            assert abs((float(np.real(pm)) / float(Sv[0])) / (Pf / P0) - 1.0) < tol, (r, pm, Pf / P0)


def test_band_spread_tells_a_band_mean_from_a_point_value():
    """`PAC.band_spread` -- the guard for the convention that cost this
    review a day.

    Far above the AM corner `S(r)·r²` is flat, so a band mean IS a point
    value and the distinction is invisible; that is the case every earlier
    gate here was written on.  It is not general.  A source behind a slow RC
    node has an in-band spectrum that is NOT `1/r²` (its `k = 0` term is
    filtered at the RC corner while the `k >= 1` terms are not), and the two
    conventions then differ by percents -- which is larger than most of the
    agreements this file asserts, so a comparison that takes a band mean on
    one side and a point value on the other measures the convention.

    Asserted on ONE orbit with the noise source MOVED, so nothing but the
    source location differs: in the tank the spread is ~1.01 and the band
    mean matches the midpoint to 0.02 %; behind the slow node the spread is
    ~1.35 and they differ by ~4 %.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0 = 6.6634

    def fixture(src):
        c = SubCircuit()
        for n in ('v', 'w', 'x'):
            c.add_node(n)
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', 'x', L=1.0)
        c['Rl'] = R('x', gnd, r=0.2)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: (u - u ** 3 / 3.0) + 0.25 * (u ** 2 - 2.0))
        c['Rs'] = R('v', 'w', r=1e2)
        c['Cs'] = C('w', gnd, c=100.0 * T0 / 1e2)
        c['n'] = IS(src, gnd, i=0.0, noisePSD=1e-6)
        return c

    out = {}
    for src in ('v', 'w'):
        cir = fixture(src)
        pss = PSS(cir, method='gear', reltol=1e-11)
        x0 = np.zeros(cir.n - 1)
        x0[0] = 2.0
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T0, timestep=T0 / 240, x0=x0, maxiterations=200)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        ov = [str(n) for n in cir.nodes].index('v')
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            spread, info = pac.band_spread(pss, ov, (0.08, 0.15), points=5,
                                           quantity='S_pm', maxsidebands=32)
        out[src] = (spread, info['mean_over_point'])
    ## the tank source: flat, so the two conventions agree and the shortcut
    ## is licensed
    assert out['v'][0] < 1.05, out['v']
    assert abs(out['v'][1] - 1.0) < 5e-3, out['v']
    ## behind the slow node: NOT flat, and the conventions differ by percents
    assert out['w'][0] > 1.2, out['w']
    assert abs(out['w'][1] - 1.0) > 0.02, out['w']
    ## and the guard must SEPARATE them, which is the whole point
    assert out['w'][0] > 1.15 * out['v'][0], out


def test_the_modal_spectrum_with_the_full_correlation_closes_on_pnoise():
    """⚠⚠ E6 CLOSED BY THE PHASE-ORBITAL CORRELATION WITH EVERY HARMONIC KEPT.

    `oscillator_spectrum(frequency_aware=False) + orbital_spectrum` over-states
    an asymmetric orbit's total by ~3x (Monte-Carlo-confirmed), and Traversa &
    Bonani's DC-harmonic cross term is ~1e-8 there.  `PAC.modal_spectrum`
    builds phase, orbital and correlation from ONE modal transfer
    (every Floquet mode, every harmonic, every input sideband), and its total
    must close on pnoise.  Measured (400 points, H = 8, 16 sidebands),
    `total/(pnoise/2)` at +3 / +10 / -10 f_amp:

        a = 0.00   1.0006  1.0006  1.0006     correlation/orbital +0.013 +0.047 -0.052
        a = 0.30   1.0145  1.0057  1.0055     correlation/orbital -1.74  -1.49  -1.38

    (the a = 0.30 excess halves on 800 points -- grid error, not the model),
    while the two-term sum reads ~3x at 10 f_amp.  Near the carrier the phase
    part is the library Lorentzian (1.00022 from 0 to 3 linewidths) and the
    correlation ~1e-6 of it.
    """
    import warnings as _w
    for a in (0.0, 0.30):
        cir, pss = _a9_vdp(cval=4.0, lval=0.25, a=a)
        pac = PAC(cir, toolkit=circuit.numeric)
        f0 = 1.0 / float(pss.period)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            _v, info = pss.ppv()
            f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
            offs = np.array([3.0, 10.0, -10.0]) * f_amp
            ms = pac.modal_spectrum(pss, offs, 0, H=8, sidebands=16)
            pn = np.array([float(np.real(pac.pnoise(pss, f0 + o, 0,
                                                    maxsidebands=16)[0]))
                           for o in offs])
            old = (np.asarray(pac.oscillator_spectrum(
                pss, offs, 0, frequency_aware=False)[0], dtype=float)
                + np.asarray(pac.orbital_spectrum(pss, offs, 0, harmonic=1, H=8),
                             dtype=float))
        parts = ms['phase'] + ms['orbital'] + ms['correlation']
        assert np.max(np.abs(parts - ms['total'])) <= 1e-12 * np.max(ms['total'])
        ratio = ms['total'] / pn
        corr_over_orb = ms['correlation'] / ms['orbital']
        if a == 0.0:
            assert np.max(np.abs(ratio - 1.0)) < 2e-3, ratio
            assert np.max(np.abs(corr_over_orb)) < 0.1, corr_over_orb
            ## near the carrier: the phase part IS the Lorentzian, no correlation
            lw = np.pi * f0 ** 2 * float(pac.diffusion_constant(pss))
            near = np.array([0.0, lw])
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                mn = pac.modal_spectrum(pss, near, 0, H=8, sidebands=16)
                lor = np.asarray(pac.oscillator_spectrum(
                    pss, near, 0, frequency_aware=False)[0], dtype=float)
            assert np.max(np.abs(mn['phase'] / lor - 1.0)) < 1e-3, mn['phase'] / lor
            assert np.max(np.abs(mn['correlation'] / mn['phase'])) < 1e-5
        else:
            assert np.max(np.abs(ratio - 1.0)) < 0.02, ratio
            ## the correlation is what closes it: large and negative
            assert np.all(corr_over_orb < -1.0), corr_over_orb
            ## and without it the modal sum is far off (a presence claim)
            old_ratio = old[1:] / pn[1:]
            assert np.all(old_ratio > 2.0), old_ratio
            np.testing.assert_allclose(
                pac.correlation_spectrum(pss, offs, 0, H=8, sidebands=16),
                ms['correlation'], rtol=0, atol=0)
    with pytest.raises(ValueError, match='harmonic must be >= 1'):
        pac.modal_spectrum(pss, np.array([f_amp]), 0, harmonic=0, H=8)

#: the correlated pair of `_xcorr_oscillator`: white `A`, 1/f `B` per node,
#: correlation `RHO` in each part (white with white, 1/f with 1/f)
_XC_A1, _XC_B1, _XC_A2, _XC_B2, _XC_RHO, _XC_FR = 1e-6, 2e-6, 0.5e-6, 1.5e-6, 0.6, 0.15


class _XCOne(Circuit):
    """One node's noise current: a white part `wa (1 + mod u^2)` that
    follows the node's voltage and a 1/f part `fb fr/f` that does not."""
    terminals = ('a', 'b')
    instparams = [Parameter(name='wa', desc='', unit='', default=0.0),
                  Parameter(name='fb', desc='', unit='', default=0.0),
                  Parameter(name='mod', desc='', unit='', default=0.0)]

    def CY(self, x, w, epar=None):
        s = _XC_FR / (abs(float(w)) / (2.0 * np.pi))
        u = float(x[0] - x[1])
        p = (1.0 + self.iparv.mod * u * u) * self.iparv.wa + self.iparv.fb * s
        return self.toolkit.array(np.array([[p, -p], [-p, p]]))


def _xc_cross(uv, ux, mod, w):
    s = _XC_FR / (abs(float(w)) / (2.0 * np.pi))
    m1, m2 = 1.0 + mod * uv * uv, 1.0 + mod * ux * ux
    return _XC_RHO * (np.sqrt(m1 * m2 * _XC_A1 * _XC_A2)
                      + np.sqrt(_XC_B1 * _XC_B2) * s)


class _XCPair(Circuit):
    """The same two currents as ONE element, correlated: the reference --
    the per-element model splits it into its white and 1/f parts."""
    terminals = ('a', 'b')
    instparams = [Parameter(name='mod', desc='', unit='', default=0.0)]

    def CY(self, x, w, epar=None):
        s = _XC_FR / (abs(float(w)) / (2.0 * np.pi))
        uv, ux, mod = float(x[0]), float(x[1]), self.iparv.mod
        p1 = (1.0 + mod * uv * uv) * _XC_A1 + _XC_B1 * s
        p2 = (1.0 + mod * ux * ux) * _XC_A2 + _XC_B2 * s
        xx = _xc_cross(uv, ux, mod, w)
        return self.toolkit.array(np.array([[p1, xx], [xx, p2]]))


class _XCSub(SubCircuit):
    """Two single-node elements, and the cross term added by a `CY`
    OVERRIDE: this circuit's `CY` is not the sum of its elements'."""
    mod = 0.0

    def CY(self, x, w, epar=circuit.defaultepar):
        out = np.array(SubCircuit.CY(self, x, w, epar))
        iv, ix = self.get_node_index('v'), self.get_node_index('x')
        xx = _xc_cross(float(x[iv]), float(x[ix]), self.mod, w)
        out[iv, ix] += xx
        out[ix, iv] += xx
        return self.toolkit.array(out)


def _xcorr_oscillator(kind, mod):
    """The asymmetric lossy LC of `_slow_node_oscillator` (no slow node),
    its noise a correlated pair on `v` and `x`: 'joint' two elements and an
    override, 'ref' one element."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0 = 6.6634
    c = _XCSub() if kind == 'joint' else SubCircuit()
    c.mod = mod
    c.add_node('v'); c.add_node('x')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', 'x', L=1.0)
    c['Rl'] = R('x', gnd, r=0.2, noisy=False)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: (u - u ** 3 / 3.0) + 0.25 * (u ** 2 - 2.0))
    if kind == 'joint':
        c['n1'] = _XCOne('v', gnd, wa=_XC_A1, fb=_XC_B1, mod=mod)
        c['n2'] = _XCOne('x', gnd, wa=_XC_A2, fb=_XC_B2, mod=mod)
    else:
        c['n'] = _XCPair('v', 'x', mod=mod)
    pss = PSS(c, method='gear', reltol=1e-11)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T0, timestep=T0 / 240, x0=x0, maxiterations=200)
    assert pss.converged
    return pss, PAC(c, toolkit=circuit.numeric), [str(n) for n in c.nodes].index('v')


def _xcorr_surfaces(pss, pac, ov, stationary):
    """Every coloured surface at the offsets the test reads (warnings go
    to the caller, which records them)."""
    f0 = 1.0 / float(pss.period)
    o = np.array([1e-3, 1e-2]) * f0
    out = {}
    ms = pac.modal_spectrum(pss, o[1:], ov)
    out['modal'] = np.array([ms['phase'][0], ms['orbital'][0],
                             ms['correlation'][0]])
    if stationary:
        return out
    out['pnoise'] = np.array([float(np.real(pac.pnoise(
        pss, o[1], ov, maxsidebands=16, cyclostationary=True)[0]))])
    out['resolved'] = np.asarray(pac.coloured_diffusion_resolved(pss, o))
    out['gamma'] = np.asarray(pac.coloured_diffusion(pss, o[:1]))
    out['lineshape'] = np.asarray(pac.oscillator_spectrum(
        pss, [0.0, o[1]], ov, fmin=1e-6 * f0)[0])
    return out


def test_noise_correlated_across_elements_is_one_joint_component():
    """A circuit whose `CY` is not the sum of its elements' -- noise
    CORRELATED across elements, added by a `CY` override -- was refused by
    every coloured surface but pnoise and the sample series until
    2026-09-26, and those two took ONE root of the whole `CY`.  Now the
    whole circuit is one element with a white and a coloured component
    (`NoiseComponents.colour_model`), so it must equal the SAME noise built as one
    element spanning both nodes, which the per-element model splits the
    same way -- bit for bit, measured, stationary and modulated.  pnoise's
    old joint root was 4.0e-5 off that reference with the white part
    following the orbit and the 1/f part not.  The correlation is not
    small here: against the same pair uncorrelated it moves pnoise +59 %,
    `coloured_diffusion` +59 %, the resolved `c` +49 / +21 %, the lineshape
    -23 / +21 % and the modal terms 5 .. 23 %.  Poisoned: the cross term
    dropped from the joint fit alone fails (first at the stationary modal
    part, 1.6e-6, which reads the model only for its white width); the
    refusals restored fail."""
    import warnings as _w
    for mod, stationary in ((0.0, True), (0.5, False)):
        got = {}
        for kind in ('joint', 'ref'):
            pss, pac, ov = _xcorr_oscillator(kind, mod)
            with _w.catch_warnings(record=True) as rec:
                _w.simplefilter('always')
                got[kind] = _xcorr_surfaces(pss, pac, ov, stationary)
            said = any('not the sum of its elements' in str(r.message)
                       for r in rec)
            assert said == (kind == 'joint'), (kind, mod)
        for k, ref in got['ref'].items():
            err = float(np.max(np.abs(got['joint'][k] - ref) / np.abs(ref)))
            assert err <= 1e-12, (mod, k, err)
