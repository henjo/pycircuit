"""Shooting tests: shooting phase noise.  Split out of test_analysis_shooting.py on
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
from pycircuit.circuit.tests._shooting_fixtures import (_lc_osc,
    _orbit_modulated_vdp,
    _raw_pair_integrals,
    _solve_vdp_noise,
    _vdp_asym,
    _vdp_ppv_method)


def test_the_diffusion_constant_matches_the_monte_carlo_measurement():
    """⚠ THE SCALAR THE WHOLE SPECTRUM IS BUILT ON, against a physical number.

    `c = (1/T) ∫ v₁ᵀ (CY/2) v₁ dt`, measured a completely different way —
    Monte Carlo on the full nonlinear circuit, phase read from
    zero-crossing timing so no PPV appears in the measurement at all —
    which gives `7.7083e-08` against `7.9516e-08` predicted, ratio 1.0316
    inside the measurement's 4.1% uncertainty at N=1200.

    ⚠ AN EARLIER VERSION OF THIS TEST PASSED AT 1.5903e-07, EXACTLY TWICE,
    against a Monte Carlo injecting `Var(i) = CY/h`. That measurement
    carried the same one-sided-as-two-sided convention the code did, so it
    confirmed the bug rather than catching it, at a ratio of 0.9965.
    `kT/C` — external to both — showed that injection reproducing 1.92× the
    right answer over ten runs. A measurement built on the assumption under
    test cannot test it.
    """
    _cir, pss, pac = _solve_vdp_noise()
    c = pac.diffusion_constant(pss)
    assert abs(c - 7.9516e-08) < 0.05 * 7.9516e-08, \
        'c = %.6e against 7.9516e-08, which a Monte Carlo with the ' \
        'CORRECT injection (Var(i) = CY/2h) measured at 7.7083e-08 -- ' \
        'ratio 1.0316, inside its 4.1%% uncertainty' % c
    ## ⚠ THE OLD VALUE HERE WAS 1.5903e-07, EXACTLY TWICE THIS, and it
    ## passed against a Monte Carlo that injected `Var(i) = CY/h`. That
    ## measurement carried the same one-sided-as-two-sided error as the
    ## code, so it confirmed the bug instead of catching it. `kT/C` --
    ## external to both -- showed that injection reproducing 1.92x the
    ## right answer over ten runs. A measurement built on the assumption
    ## under test cannot test it.
    assert abs(c * 2 - 1.5903e-07) < 0.05 * 1.5903e-07, \
        'the pre-fix value is no longer exactly twice this one; if the ' \
        'convention changed again, re-derive rather than re-fit'


def test_coloured_upconversion_needs_asymmetry_AND_loss():
    """⚠ FIVE DISCRIMINATIONS THE QUADRATIC FUNCTIONAL CANNOT PRODUCE.

    A white source contracts the MEAN OF THE SQUARE of the PPV; a coloured
    one contracts the SQUARE OF THE MEAN. Same vector, same matrix, mean
    and square exchanged — and substituting one for the other returns a
    plausible non-zero number rather than an error.

    So the gate is a pattern of ZEROS, which the quadratic functional
    cannot fake because it is large in every row:

        a      Rs      Gamma/c      c
        0.00   0.00    2.4e-22      7.95e-08
        0.00   0.20    9.7e-23      1.01e-07
        0.25   0.00    4.9e-23      8.22e-08
        0.25   0.05    2.1e-04      8.65e-08
        0.25   0.20    4.1e-03      1.20e-07

    ⚠ TWO INDEPENDENT MECHANISMS FORCE THE ZERO, and only one is the one
    designers know:

    * a SYMMETRIC waveform gives `<v> = 0` — Hajimiri & Lee, and the
      reason symmetry is the first thing reached for in a VCO;
    * a LOSSLESS LC TANK gives `<v>[0] = 0` STRUCTURALLY, whatever the
      waveform does. `v` behaves as `CᵀV₁` and `dv/dt = Gᵀv₁`, whose
      inductor row is exactly `v[0]`; periodicity of `v[1]` then forces
      `∫v[0] dt = 0`. A property of the TOPOLOGY, not of the orbit.

    That second one is why van der Pol reports zero at every asymmetry,
    and why it is useless as a positive fixture and perfect as a negative
    one. A gate built only on van der Pol would have passed an
    implementation that returns zero always.
    """
    rows = []
    for a, rs in ((0.0, 0.0), (0.0, 0.2), (0.25, 0.0),
                  (0.25, 0.05), (0.25, 0.2)):
        _cir, pss, pac = _lc_osc(a=a, rs=rs)
        c = pac.diffusion_constant(pss)
        gam = float(pac.coloured_diffusion(pss, [1.0 / pss.period])[0])
        ## ⚠ THE STRUCTURAL ZERO IS A DISCRETE IDENTITY OF THE RAW PAIR.
        ## `samples` is the pair-consistent contraction (see `ppv`), whose
        ## mean on the lossless tank is O(h^2): 8.8e-6 |v| at 240 points,
        ## 2.2e-6 at 480, so its `Gamma/c` is ~5e-10 here.  The identity
        ## is pinned on the raw pair and the consistent object is held at
        ## its measured order; both are seven decades under the rows that
        ## upconvert.
        ints, _info = _raw_pair_integrals(pss, _cir)
        vraw = np.asarray(ints) / float(pss.period)
        cy = 0.5 * np.real(pac._cy_reduced(pss, 2.0 * np.pi / pss.period))
        gam_raw = float(vraw @ cy @ vraw)
        rows.append((a, rs, c, gam, gam_raw))
        assert c > 1e-8, 'a=%r rs=%r: c = %.3e, the white functional ' \
            'should be large in EVERY row or the zeros below prove ' \
            'nothing' % (a, rs, c)
    for a, rs, c, gam, gam_raw in rows[:3]:
        assert gam_raw / c < 1e-15, \
            'a=%r rs=%r gives raw-pair Gamma/c = %.3e; with either the ' \
            'symmetry or the lossless identity intact this must vanish' \
            % (a, rs, gam_raw / c)
        assert gam / c < 1e-8, \
            'a=%r rs=%r gives Gamma/c = %.3e for the consistent object, ' \
            'whose mean here is O(h^2) -- measured 5e-10' % (a, rs, gam / c)
    for a, rs, c, gam, _graw in rows[3:]:
        assert gam / c > 1e-5, \
            'a=%r rs=%r gives Gamma/c = %.3e; breaking BOTH must ' \
            'upconvert' % (a, rs, gam / c)
    assert rows[4][3] / rows[4][2] > rows[3][3] / rows[3][2], \
        'more tank loss must upconvert more, not less'


def test_gamma_never_exceeds_c_at_the_same_density():
    """Cauchy-Schwarz on the weighted mean, asserted exactly.

    `(Σhᵢxᵢ/T)² ≤ (Σhᵢxᵢ²)/T`, so the square of the mean never exceeds
    the mean of the square. Both functionals use the SAME quadrature here
    deliberately, which makes this hold at the discrete level rather than
    only in the limit — an assertion, not an expectation. Equality would
    mean the PPV is constant over the orbit, i.e. no orbit.
    """
    for a, rs in ((0.0, 0.0), (0.25, 0.2)):
        _cir, pss, pac = _lc_osc(a=a, rs=rs)
        c = pac.diffusion_constant(pss)
        gam = float(pac.coloured_diffusion(pss, [1.0 / pss.period])[0])
        ## ⚠ PRECONDITION (2026-09-05): the inequality gam <= c is
        ## satisfied vacuously by two near-zeros, so pin c away from zero
        ## first -- otherwise a degenerate fixture would pass it.
        assert c > 1e-8, \
            'a=%r rs=%r: c = %.3e is at the floor; the inequality below ' \
            'would be two zeros agreeing' % (a, rs, c)
        assert gam <= c * (1.0 + 1e-12), \
            'a=%r rs=%r: Gamma = %.6e exceeds c = %.6e, which is ' \
            'arithmetically impossible at one CY — the two functionals ' \
            'are not sharing a quadrature' % (a, rs, gam, c)


def test_the_phase_psd_refuses_below_the_lorentzian_corner():
    """A validity boundary, not a conditioning one.

    Below the corner the excess phase is a Wiener process whose spectrum is
    singular at the origin; the finite value the real lineshape attains
    comes from the NONLINEAR phase-to-voltage map, which
    `oscillator_spectrum` carries and this does not. Returning a large
    number there would be the mistake this object invites.
    """
    _cir, pss, pac = _lc_osc()
    c = pac.diffusion_constant(pss)
    corner = np.pi * (1.0 / float(pss.period)) ** 2 * c
    with pytest.raises(ValueError, match='Lorentzian corner'):
        pac.phase_psd(pss, [corner * 0.5])
    with pytest.raises(ValueError, match='diverges at zero offset'):
        pac.phase_psd(pss, [0.0])
    ## and just above it is fine
    assert np.all(np.isfinite(pac.phase_psd(pss, [corner * 10.0])))


def test_a_flicker_source_gives_a_one_over_f_cubed_skirt():
    """⚠ THE SLOPE IS THE ASSERTION — 30 dB/decade, not 20.

    Kundert: "S_u(f) is generally pink or proportional to 1/f. Then
    S_phi(f) would be proportional to 1/f³ at low frequencies." With
    `CY ∝ 1/f` the coloured term carries one more power of `f` than the
    white one, so the skirt steepens from `1/f²` to `1/f³` below the
    flicker corner and returns to `1/f²` above it.

    ⚠ A SLOPE, NOT A STATE. Demir 1996 synthesises 1/f through a
    Lorentzian filter network at one state variable per decade because Itô
    theory admits only white driving noise — an artefact of the SDE
    formulation. No SDE is formed here, so nothing is synthesised and the
    PSS never sees a filter. The corner appearing at the right place is
    what says the frequency dependence went in correctly.
    """
    _cir, pss, pac = _lc_osc(a=0.25, rs=0.2, flicker=True,
                             fref=1.0 / 6.66)
    ## ⚠ STARTS AT 1e-5, NOT 1e-6. The first version of this test swept to
    ## 1e-6 Hz, where `2 f S_phi = 3.10` -- the skirt was carrying three
    ## times the carrier's total power. `phase_psd` now refuses there; the
    ## test was wrong, not the refusal. See the power-bound test below.
    offs = np.logspace(-5, -1, 26)
    ## the DC phase model's slopes: frequency-aware (the default since
    ## 2026-09-26) steepens the top of this sweep further, -2.16, as the
    ## amplitude mode filters above f_amp -- a different statement
    S = pac.phase_psd(pss, offs, frequency_aware=False)
    slope = np.diff(np.log10(S)) / np.diff(np.log10(offs))
    assert slope[0] < -2.9, \
        'low-offset slope is %.3f decades/decade; a 1/f source must give ' \
        '1/f^3 there, and -2 would mean the frequency dependence never ' \
        'reached CY' % slope[0]
    assert slope[-1] > -2.1, \
        'high-offset slope is %.3f; above the flicker corner the white ' \
        'term must dominate and the skirt return to 1/f^2' % slope[-1]
    ## the corner is where the two terms are equal, and it must sit inside
    ## the swept band rather than at an end of it
    assert -2.9 < np.median(slope) < -2.1, \
        'the sweep does not straddle the flicker corner (median slope ' \
        '%.3f), so neither asymptote is being tested against the other' \
        % np.median(slope)
    ## and a WHITE source on the same circuit must not steepen
    _c2, pss2, pac2 = _lc_osc(a=0.25, rs=0.2, flicker=False)
    S2 = pac2.phase_psd(pss2, offs, frequency_aware=False)
    sl2 = np.diff(np.log10(S2)) / np.diff(np.log10(offs))
    assert abs(sl2.min() + 2.0) < 1e-6 and abs(sl2.max() + 2.0) < 1e-6, \
        'a white source gave slopes in [%.4f, %.4f]; it must be exactly ' \
        '1/f^2 everywhere or the 1/f^3 above is not evidence' \
        % (sl2.min(), sl2.max())


def test_the_phase_psd_refuses_a_skirt_carrying_more_than_unit_power():
    """⚠ A SECOND FLOOR, INDEPENDENT OF THE LORENTZIAN CORNER — and for a
    coloured source it is the binding one by 306×.

    The normalised lineshape integrates to 1, and the integral over one box
    of width `df` on each side is a lower bound on it, so

        2·df·S_φ(df) ≤ 1

    is NECESSARY for the linearised skirt to be consistent with unit power.
    Vanassche, Gielen & Sansen (2003) §6 derive the same statement for a
    `1/f` input and reduce it to `df_c ≥ ε·f₀·√(2·f_1f)`; the form asserted
    here needs no assumption about the source's colour.

    ⚠ IT WAS A LIVE DEFECT, NOT A HYPOTHETICAL. The first version of
    `test_a_flicker_source_gives_a_one_over_f_cubed_skirt` swept to 1e-6 Hz
    where `2 f S_φ = 3.10` — three times the carrier's total power — and
    every assertion in it passed. The Lorentzian corner sat at 8.2e-09 Hz,
    **306× too permissive**, because it is built from `c` alone and knows
    nothing about a `Γ(f)` that grows as the offset falls.

    ⚠ AND IT IS A LOWER BOUND ON THE BREAKDOWN, NOT THE BREAKDOWN. On
    Vanassche's own example the observed flattening is at ~300 Hz, 3× the
    bound. So this refuses what is definitely invalid and admits a band
    that is already suspect — deliberately, because refusing at 3× would be
    fitting a threshold to a single example.
    """
    _cir, pss, pac = _lc_osc(a=0.25, rs=0.2, flicker=True, fref=1.0 / 6.66)
    with pytest.raises(ValueError, match='times the TOTAL power'):
        pac.phase_psd(pss, np.logspace(-6, -1, 26))
    ## and the two floors are genuinely different numbers.  ⚠ The corner
    ## is built from the WHITE functional read at the carrier -- which is
    ## what `diffusion_constant` silently returned for this 1/f source
    ## until it started refusing colour, as its docstring always said.
    c = pac._white_diffusion_at(pss, 2.0 * np.pi / float(pss.period))
    corner = np.pi * (1.0 / float(pss.period)) ** 2 * c
    ok = pac.phase_psd(pss, np.logspace(-5, -1, 26))
    assert np.all(2.0 * np.logspace(-5, -1, 26) * ok < 1.0)
    assert corner < 1e-7, \
        'the Lorentzian corner is %.3e; if it had risen to meet the power ' \
        'bound the two floors would no longer be independent' % corner


def test_the_power_bound_reproduces_vanassches_worked_example():
    """`df_c ≥ ε·f₀·√(2·f_1f)` — their closed form, from ours.

    Their §6 substitutes the traditional characteristic
    `S(df) ≈ ε²(f₀²/df²)S_n(df)` into the normalisation `∫S = 1` and bounds
    the integral below by one box each side, giving
    `1 ≥ 2ε²(f₀²/df_c)S_n(df_c)`. With `S_n = f_1f/df` that is
    `df_c ≥ ε·f₀·√(2·f_1f)`.

    ⚠ ASSERTED AS AN IDENTITY BETWEEN TWO FORMS, not as a transcription.
    `2·f·S_φ ≤ 1` is the general statement; their result is its `1/f`
    special case, and the two must agree exactly rather than approximately.
    At `ε² = 1e-19`, `f₀ = 1 GHz`, `f_1f = 50 kHz` the paper says "≥ 100 Hz"
    and this gives 100.000 Hz. Their observed flattening is ~300 Hz, which
    is the factor-of-three headroom the docstring warns about.
    """
    eps2, f0, f1f = 1e-19, 1e9, 50e3
    ## the general form: 2 f S_phi(f) = 1 with S_phi = f0^2 eps^2 f1f / f^3
    general = np.sqrt(2.0 * f0 ** 2 * eps2 * f1f)
    ## their closed form
    theirs = np.sqrt(eps2) * f0 * np.sqrt(2.0 * f1f)
    assert abs(general / theirs - 1.0) < 1e-12, \
        'the general power bound %.6g and Vanassche\'s closed form %.6g ' \
        'are not the same statement' % (general, theirs)
    assert abs(theirs - 100.0) < 1e-9, \
        'their worked example gives %.6f Hz against the ">= 100 Hz" ' \
        'printed in the paper' % theirs


class _BlueNoise(IS):
    """A source whose PSD RISES as `f²` — enough to break the power bound's
    monotonicity precondition, and nothing else in the tree can.
    """

    instparams = IS.instparams + [
        Parameter(name='fref', desc='Reference', unit='Hz', default=1.0)]

    def CY(self, x, w, epar=None):
        f = abs(float(w)) / (2.0 * np.pi)
        p = self.iparv.noisePSD * (f / self.iparv.fref) ** 4
        return self.toolkit.array([[p, -p], [-p, p]])

    ## ⚠ `f^4` RATHER THAN `f^2`, and the offsets below sit ABOVE `f0`.
    ## `S_phi ~ (c + Gamma(f))/f^2` needs `Gamma` to grow faster than `f^2`
    ## to turn the spectrum upward, AND `diffusion_constant` samples `CY`
    ## at the single frequency `f0` -- so a rising source makes `c` huge
    ## and the white term dominates every offset BELOW `f0`. The first
    ## version of this fixture swept 1e-4..1e-1 with `f0 = 0.15` and the
    ## guard never fired, because `c` was 53.3.


def test_the_power_bound_refuses_when_its_own_derivation_does_not_apply():
    """⚠ THE BOUND ADDED AN HOUR AGO HAD AN UNSTATED PRECONDITION.

    `2·Δf·S(Δf) ≤ ∫_{-Δf}^{+Δf} S ≤ 1` — and the FIRST inequality needs
    `S(f) ≥ S(Δf)` for every `|f| ≤ Δf`. The spectrum must not dip below
    its edge value anywhere further in. That holds for a monotone skirt,
    for the flattened near-carrier shape, and even with a spur (which
    *adds* power inside rather than creating a dip).

    ⚠ IT FAILS FOR A LOCKED PLL, whose phase-noise transfer function is
    HIGH-PASS: suppressed at DC, rising to the free-running level beyond
    the loop bandwidth, so it dips below its edge value everywhere inside.
    The bound is not thereby shown to be *violated* there — total power is
    still 1 — it is **no longer derived**, and a floor that is not derived
    cannot be used as one.

    ⚠ THAT IS §D SHAPE 0e A SECOND TIME, in a bound written an hour after
    shape 0e was written up: a result asserted outside the conditions its
    own derivation assumes. A correction is not self-certifying, and
    neither is a generalisation.

    Unreachable through a driven circuit today, because `phase_psd`
    refuses those — so the reachable case is a source whose density grows
    faster than `f²`, which is what `_BlueNoise` is for. The check is on
    the *shape of the returned spectrum*, so it will catch the PLL case
    when driven oscillators land, without needing to know about loops.
    """
    _cir, pss, pac = _lc_osc(a=0.25, rs=0.2, npts=240)
    ## sanity: the ordinary case is monotone and passes
    offs = np.logspace(-0.5, 1.5, 12)
    assert np.all(np.diff(pac.phase_psd(pss, offs)) < 0)

    cir2 = SubCircuit()
    cir2.add_node('v')
    cir2['C'] = C('v', gnd, c=1.0)
    cir2['B'] = BSource('v', gnd, gnd, 'v',
                        i_func=lambda u: 1.0 * (u - u ** 3 / 3.0)
                        + 0.25 * (u ** 2 - 2.0))
    cir2.add_node('x')
    cir2['L'] = L('v', 'x', L=1.0)
    cir2['Rs'] = R('x', gnd, r=0.2)
    cir2['n'] = _BlueNoise('v', gnd, i=0.0, noisePSD=1e-6, fref=1.0)
    import warnings
    pss2 = PSS(cir2, method='gear', reltol=1e-12)
    x0 = np.zeros(cir2.n - 1)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss2.solve(period=6.66, timestep=6.66 / 240, x0=x0, maxiterations=80)
    assert pss2.converged
    pac2 = PAC(cir2, toolkit=circuit.numeric)
    with pytest.raises(ValueError, match='RISES with offset'):
        pac2.phase_psd(pss2, offs)


def test_pnoise_is_phase_psd_times_the_CARRIER_POWER_not_a_psd_convention():
    """⚠⚠ THE "FACTOR OF TWO" BETWEEN `pnoise` AND `c f0^2/df^2` IS `A^2/2`,
    AND VAN DER POL MAKES IT LOOK LIKE A CONVENTION.

    Kundert section 3.5 eq (15): `L(df) = c f0^2 / df^2` for
    `f_delta << df << f0`. Our `phase_psd` matches that EXACTLY. But `pnoise`
    returns OUTPUT VOLTAGE noise, V^2/Hz, not phase noise, rad^2/Hz -- and
    the conversion is the CARRIER POWER `A^2/2`.

    ⚠ VAN DER POL'S AMPLITUDE IS 2, SO `A^2/2 = 2`, numerically
    indistinguishable from the one-sided/two-sided factor that section 0g
    caught for real. It was recorded as "a loose end, the same factor-of-two
    family". It is not: it is a fixture coincidence, section D shape 0i.

    Measured across a 36x range of carrier power -- the ratio MOVES with
    amplitude, which is what a convention factor could not do:

        A        A^2/2      pnoise/(c f0^2/df^2)   /(A^2/2)   phase_psd/L
        0.9998   0.49975    0.499334               0.999164   1.00000000
        1.9995   1.99901    1.997338               0.999164   1.00000000
        3.9990   7.99604    7.989351               0.999164   1.00000000
        5.9985   17.99108   17.976041              0.999164   1.00000000

    and the residual 0.999164 is DISCRETISATION, converging to 1 as the grid
    refines: 0.995796 / 0.999164 / 0.999870 / 1.000012 at npts = 120 / 240 /
    480 / 960. Nothing is left unexplained.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    Q = 8.0
    mu = 1.0 / (2 * np.pi * Q)

    def run(sscale, npts):
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u, _s=sscale: mu * (u - u ** 3 / (3.0 * _s * _s)))
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        pss = PSS(cir, method='gear', reltol=1e-12)
        x0 = np.zeros(cir.n - 1)
        x0[0] = 2.0 * sscale
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=2 * np.pi, timestep=2 * np.pi / npts, x0=x0,
                      maxiterations=250)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        f0 = 1.0 / float(pss.period)
        cc = pac.diffusion_constant(pss)
        X = np.asarray(pss.waveform[1], dtype=float)
        amp = float(np.max(np.abs(X[0 if 0 < pss.irefnode else 1])))
        df = f0 * 1e-6
        S_pn, _sb = pac.pnoise(pss, f0 + df, 0, sweeptype='absolute')
        S_ph = float(np.asarray(pac.phase_psd(pss, np.array([df]))).ravel()[0])
        kundert = cc * f0 * f0 / (df * df)
        return amp, S_pn, S_ph, kundert

    ## ⚠ phase_psd IS Kundert eq (15), exactly, at both amplitudes
    for sscale in (0.5, 2.0):
        amp, S_pn, S_ph, kundert = run(sscale, 480)
        assert abs(S_ph / kundert - 1.0) < 1e-9, \
            'phase_psd must equal c f0^2/df^2 exactly; got %.9f at A=%.4f' \
            % (S_ph / kundert, amp)
        ## and pnoise is that times the CARRIER POWER
        assert abs((S_pn / kundert) / (amp * amp / 2) - 1.0) < 5e-3, \
            'pnoise/(c f0^2/df^2) must be A^2/2; got %.6f against %.6f' \
            % (S_pn / kundert, amp * amp / 2)

    ## ⚠⚠ AND THE RATIO MUST MOVE WITH AMPLITUDE, or this test cannot tell a
    ## carrier power from a PSD convention -- which is the whole point
    a_lo, pn_lo, _p, k_lo = run(0.5, 480)
    a_hi, pn_hi, _p2, k_hi = run(2.0, 480)
    moved = (pn_hi / k_hi) / (pn_lo / k_lo)
    expect = (a_hi * a_hi) / (a_lo * a_lo)
    assert abs(moved / expect - 1.0) < 1e-2, \
        'the ratio must scale as A^2 (%.4f expected, %.4f seen) -- a PSD ' \
        'convention would be CONSTANT' % (expect, moved)
    assert moved > 4.0, \
        'and it must move enough to be unmistakable; got %.4f' % moved


def test_ppv_runs_under_every_integrator_and_agrees_with_gear():
    """⚠ B8 WIRED THROUGH. Building the plain transposed replay was not
    enough: `ppv` reached past `FactoredPeriod` to
    `_monodromy_matvec_transposed` DIRECTLY and refused on `fp.kind`, so
    the machinery shipped and every adjoint surface still said
    "re-solve with method='gear'". The refusal is gone and the calls go
    through the dispatcher.

    ⚠⚠ COMPARE `v[:m]`, NOT `v`. The solved-history map's state is the
    PAIR, so `ppv` returns `2m` components under gear and `m` under a
    one-step method. `norm(v)` therefore compares DIFFERENT OBJECTS across
    methods and shows a spurious 5% disagreement that converges cleanly on
    both sides — which is what makes it dangerous rather than obvious.
    `v[:m]` is the differential block under both.

    Measured: `|v[:m]|` → 0.500008 (gear/800) against 0.500811 (trap/800),
    and the phase diffusion constant `c` agrees at `O(h²)` — trap against
    gear 7.0e-2, 1.7e-2, 4.1e-3 at 200/400/800, ratios 4.2 and 4.1.
    """
    m = None
    vs, cs = {}, {}
    for method in ('gear', 'trap'):
        for npts in (200, 400, 800):
            cir, pss = _vdp_ppv_method(method, npts)
            m = cir.n - 1
            v, _info = pss.ppv()
            v = np.asarray(v, dtype=float)
            vs[(method, npts)] = float(np.linalg.norm(v[:m]))
            cs[(method, npts)] = float(PAC(cir).diffusion_constant(pss))

    ## The differential block CONVERGES to gear's; the FULL vector would
    ## not, and that is a width artifact rather than a defect.
    ##
    ## ⚠ ASSERTED AS A RATE, NOT A BOUND. The first version demanded
    ## 5e-3 at every grid and failed at 400 points with 7.65e-3 — where
    ## widening the bound would have hidden the only interesting fact,
    ## which is that the gap is second order: 7.65e-3 then 1.61e-3, ratio
    ## 4.8. A constant offset between the two maps would pass a loose
    ## bound and fail this.
    rel = [abs(vs[('trap', n)] - vs[('gear', n)]) / vs[('gear', n)]
           for n in (200, 400, 800)]
    assert rel[2] < 3e-3, \
        'trap and gear disagree on |v[:m]| by %.3e at the finest grid' \
        % rel[2]
    assert 2.5 < rel[1] / rel[2] < 8.0, \
        'the |v[:m]| gap is not closing at O(h^2) (%s); a gap that stops ' \
        'shrinking is two different objects, not two discretisations of ' \
        'one' % (['%.3e' % r for r in rel],)

    ## and `c` -- the physical quantity -- converges to gear at O(h^2)
    ref = cs[('gear', 800)]
    e = [abs(cs[('trap', n)] - ref) / abs(ref) for n in (200, 400, 800)]
    assert e[2] < 1e-2, \
        'trap\'s diffusion constant is %.3e off gear at 800 points' % e[2]
    for a, b in zip(e, e[1:]):
        assert 2.5 < a / b < 6.0, \
            'c is not converging at O(h^2) across methods (%s, ratios ' \
            '%.2f); if the rate has changed the two maps are no longer ' \
            'discretising the same object' % (e, a / b)


def test_c_agrees_between_the_ppv_form_and_the_swept_noise_path():
    """⚠ TWO INDEPENDENT CODE PATHS TO ONE PHYSICAL CONSTANT.

    This is the cross-check Kundert actually describes (*Introduction to
    RF Simulation*, p. 11–12), which is NOT a second corner formula. He
    gives the corner as `fΔ = cπf₀²` — identical to `phase_psd`'s — and
    says the small-signal sweep *"does not show the roll off"* but that
    *"it is possible to use (15) to determine fΔ"*, i.e. fit `c` from the
    far `1/Δf²` skirt and apply the same formula. So the independently
    checkable object is **`c`**, not the corner.

    The two routes share the PSS and the noise sources and nothing else:

      A. `diffusion_constant()` — the PPV contracted against `CY`.
      B. the swept `pnoise` skirt, normalised to carrier power, via
         `L(Δf) = c f₀²/Δf²`.

    Measured, van der Pol at 1600 points:

        Δf/f₀      c from the skirt
        1e-2       7.510951726e-08    ← outside the valid window
        3e-3       6.381291626e-08
        1e-3       6.263389371e-08
        3e-4       6.251214799e-08
        1e-4       6.250576334e-08    ←→ 6.250576786e-08 from path A

    ⚠⚠ **THE ESTIMATOR MUST BE THE LIMIT, NOT AN AVERAGE.** A first
    version of this took the median across all five offsets and reported
    the two paths agreeing to 0.2 % — which is not a measurement of the
    disagreement, it is a measurement of how many invalid offsets were
    included. Kundert states the window as `fΔ ≪ Δf ≪ f₀`; the outermost
    point here is 20 % high because it is outside it, not because either
    path is wrong.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 600, x0=np.array([2.0, 0.0]),
                  maxiterations=300)
    assert pss.converged

    pac = PAC(cir)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cA = float(pac.diffusion_constant(pss))
    f0 = 1.0 / pss.period

    ## carrier power in the fundamental, from the PSS waveform's own DFT.
    ## ⚠ `A²/2`, the CARRIER POWER -- not a PSD convention; see the entry
    ## on the "factor of two" that this normalisation once looked like.
    X = np.asarray(pss.waveform[1], dtype=float)[0][:-1]
    A1 = 2.0 * np.abs(np.fft.rfft(X)[1]) / len(X)
    Pc = 0.5 * A1 * A1
    assert abs(Pc - 2.0) < 1e-3, \
        'van der Pol amplitude 2 gives carrier power 2; got %.6f' % Pc

    ob = [str(nd) for nd in cir.nodes].index('v')
    got = {}
    for k in (1e-3, 1e-4):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            S, _ = pac.pnoise(pss, f0 * (1.0 + k), ob, sweeptype='absolute')
        df = k * f0
        got[k] = float(np.real(S)) / Pc * df * df / (f0 * f0)

    ## deep in the skirt the two paths must agree tightly.  ⚠ THE EXACT
    ## VALUE IS KNOWN HERE: the scipy adjoint of the exact monodromy gives
    ## c = 6.250850e-08 for this fixture; `diffusion_constant` (the
    ## pair-consistent PPV, see `ppv`) is 2.9e-5 above it and the swept
    ## path 1.4e-4 below it, each at its own discretisation, so the two
    ## sit 1.7e-4 apart and the bound is set from that, not from either.
    rel = abs(got[1e-4] - cA) / cA
    assert rel < 3e-4, \
        'the swept-noise `c` (%.9e) and diffusion_constant (%.9e) disagree ' \
        'by %.3e at Delta f/f0 = 1e-4. These are independent paths -- the ' \
        'PPV quadratic form against adjoint sideband propagation -- so a ' \
        'disagreement is a normalisation error in one of them, which is ' \
        'the class of defect a kT/C reference once caught in `c` itself' \
        % (got[1e-4], cA, rel)

    ## ⚠ AND IT MUST IMPROVE AS THE WINDOW IS ENTERED, which is what says
    ## the residual is the window rather than a constant offset
    assert abs(got[1e-4] - cA) < abs(got[1e-3] - cA), \
        'the skirt estimate does not converge toward diffusion_constant ' \
        'as the offset enters Kundert\'s window (%.9e at 1e-3, %.9e at ' \
        '1e-4, against %.9e)' % (got[1e-3], got[1e-4], cA)


def test_the_resolved_phase_diffusion_keeps_its_order_on_a_non_uniform_grid():
    """⚠ `coloured_diffusion_resolved` (and so `phase_psd`) paired PPV sample
    `j`, taken at `t_j`, with `t_{j+1}` (`tms[1:1 + n]`, weights from the
    shifted times) until 2026-09-26.  A uniform grid cannot see it (a
    common phase: white bit-identical, a Lorentzian 3.3e-16).  On a grid
    whose step varies smoothly (3:1, van der Pol a = 0.3, radau), against a
    uniform radau run at 1600 points:

        source      pairing     N=200      N=400      N=800
        white       t_{j+1}    +8.6e-5    +2.4e-5    +7.0e-6   (2nd order)
                    t_j        +1.4e-8    +5.8e-10   -3.2e-11
        Lorentzian  t_{j+1}    -9.1e-4    -5.7e-4    -3.2e-4   (below 1st)
                    t_j        +1.3e-8    +4.9e-10   -3.9e-11

    A coloured source is hit harder: its harmonics carry different weights
    `CY(|f - l f0|)`, so the per-sample phase error moves energy between
    them at first order, where for a white one it cancels in the sum.
    Gated cheaply at 200 points: the white fold against the same grid's
    `c` (Parseval), the Lorentzian against the uniform grid."""
    import warnings as _w

    def build(source, nonuniform):
        circuit.default_toolkit = circuit.numeric
        mu = 1.0 / (2.0 * np.pi * 8.0)
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.3 * u * u)
        T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
        c['n'] = (IS('v', gnd, i=0.0, noisePSD=1e-6) if source == 'white' else
                  IS('v', gnd, i=0.0, noisePSD=1e-6, noiseTau=0.3 * T))
        grid = None
        if nonuniform:
            w = 1.0 + 0.5 * np.sin(2.0 * np.pi * np.arange(200) / 200)
            grid = w / w.sum()
        pss = PSS(c, method='radau', reltol=1e-12)
        x0 = np.zeros(c.n - 1)
        x0[0] = 2.0
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 200, x0=x0, maxiterations=300,
                      grid=grid)
        assert pss.converged
        return pss, PAC(c, toolkit=circuit.numeric)

    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss, pac = build('white', True)
        f0 = 1.0 / float(pss.period)
        freqs = np.array([1e-3, 1e-2, 0.1]) * f0
        cf = pac.coloured_diffusion_resolved(pss, freqs,
                                             frequency_aware=False)  # the DC fold's timing
        c = float(pac.diffusion_constant(pss))
        assert np.max(np.abs(cf / c - 1.0)) < 1e-7, cf / c - 1.0
        out = {}
        for nonuniform in (False, True):
            pss, pac = build('lorentz', nonuniform)
            out[nonuniform] = pac.coloured_diffusion_resolved(
                pss, freqs, frequency_aware=False)
    err = out[True] / out[False] - 1.0
    assert np.max(np.abs(err)) < 1e-6, err


def test_the_dc_colour_projection_takes_a_source_that_follows_the_orbit():
    """`coloured_diffusion` -- `Gamma(f) = vbar^T (CY/2) vbar`, the colour
    that reaches the phase through the PPV's time average -- refused a
    source that follows the orbit until 2026-09-26.  It is now the l = 0
    term of the modulated fold, `(s(f)/2) |<v_1^T G>|^2` per component.
    On the asymmetric orbit (a = 0.3):

      * a SIGNED flicker against its stationary realisation (a 1/f source
        times V_v): 1.1e-15 at 1e-4 .. 1e-2 f0;
      * never above `coloured_diffusion_resolved` (Jensen): Gamma / c(f) =
        1.000 / 0.9997 / 0.997 there -- a 1/f source reaches the phase
        almost wholly through the DC term near the carrier.

    ⚠ A modulated WHITE part enters through the root of its PSD, a
    convention: `g xi` and `|g| xi` are one white process, and Gamma's
    white share reads 0.0017 of `c` that way against 0.643 through a
    signed multiplier -- both under `c`, which alone is physical for white
    noise.  Poisons: the l = 1 row for the l = 0 one; the signed amplitudes
    ignored."""
    import warnings as _w
    res = {}
    for kind in ('flicker_ref', 'flicker'):
        _c, pss, pac, ov = _orbit_modulated_vdp(kind, a=0.3, kk=0.005)
        f0 = 1.0 / float(pss.period)
        fr = np.array([1e-4, 1e-3, 1e-2]) * f0
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res[kind] = pac.coloured_diffusion(pss, fr)
            ## Jensen bounds Gamma by the DC fold (frequency-aware can sit below)
            cres = pac.coloured_diffusion_resolved(pss, fr,
                                                   frequency_aware=False)
        assert np.all(res[kind] <= cres * (1.0 + 1e-12)), (kind, res[kind] / cres)
    err = np.max(np.abs(res['flicker'] / res['flicker_ref'] - 1.0))
    assert err < 1e-9, err


def test_the_frequency_aware_fold_reads_each_band_at_f_minus_l_f0():
    """The band convention of the frequency-aware fold, MEASURED: a source
    whose density peaks at `f0 + f` (Q = 20) as an element, against its
    exact realisation -- white noise through a parallel RLC resonant there,
    into the tank through a transconductor -- whose white source cannot
    tell the bands apart (`frequency_aware_diffusion`).  At `f = 0.05 f0`,
    element / realisation = 0.9999 with the bands at `f - l f0` (the DC
    fold's), 1.083 with `f + l f0`; the DC fold 1.048.  (The slow-node
    fixtures cannot tell: there the l = 0 term carries the colour.)"""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T6 = 6.66

    def build(kind, fp, Qf=20.0, g=1e-2, Pw=1e-4, Rr=1.0):
        wp = 2.0 * np.pi * fp
        Lr, Cr = Rr / (Qf * wp), Qf / (Rr * wp)
        P = g * g * Rr * Rr * Pw

        class _Res(Circuit):
            terminals = ('p', 'n')
            instparams = []

            def CY(self, x, w, epar=None):
                ww = max(abs(float(w)), 1e-300)
                p = P / (1.0 + Qf ** 2 * (ww / wp - wp / ww) ** 2)
                return self.toolkit.array(np.array([[p, -p], [-p, p]]))
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: (u - u ** 3 / 3.0) + 0.25 * (u ** 2 - 2.0))
        c.add_node('x')
        c['L'] = L('v', 'x', L=1.0)
        c['Rs'] = R('x', gnd, r=0.2, noisy=False)
        if kind == 'element':
            c['n'] = _Res('v', gnd)
        else:
            c.add_node('r')
            c['nw'] = IS('r', gnd, i=0.0, noisePSD=Pw)
            c['Rr'] = R('r', gnd, r=Rr, noisy=False)
            c['Lr'] = L('r', gnd, L=Lr)
            c['Cr'] = C('r', gnd, c=Cr)
            c['gm'] = BSource('r', gnd, gnd, 'v', i_func=lambda u, _g=g: _g * u)
        pss = PSS(c, method='radau', reltol=1e-11)
        x0 = np.zeros(c.n - 1)
        x0[0] = 2.0
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T6, timestep=T6 / 240, x0=x0, maxiterations=200)
        assert pss.converged
        return pss, PAC(c, toolkit=circuit.numeric)
    f0 = 1.0 / 6.884048287
    fs = 0.05 * f0
    pe, pace = build('element', f0 + fs)
    pr, pacr = build('real', f0 + fs)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        ce = float(pace.coloured_diffusion_resolved(pe, [fs])[0])
        cdc = float(pace.coloured_diffusion_resolved(pe, [fs],
                                                     frequency_aware=False)[0])
        cref = float(pacr.frequency_aware_diffusion(pr, fs))
    assert abs(ce / cref - 1.0) < 1e-3, ce / cref
    assert abs(cdc / cref - 1.0) > 2e-2, cdc / cref


def test_the_harmonic_resolved_fold_is_exactly_c_for_white():
    """PARSEVAL, ASSERTED AT ROUND-OFF: `sum_l V_l^H (CY/2) V_l = c`.

    With `CY` constant the per-harmonic fold is the time average of the
    same quadratic form, and with the SAME step-weighted quadrature the
    discrete identity is exact -- which is what pins the transform's
    normalisation, the one thing a Fourier fold can silently get wrong
    by `N`, `T`, or `2 pi`.

    ⚠ AND THE DOUBLE COUNT IS MEASURED ON THE ONE FIXTURE THAT HAS IT.
    `Gamma` is exactly the `l = 0` term, so `c + Gamma - c_res == Gamma`
    to round-off -- on van der Pol that is `1e-22 c` and proves nothing
    (the inductor shorts the tank node at DC, so its PPV averages to zero
    whatever the core does: the fixture shared the claim's assumption,
    §D 0b).  `_lc_osc(a=0.25, rs=0.2)` breaks both the symmetry and the
    lossless identity and carries `Gamma/c = 4e-3`, so there the retired
    `c + Gamma` form is 0.4% high for a WHITE source, and this test would
    have failed against it.
    """
    for a, rs, floor in ((0.0, 0.0, None), (0.25, 0.2, 1e-3)):
        _cir, pss, pac = _lc_osc(a=a, rs=rs)
        f0 = 1.0 / float(pss.period)
        c = pac.diffusion_constant(pss)
        offs = np.array([1e-3, 1e-4, 3e-2]) * f0
        cres = pac.coloured_diffusion_resolved(pss, offs,
                                               frequency_aware=False)  # the DC fold's identity (2026-09-26: frequency-aware is the default)
        assert np.all(np.abs(cres / c - 1.0) < 1e-12), \
            'a=%r rs=%r: the fold is %s against c = %.6e -- Parseval ' \
            'fails, so the transform normalisation is wrong' \
            % (a, rs, cres, c)
        gam = pac.coloured_diffusion(pss, offs)
        assert np.all(np.abs((c + gam - cres) - gam) < 1e-12 * c), \
            'c + Gamma - c_res is not Gamma: the l = 0 term is not what ' \
            'Gamma computes'
        if floor is not None:
            assert np.all(gam / c > floor), \
                'a=%r rs=%r: Gamma/c = %s -- the fixture cannot see the ' \
                'double count, so the assertion above is vacuous' \
                % (a, rs, gam / c)


def test_the_ppv_samples_are_pair_consistent_and_second_order():
    """⚠⚠ THE GEAR PAIR'S FIRST BLOCK WAS A FIRST-ORDER PPV, AND `c` WITH IT.

    Fixture: van der Pol plus `0.3 u^2` at `Q = 8`.  The even term makes
    the frequency bias-sensitive (period 6.28 -> 6.73), so a kick excites
    the amplitude mode and the phase keeps accumulating while it relaxes:
    the PPV is 100x van der Pol's and `v . xdot = 1` becomes a difference
    of two O(5) terms.  That cancellation is what turns an `O(h)` rotation
    of the pair's first block into 16.6% on `c` at 400 points.

    THE REFERENCE HAS NO SHOOTING CODE IN IT.  The circuit is the explicit
    ODE `vdot = mu (v - v^3/3) + a v^2 - i_L`, `i_Ldot = v`; DOP853 at
    rtol 1e-12 gives the orbit (period 6.730654, which the shooting
    periods 6.730950 / 6.730730 / 6.730673 at 400 / 800 / 1600 points
    extrapolate to at second order), the exact monodromy by the
    variational equations, its left null vector normalised `v . f = 1`,
    and `v(t)` by transport.  A kick instrument on the same ODE agreed
    with that adjoint to 1e-4 at eleven of sixteen phases (the rest were
    an event-count artefact, exactly one period per kick).  From it:

        c_true = <v_v(t)^2> CY/2 = 5.3703e-06        (CY = 1e-6, one-sided)

    ⚠ `pnoise` gave 5.355e-06 at 400 points ALL ALONG -- it contracts in
    pair space -- and the '14% deficit' first read as a physics term was
    `c` being 16.6% high.  Gates, each against that constant:
      the corrected `c` at 400 and 800 points, and its order;
      the invariant `v(t) . xdot(t) = 1` held ALONG THE ORBIT to 1e-3
        (the first block: std 2.7e-2), with `xdot` a central difference
        of the shooting waveform so no ODE is written into the test.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    C_TRUE = 5.3703e-06
    mu = 1.0 / (2.0 * np.pi * 8.0)
    errs = []
    for npts in (400, 800):
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0)
                           + 0.3 * u ** 2)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        pss = PSS(cir, method='gear', reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=6.731, timestep=6.731 / npts,
                      x0=np.array([2.0, 0.0]), maxiterations=300)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        c = pac.diffusion_constant(pss)
        errs.append(abs(c / C_TRUE - 1.0))
        ## the invariant along the orbit, with the waveform's own tangent
        v, info = pss.ppv()
        m = cir.n - 1
        S = np.asarray(info['samples'])[:, :m]
        X = np.delete(np.asarray(pss.waveform[1], dtype=float),
                      pss.irefnode, axis=0)
        T = float(pss.period)
        h = T / npts
        Xp = X[:, :-1]
        xd = (np.roll(Xp, -1, axis=1) - np.roll(Xp, 1, axis=1)) / (2.0 * h)
        n = min(S.shape[0], xd.shape[1])
        dots = np.array([S[j] @ xd[:, j] for j in range(n)])
        ## the mean carries the central difference's own O(h^2) bias
        ## (2.1e-3 at 400 points); the STD is the invariant's test
        assert abs(dots.mean() - 1.0) < 5e-3 and dots.std() < 1e-3, \
            'npts=%d: v(t).xdot(t) = %.5f +- %.1e along the orbit; the ' \
            'phase functional must hold it at every t' \
            % (npts, dots.mean(), dots.std())
    assert errs[0] < 4e-3, \
        'c at 400 points is %.2e from the exact 5.3703e-06 -- the ' \
        'first-block PPV gave 1.7e-1 here' % errs[0]
    assert errs[1] < 1.2e-3 and errs[0] / errs[1] > 3.0, \
        'errors %.2e -> %.2e over a doubling: second order is a ratio ' \
        'of 4, the first block gave 2' % (errs[0], errs[1])


def test_the_pair_consistent_ppv_is_second_order_on_a_DAE_too():
    """⚠ THE ALGEBRAIC STATE'S SLAVED COUPLING IS O(h), AND THE FIRST BUILD
    DROPPED IT.  Series-loss tank (node `x` between `L` and `R` is
    algebraic) with an asymmetric core, so the node-`v` PPV has a real
    mean (44% of its rms).  The exact reference is the reduced ODE
    `vdot = i_B(v) - i_L`, `L i_Ldot = v - R i_L` under the scipy adjoint:

        c_true = 1.204953e-07     <v_v> = +3.137167e-02

    With the full `G` in the consistent propagation `c` was 0.60 / 0.29 /
    0.15% low at 240/480/960 points and the mean 0.37 / 0.16 / 0.08% low
    -- first order, on BOTH the raw and consistent objects, which is what
    sent the review session looking for a linear-DAE cell (empty here:
    `ppv` is autonomous-only).  With the Schur complement `G[D,NZ] -
    G[D,Z] G[A,Z]^-1 G[A,NZ]`: 2.7e-4 / 7e-5 / 2e-5 and 8.6e-4 / 2e-4 /
    5e-5.  Gates at 240 and 480 points on both, and on the order.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    C_TRUE, MEAN_TRUE = 1.204953e-07, 3.137167e-02
    errs_c, errs_m = [], []
    for npts in (240, 480):
        cir = SubCircuit()
        cir.add_node('v')
        cir.add_node('x')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', 'x', L=1.0)
        cir['Rs'] = R('x', gnd, r=0.2)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: 1.0 * (u - u ** 3 / 3.0)
                           + 0.25 * (u ** 2 - 2.0))
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        pss = PSS(cir, method='gear', reltol=1e-12)
        m = cir.n - 1
        x0 = np.zeros(m)
        x0[0] = 2.0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=6.66, timestep=6.66 / npts, x0=x0,
                      maxiterations=200)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        errs_c.append(abs(pac.diffusion_constant(pss) / C_TRUE - 1.0))
        _v, info = pss.ppv()
        S = np.asarray(info['samples_eq'])[:, :m]
        h = np.diff(np.asarray(info['times'], dtype=float))
        mean = float((S[:len(h), 0] * h).sum()) / float(pss.period)
        errs_m.append(abs(mean / MEAN_TRUE - 1.0))
    assert errs_c[0] < 6e-4 and errs_c[1] < 2e-4, \
        'c is %.2e / %.2e from the exact 1.204953e-07; the full-G ' \
        'propagation gave 6.0e-3 / 2.9e-3' % tuple(errs_c)
    assert errs_c[0] / errs_c[1] > 3.0, \
        'c errors %.2e -> %.2e: not second order' % tuple(errs_c)
    assert errs_m[0] < 2e-3 and errs_m[1] < 5e-4 and errs_m[0] / errs_m[1] > 3.0, \
        'the node-v mean is %.2e / %.2e from +3.137167e-02, or not ' \
        'second order; the full-G propagation gave 3.7e-3 / 1.6e-3' \
        % tuple(errs_m)


def test_the_consistent_propagation_names_its_index_2_boundary():
    """⚠ `G[A,Z]` NONSINGULAR *IS* THE INDEX-1 CONDITION, so the Schur
    complement in the pair-consistent propagation does not exist at index
    2.  An L-I cutset (the tank inductor split through a node that sees
    only inductors) is an autonomous index-2 oscillator that `PSS` solves;
    `ppv` must run, warn ONCE with the reason, and return finite samples.
    ⚠ AND THE FALLBACK IS PRICED WHERE IT CAN BE SEEN.  The core is the
    bias-sensitive one (`vdp + 0.3 u^2`, rows NOT in quadrature -- the
    fixture on which the first-order PPV was visible at all) with its
    tank inductor split, so the ODE and its exact `c_true = 5.3703e-06`
    are unchanged and the topology is index 2.  Against the index-1
    consistent object (0.99825 / 0.99957 / 0.99992 at 400/800/1600) the
    fallback gives 0.99822 / 0.99956 / 0.99991: second order, below 1e-5
    and below the discretisation error at every grid.  The OTHER index-2
    topology, a C-V loop (the tank capacitance split between ground and a
    DC bias rail), gives 0.99853 / 0.99964 / 0.99993 -- also second order,
    within 3e-4 of the index-1 object and six times under its own error.
    So the guard is the whole answer at index 2, both topologies, and the
    projector-chain construction comes off the roadmap (the review
    session's partition, 2026-09-05; a split of plain van der Pol had been
    blind to it, its rows being in quadrature).  Boundary named by the
    review session from the pencil: `eig(-G_red, C[D,NZ])` equals the
    finite generalised eigenvalues of `(C, G)` to 1e-12 on the series-loss
    tank, and the reduction is undefined on `li_plus_rc` and `cv_plus_rc`.
    """
    import warnings
    from pycircuit.circuit.shooting import topological_index
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    for topology in ('L-I cutset', 'C-V loop'):
        cir = SubCircuit()
        cir.add_node('v')
        if topology == 'L-I cutset':
            cir.add_node('w')
            cir['C'] = C('v', gnd, c=1.0)
            cir['L1'] = L('v', 'w', L=0.5)
            cir['L2'] = L('w', gnd, L=0.5)
        else:
            ## the tank capacitance split between ground and a DC bias
            ## rail: the loop v-C-gnd-Vb-b-C1-v is capacitors closed by a
            ## voltage source, and AC-wise C || C1 = 1 leaves the ODE alone
            cir.add_node('b')
            cir['C'] = C('v', gnd, c=0.5)
            cir['C1'] = C('v', 'b', c=0.5)
            cir['Vb'] = VS('b', gnd, v=1.0)
            cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0)
                           + 0.3 * u ** 2)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        assert topological_index(cir)[0] == 2, topology
        pss = PSS(cir, method='gear', reltol=1e-12)
        m = cir.n - 1
        x0 = np.zeros(m)
        x0[0] = 2.0
        if topology == 'C-V loop':
            ## ⚠ seed the bias node AT its source: seeded at 0 V against a
            ## 1 V source the shooting Newton did not converge at any grid,
            ## and the cell read "empty" until the seed was fixed
            x0[[str(n) for n in cir.nodes][:m].index('b')] = 1.0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=6.731, timestep=6.731 / 400, x0=x0,
                      maxiterations=300)
        assert pss.converged, topology
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            _v, info = pss.ppv()
        mine = [x for x in w if 'index > 1' in str(x.message)]
        assert len(mine) == 1, \
            '%s: the consistent propagation must warn exactly once at ' \
            'index 2; got %d' % (topology, len(mine))
        assert np.all(np.isfinite(info['samples']))
        c = PAC(cir, toolkit=circuit.numeric).diffusion_constant(pss)
        assert abs(c / 5.3703e-06 - 1.0) < 3e-3, \
            '%s: the index-2 fallback gives c %.2e from the exact ' \
            '5.3703e-06 on the fixture that can see the dropped term; the ' \
            'index-1 object gives 1.8e-3 here' \
            % (topology, abs(c / 5.3703e-06 - 1.0))


def test_B16_the_oscillator_monodromy_comes_from_the_twin_default_radau():
    """⚠⚠ THE B16 DECISION, PINNED ON THE FIXTURE THAT SHOWED IT (2026-09-05).

    Bias-sensitive core, exact `Q_lambda = 5.9083` and `c_true = 5.3703e-06`
    (scipy adjoint of the ODE).  Trapezoidal's OWN monodromy is unusable
    with either opener -- `Q` 11.1 / 28.4 / 63.9 at 400/800/1600 points
    with the default (diverging under refinement), 3086 / 12228 / 48699
    with `x0_unknown=True` (a spurious multiplier at 1) -- while its state
    and period are second order.  So "the most accurate" is per quantity:
    the state keeps the method asked for, the monodromy comes from a TWIN
    on the same grid (`PSS.monodromy_twin`).

    ⚠ THE TWIN DEFAULTS TO RADAU (2026-09-24, Andreas: "Set radau as default
    twin"; TR-BDF2 from 2026-09-05), Gear-2 and TR-BDF2 selectable.  On THIS
    fixture the Radau twin reads `Q` 5.908399 (err 1.7e-5) and `c` err
    3.3e-5 -- both at the precision of the five-digit exact values -- where
    TR-BDF2 read 5.90845 (2.6e-5) and 1.3e-4, and Gear-2 5.90942 (1.9e-4)
    and 1.7e-3.  Gates: under trap the default `ppv` reports the Radau
    twin's `Q`/`c` (tighter than TR-BDF2's `c`, asserted); `monodromy = 'gear'`
    restores the Gear-2 twin (also asserted, so the option is live); the
    native path still shows trapezoidal's own defect; the state is
    untouched; and a one-step orbit too poor to seed the twin (Euler at 400
    points: period 5% off, amplitude 55% off) gets the error with the
    reason, not a number.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)

    def build():
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0)
                           + 0.3 * u ** 2)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        return cir

    def solve(method, maxiterations=300):
        cir = build()
        pss = PSS(cir, method=method, reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=6.731, timestep=6.731 / 400,
                      x0=np.array([2.0, 0.0]), maxiterations=maxiterations)
        assert pss.converged
        return cir, pss
    cir, pss = solve('trap')
    T_state = float(pss.period)
    _v, info = pss.ppv()
    c = PAC(cir).diffusion_constant(pss)
    ## the DEFAULT twin is Radau, and it reads Q/c CLOSER to the exact than
    ## the TR-BDF2 and Gear-2 twins do on this fixture -- the tolerances
    ## below are tight enough that the Gear-2 twin's numbers (err 1.9e-4 /
    ## 1.7e-3) and the TR-BDF2 twin's `c` (1.3e-4) would FAIL them, so they
    ## encode the improvement, not just the value.
    assert info['monodromy_method'] == 'radau'
    assert abs(info['Q'] / 5.9083 - 1.0) < 5e-5, \
        'trap+radau-twin reports Q = %.5f against exact 5.9083 (err %.1e); ' \
        'the Gear-2 twin gives 5.90942, err 1.9e-4' \
        % (info['Q'], abs(info['Q'] / 5.9083 - 1.0))
    assert abs(c / 5.3703e-06 - 1.0) < 1e-4, \
        'trap+radau-twin reports c %.2e from the true; the TR-BDF2 twin ' \
        'gives 1.3e-4, the Gear-2 twin 1.7e-3' % abs(c / 5.3703e-06 - 1.0)
    assert float(pss.period) == T_state, 'the state must not move'
    assert np.asarray(pss.waveform[1]).shape == \
        np.asarray(pss.monodromy_twin().waveform[1]).shape

    ## Gear-2 is one setting away, not retired: the same run under
    ## `monodromy = 'gear'` uses the Gear-2 twin and reports its (slightly
    ## less accurate) numbers, confirming the option is live.
    cir_g, pss_g = solve('trap')
    pss_g.monodromy = 'gear'
    _vg, info_g = pss_g.ppv()
    assert info_g['monodromy_method'] == 'gear'
    assert abs(info_g['Q'] / 5.9094 - 1.0) < 3e-4, \
        'gear twin should read 5.9094 here; got %.5f' % info_g['Q']
    ## the native path: the defect, pinned so it is not rediscovered
    pss.monodromy = 'native'
    _v2, info2 = pss.ppv()
    assert info2['monodromy_method'] == 'trap' and info2['Q'] > 10.0, \
        "trap's own second multiplier read Q = %.3f; it was 11.1 here" \
        % info2['Q']
    ## Euler at this grid (orbit 55% off) is too poor to seed a twin, and the
    ## point B16 pins is the INVARIANT: the twin is REFUSED with a reason, never
    ## a plausible-wrong Q.  ⚠ THE REFUSAL MECHANISM IS ROUNDOFF-SENSITIVE HERE
    ## and is deliberately NOT pinned: this seed sits on a spurious-orbit basin
    ## boundary, so a ~1e-14 change in the step (e.g. TR-BDF2 written in the
    ## generic stage-derivative form vs the old BDF2-companion form -- the same
    ## method to 14 digits) tips the free-period Newton between two refusals --
    ## CONVERGING to a spurious limit cycle that the orbit-consistency guard
    ## then catches (Q = 1.97 against the exact 5.91, 'spurious'), or NOT
    ## CONVERGING at all ('did not converge').  Both are loud; both refuse.  The
    ## Gear-2 twin refuses the same seed (by non-convergence).  Pinning one
    ## mechanism would pin a knife-edge; the assertion accepts either refusal.
    _refused = 'spurious|too poor to seed|did not converge'
    ## ⚠ THE BUDGET IS 20, NOT 300 (2026-09-23).  The twin used to inherit the
    ## caller's `maxiterations`, and from this seed both twins refuse by NOT
    ## converging -- so at 300 the two refusals below ran their whole budget,
    ## times the solve's retry ladder: 971 + 968 traversals, 790 of this
    ## test's 807 s and the gate's single longest pole.  At 20: 42.5 + 14.6 s.
    ## (The twin is now also capped at `PSS.TWIN_MAXITER` = 40 and warns when
    ## it hits the cap -- see the twin-cap test.)  The refusal is the pinned
    ## invariant, and it holds at any budget.
    _c3, p3 = solve('euler', maxiterations=20)
    with pytest.raises(RuntimeError, match=_refused):
        p3.ppv()
    p3.monodromy = 'gear'
    with pytest.raises(RuntimeError, match=_refused):
        p3.ppv()


def test_the_pnoise_excess_over_phase_only_is_the_amplitude_mode():
    """✅ A9's open question, closed by a POSITION test (2026-09-05).

    `pnoise` exceeds the phase-only prediction `P_c f0^2 c / df^2` by 16%
    at `df/f0 = 1e-2` on van der Pol at Q = 8, and the question was what
    the excess is.  The review session proposed testing its POSITION
    rather than its size: the amplitude mode's sideband is a Lorentzian of
    half-width `f0/(2 pi Q_lambda)` and the phase part is `1/df^2`, so
    their ratio `E(df)` is a STEP with its half-rise at that corner,
    moving as `1/Q_lambda`, saturating where AM equals PM (E = 1).

    Measured over Q = 4..32 (an 8x range): the corner moves as
    `Q_lambda^-1.03`; a Lorentzian step ALONE fits badly (rms 0.05-0.13,
    corner 1.4x off), a step PLUS a term linear in `df/f0` fits to rms
    0.003 with the corner at 1.02-1.12x the prediction (converging to 1
    with Q), `E_inf = 1.01`, and a Q-INDEPENDENT linear coefficient of
    1.74 / 1.81 / 1.83 / 1.84.  A term linear in `df` against a `1/df^2`
    part is a `1/df` piece of the spectrum.  It is NOT Traversa & Bonani's
    correlation term: that coefficient goes as `Q_lambda / c`, `c` is
    exactly flat in Q here (slope -0.000), so it would scale as
    `Q_lambda`; and it persists to `df = 0.3 f0`, fifteen corners out,
    where a cross term saturates.  ✅ IT IS THE TANK'S OWN FIRST-ORDER
    ASYMMETRY, found by the review session's parity test: on the LOWER
    sideband the coefficient FLIPS SIGN (upper +1.81 / +1.84, lower -2.20
    / -2.16 at Q = 8 / 32), so it is odd in `df` -- an asymmetry of the
    response, which a correction to the even phase-only reference could
    not produce.  Its odd part is 2.003 / 2.002, and 2 is what the tank
    gives: `|Z|^2 ~ w^2/(w^2 - w0^2)^2 = (1/4k^2)(1+k)^2/(1+k/2)^2 =
    (1/4k^2)(1 + k + ...)` at `w = w0 (1+k)`, and with the far-out total
    twice the phase part (`E_inf = 1`) the linear coefficient is `2 x 1`.
    Derived, not fitted.  The even remainder (-0.19) is NOT a term: under a
    window sweep the odd part stays at 2.00 in every window and model-free
    `(E+ - E-)/2k` reads 1.99-2.00, while the even coefficient drifts with
    the window (-0.19 -> -0.32) -- absorbed step curvature plus the tank's
    even -k^2/4.  Nothing here is open.

    Gated at Q = 8 and Q = 32 on the two-term fit: the corner within 20%
    of `f0/(2 pi Q_lambda)`, `E_inf` within 10% of 1, the corner ratio
    between the two Q values within 15% of the `1/Q_lambda` prediction,
    and at Q = 8 the lower sideband's linear coefficient of the opposite
    sign with the odd part within 5% of the tank's 2.
    """
    import warnings
    from scipy.optimize import least_squares
    circuit.default_toolkit = circuit.numeric
    ks = np.logspace(-3, -0.5, 14)
    out = {}
    for Q in (8.0, 32.0):
        mu = 1.0 / (2.0 * np.pi * Q)
        T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
        npts = 400 if Q < 16 else 800
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0))
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        pss = PSS(cir, method='gear', reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / npts, x0=np.array([2.0, 0.0]),
                      maxiterations=300)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        f0 = 1.0 / float(pss.period)
        ov = [str(n) for n in cir.nodes].index('v')
        c = float(pac.diffusion_constant(pss))
        Ql = float(pss.ppv()[1]['Q'])
        X = np.asarray(pss.waveform[1], dtype=float)[ov][:-1]
        A1 = 2.0 * abs(np.fft.rfft(X)[1]) / len(X)
        Pc = 0.5 * A1 * A1
        E = np.array([float(np.real(pac.pnoise(pss, f0 * (1.0 + k), ov,
                                               sweeptype='absolute')[0]))
                      / Pc * k * k / c - 1.0 for k in ks])
        pred = 1.0 / (2.0 * np.pi * Ql)
        fit = least_squares(
            lambda q: q[0] * ks ** 2 / (ks ** 2 + q[1] ** 2) + q[2] * ks - E,
            x0=[1.0, pred, 0.5], bounds=([0, 1e-5, -10], [10, 1, 10]))
        Einf, kc, b = fit.x
        if Q == 8.0:
            ## the parity test: the LOWER sideband
            El = np.array([float(np.real(pac.pnoise(pss, f0 * (1.0 - k),
                                                    ov, sweeptype='absolute')[0]))
                           / Pc * k * k / c - 1.0 for k in ks])
            fl = least_squares(
                lambda q: q[0] * ks ** 2 / (ks ** 2 + q[1] ** 2)
                + q[2] * ks - El,
                x0=[1.0, pred, -0.5], bounds=([0, 1e-5, -10], [10, 1, 10]))
            bl = fl.x[2]
            assert bl < 0.0 < b, \
                'the linear term must be ODD in df: upper %+.2f, lower ' \
                '%+.2f' % (b, bl)
            odd = 0.5 * (b - bl)
            assert abs(odd / 2.0 - 1.0) < 0.05, \
                "the odd part is %.3f; the tank's w^2/(w^2-w0^2)^2 gives 2" \
                % odd
        rms = float(np.sqrt(np.mean(fit.fun ** 2)))
        assert rms < 0.02, 'Q=%g: the step-plus-linear form misfits E by ' \
            'rms %.3f' % (Q, rms)
        assert abs(kc / pred - 1.0) < 0.2, \
            'Q=%g: the corner sits at %.2fx f0/(2 pi Q_lambda); the ' \
            'excess is not the amplitude mode' % (Q, kc / pred)
        assert abs(Einf - 1.0) < 0.1, \
            'Q=%g: the excess saturates at %.2f, not at AM = PM' % (Q, Einf)
        out[Q] = (Ql, kc, b)
    r = (out[8.0][1] / out[32.0][1]) / (out[32.0][0] / out[8.0][0])
    assert abs(r - 1.0) < 0.15, \
        'the corner moved by %.2fx the 1/Q_lambda prediction between Q = 8 ' \
        'and Q = 32' % r
    assert abs(out[8.0][2] / out[32.0][2] - 1.0) < 0.1, \
        'the linear term is not Q-independent: %.2f vs %.2f' \
        % (out[8.0][2], out[32.0][2])


def test_the_kTC_gate_rejects_an_unscaled_CY():
    """MUTATION CHECK (P3): inject the `CY` vs `CY/2` defect that a Monte
    Carlo once confirmed rather than caught, and assert the gate fires.

    `diffusion_constant` contracts `CY/2` (one-sided to two-sided).  With
    `_cy_reduced` monkeypatched to return twice its value -- the exact
    historical bug -- `c` doubles, so a gate pinned near a reference value
    now sees 2x and must reject it.  A gate that still passed under this
    mutation would be vacuous.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=2 * np.pi, timestep=2 * np.pi / 400,
                  x0=np.array([2.0, 0.0]), maxiterations=300)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    c_ok = pac.diffusion_constant(pss)
    orig = pac._cy_reduced
    pac._cy_reduced = lambda pss_, w: 2.0 * np.asarray(orig(pss_, w))
    try:
        c_bug = pac.diffusion_constant(pss)
    finally:
        pac._cy_reduced = orig
    assert abs(c_bug / c_ok - 2.0) < 1e-9, \
        'the injected CY-vs-CY/2 defect did not double c (%.3e vs %.3e); ' \
        'the functional is not reading _cy_reduced as the gate assumes' \
        % (c_bug, c_ok)


@pytest.mark.slow
def test_the_diffusion_constants_numerical_floor_is_the_grid_not_the_tolerance():
    """How small a `c` can this stack compute before its own error dominates?

    ⚠⚠ THE PUBLISHED GATE DOES NOT TRANSFER, AND CHECKING THAT FIRST IS THE
    POINT.  Biggio et al. measure a simulator's numerical noise floor by FFT-ing
    a NOISELESS oscillator and looking between the harmonics -- whatever is
    there is the floor.  That assumes SPECTRAL ESTIMATION.  This stack is
    CLOSED FORM: `oscillator_spectrum` returns a Lorentzian scaled from
    `c = (1/T) int v^T B B^T v dt`, so a noiseless circuit has `B = 0`, `c = 0`
    and `L = -inf`.  There is no broadened spectrum to measure.

    ⚠⚠⚠ AND THE OBVIOUS REPLACEMENT IS A GATE THAT CANNOT FAIL.  Sweeping the
    source PSD and checking `c` tracks it linearly gives `c/psd` constant to
    **1.7e-16 over 28 decades** -- because `CY ~ psd` factors straight out of
    the quadratic form.  That is a STRUCTURAL IDENTITY confirming arithmetic,
    the same family as a zero-vs-zero pass.  **A measurement whose outcome is
    fixed by the algebra says nothing about the implementation.**

    The numerical error lives in the PPV and the orbit, so the knobs are the
    GRID and the TOLERANCE.  Measured on the van der Pol noise fixture::

        (1) grid, at reltol 1e-12          (2) tolerance, at npts 240
        npts   c               rel chg     reltol   c
          60   7.987354e-08    --          1e-08    8.042025661140e-08
         120   8.030800e-08    5.41e-03    1e-10    8.042025661200e-08
         240   8.042026e-08    1.40e-03    1e-12    8.042025661266e-08
         480   8.044852e-08    3.51e-04    1e-14    8.042025661208e-08
         960   8.045561e-08    8.81e-05

    **The grid change falls 4x per doubling -- O(h^2) -- and the tolerance does
    not move `c` at all past ten digits (spread ~1.5e-11).**  So the floor is
    DISCRETISATION, not the Newton tolerance, and tightening `reltol` to buy
    phase-noise accuracy buys nothing: refine the grid instead.

    At 240 points the uncertainty is ~4e-04 relative, i.e. **~0.0004 dB** on a
    reported phase noise -- far below anything that would corrupt a result.  So
    the concern is real for a spectral-estimation simulator and STRUCTURALLY
    ABSENT here; the closed-form route buys that.

    ⚠⚠ SCOPE, AND IT WAS NARROWER THAN THIS TEST CLAIMED.  `mu = 1` is not
    "moderate Q" -- van der Pol's amplitude relaxes at rate `mu`, so
    `|lambda_2| = exp(-2 pi mu)` and `Q = 1/(2 mu)`: THIS FIXTURE IS Q ~ 0.5.
    The high-Q measurement it deferred is now
    `test_the_diffusion_constant_at_high_q_has_an_analytic_reference`, and it
    changed the conclusion: the grid floor grows LINEARLY IN Q (1.8e-05 Q for
    gear), so the concern was real -- but it is a METHOD problem, and `radau`
    is six orders lower at the same grid.

    ⚠ AND THE `3.0 < ratio < 5.0` BELOW IS A `mu = 1` STATEMENT, not a
    property of `c`.  For every Q >= 5 the ratio is 8 (O(h^3)), because on an
    autonomous problem the O(h^2) error is a FREQUENCY error and the solve
    absorbs it into `T` -- measured there.  If this fixture's `mu` ever moves,
    this assertion moves with it.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    psd = 1e-6

    def vdp():
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
        c['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
        return c

    def cval(npts, reltol):
        cir = vdp()
        p = PSS(cir, method='gear', reltol=reltol)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=6.6634, timestep=6.6634 / npts,
                    x0=np.array([2.0, 0.0]), maxiterations=60)
        assert p.converged
        return float(PAC(cir, toolkit=circuit.numeric).diffusion_constant(p))

    ## (1) the grid is the knob that moves it, and it converges at O(h^2)
    cs = [cval(n, 1e-12) for n in (60, 120, 240, 480)]
    chg = [abs(cs[i + 1] - cs[i]) / abs(cs[i + 1]) for i in range(len(cs) - 1)]
    assert chg[0] > chg[1] > chg[2], 'c is not converging with the grid: %r' % chg
    ratio = chg[1] / chg[2]
    assert 3.0 < ratio < 5.0, \
        'the grid error should fall ~4x per doubling (O(h^2)), got %.2f' % ratio

    ## (2) ⚠ the TOLERANCE does not move it -- so `reltol` is the wrong dial for
    ## phase-noise accuracy, and a caller tightening it is paying for nothing.
    ct = [cval(240, rt) for rt in (1e-8, 1e-10, 1e-12, 1e-14)]
    spread = (max(ct) - min(ct)) / abs(np.mean(ct))
    assert spread < 1e-8, \
        'reltol moved c by %.3e -- if this ever becomes the limit, the floor ' \
        'story above changes and the docstring must be re-measured' % spread

    ## (3) ⚠ AND THE STRUCTURAL IDENTITY, ASSERTED SO IT IS NOT MISTAKEN FOR A
    ## GATE: c is exactly linear in the PSD, so sweeping it proves nothing.
    cir_a, cir_b = vdp(), vdp()
    cir_b['n'].ipar.noisePSD = psd * 1e-12
    cir_b.update_iparv()
    got = []
    for cc in (cir_a, cir_b):
        p = PSS(cc, method='gear', reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=6.6634, timestep=6.6634 / 240,
                    x0=np.array([2.0, 0.0]), maxiterations=60)
        got.append(float(PAC(cc, toolkit=circuit.numeric).diffusion_constant(p))
                   / float(cc['n'].ipar.noisePSD))
    assert abs(got[0] - got[1]) / abs(got[0]) < 1e-12, \
        'c/psd should be constant BY CONSTRUCTION -- if it is not, the ' \
        'quadratic form has acquired a psd dependence it should not have'


@pytest.mark.slow
def test_grid_error_measures_the_discretisation_floor_and_refuses_when_it_cannot():
    """`PSS.grid_error` VALIDATED against the analytic high-Q reference.

    The floor of this stack is discretisation, it grows linearly in `Q`, and
    it is a METHOD property -- the sibling test measures `~1.8e-05 Q` for gear
    against `~7.0e-12 Q` for radau. Shipping a PREDICTOR from those constants
    would extrapolate a fit across a regime change (gear is `O(h^3)` here and
    `O(h^2)` at `mu = 1`), so `grid_error` REFINES THE ACTUAL CIRCUIT instead.
    This test checks the instrument before anyone trusts it.

    Measured at `mu = 0.005` (`Q = 100`), 120 points refined 2x twice, against
    `c = psd/16 (1 + (11/32) mu^2)`::

        method   observed order   power_law   est rel err   true rel err
        gear         3.04           True       2.208e-04     2.257e-04
        trap         3.02           True       3.29e-06      3.320e-06
        (trap's row read `6.45 False (withheld) 4.061e-06` until the
        2026-09-14 fix -- a period-normalisation defect; see step 4)
        radau        5.03           True       2.132e-11     1.076e-11

    `gear` recovers its `O(h^3)` autonomous rate and its estimate lands within
    **2%** of the true error. `radau` recovers order 5 and over-states by 2x,
    which is the safe direction.

    ⚠⚠ `trap` IS THE REASON THE VALIDITY CHECK EXISTS. Its error changes sign
    near `Q = 100`, two terms nearly cancel, and consecutive differences then
    shrink FASTER than the error: apparent order 6.45, estimate 300x too
    small, with monotone same-signed deltas and nothing else suspicious. ⚠ A
    generic `0.5 <= order <= 8` range ACCEPTS it -- that was the first version
    of this check -- and a sign test does not catch it either (both deltas are
    negative). Only the CEILING at the method's own order rejects it, because
    a method cannot converge faster than its order.

    ⚠ And the ceiling needs the `+1.5` allowance: on an autonomous problem the
    period absorbs the leading frequency error, so `gear` (nominal 2) really
    does converge at 3.01. Without it the check would reject the shipped
    default on its own reference fixture.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    psd, T0, mu = 1e-6, 2.0 * np.pi, 0.005
    ## The `O(mu^2)` term is part of the PHYSICS, not an error: leaving it out
    ## makes radau's true error read 8.594e-06 at EVERY grid -- a constant,
    ## which is the tell that the reference and not the method is being
    ## measured. It cost one wrong reading here before it was noticed.
    analytic = psd / 16.0 * (1.0 + (11.0 / 32.0) * mu ** 2)

    def vdp():
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: mu * (u - u ** 3 / 3.0))
        c['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
        return c

    def measure(method):
        cir = vdp()
        p = PSS(cir, method=method, reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T0, timestep=T0 / 120,
                    x0=np.array([2.0, 0.0]), maxiterations=80)
            r = p.grid_error(
                lambda q: float(PAC(cir, toolkit=circuit.numeric)
                                .diffusion_constant(q)))
        true = abs(r['values'][-1] - analytic) / analytic
        return r, true

    ## 1. gear: the order is its own, and the estimate is the true error.
    r, true = measure('gear')
    assert r['power_law'], 'gear should follow a power law, order %r' % (
        r['order'],)
    assert abs(r['order'] - 3.0) < 0.3, \
        'gear order is %.2f, the autonomous O(h^3) rate is ~3' % r['order']
    assert 0.5 < r['rel_error'] / true < 2.0, \
        'gear estimate %.3e against a true error %.3e' % (
            r['rel_error'], true)

    ## 2. radau: order 5, and it must not UNDER-state.
    r_r, true_r = measure('radau')
    assert r_r['power_law'], 'radau should follow a power law'
    assert abs(r_r['order'] - 5.0) < 0.3, \
        'radau order is %.2f, want ~5' % r_r['order']
    assert r_r['rel_error'] > 0.5 * true_r, \
        'radau estimate %.3e under-states the true error %.3e' % (
            r_r['rel_error'], true_r)

    ## 3. And the six-order method gap the whole item rests on.
    assert r_r['rel_error'] < 1e-5 * r['rel_error'], \
        'radau %.3e is not far below gear %.3e' % (
            r_r['rel_error'], r['rel_error'])

    ## 4. trap: its `c` is its twin's (Radau since 2026-09-24; TR-BDF2
    ##    before, order 3.02), so it is ESTIMABLE, at the twin's order -- and
    ##    the ceiling is the higher of trap's and the twin's (5.03 here).
    ##    ⚠⚠ This step used to assert the opposite -- that `grid_error` must
    ##    REFUSE trap here, its error "changing sign near Q = 100" and a
    ##    two-grid estimate under-stating it 300x.  That sign change was a
    ##    DEFECT: `diffusion_constant` divided the twin's integral by trap's
    ##    own period (see
    ##    `test_ppv_quadratures_normalise_by_the_period_of_the_orbit_they_integrate`).
    ##    Fixed, trap reads order 3.02 and a clean estimate.
    cir = vdp()
    p = PSS(cir, method='trap', reltol=1e-12)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=T0, timestep=T0 / 120,
                x0=np.array([2.0, 0.0]), maxiterations=80)
        r_t = p.grid_error(
            lambda q: float(PAC(cir, toolkit=circuit.numeric)
                            .diffusion_constant(q)))
    true_t = abs(r_t['values'][-1] - analytic) / analytic
    assert r_t['power_law'] and abs(r_t['order'] - 5.0) < 0.3, \
        'trap (via its radau twin) should read radau\'s order 5, got ' \
        '%r' % (r_t['order'],)
    assert 0.5 < r_t['rel_error'] / true_t < 2.0, \
        'trap estimate %.3e against a true error %.3e' % (
            r_t['rel_error'], true_t)

    ## 4b. THE CEILING STILL NEEDS ITS COUNTEREXAMPLE, and the defect above
    ##     is a real one: a quantity MIXING TWO DISCRETISATIONS -- the twin's
    ##     `c` times `T_twin / T_trap`, an `O(h^3)` error plus an `O(h^2)`
    ##     one of opposite sign.  Built from real solves, it reproduces the
    ##     old apparent order 6.45 exactly, with monotone same-signed deltas;
    ##     a plain `0.5 <= order <= 8` range ACCEPTS it and only the ceiling
    ##     at the method's order refuses.  ⚠ UNDER THE TR-BDF2 TWIN: under
    ##     the radau twin (the default since 2026-09-24) the twin's `c` is
    ##     O(h^5), the O(h^2) term alone is left, and the quantity is an
    ##     honest order 2.01 -- no cancellation to refuse.
    p.monodromy = 'trbdf2'

    def mixed(q):
        c = float(PAC(cir, toolkit=circuit.numeric).diffusion_constant(q))
        return c * float(q.monodromy_twin().period) / float(q.period)
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter('always')
        r_m = p.grid_error(mixed)
    assert not r_m['power_law'], \
        'the mixed quantity (apparent order %r) was ACCEPTED; its estimate ' \
        'under-states the true error by ~300x' % (r_m['order'],)
    assert r_m['order'] is not None and r_m['order'] > 3.5 + 2.0, \
        'the counterexample must be the CEILING case, order %r' % (
            r_m['order'],)
    assert any('single power law' in str(w.message) for w in caught), \
        'the refusal must warn; got %r' % [str(w.message) for w in caught]

    ## 5. An explicit `grid` cannot be refined by a timestep, and reporting
    ##    0.0 there would be a confident lie.
    cir = vdp()
    p = PSS(cir, method='radau', reltol=1e-12)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=T0, timestep=T0 / 120, x0=np.array([2.0, 0.0]),
                maxiterations=80,
                grid=np.full(120, 1.0 / 120.0))   ## step FRACTIONS, sum 1
    try:
        p.grid_error(lambda q: 1.0)
    except ValueError as exc:
        assert 'grid' in str(exc)
    else:
        raise AssertionError('grid_error accepted an explicit grid')


def test_ppv_quadratures_normalise_by_the_period_of_the_orbit_they_integrate():
    """A PPV quadrature divides by the period of THE ORBIT ITS SAMPLES LIVE ON.

    An autonomous `trap` run reads its PPV from a TR-BDF2 twin
    (`monodromy_twin`), re-converged on the same grid, whose period differs
    from trap's by `O(h^2)`.  `diffusion_constant`, `colour_projection` and
    `coloured_diffusion_resolved` integrated the TWIN's samples over the
    twin's steps and divided by `pss.period` -- TRAP's.  So trap's `c` was
    `c_twin * T_twin / T_trap`: an `O(h^3)` positive error plus an `O(h^2)`
    period mismatch of the opposite sign.

    ⚠⚠ THAT SUM IS THE "SIGN CHANGE" OF E3 AND OF THE RADAU-DEFAULT RECORD.
    Measured at `Q = 100` against the analytic reference, signed:
    `+9.71e-05, -2.91e-06, -4.06e-06, -1.43e-06, -4.08e-07` at 120..1920
    points, reproduced from `e_trbdf2 - dT/T` to three digits at every row
    (the rows were PREDICTED before they ran), while trbdf2's own `c` falls
    monotonically at order 3.  The non-monotone row and the 300x
    under-statement `grid_error` refused were this defect, not a property of
    trapezoidal.  Every amplitude check passed, because the error is `O(h^2)`
    and SHRINKS -- only a comparison against the twin's own value sees it.

    Gate: under `trap` each quantity equals the twin's own, and the Parseval
    identity `coloured_diffusion_resolved == diffusion_constant` (white
    source) closes to round-off; before the fix they missed by exactly
    `T_twin/T_trap - 1 = -2.96e-05`, and Parseval by 1.5e-09 because the
    harmonic frequency came from the other orbit too.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    mu = 0.05

    def osc():
        ## asymmetric (0.3 u^2) so `vbar` is not a symmetry zero
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: mu * (u - u ** 3 / 3.0)
                         + 0.3 * mu * u ** 2)
        c['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        return c

    cir = osc()
    p = PSS(cir, method='trap', reltol=1e-12)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=2 * np.pi, timestep=2 * np.pi / 240,
                x0=np.array([2.0, 0.0]), maxiterations=80)
    tw = p.monodromy_twin()
    assert tw is not p, 'the fixture must exercise the twin'
    mismatch = float(tw.period) / float(p.period) - 1.0
    assert abs(mismatch) > 1e-6, \
        'the periods must differ or this gate cannot fail (%.3e)' % mismatch

    pac = PAC(cir, toolkit=circuit.numeric)
    c_host, c_twin = pac.diffusion_constant(p), pac.diffusion_constant(tw)
    assert abs(c_host / c_twin - 1.0) < 1e-12, \
        'diffusion_constant under trap is %+.3e off its twin (the period ' \
        'mismatch is %+.3e)' % (c_host / c_twin - 1.0, mismatch)

    vb_h, i_h = pac.colour_projection(p)
    vb_t, i_t = pac.colour_projection(tw)
    assert np.max(np.abs(vb_h / vb_t - 1.0)) < 1e-12, \
        'colour_projection vbar off its twin by %s' % (vb_h / vb_t - 1.0)
    assert np.max(np.abs(i_h['rms'] / i_t['rms'] - 1.0)) < 1e-12

    cr_h = pac.coloured_diffusion_resolved(p, [1e-3],
                                           frequency_aware=False)[0]  # the DC fold's identity (2026-09-26: frequency-aware is the default)
    assert abs(cr_h / c_host - 1.0) < 1e-12, \
        'Parseval under trap: resolved %.12e against c %.12e' % (cr_h, c_host)


def test_the_diffusion_constant_at_high_q_has_an_analytic_reference():
    """The floor AT HIGH Q -- the regime the original concern actually named.

    The sibling test above measures at van der Pol `mu = 1` and records high Q
    as UNTESTED.  This is that measurement, and it changes three things.

    ⚠⚠ FIRST, `mu = 1` IS NOT MODERATE Q -- IT IS Q ~ 0.5.  van der Pol's
    amplitude obeys `A' = (mu/2)(A - A^3/4)`, so linearising at `A = 2` gives a
    relaxation rate `mu` and

        |lambda_2| = exp(-mu T) ~ exp(-2 pi mu),   Q = pi / (-ln|lambda_2|)
                                                     = 1 / (2 mu).

    THAT PREDICTION IS CHECKED HERE BEFORE ANYTHING RESTS ON IT, because a
    "high Q" fixture that is not high Q would make everything below vacuous.
    Measured `|lambda_2|` against `exp(-2 pi mu)`: 0.533079/0.533488 at
    mu = 0.1, 0.881910/0.881911 at 0.02, 0.969074/0.969072 at 0.005 -- exact
    where the small-mu theory holds, and visibly WRONG at mu = 1
    (0.000859 against a predicted 0.001867), which is the right behaviour for
    an asymptotic prediction and is why mu = 1 cannot be read as high Q.
    So `mu = 0.005` is **Q = 100**.

    ⚠⚠ SECOND, AT HIGH Q THERE IS AN ANALYTIC ANSWER, so this stops being a
    self-comparison.  As `mu -> 0` the circuit is a harmonic oscillator with
    `x = 2 cos t`.  With `v = A cos(theta)` and `w = v' = -A sin(theta)`,
    `theta = atan2(-w, v)`, so a perturbation of `v` alone moves the phase by
    `dtheta/dv = w/A^2 = -sin(theta)/A`.  The PPV's v-component is therefore
    `v1 = -sin(t)/2`, and for a white current source of density `psd` into a
    1 F capacitor::

        c = <v1^2> psd     = psd/8      [two-sided CY]
                           = psd/16     [this stack's CY/2 convention]
                           = 6.25e-08   at psd = 1e-6.

    ⚠ The 1/16 rather than 1/8 IS the `CY/2` convention `diffusion_constant`
    records, so this doubles as a pin on it.

    **And the approach is O(mu^2), which is what makes the limit usable as a
    reference rather than a hope.** Measured excess over `psd/16`, radau at 480
    points::

        mu       c                excess      excess/mu^2
        0.04     6.253436656e-08   5.4986e-04   0.34366
        0.02     6.250859322e-08   1.3749e-04   0.34373
        0.01     6.250214841e-08   3.4374e-05   0.34374
        0.005    6.250053711e-08   8.5937e-06   0.343748

    A clean power law -- the ratio is 4.00 for every halving -- converging on
    `11/32 = 0.34375`.  So at Q = 100 the PHYSICS is known to five digits and
    any deviation beyond it is NUMERICS.

    ⚠ THE NEXT TERM IS MEASURED TOO, because without it this reference runs
    out before radau does.  Residual of radau at 960 points against
    `psd/16 (1 + (11/32) mu^2)`, divided by `mu^4`::

        mu       offset rel     offset/mu^4
        0.020    -8.435e-09       -0.0527
        0.010    -5.259e-10       -0.0526
        0.005    -3.172e-11       -0.0508

    Constant over a 4x sweep, so the reference extends to

        c = psd/16 (1 + (11/32) mu^2 - 0.0527 mu^4).

    ⚠⚠ THIS MATTERS FOR READING ANY HIGH-ORDER RESULT AGAINST IT.  With only
    the `mu^2` term, radau's apparent error at `mu = 0.005` reads 3.2e-11 at
    EVERY grid -- a constant, which looks like a solver floor and is not; it
    is the REFERENCE's own truncation.  `grid_error` was briefly judged to
    under-state on exactly that reading (see
    `test_grid_error_measures_the_discretisation_floor_and_refuses_when_it_cannot`).
    ⚠ A constant "error" across a grid sweep means the REFERENCE, not the
    method -- the same tell as the two failed order sweeps in the Radau
    index-2 record.  That separation is the whole point.

    ⚠⚠ THIRD, AND THE ACTIONABLE PART: THE CONCERN DOES MATERIALISE -- THE
    FLOOR GROWS LINEARLY IN Q -- AND IT IS A METHOD PROBLEM, NOT A GRID OR
    TOLERANCE ONE.  Grid uncertainty in `c` at 240 points per period
    (240-against-960), and what it is worth on a reported phase noise::

        Q      gear        trap        radau       gear in dB
         100   1.79e-03    2.63e-05    6.97e-10    0.0078
         500   9.02e-03    1.32e-04    3.48e-09    0.0392
        1000   1.82e-02    2.63e-04    6.97e-09    0.0791

    **`gear` and `radau` both scale LINEARLY IN Q** (ratios 5.04/2.02 and
    5.00/2.00 against Q ratios 5 and 2) -- so `c`'s uncertainty is
    `~1.8e-05 Q` for gear and `~7.0e-12 Q` for radau, SIX ORDERS apart.  `trap`'s
    column did not fit a clean law because it was a DEFECT (the twin's `c`
    over trap's own period, fixed 2026-09-14; it read 1.49e-06 / 1.04e-04 /
    2.36e-04); corrected it is its TR-BDF2 twin's -- linear in Q too.

    So at Q = 1000 the shipped `gear` costs 0.08 dB and by Q = 10000 it would
    cost roughly 0.7 dB -- the concern was real -- while `radau` is at 7e-08
    there, i.e. nothing.  **REFINING THE GRID IS THE EXPENSIVE ANSWER AND
    CHANGING METHOD IS THE FREE ONE.**  `reltol` remains no answer at all: 1e-8
    to 1e-14 moves `c` by 1.3e-12 at Q = 100, the same non-answer as at mu = 1.

    ⚠⚠ AND THE SIBLING TEST'S `3.0 < ratio < 5.0` IS Q-SPECIFIC, WHICH NOTHING
    SAID.  `gear` converges at O(h^2) at mu = 1 and at **O(h^3)** for every
    Q >= 5 (ratio 8.01 at the finest grids, over five doublings).  A 2nd-order
    method giving 3rd-order answers wants a mechanism, and the one measured
    here is that ON AN AUTONOMOUS PROBLEM THE PERIOD IS AN UNKNOWN: at Q = 100
    the solved PERIOD converges at order **2.00** while the waveform and `c`
    converge at **3.01**.  The O(h^2) term is a FREQUENCY error, and the
    autonomous solve absorbs it into `T` instead of leaving it in the state.
    At mu = 1 the orbit is far from harmonic, the h^2 error has a genuine
    waveform component, and `c` is order 2 again.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    psd = 1e-6
    T0 = 2.0 * np.pi
    analytic = psd / 16.0

    def vdp(mu):
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u, _m=mu: _m * (u - u ** 3 / 3.0))
        c['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
        return c

    def run(mu, npts, method='gear', reltol=1e-12):
        cir = vdp(mu)
        p = PSS(cir, method=method, reltol=reltol)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res = p.solve(period=T0, timestep=T0 / npts,
                          x0=np.array([2.0, 0.0]), maxiterations=80)
        assert p.converged, '%s mu=%g npts=%d did not converge' % (method, mu, npts)
        cval = float(PAC(cir, toolkit=circuit.numeric).diffusion_constant(p))
        per = float(np.asarray(res['period']).ravel()[0]) if 'period' in res \
            else float(getattr(p, 'period', np.nan))
        return cval, per, p

    ## (1) ⚠ THE FIXTURE IS HIGH Q, CHECKED AND NOT ASSUMED.  Without this the
    ## rest is a measurement of some other regime.
    c_g = {}
    per_g = {}
    for n in (120, 240, 480, 960):
        c_g[n], per_g[n], p_last = run(0.005, n)
    fp = p_last.factored_period()
    lam = np.sort(np.abs(np.linalg.eigvals(np.column_stack(
        [np.asarray(fp.matvec(e), float).ravel()
         for e in np.eye(fp.width)]))))[::-1]
    l2 = float(lam[1])
    pred = float(np.exp(-2.0 * np.pi * 0.005))
    assert abs(l2 - pred) < 1e-4, \
        'the Q knob is not doing what this test claims: |lambda_2| = %.6f ' \
        'against the predicted exp(-2 pi mu) = %.6f' % (l2, pred)
    Q = float(np.pi / (-np.log(l2)))
    assert Q > 90.0, 'Q = %.1f is not the high-Q regime this test is about' % Q

    ## (2) THE ANALYTIC LIMIT, and the mu^2 law that makes it usable.  radau,
    ## whose own grid error here is six orders below the physics.
    exc = {}
    for mu in (0.02, 0.01, 0.005):
        cv, _per, _p = run(mu, 480, 'radau')
        exc[mu] = (cv - analytic) / analytic
        assert exc[mu] > 0, \
            'mu=%g: c sits BELOW the harmonic limit (%.4e), which the ' \
            'amplitude correction cannot do' % (mu, exc[mu])
    for a, b in ((0.02, 0.01), (0.01, 0.005)):
        r = exc[a] / exc[b]
        assert abs(r - 4.0) < 0.05, \
            'the excess over psd/16 should fall 4x per halving of mu (an ' \
            'O(mu^2) amplitude correction); mu %g -> %g gave %.3f. If this ' \
            'is not 4 the analytic reference is wrong and every number ' \
            'below rests on nothing.' % (a, b, r)
    k = exc[0.005] / 0.005 ** 2
    assert abs(k - 0.34375) < 2e-3, \
        'the measured coefficient is %.5f, not the 11/32 on record -- the ' \
        'limit or the convention has moved' % k

    ## (3) THE SEPARATION THE ANALYTIC LIMIT BUYS: how much of each method's
    ## deviation is PHYSICS (the mu^2 term) and how much is GRID.
    phys = k * 0.005 ** 2
    c_radau, _p, _o = run(0.005, 240, 'radau')
    dev_radau = abs((c_radau - analytic) / analytic - phys) / phys
    dev_gear = abs((c_g[240] - analytic) / analytic - phys) / phys
    assert dev_radau < 1e-3, \
        "radau's deviation from the analytic limit should be the physical " \
        'mu^2 term and nothing else; it is off by %.3e of it' % dev_radau
    assert dev_gear > 100.0, \
        "gear's 240-point deviation should be dominated by the GRID (it is " \
        '%.1f times the physical term). If it is not, this fixture no longer ' \
        'separates the two and (4) below is vacuous.' % dev_gear

    ## (4) THE FLOOR ITSELF -- still negligible at Q = 100, and BOUNDED so a
    ## regression would show.  ~1.8e-3 relative is ~0.008 dB.
    floor = abs(c_g[240] - c_g[960]) / abs(c_g[960])
    assert 1e-4 < floor < 5e-3, \
        "gear's 240-point grid uncertainty at Q = 100 is %.3e; the record " \
        'says 1.8e-03 (about 0.008 dB)' % floor
    spread = dev_gear / max(dev_radau, 1e-30)
    assert spread > 1e4, \
        'the whole high-Q finding is that the METHOD dominates: radau should ' \
        'beat gear by orders here, and the ratio is only %.3g' % spread

    ## (5) `reltol` IS STILL THE WRONG DIAL, checked in the regime the concern
    ## named rather than only where it was convenient.
    ct = [run(0.005, 240, 'gear', rt)[0] for rt in (1e-8, 1e-14)]
    assert abs(ct[0] - ct[1]) / abs(ct[0]) < 1e-9, \
        'reltol moved c by %.3e at Q = 100 -- if tolerance ever becomes the ' \
        'limit, the "refine the grid" advice changes' % (
            abs(ct[0] - ct[1]) / abs(ct[0]))

    ## (6a) ⚠ THE FLOOR GROWS LINEARLY IN Q, AND THE METHOD SETS THE RATE.
    ## This is the part that says the original concern was real -- and that
    ## the answer is `radau`, not a finer grid.
    hi = {}
    for meth in ('gear', 'radau'):
        a, _p, _o = run(0.0005, 240, meth)     # Q = 1000
        b, _p, _o = run(0.0005, 960, meth)
        hi[meth] = abs(a - b) / abs(b)
    lo_gear = floor                            # Q = 100, from (4)
    c_r240, _p, _o = run(0.005, 240, 'radau')
    c_r960, _p, _o = run(0.005, 960, 'radau')
    lo_radau = abs(c_r240 - c_r960) / abs(c_r960)
    for meth, lo in (('gear', lo_gear), ('radau', lo_radau)):
        r = hi[meth] / lo
        assert 8.0 < r < 12.0, \
            '%s: the grid floor should grow LINEARLY in Q (10x from Q=100 to ' \
            'Q=1000) and grew %.2fx. The recorded rates are ~1.8e-05 Q for ' \
            'gear and ~7.0e-12 Q for radau.' % (meth, r)
    assert hi['gear'] / hi['radau'] > 1e5, \
        'the actionable finding is that the METHOD sets the rate: at Q = 1000 ' \
        'gear should be ~6 orders worse than radau and is only %.3g times' \
        % (hi['gear'] / hi['radau'])
    assert hi['gear'] * 4.3429 < 0.5, \
        'gear at Q = 1000 and 240 points is worth %.4f dB; the record says ' \
        '0.079 dB' % (hi['gear'] * 4.3429)

    ## (6b) ⚠ THE ORDER, AND ITS MECHANISM.  `c` is O(h^3) here where the
    ## sibling test asserts O(h^2) at mu = 1 -- because the O(h^2) term is a
    ## FREQUENCY error and an autonomous solve absorbs it into `T`.
    def order(v):
        d = [abs(v[i + 1] - v[i]) for i in range(len(v) - 1)]
        return [float(np.log2(d[i] / d[i + 1])) for i in range(len(d) - 1)]
    oc = order([c_g[n] for n in (120, 240, 480, 960)])
    oT = order([per_g[n] for n in (120, 240, 480, 960)])
    assert 2.8 < oc[-1] < 3.3, \
        'c converges at order %.2f at Q = 100; the record says 3.01, and the ' \
        "sibling test's 3.0 < ratio < 5.0 is a mu = 1 statement" % oc[-1]
    assert 1.8 < oT[-1] < 2.2, \
        'the solved PERIOD converges at order %.2f, not the 2.00 that makes ' \
        'the absorption story work' % oT[-1]
    assert oc[-1] - oT[-1] > 0.7, \
        'the mechanism IS the gap: the state gains an order (%.2f) over the ' \
        'period (%.2f) because the h^2 error is a frequency error the ' \
        'autonomous solve takes up. No gap, no explanation.' % (oc[-1], oT[-1])


def test_the_frequency_aware_diffusion_reaches_c_on_a_non_uniform_grid():
    """`c(f) -> c` as `f -> 0` on ANY grid: the frequency-aware PPV goes to
    the DC one, so what is left between them is the quadrature alone.

    ⚠ Until 2026-09-27 `frequency_aware_diffusion` weighted its samples with
    `np.diff(times)` -- the LEFT RECTANGLE, which `_period_weights` names as
    first order on a non-uniform grid -- while `diffusion_constant` (and so
    `c(0)`, returned directly) uses `_period_weights`.  On a smoothly varying
    3:1 grid, measured `c(0+)/c - 1` = -2.0e-3 / -1.0e-3 / -5.1e-4 at N =
    200 / 400 / 800, and 1e-15 with the shared weights.  (An ALTERNATING 3:1
    grid hides it -- 8.8e-13 -- because the rectangle's error is `(1/2)
    integral h'(t) y dt`, which averages out when `h` alternates.)
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: 0.3 * (u - u ** 3 / 3.0) + 0.1 * u * u)
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    n = 200
    w = 1.0 + 0.5 * np.sin(2.0 * np.pi * np.arange(n) / n)
    pss = PSS(cir, method='radau', reltol=1e-11)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=2 * np.pi, timestep=2 * np.pi / n, grid=w / w.sum(),
                  x0=np.array([2.0, 0.0]), maxiterations=60)
        pac = PAC(cir)
        c = pac.diffusion_constant(pss)
        cfa = pac.frequency_aware_diffusion(pss, 1e-7 / float(pss.period))
    assert abs(cfa / c - 1.0) < 1e-10, cfa / c - 1.0


def test_the_frequency_aware_diffusion_reads_a_modulated_source_on_the_ppvs_orbit():
    """A source that follows the orbit is read at each PPV sample's OWN
    state -- the orbit the PPV was computed on, `_ppv_states`.

    ⚠ Until 2026-09-27 `frequency_aware_diffusion` read it at `pss.waveform`,
    on the belief (its comment) that the frequency-aware PPV is solved on the
    solve's own orbit.  It is not, for trap and euler: their
    `factored_period()` is the monodromy TWIN's, and so is every PPV built on
    it.  The two orbits differ by the discretisation; measured c(0+)/c - 1 =
    -1.3e-4 on trap (1e-13 on gear, which has no twin), 1e-15 after.
    """
    cir, pss, pac, ov = _orbit_modulated_vdp('white', method='trap')
    assert pss.monodromy_twin() is not pss
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        c = pac.diffusion_constant(pss)
        cfa = pac.frequency_aware_diffusion(pss, 1e-7 / float(pss.period))
    assert abs(cfa / c - 1.0) < 1e-10, cfa / c - 1.0


def test_the_frequency_aware_samples_follow_a_re_solve():
    """`PAC._fa_samples` caches the frequency-aware PPV per (factored period,
    offset).  ⚠ Until 2026-09-27 the key was `id(pss.factored_period())` and
    the entry did not hold the period: re-solving the same PSS on another
    grid freed it, a new period could be born at the same address, and the
    cache handed back the OLD grid's samples (measured: 198 rows for a
    299-point solve, 2 re-solves in 6).  One PAC across re-solves must give
    what a fresh PAC gives, every time."""
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: 0.1 * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi
    pss = PSS(cir, method='gear', reltol=1e-10)
    pac = PAC(cir)
    f = 0.05
    for n in (120, 180, 120, 180, 120, 180):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / n, x0=np.array([2.0, 0.0]))
        got = pac._fa_samples(pss, f)
        want = PAC(cir)._fa_samples(pss, f)
        assert got.shape == want.shape, (n, got.shape, want.shape)
        assert np.array_equal(got, want), n


def test_gears_ppv_samples_are_second_order_on_a_non_uniform_index2_grid():
    """Item 3 of gear-first-class (2026-09-20).  On an INDEX-2 oscillator
    gear's forward solve keeps order 2 in both subspaces on any grid, its
    period and adjoint modes are second order on uniform and smooth grids --
    and its PPV SAMPLES were FIRST order on the non-uniform grid alone::

        gear, smooth grid     N=100      N=200      N=400      N=800
        c rel to radau      -5.0e-03   -2.4e-03   -1.2e-03   -5.8e-04   (2.1, 2.05, 2.0)
        gear, uniform       +1.2e-03   +3.3e-04   +8.7e-05   +2.2e-05   (3.7, 3.8, 3.9)
        trap, smooth        -1.2e-03   -2.9e-04   -7.1e-05   -1.8e-05   (4.1, 4.0, 4.0)

    The cause is `_ppv_propagate`'s index-2 fallback: with `G[A,Z]` singular
    the pair-consistent correction keeps the differential block only
    (`G[D,NZ]`), the coupling through the algebraic variables being a
    DERIVATIVE term that does not exist as a Schur complement.  The dropped
    term cancels between steps on a uniform grid -- where the fallback was
    priced "second order" -- and does not on a non-uniform one.  Fix: on a
    non-uniform solved-history grid at index >= 2 the samples come from the
    continuous adjoint's phase mode, `v_j = C_j^T q_j` normalised by
    `q_0^T C_0 xdot_0 = 1` -- the same object `floquet_modes` uses there,
    whose invariant quarters on this very fixture and grid (4.3e-4 / 1.1e-4 /
    2.9e-5 / 7.5e-6).  Measured, c relative to radau at N = 3200:
    -1.9e-4 / -2.2e-5 / +7e-8 / +8e-7 -- at the reference's floor from
    N = 400, 130x closer at 800.  AND a second, independent half: every
    period integral in `PAC` weighted its samples with `diff(times)`, the
    LEFT RECTANGLE rule -- the trapezoid rule on a uniform periodic grid and
    FIRST order on a non-uniform one by itself (`PAC._period_weights`).
    With the fixed samples alone PAC's `c` still read -1.9e-3 / -9.3e-4 /
    -4.6e-4; with periodic trapezoid weights over the same samples +4.0e-5
    / +3.6e-5 / +1.3e-5 at N = 200 / 400 / 800, the reference's floor.

    ⚠ WHY TRAP AND RADAU NEVER PAID EITHER: their `info['times']` on this
    "smooth grid" run is UNIFORM (h = T/200 everywhere, 1.1 s off the
    waveform's nodes) -- `factored_period_full` / `_dirk` and the TR-BDF2
    twin replay the converged orbit on a `linspace` grid whatever grid the
    solve used.  Their PPV, modes and noise surfaces on a non-uniform grid
    are uniform-grid replays; accurate, and not on that grid.  Recorded
    here, not changed.

    ⚠ Three instruments were rejected on the way, and the record matters
    more than the fix: cross-method comparison of the PPV WAVEFORM on this
    fixture is not one (the two kinds' sample waveforms are not related by
    any shift or sign, while their maxima and c agree), a noise VOLTAGE
    inside the C-V loop gives c = 0 exactly for every method (a voltage
    perturbation of an index-2 constraint is a differentiated input the PPV
    projection cannot represent -- recorded, not chased), and the continuous
    q pairs with CHARGE perturbations: the state-space PPV is `C^T q`, and
    normalising `q` by `q . xdot` alone lands 16.7x off with the right shape.
    ⚠ The van der Pol fixture is index 2 through a DC source INSIDE a
    capacitor loop, which also exercises the constant-source period column
    fixed in 15579a9.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric

    def fix():
        cir = SubCircuit()
        cir.add_node('v')
        cir.add_node('a')
        cir['C'] = C('v', gnd, c=4.0)
        cir['L'] = L('v', gnd, L=0.25)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['vo'] = VS('v', 'a', v=0.5)
        cir['Ca'] = C('a', gnd, c=1.0)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        return cir
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    c0 = fix()
    names = [str(x) for x in c0.nodes]
    iref = c0.get_node_index(gnd)
    red = [i for i in range(c0.n) if i != iref]
    x0 = np.zeros(c0.n - 1)
    x0[red.index(names.index('v'))] = 2.0
    x0[red.index(names.index('a'))] = 1.5

    ## it IS index 2 (the test's own projector criterion), or this proves nothing
    xf = np.zeros(c0.n)
    xf[names.index('v')], xf[names.index('a')] = 2.0, 1.5
    Cm = np.asarray(c0.C(xf), float)[np.ix_(red, red)]
    Gm = np.asarray(c0.G(xf), float)[np.ix_(red, red)]
    _U, sv, Vt = np.linalg.svd(Cm)
    d = int(np.sum(sv > len(red) * sv[0] * np.finfo(float).eps))
    Nn = Vt[d:].T
    s2 = np.linalg.svd(Nn.T @ Gm @ Nn, compute_uv=False)
    assert s2[-1] / max(s2[0], 1e-300) < 1e-10, 'the fixture is not index 2'

    def smooth(n):
        f = 1.0 + 0.5 * np.sin(2 * np.pi * np.arange(n) / n)
        return f / f.sum()

    def c_of(method, n, grid):
        cir = fix()
        ## ⚠ the period column is pinned to 'proportional' here: this test's
        ## subject is the PPV samples at a FIXED discretisation.  Under
        ## 'auto' (closing, then a proportional polish -- the seed period is
        ## 8 % low here) the two paths converge to points 0.7 ppm apart in
        ## the period, both within tolerance, and c at N = 200 moves by
        ## 1.5e-4 between them (measured +2.9e-5 vs -1.2e-4) -- the
        ## fixture's own sensitivity, not the samples
        p = PSS(cir, method=method, reltol=1e-10, period_column='proportional')
        kw = dict(period=T0, timestep=T0 / n, maxiterations=60, x0=x0)
        if grid is not None:
            kw['grid'] = grid
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            p.solve(**kw)
            assert p.converged
            p.ppv()
            c = float(np.real(PAC(cir).diffusion_constant(p)))
        return c, any('algebraic block G[A,Z] is singular' in str(w.message) for w in rec)
    cref, _ = c_of('radau', 3200, None)
    e200, w200 = c_of('gear', 200, smooth(200))
    e400, w400 = c_of('gear', 400, smooth(400))
    assert w200 and w400, 'the index-2 fallback must have fired'
    r200, r400 = abs(e200 / cref - 1.0), abs(e400 / cref - 1.0)
    assert r200 < 1e-4 and r400 < 8e-5, (r200, r400)      # was 2.4e-3 / 1.2e-3
    ## the uniform grid keeps the exact-transpose samples (second order there)
    eu, wu = c_of('gear', 400, None)
    assert wu and abs(eu / cref - 1.0) < 2e-4, eu / cref - 1.0    # measured 8.7e-5


def test_radaus_oscillator_surfaces_on_a_genuinely_non_uniform_grid_are_order_five_except_where_the_period_quadrature_caps_them():
    """Item (c) after "honour the caller's grid" (2026-09-21): what the
    one-step kinds' PPV, modes and diffusion constant do on grids they now
    actually see.  Every earlier "radau exact on the 3:1 grid" number was a
    uniform replay; these are not.  Van der Pol, radau, `c` against radau
    N = 3200 uniform, phase multiplier and mode invariant::

        grid        N=100      N=200      N=400      N=800     order   |lam0|-1  inv
        uniform    -8.9e-10   -3.3e-11   +4.3e-12   +6.0e-13   ~5      2e-11     1e-15
        alt 2:1    -2.6e-09   -9.1e-11   -2.7e-12   -4.2e-13   ~5      2e-11     1e-15
        alt 3:1    -5.2e-09   -1.6e-10   -1.1e-11   +6.3e-12   ~5      2e-11     1e-15
        smooth     +4.9e-05   +1.3e-05   +3.3e-06   +8.4e-07   2.0     2e-11     1e-15

    ⚠ THE SMOOTH ROW IS THE PERIOD QUADRATURE, NOT RADAU.  The periodic
    trapezoid rule is spectrally accurate on a uniform grid and on an
    alternating one (two interleaved uniform sums) and genuinely O(h^2) on
    a smoothly varying grid, so there `c` is second order for EVERY method
    -- `_period_quadrature`'s own "second order caps radau", measured on a
    real replay (the index-2 twin read the same: +4.7e-5 / 1.3e-5 / 3.2e-6
    / 8.1e-7).  ⚠ LIFTED THE SAME DAY (2026-09-21, Andreas: "Do as you
    recommend"): on a non-uniform EVENT-FREE grid the period weights are a
    periodic cubic spline's (`periodic_spline_weights`), and the smooth
    row reads -5.5e-9 / -2.4e-10 / -7.3e-12 / +2.6e-12 -- the reference's
    floor from N = 400.  A landed edge keeps the trapezoid (a spline rings
    through a kink), and a uniform grid never uses it.  trbdf2 (and trap
    through its twin)
    are second order on every grid with no penalty for alternation, and
    their phase multiplier leaves 1 at O(h^2) on the smooth grid only
    (1.1e-4 -> 1.7e-6), which is warned.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric

    def vdp():
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=4.0)
        cir['L'] = L('v', gnd, L=0.25)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        return cir
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    C_REF = [None]                # radau, uniform, N = 1600: within 1e-13 of N = 3200

    def run(fr, n):
        cir = vdp()
        p = PSS(cir, method='radau', reltol=1e-10)
        kw = dict(period=T0, timestep=T0 / n, maxiterations=60,
                  x0=np.array([2.0, 0.0]), break_events=False)
        if fr is not None:
            kw['grid'] = fr
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(**kw)
            assert p.converged
            p.ppv()
            c = float(np.real(PAC(cir).diffusion_constant(p)))
            fp = p.factored_period()
            assert (p._period_quadrature(fp) is not None) == (fr is not None)
            modes = p.floquet_modes(p)
        lam0 = min(abs(abs(md['lam']) - 1.0) for md in modes)
        if fr is None:
            return c, lam0
        return c / C_REF[0] - 1.0, lam0

    def alt(n, r):
        f = np.tile([r, 1.0], n // 2)
        return f / f.sum()

    def smooth(n):
        f = 1.0 + 0.5 * np.sin(2 * np.pi * np.arange(n) / n)
        return f / f.sum()
    ## ⚠ the reference is computed, not quoted: a 7-digit printed value put
    ## a flat 5.3e-8 under every row (measured before this line existed)
    C_REF[0], _ = run(None, 1600)
    ## 3:1 alternating: order 5, at the reference's floor
    e200, l200 = run(alt(200, 3.0), 200)
    e400, l400 = run(alt(400, 3.0), 400)
    assert abs(e200) < 2e-9 and abs(e400) < 2e-10, (e200, e400)
    assert l200 < 1e-9 and l400 < 1e-9, (l200, l400)
    ## smooth: was the trapezoid's second order (1.3e-5 / 3.3e-6), now the
    ## spline's -- at the reference's floor
    s200, _ = run(smooth(200), 200)
    s400, _ = run(smooth(400), 400)
    assert abs(s200) < 2e-9 and abs(s400) < 2e-10, (s200, s400)

def test_a_noise_source_on_an_index2_constraint_is_named_not_silently_zero():
    """A voltage noise in series with the DC source INSIDE a capacitor loop
    perturbs an index-2 constraint: a differentiated input whose response is
    a charge jump, which the PPV projection cannot represent -- and
    `diffusion_constant` came back EXACTLY 0 for every method, with no
    message (found 2026-09-20 while measuring gear at index 2; the same
    circuit's current noise at the node gives 3.2e-9).  `_white_diffusion_at`
    now names it once: the circuit is index >= 2 at the orbit point
    (`G[A,Z]` singular -- the PPV fallback's own test, computed here so every
    kind is covered, not only solved-history) and `CY` has power on an
    algebraic row.  Silent for a source on a differential row.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric

    def fix(noisy_vs):
        cir = SubCircuit()
        cir.add_node('v')
        cir.add_node('a')
        cir['C'] = C('v', gnd, c=4.0)
        cir['L'] = L('v', gnd, L=0.25)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['vo'] = VS('v', 'a', v=0.5, noisePSD=(1e-6 if noisy_vs else 0.0))
        cir['Ca'] = C('a', gnd, c=1.0)
        if not noisy_vs:
            cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        return cir
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    for noisy_vs, want in ((True, True), (False, False)):
        cir = fix(noisy_vs)
        names = [str(x) for x in cir.nodes]
        iref = cir.get_node_index(gnd)
        red = [i for i in range(cir.n) if i != iref]
        x0 = np.zeros(cir.n - 1)
        x0[red.index(names.index('v'))] = 2.0
        x0[red.index(names.index('a'))] = 1.5
        p = PSS(cir, method='radau', reltol=1e-10)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T0, timestep=T0 / 200, maxiterations=60, x0=x0)
        assert p.converged
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            p.ppv()
            c = float(np.real(PAC(cir).diffusion_constant(p)))
        told = any('noise source sits on an algebraic row' in str(w_.message)
                   for w_ in rec)
        assert told is want, (noisy_vs, told, c)
        if noisy_vs:
            assert c == 0.0, c
        else:
            assert c > 1e-9, c


def test_a_native_glm_ppv_and_floquet_modes_are_on_the_state():
    """Under `monodromy='native'` a GLM oscillator's `ppv()` solved the
    bordered null vector of its NORDSIECK map and returned that `r*m`-wide
    object, normalised by a tangent read off the Nordsieck right null
    vector's first block: on this van der Pol the tangent came out
    [-0.0014, 0.0057] against [0.125, 8.07], `v(0)` 1.4e3 too large and `c`
    8e3 times too large (glm3, 120 points).  `floquet_modes` refused (its
    forward sampling went through the GLM's source-free forced replay).

    Now both read the map on the STATE, ``x_0 -> x_N``
    (`_GLMPeriod.state_map`): its null vector, and the samples along the
    orbit ``S_j^T lambda_j`` -- each node's Nordsieck costate taken back
    through a startup AT that node (`_glm_node_startups`), at node 0
    exactly the null vector.  Measured against radau at 800 points: `v(0)`
    2.5e-3 / 2.2e-5 (glm2 / glm3 at 60 points, second / third order), `c`
    2e-4 / 3e-5, the second multiplier 1.1e-3 / 7e-5 of radau's.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)

    def run(method, N, mono='native'):
        c = _vdp_asym()
        c['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        p = PSS(c, method=method, reltol=1e-12)
        p.monodromy = mono
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T0, timestep=T0 / N, x0=np.array([2.0, 0.0]),
                    maxiterations=60)
            v, info = p.ppv()
            cc = PAC(c, toolkit=circuit.numeric).diffusion_constant(p)
            fm = p.floquet_modes(nmodes=2)
        return p, np.asarray(v), info, float(np.real(np.asarray(cc).ravel()[0])), fm

    _pr, vr, _ir, cr, fmr = run('radau', 800)
    for method, tv, tc, tl in (('glm2', 5e-3, 5e-4, 2e-3),
                               ('glm3', 1e-4, 1e-4, 2e-4)):
        p, v, info, cc, fm = run(method, 60)
        m = p.cir.n - 1
        assert v.shape == (m,) and np.asarray(info['samples']).shape[1] == m
        assert np.max(np.abs(v - vr)) < tv * np.max(np.abs(vr)), (method, v, vr)
        assert abs(cc / cr - 1.0) < tc, (method, cc, cr)
        lam = sorted((abs(x['lam']) for x in fm), reverse=True)
        lamr = sorted((abs(x['lam']) for x in fmr), reverse=True)
        assert abs(lam[0] - 1.0) < 1e-6 and abs(lam[1] / lamr[1] - 1.0) < tl, (method, lam, lamr)
        ## the samples' first node is `M^T v`: the null vector, up to what
        ## the border absorbs of the discrete unit multiplier's offset
        ## (measured 4.9e-7 of `|v|` for glm2 at 60 points)
        assert np.max(np.abs(np.asarray(info['samples'])[0] - v)) < 1e-5 * np.max(np.abs(v))
        ## and the frequency-aware PPV runs on the same map, reaching `ppv`
        ## as the offset vanishes
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            vf, _inf = p.frequency_aware_ppv(1e-9 / float(p.period))
        assert np.max(np.abs(np.asarray(vf) - v)) < 1e-6 * np.max(np.abs(v))


def test_ppv_gives_each_of_its_warnings_once_a_call():
    """The review's W1 (2026-10-01): `ppv()` collects its samples'
    warnings and gives each once -- but the anchor's own equation-row fill
    ran outside that record, so on an index-2 oscillator (a DC source in a
    capacitor loop) "the algebraic block G[A, Z] is singular" came TWICE
    from one call.  Now once."""
    import collections
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir.add_node('a')
    cir['C'] = C('v', gnd, c=4.0)
    cir['L'] = L('v', gnd, L=0.25)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: (u - u ** 3 / 3.0) + 0.3 * u * u)
    cir['vo'] = VS('v', 'a', v=0.5)
    cir['Ca'] = C('a', gnd, c=1.0)
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    names = [str(x) for x in cir.nodes]
    red = [i for i in range(cir.n) if i != cir.get_node_index(gnd)]
    x0 = np.zeros(cir.n - 1)
    x0[red.index(names.index('v'))] = 2.0
    x0[red.index(names.index('a'))] = 1.5
    p = PSS(cir, method='gear', reltol=1e-10, period_column='proportional')
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=T0, timestep=T0 / 200, maxiterations=60, x0=x0)
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        p.ppv()
    count = collections.Counter(str(r.message) for r in rec)
    assert any('G[A, Z] is singular' in k for k in count), list(count)
    assert max(count.values()) == 1, count
