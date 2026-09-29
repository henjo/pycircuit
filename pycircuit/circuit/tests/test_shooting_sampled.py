"""Shooting tests: shooting sampled.  Split out of test_analysis_shooting.py on
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
from pycircuit.circuit.tests._shooting_elements import (_NuModNoise,
    _Flicker08,
    _Flicker20,
    _MixedSlopeLo,
    _NuMult,
    _SgnAmpFlicker,
    _SgnPsdFlicker,
    _SwitchFlickerHdl,
    _SwitchHdl)
from pycircuit.circuit.tests._shooting_fixtures import _noise_seam
from pycircuit.circuit.tests._shooting_fixtures import (_Flicker,
    _KB,
    _ModLorentzCtl,
    _ModLorentzSigned,
    _ModLorentzThermal,
    _TwoModLorentz,
    _TEMP,
    _a9_vdp,
    _jitter_sampler,
    _mixed_exponent_rc,
    _pulse_clocked_sampler,
    _sampler_fixture,
    _sampler_fixture_method,
    _sw)


def test_the_sample_series_evaluates_a_per_band_source_only_where_it_must():
    """`sampled_variance` / `sampled_noise` with a colour that is not a power
    law (2026-09-26; Andreas: "Speeding up sampled_variance for 2b-type
    sources").  Every band frequency and sideband evaluated EVERY element
    at every injection point, again for every instant: 5.5 M leaf
    evaluations for 38 band frequencies (75 of 79 s).  Now the element
    alone (`NoiseComponents.one_element_cy`, stamped straight into the reduced matrix),
    the root cached per frequency, and the component classified once:
    STATIONARY (one point), SEPARABLE (the root per point once, times
    ``sqrt(s(w))`` from one point), else every point.  Measured on a driven
    RC, two instants, 38 band frequencies: 63.8 -> 2.3 s (stationary,
    bit-identical), 65.6 -> 1.6 s (separable, 2.2e-16), 65.2 -> 6.6 s (the
    shape moving, bit-identical).  The separable source's element is called
    7203 times against 1084050 on the general path."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    T = 1e-6

    def build(shape):
        c = SubCircuit()
        for nd in ('lo', 'out'):
            c.add_node(nd)
        c['Vlo'] = VSin('lo', gnd, va=1.0, vo=1.5, freq=1.0 / T)
        c['Ro'] = R('out', gnd, r=1e3, noisy=False)
        c['Co'] = C('out', gnd, c=0.5e-9)
        c['n'] = _ModLorentzCtl('out', gnd, 'lo', gnd, noisePSD=1e-20,
                                tau=3e-7, k=1.0, shape=shape)
        pss = PSS(c, method='radau', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 100, maxiterations=40)
        return pss, PAC(c, toolkit=circuit.numeric)

    def sv(pss, pac):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return pac.sampled_variance(pss, 1, [0.3 * T, 0.7 * T], 1e-4 / T,
                                        0.5 / T, points_per_decade=10)

    calls = [0]
    orig_cy = _ModLorentzCtl.CY

    def counting(self, x, w, epar=None):
        calls[0] += 1
        return orig_cy(self, x, w, epar)
    _ModLorentzCtl.CY = counting
    try:
        ## separable: about one evaluation per frequency, not one per point
        pss, pac = build(0.0)
        K = len(pac._stage_states(pss, pss._state_map()))
        calls[0] = 0
        fast = sv(pss, pac)
        n_fast = calls[0]
        _noise_seam(pac, separable=staticmethod(
            lambda Cs, tol=1e-9: False))                   # the general path
        calls[0] = 0
        general = sv(pss, pac)
        n_general = calls[0]
    finally:
        _ModLorentzCtl.CY = orig_cy
    assert np.max(np.abs(fast / general - 1.0)) < 1e-12, fast / general - 1.0
    ## (7203 against 1084050: the fit and the classification, then one per
    ## band frequency and sideband -- against one per injection point each)
    assert n_general > 100 * n_fast, (n_fast, n_general, K)
    ## the shape moving: the one-element stamp equals the tree walk, and the
    ## per-frequency cache changes nothing
    pss, pac = build(0.5)
    st = pac._stage_states(pss, pss._state_map())[:7]
    a = pac._noise_components(pss, st).one_element_cy(('n',), 2 * np.pi * 3e5)
    ## ⚠ each equality holds VACUOUSLY if its patch is never looked up (a
    ## refactor that stops reaching the seam through the factory): count
    walked, uncached = [], []
    walk = _noise_seam(pac, leaf_access=staticmethod(
        lambda cir, key, irn: walked.append(1)))
    b = walk(pss, st).one_element_cy(('n',), 2 * np.pi * 3e5)
    del pac._noise_components
    assert walked, 'the patched leaf_access was never consulted'
    assert np.array_equal(a, b)
    cached = sv(pss, pac)
    ## (the moving shape's cache is `perband_root`'s since 2026-09-28, which
    ## carries an element's signed columns; this element states none)
    def perband_root(self, key, mode):
        uncached.append(1)
        return lambda w: self.psd_sqrt(self.one_element_cy(key, float(w)))
    _noise_seam(pac, perband_root=perband_root)
    assert np.array_equal(sv(pss, pac), cached)
    assert uncached, 'the patched perband_root was never consulted'


def test_the_sampled_variance_is_the_covariance_at_that_instant_for_white_sources():
    """`PAC.sampled_variance` -- one seeded adjoint per (instant, series
    frequency), every sideband from its phases -- against the validated
    per-step Lyapunov route at the SAME instants.  Measured at 400 points:
    held 0.998721 kT/C both ways (1e-6), tracking 0.846606 against 0.846745
    (1.6e-4, covariance's own O(h/tau) floor).  The white series PSD is flat,
    so the band [fmin, f0/2] holds `1 - fmin/(f0/2)` of the full variance;
    fmin = 1e-6 f0 makes that 2e-6.  ⚠ Checked to fail: sources sampled one
    step early read 1.30 kT/C held; sidebands cut to N/8 read 0.77 tracking.
    """
    cir, pss, io, pac, T = _sampler_fixture(lambda c: c.__setitem__('S0', _sw()))
    ktc = _KB * _TEMP / 100e-12
    _K0, Ks = pac.covariance(pss, samples=True)
    cov = np.array([np.asarray(k, dtype=float)[io, io] for k in Ks]) / ktc
    grid = np.asarray(pss.factored_period().times, dtype=float)
    kh, kt = 149, 40
    f0 = 1.0 / T
    fmin = 1e-6 * f0
    ## a hair off the grid points: the nearest ones are used and reported
    want = [grid[kh] + 1e-12 * T, grid[kt] - 1e-12 * T]
    var = pac.sampled_variance(pss, io, want, fmin, 0.5 * f0,
                               points_per_decade=10) / ktc
    var = var / (1.0 - fmin / (0.5 * f0))
    np.testing.assert_allclose(pac.sampled_instants, [grid[kh], grid[kt]],
                               rtol=0, atol=0)
    assert abs(var[0] / cov[kh] - 1.0) < 1e-5, (var[0], cov[kh])
    assert abs(var[1] / cov[kt] - 1.0) < 5e-4, (var[1], cov[kt])
    assert abs(cov[kh] - 1.0) < 3e-3 and cov[kt] < 0.9


def test_the_sampled_series_psd_is_the_pnoise_fold_on_an_lti_circuit_white_and_one_over_f():
    """On a time-invariant circuit the sample series folds the output PSD:
    `S(f; t0) = sum_n S_out(|f + n f0|)`, whatever `t0`.  The switch with
    `gon = goff` is a noisy resistor; a 1/f source sits at the output.
    Measured: 1e-10 at 0.137 f0 and 0.01 f0 over |n| <= 10 against
    `pnoise` at the same bands (which evaluates `CY` per band itself)."""
    import warnings

    def els(c):
        c['S0'] = _sw(gon=1e-3, goff=1e-3)
        c['F0'] = _Flicker('out', gnd, i=0.0, noisePSD=1e-22, fref=1.0)
    cir, pss, io, pac, T = _sampler_fixture(els)
    f0 = 1.0 / T
    fs = np.array([0.137, 0.01]) * f0
    S = pac.sampled_noise(pss, io, [0.37 * T], fs, maxsidebands=10)[0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ref = [sum(float(np.real(pac.pnoise(pss, abs(f + k * f0), io,
                                            maxsidebands=20)[0]))
                   for k in range(-10, 11)) for f in fs]
    np.testing.assert_allclose(S, ref, rtol=1e-8, atol=0)


def test_independent_noise_sources_add_in_the_coloured_folds():
    """⚠⚠ THE DEFECT THIS FIXES (2026-09-15).  `pnoise(cyclostationary=True)`
    took ONE square root of the summed `CY(x(t), w)` per band, so independent
    sources under different modulations did not add: switch white `4kTg(t)`
    plus a constant 1/f source at the same node read +7.3 % of the total at
    0.013 f0 (+2.2 % at 0.137 f0) over white-only + flicker-only, and the
    sampled variance in hold was +9.3 % with both sources inside ONE element.
    Now one root per element x {white, coloured}: pnoise additive to 1e-15;
    the one-element sampled variance equals the two-element one to 2.3e-11.
    """
    import warnings
    ktc = _KB * _TEMP / 100e-12
    flick = lambda c: c.__setitem__('F0', _Flicker('out', gnd, i=0.0,
                                                   noisePSD=1e-22, fref=1.0))
    fx = {
        'white': _sampler_fixture(lambda c: c.__setitem__('S0', _sw())),
        'flicker': _sampler_fixture(lambda c: (c.__setitem__('S0', _sw(kb=0.0)),
                                               flick(c))),
        'both': _sampler_fixture(lambda c: (c.__setitem__('S0', _sw()), flick(c))),
    }
    f0 = 1.0 / fx['white'][4]
    p = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        for k, (_c, pss, io, pac, _T) in fx.items():
            p[k] = float(np.real(pac.pnoise(pss, 0.013 * f0, io, maxsidebands=60,
                                            cyclostationary=True)[0]))
    assert not [w for w in caught if 'not the sum of its elements' in str(w.message)]
    assert abs((p['both'] - p['white'] - p['flicker']) / p['both']) < 1e-9, p
    assert p['flicker'] / p['both'] > 0.05, p      # the flicker is not negligible

    one = _sampler_fixture(lambda c: c.__setitem__(
        'S0', _SwitchFlickerHdl('in', 'out', 'ck', gnd, gon=1e-3, goff=1e-9,
                                vth=0.0, vs=50e-3, temp=_TEMP, kb=_KB, kf=1e-20)))
    two = _sampler_fixture(lambda c: (
        c.__setitem__('S0', _sw()),
        c.__setitem__('S1', _SwitchFlickerHdl('in', 'out', 'ck', gnd, gon=1e-15,
                                              goff=1e-15, vth=0.0, vs=50e-3,
                                              temp=_TEMP, kb=0.0, kf=1e-20))))
    T = one[4]
    t_hold = 149 * T / 400
    v1 = one[3].sampled_variance(one[1], one[2], [t_hold], 1e-2 * f0, 0.5 * f0,
                                 points_per_decade=20)[0] / ktc
    v2 = two[3].sampled_variance(two[1], two[2], [t_hold], 1e-2 * f0, 0.5 * f0,
                                 points_per_decade=20)[0] / ktc
    assert abs(v1 / v2 - 1.0) < 1e-8, (v1, v2)
    assert v1 > 1.1, v1                  # the flicker is a visible share


def test_the_sampled_variance_integrates_a_power_law_between_its_points_exactly():
    """`PAC._loglog_integral` (2026-09-25): `sampled_variance`'s density is
    a power law between neighbouring points (linear in log-log), each
    interval integrated exactly -- so white, 1/f, a MOS-like 1/f^0.9 and
    1/f^2 are EXACT on any grid, where the linear trapezoid it replaced
    overestimated a 1/f band by (r - 1)^3 / 6 per point (2.4e-3 at 20 per
    decade here) and the trapezoid in ln f a white one by (ln r)^2 / 12
    (1.1e-3).  An interval with a non-positive end takes the trapezoid."""
    from scipy.integrate import trapezoid
    f = np.geomspace(1.0, 1e3, 61)                     # 20 per decade
    for p in (0.0, -1.0, -0.9, -2.0, 0.5):
        S = 3.0 * f ** p
        exact = 3.0 * (np.log(1e3) if p == -1.0
                       else (1e3 ** (p + 1.0) - 1.0) / (p + 1.0))
        got = PAC._loglog_integral(np.vstack((S, 2.0 * S)), f)
        assert np.allclose(got / (np.array([1.0, 2.0]) * exact), 1.0,
                           rtol=1e-13, atol=0), (p, got / exact)
    S = 1.0 / f
    assert trapezoid(S, f) / np.log(1e3) - 1.0 > 2e-3
    assert trapezoid(np.ones_like(f) * f, np.log(f)) / (1e3 - 1.0) - 1.0 > 1e-3
    got = PAC._loglog_integral(np.array([[1.0, 0.0, 2.0, 2.0]]),
                               np.array([1.0, 2.0, 3.0, 4.0]))
    assert abs(got[0] - (0.5 + 1.0 + 2.0)) < 1e-14


def test_the_sampled_variance_grows_as_ln_fmin_with_flicker_and_refuses_what_it_cannot_do():
    """A 1/f source makes the band variance grow by a constant per decade
    of `fmin` (measured 0.0010 kT/C per decade, flat 1e-4 -> 1e-1 f0), so
    `fmin` is required.  Refused: a band outside (0, f0/2], a series
    frequency of 0, sidebands above the grid's Nyquist, an oscillator, and
    a stage-method period map."""
    import warnings
    ktc = _KB * _TEMP / 100e-12
    cir, pss, io, pac, T = _sampler_fixture(lambda c: (
        c.__setitem__('S0', _sw(kb=0.0)),
        c.__setitem__('F0', _Flicker('out', gnd, i=0.0, noisePSD=1e-22, fref=1.0))))
    f0 = 1.0 / T
    v = [pac.sampled_variance(pss, io, [149 * T / 400], r * f0, 0.5 * f0,
                              points_per_decade=10)[0] / ktc
         for r in (1e-4, 1e-3, 1e-2)]
    inc = -np.diff(v)
    assert inc[0] > 0 and abs(inc[0] / inc[1] - 1.0) < 0.05, (v, inc)
    with pytest.raises(ValueError, match='series_fmin'):
        pac.sampled_variance(pss, io, [0.0], 0.0, 0.5 * f0)
    with pytest.raises(ValueError, match='series_fmin'):
        pac.sampled_variance(pss, io, [0.0], 1e-3 * f0, 0.6 * f0)
    with pytest.raises(ValueError, match=r'\(0, f0/2\]'):
        pac.sampled_noise(pss, io, [0.0], [0.0])
    with pytest.raises(ValueError, match='Nyquist'):
        pac.sampled_noise(pss, io, [0.0], [0.1 * f0], maxsidebands=400)
    ocir, opss = _a9_vdp()
    with pytest.raises(ValueError, match='OSCILLATOR'):
        PAC(ocir, toolkit=circuit.numeric).sampled_noise(opss, 0, [0.0], [0.01])


def _a8_buffer_chain(psd=4e-21, noisy=True, cap=2e-10, res=1e3, va=1.0, f0=1e6):
    """A8's driven buffer chain: three tanh stages into RC loads.

    ⚠⚠ THE LOAD CAPACITANCE IS NOT A FREE CHOICE, and getting it wrong cost a
    published number.  The sample series folds every sideband up to the GRID's
    Nyquist.  With `cap = 2e-12` the RC pole sits at 1/(2 pi R C) = 80 f0, far
    above the fold, so the grid -- not the circuit -- decides how much of the
    noise tail is counted: at npts = 1600 the folded variance is 3.38e-7 V^2
    at `maxsidebands = 10` against 1.83e-6 at 400, i.e. a ten-sideband read
    captures 18 % of the answer, and sigma_t climbed +14.66 / +8.71 / +4.89 %
    per grid doubling over npts 200 -> 1600.  `cap = 2e-10` puts the pole at
    0.80 f0, where the circuit's own bandwidth limits the fold and the answer
    converges (+0.24 %, then +0.14 % per doubling).
    """
    cir = SubCircuit()
    cir.add_node('in')
    cir['vs'] = VSin('in', gnd, va=va, freq=f0)
    prev = 'in'
    for k in range(3):
        nd = 'o%d' % k
        cir.add_node(nd)
        cir['B%d' % k] = BSource(prev, gnd, nd, gnd,
                                 i_func=lambda u: 5e-4 * np.tanh(3.0 * u))
        cir['R%d' % k] = R(nd, gnd, r=res, noisy=noisy)
        cir['C%d' % k] = C(nd, gnd, c=cap)
        if psd > 0:
            cir['n%d' % k] = IS(nd, gnd, i=0.0, noisePSD=psd)
        prev = nd
    return cir


def _a8_edge_variance(cir, npts, f0=1e6, method='radau', node='o2'):
    """`(variance at the rising crossing of `node`, slew there, lte_warned)`."""
    import warnings as _warnings
    T = 1.0 / f0
    pss = PSS(cir, method=method, reltol=1e-11)
    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter('always')
        pss.solve(period=T, timestep=T / npts, x0=np.zeros(cir.n - 1),
                  maxiterations=60)
    warned = any('not resolved at this accuracy' in str(w.message)
                 for w in caught)
    assert pss.converged, 'A8 chain PSS did not converge at npts = %d' % npts
    Xw = np.asarray(pss.waveform[1], dtype=float)
    grid = np.asarray(pss.factored_period().times, dtype=float)[:Xw.shape[1]]
    red = [str(n) for n in cir.nodes if str(n) != 'gnd!'].index(node)
    v = Xw[cir.get_node_index(node)]
    mid = 0.5 * (v.max() + v.min())
    cross = [k for k in range(1, len(v)) if (v[k - 1] - mid) < 0 <= (v[k] - mid)]
    assert cross, 'no rising crossing at %s (pk-pk %.5g)' % (node, float(np.ptp(v)))
    j = cross[0]
    slew = float((v[j] - v[j - 1]) / (grid[j] - grid[j - 1]))
    pac = PAC(cir, toolkit=circuit.numeric)
    ## ⚠ the gate that caught a SILENT all-zeros.  `pss.waveform` rows are the
    ## FULL MNA vector, while `_output_waveform_row` takes a REDUCED index and
    ## lifts it past irefnode itself -- so a full index quietly addressed a
    ## source's branch current and returned 0.0 for sources ON *and* OFF, a
    ## mis-addressed query wearing a clean noise floor's clothes.  The row the
    ## noise API resolves must be the node meant, and that is asserted before
    ## any number leaves this helper.
    chk = np.asarray(pac._output_waveform_row(pss, red), dtype=float)
    assert abs(float(np.ptp(chk)) - float(np.ptp(v))) <= 1e-9 * float(np.ptp(v)), \
        'output row mismatch: the API resolved pk-pk %.6f, o2 is %.6f' % (
            float(np.ptp(chk)), float(np.ptp(v)))
    var = float(np.asarray(pac.sampled_variance(pss, red, [grid[j]], 1e3,
                                                0.5 * f0,
                                                points_per_decade=5)).ravel()[0])
    return var, slew, warned


def test_the_additive_edge_jitter_is_method_independent_and_nothing_manufactures_its_floor():
    """A8's first additive number, with the three checks that decide whether
    it is an answer or an artefact.  Given the oscillator's waveform the
    buffer chain is a DRIVEN LPTV system, where `sampled_noise` applies (it
    refuses autonomous PSS by design, tested above); the oscillator's own
    share is `c`, which already ships.

    Measured 2026-09-16 on three tanh buffers, `gm R = 0.5`, 4e-21 A^2/Hz
    injected per stage, f0 = 1 MHz: sigma_t = 8.27e-11 s at o2's rising
    crossing, band [1 kHz, 500 kHz].

    ⚠ TWO GRIDS AGREEING IS NOT CONVERGENCE, so the gate is across METHOD
    FAMILIES.  gear-2 declares this fixture unresolved (local truncation error
    5.11e6 x tolerance at 400 points, 1.29e6 at 800) -- and radau, which
    reports nothing at any grid tried, lands on the same number anyway: radau
    100/200/400 give 8.249493 / 8.266431 / 8.274739e-11 s against gear 400/800
    at 8.263152 / 8.275126e-11 and trap 400 at 8.268268e-11.  That is 0.17 %
    over three families and a 4x grid range, with radau 400 and gear 800
    agreeing to 0.005 %.  The gear warning is in any case mostly the reltol
    chosen here: LTE/tol runs 5.11e6 -> 5.98e4 -> 600 as reltol goes 1e-11 ->
    1e-9 -> 1e-7 while sigma_t is unchanged to seven digits.

    ⚠⚠ WHAT THIS TEST EXISTS TO PREVENT -- three claims of mine that a
    plausible-looking harness produced and that had to be withdrawn:

      * A GRID-BOUND HEADLINE.  The first fixture's pole sat at 80 f0, so the
        grid's Nyquist set the noise bandwidth and sigma_t rose ~5-15 % per
        doubling, leaving the quoted 5.712532e-11 s about 18 % under its own
        extrapolated limit.  Hence the convergence assertion here, which that
        fixture fails outright (it moves 8.7 % over the same pair).
      * A FLOOR CONTROL THAT CONTROLLED NOTHING.  "Sources off, 242x down"
        left the RESISTORS noisy (`R.CY` is `4kT/r`), so it compared injected
        noise against thermal noise and called the remainder an instrument
        floor.  The real control needs `noisy=False` too -- and it is much
        better than the claim it replaces: the variance is then EXACTLY zero.
      * A VACUOUS FALSIFIER.  "sigma_t * slew constant" cannot fail: sigma_t
        is DEFINED as sqrt(var)/|slew|, so the product is sqrt(var)
        identically (residual measured at 1.4e-20).  It tested only whether
        the variance moved under the knob.  Likewise `PSD x 4 -> 1.99690` is
        arithmetic given exact linearity in the source PSDs, not evidence.

    ⚠ WHAT IS STILL NOT VALIDATED, stated so this is not read as more than it
    is: converting a variance to a time by the local slew is a DEFINITION
    here, not a measured jitter.  An independent time-domain check (Monte
    Carlo scatter of the actual crossing) is not done and needs machinery
    that does not exist -- see A8 in `doc/pss_roadmap_260902.md`.
    """
    psd, res, f0 = 4e-21, 1e3, 1e6

    ## the headline, on the accurate family: radau reports no resolution
    ## trouble on this fixture at any grid tried
    v_r, slew_r, warned_r = _a8_edge_variance(_a8_buffer_chain(psd, res=res), 200)
    st_radau = np.sqrt(v_r) / abs(slew_r)
    assert not warned_r, 'radau now reports this orbit unresolved: the ' \
        'fixture has moved and the number below is no longer its own'
    assert abs(st_radau / 8.266431e-11 - 1.0) < 2e-3, st_radau

    ## 1. THE ANSWER, NOT THE DISCRETISATION.  A different family at a
    ##    different grid must land on the same number.  gear-2 is the least
    ##    accurate method here and it does declare itself unresolved -- which
    ##    is exactly why agreeing with it is worth asserting, and why two gear
    ##    grids agreeing with each other would not have been.
    v_g, slew_g, warned_g = _a8_edge_variance(_a8_buffer_chain(psd, res=res),
                                              400, method='gear')
    st_gear = np.sqrt(v_g) / abs(slew_g)
    assert warned_g, 'gear no longer reports this fixture unresolved; the ' \
        'cross-family check has quietly become a comparison of two ' \
        'well-resolved solves, and proves less than it claims'
    assert abs(st_gear / st_radau - 1.0) < 3e-3, \
        'radau 200 gives %.6e s, gear 400 gives %.6e s (%.2f %% apart): the ' \
        'number depends on the integrator, so it is not the circuit\'s' % (
            st_radau, st_gear, 100.0 * (st_gear / st_radau - 1.0))

    ## the structural checks need one grid only; radau 100 is the cheapest
    ## place to run four solves against each other
    v_all, _, _ = _a8_edge_variance(_a8_buffer_chain(psd, res=res), 100)

    ## 2. nothing is manufactured: silence every source, resistors included,
    ##    and the variance must be exactly zero -- the control that the
    ##    earlier "sources off" run only appeared to be
    v_quiet, _, _ = _a8_edge_variance(
        _a8_buffer_chain(0.0, noisy=False, res=res), 100)
    assert v_quiet == 0.0, \
        'with every source silent the variance is %.6e, not 0' % v_quiet

    ## 3. the sources-off reading is RESISTOR THERMAL NOISE, and its size
    ##    follows from 4kT/R with nothing fitted
    v_res, _, _ = _a8_edge_variance(
        _a8_buffer_chain(0.0, noisy=True, res=res), 100)
    therm = 4.0 * circuit.numeric.kboltzmann * float(defaultepar.T) / res
    assert abs((v_all / v_res) / (1.0 + psd / therm) - 1.0) < 1e-4, \
        'injected/thermal variance ratio %.6f against 1 + PSD/(4kT/R) = %.6f' % (
            v_all / v_res, 1.0 + psd / therm)

    ## 4. independent sources add -- the property a coloured-fold defect once
    ##    broke (test_independent_noise_sources_add_in_the_coloured_folds)
    v_inj, _, _ = _a8_edge_variance(
        _a8_buffer_chain(psd, noisy=False, res=res), 100)
    assert abs((v_inj + v_res) / v_all - 1.0) < 1e-9, (v_inj, v_res, v_all)
    assert v_inj / v_all > 0.99, \
        'the injected sources must dominate; they are %.4f' % (v_inj / v_all)


def _a8_one_stage(psd=4e-21, noisy=True, res=1e3, cap=2e-10, gm=5e-4, va=1.0,
                  f0=1e6):
    """ONE LINEAR stage -- A8's anchor fixture.  `i = gm v_in` into `R||C`,
    no tanh anywhere, so every quantity has a closed form."""
    cir = SubCircuit()
    cir.add_node('in')
    cir.add_node('o0')
    cir['vs'] = VSin('in', gnd, va=va, freq=f0)
    cir['B0'] = BSource('in', gnd, 'o0', gnd, i_func=lambda u: gm * u)
    cir['R0'] = R('o0', gnd, r=res, noisy=noisy)
    cir['C0'] = C('o0', gnd, c=cap)
    if psd > 0:
        cir['n0'] = IS('o0', gnd, i=0.0, noisePSD=psd)
    return cir


def test_the_edge_jitter_of_a_linear_stage_matches_its_closed_form_and_the_grid_truncation_it_names():
    """A8's ANCHOR: an absolute value to hit, not just agreement between
    methods.  The cross-family test above establishes that the number belongs
    to the circuit rather than to an integrator; it cannot establish that it
    is RIGHT, because the three-tanh chain has no analytic value at any grid.
    A single LINEAR stage does.

    ⚠ THE GAP THIS CLOSES WAS POINTED OUT BY THE PEER SESSION
    (`test-pycircuit-spectre-e2`, 2026-09-16), and the point generalises: their
    sampled `kT/C` sits against an exact 1.0, so it cannot drift by tens of
    percent without announcing itself as disagreement with the PHYSICS rather
    than with another tool.  My grid-bound headline had no such anchor, which
    is precisely why only grid refinement could catch it.

    White current sources `S = S_inj + 4kT/R` into `R||C` give
    `var = S R/(4C)` over all frequencies (the resistor's own share being
    exactly `kT/C`), the flat series PSD puts `1 - fmin/(f0/2)` of that inside
    the band, and the fold reaches only the GRID's Nyquist
    `f_N = (npts/2) f0`, so a known fraction is missing:
    `captured = (2/pi) arctan(f_N/f_c)` with `f_c = 1/(2 pi R C)`.

    ⚠ THE DEFICIT WAS NAMED BEFORE IT WAS MEASURED, per the campaign's
    pre-commitment rule: `(2/pi) f_c/f_N` predicts 2.533e-03 at npts 400 and
    1.266e-03 at 800, i.e. HALVING per doubling.  Measured 2.178e-03 and
    8.976e-04 (ratio 0.412), and with the truncation term included the
    variance matches the closed form to 3.6e-04 at both grids.  The slew
    matches `2 pi f0 A_out` to 1.1e-05 and sigma_t lands 3.614654e-11 s
    against an analytic 3.618558e-11 s.

    ⚠ The truncation term is asserted to be DOING WORK rather than
    decorating: dropping it puts the coarse grid out by more than the
    tolerance the corrected form is held to.

    ⚠ Still a DEFINITION, as above: this anchors the VARIANCE and the slew,
    not the claim that crossings scatter by sigma_t.  That needs the Monte
    Carlo named in A8.
    """
    psd, res, cap, gm, va, f0 = 4e-21, 1e3, 2e-10, 5e-4, 1.0, 1e6
    fmin = 1e3
    kT = circuit.numeric.kboltzmann * float(defaultepar.T)
    fc = 1.0 / (2.0 * np.pi * res * cap)
    var_full = (psd + 4.0 * kT / res) * res / (4.0 * cap)
    band = 1.0 - fmin / (0.5 * f0)
    a_out = gm * va * res / np.sqrt(1.0 + (2.0 * np.pi * f0 * res * cap) ** 2)

    got = {}
    for npts in (400, 800):
        var, slew, warned = _a8_edge_variance(
            _a8_one_stage(psd, res=res, cap=cap, gm=gm, va=va, f0=f0), npts,
            f0=f0, node='o0')
        assert not warned, 'radau reports the linear stage unresolved at %d' % npts
        captured = (2.0 / np.pi) * np.arctan(0.5 * npts * f0 / fc)
        got[npts] = (var, slew, var_full * band * captured)

    for npts in (400, 800):
        var, slew, pred = got[npts]
        assert abs(var / pred - 1.0) < 1e-3, \
            'npts %d: sampled variance %.9e against the closed form %.9e ' \
            '(%.4f)' % (npts, var, pred, var / pred)

    ## the truncation term earns its place: without it the coarse grid misses
    ## by more than the tolerance the corrected form just passed
    d400 = 1.0 - got[400][0] / (var_full * band)
    d800 = 1.0 - got[800][0] / (var_full * band)
    assert d400 > 1.5e-3, \
        'the untruncated closed form is already within %.2e at npts 400, so ' \
        'this fixture no longer demonstrates the grid truncation' % d400
    ## and it halves with the grid, as (2/pi) f_c/f_N says it must
    assert 0.30 < d800 / d400 < 0.60, (d400, d800)

    ## the conversion factor is analytic here too
    var, slew, _ = got[400]
    assert abs(abs(slew) / (2.0 * np.pi * f0 * a_out) - 1.0) < 1e-4, slew
    st = np.sqrt(var) / abs(slew)
    assert abs(st / (np.sqrt(var_full * band) / (2.0 * np.pi * f0 * a_out))
               - 1.0) < 2e-3, st


def _a8_mc_edge_scatter(npts, nper, seed, s_tot, scale=1.0, res=1e3, cap=2e-10,
                        gm=5e-4, va=1.0, f0=1e6, burn=8):
    """Run a NOISY transient and return `(scatter of the rising crossings,
    count, mean spacing)`.

    ⚠ THE INJECTION CONVENTION IS THE WHOLE RISK, and it is not chosen here on
    reasoning.  A Monte Carlo that injects noise with the same one-sided /
    two-sided convention the analysis assumes cannot validate that analysis --
    this library has been bitten by exactly that, a MC and `diffusion_constant`
    agreeing to 0.9965 while both were 2x wrong.  So `sigma_i^2 = S/(2h)` was
    pinned against `kT/C`, which neither route can influence: backward Euler on
    R||C predicts `Var = (kT/C)/(1 + h/2tau)` = 0.993789 and measured
    1.005742 +/- 0.014891, with the rival `S/h` reading (1.987578) rejected at
    ~66 sigma.
    """
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import EulerIntegrator
    T = 1.0 / f0
    h = T / npts
    rng = np.random.default_rng(seed)
    draws = rng.normal(0.0, scale * np.sqrt(s_tot / (2.0 * h)),
                       size=int(nper * npts) + 8)
    cir = _a8_one_stage(psd=0.0, noisy=False, res=res, cap=cap, gm=gm, va=va,
                        f0=f0)
    cir['In'] = IS(gnd, 'o0', i=0.0)
    cir['In'].function.f = lambda t: float(draws[int(t / h + 0.5)])
    tran = Transient(cir, integrator=EulerIntegrator())
    sol = tran.solve(tend=nper * T, timestep=h, x0=np.zeros(cir.n),
                     fixed_timestep=True)
    v = np.asarray(sol.v('o0', gnd), dtype=float)
    t = np.arange(len(v)) * h
    ## ⚠ the threshold is ANALYTIC (0.0 -- a filtered sine with no DC path),
    ## not the midpoint of the record.  Taking it from the record let startup
    ## bias it, which moved the crossing off the steepest point and cost about
    ## half of an apparent deficit before it was caught.
    cross = []
    for k in range(1, len(v)):
        if t[k - 1] >= burn * T and v[k - 1] < 0 <= v[k]:
            frac = (0.0 - v[k - 1]) / (v[k] - v[k - 1])
            cross.append(t[k - 1] + frac * (t[k] - t[k - 1]))
    cross = np.asarray(cross)
    assert len(cross) > 2, 'only %d crossings' % len(cross)
    idx = np.round((cross - cross[0]) / T)
    resid = (cross - cross[0]) - idx * T
    slope = float(np.polyfit(idx, resid, 1)[0])
    return (float(np.std(resid - slope * idx, ddof=1)), len(cross),
            float(np.mean(np.diff(cross))))


def test_the_edge_jitter_is_the_scatter_of_a_noisy_transients_crossings_not_just_a_variance():
    """A8's LAST GAP CLOSED: sigma_t stops being a definition.

    Everything else in A8 computes `sigma_t := sqrt(var)/|slew|` and the
    roadmap said, in as many words, that converting a variance to a time this
    way is a DEFINITION and not a measured jitter -- nothing showed that a
    noisy transient's crossings actually scatter by it.  This runs the
    transient and measures the scatter.

    Measured 2026-09-16, pooled over 9 seeds x 3 grids (npts 200/400/800) and
    3528 crossings: sigma_t(MC) / closed form = 1.01099 +/- 0.01156,
    consistent with 1.000 at 0.95 sigma.  A factor of 2 in the variance is
    excluded at 26-35 sigma, a factor of 4 at 44-86 sigma.  The seed spread
    (0.03467) matches 1/sqrt(2N) for N = 392 (0.03571) to 3 %, so the scatter
    is sampling noise and the error model checks itself.

    ⚠ WHAT THIS TEST CAN AFFORD is one seed and ~192 crossings, i.e. a 5.1 %
    sampling error -- so the claim it GATES is not the 1 % agreement but that
    the conversion is not out by a FACTOR, which it places at 5.3 sigma.  The
    precise statement lives in `doc/pss_log_260902.md`; do not tighten the
    band here without buying the crossings to pay for it.  The seed is fixed,
    so this is deterministic rather than flaky (this seed sits ~1 sigma low,
    at 0.9699).

    ⚠⚠ TWO THINGS THAT WENT WRONG BUILDING THIS, both now gated below:
      * A HARNESS WHOSE CONTROLS FAILED THEIR OWN CRITERIA and one of them fed
        the measurement.  Run from `x0 = 0` with no burn-in, "crossing spacing
        = T" read 473 ppm off (a drift of 13x sigma_t per period, which would
        have swamped everything), the "floor" read 41x sigma_t, and the
        record's midpoint -- used as the crossing threshold -- was biased
        enough to move the measured slew to 0.9948 of analytic.
      * A REFINEMENT TEST TOO COARSE TO SEE WHAT IT LOOKED FOR.  Chasing a
        1.2 % deficit across npts 200/400/800 was hopeless when each grid
        carries a 2.1 % sampling error; the "deficit" was a 0.95 sigma
        fluctuation under an error bar I had underestimated by taking the sem
        from three seeds that happened to agree.  It is withdrawn.
    """
    psd, res, cap, gm, va, f0 = 4e-21, 1e3, 2e-10, 5e-4, 1.0, 1e6
    kT = circuit.numeric.kboltzmann * float(defaultepar.T)
    s_tot = psd + 4.0 * kT / res
    T = 1.0 / f0

    ## the shipped frequency-domain route, on the same fixture
    var_f, slew_f, _ = _a8_edge_variance(
        _a8_one_stage(psd, res=res, cap=cap, gm=gm, va=va, f0=f0), 200,
        f0=f0, node='o0')
    sig_freq = np.sqrt(var_f) / abs(slew_f)

    ## and the time domain, which knows nothing about it
    sig_mc, ncross, spacing = _a8_mc_edge_scatter(200, 200, 21, s_tot, res=res,
                                                  cap=cap, gm=gm, va=va, f0=f0)
    assert ncross > 150, 'only %d crossings: the run was short' % ncross
    ratio = sig_mc / sig_freq
    assert abs(ratio - 1.0) < 0.20, \
        'the crossings scatter by %.6e s, the frequency route says %.6e s ' \
        '(ratio %.4f); the sampling error here is 5.1 %%, so this is a ' \
        'factor, not noise' % (sig_mc, sig_freq, ratio)

    ## ⚠ and the alternative the whole exercise exists to exclude
    assert abs(ratio - 1.0 / np.sqrt(2.0)) > 0.15 and abs(ratio - np.sqrt(2.0)) > 0.15, \
        'ratio %.4f is consistent with a factor of 2 in the variance' % ratio

    ## the crossing ladder must be uniform, or the residuals are measured
    ## against the wrong thing
    assert abs(spacing / T - 1.0) < 1e-6, spacing

    ## POISON ALIVE: silence the source and the crossings must stop moving.
    ## Without this a dead source and a correct answer look identical.
    sig_q, nq, spacing_q = _a8_mc_edge_scatter(200, 40, 99, s_tot, scale=0.0,
                                               res=res, cap=cap, gm=gm, va=va,
                                               f0=f0)
    assert sig_q < 1e-3 * sig_mc, \
        'the noiseless harness floor is %.3e s against a signal of %.3e s' % (
            sig_q, sig_mc)
    assert abs(spacing_q / T - 1.0) < 1e-9, spacing_q


def test_the_across_period_correlation_is_the_cosine_transform_of_the_sample_series():
    """A11 step 1: `rho_k` needs NO new machinery, and in particular not
    Demir 1996.

    ⚠⚠ THIS CORRECTS A11's OWN ENTRY, which filed `rho_k` under "weeks of
    work, Demir 1996 the starting point".  That conflated two things:
    `sampled_variance` is one number at one instant and genuinely supplies
    none of it -- but `sampled_noise` returns the PSD of the SAMPLE SERIES
    `y(t0 + kT)`, and the autocovariance of a series IS the cosine transform
    of its spectrum.  So for a DRIVEN chain `rho_k`, and with it the k-cycle
    and cycle-to-cycle metrics, follows from what already ships.  Only the
    LDO's non-stationary delay modulation is still Demir's problem.

    The fixture is built so the answer is known in closed form rather than
    only cross-checked: one LINEAR stage, `tau = RC = T` exactly, and an LTI
    noise path (gm, R, C constant -- the drive only sets the waveform), so the
    sampled series is AR(1) with `rho_k = exp(-k)` exactly.

    MEASURED 2026-09-16: the transform gives 0.368067 / 0.135404 / 0.049813 /
    0.018325 against 0.367879 / 0.135335 / 0.049787 / 0.018316 -- ratio
    1.0005 at EVERY lag, a constant offset rather than one growing with k.
    Since 2026-09-25 (R_0's low end on a log grid, for 1/f sources:
    `test_jitter_metrics_integrate_a_flicker_low_end_on_a_log_grid`)
    0.368072 / 0.135410 / 0.049819 / 0.018332, ratio 1.00052 .. 1.00088:
    R_0 moved +7e-6 on this white spectrum, which the small rho_4 carries
    53x.  Without the rectangle rho_4 reads 0.9893 (was 0.9889).
    ⚠ Independently, a Monte Carlo over 1176 noisy crossings (no PSS, no
    adjoint, no spectrum) agreed within 1 sigma at every lag it can resolve:
    0.346404 and 0.131894 at k = 1, 2, i.e. 0.74 and 0.12 sigma, with
    `sigma_t` matching at 0.9901.  ⚠ A8's own fixture could NOT have tested
    this -- at `tau = 0.2 T` the first correlation is `e^-5` = 0.0067, far
    under the ~0.03 floor of that many crossings, and a quantity below the
    floor of the instrument that would check it is not a test.

    ⚠ The `sqrt(2)` / `sqrt(6)` collapse below is the check on the ALGEBRA:
    an uncorrelated series must give exactly those, and the second-difference
    coefficients `[1, -2, 1]` are easy to get wrong in a way no fixture with
    correlation would reveal.
    """
    f0 = 1e6
    T = 1.0 / f0
    fmin, fmax = f0 / 2e4, 0.5 * f0
    nfreq = 151

    def solved(cap, psd=4e-21):
        cir = _a8_one_stage(psd=psd, cap=cap, f0=f0)
        pss = PSS(cir, method='radau', reltol=1e-11)
        pss.solve(period=T, timestep=T / 400, x0=np.zeros(cir.n - 1),
                  maxiterations=60)
        assert pss.converged, 'A11 fixture PSS did not converge (cap %g)' % cap
        Xw = np.asarray(pss.waveform[1], dtype=float)
        grid = np.asarray(pss.factored_period().times, dtype=float)[:Xw.shape[1]]
        v = Xw[cir.get_node_index('o0')]
        red = [str(n) for n in cir.nodes if str(n) != 'gnd!'].index('o0')
        j = [i for i in range(1, len(v)) if v[i - 1] < 0 <= v[i]][0]
        return cir, pss, red, grid, j, v

    ## tau = RC = T  ->  rho_k = exp(-k), exactly
    cir, pss, red, grid, j, v = solved(1e-9)
    pac = PAC(cir, toolkit=circuit.numeric)
    m = pac.jitter_metrics(pss, red, grid[j], fmin, fmax, kmax=4, nfreq=nfreq,
                           dc_rectangle=True)
    for k in range(1, 5):
        assert abs(m['rho'][k - 1] / np.exp(-k) - 1.0) < 3e-3, \
            'rho_%d = %.6f against exp(-%d) = %.6f' % (
                k, m['rho'][k - 1], k, np.exp(-k))
    assert m['instant'] == grid[j], (m['instant'], grid[j])
    ## the slope is the jitter family's one estimator, a local quadratic fit
    ## at the instant (`edge_slope`; a central difference until 2026-09-29),
    ## and `kmax` defaults to 8, as `oscillator_edge_jitter`'s does
    from pycircuit.circuit.shooting._numerics import edge_slope
    assert m['slew'] == edge_slope(grid, v, m['instant'])[0]
    import inspect
    assert inspect.signature(pac.jitter_metrics).parameters['kmax'].default == 8

    ## the metrics are functions of rho and sigma_t, by a route that does not
    ## go through R -- an algebra slip in the method would not survive this
    rho = np.asarray(m['rho'], dtype=float)
    np.testing.assert_allclose(m['k_cycle'],
                               np.sqrt(2.0 * (1.0 - rho)) * m['sigma_t'],
                               rtol=1e-9)
    assert abs(m['cycle_to_cycle'] /
               (np.sqrt(6.0 - 8.0 * rho[0] + 2.0 * rho[1]) * m['sigma_t'])
               - 1.0) < 1e-9

    ## ⚠ and the flag must DO something: without the [0, fmin) rectangle the
    ## error grows with k, exactly as truncation should (measured 0.9889 of
    ## analytic at k = 4, against 1.0005 with it)
    m_raw = pac.jitter_metrics(pss, red, grid[j], fmin, fmax, kmax=4,
                               nfreq=nfreq, dc_rectangle=False)
    assert abs(m_raw['rho'][3] / np.exp(-4) - 1.0) > 5e-3, \
        'the DC rectangle no longer changes anything (%.6f); either the band ' \
        'moved or the flag stopped working' % m_raw['rho'][3]

    ## UNCORRELATED LIMIT: tau = 0.02 T, so rho -> 0 and the two metrics must
    ## collapse to sqrt(2) and sqrt(6) times sigma_t
    cir2, pss2, red2, grid2, j2, _v2 = solved(2e-11)
    pac2 = PAC(cir2, toolkit=circuit.numeric)
    m2 = pac2.jitter_metrics(pss2, red2, grid2[j2], fmin, fmax, kmax=4,
                             nfreq=nfreq, dc_rectangle=True)
    assert np.max(np.abs(m2['rho'])) < 1e-3, m2['rho']
    np.testing.assert_allclose(np.asarray(m2['k_cycle']) / m2['sigma_t'],
                               np.sqrt(2.0), rtol=1e-5)
    assert abs(m2['cycle_to_cycle'] / m2['sigma_t'] - np.sqrt(6.0)) < 1e-5

    ## refusals.  The band and kmax ones raise before any spectrum is built,
    ## so they cost nothing; the flat-instant one has to compute S first, so
    ## it is asked for the cheapest possible one.
    with pytest.raises(ValueError, match='kmax'):
        pac.jitter_metrics(pss, red, grid[j], fmin, fmax, kmax=1)
    with pytest.raises(ValueError, match='series_fmin'):
        pac.jitter_metrics(pss, red, grid[j], 0.0, fmax)
    with pytest.raises(ValueError, match='series_fmin'):
        pac.jitter_metrics(pss, red, grid[j], fmin, 0.9 * f0)
    ## ⚠⚠ THE LAST GUARD IS ON THE LINEARISATION, NOT ON "FLATNESS", and the
    ## reason is a number I misread on the way here.  At a peak the central
    ## DIFFERENCE is -8.86e-07, but the SLOPE is that over 2h = -1.767e+02 --
    ## 3.58e-04 of the steepest 4.94e+05, not 1.8e-12.  On a 400-point grid
    ## 3.58e-04 is the smallest ratio ANY instant can have, since no grid
    ## point lands exactly on the peak, so a relative-slope threshold would
    ## only ever describe the grid and would move with N.  What does not move
    ## is whether the displacement stays small against the period: at this
    ## fixture's peak sigma_t is 0.179 T, so 100x the source PSD puts it past
    ## T/2 and the first-order picture is gone.
    cir3, pss3, red3, grid3, _j3, v3 = solved(1e-9, psd=4e-19)
    pac3 = PAC(cir3, toolkit=circuit.numeric)
    with pytest.raises(ValueError, match='first-order'):
        pac3.jitter_metrics(pss3, red3, grid3[int(np.argmax(v3))], fmin, fmax,
                            kmax=2, nfreq=21)


def test_jitter_metrics_integrate_a_flicker_low_end_on_a_log_grid():
    """`jitter_metrics` with a 1/f source (2026-09-25).  Its `R_k` were a
    trapezoid on a LINEAR grid, whose first interval spans the whole 1/f
    low end: at `nfreq = 301` on this sampler (a constant 1/f current on
    the switch, noiseless otherwise) R_0..4 changed by 2.6 .. 5.8 % against
    finer grids (0.6 .. 0.8 % at the default 601).  Now ``R_k = int S df +
    int S (cos - 1) df``: the first a power law between the points of a log
    grid AND the linear one below their crossover, the plain trapezoid on
    the linear grid above it; the second on the linear grid (it vanishes as
    f^2 at low f).  R_0..4 move by <= 8.8e-5 against 2x the linear and 2x
    the log points, and R_0 meets `sampled_variance` to -2.2e-5 (its own
    40-per-decade error).  The k-cycle and cycle-to-cycle metrics depend
    on ``R_0 - R_k`` alone and did not move."""
    import warnings

    def el(cir):
        cir['S'] = _SwitchFlickerHdl('in', 'out', 'ck', gnd, vth=0.0,
                                     vs=50e-3, temp=_TEMP, kb=0.0, gon=1e-3,
                                     goff=1e-9, kf=1e-22)
    cir, pss, io, pac, T = _sampler_fixture(el, npts=100)
    f0 = 1.0 / T
    fmin, fmax = 1e-3 * f0, 0.5 * f0
    N = len(pss.factored_period().steps)
    th = float(np.asarray(pss.factored_period().times)[int(0.6 * N)])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        m = pac.jitter_metrics(pss, io, th, fmin, fmax, kmax=4, nfreq=301)
        fine = pac.jitter_metrics(pss, io, th, fmin, fmax, kmax=4, nfreq=601,
                                  points_per_decade=80)
        sv = pac.sampled_variance(pss, io, [th], fmin, fmax)[0]
    assert abs(m['R'][0] / sv - 1.0) < 1e-4, (m['R'][0], sv)
    moved = np.abs(np.asarray(m['R']) / np.asarray(fine['R']) - 1.0)
    assert np.max(moved) < 3e-4, moved


def test_the_sampled_noise_runs_natively_under_the_stage_methods():
    """radau and trbdf2 inject the source at every STAGE, so the sampled
    route reads its sensitivities at the stage abscissae and `CY` at the
    stage states (2026-09-15).  Measured on the sampler: the LTI series sum
    equals the pnoise fold to 1e-10 under both methods; the held variance is
    1.000000 (radau) / 0.999979 (trbdf2) kT/C at 400 points; the tracking
    value converges at first order (radau 0.949 -> 0.975 at 400 -> 800).
    ⚠ Checked to matter: `CY` at the end-of-step state for every stage read
    the held variance 0.867 (radau, 400 points).  ⚠ `covariance` under these
    methods read 0.876 held at 400 points until its injection moved to the
    stages as well -- see
    `test_the_stage_method_covariance_injects_at_the_stages_and_holds_kTC_across_a_switching_edge`."""
    import warnings
    ktc = _KB * _TEMP / 100e-12
    for method in ('radau', 'trbdf2'):
        def els(c):
            c['S0'] = _sw(gon=1e-3, goff=1e-3)
            c['F0'] = _Flicker('out', gnd, i=0.0, noisePSD=1e-22, fref=1.0)
        cir, pss, io, pac, T = _sampler_fixture_method(els, method)
        f0 = 1.0 / T
        f = 0.137 * f0
        S = pac.sampled_noise(pss, io, [0.37 * T], [f], maxsidebands=10)[0, 0]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ref = sum(float(np.real(pac.pnoise(pss, abs(f + k * f0), io,
                                               maxsidebands=20)[0]))
                      for k in range(-10, 11))
        assert abs(S / ref - 1.0) < 1e-8, (method, S, ref)
    held = {}
    for method, npts in (('radau', 400), ('radau', 800), ('trbdf2', 400)):
        cir, pss, io, pac, T = _sampler_fixture_method(
            lambda c: c.__setitem__('S0', _sw()), method, npts)
        f0 = 1.0 / T
        grid = np.asarray(pss.factored_period().times, dtype=float)
        N = len(pss.factored_period().steps)
        fmin = 1e-6 * f0
        v = pac.sampled_variance(pss, io, [grid[int(0.375 * N)], grid[int(0.1 * N)]],
                                 fmin, 0.5 * f0, points_per_decade=10) / ktc
        held[(method, npts)] = v / (1.0 - fmin / (0.5 * f0))
    assert abs(held[('radau', 400)][0] - 1.0) < 1e-4, held
    assert abs(held[('trbdf2', 400)][0] - 1.0) < 1e-4, held
    assert abs(held[('radau', 800)][0] - 1.0) < 1e-4, held
    e4, e8 = 1.0 - held[('radau', 400)][1], 1.0 - held[('radau', 800)][1]
    assert e4 > 0.02 and 0.4 < e8 / e4 < 0.6, (e4, e8)


def test_the_sampled_noise_runs_natively_under_a_glm():
    """A Nordsieck GLM's map on the state (2026-09-25, `_GLMStateStep`):
    `sampled_noise` passes it as a stage map, one coupling per INJECTION
    POINT -- every stage, then the substages of the startup that opens a
    step -- with `CY` at the stored stage states.  Measured on the sampler:

    * the pass's couplings at `_stage_times`, seeded with an output at node
      1, reproduce the forced replay's response there to 1e-15 (exact: one
      set of costates; step 0 and its startup's substages carry it all);
    * the LTI series sum against the pnoise fold: glm2 1.9e-5 / 2.4e-6 /
      3.0e-7, glm3 2.1e-6 / 1.1e-7 / 1.6e-8 at 100 / 200 / 400 points
      (radau 6e-15) -- NOT a coupling defect: the startup at node 0 makes
      the discretised LTI circuit periodically time-varying, so an
      instant's samples carry O(h^p)-small sidebands the stationary fold
      does not, and the gap falls at the method's order;
    * the held variance glm2 0.999825 / 0.999978 / 0.999997, glm3 0.999775
      / 0.999973 / 0.999997 kT/C at 200 / 400 / 800 points; the tracking
      value first order, as radau's.

    Pinned: the identity to 1e-12, the LTI gap falling by more than 5x
    over 100 -> 200 (glm2), the held variance within 1e-4 at 400."""
    import warnings
    ktc = _KB * _TEMP / 100e-12
    for method in ('glm2', 'glm3'):
        cir, pss, io, pac, T = _sampler_fixture_method(
            lambda c: c.__setitem__('S0', _sw()), method, 100)
        fp = pss._state_map()
        assert fp.is_glm
        rng = np.random.default_rng(2)
        m = fp.width
        f = 0.137 / T
        tinj = pac._stage_times(pss, fp)
        ## the output at NODE 1 seeds the pass: step 0 and the substages of
        ## the startup that opens it carry the whole coupling (seeded at the
        ## period's end, the sampler damps it out before step 0 and a wrong
        ## substage time passed)
        d = rng.standard_normal(m)
        u = rng.standard_normal(m) + 1j * rng.standard_normal(m)
        _g, Cp = pac._stage_pass(pss, fp, np.zeros(m), (1, d))
        assert len(tinj) == len(Cp) == len(pac._stage_states(pss, fp))
        acc = -np.sum(np.exp(2j * np.pi * f * tinj)[:, None] * Cp, axis=0) @ u
        _e, ys = pss._forced_replay(fp, f, u, y0=np.zeros(m, dtype=complex), collect=True)
        assert abs(acc - d @ ys[0]) < 1e-12 * abs(d @ ys[0]), method
    gap = []
    for npts in (100, 200):
        def els(c):
            c['S0'] = _sw(gon=1e-3, goff=1e-3)
            c['F0'] = _Flicker('out', gnd, i=0.0, noisePSD=1e-22, fref=1.0)
        cir, pss, io, pac, T = _sampler_fixture_method(els, 'glm2', npts)
        f0 = 1.0 / T
        f = 0.137 * f0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            S = pac.sampled_noise(pss, io, [0.37 * T], [f], maxsidebands=10)[0, 0]
            ref = sum(float(np.real(pac.pnoise(pss, abs(f + k * f0), io,
                                               maxsidebands=20)[0]))
                      for k in range(-10, 11))
        gap.append(abs(S / ref - 1.0))
    assert gap[0] < 1e-4 and gap[0] / gap[1] > 5.0, gap
    for method in ('glm2', 'glm3'):
        cir, pss, io, pac, T = _sampler_fixture_method(
            lambda c: c.__setitem__('S0', _sw()), method, 400)
        f0 = 1.0 / T
        grid = np.asarray(pss.factored_period().times, dtype=float)
        N = len(pss.factored_period().steps)
        fmin = 1e-6 * f0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            v = pac.sampled_variance(pss, io, [grid[int(0.375 * N)]], fmin,
                                     0.5 * f0, points_per_decade=10) / ktc
        held = float(v[0]) / (1.0 - fmin / (0.5 * f0))
        assert abs(held - 1.0) < 1e-4, (method, held)


def test_the_time_average_of_the_sampled_psd_is_the_fold_of_time_averaged_pnoise():
    """The peer's validation 2, which holds PER FREQUENCY: averaged over the
    sampling instant, the sample series samples the time-averaged
    autocorrelation at lags kT, so `mean_t0 S(f; t0) = sum_k S_avg(|f + k f0|)`
    with `S_avg` = `pnoise(cyclostationary=True)` -- over EVERY output band to
    the grid's Nyquist.  A switched capacitor buffered into an RC low-pass
    (pole f0/20) so the output is cyclostationary and band-limited.
    Measured at 60 points: 3.9e-4 with the full fold (4.4e-5 at 100 points),
    while a fold cut at |k| <= 10 is 0.25 % off (0.41 % at 100) -- so the
    identity is checked, not a truncation that happens to agree."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    fclk = 100e3
    T = 1.0 / fclk
    f0 = fclk
    cir = SubCircuit()
    for nd in ('in', 'out', 'ck', 'b', 'y'):
        cir.add_node(nd)
    cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0)
    cir['Vck'] = VSin('ck', gnd, vo=0.0, va=1.0, freq=fclk, phase=90.0)
    cir['C0'] = C('out', gnd, c=100e-12)
    cir['S0'] = _sw()
    cir['E0'] = VCVS('out', gnd, 'b', gnd, g=1.0)
    cir['Rf'] = R('b', 'y', r=1e4)
    cir['Cf'] = C('y', gnd, c=1.0 / (2 * np.pi * 1e4 * f0 / 20))
    pss = PSS(cir, method='gear', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 60, x0=np.zeros(cir.n - 1),
                  maxiterations=100)
    assert pss.converged
    iy = [str(nd) for nd in cir.nodes if str(nd) != 'gnd!'].index('y')
    pac = PAC(cir, toolkit=circuit.numeric)
    N = len(pss.factored_period().steps)
    grid = np.asarray(pss.factored_period().times, dtype=float)[:N]
    f = 0.37 * f0
    mean = float(np.mean(pac.sampled_noise(pss, iy, grid, [f])[:, 0]))
    K = N // 2
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        parts = {k: float(np.real(pac.pnoise(pss, abs(f + k * f0), iy,
                                             cyclostationary=True)[0]))
                 for k in range(-K, K + 1)}
    full = sum(parts.values())
    cut = sum(v for k, v in parts.items() if abs(k) <= 10)
    assert abs(mean / full - 1.0) < 1e-3, (mean, full)
    assert abs(mean / cut - 1.0) > 2e-3, (mean, cut)


def test_the_sampled_psd_of_a_reset_rc_is_sepke_eq_33_white_and_one_over_f():
    """Sepke et al. 2009 (TCAS-I 56(3)) eq. (33), the LTI-plus-reset case, with
    no fitted factor: when the node is reset before a window of length `Tw`
    and filters LTI within it, the one-sided sample-series PSD is
    `S(f) = sum_n |H_w(f + n f0)|^2 CY(|f + n f0|)` with the WINDOWED impulse
    response `H_w(f) = (1 - exp(-(a + j2 pi f) Tw)) / (C (a + j2 pi f))`,
    `a = G/C`, and stationary `CY` = thermal `4kTG` plus a 1/f source.
    A 1 nF / 1 MOhm node, reset by a NOISELESS switch while the clock is
    above 0.95 (10 % of the period), sampled at two instants in the window.
    Measured at 400 points (radau): worst 1.2e-3 mid-window at 10 Hz,
    1.5e-4 white and 5e-4 1/f at the end of the window (1e-5 at 1000
    points).  The controls: the reference with the window one step longer is
    1e-3..1e-2 off, the never-reset (stationary) RC 0.13x..1.87x -- so the
    window, not only the filter, is what agrees."""
    import warnings
    from pycircuit.circuit.constants import kboltzmann
    circuit.default_toolkit = circuit.numeric
    f0 = 1e3
    T = 1.0 / f0
    Cv, G, vth, L, npts = 1e-9, 1e-6, 0.95, 50, 400
    fl_psd, fl_fc = 1e-24, 1e3
    cir = SubCircuit()
    for nd in ('out', 'ck'):
        cir.add_node(nd)
    cir['Vck'] = VSin('ck', gnd, vo=0.0, va=1.0, freq=f0, phase=0.0)
    cir['C0'] = C('out', gnd, c=Cv)
    cir['R0'] = R('out', gnd, r=1.0 / G)
    cir['SW'] = _SwitchHdl('out', gnd, 'ck', gnd, gon=1e-2, goff=1e-15,
                           vth=vth, vs=1e-5, temp=_TEMP, kb=0.0)
    cir['N0'] = IS('out', gnd, i=0.0, noisePSD=fl_psd, noiseFc=fl_fc)
    pss = PSS(cir, method='radau', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=np.zeros(cir.n - 1),
                  maxiterations=60)
    assert pss.converged
    io = [str(nd) for nd in cir.nodes if str(nd) != 'gnd!'].index('out')
    pac = PAC(cir, toolkit=circuit.numeric)
    fp = pss.factored_period()
    grid = np.asarray(fp.times, dtype=float)[:len(fp.steps)]
    half = np.arccos(vth) / (2 * np.pi) * T    # reset is centred on T/4
    t_open = T / 4 + half
    kT = kboltzmann * float(circuit.defaultepar.T)
    a = G / Cv
    fs = np.array([10.0, 50.0, 100.0, 200.0, 400.0, 500.0])

    def cy(nu):
        return 4 * kT * G + fl_psd * (1.0 + fl_fc / np.abs(nu))

    def fold(h):
        return np.array([np.sum(np.abs(h(f + np.arange(-L, L + 1) * f0)) ** 2
                                * cy(f + np.arange(-L, L + 1) * f0)) for f in fs])

    def windowed(Tw):
        return lambda nu: ((1 - np.exp(-(a + 2j * np.pi * nu) * Tw))
                           / (Cv * (a + 2j * np.pi * nu)))

    for back in (0.01, 0.45):          # end of the window, and mid-window
        t0 = grid[np.argmin(np.abs(grid - (T / 4 - half - back * T) % T))]
        Tw = t0 - t_open if t0 > t_open else t0 + T - t_open
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            S = pac.sampled_noise(pss, io, [t0], fs, maxsidebands=L)[0]
        ref = fold(windowed(Tw))
        assert np.max(np.abs(S / ref - 1.0)) < 3e-3, (back, S / ref)
        assert np.min(np.abs(S / fold(windowed(Tw - T / npts)) - 1.0)) > 5e-4, back
        stationary = fold(lambda nu: 1.0 / (Cv * (a + 2j * np.pi * nu)))
        assert np.min(np.abs(S / stationary - 1.0)) > 0.1, (back, S / stationary)


def test_a_coloured_source_keeps_the_sign_of_its_scale_factor_through_the_periodic_folds():
    """⚠⚠ THE SIGN-BLIND sqrt(PSD) FOLD, FIXED WHERE THE ELEMENT STATES THE SIGN.

    A coloured source is correlated ACROSS the period, so `k(x(t)) n(t)` has
    `R(t,t') = k(t) k(t') R_n(t-t')`, and `CY = k^2 S` has lost the sign of `k`.
    The folds factored `CY` with its square root and so computed the |k|
    process.  The exact answer is available with nothing fitted: a CONSTANT
    1/f source through a signed multiplier ('A', a stationary fold) is the same
    physics as the modulated source.  pnoise at 1e-3 f0, LO = vo + sin:

        vo      sqrt(PSD) / exact     signed amplitude / exact
        1.5        1.000                1.000000000     (sign-definite: no defect)
        0.5        2.059                1.000000000
        0.0      810.6                  1.000000000     (zero-mean LO: all of it
                                                         is up-converted away)

    This is the PSP sampler's "+0.1 % flicker residual" (2026-09-17..19): Vds
    crosses zero while the switch conducts, the 1/f current follows sgn(Vds),
    and the residual survived every grid, sideband and tolerance knob because
    it was never numerical.  A commercial simulator shows a NOTCH in the
    sampled flicker at one clock amplitude, where the signed modulation
    cancels; a sign-blind fold cannot produce one.

    The contract: the scale factor OUTSIDE the noise call is the amplitude and
    carries the sign (`Element.noise_amplitudes`); a sign squared into the
    power is gone, and such a source keeps the |k| fold AND its warning.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    f0 = 1.0 / T
    k1, Rn = 0.5, 2.0

    def run(kind, vo, sampled=False):
        c = SubCircuit()
        c.add_node('lo')
        c.add_node('out')
        c['vlo'] = VSin('lo', gnd, va=1.0, vo=vo, freq=f0)
        if kind == 'A':
            c.add_node('n')
            c['xi'] = _Flicker('n', gnd, i=0.0, noisePSD=1.0, fref=1.0)
            c['Rn'] = R('n', gnd, r=Rn)
            c['M1'] = _NuMult('out', gnd, 'n', gnd, 'lo', gnd, k=k1)
        else:
            cls = _SgnAmpFlicker if kind == 'C' else _SgnPsdFlicker
            c['src'] = cls('out', gnd, 'lo', gnd, k=k1 * Rn)
        c['Ro'] = R('out', gnd, r=1.0)
        c['Co'] = C('out', gnd, c=0.2e-6)
        pss = PSS(c, method='gear', reltol=1e-10)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            pss.solve(period=T, timestep=T / 200, maxiterations=40)
            o = [str(x) for x in c.nodes].index('out')
            pac = PAC(c, toolkit=circuit.numeric)
            if sampled:
                s = pac.sampled_noise(pss, o, [0.3 * T], np.array([1e-3 * f0]),
                                      maxsidebands=16)[0, 0]
            else:
                s, _ = pac.pnoise(pss, 1e-3 * f0, o, maxsidebands=16,
                                  cyclostationary=(kind != 'A'))
        signwarn = [r for r in rec if 'Only the element knows the sign'
                    in str(r.message)]
        return float(np.real(s)), bool(signwarn)

    for vo in (0.5, 0.0):
        a, _wa = run('A', vo)
        b, wb = run('B', vo)
        c_, wc = run('C', vo)
        assert abs(c_ / a - 1.0) < 1e-6, (vo, c_ / a)
        assert wb and not wc, (vo, wb, wc)
        ## presence: the sign-blind fold is wrong by a FACTOR here
        assert b / a > (500.0 if vo == 0.0 else 1.8), (vo, b / a)
    ## sign-definite: nothing to lose, all three agree
    a, _w0 = run('A', 1.5)
    assert abs(run('B', 1.5)[0] / a - 1.0) < 1e-6
    assert abs(run('C', 1.5)[0] / a - 1.0) < 1e-6
    ## the time-sampled fold takes the same amplitudes: the two specifications
    ## agree where the sign is definite and part by a factor where it is not
    sb, sc = run('B', 1.5, sampled=True)[0], run('C', 1.5, sampled=True)[0]
    assert abs(sb / sc - 1.0) < 1e-6, sb / sc
    (sb, wsb), (sc, wsc) = run('B', 0.0, sampled=True), run('C', 0.0, sampled=True)
    assert sb / sc > 1.5, sb / sc
    ## and the sample series WARNS the |k| fold as the periodic one does
    ## (it was silent until 2026-09-29)
    assert wsb and not wsc, (wsb, wsc)

    ## the element-level contract
    el = _SgnAmpFlicker('out', gnd, 'lo', gnd, k=2.0)
    w = 2 * np.pi * 10.0
    for v in (0.7, -0.7):
        x = np.zeros(el.n)
        x[2] = v                                   # terminal `b`
        W = el.noise_amplitudes(x, w)
        assert W.shape == (el.n, 1)
        np.testing.assert_allclose(W @ W.conj().T, el.CY(x, w), rtol=1e-12,
                                   atol=1e-30)
        assert np.sign(np.real(W[0, 0])) == np.sign(v)
    assert _NuModNoise('out', gnd, 'lo', gnd).noise_amplitudes(
        np.zeros(4), w) is None                    # white only: nothing to state
    ## ⚠ and a sign squared into the POWER is not stated either: sqrt(pwr)
    ## would pass the |k| process off as signed and silence the warning
    assert _SgnPsdFlicker('out', gnd, 'lo', gnd).noise_amplitudes(
        np.zeros(4), w) is None


def _lorentz_lo_rc(kind, vo=0.0, T=1e-6):
    """A driven RC whose noise is a Lorentzian (not a power law, so it is
    evaluated per band) times V(lo), ``lo = vo + sin`` (crossing zero at
    `vo = 0`): as the ELEMENT (`kind` 'signed' or 'psd'), or REALISED -- a
    stationary Lorentzian through a multiplier by V(lo), the sign in the
    circuit, exact for every analysis.  Returns `(cir, pss, pac, out)`."""
    P, tau = 1e-20, 3e-7
    c = SubCircuit()
    for nd in ('lo', 'out'):
        c.add_node(nd)
    c['Vlo'] = VSin('lo', gnd, va=1.0, vo=vo, freq=1.0 / T)
    c['Ro'] = R('out', gnd, r=1e3, noisy=False)
    c['Co'] = C('out', gnd, c=0.5e-9)
    if kind == 'real':
        c.add_node('n')
        c['xi'] = IS('n', gnd, i=0.0, noisePSD=P, noiseTau=tau)
        c['Rn'] = R('n', gnd, r=1.0, noisy=False)
        c['M'] = _NuMult('out', gnd, 'n', gnd, 'lo', gnd, k=1.0)
    else:
        cls = _ModLorentzSigned if kind == 'signed' else _ModLorentzCtl
        c['src'] = cls('out', gnd, 'lo', gnd, noisePSD=P, tau=tau, k=1.0)
    pss = PSS(c, method='radau', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 100, maxiterations=40)
    return c, pss, PAC(c, toolkit=circuit.numeric), \
        [str(n) for n in c.nodes].index('out')


def test_the_sampler_takes_a_per_band_elements_signed_amplitudes():
    """The sample series factors a per-band (not power-law) source by the
    element's SIGNED amplitudes where it states them, as `covariance`
    already did.  ⚠ It took the root of the PSD -- the |m| process -- so
    with the modulation crossing zero `sampled_variance` read +115 % /
    +209 % of the exact answer at 0.3 T / 0.7 T, while `covariance` on the
    same element was exact (measured 2026-09-27; Andreas: "the sampler
    should use the element's signed amplitudes").  A PSD-only element keeps
    the |m| answer: the sign is not in it."""
    T = 1e-6
    sv = {}
    for kind in ('real', 'signed', 'psd'):
        _c, pss, pac, o = _lorentz_lo_rc(kind, T=T)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sv[kind] = np.asarray(pac.sampled_variance(
                pss, o, [0.3 * T, 0.7 * T], 1e-4 / T, 0.5 / T,
                points_per_decade=10), dtype=float).ravel()
    assert np.max(np.abs(sv['signed'] / sv['real'] - 1.0)) < 1e-9, sv
    assert np.min(sv['psd'] / sv['real'] - 1.0) > 1.0, sv


def test_a_moving_shape_takes_the_elements_signed_columns_on_every_surface():
    """A per-band element whose spectral SHAPE moves along the orbit -- two
    independent Lorentzians under two modulations, one crossing zero
    (`_TwoModLorentz`) -- is read per band frequency from its SIGNED
    columns in `sampled_variance`, `covariance` and the cyclostationary
    `pnoise` (Andreas, 2026-09-27: "Can covariance's own moving-shape path
    ... be fixed?").  ⚠ All three took one root of the element's PSD per
    point: SIGN-BLIND, and the two independent sources merged into one
    column.  Measured against the realisation (each Lorentzian through its
    own multiplier): +9 % / +8 % (+18 % / +12 % in pnoise) with the
    modulation crossing zero, and still +2.5 % / +5.8 % with both
    positive, where only the merging is wrong.  The PSD-only element keeps
    that answer (bit for bit: its path is unchanged).  The band is cut at
    1e-2 f0 to keep the test short; the realisation is exact on any band."""
    T = 1e-6
    P1, tau1, P2, tau2 = 1e-20, 3e-7, 2e-20, 3e-8

    def build(kind):
        c = SubCircuit()
        for nd in ('lo', 'lo2', 'out'):
            c.add_node(nd)
        c['Vlo'] = VSin('lo', gnd, va=1.0, vo=0.0, freq=1.0 / T)
        c['Vlo2'] = VSin('lo2', gnd, va=0.5, vo=1.0, freq=2.0 / T)
        c['Ro'] = R('out', gnd, r=1e3, noisy=False)
        c['Co'] = C('out', gnd, c=0.5e-9)
        if kind == 'real':
            for i, (p, tau, lo) in enumerate(((P1, tau1, 'lo'), (P2, tau2, 'lo2'))):
                n = 'n%d' % i
                c.add_node(n)
                c['xi%d' % i] = IS(n, gnd, i=0.0, noisePSD=p, noiseTau=tau)
                c['Rn%d' % i] = R(n, gnd, r=1.0, noisy=False)
                c['M%d' % i] = _NuMult('out', gnd, n, gnd, lo, gnd, k=1.0)
        else:
            c['src'] = _TwoModLorentz('out', gnd, 'lo', gnd, 'lo2', gnd,
                                      P1=P1, tau1=tau1, P2=P2, tau2=tau2)
            c['src'].signed = kind == 'signed'
        pss = PSS(c, method='radau', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 100, maxiterations=40)
        pac = PAC(c, toolkit=circuit.numeric)
        o = [str(n) for n in c.nodes].index('out')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sv = np.asarray(pac.sampled_variance(
                pss, o, [0.3 * T, 0.7 * T], 1e-2 / T, 0.5 / T,
                points_per_decade=10), dtype=float).ravel()
            _K0, Ks = pac.covariance(pss, samples=True, colour_fmin=1e-2 / T,
                                     points_per_decade=10)
            pn = np.array([float(np.real(pac.pnoise(
                pss, fr / T, o, maxsidebands=30, cyclostationary=True)[0]))
                for fr in (0.013, 0.31)])
        tms = np.asarray(pss.factored_period().times, dtype=float)
        cov = np.array([Ks[int(np.argmin(np.abs(tms[:len(Ks)] - t)))][o, o]
                        for t in (0.3 * T, 0.7 * T)])
        return np.concatenate((sv, cov, pn))
    exact, signed = build('real'), build('signed')
    assert np.max(np.abs(signed / exact - 1.0)) < 1e-9, signed / exact - 1.0


def test_signed_columns_beside_a_white_source_of_the_same_element():
    """A per-band element that states signed amplitudes for its COLOURED
    source and none for a white one beside it (`_ModLorentzThermal`: an HDL
    device with thermal and G-R noise states its coloured sources only) is
    read from its signed columns plus the root of the remainder ``CY - W
    W^H`` (`_perband_mode` 'white'); in `covariance` that remainder joins
    the white part, whose Lyapunov path covers every frequency.  ⚠ The
    amplitudes did not rebuild the CY, so the element fell back to the root
    of its PSD: +54 .. +110 % against the realisation with the Lorentzian
    crossing zero (2026-09-28); with the remainder band-limited in
    `covariance`, -1.4 % / -1.8 %."""
    T, P, tau, Pw = 1e-6, 1e-20, 3e-7, 2e-21

    def build(kind):
        c = SubCircuit()
        for nd in ('lo', 'out'):
            c.add_node(nd)
        c['Vlo'] = VSin('lo', gnd, va=1.0, vo=0.0, freq=1.0 / T)
        c['Ro'] = R('out', gnd, r=1e3, noisy=False)
        c['Co'] = C('out', gnd, c=0.5e-9)
        if kind == 'real':
            c.add_node('n')
            c['xi'] = IS('n', gnd, i=0.0, noisePSD=P, noiseTau=tau)
            c['Rn'] = R('n', gnd, r=1.0, noisy=False)
            c['M'] = _NuMult('out', gnd, 'n', gnd, 'lo', gnd, k=1.0)
            c['xw'] = IS('out', gnd, i=0.0, noisePSD=Pw)
        else:
            c['src'] = _ModLorentzThermal('out', gnd, 'lo', gnd, noisePSD=P,
                                          tau=tau, k=1.0, white=Pw)
        pss = PSS(c, method='radau', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 100, maxiterations=40)
        pac = PAC(c, toolkit=circuit.numeric)
        o = [str(n) for n in c.nodes].index('out')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sv = np.asarray(pac.sampled_variance(
                pss, o, [0.3 * T, 0.7 * T], 1e-2 / T, 0.5 / T,
                points_per_decade=10), dtype=float).ravel()
            _K0, Ks = pac.covariance(pss, samples=True, colour_fmin=1e-2 / T,
                                     points_per_decade=10)
            pn = np.array([float(np.real(pac.pnoise(
                pss, fr / T, o, maxsidebands=30, cyclostationary=True)[0]))
                for fr in (0.013, 0.31)])
        tms = np.asarray(pss.factored_period().times, dtype=float)
        cov = np.array([Ks[int(np.argmin(np.abs(tms[:len(Ks)] - t)))][o, o]
                        for t in (0.3 * T, 0.7 * T)])
        return np.concatenate((sv, cov, pn))
    exact, elem = build('real'), build('element')
    assert np.max(np.abs(elem / exact - 1.0)) < 1e-9, elem / exact - 1.0


def test_signed_columns_of_a_mixed_slope_component_are_grouped_by_exponent():
    """A power-law component whose exponent differs between its entries --
    two 1/f sources of slope 0.8 and 2.0 on two branches of one element,
    each times V(lo) crossing zero (`_MixedSlopeLo`) -- is taken by its
    SIGNED columns grouped by their own exponents (`_exponent_columns`),
    a uniform power law per group, in `sampled_variance`, `covariance` and
    the cyclostationary `pnoise`.  ⚠ It was rooted per band from its PSD,
    warned as sign-blind: 4x to 4800x the realisation (each slope
    stationary, through its own multiplier) (2026-09-28)."""
    T, k = 1e-6, 1e-20

    def build(kind):
        c = SubCircuit()
        for nd in ('lo', 'a', 'b'):
            c.add_node(nd)
        c['Vlo'] = VSin('lo', gnd, va=1.0, vo=0.0, freq=1.0 / T)
        for nd in ('a', 'b'):
            c['R' + nd] = R(nd, gnd, r=1e3, noisy=False)
            c['C' + nd] = C(nd, gnd, c=1e-9)
        if kind == 'real':
            for nd, cls in (('a', _Flicker08), ('b', _Flicker20)):
                c.add_node('n' + nd)
                c['F' + nd] = cls('n' + nd, gnd, k=k)
                c['Rn' + nd] = R('n' + nd, gnd, r=1.0, noisy=False)
                c['M' + nd] = _NuMult(nd, gnd, 'n' + nd, gnd, 'lo', gnd, k=1.0)
        else:
            c['src'] = _MixedSlopeLo('a', gnd, 'b', gnd, 'lo', gnd, k=k)
        pss = PSS(c, method='radau', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 100, maxiterations=40)
        pac = PAC(c, toolkit=circuit.numeric)
        names = [str(n) for n in c.nodes]
        oa, ob = names.index('a'), names.index('b')
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            sv = [float(np.ravel(pac.sampled_variance(
                pss, o, [0.3 * T], 1e-2 / T, 0.5 / T, points_per_decade=10))[0])
                for o in (oa, ob)]
            _K0, Ks = pac.covariance(pss, samples=True, colour_fmin=1e-2 / T,
                                     points_per_decade=10)
            pn = [float(np.real(pac.pnoise(pss, 0.013 / T, o, maxsidebands=30,
                                           cyclostationary=True)[0]))
                  for o in (oa, ob)]
        tms = np.asarray(pss.factored_period().times, dtype=float)
        j = int(np.argmin(np.abs(tms[:len(Ks)] - 0.3 * T)))
        blind = [w for w in caught if 'SIGN-BLIND' in str(w.message)]
        return np.array(sv + [Ks[j][oa, oa], Ks[j][ob, ob]] + pn), blind
    (exact, _), (elem, blind) = build('real'), build('element')
    assert not blind, [str(w.message)[:80] for w in blind]
    assert np.max(np.abs(elem / exact - 1.0)) < 1e-9, elem / exact - 1.0


def test_the_psp_sampled_flicker_has_the_notch_a_sign_blind_fold_cannot_produce():
    """⚠⚠ THE PSP SAMPLER'S FLICKER RESIDUAL, CLOSED: the 1/f current follows
    sgn(Vds), and Vds crosses zero while the switch conducts.

    `PspMosLongChannel` now contributes `sgn * flicker_noise(n_sfl)` -- the
    sign in the AMPLITUDE, `CY` unchanged -- and the sampled fold takes it.
    Swept over the clock amplitude (10 Hz, sampled mid-hold, 200 points), the
    sign-blind fold over the signed one, beside the ratio a peer session
    measured for the OLD fold against a commercial simulator:

        clock amp    0.75    0.50    0.375   0.30   0.25   0.20   0.15   0.10
        blind/signed 1.0007  1.0069  1.0602  1.420  3.585  401    5.16   2.18
        old/reference 1.0011 1.0073  1.0582  1.387  3.238  567    5.87   2.27

    so the campaign's "+0.1 %" at the fixture's 0.75 -- which survived every
    grid, sideband, tolerance and frequency-axis knob -- was never numerical.
    At 0.20 the signed modulation CANCELS through the sampler's transfer: a
    notch, reference-converged, which no fold built from a PSD can produce.
    Gated here as the notch itself (coarse grid, so as a factor).

    ⚠ One more thing the sweep found: at 0.375 alone the two folds read
    bit-identical -- one orbit sample at Vds ~ 0 carried a flicker entry at
    1e-12 of scale whose fitted exponent was rounding, and that failed the
    whole component into the per-band route.  See `uniform_exponent`.
    """
    import os
    import warnings as _w
    PDK = os.path.expanduser(
        '~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models')
    if not os.path.isdir(PDK):
        pytest.skip('IHP Open PDK not installed')
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit import psp_scaling, defaultepar
    from pycircuit.circuit.compact import PspMosLongChannel
    from pycircuit.utilities import spicecard
    T27 = 273.15 + 27.0
    was = defaultepar.T
    defaultepar.T = T27
    try:
        circuit.default_toolkit = circuit.numeric
        deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                              section='mos_tt')
        w, l = 10e-6, 1e-6
        base = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1, m=1,
                              pre_layout=1), w=w, l=l, T=T27)
        F = 1e5

        def density(amp):
            cir = SubCircuit()
            cir['Vin'] = VSin('in', gnd, vo=0.3, va=0.2, freq=F, phase=0.0)
            cir['Vck'] = VSin('ck', gnd, vo=0.75, va=amp, freq=F, phase=90.0)
            cir['M1'] = PspMosLongChannel(cm.Node('in'), cm.Node('ck'),
                                          cm.Node('out'), gnd,
                                          **dict(base, swign=0.0))
            cir['Ch'] = C('out', gnd, c=100e-12)
            pss = PSS(cir, method='gear', reltol=1e-10)
            out = {}
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                pss.solve(period=1 / F, timestep=1 / F / 100, maxiterations=60)
                assert pss.converged
                k = cir.get_node_index(cir.get_node('out'))
                k = k - 1 if k > pss.irefnode else k
                tt = np.asarray(pss.factored_period().times, float)
                tq = tt[np.argmin(abs(tt - 5e-6))]
                for tag in ('signed', 'blind'):
                    pac = PAC(cir, toolkit=circuit.numeric)
                    if tag == 'blind':
                        _noise_seam(pac, signed_amplitudes=lambda *a_, **k_: {})
                    out[tag] = float(np.real(pac.sampled_noise(
                        pss, k, [tq], np.array([10.0]), maxsidebands=40))[0, 0])
                ## a FRESH PAC: `pac` above is the deliberately blinded one
                model = PAC(cir, toolkit=circuit.numeric)._noise_components(
                    pss).model(10.0, F)
            return out, model

        ## the element states its sign, and its amplitudes rebuild the flicker
        hi, model = density(0.25)
        assert [key for key, _B, _E in model.flicker] == [('M1',)]
        assert ('M1',) in model.amplitude
        lo, _m = density(0.20)
        ## the notch: the SIGNED density collapses, the blind one does not
        assert hi['signed'] / lo['signed'] > 10.0, (hi, lo)
        assert lo['blind'] / lo['signed'] > 50.0, (hi, lo)
        assert 0.2 < hi['blind'] / lo['blind'] < 5.0, (hi, lo)
    finally:
        defaultepar.T = was


def test_the_sampled_noise_tail_closure_is_derived_and_the_resolution_limit_is_named():
    """Two limits on a held variance, told apart (2026-09-21, from a peer's
    kT/C ladder): the source spectrum beyond the covered sidebands (a
    Lorentzian's tail `1 - (2/pi) atan(F/fc)`, coefficient ONE -- the fold's
    `|n| <= L` sum is symmetric and covers |nu| < (L + 1/2) f0 whole), and
    the discretisation of the covered top sidebands at `omega h` per step.
    The peer's "coefficient 3.0" was the second at a fixed M/npts: at a
    fixed 100 sidebands the coefficient against the pure tail fell 3.03 /
    1.93 / 1.33 / 1.12 / 1.05 at 204 .. 3200 points.

    LTI limit of the switched capacitor (vck = 0: g = 0.5 mS exactly, C =
    100 pF, fc = 795.78 kHz, f0 = 100 kHz), 100 sidebands (edge 12 fc)::

        method  points  omega h   no tail    tail=True   warned
        radau   204     3.1       0.9479     0.9981      yes
        gear    204     3.1       0.8469     0.8737      yes
        gear    800     0.79      0.9326     0.9688      no

    Radau's covered bands are accurate even at omega h = 3, so its deficit
    IS the pure tail (5.2e-2 against 5.06e-2 exact) and the closure -- the
    1/nu^2 extrapolation of the two outermost covered sidebands, nothing
    fitted -- removes it to 2e-3.  Gear's deficit at 204 points is mostly
    discretisation, which no closure can remove; the warning names it.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    fclk, cval, kb, temp = 100e3, 100e-12, 1.38e-23, 300.0
    ktc = kb * temp / cval
    T = 1.0 / fclk
    gon, goff = 1e-3, 1e-9

    def build():
        cir = SubCircuit()
        cir.add_node('in')
        cir.add_node('out')
        cir.add_node('ck')
        cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0)
        cir['Vck'] = VSin('ck', gnd, vo=0.0, va=0.0, freq=fclk, phase=90.0)
        cir['S0'] = _SwitchHdl('in', 'out', 'ck', gnd, gon=gon, goff=goff,
                               vth=0.0, vs=50e-3, temp=temp, kb=kb)
        cir['C0'] = C('out', gnd, c=cval)
        return cir

    def held(method, npts, tail):
        cir = build()
        p = PSS(cir, method=method, reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T, timestep=T / npts, x0=np.zeros(cir.n - 1),
                    maxiterations=40)
        assert p.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        out = [str(n_) for n_ in cir.nodes].index('out')
        out = out if out < p.irefnode else out - 1
        tms = np.asarray(p.factored_period().times, float)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            v = float(pac.sampled_variance(p, out, [tms[npts // 2]],
                                           1e-3 * fclk, 0.5 * fclk,
                                           points_per_decade=40,
                                           maxsidebands=100, tail=tail)[0])
        warned = any('omega h' in str(w_.message) for w_ in rec)
        return v / ktc, warned
    fc = (goff + (gon - goff) / 2) / (2 * np.pi * cval)
    pure_tail = 1.0 - (2.0 / np.pi) * np.arctan(100.5 * fclk / fc)
    r0, w0 = held('radau', 204, False)
    r1, w1 = held('radau', 204, True)
    assert abs((1.0 - r0) - pure_tail) < 0.3 * pure_tail, (1.0 - r0, pure_tail)
    assert 1.0 - r1 < 5e-3, r1
    assert w0 and w1
    g0, gw0 = held('gear', 204, False)
    g8, gw8 = held('gear', 800, False)
    g8t, _ = held('gear', 800, True)
    assert gw0 and not gw8, (gw0, gw8)
    assert g0 < 0.87 and 0.92 < g8 < 0.95 and 0.96 < g8t < 0.98, (g0, g8, g8t)


def test_gears_euler_backstop_step_on_an_event_grid_replays_in_covariance_and_the_adjoint():
    """A one-step companion is a two-step companion with a zero third
    coefficient (2026-09-21, found by the driven check of the fold).
    `event_grid` lands a T/200 clock ramp inside a T/63 cell; the sliver
    from the ramp's end to the next node grows 5.4x into the following
    cell and `Gear2Integrator.check_order_drop` takes THAT step at order 1
    past the zero-stability bound -- `alphas = (1/h, -1/h)`.  The
    solved-history consumers read `alphas[2]` unguarded
    (`PAC.covariance`: `IndexError`) or refused the step as "unreachable
    through solve" (`_monodromy_matvec_transposed`: every adjoint noise
    call).  Pinned: the Euler step IS in the factored period, the forward
    covariance and the reverse `sampled_variance` -- the two recursions
    the fix touches, run independently -- agree at the hold instant to
    2e-3 where the sampled sum is resolved (N = 252; at 63 its 31
    sidebands truncate it by 5.6 %), and the held variance climbs the
    recursion's own O(h/tau)
    tracking floor with N (0.686 / 0.740 / 0.902 x kT/C at 63 / 126 / 252
    measured; the floor is the recorded gear covariance item, not this).
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-5
    ktc = 1.38e-23 * 300.0 / 100e-12

    def run(npts):
        cir = _pulse_clocked_sampler(T)
        p = PSS(cir, method='gear', reltol=1e-8)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T, timestep=T / npts, x0=np.zeros(cir.n - 1),
                    maxiterations=100)
        assert p.converged
        fp = p.factored_period()
        n_euler = sum(1 for st in fp.steps if len(st[2]) == 2)
        io = [str(n_) for n_ in cir.nodes].index('out')
        io = io if io < p.irefnode else io - 1
        pac = PAC(cir, toolkit=circuit.numeric)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            _K0, Ks = pac.covariance(p, samples=True)
            sv = pac.sampled_variance(p, io, [0.8 * T], 1.0, 50e3)
        ts = np.asarray(p.waveform[0], float)
        v = np.array([np.asarray(k, float)[io, io] for k in Ks]) / ktc
        tt = ts[:len(v)]
        at = float(np.interp(0.8 * T, tt, v))
        held = float(v[(tt > 0.6 * T) & (tt < 0.95 * T)].mean())
        return n_euler, held, at, float(np.asarray(sv, float).ravel()[0]) / ktc

    n63, h63, a63, s63 = run(63)
    ## (4 since 2026-09-28: the landed clock edges drop the order too, so
    ## the backstop step has three more Euler steps beside it in the ring)
    assert n63 >= 1, 'the Euler step must be in the ring for this to test anything: %d' % n63
    assert np.isfinite(s63) and np.isfinite(a63)      # both recursions RUN across it
    ## ⚠ the reverse-pass agreement is pinned at 252, not 63: at 63 the
    ## sampled sum covers 31 sidebands and reads 0.647 against the
    ## covariance's 0.686 -- the recorded sampled-resolution item, not
    ## the replay (0.742 / 0.740 at 126, 0.9023 / 0.9018 at 252)
    n252, h252, a252, s252 = run(252)
    assert n252 >= 1, n252
    assert abs(s252 / a252 - 1.0) < 2e-3, (s252, a252)
    ## the O(h/tau) floor, climbing.  ⚠ 0.570 / 0.781 since 2026-09-28 (were
    ## 0.686 / 0.902): the Euler steps at the landed clock edges damp the
    ## held node further -- the order drop's cost on a smooth state, here in
    ## a noise number (kT/C is 1)
    assert 0.5 < h63 < h252 < 1.0, (h63, h252)


def test_a_coloured_covariance_meets_the_sampled_variance_sign_included():
    """Coloured `covariance` against `sampled_variance` on the switched
    sampler (2026-09-25): two routes that share no integration -- here the
    FORWARD response to each input frequency over `[fmin, Nyquist]`, there
    the seeded ADJOINT and the sample series' fold over `[fmin, f0/2]`
    (they cover the same input frequencies up to holes of width `2 fmin`
    around each clock harmonic, 1e-6 f0 wide here).  The source is a 1/f
    current into the held node MODULATED by the clock, `k V_ck zeta`, whose
    sign flips twice per period; the switch is noiseless (`kb = 0`).

    Measured at 200 points, the hold instant (0.6 T), `sampled_variance` at
    its default 40 per decade: amplitude-stated (`_SgnAmpFlicker`, the
    signed process) K/sv - 1 = -3.2e-5, PSD-stated (`_SgnPsdFlicker`, the
    |m| process) the same -- `sampled_variance`'s own quadrature, second
    order in the grid ratio (-5.0e-6 at 100 per decade).  Under the linear
    trapezoid it used until 2026-09-25 the same comparison read -8.3e-5 at
    100 per decade, (r - 1)^3 / 6 per point on a 1/f band.  The two
    processes differ by 0.9 % at the hold
    (the sampled charge from before the edge, where the clock is positive,
    against the charge integrated after it, where it is negative), and
    both routes carry that difference: a sign-blind `sqrt(B)` amplitude,
    or a modulation frozen at one state, fails here.  ⚠ And both WARN the
    PSD-stated source (the |m| process where the clock changes sign), and
    neither the signed one: until 2026-09-29 both rooted it silently."""
    import warnings
    out = {}
    for kind, cls in (('amp', _SgnAmpFlicker), ('psd', _SgnPsdFlicker)):
        def elements(cir, cls=cls):
            cir['S'] = _sw(kb=0.0)
            cir['F'] = cls('out', gnd, 'ck', gnd, k=1e-11)
        cir, pss, io, pac, T = _sampler_fixture(elements, npts=200)
        f0 = 1.0 / T
        N = len(pss.factored_period().steps)
        jh = int(0.6 * N)
        th = float(np.asarray(pss.factored_period().times)[jh])
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter('always')
            _K0, seq = pac.covariance(pss, samples=True, colour_fmin=1e-6 * f0)
            sv = pac.sampled_variance(pss, io, [th], 1e-6 * f0, 0.5 * f0)[0]
        blind = sorted({str(r.message).split(':')[0] for r in rec
                        if 'touches zero along the orbit' in str(r.message)})
        assert blind == ([] if kind == 'amp' else
                         ['PAC.covariance', 'PAC.sampled_noise']), (kind, blind)
        out[kind] = seq[jh][io, io]
        assert abs(out[kind] / sv - 1.0 + 3.2e-5) < 1e-5, (kind, out[kind] / sv - 1)
    assert 5e-3 < abs(out['amp'] / out['psd'] - 1.0) < 2e-2, out


def test_a_coloured_covariance_takes_noise_correlated_across_elements():
    """The band integral (`covariance` and its kin) on the correlated
    sources of `_mixed_exponent_rc('corr')` built as two elements and a
    `CY` override: identical to the one-element build, which is pinned to
    the closed forms (refused until 2026-09-26).  The sample series, which
    took one root of the whole `CY` per band before, now takes the joint
    components: the same here (2.2e-16, measured)."""
    import warnings as _w
    K, V = {}, {}
    for kind in ('corr', 'xcorr'):
        pss, pac, names, (T, _Rv, _Cv, _k) = _mixed_exponent_rc(kind)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            K[kind] = pac.covariance(pss, colour_fmin=1e3)
            V[kind] = float(pac.sampled_variance(
                pss, names.index('o2'), [0.3 * T], 0.01 / T, 0.5 / T,
                points_per_decade=10)[0])
    assert np.max(np.abs(K['xcorr'] - K['corr'])) <= \
        1e-12 * np.max(np.abs(K['corr'])), 'the joint split is not exact'
    assert abs(V['xcorr'] / V['corr'] - 1.0) <= 1e-12, V


def test_the_sampled_series_sees_the_crossings_motion_on_a_staged_solve():
    """Item 3 of the 2026-09-22 list (2026-09-23): `sampled_noise` sits on
    the same adjoint machinery as the sideband row, and on a staged solve
    it now uses the TOTAL map, reads its sample at FIXED time and carries
    the crossings' costate injections in a third reverse pass.

    The verification is the point and it is sharp on `_jitter_sampler`:
    the HELD node has no noise source of its own -- `Ch` is noiseless, the
    switch is off, the ramps are ideal sources -- so every part of its
    sample-series variance arrives through the crossing's motion.
    Measured at `t = 0.7 T` over `[1e-3, 0.5] f0` (radau): 1.373e-8 at 100
    points and 1.500e-8 at 200 against the analytic `(s2/s1)^2 kT/C_n` =
    1.657e-8 -- 0.83 and 0.91, converging as the band the series covers
    does.  ⚠ The BORDERING itself is worth 1.15e-6 of that on radau (the
    landed grid's per-step maps already carry the switching, as
    everywhere with the compact transition) and 8.2e-3 on gear, whose
    two-step map does not -- measured both ways, which is why the code
    stays: the consumers now agree on one linearisation.  Gear's absolute
    level, 0.63 / 0.75, is its own covariance floor (see the gear test)."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    f0 = 1.0 / T
    kT = 1.380649e-23 * 300.0
    exp_h = 4.0 * kT / 1e-12
    got = {}
    for N in (100, 200):
        cir = _jitter_sampler(T)
        p = PSS(cir, method='radau', reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T, timestep=T / N, maxiterations=100, state_events=True)
        assert p.converged and p._event_columns is not None
        io = [str(n_) for n_ in cir.nodes].index('hold')
        pac = PAC(cir, toolkit=circuit.numeric)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            V = pac.sampled_variance(p, io, [0.7 * T], 1e-3 * f0, 0.5 * f0,
                                     points_per_decade=12)
            ev = p._event_columns
            p._event_columns = None
            try:
                V0 = pac.sampled_variance(p, io, [0.7 * T], 1e-3 * f0, 0.5 * f0,
                                          points_per_decade=12)
            finally:
                p._event_columns = ev
        got[N] = float(V[0]) / exp_h
        ## every part of it came through the crossing: a held, noiseless node
        assert 0.75 < got[N] < 1.02, (N, got[N])
        ## the bordering is live and small on radau's landed grid
        rel = abs(float(V[0]) - float(V0[0])) / float(V0[0])
        assert 1e-8 < rel < 1e-4, (N, rel)
    assert got[200] > got[100], got
