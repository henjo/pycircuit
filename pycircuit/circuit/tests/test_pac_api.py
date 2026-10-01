"""The PAC surfaces' contract (the review's batch 14, 2026-10-01: X6, X7,
O19, X3): a surface with diagnostics returns its value(s) and LAST an
`info` dict and leaves nothing on the PAC; one misuse has one reaction; a
refusal names the method the user called; `am_pm_noise` stops on
`pnoise`'s rule.  Each test here failed on the code before."""
import numpy as np
import pytest

from pycircuit.circuit import IS, C, L, R, SubCircuit, VSin, circuit, gnd
from pycircuit.circuit.elements import BSource
from pycircuit.circuit.shooting import PAC, PSS
from pycircuit.circuit.simwarnings import (
    AccuracyWarning,
    ModelWarning,
    UsageWarning,
)
from pycircuit.circuit.tests._warnpolicy import quiet

#: what the PAC kept between calls until 2026-10-01
_OLD_ATTRIBUTES = ('alias_stop', 'sidebands_used', 'sampled_instants',
                   'lineshape_info', 'deflated', 'matvecs', 'event_shifts',
                   'time_response', '_last_second_multiplier')


def _driven_rc(noisy=True, npts=60):
    """A driven RC low-pass: a PSS, its PAC, the output index, the period."""
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['vs'] = VSin('in', gnd, vac=1.0, va=1.0, freq=1e3)
    cir['R'] = R('in', 'out', r=1e3, noisy=noisy)
    cir['C'] = C('out', gnd, c=1e-7)
    T = 1e-3
    pss = PSS(cir, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / npts, maxiterations=40)
    assert pss.converged
    return cir, pss, PAC(cir, toolkit=circuit.numeric), 'out', T


def _vdp(psd=1e-6, npts=120):
    """A van der Pol oscillator with a white current source at its core."""
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.2 * u * u)
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
    T = 2.0 * np.pi
    pss = PSS(cir, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / npts, x0=np.array([2.0, 0.0]),
                  maxiterations=300)
    assert pss.converged
    return cir, pss, PAC(cir, toolkit=circuit.numeric)


def test_the_surfaces_return_their_diagnostics_and_leave_none_on_the_pac():
    """`(value, info)` everywhere a surface has diagnostics -- `pnoise`,
    `am_pm_noise`, `oscillator_spectrum`, the sampled family, the adjoint
    rows -- and `PAC.solve`'s `result.info`; until 2026-10-01 they were PAC
    attributes the next call overwrote (`am_pm_noise` left an earlier
    `pnoise`'s `alias_stop` standing)."""
    _cir, pss, pac, out, T = _driven_rc()
    f0 = 1.0 / T
    with quiet(AccuracyWarning):
        S, info = pac.pnoise(pss, 0.3 * f0, out)
        assert S > 0.0 and info['stop'] in ('ratio', 'bound')
        assert info['sidebands'][:3] == [0, 1, -1] and info['deflated'] is False
        am, pm, info = pac.am_pm_noise(pss, 0.3 * f0, out, sweeptype='relative')
        assert am > 0.0 and pm > 0.0 and info['bands'][0] == 0, info
        Sn, info = pac.sampled_noise(pss, out, [0.25 * T], [0.1 * f0, 0.2 * f0])
        assert Sn.shape == (1, 2) and info['instants'].shape == (1,)
        v, info = pac.sampled_variance(pss, out, [0.25 * T], 0.01 * f0, 0.5 * f0)
        assert v.shape == (1,) and 'instants' in info
        row, info = pac.adjoint_transfer_row(pss, 0.3 * f0, out)
        assert row.shape == (pss.cir.n - 1,) and info['deflated'] is False
        rows, info = pac.adjoint_sideband_row(pss, 0.3 * f0, out, [0, 1])
        assert rows.shape == (2, pss.cir.n - 1) and 'matvecs' in info
        res = pac.solve(pss, [0.3 * f0])
        assert set(res.info) == {'time_response', 'event_shifts', 'deflated',
                                 'matvecs'}, res.info
    _cir, opss, opac = _vdp()
    with quiet(ModelWarning, AccuracyWarning):
        Sv, _L, info = opac.oscillator_spectrum(opss, [1e-3], 0)
    assert Sv.shape == (1,) and info['frequency_aware'] == 'first order', info
    for p in (pac, opac):
        assert not [a for a in _OLD_ATTRIBUTES if hasattr(p, a)]


def test_a_noiseless_circuit_has_zero_noise_said_silently():
    """A density of a noiseless circuit is 0, without a word; a RATIO with no
    value raises.  Until 2026-10-01: `pnoise` ran to the Nyquist and warned
    "lower bound" about an exact zero, `orbital_spectrum` refused "no orbital
    line", `oscillator_edge_jitter` raised "is any source noisy?" in one
    place and returned 0 in another, `band_spread` returned `inf`.  (Under
    the suite's policy an unexpected warning is an error, so the silence is
    asserted by the calls themselves.)"""
    _cir, pss, pac, out, T = _driven_rc(noisy=False)
    f0 = 1.0 / T
    S, info = pac.pnoise(pss, 0.3 * f0, out)
    assert S == 0.0 and info['stop'] == 'zero', info
    am, pm, info = pac.am_pm_noise(pss, 0.3 * f0, out, sweeptype='relative')
    assert am == 0.0 and pm == 0.0 and info['stop'] == 'zero', info
    with pytest.raises(ValueError, match='undefined'):
        pac.band_spread(pss, out, (0.1, 0.2), points=3)
    _cir, opss, opac = _vdp(psd=0.0)
    with quiet(ModelWarning):
        S = opac.orbital_spectrum(opss, [1e-3, 1e-2], 0, maxharmonics=8)
    assert np.all(S == 0.0), S
    with quiet(ModelWarning):
        ej = opac.oscillator_edge_jitter(opss, 0, 0.0, kmax=3)
    assert ej['sigma_t'] == 0.0 and np.all(ej['k_cycle'] == 0.0), ej


def test_one_misuse_has_one_reaction():
    """A harmonic or a sideband is an integer -- it was truncated by `int()`
    on some surfaces and used as a float on others (`carrier_phasor`: the
    coefficient at 1.5 f0) -- and >= 1 where a carrier is read; a count a
    non-negative integer that raises above the grid everywhere (a negative
    `maxsidebands` summed nothing and returned 0; the modal one was taken as
    given); a switch a bool (`frequency_aware=None` meant True)."""
    _cir, pss, pac, out, T = _driven_rc()
    f0 = 1.0 / T
    with quiet(AccuracyWarning):
        for call in (lambda: pac.carrier_phasor(pss, out, 1.5),
                     lambda: pac.am_pm(pss, 0.3 * f0, out, harmonic=1.5),
                     lambda: pac.am_pm_noise(pss, 0.3 * f0, out, harmonic=1.5),
                     lambda: pac.adjoint_sideband_row(pss, 0.3 * f0, out, [0.5]),
                     lambda: pac.mixer_response(pss, 0.3 * f0, out, (0, 1.5)),
                     lambda: PAC.lorentzian([1.0], 1e-9, 1.0, harmonic=1.5)):
            with pytest.raises(TypeError, match='must be an integer'):
                call()
        assert pac.carrier_phasor(pss, out, 1.0) == pac.carrier_phasor(pss, out, 1)
        with pytest.raises(ValueError, match='harmonic 0 is below 1'):
            pac.am_pm_noise(pss, 0.3 * f0, out, harmonic=0)
        for call in (lambda: pac.pnoise(pss, 0.3 * f0, out, maxsidebands=-1),
                     lambda: pac.am_pm_noise(pss, 0.3 * f0, out, maxsidebands=-1),
                     lambda: pac.sampled_noise(pss, out, [0.0], [0.1 * f0],
                                               maxsidebands=-1)):
            with pytest.raises(ValueError, match='maxsidebands -1 is below 0'):
                call()
    _cir, opss, opac = _vdp()
    N = len(opss.factored_period().steps)
    with quiet(ModelWarning, AccuracyWarning):
        with pytest.raises(ValueError, match='above what the grid resolves'):
            opac.modal_spectrum(opss, [1e-2], 0, maxharmonics=4,
                                maxsidebands=2 * (N // 2 - 1) + 1)
        with pytest.raises(TypeError, match='True or False'):
            opac.oscillator_spectrum(opss, [1e-3], 0, frequency_aware=None)
        with pytest.raises(ValueError, match='kmax 0 is below 1'):
            opac.oscillator_edge_jitter(opss, 0, 0.0, kmax=0)
        with pytest.raises(TypeError, match='positional'):
            opac.oscillator_edge_jitter(opss, 0, 0.0, 3)
    ## `floquet_modes(nmodes, fp)`: the ignored first parameter is gone
    with quiet(UsageWarning):
        assert len(opss.floquet_modes(1)) == 1


def test_a_refusal_names_the_method_the_user_called():
    """The driven-circuit refusals of the oscillator surfaces, and the
    oscillator refusal of the sampled family, name the method called --
    they named the inner one (`oscillator_covariance`, `colour_projection`,
    `coloured_diffusion_resolved`, `modal_spectrum`, `sampled_noise`)."""
    _cir, pss, pac, out, _T = _driven_rc()
    for name, call in (
            ('oscillator_edge_jitter', lambda: pac.oscillator_edge_jitter(pss, out, 0.0)),
            ('orbital_mode_weights', lambda: pac.orbital_mode_weights(pss)),
            ('phase_psd', lambda: pac.phase_psd(pss, [10.0])),
            ('coloured_diffusion', lambda: pac.coloured_diffusion(pss, [10.0])),
            ('correlation_spectrum', lambda: pac.correlation_spectrum(pss, [10.0], out)),
            ('orbital_spectrum', lambda: pac.orbital_spectrum(pss, [10.0], out)),
            ('orbital_correlation', lambda: pac.orbital_correlation(pss))):
        with pytest.raises(ValueError, match=f'^PAC.{name}:'):
            call()
    _cir, opss, opac = _vdp()
    with pytest.raises(ValueError, match='^PAC.sampled_variance:'):
        opac.sampled_variance(opss, 0, [0.0], 1e-3, 0.05)


def test_am_pm_noise_stops_on_the_rule_pnoise_stops_on():
    """`am_pm_noise` summed every band pair to the cap whatever its size
    (the cost always the maximum, nothing said whether the sum converged);
    it stops on `pnoise`'s ratio test now, band pairs in `|p|` order, and
    reports the stop.  The early stop agrees with the full sum to the
    ratio's own size."""
    _cir, pss, pac = _vdp()
    f0 = 1.0 / float(pss.period)
    with quiet(ModelWarning):
        am, pm, info = pac.am_pm_noise(pss, 0.01 * f0, 0)
    with quiet(ModelWarning, AccuracyWarning):
        am_f, pm_f, info_f = pac.am_pm_noise(pss, 0.01 * f0, 0, ratio_tol=0.0)
    assert info['stop'] == 'ratio' and info_f['stop'] == 'bound', (info, info_f)
    assert len(info['bands']) < len(info_f['bands'])
    assert abs(pm / pm_f - 1.0) < 1e-8 and abs(am / am_f - 1.0) < 1e-6, \
        (pm / pm_f - 1.0, am / am_f - 1.0)
