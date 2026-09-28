"""Shooting tests: shooting pnoise.  Split out of test_analysis_shooting.py on
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
    _NuMult,
    _PllMultPd)
from pycircuit.circuit.tests._shooting_fixtures import _noise_seam
from pycircuit.circuit.tests._shooting_fixtures import (_Flicker,
    _a9_vdp,
    _adjoint_ladder,
    _coloured_vdp,
    _cos2_mixer,
    _diode_mixer,
    _sampler_fixture,
    _sw,
    _vdp_ppv_method)


def _divider():
    """A linear divider with noisy resistors, and an independent answer."""
    c = SubCircuit()
    n1, n2 = c.add_nodes('net1', 'net2')
    c['vs'] = VSin(n1, gnd, va=1.0, vac=1.0, freq=1e3)
    c['R1'] = R(n1, n2, r=9e3)
    c['R2'] = R(n2, gnd, r=1e3)
    c['C'] = C(n2, gnd, c=1e-9)
    return c


def _pnoise_at(cir, per, fout, node, npts, **kw):
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / npts, maxiterations=40)
    assert pss.converged
    irn = pss.irefnode
    k = cir.get_node_index(node)
    k = k - 1 if k > irn else k
    pac = PAC(cir, toolkit=circuit.numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        S, used = pac.pnoise(pss, fout, k, **kw)
    return S, used, pac


def test_pnoise_reduces_to_the_stationary_analysis_on_a_linear_circuit():
    """⚠ OKUMURA'S `p = 1` CASE, against a reference pnoise cannot influence.

    A linear circuit converts nothing, so every sideband but `l = 0`
    vanishes and the fold collapses to the ordinary stationary formula —
    "exactly the same as that derived for a stationary noise". That makes
    `analysis_ss.Noise` the answer: a different analysis, an `(sC + G)`
    solve, no monodromy and no period anywhere in it.

    Measured: the `l != 0` terms come back ~1e-32 of the total, and the
    ratio against `Svnout` is 1.000000.

    ⚠ THE RATE, NOT THE SIZE. Per doubling of the period grid the residual
    falls 7.80x / 7.33x / 6.75x, reaching 3.9e-09 at 400 points. That is at
    least second order — the exact exponent is not asserted here because it
    has not been established, only that the disagreement is discretisation
    and not a defect.
    """
    from pycircuit.circuit.analysis_ss import Noise
    circuit.default_toolkit = circuit.numeric
    per, fout = 1e-3, 700.0
    rels = []
    for npts in (50, 100, 200):
        cir = _divider()
        ref = complex(Noise(cir, inputsrc='vs',
                            outputnodes=(cir.get_node('net2'), gnd)
                            ).solve(fout)['Svnout']).real
        S, used, _pac = _pnoise_at(cir, per, fout, cir.get_node('net2'), npts)
        assert S > 0, 'pnoise returned %r' % S
        rels.append(abs(S - ref) / ref)
        assert 0 in used and len(used) > 1, \
            'the accumulation never looked past l=0, so the fold is untested'

    assert rels[-1] < 1e-6, \
        'pnoise disagrees with the AC noise analysis by %.3e on a LINEAR ' \
        'circuit, where the two compute the same quantity' % rels[-1]
    for a, b in zip(rels, rels[1:]):
        assert a / b > 4.0, \
            'the residual against the stationary analysis falls only %.2fx ' \
            'per doubling (%.3e -> %.3e). A constant offset would look ' \
            'small and never move' % (a / b, a, b)


def test_pnoise_folds_and_the_fold_is_not_a_rounding_term():
    """On a mixer the sidebands carry a large share of the output noise.

    62% here, so an implementation that quietly summed only `l = 0` would
    be wrong by a factor, not by a tolerance — and would still return a
    plausible-looking PSD. The linear test above cannot catch that, because
    there the fold is *supposed* to contribute nothing.
    """
    cir = _diode_mixer()
    S_all, used, _p = _pnoise_at(cir, 1e-6, 3e5, 2, 80)
    S_l0, _u, _p2 = _pnoise_at(_diode_mixer(), 1e-6, 3e5, 2, 80,
                               maxsidebands=0)
    assert S_l0 > 0 and S_all > S_l0
    share = (S_all - S_l0) / S_all
    assert share > 0.25, \
        'the sidebands contribute only %.1f%% of the output noise on a ' \
        'driven diode; either the fold is not working or this circuit ' \
        'stopped mixing' % (100 * share)
    assert max(abs(np.asarray(used))) > 5, \
        'the accumulation stopped after |l|=%d on a switching circuit' \
        % max(abs(np.asarray(used)))


def test_pnoise_says_when_the_grid_stopped_it_rather_than_the_series():
    """⚠ WHICH STOPPING RULE FIRED IS PART OF THE ANSWER.

    Ending on the ratio test means the series converged. Ending on the
    grid's Nyquist means the grid ran out first, and every sideband above
    it is MISSING rather than small — the number is a lower bound. A
    strongly switching circuit reaches that readily: measured on this
    diode, 80 and 160 points per period both end on the bound, 320 and 640
    end on the ratio test, and the totals differ by only 0.04% — so the
    warning is not proof of a bad answer, it is a statement that the
    accumulation cannot vouch for itself.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _diode_mixer()
    pss = PSS(cir, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-6, timestep=1e-6 / 80, maxiterations=40)
    irn = pss.irefnode
    k = cir.get_node_index(2)
    k = k - 1 if k > irn else k
    pac = PAC(cir, toolkit=circuit.numeric)
    with pytest.warns(RuntimeWarning, match='Nyquist'):
        pac.pnoise(pss, 3e5, k)
    assert pac.alias_stop == 'bound'

    ## the linear circuit's series does converge, and says so
    _S, _u, pac2 = _pnoise_at(_divider(), 1e-3, 700.0,
                              _divider().get_node('net2'), 100)
    assert pac2.alias_stop == 'ratio', \
        'a linear circuit folds nothing and must stop on the ratio test, ' \
        'not on the grid'


def test_pnoise_refuses_cyclostationary_sources():
    """⚠ A BIAS-DEPENDENT `CY` IS A DIFFERENT MODEL, NOT A HARDER SUM.

    Okumura's cyclostationary treatment windows each source to a single
    timestep, and the windows' Fourier coefficients CORRELATE the
    sidebands: they stop adding in power and pick up cross terms needing
    the `R_{m,n}` construction. Summing powers anyway would answer a
    different question and look entirely normal doing it.

    Every noise source in the discrete element library is bias-independent
    — a resistor's `4kT/R` does not read `x` at all — so this refusal is
    unreachable there. A compact device's `CY` does read `x`, which is why
    the check samples the orbit rather than reasoning from element types.
    """
    class _BiasNoisyR(R):
        """A resistor whose noise follows its own terminal voltage."""
        def CY(self, x, w, epar=circuit.defaultepar):
            base = 4 * self.toolkit.kboltzmann * epar.T / self.iparv.r
            iPSD = base * (1.0 + 10.0 * abs(float(np.asarray(x).ravel()[0])))
            return self.toolkit.array([[iPSD, -iPSD], [-iPSD, iPSD]])

    import warnings
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    n1, n2 = c.add_nodes('net1', 'net2')
    c['vs'] = VSin(n1, gnd, va=1.0, vac=1.0, freq=1e3)
    c['R1'] = _BiasNoisyR(n1, n2, r=9e3)
    c['R2'] = R(n2, gnd, r=1e3)
    c['C'] = C(n2, gnd, c=1e-9)
    pss = PSS(c, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=40)
    assert pss.converged
    with pytest.raises(NotImplementedError, match='cyclostationary'):
        PAC(c, toolkit=circuit.numeric).pnoise(pss, 700.0, 1)


def test_pac_refuses_an_operating_point_from_another_circuit():
    """⚠ THE DRIVEN-OSCILLATOR TRAP, WHICH IS NOT A TYPO GUARD.

    The natural way to model a driven oscillator is to solve the PSS of
    the bare oscillator and treat the injection as a perturbation. It is
    wrong: the injection DEVICE is present even when its SIGNAL is zero.
    Buonomo & Lo Schiavo — "in absence of the injection signal, the
    injection circuit affects the basic LC oscillator by CHANGING THE
    NONLINEARITY OF THE FEEDBACK LOOP … [it] can affect the start-up
    condition … OR ITS OSCILLATION AMPLITUDE, or both."

    So the free-running orbit of the circuit-with-the-device is not the
    orbit of the circuit-without-it, and every Floquet quantity built on
    the wrong one inherits the error. The analysis would converge and
    report a plausible number.

    ⚠ THE REFERENCE-NODE CHECK CANNOT CATCH THIS. Two circuits differing
    by one device have the same reference node and, as here, the same node
    count — so the existing guard passes and the answer is quietly about
    the wrong orbit.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    ## same topology, same node count, same refnode — different objects
    bare = _adjoint_ladder(3)
    other = _adjoint_ladder(3)
    pss = PSS(bare, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    assert pss.converged
    assert bare.n == other.n, 'the two circuits must look alike for this ' \
        'test to be about identity rather than shape'

    pac_other = PAC(other, toolkit=circuit.numeric)
    for call in (lambda: pac_other.solve(pss, [700.0]),
                 lambda: pac_other.adjoint_transfer_row(pss, 700.0, 1),
                 lambda: pac_other.adjoint_sideband_row(pss, 700.0, 1, 0),
                 lambda: pac_other.pnoise(pss, 700.0, 1)):
        with pytest.raises(ValueError, match='different circuit'):
            call()

    ## the matching pair is accepted
    pac_same = PAC(bare, toolkit=circuit.numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pac_same.adjoint_sideband_row(pss, 700.0, 1, 0)


class _ModulatedShot(IS):
    """Shot noise modulated by the local node voltage — `CY` reads `x`.

    The case `_cy_reduced` refuses and Hull & Meyer's construction is the
    standard treatment of. Their own worked example is shot noise
    modulated by the collector current.
    """

    def CY(self, x, w, epar=None):
        p = self.iparv.noisePSD * (
            1.0 + 0.8 * float(np.asarray(x).ravel()[0]))
        return self.toolkit.array([[p, -p], [-p, p]])


def _hm_circuit(source=None, psd=1e-18, per=1e-3):
    cir = SubCircuit()
    cir.add_node('a')
    cir.add_node('b')
    cir['vs'] = VSin('a', gnd, va=1.0, vac=1.0, freq=1.0 / per)
    cir['R'] = R('a', 'b', r=1e3)
    cir['C'] = C('b', gnd, c=1e-7)
    if source is not None:
        cir['n'] = source('b', gnd, i=0.0, noisePSD=psd)
    return cir


def _hm_pnoise(cir, freq, modulated, per=1e-3, npts=200):
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / npts, refnode=gnd,
                  maxiterations=40)
    assert pss.converged
    irn = pss.irefnode
    k = cir.get_node_index('b')
    k = k - 1 if k > irn else k
    d = np.zeros(cir.n - 1)
    d[k] = 1.0
    pac = PAC(cir, toolkit=circuit.numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        S, used = pac.pnoise(pss, freq, d, modulated=modulated)
    return S, used, np.asarray(pss.waveform[1], dtype=float)[k]


def test_the_modulated_path_reduces_to_the_stationary_one():
    """⚠ THE FIRST GATE, AND IT COSTS NOTHING: a constant `CY` must give
    the same answer through both paths, BIT FOR BIT.

    `pnoise(modulated=True)` replaces `_cy_reduced` — which samples three
    states and refuses if they differ — with a time-average over the whole
    orbit. When `CY` does not depend on `x` those are the same matrix, so
    any discrepancy is in the new quadrature rather than in the physics.

    That matters because the stationary path is gated to **1.000000**
    against `analysis_ss.Noise`; reducing to it exactly inherits that gate
    rather than asking for a new one. Measured ratio 1.000000000000 at
    1e3, 1e4 and 5e4 Hz, with the same sideband count.
    """
    cir = _hm_circuit()
    for f in (1e3, 1e4, 5e4):
        S0, u0, _x = _hm_pnoise(cir, f, False)
        S1, u1, _x = _hm_pnoise(cir, f, True)
        assert len(u0) == len(u1), \
            'f=%g: the sideband accumulation stopped differently (%d vs ' \
            '%d), so the two paths are not being compared at the same ' \
            'truncation' % (f, len(u0), len(u1))
        assert abs(S1 / S0 - 1.0) < 1e-12, \
            'f=%g: averaged %.9e against stationary %.9e; with a constant ' \
            'CY these must agree to round-off or the averaging quadrature ' \
            'is wrong' % (f, S1, S0)


def test_the_modulated_path_accepts_what_the_stationary_one_refuses():
    """⚠⚠ HULL & MEYER'S CONSTRUCTION, AND THE ROUTE TO MOS pnoise.

    `_cy_reduced` refuses a bias-dependent `CY` — correctly, because the
    stationary sum would be the wrong model. But **no physically correct
    MOS noise model has a state-independent `CY`**, so that refusal is in
    effect a blanket refusal of MOS pnoise. Hull & Meyer's answer is not
    to refuse: carry the modulation in the **response** rather than in the
    **sources** — one stationary source per device at the *cycle-averaged*
    bias, with `H_l` supplying the modulation.

    Okumura's route needs one independent source per timestep interval per
    device (~25,000 on a real circuit); this needs one. Same physics,
    `p` times cheaper, and `H_l` was already built.

    ⚠ ASSERTED AS AN IDENTITY, NOT A PLAUSIBILITY. The modulated answer
    must equal what a *frozen* source at the time-averaged `CY` gives —
    not merely lie between the frozen extremes. Measured: node b swings
    −0.847…+0.847 so the `CY` factor spans 0.323…1.677, the frozen answers
    are 7.478e-15 and 3.887e-14, and the modulated answer is **2.317e-14**
    — the frozen-at-mean value, which here is also the midpoint because
    the modulation is linear in a variable whose orbit average is zero.

    ⚠ ITS CONDITION IS CHECKABLE AND IS THE OPPOSITE OF HIGH-Q: *"none of
    the large-signal state variables may change significantly over the
    decay time of the impulse response"*, and *"high-Q filters should be
    avoided, since they cause the impulse response to ring."* So this
    degrades exactly where `λ₂ → 1` — the same boundary as everything else
    here, from a fourth direction — which makes the two constructions
    **complementary**: this one for fast-settling circuits, Okumura's
    expensive one for the high-Q case that needs it.
    """
    modc = _hm_circuit(_ModulatedShot)
    with pytest.raises(NotImplementedError, match='BIAS-DEPENDENT CY'):
        _hm_pnoise(modc, 1e4, False)

    S, _u, xs = _hm_pnoise(modc, 1e4, True)
    assert np.isfinite(S) and S > 0

    lo_fac = 1.0 + 0.8 * float(xs.min())
    hi_fac = 1.0 + 0.8 * float(xs.max())
    assert lo_fac > 0.0, 'the modulation drove CY negative; not a PSD'
    S_lo, _u1, _x1 = _hm_pnoise(_hm_circuit(IS, psd=1e-18 * lo_fac),
                                1e4, True)
    S_hi, _u2, _x2 = _hm_pnoise(_hm_circuit(IS, psd=1e-18 * hi_fac),
                                1e4, True)
    assert min(S_lo, S_hi) <= S <= max(S_lo, S_hi), \
        'S = %.4e is outside the frozen bracket [%.4e, %.4e]' \
        % (S, min(S_lo, S_hi), max(S_lo, S_hi))
    ## the stronger statement: it IS the frozen-at-mean answer
    fp = None
    import warnings
    circuit.default_toolkit = circuit.numeric
    p2 = PSS(modc, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p2.solve(period=1e-3, timestep=1e-3 / 200, refnode=gnd,
                 maxiterations=40)
    fp = p2.factored_period()
    hs = np.diff(np.asarray(fp.times, dtype=float))
    n = min(len(hs), xs.size)
    mean_fac = 1.0 + 0.8 * float((xs[:n] * hs[:n]).sum() / hs[:n].sum())
    S_mean, _u3, _x3 = _hm_pnoise(_hm_circuit(IS, psd=1e-18 * mean_fac),
                                  1e4, True)
    assert abs(S / S_mean - 1.0) < 1e-9, \
        'modulated %.6e against frozen-at-mean %.6e; the construction is ' \
        'defined as the cycle-averaged source, so these are the same ' \
        'quantity and a difference is a quadrature error' % (S, S_mean)


def test_the_sideband_sum_gives_three_eighths_and_not_one_quarter():
    """⚠⚠ THE EXTERNAL GATE FOR THE CYCLOSTATIONARY PATH — and the wrong
    answer is a specific number, not merely a wrong shape.

    Roychowdhury, Long & Feldmann (1998) follow stationary noise through a
    cascade of mixers and report that a naive stationary-only analysis
    returns **¼** of the input power where the truth is **⅜** — *"50% more
    than that predicted by the previous naive analysis"*.

    Realised here by Parseval rather than by a mixer chain, which is the
    same content in one element. A transconductance `∝ cos²(ω₀t)` has
    `H₀ = ½` and `H_±2 = ¼`, so

        Σ_l |H_l|²  =  ¼ + 1/16 + 1/16  =  3/8      (full sum)
        |H₀|²       =  ¼                            (l = 0 alone)

    and `3/8` is also `E[cos⁴]`, a two-line trigonometric identity **using
    none of the machinery it gates**. That is what makes this external in
    the strong sense: its answer cannot be influenced by the
    implementation under test. It is this problem's `kT/C`.

    MEASURED, at three output frequencies:

        maxsidebands   S/(gain²·PSD)   sidebands used
             0           0.250021      [0]
             2           0.375023      [−2 … 2]
             4           0.375023      converged
             8           0.375023      converged

    ⚠ THE RATIO IS THE ASSERTION, AND IT IS DISCRETISATION-FREE. Both
    terms carry the same `H₀` error, so `full/l₀ = 1.500005` while each
    absolute value is 2.3e-05 off `3/8` at 400 points per period. A
    tolerance on the absolute number would have been a grid test; the
    ratio is the physics.

    ⚠ AND IT PINS THE ACCUMULATION, WHICH NOTHING ELSE DID. `pnoise` was
    gated as a *ratio* against `analysis_ss.Noise` on a circuit whose
    transfer is time-INVARIANT — where every sideband but `l = 0` is zero,
    so the sum was never exercised. This is the first absolute,
    analytically known target for `H_l` on a genuinely time-varying
    transfer, and truncating at `l = 0` reproduces exactly the error
    RLF98 names.
    """
    _cir, pss, pac, d, gain2 = _cos2_mixer()
    psd = 1e-18
    for fout in (10.0, 50.0, 137.0):
        S0, u0 = pac.pnoise(pss, fout, d, maxsidebands=0)
        Sf, uf = pac.pnoise(pss, fout, d, maxsidebands=4)
        assert sorted(u0) == [0]
        assert abs(S0 / (gain2 * psd) - 0.25) < 1e-3, \
            'f=%g: l=0 alone gives %.6f, not the 1/4 a stationary-only ' \
            'analysis returns' % (fout, S0 / (gain2 * psd))
        assert abs(Sf / (gain2 * psd) - 0.375) < 1e-3, \
            'f=%g: the full sum gives %.6f, not E[cos^4] = 3/8' \
            % (fout, Sf / (gain2 * psd))
        ## the discriminating number, free of the grid
        assert abs(Sf / S0 - 1.5) < 1e-4, \
            'f=%g: full/l0 = %.6f rather than 1.5 (1.76 dB). Both share ' \
            'the same H_0 discretisation error, so this ratio is the ' \
            'physics and a tolerance on the absolute value would not be' \
            % (fout, Sf / S0)
    ## and it converges where Parseval says it must: H_l = 0 for |l| > 2
    S2, _ = pac.pnoise(pss, 50.0, d, maxsidebands=2)
    S8, _ = pac.pnoise(pss, 50.0, d, maxsidebands=8)
    assert abs(S8 / S2 - 1.0) < 1e-12, \
        'sidebands beyond l = +-2 contribute %.3e; cos^2 has exactly ' \
        'three nonzero Fourier coefficients' % abs(S8 / S2 - 1.0)


def _cs_amp(va, f0=1e6, fnt=1.0):
    """A common-source stage on a real compact MOSFET, sinusoidally driven.

    ⚠ `va` MUST BE NONZERO. At `va = 0` the source is constant, the
    circuit has no periodic excitation, and `PSS` infers an AUTONOMOUS
    problem — `pnoise` then routes into the deflated oscillator solve and
    fails inside `ppv()` three frames away. The DC limit is approached
    with a small drive, not with no drive.
    """
    from pycircuit.circuit import compact
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    for nn in ('g', 'd', 'vdd'):
        cir.add_node(nn)
    cir['vdd'] = VS('vdd', gnd, v=1.2)
    cir['vg'] = VSin('g', gnd, v=0.7, va=va, freq=f0)
    cir['rl'] = R('vdd', 'd', r=5e3)
    cir['M'] = compact.PspMosLongChannel('d', 'g', gnd, gnd, fnt=fnt)
    return cir


def _cs_solved(va, f0=1e6, npts=12):
    """One PSS on the compact model, reused for every pnoise call.

    ⚠ THE PSS DOMINATES THE COST -- 30 s at 40 points against 4 s for
    `pnoise` -- and it is the compact model's evaluation, not the grid:
    `reltol` 1e-6 to 1e-10 all take ~30 s.  So the test solves each
    operating point ONCE and calls `pnoise` twice on it.

    ⚠ AND 12 POINTS IS ENOUGH, WHICH IS NOT AN APPROXIMATION.  The
    answer is 1.320116127e-16 at 12, 20, 30 and 40 points -- identical to
    ten digits -- because the cycle average is a rectangular rule on a
    PERIODIC smooth function, which converges spectrally.  The modulation
    is still fully seen: the departure from DC is 1e-4 at `va = 2e-2` on
    the same grid.
    """
    import warnings
    cir = _cs_amp(va, f0=f0)
    pss = PSS(cir, method='gear', reltol=1e-8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1.0 / f0, timestep=1.0 / (f0 * npts), refnode=gnd,
                  maxiterations=60)
    assert pss.converged, 'va = %r did not converge' % va
    irn = pss.irefnode
    k = cir.get_node_index('d')
    k = k - 1 if k > irn else k
    d = np.zeros(cir.n - 1)
    d[k] = 1.0
    return cir, pss, PAC(cir, toolkit=circuit.numeric), d


def _cs_pn(pac, pss, d, modulated, fout=1e5):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return pac.pnoise(pss, fout, d, modulated=modulated,
                          maxsidebands=2)[0]


def test_pnoise_runs_on_a_real_mosfet_through_the_modulated_path():
    """⚠⚠ MOS pnoise, end to end, on a compact surface-potential model.

    A common-source stage on `PspMosLongChannel` with `fnt = 1`, driven at
    1 MHz. Two things this establishes that nothing else does:

    ⚠ THE BIAS-DEPENDENT REFUSAL IS REACHABLE FROM A REAL DEVICE.
    `modulated=False` raises with *"BIAS-DEPENDENT CY (varies by 0.998
    over the orbit)"* — nearly a factor of two across one cycle. That was
    recorded here as unreachable while the model was believed noiseless;
    it is not.

    ⚠⚠ AND THE MODULATED ANSWER REDUCES TO `analysis_ss.Noise`, which is a
    different analysis on a different code path. Driving the gate ever
    more weakly must recover the DC noise at the same operating point:

        va       pnoise(modulated)   ratio to analysis_ss.Noise
        1e-4     1.320248e-16        1.000000
        5e-3     1.320239e-16        0.999994
        2e-2     1.320116e-16        0.999900
        5e-2     1.319418e-16        0.999371

    ⚠ THE DEPARTURE IS QUADRATIC IN `va`, TO THREE DIGITS — 6.0e-6,
    1.00e-4, 6.29e-4 against ratios of 16 and 6.25 in `va²`. That is the
    structural signature rather than a tolerance: the first-order term
    vanishes because a sinusoid's cycle average is zero, so the leading
    correction to a cycle-averaged `CY` is second order. A linear
    departure would mean the averaging was wrong; no departure at all
    would mean the modulation was being ignored.
    """
    from pycircuit.circuit.analysis_ss import Noise

    vas = (1e-4, 5e-3, 2e-2)
    solved = [_cs_solved(va) for va in vas]

    ## the refusal, on the same solve the modulated run uses
    _c, pss, pac, d = solved[-1]
    with pytest.raises(NotImplementedError, match='BIAS-DEPENDENT CY'):
        _cs_pn(pac, pss, d, False)

    cir0 = _cs_amp(0.0)
    ref = float(np.real(Noise(cir0, inputsrc='vg',
                              outputnodes=(cir0.get_node('d'), gnd),
                              toolkit=circuit.numeric
                              ).solve(1e5, complexfreq=False)['Svnout']))
    assert ref > 0

    dev = [abs(_cs_pn(pc, ps, dd, True) / ref - 1.0)
           for _cc, ps, pc, dd in solved]
    assert dev[0] < 1e-6, \
        'at va = 1e-4 the modulated pnoise is %.3e from analysis_ss.Noise; ' \
        'the weak-drive limit must recover the DC answer' % dev[0]
    assert dev[-1] > 5e-5, \
        'the strongest drive departs by only %.3e, so the modulation is ' \
        'not being seen at all' % dev[-1]
    ## quadratic: the departure grows as va^2, not as va
    expect = (vas[2] / vas[1]) ** 2
    got = dev[2] / dev[1]
    assert abs(got / expect - 1.0) < 0.15, \
        'va %.0e -> %.0e grew the departure by %.2f against the %.2f a ' \
        'quadratic term predicts; a LINEAR departure would mean the cycle ' \
        'average is wrong' % (vas[1], vas[2], got, expect)


def test_pnoise_refuses_a_harmonic_only_when_the_sources_are_undefined_there():
    """⚠⚠ ON A HARMONIC A SIDEBAND FOLDS THE SOURCES TO DC, and some
    device models are not defined there — but "harmonics are bad" is the
    wrong rule and would refuse a valid measurement.

    Sideband `l` evaluates `CY` at `f − l·f₀`, so `f = k·f₀` evaluates it
    at **zero**. Measured on `PspMosLongChannel`:

        fnt=1, nfa=0        CY(f=0) = nan    ← a DISABLED flicker term, 0/0
        fnt=1, nfa=8e22     CY(f=0) = inf    ← the real 1/f singularity

    ⚠ THE FIRST IS THE NASTIER ONE. A caller who sets `nfa = 0` believing
    flicker is off still gets `nan` out of `pnoise`, with no exception
    anywhere — the term is disabled and its `0/f` is still `0/0`.

    ⚠⚠ AND THE FOLD TO DC IS HARMLESS WHEN THE SOURCES ARE DEFINED THERE.
    A driven divider with white sources returns 1.490351e-17 at exactly
    `f₀`, and at 2f₀ and 3f₀. So the guard checks the **sources at the
    frequency that will actually be used**, rather than refusing a
    harmonic on principle. That distinction is the whole content: a rule
    keyed on "is this a harmonic" would break a working analysis, and one
    keyed on "is the result nan" would fire after the fact without saying
    why.

    ⚠ THIS BECAME REACHABLE ONLY WHEN THE MOS NOISE MODEL WAS EXERCISED.
    The roadmap has carried the sweep-grid trap as a known hazard since
    A4d and could not test it, because every source in the discrete
    library is white and defined everywhere.
    """
    from pycircuit.circuit import compact
    import warnings
    circuit.default_toolkit = circuit.numeric

    ## the white case must keep working, at several harmonics
    cir = _divider()
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-5, refnode=gnd, maxiterations=40)
    irn = pss.irefnode
    k = cir.get_node_index('net2')
    k = k - 1 if k > irn else k
    d = np.zeros(cir.n - 1)
    d[k] = 1.0
    pac = PAC(cir, toolkit=circuit.numeric)
    for m in (1, 2, 3):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            S, _u = pac.pnoise(pss, m * 1e3, d)
        assert np.isfinite(S) and S > 0, \
            'the white divider stopped working at %d*f0 (S = %r); a fold ' \
            'to DC is harmless when the sources are defined there' % (m, S)

    ## and both MOS cases must raise, naming the mechanism
    for nfa in (0.0, 8e22):
        cirm = _cs_amp(2e-2, fnt=1.0)
        cirm['M'] = compact.PspMosLongChannel('d', 'g', gnd, gnd,
                                              fnt=1.0, nfa=nfa, ef=1.0)
        pm = PSS(cirm, method='gear', reltol=1e-8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pm.solve(period=1e-6, timestep=1e-6 / 12, refnode=gnd,
                     maxiterations=60)
        assert pm.converged
        irn = pm.irefnode
        kk = cirm.get_node_index('d')
        kk = kk - 1 if kk > irn else kk
        dd = np.zeros(cirm.n - 1)
        dd[kk] = 1.0
        pcm = PAC(cirm, toolkit=circuit.numeric)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            off, _u = pcm.pnoise(pm, 1e6 + 1e2, dd, modulated=True,
                                 maxsidebands=2)
        assert np.isfinite(off), \
            'nfa=%r: 100 Hz off the harmonic is already broken (%r), so ' \
            'the guard below would be masking a different defect' % (nfa, off)
        with pytest.raises(ValueError, match='sits on a harmonic'):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pcm.pnoise(pm, 1e6, dd, modulated=True, maxsidebands=2)


def test_the_mos_flicker_term_shows_a_one_over_f_corner_in_pnoise():
    """⚠ THE COLOURED PATH, EXERCISED FROM A REAL DEVICE.

    `pnoise` evaluates each sideband's source at `f − l·f₀`, so a
    frequency-dependent `CY` is folded at the right frequency per
    sideband. With `PspMosLongChannel`'s flicker term enabled the output
    shows a `1/f` region below the corner and the thermal plateau above:

        f (Hz)     S              local slope
        1e2        1.499166e-15   −0.7465
        1e3        2.687271e-16   −0.2659
        1e4        1.456832e-16   −0.0383
        1e5        1.333788e-16   (plateau)

    Asserted as a *shape* — steepening toward low offset and flat at high
    — rather than against a slope value, because the corner's position is
    a property of the card's `nfa` and not of the analysis. What the
    analysis has to get right is that the two regions exist and are
    ordered, which a white-only path cannot produce at all.
    """
    from pycircuit.circuit import compact
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _cs_amp(2e-2, fnt=1.0)
    cir['M'] = compact.PspMosLongChannel('d', 'g', gnd, gnd,
                                         fnt=1.0, nfa=8e22, ef=1.0)
    pss = PSS(cir, method='gear', reltol=1e-8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-6, timestep=1e-6 / 12, refnode=gnd,
                  maxiterations=60)
    assert pss.converged
    irn = pss.irefnode
    k = cir.get_node_index('d')
    k = k - 1 if k > irn else k
    d = np.zeros(cir.n - 1)
    d[k] = 1.0
    pac = PAC(cir, toolkit=circuit.numeric)
    fs = np.array([1e2, 1e3, 1e4, 1e5])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        Ss = np.array([pac.pnoise(pss, f, d, modulated=True,
                                  maxsidebands=2)[0] for f in fs])
    assert np.all(np.diff(Ss) < 0), \
        'the output noise is not falling with frequency (%s); a flicker ' \
        'term must dominate at low offset' % Ss
    slopes = np.diff(np.log10(Ss)) / np.diff(np.log10(fs))
    assert slopes[0] < -0.5, \
        'low-offset slope is %.4f; with flicker enabled the spectrum must ' \
        'steepen toward DC' % slopes[0]
    assert slopes[-1] > -0.2, \
        'high-offset slope is %.4f; above the corner the thermal term ' \
        'must dominate and the spectrum flatten' % slopes[-1]
    assert slopes[0] < slopes[-1], 'the corner is not ordered'


def test_the_coloured_source_agrees_with_its_own_realisation():
    """⚠ THE GATE FOR EVERY COLOURED PATH: the same physics built two ways.

    An `IS(noiseTau)` on the tank node and a white `IS` filtered through
    an RC into a linear `BSource` are ONE noise, and `pnoise` must not
    know which it was given.  The filter stays out of the PSS (the periods
    agree to ten digits -- A4d's 1.2e-14 reconfirmed) and the sidebands
    agree to 1.3e-4 at three offsets, with the source-side frequency a
    visible `1/(1 + (2 pi f tau)^2)` away from the output-side one.
    """
    res = {}
    for kind in ('coloured', 'filtered'):
        _c, pss, pac, ov = _coloured_vdp(kind)
        f0 = 1.0 / float(pss.period)
        res[kind] = (float(pss.period),
                     [float(np.real(pac.pnoise(pss, f0 * (1.0 + k), ov)[0]))
                      for k in (1e-2, 1e-3, 1e-4)])
    assert abs(res['coloured'][0] - res['filtered'][0]) < 1e-9 * res['coloured'][0], \
        'the filter changed the PSS: it must stay out of it'
    for a, b in zip(res['coloured'][1], res['filtered'][1]):
        assert abs(a / b - 1.0) < 1e-3, \
            'coloured %.6e vs filtered %.6e: the element and its ' \
            'realisation disagree' % (a, b)


def test_the_coloured_lineshape_transform_meets_a_closed_form():
    """`_lineshape` (2026-09-26): the normalised lineshape with a coloured
    phase, `S(f) = 2 int exp(-D/2) cos(2 pi f tau) dtau`, against a CLOSED
    FORM.  A Lorentzian-shaped colour `c_c = c1 / (1 + (nu/nu_c)^2)` beside a
    white `c_w` gives `exp(-D/2) = e^{-(a + a1) tau} e^{b (1 - e^{-g tau})}`,
    a series of Lorentzians (`b = pi c1 / nu_c`, `g = 2 pi nu_c`).  With the
    colour's corner INSIDE the white linewidth (`b = 30`, the core reshaped):
    9.8e-7 / 9.6e-7 / -1.8e-7 / -6.0e-8 / -4e-11 at 0 .. 100 linewidths (the
    corner above the linewidth, `b = 0.3`: <= 3e-7).  A 1/f colour on a
    band, against an independent reference (the structure function exact
    through Ci, the transform by mpmath's `quadosc`): <= 9e-8 at 0 ..
    100 linewidths, at flicker/white 1 and 100 (not run here: 12 s per
    offset).  The `tau` grid was the binding knob (10 per decade left
    1.8e-4 at 100 linewidths)."""
    from pycircuit.circuit.shooting import _lineshape
    cw = 1e-3 / np.pi
    c1, nuc = 3.0 * cw, 1e-4
    a, a1 = 2 * np.pi ** 2 * cw, 2 * np.pi ** 2 * c1
    b, g = np.pi * c1 / nuc, 2 * np.pi * nuc

    def ref(f):
        import mpmath as mp
        mp.mp.dps = 60
        tot, bb, w = mp.mpf(0), mp.mpf(b), 2 * mp.pi * f
        for k in range(400):
            A = mp.mpf(a) + mp.mpf(a1) + k * mp.mpf(g)
            t = (-bb) ** k / mp.factorial(k) * 2 * A / (A * A + w * w)
            tot += t
            if k > 5 and abs(t) < mp.mpf(10) ** -40 * abs(tot):
                break
        return float(mp.e ** bb * tot)

    pc, converged = _lineshape.refine(
        lambda v: c1 / (1.0 + (np.asarray(v) / nuc) ** 2), 1e-9, 1e4)
    assert converged
    shape = _lineshape.ColouredLineshape(a, pc, 4.0)
    for f in (0.0, 1e-4, 1e-3, 1e-2, 1e-1):
        assert abs(shape(f) / ref(f) - 1.0) < 3e-6, (f, shape(f), ref(f))

#: `benchmarks/lineshape_reference.py`'s values (mpmath quadosc on D exact
#: through Ci; 30 .. 130 s an offset, so held here): (offset / f0, S)
_LINESHAPE_REFERENCE = {
    ## a real oscillator's line: `_lc_osc`'s levels, linewidth 5.5e-8 f0
    'narrow': (1.7496e-8, 7.5414e-11, [
        (1e-4, 188.42067129930248), (3e-4, 3.223113074359382),
        (1e-3, 0.09357472819258551), (3e-3, 0.004740821354755711),
        (1e-2, 0.0002503903112364973), (3e-2, 2.2233264903779104e-05),
        (0.1, 1.8250151581836042e-06), (0.3, 1.9719312579969512e-07)]),
    ## a very broad line (3e-4 f0), the flicker 1x / 100x the white at 1e-3
    'broad r1': (1e-4, 1e-7, [
        (3e-4, 232.02317932733632), (1e-3, 190.9054352560955),
        (3e-3, 42.09442940352864), (1e-2, 1.1883744313501698),
        (3e-2, 0.11578136501002126), (0.1, 0.010108157518428518),
        (0.3, 0.0011149225664776597)]),
    'broad r100': (1e-4, 1e-5, [
        (3e-4, 24.81165208710349), (1e-3, 24.765433371038284),
        (3e-3, 24.362930219330035), (1e-2, 20.23369816835472),
        (3e-2, 4.230113079727684), (0.1, 0.022906421522969535),
        (0.3, 0.001499721111301821)]),
}


def _handover_against_reference(regime):
    """The module's lineshape (`_lineshape.handover`, what
    `oscillator_spectrum` returns) on a reference regime: `(offsets, true
    relative errors, estimates)`."""
    from pycircuit.circuit.shooting import _lineshape
    c_w, k, rows = _LINESHAPE_REFERENCE[regime]
    fmin, fmax = 1e-7, 0.5
    offs = np.array([r[0] for r in rows])
    ref = np.array([r[1] for r in rows])
    pc, _conv = _lineshape.refine(lambda v: k / np.asarray(v, dtype=float), fmin, fmax)
    shapes = [_lineshape.ColouredLineshape(2 * np.pi ** 2 * c_w, pc, 4.0,
                                           per_decade=2 * _lineshape.TAU_PER_DECADE,
                                           shift=sh) for sh in (True, False)]
    ctab = _lineshape.ClampedTable(
        lambda v: c_w + np.where((v >= fmin) & (v <= fmax),
                                 pc(np.clip(v, fmin, fmax)), 0.0),
        1e-3 * fmin, 1e3, per_decade=400)

    def skirt(o):
        cf = c_w + k / o
        return tuple(_lineshape.second_order_skirt(
            o, ctab, 1.0, 1e-3 * fmin, fmax, c_at_f=cf, split=sp)
            for sp in (_lineshape.SKIRT_SPLIT, _lineshape.SKIRT_SPLIT_ALT)) + (
            cf / (o * o),)
    vals, errs = _lineshape.handover(shapes, offs, skirt)
    return offs, np.array(vals) / ref - 1.0, np.array(errs)


def test_the_coloured_lineshape_meets_an_independent_reference_across_the_handover():
    """The transform-to-skirt handover (2026-09-27; Andreas: "Go with A and
    B"), against mpmath on `D` exact through Ci.

    On a real oscillator's line the output was 6.5e-5 off at 1e-2 f0, the
    "~1e-4 at the handover".  The transform there was 3.4e-7.
      * A.  Its ESTIMATE was the fault.  QUADPACK's bound sums absolute
        bounds on O(1) pieces whose tiny difference is the skirt: 100 ..
        1e4 times the true error, so the linear skirt was taken half a
        decade early.  The grid-phase move is honest (within ~2.4x), and
        QUADPACK flagged no failure anywhere; its bound now counts only
        where it does.
      * B.  The skirt is second order (`second_order_skirt`): the phase
        above a split and the core's spread, taken as corrections:
        7.8e-4 / 6.5e-5 / 6.9e-6 first order -> 3.0e-6 / 9.0e-8 / 8.1e-9
        at 3e-3 / 1e-2 / 3e-2 f0.  Its estimate is the split move plus the
        square of its own correction.  The move alone read ZERO where
        `exp(-sH2)` underflowed (the heavy 1/f line at 3e-4 f0), and the
        skirt was taken 100 % off.
    Worst true error now: 3.4e-7 (real line), 4.0e-7 (heavy flicker), and
    4.5e-6 on the broad white-dominated line.  There the expansion
    parameter is ~1e-2, B gains only ~10x, and the transform's floor near
    0.1 f0 stands."""
    ## (A alone reached 9.8e-7 on the real line, so the bar there is 5e-7)
    worst = {'narrow': 5e-7, 'broad r1': 1e-5, 'broad r100': 1e-6}
    for regime, bar in worst.items():
        offs, err, est = _handover_against_reference(regime)
        assert np.max(np.abs(err)) < bar, (regime, dict(zip(offs, err)))
        ## the estimates are honest where the error is not negligible
        honest = np.abs(err) <= np.maximum(10.0 * est, 1e-8)
        assert np.all(honest), (regime, dict(zip(offs, zip(err, est))))
    ## B pinned on its own: the second-order skirt where it is taken on the
    ## real line (8e-9 / 1e-9 / 3e-10; the linear skirt 6.9e-6 / 6.3e-7 /
    ## 7.4e-8)
    from pycircuit.circuit.shooting import _lineshape
    c_w, k, rows = _LINESHAPE_REFERENCE['narrow']
    pc, _conv = _lineshape.refine(lambda v: k / np.asarray(v, dtype=float), 1e-7, 0.5)
    ctab = _lineshape.ClampedTable(
        lambda v: c_w + np.where((v >= 1e-7) & (v <= 0.5),
                                 pc(np.clip(v, 1e-7, 0.5)), 0.0),
        1e-10, 1e3, per_decade=400)
    for o, ref in rows[5:]:
        s2 = _lineshape.second_order_skirt(o, ctab, 1.0, 1e-10, 0.5,
                                           c_at_f=c_w + k / o)
        assert abs(s2 / ref - 1.0) < 5e-8, (o, s2 / ref - 1.0)


class _DcHeldNoise(IS):
    """`CY` proportional to the voltage across the element, which is a DC
    node held at 1 V on the orbit -- constant along it, ZERO at the zero
    vector.  Exists to pin that the cyclostationarity probes lie ON the
    orbit."""

    def CY(self, x, w, epar=None):
        xv = np.asarray(x).ravel()
        p = self.iparv.noisePSD * float(xv[0] - xv[1])
        return self.toolkit.array([[p, -p], [-p, p]])


def test_the_cyclostationarity_probes_lie_on_the_orbit():
    """⚠ A FALSE POSITIVE FROM A STATE OFF THE ORBIT, found by an external
    reference-simulator cross-check (2026-09-05): `_cy_reduced` sampled `CY` at three
    states, and one of them was the ZERO VECTOR -- on the orbit only by
    accident.  A switch model reading `goff` at `v(ck) = 0` had a linear
    time-invariant RC refused as cyclostationary.  Here a noise source
    whose `CY` is proportional to a DC-held node voltage is constant along
    the orbit and zero at the origin: `pnoise` must run.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir.add_node('b')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
    cir['Vb'] = VS('b', gnd, v=1.0)
    cir['n'] = _DcHeldNoise('b', gnd, i=0.0, noisePSD=1e-6)
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(cir.n - 1)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=6.6634, timestep=6.6634 / 240, x0=x0,
                  maxiterations=60)
    assert pss.converged
    ## ⚠ PRECONDITION (the property that makes this a test of the fix,
    ## 2026-09-05): `_DcHeldNoise`'s CY must actually DIFFER between the
    ## origin and the orbit.  If it ever stopped being state-dependent, CY
    ## would be constant everywhere, probe placement could not matter, and
    ## the gate below would pass while testing nothing.
    probe = _DcHeldNoise('b', gnd, i=0.0, noisePSD=1e-6)
    cy_zero = abs(np.asarray(probe.CY(np.zeros(2), 0.0))[0, 0])
    cy_orbit = abs(np.asarray(probe.CY(np.array([1.0, 0.0]), 0.0))[0, 0])
    assert cy_orbit > 1e3 * (cy_zero + 1e-300), \
        'the fixture is not state-dependent (CY %.3e at the origin vs ' \
        '%.3e on the orbit); probe placement cannot matter and the gate ' \
        'is vacuous' % (cy_zero, cy_orbit)
    pac = PAC(cir, toolkit=circuit.numeric)
    ov = [str(n) for n in cir.nodes].index('v')
    S = pac.pnoise(pss, 1.01 / pss.period, ov)[0]   # must not refuse
    assert np.all(np.isfinite(np.real(np.asarray(S))))


def test_pnoise_over_trbdf2_matches_the_stationary_analysis_and_folds():
    """pnoise is NATIVE over TR-BDF2 -- the two-stage sideband fold.

    A TR-BDF2 step injects the source at THREE abscissae, which the ordinary
    one-injection-per-step fold cannot represent (it gave 1e11 error and 99
    spurious sidebands).  The two-vector fold (`_sideband_forced_trbdf2`)
    carries the source coupling through both stages and is verified against
    forward driven solves to machine precision.  Two end-to-end checks:

    (1) LINEAR divider: no conversion, so the fold must collapse to l=0 and
        pnoise must reduce to the AC noise analysis (a different analysis, no
        period) -- to a few ppb, and stop on the ratio test, not the grid.
    (2) DIODE MIXER (converting): sidebands carry most of the noise; TR-BDF2
        must fold them and land on the Gear-2 answer (both compute the same
        physical noise), within the discretisation gap.
    """
    import warnings
    from pycircuit.circuit.analysis_ss import Noise
    circuit.default_toolkit = circuit.numeric

    ## (1) linear divider vs the AC-noise reference
    per, fout = 1e-3, 700.0
    cir = _divider()
    ref = complex(Noise(cir, inputsrc='vs',
                        outputnodes=(cir.get_node('net2'), gnd)
                        ).solve(fout)['Svnout']).real
    cir = _divider()
    pss = PSS(cir, method='trbdf2', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 200, maxiterations=40)
    k = cir.get_node_index(cir.get_node('net2'))
    k = k - 1 if k > pss.irefnode else k
    pac = PAC(cir, toolkit=circuit.numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        S, used = pac.pnoise(pss, fout, k)
    assert abs(S - ref) / ref < 1e-6, \
        'trbdf2 pnoise disagrees with AC noise by %.2e on a LINEAR circuit' \
        % (abs(S - ref) / ref)
    assert pac.alias_stop == 'ratio', \
        'a linear circuit folds nothing; trbdf2 must stop on the ratio test'

    ## (2) diode mixer: trbdf2 folds and lands on gear
    def mix(method):
        c = _diode_mixer()
        p = PSS(c, method=method, reltol=1e-11)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p.solve(period=1e-6, timestep=1e-6 / 160, maxiterations=40)
        kk = c.get_node_index(2)
        kk = kk - 1 if kk > p.irefnode else kk
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            s, u = PAC(c, toolkit=circuit.numeric).pnoise(p, 3e5, kk)
        return s, max(abs(np.asarray(u)))
    Sg, _lg = mix('gear')
    St, lt = mix('trbdf2')
    assert lt > 5, 'trbdf2 pnoise did not fold sidebands on the mixer (max l=%d)' % lt
    assert abs(St / Sg - 1.0) < 5e-3, \
        'trbdf2 pnoise %.4e vs gear %.4e on the mixer -- they should agree' \
        % (St, Sg)


def test_pnoise_over_radau_matches_the_stationary_analysis_and_folds():
    """pnoise is NATIVE over Radau IIA(3) -- the coupled three-vector sideband
    fold.

    A Radau step injects the source at THREE abscissae through the full
    ``A (x) B`` coupling (no two-vector shortcut), which
    `_sideband_forced_radau` carries through all three stages.  Same two
    end-to-end checks as the TR-BDF2 pnoise test:

    (1) LINEAR divider: no conversion, so the fold collapses to l=0 and pnoise
        reduces to the AC noise analysis to a few ppb, stopping on the ratio
        test, not the grid.
    (2) DIODE MIXER (converting): sidebands carry most of the noise; Radau
        folds them and lands on the Gear-2 answer within the discretisation
        gap.
    """
    import warnings
    from pycircuit.circuit.analysis_ss import Noise
    circuit.default_toolkit = circuit.numeric

    ## (1) linear divider vs the AC-noise reference
    per, fout = 1e-3, 700.0
    cir = _divider()
    ref = complex(Noise(cir, inputsrc='vs',
                        outputnodes=(cir.get_node('net2'), gnd)
                        ).solve(fout)['Svnout']).real
    cir = _divider()
    pss = PSS(cir, method='radau', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 200, maxiterations=40)
    k = cir.get_node_index(cir.get_node('net2'))
    k = k - 1 if k > pss.irefnode else k
    pac = PAC(cir, toolkit=circuit.numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        S, used = pac.pnoise(pss, fout, k)
    assert abs(S - ref) / ref < 1e-6, \
        'radau pnoise disagrees with AC noise by %.2e on a LINEAR circuit' \
        % (abs(S - ref) / ref)
    assert pac.alias_stop == 'ratio', \
        'a linear circuit folds nothing; radau must stop on the ratio test'

    ## (2) diode mixer: radau folds and lands on gear
    def mix(method):
        c = _diode_mixer()
        p = PSS(c, method=method, reltol=1e-11)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p.solve(period=1e-6, timestep=1e-6 / 160, maxiterations=40)
        kk = c.get_node_index(2)
        kk = kk - 1 if kk > p.irefnode else kk
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            s, u = PAC(c, toolkit=circuit.numeric).pnoise(p, 3e5, kk)
        return s, max(abs(np.asarray(u)))
    Sg, _lg = mix('gear')
    Sr, lr = mix('radau')
    assert lr > 5, 'radau pnoise did not fold sidebands on the mixer (max l=%d)' % lr
    assert abs(Sr / Sg - 1.0) < 5e-3, \
        'radau pnoise %.4e vs gear %.4e on the mixer -- they should agree' \
        % (Sr, Sg)


def test_am_pm_noise_does_not_depend_on_where_t_equals_zero():
    """A time shift of the drive cannot move the physical AM/PM split.

    Reported by a peer session and reproduced 2026-09-14: `am_pm_noise` formed
    `a + conj(b)` / `a - conj(b)` against the TIME ORIGIN, which is the
    carrier's frame only for a cosine-phased carrier.  On a diode driven
    through 1 k (thermal noise of R1 only, 100 Hz from a 10 kHz carrier):

        drive phase   carrier phase   am/pm before   am/pm after
          0 deg         -92.34 deg       0.3047         3.3159
         90 deg          -2.34 deg       3.2816         3.3159
         37 deg         -55.34 deg       0.6508         3.3159

    AM and PM SWAPPED under a quarter-period shift.  ⚠ The pnoise identity
    `S_am + S_pm = up + lo` could not see it -- the rotation leaves `|a|` and
    `|b|` alone -- which is why every existing gate passed.  On an oscillator
    the leak is `sin^2(phi)` of the 1/df^2 PM into AM, so AM rose toward the
    carrier instead of sitting flat below the corner.
    """
    import warnings as _w
    from pycircuit.circuit.semiconductors import ZenerDiode
    F0, fm = 1e4, 100.0
    got = []
    for ph in (0.0, 90.0):
        circuit.default_toolkit = circuit.numeric
        cir = SubCircuit()
        cir['V1'] = VSin('in', gnd, vo=1.0, va=0.8, freq=F0, phase=ph)
        cir['R1'] = R('in', 'n1', r=1e3)
        cir['C1'] = C('n1', gnd, c=1e-9)
        cir['D1'] = ZenerDiode('n1', gnd, IS=1e-13)
        pss = PSS(cir, method='gear', reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=1 / F0, timestep=1 / F0 / 800)
        pac = PAC(cir, toolkit=circuit.numeric)
        full = cir.get_node_index('n1')
        k = full - 1 if full > pss.irefnode else full
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            am, pm, _ = pac.am_pm_noise(pss, fm, k, carrier=1, maxsidebands=30)
            up, _ = pac.pnoise(pss, F0 + fm, k, maxsidebands=30)
            lo, _ = pac.pnoise(pss, F0 - fm, k, maxsidebands=30)
        up, lo = float(np.real(up)), float(np.real(lo))
        assert abs((am + pm) / (up + lo) - 1.0) < 1e-6, (ph, am + pm, up + lo)
        got.append((am, pm))
    (am0, pm0), (am90, pm90) = got
    assert am0 > 2.0 * pm0, \
        'the physical split is AM-dominant here (3.32); got am/pm %.4f' % (am0 / pm0)
    for x, y, name in ((am0, am90, 'S_am'), (pm0, pm90, 'S_pm')):
        assert abs(x / y - 1.0) < 1e-9, \
            '%s moved by %.3e under a 90-degree shift of the drive: the split ' \
            'is being taken against t = 0, not the carrier' % (name, x / y - 1.0)


def test_am_pm_noise_splits_the_sideband_pair_and_obeys_its_identity():
    """`PAC.am_pm_noise` splits output noise into AM and PM parts.

    ⚠ THE GATE IS AN IDENTITY, NOT A TOLERANCE.  `pnoise` at the upper sideband
    folds exactly the bands `g = freq + p f0`, and at the lower exactly their
    negatives, so

        S_am + S_pm  ==  pnoise(carrier*f0 + freq) + pnoise(carrier*f0 - freq)

    because `|a+c|^2 + |a-c|^2 = 2|a|^2 + 2|c|^2` leaves no cross term.  A
    pairing error in the band bookkeeping -- which sideband index reaches which
    output at which sign of `g` -- breaks it.  Measured: the residual falls
    9.0e-3 -> 3.3e-11 as the sideband count goes 4 -> 64, i.e. it is TRUNCATION
    and converges away, which a wrong pairing would not do.

    ⚠⚠ THE IDENTITY ALONE IS NOT ENOUGH, AND THAT IS THE POINT OF CHECK 3.  An
    implementation that simply returned HALF the total in each of AM and PM
    would satisfy it exactly while computing nothing.  So the split is also
    required to be NON-DEGENERATE: the whole content of an AM/PM decomposition
    is that the two are unequal, which happens only because the periodic
    operating point CORRELATES the two sidebands.  Uncorrelated sidebands carry
    equal AM and PM -- the classical LTI result -- so `S_pm == S_am` is exactly
    the answer that would mean the correlation had been lost.

    ⚠ AND CHECK 4 PINS THE CONJUGATE.  The split is `a ± conj(b)`, not `a ± b`:
    the sidebands counter-rotate about the carrier.  Dropping the conjugate
    still returns two positive numbers, so only a test that computes the naive
    form and finds it DIFFERENT keeps that from rotting.
    """
    import warnings
    from pycircuit.circuit.elements import Diode
    circuit.default_toolkit = circuit.numeric

    def mixer():
        c = SubCircuit()
        c['vs'] = VSin(1, gnd, vac=1.0, va=2.0, freq=1e6, phase=20)
        c['R'] = R(1, 2, r=1e4)
        c['D'] = Diode(2, gnd)
        c['C'] = C(2, gnd, c=1e-12)
        return c

    cir = mixer()
    T = 1e-6
    f0 = 1.0 / T
    pss = PSS(cir, method='gear', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 200, maxiterations=40)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    off = 0.13 * f0

    def residual(L):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            S_am, S_pm, _ = pac.am_pm_noise(pss, off, 2, carrier=1,
                                            maxsidebands=L)
            up, _ = pac.pnoise(pss, f0 + off, 2, maxsidebands=L)
            lo, _ = pac.pnoise(pss, f0 - off, 2, maxsidebands=L)
        return abs((S_am + S_pm) - (up + lo)) / (up + lo), S_am, S_pm

    ## 1. the identity holds once the sideband sum has converged
    r64, S_am, S_pm = residual(64)
    assert r64 < 1e-9, \
        'S_am + S_pm does not equal the noise in the sideband pair it splits ' \
        '(relative residual %.3e) -- the band pairing is wrong' % r64

    ## 2. and the residual at low sideband counts is TRUNCATION: it must fall.
    ##    A pairing error leaves a residual that does not converge away.
    r8, _, _ = residual(8)
    assert r8 > r64 * 100.0, \
        'the low-order residual (%.3e) is not larger than the converged one ' \
        '(%.3e), so the agreement is not the convergence it should be' \
        % (r8, r64)

    ## 3. NON-DEGENERATE: returning half the total in each would pass (1) exactly
    assert S_am > 0.0 and S_pm > 0.0
    ratio = S_pm / S_am
    assert abs(ratio - 1.0) > 0.1, \
        'AM and PM came out equal (ratio %.4f), which is the uncorrelated-' \
        'sideband answer -- either the correlation was lost or the split is ' \
        'returning half the total twice' % ratio

    ## 4. the CONJUGATE is load-bearing: the naive `a +- b` must differ
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = pac.adjoint_sideband_row(pss, off, 2, 1)[0]
        b = pac.adjoint_sideband_row(pss, -off, 2, 1)[0]
        cy = pac._cy_reduced(pss, 2.0 * np.pi * off)
    good = float(np.real((a + np.conj(b)) @ cy @ np.conj(a + np.conj(b))))
    naive = float(np.real((a + b) @ cy @ np.conj(a + b)))
    assert abs(good - naive) > 1e-3 * abs(good), \
        'the conjugate in `a + conj(b)` made no difference here, so this ' \
        'circuit cannot pin it -- pick one whose sidebands actually rotate'


def test_a_coloured_source_is_refused_whatever_else_is_in_CY():
    """`_refuse_coloured` judged colour against ONE GLOBAL SCALE at ONE state.

    Reported by a peer session and reproduced 2026-09-15.  The guard compared
    `CY(w0)` with `CY(10 w0)` at `x_last` only, against `1e-9 * max|CY|` over
    the WHOLE matrix.  Since the PSP gate resistor carries its 4kT/rg
    (`d0e7e4e`), a low-rg device puts a white `1.27e-20` in `CY`, and the
    drain's flicker colour (`2.5e-30` on its own `1.3e-27` at that state)
    passes the global threshold `1.27e-29` -- so `covariance` folded 1/f at
    w0 and returned a held variance 4.5 % high (peer's numbers), silently.

    On a sample-and-hold (IHP lv nmos 10/1 um, rg = 1.30 ohm, flicker on,
    swign = 0; 100 kHz, clock phase 270 so the switch is OFF at t = 0):

        rg on,  clock 270    RAN          <- the defect
        rg = 0, clock 270    refused
        rg on,  clock  90    refused      (x_last is a conducting state)
        rg on,  no flicker   ran, correctly

    Two weaknesses, both fixed: the global scale (a large white entry hides a
    small coloured one), and the single state (a source white at the period
    boundary but coloured elsewhere).  Now per entry, at `x_last` and states
    spread over the orbit.
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
        deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                              section='mos_tt')
        w, l = 10e-6, 1e-6
        base = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1, m=1,
                              pre_layout=1), w=w, l=l, T=T27)
        F = 1e5

        def sampler(over, clock_phase):
            circuit.default_toolkit = circuit.numeric
            kw = dict(base, swign=0.0)
            kw.update(over)
            cir = SubCircuit()
            cir['Vin'] = VSin('in', gnd, vo=0.3, va=0.2, freq=F, phase=0.0)
            cir['Vck'] = VSin('ck', gnd, vo=0.75, va=0.75, freq=F,
                              phase=clock_phase)
            cir['M1'] = PspMosLongChannel(cm.Node('in'), cm.Node('ck'),
                                          cm.Node('out'), gnd, **kw)
            cir['Ch'] = C('out', gnd, c=100e-12)
            pss = PSS(cir, method='gear', reltol=1e-12)
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                pss.solve(period=1 / F, timestep=1 / F / 100)
            return PAC(cir, toolkit=circuit.numeric), pss

        assert base['rg'] > 1.0, 'the fixture needs the low-rg device'
        for label, over in (('rg on', {}), ('rg = 0', dict(rg=0.0))):
            pac, pss = sampler(over, 270.0)
            with pytest.raises(NotImplementedError, match='COLOURED'):
                pac._refuse_coloured(pss, 'covariance')
        ## presence control: white sources only must pass
        pac, pss = sampler(dict(nfa=0.0, nfb=0.0, nfc=0.0), 270.0)
        pac._refuse_coloured(pss, 'covariance')
    finally:
        defaultepar.T = was


def test_the_ppv_border_caches_follow_a_re_solve():
    """`PAC._deflated_solve` keeps the PPV border `(v, u)` and
    `PSS.frequency_aware_ppv` the DC PPV, once per solved orbit (2026-09-26:
    `ppv()` was 14 of `am_pm_noise`'s 27 s, recomputed for every deflated
    solve; bit-identical cached).  Both are keyed on the state map, which a
    re-solve rebuilds: the SAME objects re-solved on another grid must give
    what fresh ones give."""
    import warnings as _w
    cir, pss = _a9_vdp(a=0.3)
    pac = PAC(cir, toolkit=circuit.numeric)
    f0 = 1.0 / float(pss.period)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pac.pnoise(pss, f0 * 1.01, 0, maxsidebands=8)
        pss.frequency_aware_ppv(0.01 * f0)
        T = float(pss.period)
        pss.solve(period=T, timestep=T / 300, x0=np.array([2.0, 0.0]),
                  maxiterations=300)
        assert pss.converged
        got = (float(np.real(pac.pnoise(pss, f0 * 1.01, 0, maxsidebands=8)[0])),
               np.asarray(pss.frequency_aware_ppv(0.01 * f0)[1]['samples_eq']))
        cir2, pss2 = _a9_vdp(a=0.3)
        pss2.solve(period=T, timestep=T / 300, x0=np.array([2.0, 0.0]),
                   maxiterations=300)
        ref = (float(np.real(PAC(cir2, toolkit=circuit.numeric).pnoise(
                   pss2, f0 * 1.01, 0, maxsidebands=8)[0])),
               np.asarray(pss2.frequency_aware_ppv(0.01 * f0)[1]['samples_eq']))
    assert got[0] == ref[0], (got[0], ref[0])
    assert np.array_equal(got[1], ref[1])


def test_pnoise_cyclostationary_is_the_stationary_fold_of_the_same_physics_and_the_cycle_average_is_not():
    """2026-09-08 (Andreas: "continue with the cyclostationary construction
    on the corrected Okumura").  `pnoise(cyclostationary=True)` folds a
    bias-dependent `CY` as white noise modulated by `B(t) = sqrt(CY(x(t)))`:
    `A_p = sum_k a_{p-k} B_k`, `S = sum_p |A_p|^2`, the modulation's
    harmonics read off the PSS samples (no window count) and the sideband
    rows the stationary fold already has.  Three gates:
      1. REDUCTION: a constant `CY` gives the stationary answer to 1e-12
         (measured 4.4e-16) -- Okumura's p = 1 case;
      2. IDENTITY: the same physics written two ways must agree.  A: a
         stationary white current into R_n, then a current-mode multiplier
         k1 V_n V_lo (the stationary fold through a periodically varying
         gain, exact).  B: an HDL source at the multiplier's output with
         PSD (k1 R_n V_lo(t))^2, cyclostationary.  Both through a second
         multiplier into an RC.  Measured 9e-16 with the fold in its
         `a P a^H` form (P the DFT of CY itself).  ⚠ A first version
         built the fold from `B = sqrt(CY)` and read 2.8e-5, flat in the
         offset, the sideband count and the grid, and present only when
         the LO crosses zero (7e-16 at va = 0.2, 2.8e-5 at va = 1); the
         square-root-free form removed it, so the pin is at 1e-12.
      3. the cycle-averaged route (`modulated=True`, Hull & Meyer's
         stationary equivalent) reads 0.533 of the truth here: the power is
         right and the correlation between sidebands is gone, which is the
         whole content of the construction.
      4. COLOURED (flicker): the band-resolved fold reduces to the
         stationary fold for a stationary flicker source (4e-16) and
         matches the stationary fold of the same separable physics for a
         sign-definite modulation (1e-9 pinned; measured 1.000000); see the
         block below for what a sign-changing modulation does.
    """
    import warnings
    from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, white_noise
    from pycircuit.utilities.param import Parameter
    circuit.default_toolkit = circuit.numeric

    class Mult(Behavioural):
        params_as = 'p'
        instparams = [Parameter(name='k', desc='gain', unit='A/V^2', default=1.0)]

        @staticmethod
        def analog(p, outp, outn, a, an, b, bn):
            return Contribution(Branch(outp, outn).I, p.k * Branch(a, an).V * Branch(b, bn).V)

    class ModNoise(Behavioural):
        params_as = 'p'
        instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

        @staticmethod
        def analog(p, outp, outn, b, bn):
            return Contribution(Branch(outp, outn).I, white_noise((p.k * Branch(b, bn).V) ** 2))

    T = 1e-6; f0 = 1.0 / T; k1 = 0.5; k2 = 0.3; Rn = 2.0

    def build(kind):
        c = SubCircuit()
        for n in ('lo', 'mid', 'out'):
            c.add_node(n)
        c['vlo'] = VSin('lo', gnd, va=1.0, vo=0.3, freq=f0)
        if kind == 'A':
            c.add_node('n'); c['xi'] = IS('n', gnd, i=0.0, noisePSD=1.0); c['Rn'] = R('n', gnd, r=Rn)
            c['M1'] = Mult('mid', gnd, 'n', gnd, 'lo', gnd, k=k1)
        else:
            c['src'] = ModNoise('mid', gnd, 'lo', gnd, k=k1 * Rn)
        c['Rm'] = R('mid', gnd, r=1.0)
        c['M2'] = Mult('out', gnd, 'mid', gnd, 'lo', gnd, k=k2)
        c['Ro'] = R('out', gnd, r=1.0); c['Co'] = C('out', gnd, c=0.2e-6)
        return c

    def solve(c):
        pss = PSS(c, method='gear', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 200, maxiterations=40)
        assert pss.converged
        return pss, PAC(c, toolkit=circuit.numeric)

    cA = build('A'); pA, pacA = solve(cA); oA = [str(n) for n in cA.nodes].index('out')
    cB = build('B'); pB, pacB = solve(cB); oB = [str(n) for n in cB.nodes].index('out')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for f in (0.13 * f0, 1.37 * f0):
            s_st, _ = pacA.pnoise(pA, f, oA, maxsidebands=16)
            s_cy, _ = pacA.pnoise(pA, f, oA, maxsidebands=16, cyclostationary=True)
            assert abs(s_cy / s_st - 1.0) < 1e-12, ('reduction', f / f0, s_st, s_cy)
            sB, _ = pacB.pnoise(pB, f, oB, maxsidebands=16, cyclostationary=True)
            assert abs(sB / s_st - 1.0) < 1e-12, ('identity', f / f0, s_st, sB)
            sBm, _ = pacB.pnoise(pB, f, oB, maxsidebands=16, modulated=True)
            assert abs(sBm / s_st - 1.0) > 0.3, ('the cycle average should differ by O(1)', sBm / s_st)
        ## and the stationary fold REFUSES the bias-dependent source, as before
        with pytest.raises(NotImplementedError, match='BIAS-DEPENDENT'):
            pacB.pnoise(pB, 0.13 * f0, oB, maxsidebands=16)

    ## COLOURED (the MOS flicker case, the peer's 24x point): the band-
    ## resolved fold.  Reduction: the stationary flicker source folded
    ## cyclostationary equals the stationary fold (4e-16).  Identity: a
    ## stationary flicker source through the multiplier (A) against an HDL
    ## flicker_noise((k R V_lo)^2, 1) source (B) -- the same SEPARABLE
    ## physics -- with a SIGN-DEFINITE modulation (va = 0.2, no crossing):
    ## 1.000000; and with the gain k V_lo^2 at va = 1: 1.000000000.  ⚠ With a
    ## sign-CHANGING modulation and a coloured source the two are NOT the
    ## same physics (0.563 / 1.325, grid-independent to six digits at 200 /
    ## 400 / 800 points): m xi and |m| xi coincide for white noise and
    ## differ for a coloured one whose correlation spans the sign change --
    ## Okumura's eq. 23 in concrete form; a PSD cannot carry the sign, so
    ## the fold (and the HDL model) is the |m| one.  Two defects the
    ## coloured gates found on the way: the pair index mirrored
    ## (B_{k+l'-l} for B_{k+l-l'}; 0.49 / 0.17 on the smooth identity,
    ## invisible to the reduction) and a circular wrap pairing a harmonic
    ## with the wrong band.
    from pycircuit.circuit.hdl import flicker_noise

    class ModFlicker(Behavioural):
        params_as = 'p'
        instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

        @staticmethod
        def analog(p, outp, outn, b, bn):
            return Contribution(Branch(outp, outn).I, flicker_noise((p.k * Branch(b, bn).V) ** 2, 1))

    def build_f(kind, va):
        c = SubCircuit()
        for n in ('lo', 'mid', 'out'):
            c.add_node(n)
        c['vlo'] = VSin('lo', gnd, va=va, vo=0.3, freq=f0)
        if kind == 'A':
            c.add_node('n'); c['xi'] = _Flicker('n', gnd, i=0.0, noisePSD=1.0, fref=1.0); c['Rn'] = R('n', gnd, r=Rn)
            c['M1'] = Mult('mid', gnd, 'n', gnd, 'lo', gnd, k=k1)
        else:
            c['src'] = ModFlicker('mid', gnd, 'lo', gnd, k=k1 * Rn)
        c['Rm'] = R('mid', gnd, r=1.0)
        c['M2'] = Mult('out', gnd, 'mid', gnd, 'lo', gnd, k=k2)
        c['Ro'] = R('out', gnd, r=1.0); c['Co'] = C('out', gnd, c=0.2e-6)
        return c

    cA = build_f('A', 0.2); pA, pacA = solve(cA); oA = [str(n) for n in cA.nodes].index('out')
    cB = build_f('B', 0.2); pB, pacB = solve(cB); oB = [str(n) for n in cB.nodes].index('out')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for f in (0.13 * f0, 1.37 * f0):
            s_st, _ = pacA.pnoise(pA, f, oA, maxsidebands=16)
            s_cy, _ = pacA.pnoise(pA, f, oA, maxsidebands=16, cyclostationary=True)
            assert abs(s_cy / s_st - 1.0) < 1e-12, ('coloured reduction', f / f0, s_st, s_cy)
            sB, _ = pacB.pnoise(pB, f, oB, maxsidebands=16, cyclostationary=True)
            assert abs(sB / s_st - 1.0) < 1e-9, ('coloured identity, sign-definite modulation', f / f0, s_st, sB)
    ## and the zero-crossing flicker fixture WARNS (a PSD cannot carry the sign)
    cB1 = build_f('B', 1.0); pB1, pacB1 = solve(cB1); oB1 = [str(n) for n in cB1.nodes].index('out')
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        pacB1.pnoise(pB1, 0.13 * f0, oB1, maxsidebands=16, cyclostationary=True)
    assert any('touches zero' in str(x.message) for x in w), [str(x.message)[:80] for x in w]
    ## and the SIGN-DEFINITE squared gain (k V_lo^2, exact to nine digits)
    ## does NOT warn: its sqrt(PSD) touches zero smoothly (the order of the
    ## zero, peer) -- the warning's negative control

    class ModFlicker2(Behavioural):
        params_as = 'p'
        instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

        @staticmethod
        def analog(p, outp, outn, b, bn):
            return Contribution(Branch(outp, outn).I, flicker_noise((p.k * Branch(b, bn).V ** 2) ** 2, 1))

    cB2 = build_f('B', 1.0); cB2['src'] = ModFlicker2('mid', gnd, 'lo', gnd, k=k1 * Rn)
    pB2, pacB2 = solve(cB2); oB2 = [str(n) for n in cB2.nodes].index('out')
    with warnings.catch_warnings(record=True) as w2:
        warnings.simplefilter('always')
        pacB2.pnoise(pB2, 0.13 * f0, oB2, maxsidebands=16, cyclostationary=True)
    assert not any('touches zero' in str(x.message) for x in w2), [str(x.message)[:80] for x in w2]


def test_mos_pnoise_runs_through_the_cyclostationary_route_and_the_cycle_average_overstates_a_switched_stage():
    """ITEM 2 (Andreas, 2026-09-09): MOS pnoise through the construction it
    was built for.  An EKV NMOS driven by a 1 MHz LO from off to strong
    inversion (modulated channel noise) into node x, then a second EKV as a
    pass transistor switched by the same LO into an RC load: the modulated
    noise passes through a periodically varying transfer AFTER being
    modulated, which is where the sideband correlation lives.  Measured:
    the stationary fold REFUSES (bias-dependent CY); thermal-only the
    cyclostationary fold reads 0.376 of the cycle average at 0.1 / 1.1 /
    3.3 MHz alike (broadband transfers), at the cycle average's cost (the
    white P-form is free); with flicker at ten times thermal (kf = 1e-13,
    the EKV's kf |I|^af / f) 0.319 / 0.436 / 0.379 through the coloured
    band-resolved branch (see the cost note at the end).  The direction: the switch's own
    channel noise is largest exactly when its channel shunts it, so
    avg(|H|^2 PSD) < avg(|H|^2) avg(PSD).  ⚠ A first fixture put the noise
    at the drain of a single stage into an RC load and read cyc/avg =
    1.000 to four digits at every offset: through a time-INVARIANT
    transfer only P_0 survives and the construction cannot show -- the
    fixture, not the fold.  Pinned: refusal; 0.376 within 2 %; the flicker
    ratio at 0.1 MHz within 3 % of 0.3055 and away from the thermal one.

    ⚠⚠ 0.319 WAS THE DEFECT'S NUMBER (2026-09-15).  It came from ONE square
    root of the summed CY, and the EKV carries thermal AND flicker noise in
    one element under different modulations -- a joint root makes such
    independent sources non-additive (measured +9.3 % of a held variance
    with white + flicker in one switch element, against the same sources as
    two elements, which the split reproduces to 2.3e-11).  With one root
    per component the ratio reads 0.3055; the thermal-only 0.376 is white
    and unchanged.
    """
    import warnings
    from pycircuit.circuit import elements_hdl as eh
    circuit.default_toolkit = circuit.numeric
    EKV = dict(vto=0.5, gamma=0.7, phi=0.7, kp=1.5e-4, cox=6.9e-3, w=10e-6, l=1e-6)
    f0 = 1e6; T = 1.0 / f0

    def build(kf):
        c = SubCircuit()
        for n in ('g', 'g2', 'x', 'd', 'vdd'):
            c.add_node(n)
        c['vdd'] = VS('vdd', gnd, v=2.0)
        c['vg'] = VSin('g', gnd, vo=0.5, va=0.8, freq=f0)
        c['vg2'] = VSin('g2', gnd, vo=2.2, va=1.2, freq=f0)
        card = dict(EKV); card.update(kf=kf, af=1.0)
        c['m1'] = eh.EkvNmosHdl('x', 'g', gnd, gnd, **card)
        c['Rx'] = R('vdd', 'x', r=2e3); c['Cx'] = C('x', gnd, c=0.2e-12)
        c['m2'] = eh.EkvNmosHdl('d', 'g2', 'x', 'x', **card)
        c['RL'] = R('vdd', 'd', r=5e3); c['CL'] = C('d', gnd, c=1e-12)
        return c

    ratios = {}
    for kf in (0.0, 1e-13):
        c = build(kf)
        pss = PSS(c, method='gear', reltol=1e-9)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 200, maxiterations=60)
        assert pss.converged
        pac = PAC(c, toolkit=circuit.numeric)
        od = [str(n) for n in c.nodes].index('d')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            with pytest.raises(NotImplementedError, match='BIAS-DEPENDENT'):
                pac.pnoise(pss, 0.1e6, od, maxsidebands=16)
            sm, _ = pac.pnoise(pss, 0.1e6, od, maxsidebands=16, modulated=True)
            sc, _ = pac.pnoise(pss, 0.1e6, od, maxsidebands=16, cyclostationary=True)
        ratios[kf] = sc / sm
    assert abs(ratios[0.0] / 0.376 - 1.0) < 0.02, ratios
    assert abs(ratios[1e-13] / 0.3055 - 1.0) < 0.03, ratios
    assert abs(ratios[1e-13] - ratios[0.0]) > 0.03, ratios

    ## THE COLOURED BRANCH'S COST (2026-09-09): the profile said the cost
    ## was the circuit's CY (53 000 evaluations, 231 bands x 230 samples),
    ## not the algebra.  `_cy_colour_model` fits A + B (w1/w)^ef per entry
    ## from three frequencies, verifies at two more, and serves every band
    ## and the stop rule from the model: 12.7 s -> 1.7 s on this fixture
    ## (the cycle average's own call is 2.0 s).  Pinned: (i) the model
    ## route equals the per-band evaluation to 1e-9 -- it was 1.5e-11 --
    ## and (ii) a colour that is NOT thermal-plus-flicker (a Lorentzian
    ## term added to the fixture's CY) fails the verification, so the
    ## fold falls back to evaluating the circuit and still agrees.
    ## ⚠ Since 2026-09-15 the fit is per ELEMENT (`_cy_components_model`),
    ## so the per-band reference is forced by failing THAT fit (every
    ## element is then evaluated per band), not the whole-circuit one.
    ## ⚠⚠ AND THE REFERENCE IS NO LONGER THE SAME MODEL TO 1e-9: evaluated
    ## per band an element gets ONE root for its white and flicker parts
    ## together (there is no fit to split them), while the fitted route
    ## roots them separately.  On the EKV (thermal + flicker in one device)
    ## that in-element difference measured 4.2e-4 -- the cross-ELEMENT
    ## independence, which both routes keep, is what moved 0.319 -> 0.3055.
    ## ⚠⚠⚠ AND SINCE 2026-09-28 IT IS THE SAME MODEL AGAIN: the EKV states
    ## its flicker's signed amplitudes, and per band an element is read
    ## from them plus the root of its WHITE remainder (`_perband_mode`
    ## 'white') -- white and flicker rooted separately, as the fitted route
    ## roots them: 3.2e-13 (was 4.2e-4).  (Not bit-equal: the per-band
    ## route did run.)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert pac._noise_components(pss).colour_model(0.1e6, f0) is not None
        _noise_seam(pac, colour_fit=staticmethod(lambda Cs, ws: None))
        try:
            sfull, _ = pac.pnoise(pss, 0.1e6, od, maxsidebands=16, cyclostationary=True)
        finally:
            del pac._noise_components
    assert sc != sfull and abs(sc / sfull - 1.0) < 1e-9, (sc, sfull)

    cy_orig = c.CY
    def cy_lorentz(x, w, **kw):
        base = np.asarray(cy_orig(x, w, **kw), dtype=complex)
        n_ = base.shape[0]; bump = np.zeros_like(base)
        bump[od, od] = 1e-22 / (1.0 + (w / (2.0 * np.pi * 3.0e6)) ** 2)
        return base + bump
    c.CY = cy_lorentz
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert pac._noise_components(pss).colour_model(0.1e6, f0) is None
        sl, _ = pac.pnoise(pss, 0.1e6, od, maxsidebands=16, cyclostationary=True)
        consulted = []
        _noise_seam(pac, colour_model=lambda *a_, **k_: consulted.append(1))
        try:
            slfull, _ = pac.pnoise(pss, 0.1e6, od, maxsidebands=16, cyclostationary=True)
        finally:
            del pac._noise_components
    ## ⚠ vacuous if the patch is never looked up
    assert consulted, 'the patched colour_model was never consulted'
    assert abs(sl / slfull - 1.0) < 1e-12, (sl, slfull)
    assert abs(sl / sc - 1.0) > 1e-3, (sl, sc)


def test_the_am_corner_is_f0_over_2pi_q_lambda_not_4pi():
    """⚠ A FACTOR OF TWO IN A SHIPPED DOCSTRING, caught by the docs session
    against `am_pm_noise` itself (2026-09-09) and reproduced here before the
    text was changed.  `S_am/S_pm` on a free-running oscillator is an exact
    Lorentzian in `u = offset/f0`,

        S_am/S_pm = u^2 / (u_c^2 + u^2),   u_c = 1/(2 pi Q_lambda),

    with `Q_lambda = -1/ln|lambda_2|` from the monodromy and nothing fitted:
    0.5007 at `u_c`, 0.2007 at the previously documented `f0/(4 pi Q)`,
    0.8007 at `2 u_c` (Q = 16, 240 points; Q = 8 the same to 1e-3).  The
    old corner named a Q twice as high and overstated the in-band AM 2.5x
    at its own corner.  ⚠ Measure the corner NEAR it: far below, `S_am` is
    dominated by an O(h^2) grid term that reads like "the AM floor is an
    artefact" (the peer's near-miss; Richardson gives a finite limit).
    """
    from pycircuit.circuit.shooting import PAC
    cir, pss = _vdp_ppv_method('gear', 240, Q=16.0)[:2]
    fp = pss.factored_period()
    M = np.column_stack([fp.matvec(e) for e in np.eye(fp.width)])
    lam = np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
    q_lam = -1.0 / np.log(lam[1])
    f0 = 1.0 / float(pss.period)
    pac = PAC(cir, toolkit=circuit.numeric)
    ov = [str(n) for n in cir.nodes].index('v')
    u_c = 1.0 / (2.0 * np.pi * q_lam)
    for u, expect in ((0.5 * u_c, 0.2), (u_c, 0.5), (2.0 * u_c, 0.8)):
        am, pm, _ = pac.am_pm_noise(pss, u * f0, ov, carrier=1, maxsidebands=32)
        ratio = float(np.real(am) / np.real(pm))
        assert abs(ratio - expect) < 0.01, (u / u_c, ratio, expect)


def _a2_tone_fixture(src, asym=0.25, npts=80):
    """The A2 fixture (van der Pol tank with series loss, a `u^2` asymmetry and
    one slow RC node, tau/T = 100) with a unit-PSD noise current at `src`
    ('v' = the tank, 'w' = behind the slow node), converged at `npts`."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    T0 = 6.6634
    c = SubCircuit()
    for n in ('v', 'w', 'x'):
        c.add_node(n)
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', 'x', L=1.0)
    c['Rl'] = R('x', gnd, r=0.2)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: (u - u ** 3 / 3.0) + asym * (u ** 2 - 2.0))
    c['Rs'] = R('v', 'w', r=1e2)
    c['Cs'] = C('w', gnd, c=100.0 * T0 / 1e2)
    c['n'] = IS(src, gnd, i=0.0, noisePSD=1e-6)
    pss = PSS(c, method='gear', reltol=1e-11)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T0, timestep=T0 / npts, x0=x0, maxiterations=200)
    assert pss.converged
    return c, pss


def _forward_tone_pm(c, pss, src, tones, r=0.10, periods=250, window=100, amp=1e-4):
    """PM by quadrature of the fundamental's sidebands from FORWARD small-signal
    tone transients on the PSS's own fixed grid: a real current tone at
    `|m + r| f0` for each `m` in `tones`, driven MINUS an undriven run from the
    same state (the integrator's phase slip cancels exactly), the carrier from
    the undriven run and the sidebands `U`, `L` at `f0 (1 +- r)` from the
    difference over the last `window` periods; `PM = 0.5 |u - conj(l)|^2`
    with `u = U/C`, `l = L/C` -- `am_pm_noise`'s own combination.  Summed over
    the tones (a white source weights every band equally).  Shares the orbit
    and the integrator with `pnoise`; shares neither the adjoint nor the
    sideband assembly."""
    import warnings
    T = float(pss.period)
    npts = len(pss.factored_period().steps) + 1
    h = T / (npts - 1)
    f0 = 1.0 / T
    isrc = c.get_node_index(src)
    iv = c.get_node_index('v')
    nfull = c.n
    x0 = np.asarray(pss.waveform[1], dtype=float)[:, 0].copy()

    def run(ws):
        def inject(t):
            u = np.zeros(nfull)
            u[isrc] = amp * np.cos(ws * t)
            return u
        tr = pss._new_transient(pss._integrator_for('gear'))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tr.solve(tend=periods * T, x0=x0, timestep=h,
                           provided_function=inject if ws is not None else None,
                           fixed_timestep=True)
        t = np.asarray(res.sweep_values, dtype=float)
        X = np.asarray(res.x, dtype=float)
        X = X if X.shape[0] == nfull else X.T
        return t, X[iv]

    t0, v0 = run(None)
    N = window * (npts - 1)
    tw = t0[-N:]
    vref = v0[-N:]

    def bin_amp(x, freq):
        return np.sum(x * np.exp(-1j * 2 * np.pi * freq * tw)) * 2.0 / N

    Cc = bin_amp(vref, f0)
    total = 0.0
    for m in tones:
        t, v = run(2 * np.pi * abs(m + r) * f0)
        dv = v[-N:] - vref
        u = bin_amp(dv, (1 + r) * f0) / Cc
        l = bin_amp(dv, (1 - r) * f0) / Cc
        total += 0.5 * abs(u - np.conj(l)) ** 2
    return total


def test_pnoise_oscillator_pm_matches_a_forward_tone_transient_with_no_adjoint():
    """⚠ THE FIRST EXTERNAL GATE FOR `pnoise`'S AUTONOMOUS PM, and the
    instrument that RESOLVED a day of Monte Carlo (roadmap, 2026-09-09).

    `pnoise`'s `S_pm` is the adjoint sideband assembly folded over bands
    and split by quadrature.  A forward small-signal tone transient on the
    PSS's own grid computes the same LPTV response with no adjoint and no
    assembly: inject a real current tone at `|m + r| f0` at the source node,
    subtract the undriven run, read the fundamental's sidebands, combine as
    `0.5 |u - conj(l)|^2`, sum over the bands the white source occupies.
    The SOURCE-LOCATION double ratio -- a source behind a slow RC node over a
    source in the tank -- is normalisation-free and is exactly the quantity
    a 16-seed Monte Carlo campaign had found drifting 6 % against both
    linear constructions as the orbit's asymmetry vanished.  This route
    agreed with `pnoise` to 1.003 / 1.010 at a = 0.25 / 0 (240 points, nine
    tones) and reproduced the Monte Carlo's drift DETERMINISTICALLY from
    the estimators' definitions of phase: the drift was the estimators (a
    one-period fundamental demodulation leaks the other harmonics'
    sidebands through its boxcar; zero crossings convert every harmonic's),
    not the constructions.

    Pinned here at 80 points and 250 periods (cost: ~2 min alone; 120
    points / 400 periods took 7 min and read 0.985 too): the slow-node
    source with three tones (m = -1, 0, 1; the other bands are 1.5 % of
    its sum) and the tank source with seven (m = -3..3; its |m| = 2, 3
    bands are 15 %), against `am_pm_noise` at r = 0.10 with 32 sidebands.
    Measured 0.9855 / 0.9860 / 0.9864 at 80 / 100 / 120 points (the 1.4 %
    is the bands left out); asserted within 4 % of 1.  ⚠ The undriven subtraction
    is not optional: at 240 points the integrator's own period error slips
    the phase ~6e-3 rad per 100 periods, which leaks the carrier into the
    sideband bins at the level of the response.  ⚠ Tone m = +1 doubled in
    amplitude gave 4.000x the power: linear.
    """
    import warnings
    r = 0.10
    ratio = {}
    pm_tone = {}
    for src, tones in (('w', (-1, 0, 1)), ('v', (-3, -2, -1, 0, 1, 2, 3))):
        c, pss = _a2_tone_fixture(src)
        pac = PAC(c, toolkit=circuit.numeric)
        ov = [str(n) for n in c.nodes].index('v')
        f0 = 1.0 / float(pss.period)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _am, pm, _ = pac.am_pm_noise(pss, r * f0, ov, carrier=1, maxsidebands=32)
        ratio[src] = float(np.real(pm))
        pm_tone[src] = _forward_tone_pm(c, pss, src, tones, r=r)
    double = (pm_tone['w'] / pm_tone['v']) / (ratio['w'] / ratio['v'])
    assert abs(double - 1.0) < 0.04, (double, pm_tone, ratio)


def test_pnoise_refuses_a_coloured_source_folded_onto_dc_at_a_clock_harmonic():
    """⚠⚠ PEER REPORT, 2026-09-15.  With a 1/f source, `pnoise(cyclostationary=
    True)` at EXACTLY a clock harmonic returned 6.3e-2 V^2/Hz against 9.2e-15
    at 0.1 % either side, silently.  `1/T` rounds (99999.999999999985 Hz), so
    the folded band sits 1.5e-11 Hz from DC -- finite and enormous for 1/f --
    and the harmonic guard refused only a NON-finite CY.  Now a
    frequency-dependent CY on a harmonic refuses; white stays allowed.
    ⚠ Also fixed and pinned: the colour model was NaN at negative band
    frequencies, so the ratio stop never fired for a coloured source and
    every call ran to the Nyquist bound with a spurious warning."""
    import warnings

    def els(c):
        c['S0'] = _sw()
        c['N0'] = IS('out', gnd, i=0.0, noisePSD=1e-26, noiseFc=1e6)
    cir, pss, io, pac, T = _sampler_fixture(els, npts=200)
    f0 = 1.0 / T
    ## both forms of "on the harmonic": `k/T` lands the folded band exactly on
    ## 0 (a ZeroDivisionError from inside the guard, before the fix); the
    ## literal `k * 100e3` lands it 1.5e-11 Hz off (the finite absurd value)
    for f in (f0, 2 * f0, 100e3, 200e3):
        with pytest.raises(ValueError, match='harmonic'):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pac.pnoise(pss, f, io, maxsidebands=90, cyclostationary=True)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        v, _used = pac.pnoise(pss, 1.001 * f0, io, maxsidebands=90,
                              cyclostationary=True)
    msgs = [str(w.message) for w in caught]
    assert pac.alias_stop == 'ratio', (pac.alias_stop, msgs)
    assert not [m_ for m_ in msgs if 'invalid value' in m_ or 'Nyquist' in m_], msgs
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        v0, _u = pac.pnoise(pss, 1.001 * f0, io, maxsidebands=90,
                            cyclostationary=True, ratio_tol=0.0)
    assert np.isfinite(v) and abs(float(np.real(v)) / float(np.real(v0)) - 1.0) < 1e-4, (v, v0)

    ## white sources on the harmonic are still answered
    wcir, wpss, wio, wpac, _T = _sampler_fixture(lambda c: c.__setitem__('S0', _sw()),
                                                 npts=200)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        vw, _u = wpac.pnoise(wpss, f0, wio, maxsidebands=90, cyclostationary=True)
    assert np.isfinite(vw) and float(np.real(vw)) > 0.0


def _pll_phase_noise(K, kvco, cf, fm_list, sf=1.0, fref=1e6, npts=400):
    """`(f_c from the multiplier, [pnoise at the PHASE node])` for the locked
    loop, seeded on the STABLE branch."""
    import warnings as _w
    from pycircuit.circuit.elements_hdl import VcoHdl
    T = 1.0 / fref
    c = SubCircuit()
    for nd in ('ref', 'vco', 'ph', 'pd', 'ctl'):
        c.add_node(nd)
    c['Vref'] = VSin('ref', gnd, va=1.0, freq=fref, phase=0.0)
    c['X1'] = VcoHdl('ctl', gnd, 'vco', gnd, 'ph', f0=fref, kvco=kvco,
                     va=1.0, modulus=1.0, sf=sf)
    c['PD'] = _PllMultPd('vco', gnd, 'ref', gnd, 'pd', gnd, k=K)
    c['Rf'] = R('pd', 'ctl', r=1e3)
    c['Cf'] = C('ctl', gnd, c=cf)
    names = [str(n) for n in c.nodes if str(n) != 'gnd!']
    x0 = np.zeros(c.n - 1)
    x0[names.index('X1._state0')] = 0.25          # STABLE branch, see step 2
    x0[names.index('ph')] = 0.25
    pss = PSS(c, method='gear', reltol=1e-10)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=100)
    assert pss.converged
    lam = float(np.max(np.abs(np.linalg.eigvals(np.asarray(pss._monodromy)))))
    assert lam < 1.0, 'seeded onto the saddle; noise about an unstable orbit'
    fc = -np.log(lam) / (2.0 * np.pi * T)
    pac = PAC(c, toolkit=circuit.numeric)
    iph = names.index('ph')
    out = []
    for fm in fm_list:
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            S, _nsb = pac.pnoise(pss, fm, iph)
        out.append(float(np.real(np.asarray(S, dtype=complex).ravel()[0])))
    return fc, out


def test_a_locked_loop_shapes_its_vco_phase_noise_with_nothing_fitted():
    """A6 step 4: the closed-loop phase noise against a CLOSED FORM, with both
    of its constants predicted by other measurements and nothing fitted.

    `VcoHdl` injects white FREQUENCY noise of PSD `sf`, so the phase -- which
    the `ph` node carries in CYCLES -- is its integral:

        S_phi,free(f_m) = sf / (4 pi^2 f_m^2)

    ⚠ THE `4 pi^2` IS THE CYCLES CONVENTION AND IS NOT COSMETIC.  `phi = int f
    dt` gives `Phi = N_f / (j 2 pi f_m)`.  `VcoHdl`'s own docstring quotes
    `S_phi = sf/f_m^2`, which is the RADIAN form; applying it to this node
    reads a factor of 39.5 low, and that error presents as a CONSTANT ratio
    across every offset -- which is exactly how it was caught.  A scale error
    flat across decades is a UNITS problem; a shape error would be physics.

    A first-order loop high-passes that with corner `f_c`, and the RC's own
    pole makes the loop second order, worth `2 f_c/f_RC`:

        S_phi,closed(f_m) = sf (1 + 2 f_c/f_RC) / (4 pi^2 (f_m^2 + f_c^2))

    NOTHING HERE IS FITTED.  `f_c` is read from the FLOQUET MULTIPLIER, which
    never looks at noise (`|lam| = exp(-2 pi f_c T)`, measured 100.063 Hz
    against `kvco*K/2` = 100.000), and `2 f_c/f_RC` from R and C.  Measured
    2026-09-18 at npts 400/800/1600 -- grid-independent to eight digits --
    agreeing to ~1e-4 across three decades of offset.

    The second-pole term was found, not assumed: the first form left a
    CONSTANT 1.2586e-3 which did not move with the grid (so not
    discretisation) and was identical below and far above the corner (so not
    the loop's first pole).  It scales linearly on two independent knobs:
    `Cf` 1e-10/1e-9/1e-8 gives 1.2569e-4 / 1.2586e-3 / 1.2716e-2, and cutting
    the PD gain 100x gives 1.257e-5 -- i.e. `2 f_c/f_RC` on both.
    """
    fref, T = 1e6, 1e-6
    kvco, K, sf, rf = 2e4, 1e-2, 1.0, 1e3
    fms = [10.0, 1e3, 1e4]

    cf = 1e-9
    f_rc = 1.0 / (2.0 * np.pi * rf * cf)
    fc, vals = _pll_phase_noise(K, kvco, cf, fms, sf=sf)
    ## the corner comes from the multiplier and must match the circuit
    assert abs(fc / (kvco * K / 2.0) - 1.0) < 2e-3, (fc, kvco * K / 2.0)

    pred = [sf * (1.0 + 2.0 * fc / f_rc) / (4.0 * np.pi ** 2 * (fm ** 2 + fc ** 2))
            for fm in fms]
    ratios = [v / p for v, p in zip(vals, pred)]
    ## THE SHAPE is the strong claim: one constant across three decades,
    ## spanning the flat region, the corner and the -20 dB/decade roll-off.
    assert max(ratios) / min(ratios) - 1.0 < 1e-4, ratios
    for r in ratios:
        assert abs(r - 1.0) < 5e-4, ratios

    ## and the second-pole term is a TERM, not a fudge: ten times the filter
    ## corner divides it by ten.
    cf2 = 1e-10
    f_rc2 = 1.0 / (2.0 * np.pi * rf * cf2)
    fc2, vals2 = _pll_phase_noise(K, kvco, cf2, [1e3], sf=sf)
    base2 = sf / (4.0 * np.pi ** 2 * (1e3 ** 2 + fc2 ** 2))
    resid2 = vals2[0] / base2 - 1.0
    assert abs(resid2 / (2.0 * fc2 / f_rc2) - 1.0) < 0.02, \
        'the second-pole term must scale as 2 f_c/f_RC: %.4e vs %.4e' \
        % (resid2, 2.0 * fc2 / f_rc2)


def _nu_identity(kind, n, nonuniform):
    """pnoise of the cyclostationary IDENTITY pair (see
    `test_pnoise_cyclostationary_is_the_stationary_fold_of_the_same_physics...`):
    'A' a stationary source through a multiplier, 'B' the same physics as a
    cyclostationary source.  Optionally on a 3:1 non-uniform grid."""
    import warnings as _w
    T = 1e-6
    c = SubCircuit()
    for nd in ('lo', 'mid', 'out'):
        c.add_node(nd)
    c['vlo'] = VSin('lo', gnd, va=1.0, vo=0.3, freq=1.0 / T)
    if kind == 'A':
        c.add_node('n'); c['xi'] = IS('n', gnd, i=0.0, noisePSD=1.0)
        c['Rn'] = R('n', gnd, r=2.0)
        c['M1'] = _NuMult('mid', gnd, 'n', gnd, 'lo', gnd, k=0.5)
    else:
        c['src'] = _NuModNoise('mid', gnd, 'lo', gnd, k=1.0)
    c['Rm'] = R('mid', gnd, r=1.0)
    c['M2'] = _NuMult('out', gnd, 'mid', gnd, 'lo', gnd, k=0.3)
    c['Ro'] = R('out', gnd, r=1.0); c['Co'] = C('out', gnd, c=0.2e-6)
    grid = None
    if nonuniform:
        w = 1.0 + 0.5 * np.sin(2.0 * np.pi * np.arange(n) / n)
        grid = w / w.sum()
    pss = PSS(c, method='gear', reltol=1e-10)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T, timestep=T / n, grid=grid, maxiterations=40,
                  break_events=False)
    assert pss.converged
    o = [str(x) for x in c.nodes].index('out')
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        s, _ = PAC(c, toolkit=circuit.numeric).pnoise(
            pss, 0.05 / T, o, maxsidebands=16, cyclostationary=(kind == 'B'))
    return float(np.real(s)), pss


def test_pnoise_is_correct_on_a_non_uniform_grid_with_trapezoid_period_weights():
    """⚠⚠ `PSS.solve(grid=...)` accepts a non-uniform grid, and until 2026-09-19
    pnoise on one was WRONG BY TENS OF PERCENT AND DID NOT CONVERGE.

    Every Fourier coefficient the folds take over the period is an INTEGRAL,
    and it was evaluated as an index DFT -- `fft/N`, and `/ N` on the adjoint
    inject.  On a uniform grid that IS the periodic trapezoid rule (spectral);
    on a 3:1 grid it is not a quadrature at all: measured a flat 59 % error in
    the first harmonic, stationary pnoise converging 30 % off (2.971e-01 /
    2.963e-01 / 2.959e-01 at N = 200/400/800, rate 1.00), and the
    cyclostationary identity B/A stuck at 0.7268 / 0.7263 / 0.7261.

    FIXED with `PSS._period_quadrature`: trapezoid weights
    `(h_{n-1}+h_n)/(2T)` at the TRUE times, used by the adjoint inject (plain,
    dirk and full paths), the forward `PAC.solve`, and every harmonic of `CY`.
    Measured after: B/A = 0.99995733 / 0.99998933 / 0.99999733 (rate 4.00 --
    the gear samples' own order; a rectangle `h_n/T` measured first order and
    is the bottleneck for every method).  ⚠ Uniform grids take the ORIGINAL
    expressions (the helper returns None), so they are bit-identical.
    ⚠ A derived grid is now CORRECT, not competitive: uniform reads 3e-10
    where this reads 1e-05 at N = 800, and second order caps radau.
    ⚠ The dual-consistency identity closes for ANY weights used on both
    sides, so it could never have caught this -- only a known answer can.
    ⚠ `sampled_noise` and `covariance` were never affected: their time sums
    carry each step's own `h`.
    """
    ref, pss_u = _nu_identity('A', 400, False)
    assert pss_u._period_quadrature(pss_u.factored_period()) is None, \
        'a uniform grid must keep the index DFT (bit-identical results)'
    errs = []
    for n in (200, 400):
        a, pss_n = _nu_identity('A', n, True)
        b, _p = _nu_identity('B', n, True)
        assert pss_n._period_quadrature(pss_n.factored_period()) is not None
        ## the stationary route lands on the uniform answer (was 30 % off) ...
        assert abs(a / ref - 1.0) < 5e-4, (n, a, ref)
        ## ... and the identity holds (was 0.727, at every N)
        errs.append(abs(b / a - 1.0))
        assert errs[-1] < 1e-4, (n, b / a)
    ## second order: the quadrature reaches the integrator's order
    assert errs[0] / errs[1] > 3.0, errs


def test_a_weightless_entry_cannot_fail_a_component_into_the_sign_blind_route():
    """⚠ THE SILENT FALLBACK, twice.  A flicker component whose entries do not
    share one exponent is evaluated per band from sqrt(PSD) -- sign-blind.  The
    exponent is fitted from differences of `CY`, so its noise goes as
    1/weight, and an orbit sample at Vds ~ 0 always supplies a light entry:

        200 points,  amp 0.375:  weight 1.9e-12, exponent off by 2.2e-09
        1000 points, amp 0.25:   weight 7.1e-09, exponent off by 8.1e-09

    The first failed a 1e-9 tolerance; my repair (entries above 1e-9 vote)
    was walked through by the second, found by the peer as a fixed commit
    reproducing the OLD one to the last bit at two amplitudes.  No weight
    cut-off separates noise from signal; the COST of the wrong exponent does.
    And where the fallback does happen to a component that stated its sign,
    it now says so.
    """
    import warnings as _w
    from pycircuit.circuit.shooting._noise_components import (
        uniform_exponent, warn_signed_unused)
    B = np.zeros((3, 2, 2), dtype=complex)
    EF = np.ones((3, 2, 2))
    B[0, 0, 0] = 7.2e-20
    B[1, 0, 0] = 1.4e-31
    EF[1, 0, 0] = 1.0 - 2.2e-9                    # the 200-point offender
    assert uniform_exponent(B, EF) == 1.0
    B[2, 0, 0] = 7.1e-9 * 7.2e-20
    EF[2, 0, 0] = 1.0 + 8.1e-9                    # the 1000-point offender
    assert uniform_exponent(B, EF) == 1.0
    ## the reference is the LARGEST entry's exponent, wherever it sits
    assert uniform_exponent(B[::-1], EF[::-1]) == 1.0
    ## presence: a genuinely different exponent still fails, even when light
    EF[2, 0, 0] = 2.0
    assert uniform_exponent(B, EF) is None
    B[2, 0, 0] = 1e-12 * 7.2e-20                  # too light to cost 1e-9
    assert uniform_exponent(B, EF) == 1.0

    ## the visible fallback
    def model():
        pass
    model.amplitude = {('M1',): np.zeros((3, 2, 1))}
    B[2, 0, 0] = 1e-3 * 7.2e-20
    model.flicker = [(('M1',), B, EF)]
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        warn_signed_unused(model, 'here')
    assert len(rec) == 1 and 'SIGN-' in str(rec[0].message), rec
    EF[2, 0, 0] = 1.0
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        warn_signed_unused(model, 'here')
    assert not rec


def test_a_library_mosfets_signed_flicker_is_exact_under_a_periodic_fold():
    """⚠ THE CIRCUIT-LEVEL GATE the six converted library models did not have:
    a device whose current REVERSES under a periodic fold, against an answer
    that is exact by construction -- no reference simulator, nothing fitted.

    MOS level 1 track-and-hold; the drain current crosses zero twice per period
    while the switch conducts.  With `af = 2` the device's signed 1/f amplitude
    is EXACTLY `kappa ids(t)`, so the same physics is a CONSTANT 1/f source
    times the SENSED drain current (a CCVS on a series sense branch, into a
    multiplier) -- a stationary source through an LPTV gain, which the fold gets
    right with no sign to lose.  pnoise at 1e-3 f0, over the clock amplitude:

        ck amp   device, signed / exact    device, sign-blind / exact
        0.90        1.000113                   13.2
        0.60        1.000047                   22.9
        0.45        1.000028                   33.8
        0.30        1.000050                   58.1

    ⚠ LIVENESS IS ASSERTED, because this control nearly fooled me: a harness
    print showed the sensed current as [0, 0] (my indexing), and had that been
    true the multiplier would contribute NOTHING and "B/A = 1" would only have
    meant flicker was negligible.  Measured instead: the sensed current spans
    -7.3e-05 .. +1.2e-05 A and flicker is 44 % of the exact total.
    """
    import warnings as _w
    import pycircuit.circuit.elements_hdl as eh
    from pycircuit.circuit.elements import CCVS
    circuit.default_toolkit = circuit.numeric
    F = 1e5
    T = 1.0 / F
    mos = dict(vto=0.4, kp=2e-4, w=10e-6, l=1e-6, af=2.0)
    kf = 1e-24

    ## kappa from the element itself: W[d] = kappa ids / sqrt(f)
    el = eh.MosLevel1Hdl('d', 'g', 's', 'b', kf=kf, **mos)
    x = np.zeros(el.n)
    x[0], x[1] = 0.05, 1.5
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        wd = float(np.real(el.noise_amplitudes(x, 2 * np.pi * 10.0)[0, 0]))
        ids = float(np.real(np.asarray(el.i(x)).ravel()[0]))
    kappa = wd * np.sqrt(10.0) / ids

    def build(exact, gain=None):
        c = SubCircuit()
        for nd in ('in', 'ck', 'dd', 'out', 'x'):
            c.add_node(nd)
        c['Vin'] = VSin('in', gnd, vo=0.3, va=0.2, freq=F, phase=0.0)
        c['Vck'] = VSin('ck', gnd, vo=0.9, va=0.6, freq=F, phase=90.0)
        ## the sense branch is in BOTH circuits: they are the same orbit
        c['Hs'] = CCVS('in', 'dd', 'x', gnd, r=1.0)
        c['M1'] = eh.MosLevel1Hdl('dd', 'ck', 'out', gnd,
                                  kf=(0.0 if exact else kf), **mos)
        if exact:
            c.add_node('n')
            c['xi'] = _Flicker('n', gnd, i=0.0, noisePSD=1.0, fref=1.0)
            c['Rn'] = R('n', gnd, r=1.0)
            c['Mx'] = _NuMult('dd', 'out', 'n', gnd, 'x', gnd,
                              k=(kappa if gain is None else gain))
        c['Ch'] = C('out', gnd, c=100e-12)
        return c

    def pn(c, blind=False):
        pss = PSS(c, method='gear', reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res = pss.solve(period=T, timestep=T / 200, maxiterations=60)
            assert pss.converged
            pac = PAC(c, toolkit=circuit.numeric)
            if blind:
                _noise_seam(pac, signed_amplitudes=lambda *a_, **k_: {})
            o = [str(n) for n in c.nodes].index('out')
            s, _nsb = pac.pnoise(pss, 1e-3 * F, o, maxsidebands=40,
                                 cyclostationary=True)
        vx = np.asarray(res['tpss'].v('x'), dtype=float).ravel()
        return float(np.real(s)), vx

    a, vx = pn(build(True))
    thermal, _v = pn(build(True, gain=0.0))
    ## liveness: the current REVERSES, and the flicker is a real share of A
    assert vx.min() < -1e-5 and vx.max() > 1e-6, (vx.min(), vx.max())
    assert int(np.sum(np.diff(np.sign(vx)) != 0)) >= 2
    assert 1.0 - thermal / a > 0.2, (thermal, a)
    b, _v = pn(build(False))
    bb, _v = pn(build(False), blind=True)
    assert abs(b / a - 1.0) < 5e-4, b / a
    assert bb / a > 5.0, bb / a


def test_am_pm_noise_runs_with_its_default_sideband_count():
    """`am_pm_noise`'s default `maxsidebands=None` takes every sideband pair
    `carrier -+ p` the grid's Nyquist allows, `|p| <= N//2 - |carrier|`.  ⚠ It
    took `N//2`, asked for sideband `carrier + N//2`, and raised "above the
    grid's Nyquist" for every carrier >= 1 (found by the conventions review,
    2026-09-28); every call in the suite passed `maxsidebands`."""
    import warnings
    from pycircuit.circuit.elements import Diode
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['vs'] = VSin(1, gnd, vac=1.0, va=2.0, freq=1e6, phase=20)
    c['R'] = R(1, 2, r=1e4)
    c['D'] = Diode(2, gnd)
    c['C'] = C(2, gnd, c=1e-12)
    T = 1e-6
    pss = PSS(c, method='gear', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 40, maxiterations=40)
        pac = PAC(c, toolkit=circuit.numeric)
        N = len(pss._adjoint_host().factored_period().steps)
        default = pac.am_pm_noise(pss, 0.13 / T, 2, carrier=1)
        explicit = pac.am_pm_noise(pss, 0.13 / T, 2, carrier=1,
                                   maxsidebands=N // 2 - 1)
    assert default[:2] == explicit[:2] and default[2] == explicit[2], (default, explicit)
    assert max(abs(p_) for p_ in default[2]) == N // 2 - 1
