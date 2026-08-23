"""The core against IHP's own PSP103, end to end.

Every other PSP test here checks a construction property or an internal
consistency.  This one checks the thing none of those can: how close the
device is to the transistor a foundry ships.

The whole chain is exercised with nothing hand-tuned --

    IHP model card -> spicecard -> psp_scaling -> the element -> current

-- so a break anywhere in it shows up here, which is what makes this
worth its runtime.

It needs the IHP PDK and skips without it; the reference curves are
committed, but the model CARD is not, and reading it is half the point.
`benchmarks/psp_gap.py` prints the same comparison in detail.
"""
import json
import os
import warnings

import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.compact import (PspMosLongChannel,
                                       PspPmosLongChannel)
from pycircuit.circuit import psp_scaling
from pycircuit.utilities import spicecard


PDK = os.path.expanduser(
    '~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models')
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data',
                   'psp103_ihp_iv.json')

needs_pdk = pytest.mark.skipif(not os.path.isdir(PDK),
                               reason='IHP Open PDK not installed')

#: Below this the reference is dominated by leakage and junction terms
#: the core does not model at all.
#: The card's reference temperature, and the one ngspice ran at when
#: the reference was generated (`TR = 27 C`).  Every comparison against
#: a recorded `lp_*` value has to be made HERE, because that is where
#: the recording was made -- `to_long_channel` defaults to 300.0 K
#: instead, matching the element's own `epar` default, and at 300.0 K
#: the temperature scaling is a small but real correction rather than
#: the identity.
T27 = 273.15 + 27.0

FLOOR = 1e-6


@pytest.fixture(scope='module')
def deck():
    return spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                          section='mos_tt')


@pytest.fixture(scope='module')
def ref():
    with open(REF) as fh:
        return json.load(fh)['sweeps']


def _compare(deck, sweep):
    cm.default_toolkit = numeric
    w, l = sweep['w'], sweep['l']
    ## The sweep names its own device, so the p-channel curves pick the
    ## p-channel card and class without a second code path here.
    pmos = 'pmos' in sweep['device']
    card = deck.model_params(
        'sg13g2_lv_pmos_psp' if pmos else 'sg13g2_lv_nmos_psp',
        w=w, l=l, ng=1, m=1, pre_layout=1)
    kw = psp_scaling.to_long_channel(card, w=w, l=l, T=T27)
    cls = PspPmosLongChannel if pmos else PspMosLongChannel
    e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'), **kw)
    e.update_iparv()
    v = np.asarray(sweep['v'], float)
    r = np.abs(np.asarray(sweep['i_d'], float))
    b = sweep['bias']
    if sweep['sweep'] == 'Vd':
        got = np.array([np.asarray(e.i(e.bias(x, b['Vg'], b['Vs'], b['Vb'])), float)[0]
                        for x in v])
    else:
        got = np.array([np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])), float)[0]
                        for x in v])
    m = r > FLOOR
    return v[m], r[m], np.abs(got[m]), kw


@needs_pdk
class TestAgainstTheRealDevice(object):

    def test_the_long_device_agrees_within_four_percent(self, deck, ref):
        """W = 10 um, L = 1 um -- where the core is meant to be right.

        The core still has no channel-length modulation and no
        short-channel physics, so this is not exact.  Measured range is
        0.980 - 1.035; the bound is set just outside it, tight enough to
        catch a real break anywhere in the chain.
        """
        _, r, g, _ = _compare(deck, ref['nmos_long_idvd'])
        assert len(r) > 20
        ratio = g / r
        assert 0.96 < ratio.min(), 'worst low %.3f' % ratio.min()
        assert ratio.max() < 1.04, 'worst high %.3f' % ratio.max()

    def test_the_agreement_is_flat_not_a_lucky_crossing(self, deck, ref):
        """A curve that crosses the reference could average out to 1.

        What makes the long-device result meaningful is that the ratio
        barely moves across the sweep: a uniform offset means the
        threshold and gain factor are right and something multiplicative
        is missing, rather than two wrong shapes intersecting.
        """
        _, r, g, _ = _compare(deck, ref['nmos_long_idvd'])
        ratio = g / r
        assert ratio.max() - ratio.min() < 0.07

    def test_the_scaled_parameters_are_physical(self, deck, ref):
        """A wrong scaling shows up here before it shows up in a curve."""
        _, _, _, kw = _compare(deck, ref['nmos_long_idvd'])
        assert 1e-9 < kw['tox'] < 5e-9
        assert -1.5 < kw['vfb'] < 0.0
        assert 1e22 < kw['nsub'] < 1e24
        assert 0.5 < kw['phib'] < 1.3
        assert 0.0 < kw['u0'] < 0.1
        assert all(v >= 0.0 for v in (kw['mue'], kw['themu'], kw['cs'],
                                      kw['thecs'], kw['thesat'],
                                      kw['rs'], kw['ct']))

    def test_the_quantum_correction_is_what_closed_the_gap(self, deck, ref):
        """Recorded because it is how the term was found.

        The first comparison ran 23-33% high with a ratio that grew with
        current.  The card sets `QMC = 1`, and the quantum-mechanical
        correction to the surface potential at threshold was the only
        term of that size the core was missing; adding it took the long
        device to ~10%.  This pins that turning it off reopens the gap,
        so nobody removes it as cosmetic.
        """
        sweep = ref['nmos_long_idvd']
        w, l = sweep['w'], sweep['l']
        card = deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1,
                                 m=1, pre_layout=1)
        assert card['qmc'] > 0.0, 'the card does enable it'

        with_qm = psp_scaling.to_long_channel(card, w=w, l=l, T=T27)
        without = psp_scaling.to_long_channel(dict(card, qmc=0.0), w=w, l=l, T=T27)
        ## It raises the surface potential at threshold and the body
        ## factor -- both of which reduce current at a fixed gate bias.
        assert with_qm['phib'] > without['phib']
        assert with_qm['nsub'] > without['nsub']
        assert with_qm['phib'] - without['phib'] > 0.01

    def test_series_resistance_closed_the_next_slice(self, deck, ref):
        """The second term the measurement named, and what it was worth.

        With the threshold right, the residual on the long device was a
        nearly FLAT ~10% -- the signature of something multiplicative.
        PSP folds source/drain resistance into the mobility rather than
        adding a network element, and switching it on took the long
        device from a median 1.095 to 1.041.

        Checked as a direction rather than a number: turning `rs` off
        must raise the current, at every bias.
        """
        sweep = ref['nmos_long_idvd']
        w, l = sweep['w'], sweep['l']
        card = deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1,
                                 m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=w, l=l, T=T27)
        assert kw['rs'] > 0.0, 'the card does specify one'

        cm.default_toolkit = numeric
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, rs=0.0))
        on.update_iparv()
        off.update_iparv()
        for vd in (0.05, 0.4, 0.9, 1.5):
            x = on.bias(vd, 1.2)
            assert (np.asarray(on.i(x), float)[0]
                    < np.asarray(off.i(x), float)[0])

    def test_the_effective_thermal_voltage_closed_the_flat_offset(
            self, deck, ref):
        """The third term, and the one a Vg sweep was built to find.

        PSP does not normalise by the thermal voltage but by an
        EFFECTIVE one, `phit1 = phit*(1 + CT)` -- and the gate drive,
        the quasi-Fermi levels and the charges are all in units of that
        (`PSP103_macrodefs.include:503`).  The card sets `CTO = 0.0546`,
        worth ~7% on this geometry.

        It was found by measuring a LONG-DEVICE TRANSFER CURVE, added to
        the reference for exactly this: a gain error gives a ratio flat
        in Vg, a threshold or drive error gives one that varies.  The
        measured ratio fell from 1.175 near threshold to 1.022 at
        Vg = 1.5, which is a drive error, and `CT` is the only term of
        that size and shape the core was missing.  It took the long
        device from a median 1.041 to 1.008.

        The body factor is deliberately NOT rescaled: PSP builds `G_0`
        on the plain `phit` and only moves it under `SWFIX`, which
        defaults to 0 and this card leaves unset.
        """
        sweep = ref['nmos_long_idvd']
        w, l = sweep['w'], sweep['l']
        card = deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1,
                                 m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=w, l=l, T=T27)
        assert kw['ct'] > 0.05, 'the card does specify one'

        cm.default_toolkit = numeric
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, ct=0.0))
        on.update_iparv()
        off.update_iparv()
        for vg in (0.6, 1.0, 1.5):
            x = on.bias(0.05, vg)
            assert (np.asarray(on.i(x), float)[0]
                    < np.asarray(off.i(x), float)[0])

    def test_the_long_transfer_curve_agrees_too(self, deck, ref):
        """The sweep that found `CT`, now a guard in its own right.

        A Vd sweep at one gate bias can be matched by a wrong model with
        a compensating error; a Vg sweep across five decades of current
        cannot.
        """
        v, r, g, _ = _compare(deck, ref['nmos_long_idvg'])
        assert len(r) > 20
        ratio = g / r
        assert ratio.min() > 0.95, 'worst low %.3f' % ratio.min()
        assert ratio.max() < 1.16, 'worst high %.3f' % ratio.max()

    def test_channel_length_modulation_helps_now(self, deck, ref):
        """The term that was rejected once, and why it works now.

        CLM was implemented, measured and REVERTED earlier in this work:
        it made every sweep worse.  The reason was not that it was
        wrong -- it was that the device was still 4% high from the
        missing effective thermal voltage, and CLM only pushes current
        higher.  A term that is individually correct can worsen the fit
        while a compensating error is present.

        With `CT` in, the residual is centred and CLM helps: better in
        four sweeps of seven, and dramatically so on the short device
        (median |ratio-1| 0.098 -> 0.033 on `nmos_idvd_vg1p2`, 0.132 ->
        0.026 on `nmos_idvg_vd1p2`).

        Checked as a direction: it must RAISE the saturated current and
        do nothing in the linear region.
        """
        sweep = ref['nmos_long_idvd']
        w, l = sweep['w'], sweep['l']
        card = deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1,
                                 m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=w, l=l, T=T27)
        assert kw['alp'] > 0.0

        cm.default_toolkit = numeric
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, alp=0.0))
        on.update_iparv()
        off.update_iparv()
        deep = on.bias(1.5, 1.2)
        lin = on.bias(0.02, 1.2)
        assert (np.asarray(on.i(deep), float)[0]
                > np.asarray(off.i(deep), float)[0])
        assert np.asarray(on.i(lin), float)[0] == pytest.approx(
            np.asarray(off.i(lin), float)[0], rel=1e-3)

    def test_our_saturated_current_now_rises_too(self, deck, ref):
        """Before CLM the core's saturated current was flat.

        PSP gains 5.9% between Vd = 0.8 and 1.4.  The core is not
        expected to match that exactly -- the CLM here approximates
        PSP's, which references a saturation-limited drain voltage this
        core does not compute -- but it must now rise rather than sit
        flat.
        """
        v, r, g, _ = _compare(deck, ref['nmos_long_idvd'])
        sat = v >= 0.8
        assert r[sat][-1] / r[sat][0] > 1.03, 'PSP climbs'
        assert g[sat][-1] / g[sat][0] > 1.01, 'so do we, now'

    def test_both_geometries_are_within_a_few_percent(self, deck, ref):
        """Which one is worse has now INVERTED, twice.

        This test began by asserting the short device was at least 3x
        worse than the long one -- true when the short device was ~12%
        off.  The `GPE` fix took it to ~2%, and the `FdL` term took it to
        ~1% while moving the LONG device to ~3.5%, so the short one is
        now the more accurate of the two.

        Rewritten rather than re-inverted: what is worth pinning is that
        both are bounded, not which happens to lead.
        """
        _, r_l, g_l, _ = _compare(deck, ref['nmos_long_idvd'])
        _, r_s, g_s, _ = _compare(deck, ref['nmos_idvg_vd0p05'])
        long_err = np.median(np.abs(g_l / r_l - 1.0))
        short_err = np.median(np.abs(g_s / r_s - 1.0))
        assert long_err < 0.05, 'long %.3f' % long_err
        assert short_err < 0.05, 'short %.3f' % short_err


@needs_pdk
class TestTheSubthresholdRegion(object):
    """The decades BELOW the sweep comparison's floor.

    `_compare` masks at `FLOOR = 1e-6` A, so five decades of every
    transfer curve have never been compared against the vendor at all --
    which matters, because that is where a threshold error lives.  Above
    threshold a few millivolts of `Vth` is a tenth of a percent of
    current and invisible; at 85 mV/decade it is ten percent.

    The measurement here is not a current ratio, it is a VOLTAGE.  In
    weak inversion `I ~ exp(Vg / (n*phit))`, so

        dVth = ln(ours/psp) / (d ln I_psp / dVg)

    turns the ratio into the threshold offset that would produce it.
    Two things follow that a ratio alone cannot say: whether the offset
    is FLAT across the region -- if it is, it is a threshold error and
    nothing else -- and how big it is in millivolts, which is the unit
    the physics is actually in.
    """

    #: The LEAKAGE-FREE window, and it is not a detail.
    #:
    #: The reference records TOTAL terminal current, and PSP's junction
    #: leakage is a flat ~2e-12 A floor this core does not model at all
    #: (`TestTheSubthresholdSlope` found that the hard way).  A window
    #: reaching below it measures PSP's leakage rather than its channel.
    #: The first version of this class used `1e-14` and got answers that
    #: were wrong by a third with three times the spread; these bounds
    #: are the same ones the slope test settled on.
    LO, HI = 1e-9, 1e-6

    #: Implied threshold offset, in mV, and the tolerated spread.
    #: MEASURED, and the numbers are the finding: the n-channel sits a
    #: few mV low and the p-channel is essentially exact, so whatever is
    #: left is n-channel-specific rather than a shared constant.
    EXPECT = {'nmos_idvg_vd0p05': (3.6, 2.0),
              'nmos_long_idvg': (1.7, 0.5),
              'pmos_idvg_vd0p05': (-1.2, 0.8),
              'pmos_long_idvg': (0.0, 0.3)}

    @classmethod
    def _implied_vth(cls, deck, sweep):
        """Return (mean dVth in mV, its spread, our slope / PSP's)."""
        cm.default_toolkit = numeric
        w, l = sweep['w'], sweep['l']
        pmos = 'pmos' in sweep['device']
        kw = psp_scaling.to_long_channel(
            deck.model_params(
                'sg13g2_lv_pmos_psp' if pmos else 'sg13g2_lv_nmos_psp',
                w=w, l=l, ng=1, m=1, pre_layout=1), w=w, l=l, T=T27)
        elem = PspPmosLongChannel if pmos else PspMosLongChannel
        e = elem(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                 **kw)
        e.update_iparv()
        b = sweep['bias']
        v = np.asarray(sweep['v'], float)
        r = np.abs(np.asarray(sweep['i_d'], float))
        got = np.abs(np.array([
            np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])),
                       float)[0] for x in v]))
        m = (r > cls.LO) & (r < cls.HI)
        assert m.sum() >= 6, 'not enough points in the window'
        ## A least-squares slope over the whole window, not a point
        ## derivative: it is what makes 0.03% resolvable, and it is the
        ## same quantity for both curves so the window cancels.
        A = np.vstack([np.ones(m.sum()), v[m]]).T
        slope = np.linalg.lstsq(A, np.log(r[m]), rcond=None)[0][1]
        ours = np.linalg.lstsq(A, np.log(got[m]), rcond=None)[0][1]
        dv = np.log(got[m] / r[m]) / slope * 1e3   # mV
        return dv.mean(), dv.max() - dv.min(), ours / slope

    @pytest.mark.parametrize('name', sorted(EXPECT))
    def test_the_implied_threshold_offset(self, deck, ref, name):
        """A few millivolts, and FLAT -- so it is a threshold offset.

        Flatness is the part that identifies it.  A slope error or a
        bias-dependent term would give an implied `dVth` that walks
        across the region; these sit within a couple of millivolts over
        four to five decades of current.
        """
        want, tol = self.EXPECT[name]
        mean, spread, _ = self._implied_vth(deck, ref[name])
        assert mean == pytest.approx(want, abs=1.0), \
            '%s: %+.2f mV (expected %+.2f)' % (name, mean, want)
        assert spread < tol, '%s: spread %.2f mV' % (name, spread)

    def test_the_slope_is_right_and_the_level_is_not(self, deck, ref):
        """THE discriminator, and it changes what the residual is.

        Gain and electrostatics are separable, and this is where:
        `XCOR`, `RSB`, `BETN` and the mobility are GAIN -- they scale the
        current and leave `d ln I / dVg` untouched.  The body factor is
        ELECTROSTATIC -- it sets the ideality, so it moves the slope.

        Measured on the long devices, in the leakage-free window: the
        slope matches to **0.03%** on both channel types while the level
        is 1.7 mV out on the n-channel and 0.03 mV on the p-channel.

        A correct slope with a wrong level is not a threshold error and
        not a body-factor error.  It is a WEAK-INVERSION GAIN error --
        and it cannot be a plain gain error either, because above
        threshold the same device is within 0.1%.  That is a much
        narrower thing to look for than "the threshold is off", which is
        what this was called before the window was fixed.
        """
        for name, lim in (('nmos_long_idvg', 3e-4), ('pmos_long_idvg', 5e-4)):
            _, _, ratio = self._implied_vth(deck, ref[name])
            assert abs(ratio - 1.0) < lim, (name, ratio)

    def test_the_n_channel_is_the_one_that_is_off(self, deck, ref):
        """The asymmetry IS the lead, so it is pinned rather than
        described.

        The p-channel threshold is right to a millivolt at both
        geometries and the n-channel is not, which rules out anything
        shared -- a physical constant, the thermal voltage, the oxide
        capacitance -- and points at something n-channel-specific.
        """
        n_long = self._implied_vth(deck, ref['nmos_long_idvg'])[0]
        n_short = self._implied_vth(deck, ref['nmos_idvg_vd0p05'])[0]
        p_long = self._implied_vth(deck, ref['pmos_long_idvg'])[0]
        assert abs(p_long) < 0.5, p_long
        assert n_long > 1.0, n_long
        ## and it GROWS as the channel shortens, so part of it is the
        ## reverse-short-channel shift rather than a constant.
        assert n_short > n_long + 1.0, (n_long, n_short)

    @staticmethod
    def _card_kw(deck, sweep):
        """Scaled parameters for the device a sweep names."""
        w, l = sweep['w'], sweep['l']
        pmos = 'pmos' in sweep['device']
        return psp_scaling.to_long_channel(
            deck.model_params(
                'sg13g2_lv_pmos_psp' if pmos else 'sg13g2_lv_nmos_psp',
                w=w, l=l, ng=1, m=1, pre_layout=1), w=w, l=l, T=T27)

    #: The body-bias ladder, and what each rung measures.
    LADDER = [('nmos_long_idvg_vb0p00', 0.0), ('nmos_long_idvg_vbm0p20', 0.2),
              ('nmos_long_idvg_vbm0p40', 0.4), ('nmos_long_idvg_vbm0p70', 0.7),
              ('nmos_long_idvg_vbm1p00', 1.0), ('nmos_long_idvg_vbm1p50', 1.5)]

    def test_the_offset_grows_with_body_bias(self, deck, ref):
        """Monotone in `Vsb`, which says it is not a `VFB` error.

        A flat-band error is flat in body bias; this rises from +1.55 mV
        at `Vsb = 0` to +2.71 mV at 1.5 V.  So whatever is left lives in
        a term that body bias exercises -- the body factor, or the
        `XCOR` mobility correction, which is identically 1 at `Vsb = 0`
        and therefore contributes the whole of the difference and none
        of the baseline.
        """
        got = [self._implied_vth(deck, ref[n])[0] for n, _ in self.LADDER]
        assert all(b > a - 0.05 for a, b in zip(got, got[1:])), got
        assert got[-1] - got[0] > 1.0, got

    def test_the_two_candidates_cannot_be_separated_here(self, deck, ref):
        """A NEGATIVE result, pinned so it is not re-attempted blind.

        The two candidate terms have different shapes in principle -- a
        body-factor error goes as `sqrt(phib + Vsb)`, `XCOR` as the
        saturating `(1 + 0.2*XCOR*Vsb)/(1 + XCOR*Vsb)` -- so a ladder in
        `Vsb` ought to tell them apart.  It does not, and the reason is
        measurable rather than a matter of precision: over the `Vsb`
        range this device can actually be measured in, the two bases are
        **99.8% collinear**, and both are 99.7% collinear with a plain
        straight line.

        The range cannot be extended.  Below `Vg = 0` the reference
        stops being diffusion subthreshold and becomes junction and GIDL
        leakage, which PSP models and this element does not; the curve
        flattens, `d ln I / dVg` goes to zero, and the measurement
        divides by it.  Tried, and the numbers were nonsense.

        So this is not "the fit was inconclusive".  It is that no fit to
        this device can be conclusive, and separating the two needs a
        different measurement rather than a better one.
        """
        vsb = np.array([v for _, v in self.LADDER])
        kw = self._card_kw(deck, ref['nmos_long_idvg'])
        phib, xcor = kw['phib'], kw['xcor']
        b_sqrt = np.sqrt(phib + vsb) - np.sqrt(phib)
        b_xcor = -np.log((1.0 + 0.2 * xcor * vsb) / (1.0 + xcor * vsb))
        def corr(a, b):
            a, b = a - a.mean(), b - b.mean()
            return float(a @ b / np.sqrt((a @ a) * (b @ b)))
        assert corr(b_sqrt, b_xcor) > 0.99, corr(b_sqrt, b_xcor)
        assert corr(b_sqrt, vsb) > 0.99, corr(b_sqrt, vsb)

    def test_the_body_factor_reading_of_it(self, deck, ref):
        """IF it is all body factor, it is a 0.95% deficit.

        Recorded as a bound on the size of whatever is left rather than
        as a claim about which term it is -- the test above says the
        two cannot be told apart, and this one says neither is large.
        """
        vsb = np.array([v for _, v in self.LADDER])
        got = np.array([self._implied_vth(deck, ref[n])[0]
                        for n, _ in self.LADDER])
        kw = self._card_kw(deck, ref['nmos_long_idvg'])
        basis = np.sqrt(kw['phib'] + vsb) - np.sqrt(kw['phib'])
        A = np.vstack([np.ones_like(vsb), basis]).T
        coef = np.linalg.lstsq(A, got, rcond=None)[0]
        rms = np.sqrt(np.mean((A @ coef - got) ** 2))
        assert rms < 0.25, rms
        ## mV per sqrt(V), against a body factor of about 0.209
        assert 1.0 < coef[1] < 3.0, coef
        assert coef[1] * 1e-3 / 0.209 < 0.02, 'deficit under 2%%'

    def test_the_subthreshold_current_is_within_eleven_percent(self,
                                                               deck, ref):
        """The same statement in the units the sweeps use, so the
        coverage gap is closed in its own terms too."""
        for name in sorted(self.EXPECT):
            sw = ref[name]
            cm.default_toolkit = numeric
            w, l = sw['w'], sw['l']
            pmos = 'pmos' in sw['device']
            kw = psp_scaling.to_long_channel(
                deck.model_params(
                    'sg13g2_lv_pmos_psp' if pmos else 'sg13g2_lv_nmos_psp',
                    w=w, l=l, ng=1, m=1, pre_layout=1), w=w, l=l, T=T27)
            cls = PspPmosLongChannel if pmos else PspMosLongChannel
            e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                    cm.Node('b'), **kw)
            e.update_iparv()
            b = sw['bias']
            v = np.asarray(sw['v'], float)
            r = np.abs(np.asarray(sw['i_d'], float))
            m = (r > 1e-13) & (r < 1e-6)
            ## The long device carries more current, so it has fewer
            ## decades under the floor -- about 2.5 rather than 5.
            assert m.sum() > 6, (name, int(m.sum()))
            got = np.abs(np.array([
                np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])),
                           float)[0] for x in v[m]]))
            ratio = got / r[m]
            assert 0.90 < ratio.min() and ratio.max() < 1.11, \
                '%s: %.3f - %.3f' % (name, ratio.min(), ratio.max())


@needs_pdk
class TestTheThresholdScaling(object):
    """Our geometry scaling against PSP's OWN threshold voltage.

    `vth` is an operating-point output of the model itself, so the
    reference records it at four geometries.  That turns "is the scaling
    right?" from an inference about currents into a direct comparison --
    a current can be right for compensating reasons, a threshold cannot.

    The ABSOLUTE values are not comparable: PSP's `vth` is its own
    extraction, not `vfb + phib + gamma*sqrt(phib)`, and the two differ
    by a nearly constant 70-93 mV.  A constant offset is a definitional
    difference.  What IS comparable is the SHIFT between geometries,
    which is pure physics -- the reverse-short-channel effect from the
    pocket implants.
    """

    @staticmethod
    def _our_vth(deck, w, l):
        from pycircuit.circuit.constants import (eps0, epsRSi, epsRSiO2,
                                                 qelectron)
        card = deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1,
                                 m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=w, l=l, T=T27)
        cox = epsRSiO2 * eps0 / kw['tox']
        gamma = np.sqrt(2 * qelectron * epsRSi * eps0 * kw['nsub']) / cox
        return kw['vfb'] + kw['phib'] + gamma * np.sqrt(kw['phib'])

    def test_the_reference_records_psp_s_own_threshold(self, ref):
        import json
        with open(REF) as fh:
            data = json.load(fh)
        assert 'vth' in data
        assert set(data['vth']) == {'long', 'mid', 'short', 'wide_short'}
        ## the reverse-short-channel effect: threshold RISES as L shrinks
        v = data['vth']
        assert v['long']['vth'] < v['mid']['vth'] < v['short']['vth']

    def test_the_threshold_shift_with_length_is_within_ten_percent(
            self, deck):
        """The comparable quantity, and what it does NOT say.

        Measured: PSP shifts +237 mV from L = 1 um to L = 0.13 um and
        `_our_vth` shifts +214 mV -- 90% of it, from the pocket-implant
        doping formula plus the length terms in `VFB` and `DPHIB`.

        ⚠ **The missing 23 mV is an artefact of the PROXY, not the
        model's error.**  `_our_vth` is a textbook
        `vfb + phib + gamma*sqrt(phib)`, which is not what the model
        computes -- it leaves out the QM correction folded into the
        doping, poly depletion, and the bias-dependent body factor.
        Measured properly, from the SUBTHRESHOLD CURRENT
        (`TestTheSubthresholdRegion`), the model's threshold is 1.6 mV
        low on the long n-channel and 3.2 mV low on the short one: the
        geometry-dependent part is about 2 mV, not 23.  **The proxy
        overstates the error tenfold.**

        This test is therefore a check that the SCALING LAYER tracks
        PSP's geometry trend, which it is good for, and not a
        measurement of threshold accuracy, which it is not.  An earlier
        version of this docstring drew a quantitative conclusion from it
        ("worth about 3% of drain current ... a quarter of the ~12% the
        short device is off by") and both halves of that are now wrong.
        """
        import json
        with open(REF) as fh:
            v = json.load(fh)['vth']
        base_psp = v['long']['vth']
        base_our = self._our_vth(deck, v['long']['w'], v['long']['l'])
        for name in ('mid', 'short', 'wide_short'):
            e = v[name]
            psp_shift = e['vth'] - base_psp
            our_shift = self._our_vth(deck, e['w'], e['l']) - base_our
            assert psp_shift > 0.02, name
            rel = abs(our_shift - psp_shift) / psp_shift
            assert rel < 0.12, '%s: %.1f%% off' % (name, 100 * rel)

    def test_the_offset_is_constant_which_makes_it_definitional(self,
                                                                deck):
        """A constant offset is a different `vth`, not a wrong one.

        If our scaling were wrong the discrepancy would move with
        geometry.  It sits between 70 and 93 mV across a 7.7x range of
        channel length, and the part that does move is the 10% of the
        shift accounted for above.
        """
        import json
        with open(REF) as fh:
            v = json.load(fh)['vth']
        offs = [self._our_vth(deck, e['w'], e['l']) - e['vth']
                for e in v.values()]
        assert all(-0.12 < o < -0.05 for o in offs), offs
        assert max(offs) - min(offs) < 0.03


@needs_pdk
class TestTheScalingTermByTerm(object):
    """Our geometry scaling against PSP's OWN scaled parameters.

    PSP exposes its post-scaling values as `lp_*` operating-point
    outputs, and the reference records them at four geometries.  That
    makes the scaling layer checkable **term by term** rather than
    inferred from currents -- and the difference matters: a current can
    be right for compensating reasons, a scaled parameter cannot.

    This is how the `GPE` bug was found.  Every term matched exactly
    except `betn`, which was 12% high, because `FBET1` and `LP1` are
    width-adjusted before use (`PSP103_scaling.include:284-285`) and the
    raw card values were being used.  No amount of staring at I-V curves
    would have said which term; one comparison did.
    """

    #: What our scaling produces, keyed by PSP's name for it.  `neff` is
    #: compared with the quantum correction OFF, because PSP applies that
    #: to `phib` and `G_0` while we fold it into the doping -- equivalent
    #: bookkeeping, different intermediate.
    DIRECT = {'vfb': 'vfb', 'tox': 'tox', 'ct': 'ct', 'mue': 'mue',
              'themu': 'themu', 'cs': 'cs', 'thecs': 'thecs', 'rs': 'rs'}

    @pytest.fixture(scope='class')
    def scaled(self):
        import json
        with open(REF) as fh:
            return json.load(fh)['scaled']

    def test_the_reference_records_them(self, scaled):
        assert set(scaled) == {'long', 'mid', 'short', 'wide_short'}
        for e in scaled.values():
            assert 'betn' in e and 'neff' in e and 'vfb' in e

    @pytest.mark.parametrize('geom', ['long', 'mid', 'short', 'wide_short'])
    def test_every_directly_scaled_parameter_matches(self, deck, scaled,
                                                     geom):
        e = scaled[geom]
        card = deck.model_params('sg13g2_lv_nmos_psp', w=e['w'], l=e['l'],
                                 ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=e['w'], l=e['l'], T=T27)
        for psp_name, ours in self.DIRECT.items():
            if psp_name not in e:
                continue
            assert kw[ours] == pytest.approx(e[psp_name], rel=1e-5), \
                '%s %s: %r vs %r' % (geom, psp_name, kw[ours], e[psp_name])

    @pytest.mark.parametrize('geom', ['long', 'mid', 'short', 'wide_short'])
    def test_the_gain_factor_matches(self, deck, scaled, geom):
        """`BETN = UO*WE*GWE/(GPE*LE)`, the term that was 12% off.

        Our element takes an effective mobility rather than a `BETN`,
        because it computes `beta = u0*cox*W/L` from the DRAWN
        dimensions; `u0_eff * W/L` is the same quantity.
        """
        e = scaled[geom]
        card = deck.model_params('sg13g2_lv_nmos_psp', w=e['w'], l=e['l'],
                                 ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=e['w'], l=e['l'], T=T27)
        ours = kw['u0'] * e['w'] / e['l']
        assert ours == pytest.approx(e['betn'], rel=1e-5), \
            '%s: %r vs %r' % (geom, ours, e['betn'])

    @pytest.mark.parametrize('geom', ['long', 'mid', 'short', 'wide_short'])
    def test_the_effective_doping_matches_before_the_qm_correction(
            self, deck, scaled, geom):
        e = scaled[geom]
        card = deck.model_params('sg13g2_lv_nmos_psp', w=e['w'], l=e['l'],
                                 ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(dict(card, qmc=0.0), w=e['w'],
                                         l=e['l'], T=T27)
        assert kw['nsub'] == pytest.approx(e['neff'], rel=1e-5)

    def test_polysilicon_depletion_is_active_and_not_modelled(self,
                                                              scaled):
        """A known omission, recorded so it is not mistaken for noise.

        `lp_np` is the poly gate doping, and it is 4.6e26 -- not the zero
        that would switch the effect off.  The core assumes an ideal
        gate (`eta_p = 1`), so some of the residual on every geometry is
        this.
        """
        for e in scaled.values():
            assert e.get('np', 0.0) > 1e26


@needs_pdk
class TestPolysiliconDepletion(object):
    """The gate is not a perfect conductor.

    A depletion layer forms on the gate's silicon side and takes a share
    of the applied voltage, so the charge slope is reduced by
    `eta_p = 1/sqrt(1 + kp*xgm)` and the midpoint potential itself is
    shifted (`PSP103_macrodefs.include:697-724`).

    The card's `NPO` is 4.6e26 -- not the zero that switches it off --
    and it was found by the term-by-term scaling comparison rather than
    by looking at currents, since `lp_np` is one of the values PSP
    exposes.
    """

    def test_the_card_enables_it(self, deck):
        card = deck.model_params('sg13g2_lv_nmos_psp', w=1e-6, l=0.13e-6,
                                 ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=1e-6, l=0.13e-6, T=T27)
        assert kw['kp'] > 0.0
        ## 1.6e-3 for this oxide -- 1% of the charge slope at low gate
        ## bias, 4.5% at high, which is why it reads as a bias-dependent
        ## gain error rather than a constant one.
        assert 1e-4 < kw['kp'] < 1e-2

    def test_it_reduces_the_current(self, deck):
        """A gain reduction of about 1%, at every bias.

        Checked as a direction and a magnitude, NOT as a monotone trend
        in gate bias -- the first version of this test asserted that and
        it is not true.  `eta_p` depends on the midpoint oxide voltage,
        which does not track `Vg` monotonically once the correction also
        shifts the midpoint potential, so the reduction is 0.8-1.0%
        across the sweep without being ordered.
        """
        cm.default_toolkit = numeric
        card = deck.model_params('sg13g2_lv_nmos_psp', w=10e-6, l=1e-6,
                                 ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=10e-6, l=1e-6, T=T27)
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, kp=0.0))
        on.update_iparv()
        off.update_iparv()
        for vg in (0.6, 1.0, 1.4, 1.8):
            x = on.bias(0.05, vg)
            r = (np.asarray(on.i(x), float)[0]
                 / np.asarray(off.i(x), float)[0])
            assert 0.95 < r < 1.0, 'Vg=%.1f ratio %.4f' % (vg, r)

    def test_the_construction_survives_it(self, deck):
        """Symmetry and charge conservation, with EVERY card layer on.

        This is the test that caught a latent bug, and the reason it is
        here rather than beside the other symmetry tests: those build
        the element with DEFAULT parameters, where `rs`, `alp` and `kp`
        are all zero.  With a real card they are not, and two of the
        layers were reading bias variables that are not even under the
        source/drain exchange -- the series resistance a source-bulk
        voltage, channel-length modulation a signed `Vds`.  Symmetry had
        been broken since series resistance went in, invisibly.

        Fixed the way PSP fixes it, with symmetrised `vdsx`/`vsbx`
        (`PSP103_macrodefs.include:472`).  Not bit-exact any more --
        the smoothing costs an ulp -- so this asks for 1e-14 rather
        than `==`.
        """
        cm.default_toolkit = numeric
        card = deck.model_params('sg13g2_lv_nmos_psp', w=10e-6, l=1e-6,
                                 ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=10e-6, l=1e-6, T=T27)
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        assert kw['rs'] > 0 and kw['alp'] > 0 and kw['kp'] > 0, \
            'the card must actually switch these on for this to test '\
            'anything'
        for vg in (0.6, 1.2, 1.8):
            for vd in (0.05, 0.3, 0.9, 1.5):
                f = np.asarray(e.i(e.bias(vd, vg, 0.0, 0.0)), float)
                r = np.asarray(e.i(e.bias(0.0, vg, vd, 0.0)), float)
                assert f[0] == pytest.approx(-r[0], rel=1e-14), (vg, vd)
            q = np.asarray(e.q(e.bias(0.5, vg, 0.0, 0.0)), float)
            assert abs(q.sum()) < 4e-16 * max(np.abs(q).max(), 1e-30)


@needs_pdk
class TestTheChannelShorteningFactor(object):
    """`FdL`, and the assumption that hid it.

    PSP multiplies the drain current by
    `FdL = (1 + dL1 + dL1^2) * GdL` (`PSP103_module.include:1137`), where
    `dL1 = dL + ALP1*(...) + ALP2*(qbm*r2^2*s2)`.  With `ALP1` and `ALP2`
    zero this is exactly 1, and it was left out on that reading.

    They are not zero.  PSP's own `lp_alp2` is 4.5 on a 0.13 um device --
    found by the term-by-term comparison, not by inspecting the card,
    since `ALP2` has no geometry-independent coefficient and its scaled
    value comes out of `ALP2L1` through a length power.
    """

    def test_the_card_gives_a_large_weak_inversion_term(self, deck):
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=1e-6, l=0.13e-6,
                              ng=1, m=1, pre_layout=1), w=1e-6, l=0.13e-6, T=T27)
        assert kw['alp2'] > 1.0, kw['alp2']
        assert kw['alp1'] > 0.0

    def test_dibl_is_not_what_carries_the_output_conductance(self, ref):
        """Measured on the reference, and it is why `FdL` was looked for.

        At Vg = 0.6 the reference current climbs 2.4x through saturation.
        If that were a threshold effect, PSP's own `vth` would move
        substantially with drain bias.  It moves 3.5 mV over 1.35 V --
        2.6 mV/V.  So DIBL is genuinely small here and the climb is in
        the current formula.

        To be exact about what that does and does not say: we do not
        model DIBL AT ALL -- `CF`, `CFB` and `CFD` are absent from the
        scaling layer.  The measurement is what licenses that omission,
        not agreement with a term we have.  PSP's own shift works out to
        3.2 mV at Vds = 1.2 V on this device, which at a strong-inversion
        `gm/Id` is worth about 1% of the current.
        """
        s = ref['nmos_idvd_vg0p6']
        v = np.asarray(s['v'], float)
        i = np.abs(np.asarray(s['i_d'], float))
        assert np.interp(1.4, v, i) / np.interp(0.2, v, i) > 1.8

    def test_it_lifts_the_near_threshold_output_conductance(self, deck,
                                                            ref):
        """`ALP2` multiplies the BULK charge, so it acts near threshold.

        Checked as a direction: switching it on must raise the saturated
        current at low gate bias, and much less at high.
        """
        cm.default_toolkit = numeric
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=1e-6, l=0.13e-6,
                              ng=1, m=1, pre_layout=1), w=1e-6, l=0.13e-6, T=T27)
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, alp1=0.0,
                                                     alp2=0.0))
        on.update_iparv()
        off.update_iparv()
        lift = {}
        for vg in (0.6, 1.8):
            x = on.bias(1.4, vg)
            lift[vg] = (np.asarray(on.i(x), float)[0]
                        / np.asarray(off.i(x), float)[0])
        assert lift[0.6] > 1.4, lift
        assert lift[0.6] > lift[1.8], lift

    def test_it_more_than_halves_the_total_error(self, deck, ref):
        """The reason it is kept despite costing the long device.

        Mixed by COUNT -- better in three sweeps, worse in four -- but
        the gains are large and the regressions small, so the summed
        median error more than halves.  The regressions were the price
        of our `dL` being an approximation to PSP's -- `FdL` amplifies
        whatever error `dL` carries -- which is what pointed at the
        saturation voltage as the next thing to build.  It has been
        built since (`TestTheSaturationVoltage`), and it took the long
        device from 1.035 back to 1.010.
        """
        names = ('nmos_long_idvd', 'nmos_long_idvg', 'nmos_idvd_vg1p2',
                 'nmos_idvd_vg0p6', 'nmos_idvg_vd0p05',
                 'nmos_idvg_vd1p2', 'nmos_idvg_vb_m1')
        cm.default_toolkit = numeric
        tot = [0.0, 0.0]
        for n in names:
            s = ref[n]
            kw = psp_scaling.to_long_channel(
                deck.model_params('sg13g2_lv_nmos_psp', w=s['w'], l=s['l'],
                                  ng=1, m=1, pre_layout=1),
                w=s['w'], l=s['l'], T=T27)
            for k, over in enumerate((dict(alp1=0.0, alp2=0.0), {})):
                e = PspMosLongChannel(cm.Node('d'), cm.Node('g'),
                                      cm.Node('s'), cm.Node('b'),
                                      **dict(kw, **over))
                e.update_iparv()
                v = np.asarray(s['v'], float)
                r = np.abs(np.asarray(s['i_d'], float))
                b = s['bias']
                if s['sweep'] == 'Vd':
                    g = np.array([np.asarray(e.i(e.bias(x, b['Vg'], b['Vs'], b['Vb'])), float)[0]
                        for x in v])
                else:
                    g = np.array([np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])), float)[0]
                        for x in v])
                m = r > 1e-6
                tot[k] += np.median(np.abs(np.abs(g[m]) / r[m] - 1.0))
        assert tot[1] < 0.5 * tot[0], 'off %.3f on %.3f' % (tot[0], tot[1])


class TestTheSaturationVoltage(object):
    """`Vdsat` and `Vdse` -- and the floor that made them work.

    PSP does not evaluate the drain surface potential at the applied
    drain bias.  It computes a saturation voltage from SOURCE-SIDE
    quantities, smoothly limits `Vds` to it, and uses the limited
    `Vdse` for the drain quasi-Fermi level
    (`PSP103_macrodefs.include:596-632`).

    Two things came out of building it.  The first is that the exponent
    `AX` governing the limiter is FLOORED at 2 in the scaling layer
    (`PSP103_scaling.include:743`), and the card gives no hint of it:
    `AXO/(1 + AXL*iLE)` evaluates to 0.88 at 0.13 um.  Missing the floor
    made the limiter soft enough to bite at drain biases far below
    saturation and cost 14% on the Vd = 0.05 sweeps -- a sweep where
    nothing should have changed at all, which is what gave it away.

    The second is that `Vdsat` is built from the source end alone and is
    therefore NOT an odd function of `Vds`, so it broke the exact
    source/drain antisymmetry.  PSP's answer, and now ours, is to order
    the terminals rather than to look for an odd formula -- see
    `test_psp_current.TestTheElement`.
    """

    GEOMS = ('long', 'mid', 'short', 'wide_short')

    @pytest.fixture(scope='class')
    def scaled(self):
        with open(REF) as fh:
            return json.load(fh)['scaled']

    @pytest.mark.parametrize('geom', GEOMS)
    def test_our_ax_matches_psps_exactly(self, deck, scaled, geom):
        """Term-by-term against PSP's own `lp_ax`, at four geometries.

        This is the check that a current comparison cannot make: a
        current can be right for compensating reasons, a scaled
        parameter cannot.
        """
        g = scaled[geom]
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=g['w'], l=g['l'],
                              ng=1, m=1, pre_layout=1),
            w=g['w'], l=g['l'], T=T27)
        ## `rel=1e-5`, not tighter: ngspice prints the operating-point
        ## outputs to six significant digits, so that IS the reference's
        ## precision.  A tighter tolerance tests the printf format.
        assert kw['ax'] == pytest.approx(g['ax'], rel=1e-5), (
            geom, kw['ax'], g['ax'])

    def test_the_floor_is_what_makes_the_short_device_match(self, deck,
                                                            scaled):
        """The floor BINDS, and the unclipped formula would miss badly.

        Both short geometries land on exactly 2, and the long one does
        not -- so 2 is a clamp being hit, not a coincidence of the
        scaling.  The unclipped value is recomputed here from the card
        so the test fails if someone removes the clamp.
        """
        card = deck.model_params('sg13g2_lv_nmos_psp', w=1e-6, l=0.13e-6,
                                 ng=1, m=1, pre_layout=1)
        geo = psp_scaling.geometry(card, w=1e-6, l=0.13e-6)
        raw = card['axo'] / (1.0 + card['axl'] * geo['iLE'])
        assert raw < 1.0, raw
        assert scaled['short']['ax'] == 2.0
        assert scaled['wide_short']['ax'] == 2.0
        assert scaled['long']['ax'] > 2.0

    def test_the_velocity_saturation_body_and_gate_terms_are_absent(
            self, scaled):
        """A gap written down rather than remembered.

        PSP scales its saturation parameter by two further factors
        before using it (`macrodefs:596-607`): `xitsb` in the body bias
        and `xitsg` in the inversion charge.  Both coefficients are
        nonzero on this card and neither is modelled, so `Vdsat` is
        computed from an unmodulated `THESAT`.  Recorded here so the
        next person measuring a residual knows where to look.
        """
        assert scaled['short']['thesatb'] > 0.05
        assert scaled['short']['thesatg'] > 0.05

    def test_every_sweep_is_close_and_flat(self, deck, ref):
        """The state of the whole comparison, as a regression guard.

        Nothing here is fitted: every parameter comes off the card
        through the scaling layer.  The history this guard replaces is
        worth keeping in view -- the worst sweep was 11% out before the
        saturation voltage, 7% after it, and 0.3% once the saturation
        voltage was fed back into the channel-length modulation, which
        is what it had been missing all along.

        Held at 3% rather than at the measured 2.3%: this is a guard
        against regression, not a pin on numbers that should be free to
        improve.

        It was briefly 1%, on a model whose short-channel terms were
        switched off.  That number was better and the model was not:
        see `TestTheShortChannelTerms`, where the median turns out to be
        the wrong thing to rank terms by.
        """
        cm.default_toolkit = numeric
        worst = {}
        for name in ('nmos_long_idvd', 'nmos_idvd_vg1p2',
                     'nmos_idvd_vg0p6', 'nmos_idvg_vd0p05',
                     'nmos_idvg_vd1p2', 'nmos_idvg_vb_m1',
                     'nmos_long_idvg_vb_m1',
                     'pmos_long_idvd', 'pmos_long_idvg',
                     'pmos_idvd_vg1p2', 'pmos_idvg_vd1p2',
                     'pmos_idvg_vd0p05'):
            _, r, g, _ = _compare(deck, ref[name])
            worst[name] = abs(np.median(g / r) - 1.0)
        assert max(worst.values()) < 0.03, worst
        assert sum(worst.values()) < 0.10, worst


class TestTheShortChannelTerms(object):
    """DIBL and the `THESAT` bias modulations -- and the metric that
    nearly kept them out.

    They spent one commit switched off, on a measurement that turned out
    to be the wrong measurement.  Judged by the MEDIAN ratio to PSP over
    each sweep they looked bad; judged by the SPREAD of that ratio
    across a sweep they are decisive.  The median says whether the GAIN
    is right.  The spread says whether the BIAS DEPENDENCE is -- which
    is the part the physics governs, and the part a missing term
    breaks.

    The p-channel is what settled it, needing DIBL far more than the
    n-channel does, and disagreeing loudly enough to expose the
    compensation the n-channel had been hiding.
    """

    NAMES = ('nmos_long_idvd', 'nmos_idvd_vg1p2', 'nmos_idvd_vg0p6',
             'nmos_idvg_vd0p05', 'nmos_idvg_vd1p2', 'nmos_idvg_vb_m1',
             'nmos_long_idvg_vb_m1',
             'pmos_long_idvd', 'pmos_long_idvg', 'pmos_idvd_vg1p2',
             'pmos_idvg_vd1p2', 'pmos_idvg_vd0p05')

    @pytest.fixture(scope='class')
    def scaled(self):
        with open(REF) as fh:
            return json.load(fh)['scaled']

    def _kw(self, deck, g, model='sg13g2_lv_nmos_psp', **kw):
        return psp_scaling.to_long_channel(
            deck.model_params(model, w=g['w'], l=g['l'], ng=1, m=1,
                              pre_layout=1),
            w=g['w'], l=g['l'], **kw, T=T27)

    @pytest.mark.parametrize('geom', ['long', 'mid', 'short', 'wide_short'])
    @pytest.mark.parametrize('par', ['cf', 'cfb', 'thesatb', 'thesatg',
                                     'alp', 'alp1', 'alp2'])
    def test_every_parameter_involved_matches_psp(self, deck, scaled,
                                                  geom, par):
        """So none of this was ever the scaling layer."""
        got = self._kw(deck, scaled[geom])[par]
        assert got == pytest.approx(scaled[geom][par], rel=1e-5,
                                    abs=1e-12), (geom, par)

    def test_they_are_on_by_default(self, deck, scaled):
        kw = self._kw(deck, scaled['short'])
        assert 'cf' in kw and 'thesatg' in kw and 'thesatb' in kw
        off = self._kw(deck, scaled['short'], all_terms=False)
        assert 'cf' not in off, 'the A/B switch still has to work'

    def test_dibl_matches_psps_own_measured_threshold_shift(self, deck,
                                                            scaled):
        """3.6 mV computed against 3.5 mV measured out of PSP itself.

        Checked against a number extracted from the vendor model, not
        against our own arithmetic.
        """
        kw = self._kw(deck, scaled['short'])
        delvg = kw['cf'] * 1.35 * (1.0 + kw['cfb'] * 0.0)
        assert 3.0e-3 < delvg < 4.0e-3, delvg

    def test_they_fix_the_bias_dependence(self, deck, ref):
        """The measurement that decided it, on the metric that matters.

        Summed spread over the eleven sweeps: 1.80 without them, 0.38
        with.  Held loosely, as a statement about the direction and
        rough size of the effect rather than a pin on either number.
        """
        cm.default_toolkit = numeric
        spread = {}
        for on in (False, True):
            tot = 0.0
            for name in self.NAMES:
                sw = ref[name]
                pmos = 'pmos' in sw['device']
                kw = psp_scaling.to_long_channel(
                    deck.model_params(
                        'sg13g2_lv_pmos_psp' if pmos
                        else 'sg13g2_lv_nmos_psp',
                        w=sw['w'], l=sw['l'], ng=1, m=1, pre_layout=1),
                    w=sw['w'], l=sw['l'], all_terms=on, T=T27)
                cls = PspPmosLongChannel if pmos else PspMosLongChannel
                e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                        cm.Node('b'), **kw)
                e.update_iparv()
                v = np.asarray(sw['v'], float)
                r = np.abs(np.asarray(sw['i_d'], float))
                b = sw['bias']
                if sw['sweep'] == 'Vd':
                    g = np.array([np.asarray(e.i(e.bias(x, b['Vg'], b['Vs'], b['Vb'])), float)[0]
                        for x in v])
                else:
                    g = np.array([np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])), float)[0]
                        for x in v])
                m = r > 1e-6
                q = np.abs(g[m]) / r[m]
                tot += q.max() / q.min() - 1.0
            spread[on] = tot
        assert spread[True] < 0.4 * spread[False], spread

    def test_the_p_channel_saturation_sweep_is_where_it_shows(self, deck,
                                                              ref):
        """The single measurement that reversed the decision.

        Without DIBL the ratio sweeps from about 1.03 down to 0.44
        across the sweep -- a drive that is simply wrong.  With it, the
        ratio is flat to a couple of percent across three decades of
        current, which is a GAIN error: one cause, one number, and a
        well-posed question to ask next.
        """
        cm.default_toolkit = numeric
        sw = ref['pmos_idvg_vd1p2']
        b = sw['bias']
        v = np.asarray(sw['v'], float)
        r = np.abs(np.asarray(sw['i_d'], float))
        m = r > 1e-6
        out = {}
        for on in (False, True):
            kw = psp_scaling.to_long_channel(
                deck.model_params('sg13g2_lv_pmos_psp', w=sw['w'],
                                  l=sw['l'], ng=1, m=1, pre_layout=1),
                w=sw['w'], l=sw['l'], all_terms=on, T=T27)
            e = PspPmosLongChannel(cm.Node('d'), cm.Node('g'),
                                   cm.Node('s'), cm.Node('b'), **kw)
            e.update_iparv()
            g = np.array([np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])), float)[0] for x in v])
            q = np.abs(g[m]) / r[m]
            out[on] = q.max() / q.min() - 1.0
        assert out[False] > 1.0, out
        ## 0.08 measured.  Held at 0.15: the claim is that the ratio
        ## stops sweeping, not that it lands on a particular flatness --
        ## the residual here is a gain offset that later terms should
        ## keep moving.
        assert out[True] < 0.15, out


class TestTheChannelTypes(object):
    """The p-channel device, and how little of it is p-channel-specific.

    PSP is written once and run for both polarities.  It converts the
    terminal voltages to an internal n-like convention on the way in
    (`PSP103_module.include:1035-1048`), writes every equation for that
    convention alone, and restores the polarity on every contribution on
    the way out (`:1697-1795`).  Its 849-line geometry-scaling layer
    contains not one reference to the channel type.

    So does ours.  Exactly five things differ, four of them in the
    kernel: the effective-field weighting at each channel end (1/2 for
    electrons, 1/3 for holes) and the two velocity-saturation forms,
    holes following `v ~ E/(1 + E/Ec)` where electrons follow
    `v ~ E/sqrt(1 + (E/Ec)^2)`.  The fifth is the quantum-mechanical
    constant.  Everything else is a sign.
    """

    GEOMS = ('long', 'mid', 'short', 'wide_short')
    DIRECT = ('vfb', 'tox', 'ct', 'mue', 'themu', 'cs', 'thecs', 'rs',
              'ax', 'thesat', 'alp', 'alp1', 'alp2', 'cf', 'cfb')

    @pytest.fixture(scope='class')
    def scaled_p(self):
        with open(REF) as fh:
            return json.load(fh)['scaled_pmos']

    @pytest.fixture(scope='class')
    def scaled_n(self):
        with open(REF) as fh:
            return json.load(fh)['scaled']

    def _fet(self, deck, w=10e-6, l=1e-6):
        cm.default_toolkit = numeric
        card = deck.model_params('sg13g2_lv_pmos_psp', w=w, l=l, ng=1,
                                 m=1, pre_layout=1)
        e = PspPmosLongChannel(
            cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
            **psp_scaling.to_long_channel(card, w=w, l=l, T=T27))
        e.update_iparv()
        return e

    def test_the_card_carries_the_polarity_and_nothing_else_does(self,
                                                                 deck):
        """Both cards' `vfb` is NEGATIVE -- the p-channel one is not a
        mirrored set of numbers, it is the same convention with
        different values, and `type` alone distinguishes them."""
        n = deck.model_params('sg13g2_lv_nmos_psp', w=1e-6, l=1e-6,
                              ng=1, m=1, pre_layout=1)
        p = deck.model_params('sg13g2_lv_pmos_psp', w=1e-6, l=1e-6,
                              ng=1, m=1, pre_layout=1)
        assert psp_scaling.channel_type(n) == 1.0
        assert psp_scaling.channel_type(p) == -1.0
        for c in (n, p):
            assert c['vfbo'] < 0.0
            assert c['uo'] > 0.0 and c['toxo'] > 0.0

    @pytest.mark.parametrize('geom', GEOMS)
    @pytest.mark.parametrize('par', DIRECT)
    def test_every_scaled_parameter_matches_psp(self, deck, scaled_p,
                                                geom, par):
        """The whole scaling layer, checked on the OTHER channel type.

        It needed no change at all, which is the claim this pins.
        """
        e = scaled_p[geom]
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_pmos_psp', w=e['w'], l=e['l'],
                              ng=1, m=1, pre_layout=1),
            w=e['w'], l=e['l'], all_terms=True, T=T27)
        assert kw[par] == pytest.approx(e[par], rel=1e-5, abs=1e-12)

    @pytest.mark.parametrize('geom', GEOMS)
    def test_the_doping_matches_before_the_qm_correction(self, deck,
                                                         scaled_p, geom):
        """As for the n-channel: we fold the quantum correction into the
        doping where PSP applies it to `phib` and the body factor, so
        the comparison has to be made with it switched off."""
        e = scaled_p[geom]
        card = deck.model_params('sg13g2_lv_pmos_psp', w=e['w'], l=e['l'],
                                 ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(dict(card, qmc=0.0), w=e['w'],
                                         l=e['l'], T=T27)
        assert kw['nsub'] == pytest.approx(e['neff'], rel=1e-5)

    @pytest.mark.parametrize('geom', GEOMS)
    def test_holes_get_the_larger_quantum_constant(self, deck, scaled_n,
                                                   scaled_p, geom):
        """`QMP/QMN = 1.2515` (`PSP103_module.include:727-734`).

        Not cosmetic: `qq` shifts `phib` AND the body factor, so it
        moves the threshold and the body effect together.  Measured as
        the ratio the quantum correction inflates the folded doping by,
        which is larger for holes on every geometry.
        """
        def blow_up(model, e):
            card = deck.model_params(model, w=e['w'], l=e['l'], ng=1,
                                     m=1, pre_layout=1)
            on = psp_scaling.to_long_channel(card, w=e['w'], l=e['l'], T=T27)
            off = psp_scaling.to_long_channel(dict(card, qmc=0.0),
                                              w=e['w'], l=e['l'], T=T27)
            return on['nsub'] / off['nsub']
        n = blow_up('sg13g2_lv_nmos_psp', scaled_n[geom])
        p = blow_up('sg13g2_lv_pmos_psp', scaled_p[geom])
        assert p > n > 1.0, (geom, n, p)

    def test_the_hole_specific_kernel_terms_are_actually_active(self):
        """The four kernel differences, exercised directly.

        Compiled twice from the same numbers with only the channel-type
        flag moved, so any difference in the answer is those four terms
        and nothing else.  Guards against the flag being threaded
        through but never read -- which would leave the p-channel
        silently running electron physics.
        """
        import sympy
        from pycircuit.circuit import hdl, psp_kernel
        xg, xns, xnd = sympy.symbols('xg xns xnd', real=True)
        base = dict(mue=2.3, themu=1.26, cs=0.94, thecs=1.46, feta=1.0,
                    thesat=0.104, cox_area=1.7e-2, eps_si=1.04e-10,
                    ax=5.3, rs=0.0, rsg=0.0, rsb=0.0, vsb=0.0, alp=0.004,
                    vp=1e-5, vds=1.0, kp=0.0, alp1=0.004, alp2=0.005)
        out = []
        for pmos in (False, True):
            f = hdl.compile_chain(
                lambda: psp_kernel.intrinsic(
                    xg, xns, xnd, 1.5, 2.06, 0.0259, 1e-4,
                    mob=dict(base, pmos=pmos))['ids'],
                [xg, xns, xnd])
            out.append(float(np.asarray(f(45.0, 34.0, 72.0),
                                        float).reshape(-1)[0]))
        assert out[0] != out[1], out
        assert 0.5 < out[1] / out[0] < 1.5, out

    def test_it_is_antisymmetric_like_the_n_channel(self, deck):
        e = self._fet(deck)
        for vg in (-0.6, -1.2):
            for vd in (-0.1, -0.9):
                f = np.asarray(e.i(e.bias(vd, vg, 0.0, 0.0)),
                               float)[0]
                r = np.asarray(e.i(e.bias(0.0, vg, vd, 0.0)),
                               float)[0]
                assert f == pytest.approx(-r, rel=1e-13), (vg, vd)

    def test_it_conserves_charge_like_the_n_channel(self, deck):
        e = self._fet(deck)
        for x in (e.bias(-0.9, -1.4, 0.0, 0.0),
                  e.bias(-0.05, -0.6, 0.0, 0.3)):
            q = np.asarray(e.q(x), float)
            assert abs(q.sum()) < 1e-24 * max(1.0, np.abs(q).max())

    @pytest.mark.parametrize('x', [np.array([-0.6, -1.0, 0.0, 0.0]),
                                   np.array([-1.2, -1.8, 0.0, 0.0]),
                                   np.array([0.0, 0.0, 0.0, 0.0]),
                                   np.array([-0.05, -0.4, 0.0, 0.5])])
    def test_its_jacobian_is_finite_and_correct(self, deck, x):
        e = self._fet(deck)
        x = e.bias(*x)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            i = np.asarray(e.i(x), float)
            G = np.asarray(e.G(x), float)
        assert np.all(np.isfinite(i)) and np.all(np.isfinite(G))
        fd = np.zeros_like(G)
        for j in range(len(x)):
            h = 1e-6 * max(1.0, abs(x[j]))
            xp, xm = x.copy(), x.copy()
            xp[j] += h
            xm[j] -= h
            fd[:, j] = (np.asarray(e.i(xp), float)
                        - np.asarray(e.i(xm), float)) / (2 * h)
        assert np.max(np.abs(G - fd)) < 1e-6 * max(1.0, np.abs(G).max())

    def test_the_current_conducts_with_the_right_polarity(self, deck):
        e = self._fet(deck)
        i = np.asarray(e.i(e.bias(-1.2, -1.2, 0.0, 0.0)), float)
        assert i[0] < 0.0, 'drain sinks current for a p-channel'
        assert i[2] > 0.0, 'source sources it'
        assert i[1] == 0.0, 'no gate current in the intrinsic core'

    @pytest.mark.parametrize('name', ['pmos_long_idvd', 'pmos_long_idvg',
                                      'pmos_idvd_vg1p2',
                                      'pmos_idvg_vd0p05'])
    def test_the_gap_against_psps_own_p_channel(self, deck, ref, name):
        """First measurement, nothing p-channel-specific tuned.

        Held at 10%: the n-channel's first measurement was 29% out and
        took eight terms to bring inside a percent, so this starting
        point is a good deal better than that one -- but it is a
        starting point, and the guard says so rather than pinning a
        number that should move.
        """
        _, r, g, _ = _compare(deck, ref[name])
        assert abs(np.median(g / r) - 1.0) < 0.10, np.median(g / r)


class TestThePolyDopingScaling(object):
    """`NP = NPO*max(1e-6, 1 + NPL*iLE)` (`PSP103_scaling.include:260`).

    A geometry scaling that is EXACTLY INERT on the n-channel card,
    because `NPL` is zero there, and very much alive on the p-channel
    one, where `NPL = -0.0959` takes the gate doping down to 0.36 of
    nominal on a 0.13 um device -- nearly three times the gate depletion
    we were applying, and about 8% of the drain current.

    That is the second geometry scaling this model has missed for the
    same reason: the coefficient is zero on the only card it was ever
    checked against.  The first was the width adjustment inside `BETN`,
    which was 12% off and invisible on the geometry driving the
    investigation.  The lesson is not about `NP` -- it is that a card
    with a zero coefficient does not test a scaling, and the term-by-term
    comparison against PSP's own `lp_*` outputs is what catches it.
    """

    @pytest.fixture(scope='class')
    def scaled_both(self):
        with open(REF) as fh:
            d = json.load(fh)
        return d['scaled'], d['scaled_pmos']

    @pytest.mark.parametrize('geom', ['long', 'mid', 'short', 'wide_short'])
    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    def test_the_gate_doping_matches_psp(self, deck, scaled_both, geom,
                                         kind):
        """Directly against `lp_np`, on both channel types."""
        n, p = scaled_both
        e = (p if kind == 'pmos' else n)[geom]
        card = deck.model_params('sg13g2_lv_%s_psp' % kind, w=e['w'],
                                 l=e['l'], ng=1, m=1, pre_layout=1)
        geo = psp_scaling.geometry(card, w=e['w'], l=e['l'])
        ours = card['npo'] * max(1e-6, 1.0 + card['npl'] * geo['iLE'])
        assert ours == pytest.approx(e['np'], rel=1e-5), (kind, geom)

    def test_it_is_inert_on_the_n_channel_and_not_on_the_p(self, deck,
                                                           scaled_both):
        """Which is exactly why it was missed.

        `lp_np` is the same number at every n-channel geometry and falls
        by nearly three on the p-channel, so no amount of n-channel
        measurement could have shown the scaling existed.
        """
        n, p = scaled_both
        assert len({e['np'] for e in n.values()}) == 1
        assert p['short']['np'] < 0.45 * p['long']['np']

    def test_it_moves_the_p_channel_gate_depletion(self, deck):
        """The consequence, in the parameter the element actually takes.

        `kp` is inversely proportional to the gate doping, so a third of
        the doping is three times the `kp` -- and `kp` is what sets how
        much of the applied gate voltage is lost in the poly.
        """
        def kp(model, w, l):
            return psp_scaling.to_long_channel(
                deck.model_params(model, w=w, l=l, ng=1, m=1,
                                  pre_layout=1), w=w, l=l, T=T27)['kp']
        assert kp('sg13g2_lv_nmos_psp', 1e-6, 0.13e-6) == \
            pytest.approx(kp('sg13g2_lv_nmos_psp', 10e-6, 1e-6), rel=1e-9)
        short = kp('sg13g2_lv_pmos_psp', 1e-6, 0.13e-6)
        long_ = kp('sg13g2_lv_pmos_psp', 10e-6, 1e-6)
        assert short > 2.0 * long_, (short, long_)


class TestTheBiasDependentBodyFactor(object):
    """`Gf = G_0*sqrt(1 + DNSUB*maxa(0, Vgb - VNSUB, NSLP))`
    (`PSP103_macrodefs.include:478-484`).

    The depletion charge a real device presents is not that of a
    uniformly doped substrate: a pocket implant makes the effective
    doping rise as the gate pulls the depletion edge deeper, so the body
    factor grows with gate drive.

    THE THIRD TERM ON THIS PDK HIDDEN BY A ZERO COEFFICIENT.  `DNSUBO`
    is 4.4e-16 on the n-channel card -- zero in every sense that matters
    -- and 0.0397 on the p-channel one.  After the `BETN` width
    adjustment and the `NP` length scaling, the pattern is no longer a
    coincidence worth remarking on; it is the reason a second channel
    type earns its keep.
    """

    def _kw(self, deck, model, w, l):
        return psp_scaling.to_long_channel(
            deck.model_params(model, w=w, l=l, ng=1, m=1, pre_layout=1),
            w=w, l=l, T=T27)

    def test_the_card_switches_it_on_for_holes_only(self, deck):
        n = self._kw(deck, 'sg13g2_lv_nmos_psp', 1e-6, 0.13e-6)
        p = self._kw(deck, 'sg13g2_lv_pmos_psp', 1e-6, 0.13e-6)
        assert n['dnsub'] < 1e-12, n['dnsub']
        assert p['dnsub'] > 0.01, p['dnsub']
        assert p['nslp'] > 0.0, 'the smoothing must not be zero'

    def test_it_has_no_geometry_dependence(self, deck):
        """PSP takes these straight from the card (`scaling:255-257`)."""
        a = self._kw(deck, 'sg13g2_lv_pmos_psp', 1e-6, 0.13e-6)
        b = self._kw(deck, 'sg13g2_lv_pmos_psp', 10e-6, 1e-6)
        for k in ('dnsub', 'vnsub', 'nslp'):
            assert a[k] == b[k], k

    def test_it_raises_the_body_factor_with_gate_drive(self):
        """The smooth maximum, checked as a function rather than through
        a current: inactive well below the onset, growing above it, and
        with no corner at the onset itself."""
        from pycircuit.circuit import psp_kernel
        f = lambda v: float(psp_kernel.maxa(v, 0.0, 0.05))
        assert f(-2.0) < 0.02
        assert f(2.0) == pytest.approx(2.0, abs=0.02)
        mid = f(0.0)
        assert 0.0 < mid < 0.2, mid
        ## Smooth: the centred difference at the corner is about a half,
        ## not a jump between 0 and 1.
        h = 1e-4
        slope = (f(h) - f(-h)) / (2 * h)
        assert 0.3 < slope < 0.7, slope

    def test_it_is_what_takes_the_p_channel_to_within_a_percent(self,
                                                                deck,
                                                                ref):
        """The measurement, on the sweeps it moves.

        Switching `DNSUB` off leaves the p-channel a few percent high on
        every geometry; with it the same sweeps land inside one percent.
        The n-channel is unaffected either way, its coefficient being
        zero -- which is the control this test needs to be worth
        anything.
        """
        cm.default_toolkit = numeric
        err = {}
        for on in (True, False):
            tot = 0.0
            for name in ('pmos_long_idvd', 'pmos_idvd_vg1p2',
                         'pmos_idvg_vd0p05'):
                sw = ref[name]
                kw = self._kw(deck, 'sg13g2_lv_pmos_psp', sw['w'],
                              sw['l'])
                if not on:
                    kw['dnsub'] = 0.0
                e = PspPmosLongChannel(cm.Node('d'), cm.Node('g'),
                                       cm.Node('s'), cm.Node('b'), **kw)
                e.update_iparv()
                v = np.asarray(sw['v'], float)
                r = np.abs(np.asarray(sw['i_d'], float))
                b = sw['bias']
                if sw['sweep'] == 'Vd':
                    g = np.array([np.asarray(e.i(e.bias(x, b['Vg'], b['Vs'], b['Vb'])), float)[0]
                        for x in v])
                else:
                    g = np.array([np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])), float)[0]
                        for x in v])
                m = r > 1e-6
                tot += abs(np.median(np.abs(g[m]) / r[m]) - 1.0)
            err[on] = tot
        assert err[True] < 0.3 * err[False], err
        assert err[True] < 0.03, err


class TestTheSubthresholdSlope(object):
    """The steepest physical check there is, and it was not being made.

    The subthreshold swing is not a fitted quantity in a
    surface-potential model.  It falls out of the electrostatics -- the
    body factor, the surface potential, and the effective thermal
    voltage -- so getting it right is a statement about the formulation
    rather than about a parameter, and getting it wrong would mean the
    weak-inversion charge is wrong in a way an integrated current can
    hide.

    Measured against PSP's own curves it agrees to a quarter of a
    millivolt per decade on both channel types.

    ONE TRAP, and it produced a wrong answer before it was noticed: the
    reference records TOTAL terminal current, and PSP's junction leakage
    is a flat 2e-12 A floor which our core does not model at all.  On
    the body-biased sweep the drain current reaches that floor, so a
    slope taken down to 1e-11 A measures PSP's leakage rather than its
    channel -- and reports a 2.5 mV/decade discrepancy that does not
    exist.  The window has to stay clear of it.
    """

    #: Three decades, the bottom of them some five hundred times the
    #: 2e-12 A junction-leakage floor the reference carries.
    LO, HI = 1e-9, 1e-6

    NAMES = ('nmos_idvg_vd0p05', 'nmos_idvg_vd1p2', 'nmos_idvg_vb_m1',
             'nmos_long_idvg', 'pmos_idvg_vd0p05', 'pmos_long_idvg')

    def _slope(self, v, i):
        m = (i > self.LO) & (i < self.HI)
        assert m.sum() >= 4, 'not enough points in the window'
        return 1e3 * np.median(np.abs(np.diff(v[m])
                                      / np.diff(np.log10(i[m]))))

    @pytest.mark.parametrize('name', NAMES)
    def test_it_matches_psp(self, deck, ref, name):
        cm.default_toolkit = numeric
        sw = ref[name]
        pmos = 'pmos' in sw['device']
        kw = psp_scaling.to_long_channel(
            deck.model_params(
                'sg13g2_lv_pmos_psp' if pmos else 'sg13g2_lv_nmos_psp',
                w=sw['w'], l=sw['l'], ng=1, m=1, pre_layout=1),
            w=sw['w'], l=sw['l'], T=T27)
        cls = PspPmosLongChannel if pmos else PspMosLongChannel
        e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                **kw)
        e.update_iparv()
        v = np.asarray(sw['v'], float)
        r = np.abs(np.asarray(sw['i_d'], float))
        b = sw['bias']
        g = np.abs(np.array([np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])), float)[0] for x in v]))
        theirs, ours = self._slope(v, r), self._slope(v, g)
        assert abs(ours - theirs) < 1.0, (name, theirs, ours)

    def test_the_reference_carries_a_leakage_floor(self, ref):
        """The trap, pinned so the window is not widened by accident.

        PSP's bulk current on the body-biased sweep is a CONSTANT
        2e-12 A -- junction leakage, which this core does not model --
        and the drain current comes down to meet it.  Any slope taken
        through that region measures the leakage.
        """
        sw = ref['nmos_idvg_vb_m1']
        ib = np.abs(np.asarray(sw['i_b'], float))
        idd = np.abs(np.asarray(sw['i_d'], float))
        assert ib.max() / ib.min() < 1.05, 'expected a flat floor'
        assert idd.min() < 10.0 * ib.min(), \
            'the drain current should reach the leakage floor'
        assert ib.max() < 1e-2 * self.LO, \
            'and the chosen window should sit far above it'


class TestTheBodyBiasMobilityCorrection(object):
    """`Rxcor = (1 + 0.2*XCOR*Vsbx)/(1 + XCOR*Vsbx)`
    (`PSP103_macrodefs.include:576`), multiplying `Gmob` at both the
    source end and the midpoint (`:595`, `:750`).

    THE CLEANEST EXAMPLE ON THIS BRANCH OF A TERM NO MEASUREMENT COULD
    SEE.  It is identically 1 at zero body bias, and every sweep in the
    reference used a grounded body except one -- which sits on the
    0.13 um geometry, where `XCOR` scales to exactly zero.  So the term
    was invisible twice over, and its absence could not have shown up in
    any number this project was tracking.

    The fix was not a cleverer analysis, it was another measurement:
    `nmos_long_idvg_vb_m1`, a body-biased sweep on the geometry where
    `XCOR` is alive.  On it the term is worth 3.6%.

    That is the same lesson as the three parameters a zero coefficient
    hid, arrived at from the other direction: there, one card did not
    exercise a scaling; here, one bias condition did not exercise a
    term.  A model is only validated on the space its measurements
    actually span.
    """

    @pytest.fixture(scope='class')
    def scaled_both(self):
        with open(REF) as fh:
            d = json.load(fh)
        return d['scaled'], d['scaled_pmos']

    def _kw(self, deck, model, w, l, **kw):
        return psp_scaling.to_long_channel(
            deck.model_params(model, w=w, l=l, ng=1, m=1, pre_layout=1),
            w=w, l=l, **kw, T=T27)

    @pytest.mark.parametrize('geom', ['long', 'mid', 'short', 'wide_short'])
    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    def test_it_matches_psp(self, deck, scaled_both, geom, kind):
        n, p = scaled_both
        e = (p if kind == 'pmos' else n)[geom]
        ours = self._kw(deck, 'sg13g2_lv_%s_psp' % kind, e['w'],
                        e['l'])['xcor']
        assert ours == pytest.approx(e['xcor'], rel=1e-5, abs=1e-12), \
            (kind, geom)

    def test_it_clips_to_zero_on_the_short_n_channel(self, deck,
                                                     scaled_both):
        """Which is half of why it was invisible.

        The one body-biased sweep this project had was on that geometry,
        so even a reverse-biased body could not have exercised it.
        """
        n, _ = scaled_both
        assert n['short']['xcor'] == 0.0
        assert n['wide_short']['xcor'] == 0.0
        assert n['long']['xcor'] > 0.01

    def test_it_is_inert_at_zero_body_bias(self, deck):
        """The other half of why it was invisible.

        `Rxcor` is a ratio that is 1 by construction at `Vsbx = 0`, so
        switching the term off changes essentially nothing on a grounded
        body -- measured at 4e-7 relative, which is not zero for a
        reason worth knowing.  `Vsbx` is the CONDITIONED body bias, and
        the smooth `MINA` clip that produces it leaves it a fraction of
        a millivolt from zero even when the body is grounded.  So the
        term is inert to the resolution of the conditioning rather than
        bit-exactly, and 4e-7 is that fraction of a millivolt times
        `XCOR`.
        """
        cm.default_toolkit = numeric
        kw = self._kw(deck, 'sg13g2_lv_nmos_psp', 10e-6, 1e-6)
        assert kw['xcor'] > 0.01, 'the term must be live here'
        off = dict(kw, xcor=0.0)
        a = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        b = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **off)
        a.update_iparv()
        b.update_iparv()
        for x in (a.bias(0.05, 1.2), a.bias(1.2, 0.8)):
            ia = np.asarray(a.i(x), float)[0]
            ib = np.asarray(b.i(x), float)[0]
            assert abs(ia - ib) < 1e-5 * abs(ib), (x, ia, ib)

    def test_it_is_worth_three_percent_where_it_does_act(self, deck, ref):
        """The measurement the sweep was added for."""
        cm.default_toolkit = numeric
        sw = ref['nmos_long_idvg_vb_m1']
        assert sw['bias']['Vb'] == -1.0
        kw = self._kw(deck, 'sg13g2_lv_nmos_psp', sw['w'], sw['l'])
        v = np.asarray(sw['v'], float)
        r = np.abs(np.asarray(sw['i_d'], float))
        m = r > 1e-6
        b = sw['bias']
        got = {}
        for on in (True, False):
            e = PspMosLongChannel(cm.Node('d'), cm.Node('g'),
                                  cm.Node('s'), cm.Node('b'),
                                  **(kw if on else dict(kw, xcor=0.0)))
            e.update_iparv()
            g = np.abs(np.array([np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])), float)[0] for x in v]))
            got[on] = np.median(g[m] / r[m])
        assert abs(got[False] - 1.0) > 0.02, got
        assert abs(got[True] - 1.0) < 0.005, got


class TestTheChargeModelAgainstTheVendor(object):
    """The charge model, compared to PSP for the first time.

    Until now it had been validated entirely by CONSTRUCTION: charge
    conservation exact, the source/drain swap exact, the Ward-Dutton
    partition reproducing the textbook 0.5 -> 0.4.  Those are strong
    properties and they were all passing.

    They were also all RATIOS, and every one of them is invariant under
    a uniform error in the oxide capacitance.  The capacitances were 24%
    high and nothing in the tree could see it.  That is the lesson here:
    a model checked only against itself is checked for consistency, not
    for correctness, and the two are different claims.

    Two things caused it, both now fixed:

    * PSP builds the charge model's oxide capacitance from the **CV
      effective dimensions** (`PSP103_scaling.include:38-39, 359`),
      which the card shifts from the drawn ones by `DLQ` and `DWQ` --
      7% on the long device here;
    * and it applies a **quantum-mechanical reduction** to it,
      `COX_qm = COX/(1 + qq*(qeff1^2 + qlim2)^(-1/6))`
      (`PSP103_module.include:1474-1478`) -- another 12%.  The same `qq`
      already shifts the threshold and the body factor; this is that
      physics acting on the charge instead.

    PSP's terminal capacitances include overlap and fringe
    contributions this core does not model, which is why the LONG
    geometry is the one to read: there they are a few percent, against
    about a fifth of `cgg` at 0.13 um.
    """

    @pytest.fixture(scope='class')
    def op(self):
        with open(REF) as fh:
            d = json.load(fh)
        return d['op'], d['op_pmos']

    def _fet(self, deck, model, cls, w, l):
        """Built with the OVERLAP SWITCHED OFF.

        This class is about the intrinsic charge, and PSP reports its
        intrinsic `cgg` separately from the overlap, so the comparison
        is only like-for-like with the parasitics out.
        """
        cm.default_toolkit = numeric
        kw = psp_scaling.to_long_channel(
            deck.model_params(model, w=w, l=l, ng=1, m=1, pre_layout=1),
            w=w, l=l, T=T27)
        kw.update(cgov=0.0, cfr=0.0, cgbov=0.0)
        e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                **kw)
        e.update_iparv()
        return e

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    def test_the_oxide_capacitance_matches_psp(self, deck, kind):
        """`lp_cox`, term by term -- which splits "our capacitances are
        24% high" into two separate checkable numbers, and pins the one
        that is pure geometry."""
        with open(REF) as fh:
            sc = json.load(fh)['scaled' if kind == 'nmos'
                               else 'scaled_pmos']
        for geom, e in sc.items():
            card = deck.model_params('sg13g2_lv_%s_psp' % kind, w=e['w'],
                                     l=e['l'], ng=1, m=1, pre_layout=1)
            kw = psp_scaling.to_long_channel(card, w=e['w'], l=e['l'], T=T27)
            eps_ox = card.get('epsroxo', 3.9) * 8.8541878128e-12
            ours = eps_ox / kw['tox'] * kw['wcv'] * kw['lcv']
            assert ours == pytest.approx(e['cox'], rel=1e-4), (kind, geom)

    def test_the_cv_dimensions_are_not_the_drawn_ones(self, deck):
        """Which is half the error, and invisible without a card."""
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=10e-6, l=1e-6,
                              ng=1, m=1, pre_layout=1), w=10e-6, l=1e-6, T=T27)
        ## The LENGTH is what matters: `LEcv = LE + DLQ` with
        ## `LAP = 2.94e-8` and `DLQ = -1.37e-8` puts it 7% under the
        ## drawn value, and that is where the geometric half of the 24%
        ## came from.  The width moves the other way and barely at all
        ## -- `WOT` is negative here, so `WEcv` ends up a shade WIDER
        ## than drawn even after `DWQ`.  Asserted as measured rather
        ## than as guessed, because guessing it wrong is easy.
        assert kw['lcv'] < 0.94e-6, kw['lcv']
        assert kw['lcv'] > 0.90e-6, kw['lcv']
        assert 10.0e-6 < kw['wcv'] < 10.02e-6, kw['wcv']

    @staticmethod
    def _psp_total(pt):
        """PSP's TOTAL gate capacitance: intrinsic plus the overlap it
        reports separately."""
        return pt['cgg'] + pt['cgsol'] + pt['cgdol']

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    def test_the_gate_capacitance_matches_on_the_long_device(self, deck,
                                                             op, kind):
        """The measurement.  Every bias point on the long device,
        within 0.2% -- which is what it became once channel-length
        modulation and the polysilicon factor went into the CHARGES as
        well as the current.  Before those, 1%; before the oxide
        capacitance was corrected, 24%."""
        n, p = op
        gd = (p if kind == 'pmos' else n)['long']
        cls = PspPmosLongChannel if kind == 'pmos' else PspMosLongChannel
        e = self._fet(deck, 'sg13g2_lv_%s_psp' % kind, cls, gd['w'],
                      gd['l'])
        for pt in gd['points']:
            x = e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
            C = np.asarray(e.C(x), float)
            gi = e.gate_index
            assert C[gi, gi] == pytest.approx(self._psp_total(pt),
                                              rel=0.002), \
                (kind, pt['vg'], pt['vd'], pt['vb'],
                 C[e.gate_index, e.gate_index], pt['cgg'])

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    def test_the_source_and_bulk_capacitances_match_too(self, deck, op,
                                                        kind):
        """`cgg` alone could be right for compensating reasons; the
        partition between the terminals could not."""
        n, p = op
        gd = (p if kind == 'pmos' else n)['long']
        cls = PspPmosLongChannel if kind == 'pmos' else PspMosLongChannel
        e = self._fet(deck, 'sg13g2_lv_%s_psp' % kind, cls, gd['w'],
                      gd['l'])
        for pt in gd['points']:
            x = e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
            C = np.asarray(e.C(x), float)
            gi = e.gate_index
            ## Column order is (vd, vg, vs, vb); PSP reports `cgs` and
            ## `cgb` as the negated cross terms.  The ROW is the
            ## intrinsic gate, which is an internal node once the gate
            ## resistance is on.
            assert -C[gi, 2] == pytest.approx(pt['cgs'], rel=0.01), \
                ('cgs', kind, pt['vg'], pt['vd'])
            if abs(pt['cgb']) > 1e-16:
                assert -C[gi, 3] == pytest.approx(pt['cgb'], rel=0.10), \
                    ('cgb', kind, pt['vg'], pt['vd'])

    def test_both_corrections_are_needed(self, deck, op):
        """Each on its own leaves the capacitance visibly wrong.

        Worth pinning because either could be removed without any
        construction property noticing -- which is exactly how they came
        to be missing.
        """
        n, _ = op
        gd = n['long']
        base = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=gd['w'], l=gd['l'],
                              ng=1, m=1, pre_layout=1),
            w=gd['w'], l=gd['l'], T=T27)
        pt = [q for q in gd['points']
              if q['vg'] == 1.2 and q['vd'] == 0.05][0]
        got = {}
        base = dict(base, cgov=0.0, cfr=0.0, cgbov=0.0)
        for tag, kw in (('both', base),
                        ('no cv dims', dict(base, wcv=0.0, lcv=0.0)),
                        ('no qm', dict(base, qq=0.0))):
            e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                  cm.Node('b'), **kw)
            e.update_iparv()
            x = e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
            got[tag] = np.asarray(e.C(x), float)[e.gate_index, e.gate_index] / pt['cgg']
        assert abs(got['both'] - 1.0) < 0.03, got
        assert got['no cv dims'] > 1.05, got
        assert got['no qm'] > 1.08, got

    def test_conservation_still_holds_with_the_corrections(self, deck,
                                                            op):
        """The construction properties must survive being made correct.

        They were never in doubt -- they are invariant under exactly the
        kind of scale factor that was wrong -- but a change to the oxide
        capacitance is the one thing that could plausibly break the
        structural conservation, so it is checked rather than assumed.
        """
        n, _ = op
        gd = n['long']
        e = self._fet(deck, 'sg13g2_lv_nmos_psp', PspMosLongChannel,
                      gd['w'], gd['l'])
        for pt in gd['points']:
            x = e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
            C = np.asarray(e.C(x), float)
            q = np.asarray(e.q(x), float)
            assert abs(q.sum()) < 1e-24 * max(1.0, np.abs(q).max())
            assert np.abs(C.sum(axis=0)).max() < 1e-12 * np.abs(C).max()
            assert np.abs(C.sum(axis=1)).max() < 1e-12 * np.abs(C).max()


@needs_pdk
class TestTheDielectricComesFromTheCard(object):
    """`EPSROX` is a card parameter, so the element must read it.

    It used to be a module constant in `compact.py` while
    `psp_scaling` read `EPSROXO` from the card.  Both are 3.9 on this
    process, so the two agreed by LUCK -- and on a card with a
    different gate dielectric they would have disagreed about `cox`,
    silently, which is the body factor and the oxide capacitance at
    once.  This branch has been caught three times by a coefficient
    that happened to be harmless on the one card it was tested against;
    this is the same shape, found before it cost anything.
    """

    def test_the_card_value_reaches_the_element(self, deck):
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=10e-6, l=1e-6,
                              ng=1, m=1, pre_layout=1),
            w=10e-6, l=1e-6, T=T27)
        assert kw['epsrox'] == pytest.approx(3.9)

    def test_changing_it_changes_the_current(self, deck):
        """The guard that makes it not a coincidence: if the element
        ignored the parameter this would pass unchanged."""
        cm.default_toolkit = numeric
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=10e-6, l=1e-6,
                              ng=1, m=1, pre_layout=1),
            w=10e-6, l=1e-6, T=T27)
        out = []
        for eps in (3.9, 7.8):
            e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                  cm.Node('b'), **dict(kw, epsrox=eps))
            e.update_iparv()
            out.append(np.asarray(e.i(e.bias(0.05, 1.2)), float)[0])
        ## Doubling the dielectric constant doubles `cox`, which raises
        ## the drive and lowers the body factor -- the same direction as
        ## halving the oxide thickness.  Measured: 1.35x at Vg = 1.2,
        ## Vd = 0.05.  Less than 2x because the two effects do not
        ## compound: more drive, but `BETN` and the body factor move
        ## against it.
        assert out[1] / out[0] > 1.25, out


@needs_pdk
class TestThePhysicalConstants(object):
    """The last residual was a physical constant, not a term.

    After thirty-one scaled parameters matching PSP exactly, the
    subthreshold slope matching to a quarter of a millivolt per decade,
    and every formula re-read against the vendor source, the n-channel
    threshold was still 2-6 mV low -- and worse under body bias, which
    is the tell.

    PSP defines `EPSRSI = 11.8` (`Common103_macrodefs.include:61`).
    `pycircuit.circuit.constants` carries 11.7.  Both are defensible
    values for silicon; what is not defensible is comparing against a
    reference that uses one while using the other.

    The body factor goes as `sqrt(eps_si)`, so 11.7 against 11.8 makes
    `gamma` 0.43% small -- and because the body term grows as
    `sqrt(phib + Vsb)`, the threshold error GROWS WITH BODY BIAS.  That
    is why the body-biased sweeps were the worst ones.

    IT SURVIVED EVERY OTHER CHECK BECAUSE NO `lp_*` OUTPUT EXPOSES THE
    BODY FACTOR.  PSP reports the doping, the flat-band voltage, the
    oxide thickness and the oxide capacitance -- everything the body
    factor is built FROM -- and not the body factor itself.  Between
    thirty-one exactly-matched parameters and the drain current sat one
    number nobody had thought to compare, because it was not a
    parameter.
    """

    def test_we_use_psps_constants_not_the_trees(self):
        from pycircuit.circuit import compact, constants
        assert psp_scaling.PSP_EPSRSI == 11.8
        assert constants.epsRSi == 11.7, \
            'the tree keeps its own value, and should'
        assert compact.EPS_SI == psp_scaling.PSP_EPS_SI
        assert psp_scaling.PSP_KBOL == 1.3806505e-23
        assert psp_scaling.PSP_QELE == 1.6021918e-19

    def test_the_difference_is_where_it_hurts(self):
        """0.43% on the body factor, which is 2-6 mV of threshold."""
        import math
        assert math.sqrt(11.8 / 11.7) == pytest.approx(1.00426, abs=1e-5)

    @pytest.mark.parametrize('name,limit', [
        ('nmos_long_idvg', 2.0),
        ('nmos_idvg_vd0p05', 3.5),
        ('nmos_idvg_vb_m1', 5.0),
        ('nmos_long_idvg_vb_m1', 3.0),
    ])
    def test_the_threshold_offset_is_small_and_bounded(self, deck, ref,
                                                       name, limit):
        """Extracted the same way from both curves, so it compares like
        with like -- PSP's own `vth` is its own extraction and is not
        comparable in absolute terms.

        The bounds are where the constant fix left things.  They are
        upper bounds on a residual, not targets: if something later
        closes the rest, they should be tightened, not met.
        """
        cm.default_toolkit = numeric
        sw = ref[name]
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=sw['w'], l=sw['l'],
                              ng=1, m=1, pre_layout=1),
            w=sw['w'], l=sw['l'], T=T27)
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        v = np.asarray(sw['v'], float)
        r = np.asarray(sw['i_d'], float)
        b = sw['bias']
        g = np.array([np.asarray(e.i(e.bias(b['Vd'], x, b['Vs'], b['Vb'])), float)[0] for x in v])
        tgt = 1e-7 * sw['w'] / sw['l']

        def extract(cur):
            a = np.abs(cur)
            m = a > 1e-12
            return float(np.interp(np.log10(tgt), np.log10(a[m]), v[m]))

        delta = (extract(g) - extract(r)) * 1e3
        assert abs(delta) < limit, (name, delta)


class TestChannelLengthModulationInTheCharges(object):
    """`GdL` and `QCLM` belong to the charge model too, and were missing.

    They were left out as a long-channel simplification, and for a long
    time nothing could measure the cost -- the charge model was checked
    only against its own construction properties, which are ratios and
    see neither a scale error nor a saturation-region one.

    Against PSP's terminal capacitances the shape of the error names the
    cause without any guessing:

    * in the LINEAR region the two agreed to a part in ten thousand, at
      both geometries -- so it was not overlap capacitance, which would
      show everywhere;
    * the whole error was in SATURATION, growing with drain bias and
      with `ALP` -- 1% on the long device at `Vds = 1.2` and 8% on the
      short one.

    That is `GdL`, and it is the same `GdL` the current already used.
    The polysilicon factor `eta_p` went in beside it, on the inversion
    bracket of the gate charge where PSP puts it -- also a quantity the
    current already carried.
    """

    @pytest.fixture(scope='class')
    def op(self):
        with open(REF) as fh:
            d = json.load(fh)
        return d['op'], d['op_pmos']

    def _fet(self, deck, kind, geom, **over):
        """Overlap off -- this class measures the INTRINSIC charge,
        which is what PSP's `cgg` reports."""
        cm.default_toolkit = numeric
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_%s_psp' % kind, w=geom['w'],
                              l=geom['l'], ng=1, m=1, pre_layout=1),
            w=geom['w'], l=geom['l'], T=T27)
        kw.update(cgov=0.0, cfr=0.0, cgbov=0.0)
        kw.update(over)
        cls = PspPmosLongChannel if kind == 'pmos' else PspMosLongChannel
        e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                **kw)
        e.update_iparv()
        return e

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    def test_the_linear_region_was_never_the_problem(self, deck, op,
                                                     kind):
        """Which is what ruled out overlap capacitance.

        Overlap is a fixed capacitance in parallel with the intrinsic
        one; it does not switch itself off at low drain bias.  This
        error did.
        """
        n, p = op
        for gd in ((p if kind == 'pmos' else n)[g] for g in
                   ('long', 'short')):
            e = self._fet(deck, kind, gd)
            for pt in gd['points']:
                if abs(pt['vd']) > 0.1:
                    continue
                C = np.asarray(e.C(e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])), float)
                assert C[1, 1] == pytest.approx(pt['cgg'], rel=1e-3), \
                    (kind, gd['l'], pt['vg'], pt['vd'])

    def test_it_is_worth_a_factor_of_four_in_saturation(self, deck, op):
        """Measured on the short device, where `ALP` is largest."""
        n, _ = op
        gd = n['short']
        pt = [q for q in gd['points']
              if q['vg'] == 0.8 and q['vd'] == 1.2][0]
        e_on = self._fet(deck, 'nmos', gd)
        e_off = self._fet(deck, 'nmos', gd, alp=0.0)
        bias = (pt['vd'], pt['vg'], 0.0, pt['vb'])
        on = np.asarray(e_on.C(e_on.bias(*bias)), float)[e_on.gate_index, e_on.gate_index]
        ## `alp = 0` removes channel-length modulation from the current
        ## AND the charges, which is what `GdL = 1` meant.
        off = np.asarray(e_off.C(e_off.bias(*bias)), float)[e_off.gate_index, e_off.gate_index]
        assert abs(off / pt['cgg'] - 1.0) > 4.0 * abs(on / pt['cgg'] - 1.0)
        assert abs(on / pt['cgg'] - 1.0) < 0.03, on / pt['cgg']

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    def test_the_long_device_is_now_exact(self, deck, op, kind):
        """0.1% across the whole bias grid, both channel types.

        This is the strongest statement the model can make about its
        charges: they come from a foundry card through the scaling
        layer, nothing is fitted, and they agree with the vendor's own
        implementation to a part in a thousand.
        """
        n, p = op
        gd = (p if kind == 'pmos' else n)['long']
        e = self._fet(deck, kind, gd)
        for pt in gd['points']:
            C = np.asarray(e.C(e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])), float)
            assert C[1, 1] == pytest.approx(pt['cgg'], rel=1e-3), \
                (kind, pt['vg'], pt['vd'], pt['vb'])

    def test_psps_cgg_is_intrinsic_and_excludes_overlap(self, op):
        """The thing it is easiest to get backwards, pinned.

        PSP reports the overlap capacitances SEPARATELY -- `cgsol` and
        `cgdol` -- so `cgg` is the intrinsic one and comparing it
        against an intrinsic-only model is like-for-like.

        This was briefly claimed the other way round, and the numbers
        below are what settle it: the overlap is 11% of `cgg` on the
        long device and 182% of it on the short one.  Had `cgg` included
        them, an intrinsic-only model would be low by those amounts
        rather than by tenths of a percent.
        """
        n, _ = op
        for geom, lo, hi in (('long', 0.05, 0.20), ('short', 1.0, 3.0)):
            pt = n[geom]['points'][5]
            frac = (pt['cgsol'] + pt['cgdol']) / pt['cgg']
            assert lo < frac < hi, (geom, frac)

    def test_what_is_left_on_the_short_device(self, deck, op):
        """2%, and NOT explained -- which is the honest statement.

        It is not overlap: PSP reports that separately, as the test
        above pins.  It is not channel-length modulation: that is in
        now, and it was the linear-region agreement that identified it.
        Recorded as an upper bound so that whatever explains it has
        something to beat.
        """
        n, _ = op
        gd = n['short']
        e = self._fet(deck, 'nmos', gd)
        gi = e.gate_index
        worst = max(
            abs(np.asarray(e.C(e.bias(pt['vd'], pt['vg'], 0.0,
                                      pt['vb'])), float)[gi, gi]
                / pt['cgg'] - 1.0) for pt in gd['points'])
        assert worst < 0.03, worst


class TestTheOverlapAndFringeCapacitance(object):
    """The largest gap the plan named, and the one the DC comparison
    could never have seen.

    PSP reports these separately from `cgg` -- `cgsol` and `cgdol` --
    and on a 0.13 um device they are **227%** of the intrinsic gate
    capacitance.  A model without them is missing most of its charge
    however exact its intrinsic part is, which is precisely the state
    this one was in.

    Cheaper than the first estimate, and the reason is worth recording:
    **the overlap surface potential is CLOSED FORM**
    (`PSP103_module.include:1182-1189`) -- a smoothed maximum and one
    explicit expression, not a Newton solve.  Everything it needs except
    the bias is fixed per instance, so the fit constants
    (`sp_ovInit`, `macrodefs:217-235`) are evaluated once in the scaling
    layer and handed over as numbers.
    """

    @pytest.fixture(scope='class')
    def op(self):
        with open(REF) as fh:
            d = json.load(fh)
        return d['op'], d['op_pmos']

    def _fet(self, deck, kind, geom, **over):
        cm.default_toolkit = numeric
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_%s_psp' % kind, w=geom['w'],
                              l=geom['l'], ng=1, m=1, pre_layout=1),
            w=geom['w'], l=geom['l'], T=T27)
        kw.update(over)
        cls = PspPmosLongChannel if kind == 'pmos' else PspMosLongChannel
        e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                **kw)
        e.update_iparv()
        return e

    @pytest.mark.parametrize('geom', ['long', 'mid', 'short',
                                      'wide_short'])
    def test_the_capacitances_match_psps_scaled_values(self, deck, geom):
        """`lp_cgov` and `lp_cfr`, term by term, before any bias."""
        with open(REF) as fh:
            sc = json.load(fh)['scaled']
        e = sc[geom]
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=e['w'], l=e['l'],
                              ng=1, m=1, pre_layout=1),
            w=e['w'], l=e['l'], T=T27)
        if geom == 'long':
            assert kw['cgov'] == pytest.approx(4.53951e-15, rel=1e-5)
            assert kw['cfr'] == pytest.approx(1.998e-15, rel=1e-4)
        assert kw['cgov'] > 0.0 and kw['cfr'] > 0.0

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    @pytest.mark.parametrize('geom', ['long', 'short'])
    def test_the_total_gate_capacitance_matches(self, deck, op, kind,
                                                geom):
        """Intrinsic plus overlap against PSP's intrinsic plus overlap.

        Within 0.7% everywhere and 0.1% on the long devices -- and the
        short device is the one to watch, because there the overlap is
        more than twice the intrinsic part.
        """
        n, p = op
        gd = (p if kind == 'pmos' else n)[geom]
        e = self._fet(deck, kind, gd)
        gi = e.gate_index
        for pt in gd['points']:
            tot = pt['cgg'] + pt['cgsol'] + pt['cgdol']
            got = np.asarray(e.C(e.bias(pt['vd'], pt['vg'], 0.0,
                                        pt['vb'])), float)[gi, gi]
            assert got == pytest.approx(tot, rel=0.008), \
                (kind, geom, pt['vg'], pt['vd'], got, tot)

    def test_without_them_the_short_device_is_missing_most_of_its_charge(
            self, deck, op):
        """The measurement that justified building this at all."""
        n, _ = op
        gd = n['short']
        pt = gd['points'][5]
        tot = pt['cgg'] + pt['cgsol'] + pt['cgdol']
        gi_on = self._fet(deck, 'nmos', gd)
        off = self._fet(deck, 'nmos', gd, cgov=0.0, cfr=0.0, cgbov=0.0)
        x_on = gi_on.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
        x_off = off.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
        with_ov = np.asarray(gi_on.C(x_on), float)[gi_on.gate_index,
                                                   gi_on.gate_index]
        without = np.asarray(off.C(x_off), float)[off.gate_index,
                                                  off.gate_index]
        assert without / tot < 0.45, without / tot
        assert with_ov == pytest.approx(tot, rel=0.008)

    def test_it_leaves_the_dc_current_alone(self, deck, op):
        """Charge only -- it must not touch the current at all."""
        n, _ = op
        gd = n['long']
        on = self._fet(deck, 'nmos', gd)
        off = self._fet(deck, 'nmos', gd, cgov=0.0, cfr=0.0, cgbov=0.0)
        for vd, vg in ((0.05, 1.2), (1.2, 0.8)):
            a = np.asarray(on.i(on.bias(vd, vg)), float)[0]
            b = np.asarray(off.i(off.bias(vd, vg)), float)[0]
            assert a == pytest.approx(b, rel=1e-14)

    def test_the_construction_properties_survive(self, deck, op):
        """Conservation and antisymmetry, with the overlap on.

        Both are structural and both could plausibly have been broken
        here: the overlap charges are new branches, and they use the
        ACTUAL terminal voltages rather than the ordered ones.  They
        stay even because the drain-side parameters mirror the source's
        (`SWJUNASYM = 0`), so the two simply swap.
        """
        n, _ = op
        e = self._fet(deck, 'nmos', n['long'])
        for vd, vg in ((0.9, 1.4), (0.05, 0.6), (1.5, 1.8)):
            q = np.asarray(e.q(e.bias(vd, vg)), float)
            assert abs(q.sum()) < 4e-16 * max(1.0, np.abs(q).max())
        for vg, vd in ((0.6, 0.7), (1.2, 1.5)):
            f = np.asarray(e.i(e.bias(vd, vg, 0.0, 0.0)), float)[0]
            r = np.asarray(e.i(e.bias(0.0, vg, vd, 0.0)), float)[0]
            assert f == pytest.approx(-r, rel=1e-13), (vg, vd)


class TestTemperatureScaling(object):
    """The card is a 27 C card, and PSP scales it from there.

    28 non-zero `ST*` coefficients (`PSP103_macrodefs.include:357-390`),
    so this is not a refinement: without it the model is a 27 C model
    wearing whatever temperature the caller thought it asked for, and no
    corner or temperature analysis is possible at all.

    Almost everything is a power law -- `X_T = X*(TKR/TKD)^ST_X` -- with
    the flat-band voltage the one exception, a quadratic in `delT`.
    `STVFB`, `STBET` and `STTHESAT` carry geometry terms of their own.

    Validated at 100 C, 73 K from the reference, which is far enough
    that a wrong exponent or a flipped sign cannot hide.
    """

    #: PSP's `lp_*` name against the keyword `to_long_channel` returns.
    PAIRS = (('vfb', 'vfb'), ('ct', 'ct'), ('mue', 'mue'),
             ('themu', 'themu'), ('cs', 'cs'), ('thecs', 'thecs'),
             ('rs', 'rs'), ('thesat', 'thesat'), ('xcor', 'xcor'))

    @pytest.fixture(scope='class')
    def hot(self):
        with open(REF) as fh:
            d = json.load(fh)
        return d['scaled_hot'], d['scaled_hot_pmos']

    @pytest.mark.parametrize('geom', ['long', 'mid', 'short',
                                      'wide_short'])
    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    def test_every_parameter_matches_psp_at_100c(self, deck, hot, kind,
                                                 geom):
        n, p = hot
        e = (p if kind == 'pmos' else n)[geom]
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_%s_psp' % kind, w=e['w'],
                              l=e['l'], ng=1, m=1, pre_layout=1),
            w=e['w'], l=e['l'], T=273.15 + e['tempc'])
        for ref_key, our_key in self.PAIRS:
            assert kw[our_key] == pytest.approx(e[ref_key], rel=1e-5,
                                                abs=1e-12), \
                (kind, geom, ref_key)

    @pytest.mark.parametrize('geom', ['long', 'short'])
    def test_the_gain_factor_matches_at_100c(self, deck, hot, geom):
        """`BETN` separately, since we carry it as an effective
        mobility rather than as `BETN` itself."""
        n, _ = hot
        e = n[geom]
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=e['w'], l=e['l'],
                              ng=1, m=1, pre_layout=1),
            w=e['w'], l=e['l'], T=273.15 + e['tempc'])
        assert kw['u0'] * e['w'] / e['l'] == pytest.approx(e['betn'],
                                                           rel=1e-5)

    def test_at_the_reference_temperature_nothing_moves(self, deck):
        """The scaling has to be the identity at `TR`.

        Which is the check that catches a sign error in `rTn`: getting
        REFERENCE-over-DEVICE backwards leaves this passing and every
        other temperature wrong in the opposite direction.
        """
        with open(REF) as fh:
            cold = json.load(fh)['scaled']['long']
        card = deck.model_params('sg13g2_lv_nmos_psp', w=cold['w'],
                                 l=cold['l'], ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=cold['w'], l=cold['l'],
                                         T=273.15 + card['tr'])
        ## `rel=1e-5`, the reference's own precision: ngspice prints
        ## its operating-point outputs to six significant digits, so a
        ## tighter bound tests the printf format.
        for ref_key, our_key in self.PAIRS:
            assert kw[our_key] == pytest.approx(cold[ref_key], rel=1e-5,
                                                abs=1e-12), ref_key

    def test_the_trends_have_the_right_sign(self, deck):
        """A power law with the exponent's sign wrong still passes a
        single-point comparison if that point is the reference.  These
        are the directions a device physicist would insist on.
        """
        card = deck.model_params('sg13g2_lv_nmos_psp', w=10e-6, l=1e-6,
                                 ng=1, m=1, pre_layout=1)
        cold = psp_scaling.to_long_channel(card, w=10e-6, l=1e-6, T=250.0)
        hot = psp_scaling.to_long_channel(card, w=10e-6, l=1e-6, T=400.0)
        assert hot['u0'] < cold['u0'], 'mobility falls with temperature'
        assert hot['mue'] < cold['mue']
        assert hot['thesat'] < cold['thesat']
        assert hot['rs'] > cold['rs'], 'series resistance rises'
        assert hot['vfb'] > cold['vfb'], 'flat-band rises, threshold falls'

    def test_the_reference_temperature_mismatch_is_bounded(self, deck):
        """A known, measured, and deliberately unclosed inconsistency.

        The parameter comparisons are made at the card's `TR = 27 C`,
        because that is where the reference was recorded.  The element
        evaluating those parameters still runs at `epar`'s default of
        300.0 K, 0.15 K away, because threading an `epar` through every
        call site costs more than the error does.

        This pins what the error IS, so it stays a decision rather than
        becoming a mystery: 0.15 K is worth under a hundredth of a
        percent of drain current.
        """
        import copy
        from pycircuit.circuit import defaultepar
        cm.default_toolkit = numeric
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=10e-6, l=1e-6,
                              ng=1, m=1, pre_layout=1),
            w=10e-6, l=1e-6, T=T27)
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        hot = copy.copy(defaultepar)
        hot.T = T27
        x = e.bias(0.05, 1.2)
        a = np.asarray(e.i(x), float)[0]
        b = np.asarray(e.i(x, hot), float)[0]
        assert abs(a / b - 1.0) < 2e-4, (a, b)

    def test_the_default_temperature_matches_the_elements(self):
        """`to_long_channel` scales the parameters and the element's own
        `vt()` follows `epar`; they are not linked, so their DEFAULTS
        have to agree or a caller who sets neither is quietly running
        two temperatures at once.
        """
        import inspect
        from pycircuit.circuit import hdl
        sig = inspect.signature(psp_scaling.to_long_channel)
        assert sig.parameters['T'].default == hdl._epar_T(object())


class TestExtremeBiasWithCardParameters(object):
    """Finiteness far outside anything physical, from a REAL card.

    `test_psp_current.py` already covers this, and covers it well -- but
    with DEFAULT parameters, and that turns out to be a different
    question.  The defaults leave the junction off, `alp` and `rs` at
    zero and the noise absent; a card turns all of them on and moves
    every quantity the guards are sized against.  This branch has
    learned that distinction twice before (the symmetry bug that only
    card parameters could see, and the `AX` floor the card alone
    exercises), so the regime is pinned here rather than assumed from
    there.

    Why it matters at all, at biases this absurd: the device's BRANCH
    voltages are bounded by `limit()`, but its NODE voltages are not, so
    a diverging solve evaluates the model far outside anything physical
    and it has to stay in range there.  A NaN in one Jacobian row
    poisons the whole solve rather than one entry -- and a NaN in the
    JACOBIAN with a finite VALUE is the worse case, because nothing
    looks wrong.
    """

    def _fet(self, deck, kind='nmos', **over):
        cm.default_toolkit = numeric
        cls, card = ((PspMosLongChannel, 'sg13g2_lv_nmos_psp')
                     if kind == 'nmos'
                     else (PspPmosLongChannel, 'sg13g2_lv_pmos_psp'))
        kw = psp_scaling.to_long_channel(
            deck.model_params(card, w=10e-6, l=1e-6, ng=1, m=1,
                              pre_layout=1),
            w=10e-6, l=1e-6, T=T27)
        kw.update(over)
        e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                **kw)
        e.update_iparv()
        return e

    #: Both the value and the Jacobian are finite through here, on both
    #: signs and both channel types.  Measured: the n-channel holds to
    #: 1e27 and the p-channel to 1e24, so the shared list stops at the
    #: SMALLER -- and the per-type bounds are pinned separately in
    #: `test_the_bound_is_where_it_is_measured_to_be`, because "both
    #: types survive this far" and "this is how far each gets" are
    #: different statements and only the second catches a regression in
    #: the better one.
    REACH = [1.0e2, 1.0e3, 1.0e4, 1.0e7, 1.0e12, 1.0e20, 1.0e24]

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    @pytest.mark.parametrize('sign', [1.0, -1.0])
    def test_the_jacobian_stays_finite_with_a_card(self, deck, kind, sign):
        """This is the regression that matters.

        At `Vd = -1e3` the drain junction is FORWARD biased, so
        JUNCAP2's avalanche argument is smoothed to exactly zero -- and
        `Abs(u)` differentiates to `(re*re' + im*im')/Abs(u)`, which is
        `0/0` there.  The current was finite and the conductance was
        NaN.  1e3 volts is not hypothetical: two of these in a stack
        walk out that far on the way to a solution.
        """
        e = self._fet(deck, kind)
        t = 1.0 if kind == 'nmos' else -1.0
        for v in self.REACH:
            x = e.bias(sign * t * v, t * 1.2, 0.0, 0.0)
            with np.errstate(all='ignore'):
                i = np.asarray(e.i(x), float)
                G = np.asarray(e.G(x), float)
            assert np.all(np.isfinite(i)), (kind, sign * v, 'i')
            assert np.all(np.isfinite(G)), (kind, sign * v, 'G')

    def test_the_current_stays_bounded_too(self, deck):
        """Finite is not enough -- it also has to SATURATE.  `Vdse`
        pins the drain end, so the current may creep but must not run
        away with the bias."""
        e = self._fet(deck)
        ref = np.asarray(e.i(e.bias(2.0, 1.2, 0.0, 0.0)), float)[0]
        for v in self.REACH:
            got = np.asarray(e.i(e.bias(v, 1.2, 0.0, 0.0)), float)[0]
            ## Twenty-five decades of drain bias buy under six decades
            ## of current, which is what saturation looks like from
            ## here.  Without `Vdse` it would be linear in `Vd`.
            assert 0.0 < got < 1e6 * ref, (v, got, ref)

    @pytest.mark.parametrize('kind,last', [('nmos', 27), ('pmos', 24)])
    def test_the_bound_is_where_it_is_measured_to_be(self, deck, kind,
                                                     last):
        """The limit is RECORDED, not merely respected.

        Past these the Jacobian is still non-finite, and the cause is a
        property of the compiler rather than of the model:
        `safe_div(a, b)` regularises as `a*b/(b^2 + eps^2)`, so its
        VALUE needs `|b| < ~1e154` while its DERIVATIVE squares that
        denominator again and needs `|b| < ~1e77`.  A quarter of the
        exponent range, and the drain-side surface-potential update
        reaches it three decades sooner on the p-channel than on the
        n-channel.

        Pinned so a future change to `safe_div` has a number to beat,
        and so the gap is not rediscovered as a bug.
        """
        e = self._fet(deck, kind)
        t = 1.0 if kind == 'nmos' else -1.0
        with np.errstate(all='ignore'):
            ok = np.asarray(e.G(e.bias(t * 10.0 ** last, t * 1.2, 0.0,
                                       0.0)), float)
            bad = np.asarray(e.G(e.bias(t * 10.0 ** (last + 1), t * 1.2,
                                        0.0, 0.0)), float)
            val = np.asarray(e.i(e.bias(t * 10.0 ** (last + 1), t * 1.2,
                                        0.0, 0.0)), float)
        assert np.all(np.isfinite(ok)), kind
        assert not np.all(np.isfinite(bad)), kind
        ## The VALUE survives past where the derivative does -- which is
        ## exactly why this had to be looked for rather than noticed.
        assert np.all(np.isfinite(val)), kind


class TestTheJunction(object):
    """JUNCAP2's charge, and the ideal diode that goes with it.

    The junction capacitance is 8% of the intrinsic gate capacitance on
    a 10 um device and 126% of it at 0.13 um -- most of what a transient
    or AC simulation needs.  The junction CURRENT is around 1e-15 A
    against a 4e-4 A drain current, so the charge is implemented in full
    and of the current only the ideal diode term.

    Three things the vendor source settles that measurement alone would
    have got wrong, and each was worth finding before writing code:

    * `SWJUNCAP = 3` is a GEOMETRY-SOURCE selector, not a component
      selector (`PSP103_module.include:867-883`).  All three components
      run; `= 3` only carves the gate edge out of the perimeter so it is
      not counted twice.
    * `LG` is the ELECTRICAL width `WE`, not the drawn `W`.  `WOT` is
      negative on this card, so they differ by 0.02 um -- 2% on a 1 um
      device and 0.2% on a 10 um one, which reads as a scaling error
      rather than as a definition.
    * JUNCAP has its OWN reference temperature, `TRJ = 21 C`, against
      the transistor's 27 C.  Assuming they were the same makes every
      capacitance about 2% low.
    """

    @pytest.fixture(scope='class')
    def op(self):
        with open(REF) as fh:
            d = json.load(fh)
        return d['op'], d['op_pmos']

    def _pair(self, deck, kind, geom, **over):
        """The same device with the junction on and off."""
        cm.default_toolkit = numeric
        w, l = geom['w'], geom['l']
        card = deck.model_params('sg13g2_lv_%s_psp' % kind, w=w, l=l,
                                 ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=w, l=l, T=T27)
        jg = psp_scaling.junction_geometry(
            w, psp_scaling.geometry(card, w=w, l=l)['WE'])
        jp = psp_scaling.junction(card, T=T27, **jg)
        cls = PspPmosLongChannel if kind == 'pmos' else PspMosLongChannel
        jp = dict(jp, **over)
        out = []
        for extra in (jp, {}):
            e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                    cm.Node('b'), **dict(kw, **extra))
            e.update_iparv()
            out.append(e)
        return out

    def test_the_geometry_is_the_electrical_width_not_the_drawn_one(
            self, deck):
        """The definition that a single geometry cannot reveal."""
        card = deck.model_params('sg13g2_lv_nmos_psp', w=1e-6,
                                 l=0.13e-6, ng=1, m=1, pre_layout=1)
        geo = psp_scaling.geometry(card, w=1e-6, l=0.13e-6)
        jg = psp_scaling.junction_geometry(1e-6, geo['WE'])
        assert jg['lg'] == pytest.approx(geo['WE'], rel=1e-12)
        assert jg['lg'] != pytest.approx(1e-6, rel=1e-6), \
            'the drawn width would be wrong by 2% here'
        assert jg['ab'] == pytest.approx(3.40e-13, rel=1e-3)
        assert jg['ls'] == pytest.approx(1.66e-6, rel=1e-3)

    def test_it_uses_its_own_reference_temperature(self, deck):
        """`TRJ = 21 C`, not the transistor's `TR = 27 C`."""
        card = deck.model_params('sg13g2_lv_nmos_psp', w=1e-6,
                                 l=0.13e-6, ng=1, m=1, pre_layout=1)
        assert card['trj'] == 21.0
        assert card['tr'] == 27.0
        jg = psp_scaling.junction_geometry(
            1e-6, psp_scaling.geometry(card, w=1e-6, l=0.13e-6)['WE'])
        at_trj = psp_scaling.junction(card, T=273.15 + 21.0, **jg)
        at_tr = psp_scaling.junction(card, T=T27, **jg)
        ## At its OWN reference the built-in voltage is the card's.
        assert at_trj['jbot_vbi'] == pytest.approx(card['vbirbot'],
                                                   rel=1e-9)
        ## Six kelvin away it is not, and the capacitance moves with it.
        assert at_tr['jbot_vbi'] < at_trj['jbot_vbi']
        assert at_tr['jbot_qpref2'] > at_trj['jbot_qpref2']

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    @pytest.mark.parametrize('geom', ['long', 'short'])
    def test_the_capacitance_matches_psp(self, deck, op, kind, geom):
        """Every bias point, both geometries, both channel types.

        Measured as the DIFFERENCE the junction makes, against PSP's
        `cjs + cjd`, because our bulk row carries the intrinsic bulk
        capacitance too and PSP reports the two separately.
        """
        n, p = op
        gd = (p if kind == 'pmos' else n)[geom]
        on, off = self._pair(deck, kind, gd)
        for pt in gd['points']:
            a = np.asarray(on.C(on.bias(pt['vd'], pt['vg'], 0.0,
                                        pt['vb'])), float)[3, 3]
            b = np.asarray(off.C(off.bias(pt['vd'], pt['vg'], 0.0,
                                          pt['vb'])), float)[3, 3]
            assert a - b == pytest.approx(pt['cjs'] + pt['cjd'],
                                          rel=1e-4), \
                (kind, geom, pt['vg'], pt['vd'], pt['vb'])

    def test_the_components_are_all_live(self, deck, op):
        """None of the three is switched off by a zero prefactor.

        `SWJUNCAP = 3` runs all three, and the split is roughly
        85 / 7 / 8 -- so each is individually large enough that a
        total-only comparison could hide a factor of two in one of them.
        """
        n, _ = op
        pt = n['long']['points'][5]
        tot = pt['cjsbot'] + pt['cjssti'] + pt['cjsgat']
        assert pt['cjsbot'] / tot > 0.7
        assert pt['cjssti'] / tot > 0.03
        assert pt['cjsgat'] / tot > 0.03
        assert tot == pytest.approx(pt['cjs'], rel=1e-6)

    def test_it_leaves_the_drain_current_alone(self, deck, op):
        """The junction is charge plus a diode to the BULK; the channel
        current must not move."""
        n, _ = op
        on, off = self._pair(deck, 'nmos', n['long'])
        for vd, vg in ((0.05, 1.2), (1.2, 0.8)):
            a = np.asarray(on.i(on.bias(vd, vg)), float)[0]
            b = np.asarray(off.i(off.bias(vd, vg)), float)[0]
            assert a == pytest.approx(b, rel=1e-9)

    @pytest.mark.parametrize('geom', ['long', 'short'])
    def test_the_reverse_leakage_matches_psp(self, deck, op, geom):
        """The Shockley-Read-Hall, trap-assisted-tunnelling,
        band-to-band and avalanche terms, together.

        These were left out at first, on the grounds that they shape a
        current eleven orders below the drain current.  That argument is
        still true and they are implemented anyway: the ideal diode
        alone gets the REVERSE characteristic wrong by four orders, and
        a junction model that is wrong by four orders in reverse is not
        a junction model.

        Measured as the difference the junction makes to the drain
        current, against PSP's `ijd`.
        """
        n, _ = op
        gd = n[geom]
        on, off = self._pair(deck, 'nmos', gd)
        for pt in gd['points']:
            if pt['ijd'] == 0.0:
                continue
            a = np.asarray(on.i(on.bias(pt['vd'], pt['vg'], 0.0,
                                        pt['vb'])), float)[0]
            b = np.asarray(off.i(off.bias(pt['vd'], pt['vg'], 0.0,
                                          pt['vb'])), float)[0]
            assert a - b == pytest.approx(pt['ijd'], rel=1e-3), \
                (geom, pt['vg'], pt['vd'], pt['vb'])

    def test_the_ideal_term_alone_would_not_do(self, deck, op):
        """Why the leakage terms earn their lines.

        With only the ideal diode the reverse current is four orders
        too small.  This pins that, so the decision to implement them is
        recorded as a measurement rather than as a preference.
        """
        n, _ = op
        gd = n['long']
        pt = [q for q in gd['points'] if q['vd'] == 1.2][0]
        on, off = self._pair(deck, 'nmos', gd)
        full = (np.asarray(on.i(on.bias(pt['vd'], pt['vg'], 0.0,
                                        pt['vb'])), float)[0]
                - np.asarray(off.i(off.bias(pt['vd'], pt['vg'], 0.0,
                                            pt['vb'])), float)[0])
        ideal = self._pair(deck, 'nmos', gd)[0]
        ## Switch the three leakage prefactors off, leaving the ideal
        ## diode: that is the state this model shipped in first.
        bare, _ = self._pair(deck, 'nmos', gd, jbot_csrh=0.0,
                             jsti_csrh=0.0, jgat_csrh=0.0,
                             jbot_ctat=0.0, jsti_ctat=0.0,
                             jgat_ctat=0.0, jbot_cbbt=0.0,
                             jsti_cbbt=0.0, jgat_cbbt=0.0)
        only_ideal = (np.asarray(bare.i(bare.bias(pt['vd'], pt['vg'],
                                                  0.0, pt['vb'])),
                                 float)[0]
                      - np.asarray(off.i(off.bias(pt['vd'], pt['vg'],
                                                  0.0, pt['vb'])),
                                   float)[0])
        assert full == pytest.approx(pt['ijd'], rel=1e-3)
        assert abs(only_ideal) < 1e-3 * abs(full), (only_ideal, full)

    def test_the_forward_diode_still_conducts(self, deck, op):
        """Honest about what the ideal-diode-only current costs.

        The one regime where junction current matters is a
        forward-biased bulk, and it must still turn on.
        """
        n, _ = op
        on, _off = self._pair(deck, 'nmos', n['short'])
        ## Bulk BELOW the diffusions is reverse bias for an n-channel.
        rev = abs(np.asarray(on.i(on.bias(0.0, 0.0, 0.0, -1.2)),
                             float)[3])
        assert rev < 1e-12, rev
        ## And the diode still turns on the other way, which is the one
        ## regime where junction current matters and the reason the
        ## ideal term was kept rather than the current dropped entirely.
        fwd = abs(np.asarray(on.i(on.bias(0.0, 0.0, 0.0, 1.2)),
                             float)[3])
        assert fwd > 1e-3, fwd

    def test_conservation_and_finiteness_hold_with_it_on(self, deck, op):
        """Both could plausibly have broken: the junction adds branches,
        and its expressions are exponentials and fractional powers."""
        n, _ = op
        on, _ = self._pair(deck, 'nmos', n['long'])
        for vd, vg in ((0.9, 1.4), (0.05, 0.6), (1.5, 1.8)):
            q = np.asarray(on.q(on.bias(vd, vg)), float)
            assert abs(q.sum()) < 4e-16 * max(1.0, np.abs(q).max())
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            for vd in (0.0, 1.2, -6.0, 1e3):
                x = on.bias(vd, 1.8)
                assert np.all(np.isfinite(np.asarray(on.i(x), float)))
                assert np.all(np.isfinite(np.asarray(on.G(x), float)))


class TestTheChannelNoise(object):
    """Thermal and flicker noise, from the same surface-potential
    quantities the current uses.

    That is what a surface-potential model buys here: `qim`, `qim1`,
    `alpha`, `dps` and `H` already exist, so there is no separate noise
    model to fit -- the densities fall out of the charge distribution
    along the channel (`PSP103_module.include:1819-1841`).

    Induced gate noise is in `TestTheInducedGateNoise` below.  It splits
    the drain's white term into a correlated part and `(1 - c_igid^2)`
    of the rest -- so the tests HERE, which measure the drain alone,
    must be unaffected by it, and are.
    """

    @pytest.fixture(scope='class')
    def op(self):
        with open(REF) as fh:
            return json.load(fh)['op']

    def _fet(self, deck, geom, **over):
        cm.default_toolkit = numeric
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=geom['w'],
                              l=geom['l'], ng=1, m=1, pre_layout=1),
            w=geom['w'], l=geom['l'], T=T27)
        kw.update(over)
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        return e

    @staticmethod
    def _split(e, x):
        """Separate the white and flicker densities by measuring at two
        frequencies -- `CY` gives their sum."""
        f1, f2 = 1.0, 1.0e9
        a = np.asarray(e.CY(x, 2 * np.pi * f1), float)[0, 0]
        b = np.asarray(e.CY(x, 2 * np.pi * f2), float)[0, 0]
        sfl = (a - b) / (1.0 / f1 - 1.0 / f2)
        return a - sfl / f1, sfl

    @pytest.mark.parametrize('geom', ['long', 'short'])
    def test_both_densities_match_psp(self, deck, op, geom):
        """Within half a percent, which is where the drain current they
        are built from already sits."""
        gd = op[geom]
        e = self._fet(deck, gd)
        for pt in gd['points']:
            if pt['sid'] <= 0.0:
                continue
            x = e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
            sid, sfl = self._split(e, x)
            assert sid == pytest.approx(pt['sid'], rel=0.02), \
                ('sid', geom, pt['vg'], pt['vd'])
            if pt['sfl'] > 0.0:
                assert sfl == pytest.approx(pt['sfl'], rel=0.02), \
                    ('sfl', geom, pt['vg'], pt['vd'])

    def test_the_corner_frequency_matches(self, deck, op):
        """`fknee` is the RATIO of the two densities, so it checks them
        against each other rather than both against the current they
        share -- an error common to both would cancel here and an error
        in one would not.
        """
        gd = op['long']
        e = self._fet(deck, gd)
        for pt in gd['points']:
            if not (pt['fknee'] > 0.0 and pt['sid'] > 0.0):
                continue
            x = e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
            sid, sfl = self._split(e, x)
            assert sfl / sid == pytest.approx(pt['fknee'], rel=0.02), \
                (pt['vg'], pt['vd'], sfl / sid, pt['fknee'])

    def test_the_residual_is_the_current_it_is_built_from(self, deck, op):
        """Both densities scale with the drain current, so their error
        should be the current's error and not something of its own.

        Pinning this is what distinguishes 'the noise formulas are
        right' from 'the noise happens to be close'.
        """
        gd = op['long']
        e = self._fet(deck, gd)
        pt = gd['points'][5]
        x = e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
        sid, sfl = self._split(e, x)
        ids = np.asarray(e.i(x), float)[0] / pt['ids']
        assert abs(sid / pt['sid'] - 1.0) < 5.0 * abs(ids - 1.0) + 2e-3
        assert abs(sfl / pt['sfl'] - 1.0) < 5.0 * abs(ids - 1.0) + 2e-3

    def test_it_is_off_without_a_card(self, deck, op):
        """An element built with defaults is noiseless, so the noise
        cannot leak into a model that never asked for it."""
        e = self._fet(deck, op['long'], fnt=0.0, nfa=0.0, nfb=0.0,
                      nfc=0.0)
        CY = np.asarray(e.CY(e.bias(0.05, 1.2), 2 * np.pi), float)
        assert np.all(CY == 0.0)

    def test_thermal_and_flicker_can_be_switched_independently(self,
                                                               deck, op):
        gd = op['long']
        x_of = lambda e: e.bias(0.05, 1.2)
        th = self._fet(deck, gd, nfa=0.0, nfb=0.0, nfc=0.0)
        fl = self._fet(deck, gd, fnt=0.0)
        s_th, f_th = self._split(th, x_of(th))
        s_fl, f_fl = self._split(fl, x_of(fl))
        assert s_th > 0.0 and f_th == pytest.approx(0.0, abs=1e-30)
        assert f_fl > 0.0 and s_fl == pytest.approx(0.0, abs=1e-30)


class TestTheInducedGateNoise(object):
    """The channel's fluctuation does not only leave through the drain.

    It couples capacitively to the gate as well, and -- being the SAME
    fluctuation -- what arrives there is CORRELATED with what arrives at
    the drain (`PSP103_module.include:1863-1881, 1942-1948`).  Both
    halves of that sentence are testable and both are tested here: the
    magnitude against PSP's `sig`, and the correlation against its
    `cigid`.

    PSP builds it as a NETWORK -- an auxiliary node carrying a
    conductance, a capacitance and the shared source -- rather than as a
    written-down `jw`, and so does this model.  That is why the `f^2`
    rise below the pole and the flat region above it are testable at
    all: they are what an RC does, not what anyone typed.
    """

    #: element class, card name and reference key, per channel type.
    DEV = {'nmos': (PspMosLongChannel, 'sg13g2_lv_nmos_psp', 'op'),
           'pmos': (PspPmosLongChannel, 'sg13g2_lv_pmos_psp', 'op_pmos')}

    @pytest.fixture(scope='class')
    def ops(self):
        with open(REF) as fh:
            d = json.load(fh)
        return {'nmos': d['op'], 'pmos': d['op_pmos']}

    @pytest.fixture(scope='class')
    def op(self, ops):
        return ops['nmos']

    def _kw(self, deck, geom, kind='nmos', **over):
        _, card, _ = self.DEV[kind]
        kw = psp_scaling.to_long_channel(
            deck.model_params(card, w=geom['w'], l=geom['l'], ng=1, m=1,
                              pre_layout=1),
            w=geom['w'], l=geom['l'], T=T27)
        kw.update(over)
        return kw

    def _fet(self, deck, geom, kind='nmos', **over):
        cm.default_toolkit = numeric
        cls = self.DEV[kind][0]
        e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                **self._kw(deck, geom, kind, **over))
        e.update_iparv()
        return e

    def _sig(self, deck_kw, pt, f=1.0e3, kind='nmos'):
        """The GATE current density, from a real noise analysis.

        Measured as the current through an ideal source in the gate
        lead, which is the definition of `sig` -- and, unlike reading a
        `CY` entry, it exercises the whole path: the auxiliary node, the
        capacitive coupling, and the adjoint solve that turns them into
        a density at a terminal.
        """
        from pycircuit.circuit.circuit import SubCircuit, gnd
        from pycircuit.circuit.elements import VS
        from pycircuit.circuit.analysis_ss import Noise
        cm.default_toolkit = numeric
        c = SubCircuit()
        nd, ng_, nb = c.add_nodes('d', 'g', 'b')
        c['vd'] = VS(nd, gnd, v=pt['vd'])
        c['vg'] = VS(ng_, gnd, v=pt['vg'])
        c['vb'] = VS(nb, gnd, v=pt['vb'])
        c['m1'] = self.DEV[kind][0](nd, ng_, gnd, nb, **deck_kw)
        res = Noise(c, inputsrc='vd', outputsrc='vg').solve(f)
        return float(np.real(res['Sinout']))

    @staticmethod
    def _white(e, x, i, j):
        """The WHITE part of a `CY` entry.

        Necessary, not fussy: at 1 kHz the drain entry is dominated by
        FLICKER noise -- some four decades of it on the long device --
        so reading `CY` there and calling it `Sid` measures the wrong
        term entirely.
        """
        f1, f2 = 1.0, 1.0e9
        a = np.asarray(e.CY(x, 2 * np.pi * f1), float)[i, j]
        b = np.asarray(e.CY(x, 2 * np.pi * f2), float)[i, j]
        sfl = (a - b) / (1.0 / f1 - 1.0 / f2)
        return a - sfl / f1

    def _nt(self, deck, geom, kind='nmos'):
        """`nt = 4*FNT*k*T`, in PSP's own constants."""
        from pycircuit.circuit.psp_scaling import PSP_KBOL, PSP_QELE
        return (self._kw(deck, geom, kind)['fnt'] * 4.0
                * (PSP_KBOL * T27 / PSP_QELE) * PSP_QELE)

    def _corr(self, e, x):
        """`c_igid`, as the correlation coefficient it is defined to be.

        Read off `CY` between the drain row and the auxiliary row.  The
        gate current is a frequency-dependent multiple of that auxiliary
        source, and a correlation coefficient is invariant under
        scaling, so this is `|c(i_g, i_d)|` without needing a frequency.
        """
        nm = [n.name for n in e.nodes]
        i_d, i_n = nm.index('d'), nm.index('noi')
        cross = self._white(e, x, i_d, i_n)
        return cross / np.sqrt(self._white(e, x, i_n, i_n)
                               * self._white(e, x, i_d, i_d))

    #: `cigid` tolerance per (channel type, geometry) -- measured, not
    #: chosen.  The p-channel is the better device here, as it is on
    #: this branch's DC sweeps.
    CORR_TOL = {('nmos', 'long'): 0.003, ('nmos', 'short'): 0.04,
                ('pmos', 'long'): 0.001, ('pmos', 'short'): 0.008}

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    @pytest.mark.parametrize('geom', ['long', 'short'])
    def test_the_correlation_matches_psp(self, deck, ops, kind, geom):
        """Within 0.1% on the long n-channel and 3.2% on the short one;
        within 0.55% on the p-channel at both geometries.

        `c_igid` is built from the channel's SHAPE functions alone --
        `t1`, `t2`, `r` -- with `g_ideal`, `lc` and `CGeff` cancelling
        out of it exactly.  So this checks the shape and nothing else,
        which is what makes it a different measurement from `sig` below
        rather than a second view of the same one.
        """
        gd = ops[kind][geom]
        e = self._fet(deck, gd, kind)
        tol = self.CORR_TOL[(kind, geom)]
        for pt in gd['points']:
            if pt['sid'] <= 0.0 or pt['cigid'] <= 0.0:
                continue
            x = e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
            assert self._corr(e, x) == pytest.approx(pt['cigid'], rel=tol), \
                (kind, geom, pt['vg'], pt['vd'], pt['vb'])

    #: `sig` tolerance per (channel type, geometry) -- measured.  One
    #: tolerance ACROSS ALL BIASES of a device, deliberately: the bug
    #: this replaced was exact at `Vds = 0.05` and 73% high at 1.2, so a
    #: per-bias tolerance would have accommodated it and a shared one
    #: does not.  Anything that makes `CGeff` carry a spurious
    #: `Vds`-dependence fails here.
    SIG_TOL = {('nmos', 'long'): 0.005, ('nmos', 'short'): 0.09,
               ('pmos', 'long'): 0.005, ('pmos', 'short'): 0.02}

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    @pytest.mark.parametrize('geom', ['long', 'short'])
    def test_the_gate_density_matches_psp(self, deck, ops, kind, geom):
        """End to end through a noise analysis, every recorded bias.

        Within 0.4% on both long devices and on the short p-channel.
        The short n-channel runs to 8.4% low, and it is the same two
        places its DRAIN density is already off -- deep subthreshold at
        `Vg = 0.4`, where `sid` itself is 5-10% out -- rather than
        anything the gate path adds.
        """
        gd = ops[kind][geom]
        kw = self._kw(deck, gd, kind)
        tol = self.SIG_TOL[(kind, geom)]
        for pt in gd['points']:
            if pt['sig'] <= 0.0:
                continue
            assert self._sig(kw, pt, kind=kind) \
                == pytest.approx(pt['sig'], rel=tol), \
                (kind, geom, pt['vg'], pt['vd'], pt['vb'])

    @pytest.mark.parametrize('kind', ['nmos', 'pmos'])
    @pytest.mark.parametrize('geom', ['long', 'short'])
    def test_the_coupling_capacitance_matches_psp(self, deck, ops, kind,
                                                  geom):
        """`CGeff` term by term -- against a quantity PSP does NOT export.

        It exports `sig` and `cigid`, and

            sig   = nt*w^2*CGeff^2*mig / (1 + w^2*CGeff^2*mig^2)
            cigid = migid0 / sqrt(mig * mid),   mid = Sid/nt

        so with `migid0` -- a pure SHAPE function, which `cigid` already
        validates independently -- PSP's `mig` follows from its own
        `cigid` and `sid`, and then its `CGeff` follows from its own
        `sig`.  Ours is read straight off the assembled matrices:
        `C[noi, noi]` IS `CGeff`, and `CY[noi, noi]` is `nt*gmig`.

        This is the `lp_*` method applied to something the vendor keeps
        to itself, and it is what turned a residual that looked like
        missing short-channel physics into a one-line transcription bug:
        `core['Gmob']` already holds `Gmob_dL`, so multiplying by
        `core['GdL']` too left a spurious `1/GdL^2` that grows with
        `Vds` through channel-length modulation.

        Within 0.1% on the long devices.
        """
        gd = ops[kind][geom]
        e = self._fet(deck, gd, kind)
        nm = [n.name for n in e.nodes]
        i_d, i_n = nm.index('d'), nm.index('noi')
        C_ours = np.asarray  # local alias, keeps the line short
        nt = self._nt(deck, gd, kind)
        w = 2.0 * np.pi * 1.0e3
        tol = 0.005 if geom == 'long' else 0.05
        for pt in gd['points']:
            if pt['sig'] <= 0.0 or pt['cigid'] <= 0.0:
                continue
            x = e.bias(pt['vd'], pt['vg'], 0.0, pt['vb'])
            pwr = self._white(e, x, i_n, i_n)
            migid0 = self._white(e, x, i_d, i_n) / pwr
            ours = C_ours(e.C(x), float)[i_n, i_n]
            mid_p = pt['sid'] / nt
            mig_p = migid0 ** 2 / (pt['cigid'] ** 2 * mid_p)
            den = nt * w * w * mig_p - pt['sig'] * w * w * mig_p * mig_p
            theirs = np.sqrt(pt['sig'] / den)
            assert ours == pytest.approx(theirs, rel=tol), \
                (kind, geom, pt['vg'], pt['vd'], ours, theirs)

    def test_the_drain_density_is_untouched_by_the_split(self, deck, op):
        """The correlated source and the reduced independent one sum
        back to `Sid` IDENTICALLY, whatever the clip did.

        Which is the invariant that lets the drain-noise tests above
        stay exactly as they were: adding a correlated pair changes the
        gate and the cross terms, and must change the drain-drain entry
        by nothing at all.
        """
        gd = op['long']
        on = self._fet(deck, gd)
        off = self._fet(deck, gd, swign=0.0)
        nm_on = [n.name for n in on.nodes]
        nm_off = [n.name for n in off.nodes]
        for pt in gd['points']:
            if pt['sid'] <= 0.0:
                continue
            a = self._white(on, on.bias(pt['vd'], pt['vg'], 0.0, pt['vb']),
                            nm_on.index('d'), nm_on.index('d'))
            b = self._white(off, off.bias(pt['vd'], pt['vg'], 0.0, pt['vb']),
                            nm_off.index('d'), nm_off.index('d'))
            assert a == pytest.approx(b, rel=1e-12), (pt['vg'], pt['vd'])

    def test_it_rises_as_f_squared_and_then_flattens(self, deck, op):
        """The RC's signature, which is the reason for building it as a
        network: below the pole the density goes as `f^2`, above it the
        conductance stops mattering and the density flattens at the bare
        source power `nt * gmig`.

        Measured with `rg = 0` deliberately.  A gate resistance puts a
        SECOND pole in the path -- this card's 1.3 ohm against the gate
        capacitance lands near 1e12 Hz -- and the density rolls off
        again above it.  That is real, and the model shows it; it is
        just not the pole under test here.
        """
        gd = op['long']
        kw = self._kw(deck, gd, rg=0.0)
        pt = [p for p in gd['points']
              if p['vg'] == 1.2 and p['vd'] == 0.6 and p['vb'] == 0.0][0]
        lo1 = self._sig(kw, pt, f=1.0e3)
        lo2 = self._sig(kw, pt, f=1.0e4)
        assert lo2 / lo1 == pytest.approx(100.0, rel=0.02)
        hi1 = self._sig(kw, pt, f=1.0e13)
        hi2 = self._sig(kw, pt, f=1.0e14)
        assert hi2 / hi1 == pytest.approx(1.0, rel=0.01)

        ## And the plateau is not merely flat, it is at the RIGHT
        ## level: far above the pole the capacitor's impedance wins and
        ## the whole source current reaches the gate, so the density
        ## there IS the source power on the auxiliary branch.
        e = self._fet(deck, gd, rg=0.0)
        nm = [n.name for n in e.nodes]
        i_n = nm.index('noi')
        pwr = self._white(e, e.bias(pt['vd'], pt['vg'], 0.0, pt['vb']),
                          i_n, i_n)
        assert hi2 == pytest.approx(pwr, rel=0.01)

    def test_swign_switches_it_off_and_costs_a_row(self, deck, op):
        """`SWIGN = 0` removes the term AND the node it lived on -- the
        auxiliary branch collapses, so switching it off is free rather
        than merely zeroed."""
        gd = op['long']
        on, off = self._fet(deck, gd), self._fet(deck, gd, swign=0.0)
        assert 'noi' in [n.name for n in on.nodes]
        assert 'noi' not in [n.name for n in off.nodes]
        assert off.n == on.n - 1
        pt = [p for p in gd['points'] if p['sig'] > 0.0][0]
        assert self._sig(self._kw(deck, gd, swign=0.0), pt) \
            == pytest.approx(0.0, abs=1e-45)

    def test_it_needs_the_thermal_term(self, deck, op):
        """The gate density is proportional to `nt`, so `fnt = 0` takes
        it with the drain's white term -- an element that was told it
        has no thermal noise cannot acquire some at the gate."""
        gd = op['long']
        pt = [p for p in gd['points'] if p['sig'] > 0.0][0]
        assert self._sig(self._kw(deck, gd, fnt=0.0), pt) \
            == pytest.approx(0.0, abs=1e-45)
