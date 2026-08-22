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

import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.compact import PspMosLongChannel
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
    card = deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1, m=1,
                             pre_layout=1)
    kw = psp_scaling.to_long_channel(card, w=w, l=l)
    e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                          cm.Node('b'), **kw)
    e.update_iparv()
    v = np.asarray(sweep['v'], float)
    r = np.abs(np.asarray(sweep['i_d'], float))
    b = sweep['bias']
    if sweep['sweep'] == 'Vd':
        got = np.array([np.asarray(e.i(np.array([x, b['Vg'], b['Vs'],
                                                 b['Vb']])), float)[0]
                        for x in v])
    else:
        got = np.array([np.asarray(e.i(np.array([b['Vd'], x, b['Vs'],
                                                 b['Vb']])), float)[0]
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

        with_qm = psp_scaling.to_long_channel(card, w=w, l=l)
        without = psp_scaling.to_long_channel(dict(card, qmc=0.0), w=w, l=l)
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
        kw = psp_scaling.to_long_channel(card, w=w, l=l)
        assert kw['rs'] > 0.0, 'the card does specify one'

        cm.default_toolkit = numeric
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, rs=0.0))
        on.update_iparv()
        off.update_iparv()
        for vd in (0.05, 0.4, 0.9, 1.5):
            x = np.array([vd, 1.2, 0.0, 0.0])
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
        kw = psp_scaling.to_long_channel(card, w=w, l=l)
        assert kw['ct'] > 0.05, 'the card does specify one'

        cm.default_toolkit = numeric
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, ct=0.0))
        on.update_iparv()
        off.update_iparv()
        for vg in (0.6, 1.0, 1.5):
            x = np.array([0.05, vg, 0.0, 0.0])
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
        kw = psp_scaling.to_long_channel(card, w=w, l=l)
        assert kw['alp'] > 0.0

        cm.default_toolkit = numeric
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, alp=0.0))
        on.update_iparv()
        off.update_iparv()
        deep = np.array([1.5, 1.2, 0.0, 0.0])
        lin = np.array([0.02, 1.2, 0.0, 0.0])
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

    def test_the_short_device_is_still_the_worse_one(self, deck, ref):
        """It is, but no longer by much.

        This test used to assert the short device was at least 3x worse
        than the long one.  It stopped being true when the `GPE`
        width-adjustment bug was fixed -- the short device went from
        ~12% to ~2% and the ratio fell to about 2x.  Rewritten rather
        than relaxed, because the thing worth pinning changed: the point
        is no longer "the omissions are large" but "they are bounded".
        """
        _, r_l, g_l, _ = _compare(deck, ref['nmos_long_idvd'])
        _, r_s, g_s, _ = _compare(deck, ref['nmos_idvg_vd0p05'])
        long_err = np.median(np.abs(g_l / r_l - 1.0))
        short_err = np.median(np.abs(g_s / r_s - 1.0))
        assert short_err > long_err
        assert short_err < 0.05, 'short device %.3f' % short_err


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
        kw = psp_scaling.to_long_channel(card, w=w, l=l)
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
        """The comparable quantity, and what it says.

        Measured: PSP shifts +237 mV from L = 1 um to L = 0.13 um and we
        shift +214 mV -- 90% of it, from the pocket-implant doping
        formula plus the length terms in `VFB` and `DPHIB`.  The missing
        23 mV is worth about 3% of drain current at Vg = 1.2, which is a
        quarter of the ~12% the short device is actually off by.

        So the short-device residual is NOT predominantly a threshold
        error, and DIBL is not the next layer -- `CF` scales to 1e-7 on
        the long device and contributes under a millivolt at Vd = 0.05,
        where two of the worst sweeps sit.
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
        kw = psp_scaling.to_long_channel(card, w=e['w'], l=e['l'])
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
        kw = psp_scaling.to_long_channel(card, w=e['w'], l=e['l'])
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
                                         l=e['l'])
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
