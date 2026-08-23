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
        kw = psp_scaling.to_long_channel(card, w=1e-6, l=0.13e-6)
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
        kw = psp_scaling.to_long_channel(card, w=10e-6, l=1e-6)
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, kp=0.0))
        on.update_iparv()
        off.update_iparv()
        for vg in (0.6, 1.0, 1.4, 1.8):
            x = np.array([0.05, vg, 0.0, 0.0])
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
        kw = psp_scaling.to_long_channel(card, w=10e-6, l=1e-6)
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        assert kw['rs'] > 0 and kw['alp'] > 0 and kw['kp'] > 0, \
            'the card must actually switch these on for this to test '\
            'anything'
        for vg in (0.6, 1.2, 1.8):
            for vd in (0.05, 0.3, 0.9, 1.5):
                f = np.asarray(e.i(np.array([vd, vg, 0.0, 0.0])), float)
                r = np.asarray(e.i(np.array([0.0, vg, vd, 0.0])), float)
                assert f[0] == pytest.approx(-r[0], rel=1e-14), (vg, vd)
            q = np.asarray(e.q(np.array([0.5, vg, 0.0, 0.0])), float)
            assert abs(q.sum()) < 1e-16 * max(np.abs(q).max(), 1e-30)


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
                              ng=1, m=1, pre_layout=1), w=1e-6, l=0.13e-6)
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
                              ng=1, m=1, pre_layout=1), w=1e-6, l=0.13e-6)
        on = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                               cm.Node('b'), **kw)
        off = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                                cm.Node('b'), **dict(kw, alp1=0.0,
                                                     alp2=0.0))
        on.update_iparv()
        off.update_iparv()
        lift = {}
        for vg in (0.6, 1.8):
            x = np.array([1.4, vg, 0.0, 0.0])
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
                w=s['w'], l=s['l'])
            for k, over in enumerate((dict(alp1=0.0, alp2=0.0), {})):
                e = PspMosLongChannel(cm.Node('d'), cm.Node('g'),
                                      cm.Node('s'), cm.Node('b'),
                                      **dict(kw, **over))
                e.update_iparv()
                v = np.asarray(s['v'], float)
                r = np.abs(np.asarray(s['i_d'], float))
                b = s['bias']
                if s['sweep'] == 'Vd':
                    g = np.array([np.asarray(e.i(np.array(
                        [x, b['Vg'], b['Vs'], b['Vb']])), float)[0]
                        for x in v])
                else:
                    g = np.array([np.asarray(e.i(np.array(
                        [b['Vd'], x, b['Vs'], b['Vb']])), float)[0]
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
            w=g['w'], l=g['l'])
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

    def test_every_sweep_is_within_one_percent(self, deck, ref):
        """The state of the whole comparison, as a regression guard.

        Nothing here is fitted: every parameter comes off the card
        through the scaling layer.  The history this guard replaces is
        worth keeping in view -- the worst sweep was 11% out before the
        saturation voltage, 7% after it, and 0.3% once the saturation
        voltage was fed back into the channel-length modulation, which
        is what it had been missing all along.

        Held at 1% rather than at the measured 0.4%: this is a guard
        against regression, not a pin on numbers that should be free to
        improve.
        """
        cm.default_toolkit = numeric
        worst = {}
        for name in ('nmos_long_idvd', 'nmos_idvd_vg1p2',
                     'nmos_idvd_vg0p6', 'nmos_idvg_vd0p05',
                     'nmos_idvg_vd1p2', 'nmos_idvg_vb_m1'):
            _, r, g, _ = _compare(deck, ref[name])
            worst[name] = abs(np.median(g / r) - 1.0)
        assert max(worst.values()) < 0.01, worst
        assert sum(worst.values()) < 0.04, worst
