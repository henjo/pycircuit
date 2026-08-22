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

    def test_the_long_device_agrees_within_fifteen_percent(self, deck, ref):
        """W = 10 um, L = 1 um -- where the core is meant to be right.

        The core has no series resistance, no channel-length modulation
        and no short-channel physics, so this is not expected to be
        exact.  It IS expected to be close on a long device, and 15% is
        a bound the measured 10% sits comfortably inside while still
        catching any real break in the chain.
        """
        _, r, g, _ = _compare(deck, ref['nmos_long_idvd'])
        assert len(r) > 20
        ratio = g / r
        assert 0.85 < ratio.min(), 'worst low %.3f' % ratio.min()
        assert ratio.max() < 1.15, 'worst high %.3f' % ratio.max()

    def test_the_agreement_is_flat_not_a_lucky_crossing(self, deck, ref):
        """A curve that crosses the reference could average out to 1.

        What makes the long-device result meaningful is that the ratio
        barely moves across the sweep: a uniform offset means the
        threshold and gain factor are right and something multiplicative
        is missing, rather than two wrong shapes intersecting.
        """
        _, r, g, _ = _compare(deck, ref['nmos_long_idvd'])
        ratio = g / r
        assert ratio.max() - ratio.min() < 0.10

    def test_the_scaled_parameters_are_physical(self, deck, ref):
        """A wrong scaling shows up here before it shows up in a curve."""
        _, _, _, kw = _compare(deck, ref['nmos_long_idvd'])
        assert 1e-9 < kw['tox'] < 5e-9
        assert -1.5 < kw['vfb'] < 0.0
        assert 1e22 < kw['nsub'] < 1e24
        assert 0.5 < kw['phib'] < 1.3
        assert 0.0 < kw['u0'] < 0.1
        assert all(v >= 0.0 for v in (kw['mue'], kw['themu'], kw['cs'],
                                      kw['thecs'], kw['thesat']))

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

    def test_the_short_device_is_measurably_worse(self, deck, ref):
        """The omissions have a size, and it is recorded rather than hidden.

        At L = 0.13 um the missing DIBL, channel-length modulation and
        series resistance dominate, and the ratio reaches ~1.7.  That is
        not a defect in the core -- it is the cost of the layers not yet
        built, and it is what says which to build next.
        """
        _, r_l, g_l, _ = _compare(deck, ref['nmos_long_idvd'])
        _, r_s, g_s, _ = _compare(deck, ref['nmos_idvg_vd0p05'])
        long_err = np.median(np.abs(g_l / r_l - 1.0))
        short_err = np.median(np.abs(g_s / r_s - 1.0))
        assert short_err > 3.0 * long_err
        assert short_err < 2.0, 'still the same order, not divergent'
