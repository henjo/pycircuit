"""Shooting tests: noise devices.  Split out of test_analysis_shooting.py on
2026-09-27 (in its original order); shared helpers are in
`_shooting_fixtures.py`, the HDL elements in `_shooting_elements.py`.
"""
from pycircuit.circuit import *
from pycircuit.circuit.shooting import (PAC, algebraic_conditioning,
                                        topological_index)
import warnings
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter as _HdlParameter, white_noise)
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest
import functools as _functools


def test_the_compact_mos_models_really_are_transcapacitive():
    """The premise of the test above, pinned so it cannot go stale quietly.

    If a future change symmetrised these models' `C`, the blindness note
    would become harmless and nobody would know to remove it — and, worse,
    a real physical effect would have been lost. Ward-Dutton partition in
    saturation decouples the drain from the channel charge, so `∂q_g/∂v_d`
    is small while `∂q_d/∂v_g` is not. That asymmetry is the model being
    right, not a defect.
    """
    from pycircuit.circuit import compact
    circuit.default_toolkit = circuit.numeric
    worst = 0.0
    for name in ('PspMosLongChannel', 'PspPmosLongChannel'):
        cls = getattr(compact, name)
        inst = cls(*cls.terminals)
        for vg in np.linspace(0.0, 1.2, 7):
            for vd in np.linspace(0.0, 1.2, 7):
                Cm = np.asarray(inst.C(np.array([vd, vg, 0.0, 0.0])),
                                dtype=float)
                s = float(np.max(np.abs(Cm)))
                if s == 0.0:
                    continue
                worst = max(worst, float(np.max(np.abs(Cm - Cm.T))) / s)
    assert worst > 0.1, \
        'max|C-C^T|/max|C| = %.3e; the compact MOS models are no longer ' \
        'transcapacitive, so the transpose blind spot recorded against ' \
        'them needs re-deriving rather than deleting' % worst

    ## ⚠ AND `Cox` ITSELF, ANCHORED BY GEOMETRY RATHER THAN BY `C()`.
    ## The comparable form of McAndrew's bound is `|C_ij - C_ji| / Cox`,
    ## so `Cox` has to come from somewhere the code under test cannot
    ## move it: `eps_ox W L / tox` from the model's own parameters.
    inst = compact.PspMosLongChannel(*compact.PspMosLongChannel.terminals)
    p = inst.iparv
    cox_geom = 3.9 * 8.8541878128e-12 * p.w * p.l / p.tox
    cgg = max(abs(float(np.asarray(inst.C(np.array([0.0, vg, 0.0, 0.0])),
                                   dtype=float)[1, 1]))
              for vg in np.linspace(0.0, 2.5, 26))
    assert abs(cgg / cox_geom - 1.0) < 0.05, \
        'max Cgg = %.4e against eps_ox W L / tox = %.4e; if these have ' \
        'parted company the Cox normalisation below is no longer ' \
        'anchored by geometry' % (cgg, cox_geom)
    nonrecip = 0.0
    for vg in np.linspace(0.0, 1.2, 13):
        for vd in np.linspace(0.0, 1.2, 13):
            Cm = np.asarray(inst.C(np.array([vd, vg, 0.0, 0.0])),
                            dtype=float)
            nonrecip = max(nonrecip, float(np.max(np.abs(Cm - Cm.T))))
    assert 25.0 < nonrecip / cgg / 0.01 < 40.0, \
        'nonreciprocity is %.1fx McAndrew 1%% of Cox; the recorded 32x ' \
        'no longer describes this model' % (nonrecip / cgg / 0.01)

    ## ⚠ AND THE SAME QUANTITY AT McANDREW'S OWN CONDITION, `Vds = 0`,
    ## which is where his bound is stated and where it HOLDS. Pinned
    ## separately from the saturation figure on purpose: they are two
    ## different claims and conflating them is what produced a wrong
    ## reading of the paper.
    at_vds0 = 0.0
    for vg in np.linspace(0.0, 1.2, 13):
        Cm = np.asarray(inst.C(np.array([0.0, vg, 0.0, 0.0])), dtype=float)
        at_vds0 = max(at_vds0, float(np.max(np.abs(Cm - Cm.T))))
    assert at_vds0 / cgg < 0.01, \
        'at Vds = 0 the nonreciprocity is %.4f of Cox, above McAndrew 1%% ' \
        'where he states it. That would be the STRONGER claim -- his ' \
        'bound failing for a real compact model -- and it needs saying ' \
        'so, not silently absorbing into the saturation figure'\
        % (at_vds0 / cgg)
    assert nonrecip / at_vds0 > 20.0, \
        'the nonreciprocity is no longer created by Vds (ratio %.1f); the ' \
        'transpose gate depends on the fixture SWINGING, not on its DC ' \
        'bias' % (nonrecip / at_vds0)


def test_the_compact_mos_noise_is_off_without_a_card_not_absent():
    """⚠⚠ THIS TEST USED TO ASSERT "no noise model at all". THAT WAS WRONG.

    `PspMosLongChannel` has channel thermal and flicker noise, declared at
    `compact.py:834` as `white_noise(mult·n_sid) + flicker_noise(mult·n_sfl,
    ef)`. What is zero is the COEFFICIENTS: `fnt = 0`, `nfa = 0` by
    default, and `compact.py:1049` says why — *"`fnt = 0` switches the
    thermal term off … an element built without a card is noiseless."*

    ⚠ SO THE EARLIER MEASUREMENT WAS OF A DEFAULT-CONSTRUCTED INSTANCE AND
    THE CLAIM WAS ABOUT THE MODEL. Every sweep read `CY = 0` because the
    element had no card, and the conclusion drawn was that the feature did
    not exist. §D shape 0c: the fixture could not express the thing under
    test, and the fixture was a constructor call.

    ⚠ THE MODEL'S OWN DOCSTRING SAID SO — it lists "channel thermal and
    flicker noise" under *"Since built, and no longer absent"*, and warns
    two paragraphs later that *"a stale gap note is worse than none: it is
    trusted like a measurement and it is not one."* The note was current;
    the reader was not.

    With `fnt = 1` the element grows a noise branch (`n` goes 4 → 5) and
    `CY` is nonzero, white, and state-dependent.
    """
    from pycircuit.circuit import compact
    circuit.default_toolkit = circuit.numeric
    cls = compact.PspMosLongChannel

    bare = cls(*cls.terminals)
    assert bare.iparv.fnt == 0.0 and bare.iparv.nfa == 0.0, \
        'the noise coefficients are no longer zero by default, so ' \
        '"an element built without a card is noiseless" has changed'
    assert bare.n == 4
    for vg in (0.4, 1.2):
        cy = np.asarray(bare.CY(np.array([0.0, vg, 0.0, 0.0]),
                                2.0 * np.pi * 1e6), dtype=float)
        assert float(np.max(np.abs(cy))) == 0.0

    noisy = cls(*cls.terminals, fnt=1.0)
    assert noisy.n == 5, \
        'enabling fnt no longer adds the noise branch (n = %d)' % noisy.n
    vals = []
    for vd in (0.0, 0.6, 1.2):
        x = np.zeros(noisy.n)
        x[0], x[1] = vd, 1.0
        vals.append(float(np.asarray(
            noisy.CY(x, 2.0 * np.pi * 1e6), dtype=float)[0, 0]))
    assert min(vals) > 0.0, 'CY is still zero with fnt = 1'
    assert max(vals) / min(vals) > 1.2, \
        'CY no longer depends on the bias (%s); a state-independent MOS ' \
        'noise model would be the physically wrong one, and would also ' \
        'make pnoise(modulated=True) unnecessary' % vals


def test_the_mos_thermal_noise_satisfies_the_fluctuation_dissipation_theorem():
    """⚠⚠ THE EXTERNAL ANCHOR FOR A COMPACT MODEL'S NOISE — thermodynamics.

    At `Vds = 0` a MOSFET is in thermal equilibrium: it dissipates nothing
    and the fluctuation-dissipation theorem fixes its current noise
    completely,

        S_id  =  4 k T g_ds ,    g_ds = ∂I_d/∂V_d

    with **no model freedom whatever**. Any noise model that misses this
    is wrong regardless of what it does elsewhere, and one that hits it has
    its absolute scale anchored to thermodynamics rather than to a fit.
    This is the MOS analogue of `kT/C`.

    MEASURED at `fnt = 1`, `Vds = 0`, `T = 300 K`:

        Vg     g_ds (S)       CY[0,0]        4kT·g_ds       ratio
        0.40   5.220648e-05   8.894646e-25   8.649459e-25   1.028347
        0.80   2.180318e-04   3.690935e-24   3.612305e-24   1.021767
        1.20   3.392372e-04   5.711105e-24   5.620411e-24   1.016137
        1.50   3.918655e-04   6.579492e-24   6.492344e-24   1.013423

    ⚠ SATISFIED TO 1.3–2.8%, AND THE RESIDUAL IS STRUCTURAL RATHER THAN
    SCATTER: it has a sign and falls monotonically with `Vg`. That is a
    property of PSP's channel-integrated `n_sid` against the ideal
    `4kT·g_ds`, not a normalisation to tune.

    ⚠⚠ AND IT IS DELIBERATELY NOT TUNED. `fnt` is an exact linear scale on
    the thermal PSD — measured 0.509293 / 1.018586 / 2.037172 at
    `fnt` = 0.5/1/2 — so setting `fnt = 0.98175` would make this test read
    1.000000 and would be **fitting a physical constant to a discrepancy
    we do not understand**. The tolerance is 5%, chosen to admit the
    measured residual and to catch a factor.

    Also asserted: the thermal term is exactly WHITE (identical at 1e3,
    1e6, 1e9 Hz with the flicker coefficient off), which is what makes
    `4kT·g_ds` the right comparison at any frequency.
    """
    from pycircuit.circuit import compact
    circuit.default_toolkit = circuit.numeric
    K_B = 1.380649e-23
    T = 300.0
    inst = compact.PspMosLongChannel(
        *compact.PspMosLongChannel.terminals, fnt=1.0)
    n = inst.n

    ## white: no frequency dependence with the flicker coefficient off
    x = np.zeros(n)
    x[1] = 1.0
    vals = [float(np.asarray(inst.CY(x, 2.0 * np.pi * f), dtype=float)[0, 0])
            for f in (1e3, 1e6, 1e9)]
    assert max(vals) == min(vals), \
        'the thermal term is not white across 1e3..1e9 Hz (%s)' % vals

    ratios = []
    for vg in (0.4, 0.8, 1.2, 1.5):
        x = np.zeros(n)
        x[0], x[1] = 0.0, vg
        gds = float(np.asarray(inst.G(x), dtype=float)[0, 0])
        sid = float(np.asarray(inst.CY(x, 2.0 * np.pi * 1e6),
                               dtype=float)[0, 0])
        assert gds > 0.0
        ratios.append(sid / (4.0 * K_B * T * gds))
    for vg, r in zip((0.4, 0.8, 1.2, 1.5), ratios):
        assert abs(r - 1.0) < 0.05, \
            'at Vg = %.2f, Vds = 0 the model gives S_id = %.4f x 4kT g_ds. ' \
            'In equilibrium that ratio is fixed by thermodynamics; a ' \
            'deviation this size is a defect in the noise model, not a ' \
            'modelling choice' % (vg, r)
    ## the residual is structural: monotone in Vg, not scatter
    assert ratios == sorted(ratios, reverse=True), \
        'the FDT residual %s is no longer monotone in Vg, so it has ' \
        'stopped being the systematic effect this test documents' % ratios


## ---------------------------------------------------------------------------
## COLOURED SOURCES -- an element that has colour, and the fold that reads it
## ---------------------------------------------------------------------------

def test_the_IS_colour_is_the_named_shape():
    """`noiseTau` is white noise through an RC, `noiseFc` a flicker corner.

    Checked at the element, against the closed forms, so that every gate
    below tests the FOLD and not the source's algebra.  The Lorentzian is
    exactly realisable in-netlist and that realisation is the reference
    for the folds; flicker is not (Demir 1996: one state per decade), and
    it returns the white value at `w = 0` rather than infinity.
    """
    circuit.default_toolkit = circuit.numeric
    P, tau, fc = 1e-6, 0.3, 50.0
    lor = IS('a', gnd, i=0.0, noisePSD=P, noiseTau=tau)
    fl = IS('a', gnd, i=0.0, noisePSD=P, noiseFc=fc)
    white = IS('a', gnd, i=0.0, noisePSD=P)
    x = np.zeros(2)
    for w in (0.0, 1.0, 10.0, 1e3):
        assert np.allclose(white.CY(x, w), P * np.array([[1, -1], [-1, 1]]))
        assert np.allclose(lor.CY(x, w)[0, 0], P / (1.0 + (w * tau) ** 2),
                           rtol=1e-14)
        assert np.allclose(lor.CY(x, -w), lor.CY(x, w)), 'colour is even in w'
    assert np.allclose(fl.CY(x, 0.0)[0, 0], P), 'white at DC, not infinite'
    for w in (1.0, 10.0, 1e3):
        assert np.allclose(fl.CY(x, w)[0, 0], P * (1.0 + 2.0 * np.pi * fc / w),
                           rtol=1e-14)


def test_every_library_device_states_the_sign_of_its_flicker_current():
    """A 1/f current is a slow relative fluctuation of a conductance TIMES the
    current through it, so it follows the current's sign -- and a periodic
    noise fold needs that sign (see
    `test_a_coloured_source_keeps_the_sign_of_its_scale_factor...`).  The
    SPICE-style models wrote `flicker_noise(kf |I|^af)`, which states none and
    left them on the |m| fold; they now write `(I/|I|) * flicker_noise(...)`.

    Per device, at a bias and at its mirror: the amplitudes REBUILD the
    flicker part of `CY` (which is unchanged: the factor is +-1 to
    (1e-30/I)^2), and they change sign with the current.
    """
    import pycircuit.circuit.elements_hdl as eh
    w1, winf = 2 * np.pi * 10.0, 2 * np.pi * 1e30
    mos = [(0.5, 1.5, 0.0, 0.0), (0.0, 1.5, 0.5, 0.0)]
    cases = [('DiodeSpiceHdl', dict(kf=1e-12), [(0.7, 0.0), (-1.0, 0.0)]),
             ('GummelPoonNpnHdl', dict(kf=1e-12),
              [(1.0, 0.7, 0.0), (0.0, -0.7, 0.0)]),
             ('EkvNmosHdl', dict(kf=1e-24), mos),
             ('MosLevel1Hdl', dict(kf=1e-24), mos),
             ('MosLevel3Hdl', dict(kf=1e-24), mos),
             ('MesfetStatzHdl', dict(kf=1e-12),
              [(0.5, 0.0, 0.0), (0.0, 0.0, 0.5)])]
    for name, kw, biases in cases:
        cls = getattr(eh, name)
        with quiet():
            el = cls(*['n%d' % i for i in range(len(cls.terminals))], **kw)
            signs = []
            for bias in biases:
                x = np.zeros(el.n)
                x[:len(bias)] = bias
                W = el.noise_amplitudes(x, w1)
                assert W is not None and W.shape == (el.n, 1), name
                fl = (np.asarray(el.CY(x, w1), dtype=complex)
                      - np.asarray(el.CY(x, winf), dtype=complex))
                assert np.max(np.abs(fl)) > 0, name        # the flicker is alive
                err = np.max(np.abs(W @ W.conj().T - fl)) / np.max(np.abs(fl))
                assert err < 1e-9, (name, err)
                signs.append(np.sign(np.real(W[int(np.argmax(np.abs(W[:, 0]))), 0])))
        assert signs[0] * signs[1] == -1.0, (name, signs)


def test_the_level3_channel_noise_is_shot_noise_in_weak_inversion_and_thermal_above():
    """MosLevel3's channel noise across inversion, against the anchors that
    leave no freedom (2026-10-09; reported from a Si2302CDS card, VTO 1.17 V
    with `nfs`, where the drain noise read 32 uV/rtHz for 1 kOhm's 5.8 nV):

      * weak inversion, saturated: shot noise of the diffusion current,
        `S_id = 2 q Id` (it was `(vdsn^2/3)/1e-9`: 1e8..1e10 x that);
      * `vds -> 0`: equilibrium, `S_id = 4 k T g_ds`, weak and strong;
      * strong inversion, saturated: `(8/3) k T gm` (long channel);
      * between: monotone from the weak value (3 xn / 4 of it) to 1 --
        no spike, and nothing below the strong value (the old form
        dipped to 0.84 at Vgs 1.5)."""
    from pycircuit.circuit import elements_hdl as eh
    circuit.default_toolkit = circuit.numeric
    KB, Q, T = 1.380649e-23, 1.602176634e-19, 300.15
    m = eh.MosLevel3Hdl('d', 'g', 's', 'b', vto=1.17, kp=2e-5, w=0.273, l=2e-6,
                        nfs=1e12, tox=5e-8, nsub=1e17, gamma=0.5, phi=0.7,
                        theta=0.05, tnom=27.0)

    def at(vd, vg):
        x = np.zeros(m.n)
        x[0], x[1] = vd, vg
        G = np.asarray(m.G(x), dtype=float)
        return (float(np.asarray(m.i(x), dtype=float)[0]), G[0, 0], G[0, 1],
                float(np.asarray(m.CY(x, 2 * np.pi * 1e6), dtype=float)[0, 0]))

    for vg in (0.6, 0.9, 1.2):
        i, _gds, _gm, s = at(5.0, vg)
        assert abs(s / (2 * Q * i) - 1.0) < 1e-3, (vg, s / (2 * Q * i))
        i, gds, _gm, s = at(1e-3, vg)
        assert abs(s / (4 * KB * T * gds) - 1.0) < 0.03, (vg, s / (4 * KB * T * gds))
    i, gds, _gm, s = at(1e-3, 3.0)
    assert abs(s / (4 * KB * T * gds) - 1.0) < 0.06, s / (4 * KB * T * gds)
    ratios = []
    for vg in (1.2, 1.3, 1.4, 1.5, 1.6, 1.8):
        i, _gds, gm, s = at(5.0, vg)
        ratios.append(s / (8.0 / 3.0 * KB * T * gm))
    assert ratios == sorted(ratios, reverse=True), ratios
    assert min(ratios) > 0.97 and abs(ratios[-1] - 1.0) < 0.03, ratios
    assert ratios[0] < 4.0, ratios


@pytest.mark.parametrize('cls', ['MosLevel1Hdl', 'MosLevel3Hdl'])
def test_the_mos_flicker_forms_are_ngspice_s_nlev(cls):
    """`nlev` (2026-10-09), ngspice's flicker forms, each against its
    formula at a saturated bias, `cox` per area and `Leff = l - 2 ld`:

      0  kf |Id|^af / (f cox Leff^2)          SPICE2's; the default
      1  kf |Id|^af / (f cox W Leff)          the gate area, as a
                                              commercial simulator has it
      2  kf gm^2 / (f^af cox W Leff)          ngspice-47's default (3 too)

    `gm` is the device's own `dId/dVg` (its Jacobian, drain, source and
    bulk held); `af = 1.2` so the gm form's `f^af` is tested at two
    frequencies.  An AD8606 input pair read 5e4 x apart between forms 0
    and 2 -- `W Vov^2 / (4 L Id)` for a square law."""
    from pycircuit.circuit import elements_hdl as eh
    circuit.default_toolkit = circuit.numeric
    card = dict(vto=0.7, kp=1e-4, w=20e-6, l=2e-6, ld=0.1e-6, tox=2e-8, kf=1e-25, af=1.2)
    cox = 3.9 * 8.854187817e-12 / card['tox']
    leff, w = card['l'] - 2 * card['ld'], card['w']
    x = np.zeros(4)
    x[0], x[1] = 2.0, 1.5

    def flicker(nlev, f):
        m = getattr(eh, cls)('d', 'g', 's', 'b', nlev=nlev, **card)
        m0 = getattr(eh, cls)('d', 'g', 's', 'b', nlev=nlev, **dict(card, kf=0.0))
        s = (float(np.asarray(m.CY(x, 2 * np.pi * f), dtype=float)[0, 0])
             - float(np.asarray(m0.CY(x, 2 * np.pi * f), dtype=float)[0, 0]))
        G = np.asarray(m.G(x), dtype=float)
        return s, abs(float(np.asarray(m.i(x), dtype=float)[0])), G[0, 1]

    kf, af = card['kf'], card['af']
    for nlev, want in ((0, lambda i, gm, f: kf * i ** af / (f * cox * leff ** 2)),
                       (1, lambda i, gm, f: kf * i ** af / (f * cox * w * leff)),
                       (2, lambda i, gm, f: kf * gm ** 2 / (f ** af * cox * w * leff)),
                       (3, lambda i, gm, f: kf * gm ** 2 / (f ** af * cox * w * leff))):
        for f in (10.0, 1e3):
            s, i, gm = flicker(nlev, f)
            assert gm > 0 and i > 0
            r = s / want(i, gm, f)
            assert abs(r - 1.0) < 1e-6, (cls, nlev, f, r)
