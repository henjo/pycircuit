# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The second Phase-6 library batch: the four models chosen because each
one is the first production user of a DSL capability nothing consumed.

===================  ==========================================
model                capability it is the first real user of
===================  ==========================================
`XtalHdl`            (none -- the warm-up: a passive resonator
                     the repo could not previously build)
`FerriteBeadHdl`     ``laplace_nd`` from network coefficients
`RSkinHdl`           ``laplace_zp``
`ComparatorHdl`      ``Cross``
`MemristorHdl`       ``idt`` as a real DAE state
===================  ==========================================

Nothing here compares a model against itself.  Every number is one of

* a **closed form the model does not contain** -- the BVD impedance
  ``1/(1/(rm + s*lm + 1/(s*cm)) + s*c0)``, the parallel-RLC bead
  impedance, ``rdc*(1 + sqrt(s/ws))``, the memristor's high-frequency
  degeneration to ``V/R(ic)``, and the comparator's fold voltage
  ``vhyst*qf - atanh(qf)/gain``, which appears in the model only as a
  breakpoint hint and never in an answer;
* an **independent numerical integration** -- RK4 in numpy, from the
  Strukov/Biolek equations, for the memristor state;
* an **independent evaluation of the same rational function** -- the
  pole-zero product taken straight from the root list in numpy, against
  what ``laplace_zp``'s state-space realisation actually solves.  That
  is the check on the DSL machinery rather than on the physics;
* a **structural identity** that must hold whatever the equations are:
  the bead's ``|Z| = rdc + rp`` exactly at self-resonance, the skin
  ladder's dc floor, ``|i| <= |v|/ron`` everywhere for the memristor
  (which is what "pinched at the origin" means), the comparator's
  symmetric loop;
* a **finite difference**, for every Jacobian claim, through
  `hdl.check_jacobians`.

and each is exercised at the seams: ``rm = 0``, ``vhyst = 0``,
``gain*vhyst < 1`` (no fold at all), the memristor state driven past
both boundaries and past the ``x`` at which ``R(x)`` would be zero.
"""

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pycircuit.circuit.circuit
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.analysis_ss import AC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.elements import SubCircuit, VS, VSin, R, IS as ISource
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.hdl import Node, check_jacobians, x_layout

TWOPI = 2 * np.pi


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    """Every element in this file is built by calling its class, which
    reads ``circuit.default_toolkit`` at construction.  Set once here
    rather than at fifty call sites (`validation-design`: fix a
    convention with a fixture, not with discipline)."""
    old = pycircuit.circuit.circuit.default_toolkit
    pycircuit.circuit.circuit.default_toolkit = numeric
    yield
    pycircuit.circuit.circuit.default_toolkit = old


def _mk(cls, *nodes, **kw):
    el = cls(*[Node(nm) for nm in nodes], **kw)
    el.update_iparv()
    return el


def _quiet_transient(cir, **kw):
    tran = Transient(cir, toolkit=numeric, **{'uic': True, **kw})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(**kw.pop('solve', {}) or {})
    return tran, res


## ======================================================================
## XtalHdl -- the Butterworth-Van Dyke resonator
## ======================================================================

#: A real 32 MHz AT-cut card, not the defaults.
XTAL = dict(lm=8.0e-3, cm=3.2e-15, rm=12.0, c0=3.5e-12)


def _z_bvd(f, lm, cm, rm, c0):
    """The BVD impedance, transcribed from the network and not from
    `XtalHdl`: the motional arm in series, the holder capacitance
    across the pair."""
    s = 1j * TWOPI * np.asarray(f, dtype=complex)
    zm = rm + s * lm + 1.0 / (s * cm)
    return 1.0 / (1.0 / zm + s * c0)


def _xtal_z(card, freqs):
    """|Z| of one crystal, measured by driving it with a 1 V AC source
    and reading the source current."""
    c = SubCircuit()
    n1 = c.add_node('n1')
    c['vs'] = VS(n1, gnd, vac=1.0)
    c['x'] = eh.XtalHdl(n1, gnd, **card)
    res = AC(c).solve(freqs=np.atleast_1d(freqs))
    return -1.0 / np.asarray(res.i('vs.plus').y)


def test_crystal_matches_the_bvd_closed_form():
    """Impedance and phase against the network's own closed form, over
    eight decades and across both resonances."""
    f = np.concatenate([np.logspace(2, 9, 60),
                        np.linspace(31.4e6, 31.5e6, 40)])
    got = _xtal_z(XTAL, f)
    want = _z_bvd(f, **XTAL)
    assert_allclose(np.abs(got), np.abs(want), rtol=1e-9)
    assert_allclose(np.angle(got), np.angle(want), rtol=0, atol=1e-9)


def test_crystal_series_and_parallel_resonance():
    """``fs = 1/(2*pi*sqrt(lm*cm))`` and ``fp = fs*sqrt(1 + cm/c0)``,
    located by where the measured impedance is smallest and largest.

    The predicted frequencies are computed from the card, never from the
    sweep; the sweep only has to LAND on them.  ``|Z(fs)| = rm`` to
    within the holder capacitance's own shunting, which at 31 MHz is
    ``|1/(j*w*c0)| = 1.5 kohm`` against 12 ohm -- 0.8% and accounted
    for, not tolerated.
    """
    fs = 1.0 / (TWOPI * np.sqrt(XTAL['lm'] * XTAL['cm']))
    fp = fs * np.sqrt(1 + XTAL['cm'] / XTAL['c0'])
    ## 31.4557 MHz and 31.4701 MHz -- 14 kHz apart on 31 MHz, which is
    ## the 460 ppm pulling range of a 1100:1 capacitance ratio.
    assert 4.5e-4 < (fp - fs) / fs < 4.7e-4

    f = np.linspace(fs * (1 - 2e-4), fp * (1 + 2e-4), 4001)
    z = np.abs(_xtal_z(XTAL, f))
    assert_allclose(f[np.argmin(z)], fs, rtol=1e-5)
    assert_allclose(f[np.argmax(z)], fp, rtol=1e-5)
    ## The series minimum IS the ESR, shunted by c0: the motional arm is
    ## `rm` there and c0 is 1447 ohm at 31 MHz, so the parallel
    ## combination is 12.000 ohm to five figures.
    zc0 = 1.0 / (TWOPI * fs * XTAL['c0'])
    assert_allclose(z.min(), abs(1 / (1 / XTAL['rm'] + 1j / zc0)), rtol=1e-4)

    ## The parallel maximum.  Derived rather than quoted, because the
    ## formula that gets quoted -- lm/(rm*(c0 + cm)) -- is wrong by the
    ## capacitance ratio, 1094x here, and it looks entirely plausible.
    ## At antiresonance the motional reactance cancels c0, so
    ## X = 1/(w*c0) and the residual admittance is the arm's own
    ## conductance transformed by that reactance:
    ##
    ##     Y ~ rm/(rm^2 + X^2) ~ rm*w^2*c0^2   ->   Zp = 1/(rm*wp^2*c0^2)
    ##
    ## which is lm*cm/(rm*c0*(c0 + cm)) once wp is substituted.  Both
    ## spellings are asserted; if they ever disagree one of them is a
    ## typo rather than physics.
    zp = 1.0 / (XTAL['rm'] * (TWOPI * fp) ** 2 * XTAL['c0'] ** 2)
    assert_allclose(zp, XTAL['lm'] * XTAL['cm']
                    / (XTAL['rm'] * XTAL['c0']
                       * (XTAL['c0'] + XTAL['cm'])), rtol=1e-12)
    assert_allclose(z.max(), zp, rtol=2e-3)


def test_crystal_q_from_the_measured_bandwidth():
    """``Q = 2*pi*fs*lm/rm``, recovered from the -3 dB bandwidth of the
    bare motional arm.

    Q is a claim about the SHAPE of the resonance, and the way to test a
    shape is to measure a width: a series RLC has ``|Z| = rm*sqrt(2)``
    at ``fs*(1 +- 1/(2Q))``.  Run with ``c0 = 0`` so the closed form is
    exact rather than approximately exact -- the holder capacitance is
    tested separately, above.
    """
    card = dict(XTAL, c0=0.0)
    fs = 1.0 / (TWOPI * np.sqrt(card['lm'] * card['cm']))
    q = TWOPI * fs * card['lm'] / card['rm']
    assert 1.3e5 < q < 1.4e5
    f = np.array([fs * (1 - 1 / (2 * q)), fs, fs * (1 + 1 / (2 * q))])
    z = np.abs(_xtal_z(card, f))
    assert_allclose(z[1], card['rm'], rtol=1e-9)
    assert_allclose(z[[0, 2]], card['rm'] * np.sqrt(2.0), rtol=1e-4)


def test_crystal_degenerate_cards_stay_finite():
    """No parameter of this model is ever a divisor, so every corner of
    the card is an ordinary element and not a `Collapse` case.

    ``rm = 0`` is a lossless resonator, ``c0 = 0`` a bare motional arm,
    ``cm = 0`` an open one, ``lm = 0`` a plain resistor.  Each is
    checked for finiteness of i, q, G AND C -- the Jacobian goes
    non-finite long before the value does.
    """
    xs = [[a, 0.0, m, i] for a in (0.0, 1.0, -1e6)
          for m in (0.0, 0.5, 1e6) for i in (0.0, 1e-3, -1e6)]
    for card in (dict(XTAL, rm=0.0), dict(XTAL, c0=0.0),
                 dict(XTAL, cm=0.0), dict(XTAL, lm=0.0),
                 dict(lm=0.0, cm=0.0, rm=0.0, c0=0.0)):
        el = _mk(eh.XtalHdl, 'a', 'b', **card)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                for x in xs:
                    x = np.asarray(x, float)
                    for which in ('i', 'q', 'G', 'C'):
                        v = getattr(el, which)(x, defaultepar)
                        assert np.all(np.isfinite(np.asarray(v, float))), \
                            (card, x, which)


def test_bvd_from_fs_q_round_trips_through_a_solve():
    """`bvd_from_fs_q` takes what a datasheet prints; the element built
    from it must resonate where the datasheet said and be as sharp as it
    said.  Measured through an AC solve, not by re-doing the algebra."""
    card = eh.bvd_from_fs_q(fs=10e6, q=8e4, cm=5e-15, c0=2e-12)
    fs = 1.0 / (TWOPI * np.sqrt(card['lm'] * card['cm']))
    assert_allclose(fs, 10e6, rtol=1e-12)
    bare = dict(card, c0=0.0)
    f = np.array([fs * (1 - 1 / (2 * 8e4)), fs, fs * (1 + 1 / (2 * 8e4))])
    z = np.abs(_xtal_z(bare, f))
    assert_allclose(z[1], card['rm'], rtol=1e-9)
    assert_allclose(z[[0, 2]], card['rm'] * np.sqrt(2.0), rtol=1e-4)


def test_crystal_jacobians():
    el = _mk(eh.XtalHdl, 'a', 'b', **XTAL)
    for x in ([0.5, -0.2, 0.1, 1e-3], [0.0, 0.0, 0.0, 0.0],
              [1e6, -1e6, 3.0, -1e3]):
        chk = check_jacobians(el, x)
        assert chk.ok, str(chk)


## ======================================================================
## FerriteBeadHdl -- laplace_nd from network coefficients
## ======================================================================

#: A 600 ohm signal bead, not the class defaults (which are the same
#: numbers -- so `rdc` is moved to make sure it is not being ignored).
BEAD = dict(ls=1.2e-6, rp=600.0, cp=0.35e-12, rdc=0.5)


def _z_bead(f, ls, rp, cp, rdc):
    """``rdc`` in series with a parallel R-L-C, transcribed from the
    network.  `FerriteBeadHdl` states the same thing as one rational
    function; that they agree is the point."""
    s = 1j * TWOPI * np.asarray(f, dtype=complex)
    return rdc + 1.0 / (1.0 / rp + 1.0 / (s * ls) + s * cp)


def _drive_ac(el_factory, freqs, **kw):
    """|Z| of a two-terminal element, driven by a 1 A AC current source
    so the node voltage IS the impedance."""
    c = SubCircuit()
    n1 = c.add_node('n1')
    c['i'] = ISource(gnd, n1, iac=1.0)
    c['dut'] = el_factory(n1, gnd, **kw)
    res = AC(c).solve(freqs=np.atleast_1d(freqs))
    return np.asarray(res.v('n1').y)


def test_bead_impedance_matches_the_parallel_rlc():
    """Six decades against the closed form, magnitude and phase.

    The tolerance is 1e-9 and not 1e-3 deliberately: `laplace_nd` claims
    to realise the transfer function EXACTLY as state equations, so
    anything above solver roundoff is a defect in the realisation and
    not an approximation being tolerated.
    """
    f = np.logspace(3, 11, 120)
    got = _drive_ac(eh.FerriteBeadHdl, f, **BEAD)
    want = _z_bead(f, **BEAD)
    assert_allclose(np.abs(got), np.abs(want), rtol=1e-9)
    assert_allclose(np.angle(got), np.angle(want), rtol=0, atol=1e-9)


def test_bead_landmarks_are_the_datasheet_numbers():
    """The three things a bead datasheet actually tells you, each
    checked against its own closed form and none against the model."""
    f0 = 1.0 / (TWOPI * np.sqrt(BEAD['ls'] * BEAD['cp']))
    assert 2.4e8 < f0 < 2.5e8
    z0, zlo, zhi = _drive_ac(eh.FerriteBeadHdl,
                             [f0, f0 / 1e3, f0 * 1e3], **BEAD)
    ## AT self-resonance the reactances cancel exactly, so the impedance
    ## is rdc + rp and the phase is zero.  Both are identities, not fits.
    assert_allclose(abs(z0), BEAD['rdc'] + BEAD['rp'], rtol=1e-9)
    assert_allclose(np.angle(z0), 0.0, rtol=0, atol=1e-9)
    ## Three decades below, an inductor; three above, a capacitor -- of
    ## the part of the impedance that is NOT `rdc`.  Subtracting it is
    ## not tidying: at 245 kHz the reactance is 1.85 ohm and `rdc` is
    ## 0.5, so the series term is 3.7% of the magnitude and asserting
    ## the bare asymptote at 2e-3 fails, correctly.
    assert_allclose(abs(zlo - BEAD['rdc']),
                    TWOPI * (f0 / 1e3) * BEAD['ls'], rtol=1e-4)
    assert_allclose(np.degrees(np.angle(zlo - BEAD['rdc'])), 90.0,
                    rtol=0, atol=0.2)
    assert_allclose(abs(zhi - BEAD['rdc']),
                    1.0 / (TWOPI * (f0 * 1e3) * BEAD['cp']), rtol=1e-4)
    assert_allclose(np.degrees(np.angle(zhi - BEAD['rdc'])), -90.0,
                    rtol=0, atol=0.2)


def test_bead_is_a_short_at_dc_apart_from_rdc():
    """A real solve, not an evaluation: the dc operating point of a bead
    fed 10 mA must show ``rdc*i`` across it and nothing else, because
    the parallel inductor shorts the rest of the network out."""
    c = SubCircuit()
    n1 = c.add_node('n1')
    c['i'] = ISource(gnd, n1, i=10e-3)
    c['b'] = eh.FerriteBeadHdl(n1, gnd, **BEAD)
    res = DC(c).solve()
    assert_allclose(float(res.v(n1, gnd)), BEAD['rdc'] * 10e-3, rtol=1e-9)


def test_bead_jacobians():
    el = _mk(eh.FerriteBeadHdl, 'a', 'b', **BEAD)
    for x in ([0.5, -0.2, 1e-3, 2.0, 0.01], [0.0] * 5,
              [1e3, -1e3, -5.0, 1e4, -0.3]):
        chk = check_jacobians(el, x)
        assert chk.ok, str(chk)


## ======================================================================
## RSkinHdl -- laplace_zp, and a fractional slope
## ======================================================================

SKIN = dict(rdc=0.5, fs=1e6)


def _z_skin_exact(f, rdc, fs):
    """``rdc*(1 + sqrt(s/ws))`` -- the half-power law the ladder
    approximates, with no ladder in it."""
    return rdc * (1 + np.sqrt(1j * np.asarray(f, dtype=float) / fs))


def test_skin_matches_the_half_power_law_in_band():
    """Magnitude and phase against ``rdc*(1 + sqrt(j*f/fs))`` over the
    seven decades the ladder is advertised for.

    The band is asserted, not merely used: the point of a rational
    approximation to ``sqrt(s)`` is that it has one, and a test that
    only sampled the middle would not notice the day the ladder shrank.
    """
    f = np.logspace(2, 9, 71)          # fs*1e-4 .. fs*1e+3
    got = _drive_ac(eh.RSkinHdl, f, **SKIN)
    want = _z_skin_exact(f, **SKIN)
    ratio = np.abs(got) / np.abs(want)
    dphase = np.degrees(np.angle(got) - np.angle(want))
    assert np.max(np.abs(ratio - 1)) < 0.02, (ratio.min(), ratio.max())
    assert np.max(np.abs(dphase)) < 1.5, (dphase.min(), dphase.max())
    ## and the ripple is genuinely two-sided -- an approximation that
    ## only ever overshoots is a gain error wearing a ripple's clothes
    assert ratio.min() < 0.995 and ratio.max() > 1.005


def test_skin_slope_is_half_a_decade_per_decade():
    """The defining property, measured as a slope rather than assumed
    from the pole placement: ``d log|Z| / d log f -> 1/2`` and the phase
    tends to 45 degrees.

    Measured on ``Z - rdc``, which is the skin term alone.  The
    subtraction is the whole difference between a test that passes and
    one that measures something else: at 10*fs the constant ``rdc``
    still contributes 24% of the magnitude and drags the raw slope to
    0.46, which is neither 0.5 nor a defect -- it is the dc term, and
    the model is not claiming a half-power law for it.
    """
    f = np.logspace(7, 9, 41)          # 10*fs .. 1000*fs
    zac = _drive_ac(eh.RSkinHdl, f, **SKIN) - SKIN['rdc']
    slope = np.diff(np.log10(np.abs(zac))) / np.diff(np.log10(f))
    assert_allclose(np.mean(slope), 0.5, rtol=0, atol=0.01)
    ## The LOCAL slope ripples, because a ladder is a staircase seen
    ## from far enough away: it is 0.5 on average and swings +-0.054
    ## with a period of one cell, 10/11 of a decade.  Measured, and
    ## bounded rather than ignored -- if the ripple ever grows the
    ## placement has drifted, and the mean would not notice.
    assert np.max(np.abs(slope - 0.5)) < 0.07, slope
    assert_allclose(np.mean(np.degrees(np.angle(zac))), 45.0,
                    rtol=0, atol=0.5)


def test_skin_dc_floor_is_the_ladders_own():
    """``H(0) = 1`` for any pole-zero product, but ``sqrt(0) = 0``, so
    the ladder leaves a residue of the skin term at dc.

    It is exactly ``gain``, it is 0.3% of ``rdc`` for the default
    ladder, and it is asserted as an identity rather than bounded as
    "small": a floor that is computable is a property, and a floor that
    is merely small is a surprise waiting for the day someone changes
    the order.
    """
    zeros, poles, gain = eh._fractional_zp(0.5, eh._SKIN_LO, eh._SKIN_HI,
                                           eh._SKIN_PAIRS)
    assert 1e-3 < gain < 1e-2
    c = SubCircuit()
    n1 = c.add_node('n1')
    c['i'] = ISource(gnd, n1, i=1.0)
    c['r'] = eh.RSkinHdl(n1, gnd, **SKIN)
    res = DC(c).solve()
    assert_allclose(float(res.v(n1, gnd)), SKIN['rdc'] * (1 + gain),
                    rtol=1e-9)


def test_laplace_zp_realisation_matches_the_pole_zero_product():
    """The check on the DSL rather than on the physics.

    `laplace_zp` expands its roots into `laplace_nd` coefficients and
    `laplace_nd` builds an order-11 controllable canonical state space
    that the solver integrates.  Here the SAME root list is evaluated
    directly in numpy as ``prod(1 - s/z)/prod(1 - s/p)`` and compared
    against what the AC solve returns -- so a sign convention, a
    conjugate-pair expansion, a coefficient order or a canonical-form
    slip would all show up, and none of them is visible in the physics
    test above (which the ladder's own 2% band would hide).
    """
    zeros, poles, gain = eh._fractional_zp(0.5, eh._SKIN_LO, eh._SKIN_HI,
                                           eh._SKIN_PAIRS)
    ws = TWOPI * SKIN['fs']
    f = np.logspace(2, 10, 33)
    s = 1j * TWOPI * f
    h = np.ones_like(s)
    for k in range(0, len(zeros), 2):
        h = h * (1 - s / (zeros[k] * ws))
    for k in range(0, len(poles), 2):
        h = h / (1 - s / (poles[k] * ws))
    want = SKIN['rdc'] * (1 + gain * h)
    got = _drive_ac(eh.RSkinHdl, f, **SKIN)
    assert_allclose(np.abs(got), np.abs(want), rtol=1e-8)
    assert_allclose(np.angle(got), np.angle(want), rtol=0, atol=1e-8)


def test_skin_ladder_order_is_structural_and_buys_accuracy():
    """``pairs`` is a class property, not an instance parameter -- the
    root list's LENGTH is fixed when ``analog()`` runs.  What it buys is
    band-times-accuracy, and that is measured rather than asserted from
    the placement formula: a five-pair ladder over the same ten decades
    is visibly worse than the default eleven.
    """
    coarse = eh.skin_effect_resistor(pairs=5)
    assert 'pairs' not in [p.name for p in coarse.instparams]
    f = np.logspace(2, 9, 71)
    want = np.abs(_z_skin_exact(f, **SKIN))
    e5 = np.max(np.abs(np.abs(_drive_ac(coarse, f, **SKIN)) / want - 1))
    e11 = np.max(np.abs(np.abs(_drive_ac(eh.RSkinHdl, f, **SKIN))
                        / want - 1))
    assert e5 > 4 * e11, (e5, e11)
    assert e11 < 0.02 < e5


def test_skin_jacobians():
    el = _mk(eh.RSkinHdl, 'a', 'b', **SKIN)
    n = len(x_layout(el))
    for x in (np.zeros(n), np.linspace(-1, 1, n), np.full(n, 1e3)):
        chk = check_jacobians(el, x)
        assert chk.ok, str(chk)


## ======================================================================
## ComparatorHdl -- the first production user of Cross
## ======================================================================

COMP = dict(vref=0.0, vhyst=0.2, gain=20.0, voh=1.0, vol=-1.0)


def _fold_threshold(gain, vhyst):
    """The comparator's real trip voltage, in closed form.

    ``q = tanh(gain*(vin + vhyst*q))`` folds where ``dq/dvin`` diverges,
    i.e. where ``gain*vhyst*(1 - q^2) = 1``.  Solving that for ``q`` and
    putting it back gives the input voltage at the fold.  For
    ``gain*vhyst <= 1`` there is no fold and the threshold is ``vref``.
    """
    lg = gain * vhyst
    if lg <= 1.0:
        return 0.0
    qf = np.sqrt(1.0 - 1.0 / lg)
    return vhyst * qf - np.arctanh(qf) / gain


def _comp_dc(vin, seed, **kw):
    """Solve the comparator's dc operating point from a given latch
    seed, and return the output voltage.  Two seeds, two answers, is
    what bistability MEANS."""
    card = dict(COMP, **kw)
    c = SubCircuit()
    ni, no = c.add_node('i'), c.add_node('o')
    c['vs'] = VS(ni, gnd, v=vin)
    c['X'] = eh.ComparatorHdl(ni, gnd, no, gnd, **card)
    c['Rl'] = R(no, gnd, r=1e3)
    ## The seed is the whole experiment, so it must be a CONSISTENT
    ## starting point and not just a latch value: seeding `lat` while
    ## leaving the input node at zero hands Newton a residual dominated
    ## by the source, and it walks off the branch being asked about.
    ## With `vref = 0` that happened to work and with `vref = -0.7` it
    ## silently did not -- the two seeds converged to one answer and the
    ## element looked monostable.
    x0 = np.zeros(c.n)
    names = [n.name for n in c.nodes]
    x0[names.index('i')] = vin
    x0[names.index('X.lat')] = seed
    x0[names.index('o')] = seed
    return float(DC(c).solve(x0=x0).v(no, gnd))


def test_comparator_is_bistable_exactly_inside_the_fold_voltages():
    """Two dc solves at the same input, seeded high and low.  They must
    differ inside ``vref +- vth`` and agree outside it, and the boundary
    -- located by bisection on the SOLVE, not by evaluating a formula --
    must be the closed-form fold voltage.

    This is the strongest statement of what hysteresis is that does not
    mention time at all.
    """
    vth = _fold_threshold(COMP['gain'], COMP['vhyst'])
    assert_allclose(vth, 0.10735718591, rtol=1e-8)      # 0.537 of vhyst
    assert abs(_comp_dc(0.0, +1) - _comp_dc(0.0, -1)) > 1.9
    assert abs(_comp_dc(+0.3, +1) - _comp_dc(+0.3, -1)) < 1e-6
    assert abs(_comp_dc(-0.3, +1) - _comp_dc(-0.3, -1)) < 1e-6

    ## Bisect for the upper edge of the bistable window.
    lo, hi = 0.0, 0.3
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if abs(_comp_dc(mid, +1) - _comp_dc(mid, -1)) > 1e-6:
            lo = mid
        else:
            hi = mid
    assert_allclose(0.5 * (lo + hi), vth, rtol=1e-6)

    ## and the loop is symmetric about vref -- a structural identity of
    ## an odd nonlinearity, which no fit to one side could reveal
    lo, hi = -0.3, 0.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if abs(_comp_dc(mid, +1) - _comp_dc(mid, -1)) > 1e-6:
            hi = mid
        else:
            lo = mid
    assert_allclose(0.5 * (lo + hi), -vth, rtol=1e-6)


def test_comparator_threshold_moves_with_vref():
    """The trip window is centred on ``vref``, not on zero -- the one
    parameter a bias-dependent test could silently ignore."""
    vth = _fold_threshold(COMP['gain'], COMP['vhyst'])
    for vref in (-0.7, 0.4):
        assert abs(_comp_dc(vref, +1, vref=vref)
                   - _comp_dc(vref, -1, vref=vref)) > 1.9
        assert abs(_comp_dc(vref + 2 * vth, +1, vref=vref)
                   - _comp_dc(vref + 2 * vth, -1, vref=vref)) < 1e-6


def test_comparator_below_unity_loop_gain_has_no_hysteresis():
    """``gain*vhyst <= 1`` has no fold, so the element is a plain
    saturating comparator tripping at ``vref``.

    The seam matters twice over: the model's `Cross` expression divides
    by the loop gain and takes a square root of ``lg - 1``, and this is
    the card where both guards fire.  ``vhyst = 0`` exercises the
    divide-by-zero arm specifically.
    """
    assert _fold_threshold(4.0, 0.2) == 0.0        # lg = 0.8
    for kw in (dict(gain=4.0, vhyst=0.2), dict(gain=20.0, vhyst=0.0)):
        for vin in (-0.5, -0.05, 0.0, 0.05, 0.5):
            hi, lo = _comp_dc(vin, +1, **kw), _comp_dc(vin, -1, **kw)
            assert abs(hi - lo) < 1e-6, (kw, vin, hi, lo)
        ## and it still comparates
        assert _comp_dc(+1.0, -1, **kw) > 0.99
        assert _comp_dc(-1.0, +1, **kw) < -0.99


def test_comparator_output_swing_is_voh_vol():
    """Saturated output equals the declared rails to the tanh's own
    residual, and an asymmetric pair (0 / 3.3 V, a logic output) works
    -- the mapping from ``q`` to the output is affine, which a
    symmetric-rail test would not distinguish from a sign function."""
    card = dict(voh=3.3, vol=0.0)
    assert_allclose(_comp_dc(+1.0, +1, **card), 3.3, rtol=2e-9)
    assert_allclose(_comp_dc(-1.0, -1, **card), 0.0, rtol=0, atol=1e-8)
    ## mid-window, latched high and low: the two rails, not an average
    assert_allclose(_comp_dc(0.0, +1, **card), 3.3, rtol=1e-3)
    assert_allclose(_comp_dc(0.0, -1, **card), 0.0, rtol=0, atol=3e-3)


def _no_cross_twin():
    """`ComparatorHdl` with the two `Cross` statements removed and
    nothing else changed -- the control for the breakpoint measurement.

    Written out rather than parameterised, because a `Cross` cannot be
    switched off by an instance parameter: statements are chosen when
    ``analog()`` runs, i.e. once per class.
    """
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import Behavioural, Branch, Contribution

    class _CompNoCross(Behavioural):
        instparams = [Parameter(name=p.name, desc=p.desc, unit=p.unit,
                                default=p.default)
                      for p in eh.ComparatorHdl.instparams]

        @staticmethod
        def analog(inp, inn, outp, outn):
            bin_ = Branch(inp, inn)
            bout = Branch(outp, outn)
            lat = Node('lat')
            blat = Branch(lat, outn, 'q')
            q = blat.V
            return (Contribution(blat.V,
                                 sympy.tanh(gain * (bin_.V - vref     # noqa
                                                    + vhyst * q))),   # noqa
                    Contribution(bout.V,
                                 vol + (voh - vol) / 2 * (1 + q)))     # noqa
    return _CompNoCross


def _comp_transient(cls, freq=1e3, tend=2e-3, timestep=2e-6, **kw):
    card = dict(COMP, **kw)
    c = SubCircuit()
    ni, no = c.add_node('i'), c.add_node('o')
    c['vs'] = VSin(ni, gnd, va=1.0, freq=freq)
    c['X'] = cls(ni, gnd, no, gnd, **card)
    c['Rl'] = R(no, gnd, r=1e3)
    tran = Transient(c, toolkit=numeric, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=tend, timestep=timestep)
    return (tran.statistics,
            np.asarray(res.v('o').x[0], float),
            np.asarray(res.v('i').y, float),
            np.asarray(res.v('o').y, float))


def test_hysteresis_latches_in_a_transient_both_directions():
    """The loop, in time.  A 1 kHz sine through the comparator must
    switch high only above ``+vth`` and low only below ``-vth``, four
    times in two periods, and must NOT switch at ``vref``.

    Both directions and both thresholds, because a comparator that
    tripped at ``+vth`` in both directions would pass every test that
    looked only at the rising edge.
    """
    vth = _fold_threshold(COMP['gain'], COMP['vhyst'])
    _s, t, vi, vo = _comp_transient(eh.ComparatorHdl)
    edges = np.nonzero(np.diff(np.sign(vo)))[0]
    edges = [k for k in edges if t[k] > 0]      # skip the uic start
    assert len(edges) == 4, len(edges)
    rising = [k for k in edges if vo[k + 1] > vo[k]]
    falling = [k for k in edges if vo[k + 1] < vo[k]]
    assert len(rising) == 2 and len(falling) == 2
    ## Each switch must BRACKET its threshold.  The tolerance is 1e-5 V
    ## and it is there because `Cross` succeeds: three of the four edges
    ## land a timepoint ON the threshold to 1e-16 s, so the "sample
    ## before" and the "sample after" are the same voltage to nine
    ## digits and a strict inequality fails by 5e-9 V.  A bracket that
    ## has collapsed onto the answer is the best possible outcome here,
    ## not a failure.
    tol = 1e-5
    for k in rising:
        assert min(vi[k], vi[k + 1]) - tol <= vth <= \
            max(vi[k], vi[k + 1]) + tol, (vi[k], vth, vi[k + 1])
    for k in falling:
        assert min(vi[k], vi[k + 1]) - tol <= -vth <= \
            max(vi[k], vi[k + 1]) + tol, (vi[k], -vth, vi[k + 1])
    ## and it HOLDS its state across vref: wherever the input is well
    ## inside the loop, the output is still hard against a rail.  Stated
    ## over the whole waveform rather than as "no edge near vref",
    ## because which sample happens to bracket an edge is the step
    ## controller's business and moves with the compilation path.
    ## after the first edge, so the uic start-up is not counted: `uic`
    ## puts the latch at q = 0, its metastable point, and the first
    ## fifteen samples are it resolving off that.
    inside = (np.abs(vi) < 0.5 * vth) & (t > t[edges[0] + 1])
    assert inside.sum() >= 4, inside.sum()
    assert np.all(np.abs(vo[inside]) > 0.9), np.abs(vo[inside]).min()


def test_cross_lands_the_comparator_edges():
    """`Cross`, measured against the same element without it.

    The metric is the distance from a threshold crossing -- whose time
    is known analytically, ``asin(vth/va)/(2*pi*f)`` for a sine -- to
    the nearest timepoint the controller accepted.  That is what
    `Cross` promises and the only thing it promises.

    Measured at introduction, on a 1 ms period: without `Cross` the
    nearest timepoint to a crossing is 7.4e-7 to 1.2e-5 s away.  With
    it, THREE of the four edges are resolved to 3e-16 s and the fourth
    to 6.6e-7 -- and that scatter is a property of the predictor worth
    knowing rather than noise.  `Cross` extrapolates linearly from the
    last two accepted points and publishes one prediction; whether it
    then gets a SECOND, much better one depends on whether the
    controller polls again after landing near the corner.  Here the
    first falling edge got one refinement and the other three got two,
    so the same run places its edges ten decades apart in accuracy.
    There is no iteration to convergence, and a first-order predictor
    over a 4e-5 s step on a sine is exactly 6.6e-7 s wrong.

    It costs about 3x the accepted steps, because a breakpoint restarts
    the step controller small and it has to grow back; that cost is
    asserted as a CEILING so a future regression in it is visible, not
    hidden behind "Cross is free".
    """
    vth = _fold_threshold(COMP['gain'], COMP['vhyst'])
    va, f = 1.0, 1e3
    tc = []
    for n in range(2):
        tc += [n / f + np.arcsin(vth / va) / (TWOPI * f),
               n / f + (np.pi + np.arcsin(vth / va)) / (TWOPI * f)]

    s_off, t_off, _vi, _vo = _comp_transient(_no_cross_twin())
    s_on, t_on, _vi, _vo = _comp_transient(eh.ComparatorHdl)
    near_off = [float(np.min(np.abs(t_off - x))) for x in tc]
    near_on = [float(np.min(np.abs(t_on - x))) for x in tc]

    assert s_off.breakpoints_hit == 0
    assert s_on.breakpoints_hit >= 4              # at least one per edge
    ## the WORST edge, which is the one a step controller has to live
    ## with: 17x measured, asserted at 10x
    assert max(near_on) < max(near_off) / 10.0, (near_on, near_off)
    assert max(near_on) < 1e-6
    ## and the best edges really do land on the corner itself.  How MANY
    ## of the four do is not stable -- it depends on whether the
    ## controller polled again after landing near the corner, and that
    ## moved from three to two when the model was made chained -- so the
    ## claim is about the best edge, which is the one that shows what a
    ## converged prediction is worth.
    assert min(near_on) < 1e-15, near_on
    ## the price, bounded.  3x measured at introduction, first asserted
    ## at 4x; under the C backend (2026-08-26) the same run accepts
    ## 4.1x -- numpy's own tanh differs from libm's by an ulp, the
    ## comparator sits on tanh, and a handful of step-acceptance
    ## decisions flip.  The ceiling still binds at 4.5x.
    assert s_on.accepted_steps < 4.5 * s_off.accepted_steps
    assert s_on.rejected_steps <= s_off.rejected_steps + 2


def test_cross_predictions_use_the_fold_not_vhyst():
    """The reason the previous test works, isolated.

    ``vhyst`` is 0.2 and the real threshold is 0.107 -- they differ by a
    factor of 1.9, so a `Cross` declared at ``vhyst`` would land its
    timepoints 30 us away from the edge on this drive.  Assert that no
    accepted timepoint clusters at the ``vhyst`` crossing while several
    cluster at the fold: that distinguishes "Cross fired" from "Cross
    fired in the right place", and only the second is worth having.
    """
    va, f = 1.0, 1e3
    vth = _fold_threshold(COMP['gain'], COMP['vhyst'])
    t_fold = np.arcsin(vth / va) / (TWOPI * f)
    t_vh = np.arcsin(COMP['vhyst'] / va) / (TWOPI * f)
    assert t_vh - t_fold > 1e-5
    _s, t, _vi, _vo = _comp_transient(eh.ComparatorHdl)
    assert float(np.min(np.abs(t - t_fold))) < 1e-8
    assert float(np.min(np.abs(t - t_vh))) > 1e-7


def test_comparator_jacobians_and_finiteness():
    """Finite differences at a latched point, at the fold itself and at
    absurd bias; plus overflow/invalid/divide raised as errors over a
    sweep including the ``vhyst = 0`` and sub-unity-loop-gain cards."""
    for kw in ({}, dict(vhyst=0.0), dict(gain=4.0, vhyst=0.2)):
        el = _mk(eh.ComparatorHdl, 'a', 'b', 'c', 'd', **dict(COMP, **kw))
        for x in ([0.05, 0.0, 0.7, 0.0, 0.6, 0.0, 1e-4],
                  [_fold_threshold(COMP['gain'], COMP['vhyst']),
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                  [-3.0, 0.0, -1.0, 0.0, -0.999, 0.0, 1e-4],
                  [1e6, -1e6, 1.0, 0.0, 3.0, 0.0, -1e3]):
            chk = check_jacobians(el, x)
            assert chk.ok, (kw, x, str(chk))
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                for vin in (0.0, 1e-12, 0.5, 1e6, -1e6):
                    for q in (-3.0, -1.0, 0.0, 1.0, 3.0):
                        x = np.array([vin, 0.0, 0.5, 0.0, q, 0.0, 1e-4])
                        for which in ('i', 'q', 'G', 'C'):
                            v = np.asarray(getattr(el, which)(x, defaultepar),
                                           float)
                            assert np.all(np.isfinite(v)), (kw, x, which)


## ======================================================================
## MemristorHdl -- idt as a real DAE state
## ======================================================================

#: Biolek's own published card (Radioengineering 18(2), fig. 3).
MEM = dict(ron=100.0, roff=16e3, muv=1e-14, d=10e-9, ic=0.1, p=1.0)


def _mem_rk4(t_end, n, va, freq, ron, roff, muv, d, ic, p):
    """The Strukov/Biolek system, integrated in numpy with fixed-step
    RK4 and no reference to `elements_hdl`.

    ``dx/dt = muv*ron/d^2 * v(t)/R(x) * (1 - (x - stp(-i))^(2p))`` with
    ``R(x) = ron*x + roff*(1-x)``.  Returns ``(t, x)``.
    """
    k = muv * ron / d ** 2

    def f(tt, x):
        xc = min(max(x, 0.0), 1.0)
        v = va * np.sin(TWOPI * freq * tt)
        i = v / (ron * xc + roff * (1 - xc))
        stp = 1.0 if i < 0 else 0.0
        return k * i * (1 - abs(x - stp) ** (2 * p))

    ts = np.linspace(0.0, t_end, n + 1)
    h = ts[1] - ts[0]
    xs = np.empty(n + 1)
    xs[0] = ic
    for j in range(n):
        tt, x = ts[j], xs[j]
        k1 = f(tt, x)
        k2 = f(tt + h / 2, x + h * k1 / 2)
        k3 = f(tt + h / 2, x + h * k2 / 2)
        k4 = f(tt + h, x + h * k3)
        xs[j + 1] = x + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return ts, xs


def _mem_transient(va=1.0, freq=1.0, cycles=2, pts=2000, **kw):
    """One memristor driven by a sine, with the step size CAPPED.

    ``timestep`` alone is the opening step; the controller grows it, and
    a 2 s run at ``timestep = 2.5e-4`` takes 69 accepted steps whatever
    ``timestep`` says.  That is correct behaviour and it makes
    ``timestep`` useless as an accuracy knob -- the first version of
    `test_memristor_state_matches_an_independent_rk4` swept it over
    three decades and got the identical 1.3% error every time, which
    reads exactly like a model discrepancy.  ``timestep_max`` is the
    knob that means what ``timestep`` looks like it means.
    """
    card = dict(MEM, **kw)
    c = SubCircuit()
    n1 = c.add_node('n1')
    c['vs'] = VSin(n1, gnd, va=va, freq=freq)
    c['m'] = eh.MemristorHdl(n1, gnd, **card)
    tran = Transient(c, toolkit=numeric, uic=True,
                     timestep_max=1.0 / (freq * pts))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=cycles / freq, timestep=1.0 / (freq * pts))
    return (np.asarray(res.v('n1').x[0], float),
            np.asarray(res.v('n1').y, float),
            -np.asarray(res.i('vs.plus').y, float),
            np.asarray(res.v('m._state0').y, float))


def test_memristor_state_matches_an_independent_rk4():
    """The state is a real DAE unknown the solver integrates; check it
    against a fixed-step RK4 of the same ODE written in numpy.

    This is the test that would catch a wrong routing of the ``idt``
    state row (``q[s] = s``, ``i[s] = -arg``), a sign slip in the scratch
    node's potential contribution, or a window transcribed backwards --
    none of which the i-v loop shape would reveal, because a wrong-sign
    window still produces a loop.
    """
    tr, xr = _mem_rk4(2.0, 400000, 1.0, 1.0, **MEM)
    swing = xr.max() - xr.min()
    assert swing > 0.3                       # the reference actually moves

    errs = []
    for pts in (500, 1000, 2000):
        t, _v, _i, x = _mem_transient(va=1.0, freq=1.0, cycles=2, pts=pts)
        errs.append(np.max(np.abs(x - np.interp(t, tr, xr))))
    ## The agreement is asserted twice, and the second assertion is the
    ## one that means something: a single tolerance can be met by
    ## coincidence, but an error that falls by 4x for every halving of
    ## the step is the solver CONVERGING to this reference and to
    ## nothing else.  Gear-2 is second order and that is what is seen.
    assert errs[-1] < 1e-5 * swing, errs
    assert errs[0] / errs[1] > 3.0 and errs[1] / errs[2] > 3.0, errs


def test_memristor_is_pinched_at_the_origin():
    """``|i| <= |v|/ron`` at every timepoint, with equality only at
    ``x = 1``.

    That inequality IS the pinch: a memristor is a resistor whose value
    is remembered, so it can carry no current at zero volts however the
    state got there.  Asserted as a bound over the whole waveform rather
    than at the crossing samples, because the crossing sample is
    whichever point the controller happened to accept.
    """
    t, v, i, x = _mem_transient(va=1.2, freq=1.0, cycles=2, pts=1000)
    assert np.max(np.abs(i)) > 1e-5          # it is conducting at all
    assert np.all(np.abs(i) <= np.abs(v) / MEM['ron'] + 1e-15)
    ## and the loop really is open: the current at +0.6 V on the way up
    ## differs from the current at +0.6 V on the way down by far more
    ## than the sampling could explain
    sel = np.abs(v - 0.6) < 0.02
    assert sel.sum() >= 4, sel.sum()
    assert (i[sel].max() - i[sel].min()) / i[sel].max() > 0.05


def test_memristor_loop_collapses_at_high_frequency():
    """The signature.  Raise the drive frequency and the state has less
    time to move per cycle, so the lobes close and the device degenerates
    into the linear resistor ``R(ic)`` -- a closed form with no memristor
    physics in it at all.
    """
    r_ic = MEM['ron'] * MEM['ic'] + MEM['roff'] * (1 - MEM['ic'])
    spreads = []
    for freq in (1.0, 10.0, 1000.0):
        t, v, i, x = _mem_transient(va=1.0, freq=freq, cycles=2, pts=800)
        spreads.append(x.max() - x.min())
        if freq == 1000.0:
            ## the loop has closed onto Ohm's law at the initial state
            assert_allclose(i, v / r_ic, rtol=5e-3, atol=1e-12)
    ## and the closing is monotone in frequency
    assert spreads[0] > 3 * spreads[1] > 9 * spreads[2]
    assert spreads[2] < 1e-2 and spreads[0] > 0.3


def test_biolek_window_keeps_the_state_in_range_and_releases_it():
    """Drive it hard enough to slam into both boundaries.

    Two claims, and the second is the one Joglekar's window fails:

    * the state never leaves [0, 1] by more than the integrator's own
      local error, however hard it is driven;
    * having reached a boundary it comes BACK when the current reverses.
      A window that vanishes at both ends regardless of current sign
      would leave it stuck, and every other test in this file would
      still pass.
    """
    t, v, i, x = _mem_transient(va=3.0, freq=0.5, cycles=2, pts=1200,
                                ic=0.5)
    assert x.max() > 0.97, x.max()           # it does reach the top
    assert x.min() < 0.2, x.min()            # and comes back down
    assert -1e-6 <= x.min() and x.max() <= 1 + 1e-6, (x.min(), x.max())
    ## released: after the maximum, the state must fall again
    kmax = int(np.argmax(x))
    assert kmax < len(x) - 10
    assert x[kmax] - x[kmax:].min() > 0.5


def test_memristor_dc_pins_the_state_and_obeys_ohm_there():
    """A real dc solve.  ``idt(., ic)`` replaces the state row with
    ``s = ic`` at dc, so the operating point is Ohm's law at ``R(ic)``
    -- checked at three different ``ic`` because a model that ignored
    ``ic`` and started from zero would pass at one of them."""
    for ic in (0.0, 0.25, 1.0):
        r_ic = MEM['ron'] * ic + MEM['roff'] * (1 - ic)
        c = SubCircuit()
        n1 = c.add_node('n1')
        c['vs'] = VS(n1, gnd, v=0.4)
        c['m'] = eh.MemristorHdl(n1, gnd, **dict(MEM, ic=ic))
        res = DC(c).solve()
        assert_allclose(float(res.v('m._state0')), ic, rtol=0, atol=1e-12)
        assert_allclose(-float(res.i('vs.plus')), 0.4 / r_ic, rtol=1e-9)


def test_uic_seeds_the_state_from_ic():
    """``uic=True`` must start the transient from ``ic`` and not from
    zero.

    Pinned because it is spelling-sensitive and fails silently:
    `Transient`'s uic walk only calls a generated ``state_ic()`` for an
    element that has an instance parameter literally named ``ic``, so a
    memristor whose initial state were called ``x0`` would start at
    ``x = 0`` -- the fully undoped film, 16 kohm instead of 14.4 -- with
    no warning anywhere and a dc solve that was perfectly correct.
    """
    assert 'ic' in [p.name for p in eh.MemristorHdl.instparams]
    for ic in (0.05, 0.4, 0.9):
        t, v, i, x = _mem_transient(va=0.1, freq=1e4, cycles=1, pts=200,
                                    ic=ic)
        assert_allclose(x[0], ic, rtol=0, atol=1e-12)


def test_memristor_jacobians_and_finiteness():
    """Finite differences at ordinary bias, at both window seams, and at
    the state where ``R(x)`` would pass through zero.

    ``x = roff/(roff - ron) = 1.00629`` is six thousandths outside the
    physical range, which is well inside where a Newton iterate goes.
    The clamp in the model is what keeps it finite, so this is the point
    the test exists for.
    """
    xsing = MEM['roff'] / (MEM['roff'] - MEM['ron'])
    assert_allclose(xsing, 1.006289, rtol=1e-5)
    ## EVERY CARD, EVERY STATE, AND THE EXTREME BIASES, and until
    ## `check_jacobians` grew its UNRESOLVED verdict none of those three
    ## could be swept.  Each was excluded for a different reason and all
    ## three were properties of central differencing, not of the model:
    ##
    ## * ``x = 0`` and ``x = 1`` are the corners of the
    ##   ``minc(maxc(x, 0), 1)`` clamp, so the difference straddles a
    ##   jump in the SLOPE and returns the average of the two arms while
    ##   the Jacobian returns one of them.  No ``h`` helps; a jump has no
    ##   scale.  The instrument now detects the kink from the VALUE
    ##   alone -- the one-sided disagreement stays put as ``h`` shrinks
    ##   where a smooth function's halves -- and says so;
    ## * ``ron = 1, roff = 1e9`` has ``dR/dx = 1e9``, so a 1e-7 step is a
    ##   1% change in ``R`` and the difference carries 1e-4 of truncation
    ##   error.  Four points of this sweep FAILED before the fix (worst:
    ##   ``G = 9.998e1`` against a difference of ``9.999e1`` at
    ##   ``x = 1 - 1e-5``).  There is still no single ``h`` that serves
    ##   both this card and the default one -- measured: 1e-9 fixes this
    ##   card and breaks the default at ``x = 20`` -- which is why the
    ##   truncation term is now MEASURED per entry by Richardson rather
    ##   than tuned away;
    ## * at ``|v| = 1e4`` with ``p = 3`` the drive is 1e6 and the state
    ##   moves it by 1e-11 over the difference step, below one ulp of the
    ##   value, so the difference measured roundoff (5.8e-4 against an
    ##   analytic 5.9e-5).
    npts, resolved_pts, corners = 0, 0, 0
    for kw in ({}, dict(p=3.0), dict(ron=1.0, roff=1e9)):
        el = _mk(eh.MemristorHdl, 'a', 'b', **dict(MEM, **kw))
        for xst in (0.0, 1e-5, 0.5, 1 - 1e-5, 1.0, xsing, 1.5, -0.4, 20.0):
            for v in (0.7, -0.7, 10.0, -10.0, 1e4, -1e4):
                chk = check_jacobians(el, [v, 0.0, 1e-5, xst, 0.0])
                assert chk.ok, (kw, xst, v, str(chk))
                npts += 1
                resolved_pts += bool(chk.resolved)
                ## The clamp corners must be REPORTED as kinks, not
                ## passed quietly: on the two soft cards that is the only
                ## thing wrong with them, so the reason is pinned exactly.
                if kw.get('ron') is None and xst in (0.0, 1.0):
                    corners += 1
                    assert {u.reason for u in chk.unresolved} == {'kink'}, \
                        (kw, xst, v, str(chk))
    ## AND MOST OF THE SWEEP MUST STILL RESOLVE.  Every `chk.ok` above
    ## would also be satisfied by an instrument that had quietly stopped
    ## checking and called everything unresolved, so the count is pinned.
    ## Measured 2026-08-25: 122 of 162 points fully resolved; of the 40
    ## that are not, 24 are the clamp corners of the two soft cards, 12
    ## are the stiff card at `x = 1 - 1e-5` and `x = 1` (truncation), 3
    ## are its corner at `x = 0`, and one is `p = 3` at `x = 1.00629,
    ## v = -1e4` (roundoff -- mechanism 2, on a card that has it).
    assert npts == 162
    assert corners == 24
    assert resolved_pts >= 122, resolved_pts
    ## Finiteness is swept over a WIDER set than the finite-difference
    ## check: |v| up to 1e6 and states out to +-1e3.
    for kw in ({}, dict(p=3.0), dict(ron=1.0, roff=1e9)):
        el = _mk(eh.MemristorHdl, 'a', 'b', **dict(MEM, **kw))
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                for xst in (-1e3, -1.0, 0.0, xsing, 1.0, 1e3):
                    for v in (0.0, 1e-12, 1.0, -1.0, 1e6, -1e6):
                        x = np.array([v, 0.0, 1e-5, xst, 0.0])
                        for which in ('i', 'q', 'G', 'C'):
                            val = np.asarray(getattr(el, which)(x,
                                                                defaultepar),
                                             float)
                            assert np.all(np.isfinite(val)), (kw, x, which)


def test_the_model_has_three_kinks_and_they_are_where_the_physics_is():
    """A recorded property, not a defect.

    The model is C0 everywhere and C1 everywhere except at exactly
    three places, each of which is a physical boundary:

    * ``i = 0`` -- Biolek's ``stp(-i)`` is a hard step, so the window
      jumps from ``1 - x^(2p)`` to ``1 - (x-1)^(2p)`` as the current
      reverses.  The state equation itself is continuous, because the
      drive is ``i`` times the window and ``i`` is zero at the seam, so
      the trajectory is C1 and the simulator may integrate it;
    * ``x = 0`` and ``x = 1`` -- the corners of the resistance clamp,
      where the film is fully undoped or fully doped and the resistance
      has genuinely saturated.

    `check_jacobians` central-differences with a step of 1e-7, so at any
    of the three it straddles a jump and reports a discrepancy that is
    arithmetically correct and says nothing about the model.  Shrinking
    ``h`` does not help; a jump has no scale.

    (Until 2026-08-25 that was also the reason the sweep above skipped
    these points.  It no longer is: `check_jacobians` detects the kink
    from the VALUE -- the one-sided disagreement stays put as ``h``
    shrinks where a smooth function's halves -- and reports UNRESOLVED,
    so the corners are swept and asserted there.  This test keeps the
    stronger, model-side claim, which no verdict can make.)

    What IS checked here is the pair of claims that make the kinks
    acceptable: the value is continuous across each (the difference
    across the seam falls linearly with the step), and the one-sided
    derivatives are the two finite numbers the equations say they are.
    """
    el = _mk(eh.MemristorHdl, 'a', 'b', **MEM)
    k = MEM['muv'] * MEM['ron'] / MEM['d'] ** 2
    row = 4                       # the scratch node's V-branch row, dx/dt

    def onesided(x0, col, h):
        f0 = np.asarray(el.i(x0, defaultepar), float)
        e = np.zeros(5)
        e[col] = h
        fp = np.asarray(el.i(x0 + e, defaultepar), float)
        fm = np.asarray(el.i(x0 - e, defaultepar), float)
        return (fp - f0) / h, (f0 - fm) / h, fp - fm

    ## -- the current-reversal seam, at x = 0 where the two window arms
    ##    are 1 and 0 and so differ maximally
    x0 = np.array([0.0, 0.0, 1e-5, 0.0, 0.0])
    up, dn, _d = onesided(x0, 0, 1e-9)
    ## `i[row]` is `-(V(xs,minus) - drive)`, so d/dv is +d(drive)/dv
    assert_allclose(up[row], k / MEM['roff'], rtol=1e-6)
    assert_allclose(dn[row], 0.0, rtol=0, atol=1e-12)

    ## -- the two clamp corners, in the state
    for xc, dlo, dhi in ((0.0, 0.0, 1.0), (1.0, 1.0, 0.0)):
        xx = np.array([0.7, 0.0, 1e-5, xc, 0.0])
        up, dn, _d = onesided(xx, 3, 1e-9)
        ## below x = 0 (and above x = 1) the resistance is frozen, so the
        ## terminal current does not respond to the state at all
        r = MEM['ron'] * xc + MEM['roff'] * (1 - xc)
        slope = -0.7 * (MEM['ron'] - MEM['roff']) / r ** 2
        assert_allclose(up[0], slope * dhi, rtol=1e-5, atol=1e-12)
        assert_allclose(dn[0], slope * dlo, rtol=1e-5, atol=1e-12)
        assert abs(up[0] - dn[0]) > 0.5 * abs(slope)

    ## -- and every one of the three is CONTINUOUS in value: the jump
    ##    across the seam falls linearly with the step, which a genuine
    ##    discontinuity's would not
    for x0, col in ((np.array([0.0, 0.0, 1e-5, 0.0, 0.0]), 0),
                    (np.array([0.7, 0.0, 1e-5, 0.0, 0.0]), 3),
                    (np.array([0.7, 0.0, 1e-5, 1.0, 0.0]), 3)):
        jumps = [np.max(np.abs(onesided(x0, col, h)[2])) for h in
                 (1e-6, 1e-8, 1e-10)]
        assert jumps[0] > 0
        assert_allclose(jumps[1] / jumps[0], 1e-2, rtol=1e-3)
        assert_allclose(jumps[2] / jumps[1], 1e-2, rtol=1e-3)


def test_memristor_scratch_node_injects_nothing():
    """The idiom's correctness condition, asserted rather than argued.

    ``xs`` is written by a potential contribution and nothing else
    touches it, so KCL there forces its branch current to zero and the
    node contributes exactly nothing to the terminal currents.  If that
    were not true the memristor would be a current source as well as a
    resistor, and every i-v test above would still pass because the
    error is at the OTHER terminal.
    """
    t, v, i, x = _mem_transient(va=1.0, freq=1.0, cycles=1, pts=2000)
    c = SubCircuit()
    n1 = c.add_node('n1')
    c['vs'] = VS(n1, gnd, v=0.5)
    c['m'] = eh.MemristorHdl(n1, gnd, **MEM)
    res = DC(c).solve()
    ## the branch current of the scratch node's own branch
    el = c['m']
    names = [n.name for n in c.nodes]
    ## dc: current into the source equals V/R(ic) exactly, so nothing
    ## else is being injected at either terminal
    r_ic = MEM['ron'] * MEM['ic'] + MEM['roff'] * (1 - MEM['ic'])
    assert_allclose(-float(res.i('vs.plus')), 0.5 / r_ic, rtol=1e-12)
    ## and in transient the same identity holds sample by sample
    r_x = MEM['ron'] * np.clip(x, 0, 1) + MEM['roff'] * (1 - np.clip(x, 0, 1))
    assert_allclose(i, v / r_x, rtol=1e-9, atol=1e-15)


## ======================================================================
## What the compiler made of each of them
## ======================================================================


def test_explain_reports_the_capability_each_model_was_chosen_for():
    """`hdl.explain` is the record of what the compiler DID, and this
    pins the one line of it that each model exists to produce.

    It is a structural test and it earns its place twice: it is the only
    thing here that would notice `laplace_zp` silently building a filter
    of the wrong ORDER (the conjugate-pair trap in its own source: one
    ``(re, im)`` entry is a PAIR of roots, so listing both members
    doubles the order and every magnitude test still passes over most of
    the band), and it is how the ``ic`` DC pin is shown to be present
    rather than merely believed.
    """
    from pycircuit.circuit.hdl import explain

    def report(cls, nterm, **kw):
        el = _mk(cls, *[chr(97 + k) for k in range(nterm)], **kw)
        return explain(el, source=False, symbolic=False), x_layout(el)

    txt, lay = report(eh.XtalHdl, 2, **XTAL)
    assert 'state' not in txt.split('parameters')[0].replace('x vector', '')
    assert [k for _i, _n, k in lay].count('internal node') == 1

    txt, lay = report(eh.FerriteBeadHdl, 2, **BEAD)
    assert '2 states' in txt          # one per order of the denominator
    txt, lay = report(eh.RSkinHdl, 2, **SKIN)
    assert '%d states' % eh._SKIN_PAIRS in txt
    assert [k for _i, _n, k in lay].count('state') == eh._SKIN_PAIRS

    txt, lay = report(eh.ComparatorHdl, 4, **COMP)
    assert '2 @cross events' in txt
    assert 'state' not in txt.split('parameters')[0].replace('x vector', '')
    ## CHAINED, on purpose.  `Cross` was installed after an early return
    ## taken by every `var()`-using model, so it did nothing at all on
    ## the let-chain path until 2026-08-25 -- and `explain()` reported
    ## it as present the whole time, because it reads the compiler's
    ## record rather than the class.  So the text above is NOT evidence
    ## on its own; the class attribute is.
    assert eh.ComparatorHdl._hdl_info['chained'] is True
    el = _mk(eh.ComparatorHdl, 'a', 'b', 'c', 'd', **COMP)
    assert hasattr(el, 'accept_step') and hasattr(el, 'next_event')
    el.accept_step(0.0, np.array([-1.0, 0, -1, 0, -1, 0, 0]))
    el.accept_step(1e-4, np.array([0.0, 0, -1, 0, -1, 0, 0]))
    assert np.isfinite(el.next_event(1e-4))

    txt, lay = report(eh.MemristorHdl, 2, **MEM)
    assert '1 state' in txt
    assert any('pinned at DC by its ic' in k for _i, _n, k in lay)


def test_var_no_longer_hides_the_state_operators_from_the_compiler():
    """This test EXPIRED, which is what it was written to do.

    It used to record the reason `MemristorHdl` is written flat: the
    compiler collected `idt`, `idtmod` and `laplace_*` applications by
    walking the CONTRIBUTION STATEMENTS, so an application appearing
    only inside a `var()` definition never got a state allocated,
    survived to the printer and died there with ``Unsupported by _P:
    idt`` -- nine frames deep, naming neither `var()` nor the model.

    `$limit` had exactly the same bug and was fixed first (`hdl.py`,
    "Intermediates FIRST, and this is not tidiness"); the state
    operators now take the same treatment, so both halves assert a
    SUCCESS.  Written as a test rather than as a comment precisely
    because a comment recording a gap has no owner and does not fail
    when the gap closes -- this one failed the day it closed.
    `test_chained_first_class.py` carries the detailed pins.
    """
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       var, idt, idtmod, laplace_nd,
                                       limit_pnj, vt)

    def build(rhs_of):
        class _T(Behavioural):
            instparams = [Parameter(name='g', desc='g', unit='S',
                                    default=1e-3)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                xn = Node('xs')
                bx = Branch(xn, minus, 'd')
                return (Contribution(b.I, g * b.V * rhs_of(b, bx)),  # noqa
                        Contribution(bx.V, b.V))
        return _T

    for op, nstates in ((lambda b, bx: var(idt(bx.V, 0.5), 'x'), 1),
                        (lambda b, bx: var(idtmod(bx.V, 0.0, 1.0, 0.0),
                                           'x'), 1),
                        (lambda b, bx: var(laplace_nd(b.V, [1], [1, 1e-6]),
                                           'x'), 1)):
        cls = build(op)
        ## Compiling is not the claim -- the STATE is.  A compiler that
        ## dropped the application silently would also "compile".
        assert cls._hdl_info['chained'] is True
        assert len(cls._hdl_info['state_meta']['statenames']) == nstates

    ## and the one that was fixed
    class _Lim(Behavioural):
        instparams = [Parameter(name='IS', desc='sat', unit='A',
                                default=1e-14)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            v = var(limit_pnj(b.V, IS, vt()), 'vlim')                 # noqa
            return Contribution(b.I, IS * (sympy.exp(v / vt()) - 1))  # noqa
    assert _Lim._hdl_info['limit_spec']
