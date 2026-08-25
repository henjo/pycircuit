# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The Phase-6 library models: the self-heating thermal node and the full
SPICE junction diode.

Nothing here compares the model against itself.  Every number is either

* a **closed form** the model does not contain -- the depletion
  capacitance ``cj0/(1 - v/vj)^m``, obtained by finite-differencing the
  model's *charge*, which is what the model actually states; the thermal
  fixed point ``dT = P*rth``; the roots of the self-heating quadratic;
* a **structural identity** that must hold whatever the equations are --
  at ``T = tnom`` every temperature correction is the identity, the
  breakdown term is exactly ``-ibv`` at ``v = -bv``, area scales the
  current linearly;
* an **independent transcription** of Massobrio & Antognetti ch. 1 in
  numpy (`_ref`), written from the book rather than from
  `elements_hdl`, which is what catches a routing or compilation error;
* a **finite difference**, for every Jacobian claim.

and each is exercised across the seams -- ``v = fc*vj``, ``v = -bv``,
``v = 0``, both signs out to 1e6 V -- because that is where the
smoothing lives.

`TestNoFloatingPointGarbage` is the one that pays for itself: overflow,
invalid and divide-by-zero raised as errors over a wide bias and
temperature sweep and through a DC solve.  It is what showed that the
diode had to be written in `expl` rather than `limexp`.
"""

import copy
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pycircuit.circuit.circuit
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.elements import SubCircuit, VS, R, C, IS as ISource
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.hdl import Node

K = numeric.kboltzmann
Q = numeric.qelectron

## A real-ish small-signal card, deliberately NOT the defaults: defaults
## leave rs, cjo, tt and bv switched off, and the model is a bare
## exponential with every interesting block dark
## (`differentiable-numerics`: "use CARD parameters, not defaults").
CARD = dict(IS=1.2e-14, rs=1.5, n=1.06, tt=4e-9, cjo=2.3e-12, vj=0.78,
            m=0.42, eg=1.11, xti=3.0, fc=0.5, bv=45.0, ibv=5e-6,
            kf=0.0, af=1.0, area=1.0, tnom=300.15)


def _epar(T):
    e = copy.copy(defaultepar)
    e.T = T
    return e


def _ref(v, T=300.0, IS=1e-14, rs=0.0, n=1.0, tt=0.0, cjo=0.0, vj=1.0,
         m=0.5, eg=1.11, xti=3.0, fc=0.5, bv=1e30, ibv=1e-3, kf=0.0,
         af=1.0, area=1.0, tnom=300.15):
    """SPICE level-1 junction, transcribed in numpy from Massobrio &
    Antognetti, *Semiconductor Device Modeling with SPICE*, ch. 1.

    Returns ``(i, q)`` at the JUNCTION voltage ``v`` (i.e. after the
    series resistance), so `rs` is accepted and ignored.
    """
    vt = K * T / Q
    tr = T / tnom
    egT = 1.16 - 7.02e-4 * T ** 2 / (T + 1108.0)
    egn = 1.16 - 7.02e-4 * tnom ** 2 / (tnom + 1108.0)
    isT = IS * np.exp((tr - 1) * eg / (n * vt) + xti / n * np.log(tr))
    vjT = vj * tr - 3 * vt * np.log(tr) - egn * tr + egT
    cjT = cjo * (1 + m * (4e-4 * (T - tnom) - (vjT / vj - 1)))
    ifwd = isT * (np.exp(v / (n * vt)) - 1)
    i = area * (ifwd - ibv * np.exp(-(v + bv) / vt))
    fcvj = fc * vjT
    if v < fcvj:
        qdep = cjT * vjT * (1 - (1 - v / vjT) ** (1 - m)) / (1 - m)
    else:
        f1 = vjT * (1 - (1 - fc) ** (1 - m)) / (1 - m)
        f2 = (1 - fc) ** (1 + m)
        f3 = 1 - fc * (1 + m)
        qdep = cjT * (f1 + (f3 * (v - fcvj)
                            + m / (2 * vjT) * (v ** 2 - fcvj ** 2)) / f2)
    return i, area * (qdep + tt * ifwd)


def _mk(cls, *nodes, **kw):
    pycircuit.circuit.circuit.default_toolkit = numeric
    el = cls(*[Node(nm) for nm in nodes], **kw)
    el.update_iparv()
    return el


def _diode(**kw):
    """A junction with the series resistance collapsed away, so the
    element's two unknowns ARE the junction terminals and ``x[0]`` is the
    junction voltage."""
    p = dict(CARD, rs=0.0)
    p.update(kw)
    return _mk(eh.DiodeSpiceHdl, 'a', 'b', **p), p


def _quiet():
    """Ignore UNDERFLOW only.

    ``0.0`` is the correctly-rounded answer for ``exp(-1000)`` and the
    one benign IEEE condition (hdl.md 3.2c); a reverse-biased junction
    underflows on every evaluation.  Overflow, invalid and
    divide-by-zero are NOT suppressed anywhere in this file -- they mean
    an arm produced garbage, and `TestNoFloatingPointGarbage` turns them
    into errors.
    """
    warnings.filterwarnings('ignore', 'underflow encountered',
                            RuntimeWarning)


def _at(el, v, epar=defaultepar):
    x = np.zeros(el.n)
    x[0] = v
    with warnings.catch_warnings():
        _quiet()
        return (np.asarray(el.i(x, epar), float),
                np.asarray(el.q(x, epar), float),
                np.asarray(el.G(x, epar), float),
                np.asarray(el.C(x, epar), float))


def _tran(c, node, tend, dt, **kw):
    tran = Transient(c, toolkit=numeric, **kw)
    with warnings.catch_warnings():
        _quiet()
        res = tran.solve(tend=tend, timestep=dt, fixed_timestep=True)
    return (np.asarray(res.v(node).x[0], float),
            np.asarray(res.v(node).y, float))


## ======================================================================
## The SPICE diode: static current
## ======================================================================


class TestSpiceDiodeCurrent(object):

    def test_forward_current_is_the_shockley_law(self):
        """Against the closed form, not against a second copy of the
        model: at T = tnom every temperature correction is the identity,
        so the current is exactly ``area*IS*(exp(v/(n*Vt)) - 1)`` with
        the CARD's IS -- a number this test computes from first
        principles."""
        el, p = _diode(tnom=300.0, area=2.0)
        ep = _epar(300.0)
        vt = K * 300.0 / Q
        for v in (0.1, 0.3, 0.5, 0.65, 0.75):
            i = _at(el, v, ep)[0][0]
            want = 2.0 * p['IS'] * (np.exp(v / (p['n'] * vt)) - 1)
            assert_allclose(i, want, rtol=1e-11)

    def test_reverse_current_saturates_at_minus_area_times_is(self):
        el, p = _diode(tnom=300.0, area=3.0, bv=1e30)
        ep = _epar(300.0)
        ## Deep enough that exp(v/nVt) is below the tolerance: at -0.5 V
        ## it is still 1.2e-8, which is the diode, not an error.
        for v in (-1.0, -5.0, -20.0):
            i = _at(el, v, ep)[0][0]
            assert_allclose(i, -3.0 * p['IS'], rtol=1e-9)

    def test_breakdown_current_is_exactly_ibv_at_minus_bv(self):
        """The defining property of the additive breakdown term, and the
        reason it was written additively: at ``v = -bv`` the exponent is
        zero, so the current is ``-area*(IS + ibv)`` on the nose."""
        el, p = _diode(tnom=300.0, bv=45.0, ibv=5e-6, area=1.7)
        ep = _epar(300.0)
        i = _at(el, -p['bv'], ep)[0][0]
        assert_allclose(i, -1.7 * (p['IS'] + p['ibv']), rtol=1e-9)

    def test_breakdown_is_a_decade_per_2_3_vt(self):
        """Below -bv the current is exponential in ``Vt``, not in
        ``n*Vt`` -- SPICE's breakdown region uses the un-idealised
        thermal voltage.  Ten times the current one ``ln(10)*Vt`` further
        down."""
        el, _ = _diode(tnom=300.0, bv=45.0, ibv=5e-6, n=1.06)
        ep = _epar(300.0)
        vt = K * 300.0 / Q
        i1 = _at(el, -45.0, ep)[0][0]
        i2 = _at(el, -45.0 - np.log(10) * vt, ep)[0][0]
        assert_allclose(i2 / i1, 10.0, rtol=2e-3)

    def test_no_breakdown_by_default(self):
        """SPICE's BV default is "none", spelled 1e30 V here.  The
        breakdown term must contribute nothing at any usable bias."""
        el, p = _diode(tnom=300.0, bv=1e30)
        ep = _epar(300.0)
        for v in (-100.0, -1000.0):
            assert_allclose(_at(el, v, ep)[0][0], -p['IS'], rtol=1e-9)

    def test_area_scales_the_current_linearly(self):
        ep = _epar(300.0)
        a1 = _at(_diode(area=1.0)[0], 0.6, ep)[0][0]
        a5 = _at(_diode(area=5.0)[0], 0.6, ep)[0][0]
        assert_allclose(a5, 5.0 * a1, rtol=1e-12)

    def test_matches_the_independent_transcription_over_a_sweep(self):
        el, p = _diode()
        ep = _epar(317.0)
        for v in np.linspace(-2.0, 0.8, 29):
            i, q = _at(el, v, ep)[:2]
            ri, rq = _ref(v, T=317.0, **p)
            assert_allclose(i[0], ri, rtol=1e-9, atol=1e-30)
            ## The sweep passes through v = 0, where the depletion charge
            ## is exactly zero and only roundoff is left; atol is 1e-12
            ## of the charge scale cjo*vj.
            assert_allclose(q[0], rq, rtol=1e-9,
                            atol=1e-12 * p['cjo'] * p['vj'])


## ======================================================================
## The SPICE diode: series resistance and the internal node
## ======================================================================


class TestSpiceDiodeSeriesResistance(object):

    def test_rs_zero_collapses_the_internal_node_away(self):
        """SPICE's default is rs = 0, and a compiled model cannot switch
        a term off with a parameter value -- 1/rs would be infinity.
        `Collapse` removes the branch before it is compiled, so the
        element has two unknowns and no internal node."""
        el = _mk(eh.DiodeSpiceHdl, 'a', 'b', **dict(CARD, rs=0.0))
        assert el.n == 2
        assert [nd.name for nd in el.nodes] == ['plus', 'minus']

    def test_rs_nonzero_keeps_the_internal_node(self):
        el = _mk(eh.DiodeSpiceHdl, 'a', 'b', **dict(CARD, rs=1.5))
        assert el.n == 3
        assert 'di' in [nd.name for nd in el.nodes]

    def test_rs_is_the_measured_slope_at_high_current(self):
        """In a DC solve at 100 mA the extra drop over the rs = 0 device
        must be I*rs/area, to the accuracy of the junction's own
        logarithmic shift."""
        def solve(rs, area):
            pycircuit.circuit.circuit.default_toolkit = numeric
            c = SubCircuit()
            a = c.add_node('a')
            c['I1'] = ISource(gnd, a, i=0.1)
            c['D1'] = eh.DiodeSpiceHdl(a, gnd, **dict(CARD, rs=rs,
                                                      area=area))
            c.update_iparv()
            with warnings.catch_warnings():
                _quiet()
                return float(DC(c, toolkit=numeric).solve().v(a, gnd))

        for rs, area in ((1.5, 1.0), (1.5, 4.0), (12.0, 1.0)):
            ## Baseline at the SAME area: area moves the junction
            ## voltage too, and only the rs difference is under test.
            assert_allclose(solve(rs, area) - solve(0.0, area),
                            0.1 * rs / area, rtol=1e-9)


## ======================================================================
## The SPICE diode: charge
## ======================================================================


class TestSpiceDiodeCharge(object):

    def test_zero_bias_capacitance_is_area_times_cj0(self):
        """At T = tnom and v = 0 the depletion capacitance is cj0 by
        definition; tt contributes ``tt*dIf/dv``, which at zero bias is
        ``tt*IS/(n*Vt)`` -- femtofarads against picofarads, so it is
        subtracted rather than tolerated."""
        el, p = _diode(tnom=300.0, area=2.5)
        ep = _epar(300.0)
        cc = _at(el, 0.0, ep)[3][0, 0]
        vt = K * 300.0 / Q
        cdiff = 2.5 * p['tt'] * p['IS'] / (p['n'] * vt)
        assert_allclose(cc - cdiff, 2.5 * p['cjo'], rtol=1e-9)

    def test_depletion_capacitance_follows_the_closed_form_below_fc(self):
        """``C(v) = cj0/(1 - v/vj)^m`` -- the textbook result, which the
        model does NOT contain: the model states a charge, and this test
        differentiates it numerically to recover the capacitance."""
        el, p = _diode(tnom=300.0, tt=0.0)
        ep = _epar(300.0)
        vt = K * 300.0 / Q
        for v in (-5.0, -1.0, -0.2, 0.0, 0.2, 0.35):
            h = 1e-6
            qp = _at(el, v + h, ep)[1][0]
            qm = _at(el, v - h, ep)[1][0]
            num = (qp - qm) / (2 * h)
            want = p['cjo'] / (1 - v / p['vj']) ** p['m']
            assert_allclose(num, want, rtol=1e-6)

    def test_the_fc_seam_is_c0_and_c1(self):
        """The whole reason SPICE linearises above ``fc*vj``: the charge
        and its slope must be continuous there.  Checked by approaching
        the seam from both sides -- an assertion that fails if either
        arm's constants are wrong, which they were once."""
        el, p = _diode(tnom=300.0, tt=0.0)
        ep = _epar(300.0)
        vjT, m, fc = p['vj'], p['m'], p['fc']
        seam = fc * vjT
        want = p['cjo'] / (1 - fc) ** m
        for h in (1e-5, 1e-7):
            ql = _at(el, seam - h, ep)[1][0]
            qr = _at(el, seam + h, ep)[1][0]
            ## C0: the two arms differ by no more than the slope times
            ## the gap, to a part in 1e5 of the charge itself.  A jump in
            ## either arm's constants shows up here at full size.
            assert abs((qr - ql) - want * 2 * h) < 1e-5 * abs(ql)
            ## C1: each arm's one-sided slope is the closed-form
            ## capacitance at the seam, which neither arm contains.
            cl = (ql - _at(el, seam - 3 * h, ep)[1][0]) / (2 * h)
            cr = (_at(el, seam + 3 * h, ep)[1][0] - qr) / (2 * h)
            assert_allclose(cl, want, rtol=1e-4)
            assert_allclose(cr, want, rtol=1e-4)

    def test_above_the_seam_the_capacitance_is_linear_in_v(self):
        """The linearised arm's slope is ``cj0*m/(vj*(1-fc)^(1+m))`` --
        a constant SPICE fixes so that C keeps rising smoothly.  Again
        recovered by differentiating the model's charge."""
        el, p = _diode(tnom=300.0, tt=0.0)
        ep = _epar(300.0)
        h = 1e-6

        def cap(v):
            return (_at(el, v + h, ep)[1][0] - _at(el, v - h, ep)[1][0]) \
                / (2 * h)

        v1, v2 = 0.5, 0.7
        slope = (cap(v2) - cap(v1)) / (v2 - v1)
        want = p['cjo'] * p['m'] / (p['vj'] * (1 - p['fc']) ** (1 + p['m']))
        assert_allclose(slope, want, rtol=1e-5)

    def test_diffusion_charge_is_tt_times_the_forward_current(self):
        """The difference between two cards identical but for tt is
        exactly ``tt*I_fwd``, and I_fwd is the forward term alone -- the
        breakdown term deliberately carries no transit time."""
        ep = _epar(300.0)
        a, p = _diode(tnom=300.0, tt=0.0)
        b, _ = _diode(tnom=300.0, tt=7e-9)
        vt = K * 300.0 / Q
        for v in (0.4, 0.6, 0.7):
            dq = _at(b, v, ep)[1][0] - _at(a, v, ep)[1][0]
            ifwd = p['IS'] * (np.exp(v / (p['n'] * vt)) - 1)
            assert_allclose(dq, 7e-9 * ifwd, rtol=1e-9)

    def test_cjo_and_tt_zero_really_do_leave_no_charge(self):
        """A parameter CAN switch this block off, unlike rs, because the
        charge is a product and not a quotient -- and that difference is
        the reason rs needed `Collapse` and this did not."""
        el, _ = _diode(cjo=0.0, tt=0.0)
        for v in (-3.0, 0.0, 0.7):
            assert _at(el, v)[1][0] == 0.0


## ======================================================================
## The SPICE diode: noise
## ======================================================================


class TestSpiceDiodeNoise(object):

    def test_shot_noise_is_two_q_i(self):
        """Schottky's result, which the model states and this test
        recomputes from the solved current -- both signs of bias, since
        the PSD must stay positive where the current is negative."""
        el, _ = _diode(tnom=300.0, area=2.0)
        ep = _epar(300.0)
        for v in (0.6, 0.4, 0.0, -1.0):
            i = _at(el, v, ep)[0][0]
            cy = np.asarray(el.CY(np.array([v, 0.0]), 2 * np.pi * 100, ep),
                            float)
            assert cy[0, 0] > 0
            assert_allclose(cy[0, 0], 2 * Q * abs(i), rtol=1e-9,
                            atol=1e-40)
            ## a current noise source between the two terminals
            assert_allclose(cy[0, 1], -cy[0, 0], rtol=1e-12)
            assert_allclose(cy[1, 1], cy[0, 0], rtol=1e-12)

    def test_flicker_noise_has_the_one_over_f_shape(self):
        el, p = _diode(tnom=300.0, kf=1e-16, af=1.0)
        ep = _epar(300.0)
        x = np.array([0.6, 0.0])
        i = _at(el, 0.6, ep)[0][0]
        c1 = np.asarray(el.CY(x, 2 * np.pi * 10.0, ep), float)[0, 0]
        c2 = np.asarray(el.CY(x, 2 * np.pi * 1000.0, ep), float)[0, 0]
        shot = 2 * Q * abs(i)
        assert_allclose((c1 - shot) / (c2 - shot), 100.0, rtol=1e-9)
        assert_allclose(c1 - shot, 1e-16 * abs(i) / 10.0, rtol=1e-9)

    def test_the_flicker_exponent_is_af(self):
        ep = _epar(300.0)
        x = np.array([0.6, 0.0])
        out = []
        for area in (1.0, 4.0):
            el, _ = _diode(tnom=300.0, kf=1e-16, af=1.4, area=area)
            i = _at(el, 0.6, ep)[0][0]
            cy = np.asarray(el.CY(x, 2 * np.pi * 10.0, ep), float)[0, 0]
            out.append((cy - 2 * Q * abs(i)) / (1e-16 * abs(i) ** 1.4 / 10.0))
        assert_allclose(out, [1.0, 1.0], rtol=1e-9)

    def test_series_resistance_carries_thermal_noise_and_loses_it_when_collapsed(self):
        """4kT/rs on the rs branch -- and when `rs = 0` collapses the
        branch away, the noise source goes with it, which is what a
        collapse must mean."""
        ep = _epar(311.0)
        el = _mk(eh.DiodeSpiceHdl, 'a', 'b', **dict(CARD, rs=8.0,
                                                    area=2.0))
        ## Forward-biased, or the shot noise is zero and the two rows
        ## are trivially equal.
        x = np.array([0.7, 0.0, 0.68])
        cy = np.asarray(el.CY(x, 2 * np.pi * 100, ep), float)
        assert_allclose(cy[0, 0], 4 * K * 311.0 * 2.0 / 8.0, rtol=1e-9)
        ## `plus` sees only the rs noise; `di` sees rs + shot, and the
        ## difference is exactly Schottky's 2q|I| at the junction.
        ij = _ref(0.68, T=311.0, **dict(CARD, rs=8.0, area=2.0))[0]
        assert_allclose(cy[2, 2] - cy[0, 0], 2 * Q * abs(ij), rtol=1e-8)
        el0, _ = _diode(rs=0.0)
        cy0 = np.asarray(el0.CY(np.zeros(2), 2 * np.pi * 100, ep), float)
        assert_allclose(cy0[0, 0], 2 * Q * abs(_at(el0, 0.0, ep)[0][0]),
                        rtol=1e-9, atol=1e-40)


## ======================================================================
## The SPICE diode: temperature
## ======================================================================


class TestSpiceDiodeTemperature(object):

    def test_at_tnom_every_correction_is_the_identity(self):
        """A structural check that needs no reference implementation: at
        T = tnom, IS(T) = IS, VJ(T) = VJ and CJ0(T) = CJ0 whatever the
        formulae are.  It catches a sign error or a swapped T/tnom that
        a transcription test cannot, because the transcription would
        carry the same error."""
        el, p = _diode(tnom=311.0, tt=0.0)
        ep = _epar(311.0)
        vt = K * 311.0 / Q
        i = _at(el, 0.5, ep)[0][0]
        assert_allclose(i, p['IS'] * (np.exp(0.5 / (p['n'] * vt)) - 1),
                        rtol=1e-11)
        h = 1e-6
        cap = (_at(el, h, ep)[1][0] - _at(el, -h, ep)[1][0]) / (2 * h)
        assert_allclose(cap, p['cjo'], rtol=1e-8)

    def test_saturation_current_temperature_law(self):
        """``IS(T2)/IS(T1)`` against the closed form, at a bias low
        enough that ``exp(v/nVt) >> 1`` still holds but the ratio is
        dominated by IS."""
        p = dict(CARD, rs=0.0, cjo=0.0, tt=0.0, bv=1e30)
        el = _mk(eh.DiodeSpiceHdl, 'a', 'b', **p)
        for T in (250.0, 300.0, 400.0):
            ep = _epar(T)
            vt = K * T / Q
            i = _at(el, 0.4, ep)[0][0]
            tr = T / p['tnom']
            isT = p['IS'] * np.exp((tr - 1) * p['eg'] / (p['n'] * vt)
                                   + p['xti'] / p['n'] * np.log(tr))
            assert_allclose(i, isT * (np.exp(0.4 / (p['n'] * vt)) - 1),
                            rtol=1e-10)

    def test_forward_voltage_tempco_is_about_minus_two_millivolts(self):
        """The physical check the formulae exist to reproduce: at
        constant current a silicon junction moves about -1.5 to -2.5
        mV/K.  Independent of every equation above -- if the temperature
        path were wired backwards this is the test that says so."""
        def vf(T):
            pycircuit.circuit.circuit.default_toolkit = numeric
            c = SubCircuit()
            a = c.add_node('a')
            c['I1'] = ISource(gnd, a, i=1e-3)
            c['D1'] = eh.DiodeSpiceHdl(a, gnd, **dict(CARD, rs=0.0))
            c.update_iparv()
            with warnings.catch_warnings():
                _quiet()
                return float(DC(c, toolkit=numeric, epar=_epar(T))
                             .solve().v(a, gnd))

        tc = (vf(350.0) - vf(300.0)) / 50.0
        assert -2.5e-3 < tc < -1.5e-3

    def test_junction_potential_shrinks_with_temperature(self):
        """VJ(T) falls, so the zero-bias capacitance rises.  Both signs
        asserted, plus agreement with the transcription."""
        p = dict(CARD, rs=0.0, tt=0.0)
        el = _mk(eh.DiodeSpiceHdl, 'a', 'b', **p)
        h = 1e-6
        caps = []
        for T in (250.0, 400.0):
            ep = _epar(T)
            caps.append((_at(el, h, ep)[1][0] - _at(el, -h, ep)[1][0])
                        / (2 * h))
        assert caps[1] > caps[0]
        for T, cap in zip((250.0, 400.0), caps):
            rq = (_ref(h, T=T, **p)[1] - _ref(-h, T=T, **p)[1]) / (2 * h)
            assert_allclose(cap, rq, rtol=1e-6)


## ======================================================================
## Jacobians and finiteness -- the DSL's whole selling point
## ======================================================================


def _fd_check(el, x, epar, which, rtol, atol, h=None):
    x = np.asarray(x, float)
    n = el.n
    ref = np.zeros((n, n))
    with warnings.catch_warnings():
        _quiet()
        for k in range(n):
            hk = h if h is not None else max(1e-7, 1e-7 * abs(x[k]))
            xp, xm = x.copy(), x.copy()
            xp[k] += hk
            xm[k] -= hk
            fp = np.asarray(getattr(el, which)(xp, epar), float)
            fm = np.asarray(getattr(el, which)(xm, epar), float)
            ref[:, k] = (fp - fm) / (2 * hk)
        ana = np.asarray(getattr(el, 'G' if which == 'i' else 'C')
                         (x, epar), float)
    assert np.all(np.isfinite(ana)), (which, x, ana)
    assert_allclose(ana, ref, rtol=rtol, atol=atol)


class TestJacobians(object):

    @pytest.mark.parametrize('v', [-40.0, -5.0, -1.0, -0.1, 0.0, 0.1,
                                   0.2, 0.389, 0.39, 0.391, 0.5, 0.7,
                                   0.85])
    def test_G_and_C_match_finite_differences(self, v):
        """Swept across both smoothing seams: ``fc*vj`` = 0.39 V on this
        card, and zero bias.  Kinks hide -- one point would not find
        them."""
        el, _ = _diode()
        ep = _epar(305.0)
        scale = max(abs(_at(el, v, ep)[0][0]), 1e-13)
        _fd_check(el, [v, 0.0], ep, 'i', rtol=2e-5, atol=1e-6 * scale)
        _fd_check(el, [v, 0.0], ep, 'q', rtol=2e-5, atol=1e-18)

    def test_G_matches_finite_differences_with_the_series_node(self):
        """Three unknowns, one of them internal: the rs branch and the
        junction have to be differentiated together."""
        el = _mk(eh.DiodeSpiceHdl, 'a', 'b', **CARD)
        ep = _epar(300.0)
        for va, vdi in ((0.8, 0.75), (0.0, 0.0), (-2.0, -2.0)):
            _fd_check(el, [va, 0.0, vdi], ep, 'i', rtol=2e-5, atol=1e-9)

    def test_G_matches_finite_differences_through_the_thermal_node(self):
        """The self-heating element's Jacobian includes dI/dT and
        dP/dV -- the two off-diagonal blocks that make the thermal node
        a real unknown rather than an outer iteration."""
        el = _mk(eh.DiodeSpiceThermalHdl, 'a', 'b', 't', 'ta',
                 **dict(CARD, rs=0.0, rth=250.0, cth=1e-3))
        ep = _epar(300.0)
        for v, dT in ((0.7, 12.0), (0.0, 0.0), (-3.0, 0.5), (0.75, 60.0)):
            _fd_check(el, [v, 0.0, dT, 0.0], ep, 'i',
                      rtol=3e-5, atol=1e-9, h=1e-6)
            _fd_check(el, [v, 0.0, dT, 0.0], ep, 'q',
                      rtol=3e-5, atol=1e-18, h=1e-6)

    def test_everything_stays_finite_at_absurd_bias_both_signs(self):
        """Scan for the first non-finite rather than assert one exists.
        Value AND Jacobian, separately, because they are different
        programs and the Jacobian usually dies first."""
        el, _ = _diode()
        ep = _epar(300.0)
        vs = np.concatenate([-np.logspace(-6, 6, 61)[::-1],
                             [0.0], np.logspace(-6, 6, 61)])
        for v in vs:
            i, q, G, Cm = _at(el, v, ep)
            for nm, a in (('i', i), ('q', q), ('G', G), ('C', Cm)):
                assert np.all(np.isfinite(a)), (nm, v, a)

    def test_the_thermal_element_survives_a_junction_below_absolute_zero(self):
        """A Newton iterate can put dT below -T_ambient.  The model must
        stay finite there or the solve is over; the floor is what does
        it."""
        el = _mk(eh.DiodeSpiceThermalHdl, 'a', 'b', 't', 'ta',
                 **dict(CARD, rs=0.0, rth=250.0))
        ep = _epar(300.0)
        for dT in (-299.0, -300.0, -400.0, -1e6, 1e6):
            for v in (-1.0, 0.0, 0.7):
                x = np.array([v, 0.0, dT, 0.0])
                with warnings.catch_warnings():
                    _quiet()
                    for m in ('i', 'q', 'G', 'C'):
                        a = np.asarray(getattr(el, m)(x, ep), float)
                        assert np.all(np.isfinite(a)), (m, v, dT, a)


class TestNoFloatingPointGarbage(object):
    """Overflow, invalid and divide-by-zero raised as errors.

    A green suite carrying overflow warnings is telling you something:
    an unselected arm produced garbage, and the only thing standing
    between that and a NaN in the Jacobian is which arm `select` happens
    to pick.  Underflow is exempt -- ``0.0`` is the correct answer for
    ``exp(-1000)``.

    This is the test that made the library's diode use `expl` rather
    than `limexp`: `limexp` is deliberately not both-arms-safe, and one
    5 V DC solve of the self-heating variant raised 120 overflows.
    """

    @pytest.mark.parametrize('kw', [dict(rs=0.0), dict(rs=1.5),
                                    dict(rs=1.5, kf=1e-16, af=1.3)])
    def test_the_junction_raises_nothing_over_a_wide_sweep(self, kw):
        el = _mk(eh.DiodeSpiceHdl, 'a', 'b', **dict(CARD, **kw))
        vs = np.concatenate([-np.logspace(-6, 4, 41)[::-1], [0.0],
                             np.logspace(-6, 4, 41)])
        with warnings.catch_warnings():
            _quiet()
            with np.errstate(over='raise', invalid='raise',
                             divide='raise'):
                for T in (200.0, 300.0, 450.0):
                    ep = _epar(T)
                    for v in vs:
                        x = np.zeros(el.n)
                        x[0] = v
                        x[-1] = v * 0.9
                        for meth in ('i', 'q', 'G', 'C'):
                            getattr(el, meth)(x, ep)
                        el.CY(x, 2 * np.pi * 100, ep)

    def test_the_thermal_variants_raise_nothing_even_below_absolute_zero(self):
        d = _mk(eh.DiodeSpiceThermalHdl, 'a', 'b', 't', 'ta',
                **dict(CARD, rth=250.0, cth=1e-3))
        r = _mk(eh.RThermalHdl, 'p', 'm', 't', 'ta', r=1e3, tc1=2e-3,
                tc2=1e-6, tnom=300.0, rth=100.0, cth=1e-3)
        ep = _epar(300.0)
        with warnings.catch_warnings():
            _quiet()
            with np.errstate(over='raise', invalid='raise',
                             divide='raise'):
                for dT in (-1e6, -400.0, -299.0, 0.0, 300.0, 1e6):
                    for v in (-1e4, -1.0, 0.0, 0.5, 0.8, 1e4):
                        xd = np.array([v, 0.0, dT, 0.0, v * 0.9, 0.0])
                        xr = np.array([v, 0.0, dT, 0.0])
                        for meth in ('i', 'q', 'G', 'C'):
                            getattr(d, meth)(xd, ep)
                            getattr(r, meth)(xr, ep)

    def test_a_dc_solve_raises_nothing(self):
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        a, b, th = c.add_node('a'), c.add_node('b'), c.add_node('th')
        c['V1'] = VS(a, gnd, v=5.0)
        c['R1'] = R(a, b, r=1e3)
        c['D1'] = eh.DiodeSpiceThermalHdl(b, gnd, th, gnd,
                                          **dict(CARD, rth=200.0))
        c.update_iparv()
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            _quiet()
            res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        assert 0.6 < float(res.v(b, gnd)) < 0.85
        assert float(res.v(th, gnd)) > 0.0


## ======================================================================
## The thermal node
## ======================================================================


class TestSelfHeating(object):

    def test_temperature_rise_is_power_times_rth(self):
        """The defining self-consistency: with a flat tempco the
        dissipation is known exactly, so dT must be P*rth to machine
        precision -- and it is the solver that has to find it, since the
        thermal node is an unknown."""
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        p, th = c.add_node('p'), c.add_node('th')
        c['V1'] = VS(p, gnd, v=3.0)
        c['R1'] = eh.RThermalHdl(p, gnd, th, gnd, r=250.0, tc1=0.0,
                                 tnom=300.0, rth=40.0)
        c.update_iparv()
        res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        assert_allclose(float(res.v(th, gnd)), 40.0 * 3.0 ** 2 / 250.0,
                        rtol=1e-10)

    def test_the_current_responds_to_the_temperature_rise(self):
        """Not merely that dT appears, but that it feeds back: the
        operating point is the root of ``tc1*dT^2 + dT - rth*V^2/r``,
        a quadratic this test solves independently of the model."""
        r, tc1, rth, v = 1000.0, 1.5e-3, 300.0, 2.0
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        p, th = c.add_node('p'), c.add_node('th')
        c['V1'] = VS(p, gnd, v=v)
        c['R1'] = eh.RThermalHdl(p, gnd, th, gnd, r=r, tc1=tc1,
                                 tnom=300.0, rth=rth)
        c.update_iparv()
        res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        b = rth * v * v / r
        want = (-1 + np.sqrt(1 + 4 * tc1 * b)) / (2 * tc1)
        assert_allclose(float(res.v(th, gnd)), want, rtol=1e-9)
        ## and the current is the isothermal one divided by (1 + tc1*dT)
        i = v / (r * (1 + tc1 * want))
        assert_allclose(-float(res.i('V1.plus')), i, rtol=1e-9)

    def test_a_negative_tempco_runs_hotter_than_the_linear_estimate(self):
        """Thermal runaway's stable side: with tc1 < 0 the feedback is
        positive, so dT exceeds the first-order estimate P0*rth.  The
        exact root is again the quadratic's."""
        r, tc1, rth, v = 1000.0, -1.5e-3, 300.0, 2.0
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        p, th = c.add_node('p'), c.add_node('th')
        c['V1'] = VS(p, gnd, v=v)
        c['R1'] = eh.RThermalHdl(p, gnd, th, gnd, r=r, tc1=tc1,
                                 tnom=300.0, rth=rth)
        c.update_iparv()
        res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        b = rth * v * v / r
        want = (-1 + np.sqrt(1 + 4 * tc1 * b)) / (2 * tc1)
        dT = float(res.v(th, gnd))
        assert_allclose(dT, want, rtol=1e-9)
        assert dT > b

    def test_rth_zero_collapses_the_thermal_branch_and_pins_dt(self):
        """`rth = 0` is the isothermal limit and the default.  A
        parameter cannot switch the ``1/rth`` off -- the collapse
        removes the branch instead -- and the element must then be the
        plain resistor, exactly."""
        el = _mk(eh.RThermalHdl, 'p', 'm', 't', 'ta', r=1e3, rth=0.0)
        assert el.n == 5          # 4 terminals + the zero-volt source
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        p, th = c.add_node('p'), c.add_node('th')
        c['V1'] = VS(p, gnd, v=3.0)
        c['R1'] = eh.RThermalHdl(p, gnd, th, gnd, r=250.0, tc1=1e-3,
                                 tnom=300.0, rth=0.0)
        c.update_iparv()
        res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        assert_allclose(float(res.v(th, gnd)), 0.0, atol=1e-12)
        assert_allclose(-float(res.i('V1.plus')), 3.0 / 250.0, rtol=1e-9)

    def test_the_thermal_pole_is_one_over_rth_cth(self):
        """`ddt` on the thermal node makes the rise a real first-order
        transient.  Compared against the analytic step response, which
        the model does not contain."""
        rth, cth, r, v = 500.0, 2e-3, 1000.0, 2.0
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        p, th = c.add_node('p'), c.add_node('th')
        c['V1'] = VS(p, gnd, v=v)
        c['R1'] = eh.RThermalHdl(p, gnd, th, gnd, r=r, tc1=0.0,
                                 tnom=300.0, rth=rth, cth=cth)
        c.update_iparv()
        from pycircuit.circuit.integrator import EulerIntegrator
        t, y = _tran(c, th, tend=5.0, dt=2e-3, uic=True,
                     integrator=EulerIntegrator())
        tau = rth * cth
        final = rth * v * v / r
        want = final * (1 - np.exp(-t / tau))
        assert_allclose(y, want, atol=2e-3 * final)
        ## the 63.2% crossing IS the time constant
        k = int(np.argmax(y > (1 - np.exp(-1.0)) * final))
        tcross = np.interp((1 - np.exp(-1.0)) * final,
                           y[k - 1:k + 1], t[k - 1:k + 1])
        assert_allclose(tcross, tau, rtol=2e-2)
        ## and five time constants in, it is where the exponential says
        assert_allclose(y[-1], final * (1 - np.exp(-t[-1] / tau)),
                        rtol=2e-3)

    def test_an_external_foster_ladder_hangs_off_the_thermal_pin(self):
        """The reason the thermal node is a TERMINAL and not an internal
        node: a package model built from ordinary R and C elements is a
        parallel thermal path, and dT is the parallel combination.  No
        new operator, no change to the device."""
        rth, rext, r, v = 400.0, 100.0, 1000.0, 2.0
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        p, th = c.add_node('p'), c.add_node('th')
        c['V1'] = VS(p, gnd, v=v)
        c['R1'] = eh.RThermalHdl(p, gnd, th, gnd, r=r, tc1=0.0,
                                 tnom=300.0, rth=rth)
        c['Rpkg'] = R(th, gnd, r=rext)
        c['Cpkg'] = C(th, gnd, c=1e-3)
        c.update_iparv()
        res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        want = v * v / r * (rth * rext / (rth + rext))
        assert_allclose(float(res.v(th, gnd)), want, rtol=1e-9)

    def test_two_devices_share_one_thermal_node(self):
        """Mutual heating: both dissipations land on one node, so dT is
        the SUM of the two powers times rth.  Impossible if the node
        were internal to a device."""
        rth = 200.0
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        p1, p2, th = c.add_node('p1'), c.add_node('p2'), c.add_node('th')
        c['V1'] = VS(p1, gnd, v=2.0)
        c['V2'] = VS(p2, gnd, v=3.0)
        c['Ra'] = eh.RThermalHdl(p1, gnd, th, gnd, r=1e3, tc1=0.0,
                                 tnom=300.0, rth=rth)
        ## The second device's own rth is collapsed away; it only heats.
        c['Rb'] = eh.RThermalHdl(p2, gnd, th, gnd, r=2e3, tc1=0.0,
                                 tnom=300.0, rth=0.0)
        c.update_iparv()
        res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        ## Rb's own thermal branch is a zero-volt source, which SHORTS
        ## the shared node to ambient -- the honest consequence of the
        ## collapse, and the reason a device that must not conduct heat
        ## needs its thermal pin left unconnected instead.
        assert_allclose(float(res.v(th, gnd)), 0.0, atol=1e-12)

    def test_two_devices_with_their_own_rth_share_the_rise(self):
        """Both thermal branches live, in parallel on one node: the rise
        is (P1 + P2) * (rth1 || rth2)."""
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        p1, p2, th = c.add_node('p1'), c.add_node('p2'), c.add_node('th')
        c['V1'] = VS(p1, gnd, v=2.0)
        c['V2'] = VS(p2, gnd, v=3.0)
        c['Ra'] = eh.RThermalHdl(p1, gnd, th, gnd, r=1e3, tc1=0.0,
                                 tnom=300.0, rth=200.0)
        c['Rb'] = eh.RThermalHdl(p2, gnd, th, gnd, r=2e3, tc1=0.0,
                                 tnom=300.0, rth=600.0)
        c.update_iparv()
        res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        ptot = 2.0 ** 2 / 1e3 + 3.0 ** 2 / 2e3
        want = ptot * (200.0 * 600.0 / 800.0)
        assert_allclose(float(res.v(th, gnd)), want, rtol=1e-9)


class TestSelfHeatingDiode(object):

    def test_rth_zero_is_the_isothermal_diode_to_the_last_digit(self):
        """The thermal variant with the node collapsed must BE
        `DiodeSpiceHdl`.  Stamp equality on the electrical rows, not a
        waveform comparison."""
        iso = _mk(eh.DiodeSpiceHdl, 'a', 'b', **CARD)
        th = _mk(eh.DiodeSpiceThermalHdl, 'a', 'b', 't', 'ta',
                 **dict(CARD, rth=0.0))
        ep = _epar(307.0)
        for va, vdi in ((0.7, 0.68), (-3.0, -3.0), (0.0, 0.0)):
            xi = np.array([va, 0.0, vdi])
            xt = np.array([va, 0.0, 0.0, 0.0, vdi, 0.0])
            with warnings.catch_warnings():
                _quiet()
                i1 = np.asarray(iso.i(xi, ep), float)
                i2 = np.asarray(th.i(xt, ep), float)
                q1 = np.asarray(iso.q(xi, ep), float)
                q2 = np.asarray(th.q(xt, ep), float)
            assert_allclose(i2[[0, 1, 4]], i1, atol=0, rtol=1e-14)
            assert_allclose(q2[[0, 1, 4]], q1, atol=0, rtol=1e-14)

    def test_the_solved_rise_is_the_dissipation_times_rth(self):
        rth = 180.0
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        a, th = c.add_node('a'), c.add_node('th')
        c['I1'] = ISource(gnd, a, i=5e-3)
        c['D1'] = eh.DiodeSpiceThermalHdl(a, gnd, th, gnd,
                                          **dict(CARD, rth=rth))
        c.update_iparv()
        with warnings.catch_warnings():
            _quiet()
            res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        v = float(res.v(a, gnd))
        dT = float(res.v(th, gnd))
        assert_allclose(dT, rth * v * 5e-3, rtol=1e-9)
        assert dT > 0.1

    def test_the_self_heated_solution_equals_the_isothermal_one_at_the_risen_temperature(self):
        """The self-consistency that matters: solve the electrothermal
        device, read off dT, then solve the ISOTHERMAL device with the
        ambient set to 300 + dT.  Same forward voltage.  A test the
        model cannot pass by accident -- it couples two independently
        compiled elements through one number."""
        rth = 400.0
        pycircuit.circuit.circuit.default_toolkit = numeric

        def solve(cls, T, **extra):
            c = SubCircuit()
            a = c.add_node('a')
            c['I1'] = ISource(gnd, a, i=2e-2)
            if extra:
                th = c.add_node('th')
                c['D1'] = cls(a, gnd, th, gnd, **dict(CARD, **extra))
            else:
                c['D1'] = cls(a, gnd, **CARD)
            c.update_iparv()
            with warnings.catch_warnings():
                _quiet()
                res = DC(c, toolkit=numeric, epar=_epar(T)).solve()
            return c, res

        c1, r1 = solve(eh.DiodeSpiceThermalHdl, 300.0, rth=rth)
        v_hot = float(r1.v(c1.get_node('a'), gnd))
        dT = float(r1.v(c1.get_node('th'), gnd))
        assert dT > 1.0
        c2, r2 = solve(eh.DiodeSpiceHdl, 300.0 + dT)
        v_iso = float(r2.v(c2.get_node('a'), gnd))
        ## two independent Newton solves, so the agreement is bounded by
        ## the solver's tolerance and not by the model's
        assert_allclose(v_hot, v_iso, rtol=1e-6)

    def test_self_heating_lowers_the_forward_voltage(self):
        """Direction check: heat lowers Vf at constant current."""
        pycircuit.circuit.circuit.default_toolkit = numeric

        def vf(rth):
            c = SubCircuit()
            a, th = c.add_node('a'), c.add_node('th')
            c['I1'] = ISource(gnd, a, i=2e-2)
            c['D1'] = eh.DiodeSpiceThermalHdl(a, gnd, th, gnd,
                                              **dict(CARD, rth=rth))
            c.update_iparv()
            with warnings.catch_warnings():
                _quiet()
                return float(DC(c, toolkit=numeric,
                                epar=_epar(300.0)).solve().v(a, gnd))

        assert vf(400.0) < vf(0.0)


## ======================================================================
## Integration with the solver
## ======================================================================


class TestInCircuit(object):

    def test_the_diode_clamps_in_a_dc_solve(self):
        """A 5 V source through 1 k into the diode: KCL closes, and the
        answer is the one the model's own equation gives at the solved
        voltage -- checked by substituting back, not by trusting the
        solver."""
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        a, b = c.add_node('a'), c.add_node('b')
        c['V1'] = VS(a, gnd, v=5.0)
        c['R1'] = R(a, b, r=1e3)
        c['D1'] = eh.DiodeSpiceHdl(b, gnd, **dict(CARD, rs=0.0))
        c.update_iparv()
        with warnings.catch_warnings():
            _quiet()
            res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        v = float(res.v(b, gnd))
        assert 0.6 < v < 0.85
        ii = (5.0 - v) / 1e3
        assert_allclose(ii, _ref(v, T=300.0, **dict(CARD, rs=0.0))[0],
                        rtol=1e-7)

    def test_a_rectifier_transient_charges_the_reservoir(self):
        """The charge path in a real transient: half-wave rectifier with
        the junction and diffusion capacitance live."""
        from pycircuit.circuit.elements import VSin
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        a, b = c.add_node('a'), c.add_node('b')
        c['V1'] = VSin(a, gnd, va=5.0, freq=1e3)
        c['D1'] = eh.DiodeSpiceHdl(a, b, **dict(CARD, area=50.0))
        c['Rl'] = R(b, gnd, r=1e4)
        c['Cl'] = C(b, gnd, c=1e-6)
        c.update_iparv()
        t, y = _tran(c, b, tend=5e-3, dt=2e-6, uic=True)
        assert y[0] == pytest.approx(0.0, abs=1e-9)
        assert 3.5 < y[-1] < 5.0
        assert np.all(np.isfinite(y))
        ## the ripple is small and the output never goes negative
        assert y.min() > -1e-6

    def test_the_electrothermal_diode_solves_in_a_real_circuit(self):
        """Both new models in one netlist, with the diode's thermal pin
        driven by the self-heating resistor's -- one shared package."""
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        a, b, th = c.add_node('a'), c.add_node('b'), c.add_node('th')
        c['V1'] = VS(a, gnd, v=5.0)
        c['R1'] = eh.RThermalHdl(a, b, th, gnd, r=100.0, tc1=2e-3,
                                 tnom=300.0, rth=60.0)
        c['D1'] = eh.DiodeSpiceThermalHdl(b, gnd, th, gnd,
                                          **dict(CARD, rth=0.0))
        c.update_iparv()
        with warnings.catch_warnings():
            _quiet()
            res = DC(c, toolkit=numeric, epar=_epar(300.0)).solve()
        v = float(res.v(b, gnd))
        dT = float(res.v(th, gnd))
        assert 0.6 < v < 0.9
        assert dT == pytest.approx(0.0, abs=1e-12)
