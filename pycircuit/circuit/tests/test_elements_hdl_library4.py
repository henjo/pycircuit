# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The fourth Phase-6 library batch: the four-probe MOSFET, the opamp
macromodel, and the self-heating bipolar.

=========================  ==================================================
model                      what it is, and why it is here
=========================  ==================================================
`MosLevel1Hdl`             SPICE level 1 (Shichman-Hodges), and the **first
`MosLevel1PmosHdl`         production user of ``limit_together``.  SPICE's own
                           MOSFET limiting is four probes over four terminals
                           with a cycle in them, which the per-probe form
                           cannot compile at all -- section 1.9 asserts the
                           refusal on this very model
`OpAmpHdl`                 a Boyle-class opamp macromodel, parameterised by
                           the datasheet and tested against Boyle's own
                           closed-form relations
`GummelPoonNpnThermalHdl`  the self-heating bipolar, and thermal runaway with
                           its onset measured against the loop-gain criterion
=========================  ==================================================

**Every reference here is external to the model under test.**  In order
of how much they can catch:

* an **independent numpy transcription** of the Shichman-Hodges
  equations, written from Massobrio & Antognetti, *Semiconductor Device
  Modeling with SPICE*, 2nd ed., ch. 3 (Shichman and Hodges, *IEEE J.
  Solid-State Circuits* **3**, 285, 1968), living in `_m1_ref` below.  It
  shares no code with the DSL model;
* **relations derived by hand and written next to the number**, which is
  the strongest kind because a transcription error cannot satisfy them.
  The opamp's whole test set is of this kind: Boyle's ``GBW = gm/(2*pi*C)``,
  ``SR = I/C`` and ``f_dominant = GBW/A_ol`` are identities of the
  parameterisation, so they are asserted to parts in 1e4 off real AC and
  transient solves rather than to a band;
* the **loop-gain criterion for thermal runaway**, ``rth*dP/dT = 1``,
  solved on an independent transcription of the ideal transport BJT and
  compared against the ``rth`` at which the compiled model stops having a
  DC solution;
* **textbook limits of the channel noise** -- Nyquist at ``Vds = 0`` and
  ``(2/3)*4kT*gm`` in saturation -- which are properties of the
  Klaassen-Prins integral and not of this implementation;
* **finite differences** for every Jacobian claim, through
  `hdl.check_jacobians`, with the UNRESOLVED verdicts reported;
* **the intervention against its absence** for the limiter: the same
  model compiled with ``limiting='none'``, on the same circuit, on a
  plain Newton with ``gmin = 0``, counting Jacobian evaluations.
"""

import math
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pycircuit.circuit.circuit
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.analysis_ss import AC
from pycircuit.circuit.elements import VS, IS as ISRC, R, C
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.hdl import (Behavioural, Node, TEMP, check_jacobians,
                                   x_layout, explain, KBOLTZMANN, QELECTRON)

_KB = float(KBOLTZMANN)
_QE = float(QELECTRON)
_T0 = float(defaultepar.T)
_UT = _KB * _T0 / _QE


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    """Every element here is built by calling its class, which reads
    ``circuit.default_toolkit`` at construction."""
    old = pycircuit.circuit.circuit.default_toolkit
    pycircuit.circuit.circuit.default_toolkit = numeric
    yield
    pycircuit.circuit.circuit.default_toolkit = old


def _mk(cls, *nodes, **kw):
    el = cls(*[Node(nm) for nm in nodes], **kw)
    el.update_iparv()
    return el


## ======================================================================
## 0.  Cards, and the independent level-1 reference
## ======================================================================

#: A card with every block switched on, so that no test can pass by a
#: term being zero.  ``tnom`` is the simulator's own ambient, so the
#: temperature path is the identity everywhere except in the one test
#: that is about the temperature path.  Nothing here is fitted.
NMOS = dict(vto=0.75, kp=8e-5, gamma=0.55, phi=0.70, lambd=0.03,
            tox=2e-8, w=20e-6, l=2e-6, ld=1.5e-7,
            ## The three overlap densities are DELIBERATELY DIFFERENT.
            ## With `cgso == cgdo` a mutation that puts the drain overlap
            ## charge on the source branch is not a corruption at all,
            ## and the control silently passes -- measured, it did.
            cgso=3.5e-10, cgdo=2.5e-10, cgbo=2.0e-10,
            IS=1e-14, pb=0.85, cj=3e-4, cjsw=3e-10, mj=0.5, mjsw=0.33,
            fc=0.5, js=0.0, ad=1.6e-10, asrc=1.6e-10, pd=4.6e-5, ps=4.6e-5,
            rd=0.0, rs=0.0, rsh=0.0, nrd=2.0, nrs=2.0,
            kf=1e-27, af=1.0, tnom=_T0)

#: The same card with the two parasitic resistances on, so that the
#: internal nodes and their collapse are both exercised.
NMOS_R = dict(NMOS, rsh=25.0)

#: A minimal card: no body effect, no modulation, no charge, no
#: parasitics.  ``gamma = 0`` is a legitimate card and it is the one the
#: square-law asymptote tests use -- see `differentiable-numerics`, "the
#: dangerous value can be a PARAMETER".
NMOS_IDEAL = dict(vto=0.75, kp=8e-5, gamma=0.0, phi=0.70, lambd=0.0,
                  w=20e-6, l=2e-6, tnom=_T0)


def _m1_ref(vd, vg, vs, vb, card):
    """Shichman-Hodges, transcribed from Massobrio & Antognetti ch. 3.

    Written the way the textbook and SPICE write it -- pick the source by
    comparing the two diffusion potentials, take the threshold at THAT
    terminal's bulk bias, and return the current with the sign of the
    swap -- which is deliberately NOT how the model is written.  The
    model has one expression with no branch on the solution; if the two
    agree, the algebraic identity the model relies on is confirmed
    against the form it replaced.

    Returns ``(id_channel, i_bs, i_bd)`` in the n-channel convention.
    """
    p = dict(card)
    leff = p['l'] - 2.0 * p.get('ld', 0.0)
    beta = p['kp'] * p['w'] / leff
    ## The source is the lower diffusion; the sign is the swap.
    if vd >= vs:
        vsrc, vdrn, sgn = vs, vd, 1.0
    else:
        vsrc, vdrn, sgn = vd, vs, -1.0
    vbs = vb - vsrc
    phi = p['phi']
    gam = p.get('gamma', 0.0)
    if vbs <= 0.0:
        sarg = math.sqrt(phi - vbs)
    else:
        sarg = math.sqrt(phi) / (1.0 + vbs / (2.0 * phi))
    vth = p['vto'] + gam * (sarg - math.sqrt(phi))
    vgst = (vg - vsrc) - vth
    vds = vdrn - vsrc
    lam = p.get('lambd', 0.0)
    if vgst <= 0.0:
        ich = 0.0
    elif vds >= vgst:
        ich = 0.5 * beta * vgst ** 2 * (1.0 + lam * vds)
    else:
        ich = beta * (vgst * vds - 0.5 * vds * vds) * (1.0 + lam * vds)
    ## The two bulk junctions, in the n-channel convention.
    isb = p.get('js', 0.0) * p.get('asrc', 0.0) or p.get('IS', 1e-14)
    isd = p.get('js', 0.0) * p.get('ad', 0.0) or p.get('IS', 1e-14)
    ibs = isb * (math.exp(min((vb - vs) / _UT, 300.0)) - 1.0)
    ibd = isd * (math.exp(min((vb - vd) / _UT, 300.0)) - 1.0)
    return sgn * ich, ibs, ibd


def _m1_terminal_currents(vd, vg, vs, vb, card):
    """The four terminal currents the ELEMENT should report, built from
    `_m1_ref` and the model's declared branch orientations."""
    ich, ibs, ibd = _m1_ref(vd, vg, vs, vb, card)
    return np.array([ich - ibd, 0.0, -ich - ibs, ibs + ibd])


## ======================================================================
## 1.  MOS level 1
## ======================================================================

#: Nine bias points spanning cutoff, the linear region, saturation, both
#: signs of Vds and both signs of Vbs.  A single point cannot distinguish
#: a wrong region boundary from a wrong coefficient.
M1_BIAS = [(3.0, 2.0, 0.0, 0.0), (0.3, 2.0, 0.0, 0.0),
           (1.2, 2.0, 0.0, 0.0), (3.0, 0.5, 0.0, 0.0),
           (3.0, 2.0, 0.0, -2.0), (0.3, 2.0, 0.0, -2.0),
           (0.0, 2.0, 3.0, -2.0), (0.0, 2.0, 0.3, 0.0),
           (5.0, 5.0, 0.0, -1.0),
           ## FORWARD BULK, and these two are here because the last set
           ## did not have them: SPICE's threshold takes a different arm
           ## for `vbs > 0` -- a first-order continuation of the square
           ## root instead of the root -- and with every point at
           ## `vbs <= 0` that arm is never selected.  Measured: removing
           ## the factor of 2 from its denominator changed nothing until
           ## these two were added.
           (0.5, 2.0, 0.0, 0.55), (2.0, 2.5, 0.0, 0.65)]


def test_level1_matches_the_textbook_transcription():
    """The four terminal currents against an independent transcription of
    Shichman-Hodges, over nine bias points and both channel types.

    One tolerance for the whole set, not one per point: a model that is
    exact in saturation and wrong in the linear region is accommodated by
    a per-point band and caught by a shared one.
    """
    n = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)
    p = _mk(eh.MosLevel1PmosHdl, 'd', 'g', 's', 'b', **NMOS)
    for x in M1_BIAS:
        want = _m1_terminal_currents(*x, NMOS)
        assert_allclose(n.i(np.array(x)), want, rtol=2e-12, atol=1e-18,
                        err_msg='nmos at %r' % (x,))
        ## The p-channel is the same device with every voltage negated,
        ## so every current is too.  That is a claim about the model's
        ## polarity handling and it is checked against the SAME external
        ## reference, not against the n-channel element.
        assert_allclose(p.i(-np.array(x)), -want, rtol=2e-12, atol=1e-18,
                        err_msg='pmos at %r' % (x,))
    ## The set really does span the regions: the drain current varies by
    ## more than ten decades across it.
    ids = [abs(_m1_ref(*x, NMOS)[0]) for x in M1_BIAS]
    assert max(ids) / max(min(ids), 1e-30) > 1e10


def test_level1_is_symmetric_in_source_and_drain():
    """Exchanging the two diffusions negates the current EXACTLY.

    Not to a tolerance: the model's whole reason for being written as
    ``max(vgs - vth,0)^2 - max(vgd - vth,0)^2`` with one threshold is that
    exchanging the terminals exchanges the two squares, so the property
    is structural and holds to the last bit.

    It is the CHANNEL that is antisymmetric, so the card here switches
    the junction currents off (``IS = js = 0``, which the model floors at
    1e-30 A, twenty-seven decades below the currents compared).  With the
    default 1e-14 A junction the two diffusions sit at different bulk
    biases after the exchange and their leakages differ by 2e-14 A --
    real, correct, and not what this test is about.  Measured rather than
    assumed: the assertion below fails by exactly that on the full card.
    """
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b',
             **dict(NMOS, IS=0.0, js=0.0))
    ## Every point conducts; the property is about the CHANNEL, and a
    ## cutoff point would compare two copies of the leakage instead.
    for (vd, vg, vs, vb) in [(3.0, 2.0, 0.0, -1.0), (0.4, 1.5, 0.0, -1.0),
                             (2.0, 3.5, 1.0, -3.0), (0.0, 3.0, 2.5, 0.0)]:
        fwd = el.i(np.array([vd, vg, vs, vb]))
        rev = el.i(np.array([vs, vg, vd, vb]))
        ## The device is symmetric, so exchanging the two diffusion
        ## POTENTIALS negates the current at each pin -- and the bulk
        ## current, which sees the same pair of junctions, is unchanged.
        assert fwd[0] == -rev[0] and fwd[2] == -rev[2], (fwd, rev)
        assert fwd[3] == rev[3]
        ## and the currents being compared are not zero.
        assert abs(fwd[0]) > 1e-6


def _sign_form_gds(card, vgs, vbs):
    """``d/dvds`` of the ``sign(vds)*f(|vds|)`` formulation at
    ``vds = 0``, as the compiler would derive it.

    This is the CONTROL for the test below: the model deliberately does
    not use this form, and the reason is what this measures.
    """
    from pycircuit.circuit.hdl import sign as hsign, Abs as hAbs, maxc
    import sympy
    from pycircuit.utilities.param import Parameter

    class _Sgn(Behavioural):
        instparams = [Parameter(name='beta', desc='beta', unit='A/V^2',
                                default=1.0),
                      Parameter(name='vth', desc='threshold', unit='V',
                                default=0.5)]

        @staticmethod
        def analog(d, g, s):
            from pycircuit.circuit.hdl import Branch, Contribution, var
            bgs, bds = Branch(g, s), Branch(d, s)
            vgt = var(maxc(bgs.V - vth, 0.0), 'vgt')               # noqa
            va = var(hAbs(bds.V), 'va')
            lin = var(beta * (vgt * va - 0.5 * va * va), 'lin')    # noqa
            sat = var(0.5 * beta * vgt * vgt, 'sat')               # noqa
            mag = var(sympy.Piecewise((sat, va > vgt), (lin, True)), 'mag')
            return Contribution(bds.I, hsign(bds.V) * mag)
    el = _mk(_Sgn, 'd', 'g', 's', beta=card[0], vth=card[1])
    return float(el.G(np.array([0.0, vgs, 0.0]))[0, 0])


def test_level1_gds_at_vds_zero_is_the_on_conductance():
    """A fully-on switch at ``Vds = 0`` must have ``gds = beta*(Vgs - Vth)``,
    and the ``sign``-based formulation returns ZERO there.

    This is the reason `_mos1_analog` is written the way it is, and it is
    asserted against the alternative rather than described: the control
    below is the same physics written as ``sign(vds)*f(|vds|)``, which is
    the obvious way to express SPICE's source/drain swap in one
    expression.  `hdl.sign` has ``fdiff = 0`` and ``sign(0) = 0``, so its
    analytic ``dId/dvds`` at the origin is exactly zero -- a device Newton
    would see as an open circuit at the single most common bias in a
    digital netlist.
    """
    beta = NMOS['kp'] * NMOS['w'] / (NMOS['l'] - 2 * NMOS['ld'])
    vgs, vsb = 2.0, 1.0
    vth = NMOS['vto'] + NMOS['gamma'] * (math.sqrt(NMOS['phi'] + vsb)
                                         - math.sqrt(NMOS['phi']))
    want = beta * (vgs - vth)

    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)
    got = float(el.G(np.array([0.0, vgs, 0.0, -vsb]))[0, 0])
    ## The junction conductance at 1 V of reverse bias is 1e-30 S, thirty
    ## decades below this, so `gds` is the channel alone.
    assert_allclose(got, want, rtol=1e-12)
    assert want > 1e-4

    ## The control: the same physics, the rejected formulation, and it
    ## returns exactly 0.0.  If this ever stops being 0 the comment in
    ## `_mos1_analog` about `sign` is stale and should be re-read.
    assert _sign_form_gds((beta, vth), vgs, -vsb) == 0.0


def test_level1_body_effect_is_the_square_root_law():
    """The threshold extracted from the model's own current, against
    ``vth = vto + gamma*(sqrt(phi + vsb) - sqrt(phi))``.

    Extracted rather than read back: the current in saturation is
    ``(beta/2)*(vgs - vth)^2``, so ``sqrt(2*Id/beta) = vgs - vth`` gives
    the threshold from an observable.  A body-effect coefficient that was
    right only where it is zero would pass a test at ``vsb = 0`` and fail
    this one -- see `validation-design`, "a bias no sweep visits tests
    nothing".
    """
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **dict(NMOS, lambd=0.0))
    beta = NMOS['kp'] * NMOS['w'] / (NMOS['l'] - 2 * NMOS['ld'])
    for vsb in (0.0, 0.5, 1.0, 2.0, 4.0):
        i = el.i(np.array([4.0, 3.0, 0.0, -vsb]))[0]
        vth = 3.0 - math.sqrt(2.0 * i / beta)
        want = NMOS['vto'] + NMOS['gamma'] * (
            math.sqrt(NMOS['phi'] + vsb) - math.sqrt(NMOS['phi']))
        assert_allclose(vth, want, rtol=1e-9, err_msg='vsb=%g' % vsb)
    ## and the effect is large enough to be worth asserting: 4 V of body
    ## bias moves the threshold by more than 500 mV on this card.
    assert NMOS['gamma'] * (math.sqrt(NMOS['phi'] + 4.0)
                            - math.sqrt(NMOS['phi'])) > 0.5


def test_level1_channel_length_modulation_is_the_measured_output_slope():
    """``lambda`` is asserted as an output conductance, never read back.

    In saturation ``Id = Id0*(1 + lambda*Vds)``, so
    ``go = dId/dVds = lambda*Id/(1 + lambda*Vds)``.  That is a statement
    about the slope of a measured curve; the parameter is what it is
    fitted from.
    """
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)
    lam = NMOS['lambd']
    for vds in (2.0, 3.5, 5.0):
        x = np.array([vds, 2.0, 0.0, 0.0])
        i = el.i(x)[0]
        go = float(el.G(x)[0, 0])
        assert_allclose(go, lam * i / (1.0 + lam * vds), rtol=1e-6)
    ## With lambda = 0 the same measurement gives a flat characteristic,
    ## which is what says the term is doing the work and not something
    ## else -- a zero coefficient does not test a scaling.
    flat = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **dict(NMOS, lambd=0.0))
    assert abs(float(flat.G(np.array([3.5, 2.0, 0.0, 0.0]))[0, 0])) < 1e-16


def test_level1_gamma_and_phi_are_derived_from_nsub_and_tox_when_absent():
    """SPICE's rule, and the first model in the tree to use
    ``$param_given`` and ``$limit`` together.

    A card that gives ``gamma`` gets its ``gamma``; a card that does not
    gets ``sqrt(2*q*eps_si*Nsub)/Cox``, and ``phi`` likewise gets
    ``2*Vt*ln(Nsub/ni)``.  That distinction is what a givenness flag is
    for and what a default VALUE cannot express -- a card deliberately
    saying ``gamma = 0`` wants no body effect, and a card saying nothing
    wants the derived one, and those are two different devices.

    Both are measured from the CURRENT, through the threshold
    extraction, so the test cannot be satisfied by the parameter being
    stored correctly and used nowhere.
    """
    nsub, tox = 1e16, 2e-8
    cox = 3.9 * 8.854187817e-12 / tox
    gam_d = math.sqrt(2.0 * _QE * 11.7 * 8.854187817e-12 * nsub * 1e6) / cox
    vtnom = _KB * _T0 / _QE
    phi_d = 2.0 * vtnom * math.log(nsub * 1e6 / (1.45e10 * 1e6))
    assert 0.2 < gam_d < 0.5 and 0.5 < phi_d < 0.9, (gam_d, phi_d)

    base = dict(vto=0.75, kp=8e-5, lambd=0.0, tox=tox, nsub=nsub,
                w=20e-6, l=2e-6, tnom=_T0)
    beta = base['kp'] * base['w'] / base['l']

    def vth_of(el, vsb):
        i = float(el.i(np.array([4.0, 3.0, 0.0, -vsb]))[0])
        return 3.0 - math.sqrt(2.0 * i / beta)

    derived = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **base)
    for vsb in (0.0, 1.0, 3.0):
        want = 0.75 + gam_d * (math.sqrt(phi_d + vsb) - math.sqrt(phi_d))
        assert_allclose(vth_of(derived, vsb), want, rtol=1e-8)

    ## GIVEN wins, and the two are far enough apart that "given" and
    ## "derived" cannot be confused: gamma = 0 is a real card meaning no
    ## body effect, and it is the case a sentinel default cannot carry.
    given0 = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b',
                 **dict(base, gamma=0.0, phi=0.7))
    for vsb in (0.0, 1.0, 3.0):
        assert_allclose(vth_of(given0, vsb), 0.75, rtol=1e-9)
    ## and the derived device really did move over the same span.
    assert vth_of(derived, 3.0) - vth_of(derived, 0.0) > 0.3


def test_level1_junction_capacitance_is_the_hand_differentiated_charge():
    """``C[b][b]`` against ``dQ/dV`` of SPICE's own depletion charge,
    bottom plus sidewall, differentiated by hand.

    Both junctions are on the bulk node, so the bulk row of ``C`` is the
    sum of the two -- which is why the point below has the two junctions
    at different biases: a test with them equal would pass with the two
    grading coefficients exchanged.
    """
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)

    def cj_of(v, cz, m):
        ## SPICE's depletion capacitance with the `fc` linearisation.
        pbv, fc = NMOS['pb'], NMOS['fc']
        if v < fc * pbv:
            return cz * (1.0 - v / pbv) ** (-m)
        f2 = (1.0 - fc) ** (1.0 + m)
        f3 = 1.0 - fc * (1.0 + m)
        return cz * (f3 + m * v / pbv) / f2

    czb = NMOS['cj'] * NMOS['ad']
    czw = NMOS['cjsw'] * NMOS['pd']
    for (vbs, vbd) in [(-1.0, -3.0), (0.0, -2.0), (0.3, -0.5), (-5.0, 0.35)]:
        x = np.array([-vbd, 2.0, -vbs, 0.0])
        Cm = el.C(x)
        want = (cj_of(vbs, czb, NMOS['mj']) + cj_of(vbs, czw, NMOS['mjsw'])
                + cj_of(vbd, czb, NMOS['mj'])
                + cj_of(vbd, czw, NMOS['mjsw']))
        ## The gate overlap capacitances also land on the bulk row
        ## (`cgbo`), so they are subtracted off by name rather than by
        ## switching them off -- which would be a different device.
        got = float(Cm[3, 3]) - NMOS['cgbo'] * (NMOS['l'] - 2 * NMOS['ld'])
        assert_allclose(got, want, rtol=1e-9,
                        err_msg='vbs=%g vbd=%g' % (vbs, vbd))
    ## The two junctions differ by more than 2x over the sweep, so a
    ## model that used one bias for both would not pass.
    assert cj_of(0.3, czb, NMOS['mj']) / cj_of(-5.0, czb, NMOS['mj']) > 2.0


def test_level1_overlap_capacitances_are_constant_and_go_where_declared():
    """``cgso``/``cgdo`` scale with W, ``cgbo`` with the EFFECTIVE length,
    and none of them moves with bias.

    Read as the difference between a card with them and the same card
    without, so the intrinsic junction capacitance cancels exactly.
    """
    on = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)
    off = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b',
              **dict(NMOS, cgso=0.0, cgdo=0.0, cgbo=0.0))
    leff = NMOS['l'] - 2 * NMOS['ld']
    want_g = (NMOS['cgso'] + NMOS['cgdo']) * NMOS['w'] + NMOS['cgbo'] * leff
    for x in ([3.0, 2.0, 0.0, 0.0], [0.2, 5.0, 0.0, -2.0]):
        d = np.asarray(on.C(np.array(x))) - np.asarray(off.C(np.array(x)))
        assert_allclose(float(d[1, 1]), want_g, rtol=1e-12)
        assert_allclose(float(d[1, 0]), -NMOS['cgdo'] * NMOS['w'],
                        rtol=1e-12)
        assert_allclose(float(d[1, 2]), -NMOS['cgso'] * NMOS['w'],
                        rtol=1e-12)
        assert_allclose(float(d[1, 3]), -NMOS['cgbo'] * leff, rtol=1e-12)
    ## `cgbo` scales with the EFFECTIVE length, not the drawn one -- a
    ## 15% difference on this card, which is why `ld` is non-zero in it.
    assert abs(leff - NMOS['l']) / NMOS['l'] > 0.1


def test_level1_ohmic_resistance_collapses_or_makes_internal_nodes():
    """``rd``/``rs`` absent and ``rsh`` zero: no internal nodes at all.
    Either given: two internal nodes and a measured IR drop.

    `Collapse` is the only valid off-switch for a block whose parameter
    is a divisor -- setting the parameter to zero leaves the ``1/rd`` in
    the compiled expression.  The condition here is a CONJUNCTION
    (``rd <= 0 AND rsh*nrd <= 0``), which is what SPICE's "rd if given,
    else rsh times the squares" needs and which nothing in the tree had
    asked a `Collapse` for before.
    """
    bare = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)
    assert [e[2] for e in x_layout(bare)] == ['terminal'] * 4

    ## `rsh` alone is enough: the branch must NOT collapse.
    sheet = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS_R)
    kinds = [e[2] for e in x_layout(sheet)]
    assert kinds.count('internal node') == 2, x_layout(sheet)

    ## and the resistance is the one SPICE computes.
    rd = NMOS_R['rsh'] * NMOS_R['nrd']
    x = np.array([3.0, 2.0, 0.0, 0.0, 3.0, 0.0])
    i = sheet.i(x)
    assert_allclose(float(i[0]), 0.0, atol=1e-30)
    x[4] = 3.0 - 0.01
    assert_allclose(float(sheet.i(x)[0]), 0.01 / rd, rtol=1e-12)

    ## An explicit `rd` overrides the sheet value.
    expl_ = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **dict(NMOS_R, rd=7.0))
    x2 = np.array([3.0, 2.0, 0.0, 0.0, 3.0 - 0.01, 0.0])
    assert_allclose(float(expl_.i(x2)[0]), 0.01 / 7.0, rtol=1e-12)


def test_level1_channel_noise_hits_the_nyquist_and_two_thirds_limits():
    """The Klaassen-Prins integral for a square-law channel, checked at
    its two ends against quantities the model computes independently.

    * at ``Vds = 0`` the channel is a resistor and thermal noise must be
      ``4*k*T*gds`` EXACTLY (Nyquist), with ``gds`` taken from the
      model's own ``G``;
    * in saturation it must be ``(2/3)*4*k*T*gm``, with ``gm`` likewise.

    Neither number is written into the model; both are limits of the
    integral, so agreeing with them is evidence the integral is right.
    """
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b',
             **dict(NMOS, kf=0.0, lambd=0.0))

    def sid(x, f=1e3):
        ## `CY` takes an ANGULAR frequency; the white terms do not care
        ## and the flicker test below does.
        return float(np.real(np.asarray(
            el.CY(np.array(x), 2 * np.pi * f))[0, 0]))

    ## Nyquist, at Vds = 0.
    x0 = [1e-12, 2.5, 0.0, 0.0]
    gds = float(el.G(np.array(x0))[0, 0])
    assert_allclose(sid(x0), 4.0 * _KB * _T0 * gds, rtol=1e-6)
    assert gds > 1e-5

    ## Two thirds, deep in saturation.
    xs = [5.0, 2.5, 0.0, 0.0]
    gm = float(el.G(np.array(xs))[0, 1])
    assert_allclose(sid(xs), 2.0 / 3.0 * 4.0 * _KB * _T0 * gm, rtol=1e-9)
    assert gm > 1e-5
    ## The two limits differ by a third, so passing both is not passing
    ## one formula twice.
    assert abs(sid(x0) / sid(xs) - 1.0) > 0.2


def test_level1_flicker_noise_is_one_over_f_and_scales_with_the_card():
    """``kf*Id^af/(Cox*Leff^2*f)``, checked as a SLOPE in frequency and a
    scaling in ``kf``, with the white term subtracted off."""
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)
    white = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **dict(NMOS, kf=0.0))
    x = np.array([3.0, 2.5, 0.0, 0.0])
    idv = abs(float(el.i(x)[0]))
    cox = 3.9 * 8.854187817e-12 / NMOS['tox']
    leff = NMOS['l'] - 2 * NMOS['ld']
    for f in (1.0, 10.0, 1e3):
        w = 2 * np.pi * f
        got = (float(np.real(np.asarray(el.CY(x, w))[0, 0]))
               - float(np.real(np.asarray(white.CY(x, w))[0, 0])))
        want = NMOS['kf'] * idv ** NMOS['af'] / (cox * leff ** 2 * f)
        assert_allclose(got, want, rtol=1e-9, err_msg='f=%g' % f)
    ## and it is not swamped: at 1 Hz the flicker term is above the
    ## white one, which is what makes the subtraction meaningful.
    f1 = float(np.real(np.asarray(el.CY(x, 2 * np.pi))[0, 0]))
    f1w = float(np.real(np.asarray(white.CY(x, 2 * np.pi))[0, 0]))
    assert f1 > 2.0 * f1w


## ----------------------------------------------------------------------
## 1.9  limit_together -- the reason this model is in the batch

def _m1_variant(mode, nmos=+1):
    """`MosLevel1Hdl` compiled with a different limiter declaration.

    ``'group'`` is what ships; ``'probe'`` is SPICE's four probes declared
    per-probe, which does not compile; ``'none'`` is the unlimited
    device, which is what the limiter has to be measured against.
    """
    from pycircuit.utilities.param import Parameter          # noqa: F401
    return type('_M1_' + mode, (Behavioural,),
                dict(instparams=eh._mos1_params(),
                     analog=staticmethod(eh._mos1_analog(TEMP, nmos, mode))))


def test_the_shipped_mosfet_declares_spices_four_probes_as_one_group():
    """The production claim, read off the compiled element.

    ``fetlim`` on ``(g,si)``, ``limvds`` on ``(di,si)`` and ``pnjlim`` on
    each of ``(b,si)`` and ``(b,di)`` -- four probes, one group.  With the
    default card the two ohmic branches collapse, so the internal nodes
    ARE the terminals and the probes land on pins 1,2 / 0,2 / 3,2 / 3,0.
    """
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)
    info = el._hdl_info
    assert [(s[0], s[1]) for s in info['limit_spec']] == \
        [((1, 2), 'fet'), ((0, 2), 'vds'), ((3, 2), 'pnj'), ((3, 0), 'pnj')]
    assert info['limit_groups'] == [(False, [0, 1, 2, 3])]
    ## and `explain()` says so, which is where a designer would look.
    text = explain(el)
    assert 'fetlim on (g,s) [group 1]' in text, text
    assert text.count('[group 1]') == 4, text

    ## The p-channel declares the same four with every branch reversed,
    ## because `pnjlim` must see the voltage that has the exponential in
    ## it (`GummelPoonPnpHdl` records the same rule).
    p = _mk(eh.MosLevel1PmosHdl, 'd', 'g', 's', 'b', **NMOS)
    assert [(s[0], s[1]) for s in p._hdl_info['limit_spec']] == \
        [((2, 1), 'fet'), ((2, 0), 'vds'), ((2, 3), 'pnj'), ((0, 3), 'pnj')]


def test_the_four_probes_do_not_compile_per_probe_on_this_model():
    """The capability, asserted on the PRODUCTION model rather than on a
    toy: `test_device_limiter.py` shows the refusal on a purpose-built
    four-probe element, and this shows that the refusal is what a real
    SPICE MOSFET runs into.

    Per-probe, each probe's write-back claims a terminal.  ``fetlim``
    takes one of ``g``/``s``, ``limvds`` one of ``d``/``s``,
    ``pnjlim(vbs)`` one of ``b``/``s`` -- and ``pnjlim(vbd)`` finds both
    ``b`` and ``d`` already claimed.  There is no node it can move that
    does not undo an earlier probe, so the model is refused rather than
    silently written wrong.
    """
    with pytest.raises(ValueError, match='over-determine this device'):
        _m1_variant('probe')
    with pytest.raises(ValueError, match='limit_together'):
        _m1_variant('probe')
    ## The unlimited and the grouped forms both compile, which is what
    ## makes the refusal a statement about the DECLARATION and not about
    ## the model.
    assert _m1_variant('none') is not None
    assert _m1_variant('group') is not None


def test_level1_limiting_is_a_no_op_near_a_solution():
    """A step of a few millivolts must pass through EXACTLY.

    That property is what lets "did limiting fire?" be a convergence
    signal, and it is destroyed by a write-back that writes even when it
    did not bite (``a - (a - b)`` is not ``b``).  Asserted on the
    identity of the vector, not on a tolerance.
    """
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)
    rng = np.random.default_rng(20260826)
    x0 = np.array([1.6, 2.0, 0.0, -1.0])
    for _ in range(200):
        x = x0 + rng.uniform(-2e-3, 2e-3, 4)
        assert (el.limit(x, x0) == x).all()


def _forced_bulk(cls, ibulk, vd=0.0):
    """A current forced into the BULK, so both bulk junctions have to
    conduct it.  This is the circuit a MOSFET's ``pnjlim`` exists for, and
    the only part of a level-1 device that is exponential at all -- the
    channel is a polynomial, which is why ``fetlim`` earns much less here
    than it does on an all-region model."""
    c = SubCircuit()
    c['ib'] = ISRC(gnd, 'bulk', i=ibulk)
    c['vd'] = VS('d', gnd, v=vd)
    c['vg'] = VS('g', gnd, v=0.0)
    c['m'] = cls('d', 'g', gnd, 'bulk', **dict(NMOS, cj=0.0, cjsw=0.0))
    return c


def _plain_newton(c, maxiter=100, reltol=1e-4, vabstol=1e-6, iabstol=1e-12):
    """`StandardNewton` on the circuit, WITHOUT `DC`'s rescue chain and
    WITHOUT a gmin anchor.

    Both would rescue this circuit either way, which is exactly why
    neither is used: `validation-design` -- an intervention has to be
    measured with whatever else would paper over it turned off.
    """
    from pycircuit.circuit.nrsolver import StandardNewton
    iref = c.get_node_index(gnd)
    epar = defaultepar

    def rfunc(xr):
        x = np.insert(xr, iref, 0.0)
        F = c.i(x, epar) + c.u(0, analysis='dc', epar=epar)
        J = c.G(x, epar)
        return (np.delete(F, iref),
                np.delete(np.delete(J, iref, 0), iref, 1))

    def limiter(xr, x0r):
        x = c.limit(np.insert(xr, iref, 0.0), np.insert(x0r, iref, 0.0),
                    epar)
        return np.delete(x, iref)

    nn, nb = len(c.nodes), len(c.branches)
    abstol = np.delete(np.concatenate((iabstol * np.ones(nn),
                                       vabstol * np.ones(nb))), iref)
    xtol = np.delete(np.concatenate((vabstol * np.ones(nn),
                                     iabstol * np.ones(nb))), iref)
    x, its = StandardNewton().solve_system(
        np.zeros(c.n - 1), rfunc, numeric, reltol, abstol, xtol, maxiter,
        limiter=limiter)
    return np.insert(x, iref, 0.0), its


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_level1_limiting_rescues_a_solve_and_lands_on_the_same_answer():
    """The only test here that shows the feature doing its job.

    Same circuit, same solver, same starting point, no gmin anchor and no
    rescue ladder; the ONLY difference is whether the element declares
    the four probes.  Over four decades of forced bulk current the
    unlimited element does not converge at all and the grouped one
    converges in **at most three** Jacobian evaluations -- and lands, to
    1e-13 relative, on the point the full `DC` rescue chain finds for the
    unlimited element.  A limiter moves the PATH, never the answer.
    """
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix

    unlim, grp = _m1_variant('none'), _m1_variant('group')
    for i0 in (1e-6, 1e-4, 1e-2, 1.0):
        with pytest.raises((NoConvergenceError, SingularMatrix)):
            _plain_newton(_forced_bulk(unlim, i0))

        c = _forced_bulk(grp, i0)
        x, its = _plain_newton(c)
        assert its <= 3, (i0, its)
        vb = x[c.get_node_index('bulk')]
        ## A real forward-biased junction, not a rail or a zero.
        assert 0.4 < vb < 0.9, (i0, vb)
        ## and the same answer the rescue chain finds with no limiter.
        ref = DC(_forced_bulk(unlim, i0), toolkit=numeric).solve()
        assert_allclose(vb, float(ref.v('bulk')), rtol=1e-13)


def test_the_gmin_anchor_is_not_what_rescues_this_solve():
    """⚠ Any "with X versus without X" convergence claim has to be made
    with the anchor off, or it is measuring the anchor.

    Here the anchor is off in the measurement above (`_plain_newton` runs
    `StandardNewton` directly, which has none), and this asserts the
    other half: the anchor does not change the ANSWER either.  Four
    decades of gmin move the bulk node by nothing at all, because the
    circuit is not singular anywhere -- two forward-biased junctions are
    the best-conditioned rows in the matrix.
    """
    vals = [float(DC(_forced_bulk(_m1_variant('group'), 1e-2),
                     toolkit=numeric, gmin=g).solve().v('bulk'))
            for g in (0.0, 1e-13, 1e-12, 1e-11)]
    assert vals[1:] == vals[:-1], vals


## ----------------------------------------------------------------------
## 1.10  Jacobians, finiteness, temperature and a real circuit

@pytest.mark.parametrize('name,card,x,resolved', [
    ('saturation', NMOS, [3.0, 2.0, 0.0, -1.0], True),
    ('triode', NMOS, [0.2, 2.0, 0.0, -1.0], True),
    ## Every point where the CHANNEL is off or symmetric leaves the two
    ## reverse-biased junctions as the only thing in the drain and bulk
    ## rows, and their conductance there is ~1e-46 S against a MEASURED
    ## roundoff floor of 4.4e-23 -- twenty-three decades under.  Those
    ## entries come back UNRESOLVED, which is exactly right: a central
    ## difference at that scale carries no information about the
    ## derivative.  The verdict is `ok` (the resolvable entries all pass)
    ## and `resolved` is False, and both halves are asserted.
    ('vds exactly zero', NMOS, [1.0, 2.0, 1.0, -1.0], False),
    ('cutoff', NMOS, [3.0, 0.2, 0.0, -1.0], False),
    ('reversed', NMOS, [0.0, 2.0, 3.0, -1.0], True),
    ('forward bulk', NMOS, [0.5, 2.0, 0.0, 0.55], True),
    ('at the fc seam', NMOS, [-0.425, 2.0, 0.0, 0.0], True),
    ('right at threshold', NMOS, [2.0, 0.75, 0.0, 0.0], False),
    ('ideal card', NMOS_IDEAL, [2.0, 2.0, 0.0, 0.0], True),
    ('with parasitics', NMOS_R, [3.0, 2.0, 0.0, -1.0, 2.9, 0.05], True),
])
def test_level1_jacobians_by_finite_differences(name, card, x, resolved):
    """``G`` against ``i`` and ``C`` against ``q``, by central differences,
    at ten points chosen to sit on every seam the model has.

    ``resolved`` is asserted as well as ``ok``, in BOTH directions: seven
    of the ten points resolve completely, and the three that do not are
    named above with the mechanism and the measured floor.  Asserting
    only ``ok`` would let the instrument quietly stop resolving anything.
    """
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **card)
    res = check_jacobians(el, np.array(x, dtype=float), rtol=3e-5)
    assert res.ok, '%s\n%s' % (name, res)
    assert res.resolved == resolved, '%s: resolved=%s\n%s' % (
        name, res.resolved, res)
    ## Every unresolved entry must be unresolved for a REASON the
    ## instrument names, and on this model there are exactly two, both
    ## of them properties of the model rather than of the checker:
    ##
    ## * ROUNDOFF, on a reverse-biased junction conductance of ~1e-46 S
    ##   against a measured floor of 4.4e-23 -- the entry does not move at
    ##   any step because it is below the representable step of `i`;
    ## * TRUNCATION, at ``vgs = vth`` exactly.  The channel current is
    ##   ``(beta/2)*max(vgs - vth, 0)^2``, which is C1 but not C2 there,
    ##   so a central difference returns ``beta*h/4`` -- 4.4e-11 S at the
    ##   default step against an analytic 3.9e-13 S -- and Richardson
    ##   correctly reports that the expansion is not converging.  No `h`
    ##   fixes it: the second derivative is discontinuous.
    for e in res.unresolved:
        assert e.reason in ('roundoff', 'truncation'), (name, e)
        if e.reason == 'roundoff':
            ## the definition of the verdict: the entry is below its own
            ## MEASURED noise floor, so the difference says nothing.
            assert abs(e.ana) < e.floor, (name, e)


def test_level1_jacobian_check_fails_on_a_corrupted_derivative():
    """The instrument's own control: `check_jacobians` must FAIL when the
    Jacobian is wrong, or the ten passes above are decoration.

    Corrupted by a factor of THREE, not by a half.  At a corner between a
    slope of ``g`` and a slope of ``0`` a central difference returns
    ``g/2``, so halving the analytic entry lands exactly on it and the
    control silently passes -- `validation-design` records that
    coincidence biting twice.
    """
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **NMOS)
    x = np.array([3.0, 2.0, 0.0, -1.0])
    assert check_jacobians(el, x, rtol=3e-5).resolved

    good = el.G
    try:
        ## An INSTANCE attribute, so `check_jacobians` still finds the
        ## class's `_hdl_info` and only the derivative is wrong.
        el.G = lambda *a, **kw: 3.0 * good(*a, **kw)
        bad = check_jacobians(el, x, rtol=3e-5)
    finally:
        del el.G
    assert not bad.ok and bad.verdict == 'FAILED', bad
    ## and the same corruption on `C` is caught too, which is the half a
    ## current-only check would miss.
    goodc = el.C
    try:
        el.C = lambda *a, **kw: 3.0 * goodc(*a, **kw)
        badc = check_jacobians(el, x, rtol=3e-5)
    finally:
        del el.C
    assert not badc.ok, badc


def test_level1_stays_finite_where_no_device_belongs():
    """Value and Jacobian scanned over absurd bias, on BOTH channel types
    and on a card with every block on.

    A model can be finite everywhere and have a NaN Jacobian at one bias,
    which is the worse failure because nothing looks wrong.  Scanned for
    the boundary rather than asserted to break somewhere.
    """
    for cls, sgn in ((eh.MosLevel1Hdl, +1.0), (eh.MosLevel1PmosHdl, -1.0)):
        el = _mk(cls, 'd', 'g', 's', 'b', **NMOS)
        rng = np.random.default_rng(20260826)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for _ in range(300):
                xa = sgn * rng.uniform(-300.0, 300.0, 4)
                for f in (el.i, el.q, el.G, el.C):
                    out = np.asarray(f(xa), dtype=float)
                    assert np.isfinite(out).all(), (cls, xa, f, out)
            ## and the extremes, where every regulariser is past its
            ## seam.  The BOUNDARY is measured rather than asserted to
            ## exist: `expl` continues past its seam with exp's own
            ## third-order Taylor, so `expl(v/VT)` overflows once the
            ## cubic term does -- measured on the bulk node, finite
            ## through 1e68 V and not at 1e69 V.  1e60 is the assertion,
            ## eight decades inside the measured edge, so a regression
            ## that costs a few decades still fails it.
            for v in (1e2, 1e6, 1e12, 1e30, 1e60):
                for x in ([v, 0.0, 0.0, 0.0], [0.0, v, 0.0, 0.0],
                          [0.0, 0.0, v, 0.0], [0.0, 0.0, 0.0, v],
                          [v, v, -v, -v], [-v, -v, v, v]):
                    xa = sgn * np.array(x)
                    for f in (el.i, el.q, el.G, el.C):
                        out = np.asarray(f(xa), dtype=float)
                        assert np.isfinite(out).all(), (cls, v, x, f, out)


def test_level1_temperature_path_moves_the_threshold_and_the_mobility():
    """SPICE's ``mos1temp.c``: ``KP`` falls as ``T^-1.5``, ``PHI`` follows
    the band gap and ``VTO`` follows ``PHI``.

    Each is measured from the quantity it moves -- the mobility from the
    saturation current at a FIXED overdrive above the moved threshold,
    the threshold from the classic ``sqrt(Id)`` extrapolation -- rather
    than read back from a parameter.
    """
    ## `IS = 0` (floored at 1e-30 A) so that the extraction sees the
    ## CHANNEL alone.  With SPICE's 1e-14 A default the drain junction's
    ## own temperature path raises its leakage to 9.9e-10 A at 400 K,
    ## which is 5e-6 of the drain current and shifts the extracted
    ## threshold by 2.6 uV -- correct physics, and not what this test is
    ## about.  Measured, not assumed: with the default card the
    ## assertions below fail by exactly that.
    card = dict(NMOS_IDEAL, tnom=300.0, IS=0.0, js=0.0)
    el = _mk(eh.MosLevel1Hdl, 'd', 'g', 's', 'b', **card)

    def epar_at(T):
        e = defaultepar.copy()
        e.T = T
        return e

    def ids(T, vg, vd=5.0):
        return float(np.asarray(el.i(np.array([vd, vg, 0.0, 0.0]),
                                     epar_at(T)), float)[0])

    def vth_at(T):
        v1, v2 = 2.0, 3.0
        s1, s2 = math.sqrt(ids(T, v1)), math.sqrt(ids(T, v2))
        return v1 - s1 * (v2 - v1) / (s2 - s1)

    ## Mobility: at a fixed overdrive the saturation current goes as KP.
    for T in (250.0, 350.0, 400.0):
        r = ids(T, vth_at(T) + 2.0) / ids(300.0, vth_at(300.0) + 2.0)
        assert_allclose(r, (T / 300.0) ** -1.5, rtol=1e-9,
                        err_msg='T=%g' % T)
    ## Threshold: with gamma = 0 the whole shift is `(phi(T) - phi)/2`,
    ## which is the band-gap expression and is computed here from it.
    def phi_at(T):
        eg = 1.16 - 7.02e-4 * T * T / (T + 1108.0)
        egn = 1.16 - 7.02e-4 * 300.0 ** 2 / (300.0 + 1108.0)
        return (card['phi'] * T / 300.0
                - 3.0 * (_KB * T / _QE) * math.log(T / 300.0)
                - egn * T / 300.0 + eg)
    for T in (250.0, 400.0):
        want = card['vto'] + 0.5 * (phi_at(T) - card['phi'])
        assert_allclose(vth_at(T), want, rtol=1e-9, err_msg='T=%g' % T)
    ## and the shift is real: 100 K moves the threshold by >20 mV.
    assert abs(vth_at(400.0) - vth_at(300.0)) > 0.02


def _inverter_vout(vin, vdd, ncard, pcard):
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=vdd)
    c['vin'] = VS('vin', gnd, v=vin)
    c['mn'] = eh.MosLevel1Hdl('out', 'vin', gnd, gnd, **ncard)
    c['mp'] = eh.MosLevel1PmosHdl('out', 'vin', 'vdd', 'vdd', **pcard)
    return float(DC(c, toolkit=numeric).solve().v('out'))


def _switching_point(vdd, ncard, pcard):
    """The input at which the output sits at ``VDD/2``, by bisection on
    real DC solves.

    Defined that way on purpose: with the output pinned at ``VDD/2`` both
    devices see the SAME ``|Vds|``, so their channel-length-modulation
    factors ``(1 + lambda*|Vds|)`` are equal and cancel out of the
    current balance -- which makes the closed form below exact rather
    than an approximation with ``lambda`` neglected.
    """
    lo, hi = 0.0, vdd
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _inverter_vout(mid, vdd, ncard, pcard) > 0.5 * vdd:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def test_a_cmos_inverter_switches_where_the_textbook_says():
    """A real DC solve of two of these, against the closed-form switching
    threshold.

    With both devices in saturation the currents balance, so

        Vm = (VDD - |Vtp| + Vtn*sqrt(bn/bp)) / (1 + sqrt(bn/bp))

    -- Weste and Harris, *CMOS VLSI Design*, section 2.5.  Two ratios are
    checked, because a symmetric inverter switches at ``VDD/2`` for
    reasons that have nothing to do with the formula being right: a
    4x-wider n-channel must move the threshold DOWN by 580 mV and the
    formula must say by how much.

    **``lambda`` must be non-zero here, and that is a property of the
    circuit rather than of the model.**  With ``lambda = 0`` both devices
    are ideal current sources in saturation, so at ``Vin = Vm`` every
    output voltage in the saturation window satisfies KCL exactly and the
    operating point is not unique -- the DC solve lands wherever Newton
    stops (measured: 3.2502 V on a 5 V supply, a perfectly good solution
    of an underdetermined problem).  A finite output conductance is what
    makes the switching point a point.
    """
    vdd = 5.0
    card = dict(NMOS_IDEAL, lambd=0.03)
    wide = dict(card, w=4.0 * card['w'])

    ## Symmetric: exact by the device symmetry, and asserted tightly.
    assert_allclose(_switching_point(vdd, card, card), 0.5 * vdd, rtol=1e-9)

    ## and the curve really swings rail to rail, or "it crosses at
    ## VDD/2" would be satisfied by a straight line.
    assert _inverter_vout(4.5, vdd, card, card) < 1e-6
    assert _inverter_vout(0.5, vdd, card, card) > vdd - 1e-6

    ## Asymmetric: the formula, with nothing fitted.
    ratio = 2.0                                   # sqrt(bn/bp) = sqrt(4)
    want = (vdd - card['vto'] + card['vto'] * ratio) / (1.0 + ratio)
    got = _switching_point(vdd, wide, card)
    assert_allclose(got, want, rtol=1e-9)
    ## The two cases differ by more than half a volt, so passing both is
    ## not passing one number twice.
    assert abs(want - 0.5 * vdd) > 0.5


## ======================================================================
## 2.  The Boyle-class opamp macromodel
## ======================================================================

#: A 741-class card with the offsets zeroed, so that the gain and
#: bandwidth relations can be measured without an operating-point offset
#: saturating the open-loop amplifier.  The offsets get their own card.
OPA = dict(aol=2e5, gbw=1e6, sr=5e5, cc=30e-12, vos=0.0, ib=0.0, ios=0.0,
           rin=2e6, ricm=1e9, cmrr=1e5, psrr=1e5, vsupnom=30.0,
           rout=75.0, isc=25e-3, vdrop=1.5)

#: The same amplifier with real input imperfections.
OPA_OFF = dict(OPA, vos=1.2e-3, ib=80e-9, ios=20e-9)

VCC, VEE = 15.0, -15.0


def _opamp_tb(vp=0.0, vn=0.0, card=None, vac_p=0.0, vac_n=0.0,
              vac_cc=0.0, vac_ee=0.0, out='o'):
    """Open-loop testbench.  Every source's ``vac`` is given EXPLICITLY:
    `VS`'s default is 1.0, so a source left alone drives the AC analysis
    -- which cost an hour here, because every node in the circuit came
    back at exactly 1.0 and looked like a singular matrix."""
    c = SubCircuit()
    c['vcc'] = VS('vcc', gnd, v=VCC, vac=vac_cc)
    c['vee'] = VS('vee', gnd, v=VEE, vac=vac_ee)
    c['vp'] = VS('vp', gnd, v=vp, vac=vac_p)
    c['vn'] = VS('vn', gnd, v=vn, vac=vac_n)
    c['x'] = eh.OpAmpHdl('vp', 'vn', out, 'vcc', 'vee',
                         **(card if card is not None else OPA))
    return c


def _mid_referred(res, node):
    """The output referred to the SUPPLY MIDPOINT, which is this model's
    reference -- it has no ground pin (see `OpAmpHdl`)."""
    return (np.asarray(res.v(node).y)
            - 0.5 * (np.asarray(res.v('vcc').y)
                     + np.asarray(res.v('vee').y)))


def test_opamp_open_loop_gain_and_dominant_pole_are_boyles_relations():
    """``A_ol``, ``f_dominant = GBW/A_ol`` and ``GBW`` from one AC sweep.

    All three are IDENTITIES of Boyle's parameterisation -- ``gm1`` is
    set to ``2*pi*gbw*cc`` and ``r2`` to ``aol/gm1``, so the pole is
    ``1/(2*pi*r2*cc) = gbw/aol`` by construction -- which is why they are
    asserted to parts in 1e4 off a real solve rather than to a band.
    Reference: Boyle et al., *IEEE JSSC* **SC-9**, 353 (1974).
    """
    c = _opamp_tb(vac_p=1.0)
    fp = OPA['gbw'] / OPA['aol']
    f = np.array([1e-3, fp, 1e3, 1e5, OPA['gbw']])
    res = AC(c).solve(freqs=f)
    a = np.abs(_mid_referred(res, 'o'))
    ph = np.angle(np.asarray(res.v('o').y), deg=True)

    ## DC gain.  The common-mode term is real and is subtracted rather
    ## than ignored: driving one input alone puts half a volt of common
    ## mode in, worth `aol/(2*cmrr)` = 1.0 of output.
    assert_allclose(a[0], OPA['aol'] * (1.0 + 0.5 / OPA['cmrr']), rtol=1e-6)
    ## The dominant pole: -3 dB and -45 degrees at gbw/aol, both.
    assert_allclose(a[1] / a[0], 1.0 / math.sqrt(2.0), rtol=1e-6)
    assert_allclose(ph[1], -45.0, rtol=1e-4)
    ## Gain-bandwidth: above the pole the product is `gbw` corrected by
    ## the single-pole shape, `gbw/sqrt(1 + (fp/f)^2)`.  The correction is
    ## 1.25e-5 at 1 kHz and 1.25e-9 at 100 kHz, so the EXACT form is
    ## asserted rather than a band wide enough to hold the raw product --
    ## a band that loose would also hold a two-pole amplifier.
    for k in (2, 3, 4):
        want = OPA['gbw'] / math.sqrt(1.0 + (fp / f[k]) ** 2)
        assert_allclose(f[k] * a[k] / (1.0 + 0.5 / OPA['cmrr']),
                        want, rtol=1e-7)
    ## and unity gain lands at gbw, which is what the name means.
    assert_allclose(a[4], 1.0, rtol=1e-4)


def test_opamp_slew_rate_is_the_tail_current_over_the_compensation_cap():
    """Boyle's ``SR = I/C``, measured as the maximum slope of a real
    transient on a unity-gain buffer.

    The drive must DEMAND more slope than the amplifier can give, or the
    measurement is of the source: 8 V at 100 kHz asks for 5.0e6 V/s
    against the card's 5.0e5, ten times over.  Measured rather than
    assumed -- at 1 kHz the same buffer follows the input exactly and the
    maximum slope is the sine's own 6.3e4 V/s.
    """
    def maxslope(freq, amp):
        c = SubCircuit()
        c['vcc'] = VS('vcc', gnd, v=VCC, vac=0.0)
        c['vee'] = VS('vee', gnd, v=VEE, vac=0.0)
        c['vin'] = eh.VSinHdl('vp', gnd, vo=0.0, va=amp, freq=freq)
        c['x'] = eh.OpAmpHdl('vp', 'out', 'out', 'vcc', 'vee', **OPA)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = Transient(c, toolkit=numeric).solve(
                tend=2.0 / freq, timestep=1e-8 / (freq / 1e5),
                fixed_timestep=True)
        t = np.asarray(res.v('out').x[0], float)
        v = np.asarray(res.v('out').y, float)
        s = np.diff(v) / np.diff(t)
        return s.max(), s.min()

    up, dn = maxslope(1e5, 8.0)
    assert_allclose(up, OPA['sr'], rtol=2e-3)
    assert_allclose(-dn, OPA['sr'], rtol=2e-3)
    ## and the relation is `sr*cc`, i.e. halving `cc` at fixed `sr` is
    ## the same amplifier: the slew limit is a CURRENT limit and the
    ## parameterisation is what makes it come out as `sr`.
    slow, _ = maxslope(1e3, 8.0)
    assert_allclose(slow, 2 * np.pi * 1e3 * 8.0, rtol=5e-3)
    assert slow < 0.2 * OPA['sr']


def test_opamp_output_clamps_a_vdrop_short_of_each_rail():
    """Open-loop overdrive puts the output at ``vcc - vdrop`` and
    ``vee + vdrop``, and the clamp does NOT cost gain in the linear
    region.

    The second half is the one worth asserting: the obvious smooth clamp
    ``vs*tanh(v/vs)`` is 8% low at half the swing, which would show up as
    a gain error.  `hypsmooth` with a width relative to the swing is
    2e-6 V low at 1 V of output, so the DC gain measured at 1 V of output
    is the card's ``aol`` to seven digits.
    """
    hi = float(DC(_opamp_tb(vp=1.0), toolkit=numeric).solve().v('o'))
    lo = float(DC(_opamp_tb(vp=-1.0), toolkit=numeric).solve().v('o'))
    assert_allclose(hi, VCC - OPA['vdrop'], rtol=1e-9)
    assert_allclose(lo, VEE + OPA['vdrop'], rtol=1e-9)
    ## Linear region: 1 V out of a 13.5 V swing, gain unaffected.
    vin = 1.0 / OPA['aol']
    got = float(DC(_opamp_tb(vp=vin), toolkit=numeric).solve().v('o'))
    want = OPA['aol'] * vin * (1.0 + 0.5 / OPA['cmrr'])
    assert_allclose(got, want, rtol=2e-6)


def test_opamp_output_resistance_and_current_limit_are_one_expression():
    """``isc*tanh(dv/(rout*isc))`` has slope ``1/rout`` at the origin and
    asymptote ``isc``, so one contribution carries both numbers.

    Both are measured from DC solves with real loads: a light load reads
    the output resistance off the voltage drop, and a heavy one reads the
    current limit off the load current.
    """
    def loaded(rl, vp=0.0):
        c = _opamp_tb(vp=vp)
        c['rl'] = R('o', gnd, r=rl)
        c.update_iparv()
        return float(DC(c, toolkit=numeric).solve().v('o'))

    ## Output resistance, measured where `tanh` is linear: 1 V of output
    ## into 100 kOhm is 10 uA, four thousand times below `isc`, so the
    ## expression's slope at the origin is what is being read.  The open
    ## -circuit output is the reference, so what is measured is the DROP.
    vin = 1.0 / OPA['aol']
    vopen = loaded(1e12, vp=vin)
    for rl in (1e5, 1e4):
        v = loaded(rl, vp=vin)
        i = v / rl
        assert_allclose((vopen - v) / i, OPA['rout'], rtol=1e-3,
                        err_msg='rl=%g' % rl)
        assert i < 5e-3 * OPA['isc']
    ## Drive it to the rail and load it hard: the current saturates at
    ## `isc` and the output voltage follows the load.
    for rl in (100.0, 10.0):
        v = loaded(rl, vp=1.0)
        assert_allclose(v / rl, OPA['isc'], rtol=1e-4,
                        err_msg='rl=%g' % rl)
    ## and the limit really binds: unlimited, a 10 Ohm load at 13.5 V
    ## would be 1.35 A, fifty times the limit.
    assert (VCC - OPA['vdrop']) / 10.0 > 50 * OPA['isc']


def test_opamp_cmrr_and_psrr_are_ratios_of_measured_gains():
    """CMRR and PSRR asserted as ``Ad/Acm`` and ``Ad/Asupply``, from three
    AC solves of the same circuit with different sources driven.

    **PSRR needs the supply driven DIFFERENTIALLY** (``vcc`` up by half a
    volt, ``vee`` down by half), and this is a consequence of the model
    having no ground pin: the output is referred to the supply midpoint,
    so moving one rail alone moves the midpoint, the output follows it,
    and a naive measurement to ground reads 1.5 V of "supply gain" that
    is really the reference moving.  Measured: 1.5 against the 2.0 the
    differential drive gives.  Recorded here rather than hidden, because
    it is what a user connecting this to a single supply will meet.
    """
    def gain(**kw):
        res = AC(_opamp_tb(**kw)).solve(freqs=np.array([1e-3]))
        return abs(_mid_referred(res, 'o')[0])

    ad = gain(vac_p=1.0, vac_n=-1.0) / 2.0     # per volt of differential
    acm = gain(vac_p=1.0, vac_n=1.0)           # per volt of common mode
    assert_allclose(ad / acm, OPA['cmrr'], rtol=1e-6)

    ## PSRR: a differential supply drive of 1 V total, midpoint fixed.
    asup = gain(vac_cc=0.5, vac_ee=-0.5)
    assert_allclose(ad / asup, OPA['psrr'], rtol=1e-6)

    ## The single-rail drive, measured and named: 1.5 V of output to
    ## ground per volt of vcc, of which 0.5 is the midpoint itself.
    res = AC(_opamp_tb(vac_cc=1.0)).solve(freqs=np.array([1e-3]))
    assert_allclose(abs(np.asarray(res.v('o').y)[0]), 1.5, rtol=1e-4)


def test_opamp_input_stage_offsets_and_currents():
    """``vos``, ``ib``, ``ios``, ``rin`` and ``ricm``, each measured as
    the thing it is defined to be.

    ``vos`` is the input voltage that NULLS the output -- not a voltage
    added somewhere -- so it is found by solving for it, and the sign is
    part of the claim.
    """
    ## vos: the output is zero when the DIFFERENTIAL input equals it.
    ## Driving one input alone leaves half of `vos` as common mode, which
    ## the finite CMRR turns into 1.2 mV of output -- correct, and not
    ## what this measures.  (Measured: the single-ended null misses by
    ## exactly `aol*vos/(2*cmrr)` = 1.2e-3 V.)
    half = 0.5 * OPA_OFF['vos']
    v = float(DC(_opamp_tb(vp=half, vn=-half, card=OPA_OFF),
                 toolkit=numeric).solve().v('o'))
    assert abs(v) < 1e-9, v
    ## and it is not zero at zero input -- with this gain it is at a rail.
    v0 = float(DC(_opamp_tb(card=OPA_OFF), toolkit=numeric).solve().v('o'))
    assert_allclose(v0, VEE + OPA_OFF['vdrop'], rtol=1e-9)

    ## ib and ios: the two input currents are `ib +- ios/2`.
    ## `rin` is taken out of the way here, because a 1.2 mV differential
    ## null across 2 MOhm is 600 pA -- 0.67% of `ib`, real, and measured
    ## in its own right two assertions below.
    c = _opamp_tb(vp=half, vn=-half,
                  card=dict(OPA_OFF, rin=1e15, ricm=1e15))
    res = DC(c, toolkit=numeric).solve()
    ip = -float(res.i('vp.plus'))
    inn = -float(res.i('vn.plus'))
    assert_allclose(ip, OPA_OFF['ib'] + 0.5 * OPA_OFF['ios'], rtol=1e-6)
    assert_allclose(inn, OPA_OFF['ib'] - 0.5 * OPA_OFF['ios'], rtol=1e-6)
    assert_allclose(0.5 * (ip + inn), OPA_OFF['ib'], rtol=1e-9)
    assert_allclose(ip - inn, OPA_OFF['ios'], rtol=1e-9)

    ## rin: the differential input resistance, read as the extra current
    ## a differential input draws.
    c2 = _opamp_tb(vp=0.5, vn=-0.5, card=dict(OPA, ricm=1e15))
    r2 = DC(c2, toolkit=numeric).solve()
    assert_allclose(-float(r2.i('vp.plus')), 1.0 / OPA['rin'], rtol=1e-5)
    ## ricm: a common-mode input draws through it and not through rin.
    c3 = _opamp_tb(vp=1.0, vn=1.0, card=dict(OPA, ricm=1e7))
    r3 = DC(c3, toolkit=numeric).solve()
    assert_allclose(-float(r3.i('vp.plus')), 1.0 / 1e7, rtol=1e-5)


def test_opamp_closed_loop_inverting_gain_and_bandwidth():
    """An inverting amplifier, against the exact single-pole closed-loop
    expressions -- which is what a macromodel is FOR.

    ``A_cl = -R2/R1`` corrected by the finite loop gain, and the
    closed-loop pole at ``fp*(1 + A0*beta)``.  The second number is
    computed exactly rather than as the usual ``GBW/(1 + |A_cl|)``: that
    approximation is 0.26% off here, which is ten times the tolerance
    this asserts to.
    """
    r1, r2 = 10e3, 100e3
    c = SubCircuit()
    c['vcc'] = VS('vcc', gnd, v=VCC, vac=0.0)
    c['vee'] = VS('vee', gnd, v=VEE, vac=0.0)
    c['vin'] = VS('vin', gnd, v=0.0, vac=1.0)
    c['r1'] = R('vin', 'inn', r=r1)
    c['r2'] = R('inn', 'o', r=r2)
    ## Two parasitics are taken out of the way, and BOTH sizes are
    ## measured rather than assumed.  `rin = 2 MOhm` loads the summing
    ## node, so the feedback factor is no longer `r1/(r1 + r2)` and the
    ## closed-loop pole moves by **0.26%**; `rout = 75 Ohm` divides with
    ## the 110 kOhm feedback network and moves it by another **0.068%**.
    ## Neither is a defect -- both are what the parameters mean -- but
    ## the closed form below is about the LOOP, and with them in it is no
    ## longer exact.  The shipped card is checked against the same
    ## relation at the end, to the size those two effects predict.
    card = dict(OPA, rin=1e12, ricm=1e15, rout=1e-6)
    c['x'] = eh.OpAmpHdl(gnd, 'inn', 'o', 'vcc', 'vee', **card)

    beta = r1 / (r1 + r2)
    a0 = OPA['aol']
    fp = OPA['gbw'] / a0
    fcl = fp * (1.0 + a0 * beta)
    acl = -(r2 / r1) / (1.0 + 1.0 / (a0 * beta))

    res = AC(c).solve(freqs=np.array([1.0, fcl]))
    g = np.asarray(res.v('o').y)
    assert_allclose(abs(g[0]), abs(acl), rtol=1e-6)
    assert abs(np.angle(g[0], deg=True)) > 179.0
    ## -3 dB at the closed-loop pole, to a part in 1e4.
    assert_allclose(abs(g[1]) / abs(g[0]), 1.0 / math.sqrt(2.0), rtol=1e-4)
    ## and the loop really is fast: the closed-loop bandwidth is four
    ## decades above the open-loop pole.
    assert fcl / fp > 1e4

    ## The SHIPPED card, same circuit, same closed form: within 0.4% of
    ## the ideal relation, which is the sum of the two parasitic effects
    ## named above.  A band that would hold a two-pole amplifier is not
    ## being asserted -- 0.4% at the -3 dB point is a 0.8% error in the
    ## pole frequency.
    c['x'] = eh.OpAmpHdl(gnd, 'inn', 'o', 'vcc', 'vee', **OPA)
    c.update_iparv()
    g2 = np.asarray(AC(c).solve(freqs=np.array([1.0, fcl])).v('o').y)
    assert_allclose(abs(g2[1]) / abs(g2[0]), 1.0 / math.sqrt(2.0),
                    rtol=4e-3)


@pytest.mark.parametrize('name,card,x', [
    ('quiescent', OPA, [0.0, 0.0, 0.0, VCC, VEE, 0.0, 0.0, 0.0]),
    ('linear', OPA, [1e-5, 0.0, 2.0, VCC, VEE, 2.0, 2.0, 0.0]),
    ('saturated', OPA, [1.0, 0.0, 13.5, VCC, VEE, 200.0, 13.5, 0.0]),
    ('slewing', OPA, [0.5, 0.0, 0.0, VCC, VEE, 1.0, 0.0, 0.02]),
    ('loaded', OPA, [0.0, 0.0, -5.0, VCC, VEE, 0.0, 0.0, 0.02]),
    ('with offsets', OPA_OFF, [0.0, 0.0, 0.0, VCC, VEE, 0.0, 0.0, 0.0]),
])
def test_opamp_jacobians_by_finite_differences(name, card, x):
    """``G`` against ``i`` and ``C`` against ``q`` at six operating points
    across the amplifier's whole range, including both nonlinearities.

    The ``x`` vector is ``(inp, inn, out, vcc, vee, gain, oint, i_br)`` --
    `x_layout` says so, and the last entry is the branch current of the
    output stage's voltage contribution.
    """
    el = _mk(eh.OpAmpHdl, 'inp', 'inn', 'out', 'vcc', 'vee', **card)
    assert [e[1] for e in x_layout(el)] == [
        'inp', 'inn', 'out', 'vcc', 'vee', 'gain', 'oint', '_i_br0']
    res = check_jacobians(el, np.array(x, dtype=float), rtol=3e-5)
    assert res.ok, '%s\n%s' % (name, res)
    assert res.resolved, '%s: %d unresolved\n%s' % (
        name, len(res.unresolved), res)


def test_opamp_stays_finite_where_no_amplifier_belongs():
    """Both nonlinearities are saturating -- ``tanh`` on the output
    current and a `hypsmooth` pair on the voltage -- so the value is
    bounded by construction; what has to be checked is the JACOBIAN,
    which divides by things the value does not."""
    el = _mk(eh.OpAmpHdl, 'inp', 'inn', 'out', 'vcc', 'vee', **OPA_OFF)
    rng = np.random.default_rng(20260826)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(300):
            x = rng.uniform(-1e3, 1e3, 8)
            for f in (el.i, el.q, el.G, el.C):
                assert np.isfinite(np.asarray(f(x), float)).all(), x
        ## The BOUNDARY is measured rather than asserted to exist.  On a
        ## supply pin the Jacobian goes first, at 1e107 V, and the value
        ## follows at 1e160 -- both from `hypsmooth`, whose radicand
        ## squares its argument (documented range ~1e154) and whose
        ## smoothing width here is `1e-4` of the SWING, so it squares a
        ## quantity a decade and a half larger still.  The assertion is
        ## at 1e60, forty-seven decades inside the measured edge.
        for v in (1e3, 1e9, 1e30, 1e60):
            for k in range(8):
                x = np.zeros(8)
                x[k] = v
                for f in (el.i, el.q, el.G, el.C):
                    assert np.isfinite(np.asarray(f(x), float)).all(), (v, k)


def test_what_the_hdl_opamp_does_and_does_not_take_over_from_macromodels():
    """The roadmap's item 4 says an HDL opamp "deletes
    ``macromodels.OpAmp``'s runtime-``import jax`` lambdas".  Verified
    rather than repeated, and the answer is: the claim is TRUE and still
    live, it is NOT a drop-in replacement, and the runtime import is the
    smaller of the two things that go away.

    Measured here:

    * ``macromodels.OpAmp`` still carries ``import jax.numpy as jnp``
      INSIDE the closure it hands to `BSource`, so the import runs on
      every evaluation and on every finite-difference probe of it.  The
      *element* half of that pattern was fixed by `toolkit.derivative`;
      the *user callable* half was not, because it is user code;
    * an HDL element has no Python callable in it at all;
    * ``BSource``'s Jacobian is `toolkit.derivative`, a central
      difference at ``eps = 1e-6``.  Accurate here (1.2e-10 relative on
      this ``tanh``), so this is a structural point rather than an
      accuracy one: two extra evaluations per entry, and nothing a
      tracing backend can differentiate;
    * the bigger difference is the CLAMP.  ``macromodels.OpAmp`` limits
      with ``Vspan*tanh(v/Vspan)``, whose gain error is
      ``(v/Vspan)^2/3`` -- **0.15% at 1 V out of a 15 V swing and 12.6%
      at 8.7 V**, which is why its own test asserts the DC gain to 2e-3.
      `OpAmpHdl`'s `hypsmooth` pair is 2e-8 at the same points;
    * it is NOT a drop-in: ``macromodels.OpAmp`` is a four-terminal
      ``SubCircuit`` with a differential output and no supply pins;
      `OpAmpHdl` has five terminals, two of which are the supplies it
      takes its rails and its PSRR from.  Anything that wants a floating
      output reference should keep the old one.
    """
    import inspect
    from pycircuit.circuit import macromodels
    src = inspect.getsource(macromodels)
    assert 'import jax.numpy as jnp' in src
    ## and it is inside a function, not at module scope -- which is what
    ## makes it a per-call cost.
    fn = inspect.getsource(macromodels.OpAmp.__init__)
    assert 'import jax.numpy as jnp' in fn

    ## The HDL element has no callable parameters at all.
    el = _mk(eh.OpAmpHdl, 'inp', 'inn', 'out', 'vcc', 'vee', **OPA)
    assert not any(callable(p.default) for p in eh.OpAmpHdl.instparams)
    assert el._hdl_info['chained']

    ## Not a drop-in.
    assert macromodels.OpAmp.terminals == ('inp', 'inn', 'outp', 'outn')
    assert tuple(eh.OpAmpHdl.terminals) == ('inp', 'inn', 'out',
                                            'vcc', 'vee')

    ## The clamp, measured on both, same gain and same rails.
    aol = 1e5
    def old(vin):
        c = SubCircuit()
        c['V1'] = VS('inp', gnd, v=vin, vac=0.0)
        c['OP1'] = macromodels.OpAmp('inp', gnd, 'outp', gnd, Aol=aol,
                                     GBW=1e6, Vdd=15.0, Vss=-15.0)
        return float(DC(c, toolkit=numeric).solve().v('outp'))

    def new(vin):
        c = SubCircuit()
        c['vcc'] = VS('vcc', gnd, v=15.0, vac=0.0)
        c['vee'] = VS('vee', gnd, v=-15.0, vac=0.0)
        c['V1'] = VS('inp', gnd, v=vin, vac=0.0)
        c['vn'] = VS('inn', gnd, v=0.0, vac=0.0)
        c['x'] = eh.OpAmpHdl('inp', 'inn', 'outp', 'vcc', 'vee',
                             **dict(OPA, aol=aol, sr=1e12, cmrr=1e12,
                                    psrr=1e12, rout=50.0, vdrop=0.0,
                                    rin=1e12, ricm=1e15, isc=1.0))
        return float(DC(c, toolkit=numeric).solve().v('outp'))

    for vin, oldband in ((1e-5, 1.4e-3), (1e-4, 0.12)):
        want = vin * aol
        assert abs(old(vin) - want) / want > oldband, vin
        assert abs(new(vin) - want) / want < 1e-6, vin
    ## and both really do clamp at the rail, so the difference above is
    ## gain error in the LINEAR region and not a different rail.
    assert_allclose(old(1e-3), 15.0, rtol=1e-4)
    assert_allclose(new(1e-3), 15.0, rtol=1e-9)


## ======================================================================
## 3.  The self-heating Gummel-Poon
## ======================================================================

#: An ideal transport n-p-n: no parasitic resistance, no leakage, no
#: knees, no Early effect, no charge.  Everything that is off here is off
#: so that the independent transcription below can be short enough to be
#: read and checked by eye -- the FULL card is exercised by
#: `test_elements_hdl_library3.py`, and the one thing this batch adds is
#: the thermal node.
NPN_TH = dict(IS=1e-16, bf=100.0, nf=1.0, nr=1.0, br=1.0, eg=1.11,
              xti=3.0, tnom=300.0)

TH_VBE, TH_VCE = 0.7, 5.0


def _gp_power(dT, vbe=TH_VBE, vce=TH_VCE, card=NPN_TH):
    """The static dissipation of an ideal transport n-p-n at a junction
    temperature ``T0 + dT``, transcribed from Gummel & Poon as SPICE
    reduces it (Massobrio & Antognetti ch. 2).

    ``qb`` is exactly 1 on this card -- both Early voltages and both knee
    currents are at SPICE's "infinite" -- so the whole model is two
    exponentials, which is what makes this reference checkable by eye.
    """
    T = _T0 + dT
    vt = _KB * T / _QE
    tr = T / card['tnom']
    isT = card['IS'] * math.exp((tr - 1.0) * card['eg'] / vt
                                + card['xti'] * math.log(tr))
    ifwd = isT * (math.exp(vbe / (card['nf'] * vt)) - 1.0)
    irev = isT * (math.exp((vbe - vce) / (card['nr'] * vt)) - 1.0)
    ict = ifwd - irev
    ibe = ifwd / card['bf']
    ibc = irev / card['br']
    ## P = V(c,e)*Ic + V(b,e)*Ib, with Ic and Ib the STATIC terminal
    ## currents.  This is the same decomposition `GummelPoonNpnThermalHdl`
    ## uses, and it is written here from the topology rather than copied.
    return vce * (ict - ibc) + vbe * (ibe + ibc)


def _gp_dpdt(dT, h=1e-5):
    return (_gp_power(dT + h) - _gp_power(dT - h)) / (2.0 * h)


def _runaway_onset():
    """The saddle-node of ``dT = rth*P(dT)``, solved on the independent
    transcription alone.

    Thermal runaway is the loss of a solution, not a large number: the
    fixed point exists while the LOOP GAIN ``L = rth*dP/dT`` is below 1
    and disappears when it reaches it.  Eliminating ``rth`` between
    ``dT = rth*P`` and ``rth*P' = 1`` leaves ``dT*P'(dT) = P(dT)``, one
    equation in one unknown.  Returns ``(dT_crit, rth_crit)``.
    """
    lo, hi = 1e-6, 500.0

    def g(dT):
        return dT * _gp_dpdt(dT) - _gp_power(dT)
    assert g(lo) < 0.0 < g(hi)
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if g(m) < 0.0:
            lo = m
        else:
            hi = m
    dTc = 0.5 * (lo + hi)
    return dTc, dTc / _gp_power(dTc)


def _gp_fixed_point(rth, n=200000):
    """The self-consistent temperature rise, by damped iteration on the
    independent transcription.  Converges exactly while ``L < 1``."""
    dT = 0.0
    for _ in range(n):
        new = 0.5 * dT + 0.5 * rth * _gp_power(dT)
        if abs(new - dT) < 1e-13:
            return new
        dT = new
    return dT


def _th_circuit(rth, cth=0.0, vbe=TH_VBE, vce=TH_VCE, card=NPN_TH):
    c = SubCircuit()
    c['vc'] = VS('c', gnd, v=vce, vac=0.0)
    c['vb'] = VS('b', gnd, v=vbe, vac=0.0)
    c['q'] = eh.GummelPoonNpnThermalHdl('c', 'b', gnd, 'th', gnd,
                                        rth=rth, cth=cth, **card)
    return c


def test_thermal_bjt_at_rth_zero_is_the_isothermal_bjt():
    """``rth = 0`` collapses the thermal branch to a zero-volt source, so
    the device is `GummelPoonNpnHdl` again -- asserted on a real DC solve
    of a card with EVERY block on, to the solver's own tolerance.

    That equivalence is what makes the thermal variant a bolt-on rather
    than a fork: if it were only approximately the isothermal device, the
    library would have two Gummel-Poon models to keep in step.
    """
    full = dict(IS=1e-16, bf=120.0, nf=1.0, vaf=80.0, ikf=2e-2, ise=3e-15,
                ne=1.6, br=3.0, nr=1.0, var=20.0, ikr=5e-3, isc=1e-14,
                nc=1.9, cje=2e-13, vje=0.72, mje=0.35, cjc=1.2e-13,
                vjc=0.6, mjc=0.4, xcjc=0.6, tf=3e-11, tr=1e-9, fc=0.5,
                rb=30.0, rbm=8.0, re=1.0, rc=15.0, tnom=_T0)

    def solve(cls, **kw):
        c = SubCircuit()
        c['vc'] = VS('c', gnd, v=3.0, vac=0.0)
        c['vb'] = VS('b', gnd, v=0.78, vac=0.0)
        if cls is eh.GummelPoonNpnThermalHdl:
            c['q'] = cls('c', 'b', gnd, 'th', gnd, **dict(full, **kw))
        else:
            c['q'] = cls('c', 'b', gnd, **full)
        res = DC(c, toolkit=numeric).solve()
        return float(res.i('vc.plus')), float(res.i('vb.plus'))

    iso = solve(eh.GummelPoonNpnHdl)
    th0 = solve(eh.GummelPoonNpnThermalHdl, rth=0.0, cth=0.0)
    assert_allclose(th0, iso, rtol=1e-12)
    ## The card really is conducting, so this is not two zeros agreeing.
    assert abs(iso[0]) > 1e-4, iso
    ## and with rth ON the same device draws measurably more.
    hot = solve(eh.GummelPoonNpnThermalHdl, rth=2000.0, cth=0.0)
    assert abs(hot[0]) > 1.05 * abs(iso[0]), (hot, iso)


def test_thermal_bjt_steady_state_is_the_dissipation_times_rth():
    """``dT = rth*P`` at DC, with ``P`` rebuilt from the SOURCE currents
    of the solve rather than from anything the element reports.

    That is the whole content of the thermal node, and building ``P`` from
    the supply currents makes it a statement about the circuit: whatever
    the device dissipates has to come out of the two sources.
    """
    ## The last two rows are SATURATION -- both junctions forward biased
    ## -- and they are here because forward active does not exercise the
    ## reverse transport current at all.  Measured: deleting the `-IBC`
    ## term from the dissipation changes nothing at `Vce = 5 V` (it is
    ## 1e-16 A there) and is caught immediately at `Vce = 0.2 V`.
    for vce, vbe, rth in ((5.0, 0.7, 500.0), (5.0, 0.7, 2000.0),
                          (5.0, 0.7, 8000.0), (0.2, 0.8, 2000.0),
                          (0.05, 0.85, 500.0)):
        res = DC(_th_circuit(rth, vbe=vbe, vce=vce),
                 toolkit=numeric).solve()
        ic = -float(res.i('vc.plus'))
        ib = -float(res.i('vb.plus'))
        p = vce * ic + vbe * ib
        assert_allclose(float(res.v('th')), rth * p, rtol=1e-9,
                        err_msg='vce=%g rth=%g' % (vce, rth))
        assert p > 1e-6
    ## and the last point really is saturated, which is the physical
    ## statement that the REVERSE transport current is doing work there:
    ## the measured beta is a twentieth of `bf`, and it is the reverse
    ## current that took it there.
    res = DC(_th_circuit(500.0, vbe=0.85, vce=0.05), toolkit=numeric).solve()
    beta = abs(float(res.i('vc.plus')) / float(res.i('vb.plus')))
    assert beta < 0.1 * NPN_TH['bf'], beta


def test_thermal_bjt_matches_an_independent_fixed_point():
    """The self-consistent temperature rise against the independent
    transcription, from a fifth of the critical thermal resistance to a
    thousandth below it.

    One tolerance for the whole set.  The last point is where the loop
    gain is 0.999, i.e. where the fixed-point iteration on the reference
    takes tens of thousands of steps to converge -- and the Newton solve
    with an exact Jacobian lands on it directly, which is the difference
    between an outer iteration and a real unknown.
    """
    _, rc = _runaway_onset()
    for f in (0.2, 0.5, 0.9, 0.99, 0.999):
        got = float(DC(_th_circuit(rc * f), toolkit=numeric).solve().v('th'))
        assert_allclose(got, _gp_fixed_point(rc * f), rtol=1e-7,
                        err_msg='%g x rth_crit' % f)
    ## and the temperature rise really moves over that span.
    lo = float(DC(_th_circuit(rc * 0.2), toolkit=numeric).solve().v('th'))
    hi = float(DC(_th_circuit(rc * 0.999), toolkit=numeric).solve().v('th'))
    assert hi / lo > 10.0, (lo, hi)


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_thermal_runaway_onset_is_where_the_loop_gain_reaches_one():
    """THE POINT OF THE MODEL, and it is bracketed to 0.2%.

    Thermal runaway is the disappearance of the DC solution, and the
    criterion is ``rth*dP/dT = 1``.  ``rth_crit`` is computed here from
    the independent transcription ALONE -- the compiled model is never
    consulted -- and then the model is asked to solve at 0.999 and 1.001
    times it.  It converges at the first and has no solution at the
    second.

    The measured onset: ``dT_crit = 17.66 K``, ``rth_crit = 21190 K/W``
    at ``Vbe = 0.7 V``, ``Vce = 5 V`` on this card.
    """
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix

    dTc, rc = _runaway_onset()
    ## The criterion is satisfied to nine digits at the bracket -- that
    ## is the reference checking itself, not the model.
    assert_allclose(rc * _gp_dpdt(dTc), 1.0, rtol=1e-7)
    assert 15.0 < dTc < 20.0 and 1e4 < rc < 1e5, (dTc, rc)

    ## Below: a solution, and it is close to the critical temperature.
    below = float(DC(_th_circuit(rc * 0.999), toolkit=numeric).solve()
                  .v('th'))
    assert 0.9 * dTc < below < dTc, (below, dTc)

    ## Above: none at all.  Two decades of thermal resistance past the
    ## onset, and the rescue ladder and the gmin anchor do not
    ## manufacture one -- there is nothing to find.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for f in (1.001, 1.05, 1.5, 3.0, 100.0):
            with pytest.raises((SingularMatrix, NoConvergenceError)):
                DC(_th_circuit(rc * f), toolkit=numeric).solve()


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_thermal_runaway_is_a_transient_that_does_not_settle():
    """The same onset in TIME, which is what a designer actually sees.

    Started from ambient (``uic``), a sub-critical device settles at the
    fixed point and a super-critical one does not settle at all: it walks
    straight through the critical temperature and accelerates.

    The settling is checked as an EXPONENTIAL TAIL with the linearised
    time constant ``rth*cth/(1 - L)``, computed from the reference.  That
    is a sharper claim than "it reaches the right value": the plain
    ``rth*cth`` is 19 ms here and the measured tail is 50 ms, because the
    thermal feedback stretches it by ``1/(1 - 0.621)``.
    """
    _, rc = _runaway_onset()
    cth = 1e-6

    def run(f, tend, n=4000):
        c = _th_circuit(rc * f, cth=cth)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = Transient(c, toolkit=numeric, uic=True).solve(
                tend=tend, timestep=tend / n, fixed_timestep=True)
        return (np.asarray(res.v('th').x[0], float),
                np.asarray(res.v('th').y, float))

    ## Sub-critical: settles, from ambient, on the reference's number.
    rth = rc * 0.9
    dinf = _gp_fixed_point(rth)
    tau = rth * cth / (1.0 - rth * _gp_dpdt(dinf))
    t, y = run(0.9, 10.0 * tau)
    assert y[0] == 0.0
    assert_allclose(y[-1], dinf, rtol=1e-4)
    ## and the tail is that exponential, fitted over 3-8 time constants.
    m = (t > 3.0 * tau) & (t < 8.0 * tau)
    slope = np.polyfit(t[m], np.log(dinf - y[m]), 1)[0]
    assert_allclose(-1.0 / slope, tau, rtol=1e-2)
    ## The stretch is real: the bare thermal time constant is 2.6x
    ## shorter, so this is not `rth*cth` under another name.
    assert tau / (rth * cth) > 2.0

    ## Super-critical: no settling, and the SHAPE is the interesting
    ## part.  The device does not simply rise ever faster -- it slows
    ## down first, near where the fixed point used to be, and then
    ## accelerates.  That bottleneck is the "ghost" of the saddle-node
    ## the runaway condition destroyed, and it is what makes thermal
    ## runaway look survivable for a while and then not.
    ##
    ## Measured at 1.2x the critical thermal resistance: the rise rate
    ## bottoms out at 125 K/s at dT = 14.1 K -- 20% below the critical
    ## 17.66 K -- and is 293 K/s and still climbing at the end of the
    ## run, having passed 25.8 K.
    dTc, _ = _runaway_onset()
    tau12 = rc * 1.2 * cth
    t2, y2 = run(1.2, 6.0 * tau12)
    assert np.all(np.diff(y2) >= -1e-9)
    assert y2[-1] > 1.4 * dTc, (y2[-1], dTc)
    rate = np.diff(y2) / np.diff(t2)
    k = int(np.argmin(rate))
    assert 0 < k < len(rate) - 1, 'the minimum must be interior'
    assert 0.7 * dTc < y2[k] < 1.3 * dTc, (y2[k], dTc)
    assert np.all(np.diff(rate[k:]) > 0), 'no acceleration after the ghost'
    assert rate[-1] > 2.0 * rate[k], (rate[-1], rate[k])
    ## and the sub-critical device at the same time is going nowhere.
    assert y[-1] < 0.7 * y2[-1]


def test_the_thermal_pin_takes_an_external_foster_ladder():
    """The thermal node is a PIN, and the reason it is one is that a real
    package model is an R-C ladder somebody else builds.

    **The device's own ``rth`` is in PARALLEL with whatever is hung on
    the pin**, not in series with it, because `SelfHeating` puts its RC
    across the same ``(th, tha)`` branch the dissipation is injected
    into.  That is worth stating because the obvious way to say "no
    internal thermal resistance, use my ladder" -- ``rth = 0`` -- does
    the opposite: it collapses the branch to a zero-volt source and
    SHORTS the ladder out.  `SelfHeating`'s own docstring records that
    trap for the multi-device case; it is the same trap here.  What
    works is a large ``rth``, and the exact parallel value is what this
    asserts against.
    """
    rdev, r1, r2 = 1e6, 1500.0, 3500.0
    reff = 1.0 / (1.0 / rdev + 1.0 / (r1 + r2))
    c = SubCircuit()
    c['vc'] = VS('c', gnd, v=TH_VCE, vac=0.0)
    c['vb'] = VS('b', gnd, v=TH_VBE, vac=0.0)
    c['q'] = eh.GummelPoonNpnThermalHdl('c', 'b', gnd, 'th', gnd,
                                        rth=rdev, cth=0.0, **NPN_TH)
    c['r1'] = R('th', 'tj', r=r1)
    c['r2'] = R('tj', gnd, r=r2)
    c['c2'] = C('tj', gnd, c=1e-6)
    res = DC(c, toolkit=numeric).solve()
    lumped = float(DC(_th_circuit(reff), toolkit=numeric).solve().v('th'))
    assert_allclose(float(res.v('th')), lumped, rtol=1e-9)
    ## and the ladder's intermediate node carries its share, so the two
    ## stages really are in series with each other.
    assert_allclose(float(res.v('tj')) / float(res.v('th')),
                    r2 / (r1 + r2), rtol=1e-6)
    ## The ladder is doing the work: the device's own path carries 2.5%
    ## of the heat here, which is what `reff` being 4975 rather than
    ## 5000 says.
    assert 0.98 < reff / (r1 + r2) < 1.0


@pytest.mark.parametrize('name,rth,x', [
    ('cold', 0.0, [5.0, 0.75, 0.0, 0.0, 0.0, 0.0]),
    ('warm', 2000.0, [5.0, 0.75, 0.0, 12.0, 0.0]),
    ('near runaway', 20000.0, [5.0, 0.7, 0.0, 17.0, 0.0]),
    ('reverse', 2000.0, [-2.0, 0.3, 0.0, 1.0, 0.0]),
    ('absurd temperature', 2000.0, [5.0, 0.7, 0.0, -400.0, 0.0]),
])
def test_thermal_bjt_jacobians_by_finite_differences(name, rth, x):
    """``G`` and ``C`` by central differences, including one point below
    absolute zero -- which a Newton iterate really can visit, because the
    temperature rise is a solution variable.  That is what `SelfHeating`'s
    smooth floor exists for, and it is the point at which a hard `Max`
    would have made the Jacobian ``0/0``."""
    el = _mk(eh.GummelPoonNpnThermalHdl, 'c', 'b', 'e', 'th', 'tha',
             rth=rth, cth=1e-6, **NPN_TH)
    x = list(x)[:len(x_layout(el))]
    res = check_jacobians(el, np.array(x, dtype=float), rtol=3e-5)
    assert res.ok, '%s\n%s' % (name, res)


def test_a_self_heating_limiter_parameter_is_refused_by_name():
    """⚠ **The finding this model produced**, kept as a test because it
    is a rule a self-heating device meets and an isothermal one does not.

    `limit_pnj`'s parameters are evaluated from the card BEFORE the
    device is, so they may not reach the solution -- and with a thermal
    node the temperature-scaled saturation current DOES: ``isT`` depends
    on ``T``, which is ``$temperature + V(th,tha)``.  The natural
    spelling, the one the isothermal `GummelPoonNpnHdl` uses and the one
    roadmap 12.4 went to some trouble to make legal, is refused at
    compile time -- correctly, and by name.

    `GummelPoonNpnThermalHdl` therefore passes the AMBIENT temperature to
    the limiter (`_gp_core`'s ``tlim``).  That places the critical
    voltage against a junction that may be a hundred kelvin hotter than
    the limiter thinks; the same trade `limit_fet` records for its
    parameter-only ``vto``, and looser rather than wrong.
    """
    from pycircuit.circuit.hdl import (Branch, Contribution, var, expl,
                                       limit_pnj, vt)
    from pycircuit.utilities.param import Parameter

    def build():
        class _SelfHeatedLimiter(Behavioural):
            instparams = [Parameter(name='IS0', desc='saturation current',
                                    unit='A', default=1e-14),
                          Parameter(name='rth', desc='thermal resistance',
                                    unit='K/W', default=100.0),
                          Parameter(name='cth', desc='thermal capacitance',
                                    unit='J/K', default=0.0)]

            @staticmethod
            def analog(plus, minus, th, tha):
                heat = eh.SelfHeating(th, tha, rth, cth)       # noqa: F821
                b = Branch(plus, minus)
                ## `isT` reaches the thermal node, which is the whole
                ## point of a self-heating device -- and makes it
                ## illegal as a limiter parameter.
                isT = var(IS0 * expl(heat.T / 300.0), 'isT')   # noqa: F821
                v = var(limit_pnj(b.V, isT, vt(heat.T)), 'vd')
                i = var(isT * (expl(v / vt(heat.T)) - 1.0), 'id')
                return ((Contribution(b.I, i),)
                        + heat.dissipate(var(b.V * i, 'p')))
        return _SelfHeatedLimiter

    with pytest.raises(ValueError, match='cannot depend on the solution'):
        build()
    ## and the message names the offending terminals, which is what makes
    ## this a five-minute fix rather than a bisection.
    with pytest.raises(ValueError, match=r"reaches th \(through var"):
        build()

    ## The shipped model does compile, declares its limiter on both
    ## junctions, and its limiter is INDEPENDENT of the thermal node --
    ## which is the observable consequence of passing the ambient
    ## temperature, and the thing a future change would break silently.
    el = _mk(eh.GummelPoonNpnThermalHdl, 'c', 'b', 'e', 'th', 'tha',
             rth=1000.0, **NPN_TH)
    assert [s[1] for s in el._hdl_info['limit_spec']] == ['pnj', 'pnj']
    x0 = np.array([1.0, 0.7, 0.0, 0.0, 0.0])
    x = np.array([1.0, 5.0, 0.0, 0.0, 0.0])
    cold = el.limit(x, x0)
    for dT in (50.0, 200.0):
        xh = x.copy()
        x0h = x0.copy()
        xh[3] = x0h[3] = dT
        hot = el.limit(xh, x0h)
        assert (hot[:3] == cold[:3]).all(), (dT, hot, cold)
    ## and the limiter did bite, so this is not two no-ops agreeing.
    assert not (cold == x).all()
