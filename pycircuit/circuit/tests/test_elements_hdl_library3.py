# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The third Phase-6 library batch: the two semiconductor models the
repo could not previously build.

=====================  ======================================================
model                  what it is, and why it is here
=====================  ======================================================
`GummelPoonNpnHdl`     the SPICE Gummel-Poon bipolar transistor -- the
`GummelPoonPnpHdl`     archetype every SPICE ships, and the most
                       conspicuous absence in the tree.  Also the first
                       model with TWO ``limit_pnj`` probes sharing a
                       terminal, which is the case ``limit_spec``'s
                       write-back rule exists for
`EkvNmosHdl`           EKV 2.6's long-channel core -- one tenth of PSP's
`EkvPmosHdl`           size, all-region from ONE interpolation function,
                       the first production user of ``limit_fet``,
                       ``limit_vds`` and ``softplus``, and an
                       INDEPENDENT formulation of the physics
                       `compact.PspMosLongChannel` solves exactly
=====================  ======================================================

**Every reference here is external to the model under test.**  In
order of how much they can catch:

* an **independent numpy transcription** of the Gummel-Poon equations,
  written from Massobrio & Antognetti, *Semiconductor Device Modeling
  with SPICE*, 2nd ed., ch. 2 (the integral charge-control model of
  Gummel & Poon, BSTJ **49**, 827, 1970), living in `_gp_reference`
  below.  It shares no code with the DSL model -- not the temperature
  path, not the depletion charge, not the base-charge factor;
* **identities derived by hand from those equations and written next to
  the number**, which is the strongest kind because they cannot be
  satisfied by a transcription error that both sides share.  The two
  that carry the most weight are ``beta = BF/2`` at a collector current
  of exactly ``IKF`` (the definition of the forward knee) and
  ``Ic/go = VAF + VCE - VBE`` in forward active (the Early effect as a
  measured slope, never as a parameter read back);
* the **textbook asymptotes of the EKV interpolation function** -- a
  subthreshold slope of ``ln(10)*n*UT`` per decade below threshold, the
  square law above it, and the Ward-Dutton charge partition's 50/50 at
  ``Vds = 0`` and 40/60 in saturation, all of which are properties of
  the physics and not of this implementation;
* an **independent RK4 integration** of the two-node ODE of a
  common-emitter stage, from the same numpy transcription, for the
  transient;
* a **cross-check between two different formulations** -- EKV against
  `compact.PspMosLongChannel`, whose drain current is validated to
  1.3e-6 against a vendor's compiled PSP103.  Its result is reported in
  section 3 whatever it turns out to be, and it turns out to be
  interesting: the two agree on the SHAPE of the characteristic to a
  few percent and disagree on the absolute threshold by 160 mV between
  weak and strong inversion.  Nothing is tuned to make them agree;
* **finite differences** for every Jacobian claim, through
  `hdl.check_jacobians`.
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
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.elements import SubCircuit, VS, R
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.hdl import (Node, check_jacobians, x_layout, explain,
                                   KBOLTZMANN, QELECTRON)
from pycircuit.circuit._limiting import _pnjlim, _fetlim, _limvds

#: The tree's own constants, so the reference below computes the same
#: thermal voltage the compiled model does.  Taking them from the module
#: rather than writing 1.38e-23 here is deliberate: the reference is
#: meant to be independent of the MODEL, not of the simulator's physical
#: constants, and a mismatch there would show up as a spurious 0.9%
#: (measured, while writing this file: `defaultepar.T` is 300 K and the
#: card's `tnom` is 300.15, so a reference that assumed T = tnom was off
#: by exactly the temperature path plus the thermal voltage).
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
## 0.  The independent Gummel-Poon reference
## ======================================================================

#: A card with EVERY block switched on, so that no test can pass by a
#: term being zero.  Values are of the order a small-signal n-p-n's
#: would be; nothing here is fitted to anything.
NPN = dict(IS=1e-16, bf=120.0, nf=1.0, vaf=80.0, ikf=2e-2,
           ise=3e-15, ne=1.6, br=3.0, nr=1.0, var=20.0, ikr=5e-3,
           isc=1e-14, nc=1.9, cje=2e-13, vje=0.72, mje=0.35,
           tf=3e-11, xtf=2.0, vtf=2.0, itf=0.03,
           cjc=1.2e-13, vjc=0.6, mjc=0.4, xcjc=0.6, tr=1e-9, fc=0.5)

#: The same device with the three parasitic resistances given, so the
#: internal nodes exist and the collapse variants differ.
NPN_R = dict(NPN, rb=30.0, rbm=8.0, re=1.0, rc=15.0)

#: The IDEAL transport device: no leakage, no knee, no Early, no charge.
#: Used where a hand-derived identity needs the other blocks out of the
#: way -- and every such test says which blocks it switched off and why.
NPN_IDEAL = dict(IS=1e-16, bf=100.0, br=1.0)


def _dep_q(v, cj, vj, m, fc):
    """SPICE's depletion charge, transcribed from Massobrio &
    Antognetti eq. (2.62)-(2.64).  Scalar, branchy, and with no
    regularisation at all -- which is what makes it a reference for the
    both-arms-safe form in the model."""
    m = min(m, 0.9)
    fc = min(fc, 0.9999)
    if v < fc * vj:
        return cj * vj * (1.0 - (1.0 - v / vj) ** (1.0 - m)) / (1.0 - m)
    f1 = vj * (1.0 - (1.0 - fc) ** (1.0 - m)) / (1.0 - m)
    f2 = (1.0 - fc) ** (1.0 + m)
    f3 = 1.0 - fc * (1.0 + m)
    fcvj = fc * vj
    return cj * (f1 + (f3 * (v - fcvj)
                       + m / (2.0 * vj) * (v * v - fcvj ** 2)) / f2)


def _gp_reference(vbe, vbc, T=None, IS=1e-16, bf=100.0, nf=1.0, vaf=0.0,
                  ikf=0.0, ise=0.0, ne=1.5, br=1.0, nr=1.0, var=0.0,
                  ikr=0.0, isc=0.0, nc=2.0, rb=0.0, rbm=-1.0, re=0.0,
                  rc=0.0, cje=0.0, vje=0.75, mje=0.33, tf=0.0, xtf=0.0,
                  vtf=0.0, itf=0.0, cjc=0.0, vjc=0.75, mjc=0.33,
                  xcjc=1.0, tr=0.0, fc=0.5, xtb=0.0, eg=1.11, xti=3.0,
                  kf=0.0, af=1.0, area=1.0, tnom=300.15):
    """The Gummel-Poon model in plain numpy, from the textbook.

    Takes the two INTERNAL junction voltages (behind ``rb``/``re``/``rc``)
    and returns the terminal currents and the two stored charges.  Every
    formula is written the way the book writes it -- `exp`, not `expl`;
    `if`, not `Piecewise`; no clamps anywhere -- so that agreement with
    the compiled model is evidence about the model and not about a
    shared habit.
    """
    T = _T0 if T is None else T
    vt = _KB * T / _QE
    trat = T / tnom
    ltr = math.log(trat)
    factlog = (trat - 1.0) * eg / vt + xti * ltr
    isT = area * IS * math.exp(factlog)
    bfac = math.exp(xtb * ltr)
    bfT, brT = bf * bfac, br * bfac
    iseT = area * ise * math.exp(factlog / ne) / bfac
    iscT = area * isc * math.exp(factlog / nc) / bfac
    egT = 1.16 - 7.02e-4 * T ** 2 / (T + 1108.0)
    egn = 1.16 - 7.02e-4 * tnom ** 2 / (tnom + 1108.0)
    vjeT = vje * trat - 3.0 * vt * ltr - egn * trat + egT
    vjcT = vjc * trat - 3.0 * vt * ltr - egn * trat + egT
    cjeT = area * cje * (1.0 + mje * (4e-4 * (T - tnom)
                                      - (vjeT / vje - 1.0)))
    cjcT = area * cjc * (1.0 + mjc * (4e-4 * (T - tnom)
                                      - (vjcT / vjc - 1.0)))

    ifw = isT * (math.exp(vbe / (nf * vt)) - 1.0)
    irv = isT * (math.exp(vbc / (nr * vt)) - 1.0)
    rvaf = 0.0 if vaf <= 0 else 1.0 / vaf
    rvar = 0.0 if var <= 0 else 1.0 / var
    rikf = 0.0 if ikf <= 0 else 1.0 / (area * ikf)
    rikr = 0.0 if ikr <= 0 else 1.0 / (area * ikr)
    q1 = 1.0 / (1.0 - vbc * rvaf - vbe * rvar)
    q2 = ifw * rikf + irv * rikr
    qb = q1 * 0.5 * (1.0 + math.sqrt(1.0 + 4.0 * q2))
    ict = (ifw - irv) / qb
    ibe = ifw / bfT + iseT * (math.exp(vbe / (ne * vt)) - 1.0)
    ibc = irv / brT + iscT * (math.exp(vbc / (nc * vt)) - 1.0)
    ic = ict - ibc
    ib = ibe + ibc

    den = ifw + area * itf
    tfr = 0.0 if den == 0.0 else ifw / den
    tfx = 1.0 if vtf <= 0 else math.exp(vbc / (1.44 * vtf))
    tff = tf * (1.0 + xtf * tfr * tfr * tfx)
    qbe = tff * ifw + _dep_q(vbe, cjeT, vjeT, mje, fc)
    qjc = _dep_q(vbc, cjcT, vjcT, mjc, fc)
    rbmx = rb if rbm < 0 else rbm
    return dict(ic=ic, ib=ib, ie=-(ic + ib), ict=ict, ifwd=ifw, irev=irv,
                qb=qb, ibe=ibe, ibc=ibc, qbe=qbe,
                qbc=tr * irv + xcjc * qjc, qbx=(1.0 - xcjc) * qjc,
                rbb=(rbmx + (rb - rbmx) / qb) / area, isT=isT, vt=vt,
                cjeT=cjeT, cjcT=cjcT, vjeT=vjeT, vjcT=vjcT, tff=tff)


## ======================================================================
## 1.  The Gummel-Poon BJT
## ======================================================================

def _bjt_rows(vc, vb, ve, vbi, vci, vei, **card):
    """The six residual rows of a `GummelPoonNpnHdl` with parasitics,
    assembled BY HAND from `_gp_reference` and Kirchhoff's law.

    Written out rather than taken from the element, because the point is
    to check the compiler's term ROUTING -- which node each contribution
    lands on -- and not only the arithmetic inside it.  ``i[k]`` is the
    current LEAVING node k into the device.
    """
    r = _gp_reference(vbi - vei, vbi - vci, **card)
    area = card.get('area', 1.0)
    irc = (vc - vci) * area / card['rc']
    ire = (ve - vei) * area / card['re']
    irb = (vb - vbi) / r['rbb']
    i = np.array([irc, irb, ire,
                  -irb + r['ibe'] + r['ibc'],
                  -irc - r['ibc'] + r['ict'],
                  -ire - r['ibe'] - r['ict']])
    ## `qbx` -- the (1 - xcjc) fraction of the collector capacitance --
    ## hangs off the EXTERNAL base, so its voltage is `V(b) - V(ci)` and
    ## NOT the internal junction voltage.  Writing it with the internal
    ## one instead is a factor-of-two error at small forward bias, which
    ## is how this line got written correctly: the first version used
    ## `vbc` and the sweep found it.
    qbx = (1.0 - card.get('xcjc', 1.0)) * _dep_q(
        vb - vci, r['cjcT'], r['vjcT'], card['mjc'], card['fc'])
    q = np.array([0.0, qbx, 0.0, r['qbe'] + r['qbc'],
                  -r['qbc'] - qbx, -r['qbe']])
    return i, q, r


def test_bjt_matches_the_textbook_transcription():
    """Currents and charges against `_gp_reference`, at 40 bias points
    covering forward active, saturation, cutoff and reverse.

    The comparison is on EVERY ROW of ``i`` and ``q``, with the
    parasitic resistances given so that the internal nodes exist: that
    makes it a test of the compiler's routing -- which node each
    contribution lands on -- as well as of the arithmetic.
    """
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **NPN_R)
    assert [nm for _k, nm, _kd in x_layout(el)] == \
        ['c', 'b', 'e', 'bi', 'ci', 'ei']
    worst_i = worst_q = 0.0
    for vbe in (-0.5, 0.2, 0.5, 0.7, 0.8, 0.9):
        for vce in (-1.0, 0.05, 0.2, 0.4, 1.0, 3.0, 10.0):
            ## Internal nodes placed a little away from the terminals,
            ## so the resistor rows carry a current and cannot pass by
            ## being identically zero.
            vbi, vci, vei = vbe - 0.01, vce - 0.02, 0.005
            x = np.array([vce, vbe, 0.0, vbi, vci, vei])
            want, qw, _r = _bjt_rows(vce, vbe, 0.0, vbi, vci, vei, **NPN_R)
            got = np.asarray(el.i(x), float)
            scale = max(np.abs(want).max(), 1e-14)
            worst_i = max(worst_i, np.abs(got - want).max() / scale)
            q = np.asarray(el.q(x), float)
            qs = max(np.abs(qw).max(), 1e-20)
            worst_q = max(worst_q, np.abs(q - qw).max() / qs)
    ## 1e-12 and not "close": both sides are double-precision evaluations
    ## of the same real number through different expression trees, so
    ## anything above a few ulps of accumulated rounding is a real
    ## disagreement.  Measured worst: ~7e-13.
    assert worst_i < 1e-11, worst_i
    assert worst_q < 1e-11, worst_q


def test_bjt_area_scales_the_currents_and_not_the_voltages():
    """``area`` multiplies every current and capacitance and divides
    every resistance, which is what makes it a geometric scale factor
    rather than a fitting knob: two devices of area 1 in parallel must
    be one device of area 2, exactly."""
    one = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **NPN_R)
    two = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **dict(NPN_R, area=2.0))
    x = np.array([2.0, 0.78, 0.0, 0.77, 1.99, 0.001])
    ## The internal-node voltages are the same for both only because the
    ## resistor currents scale with area too; that is the claim.
    assert_allclose(np.asarray(two.i(x), float),
                    2.0 * np.asarray(one.i(x), float), rtol=1e-12)
    assert_allclose(np.asarray(two.q(x), float),
                    2.0 * np.asarray(one.q(x), float), rtol=1e-12)


def test_bjt_collapse_removes_the_internal_nodes():
    """With SPICE's defaults (``rb = re = rc = 0``) the three internal
    nodes are not merely shorted, they are never compiled -- so the
    element has three unknowns and the ``1/rb`` never exists."""
    bare = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **NPN)
    assert [nm for _k, nm, _kd in x_layout(bare)] == ['c', 'b', 'e']
    ## and it is the SAME device: give the parasitics a vanishing value,
    ## put the internal nodes exactly on their terminals, and the six
    ## rows fold onto the three.
    ##
    ## Folding, not slicing.  With the internal node AT the terminal the
    ## resistor carries no current, so every terminal row is exactly
    ## zero and the whole device current sits on the internal row -- the
    ## first version of this test sliced `[:3]` and compared three zeros
    ## against the right answer.  The statement that actually holds is
    ## KCL over each shorted pair.
    withr = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e',
                **dict(NPN, rb=1e-9, re=1e-9, rc=1e-9))
    x3 = np.array([2.0, 0.8, 0.0])
    x6 = np.array([2.0, 0.8, 0.0, 0.8, 2.0, 0.0])
    i6 = np.asarray(withr.i(x6), float)
    q6 = np.asarray(withr.q(x6), float)
    ## the pairs are (c, ci), (b, bi), (e, ei) -- rows 0/4, 1/3, 2/5.
    fold = np.array([i6[0] + i6[4], i6[1] + i6[3], i6[2] + i6[5]])
    foldq = np.array([q6[0] + q6[4], q6[1] + q6[3], q6[2] + q6[5]])
    assert_allclose(np.asarray(bare.i(x3), float), fold, rtol=1e-10,
                    atol=1e-18)
    assert_allclose(np.asarray(bare.q(x3), float), foldq, rtol=1e-10,
                    atol=1e-22)


def test_rbm_unset_is_rb_and_rbm_zero_collapses_with_the_base_charge():
    """SPICE's ``RBM`` defaults to ``RB``, which a default VALUE cannot
    express: a card that deliberately says ``rbm = 0`` wants a base
    resistance that falls to nothing at high injection, and a card that
    says nothing wants a constant ``rb``.

    `hdl.param_given` is the instrument for that and the model was
    written with it; it had to be replaced by a negative sentinel
    because `$param_given` and `$limit` in one model raise a
    ``TypeError`` from inside the generated ``limit()``.  See the
    parameter's comment in `elements_hdl` -- this test pins the
    behaviour, whichever mechanism delivers it.
    """
    card = dict(NPN, rb=100.0, re=0.0, rc=0.0)
    unset = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **card)
    zero = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **dict(card, rbm=0.0))
    ## Hard on: qb is large, so `rbm + (rb - rbm)/qb` -> rbm.
    x = np.array([1.0, 0.95, 0.0, 0.94, 1.0, 0.0])
    r = _gp_reference(0.94, 0.94 - 1.0, **card)
    ## High injection: `qb` is 6.4 here, which is what makes the two
    ## readings of the card differ by a factor of 6.4 rather than by
    ## rounding.
    assert r['qb'] > 5.0, r['qb']
    ## Row 1 is the base terminal: its whole current is (vb - vbi)/rbb,
    ## with `rbb = rbm + (rb - rbm)/qb`.
    i_unset = float(np.asarray(unset.i(x), float)[1])
    i_zero = float(np.asarray(zero.i(x), float)[1])
    assert_allclose(i_unset, 0.01 / 100.0, rtol=1e-12)
    assert_allclose(i_zero, 0.01 * r['qb'] / 100.0, rtol=1e-12)
    ## The ratio IS qb -- that identity is the content of the parameter.
    assert_allclose(i_zero / i_unset, r['qb'], rtol=1e-12)


## ----------------------------------------------------------------------
## 1.1  The Gummel plot

def _gummel(el, vbe, vbc=0.0, n=6):
    """``(Ic, Ib)`` of a collapsed (no-parasitic) element at a given
    junction bias.  ``vbc = 0`` is where a Gummel plot is measured."""
    x = np.array([vbe - vbc, vbe, 0.0])
    i = np.asarray(el.i(x), float)
    return float(i[0]), float(i[1])


def test_gummel_plot_beta_falls_to_half_at_a_collector_current_of_ikf():
    """The forward knee current is DEFINED by where the measured beta
    has fallen to ``BF/2``, and the Gummel-Poon algebra puts that at a
    collector current of exactly ``IKF``.  Derivation, at ``Vbc = 0``
    and with the Early voltages infinite so ``q1 = 1``:

        qb = (1 + sqrt(1 + 4*IF/IKF))/2
        Ic = IF/qb,   Ib = IF/BF,   beta = BF/qb

        beta = BF/2  <=>  qb = 2  <=>  sqrt(1 + 4*IF/IKF) = 3
                     <=>  IF = 2*IKF  <=>  Ic = IF/qb = IKF

    So the test is not "the roll-off is somewhere near ikf": the knee is
    at ``Ic = IKF`` to whatever precision the solve has.  ``ise`` and
    ``isc`` are zero here because a leakage term would add a base
    current the identity does not contain -- that is tested separately.
    """
    ikf = 5e-3
    bf = 150.0
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e',
             IS=1e-16, bf=bf, br=2.0, ikf=ikf)
    ## Locate the bias where the MEASURED beta is BF/2, by bisection on
    ## the model's own output -- nothing is read back from the card.
    def f(vbe):
        ic, ib = _gummel(el, vbe)
        return ic / ib - bf / 2.0
    lo, hi = 0.3, 1.2
    assert f(lo) > 0 > f(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    ic_knee, ib_knee = _gummel(el, 0.5 * (lo + hi))
    assert_allclose(ic_knee, ikf, rtol=1e-9)
    ## and the two ends of the same curve: beta -> BF far below the
    ## knee, and beta*sqrt(Ic) -> constant far above it (the -1/2 slope
    ## that high-level injection produces, and the reason the knee is
    ## visible at all).
    ic_lo, ib_lo = _gummel(el, 0.45)
    assert ic_lo < ikf / 1e4
    assert_allclose(ic_lo / ib_lo, bf, rtol=1e-4)
    ## Far above the knee, `qb -> sqrt(IF/IKF)` and `Ic = IF/qb`, so
    ## `IF = Ic^2/IKF`, `qb = Ic/IKF` and `beta = BF*IKF/Ic`.  The
    ## measured log-log slope of beta against COLLECTOR current is
    ## therefore -1, not the -1/2 that the slope against IF would be:
    ## a decade of collector current costs a decade of beta.  (The
    ## first version of this test asserted -1/2 and was wrong about
    ## which current the -1/2 belongs to.)
    ic1, ib1 = _gummel(el, 1.05)
    ic2, ib2 = _gummel(el, 1.15)
    assert ic1 > 30 * ikf
    slope = (math.log(ic2 / ib2) - math.log(ic1 / ib1)) \
        / (math.log(ic2) - math.log(ic1))
    assert -1.01 < slope < -0.99, slope
    ## and the asymptote itself, not only its slope.
    assert_allclose(ic2 / ib2, bf * ikf / ic2, rtol=2e-2)


def test_gummel_plot_base_current_ideality_is_ne_at_low_injection():
    """The other half of a Gummel plot: the base current has TWO slopes.

    Below about 0.4 V the ``ISE``/``NE`` recombination term dominates
    and ``d ln(Ib)/d Vbe`` is ``1/(NE*Vt)``; above it the ideal term
    ``IF/BF`` takes over and the slope is ``1/(NF*Vt)``.  The measured
    ideality factor is what a Gummel plot is READ for, and it is a
    property of the base current alone -- the collector current keeps
    the ideal slope throughout, which is the second assertion here.
    """
    ## `ise` and `ne` are chosen so the two terms cross at 0.6 V, which
    ## puts 0.30 V and 0.90 V two and a half decades either side of the
    ## crossover -- far enough that the measured ideality is the pure
    ## one to better than 0.5%, and the test says how far rather than
    ## hoping.
    ne, nf = 2.0, 1.0
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e',
             IS=1e-16, bf=100.0, ne=ne, ise=1e-13, nf=nf, br=2.0)

    def ideality(f, v, dv=2e-3):
        a, b = f(v - dv), f(v + dv)
        return 2 * dv / (_UT * math.log(b / a))

    ib = lambda v: _gummel(el, v)[1]
    ic = lambda v: _gummel(el, v)[0]
    ## The two terms never separate COMPLETELY, and neither is a pure
    ## exponential -- the ``- 1`` matters at 0.3 V, where ``exp(v/NE/Vt)``
    ## is only 331.  Differentiating ``ln(Ib)`` by hand,
    ##
    ##     Ib = ISt*(E1 - 1)/BF + ISEt*(E2 - 1)
    ##     d ln(Ib)/dv = [ISt*E1/(BF*NF) + ISEt*E2/NE] / Ib
    ##
    ## so the ideality a measurement returns is the ratio below.  That
    ## is a closed form derived from the equations, not a fit, and it
    ## is what makes the tolerance 1e-4 rather than 1%.  (The first
    ## version of this test left the ``- 1`` out, predicted 1.9935 and
    ## measured 1.9875 -- a real 0.3% that turned out to be physics.)
    card = dict(IS=1e-16, bf=100.0, ne=ne, ise=1e-13, nf=nf, br=2.0)

    def predicted(v):
        r = _gp_reference(v, 0.0, **card)
        vtT = r['vt']
        ist = r['isT']
        ## `ISEt` is `ISE*exp(factlog/NE)`, and `factlog` is recoverable
        ## from `isT/IS`.
        iset = card['ise'] * (ist / card['IS']) ** (1.0 / card['ne'])
        e1 = math.exp(v / (card['nf'] * vtT))
        e2 = math.exp(v / (card['ne'] * vtT))
        ibv = ist * (e1 - 1.0) / card['bf'] + iset * (e2 - 1.0)
        return ibv / (ist * e1 / (card['bf'] * card['nf'])
                      + iset * e2 / card['ne'])

    assert_allclose(ideality(ib, 0.30), predicted(0.30), rtol=1e-4)
    assert_allclose(ideality(ib, 0.90), predicted(0.90), rtol=1e-4)
    ## and the readings a Gummel plot is taken FOR: 2.0 at the bottom,
    ## 1.0 at the top.
    assert_allclose(ideality(ib, 0.30), ne, rtol=1.5e-2)
    assert_allclose(ideality(ib, 0.90), nf, rtol=1.5e-2)
    ## The collector current never sees either: NF at both ends.  The
    ## 2e-5 band is the ``- 1`` in ``IS*(exp(v/Vt) - 1)`` again, which
    ## at 0.30 V is one part in exp(11.6) = 1.1e5.
    assert_allclose(ideality(ic, 0.30), nf, rtol=2e-5)
    assert_allclose(ideality(ic, 0.90), nf, rtol=1e-9)
    ## And the consequence a designer cares about: beta is DEGRADED at
    ## low current, by exactly the ratio of the two base currents.
    icl, ibl = _gummel(el, 0.30)
    r = _gp_reference(0.30, 0.0, IS=1e-16, bf=100.0, ne=ne, ise=1e-13,
                      nf=nf, br=2.0)
    assert icl / ibl < 0.01 * 100.0
    assert_allclose(icl / ibl, r['ic'] / r['ib'], rtol=1e-10)


def test_early_effect_is_a_measured_output_conductance():
    """``VAF`` is verified by MEASURING ``go = dIc/dVce`` and reading the
    Early voltage back out of it, never by comparing against the card.

    In forward active with the knee and both leakage terms off,

        Ic = (IF + IS)*(1 + (Vce - Vbe)/VAF) + IS/BR
        go = (IF + IS)/VAF
        Ic/go - Vce + Vbe = VAF + (IS/BR)*VAF/(IF + IS)

    and the correction term is 8e-11 V on this card, so the identity is
    exact to every digit the difference quotient carries.  The test also
    checks the OTHER direction: with ``vaf`` infinite the same measured
    conductance collapses by six decades, which is what says the effect
    came from ``vaf`` and not from somewhere else in the model.
    """
    vaf = 80.0
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e',
             IS=1e-16, bf=120.0, br=2.0, vaf=vaf)
    for vbe in (0.6, 0.7, 0.75):
        for vce in (2.0, 5.0, 12.0):
            h = 1e-3
            icm = _gummel(el, vbe, vbe - (vce - h))[0]
            icp = _gummel(el, vbe, vbe - (vce + h))[0]
            ic = _gummel(el, vbe, vbe - vce)[0]
            go = (icp - icm) / (2 * h)
            assert_allclose(ic / go - vce + vbe, vaf, rtol=1e-6)
    ## Without vaf: the residual output conductance is the base-collector
    ## leakage alone, six decades down.
    flat = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e',
               IS=1e-16, bf=120.0, br=2.0, vaf=0.0)
    h = 1e-3
    ic = _gummel(flat, 0.7, 0.7 - 5.0)[0]
    go0 = (_gummel(flat, 0.7, 0.7 - 5.0 - h)[0]
           - _gummel(flat, 0.7, 0.7 - 5.0 + h)[0]) / (2 * h)
    assert abs(go0) < 1e-6 * ic / vaf


def test_reverse_beta_and_the_reverse_early_voltage():
    """Run the device backwards -- emitter as collector -- and the same
    two identities hold with ``BR`` and ``VAR`` in place of ``BF`` and
    ``VAF``.  That is what makes this a transport model rather than a
    forward-only one, and it is what a saturated transistor spends its
    time doing.
    """
    br, var = 3.0, 25.0
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e',
             IS=1e-16, bf=100.0, br=br, var=var)
    ## Reverse active: base-collector forward biased, base-emitter
    ## reverse biased.  Bias by (vbc, vec).
    def rev(vbc, vec):
        ## Collector at 0, base at `vbc` (forward), emitter at `vec`
        ## ABOVE the collector -- so the base-emitter junction is
        ## reverse biased by `vbc - vec` and the emitter is doing the
        ## collecting.  Getting this vector's signs wrong is the easiest
        ## mistake in the file and produced a `beta` of 78 the first
        ## time, from a device that was simply hard on in the forward
        ## direction.
        x = np.array([0.0, vbc, vec])
        i = np.asarray(el.i(x), float)
        return float(i[2]), float(i[1])      # into the emitter, base
    ie, ib = rev(0.7, 5.0)
    ## The reverse beta, MEASURED.  The Early correction is
    ## (1 + (Vec - Vbc)/VAR), so beta_measured is br times that.
    assert_allclose(ie / ib, br * (1.0 + (5.0 - 0.7) / var), rtol=2e-3)
    h = 1e-3
    go = (rev(0.7, 5.0 + h)[0] - rev(0.7, 5.0 - h)[0]) / (2 * h)
    assert_allclose(ie / go - 5.0 + 0.7, var, rtol=1e-5)


## ----------------------------------------------------------------------
## 1.2  Charge

def _cj_closed(v, cj, vj, m, fc):
    """``dQ/dv`` of SPICE's depletion charge, DIFFERENTIATED BY HAND
    rather than by differencing `_dep_q`.

    Below the seam ``C = cj*(1 - v/vj)^-m``, which is the textbook
    junction capacitance; above it the charge is a quadratic and
    ``C = cj*(f3 + m*v/vj)/f2`` is its slope, which is what the ``fc``
    linearisation exists to keep finite.  The two agree at ``v = fc*vj``
    -- that is what "C1 at the seam" means -- and that agreement is
    asserted below.
    """
    m = min(m, 0.9)
    fc = min(fc, 0.9999)
    if v < fc * vj:
        return cj * (1.0 - v / vj) ** (-m)
    f2 = (1.0 - fc) ** (1.0 + m)
    f3 = 1.0 - fc * (1.0 + m)
    return cj * (f3 + m * v / vj) / f2


def test_junction_capacitance_is_the_hand_differentiated_charge():
    """``C`` against the closed-form ``dQ/dv``, on both sides of the
    ``fc*vj`` seam and through it.

    ``xtf = 0`` here so that the transit time is a constant and the
    diffusion capacitance is exactly ``tf*d(IF)/dVbe`` -- with the
    ``xtf`` modulation on, ``tf_eff`` depends on the bias too and the
    closed form would have to reproduce the model's own arithmetic,
    which is not a reference.
    """
    card = dict(NPN, xtf=0.0, tr=0.0, xcjc=1.0, rb=0.0, re=0.0, rc=0.0)
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **card)
    for vbe in (-2.0, 0.0, 0.30, 0.35, 0.40, 0.50, 0.70, 0.85):
        for vbc in (-3.0, -0.5, 0.0, 0.25, 0.31, 0.55):
            r = _gp_reference(vbe, vbc, **card)
            x = np.array([vbe - vbc, vbe, 0.0])
            C = np.asarray(el.C(x), float)
            cbe = (_cj_closed(vbe, r['cjeT'], r['vjeT'], card['mje'],
                              card['fc'])
                   + card['tf'] * r['isT']
                   * math.exp(vbe / (card['nf'] * r['vt']))
                   / (card['nf'] * r['vt']))
            cbc = _cj_closed(vbc, r['cjcT'], r['vjcT'], card['mjc'],
                             card['fc'])
            ## Row/col 2 is the emitter: `q[e] = -qbe`, so `C[2,2]` is
            ## exactly the base-emitter capacitance.  Row/col 0 is the
            ## collector and carries the base-collector one.
            assert_allclose(C[2, 2], cbe, rtol=1e-9,
                            err_msg='vbe=%g vbc=%g' % (vbe, vbc))
            assert_allclose(C[0, 0], cbc, rtol=1e-9,
                            err_msg='vbe=%g vbc=%g' % (vbe, vbc))


def test_the_fc_linearisation_is_c1_at_the_seam():
    """The depletion charge is a different expression above ``fc*vj``,
    and the join is C1 by construction: the charge and its slope both
    continue.  Measured ACROSS the seam rather than at it, because a
    seam that is only continuous at the point itself is not what the
    linearisation claims.
    """
    card = dict(NPN, xtf=0.0, tf=0.0, tr=0.0, xcjc=1.0,
                rb=0.0, re=0.0, rc=0.0)
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **card)
    r = _gp_reference(0.0, 0.0, **card)
    seam = card['fc'] * r['vjeT']
    assert 0.3 < seam < 0.4, seam
    h = 1e-7
    qat = lambda v: float(np.asarray(
        el.q(np.array([-1.0, v, 0.0])), float)[2])
    cat = lambda v: float(np.asarray(
        el.C(np.array([-1.0, v, 0.0])), float)[2, 2])
    assert_allclose(qat(seam - h), qat(seam + h), rtol=1e-6)
    assert_allclose(cat(seam - h), cat(seam + h), rtol=1e-6)
    ## A NOTE THIS TEST WAS WRITTEN TO CATCH AND DID NOT.  The first
    ## version went on to assert that the SECOND derivative jumps, on
    ## the grounds that two different expressions cannot agree in more
    ## than they were built to.  They do: differentiating both arms,
    ##
    ##     below:  dC/dv = cj*m/vj*(1 - fc)^(-m-1)
    ##     above:  dC/dv = cj*m/(vj*(1 - fc)^(1+m))
    ##
    ## are the SAME number, so SPICE's linearisation is C2 in the
    ## charge and not merely C1.  It is C3 that breaks -- the low arm
    ## keeps curving and the high arm is a quadratic.  Measured
    ## difference in dC/dv across the seam: 1.9e-19 on 2.5e-13, i.e.
    ## nothing.  The claim below is the one that is actually true.
    assert_allclose((cat(seam + 3 * h) - cat(seam + h)) / (2 * h),
                    (cat(seam - h) - cat(seam - 3 * h)) / (2 * h),
                    rtol=1e-4)
    ## and the Piecewise is not decoration: 0.2 V above the seam the
    ## model follows the high arm and the low arm is 20% away.
    r2 = _gp_reference(0.0, 0.0, **card)
    v2 = seam + 0.35
    low_arm = r2['cjeT'] * (1.0 - v2 / r2['vjeT']) ** (-card['mje'])
    assert_allclose(cat(v2), _cj_closed(v2, r2['cjeT'], r2['vjeT'],
                                        card['mje'], card['fc']),
                    rtol=1e-9)
    assert abs(cat(v2) - low_arm) > 0.5 * cat(v2)


def test_xcjc_puts_part_of_the_collector_capacitance_outside_rb():
    """``xcjc`` distributes the base-collector capacitance along the
    base resistance: the ``xcjc`` fraction sits on the INTERNAL base and
    the rest on the external one.  With ``rb > 0`` those are different
    nodes, and the external fraction is visible as the whole of the
    external base's own charge row.
    """
    card = dict(NPN_R, xcjc=0.35, tr=0.0)
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **card)
    vb, vci = 0.75, 1.4
    x = np.array([1.5, vb, 0.0, 0.74, vci, 0.01])
    r = _gp_reference(0.74 - 0.01, 0.74 - vci, **card)
    C = np.asarray(el.C(x), float)
    ## Node b's charge is qbx and nothing else, so C[1,1] is the
    ## EXTERNAL fraction of the collector capacitance at V(b) - V(ci).
    want = (1.0 - card['xcjc']) * _cj_closed(vb - vci, r['cjcT'],
                                             r['vjcT'], card['mjc'],
                                             card['fc'])
    assert_allclose(C[1, 1], want, rtol=1e-9)
    ## and the internal fraction is on the internal base, at the
    ## INTERNAL junction voltage -- a different number, because the two
    ## capacitors see different voltages once rb carries a current.
    ## Their sum is the whole of cjc only when rb = 0.
    whole = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e',
                **dict(card, rb=0.0, re=0.0, rc=0.0))
    xw = np.array([1.5, 0.74, 0.0])
    Cw = np.asarray(whole.C(xw), float)
    rw = _gp_reference(0.74, 0.74 - 1.5, **card)
    assert_allclose(Cw[0, 0], _cj_closed(0.74 - 1.5, rw['cjcT'],
                                         rw['vjcT'], card['mjc'],
                                         card['fc']), rtol=1e-9)


## ----------------------------------------------------------------------
## 1.3  Jacobians

@pytest.mark.parametrize('name,card,x', [
    ('forward active', NPN_R, [3.0, 0.78, 0.0, 0.77, 2.98, 0.005]),
    ('saturation', NPN_R, [0.1, 0.85, 0.0, 0.83, 0.12, 0.01]),
    ('cutoff', NPN_R, [3.0, 0.0, 0.0, 0.0, 3.0, 0.0]),
    ('reverse', NPN_R, [0.0, 0.7, 3.0, 0.69, 0.01, 2.99]),
    ('collapsed', NPN, [2.0, 0.75, 0.0]),
    ## Right on the two depletion seams, where the charge changes arm.
    ('be seam', NPN, [2.0, 0.3459, 0.0]),
    ('bc seam', NPN, [0.4614, 0.75, 0.0]),
    ## High injection, where qb's square root is doing the work.
    ('high injection', NPN, [2.0, 0.95, 0.0]),
])
def test_bjt_jacobians_by_finite_differences(name, card, x):
    """``G`` against ``di/dx`` and ``C`` against ``dq/dx``, at every
    seam the model has: both depletion linearisations, the high-injection
    square root, and the four quadrants of operation."""
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **card)
    res = check_jacobians(el, np.array(x, dtype=float), rtol=2e-5)
    assert res.ok, '%s\n%s' % (name, res)


def test_bjt_stays_finite_where_no_device_belongs():
    """Every entry of ``i``, ``q``, ``G`` and ``C`` is finite at biases a
    Newton iterate reaches and a transistor never does.

    This is the property the both-arms rule exists for: the discarded
    arm of every conditional is evaluated anyway, so "the model is only
    used sensibly" is not an argument.  ``expl`` bounds the
    exponentials, `safe_pow` the depletion power and the two `maxc`
    floors the base-charge denominators.
    """
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **NPN_R)
    rng = np.random.default_rng(20260825)
    for _ in range(400):
        x = rng.uniform(-300.0, 300.0, 6)
        for f in (el.i, el.q):
            v = np.asarray(f(x), float)
            assert np.isfinite(v).all(), (x, v)
        for f in (el.G, el.C):
            v = np.asarray(f(x), float)
            assert np.isfinite(v).all(), (x, v)
    ## and the extremes, where every regulariser is past its seam.
    for s in (-1e6, -1e3, 1e3, 1e6):
        x = np.array([s, s, -s, s, -s, s], dtype=float)
        assert np.isfinite(np.asarray(el.i(x), float)).all(), s
        assert np.isfinite(np.asarray(el.G(x), float)).all(), s


## ----------------------------------------------------------------------
## 1.4  $limit on a device with TWO junctions sharing a terminal

def test_limit_pnj_write_back_on_two_junctions_sharing_a_base():
    """The case ``limit_spec``'s compile-time write-back rule exists for,
    on the first production model that has it.

    Both junctions have the internal base as their PLUS terminal.  A
    limiter that wrote every probe by moving its plus would have the
    second write undo the first, and the base-emitter voltage would come
    out unlimited; the rule therefore gives the second probe its MINUS
    terminal.  Asserted three ways: the compile-time choice, the
    arithmetic of the two sequential writes, and -- the property that
    actually matters -- that BOTH branches carry their own limited value
    at the end.
    """
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **NPN_R)
    spec = el._hdl_info['limit_spec']
    ## (rows, kind, moved row).  bi = 3, ci = 4, ei = 5.
    assert [(s[0], s[1], s[2]) for s in spec] == \
        [((3, 5), 'pnj', 3), ((3, 4), 'pnj', 4)]
    ## The parameters the compiler lambdified: (IS, VT) per probe, with
    ## VT carrying nf for the first junction and nr for the second.
    ## `hdl._args_of` builds the trailing argument list a compiled
    ## function expects; a `$limit` parameter function takes the
    ## parameters and TEMP and nothing after them.
    from pycircuit.circuit.hdl import _args_of
    allargs = _args_of(el, defaultepar)
    pargs = allargs[:len(el._hdl_info['paramnames']) + 1]
    args = [tuple(float(f(*pargs)) for f in s[3]) for s in spec]
    r = _gp_reference(0.7, 0.0, **NPN_R)
    for (isv, vtv), n in zip(args, (NPN_R['nf'], NPN_R['nr'])):
        assert_allclose(isv, r['isT'], rtol=1e-12)
        assert_allclose(vtv, n * r['vt'], rtol=1e-12)

    rng = np.random.default_rng(3)
    both_moved = 0
    for _ in range(400):
        x = rng.uniform(-20.0, 20.0, 6)
        x0 = rng.uniform(-20.0, 20.0, 6)
        out = el.limit(x, x0)
        ## WHICH node each probe moves is decided at RUNTIME now -- the
        ## terminal that drifted further from the last accepted point,
        ## with the shared base going to whichever junction is applying
        ## the larger correction.  So this asserts the property the
        ## compile-time rule existed to protect, rather than the
        ## particular pair of indices it happened to produce:
        ##
        ##   * no probe is undone -- the two corrections land on two
        ##     DIFFERENT rows, so each junction keeps its own limited
        ##     voltage;
        ##   * nothing outside this element's own three internal rows is
        ##     touched;
        ##   * and the result is bounded, which is the failure the
        ##     runtime choice was introduced to remove.
        movers = [i for i in range(6) if out[i] != x[i]]
        assert len(set(movers)) == len(movers)
        assert len(movers) <= 2, (movers, x, x0)
        assert set(movers) <= {3, 4, 5}, (movers, x, x0)
        for i in range(6):
            if i not in movers:
                assert out[i] == x[i]
        assert np.all(np.isfinite(out))
        assert np.max(np.abs(out)) <= np.max(np.abs(x)) + 40.0, (x, out)
        ## The external terminals are never this element's to move.
        assert (out[[0, 1, 2]] == x[[0, 1, 2]]).all()
        assert not np.shares_memory(out, x)
        if len(movers) == 2:
            both_moved += 1
    assert both_moved > 50, both_moved


def test_limit_pnj_is_a_no_op_near_a_solution_and_bounds_a_wild_step():
    """The two properties a limiter is FOR, measured on the element.

    A step of at most ``2*VT`` passes through exactly (`_pnjlim`'s
    escape, from Listing 1 of the PCNR paper), which is what lets
    "did limiting fire?" be a convergence signal; and a step to 40 V is
    compressed to a few tens of millivolts above the previous point,
    which is what stops Newton evaluating ``exp(40/Vt)``.
    """
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **NPN)
    ## Collapsed: b = 1, e = 2, c = 0, so the two probes are (1,2) and
    ## (1,0), moving b then c.
    assert [(s[0], s[1], s[2]) for s in el._hdl_info['limit_spec']] == \
        [((1, 2), 'pnj', 1), ((1, 0), 'pnj', 0)]
    vt = _gp_reference(0.7, 0.0, **NPN)['vt']
    for base in (0.6, 0.7, 0.75):
        x0 = np.array([2.0, base, 0.0])
        for d in (-1.5 * vt, -0.5 * vt, 0.0, 0.5 * vt, 1.5 * vt):
            x = np.array([2.0, base + d, 0.0])
            assert (el.limit(x, x0) == x).all(), (base, d)
    ## and a wild one is bounded, to SPICE's own number rather than
    ## merely to something finite: from a forward-biased point the
    ## update is compressed logarithmically,
    ## ``vlim = vold + VT*ln(1 + (vnew - vold)/VT)``, so a step to 40 V
    ## lands 190 mV above where it started.
    x0 = np.array([2.0, 0.7, 0.0])
    out = el.limit(np.array([2.0, 40.0, 0.0]), x0)
    assert_allclose(out[1], 0.7 + vt * math.log(1.0 + (40.0 - 0.7) / vt),
                    rtol=1e-12)
    assert 0.88 < out[1] < 0.90, out
    ## The same for the base-collector junction, driven hard into
    ## forward from a reverse-biased point -- the branch of `pnjlim`
    ## that starts from `vold <= 0` and lands on `VT*ln(vnew/VT)`.
    out = el.limit(np.array([-40.0, 0.7, 0.0]), x0)
    assert 0.0 < out[1] - out[0] < 1.0, out


## ----------------------------------------------------------------------
## 1.5  Noise

def test_bjt_noise_is_two_uncorrelated_shot_sources():
    """The textbook bipolar noise model: shot noise ``2*q*Ic`` in the
    collector current and ``2*q*Ib`` in the base current, UNCORRELATED,
    plus flicker noise on the base current and thermal noise on the
    three parasitic resistances.

    Uncorrelated is the substantive claim and it is testable: the two
    sources sit on different branches, ``(c,e)`` and ``(b,e)``, so the
    ``CY`` entry that couples the collector row to the base row must be
    exactly zero, while the emitter -- which carries both -- sees the
    sum.  A single source on the emitter branch would give a non-zero
    ``CY[0,1]``, so this distinguishes the right model from a plausible
    wrong one.
    """
    card = dict(NPN, kf=1e-14, af=1.2, rb=0.0, re=0.0, rc=0.0)
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **card)
    x = np.array([2.0, 0.75, 0.0])
    r = _gp_reference(0.75, 0.75 - 2.0, **card)
    f = 1e3
    CY = np.asarray(el.CY(x, 2 * np.pi * f), float)
    qe = float(QELECTRON)
    ## `ict - ibc` IS the collector current and `ibe + ibc` the base
    ## current, so the two PSDs are 2q times those.
    sc = 2 * qe * abs(r['ict'] - r['ibc'])
    sb = 2 * qe * abs(r['ibe'] + r['ibc']) \
        + card['kf'] * abs(r['ibe'] + r['ibc']) ** card['af'] / f
    assert_allclose(CY[0, 0], sc, rtol=1e-9)
    assert_allclose(CY[1, 1], sb, rtol=1e-9)
    assert CY[0, 1] == 0.0 and CY[1, 0] == 0.0
    assert_allclose(CY[2, 2], sc + sb, rtol=1e-9)
    assert_allclose(CY[0, 2], -sc, rtol=1e-9)
    assert_allclose(CY[1, 2], -sb, rtol=1e-9)
    ## The flicker term really is 1/f: halve the frequency and only the
    ## base-current PSD moves, by exactly two.
    CY2 = np.asarray(el.CY(x, np.pi * f), float)
    assert_allclose(CY2[0, 0], sc, rtol=1e-9)
    assert_allclose(CY2[1, 1] - 2 * qe * abs(r['ibe'] + r['ibc']),
                    2 * (sb - 2 * qe * abs(r['ibe'] + r['ibc'])),
                    rtol=1e-9)


def test_base_resistance_noise_follows_the_base_charge():
    """``rb`` is modulated by ``qb``, so its thermal noise is
    bias-dependent -- ``4*k*T/rbb``, not ``4*k*T/RB``.  That falls out
    of writing the noise on the resistance the model actually uses, and
    it is the reason this contribution is inside the model rather than a
    separate resistor.
    """
    card = dict(NPN, rb=200.0, rbm=20.0, re=5.0, rc=50.0, kf=0.0)
    el = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **card)
    kt4 = 4.0 * float(KBOLTZMANN) * float(defaultepar.T)
    seen = []
    for vbi in (0.4, 0.8, 0.95):
        x = np.array([2.0, vbi + 0.05, 0.0, vbi, 1.9, 0.001])
        r = _gp_reference(vbi - 0.001, vbi - 1.9, **card)
        CY = np.asarray(el.CY(x, 2 * np.pi * 1e3), float)
        ## Row 1 is the external base and its only noise source is the
        ## base resistance (the junction sources sit on internal nodes).
        assert_allclose(CY[1, 1], kt4 / r['rbb'], rtol=1e-9)
        seen.append(float(CY[1, 1]))
    ## Strictly increasing, and over this bias range the base resistance
    ## falls from RB = 200 to 44 ohm -- a factor of 4.5, which is the
    ## measured size of the effect and not a round number chosen to
    ## pass.  ``qb`` is 1.0, 1.17 and 7.4 at the three points.
    assert seen[0] < seen[1] < seen[2]
    assert seen[2] / seen[0] > 4.0, seen
    ## and the collector/emitter resistances are constant, as they must
    ## be -- they are ordinary resistors.
    for vbi in (0.4, 0.95):
        x = np.array([2.0, vbi + 0.05, 0.0, vbi, 1.9, 0.001])
        CY = np.asarray(el.CY(x, 2 * np.pi * 1e3), float)
        assert_allclose(CY[0, 0], kt4 / card['rc'], rtol=1e-9)
        assert_allclose(CY[2, 2], kt4 / card['re'], rtol=1e-9)


## ----------------------------------------------------------------------
## 1.6  Real solves

def test_dc_operating_point_solves_the_reference_equations():
    """A common-emitter stage, solved by the DC analysis, and then the
    SOLUTION checked against `_gp_reference` and Kirchhoff's law.

    This is the strongest form the check can take: it does not compare
    two evaluations of the model, it takes the voltages the solver
    returned and asks whether the independent equations are satisfied
    there.  A model that converged to a point of its OWN equations that
    is not a point of the textbook's would fail here and nowhere else.
    """
    card = dict(NPN, rb=50.0, re=2.0, rc=100.0)
    rb_ext, rc_ext, vcc, vin = 47e3, 4.7e3, 5.0, 1.0
    c = SubCircuit()
    nb = c.add_node('nb')
    nc = c.add_node('nc')
    nvi = c.add_node('nvi')
    nvcc = c.add_node('nvcc')
    c['vcc'] = VS(nvcc, gnd, v=vcc)
    c['vin'] = VS(nvi, gnd, v=vin)
    c['rb'] = R(nvi, nb, r=rb_ext)
    c['rc'] = R(nvcc, nc, r=rc_ext)
    c['q1'] = eh.GummelPoonNpnHdl(nc, nb, gnd, **card)
    res = DC(c).solve()
    vb, vcn = float(res.v(nb, gnd)), float(res.v(nc, gnd))
    ## Forward active, and not by accident -- assert it, because a
    ## solution in cutoff or in saturation would satisfy the equations
    ## just as well and would test much less of the model.
    assert 0.6 < vb < 0.8, vb
    assert vb + 0.3 < vcn < vcc, vcn
    ## The device's own sub-vector, through the circuit's node map --
    ## which is what the stamping path uses, so it cannot drift.
    sub = np.asarray(res.x, float)[c.elementnodemap['q1']]
    vc_t, vb_t, ve_t, vbi, vci, vei = [float(v) for v in sub]
    assert_allclose(vc_t, vcn, rtol=1e-12)
    assert_allclose(vb_t, vb, rtol=1e-12)
    ## Kirchhoff at the two external nodes, evaluated with the
    ## reference and the solver's own voltages.
    want, _qw, _r = _bjt_rows(vc_t, vb_t, ve_t, vbi, vci, vei, **card)
    assert_allclose((vin - vb_t) / rb_ext, want[1], rtol=1e-6)
    assert_allclose((vcc - vc_t) / rc_ext, want[0], rtol=1e-6)
    ## and the three internal rows are satisfied to the solver's
    ## tolerance: nothing external touches them, so their residual --
    ## computed from the INDEPENDENT reference -- must be zero.
    assert np.abs(want[3:]).max() < 1e-9 * abs(want[0]), want


def test_dc_current_mirror_ratio_is_the_base_current_error():
    """A current mirror, whose output error is a CLOSED FORM: two
    identical devices sharing a base take ``2*Ic/beta`` out of the input
    current, so ``Iout/Iin = beta/(beta + 2)`` before the Early effect,
    and the Early effect multiplies that by ``1 + (Vce2 - Vce1)/VAF``.

    Both terms are predicted from the card and the solved bias --
    ``beta`` is MEASURED at the input device's operating point, never
    read back from ``bf`` -- and each is several percent, so the
    prediction is doing work.
    """
    card = dict(NPN_IDEAL, bf=80.0, vaf=60.0)
    vout = 3.0
    c = SubCircuit()
    nb = c.add_node('nb')
    no = c.add_node('no')
    nvi = c.add_node('nvi')
    c['vin'] = VS(nvi, gnd, v=5.0)
    c['rin'] = R(nvi, nb, r=4.25e3)
    c['vo'] = VS(no, gnd, v=vout)
    c['q1'] = eh.GummelPoonNpnHdl(nb, nb, gnd, **card)   # diode-connected
    c['q2'] = eh.GummelPoonNpnHdl(no, nb, gnd, **card)
    res = DC(c).solve()
    vb = float(res.v(nb, gnd))
    assert 0.6 < vb < 0.85, vb
    r1 = _gp_reference(vb, 0.0, **card)             # Vce1 = Vbe
    r2 = _gp_reference(vb, vb - vout, **card)
    beta = r1['ic'] / r1['ib']
    ## Sign conventions, pinned rather than guessed: `i(<el>.plus)`
    ## is the current INTO that terminal, so the input current is
    ## positive through `rin` and the mirror's output current comes
    ## back out of the source that holds the output node.
    iin = float(res.i('rin.plus'))
    iout = float(-res.i('vo.plus'))
    assert_allclose(iin, (5.0 - vb) / 4.25e3, rtol=1e-12)
    assert 0.5e-3 < iin < 2e-3, iin
    ## The two devices, checked against the reference at the solved
    ## bias, before anything is said about their ratio.
    assert_allclose(iin, r1['ic'] + 2 * r1['ib'], rtol=1e-6)
    assert_allclose(iout, r2['ic'], rtol=1e-6)
    ## The closed form.
    predicted = (beta / (beta + 2.0)) * (1.0 + (vout - vb) / card['vaf'])
    assert_allclose(iout / iin, predicted, rtol=1e-4)
    ## and both effects are big enough to be worth predicting: ~2.5% of
    ## base-current loss and ~3.8% of Early gain, in opposite directions.
    assert 0.02 < 2.0 / beta < 0.03
    assert 0.03 < (vout - vb) / card['vaf'] < 0.05


def test_dc_differential_pair_follows_the_exponential_ratio_law():
    """The identity a bipolar differential pair is used for: with the
    emitters tied and the collectors at equal potential,

        Ic1/Ic2 = exp(Vid/Vt)

    EXACTLY -- no beta, no tail current and no device parameter appears
    in it.  It holds because both devices see the same ``qb`` and the
    same base-collector bias, so everything but the two forward
    exponentials cancels.

    The Early voltages are left "infinite" here for exactly that
    reason: with ``vaf`` finite the two ``q1`` factors differ (the bases
    sit at different potentials) and the identity acquires a
    correction.  Saying so is part of the test; picking a card that
    happens to work and not saying why is what makes a green suite
    uninformative.
    """
    card = dict(NPN_IDEAL, bf=150.0)
    for vid in (0.0, 0.010, 0.050, 0.120):
        c = SubCircuit()
        nc1 = c.add_node('nc1')
        nc2 = c.add_node('nc2')
        ne = c.add_node('ne')
        nb1 = c.add_node('nb1')
        nb2 = c.add_node('nb2')
        nvcc = c.add_node('nvcc')
        nvee = c.add_node('nvee')
        c['vcc'] = VS(nvcc, gnd, v=4.0)
        c['vee'] = VS(nvee, gnd, v=-5.0)
        c['vb1'] = VS(nb1, gnd, v=1.0 + vid / 2)
        c['vb2'] = VS(nb2, gnd, v=1.0 - vid / 2)
        ## A resistive tail rather than an ideal current source, so the
        ## solve has a real fixed point to find.
        c['rt'] = R(ne, nvee, r=3.1e3)
        c['rc1'] = R(nvcc, nc1, r=100.0)
        c['rc2'] = R(nvcc, nc2, r=100.0)
        c['q1'] = eh.GummelPoonNpnHdl(nc1, nb1, ne, **card)
        c['q2'] = eh.GummelPoonNpnHdl(nc2, nb2, ne, **card)
        res = DC(c).solve()
        ## Measured as a voltage drop across a known resistor, which
        ## is how it would be measured on a bench and avoids any
        ## question about the sign convention of a branch current.
        ic1 = (4.0 - float(res.v(nc1, gnd))) / 100.0
        ic2 = (4.0 - float(res.v(nc2, gnd))) / 100.0
        assert ic1 > 1e-9 and ic2 > 1e-9, (vid, ic1, ic2)
        assert 1e-3 < ic1 + ic2 < 4e-3, (vid, ic1 + ic2)
        assert_allclose(ic1 / ic2, math.exp(vid / _UT), rtol=1e-5)
    ## and the pair really did steer: 120 mV is more than 90 to 1.
    assert math.exp(0.120 / _UT) > 90.0


## ----------------------------------------------------------------------
## 1.7  Transient: charge storage against an independent integration

#: A card for the transient.  ``xtf = 0`` and ``xcjc = 1`` so that
#: ``qbe`` is a function of ``vbe`` alone and ``qbc`` of ``vbc`` alone --
#: with the ``xtf`` modulation on, ``tf_eff`` depends on ``vbc`` too and
#: the reference ODE below would have to reproduce the model's own
#: arithmetic instead of being derived from the equations.
TRAN_NPN = dict(IS=1e-16, bf=120.0, vaf=80.0, ikf=2e-2, br=2.0,
                cje=5e-12, vje=0.72, mje=0.35, tf=5e-10,
                cjc=2e-12, vjc=0.6, mjc=0.4, xcjc=1.0, tr=2e-9, fc=0.5)

TR_RB, TR_RC, TR_VCC = 5e3, 2e3, 5.0
TR_FREQ, TR_VO, TR_VA = 2e6, 0.78, 0.03


def _tr_vin(t):
    return TR_VO + TR_VA * math.sin(2 * math.pi * TR_FREQ * t)


def _tr_rhs(vb, vc, t, h=1e-6):
    """``(dvb/dt, dvc/dt)`` of the common-emitter stage, assembled from
    `_gp_reference` and Kirchhoff's law.

    Node b carries ``Qbe(vbe) + Qbc(vbc)`` and node c carries
    ``-Qbc(vbc)``, with ``vbe = vb`` (grounded emitter) and
    ``vbc = vb - vc``, so

        Cbe*vb' + Cbc*(vb' - vc') = (vin - vb)/Rb - Ib
                 -Cbc*(vb' - vc') = (Vcc - vc)/Rc - Ic

    a 2x2 system whose determinant is ``Cbe*Cbc``.  The two capacitances
    are central differences of the reference CHARGE -- not of the
    model's -- so nothing in this function has seen the compiled
    element.
    """
    r = _gp_reference(vb, vb - vc, **TRAN_NPN)
    cbe = (_gp_reference(vb + h, vb + h - vc, **TRAN_NPN)['qbe']
           - _gp_reference(vb - h, vb - h - vc, **TRAN_NPN)['qbe']) / (2 * h)
    cbc = (_gp_reference(vb, vb - vc + h, **TRAN_NPN)['qbc']
           - _gp_reference(vb, vb - vc - h, **TRAN_NPN)['qbc']) / (2 * h)
    a = (_tr_vin(t) - vb) / TR_RB - r['ib']
    b = (TR_VCC - vc) / TR_RC - r['ic']
    ## [[cbe+cbc, -cbc], [-cbc, cbc]] [vb', vc'] = [a, b]
    det = cbe * cbc
    return ((cbc * a + cbc * b) / det,
            (cbc * a + (cbe + cbc) * b) / det)


def _tr_rk4(vb0, vc0, tend, n):
    """Classical RK4 on the two-node system, in plain numpy."""
    dt = tend / n
    vb, vc = vb0, vc0
    out_t, out_b, out_c = [0.0], [vb], [vc]
    for k in range(n):
        t = k * dt
        k1 = _tr_rhs(vb, vc, t)
        k2 = _tr_rhs(vb + 0.5 * dt * k1[0], vc + 0.5 * dt * k1[1],
                     t + 0.5 * dt)
        k3 = _tr_rhs(vb + 0.5 * dt * k2[0], vc + 0.5 * dt * k2[1],
                     t + 0.5 * dt)
        k4 = _tr_rhs(vb + dt * k3[0], vc + dt * k3[1], t + dt)
        vb += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        vc += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        out_t.append(t + dt)
        out_b.append(vb)
        out_c.append(vc)
    return (np.array(out_t), np.array(out_b), np.array(out_c))


def _tr_run(timestep, tend):
    c = SubCircuit()
    nb = c.add_node('nb')
    ncc = c.add_node('ncn')
    nvi = c.add_node('nvi')
    nvcc = c.add_node('nvcc')
    c['vcc'] = VS(nvcc, gnd, v=TR_VCC)
    c['vin'] = eh.VSinHdl(nvi, gnd, vo=TR_VO, va=TR_VA, freq=TR_FREQ)
    c['rb'] = R(nvi, nb, r=TR_RB)
    c['rc'] = R(nvcc, ncc, r=TR_RC)
    c['q1'] = eh.GummelPoonNpnHdl(ncc, nb, gnd, **TRAN_NPN)
    tran = Transient(c, toolkit=numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=tend, timestep=timestep, fixed_timestep=True)
    t = np.asarray(res.v('nb').x[0], float)
    return t, np.asarray(res.v('nb').y, float), \
        np.asarray(res.v('ncn').y, float)


def test_bjt_transient_matches_an_independent_rk4_integration():
    """The stored charge, in time, against a fourth-order integration of
    the same two-node ODE written from the textbook equations.

    A single tolerance can be met by coincidence, so the claim is made
    twice: the waveform agrees with the reference at every sample, AND
    halving the simulator's timestep reduces the disagreement by
    approximately four, which is the second-order convergence of the
    integrator.  A wrong charge model would not converge to the right
    answer at any order.
    """
    tend = 1.5 / TR_FREQ
    errs = []
    for nstep in (300, 600, 1200):
        t, vb, vc = _tr_run(tend / nstep, tend)
        ## The reference starts from the SIMULATOR's own t = 0 point,
        ## which is its DC operating point -- so what is being compared
        ## is the trajectory, not the starting point.
        rt, rb, rc = _tr_rk4(float(vb[0]), float(vc[0]), tend, 20000)
        want_b = np.interp(t, rt, rb)
        want_c = np.interp(t, rt, rc)
        errs.append((np.abs(vb - want_b).max(), np.abs(vc - want_c).max()))
    ## The waveform really does move: the collector swings by more than
    ## 100 mV, so an "agreement" of a few mV is not agreement with a
    ## constant.
    t, vb, vc = _tr_run(tend / 1200, tend)
    assert vc.max() - vc.min() > 0.1, vc.max() - vc.min()
    ## Agreement at the finest step.
    assert errs[-1][0] < 2e-4, errs
    assert errs[-1][1] < 5e-3, errs
    ## and CONVERGENCE: each halving of the step cuts the error by
    ## between 2.5 and 6 (second order is 4; the band is wide because
    ## the error is a max over the waveform and the step controller is
    ## not exactly halving the local truncation everywhere).
    for k in (0, 1):
        ratio = errs[k][1] / errs[k + 1][1]
        assert 2.5 < ratio < 6.0, (k, errs)


def test_charge_storage_is_what_delays_the_collector():
    """The same stage with and without the charge, at a frequency where
    the transit time matters.

    With ``tf = tr = cje = cjc = 0`` the circuit is algebraic and the
    collector follows the base with no phase lag at all; with the card's
    charge it lags.  Measured as the phase of the fundamental, which is
    a property of the waveform rather than of any one sample.
    """
    tend = 3.0 / TR_FREQ
    t, vb, vc = _tr_run(tend / 2400, tend)

    def phase(sig):
        w = 2 * np.pi * TR_FREQ
        s = np.asarray(sig, float) - np.mean(sig)
        return math.atan2(float(np.sum(s * np.cos(w * t))),
                          float(np.sum(s * np.sin(w * t))))

    lag = (phase(vc) - phase(vb)) % (2 * np.pi)
    ## Inverting stage, so the "no charge" answer is exactly pi.
    c2 = SubCircuit()
    nb = c2.add_node('nb')
    ncc = c2.add_node('ncn')
    nvi = c2.add_node('nvi')
    nvcc = c2.add_node('nvcc')
    c2['vcc'] = VS(nvcc, gnd, v=TR_VCC)
    c2['vin'] = eh.VSinHdl(nvi, gnd, vo=TR_VO, va=TR_VA, freq=TR_FREQ)
    c2['rb'] = R(nvi, nb, r=TR_RB)
    c2['rc'] = R(nvcc, ncc, r=TR_RC)
    c2['q1'] = eh.GummelPoonNpnHdl(ncc, nb, gnd,
                                   **dict(TRAN_NPN, tf=0.0, tr=0.0,
                                          cje=0.0, cjc=0.0))
    tran = Transient(c2, toolkit=numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r2 = tran.solve(tend=tend, timestep=tend / 2400,
                        fixed_timestep=True)
    t2 = np.asarray(r2.v('nb').x[0], float)
    vb2 = np.asarray(r2.v('nb').y, float)
    vc2 = np.asarray(r2.v('ncn').y, float)

    def phase2(sig):
        w = 2 * np.pi * TR_FREQ
        s = np.asarray(sig, float) - np.mean(sig)
        return math.atan2(float(np.sum(s * np.cos(w * t2))),
                          float(np.sum(s * np.sin(w * t2))))

    lag0 = (phase2(vc2) - phase2(vb2)) % (2 * np.pi)
    assert abs(lag0 - np.pi) < 0.02, lag0
    ## and with the charge there is a real extra lag, whose SIZE is
    ## predicted: at 2 MHz the base-collector capacitance and the 2 kohm
    ## collector load make a pole at ``1/(2*pi*Rc*Cbc)``, so the extra
    ## phase is ``atan(w*Rc*Cbc)``.  ``Cbc`` comes from the reference
    ## charge at the measured operating point -- 0.96 pF, giving 0.0242
    ## rad -- and the measurement is 0.0227, seven percent under
    ## because Cbc also feeds back into the base.
    h = 1e-6
    vbm, vcm = float(vb.mean()), float(vc.mean())
    cbc = (_gp_reference(vbm, vbm - vcm + h, **TRAN_NPN)['qbc']
           - _gp_reference(vbm, vbm - vcm - h, **TRAN_NPN)['qbc']) / (2 * h)
    predicted = math.atan(2 * np.pi * TR_FREQ * TR_RC * cbc)
    assert 0.5e-12 < cbc < 2e-12, cbc
    extra = lag0 - lag
    assert 0.85 * predicted < extra < 1.05 * predicted, (extra, predicted)


## ----------------------------------------------------------------------
## 1.8  The p-n-p

def test_pnp_is_the_npn_with_every_voltage_and_current_reversed():
    """The p-n-p takes the SAME card -- positive ``IS``, positive
    ``cje``, positive ``vaf`` -- and conducts for ``V(e) > V(b) > V(c)``.

    Its polarity lives entirely in the order its branches are declared
    in, so the claim is exact rather than approximate: every row of
    ``i`` and ``q`` at ``-x`` must be minus the n-p-n's at ``x``, to the
    last bit.  If a sign had been carried in the arithmetic instead
    there would be rounding between the two.
    """
    npn = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **NPN_R)
    pnp = _mk(eh.GummelPoonPnpHdl, 'c', 'b', 'e', **NPN_R)
    assert pnp.terminals == ['c', 'b', 'e']
    ## The internal node order is the same, so -x is the mirrored bias.
    assert [nm for _k, nm, _kd in x_layout(pnp)] == \
        [nm for _k, nm, _kd in x_layout(npn)]
    rng = np.random.default_rng(101)
    for _ in range(120):
        x = np.concatenate([rng.uniform(-3.0, 3.0, 3),
                            rng.uniform(-3.0, 3.0, 3)])
        a = np.asarray(npn.i(x), float)
        b = np.asarray(pnp.i(-x), float)
        assert (a == -b).all(), (x, a, b)
        qa = np.asarray(npn.q(x), float)
        qb = np.asarray(pnp.q(-x), float)
        assert (qa == -qb).all(), (x, qa, qb)
        ## The Jacobians are then EQUAL, not opposite: two sign flips.
        assert (np.asarray(npn.G(x), float)
                == np.asarray(pnp.G(-x), float)).all()
        assert (np.asarray(npn.C(x), float)
                == np.asarray(pnp.C(-x), float)).all()


def test_pnp_limits_the_junction_voltage_that_has_the_exponential_in_it():
    """Why the polarity is in the branch order and not in a sign factor.

    ``limit_pnj`` bounds a branch potential in the SOLUTION vector.  For
    the p-n-p the exponential is in ``V(e) - V(b)``, so the limiter has
    to be declared on that branch; a model that kept the n-p-n's branch
    orientation and multiplied by ``-1`` inside would have `pnjlim`
    compressing the reverse excursions and letting the forward ones --
    the ones that overflow -- straight through.

    The observable: the p-n-p's probes are ``(ei,bi)`` and ``(ci,bi)``,
    with the base as their MINUS terminal, so both probes move their own
    plus and neither is undone.
    """
    pnp = _mk(eh.GummelPoonPnpHdl, 'c', 'b', 'e', **NPN_R)
    spec = pnp._hdl_info['limit_spec']
    assert [(s[0], s[1], s[2]) for s in spec] == \
        [((5, 3), 'pnj', 5), ((4, 3), 'pnj', 4)]
    ## and it does the right thing: a step that would forward-bias the
    ## p-n-p's base-emitter junction to 40 V is compressed.
    x0 = np.array([-2.0, -0.7, 0.0, -0.7, -2.0, 0.0])
    x = np.array([-2.0, -0.7, 40.0, -0.7, -2.0, 40.0])
    out = pnp.limit(x, x0)
    assert out[5] - out[3] < 1.0, out
    ## The n-p-n on the mirrored step is compressed by exactly as much.
    npn = _mk(eh.GummelPoonNpnHdl, 'c', 'b', 'e', **NPN_R)
    outn = npn.limit(-x, -x0)
    ## Both differences are the n-p-n-convention forward bias, so they
    ## are EQUAL and not opposite: the mirror is already in the branch
    ## orientation.
    assert_allclose(out[5] - out[3], outn[3] - outn[5], rtol=1e-12)


## ======================================================================
## 2.  EKV 2.6, the long-channel core
## ======================================================================

#: A card with a real body effect.  ``cox`` is given so the charge model
#: is live; W/L = 10 makes the currents comfortable to read.
EKV = dict(vto=0.5, gamma=0.7, phi=0.7, kp=1.5e-4, cox=6.9e-3,
           w=10e-6, l=1e-6)

#: The IDEAL device: no body effect, so ``n = 1`` EXACTLY and every
#: textbook asymptote is a clean number rather than a slope-factor
#: correction.  Used wherever the reference is a textbook value.
EKV_IDEAL = dict(EKV, gamma=0.0)


def _ekv_vp(vgb, card, T=None):
    """EKV's pinch-off voltage, written out from the 1995 paper.

    ``VP = VG' - PHI - GAMMA*(sqrt(VG' + (GAMMA/2)^2) - GAMMA/2)`` for
    ``VG' > 0`` and ``-PHI`` below, with
    ``VG' = VGB - VTO + PHI + GAMMA*sqrt(PHI)``.  Independent of the
    element: this is the paper's equation, transcribed.
    """
    T = _T0 if T is None else T
    phi = _ekv_phi(card, T)
    vto = card['vto'] - card.get('tcv', 0.0) * (T - card.get('tnom', 300.15))
    g = card['gamma']
    vgp = vgb - vto + phi + g * math.sqrt(phi)
    if vgp > 0.0:
        return vgp - phi - g * (math.sqrt(vgp + 0.25 * g * g) - 0.5 * g)
    return -phi


def _ekv_phi(card, T=None):
    T = _T0 if T is None else T
    tnom = card.get('tnom', 300.15)
    trat = T / tnom
    egt = 1.16 - 7.02e-4 * T ** 2 / (T + 1108.0)
    egn = 1.16 - 7.02e-4 * tnom ** 2 / (tnom + 1108.0)
    ut = _KB * T / _QE
    return card['phi'] * trat - 3.0 * ut * math.log(trat) \
        - egn * trat + egt


def _ekv_n(vp, card, T=None):
    T = _T0 if T is None else T
    ut = _KB * T / _QE
    return 1.0 + card['gamma'] / (2.0 * math.sqrt(vp + _ekv_phi(card, T)
                                                  + 4.0 * ut))


def _ids(el, vd, vg, vs=0.0, vb=0.0):
    return float(np.asarray(el.i(np.array([vd, vg, vs, vb])), float)[0])


def test_ekv_weak_inversion_is_exponential_with_the_slope_factor():
    """Below threshold ``ln(1 + e^x) -> e^x``, so
    ``ID -> Ispec*exp((VP - VS)/UT)`` and the gate characteristic is a
    straight line on a log axis whose slope is ``1/(n*UT)``.

    For the IDEAL device (``gamma = 0``) the slope factor is exactly 1,
    so the subthreshold swing is the textbook ``ln(10)*kT/q`` = 59.5
    mV/decade at 300 K -- a number with no model parameter in it at all.
    With the body effect on, the same measurement returns ``n`` from
    EKV's own formula, evaluated at the pinch-off voltage.
    """
    ## HOW FAR BELOW THRESHOLD MATTERS, and the size of the residual is
    ## known: `softplus(x) = e^x - e^2x/2 + ...`, so the exponential
    ## limit is approached with a relative error of about `e^x/2` in the
    ## SLOPE.  At the biases below, `e^x` is 2e-5 for the ideal device
    ## and 5e-4 for the real one, which is where the two tolerances
    ## come from.  (Measured at 0.10 V, four thermal voltages higher,
    ## the swing is 0.25% off -- a real effect of the interpolation
    ## function, not an error.)
    ideal = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV_IDEAL)
    h = 2e-3

    def swing_at(el, vg):
        i1 = _ids(el, 1.0, vg - h)
        i2 = _ids(el, 1.0, vg + h)
        return 2 * h / math.log10(i2 / i1)

    for vg in (-0.15, -0.10, -0.05):
        assert_allclose(swing_at(ideal, vg), math.log(10.0) * _UT,
                        rtol=1e-4)
    ## The textbook number, at 300 K: 59.5 mV per decade.
    assert_allclose(math.log(10.0) * _UT, 0.059505, rtol=1e-4)
    ## and it IS weak inversion: picoamps, twelve decades below the
    ## strong-inversion current at 3 V.
    assert _ids(ideal, 1.0, -0.05) < 1e-11
    assert _ids(ideal, 1.0, 3.0) > 1e-3

    ## WITH A BODY EFFECT THE SWING IS NOT ``ln(10)*n*UT``, and finding
    ## that out is worth more than the test that assumed it.  The
    ## subthreshold slope is set by ``dVP/dVG``, which differentiating
    ## the pinch-off equation gives in closed form:
    ##
    ##     r  = sqrt(VG' + gamma^2/4)
    ##     dVP/dVG = 1 - gamma/(2r),   so   n_p = 2r/(2r - gamma)
    ##
    ## while EKV's slope factor carries a ``+ 4*UT`` inside its root:
    ## ``n = 1 + gamma/(2*sqrt(VP + PHI + 4*UT))``.  Since
    ## ``sqrt(VP + PHI) = r - gamma/2`` exactly, the two would be
    ## RECIPROCAL without that ``4*UT`` -- with it they differ, by 10%
    ## at these biases.  The regulariser that keeps ``n`` finite at
    ## flat band is visible in the measured swing.
    ##
    ## The remaining 1.8% against ``ln(10)*UT*n_p`` is the specific
    ## current itself: ``Ispec = 2*n*beta*UT^2`` moves with the gate too,
    ## so ``d ln(ID)/dVG = 1/(n_p*UT) + d ln(n)/dVG``.  Including that
    ## term brings the prediction to 1e-4 of the measurement.
    real = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV)
    gam = EKV['gamma']
    for vg in (-0.30, -0.20, -0.10):
        phit = _ekv_phi(EKV)
        vgp = vg - EKV['vto'] + phit + gam * math.sqrt(phit)
        r = math.sqrt(vgp + 0.25 * gam * gam)
        n_p = 2 * r / (2 * r - gam)
        vp = _ekv_vp(vg, EKV)
        n = _ekv_n(vp, EKV)
        dn_dvp = -gam / (4.0 * (vp + phit + 4.0 * _UT) ** 1.5)
        dlnn_dvg = dn_dvp / (n_p * n)
        want = math.log(10.0) / (1.0 / (n_p * _UT) + dlnn_dvg)
        assert_allclose(swing_at(real, vg), want, rtol=1e-3)
        ## The two slope factors, and how far apart they are.
        assert 1.6 < n_p < 2.0, n_p
        assert 0.90 < n / n_p < 0.95, (n, n_p)
        ## and the body effect is not decoration: 98-110 mV/decade
        ## against the ideal device's 59.5.
        assert swing_at(real, vg) > 0.095


def test_ekv_weak_inversion_saturates_at_a_few_thermal_voltages():
    """The other half of the exponential limit:
    ``ID ~ exp(-VS/UT) - exp(-VD/UT)``, so the drain characteristic
    saturates once ``VD - VS`` is a few ``UT`` -- NOT at ``Vdsat``.
    That is the qualitative fact that distinguishes a subthreshold MOS
    from a square-law one, and it is a property of the interpolation
    function's weak-inversion limit alone.
    """
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV_IDEAL)
    isat = _ids(el, 1.0, 0.15)
    for k, tol in ((1.0, 0.40), (4.0, 2e-2), (8.0, 4e-4)):
        i = _ids(el, k * _UT, 0.15)
        assert_allclose(i, isat * (1.0 - math.exp(-k)), rtol=2e-3,
                        err_msg='Vds = %g*UT' % k)
        ## and the number a designer would quote: how close to saturated
        ## the device is at 1, 4 and 8 thermal voltages.
        assert abs(i / isat - 1.0) < tol, (k, i / isat)
    assert _ids(el, 4 * _UT, 0.15) > 0.98 * isat
    assert _ids(el, 1.0 * _UT, 0.15) < 0.65 * isat


def test_ekv_strong_inversion_is_the_square_law():
    """Above threshold ``ln(1 + e^x) -> x``, so

        ID -> (n*beta/2)*[(VP - VS)^2 - (VP - VD)^2]

    which is the textbook square law in triode and
    ``(n*beta/2)*(VP - VS)^2`` in saturation.  The closed form is
    evaluated here from EKV's published pinch-off and slope-factor
    formulas, not from the element.

    The residual is the interpolation function's own tail:
    ``ln(1 + e^x) = x + e^-x + ...``, so ``i`` exceeds ``x^2`` by
    ``2*x*e^-x``.  The worst case in this sweep is the REVERSE term at
    the smallest ``VP - VD`` (``x = 11.6``, giving ``2.1e-4`` on a
    difference of 318), which is 6.6e-7 relative -- hence 2e-6 and not
    "a few percent".  In deep saturation the reverse term is gone and
    the agreement is 1e-9.
    """
    card = EKV_IDEAL
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **card)
    beta = card['kp'] * card['w'] / card['l']
    for vg in (1.6, 2.2, 3.0):
        vp = _ekv_vp(vg, card)
        n = _ekv_n(vp, card)
        assert vp / (2 * _UT) > 20.0, vp
        ## Every `vd` must leave the DRAIN end in strong inversion too,
        ## or the reverse term is in moderate inversion and the square
        ## law is simply not what the model is computing there.  `vd =
        ## vp - 0.3` puts `x_r` at 2.9 and misses by 8e-5, which is the
        ## interpolation being right rather than the model being wrong.
        for vd in (0.1, 0.3, 0.5):
            want = 0.5 * n * beta * ((vp - 0.0) ** 2 - (vp - vd) ** 2)
            assert_allclose(_ids(el, vd, vg), want, rtol=2e-6)
        ## saturation: the reverse current is negligible, so only the
        ## forward square remains.
        want = 0.5 * n * beta * vp ** 2
        assert_allclose(_ids(el, vp + 1.0, vg), want, rtol=1e-9)
    ## and with the body effect on, the same closed form holds with the
    ## bias-dependent n and VP -- which is the claim that `n` is a
    ## SLOPE FACTOR and not a fitted constant.
    el2 = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV)
    for vg in (2.0, 3.0):
        for vsb in (0.0, 1.0):
            vp = _ekv_vp(vg + vsb, EKV)
            n = _ekv_n(vp, EKV)
            want = 0.5 * n * beta * (vp - vsb) ** 2
            assert_allclose(_ids(el2, vp + 1.0 + vsb, vg + vsb, vsb, 0.0),
                            want, rtol=1e-6)


def test_ekv_moderate_inversion_joins_the_two_limits_smoothly():
    """The point of the interpolation function, measured as the
    transconductance efficiency ``gm*n*UT/ID``.

    Differentiating ``ID = Ispec*ln(1 + e^x)^2`` with ``x = (VP-VS)/2UT``
    and ``dVP/dVG = 1/n``, and writing ``sp = ln(1 + e^x) = sqrt(IC)``:

        gm*n*UT/ID = (1 - exp(-sp))/sp

    which is 1 in weak inversion (``sp -> 0``) and ``1/sqrt(IC)`` in
    strong -- the two textbook limits -- and is analytic in between.
    The measured efficiency follows that closed form over SIX decades of
    inversion coefficient, monotonically, with no kink.

    A second, external comparison at the same time: the CHARGE-based
    relation ``i_f = q^2 + q`` gives ``gm*n*UT/ID = 2/(1 + sqrt(1+4*IC))``
    for the same quantity.  The two are DIFFERENT FUNCTIONS -- EKV's
    current interpolation is not the exact charge relation -- and the
    measurement says how different: they agree in both asymptotes (1 in
    weak, ``1/sqrt(IC)`` in strong) and differ by up to **12%** in
    between, because they approach the weak-inversion limit at
    different rates (``1 - sqrt(IC)/2`` against ``1 - IC``).  That is
    the price of the interpolation, it is larger than one would guess,
    and it is why EKV 2.6 uses the charge relation for the CHARGES and
    the interpolation for the CURRENT rather than one for both.
    """
    card = EKV_IDEAL
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **card)
    beta = card['kp'] * card['w'] / card['l']
    prev = None
    worst_charge = 0.0
    for vg in np.linspace(0.05, 1.6, 40):
        vp = _ekv_vp(float(vg), card)
        n = _ekv_n(vp, card)
        ispec = 2.0 * n * beta * _UT ** 2
        ## Deep saturation, so ID = Ispec*i_f and IC = ID/Ispec.
        vd = max(vp, 0.0) + 1.0
        idv = _ids(el, vd, float(vg))
        h = 1e-4
        gm = (_ids(el, vd, float(vg) + h)
              - _ids(el, vd, float(vg) - h)) / (2 * h)
        eff = gm * n * _UT / idv
        ic = idv / ispec
        sp = math.sqrt(ic)
        assert_allclose(eff, (1.0 - math.exp(-sp)) / sp, rtol=3e-4)
        ## monotone decreasing, everywhere -- no kink at any join.
        if prev is not None:
            assert eff < prev, (vg, eff, prev)
        prev = eff
        charge_form = 2.0 / (1.0 + math.sqrt(1.0 + 4.0 * ic))
        worst_charge = max(worst_charge, abs(eff / charge_form - 1.0))
    ## Six decades of inversion coefficient were covered.
    assert prev < 0.2
    ## Measured 12.2%, asserted as a band so that a change in either
    ## formulation shows up rather than being absorbed.
    assert 0.11 < worst_charge < 0.14, worst_charge


def test_ekv_is_symmetric_in_source_and_drain():
    """The model is written on voltages referred to the BULK, so
    exchanging source and drain must negate the current exactly, and
    ``Vds = 0`` must give exactly zero for any gate and body bias.

    The zero is EXACT (the forward and reverse arguments are then the
    same expression) and the antisymmetry is exact to rounding, because
    the drain-referred voltage is reconstructed as ``vds + vsb``.
    """
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV)
    for vg in (0.2, 0.8, 2.0):
        for vsb in (0.0, 0.7, 2.0):
            assert _ids(el, vsb, vg + vsb, vsb, 0.0) == 0.0
            for vds in (0.05, 0.4, 1.5):
                fwd = _ids(el, vsb + vds, vg + vsb, vsb, 0.0)
                rev = _ids(el, vsb, vg + vsb, vsb + vds, 0.0)
                assert_allclose(fwd, -rev, rtol=1e-11)
                assert fwd > 0.0


## ----------------------------------------------------------------------
## 2.2  The four terminal charges

def test_ekv_charges_sum_to_zero_structurally():
    """The gate charge is defined as minus the sum of the other three,
    so neutrality holds to the last bit at every bias -- including
    biases no device sees.  A model that computed all four
    independently would satisfy this only to the accuracy of the
    physics."""
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV)
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(300):
        x = rng.uniform(-20.0, 20.0, 4)
        q = np.asarray(el.q(x), float)
        assert np.isfinite(q).all(), x
        worst = max(worst, abs(q.sum()) / max(np.abs(q).max(), 1e-30))
    ## Not exactly zero -- it is a sum of four doubles -- but at the
    ## level of accumulated rounding rather than of physics.
    assert worst < 1e-13, worst
    ## and the same for the capacitance matrix: every row and column
    ## sums to zero, which is the statement that the device has no
    ## charge of its own.
    for x in ([1.0, 2.0, 0.0, 0.0], [0.2, 0.4, 0.1, -0.5],
              [3.0, 3.0, 0.0, -2.0]):
        C = np.asarray(el.C(np.array(x, float)), float)
        assert np.abs(C.sum(axis=0)).max() < 1e-12 * np.abs(C).max()
        assert np.abs(C.sum(axis=1)).max() < 1e-12 * np.abs(C).max()


def test_ekv_ward_dutton_partition_is_50_50_and_40_60():
    """The charge partition every quasi-static MOS model must reproduce:
    at ``Vds = 0`` the drain and source each hold half the inversion
    charge, and in saturation the split is 40/60.

    Those two numbers are properties of ``Qs = W*int (1 - x/L) Q'`` and
    ``Qd = W*int (x/L) Q'`` over a channel, not of EKV; they are what a
    charge model is checked against.  The 40/60 is asymptotic -- the
    ``- 1/2`` in the normalised charges pulls it to 0.398 at this drive
    -- so the tolerance says which.
    """
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV_IDEAL)
    ## Rows 0 and 2 of q are the drain and source charges.
    q = np.asarray(el.q(np.array([0.0, 3.0, 0.0, 0.0])), float)
    assert q[0] < 0 and q[2] < 0
    assert_allclose(q[0] / (q[0] + q[2]), 0.5, rtol=1e-12)
    ## Saturation: VP is 2.5 V, so 4 V on the drain is well past it.
    q = np.asarray(el.q(np.array([4.0, 3.0, 0.0, 0.0])), float)
    assert_allclose(q[0] / (q[0] + q[2]), 0.4, rtol=6e-3)
    ## and it really is asymptotic: at a higher drive it gets closer.
    q2 = np.asarray(el.q(np.array([9.0, 8.0, 0.0, 0.0])), float)
    e1 = abs(q[0] / (q[0] + q[2]) - 0.4)
    e2 = abs(q2[0] / (q2[0] + q2[2]) - 0.4)
    assert e2 < 0.5 * e1, (e1, e2)


def test_ekv_gate_capacitance_hits_the_textbook_limits():
    """``Cgg`` against the two numbers every MOS text prints, on the
    IDEAL device where the slope factor is exactly 1 and they are exact:

    * at ``Vds = 0`` in strong inversion the channel is uniform and the
      whole oxide capacitance is seen from the gate, ``W*L*Cox``, split
      equally between source and drain;
    * in saturation the channel is pinched off at the drain and the
      gate sees ``(2/3)*W*L*Cox``.

    Both are limits, so the test also checks that raising the gate drive
    moves the measurement TOWARDS them.
    """
    card = EKV_IDEAL
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **card)
    wlcox = card['w'] * card['l'] * card['cox']
    ## Row/col 1 is the gate.
    def cgg(vd, vg):
        return float(np.asarray(el.C(np.array([vd, vg, 0.0, 0.0])),
                                float)[1, 1])
    ## Vds = 0.
    assert_allclose(cgg(0.0, 5.0), wlcox, rtol=2e-2)
    assert_allclose(cgg(0.0, 12.0), wlcox, rtol=8e-3)
    assert abs(cgg(0.0, 12.0) - wlcox) < abs(cgg(0.0, 5.0) - wlcox)
    ## and the 50/50 split of that between source and drain.
    C = np.asarray(el.C(np.array([0.0, 12.0, 0.0, 0.0])), float)
    assert_allclose(-C[1, 0], 0.5 * wlcox, rtol=1e-2)
    assert_allclose(-C[1, 2], 0.5 * wlcox, rtol=1e-2)
    ## Saturation.
    assert_allclose(cgg(6.0, 5.0), 2.0 / 3.0 * wlcox, rtol=2e-2)
    assert_allclose(cgg(14.0, 12.0), 2.0 / 3.0 * wlcox, rtol=8e-3)
    ## The body effect changes the answer, and by a predictable amount:
    ## with `n > 1` the gate also charges the depletion region, so Cgg
    ## in saturation is `(2/(3n) + (n-1)/n)*W*L*Cox`.
    el2 = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV)
    vp = _ekv_vp(12.0, EKV)
    n = _ekv_n(vp, EKV)
    c2 = float(np.asarray(el2.C(np.array([14.0, 12.0, 0.0, 0.0])),
                          float)[1, 1])
    assert_allclose(c2, (2.0 / (3.0 * n) + (n - 1.0) / n) * wlcox,
                    rtol=2e-2)
    assert c2 > 2.0 / 3.0 * wlcox


def test_ekv_cox_zero_switches_the_charge_off_completely():
    """The charge is a PRODUCT with ``cox``, never a quotient, so the
    default card has exactly zero charge rather than a zero multiplying
    an infinity -- the trap the roadmap names and the thermal node was
    designed around."""
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **dict(EKV, cox=0.0))
    rng = np.random.default_rng(9)
    for _ in range(60):
        x = rng.uniform(-30.0, 30.0, 4)
        assert (np.asarray(el.q(x), float) == 0.0).all(), x
        assert (np.asarray(el.C(x), float) == 0.0).all(), x
    ## and the current is untouched by it.
    full = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV)
    x = np.array([1.5, 1.2, 0.0, 0.0])
    assert_allclose(np.asarray(el.i(x), float),
                    np.asarray(full.i(x), float), rtol=1e-14)


## ----------------------------------------------------------------------
## 2.3  Noise

def test_ekv_channel_noise_hits_the_nyquist_and_two_thirds_limits():
    """Channel thermal noise from the Klaassen-Prins integral,
    ``S = 4*k*T*mu*|Q_I|/L^2``, checked against the three limits it must
    reproduce -- and against the size of its departure from two of them.

    * **weak inversion, saturated**: EXACTLY ``2*q*I``, full shot noise.
      That falls out with no approximation: the normalised charges
      reduce to ``qs + qd = i_f`` there, so ``4kT*gn = 4kT*I/(2*UT)
      = 2*q*I``.  It is the textbook result for a subthreshold MOS and
      it is asserted to 1e-3;
    * **Vds = 0**: Nyquist, ``4*k*T*gds`` -- but only ASYMPTOTICALLY,
      and the gap is measured rather than tolerated.  ``gn`` is built
      from the CHARGE relation ``q = sqrt(1/4 + i) - 1/2`` while ``gds``
      comes from the CURRENT interpolation ``i = ln(1+e^x)^2``, and EKV
      2.6 uses both; the ratio is
      ``(2*sqrt(1/4+i) - 1)/(2*sqrt(i)*sigmoid(x))``, which is
      ``1 - 1/(2*sqrt(IC))`` in strong inversion.  Measured 0.932 at
      1 V of gate drive and 0.987 at 3 V -- the same
      interpolation-versus-charge inconsistency the transconductance
      test measures at 12%, seen from another side;
    * **saturation**: ``(2/3)*4*k*T*n*gm``, approached from above --
      0.668 at 2 V, 0.680 at 4 V, against 0.6667.

    Nothing here is tuned to make the limits exact.  The departures are
    a property of EKV 2.6's mixed formulation and are worth a number.
    """
    card = EKV
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **card)
    kt4 = 4.0 * _KB * _T0

    ## -- weak inversion, saturated: exactly 2*q*I.
    for vg in (-0.2, -0.1, 0.0):
        x = np.array([1.0, vg, 0.0, 0.0])
        idv = _ids(el, 1.0, vg)
        psd = float(np.asarray(el.CY(x, 2 * np.pi * 1e3), float)[0, 0])
        assert idv < 1e-8
        assert_allclose(psd, 2.0 * _QE * idv, rtol=1e-3)

    ## -- Vds = 0: Nyquist, and the exact size of the departure.
    ratios = []
    for vg in (1.0, 2.0, 3.0):
        x = np.array([0.0, vg, 0.0, 0.0])
        gds = float(np.asarray(el.G(x), float)[0, 0])
        psd = float(np.asarray(el.CY(x, 2 * np.pi * 1e3), float)[0, 0])
        assert gds > 0
        ratio = psd / (kt4 * gds)
        ratios.append(ratio)
        ## The closed form, from the two EKV relations:
        vp = _ekv_vp(vg, card)
        xx = vp / (2.0 * _UT)
        sp = math.log1p(math.exp(xx)) if xx < 30 else xx
        want = (2.0 * math.sqrt(0.25 + sp * sp) - 1.0) \
            / (2.0 * sp / (1.0 + math.exp(-xx)))
        assert_allclose(ratio, want, rtol=1e-9)
        ## and the strong-inversion asymptote of THAT: 1 - 1/(2*sqrt(IC)),
        ## which is good to 0.4% at 1 V of drive and 0.02% at 3 V --
        ## the next term is 1/(8*IC).
        assert_allclose(ratio, 1.0 - 1.0 / (2.0 * sp), rtol=5e-3)
    ## It approaches Nyquist from below, monotonically.
    assert ratios[0] < ratios[1] < ratios[2] < 1.0
    assert ratios[0] > 0.92 and ratios[2] > 0.98

    ## -- saturation: the 2/3 excess-noise factor.
    factors = []
    for vg in (2.0, 3.0, 4.0):
        vp = _ekv_vp(vg, card)
        n = _ekv_n(vp, card)
        x = np.array([vp + 1.5, vg, 0.0, 0.0])
        gm = float(np.asarray(el.G(x), float)[0, 1])
        psd = float(np.asarray(el.CY(x, 2 * np.pi * 1e3), float)[0, 0])
        assert gm > 0
        factors.append(psd / (kt4 * n * gm))
    assert all(abs(f - 2.0 / 3.0) < 0.025 for f in factors), factors
    assert 0.66 < min(factors) and max(factors) < 0.69


def test_ekv_flicker_noise_is_one_over_f():
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **dict(EKV, kf=1e-12,
                                                       af=1.1, ef=1.0))
    x = np.array([2.0, 1.5, 0.0, 0.0])
    idv = _ids(el, 2.0, 1.5)
    thermal = float(np.asarray(
        _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV).CY(
            x, 2 * np.pi * 1e3), float)[0, 0])
    for f in (1e2, 1e3, 1e4):
        psd = float(np.asarray(el.CY(x, 2 * np.pi * f), float)[0, 0])
        assert_allclose(psd - thermal, 1e-12 * abs(idv) ** 1.1 / f,
                        rtol=1e-9)
    ## and it dominates at low frequency, which is what makes it worth
    ## having: 100 Hz is three decades above the thermal floor here.
    low = float(np.asarray(el.CY(x, 2 * np.pi * 1e2), float)[0, 0])
    assert low / thermal > 100.0


## ----------------------------------------------------------------------
## 2.4  Jacobians and finiteness

@pytest.mark.parametrize('name,card,x,resolved', [
    ('weak inversion', EKV, [1.0, 0.1, 0.0, 0.0], True),
    ('moderate', EKV, [1.0, 0.55, 0.0, 0.0], True),
    ('strong saturation', EKV, [2.0, 2.0, 0.0, 0.0], True),
    ('triode', EKV, [0.05, 2.0, 0.0, 0.0], True),
    ('vds zero', EKV, [1.0, 2.0, 1.0, 0.0], True),
    ('body bias', EKV, [3.0, 2.5, 1.0, -1.0], True),
    ## Right at the `vgp = 0` seam of the pinch-off Piecewise.
    ('pinch-off seam', EKV, [1.0, -0.2856, 0.0, 0.0], True),
    ## DEEP CUTOFF IS NOT RESOLVABLE, and this is the point that bought
    ## `check_jacobians` its third verdict.  Every entry of `C` here is
    ## around 1e-25 F, so the one-band tolerance is 1e-32, and the
    ## difference at the default 1e-7 step comes back EXACTLY 0.0: the
    ## normalised charges are computed as a cancellation, so `q`'s
    ## representable step is `eps` times an INTERNAL magnitude of
    ## 1.9e-15, not times `|q|` itself.
    ##
    ## This used to be pinned with a hand-written `atol = 1e-24` and a
    ## comment proposing `max(atol, eps*|q|/h)` as the general fix.  That
    ## formula is TEN DECADES too small here (9.4e-36 against a real
    ## quantum of 4.1e-31), and even `eps*max|q|/h` misses by 2x, so the
    ## floor is now MEASURED: `check_jacobians` widens the step until the
    ## value moves clear of its own quantisation.  Doing that turns a
    ## difference of 0.0 into -1.26e-25 against an analytic -1.28e-25 --
    ## 1.5% agreement where there was none at all -- and the entries are
    ## reported UNRESOLVED rather than either FAILED or silently passed.
    ('accumulation', EKV, [1.0, -3.0, 0.0, 0.0], False),
    ('ideal device', EKV_IDEAL, [2.0, 2.0, 0.0, 0.0], True),
    ('ideal at the seam', EKV_IDEAL, [1.0, -0.2, 0.0, 0.0], False),
])
def test_ekv_jacobians_by_finite_differences(name, card, x, resolved):
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **card)
    res = check_jacobians(el, np.array(x, dtype=float), rtol=3e-5)
    assert res.ok, '%s\n%s' % (name, res)
    ## `ok` alone would be satisfied by an instrument that gave up
    ## everywhere, so the two states are pinned separately: every
    ## ordinary bias must be fully RESOLVED, and the two cutoff points
    ## must report themselves unresolved rather than pass quietly.
    assert res.resolved is resolved, '%s\n%s' % (name, res)
    if not resolved:
        assert all(u.reason == 'roundoff' for u in res.unresolved), \
            '%s\n%s' % (name, res)


def test_ekv_stays_finite_where_no_device_belongs():
    """Both arms of the pinch-off conditional are evaluated at every
    bias, and ``gamma = 0`` is a legitimate card that makes the
    discarded arm's square root sit exactly at zero -- so a
    ``0 * d(sqrt(0))`` would be NaN and would only show up on that
    card.  Both cards are swept."""
    for card in (EKV, EKV_IDEAL):
        el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **card)
        rng = np.random.default_rng(20260826)
        for _ in range(300):
            x = rng.uniform(-200.0, 200.0, 4)
            for f in (el.i, el.q, el.G, el.C):
                v = np.asarray(f(x), float)
                assert np.isfinite(v).all(), (card['gamma'], x, v)
        ## and exactly at the seam, where `vgp = 0`.
        phit = _ekv_phi(card)
        vg0 = card['vto'] - phit - card['gamma'] * math.sqrt(phit)
        for d in (-1e-15, 0.0, 1e-15, -1e-3, 1e-3):
            x = np.array([1.0, vg0 + d, 0.0, 0.0])
            for f in (el.i, el.q, el.G, el.C):
                assert np.isfinite(np.asarray(f(x), float)).all(), (card, d)


## ----------------------------------------------------------------------
## 2.5  limit_fet and limit_vds, on their first production model

def test_ekv_declares_fetlim_on_the_gate_and_limvds_on_the_drain():
    """The declaration, and the two things about it a reader needs.

    The probes are ``V(g,s)`` and ``V(d,s)`` -- SPICE's pair.  Both name
    the source, so the two write-backs must land on different rows; each
    probe then carries exactly its own limited value and the declaration
    order does not matter.

    WHICH row each takes is a runtime choice: the terminal that drifted
    further from the last accepted point.  Fixing it at compile time --
    always the gate, always the drain -- is what let a wild source drag
    a sane drain out with it, a decade per iteration.  So this asserts
    that each branch ends at its own limited value and that the body is
    never touched, not a fixed pair of indices.
    """
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV)
    assert el.terminals == ['d', 'g', 's', 'b']
    spec = el._hdl_info['limit_spec']
    assert [(s[0], s[1], s[2]) for s in spec] == \
        [((1, 2), 'fet', 1), ((0, 2), 'vds', 0)]
    assert ('2 $limit probes (fetlim on (g,s) [params at last iterate], '
            'limvds on (d,s))') in explain(el), explain(el)
    rng = np.random.default_rng(17)
    moved_source = 0
    for _ in range(300):
        x = rng.uniform(-30.0, 30.0, 4)
        x0 = rng.uniform(-30.0, 30.0, 4)
        out = el.limit(x, x0)
        ## The bulk is not a probe terminal and is never written.
        assert out[3] == x[3]
        movers = [i for i in range(4) if out[i] != x[i]]
        assert len(set(movers)) == len(movers)
        assert len(movers) <= 2 and set(movers) <= {0, 1, 2}, movers
        for i in range(4):
            if i not in movers:
                assert out[i] == x[i]
        ## Each branch ends at a bounded value -- `_limvds` clamps vds
        ## and `_fetlim` bounds the gate step, so neither branch may come
        ## out wilder than it went in.
        assert abs(out[1] - out[2]) <= max(abs(x[1] - x[2]), 60.0)
        assert abs(out[0] - out[2]) <= max(abs(x[0] - x[2]), 60.0)
        if out[2] != x[2]:
            moved_source += 1
        assert not np.shares_memory(out, x)
    ## And the source does give way sometimes -- under the old fixed
    ## rule it never did, which was the defect.
    assert moved_source > 10, moved_source


def test_fetlim_is_placed_at_the_body_biased_turn_on():
    """`limit_fet` now sees SPICE's ``von``, and the 565 mV is gone.

    SPICE passes ``fetlim`` a turn-on voltage recomputed each iteration
    from the previous iterate's bulk bias.  Until 2026-08-25 this model
    could only pass its zero-bias ``vto``, because a limiter parameter
    that reached the solution was refused at compile time; at 2 V of
    body bias that put every clamp **565 mV** below the real turn-on of
    1.06 V.  (This test used to pin that looseness as a documented
    limitation.  The limitation was a rule about ORDER -- "the limiter
    runs before the device" -- whose conclusion did not follow: the
    limiter is handed the last accepted iterate precisely so it can
    measure against it.)

    The model now writes ``von = vtoT + gamma*(sqrt(phi + vsb) -
    sqrt(phi))`` and the parameter is evaluated at ``x0``.  So an off
    device asked to jump to 100 V lands at ``von + 0.5``, where SPICE
    puts it -- and the no-op band for a small step is untouched, so
    "did limiting fire?" stays a usable convergence signal.
    """
    card = EKV
    vsb = 2.0
    phit = _ekv_phi(card)
    von = card['vto'] + card['gamma'] * (math.sqrt(phit + vsb)
                                         - math.sqrt(phit))
    assert 1.05 < von < 1.08, von
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **card)
    fet = [s for s in el._hdl_info['limit_spec'] if s[1] == 'fet'][0]
    assert fet[3][0]._wants_x is True, 'von must read the bulk bias'
    rng = np.random.default_rng(23)
    for _ in range(200):
        x0 = np.array([vsb + 1.0, vsb + rng.uniform(0.0, 1.5), vsb, 0.0])
        x = np.array([vsb + 1.0, vsb + rng.uniform(-50.0, 50.0), vsb, 0.0])
        out = el.limit(x, x0)
        assert abs(out[1] - out[2]) < 60.0
        small = np.array([x0[0], x0[1] + 0.4, vsb, 0.0])
        assert (el.limit(small, x0) == small).all()
    ## The measurement that used to record the looseness now records
    ## its absence: the 100 V jump from an off device lands at von + 0.5.
    x0 = np.array([vsb + 1.0, vsb + 0.0, vsb, 0.0])
    out = el.limit(np.array([vsb + 1.0, vsb + 100.0, vsb, 0.0]), x0)
    assert_allclose(out[1] - out[2], von + 0.5, rtol=1e-9)
    ## and at ZERO body bias von == vto (to the temperature scaling), so
    ## the isothermal, unbiased behaviour is exactly what it was.
    x0 = np.array([1.0, 0.0, 0.0, 0.0])
    out = el.limit(np.array([1.0, 100.0, 0.0, 0.0]), x0)
    assert_allclose(out[1] - out[2], card['vto'] + 0.5, rtol=1e-9)


#: The control element for the limiter experiments: `EkvNmosHdl`'s body
#: with the two ``$limit`` declarations removed and NOTHING else changed.
#: Built lazily, because compiling it costs a couple of seconds and only
#: three tests want it.
_NOLIMIT_CACHE = []


def _ekv_no_limit():
    if not _NOLIMIT_CACHE:
        from pycircuit.circuit.hdl import Behavioural

        class _EkvNoLimit(Behavioural):
            instparams = eh._ekv_params()
            analog = staticmethod(_ekv_analog_unlimited())

        _NOLIMIT_CACHE.append(_EkvNoLimit)
    return _NOLIMIT_CACHE[0]


def _count_jacobians(cir, x0=None, gmin=None):
    """Solve a circuit's DC point and count Jacobian evaluations.

    One per Newton iteration, which is the currency a limiter is judged
    in.  Returns ``(count, result)`` or ``(count, exception)`` when the
    solve failed.

    ``gmin`` reaches `DC`'s anchor (roadmap §12.3).  Pass ``0.0`` to
    measure a limiter against a bare solve: the anchor rescues an empty
    row on its own, so with it left at the default a test comparing
    "with limiting" to "without" is measuring the anchor as well.
    """
    calls = []
    orig = cir.G

    def counted(*a, **kw):
        calls.append(1)
        return orig(*a, **kw)

    cir.G = counted
    kw = {} if gmin is None else dict(gmin=gmin)
    try:
        res = DC(cir, **kw).solve(x0=x0)
    except Exception as exc:                     # noqa: BLE001
        return len(calls), exc
    return len(calls), res


LIM_CARD = dict(vto=0.5, gamma=0.7, phi=0.7, kp=1.5e-4, cox=6.9e-3,
                w=20e-6, l=2e-6)


def _forced_current_stage(cls):
    """A diode-connected EKV with a milliamp forced into it from a
    current source.  Nothing pins the drain, so a Newton step that
    overshoots takes the device into a region where the channel
    conductance underflows to EXACTLY zero and the row goes empty."""
    from pycircuit.circuit.elements import IS as ISource
    c = SubCircuit()
    nd = c.add_node('nd')
    ng = c.add_node('ng')
    c['vg'] = VS(ng, gnd, v=3.0)
    c['isrc'] = ISource(gnd, nd, i=1e-3)
    c['m'] = cls(nd, ng, gnd, gnd, **LIM_CARD)
    return c, nd


def test_fet_limiting_rescues_a_solve_that_otherwise_goes_singular():
    """What ``limit_fet``/``limit_vds`` buy, on a real solve.

    A milliamp forced into a diode-connected EKV device from a current
    source, started from the origin.  With the limiters the DC analysis
    converges in a handful of Jacobian evaluations; with the SAME model
    and only the two ``$limit`` declarations removed, a Newton step
    drives ``VP - VD`` so far negative that ``ln(1 + e^x)`` underflows
    to exactly 0.0 and the drain row becomes empty.

    **Updated 2026-08-25.** That empty row used to end the solve with a
    `SingularMatrix`; roadmap §12.3 added a gmin anchor that now rescues
    it, so the assertion is made at ``gmin = 0`` -- the state this test
    was written to measure -- and the rescue is asserted separately.
    Both halves matter: without the anchor the limiters are what save
    the solve, and with it the two routes must agree, because a
    limiter buys a shorter PATH and an anchor buys a well-posed step,
    and neither is allowed to move the answer.
    """
    n_lim, res = _count_jacobians(*[_forced_current_stage(eh.EkvNmosHdl)[0]])
    assert not isinstance(res, Exception), res
    assert n_lim < 15, n_lim

    ## Without the limiters AND without the anchor: still singular.
    c2, _nd = _forced_current_stage(_ekv_no_limit())
    _n2, res2 = _count_jacobians(c2, gmin=0.0)
    from pycircuit.circuit.analysis import SingularMatrix
    assert isinstance(res2, SingularMatrix), res2

    ## Without the limiters but WITH the anchor: rescued, and it lands
    ## on the limited model's own answer.  This is the assertion that
    ## would catch an anchor which quietly moved the solution.
    c3, nd3 = _forced_current_stage(_ekv_no_limit())
    _n3, res3 = _count_jacobians(c3)
    assert not isinstance(res3, Exception), res3
    assert_allclose(float(res3.v(nd3, gnd)),
                    float(res.v(_forced_current_stage(eh.EkvNmosHdl)[1], gnd)),
                    rtol=1e-6)
    ## and the point it converged to is a real one: 0.3 V of drain-source
    ## voltage carrying a milliamp, i.e. deep triode.
    c3, nd3 = _forced_current_stage(eh.EkvNmosHdl)
    r3 = DC(c3).solve()
    assert 0.2 < float(r3.v(nd3, gnd)) < 0.5, float(r3.v(nd3, gnd))


def _cascode_stage(cls, rleak=1e9):
    """Two stacked EKV devices from a 40 V rail.

    ``rleak`` is not decoration and its absence is a finding: with
    nothing else touching the intermediate node, a Newton iterate that
    puts the lower device in deep cutoff makes its channel conductance
    underflow to EXACTLY zero, and the node is then determined by
    nothing at all -- ``SingularMatrix: 'nm' appears in no equation``.
    EKV as published has no ``GMIN``; `compact.PspMosLongChannel` carries
    a ``GLEAK`` for exactly this reason and this model deliberately does
    not, because a 1 pS anchor is a picoamp at a volt and would sit on
    top of the weak-inversion currents this model exists to measure.
    A gigaohm to ground changes the answer by 1e-4 and removes the
    hazard; that is the testbench's job, not the model's.
    """
    c = SubCircuit()
    nd = c.add_node('nd')
    nm = c.add_node('nm')
    ng = c.add_node('ng')
    nvdd = c.add_node('nvdd')
    c['vdd'] = VS(nvdd, gnd, v=40.0)
    c['vg'] = VS(ng, gnd, v=1.0)
    c['rd'] = R(nvdd, nd, r=1e5)
    c['rl'] = R(nm, gnd, r=rleak)
    c['m1'] = cls(nd, ng, nm, gnd, **LIM_CARD)
    c['m2'] = cls(nm, ng, gnd, gnd, **LIM_CARD)
    return c, nd


def test_a_stacked_pair_without_a_dc_path_is_singular_at_the_start_not_at_the_answer():
    """The hazard `_cascode_stage` documents, pinned so that it is a
    known property rather than a surprise.

    Limiting does not REMOVE it -- a per-probe limiter bounds a step,
    and this is a row that has gone empty -- but it does change which
    circuits reach it, and that is worth recording rather than
    discovering later.  At (vdd, vg) = (5, 2.5) this used to raise, and
    now converges: the runtime write-back keeps the iterate out of the
    cutoff region where the lower device's channel conductance
    underflows to EXACTLY zero.  So the bias below is a harder one,
    chosen because it still reaches the empty row.
    """
    from pycircuit.circuit.analysis import SingularMatrix

    def _pair(vdd, vg):
        ## Built by hand, because `_cascode_stage` always adds the leak.
        cc = SubCircuit()
        nd = cc.add_node('nd')
        nm = cc.add_node('nm')
        ng = cc.add_node('ng')
        nvdd = cc.add_node('nvdd')
        cc['vdd'] = VS(nvdd, gnd, v=vdd)
        cc['vg'] = VS(ng, gnd, v=vg)
        cc['rd'] = R(nvdd, nd, r=2e3)
        cc['m1'] = eh.EkvNmosHdl(nd, ng, nm, gnd, **LIM_CARD)
        cc['m2'] = eh.EkvNmosHdl(nm, ng, gnd, gnd, **LIM_CARD)
        return cc

    ## **The name of this test used to say "structurally singular", and
    ## that was wrong.**  Measured while building the gmin anchor
    ## (roadmap §12.3): this circuit is singular at ITERATE 0 and
    ## perfectly well posed at its answer, where the `nm` conductance is
    ## 3.8e-8 S -- four decades above the anchor that gets it there.  An
    ## empty row here is a numerical fact about a starting point, not a
    ## structural fact about the netlist, and the two need different
    ## remedies.  With the anchor it now converges.
    r_anchored = DC(_pair(40.0, 0.2)).solve()
    assert np.isfinite(float(r_anchored.v('nm', gnd)))

    ## At gmin = 0 the empty row is still reachable, and the message
    ## still names it.  This is the state the test was written to pin.
    with pytest.raises(SingularMatrix, match="'nm'"):
        DC(_pair(40.0, 0.2), gmin=0.0).solve()

    ## The measurement of what the write-back fix bought, taken at
    ## gmin = 0 so that it is the limiter being measured and not the
    ## anchor: the bias this test used to use now solves with neither a
    ## leak resistor nor an anchor.
    r_easy = DC(_pair(5.0, 2.5), gmin=0.0).solve()
    assert 0.4 < float(r_easy.v('nm', gnd)) < 0.5

    ## and the hard one with a gigaohm to ground converges.
    cc = _pair(40.0, 0.2)
    cc['rl'] = R(cc.get_node('nm'), gnd, r=1e9)
    r5 = DC(cc).solve()
    assert np.isfinite(float(r5.v('nm', gnd)))
    ## The leak is a testbench artefact and not a fitting knob: three
    ## decades of it move the 40 V stage's answer by 2e-3.
    vals = []
    for rleak in (1e6, 1e9):
        cc2, _nd2 = _cascode_stage(eh.EkvNmosHdl, rleak=rleak)
        vals.append(float(DC(cc2).solve().v('nm', gnd)))
    assert_allclose(vals[0], vals[1], rtol=2e-3)
    assert 0.10 < vals[1] < 0.12, vals


def test_fet_limiting_cuts_the_iteration_count_on_a_hard_solve():
    """From the origin, the 40 V cascode takes 9 Jacobians with the
    limiters and 25 without -- a factor of 2.8, on a solve that both
    reach and reach to the same operating point.

    Both numbers are asserted as bands rather than exactly, because
    they are a property of the solver as well as of the limiter.
    """
    n_lim, res_lim = _count_jacobians(_cascode_stage(eh.EkvNmosHdl)[0])
    n_raw, res_raw = _count_jacobians(_cascode_stage(_ekv_no_limit())[0])
    assert not isinstance(res_lim, Exception), res_lim
    assert not isinstance(res_raw, Exception), res_raw
    assert_allclose(float(res_lim.v('nd', gnd)),
                    float(res_raw.v('nd', gnd)), rtol=1e-6)
    assert n_lim <= 14, n_lim
    assert n_raw >= 18, n_raw
    assert n_raw / n_lim > 1.8, (n_lim, n_raw)


def test_the_fet_write_back_moves_the_terminal_that_actually_drifted():
    """**Was a finding; is now a fix, with the old numbers kept.**

    The state-free ``$limit`` convention writes a probe's limited value
    back by moving one terminal, and that terminal used to be chosen at
    COMPILE time -- always the plus, so always the gate for
    ``limit_fet(V(g,s))`` and always the drain for ``limit_vds(V(d,s))``.

    That made the limiter a divergence GENERATOR.  The write-back is
    ``x[ra] = x[rb] + vlim``: ``vlim`` is bounded, ``x[rb]`` is not, so
    a wild node hands its magnitude to a sane one.  Elements are limited
    in instance order, so the upper device read a source node the lower
    one had not fixed yet:

        it0  nm  +5.17e+07 -> +3.20e+01   (lower device fixes its node)
        it1  nd  +4.00e+01 -> +6.14e+08   (upper device then destroys a
        it2  nd  +4.00e+01 -> +6.14e+09    perfectly good drain, a
        it3  nd  +4.00e+01 -> +5.78e+10    decade per iteration)

    Newton had the drain at a sane 40 V.  Measured then and now, same
    circuit, same start:

        limited (compile-time write-back)   225 Jacobian evaluations
        limited (runtime write-back)          8
        unlimited                            25

    and the worst single displacement of the gate -- which ``vg`` pins
    at 1.0 V -- went from **5e48 V** to a few tens of volts.

    The fix is to move whichever terminal has drifted further from the
    last accepted point: that is the node Newton is being wild about,
    and the other one is information worth keeping.
    """
    c_lim, _ = _cascode_stage(eh.EkvNmosHdl)
    x0 = np.full(c_lim.n, 10.0)
    x0[-1] = 0.0
    moves = {'n': 0, 'gate': 0, 'worst': 0.0}
    orig_limit = c_lim.limit

    def watched(x, xold, epar=defaultepar):
        before = np.array(x, dtype=float)
        out = orig_limit(x, xold, epar)
        moves['n'] += 1
        d = abs(float(out[2]) - float(before[2]))    # index 2 is 'ng'
        if d > 1e-12:
            moves['gate'] += 1
            moves['worst'] = max(moves['worst'], d)
        return out

    c_lim.limit = watched
    n_lim, res_lim = _count_jacobians(c_lim, x0=x0)
    c_raw, _ = _cascode_stage(_ekv_no_limit())
    n_raw, res_raw = _count_jacobians(c_raw, x0=np.array(x0))
    ## Both find the same answer -- limiting never moves the solution,
    ## only the path to it.
    assert not isinstance(res_lim, Exception), res_lim
    assert not isinstance(res_raw, Exception), res_raw
    assert_allclose(float(res_lim.v('nd', gnd)),
                    float(res_raw.v('nd', gnd)), rtol=1e-6)
    ## Limiting now HELPS on the case where it used to hurt 9x.
    assert n_lim < n_raw, (n_lim, n_raw)
    ## The gate may still give way -- it is a legitimate choice when the
    ## gate is what moved -- but never by a wild amount.  5e48 was the
    ## number that said the write-back was broken.
    assert moves['worst'] < 1e3, moves


def _ekv_analog_unlimited():
    """`EkvNmosHdl`'s body with the two ``$limit`` declarations removed
    and nothing else changed -- the control for the test above."""
    import sympy
    from pycircuit.circuit.hdl import (Branch, Contribution, ddt, var, vt,
                                       maxc, softplus, safe_pow, hypsmooth,
                                       TEMP)

    def analog(d, g, s, b):
        bgs, bds, bsb = Branch(g, s), Branch(d, s), Branch(s, b)
        bdb, bgb = Branch(d, b), Branch(g, b)
        T = TEMP
        ut = var(vt(T), 'ut')
        trat = var(T / tnom, 'trat')                                # noqa
        ltr = var(sympy.log(trat), 'ltrat')
        egT = var(1.16 - 7.02e-4 * T ** 2 / (T + 1108.0), 'egT')
        egn = var(1.16 - 7.02e-4 * tnom ** 2 / (tnom + 1108.0),     # noqa
                  'egtnom')
        vtoT = var(vto - tcv * (T - tnom), 'vtoT')                  # noqa
        kpT = var(kp * safe_pow(trat, bex, lo=1e-3), 'kpT')         # noqa
        phiT = var(phi * trat - 3.0 * ut * ltr - egn * trat + egT,  # noqa
                   'phiT')
        vgs = var(bgs.V, 'vgs')
        vds = var(bds.V, 'vds')
        vsb = var(bsb.V, 'vsb')
        vgb = var(vgs + vsb, 'vgb')
        vdb = var(vds + vsb, 'vdb')
        vgp = var(vgb - vtoT + phiT + gamma * sympy.sqrt(phiT), 'vgp')  # noqa
        vgr = var(sympy.sqrt(maxc(vgp, 0.0) + 0.25 * gamma ** 2     # noqa
                             + 1e-12), 'vgr')
        vp = var(sympy.Piecewise(
            (vgp - phiT - gamma * (vgr - 0.5 * gamma), vgp > 0.0),  # noqa
            (-phiT, True)), 'vp')
        nq = var(1.0 + gamma / (2.0 * sympy.sqrt(vp + phiT          # noqa
                                                 + 4.0 * ut)), 'nq')
        beta = var(kpT * w / l / (1.0 + theta * maxc(vp, 0.0)),     # noqa
                   'beta')
        ispec = var(2.0 * nq * beta * ut * ut, 'ispec')
        sf = var(softplus(var((vp - vsb) / (2.0 * ut), 'xf')), 'sf')
        sr = var(softplus(var((vp - vdb) / (2.0 * ut), 'xr')), 'sr')
        iff = var(sf * sf, 'iff')
        irr = var(sr * sr, 'irr')
        ids = var(ispec * (iff - irr), 'ids')
        sxf = var(sympy.sqrt(0.25 + iff), 'sxf')
        sxr = var(sympy.sqrt(0.25 + irr), 'sxr')
        sden = var((sxf + sxr) ** 2, 'sden')
        qdn = var(4.0 / 15.0 * (3.0 * sxr ** 3 + 6.0 * sxr ** 2 * sxf
                                + 4.0 * sxr * sxf ** 2 + 2.0 * sxf ** 3)
                  / sden - 0.5, 'qdn')
        qsn = var(4.0 / 15.0 * (3.0 * sxf ** 3 + 6.0 * sxf ** 2 * sxr
                                + 4.0 * sxf * sxr ** 2 + 2.0 * sxr ** 3)
                  / sden - 0.5, 'qsn')
        cwl = var(w * l * cox, 'cwl')                               # noqa
        qd = var(-cwl * nq * ut * qdn, 'qd')
        qs = var(-cwl * nq * ut * qsn, 'qs')
        qb = var(-cwl * gamma * sympy.sqrt(                         # noqa
            hypsmooth(phiT + vp, 1e-6 * phi))                       # noqa
            - (nq - 1.0) / nq * (qd + qs), 'qb')
        qg = var(-(qd + qs + qb), 'qg')
        return (Contribution(bds.I, ids),
                Contribution(bdb.I, ddt(qd)),
                Contribution(bgb.I, ddt(qg)),
                Contribution(bsb.I, ddt(qs)))
    return analog


## ----------------------------------------------------------------------
## 2.6  Temperature, and the p-channel

def test_ekv_temperature_path_moves_threshold_and_mobility():
    """EKV 2.6's temperature model: ``VTO(T) = VTO - TCV*(T - Tnom)``,
    ``KP(T) = KP*(T/Tnom)^BEX`` and the usual band-gap expression for
    ``PHI``.  Each is checked by MEASURING the quantity it moves --
    the threshold from an extrapolated strong-inversion intercept, the
    mobility from the saturation current at fixed overdrive -- rather
    than by reading the parameter back.
    """
    from pycircuit.circuit.circuit import ParameterDict
    card = dict(EKV_IDEAL, tcv=1.5e-3, bex=-1.5, tnom=300.0)
    el = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **card)

    def epar_at(T):
        e = defaultepar.copy()
        e.T = T
        return e

    def ids_at(T, vg, vd=6.0):
        return float(np.asarray(
            el.i(np.array([vd, vg, 0.0, 0.0]), epar_at(T)), float)[0])

    ## Threshold, by the classic extrapolation: in saturation
    ## sqrt(ID) is linear in VG, so the intercept is VT.
    def vt_extrapolated(T):
        v1, v2 = 2.0, 3.0
        s1, s2 = math.sqrt(ids_at(T, v1)), math.sqrt(ids_at(T, v2))
        return v1 - s1 * (v2 - v1) / (s2 - s1)

    ## `tcv` is a threshold FALL with temperature; over 100 K it is
    ## 150 mV, and that is what the extrapolation must see.
    d = vt_extrapolated(400.0) - vt_extrapolated(300.0)
    assert_allclose(d, -card['tcv'] * 100.0, rtol=6e-2)
    ## and it is the dominant term: with tcv = 0 the same measurement
    ## moves by an order of magnitude less (only PHI's band-gap drift).
    el0 = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b',
              **dict(card, tcv=0.0))

    def vt0(T):
        v1, v2 = 2.0, 3.0
        f = lambda vg: float(np.asarray(          # noqa: E731
            el0.i(np.array([6.0, vg, 0.0, 0.0]), epar_at(T)), float)[0])
        s1, s2 = math.sqrt(f(v1)), math.sqrt(f(v2))
        return v1 - s1 * (v2 - v1) / (s2 - s1)
    assert abs(vt0(400.0) - vt0(300.0)) < 0.25 * abs(d)

    ## Mobility: at a FIXED overdrive above the (moved) threshold the
    ## saturation current goes as KP, i.e. as (T/Tnom)^BEX.
    for T in (250.0, 350.0, 400.0):
        vt_T = vt_extrapolated(T)
        vt_0 = vt_extrapolated(300.0)
        i_T = ids_at(T, vt_T + 2.0)
        i_0 = ids_at(300.0, vt_0 + 2.0)
        assert_allclose(i_T / i_0, (T / 300.0) ** card['bex'], rtol=2e-2)


def test_ekv_pmos_is_the_nmos_with_every_voltage_and_current_reversed():
    """Same card, positive magnitudes, and the polarity carried entirely
    by the branch declarations -- so the mirror is exact rather than
    approximate."""
    n = _mk(eh.EkvNmosHdl, 'd', 'g', 's', 'b', **EKV)
    p = _mk(eh.EkvPmosHdl, 'd', 'g', 's', 'b', **EKV)
    assert p.terminals == ['d', 'g', 's', 'b']
    rng = np.random.default_rng(31)
    for _ in range(200):
        x = rng.uniform(-5.0, 5.0, 4)
        assert (np.asarray(n.i(x), float)
                == -np.asarray(p.i(-x), float)).all(), x
        assert (np.asarray(n.q(x), float)
                == -np.asarray(p.q(-x), float)).all(), x
        assert (np.asarray(n.G(x), float)
                == np.asarray(p.G(-x), float)).all(), x
    ## The limiter probes are declared on the reversed branches, so the
    ## p-channel's `fetlim` bounds `V(s,g)` -- the quantity that has the
    ## turn-on in it -- and its write-back moves the SOURCE.
    assert [(s[0], s[1], s[2]) for s in p._hdl_info['limit_spec']] == \
        [((2, 1), 'fet', 2), ((2, 0), 'vds', 0)]


def test_ekv_pmos_conducts_in_a_dc_solve():
    """A p-channel source follower from a 5 V rail: a real solve, whose
    answer is checked against the n-channel mirror of the same
    circuit."""
    card = dict(EKV, w=20e-6, l=2e-6)
    c = SubCircuit()
    ns = c.add_node('ns')
    ng = c.add_node('ng')
    nvdd = c.add_node('nvdd')
    c['vdd'] = VS(nvdd, gnd, v=5.0)
    c['vg'] = VS(ng, gnd, v=2.0)
    ## A p-channel follower needs its source PULLED UP: the drain goes
    ## to the negative rail and the load hangs from Vdd.
    c['rs'] = R(nvdd, ns, r=5e3)
    c['m'] = eh.EkvPmosHdl(gnd, ng, ns, nvdd, **card)
    vp = float(DC(c).solve().v(ns, gnd))
    ## The n-channel mirror: every rail and every node negated about
    ## the 5 V supply.
    c2 = SubCircuit()
    ns2 = c2.add_node('ns')
    ng2 = c2.add_node('ng')
    nvss = c2.add_node('nvss')
    c2['vss'] = VS(nvss, gnd, v=-5.0)
    c2['vg'] = VS(ng2, gnd, v=-2.0)
    c2['rs'] = R(nvss, ns2, r=5e3)
    c2['m'] = eh.EkvNmosHdl(gnd, ng2, ns2, nvss, **card)
    vn = float(DC(c2).solve().v(ns2, gnd))
    assert_allclose(vp, -vn, rtol=1e-9)
    ## and it really conducted: the source sits a threshold below the
    ## gate and the resistor carries a current.
    assert 2.0 < vp < 5.0, vp


## ======================================================================
## 3.  EKV against the surface-potential kernel
## ======================================================================
##
## `compact.PspMosLongChannel` solves the surface potential at each
## channel end and assembles the drain current from them; its drain
## current agrees with IHP's own compiled PSP103 to 1.3e-6 at the worst
## point of twelve sweeps.  EKV linearises the depletion charge about
## pinch-off instead.  Both are long-channel MOSFETs of the same
## physical device, so where they agree the agreement says something
## about both, and where they disagree the disagreement has a name.
##
## The comparison is set up so that NOTHING IS FITTED.  EKV's four core
## parameters are computed from PSP's physical card by the textbook
## identities:
##
##     Cox   = eps_ox/tox
##     GAMMA = sqrt(2*q*eps_si*NSUB)/Cox
##     PHI   = PHIB              (PSP's surface potential at threshold)
##     VTO   = VFB + PHI + GAMMA*sqrt(PHI)
##     KP    = U0*Cox
##
## and PSP is run with `mue = cs = thesat = 0`, which its own docstring
## records as recovering the ideal long-channel core, plus every other
## correction (CLM, DIBL, quantum, poly depletion, junction, gate
## resistance) switched off.  Importing and compiling `compact` costs
## about a minute, so everything here shares one module-scoped fixture.

PSP_CARD = dict(w=10e-6, l=10e-6, tox=2.2e-9, nsub=5e23, vfb=-0.95,
                u0=0.045, phib=0.9,
                mue=0.0, cs=0.0, thesat=0.0, ct=0.0, alp=0.0, qq=0.0,
                kp=0.0, dnsub=0.0, cf=0.0, thesatb=0.0, thesatg=0.0,
                xcor=0.0, rs=0.0, rg=0.0, mult=1.0, ax=2.0)


@pytest.fixture(scope='module')
def psp_pair():
    """One PSP element and the EKV element derived from its card.

    Module-scoped: compiling `compact` takes about a minute and there is
    nothing per-test about either element.
    """
    from pycircuit.circuit import compact
    from pycircuit.circuit.psp_scaling import (PSP_EPS_SI, PSP_QELE,
                                               PSP_EPS0)
    old = pycircuit.circuit.circuit.default_toolkit
    pycircuit.circuit.circuit.default_toolkit = numeric
    try:
        p = compact.PspMosLongChannel(Node('d'), Node('g'), Node('s'),
                                      Node('b'), **PSP_CARD)
        p.update_iparv()
        cox = 3.9 * PSP_EPS0 / PSP_CARD['tox']
        gamma = math.sqrt(2 * PSP_QELE * PSP_EPS_SI
                          * PSP_CARD['nsub']) / cox
        phi = PSP_CARD['phib']
        card = dict(vto=PSP_CARD['vfb'] + phi + gamma * math.sqrt(phi),
                    gamma=gamma, phi=phi, kp=PSP_CARD['u0'] * cox,
                    cox=cox, w=PSP_CARD['w'], l=PSP_CARD['l'])
        e = eh.EkvNmosHdl(Node('d'), Node('g'), Node('s'), Node('b'),
                          **card)
        e.update_iparv()
    finally:
        pycircuit.circuit.circuit.default_toolkit = old
    return p, e, card


def _psp_id(p, vd, vg, vs=0.0, vb=0.0):
    return float(np.asarray(p.i(p.bias(vd, vg, vs, vb)), float)[0])


def test_the_derived_ekv_card_is_the_textbook_mapping(psp_pair):
    """The parameter mapping, stated as numbers so that a reader can
    check it and so that a change in either model's constants shows up
    here rather than as a mysterious shift downstream."""
    _p, _e, card = psp_pair
    assert_allclose(card['gamma'], 0.260665, rtol=1e-4)
    assert_allclose(card['vto'], 0.197288, rtol=1e-3)
    assert_allclose(card['kp'], 7.0632e-4, rtol=1e-4)
    assert_allclose(card['cox'], 0.0156961, rtol=1e-4)


def test_ekv_and_psp_agree_on_the_subthreshold_slope_factor(psp_pair):
    """**Parameter-free comparison 1.**  The subthreshold swing depends
    on the body factor and on nothing that was mapped by hand: ``VTO``
    cancels out of a slope.

    Measured: PSP 69.6 mV/decade, EKV 72.6 -- 4.2% apart.  That gap is
    the price of EKV's linearised depletion charge against PSP's exact
    surface potential, and 4% is a good result for a formulation one
    tenth the size.
    """
    p, e, _card = psp_pair

    def swing(f, v0=0.05, v1=0.10):
        return (v1 - v0) / math.log10(f(v1) / f(v0))

    s_psp = swing(lambda vg: _psp_id(p, 1.0, vg))
    s_ekv = swing(lambda vg: _ids(e, 1.0, vg))
    ## Both are genuinely subthreshold at those gates.
    assert _psp_id(p, 1.0, 0.10) < 1e-7
    assert 0.068 < s_psp < 0.071, s_psp
    assert 0.071 < s_ekv < 0.074, s_ekv
    assert 1.03 < s_ekv / s_psp < 1.06, s_ekv / s_psp


def test_ekv_and_psp_agree_on_the_body_effect(psp_pair):
    """**Parameter-free comparison 2.**  How far the threshold moves
    with body bias is a property of the depletion charge, and it too is
    independent of ``VTO``: the measurement is a DIFFERENCE of two gate
    voltages at equal current.

    Measured, for 1.5 V of source-body bias: PSP 178 mV, EKV 163 mV --
    8.4% apart.  So the two formulations differ by 235 mV in the
    ABSOLUTE threshold (next test) and by 15 mV in how it MOVES: the
    linearised depletion charge gets the body effect nearly right and
    the threshold itself badly wrong, which is the expected signature
    and not an obvious one.
    """
    p, e, _card = psp_pair

    def vg_at(f, target, vsb):
        ## `f` takes a SOURCE-REFERRED gate voltage, so the bracket is
        ## source-referred too and nothing is added or subtracted.
        lo, hi = -1.0, 2.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if f(mid, vsb) < target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    fp = lambda vgs, vsb: _psp_id(p, 1.0 + vsb, vgs + vsb, vsb, 0.0)  # noqa
    fe = lambda vgs, vsb: _ids(e, 1.0 + vsb, vgs + vsb, vsb, 0.0)     # noqa
    ## A weak-inversion current, so what is measured is a threshold
    ## shift and not a mobility difference.
    tgt_p = fp(0.10, 0.0)
    tgt_e = fe(0.10, 0.0)
    assert tgt_p < 1e-7 and tgt_e < 1e-6
    dp = vg_at(fp, tgt_p, 1.5) - vg_at(fp, tgt_p, 0.0)
    de = vg_at(fe, tgt_e, 1.5) - vg_at(fe, tgt_e, 0.0)
    assert 0.17 < dp < 0.19, dp
    assert 0.15 < de < 0.18, de
    ## Asserted from both sides: 8.4%, so neither a regression nor an
    ## unexplained improvement passes silently.
    assert 0.05 < abs(de / dp - 1.0) < 0.12, (dp, de)


def test_ekv_and_psp_disagree_on_the_absolute_threshold(psp_pair):
    """**The result this cross-check exists to produce, and it is a
    disagreement.**

    With EKV's ``VTO`` computed from PSP's card by the textbook identity
    ``VFB + PHI + GAMMA*sqrt(PHI)``, the two models do not agree, and
    the gap is not a constant offset either:

    * in DEEP TRIODE the ratio ``I_EKV/I_PSP`` is a pure 235 mV
      threshold shift -- measured as 0.235, 0.235, 0.235, 0.241 V at
      four gate voltages a decade apart in overdrive, which is what
      makes "threshold shift" the right description rather than a
      guess;
    * in WEAK INVERSION the same measurement gives 72 mV.

    So the two formulations want thresholds that differ by 160 mV
    between the two regions.  That is exactly the linearisation: PSP's
    surface potential keeps rising above ``PHIB`` as the gate drives
    harder -- by about nine thermal voltages here -- while EKV pins it
    at ``PHI + V`` and absorbs the difference into ``VTO``.

    **The practical reading: EKV's ``VTO`` is an EXTRACTED parameter,
    not a computed one.**  Deriving it from a doping and a flat-band
    voltage, which is the obvious thing to do and is what this test
    does, is wrong by a quarter of a volt, and the error is
    region-dependent so no single extracted value removes it either.
    Nothing here is tuned to hide that.
    """
    p, e, card = psp_pair
    ## Deep triode: `I ~ n*beta*(VP - VS)*Vds`, so a threshold shift
    ## `dV` shows up as `(Vgs - Vt + dV)/(Vgs - Vt)`.
    shifts = []
    for vgs in (0.5, 0.8, 1.2, 1.8):
        r = _ids(e, 0.05, vgs) / _psp_id(p, 0.05, vgs)
        shifts.append((r - 1.0) * (vgs - card['vto']))
    assert all(0.230 < sft < 0.245 for sft in shifts), shifts
    ## The four readings agree with each other to 3%, which is what
    ## says this really is one threshold offset and not a shape error.
    assert max(shifts) / min(shifts) < 1.03, shifts

    ## Weak inversion: the shift is `n*UT*ln(ratio)` with the measured
    ## slope factor.
    ratio = _ids(e, 1.0, 0.10) / _psp_id(p, 1.0, 0.10)
    n_meas = 0.0696 / (math.log(10.0) * _UT)      # PSP's own swing
    weak_shift = n_meas * _UT * math.log(ratio)
    assert 11.0 < ratio < 12.5, ratio
    assert 0.065 < weak_shift < 0.080, weak_shift
    ## and the two shifts are far apart -- three times, 160 mV.
    assert min(shifts) - weak_shift > 0.15


def test_ekv_and_psp_agree_on_shape_once_one_threshold_is_aligned(psp_pair):
    """Having found that ``VTO`` cannot be computed, ALIGN IT ONCE and
    ask what is left.

    One degree of freedom -- the offset that makes the two currents
    equal at a single bias, ``Vgs = Vds = 1.2 V`` -- is removed, and
    nothing else is touched.  Over the strong-inversion part of a 2-D
    sweep (four gate voltages, four drain voltages, three body biases)
    the two then agree to **within 15%**, with the worst point at the
    edge of saturation.

    This is a legitimate cross-check because it is one parameter, it is
    the parameter EKV defines by extraction, and it is stated.  It is
    NOT a fit: no shape parameter is adjusted, and the residual is
    reported rather than tolerated.
    """
    p, _e, card = psp_pair
    target = _psp_id(p, 1.2, 1.2)

    def with_offset(dv):
        el = eh.EkvNmosHdl(Node('d'), Node('g'), Node('s'), Node('b'),
                           **dict(card, vto=card['vto'] + dv))
        el.update_iparv()
        return el

    lo, hi = -0.5, 1.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if _ids(with_offset(mid), 1.2, 1.2) > target:
            lo = mid
        else:
            hi = mid
    dv = 0.5 * (lo + hi)
    ## The offset needed is the strong-inversion shift, ~176 mV.
    assert 0.16 < dv < 0.19, dv
    el = with_offset(dv)

    strong, s_where = 0.0, None
    allpts, a_where = 0.0, None
    for vsb in (0.0, 0.5, 1.5):
        for vgs in (0.8, 1.2, 1.8):
            for vds in (0.05, 0.3, 1.2, 2.5):
                ip = _psp_id(p, vds + vsb, vgs + vsb, vsb, 0.0)
                ie = _ids(el, vds + vsb, vgs + vsb, vsb, 0.0)
                assert ip > 0 and ie > 0
                d = abs(math.log(ie / ip))
                if d > allpts:
                    allpts, a_where = d, (vsb, vgs, vds, ip, ie)
                ## Strong inversion only.  At Vgs = 0.8 with 1.5 V of
                ## body bias the device is in MODERATE inversion, where
                ## the effective threshold offset is somewhere between
                ## its weak and strong values -- so that corner is not
                ## testing the square law, it is testing the very
                ## region-dependence the previous test measured.
                if vgs >= 1.2:
                    if d > strong:
                        strong, s_where = d, (vsb, vgs, vds, ip, ie)
    ## 14% at the worst of the 24 strong-inversion points, and the band
    ## is asserted from BOTH sides so that an improvement is noticed too.
    assert strong < 0.17, s_where
    assert strong > 0.05, ('suspiciously good -- check the mapping',
                           s_where)
    ## and the moderate-inversion corner, reported rather than hidden:
    ## 36% at Vgs = 0.8 with 1.5 V of body bias.
    assert 0.30 < allpts < 0.45, a_where
