"""PSP103's geometry scaling, for the parameters the core model uses.

A PSP model card does not hold the numbers the equations want.  It holds
*geometry-independent* coefficients, and PSP's scaling layer
(`PSP103_scaling.include`, 849 lines) turns those plus a device's W and
L into the local parameters the physics is written in.  Without it a
card cannot be connected to a model at all.

This implements the subset `compact.PspMosLongChannel` needs -- flat-band
voltage, effective doping, surface potential at threshold, the gain
factor, mobility reduction and velocity saturation -- at nominal
temperature and with no stress or well-proximity corrections.  Each
formula cites the line it came from.

What is deliberately absent, and would matter for a short device: the
short-channel and DIBL parameters, series resistance, overlap and
fringe capacitances, gate and junction leakage, and every temperature
coefficient.  This exists to answer one question -- how far is the core
from the real device on a LONG one, where those matter least -- not to
be a complete scaling layer.

Usage::

    card = spicecard.read('cornerMOSlv.lib', section='mos_tt')
    p = card.model_params('sg13g2_lv_nmos_psp', w=10e-6, l=1e-6, ng=1)
    kw = psp_scaling.to_long_channel(p, w=10e-6, l=1e-6)
    fet = PspMosLongChannel(d, g, s, b, **kw)
"""

import math

from pycircuit.circuit.constants import (eps0, epsRSiO2,
                                         kboltzmann, qelectron)


#: PSP's normalisation lengths (`PSP103_macrodefs.include:58-59`).
LEN = 1.0e-6
WEN = 1.0e-6

#: Quantum-mechanical correction constants
#: (`PSP103_macrodefs.include:51-52`), n- and p-channel.
## PSP'S OWN PHYSICAL CONSTANTS, and they are NOT the tree's.
##
## `Common103_macrodefs.include:60-61` sets `EPSO = 8.8541878176e-12`
## and **`EPSRSI = 11.8`**, where `pycircuit.circuit.constants` carries
## 11.7.  Both values are defensible for silicon; what is not defensible
## is comparing against a reference that used one while using the other.
##
## The body factor goes as `sqrt(eps_si)`, so 11.7 against 11.8 makes
## our `gamma` 0.43% small -- and because the body term grows as
## `sqrt(phib + Vsb)`, the resulting threshold error GROWS WITH BODY
## BIAS.  Measured: 2 mV on the grounded long device, 3.6 mV on the
## grounded short one, 5.8 mV with `Vsb = 1`.  At 85 mV/decade that is
## 5-11% of the weak-inversion current, which is precisely the residual
## that survived every other check.
##
## It survived them because NO `lp_*` OUTPUT EXPOSES THE BODY FACTOR.
## Thirty-one scaled parameters matched exactly, and the one number
## between them and the current was a constant nobody had thought to
## compare.
## `KBOL` and `QELE` are on the same list (`:56-57`) and worth far less
## -- 0.047% and 0.012%, about 0.2 mV of threshold between them -- but
## they are free, and a model compared against a reference should not
## be using different physical constants from it at all.
PSP_EPS0 = 8.8541878176e-12
PSP_EPSRSI = 11.8
PSP_EPS_SI = PSP_EPSRSI * PSP_EPS0
PSP_KBOL = 1.3806505e-23
PSP_QELE = 1.6021918e-19
PSP_HBAR = 1.05457168e-34
PSP_MELE = 9.1093826e-31

QMN = 5.951993
QMP = 7.448711


def _g(card, name, default=0.0):
    """A card parameter, or its PSP default when the card omits it."""
    v = card.get(name.lower())
    return default if v is None else float(v)


def geometry(card, w, l):
    """Effective device dimensions and their reciprocals.

    `PSP103_scaling.include:28-43`.  The effective length is the drawn
    one less twice the lateral diffusion; `iLE`/`iWE` are normalised
    reciprocals, and note they are NOT small for a 1 um device -- the
    normalisation length is 1 um, so `iLE` is of order 1 and the
    length-dependent terms are fully active even on what a designer
    would call a long transistor.
    """
    iL, iW = LEN / l, WEN / w
    dl = _g(card, 'lvaro') * (1.0 + _g(card, 'lvarl') * iL) \
        * (1.0 + _g(card, 'lvarw') * iW)
    dw = _g(card, 'wvaro') * (1.0 + _g(card, 'wvarl') * iL) \
        * (1.0 + _g(card, 'wvarw') * iW)
    LE = max(l + dl - 2.0 * _g(card, 'lap'), 1.0e-9)
    WE = max(w + dw - 2.0 * _g(card, 'wot'), 1.0e-9)
    ## `Lcv`/`Wcv` are the CV dimensions WITHOUT the lateral-diffusion
    ## subtraction (`PSP103_scaling.include:40-41`), which is what the
    ## fringe and gate-bulk overlap terms scale on; `LEcv`/`WEcv` (:38-39)
    ## are the ones with it, for the oxide capacitance and the gate
    ## overlap.  Two different effective lengths, a line apart in the
    ## vendor source, and easy to swap by accident.
    dlq, dwq = _g(card, 'dlq'), _g(card, 'dwq')
    return dict(LE=LE, WE=WE, iLE=LEN / LE, iWE=WEN / WE,
                LEcv=max(LE + dlq, 1.0e-9), WEcv=max(WE + dwq, 1.0e-9),
                Lcv=max(l + dl + dlq, 1.0e-9),
                Wcv=max(w + dw + dwq, 1.0e-9))


def _neff(card, geo):
    """Effective channel doping, including the pocket implant.

    `PSP103_scaling.include:235-250`.  The three-way branch is on how the
    channel length compares with the pocket length: a long channel sees
    the pockets as a perturbation (the logarithmic form), a short one is
    all pocket.
    """
    iWE, LE = geo['iWE'], geo['LE']
    nsub0 = _g(card, 'nsubo') * max(
        1.0 + _g(card, 'nsubw') * iWE
        * math.log(1.0 + geo['WE'] / _g(card, 'wseg', 1.0)), 1.0e-3)
    npck = _g(card, 'npck') * max(
        1.0 + _g(card, 'npckw') * iWE
        * math.log(1.0 + geo['WE'] / _g(card, 'wsegp', 1.0)), 1.0e-3)
    lpck = _g(card, 'lpck') * max(
        1.0 + _g(card, 'lpckw') * iWE
        * math.log(1.0 + geo['WE'] / _g(card, 'wsegp', 1.0)), 1.0e-3)

    if LE > 2.0 * lpck:
        aa = 7.5e10
        bb = math.sqrt(nsub0 + 0.5 * npck) - math.sqrt(nsub0)
        nsub = math.sqrt(nsub0) + aa * math.log(
            1.0 + 2.0 * lpck / LE * (math.exp(bb / aa) - 1.0))
        nsub = nsub * nsub
    elif LE >= lpck:
        nsub = nsub0 + npck * lpck / LE
    else:
        nsub = nsub0 + npck * (2.0 - LE / lpck)

    iLE = geo['iLE']
    return nsub * (1.0 - _g(card, 'fol1') * iLE
                   - _g(card, 'fol2') * iLE * iLE)


def channel_type(card):
    """+1 for an n-channel card, -1 for a p-channel one.

    Deliberately NOT part of `to_long_channel`'s return value.  That
    dict is splatted straight into the element as keyword arguments, so
    anything in it has to be a parameter the element actually has, and
    the channel type is not one -- it selects the CLASS.  Returning it
    there would have turned every existing call site into a
    `TypeError`, which is a poor way to learn about a new feature.

    The card carries the polarity here and nowhere else: every parameter
    on IHP's p-channel card is a positive magnitude with the same signs
    as the n-channel one.
    """
    return -1.0 if _g(card, 'type', 1.0) < 0 else 1.0


def to_long_channel(card, w, l, T=300.0, all_terms=True):
    """Card + geometry -> `PspMosLongChannel` keyword arguments.

    `card` is what `spicecard.Deck.model_params` returned.  Returns a
    dict ready to splat into the element.

    `all_terms=False` drops two blocks, and exists only so that their
    effect can be measured:

    * **DIBL** -- `CF`, `CFB`, `CFD` (`macrodefs:473-476`);
    * the **body- and gate-bias modulation of `THESAT`** -- `THESATB`,
      `THESATG` (`macrodefs:596-607`).

    They spent one commit switched off, on a measurement that turned out
    to be the wrong measurement.  The correction is worth recording,
    because the metric was the mistake rather than the physics.

    Judged by the MEDIAN ratio to PSP over each sweep they look bad:
    summed error over the six n-channel sweeps 0.016 without them and
    0.041 with.  Judged by the SPREAD of that ratio across a sweep --
    which is what says whether the bias dependence is right, the median
    saying only whether the gain is -- they are decisive: summed spread
    over eleven sweeps 1.80 without them and 0.38 with.

    The p-channel is what made it obvious, needing DIBL far more than
    the n-channel does.  On `pmos_idvg_vd1p2` the ratio without these
    terms sweeps from 1.03 down to 0.44 -- a drive that is simply wrong
    -- and with them it is FLAT AT 1.086 across three decades of
    current, spread 1.353 -> 0.020.  A flat ratio is a gain error: one
    cause, one number, a well-posed question.  A ratio that sweeps is a
    broken bias dependence, and no amount of a favourable median makes
    that the better model.

    So the better n-channel median without them was compensating
    errors, exactly as the pattern this model has now hit three times
    predicts -- and this time the compensation was visible only because
    a second channel type disagreed.

    All twenty-one scaled parameters the reference records match PSP's
    own `lp_*` outputs exactly at four geometries, on BOTH channel
    types, so none of this is the scaling layer.  What is left is a gain
    offset -- about 2% on the n-channel and 2-9% on the p-channel -- and
    that is the next thing to chase.
    """
    ## Channel type.  The card carries the polarity HERE and nowhere
    ## else: every parameter on the p-channel card is a positive
    ## magnitude with the same signs as the n-channel one (`vfbo` is
    ## negative on both).  PSP does the same -- `PSP103_scaling.include`
    ## contains not one `CHNL_TYPE` reference, the entire scaling layer
    ## being channel-type agnostic.
    chtype = channel_type(card)

    geo = geometry(card, w, l)
    iLE, iWE, LE, WE = geo['iLE'], geo['iWE'], geo['LE'], geo['WE']

    tox = _g(card, 'toxo')
    eps_ox = _g(card, 'epsroxo', epsRSiO2) * PSP_EPS0
    cox_prime = eps_ox / tox

    ## Flat-band voltage (:230) and effective doping (:250).
    vfb = (_g(card, 'vfbo') + _g(card, 'vfbl') * iLE
           + _g(card, 'vfbw') * iWE + _g(card, 'vfblw') * iLE * iWE)
    neff = _neff(card, geo)

    ## Surface potential at threshold (`PSP103_macrodefs.include:295-308`).
    ## `phibFac` and `Eg` are the temperature-dependent parts; at 300 K
    ## they are what they are, and no temperature COEFFICIENT from the
    ## card is applied here.
    phit = PSP_KBOL * T / PSP_QELE
    eg = 1.179 - 9.025e-5 * T - 3.05e-7 * T * T
    phib_fac = max((1.045 + 4.5e-4 * T)
                   * (0.523 + 1.4e-3 * T - 1.48e-6 * T * T)
                   * T * T / 9.0e4, 1.0e-3)
    dphib = (_g(card, 'dphibo')
             + _g(card, 'dphibl') * iLE ** _g(card, 'dphiblexp', 1.0)
             + _g(card, 'dphibw') * iWE
             + _g(card, 'dphiblw') * iLE * iWE)
    phib = max(eg + dphib
               + 2.0 * phit * math.log(neff * phib_fac ** -0.75 * 4.0e-26),
               5.0e-2)

    ## QUANTUM-MECHANICAL CORRECTION
    ## (`PSP103_macrodefs.include:322-327`, switched by the card's QMC).
    ##
    ## Carriers in the inversion layer occupy quantised states a little
    ## below the surface, so the effective surface potential at threshold
    ## is higher and the body factor larger.  For a 2.2 nm oxide the
    ## shift is tens of millivolts -- small-looking, and it moves the
    ## drain current by ~30% at a fixed gate bias, which is how it was
    ## found: the first comparison against PSP103 ran that much high with
    ## a ratio that grew with current, and this is the only term of that
    ## size the core was missing.
    ##
    ## The correction is applied to `phib` and to `G_0`, so it is
    ## returned as an adjusted doping rather than an adjusted gamma --
    ## the element derives its body factor from `nsub`, and scaling that
    ## by the square of the `G_0` factor is the equivalent move.
    gamma0 = math.sqrt(2.0 * PSP_QELE * neff * PSP_EPS_SI
                       / phit) / cox_prime
    qmc = _g(card, 'qmc')
    qq = 0.0
    if qmc > 0.0:
        ## `QMN` for electrons, `QMP` for holes
        ## (`PSP103_module.include:727-734`).  A ratio of 1.25, and not
        ## cosmetic: `qq` shifts both `phib` and the body factor, so it
        ## moves the threshold and the body effect together.
        qq = (0.4 * (QMP if chtype < 0 else QMN)
              * qmc * cox_prime ** (2.0 / 3.0))
        qb0 = math.sqrt(phit * gamma0 * gamma0 * phib)
        dphibq = 0.75 * qq * qb0 ** (2.0 / 3.0)
        phib = phib + dphibq
        g_fac = 1.0 + 2.0 * (2.0 / 3.0) * dphibq / qb0
        neff = neff * g_fac * g_fac

    ## Gain factor (:286-289).  `GPE` and `GWE` are the length and width
    ## corrections to mobility; `BETN = UO*WE/(GPE*LE)*GWE`, and
    ## `BET = BETN * Cox'` (`PSP103_macrodefs.include:368`).
    ## `FBET1` and `LP1` are WIDTH-ADJUSTED before use (:284-285) -- the
    ## trailing `e` in the source is not decoration.  Using the raw card
    ## values made `GPE` 1.271 where PSP computes 1.4234, and since
    ## `BETN = UO*WE*GWE/(GPE*LE)` that is a 12% gain error, which is
    ## most of what the short device was off by.  Caught by comparing
    ## against PSP's own `lp_betn` operating-point output.
    fbet1e = _g(card, 'fbet1') * (1.0 + _g(card, 'fbet1w') * iWE)
    lp1e = _g(card, 'lp1', 1.0) * max(1.0 + _g(card, 'lp1w') * iWE, 1.0e-3)
    lp2 = _g(card, 'lp2', 1.0)
    gpe = max(1.0 + fbet1e * lp1e / LE * (1.0 - math.exp(-LE / lp1e))
              + _g(card, 'fbet2') * lp2 / LE
              * (1.0 - math.exp(-LE / lp2)), 1.0e-15)
    gwe = (1.0 + _g(card, 'betw1') * iWE
           + _g(card, 'betw2') * iWE
           * math.log(1.0 + WE / _g(card, 'wbet', 1.0)))

    ## The element computes `beta = u0 * cox * w / l` from the DRAWN
    ## dimensions, so the effective-dimension and GPE/GWE corrections are
    ## folded into an effective mobility rather than applied twice.
    u0_eff = _g(card, 'uo') * (WE / w) * (l / LE) * gwe / gpe

    ## Mobility reduction (:291-301) and velocity saturation (:309).
    mue = _g(card, 'mueo') * (1.0 + _g(card, 'muew') * iWE)
    cs = ((_g(card, 'cso')
           + _g(card, 'csl') * iLE ** _g(card, 'cslexp', 1.0))
          * (1.0 + _g(card, 'csw') * iWE)
          * (1.0 + _g(card, 'cslw') * iLE * iWE))
    thesat = ((_g(card, 'thesato')
               + _g(card, 'thesatl') * gwe / gpe
               * iLE ** _g(card, 'thesatlexp', 1.0))
              * (1.0 + _g(card, 'thesatw') * iWE)
              * (1.0 + _g(card, 'thesatlw') * iLE * iWE))

    ## POLYSILICON DEPLETION.  `kp = 2*Cox'^2*phit/(q*np*eps_si)`, with
    ## `np` floored at both `8e7/tox^2` and 5e24
    ## (`PSP103_macrodefs.include:315-319`).  The card's `NPO` is
    ## 4.6e26 -- not the zero that would switch the effect off.
    ##
    ## AND IT IS GEOMETRY-SCALED: `NP = NPO*max(1e-6, 1 + NPL*iLE)`
    ## (`PSP103_scaling.include:260`), clipped low at zero (`:705`).
    ## That went missing for the same reason the `BETN` width adjustment
    ## once did -- `NPL` is EXACTLY ZERO on the n-channel card, so the
    ## scaling is invisible on the only device it was checked against.
    ## The p-channel card sets `NPL = -0.0959`, which takes `NP` down to
    ## 0.36 of its nominal value on a 0.13 um device: nearly a factor of
    ## three more gate depletion than we were applying, and worth about
    ## 8% of the drain current there.
    np_card = _g(card, 'npo') * max(1.0e-6,
                                    1.0 + _g(card, 'npl') * iLE)
    np_card = max(np_card, 0.0)
    if np_card > 0.0:
        np_eff = max(max(np_card, 8.0e7 / (tox * tox)), 5.0e24)
        kp = (2.0 * cox_prime * cox_prime * phit
              / (PSP_QELE * np_eff * PSP_EPS_SI))
    else:
        kp = 0.0

    ## Channel-length modulation (:321, :327).  The geometry-scaled
    ## branch has NO geometry-independent term -- `ALP` is length
    ## dependence only.
    alp = (_g(card, 'alpl') * iLE ** _g(card, 'alplexp', 1.0)
           * (1.0 + _g(card, 'alpw') * iWE))
    ## `ALP1` and `ALP2` correct the channel-shortening factor in strong
    ## and WEAK inversion respectively (:323-326).  `ALP2` is 4.5 on this
    ## geometry -- assuming it zero was what left the near-threshold
    ## output conductance far too low.
    t1 = iLE ** _g(card, 'alp1lexp', 1.0)
    alp1 = (_g(card, 'alp1l1') * t1 * (1.0 + _g(card, 'alp1w') * iWE)
            / (1.0 + _g(card, 'alp1l2') * iLE * t1))
    t2 = iLE ** _g(card, 'alp2lexp', 1.0)
    alp2 = (_g(card, 'alp2l1') * t2 * (1.0 + _g(card, 'alp2w') * iWE)
            / (1.0 + _g(card, 'alp2l2') * iLE * t2))
    ## Drain-induced barrier lowering (:273-276), clipped at :714-717.
    ## `CF` carries a LENGTH POWER and nothing else -- it is a purely
    ## short-channel quantity, and on a 1 um device it scales to 1e-7,
    ## which is why omitting it cost nothing there and about a percent
    ## at 0.13 um.  PSP's companion `PSCE` block, which degrades the
    ## subthreshold slope, is all-zero on this card and is not built.
    cf = max(_g(card, 'cfl') * iLE ** _g(card, 'cflexp', 2.0)
             * (1.0 + _g(card, 'cfw') * iWE), 0.0)
    cfb = min(max(_g(card, 'cfbo'), 0.0), 1.0)
    cfd = max(_g(card, 'cfdo'), 0.0)

    ## THE CV EFFECTIVE DIMENSIONS (:38-39), which are NOT the DC ones.
    ## `LEcv = LE + DLQ`, `WEcv = WE + DWQ`, and the card sets both:
    ## `dlq = -1.37e-8`, `dwq = -1e-8` for this n-channel.  PSP builds
    ## its oxide capacitance for the CHARGE model out of these (:359),
    ## so using the drawn dimensions overstates every capacitance -- by
    ## 7% on the long device here, and the QM correction below adds
    ## another 12% on top of that.
    lecv, wecv = geo['LEcv'], geo['WEcv']

    ## OVERLAP AND FRINGE CAPACITANCE (:359-373).
    ##
    ## `CGOV` is the gate-to-source/drain overlap over the diffusion,
    ## `CFR` the outer fringe from the gate sidewall.  With
    ## `SWJUNASYM = 0` -- this card -- the drain side mirrors the source
    ## exactly (`:840-850`), which is what keeps the pair even under the
    ## source/drain exchange.  `CGBOV` is the gate-to-bulk overlap and is
    ## 4e-28 here, i.e. off; carried anyway because it is one term.
    ##
    ## The card sets no `FCGOVACC`, so the accumulation branches of the
    ## overlap charge (`module:1560-1596`) are absent entirely.
    lov = _g(card, 'lov')
    toxov = _g(card, 'toxovo', 2.0e-9)
    nov = min(max(_g(card, 'novo', 5.0e25), 1.0e23), 1.0e27)
    cgov = max(eps_ox * wecv * lov / toxov, 0.0)
    cfr = max(_g(card, 'cfrw') * geo['Wcv'] / WEN, 0.0)
    cgbov = max(_g(card, 'cgbovl') * geo['Lcv'] / LEN, 0.0)

    ## The overlap region has its own surface potential, and PSP solves
    ## it in CLOSED FORM (`module:1182-1189`) rather than iterating: one
    ## smoothed maximum and one explicit expression.  Everything the
    ## expression needs except the bias is fixed per instance, so it is
    ## computed here, in Python, and handed over as numbers -- the
    ## compiled expression stays four lines long.
    ##
    ## `sp_ovInit` (`macrodefs:217-235`) is a piecewise FIT in `1/GOV`,
    ## not a derivation; transcribed as such.
    coxov = eps_ox / toxov
    gov = math.sqrt(2.0 * PSP_QELE * nov * PSP_EPS_SI / phit) / coxov
    gov2 = gov * gov
    inv_gov = 1.0 / gov
    sp_eps = 3.1 * gov + 8.5
    sp_delta = 0.5 * sp_eps
    if inv_gov < 0.06:
        sp_a = 64.0 * inv_gov
    elif inv_gov <= 0.45:
        sp_a = 22.0 * inv_gov + 3.0
    elif inv_gov <= 1.6:
        sp_a = -7.2 * inv_gov + 15.5
    else:
        sp_a = gov
    sp_delta1 = (sp_delta + gov2 * 0.5
                 - gov * math.sqrt(sp_delta + gov2 * 0.25 + sp_a))

    ## CHANNEL NOISE (:376-384, clipped at :792-797).  The flicker
    ## coefficients carry a length correction of their own -- `Lnoi`
    ## and `ALPNOI` -- on top of the 1/(W*L) that `NFALW` already
    ## implies.  `NFC` clips to ZERO here: the card's `NFCLW` is
    ## negative and `NFC_i` is clipped low at 0 (`:796`), so the
    ## quadratic term is off however it is written.
    lnoi = max(1.0 - 2.0 * _g(card, 'lintnoi') * iLE / LEN, 1.0e-3)
    lred = 1.0 / lnoi ** _g(card, 'alpnoi', 0.0)
    nz_scale = lred * iWE * iLE
    nfa = max(nz_scale * _g(card, 'nfalw'), 0.0)
    nfb = max(nz_scale * _g(card, 'nfblw'), 0.0)
    nfc = max(nz_scale * _g(card, 'nfclw'), 0.0)
    ef = max(_g(card, 'efo', 1.0), 0.0)
    fnt = max(_g(card, 'fnto', 1.0), 0.0)
    ## `SWIGN` is a plain switch, not a scaled parameter: it is read
    ## off the card verbatim (`PSP103_parlist.include:48`).  Default 1,
    ## which is both PSP's default and what this card sets.
    swign = 1.0 if _g(card, 'swign', 1.0) >= 0.5 else 0.0

    ## GATE RESISTANCE (:604, clipped at :816).  The full expression
    ## carries a sheet-resistance term and a per-finger term; this card
    ## sets neither (`RSHG` and `RGO` are absent), leaving
    ## `RG = (RINT + RVPOLY)/(W*L)`, which reproduces PSP's own `lp_rg`
    ## of 1.3025 ohm exactly on a 10x1 um device.
    ##
    ## `RSE`, `RDE` and `RBULK` are all zero on this card -- the source
    ## and drain resistance PSP folds into the mobility instead, which
    ## this model already does.  So the gate is the only terminal
    ## resistance there is here.
    rg = max((_g(card, 'rint') + _g(card, 'rvpoly')) / (w * l), 0.0)

    ## BODY-BIAS MOBILITY CORRECTION (:299, clipped at :731).
    ## `Rxcor = (1 + 0.2*XCOR*Vsbx)/(1 + XCOR*Vsbx)` multiplies `Gmob`
    ## at both the source end and the midpoint (`macrodefs:576, 595,
    ## 750`), so it raises the mobility under reverse body bias.
    ##
    ## It is EXACTLY 1 at Vsb = 0, which is why it went unbuilt for so
    ## long: every sweep in the reference used a grounded body except
    ## one, and that one is on the geometry where `XCOR` scales to
    ## zero.  A term that no measurement touches is not a term that has
    ## been shown not to matter.
    xcor = max(_g(card, 'xcoro')
               * (1.0 + _g(card, 'xcorl') * iLE)
               * (1.0 + _g(card, 'xcorw') * iWE)
               * (1.0 + _g(card, 'xcorlw') * iLE * iWE), 0.0)

    ## Bias-dependent body factor (:255-257, no geometry term at all).
    ## `DNSUBO` is 4.4e-16 on the n-channel card -- zero in every sense
    ## that matters -- and 0.0397 on the p-channel one, so this is the
    ## THIRD term on this PDK that a zero coefficient hides from any
    ## n-channel-only measurement.
    dnsub = max(_g(card, 'dnsubo'), 0.0)
    vnsub = _g(card, 'vnsubo')
    nslp = max(_g(card, 'nslpo', 0.05), 1.0e-6)

    ## Body- and gate-bias modulation of the velocity-saturation
    ## parameter (:313-314).  Geometry-INDEPENDENT: PSP takes these
    ## straight from the card with no length or width term at all.
    ## Clips at :741-742.
    thesatb = min(max(_g(card, 'thesatbo'), -0.5), 1.0)
    thesatg = max(_g(card, 'thesatgo'), -0.5)

    ## Linear/saturation transition sharpness (:317), FLOORED AT 2
    ## (:743).  The floor is not cosmetic: this card scales `AX` to 0.88
    ## on a 0.13 um device, and an exponent below 1 makes the limiter
    ## soft enough to bite at drain biases far below saturation -- it
    ## cost 14% on the Vd = 0.05 sweeps before the clamp went in.
    axp = max(_g(card, 'axo', 18.0) / (1.0 + _g(card, 'axl') * iLE), 2.0)
    vp = _g(card, 'vpo', 0.05)

    ## INTERFACE STATES (:267).  PSP does not normalise by the thermal
    ## voltage but by an EFFECTIVE one,
    ## `phit1 = phit * (1 + CT*dCTG) * (1 + dphit1)`
    ## (`PSP103_macrodefs.include:503`), and everything downstream --
    ## the gate drive, the quasi-Fermi levels, the charges -- is in units
    ## of that.  `CTG` is zero on this card so `dCTG` is 1, and `PSCE` is
    ## absent so `dphit1` is 0, leaving `phit1 = phit*(1 + CT)`.
    ##
    ## Worth ~7% on the long device.  The body factor is NOT rescaled:
    ## `G_0` is built on the plain `phit`, and PSP only moves it to
    ## `phit1` when `SWFIX` is set, which defaults to 0 and this card
    ## leaves unset.
    ct = ((_g(card, 'cto') + _g(card, 'ctl') * iLE ** _g(card, 'ctlexp', 1.0))
          * (1.0 + _g(card, 'ctw') * iWE)
          * (1.0 + _g(card, 'ctlw') * iLE * iWE))

    ## Intrinsic series resistance (:304).  PSP applies it as a mobility
    ## term rather than a network element; see `psp_kernel`.
    rs = _g(card, 'rsw1') * iWE * (1.0 + _g(card, 'rsw2') * iWE)

    ## ---- TEMPERATURE ------------------------------------------------
    ##
    ## PSP specifies every parameter at the card's reference temperature
    ## `TR` and scales it to the simulation temperature
    ## (`PSP103_macrodefs.include:291-294, 357-390`).  With
    ## `rTn = TKR/TKD` and `delT = TKD - TKR`, almost everything is a
    ## power law -- `X_T = X * exp(ST_X * ln(rTn))`, which is
    ## `X * (TKR/TKD)^ST_X` -- and the flat-band voltage is the one
    ## exception, a quadratic in `delT`.
    ##
    ## The card carries 28 non-zero `ST*` coefficients, so this is not a
    ## refinement: without it the model is a 27 C model wearing whatever
    ## temperature the caller thought it asked for.  `STVFB`, `STBET`
    ## and `STTHESAT` are themselves geometry-scaled (`:231, 290, 312`);
    ## the rest are plain (`:270-312`).
    ##
    ## NOTE THE SIGN CONVENTION.  `rTn` is REFERENCE over DEVICE, so a
    ## positive `ST` makes the parameter FALL as the device heats -- and
    ## impact ionisation's `A2` alone takes `-STA2` (`:389`).  Easy to
    ## get backwards, and it would read as a plausible temperature
    ## coefficient of the wrong sign rather than as a bug.
    ##
    ## The caller must pass the SAME temperature the simulator will use:
    ## this scales the parameters, and the element's own `vt()` follows
    ## `epar`.  They are not linked, and disagreeing about the
    ## temperature is a quiet way to be wrong.
    tkr = 273.15 + _g(card, 'tr', 27.0)
    r_tn = tkr / T
    ln_rtn = math.log(r_tn)
    del_t = T - tkr

    def _tf(name):
        return math.exp(_g(card, 'st' + name + 'o') * ln_rtn)

    ## The three carrying geometry terms.
    st_vfb = (_g(card, 'stvfbo') + _g(card, 'stvfbl') * iLE
              + _g(card, 'stvfbw') * iWE
              + _g(card, 'stvfblw') * iLE * iWE)
    st_bet = (_g(card, 'stbeto') + _g(card, 'stbetl') * iLE
              + _g(card, 'stbetw') * iWE
              + _g(card, 'stbetlw') * iLE * iWE)
    st_thesat = (_g(card, 'stthesato') + _g(card, 'stthesatl') * iLE
                 + _g(card, 'stthesatw') * iWE
                 + _g(card, 'stthesatlw') * iLE * iWE)

    vfb = (vfb + st_vfb * del_t * (1.0 + _g(card, 'st2vfbo') * del_t)
           + _g(card, 'delvto'))
    ct = ct * _tf('ct')
    u0_eff = u0_eff * math.exp(st_bet * ln_rtn)
    themu = _g(card, 'themuo') * _tf('themu')
    mue = mue * _tf('mue')
    thecs = _g(card, 'thecso') * _tf('thecs')
    cs = cs * _tf('cs')
    xcor = xcor * _tf('xcor')
    rs = rs * _tf('rs')
    thesat = thesat * math.exp(st_thesat * ln_rtn)

    kw = dict(
        w=w, l=l, tox=tox, nsub=neff, vfb=vfb, phib=phib, u0=u0_eff,
        rs=max(rs, 0.0), rsg=_g(card, 'rsgo'), rsb=_g(card, 'rsbo'),
        ct=max(ct, 0.0), alp=max(alp, 0.0), vp=max(vp, 1.0e-6),
        alp1=max(alp1, 0.0), alp2=max(alp2, 0.0), ax=axp,
        kp=max(kp, 0.0),
        mue=max(mue, 0.0), themu=max(themu, 0.0),
        cs=max(cs, 0.0), thecs=max(thecs, 0.0),
        feta=_g(card, 'fetao', 1.0), thesat=max(thesat, 0.0),
    )
    kw.update(dnsub=dnsub, vnsub=vnsub, nslp=nslp, xcor=xcor, rg=rg,
              nfa=nfa, nfb=nfb, nfc=nfc, ef=ef, fnt=fnt, swign=swign,
              cgov=cgov, cfr=cfr, cgbov=cgbov, gov=gov, gov2=gov2,
              ov_a=sp_a, ov_d1=sp_delta1, ov_eps2=sp_eps * sp_eps,
              wcv=wecv, lcv=lecv, qq=qq)
    if all_terms:
        kw.update(thesatb=thesatb, thesatg=thesatg,
                  cf=cf, cfb=cfb, cfd=cfd)
    return kw


#: JUNCAP2's own reference temperature parameter.  It is NOT the MOSFET's
#: `TR`: this card says 21 C where the transistor says 27 C, and assuming
#: they were the same makes every junction capacitance about 2% low --
#: which reads as a geometry error rather than as a temperature one.
JUNCAP_VBILOW = 5.0e-2
JUNCAP_A = 2.0
JUNCAP_EPSCH = 0.1
JUNCAP_AERFC = 0.29214664


def junction_geometry(w, we, ng=1, z1=0.34e-6, z2=0.38e-6, wmin=0.15e-6):
    """Bottom area and the two edge lengths of one junction.

    `SWJUNCAP = 3` does not select components -- all three run -- it
    selects where the geometry comes from (`PSP103_module.include:
    867-883`): the drawn area and perimeter, with the GATE EDGE CARVED
    OUT of the perimeter so it is not counted twice.

    ``LG`` is the ELECTRICAL width `WE`, not the drawn `W`.  That is
    worth stating because it is invisible until you compare two
    geometries: `WE = W - 2*WOT` and `WOT` is negative on this card, so
    the two differ by 0.02 um -- 2% on a 1 um device and 0.2% on a 10 um
    one, which looks like a scaling error rather than like a definition.

    `z1`, `z2` and `wmin` are the PDK SUBCIRCUIT's layout constants
    (`sg13g2_moslv_mod.lib:67`), not model-card parameters, which is why
    they are arguments here rather than read from the card.
    """
    wf = max(w / ng, wmin)
    half = (ng - 1) / 2.0
    area = wf * (z1 + half * z2)
    perim = 2.0 * (wf * (half + 1.0) + z1 + half * z2)
    return dict(ab=area / ng, ls=max(perim / ng - we, 0.0), lg=we)


def junction(card, ab, ls, lg, T=300.0):
    """JUNCAP2's bias-INDEPENDENT constants, per component.

    Everything here is fixed once the card, the geometry and the
    temperature are known, so it is Python rather than symbols and the
    compiled expression carries none of it.

    Only the source-side parameter set exists: `SWJUNASYM = 0` on this
    card aliases every drain-side parameter to it
    (`JUNCAP200_InitModel.include:41-84`), so the drain differs from the
    source only through the geometry -- and for `ng = 1` not even that.

    The temperature reference is JUNCAP's OWN `TRJ`, 21 C here against
    the transistor's 27 C.
    """
    tkr = 273.15 + _g(card, 'trj', 21.0)
    tkd = max(T, 273.15 - 250.0)
    auxt = tkd / tkr
    phitr = PSP_KBOL * tkr / PSP_QELE
    phitd = PSP_KBOL * tkd / PSP_QELE
    d_phigr = -(7.02e-4 * tkr * tkr) / (1108.0 + tkr)
    d_phigd = -(7.02e-4 * tkd * tkd) / (1108.0 + tkd)

    out = {}
    geom = dict(bot=ab, sti=ls, gat=lg)
    for c in ('bot', 'sti', 'gat'):
        phig = _g(card, 'phig' + c, 1.16)
        ## `ftd` is n_i(T)/n_i(Tref); everything else keys off it.
        ftd = (auxt ** 1.5
               * math.exp(0.5 * ((phig + d_phigr) / phitr
                                 - (phig + d_phigd) / phitd)))
        vbir = _g(card, 'vbir' + c, 1.0)
        p = _g(card, 'p' + c, 0.5)
        ubi = vbir * auxt - 2.0 * phitd * math.log(ftd)
        ## A soft floor at 50 mV, written so the exponential cannot
        ## overflow when `ubi` is far below it.
        z = (JUNCAP_VBILOW - ubi) / phitd
        vbi = ubi + phitd * (z + math.log1p(math.exp(-z)) if z > 0.0
                             else math.log1p(math.exp(z)))
        cjo = _g(card, 'cjor' + c, 1.0e-3) * (vbir / vbi) ** p
        out[c] = dict(
            p=p, vbi=vbi, cjo=cjo,
            qpref=cjo * vbi / (1.0 - p), qpref2=JUNCAP_A * cjo,
            idsat=_g(card, 'idsatr' + c, 1.0e-12) * ftd * ftd,
            ftd=ftd, area=geom[c])
        ## ---- reverse-leakage constants -------------------------------
        ## `wdepnulr` uses the REFERENCE `CJOR`, not the temperature-
        ## scaled `cjo` (`JUNCAP200_InitModel.include:205-207`), and the
        ## two sidewall components carry a junction depth where the
        ## bottom one does not.
        xjun = 1.0 if c == 'bot' else _g(card, 'xjun' + c, 1.0e-7)
        cjor = _g(card, 'cjor' + c, 1.0e-3)
        wdepnulr = xjun * PSP_EPS_SI / cjor
        ## Half the bandgap, floored at the thermal voltage (`:225-227`).
        delta_e = max(0.5 * (phig + d_phigd), phitd)
        out[c].update(
            csrh=_g(card, 'csrh' + c, 100.0 if c == 'bot' else 1.0e-4),
            ctat=_g(card, 'ctat' + c, 100.0 if c == 'bot' else 1.0e-4),
            cbbt=_g(card, 'cbbt' + c,
                    1.0e-12 if c == 'bot' else 1.0e-18),
            vbir=vbir, vbirinv=1.0 / vbir,
            wdepnulr=wdepnulr, wdepnulrinv=1.0 / wdepnulr,
            omp=1.0 - p, oomp=1.0 / (1.0 - p),
            atat=delta_e / phitd,
            btatpart=math.sqrt(32.0 * _g(card, 'mefftat' + c, 0.25)
                               * PSP_MELE * PSP_QELE * delta_e ** 3)
            / (3.0 * PSP_HBAR),
            ## Linear in temperature, not Arrhenius (`:240-245`).
            fbbt=max(_g(card, 'fbbtr' + c, 1.0e9)
                     * (1.0 + _g(card, 'stfbbt' + c, -1.0e-3)
                        * (tkd - tkr)), 0.0),
            vbr=_g(card, 'vbr' + c, 10.0), pbr=_g(card, 'pbr' + c, 4.0))

    ## Instance-level constants.  `vfmin` is what guarantees the grading
    ## power never sees a non-positive base: it holds `1 - vj/vbi` above
    ## `2^(-1/pmax)` for EVERY component, which is 0.28 here.
    vbimin = min(out[c]['vbi'] for c in out)
    pmax = max(out[c]['p'] for c in out)
    vfmin = vbimin * (1.0 - JUNCAP_A ** (-1.0 / pmax))

    ## `VMAX` caps the ideal diode's exponent; a component with zero
    ## saturation current does not participate (`macrodefs:156-171`),
    ## and `idsatrgat` IS zero on this card.
    imax = _g(card, 'imax', 1.0e-3)
    vmax = []
    for c in out:
        js = out[c]['idsat'] * out[c]['area']
        vmax.append(phitd * math.log(imax / js + 1.0) if js > 0.0
                    else 1.0e8)
    ## `vbbtlim` keeps the band-to-band term's `VBIR - vbbt` positive
    ## (`macrodefs:204`); `alphaav` and the breakdown constants come
    ## from `FREV` (`InitModel:248-261`).
    alphaav = 1.0 - 1.0 / _g(card, 'frev', 1000.0)
    perfc = math.sqrt(math.pi) * JUNCAP_AERFC
    berfc = (-5.0 * JUNCAP_AERFC + 6.0 - perfc ** -2.0) / 3.0
    extra = dict(jphitr=phitr, jalphaav=alphaav,
                 jperfc=perfc, jberfc=berfc,
                 jcerfc=1.0 - JUNCAP_AERFC - berfc,
                 jvbbtlim=min(out[c]['vbir'] for c in out) - 5.0e-2)
    for c in out:
        vbr, pbr = out[c]['vbr'], out[c]['pbr']
        fstop = 1.0 / (1.0 - alphaav ** pbr)
        out[c]['fstop'] = fstop
        out[c]['vbrinv'] = 1.0 / vbr
        out[c]['slope'] = -(fstop * fstop * alphaav ** (pbr - 1.0)
                            * pbr / vbr)

    return dict(
        jphit=phitd, jvfmin=vfmin, jvch=vbimin * JUNCAP_EPSCH,
        jvbimin=vbimin, **extra,
        jvmax=min(vmax),
        ## The ideal term factors: one exponential times a total.
        jisat=sum(out[c]['idsat'] * out[c]['area'] for c in out),
        **{'j%s_%s' % (c, k): out[c][k]
           for c in out
           for k in ('p', 'vbi', 'qpref', 'qpref2', 'area', 'idsat',
                     'csrh', 'ctat', 'cbbt', 'vbir', 'vbirinv',
                     'wdepnulr', 'wdepnulrinv', 'omp', 'oomp', 'atat',
                     'btatpart', 'fbbt', 'vbr', 'pbr', 'fstop',
                     'vbrinv', 'slope', 'ftd')})
