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

from pycircuit.circuit.constants import (eps0, epsRSi, epsRSiO2,
                                         kboltzmann, qelectron)


#: PSP's normalisation lengths (`PSP103_macrodefs.include:58-59`).
LEN = 1.0e-6
WEN = 1.0e-6

#: Quantum-mechanical correction constants
#: (`PSP103_macrodefs.include:51-52`), n- and p-channel.
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
    return dict(LE=LE, WE=WE, iLE=LEN / LE, iWE=WEN / WE)


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


def to_long_channel(card, w, l, T=300.0, all_terms=False):
    """Card + geometry -> `PspMosLongChannel` keyword arguments.

    `card` is what `spicecard.Deck.model_params` returned.  Returns a
    dict ready to splat into the element.

    `all_terms=True` additionally wires up two blocks that are
    implemented, verified term by term against PSP's own scaled
    parameters, and CURRENTLY MAKE THE MEASURED FIT WORSE:

    * **DIBL** -- `CF`, `CFB`, `CFD` (`macrodefs:473-476`);
    * the **body- and gate-bias modulation of `THESAT`** -- `THESATB`,
      `THESATG` (`macrodefs:596-607`).

    Summed median error over the six n-channel sweeps: 0.016 without
    them, 0.034 with the `THESAT` modulation alone, 0.041 with both.

    They are off by default and kept in the tree deliberately, following
    the precedent this model set with channel-length modulation -- which
    was also correct, also measured worse, was left out with the
    reasoning recorded, and turned out to be the single largest accuracy
    gain the model ever took once the term it was compensating for
    arrived.  Discarding correct physics because it exposes an error
    elsewhere is how that error gets preserved.

    What the evidence says about where the error is NOT: all twenty-one
    scaled parameters the reference records -- including `cf`, `cfb`,
    `thesat`, `thesatb`, `thesatg`, `alp`, `alp1` and `alp2` -- match
    PSP's own `lp_*` outputs exactly at four geometries.  So the scaling
    is right and the discrepancy is in a formula.

    And where it probably IS: our `delVg` comes to 3.6 mV at Vds = 1.35,
    which is precisely the 3.5 mV shift measured out of PSP's own `vth`.
    In weak inversion, at this device's 85 mV/decade, 3.6 mV is a **9%**
    current change -- so DIBL is a real part of the 2.4x climb at
    Vg = 0.6 that was attributed wholly to `FdL`.  `FdL` was accepted on
    a residual that already contained this omission, and the near-
    threshold sweep is exactly where enabling DIBL now overshoots
    (1.003 -> 1.023).  That is the term to re-examine, with these on.
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
    eps_ox = _g(card, 'epsroxo', epsRSiO2) * eps0
    cox_prime = eps_ox / tox

    ## Flat-band voltage (:230) and effective doping (:250).
    vfb = (_g(card, 'vfbo') + _g(card, 'vfbl') * iLE
           + _g(card, 'vfbw') * iWE + _g(card, 'vfblw') * iLE * iWE)
    neff = _neff(card, geo)

    ## Surface potential at threshold (`PSP103_macrodefs.include:295-308`).
    ## `phibFac` and `Eg` are the temperature-dependent parts; at 300 K
    ## they are what they are, and no temperature COEFFICIENT from the
    ## card is applied here.
    phit = kboltzmann * T / qelectron
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
    gamma0 = math.sqrt(2.0 * qelectron * neff * epsRSi * eps0
                       / phit) / cox_prime
    qmc = _g(card, 'qmc')
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
    np_card = _g(card, 'npo')
    if np_card > 0.0:
        np_eff = max(max(np_card, 8.0e7 / (tox * tox)), 5.0e24)
        kp = (2.0 * cox_prime * cox_prime * phit
              / (qelectron * np_eff * epsRSi * eps0))
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

    kw = dict(
        w=w, l=l, tox=tox, nsub=neff, vfb=vfb, phib=phib, u0=u0_eff,
        rs=max(rs, 0.0), rsg=_g(card, 'rsgo'), rsb=_g(card, 'rsbo'),
        ct=max(ct, 0.0), alp=max(alp, 0.0), vp=max(vp, 1.0e-6),
        alp1=max(alp1, 0.0), alp2=max(alp2, 0.0), ax=axp,
        kp=max(kp, 0.0),
        mue=max(mue, 0.0), themu=max(_g(card, 'themuo'), 0.0),
        cs=max(cs, 0.0), thecs=max(_g(card, 'thecso'), 0.0),
        feta=_g(card, 'fetao', 1.0), thesat=max(thesat, 0.0),
    )
    if all_terms:
        kw.update(thesatb=thesatb, thesatg=thesatg,
                  cf=cf, cfb=cfb, cfd=cfd)
    return kw
