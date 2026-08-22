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


def to_long_channel(card, w, l, T=300.0):
    """Card + geometry -> `PspMosLongChannel` keyword arguments.

    `card` is what `spicecard.Deck.model_params` returned.  Returns a
    dict ready to splat into the element.
    """
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
        qq = 0.4 * QMN * qmc * cox_prime ** (2.0 / 3.0)
        qb0 = math.sqrt(phit * gamma0 * gamma0 * phib)
        dphibq = 0.75 * qq * qb0 ** (2.0 / 3.0)
        phib = phib + dphibq
        g_fac = 1.0 + 2.0 * (2.0 / 3.0) * dphibq / qb0
        neff = neff * g_fac * g_fac

    ## Gain factor (:286-289).  `GPE` and `GWE` are the length and width
    ## corrections to mobility; `BETN = UO*WE/(GPE*LE)*GWE`, and
    ## `BET = BETN * Cox'` (`PSP103_macrodefs.include:368`).
    gpe = max(1.0 + _g(card, 'fbet1') * _g(card, 'lp1', 1.0) / LE
              * (1.0 - math.exp(-LE / _g(card, 'lp1', 1.0)))
              + _g(card, 'fbet2') * _g(card, 'lp2', 1.0) / LE
              * (1.0 - math.exp(-LE / _g(card, 'lp2', 1.0))), 1.0e-15)
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

    ## Intrinsic series resistance (:304).  PSP applies it as a mobility
    ## term rather than a network element; see `psp_kernel`.
    rs = _g(card, 'rsw1') * iWE * (1.0 + _g(card, 'rsw2') * iWE)

    return dict(
        w=w, l=l, tox=tox, nsub=neff, vfb=vfb, phib=phib, u0=u0_eff,
        rs=max(rs, 0.0), rsg=_g(card, 'rsgo'), rsb=_g(card, 'rsbo'),
        mue=max(mue, 0.0), themu=max(_g(card, 'themuo'), 0.0),
        cs=max(cs, 0.0), thecs=max(_g(card, 'thecso'), 0.0),
        feta=_g(card, 'fetao', 1.0), thesat=max(thesat, 0.0),
    )
