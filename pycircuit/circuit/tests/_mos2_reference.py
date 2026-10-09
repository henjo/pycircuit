"""An independent numpy transcription of ngspice-47's MOS level 2 channel
current (`mos2temp.c`, `mos2load.c` lines 556-1010) -- the oracle the HDL
`MosLevel2Hdl` is tested against (the SPICE benchmark plan's stage 9).

Scalar and branchy, line for line with the C, ON PURPOSE: it shares no
code with the model it checks.  Only the channel current ``cdrain`` (in
the n-channel convention, drain to source, signed by the mode) is
transcribed -- the junctions, overlaps and resistors are level 3's, tested
there.  ngspice's hand derivatives (gm, gds, gmbs) are NOT transcribed:
the model's are exact derivatives of this current, checked by finite
differences.
"""
import math

#: ngspice-47's constants (`const.h`, `mos2temp.c`); `use_constants`
#: swaps in a model's own (pycircuit's are SI-2019, its eps0 CODATA-1986:
#: 3e-6 apart from ngspice's -- the floor of any comparison of the two)
NGSPICE = {'CHARGE': 1.6021766208e-19, 'CONSTboltz': 1.38064852e-23,
           'EPS0': 8.854214871e-12}
CHARGE = CONSTboltz = CONSTKoverQ = EPSSIL = EPS0 = None
REFTEMP = 300.15


def use_constants(CHARGE, CONSTboltz, EPS0):
    """Set the physical constants every function here reads."""
    g = globals()
    g.update(CHARGE=CHARGE, CONSTboltz=CONSTboltz, EPS0=EPS0,
             CONSTKoverQ=CONSTboltz / CHARGE, EPSSIL=11.7 * EPS0)


use_constants(**NGSPICE)

DEFAULTS = {'vto': 0.0, 'kp': None, 'u0': 600.0, 'gamma': 0.0, 'phi': 0.6, 'tox': 1e-7,
            'nsub': None, 'nss': 0.0, 'tpg': 1.0, 'xj': 0.0, 'ld': 0.0, 'lambd': 0.0,
            'delta': 0.0, 'ucrit': 1e4, 'uexp': 0.0, 'vmax': 0.0, 'neff': 1.0, 'nfs': 0.0,
            'w': 1e-4, 'l': 1e-4, 'tnom': 27.0, 'pb': 0.8, 'temp': 27.0, 'type': 1}


def card(**kw):
    """A level-2 card with ngspice's defaults filled in; `vto` given only
    when the caller gives it (it is computed from `nsub` otherwise)."""
    c = dict(DEFAULTS)
    c['vto_given'] = 'vto' in kw
    c['gamma_given'] = 'gamma' in kw
    c['phi_given'] = 'phi' in kw
    c['kp_given'] = 'kp' in kw and kw['kp'] is not None
    c.update(kw)
    return c


def _temp(c):
    """`mos2temp.c`: the model preprocessing and the instance's
    temperature-dependent values."""
    tnom = c['tnom'] + 273.15
    T = c['temp'] + 273.15
    typ = c['type']
    fact1 = tnom / REFTEMP
    vtnom = tnom * CONSTKoverQ
    kt1 = CONSTboltz * tnom
    egfet1 = 1.16 - (7.02e-4 * tnom * tnom) / (tnom + 1108)
    arg1 = -egfet1 / (kt1 + kt1) + 1.1150877 / (CONSTboltz * (REFTEMP + REFTEMP))
    pbfact1 = -2 * vtnom * (1.5 * math.log(fact1) + CHARGE * arg1)
    phi, gamma, vt0 = c['phi'], c['gamma'], c['vto']
    cox = 3.9 * EPS0 / c['tox']
    kp = c['kp'] if c['kp_given'] else c['u0'] * 1e-4 * cox
    xd, nsub = 0.0, c['nsub'] or 0.0
    if c['nsub'] is not None:
        if nsub * 1e6 > 1.45e16:
            if not c['phi_given']:
                phi = max(0.1, 2 * vtnom * math.log(nsub * 1e6 / 1.45e16))
            fermis = typ * 0.5 * phi
            wkfng = 3.2
            tpg = c['tpg']
            if tpg != 0:
                fermig = typ * tpg * 0.5 * egfet1
                wkfng = 3.25 + 0.5 * egfet1 - fermig
            wkfngs = wkfng - (3.25 + 0.5 * egfet1 + fermis)
            if not c['gamma_given']:
                gamma = math.sqrt(2 * 11.70 * EPS0 * CHARGE * nsub * 1e6) / cox
            if not c['vto_given']:
                vfb = wkfngs - c['nss'] * 1e4 * CHARGE / cox
                vt0 = vfb + typ * (gamma * math.sqrt(phi) + phi)
            xd = math.sqrt((EPSSIL + EPSSIL) / (CHARGE * nsub * 1e6))
        else:
            raise ValueError('Nsub < Ni')
    vt = T * CONSTKoverQ
    ratio = T / tnom
    fact2 = T / REFTEMP
    kt = T * CONSTboltz
    egfet = 1.16 - (7.02e-4 * T * T) / (T + 1108)
    arg = -egfet / (kt + kt) + 1.1150877 / (CONSTboltz * (REFTEMP + REFTEMP))
    pbfact = -2 * vt * (1.5 * math.log(fact2) + CHARGE * arg)
    ratio4 = ratio * math.sqrt(ratio)
    tkp = kp / ratio4
    phio = (phi - pbfact1) / fact1
    tphi = fact2 * phio + pbfact
    tvbi = (vt0 - typ * (gamma * math.sqrt(phi)) + 0.5 * (egfet1 - egfet)
            + typ * 0.5 * (tphi - phi))
    pbo = (c['pb'] - pbfact1) / fact1          # (pb's only use here: tBulkPot)
    tbulkpot = fact2 * pbo + pbfact
    return {'vt': vt, 'cox': cox, 'tkp': tkp, 'tphi': tphi, 'tvbi': tvbi, 'gamma': gamma,
            'xd': xd, 'nsub': nsub, 'tbulkpot': tbulkpot, 'phi': phi}


def _quartic(a1, b1, c1, d1):
    """`mos2load.c` lines 770-824: Ferrari through the resolvent cubic,
    the smallest positive root with |P| <= 1e-6, or None."""
    sig1 = (1.0, -1.0, 1.0, -1.0)
    sig2 = (1.0, 1.0, -1.0, -1.0)
    a = -b1
    b = a1 * c1 - 4.0 * d1
    c = -d1 * (a1 * a1 - 4.0 * b1) - c1 * c1
    r = -a * a / 3.0 + b
    s = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c
    r3 = r * r * r
    s2 = s * s
    p = s2 / 4.0 + r3 / 27.0
    p0 = abs(p)
    p2 = math.sqrt(p0)
    if p < 0:
        ro = math.sqrt(s2 / 4.0 + p0)
        ro = math.log(ro) / 3.0
        ro = math.exp(ro)
        fi = math.atan(-2.0 * p2 / s)
        y3 = 2.0 * ro * math.cos(fi / 3.0) - a / 3.0
    else:
        p3 = (-s / 2.0 + p2)
        p3 = math.exp(math.log(abs(p3)) / 3.0)
        p4 = (-s / 2.0 - p2)
        p4 = math.exp(math.log(abs(p4)) / 3.0)
        y3 = p3 + p4 - a / 3.0
    a3 = math.sqrt(a1 * a1 / 4.0 - b1 + y3)
    b3 = math.sqrt(y3 * y3 / 4.0 - d1)
    x4 = []
    for i in range(4):
        a4 = a1 / 2.0 + sig1[i] * a3
        b4 = y3 / 2.0 + sig2[i] * b3
        delta4 = a4 * a4 / 4.0 - b4
        if delta4 < 0:
            continue
        tmp = math.sqrt(delta4)
        x4.append(-a4 / 2.0 + tmp)
        x4.append(-a4 / 2.0 - tmp)
    xvalid, jknt = None, 0
    for x in x4:
        if x <= 0:
            continue
        poly4 = x * x * x * x + a1 * x * x * x
        poly4 = poly4 + b1 * x * x + c1 * x + d1
        if abs(poly4) > 1.0e-6:
            continue
        jknt += 1
        if jknt <= 1:
            xvalid = x
        if x > xvalid:
            continue
        xvalid = x
    return xvalid


def cdrain(c, vgs, vds, vbs, spice_band=True):
    """The channel current, drain to source in the n-channel convention
    (`mode * cdrain`), at terminal voltages already in that convention."""
    t = _temp(c)
    vt, tphi = t['vt'], t['tphi']
    w, L = c['w'], c['l']
    eff_l = L - 2 * c['ld']
    beta = t['tkp'] * w / eff_l
    unscaled_cox = t['cox'] * eff_l * w
    oxide_cap = unscaled_cox
    nfs, gamma, xd, xj = c['nfs'], t['gamma'], t['xd'], c['xj']
    vbd = vbs - vds
    vgd = vgs - vds
    mode = 1 if vds >= 0 else -1
    lvbs = vbs if mode == 1 else vbd
    lvds = mode * vds
    lvgs = vgs if mode == 1 else vgd
    phi_min_vbs = tphi - lvbs
    xlamda = c['lambd']
    if lvbs <= 0.0:
        sarg1 = math.sqrt(phi_min_vbs)
    else:
        sphi = math.sqrt(tphi)
        sarg1 = sphi / (1.0 + 0.5 * lvbs / tphi)
    if lvbs - lvds <= 0:
        barg = math.sqrt(phi_min_vbs + lvds)
    else:
        sphi = math.sqrt(tphi)
        barg = sphi / (1.0 + 0.5 * (lvbs - lvds) / tphi)
    factor = 0.125 * c['delta'] * 2.0 * math.pi * EPSSIL / unscaled_cox * eff_l
    eta = 1.0 + factor
    vbin = t['tvbi'] * c['type'] + factor * phi_min_vbs
    if gamma > 0.0 or t['nsub'] > 0.0:
        xwd = xd * barg
        xws = xd * sarg1
        argss = argsd = 0.0
        if xj > 0:
            tmp = 2.0 / xj
            args = math.sqrt(1.0 + xws * tmp)
            argd = math.sqrt(1.0 + xwd * tmp)
            tmp = 0.5 * xj / eff_l
            argss = tmp * (args - 1.0)
            argsd = tmp * (argd - 1.0)
        gamasd = gamma * (1.0 - argss - argsd)
    else:
        gamasd = gamma
    von = vbin + gamasd * sarg1
    vdsat = 0.0
    xn = 0.0
    argg = 0.0
    if nfs != 0.0 and oxide_cap != 0.0:
        cfs = CHARGE * nfs * 1e4
        ## (cdonco: the derivative chain's; the VALUE of xn uses it)
        dsrgdb = -0.5 / sarg1 if lvbs <= 0.0 else \
            -0.5 * sarg1 * (sarg1 / (tphi * math.sqrt(tphi)))
        dgddvb = _dgddvb(c, t, lvbs, lvds, sarg1, barg, eff_l)
        cdonco = -(gamasd * dsrgdb + dgddvb * sarg1) + factor
        xn = 1.0 + cfs / unscaled_cox * w * eff_l + cdonco
        tmp = vt * xn
        von = von + tmp
        argg = 1.0 / tmp
    else:
        if lvgs <= vbin:
            return 0.0
    sarg3 = sarg1 * sarg1 * sarg1
    body = barg * barg * barg - sarg3
    ## effective mobility
    if oxide_cap <= 0.0:
        ufact = 1.0
    else:
        udenom = lvgs - von
        tmp = c['ucrit'] * 100 * EPSSIL / t['cox']
        if udenom <= tmp:
            ufact = 1.0
        else:
            ufact = math.exp(c['uexp'] * math.log(tmp / udenom))
    ueff = c['u0'] * 1e-4 * ufact
    ## saturation voltage (Grove-Frohman)
    vgsx = lvgs
    gammad = gamasd / eta
    if nfs != 0 and oxide_cap != 0:
        vgsx = max(lvgs, von)
    if gammad > 0:
        gammd2 = gammad * gammad
        argv = (vgsx - vbin) / eta + phi_min_vbs
        if argv <= 0.0:
            vdsat = 0.0
        else:
            arg1 = math.sqrt(1.0 + 4.0 * argv / gammd2)
            vdsat = (vgsx - vbin) / eta + gammd2 * (1.0 - arg1) / 2.0
            vdsat = max(vdsat, 0.0)
    else:
        vdsat = (vgsx - vbin) / eta
        vdsat = max(vdsat, 0.0)
    if c['vmax'] > 0:
        gammd2 = gammad * gammad
        v1 = (vgsx - vbin) / eta + phi_min_vbs
        v2 = phi_min_vbs
        xv = c['vmax'] * eff_l / ueff
        a1 = gammad / 0.75
        b1 = -2.0 * (v1 + xv)
        c1 = -2.0 * gammad * xv
        d1 = 2.0 * v1 * (v2 + xv) - v2 * v2 - 4.0 / 3.0 * gammad * sarg3
        xvalid = _quartic(a1, b1, c1, d1)
        if xvalid is not None:
            vdsat = xvalid * xvalid - phi_min_vbs
    ## effective channel length
    bsarg = None
    if lvds != 0.0:
        gammad = gamasd
        if lvbs - vdsat <= 0:
            bsarg = math.sqrt(vdsat + phi_min_vbs)
        else:
            sphi = math.sqrt(tphi)
            bsarg = sphi / (1.0 + 0.5 * (lvbs - vdsat) / tphi)
        bodys = bsarg * bsarg * bsarg - sarg3
        if c['vmax'] <= 0:
            if t['nsub'] != 0.0 and not xlamda > 0.0:
                argv = (lvds - vdsat) / 4.0
                sargv = math.sqrt(1.0 + argv * argv)
                arg1 = math.sqrt(argv + sargv)
                xlfact = xd / (eff_l * lvds)
                xlamda = xlfact * arg1
        else:
            xdv = xd / math.sqrt(c['neff'])
            xlv = c['vmax'] * xdv / (2.0 * ueff)
            if t['nsub'] != 0.0 and not xlamda > 0.0:
                argv = max(lvds - vdsat, 0.0)
                xls = math.sqrt(xlv * xlv + argv)
                xlfact = xdv / (eff_l * lvds)
                xlamda = xlfact * (xls - xlv)
    ## punch-through
    sbiarg = math.sqrt(t['tbulkpot'])
    xwb = xd * sbiarg
    xld = eff_l - xwb
    clfact = 1.0 - xlamda * lvds
    xleff = eff_l * clfact
    deltal = xlamda * lvds * eff_l
    if t['nsub'] == 0.0:
        xwb = 0.25e-6
    if xleff < xwb:
        xleff = xwb / (1.0 + (deltal - xld) / xwb)
        clfact = xleff / eff_l
    beta1 = beta * ufact / clfact
    gammad = gamasd
    ## (`lvds <= 1e-10`: SPICE returns 0 here; the model does not -- its
    ## expressions reach 0 with vds and keep SPICE's gds as their slope.
    ## `spice_band=False` follows the model: the comparison it makes)
    if lvds <= 1.0e-10 and spice_band:
        return 0.0
    if nfs != 0 and oxide_cap != 0:
        if lvgs > von:
            return mode * _strong(beta1, lvgs, vbin, eta, lvds, vdsat, gammad, body, bodys)
    else:
        if lvgs > vbin:
            return mode * _strong(beta1, lvgs, vbin, eta, lvds, vdsat, gammad, body, bodys)
        return 0.0
    if lvgs > von:
        return mode * _strong(beta1, lvgs, vbin, eta, lvds, vdsat, gammad, body, bodys)
    ## subthreshold
    if vdsat <= 0:
        return 0.0
    vdson = min(vdsat, lvds)
    if lvds > vdsat:
        body = bodys
    cdson = beta1 * ((von - vbin - eta * vdson * 0.5) * vdson - gammad * body / 1.5)
    expg = math.exp(argg * (lvgs - von))
    return mode * cdson * expg


def _strong(beta1, lvgs, vbin, eta, lvds, vdsat, gammad, body, bodys):
    if lvds <= vdsat:
        return beta1 * ((lvgs - vbin - eta * lvds / 2.0) * lvds - gammad * body / 1.5)
    return beta1 * ((lvgs - vbin - eta * vdsat / 2.0) * vdsat - gammad * bodys / 1.5)


def _dgddvb(c, t, lvbs, lvds, sarg1, barg, eff_l):
    """`dgddvb` (lines 615-650), which the VALUE of `xn` reads through
    `cdonco`: zero unless gamma/nsub and xj are given."""
    gamma, xd, xj, tphi = t['gamma'], t['xd'], c['xj'], t['tphi']
    if not (gamma > 0.0 or t['nsub'] > 0.0) or not xj > 0:
        return 0.0
    if lvbs <= 0.0:
        dsrgdb = -0.5 / sarg1
    else:
        dsrgdb = -0.5 * sarg1 * (sarg1 / (tphi * math.sqrt(tphi)))
    if lvbs - lvds <= 0:
        dbrgdb = -0.5 / barg
    else:
        dbrgdb = -0.5 * barg * (barg / (tphi * math.sqrt(tphi)))
    tmp = 2.0 / xj
    args = math.sqrt(1.0 + xd * sarg1 * tmp)
    argd = math.sqrt(1.0 + xd * barg * tmp)
    tmp = 0.5 / eff_l
    dbargs = tmp * xd * dsrgdb / args
    dbargd = tmp * xd * dbrgdb / argd
    return -gamma * (dbargs + dbargd)
