"""JUNCAP2's junction charge, and the ideal diode that goes with it.

PSP103 does not instantiate a separate junction model -- it includes
JUNCAP200 into its own module -- but the two are different compact
models, and this keeps them apart.

WHAT IS HERE AND WHY.  The junction CAPACITANCE is 8% of the intrinsic
gate capacitance on a 10 um device and 126% of it at 0.13 um, so it is
most of what a transient or AC simulation needs.  The junction CURRENT
is around 1e-15 A against a 4e-4 A drain current -- eleven orders down,
contributing nothing to an operating point and nothing to `gm` or
`gds`.  So the charge is here in full, and of the current only the ideal
diode term, which costs three lines and is the one part that matters in
the one regime where junction current does: a forward-biased bulk.

Left out deliberately, with the vendor lines so it can be added: the
Shockley-Read-Hall recombination term (`JUNCAP200_macrodefs.include:
269-284`), trap-assisted tunnelling (`:285-303`, fifteen of the
twenty-nine lines per component and the whole reason for the `erfc`
helper), band-to-band tunnelling (`:304-311`), and avalanche breakdown
(`:312-321`).  Together they shape a reverse leakage of 2.6e-15 A at
-3 V; the ideal term alone gives 1.7e-19 A there, which is a different
number for a quantity that is eleven orders below anything else in the
model.

The bias-independent constants come from `psp_scaling.junction`, in
Python, so none of the temperature or geometry arithmetic reaches the
compiled expression.
"""

import sympy

from pycircuit.circuit import hdl


COMPONENTS = ('bot', 'sti', 'gat')

#: Junction voltages are clamped to this before anything else is done
#: with them.  A kilovolt is a hundred times the card's breakdown and
#: a thousand times its supply, so the clamp never acts on a bias any
#: real circuit reaches -- but a solver that overshoots does reach
#: absurd values, and every expression below has an exponential or a
#: fractional power in it.  Bounding the INPUT once is cheaper and
#: easier to reason about than guarding each of them.
VJUN_CLAMP = 1.0e3
JUNCAP_AERFC = 0.29214664
JUNCAP_SQRTPI = 1.77245385090551603
JUNCAP_TWOTHIRDS = 0.666666666666667


def _clamped(v, tag):
    return _v(sympy.Max(sympy.Min(v, VJUN_CLAMP), -VJUN_CLAMP),
              'jvc' + tag)


def _v(expr, name):
    return hdl.var(expr, name)


def clamp_forward(v, vfmin, vch, tag=''):
    """PSP's `hypfunction5` (`JUNCAP200_macrodefs.include:236-241`): a
    smooth clamp of `v` to `vfmin` from ABOVE, leaving `v` untouched
    below it.

    This is the guard the whole charge model rests on.  The grading
    power below raises `1 - v/vbi` to a fractional exponent, and it is
    this clamp that keeps the base positive -- `vfmin` is defined
    (`macrodefs:203`) as exactly the voltage at which the base reaches
    `2^(-1/pmax)`, which is 0.28 on this card.  Remove the clamp and the
    base goes through zero into negative territory at forward bias.

    Branchless, and its denominator is bounded below by `2*vfmin`, so
    nothing here needs guarding.
    """
    h2d = _v(v + vch * vch / vfmin, 'jh2d' + tag)
    h5 = _v(sympy.sqrt((vfmin - h2d) ** 2 + 4.0 * vch * vch), 'jh5' + tag)
    return _v(2.0 * (v * vfmin) * hdl.safe_div(1.0, vfmin + h2d + h5,
                                               eps=1e-30),
              'jvj' + tag)


def charge(v, par, tag=''):
    """The junction charge at junction voltage `v`, summed over the
    three components (`macrodefs:266-267`).

    Per component, ``Q' = qpref*(1 - (1 - vj/vbi)^(1-P)) + qpref2*(V - vj)``.

    The first term is the depletion charge, whose derivative is the
    textbook ``cjo*(1 - V/vbi)^(-P)``.  The second is what replaces the
    classical forward-bias linearisation: above `vfmin` the clamped `vj`
    stops moving, so `(V - vj)` grows linearly and pins the forward
    capacitance at `2*cjo`.  The `2` in `qpref2` and the `2` inside
    `vfmin` are the same constant, and that is not a coincidence -- it
    is what makes the capacitance continuous across the join.
    """
    v = _clamped(v, tag)
    vj = clamp_forward(v, par['jvfmin'], par['jvch'], tag)
    q = 0
    for c in COMPONENTS:
        p, vbi = par['j%s_p' % c], par['j%s_vbi' % c]
        ## The base is bounded below by `2^(-1/pmax)` by construction,
        ## so the fractional power needs no guard -- see `clamp_forward`.
        base = _v(1.0 - vj * hdl.safe_div(1.0, vbi, eps=1e-30),
                  'jbase_%s%s' % (c, tag))
        q = q + par['j%s_area' % c] * (
            par['j%s_qpref' % c] * (1.0 - base ** (1.0 - p))
            + par['j%s_qpref2' % c] * (v - vj))
    return _v(q, 'jq' + tag)


def current(v, par, tag=''):
    """The ideal diode term only (`macrodefs:268`, `:334-341`).

    `idsat*(exp(V/phit) - 1)`, with PSP's linear extrapolation above
    `VMAX` so the forward current's slope saturates rather than the
    exponential running away.  The saturation currents of the three
    components factor out of the exponential, so this is one
    exponential and one constant.

    BOTH ARMS ARE FINITE, which they are not in the vendor as written:
    its high-side arm takes `sqrt` of `1 + (V - VMAX)/phit`, and that
    quantity is -1e18 at zero bias.  PSP never evaluates it there; a DSL
    that evaluates both arms of every conditional does.  Clamping the
    excess at zero makes the unused arm harmless and changes nothing in
    the used one.
    """
    v = _clamped(v, 'i' + tag)
    phit, vmax = par['jphit'], par['jvmax']
    inv_phit = _v(hdl.safe_div(1.0, phit, eps=1e-30), 'jiphit' + tag)
    ## `exp(VMAX/phit)` is parameter-only, so it costs one evaluation
    ## per call and nothing per bias point; `expl` rather than `exp`
    ## because a card with a large `VMAX` would otherwise overflow it.
    e_vmax = _v(hdl.expl(vmax * inv_phit), 'jevmax' + tag)
    ## The exponential arm's input is CLAMPED TO ITS OWN DOMAIN, not
    ## merely selected against.  Above `VMAX` this arm is discarded --
    ## but it is still evaluated, and `expl` continues polynomially
    ## above its threshold, so at an absurd bias the cube overflows and
    ## the discarded arm poisons the derivative.  Below `VMAX` the
    ## clamp is the identity, so it costs nothing where it is used.
    lo = _v(hdl.expl(sympy.Min(v, vmax) * inv_phit), 'jelo' + tag)
    hi = _v((1.0 + sympy.Max(v - vmax, 0.0) * inv_phit) * e_vmax,
            'jehi' + tag)
    idmult = _v(sympy.Piecewise((lo, v < vmax), (hi, True)) - 1.0,
                'jidm' + tag)
    return _v(par['jisat'] * idmult, 'ji' + tag)


def _hyp2(x, x0, eps, name):
    """PSP's `hypfunction2` -- a smooth minimum with margin `eps`
    (`JUNCAP200_macrodefs.include:231`)."""
    return _v(0.5 * (x + x0 - sympy.sqrt((x - x0) ** 2 + 4.0 * eps * eps)),
              name)


def _erfc_exp(y, m, par, tag):
    """`erfc(y)*exp(m)`, PSP's rational approximation
    (`macrodefs:247-261`), written branchlessly.

    The vendor splits on the sign of `y` into `1/(1 + p*y)` and
    `1/(1 - p*y)`.  Those are the same function of `|y|`, and writing
    them as one removes a POLE: the unused arm has `1 - p*y` vanishing
    at `|y| = 1/p = 1.93`, and goes negative beyond it into a cube.  In
    a DSL that evaluates both arms that pole is not hypothetical.
    """
    perfc, berfc, cerfc = par['jperfc'], par['jberfc'], par['jcerfc']
    t = _v(hdl.safe_div(1.0, 1.0 + perfc * sympy.Abs(y), eps=1e-30),
           'jterfc' + tag)
    poly = _v((JUNCAP_AERFC * t + berfc * t * t + cerfc * t * t * t),
              'jerfp' + tag)
    ## `expl` both sides: the vendor's `expl_low` clamps only the
    ## negative tail, and `m` alone is unclamped by construction.
    pos = _v(poly * hdl.expl(-(y * y) + m), 'jerfpos' + tag)
    return _v(sympy.Piecewise(
        (pos, y > 0.0), (2.0 * hdl.expl(m) - pos, True)), 'jerfc' + tag)


def leakage(v, par, tag=''):
    """The reverse-leakage terms: Shockley-Read-Hall recombination,
    trap-assisted tunnelling, band-to-band tunnelling, and the avalanche
    multiplication that multiplies all of them
    (`JUNCAP200_macrodefs.include:269-322`).

    These shape a current of order 1e-15 A, eleven orders below the
    drain current, so nothing here moves an operating point.  What they
    do is make the junction's REVERSE characteristic right, which the
    ideal diode alone gets wrong by four orders.
    """
    v = _clamped(v, 'L' + tag)
    phitd, phitr = par['jphit'], par['jphitr']
    vmax = par['jvmax']
    inv_phit = _v(hdl.safe_div(1.0, phitd, eps=1e-30), 'jLip' + tag)

    ## `zinv = exp(V/2phit)`, with PSP's linearisation above `VMAX`.
    ## Each arm's input is clamped to its own domain: the vendor's high
    ## arm takes `sqrt(1 + (V - VMAX)/phit)`, and that radicand is -1e18
    ## at zero bias.
    zlo = _v(hdl.expl(0.5 * sympy.Min(v, vmax) * inv_phit), 'jzlo' + tag)
    zhi = _v(sympy.sqrt(sympy.Max(
        (1.0 + sympy.Max(v - vmax, 0.0) * inv_phit)
        * hdl.expl(vmax * inv_phit), 1e-300)), 'jzhi' + tag)
    zinv = _v(sympy.Piecewise((zlo, v < vmax), (zhi, True)), 'jzinv' + tag)

    ## `two_psistar`.  BOTH of the vendor's arms are kept, and the
    ## reason is worth recording because dropping one looked safe.
    ##
    ## The two are algebraically identical wherever
    ## `zinv = exp(V/2phit)` -- so it is tempting to keep only the
    ## `V <= 0` form, which needs no `z = 1/zinv` and so cannot overflow
    ## when `zinv` underflows.  But `zinv` STOPS being `exp(V/2phit)`
    ## above `VMAX`, where it is the linearisation instead, and there
    ## the identity fails: at a kilovolt forward the `V <= 0` form gives
    ## -998 where the correct arm gives +0.068, and everything
    ## downstream divides by a quantity built from it.
    ##
    ## So both arms stay, each clamped to its own domain: `z` is capped
    ## at 1, which is exact wherever the `V > 0` arm is selected (there
    ## `zinv >= 1`) and merely bounded where it is not.
    z = _v(sympy.Min(hdl.safe_div(1.0, zinv, eps=1e-30), 1.0),
           'jz' + tag)
    tps_pos = _v(2.0 * (phitd * hdl.safe_ln(
        2.0 + z + sympy.sqrt((z + 1.0) * (z + 3.0)))), 'jtpsp' + tag)
    tps_neg = _v(-v + 2.0 * (phitd * hdl.safe_ln(
        2.0 * zinv + 1.0
        + sympy.sqrt((1.0 + zinv) * (1.0 + 3.0 * zinv)))), 'jtpsn' + tag)
    tps = _v(sympy.Piecewise((tps_pos, v > 0.0), (tps_neg, True)),
             'jtps' + tag)

    vjsrh = _hyp2(v, par['jvbimin'] - tps, phitd, 'jvjsrh' + tag)
    vbbt = _hyp2(v, par['jvbbtlim'], phitr, 'jvbbt' + tag)
    vav = _hyp2(v, 0.0, 1.0e-6, 'jvav' + tag)

    total = 0
    for c in COMPONENTS:
        g = lambda k: par['j%s_%s' % (c, k)]
        ct = c + tag

        ## ---- Shockley-Read-Hall ----------------------------------
        dv = _v(sympy.Max(g('vbi') - vjsrh, 1e-30), 'jdv' + ct)
        ## Bounded into (0, 1) by the smoothing above; clamped anyway
        ## because `ln(w)` and `1/(1-w)` sit at the two ends.
        w0 = _v(sympy.Min(sympy.Max(
            1.0 - sympy.sqrt(sympy.Max(1.0 - tps * hdl.safe_div(
                1.0, dv, eps=1e-30), 0.0)), 1e-12), 1.0 - 1e-12),
            'jw0' + ct)
        dw = _v((w0 * w0 * hdl.safe_ln(w0)
                 * hdl.safe_div(1.0, 1.0 - w0, eps=1e-30) + w0)
                * (1.0 - 2.0 * g('p')), 'jdw' + ct)
        wsrh = _v(w0 + dw, 'jwsrh' + ct)
        wdep = _v(g('wdepnulr') * (dv * g('vbirinv')) ** g('p'),
                  'jwdep' + ct)
        ## `ftd` is `n_i(T)/n_i(Tref)`.  It is folded into `idsat` as a
        ## square, and SRH needs it on its own -- and `idsat` is ZERO on
        ## the gate component here, so it cannot be recovered from
        ## there.  Carried explicitly.
        asrh = _v(g('ftd') * ((zinv - 1.0) * wdep), 'jasrh' + ct)
        isrh = _v(g('csrh') * (asrh * wsrh), 'jisrh' + ct)

        ## ---- trap-assisted tunnelling -----------------------------
        btat = _v(g('btatpart') * (wdep * g('omp')
                                   * hdl.safe_div(1.0, dv, eps=1e-30)),
                  'jbtat' + ct)
        tab = _v(JUNCAP_TWOTHIRDS * g('atat')
                 * hdl.safe_div(1.0, btat, eps=1e-30), 'jtab' + ct)
        u0 = _v(tab * tab, 'ju0' + ct)
        ## A smooth `min(u0, 1)`.
        umax = _v(sympy.sqrt(u0 * u0 * hdl.safe_div(1.0, u0 * u0 + 1.0,
                                                    eps=1e-30)),
                  'jumax' + ct)
        su = _v(sympy.sqrt(sympy.Max(umax, 1e-300)), 'jsu' + ct)
        u15 = _v(umax * su, 'ju15' + ct)
        wgam = _v((1.0 + btat * u15) ** (-g('p') * g('oomp')),
                  'jwgam' + ct)
        wtat = _v(wsrh * wgam * hdl.safe_div(1.0, wsrh + wgam,
                                             eps=1e-30), 'jwtat' + ct)
        ktat = _v(sympy.sqrt(sympy.Max(0.375 * (btat * hdl.safe_div(
            1.0, su, eps=1e-30)), 1e-300)), 'jktat' + ct)
        ltat = _v(2.0 * (tab * su) - umax, 'jltat' + ct)
        mtat = _v(g('atat') * tab * su - g('atat') * umax
                  + 0.5 * (btat * u15), 'jmtat' + ct)
        xerfc = _v((ltat - 1.0) * ktat, 'jxerfc' + ct)
        gmax = _v(JUNCAP_SQRTPI * 0.5
                  * (g('atat') * _erfc_exp(xerfc, mtat, par, ct)
                     * hdl.safe_div(1.0, ktat, eps=1e-30)), 'jgmax' + ct)
        itat = _v(g('ctat') * (asrh * gmax * wtat), 'jitat' + ct)

        ## ---- band-to-band tunnelling ------------------------------
        vb = _v(sympy.Max((g('vbir') - vbbt) * g('vbirinv'), 1e-30),
                'jvb' + ct)
        fmax = _v(g('oomp') * ((g('vbir') - vbbt) * g('wdepnulrinv')
                               * hdl.safe_div(1.0, vb ** g('p'),
                                              eps=1e-30)), 'jfmax' + ct)
        ibbt = _v(g('cbbt') * (v * (fmax * fmax)
                               * hdl.expl(-g('fbbt') * hdl.safe_div(
                                   1.0, fmax, eps=1e-30))), 'jibbt' + ct)

        ## ---- avalanche multiplication -----------------------------
        ## The vendor's low arm divides by `1 - (vav/VBR)^PBR`, which is
        ## a POLE at `vav = -VBR` -- and `-VBR` is only ten volts away.
        ## The branch exists to stop short of it, but the arm is still
        ## evaluated there, so the denominator is floored at the value
        ## it has at the crossover.  That is exact inside the branch and
        ## clamps at `fstop` outside it.
        t4 = _v(sympy.Abs(vav * g('vbrinv')) ** g('pbr'), 'jt4' + ct)
        fb_lo = _v(hdl.safe_div(
            1.0, sympy.Max(1.0 - t4, 1.0 / g('fstop')), eps=1e-30),
            'jfblo' + ct)
        fb_hi = _v(g('fstop') + (vav + par['jalphaav'] * g('vbr'))
                   * g('slope'), 'jfbhi' + ct)
        fbrk = _v(sympy.Piecewise(
            (fb_lo, vav > -par['jalphaav'] * g('vbr')), (fb_hi, True)),
            'jfbrk' + ct)

        total = total + g('area') * ((isrh + itat + ibbt) * fbrk)
    return _v(total, 'jleak' + tag)
