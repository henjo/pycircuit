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
