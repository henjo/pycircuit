"""PSP's surface-potential kernel, in the Behavioural DSL.

The heart of every surface-potential compact model is an *explicit*
approximation to the root of the surface-potential equation (SPE).  PSP
does not iterate for it: it computes a good starting point and applies
two closed-form corrections, so the whole thing is a straight-line
expression the simulator can differentiate.

In normalised variables the SPE is

.. math::

    F(x) = (x_g - x)^2 - G_f^2 \\left[ e^{-x} + x - 1
           + \\delta \\left( e^{x} - x - 1
           - \\frac{x^2}{2 + x^2} \\right) \\right] = 0

with :math:`x` the normalised surface potential, :math:`x_g` the
normalised gate drive, :math:`G_f` the normalised body factor, and
:math:`\\delta = e^{-x_n}` where :math:`x_n` is the normalised
quasi-Fermi splitting.  That form is not quoted from a paper -- it is
read straight out of the shipped code, where `sp_s`'s final correction
computes exactly this residual as its ``qC`` term
(`psp103/PSP103_macrodefs.include:163`).

The form is checked for internal consistency rather than taken on trust:
the macro's ``xi1`` is exactly :math:`d\\xi_0/dx = 4x/(2+x^2)^2`, which
is what it must be for ``pC`` to be :math:`-F'(x)` and the final step to
be the Halley correction it is written as.

Translated from ``psp103/PSP103_macrodefs.include``.  **Note the PDK ships
two versions of this macro**: the PSP-based `mosvar` varactor carries a
simplified `sp_s` without the :math:`\\xi_0/\\xi_1/\\xi_2` terms, so it
solves a slightly different equation.  This module implements PSP103's,
since PSP103 is the target.

Everything here builds sympy expressions.  Inside a model being compiled
they register themselves as let-chain intermediates through `hdl.var`, so
the reuse does not explode into a tree; outside one they return plain
expressions, so they can be lambdified and tested directly.
"""

import sympy

from pycircuit.circuit import hdl
from pycircuit.circuit.psp_scaling import PSP_QELE


#: PSP's own constants, by their macro names.
X1 = 1.25
ONE_THIRD = 3.3333333333333333e-01
ONE_SIXTH = 0.1666666666666667
INV_SQRT2 = 7.0710678118654746e-01
SE05 = 2.3025850929940458e+02
KE05 = 1.0e-100

#: The magic constant in `sp_s`'s inversion-side initial guess.
_XG1_COEF = 7.324648775608221e-001


def _v(expr, name):
    """Name an intermediate -- always, never as a bare expression.

    Every intermediate here MUST become a let-chain entry.  There is no
    sensible fallback: the reuse that makes `var()` necessary inside a
    model makes the flattened expression equally unbuildable outside one.
    Measured -- building `sp_s` as one expression did not finish in ten
    minutes, where the chain takes under a second.

    So standalone use goes through `hdl.compile_chain`, which supplies
    the surrounding chain context; calling the kernel with no context at
    all raises, from `var()`, with a message saying so.
    """
    return hdl.var(expr, name)


def maxa(x, y, a):
    """PSP's ``MAXA`` -- a smooth maximum.

    ``0.5*(x + y + sqrt((x-y)^2 + a))``, the companion to `mina`
    (`PSP103_macrodefs.include:38`).  Branch-free, and the radicand is
    bounded below by ``a > 0``.
    """
    x, y = sympy.sympify(x), sympy.sympify(y)
    return 0.5 * (x + y + sympy.sqrt((x - y) ** 2 + a))


def mina(x, y, a):
    """PSP's ``MINA`` -- a smooth minimum.

    ``0.5*(x + y - sqrt((x-y)^2 + a))``.  Branch-free, and the radicand
    is bounded below by ``a > 0``, so nothing needs guarding.
    """
    x, y = sympy.sympify(x), sympy.sympify(y)
    return 0.5 * (x + y - sympy.sqrt((x - y) ** 2 + a))


def _sigma_body(a, ab, c, tau, eta, tag):
    """Shared body of ``sigma`` and ``sigma2``.

    They differ only in whether the ``a`` in two places is ``a`` or
    ``a*b``; PSP ships them as separate macros, and writing the shared
    algebra once keeps them provably the same shape.

    EVERY INTERMEDIATE IS HELD, and that is load-bearing rather than
    tidy.  Left as bare sub-expressions these get substituted into one
    another before printing, so `denom` ends up multiplied out inside
    the numerator and the compiled form evaluates products the written
    form never forms -- `nu**4` against `a*nu*tau` among them.  Written
    as a chain the same algebra is finite at arguments where the
    expanded form overflows to `inf` and then to `NaN`: measured, the
    n-channel long device with card parameters lost its Jacobian at
    ``Vd = 1e26`` and its current at ``1e34``, purely from evaluation
    ORDER.  The `var()` let-chain exists for exactly this.
    """
    nu = _v(a + c, 'sg_nu' + tag)
    mutau = _v(nu * nu + (0.5 * c * c - ab) * tau, 'sg_mutau' + tag)
    ## `mutau` appears in a denominator twice.  It is not guaranteed
    ## nonzero by anything in the algebra, and a zero here would put a
    ## NaN into the surface potential -- and so into every current the
    ## model computes.  Regularised rather than guarded, because a guard
    ## would still evaluate the division in its untaken arm.
    inv_mutau = _v(hdl.safe_div(1.0, mutau, eps=1e-30), 'sg_imu' + tag)
    denom = _v(mutau + (nu * tau * tau * inv_mutau)
               * c * (c * c * ONE_THIRD - ab), 'sg_den' + tag)
    return eta + _v(a * nu * tau, 'sg_num' + tag) \
        * hdl.safe_div(1.0, denom, eps=1e-30)


def sigma(a, c, tau, eta, tag=''):
    """PSP's ``sigma`` -- the three-argument correction."""
    a, c, tau, eta = map(sympy.sympify, (a, c, tau, eta))
    return _sigma_body(a, a, c, tau, eta, tag)


def sigma2(a, b, c, tau, eta, tag=''):
    """PSP's ``sigma2`` -- the four-argument variant.

    The vendor guards ``|tau| < 1e-120`` and returns ``eta``, because
    "sometimes tau is extremely small".  Kept, as a `Piecewise`: with
    `safe_div` underneath, the other arm is finite there anyway, but the
    guard is what the reference model computes and the point is to match
    it.
    """
    a, b, c, tau, eta = map(sympy.sympify, (a, b, c, tau, eta))
    return sympy.Piecewise((eta, sympy.Abs(tau) < 1e-120),
                           (_sigma_body(a, a * b, c, tau, eta, tag),
                            True))


def spe_residual(x, xg, xn, Gf):
    """The surface-potential equation, as a residual.

    Zero at the true surface potential.  This is what `sp_s` approximates
    the root of, and it is the honest way to check that it does: no
    vendor binary, no model card, just the equation the model's own
    correction term encodes.
    """
    x, xg, xn, Gf = map(sympy.sympify, (x, xg, xn, Gf))
    delta = sympy.exp(-xn)
    xi0 = x * x / (2.0 + x * x)
    return ((xg - x) ** 2
            - Gf ** 2 * (sympy.exp(-x) + x - 1
                         + delta * (sympy.exp(x) - x - 1 - xi0)))


def sp_s(xg, xn, delta, Gf, xi, margin=1e-5):
    """PSP's ``sp_s`` -- the surface potential, explicitly.

    Returns the normalised surface potential as a closed-form expression
    in ``xg`` (gate drive), ``xn`` (quasi-Fermi splitting), ``delta``
    (``exp(-xn)``), ``Gf`` (body factor) and ``xi``.

    Three regimes, selected by ``xg``:

    * ``|xg| <= margin`` -- a series expansion about the flat-band point,
      because the general formulae divide by quantities that vanish there;
    * ``xg < -margin`` -- accumulation, solved in ``y = -x``;
    * otherwise -- depletion and inversion.

    Those conditions are **bias-dependent**, so all three arms are
    evaluated at every bias and each must be finite everywhere.  The
    vendor's arms are not: the accumulation branch takes ``ln`` of a
    quantity that is negative in the inversion regime, the inversion
    branch takes ``sqrt`` of one that is negative in accumulation, and
    its ``exp(-eta)`` overflows there.  All are made safe here
    (hdl.md sec. 3.2c) -- in the regime where an arm is actually
    selected its arguments are well inside range, so the values are
    unchanged.

    Note this is stricter than "the selected arm must be finite".  The
    chain rule multiplies a `Piecewise`'s partial derivative -- which is
    ZERO for an unselected arm -- by that arm's own derivative, and
    ``0 * nan`` is ``nan``.  So an unselected arm that merely produces
    garbage still destroys the Jacobian, even though the value is fine.
    That is why the one-sided `expl_low`/`expl_high` are not used here
    even where the vendor uses a bare ``exp``.
    """
    xg, xn, delta, Gf, xi = map(sympy.sympify, (xg, xn, delta, Gf, xi))

    ## EACH ARM IS EVALUATED AT AN ARGUMENT INSIDE ITS OWN DOMAIN.
    ##
    ## Making every individual operation safe is not enough here, and
    ## chasing them one at a time does not converge.  Measured: with all
    ## intermediates finite, the inversion arm at xg = -77 still produced
    ## `ei ~ 1e114` -- finite, but large enough that the products in
    ## `sigma2`'s derivative overflow to inf, and inf - inf is NaN.  The
    ## Jacobian was destroyed by an arm that was never selected.
    ##
    ## Clamping the arm's INPUT fixes the whole class at once: inside the
    ## region where an arm is selected the clamp is the identity, so the
    ## values are exact; outside it the arm computes a bounded number
    ## nobody looks at.
    xg_flat = _v(sympy.Max(sympy.Min(xg, margin), -margin), 'xg_flat')
    xg_acc = _v(sympy.Min(xg, -margin), 'xg_acc')
    xg_inv = _v(sympy.Max(xg, margin), 'xg_inv')

    Gf2 = _v(Gf * Gf, 'Gf2')
    inv_Gf2 = _v(hdl.safe_div(1.0, Gf2, eps=1e-30), 'inv_Gf2')
    inv_xi = _v(hdl.safe_div(1.0, xi, eps=1e-30), 'inv_xi')

    ## ---- flat-band expansion -------------------------------------
    t1 = inv_xi * inv_xi * ONE_SIXTH * INV_SQRT2
    sp_flat = _v(xg_flat * inv_xi
                 * (1.0 + xg_flat * (1.0 - delta) * Gf * t1), 'sp_flat')

    ## ---- accumulation, in y = -x ---------------------------------
    yg = _v(-xg_acc, 'yg')
    ysub = _v(1.25 * yg * inv_xi, 'ysub')
    eta_a = _v(0.5 * (ysub + 10
                      - sympy.sqrt((ysub - 6.0) ** 2 + 64.0)), 'eta_a')
    ta = _v(yg - eta_a, 'ta')
    a_a = _v(ta * ta + Gf2 * (eta_a + 1.0), 'a_a')
    c_a = _v(2.0 * ta - Gf2, 'c_a')
    ## `a_a` is a sum of squares plus Gf2*(eta+1) and is positive in the
    ## accumulation regime, but `eta_a` goes below -1 elsewhere and the
    ## log would take a non-positive argument.
    tau_a = _v(-eta_a + hdl.safe_ln(a_a * inv_Gf2), 'tau_a')
    y0 = _v(sigma(a_a, c_a, tau_a, eta_a, '_a'), 'y0')
    d0_a = _v(hdl.expl(y0), 'd0_a')
    d1_a = _v(hdl.safe_div(1.0, d0_a, eps=1e-30), 'd1_a')
    ## xi0 = y^2/(2+y^2) and its first two derivatives, which is what
    ## makes pC and the curvature term the true -F' and F''.
    t_a = _v(1.0 / (2.0 + y0 * y0), 't_a')
    xi0_a = _v(y0 * y0 * t_a, 'xi0_a')
    xi1_a = _v(4.0 * (y0 * t_a * t_a), 'xi1_a')
    xi2_a = _v((8.0 * t_a - 12.0 * xi0_a) * t_a * t_a, 'xi2_a')
    tb = _v(yg - y0, 'tb')
    dt_a = _v(delta * d1_a, 'dt_a')
    pC_a = _v(2.0 * tb + Gf2 * (d0_a - 1.0 - dt_a
                                + delta * (1.0 - xi1_a)), 'pC_a')
    qC_a = _v(tb * tb - Gf2 * (d0_a - y0 - 1.0 + dt_a
                               + delta * (y0 - 1.0 - xi0_a)), 'qC_a')
    r_a = _v(2.0 - Gf2 * (d0_a + dt_a - delta * xi2_a), 'r_a')
    disc_a = _v(pC_a * pC_a - 2.0 * qC_a * r_a, 'disc_a')
    sp_acc = _v(-y0 - 2.0 * qC_a
                * hdl.safe_div(1.0, pC_a + hdl.safe_sqrt(disc_a),
                               eps=1e-30), 'sp_acc')

    ## ---- depletion and inversion ---------------------------------
    xg1 = _v(1.0 / (X1 + Gf * _XG1_COEF), 'xg1')
    A_fac = _v((xi * X1 * xg1 - 1.0) * xg1, 'A_fac')
    xbar = _v(xg_inv * inv_xi * (1.0 + A_fac * xg_inv), 'xbar')
    w = _v(1.0 - hdl.expl(-xbar), 'w')
    ## `xg + Gf2/4 - w` is negative in accumulation, where this arm is
    ## discarded; the sqrt has to survive being evaluated there anyway.
    x1 = _v(xg_inv + Gf2 * 0.5
            - Gf * hdl.safe_sqrt(xg_inv + Gf2 * 0.25 - w), 'x1')
    bx = _v(xn + 3.0, 'bx')
    eta_i = _v(mina(x1, bx, 5.0)
               - 0.5 * (bx - sympy.sqrt(bx * bx + 5.0)), 'eta_i')
    ti = _v(xg_inv - eta_i, 'ti')
    ei = _v(hdl.expl(-eta_i), 'ei')
    t_e = _v(1.0 / (2.0 + eta_i * eta_i), 't_e')
    xi0_e = _v(eta_i * eta_i * t_e, 'xi0_e')
    xi1_e = _v(4.0 * (eta_i * t_e * t_e), 'xi1_e')
    xi2_e = _v((8.0 * t_e - 12.0 * xi0_e) * t_e * t_e, 'xi2_e')
    a_i = _v(sympy.Max(1.0e-40,
                       ti * ti - Gf2 * (ei + eta_i - 1.0
                                        - delta * (eta_i + 1.0
                                                   + xi0_e))), 'a_i')
    b_i = _v(1.0 - 0.5 * (Gf2 * (ei - delta * xi2_e)), 'b_i')
    c_i = _v(2.0 * ti + Gf2 * (1.0 - ei
                               - delta * (1.0 + xi1_e)), 'c_i')
    tau_i = _v(xn - eta_i + hdl.safe_ln(a_i * inv_Gf2), 'tau_i')
    x0 = _v(sigma2(a_i, b_i, c_i, tau_i, eta_i, '_i'), 'x0')

    ## The vendor splits `delta0`/`delta1` three ways to keep both
    ## exp(x0) and exp(x0 - xn) inside range; `expl` does the same job in
    ## one expression and is C-3 rather than piecewise-continuous.
    ex0 = _v(hdl.expl(x0), 'ex0')
    d0_i = _v(delta * ex0, 'd0_i')
    ## A plain floored reciprocal rather than `safe_div`.  `ex0` is an
    ## exponential and therefore strictly positive, so `safe_div`'s
    ## regularisation -- which is `b/(b^2 + eps^2)` -- buys nothing here
    ## and costs half the exponent range: it SQUARES `ex0`, and `expl`
    ## is allowed to return 1e100 by construction, so the square
    ## overflows on ordinary use of the model at a bias a diverging
    ## Newton step reaches.  The floor is there only so the division
    ## cannot meet a value that underflowed to zero.
    d1_i = _v(1.0 / sympy.Max(ex0, 1e-300), 'd1_i')

    t_x = _v(1.0 / (2.0 + x0 * x0), 't_x')
    xi0_x = _v(x0 * x0 * t_x, 'xi0_x')
    xi1_x = _v(4.0 * (x0 * t_x * t_x), 'xi1_x')
    xi2_x = _v((8.0 * t_x - 12.0 * xi0_x) * t_x * t_x, 'xi2_x')
    tc = _v(xg_inv - x0, 'tc')
    pC_i = _v(2.0 * tc + Gf2 * (1.0 - d1_i + d0_i
                                - delta * (1.0 + xi1_x)), 'pC_i')
    qC_i = _v(tc * tc - Gf2 * (d1_i + x0 - 1.0 + d0_i
                               - delta * (x0 + 1.0 + xi0_x)), 'qC_i')
    r_i = _v(2.0 - Gf2 * (d1_i + d0_i - delta * xi2_x), 'r_i')
    disc_i = _v(pC_i * pC_i - 2.0 * qC_i * r_i, 'disc_i')
    sp_inv = _v(x0 + 2.0 * qC_i
                * hdl.safe_div(1.0, pC_i + hdl.safe_sqrt(disc_i),
                               eps=1e-30), 'sp_inv')

    return sympy.Piecewise((sp_flat, sympy.Abs(xg) <= margin),
                           (sp_acc, xg < -margin),
                           (sp_inv, True))


def surface_state(x, xn, delta, Gf, inv_Gf2):
    """The charge quantities at one surface potential.

    Given the normalised surface potential ``x`` at one end of the
    channel, returns ``(E, P, sq, D)``:

    * ``E = exp(-x)``, the depletion term;
    * ``P = x - 1 + E``, whose square root ``sq`` gives the bulk charge;
    * ``D = delta*(exp(x) - x - 1 - xi0)``, the inversion term -- which
      is exactly the ``delta`` group of `spe_residual`, so the current
      layer and the surface-potential solver are demonstrably built on
      the same equation.

    PSP splits each of these three ways to keep the exponentials in
    range.  Here the middle form is used with `expl`, which is C-3 and
    covers the whole line, and the SMALL-``x`` series is kept because it
    is not a range guard -- it is there because ``x - 1 + exp(-x)``
    cancels to nothing near zero, losing every significant digit.
    """
    x, xn, delta, Gf = map(sympy.sympify, (x, xn, delta, Gf))

    ## Each arm clamped to its own side: the series is valid near zero,
    ## the general form away from it (hdl.md sec. 3.2c).
    small = 1.0e-5
    xs = _v(sympy.Min(sympy.Abs(x), small), 'ss_xs')
    xb = _v(sympy.Max(x, small), 'ss_xb')

    E = _v(hdl.expl(-x), 'ss_E')
    xi0 = _v(x * x / (2.0 + x * x), 'ss_xi0')

    ## Near zero: P = x^2/2 * (1 - x/3*(1 - x/4)) to leading orders, and
    ## sq follows in closed form rather than through the cancellation.
    t_s = _v(sympy.sqrt(1.0 - ONE_THIRD * (xs * (1.0 - 0.25 * xs))), 'ss_ts')
    P_small = _v(0.5 * (xs * xs
                        * (1.0 - ONE_THIRD * (xs * (1.0 - 0.25 * xs)))),
                 'ss_Ps')
    sq_small = _v(INV_SQRT2 * (x * t_s), 'ss_sqs')
    D_small = _v(ONE_SIXTH * delta * x * x * x * (1.0 + 1.75 * x), 'ss_Ds')

    P_big = _v(xb - 1.0 + hdl.expl(-xb), 'ss_Pb')
    sq_big = _v(hdl.safe_sqrt(P_big), 'ss_sqb')
    D_big = _v(delta * (hdl.expl(xb) - xb - 1.0
                        - xb * xb / (2.0 + xb * xb)), 'ss_Db')

    near = sympy.Abs(x) < small
    P = _v(sympy.Piecewise((P_small, near), (P_big, True)), 'ss_P')
    sq = _v(sympy.Piecewise((sq_small, near), (sq_big, True)), 'ss_sq')
    D = _v(sympy.Piecewise((D_small, near), (D_big, True)), 'ss_D')
    return E, P, sq, D


#: Floor on the normalised quasi-Fermi splitting.
#:
#: `xn < 0` means a junction forward-biased past its built-in potential
#: -- outside the domain of any compact model, and numerically fatal
#: here: `delta = exp(-xn)` reaches 1e152 by `xn = -352`, and the
#: products it enters overflow to NaN.  A real device runs at xn ~ 35, so
#: this never binds in operation; it exists because Newton visits places
#: the device cannot go, and a bounded wrong answer there is worth far
#: more than a NaN.  PSP clamps its junction voltages for the same
#: reason.
XN_FLOOR = 0.0


def _bias_mod(coeff, v, name):
    """PSP's ``1 + t*v`` / ``1/(1 - t*v)`` pair.

    Used wherever PSP modulates a parameter by a bias
    (`PSP103_macrodefs.include:596-607`).  The second form is not a
    different model, it is the same first-order modulation written so
    that the factor stays POSITIVE when the coefficient is negative --
    `1 + t*v` would cross zero and change the sign of whatever it
    multiplies.  The two agree to first order in `t*v`.

    The branch is on the COEFFICIENT, which is a parameter, so the
    condition carries no bias dependence and the derivative is a plain
    `where` on a constant -- none of the `Max`-on-an-expression trouble
    documented in `intrinsic`.
    """
    coeff, v = sympy.sympify(coeff), sympy.sympify(v)
    return _v(sympy.Piecewise(
        (hdl.safe_div(1.0, 1.0 - coeff * v, eps=1e-30), coeff < 0),
        (1.0 + coeff * v, True)), name)


def _wsat(q, xitsb, thesatg, tag):
    """PSP's gate-bias modulation of the velocity-saturation parameter.

    ``wsat`` is a soft-limited inversion charge: ``100*q'/(100 + q')``
    with ``q' = q*xitsb``, so for the ordinary sub-volt charges it IS
    ``q'`` and it saturates rather than running away.  The returned
    factor multiplies ``THESAT``.
    """
    q2 = _v(q * xitsb, 'ids_qx' + tag)
    w = _v(100.0 * (q2 * hdl.safe_div(1.0, 100.0 + q2, eps=1e-30)),
           'ids_wsat' + tag)
    return _bias_mod(thesatg, w, 'ids_xitsg' + tag)


def intrinsic(xg, xn_s, xn_d, Gf, xi, phit, beta, margin=1e-5,
              mob=None):
    """The intrinsic long-channel core: current, and what charge needs.

    PSP assembles the current as
    ``Ids = BET*(FdL * qim1 * dps * Gvsatinv)``
    (`PSP103_module.include:1178`).  ``FdL`` is channel-length
    modulation and ``Gvsatinv`` velocity saturation; with both at unity
    -- the long-channel intrinsic device -- what remains is

    .. math::  I_{ds} = \\beta \\, q_{im1} \\, \\Delta\\psi

    the symmetric-linearisation charge-sheet current.  ``dps`` is the
    surface-potential difference between the channel ends and ``qim1``
    the inversion charge evaluated at the MIDPOINT potential, which is
    what makes the result exactly antisymmetric under swapping source
    and drain -- the property threshold-voltage models are famous for
    getting wrong, and the reason surface-potential models exist.

    `beta` is ``mu * Cox * W / L``; `phit` the thermal voltage; `xn_s`
    and `xn_d` the normalised quasi-Fermi levels at the two ends.

    `mob`, if given, switches on mobility reduction and velocity
    saturation: a dict of ``mue``, ``themu``, ``cs``, ``thecs``,
    ``feta``, ``thesat``, ``cox_area`` and ``eps_si``.  Left out, the
    device is the ideal long-channel one with ``Gmob = Gvsat = 1``,
    which is what the construction tests use.

    Returns a dict: ``ids`` plus the midpoint quantities
    (``qim``, ``qim1``, ``alpha``, ``dps``, ``xgm``, ``Gmob``,
    ``Gvsat``, ``zsat``) that `charges_long_channel` needs, so that
    current and charge are built from ONE evaluation of the surface
    potentials rather than two.
    """
    xg, xn_s, xn_d = map(sympy.sympify, (xg, xn_s, xn_d))
    Gf, xi, phit, beta = map(sympy.sympify, (Gf, xi, phit, beta))
    xn_s = _v(sympy.Max(xn_s, XN_FLOOR), 'ids_xns')
    xn_d = _v(sympy.Max(xn_d, XN_FLOOR), 'ids_xnd')

    Gf2 = _v(Gf * Gf, 'ids_Gf2')
    inv_Gf2 = _v(hdl.safe_div(1.0, Gf2, eps=1e-30), 'ids_invGf2')
    d_s = _v(hdl.expl(-xn_s), 'ids_ds')
    x_s = _v(sp_s(xg, xn_s, d_s, Gf, xi, margin), 'ids_xs')
    Es, Ps, sqs, Ds = surface_state(x_s, xn_s, d_s, Gf, inv_Gf2)

    ## The DRAIN end is deliberately not solved yet.  Where the
    ## saturation voltage below is active it moves `xn_d`, and solving
    ## the surface potential at the applied bias first would build a
    ## whole `sp_s` -- the most expensive thing in this file -- only for
    ## reachability pruning to throw it away again.

    ## ---- the saturation-limited drain voltage ------------------------
    ## PSP does not evaluate the drain surface potential at the applied
    ## drain bias.  It computes a saturation voltage `Vdsat` from
    ## SOURCE-SIDE quantities, smoothly limits `Vds` to it, and uses the
    ## limited `Vdse` for the drain-side quasi-Fermi level
    ## (`PSP103_macrodefs.include:596-632`).
    ##
    ## Leaving this out is what made the channel-length modulation an
    ## approximation: `dL` compares `Vds` against `Vdse`, and with no
    ## `Vdse` the second logarithm had to be dropped.  `FdL` then
    ## amplifies whatever error `dL` carries, which is why adding `FdL`
    ## cost the long device accuracy even while more than halving the
    ## total.  This closes that.
    ## Channel type.  A PLAIN PYTHON BOOL, never a symbol: PSP decides
    ## this once in `INITIAL_MODEL` (`PSP103_module.include:557-561`),
    ## it cannot vary with bias, and branching on it here keeps the
    ## compiled expression exactly the size it is for one channel type.
    pmos = bool(mob.get('pmos', False)) if mob else False
    eta_mu = ONE_THIRD if pmos else 0.5

    ## `Rxcor` (`PSP103_macrodefs.include:576`), which multiplies `Gmob`
    ## at BOTH the source end and the midpoint (`:595`, `:750`).  It is
    ## identically 1 at zero body bias, so it costs nothing where it does
    ## not act; computed once here because PSP computes it once.
    if mob is None:
        rxcor = sympy.Integer(1)
    else:
        _xc = mob.get('xcor', 0.0)
        _vsb = mob.get('vsb', 0.0)
        rxcor = _v((1.0 + 0.2 * _xc * _vsb)
                   * hdl.safe_div(1.0, 1.0 + _xc * _vsb, eps=1e-30),
                   'ids_rxcor')

    ax = 0.0 if mob is None else mob.get('ax', 0.0)
    if not isinstance(ax, sympy.Expr) and ax == 0.0:
        vdsat = None
        vdse = None
        xitsb = None
        d_d = _v(hdl.expl(-xn_d), 'ids_dd')
        x_d = _v(sp_s(xg, xn_d, d_d, Gf, xi, margin), 'ids_xd')
    else:
        xgs = _v(Gf * hdl.safe_sqrt(Ps + Ds), 'ids_xgs')
        qis = _v(Gf2 * Ds * phit
                 * hdl.safe_div(1.0, xgs + Gf * sqs, eps=1e-30), 'ids_qis')
        qbs = _v(sqs * Gf * phit, 'ids_qbs')
        alpha_s = _v(1.0 + 0.5 * (Gf * (1.0 - Es)
                                  * hdl.safe_div(1.0, sqs, eps=1e-30)),
                     'ids_alphas')
        ## Mobility at the SOURCE end, same form as the midpoint one.
        e0 = _v(1.0e-8 * mob['cox_area'] / mob['eps_si'], 'ids_ee0s')
        ## `eta_mu` is 1/2 for electrons and 1/3 for holes
        ## (`PSP103_module.include:736-743`) -- the inversion charge's
        ## weight in the effective vertical field, which differs because
        ## the two carrier types sit at different average depths below
        ## the interface.
        bs0 = _v(e0 * (qbs + eta_mu * mob['feta'] * qis) * mob['mue'],
                 'ids_mubs0')
        ## Bounded ABOVE as well as below.  `THEMU` is about 2 on this
        ## card, so a base past 1e154 overflows the power -- and a
        ## diverging Newton step, which is bounded in the device's
        ## BRANCH voltages but not in the node voltages themselves,
        ## reaches it.  The upper clamp is far outside any physical
        ## effective field; it exists so the mobility degrades to "very
        ## large" instead of to `inf`, which would go on to make an
        ## `inf - inf` somewhere downstream.
        bs_ = _v(sympy.Min(sympy.Max(bs0, 1e-300), 1e100), 'ids_mubs')
        rs0 = _v(hdl.safe_div(Ps, Ps + Ds + 1.0e-14, eps=1e-30),
                 'ids_ratios0')
        rs_ = _v(sympy.Max(rs0, 1.0e-10), 'ids_ratios')
        gmob_s = _v((1.0 + bs_ ** mob['themu']
                     + mob['cs'] * rs_ ** (0.5 * mob['thecs']))
                    * rxcor, 'ids_gmobs')
        ## The body- and gate-bias modulations of `THESAT`
        ## (`PSP103_macrodefs.include:596-607`).  `xitsb` is computed
        ## HERE, once, and reused unchanged at the midpoint below --
        ## PSP does the same, carrying it as `xitsb_dc`.  PSP guards
        ## this whole block with `Ds > ke05` and leaves `xitsb` at zero
        ## otherwise; we have no guard because we do not need one, `qis`
        ## being proportional to `Ds` and therefore already vanishing
        ## there, which makes `wsat` vanish and the factor collapse to 1.
        xitsb = _bias_mod(mob.get('thesatb', 0.0), mob.get('vsb', 0.0),
                          'ids_xitsb')
        th1 = _v(mob['thesat'] * _wsat(qis, xitsb, mob.get('thesatg', 0.0),
                                       's')
                 * hdl.safe_div(1.0, gmob_s, eps=1e-30), 'ids_th1s')
        phi_inf = _v(qis * hdl.safe_div(1.0, alpha_s, eps=1e-30) + phit,
                     'ids_phiinf')
        ysat = _v(th1 * phi_inf * INV_SQRT2, 'ids_ysat')
        if pmos:
            ## `macrodefs:609-613`.  This and the `zsat` form below are
            ## the same physical statement: electrons follow
            ## `v ~ E/sqrt(1 + (E/Ec)^2)` and holes `v ~ E/(1 + E/Ec)`,
            ## so the hole branch is one power softer in the field.
            ysat = _v(ysat * hdl.safe_div(
                1.0, sympy.sqrt(1.0 + ysat), eps=1e-30), 'ids_ysatp')
        za = _v(2.0 * hdl.safe_div(1.0, 1.0 + sympy.sqrt(1.0 + 4.0 * ysat),
                                   eps=1e-30), 'ids_za')
        t1 = _v(za * ysat, 'ids_t1s')
        Phi_0 = _v(phi_inf * za
                   * (1.0 + 0.86 * (t1 * (1.0 - t1 * za)
                                    * hdl.safe_div(1.0, 1.0 + 4.0
                                                   * (t1 * t1 * za),
                                                   eps=1e-30))), 'ids_Phi0')
        asat = _v(xgs + 0.5 * Gf2, 'ids_asat')
        Phi_2 = _v(0.98 * (Gf2 * Ds * phit
                           * hdl.safe_div(
                               1.0, asat + hdl.safe_sqrt(asat * asat
                                                         - Gf2 * Ds * 0.98),
                               eps=1e-30)), 'ids_Phi2')
        p02 = _v(Phi_0 + Phi_2, 'ids_p02')
        pp2 = _v(2.0 * (Phi_0 * Phi_2), 'ids_pp2')
        Phi_sat = _v(pp2 * hdl.safe_div(
            1.0, p02 + hdl.safe_sqrt(p02 * p02 - 1.98 * pp2), eps=1e-30),
            'ids_Phisat')
        ## NOTE THE SHAPE OF THESE TWO LINES.  `Max` must be applied to
        ## an ATOM -- a `_v` symbol -- never to an expression.  sympy
        ## differentiates `Max(f, c)` by expanding `f` into the branch
        ## conditions, so a large `f` is duplicated into every one of
        ## them; and when `f` contains a `hypsmooth`, sympy rationalises
        ## it there as `2e^2/(sqrt(z^2 + 4e^2) - z)`, whose denominator
        ## cancels to EXACTLY zero for any ordinary `z` because `4e^2`
        ## sits far below `z`'s last bit.  The value was finite and the
        ## Jacobian divided by zero.  Naming the argument first keeps it
        ## opaque and the derivative stays a one-line `where`.
        lnarg = _v(sympy.Max(
            1.0 + Phi_sat * (Phi_sat - 2.0 * asat * phit) * inv_Gf2
            * hdl.safe_div(1.0, phit * phit * Ds, eps=1e-30),
            1.0e-3), 'ids_lnarg')
        vdsat0 = _v(Phi_sat - phit * sympy.log(lnarg), 'ids_vdsat0')
        vdsat = _v(sympy.Max(vdsat0, 1.0e-6), 'ids_vdsat')

        ## `Vds` is recoverable from the arguments, since `xn_d - xn_s`
        ## IS `Vds/phit` by construction.  Limiting it and rebuilding
        ## `xn_d` keeps the caller's interface unchanged.
        vds_in = _v((xn_d - xn_s) * phit, 'ids_vdsin')
        ## Floored strictly ABOVE zero, not at it.  sympy differentiates
        ## `rat**ax` as `ax*rat**ax/rat`, so a base that reaches zero --
        ## which it does at `Vds = 0`, the most ordinary bias there is --
        ## divides by zero in the Jacobian while the value itself is
        ## perfectly finite.  Same trap as the Coulomb scattering term.
        rat0 = _v(vds_in * hdl.safe_div(1.0, vdsat, eps=1e-30), 'ids_vrat0')
        rat = _v(sympy.Max(rat0, 1e-15), 'ids_vrat')
        ## `ax` reaches here as a symbol when the caller is a compiled
        ## element, so a Python-level `ax == 0` test cannot switch this
        ## block off -- it would compile the branch anyway and divide by
        ## zero at run time.  PSP clips AX at 2 (`PSP103_scaling:743`),
        ## so there is no legitimate zero; floor it and be done.
        axc = _v(sympy.Max(ax, 0.5), 'ids_axc')

        ## `Vds*(1 + rat^ax)^(-1/ax)`, written in the two forms that are
        ## each stable on their own side of `rat = 1`.  Directly is fine
        ## below 1 and overflows above it -- `rat^6.5` leaves double
        ## range by `rat = 1e46`, which a diverging Newton step reaches
        ## -- and `huge * (1 + inf)^(-1/ax)` is `huge * 0`, i.e. NaN.
        ## Above 1, factoring `rat` out gives the algebraically identical
        ## `Vdsat*(1 + rat^-ax)^(-1/ax)`, whose base is again below 1.
        ##
        ## Each arm's input is clamped to ITS OWN domain rather than the
        ## arms merely being selected between: a Piecewise's derivative
        ## w.r.t. something used only in the discarded arm is zero, and
        ## zero times NaN is NaN.  Both are continuous at `rat = 1`,
        ## where `Vds` and `Vdsat` coincide by definition.
        rlo = _v(sympy.Max(sympy.Min(rat, 1.0), 1e-15), 'ids_rlo')
        rhi = _v(sympy.Max(hdl.safe_div(1.0, rat, eps=1e-30), 1e-15),
                 'ids_rhi')
        rhic = _v(sympy.Min(rhi, 1.0), 'ids_rhic')
        vdse = _v(sympy.Piecewise(
            (vds_in * (1.0 + rlo ** axc) ** (-1.0 / axc), rat <= 1.0),
            (vdsat * (1.0 + rhic ** axc) ** (-1.0 / axc), True)),
            'ids_vdse')
        xn_d = _v(xn_s + vdse * hdl.safe_div(1.0, phit, eps=1e-30),
                  'ids_xnd2')
        d_d = _v(hdl.expl(-xn_d), 'ids_dd2')
        x_d = _v(sp_s(xg, xn_d, d_d, Gf, xi, margin), 'ids_xd2')

    Ed, Pd, sqd, Dd = surface_state(x_d, xn_d, d_d, Gf, inv_Gf2)

    x_ds = _v(x_d - x_s, 'ids_xds')
    x_m = _v(0.5 * (x_s + x_d), 'ids_xm')

    ## The midpoint depletion term is the GEOMETRIC mean of the two ends,
    ## which is what keeps the construction symmetric.
    Em = _v(hdl.safe_sqrt(Es * Ed), 'ids_Em')
    D_bar = _v(0.5 * (Ds + Dd), 'ids_Dbar')
    Dm = _v(D_bar + 0.125 * (x_ds * x_ds * (Em - 2.0 * inv_Gf2)), 'ids_Dm')

    ## Same small-x care as `surface_state`: `x - 1 + E` cancels at zero.
    small = 1.0e-5
    xms = _v(sympy.Min(sympy.Abs(x_m), small), 'ids_xms')
    xmb = _v(sympy.Max(x_m, small), 'ids_xmb')
    near_m = sympy.Abs(x_m) < small

    t_m = _v(sympy.sqrt(1.0 - ONE_THIRD * (xms * (1.0 - 0.25 * xms))),
             'ids_tm')
    Pm_small = _v(0.5 * (xms * xms
                         * (1.0 - ONE_THIRD * (xms * (1.0 - 0.25 * xms)))),
                  'ids_Pms')
    Pm0 = _v(sympy.Piecewise((Pm_small, near_m), (xmb - 1.0 + Em, True)),
             'ids_Pm0')
    xgm0 = _v(Gf * hdl.safe_sqrt(Dm + Pm0), 'ids_xgm0')

    ## ---- polysilicon depletion --------------------------------------
    ## The gate is not a perfect conductor: a depletion layer forms on
    ## its silicon side, taking a share of the applied voltage.  PSP
    ## folds this in as `eta_p = 1/sqrt(1 + kp*xgm)` on the charge slope
    ## AND as a correction that shifts the midpoint potential itself
    ## (`PSP103_macrodefs.include:697-724`).  Implemented whole rather
    ## than as the `eta_p` factor alone -- half of a coupled correction
    ## is not a smaller version of it.
    ##
    ## `kp = 2*Cox'^2*phit/(q*np*eps_si)`, with `np` floored at both
    ## `8e7/tox^2` and 5e24 as PSP floors it (macrodefs:315-319).  For
    ## the IHP card that is 1.6e-3, worth 1% of the charge slope at low
    ## gate bias and 4.5% at high -- which is why it shows up as a
    ## bias-dependent gain error rather than a constant one.
    kp = 0.0 if mob is None else mob.get('kp', 0.0)
    if kp == 0.0 and not isinstance(kp, sympy.Expr):
        eta_p = _v(sympy.Integer(1), 'ids_etap')
        x_mc, Em_c, Dm_c, x_dsc = x_m, Em, Dm, x_ds
    else:
        eta_p = _v(hdl.safe_div(1.0, sympy.sqrt(1.0 + kp * xgm0),
                                eps=1e-30), 'ids_etap')
        d0 = _v(1.0 - Em + 2.0 * (xgm0 * inv_Gf2), 'ids_d0')
        te = _v(eta_p * hdl.safe_div(1.0, eta_p + 1.0, eps=1e-30),
                'ids_te')
        x_pm = _v(kp * (te * te * Gf2 * Dm), 'ids_xpm')
        p_pd = _v(2.0 * (xgm0 - x_pm) + Gf2 * (1.0 - Em + Dm), 'ids_ppd')
        q_pd = _v(x_pm * (x_pm - 2.0 * xgm0), 'ids_qpd')
        xi_pd = _v(1.0 - 0.5 * (Gf2 * (Em + Dm)), 'ids_xipd')
        u_pd = _v(q_pd * p_pd
                  * hdl.safe_div(1.0, p_pd * p_pd - xi_pd * q_pd,
                                 eps=1e-30), 'ids_upd')
        km = _v(hdl.expl(u_pd), 'ids_km')
        inv_km = _v(hdl.safe_div(1.0, km, eps=1e-30), 'ids_invkm')
        ## The correction applies only where the midpoint is not at flat
        ## band; near zero PSP takes the series branch and leaves the
        ## potential alone, using `eta_p` in `alpha` only.
        x_mc = _v(sympy.Piecewise((x_m, near_m), (x_m + u_pd, True)),
                  'ids_xmc')
        Em_c = _v(sympy.Piecewise((Em, near_m), (Em * inv_km, True)),
                  'ids_Emc')
        Dm_c = _v(sympy.Piecewise((Dm, near_m), (Dm * km, True)),
                  'ids_Dmc')
        Pm_b2 = _v(sympy.Max(x_mc, small) - 1.0 + Em_c, 'ids_Pmb2')
        xgm_c = _v(Gf * hdl.safe_sqrt(Dm_c + Pm_b2), 'ids_xgmc')
        km0 = _v(1.0 - Em_c + 2.0 * (xgm_c * eta_p * inv_Gf2), 'ids_km0')
        x_dsc = _v(sympy.Piecewise(
            (x_ds, near_m),
            (x_ds * km * (d0 + D_bar)
             * hdl.safe_div(1.0, km0 + km * D_bar, eps=1e-30), True)),
            'ids_xdsc')

    ## ---- midpoint charge quantities, after the correction ------------
    xmb2 = _v(sympy.Max(x_mc, small), 'ids_xmb2')
    xms2 = _v(sympy.Min(sympy.Abs(x_mc), small), 'ids_xms2')
    t_m2 = _v(sympy.sqrt(1.0 - ONE_THIRD * (xms2 * (1.0 - 0.25 * xms2))),
              'ids_tm2')
    Pm_small2 = _v(0.5 * (xms2 * xms2
                          * (1.0 - ONE_THIRD
                             * (xms2 * (1.0 - 0.25 * xms2)))), 'ids_Pms2')
    sqm_small = _v(INV_SQRT2 * (x_mc * t_m2), 'ids_sqms')
    alpha_small = _v(eta_p + INV_SQRT2 * (Gf * (1.0 - 0.5 * x_mc
                                                + ONE_SIXTH
                                                * (x_mc * x_mc))
                                          / t_m2), 'ids_as')
    Pm_big = _v(xmb2 - 1.0 + Em_c, 'ids_Pmb')
    sqm_big = _v(hdl.safe_sqrt(Pm_big), 'ids_sqmb')
    alpha_big = _v(eta_p + 0.5 * (Gf * (1.0 - Em_c)
                                  * hdl.safe_div(1.0, sqm_big, eps=1e-30)),
                   'ids_ab')

    Pm = _v(sympy.Piecewise((Pm_small2, near_m), (Pm_big, True)), 'ids_Pm')
    sqm = _v(sympy.Piecewise((sqm_small, near_m), (sqm_big, True)),
             'ids_sqm')
    alpha = _v(sympy.Piecewise((alpha_small, near_m), (alpha_big, True)),
               'ids_alpha')
    Dm = Dm_c
    x_ds = x_dsc

    xgm = _v(Gf * hdl.safe_sqrt(Dm + Pm), 'ids_xgm')
    qim = _v(phit * Gf2 * Dm
             * hdl.safe_div(1.0, xgm + Gf * sqm, eps=1e-30), 'ids_qim')
    qim1 = _v(qim + phit * alpha, 'ids_qim1')
    dps = _v(x_ds * phit, 'ids_dps')
    qbm = _v(sqm * Gf * phit, 'ids_qbm')

    ## ---- mobility reduction and velocity saturation ---------------
    ## Both depend only on MIDPOINT quantities and on `dps` SQUARED, so
    ## neither can break the source/drain antisymmetry: `Gmob` and
    ## `Gvsat` are even under the swap and the lone odd factor stays
    ## `dps`.  That is not a happy accident -- it is why PSP evaluates
    ## them at the midpoint.
    if mob is None:
        Gmob = _v(sympy.Integer(1), 'ids_Gmob')
        Gmob_dL = Gmob
        GdL = _v(sympy.Integer(1), 'ids_GdL')
        dL = _v(sympy.Integer(0), 'ids_dL')
        s1 = _v(sympy.Integer(0), 'ids_s1')
        zsat = _v(sympy.Integer(0), 'ids_zsat')
        Gvsat = _v(sympy.Integer(1), 'ids_Gvsat')
    else:
        eta_mu_m = eta_mu * mob['feta']
        qeff = _v(qbm + eta_mu_m * qim, 'ids_qeff')
        ## PSP's `E_eff0 = 1e-8 * Cox' / eps_si` (module:737).
        e_eff0 = _v(1.0e-8 * mob['cox_area'] / mob['eps_si'], 'ids_ee0')
        base0 = _v(e_eff0 * qeff * mob['mue'], 'ids_mubase0')
        ## See the note on `ids_mubs` above: clamped at both ends.
        base = _v(sympy.Min(sympy.Max(base0, 1e-300), 1e100),
                  'ids_mubase')
        ## Coulomb scattering.  PSP writes this as
        ## `cs*exp(0.5*thecs*ln(Pm/(Pm+Dm)))`, because Verilog-A's `pow`
        ## with a variable exponent is awkward -- but algebraically it is
        ## just `cs * ratio**(0.5*thecs)`, and the power form is the one
        ## to compile.  Composing `expl` with `safe_ln` nests two
        ## `Piecewise`s, and differentiating the nest emits conditions
        ## built from `logical_and.reduce` over scalars that numpy's
        ## `select` rejects outright: the Jacobian did not merely lose
        ## precision, it raised.
        ##
        ## `ratio` is floored because the exponent is below 1 (thecs =
        ## 1.18 gives 0.59), so the derivative `ratio**(0.5*thecs - 1)`
        ## diverges as the ratio goes to zero -- which it does exactly at
        ## flat band, where `Pm` is 0.  Below the floor the term is
        ## constant, which is right: no inversion charge, no scattering
        ## off it.
        ratio = _v(sympy.Max(hdl.safe_div(Pm, Pm + Dm + 1.0e-14,
                                          eps=1e-30), 1.0e-10),
                   'ids_ratio')
        coul = _v(mob['cs'] * ratio ** (0.5 * mob['thecs']), 'ids_coul')
        ## SERIES RESISTANCE.  PSP does not put a resistor in the
        ## network -- it folds the source/drain resistance into the
        ## mobility as `GR = THER*rhob*temp*qim` with
        ## `THER = 2*BET*RS` (`PSP103_macrodefs.include:381, 742`).  That
        ## keeps the device four-terminal and, more to the point, keeps
        ## it SYMMETRIC: `qim` is a midpoint quantity, so `GR` is even
        ## under the source/drain swap like everything else in `Gmob`.
        rs = mob.get('rs', 0.0)
        rsg = mob.get('rsg', 0.0)
        rsb = mob.get('rsb', 0.0)
        vsb = mob.get('vsb', 0.0)
        if rs == 0.0 and not isinstance(rs, sympy.Expr):
            gr = _v(sympy.Integer(0), 'ids_GR')
        else:
            ther = _v(2.0 * beta * rs, 'ids_ther')
            ## Body-bias dependence of the sheet resistance; the two
            ## forms are PSP's, selected by the SIGN of a parameter, so
            ## the condition is parameter-only and the topology is fixed.
            rhob = _v(sympy.Piecewise(
                (hdl.safe_div(1.0, 1.0 - rsb * vsb, eps=1e-30), rsb < 0),
                (1.0 + rsb * vsb, True)), 'ids_rhob')
            rsgf = _v(sympy.Piecewise(
                (1.0 - rsg * qim, rsg < 0),
                (hdl.safe_div(1.0, 1.0 + rsg * qim, eps=1e-30), True)),
                'ids_rsgf')
            gr = _v(ther * (rhob * rsgf * qim), 'ids_GR')

        Gmob = _v((1.0 + base ** mob['themu'] + coul + gr) * rxcor,
                  'ids_Gmob')

        ## CHANNEL-LENGTH MODULATION.  Beyond pinch-off the inversion
        ## layer stops short of the drain, so the effective channel is
        ## shorter and the current keeps rising where an ideal device's
        ## is flat.
        ##
        ## PSP writes `s1` as a RATIO of logarithms involving `Vdse`
        ## (`PSP103_macrodefs.include:753`):
        ##
        ##     s1 = ln((1 + (Vds - dps)/VP) / (1 + (Vdse - dps)/VP))
        ##
        ## The second logarithm used to be dropped, because without a
        ## saturation-limited drain voltage there was nothing to put in
        ## it -- which recovered the classic `ln(1 + (Vds - Vdsat)/VP)`
        ## with `dps` standing in for the saturation voltage, and was an
        ## APPROXIMATION to PSP rather than a translation of it.  `Vdse`
        ## exists now, so this is the real expression.
        ##
        ## The denominator matters more than it looks.  Both logarithms
        ## grow together in the linear region, where `Vdse` tracks `Vds`,
        ## so their ratio stays near 1 and `dL` stays near zero -- which
        ## is the physical statement that there is no pinch-off region to
        ## shorten.  The single-logarithm form had no way to say that and
        ## leaned on `dps` to say it instead.
        ##
        ## `hypsmooth` rather than `max(.,0)`: the excess is negative
        ## through the whole linear region and a hard clamp would put a
        ## kink in the middle of every I-V curve.
        alp = mob.get('alp', 0.0)
        if alp == 0.0 and not isinstance(alp, sympy.Expr):
            GdL = _v(sympy.Integer(1), 'ids_GdL')
            dL = _v(sympy.Integer(0), 'ids_dL')
            s1 = _v(sympy.Integer(0), 'ids_s1')
        else:
            ## `vds` arrives already symmetrised to a smooth `|Vds|`, so
            ## the excess must be measured against `|dps|` for the whole
            ## factor to stay EVEN under the source/drain exchange.
            adps = _v(2.0 * hdl.hypsmooth(dps, 1e-6) - dps, 'ids_adps')
            inv_vp = _v(hdl.safe_div(1.0, mob['vp'], eps=1e-30), 'ids_ivp')
            excess = _v(hdl.hypsmooth(mob['vds'] - adps, 1e-3), 'ids_exc')
            s1n = _v(hdl.safe_ln(1.0 + excess * inv_vp), 'ids_s1n')
            if vdse is None:
                s1 = _v(s1n, 'ids_s1')
            else:
                ## `Vdse <= Vds` always, so the denominator's excess is
                ## the smaller of the two and `s1` stays non-negative --
                ## `dL` is a channel SHORTENING and a negative one would
                ## be meaningless.  Smoothed the same way as the
                ## numerator, so neither logarithm puts a kink in the
                ## middle of an I-V curve.
                excd = _v(hdl.hypsmooth(vdse - adps, 1e-3), 'ids_excd')
                s1d = _v(hdl.safe_ln(1.0 + excd * inv_vp), 'ids_s1d')
                s1 = _v(s1n - s1d, 'ids_s1')
            dL = _v(alp * s1, 'ids_dL')
            GdL = _v(hdl.safe_div(1.0, 1.0 + dL + dL * dL, eps=1e-30),
                     'ids_GdL')
        Gmob_dL = _v(Gmob * GdL, 'ids_GmobdL')

        ## Same modulation at the midpoint, with the MIDPOINT inversion
        ## charge and the source-side `xitsb` (`macrodefs:757-765`).
        tsg = (1.0 if xitsb is None
               else _wsat(qim, xitsb, mob.get('thesatg', 0.0), 'm'))
        thesat1 = _v(mob['thesat'] * tsg
                     * hdl.safe_div(1.0, Gmob_dL, eps=1e-30), 'ids_ths')
        zsat = _v(thesat1 * thesat1 * dps * dps, 'ids_zsat')
        if pmos:
            ## `macrodefs:766-771`.  `dps >= 0` by the terminal ordering
            ## and `thesat1 >= 0`, so the denominator is at least 1 and
            ## needs no guard.
            zsat = _v(zsat * hdl.safe_div(1.0, 1.0 + thesat1 * dps,
                                          eps=1e-30), 'ids_zsatp')
        Gvsat = _v(0.5 * (Gmob_dL * (1.0 + sympy.sqrt(1.0 + 2.0 * zsat))),
                   'ids_Gvsat')

    gvinv = _v(hdl.safe_div(1.0, Gvsat, eps=1e-30), 'ids_gvinv')

    ## `FdL = (1 + dL1 + dL1^2) * GdL` (`PSP103_module.include:1137`).
    ##
    ## This was left out on the reading that `ALP1` and `ALP2` are zero,
    ## which would make `dL1 = dL` and `FdL` exactly 1.  They are NOT:
    ## PSP's own `lp_alp2` is 4.5 on a 0.13 um device.  `ALP2` multiplies
    ## the BULK charge, so it dominates near threshold -- exactly where
    ## the output conductance was worst.  At Vg = 0.6 the reference
    ## current climbs 2.5x through saturation where the core without
    ## this climbed 1.4x.
    ##
    ## And that climb is NOT a threshold effect: measured on the
    ## reference, PSP's own `vth` moves only 3.5 mV over 1.35 V of drain
    ## bias.  DIBL is genuinely small here; this is the term that
    ## carries the output conductance.
    alp1 = 0.0 if mob is None else mob.get('alp1', 0.0)
    alp2 = 0.0 if mob is None else mob.get('alp2', 0.0)
    if (not isinstance(alp1, sympy.Expr)
            and not isinstance(alp2, sympy.Expr)
            and alp1 == 0.0 and alp2 == 0.0):
        FdL = _v(sympy.Integer(1), 'ids_FdL')
    else:
        inv_qim1 = _v(hdl.safe_div(1.0, qim1, eps=1e-30), 'ids_iq1')
        r1 = _v(qim * inv_qim1, 'ids_r1')
        r2 = _v(phit * alpha * inv_qim1, 'ids_r2')
        ## `s2` takes PSP's SOFTENED `Vdsx` (`module:1135`), where `s1`
        ## above takes the ordered `Vds`.  They are different arguments
        ## and PSP spells them differently; using `|Vds|` for both makes
        ## `s2` four times too large at `Vds = 0.05`.
        s2 = _v(hdl.safe_ln(1.0 + mob.get('vdsx', mob['vds'])
                            * hdl.safe_div(1.0, mob['vp'], eps=1e-30)),
                'ids_s2')
        dL1 = _v(dL + alp1 * (inv_qim1 * r1 * s1)
                 + alp2 * (qbm * r2 * r2 * s2), 'ids_dL1')
        FdL = _v((1.0 + dL1 + dL1 * dL1) * GdL, 'ids_FdL')

    ids = _v(beta * (FdL * qim1 * dps * gvinv), 'ids')
    return dict(ids=ids, qim=qim, qim1=qim1, alpha=alpha, dps=dps,
                xgm=xgm, x_m=x_m, x_s=x_s, x_d=x_d, qbm=qbm,
                Pm=Pm, Dm=Dm, sqm=sqm, Gmob_dL=Gmob_dL, Gvsat=Gvsat,
                zsat=zsat, gvinv=gvinv, GdL=GdL, eta_p=eta_p, FdL=FdL)


def ids_long_channel(xg, xn_s, xn_d, Gf, xi, phit, beta, margin=1e-5):
    """The intrinsic long-channel drain current -- see `intrinsic`."""
    return intrinsic(xg, xn_s, xn_d, Gf, xi, phit, beta, margin)['ids']


def cox_qm(core, cox, qq, phit, pmos=False):
    """PSP's quantum-mechanical correction to the OXIDE capacitance,
    for the CHARGE model (`PSP103_module.include:1474-1478`):

        COX_qm = COX / (1 + qq*(qeff1^2 + qlim2)^(-1/6))

    with `qeff1 = qbm + eta_mu1*qim` and `qlim2 = 100*phit^2`
    (`macrodefs:322`).

    The inversion layer does not sit at the interface -- it peaks a
    little way into the silicon -- which puts a capacitance in series
    with the oxide and reduces the total.  PSP already applies the same
    `qq` to `phib` and the body factor, where it shifts the threshold;
    this is the same physics acting on the charge instead.

    Worth about 12% here, and it was invisible for as long as the charge
    model was only checked against itself: conservation, source/drain
    symmetry and the Ward-Dutton partition are all RATIOS, and a uniform
    error in the oxide capacitance divides out of every one of them.
    """
    if not isinstance(qq, sympy.Expr) and qq == 0.0:
        return cox
    eta1 = ONE_THIRD if pmos else 0.5
    qeff1 = _v(core['qbm'] + eta1 * core['qim'], 'q_qeff1')
    qlim2 = _v(100.0 * phit * phit, 'q_qlim2')
    ## The base is bounded below by `qlim2 > 0`, so the negative power
    ## needs no guard -- which is what `qlim2` is there for.
    fac = _v(1.0 + qq * (qeff1 * qeff1 + qlim2) ** (-ONE_SIXTH), 'q_qmfac')
    return _v(cox * hdl.safe_div(1.0, fac, eps=1e-30), 'q_coxqm')


def charges_long_channel(core, xg, phit, cox):
    """The intrinsic terminal charges, PSP's way.

    `core` is what `intrinsic` returned.  `cox` is the TOTAL oxide
    capacitance, ``Cox' * W * L``.  Returns ``(Qg, Qd, Qb)``; the source
    charge is not returned because it is not independent -- contributing
    these three on branches referred to the source makes
    ``Qs = -(Qg + Qd + Qb)`` exactly, which is how charge conservation
    becomes structural rather than something to verify numerically.

    The drain/source split is PSP's ``SWQPART = 0``, the linear
    distribution -- the Ward-Dutton partition, which tends to the
    textbook 40/60 ratio in saturation.  (``SWQPART = 1`` assigns
    everything to the source; not implemented, since it is not the
    default and nothing here selects it.)

    CHANNEL-LENGTH MODULATION IS IN THE CHARGES TOO.  It used to be
    left out here -- ``GdL = 1`` and ``QCLM = 0`` -- as a long-channel
    simplification, and for a long time nothing could measure what that
    cost, because the charge model was only ever checked against its own
    construction properties.  Against PSP's terminal capacitances it is
    plain: in the linear region the two agree to 1 part in 10000 at both
    geometries, and the whole error is in SATURATION, growing with drain
    bias and with ``ALP`` -- 1% on the long device at `Vds = 1.2` and 8%
    on the short one.  That is `GdL`, and it is the same `GdL` the
    current already used.

    The polysilicon depletion factor ``eta_p`` is here too, on the same
    inversion bracket PSP puts it on -- again the same quantity the
    current already carried.
    """
    qim, alpha = core['qim'], core['alpha']
    qim1, dps, xgm = core['qim1'], core['dps'], core['xgm']

    ## Voxm is the oxide voltage at the midpoint.  `H` is PSP's
    ## effective charge slope: `Gmob/Gvsat * qim1 / alpha1`, where
    ## `alpha1` carries the velocity-saturation correction
    ## (macrodefs:776-778).  With mobility off both reduce to 1 and this
    ## is the ideal `qim1/alpha`.
    Voxm = _v(xgm * phit, 'q_Voxm')
    gr = _v(core['Gmob_dL'] * core['gvinv'], 'q_gr')
    alpha1 = _v(alpha * (1.0 + 0.5 * (core['zsat'] * gr * gr)), 'q_a1')
    H = _v(gr * qim1 * hdl.safe_div(1.0, alpha1, eps=1e-30), 'q_H')
    Fj = _v(0.5 * dps * hdl.safe_div(1.0, H, eps=1e-30), 'q_Fj')
    Fj2 = _v(Fj * Fj, 'q_Fj2')
    t6 = _v(alpha * dps * ONE_SIXTH, 'q_t6')

    ## Inverted: the general expressions.  Below flat band there is no
    ## inversion charge and the gate charge is the oxide voltage alone.
    ## `PSP103_module.include:1488-1500`.  With `GdL = 1` every one of
    ## these collapses to the long-channel form it replaces: `QCLM`
    ## vanishes and the `GdL` factors drop out.
    GdL = core.get('GdL', sympy.Integer(1))
    QCLM = _v((1.0 - GdL) * (qim - 0.5 * (alpha * dps)), 'q_QCLM')
    ## `eta_p` is the polysilicon depletion factor, and PSP applies it
    ## to the whole inversion bracket of the gate charge -- the same
    ## factor the current already carries, acting on the charge.
    eta_p = core.get('eta_p', sympy.Integer(1))
    QG_on = _v(Voxm + 0.5 * (eta_p * dps
                             * (Fj * GdL * ONE_THIRD - 1.0 + GdL)),
               'q_QGon')
    QD_on = _v(0.5 * (GdL * GdL * (qim - t6 * (1.0 - Fj - 0.2 * Fj2))
                      + QCLM * (1.0 + GdL)), 'q_QDon')
    QI_on = _v(GdL * (qim + t6 * Fj) + QCLM, 'q_QIon')

    inv = xg > 0.0
    QG = _v(sympy.Piecewise((QG_on, inv), (Voxm, True)), 'q_QG')
    QD = _v(sympy.Piecewise((QD_on, inv), (0.0, True)), 'q_QD')
    QI = _v(sympy.Piecewise((QI_on, inv), (0.0, True)), 'q_QI')
    QB = _v(QG - QI, 'q_QB')

    return (_v(QG * cox, 'q_Qg'), _v(-QD * cox, 'q_Qd'),
            _v(-QB * cox, 'q_Qb'))


def noise(core, mob, beta, phit, cox, nz):
    """PSP's channel noise: thermal and flicker, as power spectral
    densities on the channel branch.

    `PSP103_module.include:1819-1841`.  Both come out of the SAME
    surface-potential quantities the current does -- `qim`, `qim1`,
    `alpha`, `dps`, `H` -- which is the point of a surface-potential
    model's noise: there is nothing to fit separately.

    Returns ``(Sid, Sfl, shared)``: the white drain-current density in
    A^2/Hz, the flicker density at 1 Hz to be divided by ``f**EF``, and
    the channel-shape intermediates `induced_gate_noise` needs -- which
    are shared rather than recomputed because they are the SAME
    description of the channel, and two copies of it could drift apart.
    """
    qim, qim1 = core['qim'], core['qim1']
    alpha, dps = core['alpha'], core['dps']

    ## `H` is the effective charge slope, the same quantity the charge
    ## model builds (`macrodefs:776-778`); recomputed here from `core`
    ## rather than threaded through, since it is three cheap lines.
    gr = _v(core['Gmob_dL'] * core['gvinv'], 'n_gr')
    alpha1 = _v(alpha * (1.0 + 0.5 * (core['zsat'] * gr * gr)), 'n_a1')
    H = _v(gr * qim1 * hdl.safe_div(1.0, alpha1, eps=1e-30), 'n_H')

    ## ---- thermal ----------------------------------------------------
    H0 = _v(qim1 * hdl.safe_div(1.0, alpha, eps=1e-30), 'n_H0')
    t1 = _v(qim * hdl.safe_div(1.0, qim1, eps=1e-30), 'n_t1')
    sqt2 = _v(0.5 * ONE_SIXTH * (dps * hdl.safe_div(1.0, H0, eps=1e-30)),
              'n_sqt2')
    t2 = _v(sqt2 * sqt2, 'n_t2')
    r = _v(H0 * hdl.safe_div(1.0, H, eps=1e-30) - 1.0, 'n_r')
    ## PSP's own floor, and it is doing real work: `lc` is squared into
    ## a denominator, so a zero there is a division by zero rather than
    ## a large number.
    lc = _v(sympy.Max(1.0 - 12.0 * (r * t2), 1.0e-20), 'n_lc')
    g_ideal = _v(beta * (core['FdL'] * qim1 * core['gvinv']), 'n_gid')
    mid0 = _v(t1 + 12.0 * t2 - 24.0 * ((1.0 + t1) * t2 * r), 'n_mid0')
    mid = _v(sympy.Max(mid0, 1.0e-40), 'n_mid')
    ## `nt = 4*FNT*k*T`, written as `4*FNT*phit*q` so the thermal
    ## voltage already in hand carries the temperature and no separate
    ## `k` or `T` is needed.
    nt = _v(nz['fnt'] * 4.0 * phit * PSP_QELE, 'n_nt')
    sid = _v(nt * (g_ideal * mid
                   * hdl.safe_div(1.0, lc * lc, eps=1e-30)), 'n_sid')

    ## ---- flicker ----------------------------------------------------
    ## Number-fluctuation flicker: the carrier density at the two
    ## channel ends and its difference, in PSP's `Cox/q` units.
    coq = _v(cox * (1.0 / PSP_QELE), 'n_coq')
    n1 = _v(coq * (alpha * phit), 'n_N1')
    nm1 = _v(coq * qim1, 'n_Nm1')
    dn1 = _v(coq * (alpha * dps), 'n_dN1')
    ## The logarithm's argument is a ratio of the two channel ends'
    ## densities.  PSP guards neither, relying on `dN1 < 2*Nm1`; floored
    ## here because the DSL differentiates it and a vanishing
    ## denominator is a pole rather than a large number.
    lo = _v(sympy.Max(nm1 - 0.5 * dn1, 1.0e-30), 'n_lo')
    ratio = _v((nm1 + 0.5 * dn1) * hdl.safe_div(1.0, lo, eps=1e-30),
               'n_rat')
    sfl0 = _v((nz['nfa'] - nz['nfb'] * n1 + nz['nfc'] * (n1 * n1))
              * hdl.safe_ln(sympy.Max(ratio, 1.0e-30))
              + (nz['nfb'] + nz['nfc'] * (nm1 - 2.0 * n1)) * dn1,
              'n_sfl0')
    prefac = _v(phit * phit * beta * hdl.safe_div(1.0, coq, eps=1e-30),
                'n_pre')
    sfl = _v(sympy.Max(prefac * (core['ids'] * core['gvinv']) * sfl0
                       * hdl.safe_div(1.0, n1, eps=1e-30), 0.0),
             'n_sfl')
    return sid, sfl, dict(t1=t1, t2=t2, sqt2=sqt2, r=r, lc=lc,
                          g_ideal=g_ideal, nt=nt)


#: Lower bound on the auxiliary node's conductance, in siemens.
#:
#: The induced-gate network below hangs on an internal node whose only
#: DC connection is that conductance, so a zero there is a FLOATING NODE
#: and a singular matrix, not merely a large impedance.  The bound is
#: safe to state rather than hope about: the gate density it guards is
#: ``Sig = w^2 CGeff^2 g nt / (g_cond^2 + w^2 CGeff^2) <= nt * g``
#: whatever ``g_cond`` is, so wherever the floor is reached at all the
#: quantity it perturbs is already below ``nt * 1e-12`` -- of order
#: 1e-33 A^2/Hz, thirteen decades under any gate noise worth reporting.
GMIG_FLOOR = 1.0e-12


def induced_gate_noise(core, shared, sid, cox_c, swign):
    """PSP's induced gate noise, as an auxiliary RC node and a
    CORRELATED pair of sources (`PSP103_module.include:1863-1881,
    1942-1948`).

    The channel's charge fluctuation does not only leave through the
    drain.  It also couples capacitively to the gate, giving a gate
    current that grows as ``f`` and -- being the SAME fluctuation -- is
    correlated with the drain noise.  Two things follow, and PSP models
    both: the gate sees a density at all, and the drain's own density is
    split into a correlated part and the remaining ``(1 - c^2)``.

    Returns a dict for the element to stamp:

    ``pwr``    the shared source power, on the auxiliary branch and,
               scaled by ``migid``, on the channel branch,
    ``gcond``  that branch's conductance,
    ``cgeff``  the coupling capacitance to the gate,
    ``migid``  the channel branch's share of the shared source (NOT
               gated by ``swign`` -- ``pwr`` already is, and gating one
               factor of a product twice would cube the switch),
    ``c2``     ``c_igid**2``, by which the channel's INDEPENDENT source
               is reduced.

    The frequency dependence is deliberately NOT in any of them.  PSP
    builds it from a real network -- a conductance and a capacitance on
    an internal node -- so the ``f^2`` rise and the roll-off above
    ``1/(2 pi mig CGeff)`` come out of the circuit rather than out of a
    hand-written ``jw``.  Reproduced here as the same network, which is
    what makes `PSP103_module.include:1969-1970`'s 1 kHz reference value
    an EXACT check rather than a low-frequency one.
    """
    t1, t2, sqt2 = shared['t1'], shared['t2'], shared['sqt2']
    r, lc, g_ideal, nt = (shared['r'], shared['lc'], shared['g_ideal'],
                          shared['nt'])

    ## `mig` is a RESISTANCE here (PSP contributes `V(NOIR)/mig`), and it
    ## carries `1/g_ideal` -- so it runs to infinity as the channel
    ## empties.  Built as its RECIPROCAL throughout: the conductance is
    ## what the matrix wants, it is the quantity that stays finite, and
    ## the source power `nt/mig` is then a product rather than a
    ## quotient of two things that both approach zero.
    mig0 = _v(t1 * (1.0 / 12.0) - t2 * (t1 + 0.2 - 12.0 * t2)
              - 1.6 * (t2 * (t1 + 1.0 - 12.0 * t2) * r), 'g_mig0')
    mig_r = _v(sympy.Max(mig0, 1.0e-40), 'g_migr')
    gmig = _v(g_ideal * (lc * lc) * hdl.safe_div(1.0, mig_r, eps=1e-45),
              'g_gmig')
    pwr = _v(swign * (nt * gmig), 'g_pwr')
    gcond = _v(sympy.Max(gmig, GMIG_FLOOR), 'g_gcond')

    ## `CGeff = Gvsat^2 COX_qm eta_p / Gmob_dL^2` (`module:1868`), so
    ## the factor is `(Gvsat/Gmob_dL)^2` -- and since
    ## `Gvsat = 0.5*Gmob_dL*(1 + sqrt(1 + 2*zsat))` (`macrodefs:771`)
    ## that ratio is a function of `zsat` ALONE and `GdL` cancels out
    ## of it exactly.
    ##
    ## The key is `Gmob_dL` and not `Gmob` because it holds
    ## `Gmob * GdL`.  It USED to be spelled `Gmob`, and that cost a
    ## real bug: reading PSP's formula literally invites
    ## `core['Gmob'] * core['GdL'] * core['gvinv']`, which carries
    ## `GdL` twice and leaves a spurious `1/GdL^2` that grows with
    ## `Vds` through channel-length modulation -- indistinguishable
    ## from a velocity-saturation error.  A key that names its contents
    ## reads wrong at the one call site that is wrong.
    gr = _v(core['Gmob_dL'] * core['gvinv'], 'g_ratio')
    cgeff = _v(swign * cox_c * core['eta_p']
               * hdl.safe_div(1.0, gr * gr, eps=1e-30), 'g_cgeff')

    ## `c_igid = migid0 / sqrt(mig * mid)` with PSP's `mid = Sid / nt`
    ## (`module:1876-1878`).  In conductance form that square root has
    ## no small denominator left in it.
    migid0 = _v(hdl.safe_div(1.0, lc * lc, eps=1e-30) * sqt2
                * (1.0 - 12.0 * t2
                   - (t1 + 19.2 * t2 - 12.0 * (t1 * t2)) * r), 'g_mid0')
    ## PSP clips `c_igid` to [0, 1] and then rebuilds `migid` from it.
    ## The clip is a bound on `migid0` itself: unclipped the rebuild
    ## returns `migid0` exactly, and the two ends of the clip are 0 and
    ## `sqid/sqig = sqrt(Sid / (nt * gmig))`.  Written that way it needs
    ## no division by the correlation coefficient.
    sid_f = _v(sympy.Max(sid, 0.0), 'g_sidf')
    gm_f = _v(sympy.Max(gmig, 1.0e-100), 'g_gmigf')
    hi = _v(sympy.sqrt(hdl.safe_div(sid_f, nt * gm_f, eps=1e-150)),
            'g_hi')
    lo = _v(sympy.Min(migid0, hi), 'g_lo')
    migid = _v(sympy.Max(lo, 0.0), 'g_migid')

    ## The reduction of the channel's INDEPENDENT source.  PSP writes it
    ## as `1 - c_igid^2`; the identity `migid^2 * pwr = c_igid^2 * Sid`
    ## is exact by construction, so the factor is read straight off the
    ## correlated source rather than from a second copy of the
    ## coefficient -- which also makes the two sources sum to `Sid`
    ## identically, whatever the clip did.
    c2r = _v(migid * migid * pwr
             * hdl.safe_div(1.0, sid_f, eps=1e-150), 'g_c2r')
    c2 = _v(sympy.Min(c2r, 1.0), 'g_c2')
    return dict(pwr=pwr, gcond=gcond, cgeff=cgeff, migid=migid, c2=c2)
