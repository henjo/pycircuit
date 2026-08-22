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


def mina(x, y, a):
    """PSP's ``MINA`` -- a smooth minimum.

    ``0.5*(x + y - sqrt((x-y)^2 + a))``.  Branch-free, and the radicand
    is bounded below by ``a > 0``, so nothing needs guarding.
    """
    x, y = sympy.sympify(x), sympy.sympify(y)
    return 0.5 * (x + y - sympy.sqrt((x - y) ** 2 + a))


def _sigma_body(a, ab, c, tau, eta):
    """Shared body of ``sigma`` and ``sigma2``.

    They differ only in whether the ``a`` in two places is ``a`` or
    ``a*b``; PSP ships them as separate macros, and writing the shared
    algebra once keeps them provably the same shape.
    """
    nu = a + c
    mutau = nu * nu + (0.5 * c * c - ab) * tau
    ## `mutau` appears in a denominator twice.  It is not guaranteed
    ## nonzero by anything in the algebra, and a zero here would put a
    ## NaN into the surface potential -- and so into every current the
    ## model computes.  Regularised rather than guarded, because a guard
    ## would still evaluate the division in its untaken arm.
    inv_mutau = hdl.safe_div(1.0, mutau, eps=1e-30)
    denom = mutau + (nu * tau * tau * inv_mutau) * c * (c * c * ONE_THIRD
                                                        - ab)
    return eta + a * nu * tau * hdl.safe_div(1.0, denom, eps=1e-30)


def sigma(a, c, tau, eta):
    """PSP's ``sigma`` -- the three-argument correction."""
    a, c, tau, eta = map(sympy.sympify, (a, c, tau, eta))
    return _sigma_body(a, a, c, tau, eta)


def sigma2(a, b, c, tau, eta):
    """PSP's ``sigma2`` -- the four-argument variant.

    The vendor guards ``|tau| < 1e-120`` and returns ``eta``, because
    "sometimes tau is extremely small".  Kept, as a `Piecewise`: with
    `safe_div` underneath, the other arm is finite there anyway, but the
    guard is what the reference model computes and the point is to match
    it.
    """
    a, b, c, tau, eta = map(sympy.sympify, (a, b, c, tau, eta))
    return sympy.Piecewise((eta, sympy.Abs(tau) < 1e-120),
                           (_sigma_body(a, a * b, c, tau, eta), True))


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
    y0 = _v(sigma(a_a, c_a, tau_a, eta_a), 'y0')
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
    x0 = _v(sigma2(a_i, b_i, c_i, tau_i, eta_i), 'x0')

    ## The vendor splits `delta0`/`delta1` three ways to keep both
    ## exp(x0) and exp(x0 - xn) inside range; `expl` does the same job in
    ## one expression and is C-3 rather than piecewise-continuous.
    ex0 = _v(hdl.expl(x0), 'ex0')
    d0_i = _v(delta * ex0, 'd0_i')
    d1_i = _v(hdl.safe_div(1.0, ex0, eps=1e-30), 'd1_i')

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
