# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Experimental symbolic toolkit built on the polynomial domain.

This toolkit behaves like :mod:`pycircuit.circuit.symbolic` but replaces the
matrix solver (and, over time, other operations) with implementations that
work in sympy's polynomial/rational *domains* (``DomainMatrix``) instead of on
free-form expressions.  The goal is aggressive performance work -- fraction-free
elimination, exact-domain arithmetic, shared-denominator assembly, lambdify
bridges -- without changing the behaviour of the stock ``symbolic`` toolkit, so
existing code and forks that depend on it are unaffected.

Select it explicitly, e.g.::

    from pycircuit.circuit import symbolic_poly
    AC(cir, toolkit=symbolic_poly).solve(sympy.Symbol('s'), complexfreq=True)

Everything not overridden here is inherited unchanged from ``symbolic``.
"""

import numpy as np
import sympy
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.exceptions import DMError

## Inherit the whole symbolic toolkit surface (flags, dot, zeros, array, det,
## conj, constants, ...); override only the aggressive pieces below.
from .symbolic import *          # noqa: F401,F403

## Marker so analyses/tests can detect this variant (symbolic is also True).
poly = True


def linearsolver(A, b):
    """Solve A x = b fraction-free over an exact polynomial/rational domain.

    Uses ``DomainMatrix.solve_den`` (Bareiss-style, fraction-free) so
    intermediate entries stay polynomials rather than the nested fractions
    ``Matrix.LUsolve`` produces -- which avoids the expression swell that makes
    LUsolve blow up on larger circuits.  The solution is returned as a numerator
    vector over a shared denominator (``x_i = numer_i / den``); the element-wise
    division is a cheap ``Mul``, not an expansion, so no swell is added.
    Canonicalisation (cancel / Poly / lambdify) is left to the caller.

    Falls back to the stock ``symbolic.linearsolver`` (LUsolve) when the entries
    do not live in a polynomial/rational domain -- e.g. transcendental functions
    land in sympy's ``EX`` domain, where the fraction-free path gives no benefit.
    """
    A = sympy.Matrix(A)
    b = sympy.Matrix(b.tolist())

    if A.shape == (1, 1):
        return np.array([(b[0] / A[0, 0])])

    try:
        dM = DomainMatrix.from_Matrix(A)
        db = DomainMatrix.from_Matrix(b)
        domain = dM.domain.unify(db.domain)
        ## Fraction-free elimination needs an *exact* domain.  Float component
        ## values land in an inexact RR[...]/CC[...] domain where solve_den is
        ## pathologically slow (and does not raise, so it would hang rather than
        ## fall back).  Convert to the corresponding exact domain (QQ[...]) first;
        ## this is both fast and exact.
        if not domain.is_Exact:
            domain = domain.get_exact()
        xnum, den = dM.convert_to(domain).solve_den(db.convert_to(domain))
        xnum = xnum.to_Matrix()
        den = domain.to_sympy(den)
        return np.array([xnum[i, 0] / den for i in range(xnum.rows)],
                        dtype=object)
    except (DMError, sympy.PolynomialError, TypeError, AttributeError):
        res = np.array((A.LUsolve(b)))
        return res.reshape((np.size(res, 0),))
