# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Toolkit classes: ``numeric``, ``symbolic`` and experimental ``symbolic_poly``.

A *toolkit* provides the math primitives (``linearsolver``, ``dot``, ``zeros``,
``conj``, ``det``, the elementary functions, physical constants, ...) that the
analyses run through, so the same analysis code works numerically or
symbolically.

Each toolkit is a thin class wrapping a backend module of primitive functions
(:mod:`pycircuit.circuit._numeric`, :mod:`pycircuit.circuit._symbolic`).  The
class adds capability flags (``symbolic``, ``poly``), an overridable method
surface, and capability introspection via :meth:`Toolkit.supports`.  Anything
not defined on the class is delegated to the backend module.

The module-level singletons ``numeric``, ``symbolic`` and ``symbolic_poly`` are
the objects passed as ``toolkit=`` to analyses and stored as
``circuit.default_toolkit``.  They are drop-in replacements for the old toolkit
*modules* -- ``toolkit.zeros(...)``, ``toolkit.symbolic`` etc. keep working.
"""

import numpy as np
import sympy
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.exceptions import DMError

from . import _numeric
from . import _symbolic


class Toolkit:
    """Base toolkit -- delegates primitive operations to a backend module."""

    symbolic = False
    poly = False

    def __init__(self, backend):
        self._backend = backend

    def __getattr__(self, name):
        ## Primitives (dot, zeros, cos, pi, kboltzmann, ...) live on the backend
        ## module.  Guard names starting with '_' to avoid recursing on
        ## self._backend before __init__ has run.
        if name.startswith('_'):
            raise AttributeError(name)
        return getattr(self._backend, name)

    def supports(self, capability):
        """Return True if this toolkit implements an optional capability.

        Optional capabilities let analyses opt into richer operations when the
        toolkit provides them, without assuming every toolkit does.
        """
        return False

    def noise_psd(self, Y, u, CY, s):
        """Return ``(transimpedance_vector, output noise PSD)``.

        Solves the reciprocal system ``Y zm = -u`` and forms the noise power
        ``zm^T CY conj(zm)``.  ``s`` is unused in this generic implementation but
        is part of the interface so toolkits that need it can override (see
        :class:`SymbolicPolyToolkit`).
        """
        Ym = self.toMatrix(Y)
        um = self.toMatrix(u)
        zm = self.linearsolver(Ym, -um)
        psd = self.dot(self.dot(zm.reshape(1, self.size(zm)), CY),
                       self.conj(zm))
        return zm, psd[0]


class NumericToolkit(Toolkit):
    """Numeric toolkit backed by numpy."""
    symbolic = False


class SymbolicToolkit(Toolkit):
    """Symbolic toolkit backed by sympy (stock LUsolve solver)."""
    symbolic = True


class SymbolicPolyToolkit(SymbolicToolkit):
    """Experimental symbolic toolkit using polynomial-domain linear algebra.

    Behaves like :class:`SymbolicToolkit` but solves linear systems
    fraction-free over sympy's polynomial/rational domain
    (``DomainMatrix.solve_den``) rather than with ``Matrix.LUsolve``.  This
    keeps intermediate entries as polynomials instead of nested fractions and
    avoids the expression swell that makes LUsolve blow up on larger circuits.
    """
    poly = True

    def supports(self, capability):
        return capability in ('num_den',)

    def linearsolver_num_den(self, A, b):
        """Solve ``A x = b`` returning ``(numerator_vector, shared_denominator)``.

        The solution is ``x_i = numerator[i] / denominator``.  In the polynomial
        domain the shared denominator is the network determinant, computed
        fraction-free so intermediate entries stay polynomials rather than the
        nested fractions ``LUsolve`` produces.  This keeps the transfer function
        as a single ``N(s)/D(s)`` without a swelling ``cancel``.

        Falls back to the stock symbolic ``LUsolve`` when the entries are not in
        a polynomial/rational domain (transcendental functions land in sympy's
        ``EX`` domain); the denominator is then 1 and the numerators are the
        LUsolve solution directly.
        """
        A = sympy.Matrix(A)
        b = sympy.Matrix(b.tolist())

        if A.shape == (1, 1):
            return np.array([b[0]], dtype=object), A[0, 0]

        try:
            dM = DomainMatrix.from_Matrix(A)
            db = DomainMatrix.from_Matrix(b)
            domain = dM.domain.unify(db.domain)
            ## Fraction-free elimination needs an *exact* domain.  Float values
            ## land in an inexact RR[...] domain where solve_den is
            ## pathologically slow (and does not raise, so it would hang rather
            ## than fall back); convert to the exact domain (QQ[...]) first.
            if not domain.is_Exact:
                domain = domain.get_exact()
            xnum, den = dM.convert_to(domain).solve_den(db.convert_to(domain))
            xnum = xnum.to_Matrix()
            den = domain.to_sympy(den)
            return (np.array([xnum[i, 0] for i in range(xnum.rows)], dtype=object),
                    den)
        except (DMError, sympy.PolynomialError, TypeError, AttributeError):
            x = np.array(self._backend.linearsolver(A, b), dtype=object)
            return x.reshape((np.size(x, 0),)), sympy.Integer(1)

    def linearsolver(self, A, b):
        """Solve ``A x = b`` fraction-free; returns the divided solution vector.

        Thin wrapper over :meth:`linearsolver_num_den` (``x_i = numer_i / den``)
        so the toolkit honours the standard solution-vector contract used by the
        analyses.
        """
        num, den = self.linearsolver_num_den(A, b)
        return np.array([ni / den for ni in num], dtype=object)

    def noise_psd(self, Y, u, CY, s):
        """Shared-denominator noise PSD.

        With ``zm = N/D`` (fraction-free), the noise power
        ``zm^T CY conj(zm)`` equals ``(N^T CY conj(N)) / (D conj(D))`` -- a
        single rational with a polynomial numerator over ``|D|^2`` rather than
        the O(n^2) sum of divided rationals the generic form builds.  Value is
        identical; only the intermediate is kept compact.

        Conjugation is ``sympy.conjugate``, which resolves to ``N(-s)``
        automatically when ``s`` is imaginary (``s = j*omega``).
        """
        N, D = self.linearsolver_num_den(self.toMatrix(Y), -self.toMatrix(u))
        N = sympy.Matrix(list(N))
        CYm = sympy.Matrix(np.asarray(CY).tolist())
        num = (N.T * CYm * N.applyfunc(sympy.conjugate))[0]
        den = D * sympy.conjugate(D)
        zm = np.array([ni / D for ni in N], dtype=object)
        return zm, num / den


## Singletons -- drop-in replacements for the old toolkit modules.
numeric = NumericToolkit(_numeric)
symbolic = SymbolicToolkit(_symbolic)
symbolic_poly = SymbolicPolyToolkit(_symbolic)
