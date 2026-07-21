# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Symbolic transfer-function post-processing.

:class:`TransferFunction` wraps a symbolic ``H(s) = numerator/denominator`` and
provides the operations an engineer wants from a symbolic small-signal result:
polynomial numerator/denominator, poles, zeros, DC gain, and a fast numeric
evaluator via :func:`sympy.lambdify`.

It is produced by :meth:`~pycircuit.circuit.analysis_ss.CircuitResultACPoly.tf`
when an analysis is run with the ``symbolic_poly`` toolkit, but it is a plain
value class and can be constructed directly from any sympy expressions.
"""

import numpy as np
import sympy


class TransferFunction:
    """A symbolic transfer function ``H(s) = num/den``.

    ``num`` and ``den`` are sympy expressions in the frequency variable ``s``
    (and possibly other symbols, treated as parameters).  They are not assumed
    to be coprime; :meth:`canonical`, :meth:`poles` and :meth:`zeros` cancel
    common factors as needed.
    """

    def __init__(self, num, den, s):
        self.num = sympy.sympify(num)
        self.den = sympy.sympify(den)
        self.s = s

    # -- expression forms ---------------------------------------------------
    def expr(self):
        """Return the raw ``num/den`` expression (no cancellation)."""
        return self.num / self.den

    def canonical(self):
        """Return ``num/den`` with common factors cancelled."""
        return sympy.cancel(self.num / self.den)

    def as_num_den(self):
        """Return ``(Poly(num, s), Poly(den, s))`` after cancellation."""
        num, den = sympy.fraction(self.canonical())
        return sympy.Poly(num, self.s), sympy.Poly(den, self.s)

    # -- poles / zeros / gain ----------------------------------------------
    def poles(self):
        """Return the poles as ``{root: multiplicity}`` (roots of the denominator)."""
        _, den = self.as_num_den()
        return sympy.roots(den)

    def zeros(self):
        """Return the zeros as ``{root: multiplicity}`` (roots of the numerator)."""
        num, _ = self.as_num_den()
        return sympy.roots(num)

    def dcgain(self):
        """Return the DC gain ``H(s=0)``."""
        return sympy.limit(self.expr(), self.s, 0)

    # -- numeric evaluation -------------------------------------------------
    def bode(self, modules='numpy'):
        """Return a fast callable ``f(s)`` for the transfer function.

        The expression must have no free symbols other than ``s`` (substitute
        parameter values first).  Evaluate on the ``s = j*omega`` axis for a
        frequency response, e.g. ``f(2j*np.pi*freqs)``.
        """
        return sympy.lambdify(self.s, self.canonical(), modules)

    def __repr__(self):
        return 'TransferFunction(%s, %s)' % (self.num, self.den)
