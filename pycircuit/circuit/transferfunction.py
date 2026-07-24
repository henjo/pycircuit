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

    >>> import sympy
    >>> from pycircuit.circuit.transferfunction import TransferFunction
    >>> s = sympy.Symbol('s')
    >>> R, C = sympy.symbols('R C', positive=True)
    >>> H = TransferFunction(1, 1 + s*R*C, s)   # a first-order low-pass
    >>> H.canonical()
    1/(C*R*s + 1)
    >>> H.poles()
    {-1/(C*R): 1}
    >>> H.zeros()
    {}
    >>> H.dcgain()
    1
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

    def _ginac_coeffs(self):
        """``(num_coeffs, den_coeffs)`` in ``s`` (low order first) via GiNaC.

        Uses :func:`pycircuit.circuit._ginac.poly_coeffs`, which does the
        cancellation and polynomial collection natively in GiNaC (~4x faster than
        sympy on the collect step, and it cancels common factors like
        ``sympy.cancel``).

        Only used when the coefficients are *symbolic* (free parameters besides
        ``s``): there the values are exact, so the rationalisation in GiNaC is
        loss-free.  For a purely numeric transfer function the coefficients are
        floats that can span a huge dynamic range (a scaled solve produces terms
        from 1e0 down to ~1e-90); GiNaC's ``limit_denominator`` would round the
        tiny high-order ones to zero and silently drop poles, so this returns
        ``None`` and the caller uses sympy/``numpy``.  ``None`` is also returned
        if the extension is unavailable or the extraction fails.
        """
        if not (self.num.free_symbols | self.den.free_symbols) - {self.s}:
            return None                          # numeric: not safe for GiNaC
        try:
            from . import _ginac
            return _ginac.poly_coeffs(self.num / self.den, self.s)
        except Exception:
            return None

    def as_num_den(self):
        """Return ``(Poly(num, s), Poly(den, s))`` after cancellation.

        When the GiNaC extension is built the coefficients are collected
        natively (faster; see :meth:`_ginac_coeffs`); otherwise sympy's
        ``cancel``/``fraction`` is used.  Either way ``num``/``den`` are coprime.
        """
        coeffs = self._ginac_coeffs()
        if coeffs is not None:
            nums, dens = coeffs
            return (sympy.Poly(list(reversed(nums)), self.s),
                    sympy.Poly(list(reversed(dens)), self.s))
        num, den = sympy.fraction(self.canonical())
        return sympy.Poly(num, self.s), sympy.Poly(den, self.s)

    # -- poles / zeros / gain ----------------------------------------------
    @staticmethod
    def _roots(poly, numeric):
        """Roots of a sympy ``Poly`` in ``s``.

        With ``numeric=False`` return the symbolic ``sympy.roots`` mapping
        ``{root: multiplicity}``.  With ``numeric=True`` require the polynomial
        to have numeric coefficients and return a numpy array of complex roots
        via ``numpy.roots`` -- fast and reliable for high-degree denominators
        where the symbolic solver has no closed form.
        """
        if not numeric:
            return sympy.roots(poly)

        coeff_symbols = set().union(*[c.free_symbols for c in poly.all_coeffs()]) \
            if poly.all_coeffs() else set()
        if coeff_symbols:
            raise ValueError(
                "numeric root-finding needs numeric coefficients; substitute "
                "values for %s first" % ', '.join(map(str, sorted(coeff_symbols, key=str))))
        return np.roots([complex(c) for c in poly.all_coeffs()])

    def poles(self, numeric=False):
        """Return the poles (roots of the denominator).

        See :meth:`_roots` for the ``numeric`` flag.
        """
        _, den = self.as_num_den()
        return self._roots(den, numeric)

    def zeros(self, numeric=False):
        """Return the zeros (roots of the numerator).

        See :meth:`_roots` for the ``numeric`` flag.
        """
        num, _ = self.as_num_den()
        return self._roots(num, numeric)

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

    def frequencyresponse(self, freqs):
        """Return the complex response ``H(s=2*pi*j*f)`` over an array of Hz.

        Convenience wrapper over :meth:`bode`; requires numeric coefficients.
        """
        return self.bode()(2j * np.pi * np.asarray(freqs, dtype=float))

    def __repr__(self):
        return 'TransferFunction(%s, %s)' % (self.num, self.den)
