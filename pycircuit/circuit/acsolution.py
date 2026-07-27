# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Representation-independent view of a solved AC system.

An AC analysis can end up holding its answer in more than one form.  The
fraction-free solver produces a numerator vector over a shared denominator; a
determinant decision diagram produces graphs that are never expanded at all.
Both can answer the same questions -- transfer functions, poles, node voltages --
but they answer them by different means.

An `ACSolution` is that common interface.  The analysis asks its toolkit for one
and hands it to the result object, so adding a representation does not add a
branch to the analysis: without this, every new backend means another
``if toolkit.supports(...)`` beside another result class, which is how analysis
code accumulates cruft.

Layering note: this module deliberately imports nothing from ``circuit`` or
``analysis`` so that ``toolkit`` can import it.  The result classes live above
both and hold an `ACSolution`; nothing here reaches upward.
"""

from abc import ABC, abstractmethod

import numpy as np
import sympy

from .transferfunction import TransferFunction


__all__ = ['ACSolution', 'NumDenSolution']


class ACSolution(ABC):
    """A solved AC system, in whatever representation the toolkit produced.

    Implementations are expected to be cheap to construct and lazy about
    expanding anything -- for large symbolic circuits an explicit ``N(s)/D(s)``
    may be astronomically large, so it must not be formed just to build a result.

    Attributes:
        s: The complex frequency symbol.
    """

    def __init__(self, s):
        self.s = s

    @abstractmethod
    def numerators(self):
        """Numerator per unknown, indexed like ``x`` (reference node included)."""

    @abstractmethod
    def denominator(self):
        """The shared denominator -- the network determinant."""

    @abstractmethod
    def node_voltages(self):
        """The solution vector ``x``, one entry per unknown."""

    @abstractmethod
    def transfer_function(self, numerator):
        """Wrap ``numerator`` over the shared denominator."""

    @abstractmethod
    def poles(self, numeric=False):
        """Roots of the denominator."""


class NumDenSolution(ACSolution):
    """The fraction-free form: a numerator vector over one shared denominator.

    This is what :class:`~pycircuit.circuit.toolkit.SymbolicPolyToolkit` and the
    GiNaC backend produce, and it is exactly the behaviour the AC analysis had
    before the representation was made pluggable.
    """

    def __init__(self, num, den, s):
        super().__init__(s)
        self._num = num
        self._den = den

    def numerators(self):
        return self._num

    def denominator(self):
        return self._den

    def node_voltages(self):
        ## Object dtype because the entries are sympy expressions, not numbers.
        return np.array([n / self._den for n in self._num], dtype=object)

    def transfer_function(self, numerator):
        return TransferFunction(numerator, self._den, self.s)

    def poles(self, numeric=False):
        return TransferFunction._roots(sympy.Poly(self._den, self.s), numeric)

    def __repr__(self):
        return '<NumDenSolution n=%d>' % len(self._num)
