# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""An `ACSolution` backed by determinant decision diagrams.

Two ways out of the same solution, and the distinction matters:

**The graph path** answers questions without ever building an expression --
poles from the s-expanded coefficients, transfer functions evaluated by walking
the diagram.  This works at sizes where an explicit ``N(s)/D(s)`` does not exist
in any useful sense.

**The compatibility path** hands back sympy expressions so that everything
already written against ``num/den`` keeps working.  That means *flattening*,
which is the one thing the representation exists to avoid, so it is guarded by
product-term count and refuses rather than hanging.  It is for small circuits.

Anything that needs a symbolic answer for a large circuit wants neither of
these directly: it wants the approximation stage, where the expression is pruned
to its dominant terms before being written down.
"""

import numpy as np
import sympy

from .acsolution import ACSolution
from .ddd import DDDSizeError, ddd_cramer, ddd_of_matrix, s_expand
from .transferfunction import TransferFunction


__all__ = ['DDDSolution']


class DDDSolution(ACSolution):
    """AC solution held as diagrams, expanded only when explicitly asked.

    Args:
        A: System matrix ``s*C + G``, reference node already removed.
        b: Right-hand side, likewise reduced.
        s: Frequency symbol.
        irefnode: Index at which the reference node is re-inserted; its voltage
            is zero, so its numerator is zero.
        order: Expansion ordering passed through to the diagram builder.
        max_terms: Guard for any operation that flattens to sympy.
    """

    def __init__(self, A, b, s, irefnode, order='min-degree', max_terms=5000):
        super().__init__(s)
        self._A = sympy.Matrix(A)
        self._b = sympy.Matrix(b).reshape(self._A.rows, 1)
        self.irefnode = irefnode
        self.order = order
        self.max_terms = max_terms

        self._den = None
        self._nums = {}
        self._sexp = None

    ## -- index bookkeeping ------------------------------------------------

    @property
    def n(self):
        """Number of unknowns *including* the reference node."""
        return self._A.rows + 1

    def _reduced(self, index):
        """Map a full index to the reduced system, or None for the ref node."""
        if index == self.irefnode:
            return None
        return index if index < self.irefnode else index - 1

    ## -- the graph path (nothing is expanded) -----------------------------

    def denominator_diagram(self):
        """The network determinant, as a diagram."""
        if self._den is None:
            self._den = ddd_of_matrix(self._A, order=self.order)
        return self._den

    def numerator_diagram(self, index):
        """Cramer numerator for one unknown, as a diagram (None for ref node)."""
        k = self._reduced(index)
        if k is None:
            return None
        if k not in self._nums:
            _, nums = ddd_cramer(self._A, self._b, indices=[k], order=self.order)
            self._nums[k] = nums[k]
        return self._nums[k]

    def s_expanded(self):
        """Denominator split by powers of ``s``, as shared diagrams."""
        if self._sexp is None:
            self._sexp = s_expand(self._A, self.s, order=self.order)
        return self._sexp

    def eval_node(self, index, env):
        """Numeric value of one unknown, by walking the diagrams."""
        num = self.numerator_diagram(index)
        if num is None:
            return 0.0 + 0.0j
        return complex(num.eval(env)) / complex(self.denominator_diagram().eval(env))

    def eval_tf(self, out_index, in_index, env):
        """Numeric ``x[out]/x[in]``, without forming the determinant.

        The shared denominator cancels in the ratio, so this needs only the two
        numerator diagrams.
        """
        num_out = self.numerator_diagram(out_index)
        num_in = self.numerator_diagram(in_index)
        if num_out is None:
            return 0.0 + 0.0j
        if num_in is None:
            raise ZeroDivisionError('input unknown is the reference node')
        return complex(num_out.eval(env)) / complex(num_in.eval(env))

    ## -- the compatibility path (expands; guarded) ------------------------

    def _flatten(self, diagram, what):
        try:
            return diagram.to_sympy(max_terms=self.max_terms)
        except DDDSizeError as exc:
            raise DDDSizeError(
                '%s cannot be written as a sympy expression: %s.  Use the '
                'diagram API instead -- poles() works without expanding, and '
                'eval_tf()/eval_node() evaluate numerically.' % (what, exc))

    def numerators(self):
        """Numerator per unknown as sympy expressions (expands -- guarded)."""
        out = []
        for p in range(self.n):
            d = self.numerator_diagram(p)
            out.append(sympy.Integer(0) if d is None
                       else self._flatten(d, 'numerator %d' % p))
        return np.array(out, dtype=object)

    def denominator(self):
        """The determinant as a sympy expression (expands -- guarded)."""
        return self._flatten(self.denominator_diagram(), 'the determinant')

    def node_voltages(self):
        den = self.denominator()
        return np.array([n / den for n in self.numerators()], dtype=object)

    def transfer_function(self, numerator):
        return TransferFunction(numerator, self.denominator(), self.s)

    def poles(self, numeric=False):
        """Circuit poles -- roots of the determinant.

        With ``numeric=True`` this uses the **graph path**: the s-expanded
        coefficients are evaluated and handed to ``numpy.roots``, so the
        determinant is never expanded and this works on circuits far beyond what
        :meth:`denominator` could return.

        With ``numeric=False`` a symbolic polynomial has to be formed, which
        expands and is therefore guarded.

        Note:
            Root-finding from expanded coefficients is ill-conditioned when they
            span many orders of magnitude, as they do for any sizeable circuit --
            a property of the polynomial form, not of the diagram.  See
            ``doc/src/circuit/ddd.rst``.
        """
        expanded = self.s_expanded()
        if numeric:
            return expanded.roots_of()

        coeffs = expanded.eval_coeffs()
        total = sum(expanded.coefficient(k).term_count()
                    for k in range(len(coeffs)))
        if total > self.max_terms:
            raise DDDSizeError(
                'the symbolic denominator would have %d product terms, above '
                'the %d limit; use poles(numeric=True), which works on the '
                'diagram without expanding' % (total, self.max_terms))
        poly = sympy.Poly(sum(c * self.s ** k for k, c in enumerate(coeffs)),
                          self.s)
        return TransferFunction._roots(poly, numeric)

    def __repr__(self):
        return '<DDDSolution n=%d order=%s>' % (self.n, self.order)
