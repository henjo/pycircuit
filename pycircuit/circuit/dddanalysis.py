# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The DDD family of analyses.

These are the ordinary analyses computed on determinant decision diagrams
instead of expanded expressions.  They exist as **separate entry points** so the
family is discoverable as a family, but they are deliberately *thin*: none of
them reimplements the analysis it is named after.

That is possible because the existing analyses already ask the toolkit for the
quantity a diagram is best at --
:class:`~pycircuit.circuit.nportanalysis.TwoPortAnalysis` reaches its parameters
through ``toolkit.det`` and ``toolkit.cofactor``, and
:class:`~pycircuit.circuit.feedback.FeedbackLoopAnalysis` forms the loop gain as
``det(Y) / det(Y_noloop)``.  Overriding those two primitives on
:class:`~pycircuit.circuit.toolkit.DDDToolkit` therefore moves the whole family
onto diagrams without a line of analysis logic being rewritten.

Duplicating the port bookkeeping or the return-ratio derivation here would leave
two implementations to keep in step, and the originals would stop being a usable
reference for checking these against.  So each class below binds the toolkit and
adds the diagram-level accessors the base class has no concept of -- and nothing
else.

The originals are untouched and remain the reference::

    from pycircuit.circuit.dddanalysis import DDDTwoPort
    from pycircuit.circuit.nportanalysis import TwoPortAnalysis

    fast = DDDTwoPort(cir, inp, inn, outp, outn).solve(freqs, complexfreq=True)
    ref = TwoPortAnalysis(cir, inp, inn, outp, outn).solve(freqs, complexfreq=True)
"""

import sympy

from .analysis_ss import TransimpedanceAnalysis
from .feedback import FeedbackLoopAnalysis
from .nportanalysis import TwoPortAnalysis
from .toolkit import ddd_toolkit


__all__ = ['DDDTwoPort', 'DDDLoopGain', 'DDDTransimpedance', 'DDDAnalysisMixin']


class DDDAnalysisMixin:
    """Binds `ddd_toolkit` and exposes the diagram behind the result.

    Everything an analysis does is inherited; this only chooses the
    representation and opens a window onto it.
    """

    def __init__(self, *args, **kvargs):
        kvargs.setdefault('toolkit', ddd_toolkit)
        super().__init__(*args, **kvargs)

    def family(self, A):
        """The shared `DDDFamily` the toolkit built for matrix ``A``.

        Useful for asking what the analysis cost -- vertex counts, how many
        product terms the answer stands for -- which the base class has no way
        to express.
        """
        return self.toolkit._ddd_family(A)

    def determinant_diagram(self, A):
        """The network determinant of ``A`` as a diagram, unexpanded."""
        return self.family(A).denominator


class DDDTwoPort(DDDAnalysisMixin, TwoPortAnalysis):
    """Two-port parameters computed on diagrams.

    Y, Z and ABCD parameters are ratios of determinants and cofactors of a
    single matrix, and the toolkit serves them all from one shared family --
    so the whole parameter set costs one construction rather than one per
    parameter.

    Examples:
        >>> from pycircuit.circuit.dddanalysis import DDDTwoPort   # doctest: +SKIP
        >>> res = DDDTwoPort(cir, n1, gnd, n2, gnd).solve(freqs)   # doctest: +SKIP
    """


class DDDLoopGain(DDDAnalysisMixin, FeedbackLoopAnalysis):
    """Loop gain and return difference computed on diagrams.

    The return difference is ``det(Y) / det(Y_noloop)``, and where the matrix
    depends linearly on the loop-forming parameter the second determinant is a
    *cofactor* of the first -- see
    :func:`~pycircuit.circuit.ddd.determinant_sensitivity`.  Both determinants
    go through the toolkit, so this needs no logic of its own.
    """

    def loop_sensitivity(self, A, parameter):
        """``d det(A)/d parameter`` as a combination of cofactors.

        The quantity a return ratio is built from, exposed directly: with
        ``det A = det A|_{k=0} + k · d det A/dk``, this is the derivative term,
        and it comes from the family the determinant already built.
        """
        from .ddd import determinant_sensitivity
        combination, _ = determinant_sensitivity(sympy.Matrix(A), parameter,
                                                 family=self.family(A))
        return combination


class DDDTransimpedance(DDDAnalysisMixin, TransimpedanceAnalysis):
    """Transimpedances computed on diagrams.

    Every transimpedance is a Cramer numerator over the shared determinant, and
    :func:`~pycircuit.circuit.ddd.ddd_cofactor_solve` produces the whole vector
    from one construction -- which is also what makes the reduced noise path
    affordable.
    """
