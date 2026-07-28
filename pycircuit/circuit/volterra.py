# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Volterra analysis based on the method described in
"Distortion analysis of analog integrated circuits" by Piet Wambacq and Willy Sansen

p129 chapter 5.2 depicts a flowchart describing the basic algorithm used in this module

"""

from numpy import array, delete, linalg, size, zeros, concatenate, pi, zeros, all, maximum
from pycircuit.circuit.analysis import Analysis, remove_row_col
from pycircuit.circuit.analysis_ss import AC
from pycircuit.circuit.dcanalysis import DC
from pycircuit.utilities import Parameter, ParameterDict, isiterable
from pycircuit.post.internalresult import InternalResultDict
from pycircuit.post import Waveform
import sympy
from sympy import Symbol, Matrix, symbols, simplify, together, factor, cancel, diff, Mul, factorial
#from symbolicelements import R, defaultepar, gnd, Diode, SubCircuit, IS, C, VCCS
from pycircuit.circuit.elements import R, Diode, IS, C, VCCS
from pycircuit.circuit.circuit import defaultepar, gnd, SubCircuit, Node, default_toolkit
from copy import copy
from pycircuit.circuit.toolkit import symbolic
from pycircuit.circuit.toolkit import numeric

class NLVCCS(VCCS):
    """Voltage-controlled current source with a cubic nonlinearity,
    ``I = gm*v + 2*v**2 + 3*v**3`` -- a minimal nonlinear element for
    exercising :func:`K`.
    """
    linear = False

    def i(self, x, epar=defaultepar):
        v = x[0]-x[1]
        I = self.ipar.gm * v + 2 * v**2 + 3 * v**3
        return array([0,0, I, -I])

    def G(self, x, epar=defaultepar):
        ## dI/dv, stamped directly: VCCS.update()/G() assume a *linear*
        ## element (I = gm*v), which would give the wrong (constant-gm)
        ## Jacobian for this cubic i() -- the same G-vs-di/dx mismatch
        ## test_element_jacobians.py polices for the shipped nonlinear
        ## elements. NLVCCS is test-only and not in that suite, but the
        ## invariant still has to hold for K() to differentiate correctly.
        v = x[0]-x[1]
        dIdv = self.ipar.gm + 4*v + 9*v**2
        return self.toolkit.array([[0, 0, 0, 0],
                                   [0, 0, 0, 0],
                                   [dIdv, -dIdv, 0, 0],
                                   [-dIdv, dIdv, 0, 0]])

    def update(self, subject):
        ## VCCS.update() precomputes and caches self._G from gm alone
        ## whenever an instance parameter changes -- correct for a linear
        ## VCCS, wrong for this one, whose Jacobian depends on x, not just
        ## parameters. G() above never reads self._G, so the cached value
        ## is never used; overriding as a no-op stops it from being computed
        ## (and crashing on a symbolic gm) in the first place.
        pass

def product(factors):
    return Mul(*factors)

def K(cir, x, ordervec, epar = defaultepar):
    """Calculate the taylor series term of the given order of i(x)

    Example:

    cir.K([1,1,0] will return the vector 1/(1!*1!) * d2I(X)/dX_0dX_1

    >>> from pycircuit.circuit.circuit import Node
    >>> default_toolkit = symbolic
    >>> d = NLVCCS(Node('plus'), Node('minus'), Node('plus'), Node('minus'), toolkit=symbolic, gm=Symbol('gm'))
    >>> epar = copy(defaultepar)
    >>> epar.T = Symbol('T')
    >>> K(d, [0,0], [1,0], epar=epar)
    array([0, 0, gm, -gm], dtype=object)
    >>> K(d, [0,0], [2,0], epar=epar)
    array([0, 0, 2, -2], dtype=object)
    >>> K(d, [0,0], [3,0], epar=epar)
    array([0, 0, 3, -3], dtype=object)
    >>> K(d, [0,0], [0,2], epar=epar)
    array([0, 0, 2, -2], dtype=object)

    """

    ## Generate a list of symbols like [x[0], x[1] ..]
    xsyms = array([Symbol('x[%d]'%i) for i in range(size(x,0))])

    ## Calculate derivative
    didx = cir.i(xsyms, epar=epar)
    for xsym, order in zip(xsyms, ordervec):
        if order > 0:
            ## sympy.diff on a plain array no longer treats a trailing
            ## integer as a derivative order (it did when this was written
            ## in 2008) -- apply it element-by-element instead.
            didx = array([diff(sympy.sympify(expr), xsym, order)
                         for expr in didx], dtype=object)

    ## Calculate taylor coefficient
    K = 1 / product(map(factorial, ordervec))

    ## Substitute x[k] with given x values
    return array([K * expr.subs(zip(xsyms, x)) for expr in didx])

class Volterra(Analysis):
    """
    Symbolic Volterra analysis class -- currently an unfinished stub.

    :func:`K` (above) is the complete, tested primitive: the Taylor-series
    term of a given order of a nonlinear element's ``i(x)``. This class was
    meant to orchestrate it across a whole circuit's nonlinear elements,
    following the flowchart in Wambacq & Sansen, *Distortion analysis of
    analog integrated circuits*, p129 ch. 5.2 -- but that orchestration was
    never written: :meth:`solve` finds the nonlinear elements and stops
    there (the call to :func:`K` was commented out, never finished, not
    later broken), and :meth:`run` returns an empty result. Completing it
    means implementing that algorithm, not fixing a bug; not attempted here
    without the source in hand -- see ``doc/architecture.md`` P14.

    The example below is deliberately linear-only. ``solve`` runs the
    circuit through ``AC(toolkit=symbolic)``, which -- like every symbolic
    toolkit -- never computes a DC operating point to linearise a nonlinear
    element around (see ``doc/architecture.md`` P11); a circuit containing a
    real nonlinear element such as :class:`NLVCCS` fails inside that ``AC``
    solve itself, before ``solve`` even reaches the point of listing
    nonlinear elements. That's a second, structural reason this class was
    never finished, not something fixed here.

    >>> cir = SubCircuit()
    >>> n1, n2 = cir.add_nodes('n1', 'n2')

    >>> cir['is'] = IS(n1, gnd, i=1e-3, iac=1)
    >>> cir['r1'] = R(n1, n2, r=Symbol('R1', positive=True))
    >>> cir['r2'] = R(n2, gnd, r=Symbol('R2', positive=True))
    >>> cir['c'] = C(n2, gnd, c=Symbol('C', positive=True))

    >>> volterra = Volterra(cir, toolkit=symbolic)
    >>> res = volterra.run()
    >>> len(res)
    0
    """
    numeric = False

    def run(self, **kvargs):
        x = self.solve(**kvargs)

        result = InternalResultDict()

        self.result = result
        return result

    def solve(self, refnode=gnd):
        x = zeros(self.cir.n)

        ac = AC(self.cir, toolkit=symbolic)
        xac = ac.solve(freqs = Symbol('s'), refnode = refnode, complexfreq = True)

        ## Find non-linear elements
        nlelements = [e for e in self.cir.xflatelements if not e.linear]
        return nlelements

if __name__ == "__main__":
    import doctest
    doctest.testmod()
