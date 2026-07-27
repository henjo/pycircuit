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
from . import _sparse_numeric
from . import _symbolic
from . import _jaxtoolkit
from .acsolution import NumDenSolution


class Toolkit:
    """Base toolkit -- delegates primitive operations to a backend module."""

    symbolic = False
    poly = False
    jax = False

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

    def ac_solution(self, A, b, s, irefnode):
        """Return an :class:`~pycircuit.circuit.acsolution.ACSolution`, or None.

        A toolkit that can express the AC answer in a compact shared form
        overrides this; returning ``None`` means "no compact form available",
        and the analysis falls back to solving frequency by frequency.

        Keeping the choice here rather than in the analysis is what lets a new
        representation be added without another branch in ``analysis_ss``.

        Args:
            A: The system matrix ``s*C + G`` with the reference node removed.
            b: Right-hand side, likewise reduced.
            s: The frequency symbol.
            irefnode: Index at which to re-insert the reference node, whose
                voltage is zero by definition.
        """
        return None

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

class SparseNumericToolkit(NumericToolkit):
    """Numeric toolkit backed by scipy.sparse."""
    symbolic = False
    poly = False
    sparse = True
    
    def build_sparse(self, data, rows, cols, shape):
        from scipy.sparse import coo_matrix
        return coo_matrix((data, (rows, cols)), shape=shape)

class JAXToolkit(Toolkit):
    """Numeric toolkit backed by JAX for auto-differentiation."""
    symbolic = False
    poly = False
    jax = True
    
    def generate_batched_eval(self, element_cls, method='i'):
        import jax
        
        cache_attr = f'_jax_batched_{method}'
        if not hasattr(element_cls, cache_attr):
            setattr(element_cls, cache_attr, {})
            
        cache = getattr(element_cls, cache_attr)
        
        def batched_eval_func(X_batch, params_batch, epar):
            if hasattr(epar, '_values'):
                key = frozenset(epar._values.items())
            else:
                key = id(epar)
                
            if key not in cache:
                @jax.jit
                def compiled_func(X_b, p_b):
                    def single_eval(x, p):
                        if method == 'i':
                            return element_cls.eval_i_pure(x, p, epar, self)
                        elif method == 'q':
                            return element_cls.eval_q_pure(x, p, epar, self)
                    
                    # --- HIGH-PERFORMANCE GPU/CPU VECTORIZATION ---
                    # Instead of evaluating thousands of Diodes/Transistors in a Python `for` loop,
                    # we use JAX's `vmap` (Vectorizing Map). It takes our `single_eval` function 
                    # (which operates on a single device) and maps it across a massive stacked array 
                    # of inputs (X_b) and parameters (p_b). 
                    # 
                    # Simultaneously, `jax.jacfwd` performs Forward-Mode Automatic Differentiation 
                    # to EXACTLY calculate the Jacobian (the Conductance / Capacitance matrix) for 
                    # every single device without manually writing calculus formulas.
                    # 
                    # JAX compiles this entire block via `@jax.jit` into a highly optimized 
                    # XLA executable that executes simultaneously in C++/GPU space.
                    val_batch = jax.vmap(single_eval)(X_b, p_b)
                    jac_batch = jax.vmap(jax.jacfwd(single_eval))(X_b, p_b)
                    return val_batch, jac_batch
                    
                cache[key] = compiled_func
                
            return cache[key](X_batch, params_batch)
            
        return batched_eval_func
        
    def generate_eval_i_and_G(self, element):
        import jax
        
        element._jax_cache_i = {}
        
        def eval_i_and_G_func(x, epar):
            if hasattr(epar, '_values'):
                key = frozenset(epar._values.items())
            else:
                key = id(epar)
                
            if key not in element._jax_cache_i:
                @jax.jit
                def compiled_func(x_inner):
                    i_vec = element.eval_i(x_inner, epar)
                    G_mat = jax.jacfwd(element.eval_i)(x_inner, epar)
                    return i_vec, G_mat
                element._jax_cache_i[key] = compiled_func
                
            return element._jax_cache_i[key](x)
            
        element.eval_i_and_G = eval_i_and_G_func

    def generate_eval_q_and_C(self, element):
        import jax
        
        element._jax_cache_q = {}
        
        def eval_q_and_C_func(x, epar):
            if hasattr(epar, '_values'):
                key = frozenset(epar._values.items())
            else:
                key = id(epar)
                
            if key not in element._jax_cache_q:
                @jax.jit
                def compiled_func(x_inner):
                    q_vec = element.eval_q(x_inner, epar)
                    C_mat = jax.jacfwd(element.eval_q)(x_inner, epar)
                    return q_vec, C_mat
                element._jax_cache_q[key] = compiled_func
                
            return element._jax_cache_q[key](x)
            
        element.eval_q_and_C = eval_q_and_C_func

    def add_at(self, lhs, indices, rhs):
        """Perform in-place addition for JAX arrays using .at[].add()"""
        # JAX arrays are immutable, so we must return the new array.
        # circuit.py does `lhs = toolkit.add_at(lhs, indices, rhs)`
        return lhs.at[indices].add(rhs)

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

    An AC solve with this toolkit returns a transfer-function-capable result
    (see :doc:`/circuit/symbolic_poly`)::

        import sympy
        from pycircuit.circuit import symbolic_poly, use_toolkit, SubCircuit, R, C, VS, gnd
        from pycircuit.circuit.analysis_ss import AC

        s = sympy.Symbol('s', imaginary=True)
        with use_toolkit(symbolic_poly):
            cir = SubCircuit()
            cir['VS'] = VS('in', gnd, vac=1)
            cir['R'] = R('in', 'out', r=sympy.Symbol('R', positive=True))
            cir['C'] = C('out', gnd, c=sympy.Symbol('C', positive=True))

        res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)
        res.tf('out', gnd).poles()      # -> {-1/(C*R): 1}
    """
    poly = True

    def supports(self, capability):
        return capability in ('num_den',)

    def ac_solution(self, A, b, s, irefnode):
        """Solve fraction-free and keep the result as ``N(s)/D(s)``."""
        num, den = self.linearsolver_num_den(A, b)
        ## Re-insert the reference node: v(refnode) = 0, so its numerator is 0.
        num = self.concatenate((num[:irefnode], self.array([0]), num[irefnode:]))
        return NumDenSolution(num, den, s)

    def linearsolver_num_den(self, A, b):
        """Solve ``A x = b`` fraction-free; returns ``(numerator_vector, denominator)``.

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


class SymengineToolkit(SymbolicPolyToolkit):
    """EXPERIMENTAL symbolic toolkit that solves through the compiled ``symengine``
    CAS, keeping the same ``N(s)/D(s)`` contract.  Correct, but **not a speed
    win** as implemented -- kept as a plug-in point for a fraction-free compiled
    backend (see ``doc/ginac_toolkit_design.md``).

    Behaves like :class:`SymbolicPolyToolkit` (assembly stays in sympy) but
    :meth:`linearsolver_num_den` solves in symengine.  The catch: symengine's
    ``LUsolve``/``det`` are *not* fraction-free -- they leave nested divisions
    that swell downstream -- so for symbolic ``N(s)/D(s)`` extraction this is
    slower than sympy's polynomial-domain ``DomainMatrix.solve_den`` (which
    :class:`SymbolicPolyToolkit` uses).  A real compiled win needs proper
    polynomial-GCD cancellation (GiNaC's ``normal()``); symengine's ``FFLU``
    fraction-free primitive does not cancel to a compact determinant here.
    Falls back to the sympy fraction-free path when symengine cannot handle the
    entries.
    """

    def linearsolver_num_den(self, A, b):
        import symengine as se
        A = sympy.Matrix(A)
        b = sympy.Matrix(b.tolist())

        if A.shape == (1, 1):
            return np.array([b[0]], dtype=object), A[0, 0]

        ## symengine drops sympy assumptions (e.g. Symbol('s', imaginary=True))
        ## on the round-trip, which would leave the result's symbols as different
        ## objects than the analysis uses (breaking poles()/conjugation).  Map the
        ## returned symbols back to the originals by name.
        orig_syms = {sym.name: sym for sym in (A.free_symbols | b.free_symbols)}

        def to_sympy(se_expr):
            e = sympy.sympify(se_expr)
            subs = {t: orig_syms[t.name] for t in e.free_symbols
                    if t.name in orig_syms}
            return e.xreplace(subs) if subs else e

        try:
            n = A.rows
            Ase = se.Matrix(n, n, [se.sympify(A[i, j])
                                   for i in range(n) for j in range(n)])
            bse = se.Matrix(n, 1, [se.sympify(b[i]) for i in range(n)])
            ## Shared denominator = network determinant; numerators = det * x
            ## (symengine cancels the det it introduced, leaving polynomials).
            den = Ase.det()
            xse = Ase.LUsolve(bse)
            num = np.array([to_sympy(xse[i] * den) for i in range(n)],
                           dtype=object)
            return num, to_sympy(den)
        except Exception:
            ## Any symengine limitation -> stock sympy fraction-free solve.
            return super().linearsolver_num_den(A, b)


class GinacToolkit(SymbolicPolyToolkit):
    """EXPERIMENTAL symbolic toolkit that solves through the compiled GiNaC
    extension.  Correct, but **not a general speed win** -- kept for evaluation.

    GiNaC's ``normal()`` performs true fraction-free cancellation (unlike
    symengine).  On a numeric-component AC system (entries ``G + s C``) the
    coefficients are first scaled to O(1) (:func:`._ginac._frequency_scale`),
    which avoids the CLN integer blow-up that otherwise stalls the solve around
    dimension 16 -- in that regime the scaled GiNaC solve is fast to high order
    (e.g. dim 24 in ~0.4 s) and beats sympy's ``DomainMatrix.solve_den``.  When
    the system cannot be scaled to O(1) (mixed numeric/symbolic entries that are
    not a single-frequency ``G + s C``), the :attr:`ginac_max_dim` guard still
    applies as a safety net; fully symbolic circuits do not blow up.

    Requires the built ``_ginac_ext`` extension (GiNaC + pybind11); the
    ``ginac_toolkit`` singleton is ``None`` when it is unavailable.  Any
    per-solve failure (parse error, singular, transcendental entries) falls back
    to sympy.
    """

    #: Fallback guard: above this dimension, a system that *cannot* be scaled to
    #: O(1) coefficients falls back to sympy to avoid the CLN blow-up.  Scalable
    #: (numeric-AC) systems are not subject to it.  Below the observed cliff ~16.
    ginac_max_dim = 12

    def linearsolver_num_den(self, A, b):
        try:
            from . import _ginac
            Am = sympy.Matrix(A)
            scalable = _ginac.detect_freq(Am) is not None
            if not scalable and Am.rows > self.ginac_max_dim:
                return super().linearsolver_num_den(A, b)
            return _ginac.linearsolver_num_den(Am, b)
        except Exception:
            return super().linearsolver_num_den(A, b)

    def eval_sweep(self, expr, params):
        """Evaluate ``expr`` over a parameter sweep via native GiNaC code.

        "Derive once, evaluate many": compiles the (possibly complex-valued)
        expression to machine code once and evaluates it per sweep point.  See
        :func:`._ginac.eval_sweep`; guarded so an oversized expression falls
        back to sympy rather than invoking the compiler on a huge source.
        """
        from . import _ginac
        return _ginac.eval_sweep(expr, params)

    def solve_native(self, A, b):
        """Solve ``A x = b`` returning a native :class:`._ginac.GinacResult`.

        Keeps the determinant-sized ``num[]``/``den`` as GiNaC ``ex`` and
        converts only the compact piece requested (transfer function,
        denominator, numeric sweep), avoiding the sympy round-trip that
        dominates a large symbolic solve.
        """
        from . import _ginac
        return _ginac.solve_native(A, b)


class DDDToolkit(SymbolicPolyToolkit):
    """Symbolic toolkit that solves via determinant decision diagrams.

    Stamping is unchanged -- circuits are built and stamped exactly as for any
    other symbolic toolkit.  What differs is the *representation* of the answer:
    instead of an expanded ``N(s)/D(s)``, the solution is a shared graph that is
    evaluated or split by powers of ``s`` without ever being written out.

    Use it like any other toolkit::

        from pycircuit.circuit.toolkit import ddd_toolkit
        res = AC(cir, toolkit=ddd_toolkit).solve(s, complexfreq=True)
        res.poles(numeric=True)      # works without expanding the determinant

    The compatibility surface (``linearsolver_num_den``, ``res.tf(...)``) has to
    expand, so it is limited by ``ddd_max_terms`` and raises rather than hanging.
    That limit is not a defect: for a fully symbolic circuit of any size the
    expanded form genuinely is unusable, which is the reason this toolkit exists.
    """

    #: Refuse to expand a diagram into sympy beyond this many product terms.
    ddd_max_terms = 5000

    def supports(self, capability):
        return capability in ('num_den', 'ddd')

    def ac_solution(self, A, b, s, irefnode):
        from .dddresult import DDDSolution
        return DDDSolution(A, b, s, irefnode, max_terms=self.ddd_max_terms)

    def noise_psd(self, Y, u, CY, s):
        """Noise PSD with every transimpedance sharing one construction.

        Noise analysis needs a transfer function from *each* noise source to the
        output -- the whole transimpedance vector ``z = Y^-1 (-u)``.  Building
        those as separate Cramer determinants re-expands minors that the network
        determinant already covered.  Expanding along the substituted column
        instead expresses each numerator through cofactors of ``Y`` itself, so
        the entire family comes out of one memo table.

        That is Shi & Tan's observation (TCAD 2001) that noise costs little more
        than a single transfer function.  Measured on an RC ladder at ``n = 15``:
        124 shared vertices for the determinant *and* all fifteen
        transimpedances, against 459 built separately, where the determinant
        alone is 40.

        The value is identical to the fraction-free form -- the shared
        denominator makes the PSD independent of any overall sign convention --
        so this is a construction saving, not a different answer.
        """
        from .ddd import DDDSizeError, ddd_cofactor_solve

        Ym = sympy.Matrix(self.toMatrix(Y))
        um = sympy.Matrix(self.toMatrix(u))
        family, nums = ddd_cofactor_solve(Ym, -um)

        def flatten(obj, what):
            n = obj.term_count()
            if n > self.ddd_max_terms:
                raise DDDSizeError(
                    '%s would expand to %d product terms, above the %d limit'
                    % (what, n, self.ddd_max_terms))
            return obj.eval()

        D = flatten(family.denominator, 'the determinant')
        N = sympy.Matrix([flatten(nums[i], 'transimpedance %d' % i)
                          for i in range(Ym.rows)])
        CYm = sympy.Matrix(np.asarray(CY).tolist())
        num = (N.T * CYm * N.applyfunc(sympy.conjugate))[0]
        den = D * sympy.conjugate(D)
        zm = np.array([ni / D for ni in N], dtype=object)
        return zm, num / den

    def linearsolver_num_den(self, A, b):
        """Solve via diagrams, then expand -- the compatibility path.

        Raises:
            DDDSizeError: If the expansion would exceed ``ddd_max_terms``.  The
                diagram API has no such limit; see :class:`DDDSolution`.
        """
        from .ddd import ddd_cramer
        Am = sympy.Matrix(A)
        den, nums = ddd_cramer(Am, b, order='min-degree')
        num = [nums[i].to_sympy(max_terms=self.ddd_max_terms)
               for i in range(Am.rows)]
        return np.array(num, dtype=object), den.to_sympy(max_terms=self.ddd_max_terms)



## Singletons -- drop-in replacements for the old toolkit modules.
numeric = NumericToolkit(_numeric)
sparse_numeric = SparseNumericToolkit(_sparse_numeric)
symbolic = SymbolicToolkit(_symbolic)
symbolic_poly = SymbolicPolyToolkit(_symbolic)
ddd_toolkit = DDDToolkit(_symbolic)

try:
    import symengine as _symengine
    symengine_toolkit = SymengineToolkit(_symbolic)
except ImportError:
    symengine_toolkit = None

try:
    from pycircuit.circuit import _ginac_ext as _ginac_ext_probe
    ginac_toolkit = GinacToolkit(_symbolic)
except ImportError:
    ginac_toolkit = None

try:
    import pycircuit.circuit._jaxtoolkit as _jaxtoolkit
    jaxtoolkit = JAXToolkit(_jaxtoolkit)
except ImportError:
    jaxtoolkit = None
