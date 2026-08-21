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
from .acsolution import NumDenSolution

## _jaxtoolkit is deliberately NOT imported here: it does `import jax` at module
## level, so importing it eagerly made JAX a hard dependency of the whole
## package -- the numeric and symbolic paths included, which never touch it.
## The guarded import at the bottom of this file is what makes it optional, and
## an unconditional import here silently defeated that guard.  See
## test_jax_optional.


class Toolkit:
    """Base toolkit -- delegates primitive operations to a backend module."""

    symbolic = False
    poly = False
    jax = False
    sparse = False

    def __init__(self, backend):
        self._backend = backend

    def __getattr__(self, name):
        ## Primitives (dot, zeros, cos, pi, kboltzmann, ...) live on the backend
        ## module.  Guard names starting with '_' to avoid recursing on
        ## self._backend before __init__ has run.
        if name.startswith('_'):
            raise AttributeError(name)
        try:
            value = getattr(self._backend, name)
        except AttributeError:
            ## Plain delegation would raise here naming the *backend module*
            ## (e.g. "module 'pycircuit.circuit._numeric' has no attribute
            ## 'sparse'"), sending a reader to the wrong file -- the toolkit
            ## class, not the backend, is what a caller actually named.
            raise AttributeError(
                "%r toolkit (backend %r) has no attribute %r"
                % (type(self).__name__, self._backend.__name__, name)) from None
        else:
            ## MEMOISE.  `__getattr__` is only consulted when normal attribute
            ## lookup fails, so writing the resolved value into the instance
            ## __dict__ means every later access resolves directly and never
            ## enters this method.  Measured 3.77 million calls in one transient
            ## benchmark run before this line existed.
            ##
            ## Only SUCCESSES are cached.  Caching failures would be the larger
            ## win -- `hasattr(toolkit, 'add_at')` on a toolkit without it is the
            ## expensive path, since it formats the message below and raises --
            ## but a negative cache would also make an attribute that appears
            ## later permanently invisible.  The hot negative probes were hoisted
            ## out of the assembly loops in stage 2b instead, which is the fix
            ## that carries no such risk.
            ##
            ## Safe against runtime mutation of the *toolkit*: assigning
            ## `toolkit.foo = ...` writes the instance __dict__ directly and so
            ## overrides any memo.  It would NOT be safe against mutating the
            ## backend module after first access; nothing in the tree does that,
            ## and it is checked rather than assumed (see the note in
            ## doc/transient_work_plan.md, item 2+.1).
            object.__setattr__(self, name, value)
            return value

    def supports(self, capability):
        """Return True if this toolkit implements an optional capability.

        Optional capabilities let analyses opt into richer operations when the
        toolkit provides them, without assuming every toolkit does.
        """
        return False

    def jacobian(self, func, x, params, epar, fallback=None):
        """Jacobian of ``func`` with respect to ``x`` -- or ``fallback``.

        An element that knows its own conductance or capacitance matrix stamps
        it directly; one differentiated automatically does not need to.  Both
        answers come from here so the *element* never has to know which kind of
        toolkit it is running under.

        Before this existed, every such element carried its own
        ``if hasattr(toolkit, 'jax') and toolkit.jax: import jax; ...`` branch --
        twenty of them in ``elements.py`` alone, each importing jax inside a
        method, which is how an optional accelerator came to be named all over
        the element layer.  Adding a differentiating backend should not mean
        touching twenty elements again.

        Args:
            func: Pure ``f(x, params, epar, toolkit)`` to differentiate.
            x: Point to differentiate at.
            params: Element parameters passed through to ``func``.
            epar: Environment parameters passed through to ``func``.
            fallback: The precomputed matrix to return when this toolkit does
                not differentiate.  A plain attribute on the element, so
                evaluating it eagerly costs nothing.  Elements whose
                non-differentiating path is a *computation* (device limiting,
                a piecewise switch model) cannot pass it eagerly -- they guard
                with ``supports('autodiff')`` instead and leave this unset.
        """
        return fallback

    def derivative(self, func, at, eps=1e-6):
        """``d func / d at`` for a scalar, user-supplied function.

        Behavioural sources take an arbitrary Python callable, so there is no
        precomputed matrix to fall back on -- the derivative has to be computed.
        This does it by central difference; a differentiating toolkit overrides
        with exact differentiation.

        The elements that need this used to write ``try: import jax ... except
        ImportError:`` inline, which chose the method by whether JAX was
        *installed* rather than by which toolkit was running.  With JAX present,
        a numeric circuit silently got JAX arrays stamped into its matrix.
        Routing through the toolkit makes the active backend the thing that
        decides, which is the only correct answer.
        """
        return (func(at + eps) - func(at - eps)) / (2 * eps)

    def nth_derivative(self, func, at, order):
        """``d^order func / d at^order`` for a scalar function of a scalar.

        Distortion analysis needs the *second and third* Taylor coefficients of
        a device characteristic, where every other analysis needs only the
        first (the Jacobian).  Rather than ask elements to declare higher
        derivatives -- a second statement of the device model, free to drift
        from ``i()``, which is the defect ``test_element_jacobians`` exists to
        catch -- they are derived from the same pure function the element
        already exposes.  See ``doc/distortion_plan.md`` section 2.3.

        The base implementation is a central difference applied ``order``
        times, with the step widened for higher orders: differencing amplifies
        round-off as ``eps**-order``, so the naive 1e-6 that is right for a
        first derivative is hopeless for a third.  A symbolic backend
        overrides this and is exact at every order.
        """
        ## eps ~ macheps**(1/(order+2)) balances truncation against round-off
        ## for a central difference; the exponent is the standard result.
        eps = 1e-16 ** (1.0 / (order + 2))

        def diff_once(f):
            return lambda v: (f(v + eps) - f(v - eps)) / (2 * eps)

        d = func
        for _ in range(order):
            d = diff_once(d)
        return d(at)

    def evaluation_groups(self, circuit):
        """Element classes this toolkit can evaluate in a batch -- none, here.

        A toolkit that evaluates elements one at a time has no use for the
        grouping, so the empty mapping means "stamp every element individually"
        and the circuit needs no branch on toolkit type.
        """
        return {}

    def batched_contributions(self, circuit, groups, methodname, x, args,
                              params_tree, ndim):
        """Stamp contributions for grouped elements, or None if not batching.

        Returning ``None`` -- the base behaviour -- tells the circuit that
        nothing was handled in bulk and every element must be visited normally.

        Args:
            circuit: The circuit being stamped.
            groups: The mapping from :meth:`evaluation_groups`.
            methodname: ``'G'``, ``'C'``, ``'i'`` or ``'q'``.
            x: State vector.
            args: Extra arguments the element method takes (``epar`` first).
            params_tree: Optional parameter override tree.
            ndim: 2 for a matrix stamp, 1 for a vector stamp.

        Returns:
            ``(values, indices)`` where ``indices`` is ``(rows, cols)`` for
            ``ndim=2`` and a flat index array for ``ndim=1``; or ``None``.
        """
        return None

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


    def matrix_from_entries(self, shape, entries):
        """A matrix from sparse stamp entries: ``[(row, col, value), ...]``.

        THE ELEMENT-STAMP CONSTRUCTOR (owner request, after the VCCS stamp
        broke twice in one day): an element that builds its stamp by
        mutating ``toolkit.zeros(...)`` embeds a backend assumption --
        in-place assignment crashes on immutable JAX arrays, and a numpy
        float intermediate rejects symbolic values -- so the construction
        belongs to the TOOLKIT, which knows its own array semantics.  The
        base implementation accumulates into nested Python lists (which
        hold floats, sympy expressions, and tracers alike) and coerces once
        through this toolkit's own ``array()``; a toolkit with a better
        native route is free to override.
        """
        n_rows, n_cols = shape
        rows = [[0] * n_cols for _ in range(n_rows)]
        for r, c, v in entries:
            rows[r][c] = rows[r][c] + v
        return self.array(rows)


class NumericToolkit(Toolkit):
    def matrix_from_entries(self, shape, entries):
        ## dtype pinned to float: the base nested-list build lets all-integer
        ## stamps (gm=1) come out int64, which changed doctest reprs and
        ## invites integer-division surprises downstream.  The symbolic
        ## toolkit keeps the base implementation -- exact integer zeros and
        ## live sympy entries are exactly what it needs.
        import numpy as _np
        n_rows, n_cols = shape
        M = _np.zeros((n_rows, n_cols))
        for r, c, v in entries:
            M[r, c] += v
        return M

    """Numeric toolkit backed by numpy."""
    symbolic = False

    def jacobian(self, func, x, params, epar, fallback=None):
        """``fallback`` if the element has one, else a central difference.

        An element that stamps its own conductance matrix still wins -- this
        only computes anything for elements that have no hand-written form,
        which is what :class:`~pycircuit.circuit.semiconductors.Semiconductor`
        relies on.  The difference used to be hand-rolled inside that class;
        having it here means there is one implementation to be wrong rather
        than one per element family.
        """
        if fallback is not None:
            return fallback

        eps = 1e-6
        f0 = func(x, params, epar, self)
        n, m = len(x), len(f0)
        J = self.zeros((m, n))
        for j in range(n):
            x_plus, x_minus = list(x), list(x)
            x_plus[j] += eps
            x_minus[j] -= eps
            f_plus = func(self.array(x_plus), params, epar, self)
            f_minus = func(self.array(x_minus), params, epar, self)
            for k in range(m):
                J[k, j] = (f_plus[k] - f_minus[k]) / (2 * eps)
        return J

class SparseNumericToolkit(NumericToolkit):
    """Numeric toolkit backed by scipy.sparse."""
    symbolic = False
    poly = False
    sparse = True
    
    def build_sparse(self, data, rows, cols, shape):
        from scipy.sparse import coo_matrix
        return coo_matrix((data, (rows, cols)), shape=shape)

class JAXToolkit(Toolkit):
    def matrix_from_entries(self, shape, entries):
        ## Same float-dtype pin as NumericToolkit, built immutably: entries
        ## may be tracers, so accumulate in a Python grid and convert once.
        import jax.numpy as _jnp
        n_rows, n_cols = shape
        rows = [[0.0] * n_cols for _ in range(n_rows)]
        for r, c, v in entries:
            rows[r][c] = rows[r][c] + v
        return _jnp.asarray(rows, dtype=_jnp.float64)

    """Numeric toolkit backed by JAX for auto-differentiation."""
    symbolic = False
    poly = False
    jax = True

    def supports(self, capability):
        ## 'autodiff' means: derivatives come from differentiating the pure
        ## element functions, so an element must hand over its smooth form and
        ## skip the manual limiting a Newton solver would otherwise apply.
        return capability in ('autodiff',)

    def jacobian(self, func, x, params, epar, fallback=None):
        """Differentiate ``func`` at ``x`` -- the fallback matrix is unused."""
        import jax
        return jax.jacfwd(func)(x, params, epar, self)

    def derivative(self, func, at, eps=None):
        """Exact scalar derivative; ``eps`` is accepted and ignored."""
        import jax
        return jax.grad(func)(at)

    def evaluation_groups(self, circuit):
        """Group a circuit's elements by class, for batched evaluation.

        Elements of one class differ only in their node map and parameters, so
        they can be evaluated in a single ``vmap`` rather than a Python loop.
        Building the grouping is meaningless for a toolkit that cannot batch,
        which is why the base returns nothing and this lives here rather than in
        :mod:`pycircuit.circuit.circuit`.
        """
        import numpy as np
        groups = {}
        for instance_name, element in circuit.elements.items():
            cls = element.__class__
            if not (hasattr(cls, 'eval_i_pure') or hasattr(cls, 'eval_q_pure')):
                continue
            group = groups.setdefault(cls, {
                'instances': [], 'nodemaps_2d': [], 'params': {},
                'rows': [], 'cols': [], '1d_indices': []})
            group['instances'].append(element)
            group['nodemaps_2d'].append(circuit._map_indices_1d[instance_name])
            rows, cols = circuit._map_indices_2d[instance_name]
            group['rows'].append(rows)
            group['cols'].append(cols)
            group['1d_indices'].append(circuit._map_indices_1d[instance_name])
            values = getattr(getattr(element, 'iparv', None), '_values', {}) or {}
            for key, value in values.items():
                ## STAGE 10.3 -- a parameter that is not a number cannot be
                ## batched, and is not a model parameter either.
                ##
                ## `ic` (initial conditions, added for uic=True) defaults to
                ## None, which no `eval_*_pure` reads and which `self.array`
                ## cannot build: JAX raises "None is not a valid value for
                ## jnp.array" and takes three vectorised-stamping tests with it.
                ## Skipping the KEY rather than the value keeps the batch
                ## rectangular -- dropping only the None entries would leave a
                ## column shorter than the instance list and misalign every
                ## parameter after it.
                if value is None:
                    group['params'].pop(key, None)
                    group.setdefault('_unbatchable', set()).add(key)
                    continue
                if key in group.get('_unbatchable', ()):
                    continue
                group['params'].setdefault(key, []).append(value)

        for group in groups.values():
            group['nodemaps_2d'] = np.array(group['nodemaps_2d'], dtype=int)
            group['rows'] = np.array(group['rows']).flatten()
            group['cols'] = np.array(group['cols']).flatten()
            group['1d_indices'] = np.array(group['1d_indices']).flatten()
            for key in group['params']:
                group['params'][key] = self.array(group['params'][key])
        return groups

    def generate_batched_eval(self, element_cls, method='i'):
        """Compile a vmapped value-and-Jacobian evaluator for one element class.

        Instead of evaluating thousands of diodes in a Python loop, ``vmap``
        maps the pure single-device function across a stacked array of inputs
        and parameters, while ``jacfwd`` differentiates it to get each device's
        conductance or capacitance matrix.  ``jit`` compiles the pair into one
        XLA executable.

        An element class supplies whichever pure forms it needs -- a diode has
        no charge, a capacitor no static current -- so the missing one
        contributes **zeros** rather than raising.  Without that, any circuit
        mixing a nonlinear element with a reactive one failed on every stamp,
        because the two land in different groups and each stamp asks every group
        for both forms.
        """
        import jax
        import jax.numpy as jnp

        cache_attr = '_jax_batched_%s' % method
        if not hasattr(element_cls, cache_attr):
            setattr(element_cls, cache_attr, {})
        cache = getattr(element_cls, cache_attr)

        wanted = 'eval_%s_pure' % ('i' if method == 'i' else 'q')
        other = 'eval_%s_pure' % ('q' if method == 'i' else 'i')

        def batched_eval_func(X_batch, params_batch, epar):
            if hasattr(epar, '_values'):
                key = frozenset(epar._values.items())
            else:
                key = id(epar)

            if key not in cache:
                @jax.jit
                def compiled_func(X_b, p_b):
                    def single_eval(x, p):
                        func = getattr(element_cls, wanted, None)
                        if func is not None:
                            return func(x, p, epar, self)
                        ## This class has no such form; contribute nothing, but
                        ## with the shape the other form defines so the stamp
                        ## still lines up.
                        shape_of = getattr(element_cls, other)
                        return jnp.zeros_like(shape_of(x, p, epar, self))

                    return (jax.vmap(single_eval)(X_b, p_b),
                            jax.vmap(jax.jacfwd(single_eval))(X_b, p_b))

                cache[key] = compiled_func
            return cache[key](X_batch, params_batch)

        return batched_eval_func

    def batched_contributions(self, circuit, groups, methodname, x, args,
                              params_tree, ndim):
        """Evaluate every grouped class in one vmapped call per class."""
        if not groups or x is None:
            return None

        method_key = 'i' if methodname in ('G', 'i') else 'q'
        want_jacobian = methodname in ('G', 'C')
        epar = args[0]

        values, rows, cols, flat = [], [], [], []
        for cls, group in groups.items():
            evaluate = self.generate_batched_eval(cls, method_key)
            params = group['params']
            if params_tree is not None and cls.__name__ in params_tree:
                params = params_tree[cls.__name__]
            value_batch, jacobian_batch = evaluate(
                x[group['nodemaps_2d']], params, epar)
            batch = jacobian_batch if want_jacobian else value_batch
            values.append(self.reshape(
                batch, (len(group['instances']), -1)).flatten())
            if ndim == 2:
                rows.append(group['rows'])
                cols.append(group['cols'])
            else:
                flat.append(group['1d_indices'])

        if not values:
            return None
        joined = self.concatenate(values) if len(values) > 1 else values[0]
        if ndim == 2:
            return joined, (np.concatenate(rows), np.concatenate(cols))
        return joined, np.concatenate(flat)

    def add_at(self, lhs, indices, rhs):
        """Perform in-place addition for JAX arrays using .at[].add()"""
        # JAX arrays are immutable, so we must return the new array.
        # circuit.py does `lhs = toolkit.add_at(lhs, indices, rhs)`
        return lhs.at[indices].add(rhs)

class SymbolicToolkit(Toolkit):
    """Symbolic toolkit backed by sympy (stock LUsolve solver)."""
    symbolic = True

    def jacobian(self, func, x, params, epar, fallback=None):
        """``fallback`` if the element has one, else an *exact* sympy derivative.

        Without this the base class returned ``fallback`` -- ``None`` for an
        element whose whole design is to be differentiated automatically -- and
        the element was left to do it itself.  ``semiconductors.py`` did, by
        central difference, which under a symbolic toolkit perturbs a *symbol*
        by ``1e-6`` and divides by ``2e-6``: a numerical approximation of an
        expression that sympy can differentiate in closed form, simultaneously
        enormous and inexact.  Differentiating here is both exact and cheaper.

        ``x`` may hold symbols or concrete values, so differentiation goes
        through fresh dummies and substitutes back -- ``diff`` needs a symbol
        to differentiate with respect to, and an operating point is often
        numeric even under a symbolic toolkit.
        """
        if fallback is not None:
            return fallback

        import sympy
        dummies = [sympy.Dummy() for _ in range(len(x))]
        f = func(self.array(dummies), params, epar, self)
        back = dict(zip(dummies, x))
        J = self.zeros((len(f), len(x)))
        for k in range(len(f)):
            for j in range(len(x)):
                J[k, j] = sympy.diff(f[k], dummies[j]).subs(back)
        return J

    def nth_derivative(self, func, at, order):
        """Exact ``d^order func / d at^order`` -- no differencing, any order.

        The base class differences ``order`` times, which loses roughly
        ``order`` significant digits; sympy differentiates in closed form, so
        the third derivative is as exact as the first.  That is the whole
        reason distortion analysis can take its Taylor coefficients from the
        model rather than from a separate declaration.
        """
        import sympy
        t = sympy.Dummy()
        return sympy.diff(func(t), t, order).subs(t, at)


class NumDenMixin:
    """The ``N(s)/D(s)`` contract shared by every fraction-free toolkit.

    Implement :meth:`linearsolver_num_den` (``A x = b`` ->
    ``(numerator_vector, shared_denominator)``, with ``x_i = numer_i / den``)
    and this mixin provides the standard divided :meth:`linearsolver` and
    declares ``supports('num_den')``.

    Deliberately *not* a base class carrying ``ac_solution``/``noise_psd`` as
    well: those two answer "what does a solve return", which
    :class:`SymbolicPolyToolkit` and :class:`DDDToolkit` do very differently
    (a flat ``NumDenSolution`` vs. a shared diagram) despite both satisfying
    this same num/den contract. Giving them a common base implied a kinship
    the code didn't use -- see problem P3 in ``doc/architecture.md``, which
    this mixin resolves.
    """
    poly = True

    def supports(self, capability):
        return capability == 'num_den'

    def linearsolver(self, A, b):
        """Solve ``A x = b`` fraction-free; returns the divided solution vector.

        Thin wrapper over :meth:`linearsolver_num_den` (``x_i = numer_i / den``)
        so the toolkit honours the standard solution-vector contract used by the
        analyses.
        """
        num, den = self.linearsolver_num_den(A, b)
        return np.array([ni / den for ni in num], dtype=object)


class SymbolicPolyToolkit(NumDenMixin, SymbolicToolkit):
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


class DDDToolkit(NumDenMixin, SymbolicToolkit):
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

    Shares the :class:`NumDenMixin` contract with :class:`SymbolicPolyToolkit`
    (both implement ``linearsolver_num_den`` and get the standard
    ``linearsolver`` for free) but is not a subclass of it: ``ac_solution`` and
    ``noise_psd`` here are diagram-aware and share nothing with
    ``SymbolicPolyToolkit``'s flat ``NumDenSolution`` form, so inheriting from
    it would have claimed a kinship this toolkit doesn't use (formerly problem
    P3 in ``doc/architecture.md``).
    """

    #: Refuse to expand a diagram into sympy beyond this many product terms.
    ddd_max_terms = 5000

    def supports(self, capability):
        return capability in ('num_den', 'ddd')

    def ac_solution(self, A, b, s, irefnode):
        from .dddresult import DDDSolution
        return DDDSolution(A, b, s, irefnode, max_terms=self.ddd_max_terms)

    def _ddd_family(self, A):
        """Shared family for a matrix, reused across calls on the same matrix.

        ``TwoPortAnalysis`` asks for a determinant and then a series of
        cofactors of the *same* matrix.  Those are all minors of one matrix, so
        rebuilding per call would throw away exactly the sharing that makes a
        diagram worth having.  One entry is enough -- an analysis works through
        one matrix at a time.
        """
        from .ddd import DDDFamily
        Am = sympy.Matrix(A)
        key = (Am.rows, Am.cols, tuple(Am))
        cached = getattr(self, '_family_cache', None)
        if cached is not None and cached[0] == key:
            return cached[1]
        family = DDDFamily(Am, sympy.zeros(Am.rows, 1))
        self._family_cache = (key, family)
        return family

    def det(self, A):
        """Determinant, built as a diagram.

        Overriding this is what lets the existing analyses run on diagrams
        without being modified: ``TwoPortAnalysis`` and ``FeedbackLoopAnalysis``
        both reach their answer through ``toolkit.det``, so they need no
        knowledge of the representation.

        The result is a sympy expression, so it expands -- guarded by
        ``ddd_max_terms``, as every compatibility path here is.
        """
        return self._ddd_family(A).denominator.to_sympy(
            max_terms=self.ddd_max_terms)

    def cofactor(self, A, i, j):
        """Cofactor, from the same shared family as :meth:`det`.

        The two-port parameters are a determinant followed by cofactors of the
        one matrix, and taking them from a single family is the difference
        between one construction and one per parameter.
        """
        from .ddd import DDDSizeError
        family = self._ddd_family(A)
        minor = family.cofactor(i, j)
        sign = -1 if (i + j) % 2 else 1
        return sign * minor.to_sympy(max_terms=self.ddd_max_terms)

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

        ## Guard first, then evaluate.  Flattening each transimpedance with its
        ## own call would re-walk the family once per unknown -- the same
        ## quadratic re-walk that cost 200x in the reduced solve.
        for i in range(Ym.rows):
            flatten(nums[i], 'transimpedance %d' % i)
        D = flatten(family.denominator, 'the determinant')
        N = sympy.Matrix([nums[i].eval() for i in range(Ym.rows)])
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
