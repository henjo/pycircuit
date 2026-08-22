# -*- coding: utf-8 -*-
# Copyright (c) 2008-2026 Pycircuit Development Team
# See LICENSE for details.
"""Behavioural: a Verilog-A-shaped element DSL compiled to MNA stamps.

An element's behaviour is written as *contribution statements* --
``Contribution(b.I, expr)`` and ``Contribution(b.V, expr)``, the DSL forms
of Verilog-A's ``I(a,b) <+ expr`` and ``V(a,b) <+ expr`` -- inside an
``analog()`` method whose arguments name the terminals.  The metaclass
compiles those equations, symbolically, into everything a pycircuit
element supplies: ``i``/``q``/``u`` vectors, exact ``G``/``C`` Jacobians,
``dudt``, the DC variants for pinned states, ``eval_i_pure``/
``eval_q_pure`` for the JAX batched path, ``state_ic`` seeding and
``periodic_states`` declarations.  See ``hdl.md`` (repository root) for
theory, the capability map against Verilog-A, and usage.

The compilation model (the same set of ideas every equation-defined
device compiler uses -- Qucs EDD, Xyce ABM, ngspice B-sources -- done
here with sympy so the Jacobians are exact and closed-form):

* the element's unknown vector is ``[terminal nodes, internal nodes,
  state nodes, branch currents]``, matching ``Circuit.update_node_map``'s
  layout exactly;
* an ``I``-contribution adds ``+expr``/``-expr`` to the two KCL rows;
* a ``V``-contribution creates a branch-current unknown ``i_b`` (added to
  the KCL rows) and the branch row ``-(v_plus - v_minus) + expr = 0``;
* each additive term routes by kind: ``ddt(arg)`` (optionally scaled by
  an x-independent factor) into the charge vector ``q``, terms free of
  circuit quantities into the source vector ``u(t)`` (they may depend on
  the time symbol ``TIME``), everything else into ``i``;
* each ``idt``/``idtmod`` application becomes a *state*: an internal node
  ``s`` with the extra equation ``ds/dt = arg`` (``q``-row ``s``,
  ``i``-row ``-arg``), the application replaced by ``s`` (``idt``) or by
  the floored wrap of ``s`` (``idtmod``, idtmod.md sec. 1.1);
* ``G = di/dx`` and ``C = dq/dx`` are sympy Jacobians of those vectors,
  compiled with ``lambdify(..., cse=True)``.

Residual convention (transient.py): ``i(x) + dq/dt + u(t) = 0``.
"""

from pycircuit.circuit import circuit
from pycircuit.circuit.circuit import defaultepar
import pycircuit.utilities.param as param

import sympy
import numpy as np

import inspect
from copy import copy


class Node(circuit.Node):
    @property
    def V(self):
        return Quantity('V', self)


class Branch(circuit.Branch):
    @property
    def V(self):
        return Quantity('V', self)

    @property
    def I(self):
        return Quantity('I', self)


## NOT a sympy.Symbol subclass, deliberately.  It used to be
## `class Parameter(param.Parameter, sympy.Symbol)`, but sympy interns
## Symbols by name: two Parameters named 'c0' in two different elements
## were THE SAME OBJECT, and whichever was declared last set `default` for
## both -- a silent cross-element parameter leak.  The mixin bought
## nothing, since `generate_code` builds its own `Symbol(name)` for each
## declared parameter and binds values from `iparv` at call time.
Parameter = param.Parameter


class ddt(sympy.Function):
    """Time derivative: routes its argument into the charge vector."""
    nargs = (1,)


class idt(sympy.Function):
    """Time integral -- Verilog-A ``idt(expr[, ic])``.

    Compiles to an extra state unknown ``s`` with equation ``ds/dt =
    expr``; the application evaluates to ``s``.  With ``ic`` given, the
    DC solve pins ``s = ic`` (the LRM: "the initial condition shall force
    the DC solution") and ``uic=True`` seeds it; without, the element has
    a DC solution only if feedback zeroes ``expr`` (idtmod.md sec. 5.1).
    """
    nargs = (1, 2)


class idtmod(sympy.Function):
    """Circular integral -- Verilog-A ``idtmod(expr, ic, modulus, offset)``.

    The state integrates unbounded within a step and evaluates through the
    floored wrap into ``[offset, offset+modulus)``; the compiled element
    declares the state through ``periodic_states()`` so the Phase-2 gauge
    shift keeps it bounded between steps (idtmod.md sec. 5.2), and seeds/
    pins from the wrapped ``ic``.
    """
    nargs = (1, 2, 3, 4)


#: The time symbol usable in source terms: ``sin(2*pi*f*TIME)`` etc.
TIME = sympy.Symbol('_hdl_time', real=True)

#: Circuit temperature in kelvin -- Verilog-A's ``$temperature``; bound from
#: ``epar.T`` at evaluation time.
TEMP = sympy.Symbol('_hdl_temp', positive=True)

#: Noise frequency in Hz inside ``flicker_noise`` PSDs; bound from the noise
#: analysis's ``w`` at evaluation time.
FREQ = sympy.Symbol('_hdl_freq', positive=True)

#: Boltzmann constant and elementary charge, for ``vt()`` and noise PSDs.
#: Taken from the numeric toolkit so a generated noise stamp equals the
#: hand-written one exactly (``elements.R.CY``) rather than nearly.
from pycircuit.circuit.constants import kboltzmann as KBOLTZMANN
from pycircuit.circuit.constants import qelectron as QELECTRON


def vt(T=TEMP):
    """Thermal voltage ``k*T/q`` -- Verilog-A's ``$vt``."""
    return KBOLTZMANN * T / QELECTRON


class _wrapfloor(sympy.Function):
    """``floor`` whose symbolic derivative is 0 -- exact almost everywhere.

    The idtmod wrap ``s - m*floor((s-o)/m)`` has slope 1 a.e. and a jump
    of exactly ``m`` at the boundary that no Jacobian can state (idtmod.md
    sec. 5.1); an unevaluated ``Derivative(floor)`` would kill lambdify
    instead.  Printed as the toolkit's ``floor``.
    """
    nargs = (1,)

    def fdiff(self, argindex=1):
        return sympy.Integer(0)


class white_noise(sympy.Function):
    """``white_noise(pwr)`` -- a flat current PSD on the contributed branch.

    Contributes ONLY to the noise correlation matrix ``CY`` (LRM: noise
    functions return zero in DC/transient); ``pwr`` may use ``TEMP``.
    """
    nargs = (1, 2)


class flicker_noise(sympy.Function):
    """``flicker_noise(pwr, exp)`` -- a ``pwr/f^exp`` current PSD."""
    nargs = (2, 3)


def limexp(x, x0=80.0):
    """Verilog-A ``limexp``: ``exp`` continued linearly above ``x0``.

    C1-continuous at the seam, so Newton keeps a bounded derivative
    instead of overflowing -- the standard convergence aid.
    """
    x = sympy.sympify(x)
    return sympy.Piecewise((sympy.exp(x), x < x0),
                           (sympy.exp(x0) * (1 + x - x0), True))


def discontinuity(degree=0):
    """Verilog-A's ``$discontinuity(degree)``: tell the solver the model
    has a discontinuity in its ``degree``-th derivative here.

    **Parsed and ignored, deliberately.**  It is advisory: it lets a
    simulator drop integration order rather than discover the corner by
    rejection, and pycircuit does that from breakpoints instead (see
    ``Cross``).  gnucap accepts it and generates an empty body for the
    same reason.  Accepting it costs nothing and lets a model written for
    another simulator compile unchanged -- 41 call sites in the vacask
    device library -- whereas rejecting it would fail those models for a
    hint they can do without.  If it ever becomes load-bearing, the
    honest implementation is to schedule a breakpoint, which is what
    ``Cross`` already does.
    """
    return sympy.Integer(0)


class _Laplace(sympy.Function):
    """Marker for a laplace_* application; see :func:`laplace_nd`."""
    nargs = (3,)


def laplace_nd(expr, num, den):
    """Verilog-A's ``laplace_nd(expr, num[], den[])``: apply the rational
    transfer function ``H(s) = sum(num[k] s^k) / sum(den[k] s^k)``.

    Coefficients ascend in powers of ``s``, as the LRM specifies.  The
    filter is realised as **state equations**, not as a convolution: an
    order-``N`` denominator introduces ``N`` unknowns in controllable
    canonical form, so the simulator integrates them with its own method
    and order and the Jacobian is exact through them.  That is the same
    construction ``idt`` uses, generalised.

    ``H`` must be proper (``len(num) <= len(den)``): an improper transfer
    function differentiates its input, which is ``ddt``'s job and needs a
    charge, not a state.
    """
    num = [sympy.sympify(c) for c in num]
    den = [sympy.sympify(c) for c in den]
    if len(den) < 2:
        raise ValueError('laplace_nd needs a denominator of order >= 1; '
                         'a constant gain is just multiplication')
    if len(num) > len(den):
        raise NotImplementedError(
            'laplace_nd requires a proper transfer function '
            '(len(num) <= len(den)); an improper one differentiates its '
            'input, which is what ddt is for')
    return _Laplace(sympy.sympify(expr), sympy.Tuple(*num),
                    sympy.Tuple(*den))


def laplace_zp(expr, zeros, poles):
    """``laplace_zp(expr, zeros[], poles[])`` -- zeros and poles instead of
    coefficients.  Each is a flat list of (real, imaginary) pairs, as the
    LRM specifies, and the pair ``(r, i)`` contributes a factor
    ``1 - s/(r + j i)`` (or ``1 - s/r`` when ``i`` is zero).

    Converted to :func:`laplace_nd` by expanding the products, so the
    realisation and everything downstream is identical.
    """
    sv = sympy.Symbol('_lap_s')

    def poly_of(pairs):
        expr_ = sympy.Integer(1)
        it = list(pairs)
        for k in range(0, len(it), 2):
            re_, im_ = sympy.sympify(it[k]), sympy.sympify(it[k + 1])
            ## `.is_zero`, NOT `== 0`: sympy's `==` is STRUCTURAL, so
            ## `Float(0.0) == 0` is False and a real pole written with a
            ## float imaginary part took the conjugate-pair branch --
            ## silently producing a filter of twice the intended order.
            if im_.is_zero:
                if re_.is_zero:
                    expr_ *= sv
                else:
                    expr_ *= (1 - sv / re_)
            else:
                ## A conjugate pair contributes a real quadratic.
                mag2 = re_ ** 2 + im_ ** 2
                expr_ *= (1 - 2 * re_ * sv / mag2 + sv ** 2 / mag2)
        return sympy.Poly(sympy.expand(expr_), sv).all_coeffs()[::-1]

    return laplace_nd(expr, poly_of(zeros), poly_of(poles))


class _Limit(sympy.Function):
    """Marker for a limited probe; see :func:`limit_pnj`."""
    nargs = (3,)


def limit_pnj(probe, IS, VT):
    """Verilog-A's ``$limit(probe, "pnjlim", ...)``: ask the simulator to
    bound this branch voltage's per-iteration excursion.

    ``probe`` must be a branch potential (``b.V``); ``IS`` and ``VT`` are
    the junction's saturation current and thermal voltage, which is what
    SPICE's ``pnjlim`` needs to place its critical voltage.  The
    expression evaluates to the probe itself, so it is written inline::

        Contribution(b.I, IS * (limexp(limit_pnj(b.V, IS, vt()) / vt()) - 1))

    The generated element implements the STATE-FREE limiter convention
    (``SubCircuit.limit``): it returns a limited copy of the solution
    sub-vector rather than keeping a private linearisation point.  That
    is what the codebase prefers and what a traced backend could ever
    support -- and it means the limiter moves the *solution*, so the
    residual and Jacobian are never taken at different points, which is
    the inconsistency classical limiting is criticised for.

    Prefer PCNR (``pcnr=True``) where the element qualifies: it makes the
    limited quantity an explicit unknown owned by this device, so two
    devices sharing a node cannot fight over it.  ``limit_pnj`` is the
    fallback for models PCNR cannot take, and for ordinary runs, where
    ``pcnr=False`` is the default.
    """
    if not (isinstance(probe, Quantity) and probe.isbranch
            and probe.quantity == 'V'):
        raise ValueError('limit_pnj limits a branch potential (b.V); '
                         'got %r' % (probe,))
    return _Limit(probe, sympy.sympify(IS), sympy.sympify(VT))


def param_given(name):
    """Verilog-A's ``$param_given(name)``: 1.0 when the instance was given
    that parameter explicitly, 0.0 when it fell back to its default.

    A *runtime* value, not a compile-time one, because the element is
    compiled once per class while givenness is a property of each
    instance.  Use it in a ``Piecewise`` exactly as Verilog-A uses it in
    an ``if``::

        rs_eff = sympy.Piecewise((rs, param_given('rs') > 0.5), (0.0, True))

    By call count this is the most-used system function in real compact
    models, which is why it is here before more glamorous operators.
    """
    return sympy.Symbol('_hdl_given_%s' % name, real=True)


def ac_stim(mag=1.0, phase=0.0):
    """Verilog-A's ``ac_stim``: a small-signal excitation, live in AC
    analysis and identically zero everywhere else.

    Contributed like any other term -- ``Contribution(b.I, ac_stim(1.0))``
    -- and routed to an AC-only source vector, so a behavioural element
    can BE a small-signal stimulus rather than merely transparent to one.
    """
    return _AcStim(sympy.sympify(mag), sympy.sympify(phase))


class _AcStim(sympy.Function):
    nargs = (2,)


def ddx(expr, probe):
    """Verilog-A ``ddx(expr, probe)``: the symbolic partial derivative of
    ``expr`` with respect to a ``V``/``I`` probe, other probes held fixed.

    Pure sympy -- ``Quantity`` is an atom, so this is a substitution and a
    ``diff``.  The probe must appear in ``expr`` as the SAME quantity
    object (``b.V`` matched against ``b.V``, not against an expanded
    node-voltage difference).
    """
    d = sympy.Dummy()
    return sympy.sympify(expr).subs(probe, d).diff(d).subs(d, probe)


class Quantity(sympy.AtomicExpr):
    """Reference to a voltage or current of a branch or node.

    A sympy atom, so it takes part in symbolic expressions (arithmetic,
    ``.subs()``, ``Matrix.jacobian()``, ...).
    """
    ## A voltage or a current is REAL.  Without this assumption sympy
    ## refuses to build `Quantity < 80` and every Piecewise over an
    ## operating point -- limexp included -- dies with "Invalid comparison
    ## of non-real".
    is_real = True
    is_commutative = True

    def __new__(cls, quantity, branch_or_node):
        if quantity not in ('V', 'I'):
            raise ValueError("quantity must be either 'V' or 'I'")
        if not isinstance(branch_or_node, (Node, Branch)):
            raise ValueError('branch_or_node must be a Branch or Node object')
        if quantity == 'I' and isinstance(branch_or_node, Node):
            raise ValueError('Current can only be taken on branches')

        obj = super().__new__(cls)
        obj.quantity = quantity
        obj.branch_or_node = branch_or_node
        return obj

    def _hashable_content(self):
        if isinstance(self.branch_or_node, Branch):
            key = ('branch', self.branch_or_node.plus.name,
                   self.branch_or_node.minus.name)
        else:
            key = ('node', self.branch_or_node.name)
        return (self.quantity,) + key

    @property
    def name(self):
        ## sympy's StrPrinter dispatches on the CLASS NAME and hits
        ## _print_Quantity (written for sympy.physics.units.Quantity), which
        ## reads .name -- without this property, printing/sorting any Mul
        ## holding a Quantity next to a Function dies at class creation.
        return repr(self)

    @property
    def isnode(self):
        return isinstance(self.branch_or_node, Node)

    @property
    def isbranch(self):
        return isinstance(self.branch_or_node, Branch)

    def __repr__(self):
        if self.isbranch:
            return '%s(%s,%s)' % (self.quantity, self.branch_or_node.plus.name,
                                  self.branch_or_node.minus.name)
        return '%s(%s)' % (self.quantity, self.branch_or_node.name)

    __str__ = __repr__


class Statement:
    pass


class Cross(Statement):
    """Verilog-A's ``@(cross(expr, direction))``: ask the transient solver
    to place a timepoint where ``expr`` crosses zero.

    ``direction`` is ``+1`` (rising), ``-1`` (falling) or ``0`` (either),
    as in the LRM.  Returned alongside contributions::

        return (Contribution(b_out.V, ...),
                Cross(b_in.V - vref, +1))

    The generated element predicts the crossing by linear extrapolation
    from the last two ACCEPTED points and publishes it through
    ``next_event``, which is the contract ``SubCircuit.next_event`` polls
    and the same first-order prediction gnucap uses.  A prediction only
    has to BRACKET the corner: the solver's break-step machinery lands on
    it and drops the integration order across it.

    Without this the controller finds a comparator's edge by rejecting
    steps until the local error is small enough -- correct, but it smears
    the corner and costs a rejection storm at every transition.
    """

    def __init__(self, expr, direction=0):
        if direction not in (-1, 0, 1):
            raise ValueError('cross direction must be -1, 0 or +1')
        self.expr = sympy.sympify(expr)
        self.direction = direction

    def nodes(self):
        nodes = set()
        for atom in self.expr.atoms():
            if isinstance(atom, Quantity):
                if atom.isbranch:
                    nodes.add(atom.branch_or_node.plus)
                    nodes.add(atom.branch_or_node.minus)
                else:
                    nodes.add(atom.branch_or_node)
        return nodes


class Contribution(Statement):
    """``Contribution(b.I, expr)`` or ``Contribution(b.V, expr)`` --
    the DSL form of Verilog-A's ``<+``."""

    def __init__(self, lhs, rhs):
        if not isinstance(lhs, Quantity) or not lhs.isbranch:
            raise ValueError(
                'the contribution target must be the I or V of a Branch '
                '(got %r); node contributions are not part of the language'
                % (lhs,))
        self.lhs = lhs
        self.rhs = sympy.sympify(rhs)

    def nodes(self):
        """The set of node objects referred to in lhs and rhs."""
        nodes = set()
        for atom in self.lhs.atoms() | self.rhs.atoms():
            if isinstance(atom, Quantity):
                if atom.isbranch:
                    nodes.add(atom.branch_or_node.plus)
                    nodes.add(atom.branch_or_node.minus)
                else:
                    nodes.add(atom.branch_or_node)
        return nodes


def _quantity_free(expr):
    """True when nothing in ``expr`` refers to a circuit quantity or state."""
    return not any(isinstance(a, Quantity) for a in expr.atoms()) and \
        not expr.atoms(_StateSymbol)


class _StateSymbol(sympy.Symbol):
    """Symbol standing for an idt/idtmod state unknown."""

    def __new__(cls, name):
        ## real=True for the same reason as Quantity: a Piecewise over a
        ## state must be constructible.
        return super().__new__(cls, name, real=True)


def _split_terms(rhs, xset):
    """Route an expression's additive terms into (i, q, u) parts.

    Called AFTER ``resolve()``, so circuit quantities and states appear as
    the x-symbols in ``xset``; a term's kind is decided by which of those
    it touches.  ``ddt(arg)`` terms -- bare, or scaled by an x-independent
    factor -- become charge; an x-DEPENDENT scaling of ``ddt`` is refused
    loudly (``g(v)*ddt(v)`` is not the derivative of any charge; move the
    factor inside the ``ddt``).  x-free terms become sources (they may
    contain ``TIME``); the rest is static current.
    """
    def xfree(e):
        return not (e.free_symbols & xset)

    rhs = sympy.expand(rhs)
    terms = rhs.args if rhs.is_Add else (rhs,)

    iterms, qterms, uterms, nterms, acterms = [], [], [], [], []
    for term in terms:
        acs = term.atoms(_AcStim)
        if acs:
            if isinstance(term, _AcStim):
                app, coeff = term, sympy.Integer(1)
            elif term.is_Mul:
                apps = [f for f in term.args if isinstance(f, _AcStim)]
                rest = [f for f in term.args if not isinstance(f, _AcStim)]
                coeff = sympy.Mul(*rest)
                if len(apps) != 1 or not xfree(coeff):
                    raise NotImplementedError(
                        'ac_stim may appear as a term or scaled by an '
                        'x-independent factor; %r is neither' % (term,))
                app = apps[0]
            else:
                raise NotImplementedError(
                    'ac_stim may only appear additively: %r' % (term,))
            mag, phase = app.args
            acterms.append(coeff * mag * sympy.exp(sympy.I * phase
                                                   * sympy.pi / 180))
            continue
        noises = term.atoms(white_noise) | term.atoms(flicker_noise)
        if noises:
            if isinstance(term, (white_noise, flicker_noise)):
                app, coeff = term, sympy.Integer(1)
            elif term.is_Mul:
                napps = [f for f in term.args
                         if isinstance(f, (white_noise, flicker_noise))]
                rest = [f for f in term.args
                        if not isinstance(f, (white_noise, flicker_noise))]
                coeff = sympy.Mul(*rest)
                if len(napps) != 1 or not xfree(coeff):
                    raise NotImplementedError(
                        'noise functions may appear as a term or scaled by '
                        'an x-independent factor; %r is neither' % (term,))
                app = napps[0]
            else:
                raise NotImplementedError(
                    'noise functions may only appear additively: %r' % (term,))
            if isinstance(app, white_noise):
                psd = coeff * app.args[0]
            else:
                psd = coeff * app.args[0] / FREQ ** app.args[1]
            nterms.append(psd)
            continue
        ddts = term.atoms(ddt)
        if ddts:
            if isinstance(term, ddt):
                qterms.append(sympy.expand(term.args[0]))
                continue
            if term.is_Mul:
                dd = [f for f in term.args if isinstance(f, ddt)]
                rest = [f for f in term.args if not isinstance(f, ddt)]
                coeff = sympy.Mul(*rest)
                if len(dd) == 1 and xfree(coeff) and not coeff.atoms(ddt):
                    qterms.append(sympy.expand(coeff * dd[0].args[0]))
                    continue
            raise NotImplementedError(
                'ddt may appear as a term or scaled by an x-independent '
                'factor; %r is neither. A state-dependent factor is not '
                'the derivative of any charge -- move it inside the ddt.'
                % (term,))
        if xfree(term):
            ## `u` carries only genuinely TIME-VARYING excitation; a
            ## constant belongs to the device's static characteristic
            ## (a diode's -IS), which keeps generated stamps identical to
            ## the hand-written ones and keeps `u` meaning one thing.
            (uterms if term.has(TIME) else iterms).append(term)
        else:
            iterms.append(term)
    return (sympy.Add(*iterms), sympy.Add(*qterms), sympy.Add(*uterms),
            nterms, sympy.Add(*acterms))


def generate_code(cls):
    """Compile the class's ``analog()`` into the element interface.

    Returns a dict consumed by :class:`BehaviouralMeta`; see the module
    docstring for the compilation model.
    """
    terminalnames = inspect.getfullargspec(cls.analog)[0]
    terminalnodes = [Node(name) for name in terminalnames]

    ## Inject parameter names as sympy symbols into the analog() globals so
    ## the body refers to instance parameters by bare name; the compiled
    ## functions take them as trailing arguments, supplied from the RESOLVED
    ## self.iparv at call time.
    analogfunc = copy(cls.analog)
    paramnames = [p.name for p in cls.instparams]
    paramsyms = [sympy.Symbol(name) for name in paramnames]
    analogfunc.__globals__.update(dict(zip(paramnames, paramsyms)))

    statements = analogfunc(*terminalnodes)
    if isinstance(statements, Statement):
        statements = (statements,)
    crossings = [st for st in statements if isinstance(st, Cross)]
    statements = [st for st in statements if not isinstance(st, Cross)]

    ## ------------------------------------------------------------------
    ## Pass 1 -- inventory: nodes, branches with current unknowns (every
    ## V-contributed or I-probed branch), and idt/idtmod states.
    nodes = set()
    vbranches = []      # branches that get a current unknown, in order
    ibranch_keys = set()

    def branch_key(branch):
        return (branch.plus.name, branch.minus.name)

    for st in crossings:
        nodes.update(st.nodes())
    for st in statements:
        nodes.update(st.nodes())
        if st.lhs.quantity == 'V':
            key = branch_key(st.lhs.branch_or_node)
            if key not in ibranch_keys:
                ibranch_keys.add(key)
                vbranches.append(st.lhs.branch_or_node)
        for atom in st.rhs.atoms(Quantity):
            if atom.isbranch and atom.quantity == 'I':
                key = branch_key(atom.branch_or_node)
                if key not in ibranch_keys:
                    ibranch_keys.add(key)
                    vbranches.append(atom.branch_or_node)

    ## ------------------------------------------------------------------
    ## NODE COLLAPSE.  `V(a,b) <+ 0` is Verilog-A's way of shorting a
    ## branch, and the standard idiom for an optional series resistance:
    ##
    ##     if (rs) I(a, ia) <+ V(a, ia)/rs;  else V(a, ia) <+ 0;
    ##
    ## Rather than carry a zero-volt source (a real branch unknown and a
    ## row for it), the two nodes become one.  Restricted, as gnucap
    ## restricts it, to UNCONDITIONAL contributions -- a collapse that
    ## depends on the operating point would change the sparsity pattern
    ## per iteration, which is a different and much larger feature.
    ##
    ## Only an INTERNAL node can be absorbed: terminals belong to the
    ## parent circuit's node map and this element cannot merge them.  A
    ## `V <+ 0` between two terminals therefore stays a zero-volt source,
    ## which is exactly what it means and costs one unknown.
    collapse = {}
    kept_statements = []
    for st in statements:
        b = st.lhs.branch_or_node
        if st.lhs.quantity == 'V' and sympy.sympify(st.rhs) == 0:
            p_, m_ = b.plus, b.minus
            p_int = p_ not in terminalnodes
            m_int = m_ not in terminalnodes
            if m_int:
                collapse[m_.name] = p_.name
                continue
            if p_int:
                collapse[p_.name] = m_.name
                continue
        kept_statements.append(st)
    if collapse:
        ## Chase chains (a -> b -> terminal) to a fixed point.
        for _ in range(len(collapse) + 1):
            collapse = {k: collapse.get(v, v) for k, v in collapse.items()}
        statements = kept_statements
        alias = {}
        for nd in list(nodes):
            tgt = collapse.get(nd.name)
            if tgt is not None:
                alias[nd] = next(o for o in nodes if o.name == tgt)
        nodes = {alias.get(nd, nd) for nd in nodes}
        ## Rewrite the surviving statements' node references.
        for st in statements:
            subs_q = {}
            for atom in st.rhs.atoms(Quantity):
                bn = atom.branch_or_node
                if isinstance(bn, Branch):
                    np_, nm_ = alias.get(bn.plus, bn.plus), \
                        alias.get(bn.minus, bn.minus)
                    if (np_, nm_) != (bn.plus, bn.minus):
                        nb = Branch(np_, nm_)
                        subs_q[atom] = Quantity(atom.quantity, nb)
                elif bn in alias:
                    subs_q[atom] = Quantity(atom.quantity, alias[bn])
            if subs_q:
                st.rhs = st.rhs.subs(subs_q)
            lb = st.lhs.branch_or_node
            np_, nm_ = alias.get(lb.plus, lb.plus), alias.get(lb.minus,
                                                              lb.minus)
            if (np_, nm_) != (lb.plus, lb.minus):
                st.lhs = Quantity(st.lhs.quantity, Branch(np_, nm_))
        ## Re-derive the branch inventory: a collapsed branch is gone.
        vbranches = [br for br in vbranches
                     if br.plus.name not in collapse
                     and br.minus.name not in collapse]
        ibranch_keys = {branch_key(br) for br in vbranches}

    ## SWITCH BRANCHES are refused, per hdl.md sec. 9 phase D.  In
    ## Verilog-A a branch may be a potential source in one operating
    ## region and a flow source in another, the two contributions being
    ## mutually exclusive under a condition.  Implementing that means
    ## re-deciding the stamp -- and the sparsity pattern -- every Newton
    ## iteration, which is gnucap's whole `_pot` re-stamping machine.
    ##
    ## Accepting BOTH unconditionally, as this compiler silently did,
    ## produces a defined but different element: a voltage source with a
    ## conductance in parallel, which is not what the model says.  A
    ## plausible wrong answer is the expensive failure, so it is named.
    _v_branches = {branch_key(st.lhs.branch_or_node)
                   for st in statements if st.lhs.quantity == 'V'}
    _i_branches = {branch_key(st.lhs.branch_or_node)
                   for st in statements if st.lhs.quantity == 'I'}
    _both = _v_branches & _i_branches
    if _both:
        raise NotImplementedError(
            'branch %s has both V and I contributions (a Verilog-A "switch '
            'branch"). Selecting between them per operating point changes '
            'the matrix structure each iteration and is not implemented; '
            'accepting both unconditionally would give a voltage source '
            'with a conductance in parallel, which is not what the model '
            'says. Split the regions into separate elements, or use a '
            'Piecewise inside ONE contribution.'
            % ', '.join('(%s,%s)' % k for k in sorted(_both)))

    internalnodes = sorted(nodes - set(terminalnodes), key=lambda n: n.name)

    ## States: each distinct idt/idtmod APPLICATION is one state.  Walk in
    ## a deterministic order.
    states = []          # (state_symbol, kind, args); kind idt/idtmod/lap
    state_subst = {}

    ## LAPLACE -> STATE EQUATIONS.  For H(s) = N(s)/D(s) with
    ## D = d0 + d1 s + ... + dN s^N, controllable canonical form gives
    ##
    ##     z1' = z2,  z2' = z3,  ...,
    ##     zN' = (u - d0 z1 - d1 z2 - ... - d_{N-1} zN) / dN
    ##     y   = n0 z1 + n1 z2 + ...
    ##
    ## i.e. N first-order states, which is exactly what this compiler
    ## already emits for `idt` -- one per integrator in the chain.  The
    ## simulator then integrates them with its own method and order, and
    ## the Jacobian is exact through them, which a convolution-based
    ## implementation cannot offer.
    for st in statements:
        for app in sorted(st.rhs.atoms(_Laplace), key=sympy.default_sort_key):
            if app in state_subst:
                continue
            u_expr, num, den = app.args
            num = list(num)
            den = list(den)
            order = len(den) - 1
            zs = []
            for _k in range(order):
                sym = _StateSymbol('_state%d' % len(states))
                states.append((sym, 'lap', None))    # rhs filled in below
                zs.append(sym)
            ## z_k' = z_{k+1} for k < N; the last row carries the input.
            rhs = []
            for k in range(order - 1):
                rhs.append(zs[k + 1])
            last = u_expr
            for k in range(order):
                last = last - den[k] * zs[k]
            rhs.append(last / den[order])
            for k in range(order):
                idx = len(states) - order + k
                sym, _kind, _ = states[idx]
                states[idx] = (sym, 'lap', (rhs[k],))
            out = sympy.Integer(0)
            for k, nk in enumerate(num):
                if k < order:
                    out += nk * zs[k]
                elif nk != 0:
                    ## num of the same order as den: the direct term
                    ## n_N/d_N, plus its correction through the chain.
                    out += nk * rhs[order - 1]
            state_subst[app] = out
    for st in statements:
        for func_cls, kind in ((idt, 'idt'), (idtmod, 'idtmod')):
            for app in sorted(st.rhs.atoms(func_cls), key=sympy.default_sort_key):
                if app in state_subst:
                    continue
                sym = _StateSymbol('_state%d' % len(states))
                states.append((sym, kind, app.args))
                if kind == 'idt':
                    state_subst[app] = sym
                else:
                    args = app.args
                    m = args[2] if len(args) > 2 else sympy.Integer(1)
                    o = args[3] if len(args) > 3 else sympy.Integer(0)
                    ## the floored wrap, idtmod.md sec. 1.1
                    state_subst[app] = \
                        sym - m * _wrapfloor((sym - o) / m)

    for _sym, _kind, args in states:
        if args is None:
            continue
        if any(a.atoms(idt) or a.atoms(idtmod) for a in args):
            raise NotImplementedError(
                'nested idt/idtmod applications are not supported yet')

    ## ------------------------------------------------------------------
    ## The unknown vector, in Circuit.update_node_map order:
    ## terminals, internal nodes, state nodes, branch currents.
    statenodenames = ['_state%d' % k for k in range(len(states))]
    xlabels = ([n.name for n in terminalnodes]
               + [n.name for n in internalnodes]
               + statenodenames
               + ['_i_br%d' % k for k in range(len(vbranches))])
    xsyms = [sympy.Symbol('_x%d' % k, real=True)
             for k in range(len(xlabels))]

    index_of = {}
    for k, node in enumerate(terminalnodes + internalnodes):
        index_of[('node', node.name)] = k
    off = len(terminalnodes) + len(internalnodes)
    for k, (sym, _kind, _args) in enumerate(states):
        index_of[('state', sym.name)] = off + k
    off += len(states)
    for k, br in enumerate(vbranches):
        index_of[('branch', branch_key(br))] = off + k

    n = len(xlabels)

    ## Substitutions: node voltages, branch-current probes, state symbols.
    subst = {}
    for node in terminalnodes + internalnodes:
        subst[Quantity('V', node)] = xsyms[index_of[('node', node.name)]]
    for br in vbranches:
        subst[Quantity('I', br)] = xsyms[index_of[('branch', branch_key(br))]]
    for sym, _kind, _args in states:
        subst[sym] = xsyms[index_of[('state', sym.name)]]

    def resolve(expr):
        ## STATES FIRST: the recorded idt/idtmod applications are in their
        ## original quantity form, and substituting branch voltages before
        ## them would rewrite the applications out from under the mapping
        ## (idt(V(a,b)) becoming idt(V(a)-V(b)), matching nothing).
        expr = expr.subs(state_subst)
        bsubs = {}
        for atom in expr.atoms(Quantity):
            if atom.isbranch and atom.quantity == 'V':
                b = atom.branch_or_node
                bsubs[atom] = Quantity('V', b.plus) - Quantity('V', b.minus)
            elif atom.isbranch and atom.quantity == 'I' and \
                    branch_key(atom.branch_or_node) not in ibranch_keys:
                raise ValueError(
                    'I(%s,%s) is probed but that branch has no current '
                    'unknown; probe only V-contributed branches.'
                    % branch_key(atom.branch_or_node))
        expr = expr.subs(bsubs)
        return expr.subs(subst)

    ## $limit: record (branch, IS, VT) for each limited probe, then
    ## replace the marker by the probe itself -- it is a request to the
    ## SOLVER, not a change to the equations.
    limits = []
    for st in statements:
        subs_l = {}
        for app in sorted(st.rhs.atoms(_Limit), key=sympy.default_sort_key):
            probe, IS_l, VT_l = app.args
            b_l = probe.branch_or_node
            key = (index_of[('node', b_l.plus.name)],
                   index_of[('node', b_l.minus.name)])
            if not any(k[0] == key for k in limits):
                limits.append((key, IS_l, VT_l))
            subs_l[app] = probe
        if subs_l:
            st.rhs = st.rhs.subs(subs_l)

    ## ------------------------------------------------------------------
    ## Pass 2 -- assemble the residual vectors.
    ivec = [sympy.Integer(0)] * n
    qvec = [sympy.Integer(0)] * n
    uvec = [sympy.Integer(0)] * n
    acvec = [sympy.Integer(0)] * n
    CYacc = {}

    for st in statements:
        b = st.lhs.branch_or_node
        p = index_of[('node', b.plus.name)]
        m_ = index_of[('node', b.minus.name)]
        if st.lhs.quantity == 'I':
            ipart, qpart, upart, npart, acpart = _split_terms(
                resolve(st.rhs), set(xsyms))
            ivec[p] += ipart; ivec[m_] -= ipart
            qvec[p] += qpart; qvec[m_] -= qpart
            uvec[p] += upart; uvec[m_] -= upart
            acvec[p] += acpart; acvec[m_] -= acpart
            for psd in npart:
                ## An I-contributed noise PSD stamps like a noisy
                ## conductance: the R.CY pattern (elements.py).
                CYacc[(p, p)] = CYacc.get((p, p), 0) + psd
                CYacc[(m_, m_)] = CYacc.get((m_, m_), 0) + psd
                CYacc[(p, m_)] = CYacc.get((p, m_), 0) - psd
                CYacc[(m_, p)] = CYacc.get((m_, p), 0) - psd
        else:
            ## V-contribution: KCL rows carry the branch current; the
            ## branch row is -(v_p - v_m) + rhs = 0 (the elements.py
            ## convention -- compare Idtmod's output row).
            bi = index_of[('branch', branch_key(b))]
            ivec[p] += xsyms[bi]
            ivec[m_] -= xsyms[bi]
            ipart, qpart, upart, npart, acpart = _split_terms(
                resolve(st.rhs), set(xsyms))
            if acpart != 0:
                raise NotImplementedError(
                    'ac_stim on a V-contributed branch is not supported '
                    'yet; contribute it with I(...)')
            if npart:
                raise NotImplementedError(
                    'noise contributions on a V-contributed branch are not '
                    'supported yet; contribute noise with I(...)')
            vp, vm = xsyms[p], xsyms[m_]
            ivec[bi] += -(vp - vm) + ipart
            qvec[bi] += qpart
            uvec[bi] += upart

    ## State equations: ds/dt = arg  ->  q-row s, i-row -arg.
    dc_pins = {}     # x-index -> (pin expression for i, seed expression)
    periodic = []    # (x-index, modulus expr, offset expr)  [idtmod only]
    for sym, kind, args in states:
        k = index_of[('state', sym.name)]
        arg_i, arg_q, arg_u, arg_n, arg_ac = _split_terms(
            resolve(args[0]), set(xsyms))
        if arg_ac != 0:
            raise NotImplementedError('ac_stim inside idt/idtmod is not '
                                      'supported')
        if arg_n:
            raise NotImplementedError('noise inside idt/idtmod is not '
                                      'supported')
        if arg_q != 0:
            raise NotImplementedError('ddt inside idt/idtmod is not supported')
        qvec[k] += subst[sym]
        ivec[k] += -arg_i
        uvec[k] += -arg_u
        if kind == 'lap':
            continue
        if len(args) > 1:
            ic = resolve(args[1])
            if not _quantity_free(ic):
                raise NotImplementedError(
                    'an idt/idtmod ic must not depend on circuit quantities')
            if kind == 'idtmod':
                m = resolve(args[2]) if len(args) > 2 else sympy.Integer(1)
                o = resolve(args[3]) if len(args) > 3 else sympy.Integer(0)
                ic = ic - m * _wrapfloor((ic - o) / m)
            dc_pins[k] = ic
        if kind == 'idtmod':
            m = resolve(args[2]) if len(args) > 2 else sympy.Integer(1)
            o = resolve(args[3]) if len(args) > 3 else sympy.Integer(0)
            if not (_quantity_free(m) and _quantity_free(o)):
                raise NotImplementedError(
                    'idtmod modulus/offset must not depend on circuit '
                    'quantities')
            periodic.append((k, m, o))

    ## The DC variant: pinned state rows become `s - ic` (fold-into-i, the
    ## elements.py convention), everything else identical.
    ivec_dc = list(ivec)
    uvec_dc = list(uvec)
    for k, ic in dc_pins.items():
        ivec_dc[k] = subst[[s for s, _k2, _a in states
                            if index_of[('state', s.name)] == k][0]] - ic
        uvec_dc[k] = sympy.Integer(0)

    ## `$param_given` symbols actually referenced, in a stable order --
    ## they become trailing arguments of every compiled function and are
    ## bound per instance from `ipar`.
    given_syms = sorted(
        {a for vec in (ivec, qvec, uvec, acvec, ivec_dc)
         for e in vec for a in e.free_symbols
         if str(a).startswith('_hdl_given_')},
        key=lambda a: str(a))
    given_names = [str(a)[len('_hdl_given_'):] for a in given_syms]

    ## ------------------------------------------------------------------
    ## Jacobians and compilation.
    G = sympy.Matrix([ivec]).jacobian(xsyms)
    C = sympy.Matrix([qvec]).jacobian(xsyms)
    G_dc = sympy.Matrix([ivec_dc]).jacobian(xsyms)
    CY = sympy.zeros(n)
    for (r_, c_), psd in CYacc.items():
        CY[r_, c_] = psd
    dudt = [sympy.diff(u_, TIME) for u_ in uvec]

    x = sympy.DeferredVector('x')
    xsubs = dict(zip(xsyms, [x[i] for i in range(n)]))

    NUMPY_MODULES = [{'_wrapfloor': np.floor}, 'numpy']

    def compile_x(expr_or_list, extra=()):
        e = expr_or_list.subs(xsubs) if hasattr(expr_or_list, 'subs') else \
            [e_.subs(xsubs) for e_ in expr_or_list]
        return sympy.lambdify(
            [x] + paramsyms + [TEMP] + given_syms + list(extra), e,
            modules=NUMPY_MODULES, cse=True)

    def compile_t(exprs):
        return sympy.lambdify([TIME] + paramsyms + [TEMP] + given_syms,
                              exprs, modules=NUMPY_MODULES, cse=True)

    funcs = dict(
        i=compile_x(ivec), i_dc=compile_x(ivec_dc),
        q=compile_x(qvec),
        G=compile_x(G), G_dc=compile_x(G_dc),
        C=compile_x(C), CY=compile_x(CY, extra=(FREQ,)),
        u=compile_t(uvec), u_dc=compile_t(uvec_dc), dudt=compile_t(dudt),
        uac=sympy.lambdify(paramsyms + [TEMP] + given_syms, acvec,
                           modules=NUMPY_MODULES, cse=True),
    )

    ## Pure single-device forms for the JAX batched path (eval_i_pure /
    ## eval_q_pure) reuse the SAME symbolic vectors, compiled lazily with
    ## sympy's jax printer on first use (jax may not be installed).
    pure_spec = dict(ivec=ivec, qvec=qvec, xsyms=xsyms, n=n)

    ## The SYMBOLIC forms, kept alongside the compiled ones.  A symbolic
    ## toolkit asked for `i(x)` with symbols in `x` was being served by a
    ## numpy lambda, which works by duck typing for plain arithmetic and
    ## fails the moment the expression contains a `floor`, a `Piecewise`
    ## or anything else numpy evaluates eagerly.  Substituting into these
    ## is exact and cannot fail that way.
    sym_spec = dict(i=ivec, q=qvec, G=list(G), C=list(C),
                    xsyms=xsyms, paramsyms=paramsyms,
                    given_syms=given_syms, shape=(n, n))

    ## ------------------------------------------------------------------
    ## PCNR participation (Aadithya/Keiter/Mei).  A device joins PCNR by
    ## making its limited quantity an explicit unknown and evaluating at
    ## ITS OWN copy, so no two devices can fight over a shared
    ## linearisation point.  The layer needs three things from the device
    ## (`pcnr.augmented_system`): which terminal pair the quantity spans,
    ## the terminal currents as a function of that quantity ALONE, and
    ## their derivative -- because it REMOVES the device's whole i/G
    ## contribution and re-stamps it at the limited voltage.
    ##
    ## So an element qualifies only if its entire current is a function of
    ## one branch voltage, it stores no charge (the layer refuses charge
    ## on a participant in transient: the charge term would have to move
    ## to the limited unknown too), and it introduces no states.  The
    ## limiter itself is generated when the current is recognisably
    ## exponential -- the case limiting exists for -- and the exponential
    ## scale is read off the expression rather than assumed.
    ## PCNR participation, detected PER CONTRIBUTION rather than by
    ## reverse-engineering the assembled `ivec`.  A device joins by making
    ## each limited quantity an explicit unknown and evaluating at ITS OWN
    ## copy, so the layer needs, per junction: the terminal pair it spans,
    ## the terminal currents as a function of that quantity ALONE, and
    ## their derivative -- it REMOVES the device's whole i/G contribution
    ## and re-stamps it at `v_lim`.
    ##
    ## So the rule is: every I-contribution must be a function of its own
    ## branch voltage alone, and at least one must be exponential (the
    ## case limiting exists for).  A device may own SEVERAL such branches
    ## -- a two-junction device is the point of this pass -- and each gets
    ## its own limited unknown, which is exactly what removes the clash
    ## PCNR is about.  Anything else (a V-contribution, a generated state,
    ## a current spanning two branch voltages) disqualifies the element:
    ## better to declare nothing than to claim a capability falsely.
    pcnr_spec = None
    if not states and not vbranches and len(terminalnodes) >= 2:
        cands, ok = [], True
        for st in statements:
            if st.lhs.quantity != 'I':
                ok = False
                break
            b = st.lhs.branch_or_node
            kp = index_of[('node', b.plus.name)]
            km = index_of[('node', b.minus.name)]
            ipart, _q, _u, _nz, _ac = _split_terms(resolve(st.rhs),
                                                   set(xsyms))
            if ipart == 0:
                continue                     # pure charge/noise: harmless
            vsym = sympy.Symbol('_pcnr_v', real=True)
            f = sympy.expand(ipart.subs(xsyms[kp], vsym + xsyms[km]))
            if f.free_symbols & set(xsyms):
                ok = False                   # not a function of V(b) alone
                break
            scales = set()
            for ex in f.atoms(sympy.exp):
                a = sympy.diff(ex.args[0], vsym)
                if a != 0 and not (a.free_symbols & {vsym}):
                    scales.add(sympy.simplify(1 / a))
            if len(scales) != 1:
                ok = False                   # no single exponential scale
                break
            VT_eff = scales.pop()
            cands.append(dict(terminals=(kp, km), vsym=vsym, f=f,
                              dfdv=sympy.diff(f, vsym), VT=VT_eff,
                              IS=sympy.simplify(
                                  VT_eff * sympy.diff(f, vsym).subs(vsym, 0))))
        if ok and cands:
            pcnr_spec = cands

    pcnr_funcs = None
    if pcnr_spec is not None:
        pcnr_funcs = []
        for spec_j in pcnr_spec:
            vs_ = spec_j['vsym']
            ## The traced backend calls these under `jit`, so they need a
            ## jax-printed twin; built lazily so importing the module does
            ## not require jax.
            _sym_j = dict(f=spec_j['f'], dfdv=spec_j['dfdv'], vsym=vs_)
            pcnr_funcs.append(dict(
                sym=_sym_j,
                terminals=spec_j['terminals'],
                i=sympy.lambdify([vs_] + paramsyms + [TEMP], spec_j['f'],
                                 modules=NUMPY_MODULES, cse=True),
                didv=sympy.lambdify([vs_] + paramsyms + [TEMP],
                                    spec_j['dfdv'], modules=NUMPY_MODULES,
                                    cse=True),
                VT=sympy.lambdify(paramsyms + [TEMP], spec_j['VT'],
                                  modules=NUMPY_MODULES),
                IS=sympy.lambdify(paramsyms + [TEMP], spec_j['IS'],
                                  modules=NUMPY_MODULES),
            ))

    state_meta = dict(
        statenames=statenodenames,
        dc_pins={k: sympy.lambdify(paramsyms, ic, modules=NUMPY_MODULES)
                 for k, ic in dc_pins.items()},
        periodic=[(k, sympy.lambdify(paramsyms, m, modules=NUMPY_MODULES),
                   sympy.lambdify(paramsyms, o, modules=NUMPY_MODULES))
                  for k, m, o in periodic],
        has_time_dependence=any(u_.has(TIME) for u_ in uvec),
    )

    ## An x-independent G/C is a CONSTANT stamp: a hand-written linear
    ## element returns its stored matrix by reference, and recomputing one
    ## per Newton iteration is the whole measured gap (benchmarks/
    ## hdl_overhead.py: 28x on a resistor's G before this).  Recorded here,
    ## cached per parameter set in the metaclass.
    xset_ = set(xsyms)
    const_G = not any((e.free_symbols & xset_) for e in G)
    const_C = not any((e.free_symbols & xset_) for e in C)

    cross_spec = None
    if crossings:
        cross_exprs = [resolve(st.expr).subs(xsubs) for st in crossings]
        cross_spec = dict(
            f=sympy.lambdify([x] + paramsyms + [TEMP] + given_syms,
                             cross_exprs, modules=NUMPY_MODULES, cse=True),
            directions=[st.direction for st in crossings])

    limit_spec = [((int(k[0]), int(k[1])),
                   sympy.lambdify(paramsyms + [TEMP], IS_l,
                                  modules=NUMPY_MODULES),
                   sympy.lambdify(paramsyms + [TEMP], VT_l,
                                  modules=NUMPY_MODULES))
                  for k, IS_l, VT_l in limits]

    branchpairs = [branch_key(br) for br in vbranches]
    internalnames = [nd.name for nd in internalnodes]

    return dict(terminalnames=terminalnames, paramnames=paramnames,
                funcs=funcs, pure_spec=pure_spec, state_meta=state_meta,
                branchpairs=branchpairs, internalnames=internalnames,
                const_G=const_G, const_C=const_C, pcnr_funcs=pcnr_funcs,
                given_names=given_names, limit_spec=limit_spec,
                cross_spec=cross_spec, sym_spec=sym_spec,
                has_ac=any(e != 0 for e in acvec))


def _params_of(self):
    ## Values from the RESOLVED iparv, plus the givenness flags from
    ## `ipar` -- givenness is a property of what the user wrote, and only
    ## `ipar` records that (see ParameterDict.update_values).
    vals = [getattr(self.iparv, name) for name in self._hdl_paramnames]
    ipar = self.ipar
    vals += [1.0 if ipar.is_given(nm) else 0.0
             for nm in self._hdl_given_names]
    return vals


def _epar_T(epar):
    return getattr(epar, 'T', 300.0)


def _symbolic_eval(self, which, x, epar):
    """Exact substitution for a symbolic toolkit -- see `sym_spec`."""
    spec = self._hdl_info['sym_spec']
    subs = dict(zip(spec['xsyms'], list(x)))
    subs.update(zip(spec['paramsyms'],
                    [getattr(self.iparv, nm)
                     for nm in self._hdl_paramnames]))
    subs[TEMP] = _epar_T(epar)
    ipar = self.ipar
    subs.update({sym: (1.0 if ipar.is_given(nm) else 0.0)
                 for sym, nm in zip(spec['given_syms'],
                                    self._hdl_given_names)})
    vec = [e.subs(subs) for e in spec[which]]
    if which in ('G', 'C'):
        n = spec['shape'][0]
        return sympy.Matrix(n, n, vec)
    return vec


def _args_of(self, epar):
    """The trailing argument list every compiled function expects:
    parameter values, then T, then the givenness flags."""
    n_par = len(self._hdl_paramnames)
    vals = _params_of(self)
    return vals[:n_par] + [_epar_T(epar)] + vals[n_par:]


def _dc(epar):
    return getattr(epar, 'analysis_kind', None) == 'dc'


class BehaviouralMeta(type):
    def __init__(cls, name, bases, dct):
        if 'analog' not in dct:
            return
        info = generate_code(cls)
        funcs = info['funcs']
        cls._hdl_paramnames = info['paramnames']
        cls._hdl_given_names = info['given_names']
        cls._hdl_info = info
        cls.terminals = info['terminalnames']
        ## Branch objects at class level, like every element: the base
        ## Circuit.__init__ counts them into n; plus/minus resolve through
        ## the terminal mapping when the branch spans terminals, and to
        ## internal nodes otherwise.
        cls.branches = tuple(circuit.Branch(circuit.Node(p), circuit.Node(m))
                             for p, m in info['branchpairs'])

        state_meta = info['state_meta']
        if state_meta['dc_pins']:
            cls.IC_KIND = 'state'

        def i(self, x, epar=defaultepar, params_tree=None):
            if getattr(self.toolkit, 'symbolic', False):
                return _symbolic_eval(self, 'i', x, epar)
            f = funcs['i_dc'] if _dc(epar) and state_meta['dc_pins'] \
                else funcs['i']
            return f(x, *_args_of(self, epar))

        const_G, const_C = info['const_G'], info['const_C']

        def update(self, subject):
            ## The iparv observer every element implements (R.update,
            ## elements.py): parameters changed, so any cached constant
            ## stamp is stale.  This is what makes the cache below a bare
            ## attribute test rather than a per-call dict key -- and it
            ## also means a late-resolved parameter expression reaches the
            ## generated code, since values are read from iparv at call
            ## time and the cache is dropped whenever they move.
            self.__dict__.pop('_hdl_Gc', None)
            self.__dict__.pop('_hdl_Cc', None)

        def G(self, x, epar=defaultepar, params_tree=None):
            if getattr(self.toolkit, 'symbolic', False):
                return _symbolic_eval(self, 'G', x, epar)
            dc = _dc(epar) and state_meta['dc_pins']
            if const_G and not dc:
                ## Constant stamp: computed once and handed back by
                ## reference, as a hand-written linear element does
                ## (benchmarks/hdl_overhead.py measured the gap this
                ## closes).  Dropped by update() on any parameter change.
                cached = self.__dict__.get('_hdl_Gc')
                if cached is None:
                    cached = np.asarray(funcs['G'](x, *_args_of(self, epar)))
                    self.__dict__['_hdl_Gc'] = cached
                return cached
            f = funcs['G_dc'] if dc else funcs['G']
            return np.asarray(f(x, *_args_of(self, epar)))

        def q(self, x, epar=defaultepar, params_tree=None):
            if getattr(self.toolkit, 'symbolic', False):
                return _symbolic_eval(self, 'q', x, epar)
            return funcs['q'](x, *_args_of(self, epar))

        def C(self, x, epar=defaultepar, params_tree=None):
            if getattr(self.toolkit, 'symbolic', False):
                return _symbolic_eval(self, 'C', x, epar)
            if const_C:
                cached = self.__dict__.get('_hdl_Cc')
                if cached is None:
                    cached = np.asarray(funcs['C'](x, *_args_of(self, epar)))
                    self.__dict__['_hdl_Cc'] = cached
                return cached
            return np.asarray(funcs['C'](x, *_args_of(self, epar)))

        def CY(self, x, w=0, epar=defaultepar):
            f_hz = np.abs(w) / (2 * np.pi)
            return np.asarray(funcs['CY'](x, *_args_of(self, epar), f_hz))

        def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
            if analysis == 'ac':
                ## ONLY `ac_stim` terms drive a small-signal analysis.
                ## The DC/transient source terms must NOT come through --
                ## injecting a device's bias constants as an AC drive is
                ## the classic ABM bias-leak (sec. 2.4).
                if not info['has_ac']:
                    return np.zeros(self.n)
                return np.asarray(funcs['uac'](*_args_of(self, epar)),
                                  dtype=complex)
            f = funcs['u_dc'] if _dc(epar) and state_meta['dc_pins'] \
                else funcs['u']
            return np.asarray(f(t, *_args_of(self, epar)), dtype=float)

        def dudt(self, t=0.0, epar=defaultepar, analysis=None,
                 params_tree=None):
            return np.asarray(funcs['dudt'](t, *_args_of(self, epar)),
                              dtype=float)

        def state_ic(self):
            out = []
            for k, icf in state_meta['dc_pins'].items():
                out.append((k, float(icf(*_params_of(self)))))
            return out

        def periodic_states(self):
            out = []
            for k, mf, of in state_meta['periodic']:
                m = float(mf(*_params_of(self)))
                if m > 0:
                    out.append((k, m, float(of(*_params_of(self)))))
            return out

        xset2 = set(info['pure_spec']['xsyms'])
        Gmat = sympy.Matrix([info['pure_spec']['ivec']]).jacobian(
            info['pure_spec']['xsyms'])
        cls.linear = (not any((e.free_symbols & xset2) for e in Gmat)
                      ## a wrap or Piecewise is discontinuous even where its
                      ## slope is constant almost everywhere
                      and not any(e.atoms(_wrapfloor) or
                                  e.atoms(sympy.Piecewise)
                                  for e in info['pure_spec']['ivec']))

        cls.update = update
        cls.i, cls.G, cls.q, cls.C, cls.CY = i, G, q, C, CY
        cls.u, cls.dudt = u, dudt
        if state_meta['dc_pins']:
            cls.state_ic = state_ic
        if state_meta['periodic']:
            cls.periodic_states = periodic_states

        ## @cross -- predict where each declared expression crosses zero
        ## and publish it as a breakpoint.  Two accepted points give the
        ## rate; that is all a first-order prediction needs, and a
        ## prediction only has to BRACKET the corner.
        cspec = info['cross_spec']
        if cspec is not None:
            def accept_step(self, t, x, epar=defaultepar, _cs=cspec):
                vals = [float(v) for v in
                        np.atleast_1d(_cs['f'](np.asarray(x, dtype=float),
                                               *_args_of(self, epar)))]
                prev = self.__dict__.get('_hdl_cross')
                self.__dict__['_hdl_cross'] = (
                    (prev[1] if prev else None), (float(t), vals))

            def reset_state(self, epar=None):
                self.__dict__.pop('_hdl_cross', None)

            def next_event(self, t, _cs=cspec):
                hist = self.__dict__.get('_hdl_cross')
                if not hist or hist[0] is None:
                    return self.toolkit.inf
                (t0, v0), (t1, v1) = hist
                if t1 <= t0:
                    return self.toolkit.inf
                best = self.toolkit.inf
                for k, direction in enumerate(_cs['directions']):
                    rate = (v1[k] - v0[k]) / (t1 - t0)
                    if rate == 0.0:
                        continue
                    if direction and np.sign(rate) != direction:
                        continue
                    ## Time from the LAST accepted point to value zero.
                    dt = -v1[k] / rate
                    if dt <= 0.0:
                        continue
                    tc = t1 + dt
                    if tc > t and tc < best:
                        best = tc
                return best

            cls.accept_step = accept_step
            cls.reset_state = reset_state
            cls.next_event = next_event

        ## $limit -- the state-free convention: return a LIMITED COPY of
        ## the sub-vector.  The limiter moves the solution, so nothing is
        ## ever evaluated at a point other than the one it linearises
        ## about, and no private iterate memory is needed (x0 is handed
        ## in).  Only generated when the model asked for it.
        lspec = info['limit_spec']
        if lspec:
            from pycircuit.circuit._limiting import _pnjlim as _pnj

            def limit(self, x, x0, epar=defaultepar, _ls=lspec):
                out = np.array(x, dtype=float, copy=True)
                x0a = np.asarray(x0, dtype=float)
                args = _args_of(self, epar)
                for (ra, rb), isf, vtf in _ls:
                    vnew = float(out[ra] - out[rb])
                    vold = float(x0a[ra] - x0a[rb])
                    vlim = _pnj(vnew, vold, float(vtf(*args)),
                                float(isf(*args)), self.toolkit)
                    ## Move the ANODE, so a device whose junctions share a
                    ## cathode does not have the second undo the first.
                    out[ra] = out[rb] + vlim
                return out

            cls.limit = limit

        ## PCNR: hand the layer the three things it needs, plus the
        ## limiter for this device's own quantity (pcnr.refine's hook).
        ## Generated only for a recognisably exponential single-branch
        ## element -- see generate_code.  With PCNR enabled the device's
        ## limiting is then PCNR's job, which is the point of the method.
        pf = info['pcnr_funcs']
        if pf is not None:
            from pycircuit.circuit._limiting import _pnjlim
            _pn_pcnr = info['paramnames']

            cls.pcnr_junctions = tuple(tuple(f['terminals']) for f in pf)

            def _pcnr_compiled(jn, which, toolkit, _pf=pf,
                               _pn=_pn_pcnr):
                """The numpy form, or its jax twin for a traced call."""
                if not getattr(toolkit, 'jax', False):
                    return _pf[jn]['i' if which == 'i' else 'didv']
                cache = _pf[jn].setdefault('_jax', {})
                if which not in cache:
                    import jax.numpy as _jnp
                    sym = _pf[jn]['sym']
                    expr = sym['f'] if which == 'i' else sym['dfdv']
                    cache[which] = sympy.lambdify(
                        [sym['vsym']]
                        + [sympy.Symbol(q) for q in _pn] + [TEMP],
                        expr, modules=[{'_wrapfloor': _jnp.floor}, 'jax'],
                        cse=True)
                return cache[which]

            def pcnr_i(v, params, epar, toolkit, jn=0, _pn=_pn_pcnr):
                f = _pcnr_compiled(jn, 'i', toolkit)
                cur = f(v, *[params[q] for q in _pn], _epar_T(epar))
                return toolkit.array([cur, -cur])

            def pcnr_didv(v, params, epar, toolkit, jn=0, _pn=_pn_pcnr):
                f = _pcnr_compiled(jn, 'didv', toolkit)
                g = f(v, *[params[q] for q in _pn], _epar_T(epar))
                return toolkit.array([g, -g])

            def pcnr_scales(params, epar, jn=0, _pf=pf, _pn=_pn_pcnr):
                """(VT, IS) for this junction -- what pnjlim needs, and
                what the traced backend cannot read off an expression."""
                args = [params[q] for q in _pn] + [_epar_T(epar)]
                return (float(_pf[jn]['VT'](*args)),
                        float(_pf[jn]['IS'](*args)))

            cls.pcnr_scales = staticmethod(pcnr_scales)

            def pcnr_limit(vnew, vold, params, epar, toolkit, jn=0,
                           _pf=pf, _pn=_pn_pcnr):
                args = [params[q] for q in _pn] + [_epar_T(epar)]
                return _pnjlim(vnew, vold, float(_pf[jn]['VT'](*args)),
                               float(_pf[jn]['IS'](*args)), toolkit)

            cls.pcnr_i = staticmethod(pcnr_i)
            cls.pcnr_didv = staticmethod(pcnr_didv)
            cls.pcnr_limit = staticmethod(pcnr_limit)

        ## eval_i_pure / eval_q_pure: compiled on first use with sympy's
        ## jax printer; the staticmethods exist only if that succeeds, so
        ## the class is admitted to the vmap groups exactly when the pure
        ## forms are real.  Params arrive as a dict on the batched path.
        cls._hdl_pure_cache = {}
        spec = info['pure_spec']
        pnames = info['paramnames']

        def _compiled_pure(which):
            cache = cls._hdl_pure_cache
            if which not in cache:
                import jax.numpy as jnp
                xv = sympy.DeferredVector('x')
                xs = dict(zip(spec['xsyms'],
                              [xv[i2] for i2 in range(spec['n'])]))
                vec = [e.subs(xs) for e in spec[which]]
                cache[which] = sympy.lambdify(
                    [xv] + [sympy.Symbol(p2) for p2 in pnames] + [TEMP],
                    vec, modules=[{'_wrapfloor': jnp.floor}, 'jax'],
                    cse=True)
            return cache[which]

        def eval_i_pure(x, params, epar, toolkit):
            f = _compiled_pure('ivec')
            return toolkit.array(f(x, *[params[p2] for p2 in pnames],
                                   _epar_T(epar)))

        def eval_q_pure(x, params, epar, toolkit):
            f = _compiled_pure('qvec')
            return toolkit.array(f(x, *[params[p2] for p2 in pnames],
                                   _epar_T(epar)))

        try:
            import sympy.printing.numpy as _sp_np
            have_jax = hasattr(_sp_np, 'JaxPrinter')
        except Exception:
            have_jax = False
        if have_jax and not state_meta['dc_pins']:
            ## DC-pinned elements are kept OFF the batched path for now:
            ## the pure form has no epar-driven variant selection, and a
            ## silently-unpinned batched DC is the failure mode item 2 of
            ## the roadmap just removed.
            cls.eval_i_pure = staticmethod(eval_i_pure)
            cls.eval_q_pure = staticmethod(eval_q_pure)


class Behavioural(circuit.Circuit, metaclass=BehaviouralMeta):
    """Behavioural circuit element -- equations in, element out.

    Write the element's terminal behaviour as contribution statements in
    an ``analog()`` staticmethod whose argument names declare the
    terminals; the metaclass compiles them into the full pycircuit element
    interface (``i``/``q``/``G``/``C``/``u``/``dudt``, DC pinning,
    ``uic`` seeding, the JAX pure forms, gauge-shift declarations).

    Example::

        class MyConductor(Behavioural):
            instparams = [Parameter(name='g', desc='Conductance', unit='S',
                                    default=1e-3)]
            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, g * b.V)

    See ``hdl.md`` for the full language: V-contributions, ``ddt``,
    ``idt``/``idtmod``, internal nodes, current probes, and time-dependent
    sources via the ``TIME`` symbol.
    """

    #: Verilog-A's ``aliasparam``: alternative accepted spellings,
    #: ``{alias: canonical}``.  Compact models carry these by the dozen
    #: (217 declarations in the vacask device library) because one device
    #: is known by different parameter names in different foundries'
    #: decks.  Declared per class::
    #:
    #:     aliasparams = {'isat': 'IS', 'js': 'IS'}
    aliasparams = {}

    def __init__(self, *args, **kvargs):
        aliases = type(self).aliasparams
        if aliases:
            for alias, canonical in aliases.items():
                if alias in kvargs:
                    if canonical in kvargs:
                        raise ValueError(
                            '%r was given both as %r and as its alias %r'
                            % (canonical, canonical, alias))
                    kvargs[canonical] = kvargs.pop(alias)
        super().__init__(*args, **kvargs)
        info = getattr(self, '_hdl_info', None)
        if info is None:
            return
        ## Compilation is keyed on `analog` appearing in the class body, so
        ## a SUBCLASS that changes only `instparams` would inherit the
        ## parent's compiled functions -- which read a different parameter
        ## list -- and answer confidently with the wrong numbers.  Caught
        ## here rather than left to produce plausible waveforms.
        declared = [p.name for p in self.instparams]
        if declared != info['paramnames']:
            raise TypeError(
                '%s changed instparams (%r) without redefining analog(), so '
                'it would run its base class\'s compiled code for %r. '
                'Redeclare analog() in this class (a one-line '
                '`analog = staticmethod(Base.analog)` is enough to trigger '
                'recompilation).'
                % (type(self).__name__, declared, info['paramnames']))
        for name in info['internalnames'] + info['state_meta']['statenames']:
            self.add_node(name)

    def next_event(self, t):
        return self.toolkit.inf


def isconstant(expr):
    for atom in expr.atoms():
        if isinstance(atom, Quantity):
            return False
    return True


if __name__ == "__main__":
    import doctest
    doctest.testmod()
