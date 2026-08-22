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

    iterms, qterms, uterms, nterms = [], [], [], []
    for term in terms:
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
            nterms)


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

    ## ------------------------------------------------------------------
    ## Pass 1 -- inventory: nodes, branches with current unknowns (every
    ## V-contributed or I-probed branch), and idt/idtmod states.
    nodes = set()
    vbranches = []      # branches that get a current unknown, in order
    ibranch_keys = set()

    def branch_key(branch):
        return (branch.plus.name, branch.minus.name)

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

    internalnodes = sorted(nodes - set(terminalnodes), key=lambda n: n.name)

    ## States: each distinct idt/idtmod APPLICATION is one state.  Walk in
    ## a deterministic order.
    states = []          # (state_symbol, kind, args) with kind 'idt'/'idtmod'
    state_subst = {}
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

    ## ------------------------------------------------------------------
    ## Pass 2 -- assemble the residual vectors.
    ivec = [sympy.Integer(0)] * n
    qvec = [sympy.Integer(0)] * n
    uvec = [sympy.Integer(0)] * n
    CYacc = {}

    for st in statements:
        b = st.lhs.branch_or_node
        p = index_of[('node', b.plus.name)]
        m_ = index_of[('node', b.minus.name)]
        if st.lhs.quantity == 'I':
            ipart, qpart, upart, npart = _split_terms(resolve(st.rhs),
                                                      set(xsyms))
            ivec[p] += ipart; ivec[m_] -= ipart
            qvec[p] += qpart; qvec[m_] -= qpart
            uvec[p] += upart; uvec[m_] -= upart
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
            ipart, qpart, upart, npart = _split_terms(resolve(st.rhs),
                                                      set(xsyms))
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
        arg_i, arg_q, arg_u, arg_n = _split_terms(resolve(args[0]),
                                                   set(xsyms))
        if arg_n:
            raise NotImplementedError('noise inside idt/idtmod is not '
                                      'supported')
        if arg_q != 0:
            raise NotImplementedError('ddt inside idt/idtmod is not supported')
        qvec[k] += subst[sym]
        ivec[k] += -arg_i
        uvec[k] += -arg_u
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
        return sympy.lambdify([x] + paramsyms + [TEMP] + list(extra), e,
                              modules=NUMPY_MODULES, cse=True)

    def compile_t(exprs):
        return sympy.lambdify([TIME] + paramsyms + [TEMP], exprs,
                              modules=NUMPY_MODULES, cse=True)

    funcs = dict(
        i=compile_x(ivec), i_dc=compile_x(ivec_dc),
        q=compile_x(qvec),
        G=compile_x(G), G_dc=compile_x(G_dc),
        C=compile_x(C), CY=compile_x(CY, extra=(FREQ,)),
        u=compile_t(uvec), u_dc=compile_t(uvec_dc), dudt=compile_t(dudt),
    )

    ## Pure single-device forms for the JAX batched path (eval_i_pure /
    ## eval_q_pure) reuse the SAME symbolic vectors, compiled lazily with
    ## sympy's jax printer on first use (jax may not be installed).
    pure_spec = dict(ivec=ivec, qvec=qvec, xsyms=xsyms, n=n)

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
    pcnr_spec = None
    if (not states and not vbranches and len(terminalnodes) >= 2
            and all(e == 0 for e in qvec)):
        nz = [k for k in range(n) if ivec[k] != 0]
        if len(nz) == 2:
            k_p, k_m = nz
            if sympy.simplify(ivec[k_p] + ivec[k_m]) == 0:
                vsym = sympy.Symbol('_pcnr_v', real=True)
                f = ivec[k_p].subs(xsyms[k_p], vsym + xsyms[k_m])
                f = sympy.simplify(sympy.expand(f))
                if not (f.free_symbols & set(xsyms)):
                    ## A function of the branch voltage alone.  Now the
                    ## exponential scale: pnjlim needs (VT, IS), and both
                    ## follow from the expression -- VT from the argument
                    ## of the exponential, IS from VT*f'(0), which is
                    ## exactly IS for the textbook IS*(exp(v/VT)-1).
                    scales = set()
                    for ex in f.atoms(sympy.exp):
                        a = sympy.diff(ex.args[0], vsym)
                        if a != 0 and not (a.free_symbols & {vsym}):
                            scales.add(sympy.simplify(1 / a))
                    if len(scales) == 1:
                        VT_eff = scales.pop()
                        IS_eff = sympy.simplify(
                            VT_eff * sympy.diff(f, vsym).subs(vsym, 0))
                        pcnr_spec = dict(
                            terminals=(k_p, k_m), vsym=vsym, f=f,
                            dfdv=sympy.diff(f, vsym),
                            VT=VT_eff, IS=IS_eff)

    pcnr_funcs = None
    if pcnr_spec is not None:
        vs_ = pcnr_spec['vsym']
        pcnr_funcs = dict(
            terminals=pcnr_spec['terminals'],
            i=sympy.lambdify([vs_] + paramsyms + [TEMP], pcnr_spec['f'],
                             modules=NUMPY_MODULES, cse=True),
            didv=sympy.lambdify([vs_] + paramsyms + [TEMP],
                                pcnr_spec['dfdv'], modules=NUMPY_MODULES,
                                cse=True),
            VT=sympy.lambdify(paramsyms + [TEMP], pcnr_spec['VT'],
                              modules=NUMPY_MODULES),
            IS=sympy.lambdify(paramsyms + [TEMP], pcnr_spec['IS'],
                              modules=NUMPY_MODULES),
        )

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

    branchpairs = [branch_key(br) for br in vbranches]
    internalnames = [nd.name for nd in internalnodes]

    return dict(terminalnames=terminalnames, paramnames=paramnames,
                funcs=funcs, pure_spec=pure_spec, state_meta=state_meta,
                branchpairs=branchpairs, internalnames=internalnames,
                const_G=const_G, const_C=const_C, pcnr_funcs=pcnr_funcs)


def _params_of(self):
    return [getattr(self.iparv, name) for name in self._hdl_paramnames]


def _epar_T(epar):
    return getattr(epar, 'T', 300.0)


def _dc(epar):
    return getattr(epar, 'analysis_kind', None) == 'dc'


class BehaviouralMeta(type):
    def __init__(cls, name, bases, dct):
        if 'analog' not in dct:
            return
        info = generate_code(cls)
        funcs = info['funcs']
        cls._hdl_paramnames = info['paramnames']
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
            f = funcs['i_dc'] if _dc(epar) and state_meta['dc_pins'] \
                else funcs['i']
            return f(x, *_params_of(self), _epar_T(epar))

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
            dc = _dc(epar) and state_meta['dc_pins']
            if const_G and not dc:
                ## Constant stamp: computed once and handed back by
                ## reference, as a hand-written linear element does
                ## (benchmarks/hdl_overhead.py measured the gap this
                ## closes).  Dropped by update() on any parameter change.
                cached = self.__dict__.get('_hdl_Gc')
                if cached is None:
                    cached = np.asarray(funcs['G'](x, *_params_of(self),
                                                   _epar_T(epar)))
                    self.__dict__['_hdl_Gc'] = cached
                return cached
            f = funcs['G_dc'] if dc else funcs['G']
            return np.asarray(f(x, *_params_of(self), _epar_T(epar)))

        def q(self, x, epar=defaultepar, params_tree=None):
            return funcs['q'](x, *_params_of(self), _epar_T(epar))

        def C(self, x, epar=defaultepar, params_tree=None):
            if const_C:
                cached = self.__dict__.get('_hdl_Cc')
                if cached is None:
                    cached = np.asarray(funcs['C'](x, *_params_of(self),
                                                   _epar_T(epar)))
                    self.__dict__['_hdl_Cc'] = cached
                return cached
            return np.asarray(funcs['C'](x, *_params_of(self),
                                         _epar_T(epar)))

        def CY(self, x, w=0, epar=defaultepar):
            f_hz = np.abs(w) / (2 * np.pi)
            return np.asarray(funcs['CY'](x, *_params_of(self),
                                          _epar_T(epar), f_hz))

        def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
            if analysis == 'ac':
                ## No behavioral AC excitation yet -- and letting the DC/
                ## transient source terms through would inject them as
                ## spurious AC drives (the classic ABM bias-leak).
                return np.zeros(self.n)
            f = funcs['u_dc'] if _dc(epar) and state_meta['dc_pins'] \
                else funcs['u']
            return np.asarray(f(t, *_params_of(self), _epar_T(epar)),
                              dtype=float)

        def dudt(self, t=0.0, epar=defaultepar, analysis=None,
                 params_tree=None):
            return np.asarray(funcs['dudt'](t, *_params_of(self),
                                            _epar_T(epar)), dtype=float)

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

        ## PCNR: hand the layer the three things it needs, plus the
        ## limiter for this device's own quantity (pcnr.refine's hook).
        ## Generated only for a recognisably exponential single-branch
        ## element -- see generate_code.  With PCNR enabled the device's
        ## limiting is then PCNR's job, which is the point of the method.
        pf = info['pcnr_funcs']
        if pf is not None:
            from pycircuit.circuit._limiting import _pnjlim
            _pn_pcnr = info['paramnames']

            cls.pcnr_junctions = (tuple(pf['terminals']),)

            def pcnr_i(v, params, epar, toolkit, _pf=pf, _pn=_pn_pcnr):
                cur = _pf['i'](v, *[params[q] for q in _pn], _epar_T(epar))
                return toolkit.array([cur, -cur])

            def pcnr_didv(v, params, epar, toolkit, _pf=pf, _pn=_pn_pcnr):
                g = _pf['didv'](v, *[params[q] for q in _pn], _epar_T(epar))
                return toolkit.array([g, -g])

            def pcnr_limit(vnew, vold, params, epar, toolkit,
                           _pf=pf, _pn=_pn_pcnr):
                args = [params[q] for q in _pn] + [_epar_T(epar)]
                return _pnjlim(vnew, vold, float(_pf['VT'](*args)),
                               float(_pf['IS'](*args)), toolkit)

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

    def __init__(self, *args, **kvargs):
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
