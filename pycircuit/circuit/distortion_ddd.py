# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Symbolic distortion on shared structure -- stage A of the DDD plan.

``doc/distortion_ddd_plan.md`` stage A, whose question is whether the cost
recorded in ``doc/distortion_mimo_plan.md`` section 8.3 is intrinsic or an
artifact of asking sympy to *simplify*.  Stage E of that plan measured
expression **size** and found it polynomial; it measured nothing about
**cost**, and ``sympy.cancel`` did not finish in 900 s at ``U^7`` on a
two-node system.

**Nothing here modifies the existing implementation.**  ``distortion.py`` is
the oracle this is checked against, and refactoring it to share code with its
own checker would destroy the independence that makes the comparison worth
anything.  The graded ring, the consistency filter and the composition
bookkeeping are reused *unchanged* -- what changes is only what a coefficient
is made of.

**The whole idea in one line:** :class:`GradedSpectrum` asks its coefficients
for ``+``, ``*``, ``conjugate()`` and ``!= 0`` and nothing else, so supplying
a different coefficient type re-runs the same recurrence over a different
algebra.  No new ring, no reimplemented recurrence.

**Why a DAG rather than a flat sum of products.**  Keeping terms as
``(coefficient, factors)`` pairs -- the shape :class:`DDDCombination` uses --
means a product of an ``n``-term and an ``m``-term operand has ``n*m`` terms,
so the expansion this is trying to avoid happens anyway.  Holding an
unexpanded :class:`Expr` graph instead means a product is *one node* pointing
at two operands.  Cost then falls on evaluation, which is memoised over the
shared structure -- the same move that measured 1000x faster than
substitution in the feasibility spike.
"""

import sympy


## Hash-consing table: structurally identical expressions become the *same*
## object, so sharing is by identity and the memo can key on it.  This is
## what buys canonical merging over sympy's structural hashing -- `a*b` and
## `b*a` intern to one node here because the arguments are sorted.
_INTERN = {}


class Expr:
    """A node in a hash-consed expression graph.

    Immutable, interned, and deliberately tiny: the only operations are the
    four the graded ring asks for, plus division, which the linear solve
    needs for its determinant.

    Attributes:
        op: ``'atom'``, ``'add'``, ``'mul'`` or ``'inv'``.
        args: Operand nodes; empty for an atom.
        value: The payload for an atom -- a number or a sympy expression.
    """

    __slots__ = ('op', 'args', 'value', '_hash')

    def __new__(cls, op, args=(), value=None):
        args = tuple(args)
        key = (op, tuple(id(a) for a in args),
               value if op == 'atom' else None)
        found = _INTERN.get(key)
        if found is not None:
            return found
        self = object.__new__(cls)
        self.op = op
        self.args = args
        self.value = value
        self._hash = hash(key)
        _INTERN[key] = self
        return self

    ## --------------------------------------------------------- building --

    @staticmethod
    def atom(value):
        """A leaf: a number, or any sympy expression treated as opaque."""
        if isinstance(value, Expr):
            return value
        return Expr('atom', (), value)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        ## Interning makes structural equality identity.  Comparison against
        ## a plain 0 has to work because GradedSpectrum drops zero terms.
        if isinstance(other, Expr):
            return self is other
        if other == 0:
            return self.op == 'atom' and self.value == 0
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        return result if result is NotImplemented else not result

    def __add__(self, other):
        other = Expr.atom(other)
        if self._is_zero():
            return other
        if other._is_zero():
            return self
        ## Flatten and order: `a + b` and `b + a` must intern to one node, or
        ## the sharing this class exists for is lost to argument order.
        if self._is_number() and other._is_number():
            ## A real circuit is mostly numeric: fold, or the graph fills
            ## with arithmetic that could have been done once at build time.
            return Expr.atom(self._number() + other._number())
        terms = self._flat('add') + other._flat('add')
        return Expr('add', _ordered(terms))

    __radd__ = __add__

    def __mul__(self, other):
        other = Expr.atom(other)
        if self._is_zero() or other._is_zero():
            return ZERO
        if self._is_one():
            return other
        if other._is_one():
            return self
        if self._is_number() and other._is_number():
            return Expr.atom(self._number() * other._number())
        factors = self._flat('mul') + other._flat('mul')
        return Expr('mul', _ordered(factors))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Expr.atom(other)
        if other._is_number():
            return self * Expr.atom(1.0 / other._number())
        return self * Expr('inv', (other,))

    def __neg__(self):
        return self * Expr.atom(-1)

    def __sub__(self, other):
        return self + (-Expr.atom(other))

    def conjugate(self):
        """Required by ``GradedSpectrum.from_phasor`` for the mirrored term."""
        if self.op == 'atom':
            value = self.value
            return Expr.atom(value.conjugate()
                             if hasattr(value, 'conjugate') else value)
        if self.op == 'inv':
            return Expr('inv', (self.args[0].conjugate(),))
        return Expr(self.op, _ordered([a.conjugate() for a in self.args]))

    ## ---------------------------------------------------------- helpers --

    def _flat(self, op):
        return list(self.args) if self.op == op else [self]

    def _is_number(self):
        if self.op != 'atom':
            return False
        value = self.value
        if isinstance(value, sympy.Expr):
            return bool(value.is_number)
        return True

    def _number(self):
        value = self.value
        return complex(value) if not isinstance(value, complex) else value

    def _is_zero(self):
        return self.op == 'atom' and self.value == 0

    def _is_one(self):
        return self.op == 'atom' and self.value == 1

    def __repr__(self):
        if self.op == 'atom':
            return 'Expr(%r)' % (self.value,)
        return 'Expr(%s, %d args)' % (self.op, len(self.args))


def _ordered(nodes):
    """Canonical argument order, so commuting operands intern together."""
    return sorted(nodes, key=id)


ZERO = Expr.atom(0)
ONE = Expr.atom(1)


def nodes_of(*roots):
    """Distinct nodes reachable from ``roots`` -- the representation's size.

    Counting the *graph* rather than the tree is the whole point: a tree
    count would report the size the representation exists to avoid.
    """
    seen, stack = set(), [r for r in roots if isinstance(r, Expr)]
    while stack:
        node = stack.pop()
        if id(node) in seen:
            continue
        seen.add(id(node))
        stack.extend(node.args)
    return len(seen)


def evaluate(roots, env):
    """Evaluate several graphs in one memoised pass over shared structure.

    The counterpart of ``ddd.eval_roots``, and for the same reason: walking
    each root separately re-walks everything they have in common, which is
    most of it, so the cost becomes quadratic in exactly the case sharing was
    meant to make cheap.

    Args:
        roots: :class:`Expr` nodes to evaluate.
        env: Mapping from sympy symbols to numbers.

    Returns:
        Dict from ``id(root)`` to its value.
    """
    memo = {}
    atoms = {}

    for root in roots:
        if not isinstance(root, Expr):
            continue
        stack = [(root, False)]
        while stack:
            node, expanded = stack.pop()
            if id(node) in memo:
                continue
            if node.op == 'atom':
                value = node.value
                key = id(node)
                if key not in atoms:
                    ## Bare symbols resolve by dict lookup; anything else
                    ## pays a substitution *once*, not once per occurrence.
                    if isinstance(value, sympy.Expr) and not value.is_number:
                        if value in env:
                            atoms[key] = env[value]
                        else:
                            atoms[key] = complex(value.subs(env))
                    else:
                        atoms[key] = complex(value)
                memo[id(node)] = atoms[key]
                continue
            if not expanded:
                stack.append((node, True))
                stack.extend((a, False) for a in node.args)
                continue
            if node.op == 'add':
                total = 0
                for a in node.args:
                    total = total + memo[id(a)]
                memo[id(node)] = total
            elif node.op == 'mul':
                total = 1
                for a in node.args:
                    total = total * memo[id(a)]
                memo[id(node)] = total
            elif node.op == 'inv':
                memo[id(node)] = 1.0 / memo[id(node.args[0])]
            else:
                raise ValueError('unknown op %r' % (node.op,))

    return {id(r): memo[id(r)] for r in roots if isinstance(r, Expr)}


def evaluate_one(root, env):
    """Convenience wrapper around :func:`evaluate` for a single graph."""
    if not isinstance(root, Expr):
        return complex(root)
    return evaluate([root], env)[id(root)]
