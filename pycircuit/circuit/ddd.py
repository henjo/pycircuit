# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Determinant Decision Diagrams -- a shared graph representation of a symbolic
determinant.

A symbolic determinant of an ``n x n`` matrix has up to ``n!`` product terms, but
those terms share an enormous amount of structure: the same *minor* is reached
along many different expansion paths.  A DDD makes that sharing explicit.  Each
vertex stands for one ``(matrix entry, minor)`` pair and contributes::

    value(v) = sign(v) * entry(v) * value(v.one_edge) + value(v.zero_edge)

so evaluating the whole determinant costs one multiply and one add per vertex --
and, notably, **no division**, so unlike an elimination-based solver there are no
pivots that can vanish at particular parameter values.

Construction follows Shi, *A Simple Implementation of Determinant Decision
Diagram* (ICCAD 2010), which is considerably easier to implement than the
original TCAD-2000 flow:

* sharing is a plain dict keyed by the minor's ``(row indices, column indices)``
  -- no BDD/ZBDD package is involved;
* a full symbol order is weakened to an *expansion order* chosen on the fly
  (which row or column to expand next), so there is no ordering pre-pass;
* signs fall out during construction, so there is no separate sign pass.

Symbol convention
-----------------
The diagram is built over **matrix entries**: one vertex carries a whole entry
such as ``g1 + g2 + s*(c1 + c2)``.  Shi & Tan (TCAD 2001) call this the *compact
symbol* representation, as against a *full symbol* one where every device gets
its own vertex; the same circuit can differ by ~3x in vertex count between the
two, so any reported size is meaningless without saying which was used.  The
entry expression is retained as the vertex payload, so component-level results
remain available to whatever consumes the graph.

See ``doc/src/circuit/ddd.rst`` for worked examples and rendered diagrams, and
``doc/ddd_conclusions.md`` for why this representation was chosen.
"""

import sympy


__all__ = ['DDD', 'DDDVertex', 'ONE', 'ZERO', 'ddd_of_matrix', 'ddd_cramer',
           'DDDSizeError']


class DDDSizeError(Exception):
    """Raised when a graph is too large for an operation that must stay bounded."""


class DDDVertex:
    """One vertex: an entry of the matrix, taken within a particular minor.

    Attributes:
        row, col: Position of the entry in the *original* matrix.
        entry: The symbolic expression stamped at that position.
        sign: ``+1`` or ``-1``, the Laplace cofactor sign, fixed at construction.
        one_edge: Diagram of the minor left after deleting ``row`` and ``col``
            (the "take this term" branch).
        zero_edge: Diagram of the remaining terms of the same expansion (the
            "skip this term" branch), or `ZERO` at the end of a group.
    """

    __slots__ = ('row', 'col', 'entry', 'sign', 'one_edge', 'zero_edge')

    is_terminal = False

    def __init__(self, row, col, entry, sign, one_edge, zero_edge):
        self.row = row
        self.col = col
        self.entry = entry
        self.sign = sign
        self.one_edge = one_edge
        self.zero_edge = zero_edge

    def __repr__(self):
        return '<DDDVertex a[%d,%d] sign=%+d>' % (self.row, self.col, self.sign)


class _Terminal:
    """The constant leaves: `ONE` (empty minor, determinant 1) and `ZERO`."""

    __slots__ = ('value', 'name')
    is_terminal = True

    def __init__(self, value, name):
        self.value = value
        self.name = name

    def __repr__(self):
        return self.name


ONE = _Terminal(1, 'ONE')
ZERO = _Terminal(0, 'ZERO')


class DDD:
    """A built determinant diagram.

    Attributes:
        root: Root vertex (or a terminal, for a trivial or structurally singular
            matrix).
        matrix: The sympy matrix it was built from.
        order: The expansion-ordering strategy used.
        n_minors: Number of distinct minors expanded -- i.e. how many entries the
            sharing table holds.
    """

    def __init__(self, root, matrix, order, n_minors):
        self.root = root
        self.matrix = matrix
        self.order = order
        self.n_minors = n_minors
        self._size = None

    def __len__(self):
        return self.size

    @property
    def size(self):
        """Vertex count excluding terminals -- ``|DDD|``.

        Doubles as an operation count: evaluating the determinant costs one
        multiply and one add per vertex.
        """
        if self._size is None:
            self._size = len(self.vertices())
        return self._size

    def vertices(self):
        """All distinct vertices, in discovery order."""
        seen, out, stack = set(), [], [self.root]
        while stack:
            v = stack.pop()
            if v.is_terminal or id(v) in seen:
                continue
            seen.add(id(v))
            out.append(v)
            stack.append(v.one_edge)
            stack.append(v.zero_edge)
        return out

    def entries(self):
        """The distinct ``(row, col)`` positions referenced by the graph."""
        return {(v.row, v.col) for v in self.vertices()}

    ## -- evaluation ------------------------------------------------------

    def eval(self, env=None):
        """Evaluate the determinant.

        Args:
            env: Mapping of sympy symbols to numbers.  Each distinct matrix
                entry is substituted once and the graph is then walked, so the
                cost is ``O(distinct entries + |DDD|)`` rather than one
                substitution per product term.

        Returns:
            The determinant -- a number when ``env`` grounds every symbol,
            otherwise a sympy expression.
        """
        env = env or {}
        values = {}
        for (i, j) in self.entries():
            e = self.matrix[i, j]
            if getattr(e, 'free_symbols', None):
                e = e.subs(env)
            try:
                values[(i, j)] = complex(e)
            except (TypeError, ValueError):
                values[(i, j)] = e                       # still symbolic

        memo = {id(ONE): 1, id(ZERO): 0}

        ## Iterative post-order.  Depth is bounded by the matrix dimension, but
        ## an explicit stack keeps evaluation independent of the recursion limit.
        stack = [(self.root, False)]
        while stack:
            node, expanded = stack.pop()
            if node.is_terminal or id(node) in memo:
                continue
            if not expanded:
                stack.append((node, True))
                stack.append((node.one_edge, False))
                stack.append((node.zero_edge, False))
                continue
            memo[id(node)] = (node.sign * values[(node.row, node.col)]
                              * memo[id(node.one_edge)]
                              + memo[id(node.zero_edge)])
        return memo[id(self.root)]

    ## -- rendering -------------------------------------------------------

    def to_dot(self, max_vertices=200, name='DDD'):
        """Return Graphviz DOT source for the diagram.

        Solid edges are 1-edges (take this term and descend into its minor);
        dashed edges are 0-edges (skip to the next term of the same expansion).

        Args:
            max_vertices: Refuse to render a graph larger than this.  The limit
                is deliberately mandatory -- a real circuit's diagram has
                thousands of vertices, would be unreadable, and would make the
                documentation build shell out to ``dot`` on something enormous.
                Rendered examples are meant to be tiny.
            name: Graph name in the DOT source.

        Raises:
            DDDSizeError: If the graph exceeds ``max_vertices``.
        """
        verts = self.vertices()
        if len(verts) > max_vertices:
            raise DDDSizeError(
                'DDD has %d vertices, above the %d limit for rendering; '
                'render a smaller example instead' % (len(verts), max_vertices))

        ids = {id(v): 'v%d' % k for k, v in enumerate(verts)}

        def ref(node):
            if node is ONE:
                return 'ONE'
            if node is ZERO:
                return 'ZERO'
            return ids[id(node)]

        lines = ['digraph %s {' % name,
                 '  rankdir=TB;',
                 '  node [shape=circle, fontsize=10];',
                 '  ONE [shape=square, label="1"];',
                 '  ZERO [shape=square, label="0"];']
        for v in verts:
            label = 'a%d%d' % (v.row, v.col)
            lines.append('  %s [label="%s%s"];'
                         % (ids[id(v)], '-' if v.sign < 0 else '', label))
        for v in verts:
            lines.append('  %s -> %s;' % (ids[id(v)], ref(v.one_edge)))
            lines.append('  %s -> %s [style=dashed];'
                         % (ids[id(v)], ref(v.zero_edge)))
        lines.append('}')
        return '\n'.join(lines)

    def __repr__(self):
        return '<DDD size=%d minors=%d order=%s>' % (
            self.size, self.n_minors, self.order)


## ------------------------------------------------------------- building --

class _Builder:
    """Layered expansion with minor-level sharing (ICCAD 2010).

    The whole mechanism is: to expand a minor, pick one of its rows or columns,
    emit one vertex per nonzero in it, and recurse on the sub-minor each of those
    leaves behind.  A dict keyed by ``(rows, cols)`` means a minor reached a
    second time is never expanded again -- that dict *is* the sharing, and it is
    why no BDD package is needed.
    """

    def __init__(self, matrix, order='min-degree'):
        self.M = matrix
        self.order = order
        self.cache = {}
        self._nonzero = {}

    def nonzero(self, i, j):
        """Structural nonzero test, memoised (sympy comparisons are not cheap)."""
        key = (i, j)
        hit = self._nonzero.get(key)
        if hit is None:
            hit = self.M[i, j] != 0
            self._nonzero[key] = hit
        return hit

    def build(self, rows, cols):
        """Diagram for the minor on ``rows`` x ``cols`` (ascending tuples)."""
        if not rows:
            return ONE                                   # empty minor: det = 1

        key = (rows, cols)
        hit = self.cache.get(key)
        if hit is not None:
            return hit

        pivot = self._select(rows, cols)
        if pivot is None:                                # a zero row or column
            self.cache[key] = ZERO
            return ZERO

        kind, pos = pivot
        if kind == 'row':
            r = rows[pos]
            terms = [(pos, q, r, cols[q]) for q in range(len(cols))
                     if self.nonzero(r, cols[q])]
        else:
            c = cols[pos]
            terms = [(p, pos, rows[p], c) for p in range(len(rows))
                     if self.nonzero(rows[p], c)]

        ## Build the sibling chain back to front, so each vertex's 0-edge is the
        ## next term of the same expansion and the last one terminates at ZERO.
        node = ZERO
        for (p, q, i, j) in reversed(terms):
            sub_rows = tuple(x for x in rows if x != i)
            sub_cols = tuple(x for x in cols if x != j)
            node = DDDVertex(i, j, self.M[i, j],
                             -1 if (p + q) % 2 else 1,
                             self.build(sub_rows, sub_cols), node)

        self.cache[key] = node
        return node

    def _select(self, rows, cols):
        """Choose which row or column of this minor to expand.

        ``min-degree`` takes whichever has fewest nonzeros, keeping the branching
        factor down on sparse circuit matrices.  ``row`` always takes the first
        row: provably optimal for a *full* matrix (Shi, TCAS-II 2010), and what
        makes the ``n * 2**(n-1)`` size identity hold exactly.

        Returns:
            ``('row'|'col', position)``, or ``None`` when some row or column is
            entirely zero -- in which case the determinant is 0.
        """
        row_deg = [sum(1 for c in cols if self.nonzero(r, c)) for r in rows]
        col_deg = [sum(1 for r in rows if self.nonzero(r, c)) for c in cols]
        if min(row_deg) == 0 or min(col_deg) == 0:
            return None
        if self.order == 'row':
            return ('row', 0)
        best_row = min(range(len(rows)), key=lambda p: row_deg[p])
        best_col = min(range(len(cols)), key=lambda q: col_deg[q])
        if row_deg[best_row] <= col_deg[best_col]:       # ties go to the row
            return ('row', best_row)
        return ('col', best_col)


def ddd_of_matrix(A, order='min-degree'):
    """Build the determinant decision diagram of a sympy matrix.

    Args:
        A: Square sympy ``Matrix``.  Entries may be arbitrary expressions; each
            becomes one vertex payload (see the module docstring on the compact
            symbol convention).
        order: ``'min-degree'`` (default) picks the sparsest row or column of
            each minor on the fly; ``'row'`` always expands the first row, which
            is optimal for full matrices.

    Returns:
        DDD: the built diagram.

    Raises:
        ValueError: If ``A`` is not square, or ``order`` is unknown.

    Example:
        >>> import sympy
        >>> a, b, c, d = sympy.symbols('a b c d')
        >>> D = ddd_of_matrix(sympy.Matrix([[a, b], [c, d]]))
        >>> sympy.simplify(D.eval() - (a*d - b*c)) == 0
        True
        >>> D.size
        4
    """
    A = sympy.Matrix(A)
    if A.rows != A.cols:
        raise ValueError('DDD needs a square matrix, got %dx%d' % (A.rows, A.cols))
    if order not in ('min-degree', 'row'):
        raise ValueError('unknown order %r' % (order,))

    builder = _Builder(A, order)
    root = builder.build(tuple(range(A.rows)), tuple(range(A.cols)))
    return DDD(root, A, order, len(builder.cache))


def ddd_cramer(A, b, indices=None, order='min-degree'):
    """Solve ``A x = b`` by Cramer's rule, as diagrams.

    ``x_i = det(A_i) / det(A)``, where ``A_i`` is ``A`` with column ``i``
    replaced by ``b``.  Every part stays a diagram: nothing is expanded into a
    sympy expression, so this is usable on systems whose explicit ``N(s)/D(s)``
    would be astronomically large.

    A transfer function between two unknowns needs no denominator at all --
    ``x_i / x_j = det(A_i) / det(A_j)`` -- since the shared ``det(A)`` cancels.

    Args:
        A: Square sympy ``Matrix``.
        b: Right-hand side (column ``Matrix`` or sequence).
        indices: Which unknowns to build numerators for.  ``None`` builds all of
            them, which is rarely what you want -- pass just the ones needed.
        order: Expansion ordering, as for :func:`ddd_of_matrix`.

    Returns:
        ``(denominator, numerators)`` -- a `DDD` and a dict mapping index to
        `DDD`.

    Example:
        >>> import sympy
        >>> a, b_, c, d = sympy.symbols('a b c d')
        >>> A = sympy.Matrix([[a, b_], [c, d]])
        >>> den, num = ddd_cramer(A, [1, 0], indices=[0])
        >>> sympy.simplify(num[0].eval() - d) == 0
        True
    """
    A = sympy.Matrix(A)
    if A.rows != A.cols:
        raise ValueError('DDD needs a square matrix, got %dx%d' % (A.rows, A.cols))
    bcol = sympy.Matrix(b).reshape(A.rows, 1)
    if indices is None:
        indices = range(A.cols)

    denominator = ddd_of_matrix(A, order=order)
    numerators = {}
    for i in indices:
        Ai = A.copy()
        Ai[:, i] = bcol
        numerators[i] = ddd_of_matrix(Ai, order=order)
    return denominator, numerators
