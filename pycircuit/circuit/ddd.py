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

import numpy as np
import sympy


__all__ = ['DDD', 'DDDVertex', 'SExpandedDDD', 'DDDFamily', 'DDDCombination',
           'NumericTerminal',
           'ONE', 'ZERO', 'ddd_of_matrix', 'ddd_cramer', 'ddd_cofactor_solve',
           's_expand', 'HierarchicalDDD', 'eval_roots',
           'reverse_cuthill_mckee', 'suppression_order', 'DDDSizeError']


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


class NumericTerminal:
    """A terminal carrying a collapsed sub-determinant.

    A plain DDD has only the ``0`` and ``1`` leaves, so every factor of every
    product term becomes a vertex -- including the ones that carry no symbol we
    care about.  Allowing terminals to hold *values* lets a whole parameter-free
    minor collapse into a single leaf, which is Pi & Shi's multi-terminal DDD
    (DAC 2000).

    It matters for more than size.  Those sub-determinants are the ones whose
    entries are actual component values, and multiplying them out symbolically is
    what turns numeric components into ever-growing exact rationals -- the
    failure that capped the GiNaC backend.  Evaluating the minor once, to a
    number or a small polynomial, means those products are never formed.

    Attributes:
        value: The sub-determinant, already evaluated.
    """

    __slots__ = ('value',)
    is_terminal = True

    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return '<NumericTerminal %s>' % (self.value,)


def _resolve(value, env):
    """Substitute ``env`` into a payload, keeping exact values exact.

    Coercing an exact Integer or Rational to machine complex would quietly turn
    an exact symbolic result into a floating-point one, so only values that are
    already inexact are converted.
    """
    if getattr(value, 'free_symbols', None):
        value = value.subs(env)
    if getattr(value, 'free_symbols', None) or not getattr(value, 'is_number', False):
        return value
    if value.has(sympy.Float):
        return complex(value)
    return value


def eval_roots(roots, env=None):
    """Evaluate several diagrams in one pass over their shared structure.

    Calling :meth:`DDD.eval` per diagram re-walks and re-substitutes everything
    they have in common, which is most of it when the diagrams come from one
    family -- the cost becomes quadratic in exactly the case sharing was meant to
    make cheap.  This resolves each payload once and memoises across all roots.

    Args:
        roots: Vertices or terminals to evaluate.
        env: Substitution.

    Returns:
        Dict mapping ``id(root)`` to its value.
    """
    env = env or {}
    memo = {}
    values = {}
    for root in roots:
        stack = [(root, False)]
        while stack:
            node, expanded = stack.pop()
            if node.is_terminal:
                if id(node) not in memo:
                    memo[id(node)] = _resolve(node.value, env)
                continue
            if id(node) in memo:
                continue
            if not expanded:
                stack.append((node, True))
                stack.append((node.one_edge, False))
                stack.append((node.zero_edge, False))
                continue
            key = id(node.entry)
            if key not in values:
                values[key] = _resolve(node.entry, env)
            memo[id(node)] = (node.sign * values[key] * memo[id(node.one_edge)]
                              + memo[id(node.zero_edge)])
    return {id(r): memo[id(r)] for r in roots}


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

    def term_count(self):
        """Number of product terms the diagram stands for (its 1-paths).

        This is the size the expression would have if expanded, so it is the
        right quantity to guard on before converting to sympy.  Computing it is
        cheap -- one pass over the DAG -- precisely because the diagram shares:
        the count can be astronomically larger than ``size``.
        """
        memo = {}
        stack = [(self.root, False)]
        while stack:
            node, expanded = stack.pop()
            if node.is_terminal:
                ## A collapsed minor is one leaf, so one term -- which is the
                ## point: it stands in for what would have been many.
                memo.setdefault(id(node), 0 if node.value == 0 else 1)
                continue
            if id(node) in memo:
                continue
            if not expanded:
                stack.append((node, True))
                stack.append((node.one_edge, False))
                stack.append((node.zero_edge, False))
                continue
            memo[id(node)] = memo[id(node.one_edge)] + memo[id(node.zero_edge)]
        return memo[id(self.root)]

    def to_sympy(self, max_terms=5000):
        """Expand the diagram into a sympy expression.

        This is the operation the representation exists to avoid, so it is
        guarded: for a fully symbolic circuit the expanded form can have more
        terms than there are atoms in anything worth counting.

        Args:
            max_terms: Refuse if the expansion would exceed this many product
                terms.

        Raises:
            DDDSizeError: If the expansion is too large.  Evaluate the diagram
                numerically, or work with the s-expanded coefficients, instead.
        """
        n = self.term_count()
        if n > max_terms:
            raise DDDSizeError(
                'expanding this diagram would give %d product terms, above the '
                '%d limit; evaluate it numerically or use the s-expanded '
                'coefficients instead of flattening it' % (n, max_terms))
        return self.eval()

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
        for v in self.vertices():
            ## Payloads are shared objects (a matrix cell, or one s-coefficient
            ## of it), so identity is a sound cache key and avoids re-substituting
            ## the same expression once per vertex that references it.
            key = id(v.entry)
            if key not in values:
                values[key] = _resolve(v.entry, env)
        memo = {}

        ## Iterative post-order.  Depth is bounded by the matrix dimension, but
        ## an explicit stack keeps evaluation independent of the recursion limit.
        stack = [(self.root, False)]
        while stack:
            node, expanded = stack.pop()
            if node.is_terminal:
                ## A collapsed minor may itself depend on the frequency, so its
                ## value needs the same substitution as a vertex payload.
                if id(node) not in memo:
                    memo[id(node)] = _resolve(node.value, env)
                continue
            if id(node) in memo:
                continue
            if not expanded:
                stack.append((node, True))
                stack.append((node.one_edge, False))
                stack.append((node.zero_edge, False))
                continue
            memo[id(node)] = (node.sign * values[id(node.entry)]
                              * memo[id(node.one_edge)]
                              + memo[id(node.zero_edge)])
        return memo[id(self.root)]

    ## -- approximation ---------------------------------------------------

    def _entry_values(self, env):
        """Numeric value of each vertex payload; raises if any stays symbolic."""
        values = {}
        for v in self.vertices():
            key = id(v.entry)
            if key not in values:
                e = v.entry
                if getattr(e, 'free_symbols', None):
                    e = e.subs(env)
                try:
                    values[key] = complex(e)
                except (TypeError, ValueError):
                    raise ValueError(
                        'entry %s is still symbolic after substitution; '
                        'dominant-term extraction ranks terms by magnitude and '
                        'so needs numeric values' % (v.entry,))
        return values

    def _best_completion(self, values):
        """Largest product magnitude reachable from each vertex to the 1 terminal.

        This is what makes the search a best-first one rather than an
        enumeration: multiplying the magnitude accumulated so far by this bound
        gives the best any completion of a partial term could achieve, so terms
        come out in decreasing order without generating the rest.
        """
        best = {}
        stack = [(self.root, False)]
        while stack:
            node, expanded = stack.pop()
            if node.is_terminal:
                best.setdefault(id(node), abs(complex(values[id(node)]))
                                if id(node) in values else abs(complex(node.value)))
                continue
            if id(node) in best:
                continue
            if not expanded:
                stack.append((node, True))
                stack.append((node.one_edge, False))
                stack.append((node.zero_edge, False))
                continue
            take = abs(values[id(node.entry)]) * best[id(node.one_edge)]
            best[id(node)] = max(take, best[id(node.zero_edge)])
        return best

    def iter_terms(self, env=None):
        """Yield product terms in decreasing order of magnitude.

        Each term is one 1-path of the diagram -- one product in the expanded
        determinant.  They are produced lazily and *in order*, so taking the
        first few costs a fraction of enumerating all of them: that is what makes
        approximation feasible on a diagram standing for millions of terms.

        Args:
            env: Substitution giving every symbol a numeric value.  May instead
                be a *sequence* of substitutions -- parameter corners -- in which
                case terms are ranked by their worst-case magnitude across all of
                them.  Ranking at a single nominal point can drop a term that
                dominates elsewhere in the design space (Yu & Sechen, 1996).

        Yields:
            ``(expression, values)`` -- the term as sympy, and its value at each
            corner (a list, in the order given).

        Raises:
            ValueError: If any entry remains symbolic after substitution.
        """
        import heapq

        envs = self._as_envs(env)
        if self.root is ZERO:
            return
        values = [self._entry_values(e) for e in envs]
        bests = [self._best_completion(v) for v in values]

        def bound(node, acc):
            return max(abs(a) * b[id(node)] for a, b in zip(acc, bests))

        start_acc = tuple(1.0 + 0j for _ in envs)
        counter = 0
        heap = [(-bound(self.root, start_acc), counter, self.root, start_acc, None)]
        while heap:
            neg, _, node, acc, taken = heapq.heappop(heap)
            if -neg == 0.0:
                continue                       # nothing reachable this way
            if node.is_terminal:
                if node.value == 0:
                    continue
                chain, factors = taken, []
                while chain is not None:
                    v, chain = chain
                    factors.append(v.sign * v.entry)
                factors.append(node.value)
                yield sympy.Mul(*reversed(factors)), list(acc)
                continue

            take = tuple(a * node.sign * v[id(node.entry)]
                         for a, v in zip(acc, values))
            b = bound(node.one_edge, take)
            if b > 0.0:
                counter += 1
                heapq.heappush(heap, (-b, counter, node.one_edge, take,
                                      (node, taken)))
            b = bound(node.zero_edge, acc)
            if b > 0.0:
                counter += 1
                heapq.heappush(heap, (-b, counter, node.zero_edge, acc, taken))

    @staticmethod
    def _as_envs(env):
        if env is None:
            return [{}]
        if isinstance(env, dict):
            return [env]
        return list(env)

    def approximate(self, env=None, tol=0.01, max_terms=200):
        """Keep the dominant terms until the value is within ``tol``.

        Terms are added in decreasing magnitude until the partial sum matches the
        exact value to a relative error of ``tol`` -- at *every* corner, if
        several were given.

        The criterion is the **error**, not the fraction of magnitude
        accumulated.  Product terms of a determinant routinely cancel, so a set
        of terms carrying most of the magnitude need not carry most of the value.

        How much this prunes depends entirely on how far apart the component
        values are; see ``doc/src/circuit/ddd.rst``.  With uniform components
        there is little to drop.

        Args:
            env: Numeric substitution, or a sequence of them (corners).
            tol: Target relative error.
            max_terms: Stop after this many terms even if ``tol`` is unmet.

        Returns:
            ``(expression, n_terms, relative_error)`` -- the error being the
            worst across corners.
        """
        envs = self._as_envs(env)
        ## Validate before evaluating, so a missing substitution reports what is
        ## actually wrong rather than failing on a float conversion later.
        for e in envs:
            self._entry_values(e)
        exact = [complex(self.eval(e)) for e in envs]
        total = [0.0 + 0j] * len(envs)
        kept = []
        for expr, values in self.iter_terms(env):
            kept.append(expr)
            total = [t + v for t, v in zip(total, values)]
            err = max(abs(t - x) / abs(x) if abs(x) > 0 else abs(t)
                      for t, x in zip(total, exact))
            if err <= tol or len(kept) >= max_terms:
                break
        else:
            err = max(abs(t - x) / abs(x) if abs(x) > 0 else abs(t)
                      for t, x in zip(total, exact)) if kept else 1.0
        return sympy.Add(*kept), len(kept), err

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

        extra = {}

        def ref(node):
            if node is ONE:
                return 'ONE'
            if node is ZERO:
                return 'ZERO'
            if node.is_terminal:
                name = 't%d' % len(extra)
                key = id(node)
                if key not in extra:
                    extra[key] = (name, node.value)
                return extra[key][0]
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
        for name, value in extra.values():
            lines.append('  %s [shape=square, label="%s"];'
                         % (name, sympy.sstr(value)[:18]))
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

    def __init__(self, matrix, order='min-degree', keep_symbolic=None,
                 collapse_max_dim=8):
        self.M = matrix
        self.order = order
        self.keep_symbolic = None if keep_symbolic is None else set(keep_symbolic)
        self.collapse_max_dim = collapse_max_dim
        self.cache = {}
        self._nonzero = {}
        self._collapsed = 0

    def collapsible(self, rows, cols):
        """True if this minor carries none of the symbols we must keep.

        Such a minor is just a value, so it can be evaluated once instead of
        contributing a vertex per factor per product term.
        """
        ## Collapsing evaluates a determinant, which is itself expensive above a
        ## certain size -- past that the collapse costs more than the vertices it
        ## saves, so large parameter-free minors are still expanded and only
        ## their sub-minors collapse.
        if self.keep_symbolic is None or not 2 <= len(rows) <= self.collapse_max_dim:
            return False
        for i in rows:
            for j in cols:
                free = getattr(self.M[i, j], 'free_symbols', None)
                if free and (free & self.keep_symbolic):
                    return False
        return True

    def collapse(self, rows, cols):
        """Evaluate a parameter-free minor into a single terminal."""
        sub = self.M.extract(list(rows), list(cols))
        value = sympy.expand(sub.det())
        self._collapsed += 1
        if value == 0:
            return ZERO
        if value == 1:
            return ONE
        return NumericTerminal(value)

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

        if self.collapsible(rows, cols):
            node = self.collapse(rows, cols)
            self.cache[key] = node
            return node

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

        if self.order == 'markowitz':
            ## Degree alone leaves many ties, and breaking them by index lets the
            ## arbitrary node numbering decide -- which is why the same circuit
            ## renumbered can give a diagram several times larger.  Score instead
            ## by what the expansion costs downstream: a row of degree d branches
            ## d ways, and each branch removes one of the intersecting columns, so
            ## weigh the row by the degrees of the lines it will consume.
            best, choice = None, None
            for p, r in enumerate(rows):
                cost = (row_deg[p] - 1) * sum(
                    col_deg[q] - 1 for q, c in enumerate(cols)
                    if self.nonzero(r, c))
                key = (cost, row_deg[p])
                if best is None or key < best:
                    best, choice = key, ('row', p)
            for q, c in enumerate(cols):
                cost = (col_deg[q] - 1) * sum(
                    row_deg[p] - 1 for p, r in enumerate(rows)
                    if self.nonzero(r, c))
                key = (cost, col_deg[q])
                if key < best:
                    best, choice = key, ('col', q)
            return choice

        best_row = min(range(len(rows)), key=lambda p: row_deg[p])
        best_col = min(range(len(cols)), key=lambda q: col_deg[q])
        if row_deg[best_row] <= col_deg[best_col]:       # ties go to the row
            return ('row', best_row)
        return ('col', best_col)


def reverse_cuthill_mckee(A):
    """A node numbering that clusters each row's nonzeros near the diagonal.

    Diagram size turns out to depend sharply on how the nodes happen to be
    numbered -- the same circuit renumbered can give a diagram several times
    larger -- because the expansion heuristic breaks ties by index.  Reordering
    to a banded form first removes that dependence: minors stay contiguous, so
    the same minors recur and the sharing table hits.

    This is the standard Cuthill-McKee sweep, reversed: start from a lowest-degree
    node, visit each level's neighbours in increasing degree, reverse the result.
    Implemented here rather than taken from scipy to keep this module free of
    anything but sympy.

    Args:
        A: Square sympy ``Matrix``; only its sparsity pattern is used.

    Returns:
        List of original indices in their new order.
    """
    n = A.rows
    adjacency = [set() for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (A[i, j] != 0 or A[j, i] != 0):
                adjacency[i].add(j)
                adjacency[j].add(i)
    degree = [len(a) for a in adjacency]

    seen = [False] * n
    order = []
    ## Several passes: a circuit matrix need not be structurally connected.
    while len(order) < n:
        remaining = [i for i in range(n) if not seen[i]]
        start = min(remaining, key=lambda i: (degree[i], i))
        seen[start] = True
        queue = [start]
        while queue:
            node = queue.pop(0)
            order.append(node)
            nbrs = sorted((j for j in adjacency[node] if not seen[j]),
                          key=lambda j: (degree[j], j))
            for j in nbrs:
                seen[j] = True
            queue.extend(nbrs)
    order.reverse()
    return order


def ddd_of_matrix(A, order='auto', keep_symbolic=None,
                  collapse_max_dim=8):
    """Build the determinant decision diagram of a sympy matrix.

    Args:
        A: Square sympy ``Matrix``.  Entries may be arbitrary expressions; each
            becomes one vertex payload (see the module docstring on the compact
            symbol convention).
        order: ``'auto'`` (default) builds both as given and band-reordered and
            keeps the smaller, which makes the result nearly independent of how
            the nodes happen to be numbered -- see
            :func:`reverse_cuthill_mckee`.  ``'min-degree'`` picks the sparsest
            row or column of each minor on the fly; ``'markowitz'`` scores ties by
            the degrees of the lines an expansion consumes; ``'row'`` always
            expands the first row, which is optimal for full matrices and is what
            makes the ``n*2**(n-1)`` identity hold exactly.
        keep_symbolic: Symbols that must stay symbolic.  When given, any minor
            containing none of them is evaluated once into a single terminal
            (see :class:`NumericTerminal`) instead of being expanded -- the
            semi-symbolic mode, where most components are numbers and only a few
            parameters are of interest.  ``None`` (default) expands everything.

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
    if order not in ('auto', 'min-degree', 'row', 'markowitz'):
        raise ValueError('unknown order %r' % (order,))

    if order == 'auto':
        ## Build both as given and band-reordered, and keep the smaller.  A
        ## favourable node numbering (nodes added in signal order, say) can beat
        ## the reordering, while a poor one is far worse than it -- so trying both
        ## keeps the good case and bounds the bad one.  Construction is
        ## milliseconds and the diagram is used many times over.
        candidates = [_build_ddd(A, 'min-degree', keep_symbolic,
                                 collapse_max_dim)]
        perm = reverse_cuthill_mckee(A)
        if perm != list(range(A.rows)):
            banded = _build_ddd(A.extract(perm, perm), 'min-degree',
                                keep_symbolic, collapse_max_dim)
            candidates.append(banded)
        return min(candidates, key=lambda d: d.size)

    return _build_ddd(A, order, keep_symbolic, collapse_max_dim)


def _build_ddd(A, order, keep_symbolic, collapse_max_dim):
    builder = _Builder(A, order, keep_symbolic=keep_symbolic,
                       collapse_max_dim=collapse_max_dim)
    root = builder.build(tuple(range(A.rows)), tuple(range(A.cols)))
    result = DDD(root, A, order, len(builder.cache))
    result.n_collapsed = builder._collapsed
    return result


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


## ------------------------------------------------- s-expanded (multiroot) --

class SExpandedDDD:
    """The determinant split by powers of the frequency variable.

    A network function is wanted as a rational polynomial in ``s``,
    ``N(s)/D(s)``, because that is what yields poles, zeros and interpretable
    coefficients.  Expanding a symbolic determinant into that form term by term
    is hopeless -- the count of ``s``-expanded product terms dwarfs even the
    determinant's own -- but the *coefficients* can be kept as diagrams.

    This is Shi & Tan's multiroot DDD (TCAD 2001): one root per power of ``s``,
    with the roots sharing subgraphs.  Sharing works because the memo key is
    ``(rows, cols, power)``, so any minor-and-power pair reached from several
    coefficients is built once.

    Attributes:
        roots: One graph per power of ``s``; ``roots[k]`` is the coefficient of
            ``s**k``.  Trailing zero coefficients are trimmed.
        matrix, s, order: Provenance.
        n_minors: Size of the sharing table.
    """

    def __init__(self, roots, matrix, s, order, n_minors):
        self.roots = roots
        self.matrix = matrix
        self.s = s
        self.order = order
        self.n_minors = n_minors
        self._size = None

    @property
    def degree(self):
        """Degree in ``s`` (``-1`` for an identically zero determinant)."""
        return len(self.roots) - 1

    def coefficient(self, k):
        """The coefficient of ``s**k`` as a `DDD`."""
        if not 0 <= k < len(self.roots):
            return DDD(ZERO, self.matrix, self.order, 0)
        return DDD(self.roots[k], self.matrix, self.order, self.n_minors)

    @property
    def size(self):
        """Distinct vertices across *all* coefficients, counting sharing once.

        Summing ``coefficient(k).size`` would double-count every shared vertex,
        which is the whole point of the representation, so it is counted over the
        union instead.
        """
        if self._size is None:
            seen, stack = set(), list(self.roots)
            while stack:
                v = stack.pop()
                if v.is_terminal or id(v) in seen:
                    continue
                seen.add(id(v))
                stack.append(v.one_edge)
                stack.append(v.zero_edge)
            self._size = len(seen)
        return self._size

    def eval_coeffs(self, env=None):
        """Evaluate every coefficient, lowest power of ``s`` first."""
        return [self.coefficient(k).eval(env) for k in range(len(self.roots))]

    def roots_of(self, env=None):
        """Numeric roots of the polynomial, i.e. the poles of ``1/det``.

        Args:
            env: Substitution grounding every symbol other than ``s``.

        Returns:
            ``numpy`` array of complex roots.

        Raises:
            ValueError: If any coefficient is still symbolic after ``env``.

        Note:
            Root-finding from *expanded* coefficients is ill-conditioned when
            they span many orders of magnitude, which they routinely do for a
            circuit -- this is a property of the polynomial form, not of the
            diagram.  Where a factored form is available it will be better
            behaved; see ``doc/ddd_conclusions.md`` §8.6.
        """
        import numpy as np

        coeffs = self.eval_coeffs(env)
        bad = [c for c in coeffs if getattr(c, 'free_symbols', None)]
        if bad:
            raise ValueError('coefficients still symbolic; substitute values '
                             'for %s first' % bad[0].free_symbols)
        ## numpy.roots wants highest power first.
        return np.roots([complex(c) for c in reversed(coeffs)])

    def approximate(self, env=None, tol=0.01, max_terms=200):
        """Approximate each coefficient of ``s`` separately.

        Approximating the determinant as a whole weights the powers of ``s`` by
        whatever the frequency happens to be, so a coefficient that dominates the
        response in one band gets discarded for another.  Approximating *after*
        the s-expansion treats every power on an equal footing, which is why
        Shi & Tan (TCAD 2001) require this order -- and it is what leaves the
        result usable across the whole frequency range rather than at one point.

        Returns:
            List of ``(expression, n_terms, relative_error)``, one per power of
            ``s``, lowest first.
        """
        return [self.coefficient(k).approximate(env, tol, max_terms)
                for k in range(len(self.roots))]

    def dominant_poles(self, env=None):
        """Estimate poles as ratios of consecutive coefficients.

        For a polynomial whose roots are widely separated -- the usual situation
        in an analog circuit, where a dominant pole is often orders of magnitude
        below the rest -- the ``k``-th pole is approximately
        ``-d[k-1] / d[k]``.

        This is what makes a big circuit *interpretable*: each estimate is a
        ratio of two symbolic coefficients, so it can be read and reasoned
        about, where the exact root of a degree-15 polynomial cannot be written
        down at all.  It is also better conditioned than root-finding, since no
        polynomial is solved.

        Args:
            env: Numeric substitution; omit for symbolic ratios.

        Returns:
            List of estimates, lowest-frequency pole first.
        """
        coeffs = self.eval_coeffs(env)
        out = []
        for k in range(1, len(coeffs)):
            lower, upper = coeffs[k - 1], coeffs[k]
            if upper == 0:
                continue
            out.append(-lower / upper)
        return out

    def __repr__(self):
        return '<SExpandedDDD degree=%d size=%d order=%s>' % (
            self.degree, self.size, self.order)


class _SBuilder(_Builder):
    """Builds coefficient diagrams, memoised on ``(rows, cols, power)``.

    The recurrence is the ordinary Laplace expansion with the entry split into
    its powers of ``s``.  With ``entry = e0 + e1*s + ...``, the coefficient of
    ``s**k`` in a minor's determinant is::

        coeff_k(M) = sum_terms sign * sum_d  e_d * coeff_(k-d)(sub-minor)

    so one vertex is emitted per ``(entry, power)`` pair rather than per entry.
    """

    def __init__(self, matrix, s, order='min-degree'):
        super().__init__(matrix, order)
        self.s = s
        self._coeffs = {}

    def entry_coeffs(self, i, j):
        """Coefficients of ``matrix[i, j]`` in ``s``, lowest power first."""
        key = (i, j)
        hit = self._coeffs.get(key)
        if hit is None:
            e = sympy.expand(self.M[i, j])
            if e.has(self.s):
                try:
                    poly = sympy.Poly(e, self.s)
                except sympy.PolynomialError as exc:
                    ## A 1/s entry -- an inductor stamped in admittance form,
                    ## typically.  The coefficient split is only defined for
                    ## polynomial entries, so say so plainly rather than letting
                    ## sympy's error surface.
                    raise ValueError(
                        'entry (%d,%d) = %s is not polynomial in %s, so it has '
                        'no s-power expansion (%s)' % (i, j, e, self.s, exc))
                hit = list(reversed(poly.all_coeffs()))
            else:
                hit = [e]
            self._coeffs[key] = hit
        return hit

    def max_degree(self):
        """Upper bound on the determinant's degree in ``s``.

        Every product term takes one entry from each row, so the degree cannot
        exceed the sum over rows of the highest degree in that row.
        """
        total = 0
        for i in range(self.M.rows):
            best = 0
            for j in range(self.M.cols):
                if self.nonzero(i, j):
                    best = max(best, len(self.entry_coeffs(i, j)) - 1)
            total += best
        return total

    def build_coeff(self, rows, cols, k):
        """Diagram for the coefficient of ``s**k`` in ``det(M[rows, cols])``."""
        if k < 0:
            return ZERO
        if not rows:
            ## An empty minor has determinant 1 -- entirely in the s^0 term.
            return ONE if k == 0 else ZERO

        key = (rows, cols, k)
        hit = self.cache.get(key)
        if hit is not None:
            return hit

        pivot = self._select(rows, cols)
        if pivot is None:
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

        ## Emit one vertex per (entry, power-of-s) pair, back to front so each
        ## 0-edge chains to the next alternative.
        node = ZERO
        for (p, q, i, j) in reversed(terms):
            sub_rows = tuple(x for x in rows if x != i)
            sub_cols = tuple(x for x in cols if x != j)
            sign = -1 if (p + q) % 2 else 1
            coeffs = self.entry_coeffs(i, j)
            for d in reversed(range(len(coeffs))):
                if coeffs[d] == 0 or d > k:
                    continue
                child = self.build_coeff(sub_rows, sub_cols, k - d)
                if child is ZERO:
                    continue                     # contributes nothing; skip it
                node = DDDVertex(i, j, coeffs[d], sign, child, node)

        self.cache[key] = node
        return node


def s_expand(A, s, order='min-degree'):
    """Split ``det(A)`` into coefficients of powers of ``s``, as shared diagrams.

    Args:
        A: Square sympy ``Matrix`` whose entries are polynomial in ``s`` (an MNA
            matrix ``G + s*C`` is degree 1 in every entry).
        s: The frequency symbol.
        order: Expansion ordering, as for :func:`ddd_of_matrix`.

    Returns:
        SExpandedDDD: one root per power of ``s``, sharing subgraphs.

    Raises:
        ValueError: If ``A`` is not square, or an entry is not polynomial in
            ``s`` (a ``1/s`` entry, say) -- the coefficient split is only defined
            for polynomial entries.

    Example:
        >>> import sympy
        >>> s, g, c = sympy.symbols('s g c')
        >>> A = sympy.Matrix([[g + s*c, -g], [-g, g + s*c]])
        >>> E = s_expand(A, s)
        >>> E.degree
        2
        >>> sympy.simplify(E.coefficient(2).eval() - c**2) == 0
        True
    """
    A = sympy.Matrix(A)
    if A.rows != A.cols:
        raise ValueError('DDD needs a square matrix, got %dx%d' % (A.rows, A.cols))

    builder = _SBuilder(A, s, order)
    ## Validate every entry up front: a clear error before any work beats one
    ## surfacing halfway through a build.
    for i in range(A.rows):
        for j in range(A.cols):
            if builder.nonzero(i, j):
                builder.entry_coeffs(i, j)

    rows = tuple(range(A.rows))
    cols = tuple(range(A.cols))
    roots = [builder.build_coeff(rows, cols, k)
             for k in range(builder.max_degree() + 1)]
    while roots and roots[-1] is ZERO:            # trim trailing zero powers
        roots.pop()
    return SExpandedDDD(roots, A, s, order, len(builder.cache))


## ------------------------------------------------- shared cofactor solve --

class DDDCombination:
    """A linear combination of diagrams: ``sum_k coeff_k * diagram_k``.

    A Cramer numerator expanded along its substituted column is a weighted sum
    of *cofactors* of the original matrix -- and cofactors are minors of that
    same matrix, so they come out of the same memo table as the determinant.
    Keeping the numerator as a combination rather than one diagram is what lets
    a whole family of transfer functions share one construction.
    """

    __slots__ = ('parts',)

    def __init__(self, parts):
        self.parts = list(parts)                 # [(coefficient, DDD), ...]

    def eval(self, env=None):
        ## The coefficients are matrix entries in their own right, so they need
        ## the same substitution as the diagrams they multiply -- evaluating the
        ## diagram alone would leave a half-substituted expression behind.
        total = 0
        for coeff, diagram in self.parts:
            total = total + _resolve(coeff, env or {}) * diagram.eval(env)
        return total

    def roots(self):
        return [d.root for d in (p[1] for p in self.parts)]

    def term_count(self):
        return sum(d.term_count() for _, d in self.parts)

    def __repr__(self):
        return '<DDDCombination parts=%d>' % len(self.parts)


def _distinct_vertices(roots):
    seen, stack = set(), list(roots)
    while stack:
        v = stack.pop()
        if v.is_terminal or id(v) in seen:
            continue
        seen.add(id(v))
        stack.append(v.one_edge)
        stack.append(v.zero_edge)
    return seen


class DDDFamily:
    """A determinant and its Cramer numerators, sharing a single memo table.

    Building each numerator as its own diagram (:func:`ddd_cramer`) re-expands
    minors that the determinant already covered, because a substituted column
    makes it a different matrix.  Expanding along that column instead expresses
    every numerator through cofactors of the *original* matrix, so the whole
    family is built once.

    This is the structure behind Shi & Tan's observation (TCAD 2001) that noise
    analysis -- which needs one transfer function per noise source -- costs
    little more than a single transfer function: the transimpedances share
    almost everything.

    Attributes:
        denominator: The network determinant, as a `DDD`.
        size: Distinct vertices across everything built, counting sharing once.
    """

    def __init__(self, A, b, order='min-degree'):
        self.matrix = sympy.Matrix(A)
        if self.matrix.rows != self.matrix.cols:
            raise ValueError('DDD needs a square matrix, got %dx%d'
                             % (self.matrix.rows, self.matrix.cols))
        self.b = sympy.Matrix(b).reshape(self.matrix.rows, 1)
        self.order = order

        self._builder = _Builder(self.matrix, order)
        self._rows = tuple(range(self.matrix.rows))
        self._cols = tuple(range(self.matrix.cols))
        self._numerators = {}
        self._cofactors = []

        self.denominator = DDD(self._builder.build(self._rows, self._cols),
                               self.matrix, order, len(self._builder.cache))

    def numerator(self, i):
        """Cramer numerator for unknown ``i``, as a `DDDCombination`.

        Expanded along the substituted column: ``sum_k b[k] * (-1)^(k+i) *
        M_ki``, where ``M_ki`` deletes row ``k`` and column ``i``.  Only the
        nonzero entries of ``b`` contribute, so a sparse right-hand side -- the
        usual case -- costs very few cofactors.
        """
        if i in self._numerators:
            return self._numerators[i]

        parts = []
        cols = tuple(c for c in self._cols if c != i)
        for k in range(self.matrix.rows):
            bk = self.b[k]
            if bk == 0:
                continue
            rows = tuple(r for r in self._rows if r != k)
            root = self._builder.build(rows, cols)
            if root is ZERO:
                continue
            sign = -1 if (k + i) % 2 else 1
            parts.append((sign * bk,
                          DDD(root, self.matrix, self.order,
                              len(self._builder.cache))))
        combination = DDDCombination(parts)
        self._numerators[i] = combination
        return combination

    def cofactor(self, k, i):
        """Diagram for the minor deleting row ``k`` and column ``i``.

        Cofactors do not depend on the right-hand side, so one family serves
        every solve against the same matrix -- which is what lets a subcircuit be
        suppressed towards all of its terminals at once.
        """
        rows = tuple(r for r in self._rows if r != k)
        cols = tuple(c for c in self._cols if c != i)
        root = self._builder.build(rows, cols)
        self._cofactors.append(root)
        return DDD(root, self.matrix, self.order, len(self._builder.cache))

    @property
    def size(self):
        roots = [self.denominator.root] + self._cofactors
        for combination in self._numerators.values():
            roots.extend(combination.roots())
        return len(_distinct_vertices(roots))

    def eval_node(self, i, env=None):
        """Value of unknown ``i`` -- numerator over the shared determinant."""
        return self.numerator(i).eval(env) / self.denominator.eval(env)

    def __repr__(self):
        return '<DDDFamily n=%d built=%d size=%d>' % (
            self.matrix.rows, len(self._numerators), self.size)


def ddd_cofactor_solve(A, b, indices=None, order='min-degree'):
    """Solve ``A x = b`` with every numerator sharing one diagram construction.

    Equivalent in value to :func:`ddd_cramer`, but far cheaper when several
    unknowns are wanted, because the numerators are expressed as cofactors of
    ``A`` rather than as determinants of ``n`` different matrices.

    Args:
        A: Square sympy ``Matrix``.
        b: Right-hand side.
        indices: Unknowns to build numerators for; ``None`` builds all.
        order: Expansion ordering.

    Returns:
        ``(family, numerators)`` -- the `DDDFamily` (carrying the shared
        denominator and the size accounting) and a dict of `DDDCombination`.

    Example:
        >>> import sympy
        >>> a, b_, c, d = sympy.symbols('a b c d')
        >>> fam, num = ddd_cofactor_solve(sympy.Matrix([[a, b_], [c, d]]), [1, 0])
        >>> sympy.simplify(num[0].eval() - d) == 0
        True
    """
    family = DDDFamily(A, b, order=order)
    if indices is None:
        indices = range(family.matrix.cols)
    return family, {i: family.numerator(i) for i in indices}


## ------------------------------------------------------------ hierarchy --

class HierarchicalDDD:
    """A determinant computed by suppressing subcircuits before the rest.

    Instead of expanding one flat determinant, eliminate a subcircuit's internal
    unknowns and stamp what is left into the remaining system -- Tan & Shi's
    hierarchical DDD (TCAD 2000), in the "symbolic stamp" formulation of Xu, Shi
    & Li (ASP-DAC 2011) that matches how a
    :class:`~pycircuit.circuit.SubCircuit` already composes.

    With unknowns partitioned into internal ``i`` and terminal ``t``, the Schur
    complement gives ``det(A) = det(A_ii)·det(A_tt - A_ti A_ii^-1 A_it)``.
    Clearing the inverse keeps everything polynomial::

        M = det(A_ii)*A_tt - A_ti * adj(A_ii) * A_it
        det(A) = det(M) / det(A_ii)**(m-1)         (m = dim of M)

    and every entry of ``adj(A_ii)`` is a cofactor of ``A_ii``, so a whole block
    comes out of one shared construction (`DDDFamily`).

    **Several blocks are suppressed in turn**, which is what makes it worth
    doing.  A single split has to choose between leaving too many terminals (the
    reduced system is nearly the original) and eliminating too much at once (the
    block's own cofactor family costs more than the flat diagram).  Suppressing a
    sequence of small blocks avoids both: no large matrix is ever built, because
    each level hands the next one something already reduced.  That is the
    structure behind the published 56x on a µA741, where the leaves are three and
    four nodes across.

    Args:
        A: Square sympy ``Matrix``.
        blocks: Either one iterable of indices (a single suppression) or a
            sequence of **disjoint** index sets, eliminated in the order given --
            leaves first, then whatever they fed.  Indices are always in the
            original numbering.
        order: Expansion ordering for every diagram built.

    Attributes:
        levels: One record per suppression, each carrying its `DDDFamily`.
        top: Diagram of whatever remained after the last suppression.
        size: Vertices across every level -- a sum, not a product.
    """

    def __init__(self, A, blocks, order='auto'):
        A = sympy.Matrix(A)
        if A.rows != A.cols:
            raise ValueError('DDD needs a square matrix, got %dx%d'
                             % (A.rows, A.cols))
        n = A.rows

        blocks = self._normalise(blocks, n)
        self.matrix = A
        self.blocks = blocks
        self.order = order
        self.levels = []

        current = A
        alive = list(range(n))                   # original index of each row
        counter = 0

        for block in blocks:
            positions = [alive.index(i) for i in block]
            rest = [p for p in range(len(alive)) if p not in set(positions)]
            if not rest:
                raise ValueError('the last block would eliminate everything; '
                                 'leave at least one unknown')
            current, family, mapping = self._suppress(
                current, positions, rest, order, counter)
            counter += 1
            self.levels.append({'family': family, 'blocks': mapping,
                                'm': current.rows})
            alive = [alive[p] for p in rest]

        self.top = ddd_of_matrix(current, order=order)
        self.reduced = current
        self.alive = alive

    @staticmethod
    def _normalise(blocks, n):
        """Accept one index set or a sequence of them; validate disjointness."""
        blocks = list(blocks)
        if blocks and all(isinstance(b, (int, np.integer)) for b in blocks):
            blocks = [blocks]
        out, seen = [], set()
        for block in blocks:
            block = tuple(sorted(set(int(i) for i in block)))
            if not block:
                raise ValueError('empty block')
            if any(i < 0 or i >= n for i in block):
                raise ValueError('block index out of range')
            if seen & set(block):
                raise ValueError('blocks must be disjoint')
            seen.update(block)
            out.append(block)
        if not out:
            raise ValueError('at least one block is required')
        if len(seen) >= n:
            raise ValueError('blocks must leave at least one unknown')
        return out

    @staticmethod
    def _suppress(M, positions, rest, order, level):
        """One Schur step: eliminate ``positions``, return the reduced matrix.

        The reduced matrix's entries are fresh symbols standing for combinations
        of the block's cofactors.  Keeping them opaque is what lets the next
        level treat this one as an ordinary matrix -- a vertex payload is never
        inspected, only evaluated.
        """
        Aii = M.extract(positions, positions)
        Ait = M.extract(positions, rest)
        Ati = M.extract(rest, positions)
        Att = M.extract(rest, rest)

        family = DDDFamily(Aii, sympy.zeros(len(positions), 1), order=order)
        if family.denominator.root is ZERO:
            raise ValueError(
                'block %d is structurally singular, so it cannot be suppressed; '
                'choose a partition whose internal sub-matrix is non-singular'
                % level)

        ni, m = len(positions), len(rest)
        cof, mapping = {}, {}
        reduced = sympy.zeros(m, m)
        for a in range(m):
            for b in range(m):
                parts = []
                if Att[a, b] != 0:
                    parts.append((Att[a, b], family.denominator))
                for p in range(ni):
                    if Ati[a, p] == 0:
                        continue
                    for k in range(ni):
                        if Ait[k, b] == 0:
                            continue
                        if (k, p) not in cof:
                            cof[(k, p)] = family.cofactor(k, p)
                        sign = -1 if (k + p) % 2 else 1
                        parts.append((-sign * Ati[a, p] * Ait[k, b], cof[(k, p)]))
                if not parts:
                    continue                     # structurally zero: stay sparse
                sym = sympy.Symbol('_lvl%d_%d_%d' % (level, a, b))
                mapping[sym] = DDDCombination(parts)
                reduced[a, b] = sym
        return reduced, family, mapping

    @property
    def size(self):
        """Vertices across every level, counting each level's sharing once."""
        return sum(lvl['family'].size for lvl in self.levels) + self.top.size

    def eval(self, env=None):
        """Evaluate ``det(A)``, one level at a time.

        Each level is resolved before the next is touched, because a level's
        diagrams are built over the previous level's stamp symbols.
        """
        env = dict(env or {})
        factors = []

        for lvl in self.levels:
            ## One pass over everything this level needs, rather than one walk
            ## per stamp entry -- they overlap almost completely.
            roots = [lvl['family'].denominator.root]
            for combination in lvl['blocks'].values():
                roots.extend(combination.roots())
            memo = eval_roots(roots, env)
            D = memo[id(lvl['family'].denominator.root)]
            if D == 0:
                raise ZeroDivisionError(
                    'an internal block is singular (its determinant is zero), '
                    'so it cannot be suppressed; choose a partition whose '
                    'internal sub-matrix is non-singular')

            ## The stamp is built as ``M = D*A_tt - ...`` because that keeps it
            ## polynomial, which is what a diagram needs.  Numerically, though,
            ## carrying that factor through and dividing by ``D**(m-1)`` at the
            ## end is ruinous: with a couple of dozen unknowns the power spans
            ## tens of orders of magnitude and the answer is lost to
            ## cancellation.  Dividing each entry by D here recovers the actual
            ## Schur complement, so each level contributes a single factor of D
            ## and nothing is ever raised to a power.
            resolved = {}
            for sym, combination in lvl['blocks'].items():
                total = 0
                for coeff, diagram in combination.parts:
                    total = total + _resolve(coeff, env) * memo[id(diagram.root)]
                resolved[sym] = total / D
            factors.append(D)
            env.update(resolved)

        value = self.top.eval(env)
        for D in factors:
            value = value * D
        return value

    def __repr__(self):
        return '<HierarchicalDDD levels=%d remaining=%d size=%d>' % (
            len(self.levels), self.top.matrix.rows, self.size)


def suppression_order(A, keep=(), limit=None):
    """Choose which unknowns to suppress, and in what order.

    Three things decide it, and only the first is a matter of taste.

    **Correctness.** Suppressing a block inverts it, so its sub-matrix must be
    non-singular.  For a single unknown that means a nonzero diagonal, which
    rules out rows that have none -- in MNA the branch equation of a voltage
    source is exactly such a row: it states ``v+ - v- = E`` and carries nothing
    on its own diagonal.  Those rows are skipped.

    **Purpose.** Whatever the answer is asked about has to survive: ports, the
    nodes a transfer function is taken between, anything to be swept or
    differentiated.  Those are named in ``keep``.

    **Cost.** Eliminating an unknown couples all of its neighbours to each other,
    so a low-degree unknown is cheap and a highly-connected one is not.  Taking
    them in increasing degree *of the pattern as it fills in* is the classic
    minimum-degree elimination ordering, and it is what keeps each block's
    cofactor family small enough to be worth building.

    Args:
        A: Square sympy ``Matrix``.
        keep: Indices that must not be suppressed.
        limit: Stop after this many; ``None`` suppresses everything eligible.

    Returns:
        List of single-index tuples, ready to pass to `HierarchicalDDD` as
        ``blocks``.

    Note:
        Fill-in is tracked on the sparsity *pattern*, which is cheap.  A pivot
        that is structurally fine can still vanish numerically for particular
        component values; `HierarchicalDDD` reports that when it happens.
    """
    A = sympy.Matrix(A)
    n = A.rows
    keep = set(int(i) for i in keep)

    neighbours = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(n):
            if i != j and (A[i, j] != 0 or A[j, i] != 0):
                neighbours[i].add(j)
                neighbours[j].add(i)

    ## A zero diagonal cannot be a 1x1 pivot.
    eligible = {i for i in range(n)
                if i not in keep and A[i, i] != 0}

    chosen = []
    while eligible and (limit is None or len(chosen) < limit):
        ## Count *all* neighbours, not just the ones still eligible: eliminating
        ## an unknown couples every neighbour to every other, and fill-in landing
        ## on a kept node is just as real as fill-in landing on an eliminable one.
        i = min(eligible, key=lambda k: (len(neighbours[k]), k))
        chosen.append((i,))
        eligible.discard(i)
        ## Eliminating i makes its neighbours mutually adjacent.
        nbrs = neighbours[i]
        for a in nbrs:
            neighbours[a].discard(i)
            for b in nbrs:
                if a != b:
                    neighbours[a].add(b)
        neighbours[i] = set()
    return chosen
