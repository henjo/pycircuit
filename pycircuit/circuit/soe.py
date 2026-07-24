"""Sequence-of-Expressions (SoE) symbolic solve and result object.

For a *fully symbolic* circuit the exploded transfer function ``N(s)/D(s)`` is
exponentially large. The SoE representation avoids that: solve ``A x = b`` by
symbolic Gaussian elimination whose intermediate results are kept as fresh
symbols, so the solution is a small, ordered list of assignments (``O(n^3)``;
linear for a tridiagonal ladder) plus a compact solution vector over those
symbols -- see ``doc/ginac_fully_symbolic.md`` (improvement #3).

The compactness comes from *sharing*: each intermediate symbol is referenced
several times but stored once. Inlining the assignments into a single expression
(:meth:`SoESolution.to_ratio`) destroys that sharing and can grow
super-linearly, so numeric evaluation resolves the assignment list *in order*
(:meth:`SoESolution.eval`) rather than expanding it. ``to_ratio`` is still useful
to bridge small results to the GiNaC ``N(s)/D(s)`` tools
(:func:`pycircuit.circuit._ginac.poly_coeffs`) for poles/zeros.
"""
import numpy as np
import sympy


def _as_matrix(M):
    if isinstance(M, sympy.MatrixBase):
        return M
    arr = np.asarray(M, dtype=object)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return sympy.Matrix(arr.tolist())


def solve_soe(A, b):
    """Solve ``A x = b`` as a sequence of expressions.

    Returns an :class:`SoESolution` holding an ordered list of ``(symbol, expr)``
    assignments (each ``expr`` small, in terms of earlier symbols and the
    original entries) and the solution vector in terms of those symbols.
    """
    A = _as_matrix(A).copy()
    b = list(_as_matrix(b))
    n = A.rows
    assigns = []
    counter = [0]

    def fresh(expr):
        # Trivial results (atoms/symbols) need no new assignment.
        if expr.is_Atom:
            return expr
        t = sympy.Symbol('t%d' % counter[0])
        counter[0] += 1
        assigns.append((t, expr))
        return t

    # Forward elimination (structural; no pivoting for the symbolic case).
    for k in range(n):
        piv = A[k, k]
        for i in range(k + 1, n):
            if A[i, k] == 0:
                continue
            m = fresh(A[i, k] / piv)
            for j in range(k, n):
                A[i, j] = fresh(A[i, j] - m * A[k, j])
            b[i] = fresh(b[i] - m * b[k])

    # Back substitution.
    x = [None] * n
    for i in range(n - 1, -1, -1):
        acc = b[i]
        for j in range(i + 1, n):
            if A[i, j] != 0 and x[j] != 0:
                acc = fresh(acc - A[i, j] * x[j])
        x[i] = fresh(acc / A[i, i])
    return SoESolution(assigns, x)


class SoESolution:
    """Solution of ``A x = b`` as a sequence of expressions (see :func:`solve_soe`).

    Attributes
    ----------
    assignments : list of (Symbol, Expr)
        Ordered intermediate assignments; resolve them in order to evaluate.
    solution : list of Expr
        Node solution ``x[i]`` in terms of the assignment symbols.
    """

    def __init__(self, assignments, solution):
        self.assignments = assignments
        self.solution = list(solution)
        self.n = len(self.solution)

    def __len__(self):
        return self.n

    # -- numeric evaluation (sequential; preserves the shared subexpressions) --

    def _symbol_table(self):
        """Name -> original symbol, matching the assumptions on the SoE symbols.

        Assignment/solution expressions carry symbols with assumptions (e.g.
        ``Symbol('R0', positive=True)``); resolving caller params by name to
        these exact symbols keeps ``lambdify`` lookups consistent whether the
        caller passes symbols or plain names.
        """
        tset = {t for t, _ in self.assignments}
        syms = set()
        for _, e in self.assignments:
            syms |= e.free_symbols
        for xi in self.solution:
            syms |= xi.free_symbols
        return {s.name: s for s in syms if s not in tset}

    def eval(self, params):
        """Evaluate every node voltage over a parameter sweep.

        ``params`` maps each original symbol (object or name) to a scalar or
        array; the arrays broadcast to a common shape ``S``. Returns a complex
        array of shape ``(n,) + S`` -- the node voltages ``x[0..n-1]``. The
        assignment list is resolved *in order* with :func:`sympy.lambdify`
        (numpy), so the compact shared form is evaluated directly rather than
        expanded.
        """
        table = self._symbol_table()
        env = {}
        for k, v in params.items():
            sym = k if isinstance(k, sympy.Symbol) and str(k) not in table \
                else table.get(str(k), sympy.Symbol(str(k)))
            env[sym] = np.asarray(v, dtype=complex)
        shape = np.broadcast_shapes(*[v.shape for v in env.values()]) \
            if env else ()

        def evaluate(expr):
            syms = sorted(expr.free_symbols, key=str)
            f = sympy.lambdify(syms, expr, modules='numpy')
            val = f(*[env[s] for s in syms])
            return np.broadcast_to(np.asarray(val, dtype=complex), shape)

        for t, expr in self.assignments:
            env[t] = evaluate(expr)
        return np.array([evaluate(xi) for xi in self.solution])

    def eval_tf(self, i, j, params):
        """Evaluate the transfer function ``x[i]/x[j]`` over a parameter sweep."""
        vals = self.eval(params)
        return vals[i] / vals[j]

    # -- symbolic bridge to the exploded N/D world (loses sharing; keep small) --

    def _inline(self, expr):
        for t, e in reversed(self.assignments):
            expr = expr.xreplace({t: e})
        return expr

    def to_ratio(self, i, j=None):
        """Inline the assignments into a single sympy expression.

        Returns ``x[i]`` (or the transfer function ``x[i]/x[j]`` when ``j`` is
        given). Inlining discards the SoE's subexpression sharing, so the result
        can grow super-linearly -- use it for small circuits or to hand a
        transfer function to the GiNaC ``N(s)/D(s)`` tools (poles/zeros). For
        numeric work prefer :meth:`eval`.
        """
        xi = self._inline(self.solution[i])
        if j is None:
            return xi
        return xi / self._inline(self.solution[j])
