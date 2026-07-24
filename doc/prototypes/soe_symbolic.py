"""Prototype: Sequence-of-Expressions (SoE) symbolic solve.

Demonstrates improvement #3 from ``doc/ginac_fully_symbolic.md``: for fully
symbolic circuits the exploded ``N(s)/D(s)`` transfer function is exponentially
large, but a *sequence of expressions* -- symbolic Gaussian elimination whose
intermediate results are kept as fresh symbols -- stays polynomial (linear for a
tridiagonal ladder) and is directly evaluable.

This is sympy-only to illustrate the representation; in production the SoE
arithmetic and (especially) the numeric evaluation of the assignment list would
run through GiNaC, which is much faster on large symbolic expressions and shares
common subexpressions in its DAG.

Run: ``python doc/prototypes/soe_symbolic.py``
"""
import time
import sympy


def solve_soe(A, b):
    """Solve ``A x = b`` as a sequence of expressions.

    Returns ``(assignments, x)`` where ``assignments`` is an ordered list of
    ``(symbol, expr)`` pairs (each ``expr`` small, in terms of earlier symbols
    and the original entries) and ``x`` is the solution in terms of those
    symbols.  Evaluate by resolving the assignments in order, then ``x``.
    """
    n = A.rows
    A = A.copy()
    b = list(b)
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
    return assigns, x


def _ladder(N):
    s = sympy.Symbol('s')
    R = [sympy.Symbol('R%d' % i, positive=True) for i in range(N)]
    C = [sympy.Symbol('C%d' % i, positive=True) for i in range(N)]
    A = sympy.zeros(N, N)
    b = sympy.zeros(N, 1)
    for i in range(N):
        A[i, i] = 1 / R[i] + s * C[i] + (1 / R[i - 1] if i > 0 else 0)
        if i + 1 < N:
            A[i, i + 1] = -1 / R[i]
            A[i + 1, i] = -1 / R[i]
    b[0] = 1 / R[0]
    return A, b, s, R, C


def main():
    print("Fully-symbolic RC ladder: direct exploded solve vs SoE")
    print("  N | direct ops | SoE ops | SoE assigns | correct")
    for N in (4, 6, 8, 10, 12):
        A, b, s, R, C = _ladder(N)
        x_direct = A.LUsolve(b)
        ops_direct = sympy.count_ops(x_direct[N - 1])
        assigns, x = solve_soe(A, b)
        ops_soe = sum(sympy.count_ops(e) for _, e in assigns)

        # Correctness: resolve the SoE numerically, compare to the direct solve.
        env = {R[i]: 100.0 * (i + 1) for i in range(N)}
        env.update({C[i]: 1e-9 * (i + 1) for i in range(N)})
        env[s] = 1j * 2 * 3.14159 * 1e6
        for t, e in assigns:
            env[t] = complex(e.subs(env))
        v_soe = complex(x[N - 1].subs(env))
        v_dir = complex(x_direct[N - 1].subs(env))
        ok = abs(v_soe - v_dir) < 1e-9 * max(1, abs(v_dir))

        print("  %2d | %6d | %5d | %6d | %s"
              % (N, ops_direct, ops_soe, len(assigns), ok))


if __name__ == '__main__':
    main()
