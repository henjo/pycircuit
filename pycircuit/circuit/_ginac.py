"""Python bridge to the GiNaC linear-solve extension (:mod:`_ginac_ext`).

Converts a sympy linear system to GiNaC-parseable strings, solves fraction-free
in GiNaC (with ``normal()`` cancellation), and converts the ``N(s)/D(s)`` result
back to sympy -- re-mapping symbol assumptions (e.g. ``Symbol('s',
imaginary=True)``) that do not survive the string round-trip.

Importing this module requires the compiled ``_ginac_ext`` extension; callers
guard for :class:`ImportError` and fall back to the pure-sympy path.
"""
import numpy as np
import sympy

from pycircuit.circuit import _ginac_ext


def _to_ginac(expr):
    # Rationalize floats first: GiNaC's normal() cancels cleanly over exact
    # rationals but not over floats (float GCD leaves the result uncancelled and
    # the solve slow), mirroring symbolic_poly's get_exact() discipline.  Then
    # translate sympy's ``**`` power to GiNaC's ``^``.
    expr = sympy.sympify(expr)
    if expr.has(sympy.Float):
        # Use the simplest rational (small denominator) rather than the exact
        # dyadic value: Rational(Float(1e-9)) has a ~90-bit denominator, and such
        # denominators multiply into astronomically large integers in the
        # determinant that stall GiNaC.  limit_denominator gives 1e-9 -> 1/1e9.
        expr = expr.replace(
            lambda e: e.is_Float,
            lambda e: sympy.Rational(e).limit_denominator(10**15))
    return str(expr).replace('**', '^')


def _to_sympy(s, orig_syms):
    e = sympy.sympify(s.replace('^', '**'))
    subs = {t: orig_syms[t.name] for t in e.free_symbols if t.name in orig_syms}
    return e.xreplace(subs) if subs else e


def linearsolver_num_den(A, b):
    """Solve ``A x = b`` in GiNaC; return ``(numerators, denominator)``.

    ``x_i = numerators[i] / denominator`` with the shared denominator equal to
    the network determinant, matching the :class:`SymbolicPolyToolkit` contract.
    """
    A = sympy.Matrix(A)
    b = sympy.Matrix(b.tolist())
    n = A.rows

    if (A.rows, A.cols) == (1, 1):
        return np.array([b[0]], dtype=object), A[0, 0]

    orig_syms = {sym.name: sym for sym in (A.free_symbols | b.free_symbols)}
    symbols = list(orig_syms.keys())
    entries = [_to_ginac(A[i, j]) for i in range(n) for j in range(n)]
    rhs = [_to_ginac(b[i]) for i in range(n)]

    num_strs, den_str = _ginac_ext.solve_numden(n, entries, rhs, symbols)

    num = np.array([_to_sympy(x, orig_syms) for x in num_strs], dtype=object)
    den = _to_sympy(den_str, orig_syms)
    return num, den
