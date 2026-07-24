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


#: Above this GiNaC-source length (characters) skip the native ``compile_ex``
#: path -- which invokes a C++ compiler at runtime -- and evaluate via sympy
#: instead.  A runaway ``compile_ex`` on a huge expression can make g++ exhaust
#: memory (it once OOM-crashed the box); keep the emitted source bounded.
MAX_COMPILE_CHARS = 200_000


def _to_ginac_eval(expr):
    """String form for the ``compile_ex`` (numeric) path.

    Unlike :func:`_to_ginac` this keeps floats as-is: ``compile_ex`` evaluates
    numerically, so rationalizing floats only bloats the source.  Only the power
    operator is translated (``**`` -> ``^``).
    """
    return str(sympy.sympify(expr)).replace('**', '^')


def eval_sweep(expr, params):
    """Evaluate ``expr`` over a parameter sweep, fast, via native GiNaC code.

    ``params`` maps each free symbol of ``expr`` (the symbol object or its name)
    to an array of values; the arrays broadcast to a common length ``N``.
    ``expr`` may be complex-valued (contain sympy ``I``, e.g. a transfer function
    evaluated at ``s = I*w``).  Returns a complex :class:`numpy.ndarray` of
    length ``N``.

    The expression is compiled to native machine code once (GiNaC
    ``compile_ex``) and evaluated per point -- the "derive once, evaluate many"
    win.  If the emitted source would exceed :data:`MAX_COMPILE_CHARS`, or the
    native path fails for any reason, it falls back to :func:`sympy.lambdify`
    (pure numpy, no compiler) so a runaway compile can never occur.
    """
    expr = sympy.sympify(expr)
    free = {sym.name: sym for sym in expr.free_symbols}

    items = [(str(k), v) for k, v in params.items()]
    names = [nm for nm, _ in items]
    missing = set(free) - set(names)
    if missing:
        raise ValueError(
            "params missing free symbols: %s" % sorted(missing))

    arrs = np.broadcast_arrays(
        *[np.atleast_1d(np.asarray(v, dtype=float)) for _, v in items])
    npts = arrs[0].shape[0] if arrs else 0

    gstr = _to_ginac_eval(expr)
    if len(gstr) <= MAX_COMPILE_CHARS:
        try:
            points = [[float(a[i]) for a in arrs] for i in range(npts)]
            res = _ginac_ext.eval_sweep(gstr, names, points)
            return np.asarray(res, dtype=complex)
        except Exception:
            pass  # fall through to the sympy path

    # Fallback: sympy lambdify -- uses the expression's real free symbols (so
    # assumptions like Symbol('s', imaginary=True) still substitute) and never
    # shells out to a compiler.
    syms = [free.get(nm, sympy.Symbol(nm)) for nm in names]
    f = sympy.lambdify(syms, expr, modules='numpy')
    out = f(*arrs)
    return np.asarray(out, dtype=complex) * np.ones(npts)


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
