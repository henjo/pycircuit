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


def _system_strings(A, b):
    """Common set-up: sympy ``(A, b)`` -> GiNaC-parseable strings + symbol map."""
    A = sympy.Matrix(A)
    b = sympy.Matrix(b.tolist())
    n = A.rows
    orig_syms = {sym.name: sym for sym in (A.free_symbols | b.free_symbols)}
    entries = [_to_ginac(A[i, j]) for i in range(n) for j in range(n)]
    rhs = [_to_ginac(b[i]) for i in range(n)]
    return n, entries, rhs, list(orig_syms.keys()), orig_syms


def linearsolver_num_den(A, b):
    """Solve ``A x = b`` in GiNaC; return ``(numerators, denominator)``.

    ``x_i = numerators[i] / denominator`` with the shared denominator equal to
    the network determinant, matching the :class:`SymbolicPolyToolkit` contract.

    This eagerly sympifies the full (determinant-sized) result.  When you only
    need a compact piece -- a transfer function, or a numeric sweep -- prefer
    :func:`solve_native`, which keeps ``num[]``/``den`` as native GiNaC ``ex``
    and converts only the small result you ask for.
    """
    Am = sympy.Matrix(A)
    if Am.shape == (1, 1):
        return np.array([sympy.Matrix(b.tolist())[0]], dtype=object), Am[0, 0]

    n, entries, rhs, symbols, orig_syms = _system_strings(A, b)
    num_strs, den_str = _ginac_ext.solve_numden(n, entries, rhs, symbols)

    num = np.array([_to_sympy(x, orig_syms) for x in num_strs], dtype=object)
    den = _to_sympy(den_str, orig_syms)
    return num, den


class GinacResult:
    """Native GiNaC solve result -- lazy, partial conversion back to sympy.

    Wraps the C++ :class:`_ginac_ext.GinacResult` (which holds the numerators
    and shared denominator as GiNaC ``ex``).  The determinant-sized full result
    is *never* round-tripped through sympy unless explicitly requested; the
    common outputs are compact and are cancelled natively in GiNaC before the
    single small piece crosses back:

    * :meth:`tf` -- a transfer function ``x_i/x_j = num_i/num_j`` (the shared
      denominator cancels), and
    * :meth:`denominator` -- the network determinant, for poles.

    Symbol assumptions that do not survive the string round-trip (e.g.
    ``Symbol('s', imaginary=True)``) are re-mapped onto the returned sympy
    expressions, as in :func:`linearsolver_num_den`.
    """

    def __init__(self, native, orig_syms):
        self._native = native
        self._orig_syms = orig_syms

    def __len__(self):
        return len(self._native)

    def numerator(self, i):
        """Raw numerator ``num[i]`` as sympy (``x_i = num[i]/den``); may be large."""
        return _to_sympy(self._native.num_str(i), self._orig_syms)

    def denominator(self):
        """Shared denominator (network determinant) as sympy -- for poles."""
        return _to_sympy(self._native.den_str(), self._orig_syms)

    def component(self, i):
        """Full solution component ``x_i = num[i]/den`` as sympy; may be large."""
        return _to_sympy(self._native.component_str(i), self._orig_syms)

    def tf(self, i, j):
        """Transfer function ``x_i/x_j`` as sympy -- cancelled natively, so compact."""
        return _to_sympy(self._native.tf_str(i, j), self._orig_syms)

    def tf_coeffs(self, i, j, var):
        """Transfer function ``x_i/x_j`` as ``(num_coeffs, den_coeffs)`` in ``var``.

        Each is a list of sympy coefficients, low order first, so that
        ``x_i/x_j = sum_k num_coeffs[k] var**k / sum_k den_coeffs[k] var**k`` --
        the canonical ``N(var)/D(var)`` form (e.g. ``var = s`` for poles/zeros or
        a ``scipy.signal`` transfer function).  Collected natively in GiNaC; only
        the coefficients are converted to sympy.
        """
        var = sympy.sympify(var)
        nums, dens = self._native.tf_coeffs(i, j, str(var))
        return ([_to_sympy(c, self._orig_syms) for c in nums],
                [_to_sympy(c, self._orig_syms) for c in dens])

    def denominator_coeffs(self, var):
        """Denominator (network determinant) as coefficients in ``var``, low first.

        The characteristic polynomial whose roots are the poles.
        """
        var = sympy.sympify(var)
        return [_to_sympy(c, self._orig_syms)
                for c in self._native.denominator_coeffs(str(var))]

    def eval_tf(self, i, j, params):
        """Numerically evaluate the transfer function ``x_i/x_j`` over a sweep.

        Composes the native-cancelled (compact) :meth:`tf` with
        :func:`eval_sweep`, so neither the large intermediate nor a large source
        ever materialises.  ``params`` is as in :func:`eval_sweep`.
        """
        return eval_sweep(self.tf(i, j), params)

    def eval_component(self, i, params):
        """Numerically evaluate the full component ``x_i = num[i]/den`` over a sweep."""
        return eval_sweep(self.component(i), params)


def solve_native(A, b):
    """Solve ``A x = b`` in GiNaC, returning a :class:`GinacResult`.

    Unlike :func:`linearsolver_num_den`, the determinant-sized numerators and
    denominator stay as native GiNaC ``ex``; only the compact pieces you ask for
    (a transfer function, the denominator, a numeric sweep) are converted back,
    avoiding the sympy round-trip that dominates a large symbolic solve.
    """
    n, entries, rhs, symbols, orig_syms = _system_strings(A, b)
    native = _ginac_ext.solve_native(n, entries, rhs, symbols)
    return GinacResult(native, orig_syms)


def _expr_strings(expr):
    """sympy ``expr`` -> (ginac string, symbol names, name->symbol map)."""
    expr = sympy.sympify(expr)
    orig_syms = {sym.name: sym for sym in expr.free_symbols}
    return _to_ginac(expr), list(orig_syms.keys()), orig_syms


def poly_coeffs(expr, var):
    """Polynomial-in-``var`` coefficients of a rational ``expr``.

    Returns ``(num_coeffs, den_coeffs)``, each a list of sympy coefficients low
    order first, so ``expr == sum_k num_coeffs[k] var**k / sum_k den_coeffs[k]
    var**k``.  The ``normal()`` + collection run in GiNaC (fast on large
    multivariate); only the coefficients are converted back to sympy.
    """
    var = sympy.sympify(var)
    gstr, symbols, orig_syms = _expr_strings(expr)
    if str(var) not in symbols:
        symbols.append(str(var))
    nums, dens = _ginac_ext.poly_coeffs(gstr, str(var), symbols)
    return ([_to_sympy(c, orig_syms) for c in nums],
            [_to_sympy(c, orig_syms) for c in dens])


def series(expr, var, order, point=0):
    """Truncated series of ``expr`` in ``var`` about ``point`` to ``order`` terms.

    Returns a sympy polynomial (the ``O(var**order)`` term dropped) -- a Laurent
    expansion for low-/high-frequency and pole/zero approximation.  Expanded in
    GiNaC (much faster than sympy on large multivariate); the compact result is
    converted back.
    """
    var = sympy.sympify(var)
    gstr, symbols, orig_syms = _expr_strings(expr)
    if str(var) not in symbols:
        symbols.append(str(var))
    s = _ginac_ext.series_expr(gstr, str(var), _to_ginac(point), int(order),
                               symbols)
    return _to_sympy(s, orig_syms)
