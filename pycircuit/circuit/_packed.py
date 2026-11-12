"""THE PACKED LAYOUT: a matrix stored on the circuit's structural pattern --
the places the stamp plan's G and C put values, and the diagonal -- from a
size up (`PACKED_FROM`), for the analyses that take it (`scope`): the DC
operating point and the multistep transient.  Below the size, and in every
other analysis, the matrices stay dense n x n numpy, as they were.

WHY.  A dense n x n float64 is 8 n^2 bytes of address space: the IBM power
grid ibmpg2t (164 898 unknowns) needs 217 GB for one, which the kernel will
not map, and its DC and transient made about 25 of them.  On the pattern a
matrix is its nonzeros' values (~1.1 M there).

THE MATRIX (`PatternMatrix`).  Its pattern (`Pattern`: n, the places in CSC
order -- column by column, rows ascending) and its values.  Every
operation the analyses make on the dense matrix is made per entry where its
order cannot matter -- the sum of two on one pattern, a scale, `abs`, a
reduction, a gather -- so each value is the dense one's bytes.  The
products the decisions read (`A @ x`, the convergence scale `|J| |x|`) are
summed in ONE order, the CSC kernel (`matvec`): the products rounded, then
added in the pattern's order, as `np.bincount` adds them -- the reference's
order at size, which the C twins repeat.  A dense product summed in BLAS's
order may round otherwise: a decision on its threshold may fall the other
way from the dense run's, nothing else moves.

LOUD, NOT SILENT.  A PatternMatrix is no ndarray: a numpy function or ufunc
it does not make raises `TypeError`, so a consumer not taught the layout
fails where it reads one -- never computes on a dense copy unseen.
`to_dense` makes one only below `DENSE_MAX` unknowns (a test's, a small
circuit's under the `PACKED_FROM` forcing); at size it raises.
"""
import threading

import numpy as np

#: the unknowns from which the DC and the multistep transient hold their
#: matrices on the pattern (`scope`); above every circuit of the gate's
#: size set (chains and meshes of 4 000 - 8 000) -- a gate forcing sets 0.
#: 65 536 while the C paths decline on the layout: ibmpg1t (54 266) keeps
#: its C path, and the grids that cannot run dense (ibmpg2t's 164 898 on)
#: take it; 32 768 once the C serves it (plan 261 L20 P2)
PACKED_FROM = 65536

#: the unknowns below which `to_dense` makes a dense copy (a fallback no
#: packed form exists for, on a small circuit); from here it raises
DENSE_MAX = 1 << 14

_SCOPE = threading.local()


class scope:
    """`with scope(cir, on):` -- inside, the stamp plan assembles `cir`'s
    G and C on the pattern (`active`) where `on` and the circuit has
    `PACKED_FROM` unknowns or more.  Nested scopes stack."""

    def __init__(self, cir, on=True):
        self.key = id(cir) if on and getattr(cir, 'n', 0) >= PACKED_FROM else None

    def __enter__(self):
        stack = getattr(_SCOPE, 'stack', None)
        if stack is None:
            stack = _SCOPE.stack = []
        stack.append(self.key)
        return self

    def __exit__(self, *exc):
        _SCOPE.stack.pop()
        return False


def active(cir):
    """Whether `cir`'s matrices are assembled on the pattern here."""
    stack = getattr(_SCOPE, 'stack', None)
    return bool(stack) and stack[-1] is not None and stack[-1] == id(cir)


def any_active():
    """Whether the innermost scope assembles on the pattern: the C paths
    that mirror the dense layout decline inside it."""
    stack = getattr(_SCOPE, 'stack', None)
    return bool(stack) and stack[-1] is not None


class Pattern:
    """A structural pattern: `n` and its places in CSC order -- `cols`,
    `rows` (int64, nnz), `Ap` (n + 1), `key` (`col * n + row`, sorted) --
    with its reductions (`reduced`) and maps onto other patterns
    (`index_of`), each made once."""

    def __init__(self, n, rows, cols):
        n = int(n)
        rows = np.asarray(rows, dtype=np.int64)
        cols = np.asarray(cols, dtype=np.int64)
        key = cols * n + rows
        order = np.argsort(key, kind='stable')
        key = key[order]
        if key.size and np.any(np.diff(key) == 0):
            raise ValueError('a pattern place given twice')
        self.n, self.key = n, key
        self.cols, self.rows = key // n, key % n
        self.nnz = int(key.size)
        self.Ap = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(np.bincount(self.cols, minlength=n), out=self.Ap[1:])
        self._red = {}
        self._rowmajor = None

    def index_of(self, rows, cols):
        """Each place's index in this pattern, -1 where it is not one."""
        k = np.asarray(cols, dtype=np.int64) * self.n + np.asarray(rows, dtype=np.int64)
        i = np.searchsorted(self.key, k)
        i = np.minimum(i, max(self.nnz - 1, 0))
        hit = self.key[i] == k if self.nnz else np.zeros(k.shape, dtype=bool)
        return np.where(hit, i, -1)

    def reduced(self, drop):
        """`(the pattern without row and column drop, the kept entries'
        indices here)` -- as `remove_row_col` cuts a dense matrix."""
        got = self._red.get(drop)
        if got is None:
            keep = np.flatnonzero((self.rows != drop) & (self.cols != drop))
            r, c = self.rows[keep], self.cols[keep]
            got = self._red[drop] = (Pattern(self.n - 1, r - (r > drop), c - (c > drop)), keep)
        return got

    def kept(self, keep):
        """`(the pattern on the rows and columns `keep` -- ascending --
        renumbered in its order, the kept entries' indices here)` -- as
        `A[np.ix_(keep, keep)]` gathers a dense matrix (the order kept: the
        entries stay in CSC order)."""
        keep = np.asarray(keep, dtype=np.int64)
        if keep.size > 1 and not (np.diff(keep) > 0).all():
            raise ValueError('kept rows and columns not ascending')
        pos = np.full(self.n, -1, dtype=np.int64)
        pos[keep] = np.arange(keep.size)
        r, c = pos[self.rows], pos[self.cols]
        sel = np.flatnonzero((r >= 0) & (c >= 0))
        return Pattern(keep.size, r[sel], c[sel]), sel

    def diagonal(self):
        """The diagonal's entries' indices, row by row (every one a place)."""
        got = self._red.get('diag')
        if got is None:
            got = self._red['diag'] = self.index_of(np.arange(self.n), np.arange(self.n))
            if (got < 0).any():
                raise ValueError('a pattern without its whole diagonal')
        return got

    def rowmajor(self):
        """The entries' order by (row, column): numpy's `nonzero` order."""
        if self._rowmajor is None:
            self._rowmajor = np.lexsort((self.cols, self.rows))
        return self._rowmajor


class PatternMatrix:
    """A square float64 matrix on a `Pattern`: zero off it (the module's
    note: THE MATRIX, LOUD, NOT SILENT)."""

    __slots__ = ('pat', 'v')
    __array_priority__ = 1000
    ndim = 2
    dtype = np.dtype(np.float64)

    def __init__(self, pat, v):
        self.pat = pat
        self.v = np.asarray(v, dtype=np.float64)

    @property
    def shape(self):
        return (self.pat.n, self.pat.n)

    @property
    def size(self):
        return self.pat.n * self.pat.n

    def __len__(self):
        return self.pat.n

    def __repr__(self):
        return f'PatternMatrix(n={self.pat.n}, nnz={self.pat.nnz})'

    def _same(self, o):
        if not (isinstance(o, PatternMatrix) and o.pat is self.pat):
            raise TypeError('PatternMatrix: an operand on another pattern, or not one')
        return o.v

    @staticmethod
    def _scalar(s):
        if isinstance(s, (int, float, np.integer, np.floating)) or (
                isinstance(s, np.ndarray) and s.ndim == 0):
            return s
        raise TypeError(f'PatternMatrix: no elementwise operation with {type(s).__name__}')

    def __add__(self, o):
        return PatternMatrix(self.pat, self.v + self._same(o))

    __radd__ = __add__

    def __sub__(self, o):
        return PatternMatrix(self.pat, self.v - self._same(o))

    def __rsub__(self, o):
        return PatternMatrix(self.pat, self._same(o) - self.v)

    def __neg__(self):
        return PatternMatrix(self.pat, -self.v)

    def __mul__(self, s):
        return PatternMatrix(self.pat, self.v * self._scalar(s))

    def __rmul__(self, s):
        return PatternMatrix(self.pat, self._scalar(s) * self.v)

    def __truediv__(self, s):
        return PatternMatrix(self.pat, self.v / self._scalar(s))

    def __abs__(self):
        return PatternMatrix(self.pat, np.abs(self.v))

    def copy(self):
        return PatternMatrix(self.pat, self.v.copy())

    def astype(self, dtype, copy=True):
        if np.dtype(dtype) != np.float64:
            raise TypeError('PatternMatrix: float64 only')
        return self.copy() if copy else self

    def reduced(self, drop):
        """Without row and column `drop`: the same bytes as the dense cut."""
        pat, keep = self.pat.reduced(drop)
        return PatternMatrix(pat, self.v[keep])

    def kept(self, keep, got=None):
        """`A[np.ix_(keep, keep)]` (`keep` ascending): the kept rows and
        columns gathered, the same bytes -- `got` the pattern's `kept(keep)`
        where a caller holds it."""
        pat, sel = got if got is not None else self.pat.kept(keep)
        return PatternMatrix(pat, self.v[sel])

    def plus_diagonal(self, d):
        """`J + eye(n) * d` for `d` finite and not negative (a scalar or a
        vector of n): every entry `+ 0.0` (a -0.0 made +0.0, as the dense
        sum makes it), `d` added on the diagonal -- `nrsolver._plus_diagonal`'s
        in-place bytes.  Any other `d` raises."""
        if isinstance(d, np.ndarray) and d.ndim == 1:
            ok = d.shape == (self.pat.n,) and bool(np.isfinite(d).all()) and not np.signbit(d).any()
        else:
            d = float(d)
            ok = np.isfinite(d) and not np.signbit(d)
        if not ok:
            raise TypeError('PatternMatrix: plus_diagonal of a negative or non-finite d')
        di = self.pat.diagonal()
        v = self.v + 0.0
        v[di] += d
        return PatternMatrix(self.pat, v)

    def matvec(self, x):
        """`A @ x` in the CSC order (THE MATRIX): `x` 1-D, n."""
        x = np.asarray(x)
        if x.ndim != 1 or x.shape[0] != self.pat.n:
            raise TypeError('PatternMatrix: matvec takes a vector of n')
        return np.bincount(self.pat.rows, weights=self.v * x[self.pat.cols], minlength=self.pat.n)

    def __matmul__(self, x):
        return self.matvec(x)

    def at(self, rows, cols):
        """The values at the places `(rows, cols)` -- 0.0 off the pattern."""
        i = self.pat.index_of(rows, cols)
        return np.where(i >= 0, self.v[np.maximum(i, 0)], 0.0)

    def __getitem__(self, k):
        if isinstance(k, tuple) and len(k) == 2:
            r, c = (np.asarray(e) for e in k)
            if r.dtype.kind in 'iu' and c.dtype.kind in 'iu':
                r, c = np.broadcast_arrays(r, c)
                return self.at(r.ravel(), c.ravel()).reshape(r.shape)
        raise TypeError('PatternMatrix: indexing by two integer arrays alone')

    def nonzero(self):
        """`np.nonzero`'s positions (row-major order) of the nonzero values."""
        o = self.pat.rowmajor()
        on = self.v[o] != 0.0
        return self.pat.rows[o][on], self.pat.cols[o][on]

    def values_rowmajor(self):
        """`(rows, cols, values)` of the nonzeros in `np.nonzero`'s order."""
        o = self.pat.rowmajor()
        v = self.v[o]
        on = v != 0.0
        return self.pat.rows[o][on], self.pat.cols[o][on], v[on]

    def max(self):
        """The largest entry, zeros off the pattern counted (a NaN propagates)."""
        m = float(np.max(self.v)) if self.v.size else 0.0
        return m if self.pat.nnz == self.size or np.isnan(m) or m >= 0.0 else 0.0

    def diagonal(self):
        """The diagonal's values, row by row."""
        return self.v[self.pat.diagonal()]

    def count_nonzero(self):
        return int(np.count_nonzero(self.v))

    def max_abs(self):
        return float(np.max(np.abs(self.v))) if self.v.size else 0.0

    def to_dense(self):
        """A dense copy, below `DENSE_MAX` unknowns alone."""
        n = self.pat.n
        if n >= DENSE_MAX:
            raise MemoryError(f'PatternMatrix: no dense copy of {n} x {n} (DENSE_MAX {DENSE_MAX})')
        out = np.zeros((n, n))
        out[self.pat.rows, self.pat.cols] = self.v
        return out

    def __array__(self, *a, **k):
        raise TypeError('PatternMatrix: not converted to a dense array (to_dense, below DENSE_MAX)')

    def __array_ufunc__(self, ufunc, method, *inputs, **kw):
        if method == '__call__' and ufunc is np.absolute and len(inputs) == 1 and not kw:
            return abs(self)
        if method == '__call__' and ufunc in (np.add, np.subtract, np.multiply, np.true_divide) and not kw:
            a, b = inputs
            op = {np.add: '__add__', np.subtract: '__sub__', np.multiply: '__mul__',
                  np.true_divide: '__truediv__'}[ufunc]
            if a is self:
                return getattr(self, op)(b)
            rop = {'__add__': '__radd__', '__sub__': '__rsub__', '__mul__': '__rmul__'}.get(op)
            if rop is not None:
                return getattr(self, rop)(a)
        raise TypeError(f'PatternMatrix: numpy {ufunc.__name__}.{method} not made on the pattern')

    def __array_function__(self, func, types, args, kwargs):
        if func is np.dot and len(args) == 2 and args[0] is self and not kwargs:
            return self.matvec(args[1])
        if func is np.count_nonzero and len(args) == 1 and not kwargs:
            return self.count_nonzero()
        if func is np.nonzero and len(args) == 1 and not kwargs:
            return self.nonzero()
        if func is np.copy and len(args) == 1:
            return self.copy()
        raise TypeError(f'PatternMatrix: numpy {func.__name__} not made on the pattern')


def own_solver(A):
    """A `linearsolver.StructuralKLU` on `A`'s own pattern, kept with the
    pattern: where a caller gives no solver and the dense default
    (`toolkit.linearsolver`, numpy's solve) has no matrix to take."""
    pat = A.pat
    sk = pat._red.get('klu')
    if sk is None:
        from pycircuit.circuit.linearsolver import StructuralKLU
        sk = pat._red['klu'] = StructuralKLU(pat.n, pat.rows, pat.cols, fallback=None)
    return sk


def choice(ls, m):
    """The solver THE PACKED LAYOUT's m x m matrices take on `ls` (an
    AutoSolver) where a C path binds one: KLU on the circuit's structure
    whatever the size (`AutoSolver.choose_packed`, `_select`'s own choice
    for a PatternMatrix), analysed as its first factor analyses it -- or
    None (no structure of that size: the C declines)."""
    ch = ls.choose_packed(m)
    if ch is None or ch is ls._dense:
        return None
    try:
        ch._analysed()
    except np.linalg.LinAlgError:
        return None
    return ch


def is_packed(a):
    return isinstance(a, PatternMatrix)

