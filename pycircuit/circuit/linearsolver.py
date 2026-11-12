"""Linear-solver strategies for the Newton iteration, in the shape the package
already uses for ``Integrator``, ``StepController``, ``Scaler`` and
``NonLinearSolver``.

STAGE 7b.  Every Newton iteration solves ``J dx = -F`` through
``toolkit.linearsolver``, which is ``numpy.linalg.solve`` on the numeric backend:
a dense LU of the whole matrix, every time.  Circuit Jacobians are extremely
sparse -- an RC ladder at n=1101 is **0.27% non-zero** -- and a sparse
factorisation of that is **74x** faster.

THE SELECTION KEYS ON FILL, NOT ON ``n``, and that is a departure from the plan.
7b specifies "DenseLU below ~200 unknowns, SuperLU above".  Measured, the
crossover follows the sparsity *pattern*: on the RC ladder SuperLU already wins
2.5x at n=31, while a dense-ish circuit at n=300 would be mis-served by a
size-only rule.  ``n`` is used only as a cheap floor below which the measurement
itself is not worth taking.

THE DEFAULT IS ``AutoSolver`` SINCE 2026-10-02 (Andreas: "Yes, set AutoSolver
default"); until then it was the old dense path, deliberately.
``numpy.linalg.solve`` and SuperLU round differently, so a circuit that switches
gets different last bits.  ``AutoSolver`` leaves anything below its thresholds
(`MIN_N_FOR_SPARSE`, 250 unknowns, the measured end-to-end crossover) on exactly
the call the code made before, and reaches for the sparse path only where the
win is worth a changed rounding: 1.17x a transient at 402 unknowns, 1.60x at
1202, the answers within 3e-17.

A further measured caution, which is why ``DenseSolver`` does NOT use
``scipy.linalg.lu_factor``: at n=31 ``numpy.linalg.solve`` takes 2e-5 s against
3.2e-4 s for ``lu_factor`` + ``lu_solve`` -- an order of magnitude *worse*, on
Python call overhead alone.  A factor/solve split only pays when the factors are
reused, and on this code path they are not: see the note on ``factor()`` below.
"""
import functools
import importlib.machinery
import importlib.util
import os
import sys
from abc import ABC, abstractmethod

import numpy

from pycircuit.circuit import _lazy

_GETRS = {}


def lu_solve(lu_and_piv, b, trans=0):
    """`scipy.linalg.lu_solve(lu_and_piv, b, trans=trans)`, the same LAPACK
    `getrs` with the same checks -- a non-finite `b` refused with the same
    message, the shapes checked, a non-zero `info` raised -- without the
    array-API wrapper around it, which was 16 of pnoise's 34 profiled
    seconds (2026-09-30: a replay solves once per step, so a small `m`
    pays the wrapper millions of times).  Bit-identical results; anything
    it does not handle (an empty or batched `b`) goes to SciPy's."""
    lu, piv = lu_and_piv
    b1 = numpy.asarray(b)
    if (b1.dtype.char in numpy.typecodes['AllFloat']
            and not numpy.isfinite(b1).all()):
        raise ValueError('array must not contain infs or NaNs')
    if b1.size == 0 or b1.ndim > 2 or lu.ndim != 2:
        from scipy.linalg import lu_solve as _sla
        return _sla(lu_and_piv, b, trans=trans)
    if lu.shape[0] != b1.shape[0]:
        raise ValueError(f'Shapes of lu {lu.shape} and b {b1.shape} are '
                         'incompatible')
    key = (lu.dtype.char, b1.dtype.char)
    getrs = _GETRS.get(key)
    if getrs is None:
        from scipy.linalg import get_lapack_funcs
        getrs, = get_lapack_funcs(('getrs',), (lu, b1))
        _GETRS[key] = getrs
    x, info = getrs(lu, piv, b1, trans=trans, overwrite_b=False)
    if info == 0:
        return x
    raise ValueError(f'illegal value in {-info}th argument of internal '
                     'gesv|posv')


_NPLU = {}


def _numpy_lapack():
    """`(ctypes, dgesv, dgetrs)` of the OpenBLAS numpy itself loaded
    (`numpy.libs/libscipy_openblas64_*.so`, ILP64), or None."""
    if 'f' not in _NPLU:
        _NPLU['f'] = None
        try:
            import ctypes

            from pycircuit.circuit import _cbuild
            lib = _cbuild.openblas(numpy)
            if lib is not None:
                gesv, getrs = lib.scipy_dgesv_64_, lib.scipy_dgetrs_64_
                gesv.restype = getrs.restype = None
                _NPLU['f'] = (ctypes, gesv, getrs)
        except (OSError, AttributeError):
            pass
    return _NPLU['f']


def _scipy_sparsetools():
    """scipy.sparse's `_sparsetools` extension: scipy's own where scipy.sparse
    is imported, else the same file loaded alone under a private name -- its
    routines the very code, without scipy.sparse's package import (0.09 s of
    a run's start-up, for this one routine; the extension imports nothing)."""
    mod = sys.modules.get('scipy.sparse._sparsetools')
    if mod is not None:
        return mod
    spec = importlib.util.find_spec('scipy')
    if spec is None or not spec.origin:
        raise ImportError('no scipy')
    d = os.path.join(os.path.dirname(spec.origin), 'sparse')
    for suffix in importlib.machinery.EXTENSION_SUFFIXES:
        path = os.path.join(d, '_sparsetools' + suffix)
        if os.path.exists(path):
            name = 'pycircuit._scipy._sparsetools'
            loader = importlib.machinery.ExtensionFileLoader(name, path)
            mod = importlib.util.module_from_spec(
                importlib.util.spec_from_file_location(name, path, loader=loader))
            loader.exec_module(mod)
            return mod
    raise ImportError('no scipy.sparse._sparsetools')


def _csc_matvec():
    """SciPy's own `csc_matvec` (`scipy.sparse._sparsetools`: what a CSC
    matrix's `dot` with a vector runs), or None."""
    if 'mv' not in _NPLU:
        _NPLU['mv'] = None
        try:
            _NPLU['mv'] = _scipy_sparsetools().csc_matvec
        except (ImportError, OSError, AttributeError):
            pass
    return _NPLU['mv']


def _csc_dot(mv, n, Ap, Ai, data, x):
    """``csc_matrix((data, Ai, Ap)).dot(x)`` for a complex128 vector `x`:
    SciPy's `_matmul_vector` -- a zero vector, then `csc_matvec` into it."""
    y = numpy.zeros(n, dtype=numpy.complex128)
    mv(n, n, Ap, Ai, data, x, y)
    return y


def _csc_of_dense(A, mv, last=None):
    """`ComplexKLUSolver.prepare`'s record of a dense square complex128 `A`
    without SciPy's construction (382 k instructions a call on the PSP
    stage): the arrays ``csc_matrix(A).astype(complex128)`` holds -- the
    nonzeros (a NaN is one, a signed zero is not) column by column, rows
    ascending, COPIED: no arithmetic -- and its `dot` as SciPy runs it
    (`_csc_dot`).  Returned with its pattern ``(nz, Ap, Ai, key)``: given
    the last one, `last`, and the nonzeros where they fell, its `Ap`, `Ai`
    and key are this record's too -- functions of the pattern alone (speed
    round 10, B3.5); the values are copied as ever."""
    n = A.shape[0]
    At = A.T
    nz = At != 0
    if last is not None and last[0].shape == nz.shape and numpy.array_equal(last[0], nz):
        nz, Ap, Ai, key = last
    else:
        Ap = numpy.zeros(n + 1, dtype=numpy.int32)
        Ap[1:] = numpy.cumsum(nz.sum(axis=1))
        Ai = numpy.ascontiguousarray(nz.nonzero()[1], dtype=numpy.int32)
        key = (n, Ap.tobytes(), Ai.tobytes())
    data = At[nz]
    return ((functools.partial(_csc_dot, mv, n, Ap, Ai, data), n, Ap, Ai,
             data.view(numpy.float64), key), (nz, Ap, Ai, key))


#: THE HELD VALUES' STATE (`StructuralKLU._held_st`, int64; csrc/pcx_held.h
#: reads and writes the same buffer): OK -- the held values are those the
#: run's numeric was last refactored on; REUSED -- the refactors `solve`
#: did not make; and a C call's checkpoint (OPEN, CK_OK, MOVED, CK_REUSED,
#: CK_NUM: the state, the count and the numeric at its open, and whether
#: the held values moved since), which the Python never opens
HELD_OK, HELD_REUSED, HELD_OPEN, HELD_CK_OK, HELD_MOVED, HELD_CK_REUSED, HELD_CK_NUM = range(7)
HELD_N = 7


def _packed_cls():
    """`_packed.PatternMatrix` (THE PACKED LAYOUT), imported where first asked."""
    from pycircuit.circuit import _packed
    return _packed.PatternMatrix


def _same_bits(a, b):
    """`a` and `b` (float64, the same length) the same bit for bit -- as C's
    memcmp: -0.0 is not 0.0, a NaN equals a NaN of its own bits."""
    return numpy.array_equal(a.view(numpy.int64), b.view(numpy.int64))


def _held_move(sk):
    """`sk`'s held state about to change (pcx_held.h's `pcx_held_move`): an
    open checkpoint keeps what it found."""
    st = sk._held_st
    if st[HELD_OPEN] and not st[HELD_MOVED]:
        if st[HELD_CK_OK]:
            nnz = sk.nnz
            sk._held_v[nnz:2 * nnz] = sk._held_v[:nnz]
        st[HELD_MOVED] = 1


def _held_lost(sk):
    """`sk`'s numeric made afresh or replaced (`pcx_held_lost`): no values
    held, and an open checkpoint's numeric, which the replacement frees, no
    longer restored."""
    _held_move(sk)
    st = sk._held_st
    st[HELD_OK] = 0
    if st[HELD_OPEN]:
        st[HELD_CK_OK] = 0


def held_open(sk):
    """A checkpoint of `sk`'s held values (pcx_held.h: THE CHECKPOINT)
    before a C call on its run numeric that the glue may hand back: 1 where
    this call owns it (none was open), else 0."""
    st = sk._held_st
    if st[HELD_OPEN]:
        return 0
    num = sk._numeric
    st[HELD_OPEN] = 1
    st[HELD_CK_OK] = int(bool(st[HELD_OK]) and bool(num))
    st[HELD_MOVED] = 0
    st[HELD_CK_REUSED] = st[HELD_REUSED]
    st[HELD_CK_NUM] = int(num or 0)
    return 1


def held_commit(sk, own):
    """The owner's checkpoint closed on a served outcome."""
    if own:
        sk._held_st[HELD_OPEN] = 0


def held_back(sk, own):
    """The owner's checkpoint closed on a hand-back (`pcx_held_back`): the
    state, the count and the numeric as it found them -- the numeric
    recorded at the open refactored on the backup again where the held
    values moved and were valid then (the same values on the same pivots:
    the same numeric), so the Python's repeat decides as a run without the
    C."""
    if not own:
        return
    st = sk._held_st
    st[HELD_REUSED] = st[HELD_CK_REUSED]
    if st[HELD_MOVED]:
        num, ok = int(st[HELD_CK_NUM]), 0
        ## (the numeric found must be the run's still: a replaced one was
        ## freed, and `_held_lost` cleared CK_OK for it)
        if st[HELD_CK_OK] and num and num == (sk._numeric or 0):
            ct, nnz = sk._ct, sk.nnz
            bak = sk._held_v[nnz:2 * nnz]
            ok = sk._lib.klu_refactor(
                sk._ptr(sk.Ap, ct.c_int), sk._ptr(sk.Ai, ct.c_int), sk._ptr(bak, ct.c_double),
                ct.c_void_p(sk._symbolic), ct.c_void_p(num), ct.byref(sk._common))
            if ok:
                sk._held_v[:nnz] = bak
        st[HELD_OK] = int(bool(ok))
    st[HELD_OPEN] = 0
    st[HELD_MOVED] = 0


def _on_pattern_pages(arr):
    """Whether `arr` is a matrix on pattern pages (`_tran_core.pattern_of`)."""
    from pycircuit.circuit import _tran_core
    return _tran_core.pattern_of(arr) is not None


def _resid_sparse(Ap, Ai, Ax, tol, b, x):
    """The refactor's residual test -- `not (max|A x - b| / max(max|b|,
    1e-300) <= tol)`, numpy's `A @ x` on the dense A -- decided on A's CSC
    (`Ap`, `Ai`, `Ax`; A zero off it) where its rounding bound decides it:
    True it fails, False it holds, None undecided (the dense test then).
    `pcx_verdict.h`'s `pcx_resid_sparse`: row i's dense sum lies within
    `e_i = 4 (cnt_i + 2) (eps t_i + 2 tiny)` of the sparse one (`t_i` the
    sum of |a| |x|, `cnt_i` the row's stored entries), so the test holds where
    it holds at the interval's upper end and fails where it fails at the
    lower; a sum past 1e300 is undecided.  Where it does not decide,
    `_resid_last`: a value not finite, or the wide pass (`pcx_verdict.h`'s
    WHERE NEITHER DECIDES)."""
    got = _resid_bound(Ap, Ai, Ax, tol, b, x)
    return got if got is not None else _resid_last(Ap, Ai, Ax, tol, b, x)


def _resid_last(Ap, Ai, Ax, tol, b, x):
    """`pcx_resid_last`: (N) where a value of x, b or A is not finite, else
    (W), the wide pass."""
    if numpy.isfinite(b).all() and numpy.isfinite(x).all() and numpy.isfinite(Ax).all():
        return _resid_wide(Ap, Ai, Ax, tol, b, x)
    return _resid_nonfinite(Ap, Ai, Ax, b, x)


def _snan(a):
    """Each NaN of `a` whose quiet bit is clear."""
    return numpy.isnan(a) & ((a.view(numpy.uint64) & numpy.uint64(1 << 51)) == 0)


def _resid_nonfinite(Ap, Ai, Ax, b, x):
    """`pcx_resid_nonfinite` (N): True (it fails) where no infinity is in A or
    x, every NaN is quiet and each row's sum of |a| |x| over its finite terms
    plus a finite |b| is under 1e300 -- the dense test then raises no flag
    and fails -- else None."""
    if (numpy.isinf(x).any() or numpy.isinf(Ax).any() or _snan(x).any() or _snan(b).any()
            or _snan(Ax).any()):
        return None
    m = len(b)
    col = numpy.repeat(numpy.arange(m), numpy.diff(Ap))
    with numpy.errstate(all='ignore'):
        xc = x[col]
        fin = ~numpy.isnan(Ax) & ~numpy.isnan(xc)
        t = numpy.bincount(Ai, weights=numpy.where(fin, numpy.abs(Ax) * numpy.abs(xc), 0.0),
                           minlength=m)
        ok = (t + numpy.where(numpy.isfinite(b), numpy.abs(b), 0.0) < 1e300).all()
    return True if ok else None


def _resid_wide(Ap, Ai, Ax, tol, b, x):
    """`pcx_resid_wide` (W) on finite values: the rows' sums in numpy's long
    double (where it holds 64 bits; else None), each row's dense value
    within `E = k u (1 + 2^-6) T + 2^-60 |b| + (2 k + 4) 2^-1074` of its
    sum -- the ends rounded to double, then the bound's test."""
    L = numpy.longdouble
    if numpy.finfo(L).nmant < 63:
        return None
    m = len(b)
    col = numpy.repeat(numpy.arange(m), numpy.diff(Ap))
    a, xc = Ax.astype(L), x[col].astype(L)
    T, S = numpy.zeros(m, dtype=L), numpy.zeros(m, dtype=L)
    numpy.add.at(T, Ai, numpy.abs(a) * numpy.abs(xc))
    k = numpy.bincount(Ai, minlength=m).astype(L)
    ab = numpy.abs(b)
    if not (T + ab.astype(L) < L(1e300)).all():
        return None
    E = (k * (L(0.5) * L(numpy.finfo(numpy.float64).eps)) * (L(1.0) + L(2.0) ** -6) * T
         + L(2.0) ** -60 * ab.astype(L) + (L(2.0) * k + L(4.0)) * L(2.0) ** -1074)
    E = E.astype(numpy.float64).astype(L)
    numpy.add.at(S, Ai, a * xc)
    bl = b.astype(L)
    u = ((S + E) - bl).astype(numpy.float64)
    v = ((S - E) - bl).astype(numpy.float64)
    au, av = numpy.abs(u), numpy.abs(v)
    h = numpy.maximum(au, av)
    same = ((u > 0.0) & (v > 0.0)) | ((u < 0.0) & (v < 0.0))
    lo_row = numpy.where(same, numpy.minimum(au, av), 0.0)
    hi = float(h.max()) if m else 0.0
    lo = float(lo_row.max()) if m else 0.0
    sc = max(float(ab.max()) if m else 0.0, 1e-300)
    if hi / sc <= tol:
        return False
    if not (lo / sc <= tol):
        return True
    return None


def _resid_bound(Ap, Ai, Ax, tol, b, x):
    """`_resid_sparse`'s first pass: the rounding bound on the double sums
    (`pcx_resid_gradual`), None where it does not decide."""
    m = len(b)
    if not (numpy.isfinite(b).all() and numpy.isfinite(x).all()):
        return None
    col = numpy.repeat(numpy.arange(m), numpy.diff(Ap))
    with numpy.errstate(all='ignore'):
        xc = x[col]
        y = numpy.bincount(Ai, weights=Ax * xc, minlength=m)
        t = numpy.bincount(Ai, weights=numpy.abs(Ax) * numpy.abs(xc), minlength=m)
        cnt = numpy.bincount(Ai, minlength=m).astype(numpy.float64)
        ab = numpy.abs(b)
        if not (t + ab < 1e300).all():
            return None
        e = 4.0 * (cnt + 2.0) * (numpy.finfo(numpy.float64).eps * t
                                 + 2.0 * numpy.finfo(numpy.float64).smallest_subnormal)
        u = (y + e) - b
        v = (y - e) - b
        au, av = numpy.abs(u), numpy.abs(v)
        h = numpy.maximum(au, av)
        same = ((u > 0.0) & (v > 0.0)) | ((u < 0.0) & (v < 0.0))
        lo_row = numpy.where(same, numpy.minimum(au, av), 0.0)
    hi = float(h.max()) if m else 0.0
    lo = float(lo_row.max()) if m else 0.0
    sc = max(float(ab.max()) if m else 0.0, 1e-300)
    if hi / sc <= tol:
        return False
    if not (lo / sc <= tol):
        return True
    return None


def _count_nonzero(arr):
    """`numpy.count_nonzero(arr)` -- at its pattern where `arr` is a matrix
    on pattern pages (`_tran_core.pattern_count`), where numpy's would read
    it whole."""
    from pycircuit.circuit import _tran_core
    return _tran_core.pattern_count(arr)


def _csc_of_packed(data, m, rows, cols, mv, last=None):
    """`_csc_of_dense(A, mv, last)` for the m x m `A` zero off the positions
    `rows`, `cols` (column by column, rows ascending), `data` its entries
    there: the same record -- the nonzeros (`!= 0`: a NaN is one, a signed
    zero is not) in that order, COPIED, their `Ap`, `Ai`, key and mask, the
    last pattern's arrays where the mask is its."""
    keep = data != 0
    r, c = rows[keep], cols[keep]
    Ai = numpy.ascontiguousarray(r, dtype=numpy.int32)
    Ap = numpy.zeros(m + 1, dtype=numpy.int32)
    numpy.cumsum(numpy.bincount(c, minlength=m), out=Ap[1:])
    key = (m, Ap.tobytes(), Ai.tobytes())
    if last is not None and last[0].shape == (m, m) and last[3] == key:
        nz, Ap, Ai, key = last
    else:
        nz = numpy.zeros((m, m), dtype=bool)
        nz[c, r] = True
    vals = data[keep]
    return ((functools.partial(_csc_dot, mv, m, Ap, Ai, vals), m, Ap, Ai,
             vals.view(numpy.float64), key), (nz, Ap, Ai, key))


def _csc_of_values(data, mv, last):
    """`_csc_of_dense`'s record of a matrix whose nonzeros fell in the
    last one's pattern, `last`, given those nonzeros in its order -- `data`,
    complex128, made elsewhere (the transform's frozen factors in C, speed
    round 12, stage 5b): the pattern's `Ap`, `Ai` and key, `data` as the
    values, its `dot` -- what `_csc_of_dense` returns where its pattern is
    the last one."""
    nz, Ap, Ai, key = last
    n = nz.shape[0]
    return ((functools.partial(_csc_dot, mv, n, Ap, Ai, data), n, Ap, Ai,
             data.view(numpy.float64), key), (nz, Ap, Ai, key))


class NumpyLU:
    """THE FACTORISATION `numpy.linalg.solve` MAKES, KEPT (speed round 10,
    B3.1): `numpy.linalg.solve(A, b)` for every real 1-D `b`, bit for bit,
    with `A` factored once.  The first `solve` IS numpy's call -- its own
    OpenBLAS's `dgesv` on a Fortran copy of `A` and a copy of `b` -- and
    keeps the factors; each later one is that library's `dgetrs` against
    them.  ⚠ NOT `dgetrf` FIRST: `dgesv` with one right-hand side factors
    on one thread below 10000 unknowns, a `dgetrf` of 100 unknowns or more
    on several, and their bits differ (210 of 360 solves from 100 to 1000
    unknowns; `dgesv` then `dgetrs`: 2920 of 2920 the same from 1 to 1000,
    threaded).  Not SciPy's `lu_factor`: SciPy's OpenBLAS rounds
    differently.  A singular `A` (`info` > 0) raises numpy's
    `LinAlgError('Singular matrix')` at every `solve`, as numpy's solve
    would.  `make` returns None where numpy's library is not found.
    `owner`: the one reader of the factors, a weak reference its caller
    keeps (None: not set) -- what decides that a kept LU may be refilled."""

    __slots__ = ('_a', '_ct', '_gesv', '_getrs', '_info', '_p', 'addr', 'n', 'owner')

    @classmethod
    def make(cls, A, reuse=None):
        """The kept LU of `A`, unfactored.  `reuse`, a kept LU of the size
        whose factors nothing reads any more, is refilled with `A` in place
        and returned -- its buffers, its LAPACK arguments and its addresses
        (`addr`, the C's) kept (speed round 10, B3.5): the values `A` copied
        as a fresh one copies them."""
        if (type(reuse) is cls and type(A) is numpy.ndarray and A.dtype == numpy.float64
                and A.shape == (reuse.n, reuse.n)):
            numpy.copyto(reuse._a[0], A)
            reuse._info = None
            return reuse
        f = _numpy_lapack()
        if f is None:
            return None
        ctypes, gesv, getrs = f
        self = cls()
        a = numpy.array(A, dtype=numpy.float64, order='F', copy=True)
        n = self.n = a.shape[0]
        ipiv = numpy.zeros(n, dtype=numpy.int64)
        nn = ctypes.c_int64(n)
        ## (each `byref` and `data_as` holds its object: the buffers live
        ## as long as these arguments)
        self._p = (ctypes.byref(ctypes.c_char(b'N')), ctypes.byref(nn),
                   ctypes.byref(ctypes.c_int64(1)), a.ctypes.data_as(ctypes.c_void_p),
                   ipiv.ctypes.data_as(ctypes.c_void_p), ctypes.c_size_t(1))
        self._a, self._info = (a, ipiv), None
        self.addr, self.owner = (a.ctypes.data, ipiv.ctypes.data), None
        self._ct, self._gesv, self._getrs = ctypes, gesv, getrs
        return self

    def solve(self, b):
        """``numpy.linalg.solve(A, b)`` for a real 1-D `b` of the size."""
        ct = self._ct
        trans, nn, one, a, ipiv, tl = self._p
        x = numpy.array(b, dtype=numpy.float64, copy=True)
        info = ct.c_int64(0)
        if self._info is None:
            self._gesv(nn, one, a, nn, ipiv, x.ctypes.data_as(ct.c_void_p), nn, ct.byref(info))
            self._info = info.value
        elif self._info == 0:
            self._getrs(trans, nn, one, a, nn, ipiv, x.ctypes.data_as(ct.c_void_p), nn,
                        ct.byref(info), tl)
        if self._info > 0:
            raise numpy.linalg.LinAlgError('Singular matrix')
        return x


class _Factored(object):
    """A reusable factorisation: whatever the solver made, plus how to use it.

    Deliberately thin.  Every solver here already has a factorisation object
    of its own (`lu_factor`'s tuple, `splu`'s object, KLU's numeric handle);
    this only gives them ONE name the callers can hold, so a caller storing
    `N` of them does not branch on which solver produced them.
    """

    __slots__ = ('_f', '_solve', '_solve_t', 'sparse')

    def __init__(self, f, solve, solve_transposed=None, sparse=False):
        self._f, self._solve = f, solve
        self._solve_t = solve_transposed
        ## a SPARSE factorisation: a caller keeping the factored matrix's
        ## companions per step (the shooting's `C`) keeps them sparse too
        self.sparse = bool(sparse)

    def solve(self, b):
        return self._solve(self._f, b)

    def solve_transposed(self, b):
        """``A^T x = b`` against the SAME factorisation, or ``None``.

        ⚠ THIS IS WHAT MAKES A REVERSE REPLAY POSSIBLE WITHOUT A REVERSE
        INTEGRATOR.  A shooting monodromy is a product of per-step solves,
        so its transpose is that product replayed in reverse order with
        transposed solves -- and every factorisation here already supports
        one (LAPACK's ``trans``, SuperLU's ``trans='T'``).  Demir &
        Roychowdhury (TCAD 22(2) 188-196) describe reverse integration as
        "often unavailable even in existing time-domain simulators and may
        require significant changes to core simulation routines"; that is
        true of a forward-only DENSE implementation and not of one that
        already stores its factors.

        `None` when the solver cannot do it, so a caller can refuse rather
        than silently transpose the wrong thing.
        """
        if self._solve_t is None:
            return None
        return self._solve_t(self._f, b)


def _split_complex(solve, b):
    """``solve(b)`` for a COMPLEX `b` against a REAL factorisation: two real
    solves.  ⚠ A real sparse factorisation (SuperLU, KLU) handed a complex
    `b` would cast it to real and drop the imaginary part; the dense
    ``lu_solve`` promotes it instead, so the dense path does not come here
    (its bits stay what they were)."""
    b = numpy.asarray(b)
    if numpy.iscomplexobj(b):
        return solve(b.real) + 1j * solve(b.imag)
    return solve(b)


class LinearSolver(ABC):
    """Solve ``A x = b`` for the Newton iteration."""

    @abstractmethod
    def solve(self, A, b, toolkit):
        """Return ``x`` such that ``A x = b``."""

    def solve_kept(self, A, b, toolkit):
        """`solve` where the caller's last solve of this very `A` left the
        solver's numeric holding its refactor (`_FrozenTransform.kept`):
        `StructuralKLU`'s makes no refactor; a solver keeping no numeric
        solves as `solve`."""
        return self.solve(A, b, toolkit)

    def chord_factor(self, A, toolkit):
        """The chord Newton's factorisation of the step's matrix: `factor`'s
        here (`StructuralKLU` keeps one from step to step)."""
        return self.factor(A, toolkit)

    def factor(self, A, toolkit):
        """Optional: return a reusable factorisation, or ``None``.

        A returned object exposes ``.solve(b)`` for the SAME ``A``.

        ⚠ STILL UNUSED ON THE TRANSIENT PATH, and the reason is worth
        stating because 7b-2 proposed otherwise.  Newton's last
        factorisation is of ``J(x_k)``; the step controller needs ``J`` at
        the *converged* ``x_{k+1}``, which nothing else factors.  Reusing
        across that boundary substitutes a different matrix -- measured at
        median 5.5e-9 and max 2.2e-5 relative on a diode circuit -- so it is
        an approximation, not a saved duplicate, and it cannot be done under
        a bit-identical gate.

        ⚠ THE SHOOTING PATH IS DIFFERENT and does use it.  Matrix-free
        variational shooting replays each stored timestep once per Krylov
        iteration against the SAME ``Jf`` -- not a nearby one -- so the
        factorisation is a saved duplicate there rather than an
        approximation, and without it the replay would refactor ``k`` times
        over.  Returning ``None`` is still valid: the caller falls back to
        ``solve`` and pays exactly that.
        """
        return None

    def solve_transposed(self, A, B, toolkit):
        """``A^T X = B`` for the `A` the last `solve` was handed (`B` m x k):
        on that solve's factorisation where the solver keeps one (KLU's
        numeric, `KLUSolver`), else ``toolkit.linearsolver`` of ``A^T`` --
        numpy's solve, its dgesv on the column-major copy, which is A's
        row-major layout (the C Newtons make the same call).  The Newton's
        rounding floor reads rows of J^-1 through it
        (`nrsolver._within_rounding`)."""
        return toolkit.linearsolver(numpy.asarray(A).T, B)


class DenseSolver(LinearSolver):
    """The historical path: whatever ``toolkit.linearsolver`` does.

    On the numeric backend that is ``numpy.linalg.solve``.  Kept as the default so
    existing results do not move, and so symbolic and other toolkits keep working
    unchanged -- this class is the only one that is toolkit-agnostic.
    """

    def solve(self, A, b, toolkit):
        return toolkit.linearsolver(A, b)

    def factor(self, A, toolkit):
        """A LAPACK LU, reusable -- or ``None`` for a symbolic matrix.

        ``toolkit.linearsolver`` has no factor/solve split to borrow, so
        this reaches for LAPACK directly.  That is safe precisely where it
        is used and nowhere else: an object-dtype matrix (the symbolic
        toolkits) returns ``None`` and the caller falls back to ``solve``,
        which stays toolkit-agnostic as this class promises.
        """
        A = numpy.asarray(A)
        if A.dtype == object:
            return None
        try:
            from scipy.linalg import lu_factor
        except ImportError:                               # pragma: no cover
            return None
        return _Factored(lu_factor(A), lambda f, b: lu_solve(f, b),
                         lambda f, b: lu_solve(f, b, trans=1))

    def __repr__(self):
        return 'DenseSolver()'


class SuperLUSolver(LinearSolver):
    """SuperLU through ``scipy.sparse.linalg.splu``.

    Discovered rather than required, the way ``symengine`` and ``_ginac_ext``
    already are: if SciPy is missing the class raises on construction and
    ``AutoSolver`` never selects it.
    """

    def __init__(self):
        try:
            import scipy.sparse as _sp
            import scipy.sparse.linalg as _spla
        except ImportError as e:                          # pragma: no cover
            raise ImportError(
                'SuperLUSolver needs scipy.sparse; install scipy or use '
                'DenseSolver()') from e
        self._sp = _sp
        self._spla = _spla

    def solve(self, A, b, toolkit):
        A = numpy.asarray(A)
        if A.dtype == object:
            ## A symbolic matrix cannot go to SuperLU; fall back rather than fail,
            ## because the selection may have been made on a numeric sibling.
            return toolkit.linearsolver(A, b)
        Acsc = self._sp.csc_matrix(A)
        return self._spla.splu(Acsc).solve(numpy.asarray(b))

    def factor(self, A, toolkit):
        """An ``splu`` object, which already IS a reusable factorisation.

        ⚠ ``splu`` takes ONE right-hand side column at a time in the shape
        the shooting replay uses, which is exactly the shape matrix-free
        wants -- one vector per Krylov iteration.  A dense multi-column
        solve is what the DENSE propagation does, and that is the
        comparison, not this.
        """
        A = numpy.asarray(A)
        if A.dtype == object:
            return None
        return _Factored(
            self._spla.splu(self._sp.csc_matrix(A)),
            lambda f, b: _split_complex(lambda r: f.solve(r), b),
            lambda f, b: _split_complex(lambda r: f.solve(r, trans='T'), b),
            sparse=True)

    def __repr__(self):
        return 'SuperLUSolver()'


## ---------------------------------------------------------------------------
## STAGE 7c -- KLU, through ctypes, discovered rather than required.
##
## KLU is the factorisation production circuit simulators use (Xyce among them)
## and the reason is the analyze/factor/REFACTOR split.  A circuit's Jacobian
## keeps its sparsity PATTERN across every Newton iteration and every timestep --
## only the values change -- so the fill-reducing ordering, which is most of the
## cost, need be computed once per circuit rather than once per solve.
##
## Measured here against `scipy`'s `splu`, which recomputes the ordering every
## call (7c's premise, confirmed: `klu_analyze` 0.000441 s against `klu_factor`
## 0.000378 s at n=3200, so the ordering is roughly half of a full factorisation,
## and `klu_refactor` at 0.000099 s is 4x cheaper than factoring again):
##
##      n      splu   klu refactor+solve   speedup
##    150  0.000163            0.000046     3.51x
##    400  0.000463            0.000060     7.75x
##    800  0.000541            0.000061     8.85x
##   1600  0.000907            0.000097     9.31x
##   3200  0.001714            0.000172     9.97x
##
## No pip dependency is added: `libklu.so` ships with SuiteSparse, and this binds
## it the way `_ginac_ext` and `symengine` are already discovered.  If the library
## is absent the class raises on construction and `AutoSolver` never selects it.
## ---------------------------------------------------------------------------

class KLUSolver(LinearSolver):
    """SuiteSparse KLU with the symbolic analysis reused across solves.

    THE REFACTOR PATH IS VALIDATED, NOT TRUSTED.  ``klu_refactor`` reuses the
    pivot ORDER chosen when the matrix was first factored.  That is the whole
    point -- and it is also its one hazard: if the values move far enough, pivots
    that were sound become unstable and the factorisation silently degrades.
    Production tools read KLU's reciprocal pivot growth to decide; reading it here
    would mean hard-coding offsets into a ``klu_common`` struct this binding
    cannot verify, so instead **the residual is checked directly** and a full
    ``klu_factor`` is redone when it is not good enough.  One sparse mat-vec
    against a factorisation is cheap, and it cannot be wrong about the layout.

    That choice is deliberate given what stage 9 and 7d kept turning up: three
    separate defects this session survived because nothing ever checked the thing
    itself.

    KLU's objects are freed where they are replaced -- `solve`'s analysis and
    factor at a new pattern, its factor at a residual fallback -- and at the
    solver's end, with `factor`'s analyses (`analyses_freed`, `factors_freed`
    count the first two).  The C units that borrow a `StructuralKLU`'s handles
    hold the solver itself and read its numeric before each call.
    """

    ## Relative residual above which a refactor is rejected and the ordering
    ## recomputed.  Loose enough not to trip on ordinary round-off, tight enough
    ## that a genuinely bad reuse is caught: a well-conditioned solve lands near
    ## 1e-16 and a degraded one is orders above.
    REFACTOR_RESIDUAL_TOL = 1e-8

    def __init__(self, libname='libklu.so.2'):
        import ctypes
        try:
            lib = ctypes.CDLL(libname)
        except OSError as e:
            raise ImportError(
                'KLUSolver needs SuiteSparse KLU (%s); install libklu / '
                'libsuitesparse, or use SuperLUSolver()' % libname) from e
        ## (scipy.sparse imported on first use: no served run reaches it)
        if importlib.util.find_spec('scipy') is None:         # pragma: no cover
            raise ImportError('KLUSolver needs scipy.sparse')
        _sp = _lazy.deferred('scipy.sparse')

        self._ct = ctypes
        self._sp = _sp
        self._lib = lib
        c_int_p = ctypes.POINTER(ctypes.c_int)
        c_dbl_p = ctypes.POINTER(ctypes.c_double)
        for nm, argt, rest in (
                ('klu_defaults', [ctypes.c_void_p], ctypes.c_int),
                ('klu_analyze', [ctypes.c_int, c_int_p, c_int_p, ctypes.c_void_p],
                 ctypes.c_void_p),
                ('klu_factor', [c_int_p, c_int_p, c_dbl_p, ctypes.c_void_p,
                                ctypes.c_void_p], ctypes.c_void_p),
                ('klu_refactor', [c_int_p, c_int_p, c_dbl_p, ctypes.c_void_p,
                                  ctypes.c_void_p, ctypes.c_void_p], ctypes.c_int),
                ('klu_solve', [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int,
                               ctypes.c_int, c_dbl_p, ctypes.c_void_p],
                 ctypes.c_int),
                ('klu_tsolve', [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int,
                                ctypes.c_int, c_dbl_p, ctypes.c_void_p],
                 ctypes.c_int),
                ('klu_free_numeric', [ctypes.c_void_p, ctypes.c_void_p],
                 ctypes.c_int),
                ('klu_free_symbolic', [ctypes.c_void_p, ctypes.c_void_p],
                 ctypes.c_int)):
            f = getattr(lib, nm)
            f.argtypes = argt
            f.restype = rest

        ## `klu_common` is ~200 bytes; over-allocating a buffer is safe and avoids
        ## committing to a struct layout that varies between SuiteSparse releases.
        self._common = ctypes.create_string_buffer(4096)
        if not lib.klu_defaults(ctypes.byref(self._common)):
            raise ImportError('klu_defaults failed; the KLU library looks unusable')

        self._pattern = None      # (n, indptr bytes, indices bytes)
        self._symbolic = None
        self._numeric = None
        self.refactors = 0        # counters, so the reuse is observable
        self.factors = 0
        self.analyses = 0
        self.residual_fallbacks = 0
        ## `factor`'s own: one analysis per sparsity pattern, one OWNED
        ## numeric factorisation per call (kept apart from `solve`'s
        ## single refactored one, whose counters the tests read)
        self._factor_symbolics = {}
        self.factor_analyses = 0
        self.numerics = 0
        ## (the analyses and factors `solve` freed where it replaced them)
        self.analyses_freed = 0
        self.factors_freed = 0

    def _ptr(self, arr, ctype):
        return arr.ctypes.data_as(self._ct.POINTER(ctype))

    def _free_numeric(self):
        """The factor `solve` holds, freed (`klu_free_numeric`): None after."""
        num, self._numeric = self._numeric, None
        if num:
            ct = self._ct
            self._lib.klu_free_numeric(ct.byref(ct.c_void_p(num)), ct.byref(self._common))
            self.factors_freed += 1

    def _free_symbolic(self):
        """The analysis `solve` holds, freed (`klu_free_symbolic`): None after."""
        sym, self._symbolic = self._symbolic, None
        if sym:
            ct = self._ct
            self._lib.klu_free_symbolic(ct.byref(ct.c_void_p(sym)), ct.byref(self._common))
            self.analyses_freed += 1

    def __del__(self):
        """KLU's objects the solver still holds, freed: `solve`'s factor and
        analysis, `factor`'s analyses -- no factorisation made from one is
        left, each holding its solver."""
        try:
            if getattr(self, '_numeric', None):
                self._free_numeric()
            if getattr(self, '_symbolic', None):
                self._free_symbolic()
            syms = getattr(self, '_factor_symbolics', None) or {}
            ct = self._ct
            for sym in syms.values():
                self._lib.klu_free_symbolic(ct.byref(ct.c_void_p(sym)), ct.byref(self._common))
            syms.clear()
        except Exception:                                     # noqa: BLE001, S110  # pragma: no cover
            pass                                              # interpreter shutdown

    def solve(self, A, b, toolkit):
        import numpy
        ctypes = self._ct
        arr = numpy.asarray(A)
        if arr.dtype == object:
            ## Symbolic matrices cannot go to KLU; fall back rather than fail.
            return toolkit.linearsolver(A, b)

        if numpy.iscomplexobj(b):
            return _split_complex(lambda r: self.solve(A, r, toolkit), b)

        Acsc = self._sp.csc_matrix(arr)
        n = Acsc.shape[0]
        Ap = numpy.ascontiguousarray(Acsc.indptr, dtype=numpy.int32)
        Ai = numpy.ascontiguousarray(Acsc.indices, dtype=numpy.int32)
        Ax = numpy.ascontiguousarray(Acsc.data, dtype=numpy.float64)
        ## ⚠ SEVERAL RIGHT-HAND SIDES: `klu_solve` takes a column-major
        ## `n x nrhs` block.  Until 2026-09-25 `b` was RAVELLED, so a matrix
        ## right-hand side (gear's dense monodromy walk hands one, `m x 2m`)
        ## crashed the PSS solve.
        rhs = numpy.asarray(b, dtype=numpy.float64)
        nrhs = 1 if rhs.ndim < 2 else rhs.shape[1]
        rhs = (numpy.ascontiguousarray(rhs.ravel()) if rhs.ndim < 2
               else numpy.asfortranarray(rhs))

        key = (n, Ap.tobytes(), Ai.tobytes())
        reused = key == self._pattern and self._symbolic and self._numeric

        if not reused:
            ## The pattern changed (or this is the first solve): a new ordering is
            ## genuinely required, and this is the cost KLU exists to pay ONCE.
            ## (The last pattern's factor and analysis freed first.)
            self._free_numeric()
            self._free_symbolic()
            self._symbolic = self._lib.klu_analyze(
                n, self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                ctypes.byref(self._common))
            self.analyses += 1
            if not self._symbolic:
                raise numpy.linalg.LinAlgError('klu_analyze failed (singular structure?)')
            self._numeric = self._lib.klu_factor(
                self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._symbolic),
                ctypes.byref(self._common))
            self.factors += 1
            if not self._numeric:
                raise numpy.linalg.LinAlgError('Singular matrix')
            self._pattern = key
        else:
            ok = self._lib.klu_refactor(
                self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._symbolic),
                ctypes.c_void_p(self._numeric), ctypes.byref(self._common))
            self.refactors += 1
            if not ok:
                raise numpy.linalg.LinAlgError('Singular matrix')

        x = rhs.copy(order='F')
        if not self._lib.klu_solve(
                ctypes.c_void_p(self._symbolic), ctypes.c_void_p(self._numeric),
                n, nrhs, self._ptr(x, ctypes.c_double), ctypes.byref(self._common)):
            raise numpy.linalg.LinAlgError('klu_solve failed')

        if reused:
            ## The validation described in the class docstring.  Only on the reuse
            ## path: a fresh factorisation chose its own pivots and needs no check.
            scale = max(float(numpy.abs(rhs).max()), 1e-300)
            resid = float(numpy.abs(Acsc.dot(x) - rhs).max()) / scale
            if not (resid <= self.REFACTOR_RESIDUAL_TOL):
                self.residual_fallbacks += 1
                self._free_numeric()
                self._numeric = self._lib.klu_factor(
                    self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                    self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._symbolic),
                    ctypes.byref(self._common))
                self.factors += 1
                if not self._numeric:
                    raise numpy.linalg.LinAlgError('Singular matrix')
                x = rhs.copy(order='F')
                if not self._lib.klu_solve(
                        ctypes.c_void_p(self._symbolic), ctypes.c_void_p(self._numeric),
                        n, nrhs, self._ptr(x, ctypes.c_double), ctypes.byref(self._common)):
                    raise numpy.linalg.LinAlgError('klu_solve failed')
        return x

    def solve_transposed(self, A, B, toolkit):
        """``A^T X = B`` on the numeric the last `solve` made of `A`
        (`klu_tsolve`, `B` m x k column-major, as `klu_solve` takes it);
        none yet: the base's own solve."""
        if not self._numeric or not self._symbolic:
            return LinearSolver.solve_transposed(self, A, B, toolkit)
        ctypes = self._ct
        x = numpy.array(B, dtype=numpy.float64, order='F', copy=True)
        nrhs = 1 if x.ndim < 2 else x.shape[1]
        if not self._lib.klu_tsolve(
                ctypes.c_void_p(self._symbolic), ctypes.c_void_p(self._numeric), x.shape[0],
                nrhs, self._ptr(x, ctypes.c_double), ctypes.byref(self._common)):
            raise numpy.linalg.LinAlgError('klu_tsolve failed')
        return x

    def factor(self, A, toolkit):
        """A reusable KLU factorisation of `A` that OWNS its numeric object
        (freed with it), with ``A x = b`` and ``A^T x = b`` (`klu_tsolve`)
        for one or several right-hand sides -- the shooting replays store
        one per step and transpose them for every adjoint surface.  The
        symbolic analysis is shared per sparsity pattern.  ``None`` for a
        symbolic matrix.

        ⚠ Until 2026-09-25 KLU had no `factor`: the replays fell back to a
        solve per replayed step with NO transposed solve, so gear's `ppv`
        and `floquet_modes` raised `AttributeError` under this solver.
        """
        ctypes = self._ct
        arr = numpy.asarray(A)
        if arr.dtype == object:
            return None
        Acsc = self._sp.csc_matrix(arr)
        n = Acsc.shape[0]
        Ap = numpy.ascontiguousarray(Acsc.indptr, dtype=numpy.int32)
        Ai = numpy.ascontiguousarray(Acsc.indices, dtype=numpy.int32)
        Ax = numpy.ascontiguousarray(Acsc.data, dtype=numpy.float64)
        key = (n, Ap.tobytes(), Ai.tobytes())
        sym = self._factor_symbolics.get(key)
        if sym is None:
            sym = self._lib.klu_analyze(
                n, self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                ctypes.byref(self._common))
            if not sym:
                raise numpy.linalg.LinAlgError('klu_analyze failed (singular structure?)')
            self._factor_symbolics[key] = sym
            self.factor_analyses += 1
        num = self._lib.klu_factor(
            self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
            self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(sym),
            ctypes.byref(self._common))
        if not num:
            raise numpy.linalg.LinAlgError('Singular matrix')
        self.numerics += 1
        fac = _KLUNumeric(self, sym, num, n)
        return _Factored(fac, lambda f, b: f.solve(b),
                         lambda f, b: f.solve_transposed(b), sparse=True)

    def __repr__(self):
        return ('KLUSolver(analyses=%d, factors=%d, refactors=%d, fallbacks=%d)'
                % (self.analyses, self.factors, self.refactors,
                   self.residual_fallbacks))


class StructuralKLU(KLUSolver):
    """KLU ON THE CIRCUIT'S STRUCTURAL PATTERN (the C stepping loop plan's
    stage 2; Andreas's decision, 2026-10-10: every sparse run).

    `KLUSolver` takes its pattern from the VALUES: each call converts the
    dense matrix (`csc_matrix`, a scan of every entry) and compares the
    pattern's bytes, and an exact zero -- a gm at cutoff -- is a new pattern
    and a new ordering.  Here the pattern is the circuit's own
    (`analysis.jacobian_structure`: every position an element can stamp,
    and the diagonal), fixed for the run: ONE analysis, and each call takes
    the values at the pattern's positions (one `take`) for a refactor
    (`solve`, validated by the residual as `KLUSolver`'s), a factor of its
    own (`factor`, with the transposed solve), or the chord's (`chord_factor`:
    `solve`'s sequence on a numeric the chord keeps from step to step).

    A matrix this pattern does not cover -- another shape, another dtype,
    a nonzero off the pattern (a monodromy, a block system the same size)
    -- is not this circuit's Jacobian: it goes to `fallback()`, today's
    choice for its size, and leaves this factorisation alone (`foreign`
    counts them).  Measured (2026-10-10, gear, 15 % of the run): rich3
    (269 unknowns) 2.71 -> 2.02 s, arom (295) 2.64 -> 1.98 s against
    SuperLU; the answers' last bits move (the plan's one re-baseline).

    THE SAME VALUES REUSE THE NUMERIC.  `solve` makes no refactor where the
    values are, bit for bit, those `_numeric` was last refactored on
    (`_held_v`; `reused` counts them): KLU's refactor of the same values on
    the same symbolic and pivot order is the same numeric, so the solve and
    its validation are what the refactor's would be.  A linear circuit's
    matrix repeats whenever its step's coefficients do -- the IBM power
    grid's 10 ns transient refactored the same values in 1 777 of its
    Newton's 2 011 refactors.  Only a refactor's values are held: a fresh
    factor's numeric (the first, a fallback's) is not a refactor's, and
    clears them.  The C units that repeat `solve`'s sequence on `_numeric`
    keep the same buffers (csrc/pcx_held.h), and one that hands back
    restores them, and the numeric, as its checkpoint found them.
    """

    def __init__(self, m, rows, cols, fallback, libname='libklu.so.2', stamped=None):
        super().__init__(libname)
        ## (a callable `(rows, cols) -> per entry, whether an element stamps
        ## it, or None`: the analysis trials KLU's order on it -- `_klu_order`)
        self._stamped = stamped
        self.m = m = int(m)
        rows = numpy.asarray(rows, dtype=numpy.intp)
        cols = numpy.asarray(cols, dtype=numpy.intp)
        self.Ai = numpy.ascontiguousarray(rows, dtype=numpy.int32)
        self.Ap = numpy.zeros(m + 1, dtype=numpy.int32)
        numpy.cumsum(numpy.bincount(cols, minlength=m), out=self.Ap[1:])
        self.take = rows * m + cols
        self.nnz = len(rows)
        self._fallback, self._fb = fallback, None
        self.foreign = 0
        self.pivot_fallbacks = 0
        ## THE CHORD'S NUMERIC (`chord_factor`): kept from step to step, its
        ## own counters (the C Newton's sparse chord writes them back)
        self._chord = None
        self.chord_refactors = 0
        self.chord_factors = 0
        self.chord_pivot_fallbacks = 0
        self.chord_residual_fallbacks = 0
        ## THE HELD VALUES (the class's note): the values `_numeric` was last
        ## refactored on, then a C checkpoint's backup of them; their state
        ## (`HELD_*`).  Buffers the C units bind for the run: never replaced
        self._held_v = numpy.zeros(2 * max(self.nnz, 1))
        self._held_st = numpy.zeros(HELD_N, dtype=numpy.int64)

    @property
    def reused(self):
        """The refactors `solve` did not make: the values were those
        `_numeric` was last refactored on."""
        return int(self._held_st[HELD_REUSED])

    def _ours(self, A):
        """`(A as an array, its values in CSC order)` where `A` is this
        circuit's reduced Jacobian -- `m x m` float64 with no nonzero off
        the pattern -- else None.  A `_packed.PatternMatrix` is its own `A`,
        its values put at their places in KLU's (`_packed_values`)."""
        from pycircuit.circuit import _packed
        if type(A) is _packed.PatternMatrix:
            return self._packed_values(A)
        arr = numpy.asarray(A)
        if arr.shape != (self.m, self.m) or arr.dtype != numpy.float64:
            return None
        Ax = arr.ravel().take(self.take)
        if _count_nonzero(arr) != numpy.count_nonzero(Ax):
            return None
        return arr, Ax

    def _packed_values(self, A):
        """`(A, KLU's values)` for a `_packed.PatternMatrix` of this size, its
        entries at their places in KLU's CSC order (zero at KLU's others), or
        None where a nonzero lies off KLU's pattern -- the map made once a
        pattern."""
        if A.shape != (self.m, self.m):
            return None
        got = self.__dict__.get('_pk_map')
        if got is None or got[0] is not A.pat:
            if self.__dict__.get('_pk_key') is None:
                kcols = numpy.repeat(numpy.arange(self.m, dtype=numpy.int64), numpy.diff(self.Ap))
                self._pk_key = (kcols * self.m + self.Ai.astype(numpy.int64), kcols)
            key = self._pk_key[0]
            want = A.pat.key
            i = numpy.minimum(numpy.searchsorted(key, want), max(key.size - 1, 0))
            idx = numpy.where(key[i] == want, i, -1) if key.size else numpy.full(want.size, -1)
            got = self._pk_map = (A.pat, idx)
        idx = got[1]
        off = idx < 0
        if off.any() and numpy.count_nonzero(A.v[off]):
            return None
        Ax = numpy.zeros(self.nnz)
        Ax[idx[~off]] = A.v[~off]
        return A, Ax

    def _packed_residual(self, Ax, rhs, x):
        """The refactor's residual test on a `_packed.PatternMatrix` --
        `not (max|A x - b| / max(max|b|, 1e-300) <= tol)`, `A x` THE PACKED
        LAYOUT's CSC kernel on KLU's own arrays (`_packed`: the reference's
        order at size), a NaN failing it as numpy's max propagates one."""
        cols = self._pk_key[1]
        if rhs.ndim < 2:
            r = numpy.bincount(self.Ai, weights=Ax * x[cols], minlength=self.m) - rhs
        else:
            r = numpy.stack([numpy.bincount(self.Ai, weights=Ax * x[cols, j], minlength=self.m)
                             for j in range(rhs.shape[1])], axis=1) - rhs
        scale = max(float(numpy.abs(rhs).max()), 1e-300)
        return not (float(numpy.abs(r).max()) / scale <= self.REFACTOR_RESIDUAL_TOL)

    def _other(self):
        self.foreign += 1
        if self._fb is None:
            self._fb = self._fallback()
        return self._fb

    def _analysed(self):
        ctypes = self._ct
        if self._symbolic is None:
            sym = self._lib.klu_analyze(
                self.m, self._ptr(self.Ap, ctypes.c_int), self._ptr(self.Ai, ctypes.c_int),
                ctypes.byref(self._common))
            self.analyses += 1
            if not sym:
                raise numpy.linalg.LinAlgError('klu_analyze failed (singular structure?)')
            ## (KLU's order tried where the pattern's always-zero places could
            ## mislead it, another kept where they do: `_klu_order`)
            got = self._order_tried(sym)
            if got is not None:
                self._lib.klu_free_symbolic(ctypes.byref(ctypes.c_void_p(sym)),
                                            ctypes.byref(self._common))
                sym = got
            self._symbolic = sym
        return self._symbolic

    def _order_tried(self, sym):
        """The symbolic to keep in place of KLU's own `sym`, or None
        (`_klu_order.better`): asked where the places an element stamps are
        known (`stamped`, given by the analysis that made the solver)."""
        if self._stamped is None:
            return None
        rows = self.Ai.astype(numpy.intp)
        cols = numpy.repeat(numpy.arange(self.m), numpy.diff(self.Ap))
        mask = self._stamped(rows, cols)
        if mask is None:
            return None
        from pycircuit.circuit import _klu_order
        return _klu_order.better(self, sym, mask)

    def _factor(self, Ax):
        ctypes = self._ct
        num = self._lib.klu_factor(
            self._ptr(self.Ap, ctypes.c_int), self._ptr(self.Ai, ctypes.c_int),
            self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._analysed()),
            ctypes.byref(self._common))
        if not num:
            raise numpy.linalg.LinAlgError('Singular matrix')
        return num

    def _free(self, num):
        ctypes = self._ct
        self._lib.klu_free_numeric(ctypes.byref(ctypes.c_void_p(num)), ctypes.byref(self._common))

    def solve(self, A, b, toolkit):
        return self._solve(A, b, toolkit, False)

    def solve_kept(self, A, b, toolkit):
        """`solve` without its refactor: the caller's last solve of this very
        `A` left the numeric holding a refactor of its values
        (`_FrozenTransform.kept`: a step's frozen real factor, refactored at
        its first solve) -- the same values on the same pivots make the same
        factors, so the refactor is not repeated, nor counted.  The residual
        check and its fresh factor stand, as in `solve`."""
        return self._solve(A, b, toolkit, True)

    def _solve(self, A, b, toolkit, kept):
        got = self._ours(A)
        if got is None:
            return self._other().solve(A, b, toolkit)
        if numpy.iscomplexobj(b):
            return _split_complex(lambda r: self._solve(A, r, toolkit, kept), b)
        ctypes = self._ct
        arr, Ax = got
        rhs = numpy.asarray(b, dtype=numpy.float64)
        nrhs = 1 if rhs.ndim < 2 else rhs.shape[1]
        rhs = (numpy.ascontiguousarray(rhs.ravel()) if rhs.ndim < 2
               else numpy.asfortranarray(rhs))
        reused = self._numeric is not None
        st = self._held_st
        if reused and not kept and st[HELD_OK] and _same_bits(Ax, self._held_v[:self.nnz]):
            ## (THE SAME VALUES: the numeric is their refactor already)
            st[HELD_REUSED] += 1
        elif reused and not kept:
            ok = self._lib.klu_refactor(
                self._ptr(self.Ap, ctypes.c_int), self._ptr(self.Ai, ctypes.c_int),
                self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._analysed()),
                ctypes.c_void_p(self._numeric), ctypes.byref(self._common))
            self.refactors += 1
            _held_move(self)
            st[HELD_OK] = int(ok)
            if ok:
                self._held_v[:self.nnz] = Ax
            else:
                _held_lost(self)
            if not ok:
                ## ⚠ THE REUSED PIVOT ORDER MET A ZERO PIVOT: the matrix may
                ## well be regular under another, so a fresh factor decides --
                ## as a failed residual does below -- and only its failure is
                ## the matrix's.  Raised here until 2026-10-11, it failed
                ## opampal's operating point (519 unknowns): the gmin and
                ## source ladders move the matrix far between solves, and a
                ## pivot the first factor chose went to zero.
                self.pivot_fallbacks += 1
                old, self._numeric = self._numeric, None
                self._free(old)
                self._numeric = self._factor(Ax)
                self.factors += 1
                reused = False
        elif not reused:
            _held_lost(self)
            self._numeric = self._factor(Ax)
            self.factors += 1
        x = rhs.copy(order='F')
        if not self._lib.klu_solve(
                ctypes.c_void_p(self._symbolic), ctypes.c_void_p(self._numeric), self.m, nrhs,
                self._ptr(x, ctypes.c_double), ctypes.byref(self._common)):
            raise numpy.linalg.LinAlgError('klu_solve failed')
        if reused:
            ## (the refactor's validation, as `KLUSolver`'s: the reused pivot
            ## order checked by the residual, a fresh factor where it fails.
            ## A matrix on pattern pages: the verdict on KLU's CSC with its
            ## rounding bound where the bound decides -- `pcx_kluseq.h`'s
            ## rule, `_resid_sparse` -- the dense test where it does not)
            fails = None
            if type(arr) is _packed_cls():
                fails = self._packed_residual(Ax, rhs, x)
            elif nrhs == 1 and _on_pattern_pages(arr):
                fails = _resid_sparse(self.Ap, self.Ai, Ax, self.REFACTOR_RESIDUAL_TOL, rhs, x)
            if fails is None:
                scale = max(float(numpy.abs(rhs).max()), 1e-300)
                resid = float(numpy.abs(arr @ x - rhs).max()) / scale
                fails = not (resid <= self.REFACTOR_RESIDUAL_TOL)
            if fails:
                self.residual_fallbacks += 1
                _held_lost(self)
                old, self._numeric = self._numeric, None
                self._free(old)
                self._numeric = self._factor(Ax)
                self.factors += 1
                x = rhs.copy(order='F')
                if not self._lib.klu_solve(
                        ctypes.c_void_p(self._symbolic), ctypes.c_void_p(self._numeric), self.m,
                        nrhs, self._ptr(x, ctypes.c_double), ctypes.byref(self._common)):
                    raise numpy.linalg.LinAlgError('klu_solve failed')
        return x

    def solve_transposed(self, A, B, toolkit):
        """``A^T X = B`` on the numeric `solve` refactored for `A` -- this
        circuit's matrix; another (`_ours`: None) went to the fallback,
        which answers it (not counted again in `foreign`), or numpy's
        solve where none was made."""
        if self._ours(A) is None:
            if self._fb is None:
                return LinearSolver.solve_transposed(self, A, B, toolkit)
            return self._fb.solve_transposed(A, B, toolkit)
        return KLUSolver.solve_transposed(self, A, B, toolkit)

    def factor(self, A, toolkit):
        """A factorisation of its own (`_KLUNumeric`, freed with it, with the
        transposed solve), on the run's one analysis."""
        got = self._ours(A)
        if got is None:
            return self._other().factor(A, toolkit)
        return self.factor_values(got[1])

    def factor_values(self, Ax):
        """`factor` of the matrix whose values in this pattern's order are
        `Ax`, for a caller that gathered them from a matrix zero off the
        pattern (`_ours`'s test made: `_tran_radau_c.final_block`)."""
        num = self._factor(numpy.ascontiguousarray(Ax, dtype=numpy.float64))
        self.numerics += 1
        fac = _KLUNumeric(self, self._symbolic, num, self.m)
        return _Factored(fac, lambda f, b: f.solve(b),
                         lambda f, b: f.solve_transposed(b), sparse=True)

    def chord_factor(self, A, toolkit):
        """THE CHORD'S FACTOR (`ChordNewton`: one matrix a step, solved again
        at every iteration): `solve`'s sequence on a numeric of the chord's
        own, kept from step to step (`_chord`) -- at its first solve a
        refactor of it on its pivot order (where a fresh factor picks the
        same pivots, the same arithmetic), a fresh factor where none is held
        or the refactor meets a zero pivot, the refactor's residual tested
        as `solve` tests it (a fresh factor where it fails); every later
        solve on the numeric as it stands.  The matrix kept for the test (it
        is the step's; the chord holds one at a time).  Another matrix
        (`_ours`: None) is the fallback's `factor`."""
        got = self._ours(A)
        if got is None:
            return self._other().factor(A, toolkit)
        f = _ChordNumeric(self, got[0], numpy.ascontiguousarray(got[1], dtype=numpy.float64))
        return _Factored(f, lambda f, b: f.solve(b), lambda f, b: f.solve_transposed(b),
                         sparse=True)

    def chord_held(self):
        """The chord's factorisation as a C chord left it (`_tran_newton_c`'s
        serve, `Transient.lte_jacobian`): its numeric as it stands (`_chord`),
        made -- a solve on it, no refactor -- in the wrapper `chord_factor`
        returns."""
        f = _ChordNumeric(self, None, None)
        f._made = True
        return _Factored(f, lambda f, b: f.solve(b), lambda f, b: f.solve_transposed(b),
                         sparse=True)

    def _chord_fresh(self, Ax):
        """The chord's numeric made afresh (the held one freed first: None
        after a singular matrix's LinAlgError, as `_numeric` is)."""
        self._free_chord()
        self._chord = self._factor(Ax)
        self.chord_factors += 1

    def _free_chord(self):
        num, self._chord = self._chord, None
        if num:
            ct = self._ct
            self._lib.klu_free_numeric(ct.byref(ct.c_void_p(num)), ct.byref(self._common))

    def __del__(self):
        try:
            if getattr(self, '_chord', None):
                self._free_chord()
        except Exception:                                     # noqa: BLE001, S110  # pragma: no cover
            pass                                              # interpreter shutdown
        KLUSolver.__del__(self)

    def walk(self):
        """A `_KLUWalk`: `factor_values` for a caller that solves with each
        block's factorisation before it hands the next one and keeps none --
        on one numeric, refactored."""
        return _KLUWalk(self)

    def __repr__(self):
        return (f'StructuralKLU(m={self.m}, nnz={self.nnz}, analyses={self.analyses}, '
                f'factors={self.factors}, refactors={self.refactors}, reused={self.reused}, '
                f'numerics={self.numerics}, '
                f'fallbacks={self.residual_fallbacks}, pivot_fallbacks={self.pivot_fallbacks}, '
                f'foreign={self.foreign}, chord={self.chord_refactors}/{self.chord_factors}/'
                f'{self.chord_pivot_fallbacks}/{self.chord_residual_fallbacks})')


class _KLUNumeric(object):
    """One KLU numeric factorisation, OWNED: freed with this object.  Holds
    its solver (the library and the `klu_common` buffer) and its symbolic
    analysis, which `klu_solve` needs alongside it."""

    def __init__(self, solver, symbolic, numeric, n):
        self._solver = solver
        self._sym, self._num, self.n = symbolic, numeric, n

    def _run(self, fn, b):
        b = numpy.asarray(b)
        if numpy.iscomplexobj(b):
            return self._run(fn, b.real) + 1j * self._run(fn, b.imag)
        x = numpy.array(b, dtype=numpy.float64, order='F', copy=True)
        nrhs = 1 if x.ndim < 2 else x.shape[1]
        ct = self._solver._ct
        if not fn(ct.c_void_p(self._sym), ct.c_void_p(self._num), self.n, nrhs,
                  x.ctypes.data_as(ct.POINTER(ct.c_double)),
                  ct.byref(self._solver._common)):
            raise numpy.linalg.LinAlgError('klu solve failed')
        return x

    def solve(self, b):
        return self._run(self._solver._lib.klu_solve, b)

    def solve_transposed(self, b):
        return self._run(self._solver._lib.klu_tsolve, b)

    def __del__(self):
        num, self._num = getattr(self, '_num', None), None
        if num:
            try:
                ct = self._solver._ct
                self._solver._lib.klu_free_numeric(
                    ct.byref(ct.c_void_p(num)), ct.byref(self._solver._common))
            except Exception:                         # pragma: no cover
                pass                                  # interpreter shutdown


class _ChordNumeric:
    """The chord's factorisation of one step's matrix (`StructuralKLU.
    chord_factor`): its first solve the sequence `solve` makes, on the
    chord's kept numeric; every later one a solve on it.  Owns nothing: the
    numeric is the solver's (`_chord`)."""

    def __init__(self, solver, arr, Ax):
        self._solver, self._arr, self._Ax, self._made = solver, arr, Ax, False

    def _make(self, rhs):
        """`solve`'s refactor path on the chord's numeric for `rhs`: the
        solution."""
        sk = self._solver
        ct, Ax = sk._ct, self._Ax
        reused = sk._chord is not None
        if reused:
            ok = sk._lib.klu_refactor(
                sk._ptr(sk.Ap, ct.c_int), sk._ptr(sk.Ai, ct.c_int), sk._ptr(Ax, ct.c_double),
                ct.c_void_p(sk._analysed()), ct.c_void_p(sk._chord), ct.byref(sk._common))
            sk.chord_refactors += 1
            if not ok:
                ## (the held pivot order met a zero pivot: a fresh factor
                ## decides, as `solve`'s)
                sk.chord_pivot_fallbacks += 1
                sk._chord_fresh(Ax)
                reused = False
        else:
            sk._chord_fresh(Ax)
        x = self._klu(sk._lib.klu_solve, rhs)
        if reused:
            scale = max(float(numpy.abs(rhs).max()), 1e-300)
            resid = float(numpy.abs(self._arr @ x - rhs).max()) / scale
            if not (resid <= sk.REFACTOR_RESIDUAL_TOL):
                sk.chord_residual_fallbacks += 1
                sk._chord_fresh(Ax)
                x = self._klu(sk._lib.klu_solve, rhs)
        return x

    def _klu(self, fn, b):
        sk = self._solver
        ct = sk._ct
        x = numpy.array(b, dtype=numpy.float64, order='F', copy=True)
        nrhs = 1 if x.ndim < 2 else x.shape[1]
        if not fn(ct.c_void_p(sk._symbolic), ct.c_void_p(sk._chord), sk.m, nrhs,
                  x.ctypes.data_as(ct.POINTER(ct.c_double)), ct.byref(sk._common)):
            raise numpy.linalg.LinAlgError('klu solve failed')
        return x

    def solve(self, b):
        b = numpy.asarray(b)
        if numpy.iscomplexobj(b):
            return self.solve(b.real) + 1j * self.solve(b.imag)
        if not self._made:
            self._made = True
            return self._make(numpy.asarray(b, dtype=numpy.float64))
        return self._klu(self._solver._lib.klu_solve, b)

    def solve_transposed(self, b):
        b = numpy.asarray(b)
        if numpy.iscomplexobj(b):
            return self.solve_transposed(b.real) + 1j * self.solve_transposed(b.imag)
        if not self._made:
            ## (the numeric made for this matrix first: on a zero right-hand
            ## side its residual test holds)
            self.solve(numpy.zeros(self._solver.m))
        return self._klu(self._solver._lib.klu_tsolve, b)


class _KLUWalk:
    """ONE NUMERIC FOR A WALK'S BLOCKS (`StructuralKLU.walk`): `factor_values`
    of each block a walk hands it, the first factored, every later one
    refactored on the same numeric -- klu_refactor on the pivot order the
    first factor chose -- and checked as `StructuralKLU.solve` checks its
    refactor, by a residual against `REFACTOR_RESIDUAL_TOL` (of a probe:
    the block times ones, solved), a fresh factor where it fails or meets a
    zero pivot.  For a caller that solves with each block's factorisation
    before it hands the next, and keeps none: a shooting walk's monodromy
    (`_pss_walks._PeriodWalks._walk_stage` without `keep`).  A factorisation
    it handed out answers until the next block's (`_KLUBorrowed`: a later
    use raises); the numeric is freed with this object.  Measured on a
    270-unknown stage block (a Howland testbench's Radau PSS): klu_refactor
    57 us, klu_factor 167 us; the solves within 6e-16 of a fresh factor's.
    `factors`, `refactors`, `fallbacks` count them (each fresh factor in the
    solver's `numerics` too)."""

    __slots__ = ('_cols', '_num', 'factors', 'fallbacks', 'gen', 'refactors', 'solver')

    def __init__(self, solver):
        self.solver = solver
        self._num = None
        self.gen = 0
        ## (each CSC entry's column, for the probe's product)
        self._cols = numpy.repeat(numpy.arange(solver.m), numpy.diff(solver.Ap))
        self.factors = self.refactors = self.fallbacks = 0

    def factor_values(self, Ax):
        sk = self.solver
        ct = sk._ct
        Ax = numpy.ascontiguousarray(Ax, dtype=numpy.float64)
        num = self._num
        self.gen += 1
        if num is not None:
            ok = sk._lib.klu_refactor(
                sk._ptr(sk.Ap, ct.c_int), sk._ptr(sk.Ai, ct.c_int), sk._ptr(Ax, ct.c_double),
                ct.c_void_p(sk._analysed()), ct.c_void_p(num), ct.byref(sk._common))
            self.refactors += 1
            if not (ok and self._holds(Ax)):
                self.fallbacks += 1
                self._num = None
                sk._free(num)
                num = None
        if num is None:
            num = self._num = sk._factor(Ax)
            self.factors += 1
            ## (a numeric of its own, as `factor_values` makes: counted there)
            sk.numerics += 1
        fac = _KLUBorrowed(self, num)
        return _Factored(fac, lambda f, b: f.solve(b),
                         lambda f, b: f.solve_transposed(b), sparse=True)

    def _holds(self, Ax):
        """The refactored numeric's check: ``b = A 1``, ``x = A^-1 b``, ``max
        |A x - b| / max(max |b|, 1e-300)`` within `REFACTOR_RESIDUAL_TOL` (a
        NaN fails); a failed solve fails."""
        sk = self.solver
        ct = sk._ct
        b = numpy.bincount(sk.Ai, weights=Ax, minlength=sk.m)
        x = b.copy()
        if not sk._lib.klu_solve(ct.c_void_p(sk._symbolic), ct.c_void_p(self._num), sk.m, 1,
                                 sk._ptr(x, ct.c_double), ct.byref(sk._common)):
            return False
        r = numpy.bincount(sk.Ai, weights=Ax * x[self._cols], minlength=sk.m) - b
        scale = max(float(numpy.abs(b).max()), 1e-300)
        return float(numpy.abs(r).max()) / scale <= sk.REFACTOR_RESIDUAL_TOL

    def __del__(self):
        num, self._num = getattr(self, '_num', None), None
        if num:
            try:
                self.solver._free(num)
            except Exception:                         # noqa: BLE001, S110  # pragma: no cover
                pass                                  # interpreter shutdown


class _KLUBorrowed(_KLUNumeric):
    """A `_KLUWalk`'s numeric as one block's factorisation: it answers while
    the walk has not refactored it for another (its generation), raises
    after, and frees nothing."""

    def __init__(self, walk, num):
        sk = walk.solver
        _KLUNumeric.__init__(self, sk, sk._symbolic, num, sk.m)
        self._walk, self._gen = walk, walk.gen

    def _run(self, fn, b):
        if self._walk.gen != self._gen or self._walk._num != self._num:
            raise RuntimeError('a walk\'s KLU factorisation used after the next block\'s')
        return _KLUNumeric._run(self, fn, b)

    def __del__(self):
        self._num = None


class ComplexKLUSolver(object):
    """SuiteSparse KLU for COMPLEX systems, through the ``klu_z_*`` family.

    The complex sibling of :class:`KLUSolver`, and the enabling dependency for
    the Radau IIA(3) cost transform: the coupled ``3m`` stage solve
    block-diagonalises through ``eig(A^{-1})`` into one REAL factorisation of
    ``(gamma_r/h) C + G`` and one COMPLEX factorisation of
    ``((alpha + i beta)/h) C + G`` (its conjugate is free), replacing an
    ``O((3m)^3)`` dense solve with two sparse ones.  Without a complex sparse
    LU the transform has no home and a real-only fallback throws the whole
    advantage away (the peer's item C); this provides it.

    ⚠ SuiteSparse's ``klu_analyze`` is PATTERN-ONLY and shared with the real
    path -- the symbolic ordering does not depend on the values being real or
    complex -- so only ``klu_z_factor``/``klu_z_refactor``/``klu_z_solve`` are
    complex.  KLU's packed-complex layout is ``[re0, im0, re1, im1, ...]``,
    which is exactly how numpy stores ``complex128``, so the values marshal by
    ``arr.view(float64)`` with no copy or interleave.  Not a `LinearSolver`
    subclass: that ABC's contract is real (`factor`/transpose semantics tuned
    for the shooting replays), and this solves complex systems for the Radau
    transform only.

    KLU's objects are freed where they are replaced -- the analysis and the
    factor at a new pattern, the factor at a residual fallback -- and at the
    solver's end (`analyses_freed`, `factors_freed` count the first two).
    Whoever holds their handles (`_tran_radau_tc`, the C stepping loop's
    transform mode) compares them with the solver's before each use.
    """

    REFACTOR_RESIDUAL_TOL = 1e-8

    def __init__(self, libname='libklu.so.2'):
        import ctypes
        import numpy
        try:
            lib = ctypes.CDLL(libname)
        except OSError as e:
            raise ImportError(
                'ComplexKLUSolver needs SuiteSparse KLU (%s); install libklu / '
                'libsuitesparse' % libname) from e
        ## (scipy.sparse imported on first use: no served run reaches it)
        if importlib.util.find_spec('scipy') is None:         # pragma: no cover
            raise ImportError('ComplexKLUSolver needs scipy.sparse')
        _sp = _lazy.deferred('scipy.sparse')
        self._ct = ctypes
        self._np = numpy
        self._sp = _sp
        self._lib = lib
        c_int_p = ctypes.POINTER(ctypes.c_int)
        c_dbl_p = ctypes.POINTER(ctypes.c_double)
        for nm, argt, rest in (
                ('klu_defaults', [ctypes.c_void_p], ctypes.c_int),
                ('klu_analyze', [ctypes.c_int, c_int_p, c_int_p, ctypes.c_void_p],
                 ctypes.c_void_p),
                ('klu_z_factor', [c_int_p, c_int_p, c_dbl_p, ctypes.c_void_p,
                                  ctypes.c_void_p], ctypes.c_void_p),
                ('klu_z_refactor', [c_int_p, c_int_p, c_dbl_p, ctypes.c_void_p,
                                    ctypes.c_void_p, ctypes.c_void_p], ctypes.c_int),
                ('klu_z_solve', [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int,
                                 ctypes.c_int, c_dbl_p, ctypes.c_void_p],
                 ctypes.c_int),
                ('klu_free_symbolic', [ctypes.c_void_p, ctypes.c_void_p], ctypes.c_int),
                ('klu_z_free_numeric', [ctypes.c_void_p, ctypes.c_void_p], ctypes.c_int)):
            f = getattr(lib, nm)
            f.argtypes = argt
            f.restype = rest
        self._common = ctypes.create_string_buffer(4096)
        if not lib.klu_defaults(ctypes.byref(self._common)):
            raise ImportError('klu_defaults failed; the KLU library looks unusable')
        self._pattern = None
        self._symbolic = None
        self._numeric = None
        ## the `prepare` record whose values the numeric holds as a REFACTOR
        ## (`solve_prepared`), or None
        self._fresh = None
        ## (the last record's pattern, `_csc_of_dense`'s)
        self._csc_last = None
        ## (the structural pattern's arrays and key, made once a structure:
        ## `prepare_structural`)
        self._skey = None
        self.analyses = 0
        self.factors = 0
        self.refactors = 0
        self.residual_fallbacks = 0
        ## (refactors on the structural pattern that met a zero pivot: a
        ## fresh factor made instead)
        self.pivot_fallbacks = 0
        ## (the analyses and factors freed where they were replaced)
        self.analyses_freed = 0
        self.factors_freed = 0

    def _ptr(self, arr, ctype):
        return arr.ctypes.data_as(self._ct.POINTER(ctype))

    def _free_numeric(self):
        """The factor held, freed (`klu_z_free_numeric`): None after."""
        num, self._numeric = self._numeric, None
        if num:
            ct = self._ct
            self._lib.klu_z_free_numeric(ct.byref(ct.c_void_p(num)), ct.byref(self._common))
            self.factors_freed += 1

    def _free_symbolic(self):
        """The analysis held, freed (`klu_free_symbolic`): None after."""
        sym, self._symbolic = self._symbolic, None
        if sym:
            ct = self._ct
            self._lib.klu_free_symbolic(ct.byref(ct.c_void_p(sym)), ct.byref(self._common))
            self.analyses_freed += 1

    def __del__(self):
        try:
            if getattr(self, '_numeric', None):
                self._free_numeric()
            if getattr(self, '_symbolic', None):
                self._free_symbolic()
        except Exception:                                     # noqa: BLE001, S110  # pragma: no cover
            pass                                              # interpreter shutdown

    def prepare(self, A):
        """``A`` as `solve` marshals it -- its product with a vector (the
        residual check's), its CSC column pointers, row indices and packed
        values, and the pattern key -- for a caller solving several
        right-hand sides against ONE matrix (`solve_prepared`: the Radau
        transform's frozen complex factor).

        A dense square complex128 `A` is marshalled by numpy
        (`_csc_of_dense`, speed round 10, B3.1): the arrays below, and the
        product as SciPy runs it.  Anything else goes through SciPy."""
        numpy = self._np
        mv = _csc_matvec()
        if (mv is not None and type(A) is numpy.ndarray and A.dtype == numpy.complex128
                and A.ndim == 2 and A.shape[0] == A.shape[1]):
            rec, self._csc_last = _csc_of_dense(A, mv, self._csc_last)
            return rec
        ## accept either a dense array or a scipy.sparse matrix
        Acsc = self._sp.csc_matrix(A).astype(numpy.complex128)
        n = Acsc.shape[0]
        Ap = numpy.ascontiguousarray(Acsc.indptr, dtype=numpy.int32)
        Ai = numpy.ascontiguousarray(Acsc.indices, dtype=numpy.int32)
        ## packed complex: view the complex128 data as 2*nnz interleaved doubles
        Ax = numpy.ascontiguousarray(Acsc.data, dtype=numpy.complex128).view(numpy.float64)
        return (Acsc.dot, n, Ap, Ai, Ax, (n, Ap.tobytes(), Ai.tobytes()))

    #: a record on an analysis solver's structural pattern
    #: (`prepare_structural`): its key's first entry
    STRUCTURAL = 'structural'

    def prepare_structural(self, A, sk):
        """`prepare(A)` on the analysis solver's STRUCTURAL pattern -- `sk`,
        a `StructuralKLU` that has made its analysis -- for the Radau
        transform's complex factor (`_RadauStages._radau_frozen`, roadmap
        259): `A`'s entries at the structure's positions (`take`), explicit
        zeros among them, in its order; the record's `Ap`, `Ai` and key the
        structure's, the same at every step -- one analysis a run, refactors
        after, where a pattern of nonzeros changes as a device's entries
        turn exactly zero.  None where `A` is not a dense square complex128
        of the structure's size, has a nonzero off it (NaN one), or SciPy's
        product is missing: the caller then takes `prepare`."""
        numpy = self._np
        mv = _csc_matvec()
        m = sk.m
        if (mv is None or type(A) is not numpy.ndarray or A.dtype != numpy.complex128
                or A.shape != (m, m)):
            return None
        sr = self._skey
        if sr is None or sr[0] is not sk or sr[1] is not sk.Ap or sr[2] is not sk.Ai:
            Ap = numpy.ascontiguousarray(sk.Ap, dtype=numpy.int32)
            Ai = numpy.ascontiguousarray(sk.Ai, dtype=numpy.int32)
            take = numpy.ascontiguousarray(sk.take, dtype=numpy.int64)
            sr = self._skey = (sk, sk.Ap, sk.Ai, Ap, Ai, take,
                               (self.STRUCTURAL, m, Ap.tobytes(), Ai.tobytes()))
        _sk, _ap, _ai, Ap, Ai, take, key = sr
        data = A.ravel()[take]
        if numpy.count_nonzero(A) != numpy.count_nonzero(data):
            return None
        return (functools.partial(_csc_dot, mv, m, Ap, Ai, data), m, Ap, Ai,
                data.view(numpy.float64), key)

    def prepare_structural_values(self, data, sk):
        """`prepare_structural(A, sk)` given `A`'s entries at the structure's
        positions, `data` (complex128, in its order: `take`'s) -- for a
        caller that made them there alone, `A` zero off the structure: the
        same record, without the gather and the count.  None where SciPy's
        product is missing or `data` is not of the structure's size."""
        numpy = self._np
        mv = _csc_matvec()
        if (mv is None or type(data) is not numpy.ndarray or data.dtype != numpy.complex128
                or data.shape != (sk.nnz,) or not data.flags.c_contiguous):
            return None
        m = sk.m
        sr = self._skey
        if sr is None or sr[0] is not sk or sr[1] is not sk.Ap or sr[2] is not sk.Ai:
            Ap = numpy.ascontiguousarray(sk.Ap, dtype=numpy.int32)
            Ai = numpy.ascontiguousarray(sk.Ai, dtype=numpy.int32)
            take = numpy.ascontiguousarray(sk.take, dtype=numpy.int64)
            sr = self._skey = (sk, sk.Ap, sk.Ai, Ap, Ai, take,
                               (self.STRUCTURAL, m, Ap.tobytes(), Ai.tobytes()))
        _sk, _ap, _ai, Ap, Ai, _take, key = sr
        return (functools.partial(_csc_dot, mv, m, Ap, Ai, data), m, Ap, Ai,
                data.view(numpy.float64), key)

    def prepare_packed(self, data, m, rows, cols):
        """`prepare(A)` given `A`'s entries `data` (complex128) at the
        positions `rows`, `cols` -- column by column, rows ascending -- `A`
        zero off them: the record `_csc_of_dense` makes of `A` (its nonzeros,
        a NaN one, a signed zero not, their pattern; the last record's
        arrays where the pattern is its), the last pattern kept as it keeps
        it.  None where SciPy's product is missing."""
        numpy = self._np
        mv = _csc_matvec()
        if mv is None or type(data) is not numpy.ndarray or data.dtype != numpy.complex128:
            return None
        rec, self._csc_last = _csc_of_packed(data, m, rows, cols, mv, self._csc_last)
        return rec

    def prepare_values(self, data, A):
        """`prepare(A)` where `A`'s nonzeros are known to fall in the last
        record's pattern, `data` holding them in its order (speed round 12,
        stage 5b: the transform's frozen factors, made in C with their
        packed values): `prepare`'s own record, without the mask and the
        gather -- and `prepare(A)` itself where it would not marshal by
        numpy or has no last pattern."""
        numpy = self._np
        mv = _csc_matvec()
        last = self._csc_last
        if (mv is None or last is None or type(A) is not numpy.ndarray
                or A.dtype != numpy.complex128 or A.ndim != 2 or A.shape[0] != A.shape[1]
                or last[0].shape != A.shape):
            return self.prepare(A)
        rec, self._csc_last = _csc_of_values(data, mv, last)
        return rec

    def solve(self, A, b):
        """Solve the complex system ``A x = b`` and return ``x`` (complex).

        Reuses the symbolic ordering across calls with the same sparsity
        pattern (``klu_z_refactor``), falling back to a full ``klu_z_factor``
        when the residual shows the reused pivots have gone bad -- the same
        validated-not-trusted rule :class:`KLUSolver` uses.
        """
        return self.solve_prepared(self.prepare(A), b)

    def solve_prepared(self, prep, b):
        """`solve` for the matrix `prepare` marshalled into `prep`.

        ONE MATRIX, SEVERAL RIGHT-HAND SIDES (speed round 10, B3.1): where
        the numeric factorisation already holds a REFACTOR of this very
        `prep`'s values (`_fresh`), the refactor `solve` would repeat is not
        made -- the same values with the same pivots give the same factors,
        so the same solution bits (measured 8000 of 8000 against a refactor
        every call).  A full factor is always followed by one refactor at
        the next call, as `solve` would make it; the residual check and its
        fallback run on every reused solve, as before.  `refactors` counts
        the refactors made."""
        numpy = self._np
        ctypes = self._ct
        dot, n, Ap, Ai, Ax, key = prep
        rhs = numpy.ascontiguousarray(numpy.asarray(b, dtype=numpy.complex128).ravel())

        reused = key == self._pattern and self._symbolic and self._numeric
        if not reused:
            ## (the last pattern's factor and analysis freed first)
            self._free_numeric()
            self._free_symbolic()
            self._symbolic = self._lib.klu_analyze(
                n, self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                ctypes.byref(self._common))
            self.analyses += 1
            self._fresh = None
            if not self._symbolic:
                raise numpy.linalg.LinAlgError('klu_analyze failed (singular structure?)')
            self._numeric = self._lib.klu_z_factor(
                self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._symbolic),
                ctypes.byref(self._common))
            self.factors += 1
            if not self._numeric:
                raise numpy.linalg.LinAlgError('Singular matrix')
            self._pattern = key
        elif self._fresh is not prep:
            self._fresh = None
            ok = self._lib.klu_z_refactor(
                self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._symbolic),
                ctypes.c_void_p(self._numeric), ctypes.byref(self._common))
            self.refactors += 1
            if not ok:
                if key[0] != self.STRUCTURAL:
                    raise numpy.linalg.LinAlgError('Singular matrix')
                ## (on the structural pattern a pivot the last fresh factor
                ## chose may be an explicit zero now: a fresh factor on the
                ## analysis decides, as `StructuralKLU`'s refactor does)
                self.pivot_fallbacks += 1
                self._factor_afresh(Ap, Ai, Ax)
            self._fresh = prep

        x = rhs.copy()
        xv = x.view(numpy.float64)
        if not self._lib.klu_z_solve(
                ctypes.c_void_p(self._symbolic), ctypes.c_void_p(self._numeric),
                n, 1, self._ptr(xv, ctypes.c_double), ctypes.byref(self._common)):
            raise numpy.linalg.LinAlgError('klu_z_solve failed')

        if reused:
            scale = max(float(numpy.abs(rhs).max()), 1e-300)
            resid = float(numpy.abs(dot(x) - rhs).max()) / scale
            if not (resid <= self.REFACTOR_RESIDUAL_TOL):
                self.residual_fallbacks += 1
                self._fresh = None
                self._factor_afresh(Ap, Ai, Ax)
                x = rhs.copy()
                xv = x.view(numpy.float64)
                if not self._lib.klu_z_solve(
                        ctypes.c_void_p(self._symbolic), ctypes.c_void_p(self._numeric),
                        n, 1, self._ptr(xv, ctypes.c_double), ctypes.byref(self._common)):
                    raise numpy.linalg.LinAlgError('klu_z_solve failed')
        return x

    def _factor_afresh(self, Ap, Ai, Ax):
        """A fresh factor of the values on the analysis held, the factor held
        freed first -- `solve_prepared`'s fallbacks: a residual the reused
        pivots failed, a structural refactor's zero pivot; numpy's error
        where the matrix is singular."""
        ctypes = self._ct
        self._free_numeric()
        self._numeric = self._lib.klu_z_factor(
            self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
            self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._symbolic),
            ctypes.byref(self._common))
        self.factors += 1
        if not self._numeric:
            raise self._np.linalg.LinAlgError('Singular matrix')

    def __repr__(self):
        return ('ComplexKLUSolver(analyses=%d, factors=%d, refactors=%d, '
                'fallbacks=%d, pivot_fallbacks=%d)' % (self.analyses, self.factors, self.refactors,
                                                       self.residual_fallbacks, self.pivot_fallbacks))


## Below this many unknowns the sparse path cannot repay its own setup, and
## measuring the fill is not worth the pass over the matrix.  Chosen from the
## measured table (see the module docstring): at n=61 the win is 1.16x, which is
## inside the noise of a Python-level dispatch.
## ⚠ RAISED 100 -> 250 (2026-10-02), END TO END: once the assembly stopped
## re-stamping constant elements (`_stamp_plan`) the transient's own time
## shrank around the solve, and SuperLU's per-call cost decides.  A gear
## transient on an RC ladder with a diode every ten sections, best of two,
## AutoSolver against the dense path: 0.83x at n = 102, 0.88x at 152, 0.92x at
## 202, 1.00x at 252, 1.05x at 302, 1.17x at 402, 1.32x at 602, 1.60x at 1202
## (answers within 3e-17).  The floor sits at the crossover.
MIN_N_FOR_SPARSE = 250

## And above this fill the matrix is dense enough that SuperLU's bookkeeping costs
## more than it saves.  The measured circuits run 0.27%-9.4%; 20% is a deliberately
## conservative bound rather than a tuned one, because no measurement here explores
## dense circuits and a threshold fitted to sparse ones would be a guess dressed as
## a number.
MAX_FILL_FOR_SPARSE = 0.20

## KLU ON THE CIRCUIT'S STRUCTURAL PATTERN (`StructuralKLU`) from this many
## unknowns, where an analysis asks for it (`structural_min_n`: the
## transient's, `Transient._STRUCTURAL_KLU_MIN_N`; every other analysis keeps
## `MIN_N_FOR_SPARSE`).  `MIN_N_FOR_SPARSE` was measured on the Python loop
## with SuperLU; the C stepping loop solves on the structure in C -- the chord
## its own factor each step, the full Newton a refactor of the run's numeric
## each iteration (newton.c's `sparse_full`, 2026-10-16) -- and its crossover
## sits far lower.  Measured (2026-10-16, the 27 gold decks, three rounds
## interleaved, medians): the chord's transient on KLU against dense, every
## deck from 40 unknowns faster -- mike2 (40) 0.97x, latch (65) 0.75x, ring
## (87) 0.78x, gm1 (129) 0.62x, mux8 (170) 0.61x, voter25 (199) 0.53x, jge
## (243) 0.64x -- the 27 decks' transients 77.0 -> ~52.6 s; and the full
## Newton's solve, timed per deck on its own Jacobians, 6.8 -> 2.8 us at 40
## unknowns, 38.5 -> 6.5 at 97, 692 -> 25 at 243 (dense LU goes threaded from
## ~100 unknowns and costs 5x its cubic trend).  Below 40 (rca, 33; gm2, 22;
## mosrect, 19) dense, where KLU was measured slower or not measured.  A
## matrix the structure does not cover (`StructuralKLU`'s fallback), and a
## run with no structure, keep `MIN_N_FOR_SPARSE`.
MIN_N_FOR_STRUCTURAL_KLU = 40

## KLU on the structure from this many unknowns for DC (and so every
## analysis's operating point) and the stage methods (TRBDF2, ESDIRK43,
## Radau); 40 was measured and refused for them (2026-10-19, roadmap 153:
## below ~100 unknowns KLU's fixed cost lost up to +31 %); the GLMs keep
## `MIN_N_FOR_SPARSE`.  Measured (2026-10-20, the five gold decks it moves,
## three rounds interleaved, best of 3): DC's operating points 927 -> 328 ms,
## every deck faster (toronto -40 %, jge -83 %), within 2e-9 V; TRBDF2's
## transients -15.1 % (each deck -3 .. -18 %), ESDIRK43's -18.9 % (-11 ..
## -23 %), no gold column changed pass/fail.  Radau kept `MIN_N_FOR_SPARSE`
## until its C transformed solve took the KLU real factor (roadmap 215);
## measured 2026-10-24 (roadmap 219, the 16 gold decks from 40 to 250
## unknowns, three rounds interleaved, medians): at 150 the five decks it
## moves 56.0 -> 45.0 s (1.25x; toronto 1.22x, mux8 1.22x, voter25 1.19x,
## e1480 1.10x, jge 1.28x), their gold metrics as before; at 40, 1.22x but
## g1310 (98) 0.94x, todd3 0.97x, latch 0.98x; nothing measured between 98
## and 130.
MIN_N_FOR_STRUCTURAL_KLU_OP_AND_STAGES = 150

## Above this many unknowns KLU's reused ordering repays its per-call cost and it
## overtakes SuperLU; below it, SuperLU's lighter setup wins.  Measured (see
## `AutoSolver._select`): KLU is 0.52x the dense path at n=152, 1.00x at n=402,
## **1.53x at n=1002** and 1.47x at n=2002.  1000 is where the tables cross, not a
## round number chosen for looks.
MIN_N_FOR_KLU = 1000


class AutoSolver(LinearSolver):
    """Pick dense or sparse per matrix, on measured fill rather than on ``n``.

    The decision is cached on the (shape, fill-bucket) of the first matrix seen,
    because a transient calls this thousands of times with the same structure and
    counting non-zeros every time would give back the saving.
    """

    def __init__(self, min_n=MIN_N_FOR_SPARSE, max_fill=MAX_FILL_FOR_SPARSE, structure=None,
                 structural_min_n=None, stamped=None):
        self.min_n = min_n
        ## (with `structure`: which of its places an element stamps --
        ## `StructuralKLU`'s `stamped`, `analysis.stamped_places`)
        self.stamped = stamped
        ## (the floor of `StructuralKLU` alone, where a structure is given:
        ## `MIN_N_FOR_STRUCTURAL_KLU` -- `min_n` unless an analysis asks)
        self.structural_min_n = min_n if structural_min_n is None else structural_min_n
        self.max_fill = max_fill
        ## a callable giving the circuit's structural pattern -- `(m, rows,
        ## cols)` or None (`analysis.jacobian_structure`): an analysis hands
        ## it in, and a sparse matrix of its size then factors on it with KLU
        ## (`StructuralKLU`; Andreas's decision, 2026-10-10)
        self.structure = structure
        self._dense = DenseSolver()
        self._sparse = None
        self._choice = None
        ## (`held`'s solvers, one per set of rows dropped)
        self._held = {}

    def _select(self, A):
        if self._choice is not None:
            return self._choice
        from pycircuit.circuit import _packed, _paths
        if type(A) is _packed.PatternMatrix:
            ## (THE PACKED LAYOUT: KLU on the circuit's structure whatever the
            ## size or the fill -- the dense solvers need the dense matrix)
            n_unknowns = A.shape[0]
            st = self.structure() if self.structure is not None else None
            if st is not None and st[0] == n_unknowns and self._structural(st, n_unknowns):
                return self._choice
            raise TypeError('PatternMatrix: no structural pattern to factor it on')
        arr = numpy.asarray(A)
        floor = (min(self.min_n, self.structural_min_n) if self.structure is not None
                 else self.min_n)
        if arr.ndim != 2 or arr.shape[0] < floor or arr.dtype == object:
            _paths.COUNTS['solver.auto:dense_small'] += 1
            self._choice = self._dense
            return self._choice
        fill = _count_nonzero(arr) / float(arr.shape[0] * arr.shape[1])
        if fill > self.max_fill:
            _paths.COUNTS['solver.auto:dense_fill'] += 1
            self._choice = self._dense
            return self._choice
        n_unknowns = arr.shape[0]
        ## WHICH sparse solver depends on `n`, and the crossover is measured.
        ##
        ## KLU wins the isolated factor+solve 3.5x-10x, but end to end it is a
        ## different story, because the solve is only a fraction of a transient and
        ## KLU pays a fixed per-call cost.  Measured on a transient, best-of-2,
        ## against the shipped dense path:
        ##
        ##      n     dense    SuperLU        KLU   SuperLU   KLU
        ##    152    3.581 s    4.378 s    6.832 s     0.82x  0.52x
        ##    402   12.413 s    9.524 s   12.444 s     1.30x  1.00x
        ##   1002   53.872 s   43.462 s   35.263 s     1.24x  1.53x
        ##   2002  161.484 s  117.111 s  110.143 s     1.38x  1.47x
        ##
        ## So KLU LOSES below ~400, ties there, and wins from ~1000 -- where it also
        ## beats SuperLU.  Selecting it everywhere would make a 152-unknown circuit
        ## twice as slow; never selecting it leaves 1.5x on the table at n >= 1000.
        ##
        ## The n >= 1002 rows could not be measured until stage 2+.5 made circuits
        ## that size buildable: construction was O(N^2.23) and a 1600-section ladder
        ## exceeded a ten-minute budget.  The earlier reading -- "KLU loses end to
        ## end" -- was correct for the sizes then reachable, and this supersedes it.
        ## THE CIRCUIT'S PATTERN, where the analysis gave one: KLU on it for
        ## every sparse run (the C stepping loop plan's stage 2) -- today's
        ## choice below for a matrix it does not cover, and where none
        ## (from `structural_min_n` -- the transient's floor is lower than
        ## `min_n`: `MIN_N_FOR_STRUCTURAL_KLU`)
        st = (self.structure() if self.structure is not None
              and n_unknowns >= self.structural_min_n else None)
        if st is not None and st[0] == n_unknowns and self._structural(st, n_unknowns):
            return self._choice
        if n_unknowns < self.min_n:
            ## (below the floor a sparse solve on the values' pattern was
            ## measured for: dense, as where no structure lowers it)
            _paths.COUNTS['solver.auto:dense_small'] += 1
            self._choice = self._dense
            return self._choice
        sp = self._value_sparse(n_unknowns)
        if sp is not None:
            self._sparse = self._choice = sp
            return self._choice
        self._choice = self._dense                          # pragma: no cover
        return self._choice

    def _structural(self, st, n_unknowns):
        """The StructuralKLU on the structure `st` chosen (its count): True,
        or False where KLU is missing."""
        from pycircuit.circuit import _paths
        try:
            self._sparse = self._choice = StructuralKLU(
                *st, fallback=lambda: (self._value_sparse(n_unknowns, count=False)
                                       if n_unknowns >= self.min_n else self._dense),
                stamped=self.stamped)
        except ImportError:
            return False
        _paths.COUNTS['solver.auto:StructuralKLU'] += 1
        return True

    def choose_on_structure(self, m):
        """`_select`'s choice for an m x m matrix, made before one is seen
        where no value of it can change the choice -- or None (the first
        matrix chooses).  Below the floor the size alone chooses: dense.
        From it, a structure of size m whose own fill is at most `max_fill`
        -- a matrix whose nonzeros lie in the structure fills no more, so
        `_select`'s fill test passes, its quotient the same form -- chooses
        the StructuralKLU on it, as `_select` makes it.  Taken by the C
        Newton (`_tran_newton_c`), which serves a first attempt only on a
        choice made."""
        if self._choice is not None:
            return self._choice
        from pycircuit.circuit import _paths
        floor = (min(self.min_n, self.structural_min_n) if self.structure is not None
                 else self.min_n)
        if m < floor:
            _paths.COUNTS['solver.auto:dense_small'] += 1
            self._choice = self._dense
            return self._choice
        if self.structure is None or m < self.structural_min_n:
            return None
        st = self.structure()
        if (st is None or st[0] != m or len(st[1]) / float(m * m) > self.max_fill
                or not self._structural(st, m)):
            return None
        return self._choice

    def _value_sparse(self, n_unknowns, count=True):
        """Today's sparse solver for `n_unknowns`, its pattern from the
        values: KLU from `MIN_N_FOR_KLU`, SuperLU below (the other where one
        is missing), or None."""
        from pycircuit.circuit import _paths
        for cls in ((KLUSolver, SuperLUSolver) if n_unknowns >= MIN_N_FOR_KLU
                    else (SuperLUSolver, KLUSolver)):
            try:
                sp = cls()
            except ImportError:
                continue
            if count:
                _paths.COUNTS['solver.auto:' + cls.__name__] += 1
            return sp
        return None

    def solve(self, A, b, toolkit):
        return self._select(A).solve(A, b, toolkit)

    def solve_kept(self, A, b, toolkit):
        """The selected solver's (`LinearSolver.solve_kept`)."""
        return self._select(A).solve_kept(A, b, toolkit)

    def solve_transposed(self, A, B, toolkit):
        """The selected solver's (`LinearSolver.solve_transposed`)."""
        return self._select(A).solve_transposed(A, B, toolkit)

    def held(self, drop):
        """The AutoSolver of this one's matrices with the rows and columns
        `drop` removed -- a DC's held nodes, which PCNR eliminates: this
        one's floors, and its structure without them, renumbered (the
        reduced matrix's own pattern).  One per set of rows, kept, as this
        one is: its choice and its factorisation outlive a solve.  Before
        it (2026-10-16) the held solve's smaller matrix was not the
        structure's: this one chose SuperLU on it and kept that choice for
        the release too, and PCNR's C declined every solve
        (`pcnr_c:linsolver`): rich3's and arom's operating points."""
        key = tuple(sorted({int(i) for i in drop}))
        got = self._held.get(key)
        if got is None:
            parent = self.structure

            def structure():
                st = parent()
                if st is None:
                    return None
                m, rows, cols = st
                keep = numpy.ones(m, dtype=bool)
                keep[list(key)] = False
                new = numpy.cumsum(keep) - 1
                rows, cols = numpy.asarray(rows), numpy.asarray(cols)
                on = keep[rows] & keep[cols]
                return int(keep.sum()), new[rows[on]], new[cols[on]]
            got = self._held[key] = AutoSolver(
                self.min_n, self.max_fill, structure=structure if parent is not None else None,
                structural_min_n=self.structural_min_n)
        return got

    def factor(self, A, toolkit):
        """The selected solver's factorisation.  ⚠ Until 2026-09-25 this
        inherited the base `None`, so the shooting replays under
        `AutoSolver` had no transposed solve (gear's `ppv` raised)."""
        return self._select(A).factor(A, toolkit)

    def chord_factor(self, A, toolkit):
        """The selected solver's chord factorisation (`StructuralKLU`'s
        kept numeric)."""
        return self._select(A).chord_factor(A, toolkit)

    def __repr__(self):
        return ('AutoSolver(min_n=%d, structural_min_n=%d, max_fill=%.2f, chose=%r)'
                % (self.min_n, self.structural_min_n, self.max_fill, self._choice))
