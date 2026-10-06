"""Small reduces and inserts as one indexed copy (speed round 11):
`analysis._reduce_small` against the slices of `_reduce_ndarray`, and
`_insert_small` against `insert_row`'s -- on drawn and enumerated arrays of
every dtype and layout, sizes across the cutoffs, every row, special values
(signed zeros, NaN payloads quiet and signalling, infinities, subnormals):
the same bytes, dtype and shape, a fresh writable C-contiguous array;
rows out of range and rows not an int take the slices, raising as they
did."""
import contextlib

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import analysis, circuit


@contextlib.contextmanager
def _slices_only():
    """The cutoffs below every size: each call takes the slices."""
    old = analysis._TAKE_2D, analysis._TAKE_1D
    analysis._TAKE_2D = analysis._TAKE_1D = -1
    try:
        yield
    finally:
        analysis._TAKE_2D, analysis._TAKE_1D = old


def _outcome(fn, *a):
    try:
        return fn(*a)
    except Exception as e:                                     # noqa: BLE001
        return type(e)


def _check(fn, src, *a):
    """`fn(src, *a)` by the indexed copy and by the slices: the same."""
    new = _outcome(fn, src, *a)
    with _slices_only():
        ref = _outcome(fn, src, *a)
    if isinstance(ref, type) or ref is None:
        assert new is ref, (new, ref)
        return
    assert type(new) is np.ndarray and new.dtype == ref.dtype and new.shape == ref.shape
    assert new.tobytes() == ref.tobytes()
    assert new.flags.c_contiguous and new.flags.writeable and new.flags.owndata
    assert not np.shares_memory(new, src)


def _bits(b):
    return float(np.frombuffer(np.uint64(b).tobytes(), dtype=np.float64)[0])


SPECIAL = (0.0, -0.0, np.inf, -np.inf, 5e-324, -5e-324, 1.7976931348623157e308,
           _bits(0x7FF8000000000000), _bits(0x7FF8000000000123), _bits(0xFFF8000000000456),
           _bits(0x7FF0000000000001), _bits(0xFFF4000000000007), 1.0, -2.5)
DTYPES = (np.float64, np.complex128, np.float32, np.int64, np.bool_, '>f8', object)


def _array(values, dtype, shape):
    a = np.array(values, dtype=np.float64).reshape(shape)
    if dtype is np.complex128:
        return a + 1j * a[::-1].copy().reshape(shape)
    if dtype is object:
        return np.array([int(v) if np.isfinite(v) else 0 for v in a.ravel()],
                        dtype=object).reshape(shape)
    return a.astype(dtype) if dtype is not np.float64 else a


def _layouts(a):
    """`a` as it is, in Fortran order, and as a strided and a transposed view."""
    yield a
    yield np.asfortranarray(a)
    big = np.zeros(tuple(2 * s for s in a.shape), dtype=a.dtype)
    big[tuple(slice(None, None, 2) for _ in a.shape)] = a
    yield big[tuple(slice(None, None, 2) for _ in a.shape)]
    if a.ndim == 2:
        yield a.T


def _reduce(A, n):
    return analysis.remove_row_col((A,), n, circuit.numeric)[0]


@settings(deadline=None, max_examples=150)
@given(data=st.data())
def test_a_drawn_reduce_is_the_slices(data):
    ndim = data.draw(st.sampled_from((1, 2)))
    N = data.draw(st.integers(1, 40) if ndim == 2 else st.integers(1, 300))
    shape = (N, N) if ndim == 2 else (N,)
    size = N * N if ndim == 2 else N
    vals = data.draw(st.lists(st.one_of(st.floats(allow_nan=False), st.sampled_from(SPECIAL)),
                              min_size=size, max_size=size))
    dtype = data.draw(st.sampled_from(DTYPES))
    n = data.draw(st.one_of(st.integers(0, N - 1), st.sampled_from((N, N + 1, -1))))
    with np.errstate(all='ignore'):
        A = _array(vals, dtype, shape)
    for a in _layouts(A):
        _check(analysis._reduce_ndarray, a, n)
        _check(analysis._reduce_ndarray, a, np.int64(n))
        if dtype is not object:
            _check(_reduce, a, n)


@settings(deadline=None, max_examples=150)
@given(data=st.data())
def test_a_drawn_insert_is_the_slices(data):
    N = data.draw(st.integers(0, 300))
    vals = data.draw(st.lists(st.one_of(st.floats(allow_nan=False), st.sampled_from(SPECIAL)),
                              min_size=N, max_size=N))
    n = data.draw(st.one_of(st.integers(0, N), st.sampled_from((N + 1, -1))))
    x = np.array(vals, dtype=np.float64)
    for v in _layouts(x):
        _check(analysis.insert_row, v, n, circuit.numeric)
        _check(analysis.insert_row, v, True if n == 1 else n, circuit.numeric)


@pytest.mark.parametrize('N', [1, 2, 3, 7, 32, 33])
def test_every_row_and_special_value_of_a_matrix(N):
    for n in range(N):
        for k, v in enumerate(SPECIAL):
            A = np.arange(N * N, dtype=np.float64).reshape(N, N) - 3.5
            A.flat[(k * 7) % (N * N)] = v
            A.flat[-1 - (k % (N * N))] = SPECIAL[-1 - k]
            _check(analysis._reduce_ndarray, A, n)


@pytest.mark.parametrize('N', [1, 2, 7, 256, 257])
def test_every_row_and_special_value_of_a_vector(N):
    rows = sorted({0, 1, N // 2, N - 1, N} & set(range(N + 1)))
    for n in rows:
        for k, v in enumerate(SPECIAL):
            x = np.linspace(-1.0, 1.0, N)
            x[k % N] = v
            if n < N:
                _check(analysis._reduce_ndarray, x, n)
            _check(analysis.insert_row, x, n, circuit.numeric)


class _CountingNumpy:
    """`numpy`, its `empty` counted (the analysis module's)."""

    def __init__(self):
        self.empties = 0

    def __getattr__(self, name):
        f = getattr(np, name)
        if name != 'empty':
            return f

        def counted(*a, **k):
            self.empties += 1
            return f(*a, **k)
        return counted


def test_a_small_reduce_and_insert_make_no_slices(monkeypatch):
    """A 7 x 7 and a 7-vector reduced, a 6-vector inserted into: one
    indexed copy each, no `empty` of the slices' (the parent: three); a
    33 x 33 and a 257-vector past the cutoffs: the slices."""
    counting = _CountingNumpy()
    monkeypatch.setattr(analysis, 'numpy', counting)
    tk = circuit.numeric
    A = np.arange(49.0).reshape(7, 7)
    analysis._reduce_ndarray(A, 3)
    analysis._reduce_ndarray(A[0], 0)
    analysis.insert_row(np.ones(6), 6, tk)
    assert counting.empties == 0
    analysis._reduce_ndarray(np.ones((33, 33)), 1)
    analysis._reduce_ndarray(np.ones(257), 1)
    analysis.insert_row(np.ones(257), 1, tk)
    assert counting.empties == 3


def test_the_cached_indices_cannot_be_written():
    analysis._reduce_ndarray(np.ones((5, 5)), 2)
    analysis.insert_row(np.ones(4), 1, circuit.numeric)
    assert analysis._TAKE_IDX
    for idx in analysis._TAKE_IDX.values():
        assert not idx.flags.writeable
