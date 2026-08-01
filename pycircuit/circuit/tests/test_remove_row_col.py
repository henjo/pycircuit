"""Tests for analysis.remove_row_col's blocked reduction (stage 7a)."""
import numpy as np
import pytest

from pycircuit.circuit import numeric
from pycircuit.circuit.analysis import remove_row_col


def _two_delete_reference(matrices, n, toolkit):
    """The form remove_row_col had before 7a: one toolkit.delete per axis.

    Kept as an independent reference rather than importing anything, so the fast
    path is checked against the semantics it replaced and not against itself.
    """
    result = []
    for A in matrices:
        for axis in range(len(A.shape)):
            A = toolkit.delete(A, [n], axis=axis)
        result.append(A)
    return tuple(result)


@pytest.mark.parametrize('n', [5, 17, 64])
@pytest.mark.parametrize('dtype', [float, complex])
def test_blocked_reduction_matches_the_two_delete_form(n, dtype):
    """Bit-identical for every removable index, both dtypes.

    Every index is swept, not just a middle one: the blocked copy writes four
    sub-blocks and two of them are EMPTY when the removed index is first or last,
    which is exactly where an off-by-one survives a spot check.
    """
    rng = np.random.default_rng(20260801)
    J = rng.random((n, n)).astype(dtype)
    f = rng.random(n).astype(dtype)
    if dtype is complex:
        J = J + 1j * rng.random((n, n))
        f = f + 1j * rng.random(n)

    for k in range(n):
        got_f, got_J = remove_row_col((f, J), k, numeric)
        exp_f, exp_J = _two_delete_reference((f, J), k, numeric)
        assert np.array_equal(got_J, exp_J), 'matrix differs at index %d' % k
        assert np.array_equal(got_f, exp_f), 'vector differs at index %d' % k
        assert got_J.dtype == exp_J.dtype


def test_non_square_and_higher_rank_fall_back():
    """The fast path is guarded; anything else keeps the general behaviour."""
    rng = np.random.default_rng(7)
    M = rng.random((6, 9))
    got, = remove_row_col((M,), 2, numeric)
    exp, = _two_delete_reference((M,), 2, numeric)
    assert np.array_equal(got, exp)
    assert got.shape == (5, 8)


def test_symbolic_toolkit_keeps_the_generic_path():
    """A sympy matrix is not a numpy float array and must not take the fast path.

    The guard is `isinstance(A, numpy.ndarray) and A.dtype != object`; a symbolic
    toolkit's arrays fail it, which is the point -- slicing does not mean the same
    thing there.
    """
    sympy = pytest.importorskip('sympy')
    from pycircuit.circuit import symbolic
    a, b, c, d = sympy.symbols('a b c d')
    M = np.array([[a, b, 0], [c, d, 0], [0, 0, 1]], dtype=object)
    got, = remove_row_col((M,), 2, symbolic)
    assert got.shape == (2, 2)
    assert got[0, 0] == a and got[1, 1] == d


def test_result_is_not_aliased_between_calls():
    """Two reductions must not share a buffer.

    A buffer cached across calls would be faster and is not safe: the transient
    keeps `J` for the step controller's LTE solve AFTER the Newton step has used
    it, so a second call would overwrite the first caller's matrix.  This is the
    test that would fail if someone adds that cache.
    """
    rng = np.random.default_rng(11)
    A = rng.random((8, 8))
    B = rng.random((8, 8))
    first, = remove_row_col((A,), 3, numeric)
    snapshot = first.copy()
    second, = remove_row_col((B,), 3, numeric)
    assert second is not first
    assert np.array_equal(first, snapshot), \
        'the first result changed when a second reduction ran -- buffers are shared'
