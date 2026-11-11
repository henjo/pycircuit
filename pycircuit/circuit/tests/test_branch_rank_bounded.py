"""THE STRUCTURAL RANK DECIDED ON THE NONZEROS (`_tran_branch._rank_bounded`):
the branch check's rounding margin (`np.linalg.norm(C)`) and KCL residual
(two n x m sums) known within bounds from `C`'s nonzeros alone, every
decision they feed taken for the bounds' whole range -- the dense passes
where a decision's edge lies inside it.  The rank the dense passes' on every
circuit; the bounds hold the dense values whatever the summation order; an
edge inside the bound, and a singular value between a band's narrowest and
widest ends, left to the dense passes."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_branch, circuit
from pycircuit.circuit._tran_branch import (
    UNDECIDED,
    _blocks_bounded,
    _norm_range,
    _rank_bounded,
    _sums_norm_range,
)
from pycircuit.circuit.elements import IS, VS, C, L, R, SubCircuit, gnd
from pycircuit.circuit.tests.test_first_attempt_c import _chain
from pycircuit.circuit.tests.test_klu_order import _layers
from pycircuit.circuit.transient import Transient


def _coupled(k=12):
    """A ladder with capacitors between nodes as well as to ground, and an
    inductor (a branch row C holds): blocks larger than one."""
    c = SubCircuit()
    c['vs'] = VS('n0', gnd, v=1.0)
    for i in range(k):
        c[f'r{i}'] = R(f'n{i}', f'n{i + 1}', r=1e3)
        c[f'cg{i}'] = C(f'n{i + 1}', gnd, c=1e-12 * (1 + i % 3))
        if i % 2:
            c[f'cc{i}'] = C(f'n{i}', f'n{i + 1}', c=3e-13)
    c['l'] = L(f'n{k}', 'x', L=1e-9)
    c['rl'] = R('x', gnd, r=50.0)
    c['is'] = IS(f'n{k // 2}', gnd, i=1e-3)
    return c


def _rank(build, bounded, monkeypatch, **kw):
    monkeypatch.setattr(_tran_branch, 'BOUNDED', bounded)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(build(), toolkit=circuit.numeric)
    for k, v in kw.items():
        setattr(tr, k, v)
    before = _paths.snapshot()
    got = tr._branch_structural_rank()
    return got, {k: v for k, v in _paths.since(before).items() if k.startswith('branch.rank')}


@pytest.mark.parametrize('build', [lambda: _layers(8), lambda: _layers(16), lambda: _chain(400),
                                   _coupled], ids=['layers8', 'layers16', 'chain400', 'coupled'])
@pytest.mark.parametrize('cert', [None, 0], ids=['cert', 'held'])
def test_the_bounded_rank_is_the_dense_one(build, cert, monkeypatch):
    """Each state decided on the nonzeros (`branch.rank:bounded`), the rank
    and scale the dense passes' -- with the certificate's generic C kept, and
    past its size (the states on one held array)."""
    kw = {} if cert is None else {'BRANCH_CERT_MAX_N': cert}
    dense, dc = _rank(build, False, monkeypatch, **kw)
    got, gc = _rank(build, True, monkeypatch, **kw)
    assert got == dense
    assert gc.get('branch.rank:bounded', 0) >= 1 and 'branch.rank:dense' not in gc, gc
    assert gc.get('branch.rank:proved', 0) == dc.get('branch.rank:proved', 0)


def _drawn(rng, n, nnz, spread):
    """A sparse n x n matrix of `nnz` entries over `spread` decades, half of
    every column's and row's mass cancelling: its dense form and nonzeros."""
    A = np.zeros((n, n))
    r, c = rng.integers(0, n, nnz), rng.integers(0, n, nnz)
    v = rng.choice([-1.0, 1.0], nnz) * 10.0 ** rng.uniform(-spread, spread, nnz)
    A[r, c] = v
    A[r, (c + 1) % n] = -v * (1.0 + rng.normal(0.0, 1e-9, nnz))
    rows, cols = np.nonzero(A)
    return A, rows, cols, A[rows, cols]


def test_the_bounds_hold_the_dense_values():
    """numpy's dense norm (BLAS's sum of n^2 squares) and its dense sums'
    norms lie inside the bounds the nonzeros give, on matrices whose sums
    cancel over many decades -- where the orders differ in their last bits."""
    rng = np.random.default_rng(11)
    moved = 0
    for trial in range(60):
        n = int(rng.integers(20, 160))
        A, rows, cols, vals = _drawn(rng, n, int(rng.integers(n, 6 * n)), float(rng.uniform(0, 12)))
        lo, hi = _norm_range(vals)
        d = float(np.linalg.norm(A))
        assert lo <= d <= hi, (trial, lo, d, hi)
        m = int(rng.integers(1, n + 1))
        for idx, sel, dense in ((cols, rows < m, A[:m].sum(axis=0)), (rows, cols < m, A[:, :m].sum(axis=1))):
            lo, hi = _sums_norm_range(idx[sel], vals[sel], n)
            d = float(np.linalg.norm(dense))
            assert lo <= d <= hi, (trial, lo, d, hi)
            moved += float(np.linalg.norm(np.bincount(idx[sel], weights=vals[sel], minlength=n))) != d
    ## (the fixtures reach the rounding: the bincount's order and numpy's
    ## parted on some, which a zero-width bound would have failed)
    assert moved > 0


def test_an_edge_inside_the_bound_is_left_to_the_dense_passes():
    """`delta + eps < tol` with `tol` inside the bounds' range: undecided --
    and the same `tol` well clear of it decided; a singular value between a
    band's widest and narrowest ends undecided, inside the narrowest None."""
    rng = np.random.default_rng(5)
    A, rows, cols, vals = _drawn(rng, 60, 240, 6.0)
    n, m, g = 60, 40, 0
    nr = _norm_range(vals)
    c = 8 * n * np.finfo(float).eps
    a = _sums_norm_range(cols[rows < m], vals[rows < m], n)
    b = _sums_norm_range(rows[cols < m], vals[cols < m], n)
    lo, hi = a[0] + b[0] + c * nr[0], a[1] + b[1] + c * nr[1]
    assert lo < hi
    assert _rank_bounded(A, 0.5 * (lo + hi), (rows, cols, vals), g, m) is UNDECIDED
    assert _rank_bounded(A, 1e3 * hi, (rows, cols, vals), g, m) is not UNDECIDED
    D = np.diag([1.0, 2.0, 3.0, 4.0])
    nz = np.nonzero(D)
    ## (2 inside the widest band (1.5, 2.7], the narrowest (2.5, 2.6] empty)
    assert _blocks_bounded(D, 2.5, (1.5, 2.5), (2.6, 2.7), nz) is UNDECIDED
    assert _blocks_bounded(D, 2.5, (0.5, 1.5), (3.5, 3.6), nz) is None
    assert _blocks_bounded(D, 2.5, (2.1, 2.2), (2.6, 2.7), nz) == 2


def test_an_undecided_state_is_the_dense_passes(monkeypatch):
    """Every state's edge inside the bound (`_rank_bounded` undecided): the
    dense passes decide, counted, the same rank."""
    dense, _ = _rank(lambda: _layers(8), False, monkeypatch)
    monkeypatch.setattr(_tran_branch, '_rank_bounded', lambda *a: UNDECIDED)
    got, gc = _rank(lambda: _layers(8), True, monkeypatch)
    assert got == dense and gc.get('branch.rank:dense', 0) == 3, gc
