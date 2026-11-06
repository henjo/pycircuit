"""THE FULL-RANK PROOF in the branch check's structural rank
(`_tran_branch._full_rank_bound`, `_unproved`; roadmap 266).  A C whose
capacitors join every node is one block of the rank's count, an SVD of
n x n at each of three states (2.5 s a state at 2 001 unknowns); a block of
`PROOF_MIN` or more is first offered a proof that every singular value
lies past the band's top -- its comparison matrix an M-matrix by two KLU
solves, so |B^-1| <= M^-1 and sigma_min >= 1 / sqrt(||M^-1||_1
||M^-1||_inf) -- counted whole where it holds.  The count is the SVD's:
the bound is a lower bound, and it must clear the top by the SVD's
rounding."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_branch, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit._tran_branch import (
    _full_rank_bound,
    _inverse_norm_bound,
    _rank_by_blocks,
)
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.transient import Transient

pytestmark = pytest.mark.skipif(_tran_branch._klu() is None, reason='libklu does not load')


def _entries(B):
    r, c = np.nonzero(B)
    return r, c, B[r, c]


def _h_block(rng, k, margin=0.2, density=0.06):
    """A nonsymmetric matrix whose |diagonal| exceeds its rows' and columns'
    other entries' sums by `margin` of them (an H-matrix), random signs."""
    N = rng.standard_normal((k, k)) * (rng.random((k, k)) < density)
    np.fill_diagonal(N, 0.0)
    s = np.maximum(np.abs(N).sum(axis=0), np.abs(N).sum(axis=1))
    d = (1.0 + margin) * s + rng.random(k) * 1e-3 + 1e-6
    return N + np.diag(d * rng.choice((-1.0, 1.0), k))


def test_the_bound_is_under_the_smallest_singular_value():
    """On H-matrices of many shapes and scales the bound is positive and no
    more than the SVD's sigma_min; a singular matrix, one missing a diagonal
    entry, one with a NaN, and one coupled more strongly than it is grounded
    get none."""
    rng = np.random.default_rng(20261106)
    for _ in range(60):
        k = int(rng.integers(2, 120))
        B = _h_block(rng, k, margin=float(rng.choice((1e-3, 0.05, 1.0)))) * 10.0 ** rng.integers(-16, 3)
        bound = _full_rank_bound(k, *_entries(B))
        smin = np.linalg.svd(B, compute_uv=False).min()
        assert 0.0 < bound <= smin * (1.0 + 1e-9)
    B = _h_block(rng, 40)
    S = B.copy()
    S[7] = S[3]
    assert _full_rank_bound(40, *_entries(S)) == 0.0
    D = B.copy()
    D[5, 5] = 0.0
    assert _full_rank_bound(40, *_entries(D)) == 0.0
    Nn = B.copy()
    Nn[2, 2] = np.nan
    assert _full_rank_bound(40, *_entries(Nn)) == 0.0
    ## (strongly coupled, weakly grounded: a chain of 2 x 2 rotations)
    W = np.eye(40) * 0.1
    for i in range(0, 40, 2):
        W[i, i + 1], W[i + 1, i] = 1.0, -1.0
    assert _full_rank_bound(40, *_entries(W)) == 0.0


def test_a_margin_within_the_sums_rounding_is_not_taken():
    """M x > 0 shown only past its sums' rounding: M = [[1, -(1 - 1e-15)],
    [-(1 - 1e-15), 1]] at x = (1, 1) has M x = 1.1e-15 each, inside the
    rounding of |M| x = 2 -- no bound, though the matrix is an M-matrix."""
    r, c = np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1])
    a = np.array([1.0, -(1.0 - 1e-15), -(1.0 - 1e-15), 1.0])
    assert _inverse_norm_bound(2, r, c, a, np.ones(2)) == np.inf
    a[[1, 2]] = -0.5
    assert _inverse_norm_bound(2, r, c, a, np.ones(2)) == pytest.approx(2.0, rel=1e-14)


def test_a_bound_within_the_svds_rounding_of_the_top_is_not_taken():
    """A connected block of 80 with one diagonal entry 1e-13 against the
    rest's 1: its bound (~1e-13) clears `tol` = 1e-14 but not the SVD's
    rounding of the top (8 k eps ||B||_F ~ 1.3e-12) -- the SVD decides,
    nothing proved."""
    k = 80
    B = np.eye(k)
    B[-1, -1] = 1e-13
    for i in range(k - 1):
        B[i, i + 1] = B[i + 1, i] = -1e-18
    assert _full_rank_bound(k, *_entries(B)) > 1e-14
    before = _paths.snapshot()
    assert _rank_by_blocks(B, 1e-14) == k
    assert _paths.since(before).get('branch.rank:proved', 0) == 0


def _blocks_matrix(rng, n):
    """Blocks of 1 to 140 under a permutation: H-matrix blocks, singular
    ones, non-dominant ones, scales across the tolerance."""
    C = np.zeros((n, n))
    perm = rng.permutation(n)
    i = 0
    while i < n:
        k = min(int(rng.choice((1, 3, 20, 64, 70, 90, 140))), n - i)
        kind = rng.integers(4)
        if kind == 0 or k < 3:
            B = _h_block(rng, k)
        elif kind == 1:
            B = _h_block(rng, k)
            B[-1] = B[0]
        elif kind == 2:
            B = rng.standard_normal((k, k)) * (rng.random((k, k)) < 0.2) + np.eye(k) * 0.01
        else:
            B = _h_block(rng, k, margin=0.01)
        idx = perm[i:i + k]
        C[np.ix_(idx, idx)] = B * 10.0 ** rng.integers(-14, -9)
        i += k
    return C


@pytest.mark.parametrize('seed', range(6))
def test_the_count_is_the_svds(seed, monkeypatch):
    """`_rank_by_blocks` with the proof gives what it gives without it
    (`PROOF_MIN` out of reach) -- its None included -- whole or through kept
    indices, its pattern given or scanned, under bands around the scales of
    the blocks; and the proof served."""
    rng = np.random.default_rng(seed)
    n = 400
    Cm = _blocks_matrix(rng, n)
    keep = np.r_[0:5, 6:n]
    red = Cm[np.ix_(keep, keep)]
    rnz = np.nonzero(red)
    cases = []
    for tol in (1e-25, 1e-15, 1e-13, 1e-11):
        for band in (None, (tol * 0.5, tol * 2.0), (tol * 0.5, tol * 1e3), (tol * 1e-3, tol * 1e5)):
            cases.append((tol, band))
    before = _paths.snapshot()
    got = [(_rank_by_blocks(Cm, t, band=b), _rank_by_blocks(Cm, t, band=b, nz=np.nonzero(Cm)),
            _rank_by_blocks(Cm, t, band=b, keep=keep), _rank_by_blocks(Cm, t, band=b, nz=rnz, keep=keep))
           for t, b in cases]
    proved = _paths.since(before).get('branch.rank:proved', 0)
    monkeypatch.setattr(_tran_branch, 'PROOF_MIN', 10 ** 9)
    want = [(_rank_by_blocks(Cm, t, band=b), _rank_by_blocks(red, t, band=b)) for t, b in cases]
    for (g0, g1, g2, g3), (w0, w1) in zip(got, want):
        assert g0 == g1 == w0
        assert g2 == g3 == w1
        assert all(g is None or type(g) is int for g in (g0, g1, g2, g3))
    assert proved > 0
    assert any(w is None for w, _ in want) and any(w is not None for w, _ in want)


def test_a_block_inside_the_band_goes_to_its_svd(monkeypatch):
    """An H-matrix block whose smallest singular value lies inside the band
    above `tol`: the proof's bound clears `tol` but not the band's top, so
    the SVD decides -- None, as without the proof."""
    rng = np.random.default_rng(7)
    k = 80
    B = _h_block(rng, k, margin=1.0)
    smin = np.linalg.svd(B, compute_uv=False).min()
    tol = smin / 50.0
    band = (tol * 0.5, tol * 100.0)
    assert _full_rank_bound(k, *_entries(B)) > tol
    assert _rank_by_blocks(B, tol, band=band) is None
    monkeypatch.setattr(_tran_branch, 'PROOF_MIN', 10 ** 9)
    assert _rank_by_blocks(B, tol, band=band) is None
    assert _rank_by_blocks(B, tol, band=(tol * 0.5, tol * 2.0)) == k


def _coupled(n):
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=2e-8, per=5e-8)
    for k in range(n):
        g = 'g0' if k == 0 else f'd{k - 1}'
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', g, gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
        c[f'cm{k}'] = C(f'd{k}', g, c=2e-15)
        if k >= 3 and k % 3 == 0:
            c[f'cl{k}'] = C(f'd{k}', f'd{k - 3}', c=1e-15)
    return c


def _rank0(monkeypatch, proof, klu=True):
    if not proof:
        monkeypatch.setattr(_tran_branch, 'PROOF_MIN', 10 ** 9)
    if not klu:
        monkeypatch.setitem(_tran_branch._KLU, 'lib', False)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_coupled(150), toolkit=circuit.numeric)
        tr.solve(tend=2e-9, timestep=1e-9)
    return tr.__dict__.get('_branch_rank0'), _paths.since(before).get('branch.rank:proved', 0)


def test_a_coupled_chains_structural_rank_is_the_svds(monkeypatch):
    """A MOS chain with a gate-drain capacitor at every stage (one block of
    its 151 nodes): the run's structural rank and scale with the proof --
    serving each of the three states -- are the SVD's; without libklu the
    SVD serves and gives the same."""
    got, proved = _rank0(monkeypatch, True)
    assert proved == 3 and got[0] > 150 and type(got[0]) is int
    with monkeypatch.context() as m:
        want, none = _rank0(m, False)
    assert got == want and none == 0
    with monkeypatch.context() as m:
        assert _rank0(m, True, klu=False) == (want, 0)
