"""The branch check's structural rank on C's pattern (`_tran_branch.
_plan_nonzeros`, `_rank_by_blocks`'s `nz` and `keep`): C's nonzeros
gathered at the stamp plan's positions -- proven all of them by one count
-- give max |C| and the blocks, and the reduced matrix is read through the
kept indices, where the check scanned and copied n x n matrices three times
a run.  The answer is the dense path's: the same rank, scale and blocks."""
import warnings
from types import SimpleNamespace

import numpy as np

from pycircuit.circuit import _tran_branch, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit._tran_branch import _BranchCheck, _plan_nonzeros, _rank_by_blocks
from pycircuit.circuit.elements import VS, C, L, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.transient import Transient


def _blocks(rng, n):
    """A matrix of random blocks (sizes 1 to 4) under a permutation, some
    entries either side of a tolerance, some blocks rank-deficient, signed
    zeros among the zeros."""
    C = np.zeros((n, n))
    perm = rng.permutation(n)
    i = 0
    while i < n:
        k = min(int(rng.integers(1, 5)), n - i)
        B = rng.standard_normal((k, k)) * 10.0 ** rng.integers(-12, 0)
        if k > 1 and rng.random() < 0.3:
            B[-1] = B[0]
        idx = perm[i:i + k]
        C[np.ix_(idx, idx)] = B
        i += k
    z = rng.random((n, n)) < 0.05
    C[z & (C == 0.0)] = -0.0
    return C


def _checker(g, m):
    return SimpleNamespace(irefnode=g, cir=SimpleNamespace(nodes=[None] * m))


def test_the_blocks_read_through_kept_indices_are_the_copys():
    """`_rank_by_blocks` on `C` through `keep` -- with its pattern given or
    scanned -- is its count on the copy `C[ix_(keep, keep)]`, the band's
    None included; on the whole `C` with the pattern given, its count
    scanned."""
    rng = np.random.default_rng(20261105)
    nones = 0
    for _ in range(200):
        n = int(rng.integers(2, 40))
        Cm = _blocks(rng, n)
        g = int(rng.integers(0, n))
        keep = np.r_[0:g, g + 1:n]
        red = Cm[np.ix_(keep, keep)]
        rows, cols = np.nonzero(Cm)
        on = (rows != g) & (cols != g)
        rnz = (rows[on] - (rows[on] > g), cols[on] - (cols[on] > g))
        assert all(np.array_equal(a, b) for a, b in zip(rnz, np.nonzero(red)))
        sv = np.linalg.svd(Cm, compute_uv=False)
        for tol in (1e-9, sv[sv > 0][-1] if np.any(sv > 0) else 1.0, 1e-30):
            for band in (None, (tol * 0.5, tol * 2.0)):
                want = _rank_by_blocks(red, tol, band=band)
                nones += want is None
                assert _rank_by_blocks(Cm, tol, band=band, nz=rnz, keep=keep) == want
                assert _rank_by_blocks(Cm, tol, band=band, keep=keep) == want
                assert (_rank_by_blocks(Cm, tol, band=band, nz=(rows, cols))
                        == _rank_by_blocks(Cm, tol, band=band))
    assert nones > 20


def test_the_reduced_rank_from_the_pattern_is_the_dense_ones():
    """`_reduced_rank` with C's nonzeros gives what it gives without them,
    on reference-embedded matrices of every size and reference."""
    rng = np.random.default_rng(20261106)
    for _ in range(200):
        n = int(rng.integers(2, 30))
        m = int(rng.integers(1, n + 1))
        g = int(rng.integers(0, m))
        Cm = _blocks(rng, n)
        Cm[g, :] = 0.0
        Cm[:, g] = 0.0
        Cm[g, :m] = -Cm[:m, :m].sum(axis=0)
        Cm[:m, g] = -Cm[:m, :m].sum(axis=1)
        rows, cols = np.nonzero(Cm)
        nz = (rows, cols, Cm[rows, cols])
        for tol in (1e-9, 1e-14, 0.3):
            eps = 8 * n * np.finfo(float).eps * float(np.linalg.norm(Cm))
            assert (_BranchCheck._reduced_rank(_checker(g, m), Cm, tol, eps, nz)
                    == _BranchCheck._reduced_rank(_checker(g, m), Cm, tol, eps))


def test_the_plans_nonzeros_are_all_of_them_or_none():
    """The plan's positions (repeats among them) gathered: the nonzeros in
    `np.nonzero`'s order, a NaN one of them, a signed zero not; a nonzero
    off the positions, positions out of range, a matrix not contiguous, or
    no plan: None."""
    Cm = np.zeros((4, 4))
    Cm[0, 0], Cm[1, 3], Cm[3, 1], Cm[2, 2] = 1.0, -2.0, np.nan, -0.0
    flat = [5, 0, 7, 13, 10, 0, 15]
    cir = SimpleNamespace(_stamp_plan=SimpleNamespace(
        methods={'C': SimpleNamespace(flat=np.asarray(flat, dtype=np.intp))}))
    rows, cols, vals = _plan_nonzeros(cir, Cm)
    r0, c0 = np.nonzero(Cm)
    assert np.array_equal(rows, r0) and np.array_equal(cols, c0)
    assert np.array_equal(vals, Cm[r0, c0], equal_nan=True)
    off = Cm.copy()
    off[3, 0] = 1e-300
    assert _plan_nonzeros(cir, off) is None
    cir._stamp_plan.methods['C'].flat = np.asarray([0, 7, 13, 16], dtype=np.intp)
    assert _plan_nonzeros(cir, Cm) is None
    cir._stamp_plan.methods['C'].flat = np.asarray(flat, dtype=np.intp)
    assert _plan_nonzeros(cir, np.asfortranarray(Cm + np.eye(4))) is None
    assert _plan_nonzeros(SimpleNamespace(), Cm) is None


def _chain(n=200):
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=4e-8, per=1e-7)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
        if k % 4 == 0:
            c[f'cc{k}'] = C(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', c=1e-15)
    c.add_node('l1')
    c['L1'] = L('d0', 'l1', L=1e-9)
    c['RL'] = R('l1', gnd, r=50.0)
    return c


def _rank0(monkeypatch, dense):
    """The structural rank and scale a run's check makes, the plan's
    nonzeros used or (`dense`) refused; and whether they served."""
    served = []
    real = _tran_branch._plan_nonzeros

    def spy(cir, Cm):
        got = None if dense else real(cir, Cm)
        served.append(got is not None)
        return got
    monkeypatch.setattr(_tran_branch, '_plan_nonzeros', spy)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(), toolkit=circuit.numeric)
        tr.solve(tend=4e-9, timestep=1e-9)
    return tr.__dict__.get('_branch_rank0'), served


def test_a_runs_structural_rank_is_the_dense_paths(monkeypatch):
    """A chain with coupling capacitors and an inductor (a branch current
    in C): the run's structural rank and scale on the plan's nonzeros are
    the dense path's, the plan serving every state."""
    got, served = _rank0(monkeypatch, dense=False)
    want, _s = _rank0(monkeypatch, dense=True)
    assert got == want and got[0] > 0
    assert served and all(served)
