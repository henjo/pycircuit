"""THE FROZEN FACTORS AT SIZE (`_RadauStages._radau_frozen_packed`; plan
261's L6b): where the transient's core is paged, the Radau transform's real
and complex factors are made at the analysis solver's structure alone --
or, before its choice, the core's pattern reduced -- from `C` and `G` at
those entries: the real one on pattern pages, the complex one's record
`prepare_structural`'s (`prepare_structural_values`) or `prepare`'s own
nonzeros (`prepare_packed`, `_csc_of_packed`).  No n x n `Cr`, `Gr` or
factor is made, and every piece is the dense path's, bit for bit."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _tran_core, circuit, linearsolver
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_radau_pattern_pages import _chain
from pycircuit.circuit.transient import Transient


def _b(a):
    return np.asarray(a).tobytes()


def _rec_bytes(rec):
    """A complex record's arrays and key (its product a partial of them)."""
    return (rec[1], _b(rec[2]), _b(rec[3]), _b(rec[4]), rec[5])


def _sparse_complex(m, rng, fill=0.06):
    """An m x m complex128 matrix: a random pattern, a few explicit (signed)
    zeros and a NaN among its entries -- and a superset of its pattern,
    column by column, rows ascending."""
    A = np.zeros((m, m), dtype=np.complex128)
    sup = rng.random((m, m)) < fill
    sup[np.arange(m), np.arange(m)] = True
    r, c = np.nonzero(sup)
    A[r, c] = rng.standard_normal(len(r)) + 1j * rng.standard_normal(len(r))
    z = rng.choice(len(r), 6, replace=False)
    A[r[z[:3]], c[z[:3]]] = 0.0
    A[r[z[3]], c[z[3]]] = complex(-0.0, -0.0)
    A[r[z[4]], c[z[4]]] = complex(np.nan, 0.0)
    order = np.lexsort((r, c))
    return A, r[order], c[order]


@pytest.mark.parametrize('seed', range(4))
def test_the_packed_record_is_the_dense_records(seed):
    """`_csc_of_packed` of a matrix's entries at a superset of its pattern is
    `_csc_of_dense`'s record of it: the nonzeros (a NaN one, a signed zero
    not), their `Ap`, `Ai`, key and mask -- and with the last pattern its
    arrays where the pattern is the same, new ones where it is not."""
    mv = linearsolver._csc_matvec()
    if mv is None:
        pytest.skip('no SciPy product')
    rng = np.random.default_rng(seed)
    m = 40
    A, rows, cols = _sparse_complex(m, rng)
    want, wlast = linearsolver._csc_of_dense(A, mv)
    got, glast = linearsolver._csc_of_packed(A[rows, cols], m, rows, cols, mv)
    assert _rec_bytes(got) == _rec_bytes(want)
    assert np.array_equal(glast[0], wlast[0]) and glast[3] == wlast[3]
    ## (the same pattern again: the last one's arrays; another: its own)
    again, alast = linearsolver._csc_of_packed(A[rows, cols] * 2.0, m, rows, cols, mv, glast)
    assert again[2] is glast[1] and again[3] is glast[2]
    assert all(a is b for a, b in zip(alast, glast))
    B = A.copy()
    B[rows[0], cols[0]] = 0.0 if B[rows[0], cols[0]] != 0 else 1.0
    wb, _ = linearsolver._csc_of_dense(B, mv, wlast)
    gb, gblast = linearsolver._csc_of_packed(B[rows, cols], m, rows, cols, mv, glast)
    assert _rec_bytes(gb) == _rec_bytes(wb) and gblast[1] is not glast[1]


def test_the_structural_values_record_is_the_structural_one():
    """`prepare_structural_values` on a matrix's entries at the structure is
    `prepare_structural`'s record of it."""
    try:
        zs = linearsolver.ComplexKLUSolver()
    except ImportError:
        pytest.skip('no complex KLU')
    rng = np.random.default_rng(5)
    m = 30
    A, rows, cols = _sparse_complex(m, rng)
    A[np.isnan(A)] = 1.0
    sk = linearsolver.StructuralKLU(m, rows, cols, fallback=lambda: None)
    want = zs.prepare_structural(A, sk)
    if want is None:
        pytest.skip('prepare_structural declines here')
    got = zs.prepare_structural_values(np.ascontiguousarray(A.ravel()[sk.take]), sk)
    assert _rec_bytes(got) == _rec_bytes(want)


@pytest.mark.parametrize('h', [1e-9, 3.7e-12, 2.5e-7])
def test_the_packed_factor_entries_are_numpys(h):
    """The factors' entries at the positions by numpy's operations on those
    entries alone are the dense factors' there -- the complex product's form
    included."""
    rng = np.random.default_rng(int(h * 1e15) % 2**31)
    m = 300
    Cr = np.zeros((m, m))
    Gr = np.zeros((m, m))
    take = np.flatnonzero(rng.random(m * m) < 0.02)
    Cr.reshape(-1)[take] = rng.standard_normal(len(take)) * 10.0 ** rng.integers(-18, -10, len(take))
    Gr.reshape(-1)[take] = rng.standard_normal(len(take)) * 10.0 ** rng.integers(-9, -2, len(take))
    Cr.reshape(-1)[take[:5]] = -0.0
    lam = np.array([3.6378342527444957 + 0j, 2.6810828736277521 + 3.0504301992474105j])
    s0, c1 = lam[0].real / h, lam[1] / h
    rf, cf = s0 * Cr + Gr, c1 * Cr + Gr
    Cp, Gp = Cr.reshape(-1)[take], Gr.reshape(-1)[take]
    assert _b(s0 * Cp + Gp) == _b(rf.reshape(-1)[take])
    assert _b(c1 * Cp + Gp) == _b(cf.reshape(-1)[take])


def _state(monkeypatch, n=300):
    """A transient run on pattern pages -- the analysis solver on its
    StructuralKLU (305 unknowns) -- and its state: x_n's C and G."""
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(n), toolkit=circuit.numeric, integrator=RadauIIA3Integrator())
        tr.solve(tend=3e-9, timestep=1e-9)
    x = np.asarray(tr._q_cache[0], dtype=float).copy()
    P = _tran_core.passes(tr, x, 'CG')
    return tr, P['C'], P['G']


def _frozen_pieces(fz):
    return (fz.m, fz.p, fz.v, _b(fz.rf), fz.rf.shape, _rec_bytes(fz.prep), fz.lu, fz.kept)


@pytest.mark.parametrize('chosen', [True, False], ids=['structure', 'unchosen'])
def test_the_packed_frozen_transform_is_the_dense_one(monkeypatch, chosen):
    """`_radau_frozen_packed` -- on the solver's StructuralKLU, and before
    its choice -- is `_radau_frozen` on the dense reduced `C` and `G`: P's
    and V's entries, the real factor's bytes, the complex record (the
    structural one, or its own nonzeros), the state each leaves."""
    tr, Cn, Gn = _state(monkeypatch)
    from pycircuit.circuit.analysis import remove_row_col
    ls = tr._get_linearsolver()
    if not chosen:
        monkeypatch.setattr(ls, '_choice', None)
    iref = tr._coupled_stage_context(np.asarray(tr._q_cache[0]), 3e-9, None).iref
    h = 1.25e-9
    pk = tr._radau_packing(iref, Cn, Gn)
    assert pk is not None and pk.ok and (pk.sk is None) == (not chosen)
    tr.__dict__.pop('_radau_Pc', None)
    zs = tr.__dict__.get('_radau_zsolver')
    last = None if zs in (None, 'unset') else zs._csc_last
    got = tr._radau_frozen_packed(pk, Cn, Gn, h)
    glast = None if zs in (None, 'unset') else zs._csc_last
    tr.__dict__.pop('_radau_Pc', None)
    if zs not in (None, 'unset'):
        zs._csc_last = last
    (Cr,) = remove_row_col((Cn,), iref, tr.toolkit)
    (Gr,) = remove_row_col((Gn,), iref, tr.toolkit)
    want = tr._radau_frozen(np.asarray(Cr, dtype=float), np.asarray(Gr, dtype=float), h)
    assert got is not None and want is not None
    assert _frozen_pieces(got) == _frozen_pieces(want)
    if not chosen and zs not in (None, 'unset'):
        assert glast[3] == zs._csc_last[3]


@pytest.mark.parametrize('c', [2.5e-10, 0.0, -0.0, -1.5, np.inf])
def test_a_combined_matrix_on_pattern_pages_is_numpys(c):
    """`combined_on_pattern`: `A + c * B`'s bytes -- on pattern pages where
    numpy's sum is zero off the pattern as they are, else numpy's."""
    n = 50
    rng = np.random.default_rng(9)
    pat = np.flatnonzero(rng.random(n * n) < 0.07)
    A = np.zeros(n * n)
    B = np.zeros(n * n)
    A[pat] = rng.standard_normal(len(pat))
    B[pat] = rng.standard_normal(len(pat))
    A[pat[:2]] = -0.0
    B[pat[2:4]] = -0.0
    A, B = A.reshape(n, n), B.reshape(n, n)
    with np.errstate(all='ignore'):
        want = A + c * B
        got = _tran_core.combined_on_pattern(A, c, B, pat, n)
    assert _b(got) == _b(want) and got.shape == want.shape


def test_the_packing_holds_every_entry_of_the_cores_reduced_pattern(monkeypatch):
    """The packing's positions -- the core's pattern reduced, column by
    column, rows ascending -- hold every reduced entry of the core's pattern;
    on the StructuralKLU, its structure's order (`take`)."""
    tr, Cn, Gn = _state(monkeypatch)
    iref = tr._coupled_stage_context(np.asarray(tr._q_cache[0]), 3e-9, None).iref
    pk = tr._radau_packing(iref, Cn, Gn)
    n = Cn.shape[0]
    full = np.zeros((n, n), dtype=bool)
    full.reshape(-1)[pk.ftake] = True
    nz = (Cn != 0) | (Gn != 0)
    nz[iref, :] = nz[:, iref] = False
    assert not np.any(nz & ~full)
    if pk.sk is not None:
        assert np.array_equal(pk.take, np.asarray(pk.sk.take, dtype=np.int64))
    else:
        assert np.all(np.diff(pk.cols * pk.m + pk.rows) > 0)
    assert pk.ok
