"""THE PACKED LAYOUT (`_packed`): from `PACKED_FROM` unknowns the DC and the
multistep transient hold their matrices on the circuit's structural pattern.
Every operation the analyses make is the dense one's bytes per entry; the
products the decisions read are summed in one defined order (the CSC
kernel); a run on the pattern is the dense run's, bit for bit, with the same
solver; a numpy function not made on the pattern raises where it meets one;
below the size nothing changes; the C paths decline inside the scope."""
import importlib
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _packed, _paths, analysis, circuit, nrsolver
from pycircuit.circuit import linearsolver as _ls
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import VS, R, SubCircuit, gnd
from pycircuit.circuit.tests.test_first_attempt_c import _chain
from pycircuit.circuit.tests.test_klu_order import _layers
from pycircuit.circuit.transient import Transient

numeric = circuit.numeric
_C_MODULES = ('_tran_newton_c', '_tran_loop_c', '_tran_lte_c', '_dc_newton_c', '_tran_branch_c', '_pcnr_c',
              '_tran_ctrl_c', '_tran_pred_c')


def _drawn(rng, n, nnz):
    """A pattern with its whole diagonal and `nnz` other places, values over
    six decades, a signed zero among them: `(PatternMatrix, its dense twin)`."""
    r = np.r_[np.arange(n), rng.integers(0, n, nnz)]
    c = np.r_[np.arange(n), rng.integers(0, n, nnz)]
    key = np.unique(c * n + r)
    pat = _packed.Pattern(n, key % n, key // n)
    v = rng.choice([-1.0, 1.0], pat.nnz) * 10.0 ** rng.uniform(-3, 3, pat.nnz)
    v[rng.integers(0, pat.nnz, 3)] = -0.0
    D = np.zeros((n, n))
    D[pat.rows, pat.cols] = v
    return _packed.PatternMatrix(pat, v), D


def test_every_operation_is_the_dense_ones_bytes():
    """Sums, scales, abs, negation, the cut, a gather, the nonzeros, the
    largest entry, the diagonal, `+ eye(n) d` -- each the dense matrix's
    bytes at every entry of the pattern (and zero off it)."""
    rng = np.random.default_rng(1)
    A, Da = _drawn(rng, 40, 120)
    B = _packed.PatternMatrix(A.pat, rng.normal(size=A.pat.nnz))
    Db = B.to_dense()

    def same(P, D):
        ## (the pattern's places by their bytes; off it the dense matrix holds
        ## zeros -- a negation's -0.0 among them -- which the pattern has not)
        D = np.asarray(D, dtype=float)
        off = np.ones(D.shape, dtype=bool)
        off[P.pat.rows, P.pat.cols] = False
        return P.v.tobytes() == D[P.pat.rows, P.pat.cols].tobytes() and not D[off].any()
    assert same(A + B, Da + Db) and same(A - B, Da - Db) and same(-A, -Da)
    assert same(A * 3.5, Da * 3.5) and same(2.0 * A, 2.0 * Da) and same(A / 7.0, Da / 7.0)
    assert same(abs(A), np.abs(Da)) and same(np.abs(A), np.abs(Da))
    for k in (0, 17, 39):
        assert same(A.reduced(k), analysis.remove_row_col((Da,), k, numeric)[0])
    r, c = rng.integers(0, 40, 50), rng.integers(0, 40, 50)
    assert A[r, c].tobytes() == Da[r, c].tobytes()
    assert all(np.array_equal(a, b) for a, b in zip(A.nonzero(), np.nonzero(Da)))
    _rr, _cc, vv = A.values_rowmajor()
    assert vv.tobytes() == Da[np.nonzero(Da)].tobytes()
    assert A.max() == Da.max() and abs(A).max() == np.abs(Da).max()
    assert A.diagonal().tobytes() == Da.diagonal().tobytes()
    d = rng.uniform(0.0, 1.0, 40)
    for g in (0.25, d):
        assert same(A.plus_diagonal(g), nrsolver._plus_diagonal(Da, g, numeric))
    assert A.count_nonzero() == np.count_nonzero(Da)


def test_the_csc_kernel_is_the_defined_order():
    """`A @ x`: each product rounded, then added in the pattern's order --
    column by column, rows ascending -- the loop written out, bit for bit."""
    rng = np.random.default_rng(2)
    A, _D = _drawn(rng, 30, 200)
    x = rng.normal(size=30) * 10.0 ** rng.uniform(-8, 8, 30)
    y = np.zeros(30)
    for j in range(30):
        for p in range(A.pat.Ap[j], A.pat.Ap[j + 1]):
            y[A.pat.rows[p]] += A.v[p] * x[j]
    assert (A @ x).tobytes() == y.tobytes() and np.dot(A, x).tobytes() == y.tobytes()


def test_nothing_dense_is_made_unseen():
    """A numpy function or ufunc not made on the pattern raises where it
    meets one; no dense copy from `DENSE_MAX` unknowns."""
    A, _D = _drawn(np.random.default_rng(3), 20, 40)
    for f in (np.asarray, np.linalg.inv, lambda a: np.sum(a), lambda a: a > 0):
        with pytest.raises(TypeError):
            f(A)
    big = _packed.PatternMatrix(_packed.Pattern(_packed.DENSE_MAX, [0], [0]), [1.0])
    with pytest.raises(MemoryError):
        big.to_dense()


def test_the_plan_assembles_on_the_pattern(monkeypatch):
    """Inside the scope the circuit's G and C are PatternMatrix, each place
    the dense assembly's bytes; the DC's cut assembly too; below the size,
    and outside the scope, dense as before."""
    cir = _chain(20)
    x = np.random.default_rng(4).uniform(-1, 1, cir.n)
    dense = {m: np.asarray(getattr(cir, m)(x, defaultepar)) for m in ('G', 'C')}
    with _packed.scope(cir):
        assert type(cir.G(x, defaultepar)) is np.ndarray
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    with _packed.scope(cir):
        got = {m: getattr(cir, m)(x, defaultepar) for m in ('G', 'C')}
        from pycircuit.circuit import _stamp_plan
        cut = _stamp_plan.assemble_matrix(cir, 'G', x, (defaultepar,), drop=0)
    for m in ('G', 'C'):
        assert type(got[m]) is _packed.PatternMatrix
        assert got[m].to_dense().tobytes() == dense[m].tobytes()
    assert cut.to_dense().tobytes() == analysis.remove_row_col((dense['G'],), 0, numeric)[0].tobytes()
    assert type(cir.G(x, defaultepar)) is np.ndarray


def _same_solver(monkeypatch):
    """StructuralKLU whatever the size and fill, as the packed layout takes:
    the dense run on the same solver."""
    init = _ls.AutoSolver.__init__

    def ainit(self, *a, **kw):
        init(self, *a, **kw)
        self.max_fill = 1.0
    monkeypatch.setattr(_ls.AutoSolver, '__init__', ainit)
    monkeypatch.setattr(Transient, '_STRUCTURAL_KLU_MIN_N', property(lambda self: 1))
    monkeypatch.setattr(DC, '_STRUCTURAL_KLU_MIN_N', 1)


def _run(build, packed, monkeypatch, c_on=False, **kw):
    for name in _C_MODULES:
        monkeypatch.setattr(importlib.import_module('pycircuit.circuit.' + name), 'ENABLED', c_on)
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0 if packed else 1 << 30)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        dc = DC(build(), toolkit=numeric).solve()
        tr = Transient(build(), toolkit=numeric)
        res = tr.solve(**kw)
    st = tr.statistics
    return ((np.asarray(dc.x, float).tobytes(), np.asarray(res.x, float).tobytes(),
             np.asarray(res.sweep_values, float).tobytes(),
             (st.accepted_steps, st.rejected_steps, st.newton_iterations)), _paths.since(before))


@pytest.mark.parametrize('build, kw', [
    (lambda: _layers(8), {'tend': 3e-10, 'timestep': 1e-11}),
    (lambda: _layers(8), {'tend': 1e-10, 'timestep': 1e-11, 'fixed_timestep': True}),
    (lambda: _chain(30), {'tend': 4e-9, 'timestep': 1e-10}),
], ids=['grid-adaptive', 'grid-fixed', 'mos-chain'])
def test_a_packed_run_is_the_dense_one(build, kw, monkeypatch):
    """The DC and the Gear-2 transient on the pattern (the threshold at 0,
    the Python's path): the dense run's bytes and statistics on the same
    solver -- every decision the same way."""
    _same_solver(monkeypatch)
    dense, _c = _run(build, False, monkeypatch, **kw)
    packed, _c = _run(build, True, monkeypatch, **kw)
    assert packed == dense


@pytest.mark.parametrize('build, kw, keys', [
    (lambda: _layers(8), {'tend': 3e-10, 'timestep': 1e-11},
     ('dcnewton_c:packed', 'newton_c:packed', 'loop_c:packed', 'core:packed')),
    (lambda: _chain(30), {'tend': 4e-9, 'timestep': 1e-10}, ('pcnr_c:packed', 'loop_c:packed')),
], ids=['grid', 'mos-chain'])
def test_the_c_paths_decline_on_the_pattern(build, kw, keys, monkeypatch):
    """With the C on, inside the scope each C entry declines, counted --
    the DC's plain Newton and PCNR, the Newton, the loop, the core: the run
    the Python's on the pattern."""
    _same_solver(monkeypatch)
    py, _c = _run(build, True, monkeypatch, **kw)
    got, counts = _run(build, True, monkeypatch, c_on=True, **kw)
    assert got == py
    for key in keys:
        assert counts.get(key, 0) > 0, (key, counts)


def test_held_nodes_and_singular_rows_on_the_pattern(monkeypatch):
    """A held node only a voltage source fixes is named the same on the
    pattern; a matrix with an empty column is reported the same way."""
    def held():
        c = SubCircuit()
        c['v'] = VS('a', gnd, v=1.0)
        c['r'] = R('a', 'b', r=1e3)
        c['r2'] = R('b', gnd, r=1e3)
        return c
    dense = DC(held(), toolkit=numeric).unholdable({'a': 0.5})
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    cir = held()
    with _packed.scope(cir):
        got = DC(cir, toolkit=numeric).unholdable({'a': 0.5})
    assert got == dense and got
    A, _D = _drawn(np.random.default_rng(5), 12, 20)
    v = A.v.copy()
    v[A.pat.cols == 5] = 0.0
    E = _packed.PatternMatrix(A.pat, v)
    assert nrsolver._structural_singularity(E, None, numeric) == \
        nrsolver._structural_singularity(E.to_dense(), None, numeric)
