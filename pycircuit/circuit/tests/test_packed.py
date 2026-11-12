"""THE PACKED LAYOUT (`_packed`): from `PACKED_FROM` unknowns the DC and the
multistep transient hold their matrices on the circuit's structural pattern.
Every operation the analyses make is the dense one's bytes per entry; the
products the decisions read are summed in one defined order (the CSC
kernel); a run on the pattern is the dense run's, bit for bit, with the same
solver; a numpy function not made on the pattern raises where it meets one;
below the size nothing changes.  The DC Newton's C serves on the pattern --
the packed Python's bytes, its record and core holding no dense matrix;
the other C paths decline inside the scope."""
import importlib
import json
import warnings

import numpy as np
import pytest

from pycircuit.circuit import (
    _dc_newton_c,
    _packed,
    _paths,
    _tran_newton_c,
    analysis,
    circuit,
    nrsolver,
)
from pycircuit.circuit import linearsolver as _ls
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import VS, C, R, SubCircuit, gnd
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
    largest entry, the diagonal, `+ eye(n) d`, the kept rows and columns
    (`np.ix_`) -- each the dense matrix's
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
    keep = np.flatnonzero(rng.uniform(size=40) < 0.7)
    assert same(A.kept(keep), Da[np.ix_(keep, keep)])


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
     ('dcnewton_c:served', 'newton_c:served', 'lte_c:served', 'loop_c:steps')),
    (lambda: _chain(30), {'tend': 4e-9, 'timestep': 1e-10},
     ('pcnr_c:packed', 'newton_c:served', 'lte_c:served', 'loop_c:steps')),
], ids=['grid', 'mos-chain'])
def test_the_c_paths_on_the_pattern(build, kw, keys, monkeypatch):
    """With the C on, inside the scope the DC Newton's C, the transient
    Newton's, the error test's and the stepping loop's serve, PCNR's
    declines, counted: the run the packed Python's, bit for bit."""
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


def test_junction_gmin_on_the_pattern_is_the_dense_loops():
    """`nrsolver._junction_gmin` on a PatternMatrix: F and J the dense
    loop's bytes, a junction pair given twice summing in the loop's order;
    a junction entry off the pattern raises."""
    rng = np.random.default_rng(6)
    n = 12
    rows = [(1, 4), (4, 7), (1, 4), (9, 2)]
    r = np.r_[np.arange(n), [a for a, b in rows] * 2, [b for a, b in rows] * 2]
    c = np.r_[np.arange(n), [a for a, b in rows] + [b for a, b in rows],
              [b for a, b in rows] + [a for a, b in rows]]
    extra = rng.integers(0, n, (2, 15))
    key = np.unique(np.r_[c, extra[1]] * n + np.r_[r, extra[0]])
    pat = _packed.Pattern(n, key % n, key // n)
    v = rng.normal(size=pat.nnz) * 10.0 ** rng.uniform(-3, 3, pat.nnz)
    A = _packed.PatternMatrix(pat, v.copy())
    D = A.to_dense()
    x = rng.normal(size=n)
    F1, F2 = rng.normal(size=n), None
    F2 = F1.copy()
    g = 1e-3
    nrsolver._junction_gmin(F1, D, x, g, rows, nrsolver._junction_tables(rows))
    nrsolver._junction_gmin(F2, A, x, g, rows, nrsolver._junction_tables(rows))
    assert F2.tobytes() == F1.tobytes()
    assert A.v.tobytes() == D[pat.rows, pat.cols].tobytes()
    off = _packed.PatternMatrix(_packed.Pattern(n, np.arange(n), np.arange(n)), np.ones(n))
    with pytest.raises(IndexError):
        nrsolver._junction_gmin(F2, off, x, g, rows, nrsolver._junction_tables(rows))


@pytest.mark.parametrize('name', ['mike2', 'gm3', 'gm19'])
def test_the_dc_newtons_c_on_the_pattern_is_the_pythons(name, monkeypatch, ordinary_dc):
    """The decks' operating points on the pattern (the threshold at 0): every
    base solve -- junction gmin, gshunt and source stepping over the base;
    refactors, fallbacks, failing and singular solves -- its outcome, the
    counters, x0 and the warnings the C's as the packed Python's, the C
    serving; gm3 and gm19 reach junction gmin on the pattern (mike2's plain
    solve converges there, on KLU)."""
    from pycircuit.circuit.tests import test_dc_newton_c as T
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    seen = []
    real = nrsolver._junction_gmin_packed

    def counted(*a):
        seen.append(1)
        return real(*a)
    monkeypatch.setattr(nrsolver, '_junction_gmin_packed', counted)
    py, _c = T._record(name, monkeypatch, False)
    jg = len(seen)
    got, counts = T._record(name, monkeypatch, True)
    assert got == py
    assert counts.get('dcnewton_c:served', 0) > 0, counts
    assert (jg > 0) is (name != 'mike2'), jg


def test_the_rounding_floor_on_the_pattern_is_the_dense_ones_and_the_cs(monkeypatch):
    """The rounding floor on the pattern (`_within_rounding` by KLU's
    transposed solve; dcnewton.c's `floor_holds` with the CSC kernel's
    scale): the packed Python's operating points the dense Python's bytes
    on the same solver, the floor deciding them (neither converges without
    it); the C's the packed Python's, served."""
    from pycircuit.circuit.tests import test_dc_newton_c as T
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    dense, _c, kind = T._floor_dc(130, False)
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    seen, orig = [], nrsolver._within_rounding

    def within(*a):
        seen.append(orig(*a))
        return seen[-1]
    py, _c, _k = T._floor_dc(130, False, within)
    got, counts, _k = T._floor_dc(130, True)
    assert kind == 'StructuralKLU' and all(type(x) is bytes for x in py)
    assert py == dense and got == py
    assert sum(seen) >= 2 and counts.get('dcnewton_c:served', 0) == 1, (seen, counts)
    assert T._floor_dc(130, False, lambda *a: False)[0] == ['NoConvergenceError'] * 2


def test_the_iteration_limit_on_the_pattern_is_the_pythons(monkeypatch):
    """An operating point whose Newton meets its iteration limit on the
    pattern: the C's `Is` for the exception's scale the CSC kernel's, every
    base solve -- the limited ones' exceptions and their text -- the packed
    Python's."""
    from pycircuit.circuit.tests import test_dc_newton_c as T
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    py, _c = T._base_solves(lambda: T._ladder(160), False, maxiter=1)
    got, counts = T._base_solves(lambda: T._ladder(160), True, maxiter=1)
    assert got == py and any(e[0] == 'raise' for e in py), py
    assert counts.get('dcnewton_c:maxiter', 0) >= 1, counts


def test_the_dc_record_on_the_pattern_holds_no_dense_matrix(monkeypatch):
    """The DC Newton's record on the pattern: its core has no n x n G bin
    and is called for its packed sets; none of the record's buffers is
    m x m; its J at KLU's entries -- and below the size the dense record as
    before."""
    from pycircuit.circuit import _tran_core
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    _same_solver(monkeypatch)
    sets = _dc_newton_c._sources()[4]['PCX_WANT_SETS']

    def record(packed):
        monkeypatch.setattr(_packed, 'PACKED_FROM', 0 if packed else 1 << 30)
        dc = DC(_layers(20), toolkit=numeric)
        before = _paths.snapshot()
        dc.solve()
        assert _paths.since(before).get('dcnewton_c:served', 0) > 0
        return dc.__dict__['_dc_newton_c'], dc.cir.n - 1
    rec, m = record(True)
    assert rec.packed and rec.core.packed and rec.core.arrays['binG'].size == 1
    assert rec.s.pattern & sets and rec.s.packed == 1
    assert max(a.size for a in rec.buf.values()) < m * m
    assert len(rec.pk[0]) == rec.s.nact == rec.klu.nnz
    dense, m = record(False)
    assert not dense.packed and dense.buf['R'].size == m * m and not dense.s.pattern & sets
    assert type(_tran_core._Core).__name__ == 'type'


def test_the_refactor_residual_on_the_pattern_is_the_pythons(monkeypatch, ordinary_dc):
    """The refactor's residual test on the pattern (`StructuralKLU.
    _packed_residual`; pcx_kluseq.h's `pcx_klu_resid_packed`, the CSC kernel
    on KLU's arrays): gm19's operating point at a tolerance near its
    refactors' residuals -- most holding, some failing to a fresh factor --
    every base solve and the counters the C's as the packed Python's, the C
    serving."""
    from pycircuit.circuit import linearsolver
    from pycircuit.circuit.tests import test_dc_newton_c as T
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    monkeypatch.setattr(linearsolver.KLUSolver, 'REFACTOR_RESIDUAL_TOL', 1e-14)
    py, _c = T._record('gm19', monkeypatch, False)
    got, counts = T._record('gm19', monkeypatch, True)
    assert got == py and counts.get('dcnewton_c:served', 0) > 0, counts
    klu = json.loads(py)['klu']
    assert any(k[3] > 0 for k in klu) and any(k[0] > 2 * k[3] for k in klu), klu



def test_a_held_solve_runs_on_the_pattern(monkeypatch, ordinary_dc):
    """A nodeset's held solve (`dcanalysis._Held`) on the pattern: the kept
    rows and columns gathered there (`PatternMatrix.kept`), KLU on the held
    matrix's own pattern -- the structure describes the unheld one, whose
    choice the held solve leaves to the unheld solve -- each state of a latch
    picked, the dense run's operating point to its tolerance."""
    from pycircuit.circuit import dcanalysis, linearsolver
    from pycircuit.circuit.tests.test_dc_held import _latch
    seen = []
    real = dcanalysis._Held.system

    def system(self, F, J):
        seen.append(type(J) is _packed.PatternMatrix)
        return real(self, F, J)
    monkeypatch.setattr(dcanalysis._Held, 'system', system)
    for hi, lo in (('q', 'qb'), ('qb', 'q')):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            monkeypatch.setattr(_packed, 'PACKED_FROM', 1 << 30)
            dense = DC(_latch(), nodeset={hi: 5.0, lo: 0.0}).solve()
            monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
            dc = DC(_latch(), nodeset={hi: 5.0, lo: 0.0})
            got = dc.solve()
        assert got.v(hi) > 4.0 and got.v(lo) < 1.0
        assert np.allclose(np.asarray(got.x), np.asarray(dense.x), rtol=1e-9, atol=1e-12)
        ch = dc._get_linearsolver()._choice
        assert type(ch) is linearsolver.StructuralKLU and ch.m == dc.cir.n - 1
    assert any(seen) and not all(seen)


def test_a_matrix_the_structure_does_not_describe_takes_its_own_klu():
    """The AutoSolver on the pattern: a matrix of the structure's size takes
    KLU on the structure (the choice, kept); one of another size -- a held
    solve's -- KLU on its own pattern, before the choice and after it, which
    it leaves as it was; each solve numpy's to rounding."""
    from pycircuit.circuit import linearsolver
    rng = np.random.default_rng(7)
    A, D = _drawn(rng, 12, 20)
    A = A.plus_diagonal(1e4)
    D = A.to_dense()
    keep = np.array([0, 1, 3, 4, 6, 7, 9, 10, 11])
    B = A.kept(keep)
    ls = linearsolver.AutoSolver(structure=lambda: (12, A.pat.rows, A.pat.cols))
    b = rng.normal(size=12)
    held = ls._select(B)
    assert ls._choice is None and held is _packed.own_solver(B) and held.m == 9
    x = ls.solve(A, b, numeric)
    assert type(ls._choice) is linearsolver.StructuralKLU and ls._choice.m == 12
    assert ls._select(B) is held and ls._choice.m == 12
    xb = ls.solve(B, b[keep], numeric)
    assert np.allclose(x, np.linalg.solve(D, b), rtol=1e-10)
    assert np.allclose(xb, np.linalg.solve(D[np.ix_(keep, keep)], b[keep]), rtol=1e-10)


def test_drawn_starts_on_the_pattern_are_the_pythons(monkeypatch, ordinary_dc):
    """gm3's operating point on the pattern from rough drawn starts: the
    ladders' rungs, the C's hand-backs on a flag or a singular factor, its
    iteration limits -- every base solve, the counters, x0 and the warnings
    the C's as the packed Python's."""
    from pycircuit.circuit.tests import test_dc_newton_c as T
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    rng = np.random.default_rng(11)
    starts = [(10.0 ** rng.uniform(1, 4), int(rng.integers(1 << 30))) for _ in range(6)]
    bails = 0
    for scale, seed in (starts[0], starts[5]):
        def x0(n, scale=scale, seed=seed):
            return np.random.default_rng(seed).standard_normal(n) * scale
        py, _c = T._record('gm3', monkeypatch, False, x0)
        got, counts = T._record('gm3', monkeypatch, True, x0)
        assert got == py and counts.get('dcnewton_c:served', 0) > 0, counts
        bails += counts.get('dcnewton_c:bail:flags', 0) + counts.get('dcnewton_c:bail:singular', 0)
    assert bails > 0


def _exp_chain(cells=40):
    """A resistor chain loaded by exponential sources (1e-14 (exp(v / 25 mV)
    - 1) A): from a wild start their conductances overflow."""
    import sympy

    from pycircuit.circuit.elements_hdl import BSourceHdl
    c = SubCircuit()
    c['vs'] = VS('s', gnd, v=1.0)
    prev = 's'
    for k in range(cells):
        c[f'r{k}'] = R(prev, f'n{k}', r=1e3)
        c[f'd{k}'] = BSourceHdl(f'n{k}', gnd, f'n{k}', gnd,
                                i_func=lambda u: 1e-14 * (sympy.exp(u / 0.025) - 1))
        prev = f'n{k}'
    return c


def test_an_overflowing_matrix_on_the_pattern_hands_back(monkeypatch, ordinary_dc):
    """THE GUARDS on the pattern: from a start where the sources' conductances
    overflow, the C hands each solve back on the flag (or the non-finite
    entry) and the packed Python's ladders run -- their outcome, x or the
    exception, the same with the C as without; from a milder start the C
    serves."""
    if _dc_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_dc_newton_c.STATUS}')
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)

    def run(on, scale):
        monkeypatch.setattr(_dc_newton_c, 'ENABLED', on)
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dc = DC(_exp_chain(), toolkit=numeric)
            try:
                x0 = np.random.default_rng(3).standard_normal(dc.cir.n) * scale
                out = np.asarray(dc.solve(x0=x0).x).tobytes()
            except Exception as e:                             # noqa: BLE001
                out = (type(e).__name__, str(e))
        return out, _paths.since(before)
    for scale, key in ((20.0, 'dcnewton_c:bail:flags'), (5.0, 'dcnewton_c:served')):
        py, _c = run(False, scale)
        got, counts = run(True, scale)
        assert got == py and counts.get(key, 0) > 0, (scale, counts)


def _vccs_chain(k=6):
    """A pulse driving k VCCS stages, each loaded by R || C: G's places
    one-way (a stage's output row, its input column) -- an unsymmetric
    pattern, whose map between the core's order and the plan's is not its
    own inverse (on a symmetric one it is)."""
    from pycircuit.circuit.elements import VCCS, VPulse
    c = SubCircuit()
    c['vin'] = VPulse('n0', gnd, v1=0.0, v2=1.0, td=2e-9, tr=1e-9, tf=1e-9, pw=5e-9, per=2e-8)
    for i in range(k):
        c[f'g{i}'] = VCCS(f'n{i}', gnd, f'n{i + 1}', gnd, gm=1e-3 * (1 + i % 3))
        c[f'r{i}'] = R(f'n{i + 1}', gnd, r=1e3)
        c[f'c{i}'] = C(f'n{i + 1}', gnd, c=1e-12)
    return c


def _loop_cases():
    from pycircuit.circuit.integrator import EulerIntegrator, TrapezoidalIntegrator
    from pycircuit.circuit.tests import test_loop_c as L
    fixed = [
        ('mos1', L._mos, None, {'tend': 60 * 2e-8, 'timestep': 2e-8}),
        ('gp-chord', L._gp, None, {'tend': 60 * 2e-8, 'timestep': 2e-8}),
        ('psp-screen', L._psp_stage, None, {'tend': 60 * 2e-8, 'timestep': 2e-8}),
        ('sparse', lambda: L._mos(250, cap=1e-13), None, {'tend': 20 * 2e-10, 'timestep': 2e-10}),
        ('pulse', L._pulsed, None, {'tend': 4e-8, 'timestep': 2.5e-10}),
        ('trap', L._pulsed, {'integrator': TrapezoidalIntegrator()},
         {'tend': 3e-8, 'timestep': 2.5e-10}),
        ('euler', L._pulsed, {'integrator': EulerIntegrator()},
         {'tend': 3e-8, 'timestep': 2.5e-10}),
        ('vccs', _vccs_chain, None, {'tend': 3e-8, 'timestep': 2.5e-10}),
    ]
    adaptive = [
        ('a-gear', L._hard_mos, None, {'tend': 4e-6, 'timestep': 2e-8}),
        ('a-pulse', L._pulsed, None, {'tend': 6e-8, 'timestep': 1e-9}),
        ('a-gp-chord', L._gp, None, {'tend': 3e-6, 'timestep': 2e-8}),
        ('a-psp-screen', L._psp_stage, None, {'tend': 3e-6, 'timestep': 2e-8}),
        ('a-sparse', lambda: L._mos(250, cap=1e-13), {'reltol': 1e-5},
         {'tend': 1e-8, 'timestep': 2e-10}),
        ('a-judge', L._pulsed, {'lte_algebraic': 'judge'}, {'tend': 4e-8, 'timestep': 1e-9}),
        ('a-charge', L._pulsed_load, None, {'tend': 4e-8, 'timestep': 1e-9}),
        ('a-vccs', _vccs_chain, None, {'tend': 4e-8, 'timestep': 1e-9}),
    ]
    return ([(nm, b, mk, kw) for nm, b, mk, kw in fixed]
            + [(nm, b, mk, dict(kw, fixed_timestep=False)) for nm, b, mk, kw in adaptive])


@pytest.mark.parametrize('name, build, make, kw', _loop_cases(), ids=[c[0] for c in _loop_cases()])
def test_the_c_loop_on_the_pattern_is_the_python_loop(name, build, make, kw, monkeypatch):
    """The stepping loop's C on the pattern (THE PACKED LAYOUT from any size)
    against the Python loop with the transient Newton's C on it: every
    attempt -- accepted, rejected and retried, the screen's, the chord's --
    the same bytes, statistics, warnings, tallies and run state (its C, Geq
    and charges as PatternMatrix on the plan's pattern) and controller; the
    loop serving."""
    from pycircuit.circuit.tests import test_loop_c as L
    if _tran_newton_c.driver() is None:
        pytest.skip('the C Newton is off')
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    ## (a short run first: a process's first run can meet its source table's
    ## stamp moved by first-use work and leave the loop -- dense too)
    L._run(build, True, monkeypatch, make=make, **dict(kw, tend=3 * kw['timestep']))
    ca = L._same(build, monkeypatch, steps_min=10, make=make, **kw)
    assert ca.get('loop_c:steps', 0) >= 10, ca


def test_the_write_backs_and_records_on_the_pattern(monkeypatch):
    """After a run the stepping loop's C served on the pattern: the cache's C,
    Geq and the companion's C PatternMatrix on the plan's own pattern (no
    `_C_pattern`: the Python on the pattern sets none); the core packed, and
    none of the Newton's record's, the error test's or the core's buffers
    an n x n or m x m one."""
    from pycircuit.circuit import _stamp_plan
    from pycircuit.circuit.tests import test_loop_c as L
    if _tran_newton_c.driver() is None:
        pytest.skip('the C Newton is off')
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    cir = L._mos(250, cap=1e-13)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(cir, toolkit=numeric)
        tr.solve(tend=20 * 2e-10, timestep=2e-10, fixed_timestep=True)
    assert _paths.since(before).get('loop_c:steps', 0) > 10
    pat = _stamp_plan._plan_for(cir).pattern
    for v in (tr._C_cache[1], tr._Geq, tr._Cmat):
        assert type(v) is _packed.PatternMatrix and v.pat is pat
    assert tr.__dict__.get('_C_pattern') is None
    n = cir.n
    core = tr.__dict__['_tran_core'][1]
    rec = tr.__dict__['_newton_c']
    lte = tr.step_controller.__dict__.get('_lte_c') if tr.step_controller else None
    assert core.packed and core.arrays['binG'].size == 1
    assert rec.s.packed == 1 and rec.s.sets == 1 and rec.pk is not None
    bufs = [*rec.buf.values(), *core.arrays.values()] + ([] if lte is None else [*lte.buf.values()])
    assert max(a.size for a in bufs) < (n - 1) * (n - 1)


def test_an_entry_bounded_by_bytes_comes_back(monkeypatch):
    """The stepping loop's outputs an entry holds bounded by their bytes
    (`_tran_loop_c.CAP_BYTES`: at size fewer rows than CAP_MAX), written
    back and entered again: the same run as at CAP_MAX."""
    from pycircuit.circuit import _tran_loop_c
    from pycircuit.circuit.tests import test_loop_c as L
    if _tran_newton_c.driver() is None:
        pytest.skip('the C Newton is off')
    assert _tran_loop_c._cap(1000) == _tran_loop_c.CAP_MAX
    assert _tran_loop_c._cap(164898) == _tran_loop_c.CAP_BYTES // (8 * 164898) < 300
    n = L._mos().n
    monkeypatch.setattr(_tran_loop_c, 'CAP_BYTES', 8 * n * 7)
    assert _tran_loop_c._cap(n) == 7
    L._same(L._mos, monkeypatch, steps_min=40, tend=60 * 2e-8, timestep=2e-8)



def test_the_refactor_residuals_on_the_stepping_path_are_the_pythons(monkeypatch):
    """The Newton's and the error test's refactor residuals on the pattern
    (newton.c `sparse_full`, lte.c, the CSC kernel deciding every one) at a
    tolerance near their values -- most holding, some failing to a fresh
    factor (the C hands those back): the loop's run the Python loop's, the
    fallbacks made."""
    from pycircuit.circuit import linearsolver
    from pycircuit.circuit.tests import test_loop_c as L
    if _tran_newton_c.driver() is None:
        pytest.skip('the C Newton is off')
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    monkeypatch.setattr(linearsolver.KLUSolver, 'REFACTOR_RESIDUAL_TOL', 1e-16)
    kw = {'tend': 6e-8, 'timestep': 1e-9, 'fixed_timestep': False}
    L._run(L._pulsed, True, monkeypatch, **dict(kw, tend=3e-9))
    ca = L._same(L._pulsed, monkeypatch, steps_min=10, **kw)
    assert ca.get('newton_c:bail:klu_residual', 0) > 0 and ca.get('lte_c:bail:residual', 0) > 0, ca
