"""The Radau transform's factors once a step (`_tran_radau._radau_frozen`,
speed round 10, B3.1).  The transform's Newton is simplified -- `C` and `G`
frozen at `x_n` -- so its two factors and their factorisations are the same
at every iteration of a step; the frozen solve makes them once.  It is the
per-iteration solve's bytes: on whole transients and a PSS with the switch
off and on, and on one captured step's factors with drawn right-hand sides,
step sizes and special values, its warnings from the same lines and its
exceptions the same; its pieces -- numpy's solve kept as an LU
(`linearsolver.NumpyLU`), the complex factor marshalled by numpy
(`ComplexKLUSolver.prepare`, `_csc_of_dense`) and refactored once
(`solve_prepared`) -- each the bytes of what it stands in for; and the
declines."""
import warnings

import numpy as np
import pytest
import scipy.sparse as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _paths, _tran_radau, circuit
from pycircuit.circuit import linearsolver as LS
from pycircuit.circuit._tran_radau import _RadauStages
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_newton_c import _mos_chain
from pycircuit.circuit.tests.test_psp_limit_c import _stage
from pycircuit.circuit.transient import Transient

SPECIAL = np.array([0.0, -0.0, np.nan, np.inf, -np.inf, 5e-324, -5e-324, 1e308, -1e308])
_NP_SOLVE = np.linalg.solve


def _klu_or_skip():
    try:
        return LS.ComplexKLUSolver()
    except ImportError as e:
        pytest.skip(f'libklu not available: {e}')


def _lapack_or_skip():
    if LS._numpy_lapack() is None:
        pytest.skip("numpy's own OpenBLAS not found")


def _others(d):
    """The counts but the frozen path's own and the once-only ones."""
    return {k: v for k, v in d.items() if not k.startswith(('radau.frozen:', 'once:'))}


def _run(build, on, make=None, **kw):
    old = _tran_radau.RADAU_FROZEN
    _tran_radau.RADAU_FROZEN = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                           **(make or {}))
            before = _paths.snapshot()
            res = tr.solve(**kw)
            d = _paths.since(before)
    finally:
        _tran_radau.RADAU_FROZEN = old
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), _others(d)), d


CASES = {
    'psp': (_stage, {}),
    'psp-transform': (_stage, {'radau_transform': True}),
    'mos-transform': (_mos_chain, {'radau_transform': True}),
    'mos-caps-transform': (lambda: _mos_chain(cap=1e-13), {'radau_transform': True}),
}


@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
@pytest.mark.parametrize('case', list(CASES))
def test_a_radau_transient_is_the_same_with_the_frozen_transform_off(case, fixed):
    _klu_or_skip()
    build, make = CASES[case]
    kw = {'tend': 6e-7, 'timestep': 2e-8}
    if fixed:
        kw['fixed_timestep'] = True
    a, da = _run(build, False, make, **dict(kw))
    b, db = _run(build, True, make, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert a[4] == b[4], ('the other counts moved', a[4], b[4])
    assert da.get('radau.frozen:off', 0) > 0, da
    assert db.get('radau.frozen:served', 0) > 0, db


def test_a_radau_pss_is_the_same_with_the_frozen_transform_off():
    """PSS's default method on a PSP circuit: the transform under 'auto'."""
    _klu_or_skip()
    from pycircuit.circuit.shooting import PSS

    def run(on):
        old = _tran_radau.RADAU_FROZEN
        _tran_radau.RADAU_FROZEN = on
        try:
            with warnings.catch_warnings(record=True) as W:
                warnings.simplefilter('always')
                warnings.simplefilter('ignore', ResourceWarning)
                p = PSS(_stage(), method='radau', reltol=1e-8)
                before = _paths.snapshot()
                p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
                d = _paths.since(before)
        finally:
            _tran_radau.RADAU_FROZEN = old
        wf = p.waveform
        return (np.asarray(wf[0], float).tobytes(), np.asarray(wf[1], float).tobytes(),
                sorted(str(w.message) for w in W), _others(d)), d
    a, _da = run(False)
    b, db = run(True)
    assert a == b
    assert db.get('radau.frozen:served', 0) > 0, db


## -- one step's factors, captured: drawn right-hand sides and step sizes --------------

def _captured(build=_stage, steps=3, **make):
    """A radau transient on the transform path and its last step's frozen
    inputs: `(tr, Cr, Gr, h)`; its solver has chosen dense by then."""
    tr = Transient(build(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(), **make)
    got = {}
    real = tr._radau_frozen

    def spy(Cr, Gr, h):
        got.update(Cr=np.array(Cr), Gr=np.array(Gr), h=h)
        return real(Cr, Gr, h)
    tr._radau_frozen = spy
    tr.solve(tend=steps * 2e-8, timestep=2e-8, fixed_timestep=True)
    del tr._radau_frozen
    assert got, 'the transform path did not run'
    return tr, got['Cr'], got['Gr'], got['h']


def _answer(fn, *args):
    """`fn(*args)`'s bytes or its exception, and its warnings with their
    lines (numpy's warnings are shown once a line)."""
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        warnings.simplefilter('ignore', ResourceWarning)
        try:
            out = np.asarray(fn(*args)).tobytes()
        except Exception as e:                             # noqa: BLE001
            out = (type(e).__name__, str(e))
    return out, [(w.category.__name__, str(w.message), w.filename, w.lineno) for w in W]


#: (only the drawn test below reads it: a captured transient is recorded
#: against the test that made it)
_CAP = {}


def _cap(kind):
    if kind not in _CAP:
        _CAP[kind] = (_captured() if kind == 'psp'
                      else _captured(_mos_chain, radau_transform=True))
    return _CAP[kind]


def _sides(tr, Cr, Gr, h, R3s):
    """Each right-hand side solved by the per-iteration solve and by one
    frozen solve, each side on a complex solver of its own from its first
    factor: the same sequence of factors and refactors but the frozen
    side's skipped ones."""
    tr._radau_zsolver = LS.ComplexKLUSolver()
    fz = tr._radau_frozen(Cr, Gr, h)
    tr._radau_zsolver = LS.ComplexKLUSolver()
    ref = [_answer(tr._radau_transform_solve, R3, Cr, Gr, h) for R3 in R3s]
    return fz, ref, (None if fz is None else [_answer(fz.solve, R3) for R3 in R3s])


@settings(deadline=None)
@given(kind=st.sampled_from(['psp', 'mos']), seed=st.integers(0, 2 ** 32 - 1),
       scale=st.integers(-40, 40), hk=st.integers(-4, 2),
       special=st.sampled_from(['none', 'none', 'one', 'some', 'zero']))
def test_the_frozen_solve_is_the_per_iteration_solve(kind, seed, scale, hk, special):
    """`_FrozenTransform.solve(R3)` against `_radau_transform_solve(R3, Cr,
    Gr, h)` on one step's factors, three right-hand sides: the bytes or the
    exception, and every warning from the same line."""
    _klu_or_skip()
    tr, Cr, Gr, h0 = _cap(kind)
    h = h0 * 10.0 ** hk
    m = Cr.shape[0]
    rng = np.random.default_rng(seed)
    R3s = []
    for _ in range(3):
        R3 = rng.standard_normal(3 * m) * 10.0 ** scale
        if special == 'one':
            R3[rng.integers(3 * m)] = rng.choice(SPECIAL)
        elif special == 'some':
            k = rng.random(3 * m) < 0.3
            R3[k] = rng.choice(SPECIAL, size=int(k.sum()))
        elif special == 'zero':
            R3[:] = 0.0
        R3s.append(R3)
    before = _paths.snapshot()
    fz, ref, got = _sides(tr, Cr, Gr, h, R3s)
    if fz is None:
        assert _paths.since(before) == {'radau.frozen:fp': 1}
        return
    assert fz.lu is not None and fz.prep is not None
    assert got == ref


def test_singular_factors_raise_as_the_per_iteration_solve():
    """A singular real factor raises numpy's error at every solve, from the
    kept LU as from numpy's solve; a non-finite factor is solved by the
    analysis solver, as before."""
    _klu_or_skip()
    tr, Cr, Gr, h = _captured()
    m = Cr.shape[0]
    R3 = np.linspace(-1.0, 1.0, 3 * m)
    sing_C, sing_G = np.zeros((m, m)), np.array(Gr)
    sing_G[1, :] = 0.0
    inf_C = np.array(Cr)
    inf_C[0, 0] = np.inf
    for C_, G_, lu in ((sing_C, sing_G, True), (inf_C, Gr, False)):
        fz, ref, got = _sides(tr, C_, G_, h, [R3, R3, -R3])
        assert fz is not None and (fz.lu is not None) == lu
        assert got == ref
    assert ref[0][0] == ('LinAlgError', 'Singular matrix') or lu is False
    assert _sides(tr, Cr, Gr, h, [R3])[1] != ref[:1]


def test_what_the_frozen_transform_declines(monkeypatch):
    """The switch, a patched or shadowed piece of the transform, a factor
    whose arithmetic would raise; a patched solver underneath is called as
    it is (the frozen solve keeps no LU or no marshalled factor for it)."""
    _klu_or_skip()
    tr, Cr, Gr, h = _captured()

    def counted(**kw):
        before = _paths.snapshot()
        fz = tr._radau_frozen(Cr, Gr, kw.get('h', h))
        return fz, {k: v for k, v in _paths.since(before).items() if k.startswith('radau.')}
    monkeypatch.setattr(_tran_radau, 'RADAU_FROZEN', False)
    assert counted() == (None, {'radau.frozen:off': 1})
    monkeypatch.undo()
    assert counted()[1] == {'radau.frozen:served': 1}
    real = _RadauStages._radau_transform_solve
    monkeypatch.setattr(_RadauStages, '_radau_transform_solve',
                        lambda self, *a: real(self, *a))
    assert counted() == (None, {'radau.frozen:patched': 1})
    monkeypatch.undo()
    tr._radau_complex_solve = tr._radau_complex_solve
    assert counted() == (None, {'radau.frozen:patched': 1})
    del tr._radau_complex_solve
    assert counted(h=1e-320) == (None, {'radau.frozen:fp': 1})
    ## the solvers underneath: called, not stood in for
    monkeypatch.setattr(LS.AutoSolver, 'solve', lambda self, A, b, tk: np.linalg.solve(A, b))
    fz, d = counted()
    assert fz.lu is None and fz.prep is not None and d == {'radau.frozen:served': 1}
    monkeypatch.undo()
    monkeypatch.setattr(LS.ComplexKLUSolver, 'solve',
                        lambda self, A, b: self.solve_prepared(self.prepare(A), b))
    fz, d = counted()
    assert fz.lu is not None and fz.prep is None and d == {'radau.frozen:served': 1}
    monkeypatch.undo()
    monkeypatch.setattr(np.linalg, 'solve', lambda A, b: _NP_SOLVE(A, b))
    assert counted()[0].lu is None
    monkeypatch.undo()
    monkeypatch.setattr(tr.toolkit, 'linearsolver', lambda A, b: _NP_SOLVE(A, b))
    assert counted()[0].lu is None


@pytest.mark.parametrize('solver', ['dense', 'superlu'])
def test_another_linear_solver_is_served(solver):
    """`linearsolver=DenseSolver()` keeps the LU; a sparse one is called."""
    _klu_or_skip()
    ls = LS.DenseSolver() if solver == 'dense' else LS.SuperLUSolver()
    a, _da = _run(_stage, False, {'linearsolver': ls}, tend=2e-7, timestep=2e-8,
                  fixed_timestep=True)
    ls = LS.DenseSolver() if solver == 'dense' else LS.SuperLUSolver()
    b, db = _run(_stage, True, {'linearsolver': ls}, tend=2e-7, timestep=2e-8,
                 fixed_timestep=True)
    assert a == b
    assert db.get('radau.frozen:served', 0) > 0, db


## -- the pieces ------------------------------------------------------------------------

@pytest.mark.parametrize('threads', [1, 4])
def test_the_kept_lu_is_numpys_solve(threads):
    """`NumpyLU(A).solve(b)` is `numpy.linalg.solve(A, b)`, bytes or error,
    for every `b` -- from 1 to 300 unknowns, on one BLAS thread and on
    four: a `dgetrf` first would differ from 100 unknowns up there (numpy's
    `dgesv` factors on one thread)."""
    _lapack_or_skip()
    tpc = pytest.importorskip('threadpoolctl')
    rng = np.random.default_rng(7 + threads)
    with tpc.threadpool_limits(limits=threads, user_api='blas'):
        for k in range(160):
            n = (1, 2, 3, 7, 30, 99, 100, 101, 160, 300)[k % 10]
            A = rng.standard_normal((n, n)) * 10.0 ** rng.integers(-9, 3, size=(n, n))
            if k % 3 == 0:
                A[rng.random((n, n)) < 0.7] = 0.0
                A[np.arange(n), np.arange(n)] += 1e-3
            if k % 7 == 0:
                A[:, rng.integers(n)] = 0.0
            lu = LS.NumpyLU.make(A)
            for _ in range(3):
                b = rng.standard_normal(n)
                assert _answer(lu.solve, b) == _answer(np.linalg.solve, A, b), n


def test_prepare_is_scipys_marshalling():
    """A dense square complex128 matrix marshalled by numpy: SciPy's CSC
    arrays, pattern key and product, bytes -- signed zeros, NaN, infinities
    and subnormals among the entries and the vectors; other inputs go
    through SciPy."""
    zs = _klu_or_skip()
    rng = np.random.default_rng(11)
    for k in range(1500):
        n = int(rng.integers(1, 30))
        A = ((rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
             * 10.0 ** rng.integers(-12, 4, size=(n, n)))
        if k % 5 >= 1:
            A[rng.random((n, n)) < 0.6] = 0.0
        if k % 5 >= 2:
            for part in (A.real, A.imag):
                s = rng.random((n, n)) < 0.1
                part[s] = rng.choice(SPECIAL, size=int(s.sum()))
        if k % 5 == 4:
            A = A.real + 0j
            A.imag[rng.random((n, n)) < 0.3] = -0.0
        p = zs.prepare(A)
        Acsc = sp.csc_matrix(A).astype(np.complex128)
        Ap = np.ascontiguousarray(Acsc.indptr, dtype=np.int32)
        Ai = np.ascontiguousarray(Acsc.indices, dtype=np.int32)
        Ax = np.ascontiguousarray(Acsc.data, dtype=np.complex128).view(np.float64)
        assert p[1] == n and p[5] == (n, Ap.tobytes(), Ai.tobytes())
        assert [(a.dtype, a.tobytes()) for a in p[2:5]] == [(a.dtype, a.tobytes())
                                                             for a in (Ap, Ai, Ax)]
        assert p[4].flags.c_contiguous
        for _ in range(2):
            x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
            if k % 3 == 0:
                s = rng.random(n) < 0.2
                x.real[s] = rng.choice(SPECIAL, size=int(s.sum()))
            with np.errstate(all='ignore'):
                assert p[0](x).tobytes() == Acsc.dot(x).tobytes()
    A = np.array([[2.0, 0.0], [1.0, 3.0]])
    for other in (A, sp.csc_matrix(A + 0j), (A + 0j)[:, :1], A.astype(np.complex64)):
        assert not hasattr(zs.prepare(other)[0], 'func'), type(other)


def test_solve_prepared_refactors_once_a_matrix():
    """Against a solver that refactors on every reused call (`solve` before
    B3.1): the same bytes on every call, one refactor a matrix where the
    same record comes back, and a residual fallback's factor followed by a
    refactor as before."""
    zs = _klu_or_skip()
    ref = LS.ComplexKLUSolver()
    rng = np.random.default_rng(5)
    n = 12
    base = (np.diag(5.0 + 3.0j + rng.random(n)) + np.diag((-1.0 - 0.5j) * np.ones(n - 1), 1)
            + np.diag((-1.0 + 0.2j) * np.ones(n - 1), -1))
    preps = [zs.prepare(base * s) for s in (1.0, 2.0 + 0.5j, 1.0)]

    def check(i, b):
        ref._fresh = None                   # (the old rule: every reused call refactors)
        a = ref.solve(base * (1.0, 2.0 + 0.5j, 1.0)[i], b)
        assert zs.solve_prepared(preps[i], b).tobytes() == a.tobytes()
    order = [0, 0, 0, 1, 1, 0, 2, 2, 0]
    for i in order:
        check(i, rng.standard_normal(n) + 1j * rng.standard_normal(n))
    ## factor, then one refactor a run of the same record: 0 | 0 0 | 1 1 | 0 | 2 2 | 0
    assert (zs.analyses, zs.factors, zs.refactors) == (1, 1, 5), zs
    assert (ref.analyses, ref.factors, ref.refactors) == (1, 1, 8), ref
    ## every reused call falls back: a factor each, then a refactor as before
    zs.REFACTOR_RESIDUAL_TOL = ref.REFACTOR_RESIDUAL_TOL = -1.0
    for i in (0, 0, 1, 0):
        check(i, rng.standard_normal(n) + 1j * rng.standard_normal(n))
    assert zs.residual_fallbacks == ref.residual_fallbacks == 4
    ## (the first call's refactor skipped: the last one was of its record)
    assert (zs.refactors - 5, ref.refactors - 8) == (3, 4), (zs, ref)
