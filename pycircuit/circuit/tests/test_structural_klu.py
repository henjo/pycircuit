"""KLU on the circuit's structural pattern (the C stepping loop plan's stage
2; Andreas's decision, 2026-10-10: every sparse run).

`analysis.jacobian_structure` must cover every nonzero the circuit's
reduced Jacobian can hold -- G's, C's, `G + a C`'s, at any state -- and
`linearsolver.StructuralKLU` must solve what the dense solver solves, hand
any matrix it does not cover to today's choice without disturbing its own
factorisation, and `AutoSolver` must pick it for a sparse run exactly when
the analysis gives it a structure.
"""
import numpy as np
import pytest

from pycircuit.circuit import elements as E
from pycircuit.circuit import gnd, linearsolver, numeric
from pycircuit.circuit.analysis import _reduce_ndarray, jacobian_structure
from pycircuit.circuit.circuit import SubCircuit, defaultepar
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.circuit.transient import Transient

klu = pytest.importorskip('pycircuit.circuit.linearsolver')
try:
    linearsolver.KLUSolver()
except ImportError:                                    # pragma: no cover
    pytest.skip('no KLU library', allow_module_level=True)


def _ladder(sections=300):
    """An RC ladder with a diode every tenth section and a pulse at its
    head: ~300 unknowns (one node a section), a sparse Jacobian."""
    c = SubCircuit()
    c['vs'] = E.VPulse('n0', gnd, v1=0.0, v2=1.0, td=1e-9, tr=1e-9, tf=1e-9, pw=2e-8, per=5e-8)
    for k in range(sections):
        a, b = f'n{k}', f'n{k + 1}'
        c[f'R{k}'] = E.R(a, b, r=100.0)
        c[f'C{k}'] = E.C(b, gnd, c=1e-12)
        if k % 10 == 5:
            c[f'D{k}'] = E.Diode(b, gnd)
    return c


def _reduced(M, iref):
    return _reduce_ndarray(np.asarray(M, dtype=float), iref)


def test_the_structure_covers_every_nonzero_of_the_jacobian():
    c = _ladder()
    iref = c.get_node_index(gnd)
    m, r, cc = jacobian_structure(c, iref)
    assert m == c.n - 1
    pat = np.zeros((m, m), dtype=bool)
    pat[r, cc] = True
    rng = np.random.default_rng(1)
    for _ in range(5):
        x = rng.normal(0.0, 0.7, c.n)
        G = _reduced(c.G(x, defaultepar), iref)
        C = _reduced(c.C(x, defaultepar), iref)
        for M in (G, C, G + 3.7e9 * C):
            assert not np.any((M != 0) & ~pat)
    ## CSC order: by column, then row; no duplicates
    key = cc * m + r
    assert np.all(np.diff(key) > 0)


def _structural(c):
    iref = c.get_node_index(gnd)
    st = jacobian_structure(c, iref)
    made = []

    def fallback():
        made.append(linearsolver.SuperLUSolver())
        return made[-1]
    return iref, st, linearsolver.StructuralKLU(*st, fallback=fallback), made


def test_it_solves_what_the_dense_solver_solves():
    c = _ladder()
    iref, st, s, _made = _structural(c)
    rng = np.random.default_rng(2)
    import numpy.linalg as la
    for _ in range(4):
        x = rng.normal(0.0, 0.5, c.n)
        J = _reduced(c.G(x, defaultepar), iref) + 1e10 * _reduced(c.C(x, defaultepar), iref)
        b = rng.normal(size=st[0])
        B = rng.normal(size=(st[0], 3))
        np.testing.assert_allclose(s.solve(J, b, numeric), la.solve(J, b), rtol=1e-9, atol=1e-12)
        np.testing.assert_allclose(s.solve(J, B, numeric), la.solve(J, B), rtol=1e-9, atol=1e-12)
        z = b + 1j * rng.normal(size=st[0])
        np.testing.assert_allclose(s.solve(J, z, numeric), la.solve(J, z), rtol=1e-9, atol=1e-12)
        f = s.factor(J, numeric)
        np.testing.assert_allclose(f.solve(b), la.solve(J, b), rtol=1e-9, atol=1e-12)
        np.testing.assert_allclose(f.solve_transposed(b), la.solve(J.T, b), rtol=1e-9, atol=1e-12)
    assert s.analyses == 1 and s.foreign == 0 and s.refactors >= 3


def test_a_matrix_it_does_not_cover_goes_to_the_fallback_and_leaves_its_factors():
    c = _ladder()
    iref, st, s, made = _structural(c)
    m = st[0]
    x = np.zeros(c.n)
    J = _reduced(c.G(x, defaultepar), iref) + 1e10 * _reduced(c.C(x, defaultepar), iref)
    b = np.ones(m)
    s.solve(J, b, numeric)
    num = s._numeric
    rng = np.random.default_rng(3)
    D = rng.normal(size=(m, m)) + m * np.eye(m)          # dense: off the pattern
    np.testing.assert_allclose(s.solve(D, b, numeric), np.linalg.solve(D, b), rtol=1e-9)
    assert s.foreign == 1 and len(made) == 1 and s._numeric == num
    f = s.factor(D, numeric)
    np.testing.assert_allclose(f.solve(b), np.linalg.solve(D, b), rtol=1e-9)
    assert s.foreign == 2 and len(made) == 1


def test_a_refactor_that_fails_its_residual_is_factored_again():
    c = _ladder()
    iref, st, s, _made = _structural(c)
    x = np.zeros(c.n)
    J = _reduced(c.G(x, defaultepar), iref) + 1e10 * _reduced(c.C(x, defaultepar), iref)
    b = np.ones(st[0])
    s.solve(J, b, numeric)
    s.REFACTOR_RESIDUAL_TOL = -1.0                       # every reuse fails its check
    y = s.solve(2.0 * J, b, numeric)
    np.testing.assert_allclose(y, np.linalg.solve(2.0 * J, b), rtol=1e-9)
    assert s.residual_fallbacks == 1 and s.factors == 2


def test_autosolver_takes_it_exactly_where_the_analysis_gives_a_structure():
    c = _ladder()
    iref = c.get_node_index(gnd)
    J = _reduced(c.G(np.zeros(c.n), defaultepar), iref) + np.eye(c.n - 1)
    b = np.ones(c.n - 1)
    a = linearsolver.AutoSolver(structure=lambda: jacobian_structure(c, iref))
    a.solve(J, b, numeric)
    assert type(a._choice) is linearsolver.StructuralKLU
    a = linearsolver.AutoSolver()
    a.solve(J, b, numeric)
    assert type(a._choice) is linearsolver.SuperLUSolver


def test_a_sparse_transient_runs_on_it_and_agrees_with_the_dense_solver():
    out = {}
    for name, ls in (('auto', None), ('dense', linearsolver.DenseSolver())):
        kw = {} if ls is None else {'linearsolver': ls}
        tr = Transient(_ladder(), reltol=1e-6, **kw)
        with quiet(AccuracyWarning):
            res = tr.solve(tend=3e-8, timestep=1e-10)
        out[name] = np.asarray(res.v('n300').y)
        if ls is None:
            assert type(tr._get_linearsolver()._choice) is linearsolver.StructuralKLU
    np.testing.assert_allclose(out['auto'], out['dense'], rtol=0, atol=1e-6)


def test_a_sparse_run_s_error_test_solves_on_its_factorisation():
    """The error test's `J^-1 Eg` on a sparse run is solved by the run's
    `StructuralKLU` (the C error test declines: its `dgesv` is not the
    reference there); a dense run's stays numpy's, and the C serves."""
    from pycircuit.circuit import _paths
    before = _paths.snapshot()
    tr = Transient(_ladder(), reltol=1e-6)
    with quiet(AccuracyWarning):
        tr.solve(tend=3e-8, timestep=1e-10)
    got = _paths.since(before)
    ls = tr._get_linearsolver()
    assert tr._sparse_linsolver() is ls and type(ls._choice) is linearsolver.StructuralKLU
    assert got.get('lte_c:sparse', 0) > 0 and ls._choice.refactors > 0
    before = _paths.snapshot()
    small = Transient(_ladder(20), reltol=1e-6)
    with quiet(AccuracyWarning):
        small.solve(tend=3e-8, timestep=1e-10)
    got = _paths.since(before)
    assert small._sparse_linsolver() is None
    assert got.get('lte_c:sparse', 0) == 0


## -- the operating point on the analysis's solver (Andreas, 2026-10-10) ----

def _dc_ladder(sections):
    """A resistive ladder driven at DC with a diode every tenth node:
    current flows, the junctions are on."""
    c = SubCircuit()
    c['vs'] = E.VS('n0', gnd, v=2.0)
    for k in range(sections):
        a, b = f'n{k}', f'n{k + 1}'
        c[f'R{k}'] = E.R(a, b, r=50.0)
        c[f'G{k}'] = E.R(b, gnd, r=5e3)
        if k % 10 == 5:
            c[f'D{k}'] = E.Diode(b, gnd)
    return c


def _op(c, dense, monkeypatch):
    from pycircuit.circuit import numeric
    from pycircuit.circuit.dcanalysis import DC
    dc = DC(c, toolkit=numeric)
    if dense:
        ## (the call before: no solver reaches the chain)
        monkeypatch.setattr(dc, '_get_linearsolver', lambda: None)
    x = np.asarray(dc.solve().x, dtype=float)
    monkeypatch.undo()
    return x, dc.__dict__.get('_default_linearsolver')


def test_a_small_operating_point_is_the_dense_call_byte_for_byte(monkeypatch):
    c = _dc_ladder(30)
    got, ls = _op(c, False, monkeypatch)
    ref, _ = _op(_dc_ladder(30), True, monkeypatch)
    assert got.tobytes() == ref.tobytes()
    assert ls._choice is ls._dense


def test_a_sparse_operating_point_factors_on_the_structure(monkeypatch):
    got, ls = _op(_dc_ladder(300), False, monkeypatch)
    ref, _ = _op(_dc_ladder(300), True, monkeypatch)
    assert type(ls._choice) is linearsolver.StructuralKLU and ls._choice.foreign == 0
    np.testing.assert_allclose(got, ref, rtol=1e-12, atol=1e-15)
