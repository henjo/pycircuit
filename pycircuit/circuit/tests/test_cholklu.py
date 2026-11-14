"""THE SPD REDUCTION ON CHOLMOD (`_cholklu`): a linear grid's transient
matrix solved through KLU's entries over CHOLMOD's supernodal Cholesky --
the voltage sources' nodes fixed, the inductors' currents eliminated, the
free nodes' symmetric positive definite system factored -- as a dense solve
makes it; the plan declining what is not a grid; a matrix CHOLMOD will not
factor KLU's; a transient on it its KLU run, its C the Python's."""
import hashlib
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _cholklu, _paths, linearsolver, numeric
from pycircuit.circuit.analysis import jacobian_structure
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import IS, VS, C, Diode, IPulse, R
from pycircuit.circuit.elements import L as Ind
from pycircuit.circuit.transient import Transient


def _grid(nx=6, ny=5):
    """An RC mesh fed by a supply through an inductor (an L from a fixed node
    to a free one), a second supply on the far corner, an inductor across
    the mesh and one to ground, a current source."""
    c = SubCircuit()
    nm = lambda i, j: f'n{i}_{j}'
    c['vdd'] = VS('s', gnd, v=1.0)
    c['lp'] = Ind('s', nm(0, 0), L=1e-10)
    c['v2'] = VS(nm(nx - 1, ny - 1), gnd, v=0.9)
    c['is'] = IS(nm(2, 2), gnd, i=1e-3)
    k = 0
    for i in range(nx):
        for j in range(ny):
            if i + 1 < nx:
                c[f'r{k}'] = R(nm(i, j), nm(i + 1, j), r=1.0 + 0.1 * (k % 7))
                k += 1
            if j + 1 < ny:
                c[f'r{k}'] = R(nm(i, j), nm(i, j + 1), r=2.0 + 0.1 * (k % 5))
                k += 1
            c[f'c{i}_{j}'] = C(nm(i, j), gnd, c=1e-12 * (1 + (i + j) % 3))
    c['lx'] = Ind(nm(1, 1), nm(3, 3), L=2e-10)
    c['lg'] = Ind(nm(4, 4), gnd, L=3e-10)
    ## (a supply made last, its inductor from an early free node: the free end
    ## the lower row, as the first supply's inductor has it the higher)
    c['v3'] = VS('late', gnd, v=0.8)
    c['ll'] = Ind(nm(1, 0), 'late', L=1.5e-10)
    return c


def _matrix(c, h):
    iref = c.get_node_index(gnd)
    x = np.zeros(c.n)
    J = np.asarray(c.G(x)) + np.asarray(c.C(x)) / h
    return np.delete(np.delete(J, iref, 0), iref, 1), iref


def _solver(c, plan=True):
    iref = c.get_node_index(gnd)
    m, rows, cols = jacobian_structure(c, iref)
    return linearsolver.StructuralKLU(m, rows, cols, fallback=linearsolver.SuperLUSolver,
                                      cholmod=_cholklu.plan(c, iref) if plan else None)


@pytest.fixture(autouse=True)
def _cholmod_here():
    if _cholklu.library(linearsolver.KLUSolver()._lib) is None:
        pytest.skip(f'no CHOLMOD here: {_cholklu.STATUS}')


def test_the_reduction_solves_as_a_dense_solve():
    """On the reduction: a factor, a refactor (h doubled), several
    right-hand sides, the transposed solve and an owned factor -- each the
    dense solve's to rounding, its backward error at the arithmetic's."""
    c = _grid()
    sk = _solver(c)
    assert sk.cholmod
    rng = np.random.default_rng(0)
    for h in (1e-12, 2e-12):
        J, _iref = _matrix(c, h)
        b = rng.standard_normal(J.shape[0])
        x = np.asarray(sk.solve(J, b, None))
        assert np.max(np.abs(x - np.linalg.solve(J, b))) <= 1e-13 * np.max(np.abs(x))
        r = J @ x - b
        assert np.max(np.abs(r)) <= 1e-15 * np.max(np.abs(J) @ np.abs(x) + np.abs(b))
        B = rng.standard_normal((J.shape[0], 3))
        X = np.asarray(sk.solve(J, B, None))
        assert np.max(np.abs(X - np.linalg.solve(J, B))) <= 1e-13 * np.max(np.abs(X))
        XT = np.asarray(sk.solve_transposed(J, B, None))
        assert np.max(np.abs(XT - np.linalg.solve(J.T, B))) <= 1e-13 * np.max(np.abs(XT))
        f = sk.factor(J, None)
        assert np.max(np.abs(np.asarray(f.solve(b)) - x)) <= 1e-13 * np.max(np.abs(x))
    assert sk.factors == 1 and sk.refactors >= 1 and sk.cholmod


def test_a_refactor_check_failing_on_the_reduction_is_kept():
    """The refactor's residual check refused (a tolerance of 0): the
    reduction's growth entry reports a sound factor (a Cholesky of a positive
    definite matrix has no pivot order to go stale), so the refactor is
    kept (THE REFACTOR KEPT ON ITS GROWTH) -- no fresh factor, counted."""
    c = _grid()
    sk = _solver(c)
    sk.REFACTOR_RESIDUAL_TOL = 0.0
    b = np.random.default_rng(3).standard_normal(sk.m)
    for h in (1e-12, 2e-12):
        J, _iref = _matrix(c, h)
        x = np.asarray(sk.solve(J, b, None))
        assert np.max(np.abs(x - np.linalg.solve(J, b))) <= 1e-13 * np.max(np.abs(x))
    assert (sk.factors, sk.refactors, sk.residual_fallbacks, sk.residual_kept) == (1, 1, 0, 1)


def test_the_plan_declines_what_is_not_a_grid():
    """A device, a voltage source between two nodes, two on one node, a
    subcircuit instance, an instance's own stamp: each named; the grid's
    rows counted by kind."""
    c = _grid()
    got = _cholklu.plan(c, c.get_node_index(gnd))
    assert list(np.bincount(got)) == [c.n - 1 - 7, 3, 4]
    cases = {}
    d = _grid()
    d['d1'] = Diode('n1_2', gnd)
    cases['element'] = d
    f = _grid()
    f['vf'] = VS('n1_2', 'n2_2', v=0.1)
    cases['vsource'] = f
    t = _grid()
    t['vt'] = VS('s', gnd, v=0.1)
    cases['vsource '] = t

    class Stage(SubCircuit):
        terminals = ('a', 'b')

        def __init__(self, *args, **kw):
            super().__init__(*args, **kw)
            self['R1'] = R('a', 'b', r=1.0)
    s = _grid()
    s['X1'] = Stage('n0_1', 'n1_1')
    cases['hierarchy'] = s
    o = _grid()
    G0 = o['r0'].G
    o['r0'].G = lambda x, *a: G0(x, *a)
    cases['instance'] = o
    for why, cir in cases.items():
        assert _cholklu.plan(cir, cir.get_node_index(gnd)) == why.strip(), why


def test_a_matrix_cholmod_refuses_is_klus(monkeypatch):
    """The inductors' diagonals zero (h infinite: the DC's matrix) -- CHOLMOD
    refuses the factor -- and a pattern a voltage source's row couples two
    nodes in: each solved by KLU, counted, the dense solve's answer."""
    c = _grid()
    sk = _solver(c)
    J, _iref = _matrix(c, np.inf)
    b = np.random.default_rng(1).standard_normal(J.shape[0])
    before = _paths.snapshot()
    x = np.asarray(sk.solve(J, b, None))
    assert _paths.since(before).get('solver.cholmod:klu') == 1
    assert np.max(np.abs(x - np.linalg.solve(J, b))) <= 1e-12 * np.max(np.abs(x))
    assert sk.cholmod and sk._ck_twin is not None and sk._ck_twin.factors == 1
    ## (a pattern not a grid's: KLU from the analysis, nothing made on CHOLMOD)
    iref = c.get_node_index(gnd)
    m, rows, cols = jacobian_structure(c, iref)
    kind = _cholklu.plan(c, iref).copy()
    kind[np.flatnonzero(kind == 0)[3]] = 1
    before = _paths.snapshot()
    bad = linearsolver.StructuralKLU(m, rows, cols, fallback=linearsolver.SuperLUSolver, cholmod=kind)
    J1, _iref = _matrix(c, 1e-12)
    x = np.asarray(bad.solve(J1, b, None))
    assert _paths.since(before).get('solver.cholmod:pattern') == 1 and not bad.cholmod
    assert np.max(np.abs(x - np.linalg.solve(J1, b))) <= 1e-12 * np.max(np.abs(x))


def test_the_threshold_the_switch_and_the_dc_keep_klu(monkeypatch):
    """Below `CHOLMOD_FROM`, with the switch off, and in the DC (its
    inductors shorts) or any analysis but the transient (`_SPD_PLAN`): KLU;
    from the threshold, a transient's matrix the reduction's
    (`solver.auto:cholmod`)."""
    def choice(make):
        c = _grid(10, 10)
        c['ip'] = IPulse('n2_3', 'n0_0', i1=0.0, i2=1e-3, td=1e-11, tr=1e-11, tf=1e-11, pw=1e-10,
                         per=4e-10)
        a = make(c)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if isinstance(a, Transient):
                a.solve(tend=5e-11, timestep=1e-11)
            else:
                a.solve()
        return a._get_linearsolver()

    def tran(c):
        return Transient(c, toolkit=numeric)

    def dc(c):
        return DC(c, toolkit=numeric)
    assert not choice(tran)._choice.cholmod
    monkeypatch.setattr(_cholklu, 'CHOLMOD_FROM', 0)
    before = _paths.snapshot()
    assert choice(tran)._choice.cholmod and _paths.since(before).get('solver.auto:cholmod') == 1
    assert choice(dc).spd_plan is None
    from pycircuit.circuit.analysis_ss import AC
    assert AC(_grid(), toolkit=numeric)._get_linearsolver().spd_plan is None
    monkeypatch.setattr(_cholklu, 'ENABLED', False)
    assert not choice(tran)._choice.cholmod


def _run(chol, c_on, monkeypatch, **make):
    from pycircuit.circuit import _tran_loop_c, _tran_lte_c, _tran_newton_c
    monkeypatch.setattr(_cholklu, 'CHOLMOD_FROM', 0 if chol else 1 << 40)
    for mod in (_tran_loop_c, _tran_lte_c, _tran_newton_c):
        monkeypatch.setattr(mod, 'ENABLED', c_on)
    c = _grid(16, 16)
    c['ip'] = IPulse('n8_8', 'n0_0', i1=0.0, i2=5e-3, td=1e-10, tr=1e-10, tf=1e-10, pw=4e-10,
                     per=2e-9)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(c, toolkit=numeric, **make)
        res = tr.solve(tend=3e-9, timestep=1e-11)
    st = tr.statistics
    x = np.asarray(res.x, float)
    sk = tr._get_linearsolver()._choice
    return (x, hashlib.sha1(x.tobytes()).hexdigest(), (st.accepted_steps, st.rejected_steps,
            st.newton_iterations), (sk.cholmod, sk.factors, sk.refactors, sk.residual_fallbacks,
            sk.chord_factors, sk.chord_refactors), _paths.since(before))


@pytest.mark.parametrize('chord', [True, False], ids=['chord', 'full'])
def test_a_grid_transient_on_cholmod_is_its_klu_run_and_its_c_the_python(chord, monkeypatch):
    """A linear grid's Gear-2 run, the threshold at 0: on the reduction the C
    Newton, error test and stepping loop serve (their KLU entries the
    reduction's), the run the Python's byte for byte; against KLU's run the
    same steps, the answers within rounding."""
    make = {'chord_jacobian': chord}
    k = _run(False, True, monkeypatch, **make)
    c = _run(True, True, monkeypatch, **make)
    p = _run(True, False, monkeypatch, **make)
    assert c[1:4] == p[1:4] and c[3][0], (c[2:4], p[2:4])
    assert c[2] == k[2] and not k[3][0]
    assert float(np.max(np.abs(c[0] - k[0]))) <= 1e-12
    d = c[4]
    assert d.get('solver.auto:cholmod') == 1 and d.get('loop_c:steps', 0) > 50, d
    assert d.get('newton_c:served', 0) > 50 and d.get('lte_c:sparse_served', 0) > 50, d
    assert not any(key.startswith(('newton_c:bail', 'lte_c:bail', 'solver.cholmod:')) for key in d), d
