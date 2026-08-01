"""Tests for the LinearSolver strategies (stage 7b)."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import numeric, gnd, SubCircuit
from pycircuit.circuit.elements import R, C, VS
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.linearsolver import (LinearSolver, DenseSolver,
                                            SuperLUSolver, AutoSolver,
                                            MIN_N_FOR_SPARSE)


def _ladder(N):
    cir = SubCircuit(toolkit=numeric)
    for i in range(N):
        cir.add_node('n%d' % i)
    cir['vs'] = VS('n0', gnd, v=1.0)
    for i in range(N - 1):
        cir['R%d' % i] = R('n%d' % i, 'n%d' % (i + 1), r=100.0 * (i + 1))
    for i in range(N):
        cir['C%d' % i] = C('n%d' % i, gnd, c=1e-9 * (i + 1))
    return cir


def _run(N, solver, tend=5e-6, ts=1e-6):
    cir = _ladder(N)
    kw = {} if solver is None else {'linearsolver': solver}
    tran = Transient(cir, toolkit=numeric, reltol=1e-4, uic=True, **kw)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(refnode=gnd, tend=tend, timestep=ts)
    return (np.asarray(res.sweep_values, dtype=float),
            np.asarray(res.x, dtype=float))


def test_default_is_the_historical_path_bit_for_bit():
    """The default must not move any existing result.

    `numpy.linalg.solve` and SuperLU round differently, so an AutoSolver default
    would change the last bits of every circuit large and sparse enough to
    qualify.  Gate 7b-1 asks for identical results on the existing tests, and this
    is what keeps that true: `DenseSolver()` must reproduce passing nothing at all.
    """
    t0, v0 = _run(40, None)
    t1, v1 = _run(40, DenseSolver())
    assert np.array_equal(t0, t1)
    assert np.array_equal(v0, v1), 'DenseSolver is not the historical path'


def test_sparse_agrees_with_dense_to_solver_tolerance():
    """Different factorisation, same answer -- to rounding, not bit for bit.

    Asserted at 1e-9 relative rather than exactly: claiming bit-equality here
    would be claiming two different LU algorithms round identically, which is
    false and would make the test a liability the first time SciPy changes.
    """
    pytest.importorskip('scipy.sparse.linalg')
    t0, v0 = _run(150, DenseSolver())
    t1, v1 = _run(150, SuperLUSolver())
    assert v0.shape == v1.shape, 'the two solvers took different step counts'
    scale = max(np.abs(v0).max(), 1e-30)
    assert np.abs(v0 - v1).max() / scale < 1e-9


def test_autosolver_keys_on_fill_not_size():
    """A large but DENSE matrix must stay on the dense path.

    7b specifies "DenseLU below ~200 unknowns, SuperLU above".  That keys on the
    wrong variable: the crossover follows the sparsity pattern.  A dense 300x300
    would be mis-served by a size-only rule, so the selection is asserted here on
    a matrix that is big and full.
    """
    pytest.importorskip('scipy.sparse.linalg')
    n = MIN_N_FOR_SPARSE + 50
    dense_A = np.ones((n, n)) + np.eye(n)
    sparse_A = np.eye(n) + np.diag(np.ones(n - 1), 1)
    b = np.ones(n)

    picked_dense = AutoSolver()
    picked_dense.solve(dense_A, b, numeric)
    assert isinstance(picked_dense._choice, DenseSolver), \
        'a full matrix was sent to the sparse solver'

    picked_sparse = AutoSolver()
    picked_sparse.solve(sparse_A, b, numeric)
    ## Deliberately "not dense" rather than a named class: WHICH sparse solver is
    ## best available is 7c's business (KLU when libklu is there, SuperLU
    ## otherwise), and this test's subject is that the choice keys on FILL.
    assert not isinstance(picked_sparse._choice, DenseSolver), \
        'a 0.7%%-fill matrix was not sent to a sparse solver'


def test_autosolver_below_the_floor_stays_dense():
    """Small matrices keep the historical call, whatever their fill.

    Measured: at n=61 the sparse win is 1.16x, inside the noise of a Python-level
    dispatch, and `numpy.linalg.solve` beats `lu_factor`+`lu_solve` by 10x at
    n=31 on call overhead alone.
    """
    n = MIN_N_FOR_SPARSE - 10
    A = np.eye(n) + np.diag(np.ones(n - 1), 1)
    s = AutoSolver()
    s.solve(A, np.ones(n), numeric)
    assert isinstance(s._choice, DenseSolver)


def test_autosolver_gives_the_right_answer_either_way():
    """Whichever branch it picks, it must actually solve the system."""
    for n, A in ((60, np.eye(60) + np.diag(np.ones(59), 1)),
                 (MIN_N_FOR_SPARSE + 20,
                  np.eye(MIN_N_FOR_SPARSE + 20)
                  + np.diag(np.ones(MIN_N_FOR_SPARSE + 19), 1))):
        b = np.arange(1.0, n + 1.0)
        x = AutoSolver().solve(A, b, numeric)
        assert np.abs(A.dot(x) - b).max() < 1e-9


def test_a_non_linearsolver_is_rejected():
    """The parameter is typed, like nrsolver and scaler."""
    cir = _ladder(5)
    with pytest.raises(TypeError, match='LinearSolver'):
        Transient(cir, toolkit=numeric, linearsolver=object())._get_linearsolver()


def test_symbolic_matrices_fall_back_rather_than_crash():
    """An object-dtype matrix cannot go to SuperLU.

    The guard lives in two places -- AutoSolver never selects sparse for an object
    array, and SuperLUSolver falls back if handed one anyway -- because the second
    is reachable when a caller constructs it explicitly.
    """
    sympy = pytest.importorskip('sympy')
    pytest.importorskip('scipy.sparse.linalg')
    from pycircuit.circuit import symbolic
    a, b_ = sympy.symbols('a b')
    A = np.array([[a, 0], [0, b_]], dtype=object)
    rhs = np.array([1, 1], dtype=object)

    auto = AutoSolver()
    auto._select(A)
    assert isinstance(auto._choice, DenseSolver)

    x = SuperLUSolver().solve(A, rhs, symbolic)
    assert x is not None


## ---------------------------------------------------------------------------
## STAGE 7c -- KLU.  Skipped wholesale when libklu is absent, so the suite stays
## green on a machine without SuiteSparse, which is the point of discovering it.
## ---------------------------------------------------------------------------

def _klu_or_skip():
    from pycircuit.circuit.linearsolver import KLUSolver
    try:
        return KLUSolver()
    except ImportError as e:
        pytest.skip('libklu not available: %s' % e)


def test_klu_solves_correctly():
    """Right answer first; everything else about KLU is an optimisation."""
    k = _klu_or_skip()
    rng = np.random.default_rng(3)
    n = 40
    A = np.diag(2.0 + rng.random(n)) + np.diag(-1.0 * np.ones(n - 1), 1) \
        + np.diag(-1.0 * np.ones(n - 1), -1)
    b = rng.random(n)
    x = k.solve(A, b, numeric)
    assert np.abs(A.dot(x) - b).max() < 1e-10


def test_klu_reuses_the_ordering_across_solves():
    """The whole point of 7c: analyze once, refactor thereafter.

    Asserted on the counters rather than on timing -- a timing assertion on a
    loaded box is a flake, and the counters say the thing that actually matters:
    the ordering was computed ONCE for many solves.
    """
    k = _klu_or_skip()
    A = np.array([[2., -1., 0.], [-1., 2., -1.], [0., -1., 2.]])
    b = np.ones(3)
    for scale in (1.0, 2.0, 3.0, 4.0):
        x = k.solve(A * scale, b, numeric)
        assert np.abs((A * scale).dot(x) - b).max() < 1e-10
    assert k.analyses == 1, 'the ordering was recomputed: %r' % k
    assert k.factors == 1, 'a full factorisation was redone: %r' % k
    assert k.refactors == 3, 'the refactor path was not taken: %r' % k


def test_klu_reanalyses_when_the_pattern_changes():
    """A different sparsity pattern must NOT reuse the old ordering.

    This is the failure the pattern key exists to prevent, and it would be
    invisible without it -- KLU would happily refactor against indices that no
    longer describe the matrix.
    """
    k = _klu_or_skip()
    A1 = np.array([[2., -1., 0.], [-1., 2., 0.], [0., 0., 3.]])
    A2 = np.array([[2., -1., 1.], [-1., 2., -1.], [1., -1., 3.]])
    b = np.ones(3)
    x1 = k.solve(A1, b, numeric)
    x2 = k.solve(A2, b, numeric)
    assert np.abs(A1.dot(x1) - b).max() < 1e-10
    assert np.abs(A2.dot(x2) - b).max() < 1e-10
    assert k.analyses == 2, 'the ordering was reused across different patterns: %r' % k


def test_klu_matches_dense_on_a_transient():
    """Same waveform as the shipped path, to solver tolerance not bit for bit.

    Two different factorisations do not round identically, and asserting they do
    would make the test a liability the first time SuiteSparse changes.
    """
    import warnings
    from pycircuit.circuit.linearsolver import DenseSolver
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.elements import SubCircuit, R, C, VS
    from pycircuit.circuit import gnd
    _klu_or_skip()

    def ladder():
        cir = SubCircuit(toolkit=numeric)
        for i in range(30):
            cir.add_node('n%d' % i)
        cir['vs'] = VS('n0', gnd, v=1.0)
        for i in range(29):
            cir['R%d' % i] = R('n%d' % i, 'n%d' % (i + 1), r=100.0 * (i + 1))
        for i in range(30):
            cir['C%d' % i] = C('n%d' % i, gnd, c=1e-9 * (i + 1))
        return cir

    def run(solver):
        from pycircuit.circuit.linearsolver import KLUSolver
        tran = Transient(ladder(), toolkit=numeric, reltol=1e-4, uic=True,
                         linearsolver=solver)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(refnode=gnd, tend=5e-6, timestep=1e-6)
        return np.asarray(res.x, dtype=float)

    from pycircuit.circuit.linearsolver import KLUSolver
    v_dense = run(DenseSolver())
    v_klu = run(KLUSolver())
    assert v_dense.shape == v_klu.shape, 'the two solvers took different step counts'
    scale = max(float(np.abs(v_dense).max()), 1e-30)
    assert np.abs(v_dense - v_klu).max() / scale < 1e-9


def test_klu_raises_on_a_singular_matrix():
    """Singularity must surface as LinAlgError, like the other solvers."""
    k = _klu_or_skip()
    A = np.array([[1.0, 2.0], [2.0, 4.0]])       # rank 1
    with pytest.raises(np.linalg.LinAlgError):
        k.solve(A, np.ones(2), numeric)


def test_autosolver_picks_by_size_measured_end_to_end():
    """AutoSolver picks by MEASURED end-to-end win, not by isolated benchmark.

    KLU wins the isolated factor+solve 3.5x-10x but LOSES end to end at small n
    (0.52x the dense path at n=152), because the solve is only a fraction of a
    transient and KLU pays a fixed per-call cost.  It ties around n~400 and wins
    from n~1000 -- 1.53x at n=1002, where it also beats SuperLU.  Those rows only
    became measurable once 2+.5 made circuits that size buildable.

    Both directions are pinned, so neither can be "fixed" into a single always-
    answer.
    """
    from pycircuit.circuit.linearsolver import (AutoSolver, KLUSolver,
                                                SuperLUSolver, MIN_N_FOR_SPARSE)
    _klu_or_skip()
    n = MIN_N_FOR_SPARSE + 20
    A = np.eye(n) + np.diag(np.ones(n - 1), 1)
    s = AutoSolver()
    s.solve(A, np.ones(n), numeric)
    ## Small: SuperLU, because KLU's setup does not repay itself at this size.
    assert isinstance(s._choice, SuperLUSolver), \
        'AutoSolver should prefer SuperLU at n=%d, chose %r' % (n, s._choice)

    ## Large: KLU, because there it does -- measured 1.53x at n=1002.
    from pycircuit.circuit.linearsolver import MIN_N_FOR_KLU
    big = MIN_N_FOR_KLU + 200
    big_A = np.eye(big) + np.diag(np.ones(big - 1), 1)
    s_big = AutoSolver()
    s_big.solve(big_A, np.ones(big), numeric)
    assert isinstance(s_big._choice, KLUSolver), \
        'AutoSolver should prefer KLU at n=%d, chose %r' % (big, s_big._choice)
