"""The C Newton serves a run's first attempt (`AutoSolver.choose_on_structure`,
`_tran_newton_c._first_choice`; roadmap 269).  The AutoSolver chose its
solver on the first matrix it saw, so the first attempt was the Python's --
its Newton on dense n x n state (J, its reductions, |J|) -- and the C
Newton declined it (`newton_c:linsolver`).  Where no value of that matrix
can change the choice -- below the floor (the size alone: dense), or a
structure of the matrix's size whose own fill is within `max_fill` (a
matrix whose nonzeros lie in it fills no more) -- the choice is made before
it, the StructuralKLU analysed as its first factor analyses it: the run's
results and its KLU counts are the Python's first attempt's."""
import warnings

import numpy as np

from pycircuit.circuit import _paths, _tran_newton_c, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.linearsolver import AutoSolver, StructuralKLU


def _structure(m, extra=0):
    """A tridiagonal m x m structure (`(m, rows, cols)`), `extra` more
    entries on the first row."""
    r = np.r_[np.arange(m), np.arange(m - 1), np.arange(1, m), np.zeros(extra, dtype=int)]
    c = np.r_[np.arange(m), np.arange(1, m), np.arange(m - 1), np.arange(2, 2 + extra)]
    return m, r, c


def test_the_choice_on_the_structure_is_the_first_matrixs():
    """Above the floor, a structure within `max_fill`: the StructuralKLU on
    it, counted once -- the class, the pattern and the count `_select` makes
    on a matrix whose nonzeros lie in the structure; below the floor dense,
    counted as `_select` counts it; a choice made is returned again."""
    m = 300
    st = _structure(m)
    a = AutoSolver(structure=lambda: st)
    b = AutoSolver(structure=lambda: st)
    J = np.zeros((m, m))
    J[st[1], st[2]] = np.random.default_rng(1).standard_normal(len(st[1]))
    before = _paths.snapshot()
    ca = a.choose_on_structure(m)
    after_a = _paths.since(before)
    cb = b._select(J)
    assert type(ca) is type(cb) is StructuralKLU
    assert np.array_equal(ca.take, cb.take) and np.array_equal(ca.Ap, cb.Ap) and ca.nnz == cb.nnz
    assert after_a.get('solver.auto:StructuralKLU') == 1
    assert a.choose_on_structure(m) is ca and a._select(J) is ca
    small = AutoSolver(structure=lambda: _structure(40))
    before = _paths.snapshot()
    assert small.choose_on_structure(40) is small._dense
    assert _paths.since(before).get('solver.auto:dense_small') == 1


def test_where_a_value_could_decide_nothing_is_chosen():
    """None, and nothing chosen, where the first matrix's values could
    change the choice: the structure's fill over `max_fill`, its size not
    the matrix's, no structure above the floor."""
    m = 300
    dense_st = _structure(m, extra=0)
    full = (m, np.repeat(np.arange(m), m), np.tile(np.arange(m), m))
    for ls in (AutoSolver(structure=lambda: full), AutoSolver(structure=lambda: dense_st, max_fill=0.001),
               AutoSolver(structure=lambda: _structure(m + 1)), AutoSolver()):
        assert ls.choose_on_structure(m) is None and ls._choice is None


def _chain(n=300):
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
    return c


def _run(monkeypatch, first):
    from pycircuit.circuit.transient import Transient
    if not first:
        monkeypatch.setattr(_tran_newton_c, '_first_choice', lambda ls, m: None)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(), toolkit=circuit.numeric)
        res = tr.solve(tend=6e-8, timestep=1e-9)
    st, ch = tr.statistics, tr._get_linearsolver()._choice
    run = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps, st.newton_iterations, st.chord_fallbacks)
    klu = tuple(getattr(ch, k) for k in ('analyses', 'factors', 'refactors', 'chord_factors',
                                          'chord_refactors', 'pivot_fallbacks', 'residual_fallbacks'))
    return run, (type(ch).__name__, klu), _paths.since(before)


def test_the_c_newton_serves_the_first_attempt_as_the_python_makes_it(monkeypatch):
    """A 302-unknown chain (a StructuralKLU): the C Newton serves every
    attempt, the first included -- no `newton_c:linsolver` decline -- and the
    run, its statistics and its KLU counts are those of the run whose first
    attempt is the Python's."""
    got, gklu, counts = _run(monkeypatch, first=True)
    with monkeypatch.context() as m:
        want, wklu, wcounts = _run(m, first=False)
    assert wcounts.get('newton_c:linsolver') == 1
    assert 'newton_c:linsolver' not in counts
    assert counts['newton_c:served'] == wcounts['newton_c:served'] + 1
    assert got == want and gklu == wklu and gklu[0] == 'StructuralKLU'
