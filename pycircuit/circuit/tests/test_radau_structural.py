"""The Radau transform's complex factor on the structural pattern (roadmap
259; `_tran_radau.STRUCTURAL`): where the analysis solver chose
`StructuralKLU`, the complex factor's record is on its pattern
(`ComplexKLUSolver.prepare_structural`) -- explicit zeros among its values,
one analysis a run, refactors after, a refactor meeting a zero pivot
factored afresh -- where a pattern of nonzeros changed, and was analysed
and factored afresh, whenever a device's entries turned exactly zero.  The
Python step, the C fold, the C transform served from Python and the C loop
alike."""
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from pycircuit.circuit import (
    _paths,
    _tran_loop_c,
    _tran_radau,
    _tran_radau_c,
    _tran_radau_tc,
    _tran_tcloop_c,
    circuit,
)
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.linearsolver import ComplexKLUSolver, StructuralKLU
from pycircuit.circuit.transient import Transient

#: the paths a run may take: (the C fold, the C transform, the C loops)
PATHS = {'python': (False, False, False), 'fold': (True, False, False),
         'served': (True, True, False), 'loop': (True, True, True)}


def _chain(n=250):
    """A pulsed MOS chain with capacitors: at 250 stages 254 unknowns --
    `StructuralKLU` -- its devices' entries turning exactly zero at cutoff."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-8, tr=1e-8, tf=1e-8, pw=2.4e-7, per=5e-7)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
    return c


def _run(monkeypatch, path='loop', structural=True, tend=3e-7, cir=None, rt=None):
    fold, tc, loops = PATHS[path]
    monkeypatch.setattr(_tran_radau, 'STRUCTURAL', structural)
    monkeypatch.setattr(_tran_radau_tc, 'FOLD', fold)
    monkeypatch.setattr(_tran_radau_tc, 'ENABLED', tc)
    monkeypatch.setattr(_tran_radau_c, 'ENABLED', tc)
    monkeypatch.setattr(_tran_tcloop_c, 'ENABLED', loops)
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', loops)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(cir if cir is not None else _chain(), toolkit=circuit.numeric,
                       integrator=RadauIIA3Integrator())
        if rt is not None:
            tr.par.reltol = rt
        res = tr.solve(tend=tend, timestep=1e-8)
    return tr, res, _paths.since(before)


def _key(tr, res):
    """The run's bytes, statistics and both factors' counts."""
    st = tr.statistics
    zs = tr.__dict__.get('_radau_zsolver')
    sk = tr._get_linearsolver()._choice
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st.accepted_steps, st.rejected_steps, st.newton_iterations, type(sk).__name__,
            None if zs is None else (zs.analyses, zs.factors, zs.refactors, zs.residual_fallbacks,
                                     zs.pivot_fallbacks),
            *(getattr(sk, k, None) for k in ('factors', 'refactors', 'pivot_fallbacks',
                                             'residual_fallbacks')))


def test_a_sparse_run_analyses_its_complex_factor_once(monkeypatch):
    """On `StructuralKLU`'s pattern the complex factor is analysed once (and
    once before the structure's own analysis is made), refactored at every
    later step; the switch off makes no structural record."""
    tr, _res, d = _run(monkeypatch, path='python', tend=1e-6)
    zs = tr._radau_zsolver
    assert type(tr._get_linearsolver()._choice) is StructuralKLU
    assert zs.analyses <= 2 and zs.refactors > 50 and d["radau.frozen:structural"] > 50
    _t0, _r0, d0 = _run(monkeypatch, path='python', structural=False)
    assert d0.get('radau.frozen:structural', 0) == 0


@pytest.mark.parametrize('path', ['fold', 'served', 'loop'])
def test_every_path_is_the_python_one(monkeypatch, path):
    """The C fold, the C transform served from Python and the C loop make
    the Python's run on the structural pattern -- bytes, statistics, both
    factors' counts, the structural records -- each serving its share."""
    a = _run(monkeypatch, path='python')
    b = _run(monkeypatch, path=path)
    assert _key(*a[:2]) == _key(*b[:2])
    assert a[2]['radau.frozen:structural'] == b[2]['radau.frozen:structural']
    served = {'fold': 'radau.frozen:fold', 'served': 'radau_tc:served',
              'loop': 'tcloop_c:attempts'}[path]
    assert b[2].get(served, 0) > 20
    ## (on the structural pattern the fold meets no record's mask)
    assert b[2]['radau.frozen:fold'] > 20 and b[2].get('radau.frozen:fold.same', 0) == 0


def test_a_dense_run_keeps_its_patterns(monkeypatch):
    """Without `StructuralKLU` (a small chain: the dense solver) nothing is
    structural -- the run the one the switch off makes."""
    a = _run(monkeypatch, cir=_chain(8))
    b = _run(monkeypatch, cir=_chain(8), structural=False)
    assert _key(*a[:2]) == _key(*b[:2])
    assert a[2].get('radau.frozen:structural', 0) == 0


def _sk(n):
    """A structure as `prepare_structural` reads one: every entry but
    (0, n - 1), CSC column by column, rows ascending, `take` each entry's
    place in the row-major matrix."""
    on = np.ones((n, n), dtype=bool)
    on[0, n - 1] = False
    Ap = np.zeros(n + 1, dtype=np.int32)
    Ap[1:] = np.cumsum(on.sum(axis=0))
    Ai = np.ascontiguousarray(on.T.nonzero()[1], dtype=np.int32)
    j = np.repeat(np.arange(n), on.sum(axis=0))
    return SimpleNamespace(m=n, Ap=Ap, Ai=Ai, take=np.asarray(Ai, dtype=np.int64) * n + j)


def test_the_structural_record_and_its_refusal():
    """The record is the structure's -- its `Ap`, `Ai`, a fixed key, the
    entries at its positions with explicit zeros -- and None for a matrix
    with a nonzero off it, or of another size or type."""
    zs = ComplexKLUSolver()
    sk = _sk(4)
    A = np.zeros((4, 4), dtype=complex)
    A[np.arange(4), np.arange(4)] = 2.0 + 1.0j
    A[2, 1] = 0.5
    rec = zs.prepare_structural(A, sk)
    assert rec[2] is zs._skey[3] and rec[3] is zs._skey[4] and rec[1] == 4
    assert rec[5][0] == ComplexKLUSolver.STRUCTURAL
    data = rec[4].view(complex)
    assert data.shape == (15,) and np.count_nonzero(data) == 5
    assert np.array_equal(data, A.ravel()[sk.take])
    assert zs.prepare_structural(A, sk)[5] is rec[5]
    B = A.copy()
    B[0, 3] = 1.0
    assert zs.prepare_structural(B, sk) is None
    assert zs.prepare_structural(A[:3, :3].copy(), sk) is None
    assert zs.prepare_structural(A.real.copy(), sk) is None


def test_a_refactor_meeting_a_zero_pivot_factors_afresh():
    """On the structural pattern a refactor whose held pivot became zero --
    a tridiagonal matrix whose middle pivot cancels, nonsingular all the
    same -- is factored afresh (`pivot_fallbacks`), the solution numpy's;
    on a pattern of nonzeros the same refactor still raises."""
    sk = SimpleNamespace(m=3, Ap=np.array([0, 2, 5, 7], dtype=np.int32),
                         Ai=np.array([0, 1, 0, 1, 2, 1, 2], dtype=np.int32),
                         take=np.array([0, 3, 1, 4, 7, 5, 8], dtype=np.int64))
    b = np.array([1.0 + 2.0j, -0.5j, 0.25])
    A = np.array([[2.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 2.0]], dtype=complex)
    B = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]], dtype=complex)
    zs = ComplexKLUSolver()
    zs.solve_prepared(zs.prepare_structural(A, sk), b)
    x = zs.solve_prepared(zs.prepare_structural(B, sk), b)
    assert zs.pivot_fallbacks == 1 and zs.analyses == 1 and zs.factors == 2
    assert np.allclose(x, np.linalg.solve(B, b), rtol=1e-14, atol=1e-15)
    zn = ComplexKLUSolver()
    zn.solve(A, b)
    with pytest.raises(np.linalg.LinAlgError):
        zn.solve(B, b)


def _err(ref, res, nodes=('d0', 'd1', 'd2', 'd3', 'd4', 'd5')):
    out = []
    for nd in nodes:
        wr, w = ref.v(nd), res.v(nd)
        tr_, yr = np.asarray(wr.x[0], float), np.asarray(wr.y, float)
        y = np.interp(tr_, np.asarray(w.x[0], float), np.asarray(w.y, float))
        out.append(np.max(np.abs(y - yr)) / (yr.max() - yr.min()))
    return max(out)


def test_a_structural_run_is_within_its_tolerance(monkeypatch):
    """Against a reltol/100 reference over a microsecond, the switch off and
    on: the answers move (a round-off change: another ordering, refactors on
    held pivots), and on the error of the first six stages is at most 1.25x
    the off run's + reltol."""
    _t, ref, _d = _run(monkeypatch, structural=False, rt=1e-6, tend=1e-6)
    t0, off, _d = _run(monkeypatch, structural=False, tend=1e-6)
    t1, on, d = _run(monkeypatch, tend=1e-6)
    assert d['radau.frozen:structural'] > 50 and _key(t0, off) != _key(t1, on)
    assert _err(ref, on) <= 1.25 * _err(ref, off) + t1.par.reltol
