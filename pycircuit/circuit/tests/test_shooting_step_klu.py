"""A fully implicit stage method's coupled stage block factored with KLU
on its structural pattern (`_InnerTransient._step_klu`) from
`SPARSE_BLOCK_MIN_N` unknowns; the answers the dense LU's to the
rounding."""
import numpy as np
import pytest

from pycircuit.circuit import _paths, circuit, gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import C, Diode, R, VSin
from pycircuit.circuit.linearsolver import DenseSolver, StructuralKLU
from pycircuit.circuit.shooting import PAC, PSS, _pss_inner
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._warnpolicy import quiet

#: sections enough for Radau's `3 m` block to reach the floor (`m =
#: sections + 2`)
SECTIONS = _pss_inner.SPARSE_BLOCK_MIN_N // 3


def _ladder(sections, diode=True):
    """A sine into an RC ladder of `sections`, a diode at its end:
    ``sections + 2`` unknowns."""
    c = SubCircuit()
    c['vs'] = VSin('n0', gnd, vac=1.0, va=1.0, freq=1e3)
    for k in range(sections):
        c[f'R{k}'] = R(f'n{k}', f'n{k + 1}', r=100.0)
        c[f'C{k}'] = C(f'n{k + 1}', gnd, c=1e-8)
    if diode:
        c['D'] = Diode(f'n{sections}', gnd)
    return c


def _pss(cir, klu, monkeypatch, npts=20, **kw):
    circuit.default_toolkit = circuit.numeric
    monkeypatch.setattr(_pss_inner, 'SHOOT_KLU', klu)
    pss = PSS(cir, method='radau', reltol=1e-10, **kw)
    with quiet(AccuracyWarning):
        res = pss.solve(period=1e-3, timestep=1e-3 / npts)
    assert pss.converged
    out = 'n{}'.format(len([e for e in cir.elements if e.startswith('R')]))
    return pss, np.asarray(res['tpss'].v(out).get_y(), dtype=float)


def test_a_large_stage_block_factors_on_its_structure_and_the_c_step_declines(monkeypatch):
    """Radau's `3 m` block from `SPARSE_BLOCK_MIN_N`, on a linear ladder
    the dense C step map (`_sens_c`) would otherwise take: the block is
    factored on its structure, and the dense map is not asked -- its LU is
    not that step's twin; the KLU step map (`_sensklu_c`) is.  (Without it
    the analysis's solver chose on the first matrix it saw, here SuperLU on
    the block.)"""
    out = {}
    for klu in (False, True):
        before = _paths.snapshot()
        pss, y = _pss(_ladder(SECTIONS, diode=False), klu, monkeypatch)
        out[klu] = (y, _paths.since(before), pss)
    (y0, d0, _), (y1, d1, pss1) = out[False], out[True]
    n = 3 * (pss1.cir.n - 1)
    assert n >= _pss_inner.SPARSE_BLOCK_MIN_N
    assert not d0.get('pss.sensklu:served')
    assert d1.get('pss.sensklu:served', 0) > 0 and not d1.get('pss.sens:served')
    klu = pss1._step_klus[n][1]
    assert type(klu) is StructuralKLU and klu.numerics > 0 and klu.foreign == 0
    assert np.allclose(y1, y0, rtol=1e-9, atol=1e-12)


def test_a_step_jacobian_and_a_small_block_stay_with_the_caller_s_solver(monkeypatch):
    monkeypatch.setattr(_pss_inner, 'SHOOT_KLU', True)
    big = PSS(_ladder(SECTIONS), method='radau')
    m = big.cir.n - 1
    assert big._step_klu(m) is None and big._step_klu(3 * m) is not None
    small = PSS(_ladder(SECTIONS - 3), method='radau')
    ms = small.cir.n - 1
    assert 3 * ms < _pss_inner.SPARSE_BLOCK_MIN_N and small._step_klu(3 * ms) is None
    own = PSS(_ladder(SECTIONS), method='radau', linearsolver=DenseSolver())
    assert own._step_klu(3 * m) is None


def test_a_factored_period_on_the_structure_answers_as_dense(monkeypatch):
    """The kept walk's block factors (`_factorise`) and their transposed
    solves -- PAC's adjoint row and its sideband response -- on the
    structure answer as without it."""
    got = {}
    for klu in (False, True):
        pss, _y = _pss(_ladder(SECTIONS), klu, monkeypatch)
        pac = PAC(pss.cir, toolkit=circuit.numeric)
        with quiet(AccuracyWarning):
            row, _info = pac.adjoint_transfer_row(pss, 300.0, f'n{SECTIONS}')
            res = pac.solve(pss, [300.0])
        got[klu] = (np.asarray(row), np.asarray(res.v(f'n{SECTIONS}', gnd)))
        if klu:
            assert pss._step_klus[3 * (pss.cir.n - 1)][1].numerics > 0
    for a, b in zip(got[True], got[False]):
        assert np.allclose(a, b, rtol=1e-8, atol=1e-14 * np.max(np.abs(b)))


def test_the_c_stage_newton_on_the_walk_s_klu_is_the_python_s(monkeypatch):
    """The walk's transient solves its stage block on the walk's KLU
    (`tr._stage_klu`, `_new_transient`): the C stage Newton makes
    StructuralKLU.solve's sequence on it, and a live stage step factors the
    C's last assembly (`_tran_radau_c.final_block`) -- the orbit and the
    monodromy the bytes of the Python Newton's and of `_stage_block`'s.
    (The linear ladder: the diode's limiter keeps a state, which the C
    declines; `test_radau_c` holds the KLU path on a nonlinear circuit.)"""
    from pycircuit.circuit import _tran_radau_c
    out = {}
    for on in (False, True):
        monkeypatch.setattr(_tran_radau_c, 'ENABLED', on)
        before = _paths.snapshot()
        pss, y = _pss(_ladder(SECTIONS, diode=False), True, monkeypatch)
        out[on] = (y.tobytes(), np.asarray(pss._monodromy).tobytes(), _paths.since(before))
        assert pss._tran.__dict__.get('_stage_klu') is pss._step_klus[3 * (pss.cir.n - 1)][1]
    assert out[True][:2] == out[False][:2]
    d1, d0 = out[True][2], out[False][2]
    assert d1.get('radau_c:on_klu', 0) > 0 and d1.get('radau_c:final', 0) > 0, d1
    assert not d0.get('radau_c:on_klu') and not d0.get('radau_c:final')


def test_a_shooting_walk_refactors_one_numeric_and_keeps_its_answer(monkeypatch):
    """A walk that keeps no step (the shooting Newton's monodromy) factors
    its first block and refactors one numeric for every later one
    (`StructuralKLU.walk`, roadmap 198): one fresh factor and `steps - 1`
    refactors a walk, none failing its check; a kept walk (PAC's factored
    period) factors each step on its own.  The orbit and the monodromy
    within the rounding of a fresh factor every step (the walk's numeric
    left out)."""
    walks = []
    orig = StructuralKLU.walk

    def walk(self):
        walks.append(orig(self))
        return walks[-1]
    monkeypatch.setattr(StructuralKLU, 'walk', walk)
    pss, y = _pss(_ladder(SECTIONS), True, monkeypatch)
    P = np.asarray(pss._monodromy)
    assert walks and all(w.factors == 1 and w.refactors == 19 and w.fallbacks == 0
                         for w in walks), [(w.factors, w.refactors, w.fallbacks) for w in walks]
    monkeypatch.setattr(StructuralKLU, 'walk', lambda self: None)
    pss0, y0 = _pss(_ladder(SECTIONS), True, monkeypatch)
    assert np.allclose(y, y0, rtol=1e-12, atol=1e-15)
    assert np.allclose(P, np.asarray(pss0._monodromy), rtol=1e-10, atol=1e-14)


def test_the_klu_step_map_in_c_is_the_python_one(monkeypatch):
    """The plain walk's step on the walk's KLU in one C call (`_sensklu_c`,
    roadmap 199): the orbit and the monodromy the bytes of the Python
    step's, every walk's numeric made as often (factors, refactors,
    fallbacks, generation) and the solver's `numerics` the same."""
    from pycircuit.circuit.shooting import _sensklu_c
    if _sensklu_c.driver() is None:
        pytest.skip(f'the KLU step map in C is off: {_sensklu_c.STATUS}')
    orig = StructuralKLU.walk
    out = {}
    for on in (False, True):
        monkeypatch.setattr(_sensklu_c, 'SENSKLU_C', on)
        walks = []

        def walk(self, walks=walks):
            walks.append(orig(self))
            return walks[-1]
        monkeypatch.setattr(StructuralKLU, 'walk', walk)
        before = _paths.snapshot()
        pss, y = _pss(_ladder(SECTIONS, diode=False), True, monkeypatch)
        sk = pss._step_klus[3 * (pss.cir.n - 1)][1]
        out[on] = (y.tobytes(), np.asarray(pss._monodromy).tobytes(),
                   [(w.factors, w.refactors, w.fallbacks, w.gen) for w in walks], sk.numerics,
                   _paths.since(before))
    assert out[True][:4] == out[False][:4]
    assert out[True][4].get('pss.sensklu:served', 0) > 30 and not out[False][4].get('pss.sensklu:served')


def _replay_runs(monkeypatch, alter=None):
    """The PSS and a PAC on it, with the converged replay as `alter` leaves
    `_stage_walk_reused` (None: as it is)."""
    from pycircuit.circuit.shooting._pss_replays import _FactoredReplays
    if alter is not None:
        monkeypatch.setattr(_FactoredReplays, '_stage_walk_reused', alter)
    before = _paths.snapshot()
    pss, y = _pss(_ladder(SECTIONS), True, monkeypatch)
    d = _paths.since(before)
    pac = PAC(pss.cir, toolkit=circuit.numeric)
    with quiet(AccuracyWarning):
        res = pac.solve(pss, [300.0])
    return y, np.asarray(res.v(f'n{SECTIONS}', gnd)), d


def test_the_converged_replay_is_the_shooting_s_last_walk(monkeypatch):
    """A stage map's converged replay takes the shooting's last walk -- the
    map at the point the solve returned, its period and grid -- instead of
    walking it again: the waveform the walked replay's to the rounding, and
    the factored period's answer (PAC, from the record) as walked."""
    y1, v1, d1 = _replay_runs(monkeypatch)
    y0, v0, d0 = _replay_runs(monkeypatch, lambda self, *a: None)
    assert d1.get('pss.replay:reused') == 1 and not d0.get('pss.replay:reused')
    assert np.allclose(y1, y0, rtol=1e-12, atol=1e-14 * np.max(np.abs(y0)))
    assert np.allclose(v1, v0, rtol=1e-10, atol=1e-14 * np.max(np.abs(v0)))


def test_a_last_walk_at_another_point_is_not_the_replay(monkeypatch):
    """The record of a walk at another point (its bytes), another period or
    grid is not taken: the replay walks, as without it."""
    from pycircuit.circuit.shooting._pss_replays import _FactoredReplays
    orig = _FactoredReplays._stage_walk_reused

    def elsewhere(self, x0, T, times, hs):
        last = self._last_stage_walk
        self._last_stage_walk = (bytes(len(last[0])),) + last[1:]
        return orig(self, x0, T, times, hs)
    y1, v1, d1 = _replay_runs(monkeypatch, elsewhere)
    monkeypatch.setattr(_FactoredReplays, '_stage_walk_reused', orig)
    y0, v0, _d0 = _replay_runs(monkeypatch, lambda self, *a: None)
    assert not d1.get('pss.replay:reused')
    assert y1.tobytes() == y0.tobytes() and v1.tobytes() == v0.tobytes()
