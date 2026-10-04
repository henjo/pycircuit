"""Nothing reads a fixed multistep step's Jacobian (`_SteppingLoop.j_unread`,
speed round 9, stage 2): the converged point is evaluated without G -- by
the C Newton's fold and by the evaluate core's 'c' -- leaving the answers,
statistics, warnings and step state of the evaluation with G; every reader
of J (the adaptive error test, a direct caller of `solve_timestep`, a
stand-in for a piece of the loop) gets its J."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_core, _tran_newton_c, circuit
from pycircuit.circuit import transient as T
from pycircuit.circuit.integrator import TrapezoidalIntegrator
from pycircuit.circuit.tests.test_newton_c import _gp_chain, _mos_chain
from pycircuit.circuit.transient import Transient


def _run(build, skip, newton=True, cls=Transient, **kw):
    """`(x bytes, time bytes, statistics, warnings, step state bytes, path
    counts)` of one transient, the skip and the C Newton as asked."""
    make = kw.pop('make', {})
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(_tran_core, 'SKIP_UNREAD_J', skip)
        mp.setattr(_tran_newton_c, 'ENABLED', newton)
        before = _paths.snapshot()
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            tr = cls(build(), toolkit=circuit.numeric, **make)
            res = tr.solve(**kw)
        d = _paths.since(before)
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    state = [tr._q_cache[1], tr._C_cache[1], tr._iq, tr._Geq, tr._Cmat,
             *tr._qlast, *tr._iqlast]
    state = b''.join(np.asarray(a, float).tobytes() for a in state if a is not None)
    assert '_j_unread' not in tr.__dict__
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), state, d)


def _same(build, newton=True, **kw):
    a = _run(build, False, newton, **dict(kw))
    b = _run(build, True, newton, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert a[4] == b[4], 'the step state moved'
    assert a[5].get('newton_c:j_unread', 0) == 0 and a[5].get('core.c:served', 0) == 0, a[5]
    return a, b


CASES = {
    'mos': (_mos_chain, {}),
    'gp-chord': (_gp_chain, {'chord_jacobian': True}),
    'trap-capacitive': (lambda: _mos_chain(cap=1e-13), {'integrator': TrapezoidalIntegrator()}),
}


@pytest.mark.parametrize('case', list(CASES))
def test_the_c_newton_folds_a_fixed_step_without_g(case):
    if _tran_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_tran_newton_c.STATUS}')
    build, make = CASES[case]
    _a, b = _same(build, tend=1e-6, timestep=2e-8, fixed_timestep=True, make=make)
    served = b[5].get('newton_c:served', 0)
    assert served >= 40 and b[5].get('newton_c:j_unread', 0) == served, b[5]


@pytest.mark.parametrize('case', list(CASES))
def test_the_core_evaluates_a_fixed_step_without_g(case):
    """The C Newton off: `_newton` in Python, the converged point through
    the core's 'c'."""
    build, make = CASES[case]
    _a, b = _same(build, newton=False, tend=1e-6, timestep=2e-8, fixed_timestep=True,
                  make=make)
    if b[5].get('core.fj:served', 0) == 0:
        pytest.skip('the core does not serve this chain here')
    assert b[5].get('core.c:served', 0) >= 40 and b[5].get('core.j:served', 0) == 0, b[5]


def test_an_adaptive_run_reads_its_j(monkeypatch):
    """The error test reads J: every attempt hands it a matrix."""
    seen = []
    real = T._LMMSteps.judge

    def judging(self, X, x_new, h, J, clamped):
        seen.append(type(J))
        return real(self, X, x_new, h, J, clamped)
    monkeypatch.setattr(T._LMMSteps, 'judge', judging)
    _a, b = _same(_mos_chain, tend=1e-6, timestep=2e-8)
    assert len(seen) > 40 and set(seen) == {np.ndarray}
    assert b[5].get('newton_c:j_unread', 0) == 0 and b[5].get('core.c:served', 0) == 0


def test_a_stand_in_for_a_piece_of_the_loop_gets_j(monkeypatch):
    """`judge` patched on the loop's class, `solve_timestep` shadowed on
    the instance, a subclass of the transient: each may read J, and each
    run evaluates it."""
    seen = []
    real = T._SteppingLoop.judge

    def judging(self):
        seen.append(type(self.J))
        return real(self)
    monkeypatch.setattr(T._SteppingLoop, 'judge', judging)
    _a, b = _same(_mos_chain, tend=4e-7, timestep=2e-8, fixed_timestep=True)
    assert len(seen) >= 20 and set(seen) == {np.ndarray}
    assert b[5].get('newton_c:j_unread', 0) == 0 and b[5].get('core.c:served', 0) == 0
    monkeypatch.undo()

    class Spied(Transient):
        pass
    b = _run(_mos_chain, True, cls=Spied, tend=4e-7, timestep=2e-8, fixed_timestep=True)
    assert b[5].get('core.c:served', 0) == 0 and b[5].get('newton_c:j_unread', 0) == 0, b[5]

    seen.clear()
    tr = Transient(_mos_chain(), toolkit=circuit.numeric)
    orig = tr.solve_timestep

    def spied(*a, **k):
        out = orig(*a, **k)
        seen.append(type(out[2]))
        return out
    tr.solve_timestep = spied
    before = _paths.snapshot()
    tr.solve(tend=4e-7, timestep=2e-8, fixed_timestep=True)
    d = _paths.since(before)
    assert len(seen) >= 20 and set(seen) == {np.ndarray}
    assert d.get('newton_c:j_unread', 0) == 0 and d.get('core.c:served', 0) == 0, d


def test_a_direct_step_after_a_fixed_run_returns_its_j():
    """The mark lives for the loop alone: `solve_timestep` called after
    the run (as the shooting walks call it) returns the Jacobian."""
    tr = Transient(_mos_chain(), toolkit=circuit.numeric)
    res = tr.solve(tend=4e-7, timestep=2e-8, fixed_timestep=True)
    x = np.ascontiguousarray(np.asarray(res.x, float)[:, -1])
    tr._dt = 2e-8
    out = tr.solve_timestep(x, 4.2e-7)
    assert type(out[2]) is np.ndarray and out[2].shape == (x.shape[0], x.shape[0])


def test_the_switch_off_folds_g():
    if _tran_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_tran_newton_c.STATUS}')
    a = _run(_mos_chain, False, tend=4e-7, timestep=2e-8, fixed_timestep=True)
    assert a[5].get('newton_c:served', 0) >= 15
    assert a[5].get('newton_c:j_unread', 0) == 0 and a[5].get('core.c:served', 0) == 0
