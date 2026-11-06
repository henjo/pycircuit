"""The step Newton lets go of its attempt (`_StepNewton._newton`, roadmap
265).  Armed by the early abort, the Newton judges its evaluations through
`judged`: it kept every Jacobian an attempt evaluated (it reads the seed's
and the count), and the chord's residual judges through it while it reads
the chord -- the chord, its residual and `judged` held each other, so each
armed attempt's n x n arrays and the chord's numeric waited for the cycle
collector, which seldom runs in a solve: +850 MB a Gear-2 hand-back at
4 000 unknowns, a 1 us run's peak 7.6 GB.  Now `judged` keeps the seed's J
alone and the chord's residual is dropped once the solve is done: the
attempt's arrays go when it returns, every result the same."""
import gc
import warnings
import weakref

import pytest

from pycircuit.circuit import (
    _paths,
    _tran_loop_c,
    _tran_lte_c,
    _tran_newton,
    _tran_newton_c,
    circuit,
)
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.nrsolver import ChordNewton
from pycircuit.circuit.transient import Transient


def _chain(n=12):
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


@pytest.fixture
def python_newton(monkeypatch):
    """Every attempt on the Python loop and the Python Newton -- where the
    early abort arms `_newton` itself."""
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', False)
    monkeypatch.setattr(_tran_newton_c, 'ENABLED', False)
    monkeypatch.setattr(_tran_lte_c, 'ENABLED', False)


def _run(**kw):
    tr = Transient(_chain(), toolkit=circuit.numeric, **kw)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr.solve(tend=1e-7, timestep=1e-9)
    return _paths.since(before)


def test_an_armed_chord_attempt_leaves_nothing_to_the_cycle_collector(python_newton):
    """The collector held off through a run of armed chord attempts: after
    it, the unreachable objects hold no chord and no arm."""
    gc.collect()
    was = gc.isenabled()
    gc.disable()
    try:
        counts = _run(chord_jacobian=True)
        gc.set_debug(gc.DEBUG_SAVEALL)
        gc.collect()
    finally:
        gc.set_debug(0)
        if was:
            gc.enable()
    left = list(gc.garbage)
    gc.garbage.clear()
    assert counts['chord:chosen'] > 20 and counts['early:checks'] > 20
    held = [o for o in left if isinstance(o, (ChordNewton, _tran_newton._EarlyArm))]
    assert not held, f'{len(held)} chords and arms left to the cycle collector'


def test_an_armed_full_newton_keeps_no_jacobian_past_the_seeds(python_newton, monkeypatch):
    """The full Newton (no chord) evaluates J every iteration: at each
    evaluation, no earlier one's J but the first's is alive -- `judged`
    holds its own copy of the seed's."""
    real = _tran_newton._StepNewton._newton
    longest, alive = [0], []

    def spy(self, func, x0, residual=None):
        refs = []

        def f(x):
            alive.append(sum(r() is not None for r in refs[1:]))
            F, J = func(x)
            refs.append(weakref.ref(J))
            longest[0] = max(longest[0], len(refs))
            return F, J
        return real(self, f, x0, residual)
    monkeypatch.setattr(_tran_newton._StepNewton, '_newton', spy)
    counts = _run(chord_jacobian=False)
    assert counts['early:checks'] > 20 and longest[0] >= 3
    assert max(alive) == 0
