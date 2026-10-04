"""The multistep Newton solve in C (`_tran_newton_c`, speed round 8, stage
4): bit for bit `_newton` and `jacobian_only`, its declines before anything
is touched, its bails without a trace."""
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _paths, _tran_newton_c, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import (
    VS,
    C,
    Diode,
    R,
    SubCircuit,
    VSin,
    gnd,
)
from pycircuit.circuit.integrator import (
    EulerIntegrator,
    Gear2Integrator,
    TrapezoidalIntegrator,
)
from pycircuit.circuit.transient import Transient


def _mos_chain(n=4, va=0.6, rl=5e3, cap=None):
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vg'] = VSin('g0', gnd, v=0.9, va=va, freq=1e6)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=rl)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', f'g{k}' if k == 0 else f'd{k - 1}', gnd, gnd)
        if cap:
            c[f'c{k}'] = C(f'd{k}', gnd, c=cap)
    return c


def _gp_chain(n=3, va=0.05):
    c = SubCircuit()
    c.add_node('vcc')
    c['vcc'] = VS('vcc', gnd, v=3.0)
    c.add_node('in')
    c['vin'] = VSin('in', gnd, v=0.75, va=va, freq=1e6)
    prev = 'in'
    for k in range(n):
        c.add_node(f'b{k}')
        c.add_node(f'c{k}')
        c[f'rb{k}'] = R(prev, f'b{k}', r=1e4)
        c[f'rc{k}'] = R('vcc', f'c{k}', r=1e3)
        c[f'Q{k}'] = eh.GummelPoonNpnHdl(f'c{k}', f'b{k}', gnd)
        prev = f'c{k}'
    return c


def _run(build, on, **kw):
    """`(x bytes, time bytes, statistics, warnings, path counts)` of one
    transient with the C solve on or off."""
    old = _tran_newton_c.ENABLED
    _tran_newton_c.ENABLED = on
    before = _paths.snapshot()
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            ## (a finalizer's warning is not the run's: a file another test left in
            ## a reference cycle closes whenever the collector runs -- inside either
            ## run, 2026-10-04)
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build(), toolkit=circuit.numeric, **kw.pop('make', {}))
            res = tr.solve(**kw)
    finally:
        _tran_newton_c.ENABLED = old
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__
           if 'seconds' not in k} if hasattr(tr.statistics, '__slots__') else \
        {k: v for k, v in vars(tr.statistics).items() if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), _paths.since(before))


def _same(build, served_min=1, **kw):
    a = _run(build, False, **dict(kw))
    b = _run(build, True, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert b[4].get('newton_c:served', 0) >= served_min, b[4]
    return a, b


@pytest.mark.parametrize('fixed', [True, False])
def test_the_plain_newton_is_the_python_one(fixed):
    """MOS level 1 (no chord): fixed steps and the default adaptive gear."""
    kw = {'tend': 2e-6, 'timestep': 2e-8}
    if fixed:
        kw['fixed_timestep'] = True
    _same(_mos_chain, **kw)


@pytest.mark.parametrize('fixed', [True, False])
def test_the_chord_is_the_python_one(fixed):
    """Gummel-Poon with the chord (SciPy's LU, called from C)."""
    kw = {'tend': 2e-6, 'timestep': 2e-8, 'make': {'chord_jacobian': True}}
    if fixed:
        kw['fixed_timestep'] = True
    a, b = _same(_gp_chain, **kw)
    assert b[4].get('chord:chosen', 0) == a[4].get('chord:chosen', 0) > 0


@pytest.mark.parametrize('integ', [Gear2Integrator, EulerIntegrator, TrapezoidalIntegrator])
def test_every_companion_and_a_capacitive_chain(integ):
    """The core's companion enum, on a chain with capacitors (the branch
    screen reads the converged point's `C` from the C solve)."""
    _same(lambda: _mos_chain(cap=1e-13), tend=1.5e-6, timestep=2e-8,
          make={'integrator': integ()})


@settings(max_examples=12, deadline=None)
@given(n=st.integers(1, 5), va=st.floats(0.05, 0.9), rl=st.floats(1e3, 2e4),
       cap=st.sampled_from([None, 1e-14, 1e-12]))
def test_drawn_chains_are_the_python_ones(n, va, rl, cap):
    _same(lambda: _mos_chain(n, va, rl, cap), served_min=0, tend=1e-6, timestep=2e-8)


def test_a_bail_leaves_no_trace(monkeypatch):
    """Every C solve runs to its end and then hands the step back (its
    entry wrapped to report `maxiter`): the Python Newton repeats each one
    -- the same answer, statistics, warnings, and every count but the
    bail's own (the source memo's entry and counts rolled back)."""
    drv = _tran_newton_c.driver()
    if drv is None:
        pytest.skip(f'the C solve is off: {_tran_newton_c.STATUS}')
    ffi, cfn, lapack = drv
    ## (once first, unmeasured: a process's first C build counts its
    ## `once:` keys in whichever run comes first)
    _run(_mos_chain, False, tend=6e-7, timestep=2e-8)
    a = _run(_mos_chain, False, tend=6e-7, timestep=2e-8)

    def bailing(s):
        cfn(s)
        return 11
    monkeypatch.setattr(_tran_newton_c, '_driver', (ffi, bailing, lapack))
    b = _run(_mos_chain, True, tend=6e-7, timestep=2e-8)
    assert a[:4] == b[:4]
    assert b[4].get('newton_c:served', 0) == 0 and b[4].get('newton_c:bail:maxiter', 0) > 10
    rest = {k: v for k, v in b[4].items() if not k.startswith('newton_c:')}
    assert rest == {k: v for k, v in a[4].items() if not k.startswith('newton_c:')}


def test_the_declines():
    """A stateful limiter (`Diode`) and a provided function: today's
    Newton, counted -- the circuit's reason once, kept with its plan."""
    def diode():
        c = _mos_chain(2)
        c['D'] = Diode('d0', gnd)
        return c
    _, b = _same(diode, served_min=0, tend=4e-7, timestep=2e-8)
    assert b[4].get('newton_c:stateful', 0) > 0 and b[4].get('newton_c:served', 0) == 0
    before = _paths.snapshot()
    Transient(_mos_chain(2), toolkit=circuit.numeric, uic=True).solve(
        tend=2e-7, timestep=2e-8, provided_function=lambda t: 0.0)
    d = _paths.since(before)
    assert d.get('newton_c:pf', 0) > 0 and d.get('newton_c:served', 0) == 0


def test_an_instance_shadow_of_the_newton_is_honoured():
    """A caller wrapping `_newton` on the instance expects its calls."""
    tr = Transient(_mos_chain(2), toolkit=circuit.numeric)
    calls = []
    real = tr._newton

    def wrapped(*a, **k):
        calls.append(1)
        return real(*a, **k)
    tr._newton = wrapped
    tr.solve(tend=2e-7, timestep=2e-8)
    assert len(calls) >= 10


def test_a_class_patch_and_a_subclass_take_the_python_newton(monkeypatch):
    """`_newton` patched on the class, or overridden by a subclass: today's
    Newton runs, every step, counted (the C stood in for both until
    2026-10-04)."""
    calls = []
    real = Transient._newton

    def counting(self, *a, **k):
        calls.append(1)
        return real(self, *a, **k)
    monkeypatch.setattr(Transient, '_newton', counting)
    before = _paths.snapshot()
    Transient(_mos_chain(2), toolkit=circuit.numeric).solve(tend=4e-7, timestep=2e-8,
                                                            fixed_timestep=True)
    d = _paths.since(before)
    assert len(calls) >= 15 and d.get('newton_c:served', 0) == 0, d
    assert d.get('newton_c:patched', 0) >= 15
    monkeypatch.undo()

    class Mine(Transient):
        def _newton(self, *a, **k):
            calls.append(2)
            return super()._newton(*a, **k)
    before = _paths.snapshot()
    Mine(_mos_chain(2), toolkit=circuit.numeric).solve(tend=4e-7, timestep=2e-8,
                                                       fixed_timestep=True)
    d = _paths.since(before)
    assert calls.count(2) >= 15 and d.get('newton_c:class', 0) >= 15, d


def test_the_c_solve_serves_the_multistep_step():
    if _tran_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_tran_newton_c.STATUS}')
    before = _paths.snapshot()
    Transient(_mos_chain(), toolkit=circuit.numeric).solve(tend=1e-6, timestep=2e-8,
                                                            fixed_timestep=True)
    d = _paths.since(before)
    assert d.get('newton_c:served', 0) >= 40
    assert d.get('newton_c:jfold', 0) == d.get('newton_c:served', 0)
