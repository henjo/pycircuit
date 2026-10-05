"""The source pass from the plan (`_stamp_plan.assemble_source`, speed round
9, stage 5; round 8's 2a): `u` and `dudt` call only the elements whose
source is not the default or a generated zero, in dict order, scattered as
the loop scatters -- bit for bit the loop (`SubCircuit.
_add_element_subvectors`), its `:called` counts with it, and the loop again
wherever a check fails."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _stamp_plan, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import IS, R, SubCircuit, VPulse, VSin, gnd
from pycircuit.circuit.tests.test_newton_c import _mos_chain
from pycircuit.circuit.tests.test_stamp_plan import mixed
from pycircuit.circuit.transient import Transient


def _pulsed():
    c = SubCircuit()
    c['vp'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-7, tr=1e-9, tf=1e-9, pw=4e-7, per=1e-6)
    c['vs'] = VSin('b', gnd, v=0.3, va=0.2, freq=2e6)
    c['ra'] = R('a', 'c', r=1e3)
    c['rb'] = R('b', 'c', r=2e3)
    c['is'] = IS('c', gnd, i=1e-4)
    c['m'] = eh.MosLevel1Hdl('c', 'b', gnd, gnd)
    return c


CIRCUITS = {'mixed': mixed, 'pulsed': _pulsed, 'mos': _mos_chain}
TIMES = (0.0, 1e-9, 1.05e-7, 3.3e-7, 9.99e-7, 2.5e-6)


def _pass(cir, m, t, analysis, on, monkeypatch):
    monkeypatch.setattr(_stamp_plan, 'SOURCE_PLAN', on)
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        warnings.simplefilter('ignore', ResourceWarning)
        v = getattr(cir, m)(t, defaultepar, analysis)
    d = _paths.since(before)
    return np.asarray(v).tobytes(), sorted(str(w.message) for w in W), d


@pytest.mark.parametrize('analysis', [None, 'tran', 'dc'])
@pytest.mark.parametrize('m', ['u', 'dudt'])
@pytest.mark.parametrize('name', list(CIRCUITS))
def test_the_plan_is_the_loop(name, m, analysis, monkeypatch):
    circuit.default_toolkit = circuit.numeric
    cir = CIRCUITS[name]()
    served = 0
    for t in TIMES:
        a = _pass(cir, m, t, analysis, False, monkeypatch)
        b = _pass(cir, m, t, analysis, True, monkeypatch)
        assert a[0] == b[0], (name, m, t, analysis)
        assert a[1] == b[1], 'the warnings moved'
        assert a[2].get(m + ':called', 0) == b[2].get(m + ':called', 0)
        served += b[2].get(f'src.{m}:served', 0)
    ## (a nested circuit's own pass is served from its own plan too)
    assert served >= len(TIMES), 'the plan did not serve'


def test_the_ac_pass_and_a_dtype_take_the_loop(monkeypatch):
    circuit.default_toolkit = circuit.numeric
    cir = _pulsed()
    before = _paths.snapshot()
    cir.u(0.0, defaultepar, 'ac')
    d = _paths.since(before)
    assert d.get('src.u:served', 0) == 0


def test_what_the_plan_declines_takes_the_loop(monkeypatch):
    """An instance shadow of `u` (a spy) set after the plan was built, or
    present when it is built, `u` patched on a class, the zero-source skip
    switched off: each takes the loop -- with its answer -- counted by its
    reason."""
    circuit.default_toolkit = circuit.numeric
    cir = _pulsed()
    ref = np.asarray(cir.u(2e-7, defaultepar, 'tran')).tobytes()

    def same(reason):
        before = _paths.snapshot()
        got = np.asarray(cir.u(2e-7, defaultepar, 'tran')).tobytes()
        d = _paths.since(before)
        assert d.get('src.u:' + reason, 0) == 1 and d.get('src.u:served', 0) == 0, d
        return got
    calls = []
    real = type(cir['vs']).u

    def spy(*a, **k):
        calls.append(1)
        return real(cir['vs'], *a, **k)
    cir['vs'].__dict__['u'] = spy
    assert same('element') == ref and calls
    _stamp_plan.invalidate(cir)
    calls.clear()
    assert same('shadowed') == ref and calls
    del cir['vs'].__dict__['u']
    _stamp_plan.invalidate(cir)
    cir.u(2e-7, defaultepar, 'tran')
    monkeypatch.setattr(type(cir['vs']), 'u', lambda self, *a, **k: real(self, *a, **k))
    assert same('patched') == ref
    monkeypatch.undo()
    monkeypatch.setattr(_stamp_plan._hdl_batch, 'SKIP_ZERO_SOURCE', False)
    assert same('gate') == ref
    monkeypatch.undo()
    before = _paths.snapshot()
    assert np.asarray(cir.u(2e-7, defaultepar, 'tran')).tobytes() == ref
    assert _paths.since(before).get('src.u:served', 0) == 1


def test_a_small_circuit_takes_the_loop(monkeypatch):
    """Under `SOURCE_MIN_ELEMENTS` the loop's visits cost less than the
    plan's checks: the plan is not asked (and serves from the threshold)."""
    circuit.default_toolkit = circuit.numeric
    small = _mos_chain(1)
    assert len(small.elements) < _stamp_plan.SOURCE_MIN_ELEMENTS
    a = _pass(small, 'u', 1e-7, 'tran', True, monkeypatch)
    assert not any(k.startswith('src.') for k in a[2]), a[2]
    assert a[0] == _pass(small, 'u', 1e-7, 'tran', False, monkeypatch)[0]
    monkeypatch.setattr(_stamp_plan, 'SOURCE_MIN_ELEMENTS', len(small.elements))
    assert _pass(small, 'u', 1e-7, 'tran', True, monkeypatch)[2].get('src.u:served') == 1


def test_an_element_added_rebuilds_the_plan():
    circuit.default_toolkit = circuit.numeric
    cir = _pulsed()
    cir.u(1e-7, defaultepar, 'tran')
    cir['vx'] = VSin('c', 'b', v=0.0, va=0.1, freq=3e6)
    before = _paths.snapshot()
    got = np.asarray(cir.u(1e-7, defaultepar, 'tran')).tobytes()
    assert _paths.since(before).get('src.u:served', 0) == 1
    _was = _stamp_plan.SOURCE_PLAN
    _stamp_plan.SOURCE_PLAN = False
    try:
        assert np.asarray(cir.u(1e-7, defaultepar, 'tran')).tobytes() == got
    finally:
        _stamp_plan.SOURCE_PLAN = _was


@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
def test_a_transient_is_the_same(fixed, monkeypatch):
    def run(on):
        monkeypatch.setattr(_stamp_plan, 'SOURCE_PLAN', on)
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(_pulsed(), toolkit=circuit.numeric)
            kw = {'tend': 1.5e-6, 'timestep': 1e-8}
            if fixed:
                kw['fixed_timestep'] = True
            res = tr.solve(**kw)
        st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__
               if 'seconds' not in k}
        return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values).tobytes(),
                st_, sorted(str(w.message) for w in W))
    assert run(False) == run(True)
