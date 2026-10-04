"""PSP's hand-written limiter in C (`_hdl_climit.PSP_LIMIT_C`, speed round 9,
stage 1): the twin is `compact.PspMosLongChannel.limit` bit for bit or it
declines; with it the walk serves PSP elements and the C Newton serves PSP
circuits, with the Python path's answers."""
import itertools
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _hdl_climit, _paths, circuit, compact
from pycircuit.circuit.elements import VS, R, SubCircuit, VSin, gnd
from pycircuit.circuit.transient import Transient

PSP = compact.PspMosLongChannel


def _stage(pmos=False):
    c = SubCircuit()
    for n in ('g', 'd', 'vdd'):
        c.add_node(n)
    c['vdd'] = VS('vdd', gnd, v=1.2)
    if pmos:
        c['vg'] = VSin('g', gnd, v=0.5, va=2e-2, freq=1e6)
        c['rl'] = R('d', gnd, r=5e3)
        c['M'] = compact.PspPmosLongChannel('d', 'g', 'vdd', 'vdd', fnt=1.0)
    else:
        c['vg'] = VSin('g', gnd, v=0.7, va=2e-2, freq=1e6)
        c['rl'] = R('vdd', 'd', r=5e3)
        c['M'] = PSP('d', 'g', gnd, gnd, fnt=1.0)
    return c


def _chain(n=3, va=2e-2):
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.2)
    c.add_node('g0')
    c['vg'] = VSin('g0', gnd, v=0.7, va=va, freq=1e6)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'M{k}'] = PSP(f'd{k}', f'g{k}' if k == 0 else f'd{k - 1}', gnd, gnd, fnt=1.0)
    return c


def _twin(el):
    cls = type(el)
    hw = _hdl_climit._handwritten(el, cls, cls._hdl_info)
    if hw is None:
        pytest.skip('PSP is not C-bound here (no compiler or store)')
    return hw.kern


def _compare(el, kern, epar, x, x0):
    """The twin's answer, None where it declines; its bytes are the law's."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ref = PSP.limit(el, x.copy(), x0.copy(), epar)
    got = kern(el, x.copy(), x0.copy(), epar)
    if got is not None:
        assert got.tobytes() == np.asarray(ref, float).tobytes(), (x, x0, ref, got)
    return got


_SPECIAL = [0.0, -0.0, 1.0, -1.0, 1.0 + 2 ** -52, 1.0 - 2 ** -53, 2.0, -0.5, 1e300, -1e300,
            5e-324, np.inf, -np.inf, np.nan]


def test_enumerated_inputs_are_the_law_or_declined():
    """Every row of x and x0 from a set holding the limit's edges, zeros of
    both signs, huge values and non-finite ones: the law's bytes, or a
    decline exactly where a value is not finite."""
    el = _stage()['M']
    kern = _twin(el)
    epar = Transient(_stage(), toolkit=circuit.numeric).epar
    n = kern.nx
    rng = np.random.default_rng(5)
    served = declined = 0
    for _ in range(6000):
        x = rng.choice(_SPECIAL, size=n).astype(float)
        x0 = rng.choice(_SPECIAL, size=n).astype(float)
        if _compare(el, kern, epar, x, x0) is None:
            declined += 1
        else:
            served += 1
    assert served > 500 and declined > 500
    ## the limit itself on each side, with finite values
    for vs, vold in itertools.product((0.0, -0.0, 0.3), (0.0, -0.0, 1.0)):
        for delta in (1.0, -1.0, np.nextafter(1.0, 2.0), np.nextafter(-1.0, -2.0), 0.0, -0.0):
            x0 = np.zeros(n)
            x0[2], x0[0], x0[1], x0[3] = vs, vs + vold, vs + vold, vs + vold
            x = x0.copy()
            x[0] = x[1] = x[3] = vs + vold + delta
            assert _compare(el, kern, epar, x, x0) is not None


@settings(deadline=None)
@given(data=st.data(), scale=st.sampled_from([0.1, 0.9, 1.0, 1.1, 5.0, 1e3]))
def test_drawn_steps_are_the_law(data, scale):
    el = _stage()['M']
    kern = _twin(el)
    epar = Transient(_stage(), toolkit=circuit.numeric).epar
    n = kern.nx
    fl = st.floats(-3.0, 3.0, allow_nan=False)
    x0 = np.array(data.draw(st.lists(fl, min_size=n, max_size=n)))
    step = np.array(data.draw(st.lists(st.floats(-scale, scale, allow_nan=False),
                                       min_size=n, max_size=n)))
    assert _compare(el, kern, epar, x0 + step, x0) is not None


def _run(build, on, monkeypatch, **kw):
    monkeypatch.setattr(_hdl_climit, 'PSP_LIMIT', on)
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        ## (a finalizer's warning is not the run's: a file another test left in
        ## a reference cycle closes whenever the collector runs -- inside either
        ## run, 2026-10-04)
        warnings.simplefilter('ignore', ResourceWarning)
        tr = Transient(build(), toolkit=circuit.numeric, **kw.pop('make', {}))
        res = tr.solve(**kw)
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), _paths.since(before))


def _same(build, monkeypatch, served_min=1, **kw):
    a = _run(build, False, monkeypatch, **dict(kw))
    b = _run(build, True, monkeypatch, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert b[4].get('newton_c:served', 0) >= served_min, b[4]
    return a, b


@pytest.mark.parametrize('pmos', [False, True], ids=['nmos', 'pmos'])
@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
def test_the_c_newton_serves_psp_with_the_python_answers(pmos, fixed, monkeypatch):
    kw = {'tend': 1e-6, 'timestep': 2e-8}
    if fixed:
        kw['fixed_timestep'] = True
    a, _b = _same(lambda: _stage(pmos), monkeypatch, **kw)
    assert a[4].get('newton_c:served', 0) == 0 and a[4].get('newton_c:limit', 0) > 0


def test_a_chain_driven_hard_enough_to_limit(monkeypatch):
    """A chain whose first steps move the gates by volts: the law limits,
    and the walk's twin limits the same rows by the same amounts."""
    _same(lambda: _chain(3, va=0.6), monkeypatch, tend=1e-6, timestep=1e-7, fixed_timestep=True)


def test_an_instance_vlimit_or_a_patched_limit_takes_the_law(monkeypatch):
    """`vlimit` set on one element, `limit` patched on the class: the Python
    law answers (called, and the same answers as with the twin off)."""
    def shadowed():
        c = _chain(3, va=0.6)
        c['M1'].vlimit = 0.4
        return c
    _a, b = _same(shadowed, monkeypatch, served_min=0, tend=4e-7, timestep=1e-7, fixed_timestep=True)
    assert b[4].get('newton_c:served', 0) == 0
    calls = []
    real = PSP.limit

    def counting(self, *a, **k):
        calls.append(1)
        return real(self, *a, **k)
    monkeypatch.setattr(PSP, 'limit', counting)
    _, b = _same(lambda: _chain(2), monkeypatch, served_min=0, tend=4e-7, timestep=1e-7,
                 fixed_timestep=True)
    assert calls and b[4].get('newton_c:served', 0) == 0


def test_a_class_vlimit_gets_its_own_twin(monkeypatch):
    """The class's `vlimit` changed: the old twin steps aside, a new walk
    builds the twin for the new value, and it is the law at that value."""
    monkeypatch.setattr(PSP, 'vlimit', 0.25)
    _same(lambda: _chain(3, va=0.6), monkeypatch, tend=6e-7, timestep=1e-7, fixed_timestep=True)
    el = _chain(1)['M0']
    assert _twin(el) is not None
