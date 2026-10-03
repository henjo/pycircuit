"""The evaluate core (`pycircuit/circuit/_tran_core.py`, speed round 7's
first commit, 2026-10-03): the transient's three evaluation sites served by
one C call must give the Python path's bytes and leave the Python path's
state, decline to it wherever they cannot serve, and step a transient and a
PSS bit for bit."""
import contextlib
import os
import subprocess
import sys
import types
import warnings

import numpy as np
import pytest

from pycircuit.circuit import PSS, _tran_core, simwarnings
from pycircuit.circuit import _hdl_cbackend as cb
from pycircuit.circuit import circuit as cm
from pycircuit.circuit.elements import VS, Diode, R, SubCircuit, VSin, gnd
from pycircuit.circuit.integrator import (
    EulerIntegrator,
    Gear2Integrator,
    ThetaIntegrator,
    TrapezoidalIntegrator,
)
from pycircuit.circuit.tests.test_hdl_batch import mos_chain
from pycircuit.circuit.transient import Transient

METHODS = ('gear', 'euler', 'trap', 'theta')
INTEGRATOR = {'gear': Gear2Integrator, 'euler': EulerIntegrator,
              'trap': TrapezoidalIntegrator, 'theta': ThetaIntegrator}


@contextlib.contextmanager
def core(on):
    was = _tran_core.CORE
    _tran_core.CORE = on
    try:
        yield
    finally:
        _tran_core.CORE = was


def _stage():
    """The harness's PSP stage: one C-bound device, a lone element."""
    from pycircuit.circuit import compact
    cm.default_toolkit = cm.numeric
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.2)
    c['vg'] = VSin('g', gnd, v=0.7, va=2e-2, freq=1e6)
    c['rl'] = R('vdd', 'd', r=5e3)
    c['M'] = compact.PspMosLongChannel('d', 'g', gnd, gnd, fnt=1.0)
    return c


def _tr(cir, method='gear', steps=6):
    """A transient a few steps in, with the step state a residual needs,
    and its last accepted state."""
    tr = Transient(cir, toolkit=cm.numeric, integrator=INTEGRATOR[method]())
    res = tr.solve(tend=steps * 2e-8, timestep=2e-8, fixed_timestep=True)
    if _tran_core.STATUS != 'c':
        pytest.skip(f'the core is off: {_tran_core.STATUS}')
    if _tran_core.core_for(tr) is None:
        from pycircuit.circuit.tests.test_hdl_batch import _why
        facts = _why(cir, 'G')
        if any(v['bound'] for v in facts.values()):
            pytest.fail(f'a bound class the core does not serve: {facts}')
        pytest.skip('the chain is not served (a class not C-bound)')
    x = np.ascontiguousarray(np.asarray(res.x, float)[:, -1])
    tr._u_memo = {}
    return tr, x


def _state(tr):
    out = [tr._iq.tobytes(), np.asarray(tr._Geq).tobytes(), np.asarray(tr._Cmat).tobytes(),
           repr(tr._companion_coeffs), tr._effective_method, type(tr.active_integrator).__name__]
    for name in ('_q_cache', '_C_cache'):
        c = getattr(tr, name, None)
        out.append(None if c is None else (c[0].tobytes(), np.asarray(c[1]).tobytes()))
    return out


def _python_j(tr, x):
    """`jacobian_only`'s body (transient.py), the Python path."""
    _iq, Geq = tr._companion_at(x)
    return tr.cir.G(x, tr.epar) + Geq


def _python_f(tr, x, t):
    """`residual_only`'s body, the Python path."""
    q = tr.cir.q(x, tr.epar)
    iq, _geq = tr.get_diff(q, tr._Cmat)
    u = tr._source_at(t, None)
    return tr.cir.i(x, tr.epar) + iq + u


def _states(x, k, seed=0):
    rng = np.random.default_rng(seed)
    out = [x]
    for _ in range(k):
        out.append(x + rng.normal(0.0, 0.3, x.size))
    out.append(np.full(x.size, 0.7))
    out.append(np.zeros(x.size))
    return out


@pytest.mark.parametrize('method', METHODS)
def test_the_residual_is_the_python_paths_bytes_and_state(method):
    tr, x = _tr(mos_chain(8), method)
    t = 1.1e-7
    for j, xt in enumerate(_states(x, 12)):
        tr._is_first_step = (j % 5 == 4)        # the order drop, where the method has one
        with core(True):
            a = tr._residual_and_jacobian(xt, t)
            sa = _state(tr)
        with core(False):
            b = tr._residual_and_jacobian(xt, t)
            sb = _state(tr)
        tr._is_first_step = False
        assert a[0].dtype == b[0].dtype and a[1].shape == b[1].shape
        assert a[0].tobytes() == b[0].tobytes(), (method, j)
        assert a[1].tobytes() == b[1].tobytes(), (method, j)
        assert sa == sb, (method, j)


@pytest.mark.parametrize('make', [lambda: mos_chain(8), _stage], ids=['chain', 'lone'])
def test_the_converged_point_and_the_chord_residual(make):
    tr, x = _tr(make())
    t = 1.1e-7
    for xt in _states(x, 6, seed=2):
        with core(True):
            none, J = _tran_core.evaluate(tr, xt, t, None, 'j')
            sa = _state(tr)
        assert none is None
        with core(False):
            Jp = _python_j(tr, xt)
            sb = _state(tr)
        assert J.tobytes() == Jp.tobytes() and sa == sb
        ## the chord's residual: the conductance from the held `_Cmat`
        with core(True):
            f = _tran_core.evaluate(tr, xt, t, None, 'f')
            sa = _state(tr)
        with core(False):
            fp = _python_f(tr, xt, t)
            sb = _state(tr)
        assert f.tobytes() == fp.tobytes() and sa == sb


def test_a_cached_capacitance_is_used_as_python_uses_it():
    """`_C_at_state` serves `C` from the memo or the converged-point cache
    before evaluating; the core asks the same lookup and skips the pass."""
    tr, x = _tr(mos_chain(6))
    fake = np.full((tr.cir.n, tr.cir.n), 1e-15)
    for on in (True, False):
        tr._C_cache = (x, fake)
        with core(on):
            _f, J = tr._residual_and_jacobian(x, 1.1e-7)
        assert tr._Cmat is fake
        if on:
            Jc = J
    assert Jc.tobytes() == J.tobytes()


def test_what_the_core_declines_takes_the_python_path():
    tr, x = _tr(mos_chain(6))
    t = 1.1e-7
    def ok():
        return _tran_core.evaluate(tr, x, t, None, 'fj') is not None
    assert ok()
    with core(False):
        assert not ok()
    assert _tran_core.evaluate(tr, x.astype(np.float32), t, None, 'fj') is None
    nan = x.copy()
    nan[2] = np.nan
    assert _tran_core.evaluate(tr, nan, t, None, 'fj') is None
    assert _tran_core.evaluate(tr, x[:-1], t, None, 'fj') is None
    ## a temperature that is not one number
    tr.epar = types.SimpleNamespace(T=np.array([300.0, 310.0]), bypasstol=-1.0)
    assert not ok()
    tr.epar = types.SimpleNamespace(T=np.array(310.0), bypasstol=-1.0)
    assert ok()                                 # (a 0-d array is one number)
    ## an integrator outside the enum
    base = tr.base_integrator
    tr.base_integrator = types.SimpleNamespace(
        check_order_drop=lambda h, hl, first: types.SimpleNamespace())
    assert not ok()
    tr.base_integrator = base
    ## an instance shadow of one of the four methods, then its removal
    el = tr.cir['M2']
    el.__dict__['q'] = lambda xx, epar=None, params_tree=None: np.zeros(el.n)
    assert not ok()
    del el.__dict__['q']
    assert ok()
    ## the switch off, and the bytes of the path it falls to
    with core(True):
        a = tr._residual_and_jacobian(x, t)
    with core(False):
        b = tr._residual_and_jacobian(x, t)
    assert a[1].tobytes() == b[1].tobytes()


def test_an_element_the_core_cannot_evaluate_leaves_the_circuit_to_python():
    cir = mos_chain(6)
    cir['D'] = Diode('d3', gnd)
    tr = Transient(cir, toolkit=cm.numeric)
    tr.solve(tend=4e-8, timestep=2e-8, fixed_timestep=True)
    assert _tran_core.core_for(tr) is None
    ## and the record is kept until the plan rebuilds
    assert tr.__dict__['_tran_core'][1] is None
    got = {}
    for on in (True, False):
        c2 = mos_chain(6)
        c2['D'] = Diode('d3', gnd)
        with core(on):
            res = Transient(c2, toolkit=cm.numeric).solve(
                tend=20 * 2e-8, timestep=2e-8, fixed_timestep=True)
        got[on] = np.asarray(res.x, float)
    assert got[True].tobytes() == got[False].tobytes()


def test_a_constant_only_circuit_is_served_with_no_kernel_at_all():
    """A circuit of hand-written constant elements has no kernel to call
    and still goes through the core: the templates, the constant
    products and the companion in C -- the parity guard's RC ladder runs
    twice as fast (401 -> 214 ms for its pairs), bit for bit."""
    from pycircuit.circuit.elements import C
    got = {}
    for on in (True, False):
        cm.default_toolkit = cm.numeric
        cir = SubCircuit()
        cir['vs'] = VSin('n0', gnd, va=2.0, freq=1e5)
        for k in range(8):
            cir[f'R{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
            cir[f'C{k}'] = C(f'n{k + 1}', gnd, c=1e-9)
        with core(on):
            tr = Transient(cir, toolkit=cm.numeric)
            res = tr.solve(tend=20e-6, timestep=1e-6)
        got[on] = np.asarray(res.x, float).tobytes()
        if on and _tran_core.STATUS == 'c':
            c_ = _tran_core.core_for(tr)
            assert c_ is not None and c_.nb == 0 and c_.ng > 0
    assert got[True] == got[False]


@pytest.mark.parametrize('method', METHODS)
def test_a_transient_steps_bit_for_bit(method):
    got = {}
    for on in (True, False):
        cir = mos_chain(6)
        with core(on):
            tr = Transient(cir, toolkit=cm.numeric, integrator=INTEGRATOR[method]())
            res = tr.solve(tend=25 * 2e-8, timestep=2e-8, fixed_timestep=True)
        got[on] = (np.asarray(res.x, float).tobytes(), tr.statistics.newton_iterations,
                   _tran_core.core_for(tr) is not None, cir)
    assert got[True][0] == got[False][0] and got[True][1] == got[False][1]
    if _tran_core.STATUS == 'c':
        from pycircuit.circuit.tests.test_hdl_batch import _why
        assert got[True][2], _why(got[True][3], 'G')


def test_the_stage_pss_is_bit_for_bit():
    got = {}
    for on in (True, False):
        ## (the harness's PSS: its grid is coarse for the accuracy it asks,
        ## and says so; the comparison is the bytes, so the warning is quiet)
        with core(on), warnings.catch_warnings():
            warnings.simplefilter('ignore', simwarnings.AccuracyWarning)
            p = PSS(_stage(), method='gear', reltol=1e-8)
            p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
        got[on] = np.asarray(p.waveform[1], float).tobytes()
    assert got[True] == got[False]


def test_a_core_with_the_wrong_accumulation_order_is_caught(monkeypatch):
    tr, x = _tr(mos_chain(6))
    src = _tran_core.CORE_C.replace('for (s = 0; s < L; s++) bin[flat[s]] += buf[s];',
                                    'for (s = L - 1; s >= 0; s--) bin[flat[s]] += buf[s];')
    assert src != _tran_core.CORE_C
    ffi, cfn, key, _cold, _secs = cb.load_kernel(src, _tran_core.CORE_CDEF)
    assert key != cb.source_key(_tran_core.CORE_C)
    drv = _tran_core.driver()
    monkeypatch.setattr(_tran_core, '_driver', (ffi, cfn, drv[2]))
    tr.__dict__.pop('_tran_core', None)
    with core(True):
        a = tr._residual_and_jacobian(x, 1.1e-7)
    with core(False):
        b = tr._residual_and_jacobian(x, 1.1e-7)
    assert a[1].tobytes() != b[1].tobytes()


def test_the_switch_is_read_from_the_environment():
    env = dict(os.environ, PYCIRCUIT_TRAN_CORE='0')
    code = 'from pycircuit.circuit import _tran_core; print(_tran_core.CORE)'
    out = subprocess.run([sys.executable, '-c', code], env=env,
                         capture_output=True, text=True, check=True)
    assert out.stdout.strip() == 'False'
