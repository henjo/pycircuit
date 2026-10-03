"""One C call per class per pass (`pycircuit/circuit/_hdl_batch.py`, speed
round 6, 2026-10-03): the stamp plan's batches of C-bound hdl elements must
give the per-element path's bytes, call no element method for an element a
batch serves, hand back what they cannot serve -- an instance shadow, a
class patch, a detached class, a temperature that is not one number -- and
mirror the element's own pack, never repack it."""
import contextlib
import os
import subprocess
import sys
import types
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _hdl_batch, _limiting, hdl
from pycircuit.circuit import _hdl_cbackend as cb
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VSin, gnd
from pycircuit.circuit.tests.test_hdl_cbackend import CHAINED, KW, _points

PASSES = ('G', 'C', 'i', 'q')


@contextlib.contextmanager
def batching(on):
    was = _hdl_batch.ENABLED
    _hdl_batch.ENABLED = on
    try:
        yield
    finally:
        _hdl_batch.ENABLED = was


@pytest.fixture
def kernel_calls(monkeypatch):
    """Every `CKernel.__call__` -- the per-element C path -- recorded (the
    limiter's kernel has a `__call__` of its own and is not)."""
    calls = []
    orig = cb.CKernel.__call__

    def spy(self, element, x, epar):
        calls.append(element)
        return orig(self, element, x, epar)
    monkeypatch.setattr(cb.CKernel, '__call__', spy)
    return calls


def mos_chain(ndev=20, make=None):
    """The harness's chain: `ndev` common-source stages, each gate on the
    previous drain, a load resistor each (benchmarks/step_machinery.py)."""
    cm.default_toolkit = cm.numeric
    make = make or (lambda d, g: eh.MosLevel1Hdl(d, g, gnd, gnd))
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['vg'] = VSin('g0', gnd, v=0.9, va=2e-2, freq=1e6)
    for k in range(ndev):
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'M{k}'] = make(f'd{k}', 'g0' if k == 0 else f'd{k - 1}')
    return c


def states(n, k=12, seed=0):
    rng = np.random.default_rng(seed)
    out = [rng.uniform(-2.0, 2.0, n) for _ in range(k)]
    out.append(np.full(n, 0.7))
    out.append(np.zeros(n))
    x = rng.uniform(-2.0, 2.0, n)
    x[1] = np.nan
    x[2] = np.inf
    out.append(x)
    return out


def _bound(cir, inst='M0'):
    """The chain's device class is C-bound here -- or the test has nothing
    to say (no compiler and no stored objects, the compile cache off)."""
    cls = type(cir[inst])
    if not cls._hdl_info.get('_c_bound'):
        pytest.skip(f'{cls.__name__} is not C-bound: {cls._hdl_backend_status}')
    return cls


def _both(cir, m, x, *args):
    with batching(True), np.errstate(all='ignore'):
        a = getattr(cir, m)(x, *args)
    with batching(False), np.errstate(all='ignore'):
        b = getattr(cir, m)(x, *args)
    assert a.dtype == b.dtype and a.shape == b.shape, m
    return a, b


def _batches(cir, m):
    return [(bt.cls, bt.n) for bt in cir.__dict__['_stamp_plan'].methods[m].batches]


def test_a_pass_is_one_call_per_class_and_the_per_element_bytes(kernel_calls):
    cir = mos_chain()
    cls = _bound(cir)
    for x in states(cir.n):
        for m in PASSES:
            with batching(False), np.errstate(all='ignore'):
                ref = getattr(cir, m)(x)
            assert len(kernel_calls) == 20, m
            kernel_calls.clear()
            with batching(True), np.errstate(all='ignore'):
                got = getattr(cir, m)(x)
            ## (a vector pass at a non-finite state is the legacy loop's on
            ## both sides, as before: the per-element path, twenty calls)
            loop = m in ('i', 'q') and not np.isfinite(x).all()
            assert len(kernel_calls) == (20 if loop else 0), m
            kernel_calls.clear()
            assert got.dtype == ref.dtype and got.shape == ref.shape
            assert got.tobytes() == ref.tobytes(), m
    plan = cir.__dict__['_stamp_plan']
    assert plan.builds == 1
    for m in PASSES:
        assert _batches(cir, m) == [(cls, 20)]
        assert plan.methods[m].calls == []       # (the rest are constant)


@pytest.mark.parametrize('name', CHAINED)
def test_every_chained_library_class_batches_bit_for_bit(name, kernel_calls):
    """Five instances of each chained library class on their own nodes: a
    C-bound class without DC pins is one call and the per-element bytes; a
    class with DC pins (its `G_dc`/`i_dc` choice is the per-element path's)
    or a class running numpy takes the loop, the bytes still the loop's."""
    cls = getattr(eh, name)
    cm.default_toolkit = cm.numeric
    cir = SubCircuit()
    nt = len(cls.terminals)
    for j in range(5):
        cir[f'e{j}'] = cls(*[f'n{j}_{t}' for t in range(nt)], **KW.get(name, {}))
    for j in range(5):
        cir[f'r{j}'] = R(f'n{j}_0', gnd, r=1e3)
    v = type(cir['e0'])
    info = v._hdl_info
    expect = (info.get('_c_bound') and not info['state_meta']['dc_pins'])
    for x in _points(cir.n, 8):
        for m in PASSES:
            a, b = _both(cir, m, x)
            assert a.tobytes() == b.tobytes(), (name, m)
    for m in PASSES:
        assert _batches(cir, m) == ([(v, 5)] if expect else []), (name, m)
    if expect:
        kernel_calls.clear()
        with batching(True), np.errstate(all='ignore'):
            for m in PASSES:
                getattr(cir, m)(_points(cir.n, 1)[0])
        assert kernel_calls == []


def test_a_mixed_chain_batches_each_class_and_leaves_the_rest(kernel_calls):
    """Two device classes interleaved with hand-written R and C and a lone
    hdl diode: a batch per class, the diode per element (one instance)."""
    def make(d, g):
        k = int(d[1:])
        if k % 2:
            return eh.MosLevel1Hdl(d, g, gnd, gnd)
        return eh.GummelPoonNpnHdl(d, g, gnd, rc=2.0, re=1.0, rb=100.0)
    cir = mos_chain(8, make)
    for k in range(8):
        cir[f'c{k}'] = C(f'd{k}', gnd, c=1e-13)
    cir['dh'] = eh.DiodeHdl('d7', gnd)
    _bound(cir, 'M1')
    _bound(cir, 'M0')
    for x in states(cir.n, 6):
        for m in PASSES:
            a, b = _both(cir, m, x)
            assert a.tobytes() == b.tobytes(), m
    kernel_calls.clear()
    with batching(True), np.errstate(all='ignore'):
        cir.G(states(cir.n, 1)[0])
    ## the diode's one kernel call is the per-element path's
    assert kernel_calls == [cir['dh']] if cir['dh'].__class__._hdl_info.get(
        '_c_bound') else kernel_calls == []
    sizes = sorted(n for _c, n in _batches(cir, 'G'))
    assert sizes == [4, 4]
    assert [inst for inst, *_ in cir.__dict__['_stamp_plan'].methods['G'].calls] == ['dh']


def test_an_instance_shadow_mid_run_is_called_and_the_rest_stays_one_call(kernel_calls):
    """PCNR installs `el.i`/`el.G` on the participating devices DURING a
    solve, with no plan rebuild: that element is called, the rest batched."""
    cir = mos_chain(6)
    _bound(cir)
    x = states(cir.n, 1)[0]
    with batching(True):
        g0 = cir.G(x)
    el = cir['M2']
    seen = []

    def shadow(xx, epar=defaultepar, params_tree=None):
        seen.append(xx.copy())
        return type(el).G(el, xx, epar) * 2.0
    el.__dict__['G'] = shadow
    kernel_calls.clear()
    with batching(True):
        a = cir.G(x)
    assert len(seen) == 1 and kernel_calls == [el]
    with batching(False):
        b = cir.G(x)
    assert len(seen) == 2 and len(kernel_calls) == 7
    assert a.tobytes() == b.tobytes()
    assert a.tobytes() != g0.tobytes()
    assert cir.__dict__['_stamp_plan'].builds == 1
    ## other passes untouched by the shadow
    kernel_calls.clear()
    with batching(True):
        cir.i(x)
    assert kernel_calls == []
    del el.__dict__['G']
    kernel_calls.clear()
    with batching(True):
        c = cir.G(x)
    assert kernel_calls == [] and c.tobytes() == g0.tobytes()


def test_a_class_patch_mid_run_takes_the_loop_for_that_method(kernel_calls, monkeypatch):
    cir = mos_chain(6)
    cls = _bound(cir)
    x = states(cir.n, 1)[0]
    with batching(True):
        cir.G(x)
        cir.i(x)
    gen = cls.G
    seen = []

    def spy(self, xx, epar=defaultepar, params_tree=None):
        seen.append(self)
        return gen(self, xx, epar)
    monkeypatch.setattr(cls, 'G', spy)
    assert not _hdl_batch.is_generated(cls, 'G') and _hdl_batch.is_generated(cls, 'i')
    kernel_calls.clear()
    with batching(True):
        a = cir.G(x)
    assert len(seen) == 6 and len(kernel_calls) == 6
    with batching(False):
        b = cir.G(x)
    assert a.tobytes() == b.tobytes()
    kernel_calls.clear()
    with batching(True):
        cir.i(x)
        cir.q(x)
        cir.C(x)
    assert kernel_calls == []
    assert cir.__dict__['_stamp_plan'].builds == 1


def test_a_detached_class_takes_the_loop_and_binds_back(kernel_calls):
    cir = mos_chain(6)
    cls = _bound(cir)
    x = states(cir.n, 1)[0]
    with batching(True):
        g0 = cir.G(x)
    hdl.set_backend('numpy', eh.MosLevel1Hdl)
    try:
        assert not cls._hdl_info.get('_c_bound')
        kernel_calls.clear()
        a, b = _both(cir, 'G', x)
        assert kernel_calls == [] and a.tobytes() == b.tobytes()
    finally:
        hdl.set_backend(None, eh.MosLevel1Hdl)
    assert cls._hdl_info.get('_c_bound')
    kernel_calls.clear()
    with batching(True):
        c = cir.G(x)
    assert kernel_calls == [] and c.tobytes() == g0.tobytes()


def test_a_parameter_change_repacks_that_element_only(kernel_calls):
    cir = mos_chain(6)
    _bound(cir)
    x = states(cir.n, 1)[0]
    with batching(True):
        g0 = cir.G(x)
    packs = {k: cir[k].__dict__['_hdl_cp'] for k in cir.elements if k.startswith('M')}
    ## (the element's own update: the circuit's notifies every element,
    ## and every pack is rebuilt on both paths alike)
    cir['M3'].ipar.vto = 0.35
    cir['M3'].update_iparv(cir.iparv)
    assert '_hdl_cp' not in cir['M3'].__dict__
    a, b = _both(cir, 'G', x)
    assert a.tobytes() == b.tobytes() and a.tobytes() != g0.tobytes()
    for k, cp in packs.items():
        if k == 'M3':
            assert cir[k].__dict__['_hdl_cp'] is not cp
        else:
            assert cir[k].__dict__['_hdl_cp'] is cp
    assert cir.__dict__['_stamp_plan'].builds == 2


def test_a_restored_state_is_read_through_its_old_pack(kernel_calls):
    """`state_restore` (the branch check's speculative Newton) puts an
    element's old `__dict__` back -- its old pack included -- with no epoch
    move: both paths read THAT pack, so the batch mirrors it and never
    repacks what the element holds."""
    cir = mos_chain(6)
    _bound(cir)
    x = states(cir.n, 1)[0]
    with batching(True):
        g0 = cir.G(x)
    old = cir['M3'].__dict__['_hdl_cp']
    snap = _limiting.state_snapshot(cir)
    cir['M3'].ipar.vto = 0.35
    cir.update_iparv()
    with batching(True):
        g1 = cir.G(x)
    assert g1.tobytes() != g0.tobytes()
    _limiting.state_restore(snap)
    assert cir['M3'].__dict__['_hdl_cp'] is old
    builds = cir.__dict__['_stamp_plan'].builds
    a, b = _both(cir, 'G', x)
    assert a.tobytes() == b.tobytes() == g0.tobytes()
    assert cir['M3'].__dict__['_hdl_cp'] is old
    assert cir.__dict__['_stamp_plan'].builds == builds


def test_the_temperature_is_written_as_the_kernel_writes_it(kernel_calls):
    cir = mos_chain(6)
    _bound(cir)
    x = states(cir.n, 1)[0]
    at = {}
    for T in (300.0, 350.0, 320, np.array(330.0), np.float64(340.0)):
        epar = types.SimpleNamespace(T=T)
        kernel_calls.clear()
        a, b = _both(cir, 'G', x, epar)
        assert a.tobytes() == b.tobytes(), T
        assert len(kernel_calls) == 6, T        # (the loop's side only)
        at[float(T)] = a.tobytes()
    assert len(set(at.values())) == len(at)     # every T told apart
    ## a temperature that is not one number: the kernel steps aside on
    ## both paths (numpy broadcasts); the batch hands the class back
    epar = types.SimpleNamespace(T=np.array([300.0, 310.0]))
    with pytest.raises(Exception), batching(False):  # noqa: B017 -- the loop's own failure
        cir.G(x, epar)


def test_a_fallback_calls_the_batched_elements_once(kernel_calls):
    """A non-real stamp from a non-batched element sends the pass to the
    legacy loop: the batched elements are then called there, once each,
    never through a batch AND the loop."""
    class _Odd(R):
        def G(self, x, epar=defaultepar):
            return R.G(self, x, epar) * (1.0 + 0.0j)
    cir = mos_chain(6)
    cir['odd'] = _Odd('d5', gnd, r=1e4)
    _bound(cir)
    x = states(cir.n, 1)[0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', np.exceptions.ComplexWarning)
        kernel_calls.clear()
        with batching(True):
            a = cir.G(x)
        assert len(kernel_calls) == 6
        kernel_calls.clear()
        with batching(False):
            b = cir.G(x)
        assert len(kernel_calls) == 6
    assert a.tobytes() == b.tobytes()
    assert _batches(cir, 'G') == [(type(cir['M0']), 6)]


def test_a_transient_steps_bit_for_bit_with_one_plan():
    from pycircuit.circuit.transient import Transient
    got = {}
    for on in (True, False):
        cir = mos_chain(6)
        _bound(cir)
        with batching(on):
            res = Transient(cir, toolkit=cm.numeric).solve(
                tend=20 * 2e-8, timestep=2e-8, fixed_timestep=True)
        got[on] = np.asarray(res.x, dtype=float)
        assert cir.__dict__['_stamp_plan'].builds == 1
    assert got[True].tobytes() == got[False].tobytes()


def test_a_driver_with_the_wrong_stride_is_caught(kernel_calls, monkeypatch):
    """The check is the bytes: a driver stepping the state by one double
    too few is not the per-element path."""
    cir = mos_chain(6)
    _bound(cir)
    x = states(cir.n, 1)[0]
    with batching(True):
        cir.G(x)
    src = _hdl_batch.PASS_C.replace('X + e * sx', 'X + e * (sx - 1)')
    assert src != _hdl_batch.PASS_C
    ffi, cfn, key, _cold, _secs = cb.load_kernel(src, _hdl_batch.PASS_CDEF)
    assert key != cb.source_key(_hdl_batch.PASS_C)
    monkeypatch.setattr(_hdl_batch, '_driver',
                        (ffi, cfn, ffi.typeof('double *'), ffi.typeof('double **')))
    a, b = _both(cir, 'G', x)
    assert a.tobytes() != b.tobytes()


def test_the_switch_is_read_from_the_environment():
    env = dict(os.environ, PYCIRCUIT_HDL_BATCH='0')
    code = 'from pycircuit.circuit import _hdl_batch; print(_hdl_batch.ENABLED)'
    out = subprocess.run([sys.executable, '-c', code],
                         env=env, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == 'False'
    assert _hdl_batch.STATUS in ('c', 'not loaded') or _hdl_batch.STATUS.startswith('off')


def test_the_generated_methods_are_told_from_their_doubles():
    gen = _hdl_batch.generated_code()
    assert {'i', 'G', 'q', 'C', 'limit', 'u', 'dudt'} <= set(gen)
    for m in PASSES:
        assert _hdl_batch.is_generated(eh.MosLevel1Hdl, m)

    class _Sub(eh.MosLevel1Hdl):
        pass
    assert _hdl_batch.is_generated(_Sub, 'G')

    class _Over(eh.MosLevel1Hdl):
        def G(self, x, epar=defaultepar, params_tree=None):
            return eh.MosLevel1Hdl.G(self, x, epar)
    assert not _hdl_batch.is_generated(_Over, 'G') and _hdl_batch.is_generated(_Over, 'i')
    assert not _hdl_batch.is_generated(R, 'G')
