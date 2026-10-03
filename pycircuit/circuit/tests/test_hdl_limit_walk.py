"""The limiter walk (`_hdl_climit.limit_walk`, speed round 6's third commit,
2026-10-03): `SubCircuit.limit`'s loop with every run of C-kernel elements
walked in one C call on the live state -- the loop's bytes in every case
the loop handles (shared nodes, duplicate rows, `x0 is x`, a hand-written
limiter mid-chain, a shadow, a decline, a rebound kernel, a nested
circuit), and the loop itself where the walk does not serve."""
import contextlib
import os
import subprocess
import sys
import types

import numpy as np
import pytest

from pycircuit.circuit import _hdl_climit as cl
from pycircuit.circuit import _limiting, hdl
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import VS, Diode, R, SubCircuit, VSin, gnd
from pycircuit.circuit.tests.test_hdl_batch import mos_chain


@contextlib.contextmanager
def walking(on):
    was = cl.WALK
    cl.WALK = on
    try:
        yield
    finally:
        cl.WALK = was


@pytest.fixture
def limit_kernel_calls(monkeypatch):
    """Every `CLimitKernel.__call__` -- the per-element kernel path."""
    calls = []
    orig = cl.CLimitKernel.__call__

    def spy(self, element, x, x0, epar):
        calls.append(element)
        return orig(self, element, x, x0, epar)
    monkeypatch.setattr(cl.CLimitKernel, '__call__', spy)
    return calls


def _bound(cir, inst='M0'):
    cls = type(cir[inst])
    if not isinstance(cls._hdl_info.get('_c_limit'), cl.CLimitKernel):
        pytest.skip(f'{cls.__name__} has no limit kernel: '
                    f'{cls._hdl_info.get("_c_limit_status")}')
    return cls


def _pairs(n, k=10, seed=1):
    """`(x, x0)`: an iterate a Newton step away from the accepted one."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(k):
        x0 = rng.uniform(-1.5, 2.0, n)
        out.append((x0 + rng.normal(0.0, 0.8, n), x0))
    return out


def _both(cir, x, x0, epar=defaultepar, same=False):
    """`cir.limit` through the walk and through the loop, on copies; the
    walk's result, the two asserted byte-equal."""
    got = []
    for on in (True, False):
        xa = x.copy()
        x0a = xa if same else x0.copy()
        with walking(on):
            r = cir.limit(xa, x0a, epar)
        assert r is xa
        got.append(xa)
    assert got[0].tobytes() == got[1].tobytes()
    return got[0]


def test_the_walk_is_one_call_and_the_loops_bytes(limit_kernel_calls):
    cir = mos_chain()
    _bound(cir)
    moved = 0
    for x, x0 in _pairs(cir.n):
        xa = x.copy()
        with walking(True):
            cir.limit(xa, x0.copy(), defaultepar)
        assert limit_kernel_calls == []
        xb = x.copy()
        with walking(False):
            cir.limit(xb, x0.copy(), defaultepar)
        assert len(limit_kernel_calls) == 20
        limit_kernel_calls.clear()
        assert xa.tobytes() == xb.tobytes()
        moved += xa.tobytes() != x.tobytes()
    assert moved > 0, 'no pair was limited: the test says nothing'
    w = cir.__dict__['_limit_walk']
    assert len(w.entries) == 20 and w.any_capable
    ## (source and bulk on gnd: duplicate rows, the last value written)
    assert all(len(set(nm.tolist())) < len(nm) for _i, _e, nm in w.entries)


def test_x0_aliasing_x_reads_the_live_state():
    cir = mos_chain(8)
    _bound(cir)
    for x, _x0 in _pairs(cir.n, 6):
        _both(cir, x, None, same=True)


def test_shared_nodes_are_read_after_the_earlier_write():
    """An inverter chain (n and p on one drain and one gate) and a cascode
    (two devices on a middle node): the next element reads what the
    previous one wrote, in dict order."""
    cm.default_toolkit = cm.numeric
    cir = SubCircuit()
    cir['vdd'] = VS('vdd', gnd, v=1.8)
    cir['vin'] = VSin('g0', gnd, v=0.9, va=0.1, freq=1e6)
    for k in range(5):
        g, d = ('g0' if k == 0 else f'd{k - 1}'), f'd{k}'
        cir[f'N{k}'] = eh.MosLevel1Hdl(d, g, gnd, gnd)
        cir[f'P{k}'] = eh.MosLevel1PmosHdl(d, g, 'vdd', 'vdd')
    cir['Ca'] = eh.MosLevel1Hdl('d4', 'g0', 'mid', gnd)
    cir['Cb'] = eh.MosLevel1Hdl('mid', 'd2', gnd, gnd)
    cir['rl'] = R('mid', 'vdd', r=1e4)
    _bound(cir, 'N0')
    _bound(cir, 'P0')
    for x, x0 in _pairs(cir.n, 8, seed=2):
        _both(cir, x, x0)
        _both(cir, x, None, same=True)


def _chain_with(mid):
    """Six stages with `mid(name)` inserted after the third, in dict order."""
    cm.default_toolkit = cm.numeric
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['vg'] = VSin('g0', gnd, v=0.9, va=2e-2, freq=1e6)
    for k in range(6):
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
        if k == 2:
            c['mid'] = mid()
    return c


def test_a_hand_written_limiter_mid_chain_stops_and_resumes_the_walk(
        limit_kernel_calls, monkeypatch):
    cir = _chain_with(lambda: Diode('d2', gnd))
    _bound(cir)
    seen = []
    orig = Diode.limit

    def spy(self, x, x0, epar=defaultepar):
        seen.append(self)
        return orig(self, x, x0, epar)
    monkeypatch.setattr(Diode, 'limit', spy)
    x, x0 = _pairs(cir.n, 1)[0]
    xa = x.copy()
    with walking(True):
        cir.limit(xa, x0.copy(), defaultepar)
    assert seen == [cir['mid']] and limit_kernel_calls == []
    xb = x.copy()
    with walking(False):
        cir.limit(xb, x0.copy(), defaultepar)
    assert len(seen) == 2 and len(limit_kernel_calls) == 6
    assert xa.tobytes() == xb.tobytes()
    w = cir.__dict__['_limit_walk']
    assert [k is None for k in w.kerns] == [False, False, False, True, False, False, False]


def test_an_instance_shadow_mid_run_is_called_and_the_walk_resumes(limit_kernel_calls):
    cir = mos_chain(6)
    _bound(cir)
    x, x0 = _pairs(cir.n, 1)[0]
    _both(cir, x, x0)
    el = cir['M2']
    seen = []

    def shadow(xx, xx0, epar=defaultepar):
        seen.append(xx.copy())
        return type(el).limit(el, xx, xx0, epar)
    el.__dict__['limit'] = shadow
    limit_kernel_calls.clear()
    xa = x.copy()
    with walking(True):
        cir.limit(xa, x0.copy(), defaultepar)
    assert len(seen) == 1 and limit_kernel_calls == [el]
    xb = x.copy()
    with walking(False):
        cir.limit(xb, x0.copy(), defaultepar)
    assert len(seen) == 2 and len(limit_kernel_calls) == 7
    assert xa.tobytes() == xb.tobytes()
    del el.__dict__['limit']
    limit_kernel_calls.clear()
    _both(cir, x, x0)
    assert limit_kernel_calls == [] + [cir[f'M{k}'] for k in range(6)]


def test_a_declined_call_is_the_closures_answer(limit_kernel_calls):
    """A NaN at a device's node: its kernel declines (Python's own sort
    order is the answer), the walk stops there, the generated `limit()`
    tries the kernel again and the closure answers."""
    cir = mos_chain(6)
    _bound(cir)
    x, x0 = _pairs(cir.n, 1)[0]
    x[list(cir.nodes).index(cir.get_node('d3'))] = np.nan
    xa = x.copy()
    with walking(True), np.errstate(all='ignore'):
        cir.limit(xa, x0.copy(), defaultepar)
    declined = list(limit_kernel_calls)
    assert 1 <= len(declined) < 6
    limit_kernel_calls.clear()
    xb = x.copy()
    with walking(False), np.errstate(all='ignore'):
        cir.limit(xb, x0.copy(), defaultepar)
    assert len(limit_kernel_calls) == 6
    assert xa.tobytes() == xb.tobytes()


def test_where_the_walk_does_not_serve_the_loop_runs(limit_kernel_calls):
    cir = mos_chain(6)
    _bound(cir)
    x, x0 = _pairs(cir.n, 1)[0]
    ## a temperature that is not one number: the kernel steps aside on
    ## both paths, the closure answers
    epar = types.SimpleNamespace(T=np.array(310.0))
    assert cl.limit_walk(cir, x.copy(), x0.copy(), epar) is None
    _both(cir, x, x0, epar)
    ## one that is: both paths the kernel's
    epar = types.SimpleNamespace(T=np.float64(310.0))
    a = _both(cir, x, x0, epar)
    assert a.tobytes() != _both(cir, x, x0).tobytes()
    ## a read-only state, a strided x0, a complex state
    ro = x.copy()
    ro.flags.writeable = False
    assert cl.limit_walk(cir, ro, x0.copy(), defaultepar) is None
    assert cl.limit_walk(cir, x.copy(), np.repeat(x0, 2)[::2], defaultepar) is None
    assert cl.limit_walk(cir, x.astype(complex), x0.copy(), defaultepar) is None
    assert cl.limit_walk(cir, x.copy()[:-1], x0.copy(), defaultepar) is None
    ## the circuit-level resolution never walks
    called = []
    with contextlib.ExitStack() as st:
        st.enter_context(pytest.MonkeyPatch.context())
        mp = pytest.MonkeyPatch()
        mp.setattr(cl, 'limit_walk', lambda *a: called.append(a) or None)
        mp.setattr(_limiting, 'CIRCUIT_LEVEL', True)
        try:
            with contextlib.suppress(Exception):
                cir.limit(x.copy(), x0.copy(), defaultepar)
        finally:
            mp.undo()
    assert called == []


def test_a_rebound_kernel_is_taken_again(limit_kernel_calls):
    cir = mos_chain(6)
    _bound(cir)
    x, x0 = _pairs(cir.n, 1)[0]
    g0 = _both(cir, x, x0)
    hdl.set_backend('numpy', eh.MosLevel1Hdl)
    try:
        assert type(cir['M0'])._hdl_info.get('_c_limit') is None
        limit_kernel_calls.clear()
        _both(cir, x, x0)
        assert limit_kernel_calls == []
        ## a circuit built while pinned: a capable class that has never had
        ## a kernel at its first walk (the second cut's fault, found by
        ## `test_device_limiter`'s PCNR grids)
        cir2 = mos_chain(4)
        x2, x02 = _pairs(cir2.n, 1, seed=5)[0]
        _both(cir2, x2, x02)
        assert limit_kernel_calls == []
        assert cir2.__dict__['_limit_walk'].any_capable
    finally:
        hdl.set_backend(None, eh.MosLevel1Hdl)
    assert isinstance(type(cir['M0'])._hdl_info.get('_c_limit'), cl.CLimitKernel)
    limit_kernel_calls.clear()
    xa = x.copy()
    with walking(True):
        cir.limit(xa, x0.copy(), defaultepar)
    ## (the new kernel taken by the walk: no per-element call, g0's bytes)
    assert limit_kernel_calls == [] and xa.tobytes() == g0.tobytes()
    assert _both(cir, x, x0).tobytes() == g0.tobytes()


def test_a_nested_circuit_walks_its_own_elements(limit_kernel_calls):
    class _Stage(SubCircuit):
        terminals = ('i', 'o', 'vdd')

        def __init__(self, *args, **kw):
            super().__init__(*args, **kw)
            nn = self.nodenames
            self['r'] = R(nn['vdd'], nn['o'], r=5e3)
            self['m'] = eh.MosLevel1Hdl(nn['o'], nn['i'], gnd, gnd)
            self['m2'] = eh.MosLevel1Hdl(nn['o'], nn['i'], gnd, gnd)
    cm.default_toolkit = cm.numeric
    cir = SubCircuit()
    cir['vdd'] = VS('vdd', gnd, v=1.8)
    cir['vg'] = VSin('g0', gnd, v=0.9, va=2e-2, freq=1e6)
    for k in range(4):
        cir[f's{k}'] = _Stage('g0' if k == 0 else f'd{k - 1}', f'd{k}', 'vdd')
        cir[f'M{k}'] = eh.MosLevel1Hdl(f'e{k}', f'd{k}', gnd, gnd)
        cir[f'r{k}'] = R('vdd', f'e{k}', r=5e3)
    _bound(cir)
    for x, x0 in _pairs(cir.n, 4, seed=3):
        _both(cir, x, x0)
    limit_kernel_calls.clear()
    with walking(True):
        cir.limit(x.copy(), x0.copy(), defaultepar)
    assert limit_kernel_calls == []
    w = cir.__dict__['_limit_walk']
    assert [k is None for k in w.kerns] == [True, False] * 4


def test_a_transient_is_bit_identical():
    from pycircuit.circuit.transient import Transient
    got = {}
    for on in (True, False):
        cir = mos_chain(6)
        _bound(cir)
        with walking(on):
            res = Transient(cir, toolkit=cm.numeric).solve(
                tend=20 * 2e-8, timestep=2e-8, fixed_timestep=True)
        got[on] = np.asarray(res.x, dtype=float)
    assert got[True].tobytes() == got[False].tobytes()


def test_the_switch_is_read_from_the_environment():
    env = dict(os.environ, PYCIRCUIT_HDL_LIMIT_WALK='0')
    code = 'from pycircuit.circuit import _hdl_climit as cl; print(cl.WALK)'
    out = subprocess.run([sys.executable, '-c', code], env=env,
                         capture_output=True, text=True, check=True)
    assert out.stdout.strip() == 'False'
    assert cl.WALK_STATUS in ('c', 'not loaded') or cl.WALK_STATUS.startswith('off')
