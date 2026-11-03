"""The limiter walk (`_hdl_climit.limit_walk`, speed round 6's third commit,
2026-10-03): `SubCircuit.limit`'s loop with every run of C-kernel elements
walked in one C call on the live state -- the loop's bytes in every case
the loop handles (shared nodes, duplicate rows, `x0 is x`, a hand-written
limiter mid-chain, a shadow, a decline, a rebound kernel, a nested
circuit, an iterate gone non-finite), and the loop itself where the walk
does not serve.  Every test runs with the walk alone (one thread) and
ahead on the pool (four, `walk_threads`): a circuit of WALK_POOL_MIN
elements or more takes the pool there, and the bytes must not move; the
pool's own cases are at the end (THE WALK AHEAD)."""
import contextlib
import glob
import os
import re
import subprocess
import sys
import textwrap
import types

import numpy as np
import pytest

from pycircuit.circuit import _cbuild, _limiting, _tran_core, hdl
from pycircuit.circuit import _hdl_climit as cl
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import VS, Diode, R, SubCircuit, VSin, gnd
from pycircuit.circuit.tests.test_hdl_batch import mos_chain


@pytest.fixture(autouse=True, params=[1, 4], ids=['alone', 'pool'])
def walk_threads(request, monkeypatch):
    """The core's thread count, which the walk takes (`_run_walk`,
    `_tran_newton_c._walk_ready`): one, the walk alone; four, the walk ahead
    on the pool where its run is long enough."""
    monkeypatch.setattr(_tran_core, 'THREADS', request.param)
    return request.param


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


def test_a_non_finite_iterate_is_walked_in_c(limit_kernel_calls):
    """An iterate gone non-finite (NaNs and infinities at random rows): on
    CPython 3.14 the walk ranks the NaN keys in C and never stops for them --
    one walk, no element on its per-call path -- with the per-element
    loop's bytes."""
    if not cl.NAN_KEYS:
        pytest.skip('off CPython 3.14 the kernels decline a NaN key')
    from pycircuit.circuit import _paths
    cir = mos_chain()
    _bound(cir)
    rng = np.random.default_rng(224)
    nan_out = 0
    for x, x0 in _pairs(cir.n, 10, seed=3):
        for r in rng.integers(0, cir.n, 3):
            x[r] = rng.choice([np.nan, np.inf, -np.inf])
        before = _paths.snapshot()
        xa = x.copy()
        with walking(True):
            cir.limit(xa, x0.copy(), defaultepar)
        d = _paths.since(before)
        assert limit_kernel_calls == [] and d.get('walk:stopped', 0) == 0
        assert d.get('walk:served', 0) == 1
        xb = x.copy()
        with walking(False):
            cir.limit(xb, x0.copy(), defaultepar)
        limit_kernel_calls.clear()
        assert xa.tobytes() == xb.tobytes()
        nan_out += bool(np.isnan(xa).any())
    assert nan_out >= 5, nan_out


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


def _chain_with(mid, stages=6):
    """`stages` stages with `mid(name)` inserted after the third, in dict
    order."""
    cm.default_toolkit = cm.numeric
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['vg'] = VSin('g0', gnd, v=0.9, va=2e-2, freq=1e6)
    for k in range(stages):
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


def test_a_nan_at_a_node_is_the_walks_on_314_else_the_closures(limit_kernel_calls):
    """A NaN at a device's node: on CPython 3.14 the walk ranks its keys in
    C and never stops; elsewhere its kernel declines (Python's own sort
    order is the answer), the walk stops there, the generated `limit()`
    tries the kernel again and the closure answers -- either way the
    per-element loop's bytes.  (The stop and resume on 3.14:
    `test_a_hand_written_limiter_mid_chain_stops_and_resumes_the_walk`.)"""
    cir = mos_chain(6)
    _bound(cir)
    x, x0 = _pairs(cir.n, 1)[0]
    x[list(cir.nodes).index(cir.get_node('d3'))] = np.nan
    xa = x.copy()
    with walking(True), np.errstate(all='ignore'):
        cir.limit(xa, x0.copy(), defaultepar)
    declined = list(limit_kernel_calls)
    if cl.NAN_KEYS:
        assert declined == []
    else:
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


## ---------------------------------------------------------------------------
## THE WALK AHEAD (walk.c's `walk_range` and `walk_ahead`): every element's
## limiter made on the pool on the state as the walk found it, then the walk
## in its order -- an element on a node an earlier one changed made again --
## the walk alone's bytes in every case.

_WALK_CAPS = _cbuild.constants(cl._WALK_CDEF)


@contextlib.contextmanager
def _threads(n):
    was = _tran_core.THREADS
    _tran_core.THREADS = n
    try:
        yield
    finally:
        _tran_core.THREADS = was


def _alone_and_ahead(cir, x, x0, same=False):
    """`cir.limit` through the walk alone (one thread) and ahead on the pool
    (four), on copies: byte-equal, and the loop's bytes.  Returns how many of
    the state's values the walk changed."""
    got = []
    for n in (1, 4):
        xa = x.copy()
        x0a = xa if same else x0.copy()
        with _threads(n), walking(True):
            assert cir.limit(xa, x0a, defaultepar) is xa
        got.append(xa)
    assert got[0].tobytes() == got[1].tobytes()
    _both(cir, x, x0, same=same)
    return int((got[0].view(np.uint64) != x.view(np.uint64)).sum())


def test_the_walk_ahead_where_writes_reach_later_elements(walk_threads):
    """Eighty stages, each gate on the previous drain: a limited drain is
    the next element's gate, so the next one is made again in the walk --
    a few nodes written (a step at three nodes), and more than
    WALK_DIRTY_MAX (a step at every node: every later element made again)
    -- x0 apart and aliasing x."""
    if walk_threads != 1:
        pytest.skip('(the test sets its own thread counts)')
    cir = mos_chain(80)
    _bound(cir)
    cap = _WALK_CAPS['WALK_DIRTY_MAX']
    rng = np.random.default_rng(232)
    few, many = [], []
    for x, x0 in _pairs(cir.n, 6, seed=5):
        _alone_and_ahead(cir, x, x0)
        _alone_and_ahead(cir, x, None, same=True)
        ## (every node a step of 3 V: every element limits)
        many.append(_alone_and_ahead(cir, x0 + rng.choice([-3.0, 3.0], cir.n), x0))
        xs = x0.copy()
        for r in rng.choice(cir.n, 3, replace=False):
            xs[r] += rng.choice([-3.0, 3.0])
        few.append(_alone_and_ahead(cir, xs, x0))
    assert any(0 < c <= cap for c in few), few
    assert max(many) > cap, many


def _diode_chain(n):
    """`n` diodes in series, each inner node shared by two: above its
    critical voltage `pnjlim`'s answer moves with its input (a FET's laws
    mostly pass a value or clamp it to a bound of the old one), so an element
    on a node an earlier one wrote answers otherwise than on the state the
    walk found."""
    cm.default_toolkit = cm.numeric
    c = SubCircuit()
    c['v'] = VS('n0', gnd, v=0.0)
    for k in range(n):
        c[f'D{k}'] = eh.DiodeSpiceHdl(f'n{k}', f'n{k + 1}' if k < n - 1 else gnd)
    return c


def test_the_walk_ahead_past_its_written_nodes_bound(walk_threads):
    """Two hundred diodes in series, every junction stepped: more than
    WALK_DIRTY_MAX nodes written in one walk, and the elements after them
    read what was written -- the walk alone's bytes, x0 apart and aliasing
    x."""
    if walk_threads != 1:
        pytest.skip('(the test sets its own thread counts)')
    cir = _diode_chain(200)
    _bound(cir, 'D0')
    rng = np.random.default_rng(232)
    most = 0
    for _ in range(6):
        x0 = np.cumsum(rng.uniform(0.0, 0.7, cir.n))[::-1].copy()
        x = x0 + rng.uniform(0.5, 3.0, cir.n) * np.sign(rng.normal(size=cir.n))
        most = max(most, _alone_and_ahead(cir, x, x0))
        _alone_and_ahead(cir, x, None, same=True)
    assert most > _WALK_CAPS['WALK_DIRTY_MAX'], most


def test_the_walk_ahead_on_shared_drains_and_gates(walk_threads):
    """Ten inverters (twenty elements: n and p on one drain and one gate),
    each reading what the previous one wrote."""
    if walk_threads != 1:
        pytest.skip('(the test sets its own thread counts)')
    cm.default_toolkit = cm.numeric
    cir = SubCircuit()
    cir['vdd'] = VS('vdd', gnd, v=1.8)
    cir['vin'] = VSin('g0', gnd, v=0.9, va=0.1, freq=1e6)
    for k in range(10):
        g, d = ('g0' if k == 0 else f'd{k - 1}'), f'd{k}'
        cir[f'N{k}'] = eh.MosLevel1Hdl(d, g, gnd, gnd)
        cir[f'P{k}'] = eh.MosLevel1PmosHdl(d, g, 'vdd', 'vdd')
    _bound(cir, 'N0')
    _bound(cir, 'P0')
    moved = 0
    for x, x0 in _pairs(cir.n, 8, seed=6):
        moved += _alone_and_ahead(cir, x, x0) > 0
        _alone_and_ahead(cir, x, None, same=True)
    assert moved > 0


def test_the_walk_ahead_stops_at_a_loop_element_and_resumes(walk_threads):
    """A loop element (a diode with a hand-written limiter, no kernel) after
    the third of twenty stages: the walk ahead stops there, the loop limits
    it, and the walk resumes ahead from the next (seventeen elements)."""
    if walk_threads != 1:
        pytest.skip('(the test sets its own thread counts)')
    cir = _chain_with(lambda: Diode('d2', gnd), stages=20)
    _bound(cir, 'M0')
    for x, x0 in _pairs(cir.n, 6, seed=7):
        _alone_and_ahead(cir, x, x0)


def test_the_walk_alone_starts_the_pool_s_threads():
    """In a fresh interpreter at four threads, a circuit's `limit` with no
    transient (the core never runs) starts the pool's workers: the walk took
    the pool."""
    script = textwrap.dedent('''
        import os
        import numpy as np
        from pycircuit.circuit.circuit import defaultepar
        from pycircuit.circuit.tests.test_hdl_batch import mos_chain

        def workers():
            out = 0
            for t in os.listdir('/proc/self/task'):
                try:
                    with open(f'/proc/self/task/{t}/comm') as fh:
                        out += fh.read().strip() == 'pcx-core'
                except OSError:
                    pass
            return out

        cir = mos_chain(20)
        x0 = np.zeros(cir.n)
        before = workers()
        cir.limit(x0 + 0.5, x0, defaultepar)
        print(before, workers())
    ''')
    env = dict(os.environ, PYCIRCUIT_CORE_THREADS='4')
    out = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True,
                         timeout=600, env=env, check=False,
                         cwd=os.path.dirname(os.path.dirname(os.path.dirname(cl.__file__))))
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.split() == ['0', '3']


def test_every_caller_clears_the_walk_s_flags_untested():
    """The walk ahead raises the union of its limiters' flags -- more than
    the walk alone raises where an element is made twice or past a stop --
    so no caller may read them: each C caller clears them before its next
    test of any flag, and the C Newton tests none."""
    calls = 0
    for path in sorted(glob.glob(os.path.join(_cbuild.CSRC_DIR, '*.c'))):
        with open(path, encoding='utf-8') as fh:
            text = re.sub(r'/\*.*?\*/', ' ', fh.read(), flags=re.DOTALL)
        for m in re.finditer(r'\bwalk(?:_fn\))?\(s->wF,', text):
            calls += 1
            after = text[m.end():]
            test, clear = after.find('fetestexcept'), after.find('feclearexcept')
            assert test == -1 or (clear != -1 and clear < test), \
                f'{os.path.basename(path)}: a flag tested after the walk before a clear'
    assert calls == 4, calls
