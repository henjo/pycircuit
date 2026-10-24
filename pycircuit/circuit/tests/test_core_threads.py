"""The evaluation core's device kernels on threads (core.c's pool,
`_tran_core.THREADS`): a call's fused-group elements made by the caller
and the pool's workers, the bins summed after in their order -- every
output byte and every floating-point flag as one thread leaves them, on
drawn calls (the kernels pushed to overflow and invalid results too) and
on whole transients (Radau's cost transform, gear2); two transients
calling at once (one on the pool, the other alone); a forked child, which
starts workers of its own; and the workers there at all."""
import ctypes
import gc
import os
import select
import signal
import threading
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_core, circuit
from pycircuit.circuit.integrator import Gear2Integrator, RadauIIA3Integrator
from pycircuit.circuit.tests.test_loop_c import _mos, _pulsed
from pycircuit.circuit.transient import Transient

pytestmark = pytest.mark.skipif(_tran_core.driver() is None,
                                reason='the evaluation core did not build')

#: (glibc's x86-64 flags: invalid, divide-by-zero, overflow, underflow,
#: inexact)
FE_ALL = 0x3d
_LIBM = ctypes.CDLL('libm.so.6')
_LIBM.fetestexcept.argtypes = _LIBM.feclearexcept.argtypes = [ctypes.c_int]
#: (the loops' own counts: where an attempt is made, never what it makes)
_OWN = ('once:', 'sources_c:', 'loop_c:', 'pred_c:', 'branch_c:', 'stageloop_c:', 'tcloop_c:',
        'radau_c:held')


def _core(stages=80):
    """A stepped transient's evaluation core: `stages` MOSFETs in its fused
    group, past the pool's minimum."""
    tr = Transient(_mos(stages, cap=1e-13), toolkit=circuit.numeric,
                   integrator=RadauIIA3Integrator())
    tr.solve(tend=4e-8, timestep=2e-8, fixed_timestep=True)
    return tr, _tran_core.core_for(tr)


def _call(core, x, bits, nthreads):
    """One core call (the passes alone) at `x` on `nthreads`: its outputs'
    bytes and the flags it raised in the caller."""
    ffi, cfn, _ = _tran_core.driver()
    n = core.n
    C, q, d = np.zeros((n, n)), np.zeros(n), np.zeros(max(n * n, 1))
    fb = ffi.from_buffer
    core.cs.nthreads = nthreads
    try:
        _LIBM.feclearexcept(FE_ALL)
        cfn(core.cs, fb('double *', x), 300.0, bits, 4, 0.0, 0.0, 0.0, 1.0, 0.0,
            fb('double *', d), fb('double *', d), fb('double *', d), fb('double *', d),
            fb('double *', d), fb('double *', C), fb('double *', q), fb('double *', d),
            fb('double *', d), fb('double *', d), fb('double *', d))
        flags = _LIBM.fetestexcept(FE_ALL)
    finally:
        core.cs.nthreads = 1
        _LIBM.feclearexcept(FE_ALL)
    return (C.tobytes(), q.tobytes(), core.arrays['binG'].tobytes(),
            core.arrays['bini'].tobytes(), core.arrays['zst'].tobytes(), flags)


def test_a_threaded_call_is_the_one_threaded_one():
    """Drawn states, small to absurd (the kernels' exp and pow overflowing,
    invalid results inside), each pass mask: on four threads every output
    byte -- the matrices, the bins, the fused staging -- and the flags
    raised as on one."""
    _tr, core = _core()
    assert core.nz and int(core.cs.nz) > 0
    rng = np.random.default_rng(220)
    n = core.n
    seen = set()
    for scale in (0.0, 0.5, 2.0, 50.0, 1e3, 1e30, 1e300):
        for bits in (0, 1, 2, 4, 7, 15):
            x = scale * rng.standard_normal(n)
            one, four = _call(core, x, bits, 1), _call(core, x, bits, 4)
            assert one == four, (scale, bits)
            seen.add(one[-1])
    ## (the drawn calls raise more than inexact: overflow and invalid among them)
    assert any(f & 0x08 for f in seen) and any(f & 0x01 for f in seen), seen


def test_a_flag_one_element_raises_on_a_worker_is_the_callers():
    """One element alone overflowing (the last stage's drain absurd, every
    other node quiet): whichever thread takes it -- the chunks are taken as
    they come, so over forty calls the workers do -- its overflow is raised
    in the caller, as on one thread."""
    tr, core = _core()
    n = core.n
    x = np.full(n, 0.3)
    x[tr.cir.get_node_index(tr.cir.get_node('d79'))] = 1e300
    one = _call(core, x, 15, 1)
    assert one[-1] & 0x08, one[-1]
    for _ in range(40):
        assert _call(core, x, 15, 4) == one


def _run(build, integ, threads, make=None, **kw):
    """One transient on `threads`: its answer, statistics, warnings and
    counts."""
    gc.collect()
    old = _tran_core.THREADS
    _tran_core.THREADS = threads
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build(), toolkit=circuit.numeric, integrator=integ(), **(make or {}))
            before = _paths.snapshot()
            res = tr.solve(**kw)
            counts = _paths.since(before)
    finally:
        _tran_core.THREADS = old
    st = tr.statistics
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            {k: getattr(st, k) for k in st.__slots__ if 'seconds' not in k},
            sorted(str(w.message) for w in W),
            {k: v for k, v in counts.items() if not k.startswith(_OWN)})


def _bytes(build, **kw):
    """A transform transient's answer as `_tran_core.THREADS` stands,
    touching no global -- for runs made at once (their caller sets the
    threads and the warnings around them)."""
    tr = Transient(build(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                   radau_transform=True)
    res = tr.solve(**kw)
    return np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes()


@pytest.mark.parametrize('stages', [60, 20])
@pytest.mark.parametrize('integ, make', [(RadauIIA3Integrator, {'radau_transform': True}),
                                         (RadauIIA3Integrator, None), (Gear2Integrator, None)],
                         ids=['radau-transform', 'radau', 'gear2'])
def test_a_transient_on_threads_is_the_one_threaded_one(integ, make, stages):
    """The pulsed chain of 60 stages and of 20 (just past the pool's
    minimum), adaptive: on four threads the bytes, the statistics, the
    warnings and every count but where the loops made their attempts (a
    process's history moves that) as on one."""
    kw = {'tend': 4e-8, 'timestep': 2.5e-10}
    one = _run(lambda: _pulsed(stages), integ, 1, make, **kw)
    four = _run(lambda: _pulsed(stages), integ, 4, make, **kw)
    assert one == four


def test_two_transients_at_once_are_their_own(monkeypatch):
    """Two transients on four threads each, from two Python threads at once
    (the core calls release the interpreter): one call at a time on the
    pool, the other's made alone -- each run as when it runs by itself."""
    kw = {'tend': 4e-8, 'timestep': 2.5e-10}
    builds = (lambda: _pulsed(60), lambda: _mos(70, cap=1e-13))
    ref = [_run(b, RadauIIA3Integrator, 1, {'radau_transform': True}, **kw) for b in builds]
    got = [None, None]
    gc.collect()
    ## (the threads and the warnings set for both runs here: each setting and
    ## restoring them would race)
    monkeypatch.setattr(_tran_core, 'THREADS', 4)

    def go(k):
        got[k] = _bytes(builds[k], **kw)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ts = [threading.Thread(target=go, args=(k,)) for k in range(2)]
        for t in ts:
            t.start()
        for t in ts:
            t.join()
    assert got[0] == ref[0][:2] and got[1] == ref[1][:2]


@pytest.mark.skipif(not hasattr(os, 'fork'), reason='no fork')
@pytest.mark.filterwarnings('ignore:This process .* is multi-threaded')
def test_a_forked_child_starts_workers_of_its_own():
    """The pool running in the parent, a fork's child runs a threaded
    transient (the parent's workers are not in it: it starts its own) --
    the parent's bytes.  (A child waiting for workers it has not would hang:
    killed after two minutes, a failure.)"""
    kw = {'tend': 4e-8, 'timestep': 2.5e-10}
    ref = _run(lambda: _pulsed(60), RadauIIA3Integrator, 4, {'radau_transform': True}, **kw)
    r, w = os.pipe()
    pid = os.fork()
    if pid == 0:
        try:
            got = _run(lambda: _pulsed(60), RadauIIA3Integrator, 4, {'radau_transform': True}, **kw)
            os.write(w, b'1' if got[:2] == ref[:2] else b'0')
        finally:
            os._exit(0)
    os.close(w)
    ready, _, _ = select.select([r], [], [], 120.0)
    ans = os.read(r, 1) if ready else b''
    os.close(r)
    if not ready:
        os.kill(pid, signal.SIGKILL)
    os.waitpid(pid, 0)
    assert ans == b'1', 'the child hung' if not ready else ans


def _workers():
    """This process's pool workers: its threads named "pcx-core"."""
    out = 0
    for t in os.listdir('/proc/self/task'):
        try:
            with open(f'/proc/self/task/{t}/comm') as f:
                out += f.read().strip() == 'pcx-core'
        except OSError:
            pass
    return out


@pytest.mark.skipif(not os.path.isdir('/proc/self/task'), reason='no /proc')
def test_the_workers_run():
    """Alive: a threaded transient has the pool's three workers running in
    this process (started by it, or by an earlier test's)."""
    _run(lambda: _pulsed(60), RadauIIA3Integrator, 4, {'radau_transform': True},
         tend=4e-8, timestep=2.5e-10)
    assert _workers() >= 3
