"""Exact invalidation (`_watch`, speed round 9, stage 4): one counter that
any change to a watched dict moves, and the readiness checks it stamps --
the evaluate core's probe and lookup, the C Newton's limiter walk.  While
the counter stands a stamped check passes unchecked; every change it
stands for, and every one no watcher sees (an instance's `__dict__`
reassigned, a class attribute), checks again -- with the answers of the
checks run in full."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _hdl_cbackend as cb
from pycircuit.circuit import _paths, _tran_core, _watch, circuit, compact
from pycircuit.circuit.tests.test_newton_c import _gp_chain, _mos_chain
from pycircuit.circuit.tests.test_psp_limit_c import _chain
from pycircuit.circuit.transient import Transient


def _on():
    if _watch.now() < 0:
        pytest.skip(f'watching is off: {_watch.STATUS}')


def test_the_counter_moves_on_a_change_and_not_on_a_read():
    _on()
    d = {'a': 1}
    before = _watch.now()
    assert _watch.arm([d], before) == before
    for change in (lambda: d.__setitem__('a', 2), lambda: d.__setitem__('b', 3),
                   lambda: d.pop('b'), lambda: d.update(c=4), lambda: d.clear()):
        e0 = _watch.now()
        _ = d.get('a'), 'a' in d, len(d), list(d)
        assert _watch.now() == e0, 'a read moved the counter'
        change()
        assert _watch.now() > e0, 'a change did not move it'


def test_a_change_during_the_check_refuses_the_stamp():
    _on()
    d = {'a': 1}
    before = _watch.now()
    _watch.arm([d], before)
    d['a'] = 2                          # (the check wrote into what it read)
    assert _watch.arm([d], before) == -1
    assert _watch.arm([d], _watch.now()) >= 0
    assert _watch.arm(['not a dict'], _watch.now()) == -1


def _served(build=_mos_chain, **kw):
    """A transient run a few steps: its core stamped."""
    tr = Transient(build(), toolkit=circuit.numeric)
    tr.solve(tend=1e-7, timestep=2e-8, fixed_timestep=True, **kw)
    core = _tran_core.core_for(tr)
    assert core is not None
    T = core.probe(tr.epar)
    assert T.__class__ is not str and core.stamp >= 0, (T, core.stamp, _watch.STATUS)
    return tr, core


def test_a_stamped_core_sees_every_change_its_check_reads(monkeypatch):
    """Each of the probe's reasons, from a stamped core: an instance shadow
    set, an element's `__dict__` replaced by one holding a shadow, a class's
    pass patched, the class detached, the pack dropped -- the probe answers
    as its full check does, and serves again once the change is undone."""
    _on()
    tr, core = _served()
    el = core.uniq_els[0]
    cls = type(el)
    ## an instance shadow (the element's dict: watched)
    el.__dict__['G'] = lambda *a, **k: None
    assert core.probe(tr.epar) == 'shadow'
    del el.__dict__['G']
    assert core.probe(tr.epar).__class__ is not str
    ## the dict REPLACED (no watcher sees it): the identity check does
    d = dict(el.__dict__)
    d['q'] = lambda *a, **k: None
    old = el.__dict__
    object.__setattr__(el, '__dict__', d)
    try:
        assert core.probe(tr.epar) == 'shadow'
    finally:
        object.__setattr__(el, '__dict__', old)
    assert core.probe(tr.epar).__class__ is not str
    ## a pass patched on the class (a class attribute: the code check)
    real = cls.G
    monkeypatch.setattr(cls, 'G', lambda self, *a, **k: real(self, *a, **k))
    assert core.probe(tr.epar) == 'generated'
    monkeypatch.undo()
    assert core.probe(tr.epar).__class__ is not str
    ## the class detached (its info dict)
    info = cls._hdl_info
    monkeypatch.setitem(info, '_c_bound', False)
    assert core.probe(tr.epar) == 'unbound'
    monkeypatch.undo()
    assert core.probe(tr.epar).__class__ is not str
    ## a parameter set: the pack dropped (the element's dict) and rebuilt,
    ## the core's mirror re-pointed at the new pack
    pack0 = el.__dict__['_hdl_cp']
    el.ipar.lambd = el.ipar.lambd * 1.5 if el.ipar.lambd else 0.01
    el.update_iparv()
    assert core.probe(tr.epar).__class__ is not str
    assert el.__dict__['_hdl_cp'] is not pack0
    assert any(m is el.__dict__['_hdl_cp'] for m in core.mirror)


def test_a_rebound_kernel_and_the_fuse_switch_are_seen(monkeypatch):
    _on()
    tr, core = _served()
    bt = core.batches[0]
    fn = bt.info['funcs'][bt.m]
    old = fn.__dict__['_hdl_c']
    new = cb.kernel_for(fn, old.nx)[0]
    monkeypatch.setattr(fn, '_hdl_c', new)
    assert core.probe(tr.epar).__class__ is not str
    assert bt.kern is new
    monkeypatch.undo()
    assert core.probe(tr.epar).__class__ is not str and bt.kern is old
    if core.nz:
        assert core.z_fn[0]
        monkeypatch.setattr(cb, 'FUSE', False)       # (the backend's module dict)
        assert core.probe(tr.epar).__class__ is not str and not core.z_fn[0]
        monkeypatch.undo()
        assert core.probe(tr.epar).__class__ is not str and core.z_fn[0]


def test_the_core_lookup_follows_the_plan(monkeypatch):
    """An element added to the circuit (its elements dict), or a parameter
    moving `ParameterDict`'s epoch: the lookup asks `_plan_for` again and
    the core is rebuilt for the new plan."""
    _on()
    tr, core = _served()
    assert _tran_core.core_for(tr) is core
    from pycircuit.circuit.elements import R
    tr.cir['rx'] = R('d0', 'd1', r=1e4)
    core2 = _tran_core.core_for(tr)
    assert core2 is not core and core2 is not None
    assert _tran_core.core_for(tr) is core2
    el = next(e for e in tr.cir.elements.values() if type(e).__name__ == 'R')
    el.ipar.r = el.ipar.r * 2.0
    el.update_iparv()
    core3 = _tran_core.core_for(tr)
    assert core3 is not core2


def test_the_walk_sees_a_class_vlimit(monkeypatch):
    """PSP's hand-written limiter twin is its class's while `vlimit` is the
    baked value -- a class attribute no dict watcher of the walk's sees:
    the stamped walk asks the twin every call."""
    _on()
    from pycircuit.circuit import _tran_newton_c
    tr = Transient(_chain(3, va=0.6), toolkit=circuit.numeric)
    before = _paths.snapshot()
    tr.solve(tend=2e-7, timestep=1e-7, fixed_timestep=True)
    assert _paths.since(before).get('newton_c:served', 0) > 0
    rec = tr.__dict__['_newton_c']
    assert rec.walk_stamp >= 0 and rec.walk_hwk
    assert _tran_newton_c._walk_ready(tr, rec)
    monkeypatch.setattr(compact.PspMosLongChannel, 'vlimit',
                        compact.PspMosLongChannel.vlimit * 2.0)
    assert not _tran_newton_c._walk_ready(tr, rec)


def _run(build, watch, monkeypatch, **kw):
    if not watch:
        monkeypatch.setattr(_watch, 'EPOCH', None)
    make = kw.pop('make', {})
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        warnings.simplefilter('ignore', ResourceWarning)
        tr = Transient(build(), toolkit=circuit.numeric, **make)
        res = tr.solve(**kw)
    monkeypatch.undo()
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), _paths.since(before))


@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
@pytest.mark.parametrize('build', [_mos_chain, _gp_chain, lambda: _chain(4, va=0.6)],
                         ids=['mos', 'gp', 'psp'])
def test_watched_and_unwatched_runs_are_the_same(build, fixed, monkeypatch):
    kw = {'tend': 1e-6, 'timestep': 2e-8}
    if fixed:
        kw['fixed_timestep'] = True
    a = _run(build, False, monkeypatch, **dict(kw))
    b = _run(build, True, monkeypatch, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    ## (a `once:` count lands on whichever run resolves a class first)
    def kept(d):
        return {k: v for k, v in d.items() if not k.startswith('once:')}
    assert kept(a[4]) == kept(b[4]), 'a path count moved'


def test_a_checker_stops_arming_past_the_cap(monkeypatch):
    """Dicts that change every call (here: a key written into an element's
    dict before each probe) re-stamp at most `MAX_ARMS` times; past that
    the core checks in full every time, unstamped."""
    _on()
    tr, core = _served()
    el = core.uniq_els[0]
    for k in range(_watch.MAX_ARMS + 4):
        el.__dict__['_churn'] = k
        assert core.probe(tr.epar).__class__ is not str
    assert core.stamp == -1 and core.arms > _watch.MAX_ARMS
