"""Radau's dense stage Newton in C (`_tran_radau_c`, speed round 9, B2): the
undamped coupled Newton of a Radau IIA(3) step in one C call.  Its stages,
the device memo and the source memo are the bytes of the Python Newton's,
with its counts and warnings -- on whole transients and on drawn seeds,
step sizes and charges at one step; every hand-back (a non-finite value, a
floating-point exception, the residual's bound, a singular system, the walk
stopping, the assembly cap) leaves the Python Newton to answer as before;
and the declines."""
import types
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _hdl_climit, _paths, _tran_core, circuit
from pycircuit.circuit import _tran_radau_c as RC
from pycircuit.circuit._tran_radau import _RadauStages
from pycircuit.circuit.elements import C, Diode, R, gnd
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_newton_c import _gp_chain, _mos_chain
from pycircuit.circuit.tests.test_psp_limit_c import _stage
from pycircuit.circuit.transient import Transient

DBL_MAX = float(np.finfo(float).max)


def _on():
    if RC.driver() is None:
        pytest.skip(f'the radau C is off: {RC.STATUS}')


def _memo(tr):
    """The device memo's two generations: keys and every record's bytes."""
    memo = getattr(tr, '_dev_memo', None) or ({}, {})
    return [sorted((k, sorted((kk, np.asarray(v).tobytes()) for kk, v in rec.items()))
                   for k, rec in g.items()) for g in memo]


def _umemo(tr):
    um = tr.__dict__.get('_u_memo') or {}
    return sorted((repr(k), np.asarray(v).tobytes()) for k, v in um.items())


def _run(build, on, make=None, **kw):
    old = RC.ENABLED
    RC.ENABLED = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                           **(make or {}))
            before = _paths.snapshot()
            res = tr.solve(**kw)
            d = _paths.since(before)
    finally:
        RC.ENABLED = old
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), _memo(tr)), d


CASES = {
    'mos': (_mos_chain, {}),
    'mos-caps': (lambda: _mos_chain(cap=1e-13), {}),
    'gp-dense': (_gp_chain, {'radau_transform': False}),
    'psp-dense': (_stage, {'radau_transform': False}),
}


@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
@pytest.mark.parametrize('case', list(CASES))
def test_a_radau_transient_is_the_same_with_the_c_off(case, fixed):
    _on()
    build, make = CASES[case]
    kw = {'tend': 6e-7, 'timestep': 2e-8}
    if fixed:
        kw['fixed_timestep'] = True
    a, _da = _run(build, False, make, **dict(kw))
    b, db = _run(build, True, make, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert a[4] == b[4], 'the device memo moved'
    assert db.get('radau_c:served', 0) > 0, db


## -- one step's Newton, captured: drawn seeds, steps and charges ---------------

def _captured(build=_mos_chain, steps=4, **make):
    """A transient whose last step's coupled solve was kept: the step's
    `(ctx, _stage_newton, _block_residual, seed0)` and copies of the device
    and source memos as its Newton found them."""
    tr = Transient(build(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(), **make)
    got = {}
    real = tr._coupled_stage_solver

    def spy(x0, t, pf):
        out = real(x0, t, pf)
        got.update(out=out, u=dict(tr.__dict__['_u_memo']),
                   dev=tuple({k: dict(v) for k, v in g.items()} for g in tr._dev_memo))
        return out
    tr._coupled_stage_solver = spy
    tr.solve(tend=steps * 2e-8, timestep=2e-8, fixed_timestep=True)
    del tr._coupled_stage_solver
    return tr, got


_CAP = {}


def _cap():
    if 'mos' not in _CAP:
        _CAP['mos'] = _captured()
    return _CAP['mos']


def _attempt(tr, got, seed, on):
    """`_stage_newton(seed)` from the captured step's memos, the C on or off:
    the stages' bytes (or the exception), the device and source memos, the
    source counts and the warnings; and every count."""
    sn = got['out'][1]
    tr.__dict__['_u_memo'] = dict(got['u'])
    tr._dev_memo = tuple({k: dict(v) for k, v in g.items()} for g in got['dev'])
    old = RC.ENABLED
    RC.ENABLED = on
    before = _paths.snapshot()
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            try:
                out = [np.asarray(y).tobytes() for y in sn(seed)]
            except Exception as e:                             # noqa: BLE001
                out = (type(e).__name__, str(e))
    finally:
        RC.ENABLED = old
    d = _paths.since(before)
    src = {k: v for k, v in d.items() if k.startswith(('umemo:', 'u:', 'src.u:'))}
    return (out, _memo(tr), _umemo(tr), src, sorted(str(w.message) for w in W)), d


def _both(seed, h=None, qn=None):
    """The captured Newton from `seed` with the C off and on (the step `h`
    and the entering charge `qn` set on its context for the pair): the two
    outcomes and the C run's counts."""
    _on()
    tr, got = _cap()
    ctx = got['out'][0]
    h0, qn0 = ctx.h, ctx.qn
    if h is not None:
        ctx.h = h
    if qn is not None:
        ctx.qn = qn
    try:
        a, da = _attempt(tr, got, seed, False)
        b, db = _attempt(tr, got, seed, True)
    finally:
        ctx.h, ctx.qn = h0, qn0
    if not db.get('radau_c:served'):
        ## a hand-back leaves no trace: every count the Python's but the C's own
        own = lambda d: {k: v for k, v in d.items() if not k.startswith('radau_c:')}
        assert own(db) == own(da), (da, db)
    return a, b, db


def _seed(dx=None):
    seed0 = _cap()[1]['out'][3]
    Y = [np.array(y, dtype=float) for y in seed0]
    if dx is not None:
        for j, d in dx.items():
            Y[j] = Y[j] + d
    return Y


def _counts(db):
    return {k: v for k, v in db.items() if k.startswith('radau_c:')}


def test_the_captured_step_is_served_and_the_pythons():
    a, b, db = _both(_seed())
    assert a == b
    assert _counts(db).get('radau_c:served') == 1, db


def _bail(db, why):
    c = _counts(db)
    assert c.get('radau_c:bail:' + why) == 1 and 'radau_c:served' not in c, c


def test_every_hand_back_answers_as_the_python():
    """Each way the C hands the solve back, planted: the Python Newton's
    answer, exception, warnings and memos, and the hand-back counted."""
    n = len(_seed()[0])
    ## a non-finite stage: `passes` declines, the Python evaluates
    bad = _seed()
    bad[1][2] = np.nan
    a, b, db = _both(bad)
    assert a == b
    _bail(db, 'nonfinite')
    ## an overflow in the residual: numpy warns
    qn = np.where(np.arange(n) % 2, DBL_MAX, -DBL_MAX)
    a, b, db = _both(_seed(), h=1e300, qn=qn)
    assert a == b and any('overflow' in w for w in a[4]), a[4]
    _bail(db, 'flags')
    ## a residual past `sum |R|`'s bound (finite, no exception)
    qn = np.full(n, 1e305)
    a, b, db = _both(_seed(), qn=qn)
    assert a == b
    _bail(db, 'rbound')
    ## a singular system (no step: J = C, singular without capacitors)
    a, b, db = _both(_seed(), h=0.0)
    assert a == b and a[0][0] == 'LinAlgError', a[0]
    _bail(db, 'dgesv')


def test_the_assembly_cap_hands_back(monkeypatch):
    _on()
    a, _da = _run(_mos_chain, False, tend=4e-7, timestep=2e-8, fixed_timestep=True)
    monkeypatch.setattr(RC, 'MAXA', 1)
    b, db = _run(_mos_chain, True, tend=4e-7, timestep=2e-8, fixed_timestep=True)
    assert a == b
    assert db.get('radau_c:bail:assemblies', 0) > 0 and 'radau_c:served' not in db, db


def test_a_stopped_walk_hands_back_and_is_taken_again():
    """The walk's first entry planted as a stop (its kernel's address 0, so
    every setup takes it): the C hands back (`walkstop`), the Python's walk
    runs the loop's own statement for that element (the same law), and once
    the address is back the next C call takes the tables again and is
    served."""
    _on()
    tr, got = _cap()
    rec = tr.__dict__['_radau_c']
    assert rec.walk is not None and rec.walk_ok
    a, _da = _attempt(tr, got, _seed(), False)
    ## (the walk the circuit holds now: both paths' setups take it)
    w = _hdl_climit._walk_for(tr.cir)
    addr0 = w.addr[0]
    w.addr[0] = 0
    w.F[0] = 0
    try:
        b, db = _attempt(tr, got, _seed(), True)
    finally:
        w.addr[0] = addr0
        w.F[0] = addr0
    assert a[0] == b[0] and a[1:3] == b[1:3] and a[4] == b[4]
    _bail(db, 'walkstop')
    assert db.get('walk:stopped', 0) > 0, db
    assert rec.walk_ok is None
    c, dc = _attempt(tr, got, _seed(), True)
    assert c == a and _counts(dc).get('radau_c:served') == 1, dc


SCALES = (0.0, 1e-12, 1e-6, 1e-3, 0.05, 0.4, 3.0, 40.0, 1e3, 1e30, 1e150, 1e300)


@settings(deadline=None, max_examples=60)
@given(data=st.data())
def test_drawn_seeds_steps_and_charges_are_the_pythons(data):
    _on()
    seed = _seed()
    n = len(seed[0])
    for j in range(3):
        sc = data.draw(st.sampled_from(SCALES))
        if sc:
            for i in data.draw(st.lists(st.integers(0, n - 1), min_size=1, max_size=4)):
                seed[j][i] += sc * data.draw(st.sampled_from((1.0, -1.0, 0.37, -2.5)))
        if data.draw(st.integers(0, 40)) == 0:
            seed[j][data.draw(st.integers(0, n - 1))] = data.draw(
                st.sampled_from((np.nan, np.inf, -np.inf, -0.0)))
    h0 = _cap()[1]['out'][0].h
    h = h0 * data.draw(st.sampled_from((1.0, 1.0, 1.0, 1e-6, 1e-3, 30.0, 1e6, 1e300)))
    qn = None
    if data.draw(st.integers(0, 5)) == 0:
        qn = np.array(_cap()[1]['out'][0].qn, dtype=float)
        qn[data.draw(st.integers(0, n - 1))] = data.draw(st.sampled_from(
            (1e-12, -3.0, 1e300, -DBL_MAX, 1e305)))
    a, b, _db = _both(seed, h=h, qn=qn)
    assert a == b


def test_drawn_seeds_reach_the_c_and_its_hand_backs():
    """Not vacuous: on a walk of drawn seeds the C serves most solves, over
    one iteration and several, and hands some back."""
    _on()
    rng = np.random.default_rng(5)
    served = handed = 0
    iters = set()
    for k in range(60):
        seed = _seed()
        sc = (1e-9, 1e-3, 0.05, 0.5, 2.0, 1e3)[k % 6]
        for j in range(3):
            seed[j][1:] += sc * rng.standard_normal(len(seed[j]) - 1)
        a, b, db = _both(seed)
        assert a == b, (k, sc)
        c = _counts(db)
        if c.get('radau_c:served'):
            served += 1
            iters.add(c['radau_c:assemblies'])
        else:
            handed += 1
    assert served > 30 and handed > 0 and len(iters) > 2, (served, handed, iters)


## -- declines -------------------------------------------------------------------

def _declines(why, build=_mos_chain, make=None, setup=None, **kw):
    _on()
    kw = {'tend': 2e-7, 'timestep': 2e-8, 'fixed_timestep': True, **kw}
    a, _da = _run(build, False, make, **dict(kw))
    if setup is not None:
        undo = setup()
    try:
        b, db = _run(build, True, make, **dict(kw))
    finally:
        if setup is not None and undo is not None:
            undo()
    assert a == b
    assert db.get('radau_c:' + why, 0) > 0 and 'radau_c:served' not in db, db


def _patch(owner, name, value):
    old = owner.__dict__[name] if isinstance(owner, type) else getattr(owner, name)

    def undo():
        setattr(owner, name, old)
    setattr(owner, name, value)
    return undo


def test_a_patched_piece_of_the_newton_declines():
    real = _RadauStages._coupled_stage_system

    def wrapped(self, *a, **k):
        return real(self, *a, **k)
    _declines('patched', setup=lambda: _patch(_RadauStages, '_coupled_stage_system', wrapped))
    real_solve = np.linalg.solve
    calls = []

    def spy(*a, **k):
        calls.append(1)
        return real_solve(*a, **k)
    _declines('patched', setup=lambda: _patch(np.linalg, 'solve', spy))
    assert calls
    real_passes = _tran_core.passes
    _declines('patched', setup=lambda: _patch(
        _tran_core, 'passes', lambda *a, **k: real_passes(*a, **k)))


def test_a_subclass_declines():
    _on()

    class Sub(Transient):
        pass
    old = RC.ENABLED
    RC.ENABLED = True
    try:
        tr = Sub(_mos_chain(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator())
        before = _paths.snapshot()
        tr.solve(tend=1e-7, timestep=2e-8, fixed_timestep=True)
        d = _paths.since(before)
    finally:
        RC.ENABLED = old
    assert d.get('radau_c:class', 0) > 0 and 'radau_c:served' not in d, d


def test_a_stateful_limiter_and_a_bypass_decline():
    def diode():
        c = _mos_chain(2)
        c.add_node('dd')
        c['rd'] = R('vdd', 'dd', r=1e3)
        c['D'] = Diode('dd', gnd)
        c['cd'] = C('dd', gnd, c=1e-12)
        return c
    _declines('stateful', build=diode)
    _declines('bypass', make={'bypass': True, 'bypasstol': 1e-9})


def test_the_switches_decline():
    _on()
    _a, da = _run(_mos_chain, False, tend=2e-7, timestep=2e-8, fixed_timestep=True)
    assert da.get('radau_c:off', 0) > 0 and 'radau_c:served' not in da, da
    _declines('core', setup=lambda: _patch(_tran_core, 'CORE_PASSES', False))


def test_a_provided_function_declines():
    _on()
    tr, got = _cap()
    ctx, _sn, _br, seed0 = got['out']
    before = _paths.snapshot()
    assert RC.solve(tr, ctx, seed0, lambda t: 0.0, lambda t: 0.0, [], True, 1e-3, 1e-6,
                    10) is None
    assert _paths.since(before) == {'radau_c:pf': 1}


def test_inputs_of_another_kind_decline():
    _on()
    tr, got = _cap()
    ctx, _sn, _br, seed0 = got['out']
    src = ctx.src
    tr.__dict__['_u_memo'] = dict(got['u'])

    def declined(why, **kw):
        args = {'tr': tr, 'ctx': ctx, 'seed': seed0, 'src': src, 'provided_function': None,
                'lims': [], 'nobypass': True, 'reltol': 1e-3, 'abstol': 1e-6, 'maxit': 10}
        args.update(kw)
        before = _paths.snapshot()
        assert RC.solve(**args) is None
        d = {k: v for k, v in _paths.since(before).items() if not k.startswith('once:')}
        assert d == {'radau_c:' + why: 1}
    declined('stateful', lims=[object()])
    declined('bypass', nobypass=False)
    declined('tol', reltol='1e-3')
    declined('tol', reltol=True)
    declined('seed', seed=[seed0[0], seed0[1]])
    declined('seed', seed=[s[:-1] for s in seed0])
    bad = types.SimpleNamespace(**vars(ctx))
    bad.h = np.float32(ctx.h)
    declined('h', ctx=bad)
    bad = types.SimpleNamespace(**vars(ctx))
    bad.qn = ctx.qn[:-1]
    declined('ctx', ctx=bad)
    tr.__dict__['_memo_put'] = tr._memo_put
    try:
        declined('shadow_tr')
    finally:
        del tr.__dict__['_memo_put']
