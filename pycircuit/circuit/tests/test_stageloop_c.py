"""The plain run's Radau steps in C (`_tran_stageloop_c`, csrc/loop.c's
stage family): a Radau IIA(3) transient, the cost transform off, with the C
stage loop on is the run with it off -- every output byte, time, statistic
and warning, the `_paths` tallies the Python attempts make, every attribute
the attempts leave on the transient (`transient._RUN_STATE`, the device
memo's two generations among them, its workspaces aside), the family's own
(its running reference, retry record, carried order) and the stage KLU's
counts -- adaptive and fixed, dense and on the stage KLU, under each
reference, retry order and `lte_algebraic` mode; an entry that fills its
outputs; its stops (a screen that fires, the estimate's solve or the
sources failing in the C, the budget spent: each after a converged Newton,
its result held for the Python's repeat; a Newton handed back); its
declines; and alive: a law or the estimate a hair off is another run.
"""
import inspect
import itertools
import math
import warnings

import numpy as np
import pytest

from pycircuit.circuit import (
    _paths,
    _tran_loop_c,
    _tran_radau_c,
    _tran_stageloop_c,
    circuit,
)
from pycircuit.circuit import transient as trmod
from pycircuit.circuit.elements import VPulse, VSin, gnd
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_loop_c import (
    _NOT_STATE,
    _OWN_COUNTS,
    _gp,
    _mos,
    _norm,
    _psp_stage,
    _pulsed,
    _pulsed_load,
)
from pycircuit.circuit.transient import Transient, _StageSteps

pytestmark = pytest.mark.skipif(_tran_loop_c.driver() is None,
                                reason='the stepping loop unit did not build')

## (the predictor's workspace for a stage step's seeds: each call copies the
## history in and reads its output, nothing carried)
_NOT_STATE_S = _NOT_STATE | {'_pred_c'}
## (counts the stage loop's own machinery makes: its declines and bails, the
## holds the Python's repeats take)
_OWN = _OWN_COUNTS + ('stageloop_c:', 'radau_c:held')
#: (the cost transform off: 'auto' leaves it on for the PSP's compiled Jacobian)
_PLAIN = {'radau_transform': False}
#: (a pulse's corners: rejections at a point, landings, order drops)
_PULSE = {'tend': 4e-8, 'timestep': 2.5e-10}
#: (... and, this tight, two at a point: the retry's order measured)
_TIGHT = {'reltol': 1e-7}


def _edge():
    """The pulsed chain behind a ten-picosecond edge: the measured order
    clamped at its floor."""
    c = _mos(6, cap=1e-13)
    c['vg'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=3e-9, tr=1e-11, tf=1e-11, pw=6e-9, per=2e-8)
    return c


def _run(build, on, monkeypatch, make=None, setup=None, **kw):
    """One run with the C stage loop on or off: `(what it left, the paths
    counted)`."""
    fams = []
    init = _StageSteps.__init__

    def kept(self, *a, **k):
        fams.append(self)
        init(self, *a, **k)
    with monkeypatch.context() as mp:
        mp.setattr(_tran_stageloop_c, 'ENABLED', on)
        mp.setattr(_StageSteps, '__init__', kept)
        before = _paths.snapshot()
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                           **dict(_PLAIN, **(make or {})))
            if setup:
                setup(tr)
            res = tr.solve(**kw)
        counts = _paths.since(before)
    st = tr.statistics
    state = {k: _norm(tr.__dict__.get(k, '<absent>')) for k in trmod._RUN_STATE
             if k not in _NOT_STATE_S}
    sk = tr.__dict__.get('_stage_klu')
    state['stage klu'] = None if sk is None else (sk.refactors, sk.factors, sk.pivot_fallbacks,
                                                  sk.residual_fallbacks, sk.numerics)
    state['family'] = _norm({k: v for k, v in vars(fams[-1]).items() if k != 'tr'})
    return {'family': fams[-1], 'x': np.asarray(res.x, float).tobytes(),
            't': np.asarray(res.sweep_values, float).tobytes(),
            'stats': {k: getattr(st, k) for k in st.__slots__ if 'seconds' not in k},
            'warnings': sorted(str(w.message) for w in W),
            'counts': {k: v for k, v in counts.items() if not k.startswith(_OWN)},
            'state': state}, counts


def _same(build, monkeypatch, steps_min=1, **kw):
    """The run with the C stage loop and without: the same; the first's
    paths (and what it left, `_same.last`)."""
    a, ca = _run(build, True, monkeypatch, **dict(kw))
    b, _cb = _run(build, False, monkeypatch, **dict(kw))
    for k in ('x', 't', 'stats', 'warnings', 'counts', 'state'):
        if k == 'state':
            assert [n for n in a[k] if a[k][n] != b[k][n]] == [], 'the state moved'
        elif k == 'counts':
            diff = {n: (a[k].get(n), b[k].get(n)) for n in set(a[k]) | set(b[k])
                    if a[k].get(n) != b[k].get(n)}
            assert not diff, (('the counts the C makes for the Python attempts\' calls '
                               '(`_tran_stageloop_c._count`) moved -- (C, Python):'), diff)
        else:
            assert a[k] == b[k], k
    assert ca.get('stageloop_c:steps', 0) >= steps_min, ca
    _same.last = a
    return ca


def _planted(monkeypatch, plant):
    """Each binding of the stage run passed to `plant(run)` before its use."""
    bind = _tran_stageloop_c._bind

    def planted(e):
        run = bind(e)
        if run is not None:
            plant(run)
        return run
    monkeypatch.setattr(_tran_stageloop_c, '_bind', planted)


@pytest.mark.parametrize('name, build, make, kw', [
    ('mos', _mos, None, {'tend': 60 * 2e-8, 'timestep': 2e-8}),
    ('mos-fixed', _mos, None, {'tend': 60 * 2e-8, 'timestep': 2e-8, 'fixed_timestep': True}),
    ('gp', _gp, None, {'tend': 60 * 2e-8, 'timestep': 2e-8}),
    ('psp-screen', _psp_stage, None, {'tend': 40 * 2e-8, 'timestep': 2e-8}),
    ('pulse', _pulsed, None, _PULSE),
    ('pulse-fixed', _pulsed, None, dict(_PULSE, fixed_timestep=True)),
    ('stage-klu', lambda: _mos(90, cap=1e-13), None, {'tend': 20 * 2e-10, 'timestep': 2e-10}),
    ('pointlocal', _pulsed, {'relref': 'pointlocal'}, _PULSE),
    ('alllocal', _pulsed, {'relref': 'alllocal'}, _PULSE),
    ('measured', _pulsed, _TIGHT, _PULSE),
    ('assumed', _pulsed, dict(_TIGHT, retry_order='assumed'), _PULSE),
    ('carried', _pulsed, dict(_TIGHT, retry_order='carried'), _PULSE),
    ('edge', _edge, {'reltol': 1e-6}, {'tend': 4e-8, 'timestep': 1e-9}),
    ('mask', _pulsed, {'lte_algebraic': 'skip'}, _PULSE),
    ('judge', _pulsed, {'lte_algebraic': 'judge'}, _PULSE),
    ('charge', _pulsed_load, None, _PULSE),
])
def test_a_run_is_the_python_one(name, build, make, kw, monkeypatch):
    """Every attempt -- accepted, rejected at a point and retried, a landing,
    an order drop -- the same bytes, statistics, warnings, tallies, state
    and family."""
    ca = _same(build, monkeypatch, steps_min=10, make=make, **kw)
    if build is _pulsed and not kw.get('fixed_timestep'):
        assert ca.get('stageloop_c:attempts', 0) > ca.get('stageloop_c:steps', 0), ca
    if name == 'stage-klu':
        assert ca.get('radau_c:on_klu', 0) > 10, ca
    if name in ('measured', 'carried', 'edge'):
        ## (two rejections at a point: the order measured -- at the floor
        ## behind the edge)
        q = _same.last['family']._q_carried
        assert q is not None and (q == 1.0) == (name == 'edge'), q


def test_an_entry_that_fills_its_outputs_comes_back(monkeypatch):
    """Seven rows an entry: the run in many entries, the same run."""
    monkeypatch.setattr(_tran_loop_c, 'CAP_MAX', 7)
    ca = _same(_pulsed, monkeypatch, steps_min=40, **_PULSE)
    assert ca.get('stageloop_c:steps', 0) > 7 * 5, ca


def test_a_screen_that_fires_hands_the_attempt_back(monkeypatch):
    """The PSP stage driven hard and the screen's tolerance raised on the
    instance, so it fires where a stage's C dips below the run's largest:
    the C stops after that attempt's Newton, its result held -- the
    Python's repeat takes it, its SVD the Python's -- the rest in C."""
    def build():
        c = _psp_stage()
        c['vg'] = VSin('g', gnd, v=0.7, va=0.5, freq=1e6)
        return c

    def setup(tr):
        tr.BRANCH_SCREEN_TOL = 0.9
    ca = _same(build, monkeypatch, steps_min=5, setup=setup, tend=1e-6, timestep=2e-8)
    assert ca.get('stageloop_c:bail:screen', 0) > 0, ca
    assert ca.get('radau_c:held', 0) == ca['stageloop_c:bail:screen'], ca


def _dgesv_failing(every):
    """A plant: the C estimate's solve (numpy's dgesv) made, then said
    singular at every `every`-th call."""
    sig = 'void(int64_t *, int64_t *, double *, int64_t *, int64_t *, double *, int64_t *, int64_t *)'

    def plant(run):
        ffi, s = run.ffi, run.s
        real = ffi.cast(sig.replace('void(', 'void(*)(', 1), s.dgesv)
        calls = []

        def solve(*args):
            real(*args)
            calls.append(1)
            if len(calls) % every == 0:
                args[-1][0] = 1
        cb = ffi.callback(sig, solve)
        run.skeep.append(cb)
        s.dgesv = ffi.cast('void *', cb)
    return plant


@pytest.mark.parametrize('name, build, kw', [
    ('dense', _pulsed, _PULSE),
    ('stage-klu', lambda: _mos(90, cap=1e-13), {'tend': 20 * 2e-10, 'timestep': 2e-10}),
])
def test_the_estimate_failing_in_c_hands_the_attempt_back(name, build, kw, monkeypatch):
    """The C estimate's solve said singular at every fifth call: the C
    stops after that attempt's Newton, its result -- on the stage KLU, the
    numeric it moved -- held and taken by the Python's repeat, whose own
    solve makes the estimate: the same run, the stage KLU's counts too."""
    _planted(monkeypatch, _dgesv_failing(5))
    ca = _same(build, monkeypatch, steps_min=10, **kw)
    assert ca.get('stageloop_c:bail:estimate', 0) > 2, ca
    assert ca.get('radau_c:held', 0) == ca['stageloop_c:bail:estimate'], ca
    if name == 'stage-klu':
        assert ca.get('radau_c:on_klu', 0) > 10, ca


def test_the_sources_failing_in_c_hand_the_attempt_back(monkeypatch):
    """The C's sources handed back at the 7th, 15th, 24th ... time it asks
    (gaps one longer each time: a fixed gap meets the same call of an
    attempt every time): before an attempt's Newton (a stage time) the
    attempt is the Python's, nothing of it made; after it (the step's end,
    its start) its result is held and taken -- the same run."""
    fails = set(itertools.accumulate(range(7, 400)))

    def plant(run):
        tab = run.tab
        tffi, real = tab.ffi, tab.cfn
        calls = []

        def entry(src):
            st = real(src)
            calls.append(1)
            return st + 100 if len(calls) in fails else st
        cb = tffi.callback('long(sources_t *)', entry)
        run.skeep.append(cb)
        run.s.src_fn = run.ffi.cast('void *', int(tffi.cast('uintptr_t', cb)))
    _planted(monkeypatch, plant)
    ca = _same(_pulsed, monkeypatch, steps_min=10, **_PULSE)
    bails, held = ca.get('stageloop_c:bail:sources', 0), ca.get('radau_c:held', 0)
    assert 0 < held < bails, ca


def test_the_budget_spent_is_the_pythons_force(monkeypatch):
    """No rejection at a point: each attempt the error test rejects spends
    the budget, and stops the C after its Newton, its result held -- the
    force-accept the Python's."""
    monkeypatch.setattr(_StageSteps, 'max_reject', 0)
    ca = _same(_pulsed, monkeypatch, steps_min=10, **_PULSE)
    assert ca.get('stageloop_c:bail:force', 0) > 0, ca
    assert ca.get('radau_c:held', 0) == ca['stageloop_c:bail:force'], ca


def test_a_newton_handed_back(monkeypatch):
    """The C stage Newton's assembly slots at three: the attempts that need
    more stop the C before their Newton is made (nothing of them held), the
    Python's Newton takes them, the C the attempts after."""
    monkeypatch.setattr(_tran_radau_c, 'MAXA', 3)
    ca = _same(_pulsed, monkeypatch, steps_min=10, **_PULSE)
    assert ca.get('stageloop_c:bail:radau', 0) > 0 and not ca.get('radau_c:held'), ca


def test_the_declines(monkeypatch):
    """A piece of the attempt patched, one shadowed on the transient, a
    constant of the law on the family's instance, numpy's sum not the C's,
    the loop's output set over its budget: the Python attempts, every one,
    the same run; switched off, nothing counted."""
    kw = dict(_PULSE)
    real = _StageSteps.judge
    with monkeypatch.context() as mp:
        mp.setattr(_StageSteps, 'judge', lambda self, *a, **k: real(self, *a, **k))
        ca = _same(_pulsed, monkeypatch, steps_min=0, **kw)
    assert ca.get('stageloop_c:patched', 0) > 0 and not ca.get('stageloop_c:steps'), ca

    def shadow(tr):
        tr._radau_error_estimate = lambda *a, **k: Transient._radau_error_estimate(tr, *a, **k)
    ca = _same(_pulsed, monkeypatch, steps_min=0, setup=shadow, **kw)
    assert ca.get('stageloop_c:shadow', 0) > 0 and not ca.get('stageloop_c:steps'), ca
    init = _StageSteps.__init__

    def safer(self, *a, **k):
        init(self, *a, **k)
        self.SAFETY = 0.8
    with monkeypatch.context() as mp:
        mp.setattr(_StageSteps, '__init__', safer)
        _a, ca = _run(_pulsed, True, monkeypatch, **kw)
    assert ca.get('stageloop_c:family', 0) > 0 and not ca.get('stageloop_c:steps'), ca
    with monkeypatch.context() as mp:
        mp.setitem(_tran_stageloop_c._MOD, 'sum', False)
        ca = _same(_pulsed, monkeypatch, steps_min=0, **kw)
    assert ca.get('stageloop_c:sum', 0) > 0 and not ca.get('stageloop_c:steps'), ca
    with monkeypatch.context() as mp:
        ## (the loop's own output set over the C stage Newton's budget)
        mp.setattr(_tran_radau_c, 'OUT_BYTES', 2 ** 10)
        ca = _same(_pulsed, monkeypatch, steps_min=0, **kw)
    assert ca.get('stageloop_c:size', 0) > 0 and not ca.get('stageloop_c:steps'), ca
    _a, cb = _run(_pulsed, False, monkeypatch, **kw)
    assert not [k for k in cb if k.startswith('stageloop_c:')], cb


def test_the_transform_and_pcnr_are_the_pythons(monkeypatch):
    """The cost transform on, or PCNR: no attempt in C, the same run."""
    ca = _same(_pulsed, monkeypatch, steps_min=0, make={'radau_transform': True}, **_PULSE)
    assert not ca.get('stageloop_c:steps'), ca
    ca = _same(_pulsed, monkeypatch, steps_min=0, make={'pcnr': True}, **_PULSE)
    assert not ca.get('stageloop_c:steps'), ca


@pytest.mark.parametrize('field', ['safety', 'gamma_r'])
def test_a_law_or_the_estimate_a_hair_off_is_another_run(field, monkeypatch):
    """Alive: the law's safety factor or the estimate's real eigenvalue a
    ulp off in the C -- the comparison sees another run."""
    a, _ = _run(_pulsed, False, monkeypatch, **_PULSE)
    load = _tran_stageloop_c._load

    def planted(loop, e, run_):
        load(loop, e, run_)
        setattr(run_.s, field, math.nextafter(getattr(run_.s, field), 0.0))
    monkeypatch.setattr(_tran_stageloop_c, '_load', planted)
    b, cb = _run(_pulsed, True, monkeypatch, **_PULSE)
    assert cb.get('stageloop_c:steps', 0) > 10 and (a['x'] != b['x'] or a['t'] != b['t'])


#: (the functions whose returned strings are decline reasons: the module's,
#: and the C loop's sources check, whose reasons it passes on)
_REASONERS = ('_eligible', '_run_ready', '_family_ready', '_state_ready', '_judge_ready',
              '_sources')


def _reasons(text):
    """Every reason `text` declines with: `_declined(loop, '<x>')` and the
    strings the reasoners return (`_REASONERS`)."""
    import ast
    out = set()
    for node in ast.walk(ast.parse(text)):
        if (isinstance(node, ast.Call) and getattr(node.func, 'id', None) == '_declined'
                and len(node.args) > 1 and isinstance(node.args[1], ast.Constant)):
            out.add(node.args[1].value)
        if isinstance(node, ast.FunctionDef) and node.name in _REASONERS:
            out |= {r.value.value for r in ast.walk(node) if isinstance(r, ast.Return)
                    and isinstance(r.value, ast.Constant) and type(r.value.value) is str}
    ## (a reasoner's '' is none: the checks passed)
    return out - {''}


def test_every_decline_says_whether_it_holds_for_the_run():
    """`_KEEPS` names every reason the stage loop declines with -- its own
    and the C loop's sources check's (a reason missing would raise at the
    decline) -- and no other."""
    found = (_reasons(inspect.getsource(_tran_stageloop_c))
             | _reasons(inspect.getsource(_tran_loop_c._sources)))
    assert found == set(_tran_stageloop_c._KEEPS), found ^ set(_tran_stageloop_c._KEEPS)
    planted = "def run(loop):\n    return _declined(loop, 'nowhere')\n" \
              "def _family_ready(loop):\n    return 'elsewhere'\n"
    assert _reasons(planted) == {'nowhere', 'elsewhere'}
