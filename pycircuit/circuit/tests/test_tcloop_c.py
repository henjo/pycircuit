"""The plain run's Radau steps on the cost transform in C (`_tran_tcloop_c`,
csrc/loop.c's stage family in its transform mode): a Radau IIA(3)
transient the transform makes ('auto' on a compiled Jacobian of 10 KB or
more) with the C loop on is the run with it off -- every output byte,
time, statistic and warning, the `_paths` tallies the Python attempts
make, every attribute the attempts leave on the transient
(`transient._RUN_STATE`: the device memo's two generations among them),
`_radau_Pc`, the family's own (its running reference, retry record,
carried order, `early_f0`) and the complex solver's counts -- fixed and
adaptive, the early abort armed (the default) or not, a new pattern of the
complex factor made in C (its analysis and fresh factor, the Python loop's
counts); an entry that fills its outputs; its stops (a screen that fires,
the estimate's solve or the sources failing in the C, a transform that
does not converge, the budget spent; a stop after a new pattern); its
declines; and alive: P or the estimate a hair off is another run.
"""
import gc
import inspect
import itertools
import math
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_loop_c, _tran_radau, _tran_tcloop_c, circuit
from pycircuit.circuit import transient as trmod
from pycircuit.circuit.elements import VSin, gnd
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_loop_c import (
    _NOT_STATE,
    _OWN_COUNTS,
    _mos,
    _norm,
    _psp_stage,
    _pulsed,
)
from pycircuit.circuit.transient import Transient, _StageSteps

pytestmark = pytest.mark.skipif(_tran_loop_c.driver() is None,
                                reason='the stepping loop unit did not build')

## (workspaces: the predictor's for a stage step's seeds, the transform
## record's buffers -- every call fills what it reads)
_NOT_STATE_T = _NOT_STATE | {'_pred_c', '_radau_tc'}
## (counts the loops' own machinery makes: their declines and bails, the
## coupled mode's holds)
_OWN = _OWN_COUNTS + ('stageloop_c:', 'tcloop_c:', 'radau_c:held')
#: (a pulse's corners: rejections at a point, landings, order drops, the
#: devices cut off -- the complex factor's pattern moving)
_PULSE = {'tend': 4e-8, 'timestep': 2.5e-10}
#: (the cost transform asked for: 'auto' leaves the MOSFET chains, their
#: Jacobians bound to C kernels, on the coupled Newton)
_TRANSFORM = {'radau_transform': True}


def _run(build, on, monkeypatch, make=None, setup=None, **kw):
    """One run with the C transform loop on or off: `(what it left, the
    paths counted)`.  (The earlier runs' objects freed first: a watched dict
    freed during a run moves the watch epoch, which the two runs would meet
    at different points.)"""
    fams = []
    init = _StageSteps.__init__

    def kept(self, *a, **k):
        fams.append(self)
        init(self, *a, **k)
    gc.collect()
    with monkeypatch.context() as mp:
        mp.setattr(_tran_tcloop_c, 'ENABLED', on)
        mp.setattr(_StageSteps, '__init__', kept)
        before = _paths.snapshot()
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                           **dict(_TRANSFORM, **(make or {})))
            if setup:
                setup(tr)
            res = tr.solve(**kw)
        counts = _paths.since(before)
    st = tr.statistics
    state = {k: _norm(tr.__dict__.get(k, '<absent>')) for k in trmod._RUN_STATE
             if k not in _NOT_STATE_T}
    zs = tr.__dict__.get('_radau_zsolver')
    state['complex solver'] = None if zs is None else (
        zs.refactors, zs.factors, zs.analyses, zs.residual_fallbacks, zs.analyses_freed,
        zs.factors_freed)
    pc = tr.__dict__.get('_radau_Pc')
    state['_radau_Pc'] = None if pc is None else (
        repr(pc[0]), np.array(pc[1], dtype=complex).tobytes())
    state['family'] = _norm({k: v for k, v in vars(fams[-1]).items() if k != 'tr'})
    return {'family': fams[-1], 'x': np.asarray(res.x, float).tobytes(),
            't': np.asarray(res.sweep_values, float).tobytes(),
            'stats': {k: getattr(st, k) for k in st.__slots__ if 'seconds' not in k},
            'warnings': sorted(str(w.message) for w in W),
            'counts': {k: v for k, v in counts.items() if not k.startswith(_OWN)},
            'state': state}, counts


def _same(build, monkeypatch, steps_min=1, **kw):
    """The run with the C transform loop and without: the same; the first's
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
                               '(`_tran_tcloop_c._count`) moved -- (C, Python):'), diff)
        else:
            assert a[k] == b[k], k
    assert ca.get('tcloop_c:steps', 0) >= steps_min, ca
    _same.last = a
    return ca


def _planted(monkeypatch, plant):
    """Each binding of the transform run passed to `plant(run)` before its
    use."""
    bind = _tran_tcloop_c._bind

    def planted(e):
        run = bind(e)
        if run is not None:
            plant(run)
        return run
    monkeypatch.setattr(_tran_tcloop_c, '_bind', planted)


def _hard_psp():
    """The PSP stage driven hard: its C moves, the screen can fire."""
    c = _psp_stage()
    c['vg'] = VSin('g', gnd, v=0.7, va=0.5, freq=1e6)
    return c


@pytest.mark.parametrize('name, build, make, kw', [
    ('psp-fixed', _psp_stage, None, {'tend': 40 * 2e-8, 'timestep': 2e-8, 'fixed_timestep': True}),
    ('psp', _psp_stage, None, {'tend': 40 * 2e-8, 'timestep': 2e-8}),
    ('psp-no-early', _psp_stage, {'radau_early_abort': 0}, {'tend': 40 * 2e-8, 'timestep': 2e-8}),
    ('mos', _mos, None, {'tend': 60 * 2e-8, 'timestep': 2e-8}),
    ('pulse', _pulsed, None, _PULSE),
    ('pulse-fixed', _pulsed, None, dict(_PULSE, fixed_timestep=True)),
    ('pulse-no-early', _pulsed, {'radau_early_abort': 0}, _PULSE),
    ('pulse-early-one', _pulsed, {'radau_early_abort': 1}, _PULSE),
    ('pulse-tight', _pulsed, {'reltol': 1e-7}, _PULSE),
    ('pulse-pointlocal', _pulsed, {'relref': 'pointlocal'}, _PULSE),
    ('pulse-skip', _pulsed, {'lte_algebraic': 'skip'}, _PULSE),
])
def test_a_run_is_the_python_one(name, build, make, kw, monkeypatch):
    """Every attempt -- accepted, rejected at a point and retried, abandoned
    at the early abort, on a new pattern of the complex factor, a landing,
    an order drop -- the same bytes, statistics, warnings, tallies, state,
    family and complex solver."""
    ca = _same(build, monkeypatch, steps_min=10, make=make, **kw)
    if build is _pulsed and not kw.get('fixed_timestep'):
        assert ca.get('tcloop_c:attempts', 0) > ca.get('tcloop_c:steps', 0), ca


def test_the_early_abort_and_a_new_pattern_are_made_in_c(monkeypatch):
    """The pulsed chain at its defaults: the early abort armed and firing,
    the complex factor's pattern moving -- both made in C (the checks and
    aborts counted as the Python's; the new patterns' analyses and fresh
    factors on the complex solver, the ones they replaced freed), the run
    the Python's."""
    ca = _same(_pulsed, monkeypatch, steps_min=10, **_PULSE)
    assert ca.get('radau.early:checks', 0) > 0 and ca.get('radau.early:aborts', 0) > 0, ca
    assert ca.get('radau_tc:klu', 0) > 0, ca
    assert not ca.get('tcloop_c:bail:fold'), ca
    ## (each pattern's analysis and factor freed where the next replaced it,
    ## by the C as by the Python: one of each left)
    _r, factors, analyses, _f, a_freed, f_freed = _same.last['state']['complex solver']
    assert (analyses - a_freed, factors - f_freed) == (1, 1), _same.last['state']['complex solver']


def test_an_entry_that_fills_its_outputs_comes_back(monkeypatch):
    """Seven rows an entry: the run in many entries, the same run."""
    monkeypatch.setattr(_tran_loop_c, 'CAP_MAX', 7)
    ca = _same(_pulsed, monkeypatch, steps_min=40, **_PULSE)
    assert ca.get('tcloop_c:steps', 0) > 7 * 5, ca


def test_a_screen_that_fires_hands_the_attempt_back(monkeypatch):
    """The PSP stage driven hard and the screen's tolerance raised on the
    instance, so it fires where a stage's C dips below the run's largest:
    the C stops after that attempt's Newton (nothing held: the Python's
    repeat refactors the same values), its SVD the Python's, the rest in
    C."""
    def setup(tr):
        tr.BRANCH_SCREEN_TOL = 0.9
    ca = _same(_hard_psp, monkeypatch, steps_min=5, setup=setup, tend=1e-6, timestep=2e-8)
    assert ca.get('tcloop_c:bail:screen', 0) > 0, ca


def test_a_transform_that_does_not_converge_is_the_pythons(monkeypatch):
    """Two Newton iterations at most: the attempts the transform does not
    converge in stop the C -- the Python's repeat falls back to the dense
    Newton -- the rest in C."""
    def setup(tr):
        tr.par.maxiter = 2
    ca = _same(_pulsed, monkeypatch, steps_min=5, setup=setup, **_PULSE)
    assert ca.get('tcloop_c:bail:transform', 0) > 0, ca


def _getrs_failing(every):
    """A plant: the C's real solves after the factor (numpy's dgetrs: the
    early check's and the estimate's) made, then said failed at every
    `every`-th call."""
    sig = ('void(char *, int64_t *, int64_t *, double *, int64_t *, int64_t *, double *, '
           'int64_t *, int64_t *, size_t)')

    def plant(run):
        tffi, ts = run.tctx.ffi, run.tctx.s
        real = tffi.cast(sig.replace('void(', 'void(*)(', 1), ts.dgetrs)
        calls = []

        def solve(*args):
            real(*args)
            calls.append(1)
            if len(calls) % every == 0:
                args[8][0] = 1
        cb = tffi.callback(sig, solve)
        run.tkeep.append(cb)
        ts.dgetrs = tffi.cast('void *', cb)
    return plant


def test_the_estimate_failing_in_c_hands_the_attempt_back(monkeypatch):
    """The C's real solves after the factor (the Newton's later ones, the
    early check's, the estimate's) said failed at every 23rd call: the C
    stops at that attempt, the Python's repeat makes it -- the same run."""
    _planted(monkeypatch, _getrs_failing(23))
    ca = _same(_pulsed, monkeypatch, steps_min=10, **_PULSE)
    stops = sum(v for k, v in ca.items() if k.startswith('tcloop_c:bail:'))
    assert stops > 2, ca


def test_the_sources_failing_in_c_hand_the_attempt_back(monkeypatch):
    """The C's sources handed back at the 7th, 15th, 24th ... time it asks
    (gaps one longer each time): before an attempt's Newton the attempt is
    the Python's; after it too (nothing held) -- the same run."""
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
        run.tkeep.append(cb)
        run.s.src_fn = run.ffi.cast('void *', int(tffi.cast('uintptr_t', cb)))
    _planted(monkeypatch, plant)
    ca = _same(_pulsed, monkeypatch, steps_min=10, **_PULSE)
    assert ca.get('tcloop_c:bail:sources', 0) > 2, ca


def test_a_stop_after_a_new_pattern_takes_it_back(monkeypatch):
    """An attempt that made a new pattern of the complex factor, stopped at
    its next source (a stage time): its analysis and fresh factor freed,
    the pattern before put back -- the Python's repeat analyses and counts
    them itself."""
    def plant(run):
        tab = run.tab
        tffi, real = tab.ffi, tab.cfn
        s = run.s
        ## (the loop's pattern buffer moved since the last source call: the
        ## attempt made a new pattern, `tc_repattern`, before its sources)
        seen = {'pcur': -1}

        def entry(src):
            st = real(src)
            if s.pcur >= 0 and s.pcur != seen['pcur']:
                seen['pcur'] = s.pcur
                return st + 100
            return st
        cb = tffi.callback('long(sources_t *)', entry)
        run.tkeep.append(cb)
        run.s.src_fn = run.ffi.cast('void *', int(tffi.cast('uintptr_t', cb)))
    _planted(monkeypatch, plant)
    ca = _same(_pulsed, monkeypatch, steps_min=10, **_PULSE)
    assert ca.get('tcloop_c:bail:sources', 0) > 0, ca


def test_the_budget_spent_is_the_pythons_force(monkeypatch):
    """No rejection at a point: each attempt the error test rejects spends
    the budget and stops the C -- the force-accept the Python's."""
    monkeypatch.setattr(_StageSteps, 'max_reject', 0)
    ca = _same(_pulsed, monkeypatch, steps_min=10, **_PULSE)
    assert ca.get('tcloop_c:bail:force', 0) > 0, ca


def test_the_declines(monkeypatch):
    """A piece of the transform's attempt patched, one shadowed on the
    transient, its frozen factors switched off, the cost transform off (the
    coupled mode's): the Python attempts, every one, the same run; switched
    off, nothing counted."""
    kw = dict(_PULSE)
    real = _tran_radau._RadauStages._radau_frozen
    with monkeypatch.context() as mp:
        mp.setattr(_tran_radau._RadauStages, '_radau_frozen',
                   lambda self, *a, **k: real(self, *a, **k))
        ca = _same(_pulsed, monkeypatch, steps_min=0, **kw)
    assert not ca.get('tcloop_c:steps'), ca

    def shadow(tr):
        tr._radau_early_check = lambda *a, **k: Transient._radau_early_check(tr, *a, **k)
    ca = _same(_pulsed, monkeypatch, steps_min=0, setup=shadow, **kw)
    assert ca.get('tcloop_c:shadow', 0) > 0 and not ca.get('tcloop_c:steps'), ca
    with monkeypatch.context() as mp:
        mp.setattr(_tran_radau, 'STAGE_FUSE', False)
        ca = _same(_pulsed, monkeypatch, steps_min=0, **kw)
    assert ca.get('tcloop_c:switch', 0) > 0 and not ca.get('tcloop_c:steps'), ca
    ca = _same(_pulsed, monkeypatch, steps_min=0, make={'radau_transform': False}, **kw)
    assert not ca.get('tcloop_c:steps'), ca
    _a, cb = _run(_pulsed, False, monkeypatch, **kw)
    assert not [k for k in cb if k.startswith('tcloop_c:')], cb


@pytest.mark.parametrize('field', ['tcm', 'lam0r'])
def test_p_or_the_factors_a_hair_off_is_another_run(field, monkeypatch):
    """Alive: P's first entry or the real factor's eigenvalue a ulp off in
    the C -- the comparison sees another run."""
    a, _ = _run(_pulsed, False, monkeypatch, **_PULSE)
    load = _tran_tcloop_c._load_tc

    def planted(e, run_, xl):
        load(e, run_, xl)
        s = run_.s
        if field == 'tcm':
            s.tcm[0] = math.nextafter(s.tcm[0], 0.0)
        else:
            s.lam0r = math.nextafter(s.lam0r, 0.0)
    monkeypatch.setattr(_tran_tcloop_c, '_load_tc', planted)
    b, cb = _run(_pulsed, True, monkeypatch, **_PULSE)
    assert cb.get('tcloop_c:steps', 0) > 10 and (a['x'] != b['x'] or a['t'] != b['t'])


#: (the functions whose returned strings are decline reasons: the module's,
#: the stage loop's it calls, and the C loop's sources check)
_REASONERS = ('_eligible', '_family_ready', '_tc_ready', '_solver_ready', '_run_ready',
              '_state_ready', '_judge_ready', '_sources')


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
    return out - {''}


def test_every_decline_says_whether_it_holds_for_the_run():
    """`_KEEPS` names every reason the transform loop declines with -- its
    own, the stage loop's checks it calls and the C loop's sources check's
    (a reason missing would raise at the decline) -- and no other."""
    from pycircuit.circuit import _tran_stageloop_c
    found = (_reasons(inspect.getsource(_tran_tcloop_c))
             | _reasons(inspect.getsource(_tran_stageloop_c._run_ready))
             | _reasons(inspect.getsource(_tran_stageloop_c._state_ready))
             | _reasons(inspect.getsource(_tran_stageloop_c._judge_ready))
             | _reasons(inspect.getsource(_tran_loop_c._sources)))
    assert found == set(_tran_tcloop_c._KEEPS), found ^ set(_tran_tcloop_c._KEEPS)
