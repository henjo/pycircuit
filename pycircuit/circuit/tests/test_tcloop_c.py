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
declines; and alive: P or the estimate a hair off is another run.  The
real factor on the circuit's StructuralKLU: every attempt's real solves on
the numeric it found, zero pivots and residual fallbacks, a nonzero off the
pattern and a failed solve stopping the C (the numerics it made freed), the
fallback the Python's; the solver's counts and its one numeric compared.
The fold on the core's pattern (a chain past the core's fill): x_n's C and
G read in place, the factors made at its positions, dense and on KLU.
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
from pycircuit.circuit.tests.test_radau_tc import (
    _crc,
    _klu,
    _klu_chain,
    _klu_state,
    _without_gm0,
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
#: (steps across the pulse's edges three Newton iterations at most: the
#: transform's frozen Jacobian does not converge where the dense Newton's
#: does -- its fallback -- on a fixed grid and, the early abort armed, on an
#: adaptive one)
_FB_FIXED = {'tend': 1.6e-7, 'timestep': 2e-9, 'fixed_timestep': True}
_FB_ADAPT = {'tend': 1.6e-7, 'timestep': 1e-9}


def _pulsed40():
    """The pulsed chain at 40 stages (45 unknowns): sparse enough for the
    core's pattern (the fold on it), on the dense real factor."""
    return _pulsed(40)


def _maxiter(k):
    def setup(tr):
        tr.par.maxiter = k
    return setup


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
    state['fallbacks'] = tr.__dict__.get('_radau_transform_fallbacks')
    state['real solver'] = _klu_state(tr)
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
    ('pulse-40', _pulsed40, None, _PULSE),
    ('pulse-40-fixed', _pulsed40, None, dict(_PULSE, fixed_timestep=True)),
    ('pulse-40-tight', _pulsed40, {'reltol': 1e-7}, _PULSE),
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
    """Two Newton iterations at most: the run's first fallback (before it the
    dense Newton's record does not exist: the Python's makes it) and those
    whose dense Newton does not converge either stop the C -- the Python's
    repeat falls back, its Newton then its ladder -- the rest in C."""
    ca = _same(_pulsed, monkeypatch, steps_min=5, setup=_maxiter(2), **_PULSE)
    assert ca.get('tcloop_c:bail:transform', 0) > 0, ca
    assert ca.get('tcloop_c:bail:radau', 0) > 0, ca


@pytest.mark.parametrize('kw, make', [(_FB_FIXED, None), (_FB_ADAPT, {'reltol': 3e-2})],
                         ids=['fixed', 'adaptive'])
def test_a_fallback_is_made_in_c(kw, make, monkeypatch):
    """The transform's Newton not converging in three iterations where the
    dense one does: the fallback made in C -- the memo rolled again, the
    dense Newton on the stage KLU, its records, its screen, the estimate on
    the fresh real factor; adaptive, after the early abort's check -- the
    run the Python's, the fallbacks counted on the transient too."""
    ca = _same(_pulsed, monkeypatch, steps_min=10, make=make, setup=_maxiter(3), **kw)
    assert ca.get('tcloop_c:fallbacks', 0) > 1, ca
    assert _same.last['state']['fallbacks'] >= ca['tcloop_c:fallbacks']
    if not kw.get('fixed_timestep'):
        assert ca.get('radau.early:checks', 0) > 0, ca


def test_after_the_cs_fallbacks_the_record_reads_its_own_inputs(monkeypatch):
    """The C's fallbacks point the transient's dense Newton record at the
    loop's buffers while they run, and back at its own after (`_call`): the
    Python's next fallback writes its charge, sources and seeds there.  Left
    on the loop's, it read the C's last attempt's -- the twin below, its
    first transient in a fresh process, counted otherwise."""
    a, ca = _run(_pulsed, True, monkeypatch, setup=_maxiter(3), **_FB_FIXED)
    assert ca.get('tcloop_c:fallbacks', 0) > 0, ca
    rec = a['family'].tr.__dict__['_radau_c']
    for k in ('qn', 'u', 'seed'):
        assert int(rec.ffi.cast('uintptr_t', getattr(rec.s, k))) == rec.buf[k].ctypes.data, k


def test_the_fallback_twin_as_a_processs_first_transient():
    """The fixed-grid fallback twin, the C's run the first transient of a
    fresh process (no earlier test warming the watch): the C and the Python
    then split the attempts otherwise -- a fallback the Python's after the
    C's -- and the run is still the Python's."""
    import subprocess
    import sys
    code = '''
import pytest
from pycircuit.circuit.tests import test_tcloop_c as T
from pycircuit.circuit.tests.test_loop_c import _pulsed
mp = pytest.MonkeyPatch()
ca = T._same(_pulsed, mp, steps_min=10, setup=T._maxiter(3), **T._FB_FIXED)
assert ca.get('tcloop_c:fallbacks', 0) > 1, ca
print('OK', ca.get('tcloop_c:fallbacks'))
'''
    r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=600,
                       check=False)
    assert r.returncode == 0 and 'OK' in r.stdout, (r.stdout[-2000:], r.stderr[-4000:])


def test_a_stop_after_a_fallbacks_newton_holds_it(monkeypatch):
    """The dense estimate's solve (numpy's dgesv, the fallback's alone) said
    failed at every call: each fallback stops after its converged Newton,
    the result held -- the Python's repeat makes the transform's part again
    and its fallback serves the held Newton.  The same run."""
    sig = 'void(int64_t *, int64_t *, double *, int64_t *, int64_t *, double *, int64_t *, int64_t *)'

    def plant(run):
        if not run.s.fb_ok:
            return
        ffi = run.ffi
        real = ffi.cast(sig.replace('void(', 'void(*)(', 1), run.s.dgesv)

        def solve(*args):
            real(*args)
            args[7][0] = 1
        cb = ffi.callback(sig, solve)
        run.tkeep.append(cb)
        run.s.dgesv = ffi.cast('void *', cb)
    _planted(monkeypatch, plant)
    ca = _same(_pulsed, monkeypatch, steps_min=10, make={'reltol': 3e-2}, setup=_maxiter(3),
               **_FB_ADAPT)
    assert ca.get('tcloop_c:bail:estimate', 0) > 0 and ca.get('radau_c:held', 0) > 0, ca


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


def test_the_fallbacks_newton_a_hair_off_is_another_run(monkeypatch):
    """Alive: the dense Newton's tableau entry A[0][0] a ulp off in the
    fallback's record (the transform's Newton keeps its own) -- the C's
    fallbacks' stages, the waveform, another run."""
    kw = dict(setup=_maxiter(3), **_FB_FIXED)
    a, _ = _run(_pulsed, False, monkeypatch, **kw)
    load = _tran_tcloop_c._load_fb

    def planted(e, run_):
        load(e, run_)
        rs = e.fb[0].s
        rs.A[0] = math.nextafter(rs.A[0], 0.0)
    monkeypatch.setattr(_tran_tcloop_c, '_load_fb', planted)
    b, cb = _run(_pulsed, True, monkeypatch, **kw)
    assert cb.get('tcloop_c:fallbacks', 0) > 0 and a['x'] != b['x']


## -- the real factor on the circuit's structural KLU ------------------------------------

#: (the KLU chain at a tight tolerance: attempts rejected at a point and
#: abandoned at the early abort)
_TIGHT = {'reltol': 1e-7}


def _setups(*fns):
    def setup(tr):
        for f in fns:
            f(tr)
    return setup


@pytest.mark.parametrize('name, make, setup, kw', [
    ('fixed', None, _klu(), dict(_PULSE, fixed_timestep=True)),
    ('adaptive', None, _klu(), _PULSE),
    ('tight', _TIGHT, _klu(), _PULSE),
    ('residual-every', _TIGHT, _klu(tol=-1.0), _PULSE),
    ('residual-some', None, _klu(tol=1e-16), dict(_PULSE, fixed_timestep=True)),
    ('residual-some-tight', _TIGHT, _klu(tol=1e-16), _PULSE),
    ('zero-pivots', _TIGHT, _klu(fail=_crc(3)), _PULSE),
])
def test_a_run_on_the_klu_real_factor_is_the_python_one(name, make, setup, kw, monkeypatch):
    """The transform's real factor on the circuit's StructuralKLU, every
    attempt in C -- accepted, rejected at a point, abandoned at the early
    abort -- its real solves (the Newton's, the early check's, the
    estimate's) on the numeric the attempt found, a zero pivot's and a
    failed residual's fresh factors, the sparse verdict and the dense one:
    the run the Python's -- bytes, statistics, warnings, tallies, state, the
    solver's counts, one numeric left and the solver's."""
    ca = _same(_klu_chain, monkeypatch, steps_min=10, make=make, setup=setup, **kw)
    real, numerics = _same.last['state']['real solver']
    assert numerics == (1, True), numerics
    if make:
        assert ca.get('tcloop_c:attempts', 0) > ca.get('tcloop_c:steps', 0), ca
        assert ca.get('radau.early:aborts', 0) > 0, ca
    if name == 'residual-every':
        assert real[3] == real[0] > 0, real
    if name.startswith('residual-some'):
        assert 0 < real[3] < real[0] and ca.get('tcloop_c:klu_dense', 0) > 0, (real, ca)
    if name == 'zero-pivots':
        assert real[2] > 0, real


def test_klu_entries_that_fill_their_outputs_hand_the_numeric_on(monkeypatch):
    """Seven rows an entry on the KLU chain, some residuals failing: the
    numeric each entry ends on the solver's, the next entry's start -- the
    same run, one numeric left."""
    monkeypatch.setattr(_tran_loop_c, 'CAP_MAX', 7)
    ca = _same(_klu_chain, monkeypatch, steps_min=40, setup=_klu(tol=1e-16), **_PULSE)
    assert ca['tcloop_c:steps'] > 7 * 5, ca
    assert _same.last['state']['real solver'][1] == (1, True)


def test_a_nonzero_off_the_klu_pattern_stops_the_c(monkeypatch):
    """The structure given without the first MOSFET's gm: while the pulse
    holds it on, the real factor has a nonzero off the pattern (`_ours`:
    the fallback solver's) -- the C stops at that attempt's factors, the
    Python's makes it; the others in C; the same run."""
    ca = _same(_klu_chain, monkeypatch, steps_min=10, setup=_klu(pattern=_without_gm0), **_PULSE)
    assert ca.get('tcloop_c:bail:fold', 0) > 0, ca
    assert _same.last['state']['real solver'][0][4] > 0


def _klu_solve_planted(after):
    """A plant: the C's KLU solves of the real factor (the Newton's, the
    early check's, the estimate's) made, then `after(calls, n, B)` -- its
    count, the size and the solution -- returning the status the C reads;
    the Python's own solves untouched."""
    sig = 'int(void *, void *, int, int, double *, void *)'

    def plant(run):
        if not run.s.rk_klu:
            return
        tffi, k = run.tctx.ffi, run.tctx.s.rk
        real = tffi.cast(sig.replace('int(', 'int(*)(', 1), k.klu_solve)
        calls = []

        def solve(sym, num, n, nrhs, B, common):
            ok = real(sym, num, n, nrhs, B, common)
            calls.append(1)
            return after(len(calls), n * nrhs, B) and ok
        cb = tffi.callback(sig, solve)
        run.tkeep.append(cb)
        k.klu_solve = tffi.cast('void *', cb)
    return plant


def test_a_klu_solve_failing_in_c_takes_the_attempts_numerics_back(monkeypatch):
    """The C's KLU solves of the real factor said failed at every 13th call,
    each refactor's residual failing (every solve a fresh factor): the C
    stops at that attempt, the numerics it made freed and the one it found
    the solver's again -- the Python's repeat makes it; the same run, one
    numeric left."""
    _planted(monkeypatch, _klu_solve_planted(lambda k, n, B: k % 13 != 0))
    ca = _same(_klu_chain, monkeypatch, steps_min=10, make=_TIGHT, setup=_klu(tol=-1.0),
               **_PULSE)
    stops = sum(v for k, v in ca.items() if k.startswith('tcloop_c:bail:'))
    assert stops > 2, ca
    assert _same.last['state']['real solver'][1] == (1, True)


def test_where_the_real_factor_is_klu_the_fallback_is_the_pythons(monkeypatch):
    """Three Newton iterations at most on the KLU chain: a transform that
    does not converge stops the C -- the fallback's dense estimate is
    numpy's solve, the Python's the analysis solver's: not bound -- and the
    Python's repeat falls back; the same run."""
    ca = _same(_klu_chain, monkeypatch, steps_min=10, setup=_setups(_klu(), _maxiter(3)),
               **_FB_FIXED)
    assert ca.get('tcloop_c:bail:transform', 0) > 0 and not ca.get('tcloop_c:fallbacks'), ca


def test_the_klu_solutions_a_hair_off_are_another_run(monkeypatch):
    """Alive: the C's KLU solutions of the real factor each a ulp up -- the
    waveform another run's."""
    a, _ = _run(_klu_chain, False, monkeypatch, setup=_klu(), **_PULSE)

    def up(k, n, B):
        for i in range(n):
            B[i] = math.nextafter(B[i], math.inf)
        return True
    _planted(monkeypatch, _klu_solve_planted(up))
    b, cb = _run(_klu_chain, True, monkeypatch, setup=_klu(), **_PULSE)
    assert cb.get('tcloop_c:steps', 0) > 10 and a['x'] != b['x']


## -- the fold on the core's pattern -----------------------------------------------------

def _fold_modes(monkeypatch):
    """Each binding's fold: `(onpat, pn, m)` as bound."""
    seen = []

    def plant(run):
        seen.append((int(run.fo.onpat), int(run.fo.pn), int(run.fo.m)))
    _planted(monkeypatch, plant)
    return seen


@pytest.mark.parametrize('build, setup, onpat', [(_pulsed, None, 0), (_pulsed40, None, 1),
                                                 (_klu_chain, _klu(), 1)],
                         ids=['small', 'dense', 'klu'])
def test_the_fold_is_on_the_cores_pattern_past_its_fill(build, setup, onpat, monkeypatch):
    """The core writes only its pattern past its fill (`_tran_core.
    PATTERN_MAX_FILL`): the loop's fold then makes the factors at its
    positions alone -- fewer than m x m -- on the dense real factor and on
    KLU; the small chain's core writes every entry and the fold is the whole
    one.  The run the Python's either way."""
    seen = _fold_modes(monkeypatch)
    _same(build, monkeypatch, steps_min=10, setup=setup, **_PULSE)
    assert seen and all(o == onpat for o, _pn, _m in seen), seen
    if onpat:
        assert all(0 < pn < m * m for _o, pn, m in seen), seen


def test_a_pattern_position_left_out_is_another_run(monkeypatch):
    """Alive: the fold on the pattern given one position fewer (the middle
    one, a stage's node: its factors' entry left zero; the last, a source
    branch's diagonal, is zero anyway) -- another run."""
    a, _ = _run(_pulsed40, False, monkeypatch, **_PULSE)

    def plant(run):
        if run.fo.onpat:
            k = int(run.fo.pn) // 2
            for name in ('pfull', 'pred', 'ptr'):
                at = run.tbuf['fold_' + name]
                at[k:-1] = at[k + 1:].copy()
            run.fo.pn -= 1
    _planted(monkeypatch, plant)
    b, cb = _run(_pulsed40, True, monkeypatch, **_PULSE)
    assert cb.get('tcloop_c:steps', 0) > 10 and (a['x'] != b['x'] or a['t'] != b['t'])


## -- the passes on the core's pattern -----------------------------------------------------

#: (the loop's buffers a pass writes: x_n's C and G, the screen's C, the end
#: passes' C and G in their two sets)
_PASS_BUFS = ('xr_C', 'xr_G', 'tcC', 'te02', 'te03', 'te12', 'te13')


@pytest.mark.parametrize('build, setup', [(_pulsed40, None), (_klu_chain, _klu())],
                         ids=['dense', 'klu'])
def test_the_pass_buffers_stay_zero_off_the_pattern(build, setup, monkeypatch):
    """The core asked for its pattern writes C there alone and G is copied
    from its bin there alone: the loop's buffers are never zeroed whole, and
    after a run they are still zero off the pattern -- as every C and G of
    the circuit is, which the passes rely on."""
    runs = []
    _planted(monkeypatch, runs.append)
    _same(build, monkeypatch, steps_min=10, setup=setup, **_PULSE)
    assert runs and all(r.s.pass_zero for r in runs)
    for r in runs:
        n = r.tctx.n
        off = np.ones(n * n, dtype=bool)
        off[np.asarray(r.tctx.core.arrays['pat'])] = False
        for k in _PASS_BUFS:
            a = r.tbuf[k].reshape(-1, n * n)
            assert not a[:, off].any(), k


def test_a_nonzero_off_the_pattern_is_another_run(monkeypatch):
    """Alive: a nonzero off the core's pattern planted in the end passes'
    two C buffers, in the pulsed input's column (its F1 moves at each edge;
    a still node's multiplies it by zero), small enough to leave the steps
    as they were (1e-15 collapses them; 1e-19 is not seen) -- nothing
    zeroes them now, and the estimate reads x_n's C whole -- another run."""
    a, _ = _run(_pulsed40, False, monkeypatch, **_PULSE)
    cir = Transient(_pulsed40(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator()).cir
    col = cir.get_node_index(cir.get_node('g0'))

    def plant(run):
        n = run.tctx.n
        off = np.ones((n, n), dtype=bool)
        off.flat[np.asarray(run.tctx.core.arrays['pat'])] = False
        off[int(run.s.iref), :] = False
        row = int(np.flatnonzero(off[:, col])[0])
        for name in ('te02', 'te12'):
            run.tbuf[name][row * n + col] = 1e-17
    _planted(monkeypatch, plant)
    b, cb = _run(_pulsed40, True, monkeypatch, **_PULSE)
    assert cb.get('tcloop_c:steps', 0) > 10 and (a['x'] != b['x'] or a['t'] != b['t'])


#: (the functions whose returned strings are decline reasons: the module's,
#: the stage loop's it calls, and the C loop's sources check)
_REASONERS = ('_eligible', '_family_ready', '_tc_ready', '_solver_ready', '_real_solver',
              '_run_ready', '_state_ready', '_judge_ready', '_sources')


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
