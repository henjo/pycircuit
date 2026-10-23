"""The shooting's plain multistep walk in C (`shooting/_lmmwalk_c`, the C
stepping loop's walk): a PSS's walks with it and without it are one
computation -- the orbit, the monodromy, the inner transient's and the walk
object's state the walks leave (the step lengths' and coefficients' types
with their values), every count but the walk's own -- under Gear-2, the
trapezoid and Euler, dense and on the run's KLU, whatever stops it: a Python
step between its calls, chunks of a step or two, a Newton handed back, a
step map handed back, a prune (the Python's views), a source that raises."""
import hashlib
import sys

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_loop_c, _tran_newton_c, circuit
from pycircuit.circuit.shooting import PSS, _lmmsens_c, _lmmwalk_c, _pss_walks, _steps
from pycircuit.circuit.simwarnings import AccuracyWarning, ConvergenceWarning
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.circuit.tests.test_newton_c import _mos_chain

_this = sys.modules[__name__]

#: the chain's drive period (1 MHz), its grid
T, STEPS = 1e-6, 40


def _b(a):
    return np.asarray(a, dtype=float).tobytes()


def _run(monkeypatch, method, on, patch=None, n=4):
    """One PSS of an `n`-stage MOSFET chain under `method`, the C walk on or
    off, `patch(mp)` applied first: `(what it left, the paths counted)` --
    an error the PSS raised among what it left."""
    with monkeypatch.context() as mp:
        mp.setattr(circuit, 'default_toolkit', circuit.numeric)
        mp.setattr(_lmmwalk_c, 'LMMWALK_C', on)
        if patch is not None:
            patch(mp)
        before = _paths.snapshot()
        pss = PSS(_mos_chain(n=n, cap=1e-13), method=method, reltol=1e-10)
        out = {}
        try:
            ## (the coarse grid's own accuracy warning; the plain map's
            ## twenty iterations under the trapezoid and Euler, which do not
            ## converge on the chain -- the same twenty either way)
            with quiet(AccuracyWarning, ConvergenceWarning):
                res = pss.solve(period=T, timestep=T / STEPS)
            out['orbit'] = _b(res['tpss'].x)
            out['monodromy'] = _b(pss._monodromy)
        except ValueError as e:
            out['raised'] = str(e)
        counts = _paths.since(before)
    tr = pss._transient()
    for k in ('_qlast', '_iqlast', '_iq', '_Cmat', '_Geq'):
        out[k] = _b(getattr(tr, k))
    out['_q_cache'] = _b(tr._q_cache[0]) + _b(tr._q_cache[1])
    out['_pred_hist'] = [(t, _b(x)) for t, x in tr._pred_hist]
    for k in ('_companion_coeffs', '_dt', '_dt_last', '_dt_last2', '_is_first_step',
              '_branch_cmax', '_effective_method'):
        out[k] = repr(getattr(tr, k, None))
    for k in ('_Jf', '_Geq', '_C'):
        out['walk' + k] = _b(getattr(pss, k))
    out['walk_coeffs'] = repr(pss._coeffs)
    klu = tr.__dict__['_newton_c'].klu
    out['klu'] = klu and (klu.refactors, klu.numerics)
    return out, counts


def _strip(counts):
    return {k: v for k, v in counts.items() if not k.startswith('pss.lmmwalk:')}


def _twin(monkeypatch, method, patch=None, n=4):
    """The PSS without the C walk and with it: what each left and every
    count but the walk's own equal; the second's paths."""
    a, ca = _run(monkeypatch, method, False, patch, n)
    b, cb = _run(monkeypatch, method, True, patch, n)
    assert [k for k in a if a[k] != b[k]] == []
    assert _strip(ca) == _strip(cb)
    return b, cb


@pytest.fixture(autouse=True)
def _warm(monkeypatch):
    """A first PSS made for every compile's once-a-process count."""
    _run(monkeypatch, 'gear', True)


@pytest.mark.parametrize('method,n', [('gear', 4), ('gear', 40), ('trap', 4), ('euler', 4)])
def test_the_c_walk_is_the_python_walk(monkeypatch, method, n):
    """A walk's steps after its first served one are the C's -- the same
    bits, states and counts -- a prune that drops a column the Python's."""
    left, counts = _twin(monkeypatch, method, n=n)
    assert counts.get('pss.lmmwalk:steps', 0) > 2 * STEPS, counts
    assert counts.get('pss.lmmwalk:walk_prune', 0) > 0, counts
    assert (left['klu'] is not None) == (n == 40)


@pytest.mark.parametrize('how', ['interleaved', 'chunks'])
def test_python_steps_and_chunks_between_the_calls(monkeypatch, how):
    """A call cut short every third step (its grid), the Python step
    between, the C taking over again: the state crosses both ways; and
    chunks of one or two steps, the state carried in the struct."""
    def patch(mp):
        if how == 'interleaved':
            orig = _lmmwalk_c.walk
            mp.setattr(_lmmwalk_c, 'walk',
                       lambda walks, j0, x, times, *a: orig(walks, j0, x, times[:j0 + 4], *a))
        else:
            mp.setattr(_lmmwalk_c, 'CHUNK_FIRST', 1)
            mp.setattr(_lmmwalk_c, 'CHUNK_MAX', 2)
    _left, counts = _twin(monkeypatch, 'gear', patch, n=40)
    assert counts.get('pss.lmmwalk:steps', 0) > STEPS, counts


def test_a_newton_handed_back_mid_walk(monkeypatch):
    """The C Newton handing back every step whose seed hashes to 0 mod 5:
    the walk stops before it (nothing of it made), the Python step takes
    it, the C the steps after."""
    def patch(mp):
        ffi, cfn, lapack = _tran_newton_c.driver()

        def entry(t):
            seed = bytes(ffi.buffer(t.xseed, 8 * t.n))
            if hashlib.sha1(seed).digest()[0] % 5 == 0:
                return _tran_newton_c._CODES['NEWTON_BAIL_X']
            return cfn(t)
        cb = ffi.callback('long(newton_t *)', entry)
        mp.setattr(_tran_newton_c, '_driver', (ffi, cb, lapack))
    _left, counts = _twin(monkeypatch, 'gear', patch, n=40)
    assert counts.get('pss.lmmwalk:newton', 0) > 3, counts
    assert counts.get('pss.lmmwalk:steps', 0) > STEPS, counts


@pytest.mark.parametrize('bail', ['LMMSENS_BAIL_SINGULAR', 'LMMSENS_BAIL_NONFINITE'])
def test_a_step_map_handed_back(monkeypatch, bail):
    """The step map handing back every step whose Jacobian hashes to 0 mod
    7: the step made, its map made as `_walk_lmm` makes it (`_lmmsens_c`'s
    decline counted, `_LiveMonodromy.step`), the C walking on after."""
    def patch(mp):
        ffi, cfn, gesv, gemm = _lmmsens_c.driver()

        def entry(t):
            raw = bytes(ffi.buffer(t.Jf, 8 * t.m * t.m))
            if hashlib.sha1(raw).digest()[0] % 7 == 0:
                return _lmmsens_c._CODES[bail]
            return cfn(t)
        cb = ffi.callback('long(lmmsens_t *)', entry)
        mp.setattr(_lmmsens_c, '_driver', (ffi, cb, gesv, gemm))
    _left, counts = _twin(monkeypatch, 'gear', patch)
    key = 'pss.lmm:' + _lmmsens_c.BAIL[_lmmsens_c._CODES[bail]]
    assert counts.get(key, 0) > 0 and counts.get('pss.lmmwalk:walk_map', 0) > 0, counts


def test_the_walks_records_are_the_pythons(monkeypatch):
    """Every plain walk made with a `record`: the C's steps recorded as the
    Python records them -- each point and time, in order."""
    seen = {}

    def patch(mp):
        orig = _pss_walks._PeriodWalks._walk_lmm
        log = seen[_lmmwalk_c.LMMWALK_C] = []

        def walk(self, *a, **k):
            if not k.get('keep') and k.get('record') is None:
                k['record'] = lambda x, t: log.append((_b(x), repr(t)))
            return orig(self, *a, **k)
        mp.setattr(_pss_walks._PeriodWalks, '_walk_lmm', walk)
    _left, counts = _twin(monkeypatch, 'gear', patch)
    assert counts.get('pss.lmmwalk:steps', 0) > 0 and seen[True] == seen[False]
    assert len(seen[True]) > 2 * STEPS


def test_a_source_that_raises_is_the_pythons_step(monkeypatch):
    """A drive whose own `f` raises past half the period: the C walk stops
    before the step whose source raises (its table's rows from there taken
    back), the Python step raises it -- the same error, state and counts
    behind it."""
    def patch(mp):
        from pycircuit.circuit import func
        orig, chain = func.Sin.f, _this._mos_chain

        def build(n=4, cap=1e-13):
            c = chain(n=n, cap=cap)
            fn = c['vg'].function

            def f(t):
                if t > T / 2:
                    raise ValueError('the drive is gone')
                return orig(fn, t)
            fn.f = f
            return c
        mp.setattr(_this, '_mos_chain', build)
    left, counts = _twin(monkeypatch, 'gear', patch)
    assert left['raised'] == 'the drive is gone' and counts.get('pss.lmmwalk:source', 0) == 1


def test_what_declines_leaves_the_python_walk(monkeypatch):
    """Switched off, a piece it stands in for patched, a walk recording its
    steps' `C` and `J` (`_record_cj`): no step in C, the same walk."""
    _left, counts = _run(monkeypatch, 'gear', False)
    assert counts.get('pss.lmmwalk:off', 0) > 0 and not counts.get('pss.lmmwalk:steps')
    orig = _steps._LiveMonodromy._prune

    def patch(mp):
        mp.setattr(_steps._LiveMonodromy, '_prune', lambda self: orig(self))
    _left, counts = _twin(monkeypatch, 'gear', patch)
    assert counts.get('pss.lmmwalk:patched', 0) > 0 and not counts.get('pss.lmmwalk:steps')

    def cj(mp):
        mp.setattr(_pss_walks._PeriodWalks, '_record_cj', True)
    _left, counts = _twin(monkeypatch, 'gear', cj)
    assert counts.get('pss.lmmwalk:walk', 0) > 0 and not counts.get('pss.lmmwalk:steps')


def test_a_one_ulp_change_in_the_c_walk_shows(monkeypatch):
    """The C's points are the walk's: one ulp planted in each call's last
    point (the walk's end, where the call ends it), and what the walks
    leave is another."""
    a, _ca = _run(monkeypatch, 'gear', True)

    def patch(mp):
        ffi, cfn = _tran_loop_c.driver()

        def planted(s):
            st = cfn(s)
            if s.walk and s.k:
                r = (s.k - 1) * s.n + (s.iref + 1) % s.n
                s.X[r] = float(np.nextafter(s.X[r], np.inf))
            return st
        mp.setattr(_tran_loop_c, '_driver', (ffi, planted))
    b, cb = _run(monkeypatch, 'gear', True, patch)
    assert cb.get('pss.lmmwalk:steps', 0) > 0 and [k for k in a if a[k] != b[k]]
