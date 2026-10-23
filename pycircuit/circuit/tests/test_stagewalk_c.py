"""The shooting's plain stage walk in C (`shooting/_stagewalk_c`): a PSS's
walks with it and without it are one computation -- the orbit, the
monodromy, the replay's records, the inner transient's state the walks
leave, every numeric's counts -- whatever stops it: a Python step between
its calls, a Newton handed back, a stop after a converged Newton (its result
held for the Python's repeat of the step), a step map not finite, a
singular one."""
import numpy as np
import pytest

from pycircuit.circuit import _paths, circuit, transient
from pycircuit.circuit import _tran_radau_c as RC
from pycircuit.circuit import linearsolver as LS
from pycircuit.circuit.shooting import PSS, _sensklu_c, _stagewalk_c
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.circuit.tests.test_newton_c import _mos_chain
from pycircuit.circuit.tests.test_shooting_step_klu import SECTIONS, _ladder

#: (the circuit, its period, the analysis' options): a linear ladder (a
#: constant `C`), a MOSFET chain (a cheap compiled Jacobian: the cost
#: transform off, as 'auto' leaves it for such circuits)
CASES = {
    'ladder': (lambda: _ladder(SECTIONS, diode=False), 1e-3, {}),
    'mos': (lambda: _mos_chain(n=SECTIONS, cap=1e-13), 1e-6, {'radau_transform': False}),
}


def _b(a):
    return np.asarray(a, dtype=float).tobytes()


def _run(monkeypatch, case, on, patch=None):
    """One PSS of `case` with the C walk on or off, `patch(monkeypatch)`
    applied first: `(what it left, the paths counted, the error it raised
    or None)`."""
    with monkeypatch.context() as mp:
        mp.setattr(circuit, 'default_toolkit', circuit.numeric)
        mp.setattr(_stagewalk_c, 'STAGEWALK_C', on)
        if patch is not None:
            patch(mp)
        walks = []
        orig = LS.StructuralKLU.walk

        def counted(self):
            walks.append(orig(self))
            return walks[-1]
        mp.setattr(LS.StructuralKLU, 'walk', counted)
        build, period, kw = CASES[case]
        before = _paths.snapshot()
        pss = PSS(build(), method='radau', reltol=1e-10, **kw)
        out, err = {}, None
        try:
            with quiet(AccuracyWarning):
                res = pss.solve(period=period, timestep=period / 20)
            out['orbit'] = _b(res['tpss'].x)
            out['monodromy'] = _b(pss._monodromy)
            lw = pss._last_stage_walk
            out['replay'] = lw and (b''.join(_b(x) for x in lw[4]),
                                    b''.join(_b(e[0]) + _b([e[1]]) + b''.join(_b(y) for y in e[2])
                                             for e in lw[5]))
        except np.linalg.LinAlgError as e:
            err = str(e)
        counts = _paths.since(before)
    tr = pss._transient()
    for k in ('_qlast', '_iqlast', '_iq', '_Cmat', '_Geq'):
        out[k] = _b(getattr(tr, k))
    out['_rk_Y'] = b''.join(_b(y) for y in tr._rk_Y)
    out['_q_cache'] = _b(tr._q_cache[0]) + _b(tr._q_cache[1])
    out['_pred_hist'] = [(t, _b(x)) for t, x in tr._pred_hist]
    for k in ('_branch_cmax', '_dt', '_dt_last', '_dt_last2', '_is_first_step'):
        out[k] = getattr(tr, k, None)
    out['walks'] = [(w.factors, w.refactors, w.fallbacks, w.gen) for w in walks]
    sk = tr.__dict__['_stage_klu']
    out['stage_klu'] = (sk.refactors, sk.factors, sk.pivot_fallbacks, sk.residual_fallbacks,
                        sk.numerics)
    return out, counts, err


def _twin(monkeypatch, case, patch=None):
    """The PSS without the C walk and with it: what each left (equal), and
    the second's paths."""
    a, _ca, ea = _run(monkeypatch, case, False, patch)
    b, cb, eb = _run(monkeypatch, case, True, patch)
    assert ea == eb
    assert [k for k in a if a[k] != b[k]] == []
    return b, cb, eb


def _sens_entry(mp, fn):
    """The step map's C entry replaced by `fn(t, cfn)` (a cffi callback):
    the Python's `_sensklu_c.block_step` and the C walk call the same."""
    ffi, cfn = _sensklu_c.driver()
    cb = ffi.callback('long(sensklu_t *)', lambda t: fn(t, cfn, ffi))
    mp.setattr(_sensklu_c, '_driver', (ffi, cb))
    return cb


@pytest.mark.parametrize('case', sorted(CASES))
def test_the_c_walk_is_the_python_walk(monkeypatch, case):
    """Every walk's steps but the first run's first (the stage KLU not yet
    analysed) are the C's: the same bits, states and counts."""
    left, counts, _e = _twin(monkeypatch, case)
    assert counts.get('pss.stagewalk:steps', 0) == 20 * len(left['walks']) - 1, counts


@pytest.mark.parametrize('how', ['interleaved', 'chunks'])
def test_python_steps_and_chunks_between_the_calls(monkeypatch, how):
    """A call stopped every third step (its grid cut short), the Python step
    between, the C taking over again: the state crosses both ways; and
    chunks of one to three steps, the state carried in the struct."""
    def patch(mp):
        if how == 'interleaved':
            orig = _stagewalk_c.walk
            mp.setattr(_stagewalk_c, 'walk',
                       lambda walks, j0, x, times, *a: orig(walks, j0, x, times[:j0 + 4], *a))
        else:
            mp.setattr(_stagewalk_c, 'CHUNK_FIRST', 1)
            mp.setattr(_stagewalk_c, 'CHUNK_MAX', 3)
    _left, counts, _e = _twin(monkeypatch, 'mos', patch)
    steps = counts.get('pss.stagewalk:steps', 0)
    assert (0 < steps < 59) if how == 'interleaved' else steps == 59, counts


def test_a_newton_handed_back_mid_walk(monkeypatch):
    """The C stage Newton's assembly cap at four: the steps that need more
    stop the walk before them (nothing of them made), the Python's Newton
    takes them, the C the steps after."""
    _left, counts, _e = _twin(monkeypatch, 'mos', lambda mp: mp.setattr(RC, 'MAXA', 4))
    assert counts.get('pss.stagewalk:bail_radau', 0) > 0, counts
    assert counts.get('pss.stagewalk:steps', 0) > 50, counts


def test_a_stop_after_the_newton_hands_its_result_to_the_python(monkeypatch):
    """Each walk's stage KLU numeric dropped at its start, so its first
    Newton factors afresh, and the step map handed back on that step (the
    walk's numeric not yet made): the C walk stops after a Newton that moved
    the stage KLU.  The Python's repeat of the step takes the held result --
    the same counts; held for other inputs, or not held, the repeat
    refactors what the C left, and the counts part."""
    def patch(mp, spoil=None):
        orig = LS.StructuralKLU.walk

        def fresh(self):
            self._numeric = None
            return orig(self)
        mp.setattr(LS.StructuralKLU, 'walk', fresh)

        def entry(t, cfn, ffi):
            if t.num == ffi.NULL:
                return _sensklu_c._CODES['SENSKLU_BAIL_FLAGS']
            return cfn(t)
        _sens_entry(mp, entry)
        if spoil is not None:
            walk = _stagewalk_c.walk

            def spoiled(walks, *a):
                got = walk(walks, *a)
                rec = walks._transient().__dict__['_radau_c']
                if rec.pre is not None:
                    rec.pre = spoil(rec.pre)
                return got
            mp.setattr(_stagewalk_c, 'walk', spoiled)
    _left, counts, _e = _twin(monkeypatch, 'mos', patch)
    assert counts.get('pss.stagewalk:bail_sens', 0) == counts.get('radau_c:held', 0) == 2
    a, _ca, _ea = _run(monkeypatch, 'mos', False, patch)
    for spoil in (lambda pre: None, lambda pre: ((pre[0][0][::-1],) + pre[0][1:],) + pre[1:]):
        b, cb, _eb = _run(monkeypatch, 'mos', True, lambda mp, s=spoil: patch(mp, s))
        assert not cb.get('radau_c:held')
        assert a['stage_klu'] != b['stage_klu']


def test_a_step_map_not_finite_on_dropped_columns(monkeypatch):
    """A NaN in the 30th step map's carried column: the step made, and
    `_LiveColumns.advance` takes the step at full width as the Python
    does, every column carried from there."""
    def patch(mp):
        calls = []

        def entry(t, cfn, ffi):
            st = cfn(t)
            calls.append(1)
            if len(calls) == 30:
                ffi.cast('double *', t.X)[2 * t.m] = float('nan')
            return st
        _sens_entry(mp, entry)
    _left, counts, _e = _twin(monkeypatch, 'mos', patch)
    assert counts.get('pss.stagewalk:nonfinite', 0) == 1 == counts.get('pss.live:nonfinite', 0)


def test_a_singular_step_map_raises_as_the_python(monkeypatch):
    """The 30th step map's block singular: the step made, the walk's
    numeric and counts as `block_step` commits them, the Python's
    LinAlgError -- the same state behind it."""
    def patch(mp):
        calls = []

        def entry(t, cfn, ffi):
            st = cfn(t)
            calls.append(1)
            return _sensklu_c._CODES['SENSKLU_SINGULAR'] if len(calls) == 30 else st
        _sens_entry(mp, entry)
    _left, counts, err = _twin(monkeypatch, 'mos', patch)
    assert err == 'Singular matrix' and counts.get('pss.stagewalk:singular', 0) == 1


def test_what_declines_leaves_the_python_walk(monkeypatch):
    """Switched off, a piece it stands in for patched, the cost transform:
    no step in C, the same walk."""
    _left, counts, _e = _run(monkeypatch, 'mos', False)
    assert counts.get('pss.stagewalk:off', 0) > 0 and not counts.get('pss.stagewalk:steps')
    orig = transient.Transient._pred_note

    def patch(mp):
        mp.setattr(transient.Transient, '_pred_note',
                   lambda self, *a, **k: orig(self, *a, **k))
    _left, counts, _e = _twin(monkeypatch, 'mos', patch)
    assert counts.get('pss.stagewalk:patched', 0) > 0 and not counts.get('pss.stagewalk:steps')
    monkeypatch.setitem(CASES, 'mos-transform', (CASES['mos'][0], 1e-6, {}))
    _left, counts, _e = _twin(monkeypatch, 'mos-transform')
    assert counts.get('pss.stagewalk:method', 0) > 0 and not counts.get('pss.stagewalk:steps')


def test_a_one_ulp_change_in_the_c_walk_shows(monkeypatch):
    """The C's points are the walk's: one ulp planted in each call's first
    point, and the orbit is another."""
    a, _ca, _ea = _run(monkeypatch, 'mos', True)

    def patch(mp):
        ffi, cfn, dgemm, dgemv = _stagewalk_c.driver()

        def planted(t):
            st = cfn(t)
            if t.j:
                x = ffi.cast('double *', t.Xout)
                x[0] = float(np.nextafter(x[0], np.inf))
            return st
        mp.setattr(_stagewalk_c, '_driver', (ffi, planted, dgemm, dgemv))
    b, cb, _eb = _run(monkeypatch, 'mos', True, patch)
    assert cb.get('pss.stagewalk:steps', 0) > 0 and a['orbit'] != b['orbit']
