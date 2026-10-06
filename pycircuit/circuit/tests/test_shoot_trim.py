"""A stage walk's per-step trims (speed round 10, B3.4): no limit sync where
no element keeps limiting state (`_InnerTransient._sync_limit_at`), the
coupled stage system in one pass (`_pss_walks._stage_block`), the step's
reduced `Jf`/`Geq`/`C` left to the walk that reads them (`_walk_lmm`), and
the reference node inserted by `analysis.insert_row`, and (speed round
11) each stage's `C` and `G` read once (`_pss_walks._stage_reads`).  Every
PSS the same
-- waveform, monodromy, warnings -- with the trims off, under every stage
method and gear; the one pass against the block loop on drawn and special
blocks, its fall-back warning from the loop's own line; the sync made
wherever a limiter keeps state; the three matrices present after a
multistep walk; the old insertion's values."""
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _limiting, _paths, circuit
from pycircuit.circuit.shooting import _pss_inner, _pss_walks
from pycircuit.circuit.tests.test_newton_c import _mos_chain
from pycircuit.circuit.tests.test_psp_limit_c import _stage


def _pss(build, method, on):
    from pycircuit.circuit.shooting import PSS
    old = _pss_inner.SHOOT_TRIM
    _pss_inner.SHOOT_TRIM = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            p = PSS(build(), method=method, reltol=1e-8)
            before = _paths.snapshot()
            p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
            d = _paths.since(before)
    finally:
        _pss_inner.SHOOT_TRIM = old
    wf = p.waveform
    mono = None if p._monodromy is None else np.asarray(p._monodromy, float).tobytes()
    return (np.asarray(wf[0], float).tobytes(), np.asarray(wf[1], float).tobytes(), mono,
            sorted(str(w.message) for w in W)), d, p


CASES = {
    'psp-radau': (_stage, 'radau', ('pss.sync:skipped', 'pss.block:served')),
    'mos-radau': (lambda: _mos_chain(cap=1e-13), 'radau', ('pss.block:served',)),
    'psp-trbdf2': (_stage, 'trbdf2', ('pss.sync:skipped',)),
    'psp-esdirk43': (_stage, 'esdirk43', ('pss.sync:skipped',)),
    'psp-gear': (_stage, 'gear', ()),
}


@pytest.mark.parametrize('case', list(CASES))
def test_a_pss_is_the_same_with_the_trims_off(case):
    build, method, keys = CASES[case]
    a, _da, pa = _pss(build, method, False)
    b, db, pb = _pss(build, method, True)
    assert a == b
    for k in keys:
        assert db.get(k, 0) > 0, (k, db)
    if not keys:
        ## (a multistep walk: nothing trimmed, its reduced matrices made)
        assert not any(k.startswith(('pss.sync:', 'pss.block:')) for k in db), db
        assert pb._Jf is not None and pb._C is not None and pb._Geq is not None
    else:
        ## (a stage walk: the matrices only `_walk_lmm` reads are not made)
        assert pb._Jf is None and pb._C is None and pb._Geq is None
        assert pa._Jf is not None
    assert not pb.__dict__.get('_skip_step_mats')


def test_the_flag_is_down_after_a_step_that_raises(monkeypatch):
    from pycircuit.circuit.shooting import PSS
    p = PSS(_stage(), method='radau', reltol=1e-8)
    p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
    seen = []

    def boom(self, *a, **k):
        seen.append(bool(self.__dict__.get('_skip_step_mats')))
        raise RuntimeError('planted')
    monkeypatch.setattr(type(p), 'solve_timestep', boom)
    with pytest.raises(RuntimeError, match='planted'):
        p._walk_stage(np.zeros(p.cir.n - 1), 1e-6, np.linspace(0.0, 1e-6, 41),
                      np.full(40, 1e-6 / 40))
    assert seen == [True]
    assert not p.__dict__.get('_skip_step_mats')


## -- the stage block --------------------------------------------------------------------

def _loop(Cs, Gs, h, A):
    """`_stage_step`'s block loop, as it stood (the oracle)."""
    s, m = A.shape[0], Cs[0].shape[0]
    Jb = np.zeros((s * m, s * m))
    for i in range(s):
        for jj in range(s):
            blk = h * A[i, jj] * Gs[jj]
            if i == jj:
                blk = Cs[i] + blk
            Jb[i * m:(i + 1) * m, jj * m:(jj + 1) * m] = blk
    return Jb


def _assembled(Cs, Gs, h, A, trim):
    """The stage system as `_stage_step` assembles it, and its warnings
    with their lines -- the one pass where it serves, else the loop."""
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        Jb = _pss_walks._stage_block(Cs, Gs, h, A) if trim else None
        if Jb is None:
            Jb = _loop(Cs, Gs, h, A)
    return Jb.tobytes(), Jb.flags.c_contiguous, [(str(w.message)) for w in W]


#: (NaNs of three payloads: where both operands of a sum are NaN the result
#: is one operand's, so only distinct payloads make the operand order seen)
_NANS = tuple(float(np.frombuffer(np.uint64(b).tobytes(), dtype=np.float64)[0])
              for b in (0x7FF8000000000000, 0x7FF8000000000123, 0xFFF8000000000456))
SPECIAL = (0.0, -0.0, np.inf, -np.inf, 5e-324, 1e-300, 1e300, -1e300, 1.0) + _NANS
_A = np.array([[0.19681547722366, -0.06553542585020, 0.02377097434822],
               [0.39442431473909, 0.29207341166523, -0.04154875212600],
               [0.37640306270047, 0.51248582618842, 0.11111111111111]])


@settings(deadline=None)
@given(data=st.data())
def test_the_one_pass_is_the_block_loop(data):
    """Drawn blocks (signed zeros, NaN, infinities, subnormals, overflowing
    products) and steps: the same system's bytes, C-contiguous, and the
    same warnings -- an operation that would warn handing the system to
    the loop, which warns from its own line as before."""
    m = data.draw(st.integers(1, 5))
    elems = st.one_of(st.floats(-1e3, 1e3), st.floats(-1e-12, 1e-12), st.sampled_from(SPECIAL))
    Cs = [np.array(data.draw(st.lists(elems, min_size=m * m, max_size=m * m)),
                   dtype=float).reshape(m, m) for _ in range(3)]
    Gs = [np.array(data.draw(st.lists(elems, min_size=m * m, max_size=m * m)),
                   dtype=float).reshape(m, m) for _ in range(3)]
    h = data.draw(st.sampled_from((1e-9, 2.5e-8, 1.0, 1e300, 0.0, -0.0, 5e-324)))
    assert _assembled(Cs, Gs, h, _A, True) == _assembled(Cs, Gs, h, _A, False)


def test_the_one_pass_declines_what_it_cannot_make():
    m = 3
    Cs = [np.eye(m)] * 3
    Gs = [np.ones((m, m))] * 3
    for args, why in (((Cs, [g.astype(np.float32) for g in Gs], 1e-9, _A), 'kind'),
                      ((Cs[:2], Gs, 1e-9, _A), 'kind'),
                      ((Cs, [np.full((m, m), np.nan)] * 3, 1e-9, _A), 'nonfinite'),
                      ((Cs, Gs, np.inf, _A), 'nonfinite'),
                      ((Cs, [np.full((m, m), 1e300)] * 3, 1e300, _A), 'fp')):
        before = _paths.snapshot()
        assert _pss_walks._stage_block(*args) is None
        assert _paths.since(before) == {'pss.block:' + why: 1}


## -- the sync and the insertion ---------------------------------------------------------

def test_the_sync_is_made_wherever_a_limiter_keeps_state(monkeypatch):
    """Skipped (counted) where no element keeps limiting state; made where
    one does, or where circuit-level limiting is on, or with the trims off."""
    from pycircuit.circuit.shooting import PSS
    p = PSS(_stage(), method='radau', reltol=1e-8)
    p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
    tr = p._transient()
    calls = []
    monkeypatch.setattr(tr.cir, 'limit', lambda x, x0, epar: calls.append(1) or x,
                        raising=False)
    x = np.zeros(tr.cir.n)

    def synced(**state):
        calls.clear()
        before = _paths.snapshot()
        with monkeypatch.context() as mp:
            for k, v in state.items():
                if k == 'lims':
                    mp.setitem(tr.__dict__, '_stateful_lims', v)
                elif k == 'level':
                    mp.setattr(_limiting, 'CIRCUIT_LEVEL', v)
                else:
                    mp.setattr(_pss_inner, 'SHOOT_TRIM', v)
            p._sync_limit_at(x)
        return len(calls), _paths.since(before).get('pss.sync:skipped', 0)
    assert synced(lims=[]) == (0, 1)
    assert synced(lims=[object()]) == (1, 0)
    assert synced(lims=[], level=True) == (1, 0)
    assert synced(lims=[], trim=False) == (1, 0)


@pytest.mark.parametrize('x', [np.array([1.0, -2.0, 3.5]), np.array([-0.0, 5e-324]),
                               np.array([1.0 + 2.0j, 3.0]), np.array([np.nan]), np.array([])],
                         ids=['float', 'zeros', 'complex', 'nan', 'empty'])
@pytest.mark.parametrize('iref', [0, 1, -1])
def test_the_insertion_is_the_old_one(x, iref):
    """`_insert_refnode` (now `analysis.insert_row`) against the old
    `concatenate`: the same dtype and bytes, a fresh array."""
    import types
    tk = circuit.numeric
    p = types.SimpleNamespace(irefnode=len(x) if iref == -1 else min(iref, len(x)), toolkit=tk)
    old = tk.concatenate((x[:p.irefnode], tk.array([0.0]), x[p.irefnode:]))
    new = _pss_inner._InnerTransient._insert_refnode(p, x)
    assert new.dtype == old.dtype and new.tobytes() == old.tobytes() and new is not x


## -- the stage reads (speed round 11) ---------------------------------------------------

def _solved(build=_stage):
    from pycircuit.circuit.shooting import PSS
    p = PSS(build(), method='radau', reltol=1e-8)
    p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
    return p


def test_the_reads_are_the_readers(monkeypatch):
    """In every coupled stage step of a radau PSS: each stage's `C` and `G`
    as `_C_at` and `_G_at` read them at the reduced stage, bit for bit, and
    the sync's skips counted as theirs are."""
    real = _pss_walks._stage_reads
    seen = []

    def both(walks, Yf):
        before = _paths.snapshot()
        got = real(walks, Yf)
        d_new = _paths.since(before)
        if got is not None:
            iref = walks.irefnode
            Ys = [walks.toolkit.concatenate((y[:iref], y[iref + 1:])) for y in Yf]
            before = _paths.snapshot()
            Cs = [np.asarray(walks._C_at(y)) for y in Ys]
            Gs = [np.asarray(walks._G_at(y)) for y in Ys]
            d_old = _paths.since(before)
            seen.append((all(a.dtype == b.dtype and a.shape == b.shape
                             and a.tobytes() == b.tobytes()
                             for a, b in zip(got[0] + got[1], Cs + Gs)),
                         d_new.get('pss.sync:skipped'), d_old.get('pss.sync:skipped')))
        return got
    monkeypatch.setattr(_pss_walks, '_stage_reads', both)
    _solved()
    assert len(seen) >= 40 and all(x == (True, 3, 3) for x in seen), seen[:3]


def test_the_reads_serve_a_radau_pss(monkeypatch):
    """The PSP stage's radau PSS: each coupled stage step's `C` and `G` read
    once (`pss.reads:served`, one a sensitivity step) -- the Python step's,
    the C step map off (`_sens_c`: where it serves it reads the memo
    itself).  (Apart from the trims' comparison above: the reads need the
    step's memo, which the other fast paths fill -- with those off, this
    fails and that compares.)"""
    from pycircuit.circuit.shooting import _sens_c
    monkeypatch.setattr(_sens_c, 'SENS_C', False)
    before = _paths.snapshot()
    _solved()
    d = _paths.since(before)
    assert d.get('pss.reads:served', 0) >= 40, d


def _declined(p, Yf):
    before = _paths.snapshot()
    assert _pss_walks._stage_reads(p, Yf) is None
    return sorted(k for k in _paths.since(before) if k.startswith('pss.reads:'))


def test_the_reads_decline_where_the_readers_read_otherwise(monkeypatch):
    """A reader patched (on the instance, on the class), a junction (PCNR's
    `G`), limiting state kept (a stateful limiter, the circuit level), a
    full stage whose reference entry is not +0.0 (not the state the readers
    rebuild), a stage the memo lacks: declined, counted -- the readers run."""
    p = _solved()
    iref = p.irefnode
    Yf = [np.array(y, dtype=float) for y in p._transient().last_step.Y]
    bad = [y.copy() for y in Yf]
    bad[0][iref] = -0.0
    assert _declined(p, bad) == ['pss.reads:state']
    odd = [y + 1.0 for y in Yf]
    for y in odd:
        y[iref] = 0.0
    assert _declined(p, odd) == ['pss.reads:memo']
    with monkeypatch.context() as mp:
        mp.setitem(p.__dict__, '_G_at', lambda x: None)
        assert _declined(p, Yf) == ['pss.reads:patched']
    with monkeypatch.context() as mp:
        mp.setattr(type(p), '_C_at', lambda self, x: None)
        assert _declined(p, Yf) == ['pss.reads:patched']
    with monkeypatch.context() as mp:
        mp.setitem(p.__dict__, '_pcnr_junctions_cache', [object()])
        assert _declined(p, Yf) == ['pss.reads:junctions']
    with monkeypatch.context() as mp:
        mp.setattr(_limiting, 'CIRCUIT_LEVEL', True)
        assert _declined(p, Yf) == ['pss.reads:limits']
    with monkeypatch.context() as mp:
        mp.setitem(p._transient().__dict__, '_stateful_lims', [object()])
        assert _declined(p, Yf) == ['pss.reads:limits']
