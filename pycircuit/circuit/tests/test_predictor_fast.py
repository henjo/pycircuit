"""The stage predictor's multistep fast path (`_tran_predictor.PRED_FAST`,
speed round 8, stage 2b): bit for bit the general path."""
import numpy as np
from hypothesis import given
from hypothesis import strategies as st

from pycircuit.circuit import _paths, _tran_predictor, circuit
from pycircuit.circuit.elements import C, R, SubCircuit, VSin, gnd
from pycircuit.circuit.transient import Transient


def _rc():
    c = SubCircuit()
    c['vs'] = VSin('in', gnd, va=1.0, freq=1e3)
    c['R'] = R('in', 'out', r=1e3)
    c['C'] = C('out', gnd, c=1e-7)
    return c


def _both(tr, ttarget, deg=None):
    """The prediction with the fast path on and off (the same instance,
    each with its own weights: the memo emptied before each)."""
    out = []
    for on in (True, False):
        old = _tran_predictor.PRED_FAST
        _tran_predictor.PRED_FAST = on
        try:
            tr.__dict__.pop('_pred_fast', None)
            tr.__dict__.pop('_pred_wmemo', None)
            out.append(tr._predict_state(ttarget, deg=deg))
        finally:
            _tran_predictor.PRED_FAST = old
    return out


def _same(a, b):
    if a is None or b is None:
        return a is None and b is None
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


_finite = st.floats(-1e6, 1e6, allow_nan=False, allow_infinity=False)


@given(times=st.lists(st.floats(0.0, 1e-3, allow_nan=False), min_size=2, max_size=8, unique=True),
       ahead=st.floats(1e-12, 1e-3), deg=st.sampled_from([1, 2, 3, 4]),
       data=st.data())
def test_the_fast_path_is_the_general_path_bit_for_bit(times, ahead, deg, data):
    """Any history behind the target (times close enough to merge in the
    1e-13 dedupe included), any degree: the same bytes, or both None."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    n = 4
    times = sorted(times)
    tr._pred_hist = [(float(t), np.array(data.draw(st.lists(_finite, min_size=n, max_size=n))))
                     for t in times]
    a, b = _both(tr, times[-1] + ahead, deg)
    assert _same(a, b)


def test_near_duplicate_times_and_a_wrapping_row():
    """Times inside the general path's 1e-13 dedupe (but outside the
    history's own 1e-14) and a periodic row: still the same bytes."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    t0 = 1e-3
    tr._pred_hist = [(t0 - 3e-6, np.array([0.1, -0.2, 0.3])),
                     (t0 - 2e-6, np.array([0.2, -0.1, 0.25])),
                     (t0 - 1e-6, np.array([0.3, 0.0, 0.2])),
                     (t0 - 1e-6 + 5e-17, np.array([0.31, 0.01, 0.21]))]
    tr._periodic_rows = [(1, 1.0, 0.0)]
    for deg in (1, 2, 3):
        a, b = _both(tr, t0, deg)
        assert _same(a, b)


def test_a_tie_in_the_distances_takes_the_general_path():
    """Two times closer than their distance's rounding -- 0 and 1.6e-45,
    behind a target 1.4e-3 away (the deep hypothesis profile's draw,
    2026-10-04) -- round to one distance: the general path's stable sort
    lists the older first and its dedupe keeps that one, so the fast path
    steps aside, counted."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    t1 = 0.00036596218663498633
    tr._pred_hist = [(0.0, np.zeros(4)), (1.5828005792270418e-45, np.array([0.0, 0.0, 0.0, 1.0])),
                     (t1, np.zeros(4))]
    before = _paths.snapshot()
    a, b = _both(tr, t1 + 0.001, 1)
    assert _same(a, b)
    assert _paths.since(before).get('pred:tie', 0) == 1


def test_the_fast_path_serves_the_multistep_step():
    """A gear transient's predictions take it: the count says so."""
    before = _paths.snapshot()
    Transient(_rc(), toolkit=circuit.numeric).solve(tend=2e-4, timestep=1e-5,
                                                    fixed_timestep=True)
    d = _paths.since(before)
    assert d.get('pred:fast', 0) > 10 and d.get('pred:general', 0) == 0


def test_an_edited_history_takes_the_general_path():
    """A history not in ascending order (a test's edit): the general path,
    counted, and the same prediction."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    tr._pred_hist = [(2e-6, np.array([0.2, 0.1, 0.0])), (1e-6, np.array([0.1, 0.0, 0.0])),
                     (3e-6, np.array([0.3, 0.2, 0.1]))]
    before = _paths.snapshot()
    a, b = _both(tr, 4e-6, 2)
    assert _same(a, b)
    assert _paths.since(before).get('pred:general', 0) == 1


def test_a_patched_fit_takes_the_general_path(monkeypatch):
    """A `_fit` patched on the class is a caller that expects its calls:
    the fast path steps aside, counted, and the patch is reached."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    tr._pred_hist = [(1e-6, np.array([0.1, 0.0, 0.0])), (2e-6, np.array([0.2, 0.1, 0.0])),
                     (3e-6, np.array([0.3, 0.2, 0.1]))]
    seen = []
    real = Transient._fit

    def counting(self, sub, tat):
        seen.append(len(sub))
        return real(self, sub, tat)
    monkeypatch.setattr(Transient, '_fit', counting)
    before = _paths.snapshot()
    tr._predict_state(4e-6, deg=2)
    assert seen == [3] and _paths.since(before).get('pred:patched', 0) == 1


## -- a coupled stage method's seeds at once, and the accept's merge (roadmap 198)

def _per_target(tr, x, ttargets):
    """`_pred_stages`'s definition: `_pred_or` per target."""
    return [tr._pred_or(np.array(x, dtype=float), t) for t in ttargets]


@given(times=st.lists(st.floats(0.0, 1e-3, allow_nan=False), min_size=1, max_size=8, unique=True),
       ahead=st.lists(st.floats(-1e-4, 1e-3), min_size=1, max_size=4), deg=st.sampled_from([1, 2, 3, 4]),
       data=st.data())
def test_the_stage_seeds_at_once_are_the_per_target_ones(times, ahead, deg, data):
    """`_pred_stages` on any history and targets -- behind or ahead of it,
    near-duplicate times included: per target the bytes `_pred_or` gives
    (its fallback a new array each), and the same path counts."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    tr.stage_predictor_degree = deg
    n = 4
    times = sorted(times)
    tr._pred_hist = [(float(t), np.array(data.draw(st.lists(_finite, min_size=n, max_size=n))))
                     for t in times]
    x = np.array(data.draw(st.lists(_finite, min_size=n, max_size=n)))
    targets = [times[-1] + a for a in ahead]
    got = []
    for fn in (tr._pred_stages, lambda x, ts: _per_target(tr, x, ts)):
        tr.__dict__.pop('_pred_fast', None)
        tr.__dict__.pop('_pred_wmemo', None)
        before = _paths.snapshot()
        out = fn(x, targets)
        got.append(([o.tobytes() for o in out], {k: v for k, v in _paths.since(before).items()
                                                  if k.startswith('pred:')}))
        assert len({id(o) for o in out}) == len(out) and all(o is not x for o in out)
    assert got[0] == got[1]


def test_a_patched_pred_or_sees_every_stage(monkeypatch):
    """A `_pred_or` patched on the class takes the per-target calls."""
    calls = []
    orig = _tran_predictor._StagePredictor._pred_or

    def counted(self, fallback, ttarget, extra=()):
        calls.append(ttarget)
        return orig(self, fallback, ttarget, extra)
    monkeypatch.setattr(_tran_predictor._StagePredictor, '_pred_or', counted)
    tr = Transient(_rc(), toolkit=circuit.numeric)
    tr.stage_predictor_degree = 2
    tr._pred_hist = [(0.0, np.zeros(4)), (1e-6, np.ones(4))]
    tr._pred_stages(np.ones(4), [1.5e-6, 2e-6, 2.5e-6])
    assert calls == [1.5e-6, 2e-6, 2.5e-6]


def _note_full(h, t, x, stages):
    """`_pred_note`'s list as the whole sort and dedupe made it."""
    h = list(h)
    for ts, ys in stages:
        h.append((float(ts), np.array(ys, dtype=float)))
    h.append((float(t), np.array(x, dtype=float)))
    h.sort(key=lambda e: e[0])
    uniq = []
    for e in h:
        if uniq and abs(e[0] - uniq[-1][0]) <= 1e-14 * max(abs(e[0]), 1.0):
            uniq[-1] = e
        else:
            uniq.append(e)
    return uniq[-Transient.PRED_HIST_MAX:]


@given(old=st.lists(st.floats(0.0, 1e-3, allow_nan=False), max_size=9, unique=True),
       new=st.lists(st.floats(0.0, 1.2e-3, allow_nan=False), min_size=1, max_size=4),
       near=st.booleans())
def test_the_accept_s_merge_is_the_whole_sort(old, new, near):
    """`_pred_note` on any history and new nodes -- after it, before it,
    among it, a new node within the dedupe's tolerance of the newest old
    one: the list the whole sort and dedupe made, entry for entry."""
    tr = Transient(_rc(), toolkit=circuit.numeric)
    hist = _note_full([], 0.0, np.zeros(2), ())[:0]
    for k, t in enumerate(sorted(old)):
        hist = _note_full(hist, t, np.full(2, float(k)), ())
    if near and hist:
        new = [hist[-1][0] * (1.0 + 1e-15)] + new
    tr._pred_hist = list(hist)
    *stages, t = new
    stages = [(ts, np.full(2, 10.0 + i)) for i, ts in enumerate(stages)]
    tr._pred_note(t, np.full(2, 99.0), stages)
    want = _note_full(hist, t, np.full(2, 99.0), stages)
    assert [(e[0], e[1].tobytes()) for e in tr._pred_hist] == [(e[0], e[1].tobytes()) for e in want]
