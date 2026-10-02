"""The stage predictor's bookkeeping trimmed bit for bit, and its
Vandermonde weights kept (speed round 4, stage D; 2026-10-02).

`_predict_state` sorted the newest nodes twice, deduplicated the node
times through a generator per pair, and solved the 3x3 Vandermonde
system at every step -- 63 us of a 864 us step in-run on the PSP stage.
The dedupe is a plain loop with the same rule and order, the newest are
sorted once, and `_fit`'s weights, a pure function of the normalised
node times, are kept on the instance keyed by their bits: the shooting
re-walks the same grid every iteration (`PRED_WEIGHT_MEMO = False` is the
old behaviour, for the comparison).
"""
import numpy as np
import pytest

from pycircuit.circuit import PSS, _tran_predictor
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import C, R, SubCircuit, VSin, gnd
from pycircuit.circuit.integrator import (
    ESDIRK43Integrator,
    Gear2Integrator,
    RadauIIA3Integrator,
)
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.transient import Transient


def _rect():
    cm.default_toolkit = numeric
    c = SubCircuit()
    for n in ('in', 'a', 'k'):
        c.add_node(n)
    c['vin'] = VSin('in', gnd, v=0.0, va=1.0, freq=1e6)
    c['rs'] = R('in', 'a', r=100.0)
    c['d'] = eh.DiodeSpiceHdl('a', 'k')
    c['rl'] = R('k', gnd, r=1e3)
    c['cl'] = C('k', gnd, c=1e-10)
    return c


def _reference_predict(tr, ttarget, extra=(), deg=None):
    """`_predict_state` as it was before stage D, written out: the pairwise
    generator dedupe, the Vandermonde solve every call, the newest sorted
    twice."""
    if tr.stage_predictor == 'off':
        return None
    nodes = [(float(tt), np.asarray(xx, dtype=float)) for tt, xx in extra]
    nodes.extend(getattr(tr, '_pred_hist', ()) or ())
    if len(nodes) < 2:
        return None
    if deg is None:
        deg = tr._pred_degree()
    nodes.sort(key=lambda e: abs(e[0] - ttarget))
    uniq = []
    for e in nodes:
        if not any(abs(e[0] - u[0]) <= 1e-13 * max(abs(e[0]), 1.0)
                   for u in uniq):
            uniq.append(e)
    nodes = uniq
    if len(nodes) < 2:
        return None
    nfit = min(int(deg) + 1, len(nodes))
    take = nodes[:nfit]

    def _fit(sub, tat):
        tv = np.array([e[0] for e in sub], dtype=float)
        scale = float(np.max(np.abs(tv - tat)))
        if not np.isfinite(scale) or scale <= 0.0:
            return None
        tau = (tv - tat) / scale
        n = len(tau)
        rhs = np.zeros(n)
        rhs[0] = 1.0
        try:
            w = np.linalg.solve(np.vander(tau, n, increasing=True).T, rhs)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(w)):
            return None
        return w @ np.array([e[1] for e in sub], dtype=float)

    newest = sorted(take, key=lambda e: -e[0])[:2]
    motion = np.abs(newest[0][1] - newest[1][1]) if len(newest) == 2 \
        else None
    if motion is None:
        return None
    pred = _fit(take, ttarget)
    if pred is None:
        return None
    newest_first = sorted(take, key=lambda e: -e[0])
    xref = newest_first[0][1]
    dt_last = abs(newest_first[0][0] - newest_first[1][0])
    ratio = (abs(ttarget - newest_first[0][0]) / dt_last) if dt_last > 0 \
        else 1.0
    w = tr.PRED_CLAMP * motion * max(ratio, 1.0)
    out = np.clip(pred, xref - w, xref + w)
    for row, _m, _o in (getattr(tr, '_periodic_rows', None) or ()):
        out[row] = xref[row]
    return out


def test_the_prediction_is_the_old_ones_bit_for_bit():
    """Random node histories (up to the eight kept), targets ahead, behind
    and amid them, `extra` nodes repeating a recorded time to the bit and
    within the dedupe's tolerance, a duplicate target: the same prediction
    as the transliterated old code, memo hits included."""
    tr = Transient(_rect(), toolkit=numeric, uic=True)
    tr.solve(tend=4e-7, timestep=2e-8, fixed_timestep=True)
    rng = np.random.default_rng(11)
    n = tr.cir.n
    for k in range(300):
        m = int(rng.integers(2, 9))
        times = np.sort(rng.uniform(0.0, 1e-6, m))
        tr._pred_hist = [(float(t), rng.standard_normal(n)) for t in times]
        tt = float(rng.uniform(-2e-7, 1.4e-6))
        extra = []
        if k % 3 == 1:
            t0 = tr._pred_hist[-1][0]
            extra = [(t0, rng.standard_normal(n)),
                     (t0 * (1.0 + 2e-14), rng.standard_normal(n))]
        elif k % 3 == 2:
            extra = [(tt, rng.standard_normal(n))]
        deg = None if k % 4 else int(rng.integers(1, 5))
        a = _reference_predict(tr, tt, extra, deg)
        b = tr._predict_state(tt, extra=extra, deg=deg)
        if a is None or b is None:
            assert a is None and b is None
        else:
            assert a.tobytes() == b.tobytes(), k
    ## and a repeated geometry is served from the memo with the same bits
    tr._pred_hist = [(1e-7, np.ones(n)), (2e-7, 2 * np.ones(n)), (3e-7, 4 * np.ones(n))]
    a = tr._predict_state(4e-7)
    assert len(tr._pred_wmemo) >= 1
    assert tr._predict_state(4e-7).tobytes() == a.tobytes()


def test_a_shooting_walk_serves_its_weights_after_the_first_period(monkeypatch):
    """The shooting re-walks the same grid every iteration: after the first
    period every `_fit` is a memo hit.  `np.vander` is the predictor's
    alone in the library, so its calls count the solves made, against the
    fits asked for (`Transient._fit`, which the parent did not have: it
    solved at every fit)."""
    fits, vanders = [0], [0]
    real_fit = Transient._fit
    real_vander = np.vander

    def counting_fit(self, sub, tat):
        fits[0] += 1
        return real_fit(self, sub, tat)

    def counting_vander(*a, **k):
        vanders[0] += 1
        return real_vander(*a, **k)
    monkeypatch.setattr(Transient, '_fit', counting_fit)
    monkeypatch.setattr(np, 'vander', counting_vander)
    p = PSS(_rect(), method='gear', reltol=1e-6)
    with quiet(AccuracyWarning):
        p.solve(period=1e-6, timestep=1e-6 / 32, maxiterations=40)
    assert fits[0] > 64                       # more than two periods' fits
    assert vanders[0] * 2 < fits[0], (fits[0], vanders[0])


@pytest.mark.parametrize('memo', [True, False])
def test_transients_and_a_pss_are_byte_identical_with_the_memo_off(memo, monkeypatch):
    """Gear, Radau and ESDIRK43 on the adaptive controller, and a shooting
    PSS: the same bytes with the weights kept and recomputed."""
    def run():
        out = []
        for integ in (Gear2Integrator, RadauIIA3Integrator, ESDIRK43Integrator):
            tr = Transient(_rect(), toolkit=numeric, integrator=integ())
            res = tr.solve(tend=1.5e-6, timestep=1e-8)
            out.append(np.asarray(res.x, float).tobytes())
            out.append(np.asarray(res.sweep_values, float).tobytes())
        p = PSS(_rect(), method='gear', reltol=1e-6)
        with quiet(AccuracyWarning):
            p.solve(period=1e-6, timestep=1e-6 / 32, maxiterations=40)
        out.append(np.asarray(p.waveform[1], float).tobytes())
        return out
    monkeypatch.setattr(_tran_predictor, 'PRED_WEIGHT_MEMO', True)
    a = run()
    monkeypatch.setattr(_tran_predictor, 'PRED_WEIGHT_MEMO', memo)
    b = run()
    assert a == b
