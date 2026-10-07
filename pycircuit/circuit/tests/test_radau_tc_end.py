"""THE STEP'S END PASSES IN THE TRANSFORM'S C CALL (speed round 12, stage
5a, `_tran_radau_tc`): once it converges the C makes the passes
`_stage_end_passes` makes next -- all four at the last stage, `C` and `G` at
the first two where the shooting reads them -- and `_stage_end_passes` takes
them for those very stages, the watch counter unmoved and `passes`'s
readiness holding as stamped, making the counters each of its calls makes.
The same results, the same counts but the hand-offs' own; every way a
hand-off is dropped answers as the calls."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_core, circuit
from pycircuit.circuit import _tran_radau_tc as TC
from pycircuit.circuit._tran_radau import _RadauStages
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_psp_limit_c import _stage
from pycircuit.circuit.tests.test_radau_tc import _on
from pycircuit.circuit.transient import Transient

#: the hand-offs' own counts: nothing else moves
OWN = ('radau_tc:end', 'radau.fuse:handed', 'radau.fuse:unhanded')


def _pss(on):
    from pycircuit.circuit.shooting import PSS
    old = TC.END
    TC.END = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            p = PSS(_stage(), method='radau', reltol=1e-8)
            before = _paths.snapshot()
            p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
            d = _paths.since(before)
    finally:
        TC.END = old
    out = (np.asarray(p.waveform[1], float).tobytes(), np.asarray(p._monodromy, float).tobytes(),
           sorted(str(w.message) for w in W))
    return out, {k: v for k, v in d.items() if not k.startswith('once:')}


def _transient(on, **kw):
    old = TC.END
    TC.END = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(_stage(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                           radau_transform=True)
            before = _paths.snapshot()
            res = tr.solve(tend=6e-7, timestep=2e-8, **kw)
            d = _paths.since(before)
    finally:
        TC.END = old
    out = (np.asarray(res.x, float).tobytes(), sorted(str(w.message) for w in W))
    return out, {k: v for k, v in d.items() if not k.startswith('once:')}


def _moved_but_own(da, db):
    keys = (set(da) | set(db)) - set(OWN)
    return {k: (da.get(k, 0), db.get(k, 0)) for k in sorted(keys) if da.get(k, 0) != db.get(k, 0)}


def test_a_radau_pss_is_the_same_with_the_end_passes_off():
    _on()
    a, da = _pss(False)
    b, db = _pss(True)
    assert a == b
    assert _moved_but_own(da, db) == {}


@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
def test_a_transform_transient_is_the_same_with_the_end_passes_off(fixed):
    _on()
    kw = {'fixed_timestep': True} if fixed else {}
    a, da = _transient(False, **kw)
    b, db = _transient(True, **kw)
    assert a == b
    assert _moved_but_own(da, db) == {}


def test_the_end_passes_are_handed_over():
    """Served: on the PSP stage's radau PSS every step's end passes come from
    the C (the parent: none -- `passes` called for each)."""
    _on()
    _b, db = _pss(True)
    assert db.get('radau_tc:end', 0) > 0, db
    assert db.get('radau.fuse:handed', 0) == db['radau_tc:end'], db
    assert db.get('radau.fuse:unhanded', 0) == 0, db


def test_no_end_passes_where_nothing_takes_them(monkeypatch):
    """The stage paths' passes switched off (`_tran_radau.STAGE_FUSE`):
    `_stage_end_passes` makes none, and the C none either -- it served, the
    hand-off never made."""
    from pycircuit.circuit import _tran_radau
    _on()
    monkeypatch.setattr(_tran_radau, 'STAGE_FUSE', False)
    _b, db = _pss(True)
    assert db.get('radau_tc:served', 0) > 0, db
    assert db.get('radau_tc:end', 0) == 0 and db.get('radau.fuse:handed', 0) == 0, db


@pytest.mark.parametrize('how', ['watch', 'held', 'token', 'partial'])
def test_a_dropped_hand_off_is_the_calls(how, monkeypatch):
    """The watch counter moved, `passes`'s readiness not as stamped, another
    stage list, a memo record holding some of the last stage's passes: the
    hand-off is dropped (counted) and the calls run -- the same results and
    counts as with the end passes off, the interference the same."""
    _on()
    orig = _RadauStages._stage_end_passes

    def interfere(self, Y):
        hand = self.__dict__.get('_tc_end')
        if hand is not None:
            if how == 'watch':
                self.__dict__['_tc_end'] = (hand[0], hand[1] - 1, *hand[2:])
            elif how == 'token':
                self.__dict__['_tc_end'] = (list(hand[0]), *hand[1:])
        if how == 'partial' and self._memo_ok():
            x = Y[-1]
            self._memo_put(x, {'q': np.asarray(self.cir.q(x, self.epar), float)})
        return orig(self, Y)
    monkeypatch.setattr(_RadauStages, '_stage_end_passes', interfere)
    if how == 'held':
        monkeypatch.setattr(_tran_core, 'passes_held', lambda tr, core: False)
    a, da = _pss(False)
    b, db = _pss(True)
    assert a == b
    assert _moved_but_own(da, db) == {}
    assert db.get('radau.fuse:unhanded', 0) > 0 and db.get('radau.fuse:handed', 0) == 0, db
