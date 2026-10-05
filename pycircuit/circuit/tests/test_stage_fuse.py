"""The stages' readers' passes in one core call a stage (speed round 10,
B3.2): at a Radau step's end the branch screen's `C` and the step end's
`q`, `i`, `G` fused (`_RadauStages._stage_end_passes`), and at the other
stages `C` and `G` together where the shooting walk reads `G` there at
once (`_stage_G_read`, set by `_walk_stage` around each step).  Every
transient and PSS the same with the switch off; the hint set around the
step alone; the declines.  And the branch screen's proxy, each reduction
made once: the old screen's decisions, scale and warnings on drawn and
special matrices (the old code is the oracle)."""
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _paths, _tran_core, _tran_radau, circuit
from pycircuit.circuit.integrator import Gear2Integrator, RadauIIA3Integrator
from pycircuit.circuit.tests.test_newton_c import _gp_chain, _mos_chain
from pycircuit.circuit.tests.test_psp_limit_c import _stage
from pycircuit.circuit.transient import Transient


def _ladder():
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import Diode, R, VSin, gnd
    c = SubCircuit()
    c.add_node('n0')
    c['vs'] = VSin('n0', gnd, va=1.2, freq=5e6)
    for k in range(3):
        c.add_node(f'n{k + 1}')
        c[f'r{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
        c[f'd{k}'] = Diode(f'n{k + 1}', gnd)
    return c


def _run(build, on, make=None, integ=RadauIIA3Integrator, **kw):
    old = _tran_radau.STAGE_FUSE
    _tran_radau.STAGE_FUSE = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build(), toolkit=circuit.numeric, integrator=integ(), **(make or {}))
            before = _paths.snapshot()
            res = tr.solve(**kw)
            d = _paths.since(before)
    finally:
        _tran_radau.STAGE_FUSE = old
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W)), d


CASES = {
    ## (the transform's own Newton: fused)
    'psp': (_stage, {}, 'served'),
    'mos-transform': (_mos_chain, {'radau_transform': True}, 'served'),
    ## (the dense Newton records its stages itself: not asked)
    'gp-dense': (_gp_chain, {'radau_transform': False}, None),
    ## (a stateful limiter on the transform: the memo does not record)
    'diode-ladder': (_ladder, {'radau_transform': True}, 'memo'),
}


@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
@pytest.mark.parametrize('case', list(CASES))
def test_a_radau_transient_is_the_same_with_the_fusion_off(case, fixed):
    build, make, key = CASES[case]
    kw = {'tend': 6e-7, 'timestep': 2e-8}
    if fixed:
        kw['fixed_timestep'] = True
    a, da = _run(build, False, make, **dict(kw))
    b, db = _run(build, True, make, **dict(kw))
    assert a == b
    assert da.get('radau.fuse:off', 0) == db.get('radau.fuse:' + key, 0) if key else True, (da, db)
    if key is None:
        assert not any(k.startswith('radau.fuse:') for k in db), db
    else:
        assert db.get('radau.fuse:' + key, 0) > 0, db


def _pss(on):
    from pycircuit.circuit.shooting import PSS
    old = _tran_radau.STAGE_FUSE
    _tran_radau.STAGE_FUSE = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            p = PSS(_stage(), method='radau', reltol=1e-8)
            before = _paths.snapshot()
            p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
            d = _paths.since(before)
    finally:
        _tran_radau.STAGE_FUSE = old
    wf = p.waveform
    ## (the monodromy too: built from the stages' `C` and `G` -- a strongly
    ## damped period map's waveform need not move where it does)
    return (np.asarray(wf[0], float).tobytes(), np.asarray(wf[1], float).tobytes(),
            np.asarray(p._monodromy, float).tobytes(),
            sorted(str(w.message) for w in W)), d, p


def test_a_radau_pss_is_the_same_and_its_stages_are_read_from_the_memo():
    """The PSS of the PSP stage by Radau: the same waveform and warnings;
    with the fusion on the shooting's stage `G` comes from the memo -- the
    circuit's own `G` evaluations drop -- and the hint is down after the
    solve."""
    a, da, _p = _pss(False)
    b, db, p = _pss(True)
    assert a == b
    assert db.get('radau.fuse:served', 0) > 0, db
    plan_g = lambda d: sum(v for k, v in d.items() if k.startswith(('plan.G:', 'batch.G:')))
    assert plan_g(db) < plan_g(da), (plan_g(da), plan_g(db))
    assert not p._transient().__dict__.get('_stage_G_read')


def test_the_hint_is_down_after_a_step_that_raises(monkeypatch):
    """`_walk_stage` sets the hint around each step alone: a step that
    raises leaves it down."""
    from pycircuit.circuit.shooting import PSS
    p = PSS(_stage(), method='radau', reltol=1e-8)
    p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
    tr = p._transient()
    seen = []

    def boom(self, *a, **k):
        seen.append(bool(tr.__dict__.get('_stage_G_read')))
        raise RuntimeError('planted')
    monkeypatch.setattr(type(p), 'solve_timestep', boom)
    with pytest.raises(RuntimeError, match='planted'):
        p._walk_stage(np.zeros(p.cir.n - 1), 1e-6, np.linspace(0.0, 1e-6, 41),
                      np.full(40, 1e-6 / 40))
    assert seen == [True]
    assert not tr.__dict__.get('_stage_G_read')


def test_what_the_fusion_declines(monkeypatch):
    """The switch, a memo that does not record (bypass), the core not
    serving: each counted, the readers evaluating as before."""
    a, _da = _run(_stage, False, tend=2e-7, timestep=2e-8, fixed_timestep=True)
    monkeypatch.setattr(_tran_core, 'CORE_PASSES', False)
    b, db = _run(_stage, True, tend=2e-7, timestep=2e-8, fixed_timestep=True)
    assert a == b and db.get('radau.fuse:core', 0) > 0, db
    monkeypatch.undo()
    _c, dc = _run(_stage, True, {'bypass': True, 'bypasstol': 1e-9}, tend=2e-7,
                 timestep=2e-8, fixed_timestep=True)
    assert dc.get('radau.fuse:memo', 0) > 0, dc


## -- the branch screen's proxy -----------------------------------------------------------

def _old_proxy(self, C):
    """The branch screen's cheap proxy as it stood before B3.2 (584411d3):
    the oracle."""
    ref = max(getattr(self, '_branch_cmax', 0.0),
              float(np.max(np.abs(C))))
    self._branch_cmax = ref
    d = np.abs(np.diag(C))
    if float(np.max(np.abs(C))) > self.BRANCH_SCREEN_TOL * ref \
            and float(np.max(d)) > self.BRANCH_SCREEN_TOL * ref:  # noqa: SIM102 (the old code, verbatim)
        if float(np.min(d[d > 0.0]) if np.any(d > 0.0) else 0.0) \
                > 1e-6 * ref:
            return False, None
    nrm = float(np.max(np.abs(C)))
    if nrm <= self.BRANCH_SCREEN_TOL * ref:
        return True, None
    return 'svd', None


def _new_screen(tr, C, monkeypatch):
    """`_branch_screen` on a planted `C`, stopped where the SVD would run."""
    with monkeypatch.context() as mp:
        mp.setattr(tr, '_C_at_state', lambda x: C)
        mp.setattr(tr, '_branch_structural_rank', lambda: (len(C), 1.0))
        mp.setattr(np.linalg, 'svd', lambda *a, **k: (_ for _ in ()).throw(_SVD()))
        try:
            return tr._branch_screen(np.zeros(len(C)))
        except _SVD:
            return 'svd', None


class _SVD(Exception):
    pass


SPECIAL = (0.0, -0.0, np.nan, np.inf, -np.inf, 5e-324, 1e-300, 1e300, 1.0, -1.0)


def _proxy_pair(C, cmax):
    """(the old proxy's decision and kept scale, the screen's) on one `C`."""
    tr = _screen_tr()
    if cmax is None:
        tr.__dict__.pop('_branch_cmax', None)
    else:
        tr._branch_cmax = cmax
    old = type('O', (), {'BRANCH_SCREEN_TOL': tr.BRANCH_SCREEN_TOL})()
    if cmax is not None:
        old._branch_cmax = cmax
    with warnings.catch_warnings(record=True) as W1:
        warnings.simplefilter('always')
        want = _old_proxy(old, C)
    from _pytest.monkeypatch import MonkeyPatch
    mp = MonkeyPatch()
    try:
        with warnings.catch_warnings(record=True) as W2:
            warnings.simplefilter('always')
            got = _new_screen(tr, C, mp)
    finally:
        mp.undo()
    return ((want[0], np.float64(old._branch_cmax).tobytes(), sorted(str(w.message) for w in W1)),
            (got[0], np.float64(tr._branch_cmax).tobytes(), sorted(str(w.message) for w in W2)))


def _enum():
    """Every branch of the proxy, enumerated: a zero beside large diagonal
    entries (the positive part's minimum), a tiny positive one, signed
    zeros, a collapsed matrix, a running scale above the matrix's, NaN and
    infinities."""
    out = []
    for diag in ((1.0, 0.0, 2.0), (1.0, -0.0, 2.0), (1.0, 1e-9, 2.0), (0.0, 0.0, 0.0),
                 (1e-300, 1.0, 1.0), (np.nan, 1.0, 1.0), (np.inf, 1.0, 1.0), (1.0, 1.0, 1.0)):
        for off in (0.0, 0.5, 1e-12, np.nan):
            C = np.full((3, 3), off)
            C[np.diag_indices(3)] = diag
            for cmax in (None, 0.0, 1.0, 1e6, np.nan):
                out.append((C, cmax))
    out.append((np.zeros((2, 2)), None))
    out.append((np.full((2, 2), -0.0), 1.0))
    return out


@pytest.mark.parametrize('k', range(len(_enum())))
def test_the_screens_proxy_is_the_old_one_on_every_branch(k):
    C, cmax = _enum()[k]
    want, got = _proxy_pair(C, cmax)
    assert got == want


@settings(deadline=None)
@given(data=st.data())
def test_the_screens_proxy_is_the_old_one(data):
    """Drawn `C` (zeros, signed zeros, NaN, infinities, subnormals, collapsed
    diagonals) and drawn running scales: the same decision, the same scale
    kept, the same warnings."""
    n = data.draw(st.integers(1, 6))
    vals = data.draw(st.lists(st.one_of(st.sampled_from(SPECIAL),
                                        st.floats(-1e3, 1e3), st.floats(-1e-12, 1e-12)),
                              min_size=n * n, max_size=n * n))
    C = np.array(vals, dtype=float).reshape(n, n)
    if data.draw(st.integers(0, 3)) == 0:
        C[np.diag_indices(n)] = data.draw(st.sampled_from((0.0, -0.0, 1e-300)))
    cmax = data.draw(st.sampled_from((None, 0.0, 1e-30, 1.0, 1e300, np.nan, np.inf)))
    tr = _screen_tr()
    if cmax is None:
        tr.__dict__.pop('_branch_cmax', None)
    else:
        tr._branch_cmax = cmax
    old = type('O', (), {'BRANCH_SCREEN_TOL': tr.BRANCH_SCREEN_TOL})()
    if cmax is not None:
        old._branch_cmax = cmax
    with warnings.catch_warnings(record=True) as W1:
        warnings.simplefilter('always')
        want = _old_proxy(old, C)
    from _pytest.monkeypatch import MonkeyPatch
    mp = MonkeyPatch()
    try:
        with warnings.catch_warnings(record=True) as W2:
            warnings.simplefilter('always')
            got = _new_screen(tr, C, mp)
    finally:
        mp.undo()
    assert (got[0], got[1] is None) == (want[0], True)
    assert np.float64(tr._branch_cmax).tobytes() == np.float64(old._branch_cmax).tobytes()
    assert sorted(str(w.message) for w in W1) == sorted(str(w.message) for w in W2)


_TR = {}


def _screen_tr():
    """A transient for the screen to run on -- constructed, never solved: the
    screen reads nothing a solve leaves (its `C` and structural rank are
    planted), and a solve here would be recorded against whichever test of
    this module ran first (the switches-off check reads it as a moved record)."""
    if 'tr' not in _TR:
        _TR['tr'] = Transient(_stage(), toolkit=circuit.numeric, integrator=Gear2Integrator())
    return _TR['tr']
