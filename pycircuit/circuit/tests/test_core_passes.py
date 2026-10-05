"""The evaluate core's passes for the stage methods (`_tran_core.passes`,
speed round 9, stage 8): the circuit's passes G, C, i, q at a point from ONE
core call (formula 4: the passes alone, through the fused kernels) -- the
bytes of `cir.G/C/i/q`, fresh arrays -- and the stage paths that take them
(Radau's full and transformed Newton, its error estimate, the DIRK stage
Newton, the stage step's finish) with the answers, statistics and warnings
of their own calls."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_core, circuit
from pycircuit.circuit.elements import Diode, gnd
from pycircuit.circuit.integrator import (
    ESDIRK43Integrator,
    RadauIIA3Integrator,
    TRBDF2Integrator,
)
from pycircuit.circuit.tests.test_newton_c import _gp_chain, _mos_chain
from pycircuit.circuit.tests.test_psp_limit_c import _stage
from pycircuit.circuit.transient import Transient

SPECIAL = (0.0, -0.0, 1e-12, -1e-12, 0.7, -0.7, 3.0, -3.0, 40.0, -40.0)


def _tr(build):
    tr = Transient(build(), toolkit=circuit.numeric)
    tr.solve(tend=1e-7, timestep=2e-8, fixed_timestep=True)
    if _tran_core.core_for(tr) is None:
        pytest.skip('the core does not serve this circuit here')
    return tr


@pytest.mark.parametrize('build', [_mos_chain, _gp_chain, _stage], ids=['mos', 'gp', 'psp'])
def test_the_passes_are_the_circuits(build):
    """Every subset of the passes, at drawn and special states: the bytes of
    the circuit's own passes, and fresh arrays."""
    tr = _tr(build)
    n = tr.cir.n
    rng = np.random.default_rng(1)
    for k in range(60):
        x = rng.standard_normal(n) * (0.3 if k % 2 else 2.0)
        if k % 3 == 0:
            x[rng.integers(0, n)] = rng.choice(SPECIAL)
        for which in ('GCiq', 'qiGC', 'q', 'qi', 'CG', 'CGi', 'iqCG', 'G'):
            before = _paths.snapshot()
            P = _tran_core.passes(tr, x, which)
            assert P is not None and _paths.since(before).get('core.p:served') == 1
            assert list(P) == list(which)
            for m in which:
                ref = np.asarray(getattr(tr.cir, m)(x, tr.epar), dtype=float)
                assert P[m].tobytes() == ref.tobytes(), (m, which, x.tolist())
        P1 = _tran_core.passes(tr, x, 'GCiq')
        P2 = _tran_core.passes(tr, x, 'GCiq')
        assert all(P1[m] is not P2[m] for m in 'GCiq')


def test_what_the_passes_decline(monkeypatch):
    tr = _tr(_mos_chain)
    x = np.zeros(tr.cir.n)

    def declined(why, xx=x):
        before = _paths.snapshot()
        assert _tran_core.passes(tr, xx, 'qi') is None
        assert _paths.since(before).get('core.p:' + why) == 1
    monkeypatch.setattr(_tran_core, 'CORE_PASSES', False)
    declined('off')
    monkeypatch.undo()
    bad = x.copy()
    bad[1] = np.nan
    declined('x', bad)
    declined('x', x.astype(np.float32))
    declined('x', x[:-1])
    tr.cir.__dict__['q'] = lambda *a, **k: None
    declined('shadow_cir')
    del tr.cir.__dict__['q']
    ## a hand-written diode: the core cannot serve the circuit
    def with_diode():
        c = _mos_chain(2)
        c['D'] = Diode('d0', gnd)
        return c
    tr2 = Transient(with_diode(), toolkit=circuit.numeric)
    tr2.solve(tend=4e-8, timestep=2e-8, fixed_timestep=True)
    before = _paths.snapshot()
    assert _tran_core.passes(tr2, np.zeros(tr2.cir.n), 'q') is None
    assert _paths.since(before).get('core.p:unservable') == 1


def _run(build, on, monkeypatch, integ, **kw):
    monkeypatch.setattr(_tran_core, 'CORE_PASSES', on)
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        warnings.simplefilter('ignore', ResourceWarning)
        tr = Transient(build(), toolkit=circuit.numeric, integrator=integ(),
                       **kw.pop('make', {}))
        before = _paths.snapshot()
        res = tr.solve(**kw)
        d = _paths.since(before)
    monkeypatch.undo()
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), d)


CASES = [
    ('radau-mos', _mos_chain, RadauIIA3Integrator, {}),
    ('radau-gp-transform', _gp_chain, RadauIIA3Integrator, {}),
    ('radau-psp', _stage, RadauIIA3Integrator, {}),
    ('radau-mos-transform', _mos_chain, RadauIIA3Integrator, {'make': {'radau_transform': True}}),
    ('trbdf2-mos', _mos_chain, TRBDF2Integrator, {}),
    ('esdirk-gp', _gp_chain, ESDIRK43Integrator, {}),
]


@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
@pytest.mark.parametrize('name, build, integ, extra', CASES, ids=[c[0] for c in CASES])
def test_a_stage_method_is_the_same_with_the_cores_passes(name, build, integ, extra, fixed,
                                                          monkeypatch):
    kw = dict(extra, tend=6e-7, timestep=2e-8)
    if fixed:
        kw['fixed_timestep'] = True
    a = _run(build, False, monkeypatch, integ, **dict(kw))
    b = _run(build, True, monkeypatch, integ, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert b[4].get('core.p:served', 0) > 0, b[4]
