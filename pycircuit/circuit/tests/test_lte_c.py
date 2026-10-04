"""The adaptive error test in C (`_tran_lte_c`, speed round 8, stage 3):
bit for bit `_charge_lte`, `_normalised` and their maximum, the running
reference with them; its declines before anything is touched, its bails
without a trace."""
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _paths, _tran_lte_c, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import integrator as ig
from pycircuit.circuit import stepcontroller as sc
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, VSin, gnd
from pycircuit.circuit.transient import Transient


def _mos_chain(n=4, va=0.9):
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vg'] = VSin('g0', gnd, v=0.9, va=va, freq=1e6)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', f'g{k}' if k == 0 else f'd{k - 1}', gnd, gnd)
    return c


def _rc_pulse():
    """Two RC sections behind a pulse: rejections at the edges, the Euler
    order drop after each breakpoint (declined: the Python chain)."""
    c = SubCircuit()
    c['vs'] = VPulse('in', gnd, v1=0.0, v2=1.0, td=1e-7, tr=1e-9, tf=1e-9, pw=4e-7, per=1e-6)
    c['R'] = R('in', 'out', r=1e3)
    c['C'] = C('out', gnd, c=1e-10)
    c['R2'] = R('out', 'o2', r=2e3)
    c['C2'] = C('o2', gnd, c=5e-11)
    return c


def _run(build, on, monkeypatch, pi=False, **kw):
    """`(x bytes, time bytes, statistics, warnings, per-attempt record,
    path counts)` of one adaptive transient with the C on or off."""
    monkeypatch.setattr(_tran_lte_c, 'ENABLED', on)
    log = []
    real = sc.StepController.evaluate_step

    def recording(self, *a, **k):
        r = real(self, *a, **k)
        run = getattr(self, '_ref_running', None)
        log.append((r, self.last_err, None if run is None else run.tobytes()))
        return r
    monkeypatch.setattr(sc.StepController, 'evaluate_step', recording)
    before = _paths.snapshot()
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            tr = Transient(build(), toolkit=circuit.numeric, **kw.pop('make', {}))
            if pi:
                tr.step_controller = sc.PIController()
            res = tr.solve(**kw)
    finally:
        monkeypatch.setattr(sc.StepController, 'evaluate_step', real)
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__
           if 'seconds' not in k} if hasattr(tr.statistics, '__slots__') else \
        {k: v for k, v in vars(tr.statistics).items() if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), log, _paths.since(before))


def _same(build, monkeypatch, served_min=1, **kw):
    a = _run(build, False, monkeypatch, **dict(kw))
    b = _run(build, True, monkeypatch, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert a[4] == b[4], 'an attempt was judged differently'
    assert b[5].get('lte_c:served', 0) >= served_min, b[5]
    return a, b


@pytest.mark.parametrize('relref', ['sigglobal', 'alllocal', 'pointlocal'])
@pytest.mark.parametrize('pi', [False, True])
def test_the_error_test_is_the_python_one(relref, pi, monkeypatch):
    """Every attempt's verdict, next step, error and running reference, for
    both controllers and every `relref`."""
    _same(_mos_chain, monkeypatch, pi=pi, make={'relref': relref}, tend=1.5e-6, timestep=2e-8)


@pytest.mark.parametrize('integ', [ig.Gear2Integrator, ig.TrapezoidalIntegrator])
def test_rejections_and_order_drops(integ, monkeypatch):
    """A pulse: rejected attempts fold into the running reference as before,
    and the Euler order drop after each edge takes the Python chain."""
    _, b = _same(_rc_pulse, monkeypatch, make={'integrator': integ()}, tend=2e-6, timestep=2e-8)
    assert b[5].get('lte_c:integrator', 0) > 0


_plain = st.one_of(st.floats(-1e3, 1e3), st.floats(-1e-12, 1e-12), st.just(0.0))
_extreme = st.one_of(_plain, st.sampled_from([1e-310, 1e300, -1e300, 1e-160]))


def _vec(data, n, elems=_plain):
    return np.array(data.draw(st.lists(elems, min_size=n, max_size=n)), dtype=float)


## (the profile's examples: 60 in the gate, 3000 at `--hypothesis-profile deep`)
@settings(deadline=None)
@given(data=st.data(), n=st.integers(2, 6), pi=st.booleans(),
       trap=st.booleans(), relref=st.sampled_from(sc.RELREF_MODES),
       h=st.floats(1e-15, 1e-3), r2=st.floats(0.01, 50.0), r3=st.floats(0.01, 50.0))
def test_the_c_is_the_chain_bit_for_bit(data, n, pi, trap, relref, h, r2, r3):
    """Drawn inputs -- huge and tiny charges (overflow, underflow), singular
    and ill-conditioned `J`, vector and scalar `abstol`, every unit split:
    the same `(err, p)` bits, running reference, warnings and exception, the
    C served or handed back.  The steps Python floats or numpy's (a period's
    grid gives numpy's)."""
    iref = data.draw(st.integers(0, n - 1))
    hs = (h, h * r2, h * r3)
    if data.draw(st.booleans()):
        hs = tuple(np.float64(v) for v in hs)
    ## (one example in five with magnitudes that overflow or underflow, one
    ## in five with a zero row in J: singular, or the reference row's)
    mag = _extreme if data.draw(st.integers(0, 4)) == 0 else _plain
    q = _vec(data, n, mag)
    ql = np.array([_vec(data, n, mag) for _ in range(3)])
    ## (a diagonal shift: drawn entries are often equal or zero, and the
    ## solve should be served more often than handed back)
    J = np.array([_vec(data, n, st.floats(-10.0, 10.0)) for _ in range(n)]) + 25.0 * np.eye(n)
    if data.draw(st.integers(0, 4)) == 0:
        J[data.draw(st.integers(0, n - 1))] = 0.0
    abstol = (_vec(data, n, st.floats(0.0, 1e-3)) if data.draw(st.booleans())
              else data.draw(st.floats(0.0, 1e-3)))
    run = (np.abs(_vec(data, n, st.floats(-2.0, 2.0)))
           if relref != 'pointlocal' and data.draw(st.booleans()) else None)
    kw = {'x_curr': _vec(data, n, st.floats(-5.0, 5.0)),
          'x_last': _vec(data, n, st.floats(-5.0, 5.0)),
          'q_curr': q, 'q_last_hist': ql, 'iq_last_hist': ql.copy(), 'h_curr': hs[0],
          'h_last': hs[1], 'h_last2': hs[2], 'no_history': False, 'J': J,
          'active_integrator': ig.TrapezoidalIntegrator() if trap else ig.Gear2Integrator(),
          'irefnode': iref, 'reltol': data.draw(st.sampled_from([1e-3, 1e-6, 0.0])),
          'abstol': abstol, 'toolkit': circuit.numeric, 'max_step': 1.0,
          'n_nodes': data.draw(st.sampled_from([None, 0, 1, n - 1, n]))}
    out = []
    for on in (False, True):
        ctrl = sc.PIController() if pi else sc.IntegralController()
        ctrl.set_relref(relref)
        if run is not None:
            ctrl._ref_running = run.copy()
        old = _tran_lte_c.ENABLED
        _tran_lte_c.ENABLED = on
        try:
            with warnings.catch_warnings(record=True) as W:
                warnings.simplefilter('always')
                try:
                    r, exc = ctrl._max_error(sc.StepLTEInputs(**kw)), None
                except Exception as e:                       # noqa: BLE001
                    r, exc = None, (type(e), str(e))
        finally:
            _tran_lte_c.ENABLED = old
        got = getattr(ctrl, '_ref_running', None)
        out.append((None if r is None else (float(r[0]).hex(), r[1]), exc,
                    None if got is None else got.tobytes(), [str(w.message) for w in W]))
    assert out[0] == out[1]


def test_a_singular_jacobian_bails_and_warns_as_before(monkeypatch):
    """The LTE solve fails: the C hands back (`dgesv`), the Python chain
    warns once and uses the charge residual -- the same `(err, p)`."""
    n = 4
    kw = {'x_curr': np.ones(n), 'x_last': np.ones(n), 'q_curr': np.arange(n, dtype=float),
          'q_last_hist': np.array([np.arange(n) * 0.5, np.arange(n) * 0.2, np.zeros(n)]),
          'iq_last_hist': np.zeros((3, n)), 'h_curr': 1e-9, 'h_last': 1e-9, 'h_last2': 1e-9,
          'no_history': False, 'J': np.zeros((n, n)), 'active_integrator': ig.Gear2Integrator(),
          'irefnode': 0, 'reltol': 1e-3, 'abstol': 1e-6, 'toolkit': circuit.numeric,
          'max_step': 1.0}
    res = []
    for on in (False, True):
        monkeypatch.setattr(_tran_lte_c, 'ENABLED', on)
        before = _paths.snapshot()
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            r = sc.IntegralController()._max_error(sc.StepLTEInputs(**kw))
        res.append((r, [str(w.message) for w in W], _paths.since(before)))
    assert res[0][:2] == res[1][:2] and len(res[1][1]) == 1
    assert res[1][2].get('lte_c:bail:dgesv', 0) == 1


def test_the_declines(monkeypatch):
    """A patched piece of the chain, an instance shadow, a controller
    subclass: the Python chain, counted, the same verdicts."""
    real = ig.third_divided_difference
    monkeypatch.setattr(ig, 'third_divided_difference', lambda *a: real(*a))
    _, b = _same(_mos_chain, monkeypatch, served_min=0, tend=6e-7, timestep=2e-8)
    assert b[5].get('lte_c:patched', 0) > 0 and b[5].get('lte_c:served', 0) == 0
    monkeypatch.undo()

    class Mine(sc.IntegralController):
        pass
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True):
        warnings.simplefilter('always')
        tr = Transient(_mos_chain(), toolkit=circuit.numeric)
        tr.step_controller = Mine()
        tr.solve(tend=6e-7, timestep=2e-8)
        d = _paths.since(before)
        calls = []
        tr = Transient(_mos_chain(), toolkit=circuit.numeric)
        tr.step_controller = ctrl = sc.IntegralController()
        ctrl._charge_lte = lambda s: (calls.append(1), sc.StepController._charge_lte(ctrl, s))[1]
        tr.solve(tend=6e-7, timestep=2e-8)
    assert d.get('lte_c:controller', 0) > 0 and d.get('lte_c:served', 0) == 0
    assert len(calls) >= 5


def test_the_c_serves_the_adaptive_step():
    if _tran_lte_c.driver() is None:
        pytest.skip(f'the C error test is off: {_tran_lte_c.STATUS}')
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True):
        warnings.simplefilter('always')
        Transient(_mos_chain(), toolkit=circuit.numeric).solve(tend=1e-6, timestep=2e-8)
    assert _paths.since(before).get('lte_c:served', 0) >= 20
