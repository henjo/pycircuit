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
    order drop after each breakpoint (served: Euler's estimate in C, 2026-10-10;
    the Python chain until then)."""
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
            ## (a finalizer's warning is not the run's: a file another test left in
            ## a reference cycle closes whenever the collector runs -- inside either
            ## run, 2026-10-04)
            warnings.simplefilter('ignore', ResourceWarning)
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
    and the Euler order drop after each edge is served too."""
    _, b = _same(_rc_pulse, monkeypatch, make={'integrator': integ()}, tend=2e-6, timestep=2e-8)
    assert b[5].get('lte_c:euler', 0) > 0 and b[5].get('lte_c:integrator', 0) == 0


_plain = st.one_of(st.floats(-1e3, 1e3), st.floats(-1e-12, 1e-12), st.just(0.0))
_extreme = st.one_of(_plain, st.sampled_from([1e-310, 1e300, -1e300, 1e-160]))


def _vec(data, n, elems=_plain):
    return np.array(data.draw(st.lists(elems, min_size=n, max_size=n)), dtype=float)


## (the profile's examples: 60 in the gate, 3000 at `--hypothesis-profile deep`)
@settings(deadline=None)
@given(data=st.data(), n=st.integers(2, 6), pi=st.booleans(),
       kind=st.sampled_from(['gear', 'trap', 'euler']), relref=st.sampled_from(sc.RELREF_MODES),
       h=st.floats(1e-15, 1e-3), r2=st.floats(0.01, 50.0), r3=st.floats(0.01, 50.0))
def test_the_c_is_the_chain_bit_for_bit(data, n, pi, kind, relref, h, r2, r3):
    """Drawn inputs -- huge and tiny charges (overflow, underflow), singular
    and ill-conditioned `J`, vector and scalar `abstol` (floored on drawn
    rows, `lte_algebraic='charge'`), every unit split:
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
    ## (Euler's: its own companion currents, and a last step that may be
    ## missing or zero -- `scale` is then 1 -- and no second one)
    iql = np.array([_vec(data, n, mag) for _ in range(3)]) if kind == 'euler' else ql.copy()
    h_last, h_last2 = hs[1], hs[2]
    if kind == 'euler':
        h_last = data.draw(st.sampled_from([h_last, h_last, None, 0.0]))
        h_last2 = data.draw(st.sampled_from([h_last2, None]))
    integ = {'gear': ig.Gear2Integrator, 'trap': ig.TrapezoidalIntegrator,
             'euler': ig.EulerIntegrator}[kind]
    kw = {'x_curr': _vec(data, n, st.floats(-5.0, 5.0)),
          'x_last': _vec(data, n, st.floats(-5.0, 5.0)),
          'q_curr': q, 'q_last_hist': ql, 'iq_last_hist': iql, 'h_curr': hs[0],
          'h_last': h_last, 'h_last2': h_last2, 'no_history': False, 'J': J,
          'active_integrator': integ(),
          'irefnode': iref, 'reltol': data.draw(st.sampled_from([1e-3, 1e-6, 0.0])),
          'abstol': abstol, 'toolkit': circuit.numeric, 'max_step': 1.0,
          'n_nodes': data.draw(st.sampled_from([None, 0, 1, n - 1, n]))}
    ## (`lte_algebraic='charge'`'s floor on drawn rows, in half the examples;
    ## chgtol / h from far below to far above `abstol`)
    if data.draw(st.booleans()):
        kw['lte_charge'] = np.array(data.draw(st.lists(st.booleans(), min_size=n, max_size=n)))
        kw['lte_chgtol'] = data.draw(st.sampled_from([1e-15, 0.0, 1e-9, 1e-3]))
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
                warnings.simplefilter('ignore', ResourceWarning)
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
            warnings.simplefilter('ignore', ResourceWarning)
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


## -- the sparse solve (the C loop plan's stage 3, step b) ------------------

def _sparse_chain():
    """A 250-stage MOS chain behind a pulse: 254 unknowns, a sparse run --
    `StructuralKLU`, its error test's `J^-1 Eg` in C."""
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['vg'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=5e-9, per=2e-8)
    for k in range(250):
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
    return c


@pytest.mark.parametrize('integ', [ig.Gear2Integrator, ig.TrapezoidalIntegrator])
def test_the_sparse_error_test_is_the_python_one(integ, monkeypatch):
    """Every attempt's verdict, error and running reference on a sparse run,
    the C solving on the run's numeric."""
    _, b = _same(_sparse_chain, monkeypatch, make={'integrator': integ(), 'reltol': 1e-5},
                 tend=1.2e-8, timestep=2e-10)
    assert b[5].get('lte_c:sparse_served', 0) > 10


def test_a_failed_residual_hands_the_solve_back(monkeypatch):
    """The refactor's validation failing on every call (the bound below
    zero): each attempt handed back -- the Python factors again -- and the
    run the Python's."""
    from pycircuit.circuit import linearsolver
    monkeypatch.setattr(linearsolver.KLUSolver, 'REFACTOR_RESIDUAL_TOL', -1.0)
    _, b = _same(_sparse_chain, monkeypatch, served_min=0, make={'reltol': 1e-5},
                 tend=6e-9, timestep=2e-10)
    assert b[5].get('lte_c:bail:residual', 0) > 5 and b[5].get('lte_c:sparse_served', 0) == 0


def _banded(m, width):
    """A banded pattern (CSC order) and its rows / cols."""
    rows, cols = [], []
    for c_ in range(m):
        for r_ in range(max(0, c_ - width), min(m, c_ + width + 1)):
            rows.append(r_)
            cols.append(c_)
    return np.array(rows), np.array(cols)


@settings(deadline=None)
@given(data=st.data(), m=st.integers(3, 9), width=st.integers(0, 2),
       kind=st.sampled_from(['gear', 'trap', 'euler']), h=st.floats(1e-12, 1e-6),
       r2=st.floats(0.1, 10.0), r3=st.floats(0.1, 10.0))
def test_the_sparse_solve_is_structural_klu_bit_for_bit(data, m, width, kind, h, r2, r3):
    """Drawn: a banded pattern and a reduced J on it (now and then a nonzero
    off it -- another matrix -- or a residual bound the refactor fails),
    drawn charges and steps: the same `(err, p)` bits, running reference,
    warnings and exception, and the solver's own counts, the C serving or
    handing back -- each side on its own solver, primed alike."""
    from pycircuit.circuit import linearsolver
    n = m + 1
    iref = data.draw(st.integers(0, n - 1))
    rows, cols = _banded(m, width)
    vals = data.draw(st.lists(st.floats(-1.0, 1.0), min_size=len(rows), max_size=len(rows)))
    Jr = np.zeros((m, m))
    Jr[rows, cols] = vals
    Jr += (2.0 + width) * np.eye(m)
    J0 = Jr.copy()
    Jr = Jr * data.draw(st.floats(0.5, 2.0))
    if data.draw(st.integers(0, 5)) == 0:
        far = [(i, j) for i in range(m) for j in range(m) if abs(i - j) > width]
        if far:
            i, j = data.draw(st.sampled_from(far))
            Jr[i, j] = 0.25
    J = np.insert(np.insert(Jr, iref, 0.0, axis=0), iref, 0.0, axis=1)
    J[iref, iref] = 1.0
    tol = data.draw(st.sampled_from([linearsolver.KLUSolver.REFACTOR_RESIDUAL_TOL] * 4 + [1e-30]))
    q = _vec(data, n)
    ql = np.array([_vec(data, n) for _ in range(3)])
    iql = np.array([_vec(data, n) for _ in range(3)])
    integ = {'gear': ig.Gear2Integrator, 'trap': ig.TrapezoidalIntegrator,
             'euler': ig.EulerIntegrator}[kind]
    out = []
    for on in (False, True):
        ls = linearsolver.AutoSolver(min_n=2, max_fill=1.0, structure=lambda: (m, rows, cols))
        ls.solve(J0, np.ones(m), circuit.numeric)          # (the run's first factor)
        sk = ls._choice
        assert type(sk) is linearsolver.StructuralKLU and sk._numeric
        sk.REFACTOR_RESIDUAL_TOL = tol
        kw = {'x_curr': _vec(data, n, st.floats(-5.0, 5.0)) if not out else out[0][-1][0],
              'x_last': _vec(data, n, st.floats(-5.0, 5.0)) if not out else out[0][-1][1],
              'q_curr': q, 'q_last_hist': ql, 'iq_last_hist': iql, 'h_curr': h,
              'h_last': h * r2, 'h_last2': h * r3, 'no_history': False, 'J': J,
              'active_integrator': integ(), 'irefnode': iref, 'reltol': 1e-3, 'abstol': 1e-9,
              'toolkit': circuit.numeric, 'max_step': 1.0, 'n_nodes': None, 'linsolver': ls}
        ctrl = sc.IntegralController()
        old = _tran_lte_c.ENABLED
        _tran_lte_c.ENABLED = on
        try:
            with warnings.catch_warnings(record=True) as W:
                warnings.simplefilter('always')
                warnings.simplefilter('ignore', ResourceWarning)
                try:
                    r, exc = ctrl._max_error(sc.StepLTEInputs(**kw)), None
                except Exception as e:                       # noqa: BLE001
                    r, exc = None, (type(e), str(e))
        finally:
            _tran_lte_c.ENABLED = old
        got = getattr(ctrl, '_ref_running', None)
        counts = (sk.analyses, sk.factors, sk.refactors, sk.residual_fallbacks, sk.foreign)
        out.append((None if r is None else (float(r[0]).hex(), r[1]), exc,
                    None if got is None else got.tobytes(), [str(w.message) for w in W], counts,
                    (kw['x_curr'], kw['x_last'])))
    assert out[0][:5] == out[1][:5]
