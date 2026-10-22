"""The dense multistep walk's monodromy step in C (`shooting/_lmmsens_c`):
bit for bit the Python step (`_LiveMonodromy.step` -- `_lmm_recursion`
through `_ColumnsOf` and numpy's solve), served on the shooting Newton's
walks, declined -- counted -- wherever it does not stand in for that step,
and handing back where numpy's solve raises or the result is not finite."""
import itertools
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, circuit, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import C, Diode, R, VSin
from pycircuit.circuit.linearsolver import DenseSolver
from pycircuit.circuit.shooting import PSS, _lmmsens_c, _steps
from pycircuit.circuit.shooting._steps import _ColumnsOf, _lmm_recursion
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._shooting_fixtures import _q20_rlc
from pycircuit.circuit.tests._warnpolicy import quiet

pytestmark = pytest.mark.skipif(_lmmsens_c.driver() is None,
                                reason='no C build: ' + _lmmsens_c.STATUS)


def _arrays(m, w, nks, seed):
    """`(P0, P1, Pq, C1, C2, Cn, Jf)`: random monodromies and capacitances
    with `nks` nonzero columns each (a -0.0 elsewhere in one, a zero to the
    column scan)."""
    rng = np.random.default_rng(seed)
    P0, P1, Pq = (rng.standard_normal((m, w)) for _ in range(3))
    Cs = []
    for nk in nks:
        Cm = np.zeros((m, m))
        cols = rng.choice(m, size=nk, replace=False)
        Cm[:, cols] = rng.standard_normal((m, nk)) * (rng.random((m, nk)) < 0.6)
        Cm[:, cols[:1]] += 1.0
        Cm[0, [c for c in range(m) if c not in cols][:1]] = -0.0
        Cs.append(np.ascontiguousarray(Cm))
    Jf = rng.standard_normal((m, m)) + 4 * np.eye(m)
    return P0, P1, Pq, Cs[0], Cs[1], Cs[2], Jf


def _python(P0, P1, Pq, C1, C2, Cn, Jf, alphas, b):
    return _lmm_recursion([P0, P1], [_ColumnsOf(C1), _ColumnsOf(C2)], Pq, _ColumnsOf(Cn),
                          alphas, b, lambda S: DenseSolver().solve(Jf, S, numeric))


CASES = list(itertools.product([(0, 0, 0), (1, 1, 1), (2, 1, 3), (5, 7, 2), (9, 9, 9)],
                               [0.0, 0.25], [2, 3]))


@pytest.mark.parametrize(('nks', 'b', 'L'), CASES,
                         ids=[f"nk{'_'.join(map(str, n))}-b{b}-L{L}" for n, b, L in CASES])
def test_the_c_step_is_the_python_step_bit_for_bit(nks, b, L):
    """Every product branch -- no column (zeros), one (numpy's own loop),
    several and all (its `cblas_dgemm`) -- with and without the companion's
    `b Pq`, a two- and a three-term companion."""
    m, w = 9, 5
    arrs = _arrays(m, w, nks, seed=sum(nks) + 7 * L)
    alphas = (1.5, -2.0, 0.5)[:L]
    want = _python(*arrs, alphas, b)
    P0, P1, Pq, C1, C2, Cn, Jf = arrs
    got = _lmmsens_c._call({}, _lmmsens_c.driver(), m, w, alphas, b, P0, P1, Pq, [C1, C2], Cn, Jf)
    assert type(got) is tuple
    for g, r in zip(got, want):
        assert g.tobytes() == np.ascontiguousarray(r).tobytes()


def test_kept_column_lists_and_pointers_answer_as_fresh_ones():
    """The ring: each call's inputs the last call's arrays moved on one
    place, their column lists and pointers kept -- the answer a fresh
    context's."""
    m, w = 9, 4
    P0, P1, Pq, C1, C2, _Cn, Jf = _arrays(m, w, (3, 4, 5), seed=3)
    rng = np.random.default_rng(4)
    alphas, wd, drv = (1.5, -2.0, 0.5), {}, _lmmsens_c.driver()
    Cs = [C1, C2]
    for _ in range(5):
        Cn = np.ascontiguousarray(Cs[0] * (1 + 0.1 * rng.standard_normal((m, m))))
        got = _lmmsens_c._call(wd, drv, m, w, alphas, 0.0, P0, P1, Pq, Cs, Cn, Jf)
        fresh = _lmmsens_c._call({}, drv, m, w, alphas, 0.0, P0, P1, Pq, Cs, Cn, Jf)
        for g, f in zip(got, fresh):
            assert g.tobytes() == f.tobytes()
        P0, P1, Pq, Cs = got[0], P0, got[1], [Cn, Cs[0]]
    assert wd['_lmmsens_c'].t.reuse == 3


def test_a_singular_jacobian_and_a_non_finite_result_are_handed_back():
    m, w = 6, 3
    P0, P1, Pq, C1, C2, Cn, Jf = _arrays(m, w, (2, 2, 2), seed=5)
    drv, alphas = _lmmsens_c.driver(), (1.5, -2.0, 0.5)
    Js = Jf.copy()
    Js[:, 2] = 0.0
    assert _lmmsens_c.BAIL[_lmmsens_c._call({}, drv, m, w, alphas, 0.0, P0, P1, Pq, [C1, C2],
                                            Cn, Js)] == 'singular'
    with pytest.raises(np.linalg.LinAlgError):
        _python(P0, P1, Pq, C1, C2, Cn, Js, alphas, 0.0)
    ## (an inf the products read: in a nonzero column of a capacitance --
    ## one in a row of P no product reads stays out of the result, the
    ## Python step's too)
    Cinf = C1.copy()
    Cinf[0, np.flatnonzero(np.any(C1 != 0, axis=0))[0]] = np.inf
    got = _lmmsens_c._call({}, drv, m, w, alphas, 0.0, P0, P1, Pq, [Cinf, C2], Cn, Jf)
    assert _lmmsens_c.BAIL[got] == 'nonfinite'
    with np.errstate(invalid='ignore', over='ignore'):
        want = _python(P0, P1, Pq, Cinf, C2, Cn, Jf, alphas, 0.0)
    assert not np.isfinite(want[0]).all()


def _clipper(two_caps=True):
    c = SubCircuit()
    c.add_nodes('a', 'b', 'd')
    c['vs'] = VSin('a', gnd, va=2.0, freq=1e3)
    c['R1'] = R('a', 'b', r=1e3)
    c['D'] = Diode('b', gnd)
    c['C1'] = C('b', gnd, c=1e-7)
    c['R2'] = R('b', 'd', r=2e3)
    c['R3'] = R('d', gnd, r=2e3)
    if two_caps:
        c['C2'] = C('d', gnd, c=1e-8)
    return c


def _walk(make, method, monkeypatch, on, **kw):
    circuit.default_toolkit = circuit.numeric
    monkeypatch.setattr(_lmmsens_c, 'LMMSENS_C', on)
    cir = make()
    m = cir.n - 1
    pss = PSS(cir, method=method, reltol=1e-9, **kw)
    pss._open_at_x0 = False
    pss.autonomous = False
    kind = pss._map_kind()
    width = 2 if kind == 'pair' else 1
    times, hs = pss._period_grid(1e-3, 60, None)
    before = _paths.snapshot()
    with quiet():
        M = pss._walk(kind, np.full(width * m, 0.1), times, hs, T=1e-3).monodromy()
    return M, _paths.since(before)


MAKERS = [(_q20_rlc, 'rlc'), (_clipper, 'clipper'), (lambda: _clipper(False), 'one cap')]


@pytest.mark.parametrize('method', ['gear2', 'trap', 'euler'])
@pytest.mark.parametrize('make', [m for m, _ in MAKERS], ids=[i for _, i in MAKERS])
def test_a_walk_s_monodromy_is_the_python_step_s(make, method, monkeypatch):
    M0, d0 = _walk(make, method, monkeypatch, False)
    M1, d1 = _walk(make, method, monkeypatch, True)
    assert M1.tobytes() == M0.tobytes()
    assert d0.get('pss.lmm:off', 0) > 0 and not d0.get('pss.lmm:served')
    ## (the one-capacitor circuit's plain map carries one live column: the
    ## Python step's)
    if make is _clipper or method == 'gear2':
        assert d1.get('pss.lmm:served', 0) > 0


def test_a_pss_on_the_c_step_is_the_python_step_s(monkeypatch):
    circuit.default_toolkit = circuit.numeric
    out = {}
    for on in (False, True):
        monkeypatch.setattr(_lmmsens_c, 'LMMSENS_C', on)
        before = _paths.snapshot()
        with quiet(), warnings.catch_warnings():
            warnings.simplefilter('ignore', AccuracyWarning)
            res = PSS(_clipper(), method='gear2').solve(period=1e-3, timestep=2e-6)
        out[on] = (np.asarray(res['tpss'].v('d').get_y(), dtype=float), _paths.since(before))
    assert out[True][0].tobytes() == out[False][0].tobytes()
    assert out[True][1].get('pss.lmm:served', 0) > 0


def test_it_declines_where_it_does_not_stand_in_for_the_python_step(monkeypatch):
    """Every column carried (`LIVE_COLUMNS` off), a caller's own solver, a
    patched piece of the Python step: declined and counted, the Python step
    run."""
    monkeypatch.setattr(_steps, 'LIVE_COLUMNS', False)
    _M, d = _walk(_clipper, 'gear2', monkeypatch, True)
    assert d.get('pss.lmm:every_column', 0) > 0 and not d.get('pss.lmm:served')
    monkeypatch.setattr(_steps, 'LIVE_COLUMNS', True)

    class Mine(DenseSolver):
        pass
    _M, d = _walk(_clipper, 'gear2', monkeypatch, True, linearsolver=Mine())
    assert d.get('pss.lmm:solver', 0) > 0 and not d.get('pss.lmm:served')
    M0, _d = _walk(_clipper, 'gear2', monkeypatch, False)
    orig = _steps._lmm_recursion

    def patched(*a, **k):
        return orig(*a, **k)
    monkeypatch.setattr(_steps, '_lmm_recursion', patched)
    M1, d = _walk(_clipper, 'gear2', monkeypatch, True)
    assert d.get('pss.lmm:patched', 0) > 0 and not d.get('pss.lmm:served')
    assert M1.tobytes() == M0.tobytes()


def test_a_kept_walk_takes_its_own_factors(monkeypatch):
    """A walk keeping its steps solves through each step's factorisation,
    not numpy's solve: the C is not asked."""
    circuit.default_toolkit = circuit.numeric
    monkeypatch.setattr(_lmmsens_c, 'LMMSENS_C', True)
    cir = _clipper()
    pss = PSS(cir, method='gear2', reltol=1e-9)
    pss._open_at_x0 = False
    pss.autonomous = False
    times, hs = pss._period_grid(1e-3, 30, None)
    before = _paths.snapshot()
    with quiet():
        pss._walk('pair', np.full(2 * (cir.n - 1), 0.1), times, hs, T=1e-3, keep=True)
    assert not any(k.startswith('pss.lmm:') for k in _paths.since(before))
