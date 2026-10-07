"""`elements_hdl.BSourceHdl`: `elements.BSource` in the HDL -- its callables
traced into a class compiled per expression.  Its stamps are `BSource`'s:
the same values to rounding (the expression is sympy's, not the
callable's order of operations), the same terminals and orientation, the
Jacobian exact where `BSource` takes a central difference.  A circuit of
such elements and constant stamps runs in C, where `BSource`'s cannot;
the van der Pol PSS agrees with `BSource`'s.  One class per expression,
its compile cached by the expression; what does not trace is refused,
saying what does."""
import warnings

import numpy as np
import pytest
import sympy

from pycircuit.circuit import _hdl_cache, _paths, _tran_core, circuit
from pycircuit.circuit.elements import IS, BSource, L, gnd
from pycircuit.circuit.elements import C as Cap
from pycircuit.circuit.elements_hdl import BSourceHdl
from pycircuit.circuit.integrator import Gear2Integrator
from pycircuit.circuit.transient import Transient

EPS = np.finfo(float).eps
MU = 1.0


def _cubic(u):
    return MU * (u - u ** 3 / 3.0)


def _charge(u):
    return 1e-9 * u + 2e-10 * u ** 2


VOLTAGES = (0.0, -0.0, 1e-3, -1e-3, 0.3, -1.0, 1.7320508075688772, -2.0, 10.0)


def _pair(**kw):
    circuit.default_toolkit = circuit.numeric
    return (BSource('a', 'b', 'c', 'd', **kw), BSourceHdl('a', 'b', 'c', 'd', **kw))


def test_the_currents_are_bsources_to_rounding_and_the_jacobian_exact():
    """`i` as `BSource`'s, to a few ulps of its terms; `G` the exact
    derivative, which `BSource`'s central difference (step 1e-6) meets to
    its own error; the terminals and their signs the same."""
    ref, hdl = _pair(i_func=_cubic)
    for v in VOLTAGES:
        for x in (np.array([v, 0.0, 0.2, -0.1]), np.array([0.5, 0.5 - v, 0.0, 0.0])):
            u = x[0] - x[1]
            a, b = hdl.i(x), ref.i(x)
            scale = abs(u) + abs(u) ** 3 / 3.0
            assert a[:2].tolist() == [0.0, 0.0] and a[3] == -a[2]
            assert abs(a[2] - b[2]) <= 8 * EPS * scale, (v, a, b)
            G, Gr = hdl.G(x), ref.G(x)
            exact = MU * (1.0 - u * u)
            assert abs(G[2, 0] - exact) <= 8 * EPS * (1.0 + u * u), (v, G[2, 0], exact)
            assert np.array_equal(G[2], [G[2, 0], -G[2, 0], 0.0, 0.0])
            assert np.array_equal(G[3], -G[2]) and not G[:2].any()
            assert abs(Gr[2, 0] - G[2, 0]) <= 1e-9 * (1.0 + abs(b[2])), (v, Gr[2, 0], G[2, 0])
            assert not hdl.q(x).any() and not hdl.C(x).any()


def test_the_charges_are_bsources_and_their_capacitance_exact():
    ref, hdl = _pair(q_func=_charge)
    for v in VOLTAGES:
        x = np.array([v, 0.0, 0.0, 0.0])
        a, b = hdl.q(x), ref.q(x)
        assert abs(a[2] - b[2]) <= 8 * EPS * (1e-9 * abs(v) + 2e-10 * v * v)
        assert a[3] == -a[2] and not a[:2].any()
        C = hdl.C(x)
        assert abs(C[2, 0] - (1e-9 + 4e-10 * v)) <= 8 * EPS * (1e-9 + 4e-10 * abs(v))
        assert abs(ref.C(x)[2, 0] - C[2, 0]) <= 1e-9 * (1e-9 + abs(b[2]))
        assert not hdl.i(x).any() and not hdl.G(x).any()


def test_sympy_functions_trace():
    """An exponential through `sympy.exp`: the value numpy's, the
    derivative exact."""
    circuit.default_toolkit = circuit.numeric
    hdl = BSourceHdl('a', 'b', 'c', 'd', i_func=lambda u: 1e-14 * (sympy.exp(u / 0.025) - 1))
    for v in (0.0, 0.3, 0.6, -0.5):
        x = np.array([v, 0.0, 0.0, 0.0])
        assert hdl.i(x)[2] == pytest.approx(1e-14 * np.expm1(v / 0.025), rel=1e-13, abs=1e-30)
        assert hdl.G(x)[2, 0] == pytest.approx(1e-14 / 0.025 * np.exp(v / 0.025), rel=1e-13)


def test_no_function_is_a_zero_element():
    circuit.default_toolkit = circuit.numeric
    z = BSourceHdl('a', 'b', 'c', 'd')
    x = np.array([0.4, -0.2, 0.1, 0.0])
    assert not z.i(x).any() and not z.G(x).any() and not z.q(x).any() and not z.C(x).any()


def test_one_class_per_expression_its_compile_keyed_by_it():
    """The same expression -- however written -- is the same class; a value
    a callable closes over is part of its expression; the compile cache
    keys each class by its expressions, so none is served another's."""
    circuit.default_toolkit = circuit.numeric
    a = BSourceHdl('a', 'b', 'c', 'd', i_func=_cubic)
    b = BSourceHdl('x', 'y', 'z', 'w', i_func=lambda v: MU * (v - v * v * v / 3.0))
    k = 2.0
    c = BSourceHdl('a', 'b', 'c', 'd', i_func=lambda u: k * (u - u ** 3 / 3.0))
    assert type(a) is type(b) and type(c) is not type(a)
    assert type(a).bsource_exprs[0] != type(c).bsource_exprs[0]
    assert not type(a)._hdl_cache_status.startswith('uncacheable'), type(a)._hdl_cache_status
    assert _hdl_cache.key_for(type(a)) != _hdl_cache.key_for(type(c))


@pytest.mark.parametrize('func, says', [
    (lambda u: np.tanh(u), 'does not trace'),
    (lambda u: u if u > 0 else 0.0, 'does not trace'),
    (lambda u: u * sympy.Symbol('k'), 'depends on k'),
    (lambda u: 1j * u, 'complex'),
], ids=['numpy', 'if', 'free-symbol', 'complex'])
def test_what_does_not_trace_is_refused_saying_what_does(func, says):
    circuit.default_toolkit = circuit.numeric
    with pytest.raises(ValueError, match=says) as e:
        BSourceHdl('a', 'b', 'c', 'd', i_func=func)
    if says == 'does not trace':
        assert 'elements.BSource' in str(e.value) and 'sympy' in str(e.value)
    with pytest.raises(TypeError, match='callable'):
        BSourceHdl('a', 'b', 'c', 'd', q_func=3.0)


def _vdp(hdl, psd=1e-6):
    """The van der Pol oscillator of `tests/_shooting_fixtures`, its
    nonlinearity either element."""
    circuit.default_toolkit = circuit.numeric
    c = circuit.SubCircuit()
    c.add_node('v')
    c['C'] = Cap('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = (BSourceHdl if hdl else BSource)('v', gnd, gnd, 'v', i_func=_cubic)
    c['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    return c, x0


def test_a_circuit_of_them_runs_in_c():
    """Served: with `BSourceHdl` the evaluate core takes the van der Pol
    circuit and the C Newton its steps; with `BSource` the core refuses
    it (a Python element) and every step is the Python Newton's."""
    _tran_core.driver()
    if _tran_core.STATUS != 'c':
        pytest.skip(f'the core is off: {_tran_core.STATUS}')
    got = {}
    for hdl in (False, True):
        cir, _x0 = _vdp(hdl)
        tr = Transient(cir, integrator=Gear2Integrator())
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr.solve(tend=2.0, timestep=0.05, fixed_timestep=True)
        got[hdl] = _paths.since(before)
    assert got[True].get('newton_c:served', 0) > 30, got[True]
    assert not got[True].get('newton_c:unservable'), got[True]
    assert got[False].get('newton_c:unservable', 0) > 30, got[False]
    assert not got[False].get('newton_c:served'), got[False]


def test_the_van_der_pol_pss_is_bsources():
    """The limit cycle, its period and its monodromy as `BSource`'s: the
    waveforms to rounding, the monodromy to the central difference's
    error in `BSource`'s Jacobian (measured 2.7e-14 and 1.1e-10)."""
    from pycircuit.circuit import PSS
    out = {}
    for hdl in (False, True):
        cir, x0 = _vdp(hdl)
        p = PSS(cir, method='gear', reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p.solve(period=6.6634, timestep=6.6634 / 200, x0=x0, maxiterations=60)
        out[hdl] = (np.asarray(p.waveform[1], float), np.asarray(p._monodromy, float),
                    float(p.period))
    (wa, ma, ta), (wb, mb, tb) = out[False], out[True]
    assert wa.shape == wb.shape and np.max(np.abs(wa - wb)) <= 1e-12 * np.max(np.abs(wa))
    assert abs(ta - tb) <= 1e-12 * ta
    assert np.max(np.abs(ma - mb)) <= 1e-8
    ea = np.sort(np.abs(np.linalg.eigvals(ma)))
    eb = np.sort(np.abs(np.linalg.eigvals(mb)))
    assert np.allclose(ea[-2:], eb[-2:], rtol=1e-8, atol=1e-12)
