"""The sources evaluated directly (`_stamp_plan.source_direct`, speed round
10, B3.7) and `Sin.f` at a scalar time (`func._sin_scalar`).  A `u` pass
whose called elements are all independent sources of the VS or IS family is
each source's value -- its parameter plus its time function at `t`, the
expression its `u` makes -- scattered as the loop's bincount scatters, in
Python floats, its checks stamped (`_watch`).  Bit for bit the loop (bytes,
dtype, warnings with their lines, the `u:called` counts) on drawn circuits,
values, times and analyses, with signed zeros and NaN payloads enumerated;
served where it should be, and the loop's answer wherever a check fails --
a method patched or shadowed, a parameter or a time function replaced, an
element added or its class changed -- with its stamp holding where a
limiter writes its element's dict every iteration.  `Sin.f`'s scalar form
is its expression's bits and type, and the expression wherever it could
warn."""
import math
import struct
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _paths, _stamp_plan, circuit, func
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import (
    IPWL,
    IS,
    VPWL,
    VS,
    C,
    Diode,
    IExp,
    IPulse,
    ISin,
    R,
    SubCircuit,
    VExp,
    VPulse,
    VSin,
    gnd,
)
from pycircuit.circuit.integrator import Gear2Integrator, RadauIIA3Integrator
from pycircuit.circuit.transient import Transient


def _nan(payload, sign=0):
    return struct.unpack('<d', struct.pack('<Q', (sign << 63) | (0x7FF8 << 48) | payload))[0]


#: values a parameter takes: ordinary, signed zeros, NaN of three payloads,
#: infinities, subnormals, the range's ends, integers
SPECIAL = (0.0, -0.0, _nan(0), _nan(0x123), _nan(0x456, 1), math.inf, -math.inf, 5e-324,
           -5e-324, 1e300, -1e300, 1.7e308, -1.7e308, 1, -2, 0)
ORDINARY = (0.7, -1.2, 2.5e-3, 1.0, -0.35)
NODES = ('a', 'b', 'c')
ANALYSES = ('tran', 'dc', None, 'noise')


def _source(kind, plus, minus, v, w):
    """A source of the family by `kind`, its value `v` (and an amplitude or
    second level `w`)."""
    if kind == 'VS':
        return VS(plus, minus, v=v)
    if kind == 'IS':
        return IS(plus, minus, i=v)
    if kind == 'VSin':
        return VSin(plus, minus, v=v, va=w, freq=2e6, td=1e-7, phase=30.0)
    if kind == 'ISin':
        return ISin(plus, minus, i=v, ia=w, freq=3e6, theta=1e5)
    if kind == 'VPulse':
        return VPulse(plus, minus, v=v, v1=0.0, v2=w, td=1e-7, tr=1e-9, tf=1e-9, pw=2e-7,
                      per=5e-7)
    if kind == 'IPulse':
        return IPulse(plus, minus, i=v, i1=w, i2=0.0, td=0.0, tr=1e-9, tf=1e-9, pw=1e-7,
                      per=3e-7)
    if kind == 'VPWL':
        return VPWL(plus, minus, v=v, tvpairs=[0.0, 0.0, 1e-7, w, 3e-7, -w])
    if kind == 'IPWL':
        return IPWL(plus, minus, i=v, tvpairs=[0.0, w, 2e-7, 0.0])
    if kind == 'VExp':
        return VExp(plus, minus, v=v, v1=0.0, v2=w, td1=1e-8, tau1=5e-8, td2=2e-7, tau2=1e-7)
    return IExp(plus, minus, i=v, i1=w, i2=0.0, td1=0.0, tau1=3e-8, td2=1e-7, tau2=5e-8)


KINDS = ('VS', 'IS', 'VSin', 'ISin', 'VPulse', 'IPulse', 'VPWL', 'IPWL', 'VExp', 'IExp')


def _pass(cir, t, an, direct):
    """`cir.u(t)`: its dtype and bytes (or the exception), its warnings with
    their lines, its `u:called` count and its paths -- by the loop where
    `direct` is False (the plan is the loop: `test_source_plan`)."""
    old = _stamp_plan.SOURCE_DIRECT, _stamp_plan.SOURCE_PLAN
    _stamp_plan.SOURCE_DIRECT, _stamp_plan.SOURCE_PLAN = direct, False
    try:
        before = _paths.snapshot()
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            try:
                v = np.asarray(cir.u(t, defaultepar, an))
                out = (v.dtype.str, v.shape, v.tobytes())
            except Exception as e:                             # noqa: BLE001
                out = (type(e).__name__, str(e))
        d = _paths.since(before)
    finally:
        _stamp_plan.SOURCE_DIRECT, _stamp_plan.SOURCE_PLAN = old
    return (out, [(w.category.__name__, str(w.message), w.filename, w.lineno) for w in W],
            d.get('u:called', 0)), d


@settings(deadline=None, max_examples=150)
@given(data=st.data())
def test_the_direct_pass_is_the_loop(data):
    """Drawn circuits: one to four sources of the family on three nodes and
    the ground (several into one node: the order of the sums), among
    resistors, capacitors, a diode and a zero-source compact model; drawn
    values (ordinary and special, of Python's and numpy's float and of int),
    times and analyses -- each pass the loop's, its warnings from the same
    lines, its `u:called` the loop's; served wherever the sources are all
    the family's."""
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    nodes = NODES + (gnd,)
    vals = st.one_of(st.sampled_from(ORDINARY), st.sampled_from(SPECIAL),
                     st.floats(allow_nan=False, allow_infinity=False, width=64))
    k = 0
    for j in range(data.draw(st.integers(1, 4))):
        for _x in range(data.draw(st.integers(0, 2))):
            kind = data.draw(st.sampled_from(('R', 'C', 'D', 'M')))
            a, b = data.draw(st.sampled_from(nodes)), data.draw(st.sampled_from(nodes))
            k += 1
            if kind == 'R':
                cir[f'x{k}'] = R(a, b, r=1e3)
            elif kind == 'C':
                cir[f'x{k}'] = C(a, b, c=1e-12)
            elif kind == 'D':
                cir[f'x{k}'] = Diode(a, b)
            else:
                cir[f'x{k}'] = eh.MosLevel1Hdl(a, b, gnd, gnd)
        kind = data.draw(st.sampled_from(KINDS))
        v = data.draw(vals)
        if data.draw(st.booleans()):
            v = np.float64(v) if isinstance(v, float) else v
        w = data.draw(st.sampled_from(ORDINARY + (0.0, -0.0, 1e300)))
        cir[f's{j}'] = _source(kind, data.draw(st.sampled_from(NODES)),
                               data.draw(st.sampled_from(nodes)), v, w)
    t = data.draw(st.sampled_from((0.0, 1e-9, 1e-7, 1.5e-7, 2.6e-7, -1e-8, 1e-3)))
    if data.draw(st.booleans()):
        t = np.float64(t)
    an = data.draw(st.sampled_from(ANALYSES))
    for _rep in range(2):
        a, _da = _pass(cir, t, an, False)
        b, db = _pass(cir, t, an, True)
        assert a == b
    ## (served, or declined as no family source -- unless a time function
    ## raised, which counts nothing on either path)
    assert (len(a[0]) == 2 or db.get('src.direct:served', 0) == 1
            or db.get('src.direct:kind', 0)), (a[0], db)


def _shared(order):
    """Sources whose values meet in one node in `order`: the sums' signed
    zeros and NaN payloads are the loop's only in its order."""
    cir = SubCircuit()
    for j, (kind, v) in enumerate(order):
        cir[f's{j}'] = IS('a', gnd, i=v) if kind == 'IS' else VS('a', 'b', v=v)
    cir['r'] = R('b', gnd, r=1.0)
    return cir


@pytest.mark.parametrize('order', [
    (('IS', 0.0), ('IS', -0.0)), (('IS', -0.0), ('IS', 0.0)), (('IS', -0.0), ('IS', -0.0)),
    (('IS', _nan(0x123)), ('IS', _nan(0x456, 1))), (('IS', _nan(0x456, 1)), ('IS', _nan(0x123))),
    (('IS', math.inf), ('IS', -math.inf)), (('IS', 1e308), ('IS', 1e308)),
    (('VS', -0.0), ('IS', 0.0), ('VS', 0.0)), (('IS', 5e-324), ('IS', -5e-324)),
    (('IS', np.float64(1.7e308)), ('IS', 1)), (('VS', np.float64(math.inf)), ('VS', 2)),
], ids=lambda o: '-'.join(f'{k}{v!r}' for k, v in o))
def test_signed_zeros_and_nan_payloads_sum_as_the_loop(order):
    circuit.default_toolkit = circuit.numeric
    cir = _shared(order)
    for an in ('tran', 'dc'):
        a, _da = _pass(cir, 1e-7, an, False)
        b, db = _pass(cir, 1e-7, an, True)
        assert a == b
        assert db.get('src.direct:served') == 1, db


def test_the_add_warns_from_u_as_the_loop():
    """A value whose add overflows or is invalid in numpy (a parameter of
    numpy's float at the range's end, an infinity against the other): the
    warning from `u`'s own line, as the loop gives it -- and the error state
    raising, raising."""
    circuit.default_toolkit = circuit.numeric
    for v, w in ((np.float64(1.7e308), 1.7e308), (np.float64(math.inf), -math.inf)):
        cir = SubCircuit()
        cir['s'] = VSin('a', gnd, v=v, va=0.0, freq=1e6, vo=0.0)
        cir['s'].function.offset = w
        cir['r'] = R('a', gnd, r=1.0)
        a, _da = _pass(cir, 1e-7, 'tran', False)
        b, db = _pass(cir, 1e-7, 'tran', True)
        assert a == b and a[1], a
        assert any(f.endswith('elements.py') for _c, _m, f, _l in a[1]), a[1]
        assert db.get('src.direct:served') == 1, db
        with np.errstate(over='raise', invalid='raise'):
            a, _da = _pass(cir, 1e-7, 'tran', False)
            b, _db = _pass(cir, 1e-7, 'tran', True)
        assert a == b and a[0][0] == 'FloatingPointError', a


def _small():
    cir = SubCircuit()
    cir['vdd'] = VS('vdd', gnd, v=1.2)
    cir['vg'] = VSin('g', gnd, v=0.7, va=2e-2, freq=1e6)
    cir['ib'] = IS('d', gnd, i=1e-5)
    cir['rl'] = R('vdd', 'd', r=5e3)
    cir['m'] = eh.MosLevel1Hdl('d', 'g', gnd, gnd)
    return cir


def _served(cir, t=2e-7, an='tran'):
    """The pass and its `src.direct` counts, against the loop's answer."""
    ref, _d = _pass(cir, t, an, False)
    got, d = _pass(cir, t, an, True)
    assert got == ref
    return {k[len('src.direct:'):]: v for k, v in d.items() if k.startswith('src.direct:')}


def test_what_the_direct_pass_declines_takes_the_loop(monkeypatch):
    """Each check that fails declines -- counted by its reason -- and the
    pass is the loop's: an instance `u` on a source or on another element,
    `u` patched on a class, a class's `function` replaced, a parameter
    shadowing its key, a source whose `u` is not the family's, an element
    whose toolkit's `array` is not the backend's, the switch -- and the pass
    is served again once the check holds; the 'ac' pass never asks it."""
    circuit.default_toolkit = circuit.numeric
    cir = _small()
    assert _served(cir) == {'served': 1}
    assert _served(cir) == {'served': 1}
    calls = []
    real = VS.u
    cir['vdd'].u = lambda *a, **k: calls.append(1) or real(cir['vdd'], *a, **k)
    assert _served(cir) == {'element': 1} and calls
    del cir['vdd'].u
    assert _served(cir) == {'served': 1}
    cir['rl'].u = lambda *a, **k: np.zeros(3)
    assert _served(cir) == {'element': 1}
    del cir['rl'].u
    assert _served(cir) == {'served': 1}
    monkeypatch.setattr(VS, 'u', lambda self, *a, **k: real(self, *a, **k))
    assert _served(cir) == {'patched': 1}
    ## (and a plan rebuilt while it is patched: the defined `u` is the one
    ## the pass stands in for, captured at its first use)
    cir['vz'] = VS('d', 'g', v=0.1)
    assert _served(cir) == {'kind': 1}
    monkeypatch.undo()
    ## (that plan holds the patched `u` -- every pass the loop's -- until it
    ## is rebuilt, as the plan's own pass: `test_source_plan`)
    assert _served(cir) == {'patched': 1}
    _stamp_plan.invalidate(cir)
    assert _served(cir) == {'served': 1}
    del cir['vz']
    assert _served(cir) == {'served': 1}
    desc = VS.__dict__['function']
    monkeypatch.setattr(VSin, 'function', property(lambda self: self.__dict__['function']),
                        raising=False)
    assert _served(cir) == {'kind': 1}
    monkeypatch.undo()
    assert VSin.__dict__.get('function') is None and VS.__dict__['function'] is desc
    assert _served(cir) == {'served': 1}
    ip = cir['ib'].iparv
    ip.__dict__['i'] = 2e-5
    assert _served(cir) == {'kind': 1}
    del ip.__dict__['i']
    assert _served(cir) == {'served': 1}
    tk2 = circuit.numeric.__class__(circuit.numeric._backend)
    tk2.array = lambda *a, **k: np.array(*a, **k)
    cir['ib'].toolkit = tk2
    assert _served(cir) == {'kind': 1}
    cir['ib'].toolkit = circuit.numeric
    assert _served(cir) == {'served': 1}
    assert _served(cir, an='ac') == {}                       # (`u` asks no direct pass)
    loop = SubCircuit._add_element_subvectors
    cir._add_element_subvectors = lambda *a, **k: loop(cir, *a, **k)
    assert _served(cir) == {'patched': 1}
    del cir._add_element_subvectors
    assert _served(cir) == {'served': 1}
    monkeypatch.setattr(SubCircuit, '_scatter_1d', staticmethod(SubCircuit._scatter_1d))
    assert _served(cir) == {'served': 1}                     # (the same function)
    monkeypatch.setattr(SubCircuit, '_scatter_1d',
                        staticmethod(lambda *a: type(cir).__mro__[1]._scatter_1d(*a)))
    assert _served(cir) == {'patched': 1}
    monkeypatch.undo()
    assert _served(cir) == {'served': 1}
    monkeypatch.setattr(_stamp_plan, 'SOURCE_DIRECT', False)
    before = _paths.snapshot()
    cir.u(2e-7, defaultepar, 'tran')
    assert _paths.since(before).get('src.direct:off') == 1
    monkeypatch.undo()
    other = _small()
    other['x'] = _UOwn('a', gnd)
    assert _served(other) == {'kind': 1}
    assert _served(other) == {'kind': 1}


class _UOwn(VS):
    """A source with a `u` of its own: not the family's expression."""

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        return self.toolkit.array([0.0, 0.0, -2.0 * t])


def test_what_changes_is_read_on_the_next_pass():
    """A parameter written, a time function replaced, an element added, a
    class changed under an element: the next pass reads it -- the loop's
    answer each time."""
    circuit.default_toolkit = circuit.numeric
    cir = _small()
    assert _served(cir) == {'served': 1}
    cir['vdd'].iparv.v = 0.9
    assert _served(cir) == {'served': 1}
    cir['vg'].function = func.Sin(0.1, 0.3, 4e6, 0.0, 0.0, 10.0)
    assert _served(cir) == {'served': 1}
    cir['vx'] = VPulse('g', 'd', v1=0.0, v2=1.0, td=1e-8, tr=1e-9, tf=1e-9, pw=1e-7, per=2e-7)
    assert _served(cir) == {'served': 1}
    cir['rl'].__class__ = _RU
    assert _served(cir) == {'element': 1}
    cir['rl'].__class__ = R
    assert _served(cir) == {'served': 1}


class _RU(R):
    def u(self, t=0.0, epar=defaultepar, analysis=None):
        return self.toolkit.array([1e-6, -1e-6])


def test_the_stamp_holds_where_a_limiter_writes_its_dict():
    """A diode writes its limiting state into its dict every iteration; the
    pass's stamp watches the sources' dicts only, and holds: a transient of
    an RC ladder ending in a diode arms it twice."""
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['vs'] = VSin('n0', gnd, va=2.0, freq=1e3)
    for k in range(6):
        cir[f'R{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
        cir[f'C{k}'] = C(f'n{k + 1}', gnd, c=1e-8)
    cir['D'] = Diode('n6', gnd)
    tr = Transient(cir, integrator=Gear2Integrator(), reltol=1e-5)
    before = _paths.snapshot()
    tr.solve(tend=3e-4, timestep=1e-5)
    d = _paths.since(before)
    sd = cir.__dict__['_stamp_plan'].methods['src:u'].direct
    assert d.get('src.direct:served', 0) > 30 and sd.arms <= 2, (d, sd.arms)


@pytest.mark.parametrize('method', ['gear', 'radau'])
def test_a_transient_and_a_pss_are_the_same(method, monkeypatch):
    """The PSP stage (a DC and a sine source): a transient and a PSS by gear
    and by radau, the direct pass and the scalar sine off and on -- the
    waveform, statistics and warnings the same."""
    from pycircuit.circuit.shooting import PSS
    from pycircuit.circuit.tests.test_psp_limit_c import _stage

    def run(on):
        monkeypatch.setattr(_stamp_plan, 'SOURCE_DIRECT', on)
        monkeypatch.setattr(func, 'SIN_SCALAR', on)
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            integ = Gear2Integrator() if method == 'gear' else RadauIIA3Integrator()
            tr = Transient(_stage(), toolkit=circuit.numeric, integrator=integ)
            res = tr.solve(tend=4e-7, timestep=2e-8)
            p = PSS(_stage(), method=method, reltol=1e-8)
            before = _paths.snapshot()
            p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
            d = _paths.since(before)
        st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__
               if 'seconds' not in k}
        return (np.asarray(res.x, float).tobytes(), st_,
                np.asarray(p.waveform[1], float).tobytes(),
                sorted(str(w.message) for w in W)), d
    a, _da = run(False)
    b, db = run(True)
    assert a == b
    assert db.get('src.direct:served', 0) > 100, db


## -- `Sin.f` at a scalar time -------------------------------------------------------------

def _sin_call(S, t, scalar):
    old = func.SIN_SCALAR
    func.SIN_SCALAR = scalar
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            try:
                v = S.f(t)
                out = (type(v).__name__, np.asarray(v).dtype.str, np.asarray(v).tobytes())
            except Exception as e:                             # noqa: BLE001
                out = (type(e).__name__, str(e))
    finally:
        func.SIN_SCALAR = old
    return out, [(w.category.__name__, str(w.message), w.filename, w.lineno) for w in W]


PARAMS = (0, 0.0, -0.0, 1, 0.7, -2.5, 1e-3, 1e5, -1e5, 1e300, -1e300, 1.7e308, 5e-324,
          math.inf, -math.inf, _nan(0), _nan(0x123, 1))


@settings(deadline=None, max_examples=400)
@given(data=st.data())
def test_the_scalar_sine_is_the_expression(data):
    """Drawn parameters (ordinary and special, of int, Python's and numpy's
    float) and times: `Sin.f` the same type and bits, and the same warnings
    from the same lines, with the scalar form on and off -- the form serving
    where nothing could warn, the expression where something could."""
    p = [data.draw(st.one_of(st.sampled_from(PARAMS), st.floats(-1e3, 1e3)))
         for _ in range(6)]
    p = [np.float64(v) if isinstance(v, float) and data.draw(st.booleans()) else v for v in p]
    with np.errstate(all='ignore'):
        S = func.Sin(p[0], p[1], p[2], p[3], p[4], p[5])
    t = data.draw(st.one_of(st.sampled_from((0.0, 1e-7, -1e-7, 1e-3, math.inf, _nan(0))),
                            st.floats(-1e-5, 1e-5)))
    if data.draw(st.booleans()):
        t = np.float64(t)
    assert _sin_call(S, t, True) == _sin_call(S, t, False)


def test_the_scalar_sine_serves_and_declines():
    """Served for a scalar time on the numeric toolkit; the expression for
    an array of times, numpy's float32, a parameter of another type, a time
    or exponent past its bounds, numpy's underflow not ignored -- whatever
    state object holds the setting."""
    S = func.Sin(0.1, 0.5, 1e6, 0.0, 0.0, 0.0)
    assert func._sin_scalar(S, 1.25e-7) is not None
    assert func._sin_scalar(S, np.float64(1.25e-7)) is not None
    assert func._sin_scalar(S, np.array([1.25e-7])) is None
    assert func._sin_scalar(S, np.float32(1.25e-7)) is None
    assert func._sin_scalar(S, 1e301) is None
    S2 = func.Sin(0.1, 0.5, 1e6, 0.0, -1e10, 0.0)
    assert func._sin_scalar(S2, 1e-6) is None              # (exp(1e4): the expression warns)
    S3 = func.Sin(np.float32(0.1), 0.5, 1e6, 0.0, 0.0, 0.0)
    assert func._sin_scalar(S3, 1e-7) is None
    with np.errstate(under='warn'):
        assert func._sin_scalar(S, 1.25e-7) is None
        with np.errstate(over='raise'):
            assert func._sin_scalar(S, 1.25e-7) is None
    ## (overflow and invalid cannot happen where it serves: their state is
    ## not read; a state set for good -- `seterr` -- is read by its value)
    with np.errstate(over='raise', invalid='raise'):
        assert func._sin_scalar(S, 1.25e-7) is not None
    old = np.seterr(divide='ignore')
    try:
        assert func._sin_scalar(S, 1.25e-7) is not None
        np.seterr(under='warn')
        assert func._sin_scalar(S, 1.25e-7) is None
    finally:
        np.seterr(**old)
    assert func._sin_scalar(S, 1.25e-7) is not None
