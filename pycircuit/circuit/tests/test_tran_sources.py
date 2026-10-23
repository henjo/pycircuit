"""The sources in C (the C stepping loop plan's stage 1, 2026-10-10):
`_tran_sources.Table` against the direct source pass it stands for
(`_stamp_plan._SourceDirect.run`) and the waveforms' own `next_event`.

Held: wherever the C answers, `u(t)` is the Python's byte for byte and the
next event the earliest of the Python's -- on every waveform (V and I) over
a dense grid and its corners, and on drawn parameters and times, one
waveform at a time (hypothesis; the `gate` profile, derandomized).  Each
decline where the Python leaves the fast float route; a table built from
inputs that changed since does not hold; a planted difference is caught.
"""
import math

import numpy as np
import pytest

hypothesis = pytest.importorskip('hypothesis')
from hypothesis import given
from hypothesis import strategies as st

from pycircuit.circuit import _paths, func, gnd
from pycircuit.circuit import _tran_sources as ts
from pycircuit.circuit import elements as E
from pycircuit.circuit.circuit import SubCircuit, defaultepar
from pycircuit.circuit.tests._warnpolicy import quiet

pytestmark = pytest.mark.skipif(ts.driver() is None, reason='the sources unit did not build')


def _direct(c):
    """The circuit's direct source pass (built by a first call)."""
    c.u(0.0, defaultepar, analysis='tran')
    return c.__dict__['_stamp_plan'].methods['src:u'].direct


def _with_loads(c, nodes):
    for k in nodes:
        c['R' + k] = E.R(k, gnd, r=1e3)
    return c


def _all_kinds():
    c = SubCircuit()
    c['v0'] = E.VS('a', gnd, v=1.5)
    c['v1'] = E.VSin('b', gnd, vo=0.1, va=1.0, freq=1e6, td=1e-7, theta=1e4, phase=30.0)
    c['v2'] = E.VPulse('c', gnd, v1=0.0, v2=1.8, td=1e-9, tr=1e-10, tf=2e-10, pw=5e-9, per=2e-8)
    c['v3'] = E.VPWL('d', gnd, tvpairs=[0.0, 0.0, 1e-9, 1.0, 3e-9, 1.0, 4e-9, -0.5])
    c['v4'] = E.VExp('e', gnd, v1=0.0, v2=1.0, td1=1e-9, tau1=2e-9, td2=6e-9, tau2=3e-9)
    c['v5'] = E.VAM('f', gnd, vo=0.1, va=0.5, fc=1e8, fm=1e7, m=0.3)
    c['v6'] = E.VSFFM('g', gnd, vo=0.0, va=1.0, fc=1e8, mdi=2.0, fm=1e7)
    c['i0'] = E.IPulse('h', gnd, i1=0.0, i2=1e-3, td=2e-9, tr=1e-10, tf=1e-10, pw=3e-9, per=1e-8)
    c['i1'] = E.ISin('k', gnd, io=1e-4, ia=1e-3, freq=5e6)
    c['i2'] = E.IPWL('m', gnd, tvpairs=[0.0, 1e-3, 2e-9, -1e-3])
    return _with_loads(c, 'abcdefghkm')


def _functions(c):
    return [el.function for _i, el in c.elements.items() if hasattr(type(el), 'function')]


def _same(c, tab, t, functions=None):
    """`tab` against the Python at `t`: 'same', 'declined' or a difference
    (raises)."""
    py = c.u(t, defaultepar, analysis='tran')
    cc = tab.u(t)
    if cc is None:
        return 'declined'
    assert cc.tobytes() == py.tobytes(), (t, py, cc)
    fns = functions if functions is not None else _functions(c)
    try:
        nxt = min(f.next_event(t) for f in fns)
    except ValueError:
        nxt = None                      # (the Python raises: a NaN time's period index)
    got = tab.next_event(t)
    if nxt is None:
        assert got is None, (t, got)
    elif got is not None:
        assert got == nxt or (math.isinf(got) and math.isinf(nxt)), (t, nxt, got)
    return 'same'


def _corners():
    """The all-kinds circuit's edges, a ulp either side, and their periods."""
    edges = [0.0, 1e-9, 1.1e-9, 6.1e-9, 6.3e-9, 2e-9, 2.1e-9, 3e-9, 4e-9, 6e-9, 1e-7, 2e-8, 4e-8]
    out = []
    for e in edges:
        out += [e, math.nextafter(e, -math.inf), math.nextafter(e, math.inf)]
    return out


def test_every_waveform_is_the_python_s_on_a_dense_grid_and_its_corners():
    c = _all_kinds()
    tab = ts.table_for(_direct(c))
    assert tab is not None
    fns = _functions(c)
    times = [float(t) for t in np.linspace(0.0, 1.2e-7, 3001)] + _corners()
    got = [_same(c, tab, t, fns) for t in times]
    assert got.count('declined') == 0, got.count('declined')


## -- drawn: one waveform at a time ------------------------------------------

_fin = st.floats(-1e3, 1e3, allow_nan=False, allow_infinity=False)
_pos = st.floats(1e-12, 1e-6, allow_nan=False)
_odd = st.sampled_from([0.0, -0.0, 5e-324, 1e-300, 1e299, -1e299, 1e301, math.inf, -math.inf,
                        math.nan])
_val = st.one_of(_fin, _odd, st.integers(-5, 5))
_time = st.one_of(st.floats(-1e-6, 1e-5, allow_nan=False), _odd)


def _one(make):
    c = SubCircuit()
    c['s'] = make()
    return _with_loads(c, 'a')


def _drawn(c, t):
    with quiet():
        tab = ts.table_for(_direct(c))
        if tab is None:
            return 'declined'
        return _same(c, tab, float(t))


@given(v1=_val, v2=_val, td=_pos, tr=_pos, tf=_pos, pw=_pos, per=st.one_of(st.just(0.0), _pos),
       t=_time, current=st.booleans())
def test_drawn_pulse(v1, v2, td, tr, tf, pw, per, t, current):
    S = (lambda: E.IPulse('a', gnd, i1=v1, i2=v2, td=td, tr=tr, tf=tf, pw=pw, per=per)) if current \
        else (lambda: E.VPulse('a', gnd, v1=v1, v2=v2, td=td, tr=tr, tf=tf, pw=pw, per=per))
    _drawn(_one(S), t)


@given(vo=_val, va=_val, freq=st.floats(0.0, 1e9), td=_pos, theta=st.floats(0.0, 1e7),
       phase=st.floats(-360.0, 360.0), t=_time)
def test_drawn_sin(vo, va, freq, td, theta, phase, t):
    _drawn(_one(lambda: E.VSin('a', gnd, vo=vo, va=va, freq=freq, td=td, theta=theta,
                               phase=phase)), t)


@given(pts=st.lists(st.tuples(st.floats(0.0, 1e-6), _val), min_size=1, max_size=6), t=_time,
       v=st.one_of(_fin, st.integers(-3, 3)))
def test_drawn_pwl(pts, t, v):
    pts = sorted(pts, key=lambda p: p[0])
    flat = [x for p in pts for x in p]
    _drawn(_one(lambda: E.VPWL('a', gnd, v=v, tvpairs=flat)), t)


@given(v1=_val, v2=_val, td1=_pos, tau1=_pos, dtd=_pos, tau2=_pos, t=_time)
def test_drawn_exp(v1, v2, td1, tau1, dtd, tau2, t):
    _drawn(_one(lambda: E.VExp('a', gnd, v1=v1, v2=v2, td1=td1, tau1=tau1, td2=td1 + dtd,
                               tau2=tau2)), t)


@given(vo=_val, va=_val, fc=st.floats(0.0, 1e9), fm=st.floats(0.0, 1e8), m=_fin, t=_time,
       sffm=st.booleans())
def test_drawn_am_and_sffm(vo, va, fc, fm, m, t, sffm):
    S = (lambda: E.VSFFM('a', gnd, vo=vo, va=va, fc=fc, mdi=m, fm=fm)) if sffm \
        else (lambda: E.VAM('a', gnd, vo=vo, va=va, fc=fc, fm=fm, m=m))
    _drawn(_one(S), t)


## -- the declines ----------------------------------------------------------

def _declined(c):
    before = dict(_paths.COUNTS)
    assert ts.table_for(_direct(c)) is None
    return {k for k, v in _paths.COUNTS.items() if k.startswith('sources_c:') and v != before.get(k, 0)}


def test_each_route_the_python_takes_otherwise_declines():
    """An int DC value (the Python adds two ints), an exponential's zero
    `tau` and a pulse's zero edge (the Python raises), a waveform subclass
    (its `f` may be another), a libm that differs from numpy's on this CPU."""
    assert _declined(_one(lambda: E.VS('a', gnd, v=2))) == {'sources_c:int'}
    assert _declined(_one(lambda: E.VExp('a', gnd, v1=0.0, v2=1.0, td1=1e-9, tau1=0.0,
                                         td2=2e-9, tau2=1e-9))) == {'sources_c:tau'}
    for edge in ('tr', 'tf'):
        ## (zeroed after the direct pass is built: building it calls the
        ## Python, which raises)
        c = _one(lambda: E.VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-9, tr=1e-10, tf=1e-10,
                                  pw=1e-9, per=0.0))
        sd = _direct(c)
        setattr(c.elements['s'].function, edge, 0.0)
        before = dict(_paths.COUNTS)
        assert ts.table_for(sd) is None
        assert {k for k, v in _paths.COUNTS.items() if k.startswith('sources_c:')
                and v != before.get(k, 0)} == {'sources_c:edge'}, edge

    class MyPulse(func.Pulse):
        pass
    c = _one(lambda: E.VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-9, tr=1e-10, tf=1e-10, pw=1e-9,
                              per=0.0))
    el = c.elements['s']
    el.function = MyPulse(0.0, 1.0, 1e-9, 1e-10, 1e-10, 1e-9, 0.0)
    assert _declined(c) == {'sources_c:class'}
    drv = ts.driver()
    try:
        ts._driver = (drv[0], drv[1], frozenset({'sin'}))
        assert _declined(_one(lambda: E.VSin('a', gnd, vo=0.0, va=1.0, freq=1e6))) \
            == {'sources_c:libm'}
    finally:
        ts._driver = drv


def test_a_time_function_of_its_own_is_the_pythons():
    """A time function with its own `f` or `next_event` on the instance -- a
    user's injected noise -- is declined: the C makes its class's waveform
    (the C loop's table once did, the run another).  A run through the C
    loop and one through the direct pass are the Python's."""
    for name in ('f', 'next_event'):
        c = _one(lambda: E.VSin('a', gnd, vo=0.0, va=1.0, freq=1e6))
        fn = c.elements['s'].function
        setattr(fn, name, getattr(fn, name))
        assert _declined(c) == {'sources_c:own'}, name
    from pycircuit.circuit import _tran_loop_c, circuit
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    from pycircuit.circuit.transient import Transient

    def run(loop, srcs, integ=None):
        c = SubCircuit()
        c['In'] = E.IS(gnd, 'o', i=0.0)
        c['In'].function.f = lambda t: 1e-3 * math.sin(2e6 * math.pi * t) + 2e-4
        c['R'] = E.R('o', gnd, r=1e3)
        c['C'] = E.C('o', gnd, c=1e-10)
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(_tran_loop_c, 'ENABLED', loop)
            mp.setattr(ts, 'ENABLED', srcs)
            tr = Transient(c, toolkit=circuit.numeric,
                           **({} if integ is None else {'integrator': integ()}))
            return np.asarray(tr.solve(tend=3e-6, timestep=1e-8).x).tobytes()
    assert run(True, True) == run(False, False)
    assert run(True, True, RadauIIA3Integrator) == run(True, False, RadauIIA3Integrator)


def test_a_call_off_the_fast_route_declines():
    """numpy's error state raising on underflow (the Python's exp would);
    a NaN time where a PWL finds no segment (the Python's None raises)."""
    c = _one(lambda: E.VExp('a', gnd, v1=0.0, v2=1.0, td1=1e-9, tau1=1e-9, td2=5e-9,
                            tau2=1e-9))
    tab = ts.table_for(_direct(c))
    with np.errstate(under='raise'):
        assert tab.u(3e-9) is None
    assert tab.u(3e-9) is not None
    c = _one(lambda: E.VPWL('a', gnd, tvpairs=[0.0, 0.0, 1e-9, 1.0]))
    tab = ts.table_for(_direct(c))
    assert tab.u(math.nan) is None
    assert _paths.COUNTS['sources_c:bail:pwl'] >= 1


def test_a_table_does_not_hold_once_an_input_changed():
    """A waveform's attribute set, a PWL's points edited in place, the DC
    value set, the function replaced: each read again by the Python, so
    the table built before does not hold."""
    c = _all_kinds()
    els = c.elements
    for change in (lambda: setattr(els['v2'].function, 'td', 2e-9),
                   lambda: els['v3'].function.values.__setitem__(1, 0.5),
                   lambda: els['v0'].iparv.__setattr__('v', 1.25),
                   lambda: setattr(els['v4'], 'function', func.TimeFunction())):
        tab = ts.table_for(_direct(c))
        assert tab is not None and tab.holds()
        change()
        assert not tab.holds()


def test_the_comparison_sees_a_planted_difference():
    """Alive: one parameter a ulp off in the table makes another vector."""
    c = _all_kinds()
    tab = ts.table_for(_direct(c))
    par = tab.arrays['par']
    k = int(tab.arrays['poff'][2])                      # the pulse's v2
    par[k + 1] = math.nextafter(par[k + 1], math.inf)
    with pytest.raises(AssertionError):
        _same(c, tab, 3e-9)
