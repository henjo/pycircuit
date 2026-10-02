"""The generated `limit()` evaluates each probe's parameter functions ONCE
per call, and keeps the ones that read no solution until the element's
parameter list changes (speed round 3, stage L1; 2026-10-02).

Until then every value was computed twice per probe per call -- the ranking
of the probes and the apply loop -- through lambdas over the model's chains:
on a 20-MosLevel1 chain the limiter was 62 % of a gear step.  The cache is
the same functions on the same objects, so the limited point, and every
result after it, is the same to the last bit (`LIMIT_PAR_CACHE = False` is
the old behaviour, for the comparison).
"""
import numpy as np
import pytest

from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VSin, gnd
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.transient import Transient


def _stage(make):
    cm.default_toolkit = numeric
    c = SubCircuit()
    for n in ('g', 'd', 'vdd'):
        c.add_node(n)
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['vg'] = VSin('g', gnd, v=0.9, va=0.4, freq=1e6)
    c['rl'] = R('vdd', 'd', r=5e3)
    c['cl'] = C('d', gnd, c=1e-13)
    c['M'] = make()
    return c


@pytest.mark.parametrize('name', ['MosLevel1Hdl', 'GummelPoonNpnHdl',
                                  'EkvNmosHdl'])
def test_a_transient_is_byte_identical_with_the_cache_off(name, monkeypatch):
    """MosLevel1 (one `limit_together` group, one parameter reading `x0`),
    Gummel-Poon (singles, none reading `x0`), EKV (singles, one reading
    `x0`): the same bytes with the values kept and with every use
    evaluating them anew."""
    cls = getattr(eh, name)
    make = (lambda: cls('d', 'g', gnd)) if name == 'GummelPoonNpnHdl' else \
        (lambda: cls('d', 'g', gnd, gnd))

    def run():
        res = Transient(_stage(make), toolkit=numeric).solve(
            tend=1e-6, timestep=1e-8, fixed_timestep=True)
        return np.asarray(res.x, float)
    a = run()
    monkeypatch.setattr(hdl, 'LIMIT_PAR_CACHE', False)
    b = run()
    assert a.tobytes() == b.tobytes()


def _counting(cls, calls):
    """Every parameter function of `cls`'s limiter wrapped to count, in the
    list the generated `limit()` reads (`info['limit_spec']`, its entries
    replaced -- the closure holds the list)."""
    spec = cls._hdl_info['limit_spec']
    saved = list(spec)

    def wrap(j, k, f):
        def g(*a):
            calls.append((j, k, a[0] if getattr(f, '_wants_x', False) else None))
            return f(*a)
        g.__dict__.update(getattr(f, '__dict__', {}))
        return g
    for j, (rows, kind, move, pfs) in enumerate(saved):
        spec[j] = (rows, kind, move,
                   tuple(wrap(j, k, f) for k, f in enumerate(pfs)))
    return saved


def test_the_values_are_kept_until_the_parameters_or_the_temperature_move():
    cls = eh.MosLevel1Hdl
    e = cls('d', 'g', gnd, gnd)
    e.update_iparv()
    calls = []
    saved = _counting(type(e), calls)
    spec = type(e)._hdl_info['limit_spec']
    try:
        x0 = np.full(e.n, 0.3)
        x1 = x0 + 0.5
        wants = [(j, k) for j, s in enumerate(spec)
                 for k, f in enumerate(s[3]) if getattr(f, '_wants_x', False)]
        pure = [(j, k) for j, s in enumerate(spec)
                for k, f in enumerate(s[3]) if not getattr(f, '_wants_x', False)]
        assert wants and pure
        e.limit(x1, x0, defaultepar)
        first = [c[:2] for c in calls]
        ## every function once, the pure ones kept from here on
        assert sorted(first) == sorted(wants + pure)
        for _ in range(3):
            del calls[:]
            e.limit(x1 + 0.1, x0, defaultepar)
            assert sorted(c[:2] for c in calls) == sorted(wants)
            ## the one reading the solution sees x0, the last accepted point
            assert all(c[2] is not None and np.array_equal(c[2], x0)
                       for c in calls)
        ## a new x0: still only the x0 readers
        del calls[:]
        e.limit(x1, x0 + 0.2, defaultepar)
        assert sorted(c[:2] for c in calls) == sorted(wants)
        ## a parameter change: everything again, once
        del calls[:]
        e.ipar.vto = float(e.iparv.vto) + 0.05
        e.update_iparv()
        e.limit(x1, x0, defaultepar)
        assert sorted(c[:2] for c in calls) == sorted(wants + pure)
        ## a temperature change: everything again, once
        del calls[:]
        from types import SimpleNamespace
        e.limit(x1, x0, SimpleNamespace(T=320.0))
        assert sorted(c[:2] for c in calls) == sorted(wants + pure)
    finally:
        spec[:] = saved


def test_a_single_probe_model_computes_each_value_once_per_call():
    cls = eh.GummelPoonNpnHdl
    e = cls('c', 'b', gnd)
    e.update_iparv()
    calls = []
    saved = _counting(type(e), calls)
    spec = type(e)._hdl_info['limit_spec']
    try:
        n = sum(len(s[3]) for s in spec)
        x0 = np.full(e.n, 0.2)
        e.limit(x0 + 0.6, x0, defaultepar)
        assert len(calls) == n           # once each, where it was twice
        del calls[:]
        e.limit(x0 + 0.7, x0, defaultepar)
        assert calls == []               # none reads x0: all kept
    finally:
        spec[:] = saved


def test_the_cache_off_evaluates_at_every_use(monkeypatch):
    monkeypatch.setattr(hdl, 'LIMIT_PAR_CACHE', False)
    cls = eh.GummelPoonNpnHdl
    e = cls('c', 'b', gnd)
    e.update_iparv()
    calls = []
    saved = _counting(type(e), calls)
    spec = type(e)._hdl_info['limit_spec']
    try:
        n = sum(len(s[3]) for s in spec)
        x0 = np.full(e.n, 0.2)
        e.limit(x0 + 0.6, x0, defaultepar)
        assert len(calls) == 2 * n       # the ranking and the apply loop
    finally:
        spec[:] = saved
