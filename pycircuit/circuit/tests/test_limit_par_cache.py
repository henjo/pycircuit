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
    """With the cache off every use evaluates: the ranking once per
    function, and the apply loop again for a probe whose limit the
    ranking's could not stand in for (another probe wrote its terminal
    first -- stage L3 reuses the ranking's limit when the inputs are the
    same bits)."""
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
        assert n <= len(calls) <= 2 * n, (len(calls), n)
        ## a second call evaluates again: nothing is kept
        del calls[:]
        e.limit(x0 + 0.6, x0, defaultepar)
        assert len(calls) >= n
    finally:
        spec[:] = saved


## ----------------------------------------------------------------------
## Stage L2: the parameter chains as twins (`_hdl_cse.optimise_limit_pars`).

@pytest.mark.parametrize('name', ['MosLevel1Hdl', 'MosLevel3Hdl',
                                  'GummelPoonNpnHdl', 'EkvNmosHdl'])
def test_the_limiter_parameter_chains_run_as_bit_identical_twins(name):
    """Every chain-compiled limiter parameter is wrapped over its CSE +
    fast-path twin: the wrapper keeps `_wants_x` and the JAX twin's
    ingredients, names the twin as its inner (with the reference on it),
    and answers the reference's bytes over a sweep."""
    from pycircuit.circuit import _hdl_cse as cs
    cls = getattr(eh, name)
    e = (cls('c', 'b', gnd) if name == 'GummelPoonNpnHdl'
         else cls('d', 'g', gnd, gnd))
    e.update_iparv()
    args = list(hdl._args_of(e, defaultepar))
    rng = np.random.default_rng(0)
    pts = [rng.uniform(-2, 2, e.n) for _ in range(30)]
    pts += [np.full(e.n, v) for v in (0.0, -0.0, 1e30, np.inf, np.nan)]
    twins = 0
    for rows, kind, move, pfs in type(e)._hdl_info['limit_spec']:
        for f in pfs:
            inner = getattr(f, '_hdl_inner', None)
            if inner is None:
                continue                      # a lambdified parameter
            assert '_hdl_ref' in inner.__dict__, (name, kind)
            assert inner._src == inner._hdl_ref._src
            assert '_hdl_limit_par' in f.__dict__
            ref = inner._hdl_ref
            wx = getattr(f, '_wants_x', False)
            for x in pts:
                with np.errstate(all='ignore'):
                    a = np.asarray(ref(x, *args) if wx else ref(*args), float)
                    b = np.asarray(f(x, *args) if wx else f(*args), float)
                assert a.tobytes() == b.tobytes(), (name, kind, x)
            twins += 1
    assert twins >= 1
    assert cs.ENABLED


## ----------------------------------------------------------------------
## Stage L3: the streamlined limiter body (2026-10-02) against a reference
## transliteration of the body as it was -- numpy-scalar branch voltages,
## every probe's limit computed twice -- over a sweep.

def _reference_limit(e, x, x0, epar=defaultepar):
    """The generated `limit()` as it was before stage L3, written out: the
    same `apply_limit` / `device_writeback`, the same canonical orders, the
    probes' parameter values from the element's own cache, and each
    probe's limit computed once for the ranking and again to apply."""
    from pycircuit.circuit._limiting import apply_limit as _lim
    from pycircuit.circuit._limiting import device_writeback as _dwb
    info = type(e)._hdl_info
    _ls = info['limit_spec']
    _lg = info.get('limit_groups') or []
    grouped = set()
    for _s, ix in _lg:
        grouped.update(ix)
    _l1 = [i for i in range(len(_ls)) if i not in grouped]
    out = np.array(x, dtype=float, copy=True)
    x0a = np.asarray(x0, dtype=float)
    args = hdl._args_of(e, epar)
    pv_of = {j: [float(f(x0a, *args)) if getattr(f, '_wants_x', False)
                 else float(f(*args)) for f in _ls[j][3]]
             for j in range(len(_ls))}
    drift = np.abs(np.asarray(out, dtype=float) - x0a)
    moved = set()

    def vlim_of(j):
        (i0, i1), k = _ls[j][0], _ls[j][1]
        vn = float(out[i0] - out[i1])
        vo = float(x0a[i0] - x0a[i1])
        return vn, _lim(k, vn, vo, pv_of[j], e.toolkit)
    for seq, idx in _lg:
        targets, shift, taken = [], {}, set()
        if seq:
            order = list(idx)
        else:
            rk = {j: vlim_of(j) for j in idx}
            order = sorted(idx, key=lambda j: (-abs(rk[j][1] - rk[j][0]),
                                               _ls[j][0][0], _ls[j][0][1]))
        for j in order:
            (ra, rb), kind, _mv, _pfs = _ls[j]
            vorig = float(out[ra] - out[rb])
            vold = float(x0a[ra] - x0a[rb])
            vin = vorig + shift.get(ra, 0.0) - shift.get(rb, 0.0)
            vlim = _lim(kind, vin, vold, pv_of[j], e.toolkit)
            if vlim != vin:
                if seq:
                    n = rb
                else:
                    n = ra if drift[ra] >= drift[rb] else rb
                    if n in taken:
                        n = rb if n == ra else ra
                taken.add(n)
                shift[n] = shift.get(n, 0.0) + (vlim - vin) * (1 if n == ra
                                                               else -1)
            targets.append((ra, rb, vorig, vlim))
        moved |= _dwb(out, targets, drift, moved)
    rank = {i: vlim_of(i) for i in _l1}
    order = sorted(_l1, key=lambda i: (-abs(rank[i][1] - rank[i][0]),
                                       _ls[i][0][0], _ls[i][0][1]))
    for i in order:
        (ra, rb), kind, move, _pfs = _ls[i]
        vnew = float(out[ra] - out[rb])
        vold = float(x0a[ra] - x0a[rb])
        vlim = _lim(kind, vnew, vold, pv_of[i], e.toolkit)
        if vlim == vnew:
            continue
        cand = ra if drift[ra] >= drift[rb] else rb
        if cand in moved:
            cand = rb if cand == ra else ra
        if cand in moved:
            cand = move
        moved.add(cand)
        if cand == ra:
            out[ra] = out[rb] + vlim
        else:
            out[rb] = out[ra] - vlim
    return out


@pytest.mark.parametrize('name', ['MosLevel1Hdl', 'MosLevel3Hdl',
                                  'GummelPoonNpnHdl', 'EkvNmosHdl',
                                  'DiodeSpiceHdl'])
def test_the_streamlined_body_answers_the_old_one_bit_for_bit(name):
    """Random states and steps of every size, ties (`x == x0`), signed
    zeros, rails: the same limited vector to the last bit, probe groups
    (MosLevel1/3), singles (the others) and a `limit_together` with a
    sequential order (through `test_device_limiter`'s own model) alike."""
    cls = getattr(eh, name)
    e = (cls('c', 'b', gnd) if name == 'GummelPoonNpnHdl'
         else cls('a', 'k') if name == 'DiodeSpiceHdl'
         else cls('d', 'g', gnd, gnd))
    e.update_iparv()
    rng = np.random.default_rng(7)
    n = e.n
    for _ in range(400):
        scale = rng.choice([1e-3, 0.1, 1.0, 10.0, 100.0])
        x0 = rng.uniform(-2.0, 2.0, n)
        x = x0 + scale * rng.standard_normal(n)
        k = rng.integers(0, 4)
        if k == 1:
            x = x0.copy()                               # a tie everywhere
        elif k == 2:
            x[rng.integers(0, n)] = -0.0                # a signed zero
            x0[rng.integers(0, n)] = 0.0
        elif k == 3:
            x[:2] = 50.0 * np.sign(x[:2])               # the rails
        a = _reference_limit(e, x, x0)
        b = e.limit(x, x0, defaultepar)
        assert a.tobytes() == b.tobytes(), (name, x, x0)
