"""The hdl limiter as a C kernel (`_hdl_climit`; speed round 5, 2026-10-03).

The kernel is proven against the Python closure it replaces, class by
class and law by law: every library class that declares `$limit` (and
its collapse variant the library instantiates), the five laws through a
probe kernel, a chained `limit_together(sequential=True)` model, the
inputs the closure accepts (lists, int vectors, a longer vector, `limit(x,
x)`), the temperatures it serves, and the calls it DECLINES -- a NaN or
an infinity among the keys it would sort -- which the closure answers.
`test_a_wrong_law_is_caught` shows the sweep can fail.
"""
import contextlib
import pickle
import types
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _hdl_cache, _limiting, hdl
from pycircuit.circuit import _hdl_cbackend as cb
from pycircuit.circuit import _hdl_climit as cl
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.hdl import Node
from pycircuit.circuit.tests.test_hdl_cbackend import _compare, numpy_backend
from pycircuit.circuit.tests.test_limit_fet import _fet
from pycircuit.circuit.toolkit import numeric

cm.default_toolkit = numeric


def _inst(cls, **kw):
    e = cls(*[Node(f'n{k}') for k in range(len(cls.terminals))], **kw)
    e.update_iparv()
    return e


def _limiting_classes():
    out = []
    for name in sorted(dir(eh)):
        cls = getattr(eh, name)
        if not (isinstance(cls, type) and issubclass(cls, hdl.Behavioural)):
            continue
        info = getattr(cls, '_hdl_info', None)
        if info and info.get('limit_spec') and info.get('chained'):
            out.append(name)
    return out


def _bound(e):
    """The element's class info when C is bound, else skip (no compiler
    and no store on this machine)."""
    info = type(e)._hdl_info
    if not info.get('_c_bound'):
        pytest.skip(f"no C backend here: {type(e)._hdl_backend_status}")
    return info


def _py_limit(e, x, x0, epar=defaultepar):
    """The closure's answer, its numpy warnings quiet (a parameter chain
    evaluates both arms of a `where`)."""
    with warnings.catch_warnings(), np.errstate(all='ignore'):
        warnings.simplefilter('ignore', RuntimeWarning)
        return e.limit(x, x0, epar)


## -- the laws ----------------------------------------------------------------

_PROBE = cl._LIMIT_C + '''
int hdl_fn(const double *x, const double *p, double *out) {
  out[0] = _lim((int)p[0], x[0], x[1], x[2], x[3]);
  return 0;
}
'''


@pytest.fixture(scope='module')
def law():
    """`law(kind, vnew, vold, p0, p1)` through a probe kernel of the
    limiter's prelude."""
    try:
        ffi, cfn, _key, _cold, _s = cb.load_kernel(_PROBE, cl.limit_cdef())
    except cb.CompileError as e:
        pytest.skip(f'no C compiler: {e}')
    dptr = ffi.typeof('double *')

    def call(kind, vnew, vold, p0=0.0, p1=0.0):
        x = np.array([vnew, vold, p0, p1], dtype=float)
        p = np.array([float(cl.KINDS[kind])])
        out = np.empty(1)
        cfn(ffi.from_buffer(dptr, x), ffi.from_buffer(dptr, p),
            ffi.from_buffer(dptr, out))
        return out[0]
    return call


_SPECIAL = [0.0, -0.0, 1e-300, -1e-300, 1e-12, 0.025, 0.05, 0.3, 0.7, 1.0,
            2.0, 3.5, 4.0, 5.5, 50.0, -50.0, 1e30, -1e30, np.inf, -np.inf,
            np.nan]


@pytest.mark.parametrize('kind', ['pnj', 'fet', 'vds', 'delta', 'id'])
def test_every_law_is_the_python_one_bit_for_bit(law, kind):
    """20 000 random points and every pair of special values (ties, signed
    zeros, NaN, infinities, the `2*VT` escape, `IS <= 0`, a denormal IS,
    `VT = 0`, `vold < 0` for vds, a NaN vto): the same bits (a NaN's sign
    bit aside, as the backend's own sweeps allow)."""
    rng = np.random.default_rng(hash(kind) & 0xffff)
    cases = []
    for _ in range(20000):
        vnew = rng.uniform(-6.0, 6.0) * 10.0 ** rng.integers(-3, 3)
        vold = vnew + rng.uniform(-3.0, 3.0) * 10.0 ** rng.integers(-3, 2)
        if rng.random() < 0.1:
            vold = vnew
        if kind == 'pnj':
            pars = (10.0 ** rng.uniform(-20, -8), rng.uniform(0.01, 0.05))
        elif kind == 'fet':
            pars = (rng.uniform(-2.0, 2.0),)
        elif kind == 'delta':
            pars = (10.0 ** rng.uniform(-2, 1),)
        else:
            pars = ()
        cases.append((vnew, vold, pars))
    for a in _SPECIAL:
        for b in _SPECIAL:
            if kind == 'pnj':
                for IS, VT in ((1e-14, 0.025), (0.0, 0.025), (-1e-14, 0.025),
                               (5e-324, 0.025), (1e-14, 0.0), (1e-14, np.nan),
                               (1e3, 0.025)):
                    cases.append((a, b, (IS, VT)))
            elif kind == 'fet':
                for vto in (0.7, -0.7, 0.0, np.nan):
                    cases.append((a, b, (vto,)))
            elif kind == 'delta':
                for vmax in (1.0, 1e-3, np.inf):
                    cases.append((a, b, (vmax,)))
            else:
                cases.append((a, b, ()))
    tally = {'equal': 0, 'nan-bits': 0, 'zero-sign': 0, 'value': 0}
    for vnew, vold, pars in cases:
        with warnings.catch_warnings(), np.errstate(all='ignore'):
            warnings.simplefilter('ignore', RuntimeWarning)
            ref = _limiting.apply_limit(kind, vnew, vold, list(pars), numeric)
        p0 = pars[0] if len(pars) > 0 else 0.0
        p1 = pars[1] if len(pars) > 1 else 0.0
        got = law(kind, vnew, vold, p0, p1)
        tally[_compare(np.array([float(ref)]), np.array([got]))] += 1
    assert tally['value'] == 0 and tally['zero-sign'] == 0, tally


## -- the whole limiter, class by class --------------------------------------

def _cases(rng, n, count):
    for _ in range(count):
        scale = rng.choice([1e-3, 0.1, 1.0, 10.0, 100.0])
        x0 = rng.uniform(-2.0, 2.0, n)
        x = x0 + scale * rng.standard_normal(n)
        k = rng.integers(0, 5)
        if k == 1:
            x = x0.copy()                               # a tie everywhere
        elif k == 2:
            x[rng.integers(0, n)] = -0.0                # a signed zero
            x0[rng.integers(0, n)] = 0.0
        elif k == 3:
            x[:2] = 50.0 * np.sign(x[:2])               # the rails
        elif k == 4:
            x[rng.integers(0, n)] = x0[rng.integers(0, n)]
        yield x, x0


@pytest.mark.parametrize('name', _limiting_classes())
def test_the_kernel_answers_the_closure_bit_for_bit(name):
    """Every library class with `$limit`: 2000 random states and steps of
    every size, ties, signed zeros, rails, a parameter change and a
    temperature change on the way -- the kernel and the closure give the
    same bytes, the kernel's a fresh array, the inputs untouched."""
    cls = getattr(eh, name)
    e = _inst(cls)
    info = _bound(e)
    kern = cl.bind(type(e), info)
    assert kern is not None, info.get('_c_limit_status')
    n = e.n
    rng = np.random.default_rng(len(name))
    epars = [defaultepar, types.SimpleNamespace(T=350.0),
             types.SimpleNamespace(T=300), types.SimpleNamespace(T=np.float64(280.5))]
    numeric_pars = [nm for nm in e.ipar._paramnames
                    if isinstance(getattr(e.ipar, nm), (int, float))
                    and not isinstance(getattr(e.ipar, nm), bool)
                    and getattr(e.ipar, nm) not in (0, 0.0)]
    for k, (x, x0) in enumerate(_cases(rng, n, 2000)):
        if k % 500 == 250 and numeric_pars:
            nm = numeric_pars[(k // 500) % len(numeric_pars)]
            setattr(e.ipar, nm, getattr(e.ipar, nm) * 1.01)
            e.update_iparv()
        epar = epars[k % len(epars)]
        xc, x0c = x.copy(), x0.copy()
        a = _py_limit(e, x, x0, epar)
        b = kern(e, x, x0, epar)
        assert b is not None, (name, k)
        assert a.tobytes() == b.tobytes(), (name, k, x, x0)
        assert not np.shares_memory(b, x) and b.dtype == np.float64
        assert x.tobytes() == xc.tobytes() and x0.tobytes() == x0c.tobytes()


def test_the_inputs_the_closure_accepts_are_served():
    """A list, an int vector, a vector longer than the class reads, the
    same array as `x` and `x0` (nothing bites: equal and fresh)."""
    e = _inst(eh.MosLevel1Hdl)
    info = _bound(e)
    kern = cl.bind(type(e), info)
    rng = np.random.default_rng(5)
    n = e.n
    x0 = rng.uniform(-1.0, 1.0, n)
    x = x0 + rng.standard_normal(n)
    for xa, x0a in ((list(x), x0), (x, list(x0)), (np.round(x * 3).astype(int), x0),
                    (np.concatenate([x, [7.0, -7.0]]), np.concatenate([x0, [1.0, 2.0]])),
                    (x, x), (x0, x0)):
        a = _py_limit(e, xa, x0a)
        b = kern(e, xa, x0a, defaultepar)
        assert b is not None and a.tobytes() == b.tobytes()
    same = kern(e, x0, x0, defaultepar)
    assert np.array_equal(same, x0) and not np.shares_memory(same, x0)
    ## what it does not serve: a 0-d state, a short one, a 2-D one
    assert kern(e, np.array(1.0), x0, defaultepar) is None
    assert kern(e, x[:n - 1], x0, defaultepar) is None
    assert kern(e, x.reshape(1, -1), x0, defaultepar) is None


def test_a_temperature_that_is_not_one_number_is_the_closures():
    e = _inst(eh.MosLevel1Hdl)
    info = _bound(e)
    kern = cl.bind(type(e), info)
    x0 = np.linspace(-0.5, 1.0, e.n)
    x = x0 + 0.8
    for T in (310.0, 310, np.float64(310.0)):
        epar = types.SimpleNamespace(T=T)
        assert kern(e, x, x0, epar).tobytes() == _py_limit(e, x, x0, epar).tobytes()
    assert kern(e, x, x0, types.SimpleNamespace(T=np.array(310.0))) is None
    assert kern(e, x, x0, types.SimpleNamespace(T=np.array([300.0, 310.0]))) is None


def test_a_nan_or_infinite_state_is_declined():
    """A NaN at a probe terminal makes a ranking key NaN, whose order only
    Python's own sort knows: the kernel declines (None) and the closure
    answers; an infinity likewise where it reaches a key."""
    e = _inst(eh.MosLevel1Hdl)
    info = _bound(e)
    kern = cl.bind(type(e), info)
    x0 = np.linspace(-0.5, 1.0, e.n)
    for bad in (np.nan, np.inf, -np.inf, 1e308):
        for row in range(e.n):
            x = x0 + 0.3
            x[row] = bad
            got = kern(e, x, x0, defaultepar)
            if got is not None:
                ## (no key reached: the same bits as the closure)
                assert got.tobytes() == _py_limit(e, x, x0).tobytes(), (bad, row)
            elif not np.isnan(bad):
                continue
            else:
                assert got is None
    x = x0 + 0.3
    x[1] = np.nan                                   # the gate row
    assert kern(e, x, x0, defaultepar) is None


def test_a_sequential_group_follows_spice_order():
    """A chained `limit_together(..., sequential=True)` model (none in the
    library): the declaration order, the minus terminal moving."""
    cls = _fet('seq', chained=True)
    e = cls('d', 'g', 's')
    e.update_iparv()
    info = _bound(e)
    assert info['limit_groups'] and info['limit_groups'][0][0] is True
    kern = cl.bind(type(e), info)
    assert kern is not None, info.get('_c_limit_status')
    rng = np.random.default_rng(9)
    for k, (x, x0) in enumerate(_cases(rng, e.n, 800)):
        a = _py_limit(e, x, x0)
        b = kern(e, x, x0, defaultepar)
        assert b is not None and a.tobytes() == b.tobytes(), k


def test_a_wrong_law_is_caught(monkeypatch):
    """The sweep can fail: a prelude whose `fetlim` clamps at `vto + 4.5`
    instead of `vto + 4.0` gives other bytes on MosLevel1."""
    e = _inst(eh.MosLevel1Hdl)
    info = _bound(e)
    wrong = cl._LIMIT_C.replace('vto + 4.0', 'vto + 4.5')
    assert wrong != cl._LIMIT_C
    monkeypatch.setattr(cl, '_LIMIT_C', wrong)
    kern = cl.bind(type(e), dict(info))
    assert kern is not None
    rng = np.random.default_rng(1)
    differ = 0
    for x, x0 in _cases(rng, e.n, 400):
        if _py_limit(e, x, x0).tobytes() != kern(e, x, x0, defaultepar).tobytes():
            differ += 1
    assert differ > 0


## -- binding, refusal, the cache ---------------------------------------------

def test_a_limiter_that_cannot_be_printed_refuses_itself_only(monkeypatch):
    e = _inst(eh.GummelPoonNpnHdl)
    info = _bound(e)
    status = type(e)._hdl_backend_status
    monkeypatch.setattr(cl, 'render', lambda info_: (_ for _ in ()).throw(
        cl.Refused('a node without a C rendering')))
    assert cl.bind(type(e), info) is None
    assert info['_c_limit'] is None
    assert info['_c_limit_status'] == 'numpy (a node without a C rendering)'
    assert type(e)._hdl_backend_status == status
    monkeypatch.undo()
    assert cl.bind(type(e), info) is not None and info['_c_limit_status'] == 'c'
    cl.unbind(info)
    assert info['_c_limit'] is None and info['_c_limit_status'].startswith('numpy')
    ## (bound again: other tests in this process read the class)
    assert cl.bind(type(e), info) is not None


def test_a_class_with_no_limit_or_no_c_source_is_refused():
    info = {'limit_spec': [], 'funcs': {}}
    with pytest.raises(cl.Refused):
        cl.render(info)
    e = _inst(eh.MosLevel1Hdl)
    info = dict(_bound(e))
    info['funcs'] = {'i': types.SimpleNamespace(_csrc=None)}
    with pytest.raises(cl.Refused):
        cl.render(info)


def test_the_object_is_keyed_by_the_limiters_prelude_too(monkeypatch):
    e = _inst(eh.DiodeSpiceHdl)
    info = _bound(e)
    csrc, _nx, _layout = cl.render(info)
    key = cb.source_key(csrc)
    monkeypatch.setattr(cl, '_LIMIT_C', cl._LIMIT_C + '\n/* moved */\n')
    assert cb.source_key(cl.render(info)[0]) != key


def test_a_bound_limiter_freezes_and_thaws_without_its_kernel():
    e = _inst(eh.EkvNmosHdl)
    info = _bound(e)
    assert cl.bind(type(e), info) is not None
    payload = _hdl_cache.freeze(info)
    pickle.dumps(payload)
    back = _hdl_cache.thaw(payload)
    assert '_c_limit' not in back and '_c_limit_status' not in back
    assert back['limit_spec'][0][3][0].__dict__.get('_hdl_limit_par') is not None


## -- bound: the closure takes the kernel (LC2, 2026-10-03) ------------------

@pytest.mark.parametrize('name', _limiting_classes())
def test_every_limiting_library_class_binds_its_kernel(name):
    """A C-bound class binds its limiter with its chain functions, and
    `explain` says so on the backend line (the digests leave that line
    out)."""
    e = _inst(getattr(eh, name))
    info = _bound(e)
    assert info.get('_c_limit') is not None, info.get('_c_limit_status')
    assert info['_c_limit_status'] == 'c'
    line = next(l_ for l_ in hdl.explain(type(e)).splitlines()
                if l_.startswith('backend:'))
    assert line.startswith('backend: c') and line.endswith('(limit: c)'), line


def test_the_closure_takes_the_kernel_and_a_numpy_pin_restores_the_closure():
    e = _inst(eh.MosLevel1Hdl)
    info = _bound(e)
    kern = info['_c_limit']
    calls = []

    def counting(el, x, x0, ep):
        calls.append(1)
        return kern(el, x, x0, ep)
    info['_c_limit'] = counting
    x0 = np.linspace(-0.5, 1.0, e.n)
    x = x0 + 0.8
    out = e.limit(x, x0)
    assert calls == [1]
    assert out.tobytes() == kern(e, x, x0, defaultepar).tobytes()
    info['_c_limit'] = kern
    hdl.set_backend('numpy', type(e))
    try:
        assert type(e)._hdl_info.get('_c_limit') is None
        assert _py_limit(e, x, x0).tobytes() == out.tobytes()
    finally:
        hdl.set_backend(None, type(e))
    ## a resolved class resolves again at once under 'auto': re-bound
    assert type(e)._hdl_info.get('_c_limit') is not None


def test_the_switch_keeps_the_closure(monkeypatch):
    """`ENABLED` off: `limit()` never consults the kernel (the chain
    functions stay on C); on again, it does."""
    e = _inst(eh.MosLevel1Hdl)
    info = _bound(e)
    kern = info['_c_limit']
    calls = []

    def counting(el, x, x0, ep):
        calls.append(1)
        return kern(el, x, x0, ep)
    info['_c_limit'] = counting
    try:
        x0 = np.linspace(-0.5, 1.0, e.n)
        monkeypatch.setattr(cl, 'ENABLED', False)
        a = _py_limit(e, x0 + 0.8, x0)
        assert calls == []
        monkeypatch.setattr(cl, 'ENABLED', True)
        b = e.limit(x0 + 0.8, x0)
        assert calls == [1] and a.tobytes() == b.tobytes()
    finally:
        info['_c_limit'] = kern


def test_a_declined_call_is_answered_by_the_closure():
    """A NaN at a probe terminal: the kernel declines, `limit()` falls
    through to the closure, whose answer is the numpy-pinned one."""
    e = _inst(eh.MosLevel1Hdl)
    info = _bound(e)
    x0 = np.linspace(-0.5, 1.0, e.n)
    x = x0 + 0.3
    x[1] = np.nan
    assert info['_c_limit'](e, x, x0, defaultepar) is None
    out = _py_limit(e, x, x0)
    with numpy_backend(type(e)):
        ref = _py_limit(e, x, x0)
    assert out.tobytes() == ref.tobytes()


def test_an_instance_on_the_jax_toolkit_never_consults_the_kernel():
    pytest.importorskip('jax')
    from pycircuit.circuit.toolkit import jaxtoolkit
    ref = _inst(eh.EkvNmosHdl)
    info = _bound(ref)
    kern = info['_c_limit']
    calls = []

    def counting(el, x, x0, ep):
        calls.append(1)
        return kern(el, x, x0, ep)
    info['_c_limit'] = counting
    try:
        e = eh.EkvNmosHdl(*[Node(f'n{k}') for k in range(4)], toolkit=jaxtoolkit)
        e.update_iparv()
        x0 = np.linspace(-0.5, 1.0, 4)
        with contextlib.suppress(Exception):      # the closure's business there
            e.limit(x0 + 0.8, x0)
        assert calls == []
        ref.limit(x0 + 0.8, x0)
        assert calls == [1]
    finally:
        info['_c_limit'] = kern


def test_a_thawed_class_binds_its_kernel_again():
    e = _inst(eh.EkvNmosHdl)
    info = _bound(e)
    back = _hdl_cache.thaw(_hdl_cache.freeze(info))
    assert back.get('_c_limit') is None
    cb.attach(type(e), back)
    cb.ensure(type(e), back, numeric)
    assert back.get('_c_limit') is not None and back['_c_limit_status'] == 'c'
    x0 = np.linspace(-0.5, 1.0, e.n)
    assert back['_c_limit'](e, x0 + 0.8, x0, defaultepar).tobytes() == \
        info['_c_limit'](e, x0 + 0.8, x0, defaultepar).tobytes()
