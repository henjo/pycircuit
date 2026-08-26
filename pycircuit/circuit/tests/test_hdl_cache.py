"""The on-disk compile cache for `Behavioural` elements (`_hdl_cache`).

The claim under test is BIT-IDENTITY: a class served from the cache
must evaluate `i`, `G`, `q`, `C`, `CY` to exactly the bytes a cold
compile produces, and `explain()` must print the same text.  Every
element in `elements_hdl` and a set of models covering the compiler's
paths (eager, chained, collapse, `$limit`, PCNR, `idtmod`, `@cross`,
`laplace`) is compiled cold, saved, rebuilt from the file and compared.

Then the things a cache gets wrong when it is wrong: it must MISS when
the analog source, a parameter default, the format number or the
compiler source changes; it must survive a corrupt, truncated or foreign
file; two processes writing the same key must not tear it; and a class
whose source cannot be recovered must compile and say so.
"""

import importlib.util
import os
import pickle
import subprocess
import sys
import textwrap
import types

import numpy as np
import pytest
import sympy

from pycircuit.circuit import _hdl_cache as hc
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import Node
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution, Node as HNode,
                                   Collapse, Cross, ddt, idtmod, var,
                                   limit_pnj, laplace_zp, expl)
from pycircuit.utilities.param import Parameter


## ----------------------------------------------------------------------
## Fixtures and helpers.

@pytest.fixture
def cache_dir(tmp_path, monkeypatch):
    """A fresh, empty cache directory, and the cache switched on."""
    d = tmp_path / 'hdlcache'
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(d))
    monkeypatch.delenv('PYCIRCUIT_HDL_CACHE', raising=False)
    monkeypatch.setattr(hc, 'ENABLED', True)
    return d


def _entries(d):
    return sorted(p for p in os.listdir(d) if p.endswith('.pkl')) \
        if os.path.isdir(d) else []


def _reclass(base, name=None):
    """A new class with `base`'s analog and parameters -- the same key,
    compiled again under whatever cache state is current."""
    return type(name or base.__name__, (Behavioural,), dict(
        analog=staticmethod(base.__dict__['analog'].__func__
                            if isinstance(base.__dict__.get('analog'),
                                          staticmethod)
                            else base.analog),
        instparams=list(base.instparams),
        ## A `params_as` class hands analog() its namespace first; a
        ## reclass that dropped it would make `p` a terminal.
        params_as=getattr(base, 'params_as', None),
        __module__=base.__module__,
        __qualname__=base.__qualname__))


def _instance(cls, **params):
    e = cls(*[Node('n%d' % k) for k in range(len(cls.terminals))], **params)
    e.update_iparv()
    return e


def _evaluate(e, rng, npoints=4):
    """`i, G, q, C, CY` at random points, as raw bytes."""
    out = []
    for _ in range(npoints):
        x = rng.uniform(-0.8, 0.8, e.n)
        with np.errstate(all='ignore'):
            out.append(tuple(np.asarray(v).tobytes() for v in (
                e.i(x), e.G(x), e.q(x), e.C(x), e.CY(x, 2 * np.pi * 1e6))))
    return out


_SURFACE = ('limit', 'pcnr_i', 'pcnr_didv', 'pcnr_limit', 'eval_i_pure',
            'eval_q_pure', 'state_ic', 'periodic_states', 'linear',
            'IC_KIND', 'terminals', 'pcnr_junctions')


def _surface(cls):
    return {a: (getattr(cls, a, None) is not None
                if callable(getattr(cls, a, None))
                else getattr(cls, a, None)) for a in _SURFACE}


def _assert_identical(cold, warm, **params):
    """The whole bit-identity claim for one class pair."""
    assert cold._hdl_cache_status == 'miss', cold._hdl_cache_status
    assert warm._hdl_cache_status == 'hit', warm._hdl_cache_status
    assert _surface(cold) == _surface(warm)
    ec, ew = _instance(cold, **params), _instance(warm, **params)
    assert ec.n == ew.n
    assert type(ec)._hdl_cache_status in ('miss', 'hit')
    ## `explain()` reads the info dict and the generated source: the
    ## strongest cheap check that the record carried everything.
    assert hdl.explain(ec, maxlines=None) == hdl.explain(ew, maxlines=None)
    rng = np.random.default_rng(20260826)
    a = _evaluate(ec, np.random.default_rng(20260826))
    b = _evaluate(ew, rng)
    assert a == b
    if getattr(cold, 'limit', None) is not None:
        x = np.full(ec.n, 0.3)
        with np.errstate(all='ignore'):
            assert repr(ec.limit(x, x * 0.5)) == repr(ew.limit(x, x * 0.5))


## ----------------------------------------------------------------------
## Models covering the compiler's paths, defined here so `inspect` can
## see their source.

class EagerDiode(Behavioural):
    """Eager path, PCNR-eligible, `sym_spec` kept."""
    instparams = [Parameter(name='IS', default=1e-14),
                  Parameter(name='n', default=1.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        vt = hdl.vt(hdl.TEMP)
        return Contribution(b.I, IS * (expl(b.V / (n * vt)) - 1.0))


class ChainedLimited(Behavioural):
    """Chained path with `$limit` and a collapse."""
    instparams = [Parameter(name='IS', default=1e-14),
                  Parameter(name='rs', default=1.0),
                  Parameter(name='cj', default=1e-12)]

    @staticmethod
    def analog(plus, minus):
        i = HNode('i')
        bj = Branch(i, minus)
        bs = Branch(plus, i)
        vt = var(hdl.vt(hdl.TEMP), 'vt')
        vj = var(limit_pnj(bj.V, IS, vt), 'vj')
        ij = var(IS * (expl(vj / vt) - 1.0), 'ij')
        return (Contribution(bj.I, ij + ddt(cj * vj)),
                Contribution(bs.I, bs.V / rs),
                Collapse(bs, rs <= 0))


class WrappedState(Behavioural):
    """`idtmod` state with a DC pin, eager path."""
    instparams = [Parameter(name='k', default=2.0),
                  Parameter(name='m', default=1.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        ph = idtmod(k * b.V, 0.0, m)
        return Contribution(b.I, 1e-3 * b.V + 1e-4 * ph)


class Crossing(Behavioural):
    """`@cross` event, chained."""
    instparams = [Parameter(name='vth', default=0.5)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        d = var(b.V - vth, 'd')
        return (Contribution(b.I, 1e-3 * d), Cross(d))


class Ladder(Behavioural):
    """`laplace_zp`, eager path, states."""
    instparams = [Parameter(name='g', default=1e-3)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, g * laplace_zp(b.V, [], [-1e3, -1e4]))


LOCAL_MODELS = [EagerDiode, ChainedLimited, WrappedState, Crossing, Ladder]


## Imported at module level ON PURPOSE: importing the library compiles
## it, and that compile must not happen inside a test whose fixture has
## pointed the cache at its own directory -- the library's entries would
## land there and the "cold" reclass below would hit.
from pycircuit.circuit import elements_hdl as eh    # noqa: E402


def _library_models():
    return [getattr(eh, n) for n in sorted(dir(eh))
            if isinstance(getattr(eh, n), type)
            and hasattr(getattr(eh, n), '_hdl_info')
            and getattr(eh, n).__module__ == eh.__name__]


## ----------------------------------------------------------------------
## Bit-identity.

@pytest.mark.parametrize('model', LOCAL_MODELS, ids=lambda m: m.__name__)
def test_local_model_identical_from_cache(cache_dir, model):
    cold = _reclass(model)
    assert len(_entries(cache_dir)) == 1
    warm = _reclass(model)
    _assert_identical(cold, warm)


def test_every_library_element_identical_from_cache(cache_dir):
    """All of `elements_hdl`, cold then from the file, one process."""
    models = _library_models()
    assert len(models) >= 26
    for model in models:
        cold = _reclass(model)
        warm = _reclass(model)
        _assert_identical(cold, warm)
    ## One entry per class, plus one per collapse variant the default
    ## parameters selected when the instances above were built.
    assert len(_entries(cache_dir)) >= len(models)


def test_identical_across_processes(cache_dir, tmp_path):
    """The cold compile in one interpreter, the hit in another: the
    values the hit produces are the bytes the cold compile wrote."""
    script = textwrap.dedent('''
        import os, sys, pickle, numpy as np
        from pycircuit.circuit import elements_hdl as eh
        from pycircuit.circuit import hdl
        from pycircuit.circuit.circuit import Node
        out = {}
        rng = np.random.default_rng(7)
        for name in ('GummelPoonNpnHdl', 'DiodeHdl', 'MosLevel1Hdl',
                     'IdtmodHdl', 'RSkinHdl', 'MemristorHdl'):
            c = getattr(eh, name)
            e = c(*[Node('n%d' % k) for k in range(len(c.terminals))])
            e.update_iparv()
            x = rng.uniform(-0.5, 0.5, e.n)
            with np.errstate(all='ignore'):
                out[name] = (c._hdl_cache_status,
                             [np.asarray(v).tobytes() for v in
                              (e.i(x), e.G(x), e.q(x), e.C(x))],
                             hdl.explain(e, maxlines=None))
        pickle.dump(out, open(sys.argv[1], 'wb'))
    ''')
    env = dict(os.environ, PYCIRCUIT_HDL_CACHE_DIR=str(cache_dir),
               PYCIRCUIT_HDL_CACHE='1')
    outs = []
    for k in range(2):
        p = tmp_path / ('run%d.pkl' % k)
        subprocess.run([sys.executable, '-c', script, str(p)], check=True,
                       env=env, timeout=600)
        with open(p, 'rb') as fh:
            outs.append(pickle.load(fh))
    for name in outs[0]:
        assert outs[0][name][0] == 'miss'
        assert outs[1][name][0] == 'hit'
        assert outs[0][name][1:] == outs[1][name][1:]


def test_collapse_variant_is_cached(cache_dir):
    """`_collapse_variant`'s per-mask recompile goes through the cache
    with its own key, and the variant's stamps are identical."""
    base1 = _reclass(ChainedLimited)
    e1 = _instance(base1, rs=0.0)
    v1 = type(e1)
    assert v1 is not base1 and v1._hdl_cache_status == 'miss'
    base2 = _reclass(ChainedLimited)
    e2 = _instance(base2, rs=0.0)
    v2 = type(e2)
    assert v2._hdl_cache_status == 'hit'
    assert e1.n == e2.n == 2
    assert _evaluate(e1, np.random.default_rng(3)) == \
        _evaluate(e2, np.random.default_rng(3))
    assert len(_entries(cache_dir)) == 2


def test_symbolic_toolkit_served_from_cache(cache_dir):
    """`sym_spec` comes back as equal sympy expressions."""
    from pycircuit.circuit.toolkit import symbolic
    cold = _reclass(EagerDiode)
    warm = _reclass(EagerDiode)
    assert warm._hdl_cache_status == 'hit'
    xs = sympy.symbols('v0 v1')
    outs = []
    for cls in (cold, warm):
        e = _instance(cls)
        e.toolkit = symbolic
        outs.append((str(e.i(xs)), str(e.G(xs))))
    assert outs[0] == outs[1]
    assert 'v0' in outs[0][0]


## ----------------------------------------------------------------------
## Invalidation.

_MODULE = '''
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, var
from pycircuit.utilities.param import Parameter


class Probe(Behavioural):
    instparams = [Parameter(name='r', default=%(default)s)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        g = var(%(scale)s / r, 'g')
        return Contribution(b.I, b.V * g)
'''


def _load_probe(tmp_path, tag, default='1e3', scale='1.0'):
    """`Probe` from a module file written for this variant.  Same
    module name and file name every time, so only the CONTENT differs
    in the key."""
    d = tmp_path / tag
    d.mkdir()
    path = d / 'probemod.py'
    path.write_text(_MODULE % dict(default=default, scale=scale))
    spec = importlib.util.spec_from_file_location('probemod', str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['probemod'] = mod
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop('probemod', None)
    return mod.Probe


def test_same_source_hits(cache_dir, tmp_path):
    a = _load_probe(tmp_path, 'a')
    b = _load_probe(tmp_path, 'b')
    assert a._hdl_cache_status == 'miss'
    assert b._hdl_cache_status == 'hit'
    assert len(_entries(cache_dir)) == 1


def test_one_character_of_analog_misses(cache_dir, tmp_path):
    _load_probe(tmp_path, 'a', scale='1.0')
    c = _load_probe(tmp_path, 'c', scale='2.0')
    assert c._hdl_cache_status == 'miss'
    assert len(_entries(cache_dir)) == 2


def test_parameter_default_misses(cache_dir, tmp_path):
    _load_probe(tmp_path, 'a', default='1e3')
    c = _load_probe(tmp_path, 'c', default='1e4')
    assert c._hdl_cache_status == 'miss'
    assert len(_entries(cache_dir)) == 2


def test_format_bump_misses(cache_dir, tmp_path, monkeypatch):
    _load_probe(tmp_path, 'a')
    monkeypatch.setattr(hc, 'CACHE_FORMAT', hc.CACHE_FORMAT + 1)
    c = _load_probe(tmp_path, 'c')
    assert c._hdl_cache_status == 'miss'
    assert len(_entries(cache_dir)) == 2


def test_compiler_source_change_misses(cache_dir, tmp_path, monkeypatch):
    """The hdl.py hash is in the key; stand in for an edit."""
    _load_probe(tmp_path, 'a')
    monkeypatch.setattr(hc, '_compiler_hash', lambda: 'edited')
    c = _load_probe(tmp_path, 'c')
    assert c._hdl_cache_status == 'miss'


def test_helper_module_change_misses(cache_dir, tmp_path):
    """A helper the analog body calls by name lives in a module whose
    source is part of the key."""
    src = '''
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution
from pycircuit.utilities.param import Parameter


def gain():
    return %s


class Probe(Behavioural):
    instparams = [Parameter(name='r', default=1e3)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, gain() * b.V / r)
'''
    def load(tag, k):
        d = tmp_path / tag
        d.mkdir()
        path = d / 'probemod2.py'
        path.write_text(src % k)
        spec = importlib.util.spec_from_file_location('probemod2', str(path))
        mod = importlib.util.module_from_spec(spec)
        sys.modules['probemod2'] = mod
        try:
            spec.loader.exec_module(mod)
        finally:
            sys.modules.pop('probemod2', None)
        return mod.Probe
    assert load('a', '1.0')._hdl_cache_status == 'miss'
    assert load('b', '1.0')._hdl_cache_status == 'hit'
    c = load('c', '2.0')
    assert c._hdl_cache_status == 'miss'
    e = _instance(c)
    assert np.allclose(e.G(np.zeros(2))[0, 0], 2e-3)


def _closure_model(which, scale):
    """Same analog TEXT for every call; the body reads `which` and
    `scale` from this enclosing scope."""
    class Closed(Behavioural):
        instparams = [Parameter(name='r', default=1e3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            if which == 'chained':
                g = var(scale / r, 'g')
            else:
                g = scale / r
            return Contribution(b.I, b.V * g)
    return Closed


def test_closure_values_are_part_of_the_key(cache_dir):
    """The suite builds classes whose analog text is identical and whose
    behaviour is decided by an enclosing variable; those must not share
    an entry."""
    a = _closure_model('eager', 1.0)
    b = _closure_model('chained', 1.0)
    c = _closure_model('eager', 2.0)
    d = _closure_model('eager', 1.0)
    assert (a._hdl_cache_status, b._hdl_cache_status, c._hdl_cache_status,
            d._hdl_cache_status) == ('miss', 'miss', 'miss', 'hit')
    assert a._hdl_info['chained'] is False and b._hdl_info['chained']
    assert np.allclose(_instance(c).G(np.zeros(2))[0, 0], 2e-3)
    assert len(_entries(cache_dir)) == 3


_MODULE_CONST = """
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution
from pycircuit.utilities.param import Parameter

SCALE = %s


class Probe(Behavioural):
    instparams = [Parameter(name='r', default=1e3)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, SCALE * b.V / r)
"""


def test_module_constant_is_part_of_the_key(cache_dir, tmp_path):
    def load(tag, k):
        d = tmp_path / tag
        d.mkdir()
        path = d / 'probemod4.py'
        path.write_text(_MODULE_CONST % k)
        spec = importlib.util.spec_from_file_location('probemod4', str(path))
        mod = importlib.util.module_from_spec(spec)
        sys.modules['probemod4'] = mod
        try:
            spec.loader.exec_module(mod)
        finally:
            sys.modules.pop('probemod4', None)
        return mod.Probe
    assert load('a', '1.0')._hdl_cache_status == 'miss'
    assert load('b', '1.0')._hdl_cache_status == 'hit'
    assert load('c', '3.0')._hdl_cache_status == 'miss'


def test_closure_container_is_keyed_by_content(cache_dir):
    """A factory that closes over a list of poles is a legitimate model
    (`skin_effect_resistor`); the list's content is in the key."""
    def make(poles):
        class Poled(Behavioural):
            instparams = [Parameter(name='g', default=1e-3)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, g * laplace_zp(b.V, [], poles))
        return Poled
    a = make([-1e3, 0, -1e4, 0])          # (real, imag) pairs
    b = make([-1e3, 0, -1e4, 0])
    c = make([-1e3, 0, -1e4, 0, -1e5, 0])
    assert (a._hdl_cache_status, b._hdl_cache_status,
            c._hdl_cache_status) == ('miss', 'hit', 'miss')
    assert _instance(a).n == 4 and _instance(c).n == 5


def test_body_that_writes_outside_itself_is_uncacheable(cache_dir):
    """A body that writes into a dict it closes over is instrumenting
    itself -- a diagnostics test does, to see the parameter injection --
    and a hit would never run the body.  Refused, with the name."""
    seen = {}

    class Peeks(Behavioural):
        instparams = [Parameter(name='r', default=1e3)]

        @staticmethod
        def analog(plus, minus):
            seen['ran'] = True
            b = Branch(plus, minus)
            return Contribution(b.I, b.V / r)
    assert Peeks._hdl_cache_status == \
        'uncacheable: analog() writes into seen, which is not its own ' \
        'local -- a cache hit would never run the body, so that write ' \
        'would not happen'
    assert seen == {'ran': True}
    assert _entries(cache_dir) == []

    class Local(Behavioural):
        """Stores into its OWN dict: not a side effect, cacheable."""
        instparams = [Parameter(name='r', default=1e3)]

        @staticmethod
        def analog(plus, minus):
            d = {}
            d['b'] = Branch(plus, minus)
            return Contribution(d['b'].I, d['b'].V / r)
    assert Local._hdl_cache_status == 'miss'


def test_key_fields(cache_dir):
    """What the key is made of, so a change to that list is deliberate."""
    key = hc.key_for(EagerDiode)
    assert key and len(key) == 64
    assert hc.key_for(EagerDiode) == key
    assert hc.key_for(ChainedLimited) != key


## ----------------------------------------------------------------------
## Damage.

def test_corrupt_entry_is_a_miss_and_repaired(cache_dir, tmp_path):
    a = _load_probe(tmp_path, 'a')
    (name,) = _entries(cache_dir)
    path = cache_dir / name
    good = path.read_bytes()
    for damage in (b'garbage', good[:len(good) // 2], b''):
        path.write_bytes(damage)
        c = _load_probe(tmp_path, 'c%d' % len(damage))
        assert c._hdl_cache_status.startswith('miss (unreadable entry'), \
            c._hdl_cache_status
        assert path.read_bytes() == good, 'the entry was rewritten'
        e = _instance(c)
        assert np.allclose(e.G(np.zeros(2)), 1e-3 * np.array([[1, -1],
                                                              [-1, 1]]))
    del a


def test_foreign_entry_is_a_miss(cache_dir, tmp_path):
    _load_probe(tmp_path, 'a')
    (name,) = _entries(cache_dir)
    path = cache_dir / name
    path.write_bytes(pickle.dumps(dict(format=hc.CACHE_FORMAT,
                                       key='someone else', payload={})))
    c = _load_probe(tmp_path, 'c')
    assert c._hdl_cache_status == 'miss (foreign entry)'


def test_entry_that_will_not_rebuild_is_a_miss(cache_dir, tmp_path):
    _load_probe(tmp_path, 'a')
    (name,) = _entries(cache_dir)
    path = cache_dir / name
    entry = pickle.loads(path.read_bytes())
    entry['payload']['funcs']['i'] = (hc._FN, 999, ('lambdify', 'def f(:',
                                                     [], []), None)
    path.write_bytes(pickle.dumps(entry))
    c = _load_probe(tmp_path, 'c')
    assert c._hdl_cache_status.startswith('miss (entry would not rebuild')
    e = _instance(c)
    assert np.allclose(e.i(np.array([1.0, 0.0])), [1e-3, -1e-3])


def test_unwritable_directory_compiles_anyway(cache_dir, tmp_path,
                                              monkeypatch):
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR',
                       str(tmp_path / 'file_not_dir'))
    (tmp_path / 'file_not_dir').write_text('in the way')
    c = _load_probe(tmp_path, 'a')
    assert c._hdl_cache_status.startswith('unsaved:'), c._hdl_cache_status
    e = _instance(c)
    assert np.allclose(e.i(np.array([1.0, 0.0])), [1e-3, -1e-3])


def test_no_temp_files_left_behind(cache_dir, tmp_path):
    _load_probe(tmp_path, 'a')
    assert all(n.endswith('.pkl') for n in os.listdir(cache_dir))


## ----------------------------------------------------------------------
## Switches and the uncacheable.

def test_env_var_disables(cache_dir, tmp_path, monkeypatch):
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE', '0')
    c = _load_probe(tmp_path, 'a')
    assert c._hdl_cache_status == 'disabled'
    assert _entries(cache_dir) == []


def test_module_flag_disables(cache_dir, tmp_path, monkeypatch):
    monkeypatch.setattr(hc, 'ENABLED', False)
    c = _load_probe(tmp_path, 'a')
    assert c._hdl_cache_status == 'disabled'
    assert _entries(cache_dir) == []


def test_source_not_recoverable_is_uncacheable(cache_dir):
    ns = {}
    exec(textwrap.dedent('''
        from pycircuit.circuit.hdl import Behavioural, Branch, Contribution
        from pycircuit.utilities.param import Parameter
        class Dyn(Behavioural):
            instparams = [Parameter(name='r', default=1e3)]
            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, b.V / r)
    '''), ns)
    c = ns['Dyn']
    assert c._hdl_cache_status == 'uncacheable: analog() source not ' \
        'recoverable'
    assert _entries(cache_dir) == []
    e = _instance(c)
    assert np.allclose(e.i(np.array([1.0, 0.0])), [1e-3, -1e-3])


def test_debug_env_warns_on_uncacheable(cache_dir, monkeypatch):
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DEBUG', '1')
    with pytest.warns(UserWarning, match='hdl cache: Dyn2'):
        exec(textwrap.dedent('''
            from pycircuit.circuit.hdl import Behavioural, Branch, Contribution
            class Dyn2(Behavioural):
                instparams = []
                @staticmethod
                def analog(plus, minus):
                    b = Branch(plus, minus)
                    return Contribution(b.I, b.V)
        '''), {})


def test_clear(cache_dir, tmp_path):
    _load_probe(tmp_path, 'a')
    assert hc.clear() == 1
    assert _entries(cache_dir) == []
    assert hc.clear() == 0


## ----------------------------------------------------------------------
## Concurrency.

def test_two_processes_same_key(cache_dir, tmp_path):
    """Both compile and write the same key at once; a third reads a
    whole, valid entry."""
    d = tmp_path / 'mod'
    d.mkdir()
    (d / 'probemod3.py').write_text(_MODULE % dict(default='1e3',
                                                   scale='1.0'))
    script = textwrap.dedent('''
        import sys, importlib.util
        spec = importlib.util.spec_from_file_location('probemod3', sys.argv[1])
        mod = importlib.util.module_from_spec(spec)
        sys.modules['probemod3'] = mod
        spec.loader.exec_module(mod)
        print(mod.Probe._hdl_cache_status)
    ''')
    env = dict(os.environ, PYCIRCUIT_HDL_CACHE_DIR=str(cache_dir),
               PYCIRCUIT_HDL_CACHE='1')
    procs = [subprocess.Popen([sys.executable, '-c', script,
                               str(d / 'probemod3.py')],
                              env=env, stdout=subprocess.PIPE, text=True)
             for _ in range(2)]
    outs = [p.communicate(timeout=300)[0].strip() for p in procs]
    assert all(p.returncode == 0 for p in procs)
    assert all(o in ('miss', 'hit') for o in outs), outs
    (name,) = _entries(cache_dir)
    entry = pickle.loads((cache_dir / name).read_bytes())
    assert entry['format'] == hc.CACHE_FORMAT
    third = subprocess.run([sys.executable, '-c', script,
                            str(d / 'probemod3.py')], env=env,
                           capture_output=True, text=True, timeout=300)
    assert third.stdout.strip() == 'hit', third.stderr


## ----------------------------------------------------------------------
## The record mechanics.

def test_bytecode_scanner_matches_dis():
    """`_global_loads` is a raw bytecode scan; `dis` is the reference."""
    n = 0

    def walk(o):
        nonlocal n
        if isinstance(o, types.FunctionType):
            n += 1
            assert hc._global_loads(o.__code__) == \
                hc._global_loads_dis(o.__code__)
            inner = getattr(o, '_hdl_inner', None)
            if inner is not None:
                walk(inner)
        elif isinstance(o, dict):
            for v in o.values():
                walk(v)
        elif isinstance(o, (list, tuple)):
            for v in o:
                walk(v)
    for model in _library_models() + LOCAL_MODELS:
        walk(model._hdl_info)
    assert n > 200


def test_record_refuses_a_rebound_global(cache_dir):
    """A function whose namespace binds a name to something the
    rebuilt namespace would not is refused, not recorded wrongly."""
    ns = hdl._chain_namespace(None)
    ns['exp'] = lambda v: v          # not numpy's
    exec('def _f(x):\n    return [exp(x)]', ns)
    fn = ns['_f']
    fn._src = 'def _f(x):\n    return [exp(x)]'
    with pytest.raises(hc.Uncacheable, match='different object'):
        hc.freeze(dict(f=fn))


def test_freeze_thaw_round_trip_in_memory():
    info = EagerDiode._hdl_info
    back = hc.thaw(pickle.loads(pickle.dumps(hc.freeze(info))))
    assert back['funcs'].keys() == info['funcs'].keys()
    args = (np.array([0.3, 0.0]), 1e-14, 1.0, 300.0)
    with np.errstate(all='ignore'):
        for k in ('i', 'G', 'q', 'C'):
            assert np.asarray(back['funcs'][k](*args)).tobytes() == \
                np.asarray(info['funcs'][k](*args)).tobytes()
    assert back['sym_spec']['i'] == info['sym_spec']['i']
    assert [str(e) for e in back['sym_spec']['G']] == \
        [str(e) for e in info['sym_spec']['G']]


def test_wants_x_survives(cache_dir):
    """`limit()` reads `_wants_x` off each limiter parameter function."""
    cold = _reclass(eh.MosLevel1Hdl)
    warm = _reclass(eh.MosLevel1Hdl)
    assert warm._hdl_cache_status == 'hit'

    def flags(o):
        if isinstance(o, types.FunctionType):
            return [getattr(o, '_wants_x', 'absent')]
        if isinstance(o, dict):
            return [v for k in sorted(o) for v in flags(o[k])]
        if isinstance(o, (list, tuple)):
            return [v for x in o for v in flags(x)]
        return []
    assert flags(cold._hdl_info['limit_spec']) == \
        flags(warm._hdl_info['limit_spec'])
    assert set(flags(warm._hdl_info['limit_spec'])) & {True, False}


def test_dc_variants_share_when_nothing_is_pinned():
    """The pre-cache waste removal: with no DC-pinned state `i_dc`,
    `G_dc`, `u_dc` ARE `i`, `G`, `u`; with one they are distinct."""
    f = ChainedLimited._hdl_info['funcs']
    assert f['i_dc'] is f['i'] and f['G_dc'] is f['G'] and \
        f['u_dc'] is f['u']
    f = WrappedState._hdl_info['funcs']
    assert WrappedState._hdl_info['state_meta']['dc_pins']
    assert f['i_dc'] is not f['i'] and f['G_dc'] is not f['G']
