"""A numeric run of cached models never imports sympy.

The modules a numeric run imports bind sympy and the HDL language's sympy
names as `_lazy` stand-ins, settled the moment sympy is imported; the compile
cache keeps an entry's sympy trees apart and loads them on first use
(`_hdl_cache._defer`).  Importing sympy was ~0.15 s of a short deck's whole
process (rca 0.55 -> 0.35 s, 2026-11-04).
"""
import json
import os
import pickle
import subprocess
import sys
import textwrap

import pytest
import sympy

from pycircuit.circuit import _hdl_cache as hc
from pycircuit.circuit import _lazy

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def _python(code):
    """`code`'s output, run by a fresh interpreter from the repository root."""
    out = subprocess.run([sys.executable, '-c', textwrap.dedent(code)], cwd=ROOT,
                         capture_output=True, text=True, timeout=900, check=False)
    assert out.returncode == 0, out.stderr[-4000:]
    return out.stdout.strip().splitlines()


#: What the two tests below import: `compact` compiles its classes as it is
#: imported -- served from the cache, once a first process has made them.
_IMPORTS = 'from pycircuit.circuit import hdl, elements_hdl, toolkit, compact, _lazy'


@pytest.fixture(scope='module')
def warm():
    """The cache holding the imports' classes (a first process compiles
    what is missing -- after an edit to the compiler, all of them)."""
    _python(_IMPORTS)


_DECKS = '''
    import sys, warnings, logging
    warnings.simplefilter('ignore'); logging.disable(logging.WARNING)
    sys.path.insert(0, 'benchmarks')
    import spice_suite
    from pycircuit.circuit import spice_import, circuit, integrator
    assert 'sympy' not in sys.modules
    circuit.default_toolkit = circuit.numeric
    for deck in ('rca', 'g1310', 'mosrect'):
        for radau in (False, True):
            imp = spice_import.import_netlist(
                spice_suite.benchdata.spice_data(spice_suite.CASES[deck][0]))
            tr, kw = imp.transient(**spice_suite.OVERRIDES.get(deck, {}))
            if radau:
                tr.par.integrator = integrator.RadauIIA3Integrator()
            tr.solve(**kw)
    print('sympy' in sys.modules)
'''


def test_a_warm_run_of_three_decks_imports_no_sympy():
    """Three gold decks -- a Gummel-Poon bipolar, MOS with diodes, a MOS
    rectifier -- imported and run under gear2 and Radau by a second process,
    every class served from the cache (the first may compile): sympy is
    never imported."""
    _python(_DECKS)
    assert _python(_DECKS)[-1] == 'False'


def test_sympy_s_import_settles_every_standin(warm):
    """Before sympy is imported the language's names, sympy and the symbolic
    backend are stand-ins; its import -- any import of it -- puts each
    target in its place in every module's namespace, leaves none pending and
    sympy its own loader, and takes the import watch down."""
    lines = _python('''
        import sys, importlib
        from pycircuit.circuit import hdl, elements_hdl, toolkit, compact, _lazy
        assert 'sympy' not in sys.modules
        names = [(hdl, 'ddt'), (hdl, 'TEMP'), (hdl, 'sympy'), (hdl, 'Quantity'), (hdl, '_CLAMPABLE'),
                 (elements_hdl, 'TEMP'), (elements_hdl, 'sympy'), (elements_hdl, '_BS_V'),
                 (toolkit, 'sympy'), (toolkit, '_symbolic'), (toolkit, 'DomainMatrix'), (compact, 'sympy')]
        assert all(_lazy.is_standin(getattr(m, n)) for m, n in names)
        assert isinstance(3, hdl.Quantity) is False and 'sympy' not in sys.modules
        import sympy
        from pycircuit.circuit import _hdl_sym, _symbolic
        left = [(k, n) for k, m in list(sys.modules.items()) if k.startswith('pycircuit')
                for n, v in vars(m).items() if _lazy.is_standin(v)]
        assert left == [], left
        assert hdl.ddt is _hdl_sym.ddt and hdl.TEMP is _hdl_sym.TEMP is elements_hdl.TEMP
        assert hdl.sympy is sympy is toolkit.sympy and toolkit._symbolic is _symbolic
        assert elements_hdl._BS_V is sympy.Symbol('v_ctrl', real=True)
        assert _lazy._PENDING == [] and _lazy._WATCH[0] is None
        assert type(sympy.__loader__).__name__ != '_Then'
        assert not any(type(f).__name__ == '_Watch' for f in sys.meta_path)
        print('ok')
    ''')
    assert lines[-1] == 'ok'


def test_every_standin_is_keyed_as_its_target(warm):
    """The cache key sees a stand-in as its target -- `_value_sig`'s text and
    the module its dependency hashes name -- so a class keyed where sympy is
    not imported and where it is has one key: every stand-in the library
    binds, made in a process without sympy, against its object here."""
    lines = _python('''
        import sys, json
        from pycircuit.circuit import hdl, elements_hdl, toolkit, compact, _lazy
        assert 'sympy' not in sys.modules
        out = []
        for s in _lazy._PENDING:
            if type(s) is _lazy.LazyModule:
                out.append(['module', s.__name__, _lazy.describe(s)])
            else:
                out.append(['standin', list(s._target), _lazy.describe(s)])
        print(json.dumps(out))
    ''')
    made = json.loads(lines[-1])
    assert len(made) >= 40
    for kind, target, (sig, dep) in made:
        if kind == 'module':
            real = __import__(target, fromlist=['_'])
        else:
            real = _lazy._value_of(tuple(target))
        assert hc._value_sig('x', real, 0, set(), False) == sig, target
        modname = getattr(real, '__module__', None)
        want = modname if isinstance(modname, str) and modname.startswith('pycircuit') else None
        assert dep == want, target


def _fresh(target):
    """A stand-in for `target`, unresolved (sympy is imported here, so
    `standin` would hand back the object itself)."""
    return _lazy.Standin(target, 'test', None, None)


def test_a_standin_forwards_as_its_target():
    """A stand-in captured before sympy's import (a closure's, a default's)
    answers as its target: arithmetic both ways, comparison, hashing,
    `sympify`, calls, attributes, `isinstance`, pickling -- the same objects
    the target gives."""
    pi, x = sympy.pi, sympy.Symbol('x')
    s = _fresh(('sympy', 'pi', False))
    for got, want in [(2 * s, 2 * pi), (s * 2, pi * 2), (s + x, pi + x), (x - s, x - pi),
                      (s / 3, pi / 3), (3 / s, 3 / pi), (s ** 2, pi ** 2), (2 ** s, 2 ** pi),
                      (-s, -pi), (abs(s), abs(pi)), (float(s), float(pi))]:
        assert got == want and type(got) is type(want)
    assert s == pi and hash(s) == hash(pi) and sympy.sympify(s) is pi
    assert sympy.srepr(sympy.sin(s * x)) == sympy.srepr(sympy.sin(pi * x))
    assert pickle.loads(pickle.dumps(s)) is pi and repr(s) == repr(pi)
    sym = _fresh(('sympy', 'Symbol', False))
    assert sym('y', real=True) is sympy.Symbol('y', real=True)
    assert isinstance(x, sym) and not isinstance(3, sym) and issubclass(sympy.Dummy, sym)
    assert sym.is_Symbol is True
    ## (a dunder is never forwarded: numpy's and copy's probes do not import sympy)
    with pytest.raises(AttributeError):
        _ = sym.__array_interface__


def test_resolve_function_puts_the_targets_in_a_closure():
    """`resolve_function` (the compiler's, before a body runs): the stand-ins
    in a function's closure cells and defaults -- and in those of a function
    a cell holds -- become their targets."""
    s, t = _fresh(('sympy', 'pi', False)), _fresh(('sympy', 'E', False))

    def outer(a=s):
        def inner():
            return t
        def body():
            return a, inner()
        return body

    body = outer()
    _lazy.resolve_function(body)
    assert body() == (sympy.pi, sympy.E)
    assert body.__closure__[0].cell_contents is sympy.pi
    assert body.__closure__[1].cell_contents.__closure__[0].cell_contents is sympy.E
    f = (lambda k=s: k)
    _lazy.resolve_function(f)
    assert f.__defaults__[0] is sympy.pi


def test_the_numpy_namespace_is_lambdify_s():
    """The namespace a lambdified record is thawed into, built without
    sympy, is lambdify's numpy namespace: the same names bound to the same
    objects (``I``, a complex, equal)."""
    import importlib
    lam = importlib.import_module('sympy.utilities.lambdify')
    lam._import('numpy')
    theirs, ours = lam.MODULES['numpy'][0], hc._numpy_namespace()
    assert set(theirs) == set(ours)
    assert [k for k in theirs if theirs[k] is not ours[k]] == ['I'] and ours['I'] == theirs['I']


def test_the_sympy_version_is_read_without_sympy(monkeypatch):
    """The key's sympy version, where sympy is not imported, read from its
    `release.py`: sympy's own."""
    want = sympy.__version__
    monkeypatch.setattr(hc, '_SYMPY_VERSION', [])
    monkeypatch.delitem(sys.modules, 'sympy')
    assert hc._sympy_version() == want and hc._SYMPY_VERSION == [want]


@pytest.fixture
def cold_library(tmp_path, monkeypatch):
    """Four library classes compiled cold into an empty cache -- a vector
    PCNR diode, MOS level 1, Gummel-Poon, EKV: `(cls, payload)` each, the
    cold info frozen and its trees apart as the cache stores it."""
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.tests.test_hdl_cache import _reclass
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path / 'hdlcache'))
    monkeypatch.delenv('PYCIRCUIT_HDL_CACHE', raising=False)
    monkeypatch.setattr(hc, 'ENABLED', True)
    out = []
    for base in (eh.DiodeSpiceHdl, eh.MosLevel1Hdl, eh.GummelPoonNpnHdl, eh.EkvNmosHdl):
        cls = _reclass(base)
        assert cls._hdl_cache_status == 'miss'
        out.append((cls, hc._defer(hc.freeze(cls._hdl_info))))
    return out


def _deferred(obj, out):
    """The `_Deferred` dicts in `obj` (its own keys: nothing loaded)."""
    if type(obj) is hc._Deferred:
        out.append(obj)
    if isinstance(obj, dict):
        for v in dict.values(obj):
            _deferred(v, out)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            _deferred(v, out)
    return out


def test_a_payload_names_no_sympy_class_outside_its_trees(cold_library):
    """A library entry's payload apart from its trees (`_TREES`) holds no
    sympy object -- nothing its pickle names comes from sympy or the
    language's sympy classes -- and thawed, its trees load on the first
    read of one, equal to the cold compile's (srepr)."""
    for cls, payload in cold_library:
        assert hc._TREES in payload, cls.__name__
        main = {k: v for k, v in payload.items() if k != hc._TREES}
        assert not hc._needs_sympy(main), cls.__name__
        blob = pickle.dumps(main, protocol=pickle.HIGHEST_PROTOCOL)
        assert b'sympy' not in blob and b'_hdl_sym' not in blob, cls.__name__
        info = hc._thaw_payload(payload)
        cold = cls._hdl_info
        for k in ('sym_spec', 'pure_spec', 'pcnr_affine', 'pcnr_funcs', 'pcnr_vector'):
            assert (cold.get(k) is None) == (info.get(k) is None), (cls.__name__, k)
        found = _deferred(info, [])
        lazies = [f._hdl_limit_par for f in _limit_fns(info, [])]
        assert found or lazies, cls.__name__
        trees = found[0]._trees if found else lazies[0]._trees
        assert trees._vals is None
        for k in ('sym_spec', 'pure_spec', 'pcnr_affine'):
            if cold.get(k) is not None:
                assert sympy.srepr(dict(info[k].items())) == sympy.srepr(cold[k]), (cls.__name__, k)
        for got, want in zip(info.get('pcnr_funcs') or (), cold.get('pcnr_funcs') or ()):
            for kk in ('sym', 'chain'):
                assert sympy.srepr(got.get(kk)) == sympy.srepr(want.get(kk)), (cls.__name__, kk)
        assert all(dict(d.items()) is not None for d in found)
        assert all(len(tuple(lp)) == 5 for lp in lazies) and trees._vals is not None


def _limit_fns(obj, out):
    """The functions carrying `_hdl_limit_par` in `obj`, in order -- a
    `_Deferred`'s own keys only (no function is deferred), nothing loaded."""
    if isinstance(obj, dict):
        for v in dict.values(obj):
            _limit_fns(v, out)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            _limit_fns(v, out)
    elif callable(obj) and '_hdl_limit_par' in getattr(obj, '__dict__', {}):
        out.append(obj)
    return out


def test_a_thawed_limiter_parameter_loads_on_its_first_use(cold_library):
    """A thawed function's `_hdl_limit_par` (the limiter parameter's sympy
    ingredients) is kept among the trees: unloaded until it is iterated or
    indexed, then the cold compile's tuple (`unlazy`); frozen again, the
    tuple."""
    seen = 0
    for cls, payload in cold_library:
        info = hc._thaw_payload(payload)
        cold, warm = _limit_fns(cls._hdl_info, []), _limit_fns(info, [])
        assert len(cold) == len(warm), cls.__name__
        for c, w in zip(cold, warm):
            lp = w._hdl_limit_par
            assert type(lp) is hc._LazyTuple and type(c._hdl_limit_par) is tuple
            assert lp._trees._vals is None or seen
            assert len(lp) == 5 and sympy.srepr(tuple(lp)) == sympy.srepr(c._hdl_limit_par)
            assert type(hc.unlazy(lp)) is tuple
            assert sympy.srepr(hc.freeze(w)[4]) == sympy.srepr(c._hdl_limit_par)
            seen += 1
    assert seen > 0
