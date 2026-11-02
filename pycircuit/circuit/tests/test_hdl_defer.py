"""The deferred compile (`hdl.defer_compile`, roadmap 231).

A `Behavioural` class a module defines between `defer_compile(name)` and
`defer_compile(name, False)` -- the device library, `elements_hdl` -- is
created as any other and its declaration checked, but compiles on its first
use: an instance of it or of a subclass, or a read of any attribute the
compile sets (`hdl._COMPILED_ATTRS`), on the class or an instance.  The
claims: the library's import compiles nothing, nor does the netlist
importer's that names eighteen of its classes (each in a fresh interpreter:
this one compiled the library at the session's start); every way
in compiles the class and then reads what an eager compile of the same model
gives; no compile sets a name `_COMPILED_ATTRS` lacks, which would read as the
bases give it before the compile; a compile that raises holds the class again;
`PYCIRCUIT_HDL_DEFER=0` compiles at creation.
"""
import ast
import gc
import importlib.util
import inspect
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

import pycircuit
from pycircuit.circuit import _hdl_cache as hc
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import Node

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(pycircuit.__file__)))


def _child(script, **extra):
    env = {k: v for k, v in os.environ.items() if k not in ('PYCIRCUIT_HDL_DEFER', 'PYCIRCUIT_GC_FREEZE')}
    env.update(extra)
    out = subprocess.run([sys.executable, '-c', textwrap.dedent(script)], cwd=_ROOT, env=env,
                         capture_output=True, text=True, timeout=600, check=False)
    assert out.returncode == 0, out.stderr[-3000:]
    return out.stdout.split()


_COUNT = '''
    import sys
    import pycircuit.circuit.elements_hdl as E
    from pycircuit.circuit import hdl, spice_import
    lib = [c for c in vars(E).values()
           if isinstance(c, type) and issubclass(c, hdl.Behavioural) and 'analog' in c.__dict__]
    print(len(lib), sum(hdl.compile_pending(c) for c in lib),
          'pycircuit.circuit._hdl_cache' in sys.modules)
'''


def test_the_library_and_the_netlist_importer_compile_nothing_at_import():
    """No class of the library compiles at import -- the compile cache is
    not even opened -- though the netlist importer holds eighteen of them in
    its tables."""
    n, pending, cache = _child(_COUNT)
    assert int(n) >= 50
    assert pending == n
    assert cache == 'False'


def test_the_switch_compiles_at_creation():
    n, pending, cache = _child(_COUNT, PYCIRCUIT_HDL_DEFER='0')
    assert int(n) >= 50 and pending == '0' and cache == 'True'


_MODEL = '''
from pycircuit.circuit import hdl
from pycircuit.circuit.hdl import (Behavioural, Branch, Collapse, Contribution, Node, ddt,
                                   expl, limit_pnj, var)
from pycircuit.utilities.param import Parameter

if {defer}:
    hdl.defer_compile(__name__)


class Limited(Behavioural):
    """Chained, `$limit`, a collapse, a charge: most of what the compile sets."""
    instparams = [Parameter(name='IS', default=1e-14),
                  Parameter(name='rs', default=1.0),
                  Parameter(name='cj', default=1e-12)]

    @staticmethod
    def analog(plus, minus):
        i = Node('i')
        bj = Branch(i, minus)
        bs = Branch(plus, i)
        vt = var(hdl.vt(hdl.TEMP), 'vt')
        vj = var(limit_pnj(bj.V, IS, vt), 'vj')
        ij = var(IS * (expl(vj / vt) - 1.0), 'ij')
        return (Contribution(bj.I, ij + ddt(cj * vj)),
                Contribution(bs.I, bs.V / rs),
                Collapse(bs, rs <= 0))
'''

_END = '''
hdl.defer_compile(__name__, False)


class Later(Limited):
    """Created after the module's definitions: compiles at its creation."""
    instparams = list(Limited.instparams)
    analog = staticmethod(Limited.__dict__['analog'].__func__)
'''

_BROKEN = '''

class Broken(Behavioural):
    """Declared well, fails to compile: its analog raises when run."""
    instparams = [Parameter(name='g', default=1e-3)]

    @staticmethod
    def analog(plus, minus):
        raise ValueError('this model does not compile')
'''


@pytest.fixture
def cache_dir(tmp_path, monkeypatch):
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path / 'hdlcache'))
    monkeypatch.delenv('PYCIRCUIT_HDL_CACHE', raising=False)
    monkeypatch.setattr(hc, 'ENABLED', True)


def _load(tmp_path, modname, defer, extra=''):
    path = tmp_path / (modname + '.py')
    path.write_text(_MODEL.format(defer=defer) + extra + (_END if defer else ''))
    spec = importlib.util.spec_from_file_location(modname, str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop(modname, None)
    return mod


def _instance(cls):
    e = cls(*[Node(f'n{k}') for k in range(len(cls.terminals))])
    e.update_iparv()
    return e


def _bytes(e):
    rng = np.random.default_rng(20261102)
    out = []
    for _ in range(4):
        x = rng.uniform(-0.6, 0.6, e.n)
        with np.errstate(all='ignore'):
            out.append(tuple(np.asarray(v).tobytes() for v in (e.i(x), e.G(x), e.q(x), e.C(x))))
    return out


_DATA = ('terminals', 'linear', 'IC_KIND', '_hdl_paramnames', '_hdl_given_names', 'pcnr_junctions',
         '_hdl_backend_status')


def _as_eager(lazy, eager):
    assert not hdl.compile_pending(lazy)
    ## (both instantiated first: the backend resolves at a class's first
    ## instance)
    assert _bytes(_instance(lazy)) == _bytes(_instance(eager))
    for name in sorted(hdl._COMPILED_ATTRS):
        assert hasattr(lazy, name) == hasattr(eager, name), name
        assert not isinstance(lazy.__dict__.get(name), hdl._Pending), name
    for name in _DATA:
        assert getattr(lazy, name, None) == getattr(eager, name, None), name


@pytest.mark.parametrize('way', ['attribute', 'instance', 'subclass', 'absent', 'all'])
def test_every_way_in_compiles_and_reads_as_an_eager_compile(cache_dir, tmp_path, way):
    eager = _load(tmp_path, f'defer_eager_{way}', False).Limited
    mod = _load(tmp_path, f'defer_lazy_{way}', True)
    lazy = mod.Limited
    assert hdl.compile_pending(lazy) and not hdl.compile_pending(eager)
    ## (the module's definitions only: a class made in its name after them
    ## compiled at its creation -- reading `Limited`'s analog off its class
    ## body did not compile `Limited`)
    assert not hdl.compile_pending(mod.Later) and isinstance(mod.Later.__dict__.get('_hdl_info'), dict)
    assert isinstance(lazy.__dict__['terminals'], hdl._Pending)
    if way == 'attribute':
        assert lazy.terminals == eager.terminals
    elif way == 'instance':
        _instance(lazy)
    elif way == 'subclass':
        sub = type('Sub', (lazy,), {'__module__': lazy.__module__})
        assert hdl.compile_pending(lazy)
        _instance(sub)
    elif way == 'absent':
        ## a name this model's compile does not set reads as the bases give it
        assert hasattr(lazy, 'state_ic') == hasattr(eager, 'state_ic')
    else:
        assert hdl.compile_deferred(mod.__name__) == 1
    _as_eager(lazy, eager)


def test_a_compile_that_raises_holds_the_class_again(cache_dir, tmp_path):
    mod = _load(tmp_path, 'defer_broken', True, _BROKEN)
    for _ in range(2):
        with pytest.raises(ValueError, match='does not compile'):
            mod.Broken.terminals  # noqa: B018
        assert hdl.compile_pending(mod.Broken)
        assert isinstance(mod.Broken.__dict__['i'], hdl._Pending)
    with pytest.raises(ValueError, match='does not compile'):
        mod.Broken(Node('a'), Node('b'))
    ## (and the model's other class is untouched)
    assert hdl.compile_pending(mod.Limited)
    del mod
    gc.collect()


def test_no_compile_sets_a_name_the_descriptors_lack():
    """Every attribute the compile assigns on the class is held while the
    compile waits: statically, each `cls.<name> = ...` target of the
    metaclass's compile; dynamically, every name any compile of this process
    set (the library's, at the session's start, and every test's model)."""
    tree = ast.parse(textwrap.dedent(inspect.getsource(hdl.BehaviouralMeta._hdl_compile)))
    targets = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
            for t in (node.targets if isinstance(node, ast.Assign) else [node.target]):
                for e in (t.elts if isinstance(t, ast.Tuple) else [t]):
                    if (isinstance(e, ast.Attribute) and isinstance(e.value, ast.Name)
                            and e.value.id == 'cls'):
                        targets.add(e.attr)
    assert targets and targets <= hdl._COMPILED_ATTRS, sorted(targets - hdl._COMPILED_ATTRS)
    from pycircuit.circuit import elements_hdl
    hdl.compile_deferred(elements_hdl.__name__)
    assert hdl._UNLISTED == set()
