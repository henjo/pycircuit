"""The C backend as the default: 'auto' (the speed plan round 2, stage D;
Andreas, 2026-10-02: every chained model, built at a class's first
instance).

'auto' runs a chained class on C where C can be served -- the kernels
printed, cffi installed, the compile cache on, a compiler or the stored
objects -- and on numpy otherwise, quietly, with the reason in its status.
It resolves at the class's FIRST INSTANCE: nothing builds at import or at
class creation, every instance runs one backend from its first evaluation,
and what the analyses read about the backend at setup (the 'auto' Newton
options' `compiled_jacobian_size`) cannot change under them.
"""
import os
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _hdl_cbackend as cb
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import Node, defaultepar
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, var
from pycircuit.utilities.param import Parameter

_CC = cb.find_compiler()[0]
needs_cc = pytest.mark.skipif(_CC is None, reason='no C compiler')

_MODEL = '''
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, var
from pycircuit.utilities.param import Parameter


class Fresh(Behavioural):
    instparams = [Parameter(name='gg', desc='g', unit='S', default=%s)]

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        u = var(b.V, 'u')
        return Contribution(b.I, gg * u * u + u)               # noqa: F821
'''


@pytest.fixture
def default_backend():
    """Nothing pins a backend for the test -- neither the environment's
    variable nor the module flag.  Yields `touch(*classes)`: a class that
    already exists was attached under the run's own choice (a suite run
    with `PYCIRCUIT_HDL_BACKEND=numpy` attached the library to numpy at
    import), so a test re-attaches the classes it uses under the default --
    and afterwards they are attached under the run's choice again, or C
    would leak into the rest of a numpy run."""
    saved_env = os.environ.pop('PYCIRCUIT_HDL_BACKEND', None)
    saved_flag = hdl.BACKEND
    hdl.BACKEND = None
    touched = []

    def touch(*classes):
        for c in classes:
            hdl.set_backend(None, c)
            touched.append(c)
    try:
        yield touch
    finally:
        if saved_env is not None:
            os.environ['PYCIRCUIT_HDL_BACKEND'] = saved_env
        hdl.BACKEND = saved_flag
        for c in touched:
            hdl.set_backend(None, c)


def _fresh(tmp_path, tag, default='1.5'):
    """A chained class compiled from a module of its own (a new class, so
    unresolved); the same `default`, the same text and keys."""
    import importlib.util
    import sys
    d = tmp_path / tag
    d.mkdir()
    path = d / 'freshmod.py'
    path.write_text(_MODEL % default)
    spec = importlib.util.spec_from_file_location('freshmod', str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['freshmod'] = mod
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop('freshmod', None)
    return mod.Fresh


def _sos(d):
    return sorted(n for n in os.listdir(d) if n.endswith('.so'))


def test_the_default_backend_is_auto(default_backend):
    class M(Behavioural):
        instparams = [Parameter(name='gg', desc='g', unit='S', default=1.0)]

        @staticmethod
        def analog(p, m):
            b = Branch(p, m)
            return Contribution(b.I, gg * var(b.V, 'u'))       # noqa: F821
    assert hdl._backend_requested(M) == 'auto'


@needs_cc
def test_auto_builds_at_the_first_instance_not_at_class_creation(
        default_backend, tmp_path, monkeypatch):
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path / 'store'))
    ## (every `Fresh` prints the same C -- its parameters travel at run time
    ## -- so a kernel another test loaded would be served from memory)
    monkeypatch.setattr(cb, '_loaded', {})
    cls = _fresh(tmp_path, 't5', '1.25')
    assert cls._hdl_backend_status == 'auto (resolved at the first instance)'
    store = tmp_path / 'store'
    assert not store.exists() or _sos(store) == []
    assert not cls._hdl_info['funcs']['i'].__dict__.get('_hdl_c')
    e = cls(Node('a'), Node('b'))
    e.update_iparv()
    assert cls._hdl_backend_status == 'c'
    assert _sos(store)
    ## the C kernel answers the numpy function's bytes
    x = np.array([0.37, -0.11])
    args = list(hdl._args_of(e, defaultepar))
    for k in ('i', 'G'):
        f = cls._hdl_info['funcs'][k]
        assert getattr(e, k)(x).tobytes() == \
            np.asarray(f(x, *args), float).tobytes()


@needs_cc
def test_lifting_a_pin_restores_what_the_instances_ran(default_backend,
                                                      tmp_path, monkeypatch):
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path / 'store'))
    cls = _fresh(tmp_path, 't6', '1.5')
    cls(Node('a'), Node('b'))
    assert cls._hdl_backend_status == 'c'
    hdl.set_backend('numpy', cls)
    assert cls._hdl_backend_status == 'numpy'
    assert not cls._hdl_info.get('_c_bound')
    hdl.set_backend(None, cls)
    assert cls._hdl_backend_status == 'c'
    assert cls._hdl_info['_c_bound']


@needs_cc
@pytest.mark.parametrize('what', ['cffi', 'compiler', 'cache'])
def test_auto_runs_numpy_quietly_where_c_cannot_be_served(
        default_backend, tmp_path, monkeypatch, what):
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path / 'store'))
    if what == 'cffi':
        real = cb.importlib.util.find_spec
        monkeypatch.setattr(cb.importlib.util, 'find_spec',
                            lambda n, *a: None if n == 'cffi'
                            else real(n, *a))
        want = 'numpy (auto: cffi not installed)'
    elif what == 'compiler':
        monkeypatch.setattr(cb, '_compiler', (None, 'none here'))
        want = 'numpy (auto: no C compiler)'
    else:
        cls0 = _fresh(tmp_path, 't7', '1.75')   # compiled, the cache on
        monkeypatch.setenv('PYCIRCUIT_HDL_CACHE', '0')
        want = 'numpy (auto: the compile cache is off)'
    cls = cls0 if what == 'cache' else _fresh(tmp_path, 't8', '2.25')
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        e = cls(Node('a'), Node('b'))
        e.update_iparv()
        e.i(np.array([0.5, 0.0]))
    assert cls._hdl_backend_status == want
    assert not cls._hdl_info.get('_c_bound')


@needs_cc
def test_a_stored_object_serves_auto_without_a_compiler(
        default_backend, tmp_path, monkeypatch):
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path / 'store'))
    monkeypatch.setattr(cb, '_loaded', {})
    cls = _fresh(tmp_path, 't9', '2.5')
    cls(Node('a'), Node('b'))
    assert cls._hdl_backend_status == 'c'
    monkeypatch.setattr(cb, '_compiler', (None, 'none here'))
    monkeypatch.setattr(cb, '_loaded', {})
    again = _fresh(tmp_path, 't9b', '2.5')       # same text, same keys
    again(Node('a'), Node('b'))
    assert again._hdl_backend_status == 'c'


@needs_cc
def test_a_collapse_variant_resolves_at_its_own_first_instance(
        default_backend):
    default_backend(eh.DiodeSpiceHdl)
    e = eh.DiodeSpiceHdl(Node('a'), Node('k'), rs=10.0)
    var_ = type(e)
    assert var_._hdl_backend_status == 'c'
    assert var_._hdl_info['_c_bound']


@needs_cc
def test_every_chained_library_class_resolves_to_c(default_backend):
    """Every chained class's kernels build (or load) and serve its first
    instance: a C rendering the compiler refused would show here."""
    import inspect
    bad = []
    for name, c in sorted(vars(eh).items()):
        if not (inspect.isclass(c) and issubclass(c, Behavioural)
                and c.__module__ == eh.__name__
                and c._hdl_info.get('chained')):
            continue
        default_backend(c)
        e = c(*[Node(f'n{k}') for k in range(len(c.terminals))])
        if type(e)._hdl_backend_status != 'c':
            bad.append((name, type(e)._hdl_backend_status))
    assert not bad, bad


@needs_cc
def test_the_newton_options_read_the_backend_fixed_at_construction(
        default_backend):
    """'auto' resolves before any analysis: `compiled_jacobian_size` -- what
    the 'auto' Newton options read at a transient's setup -- is the same
    before the first evaluation and after a run."""
    from pycircuit.circuit import circuit
    from pycircuit.circuit._tran_newton import compiled_jacobian_size
    from pycircuit.circuit.elements import VS, R, SubCircuit, VSin, gnd
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    for n in ('g', 'd', 'vdd'):
        c.add_node(n)
    c['vdd'] = VS('vdd', gnd, v=1.2)
    c['vg'] = VSin('g', gnd, v=0.7, va=2e-2, freq=1e6)
    c['rl'] = R('vdd', 'd', r=5e3)
    default_backend(eh.MosLevel1Hdl)
    c['M'] = eh.MosLevel1Hdl('d', 'g', gnd, gnd)
    before = compiled_jacobian_size(c)
    assert type(c['M'])._hdl_backend_status == 'c'
    Transient(c, toolkit=circuit.numeric).solve(
        tend=2e-7, timestep=2e-8, fixed_timestep=True)
    assert compiled_jacobian_size(c) == before


@needs_cc
def test_a_jax_instance_leaves_the_class_to_its_jax_twin(default_backend,
                                                        tmp_path,
                                                        monkeypatch):
    jax = pytest.importorskip('jax')
    import jax.numpy as jnp

    from pycircuit.circuit.toolkit import jaxtoolkit
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path / 'store'))
    cls = _fresh(tmp_path, 't3', '2.75')
    ej = cls(Node('a'), Node('b'), toolkit=jaxtoolkit)
    ej.update_iparv()
    assert cls._hdl_backend_status == 'auto (resolved at the first instance)'
    got = jax.jit(lambda v: ej.i(v))(jnp.asarray([0.4, 0.1]))
    en = cls(Node('a'), Node('b'))
    en.update_iparv()
    assert cls._hdl_backend_status == 'c'
    assert np.allclose(np.asarray(got), en.i(np.array([0.4, 0.1])),
                       rtol=1e-12)
