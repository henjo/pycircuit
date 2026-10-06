"""Parameter reads where they are many (speed round 12, stage 4): a read
through `ParameterDict.__getattr__` costs ~2.3 k instructions -- the failed
lookup and the hook, which no leaner body avoids (measured: a helper reading
the values dict 2.1 k, inlined 1.0 k).  So the reads are made fewer where
they are many: `BSource` reads its function once a call (twice before), and
an element's pack (`hdl._params_of`) reads each value from the values dict
exactly where `getattr` would reach it -- `getattr` itself everywhere else."""
import cProfile
import pstats

import numpy as np
import pytest

from pycircuit.circuit import circuit, hdl
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import Node
from pycircuit.circuit.elements import BSource, gnd
from pycircuit.utilities.param import ParameterDict


def _getattr_calls(fn):
    """`fn()`'s calls of `ParameterDict.__getattr__` (a profile: nothing
    patched, so nothing a genuineness check could see)."""
    pr = cProfile.Profile()
    pr.enable()
    out = fn()
    pr.disable()
    n = sum(v[1] for (path, _l, f), v in pstats.Stats(pr).stats.items()
            if f == '__getattr__' and path.endswith('param.py'))
    return out, n


def _mos():
    e = eh.MosLevel1Hdl(*[Node(f'n{k}') for k in range(4)], vto=0.5, lambd=0.02)
    e.update_iparv()
    return e


def _reference(e):
    vals = [getattr(e.iparv, n) for n in e._hdl_paramnames]
    return vals + [1.0 if e.ipar.is_given(n) else 0.0 for n in e._hdl_given_names]


def test_the_pack_reads_every_value_as_getattr_does():
    """The very objects `getattr` returns (the givenness flags after them:
    equal floats -- each function's `1.0` is its own constant)."""
    e = _mos()
    got, ref = hdl._params_of(e), _reference(e)
    k = len(e._hdl_paramnames)
    assert len(got) == len(ref)
    assert all(a is b for a, b in zip(got[:k], ref[:k], strict=True))
    assert got[k:] == ref[k:] and all(type(a) is float for a in got[k:])


def test_the_pack_reads_no_value_through_getattr():
    """The values dict serves every name (the parent: one `__getattr__` a
    name, ~40 for a MOSFET)."""
    e = _mos()
    _vals, n = _getattr_calls(lambda: hdl._params_of(e))
    assert n == 0, n


def test_switched_off_the_pack_reads_through_getattr(monkeypatch):
    """`PYCIRCUIT_PARAM_DIRECT=0`: `getattr` a name, as the parent."""
    e = _mos()
    monkeypatch.setattr(hdl, 'PARAM_DIRECT', False)
    vals, n = _getattr_calls(lambda: hdl._params_of(e))
    assert n == len(e._hdl_paramnames) and vals == _reference(e)


def test_a_shadow_a_subclass_and_a_patched_getattr_are_getattrs():
    """Where `getattr` would not reach the values dict -- an instance
    attribute shadowing a parameter, a `ParameterDict` subclass, a patched
    `__getattr__` -- the pack calls `getattr`, and gets what it gets."""
    e = _mos()
    name = e._hdl_paramnames[0]
    e.iparv.__dict__[name] = 'shadow'
    try:
        assert hdl._params_of(e)[0] == 'shadow'
        assert hdl._params_of(e) == _reference(e)
    finally:
        del e.iparv.__dict__[name]

    class Sub(ParameterDict):
        def __getattr__(self, key):
            if key in self.__dict__.get('_parameters', ()):
                return 'sub'
            return ParameterDict.__getattr__(self, key)
    e2 = _mos()
    sub = Sub(*[e2.iparv._parameters[n] for n in e2.iparv._paramnames])
    sub._values.update(e2.iparv._values)
    e2.iparv = sub
    assert hdl._params_of(e2)[:len(e2._hdl_paramnames)] == ['sub'] * len(e2._hdl_paramnames)

    e3 = _mos()
    real = ParameterDict.__getattr__
    try:
        ParameterDict.__getattr__ = lambda self, key: ('patched' if key in
                                                       self.__dict__.get('_parameters', ())
                                                       else real(self, key))
        assert hdl._params_of(e3)[0] == 'patched'
    finally:
        ParameterDict.__getattr__ = real
    assert hdl._params_of(e3) == _reference(e3)


def test_a_lookup_of_the_class_own_is_getattr(monkeypatch):
    """A `__getattribute__` on the class (none today: the default lookup)
    answers every read -- the pack calls `getattr` and gets its answer."""
    e = _mos()
    ga = object.__getattribute__

    def lookup(self, key):
        d = ga(self, '__dict__')
        if key != '_parameters' and key in d.get('_parameters', ()):
            return 'looked up'
        return ga(self, key)
    monkeypatch.setattr(ParameterDict, '__getattribute__', lookup, raising=False)
    k = len(e._hdl_paramnames)
    assert hdl._params_of(e)[:k] == ['looked up'] * k
    monkeypatch.undo()
    assert hdl._params_of(e) == _reference(e)


@pytest.mark.parametrize('m', ['i', 'G', 'q', 'C'])
def test_a_behavioural_source_reads_its_function_once_a_call(m):
    """`BSource` asked its parameters for the function twice a call (the
    parent: 2 reads; now 1) -- the same function, the same answer."""
    circuit.default_toolkit = circuit.numeric
    b = BSource(Node('p'), Node('n'), gnd, Node('c'),
                i_func=lambda u: 1e-3 * u ** 3, q_func=lambda u: 1e-9 * u ** 2)
    b.update_iparv()
    x = np.array([0.7, 0.1, 0.0, 0.0])
    _out, n = _getattr_calls(lambda: getattr(b, m)(x))
    assert n == 1, n
    ## (and without a function: the zeros, one read)
    b0 = BSource(Node('p'), Node('n'), gnd, Node('c'))
    b0.update_iparv()
    z, n0 = _getattr_calls(lambda: getattr(b0, m)(x))
    assert n0 == 1 and not np.any(z)
