"""numpy's module reductions at their floor, exactly (speed round 12, stage
3): the numeric toolkit's `alltrue` -- the array's own `all` for an exact
ndarray -- the C lookup decided at the first differing entry, the memo's
keys without `asarray` for an exact float64 array.  Each served path is the
one it stands in for, result and type; the tests that count fail on the
parent tree."""
import numpy as np
import pytest

from pycircuit.circuit import _numeric, circuit
from pycircuit.circuit.elements import VS, C, R, SubCircuit, gnd
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.transient import Transient

ARRAYS = [np.array([True, False]), np.array([True, True]), np.array([], bool), np.array(True),
          np.array(False), np.array([1.0, np.nan]), np.array([0.0, -0.0]),
          np.array([[1, 2], [3, 0]]), np.zeros((0, 3)), np.array([1 + 0j, 0j]),
          np.array(['a', ''], dtype=object), np.array([None, 1], dtype=object),
          np.array([np.inf, -np.inf]), np.arange(10) > -1, np.ones((3, 4), bool)[:, ::2],
          np.array([1.0, 2.0])[::-1], np.asfortranarray(np.ones((2, 2)))]
OTHERS = [[True, False], (1, 1), True, False, np.True_, np.float64(0.0), 3,
          np.ma.masked_array([1, 0], mask=[False, True])]


def test_alltrue_is_numpys_all_on_every_kind_of_argument():
    for a in ARRAYS + OTHERS:
        with np.errstate(all='ignore'):
            r1, r2 = np.all(a), _numeric.alltrue(a)
        assert type(r1) is type(r2) and r1 == r2, (a, r1, r2)
    a = np.array([[True, False], [True, True]])
    for kw in ({'axis': 0}, {'axis': 1, 'keepdims': True}):
        assert np.array_equal(_numeric.alltrue(a, **kw), np.all(a, **kw))
    out = np.empty(2, bool)
    ## (numpy's positional order: a, axis, out)
    assert _numeric.alltrue(a, 0, out) is out


def test_an_exact_array_takes_its_own_method(monkeypatch):
    """The module function's dispatch is not entered for an exact ndarray
    (the parent: once a call); anything else is numpy's own call."""
    from numpy._core import fromnumeric
    calls = []
    real = fromnumeric._wrapreduction_any_all

    def counting(*a, **k):
        calls.append(1)
        return real(*a, **k)
    monkeypatch.setattr(fromnumeric, '_wrapreduction_any_all', counting)
    assert numeric.alltrue(np.array([True, True])) and calls == []
    assert numeric.alltrue([True, True]) and calls == [1]


def _transient():
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['v'] = VS('a', gnd, v=1.0)
    c['r'] = R('a', 'b', r=1e3)
    c['c'] = C('b', gnd, c=1e-9)
    return Transient(c, toolkit=circuit.numeric)


def test_a_lookup_that_misses_is_decided_at_the_first_entry(monkeypatch):
    """A first entry that differs is a miss without the full comparison (the
    parent: the comparison every time); NaN there never equal; a signed zero
    there equal, the full comparison deciding; a difference further on, the
    full comparison's miss; an equal copy, its hit."""
    tr = _transient()
    n = tr.cir.n
    x0 = np.linspace(0.1, 0.4, n)
    Cm = np.eye(n)
    tr._C_cache = (x0, Cm)
    calls = []
    real = numeric.alltrue

    def counting(a):
        calls.append(1)
        return real(a)
    monkeypatch.setattr(tr.toolkit, 'alltrue', counting)
    x = x0.copy()
    x[0] += 1.0
    assert tr._C_lookup(x) is None and calls == []
    x = x0.copy()
    x[0] = np.nan
    assert tr._C_lookup(x) is None and calls == []
    assert tr._C_lookup(x0.copy()) is Cm and calls == [1]
    z = x0.copy()
    z[0] = 0.0
    tr._C_cache = (z, Cm)
    zm = z.copy()
    zm[0] = -0.0
    assert tr._C_lookup(zm) is Cm and calls == [1, 1]
    x = z.copy()
    x[-1] += 1.0
    assert tr._C_lookup(x) is None and calls == [1, 1, 1]


@pytest.mark.parametrize('key', ['array', 'list', 'float32', 'strided', 'big_endian'])
def test_the_memo_keys_an_array_and_its_other_forms_alike(key):
    """The key is the state's float64 bytes however the state comes: an
    exact array is its own `asarray`, everything else converted as before."""
    tr = _transient()
    tr._memo_step()
    ## (values a float32 holds exactly: its conversion is the same state)
    x = (np.arange(tr.cir.n) + 1.0) * 0.25
    form = {'array': x.copy(), 'list': list(x), 'float32': x.astype(np.float32),
            'strided': np.repeat(x, 2)[::2], 'big_endian': x.astype('>f8')}[key]
    rec = {'C': np.eye(tr.cir.n)}
    tr._memo_put(form, rec)
    got = tr._memo_get(x)
    assert got is not None and got['C'] is rec['C']
    tr._memo_put(x, {'G': 1})
    assert tr._memo_get(form)['G'] == 1
