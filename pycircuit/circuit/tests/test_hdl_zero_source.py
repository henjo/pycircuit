"""The zero-source skip (`_hdl_batch.zero_source`, speed round 6's second
commit, 2026-10-03): an hdl element whose generated `u`/`dudt` returns
literal zeros adds nothing to a bin that starts at +0.0, so the assembly
skips its call as it skips `Circuit.u`'s default -- the bytes the loop's;
an element with a source term, a shadow, an override and the 'ac' path
still called."""
import contextlib

import numpy as np
import pytest
import sympy

from pycircuit.circuit import _hdl_batch, hdl
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import R, SubCircuit, gnd
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, Parameter, var
from pycircuit.circuit.tests.test_hdl_batch import mos_chain
from pycircuit.circuit.tests.test_hdl_cbackend import CHAINED


@contextlib.contextmanager
def skipping(on):
    was = _hdl_batch.SKIP_ZERO_SOURCE
    _hdl_batch.SKIP_ZERO_SOURCE = on
    try:
        yield
    finally:
        _hdl_batch.SKIP_ZERO_SOURCE = was


@pytest.fixture
def source_calls(monkeypatch):
    """Every generated `u`/`dudt` body that reached its compiled function:
    `hdl._args_of` is what each calls first (a kernel never does)."""
    calls = []
    orig = hdl._args_of

    def spy(self, epar):
        calls.append(self)
        return orig(self, epar)
    monkeypatch.setattr(hdl, '_args_of', spy)
    return calls


class _Driven(Behavioural):
    """A conductance with a source term: `u` is not zeros."""
    instparams = [Parameter(name='gg', desc='g', unit='S', default=1e-3),  # noqa: RUF012 -- the DSL's form
                  Parameter(name='ia', desc='amplitude', unit='A', default=2e-3)]

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        u = var(b.V, 'u')
        return Contribution(b.I, gg * u + ia * sympy.sin(1e6 * hdl.TIME))  # noqa: F821


def _both(cir, m, *args):
    with skipping(True):
        a = getattr(cir, m)(*args)
    with skipping(False):
        b = getattr(cir, m)(*args)
    assert a.dtype == b.dtype and a.shape == b.shape
    assert a.tobytes() == b.tobytes(), m
    return a


def test_a_zero_source_is_skipped_and_the_vector_is_the_loops(source_calls):
    cir = mos_chain(6)
    for analysis in ('tran', 'dc', None):
        for m in ('u', 'dudt'):
            with skipping(True):
                a = getattr(cir, m)(1e-7, defaultepar, analysis)
            assert source_calls == [], (m, analysis)
            with skipping(False):
                b = getattr(cir, m)(1e-7, defaultepar, analysis)
            ## (outside the time domain the generated method answers zeros
            ## before its compiled function: no call on either side)
            assert len(source_calls) == (6 if analysis in cm.timedomain_analyses else 0)
            source_calls.clear()
            assert a.dtype == b.dtype == np.float64 and a.shape == b.shape
            assert a.tobytes() == b.tobytes(), (m, analysis)


def test_a_source_term_is_still_called(source_calls):
    cm.default_toolkit = cm.numeric
    cir = SubCircuit()
    for k in range(4):
        cir[f's{k}'] = _Driven(f'n{k}', gnd)
        cir[f'r{k}'] = R(f'n{k}', gnd, r=1e3)
        cir[f'M{k}'] = eh.MosLevel1Hdl(f'n{k}', gnd, gnd, gnd)
    assert not _hdl_batch.zero_source(type(cir['s0']), 'u')
    assert _hdl_batch.zero_source(type(cir['M0']), 'u')
    with skipping(True):
        a = cir.u(2.5e-7, defaultepar, 'tran')
    assert source_calls == [cir[f's{k}'] for k in range(4)]
    assert a.any()
    source_calls.clear()
    with skipping(False):
        b = cir.u(2.5e-7, defaultepar, 'tran')
    assert len(source_calls) == 8
    assert a.tobytes() == b.tobytes()


def test_the_ac_path_keeps_its_calls_and_its_dtype(source_calls):
    cir = mos_chain(6)
    with skipping(True):
        a = cir.u(0.0, defaultepar, 'ac')
    with skipping(False):
        b = cir.u(0.0, defaultepar, 'ac')
    assert a.dtype == b.dtype == np.complex128
    assert a.tobytes() == b.tobytes()


def test_a_shadow_and_an_override_are_called(source_calls):
    cir = mos_chain(6)
    el = cir['M2']
    el.__dict__['u'] = lambda t=0.0, epar=defaultepar, analysis=None, params_tree=None: (
        np.full(el.n, 1e-3))
    a = _both(cir, 'u', 1e-7, defaultepar, 'tran')
    assert a.any()
    del el.__dict__['u']

    class _Over(eh.MosLevel1Hdl):
        def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
            return np.full(self.n, 2e-3)
    assert not _hdl_batch.zero_source(_Over, 'u')
    assert _hdl_batch.zero_source(eh.MosLevel1Hdl, 'u')       # (the shared info)
    cir['over'] = _Over('d1', 'd2', gnd, gnd)
    a = _both(cir, 'u', 1e-7, defaultepar, 'tran')
    assert a.any()


@pytest.mark.parametrize('name', CHAINED)
def test_every_chained_library_class_has_a_zero_source(name):
    cls = getattr(eh, name)
    assert _hdl_batch.zero_source(cls, 'u') and _hdl_batch.zero_source(cls, 'dudt')
    assert '_u_zero' in cls._hdl_info


def test_the_source_reader():
    def fn(src):
        f = lambda: None
        f._src = src
        return f
    assert _hdl_batch._returns_zeros(fn('def _f(t, a):\n    return [0, 0.0, -0.0]'))
    assert not _hdl_batch._returns_zeros(fn('def _f(t, a):\n    return [0, 1e-3]'))
    assert not _hdl_batch._returns_zeros(fn('def _f(t, a):\n    x = 0\n    return [x]'))
    assert not _hdl_batch._returns_zeros(fn('def _f(t, a):\n    return 0'))
    assert not _hdl_batch._returns_zeros(lambda: None)
    assert not _hdl_batch.zero_source(R, 'u')


def test_a_transient_is_bit_identical():
    from pycircuit.circuit.transient import Transient
    got = {}
    for on in (True, False):
        cir = mos_chain(6)
        with skipping(on):
            res = Transient(cir, toolkit=cm.numeric).solve(
                tend=20 * 2e-8, timestep=2e-8, fixed_timestep=True)
        got[on] = np.asarray(res.x, dtype=float)
    assert got[True].tobytes() == got[False].tobytes()
