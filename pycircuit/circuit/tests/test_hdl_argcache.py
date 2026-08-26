# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""`hdl._args_of` is memoised on the instance (roadmap sec. 23).

It builds the trailing argument list every compiled function is called
with, so it runs once per ``i``/``G``/``q``/``C``/``u`` call -- once per
element per Newton iteration per timestep.  Profiling a 49-element
transient counted **125 029 calls**, each rebuilding an identical list:
794 k ``getattr`` and 421 k ``ParameterDict.__getattr__``.  Parameters do
not change during a transient, so the work was pure repetition; caching
it measured 1.16x end-to-end at 121 devices with bit-identical answers.

A value cache is only ever as good as its invalidation, and that is what
this file is mostly about.  `_args_of` reads TWO parameter dicts --
values from ``iparv``, givenness flags from ``ipar`` -- and the cache
hangs off the ``iparv`` observer, so the tests below walk every route a
parameter can move by:

* ``iparv`` written directly              -> ``notify`` -> ``update``
* ``ipar`` written directly               -> ``_ipar_changed`` ->
  ``update_iparv`` -> ``iparv.update_values`` -> ``notify`` -> ``update``
* ``ipar.set()`` and an explicit ``update_iparv()``

`T` is in the cache KEY rather than the cached list: ``epar.T`` varies
within a run without touching either dict, so it has to select between
entries instead of being stored beside them.

`test_the_cache_is_actually_used` is the one that gives the rest their
meaning -- every other test in this file passes with the memoisation
deleted.
"""

import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd, elements_hdl as eh
from pycircuit.circuit.circuit import SubCircuit, defaultepar
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import hdl

cm.default_toolkit = numeric


class _Epar(object):
    """A minimal epar carrying just a temperature."""

    def __init__(self, T):
        self.T = T


def _r_index(el):
    return el._hdl_paramnames.index('r')


def _r_of(el, T=300.0):
    return hdl._args_of(el, _Epar(T))[_r_index(el)]


def test_the_cache_is_actually_used():
    ## Without this, every other test in the file passes with the
    ## memoisation removed -- they would only be re-testing that
    ## `_args_of` reads live parameters, which it always did.
    el = eh.RHdl('a', gnd, r=1e3)
    first = hdl._args_of(el, _Epar(300.0))
    second = hdl._args_of(el, _Epar(300.0))
    assert first is second, 'the argument list is being rebuilt'
    assert el.__dict__.get('_hdl_args') is not None


def test_a_write_through_iparv_invalidates():
    el = eh.RHdl('a', gnd, r=1e3)
    assert _r_of(el) == 1e3
    el.iparv.r = 7e3
    assert _r_of(el) == 7e3


def test_a_write_through_ipar_invalidates():
    ## The indirect route, and the one a cache on `iparv` could plausibly
    ## miss: `ipar` is a DIFFERENT dict, and `_args_of` reads the
    ## givenness flags from it.
    el = eh.RHdl('a', gnd, r=1e3)
    assert _r_of(el) == 1e3
    el.ipar.r = 2e3
    assert _r_of(el) == 2e3


def test_ipar_set_and_update_iparv_invalidate():
    el = eh.RHdl('a', gnd, r=1e3)
    assert _r_of(el) == 1e3
    el.ipar.set(r=4e3)
    assert _r_of(el) == 4e3
    el.iparv.r = 5e3
    el.update_iparv()
    assert _r_of(el) == 4e3, 'update_iparv should restore from ipar'


def test_the_stamp_follows_the_parameter_not_the_cache():
    ## The failure this whole file guards against, at the level a user
    ## would meet it: a resistor whose conductance keeps its old value.
    el = eh.RHdl('a', gnd, r=1e3)
    x = np.zeros(2)
    g1 = float(np.asarray(el.G(x, defaultepar))[0, 0])
    el.ipar.r = 1e6
    g2 = float(np.asarray(el.G(x, defaultepar))[0, 0])
    assert g1 == pytest.approx(1e-3, rel=1e-12)
    assert g2 == pytest.approx(1e-6, rel=1e-12)


def test_temperature_selects_rather_than_collides():
    el = eh.RHdl('a', gnd, r=1e3)
    names = el._hdl_paramnames
    a300 = hdl._args_of(el, _Epar(300.0))
    a400 = hdl._args_of(el, _Epar(400.0))
    assert a300[len(names)] == 300.0
    assert a400[len(names)] == 400.0
    ## ...and the parameter values survive the temperature change
    assert a300[_r_index(el)] == a400[_r_index(el)] == 1e3
    ## going back to the first temperature must still give that T
    assert hdl._args_of(el, _Epar(300.0))[len(names)] == 300.0


def test_a_non_scalar_temperature_is_not_cached():
    ## An array or symbolic T would make the key comparison ambiguous.
    ## Such a call must still WORK, just not be memoised.
    el = eh.RHdl('a', gnd, r=1e3)
    T = np.array([300.0, 310.0])
    out = hdl._args_of(el, _Epar(T))
    assert np.array_equal(np.asarray(out[len(el._hdl_paramnames)]), T)
    assert el.__dict__.get('_hdl_args') is None


def test_a_transient_is_bit_identical_to_an_uncached_run():
    ## The acceptance measurement of roadmap sec. 23, shrunk to a test:
    ## the cache may not move a single bit of the answer.
    def build():
        c = SubCircuit()
        c['vs'] = eh.VSinHdl('n0', gnd, va=5.0, freq=1e3)
        for k in range(3):
            c['D%d' % k] = eh.DiodeHdl('n%d' % k, 'n%d' % (k + 1))
            c['R%d' % k] = eh.RHdl('n%d' % (k + 1), gnd, r=1e3)
            c['C%d' % k] = eh.CHdl('n%d' % (k + 1), gnd, c=1e-7)
        return c

    def solve():
        r = Transient(build()).solve(refnode=gnd, tend=2e-4, timestep=1e-5)
        return np.asarray(r.v('n3', gnd))

    cached = solve()

    ## Same circuit with the memoisation bypassed, in the same process.
    real = hdl._args_of

    def uncached(self, epar):
        n_par = len(self._hdl_paramnames)
        vals = hdl._params_of(self)
        return vals[:n_par] + [hdl._epar_T(epar)] + vals[n_par:]

    hdl._args_of = uncached
    try:
        plain = solve()
    finally:
        hdl._args_of = real

    assert np.array_equal(cached, plain), 'the cache changed the answer'


def test_a_sweep_over_a_parameter_sees_every_value():
    ## The shape that would break a naive cache in production: one
    ## element re-parameterised in a loop, which is what a sweep does.
    el = eh.RHdl('a', gnd, r=1e3)
    x = np.zeros(2)
    seen = []
    for r in (1e2, 1e3, 1e4, 1e5):
        el.ipar.r = r
        seen.append(float(np.asarray(el.G(x, defaultepar))[0, 0]))
    assert seen == [pytest.approx(1.0 / r, rel=1e-12)
                    for r in (1e2, 1e3, 1e4, 1e5)]


def test_an_int_temperature_does_not_answer_a_float_one():
    ## `defaultepar` carries T as an INT (300), and ``300 == 300.0``.
    ## Matching the key on value alone would return the int-T argument
    ## list to a float-T caller -- numerically harmless here, but a
    ## divergence from the uncached path, which passed the float
    ## through.  The type is therefore part of the key.
    el = eh.RHdl('a', gnd, r=1e3)
    iT = len(el._hdl_paramnames)
    a_int = hdl._args_of(el, _Epar(300))
    assert type(a_int[iT]) is int
    a_float = hdl._args_of(el, _Epar(300.0))
    assert type(a_float[iT]) is float, 'int entry answered a float call'
    ## ...and going back to the int must still give an int
    assert type(hdl._args_of(el, _Epar(300))[iT]) is int


def test_the_identity_fast_path_is_used():
    ## Discriminating on NaN, because the obvious test does not
    ## discriminate at all: asserting that two calls return the same list
    ## passes whether the hit came from `is` or from the value+type
    ## comparison, since both return the cached object.  (Written that way
    ## first, and it passed with the fast path deleted.)
    ##
    ## A NaN temperature separates them: `nan is nan` is True for one
    ## object, while `nan == nan` is False.  So the entry can only be
    ## reused through the identity check -- without it every call misses,
    ## rebuilds, and returns a NEW list.
    el = eh.RHdl('a', gnd, r=1e3)
    e = _Epar(float('nan'))
    first = hdl._args_of(el, e)
    assert hdl._args_of(el, e) is first, 'the identity fast path is gone'


def test_a_caller_that_rebuilds_T_still_gets_the_right_answer():
    ## Identity is STRICTER than value equality, so a caller handing a
    ## freshly built float each time simply falls through to the value
    ## test.  It must not get a wrong answer, and must not get a stale
    ## one after a parameter changes.
    el = eh.RHdl('a', gnd, r=1e3)
    iT = len(el._hdl_paramnames)
    for _ in range(3):
        a = hdl._args_of(el, _Epar(float('3e2')))   # a new object each time
        assert a[iT] == 300.0
        assert a[_r_index(el)] == 1e3
    el.ipar.r = 8e3
    assert hdl._args_of(el, _Epar(float('3e2')))[_r_index(el)] == 8e3


def test_identity_does_not_defeat_invalidation():
    ## The dangerous combination: the SAME epar object (so the identity
    ## key hits) across a parameter change.  `update()` must still have
    ## dropped the entry.
    el = eh.RHdl('a', gnd, r=1e3)
    e = _Epar(300.0)
    assert hdl._args_of(el, e)[_r_index(el)] == 1e3
    el.ipar.r = 5e3
    assert hdl._args_of(el, e)[_r_index(el)] == 5e3, 'identity hid a change'
