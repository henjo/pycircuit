"""The C backend for chained `Behavioural` elements (`_hdl_cbackend`).

The claim under test is BIT-IDENTITY with the numpy path: an element
whose class runs `backend='c'` must return, from `i`, `G`, `q` and `C`,
the same BYTES the generated numpy source returns, over sweeps that
include the extremes the kernel's safety primitives exist for.  Three
exceptions, each measured, each named, each pinned by its own test
rather than hidden in a tolerance:

* **tanh** -- numpy ships its own vectorised `tanh`, which differs from
  libm's by an ulp in ~30% of arguments.  A chain using `tanh` agrees
  with the numpy path only to that ulp (amplified where the model
  cancels), and agrees BITWISE with the same numpy source run with
  libm's `tanh` -- which is what `test_tanh_is_the_whole_difference`
  asserts, converting the tolerance into a named cause.
* **the sign of exact zeros** -- the numpy path computes integer-typed
  subchains (`numpy.where(c, 1, 0)` is int64, and an integer zero
  carries no sign) where C is all doubles, so a C zero can be `-0.0`
  where numpy's is `+0.0`.  Values compare equal; only bytes differ.
* **a parameter set the numpy path cannot evaluate** -- an
  all-parameter denominator of exactly zero (`area/rc` at the default
  `rc=0`) raises `ZeroDivisionError` from Python where C returns
  `inf`.  That is a fact about the parameter set (`param = 0` is not
  an off-switch), so the sweeps here use parameter sets the numpy path
  accepts.

Then the machinery around the kernel: selection and status (a request
for C that cannot be served must SAY so, loudly, and run numpy), the
`.so` store (keyed by source + compiler + flags; corrupt objects
rebuilt; warm objects served with no compiler on the PATH), the compile
cache carrying the C text across interpreters, solver-level parity on
DC and transient, and the mutation checks that prove these tests can
fail: a broken kernel helper is caught by the sweep, and every key
ingredient actually changes the key.
"""

import contextlib
import itertools
import math
import os
import re
import shutil
import subprocess
import sys
import textwrap

import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit
import pycircuit.circuit.circuit as cm
from pycircuit.circuit import _hdl_cache as hc
from pycircuit.circuit import _hdl_cbackend as cb
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import Node, defaultepar
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Collapse, var, maxc, minc)
from pycircuit.utilities.param import Parameter

_CC = cb.find_compiler()[0]
needs_cc = pytest.mark.skipif(_CC is None, reason='no C compiler')

PDK = os.path.expanduser(
    '~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models')
needs_pdk = pytest.mark.skipif(not os.path.isdir(PDK),
                               reason='IHP Open PDK not present')


## ----------------------------------------------------------------------
## Helpers.

#: Parameter overrides that turn ON the parasitic branches these classes
#: collapse away at their defaults, so the sweep exercises them.
#:
#: NOT a defect list, and an earlier version of this comment said it was
#: ("classes whose DEFAULTS the numpy path cannot evaluate").  Measured
#: 2026-08-26 across all 37 library classes at 28 bias points and four
#: methods: no class raises at its defaults.  What raises is the
#: UNCOLLAPSED BASE function called directly -- which no instance runs,
#: because `Collapse` retargets every instance to a variant with the
#: dividing branch removed.  That is what `rc = 0` is declared to mean.
KW = {'GummelPoonNpnHdl': dict(rc=2.0, re=1.0, rb=100.0),
      'GummelPoonPnpHdl': dict(rc=2.0, re=1.0, rb=100.0),
      'GummelPoonNpnThermalHdl': dict(rc=2.0, re=1.0, rb=100.0,
                                      rth=200.0),
      'RThermalHdl': dict(rth=40.0),
      'PhotodiodeHdl': dict(rsh=1e6, rs=1.0),
      'LedHdl': dict(rs=2.0),
      'MesfetStatzHdl': dict(rs=1.0, rd=1.0),
      'DiodeSpiceThermalHdl': dict(rth=200.0, rs=1.0),
      'DiodeSpiceHdl': dict(IS=1e-14, n=1.6, rs=2.0, cjo=1e-12, tt=1e-9)}


def chained_classes():
    out = []
    for name, c in sorted(vars(eh).items()):
        if isinstance(c, type) and issubclass(c, Behavioural) and \
                c is not Behavioural and '_hdl_info' in c.__dict__:
            out.append(name)
    return out


ALL_CLASSES = chained_classes()
CHAINED = [n for n in ALL_CLASSES if getattr(eh, n)._hdl_info['chained']]
EAGER = [n for n in ALL_CLASSES if n not in CHAINED]


def _instance(cls, **kw):
    e = cls(*[Node('n%d' % k) for k in range(len(cls.terminals))], **kw)
    e.update_iparv()
    return e


@contextlib.contextmanager
def c_backend(cls):
    """The class pinned to C for the block, restored (to whatever the
    environment says) afterwards -- a leaked pin would silently turn
    the REST of the suite into a C-backend run."""
    hdl.set_backend('c', cls)
    try:
        yield
    finally:
        hdl.set_backend(None, cls)


@contextlib.contextmanager
def numpy_backend(cls):
    """The class pinned to numpy for the block -- the REFERENCE side of a
    comparison: under the default ('auto') a class with an instance runs
    C, so a reference taken through the element methods unpinned would
    compare C with C -- restored afterwards."""
    hdl.set_backend('numpy', cls)
    try:
        yield
    finally:
        hdl.set_backend(None, cls)


def _points(n, count=50, seed=0):
    rng = np.random.default_rng(seed)
    pts = [np.ascontiguousarray(p) for p in rng.uniform(-2, 2, (count, n))]
    pts += [np.full(n, v) for v in (1e30, -1e30, 100.0, -100.0, 0.7, 0.0)]
    return pts


def _c_funcs(cls):
    """The distinct chained functions of `cls` that carry a C kernel."""
    funcs = cls._hdl_info['funcs']
    seen, out = set(), []
    for name in cb.C_FUNCS:
        f = funcs.get(name)
        if f is None or id(f) in seen:
            continue
        seen.add(id(f))
        if f.__dict__.get('_hdl_c') is not None:
            out.append((name, f))
    return out


def _compare(ref, out):
    """`'equal'`, `'nan-bits'`, `'zero-sign'` or `'value'`.

    ⚠ **`'zero-sign'` used to mean two different things.**  It was
    returned whenever the values all compared equal but the bytes did
    not -- which is true of a signed zero AND of two NaNs with different
    sign bits, and those are not the same finding.  Split 2026-08-27
    (roadmap sec. 37), after a change was refused on the strength of "16
    signed zeros in PSP's G" that turned out to be **240 NaN sign bits**
    at 16 bias points, all of them points where `G` is already NaN in
    both paths (`fff8000000000000` against `7ff8000000000000`).

    The distinction matters because the two carry different weight.  A
    signed zero is a real, observable value: `1/-0.0` is `-inf`.  **A
    NaN's sign bit is not specified by IEEE-754** for the operations
    that produce one, is not preserved through arithmetic, and cannot be
    observed by any consumer -- every use of a NaN yields a NaN.
    Requiring the two backends to agree on it is requiring something the
    standard does not define.
    """
    ref = np.asarray(ref, float)
    out = np.asarray(out, float)
    if ref.tobytes() == out.tobytes():
        return 'equal'
    both_nan = np.isnan(ref) & np.isnan(out)
    eq = (ref == out) | both_nan
    if not bool(eq.all()):
        return 'value'
    ## Every differing byte is inside a NaN?  Then the disagreement is
    ## about a bit nothing defines.
    rb = ref.view(np.uint64).reshape(ref.shape)
    ob = out.view(np.uint64).reshape(out.shape)
    differing = rb != ob
    return 'nan-bits' if bool((differing & ~both_nan).sum() == 0) else 'zero-sign'


def _sweep(e, cls, pts):
    """{func name: {'equal': n, 'zero-sign': n, 'value': n}} comparing
    the C kernel against the numpy function it was printed from."""
    args = [float(v) for v in hdl._args_of(e, defaultepar)]
    out = {}
    for name, f in _c_funcs(cls):
        kern = f.__dict__['_hdl_c']
        tally = {'equal': 0, 'nan-bits': 0, 'zero-sign': 0,
                 'value': 0}
        for x in pts:
            with np.errstate(all='ignore'):
                ref = np.asarray(f(x, *args), float)
            tally[_compare(ref, kern(e, x, defaultepar))] += 1
        out[name] = tally
    return out


def _libm_tanh_twin(f):
    """`f._src` re-executed with `numpy.tanh` swapped for libm's --
    the ONLY difference from the real numpy path."""
    import types
    proxy = types.ModuleType('numpy')
    proxy.__dict__.update(np.__dict__)
    proxy.tanh = lambda v: np.float64(math.tanh(v))
    ns = hdl._chain_namespace(dict(hdl._KERNEL_NUMPY, _wrapfloor=np.floor))
    ns['numpy'] = proxy
    exec(compile(f._src, '<libm-tanh>', 'exec'), ns)
    return ns['_f']


## ----------------------------------------------------------------------
## Bit identity across the library.

@needs_cc
class TestLibraryBitIdentity(object):
    """Every chained library class, every per-iteration function, 56
    points including +-1e30: the C kernel returns the numpy bytes, with
    only the two named exceptions."""

    @pytest.mark.parametrize('name', CHAINED)
    def test_class_bitwise(self, name):
        cls = getattr(eh, name)
        e = _instance(cls, **KW.get(name, {}))
        n = len(hdl.x_layout(cls))
        args = [float(v) for v in hdl._args_of(e, defaultepar)]
        with c_backend(cls):
            assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
            found = _c_funcs(cls)
            assert found, 'no C kernels attached'
            uses_tanh = any('tanh' in f._src for _nm, f in found)
            if not uses_tanh:
                tallies = _sweep(e, cls, _points(n))
                for fname, t in tallies.items():
                    ## Value-identical everywhere; bytes identical bar
                    ## the zero-sign exception -- and a function that is
                    ## MOSTLY zero-sign would mean something else is
                    ## wrong, so require byte-equality to dominate.
                    ## `nan-bits` is counted separately and NOT bounded:
                    ## a NaN's sign bit is not defined by IEEE-754 (see
                    ## `_compare`).  A NaN that should not be there shows
                    ## up in the twin finiteness tests, not here.
                    assert t['value'] == 0, (fname, t)
                    assert t['equal'] >= t['zero-sign'], (fname, t)
                return
            ## tanh classes: the ulp of numpy's own tanh, amplified
            ## where the model cancels, is the ONLY allowed deviation
            ## from the true numpy path (the twin test pins that it IS
            ## tanh).  One shared band across all points and functions,
            ## absolute against the function's own scale.
            for fname, f in found:
                kern = f.__dict__['_hdl_c']
                for x in _points(n):
                    with np.errstate(all='ignore'):
                        ref = np.asarray(f(x, *args), float)
                    out = kern(e, x, defaultepar)
                    eq = (ref == out) | (np.isnan(ref) & np.isnan(out))
                    if bool(eq.all()):
                        continue
                    scale = np.nanmax(np.abs(np.where(np.isfinite(ref),
                                                      ref, 0.0)))
                    assert np.allclose(
                        np.where(eq, 0.0, ref), np.where(eq, 0.0, out),
                        rtol=1e-6, atol=1e-9 * max(1.0, scale)), \
                        (fname, x[:4])

    @pytest.mark.parametrize('name', sorted(
        n for n in CHAINED
        if any('tanh' in f._src for _nm, f in
               [(k, getattr(eh, n)._hdl_info['funcs'][k])
                for k in ('i', 'G')])))
    def test_tanh_is_the_whole_difference(self, name):
        """For every tanh-using class, the C kernel agrees BITWISE with
        the numpy source run with libm's tanh: the ulp against the real
        numpy path is numpy's own tanh, nothing else."""
        cls = getattr(eh, name)
        e = _instance(cls, **KW.get(name, {}))
        n = len(hdl.x_layout(cls))
        args = [float(v) for v in hdl._args_of(e, defaultepar)]
        with c_backend(cls):
            for fname, f in _c_funcs(cls):
                twin = _libm_tanh_twin(f)
                kern = f.__dict__['_hdl_c']
                for x in _points(n):
                    with np.errstate(all='ignore'):
                        ref = np.asarray(twin(x, *args), float)
                    got = _compare(ref, kern(e, x, defaultepar))
                    assert got != 'value', (fname, x[:4], got)

    @pytest.mark.parametrize('name', EAGER)
    def test_eager_classes_stay_numpy_and_say_so(self, name):
        cls = getattr(eh, name)
        with c_backend(cls):
            assert cls._hdl_backend_status.startswith('numpy (eager path')
        ## and the request did not break evaluation
        e = _instance(cls, **KW.get(name, {}))
        n = len(hdl.x_layout(cls))
        with np.errstate(all='ignore'):
            e.i(np.zeros(n))


## ----------------------------------------------------------------------
## Both arms.

class _OverflowArms(Behavioural):
    """The discarded arm overflows: `exp(u)` is inf past u ~ 710, and
    the selection must still return the finite arm -- the property
    `differentiable-numerics` protects, kept in C because `_sel` is a
    CALL whose arguments are both evaluated first."""
    instparams = [Parameter(name='gg', desc='scale', unit='A', default=1e-3)]

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        u = var(b.V, 'u')
        big = var(sympy.exp(u), 'big')
        kept = sympy.Piecewise((big, u < 20.0),
                               (sympy.exp(20.0) * (u - 19.0), True))
        return Contribution(b.I, gg * kept)                    # noqa: F821


@needs_cc
class TestBothArmsPreserved(object):

    def test_finite_where_the_discarded_arm_overflows(self):
        e = _instance(_OverflowArms)
        x = np.array([1000.0, 0.0])   # exp(1000) = inf in the dead arm
        with numpy_backend(_OverflowArms), np.errstate(all='ignore'):
            ref_i = e.i(x).copy()
            ref_G = e.G(x).copy()
        ## The VALUE survives the overflowing dead arm; the Jacobian at
        ## this bias is NaN on the numpy path (the dead arm's derivative
        ## poisons it -- exactly the `differentiable-numerics` case of
        ## an unclamped arm input), and the backend must NOT "fix" that:
        ## same bytes, NaN for NaN.
        assert np.isfinite(ref_i).all()
        assert not np.isfinite(ref_G).all()
        with c_backend(_OverflowArms):
            assert _OverflowArms._hdl_backend_status == 'c'
            with np.errstate(all='ignore'):
                ci, cG = e.i(x), e.G(x)
        assert ref_i.tobytes() == ci.tobytes()
        assert np.isfinite(ci).all()
        assert ref_G.tobytes() == cG.tobytes()

    def test_both_finite_below_the_seam(self):
        e = _instance(_OverflowArms)
        for v in (15.0, 25.0, 700.0):
            x = np.array([v, 0.0])
            with numpy_backend(_OverflowArms), np.errstate(all='ignore'):
                ref_i, ref_G = e.i(x).copy(), e.G(x).copy()
            with c_backend(_OverflowArms):
                with np.errstate(all='ignore'):
                    ci, cG = e.i(x), e.G(x)
            assert ref_i.tobytes() == ci.tobytes()
            assert ref_G.tobytes() == cG.tobytes()
            assert np.isfinite(ci).all()


## ----------------------------------------------------------------------
## The pow sentinel.

def _pow_sentinel():
    """A double where glibc's `pow(x, 2.0)` differs from `x*x` -- the
    value that catches a compiler folding `pow` (measured rate ~1 in
    1200, so 50 000 draws miss with probability ~1e-18)."""
    rng = np.random.default_rng(7)
    for x in rng.uniform(0.5, 2.0, 50000):
        if math.pow(x, 2.0) != x * x:
            return float(x)
    return None


class _SquareModel(Behavioural):
    instparams = [Parameter(name='gg', desc='scale', unit='A', default=1.0)]

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        u = var(b.V, 'u')            # an intermediate, so ** stays **
        return Contribution(b.I, gg * u ** 2)                  # noqa: F821


@needs_cc
class TestPowSentinel(object):
    """numpy's scalar `x ** 2` is glibc `pow`; gcc folds `pow(x, 2.0)`
    to the (more accurate) `x*x` at -O1+ unless told not to.  The
    sentinel pins the `-fno-builtin-pow` contract, and its mutation arm
    proves the sentinel can fail."""

    def _x(self):
        s = _pow_sentinel()
        if s is None:                            # pragma: no cover
            pytest.skip('no pow/mul-differing double found in 50k draws')
        return np.array([s, 0.0])

    def test_emitted_pow_is_a_real_pow(self):
        x = self._x()
        e = _instance(_SquareModel)
        assert 'pow(' in _SquareModel._hdl_info['funcs']['i']._csrc
        with numpy_backend(_SquareModel):
            ref = e.i(x).copy()
        with c_backend(_SquareModel):
            got = e.i(x)
        assert ref.tobytes() == got.tobytes()
        ## and the point is a real sentinel: pow and mul DO differ here
        s = float(x[0])
        assert math.pow(s, 2.0) != s * s

    def test_without_the_flag_the_sentinel_fires(self, tmp_path,
                                                 monkeypatch):
        """Drop `-fno-builtin-pow` and the same source, same key
        machinery, produces DIFFERENT bytes at the sentinel -- the
        test's power, measured."""
        x = self._x()
        e = _instance(_SquareModel)
        with numpy_backend(_SquareModel):
            ref = e.i(x).copy()
        ## ⚠ THE CLASS IS RE-BOUND AFTER THE FLAGS ARE RESTORED (2026-10-03):
        ## `c_backend`'s un-pin runs while the flags are still patched, so it
        ## re-bound `_SquareModel` to kernels compiled WITHOUT
        ## `-fno-builtin-pow`, and the module's later tests ran them (the
        ## leak detector's report).  The patches end in this test now, and
        ## the class is re-resolved under the real flags.
        try:
            with pytest.MonkeyPatch.context() as mp:
                mp.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))
                mp.setattr(cb, 'CFLAGS', tuple(
                    f for f in cb.CFLAGS if f != '-fno-builtin-pow'))
                mp.setattr(cb, '_loaded', {})
                with c_backend(_SquareModel):
                    assert _SquareModel._hdl_backend_status == 'c'
                    got = e.i(x)
        finally:
            hdl.set_backend(None, _SquareModel)
        assert ref.tobytes() != got.tobytes()
        assert np.allclose(ref, got, rtol=1e-14)   # one ulp, not garbage


## ----------------------------------------------------------------------
## Selection, status, fallback.

_SMALL = """
import numpy as np
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, var
from pycircuit.circuit.circuit import Node
from pycircuit.utilities.param import Parameter
import sympy

class M(Behavioural):
    instparams = [Parameter(name='gg', desc='g', unit='S', default=2.0)]
    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        u = var(b.V, 'u')
        return Contribution(b.I, gg * sympy.tanh(u) + gg * u ** 2)

e = M(Node('a'), Node('b'))
e.update_iparv()
x = np.array([0.625, 0.0])
print('STATUS', M._hdl_backend_status)
print('BYTES', e.i(x).tobytes().hex(), e.G(x).tobytes().hex())
"""


def _run(code, env=None, check=True, script=None):
    """`code` in a fresh interpreter.  With `script` (a path), the code
    runs from a real file -- `inspect.getsource` works, so the class is
    compile-cacheable; `-c` classes never are."""
    full = dict(os.environ)
    full.update(env or {})
    if script is not None:
        with open(script, 'w') as fh:
            fh.write(code)
        cmd = [sys.executable, str(script)]
    else:
        cmd = [sys.executable, '-c', code]
    p = subprocess.run(cmd, capture_output=True, text=True, env=full)
    if check:
        assert p.returncode == 0, p.stderr
    return p


@needs_cc
class TestSelectionAndFallback(object):

    def test_env_var_selects_c(self, tmp_path):
        p = _run(_SMALL, env={'PYCIRCUIT_HDL_BACKEND': 'c',
                              'PYCIRCUIT_HDL_CACHE_DIR': str(tmp_path)})
        assert 'STATUS c' in p.stdout

    def test_env_var_typo_is_loud(self, tmp_path):
        p = _run(_SMALL, env={'PYCIRCUIT_HDL_BACKEND': 'C99',
                              'PYCIRCUIT_HDL_CACHE_DIR': str(tmp_path)},
                 check=False)
        assert p.returncode != 0
        assert 'unknown HDL backend' in p.stderr

    def test_no_compiler_falls_back_and_says_why(self, tmp_path):
        """PATH stripped: `backend='c'` yields the numpy result, the
        status names the missing compiler, a warning is issued and
        NOTHING raises."""
        ref = _run(_SMALL, env={'PYCIRCUIT_HDL_BACKEND': 'numpy',
                                'PYCIRCUIT_HDL_CACHE_DIR': str(tmp_path)})
        code = ('import warnings\n'
                'with warnings.catch_warnings(record=True) as w:\n'
                '    warnings.simplefilter("always")\n' +
                textwrap.indent(_SMALL, '    ') +
                '\nprint("WARNED", any("no C compiler" in str(x.message)'
                ' for x in w))\n')
        p = _run(code, env={'PYCIRCUIT_HDL_BACKEND': 'c',
                            'PYCIRCUIT_HDL_CACHE_DIR': str(tmp_path),
                            'PATH': '/nonexistent'})
        assert 'STATUS numpy (compile failed: no C compiler' in p.stdout
        assert 'WARNED True' in p.stdout
        assert [ln for ln in p.stdout.splitlines() if
                ln.startswith('BYTES')] == \
            [ln for ln in ref.stdout.splitlines() if ln.startswith('BYTES')]

    def test_per_class_attribute_and_explain(self, tmp_path, monkeypatch):
        monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))

        class M(Behavioural):
            hdl_backend = 'c'
            instparams = [Parameter(name='gg', desc='g', unit='S',
                                    default=1.0)]

            @staticmethod
            def analog(p, m):
                b = Branch(p, m)
                u = var(b.V, 'u')
                return Contribution(b.I, gg * u ** 3)          # noqa: F821

        assert M._hdl_backend_status == 'c'
        assert 'backend: c' in hdl.explain(M, source=False, symbolic=False)
        hdl.set_backend('numpy', M)
        assert M._hdl_backend_status == 'numpy'
        assert 'backend: numpy' in hdl.explain(M, source=False,
                                               symbolic=False)

    def test_set_backend_default_applies_to_new_classes(self, tmp_path,
                                                        monkeypatch):
        monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))
        hdl.set_backend('c')
        try:
            class M(Behavioural):
                instparams = [Parameter(name='gg', desc='g', unit='S',
                                        default=1.0)]

                @staticmethod
                def analog(p, m):
                    b = Branch(p, m)
                    u = var(b.V, 'u')
                    return Contribution(b.I, gg * u * u)       # noqa: F821

            assert M._hdl_backend_status == 'c'
        finally:
            hdl.set_backend(None)

    def test_collapse_variant_follows_the_pin(self):
        """A collapsing model RUNS as a compiled variant subclass;
        pinning the base must reach it, or `set_backend('c', cls)`
        silently leaves every existing instance on numpy (found by the
        benchmark: MosLevel3 at exactly 1x)."""
        cls = eh.MosLevel3Hdl
        e = _instance(cls)
        var_cls = type(e)
        assert var_cls is not cls, 'expected a collapse variant'
        n = len(hdl.x_layout(var_cls))
        x = np.linspace(-0.3, 0.3, n)
        with numpy_backend(cls), np.errstate(all='ignore'):
            assert var_cls._hdl_backend_status == 'numpy'
            ref = e.G(x).copy()
        with c_backend(cls):
            assert var_cls._hdl_backend_status == 'c', \
                var_cls._hdl_backend_status
            with np.errstate(all='ignore'):
                got = e.G(x)
        assert ref.tobytes() == got.tobytes()
        ## unpinned: the variant -- which has an instance -- runs whatever
        ## the environment's default gives one ('auto': C)
        want = {'numpy': 'numpy', 'c': 'c', 'auto': 'c'}[
            hdl._backend_requested(var_cls)]
        assert var_cls._hdl_backend_status == want

    def test_results_identical_through_the_element_methods(self):
        """The same instance, backend toggled around it: `i/G/q/C`
        bytes unchanged, and the packed-parameter cache follows a
        parameter update."""
        cls = eh.GummelPoonNpnHdl
        e = _instance(cls, **KW['GummelPoonNpnHdl'])
        n = len(hdl.x_layout(cls))
        x = np.linspace(-0.4, 0.4, n)
        with numpy_backend(cls), np.errstate(all='ignore'):
            ref = [getattr(e, m)(x).copy() for m in ('i', 'G', 'q', 'C')]
        with c_backend(cls):
            with np.errstate(all='ignore'):
                got = [getattr(e, m)(x) for m in ('i', 'G', 'q', 'C')]
            assert all(a.tobytes() == b.tobytes()
                       for a, b in zip(ref, got))
            ## a parameter change must invalidate the packed vector
            e.ipar.rc = 4.0
            e.update_iparv()
            with np.errstate(all='ignore'):
                after_c = e.i(x).copy()
        with numpy_backend(cls), np.errstate(all='ignore'):
            after_np = e.i(x)
        assert after_c.tobytes() == after_np.tobytes()
        assert after_c.tobytes() != ref[0].tobytes()


## ----------------------------------------------------------------------
## The .so store.

@needs_cc
class TestSoStore(object):

    def _model_source(self, tag):
        return _SMALL.replace("default=2.0", "default=2.%d" % tag)

    def test_round_trip_second_interpreter_needs_no_compiler(self,
                                                             tmp_path):
        """Interpreter one compiles; interpreter two, with the PATH
        stripped, must still run backend 'c' from the stored `.so`
        (dlopen only) and return the same bytes."""
        env = {'PYCIRCUIT_HDL_BACKEND': 'c',
               'PYCIRCUIT_HDL_CACHE_DIR': str(tmp_path)}
        p1 = _run(_SMALL, env=env)
        assert 'STATUS c' in p1.stdout
        sos = [f for f in os.listdir(tmp_path) if f.endswith('.so')]
        assert sos, 'no shared object was stored'
        p2 = _run(_SMALL, env=dict(env, PATH='/nonexistent'))
        assert 'STATUS c' in p2.stdout
        b1 = [ln for ln in p1.stdout.splitlines() if ln.startswith('BYTES')]
        b2 = [ln for ln in p2.stdout.splitlines() if ln.startswith('BYTES')]
        assert b1 == b2 and b1

    def test_corrupt_so_is_rebuilt(self, tmp_path):
        env = {'PYCIRCUIT_HDL_BACKEND': 'c',
               'PYCIRCUIT_HDL_CACHE_DIR': str(tmp_path)}
        p1 = _run(_SMALL, env=env)
        for f in os.listdir(tmp_path):
            if f.endswith('.so'):
                with open(os.path.join(tmp_path, f), 'wb') as fh:
                    fh.write(b'not an ELF object')
        p2 = _run(_SMALL, env=env)
        assert 'STATUS c' in p2.stdout
        assert [ln for ln in p1.stdout.splitlines()
                if ln.startswith('BYTES')] == \
               [ln for ln in p2.stdout.splitlines()
                if ln.startswith('BYTES')]

    def test_every_key_ingredient_changes_the_key(self, monkeypatch):
        """Source, kernel prelude, flags, compiler identity: each must
        move the key, or a stale binary would be served."""
        csrc = 'void hdl_fn(const double *x, const double *p, '\
               'double *out) { out[0] = x[0]; }\n'
        k0 = cb.source_key(csrc)
        assert cb.source_key(csrc.replace('x[0]', 'p[0]')) != k0
        monkeypatch.setattr(cb, 'CFLAGS', cb.CFLAGS + ('-DX',))
        k_flags = cb.source_key(csrc)
        monkeypatch.undo()
        assert k_flags != k0
        monkeypatch.setattr(hdl, '_KERNEL_C',
                            hdl._KERNEL_C + '/* mutated */\n')
        k_kernel = cb.source_key(csrc)
        monkeypatch.undo()
        assert k_kernel != k0
        ## the compiler's identity lives in the FILENAME, so that a
        ## warm store can be served with no compiler at all
        p0 = os.path.basename(cb.so_path(k0))
        monkeypatch.setattr(cb, '_compiler', ('/usr/bin/cc', 'other 1.0'))
        p_cc = os.path.basename(cb.so_path(cb.source_key(csrc)))
        monkeypatch.undo()
        assert p_cc != p0
        assert cb.source_key(csrc) == k0

    def test_broken_kernel_helper_is_caught_by_the_sweep(self, tmp_path,
                                                         monkeypatch):
        """Mutation check: `_npmax` with numpy's NaN rule REMOVED must
        make the bit-identity sweep fail -- a sweep that would still
        pass would be no evidence."""
        monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))

        class M(Behavioural):
            instparams = [Parameter(name='gg', desc='g', unit='S',
                                    default=1.0)]

            @staticmethod
            def analog(p, m):
                b = Branch(p, m)
                u = var(b.V, 'u')
                return Contribution(b.I, gg * maxc(u, 0.5))    # noqa: F821

        e = _instance(M)
        line = next(ln for ln in hdl._KERNEL_C.splitlines()
                    if ln.startswith('static inline double _npmax('))
        broken = hdl._KERNEL_C.replace(line, line.replace(' || a != a', ''))
        assert broken != hdl._KERNEL_C
        monkeypatch.setattr(hdl, '_KERNEL_C', broken)
        with c_backend(M):
            assert M._hdl_backend_status == 'c'
            tallies = _sweep(e, M, [np.array([np.nan, 0.0])])
        assert any(t['value'] for t in tallies.values()), tallies

    def test_compile_error_reports_and_runs_numpy(self, tmp_path,
                                                  monkeypatch):
        monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))

        class M(Behavioural):
            instparams = [Parameter(name='gg', desc='g', unit='S',
                                    default=1.0)]

            @staticmethod
            def analog(p, m):
                b = Branch(p, m)
                u = var(b.V, 'u')
                return Contribution(b.I, gg * u)               # noqa: F821

        e = _instance(M)
        with numpy_backend(M):
            ref = e.i(np.array([0.5, 0.0])).copy()
        monkeypatch.setattr(hdl, '_KERNEL_C',
                            hdl._KERNEL_C + 'this is not C\n')
        ## (a `CostWarning` since 2026-10-01; the default `UserWarning`
        ## before)
        from pycircuit.circuit.simwarnings import CostWarning
        with pytest.warns(CostWarning, match='build failed'):
            hdl.set_backend('c', M)
        try:
            assert M._hdl_backend_status.startswith(
                'numpy (compile failed:')
            got = e.i(np.array([0.5, 0.0]))
            assert got.tobytes() == ref.tobytes()
        finally:
            ## (numpy, not the default: 'auto' would build again against
            ## the broken prelude, still patched in)
            hdl.set_backend('numpy', M)


## ----------------------------------------------------------------------
## The compile cache carries the C text.

@needs_cc
class TestCompileCacheCarriesC(object):

    def test_freeze_thaw_preserves_the_c_rendering(self):
        info = eh.GummelPoonNpnHdl._hdl_info
        thawed = hc.thaw(hc.freeze(info))
        for k in ('i', 'G'):
            a, b = info['funcs'][k], thawed['funcs'][k]
            assert getattr(a, '_csrc', None) is not None
            assert a._csrc == b._csrc
            assert tuple(a._cshape) == tuple(b._cshape)
            assert tuple(a._clayout) == tuple(b._clayout)

    def test_cache_hit_still_serves_backend_c(self, tmp_path):
        """Second interpreter: pickle hit AND `.so` hit -- no sympy, no
        compiler run -- same bytes."""
        env = {'PYCIRCUIT_HDL_BACKEND': 'c',
               'PYCIRCUIT_HDL_CACHE': '1',
               'PYCIRCUIT_HDL_CACHE_DIR': str(tmp_path)}
        code = _SMALL + "\nfrom pycircuit.circuit import _hdl_cache\n" \
                        "print('CACHE', M._hdl_cache_status)\n"
        script = str(tmp_path / 'small_model_script.py')
        p1 = _run(code, env=env, script=script)
        assert 'CACHE miss' in p1.stdout, p1.stdout
        p2 = _run(code, env=dict(env, PATH='/nonexistent'), script=script)
        assert 'CACHE hit' in p2.stdout
        assert 'STATUS c' in p2.stdout
        assert [ln for ln in p1.stdout.splitlines()
                if ln.startswith('BYTES')] == \
               [ln for ln in p2.stdout.splitlines()
                if ln.startswith('BYTES')]


## ----------------------------------------------------------------------
## Solver parity.

@needs_cc
class TestSolverParity(object):
    """The Newton path sees only the returned arrays, and those are
    byte-identical -- so DC and transient answers must match the numpy
    backend to solver precision."""

    def _bjt_circuit(self):
        from pycircuit.circuit.elements import SubCircuit, VS, R
        from pycircuit.circuit import gnd
        pycircuit.circuit.circuit.default_toolkit = \
            __import__('pycircuit.circuit.toolkit',
                       fromlist=['numeric']).numeric
        c = SubCircuit()
        vc, vb, out = c.add_node('vc'), c.add_node('vb'), c.add_node('out')
        c['VC'] = VS(vc, gnd, v=5.0)
        c['VB'] = VS(vb, gnd, v=0.7)
        c['RB'] = R(vb, 'b', r=1e4)
        c['RC'] = R(vc, out, r=2e3)
        c['Q1'] = eh.GummelPoonNpnHdl(out, 'b', gnd,
                                      **KW['GummelPoonNpnHdl'])
        c.update_iparv()
        return c, out

    def test_bjt_dc(self):
        from pycircuit.circuit.dcanalysis import DC
        from pycircuit.circuit.toolkit import numeric
        from pycircuit.circuit import gnd
        c, out = self._bjt_circuit()
        with numpy_backend(eh.GummelPoonNpnHdl):
            ref = float(DC(c, toolkit=numeric).solve().v(out, gnd))
        with c_backend(eh.GummelPoonNpnHdl):
            got = float(DC(c, toolkit=numeric).solve().v(out, gnd))
        assert abs(got - ref) <= 1e-12 * max(1.0, abs(ref))

    def test_bjt_transient(self):
        from pycircuit.circuit.transient import Transient
        from pycircuit.circuit.elements import SubCircuit, VSin, R, C as Cap
        from pycircuit.circuit.toolkit import numeric
        from pycircuit.circuit import gnd
        pycircuit.circuit.circuit.default_toolkit = numeric

        def build():
            c = SubCircuit()
            vc, vb, out = (c.add_node('vc'), c.add_node('vb'),
                           c.add_node('out'))
            c['VC'] = __import__('pycircuit.circuit.elements',
                                 fromlist=['VS']).VS(vc, gnd, v=5.0)
            c['VB'] = VSin(vb, gnd, vo=0.65, va=0.05, freq=1e6)
            c['RB'] = R(vb, 'b', r=1e4)
            c['RC'] = R(vc, out, r=2e3)
            c['CL'] = Cap(out, gnd, c=1e-11)
            c['Q1'] = eh.GummelPoonNpnHdl(out, 'b', gnd,
                                          **KW['GummelPoonNpnHdl'])
            c.update_iparv()
            return c

        def wave():
            ## ONE Newton on both backends: the 'auto' `chord_jacobian`
            ## decides by the Jacobian's cost, and the C-bound model counts a
            ## hundredth of the numpy one (`C_KERNEL_SHARE`: 57 KB against
            ## 0.6 KB), so left to it the two runs take different Newtons
            ## and agree to the Newton tolerance (9.4e-10), not to rounding
            res = Transient(build(), toolkit=numeric,
                            chord_jacobian=False).solve(
                tend=2e-6, timestep=2e-8, fixed_timestep=True)
            return np.asarray(res.v('out', gnd), float)

        with numpy_backend(eh.GummelPoonNpnHdl):
            ref = wave()
        with c_backend(eh.GummelPoonNpnHdl):
            got = wave()
        assert got.shape == ref.shape
        assert np.allclose(got, ref, rtol=1e-12, atol=1e-12)

    @needs_pdk
    def test_psp_dc_and_transient(self):
        from pycircuit.circuit.dcanalysis import DC
        from pycircuit.circuit.transient import Transient
        from pycircuit.circuit.elements import SubCircuit, VS, VSin, R
        from pycircuit.circuit.toolkit import numeric
        from pycircuit.circuit import gnd, psp_scaling
        from pycircuit.utilities import spicecard
        from pycircuit.circuit.compact import PspMosLongChannel
        pycircuit.circuit.circuit.default_toolkit = numeric
        deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                              section='mos_tt')
        w, l = 10e-6, 1e-6
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1, m=1,
                              pre_layout=1), w=w, l=l, T=273.15 + 27)

        def build(ac):
            c = SubCircuit()
            vdd, g, d = c.add_node('vdd'), c.add_node('g'), c.add_node('d')
            c['VDD'] = VS(vdd, gnd, v=1.2)
            c['VG'] = VSin(g, gnd, vo=0.8, va=0.1, freq=1e6) if ac \
                else VS(g, gnd, v=0.8)
            c['RL'] = R(vdd, d, r=10e3)
            c['M1'] = PspMosLongChannel(d, g, gnd, gnd, **kw)
            c.update_iparv()
            return c

        def dc():
            return float(DC(build(False), toolkit=numeric)
                         .solve().v('d', gnd))

        def tran():
            res = Transient(build(True), toolkit=numeric).solve(
                tend=1e-6, timestep=2e-8, fixed_timestep=True)
            return np.asarray(res.v('d', gnd), float)

        with numpy_backend(PspMosLongChannel):
            ref_dc, ref_tr = dc(), tran()
        with c_backend(PspMosLongChannel):
            assert PspMosLongChannel._hdl_backend_status == 'c'
            got_dc, got_tr = dc(), tran()
        assert abs(got_dc - ref_dc) <= 1e-12 * max(1.0, abs(ref_dc))
        assert np.allclose(got_tr, ref_tr, rtol=1e-12, atol=1e-12)


## ----------------------------------------------------------------------
## PSP bit identity (the spike's sweep, as a test).

@needs_cc
@needs_pdk
@pytest.mark.slow
class TestPspBitIdentity(object):
    """>= 1000 points of the spike's sweep (extreme grid + reference
    biases): `i`, `G`, `q` byte-identical; `C` value-identical with
    only zero-sign byte differences (the integer-lattice exception --
    PSP's charge Jacobian has integer-valued zero cells)."""

    @pytest.fixture(scope='class')
    def psp(self):
        from pycircuit.circuit import psp_scaling
        from pycircuit.utilities import spicecard
        from pycircuit.circuit.toolkit import numeric
        pycircuit.circuit.circuit.default_toolkit = numeric
        deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                              section='mos_tt')
        w, l = 10e-6, 1e-6
        kw = psp_scaling.to_long_channel(
            deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1, m=1,
                              pre_layout=1), w=w, l=l, T=273.15 + 27)
        from pycircuit.circuit.compact import PspMosLongChannel
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        return e

    def _points(self, e):
        ext = (-1e30, -100.0, -1.0, 0.0, 0.7, 1.2, 100.0, 1e30)
        combos = list(itertools.product(ext, repeat=4))
        rng = np.random.default_rng(0)
        keep = rng.choice(len(combos), 700, replace=False)
        raw = [combos[k] for k in sorted(keep)]
        ## reference-sweep biases, without importing the benchmark
        for vg in np.linspace(0.0, 1.2, 61):
            raw.append((0.05, vg, 0.0, 0.0))
            raw.append((1.2, vg, 0.0, 0.0))
            raw.append((0.05, vg, 0.0, -0.6))
        for vd in np.linspace(0.0, 1.2, 61):
            raw.append((vd, 0.6, 0.0, 0.0))
            raw.append((vd, 1.2, 0.0, 0.0))
        with np.errstate(all='ignore'):
            return [np.ascontiguousarray(e.bias(*c), float) for c in raw]

    def test_the_spikes_sweep(self, psp):
        e = psp
        cls = type(e)
        pts = self._points(e)
        assert len(pts) >= 1000
        with c_backend(cls):
            assert cls._hdl_backend_status == 'c'
            tallies = _sweep(e, cls, pts)
        assert set(tallies) == {'i', 'G', 'q', 'C'}
        for k in ('i', 'G', 'q'):
            t = tallies[k]
            ## Strict: no value may differ and no ZERO may change sign.
            ##
            ## `nan-bits` is deliberately not required to be zero, and
            ## that is a 2026-08-27 correction rather than a relaxation.
            ## It was folded into `zero-sign` before, and on that reading
            ## a change was refused for "16 signed zeros in G" that were
            ## really 240 NaN SIGN BITS at 16 bias points -- every one of
            ## them a point where `G` is already NaN in both paths
            ## (biases like +-1e30, 20 of 36 cells NaN).  IEEE-754 does
            ## not define that bit, arithmetic does not preserve it, and
            ## no consumer can see it.  Requiring it made the suite
            ## sensitive to something that is not a computation.
            ##
            ## What still fails here: any differing VALUE, and any zero
            ## whose sign flips -- `1/-0.0` is `-inf`, so that one is
            ## observable and stays banned.
            assert t['value'] == 0 and t['zero-sign'] == 0, (k, t)
        assert tallies['C']['value'] == 0, tallies['C']
        ## The NaN-bit count is asserted as a RECORDED number rather than
        ## left free, so a change that starts producing NaNs somewhere
        ## new is still caught here.  Measured 2026-08-27: 16 points on
        ## `G`, none on `i` or `q`.
        assert tallies['G']['nan-bits'] <= 32, tallies['G']
        assert tallies['i']['nan-bits'] == 0, tallies['i']
        assert tallies['q']['nan-bits'] == 0, tallies['q']


## ----------------------------------------------------------------------
## The kept libm calls, declared const (F1b, 2026-10-02).

_KEPT = tuple(f[len('-fno-builtin-'):] for f in cb.CFLAGS
              if f.startswith('-fno-builtin-'))


def _const_declared(prelude):
    return re.findall(r'double (\w+)\([^)]*\) __attribute__\(\(const\)\);',
                      prelude)


def test_every_kept_libm_call_is_declared_const():
    """The prelude declares `const` exactly the functions `CFLAGS` keeps
    real libm calls: one the flags do not keep, gcc may already fold its
    own way; one they keep and the prelude misses is computed again at
    every repeat."""
    declared = _const_declared(hdl._KERNEL_C)
    assert len(declared) == len(set(declared))
    assert _KEPT and set(declared) == set(_KEPT)


class _RepeatedPow(Behavioural):
    instparams = [Parameter(name='gg', desc='scale', unit='A', default=1.0)]

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        u = var(b.V, 'u')
        w = var(u ** 2.5, 'w')
        return Contribution(b.I, gg * w * w + w)                # noqa: F821


@needs_cc
@pytest.mark.skipif(not shutil.which('objdump'),
                    reason='no objdump')
def test_const_merges_the_repeated_calls_and_keeps_the_bytes(tmp_path,
                                                             monkeypatch):
    """`G` prints the local partial `pow(u, 1.5)` once per unknown: with
    the declarations gcc calls it once, without them twice -- the
    declaration's power, measured on the object -- and both kernels
    return the numpy path's bytes."""
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))
    monkeypatch.setattr(cb, '_loaded', {})
    e = _instance(_RepeatedPow)
    fn = _RepeatedPow._hdl_info['funcs']['G']
    assert fn._csrc.count('pow(') >= 3
    xs = [np.array([v, 0.0]) for v in (0.3, 1.7, 2.0, 1e3)]
    with numpy_backend(_RepeatedPow):
        ref = [e.G(x).copy() for x in xs]

    def pow_calls():
        with c_backend(_RepeatedPow):
            assert _RepeatedPow._hdl_backend_status == 'c'
            assert all(e.G(x).tobytes() == r.tobytes()
                       for x, r in zip(xs, ref))
            so = cb.so_path(fn.__dict__['_hdl_c'].key)
        out = subprocess.run(['objdump', '-d', so], capture_output=True,
                             text=True, check=True).stdout
        return sum('<pow@plt>' in ln and 'call' in ln
                   for ln in out.splitlines())

    merged = pow_calls()
    plain = hdl._KERNEL_C
    for name in _KEPT:
        plain = '\n'.join(ln for ln in plain.splitlines()
                          if not ln.startswith(f'double {name}('))
    assert not _const_declared(plain)
    monkeypatch.setattr(hdl, '_KERNEL_C', plain)
    monkeypatch.setattr(cb, '_loaded', {})
    assert pow_calls() > merged >= 1


## ----------------------------------------------------------------------
## The C backend's defects and gaps, closed before it became the default
## (2026-10-02).

class _TieModel(Behavioural):
    instparams = [Parameter(name='gg', desc='scale', unit='A', default=1.0)]

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        u = var(b.V, 'u')
        return Contribution(b.I, gg / maxc(u, -u) + gg / minc(u, -u))  # noqa: F821


def test_the_prelude_carries_numpys_signed_zero_tie_rule():
    """`_npmax`/`_npmin` return the operand numpy returns on a tie of
    signed zeros -- measured, because it follows the hardware (numpy 2.5 on
    x86: the second)."""
    rule = hdl._numpy_tie_rule()
    assert rule == hdl._TIE_RULE
    assert hdl._NPMAXMIN_C[rule] in hdl._KERNEL_C
    if rule is not None:
        a, b = np.float64(-0.0), 0.0
        r = np.maximum(a, b)
        assert np.signbit(r) == np.signbit(a if rule == 'first' else b)


@needs_cc
def test_a_signed_zero_tie_has_numpys_sign_on_c(tmp_path, monkeypatch):
    """`gg / maxc(u, -u) + gg / minc(u, -u)` at u = +-0: the tie's zero
    decides the infinity's sign.  The prelude returned the FIRST operand
    where numpy returns the second -- +inf against -inf (DEFECT, fixed
    2026-10-02; no byte sweep had met such a tie)."""
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))
    e = _instance(_TieModel)
    xs = [np.array(v) for v in ([0.0, 0.0], [-0.0, 0.0], [0.0, -0.0])]
    with numpy_backend(_TieModel), np.errstate(all='ignore'):
        ref = [e.i(x).copy() for x in xs]
    assert all(np.all(np.isinf(r)) for r in ref)
    with c_backend(_TieModel):
        assert _TieModel._hdl_backend_status == 'c'
        got = [e.i(x) for x in xs]
    assert [g.tobytes() for g in got] == [r.tobytes() for r in ref]


@needs_cc
def test_a_c_bound_class_under_the_jax_backend_runs_its_jax_twin():
    """Under a JAX toolkit `x` is a tracer no C kernel can read: the jax
    twin answers, whatever the class's backend.  The C kernel was consulted
    first and failed on the tracer (DEFECT, fixed 2026-10-02)."""
    jax = pytest.importorskip('jax')
    import jax.numpy as jnp

    from pycircuit.circuit.toolkit import jaxtoolkit
    nodes = [Node(f'n{k}') for k in range(len(eh.EkvNmosHdl.terminals))]
    ref_el = eh.EkvNmosHdl(*nodes)
    ref_el.update_iparv()
    x = 0.1 + 0.2 * np.arange(len(ref_el.nodes), dtype=float)
    with numpy_backend(eh.EkvNmosHdl):
        ref = ref_el.i(x)
    with c_backend(eh.EkvNmosHdl):
        el = eh.EkvNmosHdl(*nodes, toolkit=jaxtoolkit)
        el.update_iparv()
        assert type(el)._hdl_backend_status == 'c'
        got = jax.jit(lambda v: el.i(v))(jnp.asarray(x))
    assert np.allclose(np.asarray(got), ref, rtol=1e-10, atol=0.0)


@needs_cc
def test_a_temperature_array_runs_the_numpy_function(tmp_path):
    """A kernel takes one temperature; an array of them (a temperature
    sweep's epar) is served by the numpy function, which broadcasts --
    the kernel raised on it (until 2026-10-02)."""
    from types import SimpleNamespace as Epar
    e = _instance(eh.DiodeSpiceHdl)
    assert eh.DiodeSpiceHdl._hdl_info['funcs']['i']._clayout[1] is not None
    x = np.array([0.6, 0.0])
    ep = Epar(T=np.array([280.0, 300.0, 320.0]))
    with numpy_backend(eh.DiodeSpiceHdl):
        ref = e.i(x, ep).copy()
        ref_one = e.i(x, Epar(T=300.0)).copy()
    with c_backend(eh.DiodeSpiceHdl):
        assert type(e)._hdl_backend_status == 'c'
        got = e.i(x, ep)
        one = e.i(x, Epar(T=300.0))
    assert got.tobytes() == ref.tobytes()
    assert one.tobytes() == ref_one.tobytes()


@needs_cc
def test_the_kernel_takes_a_read_only_state_and_any_one_number_as_t():
    """The buffer-protocol call (2026-10-02): a read-only `x` is read, not
    written; a temperature that is a numpy scalar or a 0-d array is one
    number and taken; an int is taken; the bytes are the numpy function's
    in every case."""
    from types import SimpleNamespace as Epar
    e = _instance(eh.DiodeSpiceHdl)
    f = type(e)._hdl_info['funcs']['i']
    x = np.array([0.6, 0.0])
    x.setflags(write=False)
    temps = (300.0, 310, np.float64(320.0), np.array(330.0))
    with numpy_backend(eh.DiodeSpiceHdl):
        ref = [e.i(x, Epar(T=T)).copy() for T in temps]
    with c_backend(eh.DiodeSpiceHdl):
        assert type(e)._hdl_backend_status == 'c'
        assert f.__dict__['_hdl_c'] is not None
        for T, r in zip(temps, ref, strict=True):
            got = e.i(x, Epar(T=T))
            assert got.tobytes() == r.tobytes(), T
            assert got is not r and got.flags.writeable
    assert not x.flags.writeable


@needs_cc
def test_without_cffi_a_request_for_c_runs_numpy_and_says_why(tmp_path,
                                                              monkeypatch):
    """The kernels load through cffi; without it a request for C ran into
    an ImportError out of class creation."""
    ## ⚠ THE UN-PIN RUNS AFTER cffi IS VISIBLE AGAIN (2026-10-03): it used to
    ## run in this test's `finally`, while `monkeypatch` still hid cffi (it
    ## is undone at fixture teardown, later), so `_SquareModel` re-resolved
    ## to numpy and stayed there for the module's later tests -- the leak
    ## detector's report; the same teardown-order shape as the one-worker
    ## anomaly (test_hdl_cse's reference swap).
    from pycircuit.circuit.simwarnings import CostWarning
    try:
        with pytest.MonkeyPatch.context() as mp:
            mp.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))
            real = cb.importlib.util.find_spec
            mp.setattr(cb.importlib.util, 'find_spec',
                       lambda n, *a: None if n == 'cffi' else real(n, *a))
            e = _instance(_SquareModel)
            with pytest.warns(CostWarning, match='cffi is not installed'):
                hdl.set_backend('c', _SquareModel)
            assert _SquareModel._hdl_backend_status == \
                'numpy (cffi not installed)'
            assert not _SquareModel._hdl_info.get('_c_bound')
            e.i(np.array([0.5, 0.0]))
    finally:
        hdl.set_backend(None, _SquareModel)


@needs_cc
def test_a_c_bound_class_is_flagged_and_never_fused():
    """`info['_c_bound']` follows attach / detach, and an evaluation
    session hands a C-bound class's calls to its separate kernels."""
    from pycircuit.circuit import _evalhint, _hdl_cse
    e = _instance(eh.MosLevel1Hdl)
    info = type(e)._hdl_info
    x = 0.1 + 0.2 * np.arange(len(e.nodes), dtype=float)
    with c_backend(eh.MosLevel1Hdl):
        assert info['_c_bound'] is True
        with _evalhint.evaluating('i', 'G'):
            assert _hdl_cse.take(e, 'i', x, defaultepar, info,
                                 hdl._args_of) is None
    with numpy_backend(eh.MosLevel1Hdl):
        assert not info.get('_c_bound')


def test_an_exponent_that_is_a_where_value_has_no_c_rendering():
    """numpy raises to a 0-d ARRAY exponent through its fast power path
    (2, 1/2, -1 squared / rooted / reciprocated, not `pow`), so no one C
    form agrees: the printer refuses it, and the class runs numpy."""

    class PowWhere(Behavioural):
        instparams = [Parameter(name='gg', desc='g', unit='S', default=1.0)]

        @staticmethod
        def analog(p, m):
            b = Branch(p, m)
            u = var(b.V, 'u')
            k = var(sympy.Piecewise((2, u > 0), (0.5, True)), 'k')
            return Contribution(b.I, gg * (u * u + 1) ** k)     # noqa: F821

    fn = PowWhere._hdl_info['funcs']['i']
    assert getattr(fn, '_csrc', None) is None
    assert 'exponent that is a numpy.where value' in fn._creason


@needs_cc
def test_a_key_being_built_elsewhere_is_left_to_its_builder(tmp_path,
                                                           monkeypatch):
    """The build lock: while another holder has a key, the parallel
    pre-build leaves it alone (no duplicate `cc`), and `kernel_for` builds
    it once the lock is free."""
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))
    monkeypatch.setattr(cb, '_loaded', {})
    fn = _SquareModel._hdl_info['funcs']['i']
    key = cb.source_key(fn._csrc)
    assert cb.find_so(key) is None
    with cb._key_lock(key) as mine:
        assert mine
        with cb._key_lock(key, blocking=False) as other:
            assert not other
        cb._build_missing_parallel([fn])
        assert cb.find_so(key) is None
    _kern, cold = cb.kernel_for(fn, fn._cshape[0])
    assert cold and cb.find_so(key) is not None


if __name__ == '__main__':                        # pragma: no cover
    sys.exit(pytest.main([__file__, '-v']))
