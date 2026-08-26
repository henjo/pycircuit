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
                                   Collapse, var, maxc)
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
    """`'equal'`, `'zero-sign'` (values equal, bytes differ only in the
    sign of exact zeros) or `'value'`."""
    ref = np.asarray(ref, float)
    out = np.asarray(out, float)
    if ref.tobytes() == out.tobytes():
        return 'equal'
    eq = (ref == out) | (np.isnan(ref) & np.isnan(out))
    return 'zero-sign' if bool(eq.all()) else 'value'


def _sweep(e, cls, pts):
    """{func name: {'equal': n, 'zero-sign': n, 'value': n}} comparing
    the C kernel against the numpy function it was printed from."""
    args = [float(v) for v in hdl._args_of(e, defaultepar)]
    out = {}
    for name, f in _c_funcs(cls):
        kern = f.__dict__['_hdl_c']
        tally = {'equal': 0, 'zero-sign': 0, 'value': 0}
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
        with np.errstate(all='ignore'):
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
            with np.errstate(all='ignore'):
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
        ref = e.i(x).copy()
        monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))
        monkeypatch.setattr(cb, 'CFLAGS', tuple(
            f for f in cb.CFLAGS if f != '-fno-builtin-pow'))
        monkeypatch.setattr(cb, '_loaded', {})
        with c_backend(_SquareModel):
            assert _SquareModel._hdl_backend_status == 'c'
            got = e.i(x)
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
        with np.errstate(all='ignore'):
            ref = e.G(x).copy()
        with c_backend(cls):
            assert var_cls._hdl_backend_status == 'c', \
                var_cls._hdl_backend_status
            with np.errstate(all='ignore'):
                got = e.G(x)
        assert ref.tobytes() == got.tobytes()
        ## unpinned: the variant follows whatever the environment says,
        ## exactly as the base does
        assert var_cls._hdl_backend_status == cls._hdl_backend_status

    def test_results_identical_through_the_element_methods(self):
        """The same instance, backend toggled around it: `i/G/q/C`
        bytes unchanged, and the packed-parameter cache follows a
        parameter update."""
        cls = eh.GummelPoonNpnHdl
        e = _instance(cls, **KW['GummelPoonNpnHdl'])
        n = len(hdl.x_layout(cls))
        x = np.linspace(-0.4, 0.4, n)
        with np.errstate(all='ignore'):
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
        with np.errstate(all='ignore'):
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
        broken = hdl._KERNEL_C.replace(
            'return (a >= b || a != a) ? a : b;', 'return (a >= b) ? a : b;')
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
        ref = e.i(np.array([0.5, 0.0])).copy()
        monkeypatch.setattr(hdl, '_KERNEL_C',
                            hdl._KERNEL_C + 'this is not C\n')
        with pytest.warns(UserWarning, match='build failed'):
            hdl.set_backend('c', M)
        try:
            assert M._hdl_backend_status.startswith(
                'numpy (compile failed:')
            got = e.i(np.array([0.5, 0.0]))
            assert got.tobytes() == ref.tobytes()
        finally:
            hdl.set_backend(None, M)


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
            res = Transient(build(), toolkit=numeric).solve(
                tend=2e-6, timestep=2e-8, fixed_timestep=True)
            return np.asarray(res.v('out', gnd), float)

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
            assert t['value'] == 0 and t['zero-sign'] == 0, (k, t)
        assert tallies['C']['value'] == 0, tallies['C']


if __name__ == '__main__':                        # pragma: no cover
    sys.exit(pytest.main([__file__, '-v']))
