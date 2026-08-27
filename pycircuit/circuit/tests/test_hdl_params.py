# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The parameter namespace (roadmap S5, second half): ``params_as = 'p'``
hands ``analog()`` a `hdl.ParamNamespace` whose attributes ARE the
parameter symbols.

The claim is BIT-IDENTITY: a model written with ``p`` compiles to the
same generated source, the same ``i``/``G``/``q``/``C``/``CY`` bytes and
the same ``explain()`` text as the same model written with bare names.
Then the three things the namespace exists for -- a helper shared by
several classes, a keyword-named parameter, and a typo caught at compile
time with the class and the declared names in the message -- and the
cache key, which must tell a ``params_as`` class from a bare one whose
``analog()`` source is the same text.

`elements_hdl`'s five ``params_as`` adopters (the three Gummel-Poon
classes, whose helper used to be rebound by hand, and the two SPICE
diodes, whose helper took nineteen arguments) are pinned to hashes
RECORDED from the pre-change source: the same card, the same three
random points, the same ``i/G/q/C`` byte string, the same ``explain()``
text.  The recipe is in `_digest`; the values were taken on 2026-08-26
from commit cb99b18 with the cache disabled, and re-taken from the
changed source before being written down here.  When ``explain()``
grew its ``backend:`` line (the C backend, 2026-08-26) that line was
excluded from the hashed text rather than folded into new digests --
it names the EVALUATION backend, which the suite legitimately varies
-- and the stripped text still hashes to the ORIGINAL recorded
values, so everything this record pins is intact.
"""

import hashlib
import importlib.util
import sys
import textwrap

import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import hdl
from pycircuit.circuit import _hdl_cache as hc
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution, var,
                                   param_given, ParamNamespace, explain)
from pycircuit.utilities.param import Parameter


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    old = cm.default_toolkit
    cm.default_toolkit = numeric
    yield
    cm.default_toolkit = old


def _P(name, default):
    return Parameter(name=name, desc=name, unit='', default=default)


def _make(name, analog, params=(('g', 1e-3),), **extra):
    body = dict(instparams=[_P(n, d) for n, d in params],
                analog=staticmethod(analog))
    body.update(extra)
    return hdl.BehaviouralMeta(name, (Behavioural,), body)


def _bytes(el, x):
    out = []
    for m in ('i', 'G', 'q', 'C', 'CY'):
        f = getattr(el, m, None)
        if f is None:
            continue
        v = f(x) if m != 'CY' else f(x, 1.0)
        out.append(np.ascontiguousarray(np.asarray(v, float)).tobytes())
    return b''.join(out)


def _srcs(cls):
    """The generated source of every compiled function: the let-chain's
    `_src`, or the text lambdify registered with `linecache`."""
    import inspect
    out = {}
    for k, f in cls._hdl_info['funcs'].items():
        if hasattr(f, '_src'):
            out[k] = f._src
        elif callable(f):
            try:
                out[k] = inspect.getsource(f)
            except (OSError, TypeError):
                pass
    return out


## ----------------------------------------------------------------------
## Bit-identity between the two spellings.

class TestSameModelBothWays(object):
    """One model, written bare and written with ``p``; the compiled
    element must not be able to tell which."""

    @staticmethod
    def _eager_pair():
        def bare(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V ** 3 + k * b.V)   # noqa: F821

        def named(p, plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, p.g * b.V ** 3 + p['k'] * b.V)

        params = (('g', 1e-3), ('k', 2e-4))
        return (_make('Eager', bare, params),
                _make('Eager', named, params, params_as='p'))

    @staticmethod
    def _chained_pair():
        def bare(plus, minus):
            b = Branch(plus, minus)
            gg = var(g * k, 'gg')                              # noqa: F821
            sel = sympy.Piecewise((g, param_given('k') > 0.5),  # noqa: F821
                                  (0.0, True))
            return Contribution(b.I, gg * b.V + sel + hdl.ddt(c * b.V))  # noqa

        def named(p, plus, minus):
            b = Branch(plus, minus)
            gg = var(p.g * p.k, 'gg')
            sel = sympy.Piecewise((p.g, p.given('k') > 0.5), (0.0, True))
            return Contribution(b.I, gg * b.V + sel + hdl.ddt(p.c * b.V))

        params = (('g', 1e-3), ('k', 2.0), ('c', 1e-12))
        return (_make('Chained', bare, params),
                _make('Chained', named, params, params_as='p'))

    @pytest.mark.parametrize('pair', ['_eager_pair', '_chained_pair'])
    def test_generated_source_bytes_and_explain_are_identical(self, pair):
        bare, named = getattr(self, pair)()
        assert list(named.terminals) == ['plus', 'minus']
        assert _srcs(bare) and _srcs(bare) == _srcs(named)
        assert bare._hdl_info['paramnames'] == named._hdl_info['paramnames']
        a = bare('a', 'b', k=3.0)
        b = named('a', 'b', k=3.0)
        a.update_iparv()
        b.update_iparv()
        x = np.array([0.37, -0.11])
        assert _bytes(a, x) == _bytes(b, x)
        assert a.explain(maxlines=None) == b.explain(maxlines=None)
        assert explain(bare, maxlines=None) == explain(named, maxlines=None)

    def test_the_comparison_can_fail(self):
        """Mutation check: the identity above is not vacuous.  Swap one
        parameter for another in the ``p`` spelling and every one of the
        three comparisons breaks."""
        bare, _ = self._eager_pair()

        def mutated(p, plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, p.k * b.V ** 3 + p['g'] * b.V)

        named = _make('Eager', mutated, (('g', 1e-3), ('k', 2e-4)),
                      params_as='p')
        assert _srcs(bare) != _srcs(named)
        a = bare('a', 'b')
        b = named('a', 'b')
        a.update_iparv()
        b.update_iparv()
        x = np.array([0.37, -0.11])
        assert _bytes(a, x) != _bytes(b, x)
        assert a.explain(maxlines=None) != b.explain(maxlines=None)

    def test_the_attributes_are_the_symbols(self):
        ns = ParamNamespace('X', ['g', 'lambda'],
                            [hdl._param_symbol('g'),
                             hdl._param_symbol('lambda')])
        assert ns.g is sympy.Symbol('g')
        assert ns['g'] is ns.g
        assert ns['lambda'] is sympy.Symbol('_hdl_kw_lambda')
        assert ns.names == ('g', 'lambda')
        assert list(ns) == ['g', 'lambda'] and len(ns) == 2
        assert 'g' in ns and 'gg' not in ns
        assert 'X' in repr(ns) and 'lambda' in repr(ns)
        assert ns.given('g') == param_given('g')
        with pytest.raises(AttributeError, match='read-only'):
            ns.g = 2
        ## numpy's and sympy's duck-type probes must see a plain
        ## AttributeError, not the model-facing message.
        assert not hasattr(ns, '__array__')
        assert not hasattr(ns, '_sympy_')

    def test_a_pin_may_still_be_called_p(self):
        """Without ``params_as`` the first argument is a terminal
        whatever it is called; ``analog(p, n)`` is a spelling several
        test models use, which is why opting in is EXPLICIT rather than
        keyed on the name."""
        def analog(p, n):
            b = Branch(p, n)
            return Contribution(b.I, g * b.V)                  # noqa: F821

        cls = _make('PinP', analog)
        assert list(cls.terminals) == ['p', 'n']
        el = cls('a', 'b', g=2e-3)
        el.update_iparv()
        assert np.allclose(el.i(np.array([1.0, 0.0])), [2e-3, -2e-3])


## ----------------------------------------------------------------------
## What the namespace makes possible.

def _shared_helper(p, plus, minus):
    """A module-level helper two classes call with different cards; it
    reads the card through ``p`` and has no idea which class called."""
    b = Branch(plus, minus)
    gg = var(p.g * p.scale, 'gg')
    return Contribution(b.I, gg * b.V)


class TestAHelperSharedByTwoClasses(object):

    def test_each_class_sees_its_own_card(self):
        def analog_a(p, plus, minus):
            return _shared_helper(p, plus, minus)

        def analog_b(p, plus, minus):
            return _shared_helper(p, plus, minus)

        A = _make('A', analog_a, (('g', 1e-3), ('scale', 2.0)),
                  params_as='p')
        B = _make('B', analog_b, (('g', 5e-3), ('scale', 10.0)),
                  params_as='p')
        a = A('a', 'b')
        b = B('a', 'b')
        a.update_iparv()
        b.update_iparv()
        x = np.array([1.0, 0.0])
        assert np.allclose(a.i(x), [2e-3, -2e-3])
        assert np.allclose(b.i(x), [5e-2, -5e-2])
        ## And the card given on the instance, not just the default.
        c = B('a', 'b', scale=1.0)
        c.update_iparv()
        assert np.allclose(c.i(x), [5e-3, -5e-3])

    def test_the_helper_cannot_be_reached_bare(self):
        """The control: the same helper written with bare names does
        not resolve them, because its ``__globals__`` are this module's
        and the compiler injects into a private copy of analog()'s."""
        def helper(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                  # noqa: F821

        def analog(plus, minus):
            return helper(plus, minus)

        with pytest.raises(NameError, match=r"Bare\.analog\(\) uses the "
                                            r"name 'g'.*helper"):
            _make('Bare', analog)

    def test_the_library_helpers_take_p(self):
        """`_gp_core` and `_spice_diode` are the two helpers this was
        built for; ``_with_params`` is gone."""
        import inspect
        assert inspect.getfullargspec(eh._gp_core).args[0] == 'p'
        ## `junction` (2026-08-26, fifth batch): the optical devices'
        ## hook for an extra current across the helper's own junction
        ## branch, optional and None by default.
        assert inspect.getfullargspec(eh._spice_diode).args == \
            ['p', 'a', 'c', 'T', 'junction']
        assert not hasattr(eh, '_with_params')
        for name in ('GummelPoonNpnHdl', 'GummelPoonPnpHdl',
                     'GummelPoonNpnThermalHdl', 'DiodeSpiceHdl',
                     'DiodeSpiceThermalHdl'):
            assert getattr(eh, name).params_as == 'p'


class TestKeywordNamedParameters(object):

    def test_lambda_and_as_are_reachable_through_p(self):
        def analog(p, plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, p.g * (1 + p['lambda'] * p['as'])
                                * b.V)

        cls = _make('Kw', analog, (('g', 1e-3), ('lambda', 0.5),
                                   ('as', 2.0)), params_as='p')
        assert cls._hdl_info['paramnames'] == ['g', 'lambda', 'as']
        el = cls('a', 'b', **{'lambda': 1.0})
        el.update_iparv()
        x = np.array([1.0, 0.0])
        assert np.allclose(el.i(x), [3e-3, -3e-3])
        assert np.allclose(el.G(x)[0, 0], 3e-3)
        ## The generated code cannot spell `lambda`; the symbol is
        ## mangled, and explain() shows the mangled name in the source
        ## and the declared name in the parameter list.
        text = explain(cls, maxlines=None)
        assert '  lambda\n' in text and '_hdl_kw_lambda' in text
        assert 'def _lambdifygenerated(x, g, _hdl_kw_lambda, _hdl_kw_as' \
            in text

    def test_lambda_in_a_let_chain(self):
        """The chained printer emits argument names VERBATIM, which is
        where an unmangled keyword would have been a SyntaxError."""
        def analog(p, plus, minus):
            b = Branch(plus, minus)
            clm = var(1 + p['lambda'] * b.V, 'clm')
            return Contribution(b.I, p.g * clm * b.V)

        cls = _make('KwChain', analog, (('g', 1e-3), ('lambda', 0.5)),
                    params_as='p')
        assert cls._hdl_info['chained']
        el = cls('a', 'b', **{'lambda': 2.0})
        el.update_iparv()
        x = np.array([1.0, 0.0])
        assert np.allclose(el.i(x), [3e-3, -3e-3])

    def test_the_bare_style_refuses_a_keyword_name(self):
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                  # noqa: F821

        with pytest.raises(TypeError,
                           match=r"KwBare declares the parameter 'lambda', "
                                 r"which cannot be read as bare name inside "
                                 r"analog\(\): a Python keyword is not a "
                                 r"valid name\. Set params_as = 'p' on the "
                                 r"class and read it as p\['lambda'\]"):
            _make('KwBare', analog, (('g', 1e-3), ('lambda', 0.5)))
        with pytest.raises(TypeError, match=r"parameters 'lambda', 'as'"):
            _make('KwBare2', analog, (('g', 1e-3), ('lambda', 0.5),
                                      ('as', 1.0)))
        with pytest.raises(TypeError, match=r"'my-p'.*not a valid Python "
                                            r"identifier"):
            _make('KwBare3', analog, (('g', 1e-3), ('my-p', 0.5)))

    def test_mos_level1_accepts_spice_spellings_as_aliases(self):
        """`MosLevel1Hdl` keeps ``lambd``/``asrc`` as its canonical names
        (it reads them bare); SPICE's ``LAMBDA`` and ``AS`` are aliases
        on the instance.  `Parameter` itself has no alias support --
        `Behavioural.aliasparams` is the mechanism."""
        for cls in (eh.MosLevel1Hdl, eh.MosLevel1PmosHdl):
            assert cls.aliasparams == {'lambda': 'lambd', 'as': 'asrc'}
            a = cls('d', 'g', 's', 'b', **{'lambda': 0.02, 'as': 3e-12})
            b = cls('d', 'g', 's', 'b', lambd=0.02, asrc=3e-12)
            a.update_iparv()
            b.update_iparv()
            assert a.iparv.lambd == b.iparv.lambd == 0.02
            assert a.iparv.asrc == b.iparv.asrc == 3e-12
            with pytest.raises(ValueError, match='alias'):
                cls('d', 'g', 's', 'b', lambd=0.02, **{'lambda': 0.02})


class TestATypoIsCaughtAtCompileTime(object):

    MSG = (r"Typo\.analog\(\) reads parameter 'gg', which Typo does not "
           r"declare \(its instparams are 'g', 'k'\)\. Declare it: "
           r"instparams = \[\.\.\., Parameter\(name='gg', desc=\.\.\., "
           r"unit=\.\.\., default=\.\.\.\)\]\.")

    def test_attribute_item_and_given(self):
        params = (('g', 1e-3), ('k', 2.0))

        def attr(p, plus, minus):
            return Contribution(Branch(plus, minus).I,
                                p.gg * Branch(plus, minus).V)

        def item(p, plus, minus):
            return Contribution(Branch(plus, minus).I,
                                p['gg'] * Branch(plus, minus).V)

        def given(p, plus, minus):
            return Contribution(Branch(plus, minus).I,
                                p.given('gg') * Branch(plus, minus).V)

        with pytest.raises(AttributeError, match=self.MSG):
            _make('Typo', attr, params, params_as='p')
        with pytest.raises(KeyError, match=self.MSG):
            _make('Typo', item, params, params_as='p')
        with pytest.raises(KeyError, match=self.MSG):
            _make('Typo', given, params, params_as='p')

    def test_the_control_compiles(self):
        """Mutation check for the test above: the same three bodies with
        the name spelled right compile and evaluate."""
        params = (('g', 1e-3), ('k', 2.0))

        def ok(p, plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, (p.g + p['g'] * p.given('k')) * b.V)

        cls = _make('Typo', ok, params, params_as='p')
        el = cls('a', 'b', k=1.0)
        el.update_iparv()
        assert np.allclose(el.i(np.array([1.0, 0.0])), [2e-3, -2e-3])

    def test_a_bare_name_in_a_params_as_class_says_where_it_went(self):
        def analog(p, plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                  # noqa: F821

        with pytest.raises(NameError,
                           match=r"BareInNs\.analog\(\) uses the name 'g'.*"
                                 r"sets params_as = 'p'.*NOT bound as bare "
                                 r"names -- they are read as p\.g or "
                                 r"p\['g'\]"):
            _make('BareInNs', analog, params_as='p')

    def test_the_declaration_is_checked(self):
        def no_p(plus, minus):
            return Contribution(Branch(plus, minus).I, 0.0)

        with pytest.raises(TypeError,
                           match=r"WrongFirst sets params_as = 'p', so "
                                 r"analog\(\)'s FIRST argument must be named "
                                 r"'p'.*it is analog\(plus, minus\)\. Write "
                                 r"\"def analog\(p, plus, minus\)\""):
            _make('WrongFirst', no_p, params_as='p')

        def only_p(p):
            return None

        with pytest.raises(TypeError,
                           match=r"NoTerms\.analog\(\) declares no terminal "
                                 r"arguments.*\"def analog\(p, plus, "
                                 r"minus\)\""):
            _make('NoTerms', only_p, params_as='p')
        with pytest.raises(TypeError,
                           match=r"BadName sets params_as = 'lambda', which "
                                 r"is not a Python identifier"):
            _make('BadName', only_p, params_as='lambda')

        def selfy(p, self, minus):
            return None

        with pytest.raises(TypeError, match=r"'self'.*FIRST TERMINAL"):
            _make('Selfy', selfy, params_as='p')

    def test_any_identifier_may_name_the_namespace(self):
        def analog(par, plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, par.g * b.V)

        cls = _make('Par', analog, params_as='par')
        assert list(cls.terminals) == ['plus', 'minus']
        el = cls('a', 'b', g=4e-3)
        el.update_iparv()
        assert np.allclose(el.i(np.array([1.0, 0.0])), [4e-3, -4e-3])

    def test_given_is_param_given(self):
        """`p.given(x)` is `param_given(x)` with the name checked: the
        symbol is the same object, so the runtime is the same path."""
        seen = {}

        def analog(p, plus, minus):
            seen['g'] = p.given('k')
            b = Branch(plus, minus)
            return Contribution(b.I, sympy.Piecewise(
                (p.k, p.given('k') > 0.5), (p.g, True)) * b.V)

        cls = _make('Given', analog, (('g', 1e-3), ('k', 5e-3)),
                    params_as='p')
        assert seen['g'] == param_given('k')
        assert cls._hdl_info['given_names'] == ['k']
        a = cls('a', 'b')
        b = cls('a', 'b', k=5e-3)
        a.update_iparv()
        b.update_iparv()
        x = np.array([1.0, 0.0])
        assert np.allclose(a.i(x), [1e-3, -1e-3])
        assert np.allclose(b.i(x), [5e-3, -5e-3])

    def test_a_collapse_variant_keeps_the_namespace(self):
        """Collapse variants are compiled subclasses; `params_as` is
        inherited and the variant's analog() is handed ``p`` again."""
        def analog(p, plus, minus):
            mid = hdl.Node('mid')
            b1, b2 = Branch(plus, mid), Branch(mid, minus)
            return (Contribution(b1.I, b1.V / p.r1),
                    Contribution(b2.I, b2.V / p.r2),
                    hdl.Collapse(b1, p.r1 <= 0))

        cls = _make('Coll', analog, (('r1', 0.0), ('r2', 1e3)),
                    params_as='p')
        a = cls('a', 'b')
        a.update_iparv()
        assert type(a) is not cls and type(a).params_as == 'p'
        assert np.allclose(a.i(np.array([1.0, 0.0])), [1e-3, -1e-3])
        b = cls('a', 'b', r1=1e3)
        b.update_iparv()
        assert np.allclose(b.i(np.array([1.0, 0.0, 0.5])),
                           [5e-4, -5e-4, 0.0])


## ----------------------------------------------------------------------
## The cache key.

_MODULE = '''
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, var
from pycircuit.utilities.param import Parameter


class Probe(Behavioural):
    params_as = 'p'
    instparams = [Parameter(name='r', default=%(default)s)]

    @staticmethod
    def analog(p, plus, minus):
        b = Branch(plus, minus)
        g = var(%(scale)s / p.r, 'g')
        return Contribution(b.I, b.V * g)
'''


@pytest.fixture
def cache_dir(tmp_path, monkeypatch):
    d = tmp_path / 'hdlcache'
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(d))
    monkeypatch.delenv('PYCIRCUIT_HDL_CACHE', raising=False)
    monkeypatch.setattr(hc, 'ENABLED', True)
    return d


def _load_probe(tmp_path, tag, default='1e3', scale='1.0'):
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


class TestTheCacheKey(object):

    def test_same_source_hits_and_is_identical(self, cache_dir, tmp_path):
        a = _load_probe(tmp_path, 'a')
        b = _load_probe(tmp_path, 'b')
        assert a._hdl_cache_status == 'miss'
        assert b._hdl_cache_status == 'hit'
        ea, eb = a('x', 'y', r=2e3), b('x', 'y', r=2e3)
        ea.update_iparv()
        eb.update_iparv()
        x = np.array([0.3, -0.2])
        assert _bytes(ea, x) == _bytes(eb, x)
        assert ea.explain(maxlines=None) == eb.explain(maxlines=None)

    def test_body_and_default_changes_miss(self, cache_dir, tmp_path):
        _load_probe(tmp_path, 'a')
        assert _load_probe(tmp_path, 'b', scale='2.0')._hdl_cache_status \
            == 'miss'
        assert _load_probe(tmp_path, 'c', default='1e4')._hdl_cache_status \
            == 'miss'

    def test_params_as_is_in_the_key(self):
        """Two classes with the SAME analog() source text, one bare and
        one `params_as`, must not share an entry: the bare one's first
        argument is a terminal, so the compiled element has a different
        pin count.  Keyed on plain classes, because the bare one does
        not compile."""
        src = textwrap.dedent('''
            def analog(p, plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, p.g * b.V)
        ''')
        ns = dict(Branch=Branch, Contribution=Contribution)
        exec(compile(src, __file__, 'exec'), ns)
        body = dict(instparams=[_P('g', 1e-3)], analog=staticmethod(ns['analog']),
                    __module__=__name__)
        bare = type('Same', (), dict(body))
        named = type('Same', (), dict(body, params_as='p'))
        assert hc.key_for(bare) != hc.key_for(named)
        assert hc.key_for(named) == hc.key_for(
            type('Same', (), dict(body, params_as='p')))


## ----------------------------------------------------------------------
## The library adopters, against hashes recorded before the change.

_BJT = dict(IS=2e-15, bf=150.0, br=4.0, nf=1.02, nr=1.05, vaf=80.0, var=12.0,
            ikf=0.05, ikr=0.01, ise=3e-14, ne=1.6, isc=1e-13, nc=1.8,
            rb=120.0, rbm=20.0, re=0.8, rc=2.5, cje=1.2e-12, vje=0.75,
            mje=0.35, cjc=0.5e-12, vjc=0.6, mjc=0.3, xcjc=0.7, tf=3e-10,
            xtf=2.0, vtf=5.0, itf=0.02, tr=5e-8, xtb=1.5, eg=1.12, xti=3.2,
            tnom=300.15, area=1.3, rth=150.0, cth=1e-4)
_DIO = dict(IS=3e-15, rs=1.7, n=1.05, tt=4e-9, cjo=2e-12, vj=0.7, m=0.4,
            eg=1.1, xti=3.0, fc=0.5, bv=50.0, ibv=1e-4, kf=1e-15, af=1.1,
            area=2.0, tnom=300.15, rth=150.0, cth=1e-4)

#: (class, card, pins, x-length, three point digests, explain digest),
#: recorded from cb99b18's `elements_hdl` -- the `_with_params` /
#: nineteen-argument spelling -- with the cache disabled.
_RECORDED = [
    ## ⚠ Re-recorded again 2026-08-27 (roadmap sec. 43) for the thermal
    ## node and the diode's series branch becoming probes.  Same story
    ## as below: EXPLAIN text only, POINT digests untouched, and
    ## `i`/`G`/`q`/`C` measured bit-identical across every affected
    ## model and card.
    ##
    ## ⚠ The three Gummel-Poon EXPLAIN digests were re-recorded
    ## 2026-08-27 (roadmap sec. 42): the base-resistance branch is now a
    ## `limit_identity` PCNR probe, so `explain()` lists one more probe.
    ##
    ## Their POINT digests are UNCHANGED, and that is the evidence that
    ## matters -- this test separates the two on purpose, and only the
    ## text moved.  Independently measured: `i`/`G`/`q`/`C` are
    ## BIT-identical across 25 random biases for all three classes, max
    ## relative difference exactly 0.000e+00.  An identity probe applies
    ## no law, so it must not move a number, and it did not.
    ('GummelPoonNpnHdl', _BJT, ('c', 'b', 'e'), 6,
     ['4d746c561ac7b430', '924f3d39a2c404df', '5ea1b16f597dd226'],
     '3b86779394634e8e'),
    ('GummelPoonPnpHdl', _BJT, ('c', 'b', 'e'), 6,
     ['2c024addadbe4ba4', '38dfed3011c5a3f5', 'e21e0a40a85e01bf'],
     '2cdbb30b8f2ec50d'),
    ('GummelPoonNpnThermalHdl', _BJT, ('c', 'b', 'e', 'th', 'tha'), 8,
     ['d2bcfb160947baec', 'fde3cc1f8ae9d6db', '943ddcd06908762d'],
     '651b1d49cd49caf9'),
    ## ⚠ The two SPICE-diode `explain` digests were RE-RECORDED
    ## 2026-08-27 for `_autohold` (roadmap sec. 36): the regularisers now
    ## hold their own arguments, so the chain carries more named
    ## intermediates and `explain()` lists them.
    ##
    ## `DiodeSpiceHdl`'s three POINT digests are UNCHANGED -- its numbers
    ## did not move at all, only the text.  `DiodeSpiceThermalHdl`'s did
    ## move, and were measured before re-recording: 2880 samples, 92 of
    ## 920 comparable differ, **max relative 1.213e-15, median 1.568e-16**
    ## -- reassociation at the last bit, because holding stops sympy
    ## flattening across the boundary.
    ##
    ## Old explain digests: DiodeSpiceHdl 4faadd376cafc95a,
    ## DiodeSpiceThermalHdl b2230da8487...dbfb58d3a752 (see git).
    ('DiodeSpiceHdl', _DIO, ('a', 'c'), 3,
     ['eadd075b306b20e0', 'ac934e5f1d3b6a30', 'de1b14f3dc3cb89c'],
     '149445302b3077ef'),
    ## Old: points b2230da8487a876f / cda54034e9a55933 / 8361dbfb58d3a752,
    ## explain 0dcf03f6b937d283.  Moved by 1.213e-15 relative -- see the
    ## note above.
    ('DiodeSpiceThermalHdl', _DIO, ('a', 'c', 'th', 'tha'), 5,
     ['2779caa138b9b423', 'a630e803c927764a', '5452a4a0006fbacb'],
     'f2ec390e8b1e6ad9'),
]


def _digest(data):
    return hashlib.sha256(data).hexdigest()[:16]


def _points(el, n):
    """Three points, `RandomState(7)`, `0.7*uniform(-1, 1, n)` each;
    the digest is over the concatenated `i, G, q, C` bytes."""
    rng = np.random.RandomState(7)
    out = []
    for _ in range(3):
        x = 0.7 * rng.uniform(-1, 1, n)
        out.append(_digest(b''.join(
            np.ascontiguousarray(np.asarray(getattr(el, m)(x),
                                            float)).tobytes()
            for m in ('i', 'G', 'q', 'C'))))
    return out


@pytest.mark.parametrize('name,card,pins,n,points,expl', _RECORDED,
                         ids=[r[0] for r in _RECORDED])
def test_library_adopter_is_bit_identical_to_the_record(name, card, pins, n,
                                                        points, expl):
    cls = getattr(eh, name)
    declared = [p.name for p in cls.instparams]
    el = cls(*pins, **{k: v for k, v in card.items() if k in declared})
    el.update_iparv()
    assert len(hdl.x_layout(el)) == n
    assert _points(el, n) == points
    ## The `backend:` line reports which EVALUATION backend is bound
    ## (numpy or C, per $PYCIRCUIT_HDL_BACKEND) -- deliberately outside
    ## this record, which pins the COMPILE: the suite runs in both
    ## backend states and the compile is the same in both.
    text = '\n'.join(ln for ln in el.explain(maxlines=None).splitlines()
                     if not ln.startswith('backend:'))
    assert _digest(text.encode()) == expl


def test_the_record_can_fail():
    """Mutation check: a one-parameter change to the card moves every
    digest, so an element that silently read a different parameter
    (the failure a mis-bound helper produces) would not match."""
    cls = eh.GummelPoonNpnHdl
    card = dict(_BJT, bf=151.0)
    declared = [p.name for p in cls.instparams]
    el = cls('c', 'b', 'e', **{k: v for k, v in card.items()
                              if k in declared})
    el.update_iparv()
    assert _points(el, 6) != _RECORDED[0][4]
