# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The HDL's syntactic surface, its parameter namespace, and the two
instruments (`explain`, `check_jacobians`).

Three of the mistakes checked here used to be **silent** -- an argument
named ``self`` became a terminal, a class-body ``terminals`` was
discarded, an internal ``Node('out')`` merged into the ``out`` pin -- and
a silent wrong answer is worse than any message.  So each of those tests
does two things: it asserts the error now fires, and it demonstrates the
mechanism that used to swallow it *still works for every other name*, so
the check is the only thing standing between the author and a differently
wired element.

The rest pin message CONTENT with ``match=``, not just the exception
type: an error that names no class and offers no fix is the failure this
whole file exists to remove, and only the text distinguishes the two.

Every check is paired with a control -- the same model written correctly,
compiling and evaluating -- because a test that only asserts a raise
passes just as well on a compiler that refuses everything.
"""

import math
import re
import types

import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import hdl
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.hdl import (Behavioural, Branch, Node, Contribution,
                                   ddt, idt, var, explain, check_jacobians,
                                   x_layout)
from pycircuit.utilities.param import Parameter


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    """Every element here is evaluated numerically.  Set centrally, so a
    test added later cannot forget it."""
    old = cm.default_toolkit
    cm.default_toolkit = numeric
    yield
    cm.default_toolkit = old


def _P(name, default):
    return Parameter(name=name, desc=name, unit='', default=default)


def _make(name, analog, params=(('g', 1e-3),), **extra):
    """Compile a class from a plain function, so a model that must FAIL
    to compile can live inside a test rather than at import time."""
    body = dict(instparams=[_P(n, d) for n, d in params],
                analog=staticmethod(analog))
    body.update(extra)
    return hdl.BehaviouralMeta(name, (Behavioural,), body)


## A module global, used to prove that ordinary names still resolve
## inside analog() after the parameter injection stopped mutating this
## module's dict.
SCALE = 3.0


class TestTheSyntacticChecks(object):

    def test_a_correct_model_compiles_and_evaluates(self):
        """The control for everything below.  A compiler that refused
        every input would pass all the `raises` tests in this class."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                    # noqa: F821

        el = _make('Good', analog)('a', 'b', g=2e-3)
        el.update_iparv()
        assert list(el.terminals) == ['plus', 'minus']
        assert np.allclose(el.i(np.array([1.0, 0.0])), [2e-3, -2e-3])

    def test_forgetting_the_return_names_the_class_and_the_fix(self):
        """Verilog-A has no `return`, so this is the first thing a
        newcomer writes.  It used to be `TypeError: 'NoneType' object is
        not iterable`, raised from inside the compiler."""
        def analog(plus, minus):
            Contribution(Branch(plus, minus).I,
                         g * Branch(plus, minus).V)              # noqa: F821

        with pytest.raises(TypeError, match=r'NoReturn\.analog\(\) returned '
                                            r'None.*must RETURN'):
            _make('NoReturn', analog)

    def test_a_bare_expression_is_refused(self):
        def analog(plus, minus):
            return g * Branch(plus, minus).V                     # noqa: F821

        with pytest.raises(TypeError,
                           match=r'Bare\.analog\(\) returned .*not a '
                                 r'contribution statement'):
            _make('Bare', analog)

    def test_a_self_argument_is_refused_rather_than_made_a_terminal(self):
        """It used to become the element's FIRST TERMINAL, shifting every
        connection by one pin, without a word.

        The second half of this test is the evidence that the silence is
        gone rather than merely relocated: rename that first argument to
        anything else and the compiler DOES make it a terminal, quietly,
        exactly as it used to for `self`.
        """
        def analog(self, plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                    # noqa: F821

        with pytest.raises(TypeError,
                           match=r"WithSelf\.analog\(\) takes 'self' as its "
                                 r'first argument.*FIRST TERMINAL'):
            _make('WithSelf', analog)

        def analog2(bulk, plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                    # noqa: F821

        silent = _make('ThreePin', analog2)
        assert list(silent.terminals) == ['bulk', 'plus', 'minus']

    def test_a_class_body_terminals_that_disagrees_is_refused(self):
        """`terminals` in the class body is DISCARDED -- the metaclass
        overwrites it with analog()'s argument names.  Refused only when
        the two disagree, which is the case where the discarding changes
        the element; a declaration that restates the signature is
        redundant, not wrong, and existing models carry one."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                    # noqa: F821

        with pytest.raises(TypeError,
                           match=r'Reordered declares terminals = .*'
                                 r'DISCARDED'):
            _make('Reordered', analog, terminals=('minus', 'plus'))

        same = _make('Restated', analog, terminals=('plus', 'minus'))
        assert list(same.terminals) == ['plus', 'minus']

    def test_an_internal_node_that_collides_with_a_terminal_is_refused(self):
        """Nodes are identified by name, so `Node('out')` in an element
        with an `out` terminal used to become the terminal -- one unknown
        fewer and the internal equation added to the pin's.

        The renamed control is what shows the size of what was lost: the
        same model with a non-colliding name has THREE unknowns, so the
        collision was not a naming nicety, it deleted a node.
        """
        def analog(plus, out):
            mid = Node('out')
            return (Contribution(Branch(plus, mid).I, g * Branch(plus, mid).V),
                    Contribution(Branch(mid, out).I,
                                 g * Branch(mid, out).V))        # noqa: F821

        with pytest.raises(ValueError,
                           match=r"Collide\.analog\(\) builds its own "
                                 r"Node\('out'\).*already a terminal"):
            _make('Collide', analog)

        def analog2(plus, out):
            mid = Node('mid')
            return (Contribution(Branch(plus, mid).I, g * Branch(plus, mid).V),
                    Contribution(Branch(mid, out).I,
                                 g * Branch(mid, out).V))        # noqa: F821

        el = _make('Renamed', analog2)('a', 'b')
        el.update_iparv()
        assert el.n == 3
        assert [e[2] for e in x_layout(el)] == ['terminal', 'terminal',
                                                'internal node']

    def test_a_python_if_on_a_circuit_quantity_names_the_class(self):
        """Raw sympy said `cannot determine truth value of Relational`
        and nothing about where it came from or what to write."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            if b.V > 0:
                return Contribution(b.I, g * b.V)                # noqa: F821
            return Contribution(b.I, 0)

        with pytest.raises(TypeError,
                           match=r'IfOnV\.analog\(\).*where Python wanted a '
                                 r'bool.*Piecewise'):
            _make('IfOnV', analog)

        def analog2(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, sympy.Piecewise((g * b.V, b.V > 0),
                                                     (0.0, True)))  # noqa
        el = _make('PiecewiseOnV', analog2)('a', 'b', g=1e-3)
        el.update_iparv()
        assert np.allclose(el.i(np.array([1.0, 0.0])), [1e-3, -1e-3])
        assert np.allclose(el.i(np.array([-1.0, 0.0])), [0.0, 0.0])

    @pytest.mark.parametrize('name,fn', [('math', math.exp),
                                         ('numpy', np.exp)])
    def test_a_float_math_function_on_a_quantity_names_the_class(self, name,
                                                                 fn):
        """`math.exp` raised `Cannot convert expression to float`;
        `numpy.exp` raised a ufunc message about a `Mul` with no callable
        `exp` method.  Two different texts, one mistake, neither naming
        the model or the replacement."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * fn(b.V / 0.026))        # noqa: F821

        with pytest.raises(TypeError,
                           match=r'FloatMath\.analog\(\): a math or numpy '
                                 r'function.*sympy\.exp'):
            _make('FloatMath', analog)

        def analog2(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * sympy.exp(b.V / 0.026))  # noqa: F821

        el = _make('SympyMath' + name, analog2)('a', 'b', g=1e-3)
        el.update_iparv()
        assert np.isfinite(np.asarray(el.i(np.array([0.1, 0.0])), float)).all()

    def test_star_args_cannot_become_terminals(self):
        def analog(*terminals):
            b = Branch(terminals[0], terminals[1])
            return Contribution(b.I, g * b.V)                    # noqa: F821

        with pytest.raises(TypeError,
                           match=r'Star\.analog\(\) declares \*terminals'):
            _make('Star', analog)

    def test_the_wrong_number_of_connection_nodes_is_refused(self):
        """`Circuit.__init__` zips terminals against the positional
        arguments, and zip is silent about a length mismatch."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                    # noqa: F821

        cls = _make('TwoPin', analog)
        with pytest.raises(TypeError,
                           match=r'TwoPin has 2 terminals \(plus, minus\) but '
                                 r'was given 1 connection node'):
            cls('a')
        with pytest.raises(TypeError, match=r'was given 3 connection nodes'):
            cls('a', 'b', 'c')
        ## Zero is still legal: that is the add_instance(name, el, plus=n)
        ## form, where the connection is made by name afterwards.
        assert cls().n == 2

    def test_an_unknown_parameter_names_the_class_and_the_valid_names(self):
        """It used to be a bare
        `KeyError: 'parameter R not in parameter dictionary'`."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, b.V / r)                    # noqa: F821

        cls = _make('Res', analog, params=(('r', 1e3),))
        with pytest.raises(TypeError,
                           match=r"Res has no parameter 'R'\. Parameter names "
                                 r"are case sensitive: did you mean 'r'\?"):
            cls('a', 'b', R=1e3)
        with pytest.raises(TypeError,
                           match=r"Res has no parameter 'rr'\. It accepts: r"):
            cls('a', 'b', rr=1e3)
        with pytest.raises(TypeError,
                           match=r"'plus' is a TERMINAL of this element"):
            cls('a', 'b', plus=1.0)
        el = cls('a', 'b', r=1e3)
        el.update_iparv()
        assert np.allclose(el.i(np.array([1.0, 0.0])), [1e-3, -1e-3])


class TestTheParameterNamespace(object):
    """Item B: the injection used to be
    ``analogfunc.__globals__.update(...)``, and `copy()` of a function
    returns the same object, so ``__globals__`` IS the defining module's
    dict."""

    def test_compiling_a_model_does_not_write_into_this_module(self):
        """The leak, reproduced in the smallest form that shows it: a
        parameter named `vt` used to leave `Symbol('vt')` bound in the
        module that defined the class, permanently, shadowing anything
        of that name."""
        assert 'vt' not in globals()

        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, IS * (sympy.exp(b.V / vt) - 1))  # noqa

        _make('LeakProbe', analog, params=(('IS', 1e-14), ('vt', 0.026)))
        assert 'vt' not in globals()
        assert 'IS' not in globals()
        assert isinstance(hdl.vt, types.FunctionType)

    def test_elements_hdl_keeps_the_helper_the_leak_used_to_replace(self):
        """`elements_hdl.DiodeHdl` declares a parameter named `vt`, and
        `elements_hdl` imports `hdl.vt`, the thermal-voltage FUNCTION.
        After importing the module, `vt` in it used to be a Symbol."""
        assert not isinstance(getattr(eh, 'vt', None), sympy.Symbol)
        assert eh._vt is hdl.vt
        assert callable(eh._vt)

    def test_a_parameter_of_another_class_is_not_in_scope(self):
        """The expensive half of the leak: class B could use a name only
        class A declared, compile without a word, and fail at the first
        evaluation as a bare `NameError` from inside generated code.  Now
        it fails at COMPILE time, naming B and the name."""
        def analog_a(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, gain * b.V)                 # noqa: F821

        _make('Donor', analog_a, params=(('gain', 2e-3),))

        def analog_b(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, gain * b.V)                 # noqa: F821

        with pytest.raises(NameError,
                           match=r"Borrower\.analog\(\) uses the name 'gain'"
                                 r".*instparams"):
            _make('Borrower', analog_b, params=(('g', 1e-3),))

    def test_bare_parameter_names_still_resolve(self):
        """The convention 63 `# noqa: F821` sites depend on.  Breaking
        the leak must not break this."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                    # noqa: F821

        el = _make('BareName', analog)('a', 'b', g=4e-3)
        el.update_iparv()
        assert np.allclose(el.i(np.array([1.0, 0.0])), [4e-3, -4e-3])

    def test_ordinary_module_globals_are_still_visible(self):
        """The private namespace is a COPY of the defining module's, so
        helpers, constants and imports resolve exactly as before."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, SCALE * g * b.V)            # noqa: F821

        el = _make('UsesGlobal', analog)('a', 'b', g=1e-3)
        el.update_iparv()
        assert np.allclose(el.i(np.array([1.0, 0.0])), [3e-3, -3e-3])

    def test_analog_can_read_its_own_parameters_out_of_globals(self):
        """Models ask the compiler for their parameter symbols by reading
        the injection back out of `globals()`; the private namespace has
        to answer the same way."""
        seen = {}

        def analog(plus, minus):
            seen['g'] = globals()['g']
            b = Branch(plus, minus)
            return Contribution(b.I, seen['g'] * b.V)

        _make('ReadsGlobals', analog)
        assert seen['g'] == sympy.Symbol('g')


class TestExplain(object):

    @staticmethod
    def _states_model():
        def analog(plus, minus):
            mid = Node('mid')
            b_in = Branch(plus, mid)
            b_out = Branch(mid, minus)
            return (Contribution(b_in.V, idt(b_in.V) + g * b_out.V),
                    Contribution(b_out.I, g * b_out.V))          # noqa: F821
        return _make('Layout', analog)

    def test_x_layout_puts_states_between_internal_nodes_and_currents(self):
        """The motivating failure: the x-vector order was reconstructed
        from `el.nodes` and `el.branches`, the state node was missed, and
        the result was `IndexError: index 5 is out of bounds for axis 0
        with size 5`."""
        el = self._states_model()('a', 'b')
        el.update_iparv()
        kinds = [k for _i, _n, k in x_layout(el)]
        assert kinds[:2] == ['terminal', 'terminal']
        assert kinds[2] == 'internal node'
        assert kinds[3].startswith('state')
        assert kinds[4].startswith('branch current')
        assert len(kinds) == el.n
        ## and the layout is the one i() actually takes
        assert len(np.asarray(el.i(np.zeros(el.n)), float)) == el.n

    def test_explain_reports_the_layout_and_the_compilation_path(self):
        el = self._states_model()('a', 'b')
        el.update_iparv()
        text = explain(el, source=False, symbolic=False)
        assert 'Layout: 2 terminals, 5 unknowns' in text
        assert 'x[2]  mid' in text
        assert 'flat lambdify expressions' in text
        assert '1 state' in text

        def chained(plus, minus):
            b = Branch(plus, minus)
            i = var(g * b.V, 'i')                                # noqa: F821
            return Contribution(b.I, i)
        el2 = _make('Chained', chained)('a', 'b')
        el2.update_iparv()
        assert 'let-chain' in explain(el2, source=False, symbolic=False)

    def test_explain_shows_the_symbolic_vectors_and_the_generated_source(self):
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                    # noqa: F821

        el = _make('Src', analog)('a', 'b', g=1e-3)
        el.update_iparv()
        text = explain(el)
        assert re.search(r'i\[0\] = .*g', text)
        assert re.search(r'G\[0,0\] = g', text)
        assert 'generated source for i(x):' in text
        assert 'def _lambdifygenerated' in text

    def test_explain_flags_a_parameter_expression_that_is_unresolved(self):
        """The failure it exists to catch: a hierarchical expression
        (`r='2*rbase'`) that nothing has resolved leaves the element
        running on the DEFAULT, and every other reading of the element
        looks normal.  A plain value resolves immediately and must NOT be
        flagged, or the flag is noise."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, b.V / r)                    # noqa: F821

        cls = _make('Stale', analog, params=(('r', 1e3),))
        text = explain(cls('a', 'b', r='2*rbase'), source=False,
                       symbolic=False)
        assert "r            = 1000.0" in text
        assert "ipar says '2*rbase'" in text

        plain = cls('a', 'b', r=7.0)
        plain.update_iparv()
        assert 'ipar says' not in explain(plain, source=False, symbolic=False)

    def test_explain_refuses_something_that_is_not_a_compiled_element(self):
        from pycircuit.circuit.elements import R
        with pytest.raises(TypeError, match='not a compiled HDL element'):
            explain(R)


class TestCheckJacobians(object):

    @staticmethod
    def _diode():
        el = eh.DiodeHdl('p', 'n', IS=1e-13)
        el.update_iparv()
        return el

    def test_it_passes_where_the_derivative_is_right(self):
        res = check_jacobians(self._diode(), [0.42, 0.0])
        assert res.ok, repr(res)
        assert 'ok' in repr(res)
        assert res.results['G']['err'] < 1e-4 * res.results['G']['scale']

    def test_it_fails_on_a_G_that_does_not_match_i(self):
        """The test that makes the instrument evidence: a *deliberately*
        halved `G`.  A checker that cannot fail on this one is
        decoration."""
        class HalfG(eh.DiodeHdl):
            def G(self, x, epar=defaultepar, params_tree=None):
                return 0.5 * np.asarray(super().G(x, epar), float)

        el = HalfG('p', 'n', IS=1e-13)
        el.update_iparv()
        res = check_jacobians(el, [0.42, 0.0])
        assert not res.ok
        assert 'FAILED' in repr(res)
        assert res.results['G']['worst'] == (0, 0)
        assert 'FAILS rtol' in repr(res)
        ## and the pairing is real: C against q is untouched and passes
        assert res.results['C']['ok']

    def test_it_fails_on_a_C_that_does_not_match_q(self):
        """A charge-carrying element, halved `C`.  Written from a local
        model rather than `DiodeSpiceHdl` because an element with a
        `Collapse` retargets its own `__class__` to a compiled variant at
        construction, which discards a subclass's overrides."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, IS * (sympy.exp(b.V / 0.026) - 1)
                                + ddt(cj * b.V))                 # noqa: F821

        base = _make('Charged', analog,
                     params=(('IS', 1e-14), ('cj', 1e-12)))

        class HalfC(base):
            def C(self, x, epar=defaultepar, params_tree=None):
                return 0.5 * np.asarray(super().C(x, epar), float)

        el = HalfC('p', 'n')
        el.update_iparv()
        res = check_jacobians(el, [0.5, 0.0])
        assert not res.ok
        assert res.results['G']['ok']
        assert not res.results['C']['ok']

    def test_it_reports_a_non_finite_jacobian_rather_than_raising(self):
        """A NaN Jacobian is the failure the documentation warns about
        and had no instrument for.  It must come back as a REPORT: the
        model still evaluates, so nothing else in the run complains."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g / b.V)                    # noqa: F821

        el = _make('Singular', analog)('a', 'b')
        el.update_iparv()
        res = check_jacobians(el, [0.0, 0.0])
        assert not res.ok
        assert any(what == 'G' for what, _idx in res.nonfinite)
        assert 'NON-FINITE' in repr(res)
        ## the same element away from the singularity is clean
        assert check_jacobians(el, [1.0, 0.0]).ok

    def test_a_wrong_length_x_names_the_layout(self):
        el = self._diode()
        with pytest.raises(ValueError,
                           match=r'(?s)DiodeHdl takes an x of length 2, got '
                                 r'\[3\].*x\[0\]'):
            check_jacobians(el, [0.1, 0.0, 0.0])

    def test_the_element_method_is_the_same_check(self):
        el = self._diode()
        assert el.check_jacobians([0.3, 0.0]).ok
        assert 'DiodeHdl' in el.explain(source=False, symbolic=False)


class TestThePcnrEligibilityClaim(object):
    """`limexp`'s docstring tells a model author when to prefer `expl`.
    The condition it states is PCNR eligibility, and these pin the two
    halves of that claim so the docstring stays a measurement."""

    @pytest.mark.parametrize('card, qualifies', [
        (dict(), True),
        (dict(cjo=0.0, tt=0.0, rs=0.0), True),
        (dict(cjo=1e-12, tt=1e-9, rs=2.0), False),
    ])
    def test_the_full_spice_diode_qualifies_only_without_rs(self, card,
                                                            qualifies):
        """What decides it is `rs`, and NOT for the stated reason.

        This test used to read `never_qualifies_for_pcnr`, because
        `var()` disqualified every chained model; roadmap 10.2 removed
        that, and what is left is the series resistance.  Roadmap 10.2
        predicted the outcome and gave the wrong reason for it -- it says
        the rs diode's "current spans two branch voltages and the
        free-symbol test still refuses it".  Measured, EACH of the two
        contributions is a function of its own branch voltage alone
        (that is what the internal node buys); what refuses the device is
        that the series resistor's current is not exponential, and the
        rule is every current, not some current.  With `rs = 0` the
        collapse removes that contribution entirely and the diode
        qualifies -- charge, chain and all.
        """
        el = eh.DiodeSpiceHdl('a', 'b', **card)
        el.update_iparv()
        assert hasattr(el, 'pcnr_i') is qualifies
        assert el._hdl_info['chained']

    def test_charge_alone_does_not_disqualify_a_device(self):
        """Recorded because the opposite is written down in more than one
        place: `pcnr.augmented_system` used to refuse a charge-carrying
        participant and no longer does."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, IS * (sympy.exp(b.V / 0.026) - 1)
                                + ddt(cj * b.V))                 # noqa: F821

        el = _make('ChargedJunction', analog,
                   params=(('IS', 1e-14), ('cj', 1e-12)))('a', 'b')
        el.update_iparv()
        assert hasattr(el, 'pcnr_i')
        assert 'PCNR: 1 junction' in explain(el, source=False, symbolic=False)

    def test_a_var_chain_qualifies_for_pcnr(self):
        """The other half, and it reversed with roadmap 10.2.

        The detector walks the let-chain instead of reading the assembled
        expression, so a junction behind a `var()` recovers the same
        `IS`/`VT` as its flat spelling.  `test_chained_first_class.py`
        pins the numbers; this one pins the claim the `limexp` docstring
        makes.
        """
        def analog(plus, minus):
            b = Branch(plus, minus)
            i = var(IS * (sympy.exp(b.V / 0.026) - 1), 'i')      # noqa: F821
            return Contribution(b.I, i)

        el = _make('ChainedJunction', analog,
                   params=(('IS', 1e-14),))('a', 'b')
        el.update_iparv()
        assert el._hdl_info['chained']
        assert hasattr(el, 'pcnr_i')
        assert 'PCNR: 1 junction' in explain(el, source=False,
                                             symbolic=False)
