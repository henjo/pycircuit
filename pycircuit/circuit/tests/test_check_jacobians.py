"""`check_jacobians`'s third verdict, and `$limit` parameters written
with `var()` -- roadmap `hdl_roadmap_260824.md` sections 12.2 and 12.4.

**12.2.** A finite difference is not an oracle.  Three separate
mechanisms make it report a large discrepancy on a model that is right,
and each was hit by a real model in this tree before it was understood:

1. a **kink** -- a window seam or a clamp corner.  Central differencing
   returns the average of the two one-sided slopes while the Jacobian
   returns one of them, and no ``h`` helps because a jump has no scale;
2. **roundoff** -- an entry whose own signal is below the representable
   step of the value it differentiates, so the difference comes back
   exactly ``0.0``;
3. **truncation** -- a stiff card where the ``h^2`` term is not small.

Each of the three is built here DELIBERATELY, in an element of a dozen
lines, so the instrument is tested against a mechanism rather than
against whichever model happened to expose it.

**And the half that matters is the second one.**  A change that marked
everything unresolved would satisfy every "is it ok?" assertion in this
tree and destroy the instrument, so every element below is also run with
a deliberately corrupted derivative at the SAME point, and must still be
caught.  `TestTheCorruptedJacobianIsStillCaught` is the reason this file
is evidence.

**12.4.** `limit_pnj`'s parameter expressions used to be lambdified over
the parameter namespace alone, so an `_IntermediateSymbol` from `var()`
was not in scope and the model died with a bare ``NameError`` on the
first Newton iteration.  They are now resolved against the let-chain.
"""
import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import hdl
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution, ddt,
                                   var, maxc, limexp, vt, limit_pnj,
                                   limit_fet, check_jacobians,
                                   JAC_PASS, JAC_UNRESOLVED, JAC_FAIL)
from pycircuit.utilities.param import Parameter


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    """Every element here is evaluated numerically.  Set centrally, so a
    test added later cannot forget it."""
    old = cm.default_toolkit
    cm.default_toolkit = numeric
    yield
    cm.default_toolkit = old


def _make(name, analog, params=(('g', 1e-3),)):
    """Compile a class from a plain function, so an element built to
    misbehave lives inside a test rather than at import time."""
    return hdl.BehaviouralMeta(
        name, (Behavioural,),
        dict(instparams=[Parameter(n, default=d) for n, d in params],
             analog=staticmethod(analog)))


def _inst(cls, **kw):
    el = cls('a', 'b', **kw)
    el.update_iparv()
    return el


## ----------------------------------------------------------------------
## The three elements, one per mechanism.

def _kinked():
    """``i = g * max(v, 0)``.  C0 everywhere, C1 everywhere but ``v = 0``.

    This is the memristor's clamp corner and its window seam reduced to
    one line: at the corner the forward slope is ``g``, the backward one
    is ``0``, a central difference returns ``g/2``, and the analytic
    Jacobian returns one of the arms.
    """
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, g * maxc(b.V, 0.0))         # noqa: F821
    return _inst(_make('Kinked', analog))


def _frozen(A=1e11, c=1e-25):
    """A charge computed as a CANCELLATION, EKV's deep-cutoff shape.

    ``big`` is held at ``1e11`` and the charge is ``c * (big - 1e11)``,
    so ``dq/dv`` is exactly ``c`` while ``q``'s representable step is
    ``eps * 1e11 = 1.5e-5`` times ``c``.  A 1e-7 step therefore does not
    move ``q`` AT ALL and the difference is exactly zero -- which is what
    EKV's normalised charges do at ``vgs = -3 V``.

    The cancellation has to live behind a `var()` or sympy simplifies it
    away, which is also why the real one is invisible from the outputs:
    the magnitude that sets the quantum is INTERNAL.
    """
    def analog(plus, minus):
        b = Branch(plus, minus)
        big = var(A + b.V, 'big')
        return Contribution(b.I, g * b.V                     # noqa: F821
                            + ddt(c * (big - A)))
    return _inst(_make('Frozen', analog))


def _stiff(K=1e5):
    """``i = g * exp(K*v)`` with ``K = 1e5``: a stiff card.

    The memristor at ``ron = 1, roff = 1e9`` in one line.  At ``v = 0``
    the third derivative is ``g*K^3 = 1e12``, so a 1e-7 central
    difference carries ``h^2/6 * g*K^3`` of truncation -- 1.7e-3 against
    a band of 1e-2 * 1e-7 -- and the analytic Jacobian gets blamed.
    """
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, g * sympy.exp(K * b.V))     # noqa: F821
    return _inst(_make('Stiff', analog))


def _diode():
    """The control: an ordinary diode at an ordinary bias, where every
    noise floor must land well below the tolerance band."""
    el = eh.DiodeHdl('p', 'n', IS=1e-13)
    el.update_iparv()
    return el


def _scaled(el, which, factor):
    """The same element with ``G`` (or ``C``) deliberately wrong by
    ``factor``.

    Built by subclassing the element's OWN class, so the corrupted twin
    differs from the original in exactly one multiplication.
    """
    base = type(el)

    def bad(self, x, epar=defaultepar, params_tree=None, _b=base, _w=which):
        return factor * np.asarray(getattr(_b, _w)(self, x, epar), float)

    cls = type('Wrong' + which + base.__name__, (base,), {which: bad})
    out = cls('a', 'b')
    out.iparv = el.iparv
    return out


## ----------------------------------------------------------------------
## 1.  Each mechanism is reported as UNRESOLVED, and named.

class TestTheThreeMechanisms(object):

    def test_a_kink_is_unresolved_and_only_at_the_kink(self):
        """A jump in the slope, detected from the VALUE alone.

        The discriminator is that the one-sided disagreement
        ``|f(x+h) - 2f(x) + f(x-h)| / h`` STAYS PUT as ``h`` shrinks,
        where a smooth function's halves.  It never looks at the
        analytic Jacobian, which is what makes it a statement about the
        model rather than an excuse for the comparison.
        """
        el = _kinked()
        res = check_jacobians(el, [0.0, 0.0])
        assert res.ok, res
        assert res.verdict == 'UNRESOLVED'
        assert not res.resolved
        assert {u.reason for u in res.unresolved} == {'kink'}
        ## the size is the mechanism: the difference returns the average
        ## of the two arms, so it is out by half the jump.
        u = res.unresolved[0]
        assert_close_ratio(u.err, 0.5 * 1e-3, 1e-6)
        ## and one step either side of the corner it is an ordinary
        ## linear element again, fully resolved.
        for v in (0.5, -0.5):
            clean = check_jacobians(el, [v, 0.0])
            assert clean.ok and clean.resolved, clean

    def test_a_signal_below_the_values_quantum_is_unresolved(self):
        """The difference comes back EXACTLY zero, and the floor for it
        cannot be computed -- it has to be measured.

        The recorded proposal for this was ``max(atol, eps*|q|/h)``.
        Measured on the model that motivated it, that formula is TEN
        DECADES too small: EKV's ``q`` is 4.2e-27 there and its true quantum
        is 4.1e-31, because the cancellation happens on an INTERNAL
        magnitude of 1.9e-15 that no output shows.  So `check_jacobians`
        widens the step until the value clears its own quantisation and
        reports what it then measures.
        """
        el = _frozen()
        res = check_jacobians(el, [0.0, 0.0])
        assert res.ok, res
        assert not res.resolved
        assert {u.reason for u in res.unresolved} == {'roundoff'}
        ## The plain difference at the default step is literally zero...
        q = np.asarray(el.q([0.0, 0.0], defaultepar), float)
        qp = np.asarray(el.q([1e-7, 0.0], defaultepar), float)
        assert np.all(q == qp)
        ## ...and the widened one is a real number that corroborates the
        ## analytic entry to better than 2%, which "FAILS, fd = 0.0"
        ## never did.
        u = [e for e in res.unresolved if (e.row, e.col) == (0, 0)][0]
        assert u.fd != 0.0
        assert_close_ratio(u.fd, u.ana, 0.02)

    def test_a_stiff_card_is_unresolved_by_truncation(self):
        """``h^2/6 * f'''`` measured per entry by RICHARDSON -- the same
        column differenced at ``h`` and at ``2h``, whose disagreement is
        three times the ``h^2`` term -- rather than bounded a priori.

        There is no single ``h`` that serves this and an ordinary card:
        measured on the memristor, 1e-9 fixes the stiff one and breaks
        the default at ``x = 20``.
        """
        el = _stiff()
        res = check_jacobians(el, [0.0, 0.0])
        assert res.ok, res
        assert not res.resolved
        assert {u.reason for u in res.unresolved} == {'truncation'}
        ## the size is h^2/6 * g * K^3, within a factor of two
        assert_close_ratio(res.unresolved[0].err,
                           1e-14 / 6.0 * 1e-3 * 1e15, 0.5)

    def test_an_ordinary_model_is_fully_resolved(self):
        """The control, and it is not decoration: every floor above must
        land BELOW the tolerance band on a well-conditioned element, or
        the instrument has been widened rather than taught.

        Measured on a diode at 0.42 V: band 4.4e-12, roundoff floor
        2.5e-15, truncation floor 7.1e-15, no kink.
        """
        res = check_jacobians(_diode(), [0.42, 0.0])
        assert res.ok and res.resolved, res
        assert res.verdict == 'ok'
        assert res.unresolved == []
        g = res.results['G']
        assert np.all(g['status'] == JAC_PASS)
        assert np.max(g['roundoff']) < 0.1 * g['atol']
        assert np.max(g['truncation']) < 0.1 * g['atol']
        assert not np.any(g['kink'])


## ----------------------------------------------------------------------
## 2.  THE HALF THAT MATTERS.

class TestTheCorruptedJacobianIsStillCaught(object):
    """Every element of section 1, at the SAME point, with its analytic
    Jacobian deliberately wrong.

    An UNRESOLVED verdict is a statement that the finite difference is
    noisier than the discrepancy.  A factor of two or three is not a
    discrepancy of that size, so all of these must still FAIL -- and if
    any of them stops failing, the third verdict has become a way of not
    checking.

    Measured on the memristor's own sweep with `G` halved throughout:
    **156 of 162 points still caught it**, and the six that did not are
    one point of one card where a 1e-7 difference carries no information
    about the derivative at all.
    """

    ## `factor` is 0.5 everywhere but the kink, and THAT EXCEPTION IS
    ## ARITHMETIC, not a weakness of the kink verdict.  At a corner
    ## between a slope of `g` and a slope of `0` the central difference
    ## sits at `g/2`, so a HALVED analytic Jacobian lands exactly on it:
    ## `assert 0.5*g != g/2` is not a test anyone can pass.  The same
    ## coincidence appears on the memristor's stiff card at `x = 1`,
    ## where the forward arm is likewise flat.  A factor of three is
    ## outside the jump and is caught.
    @pytest.mark.parametrize('name,el,x,which,factor', [
        ('kink', 'kinked', [0.0, 0.0], 'G', 3.0),
        ('frozen', 'frozen', [0.0, 0.0], 'C', 0.5),
        ('stiff', 'stiff', [0.0, 0.0], 'G', 0.5),
        ('diode', 'diode', [0.42, 0.0], 'G', 0.5),
    ])
    def test_a_wrong_jacobian_is_caught_at_the_awkward_point(
            self, name, el, x, which, factor):
        good = {'kinked': _kinked, 'frozen': _frozen, 'stiff': _stiff,
                'diode': _diode}[el]()
        assert check_jacobians(good, x).ok            # the control
        bad = _scaled(good, which, factor)
        res = check_jacobians(bad, x)
        assert not res.ok, res
        assert res.verdict == 'FAILED'
        assert any(f.which == which for f in res.failures), res
        assert 'FAILED' in repr(res)

    def test_the_kinks_blind_spot_is_the_jump_and_is_measured(self):
        """Named so it expires if the boundary moves.

        At the `maxc(v, 0)` corner the analytic Jacobian is one arm and
        the difference is their average, so an error of half the jump is
        exactly what the point already has and cannot be separated from
        it.  This pins that the blind spot is that size and no larger --
        which is the honest content of the kink verdict.
        """
        el = _kinked()
        ## the difference reports 0.5*g (the average of `g` and `0`) and
        ## the floor is 0.75 of the jump, so the analytic entry may be
        ## anywhere in [-0.25*g, 1.25*g] without being separable from
        ## the corner.  Asserted as a boundary, both sides.
        assert check_jacobians(_scaled(el, 'G', 1.2), [0.0, 0.0]).ok
        assert not check_jacobians(_scaled(el, 'G', 1.5), [0.0, 0.0]).ok
        assert not check_jacobians(_scaled(el, 'G', -0.5), [0.0, 0.0]).ok

    def test_the_ladder_does_not_excuse_a_value_that_never_moves(self):
        """The widening probe's discriminator, tested against its own
        failure mode.

        A charge that is genuinely INDEPENDENT of the branch voltage
        stays frozen at every step up to the cap, and a non-zero ``C``
        against it is a real defect, not a resolution limit.  Without
        this the ladder would be a licence to pass anything whose
        difference came back zero.
        """
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V)                # noqa: F821

        el = _inst(_make('Flat', analog))

        class Invented(type(el)):
            def C(self, x, epar=defaultepar, params_tree=None):
                out = np.zeros((2, 2))
                out[0, 0] = 1e-13                # a capacitance that is not
                return out                       # in the charge at all

        liar = Invented('a', 'b')
        liar.update_iparv()
        res = check_jacobians(liar, [0.3, 0.0])
        assert not res.ok, res
        assert [(f.which, f.row, f.col) for f in res.failures] == [('C', 0, 0)]

    def test_a_wrong_entry_larger_than_the_kink_still_fails(self):
        """A kink excuses an error up to the size of the jump and no
        further, and that limit is asserted rather than described.

        This is the honest boundary of the kink verdict: at a corner the
        derivative is genuinely two-valued, so an error SMALLER than the
        jump cannot be separated from the ambiguity.  An error larger
        than it can be, and is.
        """
        el = _kinked()
        ## The difference sits at 0.5*g and the floor is 0.75 of the
        ## jump, so a factor of 0.6 is 0.1*g away and inside; a factor
        ## of 3 is 2.5*g away and well outside.
        assert check_jacobians(_scaled(el, 'G', 0.6), [0.0, 0.0]).ok
        assert not check_jacobians(_scaled(el, 'G', 3.0), [0.0, 0.0]).ok


## ----------------------------------------------------------------------
## 3.  The verdict surface itself.

class TestTheVerdictSurface(object):

    def test_ok_and_resolved_are_different_claims(self):
        res = check_jacobians(_kinked(), [0.0, 0.0])
        assert res.ok is True
        assert res.resolved is False
        assert res.verdict == 'UNRESOLVED'
        ## `ok` is the one existing tests assert, so an unresolved entry
        ## must not fail them; `resolved` is the stronger statement and
        ## has to be asked for.
        assert res.failures == []
        assert len(res.unresolved) > 0

    def test_the_repr_says_which_entries_and_why(self):
        text = repr(check_jacobians(_frozen(), [0.0, 0.0]))
        assert 'not resolvable here' in text
        assert 'UNRESOLVED' in text
        assert 'roundoff' in text
        assert 'noise floor' in text
        ## and it names the row and column by their layout labels
        assert '(plus / plus)' in text

    def test_entries_carry_the_numbers_not_just_a_flag(self):
        res = check_jacobians(_stiff(), [0.0, 0.0])
        u = res.unresolved[0]
        assert u.which == 'G'
        assert u.err == pytest.approx(abs(u.ana - u.fd))
        assert u.err <= u.floor
        assert u.reason == 'truncation'
        ## the per-entry status matrix is the programmatic form
        st = res.results['G']['status']
        assert set(np.unique(st)) <= {JAC_PASS, JAC_UNRESOLVED, JAC_FAIL}
        assert np.count_nonzero(st == JAC_UNRESOLVED) == len(res.unresolved)

    def test_a_non_finite_jacobian_is_not_quietly_unresolved(self):
        """NOT COMPARABLE outranks UNRESOLVED.  A NaN is not a
        resolution limit and must not be absorbed by one."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g / b.V)                # noqa: F821

        el = _inst(_make('Singular', analog))
        res = check_jacobians(el, [0.0, 0.0])
        assert not res.ok
        assert res.verdict == 'NOT COMPARABLE'
        assert any(what == 'G' for what, _ in res.nonfinite)

    def test_an_explicit_atol_still_wins_where_it_is_wider(self):
        """The `atol` argument is untouched: it is still ONE band for
        the whole matrix, and the per-entry floors sit under it."""
        el = _stiff()
        assert check_jacobians(el, [0.0, 0.0], atol=1.0).resolved


## ----------------------------------------------------------------------
## 4.  Roadmap 12.4 -- `$limit` parameters may be written with `var()`.

def _pnj_pair():
    """The same diode twice: once with the limiter's saturation current
    named by `var()`, once with the arithmetic spelled out again.

    This is the Gummel-Poon BJT's duplicate `isT`, reduced to the
    smallest model that has it.
    """
    def with_var(plus, minus):
        b = Branch(plus, minus)
        isT = var(IS * 2.0, 'isT')                           # noqa: F821
        v = limit_pnj(b.V, isT, vt())
        return Contribution(b.I, isT * (limexp(v / vt()) - 1))

    def spelled_twice(plus, minus):
        b = Branch(plus, minus)
        isT = var(IS * 2.0, 'isT')                           # noqa: F821
        v = limit_pnj(b.V, IS * 2.0, vt())                   # noqa: F821
        return Contribution(b.I, isT * (limexp(v / vt()) - 1))

    p = (('IS', 1e-14),)
    return (_inst(_make('PnjVar', with_var, p)),
            _inst(_make('PnjTwice', spelled_twice, p)))


class TestLimitParametersMayUseVar(object):

    def test_both_spellings_give_the_same_limiter_parameters(self):
        a, b = _pnj_pair()
        args_a = hdl._args_of(a, defaultepar)
        args_b = hdl._args_of(b, defaultepar)
        pa = [float(f(*args_a)) for f in a._hdl_info['limit_spec'][0][3]]
        pb = [float(f(*args_b)) for f in b._hdl_info['limit_spec'][0][3]]
        assert pa == pb
        assert pa[0] == 2e-14                       # IS * 2, through the chain

    def test_both_spellings_limit_identically(self):
        """The parameters agreeing is not the claim; the LIMITER
        behaving identically is.  Swept over both branches of `pnjlim`
        and over its no-op band."""
        a, b = _pnj_pair()
        rng = np.random.default_rng(20260825)
        worst = 0.0
        for _ in range(200):
            x = np.array([rng.uniform(-5.0, 1.2), 0.0])
            x0 = np.array([rng.uniform(-5.0, 1.2), 0.0])
            la = a.limit(x, x0, defaultepar)
            lb = b.limit(x, x0, defaultepar)
            worst = max(worst, float(np.max(np.abs(la - lb))))
        assert worst == 0.0
        ## and the pair really is exercising the limiter, not agreeing
        ## because neither ever bites.
        bites = a.limit(np.array([2.0, 0.0]), np.array([0.6, 0.0]),
                        defaultepar)
        assert bites[0] < 2.0

    def test_it_works_for_limit_fet_too(self):
        """The resolution is in `limit_spec`'s compile, not in
        `limit_pnj`, so every kind gets it."""
        def analog(d, s):
            b = Branch(d, s)
            vth = var(VTO + 0.25, 'vth')                     # noqa: F821
            vgs = limit_fet(b.V, vth)
            return Contribution(b.I, g * maxc(vgs - VTO, 0.0))  # noqa: F821

        el = _inst(_make('FetVar', analog, (('g', 1e-3), ('VTO', 0.5))))
        (_ra, _rb), kind, _move, pars = el._hdl_info['limit_spec'][0]
        assert kind == 'fet'
        assert float(pars[0](*hdl._args_of(el, defaultepar))) == 0.75

    def test_a_var_that_reads_the_solution_is_evaluated_at_the_last_iterate(self):
        """A limiter parameter may read the solution, and is then taken
        at the LAST ACCEPTED iterate -- SPICE's `von` semantics.

        This used to assert a compile-time refusal ("evaluated BEFORE the
        device, so it cannot depend on the solution").  The order was
        right and the conclusion did not follow: a limiter is handed
        `x0` precisely so it can measure against it, and a parameter
        evaluated there is as well-defined as `vold` is.  Inverted on
        2026-08-25.

        The check is DISCRIMINATING by construction: `IS` falls thirteen
        decades per volt, so `pnjlim`'s critical voltage is ~0.73 V at
        the last iterate (0 V) and ~1.4 V at the proposed one (0.9 V).
        Limiting fires only if the parameter was taken at `x0`.
        """
        from pycircuit.circuit.hdl import expl

        def analog(plus, minus):
            b = Branch(plus, minus)
            isv = var(IS * expl(-0.7755 * b.V / vt()), 'isv')  # noqa: F821
            v = limit_pnj(b.V, isv, vt())
            return Contribution(b.I, isv * (limexp(v / vt()) - 1))

        el = _inst(_make('SolPnj', analog, (('IS', 1e-14),)))
        fn = el._hdl_info['limit_spec'][0][3][0]
        assert fn._wants_x is True
        out = el.limit(np.array([0.9, 0.0]), np.array([0.0, 0.0]))
        assert out[0] < 0.9, out            # it bit: parameter taken at x0
        ## and the constant parameter beside it is still a plain callable
        assert el._hdl_info['limit_spec'][0][3][1]._wants_x is False

    def test_the_chain_is_read_not_inlined(self):
        """`var()` exists to stop an expression doubling at every level,
        and resolving it into a limiter parameter must not undo that.

        A chain 40 deep, each level mentioning the previous one twice,
        is 2**40 tree occurrences if inlined.  It compiles here because
        `_chain_compile` walks the DAG instead.
        """
        def analog(plus, minus):
            b = Branch(plus, minus)
            s = var(IS * 1.0, 'lvl0')                        # noqa: F821
            for k in range(40):
                s = var(s * s / (s + 1e-30), 'lvl%d' % (k + 1))
            v = limit_pnj(b.V, s, vt())
            return Contribution(b.I, s * (limexp(v / vt()) - 1))

        el = _inst(_make('DeepPnj', analog, (('IS', 1e-14),)))
        val = float(el._hdl_info['limit_spec'][0][3][0](
            *hdl._args_of(el, defaultepar)))
        assert np.isfinite(val) and val > 0.0


def assert_close_ratio(got, want, rtol):
    assert abs(got - want) <= rtol * abs(want), (got, want)
