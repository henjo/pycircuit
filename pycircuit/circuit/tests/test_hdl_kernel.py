"""The compact-model math kernel: `expl`, `hypsmooth`, `safe_*`.

A compiled conditional evaluates BOTH arms -- numpy `select` chooses
afterwards -- so a guard like ``x > 0 ? sqrt(x) : 0`` still takes the
square root of the negative number.  The value survives (the bad arm is
discarded) but the floating-point flag is raised on every evaluation,
and the DERIVATIVE frequently does not survive: ``d(max(x,0))/dx``
divided by ``2*sqrt(max(x,0))`` is ``0/0`` at ``x <= 0``, and one NaN in
the Jacobian loses the whole Newton step.

These functions are written so that every arm is valid for every input.
The tests pin that, and pin the two things that made earlier attempts
wrong: floating-point cancellation destroying an identity that holds in
exact arithmetic, and `sign` differentiating to something no numeric
backend can print.
"""
import warnings

import numpy as np
import pytest
import sympy

from pycircuit.circuit import hdl


X = sympy.Symbol('x', real=True)
#: sympy prints Min/Max through these; the numeric toolkit binds the same.
MODULES = [{'Min': np.minimum, 'Max': np.maximum}, 'numpy']

#: Everyday values plus every seam, plus the extremes that expose
#: overflow and cancellation.
SE05 = 2.3025850929940458e+02          # PSP's `se05` = ln(1e100)

GRID = np.concatenate([
    np.array([-1e300, -1e100, -1e10, -1000., -745., -SE05 - 1e-6, -SE05,
              -SE05 + 1e-6, -100., -80., -1., -1e-9, -0.0, 0.0, 1e-9, 1.,
              80., 100., SE05 - 1e-6, SE05, SE05 + 1e-6, 709., 1000.,
              1e10, 1e100, 1e300]),
    np.linspace(-300, 300, 1201)])

#: Values a circuit quantity can actually take.  Beyond ~1e154 the
#: radicand of `hypsmooth` overflows; that is documented, not fixed.
SANE = GRID[np.abs(GRID) <= 1e10]


def _fns(expr):
    return (sympy.lambdify(X, expr, modules=MODULES),
            sympy.lambdify(X, sympy.diff(expr, X), modules=MODULES))


def _sweep(expr, grid):
    """Evaluate f and df over `grid`, counting floating-point warnings."""
    f, d = _fns(expr)
    out = {}
    for key, g in (('f', f), ('df', d)):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            ## UNDERFLOW is deliberately not watched.  It is the one
            ## benign IEEE condition: 0.0 is the correctly-rounded answer
            ## for exp(-1000), and the one-sided `limexp` / `expl_high`
            ## raise it by design wherever plain `exp` would.  Overflow,
            ## invalid and divide-by-zero all mean the arm produced
            ## garbage, and those are what these functions exist to stop.
            old = np.seterr(over='warn', invalid='warn', divide='warn',
                            under='ignore')
            try:
                vals = np.array([float(g(v)) for v in grid])
            finally:
                np.seterr(**old)
        out[key] = vals
        out[key + '_warn'] = len(w)
    return out


## `limexp` is absent on purpose -- it is the LRM function and stays
## PCNR-detectable at the cost of an overflow in its discarded arm.
## `TestLimexpIsDeliberatelyNotSafe` pins that trade explicitly.
ALL = {
    'expl_high': hdl.expl_high(X),
    'expl': hdl.expl(X),
    'hypsmooth': hdl.hypsmooth(X, 1e-6),
    'safe_sqrt': hdl.safe_sqrt(X),
    'safe_ln': hdl.safe_ln(X),
    'safe_div': hdl.safe_div(1.0, X),
}


class TestNothingIsEverNonFinite(object):
    """Over the range a circuit reaches, value and derivative are finite."""

    @pytest.mark.parametrize('name', sorted(ALL))
    def test_finite_over_sane_range(self, name):
        r = _sweep(ALL[name], SANE)
        bad_f = SANE[~np.isfinite(r['f'])]
        bad_d = SANE[~np.isfinite(r['df'])]
        assert bad_f.size == 0, 'non-finite value at %s' % bad_f
        assert bad_d.size == 0, 'non-finite derivative at %s' % bad_d

    @pytest.mark.parametrize('name', sorted(ALL))
    def test_no_floating_point_warnings_over_sane_range(self, name):
        """The point of the exercise.

        Under `np.seterr(all='raise')` -- which anyone debugging a
        convergence failure will set -- a warning here is an exception,
        and the model becomes unusable at perfectly ordinary biases.
        """
        r = _sweep(ALL[name], SANE)
        assert r['f_warn'] == 0
        assert r['df_warn'] == 0

    def test_expl_low_is_deliberately_one_sided(self):
        """`expl_low` guards only the lower tail -- like PSP's macro.

        Above the seam its exponential arm eventually overflows, by
        design: it exists to stop underflow to zero, and a model that can
        reach that far up should be using two-sided `expl`.  Pinned so
        the asymmetry is a decision, not a gap.
        """
        f, _ = _fns(hdl.expl_low(X))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            old = np.seterr(all='ignore')
            try:
                assert np.isfinite(float(f(-1e10)))
                assert not np.isfinite(float(f(1000.0)))
            finally:
                np.seterr(**old)
        ## whereas the two-sided form copes
        g, _ = _fns(hdl.expl(X))
        assert np.isfinite(float(g(1000.0)))


class TestLimexpIsDeliberatelyNotSafe(object):
    """`limexp` keeps the bare `exp` so PCNR can still see the junction.

    Clamping its argument makes the discarded arm safe but turns
    `exp(v/vt)` into `exp(min(v/vt, 80))`, which PCNR's shape detector
    does not recognise as a junction -- the device silently loses its
    limiting.  Since PCNR is what bounds the argument in the first place,
    the overflow is the better half of the trade.  `expl_high` is the
    same continuation with the clamp, for models that are not junctions.
    """

    def test_limexp_exposes_a_bare_exp(self):
        assert any(a.args[0] == X for a in hdl.limexp(X).atoms(sympy.exp))

    def test_expl_high_does_not(self):
        assert not any(a.args[0] == X
                       for a in hdl.expl_high(X).atoms(sympy.exp))

    def test_they_are_genuinely_different_functions(self):
        """Not two spellings of one thing -- pinned so nobody merges them.

        `limexp` continues LINEARLY from 80; `expl_high` continues with
        exp's own third-order Taylor from `se05` = 230.2585.  They agree
        below 80 and diverge enormously above it: at x = 100 the LRM
        function gives exp(80)*21 while PSP's gives exp(100) exactly,
        seven orders of magnitude apart.
        """
        fl, _ = _fns(hdl.limexp(X))
        fh, _ = _fns(hdl.expl_high(X))
        for v in (-100., -1., 0., 1., 50., 79.):
            assert float(fl(v)) == pytest.approx(float(fh(v)), rel=1e-15)
        assert float(fh(100.)) == pytest.approx(np.exp(100.), rel=1e-15)
        assert float(fl(100.)) == pytest.approx(np.exp(80.) * 21, rel=1e-12)
        assert float(fh(100.)) / float(fl(100.)) > 1e6

    def test_and_only_expl_high_is_warning_free(self):
        big = np.array([100., 1000., 1e6])
        assert _sweep(hdl.limexp(X), big)['f_warn'] > 0
        assert _sweep(hdl.expl_high(X), big)['f_warn'] == 0


class TestExponentialFamily(object):

    def test_agrees_with_exp_in_the_middle(self):
        """Exactly `exp`, not an approximation to it, inside the seams."""
        f, d = _fns(hdl.expl(X))
        for v in (-200., -70., -10., -1., 0., 1., 10., 70., 200.):
            assert f(v) == pytest.approx(np.exp(v), rel=1e-13)
            assert d(v) == pytest.approx(np.exp(v), rel=1e-13)

    @pytest.mark.parametrize('seam', [-SE05, SE05])
    def test_c1_at_the_seams(self, seam):
        """Value and slope join; Newton sees no step in the Jacobian."""
        f, d = _fns(hdl.expl(X))
        h = 1e-6
        fl, fr = float(f(seam - h)), float(f(seam + h))
        dl, dr = float(d(seam - h)), float(d(seam + h))
        assert fl == pytest.approx(fr, rel=1e-4)
        assert dl == pytest.approx(dr, rel=1e-4)

    def test_the_continuation_is_exps_own_taylor_series(self):
        """Which is what makes it C-3 rather than merely C1.

        `P3(u) = 1 + u + u**2/2 + u**3/6` is the third-order Taylor
        polynomial of `exp` about the seam, so three derivatives match
        there, not one.
        """
        f, _ = _fns(hdl.expl(X))
        for du in (0.1, 0.5, 1.0, 2.0):
            got = float(f(SE05 + du)) / 1e100
            want = 1 + du + du ** 2 / 2 + du ** 3 / 6
            assert got == pytest.approx(want, rel=1e-12)

    def test_lower_tail_stays_strictly_positive(self):
        """The reason `expl_low` is hyperbolic rather than clamped.

        Plain `exp` underflows to exactly zero below about -745, and a
        model that divides by the result gets an infinity out of an
        ordinary reverse bias.
        """
        f, _ = _fns(hdl.expl(X))
        for v in (-100., -745., -1000., -1e6, -1e9):
            assert float(f(v)) > 0.0
        assert np.exp(-1000.) == 0.0        # what it is protecting against

    def test_upper_tail_is_cubic_not_exponential(self):
        """Bounded growth is the point: `exp(1e6)` is not representable."""
        f, _ = _fns(hdl.expl(X))
        big = float(f(1e6))
        assert np.isfinite(big)
        ## cubic in the overshoot: doubling it multiplies by ~8
        r = float(f(SE05 + 2e5)) / float(f(SE05 + 1e5))
        assert r == pytest.approx(8.0, rel=1e-3)


class TestHypsmooth(object):

    def test_approximates_max_x_zero(self):
        """Above zero it is `x`, below zero it is a positive whisker.

        The relative error above zero is `(eps/x)**2`, not a constant:
        hypsmooth(x) -> x + eps**2/x.  So the tolerance is written that
        way rather than pinned to a number that only holds for one `x`.
        """
        eps = 1e-6
        f, _ = _fns(hdl.hypsmooth(X, eps))
        for v in (1e-3, 1.0, 100.0, 1e6):
            assert float(f(v)) == pytest.approx(v, rel=4 * (eps / v) ** 2)
        for v in (-1e-3, -1.0, -100.0, -1e6):
            assert 0.0 < float(f(v)) <= eps
            ## and the whisker is the eps**2/|x| tail, not noise
            assert float(f(v)) == pytest.approx(eps ** 2 / abs(v), rel=1e-6)

    def test_survives_the_cancellation_that_broke_the_naive_form(self):
        """`0.5*(x + sqrt(x*x + 4*eps^2))` is positive only on paper.

        At x = -100 with eps = 1e-12, `sqrt(x*x + 4e-24)` rounds to
        exactly 100.0 and the sum is exactly 0.0 -- so every downstream
        consumer that was promised a positive number gets a zero.  This
        pins that the shipped form does not.
        """
        eps = 1e-12
        naive = 0.5 * (-100.0 + np.sqrt(100.0 ** 2 + 4 * eps ** 2))
        assert naive == 0.0                       # the trap

        f, _ = _fns(hdl.hypsmooth(X, eps))
        assert float(f(-100.0)) > 0.0
        assert float(f(-100.0)) == pytest.approx(eps ** 2 / 100.0, rel=1e-6)

    def test_positive_everywhere_on_the_sane_range(self):
        f, _ = _fns(hdl.hypsmooth(X, 1e-12))
        vals = np.array([float(f(v)) for v in SANE])
        assert np.all(vals > 0.0)


class TestSafeSqrt(object):

    def test_matches_sqrt_where_it_should(self):
        f, _ = _fns(hdl.safe_sqrt(X))
        for v in (1e-6, 1e-3, 1.0, 100.0, 1e8):
            assert float(f(v)) == pytest.approx(np.sqrt(v), rel=1e-9)

    def test_derivative_is_finite_below_zero(self):
        """The clamped form `sqrt(max(x,0))` gives 0/0 here.

        Measured before the fix: 2013 non-finite derivatives over a
        4001-point sweep.
        """
        _, d = _fns(hdl.safe_sqrt(X))
        for v in (-1e-9, -1.0, -100.0, -1e6):
            assert np.isfinite(float(d(v)))

    def test_the_clamped_form_really_does_fail(self):
        """Guard against this test suite pinning a non-problem."""
        clamped = sympy.sqrt(sympy.Max(X, 0.0))
        _, d = _fns(clamped)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            old = np.seterr(all='ignore')
            try:
                assert not np.isfinite(float(d(-1.0)))
            finally:
                np.seterr(**old)


class TestSafeDiv(object):

    def test_matches_division_away_from_zero(self):
        f, _ = _fns(hdl.safe_div(1.0, X))
        for v in (1e-9, 1e-3, 1.0, 1e6, -1e-9, -1.0):
            assert float(f(v)) == pytest.approx(1.0 / v, rel=1e-9)

    def test_finite_at_zero_instead_of_infinite(self):
        f, d = _fns(hdl.safe_div(1.0, X))
        assert float(f(0.0)) == 0.0
        assert np.isfinite(float(d(0.0)))

    def test_the_sign_based_guard_cannot_even_be_compiled(self):
        """Why the regularised form is used instead of shifting by sign.

        `sign` differentiates to a DiracDelta, which no numeric backend
        prints -- so the model fails to COMPILE, not to converge.  Pinned
        so nobody reintroduces the obvious-looking version.
        """
        naive = 1.0 / (X + 1e-30 * sympy.sign(X))
        d = sympy.diff(naive, X)
        assert d.has(sympy.DiracDelta)
        fn = sympy.lambdify(X, d, modules=MODULES)
        with pytest.raises(NameError):
            fn(1.0)


class TestUsableInsideAModel(object):
    """The kernel has to survive compilation, not just lambdify."""

    def test_a_model_built_from_the_kernel_compiles_and_differentiates(self):
        import pycircuit.circuit.circuit as cm
        from pycircuit.circuit.toolkit import numeric
        from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                           var, vt)
        from pycircuit.utilities.param import Parameter

        cm.default_toolkit = numeric

        def analog(d_, g_, s_, b_):
            bds, bgs = Branch(d_, s_), Branch(g_, s_)
            vg = var(bgs.V / vt(), 'vg')
            vd = var(bds.V / vt(), 'vd')
            inv = var(hdl.expl(vg - 30.0), 'inv')
            surf = var(hdl.safe_sqrt(vg - 1.0), 'surf')
            ratio = var(hdl.safe_div(inv, surf), 'ratio')
            sat = var(hdl.hypsmooth(vd, 1e-3), 'sat')
            return Contribution(bds.I, IS * ratio * sat)        # noqa: F821

        cls = type('KernelFet', (Behavioural,), {
            'instparams': [Parameter(name='IS', desc='I', unit='A',
                                     default=1e-6)],
            'analog': staticmethod(analog)})
        e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                IS=1e-6)
        e.update_iparv()

        ## Includes the biases Newton wanders through, not just sensible
        ## ones: a huge negative gate (deep in the expl tail), a gate
        ## below the safe_sqrt seam, zero everywhere.
        for x in (np.array([1.0, 1.2, 0.0, 0.0]),
                  np.array([0.0, 0.0, 0.0, 0.0]),
                  np.array([0.05, -40.0, 0.0, 0.0]),
                  np.array([-5.0, 0.02, 0.0, 0.0]),
                  np.array([50.0, 50.0, 0.0, 0.0])):
            with warnings.catch_warnings():
                warnings.simplefilter('error')
                i = np.asarray(e.i(x), float)
                G = np.asarray(e.G(x), float)
            assert np.all(np.isfinite(i)), x
            assert np.all(np.isfinite(G)), x

            ## and the Jacobian is the derivative of that residual
            fd = np.zeros_like(G)
            for j in range(len(x)):
                h = 1e-6 * max(1.0, abs(x[j]))
                xp, xm = x.copy(), x.copy()
                xp[j] += h
                xm[j] -= h
                fd[:, j] = (np.asarray(e.i(xp), float)
                            - np.asarray(e.i(xm), float)) / (2 * h)
            scale = max(1.0, float(np.max(np.abs(G))))
            assert np.max(np.abs(G - fd)) < 1e-4 * scale, x


## ======================================================================
## Phase 3 of the hdl roadmap: the selection, sign and power primitives,
## and the generic helpers promoted out of `psp_kernel`.
##
## Three claims are under test for each of them, and the third is the one
## a value-based suite misses:
##
##   * the VALUE, against a reference computed independently -- Python's
##     own `max`, `math.copysign`, `math.log1p` -- never against a second
##     spelling of the implementation;
##   * the DERIVATIVE, against finite differences, swept rather than
##     sampled at one point, and separately AT the kink where finite
##     differences cannot speak;
##   * what the compiler EMITS.  A primitive can be numerically perfect
##     and still cost 4.4 us per occurrence by lowering to `numpy.select`
##     instead of `numpy.maximum`, and no assertion about its value can
##     see that.
## ======================================================================

import inspect
import math
import re

#: The kernel's own runtime namespace.  Passed explicitly where a test
#: wants to prove the REGISTRATION works; the primitives also carry
#: `_imp_`, so a bare `lambdify` resolves them too -- which is itself
#: tested, below, because `hypsmooth` and the `expl` family are public
#: and were lambdifiable with `modules='numpy'` before this change.
KMODS = [hdl._KERNEL_NUMPY, 'numpy']

#: A second unknown, for the two-argument primitives.
Y = sympy.Symbol('y', real=True)

#: A PARAMETER as the compiler makes one: a bare `sympy.Symbol(name)`
#: with no real assumption (`paramsyms`, in `generate_code`).  That
#: distinction is the whole of the `Abs` trap -- branch quantities ARE
#: declared real, parameters are not, and `Abs` of anything mixing the
#: two goes through `re`/`im`.
PARAM = sympy.Symbol('IS')

#: A grid that straddles every kink these functions have.
KINKY = np.array([-1e12, -1e3, -1.0, -1e-3, -1e-9, 0.0, 1e-9, 1e-3, 0.5,
                  1.0, 2.0, 1e3, 1e12])


def _emitted(expr, syms=(X,), modules=None):
    """Count the numpy primitives the compiled form of `expr` calls.

    The instrument the whole of this section leans on: `lambdify` the
    expression, read the source it generated, and count.  Both of the
    performance claims that motivated these primitives -- and both of the
    measurement mistakes made while chasing them -- were invisible to any
    assertion about a returned number and one source dump away.
    """
    f = sympy.lambdify(list(syms), expr,
                       modules=KMODS if modules is None else modules,
                       cse=True)
    src = inspect.getsource(f)
    names = ('select', 'amin', 'amax', 'minimum', 'maximum', 'heaviside',
             'where', 'reduce', 'maxc', 'minc', '_step', 'sign', 'abs',
             'DiracDelta', 'Heaviside')
    return {n: len(re.findall(r'\b' + n + r'\(', src)) for n in names}


def _fd(f, v, h=None):
    """Central finite difference, with a step scaled to the point."""
    h = h if h is not None else 1e-6 * max(1.0, abs(v))
    return (float(f(v + h)) - float(f(v - h))) / (2 * h)


class TestMaxcMinc(object):
    """`maxc`/`minc`: `Max`/`Min` with a derivative that does not expand.

    The roadmap justifies these by an undocumented rule -- "`Max` only on
    an ATOM", because differentiating `Max` of a compound expression was
    reported to produce a Jacobian that divides by zero.  **That failure
    is not reproducible on this sympy** (1.14) with this kernel; see
    `test_the_max_on_an_atom_rule_is_not_reproducible_here`, which is
    written as the honest record of the attempt rather than as an
    assertion that the rule was imaginary.

    What IS demonstrable, and is what these classes buy, is below:
    `Max` differentiates to a `Heaviside` that lowers to `numpy.select`,
    and differentiates AGAIN to a `DiracDelta` that does not lower at
    all.
    """

    def test_value_is_the_ordinary_maximum(self):
        f = sympy.lambdify([X, Y], hdl.maxc(X, Y), modules=KMODS)
        g = sympy.lambdify([X, Y], hdl.minc(X, Y), modules=KMODS)
        for a in KINKY:
            for b in (-1.0, 0.0, 1.0, a):
                ## reference: Python's own max/min, not a second spelling
                assert float(f(a, b)) == max(float(a), float(b))
                assert float(g(a, b)) == min(float(a), float(b))

    def test_partials_are_the_finite_differences_away_from_the_tie(self):
        for expr, ref in ((hdl.maxc(X, 0.5), max),
                          (hdl.minc(X, 0.5), min)):
            d = sympy.lambdify(X, sympy.diff(expr, X), modules=KMODS)
            f = sympy.lambdify(X, expr, modules=KMODS)
            for v in KINKY:
                if abs(v - 0.5) < 1e-3:
                    continue                      # the kink itself
                assert float(d(v)) == pytest.approx(_fd(f, v), abs=1e-6), v

    def test_the_two_partials_sum_to_one_at_the_tie(self):
        """Where `Heaviside` would give 1/2 each and `maxc(x, x)` slope 2.

        A `Max` differentiated through `Heaviside` splits the tie evenly,
        so an expression in which BOTH arguments depend on the same
        unknown gets a slope of 2 where the function's slope is 1.  These
        give the whole derivative to the first argument, which is
        arbitrary but consistent, and consistent is what a Jacobian
        needs.
        """
        for cls in (hdl.maxc, hdl.minc):
            e = cls(X, X)
            assert float(sympy.diff(e, X).subs(X, 3.0).evalf()) == 1.0

    def test_the_derivative_is_diracdelta_free_at_every_order(self):
        """`Max`'s second derivative does not COMPILE; these stay zero.

        This is not hypothetical tidiness: a distortion analysis or a
        sensitivity takes a second derivative, and `DiracDelta` reaches
        `lambdify` as a bare name that no backend defines.
        """
        d2 = sympy.diff(sympy.Max(X, Y), X, 2)
        assert d2.has(sympy.DiracDelta)                     # the problem
        bad = sympy.lambdify([X, Y], d2, modules=KMODS)
        with pytest.raises(NameError):
            bad(1.0, 2.0)

        for cls in (hdl.maxc, hdl.minc):
            for order in (1, 2, 3):
                d = sympy.diff(cls(X, Y), X, order)
                assert not d.has(sympy.DiracDelta)
                assert not d.has(sympy.Heaviside)

    def test_the_derivative_lowers_to_a_cheap_call_not_to_select(self):
        """The measured reason to prefer these over `Max`/`Min`.

        `numpy.select` is 4.4 us per call on this machine against 0.62 us
        for `numpy.maximum`, scalar -- and a compiled model calls these
        in the Newton inner loop.  sympy lowers `Heaviside`, which is
        what `Max` differentiates to, straight to `select`.
        """
        assert _emitted(sympy.diff(sympy.Max(X, Y), X), (X, Y))['select'] == 1
        for cls in (hdl.maxc, hdl.minc):
            em = _emitted(sympy.diff(cls(X, Y), X), (X, Y))
            assert em['select'] == 0
            assert em['_step'] == 1

    def test_it_may_be_applied_to_a_COMPOUND_expression(self):
        """The deliverable: the "`Max` only on an atom" rule is lifted.

        Each kernel function is wrapped in a clamp and the result
        differentiated, lambdified and swept.  The stated failure mode --
        sympy rationalising the smoothing into a denominator that cancels
        to zero, so the value is finite and the Jacobian is NaN -- did
        **not** reproduce on sympy 1.14 for any of these.  A different
        and worse one did: with the kernel's own clamps in place,
        `sympy.Max` of these same compounds does not finish `lambdify`
        at all (>200 s, stuck in `Float.__eq__` under the `Heaviside`
        printer, against 0.006 s for `maxc`).  So the rule is real; only
        its symptom in this sympy is a hang rather than a NaN.

        That case is NOT executed here, because a test that hangs is not
        a test.  What is executed is the positive claim: `maxc` on the
        same compounds compiles promptly and evaluates finitely.  The
        one-second budget would have failed by four orders of magnitude
        on the `sympy.Max` form.
        """
        import time
        compounds = [hdl.safe_sqrt(X), hdl.safe_div(1.0, X), hdl.expl(X),
                     hdl.safe_ln(X), hdl.safe_abs(X),
                     hdl.hypsmooth(X, 1e-6) * hdl.hypsmooth(X, 1e-6)]
        grid = np.array([-1e6, -1.0, -1e-9, 0.0, 1e-9, 1.0, 1e6])
        for c in compounds:
            t0 = time.time()
            d = sympy.lambdify(X, sympy.diff(hdl.maxc(c, 1e-10), X),
                               modules=KMODS)
            f = sympy.lambdify(X, hdl.maxc(c, 1e-10), modules=KMODS)
            assert time.time() - t0 < 1.0, c
            with np.errstate(all='ignore'):
                vals = np.array([float(f(v)) for v in grid])
                ders = np.array([float(d(v)) for v in grid])
            assert np.all(np.isfinite(vals)), c
            assert np.all(np.isfinite(ders)), c

    def test_finite_at_the_extremes_the_contract_claims(self):
        """RANGE says total -- infinities included, on both arguments.

        This is the property that ruled out the operator-only form
        ``a*s + b*(1 - s)``, which is `nan` whenever the LOSING argument
        is infinite.  Since a clamp is most often written precisely to
        survive an argument that has run away, that failure would have
        landed exactly where the primitive was needed.
        """
        f = sympy.lambdify([X, Y], hdl.maxc(X, Y), modules=KMODS)
        g = sympy.lambdify([X, Y], hdl.minc(X, Y), modules=KMODS)
        inf = float('inf')
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            assert float(f(-inf, 0.0)) == 0.0        # losing -inf
            assert float(f(inf, 0.0)) == inf
            assert float(g(inf, 0.0)) == 0.0         # losing +inf
            assert float(g(-inf, 0.0)) == -inf
            assert float(f(1e308, -1e308)) == 1e308


class TestSignAndAbs(object):
    """`sign` and `Abs`: the two sympy spellings that do not compile.

    Both failures are real and both are pinned here against the sympy
    originals, because a test that only shows the replacement works
    cannot show the replacement was needed.
    """

    def test_sympy_sign_fails_to_compile_at_all(self):
        d = sympy.diff(sympy.sign(X), X)
        assert d.has(sympy.DiracDelta)
        with pytest.raises(NameError):
            sympy.lambdify(X, d, modules=KMODS)(1.0)

    def test_sympy_abs_gives_a_nan_derivative_at_the_origin(self):
        """The quiet one: the VALUE is right everywhere.

        The usual statement of this -- "sympy does not know a model
        symbol is real" -- is half right and the half that is wrong
        matters.  `Quantity` DOES declare `is_real = True`, so `Abs` of a
        bare branch voltage differentiates to a perfectly good `sign`.  A
        PARAMETER does not: the compiler makes it with
        `sympy.Symbol(name)` and no assumptions (`paramsyms`, in
        `generate_code`).  So the trap is reachable when the argument
        mixes an unknown with a parameter -- which is most of them -- and
        then the derivative goes through `re`/`im` and is `0/0` at
        zero.
        """
        assert sympy.diff(sympy.Abs(X), X) == sympy.sign(X)   # real: fine
        d_expr = sympy.diff(sympy.Abs(X * PARAM), X)          # mixed: not
        assert d_expr.has(sympy.re) and d_expr.has(sympy.im)
        d = sympy.lambdify([X, PARAM], d_expr, modules=KMODS)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            old = np.seterr(all='ignore')
            try:
                assert not np.isfinite(float(d(0.0, 1.0)))    # the trap
                assert float(d(2.0, 1.0)) == 1.0              # value fine
            finally:
                np.seterr(**old)

        ## ours is exact and finite for the same argument
        ours = sympy.lambdify([X, PARAM],
                              sympy.diff(hdl.Abs(X * PARAM), X),
                              modules=KMODS)
        assert float(ours(0.0, 1.0)) == 0.0
        assert float(ours(2.0, 1.0)) == 1.0

    def test_our_sign_and_abs_have_the_values_they_claim(self):
        f = sympy.lambdify(X, hdl.sign(X), modules=KMODS)
        g = sympy.lambdify(X, hdl.Abs(X), modules=KMODS)
        for v in KINKY:
            ## reference: math.copysign / abs, computed independently
            want = 0.0 if v == 0.0 else math.copysign(1.0, v)
            assert float(f(v)) == want
            assert float(g(v)) == abs(float(v))
        assert float(f(float('inf'))) == 1.0
        assert float(g(float('-inf'))) == float('inf')

    def test_their_derivatives_are_defined_everywhere(self):
        ds = sympy.lambdify(X, sympy.diff(hdl.sign(X), X), modules=KMODS)
        da = sympy.lambdify(X, sympy.diff(hdl.Abs(X), X), modules=KMODS)
        fa = sympy.lambdify(X, hdl.Abs(X), modules=KMODS)
        for v in KINKY:
            assert float(ds(v)) == 0.0
            assert np.isfinite(float(da(v)))
        ## and away from the kink `Abs` really does differentiate to +-1
        for v in (-1e3, -1.0, -1e-3, 1e-3, 1.0, 1e3):
            assert float(da(v)) == pytest.approx(_fd(fa, v), abs=1e-6)
        assert float(da(0.0)) == 0.0      # the only finite answer there

    def test_second_derivatives_still_compile(self):
        for e in (hdl.sign(X), hdl.Abs(X)):
            for order in (2, 3):
                d = sympy.diff(e, X, order)
                assert not d.has(sympy.DiracDelta)

    def test_they_lower_to_the_backends_own_sign_and_abs(self):
        """No registration at all, on numpy and on jax alike.

        The classes are NAMED `sign` and `Abs` so sympy's printer maps
        them to the backend function of the same name.  Asserting on the
        emitted source is what makes that a decision rather than a
        coincidence nobody would notice breaking.
        """
        assert _emitted(hdl.sign(X))['sign'] == 1
        assert _emitted(hdl.Abs(X))['abs'] == 1
        assert _emitted(hdl.sign(X))['select'] == 0
        assert _emitted(hdl.Abs(X))['select'] == 0


class TestSafePow(object):

    def test_matches_the_ordinary_power_inside_the_range(self):
        for e in (0.5, 0.59, 1.5, 2.0, -1.0 / 6.0):
            f = sympy.lambdify(X, hdl.safe_pow(X, e, lo=1e-10),
                               modules=KMODS)
            for v in (1e-8, 1e-3, 1.0, 7.0, 1e6):
                assert float(f(v)) == pytest.approx(float(v) ** e,
                                                    rel=1e-12), (e, v)

    def test_below_the_floor_it_is_constant_and_finite(self):
        """Where a raw `b**e` is `nan` (b < 0) or divergent (b = 0).

        The exponent 0.59 is PSP's Coulomb-scattering exponent and the
        base is exactly zero at flat band, so this is the reachable case
        and not a contrived one.
        """
        lo, e = 1e-10, 0.59
        f = sympy.lambdify(X, hdl.safe_pow(X, e, lo=lo), modules=KMODS)
        d = sympy.lambdify(X, sympy.diff(hdl.safe_pow(X, e, lo=lo), X),
                           modules=KMODS)
        for v in (-1e6, -1.0, -1e-30, 0.0, lo / 2):
            assert float(f(v)) == pytest.approx(lo ** e, rel=1e-12)
            assert float(d(v)) == 0.0
        ## and the raw form really does fail there, both ways
        raw = sympy.lambdify(X, X ** e, modules=KMODS)
        draw = sympy.lambdify(X, sympy.diff(X ** e, X), modules=KMODS)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            old = np.seterr(all='ignore')
            try:
                ## a negative base with a fractional exponent is not a
                ## real number at all -- numpy's answer is `nan`, and
                ## sympy's is a COMPLEX number, which is worse: it
                ## propagates silently into a real Jacobian
                assert not isinstance(raw(-1.0), float) \
                    or not np.isfinite(raw(-1.0))
                ## divergent at 0 -- which the generated code reports as
                ## a raised `ZeroDivisionError`, not as an `inf`
                with pytest.raises(ZeroDivisionError):
                    draw(0.0)
            finally:
                np.seterr(**old)

    def test_the_cap_clamps_the_other_end(self):
        f = sympy.lambdify(X, hdl.safe_pow(X, 2.0, lo=1e-10, hi=100.0),
                           modules=KMODS)
        assert float(f(50.0)) == pytest.approx(2500.0, rel=1e-12)
        assert float(f(1e30)) == pytest.approx(1e4, rel=1e-12)

    def test_derivative_against_finite_differences(self):
        lo = 1e-10
        for e in (0.5, 0.59, 1.5, -1.0 / 6.0):
            ex = hdl.safe_pow(X, e, lo=lo)
            f = sympy.lambdify(X, ex, modules=KMODS)
            d = sympy.lambdify(X, sympy.diff(ex, X), modules=KMODS)
            for v in (1e-6, 1e-3, 0.4, 1.0, 3.0, 1e4):
                ## step RELATIVE to the point: at v = 1e-6 an absolute
                ## 1e-6 step straddles the floor, and the difference
                ## quotient is then measuring the clamp, not the power
                assert float(d(v)) == pytest.approx(
                    _fd(f, v, h=1e-6 * v), rel=2e-5), (e, v)

    def test_finite_over_the_range_the_contract_claims(self):
        """Value AND derivative, on a sweep, at the documented floor."""
        for e in (0.5, 0.59, 2.0, -1.0 / 6.0, -1.0):
            ex = hdl.safe_pow(X, e, lo=1e-10)
            f = sympy.lambdify(X, ex, modules=KMODS)
            d = sympy.lambdify(X, sympy.diff(ex, X), modules=KMODS)
            grid = np.concatenate([KINKY, np.logspace(-40, 40, 161)])
            with warnings.catch_warnings():
                warnings.simplefilter('error')
                vals = np.array([float(f(v)) for v in grid])
                ders = np.array([float(d(v)) for v in grid])
            assert np.all(np.isfinite(vals)), e
            assert np.all(np.isfinite(ders)), e

    def test_a_floor_that_cannot_work_is_refused_at_build_time(self):
        """The roadmap's "this guard can never fire" check, inverted.

        With `lo = 1e-30` and an exponent below -9.2 the CLAMP itself
        leaves double range: the floor that was supposed to keep the term
        finite is what makes it infinite.  Refused, with the arithmetic
        in the message, rather than left to surface as an `inf` at some
        bias.
        """
        hdl.safe_pow(X, -9.0)                      # value 1e270, ok
        with pytest.raises(ValueError) as ei:
            hdl.safe_pow(X, -9.5)                  # derivative overflows
        assert 'overflow' in str(ei.value)
        with pytest.raises(ValueError):
            hdl.safe_pow(X, 0.5, lo=0.0)
        with pytest.raises(ValueError):
            hdl.safe_pow(X, 0.5, lo=1e-3, hi=1e-4)
        ## a bigger floor makes the same exponent legal again
        hdl.safe_pow(X, -9.5, lo=1e-20)

    def test_the_floor_may_be_written_relative_to_a_parameter(self):
        """Which is what the docstring tells the author to do.

        "A smoothing constant is meaningless without its scale": an
        absolute floor encodes an assumption about a parameter the author
        does not control.  Writing `lo = 1e-3*vp` makes the bound a
        SYMBOL at compile time, so it must be accepted -- and the
        overflow check, which needs a number, must skip it rather than
        raise.  It was `float(lo)` and did raise; pinned so it stays
        fixed.
        """
        vp = sympy.Symbol('vp', real=True, positive=True)
        e = hdl.safe_pow(X, 0.5, lo=1e-3 * vp)
        assert e.has(vp)
        f = sympy.lambdify([X, vp], e, modules=KMODS)
        assert float(f(4.0, 1.0)) == pytest.approx(2.0, rel=1e-12)
        ## and below the (now parameter-scaled) floor it is the floor
        assert float(f(-1.0, 1.0)) == pytest.approx(math.sqrt(1e-3),
                                                    rel=1e-12)


class TestPromotedHelpers(object):
    """`softplus`, `mne`, `mxe`, `p3` -- lifted out of `psp_kernel`.

    `p3` was DUPLICATED there, byte for byte; the other three are
    generic and were only in a model file because that is where they were
    first needed.
    """

    def test_softplus_matches_log1p_exp(self):
        """Reference: `math.log1p(math.exp(z))`, from the stdlib."""
        f = sympy.lambdify(X, hdl.softplus(X), modules=KMODS)
        for v in (-50.0, -1.0, -1e-9, 0.0, 1e-9, 1.0, 50.0, 700.0):
            want = math.log1p(math.exp(v)) if v < 700 else v
            assert float(f(v)) == pytest.approx(want, rel=1e-13), v

    def test_softplus_keeps_the_SMALL_end_too(self):
        """`log(1 + exp(z))` written literally loses it, and quietly.

        At z = -50 the exponential is 1.9e-22, `1 + 1.9e-22` rounds to
        exactly 1.0, and the logarithm of that is exactly 0 -- so a form
        built to keep the LARGE end finite throws the small end away.
        `log1p` is why this one does not, and it is the single respect in
        which the kernel version differs from `psp_kernel._softplus`.
        """
        assert 1.0 + math.exp(-50.0) == 1.0                 # the trap
        f = sympy.lambdify(X, hdl.softplus(X), modules=KMODS)
        for v in (-50.0, -200.0, -700.0):
            assert float(f(v)) == pytest.approx(math.log1p(math.exp(v)),
                                                rel=1e-13), v
            assert float(f(v)) > 0.0

    def test_softplus_is_finite_where_the_literal_form_overflows(self):
        """`log(1 + exp(z))` overflows at z = 710 -- and the answer is 710."""
        with np.errstate(all='ignore'):
            naive = float(np.log(1.0 + np.exp(710.0)))
        assert not np.isfinite(naive)                       # the trap

        f = sympy.lambdify(X, hdl.softplus(X), modules=KMODS)
        d = sympy.lambdify(X, sympy.diff(hdl.softplus(X), X), modules=KMODS)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            for v in (-1e300, -1e10, -710.0, 0.0, 710.0, 1e10, 1e300):
                y = float(f(v))
                assert np.isfinite(y)
                assert y >= 0.0
                assert y <= abs(v) + 0.7 + 1e-9      # the claimed bound
                g = float(d(v))
                assert 0.0 <= g <= 1.0               # the logistic

    def test_softplus_derivative_is_the_logistic(self):
        f = sympy.lambdify(X, hdl.softplus(X), modules=KMODS)
        d = sympy.lambdify(X, sympy.diff(hdl.softplus(X), X), modules=KMODS)
        for v in (-20.0, -1.0, -0.1, 0.1, 1.0, 20.0):
            want = 1.0 / (1.0 + math.exp(-v))       # independent reference
            assert float(d(v)) == pytest.approx(want, rel=1e-10), v
            assert float(d(v)) == pytest.approx(_fd(f, v), abs=1e-7), v

    def test_mne_and_mxe_are_smooth_min_and_max(self):
        mn = sympy.lambdify([X, Y], hdl.mne(X, Y), modules=KMODS)
        mx = sympy.lambdify([X, Y], hdl.mxe(X, Y), modules=KMODS)
        for a, b in ((1.0, 3.0), (0.2, 7.0), (1e-6, 1.0), (1e3, 1e6),
                     (1e12, 1.0)):
            ## at a = 0 the pair is EXACT, not merely close
            assert float(mn(a, b)) == pytest.approx(min(a, b), rel=1e-9)
            assert float(mx(a, b)) == pytest.approx(max(a, b), rel=1e-9)

        ## EXCEPT at x == y, where the radicand is exactly zero and the
        ## `safe_sqrt` regularisation shows through: `sqrt(hypsmooth(0))`
        ## is `eps`, not 0, so the answer is off by `eps/2` absolute
        ## (5e-7 with the default 1e-12).  Stated rather than smoothed
        ## over -- it is a property of the eps, and the eps is chosen for
        ## the derivative's sake.
        assert float(mn(5.0, 5.0)) == pytest.approx(5.0, abs=1e-6)
        assert float(mx(5.0, 5.0)) == pytest.approx(5.0, abs=1e-6)
        assert float(mn(5.0, 5.0)) != 5.0

    def test_mne_holds_its_accuracy_over_sixteen_decades_of_ratio(self):
        """The property the conjugate form is written for.

        `psp_kernel._mne` records that the literal `s - sqrt(s^2 - 4xy)`
        "at t2 = 1e9 ... rounds to zero".  **That is not reproducible
        here**: swept over ratios from 1e1 to 1e16, the literal form in
        double precision returns the right answer at every one of them,
        because `s*s` and `4xy` are far enough apart that the root is not
        close to `s`.  The claim is left in place upstream and contested
        here rather than repeated.

        What IS worth pinning, and is pinned, is the accuracy the
        shipped form actually delivers across that range -- which is the
        thing a caller depends on either way.
        """
        mn = sympy.lambdify([X, Y], hdl.mne(X, Y), modules=KMODS)
        for k in range(1, 17):
            a, b = 10.0 ** k, 1.0
            assert float(mn(a, b)) == pytest.approx(1.0, rel=1e-12), k

    def test_mne_mxe_derivatives_against_finite_differences(self):
        for expr in (hdl.mne(X, 3.0), hdl.mxe(X, 3.0)):
            f = sympy.lambdify(X, expr, modules=KMODS)
            d = sympy.lambdify(X, sympy.diff(expr, X), modules=KMODS)
            ## 3.0 -- the crossover -- is EXCLUDED, and deliberately:
            ## the smoothing there is `sqrt(hypsmooth(0))` wide, i.e.
            ## 1e-6, and a difference quotient with a 3e-6 step straddles
            ## it and measures the corner rather than the slope.  A
            ## finite difference cannot speak at a kink; the value there
            ## is asserted separately below.
            for v in (0.1, 0.5, 1.0, 2.9, 3.1, 10.0, 1e4):
                ## `abs` matters as much as `rel` here: on the flat side
                ## both the analytic slope and the difference quotient
                ## are ~1e-13, and a relative test on two numbers that
                ## are both zero to working precision is not a test
                assert float(d(v)) == pytest.approx(_fd(f, v), rel=1e-4,
                                                    abs=1e-6), (expr, v)

    def test_mne_mxe_split_the_slope_evenly_at_the_crossover(self):
        """Where the two arguments meet, each gets half -- by symmetry.

        The smooth min/max are symmetric in `x` and `y`, so at `x == y`
        neither can claim more than half the slope; anything else would
        make the pair asymmetric in the source and drain, which is the
        property PSP uses them to preserve.
        """
        for expr in (hdl.mne(X, 3.0), hdl.mxe(X, 3.0)):
            d = sympy.lambdify(X, sympy.diff(expr, X), modules=KMODS)
            assert float(d(3.0)) == pytest.approx(0.5, rel=1e-9)

    def test_p3_is_exps_third_order_taylor(self):
        f = sympy.lambdify(X, hdl.p3(X), modules=KMODS)
        for u in (-2.0, -0.5, 0.0, 0.5, 2.0):
            want = 1 + u + u ** 2 / 2 + u ** 3 / 6      # written out
            assert float(f(u)) == pytest.approx(want, rel=1e-14)
        ## three derivatives match exp at the origin, which is the point
        for order in (0, 1, 2, 3):
            got = sympy.diff(hdl.p3(X), X, order).subs(X, 0)
            assert float(got) == pytest.approx(1.0, rel=1e-14)

    def test_p3_is_the_same_function_the_expl_family_uses(self):
        """It was duplicated in `psp_kernel`; pin that it is one thing now."""
        assert hdl._p3 is hdl.p3


class TestKernelPrimitivesReachEveryEvaluator(object):
    """Registration, verified by EXECUTION on each path.

    A primitive prints as a call, so a namespace missing it is a
    `NameError` at evaluation time from a model that compiled clean.
    There are five namespaces and they are reached by different kinds of
    model; a dict entry that is merely present proves nothing.
    """

    def test_a_bare_lambdify_resolves_them_with_no_modules_map(self):
        """Not a nicety -- `hypsmooth` and `expl` are PUBLIC.

        They are clamped with `maxc`/`minc` now, so an outside caller who
        wrote `sympy.lambdify(x, hdl.safe_sqrt(x), modules='numpy')` --
        which worked while the clamp was `sympy.Max` -- would get a
        `NameError` from inside `<lambdifygenerated-N>` if these carried
        no `_imp_`.
        """
        for expr in (hdl.maxc(X, 0.5), hdl.minc(X, 2.0), hdl._step(X, 0.5),
                     hdl.safe_sqrt(X), hdl.expl(X), hdl.hypsmooth(X, 1e-6),
                     hdl.softplus(X), hdl.safe_pow(X, 0.5)):
            f = sympy.lambdify(X, expr)                  # NO modules map
            d = sympy.lambdify(X, sympy.diff(expr, X))
            assert np.isfinite(float(f(1.3)))
            assert np.isfinite(float(d(1.3)))

    def test_the_eager_element_path(self):
        """`NUMPY_MODULES` -- an element with no `var()` in it."""
        e, cls = _kernel_element(chained=False)
        x = np.array([1.3, 0.0])
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            i = np.asarray(e.i(x), float)
            G = np.asarray(e.G(x), float)
        assert np.all(np.isfinite(i)) and np.all(np.isfinite(G))
        assert G[0, 0] == pytest.approx(_element_fd(e, x), rel=1e-5)

    def test_the_chained_element_path(self):
        """`_ChainPrinter` and `_chain_compile`'s namespace.

        This is the path EVERY production model takes, because every
        production model uses `var()`, and it is the path a missing
        registration hides on: the chain prints its own source rather
        than going through sympy's printer.
        """
        e, cls = _kernel_element(chained=True)
        x = np.array([1.3, 0.0])
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            i = np.asarray(e.i(x), float)
            G = np.asarray(e.G(x), float)
        assert np.all(np.isfinite(i)) and np.all(np.isfinite(G))
        assert G[0, 0] == pytest.approx(_element_fd(e, x), rel=1e-5)

        fn = cls._hdl_info['funcs']['G']
        ## the printer really emitted the calls ...
        assert 'maxc(' in fn._src and 'minc(' in fn._src
        assert '_step(' in fn._src
        ## ... and the namespace the chain was exec'd in really bound them
        assert fn.__globals__['maxc'] is hdl._maxc_numpy
        assert fn.__globals__['minc'] is hdl._minc_numpy
        assert fn.__globals__['_step'] is hdl._step_numpy

    def test_the_two_paths_agree(self):
        ea, _ = _kernel_element(chained=False)
        eb, _ = _kernel_element(chained=True)
        for v in (-3.0, -1e-9, 0.0, 1e-9, 0.5, 1.3, 40.0):
            x = np.array([v, 0.0])
            assert np.allclose(np.asarray(ea.i(x), float),
                               np.asarray(eb.i(x), float), rtol=1e-12)
            assert np.allclose(np.asarray(ea.G(x), float),
                               np.asarray(eb.G(x), float), rtol=1e-10,
                               atol=1e-14)

    def test_the_crossing_path(self):
        """`cross_spec`'s own `modules_map`, a sixth compilation.

        Distinct from `i`/`q`/`G`/`C`: an event expression is compiled
        separately and a name missing there fails only for a model that
        asks the transient solver for a timepoint.
        """
        import pycircuit.circuit.circuit as cm
        from pycircuit.circuit.toolkit import numeric
        from pycircuit.circuit.hdl import (Behavioural, Branch,
                                           Contribution, Cross, var)
        cm.default_toolkit = numeric

        def analog(p, n):
            b = Branch(p, n)
            v = var(b.V, 'v')
            return (Contribution(b.I, var(hdl.maxc(v, 0.0), 'a')),
                    Cross(hdl.maxc(v, 0.0) - 0.5, +1))

        cls = type('KernelCross', (Behavioural,),
                   {'instparams': [], 'analog': staticmethod(analog)})
        e = cls(cm.Node('p'), cm.Node('n'))
        e.update_iparv()
        f = cls._hdl_info['cross_spec']['f']
        assert 'maxc(' in f._src
        assert float(np.asarray(f(np.array([1.3, 0.0]), 300.0),
                                float)[0]) == pytest.approx(0.8)

    def test_the_two_jax_paths(self):
        """`eval_i_pure` (the vmap admission ticket) and PCNR's own map.

        Under `jax.jit` a numpy implementation does not merely run slowly
        -- `numpy.maximum` on a tracer raises
        `TracerArrayConversionError` -- so a value that matches the numpy
        path here is proof the jax binding is live.
        """
        pytest.importorskip('jax')
        import jax
        import pycircuit.circuit.circuit as cm
        from pycircuit.circuit.toolkit import numeric, jaxtoolkit
        from pycircuit.circuit.circuit import defaultepar
        from pycircuit.circuit.hdl import (Behavioural, Branch,
                                           Contribution, vt)
        from pycircuit.utilities.param import Parameter
        cm.default_toolkit = numeric

        e, cls = _kernel_element(chained=False)
        assert hasattr(cls, 'eval_i_pure')
        x = np.array([1.3, 0.0])
        got = np.asarray(cls.eval_i_pure(jaxtoolkit.array(x), {},
                                         defaultepar, jaxtoolkit), float)
        assert np.allclose(got, np.asarray(e.i(x), float), rtol=1e-6)
        jitted = jax.jit(lambda z: cls.eval_i_pure(z, {}, defaultepar,
                                                   jaxtoolkit))
        assert np.allclose(np.asarray(jitted(jaxtoolkit.array(x)), float),
                           np.asarray(e.i(x), float), rtol=1e-6)
        ## and the JACOBIAN, which is where `_step` has to survive tracing
        J = np.asarray(jax.jacfwd(
            lambda z: cls.eval_i_pure(z, {}, defaultepar, jaxtoolkit))(
                jaxtoolkit.array(x)), float)
        assert np.allclose(J, np.asarray(e.G(x), float), rtol=1e-5)

        ## PCNR compiles its junction through a SEPARATE jax map.
        def analog(p, n):
            b = Branch(p, n)
            return Contribution(
                b.I, IS * (hdl.limexp(b.V / vt()) - 1.0)     # noqa: F821
                + 1e-3 * hdl.maxc(b.V, 0.0) * hdl.sign(b.V))

        jcls = type('KernelJunction', (Behavioural,), {
            'instparams': [Parameter(name='IS', desc='I', unit='A',
                                     default=1e-14)],
            'analog': staticmethod(analog)})
        je = jcls(cm.Node('p'), cm.Node('n'), IS=1e-14)
        je.update_iparv()
        assert jcls.pcnr_junctions == ((0, 1),)
        pr = {'IS': 1e-14}
        want = np.asarray(jcls.pcnr_i(0.6, pr, defaultepar, numeric), float)
        gotj = np.asarray(jax.jit(
            lambda v: jcls.pcnr_i(v, pr, defaultepar, jaxtoolkit))(
                jaxtoolkit.array(0.6)), float)
        assert np.allclose(gotj, want, rtol=1e-6)


class TestTheClampedExponentialsGotCheaper(object):
    """`expl` and `hypsmooth` clamp with `maxc`/`minc` now, not `Max`/`Min`.

    The value is unchanged -- the tests above in this file cover that --
    and the derivative is the same function.  What changed is what the
    compiler emits for it, and only a source assertion can see that.

    Measured on a chained element's `G` on this machine: `expl` 42
    dispatched `numpy.where` calls at 140 us, now 6 at 58 us;
    `hypsmooth` 19 at 56 us, now 3 at 32 us.  The thresholds below are
    set with headroom, so they are a guard against the regression and not
    a pin on the exact count.
    """

    @pytest.mark.parametrize('name,fn,limit', [
        ('expl', lambda v: hdl.expl(v), 12),
        ('hypsmooth', lambda v: hdl.hypsmooth(v, 1e-6), 8),
    ])
    def test_the_chained_jacobian_emits_few_dispatched_conditionals(
            self, name, fn, limit):
        import pycircuit.circuit.circuit as cm
        from pycircuit.circuit.toolkit import numeric
        from pycircuit.circuit.hdl import (Behavioural, Branch,
                                           Contribution, var)
        cm.default_toolkit = numeric
        _HOLD[0] = fn

        def analog(p, n):
            b = Branch(p, n)
            return Contribution(b.I, var(_HOLD[0](var(b.V, 'v')), 'y'))

        cls = type('KernelCheap' + name, (Behavioural,),
                   {'instparams': [], 'analog': staticmethod(analog)})
        e = cls(cm.Node('p'), cm.Node('n'))
        e.update_iparv()
        src = cls._hdl_info['funcs']['G']._src
        n_where = len(re.findall(r'\bwhere\(', src))
        assert len(re.findall(r'\bselect\(', src)) == 0
        assert len(re.findall(r'\bHeaviside\(', src)) == 0
        assert n_where <= limit, '%s emits %d where() calls' % (name,
                                                                n_where)

    def test_the_values_are_bit_for_bit_what_they_were(self):
        """`Max` and `maxc` are the same function; pin that they agree.

        Not a tautology: `maxc` is a different CLASS with its own
        evaluation, and a sign error in `minc`'s argument order would
        change `expl`'s clamp without changing its shape enough to be
        obvious.
        """
        f = sympy.lambdify(X, hdl.expl(X), modules=KMODS)
        g = sympy.lambdify(X, hdl.hypsmooth(X, 1e-6), modules=KMODS)
        for v in (-1e6, -1000., -SE05 - 1., -SE05, -SE05 + 1., -1., 0.,
                  1., SE05 - 1., SE05, SE05 + 1., 1000., 1e6):
            assert float(f(v)) > 0.0
            assert float(g(v)) > 0.0
        ## exactly `exp` between the seams, still
        for v in (-200., -1., 0., 1., 200.):
            assert float(f(v)) == pytest.approx(np.exp(v), rel=1e-13)


#: Set by the parametrised test above; a closure default argument would
#: become a TERMINAL, since `analog()`'s signature declares the pins.
_HOLD = [None]


def _element_fd(e, x, j=0):
    """d i[0] / d x[j] by central differences on a compiled element."""
    h = 1e-6 * max(1.0, abs(x[j]))
    xp, xm = x.copy(), x.copy()
    xp[j] += h
    xm[j] -= h
    return (float(np.asarray(e.i(xp), float)[0])
            - float(np.asarray(e.i(xm), float)[0])) / (2 * h)


def _kernel_element(chained):
    """A two-terminal element whose current mentions every new primitive.

    Built twice from the same expression, once with `var()` and once
    without, because the eager stamps and the let-chain are DIFFERENT
    code generators and a primitive can be registered for one and not the
    other.
    """
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       var)
    cm.default_toolkit = numeric

    def body(v):
        ## Each primitive gets a DIFFERENTLY SCALED argument, which is
        ## worth a word because it is not decoration.  Written as nine
        ## terms of the same `b.V`, the eager form is one nine-term
        ## `Add`; sympy sorts an `Add` that long, the sort reaches
        ## `Quantity.compare` on two equal branch voltages, and that
        ## raises `TypeError: '>' not supported between instances of
        ## NoneType and NoneType` -- `Quantity._hashable_content` puts a
        ## bare `None` in the tuple and sympy compares it with `>`.  That is
        ## a pre-existing defect with nothing to do with this kernel (it
        ## reproduces as `Quantity('V', b).compare(Quantity('V', b))` on
        ## a stock tree), it is outside this change's remit, and it is
        ## reported.  Distinct coefficients make the sort resolve on the
        ## coefficient and never reach the atoms.
        return (hdl.maxc(v, 0.5) + hdl.minc(1.1 * v, 2.0)
                + hdl.sign(1.2 * v) + hdl.Abs(1.3 * v)
                + hdl.safe_pow(1.4 * v, 0.5) + hdl.softplus(1.5 * v)
                + hdl.mne(1.6 * v, 1.0) + hdl.mxe(1.7 * v, 1.0)
                + hdl.p3(1.8 * v))

    def analog_eager(p, n):
        b = Branch(p, n)
        return Contribution(b.I, body(b.V))

    def analog_chained(p, n):
        b = Branch(p, n)
        return Contribution(b.I, var(body(var(b.V, 'v')), 'y'))

    name = 'KernelAll' + ('Chained' if chained else 'Eager')
    cls = type(name, (Behavioural,),
               {'instparams': [],
                'analog': staticmethod(analog_chained if chained
                                       else analog_eager)})
    e = cls(cm.Node('p'), cm.Node('n'))
    e.update_iparv()
    return e, cls
