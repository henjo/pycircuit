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
