"""PSP's surface-potential kernel, checked against the equation it solves.

`psp_kernel.sp_s` is a translation of PSP103's explicit surface-potential
approximation -- the heart of the model, and the piece whose correctness
everything else rests on.

It is checked WITHOUT a vendor binary or a model card, by evaluating the
surface-potential equation at the returned root.  That residual is not
taken from a paper: it is the expression PSP's own final correction step
computes as its ``qC`` term, so the test and the code are checked against
the same source rather than against each other.

A second, independent check: the DSL's forward-accumulated ``d(sp)/d(xg)``
is compared against the EXACT implicit derivative ``-F_xg / F_x`` of that
equation.  Those two agreeing means both the translation and the
compiler's Jacobian are right, by different routes.
"""
import math

import numpy as np
import pytest
import sympy

from pycircuit.circuit import hdl
from pycircuit.circuit import psp_kernel as K


XG, XN, GF, XI = sympy.symbols('xg xn Gf xi', real=True)
X = sympy.Symbol('x', real=True)
MODULES = [{'Min': np.minimum, 'Max': np.maximum}, 'numpy']

#: PSP's own definition, `mosvar.va:584`: xi = 1 + Gf/sqrt(2).
INV_SQRT2 = 7.0710678118654746e-01


def _xi(gf):
    return 1.0 + gf * INV_SQRT2


#: The physically meaningful domain.
#:
#: `xn` is the normalised quasi-Fermi splitting: 2*phi_F/vt is about 35
#: for a real MOSFET at room temperature.  `Gf` is the normalised body
#: factor, gamma/sqrt(vt), which for gamma = 0.3..1.0 V^0.5 lands around
#: 2..6.  The kernel is checked over a generous envelope of that, and its
#: behaviour OUTSIDE it is pinned separately as a known limit rather than
#: quietly avoided.
GF_PHYS = (0.3, 1.0, 3.0, 5.0)
XN_PHYS = (15.0, 25.0, 35.0, 45.0)


@pytest.fixture(scope='module')
def fns():
    """Compile `sp_s`, its xg-derivative, and the residual."""
    build = lambda: K.sp_s(XG, XN, sympy.exp(-XN), GF, XI)
    val = hdl.compile_chain(build, [XG, XN, GF, XI])
    jac = hdl.compile_chain(build, [XG, XN, GF, XI], wrt=[XG])
    res = K.spe_residual(X, XG, XN, GF)
    return dict(
        sp=lambda *a: float(val(*a)[0]),
        dsp=lambda *a: float(jac(*a)[0][0]),
        F=sympy.lambdify((X, XG, XN, GF), res, MODULES),
        F_x=sympy.lambdify((X, XG, XN, GF), sympy.diff(res, X), MODULES),
        F_xg=sympy.lambdify((X, XG, XN, GF), sympy.diff(res, XG), MODULES),
    )


def _xg_sweep():
    return np.concatenate([np.linspace(-100.0, -1e-4, 40),
                           np.array([-1e-6, 0.0, 1e-6]),
                           np.linspace(1e-4, 300.0, 80)])


def _scale(xg, s, gf):
    """Size of the terms that cancel in F, for a relative measure."""
    return max(abs(xg - s) ** 2, gf ** 2 * max(abs(s), 1.0), 1.0)


class TestItSolvesTheEquation(object):

    @pytest.mark.parametrize('gf', GF_PHYS)
    @pytest.mark.parametrize('xn', XN_PHYS)
    def test_residual_is_negligible_over_the_physical_domain(
            self, fns, gf, xn):
        """The whole claim, in one assertion.

        Measured worst over this envelope is 6.7e-8 and typically ~2e-9,
        which is the accuracy PSP's explicit approximation is built for.
        """
        xi = _xi(gf)
        worst, at = 0.0, None
        for xg in _xg_sweep():
            s = fns['sp'](xg, xn, gf, xi)
            r = abs(float(fns['F'](s, xg, xn, gf))) / _scale(xg, s, gf)
            if r > worst:
                worst, at = r, xg
        assert worst < 1e-6, 'worst %.2e at xg=%.4g' % (worst, at)

    def test_it_is_monotonic_in_the_gate_drive(self, fns):
        """Surface potential rises with gate drive, everywhere."""
        xi = _xi(1.0)
        xs = _xg_sweep()
        sp = np.array([fns['sp'](x, 35.0, 1.0, xi) for x in xs])
        assert np.all(np.diff(sp) > 0)

    def test_it_is_finite_far_outside_any_sane_bias(self, fns):
        """Newton visits absurd points on its way to the solution."""
        xi = _xi(1.0)
        for xg in (-1e4, -1e3, 1e3, 1e4, 1e6):
            s = fns['sp'](xg, 35.0, 1.0, xi)
            assert np.isfinite(s), xg

    def test_flat_band_gives_zero_surface_potential(self, fns):
        assert fns['sp'](0.0, 35.0, 1.0, _xi(1.0)) == pytest.approx(0.0,
                                                                    abs=1e-12)

    def test_the_three_regimes_join_continuously(self, fns):
        """`|xg| <= margin` switches arms; the seam must not show.

        Continuity means the gap across the seam shrinks WITH the
        interval.  Measured, it is ~0.586*2h -- exactly the function's
        own slope, so there is no step at all.  (A fixed step would stay
        constant as h shrinks; that is what this distinguishes.)
        """
        xi, m = _xi(1.0), 1e-5
        for seam in (-m, m):
            prev = None
            for h in (1e-6, 1e-7, 1e-8, 1e-9):
                gap = abs(fns['sp'](seam + h, 35.0, 1.0, xi)
                          - fns['sp'](seam - h, 35.0, 1.0, xi))
                assert gap < 3.0 * h, 'seam %g h %g gap %.3e' % (seam, h,
                                                                 gap)
                if prev is not None:
                    ## shrinking by ~10x per decade of h, i.e. linear
                    assert gap < 0.2 * prev
                prev = gap

    def test_the_slope_is_continuous_across_the_seam(self, fns):
        """C1, checked just OUTSIDE the seam on each side.

        Not AT it: `sp` is built from `Max`/`Min` clamps, and sympy
        differentiates those through `Heaviside`, whose value at zero is
        1/2 by convention.  So the derivative evaluated exactly at a
        clamp point is the AVERAGE of the one-sided derivatives -- an
        artifact of asking at the corner, not a property of the model.
        Pinned here because it cost a debugging session.
        """
        xi, m = _xi(1.0), 1e-5
        for seam in (-m, m):
            h = 1e-8
            lo = fns['dsp'](seam - h, 35.0, 1.0, xi)
            hi = fns['dsp'](seam + h, 35.0, 1.0, xi)
            assert lo == pytest.approx(hi, rel=1e-6)

    def test_the_derivative_exactly_at_a_clamp_is_the_average(self, fns):
        """The `Heaviside(0) = 1/2` artifact itself, pinned."""
        xi, m = _xi(1.0), 1e-5
        h = 1e-9
        one_sided = fns['dsp'](-m + h, 35.0, 1.0, xi)
        at_seam = fns['dsp'](-m, 35.0, 1.0, xi)
        assert at_seam == pytest.approx(0.5 * one_sided, rel=1e-3)


class TestTheDerivativeIsRight(object):
    """Forward accumulation vs the exact implicit derivative."""

    @pytest.mark.parametrize('gf', GF_PHYS)
    def test_matches_implicit_differentiation(self, fns, gf):
        """`dx/dxg = -F_xg / F_x`, exactly, wherever F_x is nonzero.

        Skipped near flat band, where F has a DOUBLE root: F_x vanishes
        and the implicit derivative is genuinely undefined there, so a
        disagreement would say nothing about the code.
        """
        xi, xn = _xi(gf), 35.0
        worst, at, n = 0.0, None, 0
        for xg in _xg_sweep():
            s = fns['sp'](xg, xn, gf, xi)
            fx = float(fns['F_x'](s, xg, xn, gf))
            if abs(fx) <= 1e-6 * max(abs(xg - s), gf ** 2, 1.0):
                continue
            want = -float(fns['F_xg'](s, xg, xn, gf)) / fx
            err = abs(fns['dsp'](xg, xn, gf, xi) - want) / max(abs(want),
                                                               1e-3)
            n += 1
            if err > worst:
                worst, at = err, xg
        assert n > 100, 'too few usable points (%d)' % n
        assert worst < 1e-5, 'worst %.2e at xg=%.4g' % (worst, at)

    def test_it_also_matches_finite_differences(self, fns):
        """A cruder check, but one that shares no algebra with the above."""
        xi, xn, gf = _xi(1.0), 35.0, 1.0
        for xg in (-50., -5., 1., 10., 50., 150.):
            h = 1e-5 * max(1.0, abs(xg))
            fd = (fns['sp'](xg + h, xn, gf, xi)
                  - fns['sp'](xg - h, xn, gf, xi)) / (2 * h)
            assert fns['dsp'](xg, xn, gf, xi) == pytest.approx(fd, rel=1e-5)


class TestKnownLimits(object):
    """Where the approximation stops being accurate, pinned as a fact.

    Both are outside the domain PSP is fitted for, and both degrade
    smoothly rather than failing.  Recorded so that a future measurement
    landing here is recognised as the approximation's edge and not as a
    regression.
    """

    def test_it_degrades_as_the_body_factor_vanishes(self, fns):
        """`Gf -> 0` makes the equation degenerate.

        `F = (xg - x)^2 - Gf^2*[...]`, so as `Gf -> 0` the root becomes a
        double root at `xg` and the problem is ill-conditioned.  A real
        body factor is 2..6; 0.1 means gamma ~ 0.016 V^0.5.
        """
        got = {}
        for gf in (5.0, 1.0, 0.3, 0.1):
            xi, worst = _xi(gf), 0.0
            for xg in _xg_sweep():
                s = fns['sp'](xg, 35.0, gf, xi)
                worst = max(worst, abs(float(fns['F'](s, xg, 35.0, gf)))
                            / _scale(xg, s, gf))
            got[gf] = worst
        assert got[5.0] < 1e-6 and got[1.0] < 1e-6 and got[0.3] < 1e-6
        assert got[0.1] > 1e-6           # the edge, and it is smooth
        assert got[0.1] < 1e-3

    def test_it_degrades_at_zero_quasi_fermi_splitting(self, fns):
        """`xn = 0` means delta = 1: no inversion barrier at all.

        A real device has xn ~ 35.  At xn = 0 the residual reaches ~1e-5
        where the physical domain gives ~1e-9.
        """
        xi = _xi(1.0)
        worst = max(abs(float(fns['F'](fns['sp'](xg, 0.0, 1.0, xi),
                                       xg, 0.0, 1.0)))
                    / _scale(xg, fns['sp'](xg, 0.0, 1.0, xi), 1.0)
                    for xg in _xg_sweep())
        assert worst > 1e-7
        assert worst < 1e-3


class TestTheHelperFunctions(object):

    def test_mina_is_a_smooth_minimum(self):
        x, y = sympy.symbols('x y', real=True)
        f = sympy.lambdify((x, y), K.mina(x, y, 5.0), MODULES)
        ## well away from the corner it is the true minimum
        assert float(f(0.0, 100.0)) == pytest.approx(0.0, abs=0.02)
        assert float(f(100.0, 0.0)) == pytest.approx(0.0, abs=0.02)
        ## and it never exceeds it
        for a, b in ((1.0, 2.0), (-3.0, 5.0), (0.0, 0.0), (7.0, 7.0)):
            assert float(f(a, b)) <= min(a, b) + 1e-12

    def test_sigma_and_sigma2_agree_when_b_is_one(self):
        a, c, tau, eta = sympy.symbols('a c tau eta', real=True)
        f1 = sympy.lambdify((a, c, tau, eta), K.sigma(a, c, tau, eta),
                            MODULES)
        f2 = sympy.lambdify((a, c, tau, eta),
                            K.sigma2(a, 1.0, c, tau, eta), MODULES)
        for args in ((2.0, 3.0, 0.5, 1.0), (1e-3, 0.2, -2.0, 0.1),
                     (5.0, -1.0, 3.0, -2.0)):
            assert float(f1(*args)) == pytest.approx(float(f2(*args)),
                                                     rel=1e-12)

    def test_sigma2_takes_the_vendors_tiny_tau_shortcut(self):
        a, b, c, tau, eta = sympy.symbols('a b c tau eta', real=True)
        f = sympy.lambdify((a, b, c, tau, eta),
                           K.sigma2(a, b, c, tau, eta), MODULES)
        assert float(f(2.0, 1.0, 3.0, 1e-200, 7.0)) == 7.0

    def test_the_residual_is_zero_at_flat_band(self):
        f = sympy.lambdify((X, XG, XN, GF),
                           K.spe_residual(X, XG, XN, GF), MODULES)
        assert float(f(0.0, 0.0, 35.0, 1.0)) == pytest.approx(0.0, abs=1e-18)


class TestUsableInsideAModel(object):
    """The kernel has to survive the element compiler, not just lambdify."""

    def test_a_model_using_sp_s_compiles_and_differentiates(self):
        import pycircuit.circuit.circuit as cm
        from pycircuit.circuit.toolkit import numeric
        from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                           Node, var, vt)
        from pycircuit.utilities.param import Parameter

        cm.default_toolkit = numeric

        class SurfacePotentialProbe(Behavioural):
            """A toy device whose current is the surface potential.

            Not physics -- a harness. The point is that `sp_s` compiles
            inside `analog()`, lands in the stamps, and differentiates.
            """
            instparams = [Parameter(name='gf', desc='body factor', unit='',
                                    default=1.0),
                          Parameter(name='xnq', desc='xn', unit='',
                                    default=35.0),
                          Parameter(name='k', desc='scale', unit='A',
                                    default=1e-6)]

            @staticmethod
            def analog(g, b):
                br = Branch(g, b)
                xg = var(br.V / vt(), 'xg')
                xi = var(1.0 + gf * INV_SQRT2, 'xi')          # noqa: F821
                sp = K.sp_s(xg, xnq, sympy.exp(-xnq), gf, xi)  # noqa: F821
                return Contribution(br.I, k * sp)              # noqa: F821

        e = SurfacePotentialProbe(Node('g'), Node('b'), gf=1.0, xnq=35.0,
                                  k=1e-6)
        e.update_iparv()
        assert SurfacePotentialProbe._hdl_info['chained'] is True

        for v in (-2.0, -0.5, 0.0, 0.5, 1.0, 2.0):
            x = np.array([v, 0.0])
            i = np.asarray(e.i(x), float)
            G = np.asarray(e.G(x), float)
            assert np.all(np.isfinite(i)), v
            assert np.all(np.isfinite(G)), v
            ## and G is the derivative of that i
            h = 1e-6
            xp, xm = np.array([v + h, 0.0]), np.array([v - h, 0.0])
            fd = (np.asarray(e.i(xp), float)
                  - np.asarray(e.i(xm), float)) / (2 * h)
            assert G[:, 0] == pytest.approx(fd, rel=1e-5, abs=1e-14), v

    def test_the_kernel_refuses_to_build_without_a_chain(self):
        """No silent fallback to a flattened expression.

        Building `sp_s` as one expression does not finish in ten minutes,
        so a fallback would look like a hang rather than an error.
        """
        with pytest.raises(RuntimeError):
            K.sp_s(XG, XN, sympy.exp(-XN), GF, XI)
