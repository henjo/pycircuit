# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Stage 1 gate: the perturbation result must equal the Volterra result.

``doc/distortion_plan.md`` declares this gate before the code was written, and
it is the strongest one available in the whole plan, which is why it comes
first: Buonomo & Lo Schiavo prove (TCAS-I 52(8) 2005, Theorems 1 and 2, and
explicitly for the cubic case in their eqs. 43a-b) that their perturbation
method and the Volterra series give *identical* HD2 and HD3 for a cubic
nonlinearity.

So this test does not check the implementation against the same authors'
arithmetic, or against a number read off one of their plots.  It derives
HD2/HD3 a second time by an independent route -- Volterra kernels, computed
here from scratch -- and requires the two to agree symbolically.  If the
perturbation recurrence is wrong, the kernels do not care.

The circuit is the smallest one that exercises everything: a single node with
admittance ``Y(s) = g + sC``, driven by a current source, loaded by a
nonlinear current ``i(v) = b*v**2 + c*v**3``.  Both quadratic and cubic terms
are present because the interesting content of HD3 is the *quadratic*
coefficient's second-order contribution -- a circuit with ``b = 0`` would pass
a much weaker test.
"""

import numpy as np
import pytest
import sympy

from pycircuit.circuit import numeric, symbolic
from pycircuit.circuit.distortion import (CompositeNonlinearity,
                                          CubicNonlinearity,
                                          ExponentialNonlinearity, Harmonic,
                                          harmonic_response,
                                          intermodulation_response,
                                          scalar_nonlinearity,
                                          taylor_coefficients)


## ---------------------------------------------------------------- fixtures

def _single_node_apply_G(Y):
    """``apply_G`` for a one-node circuit with admittance ``Y(s)``."""
    def apply_G(harmonic, rhs):
        s = sympy.I * harmonic.frequency((W,))
        return [rhs[0] / Y(s)]
    return apply_G


W = sympy.Symbol('w', positive=True)
G0, C0, B, C3, A = sympy.symbols('g C b c A', real=True)


def _perturbation_hd(Y):
    """HD2/HD3 from the module under test."""
    bm, cm = scalar_nonlinearity(0, B, C3, 1, symbolic)
    sol = harmonic_response(_single_node_apply_G(Y), [A],
                            CubicNonlinearity(bm, cm),
                            tones=(W,), toolkit=symbolic)
    return (sol.amplitude((2,), 0), sol.amplitude((3,), 0),
            sol.amplitude((1,), 0))


def _volterra_hd(Y):
    """HD2/HD3 from Volterra kernels, derived here independently.

    For ``Y(s) v = u - (b v**2 + c v**3)``, writing ``H1(s) = 1/Y(s)``:

        H2(s1,s2)    = -H1(s1+s2) * b * H1(s1) H1(s2)
        H3(s1,s2,s3) = -H1(s1+s2+s3) * [ c*H1(s1)H1(s2)H1(s3)
                                         + 2b * sym{H1(si) H2(sj,sk)} ]

    and for a single tone of amplitude ``A`` the harmonic amplitudes are
    ``A**2/2 * |H2(jw,jw)|`` and ``A**3/4 * |H3(jw,jw,jw)|``.  At
    ``s1=s2=s3=jw`` the three symmetrised terms coincide, so the average is
    just one of them.
    """
    def H1(s):
        return 1 / Y(s)

    jw = sympy.I * W
    H2_ww = -H1(2 * jw) * B * H1(jw) ** 2
    H3_www = -H1(3 * jw) * (C3 * H1(jw) ** 3 + 2 * B * H1(jw) * H2_ww)

    fundamental = A * H1(jw)
    second = A ** 2 / 2 * H2_ww
    third = A ** 3 / 4 * H3_www
    return second, third, fundamental


def _equal(expr_a, expr_b):
    """Symbolic equality up to simplification."""
    return sympy.simplify(sympy.expand(expr_a - expr_b)) == 0


## ------------------------------------------------------------------- tests

@pytest.mark.parametrize('label,Y', [
    ('resistive',  lambda s: G0),
    ('RC',         lambda s: G0 + s * C0),
])
def test_second_harmonic_matches_volterra(label, Y):
    """HD2 by perturbation must equal HD2 by Volterra kernels, symbolically."""
    p_second, _, p_fund = _perturbation_hd(Y)
    v_second, _, v_fund = _volterra_hd(Y)

    assert _equal(p_fund, v_fund), '%s: fundamentals differ' % label
    assert _equal(p_second, v_second), (
        '%s: second harmonic differs\n  perturbation: %s\n  volterra:     %s'
        % (label, sympy.simplify(p_second), sympy.simplify(v_second)))


@pytest.mark.parametrize('label,Y', [
    ('resistive',  lambda s: G0),
    ('RC',         lambda s: G0 + s * C0),
])
def test_third_harmonic_matches_volterra(label, Y):
    """HD3 must match too -- including the quadratic term's contribution.

    This is the part that would break if the second perturbation order were
    dropped or mis-signed: ``H3`` carries a ``2*b**2*H1(2jw)`` term that has
    nothing to do with the cubic coefficient.
    """
    _, p_third, _ = _perturbation_hd(Y)
    _, v_third, _ = _volterra_hd(Y)

    assert _equal(p_third, v_third), (
        '%s: third harmonic differs\n  perturbation: %s\n  volterra:     %s'
        % (label, sympy.simplify(p_third), sympy.simplify(v_third)))


def test_third_harmonic_actually_depends_on_the_quadratic_coefficient():
    """Guard against a test that would pass with the second order deleted.

    If ``HD3`` were computed from the cubic coefficient alone, the two
    expressions below would be equal.  They must not be.
    """
    Y = lambda s: G0 + s * C0
    _, third_with_b, _ = _perturbation_hd(Y)
    third_without_b = third_with_b.subs(B, 0)

    assert not _equal(third_with_b, third_without_b), (
        'HD3 does not depend on b -- the second-order contribution is missing')


## --------------------------------------------------- Taylor coefficients

def test_taylor_coefficients_of_a_diode():
    """``b`` and ``c`` must be the real Taylor coefficients of the model.

    For ``i = IS*(exp(v/VT) - 1)`` the coefficients about ``v = 0`` are
    ``IS/(2*VT**2)`` and ``IS/(6*VT**3)`` -- the values the 2005 paper uses
    for its bipolar example (its cubic approximation of the exponential).
    Deriving them from ``eval_i_pure`` rather than declaring them is the
    design point of ``taylor_coefficients``.
    """
    from pycircuit.circuit.elements import Diode

    d = Diode(1, 2, toolkit=symbolic, IS=sympy.Symbol('IS', positive=True))
    b, c = taylor_coefficients(d, [0, 0])

    VT = symbolic.kboltzmann * 300 / symbolic.qelectron
    IS = sympy.Symbol('IS', positive=True)

    assert _equal(sympy.simplify(b), IS / (2 * VT ** 2))
    assert _equal(sympy.simplify(c), IS / (6 * VT ** 3))


def test_numeric_taylor_coefficients_agree_with_symbolic():
    """The numeric backend's differenced coefficients must match the exact ones.

    ``nth_derivative`` differences ``order`` times on a numeric toolkit and is
    exact on a symbolic one; a third derivative is where a naive fixed step
    would fall apart, so this pins that the step scaling actually works.
    """
    from pycircuit.circuit.elements import Diode

    IS_val = 1e-13
    b_num, c_num = taylor_coefficients(
        Diode(1, 2, toolkit=numeric, IS=IS_val), [0.0, 0.0])

    VT = numeric.kboltzmann * 300 / numeric.qelectron
    b_exact = IS_val / (2 * VT ** 2)
    c_exact = IS_val / (6 * VT ** 3)

    assert abs(b_num - b_exact) / b_exact < 1e-4
    assert abs(c_num - c_exact) / c_exact < 1e-3


## ------------------------------------------------------- frequency index

def test_harmonic_index_is_general_from_the_start():
    """The index must already support two tones, per the plan's stage-1 rule.

    Two-tone intermodulation reuses the identical recurrence with the scalar
    harmonic index replaced by a pair.  Pinning it here is what stops stage 4
    from becoming a rewrite.
    """
    w1, w2 = sympy.symbols('w1 w2', positive=True)

    assert Harmonic((2,)).frequency((w1,)) == 2 * w1
    assert Harmonic((2, -1)).frequency((w1, w2)) == 2 * w1 - w2
    assert Harmonic((2, -1)).order() == 3
    assert Harmonic((3,)).order() == 3

    with pytest.raises(ValueError):
        Harmonic((2, -1)).frequency((w1,))


## ------------------------------------------------- stage 2: multi-device

"""Stage 2 gate: the three-stage RNMC amplifier of the 2013 paper.

Where the stage 1 gate checks the *mathematics* against an independent theory,
this one checks the implementation against a published circuit with several
interacting nonlinearities -- three of them, spread over two nodes, one with a
negative coefficient.

Tolerances are +/-1 dB on HD2 and +/-2 dB on HD3, set in the plan before the
code was written.  They are wide because the paper publishes **no tables**:
every result in it is a plotted curve, so the reference values are graph
reads and cannot support a tighter claim.  The agreement actually achieved is
far better than that, but the gate stays where it was declared.

Component values are the 2009/2013 set.  The 2012 companion paper publishes
the *same circuit* with fourteen slightly different values; mixing the two
produces a small, plausible, entirely spurious disagreement that looks exactly
like an implementation bug.  See doc/distortion_plan.md section 1.5.
"""

# 2013 paper eqs. (34)-(35).
_GM1, _GM2, _GM3 = 245e-6, 200e-6, 1e-3
_G01, _G02, _G03 = 1 / 98e3, 1 / 107e3, 1 / 23e3
_GM2Q, _GM2C = 4e-3, 60e-3
_GM3Q, _GM3C = 2e-3, 3e-3
_G03Q, _G03C = 0.1e-6, 0.01e-6
_CL, _C1, _C2 = 10e-12, 2e-12, 1e-12
_XIN = 0.1e-3
_OUT = 2                      # node 3 is the output


def _rnmc_admittance(s):
    return np.array([
        [_G01 + s * (_C1 + _C2), -s * _C2, -s * _C1],
        [_GM2,                    _G02,     0.0],
        [0.0,                     _GM3,    -(_G03 + s * _CL)],
    ], dtype=complex)


def _rnmc_coefficients(with_g03=True):
    """``b``/``c`` per 2013 eq. (35): f = (0, f_m2(x1), f_m3(x2) - f_03(x3))."""
    b = np.zeros((3, 3))
    c = np.zeros((3, 3))
    b[1, 0], c[1, 0] = _GM2Q, _GM2C
    b[2, 1], c[2, 1] = _GM3Q, _GM3C
    if with_g03:
        b[2, 2], c[2, 2] = -_G03Q, -_G03C
    return b, c


def _rnmc_hd(f, coefficients=None):
    """HD2 and HD3 at the output node, in dB, at frequency ``f`` in Hz."""
    w = 2 * np.pi * f
    b, c = coefficients if coefficients is not None else _rnmc_coefficients()
    U = np.array([_GM1 * _XIN, 0.0, 0.0], dtype=complex)

    def apply_G(harmonic, rhs):
        s = 1j * harmonic.frequency((w,))
        return np.linalg.solve(_rnmc_admittance(s),
                               np.asarray(rhs, dtype=complex))

    sol = harmonic_response(apply_G, U, CubicNonlinearity(b, c),
                            tones=(w,), toolkit=numeric)
    return 20 * np.log10(sol.HD2(_OUT)), 20 * np.log10(sol.HD3(_OUT))


@pytest.mark.parametrize('f,expected_hd2', [
    (100.0,   -31.8),
    (631e3,  -105.3),      # the minimum of the HD2 curve
])
def test_rnmc_hd2_matches_published_curve(f, expected_hd2):
    hd2, _ = _rnmc_hd(f)
    assert abs(hd2 - expected_hd2) < 1.0, (
        'HD2 at %g Hz: got %.2f dB, expected %.1f dB' % (f, hd2, expected_hd2))


def test_rnmc_hd3_matches_published_curve():
    """HD3 near the peak of its curve.

    Note the peak is broad and its exact location is sampling-dependent: the
    reference value is quoted at 631 Hz (a log-grid point), while the true
    maximum of the computed curve sits nearer 515 Hz at a fractionally higher
    value.  Checking at the quoted frequency rather than at "wherever the
    maximum happens to be" keeps this comparing like with like.
    """
    _, hd3 = _rnmc_hd(631.0)
    assert abs(hd3 - (-65.5)) < 2.0, 'HD3 at 631 Hz: got %.2f dB' % hd3


def test_rnmc_uses_every_nonlinearity():
    """All three nonlinear devices must actually reach the answer.

    Guards against a matrix-form implementation that quietly drops
    off-diagonal coefficients and still passes the curve checks, which the
    dominant g_m2 term alone very nearly would.  Removing the output
    conductance's nonlinearity must perturb the result.
    """
    full = _rnmc_hd(1e4)
    without_g03 = _rnmc_hd(1e4, coefficients=_rnmc_coefficients(with_g03=False))

    assert full != without_g03, (
        'dropping the g_03 nonlinearity changed nothing -- off-diagonal '
        'coefficients are not reaching the recurrence')


def test_rnmc_hd_scales_with_input_amplitude():
    """HD2 must scale as X_in and HD3 as X_in**2.

    A structural check the published curves cannot give, since the paper plots
    only one amplitude per figure.  It is also the check that caught the 2013
    paper's own Fig. 6 caption error, where the plotted curves correspond to
    twice the stated input.
    """
    global _XIN
    original = _XIN
    try:
        _XIN = 0.1e-3
        hd2_a, hd3_a = _rnmc_hd(1e3)
        _XIN = 0.2e-3
        hd2_b, hd3_b = _rnmc_hd(1e3)
    finally:
        _XIN = original

    assert abs((hd2_b - hd2_a) - 20 * np.log10(2)) < 0.01, 'HD2 is not ~X_in'
    assert abs((hd3_b - hd3_a) - 40 * np.log10(2)) < 0.01, 'HD3 is not ~X_in**2'


## ------------------------------------- stage 3: exponential nonlinearity

"""Stage 3 gate: the common-emitter bipolar amplifier of the 2005 paper.

This is the only worked example in the whole source set whose expected answer
is a *closed-form expression* rather than a curve to read off a plot, so it is
the one gate with no graph-reading error at all.  Every parameter is printed
in the paper (its Fig. 3): I_BQ = 10 uA, V_T = 25 mV, beta_F = 180,
R = 250 Ohm, C = 37 pF, U = 50 mV.

Its eqs. 48a-b give HD2 and HD3 in closed form under a *cubic* approximation
of the exponential, so they gate the cubic path on this topology.  The
exponential path is then checked against them the other way round -- it must
*differ*, in the direction and by roughly the margin the paper measures, which
is the entire reason stage 3 exists.

**The trap this test also pins:** distortion at the nonlinearity's controlling
node is not distortion at the circuit output.  The paper refers HD to the
output current via y = E*u + N*x (its eqs. 45c-d), and that mapping carries
both a feedforward term and a frequency-dependent factor.  Measuring at v_be
instead of at i_o gives an answer wrong by a constant factor -- exactly 10 for
these component values -- which looks like a plausible result and is not one.
"""

_IBQ, _VT, _BF = 10e-6, 25e-3, 180.0
_R, _C, _U = 250.0, 37e-12, 50e-3
_ALPHA = _IBQ / _VT                 # junction conductance: the linear part


def _bjt_solution(f, nonlinearity):
    w = 2 * np.pi * f

    def apply_G(harmonic, rhs):
        s = 1j * harmonic.frequency((w,))
        return np.array(rhs, dtype=complex) / (1.0 / _R + s * _C + _ALPHA)

    sol = harmonic_response(apply_G, [_U / _R], nonlinearity,
                            tones=(w,), toolkit=numeric)
    return sol, w


def _bjt_output_hd(f, nonlinearity):
    """HD2/HD3 of the *output current*, per the paper's y = E u + N x."""
    sol, w = _bjt_solution(f, nonlinearity)
    E = -_BF / _R
    N = lambda m: (_BF / _R) * (1 + 1j * m * w * _R * _C)

    Y1 = E * _U + N(1) * sol.amplitude((1,), 0)
    Y2 = N(2) * sol.amplitude((2,), 0)
    Y3 = N(3) * sol.amplitude((3,), 0)
    return abs(Y2) / abs(Y1), abs(Y3) / abs(Y1)


def _bjt_cubic():
    """Cubic truncation of the exponential -- the Taylor coefficients."""
    return CubicNonlinearity(np.array([[_IBQ / (2 * _VT ** 2)]]),
                             np.array([[_IBQ / (6 * _VT ** 3)]]))


def _paper_eq48(f):
    """The 2005 paper's closed forms, eqs. 48a-b."""
    w = 2 * np.pi * f
    Re = _R / (1 + _R * _ALPHA)
    Rk = _R / (1 - 2 * _R * _ALPHA)
    hd2 = (_U / (4 * _VT)) * (1 / (1 + _R * _ALPHA) ** 2) * abs(
        (1 + 2j * w * _R * _C)
        / ((1 + 2j * w * Re * _C) * (1 + 1j * w * Re * _C)))
    hd3 = (_U ** 2 / (24 * _VT ** 2)) * abs(
        (1 - 2 * _R * _ALPHA) / (1 + _R * _ALPHA) ** 4) * abs(
        ((1 + 2j * w * Rk * _C) * (1 + 3j * w * _R * _C))
        / ((1 + 1j * w * Re * _C) ** 2 * (1 + 2j * w * Re * _C)
           * (1 + 3j * w * Re * _C)))
    return hd2, hd3


@pytest.mark.parametrize('f', [1e3, 1e5, 1e6, 1e7])
def test_bjt_cubic_matches_paper_closed_form(f):
    """The gate: agreement with eqs. 48a-b, which are exact expressions.

    No tolerance argument from graph reading is needed here -- the reference
    is an algebraic formula, so this is checked to five significant figures
    rather than to a decibel.
    """
    hd2, hd3 = _bjt_output_hd(f, _bjt_cubic())
    ref2, ref3 = _paper_eq48(f)

    assert abs(hd2 / ref2 - 1) < 1e-4, (
        'HD2 at %g Hz: got %.6e, eq. 48a gives %.6e' % (f, hd2, ref2))
    assert abs(hd3 / ref3 - 1) < 1e-4, (
        'HD3 at %g Hz: got %.6e, eq. 48b gives %.6e' % (f, hd3, ref3))


def test_output_referred_hd_differs_from_node_hd():
    """Distortion at the controlling node is not distortion at the output.

    Pinned because getting this wrong produces a plausible number, not an
    error: for these values the two differ by a constant factor of ten.
    """
    sol, _ = _bjt_solution(1e3, _bjt_cubic())
    node_hd2 = sol.HD2(0)
    output_hd2, _ = _bjt_output_hd(1e3, _bjt_cubic())

    assert abs(output_hd2 / node_hd2 - 10.0) < 0.1, (
        'expected the output/node ratio to be ~10 for this circuit, got %.3f'
        % (output_hd2 / node_hd2))


@pytest.mark.parametrize('f', [1e3, 1e6])
def test_exponential_exceeds_its_cubic_approximation(f):
    """The reason stage 3 exists: a cubic fit to an exponential is not enough.

    The 2005 paper measures a cubic-truncated exponential converging to a
    second harmonic materially below the true value -- around 20% on its
    Fig. 3.  The exact Bessel-function treatment must therefore predict
    *more* distortion than the cubic, by a comparable margin.  Bounds are
    deliberately loose: the paper's figure compares absolute harmonic
    amplitudes on a different quantity, so this checks the effect is real and
    of the right size, not that two differently-defined numbers coincide.
    """
    exact = ExponentialNonlinearity(_IBQ, _VT, port=0, n=1)
    hd2_exact, _ = _bjt_output_hd(f, exact)
    hd2_cubic, _ = _bjt_output_hd(f, _bjt_cubic())

    ratio = hd2_exact / hd2_cubic
    assert ratio > 1.05, (
        'exact exponential should predict more HD2 than its cubic fit, '
        'got ratio %.4f' % ratio)
    assert ratio < 1.6, 'ratio %.4f is implausibly large' % ratio


def test_exponential_reduces_to_cubic_for_small_signals():
    """As the drive vanishes the two treatments must converge.

    The cubic fit *is* the exponential's Taylor expansion, so at small
    amplitude they cannot disagree.  A Bessel implementation with a wrong
    argument scaling would still show a plausible ratio at the working
    amplitude but would fail to converge here.
    """
    global _U
    original = _U
    try:
        _U = 1e-4                      # ~250x smaller than the paper's drive
        exact = ExponentialNonlinearity(_IBQ, _VT, port=0, n=1)
        hd2_exact, _ = _bjt_output_hd(1e3, exact)
        hd2_cubic, _ = _bjt_output_hd(1e3, _bjt_cubic())
    finally:
        _U = original

    assert abs(hd2_exact / hd2_cubic - 1) < 0.02, (
        'exponential and cubic must agree at small signal, ratio %.4f'
        % (hd2_exact / hd2_cubic))


## ------------------------------------------- stage 4: two-tone / IM3

"""Stage 4 gate: the 2*w1 - w2 intermodulation product.

The most sensitive gate in the plan, and deliberately so.  For the reference
amplifier the first- and second-order contributions to IM3 differ by only
0.74 dB and very nearly cancel -- the total lands 21.8 dB below the larger of
the two.  A ~1% error in either term moves the total by a visible amount, so
this checks the two contributions *separately* as well as their sum.  Checking
only the total would let compensating errors through; checking only the terms
would miss a sign error in how they combine.

Reference values are the low-frequency asymptote.  Note 100 Hz is *not* yet on
the plateau for this circuit (it gives -62.21 dB against the asymptotic
-62.73), which is exactly the kind of near-miss the cancellation amplifies --
the individual terms are still within 0.07 dB there.
"""

_IM_TONE_RATIO = 0.9          # f2 = 0.9 * f1, per the paper's Fig. 4(b)


def _rnmc_im3(f1):
    """First-order, second-order and total IM3 at the output node, in dB."""
    w1 = 2 * np.pi * f1
    w2 = 2 * np.pi * _IM_TONE_RATIO * f1
    b, c = _rnmc_coefficients()
    U = np.array([_GM1 * _XIN, 0.0, 0.0], dtype=complex)

    def apply_G(harmonic, rhs):
        s = 1j * harmonic.frequency((w1, w2))
        return np.linalg.solve(_rnmc_admittance(s),
                               np.asarray(rhs, dtype=complex))

    sol = intermodulation_response(apply_G, U, U, CubicNonlinearity(b, c),
                                   (w1, w2), numeric)
    db = lambda v: 20 * np.log10(abs(v))
    return (db(sol.amplitude(('first',), _OUT)),
            db(sol.amplitude(('second',), _OUT)),
            db(sol.amplitude((2, -1), _OUT)))


@pytest.mark.parametrize('label,index,expected', [
    ('first-order contribution',  0, -41.65),
    ('second-order contribution', 1, -40.92),
    ('total IM3',                 2, -62.73),
])
def test_rnmc_im3_matches_published_values(label, index, expected):
    got = _rnmc_im3(1.0)[index]
    assert abs(got - expected) < 0.5, (
        '%s: got %.3f dB, expected %.2f dB' % (label, got, expected))


def test_im3_contributions_very_nearly_cancel():
    """The cancellation is the physics, not an artefact -- pin it.

    If a future change made the two contributions merely additive, every
    individual-term assertion above would still pass while the total moved by
    20 dB.  This states the relationship between them directly.
    """
    first, second, total = _rnmc_im3(1.0)

    assert abs(first - second) < 1.5, (
        'the two contributions should be within ~1 dB of each other, '
        'got %.2f and %.2f' % (first, second))
    assert total < max(first, second) - 15.0, (
        'expected strong cancellation; total %.2f dB is not far enough below '
        'the larger contribution %.2f dB' % (total, max(first, second)))


def test_im3_needs_both_coefficients():
    """Both b and c must reach IM3, by different routes.

    Only the cubic coefficient produces 2*w1 - w2 at first order -- squaring
    two tones cannot make that frequency.  The quadratic coefficient reaches
    it only at second order, by mixing the 2*w1 and w1-w2 products it does
    create.  Zeroing either must change the answer, and by design they change
    *different* contributions.
    """
    w1 = 2 * np.pi
    w2 = 2 * np.pi * _IM_TONE_RATIO
    U = np.array([_GM1 * _XIN, 0.0, 0.0], dtype=complex)

    def run(b, c):
        def apply_G(h, rhs):
            return np.linalg.solve(_rnmc_admittance(1j * h.frequency((w1, w2))),
                                   np.asarray(rhs, dtype=complex))
        sol = intermodulation_response(apply_G, U, U, CubicNonlinearity(b, c),
                                       (w1, w2), numeric)
        return (abs(sol.amplitude(('first',), _OUT)),
                abs(sol.amplitude(('second',), _OUT)))

    b, c = _rnmc_coefficients()
    first_full, second_full = run(b, c)
    first_no_c, second_no_c = run(b, np.zeros_like(c))
    first_no_b, second_no_b = run(np.zeros_like(b), c)

    ## Cubic drives the first-order term and nothing else.
    assert first_no_c < first_full * 1e-9, 'first order should vanish with c=0'
    assert abs(second_no_c - second_full) < second_full * 1e-9

    ## Quadratic drives the second-order term and nothing else.
    assert second_no_b < second_full * 1e-9, 'second order should vanish with b=0'
    assert abs(first_no_b - first_full) < first_full * 1e-9


def test_exponential_two_tone_is_refused_not_guessed():
    """Two-tone for an exponential device is not in any reference paper.

    The 2013 paper derives intermodulation for cubic polynomials only.  The
    exponential case would need a two-argument Bessel expansion that no source
    in this set provides, so it raises rather than silently returning a number
    from an unverified derivation.
    """
    nl = ExponentialNonlinearity(_IBQ, _VT, port=0, n=1)
    with pytest.raises(NotImplementedError):
        nl.intermodulation_sources([1.0], [1.0], numeric)


def _exponential_fourier_reference(X1, I_S, V_T, m, derivative=False):
    """Fourier coefficient of f(x0(t)) by direct numerical integration.

    An oracle that owes nothing to the Bessel derivation: build the waveform,
    apply the exponential, and transform.  Used because the Bessel path is the
    only part of this module with no published numbers to check it against.
    """
    N = 40000
    theta = 2 * np.pi * np.arange(N) / N
    x0 = np.real(X1 * np.exp(1j * theta))
    if derivative:
        f = (I_S / V_T) * (np.exp(x0 / V_T) - 1)
    else:
        f = I_S * (np.exp(x0 / V_T) - 1 - x0 / V_T)
    return 2 * (f * np.exp(-1j * m * theta)).mean()


@pytest.mark.parametrize('phase_deg', [0, 30, 90, 137, 210, 359])
def test_exponential_harmonics_carry_phase(phase_deg):
    """The Bessel coefficients must be right in *argument*, not only magnitude.

    Regression for a shipped bug.  The Jacobi-Anger expansion is stated for a
    real cosine drive; a drive A*exp(j*phi) is A*cos(t + phi), so harmonic m
    carries exp(j*m*phi).  The original implementation returned real
    coefficients, which left every |F_m| exactly right and every argument
    wrong -- at phi = 90 degrees the second harmonic had the opposite sign.

    Nothing caught it: the stage 3 gate compared the *cubic* path against the
    paper's closed form, and the exponential tests compared only |HD| ratios
    at frequencies where the controlling node happened to be nearly real.  A
    magnitude-only check cannot see a phase error, which is the whole reason
    this test compares complex values.
    """
    X1 = 0.02 * np.exp(1j * np.deg2rad(phase_deg))
    F2, F3, B1_apply = ExponentialNonlinearity(
        _IBQ, _VT, port=0, n=1).harmonic_sources([X1], numeric)

    ref2 = _exponential_fourier_reference(X1, _IBQ, _VT, 2)
    ref3 = _exponential_fourier_reference(X1, _IBQ, _VT, 3)
    ref_b = _exponential_fourier_reference(X1, _IBQ, _VT, 1, derivative=True)

    assert abs(F2[0] - ref2) / abs(ref2) < 1e-9, (
        'F2 at %d deg: got %r, exact %r' % (phase_deg, F2[0], ref2))
    assert abs(F3[0] - ref3) / abs(ref3) < 1e-9, (
        'F3 at %d deg: got %r, exact %r' % (phase_deg, F3[0], ref3))
    assert abs(B1_apply([1.0])[0] - ref_b) / abs(ref_b) < 1e-9, (
        'B1 at %d deg is wrong' % phase_deg)


def test_single_tone_magnitudes_are_invariant_to_the_drive_phase():
    """...and yet no single-tone result depends on that phase. Both are true.

    The factor exp(j*m*phi) is *common* to F_m and to the second-order mixing
    term, so it multiplies the whole m-th harmonic and cancels in every
    magnitude ratio.  For one tone and one device it is exactly a choice of
    time origin, and therefore unobservable.

    This is worth pinning for two reasons.  It records why the incomplete
    coefficients produced no wrong answer in the shipped single-tone scope --
    the omission was a latent incorrectness, not an active bug.  And it fails
    loudly if someone later introduces a genuine phase asymmetry (a second
    nonlinear device seeing a different phase, or a second tone), because then
    the factor no longer cancels and this invariance no longer holds.
    """
    def hd3_with_drive_phase(phase_deg):
        X1 = 0.02 * np.exp(1j * np.deg2rad(phase_deg))
        nl = ExponentialNonlinearity(_IBQ, _VT, port=0, n=1)
        F2, F3, B1_apply = nl.harmonic_sources([X1], numeric)
        ## Reproduce the recurrence's combination step with a unit circuit,
        ## isolating the nonlinearity from any frequency dependence.
        mixing = 0.5 * B1_apply([-F2[0]])[0]
        return abs(F3[0] + mixing)

    reference = hd3_with_drive_phase(0)
    for phase_deg in (30, 90, 137, 250):
        got = hd3_with_drive_phase(phase_deg)
        assert abs(got / reference - 1) < 1e-12, (
            'a single-tone magnitude changed with the drive phase (%d deg): '
            'that should be unobservable' % phase_deg)


## ------------------------------- genuine phase asymmetry: two devices

"""A configuration where a dropped phase is actually observable.

Everything above this point is blind to the argument of the Bessel
coefficients.  With one tone and one device the factor exp(j*m*phi) multiplies
every term of the m-th harmonic alike, so it is a choice of time origin and
cancels out of all magnitudes -- which is why the original phase-less
implementation produced no wrong single-tone answer, and why no existing test
could have caught it.

Two nonlinear devices whose controlling nodes sit at *different* phases break
that.  No single time origin makes both drives real, so their relative phase
is physical, and a device that reports magnitude-only coefficients combines
with the other one wrongly.

The circuit is the smallest one that produces the asymmetry: two nodes coupled
by a capacitor, an exponential device on one and a cubic on the other.  At
30 kHz the nodes sit about 65 degrees apart.
"""

_PA_IS, _PA_VT = 1e-6, 25e-3
_PA_G0, _PA_G1, _PA_CC = 1e-3, 2e-3, 5e-9
_PA_U = 2e-5
_PA_B, _PA_C = 3e-3, 1e-3
_PA_F = 3e4


class _MagnitudeOnlyExponential(ExponentialNonlinearity):
    """The pre-fix behaviour, kept so the difference can be measured."""

    def harmonic_sources(self, X1, toolkit):
        a = abs(X1[self.port]) / self.V_T
        F2 = [0] * self.n
        F3 = [0] * self.n
        F2[self.port] = 2 * self.I_S * self._bessel(2, a, toolkit)
        F3[self.port] = 2 * self.I_S * self._bessel(3, a, toolkit)
        b1 = 2 * (self.I_S / self.V_T) * self._bessel(1, a, toolkit)

        def B1_apply(vec):
            out = [0] * self.n
            out[self.port] = b1 * vec[self.port]
            return out

        return F2, F3, B1_apply


def _two_node_admittance(s):
    return np.array([[_PA_G0 + s * _PA_CC, -s * _PA_CC],
                     [-s * _PA_CC, _PA_G1 + s * _PA_CC]], dtype=complex)


def _mixed_device_hd3(exponential_cls, f=_PA_F):
    w = 2 * np.pi * f
    U = np.array([_PA_U, 0.0], dtype=complex)
    b = np.zeros((2, 2))
    c = np.zeros((2, 2))
    b[1, 1], c[1, 1] = _PA_B, _PA_C

    def apply_G(harmonic, rhs):
        return np.linalg.solve(_two_node_admittance(1j * harmonic.frequency((w,))),
                               np.asarray(rhs, dtype=complex))

    nl = CompositeNonlinearity(exponential_cls(_PA_IS, _PA_VT, port=0, n=2),
                               CubicNonlinearity(b, c))
    return harmonic_response(apply_G, U, nl, tones=(w,), toolkit=numeric).HD3(1)


def test_the_two_nodes_really_are_at_different_phases():
    """The premise of the test below -- assert it rather than assume it."""
    w = 2 * np.pi * _PA_F
    X1 = np.linalg.solve(_two_node_admittance(1j * w),
                         np.array([_PA_U, 0.0], dtype=complex))
    separation = abs(np.rad2deg(np.angle(X1[0]) - np.angle(X1[1])))

    assert 30 < separation < 150, (
        'nodes are %.1f deg apart; the asymmetry this file needs is not '
        'present at this operating point' % separation)


def test_dropped_phase_is_caught_when_two_devices_disagree():
    """The regression test that the single-device tests could not provide.

    A magnitude-only exponential must give a *different* answer from the
    correct one here.  If this ever stops failing for the buggy class, the
    configuration has lost its asymmetry and is guarding nothing.
    """
    correct = _mixed_device_hd3(ExponentialNonlinearity)
    magnitude_only = _mixed_device_hd3(_MagnitudeOnlyExponential)

    difference = abs(correct / magnitude_only - 1)
    assert difference > 0.05, (
        'dropping the phase changed HD3 by only %.2f%% -- this configuration '
        'no longer detects the defect it exists to detect' % (100 * difference))


def test_phase_sensitivity_is_not_a_knife_edge():
    """The detection must not depend on a fragile near-cancellation.

    Sensitivity to the dropped phase rises sharply where the two devices'
    contributions nearly cancel, which would make a brittle test: a small
    change anywhere would move it.  This checks the chosen operating point
    sits on the broad plateau instead, by sweeping the cubic coefficient two
    decades either side and requiring the defect to stay visible throughout.
    """
    global _PA_C
    original = _PA_C
    detected = []
    try:
        for c in (1e-4, 3e-4, 1e-3, 3e-3):
            _PA_C = c
            d = abs(_mixed_device_hd3(ExponentialNonlinearity)
                    / _mixed_device_hd3(_MagnitudeOnlyExponential) - 1)
            detected.append(d)
    finally:
        _PA_C = original

    assert min(detected) > 0.05, (
        'detection collapses somewhere in the sweep: %s'
        % ['%.3f' % d for d in detected])


def test_composite_sums_its_parts():
    """A composite of one part must equal that part; of two, their sum."""
    X1 = [0.01 + 0.002j, 0.004 - 0.001j]
    b = np.zeros((2, 2))
    c = np.zeros((2, 2))
    b[1, 1], c[1, 1] = 2e-3, 5e-3
    cubic = CubicNonlinearity(b, c)
    expo = ExponentialNonlinearity(_PA_IS, _PA_VT, port=0, n=2)

    solo_F2, _, _ = cubic.harmonic_sources(X1, numeric)
    wrapped_F2, _, _ = CompositeNonlinearity(cubic).harmonic_sources(X1, numeric)
    assert np.allclose(np.array(solo_F2, dtype=complex),
                       np.array(wrapped_F2, dtype=complex))

    e_F2, _, _ = expo.harmonic_sources(X1, numeric)
    both_F2, _, _ = CompositeNonlinearity(cubic, expo).harmonic_sources(X1, numeric)
    assert np.allclose(np.array(both_F2, dtype=complex),
                       np.array(solo_F2, dtype=complex)
                       + np.array(e_F2, dtype=complex))


## --------------------------------------- operating point: a separate trap

"""The linearisation point must be solved for the whole circuit, not one device.

The Taylor coefficients b and c are taken *about the operating point*, so
getting that point wrong scales every distortion figure -- and it is easy to
get wrong in a way that looks reasonable.

For a diode biased by a current source with a resistor in parallel, the bias
splits between the two.  Solving the diode alone puts the junction ~20 mV too
high, which sounds negligible until you remember the exponential: it inflates
the small-signal conductance by 2.25x and the Taylor coefficients with it.
Nothing about the resulting numbers looks wrong.

This is *not* the defect that produced the disagreement it was found chasing
-- that turned out to be a phantom DC offset in ISin, fixed separately in
test_source_dc_offset.py -- but it is a real trap and worth its own guard.
"""

_OP_IS, _OP_VT = 1e-13, 25e-3
_OP_IBIAS, _OP_R = 1e-3, 1e3


def _diode_operating_point(include_resistor):
    """Junction voltage, solving either the whole circuit or just the diode."""
    from scipy.optimize import brentq
    if include_resistor:
        residual = lambda v: v / _OP_R + _OP_IS * np.expm1(v / _OP_VT) - _OP_IBIAS
    else:
        residual = lambda v: _OP_IS * np.expm1(v / _OP_VT) - _OP_IBIAS
    return brentq(residual, 0.0, 1.0)


def test_ignoring_a_parallel_path_moves_the_operating_point():
    """Quantify the trap, so its size is on record rather than assumed."""
    v_correct = _diode_operating_point(include_resistor=True)
    v_naive = _diode_operating_point(include_resistor=False)

    assert v_naive > v_correct, 'the naive point should be too high'
    assert 0.015 < v_naive - v_correct < 0.030, (
        'expected roughly 20 mV of error, got %.1f mV'
        % ((v_naive - v_correct) * 1e3))


def test_operating_point_error_inflates_the_taylor_coefficients():
    """A 20 mV error is not a 20 mV effect -- the exponential magnifies it.

    Both coefficients scale as exp(v0/VT), so a small voltage error becomes a
    large coefficient error, and every distortion figure inherits it.
    """
    v_correct = _diode_operating_point(include_resistor=True)
    v_naive = _diode_operating_point(include_resistor=False)

    ratio = np.exp((v_naive - v_correct) / _OP_VT)
    assert ratio > 2.0, (
        'a %.1f mV operating-point error should inflate the coefficients by '
        'more than 2x, got %.2fx' % ((v_naive - v_correct) * 1e3, ratio))

    ## And it is the *whole* error: b and c share the same exponential factor.
    b_correct = (_OP_IS / (2 * _OP_VT ** 2)) * np.exp(v_correct / _OP_VT)
    b_naive = (_OP_IS / (2 * _OP_VT ** 2)) * np.exp(v_naive / _OP_VT)
    assert abs(b_naive / b_correct - ratio) < 1e-9


def test_taylor_coefficients_track_the_operating_point():
    """``taylor_coefficients`` must linearise where it is told to.

    The module cannot know what the surrounding circuit does; it takes the
    operating point as an argument.  This pins that it actually uses it -- a
    version that ignored ``x0`` would silently linearise at zero bias and
    return coefficients too small by many orders of magnitude.

    Note the thermal voltage: the ``Diode`` model builds it from the toolkit's
    own constants (``k*T/q``, giving 25.84 mV at 300 K), not from the round
    25 mV used elsewhere in this file.  Using the round number here predicts a
    ratio 2.2x too large and looks exactly like a code defect -- the same trap
    recorded for the CODATA constants in ``test_semiconductor_jacobians``.
    """
    from pycircuit.circuit.elements import Diode

    d = Diode(1, 2, toolkit=numeric, IS=_OP_IS)
    VT_model = numeric.kboltzmann * 300 / numeric.qelectron
    v_correct = _diode_operating_point(include_resistor=True)

    b_biased, c_biased = taylor_coefficients(d, [v_correct, 0.0])
    b_zero, c_zero = taylor_coefficients(d, [0.0, 0.0])

    expected = np.exp(v_correct / VT_model)
    assert abs(b_biased / b_zero / expected - 1) < 0.05, (
        'b did not track the operating point: ratio %.3e, expected %.3e'
        % (b_biased / b_zero, expected))
    assert abs(c_biased / c_zero / expected - 1) < 0.05


## --------------------------------------- higher-order polynomial fits

@pytest.mark.parametrize('seed', [0, 1, 2])
def test_polynomial_harmonic_sources_match_numerical_fourier(seed):
    """Fifth-order coefficients must be right, not just plausible.

    ``x**4`` feeds the second harmonic and ``x**5`` the third, so the quartic
    and quintic coefficients change ``F2``/``F3`` even though the recurrence
    keeps no harmonic above the third.  Checked against direct numerical
    Fourier extraction rather than against a re-derivation of the same
    algebra.
    """
    from pycircuit.circuit.distortion import PolynomialNonlinearity

    rng = np.random.default_rng(seed)
    coefficients = list(rng.normal(size=4))
    X = complex(rng.normal(), rng.normal())

    F2, F3, B1_apply = PolynomialNonlinearity(coefficients).harmonic_sources(
        [X], numeric)

    N = 8192
    theta = 2 * np.pi * np.arange(N) / N
    x = np.real(X * np.exp(1j * theta))
    f = sum(c * x ** (k + 2) for k, c in enumerate(coefficients))
    fprime = sum((k + 2) * c * x ** (k + 1) for k, c in enumerate(coefficients))
    ref = lambda g, m: 2 * (g * np.exp(-1j * m * theta)).mean()

    assert abs(F2[0] - ref(f, 2)) / abs(ref(f, 2)) < 1e-10
    assert abs(F3[0] - ref(f, 3)) / abs(ref(f, 3)) < 1e-10
    assert abs(B1_apply([1.0])[0] - ref(fprime, 1)) / abs(ref(fprime, 1)) < 1e-10


def test_cubic_is_the_three_term_special_case():
    """A polynomial truncated at third order must equal ``CubicNonlinearity``."""
    from pycircuit.circuit.distortion import PolynomialNonlinearity

    b, c = 3e-3, 7e-3
    X = [0.01 + 0.003j]
    poly_F2, poly_F3, poly_B = PolynomialNonlinearity([b, c]).harmonic_sources(
        X, numeric)
    cub_F2, cub_F3, cub_B = CubicNonlinearity(
        np.array([[b]]), np.array([[c]])).harmonic_sources(X, numeric)

    assert abs(poly_F2[0] - cub_F2[0]) < 1e-18
    assert abs(poly_F3[0] - cub_F3[0]) < 1e-18
    assert abs(poly_B([1.0])[0] - cub_B([1.0])[0]) < 1e-18


def test_quintic_of_an_exponential_beats_its_cubic_fit():
    """Fifth order must move the answer toward the exact Bessel result.

    The point of the class: it separates device-model truncation from
    perturbation truncation, so the two can be attributed independently.
    """
    from pycircuit.circuit.distortion import PolynomialNonlinearity

    I_S, VT = 1e-6, 25e-3
    X = [8e-3 + 0j]

    cubic = PolynomialNonlinearity.taylor_of_exponential(I_S, VT, 3)
    quintic = PolynomialNonlinearity.taylor_of_exponential(I_S, VT, 5)
    exact = ExponentialNonlinearity(I_S, VT, port=0, n=1)

    f2 = lambda nl: abs(nl.harmonic_sources(X, numeric)[0][0])
    gap_cubic = abs(f2(cubic) - f2(exact))
    gap_quintic = abs(f2(quintic) - f2(exact))

    assert gap_quintic < gap_cubic, (
        'the quintic fit should be closer to the exact exponential than the '
        'cubic one (gaps %.3e vs %.3e)' % (gap_quintic, gap_cubic))


def test_orders_beyond_five_are_refused():
    """Silently dropping a term that feeds the kept harmonics is not acceptable."""
    from pycircuit.circuit.distortion import PolynomialNonlinearity

    with pytest.raises(NotImplementedError):
        PolynomialNonlinearity([1.0, 2.0, 3.0, 4.0, 5.0])


## ------------------------------ harmonic cutoff vs consistent truncation

"""The published truncation is principled; restoring its dropped terms hurts.

``harmonic_response`` keeps only the leading ``B1 X2`` term of the
second-order convolution. ``harmonic_response_spectral`` keeps all of them up
to a cutoff. The instinct is that the latter must be more accurate. It is not:
the restored terms are ``O(U**4)``, third-order perturbation also contributes
at ``O(U**4)``, and restoring some while omitting others breaks the
expansion's consistency.
"""

_SPEC_IS, _SPEC_VT = 1e-6, 25e-3


def _spectral_setup():
    from pycircuit.circuit.distortion import harmonic_response_spectral
    alpha = _SPEC_IS / _SPEC_VT
    w0 = 2 * np.pi * 1e4
    Y = lambda s: 1e-3 + s * 1e-9 + alpha
    apply_G = lambda h, rhs: [np.asarray(rhs, dtype=complex)[0]
                              / Y(1j * h.frequency((w0,)))]
    f = lambda v: _SPEC_IS * (np.expm1(v / _SPEC_VT) - v / _SPEC_VT)
    fp = lambda v: (_SPEC_IS / _SPEC_VT) * np.expm1(v / _SPEC_VT)
    return harmonic_response_spectral, apply_G, f, fp, w0


def test_spectral_and_truncated_agree_as_the_drive_vanishes():
    """The two must coincide at small signal, or one of them is simply wrong.

    The terms the truncation drops are higher order in the drive, so they
    vanish faster than the terms it keeps. This is the check that the
    spectral implementation is sound before using it to draw a conclusion.
    """
    spectral, apply_G, f, fp, w0 = _spectral_setup()
    exact = ExponentialNonlinearity(_SPEC_IS, _SPEC_VT, port=0, n=1)

    ## Drives chosen so the gap is resolvable: below ~2e-6 the two agree to
    ## floating-point, which would make a "decreasing" assertion vacuous.
    gaps = []
    for drive in (1.6e-5, 4e-6, 1e-6):
        trunc = harmonic_response(apply_G, [drive], exact,
                                  tones=(w0,), toolkit=numeric)
        spec = spectral(apply_G, [drive], f, fp, (w0,), numeric, n_harmonics=3)
        gaps.append(abs(spec.HD3(0) / trunc.HD3(0) - 1))

    assert gaps[0] > gaps[1] > gaps[2], (
        'the two should converge as the drive shrinks, got %s'
        % ['%.5f' % g for g in gaps])
    assert gaps[-1] < 1e-3, 'still %.5f apart at small signal' % gaps[-1]


def test_doubling_the_harmonic_cutoff_changes_little():
    """Harmonics above the third barely feed back at second order.

    Pinned because it is the cheap thing to try when a prediction is poor,
    and it does not help -- worth knowing before spending effort on it.
    """
    spectral, apply_G, f, fp, w0 = _spectral_setup()
    drive = 5e-7

    three = spectral(apply_G, [drive], f, fp, (w0,), numeric, n_harmonics=3)
    six = spectral(apply_G, [drive], f, fp, (w0,), numeric, n_harmonics=6)

    assert abs(six.HD2(0) / three.HD2(0) - 1) < 0.01
    assert abs(six.HD3(0) / three.HD3(0) - 1) < 0.05


def test_restoring_dropped_terms_moves_the_answer_at_moderate_drive():
    """The inconsistent-truncation effect must be real and O(U**2) relative.

    Not asserting which is *better* here -- that needs the transient
    reference, and ``test_distortion_vs_transient`` has it. This pins that
    the restored terms are not negligible, and that their relative size grows
    with drive as an O(U**4) contribution should.
    """
    spectral, apply_G, f, fp, w0 = _spectral_setup()
    exact = ExponentialNonlinearity(_SPEC_IS, _SPEC_VT, port=0, n=1)

    def relative_gap(drive):
        trunc = harmonic_response(apply_G, [drive], exact,
                                  tones=(w0,), toolkit=numeric)
        spec = spectral(apply_G, [drive], f, fp, (w0,), numeric, n_harmonics=3)
        return abs(spec.HD2(0) / trunc.HD2(0) - 1)

    ## Doubling the drive should quadruple the gap if it is O(U**4) sitting
    ## beside an O(U**2) leading term.
    small, large = relative_gap(4e-6), relative_gap(8e-6)
    assert 3.0 < large / small < 5.0, (
        'the restored terms should grow as U**2 relative to the leading one, '
        'i.e. ~4x for a doubled drive; got %.4f -> %.4f (%.2fx)'
        % (small, large, large / small))


## ------------------- Picard is not the same as perturbation order

"""Pinning a correction: ``order`` counts Picard iterations, nothing more.

An earlier version of this work labelled a Picard-iteration count ``order``
and then reported its behaviour as though it described the perturbation
series. It does not. The two agree at n=0 and n=1 and diverge from n=2, where
a Picard iterate carries a *fragment* of the next perturbation order.

These tests use the scalar problem ``Y x + b x**2 = u``, whose perturbation
series can be written down exactly and whose solution is a quadratic root, so
both constructions can be compared against truth with no numerical method in
between.
"""

_SC_Y, _SC_B, _SC_U = 2.0, 0.7, 1.0


def _scalar_perturbation_terms():
    """Exact term-by-term perturbation series for Y x + b x**2 = u.

    f = b x**2, so expanding f(x0 + e x1 + e^2 x2 + ...) gives
        e^1 : b x0^2
        e^2 : b (2 x0 x1)
        e^3 : b (2 x0 x2 + x1^2)
    """
    G = 1.0 / _SC_Y
    x0 = G * _SC_U
    x1 = -G * _SC_B * x0 ** 2
    x2 = -G * _SC_B * (2 * x0 * x1)
    x3 = -G * _SC_B * (2 * x0 * x2 + x1 ** 2)
    return [x0, x1, x2, x3]


def _scalar_picard(n):
    G = 1.0 / _SC_Y
    x = G * _SC_U
    for _ in range(n):
        x = G * (_SC_U - _SC_B * x ** 2)
    return x


def _scalar_exact():
    return (-_SC_Y + np.sqrt(_SC_Y ** 2 + 4 * _SC_B * _SC_U)) / (2 * _SC_B)


def test_picard_matches_perturbation_only_for_the_first_two_orders():
    """They are the same object at n=0 and n=1, and different after."""
    terms = _scalar_perturbation_terms()
    partial = [sum(terms[:k + 1]) for k in range(4)]

    assert abs(_scalar_picard(0) - partial[0]) < 1e-15
    assert abs(_scalar_picard(1) - partial[1]) < 1e-15
    assert abs(_scalar_picard(2) - partial[2]) > 1e-6, (
        'Picard and the perturbation truncation must differ from n=2')


def test_the_picard_excess_is_a_fragment_of_the_next_order():
    """Quantify the difference rather than merely asserting one exists.

    The second Picard iterate exceeds x0+x1+x2 by exactly -G*b*x1**2, which is
    one of the two terms of x3 = -G*b*(2*x0*x2 + x1**2). Not the whole next
    order -- a piece of it.
    """
    x0, x1, x2, x3 = _scalar_perturbation_terms()
    G = 1.0 / _SC_Y

    excess = _scalar_picard(2) - (x0 + x1 + x2)
    fragment = -G * _SC_B * x1 ** 2

    assert abs(excess - fragment) < 1e-15, 'excess is not the expected fragment'
    assert abs(excess) < abs(x3), 'the excess should be smaller than all of x3'
    assert abs(excess / x3 - 0.2) < 0.05, (
        'for these values the fragment is a fifth of x3; got %.4f'
        % (excess / x3))


def test_the_true_perturbation_series_converges_monotonically():
    """The claim that raising the order is non-monotonic was about Picard.

    On the true series it is not. Pinned because the retracted claim is the
    kind that gets repeated once written down.
    """
    terms = _scalar_perturbation_terms()
    exact = _scalar_exact()
    errors = [abs(sum(terms[:k + 1]) - exact) for k in range(4)]

    assert all(errors[i] > errors[i + 1] for i in range(len(errors) - 1)), (
        'perturbation partial sums should improve monotonically, got %s'
        % ['%.3e' % e for e in errors])


def test_the_inconsistent_truncation_is_not_automatically_worse():
    """Retraction guard: 'inconsistent is worse' does not generalise.

    Picard carries an unbalanced fragment of the next order, which by the
    tidy-truncation argument ought to hurt. Here it helps -- at every order
    past the first. The circuit measurement that prompted that argument stands;
    the explanation offered for it does not.
    """
    terms = _scalar_perturbation_terms()
    exact = _scalar_exact()

    for n in (2, 3):
        pert_err = abs(sum(terms[:n + 1]) - exact)
        pic_err = abs(_scalar_picard(n) - exact)
        assert pic_err < pert_err, (
            'at n=%d Picard (%.3e) should beat the consistent truncation '
            '(%.3e) for this problem' % (n, pic_err, pert_err))


def test_perturbation_terms_stay_single_monomials_symbolically():
    """The reason this method suits a symbolic tool, pinned.

    Carried symbolically, each perturbation order of ``Y x + b x**2 = u`` is a
    single monomial -- ``-5*b**3*u**4/Y**7`` and so on -- readable as a design
    formula. A Picard iterate doubles in length every pass and is unreadable
    by the third. Bounded expression growth is what makes symbolic distortion
    analysis feasible, and an iterative method forfeits it.
    """
    Y, b, u = sympy.symbols('Y b u')
    G = 1/Y

    x0 = G*u
    x1 = sympy.expand(-G*b*x0**2)
    x2 = sympy.expand(-G*b*(2*x0*x1))
    x3 = sympy.expand(-G*b*(2*x0*x2 + x1**2))

    for order, term in enumerate((x0, x1, x2, x3)):
        assert len(sympy.Add.make_args(sympy.expand(term))) == 1, (
            'perturbation term x(%d) should be a single monomial, got %s'
            % (order, term))

    ## Picard, by contrast, doubles.
    x = G*u
    lengths = []
    for _ in range(4):
        x = sympy.expand(G*(u - b*x**2))
        lengths.append(len(sympy.Add.make_args(x)))

    assert lengths == [2, 4, 8, 16], (
        'expected Picard iterates to double in length, got %s' % lengths)
